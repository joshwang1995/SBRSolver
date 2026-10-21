#!/usr/bin/env python3
"""Generate a dense, geometry-aware receiver layout for an SBRSolver scene.

Reads the scene's STL file(s), lays a regular grid (or a line) of receivers and
keeps only the points that make sense for a coverage map:

* inside the region asked for (``--region``):
  ``reachable`` (default)  the area the transmitter's rays can reach without
                            passing through a surface: a flood fill on a raster
                            of the walls, seeded at the transmitter;
  ``enclosure``            everything inside the outermost walls, enclosed
                            rooms included (a flood fill seeded outside the
                            scene marks the exterior);
  ``bbox``                 no wall test, only the box given with ``--bbox``
                            (open scenes such as a ground plane or a corner);
* at least ``--margin`` metres (2D) away from every wall that spans the
  receiver height, so no receiver sits on or inside a surface;
* at least ``--tx-clearance`` metres (3D) away from the transmitter.

The grid is offset by half a spacing from the scene's minimum corner, so a
wall on a round coordinate (x = 2.0, y = 3.0 ...) never coincides with a
receiver. ``--line`` places receivers along a segment instead of a grid.

Output: an SBRSolver receiver file (title line, count, then ``x y z`` per
line) for ``[receivers] layout = file``. ``--plot`` also writes a top-down PNG
of the walls, the transmitter and the receivers (needs Pillow).

Examples (run from the repository root):

  python3 tools/generate_receivers.py --stl SBRSolver_rev2/data/stl_files/office.stl \\
      --stl SBRSolver_rev2/data/stl_files/office_glass_panel.stl \\
      --tx 0.835938 4.53906 2.5 --height 1.5 --spacing 0.25 --region enclosure \\
      --out SBRSolver_rev2/data/fine/RX_Office_fine.dat --plot office_fine.png

  python3 tools/generate_receivers.py --stl SBRSolver_rev2/data/stl_files/ground.stl \\
      --tx 0 0 5 --line 0 1 3 0 999 3 --spacing 0.1 --region bbox \\
      --out SBRSolver_rev2/data/fine/RX_Ground_fine.dat
"""
import argparse
import re
import struct
import sys
from collections import deque

import numpy as np


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
def read_stl(path):
    """Triangles of a binary or ASCII STL as an (n, 3, 3) float array."""
    data = open(path, 'rb').read()
    if data[:5] == b'solid' and b'facet' in data[:512]:
        text = data.decode('ascii', 'replace')
        tris = []
        for m in re.finditer(r'outer loop(.*?)endloop', text, re.S):
            vs = re.findall(r'vertex\s+(\S+)\s+(\S+)\s+(\S+)', m.group(1))
            if len(vs) == 3:
                tris.append([[float(c) for c in v] for v in vs])
        return np.array(tris, dtype=float).reshape(-1, 3, 3)
    n = struct.unpack_from('<I', data, 80)[0]
    rec = np.dtype([('n', '<3f4'), ('v', '<9f4'), ('a', '<u2')])
    arr = np.frombuffer(data[84:84 + 50 * n], dtype=rec)
    return arr['v'].reshape(n, 3, 3).astype(float)


def wall_segments(tris, height):
    """2D segments (x0, y0, x1, y1) of the vertical faces that span `height`.

    Every edge of a vertical triangle is projected to the plane; horizontal
    edges give the wall's footprint, the vertical ones degenerate to points.
    Faces entirely below or above the receiver height are ignored (a desk, a
    half-height wall under the receivers), faces crossing it count.
    """
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    norm = np.linalg.norm(n, axis=1)
    ok = norm > 0
    vertical = ok & (np.abs(n[:, 2]) <= 1e-6 * np.maximum(norm, 1e-300))
    zmin = tris[:, :, 2].min(axis=1)
    zmax = tris[:, :, 2].max(axis=1)
    spans = vertical & (zmin <= height + 1e-9) & (zmax >= height - 1e-9)
    segs = []
    for t in tris[spans]:
        for a, b in ((0, 1), (1, 2), (2, 0)):
            p, q = t[a], t[b]
            if abs(p[0] - q[0]) > 1e-9 or abs(p[1] - q[1]) > 1e-9:
                segs.append((p[0], p[1], q[0], q[1]))
    return np.array(segs, dtype=float).reshape(-1, 4)


def distance_to_segments(points, segs):
    """Shortest 2D distance from each point (n, 2) to any segment (m, 4)."""
    if len(segs) == 0:
        return np.full(len(points), np.inf)
    best = np.full(len(points), np.inf)
    a = segs[:, 0:2]
    d = segs[:, 2:4] - a
    dd = np.maximum((d * d).sum(axis=1), 1e-300)
    chunk = max(1, 2_000_000 // max(len(segs), 1))
    for start in range(0, len(points), chunk):
        p = points[start:start + chunk]
        ap = p[:, None, :] - a[None, :, :]                     # (k, m, 2)
        t = np.clip((ap * d[None, :, :]).sum(axis=2) / dd[None, :], 0.0, 1.0)
        closest = a[None, :, :] + t[:, :, None] * d[None, :, :]
        dist = np.linalg.norm(p[:, None, :] - closest, axis=2).min(axis=1)
        best[start:start + chunk] = dist
    return best


# ---------------------------------------------------------------------------
# Raster flood fill
# ---------------------------------------------------------------------------
class Raster:
    """Occupancy raster of the walls; `cell` metres per cell, one cell of
    padding around the scene so the exterior is connected."""

    def __init__(self, bounds, cell):
        self.cell = cell
        self.x0 = bounds[0] - 2 * cell
        self.y0 = bounds[1] - 2 * cell
        self.nx = int(np.ceil((bounds[2] - self.x0) / cell)) + 3
        self.ny = int(np.ceil((bounds[3] - self.y0) / cell)) + 3
        self.wall = np.zeros((self.ny, self.nx), dtype=bool)

    def index(self, x, y):
        return int((x - self.x0) / self.cell), int((y - self.y0) / self.cell)

    def rasterise(self, segs):
        # 8-connected lines (a 4-connected flood fill cannot cross them)
        for x0, y0, x1, y1 in segs:
            i0, j0 = self.index(x0, y0)
            i1, j1 = self.index(x1, y1)
            steps = max(abs(i1 - i0), abs(j1 - j0), 1)
            for s in range(steps + 1):
                i = i0 + (i1 - i0) * s // steps
                j = j0 + (j1 - j0) * s // steps
                if 0 <= i < self.nx and 0 <= j < self.ny:
                    self.wall[j, i] = True

    def flood(self, seeds):
        reached = np.zeros_like(self.wall)
        queue = deque()
        for x, y in seeds:
            i, j = self.index(x, y)
            if not (0 <= i < self.nx and 0 <= j < self.ny):
                raise SystemExit(f'seed ({x}, {y}) lies outside the scene')
            if self.wall[j, i]:
                raise SystemExit(f'seed ({x}, {y}) lies on a wall; move it or lower --cell')
            reached[j, i] = True
            queue.append((i, j))
        while queue:
            i, j = queue.popleft()
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                a, b = i + di, j + dj
                if 0 <= a < self.nx and 0 <= b < self.ny and not reached[b, a] and not self.wall[b, a]:
                    reached[b, a] = True
                    queue.append((a, b))
        return reached

    def lookup(self, mask, points):
        out = np.zeros(len(points), dtype=bool)
        for k, (x, y) in enumerate(points):
            i, j = self.index(x, y)
            out[k] = 0 <= i < self.nx and 0 <= j < self.ny and mask[j, i]
        return out


# ---------------------------------------------------------------------------
# Candidate points
# ---------------------------------------------------------------------------
def grid_points(bounds, spacing, height, offset):
    xs = np.arange(bounds[0] + offset, bounds[2] - offset / 2 + 1e-9, spacing)
    ys = np.arange(bounds[1] + offset, bounds[3] - offset / 2 + 1e-9, spacing)
    X, Y = np.meshgrid(xs, ys)
    pts = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, height)])
    return pts, (len(xs), len(ys))


def line_points(line, spacing):
    p0 = np.array(line[0:3], dtype=float)
    p1 = np.array(line[3:6], dtype=float)
    length = np.linalg.norm(p1 - p0)
    count = int(round(length / spacing)) + 1
    t = np.linspace(0.0, 1.0, count)
    return p0[None, :] + t[:, None] * (p1 - p0)[None, :]


def write_receivers(path, points, title):
    with open(path, 'w') as f:
        f.write(f'{title}\n{len(points)}\n')
        for x, y, z in points:
            f.write(f'{x:.6g} {y:.6g} {z:.6g}\n')


def plot(path, segs, bounds, tx, points, rejected, scale):
    try:
        from PIL import Image, ImageDraw
    except ImportError:
        print('Pillow is not installed, no plot written', file=sys.stderr)
        return
    margin = 1.0
    x0, y0 = bounds[0] - margin, bounds[1] - margin
    x1, y1 = bounds[2] + margin, bounds[3] + margin
    W, H = int((x1 - x0) * scale) + 1, int((y1 - y0) * scale) + 1
    img = Image.new('RGB', (W, H), 'white')
    draw = ImageDraw.Draw(img)

    def px(p):
        return ((p[0] - x0) * scale, (y1 - p[1]) * scale)

    for p in rejected:
        x, y = px(p)
        draw.ellipse([x - 1, y - 1, x + 1, y + 1], fill=(220, 220, 220))
    for p in points:
        x, y = px(p)
        draw.ellipse([x - 1.2, y - 1.2, x + 1.2, y + 1.2], fill=(0, 120, 255))
    for sx0, sy0, sx1, sy1 in segs:
        draw.line([px((sx0, sy0)), px((sx1, sy1))], fill=(0, 0, 0), width=2)
    x, y = px(tx)
    draw.ellipse([x - 4, y - 4, x + 4, y + 4], fill=(255, 0, 0))
    img.save(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--stl', action='append', required=True, help='scene STL file (repeat for several)')
    ap.add_argument('--tx', nargs=3, type=float, required=True, metavar=('X', 'Y', 'Z'), help='transmitter position')
    ap.add_argument('--height', type=float, default=1.5, help='receiver height for the grid [m] (default 1.5)')
    ap.add_argument('--spacing', type=float, default=0.25, help='grid or line spacing [m] (default 0.25)')
    ap.add_argument('--offset', type=float, default=None, help='grid offset from the minimum corner (default spacing/2)')
    ap.add_argument('--margin', type=float, default=0.25, help='minimum 2D distance from any wall [m] (default 0.25)')
    ap.add_argument('--tx-clearance', type=float, default=0.5, help='minimum 3D distance from the transmitter [m] (default 0.5)')
    ap.add_argument('--region', choices=('reachable', 'enclosure', 'bbox'), default='reachable')
    ap.add_argument('--bbox', nargs=4, type=float, metavar=('XMIN', 'YMIN', 'XMAX', 'YMAX'),
                    help='restrict the grid to this box (default: the scene extents)')
    ap.add_argument('--line', nargs=6, type=float, metavar=('X0', 'Y0', 'Z0', 'X1', 'Y1', 'Z1'),
                    help='receivers along this segment instead of a grid')
    ap.add_argument('--cell', type=float, default=0.05, help='raster cell of the wall flood fill [m] (default 0.05)')
    ap.add_argument('--out', required=True, help='receiver file to write')
    ap.add_argument('--title', default=None, help='first line of the receiver file')
    ap.add_argument('--plot', default=None, help='write a top-down PNG of walls, transmitter and receivers')
    ap.add_argument('--plot-scale', type=float, default=20.0, help='pixels per metre in the plot')
    args = ap.parse_args()

    tris = np.concatenate([read_stl(s) for s in args.stl])
    if len(tris) == 0:
        raise SystemExit('no triangles read')
    allpts = tris.reshape(-1, 3)
    extents = (allpts[:, 0].min(), allpts[:, 1].min(), allpts[:, 0].max(), allpts[:, 1].max())
    bounds = tuple(args.bbox) if args.bbox else extents
    tx = np.array(args.tx, dtype=float)

    if args.line:
        candidates = line_points(args.line, args.spacing)
        height = float(np.mean(candidates[:, 2]))
        shape = (len(candidates), 1)
    else:
        offset = args.spacing / 2.0 if args.offset is None else args.offset
        candidates, shape = grid_points(bounds, args.spacing, args.height, offset)
        height = args.height
    segs = wall_segments(tris, height)
    reasons = {}
    keep = np.ones(len(candidates), dtype=bool)

    # 1. region
    if args.region != 'bbox':
        raster = Raster(extents, args.cell)
        raster.rasterise(segs)
        if args.region == 'reachable':
            reached = raster.flood([(tx[0], tx[1])])
            inside = raster.lookup(reached, candidates[:, :2])
        else:
            outside = raster.flood([(extents[0] - args.cell, extents[1] - args.cell)])
            inside = ~raster.lookup(outside, candidates[:, :2])
        reasons['outside the %s region' % args.region] = int((~inside & keep).sum())
        keep &= inside
    # 2. margin from walls
    dist = distance_to_segments(candidates[:, :2], segs)
    near = dist < args.margin
    reasons['within %g m of a wall' % args.margin] = int((near & keep).sum())
    keep &= ~near
    # 3. clearance from the transmitter
    close = np.linalg.norm(candidates - tx[None, :], axis=1) < args.tx_clearance
    reasons['within %g m of the transmitter' % args.tx_clearance] = int((close & keep).sum())
    keep &= ~close

    points = candidates[keep]
    if len(points) == 0:
        raise SystemExit('no receiver survived the filters')
    title = args.title or ('receiver positions (%s, spacing %g m, height %g m, margin %g m, %s)'
                           % (args.region, args.spacing, height, args.margin, ' + '.join(args.stl)))
    write_receivers(args.out, points, title)
    print('%s: %d receivers kept of %d candidates (%d x %d)' % (args.out, len(points), len(candidates), shape[0], shape[1]))
    for reason, count in reasons.items():
        print('  rejected %6d %s' % (count, reason))
    print('  extents x [%g, %g] y [%g, %g] z [%g, %g]' % (points[:, 0].min(), points[:, 0].max(), points[:, 1].min(),
                                                          points[:, 1].max(), points[:, 2].min(), points[:, 2].max()))
    print('  closest wall %.3f m, closest to the transmitter %.3f m'
          % (dist[keep].min() if np.isfinite(dist[keep]).any() else float('inf'),
             np.linalg.norm(points - tx[None, :], axis=1).min()))
    if args.plot:
        plot(args.plot, segs, bounds, tx, points, candidates[~keep], args.plot_scale)
        print('  plot written to', args.plot)


if __name__ == '__main__':
    main()
