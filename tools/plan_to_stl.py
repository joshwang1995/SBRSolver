#!/usr/bin/env python3
"""Convert an architectural floor plan (vector PDF) into SBRSolver STL files.

The plan's line work gives the walls. A new wall is drawn as two parallel
lines a wall's thickness apart. A glazed partition is drawn as one line with
a jamb block at an end. The existing walls, the concrete core and the
columns are drawn as filled black polygons. The tool reads those vectors with
pdfplumber and pairs the parallel lines into wall centre lines. A wall drawn
with three or more lines keeps one centre line per band of lines a wall's
thickness wide. A rectangle wider than a wall outlines an object (a desk, a
table, a cabinet) and gives no lines. The edges of the filled polygons are
not wall lines, because the fill models that wall. A line that runs along a
fill and past its end keeps the part past the end. The tool removes the door
leaves (a short line from the hinge of a door-swing arc), stair treads with
their landing edges and glyph-sized fills. A small fill on the edge of a
wall fill is a chamfer or a jamb and stays. It then sorts the walls into
classes that decide their height and material:

* ``exterior``: walls on the outline of the plate, modelled full height as
  glass. A ray from such a wall leaves the plan without crossing anything.
  A long wall off the 90 degree grid is also exterior.
* ``glazing``: a single-line wall, full height as glass.
* ``partition``: walls that join the wall network or the core, full height.
  Also a free-standing piece in line with such walls within 1.2 m of both
  ends (the pier between two doors).
* ``booth``: a closed box of 1 m or more that joins nothing (a phone booth).
* ``partition`` also: a free-standing group of 4 m or more with a corner.
* ``screen``: a long free-standing wall that joins nothing (a desk screen).
  Also a desk row: a piece of 1 m to 7.5 m with one free end and no door
  swing at either end. It has a parallel twin of about the same length
  within 6 m that also has a free end.
* The tool drops short free-standing pieces (bench ends, symbols) and short
  pieces that hug the curtain wall (sills, mullions, column faces). It also
  drops short pieces off the 90 degree grid (break lines, symbols).
* The tool also drops furniture outlines. A piece under 2 ft with a free
  end and no door at either end is one (a seat, a cabinet, a locker). So
  is a piece under 1 m with a free end and no door beside a parallel piece
  of 0.5 to 1.5 m.

A filled polygon 8 in or thicker is ``core`` (concrete). A thinner one is an
``existing`` wall (plaster board, full height). The tool drops the dashes of
dotted lines and fills narrower than 1.5 in.

Walls and screens become single planes at the centre line (one quad each).
The polygons become extruded prisms (side faces only). Two fills that abut
get no face between them, and a column that overlaps a wall fill writes the
shared face once. The floor and the suspended ceiling cover the
plate outline. The tool writes one STL per surface group: partitions,
curtain wall, glazed partitions, core, screens, floor and ceiling. It also
writes ``walls.json`` (the walls and the polygons, in metres) and
``walls.csv`` (every wall with its class, length, thickness, height and
material). With ``--overlay`` it writes a PNG of the plan with the classes
drawn on it.

Coordinates are in metres with the origin at the south-west corner of the
plan region (``--region``). x points east, y north and z up. The floor is at
z = 0. Give ``--outline`` in this frame. To find the values, run the tool once
without it. Then read the exterior wall ends from ``walls.json``.

Example (a large office plan, 1/8" = 1'-0" on a 48 x 36 in sheet):

  python3 -I tools/plan_to_stl.py layout.pdf --out SBRSolver_rev2/data/large_office \\
      --region 400 250 1560 2480 --overlay overlay.png \\
      --outline "0.38,9.93 27.63,9.92 28.3,9.22 35.2,3.13 37.44,3.26 37.44,9.75 38.3,9.9 \\
                 38.3,65.63 33.47,65.79 33.48,74.74 0.48,74.74 1.14,65.79 0.38,65.63"

Needs pdfplumber, numpy and Pillow.
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import math
import struct
from pathlib import Path

import numpy as np

Point = tuple[float, float]


class PlanToStl:
    """Reads one PDF page and writes the scene files for it."""

    # wall thickness in inches that two parallel lines may stand apart
    THICKNESS_IN = (2.5, 13.0)
    OVERLAY_SCALE = 2           # overlay pixels per PDF point
    MATERIAL = {'exterior': 'glass', 'glazing': 'glass', 'partition': 'plaster_board', 'booth': 'plaster_board',
                'existing': 'plaster_board', 'core': 'concrete', 'screen': 'wood'}

    def __init__(self) -> None:
        parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('pdf', type=Path, help='the plan, a vector PDF (the first page is read)')
        parser.add_argument('--out', type=Path, required=True, help='directory for the STL files and walls.json')
        parser.add_argument('--feet-per-point', type=float, default=1.0 / 9.0,
                            help='plan scale: real feet per PDF point (default 1/9: a 1/8" = 1\'-0" drawing)')
        parser.add_argument('--region', nargs=4, type=float, metavar=('X0', 'Y0', 'X1', 'Y1'),
                            help='part of the sheet that holds the plan, in PDF points (default: the whole page)')
        parser.add_argument('--ceiling', type=float, default=2.7, help='ceiling height and full wall height [m] (default 2.7)')
        parser.add_argument('--screen-height', type=float, default=1.2, help='height of free-standing screens [m] (default 1.2)')
        parser.add_argument('--booth-height', type=float, default=2.2, help='height of free-standing booths [m] (default 2.2)')
        parser.add_argument('--outline', type=str, default=None,
                            help='plate outline for the floor and ceiling as "x,y x,y ..." in metres, counter-clockwise '
                                 '(default: the convex hull of the exterior walls)')
        parser.add_argument('--overlay', type=Path, default=None, help='write a PNG of the plan with the wall classes')
        self.args = parser.parse_args()
        self.pt_per_ft = 1.0 / self.args.feet_per_point
        self.pt_per_m = self.pt_per_ft / 0.3048
        self.join_tol = 0.5 * self.pt_per_ft           # 6 in: ends closer than this are joined
        self.min_gap = self.THICKNESS_IN[0] / 12 * self.pt_per_ft
        self.max_gap = self.THICKNESS_IN[1] / 12 * self.pt_per_ft
        self.polys: list[list[Point]] = []
        self.poly_class: list[str] = []                                # core, existing or drop, per polygon
        self.segments: list[tuple[float, float, float, float]] = []   # angle, offset, lo, hi
        self.singles: list[tuple[float, float, float, float]] = []    # long lines with a jamb block at an end
        self.arcs: list[tuple[Point, float, tuple]] = []               # hinge, radius, bounding box
        self.walls: list[list[float]] = []                             # angle, offset, lo, hi, gap (0: a single line)
        self.records: list[dict] = []
        self.closures: list[tuple[Point, Point]] = []
        self.bounds = (0.0, 0.0, 0.0, 0.0)
        self.plan_image = None

    # ------------------------------------------------------------------ input
    def read_pdf(self) -> None:
        """Loads the filled polygons, straight line pieces and door arcs."""
        import pdfplumber
        with pdfplumber.open(self.args.pdf) as pdf:
            page = pdf.pages[0]
            region = tuple(self.args.region) if self.args.region else (0.0, 0.0, page.width, page.height)
            objects = list(page.lines) + list(page.rects) + list(page.curves)
            if self.args.overlay:
                self.plan_image = page.crop(region).to_image(resolution=72 * self.OVERLAY_SCALE).original
        def black(colour) -> bool:
            return str(colour) in ('0.0', '0', '[0.0]', '[0]')
        inch = self.pt_per_ft / 12
        small: list[list[Point]] = []
        jambs: list[Point] = []
        for o in objects:
            cx, cy = (o['x0'] + o['x1']) / 2, (o['top'] + o['bottom']) / 2
            if not (region[0] <= cx <= region[2] and region[1] <= cy <= region[3]):
                continue
            if o['object_type'] == 'rect':
                pts = [(o['x0'], o['top']), (o['x1'], o['top']), (o['x1'], o['bottom']), (o['x0'], o['bottom']), (o['x0'], o['top'])]
            elif o.get('pts'):
                pts = [(float(p[0]), float(p[1])) for p in o['pts']]
            else:
                pts = [(o['x0'], o['top']), (o['x1'], o['bottom'])]
            if o.get('fill') and black(o.get('non_stroking_color')) and len(pts) >= 3:
                xs = [p[0] for p in pts]
                ys = [p[1] for p in pts]
                poly = pts if pts[0] != pts[-1] else pts[:-1]
                # glyphs and symbols are under 1 ft across, a column is 1.5 ft
                if max(max(xs) - min(xs), max(ys) - min(ys)) < self.pt_per_ft:
                    small.append(poly)
                else:
                    self.polys.append(poly)
            if not o.get('stroke') or not black(o.get('stroking_color')) or (o.get('linewidth') or 0) >= 1.5:
                continue
            if o['object_type'] == 'curve' and self._door_arc(pts):
                continue
            if o['object_type'] == 'rect':
                # a jamb block, the glazing symbol at the end of a single line
                sides = sorted((o['x1'] - o['x0'], o['bottom'] - o['top']))
                if 0.08 * self.pt_per_m <= sides[1] <= 0.13 * self.pt_per_m and sides[0] <= 0.07 * self.pt_per_m:
                    jambs.append((cx, cy))
                # a rect wider than a wall outlines an object (a desk, a
                # table, a cabinet): its edges are not wall lines
                if sides[0] > self.max_gap:
                    continue
            for (ax, ay), (bx, by) in zip(pts, pts[1:]):
                dx, dy = bx - ax, by - ay
                if math.hypot(dx, dy) < self.pt_per_ft:
                    continue
                angle = math.degrees(math.atan2(dy, dx)) % 180.0
                u, n = self._axes(angle)
                ta, tb = ax * u[0] + ay * u[1], bx * u[0] + by * u[1]
                self.segments.append((angle, ax * n[0] + ay * n[1], min(ta, tb), max(ta, tb)))
        self.bounds = region
        # a small fill with up to five corners on the edge of a wall fill is
        # part of that wall: a chamfered corner or a jamb
        edges = [(a, b) for pl in self.polys if self._caliper(pl)[0] >= 1.5 * inch for a, b in zip(pl, pl[1:] + pl[:1])]
        bits = [pl for pl in small if len(pl) <= 5 and any(self._dist(v, a, b) <= inch for v in pl for a, b in edges)]
        self.polys += bits
        self.merge_triangles()
        self.drop_polygon_edges()
        self.singles = [s for s in self.segments if s[3] - s[2] >= 1.5 * self.pt_per_m
                        and any(math.hypot(e[0] - j[0], e[1] - j[1]) <= 0.12 * self.pt_per_m for e in self.endpoints(s) for j in jambs)]
        print(f'{len(self.polys)} filled polygons ({len(small) - len(bits)} glyph-sized dropped), {len(self.segments)} line pieces, '
              f'{len(self.arcs)} door arcs, {len(jambs)} jamb blocks')

    def _door_arc(self, pts: list[Point]) -> bool:
        """Records a quarter-circle door swing (2 to 5 ft) and returns True for one.

        The arc's end points sit on two diagonal corners of its bounding box.
        The hinge is the remaining corner away from the curve.
        """
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        w, h = max(xs) - min(xs), max(ys) - min(ys)
        if not (2 * self.pt_per_ft <= w <= 5 * self.pt_per_ft and 0.8 <= w / max(h, 1e-9) <= 1.25):
            return False
        corners = [(min(xs), min(ys)), (max(xs), min(ys)), (max(xs), max(ys)), (min(xs), max(ys))]
        def corner(p: Point) -> int | None:
            for i, c in enumerate(corners):
                if math.hypot(p[0] - c[0], p[1] - c[1]) < 1.0:
                    return i
            return None
        a, b = corner(pts[0]), corner(pts[-1])
        if a is None or b is None or (a - b) % 4 != 2:
            return False
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
        hinge = max((i for i in range(4) if i not in (a, b)), key=lambda i: math.hypot(corners[i][0] - mx, corners[i][1] - my))
        self.arcs.append((corners[hinge], (w + h) / 2, (corners[0][0], corners[0][1], corners[2][0], corners[2][1])))
        return True

    def merge_triangles(self) -> None:
        """Joins filled triangles to the polygons they share an edge with.

        The exporter splits each fill into triangles. A triangle that shares
        an edge with another polygon adds its apex to that polygon. The
        passes repeat until no triangle joins. A fill drawn twice gives one
        polygon: identical triangles are read once, and a polygon inside
        another one (within 1 in) is dropped.
        """
        def key(p: Point) -> tuple[float, float]:
            return (round(p[0], 2), round(p[1], 2))
        seen: set[frozenset] = set()
        polys = []
        for pl in self.polys:
            k = frozenset(key(p) for p in pl)
            if k not in seen:
                seen.add(k)
                polys.append(pl)
        while True:
            by_edge: dict[tuple, list[int]] = collections.defaultdict(list)
            for i, pl in enumerate(polys):
                for a, b in zip(pl, pl[1:] + pl[:1]):
                    by_edge[tuple(sorted((key(a), key(b))))].append(i)
            gone: set[int] = set()
            for i, pl in enumerate(polys):
                if len(pl) != 3:
                    continue
                for a, b in zip(pl, pl[1:] + pl[:1]):
                    twins = [j for j in by_edge[tuple(sorted((key(a), key(b))))] if j != i and j not in gone]
                    if not twins:
                        continue
                    other = polys[twins[0]]
                    k = next((k for k in range(len(other)) if {key(other[k]), key(other[(k + 1) % len(other)])} == {key(a), key(b)}), None)
                    if k is None:
                        continue
                    apex = next(p for p in pl if key(p) not in (key(a), key(b)))
                    polys[twins[0]] = other[:k + 1] + [apex] + other[k + 1:]
                    gone.add(i)
                    break
            if not gone:
                break
            polys = [pl for i, pl in enumerate(polys) if i not in gone]
        inch = self.pt_per_ft / 12
        order = sorted(range(len(polys)), key=lambda i: -abs(self._area(polys[i])))
        covered: set[int] = set()
        for n, i in enumerate(order):
            for j in order[:n]:
                if j not in covered and all(self._inside(v, polys[j]) or self._dist_to(v, polys[j]) <= inch for v in polys[i]):
                    covered.add(i)
                    break
        self.polys = [pl for i, pl in enumerate(polys) if i not in covered]

    def _inside(self, pt: Point, poly: list[Point]) -> bool:
        x, y = pt
        ins = False
        for (ax, ay), (bx, by) in zip(poly, poly[1:] + poly[:1]):
            if (ay > y) != (by > y) and x < ax + (y - ay) * (bx - ax) / (by - ay):
                ins = not ins
        return ins

    def _dist_to(self, pt: Point, poly: list[Point]) -> float:
        return min(self._dist(pt, a, b) for a, b in zip(poly, poly[1:] + poly[:1]))

    def drop_polygon_edges(self) -> None:
        """Cuts the parts of the line pieces that lie on the edges of the filled polygons.

        A line often runs along a filled wall and continues past its end as
        the face of a hollow wall. The part past the end stays. The end cap
        of a long fill is not a wall face: a line along it belongs to the
        wall the fill abuts, so the cap does not cut it.
        """
        inch = self.pt_per_ft / 12
        edges = []
        for pl in self.polys:
            width, length = self._caliper(pl)[:2]
            if width < 1.5 * inch:
                continue
            for a, b in zip(pl, pl[1:] + pl[:1]):
                if length >= 0.5 * self.pt_per_m and math.hypot(b[0] - a[0], b[1] - a[1]) < self.pt_per_ft:
                    continue
                edges.append((a, b))
        if not edges:
            return
        A = np.array([e[0] for e in edges], dtype=float)
        D = np.array([e[1] for e in edges], dtype=float) - A
        L = np.hypot(D[:, 0], D[:, 1])
        A, D, L = A[L > 1e-9], D[L > 1e-9], L[L > 1e-9]
        U = D / L[:, None]
        N = np.stack([-U[:, 1], U[:, 0]], axis=1)
        angle = np.degrees(np.arctan2(D[:, 1], D[:, 0])) % 180
        offset = (A * N).sum(axis=1)
        kept, cut = [], 0
        for seg in self.segments:
            p, q = self.endpoints(seg)
            u = self._axes(seg[0])[0]
            parallel = np.abs(((angle - seg[0]) + 90) % 180 - 90) <= 1.5
            near = (np.abs(p[0] * N[:, 0] + p[1] * N[:, 1] - offset) <= 1.5) & (np.abs(q[0] * N[:, 0] + q[1] * N[:, 1] - offset) <= 1.5)
            ta = A[:, 0] * u[0] + A[:, 1] * u[1]
            tb = ta + D[:, 0] * u[0] + D[:, 1] * u[1]
            on_edge = parallel & near & (np.minimum(np.maximum(ta, tb), seg[3]) - np.maximum(np.minimum(ta, tb), seg[2]) > 0)
            spans = [(float(min(a, b)), float(max(a, b))) for a, b in zip(ta[on_edge], tb[on_edge])]
            runs = self._uncovered(seg[2], seg[3], spans)
            cut += runs != [(seg[2], seg[3])]
            kept += [(seg[0], seg[1], lo, hi) for lo, hi in runs]
        print(f'{cut} line pieces on polygon edges removed or cut back')
        self.segments = kept

    def _uncovered(self, lo: float, hi: float, spans: list[tuple[float, float]], least: float | None = None) -> list[tuple[float, float]]:
        """The parts of [lo, hi] that no span covers, 1 ft or longer unless ``least`` says otherwise."""
        least = self.pt_per_ft if least is None else least
        out = []
        for a, b in sorted(spans):
            if a > lo:
                out.append((lo, min(a, hi)))
            lo = max(lo, b)
            if lo >= hi:
                break
        if lo < hi:
            out.append((lo, hi))
        return [(a, b) for a, b in out if b - a >= least]

    # ------------------------------------------------------------ wall finding
    def _leaf(self, p: Point, q: Point, length: float, tol: float) -> bool:
        """True for a door leaf: a line from the hinge of a door-swing arc into the arc's box."""
        for hinge, r, box in self.arcs:
            if not 0.7 * r <= length <= 1.2 * r:
                continue
            for near, far in ((p, q), (q, p)):
                if (math.hypot(near[0] - hinge[0], near[1] - hinge[1]) <= tol
                        and box[0] - tol <= far[0] <= box[2] + tol and box[1] - tol <= far[1] <= box[3] + tol):
                    return True
        return False

    def pair_lines(self) -> None:
        """Pairs parallel line pieces a wall's thickness apart into centre lines.

        The door leaves go first: their edges would pair with the wall faces
        next to the hinge. A single line with a jamb block at an end is a
        glazed partition. It becomes a centre line with gap 0 where no pair
        covers it.
        """
        leaves = [s for s in self.segments if self._leaf(*self.endpoints(s), s[3] - s[2], 3.0)]
        self.segments = [s for s in self.segments if s not in leaves]
        bins: dict[int, list] = collections.defaultdict(list)
        for s in self.segments:
            bins[round(s[0]) % 180].append(s)
        seen = set()
        for key, group in bins.items():
            if key in (45, 135):            # the hatch direction of a plan
                continue
            candidates = group + bins.get((key + 1) % 180, []) + bins.get((key - 1) % 180, [])
            for a in group:
                a_off, a_lo, a_hi = self._in_frame(a, a[0])
                for b in candidates:
                    if b is a or not self._parallel(a[0], b[0]):
                        continue
                    b_off, b_lo, b_hi = self._in_frame(b, a[0])
                    gap = abs(b_off - a_off)
                    lo, hi = max(a_lo, b_lo), min(a_hi, b_hi)
                    if self.min_gap <= gap <= self.max_gap and hi - lo >= self.pt_per_ft:
                        k = tuple(sorted((a, b)))
                        if k not in seen:
                            seen.add(k)
                            self.walls.append([a[0], (a_off + b_off) / 2, lo, hi, gap])
        self.walls = self._merge(self._merge(self.walls))
        glazing = []
        for s in self._merge([[s[0], s[1], s[2], s[3], 0.0] for s in self.singles]):
            spans = []
            for w in self.walls:
                off, lo, hi = self._in_frame(w, s[0])
                if self._parallel(w[0], s[0]) and abs(off - s[1]) <= w[4] / 2 + 1.0:
                    spans.append((lo, hi))
            glazing += [[s[0], s[1], lo, hi, 0.0] for lo, hi in self._uncovered(s[2], s[3], spans)]
        self.walls += glazing
        print(f'{len(leaves)} door leaves removed, {len(self.walls) - len(glazing)} wall centre lines, '
              f'{len(glazing)} single-line glazed partitions')

    def _axes(self, angle: float) -> tuple[Point, Point]:
        """Unit vectors along and across a direction in degrees."""
        r = math.radians(angle)
        u = (math.cos(r), math.sin(r))
        return u, (-u[1], u[0])

    def _parallel(self, a: float, b: float) -> bool:
        return abs(((b - a) + 90) % 180 - 90) <= 1.5

    def _in_frame(self, seg: tuple, angle: float) -> tuple[float, float, float]:
        """Offset and extent of a line piece in the frame of another direction."""
        u, n = self._axes(angle)
        p, q = self.endpoints(seg)
        off = ((p[0] + q[0]) / 2) * n[0] + ((p[1] + q[1]) / 2) * n[1]
        ta, tb = p[0] * u[0] + p[1] * u[1], q[0] * u[0] + q[1] * u[1]
        return off, min(ta, tb), max(ta, tb)

    def _merge(self, walls: list[list[float]]) -> list[list[float]]:
        """Joins collinear overlapping centre lines."""
        out: list[list[float]] = []
        by_angle: dict[int, list] = collections.defaultdict(list)
        for w in walls:
            by_angle[round(w[0])].append(w)
        for group in by_angle.values():
            group.sort(key=lambda w: (round(w[1]), w[2]))
            merged: list[list[float]] = []
            for w in group:
                for m in merged:
                    if abs(m[1] - w[1]) <= 1.0 and w[2] <= m[3] + 1.0 and w[3] >= m[2] - 1.0:
                        m[2], m[3], m[4] = min(m[2], w[2]), max(m[3], w[3]), max(m[4], w[4])
                        break
                else:
                    merged.append(list(w))
            out += merged
        return out

    def endpoints(self, w: list[float] | tuple) -> tuple[Point, Point]:
        u, n = self._axes(w[0])
        return ((w[2] * u[0] + w[1] * n[0], w[2] * u[1] + w[1] * n[1]),
                (w[3] * u[0] + w[1] * n[0], w[3] * u[1] + w[1] * n[1]))

    def drop_non_walls(self) -> None:
        """Removes door leaves, the parts inside the filled polygons, stair treads and landing edges."""
        boxes = [self._bbox(pl) for pl in self.polys]
        step = 0.25 * self.pt_per_ft
        kept, leaves, in_core, clipped = [], 0, 0, 0
        for w in self.walls:
            p, q = self.endpoints(w)
            length = w[3] - w[2]
            # a leaf drawn as two lines pairs with itself: the centre line
            # then lies half the gap from the hinge
            if self._leaf(p, q, length, 3.0 + w[4] / 2):
                leaves += 1
                continue
            # sample the centre line every 3 in: the parts over a filled
            # polygon are modelled by its prism, the others stay walls
            box = (min(p[0], q[0]) - 1, min(p[1], q[1]) - 1, max(p[0], q[0]) + 1, max(p[1], q[1]) + 1)
            near = [pl for pl, b in zip(self.polys, boxes) if b[0] <= box[2] and b[2] >= box[0] and b[1] <= box[3] and b[3] >= box[1]]
            count = max(1, math.ceil(length / step))
            outside = []
            for k in range(count):
                t = (k + 0.5) / count
                outside.append(not any(self._inside((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])), pl) for pl in near))
            runs = [(k0, k1) for k0, k1 in self._runs(outside) if (k1 - k0) * length / count >= self.pt_per_ft]
            if not runs:
                in_core += 1
            elif runs == [(0, count)]:
                kept.append(w)
            else:
                clipped += 1
                kept += [[w[0], w[1], w[2] + k0 * length / count, w[2] + k1 * length / count, w[4]] for k0, k1 in runs]
        # treads: a run of four or more short parallel pieces, each 6 in to
        # 2 ft from the next, that overlap along the flight and span 2 ft or
        # more across it (a wall drawn with five lines spans less)
        short = 7 * self.pt_per_ft
        links = []
        for i, w in enumerate(kept):
            for j in range(i + 1, len(kept)):
                v = kept[j]
                if w[3] - w[2] > short or v[3] - v[2] > short or round(v[0]) != round(w[0]):
                    continue
                if (0.5 * self.pt_per_ft <= abs(v[1] - w[1]) <= 2 * self.pt_per_ft
                        and min(v[3], w[3]) - max(v[2], w[2]) > 0.5 * min(v[3] - v[2], w[3] - w[2])):
                    links.append((i, j))
        roots = self._roots(len(kept), links)
        runs: dict[int, list[int]] = collections.defaultdict(list)
        for i, w in enumerate(kept):
            if w[3] - w[2] <= short:
                runs[roots[i]].append(i)
        flights = [idx for idx in runs.values() if len(idx) >= 4
                   and max(kept[j][1] for j in idx) - min(kept[j][1] for j in idx) >= 2 * self.pt_per_ft]
        treads = {i for idx in flights for i in idx}
        # a landing edge: a piece with both ends within 1 ft of a flight's box
        ends = [self.endpoints(w) for w in kept]
        for idx in flights:
            x0, y0, x1, y1 = self._bbox([e for j in idx for e in ends[j]])
            treads.update(i for i in range(len(kept)) if all(x0 - self.pt_per_ft <= e[0] <= x1 + self.pt_per_ft
                                                             and y0 - self.pt_per_ft <= e[1] <= y1 + self.pt_per_ft for e in ends[i]))
        self.walls = [w for i, w in enumerate(kept) if i not in treads]
        print(f'dropped {leaves} door leaves, {in_core} pieces inside the filled polygons ({clipped} clipped), '
              f'{len(treads)} stair treads and landing edges: {len(self.walls)} walls')

    def _roots(self, n: int, links: list[tuple[int, int]]) -> list[int]:
        """The group root of each of n items after joining the linked pairs."""
        parent = list(range(n))
        def find(i: int) -> int:
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i
        for i, j in links:
            parent[find(i)] = find(j)
        return [find(i) for i in range(n)]

    def _runs(self, flags: list[bool]) -> list[tuple[int, int]]:
        """Index ranges [start, end) of the consecutive True values."""
        runs = []
        start = None
        for k, flag in enumerate(flags + [False]):
            if flag and start is None:
                start = k
            elif not flag and start is not None:
                runs.append((start, k))
                start = None
        return runs

    def drop_stacked(self) -> None:
        """Keeps one centre line where a wall is drawn with three or more lines.

        Every pair of the lines gives a centre line. The longest one (to the
        foot, the thinnest on a tie) stays. A parallel one whose thickness
        band meets the kept one's band is cut back to the parts the kept one
        does not cover. A single-line wall (gap 0) covers nothing.
        """
        order = sorted(range(len(self.walls)), key=lambda i: (-round((self.walls[i][3] - self.walls[i][2]) / self.pt_per_ft), self.walls[i][4]))
        kept: list[tuple[int, list[float]]] = []
        dropped = 0
        for i in order:
            w = self.walls[i]
            covered = [False] * int(math.ceil(w[3] - w[2]))
            for _, v in kept:
                if v[4] == 0 or not self._parallel(v[0], w[0]):
                    continue
                off, lo, hi = self._in_frame(v, w[0])
                if abs(off - w[1]) > (v[4] + w[4]) / 2 + 1.0:
                    continue
                for k in range(max(0, int(lo - w[2])), min(len(covered), int(math.ceil(hi - w[2])))):
                    covered[k] = True
            if not any(covered):
                kept.append((i, w))
                continue
            dropped += 1
            for k0, k1 in self._runs([not c for c in covered]):
                if k1 - k0 >= self.pt_per_ft:
                    kept.append((i, [w[0], w[1], w[2] + k0, min(w[3], w[2] + k1), w[4]]))
        self.walls = [w for _, w in sorted(kept, key=lambda item: item[0])]
        print(f'{dropped} stacked centre lines removed or cut back: {len(self.walls)} walls')

    # ----------------------------------------------------- filled polygons
    def sort_polys(self) -> None:
        """Sorts the filled polygons into core, existing wall and drop.

        The core and the columns are 8 in and thicker. A thinner fill is an
        existing wall. A fill under 1.5 in wide is a drawn line. A short thin
        fill that touches no wall, no long fill and no kept fill is a dash of
        a dotted line.
        """
        inch = self.pt_per_ft / 12
        ends = [self.endpoints(w) for w in self.walls]
        sizes = [self._caliper(pl)[:2] for pl in self.polys]
        def touches(i: int) -> bool:
            for v in self.polys[i]:
                if any(self._dist(v, a, b) <= self.join_tol for a, b in ends):
                    return True
                for j, pl in enumerate(self.polys):
                    if j == i or sizes[j][1] < 0.5 * self.pt_per_m or sizes[j][0] < 1.5 * inch:
                        continue
                    if any(self._dist(v, a, b) <= self.join_tol for a, b in zip(pl, pl[1:] + pl[:1])):
                        return True
            return False
        self.poly_class = []
        for i, pl in enumerate(self.polys):
            width, length = sizes[i]
            if width < 1.5 * inch or abs(self._area(pl)) < 1e-6:
                cls = 'drop'
            elif width >= 8 * inch:
                cls = 'core'
            elif length < 0.5 * self.pt_per_m and not touches(i):
                cls = 'drop'
            else:
                cls = 'existing'
            self.poly_class.append(cls)
        # a short fill on the edge of a kept fill is part of that wall (a
        # jamb on a chamfer), so a chain of short fills from a wall stays
        changed = True
        while changed:
            kept = [(a, b) for pl, cls in zip(self.polys, self.poly_class) if cls != 'drop' for a, b in zip(pl, pl[1:] + pl[:1])]
            changed = False
            for i, pl in enumerate(self.polys):
                if self.poly_class[i] != 'drop' or sizes[i][0] < 1.5 * inch or abs(self._area(pl)) < 1e-6:
                    continue
                if any(self._dist(v, a, b) <= inch for v in pl for a, b in kept):
                    self.poly_class[i] = 'existing'
                    changed = True
        counts = collections.Counter(self.poly_class)
        print('filled polygons: ' + ', '.join(f'{k} {v}' for k, v in sorted(counts.items())))

    def _bbox(self, pl: list[Point]) -> tuple[float, float, float, float]:
        xs = [p[0] for p in pl]
        ys = [p[1] for p in pl]
        return (min(xs), min(ys), max(xs), max(ys))

    def _caliper(self, pl: list[Point]) -> tuple[float, float, Point, Point]:
        """Width, length and the two ends of the long axis of a polygon.

        The width is the smallest extent across any edge direction. The
        length is the extent along that edge.
        """
        best = None
        for (ax, ay), (bx, by) in zip(pl, pl[1:] + pl[:1]):
            L = math.hypot(bx - ax, by - ay)
            if L < 1e-9:
                continue
            u = ((bx - ax) / L, (by - ay) / L)
            n = (-u[1], u[0])
            along = [x * u[0] + y * u[1] for x, y in pl]
            across = [x * n[0] + y * n[1] for x, y in pl]
            width = max(across) - min(across)
            if best is None or width < best[0]:
                mid = (max(across) + min(across)) / 2
                best = (width, max(along) - min(along),
                        (min(along) * u[0] + mid * n[0], min(along) * u[1] + mid * n[1]),
                        (max(along) * u[0] + mid * n[0], max(along) * u[1] + mid * n[1]))
        return best if best else (0.0, 0.0, pl[0], pl[0])

    # ------------------------------------------------------------- classes
    def classify(self) -> None:
        """Assigns exterior / glazing / partition / booth / screen / drop to every wall."""
        ends = [self.endpoints(w) for w in self.walls]
        poly_edges = [(a, b) for pl, c in zip(self.polys, self.poly_class) if c != 'drop' for a, b in zip(pl, pl[1:] + pl[:1])]
        touches_core = [False] * len(self.walls)
        free_ends = [2] * len(self.walls)
        links = []
        for i, (p, q) in enumerate(ends):
            for e in (p, q):
                joined = False
                for j, (r, s) in enumerate(ends):
                    if i != j and self._dist(e, r, s) <= self.join_tol:
                        links.append((i, j))
                        joined = True
                if any(self._dist(e, a, b) <= self.join_tol for a, b in poly_edges):
                    touches_core[i] = True
                    joined = True
                if joined:
                    free_ends[i] -= 1
        roots = self._roots(len(self.walls), links)
        comp: dict[int, list[int]] = collections.defaultdict(list)
        for i in range(len(self.walls)):
            comp[roots[i]].append(i)
        def length(idx: list[int]) -> float:
            return sum(self.walls[i][3] - self.walls[i][2] for i in idx)
        main = {i for idx in comp.values() if any(touches_core[i] for i in idx) or length(idx) > 200 * self.pt_per_ft for i in idx}

        all_edges = ends + poly_edges
        def exterior(i: int) -> bool:
            w = self.walls[i]
            p, q = ends[i]
            if round(w[0]) % 90 != 0 and w[3] - w[2] > 15 * self.pt_per_ft:
                return True                 # a long oblique wall is the curtain wall
            if i not in main or w[3] - w[2] < 4 * self.pt_per_ft:
                return False
            mid = ((p[0] + q[0]) / 2, (p[1] + q[1]) / 2)
            n = self._axes(w[0])[1]
            for sign in (1, -1):
                start = (mid[0] + sign * 3 * n[0], mid[1] + sign * 3 * n[1])
                far = (mid[0] + sign * 3000 * n[0], mid[1] + sign * 3000 * n[1])
                if not any(self._crosses(start, far, a, b) for a, b in all_edges if (a, b) != ends[i]):
                    return True
            return False

        # desk rows: one free end, no door at either end and a parallel twin
        # of about the same length within 6 m. A room wall joins at both ends
        # or has a door swing at its open end, up to two leaves from the
        # hinge (a double door).
        exterior_set = {i for i in range(len(self.walls)) if self.walls[i][4] > 0 and exterior(i)}
        def door_at(e: Point) -> bool:
            return any(math.hypot(e[0] - hinge[0], e[1] - hinge[1]) <= 2.2 * r for hinge, r, _ in self.arcs)
        def desk_row(i: int) -> bool:
            w = self.walls[i]
            length = w[3] - w[2]
            if free_ends[i] != 1 or not self.pt_per_m <= length <= 7.5 * self.pt_per_m or touches_core[i]:
                return False
            if any(door_at(e) for e in ends[i]):
                return False
            for j, v in enumerate(self.walls):
                if j == i or free_ends[j] == 0 or round(v[0]) != round(w[0]) or abs(v[1] - w[1]) > 6 * self.pt_per_m:
                    continue
                if 0.7 <= (v[3] - v[2]) / length <= 1.4 and min(v[3], w[3]) - max(v[2], w[2]) > 0.5 * length:
                    return True
            return False
        def bench(i: int) -> bool:
            """A short piece with a free end and no door: a seat, a cabinet or a locker outline.

            It is under 2 ft, or under 1 m with a parallel piece of 0.5 to
            1.5 m within 1.5 m beside it (a row of furniture).
            """
            w = self.walls[i]
            length = w[3] - w[2]
            if free_ends[i] != 1 or length >= self.pt_per_m or any(door_at(e) for e in ends[i]):
                return False
            if length < 2 * self.pt_per_ft:
                return True
            for j, v in enumerate(self.walls):
                if j == i or round(v[0]) != round(w[0]) or not 0 < abs(v[1] - w[1]) <= 1.5 * self.pt_per_m:
                    continue
                if 0.5 * self.pt_per_m <= v[3] - v[2] <= 1.5 * self.pt_per_m and min(v[3], w[3]) > max(v[2], w[2]):
                    return True
            return False
        def corner(idx: list[int]) -> bool:
            return any(not self._parallel(self.walls[i][0], self.walls[j][0]) for i in idx for j in idx)
        def wall_detail(i: int) -> bool:
            """A short piece that hugs the curtain wall: a sill, a mullion or a column face."""
            w = self.walls[i]
            if w[3] - w[2] > 2 * self.pt_per_m:
                return False
            return any(round(self.walls[j][0]) == round(w[0]) and abs(self.walls[j][1] - w[1]) <= 0.3 * self.pt_per_m
                       and min(self.walls[j][3], w[3]) > max(self.walls[j][2], w[2]) for j in exterior_set)
        def oblique(i: int) -> bool:
            """A short piece off the 90 degree grid that is not along the curtain wall: a symbol."""
            w = self.walls[i]
            if w[3] - w[2] >= 2 * self.pt_per_m or round(w[0]) % 90 == 0:
                return False
            return not any(self._parallel(w[0], self.walls[j][0]) for j in exterior_set)
        def pier(i: int) -> bool:
            """A free-standing piece in line with main walls within a door width of both ends."""
            w = self.walls[i]
            for end, sign in ((w[2], -1), (w[3], 1)):
                for j in main:
                    off, lo, hi = self._in_frame(self.walls[j], w[0])
                    if self._parallel(self.walls[j][0], w[0]) and abs(off - w[1]) <= self.join_tol \
                            and 0 < sign * ((lo if sign > 0 else hi) - end) <= 1.2 * self.pt_per_m:
                        break
                else:
                    return False
            return True
        self.records = []
        for i, w in enumerate(self.walls):
            idx = comp[roots[i]]
            if w[4] == 0:
                cls = 'glazing'             # a single line with a jamb block
            elif i in exterior_set:
                cls = 'exterior'
            elif wall_detail(i) or oblique(i) or bench(i):
                cls = 'drop'
            elif i in main and desk_row(i):
                cls = 'screen'
            elif i in main:
                cls = 'partition'
            elif len(idx) >= 4 and all(free_ends[j] == 0 for j in idx) and self._span(idx, ends) >= self.pt_per_m:
                cls = 'booth'               # a closed box of 1 m or more, a smaller one is furniture
            elif corner(idx) and length(idx) >= 4 * self.pt_per_m:
                cls = 'partition'           # a free-standing room corner, desk screens do not turn
            elif pier(i):
                cls = 'partition'           # the wall between two doorways
            elif w[3] - w[2] >= 8 * self.pt_per_ft:
                cls = 'screen'
            else:
                cls = 'drop'
            self.records.append({'class': cls, 'angle': w[0], 'length_m': (w[3] - w[2]) / self.pt_per_m,
                                 'thickness_in': w[4] / self.pt_per_ft * 12, 'p': ends[i][0], 'q': ends[i][1],
                                 'free_ends': free_ends[i]})
        # close the envelope: a free end of a long exterior wall joins the
        # nearest free end of another one within 2 m (the corner columns),
        # or else the nearest corner of a core polygon within 2 m (a pier)
        ext = [i for i, r in enumerate(self.records) if r['class'] == 'exterior' and r['length_m'] >= 2.0]
        free = [(i, e) for i in ext for e in ends[i] if not any(self._dist(e, *ends[j]) <= self.join_tol for j in ext if j != i)]
        corners = [v for pl, c in zip(self.polys, self.poly_class) if c == 'core' for v in pl]
        used: set = set()
        for a_i, a in free:
            if (a_i, a) in used:
                continue
            near = [(math.hypot(a[0] - b[0], a[1] - b[1]), b_i, b) for b_i, b in free if b_i != a_i and (b_i, b) not in used]
            near = [c for c in near if c[0] <= 2 * self.pt_per_m]
            if near:
                dist, b_i, b = min(near)
                used.update({(a_i, a), (b_i, b)})
                self.closures.append((a, b))
                continue
            near = [(math.hypot(a[0] - v[0], a[1] - v[1]), v) for v in corners]
            near = [c for c in near if self.join_tol < c[0] <= 2 * self.pt_per_m]
            if near:
                used.add((a_i, a))
                self.closures.append((a, min(near)[1]))
        counts = collections.Counter(r['class'] for r in self.records)
        print('classes: ' + ', '.join(f'{k} {v}' for k, v in sorted(counts.items())) + f', {len(self.closures)} envelope closures')

    def _span(self, idx: list[int], ends: list[tuple[Point, Point]]) -> float:
        """The shorter side of the bounding box of the walls in idx."""
        xs = [p[0] for i in idx for p in ends[i]]
        ys = [p[1] for i in idx for p in ends[i]]
        return min(max(xs) - min(xs), max(ys) - min(ys))

    def _dist(self, pt: Point, a: Point, b: Point) -> float:
        dx, dy = b[0] - a[0], b[1] - a[1]
        l2 = dx * dx + dy * dy
        t = 0.0 if l2 == 0 else max(0.0, min(1.0, ((pt[0] - a[0]) * dx + (pt[1] - a[1]) * dy) / l2))
        return math.hypot(pt[0] - (a[0] + t * dx), pt[1] - (a[1] + t * dy))

    def _crosses(self, p1: Point, p2: Point, p3: Point, p4: Point) -> bool:
        d = (p2[0] - p1[0]) * (p4[1] - p3[1]) - (p2[1] - p1[1]) * (p4[0] - p3[0])
        if abs(d) < 1e-12:
            return False
        t = ((p3[0] - p1[0]) * (p4[1] - p3[1]) - (p3[1] - p1[1]) * (p4[0] - p3[0])) / d
        u = ((p3[0] - p1[0]) * (p2[1] - p1[1]) - (p3[1] - p1[1]) * (p2[0] - p1[0])) / d
        return 0 < t < 1 and 0 < u < 1

    # -------------------------------------------------------------- output
    def to_m(self, pt: Point) -> Point:
        return ((pt[0] - self.bounds[0]) / self.pt_per_m, (self.bounds[3] - pt[1]) / self.pt_per_m)

    def write(self) -> None:
        """Writes the STL files, walls.json, walls.csv and the overlay."""
        args = self.args
        args.out.mkdir(parents=True, exist_ok=True)
        height = {'exterior': args.ceiling, 'glazing': args.ceiling, 'partition': args.ceiling, 'booth': args.booth_height,
                  'screen': args.screen_height}
        walls = [dict(r, p=self.to_m(r['p']), q=self.to_m(r['q']), height_m=height[r['class']])
                 for r in self.records if r['class'] != 'drop']
        walls += [{'class': 'exterior', 'angle': math.degrees(math.atan2(b[1] - a[1], b[0] - a[0])) % 180.0,
                   'length_m': math.hypot(a[0] - b[0], a[1] - b[1]) / self.pt_per_m,
                   'thickness_in': 0.0, 'p': self.to_m(a), 'q': self.to_m(b), 'height_m': args.ceiling, 'closure': True}
                  for a, b in self.closures]
        core = [[self.to_m(p) for p in pl] for pl, c in zip(self.polys, self.poly_class) if c == 'core']
        existing = [[self.to_m(p) for p in pl] for pl, c in zip(self.polys, self.poly_class) if c == 'existing']
        groups: dict[str, list] = collections.defaultdict(list)
        for w in walls:
            groups[w['class']] += self._quad(w['p'], w['q'], 0.0, w['height_m'])
        if args.outline:
            outline = [tuple(float(v) for v in p.split(',')) for p in args.outline.split()]
        else:
            outline = self._hull([w['p'] for w in walls if w['class'] == 'exterior'] + [w['q'] for w in walls if w['class'] == 'exterior'])
        floor, ceiling = [], []
        for a, b, c in self._ear_clip(outline):
            floor.append(((a[0], a[1], 0.0), (b[0], b[1], 0.0), (c[0], c[1], 0.0)))
            ceiling.append(((a[0], a[1], args.ceiling), (c[0], c[1], args.ceiling), (b[0], b[1], args.ceiling)))
        files = {'walls_partition.stl': groups['partition'] + groups['booth'] + self._prisms(existing, existing + core, args.ceiling),
                 'walls_exterior.stl': groups['exterior'], 'glazing.stl': groups['glazing'],
                 'core.stl': self._prisms(core, existing + core, args.ceiling),
                 'screens.stl': groups['screen'], 'floor.stl': floor, 'ceiling.stl': ceiling}
        for name, tris in files.items():
            self._write_stl(args.out / name, tris)
            print(f'{name}: {len(tris)} triangles')
        with open(args.out / 'walls.json', 'w') as f:
            json.dump({'walls': walls, 'core': core, 'existing': existing, 'outline': outline, 'ceiling': args.ceiling}, f)
        self._write_csv(args.out / 'walls.csv', walls, existing, core, args.ceiling)
        if args.overlay:
            self._overlay(args.overlay)

    def _write_csv(self, path: Path, walls: list[dict], existing: list, core: list, ceiling: float) -> None:
        """Writes one row per wall. A polygon is listed as the segment along its long axis."""
        rows = [(r['class'], r['p'], r['q'], r['thickness_in'], r['height_m'],
                 'closure' if r.get('closure') else {2: 'no', 1: 'one end', 0: 'both ends'}[r.get('free_ends', 0)]) for r in walls]
        for cls, polys in (('existing', existing), ('core', core)):
            for pl in polys:
                width, _, p, q = self._caliper(pl)
                rows.append((cls, p, q, width / 0.0254, ceiling, 'yes'))
        with open(path, 'w', newline='') as f:
            out = csv.writer(f, lineterminator='\n')
            out.writerow(['id', 'class', 'x0', 'y0', 'x1', 'y1', 'length_m', 'thickness_in', 'height_m', 'material', 'joined'])
            for i, (cls, p, q, thickness, height, joined) in enumerate(rows):
                out.writerow([i, cls, *(f'{v:.2f}' for v in (*p, *q)), f'{math.hypot(q[0] - p[0], q[1] - p[1]):.2f}',
                              f'{thickness:.1f}', height, self.MATERIAL[cls], joined])

    def _prisms(self, polys: list, others: list, height: float) -> list:
        """Side faces of the polygons extruded from the floor to height.

        Where an edge lies on an edge of another polygon in ``others`` (two
        abutting fills of one wall) that part gets no face: a ray crosses
        the wall once. Where two fills overlap on the same side of the edge
        (a column on a wall strip) the later one in ``others`` writes the
        face.
        """
        edges = [(k, a, b, 1.0 if self._area(pl) > 0 else -1.0) for k, pl in enumerate(others)
                 if abs(self._area(pl)) > 1e-4 for a, b in zip(pl, pl[1:] + pl[:1])]
        owner = np.array([e[0] for e in edges])
        A = np.array([e[1] for e in edges], dtype=float)
        D = np.array([e[2] for e in edges], dtype=float) - A
        L = np.hypot(D[:, 0], D[:, 1])
        # the edge direction with the fill on its left: fills on opposite
        # sides of a shared edge run opposite ways, overlapping fills the same way
        W = D * (np.array([e[3] for e in edges]) / np.maximum(L, 1e-9))[:, None]
        position = {id(pl): k for k, pl in enumerate(others)}
        tris: list = []
        for pl in polys:
            if abs(self._area(pl)) <= 1e-4:
                continue
            mine = position.get(id(pl), -1)
            turn = 1.0 if self._area(pl) > 0 else -1.0
            for a, b in zip(pl, pl[1:] + pl[:1]):
                length = math.hypot(b[0] - a[0], b[1] - a[1])
                if length < 1e-9:
                    continue
                u = ((b[0] - a[0]) / length, (b[1] - a[1]) / length)
                n = (-u[1], u[0])
                along = (A[:, 0] - a[0]) * u[0] + (A[:, 1] - a[1]) * u[1]
                across = (A[:, 0] - a[0]) * n[0] + (A[:, 1] - a[1]) * n[1]
                heading = (W[:, 0] * u[0] + W[:, 1] * u[1]) * turn
                collinear = (owner != mine) & (L > 1e-9) & (np.abs(across) <= 0.005)
                same = collinear & ((heading <= -0.9999) | ((heading >= 0.9999) & (owner > mine)))
                spans = [(float(min(t, t + s)), float(max(t, t + s))) for t, s in
                         zip(along[same], (D[same, 0] * u[0] + D[same, 1] * u[1]))]
                for lo, hi in self._uncovered(0.0, length, spans, 0.01):
                    tris += self._quad((a[0] + lo * u[0], a[1] + lo * u[1]), (a[0] + hi * u[0], a[1] + hi * u[1]), 0.0, height)
        return tris

    def _quad(self, p: Point, q: Point, z0: float, z1: float) -> list:
        a, b = (p[0], p[1], z0), (q[0], q[1], z0)
        c, d = (q[0], q[1], z1), (p[0], p[1], z1)
        return [(a, b, c), (a, c, d)]

    def _area(self, pl: list[Point]) -> float:
        return 0.5 * sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(pl, pl[1:] + pl[:1]))

    def _hull(self, points: list[Point]) -> list[Point]:
        pts = sorted(set(points))
        def half(seq):
            h: list[Point] = []
            for p in seq:
                while len(h) >= 2 and (h[-1][0] - h[-2][0]) * (p[1] - h[-2][1]) - (h[-1][1] - h[-2][1]) * (p[0] - h[-2][0]) <= 0:
                    h.pop()
                h.append(p)
            return h
        return half(pts)[:-1] + half(reversed(pts))[:-1]

    def _ear_clip(self, poly: list[Point]) -> list[tuple[Point, Point, Point]]:
        """Triangulates a simple counter-clockwise polygon."""
        def cross(o: Point, a: Point, b: Point) -> float:
            return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
        def inside(p: Point, a: Point, b: Point, c: Point) -> bool:
            return cross(a, b, p) >= 0 and cross(b, c, p) >= 0 and cross(c, a, p) >= 0
        pts = list(poly)
        tris = []
        while len(pts) > 3:
            for i in range(len(pts)):
                a, b, c = pts[i - 1], pts[i], pts[(i + 1) % len(pts)]
                if cross(a, b, c) > 0 and not any(inside(p, a, b, c) for p in pts if p not in (a, b, c)):
                    tris.append((a, b, c))
                    del pts[i]
                    break
            else:
                raise SystemExit('the outline is not a simple counter-clockwise polygon')
        tris.append((pts[0], pts[1], pts[2]))
        return tris

    def _write_stl(self, path: Path, tris: list) -> None:
        with open(path, 'wb') as f:
            f.write(path.name.encode().ljust(80, b'\0'))
            f.write(struct.pack('<I', len(tris)))
            for a, b, c in tris:
                n = np.cross(np.subtract(b, a), np.subtract(c, a))
                n = n / (np.linalg.norm(n) or 1.0)
                f.write(struct.pack('<3f', *n))
                for p in (a, b, c):
                    f.write(struct.pack('<3f', *p))
                f.write(struct.pack('<H', 0))

    def _overlay(self, path: Path) -> None:
        """Draws the polygon and wall classes over a faded render of the plan region."""
        from PIL import Image, ImageDraw
        x0, y0 = self.bounds[:2]
        img = Image.blend(self.plan_image, Image.new('RGB', self.plan_image.size, 'white'), 0.6)
        draw = ImageDraw.Draw(img)
        def shift(p: Point) -> Point:
            return ((p[0] - x0) * self.OVERLAY_SCALE, (p[1] - y0) * self.OVERLAY_SCALE)
        fill = {'core': (0, 0, 0), 'existing': (230, 120, 0), 'drop': (190, 190, 190)}
        for pl, c in zip(self.polys, self.poly_class):
            draw.polygon([shift(p) for p in pl], fill=fill[c])
        colour = {'exterior': (0, 140, 200), 'glazing': (0, 170, 170), 'partition': (220, 0, 0), 'booth': (160, 0, 200),
                  'screen': (0, 90, 255), 'drop': (190, 190, 190)}
        for r in self.records:
            draw.line([shift(r['p']), shift(r['q'])], fill=colour[r['class']], width=3)
        for a, b in self.closures:
            draw.line([shift(a), shift(b)], fill=(0, 140, 200), width=3)
        img.save(path)

    def run(self) -> None:
        self.read_pdf()
        self.pair_lines()
        self.drop_non_walls()
        self.drop_stacked()
        self.sort_polys()
        self.classify()
        self.write()


if __name__ == '__main__':
    PlanToStl().run()
