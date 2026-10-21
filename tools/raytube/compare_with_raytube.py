#!/usr/bin/env python3
"""Compare SBRSolver's ElectricField.csv with RayTubeTracing's.

SBRSolver writes ``X,Y,Z,Ex,Ey,Ez`` with complex numbers such as
``-2.52878846646+0.000900257848066i``. RayTubeTracing.exe writes only
``E_x, E_y, E_z`` (``+0.0097240028 -0.003705218I`` in ElectricField.csv,
``+0.0097240028-0.003705218i`` in ElectricFieldMatlab.csv) and its receiver
positions are taken from the ReceivedPowers.csv next to it. Either layout is
accepted for either argument.

A reference without any field or path at the compared receivers fails the
comparison (nothing would be checked) unless --allow-empty is given.

Receivers are matched by row order and their positions must agree. Each
field component must satisfy |a - b| <= atol * max|E| + rtol * |b| where the
maximum is taken over the reference file. RayTubeTracing prints 8
significant digits, so the default rtol of 1e-7 (one unit in the last
printed digit) demands agreement to every digit it prints; the scaled atol
covers components that are zero up to rounding (for example the cross-polar
component over a flat ground).

Usage: compare_with_raytube.py <sbrsolver.csv> <raytube.csv> [--rtol R] [--atol A]
       [--worst N] [--paths sbr_RayPath.vtk raytube_ReceivedPowers.csv]
Exit status 0 when every receiver agrees, 1 otherwise.
"""
import argparse
import math
import os
import re
import sys

COMPLEX = re.compile(r"^\s*([+-]?[^+-]+(?:[eE][+-]?\d+)?)\s*([+-]\s*[^+-]+(?:[eE][+-]?\d+)?)\s*[iIjJ]\s*$")


def parse_complex(text):
    text = text.strip().replace("nan(ind)", "nan")
    m = COMPLEX.match(text)
    if not m:
        raise ValueError("cannot parse complex number '%s'" % text)
    return complex(float(m.group(1)), float(m.group(2).replace(" ", "")))


def read_field(path):
    """Rows of (position or None, (Ex, Ey, Ez)).

    RayTubeTracing.exe writes only the three field columns to
    ElectricField.csv; its receiver positions then come from the
    ReceivedPowers.csv next to it when that file exists.
    """
    rows = []
    with open(path) as fh:
        header = fh.readline()
        for number, line in enumerate(fh, 2):
            line = line.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) == 6:
                pos = tuple(float(p) for p in parts[:3])
                field = tuple(parse_complex(p) for p in parts[3:])
            elif len(parts) == 3:
                pos = None
                field = tuple(parse_complex(p) for p in parts)
            else:
                raise ValueError("%s:%d: expected 3 or 6 columns, got %d" % (path, number, len(parts)))
            rows.append((pos, field))
    if rows and rows[0][0] is None:
        powers = os.path.join(os.path.dirname(os.path.abspath(path)), "ReceivedPowers.csv")
        if os.path.exists(powers):
            positions = []
            with open(powers) as fh:
                fh.readline()
                for line in fh:
                    parts = [p.strip() for p in line.split(",")]
                    if len(parts) >= 3:
                        positions.append(tuple(float(p) for p in parts[:3]))
            if len(positions) == len(rows):
                rows = [(positions[i], rows[i][1]) for i in range(len(rows))]
    return header.strip(), rows


def read_raytube_path_counts(path):
    """nDuplicateRemoved column of ReceivedPowers.csv, per receiver."""
    counts = []
    with open(path) as fh:
        fh.readline()
        for line in fh:
            parts = [p.strip() for p in line.split(",")]
            if len(parts) >= 9:
                counts.append(int(parts[8]))
    return counts


def read_sbr_path_counts(path, receivers):
    """Number of polylines per receiver in SBRSolver's RayPath.vtk."""
    counts = [0] * receivers
    with open(path) as fh:
        lines = fh.read().split("\n")
    for i, line in enumerate(lines):
        tokens = line.split()
        if len(tokens) == 4 and tokens[0] == "ReceiverID":
            n = int(tokens[2])
            values = []
            j = i + 1
            while len(values) < n and j < len(lines):
                values.extend(int(v) for v in lines[j].split())
                j += 1
            for r in values:
                if 0 <= r < receivers:
                    counts[r] += 1
            break
    return counts


def fmt(c):
    return "%+.8g%+.8gi" % (c.real, c.imag)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sbrsolver")
    parser.add_argument("raytube")
    parser.add_argument("--rtol", type=float, default=1e-7)
    parser.add_argument("--atol", type=float, default=1e-9, help="scaled by the largest |E| component of the reference")
    parser.add_argument("--worst", type=int, default=5, help="how many of the worst receivers to print")
    parser.add_argument("--paths", nargs=2, metavar=("SBR_RAYPATH_VTK", "RAYTUBE_RECEIVEDPOWERS_CSV"),
                        help="also compare the number of paths found per receiver")
    parser.add_argument("--receivers", metavar="FIRST-LAST",
                        help="compare only this range of receiver indices (inclusive), e.g. 0-1000")
    parser.add_argument("--allow-empty", action="store_true",
                        help="pass even when the reference has no field at any compared receiver")
    args = parser.parse_args()

    _, sbr = read_field(args.sbrsolver)
    _, ref = read_field(args.raytube)
    ok = True
    if len(sbr) != len(ref):
        print("FAIL: %d receivers in %s but %d in %s" % (len(sbr), args.sbrsolver, len(ref), args.raytube))
        return 1
    first, last = 0, len(ref) - 1
    if args.receivers:
        m = re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*", args.receivers)
        if not m or int(m.group(1)) > int(m.group(2)) or int(m.group(2)) >= len(ref):
            print("FAIL: --receivers must be FIRST-LAST within 0-%d" % (len(ref) - 1))
            return 1
        first, last = int(m.group(1)), int(m.group(2))
        print("comparing receivers %d to %d of %d" % (first, last, len(ref)))
    sbr = sbr[first:last + 1]
    ref = ref[first:last + 1]

    scale = max((abs(c) for _, f in ref for c in f), default=0.0)
    if not math.isfinite(scale):
        print("FAIL: the reference contains non-finite values")
        return 1
    atol = args.atol * scale
    if scale == 0.0 and not args.allow_empty:
        print("FAIL: the reference field is zero at every compared receiver, so nothing would be compared "
              "(move the receivers into the transmitter's coverage, or pass --allow-empty)")
        return 1

    worst = []
    bad = 0
    max_rel = 0.0
    max_abs = 0.0
    for i, ((ps, fs), (pr, fr)) in enumerate(zip(sbr, ref), first):
        if ps is None or pr is None:
            pr = ps if pr is None else pr
        elif any(abs(a - b) > 1e-6 * max(1.0, abs(b)) for a, b in zip(ps, pr)):
            print("FAIL: receiver %d is at %s in %s but at %s in %s" % (i, ps, args.sbrsolver, pr, args.raytube))
            return 1
        err = 0.0
        rel = 0.0
        for a, b in zip(fs, fr):
            if not (math.isfinite(a.real) and math.isfinite(a.imag) and math.isfinite(b.real) and math.isfinite(b.imag)):
                err = float("inf")
                rel = float("inf")
                continue
            d = abs(a - b)
            err = max(err, d)
            rel = max(rel, d / (atol + args.rtol * abs(b)) if (atol + args.rtol * abs(b)) > 0 else (0.0 if d == 0 else float("inf")))
        max_abs = max(max_abs, err)
        max_rel = max(max_rel, err / scale if scale > 0 else err)
        if rel > 1.0:
            bad += 1
        worst.append((rel, i, pr, fs, fr))
    worst.sort(key=lambda w: -w[0])

    print("%d receivers, max |dE| = %.3e (%.3e of the largest reference component %.6g)" % (len(ref), max_abs, max_rel, scale))
    if args.paths:
        if not os.path.exists(args.paths[0]):
            print("note: %s not found, path counts not compared (set [outputs] ray_paths = true)" % args.paths[0])
        else:
            total = first + len(ref) if not args.receivers else None
            sbr_all = read_sbr_path_counts(args.paths[0], max(total or 0, last + 1))
            ref_all = read_raytube_path_counts(args.paths[1])
            if len(ref_all) <= last:
                print("FAIL: %s lists %d receivers" % (args.paths[1], len(ref_all)))
                ok = False
            else:
                sbr_counts = sbr_all[first:last + 1]
                ref_counts = ref_all[first:last + 1]
                if sum(ref_counts) == 0 and not args.allow_empty:
                    print("FAIL: the reference found no path at any compared receiver")
                    ok = False
                differing = [i for i in range(len(ref_counts)) if sbr_counts[i] != ref_counts[i]]
                print("paths per receiver: sbrsolver %d, raytube %d, receivers with a different count: %d%s"
                      % (sum(sbr_counts), sum(ref_counts), len(differing),
                         (" (first: " + ", ".join("%d: %d vs %d" % (first + i, sbr_counts[i], ref_counts[i]) for i in differing[:8]) + ")") if differing else ""))
                if differing:
                    ok = False
    if bad:
        ok = False
        print("FAIL: %d receivers exceed the tolerance (rtol %g, atol %g)" % (bad, args.rtol, atol))
    for rel, i, pos, fs, fr in worst[:args.worst]:
        where = "at (%.6g, %.6g, %.6g)" % pos if pos is not None else ""
        print("  receiver %d %s: %.2fx tolerance" % (i, where, rel))
        for k, name in enumerate("xyz"):
            print("    E%s sbrsolver %-28s raytube %s" % (name, fmt(fs[k]), fmt(fr[k])))
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
