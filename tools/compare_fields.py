#!/usr/bin/env python3
"""Compare two ElectricField.csv files produced by sbrsolver.

Each row holds X,Y,Z followed by the complex Ex,Ey,Ez components written
with 12 significant digits (e.g. "-2.52878846646+0.000900257848066i"). The
comparison is tolerance based so that results produced on different platforms
or compilers (which differ in the last bits of libm functions) still count as
equal, while real regressions are caught: the default relative tolerance of
2e-5 is far above the printing granularity and absorbs that noise.
Non-finite values (nan, inf) in either file always count as a mismatch.

Usage: compare_fields.py <reference.csv> <candidate.csv> [--rtol R] [--atol A]
Exit status is 0 when the files agree, 1 otherwise.
"""
import argparse
import csv
import math
import re
import sys


def parse_complex(text):
    text = text.strip()
    if text.endswith("i"):
        text = text[:-1] + "j"
    # MSVC spells NaN as "nan(ind)" or "-nan(ind)"
    text = re.sub(r"nan\([^)]*\)", "nan", text)
    return complex(text)


def is_finite(value):
    return math.isfinite(value.real) and math.isfinite(value.imag)


def read_rows(path):
    rows = []
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        for line in reader:
            if not line:
                continue
            pos = tuple(float(v) for v in line[:3])
            field = tuple(parse_complex(v) for v in line[3:6])
            rows.append((pos, field))
    return header, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reference")
    parser.add_argument("candidate")
    parser.add_argument("--rtol", type=float, default=2e-5,
                        help="relative tolerance per component (default 2e-5, platform and compiler noise)")
    parser.add_argument("--atol", type=float, default=1e-9,
                        help="absolute tolerance, scaled by the largest field magnitude in the reference")
    args = parser.parse_args()

    ref_header, ref = read_rows(args.reference)
    cand_header, cand = read_rows(args.candidate)

    if ref_header != cand_header:
        print(f"header mismatch: {ref_header} vs {cand_header}")
        return 1
    if len(ref) != len(cand):
        print(f"row count mismatch: {len(ref)} reference rows vs {len(cand)} candidate rows")
        return 1

    scale = max((abs(c) for _, field in ref for c in field), default=1.0)
    atol = args.atol * scale
    worst = 0.0
    bad = 0
    for idx, ((rpos, rfield), (cpos, cfield)) in enumerate(zip(ref, cand), start=2):
        if any(not math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-9) for a, b in zip(rpos, cpos)):
            print(f"line {idx}: receiver position differs: {rpos} vs {cpos}")
            bad += 1
            continue
        for comp, (a, b) in zip("xyz", zip(rfield, cfield)):
            if not (is_finite(a) and is_finite(b)):
                bad += 1
                if bad <= 10:
                    print(f"line {idx}: E{comp} is not finite: {a} vs {b}")
                continue
            err = abs(a - b)
            tol = atol + args.rtol * max(abs(a), abs(b))
            worst = max(worst, err / tol if tol else err)
            if err > tol:
                bad += 1
                if bad <= 10:
                    print(f"line {idx}: E{comp} differs: {a} vs {b} (|diff|={err:.3e}, tol={tol:.3e})")

    if bad:
        print(f"FAILED: {bad} mismatching values (worst error / tolerance = {worst:.3g})")
        return 1
    print(f"OK: {len(ref)} receivers match (worst error / tolerance = {worst:.3g})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
