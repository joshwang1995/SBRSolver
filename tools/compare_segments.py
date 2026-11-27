#!/usr/bin/env python3
"""Compare the traced segments of two sbrsolver runs, launch by launch.

    compare_segments.py REF.csv CAND.csv [--max-differing N]
                        [--max-differing-fraction F] [--rtol R] [--allow FILE]

Both files are Segments.csv files (written with --segments-csv or
[outputs] segments_csv = true). A launch matches when both runs list the
same sequence of (Index, Parent, Kind, TriangleID, CoplanarID, MaterialID)
and every Length and UnfoldedStart agrees to --rtol (default 1e-9) relative
to the larger value. The check passes when the launch counts are equal and
at most --max-differing launches (default 0), or at most the fraction
--max-differing-fraction of all launches, differ; launches listed in
--allow (one id per line, a reason after '#') do not count.

The verdict line is the interface for scripts:

    RESULT segments: PASS (launches=L, differing=D, allowed=A, max_rel_err=E)

Exit status: 0 on pass, 1 on fail, 2 on a usage error, 77 when CAND.csv
does not exist (the candidate run was skipped).
"""
import argparse
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

SKIP = 77
HEADER = ["LaunchID", "Index", "Parent", "Kind", "TriangleID", "CoplanarID", "MaterialID",
          "Length", "UnfoldedStart", "Ox", "Oy", "Oz", "Dx", "Dy", "Dz"]

Key = Tuple[int, int, int, int, int, int]
Row = Tuple[Key, float, float]


def read_launches(path: Path) -> Dict[int, List[Row]]:
    """Returns the rows of every launch, in file order."""
    launches: Dict[int, List[Row]] = {}
    with open(path) as fh:
        header = next(fh).strip().split(",")
        if header != HEADER:
            raise SystemExit(f"{path}: header is {header}, expected {HEADER}")
        for line in fh:
            if not line.strip():
                continue
            cells = line.strip().split(",")
            launch = int(cells[0])
            key = tuple(int(c) for c in cells[1:7])
            launches.setdefault(launch, []).append((key, float(cells[7]), float(cells[8])))
    return launches


def read_allow(path: Path) -> Set[int]:
    allowed = set()
    with open(path) as fh:
        for line in fh:
            text = line.split("#", 1)[0].strip()
            if text:
                allowed.add(int(text))
    return allowed


def relative_error(a: float, b: float) -> float:
    scale = max(abs(a), abs(b))
    return abs(a - b) / scale if scale > 0.0 else 0.0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reference", type=Path, help="Segments.csv of the reference run")
    parser.add_argument("candidate", type=Path, help="Segments.csv of the run under test")
    parser.add_argument("--max-differing", type=int, default=0, help="launches that may differ")
    parser.add_argument("--max-differing-fraction", type=float, default=None,
                        help="fraction of all launches that may differ, in place of --max-differing")
    parser.add_argument("--rtol", type=float, default=1e-9, help="relative tolerance on Length and UnfoldedStart")
    parser.add_argument("--allow", type=Path, default=None, help="file of launch ids that may differ, one per line")
    args = parser.parse_args()

    if not args.candidate.exists():
        print(f"SKIP: {args.candidate} not found (the candidate run was skipped)")
        return SKIP
    reference = read_launches(args.reference)
    candidate = read_launches(args.candidate)
    if len(reference) != len(candidate) or set(reference) != set(candidate):
        print(f"RESULT segments: FAIL (launches={len(reference)} vs {len(candidate)})")
        return 1
    allowed = read_allow(args.allow) if args.allow else set()

    differing: List[Tuple[int, str]] = []
    allowed_hits = 0
    max_rel_err = 0.0
    for launch in sorted(reference):
        ref_rows = reference[launch]
        cand_rows = candidate[launch]
        reason = ""
        if len(ref_rows) != len(cand_rows):
            reason = f"{len(ref_rows)} vs {len(cand_rows)} segments"
        else:
            for (rkey, rlength, rstart), (ckey, clength, cstart) in zip(ref_rows, cand_rows):
                if rkey != ckey:
                    reason = f"segment {rkey[0]}: {rkey} vs {ckey}"
                    break
                err = max(relative_error(rlength, clength), relative_error(rstart, cstart))
                max_rel_err = max(max_rel_err, err)
                if err > args.rtol:
                    reason = f"segment {rkey[0]}: length {rlength} vs {clength}, start {rstart} vs {cstart}"
                    break
        if reason:
            if launch in allowed:
                allowed_hits += 1
            else:
                differing.append((launch, reason))

    for launch, reason in differing[:10]:
        print(f"launch {launch}: {reason}")
    limit = args.max_differing
    if args.max_differing_fraction is not None:
        limit = int(args.max_differing_fraction * len(reference))
    ok = len(differing) <= limit
    verdict = "PASS" if ok else "FAIL"
    print(f"RESULT segments: {verdict} (launches={len(reference)}, differing={len(differing)}, allowed={allowed_hits}, "
          f"max_rel_err={max_rel_err:.3g})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
