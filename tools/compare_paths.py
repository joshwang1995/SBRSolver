#!/usr/bin/env python3
"""Compare the paths and fields of two sbrsolver runs, receiver by receiver.

    compare_paths.py REF_DIR CAND_DIR [--max-differing N] [--field-rtol R]
                     [--max-db D] [--allow FILE]

Both directories hold a Paths.csv (written with --paths-csv or
[outputs] paths_csv = true) and an ElectricField.csv. A path is identified by
its key: the transmitter and the sequence of surfaces and interaction types,
the Interactions column. The check passes when:

* both runs have the same receivers at the same positions;
* at most --max-differing receivers (default 0) have different key sets,
  not counting the receivers listed in --allow (one id per line, a reason
  after '#');
* at every receiver with equal key sets, each field component agrees to
  --field-rtol (default 1e-9) relative to the larger magnitude, with an
  absolute floor of 1e-15 V/m;
* at every differing receiver, |E| differs by at most --max-db dB when that
  bound is given.

The verdict line is the interface for scripts:

    RESULT paths: PASS (receivers=R, differing=D, allowed=A,
                        max_field_rel_err=E, max_db_diff=X)

Exit status: 0 on pass, 1 on fail, 2 on a usage error, 77 when CAND_DIR has
no ElectricField.csv (the candidate run was skipped).
"""
import argparse
import math
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

SKIP = 77
HEADER = ["ReceiverID", "LaunchID", "TransmitterID", "Interactions", "Ex", "Ey", "Ez"]

Field = Tuple[complex, complex, complex]
Key = Tuple[int, str]


def parse_complex(text: str) -> complex:
    text = text.strip()
    if text.endswith("i"):
        text = text[:-1] + "j"
    return complex(text)


def read_fields(path: Path) -> List[Tuple[Tuple[float, ...], Field]]:
    """Returns (position, (Ex, Ey, Ez)) per receiver, in file order."""
    rows = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            if not line.strip():
                continue
            cells = line.strip().split(",")
            position = tuple(float(v) for v in cells[:3])
            field = tuple(parse_complex(v) for v in cells[3:6])
            rows.append((position, field))
    return rows


def read_keys(path: Path, receiver_count: int) -> List[Set[Key]]:
    """Returns the set of path keys of every receiver."""
    keys: List[Set[Key]] = [set() for _ in range(receiver_count)]
    with open(path) as fh:
        header = next(fh).strip().split(",")
        if header != HEADER:
            raise SystemExit(f"{path}: header is {header}, expected {HEADER}")
        for line in fh:
            if not line.strip():
                continue
            cells = line.strip().split(",")
            receiver = int(cells[0])
            if receiver < 0 or receiver >= receiver_count:
                raise SystemExit(f"{path}: receiver {receiver} is outside 0..{receiver_count - 1}")
            keys[receiver].add((int(cells[2]), cells[3]))
    return keys


def read_allow(path: Path) -> Set[int]:
    allowed = set()
    with open(path) as fh:
        for line in fh:
            text = line.split("#", 1)[0].strip()
            if text:
                allowed.add(int(text))
    return allowed


def magnitude(field: Field) -> float:
    return math.sqrt(sum(abs(c) ** 2 for c in field))


def db_difference(reference: Field, candidate: Field) -> float:
    ref_mag = magnitude(reference)
    cand_mag = magnitude(candidate)
    if ref_mag > 0.0 and cand_mag > 0.0:
        return abs(20.0 * math.log10(cand_mag / ref_mag))
    if ref_mag == 0.0 and cand_mag == 0.0:
        return 0.0
    return math.inf


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ref_dir", type=Path, help="output directory of the reference run")
    parser.add_argument("cand_dir", type=Path, help="output directory of the run under test")
    parser.add_argument("--max-differing", type=int, default=0, help="receivers that may have different key sets")
    parser.add_argument("--field-rtol", type=float, default=1e-9, help="relative tolerance per field component on equal receivers")
    parser.add_argument("--max-db", type=float, default=None, help="largest |E| difference in dB at a differing receiver")
    parser.add_argument("--allow", type=Path, default=None, help="file of receiver ids that may differ, one per line")
    args = parser.parse_args()

    cand_field_path = args.cand_dir / "ElectricField.csv"
    if not cand_field_path.exists():
        print(f"SKIP: {cand_field_path} not found (the candidate run was skipped)")
        return SKIP
    ref_fields = read_fields(args.ref_dir / "ElectricField.csv")
    cand_fields = read_fields(cand_field_path)
    if len(ref_fields) != len(cand_fields):
        print(f"RESULT paths: FAIL (receivers={len(ref_fields)} vs {len(cand_fields)})")
        return 1
    receivers = len(ref_fields)
    for i, ((rpos, _), (cpos, _)) in enumerate(zip(ref_fields, cand_fields)):
        if any(not math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-9) for a, b in zip(rpos, cpos)):
            print(f"RESULT paths: FAIL (receiver {i} is at {rpos} in one run and {cpos} in the other)")
            return 1

    ref_keys = read_keys(args.ref_dir / "Paths.csv", receivers)
    cand_keys = read_keys(args.cand_dir / "Paths.csv", receivers)
    allowed = read_allow(args.allow) if args.allow else set()

    differing: List[int] = []
    allowed_hits = 0
    max_rel_err = 0.0
    max_db = 0.0
    field_failures = 0
    for i in range(receivers):
        rfield = ref_fields[i][1]
        cfield = cand_fields[i][1]
        if ref_keys[i] != cand_keys[i]:
            if i in allowed:
                allowed_hits += 1
            else:
                differing.append(i)
            max_db = max(max_db, db_difference(rfield, cfield))
            continue
        for comp, (a, b) in zip("xyz", zip(rfield, cfield)):
            err = abs(a - b)
            scale = max(abs(a), abs(b))
            rel = err / scale if scale > 0.0 else 0.0
            max_rel_err = max(max_rel_err, rel)
            if err > 1e-15 and rel > args.field_rtol:
                field_failures += 1
                if field_failures <= 10:
                    print(f"receiver {i}: E{comp} differs: {a} vs {b} (relative {rel:.3g})")

    for i in differing[:10]:
        only_ref = sorted(ref_keys[i] - cand_keys[i])
        only_cand = sorted(cand_keys[i] - ref_keys[i])
        print(f"receiver {i}: only in reference: {only_ref}; only in candidate: {only_cand}")

    ok = len(differing) <= args.max_differing and field_failures == 0
    if args.max_db is not None and max_db > args.max_db:
        ok = False
    db_text = "inf" if math.isinf(max_db) else f"{max_db:.3g}"
    verdict = "PASS" if ok else "FAIL"
    print(f"RESULT paths: {verdict} (receivers={receivers}, differing={len(differing)}, allowed={allowed_hits}, "
          f"max_field_rel_err={max_rel_err:.3g}, max_db_diff={db_text})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
