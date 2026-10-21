#!/usr/bin/env python3
"""Sanity check of an SBRSolver field file against its receiver file.

    check_fields.py ElectricField.csv RX_file.dat [--all-nonzero] [--max-zero N]

Passes (exit 0) when the CSV holds exactly one row per receiver of the
receiver file, at the receiver's position, with finite complex components.
``--all-nonzero`` also requires a non-zero field at every receiver (scenes
where every receiver has a line of sight); ``--max-zero N`` allows at most N
receivers without any field. Only the standard library is used, so the
check runs wherever ctest finds a Python interpreter.
"""
import math
import sys


def read_receivers(path):
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip()]
    count = int(lines[1])
    points = [tuple(float(v) for v in l.split()[:3]) for l in lines[2:2 + count]]
    if len(points) != count:
        raise SystemExit('%s: header announces %d receivers, %d found' % (path, count, len(points)))
    return points


def parse_complex(text):
    # SBRSolver writes "+1.23-0.45i" (no space, 'i' suffix)
    return complex(text.strip().replace('i', 'j'))


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 2
    csv_path, rx_path = argv[1], argv[2]
    all_nonzero = '--all-nonzero' in argv
    max_zero = None
    if '--max-zero' in argv:
        max_zero = int(argv[argv.index('--max-zero') + 1])
    receivers = read_receivers(rx_path)
    with open(csv_path) as f:
        rows = [l.strip().split(',') for l in f if l.strip()]
    header, rows = rows[0], rows[1:]
    if len(rows) != len(receivers):
        print('FAIL: %d rows for %d receivers' % (len(rows), len(receivers)))
        return 1
    errors = 0
    zeros = 0
    for i, (row, rx) in enumerate(zip(rows, receivers)):
        if len(row) < 6:
            print('FAIL: row %d has %d columns' % (i, len(row)))
            return 1
        pos = tuple(float(v) for v in row[:3])
        if max(abs(a - b) for a, b in zip(pos, rx)) > 1e-6 * max(1.0, max(abs(c) for c in rx)):
            print('FAIL: row %d is at %s, receiver %d is at %s' % (i, pos, i, rx))
            errors += 1
            if errors > 10:
                return 1
            continue
        field = [parse_complex(c) for c in row[3:6]]
        magnitude = math.sqrt(sum(abs(c) ** 2 for c in field))
        if not math.isfinite(magnitude):
            print('FAIL: non-finite field at receiver %d %s' % (i, rx))
            errors += 1
        elif magnitude == 0.0:
            zeros += 1
    if errors:
        return 1
    print('%d receivers, %d without field' % (len(receivers), zeros))
    if all_nonzero and zeros:
        print('FAIL: --all-nonzero but %d receivers have no field' % zeros)
        return 1
    if max_zero is not None and zeros > max_zero:
        print('FAIL: %d receivers have no field, at most %d allowed' % (zeros, max_zero))
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
