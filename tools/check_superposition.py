#!/usr/bin/env python3
"""Check that a run with several transmitters is the sum of the single-transmitter runs.

    check_superposition.py SBRSOLVER CONFIG.ini FULL_DIR WORK_DIR [--rtol R] [--timeout S]

``FULL_DIR`` holds the results of ``sbrsolver CONFIG.ini FULL_DIR``, where
CONFIG.ini has two or more ``[transmitter]`` sections. The script writes one
configuration per transmitter into WORK_DIR: that section alone, with the
first section's values for the keys it leaves out (never its name), and every
file name made absolute. It runs SBRSOLVER on each into ``WORK_DIR/tx<index>``
and checks:

* ``ElectricField.csv``: every run has the same receivers, and at each one the
  field of the full run equals the sum of the single-transmitter fields,
  within R (default 1e-9) times the largest single field at that receiver.
  Every transmitter must reach at least one receiver.
* ``RayPath.vtk``: every path of the full run carries a ``TransmitterID``, and
  for every transmitter the full run holds as many paths with its id as the
  single run found, at least one.

Every check prints a ``RESULT`` line. The exit status is 0 only when every
check passes. The generated configurations and results stay in WORK_DIR.
"""
import argparse
import os
import re
import subprocess
import sys

PATH_KEYS = ('stl_file', 'file', 'pattern_file')
ERROR_MARKERS = re.compile(r'^(Error\b|Assertion)', re.M)


def read_sections(path):
    """The configuration as a list of (name, [(key, value)]) in file order.

    Section names and keys are lower case, as the solver reads them.
    """
    sections = []
    with open(path, encoding='utf-8-sig') as f:
        for raw in f:
            line = raw.strip()
            # Like the solver: a whole line is a comment when it starts with
            # '#' or ';', and an inline comment starts with '#' after
            # whitespace, so ';' stays usable as the matrix row separator.
            if not line or line[0] in '#;':
                continue
            text = re.split(r'\s#', line, maxsplit=1)[0].strip()
            if not text:
                continue
            if text.startswith('['):
                sections.append((text.strip('[]').strip().lower(), []))
                continue
            if '=' not in text or not sections:
                raise SystemExit('%s: cannot parse line %r' % (path, line))
            key, value = (s.strip() for s in text.split('=', 1))
            sections[-1][1].append((key.lower(), value))
    return sections


def transmitter_names(transmitters):
    """The name of every transmitter as the solver assigns it.

    A section without a name gets the first unused of tx0, tx1, ...
    """
    names = []
    for entries in transmitters:
        name = dict(entries).get('name')
        if name is None:
            n = 0
            while 'tx%d' % n in names:
                n += 1
            name = 'tx%d' % n
        names.append(name)
    return names


def absolute(key, value, base):
    if key in PATH_KEYS and not os.path.isabs(value):
        return os.path.normpath(os.path.join(base, value))
    return value


def write_single(sections, index, name, base, path):
    """Writes the configuration with only transmitter `index`, named `name`."""
    transmitters = [entries for section, entries in sections if section == 'transmitter']
    merged = {key: value for key, value in transmitters[0] if key != 'name'}
    merged.update(dict(transmitters[index]))
    merged['name'] = name
    with open(path, 'w', encoding='utf-8') as f:
        f.write('# transmitter %d (%s) of %s alone, written by check_superposition.py\n'
                % (index, name, os.path.basename(base)))
        written = False
        for section, entries in sections:
            if section == 'transmitter':
                if written:
                    continue
                written = True
                entries = list(merged.items())
            f.write('\n[%s]\n' % section)
            for key, value in entries:
                f.write('%s = %s\n' % (key, absolute(key, value, os.path.dirname(base))))


def require_file(path):
    if not os.path.isfile(path):
        raise SystemExit('FAIL: %s: file not found (the solver writes it when [outputs] keeps the default true)' % path)


def parse_complex(text):
    text = text.strip()
    if text.endswith('i'):
        text = text[:-1] + 'j'
    return complex(text)


def read_fields(path):
    """[(x, y, z), (Ex, Ey, Ez)] per receiver."""
    require_file(path)
    rows = []
    with open(path) as f:
        header = f.readline()
        if not header.startswith('X,Y,Z'):
            raise SystemExit('%s: unexpected header %r' % (path, header.strip()))
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split(',')
            if len(parts) != 6:
                raise SystemExit('%s: expected 6 columns, got %r' % (path, line.strip()))
            rows.append((tuple(float(v) for v in parts[:3]), tuple(parse_complex(v) for v in parts[3:])))
    if not rows:
        raise SystemExit('%s: no receivers' % path)
    return rows


def read_path_tags(path):
    """The number of paths (LINES) and the list of TransmitterID values, if present."""
    require_file(path)
    lines = None
    tags = []
    with open(path) as f:
        for line in f:
            if line.startswith('LINES '):
                lines = int(line.split()[1])
            elif line.startswith('TransmitterID '):
                count = int(line.split()[2])
                tags = [int(f.readline()) for _ in range(count)]
    if lines is None:
        raise SystemExit('%s: no LINES block' % path)
    return lines, tags


def run_solver(solver, config, out_dir, timeout):
    try:
        result = subprocess.run([solver, config, out_dir], capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        raise SystemExit('FAIL: %s %s did not finish within %g s' % (solver, config, timeout))
    output = result.stdout + result.stderr
    if result.returncode != 0:
        sys.stderr.write(output[-4000:])
        raise SystemExit('FAIL: %s %s exited with %d' % (solver, config, result.returncode))
    marker = ERROR_MARKERS.search(output)
    if marker:
        sys.stderr.write(output[-4000:])
        raise SystemExit('FAIL: %s %s reported %r' % (solver, config, marker.group(0)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('solver')
    ap.add_argument('config')
    ap.add_argument('full_dir')
    ap.add_argument('work_dir')
    ap.add_argument('--rtol', type=float, default=1e-9, help='relative tolerance of the field sum (default 1e-9)')
    ap.add_argument('--timeout', type=float, default=1800, help='seconds allowed per solver run (default 1800)')
    args = ap.parse_args()

    config = os.path.abspath(args.config)
    print('check_superposition: solver=%s config=%s full=%s work=%s rtol=%g timeout=%gs'
          % (args.solver, config, args.full_dir, args.work_dir, args.rtol, args.timeout))
    sections = read_sections(config)
    transmitters = [entries for section, entries in sections if section == 'transmitter']
    if len(transmitters) < 2:
        raise SystemExit('FAIL: %s has %d [transmitter] section(s); the check needs at least two'
                         % (config, len(transmitters)))
    names = transmitter_names(transmitters)
    os.makedirs(args.work_dir, exist_ok=True)

    singles = []
    for index, name in enumerate(names):
        ini = os.path.join(args.work_dir, 'tx%d.ini' % index)
        write_single(sections, index, name, config, ini)
        out = os.path.join(args.work_dir, 'tx%d' % index)
        run_solver(args.solver, ini, out, args.timeout)
        singles.append((name, read_fields(os.path.join(out, 'ElectricField.csv')),
                        read_path_tags(os.path.join(out, 'RayPath.vtk'))[0]))
    full = read_fields(os.path.join(args.full_dir, 'ElectricField.csv'))
    for name, rows, _ in singles:
        if len(rows) != len(full):
            raise SystemExit('FAIL: %s has %d receivers, the full run %d' % (name, len(rows), len(full)))
    failures = 0

    # 1. the field of the full run is the sum of the single fields, and every
    # transmitter reaches at least one receiver
    worst = 0.0
    nonzero = 0
    reached = [0] * len(singles)
    for k, (position, field) in enumerate(full):
        total = [0j, 0j, 0j]
        scale = 0.0
        for i, (name, rows, _) in enumerate(singles):
            if rows[k][0] != position:
                raise SystemExit('FAIL: receiver %d of %s is at %s, in the full run at %s'
                                 % (k, name, rows[k][0], position))
            total = [t + f for t, f in zip(total, rows[k][1])]
            largest = max(abs(c) for c in rows[k][1])
            if largest > 0.0:
                reached[i] += 1
            scale = max(scale, largest)
        error = max(abs(f - t) for f, t in zip(field, total))
        if scale > 0.0:
            nonzero += 1
            worst = max(worst, error / scale)
        elif error > 0.0:
            worst = float('inf')
    detail = ', '.join('%s=%d' % (singles[i][0], reached[i]) for i in range(len(singles)))
    unreached = [singles[i][0] for i in range(len(singles)) if reached[i] == 0]
    if unreached:
        print('RESULT superposition: FAIL (no field at any receiver from transmitter %s; receivers reached: %s)'
              % (', '.join(unreached), detail))
        failures += 1
    elif worst > args.rtol:
        print('RESULT superposition: FAIL (receivers=%d, nonzero=%d, max_rel_err=%.3g, rtol=%.3g, reached: %s)'
              % (len(full), nonzero, worst, args.rtol, detail))
        failures += 1
    else:
        print('RESULT superposition: PASS (receivers=%d, nonzero=%d, max_rel_err=%.3g, rtol=%.3g, reached: %s)'
              % (len(full), nonzero, worst, args.rtol, detail))

    # 2. every path of the full run carries its transmitter, and each
    # transmitter has as many as its single run, at least one
    lines, tags = read_path_tags(os.path.join(args.full_dir, 'RayPath.vtk'))
    if len(tags) != lines:
        print('RESULT path_tags: FAIL (full run has %d paths but %d TransmitterID values)' % (lines, len(tags)))
        failures += 1
    else:
        counts = [tags.count(i) for i in range(len(singles))]
        expected = [paths for _, _, paths in singles]
        detail = ', '.join('%s=%d/%d' % (singles[i][0], counts[i], expected[i]) for i in range(len(singles)))
        if counts != expected or min(expected) == 0:
            print('RESULT path_tags: FAIL (full/single paths per transmitter: %s)' % detail)
            failures += 1
        else:
            print('RESULT path_tags: PASS (full/single paths per transmitter: %s)' % detail)

    if failures:
        print('%d check(s) failed. The generated configurations and results are in %s' % (failures, args.work_dir))
        sys.exit(1)


if __name__ == '__main__':
    main()
