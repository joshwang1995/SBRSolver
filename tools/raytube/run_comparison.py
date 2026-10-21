#!/usr/bin/env python3
"""Compare one sbrsolver run with the RayTubeTracing reference solver.

Steps (all paths are created under --work):
 1. convert the INI into a RayTubeTracing run directory (ini_to_raytube.py),
    checking the per-face materials against the Faces.vtk sbrsolver wrote;
 2. either run RayTubeTracing.exe there (--exe, optionally through
    --launcher such as wine64 on Linux) or take the committed reference
    output (--reference DIR holding ElectricField.csv and ReceivedPowers.csv);
 3. compare sbrsolver's ElectricField.csv and path counts with the
    reference (compare_with_raytube.py), over the receiver range the
    manifest (cases.json) names for the scene;
 4. with --exe and --reference, also check that the freshly generated
    reference still equals the committed one (so a stale committed file is
    noticed on the platform that can run the executable);
 5. with --save-reference DIR, copy the fresh reference output there to
    (re)generate the committed files.

Usage: run_comparison.py <scene.ini> --sbrsolver-output DIR --work DIR
           [--reference DIR] [--exe RayTubeTracing.exe] [--launcher CMD]
           [--manifest cases.json] [--save-reference DIR] [--rtol R]
Exit status 0 when everything agrees, 1 otherwise.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PYTHON = [sys.executable, "-I"]


def run(cmd, cwd=None, log=None):
    print("+ " + " ".join(cmd), flush=True)
    if log:
        with open(log, "w") as fh:
            result = subprocess.run(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT)
    else:
        result = subprocess.run(cmd, cwd=cwd)
    return result.returncode


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ini")
    parser.add_argument("--sbrsolver-output", required=True, help="directory sbrsolver wrote for this scene")
    parser.add_argument("--work", required=True, help="scratch directory for the converted scene and the reference run")
    parser.add_argument("--reference", help="committed reference directory (ElectricField.csv, ReceivedPowers.csv)")
    parser.add_argument("--exe", help="RayTubeTracing.exe to run; without it the committed reference is used")
    parser.add_argument("--launcher", help="program that runs the executable, e.g. wine64 on Linux")
    parser.add_argument("--manifest", default=os.path.join(HERE, "cases.json"))
    parser.add_argument("--save-reference", metavar="DIR", help="copy the fresh reference output to DIR")
    parser.add_argument("--rtol", type=float, default=1e-7)
    args = parser.parse_args()

    if not args.exe and not args.reference:
        print("run_comparison: give --exe to run the reference solver or --reference for committed output")
        return 2

    with open(args.manifest) as fh:
        manifest = json.load(fh)
    name = os.path.splitext(os.path.basename(args.ini))[0]
    case = manifest.get(name, {})
    tesselate = int(case.get("tesselate", 0))
    receivers = case.get("receivers")
    threads = int(case.get("threads", 0))

    work = os.path.abspath(args.work)
    ref_dir = os.path.join(work, "raytube")
    sbr_field = os.path.join(args.sbrsolver_output, "ElectricField.csv")
    sbr_paths = os.path.join(args.sbrsolver_output, "RayPath.vtk")
    sbr_faces = os.path.join(args.sbrsolver_output, "Faces.vtk")
    if not os.path.exists(sbr_field):
        print("run_comparison: %s not found (run sbrsolver first)" % sbr_field)
        return 1

    # 1. convert
    convert = PYTHON + [os.path.join(HERE, "ini_to_raytube.py"), args.ini, ref_dir, "--threads", str(threads)]
    if tesselate > 0:
        convert += ["--tesselate", str(tesselate)]
    if os.path.exists(sbr_faces):
        convert += ["--check-faces", sbr_faces]
    else:
        print("note: %s not found, per-face materials not cross-checked (set [scene] save_faces = true)" % sbr_faces)
    if run(convert) != 0:
        return 1

    # 2. reference output
    if args.exe:
        exe = os.path.abspath(args.exe)
        cmd = ([args.launcher] if args.launcher else []) + [exe]
        log = os.path.join(work, "raytube.log")
        status = run(cmd, cwd=ref_dir, log=log)
        fresh = os.path.join(ref_dir, "ElectricField.csv")
        if status != 0 or not os.path.exists(fresh):
            print("run_comparison: the reference solver failed (status %d), see %s" % (status, log))
            with open(log, errors="replace") as fh:
                print(fh.read()[-2000:])
            return 1
        reference = ref_dir
        if args.save_reference:
            os.makedirs(args.save_reference, exist_ok=True)
            for f in ("ElectricField.csv", "ReceivedPowers.csv"):
                shutil.copyfile(os.path.join(ref_dir, f), os.path.join(args.save_reference, f))
            print("saved the reference output to %s" % args.save_reference)
    else:
        reference = os.path.abspath(args.reference)

    # 3. compare sbrsolver with the reference
    compare = PYTHON + [os.path.join(HERE, "compare_with_raytube.py"), sbr_field, os.path.join(reference, "ElectricField.csv"),
                        "--rtol", str(args.rtol), "--paths", sbr_paths, os.path.join(reference, "ReceivedPowers.csv")]
    if receivers:
        compare += ["--receivers", receivers]
    ok = run(compare) == 0

    # 4. fresh reference against the committed one
    if args.exe and args.reference:
        print("checking the fresh reference output against the committed %s" % args.reference)
        drift = PYTHON + [os.path.join(HERE, "compare_with_raytube.py"), os.path.join(ref_dir, "ElectricField.csv"),
                          os.path.join(os.path.abspath(args.reference), "ElectricField.csv"), "--rtol", str(args.rtol), "--worst", "3"]
        if run(drift) != 0:
            print("run_comparison: the committed reference differs from what RayTubeTracing produces now; "
                  "regenerate it with --save-reference")
            ok = False

    print("run_comparison: %s %s" % (name, "PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
