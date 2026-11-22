"""Cross-language validation: Python vs C++ (and MATLAB when its CSVs exist).
Usage: python tools/cross_validate.py [--tol 1e-10]
Expects the C++ binary at build/nfs_run. MATLAB outputs (from MATLAB/tests/run_validation_cases.m)
are read from results/validation/matlab/<case>.csv if present; otherwise MATLAB is reported NOT RUN."""
import argparse, glob, os, subprocess, sys, datetime, platform
import numpy as np
ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(ROOT, "Python"))
from nanosim.core.case import read_case, run_case, write_csv

ap = argparse.ArgumentParser(); ap.add_argument("--tol", type=float, default=1e-10)
tol = ap.parse_args().tol
cases = sorted(glob.glob(os.path.join(ROOT, "simulations", "validation_cases", "*.txt")))
outdir = os.path.join(ROOT, "results", "validation"); os.makedirs(outdir, exist_ok=True)
for sub in ("python", "cpp", "matlab"): os.makedirs(os.path.join(outdir, sub), exist_ok=True)
exe = os.path.join(ROOT, "build", "nfs_run")

def rows_minA_ok(P): return P[:, 3].min() >= -1e-12

rows, ok_all = [], True
for cf in cases:
    name = os.path.splitext(os.path.basename(cf))[0]
    pyp, cpp, mat = (os.path.join(outdir, s, name + ".csv") for s in ("python", "cpp", "matlab"))
    write_csv(pyp, run_case(read_case(cf)))
    subprocess.run([exe, cf, cpp], check=True)
    P = np.loadtxt(pyp, delimiter=",", skiprows=1); C = np.loadtxt(cpp, delimiter=",", skiprows=1)
    dpc = np.max(np.abs(P[:, 1:] - C[:, 1:])); dwl = np.max(np.abs(P[:, 0] - C[:, 0]))
    dpm = None
    if os.path.exists(mat):
        Mx = np.loadtxt(mat, delimiter=",", skiprows=1); dpm = np.max(np.abs(P[:, 1:] - Mx[:, 1:]))
    ok = rows_minA_ok(P) and dpc <= tol and dwl <= 1e-9 and (dpm is None or dpm <= tol); ok_all &= ok
    rows.append((name, len(P), dwl, dpc, dpm, P[:, 3].min(), ok))

with open(os.path.join(outdir, "cross_language_report.md"), "w") as f:
    f.write(f"# Cross-language validation report\n\nGenerated {datetime.datetime.now():%Y-%m-%d %H:%M} on {platform.platform()}, "
            f"Python {platform.python_version()}, NumPy {np.__version__}.\n\nTolerance (max |dR|,|dT|,|dA|): {tol:g}.\n"
            "MATLAB column 'NOT RUN' means MATLAB was unavailable and no MATLAB CSVs were supplied.\n\n"
            "| case | points | max d(wl) nm | max abs diff Py-C++ | max abs diff Py-MATLAB | min A (passivity: >= -1e-12) | status |\n|---|---|---|---|---|---|---|\n")
    for n, N, dw, d1, d2, e, ok in rows:
        f.write(f"| {n} | {N} | {dw:.1e} | {d1:.2e} | {'NOT RUN' if d2 is None else f'{d2:.2e}'} | {e:.2e} | {'PASS' if ok else 'FAIL'} |\n")
    f.write(f"\n**Overall: {'PASS' if ok_all else 'FAIL'}** ({len(rows)} cases)\n")
print(open(os.path.join(outdir, "cross_language_report.md")).read())
sys.exit(0 if ok_all else 1)
