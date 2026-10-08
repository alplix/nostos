"""Checks scripts/make_work_units.py (the data-preparation half of a future BOINC work generator)
against the rest of the pipeline:
  1. its 'Oumuamua element set is exactly src/iso_solutions.h's OUMUAMUA_ELEMENTS_BARY_1600 -- the
     barycentric 1600 elements whose exact inbound asymptote the BOINC app starts from (2026-10-08
     fix); a stale heliocentric set here would silently reintroduce the fixed bug in every WU;
  2. its velocity covariance is positive definite with the 7c per-axis sigmas (~0.020/0.020/0.022
     km/s) that tools/oumuamua_velocity_covariance.c prints;
  3. the parallax zero-point: --zpt none/global/lindegren write raw, raw + 0.017 and raw - Z mas,
     with Z from the official gaiadr3_zeropoint package (HIP 3757: Z = -0.024473 mas, the same value
     tests/test_parallax_zeropoint.c pins for the C port). Skipped, not failed, if the package is
     not installed (pip install gaiadr3-zeropoint).
Run: python tests/test_make_work_units.py   (part of `make test`)
"""
import csv
import math
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import make_work_units as mwu  # noqa: E402

failures = 0


def check(ok, what):
    global failures
    print(("  ok:   " if ok else "  FAIL: ") + what)
    if not ok:
        failures += 1


print("test_make_work_units:")

hdr = open(os.path.join(ROOT, "src", "iso_solutions.h")).read()
m = re.search(r"OUMUAMUA_ELEMENTS_BARY_1600 = \{(.*?)\};", hdr, re.S)
c_el = {k: float(v) for k, v in re.findall(r"\.(\w+) = ([-0-9.eE+]+)", m.group(1))}
check(c_el == mwu.OUMUAMUA_ELEMENTS,
      "make_work_units.py elements == iso_solutions.h OUMUAMUA_ELEMENTS_BARY_1600 (exact)")

c = mwu.OUMUAMUA_VEL_COV_KPCMYR2
cov = [[c[0], c[1], c[2]], [c[1], c[3], c[4]], [c[2], c[4], c[5]]]
kms = 977.7922216807892  # NOSTOS_KMS_PER_KPCMYR
sig = [math.sqrt(cov[i][i]) * kms for i in range(3)]
check(all(abs(s - r) < 2e-4 for s, r in zip(sig, (0.020198, 0.019906, 0.021668))),
      f"covariance per-axis sigma {sig[0]:.6f} {sig[1]:.6f} {sig[2]:.6f} km/s == 7c Table-1 values")
d1 = cov[0][0]
d2 = cov[0][0] * cov[1][1] - cov[0][1] ** 2
d3 = (cov[0][0] * (cov[1][1] * cov[2][2] - cov[1][2] ** 2) - cov[0][1] * (cov[0][1] * cov[2][2] - cov[1][2] * cov[0][2])
      + cov[0][2] * (cov[0][1] * cov[1][2] - cov[1][1] * cov[0][2]))
check(d1 > 0 and d2 > 0 and d3 > 0, "covariance positive definite (all leading minors > 0)")

# A two-star CSV with the zero-point inputs: HIP 3757 (real DR3 values) and a 6p star whose
# recipe input is null (pseudocolour missing) -> global fallback.
cols = ["source_id", "ra", "dec", "parallax", "pmra", "pmdec", "radial_velocity", "phot_g_mean_mag",
        "nu_eff_used_in_astrometry", "pseudocolour", "ecl_lat", "astrometric_params_solved"]
rows = [["2525688198820543360", "12.056371339593694", "-5.136005920741353", "40.22675666148508",
         "181.89740146004812", "-120.13173095843229", "26.607231", "11.169264", "1.3237802", "null",
         "-9.490552045804828", "31"],
        ["1", "10.0", "10.0", "5.0", "1.0", "1.0", "10.0", "15.0", "null", "null", "20.0", "95"]]
try:
    import zero_point  # noqa: F401
    have_pkg = True
except ImportError:
    have_pkg = False

with tempfile.TemporaryDirectory() as td:
    path = os.path.join(td, "stars.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        w.writerows(rows)
    modes = ["none", "global"] + (["lindegren"] if have_pkg else [])
    plx = {}
    for mode in modes:
        out = os.path.join(td, mode)
        r = subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "make_work_units.py"), path, out,
                            "--zpt", mode], capture_output=True, text=True)
        check(r.returncode == 0, f"--zpt {mode}: script runs")
        lines = open(os.path.join(out, "wu_000000_in")).read().splitlines()
        check(lines[0].split() == [repr(mwu.OUMUAMUA_ELEMENTS[k]) for k in ("e", "q_au", "i_deg", "node_deg", "peri_deg")],
              f"--zpt {mode}: WU line 1 carries the barycentric 1600 elements")
        plx[mode] = {ln.split()[0]: float(ln.split()[3]) for ln in lines[4:]}
    check(plx["none"]["2525688198820543360"] == 40.22675666148508, "--zpt none writes the raw parallax")
    check(abs(plx["global"]["2525688198820543360"] - (40.22675666148508 + 0.017)) < 1e-12,
          "--zpt global writes parallax + 0.017 mas")
    if have_pkg:
        check(abs(plx["lindegren"]["2525688198820543360"] - (40.22675666148508 + 0.024473)) < 1e-6,
              "--zpt lindegren: HIP 3757 parallax - Z with Z = -0.024473 mas (official package)")
        check(abs(plx["lindegren"]["1"] - 5.017) < 1e-12, "--zpt lindegren: null recipe input -> global fallback")
    else:
        print("  SKIP: gaiadr3_zeropoint not installed -- lindegren mode not tested")

if failures:
    print(f"test_make_work_units: {failures} FAILURE(S)")
    sys.exit(1)
print("test_make_work_units: all passed")
