"""Reference values for tests/test_parallax_zeropoint.c from the Gaia collaboration's OFFICIAL
parallax zero-point package (`pip install gaiadr3-zeropoint`, version 0.1.0; Lindegren et al. 2021).

Prints one C initializer line per case: {G, nu_eff, pseudocolour, ecl_lat, params_solved, Z_mas}.
The cases cover 5p and 6p solutions, every G interval of the coefficient tables (including the
6/21 clamps), colours inside and outside the calibrated range, both ecliptic hemispheres, and
three real Gaia DR3 stars from the 'Oumuamua shortlist (HIP 3757 among them), if the shortlist
CSV is available.

Usage: python scripts/zeropoint_crosscheck.py [shortlist_rows.csv]
"""
import csv
import math
import sys
import warnings

import numpy as np
from zero_point import zpt

zpt.load_tables()
warnings.simplefilter("ignore")  # out-of-range cases are deliberate here


def get_zpt(g, nu, pc, b, ps):
    # Array call: the package's scalar path uses np.can_cast on Python floats, which NumPy >= 2
    # rejects (NEP 50); the array path is the same arithmetic. NaN = column null in the catalogue.
    return float(zpt.get_zpt(np.array([g]), np.array([nu]), np.array([pc]), np.array([b]),
                             np.array([ps]))[0])

cases = [
    # G, nu_eff, pseudocolour, ecl_lat, params_solved
    (5.0, 1.50, math.nan, 10.0, 31),
    (8.3, 1.55, math.nan, -35.0, 31),
    (10.95, 1.40, math.nan, 60.0, 31),
    (11.5, 1.62, math.nan, -80.0, 31),
    (12.5, 1.30, math.nan, 5.0, 31),
    (13.0, 1.75, math.nan, 45.0, 31),
    (14.2, 1.48, math.nan, -20.0, 31),
    (16.0, 1.20, math.nan, 70.0, 31),
    (17.0, 1.85, math.nan, -50.0, 31),
    (18.3, 1.58, math.nan, 0.0, 31),
    (19.5, 1.45, math.nan, 30.0, 31),
    (20.5, 1.35, math.nan, -65.0, 31),
    (21.5, 1.50, math.nan, 15.0, 31),
    (15.5, 2.00, math.nan, 25.0, 31),   # nu_eff out of range
    (9.0, math.nan, 1.50, 40.0, 95),
    (12.0, math.nan, 1.30, -10.0, 95),
    (14.0, math.nan, 1.65, 75.0, 95),
    (16.8, math.nan, 1.20, -40.0, 95),  # pseudocolour below range
    (19.2, math.nan, 1.80, 55.0, 95),   # pseudocolour above range
    (20.7, math.nan, 1.45, -5.0, 95),
]

real = []
if len(sys.argv) > 1:
    wanted = {"2525688198820543360", "2618904962125653888", "4399766976222478336"}
    with open(sys.argv[1], newline="") as f:
        rd = csv.reader(f)
        hdr = next(rd)
        col = {n: i for i, n in enumerate(hdr)}
        for r in rd:
            if r[col["source_id"]] in wanted:
                def fv(n):
                    x = r[col[n]]
                    return math.nan if x in ("", "null") else float(x)
                real.append((r[col["source_id"]], fv("phot_g_mean_mag"), fv("nu_eff_used_in_astrometry"),
                             fv("pseudocolour"), fv("ecl_lat"), int(r[col["astrometric_params_solved"]]),
                             fv("parallax")))

for g, nu, pc, b, ps in cases:
    z = get_zpt(g, nu, pc, b, ps)
    print(f"    {{ {g!r}, {nu!r}, {pc!r}, {b!r}, {ps}, {z!r} }},")
print("    // real Gaia DR3 stars")
for sid, g, nu, pc, b, ps, plx in real:
    z = get_zpt(g, nu, pc, b, ps)
    print(f"    {{ {g!r}, {nu!r}, {pc!r}, {b!r}, {ps}, {z!r} }}, // {sid}, parallax {plx}")
