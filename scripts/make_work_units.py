"""Chunks a downloaded Gaia bulk CSV (scripts/download_filter_6d.py's output) into BOINC work-unit
input files matching src/main_nostos_boinc.cpp / src/nostos_io.h's exact format -- the data-
preparation half of a future work generator.

Deliberately does NOT call BOINC's create_work or touch any database -- there's no BOINC project
instance set up for Nostos on any server yet (unlike bitboinc-app, a mature production project this
repo shares infrastructure with but is NOT the same project), and that's a real deployment decision
for the project owner to make explicitly, not something to build blind ahead of time (see
README.md's "BOINC app" roadmap item). This script only produces the plain input files themselves,
which can be fed to the compiled nostos_cpu_linux_x86_64 binary directly (as already verified
end-to-end on the bitboinc-app server, see NOTES_method_selection.md) to confirm they're read
correctly -- the actual registration/distribution step is left for later.

Usage:
    python scripts/make_work_units.py <gaia_csv> <out_dir> [--chunk-size N] [--max-rows N]
                                      [--zpt lindegren|global|none]

Gaia DR3 parallax zero-point (2026-10-08, default ON, same as src/gaia_scan.c and
tools/mc_shortlist.c): the parallax written to each work unit is parallax - Z, with Z from the Gaia
collaboration's official `gaiadr3_zeropoint` package (Lindegren et al. 2021; `pip install
gaiadr3-zeropoint`). The BOINC app itself never sees the recipe inputs, so the correction has to
happen here. Stars whose recipe inputs are null/invalid get the global -0.017 mas offset (counted);
--zpt global applies -0.017 mas to every star; --zpt none writes the raw catalogue parallax.
"""
import argparse
import csv
import math
import os
import sys
import warnings

ZPT_GLOBAL_MAS = -0.017  # src/parallax_zeropoint.h NOSTOS_ZPT_GLOBAL_MAS

# 1I/'Oumuamua: JPL Horizons BARYCENTRIC osculating elements at 1600-01-01 TDB (solution JPL#16 =
# 7c), the same set as src/iso_solutions.h's OUMUAMUA_ELEMENTS_BARY_1600. src/main_nostos_boinc.cpp
# starts the ISO on the exact inbound asymptote of these elements with the total solar-system GM
# (2026-10-08 fix, NOTES_method_selection.md "Residual distance gap"; the heliocentric 2017-epoch
# set used before pointed 0.077 deg away and, evaluated 9.08 yr before perihelion, was 2.3% too fast).
OUMUAMUA_ELEMENTS = {
    "e": 1.198093111333325,
    "q_au": 0.2526436572799559,
    "i_deg": 122.7927648087493,
    "node_deg": 24.24788585536907,
    "peri_deg": 241.8944051367842,
}

# 'Oumuamua's asymptotic-velocity covariance (kpc/Myr)^2 in this project's Galactic axes: Bailer-Jones
# et al. 2018 (AJ 156, 205) Table 1, solution 7c, rotated from barycentric ICRF by
# tools/oumuamua_velocity_covariance.c (2026-10-08). Per-axis 1-sigma ~0.020, 0.020, 0.022 km/s. This
# replaces the JPL SBDB osculating-element covariance propagated to the velocity (~0.001 km/s), which
# left out the non-gravitational parameters' uncertainty -- a real ~20-40x underestimate. The matrix
# is nearly singular (published correlations ~0.99) but positive definite (smallest eigenvalue
# 1.6e-15), so mc_search_stage1()'s Cholesky accepts it.
OUMUAMUA_VEL_COV_KPCMYR2 = (
    4.2671045967e-10, 4.1958328879e-10, 4.5772944273e-10,
    4.1446922995e-10, 4.4977060095e-10,
    4.9105925054e-10,
)  # cov00 cov01 cov02 cov11 cov12 cov22

DEFAULT_SEED_BASE = 1  # a fixed, reproducible default -- a real generator would pick this per
                        # work unit (e.g. derived from the ISO + chunk index) for genuinely
                        # independent Monte Carlo draws across work units; left as a single fixed
                        # value here since this script doesn't produce more than one real ISO's
                        # worth of work units yet.
DEFAULT_TOTAL_TIME_MYR = -100.0  # matches NOTES_method_selection.md's settled 100 Myr backward window
DEFAULT_N_STEPS = 100000  # dt=1e-3 Myr=1000yr over the full 100 Myr window
DEFAULT_DISTANCE_THRESHOLD_PC = 50.0  # matches gaia_scan.c's own widened threshold -- see its
                                       # comment for why 5pc (calibrated on a different model's
                                       # scale) was too tight for this project's own numbers


def write_work_unit(path, stars_chunk):
    with open(path, "w") as f:
        el = OUMUAMUA_ELEMENTS
        f.write(f"{el['e']} {el['q_au']} {el['i_deg']} {el['node_deg']} {el['peri_deg']}\n")
        f.write(" ".join(str(x) for x in OUMUAMUA_VEL_COV_KPCMYR2) + "\n")
        f.write(f"{DEFAULT_SEED_BASE}\n")
        f.write(f"{DEFAULT_TOTAL_TIME_MYR} {DEFAULT_N_STEPS} {DEFAULT_DISTANCE_THRESHOLD_PC}\n")
        for star in stars_chunk:
            f.write(
                f"{star['source_id']} {star['ra']} {star['dec']} {star['parallax']} "
                f"{star['pmra']} {star['pmdec']} {star['rv']}\n"
            )


def _float_or_nan(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def zero_points(rows, mode):
    """Z (mas) for each raw CSV row dict; returns (list_of_Z, n_global_fallback)."""
    if mode == "none":
        return [0.0] * len(rows), 0
    if mode == "global":
        return [ZPT_GLOBAL_MAS] * len(rows), 0
    import numpy as np
    from zero_point import zpt as _zpt
    if not getattr(zero_points, "_loaded", False):
        _zpt.load_tables()
        zero_points._loaded = True
    g = np.array([_float_or_nan(r.get("phot_g_mean_mag")) for r in rows])
    nu = np.array([_float_or_nan(r.get("nu_eff_used_in_astrometry")) for r in rows])
    pc = np.array([_float_or_nan(r.get("pseudocolour")) for r in rows])
    b = np.array([_float_or_nan(r.get("ecl_lat")) for r in rows])
    ps = np.array([_float_or_nan(r.get("astrometric_params_solved")) for r in rows])
    ok = (ps == 31) | (ps == 95)
    z = np.full(len(rows), np.nan)
    if ok.any():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # out-of-range G/colour: evaluated anyway, like the C port
            # Array call: the package's scalar path breaks under NumPy >= 2 (NEP 50).
            z[ok] = _zpt.get_zpt(g[ok], nu[ok], pc[ok], b[ok], ps[ok])
    bad = ~np.isfinite(z)
    z[bad] = ZPT_GLOBAL_MAS
    return [float(x) for x in z], int(bad.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gaia_csv")
    ap.add_argument("out_dir")
    # Default picked from a real timing measurement (2026-10-06, on the bitboinc-app production
    # server -- a capable machine, likely faster than a typical volunteer host): 10 stars x
    # MC_N_ISO_CLONES=100 clones x DEFAULT_N_STEPS=100000 steps took 23.05s wall-clock, i.e.
    # ~2.3s/star. 150 stars at that rate is ~5.75 minutes on THIS hardware -- a reasonable BOINC
    # task duration target even accounting for slower volunteer machines, and far short of the
    # ~38 minutes 1000 stars/WU (an earlier, un-benchmarked default) was measured to actually take.
    # Chunk size is the right lever to tune here, not MC_N_ISO_CLONES (a scientifically-justified
    # constant, see NOTES_method_selection.md) or DEFAULT_N_STEPS (a resolution/precision choice,
    # not a free performance knob) -- it doesn't change what's being computed per star at all, only
    # how many stars one work unit covers.
    ap.add_argument("--chunk-size", type=int, default=150)
    ap.add_argument("--max-rows", type=int, default=None)
    ap.add_argument("--zpt", choices=["lindegren", "global", "none"], default="lindegren")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    chunk = []
    n_rows = 0
    n_written = 0
    n_chunks = 0
    n_zpt_fallback = 0
    pending = []  # raw rows, zero-pointed in batches (one vectorised call per batch)

    def flush_pending():
        nonlocal chunk, n_written, n_chunks, n_zpt_fallback, pending
        zs, nfb = zero_points(pending, args.zpt)
        n_zpt_fallback += nfb
        for row, z in zip(pending, zs):
            try:
                # Same order as src/gaia_scan.c: correct first, then reject non-positive parallaxes.
                parallax = float(row["parallax"]) - z
                if parallax <= 0:
                    continue
                star = {
                    "source_id": int(row["source_id"]),
                    "ra": float(row["ra"]),
                    "dec": float(row["dec"]),
                    "parallax": parallax,
                    "pmra": float(row["pmra"]),
                    "pmdec": float(row["pmdec"]),
                    "rv": float(row["radial_velocity"]),
                }
            except (ValueError, KeyError, TypeError):
                continue
            chunk.append(star)
            n_written += 1
            if len(chunk) >= args.chunk_size:
                write_work_unit(os.path.join(args.out_dir, f"wu_{n_chunks:06d}_in"), chunk)
                n_chunks += 1
                chunk = []
        pending = []

    with open(args.gaia_csv, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            n_rows += 1
            if args.max_rows and n_rows > args.max_rows:
                break
            pending.append(row)
            if len(pending) >= 20000:
                flush_pending()
    flush_pending()

    if chunk:
        write_work_unit(os.path.join(args.out_dir, f"wu_{n_chunks:06d}_in"), chunk)
        n_chunks += 1

    print(f"{n_rows} rows read, {n_written} stars written across {n_chunks} work units "
          f"(chunk size {args.chunk_size}) in {args.out_dir}; parallax zero-point: {args.zpt}"
          f"{f' ({n_zpt_fallback} stars on the global -0.017 mas fallback)' if args.zpt == 'lindegren' else ''}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
