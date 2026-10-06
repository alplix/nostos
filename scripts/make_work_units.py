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
"""
import argparse
import csv
import os
import sys

# 1I/'Oumuamua, same real JPL Horizons elements used throughout this project's tests (see
# tests/test_horizons_crosscheck.c, tests/test_pipeline_oumuamua_hip3757.c).
OUMUAMUA_ELEMENTS = {
    "e": 1.200954455284624,
    "q_au": 0.2559200355462156,
    "i_deg": 122.7402198359025,
    "node_deg": 24.59683473134646,
    "peri_deg": 241.8043518577817,
}

# Placeholder velocity-uncertainty covariance (kpc/Myr)^2 -- NOT yet derived from 'Oumuamua's real
# published orbital-solution covariance matrix (that conversion -- from whatever units/parametrization
# the JPL/MPC solution's covariance is quoted in, into this project's internal kpc/Myr^2 velocity
# covariance -- is itself unfinished work, a real TODO, not a decision made here). Small enough that
# mc_search_stage1()'s 100 clones will cluster tightly around the best-fit trajectory until a real
# covariance is used.
PLACEHOLDER_VEL_COV_KPCMYR2 = (1e-6, 0.0, 0.0, 1e-6, 0.0, 1e-6)  # cov00 cov01 cov02 cov11 cov12 cov22

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
        f.write(" ".join(str(x) for x in PLACEHOLDER_VEL_COV_KPCMYR2) + "\n")
        f.write(f"{DEFAULT_SEED_BASE}\n")
        f.write(f"{DEFAULT_TOTAL_TIME_MYR} {DEFAULT_N_STEPS} {DEFAULT_DISTANCE_THRESHOLD_PC}\n")
        for star in stars_chunk:
            f.write(
                f"{star['source_id']} {star['ra']} {star['dec']} {star['parallax']} "
                f"{star['pmra']} {star['pmdec']} {star['rv']}\n"
            )


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
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    chunk = []
    n_rows = 0
    n_written = 0
    n_chunks = 0

    with open(args.gaia_csv, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            n_rows += 1
            if args.max_rows and n_rows > args.max_rows:
                break
            try:
                parallax = float(row["parallax"])
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
            except (ValueError, KeyError):
                continue
            chunk.append(star)
            n_written += 1
            if len(chunk) >= args.chunk_size:
                write_work_unit(os.path.join(args.out_dir, f"wu_{n_chunks:06d}_in"), chunk)
                n_chunks += 1
                chunk = []

    if chunk:
        write_work_unit(os.path.join(args.out_dir, f"wu_{n_chunks:06d}_in"), chunk)
        n_chunks += 1

    print(f"{n_rows} rows read, {n_written} stars written across {n_chunks} work units "
          f"(chunk size {args.chunk_size}) in {args.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
