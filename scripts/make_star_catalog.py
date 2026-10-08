"""Writes the compact star catalogue the BOINC work generator cuts work units from.

One line per usable Gaia DR3 6D star, in catalogue order, gzip-compressed:
    source_id ra_deg dec_deg parallax_mas pmra_masyr pmdec_masyr rv_kms
with the parallax already zero-point corrected (Lindegren et al. 2021 by default) -- exactly the
star lines of a work-unit input file (src/main_nostos_boinc.cpp, src/nostos_io.h), so the server
only has to prepend an ISO header to a slice of this file. Same row filtering and correction as
scripts/make_work_units.py (it reuses that script's code); the precision kept (1e-9 deg, 1e-6 mas,
1e-5 mas/yr, 1e-4 km/s) is far below every Gaia uncertainty.

Usage:
    python scripts/make_star_catalog.py <gaia_csv> <out.txt.gz> [--zpt lindegren|global|none] [--max-rows N]
"""
import argparse
import csv
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_work_units import zero_points  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gaia_csv")
    ap.add_argument("out_gz")
    ap.add_argument("--zpt", choices=["lindegren", "global", "none"], default="lindegren")
    ap.add_argument("--max-rows", type=int, default=None)
    args = ap.parse_args()

    n_rows = n_written = n_fallback = 0
    pending = []
    tmp = args.out_gz + ".part"
    with gzip.open(tmp, "wt", compresslevel=6, newline="\n") as out:
        def flush():
            nonlocal pending, n_written, n_fallback
            zs, nfb = zero_points(pending, args.zpt)
            n_fallback += nfb
            lines = []
            for row, z in zip(pending, zs):
                try:
                    parallax = float(row["parallax"]) - z
                    if parallax <= 0:
                        continue
                    lines.append("%d %.9f %.9f %.6f %.5f %.5f %.4f\n" % (
                        int(row["source_id"]), float(row["ra"]), float(row["dec"]), parallax,
                        float(row["pmra"]), float(row["pmdec"]), float(row["radial_velocity"])))
                except (ValueError, KeyError, TypeError):
                    continue
            out.write("".join(lines))
            n_written += len(lines)
            pending = []

        with open(args.gaia_csv, newline="") as f:
            for row in csv.DictReader(f):
                n_rows += 1
                if args.max_rows and n_rows > args.max_rows:
                    break
                pending.append(row)
                if len(pending) >= 50000:
                    flush()
                    if n_rows % 1000000 == 0:
                        print(f"{n_rows} rows, {n_written} stars", file=sys.stderr, flush=True)
        flush()
    os.replace(tmp, args.out_gz)
    print(f"done: {n_rows} rows read, {n_written} stars written to {args.out_gz} "
          f"(zero-point {args.zpt}, {n_fallback} on the global fallback)", file=sys.stderr)


if __name__ == "__main__":
    main()
