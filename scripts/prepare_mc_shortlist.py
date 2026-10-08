"""Turn the shortlisted stars' raw Gaia DR3 rows into tools/mc_shortlist.c's input format, plus a
data-quality flag table.

The rows themselves come from ONE streaming pass over the 43 GB catalog CSV -- done with a fixed-
string grep for ",<source_id>," at idle CPU priority (see NOTES_method_selection.md, "Stage-2
Monte Carlo of the 99-star shortlist"), which can in principle also hit a different numeric
column. This script therefore re-checks the actual source_id column of every matched row against
the header, and collapses exact duplicate rows (the known ~0.4% duplicate-row issue in the
downloaded CSV) -- but refuses to silently pick between NON-identical rows for the same source.

Usage:
  python scripts/prepare_mc_shortlist.py <gaia_csv_for_header> <matched_rows.csv> <ids.txt>
                                         <stars_out.txt> <quality_out.tsv>
"""
import csv
import sys

ASTRO = ["ra", "dec", "parallax", "pmra", "pmdec", "radial_velocity"]
ERRS = ["ra_error", "dec_error", "parallax_error", "pmra_error", "pmdec_error", "radial_velocity_error"]
CORRS = ["ra_dec_corr", "ra_parallax_corr", "ra_pmra_corr", "ra_pmdec_corr", "dec_parallax_corr",
         "dec_pmra_corr", "dec_pmdec_corr", "parallax_pmra_corr", "parallax_pmdec_corr", "pmra_pmdec_corr"]
# Gaia DR3 parallax zero-point inputs (Lindegren et al. 2021), appended as five extra columns since
# 2026-10-08; tools/mc_shortlist.c applies the correction itself (src/parallax_zeropoint.c). A null
# is written as "nan" (e.g. pseudocolour for every 5-parameter solution).
ZPT_INPUTS = ["phot_g_mean_mag", "nu_eff_used_in_astrometry", "pseudocolour", "ecl_lat",
              "astrometric_params_solved"]
QUALITY = ["designation", "phot_g_mean_mag", "bp_rp", "parallax_over_error", "ruwe",
           "astrometric_params_solved", "astrometric_excess_noise_sig", "ipd_frac_multi_peak",
           "ipd_gof_harmonic_amplitude", "duplicated_source", "non_single_star",
           "radial_velocity_error", "rv_method_used", "rv_nb_transits", "rv_chisq_pvalue",
           "rv_renormalised_gof", "rv_amplitude_robust", "rv_template_teff", "grvs_mag",
           "phot_variable_flag", "teff_gspphot", "logg_gspphot", "distance_gspphot"]


def main():
    header_csv, matched, ids_path, stars_out, quality_out = sys.argv[1:6]
    with open(header_csv, newline="") as f:
        header = next(csv.reader(f))
    col = {name: i for i, name in enumerate(header)}
    wanted = [line.strip() for line in open(ids_path) if line.strip()]
    wanted_set = set(wanted)

    rows = {}
    n_lines = n_false_hits = n_exact_dups = 0
    with open(matched, newline="") as f:
        for r in csv.reader(f):
            n_lines += 1
            if len(r) != len(header):
                raise SystemExit(f"row with {len(r)} fields, header has {len(header)}")
            sid = r[col["source_id"]]
            if sid not in wanted_set:
                n_false_hits += 1
                continue
            if sid in rows:
                if rows[sid] == r:
                    n_exact_dups += 1
                    continue
                raise SystemExit(f"non-identical duplicate rows for {sid} -- inspect by hand")
            rows[sid] = r

    missing = [s for s in wanted if s not in rows]
    print(f"{n_lines} matched lines, {n_false_hits} false hits (id in another column), "
          f"{n_exact_dups} exact duplicate rows collapsed, {len(rows)} unique stars, "
          f"{len(missing)} missing: {missing}")

    def val(r, name):
        x = r[col[name]]
        return None if x in ("", "null", "NULL", "nan") else x

    n_null_corr = 0
    with open(stars_out, "w") as out, open(quality_out, "w") as q:
        out.write("# source_id " + " ".join(ASTRO + ERRS + CORRS + ZPT_INPUTS) + "\n")
        q.write("source_id\t" + "\t".join(QUALITY) + "\tnull_fields\n")
        for sid in wanted:
            if sid not in rows:
                continue
            r = rows[sid]
            nulls = []
            fields = [sid]
            for name in ASTRO + ERRS:
                x = val(r, name)
                if x is None:
                    raise SystemExit(f"{sid}: required field {name} is null")
                fields.append(x)
            for name in CORRS:
                x = val(r, name)
                if x is None:
                    nulls.append(name)
                    n_null_corr += 1
                    x = "0"  # documented assumption: missing correlation treated as uncorrelated
                fields.append(x)
            for name in ZPT_INPUTS:
                fields.append(val(r, name) or "nan")
            out.write(" ".join(fields) + "\n")
            q.write(sid + "\t" + "\t".join(val(r, n) or "null" for n in QUALITY) + "\t"
                    + (",".join(nulls) or "-") + "\n")
    print(f"{n_null_corr} null correlation values set to 0")


if __name__ == "__main__":
    main()
