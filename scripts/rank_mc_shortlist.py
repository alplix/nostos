"""Merge tools/mc_shortlist.c's per-star summary (and per-sample dump) with the full-catalog
scan's nominal numbers and the Gaia quality flags from scripts/prepare_mc_shortlist.py, and rank
the stars by Monte Carlo robustness.

Robustness ranking (a triage order, not a probability of origin): primarily the fraction of MC
samples passing within 1 pc, then within 2 pc, then the median distance. Also reports the joint
"close AND slow" fraction (d < 1 pc and v < 20 km/s, computed per sample).

Quality flags (any one is a reason for extra skepticism, not automatic rejection):
  RUWE>1.4 (poor single-star astrometric fit -- possible unresolved binary), DUP
  (duplicated_source), NSS (Gaia DR3 non_single_star != 0), RVERR>5 (radial_velocity_error
  > 5 km/s), RVVAR (rv_chisq_pvalue < 0.01: RV epochs inconsistent with constant velocity --
  possible spectroscopic binary), IPD (ipd_frac_multi_peak > 10: multiple-peak images, likely a
  close double), PLX<10 (parallax_over_error < 10).

Usage: python scripts/rank_mc_shortlist.py <summary.txt> <samples.txt> <quality.tsv>
                                           <scan_full.txt> <out_table.txt>
"""
import csv
import sys
from collections import defaultdict


def flags_for(q):
    def f(x):
        try:
            return float(x)
        except ValueError:
            return None
    out = []
    if (f(q["ruwe"]) or 0) > 1.4:
        out.append("RUWE>1.4")
    if q["duplicated_source"] in ("True", "true", "1"):
        out.append("DUP")
    if q["non_single_star"] not in ("0", "null"):
        out.append("NSS")
    if (f(q["radial_velocity_error"]) or 0) > 5:
        out.append("RVERR>5")
    p = f(q["rv_chisq_pvalue"])
    if p is not None and p < 0.01:
        out.append("RVVAR")
    if (f(q["ipd_frac_multi_peak"]) or 0) > 10:
        out.append("IPD")
    if (f(q["parallax_over_error"]) or 1e9) < 10:
        out.append("PLX<10")
    return out


def main():
    summary_path, samples_path, quality_path, scan_path, out_path = sys.argv[1:6]
    cols = None
    stars = []
    for line in open(summary_path):
        if line.startswith("# source_id"):
            cols = line[2:].split()
            continue
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        stars.append(dict(zip(cols, parts)))

    joint = defaultdict(lambda: [0, 0])
    for line in open(samples_path):
        if line.startswith("#"):
            continue
        sid, _, d, t, v, _edge = line.split()
        joint[sid][1] += 1
        if float(d) < 1.0 and float(v) < 20.0:
            joint[sid][0] += 1

    quality = {r["source_id"]: r for r in csv.DictReader(open(quality_path), delimiter="\t")}

    wanted = {s["source_id"] for s in stars}
    scan = {}
    for line in open(scan_path):
        sid = line.split(" ", 1)[0]
        if sid in wanted and sid not in scan:
            kv = dict(p.split("=") for p in line.split()[1:] if "=" in p)
            scan[sid] = (float(kv["gravity_d_pc"]), float(kv["gravity_t_myr"]), float(kv["rel_speed_kms"]))

    for s in stars:
        for k in list(s):
            if k != "source_id":
                s[k] = float(s[k])
        n = joint[s["source_id"]]
        s["f_close_slow"] = n[0] / n[1] if n[1] else 0.0
        s["flags"] = flags_for(quality[s["source_id"]])
        s["designation"] = quality[s["source_id"]]["designation"]
        s["G"] = quality[s["source_id"]]["phot_g_mean_mag"]
        s["ruwe"] = quality[s["source_id"]]["ruwe"]
        s["rv_err"] = quality[s["source_id"]]["radial_velocity_error"]
        s["scan"] = scan.get(s["source_id"])

    stars.sort(key=lambda s: (-s["f_d_lt_1"], -s["f_d_lt_2"], s["d_p50"]))

    with open(out_path, "w") as out:
        out.write(f"# Stage-2 Monte Carlo of a {len(stars)}-star 'Oumuamua shortlist (tools/mc_shortlist.c), ranked by\n"
                  "# robustness: fraction of samples with d<1pc, then d<2pc, then median d. d in pc, t in Myr,\n"
                  "# v in km/s. 'scan' = the full-catalog gaia_scan.c nominal; 'nom' = this run's own nominal\n"
                  "# (same code path, should match). f_cs = fraction with d<1pc AND v<20km/s.\n"
                  "# Measurement uncertainty only (Gaia star + 'Oumuamua asymptote); no Galactic-potential\n"
                  "# uncertainty and no Gaia parallax zero-point correction. Since the 2026-10-08 starting-state\n"
                  "# fix the pipeline reproduces published encounter distances to a few hundredths of a pc\n"
                  "# (NOTES_method_selection.md, \"Residual distance gap\").\n")
        out.write("# rank source_id n_valid scan_d nom_d d_p05 d_p16 d_p50 d_p84 d_p95 "
                  "t_p16 t_p50 t_p84 v_p16 v_p50 v_p84 f_d<0.5 f_d<1 f_d<2 f_cs G ruwe rv_err flags\n")
        for i, s in enumerate(stars, 1):
            sd = s["scan"][0] if s["scan"] else float("nan")
            out.write(f"{i:3d} {s['source_id']:>20s} {int(s['n_valid']):5d} {sd:7.4f} {s['nom_d_pc']:7.4f} "
                      f"{s['d_p05']:7.3f} {s['d_p16']:7.3f} {s['d_p50']:7.3f} {s['d_p84']:7.3f} {s['d_p95']:7.3f} "
                      f"{s['t_p16']:8.3f} {s['t_p50']:8.3f} {s['t_p84']:8.3f} "
                      f"{s['v_p16']:6.2f} {s['v_p50']:6.2f} {s['v_p84']:6.2f} "
                      f"{s['f_d_lt_0.5']:.3f} {s['f_d_lt_1']:.3f} {s['f_d_lt_2']:.3f} {s['f_close_slow']:.3f} "
                      f"{float(s['G']):5.2f} {float(s['ruwe']):6.2f} {float(s['rv_err']):5.2f} "
                      f"{','.join(s['flags']) or '-'}\n")

    # Consistency check: this run's nominal must reproduce the full-catalog scan.
    worst = max(abs(s["nom_d_pc"] - s["scan"][0]) / s["scan"][0] for s in stars if s["scan"])
    print(f"max relative |nominal - scan| distance difference: {worst:.2e}")
    print(f"window-limited samples (total): {int(sum(s['n_window_limited'] for s in stars))}")
    print(f"rejected samples (total): {int(sum(s['n_rejected'] for s in stars))}")


if __name__ == "__main__":
    main()
