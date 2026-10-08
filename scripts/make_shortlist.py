"""Build the stage-2 shortlist from one or more src/gaia_scan.c candidate files (e.g. the two
--byte-range halves of one scan): deduplicate by source_id (keep the first occurrence -- the known
~0.4% duplicate-row issue of the downloaded CSV), drop WINDOW_LIMITED rows, and keep
gravity_d_pc < --dmax AND rel_speed_kms < --vmax (default 2 pc / 50 km/s, the Bailer-Jones-style
plausibility cut used since 2026-10-07).

Writes the human-readable list (same format as scan_oumuamua_full_asym_most_promising.txt: a count
line, a blank line, then one star per line sorted by distance) and a plain ID list. With
--union-ids, the ID list is the union (new shortlist first, then any extra IDs not already in it),
so an older shortlist can be followed star by star; with --have-rows, also writes the IDs whose
DR3 rows are not yet in the given CSV (they need one more grep pass over the catalogue).

Replaces the scratch-only builder used for the 163-star run (NOTES_method_selection.md). Parsing
trap it avoids: the count line's first token is numeric, so only lines whose first token is an
integer AND that contain key=value fields count as candidates.

Usage:
  python scripts/make_shortlist.py <out_promising.txt> <out_ids.txt> <scan.txt> [<scan2.txt> ...]
         [--dmax 2.0] [--vmax 50.0] [--union-ids old_ids_or_promising.txt]
         [--have-rows rows.csv --missing-ids out_missing.txt]
"""
import argparse
import csv


def parse_scan(paths):
    rows, n_dup, n_wl = {}, 0, 0
    for path in paths:
        for line in open(path):
            parts = line.split()
            if not parts or not parts[0].isdigit():
                continue
            kv = dict(p.split("=", 1) for p in parts[1:] if "=" in p)
            if "gravity_d_pc" not in kv:
                continue
            if parts[0] in rows:
                n_dup += 1
                continue
            wl = "WINDOW_LIMITED" in parts
            n_wl += wl
            rows[parts[0]] = (float(kv["gravity_d_pc"]), float(kv["rel_speed_kms"]),
                              float(kv["gravity_t_myr"]), wl, kv.get("zpt_mas"))
    return rows, n_dup, n_wl


def read_ids(path):
    ids = []
    for line in open(path):
        p = line.split()
        if p and p[0].isdigit() and (len(p) == 1 or "d=" in line):
            ids.append(p[0])
    return ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out_promising")
    ap.add_argument("out_ids")
    ap.add_argument("scan", nargs="+")
    ap.add_argument("--dmax", type=float, default=2.0)
    ap.add_argument("--vmax", type=float, default=50.0)
    ap.add_argument("--union-ids")
    ap.add_argument("--have-rows")
    ap.add_argument("--missing-ids")
    a = ap.parse_args()

    rows, n_dup, n_wl = parse_scan(a.scan)
    short = sorted((v[0], sid, v) for sid, v in rows.items()
                   if v[0] < a.dmax and v[1] < a.vmax and not v[3])
    with open(a.out_promising, "w") as o:
        o.write(f"{len(short)} candidates with gravity_d_pc < {a.dmax} AND rel_speed_kms < {a.vmax}, "
                f"sorted by distance\n\n")
        for d, sid, v in short:
            z = f"  zpt={float(v[4]):+.4f} mas" if v[4] is not None else ""
            o.write(f"{sid:>22}  d={d:9.4f} pc  v={v[1]:7.2f} km/s  t={v[2]:9.3f} Myr{z}\n")
    ids = [sid for _, sid, _ in short]
    print(f"scan rows (dedup): {len(rows)} ({n_dup} duplicate rows skipped), window-limited: {n_wl}")
    print(f"shortlist (d < {a.dmax} pc, v < {a.vmax} km/s): {len(ids)}")
    if a.union_ids:
        old = read_ids(a.union_ids)
        print(f"extra list: {len(old)}; in new shortlist: {len(set(old) & set(ids))}; "
              f"dropped: {len(set(old) - set(ids))}; new entries: {len(set(ids) - set(old))}")
        ids = list(dict.fromkeys(ids + old))
        print(f"union: {len(ids)}")
    open(a.out_ids, "w").write("\n".join(ids) + "\n")
    if a.have_rows:
        have = set()
        with open(a.have_rows, newline="") as f:
            r = csv.reader(f)
            i = next(r).index("source_id")
            for row in r:
                have.add(row[i])
        missing = [s for s in ids if s not in have]
        open(a.missing_ids, "w").write("\n".join(missing) + ("\n" if missing else ""))
        print(f"rows still to extract: {len(missing)}")


if __name__ == "__main__":
    main()
