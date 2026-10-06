"""Post-processes src/gaia_scan.c's raw candidate output into a sorted, annotated summary --
NOT a discovery claim, a triage tool. Flags each candidate's plausibility the same way
Bailer-Jones et al. 2018 do: favoring LOW relative velocity over raw closest-approach distance
(their own most-favored 'Oumuamua candidate, HD 292249, is their LOWEST relative-velocity one, not
their closest-approach one -- see NOTES_method_selection.md), and separately flags anything with an
implausibly high relative speed as a likely data-quality issue rather than a real encounter
candidate at all (normal disk-star relative velocities are tens of km/s, not hundreds). Also
separates out anything gaia_scan.c itself flagged as WINDOW_LIMITED -- its gravity-confirmation
search window still hadn't converged after every retry attempt, so that candidate's reported
distance/time shouldn't be trusted the way a converged one can be (see gaia_scan.c's own comment
on the low-relative-velocity window-sizing bug this flag exists to catch).
"""
import sys

IMPLAUSIBLE_SPEED_KMS = 150.0  # well above normal disk-star kinematics -- see this file's
                                # top comment; not a hard physical limit, just a "treat with
                                # real skepticism, probably a bad radial_velocity measurement"
                                # flag, consistent with how Bailer-Jones et al. 2018 use
                                # relative velocity as a plausibility signal, not distance alone.


def parse_line(line):
    parts = line.split()
    source_id = parts[0]
    fields = {}
    window_limited = False
    for p in parts[1:]:
        if "=" not in p:
            if p == "WINDOW_LIMITED":
                window_limited = True
            continue
        key, val = p.split("=")
        fields[key] = float(val)
    fields["window_limited"] = window_limited
    return source_id, fields


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "D:/nostos_data/scan_oumuamua_candidates.txt"
    rows = []
    seen_source_ids = set()
    n_duplicates = 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            source_id, fields = parse_line(line)
            if source_id in seen_source_ids:
                # A real, if small (~0.4% in one 14M-row preview run), data-quality issue: the
                # same source_id's row appears more than once in gaia_scan.c's output, with
                # IDENTICAL computed values both times -- consistent with the same star's row
                # appearing twice in the downloaded Gaia CSV itself (gaia_scan.c just screens
                # whatever rows the file hands it), not a bug in the screening logic. Root cause
                # not yet tracked down (possibly overlapping/retried chunks in
                # download_filter_6d.py's resume logic) -- deduplicating here, by keeping only
                # the first occurrence, is a pragmatic guard against double-counting a candidate
                # rather than a fix for whatever is producing the duplicate row upstream.
                n_duplicates += 1
                continue
            seen_source_ids.add(source_id)
            rows.append((source_id, fields))

    print(f"{len(rows)} candidates total" +
          (f" ({n_duplicates} duplicate source_id rows skipped)" if n_duplicates else "") + "\n")

    unreliable = [r for r in rows if r[1]["window_limited"]]
    reliable = [r for r in rows if not r[1]["window_limited"]]
    print(f"-- {len(unreliable)} of {len(rows)} candidates are WINDOW_LIMITED (gravity "
          f"confirmation didn't converge -- distance/time NOT trustworthy, excluded below) --\n")

    plausible = [r for r in reliable if r[1]["rel_speed_kms"] < IMPLAUSIBLE_SPEED_KMS]
    flagged = [r for r in reliable if r[1]["rel_speed_kms"] >= IMPLAUSIBLE_SPEED_KMS]

    print(f"-- {len(plausible)} plausible-velocity, window-converged candidates "
          f"(< {IMPLAUSIBLE_SPEED_KMS} km/s), sorted by relative speed (lowest = most plausible "
          f"origin, per Bailer-Jones et al. 2018) --")
    for source_id, f in sorted(plausible, key=lambda r: r[1]["rel_speed_kms"])[:60]:
        print(f"  {source_id:>22s}  v={f['rel_speed_kms']:7.2f} km/s  "
              f"d={f['gravity_d_pc']:9.3f} pc  t={f['gravity_t_myr']:9.3f} Myr")

    print(f"\n-- {len(flagged)} window-converged candidates with implausibly high relative speed "
          f"(>= {IMPLAUSIBLE_SPEED_KMS} km/s) -- likely data-quality artifacts, not real "
          f"encounters, not investigated further here --")
    for source_id, f in sorted(flagged, key=lambda r: -r[1]["rel_speed_kms"])[:20]:
        print(f"  {source_id:>22s}  v={f['rel_speed_kms']:7.2f} km/s  "
              f"d={f['gravity_d_pc']:9.3f} pc  t={f['gravity_t_myr']:9.3f} Myr")

    print(f"\n-- top 10 by distance alone (gravity-confirmed, window-converged only), for reference --")
    for source_id, f in sorted(reliable, key=lambda r: r[1]["gravity_d_pc"])[:10]:
        flag = " [HIGH-V]" if f["rel_speed_kms"] >= IMPLAUSIBLE_SPEED_KMS else ""
        print(f"  {source_id:>22s}  d={f['gravity_d_pc']:9.3f} pc  t={f['gravity_t_myr']:9.3f} Myr  "
              f"v={f['rel_speed_kms']:7.2f} km/s{flag}")


if __name__ == "__main__":
    main()
