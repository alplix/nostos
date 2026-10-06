"""Post-processes src/gaia_scan.c's raw candidate output into a sorted, annotated summary --
NOT a discovery claim, a triage tool. Flags each candidate's plausibility the same way
Bailer-Jones et al. 2018 do: favoring LOW relative velocity over raw closest-approach distance
(their own most-favored 'Oumuamua candidate, HD 292249, is their LOWEST relative-velocity one, not
their closest-approach one -- see NOTES_method_selection.md), and separately flags anything with an
implausibly high relative speed as a likely data-quality issue rather than a real encounter
candidate at all (normal disk-star relative velocities are tens of km/s, not hundreds).
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
    for p in parts[1:]:
        key, val = p.split("=")
        fields[key] = float(val)
    return source_id, fields


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "D:/nostos_data/scan_oumuamua_candidates.txt"
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            source_id, fields = parse_line(line)
            rows.append((source_id, fields))

    print(f"{len(rows)} candidates total\n")

    plausible = [r for r in rows if r[1]["rel_speed_kms"] < IMPLAUSIBLE_SPEED_KMS]
    flagged = [r for r in rows if r[1]["rel_speed_kms"] >= IMPLAUSIBLE_SPEED_KMS]

    print(f"-- {len(plausible)} plausible-velocity candidates (< {IMPLAUSIBLE_SPEED_KMS} km/s), "
          f"sorted by relative speed (lowest = most plausible origin, per Bailer-Jones et al. 2018) --")
    for source_id, f in sorted(plausible, key=lambda r: r[1]["rel_speed_kms"]):
        print(f"  {source_id:>22s}  v={f['rel_speed_kms']:7.2f} km/s  "
              f"d={f['gravity_d_pc']:6.3f} pc  t={f['gravity_t_myr']:8.3f} Myr")

    print(f"\n-- {len(flagged)} candidates with implausibly high relative speed "
          f"(>= {IMPLAUSIBLE_SPEED_KMS} km/s) -- likely data-quality artifacts, not real "
          f"encounters, not investigated further here --")
    for source_id, f in sorted(flagged, key=lambda r: -r[1]["rel_speed_kms"]):
        print(f"  {source_id:>22s}  v={f['rel_speed_kms']:7.2f} km/s  "
              f"d={f['gravity_d_pc']:6.3f} pc  t={f['gravity_t_myr']:8.3f} Myr")

    print(f"\n-- top 10 by distance alone (gravity-confirmed), for reference --")
    for source_id, f in sorted(rows, key=lambda r: r[1]["gravity_d_pc"])[:10]:
        flag = " [HIGH-V]" if f["rel_speed_kms"] >= IMPLAUSIBLE_SPEED_KMS else ""
        print(f"  {source_id:>22s}  d={f['gravity_d_pc']:6.3f} pc  t={f['gravity_t_myr']:8.3f} Myr  "
              f"v={f['rel_speed_kms']:7.2f} km/s{flag}")


if __name__ == "__main__":
    main()
