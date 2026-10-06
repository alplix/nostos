// Standalone CLI research tool (NOT the BOINC app -- see main_nostos_boinc.cpp for that): streams
// through a downloaded Gaia bulk CSV (scripts/download_filter_6d.py's output,
// D:/nostos_data/gaia_dr3_6d.csv) one line at a time (never loads the whole ~33M-row file into
// memory) and screens every star against a given ISO for a close encounter.
//
// Two-stage screening within this single tool, matching the project's established two-stage
// design (NOTES_method_selection.md) but adapted for a fast LOCAL scan rather than distributed
// BOINC work units:
//   1. Cheap analytic straight-line closest approach (O(1) per star, no integration loop at all --
//      exact textbook formula, same one test_encounter.c's own exact case validates) for EVERY
//      star in the file. This is what actually makes scanning ~33M rows in a single local pass
//      practical.
//   2. Any star whose straight-line closest approach falls under CANDIDATE_THRESHOLD_PC gets
//      re-checked with the REAL encounter search (full Allen & Santillan gravity,
//      find_closest_approach()) for a more trustworthy number -- gravity was found to barely
//      matter for the two published candidates checked so far
//      (tests/test_pipeline_oumuamua_hip3757.c), but there's no reason to assume that holds for
//      every star in the catalog, so cheap-filter-then-confirm is the right shape even though both
//      stages currently tend to agree closely.
//
// No Monte Carlo here (that's mc_search_stage1(), a separate, slower concern for whatever
// shortlist this produces) -- this is the single-best-fit scan meant to find out, for the first
// time, which real Gaia stars are even in the right ballpark at all.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "coords.h"
#include "orbital_elements.h"
#include "encounter.h"
#include "potential.h"
#include "units.h"
#include "scan_logic.h"

#define MAX_LINE 65536
#define MAX_COLS 200
// Widened from an initial 5.0 pc, 2026-10-06: under THIS project's own model, Bailer-Jones et
// al. 2018's own two best-documented candidates (HIP 3757, HD 292249) come out at ~24-27 pc in
// the straight-line screen -- nowhere near their OWN reported 0.6-1.6 pc (the still-unresolved
// discrepancy documented at length in NOTES_method_selection.md's "Pipeline validation" section).
// A 5 pc threshold calibrated assuming agreement with their absolute scale would have silently
// excluded the exact kind of candidate this project's own validated test cases turn out to look
// like under its own model. Widened to a threshold that's actually self-consistent with how this
// project's own numbers behave, not borrowed from a different model's scale.
#define CANDIDATE_THRESHOLD_PC 50.0

typedef struct {
    int source_id, ra, dec, parallax, pmra, pmdec, rv;
} ColumnIndices;

static int find_column(char **headers, int n_headers, const char *name) {
    for (int i = 0; i < n_headers; i++) {
        if (strcmp(headers[i], name) == 0) return i;
    }
    return -1;
}

static int split_line(char *line, char **out_fields, int max_fields) {
    int n = 0;
    char *p = line;
    out_fields[n++] = p;
    while (*p && n < max_fields) {
        if (*p == ',') {
            *p = '\0';
            out_fields[n++] = p + 1;
        }
        p++;
    }
    return n;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s <gaia_csv_path> [max_rows]\n", argv[0]);
        return 1;
    }
    const char *csv_path = argv[1];
    long max_rows = (argc >= 3) ? atol(argv[2]) : -1; // -1 = no limit, scan the whole file

    FILE *f = fopen(csv_path, "r");
    if (!f) {
        fprintf(stderr, "can't open %s\n", csv_path);
        return 1;
    }

    char *line = malloc(MAX_LINE);
    char *fields[MAX_COLS];
    char *headers_storage = malloc(MAX_LINE);
    char *headers[MAX_COLS];

    if (!fgets(line, MAX_LINE, f)) {
        fprintf(stderr, "empty file or read error\n");
        return 1;
    }
    strcpy(headers_storage, line);
    int n_headers = split_line(headers_storage, headers, MAX_COLS);
    for (int i = 0; i < n_headers; i++) {
        size_t len = strlen(headers[i]);
        while (len > 0 && (headers[i][len-1] == '\n' || headers[i][len-1] == '\r')) {
            headers[i][--len] = '\0';
        }
    }

    ColumnIndices col;
    col.source_id = find_column(headers, n_headers, "source_id");
    col.ra = find_column(headers, n_headers, "ra");
    col.dec = find_column(headers, n_headers, "dec");
    col.parallax = find_column(headers, n_headers, "parallax");
    col.pmra = find_column(headers, n_headers, "pmra");
    col.pmdec = find_column(headers, n_headers, "pmdec");
    col.rv = find_column(headers, n_headers, "radial_velocity");
    if (col.source_id < 0 || col.ra < 0 || col.dec < 0 || col.parallax < 0 ||
        col.pmra < 0 || col.pmdec < 0 || col.rv < 0) {
        fprintf(stderr, "couldn't find all required columns in the header "
                "(source_id=%d ra=%d dec=%d parallax=%d pmra=%d pmdec=%d radial_velocity=%d)\n",
                col.source_id, col.ra, col.dec, col.parallax, col.pmra, col.pmdec, col.rv);
        return 1;
    }
    fprintf(stderr, "columns found: source_id=%d ra=%d dec=%d parallax=%d pmra=%d pmdec=%d rv=%d\n",
            col.source_id, col.ra, col.dec, col.parallax, col.pmra, col.pmdec, col.rv);

    // 1I/'Oumuamua, same real JPL Horizons elements as tests/test_horizons_crosscheck.c and
    // tests/test_pipeline_oumuamua_hip3757.c, evaluated at the same ~9-years-post-perihelion
    // epoch that test uses (see its comment for why t=0/perihelion is the wrong choice).
    OrbitalElements oumuamua_el = {
        .e = 1.200954455284624, .q_au = 0.2559200355462156,
        .i_deg = 122.7402198359025, .node_deg = 24.59683473134646, .peri_deg = 241.8043518577817,
    };
    double t_now_myr = 9.08 * 365.25 * NOSTOS_MYR_PER_DAY;
    PhaseSpace iso_ps;
    orbital_elements_to_galactocentric(oumuamua_el, t_now_myr, &solar_params_default, &iso_ps);

    Potential mw = milkyway_potential(&milkyway_potential_default);

    long n_rows = 0, n_parsed = 0, n_candidates = 0;
    while (fgets(line, MAX_LINE, f)) {
        n_rows++;
        if (max_rows >= 0 && n_rows > max_rows) break;
        if (n_rows % 2000000 == 0) {
            fprintf(stderr, "  ...%ld rows scanned, %ld candidates so far\n", n_rows, n_candidates);
        }

        char line_copy[MAX_LINE];
        strncpy(line_copy, line, MAX_LINE - 1);
        line_copy[MAX_LINE - 1] = '\0';
        int n_fields = split_line(line_copy, fields, MAX_COLS);
        if (n_fields <= col.rv) continue; // malformed/short row, skip

        // Any of these being "null" shouldn't happen for ra/dec/parallax/pmra/pmdec/rv in this
        // already-RV-filtered file, but check defensively rather than crash on a bad row.
        int bad = 0;
        for (int c = 0; c < 7; c++) {
            int idx = (c==0?col.source_id:c==1?col.ra:c==2?col.dec:c==3?col.parallax:c==4?col.pmra:c==5?col.pmdec:col.rv);
            if (strcmp(fields[idx], "null") == 0 || fields[idx][0] == '\0') { bad = 1; break; }
        }
        if (bad) continue;

        unsigned long long source_id = strtoull(fields[col.source_id], NULL, 10);
        double ra = atof(fields[col.ra]);
        double dec = atof(fields[col.dec]);
        double parallax = atof(fields[col.parallax]);
        double pmra = atof(fields[col.pmra]);
        double pmdec = atof(fields[col.pmdec]);
        double rv = atof(fields[col.rv]);
        if (parallax <= 0) continue; // non-physical / unusable distance

        n_parsed++;

        Astrometry a = { ra, dec, parallax, pmra, pmdec, rv };
        PhaseSpace star_ps;
        astrometry_to_galactocentric(a, &solar_params_default, &star_ps);

        double t_star, d_star_pc;
        if (!straight_line_closest_approach(iso_ps, star_ps, &t_star, &d_star_pc)) continue; // degenerate, never happens for real data

        // Only PAST encounters matter for an origin search (t_star>0 means the straight-line
        // closest approach is in the future, irrelevant to "where did the ISO come from" -- real
        // bug found while building this: the gravity-confirmation step below always searched
        // backward regardless of t_star's sign, so a future-candidate star got "confirmed" against
        // an unrelated backward trajectory and produced a wildly wrong, meaningless distance).
        if (t_star < 0 && d_star_pc < CANDIDATE_THRESHOLD_PC) {
            // Confirm with the real encounter search (full gravity). The straight-line estimate's
            // OWN time is only a straight-line-model prediction of where to center the search
            // window -- over tens of Myr (this project's backward-integration window goes out to
            // 100 Myr, see NOTES_method_selection.md), gravity's cumulative effect can shift the
            // true minimum meaningfully away from that naive guess, which is part of why the
            // gravity confirmation exists at all. Real bug found while building this: the first
            // window (1.5x the straight-line t_star) was often too narrow for large |t_star|,
            // silently reporting the window's own EDGE (not a true minimum) as the answer for
            // several candidates -- visible as a far LARGER "confirmed" distance than the
            // straight-line screen found, the opposite of what a real confirmation should ever do.
            // Fixed with a generous, still-bounded widening (2.5x + a larger fixed pad) AND an
            // explicit edge check that re-runs with a wider window if the result lands within one
            // step of either boundary, rather than trusting a result that might just be truncated.
            // Both this and the straight-line screen above were pulled out into scan_logic.c/.h
            // after these two bugs were found by ad-hoc testing with zero automated coverage --
            // see tests/test_scan_logic.c, which now exercises both independently.
            double initial_window_myr = fabs(t_star) * 2.5 + 2.0;
            double window_myr;
            EncounterResult r = gravity_confirm(&mw, star_ps, iso_ps, initial_window_myr, 2000, 4, &window_myr);
            double rel_speed_kms = nostos_kpcmyr_to_kms(v3_norm(v3_sub(r.star_vel_at_min, r.iso_vel_at_min)));

            printf("%llu straight_line_d_pc=%.4f straight_line_t_myr=%.4f "
                   "gravity_d_pc=%.4f gravity_t_myr=%.4f rel_speed_kms=%.2f search_window_myr=%.2f\n",
                   source_id, d_star_pc, t_star, r.min_distance_kpc * 1000.0, r.time_at_min_myr,
                   rel_speed_kms, window_myr);
            fflush(stdout);
            n_candidates++;
        }
    }

    fprintf(stderr, "\nDone: %ld rows read, %ld parsed/used, %ld candidates under %.1f pc "
            "(straight-line screen)\n", n_rows, n_parsed, n_candidates, CANDIDATE_THRESHOLD_PC);

    fclose(f);
    free(line);
    free(headers_storage);
    return 0;
}
