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
#include "iso_solutions.h"
#include "parallax_zeropoint.h"
#include "csv_range.h"

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
    // Parallax zero-point inputs (Lindegren et al. 2021; src/parallax_zeropoint.h).
    int gmag, nu_eff, pseudocolour, ecl_lat, params_solved;
} ColumnIndices;

static double field_or_nan(const char *s) {
    if (s[0] == '\0' || strcmp(s, "null") == 0 || strcmp(s, "NULL") == 0) return NAN;
    return atof(s);
}

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

static void usage(const char *prog) {
    fprintf(stderr,
            "usage: %s <gaia_csv_path> [max_rows] [--iso=oumuamua|borisov] "
            "[--zpt=lindegren|global|none] [--byte-range=START:END]\n"
            "  --iso         interstellar object to screen against (default oumuamua)\n"
            "  --zpt         Gaia DR3 parallax zero-point: the full Lindegren et al. 2021 recipe\n"
            "                (default), the global -0.017 mas offset, or none (pre-2026-10-08)\n"
            "  --byte-range  only rows whose first byte lies in [START, END) of the file, so a scan\n"
            "                can be split across processes (the header is always read at offset 0)\n",
            prog);
}

int main(int argc, char **argv) {
    const char *csv_path = NULL;
    long max_rows = -1; // -1 = no limit, scan the whole file
    const char *iso_name = "oumuamua";
    ZptMode zpt_mode = ZPT_MODE_LINDEGREN;
    long long byte_start = 0, byte_end = -1;
    int n_pos = 0;
    for (int i = 1; i < argc; i++) {
        if (strncmp(argv[i], "--iso=", 6) == 0) {
            iso_name = argv[i] + 6;
        } else if (strncmp(argv[i], "--zpt=", 6) == 0) {
            if (gaia_zpt_mode_parse(argv[i] + 6, &zpt_mode) != 0) { usage(argv[0]); return 1; }
        } else if (strncmp(argv[i], "--byte-range=", 13) == 0) {
            if (sscanf(argv[i] + 13, "%lld:%lld", &byte_start, &byte_end) != 2 || byte_start < 0 ||
                byte_end <= byte_start) { usage(argv[0]); return 1; }
        } else if (argv[i][0] == '-' && argv[i][1] == '-') {
            usage(argv[0]); return 1;
        } else if (n_pos == 0) {
            csv_path = argv[i]; n_pos++;
        } else if (n_pos == 1) {
            max_rows = atol(argv[i]); n_pos++;
        } else {
            usage(argv[0]); return 1;
        }
    }
    if (!csv_path) { usage(argv[0]); return 1; }

    FILE *f = fopen(csv_path, "rb"); // binary: --byte-range offsets must be exact byte positions
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
    col.gmag = find_column(headers, n_headers, "phot_g_mean_mag");
    col.nu_eff = find_column(headers, n_headers, "nu_eff_used_in_astrometry");
    col.pseudocolour = find_column(headers, n_headers, "pseudocolour");
    col.ecl_lat = find_column(headers, n_headers, "ecl_lat");
    col.params_solved = find_column(headers, n_headers, "astrometric_params_solved");
    int have_zpt_cols = col.gmag >= 0 && col.nu_eff >= 0 && col.pseudocolour >= 0 &&
                        col.ecl_lat >= 0 && col.params_solved >= 0;
    if (zpt_mode == ZPT_MODE_LINDEGREN && !have_zpt_cols) {
        // Documented fallback (NOTES_method_selection.md, "Gaia DR3 parallax zero-point"): without
        // the five recipe inputs, apply the global QSO-based offset rather than nothing.
        fprintf(stderr, "WARNING: zero-point columns missing (phot_g_mean_mag=%d nu_eff=%d "
                "pseudocolour=%d ecl_lat=%d astrometric_params_solved=%d); using --zpt=global\n",
                col.gmag, col.nu_eff, col.pseudocolour, col.ecl_lat, col.params_solved);
        zpt_mode = ZPT_MODE_GLOBAL;
    }
    int max_col = 0;
    {
        int need[12] = { col.source_id, col.ra, col.dec, col.parallax, col.pmra, col.pmdec, col.rv,
                         col.gmag, col.nu_eff, col.pseudocolour, col.ecl_lat, col.params_solved };
        int n_need = (zpt_mode == ZPT_MODE_LINDEGREN) ? 12 : 7;
        for (int k = 0; k < n_need; k++) if (need[k] > max_col) max_col = need[k];
    }

    // 1I/'Oumuamua's starting state: a point on its exact INBOUND asymptote at the Gaia DR3 epoch
    // (J2016.0), from JPL's barycentric 1600 elements with the total solar-system GM (see
    // src/iso_solutions.h and orbital_elements.h). History, all in NOTES_method_selection.md:
    //  - 2026-10-07 "Found the actual bug": the OUTBOUND (post-perihelion) state was used, whose
    //    direction differs by the full turning angle (112.7 deg) -- distances off 15-200x.
    //  - 2026-10-08 "Residual distance gap": the replacement, the Keplerian state 9.08 yr BEFORE
    //    perihelion, had the right direction but was still 2.3% (0.61 km/s) too fast (vis-viva
    //    converges only like 1/r at ~55 AU) -- the residual ~1.4-2x distance gap; and heliocentric
    //    2017-epoch elements point 0.077 deg away from the true barycentric asymptote.
    // 2I/Borisov (--iso=borisov, 2026-10-08): the same construction from JPL#54's barycentric 1600
    // elements (src/iso_solutions.h), validated against Ross 573 in
    // tests/test_pipeline_borisov_ross573.c.
    OrbitalElements iso_el;
    double t_start_days;
    if (strcmp(iso_name, "oumuamua") == 0) {
        iso_el = OUMUAMUA_ELEMENTS_BARY_1600; t_start_days = OUMUAMUA_T_START_DR3_DAYS;
    } else if (strcmp(iso_name, "borisov") == 0) {
        iso_el = BORISOV_ELEMENTS_BARY_1600; t_start_days = BORISOV_T_START_DR3_DAYS;
    } else {
        fprintf(stderr, "unknown --iso=%s (oumuamua|borisov)\n", iso_name);
        return 1;
    }
    double t_start_myr = t_start_days * NOSTOS_MYR_PER_DAY;
    PhaseSpace iso_ps;
    orbital_elements_inbound_asymptote_gm(iso_el, NOSTOS_GM_SSB, t_start_myr,
                                          &solar_params_default, &iso_ps);
    fprintf(stderr, "iso=%s zpt=%s byte_range=[%lld, %lld)\n", iso_name, gaia_zpt_mode_name(zpt_mode),
            byte_start, byte_end);

    // --byte-range: rows whose FIRST byte lies in [byte_start, byte_end), so adjacent ranges cover
    // every row exactly once (src/csv_range.h, tests/test_csv_range.c).
    if (!csv_range_begin(f, byte_start, line, MAX_LINE)) {
        fprintf(stderr, "nothing at or after offset %lld\n", byte_start);
        return 0;
    }

    Potential mw = milkyway_potential(&milkyway_potential_default);

    long n_rows = 0, n_parsed = 0, n_candidates = 0;
    long n_zpt_out_of_range = 0, n_zpt_fallback = 0;
    while (csv_range_next(f, byte_end, line, MAX_LINE)) {
        n_rows++;
        if (max_rows >= 0 && n_rows > max_rows) break;
        if (n_rows % 2000000 == 0) {
            fprintf(stderr, "  ...%ld rows scanned, %ld candidates so far\n", n_rows, n_candidates);
        }

        char line_copy[MAX_LINE];
        strncpy(line_copy, line, MAX_LINE - 1);
        line_copy[MAX_LINE - 1] = '\0';
        int n_fields = split_line(line_copy, fields, MAX_COLS);
        if (n_fields <= max_col) continue; // malformed/short row, skip

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
        // Gaia DR3 parallax zero-point (2026-10-08, default ON): parallax_corrected = parallax - Z,
        // Z ~ -0.01..-0.05 mas from Lindegren et al. 2021's recipe (src/parallax_zeropoint.h).
        ZptStatus zst = ZPT_OK;
        double zpt;
        if (zpt_mode == ZPT_MODE_LINDEGREN) {
            zpt = gaia_dr3_parallax_zeropoint(field_or_nan(fields[col.gmag]), field_or_nan(fields[col.nu_eff]),
                                              field_or_nan(fields[col.pseudocolour]),
                                              field_or_nan(fields[col.ecl_lat]),
                                              atoi(fields[col.params_solved]), &zst);
        } else {
            zpt = gaia_zeropoint_for_mode(zpt_mode, NAN, NAN, NAN, NAN, 0, &zst);
        }
        parallax -= zpt;
        if (parallax <= 0) continue; // non-physical / unusable distance
        if (zst == ZPT_OUT_OF_RANGE) n_zpt_out_of_range++;
        if (zst == ZPT_FALLBACK) n_zpt_fallback++;

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
            //
            // Second real bug found, in the SAME family, while running the first full scan at a
            // widened 50pc threshold: for LOW relative-velocity candidates specifically (exactly
            // the most scientifically interesting ones, per Bailer-Jones et al. 2018's own
            // preference), the straight-line t_star is a poor predictor of where to center the
            // window at all -- low relative velocity means the straight-line model's own
            // implied timescale is long, which is exactly the regime where gravity's cumulative
            // effect over that long a time can move the true minimum far from a naive 2.5x*t_star
            // guess, no matter how many times that guess gets doubled. Symptom: several
            // "confirmed" distances came back wildly larger than the 50pc straight-line screen
            // that was supposed to be getting refined (e.g. one case: d*=50pc straight-line vs a
            // "confirmed" 5000pc). Fixed by never starting the window below this project's own
            // settled 100 Myr backward-integration trust limit (NOTES_method_selection.md) --
            // cheap to do since this is a single confirmation trajectory, not 100 Monte Carlo
            // clones -- AND by explicitly checking the FINAL result for the same at-the-edge
            // condition gravity_confirm() itself checks internally, flagging it in the output
            // rather than silently trusting a result that still might be window-limited after
            // every retry attempt.
            double initial_window_myr = fabs(t_star) * 2.5 + 2.0;
            if (initial_window_myr < 100.0) initial_window_myr = 100.0;
            double window_myr;
            EncounterResult r = gravity_confirm(&mw, star_ps, iso_ps, initial_window_myr, 2000, 6, &window_myr);
            double rel_speed_kms = nostos_kpcmyr_to_kms(v3_norm(v3_sub(r.star_vel_at_min, r.iso_vel_at_min)));
            double dt_final = window_myr / (int)(window_myr * 2000.0);
            int still_at_edge = (fabs(fabs(r.time_at_min_myr) - window_myr) < 2.0 * dt_final);

            printf("%llu straight_line_d_pc=%.4f straight_line_t_myr=%.4f "
                   "gravity_d_pc=%.4f gravity_t_myr=%.4f rel_speed_kms=%.2f search_window_myr=%.2f "
                   "zpt_mas=%.6f%s%s\n",
                   source_id, d_star_pc, t_star, r.min_distance_kpc * 1000.0, r.time_at_min_myr,
                   rel_speed_kms, window_myr, zpt, still_at_edge ? " WINDOW_LIMITED" : "",
                   zst == ZPT_OUT_OF_RANGE ? " ZPT_OUT_OF_RANGE" : zst == ZPT_FALLBACK ? " ZPT_FALLBACK" : "");
            fflush(stdout);
            n_candidates++;
        }
    }

    fprintf(stderr, "\nDone: %ld rows read, %ld parsed/used, %ld candidates under %.1f pc "
            "(straight-line screen)\n", n_rows, n_parsed, n_candidates, CANDIDATE_THRESHOLD_PC);
    fprintf(stderr, "zero-point (%s): %ld used stars outside the recipe's calibrated G/colour range "
            "(evaluated anyway, flagged ZPT_OUT_OF_RANGE), %ld on the global fallback\n",
            gaia_zpt_mode_name(zpt_mode), n_zpt_out_of_range, n_zpt_fallback);

    fclose(f);
    free(line);
    free(headers_storage);
    return 0;
}
