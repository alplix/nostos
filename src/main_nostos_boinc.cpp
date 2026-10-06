// Nostos, wrapped as a real BOINC app -- CPU first, matching this project family's usual order
// (see README.md's roadmap). Mirrors the established wrapper shape from the sibling bitboinc-app
// project's own main_range_*_boinc.cpp files (boinc_init/checkpoint/fraction_done), even though
// Nostos is a separate repo -- same infrastructure (README.md: "sharing bitboinc's BOINC
// infrastructure"), same conventions.
//
// One work unit = one interstellar object (ISO) against one CHUNK of Gaia stars. Runs
// mc_search_stage1() (src/mc_search.h -- the exhaustive, cheap Monte Carlo pass described in
// NOTES_method_selection.md's two-stage design) for every star in the chunk, and reports only
// the stars whose closest approach falls under the work unit's own distance threshold -- matching
// this project's stated approach (README.md: "Close encounters below a threshold distance are
// reported back; everything else is discarded to keep result sizes small").
//
// Input file format (plain text, one work unit's data):
//   Line 1: ISO orbital elements and search parameters, space-separated:
//     e q_au i_deg node_deg peri_deg
//     cov00 cov01 cov02 cov11 cov12 cov22   (upper triangle of the 3x3 velocity covariance,
//                                             (kpc/Myr)^2 -- symmetric, so only 6 numbers needed)
//     seed_base  (uint64, XORed with each star's own source_id for a reproducible-but-
//                 per-star-independent Monte Carlo seed -- see mc_search.h's own seed comment)
//     total_time_myr  n_steps  distance_threshold_pc
//   Remaining lines: one star per line, space-separated:
//     source_id ra_deg dec_deg parallax_mas pmra_masyr pmdec_masyr rv_kms
//
// Output file format (plain text): one line per star whose mc_search_stage1() result's
// min_distance_kpc fell under the threshold:
//   source_id min_distance_pc mean_min_distance_pc time_at_min_myr
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>

#include "boinc_api.h"
#include "filesys.h"

#include "coords.h"
#include "orbital_elements.h"
#include "potential.h"
#include "mc_search.h"
#include "units.h"
#include "nostos_io.h"

int main(int argc, char **argv) {
    BOINC_OPTIONS options;
    boinc_options_defaults(options);
    int retval = boinc_init_options(&options);
    if (retval) {
        fprintf(stderr, "boinc_init_options failed: %d\n", retval);
        exit(retval);
    }

    char input_path[512], output_path[512];
    if (boinc_resolve_filename("in", input_path, sizeof(input_path))) {
        fprintf(stderr, "can't resolve input file 'in'\n");
        boinc_finish(1);
        return 1;
    }
    if (boinc_resolve_filename("out", output_path, sizeof(output_path))) {
        fprintf(stderr, "can't resolve output file 'out'\n");
        boinc_finish(1);
        return 1;
    }

    OrbitalElements el;
    double vel_cov[3][3];
    uint64_t seed_base;
    double total_time_myr, distance_threshold_pc;
    int n_steps;
    StarRecord *stars = NULL;
    size_t n_stars = 0;
    FILE *inf = boinc_fopen(input_path, "r");
    bool parsed = inf && parse_input_stream(inf, &el, vel_cov, &seed_base, &total_time_myr,
                                             &n_steps, &distance_threshold_pc, &stars, &n_stars);
    if (inf) fclose(inf);
    if (!parsed) {
        fprintf(stderr, "can't open or parse input file %s\n", input_path);
        boinc_finish(1);
        return 1;
    }

    PhaseSpace iso_ps;
    // Evaluated ~9 years BEFORE perihelion (the INBOUND asymptotic state), not after it and not
    // at t=0 -- see NOTES_method_selection.md's "Found the actual bug" section. t=0/perihelion
    // gives a velocity still mid-deflection by the Sun's gravity; the OUTBOUND (post-perihelion)
    // asymptote, used here until this fix, points in the WRONG direction for this long Milky-Way-
    // only backward integration -- it differs from the correct inbound direction by the orbit's
    // full hyperbolic turning angle (112.7 degrees for 'Oumuamua's own orbit, not a small
    // correction), since this integrator never re-simulates the Sun's actual deflection near
    // perihelion at all. Using the inbound asymptote's direction from the start sidesteps that
    // entirely. The generator that writes this work unit's input file is responsible for having
    // given `el` consistent with that same convention; this wrapper evaluates a further FIXED
    // negative offset from perihelion to make sure the asymptotic regime (not the still-curving
    // near-Sun part of the orbit) is reached, same offset every time for consistency.
    double t_eval_myr = -9.08 * 365.25 * NOSTOS_MYR_PER_DAY;
    if (!orbital_elements_to_galactocentric(el, t_eval_myr, &solar_params_default, &iso_ps)) {
        fprintf(stderr, "invalid orbital elements (not hyperbolic, or q<=0)\n");
        free(stars);
        boinc_finish(1);
        return 1;
    }

    Potential mw = milkyway_potential(&milkyway_potential_default);

    const char *checkpoint_name = "nostos_checkpoint";
    size_t start_index = 0;
    FILE *cp = boinc_fopen(checkpoint_name, "r");
    if (cp) {
        if (fscanf(cp, "%zu", &start_index) != 1) start_index = 0;
        fclose(cp);
        fprintf(stderr, "resuming from checkpoint: star %zu of %zu\n", start_index, n_stars);
    }

    FILE *outf = boinc_fopen(output_path, start_index > 0 ? "a" : "w");
    if (!outf) {
        fprintf(stderr, "can't open output file %s\n", output_path);
        free(stars);
        boinc_finish(1);
        return 1;
    }

    double threshold_kpc = distance_threshold_pc / 1000.0;

    for (size_t i = start_index; i < n_stars; i++) {
        const StarRecord &star = stars[i];
        PhaseSpace star_ps;
        astrometry_to_galactocentric(star.astrometry, &solar_params_default, &star_ps);

        uint64_t seed = combine_seed(seed_base, star.source_id);
        McEncounterSummary summary = mc_search_stage1(
            &mw, iso_ps, vel_cov, star_ps, total_time_myr, n_steps, seed
        );

        if (summary.n_clones_run > 0 && summary.min_distance_kpc < threshold_kpc) {
            fprintf(outf, "%llu %.6f %.6f %.6f\n",
                    (unsigned long long)star.source_id,
                    summary.min_distance_kpc * 1000.0,
                    summary.mean_min_distance_kpc * 1000.0,
                    summary.time_at_min_myr);
        }

        if (boinc_time_to_checkpoint()) {
            FILE *cpw = boinc_fopen(checkpoint_name, "w");
            if (cpw) {
                fprintf(cpw, "%zu\n", i + 1);
                fclose(cpw);
            }
            boinc_checkpoint_completed();
        }

        if (n_stars > 0) {
            boinc_fraction_done((double)(i + 1) / (double)n_stars);
        }
    }

    fclose(outf);
    free(stars);
    boinc_finish(0);
    return 0;
}
