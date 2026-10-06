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
#include <vector>

#include "boinc_api.h"
#include "filesys.h"

#include "coords.h"
#include "orbital_elements.h"
#include "potential.h"
#include "mc_search.h"
#include "units.h"

struct StarRecord {
    uint64_t source_id;
    Astrometry astrometry;
};

static bool read_input(
    const char *path, OrbitalElements *el, double vel_cov[3][3], uint64_t *seed_base,
    double *total_time_myr, int *n_steps, double *distance_threshold_pc,
    std::vector<StarRecord> *stars
) {
    FILE *f = boinc_fopen(path, "r");
    if (!f) {
        fprintf(stderr, "can't open input file %s\n", path);
        return false;
    }

    double cov00, cov01, cov02, cov11, cov12, cov22;
    if (fscanf(f, "%lf %lf %lf %lf %lf", &el->e, &el->q_au, &el->i_deg, &el->node_deg, &el->peri_deg) != 5) {
        fprintf(stderr, "malformed input file: can't read orbital elements\n");
        fclose(f);
        return false;
    }
    if (fscanf(f, "%lf %lf %lf %lf %lf %lf", &cov00, &cov01, &cov02, &cov11, &cov12, &cov22) != 6) {
        fprintf(stderr, "malformed input file: can't read velocity covariance\n");
        fclose(f);
        return false;
    }
    vel_cov[0][0] = cov00; vel_cov[0][1] = cov01; vel_cov[0][2] = cov02;
    vel_cov[1][0] = cov01; vel_cov[1][1] = cov11; vel_cov[1][2] = cov12;
    vel_cov[2][0] = cov02; vel_cov[2][1] = cov12; vel_cov[2][2] = cov22;

    unsigned long long seed_tmp;
    if (fscanf(f, "%llu", &seed_tmp) != 1) {
        fprintf(stderr, "malformed input file: can't read seed_base\n");
        fclose(f);
        return false;
    }
    *seed_base = (uint64_t)seed_tmp;

    if (fscanf(f, "%lf %d %lf", total_time_myr, n_steps, distance_threshold_pc) != 3) {
        fprintf(stderr, "malformed input file: can't read search parameters\n");
        fclose(f);
        return false;
    }

    StarRecord rec;
    unsigned long long sid;
    while (fscanf(f, "%llu %lf %lf %lf %lf %lf %lf",
                  &sid, &rec.astrometry.ra_deg, &rec.astrometry.dec_deg, &rec.astrometry.parallax_mas,
                  &rec.astrometry.pmra_masyr, &rec.astrometry.pmdec_masyr, &rec.astrometry.rv_kms) == 7) {
        rec.source_id = (uint64_t)sid;
        stars->push_back(rec);
    }

    fclose(f);
    return true;
}

// Mixes a per-work-unit seed with a per-star source_id into one independent-looking 64-bit
// seed -- splitmix64's own mixing step (see montecarlo.c), applied here once up front so
// neighboring source_ids (which, within one Gaia catalog chunk, are often numerically close)
// don't produce visibly-correlated Monte Carlo draws between stars.
static uint64_t combine_seed(uint64_t seed_base, uint64_t source_id) {
    uint64_t z = seed_base ^ (source_id + 0x9E3779B97F4A7C15ULL + (seed_base << 6) + (seed_base >> 2));
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

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
    std::vector<StarRecord> stars;
    if (!read_input(input_path, &el, vel_cov, &seed_base, &total_time_myr, &n_steps,
                     &distance_threshold_pc, &stars)) {
        boinc_finish(1);
        return 1;
    }

    PhaseSpace iso_ps;
    // Evaluated ~9 years after perihelion, NOT at t=0 -- see
    // tests/test_pipeline_oumuamua_hip3757.c's comment and NOTES_method_selection.md for the real
    // bug this choice fixes (t=0/perihelion gives a velocity still mid-deflection by the Sun's
    // gravity, pointing the wrong direction for the long backward integration). The generator
    // that writes this work unit's input file is responsible for having already advanced `el` to
    // represent the ISO's elements at whatever epoch is ~years-not-decades past its own
    // perihelion -- this wrapper just evaluates slightly further forward from THAT epoch to reach
    // the asymptotic regime, same fixed offset every time for consistency.
    double t_eval_myr = 9.08 * 365.25 * NOSTOS_MYR_PER_DAY;
    if (!orbital_elements_to_galactocentric(el, t_eval_myr, &solar_params_default, &iso_ps)) {
        fprintf(stderr, "invalid orbital elements (not hyperbolic, or q<=0)\n");
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
        fprintf(stderr, "resuming from checkpoint: star %zu of %zu\n", start_index, stars.size());
    }

    FILE *outf = boinc_fopen(output_path, start_index > 0 ? "a" : "w");
    if (!outf) {
        fprintf(stderr, "can't open output file %s\n", output_path);
        boinc_finish(1);
        return 1;
    }

    double threshold_kpc = distance_threshold_pc / 1000.0;

    for (size_t i = start_index; i < stars.size(); i++) {
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

        if (!stars.empty()) {
            boinc_fraction_done((double)(i + 1) / (double)stars.size());
        }
    }

    fclose(outf);
    boinc_finish(0);
    return 0;
}
