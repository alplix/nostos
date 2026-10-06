#include "mc_search.h"
#include <math.h>

McEncounterSummary mc_search_stage1(
    const Potential *pot, PhaseSpace iso_mean, const double vel_cov_kpcmyr2[3][3],
    PhaseSpace star, double total_time_myr, int n_steps, uint64_t seed
) {
    McEncounterSummary summary = {0};

    double cov6[MC_MAX_DIM][MC_MAX_DIM] = {{0}};
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            cov6[i][j] = vel_cov_kpcmyr2[i][j];
        }
    }
    double L[MC_MAX_DIM][MC_MAX_DIM];
    if (mc_cholesky(cov6, 3, L) != 0) {
        return summary; // n_clones_run stays 0 -- caller's signal that this input was rejected
    }

    double mean6[MC_MAX_DIM] = {iso_mean.vel.x, iso_mean.vel.y, iso_mean.vel.z, 0, 0, 0};

    McRng rng;
    mc_rng_seed(&rng, seed);

    double min_d = INFINITY;
    double sum_min = 0.0;
    double time_at_min = 0.0;

    for (int c = 0; c < MC_N_ISO_CLONES; c++) {
        double drawn[MC_MAX_DIM];
        mc_sample_multivariate_gaussian(&rng, mean6, L, 3, drawn);

        State iso_clone = {iso_mean.pos, {drawn[0], drawn[1], drawn[2]}};
        State star_state = {star.pos, star.vel};

        EncounterResult r = find_closest_approach(pot, &star_state, &iso_clone, total_time_myr, n_steps);
        sum_min += r.min_distance_kpc;
        if (r.min_distance_kpc < min_d) {
            min_d = r.min_distance_kpc;
            time_at_min = r.time_at_min_myr;
        }
    }

    summary.min_distance_kpc = min_d;
    summary.mean_min_distance_kpc = sum_min / MC_N_ISO_CLONES;
    summary.time_at_min_myr = time_at_min;
    summary.n_clones_run = MC_N_ISO_CLONES;
    return summary;
}
