#include "mc_search.h"
#include <math.h>
#include <stdlib.h>

// The two helpers below exist for speed only and leave every result bit-identical (2026-11-04, docs/BENCHMARK.md,
// "Results"). They live here and not in integrator.h or encounter.c because those files are copies shared with
// phygas (NOS-020), and the app's only path to them is this function.

// leapfrog_step() of integrator.h with the acceleration carried from one step to the next: the a1 a step computes at
// its new position is the a0 the next step would compute at that same position, and the acceleration is a pure
// function of the position, so passing it on saves one of the two calls per step and changes no bit. The arithmetic
// is written exactly as in leapfrog_step(), so a compiler that fuses multiply-adds fuses the same ones.
static inline void leapfrog_step_carry(const Potential *pot, State *state, double dt, Vec3 *accel) {
    Vec3 v_half = v3_add(state->vel, v3_scale(*accel, dt * 0.5));
    state->pos = v3_add(state->pos, v3_scale(v_half, dt));
    *accel = potential_acceleration(pot, state->pos);
    state->vel = v3_add(v_half, v3_scale(*accel, dt * 0.5));
}

// find_closest_approach() for one clone against a star trajectory computed beforehand: the star's state after each
// step is the same for all MC_N_ISO_CLONES clones (the bodies do not attract each other), so it is integrated once
// per star instead of once per clone. Same dt, same separations, same comparisons as find_closest_approach(); only
// the distance and the time of the minimum are kept, which is all mc_search_stage1() reads.
static void closest_approach_to_track(
    const Potential *pot, Vec3 star_pos0, const Vec3 *star_track, State iso, double total_time, int n_steps,
    double *min_d_out, double *t_min_out
) {
    double dt = total_time / n_steps;
    double min_d = v3_norm(v3_sub(star_pos0, iso.pos));
    double t_min = 0.0;
    Vec3 a = potential_acceleration(pot, iso.pos);
    for (int i = 0; i < n_steps; i++) {
        leapfrog_step_carry(pot, &iso, dt, &a);
        double d = v3_norm(v3_sub(star_track[i], iso.pos));
        if (d < min_d) {
            min_d = d;
            t_min = (i + 1) * dt;
        }
    }
    *min_d_out = min_d;
    *t_min_out = t_min;
}

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

    // The star's positions after each step, shared by all clones (closest_approach_to_track()): 4 acceleration calls
    // per clone-step become about 1. Measured 2026-11-04 (docs/BENCHMARK.md, Results): 2.8x faster on an Apple M4,
    // 2.7x on Linux arm64 with GCC, the same 17-digit output. 24 bytes per step (2.4 MB at the usual 100,000 steps).
    // Without the memory the old per-clone path below runs, with the same results.
    Vec3 *star_track = n_steps > 0 ? (Vec3 *)malloc(sizeof(Vec3) * (size_t)n_steps) : NULL;
    if (star_track) {
        State s = {star.pos, star.vel};
        double dt = total_time_myr / n_steps;
        Vec3 a = potential_acceleration(pot, s.pos);
        for (int i = 0; i < n_steps; i++) {
            leapfrog_step_carry(pot, &s, dt, &a);
            star_track[i] = s.pos;
        }
    }

    for (int c = 0; c < MC_N_ISO_CLONES; c++) {
        double drawn[MC_MAX_DIM];
        mc_sample_multivariate_gaussian(&rng, mean6, L, 3, drawn);

        State iso_clone = {iso_mean.pos, {drawn[0], drawn[1], drawn[2]}};
        double d_c, t_c;
        if (star_track) {
            closest_approach_to_track(pot, star.pos, star_track, iso_clone, total_time_myr, n_steps, &d_c, &t_c);
        } else {
            State star_state = {star.pos, star.vel};
            EncounterResult r = find_closest_approach(pot, &star_state, &iso_clone, total_time_myr, n_steps);
            d_c = r.min_distance_kpc;
            t_c = r.time_at_min_myr;
        }
        sum_min += d_c;
        if (d_c < min_d) {
            min_d = d_c;
            time_at_min = t_c;
        }
    }
    free(star_track);

    summary.min_distance_kpc = min_d;
    summary.mean_min_distance_kpc = sum_min / MC_N_ISO_CLONES;
    summary.time_at_min_myr = time_at_min;
    summary.n_clones_run = MC_N_ISO_CLONES;
    return summary;
}
