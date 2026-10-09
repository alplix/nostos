// Server-side follow-up for the stars the volunteers found (bin/nostos_detail.py runs it).
//
// The volunteer app reports only the closest and mean closest-approach distance and the time of
// the closest clone. For the stars that pass within the threshold (a few thousand per object out of
// 33 million) this tool re-runs exactly the same Monte Carlo (same seed, same 100 ISO clones, same
// integrator) and keeps every clone's result, which gives what the literature reports for a
// candidate (Bailer-Jones et al. 2018/2020): 5% / 50% / 95% points of the closest-approach distance,
// of its time and of the relative speed at closest approach. A slow encounter is the one that can
// matter for an origin; distance alone can't say that.
//
// It also writes the nominal paths (the ISO's mean velocity, the star's catalogue values) relative
// to the Sun, which the 3D view on nostos_3d.php draws.
//
// Input (stdin): a work unit in the volunteer app's format (header + star lines). Output (stdout):
//   I x,y,z;x,y,z;...            ISO path, pc, one point every TRACK_STEP_MYR from now
//   S id dmin dmean d05 d50 d95 t05 t50 t95 v05 v50 v95     pc, Myr (negative = past), km/s
//   P id x,y,z;...               the star's path, same spacing, until a little after the encounter
//   C id x,y,z,vx,vy,vz;...      per clone at its closest approach: the object's position (pc) and
//                                velocity (km/s) relative to the star, Galactic axes (close-up view)
//   G x,y,z;...                  the Sun's own Galactocentric path (pc), for the Galaxy view
//   # done N
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "coords.h"
#include "orbital_elements.h"
#include "potential.h"
#include "encounter.h"
#include "mc_search.h"
#include "montecarlo.h"
#include "nostos_io.h"
#include "units.h"

#define TRACK_STEP_MYR 0.5

static int cmp_double(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : x > y;
}

// q-th quantile of a sorted array, linear interpolation
static double quantile(const double *v, int n, double q) {
    double pos = q * (n - 1);
    int i = (int)pos;
    if (i >= n - 1) return v[n - 1];
    return v[i] + (pos - i) * (v[i + 1] - v[i]);
}

// Positions every `every` steps (including t=0), n_out points; Galactocentric, kpc.
static void sample_path(const Potential *pot, State s, double total_time, int n_steps, int every,
                        int n_out, Vec3 *out) {
    double dt = total_time / n_steps;
    out[0] = s.pos;
    for (int k = 1; k < n_out; k++) {
        for (int i = 0; i < every; i++) leapfrog_step(pot, &s, dt);
        out[k] = s.pos;
    }
}

static void print_path(const Vec3 *p, const Vec3 *sun, int n) {
    for (int k = 0; k < n; k++) {
        printf("%s%.1f,%.1f,%.1f", k ? ";" : "", (p[k].x - sun[k].x) * 1000.0,
               (p[k].y - sun[k].y) * 1000.0, (p[k].z - sun[k].z) * 1000.0);
    }
    printf("\n");
}

int main(void) {
    OrbitalElements el;
    double vel_cov[3][3], total_time, threshold_pc;
    uint64_t seed_base;
    int n_steps;
    StarRecord *stars = NULL;
    size_t n_stars = 0;
    if (!parse_input_stream(stdin, &el, vel_cov, &seed_base, &total_time, &n_steps, &threshold_pc,
                            &stars, &n_stars)) {
        fprintf(stderr, "can't parse input\n");
        return 1;
    }
    PhaseSpace iso;
    if (!orbital_elements_inbound_asymptote_gm(el, NOSTOS_GM_SSB, 0.0, &solar_params_default, &iso)) {
        fprintf(stderr, "invalid orbital elements\n");
        return 1;
    }
    Potential mw = milkyway_potential(&milkyway_potential_default);

    // same Cholesky factor / clone draws as mc_search_stage1()
    double cov6[MC_MAX_DIM][MC_MAX_DIM] = {{0}}, L[MC_MAX_DIM][MC_MAX_DIM];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) cov6[i][j] = vel_cov[i][j];
    if (mc_cholesky(cov6, 3, L) != 0) {
        fprintf(stderr, "velocity covariance is not positive definite\n");
        return 1;
    }
    double mean6[MC_MAX_DIM] = {iso.vel.x, iso.vel.y, iso.vel.z, 0, 0, 0};

    double dt = total_time / n_steps;
    int every = (int)lround(TRACK_STEP_MYR / fabs(dt));
    if (every < 1) every = 1;
    int n_track = n_steps / every + 1;
    Vec3 *sun = malloc(sizeof(Vec3) * n_track), *path = malloc(sizeof(Vec3) * n_track);
    Vec3 zero = {0, 0, 0};
    State sun_state = {helio_equatorial_to_galactocentric_pos(zero, &solar_params_default),
                       helio_equatorial_to_galactocentric_vel(zero, &solar_params_default)};
    sample_path(&mw, sun_state, total_time, n_steps, every, n_track, sun);
    State iso_state = {iso.pos, iso.vel};
    sample_path(&mw, iso_state, total_time, n_steps, every, n_track, path);
    printf("I ");
    print_path(path, sun, n_track);
    printf("G ");
    for (int k = 0; k < n_track; k++) {
        printf("%s%.0f,%.0f,%.0f", k ? ";" : "", sun[k].x * 1000.0, sun[k].y * 1000.0, sun[k].z * 1000.0);
    }
    printf("\n");

    double d[MC_N_ISO_CLONES], t[MC_N_ISO_CLONES], v[MC_N_ISO_CLONES];
    Vec3 rel_pos[MC_N_ISO_CLONES], rel_vel[MC_N_ISO_CLONES];
    for (size_t s = 0; s < n_stars; s++) {
        PhaseSpace star;
        astrometry_to_galactocentric(stars[s].astrometry, &solar_params_default, &star);
        McRng rng;
        mc_rng_seed(&rng, combine_seed(seed_base, stars[s].source_id));
        double dmin = INFINITY, dsum = 0;
        for (int c = 0; c < MC_N_ISO_CLONES; c++) {
            double drawn[MC_MAX_DIM];
            mc_sample_multivariate_gaussian(&rng, mean6, L, 3, drawn);
            State iso_clone = {iso.pos, {drawn[0], drawn[1], drawn[2]}};
            State star_state = {star.pos, star.vel};
            EncounterResult r = find_closest_approach(&mw, &star_state, &iso_clone, total_time, n_steps);
            d[c] = r.min_distance_kpc * 1000.0;
            t[c] = r.time_at_min_myr;
            v[c] = nostos_kpcmyr_to_kms(v3_norm(v3_sub(r.star_vel_at_min, r.iso_vel_at_min)));
            rel_pos[c] = v3_scale(v3_sub(r.iso_pos_at_min, r.star_pos_at_min), 1000.0);
            rel_vel[c] = v3_scale(v3_sub(r.iso_vel_at_min, r.star_vel_at_min), NOSTOS_KMS_PER_KPCMYR);
            dsum += d[c];
            if (d[c] < dmin) dmin = d[c];
        }
        qsort(d, MC_N_ISO_CLONES, sizeof(double), cmp_double);
        qsort(t, MC_N_ISO_CLONES, sizeof(double), cmp_double);
        qsort(v, MC_N_ISO_CLONES, sizeof(double), cmp_double);
        double t50 = quantile(t, MC_N_ISO_CLONES, 0.5);
        printf("S %llu %.6f %.6f %.4f %.4f %.4f %.4f %.4f %.4f %.3f %.3f %.3f\n",
               (unsigned long long)stars[s].source_id, dmin, dsum / MC_N_ISO_CLONES,
               quantile(d, MC_N_ISO_CLONES, 0.05), quantile(d, MC_N_ISO_CLONES, 0.5),
               quantile(d, MC_N_ISO_CLONES, 0.95),
               quantile(t, MC_N_ISO_CLONES, 0.05), t50, quantile(t, MC_N_ISO_CLONES, 0.95),
               quantile(v, MC_N_ISO_CLONES, 0.05), quantile(v, MC_N_ISO_CLONES, 0.5),
               quantile(v, MC_N_ISO_CLONES, 0.95));
        // star path until a little past the encounter (at least 2 Myr, 20% more)
        double t_end = fabs(t50) + fmax(2.0, 0.2 * fabs(t50));
        int n = (int)ceil(t_end / TRACK_STEP_MYR) + 1;
        if (n > n_track) n = n_track;
        State st = {star.pos, star.vel};
        sample_path(&mw, st, total_time, n_steps, every, n, path);
        printf("P %llu ", (unsigned long long)stars[s].source_id);
        print_path(path, sun, n);
        printf("C %llu ", (unsigned long long)stars[s].source_id);
        for (int c = 0; c < MC_N_ISO_CLONES; c++) {
            printf("%s%.2f,%.2f,%.2f,%.1f,%.1f,%.1f", c ? ";" : "", rel_pos[c].x, rel_pos[c].y, rel_pos[c].z,
                   rel_vel[c].x, rel_vel[c].y, rel_vel[c].z);
        }
        printf("\n");
        fflush(stdout);
    }
    printf("# done %zu\n", n_stars);
    free(sun); free(path); free(stars);
    return 0;
}
