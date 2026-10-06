// Validates mc_search_stage1()'s Monte Carlo aggregation logic -- NOT
// orbital dynamics (that's test_kepler.c/test_encounter.c's job, already
// validated) or the sampler itself (test_montecarlo.c's job). Uses the
// zero-gravity Kepler potential (GM=0, same as test_encounter.c's own exact
// case) so any discrepancy found here is attributable to mc_search.c's own
// logic, not to orbit-integration error.
#include <stdio.h>
#include <math.h>
#include "../src/potential.h"
#include "../src/mc_search.h"

static int g_failures = 0;

static void check(int ok, const char *what) {
    if (!ok) { printf("  FAIL: %s\n", what); g_failures++; }
    else { printf("  ok:   %s\n", what); }
}

static void test_rejects_bad_covariance(void) {
    printf("test_rejects_bad_covariance:\n");
    KeplerParams kp = { .GM = 0.0 };
    Potential pot = kepler_potential(&kp);
    PhaseSpace iso = { .pos = {15,11,0}, .vel = {0.5, 0.1, 0} };
    PhaseSpace star = { .pos = {10,10,0}, .vel = {1.0, 0.3, 0} };
    // Not positive-definite (same construction as test_montecarlo.c's own bad-matrix case).
    double bad_cov[3][3] = { {1.0, 2.0, 0.0}, {2.0, 1.0, 0.0}, {0.0, 0.0, 1.0} };
    McEncounterSummary s = mc_search_stage1(&pot, iso, bad_cov, star, 20.0, 2000, 42ULL);
    check(s.n_clones_run == 0, "a non-positive-definite covariance matrix is rejected (n_clones_run=0)");
}

static void test_near_deterministic_matches_single_trajectory(void) {
    printf("test_near_deterministic_matches_single_trajectory:\n");
    KeplerParams kp = { .GM = 0.0 };
    Potential pot = kepler_potential(&kp);
    PhaseSpace iso = { .pos = {15,11,0}, .vel = {0.5, 0.1, 0} };
    PhaseSpace star = { .pos = {10,10,0}, .vel = {1.0, 0.3, 0} };

    // A tiny (but still positive-definite) covariance should make every clone's trajectory
    // essentially identical to the unperturbed one -- so the Monte Carlo summary should reduce
    // to the plain single-trajectory answer, which test_encounter.c already validated exactly.
    double tiny_cov[3][3] = { {1e-16,0,0}, {0,1e-16,0}, {0,0,1e-16} };
    McEncounterSummary s = mc_search_stage1(&pot, iso, tiny_cov, star, 20.0, 2000000, 1ULL);
    check(s.n_clones_run == MC_N_ISO_CLONES, "runs all MC_N_ISO_CLONES clones for a valid covariance");

    State a = { .pos = iso.pos, .vel = iso.vel };
    State b = { .pos = star.pos, .vel = star.vel };
    EncounterResult direct = find_closest_approach(&pot, &b, &a, 20.0, 2000000);
    printf("  direct (unperturbed): d=%.9f  mc (tiny cov): min=%.9f mean=%.9f\n",
           direct.min_distance_kpc, s.min_distance_kpc, s.mean_min_distance_kpc);
    check(fabs(s.min_distance_kpc - direct.min_distance_kpc) / direct.min_distance_kpc < 1e-4,
          "with ~zero velocity uncertainty, MC min distance matches the single deterministic trajectory");
    check(fabs(s.mean_min_distance_kpc - direct.min_distance_kpc) / direct.min_distance_kpc < 1e-4,
          "with ~zero velocity uncertainty, MC mean distance also matches (all clones collapse together)");
}

static void test_determinism_and_spread(void) {
    printf("test_determinism_and_spread:\n");
    KeplerParams kp = { .GM = 0.0 };
    Potential pot = kepler_potential(&kp);
    PhaseSpace iso = { .pos = {15,11,0}, .vel = {0.5, 0.1, 0} };
    PhaseSpace star = { .pos = {10,10,0}, .vel = {1.0, 0.3, 0} };
    // A real, non-tiny velocity covariance (kpc/Myr)^2 -- plausible scale
    // for an ISO orbital solution's velocity uncertainty.
    double cov[3][3] = { {1e-4,0,0}, {0,1e-4,0}, {0,0,1e-4} };

    McEncounterSummary s1 = mc_search_stage1(&pot, iso, cov, star, 20.0, 20000, 777ULL);
    McEncounterSummary s2 = mc_search_stage1(&pot, iso, cov, star, 20.0, 20000, 777ULL);
    printf("  same seed: run1 min=%.12f mean=%.12f | run2 min=%.12f mean=%.12f\n",
           s1.min_distance_kpc, s1.mean_min_distance_kpc, s2.min_distance_kpc, s2.mean_min_distance_kpc);
    check(s1.min_distance_kpc == s2.min_distance_kpc && s1.mean_min_distance_kpc == s2.mean_min_distance_kpc,
          "same seed produces bit-identical results (required for BOINC quorum validation)");

    McEncounterSummary s3 = mc_search_stage1(&pot, iso, cov, star, 20.0, 20000, 888ULL);
    printf("  different seed: mean=%.12f (vs %.12f for seed 777)\n", s3.mean_min_distance_kpc, s1.mean_min_distance_kpc);
    check(s1.mean_min_distance_kpc != s3.mean_min_distance_kpc,
          "a different seed produces a different result (confirms randomness is actually used)");

    check(s1.min_distance_kpc <= s1.mean_min_distance_kpc + 1e-12,
          "the minimum across clones is never larger than their mean (basic sanity on the aggregation)");
}

int main(void) {
    test_rejects_bad_covariance();
    test_near_deterministic_matches_single_trajectory();
    test_determinism_and_spread();

    if (g_failures) {
        printf("\n%d check(s) FAILED\n", g_failures);
        return 1;
    }
    printf("\nAll checks passed.\n");
    return 0;
}
