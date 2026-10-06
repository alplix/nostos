// Validates find_closest_approach() against an EXACT analytic case, deliberately decoupled from
// gravity entirely: with GM=0 (see potential.c's Kepler potential, acceleration is identically
// zero when GM=0) both bodies move in straight lines at constant velocity, and the time/distance
// of closest approach between two uniformly-moving points has a textbook closed-form answer
// (minimize |r0 + v_rel*t|^2, a parabola in t) to check the search loop against independent of
// any orbit-integration correctness question (that's test_kepler.c's job).
#include <stdio.h>
#include <math.h>
#include "../src/potential.h"
#include "../src/encounter.h"

static int g_failures = 0;

static void check(int ok, const char *what) {
    if (!ok) { printf("  FAIL: %s\n", what); g_failures++; }
    else { printf("  ok:   %s\n", what); }
}

static void test_zero_gravity_exact(void) {
    printf("test_zero_gravity_exact:\n");
    KeplerParams kp = { .GM = 0.0 };
    Potential pot = kepler_potential(&kp);

    // Both bodies start well away from the origin and keep moving away from it for the whole
    // window (the Kepler potential's acceleration formula divides by r^3, which is 0/0 = NaN
    // exactly AT the origin even though GM=0 makes the force zero everywhere else -- avoid that
    // coordinate singularity rather than relying on it happening to cancel out). Velocities are
    // picked so the analytic closest approach falls partway through the window, not right at an
    // endpoint (which would make the test too easy to pass by accident).
    State a = { .pos = {15.0, 11.0, 0}, .vel = {0.5, 0.1, 0} };
    State b = { .pos = {10.0, 10.0, 0}, .vel = {1.0, 0.3, 0} };

    Vec3 r0 = v3_sub(a.pos, b.pos);
    Vec3 vrel = v3_sub(a.vel, b.vel);
    double t_exact = -v3_dot(r0, vrel) / v3_dot(vrel, vrel);
    Vec3 sep_exact = v3_add(r0, v3_scale(vrel, t_exact));
    double d_exact = v3_norm(sep_exact);
    printf("  analytic: t*=%.6f d*=%.6f\n", t_exact, d_exact);

    int n_steps = 2000000;
    double total_time = 20.0; // comfortably spans t_exact
    EncounterResult r = find_closest_approach(&pot, &a, &b, total_time, n_steps);
    printf("  numeric:  t=%.6f d=%.6f\n", r.time_at_min_myr, r.min_distance_kpc);

    double dt = total_time / n_steps;
    check(fabs(r.time_at_min_myr - t_exact) < 2 * dt, "time of closest approach matches analytic, within one step");
    check(fabs(r.min_distance_kpc - d_exact) / d_exact < 1e-4, "minimum distance matches analytic to < 1e-4 relative");
}

static void test_gravity_sanity(void) {
    printf("test_gravity_sanity:\n");
    KeplerParams kp = { .GM = 1.0 };
    Potential pot = kepler_potential(&kp);

    // Two bodies on the same circular orbit, offset in phase -- their separation oscillates but
    // is bounded; just check the function runs, finds SOME minimum no larger than the starting
    // separation, and reports a time within the integration window. Not an exact check (that's
    // the zero-gravity test's job) -- just confirms nothing breaks once real dynamics are involved.
    double r0 = 2.0;
    double v_circ = sqrt(kp.GM / r0);
    State a = { .pos = {r0, 0, 0}, .vel = {0, v_circ, 0} };
    State b = { .pos = {0, r0, 0}, .vel = {-v_circ, 0, 0} }; // same orbit, 90 degrees ahead

    double d0 = v3_norm(v3_sub(a.pos, b.pos));
    double total_time = 2 * 3.14159265358979323846 * sqrt(r0*r0*r0/kp.GM); // one period
    EncounterResult r = find_closest_approach(&pot, &a, &b, total_time, 200000);
    printf("  d0=%.6f min_d=%.6f at t=%.6f (window=%.6f)\n", d0, r.min_distance_kpc, r.time_at_min_myr, total_time);

    check(r.min_distance_kpc <= d0 + 1e-9, "found minimum is no larger than the starting separation");
    check(r.time_at_min_myr >= 0 && r.time_at_min_myr <= total_time, "reported time falls within the integration window");
}

int main(void) {
    test_zero_gravity_exact();
    test_gravity_sanity();

    if (g_failures) {
        printf("\n%d check(s) FAILED\n", g_failures);
        return 1;
    }
    printf("\nAll checks passed.\n");
    return 0;
}
