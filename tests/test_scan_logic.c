// Validates scan_logic.c's two functions -- pulled out of gaia_scan.c and given real test
// coverage after TWO real bugs were found in that logic by ad-hoc testing alone (a past-vs-future
// sign bug, and a search-window-truncation bug -- see NOTES_method_selection.md for both). These
// tests are specifically shaped to catch a regression of either one.
#include <stdio.h>
#include <math.h>
#include "../src/scan_logic.h"
#include "../src/potential.h"

static int g_failures = 0;

static void check(int ok, const char *what) {
    if (!ok) { printf("  FAIL: %s\n", what); g_failures++; }
    else { printf("  ok:   %s\n", what); }
}

static void test_straight_line_matches_analytic(void) {
    printf("test_straight_line_matches_analytic:\n");
    // Same exact scenario as test_encounter.c's test_zero_gravity_exact -- known analytic answer
    // (t*=9.310345, d*=0.928477), cross-validating straight_line_closest_approach() against an
    // already-independently-validated reference rather than re-deriving a new one.
    PhaseSpace a = { .pos = {15.0, 11.0, 0}, .vel = {0.5, 0.1, 0} };
    PhaseSpace b = { .pos = {10.0, 10.0, 0}, .vel = {1.0, 0.3, 0} };
    double t_star, d_star_pc;
    bool ok = straight_line_closest_approach(a, b, &t_star, &d_star_pc);
    check(ok, "returns true for a non-degenerate pair");
    printf("  t*=%.6f (expect 9.310345), d*=%.6f kpc -> %.3f pc (expect 928.477 pc)\n",
           t_star, d_star_pc / 1000.0, d_star_pc);
    check(fabs(t_star - 9.310345) < 1e-4, "time matches the known analytic value");
    check(fabs(d_star_pc - 928.477) < 1.0, "distance (in pc) matches the known analytic value");
}

static void test_straight_line_degenerate(void) {
    printf("test_straight_line_degenerate:\n");
    // Identical velocities -- vrel=0, the function's own documented degenerate case.
    PhaseSpace a = { .pos = {1,2,3}, .vel = {0.5, 0.5, 0.5} };
    PhaseSpace b = { .pos = {4,5,6}, .vel = {0.5, 0.5, 0.5} };
    double t_star, d_star_pc;
    bool ok = straight_line_closest_approach(a, b, &t_star, &d_star_pc);
    check(!ok, "returns false for identical velocities (would divide by zero otherwise)");
}

static void test_sign_past_vs_future(void) {
    printf("test_sign_past_vs_future:\n");
    // Two scenarios, mirror images in time: one where the closest approach is in the past
    // (t_star<0) and one in the future (t_star>0) -- this is exactly the distinction the real
    // sign bug in gaia_scan.c got wrong (it searched backward regardless of this sign).
    PhaseSpace a_past = { .pos = {15.0, 11.0, 0}, .vel = {0.5, 0.1, 0} };
    PhaseSpace b_past = { .pos = {10.0, 10.0, 0}, .vel = {1.0, 0.3, 0} };
    double t1, d1;
    straight_line_closest_approach(a_past, b_past, &t1, &d1);
    check(t1 > 0, "sanity: this reference scenario's closest approach is forward in time from the start point");

    // Negate both velocities -- time-reversal symmetry means the closest approach should now be
    // at -t1 instead of +t1 (same geometry, motion reversed).
    PhaseSpace a_future = { .pos = a_past.pos, .vel = {-a_past.vel.x, -a_past.vel.y, -a_past.vel.z} };
    PhaseSpace b_future = { .pos = b_past.pos, .vel = {-b_past.vel.x, -b_past.vel.y, -b_past.vel.z} };
    double t2, d2;
    straight_line_closest_approach(a_future, b_future, &t2, &d2);
    printf("  forward case: t=%.6f | time-reversed case: t=%.6f (expect -%.6f)\n", t1, t2, t1);
    check(fabs(t2 - (-t1)) < 1e-9, "time-reversing both velocities flips the sign of t_star exactly");
    check(fabs(d2 - d1) < 1e-9, "distance is unchanged by time-reversal (same geometry)");
}

static void test_gravity_confirm_widens_when_needed(void) {
    printf("test_gravity_confirm_widens_when_needed:\n");
    // Same zero-gravity exact scenario (true t*=9.310345) but with a deliberately too-small
    // initial window (2.0 Myr) -- this is exactly the bug found in gaia_scan.c: without the
    // widening-retry logic, this would silently report the window's own edge (at t=-2.0) instead
    // of the true minimum far outside it.
    KeplerParams kp = { .GM = 0.0 };
    Potential pot = kepler_potential(&kp);
    PhaseSpace iso = { .pos = {15.0, 11.0, 0}, .vel = {0.5, 0.1, 0} };
    PhaseSpace star = { .pos = {10.0, 10.0, 0}, .vel = {1.0, 0.3, 0} };
    // Flip signs so the true closest approach is in the PAST (negative time), matching how
    // gaia_scan.c actually calls this (it only ever searches backward).
    PhaseSpace iso_rev = { .pos = iso.pos, .vel = {-iso.vel.x, -iso.vel.y, -iso.vel.z} };
    PhaseSpace star_rev = { .pos = star.pos, .vel = {-star.vel.x, -star.vel.y, -star.vel.z} };

    double final_window;
    EncounterResult r = gravity_confirm(&pot, star_rev, iso_rev, 2.0, 200000, 4, &final_window);
    printf("  found: t=%.6f d=%.6f (true: t=-9.310345 d=0.928477), final_window=%.2f Myr "
           "(started at 2.0)\n", r.time_at_min_myr, r.min_distance_kpc, final_window);
    check(final_window > 2.0, "the window actually widened beyond its too-small starting value");
    check(fabs(r.time_at_min_myr - (-9.310345)) < 1e-3, "found the TRUE minimum, not the original window's edge");
    check(fabs(r.min_distance_kpc - 0.928477) / 0.928477 < 1e-3, "distance matches the true analytic value");
}

static void test_gravity_confirm_no_widening_needed(void) {
    printf("test_gravity_confirm_no_widening_needed:\n");
    // Same scenario, but this time the initial window is already generous enough -- confirms the
    // function does NOT needlessly widen when it doesn't have to (final_window == initial).
    KeplerParams kp = { .GM = 0.0 };
    Potential pot = kepler_potential(&kp);
    PhaseSpace iso_rev = { .pos = {15.0, 11.0, 0}, .vel = {-0.5, -0.1, 0} };
    PhaseSpace star_rev = { .pos = {10.0, 10.0, 0}, .vel = {-1.0, -0.3, 0} };

    double final_window;
    EncounterResult r = gravity_confirm(&pot, star_rev, iso_rev, 20.0, 200000, 4, &final_window);
    printf("  found: t=%.6f d=%.6f, final_window=%.2f Myr (started at 20.0)\n",
           r.time_at_min_myr, r.min_distance_kpc, final_window);
    check(final_window == 20.0, "window stays at its initial value when already wide enough");
    check(fabs(r.time_at_min_myr - (-9.310345)) < 1e-3, "still finds the correct minimum");
}

int main(void) {
    test_straight_line_matches_analytic();
    test_straight_line_degenerate();
    test_sign_past_vs_future();
    test_gravity_confirm_widens_when_needed();
    test_gravity_confirm_no_widening_needed();

    if (g_failures) {
        printf("\n%d check(s) FAILED\n", g_failures);
        return 1;
    }
    printf("\nAll checks passed.\n");
    return 0;
}
