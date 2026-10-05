// Validates the leapfrog integrator itself, completely independent of any
// Milky-Way-potential parameter uncertainty: a Kepler (point-mass) orbit has
// an exact analytic solution and exact conservation laws to check against.
//
//   1. Energy and z-angular-momentum stay bounded (don't secularly drift)
//      over many orbital periods.
//   2. A circular orbit returns to within a small tolerance of its starting
//      position after exactly one period.
//   3. An eccentric orbit's position at several times matches the classic
//      Kepler's-equation analytic solution.
//
// If any of these fail, nothing built on top of integrator.h can be trusted,
// so this must pass before potential.c's Milky-Way model is worth testing.
#include <stdio.h>
#include <math.h>
#include "../src/potential.h"
#include "../src/integrator.h"

#define TEST_PI 3.14159265358979323846

static int g_failures = 0;

static void check(int ok, const char *what) {
    if (!ok) {
        printf("  FAIL: %s\n", what);
        g_failures++;
    } else {
        printf("  ok:   %s\n", what);
    }
}

// --- Test 1: energy + angular momentum conservation over many periods ----
static void test_conservation(void) {
    printf("test_conservation:\n");
    KeplerParams kp = { .GM = 1.0 };
    Potential pot = kepler_potential(&kp);

    double r0 = 2.0;
    double v_circ = sqrt(kp.GM / r0);
    State s = { .pos = {r0, 0, 0}, .vel = {0, v_circ * 0.6, v_circ * 0.3} };
    // (deliberately NOT a perfect circle -- an eccentric orbit stresses the
    // integrator's conservation properties more than a circle would)

    double r_for_E = v3_norm(s.pos);
    double v2 = v3_dot(s.vel, s.vel);
    double E0 = 0.5 * v2 + potential_energy(&pot, s.pos);
    Vec3 L0v = v3_cross(s.pos, s.vel);
    double Lz0 = L0v.z;
    (void)r_for_E;

    double period = 2 * TEST_PI * sqrt(r0 * r0 * r0 / kp.GM) * 3.0; // a few periods' worth of time unit
    int n_steps = 200000;
    double dt = period / n_steps;

    double max_dE = 0, max_dLz = 0;
    for (int i = 0; i < n_steps; i++) {
        leapfrog_step(&pot, &s, dt);
        double v2i = v3_dot(s.vel, s.vel);
        double Ei = 0.5 * v2i + potential_energy(&pot, s.pos);
        Vec3 Liv = v3_cross(s.pos, s.vel);
        double dE = fabs((Ei - E0) / E0);
        double dLz = fabs((Liv.z - Lz0) / Lz0);
        if (dE > max_dE) max_dE = dE;
        if (dLz > max_dLz) max_dLz = dLz;
    }
    printf("  max relative energy drift: %.3e\n", max_dE);
    printf("  max relative Lz drift:     %.3e\n", max_dLz);
    check(max_dE < 1e-4, "energy conserved to < 1e-4 relative over 3 periods");
    check(max_dLz < 1e-8, "Lz conserved to < 1e-8 relative (exact for an axisymmetric force law)");
}

// --- Test 2: circular orbit returns to start after one period ------------
static void test_circular_period(void) {
    printf("test_circular_period:\n");
    KeplerParams kp = { .GM = 1.0 };
    Potential pot = kepler_potential(&kp);

    double r0 = 3.0;
    double v_circ = sqrt(kp.GM / r0);
    State s = { .pos = {r0, 0, 0}, .vel = {0, v_circ, 0} };

    double period = 2 * TEST_PI * sqrt(r0 * r0 * r0 / kp.GM);
    int n_steps = 20000;
    integrate(&pot, &s, period, n_steps, NULL, NULL);

    double dx = s.pos.x - r0, dy = s.pos.y - 0, dz = s.pos.z - 0;
    double drift = sqrt(dx * dx + dy * dy + dz * dz);
    printf("  position drift after 1 period: %.3e kpc (r0=%.1f)\n", drift, r0);
    check(drift / r0 < 1e-5, "circular orbit returns to start within 1e-5*r0 after one period");
}

// --- Test 3: eccentric orbit vs. the analytic Kepler-equation solution ---
static double solve_kepler_E(double M, double e) {
    double E = M; // initial guess
    for (int i = 0; i < 100; i++) {
        double f = E - e * sin(E) - M;
        double fp = 1 - e * cos(E);
        double dE = f / fp;
        E -= dE;
        if (fabs(dE) < 1e-14) break;
    }
    return E;
}

static void test_eccentric_vs_analytic(void) {
    printf("test_eccentric_vs_analytic:\n");
    double GM = 1.0;
    double a = 2.0, e = 0.5;
    KeplerParams kp = { .GM = GM };
    Potential pot = kepler_potential(&kp);

    // Start exactly at periapsis: r_p = a(1-e), moving purely in +y, with
    // vis-viva speed v_p = sqrt(GM*(1+e)/(a*(1-e))).
    double r_p = a * (1 - e);
    double v_p = sqrt(GM * (1 + e) / (a * (1 - e)));
    State s = { .pos = {r_p, 0, 0}, .vel = {0, v_p, 0} };

    double n = sqrt(GM / (a * a * a)); // mean motion
    double period = 2 * TEST_PI / n;

    int n_steps = 500000;
    double dt = period / n_steps;

    double max_rel_err = 0;
    int check_points[] = {1, n_steps / 8, n_steps / 4, n_steps / 2, (n_steps * 3) / 4};
    int next_check = 0;
    for (int i = 1; i <= n_steps && next_check < 5; i++) {
        leapfrog_step(&pot, &s, dt);
        if (i == check_points[next_check]) {
            double t = i * dt;
            double M = n * t;
            double E = solve_kepler_E(M, e);
            double x_analytic = a * (cos(E) - e);
            double y_analytic = a * sqrt(1 - e * e) * sin(E);
            double dx = s.pos.x - x_analytic, dy = s.pos.y - y_analytic;
            double err = sqrt(dx * dx + dy * dy);
            double rel_err = err / a;
            printf("  t=%.4f: integrator=(%.6f,%.6f) analytic=(%.6f,%.6f) rel_err=%.3e\n",
                   t, s.pos.x, s.pos.y, x_analytic, y_analytic, rel_err);
            if (rel_err > max_rel_err) max_rel_err = rel_err;
            next_check++;
        }
    }
    check(max_rel_err < 1e-4, "eccentric-orbit position matches analytic Kepler solution to < 1e-4*a");
}

int main(void) {
    test_conservation();
    test_circular_period();
    test_eccentric_vs_analytic();

    if (g_failures) {
        printf("\n%d check(s) FAILED\n", g_failures);
        return 1;
    }
    printf("\nAll checks passed.\n");
    return 0;
}
