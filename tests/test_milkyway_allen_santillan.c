// Validates the parts of Allen & Santillan (1991)'s Milky Way halo that
// tests/test_milkyway_galpy.c can't reach: there's no stock galpy potential
// class matching its generalized truncated enclosed-mass law (unlike the
// earlier NFW placeholder, which galpy supports directly), so this checks
// internal self-consistency instead -- specifically the two places a sign
// or continuity mistake in potential.c's hand-derived halo energy() formula
// would most likely show up:
//
//   1. Energy and acceleration must agree: a numerical derivative of
//      milkyway_energy() along a radial direction should match
//      milkyway_accel()'s analytic output (accel = -dPhi/dr). This is a
//      real, independent check of the energy formula -- it was integrated
//      by hand from the acceleration law (see potential.c's comment), so a
//      mistake in that integration is exactly a mistake this would catch.
//   2. The halo's two branches (r<Lambda generalized power law, r>=Lambda
//      point-mass exterior) must meet continuously AT r=Lambda, in both
//      acceleration and energy -- a real risk point for any piecewise
//      formula.
//
// Also prints the resulting total rotation curve (same format as the old
// test_milkyway.c) as a final by-eye sanity check that this is still
// recognizably "Milky Way shaped" (~220-235 km/s near the solar radius).
#include <stdio.h>
#include <math.h>
#include "../src/potential.h"
#include "../src/units.h"

static int g_failures = 0;

static void check(int ok, const char *what) {
    if (!ok) { printf("  FAIL: %s\n", what); g_failures++; }
    else { printf("  ok:   %s\n", what); }
}

// Isolate the halo alone (zero the other two components) so this file's
// checks are specifically about the halo's own formula, not the sum.
static Potential halo_only(void) {
    static MilkyWayParams p;
    p = milkyway_potential_default;
    p.bulge_GM = 0;
    p.disk_GM = 0;
    return milkyway_potential(&p);
}

static void test_energy_matches_accel_numerically(void) {
    printf("test_energy_matches_accel_numerically (halo only):\n");
    Potential halo = halo_only();
    // One point well inside Lambda, one well outside -- both branches.
    double test_r[2] = {20.0, 150.0};
    for (int i = 0; i < 2; i++) {
        double r = test_r[i];
        Vec3 pos = {r, 0.0, 0.0};
        double h = 1e-4; // kpc, small relative to r=20..150
        double phi_plus = potential_energy(&halo, (Vec3){r + h, 0.0, 0.0});
        double phi_minus = potential_energy(&halo, (Vec3){r - h, 0.0, 0.0});
        double dphi_dr_numeric = (phi_plus - phi_minus) / (2.0 * h);
        Vec3 a = potential_acceleration(&halo, pos);
        // a.x is the full radial acceleration here (pos.y=pos.z=0), and
        // accel = -dPhi/dr, so dPhi/dr should equal -a.x.
        double rel_err = fabs(dphi_dr_numeric - (-a.x)) / fabs(a.x);
        printf("  r=%.1f kpc: dPhi/dr(numeric)=%.6e  -accel=%.6e  relerr=%.2e\n",
               r, dphi_dr_numeric, -a.x, rel_err);
        char label[64];
        snprintf(label, sizeof(label), "energy/acceleration agree at r=%.1f kpc", r);
        check(rel_err < 1e-6, label);
    }
}

static void test_continuity_at_cutoff(void) {
    printf("test_continuity_at_cutoff:\n");
    Potential halo = halo_only();
    double Lambda = milkyway_potential_default.halo_cutoff;
    double eps = 1e-6; // kpc
    Vec3 pos_in = {Lambda - eps, 0.0, 0.0};
    Vec3 pos_out = {Lambda + eps, 0.0, 0.0};

    double phi_in = potential_energy(&halo, pos_in);
    double phi_out = potential_energy(&halo, pos_out);
    double phi_rel_err = fabs(phi_in - phi_out) / fabs(phi_in);
    printf("  energy: just inside=%.10e  just outside=%.10e  relerr=%.2e\n", phi_in, phi_out, phi_rel_err);
    check(phi_rel_err < 1e-5, "energy is continuous across r=Lambda");

    Vec3 a_in = potential_acceleration(&halo, pos_in);
    Vec3 a_out = potential_acceleration(&halo, pos_out);
    double a_rel_err = fabs(a_in.x - a_out.x) / fabs(a_in.x);
    printf("  accel:  just inside=%.10e  just outside=%.10e  relerr=%.2e\n", a_in.x, a_out.x, a_rel_err);
    check(a_rel_err < 1e-5, "acceleration is continuous across r=Lambda");
}

static void print_rotation_curve(void) {
    printf("Rotation curve (z=0 plane), Allen & Santillan (1991):\n");
    printf("  R (kpc)   v_c (km/s)\n");
    double radii[] = {1.0, 2.0, 4.0, 6.0, 8.0, 8.122, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0};
    Potential mw = milkyway_potential(&milkyway_potential_default);
    for (size_t i = 0; i < sizeof(radii)/sizeof(radii[0]); i++) {
        double R = radii[i];
        Vec3 pos = {R, 0.0, 0.0};
        Vec3 a = potential_acceleration(&mw, pos);
        double a_in = -a.x; // inward-pointing magnitude at this point (y=z=0)
        double v_c_kpcmyr = sqrt(a_in * R);
        double v_c_kms = nostos_kpcmyr_to_kms(v_c_kpcmyr);
        printf("  %7.3f   %9.2f\n", R, v_c_kms);
    }
    printf("\nFor reference: the real Milky Way's rotation curve is close to\n"
           "flat at roughly 220-235 km/s from a few kpc out to beyond the solar\n"
           "radius (~8.1-8.2 kpc) -- this is the model Allen & Santillan (1991)\n"
           "themselves fit to reproduce, so this is a by-eye sanity check, not\n"
           "an independent validation (that's what the checks above are for).\n");
}

int main(void) {
    test_energy_matches_accel_numerically();
    test_continuity_at_cutoff();
    print_rotation_curve();

    if (g_failures) {
        printf("\n%d check(s) FAILED\n", g_failures);
        return 1;
    }
    printf("\nAll checks passed.\n");
    return 0;
}
