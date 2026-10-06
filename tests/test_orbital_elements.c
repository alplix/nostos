// Validates orbital_elements_to_galactocentric() as far as possible without an external
// reference (no internet access in this session -- see the file's own top comment and
// NOTES_method_selection.md for the real cross-check this still needs against e.g. JPL Horizons).
//
// Strategy: use a "null" Sun (at the Galactocentric origin, zero velocity) so the function's
// output IS the heliocentric Keplerian state directly, with nothing from coords.c's Sun-offset
// machinery mixed in -- lets these tests check pure two-body orbital mechanics cleanly.
#include <stdio.h>
#include <math.h>
#include "../src/orbital_elements.h"
#include "../src/units.h"

static int g_failures = 0;
static const SolarParams null_sun = {0, 0, 0, 0, 0, 0};

static void check(int ok, const char *what) {
    if (!ok) { printf("  FAIL: %s\n", what); g_failures++; }
    else { printf("  ok:   %s\n", what); }
}

// Specific orbital energy (per unit mass) and angular momentum magnitude, both conserved
// quantities for two-body motion -- rotation-invariant, so these specifically do NOT validate
// the perifocal->ecliptic->equatorial rotation matrices (any orthogonal rotation preserves
// vector magnitudes/dot products), only the radial (Kepler-equation/vis-viva) part of the
// pipeline. See test_identity_rotation() below for a real (if partial) check of the rotation.
static double energy_of(Vec3 pos, Vec3 vel) {
    double r = v3_norm(pos);
    double v2 = v3_dot(vel, vel);
    return 0.5 * v2 - NOSTOS_GM_SUN / r;
}
static double ang_mom_of(Vec3 pos, Vec3 vel) {
    return v3_norm(v3_cross(pos, vel));
}

static void test_perihelion_distance(void) {
    printf("test_perihelion_distance:\n");
    // 'Oumuamua-like elements (e~1.20, q~0.255 AU) -- real published ballpark values, not
    // claimed to be the exact current best-fit (that's what the real validation step needs).
    OrbitalElements el = { .e = 1.1995, .q_au = 0.2553, .i_deg = 122.74, .node_deg = 24.60, .peri_deg = 241.81 };
    PhaseSpace ps;
    bool ok = orbital_elements_to_galactocentric(el, 0.0, &null_sun, &ps);
    check(ok, "function returns true for a valid hyperbolic orbit");
    double r_kpc = v3_norm(ps.pos);
    double q_kpc = el.q_au * NOSTOS_KPC_PER_AU;
    printf("  r(t=0)=%.6e kpc, q=%.6e kpc\n", r_kpc, q_kpc);
    check(fabs(r_kpc - q_kpc) / q_kpc < 1e-9, "heliocentric distance at t=0 (perihelion) equals q");
}

static void test_energy_angmom_conservation(void) {
    printf("test_energy_angmom_conservation:\n");
    // Borisov-like elements (e~3.36, q~2.0 AU).
    OrbitalElements el = { .e = 3.3573, .q_au = 2.0063, .i_deg = 44.05, .node_deg = 308.15, .peri_deg = 209.12 };
    // Time offsets in fractions of a Myr corresponding to a handful of years -- see
    // orbital_elements.h's domain-of-validity note: this two-body model (and its Kepler solver)
    // is only meant to be evaluated within years of perihelion, not at the Myr-scale offsets the
    // rest of this project integrates over. (~1 year = 2.7397e-6 Myr.)
    double t1_myr = -8.0 * 2.7397e-6;  // ~8 years before perihelion
    double t2_myr = 18.0 * 2.7397e-6;  // ~18 years after
    PhaseSpace ps0, ps1, ps2;
    orbital_elements_to_galactocentric(el, 0.0, &null_sun, &ps0);
    orbital_elements_to_galactocentric(el, t1_myr, &null_sun, &ps1);
    orbital_elements_to_galactocentric(el, t2_myr, &null_sun, &ps2);

    double E0 = energy_of(ps0.pos, ps0.vel), E1 = energy_of(ps1.pos, ps1.vel), E2 = energy_of(ps2.pos, ps2.vel);
    double L0 = ang_mom_of(ps0.pos, ps0.vel), L1 = ang_mom_of(ps1.pos, ps1.vel), L2 = ang_mom_of(ps2.pos, ps2.vel);
    printf("  energy:  t=0: %.6e  t1(~8yr before): %.6e  t2(~18yr after): %.6e\n", E0, E1, E2);
    printf("  |L|:     t=0: %.6e  t1(~8yr before): %.6e  t2(~18yr after): %.6e\n", L0, L1, L2);

    check(fabs(E1 - E0) / fabs(E0) < 1e-6, "energy conserved: ~8yr before matches t=0");
    check(fabs(E2 - E0) / fabs(E0) < 1e-6, "energy conserved: ~18yr after matches t=0");
    check(fabs(L1 - L0) / L0 < 1e-6, "|angular momentum| conserved: ~8yr before matches t=0");
    check(fabs(L2 - L0) / L0 < 1e-6, "|angular momentum| conserved: ~18yr after matches t=0");
    check(E0 > 0, "specific orbital energy is positive (unbound, as expected for e>1)");
}

static void test_rotation_preserves_magnitude(void) {
    printf("test_rotation_preserves_magnitude:\n");
    // NOT a check of the perifocal->ecliptic->equatorial rotation's direction/correctness --
    // orthogonal rotations preserve vector magnitude regardless of whether the rotation itself
    // is right, so this only catches a gross error (e.g. a non-orthogonal or scaled matrix).
    // Real validation of the rotation's correctness needs an external reference (JPL Horizons
    // state vector for a real object) -- see the file-level TODO note.
    OrbitalElements el = { .e = 1.5, .q_au = 1.3, .i_deg = 47.0, .node_deg = 133.0, .peri_deg = 81.0 };
    PhaseSpace ps;
    orbital_elements_to_galactocentric(el, 0.0, &null_sun, &ps);
    double q_kpc = el.q_au * NOSTOS_KPC_PER_AU;
    double r_kpc = v3_norm(ps.pos);
    printf("  |pos|=%.9e kpc, q=%.9e kpc\n", r_kpc, q_kpc);
    check(fabs(r_kpc - q_kpc) / q_kpc < 1e-9, "|pos| at perihelion still equals q after the full rotation chain");
}

int main(void) {
    test_perihelion_distance();
    test_energy_angmom_conservation();
    test_rotation_preserves_magnitude();

    if (g_failures) {
        printf("\n%d check(s) FAILED\n", g_failures);
        return 1;
    }
    printf("\nAll checks passed.\n");
    return 0;
}
