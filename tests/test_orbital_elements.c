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
#include "../src/iso_solutions.h"

static int g_failures = 0;
static const SolarParams null_sun = {0, 0, 0, 0, 0, 0};
static const double PI_D = 3.14159265358979323846;

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


// orbital_elements_inbound_asymptote(): the starting state every backward Galactic integration in
// this project now uses (see orbital_elements.h). Checked against the full Keplerian solution far
// out on the inbound branch (where the Sun's pull really is negligible), plus the closed-form
// geometry. Also pins the size of the bug it replaced: the Keplerian state 9.08 yr before
// perihelion -- the old "asymptotic" starting point -- is 2.3% too fast for 'Oumuamua.
static void test_inbound_asymptote(void) {
    printf("test_inbound_asymptote:\n");
    const double yr = 365.25 * NOSTOS_MYR_PER_DAY;
    OrbitalElements cases[2] = {
        { .e = 1.200954455284624, .q_au = 0.2559200355462156, .i_deg = 122.7402198359025,
          .node_deg = 24.59683473134646, .peri_deg = 241.8043518577817 },  // 1I/'Oumuamua
        { .e = 3.356475782676596, .q_au = 2.006520878500843, .i_deg = 44.05264247909138,
          .node_deg = 308.1477292269942, .peri_deg = 209.1236864378081 },  // 2I/Borisov
    };
    for (int c = 0; c < 2; c++) {
        OrbitalElements el = cases[c];
        double a_kpc = el.q_au * NOSTOS_KPC_PER_AU / (el.e - 1.0);
        double v_inf = sqrt(NOSTOS_GM_SUN / a_kpc);

        PhaseSpace asym0, asymT, far_in, far_out, near9;
        double T = -20000.0 * yr; // ~0.5-0.7 million AU out on the inbound leg
        check(orbital_elements_inbound_asymptote(el, 0.0, &null_sun, &asym0), "asymptote: valid orbit accepted");
        orbital_elements_inbound_asymptote(el, T, &null_sun, &asymT);
        orbital_elements_to_galactocentric(el, T, &null_sun, &far_in);
        orbital_elements_to_galactocentric(el, -T, &null_sun, &far_out);
        double old_t_yr = (c == 0) ? -9.08 : -3.0; // the old "asymptotic" evaluation points
        orbital_elements_to_galactocentric(el, old_t_yr * yr, &null_sun, &near9);

        double speed = v3_norm(asym0.vel);
        Vec3 u = v3_scale(asym0.vel, 1.0 / speed);
        double r_far = v3_norm(far_in.pos);
        // vis-viva: the real orbit at r_far must have v^2 = v_inf^2 + 2GM/r exactly.
        double v2_expected = v_inf * v_inf + 2.0 * NOSTOS_GM_SUN / r_far;
        double ang_far = asin(v3_norm(v3_cross(u, far_in.vel)) / v3_norm(far_in.vel));
        double turn = acos(fmax(-1.0, fmin(1.0, v3_dot(u, far_out.vel) / v3_norm(far_out.vel))));
        // perpendicular distance of the real far-out position from the asymptotic line
        Vec3 rel = v3_sub(far_in.pos, asym0.pos);
        Vec3 perp = v3_sub(rel, v3_scale(u, v3_dot(rel, u)));
        double b_kpc = a_kpc * sqrt(el.e * el.e - 1.0);
        double near_excess = v3_norm(near9.vel) / v_inf - 1.0;
        printf("  %s: v_inf=%.4f km/s, angle(asym, Kepler at -20 kyr)=%.2e rad, turning angle=%.4f deg "
               "(2 asin(1/e)=%.4f), |B|/b-1=%.1e, perp. offset at -20 kyr=%.2e AU, "
               "old %.2f yr state too fast by %.2f%%\n",
               c == 0 ? "1I" : "2I", nostos_kpcmyr_to_kms(v_inf), ang_far, turn * 180.0 / PI_D,
               2.0 * asin(1.0 / el.e) * 180.0 / PI_D, v3_norm(asym0.pos) / b_kpc - 1.0,
               v3_norm(perp) * NOSTOS_AU_PER_KPC, old_t_yr, 100.0 * near_excess);

        check(fabs(speed / v_inf - 1.0) < 1e-12, "asymptote: |v| = sqrt(GM/|a|)");
        check(fabs(v3_dot(far_in.vel, far_in.vel) / v2_expected - 1.0) < 1e-9,
              "asymptote: far-out Kepler speed obeys vis-viva with this v_inf");
        check(ang_far < 1e-6, "asymptote: direction matches the Kepler velocity far out on the inbound leg");
        check(fabs(turn - 2.0 * asin(1.0 / el.e)) < 1e-4,
              "asymptote: angle to the outbound direction is the turning angle 2*asin(1/e)");
        check(fabs(v3_norm(asym0.pos) / b_kpc - 1.0) < 1e-12 && fabs(v3_dot(asym0.pos, u)) < 1e-12 * b_kpc * 1e3,
              "asymptote: t=0 point is the impact-parameter point (|B|=b, B perpendicular to v)");
        check(v3_norm(perp) * NOSTOS_AU_PER_KPC < 1e-3,
              "asymptote: the real orbit converges onto this line far out (< 0.001 AU off)");
        check(fabs(v3_norm(v3_sub(asymT.pos, asym0.pos)) / (fabs(T) * v_inf) - 1.0) < 1e-12,
              "asymptote: position advances along the line at v_inf");
        check(near_excess > 0.01,
              "regression: the old near-Sun (-9.08 yr for 1I, -3.0 yr for 2I) starting state is NOT asymptotic (>1% too fast)");
    }
}

// The pipeline's actual ISO starting velocity (src/iso_solutions.h): the inbound asymptote of JPL's
// barycentric 1600 elements with the total solar-system GM must reproduce the asymptote Bailer-Jones
// et al. 2018 published for the same orbit solution (7c = JPL#16; their Table 1: radiant RA 279.4987,
// Dec 33.8649 deg, v_inf 26.391543 km/s, barycentric ICRF) -- an end-to-end external check of the
// asymptote construction AND the frame/GM choice. The heliocentric 2017-epoch elements, by
// contrast, give a radiant 0.077 deg away (RA 279.5752, Dec 33.9079).
static void test_barycentric_asymptote_matches_paper(void) {
    printf("test_barycentric_asymptote_matches_paper:\n");
    const double r2d = 180.0 / PI_D;
    PhaseSpace ps;
    check(orbital_elements_inbound_asymptote_gm(OUMUAMUA_ELEMENTS_BARY_1600, NOSTOS_GM_SSB, 0.0, &null_sun, &ps),
          "barycentric 1600 elements accepted");
    // null_sun output is in Galactic axes; project onto the images of the ICRS basis vectors.
    Vec3 ex = helio_equatorial_to_galactocentric_vel((Vec3){1, 0, 0}, &null_sun);
    Vec3 ey = helio_equatorial_to_galactocentric_vel((Vec3){0, 1, 0}, &null_sun);
    Vec3 ez = helio_equatorial_to_galactocentric_vel((Vec3){0, 0, 1}, &null_sun);
    Vec3 v = { v3_dot(ex, ps.vel), v3_dot(ey, ps.vel), v3_dot(ez, ps.vel) };
    double speed = v3_norm(v);
    double ra = atan2(-v.y, -v.x) * r2d; // radiant = direction it comes FROM
    if (ra < 0) ra += 360.0;
    double dec = asin(-v.z / speed) * r2d;
    double v_kms = nostos_kpcmyr_to_kms(speed);
    printf("  radiant RA=%.5f Dec=%.5f v_inf=%.6f km/s (paper: 279.4987, 33.8649, 26.391543)\n", ra, dec, v_kms);
    check(fabs(ra - 279.4987) < 2e-4 && fabs(dec - 33.8649) < 2e-4,
          "radiant matches Bailer-Jones et al. 2018 Table 1 (7c) to its printed precision");
    // 1e-4 km/s (4e-6 relative): NOSTOS_G itself is only quoted to ~1e-6 (units.h), which alone
    // shifts v_inf by ~2e-5 km/s.
    check(fabs(v_kms - 26.391543) < 1e-4, "v_inf matches Bailer-Jones et al. 2018 Table 1 (7c) to 1e-4 km/s");
}

int main(void) {
    test_perihelion_distance();
    test_energy_angmom_conservation();
    test_rotation_preserves_magnitude();
    test_inbound_asymptote();
    test_barycentric_asymptote_matches_paper();

    if (g_failures) {
        printf("\n%d check(s) FAILED\n", g_failures);
        return 1;
    }
    printf("\nAll checks passed.\n");
    return 0;
}
