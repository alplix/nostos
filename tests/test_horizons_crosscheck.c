// Cross-checks orbital_elements_to_galactocentric()'s full Kepler-solve + perifocal->ecliptic->
// equatorial rotation chain against a REAL external reference: a JPL Horizons state vector for
// 1I/'Oumuamua, fetched live from ssd.jpl.nasa.gov's Horizons API (not a value typed from memory
// or a textbook example). This is the cross-check flagged as a real TODO in orbital_elements.c's
// top comment, coords.h's top comment, and NOTES_method_selection.md -- until now, every other
// test in this project only checked internal self-consistency (energy/angular-momentum
// conservation, round-trip transforms), never an independent source of truth.
//
// Reference data (JPL Horizons, Sun-centered, ICRF frame, AU / AU-day units), both queried for
// the SAME epoch (2026-10-06):
//
//   Osculating elements (COMMAND='A/2017 U1', EPHEM_TYPE=ELEMENTS, epoch JD 2458041.5):
//     e  = 1.200954455284624
//     q  = 0.2559200355462156 AU
//     i  = 122.7402198359025 deg
//     OM = 24.59683473134646 deg  (longitude of ascending node)
//     W  = 241.8043518577817 deg (argument of perihelion)
//     Tp = JD 2458006.000109129120 (time of perihelion passage)
//
//   State vector (same COMMAND, EPHEM_TYPE=VECTORS, same epoch JD 2458041.5 = 2017-10-15 TDB):
//     pos (AU)     = ( 0.9946052552514761,  0.4838636943944900,  0.1409054524559812)
//     vel (AU/day) = ( 0.0256831488176475,  0.0021047028221171,  0.0099543489726916)
//
// Strategy: call orbital_elements_to_galactocentric() with a "null" Sun (same trick as
// test_orbital_elements.c) so its output is the raw heliocentric state in Galactic-aligned axes
// (coords.c's helio_equatorial_to_galactocentric_{pos,vel}() always applies the ICRS->Galactic
// rotation A_G, independent of the Sun's position/velocity) -- then rotate back by A_G's
// transpose (A_G is orthogonal, so its transpose is its inverse) to recover the ICRS-aligned
// state that's directly comparable to the Horizons vector above.
#include <stdio.h>
#include <math.h>
#include "../src/orbital_elements.h"
#include "../src/units.h"

static int g_failures = 0;
static const SolarParams null_sun = {0, 0, 0, 0, 0, 0};

// Same matrix as coords.c's A_G (duplicated here deliberately -- this test's whole point is to
// be an independent check, so it shouldn't share a static function with the code under test;
// see coords.c's own top comment for this constant's provenance).
static const double A_G[3][3] = {
    {-0.0548755604162154, -0.8734370902348850, -0.4838350155487132},
    { 0.4941094278755837, -0.4448296299600112,  0.7469822444972189},
    {-0.8676661490190047, -0.1980763734312015,  0.4559837761750669},
};

static Vec3 mat_transpose_vec(const double m[3][3], Vec3 v) {
    return (Vec3){
        m[0][0]*v.x + m[1][0]*v.y + m[2][0]*v.z,
        m[0][1]*v.x + m[1][1]*v.y + m[2][1]*v.z,
        m[0][2]*v.x + m[1][2]*v.y + m[2][2]*v.z,
    };
}

static void check(int ok, const char *what) {
    if (!ok) { printf("  FAIL: %s\n", what); g_failures++; }
    else { printf("  ok:   %s\n", what); }
}

int main(void) {
    printf("test_horizons_crosscheck:\n");

    OrbitalElements el = {
        .e = 1.200954455284624,
        .q_au = 0.2559200355462156,
        .i_deg = 122.7402198359025,
        .node_deg = 24.59683473134646,
        .peri_deg = 241.8043518577817,
    };
    double t_peri_jd = 2458006.000109129120;
    double t_eval_jd = 2458041.500000000;
    double t_myr = (t_eval_jd - t_peri_jd) * NOSTOS_MYR_PER_DAY;

    PhaseSpace ps;
    bool ok = orbital_elements_to_galactocentric(el, t_myr, &null_sun, &ps);
    check(ok, "function returns true for 'Oumuamua's real hyperbolic elements");

    // ps.pos/ps.vel are in Galactic-aligned axes (A_G already applied) -- rotate back to
    // ICRS-aligned axes for a direct comparison against the Horizons vector.
    Vec3 pos_icrs_kpc = mat_transpose_vec(A_G, ps.pos);
    Vec3 vel_icrs_kpcmyr = mat_transpose_vec(A_G, ps.vel);

    Vec3 pos_au = {
        pos_icrs_kpc.x * NOSTOS_AU_PER_KPC,
        pos_icrs_kpc.y * NOSTOS_AU_PER_KPC,
        pos_icrs_kpc.z * NOSTOS_AU_PER_KPC,
    };
    Vec3 vel_au_day = {
        vel_icrs_kpcmyr.x * NOSTOS_AU_PER_KPC * NOSTOS_MYR_PER_DAY,
        vel_icrs_kpcmyr.y * NOSTOS_AU_PER_KPC * NOSTOS_MYR_PER_DAY,
        vel_icrs_kpcmyr.z * NOSTOS_AU_PER_KPC * NOSTOS_MYR_PER_DAY,
    };

    Vec3 pos_ref = {0.9946052552514761, 0.4838636943944900, 0.1409054524559812};
    Vec3 vel_ref = {0.0256831488176475, 0.0021047028221171, 0.0099543489726916};

    double pos_err = sqrt(
        (pos_au.x-pos_ref.x)*(pos_au.x-pos_ref.x) +
        (pos_au.y-pos_ref.y)*(pos_au.y-pos_ref.y) +
        (pos_au.z-pos_ref.z)*(pos_au.z-pos_ref.z)
    ) / sqrt(pos_ref.x*pos_ref.x + pos_ref.y*pos_ref.y + pos_ref.z*pos_ref.z);
    double vel_err = sqrt(
        (vel_au_day.x-vel_ref.x)*(vel_au_day.x-vel_ref.x) +
        (vel_au_day.y-vel_ref.y)*(vel_au_day.y-vel_ref.y) +
        (vel_au_day.z-vel_ref.z)*(vel_au_day.z-vel_ref.z)
    ) / sqrt(vel_ref.x*vel_ref.x + vel_ref.y*vel_ref.y + vel_ref.z*vel_ref.z);

    printf("  computed pos (AU):     (%.10f, %.10f, %.10f)\n", pos_au.x, pos_au.y, pos_au.z);
    printf("  horizons pos (AU):     (%.10f, %.10f, %.10f)\n", pos_ref.x, pos_ref.y, pos_ref.z);
    printf("  relative pos error:    %.3e\n", pos_err);
    printf("  computed vel (AU/day): (%.10f, %.10f, %.10f)\n", vel_au_day.x, vel_au_day.y, vel_au_day.z);
    printf("  horizons vel (AU/day): (%.10f, %.10f, %.10f)\n", vel_ref.x, vel_ref.y, vel_ref.z);
    printf("  relative vel error:    %.3e\n", vel_err);

    check(pos_err < 1e-4, "position matches JPL Horizons to < 1e-4 relative error");
    check(vel_err < 1e-4, "velocity matches JPL Horizons to < 1e-4 relative error");

    if (g_failures) {
        printf("\n%d check(s) FAILED\n", g_failures);
        return 1;
    }
    printf("\nAll checks passed -- orbital_elements.c's full rotation chain is now independently\n"
           "validated against a real JPL Horizons state vector, not just internal self-consistency.\n");
    return 0;
}
