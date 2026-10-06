// Closes the real external-validation gap coords.h's top comment used to flag: the round-trip
// self-consistency check in test_coords.c can catch a sign error or non-orthogonal matrix, but
// NOT an internally-consistent-but-wrong convention (e.g. a transposed or differently-defined
// A_G rotation matrix would round-trip with itself just fine). This checks
// astrometry_to_galactocentric() against astropy (a real, independent, widely-used
// implementation of the exact same ICRS<->Galactic transform) for a real star -- HIP 3757,
// Bailer-Jones et al. 2018's closest 'Oumuamua-encounter candidate (same star
// tests/test_pipeline_oumuamua_hip3757.c uses).
//
// Reference values from astropy (see scripts/astropy_crosscheck.py), using the SAME real Gaia
// DR2 astrometry: ra=12.05634594268451 deg, dec=-5.1359892365165924 deg,
// parallax=40.20064272792413 mas, pmra=182.16070593050424 mas/yr, pmdec=-120.15964959812227
// mas/yr, rv=26.880687935465836 km/s.
//
// astropy's SkyCoord.transform_to(Galactic()) gives the Galactic-ALIGNED HELIOCENTRIC frame --
// i.e. the same intermediate quantity as coords.c's internal `helio_gal` (rotated by A_G, but
// NOT yet shifted to the Sun's Galactocentric position or had the Sun's own velocity added). This
// test uses a "null" Sun (position/velocity both zero) so astrometry_to_galactocentric()'s raw
// output IS that same intermediate quantity, directly comparable -- same technique
// test_orbital_elements.c uses for its own null-Sun checks.
#include <stdio.h>
#include <math.h>
#include "../src/coords.h"
#include "../src/units.h"

static int g_failures = 0;

static void check_close(double got, double ref, double tol, const char *what) {
    double rel_err = (ref != 0.0) ? fabs(got - ref) / fabs(ref) : fabs(got);
    int ok = rel_err < tol;
    if (!ok) { printf("  FAIL: %s (got=%.6f ref=%.6f relerr=%.2e)\n", what, got, ref, rel_err); g_failures++; }
    else { printf("  ok:   %s (relerr=%.2e)\n", what, rel_err); }
}

int main(void) {
    printf("test_coords_astropy_crosscheck:\n");

    Astrometry hip3757 = {
        .ra_deg = 12.05634594268451,
        .dec_deg = -5.1359892365165924,
        .parallax_mas = 40.20064272792413,
        .pmra_masyr = 182.16070593050424,
        .pmdec_masyr = -120.15964959812227,
        .rv_kms = 26.880687935465836,
    };
    static const SolarParams null_sun = {0, 0, 0, 0, 0, 0};
    PhaseSpace ps;
    astrometry_to_galactocentric(hip3757, &null_sun, &ps);

    // Position: astropy gives pc, convert this project's kpc output to pc for comparison.
    double pos_pc[3] = { ps.pos.x * 1000.0, ps.pos.y * 1000.0, ps.pos.z * 1000.0 };
    double pos_ref_pc[3] = { -4.772120097707575, 8.006357763713293, -23.06299813752978 };
    printf("  computed pos (pc): (%.6f, %.6f, %.6f)\n", pos_pc[0], pos_pc[1], pos_pc[2]);
    printf("  astropy  pos (pc): (%.6f, %.6f, %.6f)\n", pos_ref_pc[0], pos_ref_pc[1], pos_ref_pc[2]);
    check_close(pos_pc[0], pos_ref_pc[0], 1e-3, "position.x matches astropy");
    check_close(pos_pc[1], pos_ref_pc[1], 1e-3, "position.y matches astropy");
    check_close(pos_pc[2], pos_ref_pc[2], 1e-3, "position.z matches astropy");

    // Velocity: astropy gives km/s, convert this project's kpc/Myr output to km/s.
    double vel_kms[3] = {
        nostos_kpcmyr_to_kms(ps.vel.x), nostos_kpcmyr_to_kms(ps.vel.y), nostos_kpcmyr_to_kms(ps.vel.z)
    };
    double vel_ref_kms[3] = { -16.131078667501903, -13.946210645057471, -30.49655736230923 };
    printf("  computed vel (km/s): (%.6f, %.6f, %.6f)\n", vel_kms[0], vel_kms[1], vel_kms[2]);
    printf("  astropy  vel (km/s): (%.6f, %.6f, %.6f)\n", vel_ref_kms[0], vel_ref_kms[1], vel_ref_kms[2]);
    check_close(vel_kms[0], vel_ref_kms[0], 1e-2, "velocity.x matches astropy");
    check_close(vel_kms[1], vel_ref_kms[1], 1e-2, "velocity.y matches astropy");
    check_close(vel_kms[2], vel_ref_kms[2], 1e-2, "velocity.z matches astropy");

    if (g_failures) {
        printf("\n%d check(s) FAILED\n", g_failures);
        return 1;
    }
    printf("\nAll checks passed -- astrometry_to_galactocentric()'s ICRS->Galactic rotation and\n"
           "velocity construction are independently validated against astropy, for a real star,\n"
           "not just round-trip self-consistency.\n");
    return 0;
}
