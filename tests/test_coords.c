// Round-trip self-consistency check for coords.c: astrometry -> Galactocentric phase space ->
// astrometry should recover the original input (to floating-point precision). This is NOT a
// check against any external reference (no internet access when this was written) -- it only
// proves the forward and inverse transforms are mutually consistent, which would catch a sign
// error or a non-orthogonal matrix, but would NOT catch, say, a consistently-applied wrong
// convention. See coords.h's top comment -- cross-checking against astropy/galpy is still a
// real TODO before trusting this for actual star positions.
#include <stdio.h>
#include <math.h>
#include "../src/coords.h"

static int g_failures = 0;

static void check_close(double a, double b, double tol, const char *what) {
    double diff = fabs(a - b);
    double scale = fabs(a) > fabs(b) ? fabs(a) : fabs(b);
    double rel = scale > 1e-12 ? diff / scale : diff;
    if (rel > tol) {
        printf("  FAIL: %s (got %.10f, want %.10f, rel_err=%.3e)\n", what, b, a, rel);
        g_failures++;
    } else {
        printf("  ok:   %s (rel_err=%.3e)\n", what, rel);
    }
}

static void roundtrip_case(const char *label, Astrometry a) {
    printf("%s:\n", label);
    PhaseSpace ps;
    astrometry_to_galactocentric(a, &solar_params_default, &ps);
    Astrometry back = galactocentric_to_astrometry(ps, &solar_params_default);

    check_close(a.ra_deg, back.ra_deg, 1e-9, "ra_deg");
    check_close(a.dec_deg, back.dec_deg, 1e-9, "dec_deg");
    check_close(a.parallax_mas, back.parallax_mas, 1e-9, "parallax_mas");
    check_close(a.pmra_masyr, back.pmra_masyr, 1e-9, "pmra_masyr");
    check_close(a.pmdec_masyr, back.pmdec_masyr, 1e-9, "pmdec_masyr");
    check_close(a.rv_kms, back.rv_kms, 1e-9, "rv_kms");
}

int main(void) {
    // A handful of made-up but physically reasonable cases spanning different sky positions
    // (including near the poles, where ra becomes degenerate/sensitive, and near ra=0/360 wrap).
    roundtrip_case("generic mid-latitude star", (Astrometry){
        .ra_deg = 123.456, .dec_deg = 34.5, .parallax_mas = 5.2,
        .pmra_masyr = 12.3, .pmdec_masyr = -7.8, .rv_kms = 23.4
    });
    roundtrip_case("near ra=0 wrap", (Astrometry){
        .ra_deg = 0.5, .dec_deg = -10.0, .parallax_mas = 2.0,
        .pmra_masyr = -40.0, .pmdec_masyr = 15.0, .rv_kms = -120.0
    });
    roundtrip_case("near north celestial pole", (Astrometry){
        .ra_deg = 200.0, .dec_deg = 89.5, .parallax_mas = 1.0,
        .pmra_masyr = 3.0, .pmdec_masyr = 3.0, .rv_kms = 5.0
    });
    roundtrip_case("near south celestial pole", (Astrometry){
        .ra_deg = 310.0, .dec_deg = -89.5, .parallax_mas = 8.0,
        .pmra_masyr = -3.0, .pmdec_masyr = 2.0, .rv_kms = 50.0
    });
    roundtrip_case("distant, near-zero proper motion", (Astrometry){
        .ra_deg = 270.0, .dec_deg = 0.0, .parallax_mas = 0.3,
        .pmra_masyr = 0.1, .pmdec_masyr = -0.05, .rv_kms = -300.0
    });

    if (g_failures) {
        printf("\n%d check(s) FAILED\n", g_failures);
        return 1;
    }
    printf("\nAll checks passed.\n");
    return 0;
}
