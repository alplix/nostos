// Pins src/parallax_zeropoint.c (C port of the Gaia collaboration's gaiadr3_zeropoint 0.1.0,
// Lindegren et al. 2021) against reference values produced by the OFFICIAL Python package itself
// (scripts/zeropoint_crosscheck.py; the package rounds to 1e-6 mas, hence the tolerance). Cases
// cover 5p and 6p solutions, every G interval of the coefficient tables including the clamps below
// G=6 and above G=21, colours inside and outside the calibrated range, both ecliptic hemispheres,
// and three real DR3 'Oumuamua-shortlist stars (HIP 3757 among them). Also checks the status
// flags, the global fallback and the --zpt mode parsing.
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../src/parallax_zeropoint.h"

static int g_failures = 0;

static void check(int ok, const char *what) {
    if (!ok) { printf("  FAIL: %s\n", what); g_failures++; }
    else printf("  ok:   %s\n", what);
}

typedef struct { double g, nu, pc, b; int ps; double z_ref; } Case;

int main(void) {
    printf("test_parallax_zeropoint:\n");
    const double nan = NAN;
    static const Case cases[] = {
        { 5.0, 1.5, NAN, 10.0, 31, -0.03746 },
        { 8.3, 1.55, NAN, -35.0, 31, -0.023525 },
        { 10.95, 1.4, NAN, 60.0, 31, -0.025763 },
        { 11.5, 1.62, NAN, -80.0, 31, -0.030354 },
        { 12.5, 1.3, NAN, 5.0, 31, -0.011411 },
        { 13.0, 1.75, NAN, 45.0, 31, -0.029943 },
        { 14.2, 1.48, NAN, -20.0, 31, -0.042854 },
        { 16.0, 1.2, NAN, 70.0, 31, -0.052891 },
        { 17.0, 1.85, NAN, -50.0, 31, 0.001765 },
        { 18.3, 1.58, NAN, 0.0, 31, -0.022657 },
        { 19.5, 1.45, NAN, 30.0, 31, -0.015368 },
        { 20.5, 1.35, NAN, -65.0, 31, 0.000517 },
        { 21.5, 1.5, NAN, 15.0, 31, -0.016563 },
        { 15.5, 2.0, NAN, 25.0, 31, 0.004561 },
        { 9.0, NAN, 1.5, 40.0, 95, -0.029822 },
        { 12.0, NAN, 1.3, -10.0, 95, -0.012767 },
        { 14.0, NAN, 1.65, 75.0, 95, -0.024085 },
        { 16.8, NAN, 1.2, -40.0, 95, -0.059035 },
        { 19.2, NAN, 1.8, 55.0, 95, -0.005032 },
        { 20.7, NAN, 1.45, -5.0, 95, -0.017811 },
        // real Gaia DR3 stars
        { 11.169264, 1.3237802, NAN, -9.490552045804828, 31, -0.024473 }, // 2525688198820543360 = HIP 3757
        { 12.373047, 1.4066502, NAN, 2.3069070166060857, 31, -0.011996 }, // 2618904962125653888
        { 11.13299, 1.4818748, NAN, 16.123400126834508, 31, -0.035144 },  // 4399766976222478336
    };
    int n = (int)(sizeof cases / sizeof cases[0]);
    double max_err = 0.0;
    for (int i = 0; i < n; i++) {
        ZptStatus st;
        double z = gaia_dr3_parallax_zeropoint(cases[i].g, cases[i].nu, cases[i].pc, cases[i].b,
                                               cases[i].ps, &st);
        double err = fabs(z - cases[i].z_ref);
        if (err > max_err) max_err = err;
        if (err > 6e-7) {
            printf("  FAIL: case %d G=%.3f: got %.7f, official package %.6f\n", i, cases[i].g, z,
                   cases[i].z_ref);
            g_failures++;
        }
    }
    char msg[160];
    snprintf(msg, sizeof msg, "%d cases match the official gaiadr3_zeropoint package (max |diff| %.1e mas)",
             n, max_err);
    check(max_err <= 6e-7, msg);

    ZptStatus st;
    gaia_dr3_parallax_zeropoint(14.2, 1.48, nan, -20.0, 31, &st);
    check(st == ZPT_OK, "in-range 5p star -> ZPT_OK");
    gaia_dr3_parallax_zeropoint(5.0, 1.5, nan, 10.0, 31, &st);
    check(st == ZPT_OUT_OF_RANGE, "G < 6 -> ZPT_OUT_OF_RANGE");
    gaia_dr3_parallax_zeropoint(15.5, 2.0, nan, 25.0, 31, &st);
    check(st == ZPT_OUT_OF_RANGE, "nu_eff > 1.9 -> ZPT_OUT_OF_RANGE");
    gaia_dr3_parallax_zeropoint(16.8, nan, 1.2, -40.0, 95, &st);
    check(st == ZPT_OUT_OF_RANGE, "pseudocolour < 1.24 -> ZPT_OUT_OF_RANGE");
    // The 6p recipe ignores nu_eff, the 5p one ignores pseudocolour.
    double z6a = gaia_dr3_parallax_zeropoint(12.0, 1.1, 1.3, -10.0, 95, NULL);
    double z6b = gaia_dr3_parallax_zeropoint(12.0, 1.7, 1.3, -10.0, 95, NULL);
    check(z6a == z6b, "6p result independent of nu_eff");

    double zf = gaia_dr3_parallax_zeropoint(14.0, nan, 1.5, 20.0, 31, &st);
    check(st == ZPT_FALLBACK && zf == NOSTOS_ZPT_GLOBAL_MAS, "5p star with null nu_eff -> global fallback");
    zf = gaia_dr3_parallax_zeropoint(14.0, 1.5, 1.5, 20.0, 3, &st);
    check(st == ZPT_FALLBACK && zf == NOSTOS_ZPT_GLOBAL_MAS, "2-parameter solution -> global fallback");
    zf = gaia_dr3_parallax_zeropoint(nan, 1.5, nan, 20.0, 31, &st);
    check(st == ZPT_FALLBACK, "null G -> global fallback");

    ZptMode m;
    check(gaia_zpt_mode_parse("lindegren", &m) == 0 && m == ZPT_MODE_LINDEGREN, "parse 'lindegren'");
    check(gaia_zpt_mode_parse("global", &m) == 0 && m == ZPT_MODE_GLOBAL, "parse 'global'");
    check(gaia_zpt_mode_parse("none", &m) == 0 && m == ZPT_MODE_NONE, "parse 'none'");
    check(gaia_zpt_mode_parse("bogus", &m) != 0, "reject unknown mode");
    check(gaia_zeropoint_for_mode(ZPT_MODE_NONE, 14.2, 1.48, nan, -20.0, 31, NULL) == 0.0, "mode none -> 0");
    check(gaia_zeropoint_for_mode(ZPT_MODE_GLOBAL, 14.2, 1.48, nan, -20.0, 31, NULL) == NOSTOS_ZPT_GLOBAL_MAS,
          "mode global -> -0.017 mas");
    check(fabs(gaia_zeropoint_for_mode(ZPT_MODE_LINDEGREN, 14.2, 1.48, nan, -20.0, 31, NULL) + 0.042854) < 6e-7,
          "mode lindegren -> full recipe");
    check(strcmp(gaia_zpt_mode_name(ZPT_MODE_LINDEGREN), "lindegren") == 0, "mode name round trip");

    if (g_failures) { printf("test_parallax_zeropoint: %d FAILURE(S)\n", g_failures); return 1; }
    printf("test_parallax_zeropoint: all passed\n");
    return 0;
}
