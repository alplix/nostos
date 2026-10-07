// Validates src/star_mc.c (star-side Monte Carlo resampling from Gaia's own uncertainties).
// Three kinds of checks:
//   1. Exact: gaia_build_cov6() puts sigma^2 on the diagonal, corr*sigma_i*sigma_j in the right
//      (symmetric) off-diagonal slots in Gaia's own column order, and leaves rv uncorrelated.
//   2. Exact: gaia_apply_offsets() converts an alpha*/delta offset in MAS to a degree shift with
//      the cos(dec) factor on RA -- the unit trap that tools/mc_vs_bestfit_ross573.c fell into
//      (mas sigma added straight to a degree mean), now pinned down by a test.
//   3. Statistical: many draws reproduce the requested per-parameter sigmas (in Gaia's native
//      units, measured back from the drawn degrees) and one requested correlation, within a
//      generous multiple of the expected sampling error.
#include <stdio.h>
#include <math.h>
#include "../src/star_mc.h"

static int g_failures = 0;

static void check(int ok, const char *what) {
    if (!ok) { printf("  FAIL: %s\n", what); g_failures++; }
    else { printf("  ok:   %s\n", what); }
}

static GaiaStarUncertainty example_star(void) {
    GaiaStarUncertainty s = {
        .source_id = 1,
        .mean = { .ra_deg = 120.0, .dec_deg = 60.0, .parallax_mas = 10.0,
                  .pmra_masyr = 50.0, .pmdec_masyr = -30.0, .rv_kms = 12.0 },
        .ra_error_mas = 0.02, .dec_error_mas = 0.03, .parallax_error_mas = 0.04,
        .pmra_error_masyr = 0.05, .pmdec_error_masyr = 0.06, .rv_error_kms = 1.5,
        .corr = { 0.1, -0.1, 0.2, -0.1, 0.1, -0.15, 0.2, 0.5, -0.1, 0.15 },
    };
    return s;
}

static void test_cov_layout(void) {
    printf("test_cov_layout:\n");
    GaiaStarUncertainty s = example_star();
    double cov[MC_MAX_DIM][MC_MAX_DIM];
    gaia_build_cov6(&s, cov);
    double sig[6] = {0.02, 0.03, 0.04, 0.05, 0.06, 1.5};
    int diag_ok = 1;
    for (int i = 0; i < 6; i++) if (fabs(cov[i][i] - sig[i] * sig[i]) > 1e-15) diag_ok = 0;
    check(diag_ok, "diagonal = sigma^2 for all six parameters");
    check(fabs(cov[0][1] - 0.1 * 0.02 * 0.03) < 1e-15 && cov[0][1] == cov[1][0], "ra_dec_corr -> [0][1], symmetric");
    check(fabs(cov[2][3] - 0.5 * 0.04 * 0.05) < 1e-15, "parallax_pmra_corr -> [2][3]");
    check(fabs(cov[3][4] - 0.15 * 0.05 * 0.06) < 1e-15, "pmra_pmdec_corr -> [3][4]");
    check(fabs(cov[1][4] - 0.2 * 0.03 * 0.06) < 1e-15, "dec_pmdec_corr -> [1][4]");
    int rv_ok = 1;
    for (int i = 0; i < 5; i++) if (cov[i][5] != 0.0 || cov[5][i] != 0.0) rv_ok = 0;
    check(rv_ok, "rv uncorrelated with the astrometric parameters");
    double L[MC_MAX_DIM][MC_MAX_DIM];
    check(mc_cholesky(cov, 6, L) == 0, "a realistic Gaia covariance is positive-definite");
}

static void test_offset_units(void) {
    printf("test_offset_units:\n");
    Astrometry m = { .ra_deg = 120.0, .dec_deg = 60.0, .parallax_mas = 10.0,
                     .pmra_masyr = 50.0, .pmdec_masyr = -30.0, .rv_kms = 12.0 };
    double off[MC_MAX_DIM] = { 3.6e6, 1.8e6, 0.5, 1.0, -2.0, 3.0 };
    Astrometry a = gaia_apply_offsets(&m, off);
    // 3.6e6 mas of alpha* = 1 deg of great-circle arc = 1/cos(60deg) = 2 deg of RA.
    check(fabs(a.ra_deg - 122.0) < 1e-9, "1 deg of alpha* at dec=60 -> +2 deg of RA");
    check(fabs(a.dec_deg - 60.5) < 1e-12, "1.8e6 mas of delta -> +0.5 deg of Dec");
    check(fabs(a.parallax_mas - 10.5) < 1e-12 && fabs(a.pmra_masyr - 51.0) < 1e-12 &&
          fabs(a.pmdec_masyr + 32.0) < 1e-12 && fabs(a.rv_kms - 15.0) < 1e-12,
          "parallax/pm/rv offsets added directly in native units");
    double small[MC_MAX_DIM] = { 0.03, 0, 0, 0, 0, 0 };
    Astrometry b = gaia_apply_offsets(&m, small);
    check(fabs(b.ra_deg - m.ra_deg) < 1e-7, "a 0.03 mas offset moves RA by ~1e-8 deg, NOT 0.03 deg "
          "(the tools/mc_vs_bestfit_ross573.c unit trap)");
}

static void test_sampling_statistics(void) {
    printf("test_sampling_statistics:\n");
    GaiaStarUncertainty s = example_star();
    double cov[MC_MAX_DIM][MC_MAX_DIM], L[MC_MAX_DIM][MC_MAX_DIM];
    gaia_build_cov6(&s, cov);
    mc_cholesky(cov, 6, L);
    McRng rng;
    mc_rng_seed(&rng, 777ULL);
    const int N = 400000;
    double sum[6] = {0}, sum2[6] = {0}, sum_plx_pmra = 0;
    double cosd = cos(60.0 * 3.14159265358979323846 / 180.0);
    int n_ok = 0;
    for (int k = 0; k < N; k++) {
        Astrometry a;
        if (!gaia_sample_astrometry(&rng, &s, L, &a)) continue;
        double x[6] = {
            (a.ra_deg - s.mean.ra_deg) * cosd * 3.6e6, (a.dec_deg - s.mean.dec_deg) * 3.6e6,
            a.parallax_mas - s.mean.parallax_mas, a.pmra_masyr - s.mean.pmra_masyr,
            a.pmdec_masyr - s.mean.pmdec_masyr, a.rv_kms - s.mean.rv_kms,
        };
        for (int i = 0; i < 6; i++) { sum[i] += x[i]; sum2[i] += x[i] * x[i]; }
        sum_plx_pmra += x[2] * x[3];
        n_ok++;
    }
    check(n_ok == N, "no draws rejected for a parallax 250 sigma from zero");
    double sig[6] = {0.02, 0.03, 0.04, 0.05, 0.06, 1.5};
    int sig_ok = 1;
    for (int i = 0; i < 6; i++) {
        double sd = sqrt(sum2[i] / n_ok - (sum[i] / n_ok) * (sum[i] / n_ok));
        printf("  param %d: empirical sigma %.5g vs requested %.5g\n", i, sd, sig[i]);
        // relative sampling error of a std estimate ~ 1/sqrt(2N) ~ 1.1e-3; allow 1%.
        if (fabs(sd / sig[i] - 1.0) > 0.01) sig_ok = 0;
    }
    check(sig_ok, "empirical sigmas match requested sigmas (in native Gaia units) to 1%");
    double corr = (sum_plx_pmra / n_ok) / (0.04 * 0.05);
    printf("  parallax-pmra empirical corr %.4f vs requested 0.5\n", corr);
    check(fabs(corr - 0.5) < 0.01, "requested parallax_pmra correlation reproduced");
}

static void test_rejects_nonphysical_parallax(void) {
    printf("test_rejects_nonphysical_parallax:\n");
    GaiaStarUncertainty s = example_star();
    s.mean.parallax_mas = 0.01;
    s.parallax_error_mas = 1.0; // ~half of all draws go negative
    double cov[MC_MAX_DIM][MC_MAX_DIM], L[MC_MAX_DIM][MC_MAX_DIM];
    gaia_build_cov6(&s, cov);
    mc_cholesky(cov, 6, L);
    McRng rng;
    mc_rng_seed(&rng, 5ULL);
    int n_rej = 0, bad_out = 0;
    for (int k = 0; k < 10000; k++) {
        Astrometry a;
        if (!gaia_sample_astrometry(&rng, &s, L, &a)) n_rej++;
        else if (a.parallax_mas <= 0) bad_out++;
    }
    printf("  %d of 10000 rejected\n", n_rej);
    check(n_rej > 4000 && n_rej < 6000, "roughly half of draws rejected when parallax ~ 0 +- 1");
    check(bad_out == 0, "no accepted draw ever has parallax <= 0");
}

int main(void) {
    test_cov_layout();
    test_offset_units();
    test_sampling_statistics();
    test_rejects_nonphysical_parallax();
    if (g_failures) { printf("\n%d FAILURE(S)\n", g_failures); return 1; }
    printf("\nall star_mc tests passed\n");
    return 0;
}
