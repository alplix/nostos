// Validates montecarlo.c's PRNG, Gaussian sampler, and correlated
// multivariate sampling. Two kinds of checks:
//   1. Exact (not statistical): the Cholesky factor L actually satisfies
//      L*L^T = cov, to floating-point precision -- a pure linear-algebra
//      check, no randomness involved.
//   2. Statistical: draw many samples and check the EMPIRICAL mean/
//      covariance converges to the requested mean/covariance within the
//      expected sampling error (~1/sqrt(N) for N draws) -- this is how you
//      validate a random sampler (there's no "exact" answer to check a
//      single draw against), with a generous tolerance so a correct
//      sampler essentially never fails this by bad luck, while a real bug
//      (wrong mean, wrong variance, uncorrelated when correlation was
//      requested) would still show up as an obvious, large discrepancy.
#include <stdio.h>
#include <math.h>
#include "../src/montecarlo.h"

static int g_failures = 0;

static void check(int ok, const char *what) {
    if (!ok) { printf("  FAIL: %s\n", what); g_failures++; }
    else { printf("  ok:   %s\n", what); }
}

static void test_cholesky_exact(void) {
    printf("test_cholesky_exact:\n");
    // A real (not diagonal) 3x3 covariance matrix -- symmetric,
    // positive-definite (diagonally dominant, a sufficient condition).
    double cov[MC_MAX_DIM][MC_MAX_DIM] = {
        {4.0, 1.0, 0.5, 0,0,0},
        {1.0, 9.0, 2.0, 0,0,0},
        {0.5, 2.0, 1.0, 0,0,0},
        {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0},
    };
    double L[MC_MAX_DIM][MC_MAX_DIM];
    int rc = mc_cholesky(cov, 3, L);
    check(rc == 0, "cholesky succeeds on a positive-definite matrix");

    double max_err = 0.0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double sum = 0.0;
            for (int k = 0; k < 3; k++) sum += L[i][k] * L[j][k];
            double err = fabs(sum - cov[i][j]);
            if (err > max_err) max_err = err;
        }
    }
    printf("  max |L*L^T - cov| = %.3e\n", max_err);
    check(max_err < 1e-12, "L*L^T reconstructs the original covariance matrix");

    // A matrix that is NOT positive-definite (negative eigenvalue) should
    // be rejected, not silently produce garbage.
    double bad_cov[MC_MAX_DIM][MC_MAX_DIM] = {
        {1.0, 2.0, 0,0,0,0},
        {2.0, 1.0, 0,0,0,0},
        {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0},
    };
    double L2[MC_MAX_DIM][MC_MAX_DIM];
    int rc2 = mc_cholesky(bad_cov, 2, L2);
    check(rc2 != 0, "cholesky rejects a non-positive-definite matrix");
}

static void test_gaussian_statistics(void) {
    printf("test_gaussian_statistics:\n");
    McRng rng;
    mc_rng_seed(&rng, 12345ULL);
    const int N = 2000000;
    double sum = 0.0, sum2 = 0.0;
    for (int i = 0; i < N; i++) {
        double x = mc_gaussian(&rng);
        sum += x;
        sum2 += x * x;
    }
    double mean = sum / N;
    double var = sum2 / N - mean * mean;
    printf("  N=%d: empirical mean=%.5f (expect ~0), empirical var=%.5f (expect ~1)\n", N, mean, var);
    // Standard error of the mean for N(0,1) over 2e6 draws is ~1/sqrt(2e6)
    // =~ 7e-4; 0.01 is a generous ~14-sigma margin, essentially never fails
    // on a correct sampler but still catches a real bug (e.g. a scale or
    // sign error) immediately.
    check(fabs(mean) < 0.01, "empirical mean of standard normal draws is close to 0");
    check(fabs(var - 1.0) < 0.02, "empirical variance of standard normal draws is close to 1");
}

static void test_multivariate_statistics(void) {
    printf("test_multivariate_statistics:\n");
    double mean[MC_MAX_DIM] = {10.0, -5.0, 2.0, 0,0,0};
    double cov[MC_MAX_DIM][MC_MAX_DIM] = {
        {4.0, -1.5, 0.0, 0,0,0},
        {-1.5, 2.0, 0.5, 0,0,0},
        {0.0, 0.5, 1.0, 0,0,0},
        {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0},
    };
    double L[MC_MAX_DIM][MC_MAX_DIM];
    int rc = mc_cholesky(cov, 3, L);
    check(rc == 0, "cholesky succeeds on the test covariance");

    McRng rng;
    mc_rng_seed(&rng, 987654321ULL);
    const int N = 2000000;
    double sum[3] = {0,0,0};
    double sum2[3][3] = {{0}};
    for (int i = 0; i < N; i++) {
        double out[MC_MAX_DIM];
        mc_sample_multivariate_gaussian(&rng, mean, L, 3, out);
        for (int a = 0; a < 3; a++) sum[a] += out[a];
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                sum2[a][b] += (out[a] - mean[a]) * (out[b] - mean[b]);
            }
        }
    }
    double emp_mean[3];
    for (int a = 0; a < 3; a++) emp_mean[a] = sum[a] / N;
    printf("  empirical mean: (%.4f, %.4f, %.4f) expect (%.1f, %.1f, %.1f)\n",
           emp_mean[0], emp_mean[1], emp_mean[2], mean[0], mean[1], mean[2]);
    check(fabs(emp_mean[0] - mean[0]) < 0.02, "component 0 mean matches");
    check(fabs(emp_mean[1] - mean[1]) < 0.02, "component 1 mean matches");
    check(fabs(emp_mean[2] - mean[2]) < 0.02, "component 2 mean matches");

    printf("  empirical covariance vs requested:\n");
    int cov_ok = 1;
    for (int a = 0; a < 3; a++) {
        for (int b = 0; b < 3; b++) {
            double emp_cov = sum2[a][b] / N;
            printf("    [%d][%d]: empirical=%.4f requested=%.4f\n", a, b, emp_cov, cov[a][b]);
            if (fabs(emp_cov - cov[a][b]) > 0.03) cov_ok = 0;
        }
    }
    check(cov_ok, "empirical covariance matrix matches the requested one (diagonal AND off-diagonal, i.e. the requested correlation is really reproduced)");
}

int main(void) {
    test_cholesky_exact();
    test_gaussian_statistics();
    test_multivariate_statistics();

    if (g_failures) {
        printf("\n%d check(s) FAILED\n", g_failures);
        return 1;
    }
    printf("\nAll checks passed.\n");
    return 0;
}
