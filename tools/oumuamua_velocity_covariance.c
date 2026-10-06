// One-off computation, not a regression test: propagates 'Oumuamua's REAL orbit-fit covariance
// matrix (JPL Small-Body Database, Solution 7c / Micheli et al. 2018, fetched 2026-10-06 via
// https://ssd-api.jpl.nasa.gov/sbdb.api?sstr=1I&full-prec=true&cov=mat) through
// orbital_elements_to_galactocentric() to produce a real Galactocentric velocity covariance in
// (kpc/Myr)^2 -- the quantity scripts/make_work_units.py needs for mc_search_stage1()'s Monte
// Carlo clones, and which that script had been using a placeholder (1e-6 (kpc/Myr)^2 on the
// diagonal, explicitly flagged as such) for instead.
//
// Method: standard linear uncertainty propagation. The orbit solution's 9x9 covariance is over
// (e, q, tp, node, peri, i, A1, A2, A3) -- the first six are the elements
// orbital_elements_to_galactocentric() actually takes (tp, time of perihelion passage, enters
// through this project's t_myr argument rather than as a struct field -- see below); A1/A2/A3
// are non-gravitational-acceleration model parameters this project's two-body Keplerian model
// doesn't use at all. Marginalizing them out is NOT an approximation here: the marginal
// covariance of a subset of jointly-Gaussian variables is exactly the corresponding principal
// submatrix, so taking the top-left 6x6 block of the 9x9 is exact, not a simplification.
//
// tp's role: OrbitalElements has no tp field -- orbital_elements_to_galactocentric() takes
// t_myr directly, "time relative to perihelion passage". Evaluating at a fixed absolute moment
// while tp varies is exactly equivalent to perturbing t_myr itself by -(delta_tp in Myr), so tp's
// uncertainty is propagated as a perturbation of the t_myr argument, not a struct field.
//
// Evaluation point: this project's own pipeline (tests/test_pipeline_oumuamua_hip3757.c)
// evaluates at t_myr corresponding to ~9.08 years post-perihelion, the point found to reach the
// asymptotic hyperbolic-excess-velocity regime (velocity has stopped changing fast -- see that
// test's own comment). The same t_myr is used here, since the velocity covariance is only
// meaningful once the velocity itself has stopped evolving quickly with time.
//
// Caveat, stated plainly rather than glossed over: the SBDB solution returned by default is
// quoted at epoch JD 2458080.5, while this project's hardcoded OUMUAMUA_ELEMENTS (used
// elsewhere, e.g. scripts/make_work_units.py, tests/test_pipeline_oumuamua_hip3757.c) came from
// a JPL Horizons fetch at epoch JD 2458041.5 -- both are the same underlying orbit determination
// (data_arc 80 days, first_obs 2017-10-14, last_obs 2018-01-02, "Sol. 7c, Micheli et al. 2018"),
// just propagated/reported at slightly different epochs, and the element values agree to
// ~1-2e-4 relative. That's close enough that evaluating the Jacobian at THIS solution's own
// best-fit point (rather than re-deriving a covariance at the other epoch, which SBDB's API
// doesn't offer on demand) is a reasonable, honestly-caveated approximation -- not treated as
// exact.
#include <stdio.h>
#include <math.h>
#include "../src/orbital_elements.h"
#include "../src/coords.h"
#include "../src/units.h"

// Best-fit elements at the covariance's own epoch (JD 2458080.5), read directly off the SBDB
// response fetched 2026-10-06 -- NOT this project's other hardcoded OUMUAMUA_ELEMENTS (see the
// top-of-file caveat on why these differ slightly, and why that's fine here).
#define E_BEST     1.201133796102373
#define Q_AU_BEST  0.2559115812959116
#define I_DEG_BEST 122.7417062847286
#define NODE_DEG_BEST 24.59690955523242
#define PERI_DEG_BEST 241.8105360304898

// Evaluation time: ~9.08 years post-perihelion, this project's own settled asymptotic-regime
// point (tests/test_pipeline_oumuamua_hip3757.c).
#define T_MYR_EVAL (9.08 * 365.25 * NOSTOS_MYR_PER_DAY)

// 6x6 covariance submatrix (e, q[au], tp[day], node[deg], peri[deg], i[deg]) -- exactly the
// top-left 6x6 block of the SBDB response's 9x9 "covariance.data", in its own native element
// units (JPL reports covariances directly in each element's native unit, not normalized).
static const double COV6[6][6] = {
    { 4.437021723735803e-10, 1.402645684430185e-10, -2.290597192603328e-9, -4.866487927086338e-9,  2.630680330585727e-8,  3.783694880329765e-9 },
    { 1.402645684430185e-10, 4.44023866269166e-11,  -7.704816361638409e-10,-1.55813687037147e-9,   8.310783593401902e-9,  1.153616408668474e-9 },
    {-2.290597192603328e-9, -7.704816361638409e-10,  6.982193721062611e-8,  5.061551619658835e-8, -1.267295484829105e-7,  3.435990038007021e-8 },
    {-4.866487927086338e-9, -1.55813687037147e-9,    5.061551619658835e-8,  6.462880452547438e-8, -2.844751574251313e-7, -1.765593440441589e-8 },
    { 2.630680330585727e-8,  8.310783593401902e-9,  -1.267295484829105e-7, -2.844751574251313e-7,  1.561285075915765e-6,  2.328863207568894e-7 },
    { 3.783694880329765e-9,  1.153616408668474e-9,   3.435990038007021e-8, -1.765593440441589e-8,  2.328863207568894e-7,  8.309282984005662e-8 },
};
// Row/column order: e, q, tp, node, peri, i -- matches SBDB's own "labels" field.

static void eval_vel(double e, double q_au, double tp_day_offset, double node_deg,
                      double peri_deg, double i_deg, const SolarParams *sun, Vec3 *vel_out) {
    OrbitalElements el = { e, q_au, i_deg, node_deg, peri_deg };
    double t_myr = T_MYR_EVAL - tp_day_offset * NOSTOS_MYR_PER_DAY;
    PhaseSpace ps;
    if (!orbital_elements_to_galactocentric(el, t_myr, sun, &ps)) {
        fprintf(stderr, "orbital_elements_to_galactocentric failed for a perturbed point\n");
        *vel_out = (Vec3){0, 0, 0};
        return;
    }
    *vel_out = ps.vel;
}

int main(void) {
    const SolarParams *sun = &solar_params_default;

    // Step sizes for central finite differences -- each a small fraction of that element's own
    // sigma (sqrt of its diagonal COV6 entry), so the perturbation is well within the locally-
    // linear regime the Jacobian approximation needs, but large enough to stay well above
    // double-precision floor given orbital_elements_to_galactocentric()'s own trig/hyperbolic
    // chain.
    double sigma[6];
    for (int i = 0; i < 6; i++) sigma[i] = sqrt(COV6[i][i]);
    double h[6];
    for (int i = 0; i < 6; i++) h[i] = sigma[i] * 0.1;

    double best[6] = { E_BEST, Q_AU_BEST, 0.0, NODE_DEG_BEST, PERI_DEG_BEST, I_DEG_BEST };
    // (index 2, tp, is an OFFSET in days from the best-fit tp, so its best value is 0)

    double J[6][3]; // J[i][k] = d(vel_k)/d(x_i)

    for (int i = 0; i < 6; i++) {
        double xp[6], xm[6];
        for (int j = 0; j < 6; j++) { xp[j] = best[j]; xm[j] = best[j]; }
        xp[i] += h[i];
        xm[i] -= h[i];

        Vec3 vp, vm;
        eval_vel(xp[0], xp[1], xp[2], xp[3], xp[4], xp[5], sun, &vp);
        eval_vel(xm[0], xm[1], xm[2], xm[3], xm[4], xm[5], sun, &vm);

        J[i][0] = (vp.x - vm.x) / (2.0 * h[i]);
        J[i][1] = (vp.y - vm.y) / (2.0 * h[i]);
        J[i][2] = (vp.z - vm.z) / (2.0 * h[i]);
    }

    fprintf(stderr, "Jacobian d(vel_kpcmyr)/d(element) [rows: e,q_au,tp_day,node_deg,peri_deg,i_deg]:\n");
    for (int i = 0; i < 6; i++) {
        fprintf(stderr, "  [%d] %+.6e %+.6e %+.6e\n", i, J[i][0], J[i][1], J[i][2]);
    }

    // Cov_vel (3x3) = J^T * COV6 * J
    double tmp[6][3]; // tmp = COV6 * J
    for (int i = 0; i < 6; i++) {
        for (int k = 0; k < 3; k++) {
            double s = 0.0;
            for (int j = 0; j < 6; j++) s += COV6[i][j] * J[j][k];
            tmp[i][k] = s;
        }
    }
    double covvel[3][3]; // J^T * tmp
    for (int a = 0; a < 3; a++) {
        for (int k = 0; k < 3; k++) {
            double s = 0.0;
            for (int i = 0; i < 6; i++) s += J[i][a] * tmp[i][k];
            covvel[a][k] = s;
        }
    }

    fprintf(stderr, "\nResulting Galactocentric velocity covariance, (kpc/Myr)^2:\n");
    for (int a = 0; a < 3; a++) {
        fprintf(stderr, "  %+.6e %+.6e %+.6e\n", covvel[a][0], covvel[a][1], covvel[a][2]);
    }

    double sigma_kms[3];
    for (int a = 0; a < 3; a++) sigma_kms[a] = sqrt(covvel[a][a]) * NOSTOS_KMS_PER_KPCMYR;
    fprintf(stderr, "\nPer-axis 1-sigma velocity uncertainty: %.6f %.6f %.6f km/s\n",
            sigma_kms[0], sigma_kms[1], sigma_kms[2]);

    printf("PLACEHOLDER_VEL_COV_KPCMYR2 replacement (cov00 cov01 cov02 cov11 cov12 cov22):\n");
    printf("%.10e %.10e %.10e %.10e %.10e %.10e\n",
           covvel[0][0], covvel[0][1], covvel[0][2], covvel[1][1], covvel[1][2], covvel[2][2]);

    return 0;
}
