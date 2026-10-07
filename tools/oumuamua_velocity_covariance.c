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
// evaluates at t_myr corresponding to ~9.08 years BEFORE perihelion -- the INBOUND asymptote,
// the velocity direction 'Oumuamua actually had while approaching through the Galaxy, which is
// what the long Milky-Way-potential-only backward integration needs (see that test file's top
// comment for the real bug this matters for: the outbound, post-perihelion asymptote points in
// a materially different direction, by the orbit's full hyperbolic turning angle, and silently
// carrying that wrong direction through a backward integration that never re-simulates the
// Sun's actual deflection was the dominant source of this project's long-open encounter-distance
// discrepancy). The magnitude (9.08 years) is where the velocity reaches the asymptotic
// hyperbolic-excess regime (stops changing fast) -- same convergence rate on the inbound branch
// as the outbound one, by the two-body Kepler problem's time symmetry. The same t_myr is used
// here, since the velocity covariance is only meaningful once the velocity itself has stopped
// evolving quickly with time.
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
#include "../src/iso_solutions.h"

// Best-fit elements at the covariance's own epoch (JD 2458080.5), read directly off the SBDB
// response fetched 2026-10-06 -- NOT this project's other hardcoded OUMUAMUA_ELEMENTS (see the
// top-of-file caveat on why these differ slightly, and why that's fine here).
#define E_BEST     1.201133796102373
#define Q_AU_BEST  0.2559115812959116
#define I_DEG_BEST 122.7417062847286
#define NODE_DEG_BEST 24.59690955523242
#define PERI_DEG_BEST 241.8105360304898

// Evaluation time: ~9.08 years BEFORE perihelion -- the inbound asymptote, this project's own
// settled correct starting point (tests/test_pipeline_oumuamua_hip3757.c).
// 2026-10-08: now the EXACT inbound asymptotic velocity (orbital_elements_inbound_asymptote()), not
// the Keplerian velocity at -9.08 yr, which was still 2.3% too fast (orbital_elements.h). On the
// asymptote the velocity no longer depends on tp at all (tp only slides the position along the
// line), so the tp row of the Jacobian is now exactly zero.
#define T_MYR_EVAL (OUMUAMUA_T_START_DR3_DAYS * NOSTOS_MYR_PER_DAY)

// 6x6 covariance submatrix (e, q[au], tp[day], node[deg], peri[deg], i[deg]) -- the real SBDB
// orbit-fit covariance block, now kept in one shared place (src/iso_solutions.h) so this tool and
// tools/mc_shortlist.c cannot drift apart. Row/column order: e, q, tp, node, peri, i.
#define COV6 OUMUAMUA_SBDB_COV6

static void eval_vel(double e, double q_au, double tp_day_offset, double node_deg,
                      double peri_deg, double i_deg, const SolarParams *sun, Vec3 *vel_out) {
    OrbitalElements el = { e, q_au, i_deg, node_deg, peri_deg };
    double t_myr = T_MYR_EVAL - tp_day_offset * NOSTOS_MYR_PER_DAY;
    PhaseSpace ps;
    if (!orbital_elements_inbound_asymptote(el, t_myr, sun, &ps)) {
        fprintf(stderr, "orbital_elements_inbound_asymptote failed for a perturbed point\n");
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

    printf("SBDB osculating-element covariance propagated to the asymptotic velocity, (kpc/Myr)^2 "
           "(cov00 cov01 cov02 cov11 cov12 cov22):\n");
    printf("%.10e %.10e %.10e %.10e %.10e %.10e\n",
           covvel[0][0], covvel[0][1], covvel[0][2], covvel[1][1], covvel[1][2], covvel[2][2]);

    // 2026-10-08: the covariance scripts/make_work_units.py should actually use. The SBDB block above
    // leaves out the non-gravitational parameters (A1-A3), whose uncertainty propagates into the
    // inbound asymptote through the whole passage; Bailer-Jones et al. 2018's Table 1 publishes
    // the full asymptotic-velocity covariance for this solution (7c), ~40x larger. Rotate it from
    // barycentric ICRF km/s into this project's Galactic axes: C_gal = A C_icrs A^T.
    static const SolarParams no_sun = { 0, 0, 0, 0, 0, 0 };
    double c_icrs[3][3];
    {
        const double *sg = OUMUAMUA_7C_VINF_SIGMA_KMS, *rc = OUMUAMUA_7C_VINF_CORR;
        const double rho[3][3] = { { 1.0, rc[0], rc[1] }, { rc[0], 1.0, rc[2] }, { rc[1], rc[2], 1.0 } };
        for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) c_icrs[a][b] = rho[a][b] * sg[a] * sg[b];
    }
    Vec3 col[3]; // A applied to the ICRS basis vectors = the columns of A
    col[0] = helio_equatorial_to_galactocentric_vel((Vec3){1, 0, 0}, &no_sun);
    col[1] = helio_equatorial_to_galactocentric_vel((Vec3){0, 1, 0}, &no_sun);
    col[2] = helio_equatorial_to_galactocentric_vel((Vec3){0, 0, 1}, &no_sun);
    double A[3][3];
    for (int j = 0; j < 3; j++) { A[0][j] = col[j].x; A[1][j] = col[j].y; A[2][j] = col[j].z; }
    double c_gal[3][3];
    for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) {
        double s = 0.0;
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) s += A[a][i] * c_icrs[i][j] * A[b][j];
        c_gal[a][b] = s / (NOSTOS_KMS_PER_KPCMYR * NOSTOS_KMS_PER_KPCMYR);
    }
    fprintf(stderr, "\nBailer-Jones et al. 2018 Table 1 (7c) covariance, Galactic axes, per-axis 1-sigma: "
            "%.6f %.6f %.6f km/s\n", sqrt(c_gal[0][0]) * NOSTOS_KMS_PER_KPCMYR,
            sqrt(c_gal[1][1]) * NOSTOS_KMS_PER_KPCMYR, sqrt(c_gal[2][2]) * NOSTOS_KMS_PER_KPCMYR);
    printf("OUMUAMUA_VEL_COV_KPCMYR2 for scripts/make_work_units.py (Bailer-Jones 2018 Table 1, 7c; "
           "cov00 cov01 cov02 cov11 cov12 cov22):\n");
    printf("%.10e %.10e %.10e %.10e %.10e %.10e\n",
           c_gal[0][0], c_gal[0][1], c_gal[0][2], c_gal[1][1], c_gal[1][2], c_gal[2][2]);

    return 0;
}
