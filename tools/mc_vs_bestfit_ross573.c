// Diagnostic, not production code: directly tests the one remaining untested hypothesis for the
// "open discrepancy" documented at length in NOTES_method_selection.md and
// tests/test_pipeline_borisov_ross573.c / tests/test_pipeline_oumuamua_hip3757.c -- that
// Bailer-Jones et al.'s reported numbers are MEDIANS of ~2000 Monte Carlo surrogates (resampling
// BOTH the star's and the ISO's full phase-space uncertainty), while this project's pipeline
// tests so far have only ever reported the SINGLE deterministic best-fit trajectory. Those two
// quantities are not the same thing in general, and for exactly the geometric reason explored
// here: a deterministic trajectory has one particular (not necessarily minimal) miss distance,
// while resampling over the full uncertainty ellipsoid of BOTH bodies explores many nearby
// trajectories, some of which pass much closer purely by chance -- so the MEDIAN over many
// resampled trajectories can plausibly sit well below the single deterministic value, especially
// when position uncertainty (not just velocity) is included, which mc_search_stage1() deliberately
// does NOT do for the ISO (see mc_search.h's own comment on why ISO position isn't resampled --
// this tool, being a one-off diagnostic rather production BOINC-quorum code, isn't bound by that
// same cost constraint).
//
// Method: resample BOTH bodies' FULL 6D state (not just velocity) via real uncertainty data:
//   - Ross 573 (the star): Gaia DR2's own 5x5 correlated astrometric covariance (ra, dec,
//     parallax, pmra, pmdec) plus its independent radial-velocity error -- queried directly from
//     gaiadr2.gaia_source (source_id 5162123155863791744), same source as
//     tests/test_pipeline_borisov_ross573.c's mean values.
//   - 2I/Borisov (the ISO): JPL Small-Body Database's own 6x6 marginal covariance over
//     (e, q, tp, node, peri, i) -- same source/method as
//     tools/oumuamua_velocity_covariance.c's 'Oumuamua covariance, just for Borisov's solution 54
//     (fetched 2026-10-06 via https://ssd-api.jpl.nasa.gov/sbdb.api?sstr=2I&full-prec=true&cov=mat).
// For each of 2000 draws (matching Bailer-Jones et al.'s own count), sample a full star state and
// a full ISO state independently, run the SAME find_closest_approach() this project's pipeline
// tests use, and record the resulting minimum distance/time/speed. Report the median and a 90%
// interval across the 2000 draws, for direct comparison against both the deterministic best-fit
// (tests/test_pipeline_borisov_ross573.c: d=1.06 pc, t=-0.893 Myr, v=23.0 km/s) and
// Bailer-Jones et al. 2020's reported MC median (d=0.068 pc, t=-909 kyr, v=22.6 km/s).
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/coords.h"
#include "../src/orbital_elements.h"
#include "../src/encounter.h"
#include "../src/potential.h"
#include "../src/units.h"
#include "../src/montecarlo.h"
#include "../src/star_mc.h"

#define N_DRAWS 2000
#define SEARCH_WINDOW_MYR 5.0
#define N_STEPS (int)(SEARCH_WINDOW_MYR * 100000)

// Ross 573, Gaia DR2 source_id 5162123155863791744 -- means (same as
// tests/test_pipeline_borisov_ross573.c) plus the full 5x5 correlated astrometric error matrix
// and independent RV error, queried 2026-10-06.
static const double STAR_MEAN[MC_MAX_DIM] = {
    52.33270555516616, -11.679678598930318, 47.536273832249684,
    54.42974631770208, -304.38948045019396, 13.060251980364459,
}; // ra_deg, dec_deg, parallax_mas, pmra_masyr, pmdec_masyr, rv_kms

// Raw Gaia-archive-queried per-parameter standard errors (ra, dec, parallax, pmra, pmdec, rv;
// ra/dec errors in MAS of alpha*/delta as Gaia reports them, the rest in STAR_MEAN's units) and the 10 pairwise correlation coefficients among the first five
// (rv is Gaia's own independent spectroscopic measurement, uncorrelated with the astrometric
// solution) -- the covariance is built from these by src/star_mc.c's gaia_build_cov6() as
// cov[i][j] = corr[i][j] * sigma[i] * sigma[j] (diagonal: sigma[i]^2), rather than hand-computing
// and transcribing 15 products, which is exactly the kind of arithmetic that's easy to get
// subtly wrong without it showing up as an obvious failure.
static const double STAR_SIGMA[MC_MAX_DIM] = {
    0.02758606936747954, 0.023243752052676182, 0.034845658889076744,
    0.059758726384874165, 0.06374894646745397, 0.401853480238817,
};
// order: ra_dec, ra_parallax, ra_pmra, ra_pmdec, dec_parallax, dec_pmra, dec_pmdec,
// parallax_pmra, parallax_pmdec, pmra_pmdec
static const double STAR_CORR[10] = {
    -0.079618886, 0.29959196, -0.06498557, 0.09616806,
    -0.118148394, 0.17030483, 0.26450762,
    0.048970062, -0.107570626,
    0.15593974,
};
// Bug found 2026-10-07 (NOTES_method_selection.md, "Stage-2 Monte Carlo of the 99-star
// shortlist"): this tool originally built a 6x6 covariance in (ra_deg, dec_deg, ...) space from
// the sigmas above and sampled ra/dec directly -- but Gaia's ra_error/dec_error are in MAS (and
// ra_error is for alpha* = ra*cos(dec)), so every RA/Dec draw was perturbed by ~0.03 DEGREES
// instead of ~0.03 mas (~100 arcsec, ~0.01 pc at Ross 573's distance). Now delegated to
// src/star_mc.c, which keeps the covariance in Gaia's native units and converts the drawn
// offsets to degrees correctly (tested in tests/test_star_mc.c).
static GaiaStarUncertainty star_unc(void) {
    GaiaStarUncertainty s = {
        .source_id = 5162123155863791744ULL,
        .mean = { STAR_MEAN[0], STAR_MEAN[1], STAR_MEAN[2], STAR_MEAN[3], STAR_MEAN[4], STAR_MEAN[5] },
        .ra_error_mas = STAR_SIGMA[0], .dec_error_mas = STAR_SIGMA[1],
        .parallax_error_mas = STAR_SIGMA[2], .pmra_error_masyr = STAR_SIGMA[3],
        .pmdec_error_masyr = STAR_SIGMA[4], .rv_error_kms = STAR_SIGMA[5],
    };
    for (int k = 0; k < GAIA_N_CORR; k++) s.corr[k] = STAR_CORR[k];
    return s;
}

// 2I/Borisov, JPL SBDB Solution 54 (same elements as tests/test_pipeline_borisov_ross573.c) plus
// its 6x6 marginal covariance over (e, q, tp, node, peri, i) -- top-left block of the real 10x10
// SBDB covariance (the other four parameters, A1/A2/A3/DT, are non-gravitational-model
// parameters this project's two-body Keplerian model doesn't use; marginalizing them out by
// taking the submatrix is exact, same reasoning as tools/oumuamua_velocity_covariance.c).
static const double ISO_MEAN[MC_MAX_DIM] = {
    3.356475782676596, 2.006520878500843, 0.0, 308.1477292269942, 209.1236864378081, 44.05264247909138,
}; // e, q_au, tp_offset_day (0 = best fit), node_deg, peri_deg, i_deg
static const double ISO_COV[MC_MAX_DIM][MC_MAX_DIM] = {
    { 3.624937870909242e-10, 4.794054472467175e-11, -2.464202390944168e-10, 5.196830156130616e-10, 7.456882531816074e-11, -1.620775619495585e-10 },
    { 4.794054472467175e-11, 1.082666613279811e-11, 2.151068454174362e-10, 1.360371010275828e-10, 2.105858394578639e-10, -2.130754641675636e-11 },
    { -2.464202390944168e-10, 2.151068454174362e-10, 2.018072950651087e-8, 3.983156791931724e-9, 1.534464172414489e-8, 5.375873483794185e-10 },
    { 5.196830156130616e-10, 1.360371010275828e-10, 3.983156791931724e-9, 2.026904186820056e-9, 3.416233822860535e-9, -1.456300667469856e-10 },
    { 7.456882531816074e-11, 2.105858394578639e-10, 1.534464172414489e-8, 3.416233822860535e-9, 1.21260867833924e-8, 2.316543101411908e-10 },
    { -1.620775619495585e-10, -2.130754641675636e-11, 5.375873483794185e-10, -1.456300667469856e-10, 2.316543101411908e-10, 2.779994339873895e-10 },
}; // order: e, q, tp, node, peri, i -- matches SBDB's own "labels"

static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    return (da > db) - (da < db);
}

static double percentile(double *sorted, int n, double p) {
    double idx = p * (n - 1);
    int lo = (int)floor(idx), hi = (int)ceil(idx);
    if (lo == hi) return sorted[lo];
    double frac = idx - lo;
    return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
}

int main(void) {
    GaiaStarUncertainty star = star_unc();
    double star_cov[MC_MAX_DIM][MC_MAX_DIM];
    gaia_build_cov6(&star, star_cov);

    double star_L[MC_MAX_DIM][MC_MAX_DIM], iso_L[MC_MAX_DIM][MC_MAX_DIM];
    if (mc_cholesky(star_cov, 6, star_L) != 0) { fprintf(stderr, "star covariance not PD\n"); return 1; }
    if (mc_cholesky(ISO_COV, 6, iso_L) != 0) { fprintf(stderr, "ISO covariance not PD\n"); return 1; }

    McRng star_rng, iso_rng;
    mc_rng_seed(&star_rng, 1001);
    mc_rng_seed(&iso_rng, 2002);

    Potential mw = milkyway_potential(&milkyway_potential_default);
    double t_now_myr = -3.0 * 365.25 * NOSTOS_MYR_PER_DAY; // the INBOUND asymptote, same
                                                             // evaluation point as
                                                             // tests/test_pipeline_borisov_ross573.c
                                                             // (see its top comment for why
                                                             // inbound, not outbound, matters)

    double *dist = malloc(N_DRAWS * sizeof(double));
    double *time = malloc(N_DRAWS * sizeof(double));
    double *speed = malloc(N_DRAWS * sizeof(double));
    int n_valid = 0;

    for (int k = 0; k < N_DRAWS; k++) {
        double iso_draw[MC_MAX_DIM];
        Astrometry star_a;
        int star_ok = gaia_sample_astrometry(&star_rng, &star, star_L, &star_a);
        mc_sample_multivariate_gaussian(&iso_rng, ISO_MEAN, iso_L, 6, iso_draw);

        if (!star_ok || iso_draw[0] <= 1.0 || iso_draw[1] <= 0.0) continue; // reject a
            // nonphysical draw (negative parallax or non-hyperbolic/negative-q elements) rather
            // than silently feeding it through -- rare at these tight uncertainties, but a real
            // possibility for a 3-sigma-tail draw, and not handled by either underlying function.

        PhaseSpace star_ps;
        astrometry_to_galactocentric(star_a, &solar_params_default, &star_ps);

        OrbitalElements iso_el = {
            .e = iso_draw[0], .q_au = iso_draw[1], .i_deg = iso_draw[5],
            .node_deg = iso_draw[3], .peri_deg = iso_draw[4],
        };
        double t_myr = t_now_myr - iso_draw[2] * NOSTOS_MYR_PER_DAY;
        PhaseSpace iso_ps;
        if (!orbital_elements_to_galactocentric(iso_el, t_myr, &solar_params_default, &iso_ps)) continue;

        State star_state = {star_ps.pos, star_ps.vel};
        State iso_state = {iso_ps.pos, iso_ps.vel};
        EncounterResult r = find_closest_approach(&mw, &star_state, &iso_state, -SEARCH_WINDOW_MYR, N_STEPS);

        dist[n_valid] = r.min_distance_kpc * 1000.0;
        time[n_valid] = r.time_at_min_myr;
        speed[n_valid] = nostos_kpcmyr_to_kms(v3_norm(v3_sub(r.star_vel_at_min, r.iso_vel_at_min)));
        n_valid++;
    }

    printf("%d of %d draws valid\n", n_valid, N_DRAWS);
    qsort(dist, n_valid, sizeof(double), cmp_double);
    qsort(time, n_valid, sizeof(double), cmp_double);
    qsort(speed, n_valid, sizeof(double), cmp_double);

    printf("\ndistance (pc): median=%.4f  90%% CI=[%.4f, %.4f]\n",
           percentile(dist, n_valid, 0.5), percentile(dist, n_valid, 0.05), percentile(dist, n_valid, 0.95));
    printf("time (Myr):    median=%.4f  90%% CI=[%.4f, %.4f]\n",
           percentile(time, n_valid, 0.5), percentile(time, n_valid, 0.05), percentile(time, n_valid, 0.95));
    printf("speed (km/s):  median=%.4f  90%% CI=[%.4f, %.4f]\n",
           percentile(speed, n_valid, 0.5), percentile(speed, n_valid, 0.05), percentile(speed, n_valid, 0.95));

    printf("\nfor comparison:\n");
    printf("  this project's deterministic best fit: d=1.06 pc, t=-0.893 Myr, v=23.0 km/s\n");
    printf("  Bailer-Jones et al. 2020 (MC median):  d=0.068 pc, t=-0.909 Myr, v=22.6 km/s\n");

    free(dist); free(time); free(speed);
    return 0;
}
