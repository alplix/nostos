// Star-side Monte Carlo resampling from a Gaia source's own reported uncertainties -- the
// "stage 2" half of NOTES_method_selection.md's two-stage Monte Carlo design (Bailer-Jones et
// al. 2019: clones of each candidate STAR, Gaussian in the observed astrometric quantities, then
// propagated to phase space the normal way). Stage 1 (src/mc_search.c) only resamples the ISO.
//
// Units, stated explicitly because this is exactly where a real bug was found (2026-10-07, in
// tools/mc_vs_bestfit_ross573.c, which built its covariance inline before this module existed):
// Gaia's ra_error/dec_error are in MILLIARCSECONDS (and ra_error is sigma of ra*cos(dec), i.e.
// alpha*, not of ra itself), while ra/dec are in DEGREES. That tool added a mas-sized sigma
// directly to a degree-valued mean -- a ~3.6e6x-too-large position perturbation (0.028 "deg"
// instead of 0.028 mas). This module keeps the covariance in Gaia's own native units
// (mas for alpha*/delta offsets, mas, mas/yr, km/s) and converts the drawn alpha*/delta OFFSETS
// to degrees only when applying them (gaia_apply_offsets()), which is the one place the
// cos(dec) and mas->deg factors appear.
#ifndef NOSTOS_STAR_MC_H
#define NOSTOS_STAR_MC_H

#include <stdint.h>
#include <stdbool.h>
#include "coords.h"
#include "montecarlo.h"

// Index order of the 10 Gaia correlation coefficients, exactly the Gaia archive's own column
// order: ra_dec, ra_parallax, ra_pmra, ra_pmdec, dec_parallax, dec_pmra, dec_pmdec,
// parallax_pmra, parallax_pmdec, pmra_pmdec.
#define GAIA_N_CORR 10

typedef struct {
    uint64_t source_id;
    Astrometry mean;          // ra/dec degrees, parallax mas, pm mas/yr, rv km/s
    double ra_error_mas;      // sigma of alpha* = ra*cos(dec), mas (Gaia's own convention)
    double dec_error_mas;
    double parallax_error_mas;
    double pmra_error_masyr;
    double pmdec_error_masyr;
    double rv_error_kms;
    double corr[GAIA_N_CORR]; // see GAIA_N_CORR's comment for the order; 0 for a missing value
} GaiaStarUncertainty;

// Builds the 6x6 covariance in native Gaia units, ordered (alpha* offset mas, delta offset mas,
// parallax mas, pmra mas/yr, pmdec mas/yr, rv km/s): cov[i][j] = corr_ij * sigma_i * sigma_j for
// the 5 astrometric parameters, and rv uncorrelated with all of them (Gaia's RV is an independent
// spectroscopic measurement; DR3 publishes no astrometry-RV correlation).
void gaia_build_cov6(const GaiaStarUncertainty *s, double cov[MC_MAX_DIM][MC_MAX_DIM]);

// Applies one drawn offset vector (same order/units as gaia_build_cov6()) to the star's mean
// astrometry: ra += d_alpha*/cos(dec)/3.6e6 deg, dec += d_delta/3.6e6 deg, the rest added
// directly.
Astrometry gaia_apply_offsets(const Astrometry *mean, const double offset[MC_MAX_DIM]);

// Draws one resampled astrometry. L must be mc_cholesky() of gaia_build_cov6()'s output.
// Returns false for a non-physical draw (parallax <= 0) -- callers should skip it and count it,
// not silently feed it through astrometry_to_galactocentric().
bool gaia_sample_astrometry(McRng *rng, const GaiaStarUncertainty *s,
                            const double L[MC_MAX_DIM][MC_MAX_DIM], Astrometry *out);

#endif
