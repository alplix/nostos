// Shared published orbit solutions for interstellar objects, so the numbers live in ONE place for
// any tool that needs both a best-fit element set and its covariance (first user:
// tools/mc_shortlist.c; tools/oumuamua_velocity_covariance.c takes its covariance from here too).
// Some tests still carry their own literal copies of the heliocentric best-fit elements.
#ifndef NOSTOS_ISO_SOLUTIONS_H
#define NOSTOS_ISO_SOLUTIONS_H

#include "orbital_elements.h"

// 1I/'Oumuamua heliocentric best-fit elements: JPL Horizons, epoch JD 2458041.5 -- the set the
// catalogue scans used until 2026-10-08 (via the near-Sun Keplerian state; superseded by the
// barycentric 1600 set below for the starting state).
static const OrbitalElements OUMUAMUA_ELEMENTS_HORIZONS = {
    .e = 1.200954455284624, .q_au = 0.2559200355462156,
    .i_deg = 122.7402198359025, .node_deg = 24.59683473134646, .peri_deg = 241.8043518577817,
};

// THE PIPELINE'S STARTING STATE (2026-10-08, see NOTES_method_selection.md "Residual distance gap").
// Every backward Galactic integration starts the ISO on its exact inbound asymptote,
// orbital_elements_inbound_asymptote_gm(<*_ELEMENTS_BARY_1600>, NOSTOS_GM_SSB, t, ...), with t =
// (star-catalogue epoch - barycentric tp), i.e. at the point of the asymptotic line at the moment
// the star's astrometry refers to (the along-line position is a negligible detail, tens of AU).
//
// The element sets are BARYCENTRIC osculating elements (ecliptic J2000) at 1600-01-01 TDB from JPL
// Horizons (centre @0, fetched 2026-10-08), i.e. AFTER JPL's own integration of the whole passage
// backward with planets and the solution's non-gravitational model, at r ~ 2300 AU (1I) /
// 2800 AU (2I), where the remaining solar-system deflection is < 1e-6 rad. Same solutions as
// before: 1I = JPL#16 ("Sol. 7c, Micheli et al. 2018"), 2I = JPL#54. For 1I this reproduces
// Bailer-Jones et al. 2018's Table 1 asymptote for solution 7c to every printed digit (radiant
// RA 279.4987, Dec 33.8649, v_inf 26.391543 km/s; tests/test_orbital_elements.c).
static const OrbitalElements OUMUAMUA_ELEMENTS_BARY_1600 = {
    .e = 1.198093111333325, .q_au = 0.2526436572799559,
    .i_deg = 122.7927648087493, .node_deg = 24.24788585536907, .peri_deg = 241.8944051367842,
};
static const OrbitalElements BORISOV_ELEMENTS_BARY_1600 = {
    .e = 3.362344839910536, .q_au = 2.013421232409241,
    .i_deg = 44.04464609191959, .node_deg = 308.1044985515907, .peri_deg = 209.1125494955898,
};
#define OUMUAMUA_TP_BARY_JD 2458005.600382605568
#define BORISOV_TP_BARY_JD 2458826.211501355749

// Heliocentric perihelion times of the 2017/2020-epoch heliocentric solutions (kept for tools that
// still use those elements, e.g. as the centre of the SBDB element covariance).
#define OUMUAMUA_TP_JD 2458006.000109129120   // JPL Horizons, epoch JD 2458041.5 solution above
#define BORISOV_TP_JD 2458826.052845906059    // JPL SBDB solution 54 (epoch JD 2458853.5)
#define GAIA_DR2_EPOCH_JD 2457206.375         // J2015.5
#define GAIA_DR3_EPOCH_JD 2457389.0           // J2016.0
// Time of the Gaia DR3 epoch relative to 'Oumuamua's (barycentric) perihelion, days (~ -616.6).
#define OUMUAMUA_T_START_DR3_DAYS (GAIA_DR3_EPOCH_JD - OUMUAMUA_TP_BARY_JD)
// Same for 2I/Borisov (~ -1437.2 days).
#define BORISOV_T_START_DR3_DAYS (GAIA_DR3_EPOCH_JD - BORISOV_TP_BARY_JD)

// 'Oumuamua's asymptotic-velocity uncertainty for solution 7c as published by Bailer-Jones et al.
// 2018 (AJ 156, 205), Table 1: 1-sigma of the barycentric ICRF velocity components (km/s) and their
// correlations (vx,vy), (vx,vz), (vy,vz). Unlike OUMUAMUA_SBDB_COV6 below (osculating elements
// only), this includes the non-gravitational-parameter uncertainty propagated through the whole
// passage, so it is ~40x larger (0.01-0.03 km/s) -- the right ISO uncertainty for the stage-2
// Monte Carlo. Values transcribed from the paper's Table 1 via the arXiv HTML rendering.
static const double OUMUAMUA_7C_VINF_SIGMA_KMS[3] = { 0.0101308, 0.0307732, 0.0149700 };
static const double OUMUAMUA_7C_VINF_CORR[3] = { 0.994705, -0.989469, -0.999101 };

// 2I/Borisov's asymptotic-velocity uncertainty as published by Bailer-Jones et al. 2020 (A&A 634,
// A14), Table 1, for JPL solution 47 -- the solution all of that paper's encounter results use
// (radiant RA 32.79720, Dec 59.44014 deg, v_inf 32.286894 km/s): 1-sigma of the barycentric ICRF
// components (km/s) and correlations (vx,vy), (vx,vz), (vy,vz). Includes the non-gravitational-
// parameter uncertainty of the passage. The pipeline's NOMINAL Borisov state is JPL#54's
// barycentric asymptote (BORISOV_ELEMENTS_BARY_1600; a later fit with more data, radiant within
// ~0.02 deg of #47), so tools/mc_shortlist.c --iso=borisov applies this covariance around #54 --
// the closest published stand-in for #54's own (not published as an asymptote covariance). At
// 0.001-0.004 km/s it is negligible next to any Gaia star's RV error. Transcribed 2026-10-08 from
// the arXiv rendering of Table 1 (arXiv:1912.10213, re-checked 2026-10-08; "we use only solution 47").
static const double BORISOV_47_VINF_SIGMA_KMS[3] = { 0.0024122, 0.0011258, 0.0041225 };
static const double BORISOV_47_VINF_CORR[3] = { 0.9282, 0.9886, 0.9500 };

// 6x6 marginal covariance over (e, q[au], tp[day], node[deg], peri[deg], i[deg]) -- the top-left
// block of JPL SBDB's 9x9 orbit-fit covariance (Solution 7c / Micheli et al. 2018, quoted at epoch
// JD 2458080.5; fetched 2026-10-06 via
// https://ssd-api.jpl.nasa.gov/sbdb.api?sstr=1I&full-prec=true&cov=mat). A1/A2/A3 (non-grav
// parameters) marginalized out by taking the principal submatrix, which is exact. Same epoch
// caveat as tools/oumuamua_velocity_covariance.c: applied around the Horizons-epoch best fit
// above, whose elements agree with SBDB's to ~1-2e-4 relative.
static const double OUMUAMUA_SBDB_COV6[6][6] = {
    { 4.437021723735803e-10, 1.402645684430185e-10, -2.290597192603328e-9, -4.866487927086338e-9,  2.630680330585727e-8,  3.783694880329765e-9 },
    { 1.402645684430185e-10, 4.44023866269166e-11,  -7.704816361638409e-10,-1.55813687037147e-9,   8.310783593401902e-9,  1.153616408668474e-9 },
    {-2.290597192603328e-9, -7.704816361638409e-10,  6.982193721062611e-8,  5.061551619658835e-8, -1.267295484829105e-7,  3.435990038007021e-8 },
    {-4.866487927086338e-9, -1.55813687037147e-9,    5.061551619658835e-8,  6.462880452547438e-8, -2.844751574251313e-7, -1.765593440441589e-8 },
    { 2.630680330585727e-8,  8.310783593401902e-9,  -1.267295484829105e-7, -2.844751574251313e-7,  1.561285075915765e-6,  2.328863207568894e-7 },
    { 3.783694880329765e-9,  1.153616408668474e-9,   3.435990038007021e-8, -1.765593440441589e-8,  2.328863207568894e-7,  8.309282984005662e-8 },
};

#endif
