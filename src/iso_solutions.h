// Shared published orbit solutions for interstellar objects, so the numbers live in ONE place for
// any tool that needs both a best-fit element set and its covariance (first user:
// tools/mc_shortlist.c; tools/oumuamua_velocity_covariance.c takes its covariance from here too).
// The pipeline tests and src/gaia_scan.c still carry their own literal copies of the best-fit
// elements (identical digits) -- deliberately left alone rather than churned.
#ifndef NOSTOS_ISO_SOLUTIONS_H
#define NOSTOS_ISO_SOLUTIONS_H

#include "orbital_elements.h"

// 1I/'Oumuamua best-fit elements: JPL Horizons, epoch JD 2458041.5 -- the exact set src/gaia_scan.c
// screened the whole catalog with (so a Monte Carlo run centred on it is directly comparable to
// that scan's nominal numbers).
static const OrbitalElements OUMUAMUA_ELEMENTS_HORIZONS = {
    .e = 1.200954455284624, .q_au = 0.2559200355462156,
    .i_deg = 122.7402198359025, .node_deg = 24.59683473134646, .peri_deg = 241.8043518577817,
};

// Inbound-asymptote evaluation point, years before perihelion (see NOTES_method_selection.md,
// "Found the actual bug").
#define OUMUAMUA_T_EVAL_YEARS (-9.08)

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
