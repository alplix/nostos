// Work-unit input parsing for main_nostos_boinc.cpp, pulled into its own plain-C module (not
// C++, no std::vector) so it can be unit-tested with the same plain gcc toolchain as the rest of
// this project, without needing BOINC headers/libs locally -- see tests/test_nostos_io.c. The
// BOINC-specific parts (boinc_fopen, resolve_filename, checkpointing) stay in
// main_nostos_boinc.cpp itself; this module only knows about a plain FILE*.
//
// Input file format (plain text, one work unit's data) -- see main_nostos_boinc.cpp's own top
// comment for the full field-by-field description:
//   Line 1: e q_au i_deg node_deg peri_deg
//   Line 2: cov00 cov01 cov02 cov11 cov12 cov22 (upper triangle of the 3x3 velocity covariance)
//   Line 3: seed_base
//   Line 4: total_time_myr n_steps distance_threshold_pc
//   Remaining lines: source_id ra_deg dec_deg parallax_mas pmra_masyr pmdec_masyr rv_kms
#ifndef NOSTOS_IO_H
#define NOSTOS_IO_H

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include "coords.h"
#include "orbital_elements.h"

typedef struct {
    uint64_t source_id;
    Astrometry astrometry;
} StarRecord;

// On success, *stars_out is a malloc'd array of *n_stars_out StarRecords -- caller owns it and
// must free() it. On failure, *stars_out/*n_stars_out are left unset; caller must not free them.
bool parse_input_stream(
    FILE *f, OrbitalElements *el, double vel_cov[3][3], uint64_t *seed_base,
    double *total_time_myr, int *n_steps, double *distance_threshold_pc,
    StarRecord **stars_out, size_t *n_stars_out
);

// Mixes a per-work-unit seed with a per-star source_id into one independent-looking 64-bit seed
// (splitmix64's own mixing step, see montecarlo.c) -- used so neighboring source_ids (often
// numerically close within one Gaia catalog chunk) don't produce visibly-correlated Monte Carlo
// draws between stars.
uint64_t combine_seed(uint64_t seed_base, uint64_t source_id);

#endif
