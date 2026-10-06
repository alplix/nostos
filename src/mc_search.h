// Stage-1 (exhaustive, cheap) Monte Carlo encounter search: one star's
// single best-fit 6D state against MC_N_ISO_CLONES Monte Carlo realizations
// of an interstellar object's velocity uncertainty. This is the pass meant
// to run across the ENTIRE Gaia catalog (tens of millions of stars) in the
// distributed BOINC work units -- see NOTES_method_selection.md's "Monte
// Carlo sampling scheme" section for why this is a two-stage design (full
// star-side Monte Carlo too, matching Bailer-Jones et al. 2019 exactly, is
// reserved for a stage-2 follow-up on whatever stars survive stage 1 with a
// close-enough approach -- not implemented here, since it only ever runs on
// a short, already-shortlisted list rather than the full catalog).
#ifndef NOSTOS_MC_SEARCH_H
#define NOSTOS_MC_SEARCH_H

#include <stdint.h>
#include "coords.h"
#include "encounter.h"
#include "montecarlo.h"

// Bailer-Jones et al. 2019's own number (arXiv:1911.02473) -- see
// NOTES_method_selection.md, not picked arbitrarily.
#define MC_N_ISO_CLONES 100

typedef struct {
    double min_distance_kpc;      // smallest closest-approach distance across all clones
    double mean_min_distance_kpc; // average closest-approach distance across clones -- a
                                   // large gap between this and min_distance_kpc signals a
                                   // wide/uncertain result, not a tight one
    double time_at_min_myr;       // time (since integration start) of the single closest
                                   // pairing found, i.e. from the clone that produced
                                   // min_distance_kpc
    int n_clones_run;             // always MC_N_ISO_CLONES unless the covariance matrix was
                                   // rejected as non-positive-definite, in which case 0 (and
                                   // every other field is left at its zero-initialized value --
                                   // callers should check this before trusting the result)
} McEncounterSummary;

// iso_mean: the ISO's best-fit starting state (position treated as exactly
// known -- see this file's top comment on why position isn't resampled).
// vel_cov_kpcmyr2: the 3x3 covariance matrix of the ISO's velocity
// uncertainty, in (kpc/Myr)^2 -- convert from whatever units the published
// orbital solution's covariance is quoted in (typically km/s or km^2/s^2)
// using units.h's helpers before calling this.
// star: the Gaia star's single best-fit 6D state (no Monte Carlo on the
// star side in this stage-1 pass).
// seed: must differ per (ISO, star) pair for genuinely independent draws
// across the catalog -- e.g. derived from the star's Gaia source_id -- but
// must also be exactly REPRODUCIBLE from the same inputs (not time-based),
// since BOINC quorum validation re-runs the same work unit on a different
// host and expects the same answer.
McEncounterSummary mc_search_stage1(
    const Potential *pot, PhaseSpace iso_mean, const double vel_cov_kpcmyr2[3][3],
    PhaseSpace star, double total_time_myr, int n_steps, uint64_t seed
);

#endif
