// The actual per-(star, ISO-realization) search: integrate both bodies backward through the same
// background potential (no mutual gravity between them -- each trajectory is independent in a
// smooth potential, per the project's own stated approach) and find their closest approach.
#ifndef NOSTOS_ENCOUNTER_H
#define NOSTOS_ENCOUNTER_H

#include "integrator.h"

typedef struct {
    double min_distance_kpc;
    double time_at_min_myr; // time since the start of this integration call (same sign as total_time)
    Vec3 star_pos_at_min;   // Galactocentric position of each body when the minimum occurred --
    Vec3 iso_pos_at_min;    // useful for diagnostics/plotting, not needed for the distance itself
} EncounterResult;

// Integrates `star` and `iso` together, step by step, over `total_time` Myr (negative = backward,
// the normal case for this project) using `n_steps` fixed leapfrog steps each, tracking the
// minimum separation between them at every step boundary (including the starting point, step 0).
// Both State structs are advanced in place to their final integrated position -- pass copies if
// the caller still needs the original (e.g. unperturbed) state afterward.
//
// Resolution note: the reported minimum is the smallest separation among the n_steps+1 SAMPLED
// points (start plus every step), not a continuously-refined true minimum -- the true closest
// approach between two sample points could in principle be slightly smaller. For how finely this
// matters in practice (dt vs. typical encounter velocities ~tens of km/s), see
// NOTES_method_selection.md.
EncounterResult find_closest_approach(
    const Potential *pot, State *star, State *iso, double total_time, int n_steps
);

#endif
