// Reusable screening logic for gaia_scan.c, pulled out into its own testable module after two
// real bugs were found in it by ad-hoc testing alone (a past-vs-future sign bug, and a search-
// window-truncation bug -- see NOTES_method_selection.md) with zero automated test coverage at
// the time. Both functions here are now independently tested (tests/test_scan_logic.c).
#ifndef NOSTOS_SCAN_LOGIC_H
#define NOSTOS_SCAN_LOGIC_H

#include <stdbool.h>
#include "coords.h"
#include "encounter.h"
#include "potential.h"

// Cheap analytic straight-line closest approach between two phase-space points (O(1), no
// integration) -- the same exact textbook formula test_encounter.c's own exact case validates.
// Returns false (and leaves *t_star_myr/*d_star_pc unchanged) if the two points have identical
// velocity (degenerate, never happens for real distinct objects, but checked rather than
// dividing by zero).
bool straight_line_closest_approach(
    PhaseSpace a, PhaseSpace b, double *t_star_myr, double *d_star_pc
);

// Confirms a straight-line candidate with the real encounter search (full gravity), searching
// backward from "now" (t=0) with an initial window of initial_window_myr. If the found minimum
// falls within 2 integration steps of either end of the search window (meaning it might just be
// the window's own edge, not a true minimum -- the bug this function exists to avoid repeating),
// doubles the window and retries, up to max_attempts times. Returns the final (possibly widened)
// window actually used in *final_window_myr, so a caller can tell whether/how much widening
// happened.
EncounterResult gravity_confirm(
    const Potential *pot, PhaseSpace star, PhaseSpace iso,
    double initial_window_myr, int steps_per_myr, int max_attempts,
    double *final_window_myr
);

#endif
