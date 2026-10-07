// Second interstellar object for the "Pipeline validation" roadmap item (the first,
// 1I/'Oumuamua, is in test_pipeline_oumuamua_hip3757.c -- read that file's top comment first for
// the full story of the inbound/outbound-asymptote bug this file's history also went through).
// Target: Ross 573, the single closest encounter Bailer-Jones, Farnocchia, Ye, Meech & Micheli
// 2020 (A&A, arXiv:1912.10213) found for 2I/Borisov across their entire 7.4-million-star Gaia
// DR2 search -- their headline result for this object, not a hand-picked easy case:
//   t = -909 kyr (median of 2000 Monte Carlo surrogates, 90% CI -937..-882 kyr)
//   d = 0.068 pc (median, 90% CI 0.053..0.091 pc)
//   v = 22.6 km/s (median, 90% CI 22.0..23.3 km/s)
// (That paper reports Monte Carlo medians, not a single best-fit trajectory, same caveat as
// test_pipeline_oumuamua_hip3757.c already notes for the 'Oumuamua comparisons.)
//
// Ross 573's identity (Gaia DR2 source_id 5162123155863791744) confirmed independently via
// SIMBAD: parallax 47.5505 mas and proper motion (54.471, -304.357) mas/yr there match this
// test's Gaia-archive-queried values (47.536, (54.430, -304.389)) to 3-4 significant figures.
//
// 2I/Borisov's orbital elements: real JPL Small-Body Database values (Solution 54, epoch
// JD 2458853.5, fetched 2026-10-06 via
// https://ssd-api.jpl.nasa.gov/sbdb.api?sstr=2I&full-prec=true) -- same source and method as
// 'Oumuamua's elements (test_horizons_crosscheck.c), just a different interstellar object.
//
// Evaluation epoch: 3.0 years BEFORE perihelion (Borisov's perihelion was JD 2458826.05, i.e.
// 2019-12-08) -- the INBOUND asymptote, not the outbound one. This file originally used +3.0
// years (outbound) and found a ~200x-too-large closest-approach distance; see
// test_pipeline_oumuamua_hip3757.c's top comment for why that's wrong (the outbound asymptote's
// direction differs from the correct inbound one by the orbit's full hyperbolic turning angle,
// 2*asin(1/e) = 34.7 degrees for this orbit, and this backward integration never re-simulates
// the Sun's actual deflection to correct for that). The 3.0-year MAGNITUDE itself was chosen by
// checking where the computed velocity converges to its asymptotic (hyperbolic-excess) value --
// at 3.0 years the computed speed is within ~0.03% of its value at 10.0 years (233.3 km/s at
// 0.5yr vs 235.30 km/s at 3.0yr vs 235.46 km/s at 10.0yr, all on the outbound branch where this
// check was originally done; by the time-symmetry of the two-body Kepler problem the inbound
// branch converges at the same rate), and matches this orbit's analytic hyperbolic-excess speed
// sqrt(-GM_sun/a) = 32.28 km/s (heliocentric) well within that same margin.
//
// CORRECTION, 2026-10-08: that convergence claim was wrong. The 0.03% figure compared
// GALACTOCENTRIC speeds (235 km/s, dominated by the Sun's own motion), which hides the
// heliocentric error: at 3.0 yr before perihelion (r = 22 AU) the heliocentric speed is
// 33.50 km/s, 3.8% (1.22 km/s) above v_inf = 32.28 km/s (vis-viva v^2 = v_inf^2 + 2GM/r converges
// only like 1/r). Over the ~0.9 Myr look-back that excess is ~1 pc of along-track error -- the
// whole "residual ~15x distance gap" (1.06 pc vs 0.068 pc). Now started from the EXACT inbound
// asymptote at the Gaia DR2 epoch: t=-0.910 Myr, d=0.074 pc, v=22.62 km/s from the heliocentric
// 2020-epoch elements, and t=-0.9096, d=0.0802, v=22.63 from JPL's barycentric 1600 elements (the
// pipeline convention, src/iso_solutions.h) -- inside the paper's own 90% intervals for all three
// quantities. With the paper's own inputs (its Table 1 asymptote for solution 47, parallax
// zero-point) the nominal is 0.084 pc and a 2000-surrogate MC median 0.086 [0.075,0.102] pc vs the
// published 0.068 [0.053,0.091]: ~1 sigma high. The leftover ~0.015 pc is below the ~0.025-0.03 pc
// by which the Galactic tide alone moves this tiny miss vector (zero potential: 0.064 pc), so it is
// at the level of potential-model / implementation details; see NOTES_method_selection.md.
#include <stdio.h>
#include <math.h>
#include "../src/coords.h"
#include "../src/orbital_elements.h"
#include "../src/encounter.h"
#include "../src/potential.h"
#include "../src/units.h"
#include "../src/iso_solutions.h"

static int g_failures = 0;

static void check(int ok, const char *what) {
    if (!ok) { printf("  FAIL: %s\n", what); g_failures++; }
    else { printf("  ok:   %s\n", what); }
}

int main(void) {
    printf("test_pipeline_borisov_ross573:\n");

    // Exact inbound asymptote of JPL's barycentric 1600 elements for the same solution (JPL#54,
    // src/iso_solutions.h) with the total solar-system GM, placed at the Gaia DR2 epoch (J2015.5)
    // of the star data below -- the pipeline's own starting-state convention.
    double t_inbound_myr = (GAIA_DR2_EPOCH_JD - BORISOV_TP_BARY_JD) * NOSTOS_MYR_PER_DAY;
    PhaseSpace iso_ps;
    bool ok = orbital_elements_inbound_asymptote_gm(BORISOV_ELEMENTS_BARY_1600, NOSTOS_GM_SSB,
                                                    t_inbound_myr, &solar_params_default, &iso_ps);
    check(ok, "2I/Borisov's orbital elements convert to a valid Galactocentric state");

    // Ross 573, real Gaia DR2 astrometry, source_id 5162123155863791744 (queried 2026-10-06 from
    // gaiadr2.gaia_source, identity cross-confirmed via SIMBAD -- see this file's top comment).
    Astrometry ross573 = {
        .ra_deg = 52.33270555516616,
        .dec_deg = -11.679678598930318,
        .parallax_mas = 47.536273832249684,
        .pmra_masyr = 54.42974631770208,
        .pmdec_masyr = -304.38948045019396,
        .rv_kms = 13.060251980364459,
    };
    PhaseSpace star_ps;
    astrometry_to_galactocentric(ross573, &solar_params_default, &star_ps);

    // Separation at the Gaia epoch between Ross 573 and the point on Borisov's inbound asymptote.
    // Printed only as a sanity-check starting distance for the search below.
    double d0_pc = v3_norm(v3_sub(star_ps.pos, iso_ps.pos)) * 1000.0;
    printf("  separation at the inbound starting point: %.2f pc\n", d0_pc);

    Potential mw = milkyway_potential(&milkyway_potential_default);
    State star_state = {star_ps.pos, star_ps.vel};
    State iso_state = {iso_ps.pos, iso_ps.vel};
    double search_window_myr = 2.0; // reported encounter at -0.909 Myr; same margin convention
                                     // as test_pipeline_oumuamua_hip3757.c's HIP 3757 case.
                                     // Checked manually at 5.0 Myr too (same HD 292249-style
                                     // widening sanity check) -- found the same minimum, so this
                                     // isn't a window-truncation artifact.
    int n_steps = (int)(search_window_myr * 100000); // ~10 yr steps
    EncounterResult r = find_closest_approach(&mw, &star_state, &iso_state, -search_window_myr, n_steps);

    double min_d_pc = r.min_distance_kpc * 1000.0;
    Vec3 rel_vel = v3_sub(r.star_vel_at_min, r.iso_vel_at_min);
    double rel_speed_kms = nostos_kpcmyr_to_kms(v3_norm(rel_vel));

    printf("  found: t=%.4f Myr, distance=%.4f pc, rel. speed=%.2f km/s\n",
           r.time_at_min_myr, min_d_pc, rel_speed_kms);
    printf("  Bailer-Jones et al. 2020: t=-0.909 Myr, distance=0.068 pc, rel. speed=22.6 km/s\n");

    // Asserted against the paper's own 90% intervals (t -937..-882 kyr, d 0.053..0.091 pc,
    // v 22.0..23.3 km/s). Before the 2026-10-08 exact-asymptote fix: d=1.06 pc, t=-0.893, v=23.01 --
    // the distance check would fail by a factor ~12 (time and speed were already inside).
    check(r.time_at_min_myr < 0, "closest approach found in the past (not the future)");
    check(r.time_at_min_myr > -0.937 && r.time_at_min_myr < -0.882,
          "closest-approach time inside Bailer-Jones et al. 2020's 90% interval");
    check(rel_speed_kms > 22.0 && rel_speed_kms < 23.3,
          "relative speed inside Bailer-Jones et al. 2020's 90% interval");
    check(min_d_pc > 0.053 && min_d_pc < 0.091,
          "closest-approach distance inside Bailer-Jones et al. 2020's 90% interval");
    printf("\n");

    if (g_failures) {
        printf("%d check(s) FAILED\n", g_failures);
        return 1;
    }
    printf("Validation passed (time, distance and speed all inside Bailer-Jones et al. 2020's "
           "90%% intervals).\n");
    return 0;
}
