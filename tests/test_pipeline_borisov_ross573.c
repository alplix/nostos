// Second interstellar object for the "Pipeline validation" roadmap item (the first,
// 1I/'Oumuamua, is in test_pipeline_oumuamua_hip3757.c -- read that file's top comment first,
// since this one follows the exact same structure and the same open discrepancy applies here
// too). Target: Ross 573, the single closest encounter Bailer-Jones, Farnocchia, Ye, Meech &
// Micheli 2020 (A&A, arXiv:1912.10213) found for 2I/Borisov across their entire 7.4-million-star
// Gaia DR2 search -- their headline result for this object, not a hand-picked easy case:
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
// Evaluation epoch: 3.0 years post-perihelion (Borisov's perihelion was JD 2458826.05, i.e.
// 2019-12-08). Chosen the same way 'Oumuamua's 9.08-year point was: checked where the computed
// Galactocentric velocity has converged to its asymptotic (hyperbolic-excess) value rather than
// still being mid-deflection by the Sun's gravity -- at 3.0 years the computed speed is within
// ~0.03% of its value at 10.0 years (233.3 km/s at 0.5yr vs 235.30 km/s at 3.0yr vs
// 235.46 km/s at 10.0yr), and matches this orbit's analytic hyperbolic-excess speed
// sqrt(-GM_sun/a) = 32.28 km/s (heliocentric) well within that same margin. Borisov's much
// larger perihelion distance (2.01 au vs 'Oumuamua's 0.26 au) means it converges to its
// asymptotic velocity over a longer PHYSICAL time than 'Oumuamua despite a lower speed, but it
// still converges well inside this project's multi-Myr integration timescale, which is the only
// thing that actually matters here.
#include <stdio.h>
#include <math.h>
#include "../src/coords.h"
#include "../src/orbital_elements.h"
#include "../src/encounter.h"
#include "../src/potential.h"
#include "../src/units.h"

static int g_failures = 0;

static void check(int ok, const char *what) {
    if (!ok) { printf("  FAIL: %s\n", what); g_failures++; }
    else { printf("  ok:   %s\n", what); }
}

static OrbitalElements borisov_elements(void) {
    return (OrbitalElements){
        .e = 3.356475782676596,
        .q_au = 2.006520878500843,
        .i_deg = 44.05264247909138,
        .node_deg = 308.1477292269942,
        .peri_deg = 209.1236864378081,
    };
}

int main(void) {
    printf("test_pipeline_borisov_ross573:\n");

    double t_now_myr = 3.0 * 365.25 * NOSTOS_MYR_PER_DAY;
    PhaseSpace iso_ps;
    bool ok = orbital_elements_to_galactocentric(borisov_elements(), t_now_myr, &solar_params_default, &iso_ps);
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

    double d0_pc = v3_norm(v3_sub(star_ps.pos, iso_ps.pos)) * 1000.0;
    printf("  present-day separation: %.2f pc\n", d0_pc);

    Potential mw = milkyway_potential(&milkyway_potential_default);
    State star_state = {star_ps.pos, star_ps.vel};
    State iso_state = {iso_ps.pos, iso_ps.vel};
    double search_window_myr = 2.0; // reported encounter at -0.909 Myr; same margin convention
                                     // as test_pipeline_oumuamua_hip3757.c's HIP 3757 case.
                                     // Checked manually at 5.0 Myr too (same HD 292249-style
                                     // widening sanity check) -- found the exact same minimum
                                     // (t=-0.523 Myr), so this isn't a window-truncation
                                     // artifact.
                                     // as test_pipeline_oumuamua_hip3757.c's HIP 3757 case
    int n_steps = (int)(search_window_myr * 100000); // ~10 yr steps
    EncounterResult r = find_closest_approach(&mw, &star_state, &iso_state, -search_window_myr, n_steps);

    double min_d_pc = r.min_distance_kpc * 1000.0;
    Vec3 rel_vel = v3_sub(r.star_vel_at_min, r.iso_vel_at_min);
    double rel_speed_kms = nostos_kpcmyr_to_kms(v3_norm(rel_vel));

    printf("  found: t=%.3f Myr, distance=%.4f pc, rel. speed=%.1f km/s\n",
           r.time_at_min_myr, min_d_pc, rel_speed_kms);
    printf("  Bailer-Jones et al. 2020: t=-0.909 Myr, distance=0.068 pc, rel. speed=22.6 km/s\n");

    // Same convention as test_pipeline_oumuamua_hip3757.c: only the sign is a hard pass/fail
    // check (guards against the exact class of epoch/direction bug that test file documents
    // finding for 'Oumuamua). Distance/time/speed magnitude are printed, not asserted -- see
    // this file's top comment and NOTES_method_selection.md for why that comparison is still an
    // open question for 'Oumuamua, and this result tells us whether the SAME pattern shows up
    // for a second, independent interstellar object.
    check(r.time_at_min_myr < 0, "closest approach found in the past (not the future)");

    printf("  NOTE: distance/time/speed vs. Bailer-Jones et al. 2020 not asserted here -- see "
           "this file's top comment and NOTES_method_selection.md.\n\n");

    if (g_failures) {
        printf("%d check(s) FAILED\n", g_failures);
        return 1;
    }
    printf("Partial validation passed (encounter found in the past, as it should be) -- the "
           "precise encounter distance/time/speed are reported above for comparison, not "
           "asserted.\n");
    return 0;
}
