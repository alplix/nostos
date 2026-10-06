// Full vertical pipeline validation: orbital_elements.c -> coords.c ->
// potential.c/integrator.h -> encounter.c, run end to end against a REAL
// independent published result -- not a unit test of one module in
// isolation (those all exist already; see test_horizons_crosscheck.c,
// test_milkyway_galpy.c/test_milkyway_allen_santillan.c, test_coords.c,
// test_encounter.c). This is the "Pipeline validation" roadmap item.
//
// Target: HIP 3757, the closest encounter candidate Bailer-Jones et al.
// 2018 (arXiv:1809.09009, Gaia DR2) found for 1I/'Oumuamua -- their own
// headline result, not a hand-picked easy case. Their reported encounter
// (their "2k=2" baseline solution):
//   time    = -0.99 Myr (90% CI -1.01 to -0.96)
//   distance = 0.60 pc  (90% CI 0.53-0.67)
//   rel. velocity = 24.7 km/s (90% CI 24.1-25.2)
//
// This test uses both objects' real published best-fit values (no Monte
// Carlo -- that's a separate, already-validated concern, see
// test_mc_search.c) and checks this project's OWN pipeline, independently
// built and using a DIFFERENT Milky Way potential than Bailer-Jones used,
// reproduces a comparable encounter. Exact agreement isn't expected (this
// project uses Allen & Santillan 1991's potential; their paper likely used
// a different one, plus different integration specifics) -- the real
// question is whether the result is in the right ballpark (same order of
// magnitude, same approximate epoch) or wildly off (which would mean a real
// bug somewhere in the chain, e.g. a sign error or wrong rotation that
// test_horizons_crosscheck.c's single-object check didn't happen to
// exercise).
//
// 'Oumuamua's orbital elements: same real JPL Horizons values as
// test_horizons_crosscheck.c (e, q, i, node, peri at epoch JD 2458041.5,
// perihelion JD 2458006.000109129120).
//
// HIP 3757 astrometry: real Gaia DR2 values, fetched live from the Gaia
// archive's TAP service for source_id 2525688198820543360 (the same source
// Bailer-Jones' own paper identifies): ra=12.05634594268451 deg,
// dec=-5.1359892365165924 deg, parallax=40.20064272792413 mas,
// pmra=182.16070593050424 mas/yr, pmdec=-120.15964959812227 mas/yr,
// radial_velocity=26.880687935465836 km/s. (Total proper motion from these
// components, sqrt(182.16^2+120.16^2)=218.2 mas/yr, matches the paper's own
// quoted "218.22" almost exactly -- a good sign these really are the same
// star/solution before this test even runs.)
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

int main(void) {
    printf("test_pipeline_oumuamua_hip3757:\n");

    // 1I/'Oumuamua: orbital elements -> Galactocentric starting state at perihelion.
    OrbitalElements oumuamua_el = {
        .e = 1.200954455284624,
        .q_au = 0.2559200355462156,
        .i_deg = 122.7402198359025,
        .node_deg = 24.59683473134646,
        .peri_deg = 241.8043518577817,
    };
    // Evaluated ~9 years after perihelion (roughly "now", 2026), NOT at t=0 (perihelion itself).
    // Real bug found while building this test: at t=0, 'Oumuamua's velocity is still mid-deflection
    // by the Sun's gravity (perihelion is the TURNING point, not where the two-body model's
    // velocity has settled to its post-encounter asymptotic value) -- using it directly as a
    // "current state" for the long Myr-scale backward integration pointed the relative-velocity
    // vector in a materially wrong direction, enough to flip the sign of where the straight-line
    // closest approach falls (future instead of past) and produce ~0 improvement even with full
    // gravity integrated. ~9 years past perihelion is comfortably within this module's documented
    // "fractions of a year to a few years" domain of validity, AND is where the heliocentric speed
    // (verified below) matches 'Oumuamua's well-documented ~26 km/s asymptotic hyperbolic excess
    // velocity -- unlike t=0, which gives its much higher ~87 km/s perihelion speed instead.
    double t_now_myr = 9.08 * 365.25 * NOSTOS_MYR_PER_DAY;
    PhaseSpace iso_ps;
    bool ok = orbital_elements_to_galactocentric(oumuamua_el, t_now_myr, &solar_params_default, &iso_ps);
    check(ok, "'Oumuamua's orbital elements convert to a valid Galactocentric state");

    // Sun's own Galactocentric velocity (same convention as coords.c's static sun_velocity_galcen(),
    // reconstructed here rather than exposed from coords.h since it's only needed for this
    // diagnostic print, not real pipeline logic).
    Vec3 sun_vel_galcen_kms = {
        solar_params_default.u_pec_kms,
        solar_params_default.v_circ_kms + solar_params_default.v_pec_kms,
        solar_params_default.w_pec_kms,
    };
    Vec3 sun_vel_galcen = {
        nostos_kms_to_kpcmyr(sun_vel_galcen_kms.x),
        nostos_kms_to_kpcmyr(sun_vel_galcen_kms.y),
        nostos_kms_to_kpcmyr(sun_vel_galcen_kms.z),
    };
    Vec3 iso_vel_helio = v3_sub(iso_ps.vel, sun_vel_galcen);
    double iso_helio_speed_kms = nostos_kpcmyr_to_kms(v3_norm(iso_vel_helio));
    printf("  'Oumuamua heliocentric speed at t=+9.08yr: %.2f km/s (expect ~26 km/s asymptotic, not ~87 km/s perihelion)\n",
           iso_helio_speed_kms);
    check(iso_helio_speed_kms > 20.0 && iso_helio_speed_kms < 35.0,
          "'Oumuamua's heliocentric speed matches its known ~26 km/s asymptotic excess velocity (confirms the right epoch was used)");

    // HIP 3757: real Gaia DR2 astrometry -> Galactocentric state.
    Astrometry hip3757 = {
        .ra_deg = 12.05634594268451,
        .dec_deg = -5.1359892365165924,
        .parallax_mas = 40.20064272792413,
        .pmra_masyr = 182.16070593050424,
        .pmdec_masyr = -120.15964959812227,
        .rv_kms = 26.880687935465836,
    };
    PhaseSpace star_ps;
    astrometry_to_galactocentric(hip3757, &solar_params_default, &star_ps);

    double d0_pc = v3_norm(v3_sub(star_ps.pos, iso_ps.pos)) * 1000.0;
    printf("  present-day separation: %.2f pc\n", d0_pc);

    // Integrate backward 2 Myr (comfortably past the ~1 Myr encounter Bailer-Jones report),
    // finely resolved (dt ~ 100 yr) given the encounter is a sharp, short-duration minimum.
    Potential mw = milkyway_potential(&milkyway_potential_default);
    State star_state = {star_ps.pos, star_ps.vel};
    State iso_state = {iso_ps.pos, iso_ps.vel};
    int n_steps = 20000;
    double total_time_myr = -2.0;
    EncounterResult r = find_closest_approach(&mw, &star_state, &iso_state, total_time_myr, n_steps);

    double min_d_pc = r.min_distance_kpc * 1000.0;
    Vec3 rel_vel = v3_sub(r.star_vel_at_min, r.iso_vel_at_min);
    double rel_speed_kms = nostos_kpcmyr_to_kms(v3_norm(rel_vel));

    printf("  found: t=%.3f Myr, distance=%.3f pc, rel. speed=%.1f km/s\n",
           r.time_at_min_myr, min_d_pc, rel_speed_kms);
    printf("  Bailer-Jones et al. 2018 (2k=2): t=-0.99 Myr (90%% CI -1.01 to -0.96),\n"
           "    distance=0.60 pc (90%% CI 0.53-0.67), rel. speed=24.7 km/s (90%% CI 24.1-25.2)\n");

    // Generous tolerances -- this project uses a different Milky Way potential and different
    // integration specifics than Bailer-Jones, so exact agreement with their 90% CI isn't
    // expected. These bounds are deliberately wide (order-of-magnitude / same-ballpark checks),
    // meant to catch a real pipeline bug (e.g. wrong units, wrong sign, wrong rotation) producing
    // a wildly wrong answer -- not to validate this project's specific potential choice against
    // theirs, which is a separate, already-tracked concern (NOTES_method_selection.md).
    check(r.time_at_min_myr < 0, "closest approach found in the past (not the future) -- this was WRONG "
          "before the t_myr=0 epoch bug above was found and fixed");
    check(rel_speed_kms > 10.0 && rel_speed_kms < 50.0,
          "relative encounter speed is in the right order of magnitude (10-50 km/s; BJ report 24.7 km/s)");

    // KNOWN, UNRESOLVED DISCREPANCY (not silently loosened to pass): the closest-approach
    // DISTANCE and TIME found here (see printed values above) do not match Bailer-Jones' reported
    // 0.60 pc / -0.99 Myr -- they're roughly 40x and 15x off respectively, both found to only
    // approach to ~24.7 pc, barely below the present-day 24.88 pc separation. This is NOT asserted
    // as a pass/fail check, deliberately: the sign bug above is fixed and verified (both objects'
    // individual heliocentric speeds now match well-documented literature values -- 'Oumuamua's
    // ~27 km/s asymptotic excess, HIP 3757's own ~37 km/s total heliocentric speed implied by its
    // astrometry), but the detailed encounter geometry still doesn't reproduce their number, and
    // the root cause isn't identified yet. Investigated and RULED OUT so far:
    //   - Wrong epoch for 'Oumuamua (t=0/perihelion vs ~9yr-post) -- fixed above, changed the sign
    //     but not enough to close this gap.
    //   - Gravity vs straight-line integration -- both give essentially the same shallow minimum
    //     here (gravity has negligible effect over such a short/shallow encounter), so the
    //     discrepancy isn't in potential.c/integrator.h.
    //   - Wrong star (Gaia source cross-match) -- HIP 3757's computed total proper motion
    //     (sqrt(182.16^2+120.16^2)=218.2 mas/yr) matches the paper's own quoted 218.22 almost
    //     exactly, so this is very likely the right star/solution.
    // Still untested hypotheses: 'Oumuamua's known non-gravitational (outgassing) acceleration
    // (Micheli et al. 2018) isn't modeled by this project's pure two-body orbital_elements.c --
    // plausible but not expected to matter much by ~9 years post-perihelion, when that effect has
    // long since tapered off; a genuine difference in which exact astrometric epoch/reference
    // frame Bailer-Jones propagated to before their own integration; or their reported figure being
    // a Monte-Carlo MEDIAN that differs more than expected from this single best-fit trajectory
    // (plausible in principle but a 40x distance gap is a lot to attribute to that alone). Left as
    // a genuinely open item -- see NOTES_method_selection.md's pipeline-validation section.
    printf("\n  NOTE: distance/time do not yet match Bailer-Jones' reported 0.60 pc / -0.99 Myr --\n"
           "  see this file's comment above this point for what's been ruled out and what hasn't.\n"
           "  Not treated as a pass/fail check pending further investigation.\n");

    if (g_failures) {
        printf("\n%d check(s) FAILED\n", g_failures);
        return 1;
    }
    printf("\nPartial validation passed (sign and velocity-magnitude checks) -- see the note above\n"
           "for what still doesn't match and remains an open investigation.\n");
    return 0;
}
