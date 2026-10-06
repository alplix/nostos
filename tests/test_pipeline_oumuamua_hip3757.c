// Full vertical pipeline validation: orbital_elements.c -> coords.c ->
// potential.c/integrator.h -> encounter.c, run end to end against REAL
// independent published results -- not a unit test of one module in
// isolation (those all exist already; see test_horizons_crosscheck.c,
// test_coords_astropy_crosscheck.c, test_milkyway_galpy.c/
// test_milkyway_allen_santillan.c, test_coords.c, test_encounter.c). This
// is the "Pipeline validation" roadmap item.
//
// Targets: the two most distinct candidates Bailer-Jones et al. 2018
// (arXiv:1809.09009, Gaia DR2) found for 1I/'Oumuamua's origin -- their own
// headline results, not hand-picked easy cases:
//   HIP 3757  (closest approach, highest relative speed): t=-0.99 Myr,
//             d=0.60 pc, v=24.7 km/s (their "2k=2" solution)
//   HD 292249 (more distant, lowest relative speed, their most PLAUSIBLE
//             candidate): t=-3.765 Myr, d=1.60 pc, v=10.7 km/s
// Both stars' identities confirmed independently via SIMBAD (matches the
// Gaia DR2 source_id in each case) and, for HIP 3757, via its total proper
// motion (sqrt(182.16^2+120.16^2)=218.2 mas/yr matching the paper's own
// quoted 218.22 to 4 significant figures).
//
// Status: this is a genuinely OPEN validation, documented honestly rather
// than adjusted to pass. BOTH stars show the exact same pattern -- found
// encounter in the right temporal direction (past) and right
// order-of-magnitude relative speed, but a MUCH shallower/earlier minimum
// than reported (tens of pc instead of ~1 pc). Despite extensive
// investigation (see NOTES_method_selection.md's full writeup), every
// individual pipeline component has now been independently validated and
// checks out:
//   - 'Oumuamua's orbital_elements.c propagation: matches real JPL Horizons
//     data to ~6e-4 relative error even 9 years post-perihelion (not just
//     the ~1e-6 at +36 days test_horizons_crosscheck.c already checks).
//   - HIP 3757's astrometry_to_galactocentric() output: matches astropy's
//     independent ICRS->Galactic implementation to ~1e-7
//     (test_coords_astropy_crosscheck.c).
//   - Full gravity (Allen & Santillan potential) vs straight-line give
//     essentially the SAME shallow minimum for both stars, even searched
//     out to -6 Myr for HD 292249 (well past its reported -3.765 Myr) --
//     gravity genuinely isn't the differentiator here.
//   - 'Oumuamua's known non-gravitational (outgassing) acceleration
//     (Micheli et al. 2018, confirmed present in the JPL Horizons solution
//     used here) is quantitatively too small (~0.05 km/s velocity change
//     from its few remaining active months after this test's reference
//     epoch) to explain a 15-40x distance/time discrepancy.
// The systematic nature (same qualitative mismatch for two very different
// stars) with every individual component independently validated is itself
// a real, useful finding -- it rules out a component-level bug and points
// at something about the OVERALL comparison (Bailer-Jones' exact method,
// reference epoch, or Monte-Carlo-median-vs-single-best-fit) not yet
// understood. Left open rather than papered over.
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

// 'Oumuamua's orbital elements, shared by both cases -- real JPL Horizons values, same as
// test_horizons_crosscheck.c (e, q, i, node, peri at epoch JD 2458041.5, perihelion
// JD 2458006.000109129120).
static OrbitalElements oumuamua_elements(void) {
    return (OrbitalElements){
        .e = 1.200954455284624,
        .q_au = 0.2559200355462156,
        .i_deg = 122.7402198359025,
        .node_deg = 24.59683473134646,
        .peri_deg = 241.8043518577817,
    };
}

static void run_case(
    const char *star_name, Astrometry star_astrometry,
    double bj_time_myr, double bj_distance_pc, double bj_speed_kms,
    double search_window_myr
) {
    printf("%s:\n", star_name);

    // Evaluated ~9 years after perihelion (roughly "now"), NOT at t=0 (perihelion itself). Real
    // bug found while building this test (see NOTES_method_selection.md): at t=0, 'Oumuamua's
    // velocity is still mid-deflection by the Sun's gravity, not yet at its post-encounter
    // asymptotic value -- using it directly pointed the relative-velocity vector in a materially
    // wrong direction, enough to flip whether the closest approach falls in the future or past.
    double t_now_myr = 9.08 * 365.25 * NOSTOS_MYR_PER_DAY;
    PhaseSpace iso_ps;
    bool ok = orbital_elements_to_galactocentric(oumuamua_elements(), t_now_myr, &solar_params_default, &iso_ps);
    check(ok, "'Oumuamua's orbital elements convert to a valid Galactocentric state");

    PhaseSpace star_ps;
    astrometry_to_galactocentric(star_astrometry, &solar_params_default, &star_ps);

    double d0_pc = v3_norm(v3_sub(star_ps.pos, iso_ps.pos)) * 1000.0;
    printf("  present-day separation: %.2f pc\n", d0_pc);

    Potential mw = milkyway_potential(&milkyway_potential_default);
    State star_state = {star_ps.pos, star_ps.vel};
    State iso_state = {iso_ps.pos, iso_ps.vel};
    int n_steps = (int)(search_window_myr * 100000); // ~10 yr steps
    EncounterResult r = find_closest_approach(&mw, &star_state, &iso_state, -search_window_myr, n_steps);

    double min_d_pc = r.min_distance_kpc * 1000.0;
    Vec3 rel_vel = v3_sub(r.star_vel_at_min, r.iso_vel_at_min);
    double rel_speed_kms = nostos_kpcmyr_to_kms(v3_norm(rel_vel));

    printf("  found: t=%.3f Myr, distance=%.3f pc, rel. speed=%.1f km/s\n",
           r.time_at_min_myr, min_d_pc, rel_speed_kms);
    printf("  Bailer-Jones et al. 2018: t=%.3f Myr, distance=%.2f pc, rel. speed=%.1f km/s\n",
           bj_time_myr, bj_distance_pc, bj_speed_kms);

    // Only the sign is asserted as pass/fail -- this directly guards against a regression of the
    // real bug this test found and fixed (evaluating 'Oumuamua at the wrong epoch flipped this).
    // Distance, time, and speed magnitude are all part of the same open, documented discrepancy
    // (see this file's top comment) -- printed for visibility, not asserted, since HD 292249's
    // found speed (51.7 km/s vs their reported 10.7 km/s) is itself part of that same open
    // mismatch, not a separate thing to tolerance-tune around.
    check(r.time_at_min_myr < 0, "closest approach found in the past (not the future)");

    printf("  NOTE: distance/time/speed do not match Bailer-Jones' reported values -- a known, "
           "unresolved discrepancy (not a pass/fail check here, beyond the sign). See this file's "
           "and NOTES_method_selection.md's top comments.\n\n");
}

int main(void) {
    // HIP 3757: real Gaia DR2 astrometry, source_id 2525688198820543360.
    Astrometry hip3757 = {
        .ra_deg = 12.05634594268451,
        .dec_deg = -5.1359892365165924,
        .parallax_mas = 40.20064272792413,
        .pmra_masyr = 182.16070593050424,
        .pmdec_masyr = -120.15964959812227,
        .rv_kms = 26.880687935465836,
    };
    run_case("test_pipeline_hip3757", hip3757, -0.99, 0.60, 24.7, 2.0);

    // HD 292249: real Gaia DR2 astrometry, source_id 3107000885484340224 (identity confirmed via
    // SIMBAD, which independently agrees on this source_id and its ra/dec/parallax/pm).
    Astrometry hd292249 = {
        .ra_deg = 101.09513645394681,
        .dec_deg = -1.8534936165140703,
        .parallax_mas = 24.124878291405853,
        .pmra_masyr = -5.327852712260631,
        .pmdec_masyr = -70.91032235521574,
        .rv_kms = 33.124253095329394,
    };
    run_case("test_pipeline_hd292249", hd292249, -3.765, 1.60, 10.7, 6.0);

    if (g_failures) {
        printf("%d check(s) FAILED\n", g_failures);
        return 1;
    }
    printf("Partial validation passed for both stars (encounter found in the past, as it should\n"
           "be) -- the precise encounter distance/time/speed remain an open, documented\n"
           "discrepancy.\n");
    return 0;
}
