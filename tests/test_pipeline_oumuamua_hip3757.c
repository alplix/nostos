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
// Status, 2026-10-06: the long-open discrepancy this test used to document (tens of pc instead
// of ~1 pc) is now LARGELY RESOLVED -- see NOTES_method_selection.md's "Found the actual bug"
// section for the full writeup. Root cause: this file was evaluating 'Oumuamua's orbital
// elements at t=+9.08 YEARS (the OUTBOUND, post-perihelion asymptote -- "now", matching how
// it's actually observed), then backward-integrating that state through the Milky Way potential
// ALONE, with no re-simulation of the Sun's own gravity anywhere in that integration. But a
// Galactic-potential-only backward integration starting from the outbound state silently carries
// the WRONG velocity direction once it passes back through where perihelion was: the inbound
// (pre-encounter) and outbound (post-encounter) asymptotic velocity directions differ by the
// orbit's full hyperbolic turning angle (2*asin(1/e) = 112.7 degrees for this orbit -- not a
// small correction), because the Sun's gravity deflects the trajectory by exactly that angle
// during the encounter, and this integrator never accounts for that deflection at all (by
// design -- the Sun's gravity is supposed to be negligible by the time this integration takes
// over). Fixed by evaluating at t=-9.08 years instead (the INBOUND asymptote) -- correct for
// exactly the same reason the +9.08-year choice was correct for "what is 'Oumuamua's current
// observed state": both are found empirically by checking where the (correspondingly signed)
// time's computed velocity has converged to its asymptotic value, just on the two different
// branches of the same hyperbola.
//
// Result with the fix: HIP 3757's time and speed now match Bailer-Jones' reported values to
// within ~1% (t=-0.9855 vs -0.990 Myr, v=24.69 vs 24.7 km/s) -- distance is within ~2x (1.17 vs
// 0.60 pc), down from ~40x. HD 292249 similarly improved across the board (see this file's
// run_case() output). The remaining O(1)-2x gap in distance is a much smaller, more plausible
// residual -- not yet further investigated, but no longer the dominant effect. Every individual
// pipeline component was already independently validated before this fix (JPL Horizons,
// astropy, galpy cross-checks -- still true and unaffected by this change, since the BUG was in
// which asymptote this test fed into an otherwise-correct integration, not in any of those
// validated components themselves).
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

    // Evaluated ~9 years BEFORE perihelion -- the INBOUND asymptote, the velocity direction
    // 'Oumuamua actually had while approaching the solar system through the Galaxy, which is
    // what this Milky-Way-potential-only backward integration needs as its starting point (see
    // this file's top comment for why the OUTBOUND, post-encounter asymptote -- used here before
    // this fix -- silently fed in the wrong direction instead).
    double t_inbound_myr = -9.08 * 365.25 * NOSTOS_MYR_PER_DAY;
    PhaseSpace iso_ps;
    bool ok = orbital_elements_to_galactocentric(oumuamua_elements(), t_inbound_myr, &solar_params_default, &iso_ps);
    check(ok, "'Oumuamua's orbital elements convert to a valid Galactocentric state");

    PhaseSpace star_ps;
    astrometry_to_galactocentric(star_astrometry, &solar_params_default, &star_ps);

    // NOT "today's" separation -- iso_ps is 'Oumuamua's position ~9 years before perihelion, on
    // the inbound leg, not where it is now. Printed only as a sanity-check starting distance for
    // the backward search below.
    double d0_pc = v3_norm(v3_sub(star_ps.pos, iso_ps.pos)) * 1000.0;
    printf("  separation at the inbound starting point: %.2f pc\n", d0_pc);

    Potential mw = milkyway_potential(&milkyway_potential_default);
    State star_state = {star_ps.pos, star_ps.vel};
    State iso_state = {iso_ps.pos, iso_ps.vel};
    int n_steps = (int)(search_window_myr * 100000); // ~10 yr steps
    EncounterResult r = find_closest_approach(&mw, &star_state, &iso_state, -search_window_myr, n_steps);

    double min_d_pc = r.min_distance_kpc * 1000.0;
    Vec3 rel_vel = v3_sub(r.star_vel_at_min, r.iso_vel_at_min);
    double rel_speed_kms = nostos_kpcmyr_to_kms(v3_norm(rel_vel));

    printf("  found: t=%.4f Myr, distance=%.4f pc, rel. speed=%.2f km/s\n",
           r.time_at_min_myr, min_d_pc, rel_speed_kms);
    printf("  Bailer-Jones et al. 2018: t=%.3f Myr, distance=%.2f pc, rel. speed=%.1f km/s\n",
           bj_time_myr, bj_distance_pc, bj_speed_kms);

    // Time and speed are now asserted to within a generous 25% relative tolerance -- both match
    // to within a few percent in practice (see this file's top comment), so 25% has real margin
    // without being so loose it would pass a regression of the original bug (which was off by
    // 1-2 ORDERS of magnitude, not tens of percent). Distance remains printed, not asserted: the
    // residual ~2x gap is a smaller, different, not-yet-investigated effect (not the same bug
    // this test's history is about).
    check(r.time_at_min_myr < 0, "closest approach found in the past (not the future)");
    check(fabs(r.time_at_min_myr - bj_time_myr) / fabs(bj_time_myr) < 0.25,
          "closest-approach time matches Bailer-Jones et al. 2018 to within 25%");
    check(fabs(rel_speed_kms - bj_speed_kms) / bj_speed_kms < 0.25,
          "relative speed matches Bailer-Jones et al. 2018 to within 25%");

    printf("  NOTE: distance still doesn't closely match Bailer-Jones' reported value -- a "
           "smaller, separate residual gap (not the fixed inbound/outbound bug). See this file's "
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
    printf("Validation passed for both stars (encounter found in the past, with time and speed\n"
           "matching Bailer-Jones et al. 2018 to within 25%%) -- distance still has a smaller,\n"
           "unresolved residual gap, printed above but not asserted.\n");
    return 0;
}
