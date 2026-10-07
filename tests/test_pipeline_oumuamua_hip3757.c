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
//
// 2026-10-08, the residual gap closed (NOTES_method_selection.md, "Residual distance gap"): the
// "-9.08 yr" state was asymptotic in DIRECTION but not in SPEED (2.3% too fast at ~55 AU, vis-viva
// converging only like 1/r). The pipeline state is now the exact inbound asymptote of JPL's
// barycentric 1600 elements (src/iso_solutions.h). Published values (Bailer-Jones et al. 2018,
// Table 2, solution "2 k=2", medians and 5%/95% of 2000 surrogates):
//   HIP 3757   t=-986 [-1010,-964] kyr  d=0.599 [0.527,0.669] pc  v=24.7 [24.1,25.2] km/s
//   HD 292249  t=-3765 [-4361,-3300]    d=1.605 [1.377,1.868]     v=10.7 [9.3,12.2]
// Two checks per star:
//   (1) the PIPELINE (solution 7c via JPL Horizons, raw Gaia DR2 astrometry) must land inside the
//       published 90% intervals for t, d and v. Before the fix d was 1.17 / 2.31 pc (outside).
//   (2) with the PAPER'S OWN INPUTS -- its Table 1 asymptote for solution 2 k=2 placed at the
//       barycentre at J2015.5, and DR2 parallaxes corrected by its zero-point (+0.029 mas) -- the
//       nominal must match the published MEDIANS closely: d to 3% (found 0.6023 vs 0.599 and
//       1.6431 vs 1.605), t and v to 1%. The remaining ~2% for HD 292249 is the Galactic-potential
//       model (Allen & Santillan here vs their Dauphole & Colin model; that swap alone gives
//       1.627 pc) plus MC median-vs-nominal.
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

typedef struct { double med, lo, hi; } Interval; // published median and 5%/95% bounds

// Bailer-Jones et al. 2018 Table 1, solution "2 k=2": radiant (RA, Dec, deg) and v_inf (km/s) of
// 'Oumuamua's barycentric inbound asymptote. The ISO is placed at the barycentre at t=0=J2015.5,
// as in the paper.
static PhaseSpace paper_iso_2k2(void) {
    const double d2r = 3.14159265358979323846 / 180.0, ra = 279.4752 * d2r, dec = 33.8595 * d2r;
    Vec3 radiant = { cos(dec) * cos(ra), cos(dec) * sin(ra), sin(dec) };
    PhaseSpace ps;
    ps.pos = helio_equatorial_to_galactocentric_pos((Vec3){ 0, 0, 0 }, &solar_params_default);
    ps.vel = helio_equatorial_to_galactocentric_vel(v3_scale(radiant, -nostos_kms_to_kpcmyr(26.420410)),
                                                   &solar_params_default);
    return ps;
}

static EncounterResult encounter(PhaseSpace star_ps, PhaseSpace iso_ps, double window_myr, double *v_kms) {
    Potential mw = milkyway_potential(&milkyway_potential_default);
    State star_state = {star_ps.pos, star_ps.vel};
    State iso_state = {iso_ps.pos, iso_ps.vel};
    EncounterResult r = find_closest_approach(&mw, &star_state, &iso_state, -window_myr,
                                              (int)(window_myr * 100000)); // ~10 yr steps
    *v_kms = nostos_kpcmyr_to_kms(v3_norm(v3_sub(r.star_vel_at_min, r.iso_vel_at_min)));
    return r;
}

static void run_case(
    const char *star_name, Astrometry star_astrometry,
    Interval bj_t, Interval bj_d, Interval bj_v, double search_window_myr
) {
    printf("%s:\n", star_name);

    // The pipeline's own starting state: 'Oumuamua on its exact inbound asymptote (JPL barycentric
    // 1600 elements, total solar-system GM -- src/iso_solutions.h), at the Gaia DR2 epoch J2015.5
    // that the star astrometry below refers to. (Before 2026-10-07: the OUTBOUND state; before
    // 2026-10-08: the Keplerian state 9.08 yr before perihelion, 2.3% too fast. See the top
    // comment.)
    double t_inbound_myr = (GAIA_DR2_EPOCH_JD - OUMUAMUA_TP_BARY_JD) * NOSTOS_MYR_PER_DAY;
    PhaseSpace iso_ps;
    bool ok = orbital_elements_inbound_asymptote_gm(OUMUAMUA_ELEMENTS_BARY_1600, NOSTOS_GM_SSB,
                                                    t_inbound_myr, &solar_params_default, &iso_ps);
    check(ok, "'Oumuamua's orbital elements convert to a valid Galactocentric state");

    PhaseSpace star_ps;
    astrometry_to_galactocentric(star_astrometry, &solar_params_default, &star_ps);

    // Separation at the Gaia epoch between the star and the point on 'Oumuamua's inbound asymptote
    // (within ~1 AU of the Sun). Printed only as a sanity-check starting distance for the backward
    // search below.
    double d0_pc = v3_norm(v3_sub(star_ps.pos, iso_ps.pos)) * 1000.0;
    printf("  separation at the inbound starting point: %.2f pc\n", d0_pc);

    double rel_speed_kms;
    EncounterResult r = encounter(star_ps, iso_ps, search_window_myr, &rel_speed_kms);
    double min_d_pc = r.min_distance_kpc * 1000.0;

    printf("  pipeline: t=%.4f Myr, distance=%.4f pc, rel. speed=%.3f km/s\n",
           r.time_at_min_myr, min_d_pc, rel_speed_kms);
    printf("  Bailer-Jones et al. 2018: t=%.3f [%.3f,%.3f] Myr, d=%.3f [%.3f,%.3f] pc, v=%.1f [%.1f,%.1f] km/s\n",
           bj_t.med, bj_t.lo, bj_t.hi, bj_d.med, bj_d.lo, bj_d.hi, bj_v.med, bj_v.lo, bj_v.hi);

    check(r.time_at_min_myr < 0, "closest approach found in the past (not the future)");
    check(r.time_at_min_myr > bj_t.lo && r.time_at_min_myr < bj_t.hi,
          "pipeline: time inside the published 90% interval");
    check(rel_speed_kms > bj_v.lo && rel_speed_kms < bj_v.hi,
          "pipeline: relative speed inside the published 90% interval");
    check(min_d_pc > bj_d.lo && min_d_pc < bj_d.hi,
          "pipeline: distance inside the published 90% interval");

    // Same star, the paper's own inputs (see the top comment).
    Astrometry zp = star_astrometry;
    zp.parallax_mas += 0.029; // Lindegren et al. 2018 DR2 zero-point, as applied in the paper
    PhaseSpace star_zp;
    astrometry_to_galactocentric(zp, &solar_params_default, &star_zp);
    double v_p;
    EncounterResult rp = encounter(star_zp, paper_iso_2k2(), search_window_myr, &v_p);
    double d_p = rp.min_distance_kpc * 1000.0;
    printf("  paper inputs (2 k=2 asymptote, parallax zero-point): t=%.4f Myr, d=%.4f pc, v=%.3f km/s\n",
           rp.time_at_min_myr, d_p, v_p);
    check(fabs(rp.time_at_min_myr / bj_t.med - 1.0) < 0.01, "paper inputs: time within 1% of the published median");
    check(fabs(v_p / bj_v.med - 1.0) < 0.01, "paper inputs: speed within 1% of the published median");
    check(fabs(d_p / bj_d.med - 1.0) < 0.03, "paper inputs: distance within 3% of the published median");
    printf("\n");
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
    run_case("test_pipeline_hip3757", hip3757, (Interval){ -0.986, -1.010, -0.964 },
             (Interval){ 0.599, 0.527, 0.669 }, (Interval){ 24.7, 24.1, 25.2 }, 2.0);

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
    run_case("test_pipeline_hd292249", hd292249, (Interval){ -3.765, -4.361, -3.300 },
             (Interval){ 1.605, 1.377, 1.868 }, (Interval){ 10.7, 9.3, 12.2 }, 6.0);

    if (g_failures) {
        printf("%d check(s) FAILED\n", g_failures);
        return 1;
    }
    printf("Validation passed for both stars (pipeline inside Bailer-Jones et al. 2018's 90%%\n"
           "intervals; with the paper's own inputs, medians reproduced to <= 3%%).\n");
    return 0;
}
