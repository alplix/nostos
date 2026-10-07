// Stage-2 Monte Carlo for a shortlist of candidate stars (NOTES_method_selection.md's two-stage
// design: "any star whose stage-1 closest approach falls under some threshold distance gets the
// full Bailer-Jones treatment"). For every star in the input file:
//   1. NOMINAL: the star's best-fit Gaia astrometry against 'Oumuamua's best-fit inbound state,
//      run through EXACTLY the same gravity_confirm() call src/gaia_scan.c uses (same window rule,
//      same 2000 steps/Myr, same 6 widening attempts) -- so this column should reproduce the
//      full-catalog scan's own numbers, a built-in consistency check.
//   2. MONTE CARLO: n_samples independent pairings of a resampled star (full Gaia 5x5 correlated
//      astrometric covariance + independent RV error, src/star_mc.c) with a resampled 'Oumuamua
//      (asymptotic velocity drawn from Bailer-Jones et al. 2018's Table 1 covariance for solution
//      7c, src/iso_solutions.h; until 2026-10-08 the JPL SBDB osculating-element covariance --
//      i.e. Bailer-Jones et al.'s "surrogates" resampling BOTH bodies). Each pairing again goes through gravity_confirm(), so the
//      definition of "the encounter" (global minimum separation within the backward window) is
//      identical to the nominal/scan one.
//
// Usage: mc_shortlist <stars.txt> <n_samples> <summary_out> [samples_out|-] [min_window_myr]
//   stars.txt: one star per line (lines starting with '#' ignored), whitespace-separated:
//     source_id ra_deg dec_deg parallax_mas pmra pmdec rv_kms
//     ra_error_mas dec_error_mas parallax_error_mas pmra_error pmdec_error rv_error_kms
//     ra_dec_corr ra_parallax_corr ra_pmra_corr ra_pmdec_corr dec_parallax_corr dec_pmra_corr
//     dec_pmdec_corr parallax_pmra_corr parallax_pmdec_corr pmra_pmdec_corr
//   (scripts/prepare_mc_shortlist.py writes this from the matched Gaia CSV rows.)
//   min_window_myr: floor on the backward search window (default 100, gaia_scan.c's own value).
//
// Summary output: one line per star (header line documents the columns). Distances pc, times
// Myr (negative = past), speeds km/s; percentiles over the valid MC samples.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "coords.h"
#include "orbital_elements.h"
#include "encounter.h"
#include "potential.h"
#include "units.h"
#include "montecarlo.h"
#include "scan_logic.h"
#include "star_mc.h"
#include "iso_solutions.h"

#define STEPS_PER_MYR 2000
#define MAX_ATTEMPTS 6

static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    return (da > db) - (da < db);
}

static double percentile(const double *sorted, int n, double p) {
    if (n <= 0) return NAN;
    double idx = p * (n - 1);
    int lo = (int)floor(idx), hi = (int)ceil(idx);
    if (lo == hi) return sorted[lo];
    double frac = idx - lo;
    return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
}

static uint64_t mix64(uint64_t z) {
    z += 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

typedef struct { double d_pc, t_myr, v_kms; int at_edge; } Enc;

static Enc run_encounter(const Potential *mw, PhaseSpace star, PhaseSpace iso, double min_window) {
    // Same window rule as src/gaia_scan.c: straight-line t_star only sizes the window, floored at
    // min_window (100 Myr there).
    double t_star = 0.0, d_star = 0.0;
    straight_line_closest_approach(iso, star, &t_star, &d_star);
    double w0 = fabs(t_star) * 2.5 + 2.0;
    if (w0 < min_window) w0 = min_window;
    double window;
    EncounterResult r = gravity_confirm(mw, star, iso, w0, STEPS_PER_MYR, MAX_ATTEMPTS, &window);
    double dt_final = window / (int)(window * STEPS_PER_MYR);
    Enc e;
    e.d_pc = r.min_distance_kpc * 1000.0;
    e.t_myr = r.time_at_min_myr;
    e.v_kms = nostos_kpcmyr_to_kms(v3_norm(v3_sub(r.star_vel_at_min, r.iso_vel_at_min)));
    e.at_edge = (fabs(fabs(r.time_at_min_myr) - window) < 2.0 * dt_final);
    return e;
}

static int parse_star(const char *line, GaiaStarUncertainty *s) {
    unsigned long long sid;
    double v[23];
    int n = sscanf(line,
        "%llu %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
        "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
        &sid, &v[0], &v[1], &v[2], &v[3], &v[4], &v[5], &v[6], &v[7], &v[8], &v[9], &v[10], &v[11],
        &v[12], &v[13], &v[14], &v[15], &v[16], &v[17], &v[18], &v[19], &v[20], &v[21]);
    if (n != 23) return 0;
    s->source_id = sid;
    s->mean = (Astrometry){ v[0], v[1], v[2], v[3], v[4], v[5] };
    s->ra_error_mas = v[6]; s->dec_error_mas = v[7]; s->parallax_error_mas = v[8];
    s->pmra_error_masyr = v[9]; s->pmdec_error_masyr = v[10]; s->rv_error_kms = v[11];
    for (int k = 0; k < GAIA_N_CORR; k++) s->corr[k] = v[12 + k];
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: %s <stars.txt> <n_samples> <summary_out> [samples_out|-] [min_window_myr]\n", argv[0]);
        return 1;
    }
    int n_samples = atoi(argv[2]);
    const char *samples_path = (argc >= 5 && strcmp(argv[4], "-") != 0) ? argv[4] : NULL;
    double min_window = (argc >= 6) ? atof(argv[5]) : 100.0;
    if (n_samples < 1) { fprintf(stderr, "n_samples must be >= 1\n"); return 1; }

    FILE *in = fopen(argv[1], "r");
    if (!in) { fprintf(stderr, "can't open %s\n", argv[1]); return 1; }
    FILE *out = fopen(argv[3], "w");
    if (!out) { fprintf(stderr, "can't open %s\n", argv[3]); return 1; }
    FILE *samp = samples_path ? fopen(samples_path, "w") : NULL;
    if (samples_path && !samp) { fprintf(stderr, "can't open %s\n", samples_path); return 1; }

    fprintf(out, "# mc_shortlist: n_samples=%d min_window_myr=%.1f steps_per_myr=%d\n",
            n_samples, min_window, STEPS_PER_MYR);
    fprintf(out, "# source_id n_valid n_rejected n_window_limited nom_d_pc nom_t_myr nom_v_kms "
                 "d_p05 d_p16 d_p50 d_p84 d_p95 t_p05 t_p16 t_p50 t_p84 t_p95 "
                 "v_p05 v_p16 v_p50 v_p84 v_p95 f_d_lt_0.5 f_d_lt_1 f_d_lt_2\n");
    if (samp) fprintf(samp, "# source_id sample d_pc t_myr v_kms window_limited\n");

    Potential mw = milkyway_potential(&milkyway_potential_default);
    // ISO start: exact inbound asymptote at the Gaia DR3 epoch from the barycentric 1600 elements,
    // same as src/gaia_scan.c (src/iso_solutions.h). ISO uncertainty (2026-10-08): the asymptotic
    // velocity is drawn from Bailer-Jones et al. 2018's own Table 1 covariance for solution 7c
    // (barycentric ICRF km/s), which -- unlike the SBDB osculating-element covariance used before --
    // includes the non-gravitational-parameter uncertainty of the whole passage (~0.01-0.03 km/s).
    const double t_eval_myr = OUMUAMUA_T_START_DR3_DAYS * NOSTOS_MYR_PER_DAY;
    PhaseSpace iso_nom;
    if (!orbital_elements_inbound_asymptote_gm(OUMUAMUA_ELEMENTS_BARY_1600, NOSTOS_GM_SSB, t_eval_myr,
                                               &solar_params_default, &iso_nom)) {
        fprintf(stderr, "bad nominal ISO elements\n"); return 1;
    }
    const double iso_mean[MC_MAX_DIM] = { 0.0, 0.0, 0.0 };
    double iso_cov[MC_MAX_DIM][MC_MAX_DIM] = {{0}};
    {
        const double *sg = OUMUAMUA_7C_VINF_SIGMA_KMS, *rc = OUMUAMUA_7C_VINF_CORR;
        const double rho[3][3] = { { 1.0, rc[0], rc[1] }, { rc[0], 1.0, rc[2] }, { rc[1], rc[2], 1.0 } };
        for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) iso_cov[i][j] = rho[i][j] * sg[i] * sg[j];
    }
    double iso_L[MC_MAX_DIM][MC_MAX_DIM];
    if (mc_cholesky(iso_cov, 3, iso_L) != 0) {
        // The published correlations are rounded to 6 digits and the matrix is nearly singular;
        // a 1e-6 relative diagonal jitter restores positive-definiteness without changing sigmas
        // meaningfully.
        for (int i = 0; i < 3; i++) iso_cov[i][i] *= 1.0 + 1e-6;
        if (mc_cholesky(iso_cov, 3, iso_L) != 0) { fprintf(stderr, "ISO covariance not PD\n"); return 1; }
    }
    static const SolarParams no_sun = { 0, 0, 0, 0, 0, 0 }; // rotation ICRS -> Galactic axes only

    double *d = malloc(n_samples * sizeof(double));
    double *t = malloc(n_samples * sizeof(double));
    double *v = malloc(n_samples * sizeof(double));

    char line[8192];
    int star_no = 0;
    while (fgets(line, sizeof line, in)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        GaiaStarUncertainty s;
        if (!parse_star(line, &s)) { fprintf(stderr, "skipping unparseable line: %.60s\n", line); continue; }
        star_no++;

        PhaseSpace star_nom;
        astrometry_to_galactocentric(s.mean, &solar_params_default, &star_nom);
        Enc nom = run_encounter(&mw, star_nom, iso_nom, min_window);

        double cov[MC_MAX_DIM][MC_MAX_DIM], L[MC_MAX_DIM][MC_MAX_DIM];
        gaia_build_cov6(&s, cov);
        if (mc_cholesky(cov, 6, L) != 0) {
            fprintf(stderr, "%llu: star covariance not positive-definite, skipping MC\n",
                    (unsigned long long)s.source_id);
            fprintf(out, "%llu 0 0 0 %.4f %.4f %.3f NOT_PD\n", (unsigned long long)s.source_id,
                    nom.d_pc, nom.t_myr, nom.v_kms);
            fflush(out);
            continue;
        }

        McRng star_rng, iso_rng;
        mc_rng_seed(&star_rng, mix64(s.source_id ^ 0x5354415253454544ULL));
        mc_rng_seed(&iso_rng, mix64(s.source_id ^ 0x49534f5345454421ULL));

        int n_valid = 0, n_rej = 0, n_edge = 0;
        for (int k = 0; k < n_samples; k++) {
            Astrometry a;
            double iso_draw[MC_MAX_DIM];
            int ok_star = gaia_sample_astrometry(&star_rng, &s, L, &a);
            mc_sample_multivariate_gaussian(&iso_rng, iso_mean, iso_L, 3, iso_draw);
            if (!ok_star) { n_rej++; continue; }
            Vec3 dv_icrs = { nostos_kms_to_kpcmyr(iso_draw[0]), nostos_kms_to_kpcmyr(iso_draw[1]),
                             nostos_kms_to_kpcmyr(iso_draw[2]) };
            PhaseSpace iso_ps = iso_nom, star_ps;
            iso_ps.vel = v3_add(iso_nom.vel, helio_equatorial_to_galactocentric_vel(dv_icrs, &no_sun));
            astrometry_to_galactocentric(a, &solar_params_default, &star_ps);
            Enc e = run_encounter(&mw, star_ps, iso_ps, min_window);
            if (e.at_edge) n_edge++;
            d[n_valid] = e.d_pc; t[n_valid] = e.t_myr; v[n_valid] = e.v_kms;
            if (samp) fprintf(samp, "%llu %d %.5f %.5f %.4f %d\n", (unsigned long long)s.source_id,
                              k, e.d_pc, e.t_myr, e.v_kms, e.at_edge);
            n_valid++;
        }

        int n05 = 0, n1 = 0, n2 = 0;
        for (int k = 0; k < n_valid; k++) {
            if (d[k] < 0.5) n05++;
            if (d[k] < 1.0) n1++;
            if (d[k] < 2.0) n2++;
        }
        qsort(d, n_valid, sizeof(double), cmp_double);
        qsort(t, n_valid, sizeof(double), cmp_double);
        qsort(v, n_valid, sizeof(double), cmp_double);
        static const double P[5] = {0.05, 0.16, 0.50, 0.84, 0.95};
        fprintf(out, "%llu %d %d %d %.4f %.4f %.3f", (unsigned long long)s.source_id,
                n_valid, n_rej, n_edge, nom.d_pc, nom.t_myr, nom.v_kms);
        for (int i = 0; i < 5; i++) fprintf(out, " %.4f", percentile(d, n_valid, P[i]));
        for (int i = 0; i < 5; i++) fprintf(out, " %.4f", percentile(t, n_valid, P[i]));
        for (int i = 0; i < 5; i++) fprintf(out, " %.3f", percentile(v, n_valid, P[i]));
        double nv = n_valid > 0 ? n_valid : 1;
        fprintf(out, " %.4f %.4f %.4f\n", n05 / nv, n1 / nv, n2 / nv);
        fflush(out);
        if (samp) fflush(samp);
        fprintf(stderr, "[%d] %llu nominal d=%.4f pc | MC median d=%.4f pc, f(<1pc)=%.3f, %d valid, %d edge\n",
                star_no, (unsigned long long)s.source_id, nom.d_pc, percentile(d, n_valid, 0.5),
                n1 / nv, n_valid, n_edge);
    }

    free(d); free(t); free(v);
    fclose(in); fclose(out);
    if (samp) fclose(samp);
    return 0;
}
