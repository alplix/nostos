// Gaia (E)DR3 parallax zero-point, Lindegren et al. 2021 (A&A 649, A4), as implemented by the
// Gaia collaboration's official Python package `gaiadr3_zeropoint` (P. Ramos & A. Brown, after
// L. Lindegren's Matlab code; coefficient tables z5_200720 / z6_200720, version 0.1.0). This is a
// line-by-line C port so the C tools (src/gaia_scan.c, tools/mc_shortlist.c) can apply the
// correction without Python; tests/test_parallax_zeropoint.c pins it against reference values
// produced by the official package itself (scripts/zeropoint_crosscheck.py).
//
// Usage: parallax_corrected = parallax_catalogue - Z, where Z = the returned value (mas; Z is
// typically -0.02 to -0.04 mas, so the correction INCREASES the parallax, i.e. brings stars closer).
//
// Inputs are the Gaia DR3 columns phot_g_mean_mag, nu_eff_used_in_astrometry (5-parameter
// solutions, astrometric_params_solved = 31), pseudocolour (6-parameter solutions, = 95), ecl_lat
// (degrees) and astrometric_params_solved. The recipe is calibrated for 6 < G < 21,
// 1.1 < nu_eff < 1.9 (5p), 1.24 < pseudocolour < 1.72 (6p). Outside those ranges this function does
// what the official package does with warnings enabled: G is clamped to the table ends (6 / 21) and
// the colour basis functions are evaluated as-is (the 6p ones saturate at 1.24/1.72 by construction);
// the status is then ZPT_OUT_OF_RANGE so callers can count/flag such stars.
//
// When the recipe cannot be evaluated at all (2-parameter solution, a NaN/missing input), the
// status is ZPT_FALLBACK and the global value NOSTOS_ZPT_GLOBAL_MAS is returned instead -- the
// QSO-based median offset Lindegren et al. 2021 report for the whole catalogue (-17 uas), the
// value most papers apply when they don't use the full recipe.
#ifndef NOSTOS_PARALLAX_ZEROPOINT_H
#define NOSTOS_PARALLAX_ZEROPOINT_H

#define NOSTOS_ZPT_GLOBAL_MAS (-0.017)

typedef enum {
    ZPT_OK = 0,           // inside the recipe's calibrated range
    ZPT_OUT_OF_RANGE = 1, // evaluated, but G or colour outside the calibrated range
    ZPT_FALLBACK = 2      // recipe not applicable; NOSTOS_ZPT_GLOBAL_MAS returned
} ZptStatus;

// Which correction a tool applies (command-line option --zpt=lindegren|global|none).
typedef enum {
    ZPT_MODE_NONE = 0,      // raw catalogue parallax (the behaviour before 2026-10-08)
    ZPT_MODE_GLOBAL = 1,    // subtract NOSTOS_ZPT_GLOBAL_MAS for every star
    ZPT_MODE_LINDEGREN = 2  // full Lindegren et al. 2021 recipe (default)
} ZptMode;

// Returns Z in mas. *status may be NULL.
double gaia_dr3_parallax_zeropoint(double phot_g_mean_mag, double nu_eff_used_in_astrometry,
                                   double pseudocolour, double ecl_lat_deg,
                                   int astrometric_params_solved, ZptStatus *status);

// Z to subtract under a given mode (0 for ZPT_MODE_NONE; the global value for ZPT_MODE_GLOBAL).
double gaia_zeropoint_for_mode(ZptMode mode, double phot_g_mean_mag,
                               double nu_eff_used_in_astrometry, double pseudocolour,
                               double ecl_lat_deg, int astrometric_params_solved,
                               ZptStatus *status);

// Parses "lindegren" / "global" / "none"; returns 0 on success, -1 on an unknown name.
int gaia_zpt_mode_parse(const char *name, ZptMode *mode);
const char *gaia_zpt_mode_name(ZptMode mode);

#endif
