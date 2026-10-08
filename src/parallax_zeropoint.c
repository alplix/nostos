// C port of gaiadr3_zeropoint 0.1.0 (zero_point/zpt.py, _calc_zpt + get_zpt); see the header.
#include "parallax_zeropoint.h"
#include <math.h>
#include <string.h>

#define ZPT_NG 13

// G-magnitude knots shared by both tables.
static const double ZPT_G[ZPT_NG] = {
    6.0, 10.8, 11.2, 11.8, 12.2, 12.9, 13.1, 15.9, 16.1, 17.5, 19.0, 20.0, 21.0
};

// z5_200720.txt: columns are (j,k) = (0,0) (0,1) (0,2) (1,0) (1,1) (2,0) (3,0) (4,0); uas.
#define Z5_M 8
static const int Z5_J[Z5_M] = { 0, 0, 0, 1, 1, 2, 3, 4 };
static const int Z5_K[Z5_M] = { 0, 1, 2, 0, 1, 0, 0, 0 };
static const double Z5_Q[ZPT_NG][Z5_M] = {
    { -26.98,  -9.62,  27.40,  -25.1,  -0.0, -1257,     0,     0 },
    { -27.23,  -3.07,  23.04,   35.3,  15.7, -1257,     0,     0 },
    { -30.33,  -9.23,   9.08,  -88.4, -11.8, -1257,     0,     0 },
    { -33.54, -10.08,  13.28, -126.7,  11.6, -1257,     0,     0 },
    { -13.65,  -0.07,   9.35, -111.4,  40.6, -1257,     0,     0 },
    { -19.53,  -1.64,  15.86,  -66.8,  20.6, -1257,     0,     0 },
    { -37.99,   2.63,  16.14,   -5.7,  14.0, -1257, 107.9, 104.3 },
    { -38.33,   5.61,  15.42,      0,  18.7, -1189, 243.8, 155.2 },
    { -31.05,   2.83,   8.59,      0,  15.5, -1404, 105.5, 170.7 },
    { -29.18,  -0.09,   2.41,      0,  24.5, -1165, 189.7, 325.0 },
    { -18.40,   5.98,  -6.46,      0,   5.5,     0,     0, 276.6 },
    { -12.65,  -4.57,  -7.46,      0,  97.9,     0,     0,     0 },
    { -18.22, -15.24, -18.54,      0, 128.2,     0,     0,     0 },
};

// z6_200720.txt: (j,k) = (0,0) (0,1) (0,2) (1,0) (1,1) (1,2) (2,0); uas.
#define Z6_M 7
static const int Z6_J[Z6_M] = { 0, 0, 0, 1, 1, 1, 2 };
static const int Z6_K[Z6_M] = { 0, 1, 2, 0, 1, 2, 0 };
static const double Z6_Q[ZPT_NG][Z6_M] = {
    { -27.85,  -7.78,  27.47,  -32.1,  14.4,    9.5,   -67 },
    { -28.91,  -3.57,  22.92,    7.7,  12.6,    1.6,  -572 },
    { -26.72,  -8.74,   9.36,  -30.3,   5.6,   17.2, -1104 },
    { -29.04,  -9.69,  13.63,  -49.4,  36.3,   17.7, -1129 },
    { -12.39,  -2.16,  10.23,  -92.6,  19.8,   27.6,  -365 },
    { -18.99,  -1.93,  15.90,  -57.2,  -8.0,   19.9,  -554 },
    { -38.29,   2.59,  16.20,  -10.5,   1.4,    0.4,  -960 },
    { -36.83,   4.20,  15.76,   22.3,  11.1,   10.0, -1367 },
    { -28.37,   1.99,   9.28,   50.4,  17.2,   13.7, -1351 },
    { -24.68,  -1.37,   3.52,   86.8,  19.8,   21.3, -1380 },
    { -15.32,   4.01,  -6.03,   29.2,  14.1,    0.4,  -563 },
    { -13.73, -10.92,  -8.30,  -74.4, 196.4,  -42.0,   536 },
    { -29.53, -20.34, -18.74,  -39.5, 326.8, -262.3,  1598 },
};

static double dmin(double a, double b) { return a < b ? a : b; }
static double dmax(double a, double b) { return a > b ? a : b; }

double gaia_dr3_parallax_zeropoint(double g, double nu_eff, double pseudocolour, double ecl_lat_deg,
                                   int params_solved, ZptStatus *status) {
    ZptStatus st = ZPT_OK;
    double colour;
    int m;
    const int *J, *K;
    const double (*Q)[8];
    const double (*Q6)[7];
    int is5 = (params_solved == 31), is6 = (params_solved == 95);
    if (is5) colour = nu_eff;
    else if (is6) colour = pseudocolour;
    else { if (status) *status = ZPT_FALLBACK; return NOSTOS_ZPT_GLOBAL_MAS; }
    if (isnan(g) || isnan(colour) || isnan(ecl_lat_deg)) {
        if (status) *status = ZPT_FALLBACK;
        return NOSTOS_ZPT_GLOBAL_MAS;
    }
    if (g >= 21.0 || g <= 6.0) st = ZPT_OUT_OF_RANGE;
    if (is5 && (colour >= 1.9 || colour <= 1.1)) st = ZPT_OUT_OF_RANGE;
    if (is6 && (colour >= 1.72 || colour <= 1.24)) st = ZPT_OUT_OF_RANGE;

    // Colour basis functions c_j and ecliptic-latitude basis functions b_k (zpt.py, _calc_zpt).
    double c[5];
    c[0] = 1.0;
    c[1] = dmax(-0.24, dmin(0.24, colour - 1.48));
    c[2] = pow(dmin(0.24, dmax(0.0, 1.48 - colour)), 3);
    c[3] = dmin(0.0, colour - 1.24);
    c[4] = dmax(0.0, colour - 1.72);
    double sb = sin(ecl_lat_deg * (3.14159265358979323846 / 180.0));
    double b[3] = { 1.0, sb, sb * sb - 1.0 / 3.0 };

    // G bin: np.digitize(g, G, right=False) - 1, clamped to [0, NG-2]; then the interpolation
    // weight h clamped to [0, 1] (so G outside 6-21 takes the end values).
    int ig = 0;
    while (ig < ZPT_NG && g >= ZPT_G[ig]) ig++;
    ig -= 1;
    if (ig < 0) ig = 0;
    if (ig > ZPT_NG - 2) ig = ZPT_NG - 2;
    double h = dmax(0.0, dmin(1.0, (g - ZPT_G[ig]) / (ZPT_G[ig + 1] - ZPT_G[ig])));

    double z_uas = 0.0;
    if (is5) {
        m = Z5_M; J = Z5_J; K = Z5_K; Q = Z5_Q;
        for (int i = 0; i < m; i++)
            z_uas += ((1.0 - h) * Q[ig][i] + h * Q[ig + 1][i]) * c[J[i]] * b[K[i]];
    } else {
        m = Z6_M; J = Z6_J; K = Z6_K; Q6 = Z6_Q;
        for (int i = 0; i < m; i++)
            z_uas += ((1.0 - h) * Q6[ig][i] + h * Q6[ig + 1][i]) * c[J[i]] * b[K[i]];
    }
    if (status) *status = st;
    return z_uas * 0.001;
}

double gaia_zeropoint_for_mode(ZptMode mode, double g, double nu_eff, double pseudocolour,
                               double ecl_lat_deg, int params_solved, ZptStatus *status) {
    switch (mode) {
    case ZPT_MODE_NONE:
        if (status) *status = ZPT_OK;
        return 0.0;
    case ZPT_MODE_GLOBAL:
        if (status) *status = ZPT_OK;
        return NOSTOS_ZPT_GLOBAL_MAS;
    case ZPT_MODE_LINDEGREN:
    default:
        return gaia_dr3_parallax_zeropoint(g, nu_eff, pseudocolour, ecl_lat_deg, params_solved, status);
    }
}

int gaia_zpt_mode_parse(const char *name, ZptMode *mode) {
    if (strcmp(name, "lindegren") == 0) { *mode = ZPT_MODE_LINDEGREN; return 0; }
    if (strcmp(name, "global") == 0) { *mode = ZPT_MODE_GLOBAL; return 0; }
    if (strcmp(name, "none") == 0) { *mode = ZPT_MODE_NONE; return 0; }
    return -1;
}

const char *gaia_zpt_mode_name(ZptMode mode) {
    switch (mode) {
    case ZPT_MODE_NONE: return "none";
    case ZPT_MODE_GLOBAL: return "global";
    default: return "lindegren";
    }
}
