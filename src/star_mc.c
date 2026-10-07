#include "star_mc.h"
#include <math.h>

#define MAS_PER_DEG 3.6e6
#define DEG2RAD (3.14159265358979323846 / 180.0)

void gaia_build_cov6(const GaiaStarUncertainty *s, double cov[MC_MAX_DIM][MC_MAX_DIM]) {
    double sigma[6] = {
        s->ra_error_mas, s->dec_error_mas, s->parallax_error_mas,
        s->pmra_error_masyr, s->pmdec_error_masyr, s->rv_error_kms,
    };
    for (int i = 0; i < MC_MAX_DIM; i++)
        for (int j = 0; j < MC_MAX_DIM; j++)
            cov[i][j] = 0.0;
    for (int i = 0; i < 6; i++) cov[i][i] = sigma[i] * sigma[i];
    static const int ii[GAIA_N_CORR] = {0, 0, 0, 0, 1, 1, 1, 2, 2, 3};
    static const int jj[GAIA_N_CORR] = {1, 2, 3, 4, 2, 3, 4, 3, 4, 4};
    for (int k = 0; k < GAIA_N_CORR; k++) {
        double c = s->corr[k] * sigma[ii[k]] * sigma[jj[k]];
        cov[ii[k]][jj[k]] = c;
        cov[jj[k]][ii[k]] = c;
    }
}

Astrometry gaia_apply_offsets(const Astrometry *mean, const double offset[MC_MAX_DIM]) {
    Astrometry a = *mean;
    double cosd = cos(mean->dec_deg * DEG2RAD);
    a.ra_deg += offset[0] / cosd / MAS_PER_DEG;
    a.dec_deg += offset[1] / MAS_PER_DEG;
    a.parallax_mas += offset[2];
    a.pmra_masyr += offset[3];
    a.pmdec_masyr += offset[4];
    a.rv_kms += offset[5];
    return a;
}

bool gaia_sample_astrometry(McRng *rng, const GaiaStarUncertainty *s,
                            const double L[MC_MAX_DIM][MC_MAX_DIM], Astrometry *out) {
    static const double zero[MC_MAX_DIM] = {0};
    double off[MC_MAX_DIM];
    mc_sample_multivariate_gaussian(rng, zero, L, 6, off);
    Astrometry a = gaia_apply_offsets(&s->mean, off);
    if (a.parallax_mas <= 0.0) return false;
    *out = a;
    return true;
}
