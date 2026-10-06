#include "nostos_io.h"
#include <stdlib.h>

bool parse_input_stream(
    FILE *f, OrbitalElements *el, double vel_cov[3][3], uint64_t *seed_base,
    double *total_time_myr, int *n_steps, double *distance_threshold_pc,
    StarRecord **stars_out, size_t *n_stars_out
) {
    double cov00, cov01, cov02, cov11, cov12, cov22;
    if (fscanf(f, "%lf %lf %lf %lf %lf", &el->e, &el->q_au, &el->i_deg, &el->node_deg, &el->peri_deg) != 5) {
        fprintf(stderr, "malformed input: can't read orbital elements\n");
        return false;
    }
    if (fscanf(f, "%lf %lf %lf %lf %lf %lf", &cov00, &cov01, &cov02, &cov11, &cov12, &cov22) != 6) {
        fprintf(stderr, "malformed input: can't read velocity covariance\n");
        return false;
    }
    vel_cov[0][0] = cov00; vel_cov[0][1] = cov01; vel_cov[0][2] = cov02;
    vel_cov[1][0] = cov01; vel_cov[1][1] = cov11; vel_cov[1][2] = cov12;
    vel_cov[2][0] = cov02; vel_cov[2][1] = cov12; vel_cov[2][2] = cov22;

    unsigned long long seed_tmp;
    if (fscanf(f, "%llu", &seed_tmp) != 1) {
        fprintf(stderr, "malformed input: can't read seed_base\n");
        return false;
    }
    *seed_base = (uint64_t)seed_tmp;

    if (fscanf(f, "%lf %d %lf", total_time_myr, n_steps, distance_threshold_pc) != 3) {
        fprintf(stderr, "malformed input: can't read search parameters\n");
        return false;
    }

    size_t capacity = 64, n = 0;
    StarRecord *stars = (StarRecord *)malloc(capacity * sizeof(StarRecord));
    if (!stars) {
        fprintf(stderr, "out of memory allocating star array\n");
        return false;
    }

    unsigned long long sid;
    StarRecord rec;
    while (fscanf(f, "%llu %lf %lf %lf %lf %lf %lf",
                  &sid, &rec.astrometry.ra_deg, &rec.astrometry.dec_deg, &rec.astrometry.parallax_mas,
                  &rec.astrometry.pmra_masyr, &rec.astrometry.pmdec_masyr, &rec.astrometry.rv_kms) == 7) {
        rec.source_id = (uint64_t)sid;
        if (n >= capacity) {
            capacity *= 2;
            StarRecord *grown = (StarRecord *)realloc(stars, capacity * sizeof(StarRecord));
            if (!grown) {
                free(stars);
                fprintf(stderr, "out of memory growing star array\n");
                return false;
            }
            stars = grown;
        }
        stars[n++] = rec;
    }

    *stars_out = stars;
    *n_stars_out = n;
    return true;
}

uint64_t combine_seed(uint64_t seed_base, uint64_t source_id) {
    uint64_t z = seed_base ^ (source_id + 0x9E3779B97F4A7C15ULL + (seed_base << 6) + (seed_base >> 2));
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
