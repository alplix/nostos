#include "encounter.h"

EncounterResult find_closest_approach(
    const Potential *pot, State *star, State *iso, double total_time, int n_steps
) {
    double dt = total_time / n_steps;

    EncounterResult result;
    Vec3 sep0 = v3_sub(star->pos, iso->pos);
    result.min_distance_kpc = v3_norm(sep0);
    result.time_at_min_myr = 0.0;
    result.star_pos_at_min = star->pos;
    result.iso_pos_at_min = iso->pos;
    result.star_vel_at_min = star->vel;
    result.iso_vel_at_min = iso->vel;

    for (int i = 0; i < n_steps; i++) {
        leapfrog_step(pot, star, dt);
        leapfrog_step(pot, iso, dt);

        Vec3 sep = v3_sub(star->pos, iso->pos);
        double d = v3_norm(sep);
        if (d < result.min_distance_kpc) {
            result.min_distance_kpc = d;
            result.time_at_min_myr = (i + 1) * dt;
            result.star_pos_at_min = star->pos;
            result.iso_pos_at_min = iso->pos;
            result.star_vel_at_min = star->vel;
            result.iso_vel_at_min = iso->vel;
        }
    }

    return result;
}
