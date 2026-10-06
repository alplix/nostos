#include "scan_logic.h"
#include <math.h>

bool straight_line_closest_approach(PhaseSpace a, PhaseSpace b, double *t_star_myr, double *d_star_pc) {
    Vec3 r0 = v3_sub(b.pos, a.pos);
    Vec3 vrel = v3_sub(b.vel, a.vel);
    double vrel2 = v3_dot(vrel, vrel);
    if (vrel2 <= 0) return false;
    double t_star = -v3_dot(r0, vrel) / vrel2;
    Vec3 sep = v3_add(r0, v3_scale(vrel, t_star));
    *t_star_myr = t_star;
    *d_star_pc = v3_norm(sep) * 1000.0;
    return true;
}

EncounterResult gravity_confirm(
    const Potential *pot, PhaseSpace star, PhaseSpace iso,
    double initial_window_myr, int steps_per_myr, int max_attempts,
    double *final_window_myr
) {
    double window_myr = initial_window_myr;
    EncounterResult r;
    for (int attempt = 0; attempt < max_attempts; attempt++) {
        int n_steps = (int)(window_myr * steps_per_myr);
        if (n_steps < 1000) n_steps = 1000; // floor, same safety margin as the original gaia_scan.c logic
        State star_state = {star.pos, star.vel};
        State iso_state = {iso.pos, iso.vel};
        r = find_closest_approach(pot, &star_state, &iso_state, -window_myr, n_steps);
        double dt = window_myr / n_steps;
        bool at_edge = (fabs(fabs(r.time_at_min_myr) - window_myr) < 2.0 * dt);
        if (!at_edge || attempt == max_attempts - 1) break;
        window_myr *= 2.0;
    }
    *final_window_myr = window_myr;
    return r;
}
