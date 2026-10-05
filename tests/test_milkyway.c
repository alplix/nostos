// Not a pass/fail correctness test (there's nothing independent to check the
// Milky Way potential's parameters against yet -- see potential.h's caveat).
// This just prints the rotation curve v_c(R) = sqrt(R * |a_R(R)|) the current
// milkyway_potential_default parameters produce, in km/s, so they can be
// eyeballed/tuned against the real Milky Way's well-measured ~220-235 km/s
// near the solar radius (~8.1-8.2 kpc) before anything is built on top of it.
#include <stdio.h>
#include <math.h>
#include "../src/potential.h"
#include "../src/units.h"

int main(void) {
    Potential pot = milkyway_potential(&milkyway_potential_default);

    printf("Rotation curve (z=0 plane):\n");
    printf("  R (kpc)   v_c (km/s)\n");
    double radii[] = {1, 2, 4, 6, 8, 8.122, 10, 12, 15, 20, 25, 30};
    for (size_t i = 0; i < sizeof(radii) / sizeof(radii[0]); i++) {
        double R = radii[i];
        Vec3 pos = {R, 0, 0};
        Vec3 a = potential_acceleration(&pot, pos);
        double a_R = fabs(a.x); // a.x is the only nonzero component at y=z=0
        double v_c_internal = sqrt(R * a_R);
        double v_c_kms = nostos_kpcmyr_to_kms(v_c_internal);
        printf("  %7.3f   %8.2f\n", R, v_c_kms);
    }

    printf("\nFor reference: the real Milky Way's rotation curve is close to\n");
    printf("flat at roughly 220-235 km/s from a few kpc out to beyond the solar\n");
    printf("radius (~8.1-8.2 kpc). If the numbers above don't look like that,\n");
    printf("milkyway_potential_default's masses/scales in potential.c need\n");
    printf("adjusting (or replacing with a published model's parameters).\n");
    return 0;
}
