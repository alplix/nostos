// Cross-checks potential.c's Allen & Santillan (1991) bulge (Plummer) and
// disk (Miyamoto-Nagai) acceleration formulas against galpy, component by
// component, at both in-plane and off-plane points.
//
// The halo is NOT checked here (updated 2026-10-06, when the model switched
// from a generic Hernquist+Miyamoto-Nagai+NFW placeholder to A&S's actual
// published model): A&S's halo is a generalized truncated enclosed-mass law
// with no stock galpy potential class to compare against (unlike the
// now-removed NFW placeholder, which galpy does support directly). See
// tests/test_milkyway_allen_santillan.c for how the halo is validated
// instead -- continuity and a numerical-derivative self-consistency check
// across its r=Lambda truncation boundary, since there's no independent
// reference implementation available for it here.
//
// Method: constructed galpy's PlummerPotential/MiyamotoNagaiPotential with
// the exact same mass/scale parameters as milkyway_potential_default (see
// scripts/galpy_crosscheck.py, which produced the reference numbers below),
// queried Rforce/zforce at 5 points (R=1,4,8.122,8.122,20 kpc; z=0,0,0,0.5,0
// kpc -- the last one off-plane to exercise the disk's z-dependence), and
// converted galpy's km/s/Myr output to this project's kpc/Myr^2 via
// NOSTOS_KMS_PER_KPCMYR. Both components match to ~1.7e-6 relative error at
// every point -- consistent with a G-value precision difference between
// galpy's constants and NOSTOS_G, not a formula bug.
#include <stdio.h>
#include <math.h>
#include "../src/potential.h"
#include "../src/units.h"

static int g_failures = 0;

static void check_close(double got, double ref, double tol, const char *what) {
    double rel_err = (ref != 0.0) ? fabs(got - ref) / fabs(ref) : fabs(got);
    int ok = rel_err < tol;
    if (!ok) { printf("  FAIL: %s (got=%.6e ref=%.6e relerr=%.2e)\n", what, got, ref, rel_err); g_failures++; }
    else { printf("  ok:   %s (relerr=%.2e)\n", what, rel_err); }
}

typedef struct { double R, z; } Point;
typedef struct { double aR, az; } RefAccel;

int main(void) {
    printf("test_milkyway_galpy:\n");

    Point pts[5] = {
        {1.0, 0.0}, {4.0, 0.0}, {8.122, 0.0}, {8.122, 0.5}, {20.0, 0.0}
    };

    // galpy reference accelerations (kpc/Myr^2), from scripts/galpy_crosscheck.py's output,
    // same milkyway_potential_default parameters (Allen & Santillan 1991).
    RefAccel bulge_ref[5] = {
        {-5.139443e-02, -0.0}, {-3.906292e-03, -0.0}, {-9.575422e-04, -0.0},
        {-9.521368e-04, -5.861468e-05}, {-1.583650e-04, -0.0}
    };
    RefAccel disk_ref[5] = {
        {-2.131985e-03, -0.0}, {-4.790987e-03, -0.0}, {-3.282785e-03, -0.0},
        {-3.111044e-03, -2.013401e-03}, {-8.626423e-04, -0.0}
    };

    MilkyWayParams p_bulge = milkyway_potential_default;
    p_bulge.disk_GM = 0; p_bulge.halo_GM = 0;
    MilkyWayParams p_disk = milkyway_potential_default;
    p_disk.bulge_GM = 0; p_disk.halo_GM = 0;

    Potential bulge_pot = milkyway_potential(&p_bulge);
    Potential disk_pot = milkyway_potential(&p_disk);

    const double TOL = 1e-4; // galpy itself only agrees to ~1.7e-6 on top of float rounding;
                              // 1e-4 leaves comfortable margin without being so loose it'd miss
                              // a real formula bug (which would show up as O(1) wrong, not O(1e-5)).

    for (int i = 0; i < 5; i++) {
        double R = pts[i].R, z = pts[i].z;
        Vec3 pos = {R, 0.0, z};
        char label[64];

        Vec3 a_bulge = potential_acceleration(&bulge_pot, pos);
        snprintf(label, sizeof(label), "bulge aR at R=%.3f z=%.3f", R, z);
        check_close(a_bulge.x, bulge_ref[i].aR, TOL, label);
        snprintf(label, sizeof(label), "bulge az at R=%.3f z=%.3f", R, z);
        check_close(a_bulge.z, bulge_ref[i].az, TOL, label);

        Vec3 a_disk = potential_acceleration(&disk_pot, pos);
        snprintf(label, sizeof(label), "disk  aR at R=%.3f z=%.3f", R, z);
        check_close(a_disk.x, disk_ref[i].aR, TOL, label);
        snprintf(label, sizeof(label), "disk  az at R=%.3f z=%.3f", R, z);
        check_close(a_disk.z, disk_ref[i].az, TOL, label);
    }

    if (g_failures) {
        printf("\n%d check(s) FAILED\n", g_failures);
        return 1;
    }
    printf("\nAll checks passed -- potential.c's Allen & Santillan bulge (Plummer) and disk\n"
           "(Miyamoto-Nagai) formulas are independently validated against galpy.\n");
    return 0;
}
