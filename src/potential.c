#include "potential.h"
#include "units.h"
#include <math.h>

// ---- Kepler -------------------------------------------------------------

static void kepler_accel(const Potential *pot, Vec3 pos, Vec3 *out) {
    const KeplerParams *p = (const KeplerParams *)pot->params;
    double r = v3_norm(pos);
    double f = -p->GM / (r * r * r); // a = -GM/r^3 * r_vec
    *out = v3_scale(pos, f);
}

static double kepler_energy(const Potential *pot, Vec3 pos) {
    const KeplerParams *p = (const KeplerParams *)pot->params;
    double r = v3_norm(pos);
    return -p->GM / r;
}

Potential kepler_potential(const KeplerParams *params) {
    Potential pot;
    pot.acceleration = kepler_accel;
    pot.energy = kepler_energy;
    pot.params = params;
    return pot;
}

// ---- Milky Way: Hernquist bulge + Miyamoto-Nagai disk + NFW halo --------

// See potential.h's big comment: these values are an unvalidated first-pass
// placeholder, not yet checked against a published model or galpy.
// *_GM fields are G_internal * mass_in_Msun (1e10/6.8e10/5.4e11 Msun for
// bulge/disk/halo respectively) -- computed here rather than left as raw
// masses so callers never need to know the unit-conversion constant.
const MilkyWayParams milkyway_potential_default = {
    .bulge_GM = NOSTOS_G * 1.0e10,
    .bulge_c = 0.7,
    .disk_GM = NOSTOS_G * 6.8e10,
    .disk_a = 3.0,
    .disk_b = 0.28,
    .halo_GM = NOSTOS_G * 5.4e11,
    .halo_rs = 16.0,
};

static void milkyway_accel(const Potential *pot, Vec3 pos, Vec3 *out) {
    const MilkyWayParams *p = (const MilkyWayParams *)pot->params;

    double R2 = pos.x * pos.x + pos.y * pos.y;
    double R = sqrt(R2);
    double r = v3_norm(pos);

    Vec3 a_total = {0, 0, 0};

    // Hernquist bulge: Phi = -GM/(r+c); accel = -GM/((r+c)^2 * r) * pos
    {
        double denom = (r + p->bulge_c);
        double f = -p->bulge_GM / (denom * denom * (r > 0 ? r : 1e-12));
        a_total = v3_add(a_total, v3_scale(pos, f));
    }

    // Miyamoto-Nagai disk.
    {
        double zeta = sqrt(pos.z * pos.z + p->disk_b * p->disk_b);
        double azeta = p->disk_a + zeta;
        double D2 = R2 + azeta * azeta;
        double D = sqrt(D2);
        double D3 = D2 * D;
        double aR = -p->disk_GM * R / D3; // acceleration in the R direction
        double az = -p->disk_GM * pos.z * azeta / (zeta > 0 ? zeta : 1e-12) / D3;
        if (R > 1e-12) {
            a_total.x += aR * pos.x / R;
            a_total.y += aR * pos.y / R;
        }
        a_total.z += az;
    }

    // NFW halo: enclosed-mass form, M(<r) = M_s*[ln(1+r/rs) - r/(r+rs)],
    // accel = -G*M(<r)/r^2 * r_hat.
    {
        double rs = p->halo_rs;
        double x = (r > 1e-12 ? r : 1e-12) / rs;
        double menc_factor = log(1.0 + x) - x / (1.0 + x); // dimensionless, * M_s = enclosed mass
        double r_safe = (r > 1e-12 ? r : 1e-12);
        double f = -p->halo_GM * menc_factor / (r_safe * r_safe * r_safe);
        a_total = v3_add(a_total, v3_scale(pos, f));
    }

    *out = a_total;
}

static double milkyway_energy(const Potential *pot, Vec3 pos) {
    const MilkyWayParams *p = (const MilkyWayParams *)pot->params;

    double R2 = pos.x * pos.x + pos.y * pos.y;
    double r = v3_norm(pos);

    double phi_bulge = -p->bulge_GM / (r + p->bulge_c);

    double zeta = sqrt(pos.z * pos.z + p->disk_b * p->disk_b);
    double azeta = p->disk_a + zeta;
    double D = sqrt(R2 + azeta * azeta);
    double phi_disk = -p->disk_GM / D;

    double rs = p->halo_rs;
    double x = (r > 1e-12 ? r : 1e-12) / rs;
    // Phi_halo(r) = -(G*M_s/r) * ln(1+r/rs)
    double phi_halo = -p->halo_GM / (r > 1e-12 ? r : 1e-12) * log(1.0 + x);

    return phi_bulge + phi_disk + phi_halo;
}

Potential milkyway_potential(const MilkyWayParams *params) {
    Potential pot;
    pot.acceleration = milkyway_accel;
    pot.energy = milkyway_energy;
    pot.params = params;
    return pot;
}
