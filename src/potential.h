// Generic Galactic-potential interface: an acceleration law (per unit mass, so
// it applies equally to a Gaia star or an interstellar object -- orbits in a
// smooth background potential don't depend on the test particle's own mass)
// plus a scalar potential energy, used for energy-conservation checks.
//
// All physical quantities are in the internal unit system from units.h:
// kpc, Myr, kpc/Myr, Msun.
#ifndef NOSTOS_POTENTIAL_H
#define NOSTOS_POTENTIAL_H

#include "vec3.h"

typedef struct Potential Potential;

struct Potential {
    void (*acceleration)(const Potential *pot, Vec3 pos, Vec3 *accel_out);
    double (*energy)(const Potential *pot, Vec3 pos);
    const void *params;
};

static inline Vec3 potential_acceleration(const Potential *pot, Vec3 pos) {
    Vec3 a;
    pot->acceleration(pot, pos, &a);
    return a;
}

static inline double potential_energy(const Potential *pot, Vec3 pos) {
    return pot->energy(pot, pos);
}

// -- Point-mass (Kepler) potential --------------------------------------
// Phi(r) = -GM/r. Exists purely to validate the integrator itself (closed
// orbits, exact conservation, exact analytic position-vs-time via Kepler's
// equation) completely independent of any Milky-Way parameter uncertainty.
typedef struct {
    double GM; // kpc^3 / Myr^2
} KeplerParams;

Potential kepler_potential(const KeplerParams *params);

// -- Milky Way potential (placeholder, needs calibration) ---------------
// Standard 3-component Galactic mass model: Hernquist bulge + Miyamoto-Nagai
// disk + NFW halo. This is the textbook SHAPE used by essentially every
// Milky Way mass model in the literature (incl. galpy's MWPotential2014,
// which uses the same three ingredients with different amplitudes/scales).
//
// The parameter VALUES below (milkyway_potential_default) are still a
// first-pass placeholder, not published parameters like Allen & Santillan
// 1991 or galpy's actual MWPotential2014 amplitudes -- picked only to give a
// roughly plausible rotation curve (~220-235 km/s near the solar radius).
// Per the project roadmap's "method selection" step, swapping in a
// published parameter set is still worth doing before trusting a real
// encounter result scientifically.
//
// What IS now validated (2026-10-06, tests/test_milkyway_galpy.c): the
// ACCELERATION FORMULAS themselves -- Hernquist, Miyamoto-Nagai, and the NFW
// enclosed-mass law, in-plane and off-plane -- match galpy (constructed with
// these same placeholder parameters) to ~1.7e-6 relative error, component by
// component. So the code correctly implements the physics these three
// standard forms describe; what's still open is only the choice of
// amplitude/scale NUMBERS, not whether the formulas are right.
// tests/test_milkyway.c's rotation-curve print remains useful for an
// empirical by-eye sanity check while tuning those numbers.
typedef struct {
    double bulge_GM;   // kpc^3/Myr^2
    double bulge_c;    // kpc, Hernquist scale radius
    double disk_GM;    // kpc^3/Myr^2
    double disk_a;     // kpc, Miyamoto-Nagai radial scale length
    double disk_b;     // kpc, Miyamoto-Nagai vertical scale height
    double halo_GM;     // kpc^3/Myr^2, NFW characteristic mass scale (G*M_s)
    double halo_rs;    // kpc, NFW scale radius
} MilkyWayParams;

extern const MilkyWayParams milkyway_potential_default;

Potential milkyway_potential(const MilkyWayParams *params);

#endif
