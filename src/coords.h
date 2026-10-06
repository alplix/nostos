// Observed astrometry (ICRS RA/Dec/parallax/proper motion/radial velocity -- what Gaia and the
// ISO orbital-element solutions actually give us) <-> Galactocentric Cartesian phase space (what
// the integrator in integrator.h/potential.h needs). See the big caveat comment on
// icrs_to_galactic() below before trusting this for real science -- the rotation matrix is a
// well-known published constant, but the full 6D (position+velocity) pipeline has not yet been
// cross-checked against an independent reference (astropy/galpy) the way potential.c's Milky Way
// model is flagged as needing. The round-trip self-consistency test in
// tests/test_coords.c (forward transform then exact inverse, recovering the original input to
// float precision) is the only validation so far -- real but partial: it would catch a sign error
// or a non-orthogonal matrix, but NOT catch, say, an internally-consistent but wrong convention.
#ifndef NOSTOS_COORDS_H
#define NOSTOS_COORDS_H

#include "vec3.h"

typedef struct {
    double ra_deg;       // ICRS right ascension, degrees
    double dec_deg;      // ICRS declination, degrees
    double parallax_mas; // parallax, milliarcseconds (must be > 0)
    double pmra_masyr;   // proper motion in RA, mas/yr (already *cos(dec), Gaia's own convention)
    double pmdec_masyr;  // proper motion in Dec, mas/yr
    double rv_kms;       // radial velocity, km/s (positive = receding)
} Astrometry;

typedef struct {
    Vec3 pos; // kpc
    Vec3 vel; // kpc/Myr
} PhaseSpace;

// Sun's assumed position and velocity in the Galactocentric frame this module defines (see
// coords.c's top comment for the exact frame convention: right-handed, Galactic center at
// origin, Sun at (-r0, 0, z0), rotation in the +y direction at the Sun's location).
// Defaults are commonly-cited modern values (GRAVITY Collaboration 2018 for r0; Schonrich et al.
// 2010 for the Sun's peculiar motion) -- like the Milky Way potential, NOT yet cross-checked
// against a reference implementation, just standard literature values.
typedef struct {
    double r0_kpc;     // Sun's Galactocentric radius
    double z0_kpc;      // Sun's height above the Galactic plane
    double v_circ_kms;  // circular rotation speed at r0
    double u_pec_kms;   // Sun's peculiar motion relative to the LSR, towards the Galactic center
    double v_pec_kms;   // ... in the direction of Galactic rotation
    double w_pec_kms;   // ... towards the North Galactic Pole
} SolarParams;

extern const SolarParams solar_params_default;

// Forward: observed astrometry -> Galactocentric phase space.
void astrometry_to_galactocentric(Astrometry a, const SolarParams *sun, PhaseSpace *out);

// Inverse: Galactocentric phase space -> observed astrometry (as seen from the Sun's assumed
// position/velocity). Exists mainly so tests/test_coords.c can round-trip a made-up Astrometry
// through both functions and check it comes back unchanged -- the main self-consistency check
// this module has.
Astrometry galactocentric_to_astrometry(PhaseSpace ps, const SolarParams *sun);

#endif
