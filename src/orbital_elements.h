// Heliocentric orbital elements (the form JPL Small-Body Database / Minor Planet Center publish
// solar-system body orbits in) -> Galactocentric phase space. A completely different starting
// point from coords.h's astrometry (observed sky position/motion) -- this is "here's the shape
// and orientation of the orbit", not "here's where it appeared on the sky and how fast it moved".
//
// Only hyperbolic orbits (eccentricity > 1) are implemented -- the physically relevant case for
// interstellar objects, which are unbound by definition. An elliptical or parabolic caller would
// need a different anomaly-solving branch this file doesn't have.
//
// Domain of validity, found the hard way while writing tests/test_orbital_elements.c: this
// two-body (Sun + ISO only) Keplerian model is only meaningful while the ISO is actually within
// the Sun's gravitational influence -- in practice a matter of YEARS around perihelion for a
// solar-system-scale hyperbolic flyby (e.g. q ~ 1 AU, e ~ 1-4), not the Myr timescales the rest
// of this project integrates over. Calling this with a multi-Myr `t_myr` produces a mean anomaly
// enormous enough to overflow sinh/cosh in the Kepler solver (hit this directly: -nan from
// exactly that). The intended use is to evaluate this near t_myr=0 (perihelion, or whatever
// epoch the published elements are referenced to) to get a STARTING Galactocentric state, which
// the long Myr-scale backward integration then continues under the Galactic potential alone
// (integrator.h / potential.c) -- the Sun's own gravity is negligible at that point and isn't
// part of that integration at all.
#ifndef NOSTOS_ORBITAL_ELEMENTS_H
#define NOSTOS_ORBITAL_ELEMENTS_H

#include <stdbool.h>
#include "coords.h"

typedef struct {
    double e;        // eccentricity, > 1 for a hyperbolic orbit
    double q_au;      // perihelion distance, AU
    double i_deg;     // inclination, degrees, in the ecliptic J2000 frame
    double node_deg;  // longitude of ascending node (capital-omega), degrees
    double peri_deg;  // argument of perihelion (lowercase-omega), degrees
} OrbitalElements;

// Returns the heliocentric-then-Galactocentric phase space state at `t_myr` Myr relative to
// perihelion passage (negative = before perihelion, the normal case when backward-integrating
// from a currently-observed position towards where the object came from). Returns false (and
// leaves *out unchanged) if el.e <= 1 (not hyperbolic) or el.q_au <= 0.
bool orbital_elements_to_galactocentric(
    OrbitalElements el, double t_myr, const SolarParams *sun, PhaseSpace *out
);

#endif
