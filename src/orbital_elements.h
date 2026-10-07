// Heliocentric orbital elements (the form JPL Small-Body Database / Minor Planet Center publish
// solar-system body orbits in) -> Galactocentric phase space. A completely different starting
// point from coords.h's astrometry (observed sky position/motion) -- this is "here's the shape
// and orientation of the orbit", not "here's where it appeared on the sky and how fast it moved".
//
// Only hyperbolic orbits (eccentricity > 1) are implemented -- the physically relevant case for
// interstellar objects, which are unbound by definition. An elliptical or parabolic caller would
// need a different anomaly-solving branch this file doesn't have.
//
// Cross-checked 2026-10-06 against a real JPL Horizons state vector for 1I/'Oumuamua -- see
// tests/test_horizons_crosscheck.c and orbital_elements.c's own comment. Position and velocity
// both match to ~1e-6 relative error.
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

// The ISO's PRE-ENCOUNTER state for a Galactic-potential-only backward integration: a point on the
// exact INBOUND ASYMPTOTE of the hyperbola, moving with the exact inbound asymptotic velocity
// (|v| = v_inf = sqrt(GM_sun/|a|), direction = limit nu -> -nu_inf). Position is
// B + v_inf_vec * t_myr, where B is the impact-parameter point (the asymptote's closest point to
// the Sun, |B| = |a| sqrt(e^2-1)) and t_myr is time relative to perihelion (pass the star data's
// epoch minus the perihelion time; the along-track ambiguity of a straight line vs the real orbit
// is only a logarithmic ~|a| ln(...) term, tens of AU at most).
//
// USE THIS, not orbital_elements_to_galactocentric() at some "large enough" pre-perihelion time,
// to start the Myr-scale backward integration. Real bug (2026-10-08, see
// NOTES_method_selection.md "Residual distance gap: the ISO speed was never asymptotic"): the
// pipeline used to evaluate the Keplerian state 9.08 yr (1I) / 3.0 yr (2I) before perihelion and
// treat it as asymptotic. The DIRECTION there is converged (~1e-4 rad, it converges like 1/r^2),
// but the SPEED is not -- vis-viva v^2 = v_inf^2 + 2GM/r converges only like 1/r: 'Oumuamua was
// 2.3% (0.61 km/s) and Borisov 3.8% (1.22 km/s) too fast, a ~0.6-1.1 pc along-track error per
// Myr of look-back -- exactly the "time and speed right, distance off 1.4-15x" residual.
// Same validity conditions as above (hyperbolic, q > 0); returns false otherwise.
bool orbital_elements_inbound_asymptote(
    OrbitalElements el, double t_myr, const SolarParams *sun, PhaseSpace *out
);

// Same, with an explicit central GM (internal units, kpc^3/Myr^2). The variant above uses GM_sun,
// right for HELIOCENTRIC elements. For the pipeline, prefer BARYCENTRIC osculating elements taken
// far out on the inbound leg (JPL Horizons, centre @0, e.g. epoch 1600 -- these already include
// planetary and non-gravitational perturbations of the whole passage) together with the total
// solar-system GM (NOSTOS_GM_SSB): that is exactly how Bailer-Jones et al. 2018/2020 obtain the
// ISO's incoming velocity ("integrate back to 3000 BCE, then extrapolate to the Keplerian
// inbound asymptote", barycentric), and it reproduces their Table 1 asymptote for 1I solution 7c
// to every printed digit (tests/test_orbital_elements.c). Heliocentric 2017-epoch elements instead
// point 0.077 deg (1I) / 0.021 deg (2I) away -- ~0.1 pc at HD 292249's -3.8 Myr.
bool orbital_elements_inbound_asymptote_gm(
    OrbitalElements el, double gm, double t_myr, const SolarParams *sun, PhaseSpace *out
);

#endif
