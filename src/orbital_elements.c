// See orbital_elements.h for the public interface. Implementation path:
//   1. Solve the hyperbolic Kepler equation for the hyperbolic anomaly H at the requested time.
//   2. Get true anomaly nu and radius r from H (via the standard, unambiguous tan(nu/2) relation
//      and the orbit equation r = p/(1+e*cos(nu)) -- deliberately NOT using the more commonly
//      quoted x=a(e-cosh H)-style perifocal formulas directly, since those have had real,
//      recurring sign-convention disagreements across textbooks; going through true anomaly and
//      the orbit equation is unambiguous and easy to sanity-check (r=q exactly at nu=0,
//      increasing monotonically away from perihelion).
//   3. Perifocal position/velocity from r, nu (standard conic formulas, valid for any
//      eccentricity).
//   4. Rotate perifocal -> ecliptic via the classic 3-1-3 Euler sequence (node, inclination,
//      argument of perihelion).
//   5. Rotate ecliptic -> ICRS-equatorial-aligned (a single rotation about the x-axis by the
//      J2000 mean obliquity).
//   6. Hand off to coords.c's shared helio_equatorial_to_galactocentric_{pos,vel}() -- the exact
//      same rotate-to-Galactic-then-shift-to-Galactocentric logic astrometry_to_galactocentric()
//      uses, so this module and coords.c can never disagree about the Sun's assumed position/
//      velocity or the ICRS->Galactic rotation.
//
// Steps 2-5 are standard orbital mechanics (e.g. Curtis, "Orbital Mechanics for Engineering
// Students"; Vallado, "Fundamentals of Astrodynamics and Applications").
//
// Cross-checked 2026-10-06 against a real JPL Horizons state vector for 1I/'Oumuamua (fetched
// live from ssd.jpl.nasa.gov's Horizons API, not a textbook example) -- see
// tests/test_horizons_crosscheck.c. Both position and velocity match to ~1e-6 relative error,
// which validates not just the Kepler-solve/radial part (already covered by
// tests/test_orbital_elements.c's energy/angular-momentum checks) but the actual DIRECTION of
// the full perifocal->ecliptic->equatorial->Galactic rotation chain -- the thing the older
// energy/angular-momentum/magnitude-only tests explicitly could not catch, since any orthogonal
// rotation preserves those quantities regardless of whether the rotation itself is right.
#include "orbital_elements.h"
#include "units.h"
#include <math.h>

#define DEG2RAD (3.14159265358979323846 / 180.0)

// J2000 mean obliquity of the ecliptic (IAU-standard value, degrees).
#define OBLIQUITY_DEG 23.43929111

static double solve_hyperbolic_kepler(double M, double e) {
    // H=M as the initial guess (as in potential.c's analogous elliptical solver) diverges for
    // large |M| -- e*sinh(H) grows like (e/2)*exp(|H|), so a naive Newton step from a starting
    // point far from the true root overflows long before converging. Real bug, hit while testing
    // this function: calling it with a multi-Myr time offset (nonsensical for a solar-system-
    // scale hyperbolic flyby anyway -- see this module's own domain-of-validity note in
    // orbital_elements.h) produced -nan from exactly this overflow. Fixed with the standard
    // robust starting guess (e.g. Vallado): for M=0, H=0 exactly; otherwise H = sign(M) *
    // ln(2|M|/e + 1.8), the asymptotic inverse of e*sinh(H)-H for large |H|, which stays close
    // to the true root across the whole range this project actually needs (fractions of a year
    // to a few years around perihelion -- not millions of years, see this file's domain note).
    double H;
    if (M == 0.0) {
        H = 0.0;
    } else {
        double sign = (M > 0) ? 1.0 : -1.0;
        H = sign * log(2.0 * fabs(M) / e + 1.8);
    }
    for (int iter = 0; iter < 200; iter++) {
        double f = e * sinh(H) - H - M;
        double fp = e * cosh(H) - 1.0;
        double dH = f / fp;
        H -= dH;
        if (fabs(dH) < 1e-14) break;
    }
    return H;
}

// Perifocal (x towards perihelion, z = orbit normal) -> heliocentric ICRS-equatorial-aligned axes:
// the classic 3-1-3 Euler rotation into the J2000 ecliptic (node Omega about z, inclination i about
// the new x, argument of perihelion omega about the original z -- applied as the single combined
// matrix, standard form e.g. Curtis eq. 4.36), then a rotation about x by the J2000 mean obliquity.
// Shared by orbital_elements_to_galactocentric() and orbital_elements_inbound_asymptote() so the
// two can never disagree about orientation (the chain is cross-checked against JPL Horizons in
// tests/test_horizons_crosscheck.c).
static Vec3 perifocal_to_icrs(OrbitalElements el, Vec3 v) {
    double Om = el.node_deg * DEG2RAD;
    double inc = el.i_deg * DEG2RAD;
    double w = el.peri_deg * DEG2RAD;
    double cO = cos(Om), sO = sin(Om), ci = cos(inc), si = sin(inc), cw = cos(w), sw = sin(w);
    double R[3][3] = {
        { cO*cw - sO*sw*ci, -cO*sw - sO*cw*ci,  sO*si },
        { sO*cw + cO*sw*ci, -sO*sw + cO*cw*ci, -cO*si },
        { sw*si,             cw*si,             ci    },
    };
    Vec3 ecl = {
        R[0][0]*v.x + R[0][1]*v.y + R[0][2]*v.z,
        R[1][0]*v.x + R[1][1]*v.y + R[1][2]*v.z,
        R[2][0]*v.x + R[2][1]*v.y + R[2][2]*v.z,
    };
    double eps = OBLIQUITY_DEG * DEG2RAD;
    double ce = cos(eps), se = sin(eps);
    return (Vec3){ ecl.x, ecl.y*ce - ecl.z*se, ecl.y*se + ecl.z*ce };
}

bool orbital_elements_to_galactocentric(
    OrbitalElements el, double t_myr, const SolarParams *sun, PhaseSpace *out
) {
    if (el.e <= 1.0 || el.q_au <= 0.0) return false;

    double q_kpc = el.q_au * NOSTOS_KPC_PER_AU;
    double a_kpc = q_kpc / (el.e - 1.0); // magnitude of the semi-major axis (positive)
    double p_kpc = q_kpc * (1.0 + el.e); // semi-latus rectum -- sign-unambiguous, same formula
                                          // works for ellipse or hyperbola

    double n = sqrt(NOSTOS_GM_SUN / (a_kpc * a_kpc * a_kpc)); // mean motion, 1/Myr
    double M = n * t_myr; // mean anomaly; unbounded (not mod 2*pi) for a hyperbolic orbit
    double H = solve_hyperbolic_kepler(M, el.e);

    double nu = 2.0 * atan2(sqrt(el.e + 1.0) * sinh(H / 2.0), sqrt(el.e - 1.0) * cosh(H / 2.0));
    double r_kpc = p_kpc / (1.0 + el.e * cos(nu));

    // Perifocal frame (x towards perihelion, z = orbit normal). Standard conic position/velocity
    // in terms of true anomaly -- valid for any eccentricity, no sign ambiguity.
    double h = sqrt(NOSTOS_GM_SUN * p_kpc); // specific angular momentum, kpc^2/Myr
    Vec3 pos_peri = { r_kpc * cos(nu), r_kpc * sin(nu), 0.0 };
    Vec3 vel_peri = {
        -(NOSTOS_GM_SUN / h) * sin(nu),
        (NOSTOS_GM_SUN / h) * (el.e + cos(nu)),
        0.0,
    };

    Vec3 pos_icrs = perifocal_to_icrs(el, pos_peri);
    Vec3 vel_icrs = perifocal_to_icrs(el, vel_peri);

    out->pos = helio_equatorial_to_galactocentric_pos(pos_icrs, sun);
    out->vel = helio_equatorial_to_galactocentric_vel(vel_icrs, sun);
    return true;
}

bool orbital_elements_inbound_asymptote(
    OrbitalElements el, double t_myr, const SolarParams *sun, PhaseSpace *out
) {
    return orbital_elements_inbound_asymptote_gm(el, NOSTOS_GM_SUN, t_myr, sun, out);
}

bool orbital_elements_inbound_asymptote_gm(
    OrbitalElements el, double gm, double t_myr, const SolarParams *sun, PhaseSpace *out
) {
    if (el.e <= 1.0 || el.q_au <= 0.0 || gm <= 0.0) return false;

    double e = el.e;
    double q_kpc = el.q_au * NOSTOS_KPC_PER_AU;
    double a_kpc = q_kpc / (e - 1.0);        // |a|
    double s = sqrt(e * e - 1.0);
    double v_inf = sqrt(gm / a_kpc); // hyperbolic excess speed (vis-viva, r -> infinity)

    // Perifocal frame, focus (Sun, or the barycentre for barycentric elements) at the origin, perihelion on +x, motion counter-clockwise
    // (+y at perihelion). The inbound branch comes in from true anomaly nu -> -nu_inf, with
    // cos(nu_inf) = -1/e; the limiting velocity there, from the conic formula
    // v = (GM/h) * (-sin nu, e + cos nu), is v_inf * (1/e, sqrt(e^2-1)/e) -- a unit vector, since
    // (1 + e^2 - 1)/e^2 = 1. (The outbound asymptote is v_inf * (-1/e, sqrt(e^2-1)/e); the angle
    // between the two is the turning angle 2*asin(1/e).)
    Vec3 u_peri = { 1.0 / e, s / e, 0.0 };
    // The asymptote passes through the hyperbola's centre, (|a| e, 0) from the focus. Its point of
    // closest approach to the Sun -- the impact-parameter vector, |B| = |a| sqrt(e^2-1) = b -- is
    // C - (C.u) u with C.u = |a|: B = |a| ((e^2-1)/e, -sqrt(e^2-1)/e).
    Vec3 b_peri = { a_kpc * (e * e - 1.0) / e, -a_kpc * s / e, 0.0 };
    Vec3 pos_peri = v3_add(b_peri, v3_scale(u_peri, v_inf * t_myr));
    Vec3 vel_peri = v3_scale(u_peri, v_inf);

    out->pos = helio_equatorial_to_galactocentric_pos(perifocal_to_icrs(el, pos_peri), sun);
    out->vel = helio_equatorial_to_galactocentric_vel(perifocal_to_icrs(el, vel_peri), sun);
    return true;
}
