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

    // Perifocal -> ecliptic: classic 3-1-3 Euler rotation (node Omega about z, inclination i
    // about the new x, argument of perihelion omega about the original z -- applied here as the
    // single combined matrix, standard form e.g. Curtis eq. 4.36).
    double Om = el.node_deg * DEG2RAD;
    double inc = el.i_deg * DEG2RAD;
    double w = el.peri_deg * DEG2RAD;
    double cO = cos(Om), sO = sin(Om), ci = cos(inc), si = sin(inc), cw = cos(w), sw = sin(w);
    double R[3][3] = {
        { cO*cw - sO*sw*ci, -cO*sw - sO*cw*ci,  sO*si },
        { sO*cw + cO*sw*ci, -sO*sw + cO*cw*ci, -cO*si },
        { sw*si,             cw*si,             ci    },
    };
    Vec3 pos_ecl = {
        R[0][0]*pos_peri.x + R[0][1]*pos_peri.y + R[0][2]*pos_peri.z,
        R[1][0]*pos_peri.x + R[1][1]*pos_peri.y + R[1][2]*pos_peri.z,
        R[2][0]*pos_peri.x + R[2][1]*pos_peri.y + R[2][2]*pos_peri.z,
    };
    Vec3 vel_ecl = {
        R[0][0]*vel_peri.x + R[0][1]*vel_peri.y + R[0][2]*vel_peri.z,
        R[1][0]*vel_peri.x + R[1][1]*vel_peri.y + R[1][2]*vel_peri.z,
        R[2][0]*vel_peri.x + R[2][1]*vel_peri.y + R[2][2]*vel_peri.z,
    };

    // Ecliptic -> ICRS-equatorial-aligned: rotation about x by the mean obliquity.
    double eps = OBLIQUITY_DEG * DEG2RAD;
    double ce = cos(eps), se = sin(eps);
    Vec3 pos_icrs = { pos_ecl.x, pos_ecl.y*ce - pos_ecl.z*se, pos_ecl.y*se + pos_ecl.z*ce };
    Vec3 vel_icrs = { vel_ecl.x, vel_ecl.y*ce - vel_ecl.z*se, vel_ecl.y*se + vel_ecl.z*ce };

    out->pos = helio_equatorial_to_galactocentric_pos(pos_icrs, sun);
    out->vel = helio_equatorial_to_galactocentric_vel(vel_icrs, sun);
    return true;
}
