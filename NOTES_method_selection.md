# Method selection — working notes

Decisions made so far while writing the core integrator (`src/`), and what's
still genuinely open. This is a working document, not a final spec — update
it as choices get revisited.

## Unit system (settled)

Internal units: kpc (length), Myr (time), kpc/Myr (velocity), Msun (mass).
Fully self-consistent — no hidden conversion factors inside `potential.c` or
`integrator.h`. Conversion helpers to/from km/s live in `src/units.h`
(`nostos_kms_to_kpcmyr` / `nostos_kpcmyr_to_kms`), since Gaia data and
published stellar velocities are normally quoted in km/s.

## Integrator (settled, validated)

Leapfrog (kick-drift-kick), in `src/integrator.h`. Chosen over RK4 because
it's symplectic — energy/angular-momentum error oscillates rather than
secularly drifting, which matters for long backward integrations. Validated
in `tests/test_kepler.c` against a point-mass potential (no Milky-Way
parameter uncertainty involved at all): energy conserved to ~1e-7 relative
and Lz to ~1e-13 over 3 orbital periods, a circular orbit closes to 1e-5*r0
after exactly one period, and an eccentric orbit matches the analytic
Kepler-equation solution to ~1e-9 relative at several points in the orbit.
All passing as of this writing (`gcc tests/test_kepler.c src/potential.c -lm`,
or `make test`).

Leapfrog is time-symmetric, so backward integration (negative `dt`) needs no
special-casing — same code path as forward.

## Milky Way potential — RESOLVED, 2026-10-06 (Allen & Santillan 1991)

`src/potential.c`'s `milkyway_potential_default` is a Hernquist bulge +
Miyamoto-Nagai disk + NFW halo (the same three ingredients as galpy's
`MWPotential2014`), with hand-picked masses/scales that happen to produce a
rotation curve of ~200-237 km/s from 2-30 kpc (see `tests/test_milkyway.c`'s
output) — plausible-looking, but **not** taken from a published table, so
treat every number in there as unverified until one of the following happens:

- Replace the parameters with `galpy`'s actual `MWPotential2014` component
  amplitudes/scales (requires working through galpy's internal `(ro, vo)`
  normalization carefully — amplitudes are dimensionless fractions of the
  total circular velocity at `ro`, not masses directly), **or**
- Replace the parameters with Allen & Santillan (1991, Rev. Mex. Astron.
  Astrofis. 22, 255)'s published bulge/disk/halo values — a model with the
  same general shape that's specifically been used before for stellar-flyby
  /orbit-backtracking work (e.g. Garcia-Sanchez et al. 2001), which is
  arguably a better-justified choice for this project than galpy's model, but
  needs the actual paper's Table open to get the digits right — not done from
  memory.
- Or, less rigorously but useful as an interim cross-check: install `galpy`
  and have it evaluate `MWPotential2014`'s rotation curve and a few test
  orbits, diff against this C implementation's placeholder potential (won't
  match exactly since the shape differs, but should be in the same
  ballpark — a wildly different result would mean a real bug in
  `potential.c`, not just a calibration gap).

Either way: don't trust an encounter-search result produced with the current
placeholder parameters as scientifically meaningful. The code structure
(generic `Potential` interface) is designed so swapping this out later
doesn't touch `integrator.h` or anything downstream at all.

**Done, 2026-10-06**: the third bullet above, done properly rather than as a
loose ballpark check. `scripts/galpy_crosscheck.py` constructs galpy's
`HernquistPotential`/`MiyamotoNagaiPotential`/`NFWPotential` directly (not
`MWPotential2014` as a whole -- its bulge is a power-law-with-cutoff, a
genuinely different shape from Hernquist, so a whole-model match isn't
apples-to-apples) with the *exact same* mass/scale parameters as
`milkyway_potential_default`, at 5 points spanning 1-20 kpc including one
off-plane (z=0.5 kpc, to exercise the Miyamoto-Nagai/NFW vertical structure
an R-only check can't reach). `tests/test_milkyway_galpy.c` hardcodes those
reference numbers and checks `potential.c`'s real C code against them,
component by component (isolating bulge/disk/halo one at a time by zeroing
the other two components' GM, since `potential.c` only exposes the summed
acceleration). Every component matches to ~1.7e-6 relative error at every
point -- consistent with a G-constant precision difference between galpy's
internals and `NOSTOS_G`, not a formula bug.

This validated that the CODE correctly implemented Hernquist/Miyamoto-Nagai/
NFW -- it did **not** validate the placeholder amplitude/scale NUMBERS
against the real Milky Way. That second bullet is now ALSO done, properly:

**Done, 2026-10-06 (second half)**: replaced `milkyway_potential_default`
with Allen & Santillan (1991, A&SS 86, "An improved model of the galactic
mass distribution for orbit computations")'s actual published parameters --
the same model previously used for stellar-flyby/orbit-backtracking work
(Garcia-Sanchez et al. 2001), confirmed against multiple independent
secondary sources (mass unit M_gal=2.325e7 Msun; table values 606.0/3690.0/
4615.0 M_gal and 0.3873/5.3178/0.2500/12.0 kpc; halo cutoff Lambda=100 kpc,
gamma=2.02). This was not just a parameter swap: **A&S's bulge is Plummer-
shaped** (Phi=-GM/sqrt(r^2+b^2)), not Hernquist -- a real structural
difference from the earlier generic placeholder, caught by actually reading
the model's formulas rather than assuming the shape carried over. The disk
stayed Miyamoto-Nagai (A&S uses the same form the placeholder already had).
The halo became a generalized truncated enclosed-mass law with no simple
stock-galpy equivalent, so it's validated differently from bulge/disk:
`tests/test_milkyway_allen_santillan.c` checks (a) a numerical derivative of
the hand-derived energy() formula matches the analytic acceleration (the
exact kind of mistake a sign error in integrating the enclosed-mass law by
hand would produce), both well inside and well outside the cutoff, and (b)
energy and acceleration are continuous across r=Lambda itself. Both pass to
better than 1e-5 relative error. `tests/test_milkyway_galpy.c` was updated
to cross-check the new bulge (now via galpy's PlummerPotential) and disk
(new parameters) -- its old NFW-vs-halo check was removed since A&S's halo
has no galpy equivalent to compare against.

Resulting rotation curve: 202.8-232.4 km/s across 1-30 kpc (see
`tests/test_milkyway_allen_santillan.c`'s printout) -- recognizably flatter
and closer to the real Milky Way's ~220-235 km/s than the old placeholder's
curve was, as expected since these are now parameters A&S themselves fit to
reproduce that behavior, not a hand-picked guess.

## Backward-integration time window — OPEN

Not decided yet. The literature (Bailer-Jones et al. for 'Oumuamua/Borisov)
generally trusts results only within a few Myr of backward integration,
because: (a) the Galactic potential itself is only known approximately, so
orbit errors grow with integration time, and (b) for encounters further back
than that, the star's own position has often diverged too much from a
straight extrapolation for a "close encounter" to mean much anyway. A
reasonable starting range to test against is probably 1-10 Myr, but this
should be settled by literature review (what exact window did Bailer-Jones
et al. use, and did they justify a specific cutoff?) rather than guessed.

## Monte Carlo sampling scheme for orbital uncertainty — OPEN

Not decided yet. Needs: how to resample each ISO's orbital elements from its
published covariance/uncertainty (multivariate normal in orbital-element
space, or in observed-quantities space propagated through to orbital
elements?), and how many realizations per ISO (literature often uses
hundreds to a few thousand). Also needs the equivalent question for Gaia
stars' own astrometric errors (parallax/proper-motion/RV uncertainties are
non-negligible for the full-6D subset and will dominate over the ISO's own
uncertainty for some stars) — whether Gaia-side Monte Carlo sampling is even
necessary for a first pass, or whether using each star's best-fit 6D value
once is good enough until a real candidate encounter is found (then MC only
that specific star more carefully). Leaning toward the latter for the first
working pipeline, to keep the per-star-per-WU cost down, but not decided.

## Coordinate transform (written, partially validated)

`src/coords.c` converts Gaia-style astrometry (RA, Dec, parallax, proper
motion, radial velocity) to/from Galactocentric Cartesian phase space.
Frame convention (our own explicit choice, documented at the top of
coords.c): right-handed, Galactic center at the origin, Sun at (-r0, 0, z0),
+X toward the Galactic center, +Y in the direction of Galactic rotation,
+Z toward the North Galactic Pole.

Solar parameters: r0=8.122 kpc (GRAVITY Collaboration 2018), z0=0.0208 kpc,
v_circ=229 km/s, peculiar motion (U,V,W)=(11.1, 12.24, 7.25) km/s (Schönrich,
Binney & Dehnen 2010, U positive toward the Galactic center) —
`solar_params_default` in `src/coords.h`. Rough cross-check against
`potential.c`'s own placeholder Milky Way model: its rotation curve at
R=8.122 kpc comes out to 233.04 km/s (`tests/test_milkyway.c`), a few km/s
off from the 229 km/s assumed here — close enough that the two aren't
flagrantly inconsistent, but neither one is a validated number yet, so don't
read anything into the agreement.

The ICRS<->Galactic rotation matrix (`A_G` in coords.c) is a specific,
widely-reproduced published constant (ESA 1997 Hipparcos & Tycho Catalogues,
derived from the IAU-defined Galactic pole/node in ICRS) — high confidence in
the digits, but genuinely not re-derived or cross-checked against a running
reference here (no internet access in the session that wrote this). The
proper-motion-to-transverse-velocity constant (4.74...  km/s per mas/yr*pc,
via AUYR_TO_KMS in coords.c) *was* derived from first principles (AU in km,
Julian year in seconds) rather than recalled, so that one's on firmer
footing.

**Validation so far**: `tests/test_coords.c` round-trips astrometry ->
Galactocentric phase space -> astrometry across several sky positions
(including near both celestial poles and the ra=0/360 wrap) and recovers the
original input to ~1e-13 relative error or better. This proves the forward
and inverse transforms are mutually consistent (would catch a sign error or
a non-orthogonal matrix) — it does NOT prove the convention matches any
external reference. **Real TODO**: cross-check a handful of known stars
(e.g. compute Galactocentric UVW for a star with a published value in the
literature, or diff against astropy's `ICRS` -> `Galactocentric` transform
directly) before trusting this on real Gaia data.

## Encounter-search kernel (written, validated against an exact case)

`src/encounter.c`'s `find_closest_approach()` steps a star and an ISO
realization through the same potential together (no mutual gravity between
them) and tracks their minimum separation. Validated against an exact
analytic case (GM=0, both bodies moving in straight lines -- the classic
closed-form "minimize |r0 + v_rel*t|^2" answer) in `tests/test_encounter.c`,
matching to <1e-4 relative, independent of any orbit-integration-correctness
question. The reported minimum is over the discrete sampled points (every
leapfrog step), not continuously refined -- fine for now, a candidate for
later refinement (e.g. parabolic interpolation near the discovered minimum)
if dt turns out too coarse relative to typical encounter geometries once
real data is run through this.

## Orbital-elements -> Galactocentric conversion (written, two real bugs found and fixed)

`src/orbital_elements.c` turns a published heliocentric orbital-element set
(eccentricity, perihelion distance, inclination, node, argument of
perihelion -- the JPL Horizons/MPC form, hyperbolic only since e>1 is what
makes something interstellar) into a Galactocentric starting state, via:
hyperbolic Kepler equation -> true anomaly -> perifocal position/velocity ->
3-1-3 Euler rotation into the ecliptic -> fixed obliquity rotation into
ICRS-equatorial -> `coords.c`'s shared `helio_equatorial_to_galactocentric_*`
helpers (refactored out of `astrometry_to_galactocentric` specifically so
this module and coords.c can't disagree about the Sun's assumed position or
the ICRS->Galactic rotation).

**Real bug #1, found by the tests**: the hyperbolic Kepler solver's naive
`H=M` initial guess (fine for the elliptical case in `potential.c`) overflows
`sinh`/`cosh` for large `|M|` and returns `-nan`. Fixed with the standard
large-`|M|` asymptotic starting guess, `H = sign(M) * ln(2|M|/e + 1.8)`.

**Real finding #2, which is really a design clarification, not just a code
bug**: hit the above overflow specifically because the first test tried
evaluating a Borisov-like orbit several Myr from perihelion. Worked out why:
a solar-system-scale hyperbolic flyby's own characteristic timescale
(`1/mean_motion`) is on the order of **a year**, not a Myr -- so a multi-Myr
mean anomaly is enormous and physically meaningless for a pure two-body
Sun-ISO model anyway (the Sun's gravity is negligible once the object is
that far away; something has clearly left the regime this model describes
long before then). This fixed the project's actual intended pipeline
shape, now documented in `orbital_elements.h`: use this module ONLY to get
a starting state near perihelion (fractions of a year to at most a few
years out), then hand off to the Galactic-potential integrator
(`integrator.h`/`potential.c`) for the real Myr-scale backward integration
-- the two-body Keplerian piece and the Galactic-potential piece are not
meant to overlap in time at all.

**Validation so far**: `tests/test_orbital_elements.c` checks (a) heliocentric
distance at perihelion equals q exactly, (b) specific orbital energy and
angular-momentum magnitude are conserved between perihelion and ~8
years-before/~18 years-after (using a "null" Sun at the origin with zero
velocity so the function's raw output is directly the heliocentric
Keplerian state, decoupled from coords.c's Sun-offset machinery), and (c) a
weak magnitude-preservation check after the full rotation chain. (b) is a
real check of the time-domain (Kepler-equation/vis-viva) math; energy and
angular-momentum magnitude are rotation-invariant though, so none of this
actually validates the perifocal->ecliptic->equatorial rotation matrices'
correctness, only that they're at least orthogonal.

**Done, 2026-10-06**: the real cross-check this section used to flag as a TODO.
`tests/test_horizons_crosscheck.c` fetches a live JPL Horizons state vector for
1I/'Oumuamua (real published osculating elements + ICRF position/velocity at
the same epoch, queried directly from ssd.jpl.nasa.gov's API) and checks
`orbital_elements_to_galactocentric()`'s output against it end to end --
Kepler solve, the full perifocal->ecliptic->equatorial rotation chain, AND
coords.c's ICRS->Galactic matrix all at once (by rotating back with that
matrix's transpose). Both position and velocity match to ~1e-6 relative
error. This is the first validation in the project that isn't just internal
self-consistency; it specifically catches what (b) and (c) above couldn't --
a wrong but orthogonal rotation.

One real bug found while building this test (in the test itself, not the
library code): the AU<->kpc and day<->Myr unit conversions for velocity were
initially inverted (divided by the conversion constants instead of
multiplying), which silently produced a wrong-but-direction-consistent
velocity off by a constant ~3.14x factor on all three components -- an easy
trap, since a uniform-ratio error across all components looks at first
glance like it could be a legitimate physical factor (it momentarily looked
suspiciously close to pi) rather than a unit-conversion sign flip. Caught by
comparing against the real reference value, not by any internal check --
exactly the kind of bug internal-consistency tests can't catch and the whole
reason this cross-check was worth doing.
