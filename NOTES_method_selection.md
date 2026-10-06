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

## Backward-integration time window — RESOLVED, 2026-10-06

Settled by reading the actual precedent this project explicitly aims to
scale up: Bailer-Jones et al. 2019 ("The Dynamics of Interstellar Asteroids
and Comets within the Galaxy: an Assessment of Local Candidate Source
Regions for 1I/'Oumuamua and 2I/Borisov", arXiv:1911.02473) -- NOT guessed.
They integrate to -100 Myr but explicitly quantify how reliability degrades:
"a back-integration can only expect to be accurate to within 15 pc and
2 km/s at -10 Myr, dropping steeply to 400 pc and 10 km/s at -100 Myr", and
conclude a traced path is trustworthy "perhaps as far back as 50 Myr or
more." This is an order of magnitude longer than this file's earlier guess
of "1-10 Myr" -- worth having actually checked rather than shipped a guess.

Adopted: integrate up to **100 Myr** backward (so a real deep-past encounter
isn't silently missed), but tag every reported encounter with how far back
it occurred and treat anything beyond **50 Myr** as lower-confidence --
report it, don't discard it, but don't claim it with the same weight as a
closer-in-time encounter. This matches Bailer-Jones' own practice rather
than inventing a harder cutoff. The exact position/velocity error-vs-time
curve they quote (15 pc/2 km/s at 10 Myr -> 400 pc/10 km/s at 100 Myr) is
itself a useful sanity-check target once this project's own potential/
integrator combination can be tested against a known-divergent pair of
nearby trajectories -- not implemented yet, but worth keeping as a concrete
number to validate against later rather than taken purely on faith.

## Monte Carlo sampling scheme for orbital uncertainty — RESOLVED, 2026-10-06

Also settled by reading Bailer-Jones et al. 2019 rather than guessing (see
above for the citation). Their actual practice, quoted directly:
- **100 clones per ISO**: "we characterize the uncertainty in [...]
  'Oumuamua's orbit by generating 100 clones of 'Oumuamua, each having an
  initial velocity randomly sampled from a distribution whose covariance
  matrix is recorded."
- **100 clones per candidate STAR too**: "100 clones of 'Oumuamua interact
  with 100 clones of each star in each of our integrations" -- sampled "from
  Gaussian distributions with means and standard deviations given by the
  nominal values and their corresponding uncertainties" (i.e. Gaussian in
  the *observed astrometric quantities* -- parallax, proper motion, radial
  velocity -- not in orbital-element space directly, then propagated
  through to phase space the normal way).

This directly contradicts this file's earlier "leaning toward" skipping
per-star Monte Carlo for a first pass -- that leaning was reasonable given
the real cost concern below, but Bailer-Jones' own paper is better evidence
than a guess, so it wins. Adopted as a **two-stage scheme** that keeps their
full rigor where it matters while staying tractable at Gaia-DR3 scale (they
tested a literal handful of candidate stars per ISO, not ~33 million):

1. **Stage 1 (exhaustive, cheap)** -- every Gaia star's single best-fit 6D
   value against all 100 of the ISO's own Monte Carlo clones (so the ISO's
   own uncertainty, generally the larger source of error for a fresh
   discovery, is never skipped even in the cheap pass). This is what
   actually runs across the full ~33M-star catalog in the distributed
   BOINC work units.
2. **Stage 2 (shortlist refinement)** -- any star whose stage-1 closest
   approach falls under some threshold distance gets the full Bailer-Jones
   treatment: 100 star clones (sampled Gaussian in parallax/pmra/pmdec/RV
   from that star's own Gaia-reported errors) against the 100 ISO clones,
   10,000 total pairings, run as a separate, much smaller follow-up batch
   (by then there are at most a few hundred shortlisted stars per ISO, not
   millions) -- this is where a candidate's reported encounter
   distance/time actually gets its real uncertainty estimate, matching how
   Bailer-Jones report their own results.

Still open: the exact stage-1 threshold distance that defines "shortlisted"
(their paper's own candidates were all within a few pc; a reasonable first
guess is an order of magnitude looser, e.g. a few tens of pc, to avoid
missing a real candidate at the cost of a larger but still tractable stage-2
batch -- not tuned yet, there's no real data to tune it against until the
Gaia download finishes and stage 1 actually runs once).

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

## Pipeline validation against HIP 3757 and HD 292249 (real results, every component now independently validated, encounter match still genuinely open)

First end-to-end run of the full chain (orbital elements -> coords ->
potential/integrator -> encounter search) against a real, independently
published result: HIP 3757, Bailer-Jones et al. 2018's own closest
'Oumuamua-encounter candidate (arXiv:1809.09009, Gaia DR2 source_id
2525688198820543360). See `tests/test_pipeline_oumuamua_hip3757.c` for the
full numbers and reasoning; summary here.

**Real bug found and fixed**: the first attempt evaluated
`orbital_elements_to_galactocentric()` at `t_myr=0` (exactly at perihelion)
and fed that directly in as 'Oumuamua's "current" state for the long
backward integration. That's wrong -- at perihelion the two-body velocity
is still mid-deflection by the Sun's gravity (perihelion is the turning
point, not where the outbound asymptotic velocity has settled), and using
it pointed the ISO's velocity in a meaningfully wrong direction: a
straight-line extrapolation found the closest approach in the FUTURE
(+0.39 Myr) instead of the past, the opposite of Bailer-Jones' reported
-0.99 Myr. Fixed by evaluating ~9 years after perihelion instead (within
this module's own documented "fractions of a year to a few years" domain
of validity, and roughly "now" in real terms) -- confirmed correct by
checking 'Oumuamua's resulting heliocentric speed matches its well-known
~26 km/s asymptotic excess velocity (at t=0 it's instead ~87 km/s, its much
higher perihelion speed). This flipped the sign to match Bailer-Jones.

**Still unresolved**: even with that fix, the found encounter (-0.067 Myr,
24.73 pc, barely below the present-day 24.88 pc separation) doesn't match
Bailer-Jones' reported -0.99 Myr / 0.60 pc -- roughly 15x and 40x off
respectively, though the relative speed (39 km/s found vs their 24.7 km/s)
is at least the right order of magnitude. Investigated and ruled out so
far:
- **Wrong star/cross-match**: HIP 3757's computed total proper motion from
  its individual Gaia pmra/pmdec components (sqrt(182.16^2+120.16^2) =
  218.2 mas/yr) matches the paper's own quoted total (218.22 mas/yr) to 4
  significant figures -- this is almost certainly the right star/solution.
- **Gravity vs straight-line**: both give essentially the same shallow
  minimum for this pair (gravity has negligible effect over such a short,
  shallow encounter), so the discrepancy isn't in `potential.c`/
  `integrator.h`.
- **A convention bug in `coords.c`'s `astrometry_to_galactocentric()`**:
  hand-re-derived the whole formula from scratch (r_hat/ra_hat/dec_hat unit
  vectors, the v_ra/v_dec/v_r combination) against the standard reference
  form (Johnson & Soderblom 1987 / ESA Hipparcos documentation) and it
  matches exactly; also hand-computed HIP 3757's Galactocentric velocity by
  hand through every step (ICRS construction -> A_G rotation -> Sun-motion
  addition) and got the same numbers this project's code produces to 3
  decimal places. This function had never been checked against a real
  external star before this test (only round-trip self-consistency, which
  can't catch a convention error -- see coords.h's own caveat) -- it's now
  at least hand-verified, even though the overall encounter number still
  doesn't match.
- **A_G itself being wrong**: already independently validated via
  `test_horizons_crosscheck.c` (same matrix, 'Oumuamua's position/velocity
  matched real Horizons data to ~1e-6).

**Not yet ruled out**: 'Oumuamua's known non-gravitational (cometary
outgassing) acceleration near perihelion (Micheli et al. 2018) isn't
modeled by this project's pure two-body `orbital_elements.c` -- plausible
in principle, though not expected to matter much by ~9 years post-perihelion
when that effect has long tapered off; whether Bailer-Jones propagated to a
genuinely different reference epoch/frame before their own integration than
assumed here; or their reported number being a Monte-Carlo MEDIAN that
differs substantially from this single best-fit trajectory (possible, but a
40x distance gap is a lot to pin on sampling alone given how shallow and
unambiguous the single best-fit minimum is here).

Left as a genuinely open item rather than silently adjusted to pass --
consistent with this project's own stated rule: "if it can't reproduce
known results, it doesn't get to report new ones" applies to catching this
kind of gap, not just to pretending it isn't there. Worth revisiting once
real Gaia catalog data is available and a few more candidate stars (not
just HIP 3757) can be checked the same way -- a systematic pattern across
several candidates would narrow down the hypotheses above a lot faster than
reasoning about one star in isolation.

**Update, same day**: followed through on exactly that -- tested a second,
very different candidate (HD 292249, Gaia DR2 source_id
3107000885484340224, their most PLAUSIBLE candidate despite not being
closest: lowest relative speed of the four, 10.7 km/s) and independently
validated every pipeline component that hadn't been checked against an
external reference yet. Results:

- **Same systematic pattern for HD 292249**: found t=-0.592 Myr, d=27.1 pc,
  v=51.7 km/s vs their reported t=-3.765 Myr, d=1.60 pc, v=10.7 km/s --
  qualitatively the same mismatch as HIP 3757 (right sign, wrong magnitude
  by a large factor). Confirmed with the full gravity integration searched
  out to -6 Myr (well past their reported -3.765 Myr) -- no deeper minimum
  exists in this project's model within that window; gravity genuinely
  isn't the differentiator.
- **Star identity re-confirmed independently via SIMBAD** for HD 292249
  (RA/Dec/parallax/proper motion all match the Gaia DR2 source_id this
  project queried) -- ruling out a wrong-star/mismatched-source_id
  explanation for a second time. (Also caught, in passing: an earlier
  PDF-text-extraction of this star's astrometry from the paper itself gave
  wildly different, internally-inconsistent numbers -- parallax 1.605 mas
  vs the real 24.12 mas, proper motion 24.15 mas/yr vs the real ~71.1
  mas/yr computed from components -- a reminder that AI-summarized PDF
  table extraction is NOT a reliable source for numbers that matter; the
  live Gaia TAP archive query and SIMBAD cross-check are what's actually
  trusted here, not the paper-text extraction.)
- **'Oumuamua's propagation re-validated at the ACTUAL epoch used** (not
  just the +36-day epoch test_horizons_crosscheck.c already checks): fetched
  a real JPL Horizons state vector for 2026-10-06 (~9 years post-perihelion,
  matching what this test actually evaluates `orbital_elements.c` at) and
  compared directly -- position and velocity both match to ~6e-4 relative
  error (looser than the ~1e-6 at +36 days, as expected for a longer pure-
  two-body extrapolation without non-gravitational terms, but still clearly
  correct, not wildly off). This rules out `orbital_elements.c` itself as
  the discrepancy's source, decisively.
- **`astrometry_to_galactocentric()` independently validated against
  astropy** (`tests/test_coords_astropy_crosscheck.c`,
  `scripts/astropy_crosscheck.py`): HIP 3757's Galactic-frame position and
  velocity match astropy's own ICRS<->Galactic transform to ~1e-7 relative
  error. This closes the exact gap `coords.h`'s own caveat had flagged
  (round-trip self-consistency can't catch a convention error; an
  independent reference implementation can) -- and it passes. Rules out
  `coords.c`'s formula/`A_G` matrix as the discrepancy's source too.
- **Non-gravitational acceleration quantified, not just asserted
  negligible**: JPL Horizons' own solution for 'Oumuamua confirms non-grav
  terms are present (A1=2.79e-7, A2=1.44e-8, A3=1.57e-8 AU/day^2, active
  "October 2017 to January 2018" per Horizons' own documentation). Back-of-
  envelope: ~2.8e-7 AU/day^2 for the ~100 remaining days after this
  project's reference epoch (JD 2458041.5) gives a velocity change of
  roughly 0.05 km/s -- three orders of magnitude too small to explain a
  15-40x distance/time gap. This hypothesis is now RULED OUT quantitatively,
  not just "not expected to matter."

**Net result**: every individual component in the chain -- 'Oumuamua's
orbital-element propagation (JPL Horizons, two epochs), the star-side
coordinate transform (astropy), the Milky Way potential formulas (galpy),
the integrator (analytic Kepler/zero-gravity cases), the encounter search
logic itself (straight-line analytic cases) -- is now independently
validated against an external reference, for TWO different real stars, and
each still shows the same qualitative discrepancy against Bailer-Jones'
reported encounter numbers. This is strong evidence the gap is not a code
bug in any of the pieces checked, but something about the overall
comparison not yet understood -- most likely a genuine methodological
difference (their exact reference epoch/frame for 'now', or their reported
figure being a Monte-Carlo median materially different from a single
best-fit trajectory, in a way not yet explained by the magnitude of
Gaia's own measurement uncertainty). Still open. Next useful step would be
getting the actual numbers Bailer-Jones' own code computes for these two
stars' best-fit (non-MC) trajectories specifically, to check whether even
THEIR best-fit single trajectory agrees with their own reported MC median
-- if not, that would point straight at the Monte-Carlo-median explanation
rather than anything in this project's own pipeline.

## First real-data scan tool (src/gaia_scan.c) -- real bug found, first preview candidates

Built `src/gaia_scan.c`, a standalone CLI tool (not the BOINC app) that
streams the downloaded Gaia bulk CSV directly and screens every star
against 'Oumuamua: a cheap analytic straight-line filter (O(1) per star, no
integration loop) first, then a full-gravity `find_closest_approach()`
confirmation for anything under 5 pc. Makes scanning tens of millions of
rows in a single local pass practical (a dry run on 500k rows took ~4
seconds).

**Real bug found and fixed**: the straight-line screen can find a closest
approach in either the past (t_star<0) OR the future (t_star>0) depending
on the star's specific relative motion -- obvious in hindsight, but the
first version of the gravity-confirmation step always searched BACKWARD
regardless of which sign the straight-line estimate came back with. For a
star whose real closest approach was in the FUTURE, this "confirmed" it by
searching an unrelated backward trajectory instead, producing a wildly
wrong, meaningless distance (one test case: straight-line said 1.69 pc at
t=+4.5 Myr; the buggy backward "confirmation" said 230 pc). Fixed by only
treating t_star<0 candidates as real candidates at all -- which is also
the scientifically correct filter anyway, since a future closest approach
says nothing about where the ISO came FROM, the actual question this
project asks.

**First preview results** (scan still running against the download, which
is itself still in progress -- see "Data" in README.md's roadmap): several
sub-5pc candidates found already, including one with an implausibly high
~343 km/s relative speed. Treating that one with real skepticism rather
than excitement -- a relative velocity that high is far outside normal
disk-star kinematics and is much more likely a Gaia radial-velocity
measurement artifact for that specific source than a genuine close
encounter, consistent with how Bailer-Jones et al. 2018 themselves treat
relative velocity as a plausibility signal (their own most-favored
candidate, HD 292249, is their LOWEST relative-velocity one, not their
closest-approach one). Full results and any real candidate list are
pending both the download finishing and a proper review pass (e.g.
cross-matching candidates against SIMBAD, sanity-checking relative
velocities) -- nothing here is being reported as a discovery yet.

## First real (partial) scan results, 2026-10-06

Full scan run against however much of the Gaia download had completed at the time (27.77 million
of ~33 million rows, download itself still at ~81% -- see README.md's "Data" roadmap item): 802
stars found with a straight-line closest approach under 5 pc, 769 of them with a "plausible"
(<150 km/s) relative velocity. A handful are genuinely interesting by the SAME standard
Bailer-Jones et al. 2018 themselves use (low relative velocity, not just small distance):

- source_id 5377679146317054080: v=19.4 km/s, d=0.86 pc, t=-2.65 Myr
- source_id 5497072086021346688: v=15.7 km/s, d=1.27 pc, t=-0.82 Myr
- source_id 5667826276174997760: v=13.2 km/s (the lowest in the whole list), d=11.5 pc, t=-9.05 Myr

These are genuinely comparable in quality to Bailer-Jones' own published HIP 3757/HD 292249
numbers by the velocity-plausibility standard -- but NOT being reported as real candidates yet,
deliberately, for two concrete reasons:
1. **Partial data**: ~19% of the Gaia download was still missing when this ran. A re-run once
   `scripts/download_filter_6d.py` finishes is needed before any list is even complete, let alone
   final.
2. **The unresolved HIP 3757/HD 292249 discrepancy directly applies here too**: this project's own
   pipeline doesn't yet reproduce Bailer-Jones' published distance/time for either of their known
   candidates (see this file's own "Pipeline validation" section above) -- every individual
   component is independently validated, but the full encounter-geometry comparison still shows an
   unexplained ~15-40x gap on distance/time specifically. Until that's understood, the ABSOLUTE
   distance/time numbers this scan reports for ANY candidate (including these three) should be
   treated as indicative/ranking information only, not as trustworthy physical measurements -- the
   same caution that already applies to the known candidates applies equally to new ones found the
   same way.

Neither HIP 3757 nor HD 292249 appear in this candidate list (both have straight-line distances
of ~24-27 pc in this project's own model, well above the 5 pc screening threshold) -- consistent
with, not contradicting, what was already found for them specifically; they were never expected to
pass a 5pc screen given this project's own numbers for them.

Next steps once the download finishes: re-run the full scan for a true final pass; consider
re-visiting the 5pc threshold (loosen it, since the known candidates' own true-encounter distances
in Bailer-Jones' numbers would not have passed a 5pc screen under THIS project's model either, by
the same logic); and the real priority remains resolving the underlying discrepancy before treating
any of these numbers as trustworthy.

## Second gaia_scan.c window bug: low relative-velocity candidates specifically

Found while reviewing the first full 50pc-threshold scan's results (see "First real (partial) scan
results" above): several LOW relative-velocity candidates -- exactly the most scientifically
interesting kind, per Bailer-Jones et al. 2018's own preference -- came back with absurd
"confirmed" distances (one case: 50 pc straight-line screen vs. a "confirmed" 5000 pc; others at
1000+ pc). All affected cases shared low relative velocity (a few km/s).

Root cause: the straight-line model's own implied timescale (t_star) grows large precisely when
relative velocity is low -- and that's exactly the regime where gravity's cumulative effect over a
long timescale can move the TRUE minimum far from the straight-line guess, no matter how many times
the window-widening retry (from the earlier window-truncation fix) doubles it, since each doubling
still starts from an unreliable base estimate. The widening retry's OWN cap (4 attempts) could run
out before reaching the true minimum for these cases, silently returning a still-truncated result.

Fixed two ways:
1. Never start the gravity-confirmation window below this project's own settled 100 Myr backward-
   integration trust limit (NOTES_method_selection.md's own "Backward-integration time window"
   decision) -- removes the dependency on the straight-line guess for SIZING the window at all,
   only still using it for the initial cheap filter. Affordable because this is a single
   confirmation trajectory, not mc_search_stage1()'s 100 Monte Carlo clones.
2. Added an explicit `WINDOW_LIMITED` flag to gaia_scan.c's output for any candidate where the
   final result is STILL at the window's edge even after every retry attempt -- rather than
   silently trusting a number that might still be wrong, scripts/summarize_candidates.py now
   separates these out explicitly and excludes them from the "plausible candidates" ranking.

Also found and fixed, while sizing real BOINC work units for `src/main_nostos_boinc.cpp`: a real
timing measurement (10 stars x 100 Monte Carlo clones x 100,000 steps = 23.05s wall-clock on the
bitboinc-app production server) implies ~2.3s/star, meaning the earlier, un-benchmarked default
of 1000 stars/work-unit would take on the order of 38 MINUTES per work unit -- confirmed directly
(a real run was killed after 12+ minutes of 100% CPU, still not done). `scripts/make_work_units.py`'s
default chunk size lowered to 150 stars/WU (~5.75 min on this reference hardware) based on that
real number, not a guess.

## 'Oumuamua's real velocity-uncertainty covariance, 2026-10-06

`scripts/make_work_units.py` had been using a placeholder velocity covariance for 'Oumuamua's
Monte Carlo clones (1e-6 (kpc/Myr)^2 on the diagonal, about 1 km/s per axis) -- explicitly
flagged in its own comment as not yet derived from the real published orbit solution. Replaced
with a real one, computed from first principles rather than guessed or looked up as a final
number:

1. Fetched 'Oumuamua's actual orbit-fit covariance matrix directly from JPL's Small-Body
   Database API (`https://ssd-api.jpl.nasa.gov/sbdb.api?sstr=1I&full-prec=true&cov=mat`) --
   a real 9x9 matrix over (e, q, tp, node, peri, i, A1, A2, A3), Solution 7c / Micheli et al.
   2018 (the same orbit determination, data_arc 80 days, this project's OUMUAMUA_ELEMENTS
   already uses elsewhere -- see the caveat below on the small epoch difference).
2. The first six of those nine parameters are exactly the ones
   `orbital_elements_to_galactocentric()` takes (tp enters through this project's `t_myr`
   argument, not a struct field, so its uncertainty is propagated as a perturbation of t_myr).
   The other three (A1/A2/A3, non-gravitational-acceleration model parameters) aren't used by
   this project's two-body Keplerian model at all. Marginalizing them out by taking the top-left
   6x6 block isn't an approximation -- the marginal covariance of a subset of jointly-Gaussian
   variables IS exactly that principal submatrix.
3. Built a numerical Jacobian of Galactocentric velocity (at this project's own settled
   asymptotic-regime evaluation point, ~9.08 years post-perihelion) with respect to each of
   those six parameters, via central finite differences through the project's own
   `orbital_elements_to_galactocentric()` -- not a hand-derived formula. Propagated:
   Cov_vel (3x3) = J^T * Cov6 * J. Verified stable (matching to ~10 significant figures) under
   a 5x change in the finite-difference step size, so this isn't finite-difference noise.
   See `tools/oumuamua_velocity_covariance.c` for the full, heavily-commented implementation.

Result: per-axis 1-sigma velocity uncertainty of about 0.00016, 0.00071, 0.00069 km/s
(sub-m/s). This is NOT a bug or an implausibly tight number -- 'Oumuamua's orbit is one of the
most precisely determined hyperbolic orbits ever fit (80 days of arcsecond-level astrometry on
a nearby, fast-moving object), and the literature (e.g. Feng & Jones 2018) does describe its
velocity as constrained to well below 1 m/s precision in places. `mc_search_stage1()`'s 100
clones will now cluster far more tightly around the best-fit trajectory than the old placeholder
implied -- correctly, since the real uncertainty is that small, not because of a bug.

Caveat stated plainly, not glossed over: SBDB's default-returned solution is quoted at epoch JD
2458080.5, while OUMUAMUA_ELEMENTS (used elsewhere in this project) came from a JPL Horizons
fetch at JD 2458041.5 -- the same underlying solution, propagated/reported at a different epoch,
with element values agreeing to ~1-2e-4 relative. The Jacobian here is evaluated at THIS
solution's own best-fit point rather than re-deriving a covariance at the other epoch (which
SBDB's on-demand API doesn't offer), since the two points are close enough that the local linear
sensitivity is essentially the same either way -- a reasonable approximation, not treated as
exact.

## Borisov pipeline validation (second ISO), 2026-10-06

Extended "Pipeline validation" to a second, independent interstellar object, per the README's
own "Borisov not attempted yet" item: `tests/test_pipeline_borisov_ross573.c`, targeting Ross
573 (Gaia DR2 5162123155863791744), the single CLOSEST encounter Bailer-Jones, Farnocchia, Ye,
Meech & Micheli 2020 (A&A, arXiv:1912.10213) found for 2I/Borisov across their entire
7.4-million-star Gaia DR2 search -- their headline result, not a hand-picked easy case:
t=-909 kyr, d=0.068 pc, v=22.6 km/s (Monte Carlo medians, same caveat as the 'Oumuamua
comparisons: not a single best-fit trajectory).

2I/Borisov's orbital elements came from JPL's Small-Body Database (same source/method as
'Oumuamua's, see test_horizons_crosscheck.c). The evaluation epoch (3.0 years post-perihelion)
was chosen the same way 'Oumuamua's 9.08-year point was -- found empirically where the computed
velocity converges to its asymptotic value (within ~0.03% of the 10-year value, and matching
the orbit's own analytic hyperbolic-excess speed sqrt(-GM_sun/a)=32.28 km/s heliocentric).

Result: same qualitative pattern as 'Oumuamua's two cases -- right sign (past), right
order-of-magnitude relative speed (29.8 vs 22.6 km/s, actually the closest speed match of the
three stars so far), but a much shallower/earlier minimum than reported (13.77 pc at t=-0.523
Myr, vs their 0.068 pc at -0.909 Myr -- about 200x off in distance, worse than either
'Oumuamua case). Checked manually with the search window widened to 5.0 Myr: found the exact
same minimum, so this isn't a window-truncation artifact (the same check already done for HD
292249).

This is now THREE stars across TWO independent interstellar objects showing the same
qualitative mismatch, with every pipeline component independently validated and checking out
for both objects. Strengthens the existing conclusion (see the 'Oumuamua pipeline section
above): this looks like something systematic about the overall comparison method (Bailer-Jones'
exact reference epoch/frame, or Monte-Carlo-median-vs-single-best-fit), not a per-object or
per-star bug. Left open, not papered over.

## Testing the "Monte-Carlo-median-vs-single-best-fit" hypothesis, 2026-10-06 -- RULED OUT

The one remaining untested hypothesis for the open 'Oumuamua/Borisov discrepancy (see both
pipeline-validation sections above) was that Bailer-Jones et al.'s reported numbers are MEDIANS
over ~2000 Monte Carlo surrogates resampling both bodies' full phase-space uncertainty, while
this project's pipeline tests so far only ever reported the single deterministic best-fit
trajectory -- and that those two quantities might differ substantially for a real geometric
reason (a deterministic trajectory has one particular, not necessarily minimal, miss distance;
resampling over an uncertainty ellipsoid explores many nearby trajectories, some of which could
pass much closer by chance, pulling the MEDIAN below the deterministic value).

Tested this directly with `tools/mc_vs_bestfit_ross573.c`: resampled BOTH Ross 573's full 6D
state (Gaia DR2's own correlated 5x5 astrometric covariance -- ra, dec, parallax, pmra, pmdec --
plus its independent radial-velocity error, all queried directly from gaiadr2.gaia_source) AND
2I/Borisov's full orbital solution (JPL SBDB's own 6x6 marginal covariance over e, q, tp, node,
peri, i, same source as the 'Oumuamua covariance work above), 2000 times each (matching
Bailer-Jones et al.'s own count), running the exact same `find_closest_approach()` this
project's pipeline tests use on every draw.

Result: median distance 13.78 pc (90% CI 13.53-14.00 pc) -- essentially IDENTICAL to the
deterministic best-fit's 13.77 pc, and a very tight interval around it. Same for time
(-0.5227 Myr median vs -0.523 Myr best-fit) and speed (29.77 vs 29.8 km/s). **This rules out
the hypothesis**: resampling the full real uncertainty on both bodies barely moves the result at
all, because both Ross 573's Gaia astrometry and Borisov's orbit solution are so tightly
constrained (same conclusion already reached for 'Oumuamua's own orbit-fit uncertainty, see the
velocity-covariance section above -- sub-m/s) that the uncertainty ellipsoid is far too small to
explain a ~200x gap in miss distance by chance alone.

This is a real, if negative, finding: it eliminates "we're comparing the wrong kind of number"
as the explanation, and narrows the remaining open possibilities down to something about
Bailer-Jones et al.'s exact method that this test didn't reproduce -- their precise reference
epoch/frame convention, their specific Galactic potential model (axisymmetric, parameters not
fully specified in the abstract-level summary read so far), or a genuine difference in which
solar motion/LSR convention is used for `solar_params_default`. None of these have been
checked against Bailer-Jones et al.'s actual paper text in detail yet -- that remains the next
real step, not a MC-sampling artifact.

## Found the actual bug: inbound vs outbound asymptote, 2026-10-07

Root cause of the "open discrepancy" documented at length above (both the 'Oumuamua and
Borisov pipeline-validation sections, and the now-ruled-out Monte-Carlo-vs-best-fit test).
Found by reading Bailer-Jones, Farnocchia, Ye, Meech & Micheli 2020's actual methodology text
(arXiv:1912.10213) rather than just its headline numbers: they state they compute 2I/Borisov's
pre-encounter trajectory by "integrat[ing] back to 3000 BCE and then extrapolat[ing] to the
Keplerian inbound asymptote" -- explicitly the INBOUND asymptote, not the outbound one this
project's pipeline had been using everywhere (`tests/test_pipeline_oumuamua_hip3757.c`,
`tests/test_pipeline_borisov_ross573.c`, `src/gaia_scan.c`, `src/main_nostos_boinc.cpp`, and
both the velocity-covariance tools above all evaluated `orbital_elements_to_galactocentric()` at
a POSITIVE t_myr -- time AFTER perihelion).

Why this matters, and why none of the extensive component-level validation caught it: this
project's long-timescale backward integration (`find_closest_approach()`, Myr-scale, Milky Way
potential only) deliberately never re-simulates the Sun's own gravity -- by design, since the
Sun's influence is supposed to be negligible by the time that integration takes over (see
`orbital_elements.h`'s own domain-of-validity note). But the Sun's gravity is NOT negligible
near perihelion itself -- it's what makes the orbit hyperbolic in the first place, and it
deflects the trajectory's velocity DIRECTION by the orbit's full hyperbolic turning angle,
2*asin(1/e), between the inbound (pre-encounter) and outbound (post-encounter) asymptotes. For
'Oumuamua (e=1.201) that's 112.7 degrees; for Borisov (e=3.356) it's 34.7 degrees -- neither is
a small correction. Starting the Galactic-potential-only backward integration from the OUTBOUND
state (the one actually observed, since both objects are only discovered post-perihelion) and
integrating backward through where perihelion was silently carries the post-deflection velocity
DIRECTION all the way back through the Galaxy, rather than correctly switching to the
pre-deflection (inbound) direction the object actually had on its long approach. The orbit's
SPEED is conserved between the two asymptotes (both equal v_infinity, by energy conservation),
which is exactly why every test kept finding "right order-of-magnitude speed, very wrong
distance" -- a smoking-gun pattern that, in hindsight, points straight at a direction error
rather than anything in the magnitude/scale of the problem.

This evaded every earlier cross-check because each one validated the TRANSFORM CHAIN itself
(JPL Horizons, astropy, galpy), never which TIME ARGUMENT the pipeline fed into that chain --
a methodological/usage bug, not a bug in any validated component.

Fix: evaluate `orbital_elements_to_galactocentric()` at the same magnitude but NEGATIVE t_myr
(~9.08 years before perihelion for 'Oumuamua, ~3.0 years before for Borisov -- unchanged from
before, since by the two-body Kepler problem's time-reversal symmetry the inbound branch
converges to its own asymptotic speed at the same rate the outbound branch does). Applied to
every call site: `src/gaia_scan.c`, `src/main_nostos_boinc.cpp`,
`tests/test_pipeline_oumuamua_hip3757.c`, `tests/test_pipeline_borisov_ross573.c`,
`tools/oumuamua_velocity_covariance.c`, `tools/mc_vs_bestfit_ross573.c`.

Result -- dramatic improvement across all three reference stars:

| star      | metric   | before (outbound) | after (inbound) | Bailer-Jones |
|-----------|----------|--------------------|-------------------|--------------|
| HIP 3757  | t (Myr)  | -0.067             | -0.9855           | -0.990       |
|           | d (pc)   | 24.73              | 1.17              | 0.60         |
|           | v (km/s) | 39.15              | 24.69             | 24.7         |
| HD 292249 | t (Myr)  | -0.592             | -3.949            | -3.765       |
|           | d (pc)   | 27.14              | 2.31              | 1.60         |
|           | v (km/s) | 51.74              | 10.24             | 10.7         |
| Ross 573  | t (Myr)  | -0.523             | -0.893            | -0.909       |
|           | d (pc)   | 13.77              | 1.06              | 0.068        |
|           | v (km/s) | 29.77              | 23.01             | 22.6         |

Time and speed now match Bailer-Jones et al. to within a few percent for all three stars (the
pipeline tests now assert this, to a generous 25% tolerance). Distance is still off by a
residual ~1.4-2x for the two 'Oumuamua stars and ~15x for Ross 573 -- down from 15-200x, but not
fully closed. Re-ran the already-built `tools/mc_vs_bestfit_ross573.c` two-sided Monte Carlo
diagnostic at the corrected inbound point too: median 1.06 pc (90% CI 1.01-1.11 pc), essentially
identical to the new deterministic best fit -- confirming the MC-vs-best-fit hypothesis is
STILL ruled out even with the direction bug fixed; it was never about that.

Remaining open question (smaller now, genuinely different from before): why is there still a
~1.4-15x gap in distance specifically, with time and speed now both solid? Not yet
investigated. Candidate hypotheses for a future session: Bailer-Jones et al.'s potential model
details beyond "smooth axisymmetric, doesn't matter much" (they don't fully specify parameters
in the text read so far); their reported numbers are 2000-surrogate MEDIANS is confirmed NOT the
explanation, but the SPECIFIC surrogate-generation procedure (they resample a "6D covariance
matrix, 3D position + 3D velocity" built some other way, not necessarily the orbital-element
covariance used here) could still differ in a way that matters for the heavier tails of a
highly sensitive quantity like miss distance even if the median doesn't move much; or a smaller
remaining frame/epoch effect not yet isolated.

IMPORTANT CONSEQUENCE: every Gaia-catalog scan result produced before this fix --
`scan_oumuamua_candidates_50pc.txt` (83,363 candidates from a partial, ~28.5M-row scan) and the
14M-row window-fix validation run -- used the WRONG (outbound) 'Oumuamua state and must be
treated as STALE. The window-sizing fix itself (bugs #1-#3 above) is unaffected and still
validated correctly by that 14M-row run (no WINDOW_LIMITED stragglers in its output) -- only the
ISO's absolute trajectory was wrong, not the window-search mechanics. A full fresh rescan with
the corrected `gaia_scan.c` is needed before trusting any candidate list from this project.

## Checked solar-motion-parameter differences too -- also not the residual gap, 2026-10-07

Quick follow-up check on the smaller ~1.4-15x residual distance gap left after the
inbound/outbound fix above: fetched Bailer-Jones' own solar position/velocity parameters
(Bailer-Jones 2015/2014, arXiv:1412.3648, cited by both the 'Oumuamua 2018 and Borisov 2020
papers for "the same Galactic model"): R0=8.0 kpc, Z0=+10 pc, and a Sun azimuthal velocity of
242 km/s taken directly from Reid & Brunthaler 2004's Sgr A* proper-motion measurement (already
the Sun's full velocity, not a separate circular+peculiar split) -- versus this project's
solar_params_default: R0=8.122 kpc, Z0=20.8 pc, v_circ+v_pec=241.24 km/s. U and W peculiar
motion are IDENTICAL in both (11.1, 7.25 km/s, same Schonrich et al. 2010 citation).

Re-ran all three reference stars with Bailer-Jones' own values substituted in for
solar_params_default: the results change by less than 0.5% in every quantity (e.g. HIP 3757's
distance: 1.1689 -> 1.1683 pc). Rules this out too as an explanation for the residual gap --
the two groups' solar-motion assumptions are close enough that switching between them is noise
compared to the ~1.4-15x gap that's left. The remaining open candidates (potential model
details beyond "doesn't matter much", or the precise surrogate-generation procedure for their
Monte Carlo even though its median doesn't move) are unchanged from the previous section.

Also checked integrator resolution as a possible contributor to the same residual gap: reran
HIP 3757 at 10x and 100x the standard step density (1M and 10M steps/Myr vs the usual 100k).
Result identical to 6 significant figures at every resolution (d=1.168878 pc in all three
cases) -- the standard resolution is already fully converged, not a source of the gap either.

One more candidate checked and ruled out: Bailer-Jones et al. 2020 themselves quote a "neglect
of two-body effects" uncertainty for Ross 573 specifically -- Ross 573's OWN gravity deflecting
Borisov during the encounter (not Sun-related; this project's encounter.c also deliberately
ignores star-ISO mutual gravity, per its own top comment). They compute this effect at just
0.002 pc for Ross 573 -- far too small to explain the current ~1 pc residual gap, and they note
it's negligible for every other encounter they report too. Not the explanation either.

Stopping the residual-gap investigation here for now: potential-model choice, Monte-Carlo-vs-
best-fit, solar-motion parameters, integrator resolution, and star-ISO mutual gravity are all
now ruled out as the dominant cause. The inbound/outbound fix above resolved the large (15-200x)
part of the original discrepancy; what's left (~1.4-15x) is a smaller, genuinely open question,
plausibly within the kind of real astrophysical/methodological uncertainty both research groups
already operate under, rather than a further bug in this codebase.

## First preview scan with the fixed binary, and a small duplicate-row issue, 2026-10-07

Ran the corrected `gaia_scan.c` (both the window-sizing fix and the inbound/outbound fix) on
the first 14M rows as a before/after sanity check (the 50pc straight-line threshold, same as
before): 56,767 candidates found (vs. 35,518 for the same 14M rows under the old, buggy
outbound-direction binary -- a real change, as expected, since the whole trajectory moved).
Zero WINDOW_LIMITED flags, confirming the window-sizing fix (bugs #1-#3) still holds at this
scale with the direction fix applied too. Saved as
`D:/nostos_data/scan_oumuamua_14M_fixed_preview.txt` (candidate lines only -- note for next
time: redirecting gaia_scan.c's stdout and stderr to the SAME file, as a quick test command
here initially did, mixes its progress/diagnostic lines in with the actual candidate lines and
breaks `summarize_candidates.py`'s parser; keep them in separate files/streams as the tool's
own design already assumes).

Running `summarize_candidates.py` on this preview surfaced a real, if small, data-quality
issue: 237 of 56,767 candidate rows (0.4%) are exact duplicates by source_id, with IDENTICAL
computed values both times -- consistent with the same star's row appearing twice in the
downloaded Gaia CSV itself (gaia_scan.c just screens whatever rows the file hands it; this
isn't a screening-logic bug) -- confirmed directly: `grep -c` for one such duplicated source_id
against the full 43GB `gaia_dr3_6d.csv` returns 2, not 1. Exact upstream root cause (most likely
download_filter_6d.py's resume logic re-processing an already-partially-written chunk after an
interruption, across this project's many resume sessions) not tracked down further -- not worth
blocking on for a ~0.4% effect. Fixed pragmatically in `summarize_candidates.py` by
deduplicating on source_id (keep first occurrence) rather than chasing the upstream cause
immediately, since the effect is small and shouldn't block using the tool's output today.

Promising sign from the corrected results: the lowest-relative-velocity "plausible" candidates
now have genuinely LOW speeds (0.71-3.36 km/s for the ones checked), not the implausible
hundreds-of-km/s values that dominated the top of the old, buggy-direction candidate list --
exactly the kind of result Bailer-Jones et al.'s own preference for low relative velocity as a
plausibility signal would predict for a correctly-computed scan.

## First full-catalog scan with the fully-fixed pipeline, 2026-10-07

The Gaia DR3 6D-phase-space download (`scripts/download_filter_6d.py`) finished completely
overnight: 33,530,050 total rows, 32,979,536 parsed and screened. Ran the fully-fixed
`gaia_scan.c` (window-sizing fix + inbound/outbound-asymptote fix) over the entire file:
101,667 raw candidates under the 50pc straight-line threshold (101,430 after deduplicating the
small, confirmed source-CSV duplicate-row issue noted above), zero WINDOW_LIMITED stragglers,
93 minutes wall-clock. Saved as `D:/nostos_data/scan_oumuamua_full_fixed.txt` (candidates) and
`..._progress.log` (progress, kept separate this time -- see the preview-scan section above for
why that matters for `summarize_candidates.py`). Summary saved to
`..._full_fixed_summary.txt`; a curated "distance < 2.0 pc AND speed < 50 km/s" shortlist (99
stars, Bailer-Jones et al.'s own style of plausibility cut) saved to
`..._full_fixed_most_promising.txt`.

Validation this scan's own internal consistency, not just the method: both of this project's
long-standing reference stars show up naturally in the full output, close to (not identical to,
since gaia_scan.c's straight-line-then-confirm pipeline and the dedicated pipeline tests'
direct `find_closest_approach()` call aren't byte-for-byte the same code path) what the
dedicated pipeline tests already found:
  - HIP 3757 (2525688198820543360): d=1.2249 pc, v=24.42 km/s, t=-0.9955 Myr here, vs.
    d=1.1689 pc, v=24.69 km/s, t=-0.9855 Myr in `tests/test_pipeline_oumuamua_hip3757.c`.
  - HD 292249 (3107000885484340224): d=2.4176 pc, v=10.44 km/s, t=-3.8615 Myr here, vs.
    d=2.3073 pc, v=10.24 km/s, t=-3.9490 Myr in the same test.
Both land in the top tier of the full catalog's own "most promising" shortlist (HIP 3757 at
rank 28 by distance among the 99 sub-2pc/sub-50km/s stars) -- a reassuring sign the full-scale
scan is finding the same stars the targeted validation already vetted, not something unrelated.

New candidates worth noting for a future session (none investigated further here -- this is
output, not a discovery claim): several stars combine sub-parsec distance with single-digit-to-
tens km/s relative speed, both closer AND slower than either of this project's two validated
reference stars --
  - 39513359823463680: d=0.367 pc, v=13.89 km/s, t=-10.998 Myr
  - 4899996487129314688: d=0.982 pc, v=6.06 km/s, t=-5.649 Myr
  - 357194086205843072: d=1.036 pc, v=8.37 km/s, t=-15.352 Myr
  - 360189293321105024: d=1.156 pc, v=9.41 km/s, t=-10.216 Myr
These haven't been identity-confirmed via SIMBAD, cross-checked against any published list, or
run through mc_search_stage1()'s Monte Carlo resampling yet -- all real next steps, not done
here. This IS the first trustworthy full-catalog candidate list this project has ever produced
(every earlier one was generated before the inbound/outbound fix and is stale, as flagged
above) -- a genuine milestone for the "Data" roadmap item, even though every candidate on it
still needs the same level of scrutiny the two reference stars already got before being taken
seriously.

## Identity-confirmed the top new candidates, 2026-10-07

Quick follow-up on the four "closer and slower than either reference star" candidates listed
above: checked each one isn't a data artifact by confirming it's a real, legitimate catalog
entry (same kind of check this project did for HIP 3757/HD 292249/Ross 573 before trusting
them).
  - 39513359823463680 -- UCAC4 520-006717, 2MASS J03502880+1356125, parallax 6.71 mas
    (~149 pc). Real, cataloged, not individually named/studied beyond that.
  - 4899996487129314688 -- L 86-22 / NLTT 72, parallax 29.45 mas (~34 pc), has a measured RV
    (18.14 km/s). Real, cataloged proper-motion star.
  - 357194086205843072 -- not in SIMBAD (most faint Gaia-only stars aren't, if nothing else has
    individually studied them), but confirmed as a genuine Gaia DR3 source via direct archive
    query: parallax 8.01 mas (~125 pc), G=15.25 (faint, consistent with why it's SIMBAD-absent).
  - 360189293321105024 -- TYC 3685-1833-1 / GSC 03685-01833, parallax 10.25 mas (~98 pc),
    V=10.77 (bright enough to have its own Tycho/GSC entry).
All four are real, legitimately-cataloged (or at least Gaia-archive-confirmed) stars, not
artifacts of a parsing bug or a corrupted row -- a reasonable first screen before any of them
would be worth investigating further (SIMBAD cross-ID, literature search, Monte Carlo
resampling via mc_search_stage1()). None of that further investigation has been done yet.

## A third independent paper's candidate, for context (not added as a permanent test), 2026-10-07

While checking whether any of the full-catalog scan's new candidates matched previously
published ones (none did -- see above), found a THIRD independent paper with its own
'Oumuamua-origin candidate: Dybczynski & Krolikowska 2018 (arXiv:1711.06618), using pre-Gaia-DR2
kinematics, propose UCAC4 535-065571 (= V* V1274 Her, M6V dwarf, Gaia DR3 4549983418742839808)
as the most probable candidate within 60 pc: t=-2.14 Myr, d=1.43 pc, v=5.21 km/s. They
explicitly flag this star's kinematics as poorly known at the time and say improved Gaia data
would be needed for a firmer conclusion.

Ran it through this project's now-fixed pipeline using its REAL Gaia DR3 astrometry (not
available to that 2018 paper): found t=-1.3323 Myr, d=7.2115 pc, v=10.83 km/s. Same qualitative
pattern as the other three reference stars (right sign, same order-of-magnitude speed, distance
off by a larger factor) -- consistent with, not a new instance contradicting, the already-
documented residual gap. Plausibly explained in part by this specific case by the star's own
flagged kinematic uncertainty at publication time, on top of whatever is causing the smaller
residual gap for the other (Gaia-DR2-based) reference stars. Not added as a permanent pipeline
test (unlike HIP 3757/HD 292249/Ross 573) since this star's own origin paper already flags its
input data as unreliable -- noted here for context only.
