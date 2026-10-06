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
