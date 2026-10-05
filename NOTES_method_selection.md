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

## Milky Way potential (NOT settled — placeholder in place)

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

## Coordinate transform (not started)

Converting Gaia's (RA, Dec, parallax, proper motion, radial velocity) and an
ISO's heliocentric orbital elements into a common Galactocentric Cartesian
frame isn't written yet. Needs: solar position (R0, z0) and velocity
(circular + peculiar motion) constants — common modern choices are
R0=8.122 kpc (GRAVITY Collaboration 2018) and peculiar motion
(U,V,W)=(11.1, 12.24, 7.25) km/s (Schönrich et al. 2010), but not yet decided
for this project specifically (should match whatever the final Milky Way
potential model assumes for its own R0/v_c, to stay self-consistent).
