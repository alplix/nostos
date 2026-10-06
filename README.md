# Nostos

**Nostos** (νόστος, "homecoming") is a volunteer-computing project, sharing [bitboinc](https://bitboinc.athena.org.tr)'s BOINC infrastructure, that searches for the **origin stars of interstellar objects** — asteroids and comets like 1I/'Oumuamua and 2I/Borisov that passed through the Solar System on hyperbolic trajectories, having come from somewhere else in the Milky Way.

## The idea

An interstellar object's current position and velocity, projected backward through the Galaxy's gravitational potential, traces a path through where it's been. Somewhere along that path, it may have had a close encounter with a star — close enough, long enough ago, to plausibly be where it was ejected from.

A handful of papers (most notably Bailer-Jones et al. for 'Oumuamua and Borisov) have done this kind of search, but against a hand-picked shortlist of candidate stars — a few dozen, chosen because they looked promising. Gaia's full catalog has tens of millions of stars with complete 6D phase-space data (position + velocity). Nobody has run the exhaustive version: every known interstellar object, against *every* Gaia star with full 6D data, with proper Monte Carlo propagation of each object's own orbital uncertainty. That's the computational gap this project fills — and it isn't a one-off: Rubin Observatory/LSST is expected to find many more interstellar objects in the coming years, so this is lasting infrastructure, not a single analysis.

## How it works (planned)

1. Each known interstellar object's orbital elements (public, from JPL Small-Body Database / Minor Planet Center) are resampled into many Monte Carlo realizations spanning its measurement uncertainty.
2. Gaia DR3's ~33 million full-6D-phase-space stars are split into chunks and distributed to volunteers as BOINC work units.
3. Each work unit integrates its chunk of Gaia stars *and* the interstellar object's Monte Carlo realizations backward through a standard Milky Way gravitational potential (no star-star N-body needed — each trajectory is independent in a smooth background potential), looking for closest-approach distance and time.
4. Close encounters below a threshold distance are reported back; everything else is discarded to keep result sizes small.
5. The pipeline is validated against the *already-published* candidate origin stars for 'Oumuamua and Borisov before being trusted on anything new — if it can't reproduce known results, it doesn't get to report new ones.

## Status

Gaia DR3 6D-phase-space download in progress (see `scripts/download_filter_6d.py`,
resumable via `scripts/progress.txt`). Core orbit integrator, the
astrometry<->Galactocentric coordinate transform, the closest-approach
search kernel, and the orbital-elements starting-state converter are all
written and validated (see `src/` and `tests/`, and `NOTES_method_selection.md`
for the choices made so far, two real bugs found and fixed along the way,
and what's still open). A CPU BOINC worker app exists and has been
compiled and functionally tested (`src/main_nostos_boinc.cpp`) -- not yet
deployed; no work generator, validator, or server-side registration yet.

## Roadmap

- [x] **Method selection** — integrator, coordinate transform, orbital-elements conversion, a validated Milky Way potential (Allen & Santillan 1991, published parameters), the backward-integration time window (100 Myr, with results beyond 50 Myr flagged lower-confidence), and the Monte Carlo resampling scheme (100 clones per ISO, two-stage design for the star side) are all settled -- the last two by reading Bailer-Jones et al. 2019's actual published methodology rather than guessing. See `NOTES_method_selection.md`.
- [ ] **Data** — pull the Gaia DR3 full-6D-phase-space subset (sources with measured radial velocity) via the Gaia archive TAP service; keep only the columns actually needed (position, parallax, proper motion, radial velocity, and their errors).
- [x] **Core integrator** — `src/integrator.h` (leapfrog) + `src/potential.c` written and validated against exact conservation laws and the analytic Kepler solution (`tests/test_kepler.c`, all passing). The Milky Way potential is now Allen & Santillan (1991)'s actual published model (Plummer bulge + Miyamoto-Nagai disk + a truncated generalized-enclosed-mass halo), replacing the earlier generic Hernquist+MN+NFW placeholder -- bulge and disk are cross-checked against `galpy` component by component (`tests/test_milkyway_galpy.c`, ~1.7e-6 relative error, in-plane and off-plane); the halo (no stock galpy equivalent) is validated via energy/acceleration self-consistency and continuity across its truncation radius (`tests/test_milkyway_allen_santillan.c`). Resulting rotation curve: 203-232 km/s across 1-30 kpc. `src/coords.c` (astrometry <-> Galactocentric phase space) and `src/orbital_elements.c` (ISO orbital elements -> Galactocentric starting state) are written and tested (`tests/test_coords.c`, `tests/test_orbital_elements.c`, all passing) and cross-checked end to end against a real JPL Horizons state vector for 1I/'Oumuamua (`tests/test_horizons_crosscheck.c`) -- position and velocity both match to ~1e-6 relative error, confirming the full rotation chain's direction/sign convention, not just internal self-consistency.
- [x] **Encounter search kernel** — `src/encounter.c`'s `find_closest_approach()` integrates a star and an ISO realization together and tracks their minimum separation. Validated against an exact analytic (GM=0, straight-line motion) case (`tests/test_encounter.c`, all passing). `src/montecarlo.c` (a seedable PRNG + Cholesky-based correlated Gaussian sampler, validated statistically in `tests/test_montecarlo.c`) and `src/mc_search.c` (the stage-1 Monte Carlo wrapper around the above -- 100 ISO clones per star, per Bailer-Jones et al. 2019's own number, deterministic per seed for BOINC quorum validation; `tests/test_mc_search.c`) implement the Monte Carlo scheme decided in `NOTES_method_selection.md`.
- [~] **Pipeline validation** — tested against two of Bailer-Jones et al. 2018's own 'Oumuamua-encounter candidates, HIP 3757 and HD 292249 (`tests/test_pipeline_oumuamua_hip3757.c`). Found and fixed a real bug along the way (evaluating 'Oumuamua's orbital elements at perihelion instead of near "now" pointed its velocity in a meaningfully wrong direction, flipping whether the closest approach falls in the past or future). Every individual pipeline component is now independently validated against an external reference (galpy for the potential, astropy for the coordinate transform, JPL Horizons at two epochs for 'Oumuamua's propagation) and checks out -- yet both stars still show the same qualitative mismatch against the reported encounter distance/time (right sign and order-of-magnitude speed, ~15-40x off on distance/time). A genuinely open discrepancy, several hypotheses ruled out (see `NOTES_method_selection.md`). Borisov not attempted yet.
- [~] **BOINC app** — CPU worker app written (`src/main_nostos_boinc.cpp`): one work unit = one ISO against one chunk of Gaia stars, runs `mc_search_stage1()` per star, reports only encounters under the work unit's own distance threshold. Compiled and functionally tested end-to-end (statically linked against the sibling bitboinc-app project's own BOINC libs, real star data in, correct filtered output out) -- not yet deployed, no work generator or validator yet, no server-side registration (deliberately: those are production decisions, not made without the user's own review).
- [ ] **GPU acceleration** — the integration step parallelizes naturally across stars; port to CUDA/OpenCL/Metal once the CPU path is proven correct.
- [ ] **Publication** — results (and, eventually, any genuinely interesting close encounter) published openly, data releases via Zenodo.

## Why "Nostos"

The Odyssey's word for a hero's long journey home. Fitting for a project whose entire point is finding out where a wandering visitor from outside the Solar System actually came from.
