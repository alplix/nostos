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

Just started. No code yet. See the roadmap below.

## Roadmap

- [ ] **Method selection** — settle on a specific Galactic potential model (e.g. `MWPotential2014`), the backward-integration time window (existing literature stays within a few Myr, for good reason — orbits become chaotic/unreliable further back), and the Monte Carlo sampling scheme for orbital uncertainty.
- [ ] **Data** — pull the Gaia DR3 full-6D-phase-space subset (sources with measured radial velocity) via the Gaia archive TAP service; keep only the columns actually needed (position, parallax, proper motion, radial velocity, and their errors).
- [ ] **Core integrator** — a small, fast C/C++ orbit integrator: analytic potential force law + a standard integrator (leapfrog or RK4). Differentially tested against a reference implementation (e.g. Python's `galpy`) on known test orbits before anything else is built on top of it.
- [ ] **Encounter search kernel** — the actual per-(Gaia star, ISO Monte Carlo sample) closest-approach calculation. Validated against hand-computed small test cases.
- [ ] **Pipeline validation** — run the full thing against 'Oumuamua and Borisov and confirm it independently finds the same candidate origin stars the published literature already identified.
- [ ] **BOINC app** — CPU first (matching this project family's usual order), work generator that chunks the Gaia catalog, assimilator that collects genuine close encounters.
- [ ] **GPU acceleration** — the integration step parallelizes naturally across stars; port to CUDA/OpenCL/Metal once the CPU path is proven correct.
- [ ] **Publication** — results (and, eventually, any genuinely interesting close encounter) published openly, data releases via Zenodo.

## Why "Nostos"

The Odyssey's word for a hero's long journey home. Fitting for a project whose entire point is finding out where a wandering visitor from outside the Solar System actually came from.
