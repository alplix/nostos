#!/usr/bin/env python3
"""Writes the input file `in` of a benchmark profile (scripts/bench/profiles/*.wu), the same bytes on every machine.

usage: make_bench_wu.py <profile.wu> <out_file>

The Gaia data of a real work unit is not in the repository, so the stars are made here. Time per star does not depend
on the star (every star runs all MC_N_ISO_CLONES clones for every one of the n_steps steps, with no early exit), so a
made-up star costs what a real one costs; what the data does decide is which stars come under the distance threshold
and are written to `out`. That is why the profiles carry two real stars (HIP 3757 and HD 292249, whose astrometry is in
tests/test_pipeline_oumuamua_hip3757.c): at the production threshold they are the rows of `out`.

Layout of `in`: scripts/make_work_units.py's own writer (write_work_unit), with its 'Oumuamua elements and velocity
covariance, so a change to the production file format or elements reaches the profiles; run_bench.sh records the md5 of
the file it generated, which shows such a change. The profile file has one `key value` per line, `#` starts a comment:

  iso oumuamua        the only ISO scripts/make_work_units.py knows
  stars N             stars in the file, the real ones included
  anchors K           the first K stars are the real ones (0, 1 or 2); the rest are made up
  star_seed S         seed of the made-up stars; the same seed in two profiles makes one the first stars of the other
  seed_base B         line 3 of `in`, mixed with each source_id into the clone seeds (src/nostos_io.c combine_seed)
  time_myr T          integration window, negative = backwards (default: make_work_units.py's)
  steps N             leapfrog steps (default: make_work_units.py's)
  threshold_pc X      a star is written to `out` when its closest clone comes under X pc
  apps cpu|gpu        which apps the profile is meant for; read by mac_ab.sh, ignored here

The made-up stars use integer arithmetic only (splitmix64 and micro-units), so no libm or Python version can change them:
ra 0..360 deg, dec -80..80 deg, parallax 1..20 mas, proper motions -80..80 mas/yr, radial velocity -100..100 km/s,
source_ids from 4e18 up. Every star takes 7 draws, so star k is the same in every profile with the same star_seed.
"""
import os
import sys

sys.dont_write_bytecode = True  # importing make_work_units must not leave a __pycache__ in scripts/

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import make_work_units as mwu  # noqa: E402

MASK = (1 << 64) - 1

# Gaia DR2 astrometry copied from tests/test_pipeline_oumuamua_hip3757.c (raw, no zero-point). Both stars come within a
# few pc of 'Oumuamua's path (the 6-star run in NOTES_method_selection.md gave 0.56 and 1.6 pc).
ANCHORS = [
    dict(source_id=2525688198820543360, ra=12.05634594268451, dec=-5.1359892365165924,
         parallax=40.20064272792413, pmra=182.16070593050424, pmdec=-120.15964959812227, rv=26.880687935465836),  # HIP 3757
    dict(source_id=3107000885484340224, ra=101.09513645394681, dec=-1.8534936165140703,
         parallax=24.124878291405853, pmra=-5.327852712260631, pmdec=-70.91032235521574, rv=33.124253095329394),  # HD 292249
]

KEYS = {"iso": str, "stars": int, "anchors": int, "star_seed": int, "seed_base": int, "time_myr": float, "steps": int,
        "threshold_pc": float, "apps": str}


class SplitMix64:
    """The PRNG of src/montecarlo.h (splitmix64), here only as a source of numbers that is the same on every machine."""

    def __init__(self, seed):
        self.state = seed & MASK

    def next(self):
        self.state = (self.state + 0x9E3779B97F4A7C15) & MASK
        z = self.state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK
        return z ^ (z >> 31)

    def below(self, n):
        return self.next() % n


def synthetic_star(rng, k):
    # source_ids rise with k (k * 1000003 + a draw below 1e6), so they are unique
    sid = 4000000000000000000 + k * 1000003 + rng.below(1000000)
    ra = rng.below(360000000)                    # micro-degrees
    dec = rng.below(160000000) - 80000000
    plx = 1000000 + rng.below(19000000)          # micro-mas
    pmra = rng.below(160000000) - 80000000       # micro-mas/yr
    pmdec = rng.below(160000000) - 80000000
    rv = rng.below(200000000) - 100000000        # micro-km/s
    return dict(source_id=sid, ra=ra / 1e6, dec=dec / 1e6, parallax=plx / 1e6, pmra=pmra / 1e6, pmdec=pmdec / 1e6,
                rv=rv / 1e6)


def read_profile(path):
    prof = {}
    with open(path) as f:
        for n, line in enumerate(f, 1):
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 2 or parts[0] not in KEYS:
                sys.exit("%s:%d: cannot read '%s' (keys: %s)" % (path, n, line, " ".join(sorted(KEYS))))
            try:
                prof[parts[0]] = KEYS[parts[0]](parts[1])
            except ValueError:
                sys.exit("%s:%d: bad value in '%s'" % (path, n, line))
    return prof


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 2
    prof = read_profile(argv[1])
    for key in ("iso", "stars", "star_seed", "threshold_pc"):
        if key not in prof:
            sys.exit("%s: needs '%s'" % (argv[1], key))
    if prof["iso"] != "oumuamua":
        sys.exit("%s: iso must be oumuamua (scripts/make_work_units.py has no other ISO)" % argv[1])
    n = prof["stars"]
    anchors = prof.get("anchors", 0)
    if not 0 <= anchors <= len(ANCHORS) or anchors > n or n < 1:
        sys.exit("%s: stars must be at least 1 and anchors between 0 and min(%d, stars)" % (argv[1], len(ANCHORS)))

    stars = list(ANCHORS[:anchors])
    rng = SplitMix64(prof["star_seed"])
    for k in range(n - anchors):
        stars.append(synthetic_star(rng, k))

    # write_work_unit() reads these module constants for lines 3 and 4 of the file
    mwu.DEFAULT_SEED_BASE = prof.get("seed_base", mwu.DEFAULT_SEED_BASE)
    mwu.DEFAULT_TOTAL_TIME_MYR = prof.get("time_myr", mwu.DEFAULT_TOTAL_TIME_MYR)
    mwu.DEFAULT_N_STEPS = prof.get("steps", mwu.DEFAULT_N_STEPS)
    mwu.DEFAULT_DISTANCE_THRESHOLD_PC = prof["threshold_pc"]
    mwu.write_work_unit(argv[2], stars)
    print("%s: %d stars (%d real), seed_base %d, %d steps over %g Myr, threshold %g pc"
          % (argv[2], n, anchors, mwu.DEFAULT_SEED_BASE, mwu.DEFAULT_N_STEPS, mwu.DEFAULT_TOTAL_TIME_MYR,
             mwu.DEFAULT_DISTANCE_THRESHOLD_PC), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
