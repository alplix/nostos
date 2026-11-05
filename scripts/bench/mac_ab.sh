#!/bin/bash
# One-command A/B on a Mac (Apple Silicon or Intel; it also runs on Linux): BASE_REF built the way the repository documents
# its CPU app (scripts/build_boinc_all.sh: -O2, -std=c99 modules), against this working tree, on the profiles in
# scripts/bench/profiles. Both sides run against the fake BOINC runtime of scripts/bench/boinc_stub.
# usage: scripts/bench/mac_ab.sh [rounds per profile, default 3] [profiles, default "quick short"]
#   env: BASE_REF  the baseline (default upstream/main if that remote branch exists, else origin/main)
#        NEW_OPT   compiler flags of the working-tree build (default -O2, the documented ones: with an unchanged tree
#                  the table is then an A/A comparison, which shows how much two identical builds differ here)
#        EXACT     1 (default here) prints 17 digits instead of 6 in `out`, so "same output" means the same doubles;
#                  0 builds both sides with the production format
#        EXPERIMENTS=1  also builds BASE_REF with -O3 and with -O2 -flto, the same-results flag candidates of
#                  docs/BENCHMARK.md (they stay "same output" only if the compiler fuses nothing differently)
#        NOSTOS_WORK (default ~/nostos-work): builds in build/ and bins/, results in bench-runs/<date>-<profile>/
# Needs the Xcode command line tools (or gcc, g++), python3 and git; no BOINC checkout. The CPU app only: there is no
# Metal app, and the CUDA and OpenCL apps need an NVIDIA or AMD GPU with fp64 (profile gpu-batch, build_variant.sh on
# Linux). The profiles with `apps gpu` are skipped. About 1 minute of building, then (quick + short, estimated from the
# 2.3-2.6 s per star of NOTES_method_selection.md, not measured) 2 labels x 3 rounds x about 2 minutes = 13 minutes. The
# production-sized profile is cpu-wu: `scripts/bench/mac_ab.sh 1 cpu-wu` takes about 2 hours.
#
# Binaries (labels):
#   cpu-base  BASE_REF, the documented build (-O2)
#   cpu-new   this working tree, NEW_OPT
#   cpu-O3, cpu-lto  with EXPERIMENTS=1: BASE_REF with -O3, with -O2 -flto
set -e -o pipefail
here="$(cd "$(dirname "$0")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
ROUNDS=${1:-3}
PROFILES=${2:-"quick short"}
W=${NOSTOS_WORK:-$HOME/nostos-work}
export NOSTOS_WORK=$W
export EXACT=${EXACT:-1}
if [ -z "${BASE_REF:-}" ]; then
  if git -C "$repo" rev-parse -q --verify upstream/main > /dev/null; then BASE_REF=upstream/main; else BASE_REF=origin/main; fi
fi
git -C "$repo" rev-parse -q --verify "$BASE_REF^{commit}" > /dev/null || { echo "no ref $BASE_REF (BASE_REF=...)"; exit 2; }
command -v python3 > /dev/null || { echo "python3 is needed"; exit 2; }
B=$W/bins

REF=$BASE_REF OPT="-O2" BACKENDS=cpu bash "$here/build_variant.sh" base
OPT="${NEW_OPT:--O2}" BACKENDS=cpu bash "$here/build_variant.sh" new
LABELS="cpu-base=$B/base/cpu cpu-new=$B/new/cpu"
if [ "${EXPERIMENTS:-0}" = 1 ]; then
  REF=$BASE_REF OPT="-O3" BACKENDS=cpu bash "$here/build_variant.sh" o3
  REF=$BASE_REF OPT="-O2 -flto" BACKENDS=cpu bash "$here/build_variant.sh" lto
  LABELS="$LABELS cpu-O3=$B/o3/cpu cpu-lto=$B/lto/cpu"
fi
set -- $LABELS
NLINES=$(( $# + 1 ))
for prof in $PROFILES; do
  apps=$(awk '$1 == "apps" {print $2}' "$here/profiles/$prof.wu" 2> /dev/null)
  if [ "${apps:-cpu}" != cpu ]; then
    echo "== $prof is for the GPU apps; skipped (nothing to run it with here)"
    continue
  fi
  echo "== $prof, $ROUNDS rounds, BASE_REF $BASE_REF, EXACT=$EXACT"
  bash "$here/run_bench.sh" "$prof" "$ROUNDS" $LABELS | tail -n "$NLINES"
done
