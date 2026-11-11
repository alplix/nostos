#!/bin/bash
# Linux A/B of the CPU app: BASE_REF built the documented way against this working tree, both with GCC (the compiler of
# the Linux release builds, scripts/build_boinc_all.sh), natively on Linux or in a Debian 12 container on a Mac (docker,
# or DOCKER=podman). On Apple Silicon the container is arm64 and runs at native speed, but it shares the Mac's cores and
# load: it gives ratios between Linux builds, not Linux-against-macOS numbers.
# usage: scripts/bench/linux_ab.sh [rounds, default 3] [profiles, default "quick short"]
#   env: NOSTOS_WORK (default ~/nostos-work; from a Mac the container works in $NOSTOS_WORK/linux), BASE_REF (default
#        upstream/main, else origin/main), EXACT (default 1: 17 digits, see build_variant.sh), CLANG=1 also builds both
#        sides with clang, EXPERIMENTS=1 also BASE_REF with -O3 and -O2 -flto, DOCKER (default docker)
# Labels: gcc-base, gcc-new (and clang-base, clang-new, gcc-O3, gcc-lto). Builds of one compiler must write the same
# `out`; GCC and clang may differ at rounding level (NOS-013).
# Needs on Linux: gcc, g++, git, python3, GNU time (and clang). On a Mac: docker or podman, network once for the image.
set -e -o pipefail
here="$(cd "$(dirname "$0")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
ROUNDS=${1:-3}
PROFILES=${2:-quick short}
if [ "$(uname)" = Darwin ]; then
  D=${DOCKER:-docker}
  W=${NOSTOS_WORK:-$HOME/nostos-work}/linux
  IMAGE=bitboinc-bench:linux-1      # bump the tag when the package list changes
  mkdir -p "$W/image"
  if ! $D image inspect $IMAGE > /dev/null 2>&1; then
    cat > "$W/image/Dockerfile" <<'EOF'
FROM python:3.12-slim-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends gcc g++ clang git patch time libc6-dev \
    libjemalloc-dev libmimalloc-dev && rm -rf /var/lib/apt/lists/* && git config --global --add safe.directory '*'
EOF
    $D build -q -t $IMAGE "$W/image"
  fi
  if [ -z "${BASE_REF:-}" ]; then
    if git -C "$repo" rev-parse -q --verify upstream/main > /dev/null; then BASE_REF=upstream/main; else BASE_REF=origin/main; fi
  fi
  exec $D run --rm -v "$repo":/src:ro -v "$W":/work -e NOSTOS_WORK=/work -e BASE_REF="$BASE_REF" -e EXACT="${EXACT:-1}" \
    -e CLANG="${CLANG:-0}" -e EXPERIMENTS="${EXPERIMENTS:-0}" $IMAGE bash /src/scripts/bench/linux_ab.sh "$ROUNDS" "$PROFILES"
fi

W=${NOSTOS_WORK:-$HOME/nostos-work}
export NOSTOS_WORK=$W EXACT=${EXACT:-1}
B=$W/bins
if [ -z "${BASE_REF:-}" ]; then
  if git -C "$repo" rev-parse -q --verify upstream/main > /dev/null; then BASE_REF=upstream/main; else BASE_REF=origin/main; fi
fi
echo "== $(gcc --version | head -1), $(ldd --version 2>&1 | head -1), $(uname -m)"
REF=$BASE_REF CC=gcc CXX=g++ BACKENDS=cpu bash "$here/build_variant.sh" lx-base
CC=gcc CXX=g++ BACKENDS=cpu bash "$here/build_variant.sh" lx-new
LABELS="gcc-base=$B/lx-base/cpu gcc-new=$B/lx-new/cpu"
if [ "${EXPERIMENTS:-0}" = 1 ]; then
  REF=$BASE_REF CC=gcc CXX=g++ OPT="-O3" BACKENDS=cpu bash "$here/build_variant.sh" lx-o3
  REF=$BASE_REF CC=gcc CXX=g++ OPT="-O2 -flto" BACKENDS=cpu bash "$here/build_variant.sh" lx-lto
  LABELS="$LABELS gcc-O3=$B/lx-o3/cpu gcc-lto=$B/lx-lto/cpu"
fi
if [ "${CLANG:-0}" = 1 ]; then
  echo "== $(clang --version | head -1)"
  REF=$BASE_REF CC=clang CXX=clang++ BACKENDS=cpu bash "$here/build_variant.sh" lx-clang-base
  CC=clang CXX=clang++ BACKENDS=cpu bash "$here/build_variant.sh" lx-clang-new
  LABELS="$LABELS clang-base=$B/lx-clang-base/cpu clang-new=$B/lx-clang-new/cpu"
fi
for prof in $PROFILES; do
  echo "== $prof, $ROUNDS rounds, BASE_REF $BASE_REF, EXACT=$EXACT"
  bash "$here/run_bench.sh" "$prof" "$ROUNDS" $LABELS
done
