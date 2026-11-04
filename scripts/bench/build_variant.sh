#!/bin/bash
# Builds one variant of the Nostos BOINC apps against the fake BOINC runtime (scripts/bench/boinc_stub), for run_bench.sh.
# usage: build_variant.sh <name> [patch ...]
#   Copies the repository's tracked files (or, with REF, a git ref) to $NOSTOS_WORK/build/<name>, applies the patches (git
#   diff format, -p1), compiles the way scripts/build_boinc_all.sh (CPU) and scripts/build_gpu_server.sh (CUDA, OpenCL) do,
#   and writes bins/<name>/cpu (and cuda, opencl) with a <backend>.info file that records how it was built.
#   env: OPT       flags of every compile and of the link (default -O2, the documented flag), e.g. "-O3", "-O2 -flto",
#                  "-O2 -ffp-contract=off"; the C modules always get -std=c99 as well
#        CC, CXX   compilers (default clang and clang++ on macOS, gcc and g++ elsewhere)
#        BACKENDS  default "cpu", plus "cuda" when nvcc is in PATH on Linux, plus "opencl" on Linux when CL/cl.h is found
#                  (in /usr/include or in $OPENCL_INC). There is no Metal app, and the OpenCL app needs fp64, which Apple's
#                  GPUs lack, so a Mac builds the CPU app only; naming a backend that cannot be built is an error.
#        EXACT  1 prints the three numbers of every row of `out` with %.17g instead of %.6f (the copy's sources are edited,
#               nothing else changes), so that the hash of `out` covers every digit of the doubles. Build both sides with
#               it when a change touches arithmetic; the production format is the default.
#        ALLOC  system (default), jemalloc (JEMALLOC = static library or linker flags; default Homebrew's libjemalloc.a on
#               macOS, -ljemalloc elsewhere) or mimalloc (MIMALLOC, default -lmimalloc). The hot loops allocate nothing,
#               so no gain is expected; the knob is there to show that.
#        STATIC 1 links the CPU app with -static like the Linux release builds (default 1 on Linux, 0 on macOS, where it
#               cannot be done); 0 if the static libraries are missing
#        CXXSTD language level of the CPU wrapper (default none: the compiler's own, as build_boinc_all.sh does)
#        NVCCFLAGS (default "-O3 -std=c++17 --fmad=false", build_gpu_server.sh), CUDA_GENCODE (default -arch=native; the
#               release builds list sm_50 to sm_120), OCL_VENDOR (default none: any GPU; releases use amd or nvidia)
#        REF    build a git ref of this repository (e.g. upstream/main) instead of its working tree
#        NOSTOS_WORK (default ~/nostos-work), NOSTOS_SRC (the repository to copy, default this one)
#   Only tracked files are copied: `git add -N` a new file of the candidate, or it is not built.
#   The fake BOINC runtime always comes from this script's own directory, so a REF that predates scripts/bench builds too.
set -e -o pipefail
here="$(cd "$(dirname "$0")" && pwd)"
repo=${NOSTOS_SRC:-$(cd "$here/../.." && pwd)}
[ $# -ge 1 ] || { awk 'NR>1 && /^#/ {print; next} NR>1 {exit}' "$0"; exit 2; }
name=$1; shift
W=${NOSTOS_WORK:-$HOME/nostos-work}
OPT=${OPT:--O2}; ALLOC=${ALLOC:-system}; EXACT=${EXACT:-0}; CXXSTD=${CXXSTD:-}
OS=$(uname)
if [ "$OS" = Darwin ]; then
  CC=${CC:-clang}; CXX=${CXX:-clang++}; STATIC=${STATIC:-0}
else
  CC=${CC:-gcc}; CXX=${CXX:-g++}; STATIC=${STATIC:-1}
fi
if [ -z "${BACKENDS:-}" ]; then
  BACKENDS=cpu
  if [ "$OS" = Linux ]; then
    if command -v nvcc > /dev/null 2>&1; then BACKENDS="$BACKENDS cuda"; fi
    if [ -f /usr/include/CL/cl.h ] || [ -f "${OPENCL_INC:-/nonexistent}/CL/cl.h" ]; then BACKENDS="$BACKENDS opencl"; fi
  fi
fi
T=$W/build/$name; out=$W/bins/$name
rm -rf "$T"; mkdir -p "$T" "$out"
log=$T/build.log
if [ -n "${REF:-}" ]; then
  git -C "$repo" rev-parse -q --verify "$REF^{commit}" > /dev/null || { echo "no such ref in $repo: $REF"; exit 2; }
  (cd "$repo" && git archive "$REF") | tar -xf - -C "$T"
  source_desc="$REF ($(git -C "$repo" rev-parse --short "$REF"))"
else
  (cd "$repo" && git ls-files -z | tar --null -T - -cf -) | tar -xf - -C "$T"
  source_desc="working tree at $(git -C "$repo" rev-parse --short HEAD)"
  if [ -n "$(git -C "$repo" status --porcelain --untracked-files=no | head -1)" ]; then source_desc="$source_desc plus uncommitted changes"; fi
fi
for p in "$@"; do
  [ -f "$p" ] || { echo "no patch file $p"; exit 2; }
  a="$(cd "$(dirname "$p")" && pwd)/$(basename "$p")"
  (cd "$T" && patch -p1 -s < "$a") || { echo "$p does not apply to $source_desc"; exit 1; }
done
[ -f "$T/src/main_nostos_boinc.cpp" ] || { echo "the copy has no src/main_nostos_boinc.cpp: wrong ref or repository?"; exit 1; }
# the output format of both wrappers, in a form the hash of out can tell digit for digit; sed -i differs between BSD and
# GNU, so write a new file and move it
if [ "$EXACT" = 1 ]; then
  for f in main_nostos_boinc.cpp main_nostos_gpu_boinc.cpp; do
    sed 's/%\.6f/%.17g/g' "$T/src/$f" > "$T/src/$f.exact" && mv "$T/src/$f.exact" "$T/src/$f"
    grep -q '%\.17g' "$T/src/$f" || { echo "EXACT=1: no %.6f in src/$f, has the output format changed?"; exit 1; }
  done
fi
case $ALLOC in
  system) L="" ;;
  jemalloc)
    if [ "$OS" = Darwin ]; then L=${JEMALLOC:-$(brew --prefix 2> /dev/null || echo /opt/homebrew)/lib/libjemalloc.a}
    else L=${JEMALLOC:--ljemalloc}; fi ;;
  mimalloc) L=${MIMALLOC:--lmimalloc} ;;
  *) echo "ALLOC must be system, jemalloc or mimalloc"; exit 2 ;;
esac
case $L in /*) [ -f "$L" ] || { echo "no $L (set JEMALLOC or MIMALLOC)"; exit 2; } ;; esac

# one compile or link step: on failure the error lines of the log are shown and the script stops
step() {
  echo "+ $*" >> "$log"
  "$@" >> "$log" 2>&1 || { grep -E 'error|undefined|ld:|cannot find' "$log" | tail -8; echo "FAILED: $* (see $log)"; exit 1; }
}
# a step that must have written its output
have() { [ -s "$1" ] || { echo "$1 was not written (see $log)"; exit 1; }; }

obj=$T/obj; mkdir -p "$obj"
CORE="coords orbital_elements potential encounter montecarlo mc_search nostos_io"   # build_boinc_all.sh's CSRC
for m in $CORE; do
  step $CC $OPT -std=c99 -I "$T/src" -c "$T/src/$m.c" -o "$obj/$m.o"
  have "$obj/$m.o"
done
step $CXX $OPT -I "$here/boinc_stub" -c "$here/boinc_stub/boinc_stub.cpp" -o "$obj/boinc_stub.o"
have "$obj/boinc_stub.o"

# the GPU wrappers need only these (build_gpu_server.sh's CSRC)
gpu_objs=()
for m in coords orbital_elements potential montecarlo nostos_io; do gpu_objs+=("$obj/$m.o"); done

for b in $BACKENDS; do
  rm -f "$out/$b" "$out/$b.info"
  case $b in
    cpu)
      cpu_objs=()
      for m in $CORE; do cpu_objs+=("$obj/$m.o"); done
      step $CXX $OPT $CXXSTD -I "$T/src" -I "$here/boinc_stub" -c "$T/src/main_nostos_boinc.cpp" -o "$obj/main_cpu.o"
      have "$obj/main_cpu.o"
      link_static=""; if [ "$STATIC" = 1 ]; then link_static=-static; fi
      step $CXX $OPT $link_static "$obj/main_cpu.o" "$obj/boinc_stub.o" "${cpu_objs[@]}" $L -lpthread -lm -o "$out/cpu" ;;
    cuda)
      command -v nvcc > /dev/null 2>&1 || { echo "no nvcc: the CUDA app cannot be built here"; exit 1; }
      step nvcc ${NVCCFLAGS:--O3 -std=c++17 --fmad=false} ${CUDA_GENCODE:--arch=native} -I "$T/src" -c "$T/src/nostos_cuda.cu" -o "$obj/nostos_cuda.o"
      have "$obj/nostos_cuda.o"
      step $CXX $OPT -std=c++17 -DNG_CUDA -I "$T/src" -I "$here/boinc_stub" -c "$T/src/main_nostos_gpu_boinc.cpp" -o "$obj/main_cuda.o"
      have "$obj/main_cuda.o"
      # nvcc links, so that cudart and the host libraries are found without knowing the toolkit's directories
      step nvcc ${CUDA_GENCODE:--arch=native} "$obj/main_cuda.o" "$obj/nostos_cuda.o" "$obj/boinc_stub.o" "${gpu_objs[@]}" -lpthread -o "$out/cuda" ;;
    opencl)
      cl_inc=""; if [ -n "${OPENCL_INC:-}" ]; then cl_inc="-I $OPENCL_INC"; fi
      posix=""; if [ "$OS" = Linux ]; then posix=-D_POSIX_C_SOURCE=200809L; fi    # clock_gettime under -std=c99 (291f413)
      vendor=""; if [ -n "${OCL_VENDOR:-}" ]; then vendor="-DEXPECTED_GPU_VENDOR=\"$OCL_VENDOR\""; fi
      step $CC $OPT -std=c99 $posix $cl_inc -DCL_TARGET_OPENCL_VERSION=120 -I "$T/src" -c "$T/src/cl_dynload.c" -o "$obj/cl_dynload.o"
      step $CC $OPT -std=c99 $posix $cl_inc -DCL_TARGET_OPENCL_VERSION=120 $vendor -I "$T/src" -c "$T/src/nostos_opencl.c" -o "$obj/nostos_opencl.o"
      step $CXX $OPT -std=c++17 -DNG_OPENCL -I "$T/src" -I "$here/boinc_stub" -c "$T/src/main_nostos_gpu_boinc.cpp" -o "$obj/main_opencl.o"
      have "$obj/main_opencl.o"
      libdl=""; if [ "$OS" = Linux ]; then libdl=-ldl; fi
      step $CXX $OPT "$obj/main_opencl.o" "$obj/nostos_opencl.o" "$obj/cl_dynload.o" "$obj/boinc_stub.o" "${gpu_objs[@]}" \
        $libdl -lpthread -lm -o "$out/opencl" ;;
    *) echo "unknown backend $b (cpu, cuda, opencl)"; exit 2 ;;
  esac
  [ -x "$out/$b" ] || { echo "build of $b failed: no executable (see $log)"; exit 1; }
  {
    echo "name $name, backend $b, built $(date '+%Y-%m-%d %H:%M:%S') on $(uname -sm)"
    echo "source $source_desc; patches: ${*:-none}"
    echo "OPT=$OPT EXACT=$EXACT ALLOC=$ALLOC STATIC=$STATIC CXXSTD=${CXXSTD:-default}"
    echo "CC $($CC --version 2>&1 | head -1)"
    echo "CXX $($CXX --version 2>&1 | head -1)"
  } > "$out/$b.info"
done
echo "built $out: $(cd "$out" && ls | grep -v '\.info$' | tr '\n' ' ')($source_desc, OPT=$OPT EXACT=$EXACT ALLOC=$ALLOC)"
