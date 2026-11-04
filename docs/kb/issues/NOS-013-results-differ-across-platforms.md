---
id: NOS-013
title: Results differ across platforms in the last digits; does the quorum compare tolerate that and threshold-edge rows?
type: question
status: new
severity: medium
platforms: [aarch64, android, windows, clang, cuda, opencl]
area: determinism
next: Decide: does the quorum compare use a numeric tolerance and accept a row missing at the threshold? Then diff ARM64/Android builds
pattern: PAT-018
aliases: []
reported: 2026-10-27
source: audit 2026-10-27; commit messages of 1ecc658 and dea832d; question 3 of the 2026-10-21 plan (validator tolerance)
where: scripts/build_boinc_all.sh:23, 43, scripts/build_boinc_android.sh:11-17, src/nostos_cuda.cu:1-2, src/main_nostos_boinc.cpp:135 @ eb7029e
fixed_in:
---

## Symptom

Measured by the maintainer (commit message of 1ecc658): on a 3-star test work unit "every build except Windows ARM64
and Android gave the same output, to the last printed digit (native, qemu, wine, Mac arm64 and x86_64)". The GPU apps
gave the same output as the Linux CPU app on 3 and 500 stars (dea832d).

## Cause

Partly verified by reading the build scripts and sources at eb7029e, and by counting FMA instructions with the two
compilers available here (below). The two differing compilers themselves were not run.

* FMA contraction. Both differing builds are compiled with clang for aarch64: llvm-mingw for Windows ARM64
  (`scripts/build_boinc_all.sh:43`) and the Android NDK r27c (`scripts/build_boinc_android.sh:11-17`). Upstream clang
  fuses `a*b+c` within one expression by default (`-ffp-contract=on`; Apple clang does not, see below), and the orbit code is full of such expressions
  (`src/potential.c:63, 71, 80`; `src/montecarlo.c:52`, `sum += L[i][j] * z[j]`). GCC with `-std=c99`, which the
  other Linux and Windows builds use (`scripts/build_boinc_all.sh:23`), does not contract. So the comment in
  `src/nostos_cuda.cu:1-2`, "the CPU builds do not contract a*b+c into FMAs", is not true for these two builds.
  Measured on 2026-10-27 on macOS arm64 at eb7029e (`-O2 -S`, counting `fmadd`/`fmsub`/`fnmadd`/`fnmsub` in
  `potential.c`, `montecarlo.c`, `coords.c`, `encounter.c`, `orbital_elements.c`): GCC 13.2 with `-std=c99` emits 0
  in every file (with `-std=gnu99`: 16, 2, 41, 10, 19); Apple clang 14 emits 0 by default and 15, 2, 44, 4, 48 with
  `-ffp-contract=on`, the upstream clang default. So the code has many contractible expressions, the GCC builds do
  not fuse them, and Apple clang does not either, which fits the macOS arm64 build matching (its build script is not
  in the repository, NOS-018). The NDK r27c and llvm-mingw clangs were not available here; that they fuse is the
  likely explanation, still to measure. The C++ wrappers, which g++ compiles with contraction on, contain no
  `a*b+c` on result values.
* libm. `pow` (`src/potential.c:57, 71`), `log` and `cos` (`src/montecarlo.c:15-16`, the clone draws) and the
  trigonometry of `src/coords.c` come from glibc, the Windows C runtime, bionic, Apple's libm, or the CUDA and OpenCL
  built-ins; their last bits can differ.
* Thresholds. A row is written only when the closest clone is under the threshold (`src/main_nostos_boinc.cpp:135`,
  `src/main_nostos_gpu_boinc.cpp:168`), with six decimals in pc. A last-bit difference can change the last printed
  digit, or add or drop the whole row of a star at the threshold.

See [PAT-018](../shared/patterns/PAT-018-floating-point-differs-across-platforms.md). The quorum compare
(`server/nostos_compare.sh`, named in `src/nostos_gpu.h:10`) is not in this repository.

Questions for the maintainers: does the compare match lines exactly or numbers with a tolerance, and what does it do
when one result has a row that the other lacks? Are CPU and GPU versions, and the Android and Windows ARM64 builds,
mixed in one quorum?

## Fix

Measure first: run one work unit on the Linux x86_64 build and on the Android and Windows ARM64 builds and diff; then
rebuild those two with `-ffp-contract=off` and diff again. If that removes the difference, add `-ffp-contract=off` to
every C compile and correct the comment in `src/nostos_cuda.cu`. A numeric tolerance in the compare covers the libm
part; a row present in only one result can be accepted when the other side's star is within the tolerance of the
threshold. Do not change the science results silently: record the measured effect here.

## Test

A cross-platform diff over a few small work units (native, qemu, Wine), run whenever compiler flags or math code
change.

## Log

- 2026-10-27: recorded as a question; build flags and code read at eb7029e, the measurement taken from the 1ecc658
  commit message.
- 2026-10-27: re-checked after the pattern was corrected (count FMA instructions instead of reasoning from compiler
  versions): counted with GCC 13.2 and Apple clang 14 on arm64; the Apple clang result explains why the macOS build
  matched. Status unchanged.
