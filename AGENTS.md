# AGENTS.md: nostos

Orientation for coding agents and new contributors. The science, the method and its validation are in
[README.md](README.md) and [NOTES_method_selection.md](NOTES_method_selection.md); this file covers the code. Line
numbers refer to commit eb7029e.

## What one work unit is

One work unit is one interstellar object (ISO) against one chunk of Gaia DR3 stars. The input file `in` is plain text:
line 1 the ISO's barycentric orbital elements (`e q_au i_deg node_deg peri_deg`, epoch 1600, far out on the inbound
leg), line 2 the upper triangle of its 3x3 velocity covariance in (kpc/Myr)^2, line 3 a 64-bit `seed_base`, line 4
`total_time_myr n_steps distance_threshold_pc`, then one star per line: `source_id ra_deg dec_deg parallax_mas
pmra_masyr pmdec_masyr rv_kms` (parallax already zero-point corrected). The app starts the ISO on its inbound
asymptote, and for every star draws 100 ISO velocity clones (seed `combine_seed(seed_base, source_id)`), integrates
star and clone backwards through the Allen & Santillan (1991) Milky Way potential with a leapfrog (defaults in
`scripts/make_work_units.py`: -100 Myr, 100,000 steps, 50 pc), and keeps each clone's closest approach. The output
file `out` has one line `source_id min_pc mean_pc t_myr` (closest clone, mean over clones, time of the closest clone)
for every star whose closest clone came under the threshold, and ends with `# done <n_stars>`. The server's generator
(not here) makes CPU work units of 1500 stars (about 2.5 s per star, so about an hour) and GPU work units of 20,000
stars (commit messages of 1ecc658 and dea832d).

## Layout

| Path | What it holds |
|---|---|
| `src/main_nostos_boinc.cpp` | CPU BOINC app: read `in`, `mc_search_stage1()` per star, checkpoint, write `out` |
| `src/main_nostos_gpu_boinc.cpp` | GPU BOINC app (`-DNG_CUDA` or `-DNG_OPENCL`): same files; clone draws and per-star reduction on the host, integration on the GPU in batches of 2000 stars |
| `src/nostos_gpu.h` | GPU backend interface; adaptive kernel slices (`ng_next_slice`, about 0.25 s per launch) |
| `src/nostos_cuda.cu`, `src/nostos_opencl.c` | CUDA and OpenCL backends (device choice, fp64 check, launches) |
| `src/nostos_kernel_body.h` | device code of both backends: operation-by-operation copy of `milkyway_accel()` and `leapfrog_step()` |
| `src/nostos_kernel_cl.inc` | generated from `nostos_kernel_body.h` by `scripts/embed_kernel.py`; do not edit |
| `src/cl_dynload.c`, `.h` | run-time OpenCL loader and vendor matching (copy, see Invariants) |
| `src/potential.c`, `.h`, `src/integrator.h`, `src/vec3.h`, `src/units.h` | orbit core: Allen & Santillan potential, leapfrog, vectors, units (kpc, Myr, kpc/Myr) (copies) |
| `src/coords.c`, `.h` | ICRS astrometry to Galactocentric phase space, solar parameters (copy) |
| `src/montecarlo.c`, `.h` | seedable RNG, Box-Muller, Cholesky (copy) |
| `src/orbital_elements.c`, `.h`, `src/iso_solutions.h` | hyperbolic elements to Galactocentric state, inbound asymptote; published solutions of 1I, 2I, 3I |
| `src/encounter.c`, `.h`, `src/mc_search.c`, `.h` | closest approach of one star and one ISO; stage 1 (100 ISO clones per star) |
| `src/nostos_io.c`, `.h` | work-unit parser `parse_input_stream()` and `combine_seed()` |
| `src/star_mc.c`, `.h` | stage 2: star-side resampling from Gaia errors (used by `tools/`) |
| `src/gaia_scan.c`, `src/scan_logic.c`, `src/csv_range.c`, `src/parallax_zeropoint.c` (and `.h`) | standalone full-catalogue scan tool, not the BOINC app: straight-line screen, gravity confirmation, byte-range split, DR3 parallax zero-point |
| `src/nostos_detail.c` | server-side follow-up: re-runs stage 1 of the found stars with the same seeds, writes percentiles and paths (built on the server; no Makefile target) |
| `tools/` | stage-2 MC of a shortlist (`mc_shortlist.c`) and two diagnostics; `make tools` |
| `tests/` | C tests per module and against published encounters, `test_make_work_units.py`; `make test` |
| `scripts/build_*.sh`, `scripts/build_*.bat` | builds, see Building |
| `scripts/list_gaia_bulk.py`, `scripts/download_filter_6d.py` | list and download the Gaia DR3 bulk files, keep rows with a radial velocity |
| `scripts/make_work_units.py`, `scripts/make_star_catalog.py` | work-unit input files; the compact star catalogue the server cuts work units from (both apply the parallax zero-point) |
| `scripts/embed_kernel.py` | regenerates `src/nostos_kernel_cl.inc` |
| other `scripts/` | scan post-processing (`summarize_candidates.py`, `make_shortlist.py`, `prepare_mc_shortlist.py`, `rank_mc_shortlist.py`, `scan_worker.sh`) and cross-checks against astropy, galpy and the official zero-point package |
| `docs/kb/` | knowledge base (below) |

There are no vendored third-party trees.

## Invariants

Break one of these and results become silently wrong or stop matching the server.

* **Input format.** Written by `scripts/make_work_units.py` and by the server's generator (header plus a slice of the
  `make_star_catalog.py` catalogue); parsed only by `parse_input_stream()` (`src/nostos_io.c`), which both apps and
  `src/nostos_detail.c` use. The parser stops at the first star line that does not parse; `# done N` counts what was
  read.
* **Output format.** Both apps write `"%llu %.6f %.6f %.6f\n"` and `"# done %zu\n"`
  (`src/main_nostos_boinc.cpp:136-140, 159`, `src/main_nostos_gpu_boinc.cpp:169-170, 179`). The server's validator
  checks the trailer, the star count and value ranges, and a compare script checks the two results of a quorum
  (commit messages of 1ecc658 and dea832d). Change the format only together with the server.
* **Seeds.** Each star's RNG is seeded with `combine_seed(seed_base, source_id)` (`src/nostos_io.c:65-70`) and the
  clones are drawn with `mc_sample_multivariate_gaussian()` in clone order. Changing the mixing, the RNG, the clone
  count (`MC_N_ISO_CLONES`) or the order changes every result, breaks quorums with older app versions and the
  server's `nostos_detail` re-run, whose distances must equal the volunteers' (728ea35).
* **CPU and GPU give the same numbers.** `src/nostos_kernel_body.h` repeats `milkyway_accel()` and `leapfrog_step()`
  operation by operation, and `src/main_nostos_gpu_boinc.cpp:104-172` repeats `mc_search_stage1()` (Cholesky, draws,
  minimum and mean in clone order). Change both sides together, then run `python scripts/embed_kernel.py`; the
  potential's parameters go to the kernel in the order of the `NG_*` enum (`src/nostos_gpu.h:19`). At eb7029e the
  generated `.inc` matches the body (checked).
* **Floating point.** C sources are compiled with `-std=c99` (GCC then does not fuse multiply-adds), CUDA with
  `--fmad=false`, the OpenCL kernel with `#pragma OPENCL FP_CONTRACT OFF`, 32-bit x86 with `-msse2 -mfpmath=sse`.
  GCC with `-std=c99` and Apple clang 14 emit no fused multiply-adds for this code (counted on 2026-10-27); the
  Android and Windows ARM64 builds (upstream clang) gave different output, probably because they fuse (NOS-013).
* **ISO elements.** The work unit must carry barycentric elements of the inbound leg (`src/iso_solutions.h`,
  `*_BARY_1600`); `scripts/make_work_units.py` must hold the same numbers (`tests/test_make_work_units.py`).
* **Copies.** `src/coords.*`, `src/potential.*`, `src/montecarlo.*`, `src/integrator.h`, `src/vec3.h`, `src/units.h`
  and `src/cl_dynload.*` are byte-identical to phygas's (`cl_dynload.c` differs from planodia's, NOS-016). They are
  listed in `docs/kb/kb.conf`; change all copies together (NOS-020).

## Checkpoint and resume

* The file is `nostos_checkpoint` in the slot directory: one decimal number, the index of the next star. It is the
  only authority; `out` is reopened with `"a"` when that number is above 0 and created with `"w"` when it is 0, missing
  or unreadable.
* CPU app: after each star, when `boinc_time_to_checkpoint()` returns 1, the file is rewritten in place inside that
  critical section, then `boinc_checkpoint_completed()` (`src/main_nostos_boinc.cpp:143-150`). `out` is not flushed
  first. Progress is reported after every star.
* GPU app: after each batch (`NOSTOS_BATCH_STARS`, default 2000): the batch's rows, `fflush(out)`, the file rewritten
  in place, `boinc_checkpoint_completed()`, with no timer and no critical section (`src/main_nostos_gpu_boinc.cpp:161-177`).
  Progress is reported after every kernel slice.
* A resume must reproduce the rows of the remaining stars byte for byte; per-star seeds make a recomputed star give
  the same row on the same platform.
* Known defects at eb7029e:
  * Rows written after the last checkpoint are kept and written again, and the CPU app can count rows that never
    reached the disk (NOS-001).
  * A restart after the `# done` trailer was written but before `boinc_finish()` completed appends a second trailer
    (NOS-001).
  * The checkpoint is rewritten in place, its write is not checked, and completion is reported even when it failed
    (NOS-002).
  * The GPU app writes the checkpoint after every batch, without `boinc_time_to_checkpoint()` (NOS-003).
  * A checkpoint that cannot be opened restarts the task from star 0 and truncates `out` (NOS-010).
  * Write errors on `out` are not checked before `boinc_finish(0)` (NOS-011).
  * Windows builds write `out` with CRLF line ends (NOS-014).

## Building

Tests and tools (any Unix with a C compiler; `gcc` may be clang):

```
make test                    # builds build/test_*, runs them, then tests/test_make_work_units.py
make test PYTHON=python3     # where there is no `python`
make tools                   # build/gaia_scan, build/mc_shortlist, build/mc_vs_bestfit_ross573, build/oumuamua_velocity_covariance
make clean
```

CPU BOINC app, on the bitboinc build server, from a directory that contains `src/`:
`bash scripts/build_boinc_all.sh` (cross compilers, BOINC libraries per target in `/opt/bitboinc-src/boinc-<target>`,
llvm-mingw in `/opt/llvm-mingw`). Output `out/nostos_v${NOSTOS_VERSION:-1}_<os>_<arch>[.exe]` for Linux x86_64,
aarch64, armhf, i686, ppc64le, riscv64 and Windows x64, x86, ARM64. Per target it runs:

```
$CC  -O2 -std=c99 -I src -c src/<m>.c     # coords orbital_elements potential encounter montecarlo mc_search nostos_io
$CXX -O2 -I src -I $B/api -I $B/lib -c src/main_nostos_boinc.cpp
$CXX -static obj/*.o $B/api/libboinc_api.a $B/lib/libboinc.a -lpthread -lm     # Windows: -lpsapi
```

Linux builds are static; the script prints the number of `GLIBC_` symbols (expected 0). Other CPU builds:
`scripts/build_boinc_android.sh` (aarch64, WSL as root, NDK r27c, API 21, BOINC built with the NDK,
`-static-libstdc++`), `scripts/build_boinc_linux.sh` (local x86_64 static build in WSL against BOINC 8.2.8,
`BOINC_SRC`, `BOINC_BUILD`). Both have fixed developer paths (NOS-015). The macOS script named in
`build_boinc_all.sh` is not in the repository (NOS-018).

GPU apps: `bash scripts/build_gpu_server.sh` on the build server (sudo Docker; BOINC source in
`/opt/bitboinc-src/boinc`, built inside each container): Linux CUDA in `nvidia/cuda:12.8.1-devel-ubuntu22.04`
(`nvcc -O3 -std=c++17 --fmad=false`, sm_50 to sm_120, static cudart), Linux OpenCL in `ubuntu:22.04`, Windows OpenCL
with the host's MinGW-w64; one OpenCL binary per vendor (`-DEXPECTED_GPU_VENDOR="amd"` tagged `ati`, `"nvidia"`).
The Linux binaries need glibc 2.34 (NOS-005). Windows CUDA: `scripts\build_gpu_windows.bat` (MSVC 2022 Build Tools,
CUDA 13.3, sm_75 and newer, NOS-021); `scripts\build_gpu_opencl_wintest.bat` is a local MSVC test build. After editing
`src/nostos_kernel_body.h`: `python scripts/embed_kernel.py`.

`NOSTOS_VERSION` (default 1) only sets the file names. GPU run-time knobs: `NOSTOS_BATCH_STARS` (2000),
`NOSTOS_STEPS_PER_LAUNCH` (upper bound of a slice, 5000). `--device` applies only when `init_data.xml` assigns no
device; `--cl-platform` (default 0) is tried first and kept only if it has GPUs of the build's vendor.

Data (developer machine): `python scripts/list_gaia_bulk.py`, `python scripts/download_filter_6d.py` (fixed
`D:/nostos_data`, NOS-009), then `python scripts/make_star_catalog.py <gaia_csv> <out.txt.gz>` or
`python scripts/make_work_units.py <gaia_csv> <out_dir> [--chunk-size N] [--zpt lindegren|global|none]`. The download
marks a bulk file done after five failed attempts, so its rows are lost for good (NOS-006), and writes a file's rows
twice when it stops between writing them and marking the file done (NOS-008).

## Tests

`make test` runs 17 C tests (orbit integrator against Kepler, the potential against galpy and its own consistency,
coordinates against astropy and JPL Horizons, Monte Carlo statistics, the parser, the scan logic, the zero-point port,
and three published encounters of 'Oumuamua and Borisov) and `tests/test_make_work_units.py`. Run on 2026-10-27 at
eb7029e on macOS arm64 (Apple clang 14 as `gcc`, Python 3.11) as `make test PYTHON=python3` in a scratch copy: all
passed in about 13 s; the Lindegren zero-point case was skipped because `gaiadr3_zeropoint` was not installed.
`make tools` also built.

Not covered: the BOINC wrappers (no BOINC stub, no kill/resume harness), the GPU backends (CPU/GPU agreement was
checked by hand, dea832d and 291f413), the data scripts and the build scripts.

What to run for which change:

* potential, integrator, encounter: `make test`; mirror the change in `nostos_kernel_body.h`, run `embed_kernel.py`,
  and compare one work unit between the CPU and a GPU build (`cmp` of `out`).
* coords, orbital elements, `iso_solutions.h`: `make test` (coordinate, Horizons and pipeline tests,
  `test_make_work_units.py`).
* montecarlo, mc_search, nostos_io: `make test`; any change of numbers also affects `nostos_detail.c` and quorums with
  older versions.
* `main_*_boinc.cpp`: build against the BOINC libraries, run a small work unit standalone (a directory with the input
  named `in`; the app writes `out` and `boinc_finish_called`, and `nostos_checkpoint` once one was taken) and `cmp`
  `out` with the previous build's. For checkpoint changes, stop and restart with a quit that flushes (a BOINC client, or a stub that calls
  `fflush(NULL)`), not with SIGKILL (PAT-002).
* shared copies: the same change in phygas (and planodia for `cl_dynload`).

## Environment notes

* macOS: `make test` works with Apple clang; pass `PYTHON=python3` if `python` is missing.
* Release builds need the bitboinc build server (cross compilers, per-target BOINC libraries, Docker) or a Windows PC
  with WSL and MSVC; the scripts cannot run elsewhere without edits (NOS-015).
* The kernel needs fp64: no Metal, Intel or Adreno builds; OpenCL plan classes `opencl_ati_fp64` and
  `opencl_nvidia_fp64` (dea832d); the app refuses a device without fp64 (`src/nostos_opencl.c:52-55`).
* GPU errors: initialisation and kernel failures end the task with `boinc_finish(1)`, a computation error
  (NOS-004, not confirmed yet), and the OpenCL launch loop ignores the results of `clSetKernelArg` (NOS-022).
* Linux CUDA uses CUDA 12.8.1 (sm_50 and newer), Windows CUDA uses CUDA 13.3 (sm_75 and newer).
* Linux OpenCL needs `-D_POSIX_C_SOURCE=200809L` for `clock_gettime` (291f413).
* Build inputs are not pinned: container images by tag, unversioned `apt-get` packages, BOINC from the build server
  without a printed version (NOS-023).
* MSVC's C++ rejects compound literals, so `src/vec3.h` uses `v3_make()`; the C headers need `extern "C"` in the C++
  wrappers (the 2026-10-08 link failure).
* The Gaia CSV is 43 GB. Python packages: `requests` (download), `numpy` and `gaiadr3-zeropoint` (default `--zpt
  lindegren`), `astropy`, `galpy`, `astroquery` (cross-checks).

## Conventions for changes

* Small branches, one fix per commit, each with the test that fails before it where possible.
* Commit messages name the issue (`KB: NOS-NNN`) and the pattern if any (`Pattern: PAT-NNN`).
* Measured claims only: say what was run, on which platform, at which commit. Date comments that record a finding
  (the code already does: "2026-10-09: ...").
* Keep the CPU app, the GPU host code and both GPU backends symmetric; a change to the orbit code goes to the kernel
  body and to every shared copy.
* Never change science results silently: record the measured effect in the issue and in NOTES_method_selection.md.

## Not in this repository

* The server side named in the commit messages (a `server/` directory elsewhere): database schema (ISOs, encounters,
  details), work generator, validator, quorum compare (`server/nostos_compare.sh`), assimilator, the job that runs
  `nostos_detail`, BOINC templates, plan classes, and the web pages (`nostos.php`, `nostos_3d.php`, admin and export
  pages).
* The data: the Gaia CSV, `progress.txt`, scan results (under `D:/nostos_data` on the developer machine).
* The BOINC libraries and the macOS build script (NOS-018).
* The other copies of the shared files: phygas and planodia.

## Knowledge base

`docs/kb/` holds this repository's known bugs and platform traps, and the bug patterns shared with the other bitboinc
applications. Workflow: [docs/kb/shared/PROCESS.md](docs/kb/shared/PROCESS.md).

- Ready to fix, with next steps: `python3 docs/kb/shared/kb.py list --open --action fix` (also `--area`,
  `--platform`). Then `kb.py show NOS-NNN` for the issue and `kb.py show PAT-NNN` for the pattern it links, with the
  BOINC API facts the pattern cites. In the code, read the functions that `where:` and the Cause name.
- All repositories checked out next to this one: `kb.py todo`, then `kb.py todo --pattern PAT-NNN` for the next step
  and place of each issue; `kb.py check PAT-NNN` checks a pattern everywhere.
- Before debugging a symptom: `grep -ril "<words>" ../*/docs/kb/issues docs/kb/shared/patterns`, then `kb.py show`
  the ids it finds. Read other repositories through the shell (`kb.py show`, `git -C ../<repo> grep -n`, `sed -n`):
  opening one of their files in an editor tool can load their whole AGENTS.md as well.
- New report: `kb.py new issue "symptom in a few words" area=... next=...`.
- After a fix: `kb.py set NOS-NNN status=fix-proposed --log "..."` (`status=fixed fixed_in=<commit>` once merged),
  and update the passages of this file it prints. The commit names the issue (`KB: NOS-NNN`) and the pattern, if any
  (`Pattern: PAT-NNN`).
- Before a pull request: `kb.py index && kb.py lint --strict`.
