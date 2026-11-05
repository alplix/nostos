# Benchmarking Nostos (`scripts/bench/`)

How to measure whether a change makes Nostos faster and whether it keeps the results, the scripts that do it, and a
place for the results. The method is the one the Kampylos app uses for its speed work (its document lives in that
repository, so the rules are repeated here, adapted to this app). The scripts build the real BOINC wrappers
(`src/main_nostos_boinc.cpp`, `src/main_nostos_gpu_boinc.cpp`) against a small fake BOINC runtime
(`scripts/bench/boinc_stub`) and run them on fixed work units. Written 2026-11-04.

**Status (2026-11-04):** the suite was written on an Apple M1 and first run on an Apple M4 Pro the same day, where it
built and ran without a fix. It now also runs Linux arm64 in a container (`linux_ab.sh`) and Windows and Linux x86_64
builds under emulation (`cross_check.sh`). Measured shares replace the estimates in "Where the time goes", and the
first same-results speed-up (candidates 1 and 2, `src/mc_search.c`) is measured in Results. The GPU apps are still
unmeasured.

## Method

### Two kinds of speed change

* **Same results**: build flags, link mode, allocator, GPU batch and launch sizes, and rewrites that do the same
  floating-point operations on the same values in the same order, only not twice (using a value that was already
  computed). The result file `out` must stay byte-identical to the old build's on the same toolchain. Such a change
  needs a measured gain, the output check and the correctness checks below, nothing more. `out` prints distances with
  six decimals, far fewer digits than a double has, so a change that touches arithmetic is checked with `EXACT=1`
  (17 digits) on the profile `quick`.
* **Result-changing**: anything that changes a number or a row: other or fewer `pow` calls or a rewritten formula,
  contraction into fused multiply-adds or fast-math, another step count, window or clone count, another RNG, seed or
  draw order, skipping stars or clones after a screen. Measure it the same way and treat it as a science change that the
  maintainers decide on: the server compares the two results of a quorum (`server/nostos_compare.sh`, not in this
  repository; whether it tolerates last-digit differences is the open question NOS-013), may mix app versions, and re-runs
  every found star with `nostos_detail`, whose distances must equal the volunteers' (728ea35). Until it is decided,
  keep such a change behind a compile-time switch that is off by default, and record the measured effect in the issue
  and in `NOTES_method_selection.md` (AGENTS.md, "Conventions for changes").

### Rules for every measurement

1. **Build both sides the same way, from clean copies.** The baseline is a git ref (`REF=upstream/main`) built with the
   documented build lines (`scripts/build_boinc_all.sh`: `-O2`, `-std=c99` for the C modules); the candidate is the
   working tree, or the baseline plus a patch. `build_variant.sh` copies the tracked files, so a stale object file or a
   scratch file in the working tree cannot leak in. (A new file of the candidate must be `git add -N`ed, or it is not
   copied.) Both sides link the same fake BOINC runtime.
2. **Find out where the time goes first.** Sample a run (`sample <pid> 10` on macOS, `perf record -g` on Linux) and, for
   a GPU app, split the kernel time from the host time. A change saves at most the share of time it touches: a change
   to 10 % of the time cannot gain more than 1.11x.
3. **Use a profile that exercises the change.** A change that skips work on far stars does nothing on `quick` (every
   star is written to `out` there, so none is skipped) and shows on `short` and `cpu-wu`; a change to arithmetic is
   checked on `quick` with `EXACT=1`; a change to GPU batching needs `gpu-batch` or more stars. Confirm the final
   number on a whole work unit (`cpu-wu`).
4. **Interleave and repeat.** `run_bench.sh` runs A B C A B C ... in fresh directories; use at least 3 rounds. The load
   from other programs changes over minutes, and interleaving makes it hit every build alike. Compare medians.
5. **Accept by the rule, not by eye.** A candidate is faster when its median is at least 5 % below the baseline's and
   its slowest run beats the baseline's fastest, or when its median is at least 10 % below. Anything else is "no clear
   gain" and is reported as that.
6. **Check the output of every run.** Each run's whole `out` is hashed: the rows `source_id min_pc mean_pc t_myr` and
   the trailer `# done <n_stars>`, nothing in it is a time, a host or a version. `same output` means every run of a
   label matches the baseline; a run whose last line is not `# done <n_stars>` fails. Compare like with like: other
   compilers, operating systems and libm (`pow`, `log`, `cos` of glibc, Apple, the CUDA and OpenCL built-ins) differ at
   the last digits (NOS-013), so a Linux GCC build is compared with a Linux GCC build and a CUDA build with a CUDA build.
7. **Measure on a quiet machine.** Pause the BOINC client (CPU and GPU tasks, Nostos tasks included) and other heavy
   programs, keep a laptop on power, and look at `uptime` (the scripts record the load). On a loaded machine an
   interleaved comparison is still fair, but its absolute times and spreads are not worth quoting.
8. **Write it down.** Date, machine, OS, compiler, load, profile, rounds, medians and spread, verdict and output go into
   the results below. A number in a commit message, a README or a code comment comes from there, and a platform that
   was not measured is called unmeasured.

### Steps for one change

1. Sample the current build on the profile closest to production (`short`, or `cpu-wu`) and note where the time goes.
2. Make the change in the working tree, or as a patch (`build_variant.sh <name> <patch>`; `git diff > x.patch` makes one);
   build the baseline with `REF=`.
3. Run one round on `quick` with `EXACT=1` (`scripts/bench/mac_ab.sh 1 quick`) to see that the output stays the same.
4. Run 3 or more rounds on a quiet machine on each profile the change touches (`mac_ab.sh 3 short`), and a whole work
   unit (`cpu-wu`) for the final number.
5. Run the correctness checks: `make test PYTHON=python3` and the tests that AGENTS.md routes the touched files to; for
   a change in the orbit code also the kernel mirror (`src/nostos_kernel_body.h`, then `python3 scripts/embed_kernel.py`)
   and a `cmp` of `out` between a CPU and a GPU build. There is no kill-and-resume test for the wrappers (see "What the
   suite does not cover").
6. Measure other platforms where possible: the same two scripts on Linux (GCC, static, as the release build), and the
   CUDA and OpenCL builds on a machine with the GPU. Name the platforms that were not measured.
7. Record the results here, and update the build lines in AGENTS.md and `scripts/build_*.sh` when the way to build changes.

### Traps carried over from the Kampylos measurements

These were seen with Kampylos on the same kind of Mac; Nostos has not been measured yet, but the causes are not specific
to either app.

* Other programs: on an M1 the same build varied by up to 40 % between rounds at load 5 to 30, and a BOINC task on the
  same GPU doubled a GPU run.
* Efficiency cores: `taskpolicy -b` (background QoS) puts a process on the efficiency cores of Apple Silicon, `nice` does
  not, and the BOINC client's way of lowering priority decides which a volunteer's Mac does. `RUN_PREFIX` compares both.
* Containers on a Mac run on the same cores and share their load: good for ratios between Linux builds, not for
  comparing Linux with macOS.

## One command on a Mac

```
scripts/bench/mac_ab.sh              # 3 rounds of the profiles quick and short (M4 Pro: about 6 minutes)
scripts/bench/mac_ab.sh 5 short      # 5 rounds of one profile
scripts/bench/mac_ab.sh 1 cpu-wu     # one production-sized work unit per build (M4 Pro: about 35 min base, 12 new)
scripts/bench/linux_ab.sh            # the same with GCC on Linux, in a Debian 12 container on a Mac (CLANG=1: also clang)
scripts/bench/cross_check.sh         # Windows (MinGW-w64, Wine) and Linux x86_64 builds, emulated: compile, run, which out
EXPERIMENTS=1 scripts/bench/mac_ab.sh 3 short   # also BASE_REF with -O3 and with -O2 -flto
NEW_OPT="-O2 -flto" scripts/bench/mac_ab.sh     # the working tree with other flags
EXACT=0 scripts/bench/mac_ab.sh      # production output format (6 decimals) instead of 17 digits
```

It needs the Xcode command line tools (or `gcc` and `g++`), `python3` and `git`; no BOINC checkout. It builds, in
`~/nostos-work`, the baseline `upstream/main` (or `origin/main`; `BASE_REF=` picks another ref) the documented way
(`-O2`, the build lines of AGENTS.md "Building" and `scripts/build_boinc_all.sh`) and this working tree with `NEW_OPT`
(default the same `-O2`), then runs each profile. On an unchanged tree the two builds are identical, so the table is an
A/A comparison: it shows the noise of this machine, and what a "no clear gain" looks like. Each profile ends with a table
like this one (`~/nostos-work/bench-runs/<date>-<profile>/summary.txt`); the numbers are placeholders for the format,
not a measurement:

```
label        runs  median s    min s    max s    cpu s  RSS MB  GPU%   vs 1st  verdict / output
cpu-base        3       M0       ..       ..       ..      ..     -   1.000x  baseline; same output
cpu-new         3       M1       ..       ..       ..      ..     -   M0/M1x  <ACCEPT | no clear gain | slower> (+x %); same output
```

`verdict` applies the rule of Rule 5. `same output` means every run of that label wrote the same `out` as the baseline,
`same output as X` the same as an earlier label (a CUDA build after CPU builds is compared with the CUDA build before
it, not with the CPU ones); otherwise the column shows how far the closest `min_pc` of the work unit moved.
`GPU%` is only filled with `GPU_SAMPLE=1`.

On Linux, and for the GPU apps, use the two scripts directly (`build_variant.sh` builds `cpu`, plus `cuda` when `nvcc` is
in `PATH` and `opencl` when `CL/cl.h` is installed):

```
REF=upstream/main scripts/bench/build_variant.sh base       # ~/nostos-work/bins/base/{cpu,cuda,opencl}
scripts/bench/build_variant.sh new                          # this working tree
scripts/bench/run_bench.sh short 3 base=$HOME/nostos-work/bins/base/cpu new=$HOME/nostos-work/bins/new/cpu
scripts/bench/run_bench.sh gpu-batch 3 base=$HOME/nostos-work/bins/base/cuda new=$HOME/nostos-work/bins/new/cuda
```

## The pieces

| file | what it does |
|---|---|
| `build_variant.sh <name> [patch ...]` | copies the tracked files (or `REF=<git ref>`), applies patches, compiles the C modules with `-std=c99` and the wrapper as `scripts/build_boinc_all.sh` does, links the fake BOINC runtime, and writes `bins/<name>/cpu` (and `cuda`, `opencl`: see below) with a `.info` file of how it was built. `OPT`, `CC`, `CXX`, `ALLOC`, `STATIC`, `EXACT` are the knobs (`build_variant.sh` without arguments lists them) |
| `run_bench.sh <profile> <rounds> label=binary ...` | makes the input, then interleaved runs in fresh directories: wall, user and system time, peak memory, a hash of `out`, and the closest `min_pc` of the run. `RUN_PREFIX` runs every binary through a command (`taskpolicy -b`, `nice -n 19`) |
| `summarize_bench.py <results.tsv>` | the table above |
| `make_bench_wu.py <profile.wu> <file>` | writes the input `in` of a profile; used by `run_bench.sh`, the md5 of its result goes into the results file |
| `profiles/*.wu` | number of stars, how many are real, seeds, window, steps, threshold of one work unit |
| `boinc_stub/` | the fake BOINC runtime: `boinc_api.h`, `filesys.h`, `boinc_stub.cpp` (what differs from the real library is in the comment at the top of the `.cpp`) |
| `mac_ab.sh` | everything above in one command |
| `linux_ab.sh` | `BASE_REF` against this tree with GCC (and clang with `CLANG=1`) on Linux, natively or in a Debian 12 container on a Mac (arm64 on Apple Silicon: native speed, shared cores) |
| `cross_check.sh` | the CPU app for Windows (MinGW-w64 GCC, static, run by Wine) and Linux x86_64 (GCC, static), `BASE_REF` and this tree, in an x86_64 Debian 12 container: compiles, runs, which `out`. Under Rosetta on Apple Silicon: no speed result |
| `patches/` | result-changing experiments, off unless their `-D` is given: `fast-pow.patch` (`-DNOSTOS_FAST_POW`), `skip-far.patch` (`-DNOSTOS_SKIP_MARGIN_PC=<pc>`); see Results |

The backends. **CPU**: built and run on macOS arm64 and on Linux. **CUDA**: built only where `nvcc` is found (Linux); the
line is `scripts/build_gpu_server.sh`'s, for the card's own architecture (`-arch=native`; `CUDA_GENCODE` overrides).
**OpenCL**: built only where `CL/cl.h` is found (Linux, or `BACKENDS=opencl` by hand); it runs only on a GPU with fp64.
There is **no Metal app**, and a Mac cannot run the GPU apps: Apple's GPUs have no fp64, and the apps refuse a device
without it (`src/nostos_opencl.c:52-55`). So on the M1 only the CPU app is measured, and the GPU backends are a gap:
`build_variant.sh` and the `gpu-batch` profile are written for them from the release build scripts and have not been
built or run anywhere.

The work units. The Gaia data of a real work unit is not in the repository, so `make_bench_wu.py` makes the stars, from a
fixed seed with integer arithmetic only (the same bytes on every machine and Python version), and puts two real stars in
front: HIP 3757 and HD 292249, whose Gaia astrometry is in `tests/test_pipeline_oumuamua_hip3757.c`. The header
(elements, covariance, window, steps) is written by `scripts/make_work_units.py`'s own writer, so the profiles follow
that file. Time per star does not depend on the star (every star runs 100 clones for 100,000 steps, without an early
exit), so made-up stars cost what real ones cost. What the data does decide is which stars come under the distance
threshold and are written to `out`: with the production threshold of 50 pc a random star never does, and the two real
ones do, as in a real work unit where most chunks return a trailer and a few rows. All profiles use 'Oumuamua's elements
and covariance (the only ISO `make_work_units.py` has) and `seed_base` 1; they share `star_seed`, so each profile
starts with the stars of the shorter ones. The times below are estimates from the 2.3-2.6 s per star of the 2026-10
reports (`NOTES_method_selection.md`: 23.05 s for 10 stars on the reference server, 15.4 s for 6 stars on a PC), not
measured on an M1.

| profile | stars | threshold | apps | estimated run | why |
|---|---|---|---|---|---|
| `quick` | 12 (2 real) | 100000 pc: every star is written | CPU | about 30 s on an M1 | smoke test, and the output check for arithmetic changes: the hash covers all numbers of all 12 stars (with `EXACT=1`, every digit) |
| `short` | 40 (2 real) | 50 pc | CPU | about 100 s on an M1 | the standard A/B: long enough for a median, the production output shape (2 rows and the trailer), shows a gain from skipping far stars |
| `cpu-wu` | 1500 (2 real) | 50 pc | CPU | about 1 hour on an M1 | one production CPU work unit (AGENTS.md: 1500 stars, about 2.5 s per star); the final number, once per candidate |
| `gpu-batch` | 2000 (2 real) | 50 pc | CUDA, OpenCL | a minute or two on a mid-range GPU (code comment) | one full batch of the GPU apps (`NOSTOS_BATCH_STARS` 2000, 200,000 star-clone pairs); the release GPU work unit is ten of these |

## Where the time goes

**Measured on 2026-11-04** (Apple M4 Pro, `sample` over 3 s of `quick`, EXACT=1 builds): in the baseline `pow` took
72 % of the samples (`pow` 63 % and its stub 9 %), `milkyway_accel()` without `pow` 22 %, `find_closest_approach()`
6 %. After candidates 1 and 2 (this tree): `pow` 55 %, `mc_search_stage1()` with the inlined clone loop 27 % (the
leapfrog arithmetic, the `sqrt` of each step's distance and the stored star positions), `milkyway_accel()` 18 %. One
star took 1.26 s in the baseline and 0.46 s now (`short`, 40 stars). The estimate below got the order right (the three
`pow` the largest share, the rest invisible) and the per-star time wrong by 2x (it was the M1's 2.3-2.6 s).

The estimate of 2026-11-04, read from the code before anything was run:

* **CPU app: `milkyway_accel()`, almost all of it.** One star is 100 clones x 100,000 steps; one step is two
  `leapfrog_step()` calls (star, clone) and each calls `potential_acceleration()` twice (`src/integrator.h:21`,
  `src/encounter.c`). That is 4 x 100,000 x 100 = 4.0e7 acceleration calls per star, plus one `sqrt` for the distance per
  step. At the reported 2.3-2.6 s per star that is 58-65 ns per call, which fits what one call does (`src/potential.c`):
  3 `pow` (the bulge's `pow(.., 1.5)` at line 71, two in the halo's `halo_enclosed_GM` at line 57), 4 `sqrt` (`R`, `r`,
  the disk's `zeta` and `D`) and 9 divisions, called through a function pointer from another translation unit. Expected
  split of a call, to be checked: the three `pow` the largest share, then the divisions and square roots, the indirect
  call and `Vec3` copies a small share.
* **Everything else is invisible on the CPU app.** Per star: 300 Gaussian draws (Box-Muller: `log`, `sqrt`, `cos`), one
  3x3 Cholesky and one coordinate conversion, tens of microseconds against seconds. Per task: parsing `in`
  (`fscanf` per star), the ISO's start state, a checkpoint file write once per client period, a few output rows. The
  allocator is not involved (one `realloc`-grown star array). Peak memory is a few hundred KB.
* **GPU apps: the kernel, `ng_steps()` in `src/nostos_kernel_body.h`**, the same arithmetic in fp64, one thread per
  (star, clone) pair, so the same 4 calls per step. Consumer GPUs run fp64 at a small fraction of their fp32 rate, so
  the batch is bound by fp64 `pow`, `sqrt` and division. The host side per batch of 2000 stars is the draw of 200,000
  clones (600,000 Gaussians, tens of milliseconds, single-threaded), filling and uploading about 20 MB, and the reduction
  per star; the kernel takes a minute or two. Slices of about 0.25 s each end with a synchronisation
  (`ng_next_slice()`), at most 5000 steps per launch, so at least 20 launches per batch and an idle gap of tens of
  microseconds between them: well under 1 %. The host thread probably spins while it waits (CUDA's default
  `cudaDeviceSynchronize`), so a GPU run's CPU time is not the host's work.
* **The star's own trajectory is computed 100 times per star**, once per clone, though it is the same every time (the
  bodies do not attract each other). That is half of all acceleration calls. And inside a step the acceleration at the
  end of one step is computed again at the start of the next. Together this is the source of the first two candidates
  below.

## Candidate speed-ups, in the order to measure them

Ranked by expected gain and by how little a decision they need. Each says whether it keeps the results.

1. **Integrate the star once per star, not once per clone. Same results.** `find_closest_approach()`
   (`src/encounter.c`) advances the star and the clone together; the star's state is identical for all 100 clones.
   Store the star's positions (100,001 x 3 doubles, 2.4 MB) when the first clone runs, and have the other 99 compare
   against them. The arithmetic and its order are unchanged, so `out` is byte-identical. Expected: about 2x on the CPU
   app (2 of the 4 calls per step per clone go away). `mc_search_stage1()` is the only caller on the app's path; the
   other callers (`scan_logic.c`, `nostos_detail.c`, tools) keep the old function. Not directly available to the GPU
   apps, where each thread owns one pair: sharing a trajectory across the 100 threads of a star needs another
   decomposition of the kernel, a design question, not a quick patch.
2. **Carry the end-of-step acceleration into the next step. Same results.** `leapfrog_step()` computes `a1` at the
   new position and the next step computes `a0` at the same position: the same function of the same value, so reuse
   cannot change a bit (on builds that do not fuse operations differently; the hash decides). Expected: about 2x on the
   CPU app, and with item 1 about 4x (1 call per clone-step instead of 4). The same change applies to the kernel body
   (`ng_leapfrog`), which halves the GPU's work: keep the acceleration in registers within a slice and recompute it
   once at the start of each slice, so the checkpointed state does not grow. `src/integrator.h` is a copy shared with
   other apps (NOS-020): do it in `find_closest_approach()` or change every copy together. The kernel must mirror the CPU
   arithmetic operation by operation (AGENTS.md "Invariants"), then `embed_kernel.py`, then a CPU against GPU `cmp`.
3. **Let the compiler see the potential: `-flto`, or one translation unit for the C modules. Same results if nothing
   is fused differently.** `milkyway_accel()` is called through a function pointer from another file, so without
   link-time optimisation it cannot be inlined into the step loop. `EXPERIMENTS=1 mac_ab.sh` measures `-flto`.
   Expected: 5-15 % on the CPU app, and relatively the same after items 1 and 2.
4. **Other flags: `-O3`, a newer instruction set. Same results only if the compiler does not contract.** GCC with
   `-std=c99` and Apple clang 14 emit no fused multiply-adds for this code, upstream clang does by default
   (`-ffp-contract=on`) and then gave other output (NOS-013); `-O2 -ffp-contract=off` keeps it safe. A build for the
   host's own CPU (`-march=native`) cannot be shipped, since volunteers' CPUs differ. Expected: under 5 % (scalar code
   around libm calls); listed because it is one line to try. `EXPERIMENTS=1` measures `-O3`.
5. **Fewer `pow` calls. Result-changing at the last bits.** `pow(t, 1.5)` equals `t * sqrt(t)` mathematically and the
   halo's `x^gamma` equals `x * x^(gamma-1)`, so the three `pow` of a call could become one. Expected: 1.3-1.7x per call
   if the `pow` are the largest share of a call, to be checked by sampling. Every acceleration then differs in its last bits, so
   printed digits of some rows move and rows at the 50 pc threshold can appear or vanish (NOS-013), `nostos_detail`'s
   re-run and quorums with older app versions are affected, and the GPU kernel has to follow. Needs the maintainers'
   decision; measure how many rows of `cpu-wu` change and by how much, behind a compile switch.
6. **Skip the clones of stars that cannot be reported. Result-changing in principle.** The app runs all 100 clones for
   every star, but writes only stars whose closest clone is under the threshold, which almost no star is. Running one
   clone first and the rest only when it comes within a margin of the threshold (or the straight-line screen of
   `src/gaia_scan.c`) would skip about 99 % of the work for a typical star. The rows written stay the same as long as the
   margin is safe, which has to be shown on real chunks (the spread of the 100 clones is set by the covariance, about
   0.02 km/s per axis, i.e. about 2 pc over 100 Myr, small against 50 pc, but not proven for every star), and a wrongly
   skipped star is a missing row. The largest gain of this list and a method decision: measure on `short` and `cpu-wu`
   only, since `quick` reports every star and would show no gain.
7. **Vector instructions across clones. Result-changing unless `pow` stays scalar.** Add, multiply, divide and `sqrt`
   give identical results in vector form; libm has no vector `pow` that matches the scalar one bit for bit. Little to gain
   while `pow` is a large share; revisit after item 5.
8. **GPU scheduling: `NOSTOS_BATCH_STARS`, `NOSTOS_STEPS_PER_LAUNCH`, the block size (128 CUDA, 64 OpenCL).
   Same results** (the code says that splitting steps does not change a result). Expected to matter little because the
   kernel dominates, except where registers spill (a `double` array per thread) or a GPU is shared. GPU machine only.
9. **Allocator, static against dynamic linking, `CXXSTD`. Same results.** Expected: nothing. `ALLOC=jemalloc` and `STATIC=0`
   exist to show that.

Not a candidate for volunteers' throughput: threads inside one task. BOINC runs one task per core, so a second
thread per task shortens one task and not the work done per hour.

## Traps specific to Nostos

* **Determinism.** A run of one binary on one work unit writes the same `out` every time, so any difference between two
  runs of a label is a real finding (the summary says "results differ between its own runs"). Between builds it holds
  only on the same compiler, flags, libm and OS (NOS-013); a Linux hash is never compared with a macOS hash. `out` is
  six decimals of pc: a change of 1e-12 relative in a distance almost never changes it, so a "same output" at the
  production format does not prove "same arithmetic". Use `EXACT=1` for that, on `quick`, where every star's three
  numbers are in `out`. On `short` and `cpu-wu` only the two real stars are, and a change that alters other stars'
  numbers is invisible there.
* **The threshold hides or shows work.** The profile threshold decides which rows exist. `quick` raises it to 100,000 pc so
  that every star is written (a star's closest approach is never larger than its separation from the ISO at the start,
  and the made-up parallaxes of 1-20 mas put every star within a kiloparsec); the other profiles use the production 50 pc. Do not compare a `quick` time with a `short` time
  as if the work differed: it does not.
* **RNG seeds.** Every clone is drawn from `combine_seed(seed_base, source_id)` in clone order (`src/nostos_io.c`,
  `src/montecarlo.c`), so a change of `seed_base`, of the star list or of the order of draws changes every result. The
  profiles fix both; the md5 of the generated `in` is in the header of every `results.tsv`. If two sessions show
  different md5s for one profile, `scripts/make_work_units.py` (elements, covariance, format) or the generator changed,
  and their results are not comparable. A change of the clone count or of the draw order is a result-changing change.
* **Checkpoints.** The CPU app writes `nostos_checkpoint` after a star when the stub says it is time (every 60 s by default,
  `NOSTOS_STUB_CHECKPOINT_SECONDS`), the GPU app after every batch whatever the clock says (NOS-003). Every run starts in
  a fresh directory without a checkpoint, so it never resumes; a leftover `nostos_checkpoint` would make the app skip
  the first stars and append to an old `out`. A change to checkpointing (NOS-001 to NOS-003, NOS-010) cannot be judged
  here: the stub does not kill, and a kill would not flush stdio the way the real client's quit does
  (`docs/kb/shared/boinc-api.md`, PAT-002).
* **GPU sharing.** The kernel slices (`ng_next_slice`) follow the time one launch takes, so a GPU that another task or the
  desktop also uses runs with shorter slices and the whole batch takes longer, without changing a result. Measure GPU apps
  on a GPU that nothing else uses (pause the BOINC client's GPU tasks, no display on it where the driver has a watchdog,
  NOS-012). `GPU_SAMPLE=1` records the utilisation as a check of that. A host thread that spins in `cudaDeviceSynchronize`
  makes the CPU column of a GPU run look like work.
* **Thread counts.** Both apps are single-threaded and the stub reports no CPU budget; nothing here depends on core
  count. The client runs one task per core, so on a machine with simultaneous multithreading two tasks may share one
  core's floating-point units and each run slower than alone; `run_bench.sh` runs one at a time. A gain in the number of
  `pow` calls may weigh differently there.
* **The task's fixed costs.** The real `boinc_finish()` waits 2 s (`docs/kb/shared/boinc-api.md`); the stub does not, so
  a run here is 2 s shorter than a task, about 7 % of `quick` and 0.06 % of `cpu-wu`.
* **Long runs.** A fanless laptop may throttle during `cpu-wu`; interleaving spreads that over both sides but does not
  remove it. Look at the spread, not only the median.

## What the suite does not cover

* GPU apps on a Mac (not possible: no fp64), and any GPU measurement at all until someone runs `build_variant.sh` and
  `gpu-batch` on a machine with a card. The CUDA and OpenCL build lines in `build_variant.sh` are untested, and the CUDA
  one builds for one architecture, not the release list of sm_50 to sm_120 (a binary for the card in use runs the same
  code).
* Windows on real hardware, Android, 32-bit x86, Windows ARM64 (llvm-mingw, a clang) and the other release targets;
  the release builds themselves (`scripts/build_boinc_all.sh` needs the build server, NOS-015). `cross_check.sh` builds
  the Windows x86_64 CPU app with MinGW-w64 GCC (the release uses llvm-mingw for some targets) and runs it under Wine.
* The real BOINC library: the fake runtime has no timer thread, no suspend, quit or abort, no low-priority setting, and
  `boinc_finish()` exits at once (comment in `boinc_stub/boinc_stub.cpp`).
* Checkpoint and resume (no kill harness exists for the wrappers), the quorum compare, the validator, `nostos_detail`
  and the other server-side programs (`src/nostos_detail.c` has no build target), the data scripts, `tools/`, and
  `src/gaia_scan.c`, which has its own, different hot loop.
* Many tasks at once, as on a volunteer's machine (shared caches, simultaneous multithreading, thermal limits).
* The science: that a faster build finds the same encounters is shown only by `out` on these work units.
* Real Gaia chunks: the made-up stars cost what real ones cost but are not the distribution of a real chunk, so a
  change that skips work by screening (candidate 6) needs a run on real chunks before anyone relies on its gain.

## Results

### 2026-11-04: the star integrated once per star, the acceleration carried over (Apple M4 Pro)

Machine: Apple M4 Pro (8 performance + 4 efficiency cores, 24 GB), macOS 27.0, Apple clang 21.0.0 (Command Line Tools),
on power, BOINC client's tasks suspended, load 1.6-2.9 (the benchmark itself 1, the window server and BOINC Manager about
0.7). Builds: base = `upstream/main` (2fd4153), `-O2`, the documented lines; new = this working tree: `src/mc_search.c`
integrates the star once per star (candidate 1) and carries the end-of-step acceleration into the next step
(candidate 2), `-O2`. EXACT=1 (17 digits) unless said otherwise. 3 interleaved rounds, medians.

| profile | build | median s | min-max s | vs base | output |
|---|---|---|---|---|---|
| `quick` (12 stars, all written) | base | 15.47 | 15.21-15.68 | 1.00x | |
| | new | 5.50 | 5.46-5.51 | **2.81x**, ACCEPT | same, 17 digits |
| | base `-O3` | 15.57 | 15.36-15.78 | 0.99x, no clear gain | same |
| | base `-O2 -flto` | 15.63 | 15.45-15.72 | 0.99x, no clear gain | same |
| `short` (40 stars, 50 pc, 8 rows) | base | 52.62 | 52.00-52.82 | 1.00x | |
| | new | 18.45 | 18.37-18.48 | **2.85x**, ACCEPT | same |

The rows are byte-identical because the arithmetic is: the star's state after each step is the same for all 100 clones
(the bodies do not attract each other), and the acceleration a step computes at its new position is the one the next
step would compute there. `make test` (all 17 C tests and `test_make_work_units.py`) passes with the change. Expected
was about 4x; the remaining cost is the `pow` calls (one acceleration per clone-step is left) and the loop itself.

Other platforms, the same two builds (`linux_ab.sh`, `cross_check.sh`):

| platform | compiler | base | new | output |
|---|---|---|---|---|
| Linux arm64 container, 3 rounds, `quick` / `short` | GCC 12.2 | 14.16 s / 47.33 s | 5.32 s / 17.73 s (**2.66x / 2.67x**) | same |
| Linux arm64 container, 3 rounds | Debian clang 14.0.6 | 13.63 s / 45.47 s | 5.25 s / 17.54 s (2.70x) | same as clang base, other digits than GCC |
| Linux x86_64 container under Rosetta, 1 round, `quick` | GCC 12.2, static | 33.5 s | 8.4 s | same |
| Windows x86_64, MinGW-w64 GCC 16.2, static, Wine 8.0 under Rosetta, 1 round, `quick` | | 28.5 s | 8.1 s | same |

The emulated times only show that the builds run; they are not speed results.

**Results across platforms (NOS-013), measured.** On `quick` at 17 digits, every GCC build writes the same `out`:
Linux arm64 GCC 12, Linux x86_64 GCC 12 and Windows MinGW-w64 GCC 16 (hash f24ba714e9c3, base and new alike). The
clang builds differ from them: Apple clang 21 (6d11025beb46) and Debian clang 14 (fe54d8d5ac39). Apple clang 21 with
`-ffp-contract=off` writes exactly the GCC `out` (f24ba714e9c3; 74 of 200 rows of the 200-star check below differ from
clang's default build), at no speed cost (93.6 s against 93.0 s, one run). So the difference of the clang builds is
fused multiply-adds, which clang performs by default even with `-std=c99`; `-ffp-contract=off` on the clang targets
(macOS, Windows ARM64 and the other llvm-mingw builds, Android) would make every CPU build give the same rows. Not
applied: the release build scripts need the build server, and the change moves those builds' rows (NOS-013).

### 2026-11-04: result-changing experiments (Apple M4 Pro)

`scripts/bench/patches/` holds two experiments, each off unless its `-D` is given, built on top of this tree:

* `fast-pow.patch`, `-DNOSTOS_FAST_POW` (candidate 5): `(r^2+b^2)^1.5` as `d * sqrt(d)` and the halo's `x^gamma` as
  `x * x^(gamma-1)`, one `pow` instead of three per acceleration. On this Mac's libm each form differs from the
  original in the last bit for about 27 % of arguments.
* `skip-far.patch`, `-DNOSTOS_SKIP_MARGIN_PC=<pc>` (candidate 6): the star's closest approach to the ISO's mean orbit is
  computed first (one extra integration), and a star farther than threshold + margin skips its 100 clones and is not
  written.

| profile | new | fast-pow | skip-far, 100 pc margin | output |
|---|---|---|---|---|
| `short`, 3 rounds | 18.45 s | 15.79 s (1.17x) | 13.66 s (1.35x) | all three the same 8 rows, 17 digits |
| `quick`, 3 rounds (nothing is skipped: threshold 100,000 pc) | 5.58 s | 4.78 s (1.17x) | | same 12 rows, 17 digits |
| 200 stars, all written (`quick` with `stars 200`), 1 run | 93.0 s | 79.1 s (1.18x) | | 1 of 200 rows differs |

Fast `pow` changes the arithmetic, yet the rows rarely move: the half-step kick `a * dt/2` is about a millionth of the
velocity, so a last-bit change of the acceleration is usually lost when the velocity is rounded. The one differing row
of 200 differs in the mean distance's 16th digit (73.008710699751688 against 73.00871069975166), which the production
format (`%.6f`) does not print. It is still a science change in principle (and the GPU kernel would have to follow),
so it stays a switch: the maintainers decide, with this number.

Skipping far stars kept the rows of `short`, but the made-up stars of the profiles all pass within a kiloparsec of the
ISO, most within the margin, so the 1.35x says little. A skipped star costs 2 integrations instead of 101 (about 50x
less); the gain on a real chunk is the share of its stars that stay farther than threshold + margin, which needs a real
chunk to measure, and the margin needs the spread of the clones (about 2 pc per 100 Myr from the covariance, nonlinear
effects not bounded). A wrongly skipped star is a missing row.

### Next steps by platform (plan of 2026-11-04)

| platform | app | next step | expected | evidence |
|---|---|---|---|---|
| every CPU build | CPU | release this tree (candidates 1 and 2) | about 2.7-2.9x, same rows | above |
| CUDA, OpenCL | GPU | carry the acceleration in `ng_leapfrog` (`src/nostos_kernel_body.h`), then `embed_kernel.py` | up to 2x on the kernel (half the acceleration calls), same rows | not measurable here (needs an fp64 GPU) |
| CUDA, OpenCL | GPU | integrate the star once per star: another kernel decomposition (a star's 100 clones share its track) | up to 2x more | design question |
| clang targets (macOS, Windows ARM64, Android) | CPU | `-ffp-contract=off`, decide (NOS-013) | no speed change; the rows of every platform agree | above |
| all | CPU, GPU | fast `pow`, decide | 1.17x on the CPU after candidates 1 and 2; 1 row in 200 moves in the 16th digit | above |
| all | CPU, GPU | skip far stars, decide, after a run on a real chunk | up to ~50x per far star | above |
| Linux | CPU | nothing else: `-O3` and LTO gained nothing; the allocator is not used in the loop | | above |

One section per measurement session, newest first, in this form.

```
### <date>: <change or question> (<machine>)

Machine: <chip, cores, RAM>, <OS version>, <compiler --version>, power, BOINC client paused (yes/no), load before the
run (the `# load:` line of results.tsv).
Builds: base = <ref and commit, OPT>; new = <working tree at commit, patch, OPT>; EXACT=<0|1>; the `.info` files of
bins/<name>/.
Profile: <name>, rounds <n>, RUN_PREFIX <none|taskpolicy -b|nice -n 19>, md5 of `in` (from results.tsv).

<the table of summary.txt>

Verdict: <by the rule>. Output: <same | other, how far>. Where the time went: <sample result>. Not measured: <platforms,
backends, profiles>.
```
