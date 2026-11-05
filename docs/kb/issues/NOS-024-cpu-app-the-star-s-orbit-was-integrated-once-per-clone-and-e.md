---
id: NOS-024
title: CPU app: the star's orbit was integrated once per clone and each acceleration computed twice
type: bug
status: fix-proposed
severity: low
platforms: [all]
area: performance
action: review
next: Port the acceleration carry to ng_leapfrog in src/nostos_kernel_body.h (GPU, needs an fp64 GPU to measure)
blocked_by: []
pattern:
aliases: []
reported: 2026-11-04
source: M4 benchmark 2026-11-04 (docs/BENCHMARK.md, Results)
where: src/mc_search.c:28-50, src/encounter.c:17-31 @ 4bce441
fixed_in:
---

## Symptom

Wasted time, no wrong result: a CPU task spent about 2.8 times the CPU time it needs (1.26 s per star on an Apple M4
Pro, 0.46 s with the fix, same output to 17 digits).

## Cause

`mc_search_stage1()` called `find_closest_approach()` once per clone (`src/mc_search.c` at 4bce441), which advanced the
star and the clone together (`src/encounter.c:17-31`). The star's state after each step is the same for all 100 clones,
because the bodies do not attract each other, so 99 of its 100 integrations repeated the first. And `leapfrog_step()`
(`src/integrator.h:21-29`) computed the acceleration at the start of a step that the previous step had just computed
at the same position. Together: 4 acceleration calls per clone-step where about 1 is needed. Sampled: `pow` (three per
acceleration) took 72 % of the time.

## Fix

`src/mc_search.c` (this branch): the star is integrated once per star and its positions are stored (2.4 MB at 100,000
steps; without the memory the old path runs), and each clone is compared with them; both integrations carry the
end-of-step acceleration into the next step (`leapfrog_step_carry()`). Same arithmetic in the same order, so the rows
are bit-identical. Kept out of `integrator.h` and `encounter.c`, which are copies shared with phygas (NOS-020) and used
by other callers. The GPU kernel is unchanged and still gives the same numbers as the CPU app; the carry would halve its
acceleration calls too (`next:`), and sharing the star's track among a star's clones needs another kernel decomposition.

## Test

`make test` passes. `EXACT=1 scripts/bench/mac_ab.sh 3 "quick short"`: same `out` at 17 digits, `quick` 15.47 s ->
5.50 s, `short` 52.62 s -> 18.45 s on an Apple M4 Pro; `linux_ab.sh`: 2.66-2.70x with GCC 12 and clang 14;
`cross_check.sh`: Windows (MinGW-w64) and Linux x86_64 builds write the same `out` as before.

## Log

- 2026-11-04: created.
- 2026-11-04: fix written and measured (docs/BENCHMARK.md, Results).
