---
id: NOS-003
title: GPU app writes the checkpoint after every batch without boinc_time_to_checkpoint()
type: bug
status: confirmed
severity: low
platforms: [cuda, opencl]
area: checkpoint
next: GPU app: write the per-batch checkpoint inside begin/end_critical_section(); call boinc_checkpoint_completed() only after a 1
pattern: PAT-004
aliases: []
reported: 2026-10-27
source: audit 2026-10-27
where: src/main_nostos_gpu_boinc.cpp:121, 132-178 @ eb7029e
fixed_in:
---

## Symptom

The volunteer's "request tasks to checkpoint at most every N seconds" setting has no effect on the GPU app, and its
checkpoint write is not protected against a quit.

## Cause

Verified by reading `src/main_nostos_gpu_boinc.cpp` at eb7029e. After every batch of `NOSTOS_BATCH_STARS` stars
(default 2000, `:121`) the app writes `nostos_checkpoint` and calls `boinc_checkpoint_completed()` (`:173-176`). It
never calls `boinc_time_to_checkpoint()`, so the write is outside any critical section and the completion is reported
without a preceding 1 (see [PAT-004](../shared/patterns/PAT-004-checkpoint-not-gated-by-time-to-checkpoint.md)).

Cadence and size: dea832d's commit message gives 500 stars in 27 s on an RTX 5070 Ti, so one batch is about 110 s
there and longer on slower GPUs. The file is one number. So the write is neither large nor frequent; what makes this
row `affected` is that the write is not atomic (NOS-002) and a quit can cut it. The CPU app asks the timer correctly
(`src/main_nostos_boinc.cpp:143`).

## Fix

Proposed. The batch is the natural checkpoint unit of the GPU app, and the client quits GPU tasks on every suspend,
so a checkpoint per batch is worth keeping. The pattern's `not-affected` form fits: keep the own cadence, write the
file atomically (NOS-002) between `boinc_begin_critical_section()` and `boinc_end_critical_section()`, and call
`boinc_checkpoint_completed()` only when `boinc_time_to_checkpoint()` (asked once per batch, PAT-001) returned 1.
The alternative, gating the write on the timer, loses up to a whole extra batch per quit for volunteers with a long
checkpoint interval.

## Test

A stub BOINC library that logs `boinc_time_to_checkpoint()` results, critical-section entries and checkpoint writes:
every write must be inside a begin/end pair, and `boinc_checkpoint_completed()` must follow only a 1.

## Log

- 2026-10-27: found while checking PAT-004; verified by reading the code at eb7029e.
