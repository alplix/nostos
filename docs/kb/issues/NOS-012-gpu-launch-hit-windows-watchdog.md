---
id: NOS-012
title: GPU kernel launches of a fixed 500 steps ran past the Windows driver watchdog on GPUs with weak fp64
type: bug
status: fixed
severity: medium
platforms: [windows, opencl, cuda]
area: gpu-runtime
pattern: PAT-016
aliases: []
reported: 2026-10-27
source: commit message of 291f413 (2026-10-09)
where: src/nostos_gpu.h:33-67, src/nostos_opencl.c:100-121, src/nostos_cuda.cu:46-58 @ eb7029e
fixed_in: 291f413
---

## Symptom

OpenCL error -36 (`CL_INVALID_COMMAND_QUEUE`, after the driver reset) and a computation error; one host with a
GTX 750 Ti failed 33 tasks (291f413 commit message).

## Cause

The first GPU version (dea832d) ran the 100,000 integration steps in launches of a fixed 500 steps
(`NOSTOS_STEPS_PER_LAUNCH`, default 500) over all 200,000 (star, clone) pairs of a batch. On a GPU with weak double
precision one such launch took over 2 s, past the Windows TDR limit. See
[PAT-016](../shared/patterns/PAT-016-long-gpu-launch-hits-watchdog.md).

## Fix

291f413 (2026-10-09): every launch is timed (`ng_now()`), the first slice is 20 steps (`NG_FIRST_SLICE`), and each next
slice is sized to take about 0.25 s (`NG_TARGET_LAUNCH_S`), growing at most 4x per launch and capped by
`NOSTOS_STEPS_PER_LAUNCH` (now default 5000) (`src/nostos_gpu.h:33-67`). Both backends use it
(`src/nostos_opencl.c:100-121`, `src/nostos_cuda.cu:46-58`). The commit message says the output is unchanged (checked
against v1 and the CPU app): splitting the steps does not change the arithmetic. Verified by reading the code at
eb7029e. Remaining worst case is the first slice: from the figure above (500 steps over 2 s), 20 steps take about
0.1 s on that GPU. Not measured here.

## Test

Each batch prints `nostos_opencl: N pairs, M launches, last slice K steps` (or `nostos_cuda: ...`) to stderr. The
duration of each launch is not logged; logging it once per batch on the slowest supported GPU with a full batch would
show the margin.

## Log

- 2026-10-09: fixed upstream in 291f413 after the failures on one host.
- 2026-10-27: recorded while checking PAT-016; fix verified by reading the code at eb7029e.
