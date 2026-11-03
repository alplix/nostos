---
id: NOS-004
title: GPU initialisation and kernel errors end the task with boinc_finish(1)
type: question
status: new
severity: medium
platforms: [cuda, opencl]
area: gpu-runtime
next: Decide: may GPU init and runtime errors end in boinc_temporary_exit() with a cap of about 5 retries without progress?
pattern: PAT-013
aliases: []
reported: 2026-10-27
source: audit 2026-10-27; question 4 of the 2026-10-21 plan (temporary exit on GPU errors)
where: src/main_nostos_gpu_boinc.cpp:73, 78, 159 @ eb7029e
fixed_in:
---

## Symptom

Any GPU error ends the task as a computation error: the GPU work of the task is lost although a checkpoint after the
last finished batch exists, and the host's error count goes up. Example: before 291f413 the Windows driver watchdog
killed long kernels (OpenCL error -36) and one host with a GTX 750 Ti failed 33 tasks (291f413 commit message, NOS-012).

## Cause

Verified by reading the code at eb7029e. Every error exit after start-up, classified as in
[PAT-013](../shared/patterns/PAT-013-transient-gpu-error-ends-task.md):

1. No usable device or initialisation failed: `ng_cuda_init()` (`src/nostos_cuda.cu:18-29`: no CUDA device,
   `cudaGetDeviceCount`/`cudaSetDevice` errors, which include a driver that is too old) and `ng_cl_init()`
   (`src/nostos_opencl.c:23-84`: no OpenCL library, no GPU of the build's vendor, no fp64, context, queue, program or
   kernel creation failed). Both end in `boinc_finish(1)` (`src/main_nostos_gpu_boinc.cpp:73`, `:78`).
2. Runtime error during the computation: `ng_cuda_run()` (`CK()` on allocation, copy, launch and synchronize) and
   `ng_cl_run()` (buffer allocation, kernel error, result read). They end in `boinc_finish(1)` (`:159`).
3. Bad input: unreadable or malformed `in` (`:97`), invalid orbital elements (`:102`). A computation error is right
   here.

Some causes in group 1 are permanent on that host (no fp64: the plan classes `opencl_*_fp64` should prevent it; a
kernel that does not build), others are transient (driver reset or update, device missing after a wake). Group 2 is
mostly transient (watchdog, out of memory while another program uses the GPU).

## Fix

A decision for the maintainers: should groups 1 (transient causes) and 2 use `boinc_temporary_exit()` with a retry
counter in the slot directory, reset when a batch finishes, and `boinc_finish(1)` after about five temporary exits
without progress, as the pattern proposes? The per-batch checkpoint already keeps the progress, so a retry would lose
at most one batch. Nothing to implement before the answer.

## Test

A test build where an environment variable makes the Nth launch fail: the run must end with a temporary exit and an
intact checkpoint, and repeated failures without progress must end in `boinc_finish(1)` after the limit.

## Log

- 2026-10-27: recorded as a question for the maintainers; exits verified by reading the code at eb7029e.
