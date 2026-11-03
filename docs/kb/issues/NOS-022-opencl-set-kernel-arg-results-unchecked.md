---
id: NOS-022
title: OpenCL launch loop ignores the results of clSetKernelArg
type: bug
status: confirmed
severity: low
platforms: [opencl]
area: gpu-runtime
next: ng_cl_run(): check the nine clSetKernelArg results and treat a failure like a failed launch (goto out, rc = 1)
pattern: PAT-027
aliases: []
reported: 2026-10-27
source: audit 2026-10-27
where: src/nostos_opencl.c:102-121 @ eb7029e
fixed_in:
---

## Symptom

None observed. If a `clSetKernelArg` call failed after the first slice, the next launch would run with that
argument's previous value and the batch would produce rows that look valid: a stale `step0` (argument 6) gives a
wrong time of closest approach (`t_min = (j + 1) * dt` in `src/nostos_kernel_body.h:84`), a stale step count
(argument 7) integrates the wrong number of steps and so gives wrong distances.

## Cause

Verified by reading the code at eb7029e. In every slice `ng_cl_run()` sets the nine kernel arguments and ignores the
nine results (`src/nostos_opencl.c:105-113`); the launch, `clFinish` and the result read are checked (`:114-116`,
`:124`). On the first slice an argument that was never set makes the launch fail with `CL_INVALID_KERNEL_ARGS`, which
is caught; on later slices the old value is still set. The host result buffer is freshly allocated per batch, so a
failed read cannot pass old data. The CUDA backend checks every call (`CK()`, `src/nostos_cuda.cu:36-67`). See
[PAT-027](../shared/patterns/PAT-027-unchecked-gpu-transfers-use-stale-data.md).

For a valid kernel, index and size, `clSetKernelArg` can practically fail only with `CL_OUT_OF_RESOURCES` or
`CL_OUT_OF_HOST_MEMORY`, hence low.

## Fix

Proposed: combine the nine results (`err |= ...` or check each) and treat any failure like a failed launch (`goto
out`, `rc = 1`). The six arguments that do not change within a batch can be set once before the loop.

## Test

A test build that passes a wrong size for argument 6 on the second slice must print "kernel failed" (or a new
"clSetKernelArg failed") and end without writing rows for that batch.

## Log

- 2026-10-27: found while checking PAT-027; verified by reading the code at eb7029e.
