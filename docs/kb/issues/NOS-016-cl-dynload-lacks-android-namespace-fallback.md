---
id: NOS-016
title: cl_dynload.c lacks planodia's Android linker-namespace fallback (3f6c1b7)
type: bug
status: confirmed
severity: info
platforms: [android, opencl]
area: device-selection
next: Copy planodia's src/cl_dynload.c (3f6c1b7) here and to phygas so the registered copies are identical again
pattern: PAT-017
aliases: []
reported: 2026-10-27
source: audit 2026-10-27
where: src/cl_dynload.c:26-60 @ eb7029e
fixed_in:
---

## Symptom

None today. nostos has no Android OpenCL build: `scripts/build_boinc_android.sh` builds only the CPU app, and dea832d
states that no Adreno build is made because those GPUs lack fp64. If an Android OpenCL build were added with this
loader, many phones would report "no OpenCL ICD found" inside the BOINC client (see
[PAT-017](../shared/patterns/PAT-017-android-opencl-library-not-found.md)).

## Cause

Checked with `cmp` on 2026-10-27: `src/cl_dynload.c` is byte-identical to phygas's copy (phygas 0a3eda2) and differs
from planodia's, where 3f6c1b7 added `android_dlopen_ext()` with the vendor linker namespaces, vendor library paths and
the `dlerror()` text in the error message. `src/cl_dynload.h` is identical in all three. The POSIX branch here only
calls plain `dlopen()` on `libOpenCL.so`, `libOpenCL.so.1` and `libmali.so` (`:45-56`). `kb.py status ../planodia`
reports this as DRIFT.

## Fix

Proposed: copy planodia's `src/cl_dynload.c` (3f6c1b7) here and to phygas, so the registered copies are identical
again. Outside Android it changes only the error message (it adds `dlerror()`), so Linux and Windows behaviour stays
the same.

## Test

`python3 docs/kb/shared/kb.py status ../phygas ../planodia` shows no drift for `cl_dynload.c`; on Android only if an
Android OpenCL build is ever added (the test in PAT-017).

## Log

- 2026-10-27: found while checking PAT-017; copies compared at eb7029e (planodia 3f6c1b7, phygas 0a3eda2).
