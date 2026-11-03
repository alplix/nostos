---
id: NOS-021
title: Windows CUDA build has no code for GPUs below compute capability 7.5; does the plan class exclude them?
type: question
status: new
severity: medium
platforms: [windows, cuda]
area: build
next: Decide: require compute capability 7.5 and a CUDA 13 driver in the Windows CUDA plan class, or build it with CUDA 12.8?
pattern: PAT-024
aliases: []
reported: 2026-10-27
source: audit 2026-10-27
where: scripts/build_gpu_windows.bat:5, 10, 15-17, scripts/build_gpu_server.sh:31, 59 @ eb7029e
fixed_in:
---

## Symptom

Expected, not observed: on a Windows host with a Maxwell, Pascal or Volta NVIDIA GPU (for example GTX 750 Ti, GTX 9xx,
GTX 10xx, Titan V), or with a driver too old for CUDA 13, the CUDA app fails at once ("no kernel image is available
for execution on the device" at the first launch, or an insufficient-driver error in `ng_cuda_init()`) and ends with
`boinc_finish(1)` (NOS-004). The same GPU works with the Linux CUDA app.

## Cause

Verified by reading the scripts at eb7029e. The Windows CUDA app is built with CUDA 13.3 (`build_gpu_windows.bat:5`)
and code for sm_75 to sm_120 plus compute_120 PTX (`:10`); the PTX cannot be compiled for older GPUs, and CUDA 13 no
longer supports building for sm_50 to sm_72. The Linux CUDA app is built with CUDA 12.8.1 for sm_50 and newer
(`build_gpu_server.sh:31, 59`). The runtime is linked statically (nvcc default), so the Windows app also needs a driver
that supports CUDA 13 (the r580 series or newer).

The plan classes are on the server and not in this repository, so whether such hosts get Windows CUDA tasks is not
known here.

## Fix

A decision for the maintainers: require compute capability 7.5 and a CUDA 13 driver in the Windows CUDA plan class, or
build the Windows app with CUDA 12.8 and the same gencode list as Linux.

## Test

Check the plan class on the server; run the Windows binary standalone on a Pascal GPU and read stderr.

## Log

- 2026-10-27: found while reading the build scripts at eb7029e.
- 2026-10-27: origin of pattern PAT-024; linked.
