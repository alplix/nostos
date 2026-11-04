---
id: NOS-015
title: Android, local Linux and Windows build scripts hard-code paths of one developer machine
type: bug
status: confirmed
severity: low
platforms: [all]
area: build
next: Build scripts: find the repo from the script's own path; take BOINC, NDK, CUDA and OpenCL-header paths from env vars
pattern: PAT-021
aliases: []
reported: 2026-10-27
source: audit 2026-10-27
where: scripts/build_boinc_android.sh:8-14, scripts/build_boinc_linux.sh:13, scripts/build_gpu_windows.bat:5, 9, scripts/build_gpu_opencl_wintest.bat:6, 11-12 @ eb7029e
fixed_in:
---

## Symptom

The Android build, the local Linux build and the two Windows batch files work only on one Windows PC with WSL and a
particular folder layout. Elsewhere they fail, or they quietly use headers or libraries from another checkout.

## Cause

Verified by reading the scripts at eb7029e (`git grep` as in
[PAT-021](../shared/patterns/PAT-021-build-depends-on-developer-paths.md)):

* `scripts/build_boinc_android.sh:8`: the repository is a fixed WSL path under `/mnt/c/Users/` (one user's Documents
  folder), not overridable; the NDK and the BOINC libraries are fixed paths under `/root` (`:10-14`).
* `scripts/build_boinc_linux.sh:13`: the same fixed repository path; `BOINC_SRC` and `BOINC_BUILD` can be
  overridden (`:14-15`).
* `scripts/build_gpu_windows.bat:9` and `scripts/build_gpu_opencl_wintest.bat:6`: the BOINC MSVC libraries are taken
  from a fixed folder under `C:\Users\` belonging to another project's checkout; `build_gpu_windows.bat:5` fixes the
  CUDA 13.3 install path; `build_gpu_opencl_wintest.bat:11-12` takes the OpenCL headers from `C:\temp`.

The build-server scripts (`scripts/build_boinc_all.sh`, `scripts/build_gpu_server.sh`) use `/opt/bitboinc-src` and
`/opt/llvm-mingw`; that is the build server's layout, stated in their headers, and is not counted here. The data
scripts' `D:/nostos_data` is NOS-009.

## Fix

Proposed: find the repository from the script's own location (`$(dirname "$0")/..`, `%~dp0..`), and take the BOINC,
NDK, CUDA and OpenCL-header locations from environment variables with today's values as defaults; list the variables
in AGENTS.md under Building.

## Test

Run each script from a fresh clone in another directory (or a container) with the variables set; it must build
without editing.

## Log

- 2026-10-27: found while checking PAT-021; verified by reading the scripts at eb7029e.
