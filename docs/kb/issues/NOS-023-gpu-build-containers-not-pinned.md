---
id: NOS-023
title: GPU build containers use moving image tags and unversioned packages, and no build records the BOINC version
type: bug
status: confirmed
severity: low
platforms: [linux, windows, cuda, opencl]
area: build
next: Pull the GPU build images by digest (with the NOS-005 switch) and print gcc, nvcc and BOINC versions into the build log
pattern: PAT-029
aliases: []
reported: 2026-10-27
source: audit 2026-10-27
where: scripts/build_gpu_server.sh:13, 18, 59, 61, scripts/build_boinc_all.sh:17-26 @ eb7029e
fixed_in:
---

## Symptom

Two builds of the same commit can differ (another compiler patch release, another statically linked libstdc++, other
BOINC libraries), and an older release cannot be rebuilt exactly to find out why its results differ.

## Cause

Verified by reading the scripts at eb7029e (`git grep` as in
[PAT-029](../shared/patterns/PAT-029-build-dependencies-not-pinned.md)):

* The GPU builds run in `nvidia/cuda:12.8.1-devel-ubuntu22.04` and `ubuntu:22.04` (`scripts/build_gpu_server.sh:59,
  61`), pulled by tag; both tags are re-published with distribution updates, and no digest is given.
* Inside them, `apt-get install -y build-essential m4 pkg-config file opencl-headers` installs unversioned packages
  (`:18`); the GCC, libstdc++ and libgcc that it brings are linked statically into the shipped binaries (`:37`, `:53`).
* The BOINC source is the build server's `/opt/bitboinc-src/boinc` (`:13`), and the CPU builds link the server's
  per-target BOINC libraries and cross compilers (`scripts/build_boinc_all.sh:17-26`); no script prints their
  versions. The local WSL build names BOINC 8.2.8 only in a comment (`scripts/build_boinc_linux.sh:3`); the Windows
  CUDA build fixes the toolkit by its install path (CUDA 13.3).

No script clones a repository or installs Python packages.

## Fix

Proposed: pull the images by digest (or a dated tag), together with the move to the Ubuntu 20.04 images of NOS-005;
print `gcc --version`, `nvcc --version` and the BOINC version (`git -C <boinc> describe`, or the version in its
`configure.ac`) into the build log next to the binaries in `out/`.

## Test

Build the same commit twice in clean containers some days apart and compare the printed versions and the binaries'
output on a small work unit.

## Log

- 2026-10-27: found while checking PAT-029; verified by reading the scripts at eb7029e.
