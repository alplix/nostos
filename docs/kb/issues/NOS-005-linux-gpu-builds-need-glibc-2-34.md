---
id: NOS-005
title: Linux GPU builds are linked against glibc 2.35 and need GLIBC_2.34, so they do not start on Debian 11 or Ubuntu 20.04
type: bug
status: confirmed
severity: high
platforms: [linux, x86_64, cuda, opencl]
area: build
next: build_gpu_server.sh: use the ubuntu20.04 images as phygas does and fail the build when the newest GLIBC_ is above 2.31
pattern: PAT-022
aliases: []
reported: 2026-10-27
source: audit 2026-10-27
where: scripts/build_gpu_server.sh:3-7, 36-38, 52-54, 59, 61, 78-80 @ eb7029e
fixed_in:
---

## Symptom

On Linux hosts whose glibc is older than 2.34 (Debian 11, Ubuntu 20.04, RHEL 8 and its rebuilds) the CUDA and OpenCL
apps do not start: the dynamic loader prints "version `GLIBC_2.34' not found" and every task errors out within a
second. No volunteer report recorded yet.

## Cause

Verified by reading `scripts/build_gpu_server.sh` at eb7029e. The Linux CUDA app is built in
`nvidia/cuda:12.8.1-devel-ubuntu22.04` (`:59`) and the Linux OpenCL apps in `ubuntu:22.04` (`:61`); Ubuntu 22.04 has
glibc 2.35. The link lines (`:36-38`, `:52-54`) link libstdc++, libgcc and cudart statically but glibc dynamically.
An executable linked against glibc 2.34 or newer references `__libc_start_main@GLIBC_2.34` (and the pthread
functions, merged into libc in 2.34, get `GLIBC_2.34` versions), so it needs glibc 2.34 at least. The commit that
added these builds says so itself: dea832d, "Linux CUDA + OpenCL (amd, nvidia) in Ubuntu 22.04 containers (GLIBC
2.34)". The script's comment "The Ubuntu 22.04 containers keep the GLIBC requirement low" (`:7`) is therefore wrong
for the bitboinc baseline: phygas's `scripts/build_gpu_server.sh`, which has the same structure, builds its GPU apps in
the Ubuntu 20.04 images for glibc 2.31 (Debian 11 / Ubuntu 20.04 and newer); apart from the app names and sources the
two scripts differ only in the image tags (diff of both at eb7029e and phygas 0a3eda2).

The CPU Linux builds are fully static (`scripts/build_boinc_all.sh:26`) and are not affected.

Not measured here (the Docker build was not run): the newest `GLIBC_` version in the shipped binaries. The script
prints it at the end (`:78-80`); to check a built binary:

```
objdump -T out/nostos_cuda_v1_linux_x86_64 | grep -oE 'GLIBC_[0-9.]+' | sort -V | tail -1     # expected: GLIBC_2.34
docker run --rm -v "$PWD/out:/o" debian:11 /o/nostos_opencl_v1_linux_x86_64_nvidia            # expected: GLIBC_2.34 not found
```

## Fix

Proposed: use `nvidia/cuda:12.8.1-devel-ubuntu20.04` and `ubuntu:20.04` as phygas does, correct the comment, and
make the script fail when the newest `GLIBC_` version of a Linux binary is above 2.31. It needs a new app version on
the server. Volunteers with glibc 2.31 to 2.33 can then run the GPU app.

## Test

After the build, the `objdump` line above must print `GLIBC_2.31` or lower, and the `debian:11` run must get past the
dynamic loader (it may then stop because no OpenCL or CUDA driver is present in the container, which is expected).

## Log

- 2026-10-27: verified by reading scripts/build_gpu_server.sh at eb7029e and the dea832d commit message; build not run.
