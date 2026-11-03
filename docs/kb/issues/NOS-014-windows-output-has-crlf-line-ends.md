---
id: NOS-014
title: Windows builds write the result file with CRLF line ends
type: bug
status: confirmed
severity: low
platforms: [windows, mingw, msvc]
area: io
next: Open out with "wb" and "ab" in src/main_nostos_boinc.cpp:115 and src/main_nostos_gpu_boinc.cpp:129
pattern: PAT-019
aliases: []
reported: 2026-10-27
source: audit 2026-10-27; commit message of dea832d
where: src/main_nostos_boinc.cpp:115, src/main_nostos_gpu_boinc.cpp:129 @ eb7029e
fixed_in:
---

## Symptom

On Windows `out` ends every line with `\r\n`, so its bytes differ from the same result on every other platform. The
commit message of dea832d says the server's validator, compare and assimilator strip CR for this reason: otherwise
"every Windows result" would have been invalid.

## Cause

Verified by reading the code at eb7029e. Both apps open `out` with `"w"` or `"a"` and the checkpoint with `"w"` (text
mode), and the Windows C runtime turns `\n` into `\r\n`. Windows builds: the CPU app with MinGW-w64 and llvm-mingw
(`scripts/build_boinc_all.sh:41-43`), the OpenCL app with MinGW-w64 (`scripts/build_gpu_server.sh:64-76`), the CUDA
app with MSVC (`scripts/build_gpu_windows.bat`). The checkpoint is private to the app and holds one number read back
with `fscanf("%zu")`, so its CRLF is harmless (the corrected pattern counts such files as not affected); the input
parser accepts CRLF input (`fscanf` skips it as white space). Only `out` matters. See
[PAT-019](../shared/patterns/PAT-019-text-mode-output-crlf-on-windows.md).

## Fix

Proposed: open `out` with `"wb"` and `"ab"` (the checkpoint may follow for uniformity). It is needed anyway if the resume trims `out` by byte
length (NOS-001). The server's CR stripping can stay for older app versions.

## Test

Build with MinGW-w64, run a small work unit under Wine: `grep -c $'\r' out` must print 0 and `cmp` with the Linux
output must succeed (apart from NOS-013 differences, which do not occur on x86_64).

## Log

- 2026-10-27: found while checking PAT-019; verified by reading the code at eb7029e.
- 2026-10-27: re-checked after the pattern was corrected (private text files): the checkpoint is private and no longer
  part of the issue; title narrowed to the result file. `out` is still text mode, so still affected.
