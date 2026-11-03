---
id: NOS-011
title: Write errors on the result file are not checked before boinc_finish(0)
type: bug
status: confirmed
severity: low
platforms: [all]
area: io
next: Both apps: check fflush, ferror and fclose of out before boinc_finish(0); on failure print the error and boinc_finish(1)
pattern: PAT-006
aliases: []
reported: 2026-10-27
source: audit 2026-10-27
where: src/main_nostos_boinc.cpp:136-162, src/main_nostos_gpu_boinc.cpp:169-182 @ eb7029e
fixed_in:
---

## Symptom

A full disk, a quota or an I/O error while `out` is written goes unnoticed and the task reports success.

## Cause

Verified by reading the code at eb7029e. Neither app checks the result of `fprintf`, `fflush` (GPU `:173`) or
`fclose`, nor calls `ferror()`, before `boinc_finish(0)` (`src/main_nostos_boinc.cpp:136-140, 159-162`,
`src/main_nostos_gpu_boinc.cpp:169-173, 179-182`). See
[PAT-006](../shared/patterns/PAT-006-write-errors-unchecked-before-finish.md).

What limits the damage: the server's validator checks the `# done <n_stars>` trailer and the star count (commit
message of 1ecc658), so a file cut at the end is rejected. A row lost in the middle while the trailer is still
written (space freed in between) would pass.

## Fix

Proposed: before `boinc_finish(0)`, `if (fflush(outf) != 0 || ferror(outf) || fclose(outf) != 0)` print the error and
end with `boinc_finish(1)` (or a temporary exit once the resume can repair `out`, NOS-001). In the checkpoint path,
check `fflush(outf)` too, so a checkpoint never counts rows that failed to write.

## Test

On Linux, in a standalone run directory: `ln -s /dev/full out`, run the app: it must not end with status 0.

## Log

- 2026-10-27: found while checking PAT-006; verified by reading the code at eb7029e.
