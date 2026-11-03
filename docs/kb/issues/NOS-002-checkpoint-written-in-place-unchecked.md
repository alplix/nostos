---
id: NOS-002
title: Checkpoint rewritten in place, write results unchecked, completion reported even when the write failed
type: bug
status: confirmed
severity: medium
platforms: [all]
area: checkpoint
next: Write nostos_checkpoint via a .tmp file with checked fprintf/fflush/fclose and boinc_rename(); report completed() only on success
pattern: PAT-003
aliases: []
reported: 2026-10-27
source: audit 2026-10-27
where: src/main_nostos_boinc.cpp:143-150, src/main_nostos_gpu_boinc.cpp:174-176 @ eb7029e
fixed_in:
---

## Symptom

After a crash, a power loss or (GPU app) a quit at the wrong moment, `nostos_checkpoint` is empty. The next start
does not print "resuming from checkpoint", begins at star 0 and truncates `out`: all work done so far is computed
again. When the write fails (disk full, a file held by another program on Windows), the client is still told that a
checkpoint was made.

## Cause

Verified by reading the code at eb7029e. Of the four cases of
[PAT-003](../shared/patterns/PAT-003-checkpoint-not-written-atomically.md):

* (a) both apps open the only copy with `boinc_fopen("nostos_checkpoint", "w")`, which truncates it
  (`src/main_nostos_boinc.cpp:144`, `src/main_nostos_gpu_boinc.cpp:174`);
* (b) the results of `fprintf` and `fclose` are not checked (`:146-147`, `:175`);
* (c) does not apply: no rename is used;
* (d) `boinc_checkpoint_completed()` is called whether or not the file could be opened (`:149`, `:176`).

The CPU app writes inside the critical section that `boinc_time_to_checkpoint()` opened, so a client quit waits
until the write is done; only a crash, a kill or a power loss can cut it. The GPU app has no critical section: a quit
between the `fopen` and the `fprintf` leaves an empty file. (A quit after the `fprintf` is harmless: the library's
`fflush(NULL)` completes the file.)

An empty or unreadable checkpoint makes `fscanf` fail, the start becomes 0 and `out` is opened with `"w"`
(`src/main_nostos_boinc.cpp:110-115`, `src/main_nostos_gpu_boinc.cpp:126-129`). Output and checkpoint stay consistent,
so no row is wrong; the loss is work, up to a whole task (a CPU work unit of 1500 stars is about an hour; a GPU work
unit has 20,000 stars).

With (d), a failed write keeps the previous checkpoint while the client records a new one; the next resume then
starts at the older star and the rows written since are appended again (NOS-001).

## Fix

Proposed: one helper for both apps, as in the pattern: write `nostos_checkpoint.tmp` with `boinc_fopen(..., "wb")`,
check `fprintf`, `fflush` and `fclose`, then `boinc_rename()` it over `nostos_checkpoint`; call
`boinc_checkpoint_completed()` only when all of that succeeded, otherwise `boinc_end_critical_section()` in the CPU
path (PAT-001). Do it together with NOS-001, whose fix puts more into the checkpoint. No change for volunteers except
that an interrupted write no longer costs the task's progress.

## Test

With a stub BOINC library: make the write fail (a read-only slot directory, or `nostos_checkpoint.tmp` as a symbolic
link to `/dev/full` on Linux) and check that the old checkpoint is unchanged and no completion was reported. Kill the
process in a loop while it writes and check that the checkpoint always parses.

## Log

- 2026-10-27: found while checking PAT-003; verified by reading the code at eb7029e.
