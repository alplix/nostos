---
id: NOS-010
title: A checkpoint that cannot be opened restarts the task from star 0 and truncates the output
type: bug
status: confirmed
severity: low
platforms: [all]
area: checkpoint
next: Both apps: start from 0 only if !boinc_file_exists(nostos_checkpoint); any other open failure -> boinc_temporary_exit(300)
pattern: PAT-005
aliases: []
reported: 2026-10-27
source: audit 2026-10-27
where: src/main_nostos_boinc.cpp:106-115, src/main_nostos_gpu_boinc.cpp:124-129 @ eb7029e
fixed_in:
---

## Symptom

After a restart the task begins again at star 0 although it had a checkpoint (stderr has no "resuming from
checkpoint" line), and every finished star is computed again. Most likely on Windows hosts whose antivirus or backup
program holds the file at that moment.

## Cause

Verified by reading the code at eb7029e. Both apps treat any NULL from `boinc_fopen("nostos_checkpoint", "r")` as
"no checkpoint": the start stays 0 and `out` is opened with `"w"`, which truncates it
(`src/main_nostos_boinc.cpp:108-115`, `src/main_nostos_gpu_boinc.cpp:125-129`). NULL also comes back when the file
exists but cannot be opened (on Windows a sharing violation after 5 s of retries; EACCES; EMFILE), see
[PAT-005](../shared/patterns/PAT-005-open-failure-treated-as-fresh-start.md) and
[boinc-api.md](../shared/boinc-api.md#files). Output and start stay consistent, so the result is still correct; the
loss is the work done so far, up to a whole task.

A checkpoint that opens but does not parse (empty after NOS-002) also gives start 0 with a truncated `out`; the pattern
accepts that, because the two stay consistent.

## Fix

Proposed: a missing file (`!boinc_file_exists("nostos_checkpoint")`) means a fresh start; any other open failure ends
with `boinc_temporary_exit(300, "cannot open the checkpoint, will retry", false)` and leaves both files alone. Put it
in the checkpoint helper shared by both apps (NOS-001, NOS-002).

## Test

On Unix: run a standalone work unit until a checkpoint exists, `chmod 000 nostos_checkpoint`, start again: `out` must
not be truncated and the run must end with a temporary exit (with a stub BOINC library, or the real one in standalone
mode, which writes `temporary_exit`).

## Log

- 2026-10-27: found while checking PAT-005; verified by reading the code at eb7029e.
