---
id: NOS-001
title: Output rows written after the last checkpoint are kept on resume and written a second time
type: bug
status: confirmed
severity: medium
platforms: [all]
area: checkpoint
next: Both apps: keep the found rows in nostos_checkpoint (next star + rows, written atomically) and write out once at the end
pattern: PAT-002
aliases: []
reported: 2026-10-27
source: audit 2026-10-27
where: src/main_nostos_boinc.cpp:106-150, src/main_nostos_gpu_boinc.cpp:124-177 @ eb7029e
fixed_in:
---

## Symptom

After a task was stopped and started again, its result file `out` holds the same `source_id` line twice, or a
torn line followed by a full copy of it, or two `# done` lines. The CPU app can also lose rows: the checkpoint counts a
star whose row never reached the disk. No volunteer report yet; found by reading the code.

## Cause

Verified by reading `src/main_nostos_boinc.cpp` and `src/main_nostos_gpu_boinc.cpp` at eb7029e, with the BOINC
facts in [boinc-api.md](../shared/boinc-api.md#exit-paths). The general mechanism is
[PAT-002](../shared/patterns/PAT-002-output-after-checkpoint-kept-on-resume.md).

CPU app. The resume position comes only from `nostos_checkpoint` (`:108-113`). When it is above 0, `out` is opened
with `"a"` (`:115`) and is not cut back. A row is written with `fprintf` when its star finishes (`:135-141`) and stays
in the stdio buffer. The checkpoint is written only when `boinc_time_to_checkpoint()` returns 1 (`:143-150`). When the
client quits the task, the BOINC library calls `fflush(NULL)` before it exits, so every row written since the last
checkpoint reaches `out`. The resume starts at the checkpointed star again and appends those rows a second time. The
seed of each star is `combine_seed(seed_base, source_id)`, so on the same host the second copy is byte-identical. A
quit in the middle of an `fprintf` (the flush runs from the SIGALRM handler) leaves a torn line before the copy.

The CPU app also writes the checkpoint without `fflush(outf)` first. If the process ends without the library's flush
(a crash, SIGKILL, the operating system ending processes at shutdown), the checkpoint counts stars whose rows were
still in the buffer, and those rows are missing from the final result. The GPU app flushes first (`:173`).

GPU app. The same resume logic (`:125-129`). After each batch it writes the batch's rows (`:161-172`), flushes `out`
(`:173`) and then rewrites the checkpoint (`:174-175`), outside any critical section. A quit between the first
`fprintf` of a batch and the `fclose` of the checkpoint keeps those rows and the old checkpoint, so the whole batch
(2000 stars by default) is computed and its rows appended again. The window is short (host-side work only), but the
client quits GPU tasks on every suspend, so the window is crossed often.

End marker. Both apps write `# done <n_stars>` after the loop without looking for one (`src/main_nostos_boinc.cpp:159`,
`src/main_nostos_gpu_boinc.cpp:179`). If the process ends after the trailer reached the file but before
`boinc_finish()` wrote `boinc_finish_called` (a quit in the short window before `boinc_finish()` sets `finishing`, or a
kill or power cut during its 2 s sleep), the client starts the task again. The GPU app's last checkpoint is `n_stars`
(written after the last batch), so the restart skips the loop and appends a second `# done` line; the CPU app resumes
at its last checkpoint and appends the rows after it and a second trailer. The window is small; noted because the
validator checks the trailer (commit message of 1ecc658).

Why medium and not the pattern's high: rows are sparse. In the 2026-10-07 catalogue scan about 0.3% of the screened
stars had a straight-line closest approach under 50 pc (the default threshold of `scripts/make_work_units.py`; the
server's threshold is not in this repository), so most quits find no row since the last checkpoint. With about 2.5 s
per star on the CPU and the default 60 s checkpoint interval, about 24 stars pass between two checkpoints. What the
server's quorum compare and assimilator do with a duplicated line is not visible here.

## Fix

Proposed, not implemented. Rows are few (a CPU work unit of 1500 stars gives a handful), so the simplest fix keeps
them in the checkpoint: the checkpoint file holds the next star index and every row so far, written atomically
(NOS-002), and `out` is written once, at the end, from memory. Then nothing is appended on resume and no trim is
needed.

The alternative is the pattern's option (a): at each checkpoint `fflush(outf)`, store the next star index and the
byte length of `out`; on resume copy the first that many bytes of `out` to a temporary file and `boinc_rename()` it
over `out` (a checkpoint that claims more bytes than `out` has means restart from 0), then append. This needs binary
mode for `out` (NOS-014), otherwise byte lengths differ on Windows. Option (b) of the pattern (the output is the
authority) does not fit: rows are written only for some stars, so the row count does not give the star index.

Either way, flush `out` before the checkpoint, and write `# done` only when it is not there yet (with the first fix,
`out` is written once at the end, so this comes for free). One helper used by both apps. Volunteers: no change to results; resumed tasks no longer carry
duplicate rows.

## Test

Needs a stub BOINC library (none in this repository yet) whose quit does what the real one does: `fflush(NULL)` and
`exit(0)` at a random star, not SIGKILL. Run a small work unit with a large threshold (so most stars give a row) to the
end, then again with repeated stub quits and restarts, and `cmp` the two outputs.

## Log

- 2026-10-27: found while checking PAT-002; verified by reading the code at eb7029e.
- 2026-10-27: re-checked after the pattern was corrected (flush before checkpoint, end markers): the missing flush in
  the CPU app was already described; added the second `# done` trailer case. Status and severity unchanged.
