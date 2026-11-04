---
id: NOS-008
title: Gaia download writes a bulk file's rows twice when it stops between writing them and marking the file done
type: bug
status: confirmed
severity: medium
platforms: [all]
area: data-pipeline
next: download_filter_6d.py: store the CSV size with each done key and cut back to it on start; dedupe source_id in make_star_catalog.py
pattern: PAT-023
aliases: [N3]
reported: 2026-10-21
source: audit 2026-10-21 (plan section 6.5, N3); NOTES_method_selection.md, 2026-10-07
where: scripts/download_filter_6d.py:107-132 @ eb7029e
fixed_in:
---

## Symptom

The downloaded Gaia CSV holds some rows twice. Observed by the maintainer (NOTES_method_selection.md, "First preview
scan with the fixed binary, and a small duplicate-row issue, 2026-10-07"): 237 of 56,767 candidate rows (0.4%) were
exact duplicates, and `grep -c` for one such `source_id` in the 43 GB CSV returned 2. The scan scripts deduplicate by
`source_id` (`scripts/summarize_candidates.py`, `scripts/make_shortlist.py`, `scripts/prepare_mc_shortlist.py`), but
`scripts/make_star_catalog.py` and `scripts/make_work_units.py` do not, so the star catalogue the BOINC work units are
cut from repeats those stars: they are searched twice and their encounter rows reach the server twice.

## Cause

Verified by reading the script at eb7029e. For each finished bulk file the rows are written and flushed
(`:125-131`), and only then is its key appended to `progress.txt` (`:132`). A stop between the two (Ctrl-C, a crash,
a reboot) leaves the rows in the CSV but not the key in `progress.txt`, so the next run downloads that file again and
appends all its rows again (the CSV is opened with `"a"` when it exists, `:107-110`). A stop inside the `writerow`
loop has the same effect for the part of the file that the write buffer had already flushed, and can leave a torn last
line that the next run's first row is glued to. The NOTES entry named this resume path as the most likely cause
without tracing it; the code reading confirms that the path exists. Not reproduced here.

## Fix

Proposed. The plan's order (write, flush, `os.fsync`, then mark done) narrows the window but cannot close it with
two files. Robust: record the CSV's byte length with each key in `progress.txt` and, at start, truncate the CSV to the
last recorded length before appending (this also drops torn lines); or write each bulk file's rows to a file of its
own and concatenate at the end. For the existing CSV: deduplicate by `source_id` in `scripts/make_star_catalog.py`
(keep the first occurrence, as the scan scripts do) and say so in its docstring.

## Test

A pytest test with a mocked download and a fault injected between the write and `mark_done()`: after a second run
every `source_id` appears once and every line parses. No tests exist for the data scripts yet.

## Log

- 2026-10-07: duplicates observed and worked around downstream (NOTES_method_selection.md).
- 2026-10-21: reported in the audit plan as N3 (against ec7aefe).
- 2026-10-27: recorded from the 2026-10-21 audit, verified by reading the code at eb7029e.
- 2026-10-27: instance of pattern PAT-023 (case b); linked.
