---
id: NOS-006
title: Gaia download marks a bulk file done after five failed attempts, so its rows are lost for good
type: bug
status: confirmed
severity: high
platforms: [all]
area: data-pipeline
next: Find bulk files with no rows in the CSV (source_id HEALPix ranges) and fetch them; stop marking failed files done
pattern: PAT-023
aliases: [N1]
reported: 2026-10-21
source: audit 2026-10-21 (plan sections 1 and 6.5, N1)
where: scripts/download_filter_6d.py:85-97, 117-132 @ eb7029e
fixed_in:
---

## Symptom

The output CSV silently lacks every row of a Gaia bulk file whose download or filtering failed five times in a row.
A second run does not retry it, because `progress.txt` lists it as done. The only trace is a line
`GIVING UP on <key> after 5 attempts` in the console output of the run.

Lead, not verified: the README gives 33,530,050 rows for the finished download, while Gaia DR3 publishes a radial
velocity for 33,812,183 sources. If the CSV was meant to hold all of them, 282,133 rows (0.8%) are missing, about the
content of 28 bulk files if the 3386 bulk files hold about 10,000 radial-velocity rows each on average; and the CSV also holds some
rows twice (NOS-008), so the shortfall of distinct stars is larger. Lost files would explain it; another cause is
possible. Those stars are never searched by the BOINC work units cut from this CSV.

## Cause

Verified by reading the script at eb7029e (the lines are unchanged since ec7aefe; the only later commit, 0cfdedd,
touched other lines). `process_one()` retries five times and then returns `(key, None, [], 0)` (`:92-96`). `main()`
then treats it like a success: the empty row list is "written" and `mark_done(key)` appends the key to `progress.txt`
(`:120-132`).

## Fix

Proposed (plan N1): return a failure flag instead of an empty result, do not call `mark_done()` for a failed file,
list the failed keys at the end and exit with a non-zero status, so a rerun retries exactly those.

To repair the existing CSV: find the bulk files that contributed no rows. A Gaia DR3 `source_id` carries the HEALPix
index of the source (`source_id // 2**35` is the level-12 index), and the bulk file names carry a HEALPix range (check
the level in the Gaia DR3 bulk-download documentation), so one pass over the CSV's `source_id` column shows the empty
ranges without downloading anything. Download and filter those files and append them (after the deduplication of
NOS-008).

## Test

A pytest test with a mocked download that always fails for one key: that key must not appear in `progress.txt`, the
run must report it and exit non-zero, and a second run must try it again. No tests exist for the data scripts yet.

## Log

- 2026-10-21: reported in the audit plan as N1 (against ec7aefe), "to confirm".
- 2026-10-27: recorded from the 2026-10-21 audit, verified by reading the code at eb7029e.
- 2026-10-27: origin of pattern PAT-023 (case a); linked.
