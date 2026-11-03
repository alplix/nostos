---
id: NOS-007
title: Gaia download kept every finished future and its rows in memory until the end of the run
type: bug
status: fixed
severity: medium
platforms: [all]
area: data-pipeline
pattern:
aliases: [N2]
reported: 2026-10-21
source: audit 2026-10-21 (plan sections 1 and 6.5, N2)
where: scripts/download_filter_6d.py:26, 117-124 @ eb7029e
fixed_in: 0cfdedd
---

## Symptom

The memory use of `scripts/download_filter_6d.py` grew with every finished bulk file until the machine ran out
(the README's "Data" item: "real RAM-leak bug found and fixed along the way").

## Cause

The dict `futures` maps every submitted future to its key. A finished future holds its result, which is the header
and every kept row of that bulk file, and nothing removed the entries, so all results stayed in memory for the whole
run (over a thousand bulk files, as the code comment at `:121-124` says).

## Fix

Fixed upstream in 0cfdedd (2026-10-06), before the audit plan was written; the plan reviewed ec7aefe, the commit
before it. The loop deletes each future from the dict as soon as its rows are written (`:121`), and `N_WORKERS` went
from 16 to 8 (`:26`). Verified by reading at eb7029e: after the deletion only the loop variable refers to the
finished future, until the next iteration. All futures are still submitted at once (`:118`), but a future that has not
run yet holds only its key, so the bounded in-flight window the plan proposed is not needed for memory. Not measured
here. Found with `git log --oneline ec7aefe..HEAD -- scripts/download_filter_6d.py`.

## Test

None. A test would run `main()` with a mocked download that returns large row lists and check with `tracemalloc` that
the peak stays near `N_WORKERS` files' worth.

## Log

- 2026-10-21: reported in the audit plan as N2 (against ec7aefe).
- 2026-10-27: recorded from the 2026-10-21 audit; found fixed in 0cfdedd, verified by reading the code at eb7029e.
