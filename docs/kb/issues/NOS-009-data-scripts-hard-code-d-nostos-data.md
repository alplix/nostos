---
id: NOS-009
title: Gaia data scripts hard-code D:/nostos_data and disagree on where the file list is
type: bug
status: confirmed
severity: low
platforms: [all]
area: data-pipeline
next: Read NOSTOS_DATA_DIR (default D:/nostos_data) in download_filter_6d.py and list_gaia_bulk.py; note it in AGENTS.md
pattern: PAT-021
aliases: [N4]
reported: 2026-10-21
source: audit 2026-10-21 (plan sections 3.3 and 6.5, N4)
where: scripts/download_filter_6d.py:21-24, scripts/list_gaia_bulk.py:34 @ eb7029e
fixed_in:
---

## Symptom

On a machine without a `D:` drive, `scripts/download_filter_6d.py` stops at once (it cannot open
`D:/nostos_data/gaia_file_list.txt`). `scripts/list_gaia_bulk.py` writes `gaia_file_list.txt` into the current
directory, so the two scripts only work together after the file is moved by hand.

## Cause

Verified by reading the scripts at eb7029e. `DATA_DIR = "D:/nostos_data"` is a constant with no override
(`download_filter_6d.py:21`); the file list, the output CSV and `progress.txt` are derived from it (`:22-24`).
`list_gaia_bulk.py:34` writes `gaia_file_list.txt` relative to the current directory. `.gitignore` lists those file
names as if they lived in the repository root. Other defaults are overridable and do not count:
`scripts/summarize_candidates.py:39` (command-line argument), `src/gaia_scan.c` (the CSV path is an argument;
`D:/nostos_data` appears only in a comment). Build-script paths are NOS-015.

## Fix

Proposed (plan N4): a `NOSTOS_DATA_DIR` environment variable (or a `--data-dir` argument) read by both scripts, with
`D:/nostos_data` as the unchanged default, and a line in AGENTS.md.

## Test

Run both scripts with `NOSTOS_DATA_DIR` pointing to a temporary directory and the HTTP calls mocked: the list, the
CSV and `progress.txt` must all land there.

## Log

- 2026-10-21: reported in the audit plan as N4 (against ec7aefe).
- 2026-10-27: recorded from the 2026-10-21 audit, verified by reading the code at eb7029e.
