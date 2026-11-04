---
id: NOS-019
title: README named MWPotential2014 as the Galactic potential while the code uses Allen & Santillan (1991)
type: bug
status: fixed
severity: info
platforms: [all]
area: docs
pattern:
aliases: [plan 3.3 nostos gotcha]
reported: 2026-10-21
source: audit 2026-10-21 (plan section 3.3)
where: README.md:44, 46 @ eb7029e
fixed_in: cad9a30
---

## Symptom

The plan of 2026-10-21 noted for nostos: "The README mentions MWPotential2014, but phygas uses Allen & Santillan
1991." A reader would expect the wrong Milky Way model.

## Cause

At ec7aefe (the commit the plan reviewed) the README's roadmap said only that a Galactic potential still had to be
chosen, "e.g. `MWPotential2014`"; there was no code yet. `git log -S MWPotential2014 -- README.md` shows that cad9a30
(2026-10-06) removed that text, and 0c43404 then replaced the placeholder potential by Allen & Santillan (1991). At
eb7029e the README names Allen & Santillan (`README.md:44, 46`), `src/potential.c:28-45` implements it with the
published parameters, and `NOTES_method_selection.md` mentions `MWPotential2014` only in the history of that choice
(`:35-66`), as an option that was not taken. `src/potential.c` is identical to phygas's copy, so both applications use
the same potential.

## Fix

Nothing left to do; recorded so the plan's note is not followed up again.

## Test

`git grep -n MWPotential` finds only `NOTES_method_selection.md`, in the history of the choice.

## Log

- 2026-10-21: noted in the audit plan (section 3.3) against ec7aefe.
- 2026-10-27: recorded from the 2026-10-21 audit; already fixed in cad9a30, verified at eb7029e.
