---
id: NOS-017
title: README status and the make_work_units.py docstring predate the deployment and the GPU apps
type: bug
status: confirmed
severity: low
platforms: [all]
area: docs
next: Update README Status and Roadmap (deployed CPU and GPU apps), the make_work_units.py docstring and the cl_dynload.h comment
pattern:
aliases: []
reported: 2026-10-27
source: audit 2026-10-27
where: README.md:39-40, 49-50, scripts/make_work_units.py:1-12, src/cl_dynload.h:2, 18 @ eb7029e
fixed_in:
---

## Symptom

A reader of the README learns that the CPU app is "not yet deployed; no work generator, validator, or server-side
registration yet" (`README.md:39-40`, also the "BOINC app" roadmap item, `:49`) and that GPU acceleration is not
started (`:50`, unchecked box). `scripts/make_work_units.py:1-12` says that no BOINC project instance exists for
Nostos and calls itself the data half of a "future" work generator. `src/cl_dynload.h:2, 18` refer to
`opencl_search.c`, which does not exist in this repository (phygas's copy has the same text).

## Cause

The texts were not updated by 1ecc658 (2026-10-08, "Deploy Nostos on bitboinc: builds for every platform, server
side"), dea832d and 291f413 (the CUDA and OpenCL apps) or 17c83eb (3I/ATLAS added to the search). The commit messages
describe a server side (work generator, validator, assimilator, pages) that lives outside this repository. Verified by
reading the files at eb7029e.

## Fix

Proposed: update the README's Status and Roadmap (CPU and GPU apps, platforms, where the server side lives, which
objects are searched), the docstring of `scripts/make_work_units.py` (the server's generator cuts work units from
`scripts/make_star_catalog.py`'s catalogue), and the `cl_dynload.h` comment in every copy.

## Test

None (documentation).

## Log

- 2026-10-27: found while reading the repository; verified at eb7029e.
