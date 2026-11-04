---
id: NOS-020
title: The orbit core and the OpenCL loader are copies of files in phygas and planodia
type: design-note
status: documented
severity: info
platforms: [all]
area: build
pattern:
aliases: [plan 3.3 nostos gotcha, plan question 7 (make_work_units.py)]
reported: 2026-10-21
source: audit 2026-10-21 (plan sections 3.3 and 8)
where: src/coords.c, src/potential.c, src/montecarlo.c, src/integrator.h, src/cl_dynload.c @ eb7029e
fixed_in:
---

## Symptom

A fix to the Milky Way potential, the coordinate transform, the random numbers, the integrator or the OpenCL loader
in one repository does not reach the other one by itself.

## Cause

Checked with `cmp` on 2026-10-27 against phygas 0a3eda2 and planodia 3f6c1b7: `src/coords.c`, `src/coords.h`,
`src/potential.c`, `src/potential.h`, `src/montecarlo.c`, `src/montecarlo.h`, `src/integrator.h`, `src/vec3.h`,
`src/units.h`, `src/cl_dynload.c` and `src/cl_dynload.h` are byte-identical to phygas's; `src/cl_dynload.h` is also
identical to planodia's, while `src/cl_dynload.c` differs from planodia's (NOS-016). They are registered as
`shared-code` in `docs/kb/kb.conf`. `scripts/build_gpu_server.sh` and `scripts/embed_kernel.py` have the same layout
in both repositories but are not copies.

The plan of 2026-10-21 said the "Nostos core" lived only in phygas: it reviewed nostos at ec7aefe, before the core was
committed here (cad9a30 onwards). Its question "where does make_work_units.py live?" is answered by the code:
`scripts/make_work_units.py` is in this repository since 2389f5b, and `scripts/make_star_catalog.py` (which imports
`zero_points` from it) since 1ecc658. phygas's own `scripts/make_star_catalog.py` imports the same function from a
nostos checkout next to it (default `../../nostos/scripts` relative to its `scripts/`, option `--nostos-scripts`).
Whether the core should live in one repository only is still a question for the maintainers.

## Fix

Nothing to change now. A change to one copy goes to the other copies in the same series of pull requests;
`python3 docs/kb/shared/kb.py status ../phygas ../planodia` reports `DRIFT` when they differ.

## Test

`kb.py status ../phygas ../planodia`: no drift except the one recorded in NOS-016.

## Log

- 2026-10-21: noted in the audit plan (sections 3.3 and 8) against ec7aefe.
- 2026-10-27: recorded from the 2026-10-21 audit; copies compared at eb7029e.
