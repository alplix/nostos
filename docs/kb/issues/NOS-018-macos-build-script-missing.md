---
id: NOS-018
title: scripts/build_boinc_mac.sh is referenced but not in the repository
type: bug
status: confirmed
severity: low
platforms: [macos]
area: build
action: decide
next: Maintainers: commit the macOS build script used for the Mac builds, or correct the reference in build_boinc_all.sh:4
pattern:
aliases: []
reported: 2026-10-27
source: audit 2026-10-27
where: scripts/build_boinc_all.sh:4 @ eb7029e
fixed_in:
---

## Symptom

The macOS CPU app (the 1ecc658 commit message reports Mac arm64 and x86_64 builds in its cross-platform check) cannot
be rebuilt from this repository, and its compiler flags cannot be reviewed (they matter for NOS-013).

## Cause

`scripts/build_boinc_all.sh:4` says "Mac (scripts/build_boinc_mac.sh) and Android (scripts/build_boinc_android.sh)
are built on their own hosts". There is no `scripts/build_boinc_mac.sh` at eb7029e, and
`git log --all -- scripts/build_boinc_mac.sh` finds no commit that ever had it.

## Fix

Proposed: commit the macOS build script (without developer paths, NOS-015), or correct the comment.

## Test

None.

## Log

- 2026-10-27: found while reading the build scripts at eb7029e.
- 2026-10-27: action set to decide: only the maintainers have the macOS build script, so no one else can fix this.
