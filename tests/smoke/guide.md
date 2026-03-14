# Smoke Fixtures Guide

This directory holds executable smoke assets for the canonical `make smoke` family. Smoke tests are intentionally small but they run the real binaries and runtime pathways, so they provide stronger confidence than dry-run-only checks.

## What The Smoke Runner Verifies

- `bin/simulator` and `bin/postprocessor` launch and respond to `-help`.
- `bin/picurv init` produces a self-contained case with copied binaries and origin metadata.
- Template matrix init/validate/dry-run coverage across `flat_channel`, `bent_channel`, and `brownian_motion`.
- `picurv run --dry-run --format json` emits valid solve/post plans.
- Restart planning resolves `run_control.restart_from_run_dir` correctly.
- Tiny real solve+post paths for flat, bent, and particle-enabled flat cases.
- Restart branch coverage (`load` and `init`) including restart-equivalence checks.
- Analytical Brownian runtime path including particle output and MSD statistics.
- Multi-rank tiny runtime flows for flat/bent and particle restart branches.

## Why This Matters

These checks protect high-risk integration boundaries: case initialization, artifact generation, restart continuity, postprocessing compatibility, and MPI rank behavior.

When a smoke check fails, assume a user-visible workflow regression until proven otherwise.

## Canonical Commands

- `make smoke`
- `make smoke-mpi`
- `make smoke-mpi-matrix`

## Runtime Sequence Labels

The documentation pages refer to smoke sequences using these labels:

- `S0`: template matrix init/validate/dry-run checks
- `S1`: tiny flat solve+post
- `S1b`: tiny bent solve+post
- `S2`: tiny flat particle solve+post (default Trilinear interpolation)
- `S2b`: tiny flat particle solve+post with `CornerAveraged` interpolation (legacy path regression)
- `S3`: particle restart branches (`load`, `init`)
- `S4`: analytical Brownian path with particle VTP + MSD CSV
- `S5`: multi-rank tiny runtime paths (flat/bent + particle restart branches)
- `S6`: restart-equivalence continuity check

## Useful Environment Knobs

- `KEEP_SMOKE_TMP=1`: keep temp workspace for debugging.
- `SMOKE_MPI_NPROCS=<n>`: rank count for `make smoke-mpi`.
- `SMOKE_MPI_MATRIX_NPROCS="<n1> <n2> ..."`: rank list for `make smoke-mpi-matrix`.

## Maintenance Guidance

- Keep smoke runtime small enough for frequent local execution.
- Prioritize deterministic assertions over broad but noisy checks.
- Update fixture docs whenever smoke scope is extended.
