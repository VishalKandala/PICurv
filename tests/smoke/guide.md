# Smoke Fixtures Guide

This directory holds executable smoke assets for the canonical `make smoke` target.

The smoke runner verifies:

- `bin/simulator` launches and responds to `-help`
- `bin/postprocessor` launches and responds to `-help`
- `bin/picurv init` creates a self-contained case with copied binaries and origin metadata
- template matrix init/validate/dry-run coverage across `flat_channel`, `bent_channel`, and `brownian_motion`
- `picurv run --dry-run --format json` produces a valid solve/post plan
- restart planning resolves `run_control.restart_from_run_dir` into the expected restart source directory
- real end-to-end solve+post flow on a tiny flat-channel case
- real end-to-end solve+post flow on a tiny bent-channel case
- real end-to-end solve+post flow on a tiny particle-enabled flat-channel case
- restart flows for both particle restart branches (`load` and `init`) against a generated base run
- restart equivalence check: continuous tiny flat run vs split tiny flat restart (same final-step continuity metric tolerance)
- real end-to-end analytical Brownian case including particle VTP output and MSD statistics CSV
- multi-rank tiny solve+post runtime flows for flat and bent channels
- multi-rank flat-channel particle runtime flow with restart `load` and `init` branches

These checks are intentionally tiny, but they execute the real runtime path (solver, restart, and postprocessor),
not just dry-run wiring.

Canonical commands:

- `make smoke`
- `make smoke-mpi`
- `make smoke-mpi-matrix`
