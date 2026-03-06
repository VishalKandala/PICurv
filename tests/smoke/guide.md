# Smoke Fixtures Guide

This directory holds executable smoke assets for the canonical `make smoke` target.

The smoke runner verifies:

- `bin/simulator` launches and responds to `-help`
- `bin/postprocessor` launches and responds to `-help`
- `bin/picurv init` creates a self-contained case with copied binaries and origin metadata
- `picurv run --dry-run --format json` produces a valid solve/post plan
- restart planning resolves `run_control.restart_from_run_dir` into the expected restart source directory
- real end-to-end solve+post flow on a tiny flat-channel case
- real end-to-end solve+post flow on a tiny particle-enabled flat-channel case
- restart flows for both particle restart branches (`load` and `init`) against a generated base run
- real end-to-end analytical Brownian case including particle VTP output and MSD statistics CSV

These checks are intentionally tiny, but they execute the real runtime path (solver, restart, and postprocessor),
not just dry-run wiring.
