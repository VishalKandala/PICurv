# Smoke Fixtures Guide

This directory holds executable smoke assets for the canonical `make smoke` target.

The smoke runner verifies:

- `bin/simulator` launches and responds to `-help`
- `bin/postprocessor` launches and responds to `-help`
- `bin/picurv init` creates a self-contained case with copied binaries and origin metadata
- `picurv run --dry-run --format json` produces a valid solve/post plan
- restart planning resolves `run_control.restart_from_run_dir` into the expected restart source directory

These checks stay lightweight (no full CFD solve/post run), but now cover end-to-end workflow orchestration.
