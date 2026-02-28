# Examples Guide

This directory contains runnable case templates used by `./scripts/pic.flow init` and reference configs for sweeps/clusters.

## Subdirectories

- `flat_channel/`: programmatic-grid baseline case.
- `bent_channel/`: file-grid baseline case.
- `brownian_motion/`: analytical zero-flow particle-diffusion verification case.
- `master_template/`: fully annotated reference templates for all config roles.

## How To Use

- Initialize a starter case:
  - `./scripts/pic.flow init flat_channel --dest my_case`
  - `./scripts/pic.flow init bent_channel --dest my_case`
- Validate copied configs before running:
  - `./scripts/pic.flow validate --case my_case/<case>.yml --solver my_case/<solver>.yml --monitor my_case/<monitor>.yml --post my_case/<post>.yml`
- `init` links/copies compiled binaries (`picsolver`, `postprocessor`) when `bin/` exists.
- Use master templates to discover full schema/options before trimming for production runs.
- The config roles are modular: you can often reuse a `solver.yml`, `monitor.yml`, or `post.yml` from one example with a different `case.yml` when the contract still makes sense.

## Related Docs

- https://vishalkandala.me/picurv-docs/02_Tutorial_Programmatic_Grid.html
- https://vishalkandala.me/picurv-docs/03_Tutorial_File-Based_Grid.html
- https://vishalkandala.me/picurv-docs/49_Workflow_Recipes_and_Config_Cookbook.html
- https://vishalkandala.me/picurv-docs/37_Sweep_Studies_Guide.html
