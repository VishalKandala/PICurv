# Examples Guide

This directory contains runnable case templates used by `./bin/picurv init` and reference configs for sweeps/clusters.

## Subdirectories

- `flat_channel/`: programmatic-grid baseline case.
- `bent_channel/`: file-grid baseline case.
- `brownian_motion/`: analytical zero-flow particle-diffusion verification case.
- `master_template/`: fully annotated reference templates for all config roles.

## How To Use

- Initialize a starter case:
  - `./bin/picurv init flat_channel --dest my_case`
  - `./bin/picurv init bent_channel --dest my_case`
- Validate copied configs before running:
  - `./bin/picurv validate --case my_case/<case>.yml --solver my_case/<solver>.yml --monitor my_case/<monitor>.yml --post my_case/<post>.yml`
- If project binaries are already built, `init` copies the available executables from `bin/`
  into the new case directory, including `picurv`, so the case can run self-contained.
- `init` also writes `.picurv-origin.json`, so from inside that case you can later run
  `./picurv status-source`, `./picurv build`, `./picurv pull-source`, `./picurv sync-binaries`,
  and `./picurv sync-config`
  against the original source repository.
- Use master templates to discover full schema/options before trimming for production runs.
- The config roles are modular: you can often reuse a `solver.yml`, `monitor.yml`, or `post.yml` from one example with a different `case.yml` when the contract still makes sense.

## Related Docs

- https://vishalkandala.me/picurv-docs/02_Tutorial_Programmatic_Grid.html
- https://vishalkandala.me/picurv-docs/03_Tutorial_File-Based_Grid.html
- https://vishalkandala.me/picurv-docs/49_Workflow_Recipes_and_Config_Cookbook.html
- https://vishalkandala.me/picurv-docs/37_Sweep_Studies_Guide.html
