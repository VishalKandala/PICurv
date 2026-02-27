# Examples Guide

This directory contains runnable case templates used by `./scripts/pic.flow init` and reference configs for sweeps/clusters.

## Subdirectories

- `flat_channel/`: programmatic-grid baseline case.
- `bent_channel/`: file-grid baseline case.
- `master_template/`: fully annotated reference templates for all config roles.

## How To Use

- Initialize a starter case:
  - `./scripts/pic.flow init flat_channel --dest my_case`
  - `./scripts/pic.flow init bent_channel --dest my_case`
- Validate copied configs before running:
  - `./scripts/pic.flow validate --case my_case/<case>.yml --solver my_case/<solver>.yml --monitor my_case/<monitor>.yml --post my_case/<post>.yml`
- `init` links/copies compiled binaries (`picsolver`, `postprocessor`) when `bin/` exists.
- Use master templates to discover full schema/options before trimming for production runs.

## Related Docs

- https://vishalkandala.me/picurv-docs/02_Tutorial_Programmatic_Grid.html
- https://vishalkandala.me/picurv-docs/03_Tutorial_File-Based_Grid.html
- https://vishalkandala.me/picurv-docs/37_Sweep_Studies_Guide.html
