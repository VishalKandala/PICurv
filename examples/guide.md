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
- Use master templates to discover full schema/options before trimming for production runs.

## Related Docs

- `docs/pages/02_Tutorial_Programmatic_Grid.md`
- `docs/pages/03_Tutorial_File-Based_Grid.md`
- `docs/pages/37_Sweep_Studies_Guide.md`
