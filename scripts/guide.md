# Scripts Guide

This directory contains Python automation and utility scripts used for build/run orchestration and developer tooling.

## Primary Scripts

- `picurv`: main CLI conductor (`init`, `build`, `sync-binaries`, `sync-config`, `pull-source`, `status-source`, `validate`, `run`, `sweep`).
- `grid.gen`: structured-grid generation utility.
- `audit_ingress.py`: checks YAML-to-C option ingress consistency.
- `check_markdown_links.py`: CI/local markdown link validator for docs.

## Typical Usage

- bootstrap/install:
  - `./scripts/picurv build`
- local run lifecycle after build:
  - `./bin/picurv init ...`
  - `cd my_case && ./picurv status-source && ./picurv build && ./picurv sync-binaries`
  - `./bin/picurv validate ...`
  - `./bin/picurv run ...`
- study orchestration:
  - `./bin/picurv sweep --study ... --cluster ...`
- direct grid generation:
  - `python3 scripts/grid.gen --config config/grids/coarse_square_tube_curved.cfg pipe`

## Maintenance Rules

1. Keep CLI contract aligned with docs pages (`05`, `07`, `08`, `09`, `10`, `37`).
2. Prefer explicit validation errors with stable `ERROR <CODE>` formatting.
3. If new YAML keys map to C flags, update ingress docs and tests in same change.

## Related Docs

- https://vishalkandala.me/picurv-docs/05_The_Conductor_Script.html
- https://vishalkandala.me/picurv-docs/48_Grid_Generator_Guide.html
- https://vishalkandala.me/picurv-docs/15_Config_Ingestion_Map.html
- https://vishalkandala.me/picurv-docs/49_Workflow_Recipes_and_Config_Cookbook.html
- `tests/guide.md`
