# Scripts Guide

This directory contains Python automation and utility scripts used for build/run orchestration and developer tooling.

## Primary Scripts

- `pic.flow`: main CLI conductor (`init`, `build`, `validate`, `run`, `sweep`).
- `grid.gen`: structured-grid generation utility.
- `audit_ingress.py`: checks YAML-to-C option ingress consistency.
- `check_markdown_links.py`: CI/local markdown link validator for docs.

## Typical Usage

- local run lifecycle:
  - `./scripts/pic.flow init ...`
  - `./scripts/pic.flow validate ...`
  - `./scripts/pic.flow run ...`
- study orchestration:
  - `./scripts/pic.flow sweep --study ... --cluster ...`

## Maintenance Rules

1. Keep CLI contract aligned with docs pages (`05`, `07`, `08`, `09`, `10`, `37`).
2. Prefer explicit validation errors with stable `ERROR <CODE>` formatting.
3. If new YAML keys map to C flags, update ingress docs and tests in same change.

## Related Docs

- `docs/pages/05_The_Conductor_Script.md`
- `docs/pages/15_Config_Ingestion_Map.md`
- `tests/guide.md`
