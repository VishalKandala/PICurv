# Scripts Guide

This directory is the Python/shell control plane for PICurv. It covers user workflow orchestration (`picurv`), mesh generation, docs hygiene, ingress consistency checks, and coverage gates used by `make coverage-*`.

## Script Inventory and Purpose

- `picurv`
  - main user CLI: `init`, `build`, `sync-binaries`, `sync-config`, `pull-source`, `status-source`, `validate`, `run`, `sweep`
  - validates YAML role contracts, emits structured errors, generates solver/post runtime artifacts, and launches local or Slurm workflows
- `grid.gen`
  - standalone structured-grid utility for `grid.mode: grid_gen` and legacy 1D-axis conversion
  - writes PICGRID plus optional `.info`/`.vts` outputs
  - includes `legacy1d` subcommand for headerless legacy payload -> canonical PICGRID conversion
- `audit_ingress.py`
  - scans `src/setup.c` and `src/io.c` for `PetscOptionsGet*`/`PetscOptionsHasName` usage
  - compares discovered flags with `scripts/audit_ingress_manifest.json`
- `check_markdown_links.py`
  - local markdown link checker for `README.md`, `docs/**/*.md`, and `examples/**/*.md`
- `python_coverage_gate.py`
  - stdlib-trace-based Python line coverage gate used by `make coverage-python`
  - default target file is `scripts/picurv`
- `c_coverage_gate.py`
  - gcov summarizer and threshold gate used by `make coverage-c`
  - computes weighted line coverage over `src/*.c`
- `bootstrap_install.sh`
  - Debian/Ubuntu bootstrap helper for base packages, Python deps, optional PETSc build, and PICurv binary build
- `generate_doxygen_fallback_indexes.py`
  - helper for generating structured/fallback Doxygen index pages (`files*.html`, `annotated*.html`) when default indexes are incomplete
- `convert_grid_from_legacy_to_picgrid.py`
  - backward-compatible wrapper that delegates to `grid.gen legacy1d`
  - retained for transition workflows; canonical implementation now lives in `grid.gen`

## High-Value CLI Contracts (`scripts/picurv`)

- structured errors always use:
  - `ERROR <CODE> | key=... | file=... | message=... | hint=...`
- `run --dry-run --format json` is the safest way to inspect planned commands/artifacts
- in cluster mode:
  - solver rank count is derived from `cluster.yml` resources
  - post stage remains single-task by default
- `validate --strict` adds additional filesystem and study-base-config checks

## Typical Usage Flows

- build binaries:
  - `./scripts/picurv build`
- initialize a case and keep it synchronized:
  - `./bin/picurv init flat_channel --dest my_case`
  - `cd my_case && ./picurv status-source --format json`
  - `./picurv pull-source && ./picurv build && ./picurv sync-binaries`
  - `./picurv sync-config --prune`
- preflight before execution:
  - `./bin/picurv validate --case ... --solver ... --monitor ... --post ... --strict`
  - `./bin/picurv run --solve --post-process ... --dry-run --format json`
- cluster workflow:
  - `./bin/picurv run --solve --post-process ... --cluster cluster.yml --no-submit`
- sweep workflow:
  - `./bin/picurv sweep --study study.yml --cluster cluster.yml --no-submit`

## Maintenance Rules

1. Keep `picurv` CLI behavior and docs synchronized, especially option preconditions and defaulting rules.
2. Any new YAML-to-C option ingress must update all of:
   - `scripts/audit_ingress_manifest.json`
   - docs ingestion map (`docs/pages/15_Config_Ingestion_Map.md`)
   - regression tests (`tests/test_config_regressions.py` and/or CLI smoke tests)
3. Keep error-code semantics stable for automation.
4. Prefer deterministic dry-run and no-submit behaviors for new orchestration features.

## Debugging Guidance

- for run-plan issues:
  - use `run --dry-run --format json` and inspect `stages.*.launch_command_string`
- for cluster issues:
  - compare generated `scheduler/*.sbatch` with `cluster.yml`
- for docs drift:
  - run `python3 scripts/check_markdown_links.py`
- for ingress drift:
  - run `python3 scripts/audit_ingress.py --show-scanned`
- for coverage gate failures:
  - inspect `coverage/python/summary.txt` and `coverage/c/summary.txt`

## Related Docs

- https://vishalkandala.me/picurv-docs/05_The_Conductor_Script.html
- https://vishalkandala.me/picurv-docs/15_Config_Ingestion_Map.html
- https://vishalkandala.me/picurv-docs/40_Testing_and_Quality_Guide.html
- https://vishalkandala.me/picurv-docs/48_Grid_Generator_Guide.html
- `tests/guide.md`
