# Scripts Guide

This directory is the Python/shell control plane for PICurv. It covers user workflow orchestration (`picurv`), mesh generation, docs hygiene, ingress consistency checks, and coverage gates used by `make coverage-*`.

## Script Inventory and Purpose

- `picurv`
  - main user CLI: `init`, `build`, `sync-binaries`, `sync-config`, `pull-source`, `status-source`, `validate`, `precompute`, `run`, `sweep`
  - validates YAML role contracts, emits structured errors, generates solver/post runtime artifacts, and launches local or Slurm workflows
- `grid.gen`
  - standalone structured-grid utility for `grid.mode: grid_gen` and legacy 1D-axis conversion
  - writes PICGRID plus optional `.info`/`.vts` outputs
  - includes `legacy1d` subcommand for headerless legacy payload -> canonical PICGRID conversion
- `profile.gen`
  - standalone dimensional `PICSLICE` profile generator used by `picurv` generated and field-sliced profile orchestration
  - currently supports `square_duct_poiseuille`
- `ic.gen`
  - standalone expression-based `Ucat`/`Ucont` initial-condition generator
  - evaluates fields on staged nondimensional PICGRID cell or face geometry and writes PETSc binary vectors
- `plot.gen`
  - standalone renderer for versioned normalized time-history JSON requests
  - used by `picurv summarize --plot`; owns matplotlib styling, interactive display, saving, and headless fallback
  - requires `matplotlib`, which `bootstrap_install.sh` installs by default
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
- `audit_function_docs.py`
  - repository-wide function-doc audit for C headers/sources/tests and Python product/test code
  - used by repo-consistency tests to prevent doc-coverage regressions

## CLI Help and Option References

Use script-local `--help` as the first source of truth:

- `./scripts/picurv --help`
- `./scripts/picurv run --help`
- `./scripts/picurv precompute --help`
- `./scripts/picurv validate --help`
- `./scripts/picurv sweep --help`
- `python3 scripts/grid.gen --help`
- `python3 scripts/grid.gen legacy1d --help`
- `python3 scripts/profile.gen --help`
- `python3 scripts/ic.gen --help`
- `python3 scripts/plot.gen --help`
- `python3 scripts/audit_ingress.py --help`
- `python3 scripts/check_markdown_links.py --help`
- `python3 scripts/python_coverage_gate.py --help`
- `python3 scripts/c_coverage_gate.py --help`
- `python3 scripts/generate_doxygen_fallback_indexes.py --help`
- `python3 scripts/convert_grid_from_legacy_to_picgrid.py --help`
- `python3 scripts/audit_function_docs.py`
- `bash scripts/bootstrap_install.sh --help`

Detailed option coverage lives in:

- `picurv`: `docs/pages/05_The_Conductor_Script.md` (`Full Command and Option Matrix`)
- `grid.gen`: `docs/pages/48_Grid_Generator_Guide.md` (`Option Families`)
- testing/coverage flow context: `docs/pages/40_Testing_and_Quality_Guide.md`

Helper script option summary:

- `audit_ingress.py`
  - `--manifest` to select manifest JSON path.
  - `--show-scanned` to print discovered PETSc flags before drift comparison.
- `check_markdown_links.py`
  - `--repo-root`, `--docs-dir`, `--examples-dir` to scope scanning.
  - `--no-readme` to exclude root `README.md`.
- `python_coverage_gate.py`
  - `--target` (repeatable) to select measured files.
  - `--pytest-args` to forward pytest arguments.
  - `--min-line`, `--output-dir` to control gate threshold and artifacts.
- `c_coverage_gate.py`
  - `--src-dir`, `--obj-dir`, `--output-dir` for coverage input/output paths.
  - `--min-line` for weighted line threshold.
- `generate_doxygen_fallback_indexes.py`
  - `--repo-root`, `--html-dir` for fallback index generation scope.
- `convert_grid_from_legacy_to_picgrid.py`
  - `--input`, `--output`, `--axis-columns`, `--allow-trailing` for legacy payload conversion.
- `bootstrap_install.sh`
  - `--install-petsc` to build PETSc in the bootstrap path.
  - `--petsc-version`, `--petsc-prefix`, `--petsc-arch` for PETSc installation details.
  - `--python-bin` to pick interpreter for dependency install/build checks.
  - `--skip-system-deps` to skip apt package installation.

## High-Value CLI Contracts (`scripts/picurv`)

- structured errors always use:
  - `ERROR <CODE> | key=... | file=... | message=... | hint=...`
- `run --dry-run --format json` is the safest way to inspect planned commands/artifacts
- `precompute --case ...` writes deterministic grid/profile artifacts without launching solver/post stages
- `run --no-submit` stages run artifacts without starting local execution or Slurm submission; `submit --run-dir` consumes the staged package later
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
- staged run workflow:
  - `./bin/picurv run --solve --post-process ... --no-submit`
  - `./bin/picurv submit --run-dir runs/<run_id>`
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
5. Keep Python function docstrings Doxygen-compatible so `scripts/audit_function_docs.py` stays green.

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
