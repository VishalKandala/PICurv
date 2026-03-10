@page 05_The_Conductor_Script The Conductor Script: picurv

@anchor _The_Conductor_Script

`picurv` is the workflow orchestrator for PICurv. It validates YAML inputs, generates C runtime artifacts, and runs or schedules solver/postprocessing stages.
It is also the primary user-facing contract layer: many defaults, aliases, and translation rules are enforced here before the C solver starts.

@tableofcontents

@section p05_usage_sec 1. General Usage

```bash
./bin/picurv [COMMAND] [ARGS...]
```

If `bin/picurv` does not exist yet, use `./scripts/picurv build` once to install it.

Primary commands:
- `init`
- `build`
- `sync-binaries`
- `sync-config`
- `pull-source`
- `status-source`
- `run`
- `sweep`
- `validate`

Core idea:

- `case.yml`, `solver.yml`, `monitor.yml`, and `post.yml` are modular profiles.
- You do not need to rewrite all four every time.
- In normal use, you mix and match a case definition with reusable solver, monitor, and post profiles as long as the combination is contract-compatible.
- Example: the same `case.yml` can be paired with a quieter `monitor.yml`, a different solver strategy, or a different post recipe without changing the physical setup file.

Help:
```bash
./bin/picurv --help
./bin/picurv init --help
./scripts/picurv build --help
./bin/picurv sync-binaries --help
./bin/picurv sync-config --help
./bin/picurv pull-source --help
./bin/picurv status-source --help
./bin/picurv run --help
./bin/picurv sweep --help
./bin/picurv validate --help
```

@section p05_init_sec 2. init: Create A New Case Directory

```bash
./bin/picurv init <template_name> [--dest <new_dir>]
```

Behavior:

- copies `examples/<template_name>/` into a new working directory,
- optionally renames the destination via `--dest`,
- writes `.picurv-origin.json` with the source repo path and template name,
- copies the full executable set from `bin/` into the new case directory so the case is self-contained,
- includes `picurv`, `simulator`, and `postprocessor` when they exist in `bin/`,
- copies built executables when they exist in the source repo `bin/`.

Examples:

```bash
./bin/picurv init flat_channel --dest my_first_case
./bin/picurv init bent_channel --dest my_bent_case
```

Use `init` when you want a runnable starting point that you can `cd` into and launch with `./picurv ...`.
Use reusable profile libraries directly when you already have your own case directory layout.

@section p05_build_sec 3. build: Build Project Executables

```bash
./scripts/picurv build [--source-root <repo>] [--case-dir <case>] [MAKE_ARGS...]
```

Behavior:

- calls the top-level `Makefile` via `make`,
- resolves the source repo from `.picurv-origin.json` when run from an initialized case,
- passes any trailing arguments directly through to the Make/build layer,
- can rebuild or clean the source repo without leaving a copied case directory,
- is the recommended command for normal users instead of invoking `make` manually.

Examples:

```bash
./scripts/picurv build
./scripts/picurv build clean-project
./scripts/picurv build SYSTEM=cluster
./scripts/picurv build postprocessor
./my_case/picurv build clean-project
```

@section p05_sync_sec 3b. sync-binaries / sync-config / pull-source / status-source

These commands are intended for self-contained case directories created by `init`.

`sync-binaries` refreshes the case-local copies of `picurv`, `simulator`, and `postprocessor`
from the source repo `bin/`:

```bash
./my_case/picurv sync-binaries
```

`sync-config` refreshes files from `examples/<template_name>/` into the case directory.
By default it preserves user-modified files and only copies missing files:

```bash
./my_case/picurv sync-config
./my_case/picurv sync-config --overwrite
./my_case/picurv sync-config --prune
```

`--prune` is conservative: it removes only files previously recorded as template-managed
that no longer exist in the source template. User-created case files are not pruned.

`pull-source` runs `git pull --rebase` in the original source repo so you can update code
without leaving the case directory:

```bash
./my_case/picurv pull-source
./my_case/picurv pull-source --no-rebase
./my_case/picurv pull-source --remote origin --branch main
```

`status-source` inspects source commit drift, copied binary drift, and template-file drift
before you decide what to sync:

```bash
./my_case/picurv status-source
./my_case/picurv status-source --format json
```

For older cases that do not yet have `.picurv-origin.json`, pass `--source-root /path/to/PICurv`.
For `sync-config`, also pass `--template-name <example_name>` if the template cannot be inferred.

`status-source --format json` payload highlights:

- `source_repo_root`, `case_dir`
- `last_known_source_git_commit`, `current_source_git_commit`, `source_commit_changed`
- `binaries`:
  - `case_bin_current`, `case_bin_different`, `case_bin_missing`
- `config`:
  - `case_current_files`, `case_modified_files`, `case_missing_files`
  - `template_removed_since_last_sync` (when template-managed tracking is available)

@section p05_run_sec 4. run: Single-Case Workflow

```bash
./bin/picurv run [STAGES] [INPUTS] [OPTIONS]
```

Stages:
- `--solve`
- `--post-process`

Inputs for `--solve`:
- `--case <case.yml>`
- `--solver <solver.yml>`
- `--monitor <monitor.yml>`

Inputs for `--post-process`:
- `--post <post.yml>`
- either same invocation with `--solve`, or `--run-dir <existing_run_dir>`

MPI/local options:
- `-n, --num-procs`
  - applies to solver stage launch sizing.
  - post-processing defaults to one rank/task.
  - local multi-rank runs resolve launcher overrides in this order: `PICURV_MPI_LAUNCHER`, `MPI_LAUNCHER`, nearest `.picurv-execution.yml`, nearest legacy `.picurv-local.yml`, then default `mpiexec`.

Cluster/Slurm options:
- `--cluster <cluster.yml>`
- `--scheduler slurm` (optional explicit selector)
- `--no-submit` (generate scripts/manifests only)

Preflight options:
- `--dry-run` (resolve and print launch/artifact plan only)
- `--format json` (machine-readable output for `--dry-run`)

Local example:
```bash
./bin/picurv run --solve --post-process -n 8 \
  --case my_case/case.yml \
  --solver my_case/solver.yml \
  --monitor my_case/monitor.yml \
  --post my_case/post.yml
```
In this command, solver runs with 8 ranks; post-processing still defaults to 1 rank.

Slurm example (generate + submit):
```bash
./bin/picurv run --solve --post-process \
  --case my_case/case.yml \
  --solver my_case/solver.yml \
  --monitor my_case/monitor.yml \
  --post my_case/post.yml \
  --cluster my_case/cluster.yml
```

Slurm example (generate only):
```bash
./bin/picurv run --solve --post-process \
  --case my_case/case.yml \
  --solver my_case/solver.yml \
  --monitor my_case/monitor.yml \
  --post my_case/post.yml \
  --cluster my_case/cluster.yml \
  --no-submit
```

Dry-run example (no file writes):
```bash
./bin/picurv run --solve --post-process \
  --case my_case/case.yml \
  --solver my_case/solver.yml \
  --monitor my_case/monitor.yml \
  --post my_case/post.yml \
  --dry-run --format json
```

Common `run` use cases:

- first full local run: `--solve --post-process`
- solver-only run: `--solve`
- post-only rerun on an existing run directory: `--post-process --run-dir ...`
- cluster script generation without submit: `--cluster ... --no-submit`
- planning and CI-style checks: `--dry-run`

@section p05_sweep_sec 5. sweep: Parameter Study via Slurm Arrays

```bash
./bin/picurv sweep \
  --study my_study/study.yml \
  --cluster my_study/cluster.yml [--no-submit]
```

Behavior:
- expands parameter matrix from `study.yml`
- materializes case directories under `studies/<study_id>/cases/`
- generates `solver_array.sbatch` and `post_array.sbatch`
- submits post array with `afterok:<solver_jobid>` dependency (unless `--no-submit`)
- aggregates metrics and emits plots in `studies/<study_id>/results/`

@section p05_validate_sec 6. validate: Config-Only Checks

```bash
./bin/picurv validate \
  --case my_case/case.yml \
  --solver my_case/solver.yml \
  --monitor my_case/monitor.yml \
  --post my_case/post.yml
```

`validate` does not launch solver/post and does not create run/study artifacts.

What `validate` is for:

- check a new profile combination before running,
- confirm a modified template still satisfies the current schema,
- catch mode-dependent contract errors before the C runtime,
- inspect warnings where `picurv` preserves a C-side default intentionally.

@section p05_command_matrix_sec 7. Full Command and Option Matrix

This section is intentionally exhaustive and mirrors the current `argparse` contract in `scripts/picurv`.
Use it as the authoritative option reference when writing docs, examples, wrappers, or CI jobs.

`run`:
- stages:
  - `--solve` (requires `--case`, `--solver`, `--monitor`)
  - `--post-process` (requires `--post`; requires `--run-dir` when `--solve` is not selected)
- solver inputs:
  - `--case <path>`
  - `--solver <path>`
  - `--monitor <path>`
- post inputs:
  - `--post <path>`
  - `--run-dir <path>`
- launch controls:
  - `-n, --num-procs <int>` (solver stage only; post remains 1 rank/task)
  - `--cluster <cluster.yml>` (enables Slurm mode)
  - `--scheduler <name>` (must be used with `--cluster`; must match `cluster.yml:scheduler.type`)
  - `--no-submit` (render scripts/manifests without `sbatch`)
  - `--dry-run` (no file writes; plan only)
  - `--format {text,json}` (dry-run output format)

`validate`:
- role selectors:
  - `--case <path>`
  - `--solver <path>`
  - `--monitor <path>`
  - `--post <path>`
  - `--cluster <path>`
  - `--study <path>`
- stricter policy:
  - `--strict` (adds additional checks for selected roles; documented below)

`sweep`:
- required:
  - `--study <study.yml>`
  - `--cluster <cluster.yml>`
- optional:
  - `--no-submit`

`init`:
- positional:
  - `<template_name>`
- optional:
  - `--dest <dir>`
  - `--source-root <repo>`

`build`:
- optional:
  - `--source-root <repo>`
  - `--case-dir <case_dir>`
- passthrough:
  - trailing `MAKE_ARGS...` are passed directly to `make`

`sync-binaries`:
- `--case-dir <case_dir>`
- `--source-root <repo>`

`sync-config`:
- `--case-dir <case_dir>`
- `--source-root <repo>`
- `--template-name <template>`
- `--overwrite`
- `--prune`

`pull-source`:
- `--case-dir <case_dir>`
- `--source-root <repo>`
- `--remote <git_remote>`
- `--branch <git_branch>`
- `--no-rebase`

`status-source`:
- `--case-dir <case_dir>`
- `--source-root <repo>`
- `--template-name <template>`
- `--format {text,json}`

@section p05_strict_validate_sec 8. `validate --strict`: Additional Checks

`--strict` does not change baseline schema validation, but it adds file-system consistency checks for selected roles:

- with `--post`:
  - if `source_data.directory` is not `<solver_output_dir>`, the resolved directory must exist.
- with `--study`:
  - base configs listed in `study.base_configs` are loaded and revalidated as real case/solver/monitor/post bundles.
  - this catches study files that are syntactically valid but point to invalid base configurations.

Use `--strict` in CI/pre-submit checks when validating reusable profile libraries or study manifests.

@section p05_dry_run_schema_sec 9. Dry-Run JSON Plan Schema

`picurv run --dry-run --format json` emits a deterministic plan payload with these top-level keys:

- `mode` (`"dry-run"`)
- `created_at` (ISO timestamp)
- `launch_mode` (`local` or `slurm`)
- `warnings` (list)
- `inputs` (resolved absolute paths)
- `stages` (stage-specific launch plans)
- `artifacts` (predicted files/directories, deduplicated)
- `run_id_preview` / `run_dir_preview` (when known)
- `solver_num_procs_effective`
- `post_num_procs_effective`
- `num_procs_effective` (currently mirrors solver count)

Stage entries under `stages.solve` and `stages.post-process` include:

- `mode` (`local` or `slurm`)
- `num_procs_effective`
- `launch_command` (tokenized command list)
- `launch_command_string` (shell-ready display string)

Additional stage fields:
- `script` (Slurm script path in cluster mode)
- `source_data_directory` (post stage source directory resolution)
- `restart_source_directory` (solve stage, when restart source is resolved)

Dry-run guarantees:
- no run directory creation,
- no control/post recipe writes,
- no scheduler script writes,
- no job submission.

@section p05_structured_errors_sec 10. Structured Error Output Contract

Validation and CLI usage errors are emitted as one-line, machine-parseable records:

```text
ERROR <CODE> | key=<yaml_or_cli_key> | file=<path_or_dash> | message=<summary> | hint=<actionable_hint>
```

Current normalized error code set:

- `CLI_USAGE_INVALID`
- `CFG_MISSING_SECTION`
- `CFG_MISSING_KEY`
- `CFG_INVALID_TYPE`
- `CFG_INVALID_VALUE`
- `CFG_FILE_NOT_FOUND`
- `CFG_GRID_PARSE`
- `CFG_INCONSISTENT_COMBO`

This contract is exercised by Python tests and should remain stable for wrappers and CI parsers.

@section p05_modular_sec 11. Modular Profile Strategy

PICurv is intended to be used with reusable profile libraries.

Typical pattern:

1. keep `case.yml` focused on physics, grid, BCs, and run duration,
2. keep `solver.yml` focused on numerical strategy,
3. keep `monitor.yml` focused on logging and I/O cadence,
4. keep `post.yml` focused on analysis outputs,
5. recombine them as needed for different studies.

Examples:

- same `case.yml` + a lighter `monitor.yml` for fast debug runs,
- same `case.yml` + a stricter `solver.yml` for convergence checks,
- same `case.yml` + multiple `post.yml` recipes for different analysis outputs,
- same `solver.yml` reused across many cases when the discretization strategy is stable.

This is why `picurv` treats these roles as separate inputs instead of one monolithic file.

For prebuilt reusable profiles, also see the local guides under:

- `config/solvers/`
- `config/monitors/`
- `config/postprocessors/`
- `config/schedulers/`
- `examples/master_template/`

@section p05_artifacts_sec 12. Generated Runtime Artifacts

Single run (`run`):
- `runs/<run_id>/config/*.control`, `bcs*.run`, `post.run`, plus optional `whitelist.run` / `profile.run` sidecars when enabled
- `runs/<run_id>/scheduler/solver.sbatch`, `post.sbatch` (cluster mode)
- `runs/<run_id>/scheduler/submission.json` (cluster mode)
- `runs/<run_id>/manifest.json`

Sweep (`sweep`):
- `studies/<study_id>/cases/<case_i>/...`
- `studies/<study_id>/scheduler/case_index.tsv`
- `studies/<study_id>/scheduler/solver_array.sbatch`
- `studies/<study_id>/scheduler/post_array.sbatch`
- `studies/<study_id>/scheduler/submission.json`
- `studies/<study_id>/results/metrics_table.csv`
- `studies/<study_id>/results/plots/*.png` (if plotting enabled and matplotlib available)
- `studies/<study_id>/study_manifest.json`

@section p05_next_steps_sec 13. Next Steps

- Config contract: **@subpage 14_Config_Contract**
- User workflows: **@subpage 11_User_How_To_Guides**
- Extensibility: **@subpage 17_Workflow_Extensibility**
- First-run path: **@subpage 02_Tutorial_Programmatic_Grid**
- Grid generator details: **@subpage 48_Grid_Generator_Guide**
- Modular examples and recipes: **@subpage 49_Workflow_Recipes_and_Config_Cookbook**
- Troubleshooting map: **@subpage 39_Common_Fatal_Errors**
- Smoke/CI verification details: **@subpage 40_Testing_and_Quality_Guide**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **The Conductor Script: picurv** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

Treat this page as both a conceptual reference and a runbook. If you are debugging, pair the method/procedure described here with monitor output, generated runtime artifacts under `runs/<run_id>/config`, and the associated solver/post logs so numerical intent and implementation behavior stay aligned.

### What To Extract Before Changing A Case

- Identify which YAML role or runtime stage this page governs.
- List the primary control knobs (tolerances, cadence, paths, selectors, or mode flags).
- Record expected success indicators (convergence trend, artifact presence, or stable derived metrics).
- Record failure signals that require rollback or parameter isolation.

### Practical CFD Troubleshooting Pattern

1. Reproduce the issue on a tiny case or narrow timestep window.
2. Change one control at a time and keep all other roles/configs fixed.
3. Validate generated artifacts and logs after each change before scaling up.
4. If behavior remains inconsistent, compare against a known-good baseline example and re-check grid/BC consistency.
