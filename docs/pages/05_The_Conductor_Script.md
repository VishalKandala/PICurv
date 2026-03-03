@page 05_The_Conductor_Script The Conductor Script: picurv

`picurv` is the workflow orchestrator for PICurv. It validates YAML inputs, generates C runtime artifacts, and runs or schedules solver/postprocessing stages.
It is also the primary user-facing contract layer: many defaults, aliases, and translation rules are enforced here before the C solver starts.

@tableofcontents

@section usage_sec 1. General Usage

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

@section init_sec 2. init: Create A New Case Directory

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

@section build_sec 3. build: Build Project Executables

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

@section sync_sec 3b. sync-binaries / sync-config / pull-source / status-source

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

@section run_sec 4. run: Single-Case Workflow

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

@section sweep_sec 5. sweep: Parameter Study via Slurm Arrays

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

@section validate_sec 6. validate: Config-Only Checks

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

@section modular_sec 7. Modular Profile Strategy

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

@section artifacts_sec 8. Generated Runtime Artifacts

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

@section next_steps_sec 9. Next Steps

- Config contract: **@subpage 14_Config_Contract**
- User workflows: **@subpage 11_User_How_To_Guides**
- Extensibility: **@subpage 17_Workflow_Extensibility**
- First-run path: **@subpage 02_Tutorial_Programmatic_Grid**
- Grid generator details: **@subpage 48_Grid_Generator_Guide**
- Modular examples and recipes: **@subpage 49_Workflow_Recipes_and_Config_Cookbook**
- Troubleshooting map: **@subpage 39_Common_Fatal_Errors**
- Smoke/CI verification details: **@subpage 40_Testing_and_Quality_Guide**
