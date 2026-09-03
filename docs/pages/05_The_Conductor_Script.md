@page 05_The_Conductor_Script The Conductor Script: picurv

@anchor _The_Conductor_Script

`picurv` is the workflow orchestrator for PICurv. It validates YAML inputs, generates C runtime artifacts, and runs or schedules solver/postprocessing stages.
It is also the primary user-facing contract layer: many defaults, aliases, and translation rules are enforced here before the C solver starts.

@tableofcontents

@section p05_usage_sec 1. General Usage

```bash
picurv [COMMAND] [ARGS...]
```

After `make all`, `bin/picurv` is a launcher for the stable `picurv_cli/picurv`
source-tree entrypoint. The implementation lives in `picurv_cli/`. The launcher
uses `.picurv-venv/bin/python` when bootstrap
created a managed environment, then a bootstrap-recorded or user-provided
`PICURV_PYTHON`, and finally `python3`.
Source `etc/picurv.sh` to add `bin/` to your PATH and expose `picurv_cli/` as a
fallback so `picurv` works from any directory even if `bin/picurv` is
temporarily absent before `make conductor` recreates the launcher.
If `bin/picurv` does not exist yet, run `./picurv_cli/picurv build` or `make conductor`.

Primary commands:
- `init`
- `build`
- `version`, `versions`, and `source`
- `inputs`
- `sync-config`
- `run`
- `precompute`
- `storage`
- `summarize`
- `submit`
- `cancel`
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
./picurv_cli/picurv build --help
./bin/picurv version --help
./bin/picurv versions --help
./bin/picurv source --help
./bin/picurv inputs --help
./bin/picurv sync-config --help
./bin/picurv run --help
./bin/picurv precompute --help
./bin/picurv summarize --help
./bin/picurv submit --help
./bin/picurv cancel --help
./bin/picurv sweep --help
./bin/picurv validate --help
```

@section p05_init_sec 2. init: Create A New Case Directory

```bash
picurv init <template_name> [--dest <new_dir>] [--pin-binaries]
```

Behavior:

- copies `examples/<template_name>/` into a new working directory,
- creates the uniform `<workspace>/config/`, `<workspace>/inputs/`, `<workspace>/assets/`,
  `<workspace>/runs/`, and `<workspace>/studies/` tree,
- classifies and relocates template YAML to canonical roles under `<workspace>/config/`,
  preserving colliding variants below `<workspace>/config/variants/`,
- writes `.picurv-workspace.yml`; `software.picurv` is empty by default and therefore
  means the latest active installation,
- optionally renames the destination via `--dest`,
- writes `.picurv-origin.json` with the source repo path and template name,
- writes `.picurv-execution.yml` for optional site-specific launcher overrides,
- seeds that file from a repo-root `.picurv-execution.yml` when the source clone already has one, otherwise from inert defaults,
- does **not** copy binaries by default; runtime executables are resolved from the project `bin/` directory via PATH.

Binary pinning (`--pin-binaries`):

- when `--pin-binaries` is passed, `simulator` and `postprocessor` are copied into the case directory,
- case-local copies take precedence over `bin/` originals at runtime (`resolve_runtime_executable` checks the invocation directory first),
- use this when submitting Slurm jobs and you may rebuild the repo before the job runs,
- `picurv` itself is never copied — it is always used from PATH and is safe to update mid-run since it only launches the C binaries.

Equivalent manual step (after init): `picurv sync-binaries --case-dir <case>`.

Examples:

```bash
picurv init flat_channel --dest my_first_case
picurv init bent_channel --dest my_bent_case --pin-binaries
```

@section p05_build_sec 3. build: Build Project Executables

```bash
./picurv_cli/picurv build [--source-root <repo>] [--case-dir <case>] [MAKE_ARGS...]
```

Behavior:

- calls the top-level `Makefile` via `make`,
- defaults to `make all` when you do not provide an explicit make target,
- resolves the source repo from `.picurv-origin.json` when run from an initialized case,
- passes any trailing arguments directly through to the Make/build layer,
- can rebuild or clean the source repo without leaving a copied case directory,
- writes the streamed build output to `<repo>/logs/build.log` in the source repo,
- is the recommended command for normal users instead of invoking `make` manually.

Direct `make all` keeps the traditional stdout-only behavior. When you want the
same source-repo build plus a warnings-only artifact, use:

```bash
make audit-build
```

This writes:

- `<repo>/logs/build.log`
- `<repo>/logs/build.warnings.log`

Examples:

```bash
./picurv_cli/picurv build
./picurv_cli/picurv build clean-project
./picurv_cli/picurv build SYSTEM=cluster
./picurv_cli/picurv build postprocessor
./my_case/picurv build clean-project
make audit-build
```

@section p05_sync_sec 3b. Source, Version, And Legacy Maintenance Commands

`picurv version` prints the release, Git commit, and dirty-tree build identity shared
by the Python conductor, `simulator --version`, and `postprocessor --version`.
`picurv source update` fetches branches and tags without changing the active checkout.
`picurv versions list` reports tags; `versions install <version>` or
`versions activate [version]` requires a clean source checkout, selects the exact Git
version, and rebuilds. With no argument, `activate` reads the workspace's exact
`software.picurv` value and refuses ranges.

`sync-binaries`, `pull-source`, and `status-source` remain legacy maintenance surfaces.
New workspaces normally use the active installation and the version commands instead
of copying executables into every case.

`sync-binaries` copies `simulator` and `postprocessor` from the source repo `bin/`
into a case directory for version-pinning (optional — normally binaries are resolved via PATH):

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
`sync-config` does not copy `execution.example.yml` into a case; instead it creates
`.picurv-execution.yml` only when the case does not already have one.

`pull-source` refreshes every local branch with a configured upstream, then restores the
branch you started on, so you can update code without leaving the case directory:

```bash
./my_case/picurv pull-source
./my_case/picurv pull-source --current-branch-only
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
- `--monitor <monitor.yml>` (logging, diagnostics, profiling, output cadence, directories)

Inputs for `--post-process`:
- `--post <post.yml>`
- either same invocation with `--solve`, or `--run-dir <existing_run_dir>`

MPI/local options:
- `-n, --num-procs`
  - applies to solver and field postprocessor launch sizing.
  - PICurv strips conflicting MPI size flags from post launches and rewrites them to the requested rank count.
  - local multi-rank runs resolve launcher overrides in this order: `PICURV_MPI_LAUNCHER`, `MPI_LAUNCHER`, nearest `.picurv-execution.yml`, nearest legacy `.picurv-local.yml`, then default `mpiexec`.

Staging/execution options:
- `--cluster <cluster.yml>`
- `--scheduler slurm` (optional explicit selector)
- `--no-submit` (stage run artifacts without starting the local or Slurm backend)

Preflight options:
- `--dry-run` (resolve and print launch/artifact plan only, including PETSc diagnostic flags and log paths)
- `--format json` (machine-readable output for `--dry-run`)

Local example:
```bash
./bin/picurv run --solve --post-process -n 8 \
  --case my_case/case.yml \
  --solver my_case/solver.yml \
  --monitor my_case/monitor.yml \
  --post my_case/post.yml
```
In this command, the solver and field postprocessor both run with 8 ranks.

Slurm example (generate + submit):
```bash
./bin/picurv run --solve --post-process \
  --case my_case/case.yml \
  --solver my_case/solver.yml \
  --monitor my_case/monitor.yml \
  --post my_case/post.yml \
  --cluster my_case/cluster.yml
```

Stage artifacts without starting execution:
```bash
./bin/picurv run --solve --post-process \
  --case my_case/case.yml \
  --solver my_case/solver.yml \
  --monitor my_case/monitor.yml \
  --post my_case/post.yml \
  --no-submit
```
Add `--cluster my_case/cluster.yml` to stage Slurm scripts instead of local commands.

Follow-up execution/submission from existing artifacts:
```bash
./bin/picurv submit --run-dir runs/<run_id>
```

Post-only continuation examples:

- Catch up a live run without editing `post.yml.start_step`. Keep the full desired analysis window in `post.yml`, then let `--continue` move the launch cursor inside that window:

```bash
./bin/picurv run --post-process --continue \
  --run-dir runs/search_robustness_20260322-073415 \
  --post search_robustness_analysis.yml
```

- If the solver has only written source data through step `420`, PICurv launches only the fully available prefix in the requested stride. A later `--continue` run picks up the newer steps after the solver produces them.
- If the same recipe already post-processed the requested window, PICurv skips the launch and reports that the run is already caught up.
- If you change the recipe itself, for example by adding `Qcrit` or changing the statistics output prefix, PICurv treats that as a new recipe lineage and starts again from the configured `start_step`.
- PICurv allows only one post writer per run directory. If a second post job targets the same run, it is refused immediately instead of racing on `<run.visualization>/` or the statistics output beneath it.

Graceful shutdown note:

- generated Slurm solver jobs enable a runtime walltime guard by default; after 10 completed warmup steps, PICurv estimates timestep cost and requests the same graceful final-write path before remaining walltime gets too tight.
- if `cluster.yml -> execution.extra_sbatch.signal` requests an early warning signal, PICurv also traps `SIGUSR1`, `SIGTERM`, and `SIGINT`, then writes one last snapshot at the next safe checkpoint even when the normal recording interval has not been reached.
- tune the automatic estimator in `cluster.yml -> execution.walltime_guard`; keep `execution.extra_sbatch.signal` as fallback protection for preemption/manual termination or jobs that may not reach the warmup window.
- use `signal: "USR1@300"` for `srun`-launched jobs, or `signal: "B:USR1@300"` plus `exec mpirun ...` for direct `mpirun` batch launches.

Runtime stream logs:

- C-managed logs remain under `<run.runtime_logs>/`.
- wrapper stdout/stderr stream logs now live under `<run.scheduler>/` for both local and Slurm launches.
- this avoids collisions with solver startup, which recreates the C log directory.

@section p05_summarize_sec 5. summarize: Read-Only Run Configuration and Health Summary

```bash
./bin/picurv summarize --run-dir <run_dir> \
  [--overview] [--case] [--solver] [--monitor] \
  [--list-plot-series | --plot <qualified-series>] [--last <n>] \
  [--latest | --max-step | --step <n>] [--format json]
```

Behavior:

- `--overview` reports run metadata plus curated case, solver, and monitor summaries,
- `--case`, `--solver`, and `--monitor` are additive selectors for individual copied configs,
- config-only selectors work before runtime logs exist,
- timestep selectors read existing runtime artifacts and still require summary-capable logs,
- plain `summarize --run-dir ...` preserves the implicit latest-step health behavior,
- `--list-plot-series` discovers numeric scalar histories available to `--plot`,
- `--plot <qualified-series>` delegates time-history rendering to `generators/plot.gen`,
- plotting requires `matplotlib`, which bootstrap installs in the managed PICurv environment by default,
- continued runs write `# Continuation from step <N>` separators into append-only logs,
- plots use the latest continuation segment by default; `--last <n>` keeps its last N chronological records per line,
- `--plot-output <path>` saves without opening a window; headless interactive requests save under `<run.analysis>/plots/`,
- builds a best-effort step summary without changing solver output,
- works for active runs and completed runs,
- reports unavailable sections when a source log is missing or disabled.

Examples:

```bash
./bin/picurv summarize --run-dir runs/my_case_20260310-120000 --overview
./bin/picurv summarize --run-dir runs/my_case_20260310-120000 --case --solver
./bin/picurv summarize --run-dir runs/my_case_20260310-120000 --latest
./bin/picurv summarize --run-dir runs/my_case_20260310-120000 --monitor --step 500
./bin/picurv summarize --run-dir runs/my_case_20260310-120000 --list-plot-series
./bin/picurv summarize --run-dir runs/my_case_20260310-120000 --plot momentum.residual_norm --last 100
./bin/picurv summarize --run-dir runs/my_case_20260310-120000 --plot memory.process_peak_mb_max --plot-output peak.png
./bin/picurv summarize --run-dir runs/my_case_20260310-120000 --latest --format json
```

Configuration summaries include configured values plus useful normalized and
derived values such as Reynolds number, nondimensional timestep, normalized
solver selections, effective monitoring cadence, and enabled diagnostics. Text
output uses glanceable dashboards with aligned field groups and compact tables;
use `--format json` for structured automation output.

Plotting supports numeric scalar histories from continuity, particle metrics,
momentum, Poisson, profiling, runtime memory, and solution-convergence logs.
Block-based histories render one line per block and profiling histories render
one line per function. Residual/norm series use logarithmic scaling when all
values are positive; `--linear-y` forces linear scaling. Plot mode is standalone
and does not combine with overview/config/selected-step selectors.

`generators/plot.gen` is the standalone rendering layer. It accepts the versioned,
normalized JSON request produced by `picurv` through stdin, and can also render
a saved request directly:

```bash
python3 generators/plot.gen --input request.json
cat request.json | python3 generators/plot.gen --input -
python3 generators/plot.gen --input request.json --output history.png
```

The normalized request contains plot labels, scale selection, window metadata,
and ordered labeled lines of numeric `[timestep, value]` points. `plot.gen`
does not parse PICurv run directories or logs.

Typical sources:

- `<run.runtime_logs>/Continuity_Metrics.log`
- `<run.runtime_logs>/Particle_Metrics.log` (per-step `Lost` plus run-local `Lost Total` when available)
- `<run.runtime_logs>/Momentum_Solver_Convergence_History_Block_*.log`
- `<run.runtime_logs>/Poisson_Solver_Convergence_History_Block_*.log`
- `<run.runtime_logs>/solution_convergence.log` (mode-specific speed/KE drift and L2 norms)
- `<run.runtime_logs>/Profiling_Timestep_Summary.csv` when enabled
- `<run.runtime_logs>/Runtime_Memory.log` when `monitor.diagnostics.runtime_memory_log.enabled` is true
- `<run.runtime_logs>/les_coefficient.csv` when `case.yml -> models.physics.turbulence.les.diagnostics.enabled` is true; the `les.*` series carry the effective coefficient, its spread, eddy-viscosity levels, the modelled subgrid energy, and the pre-clipping backscattering and limited fractions
- `<run.runtime_logs>/PETSc_*_Solver.log` / `<run.runtime_logs>/PETSc_*_PostProcessor.log` when file-backed PETSc diagnostics are enabled
- `<run.scheduler>/*_solver.log` or `<run.scheduler>/solver_*.out` for sampled particle snapshot previews

When particle snapshots are available, `summarize` reports sampled diagnostics such as:

- speed min/mean/max/std and stagnant-count
- sampled position bounds and centroid
- sampled rank counts and duplicate sampled cells
- sampled weight min/max by component
- sanity checks for duplicate PIDs and non-finite/zero/negative weights
- top sampled speeds and sampled deltas versus the previous snapshot when matching PIDs exist

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
- staged local/cluster run without execution: `--no-submit`
- delayed submission of an already-staged run: `submit --run-dir ...`
- planning and CI-style checks: `--dry-run`

@section p05_precompute_sec 5a. precompute: Generate Deterministic Artifacts

```bash
picurv precompute --case config/case.yml [--only grid,initial-condition,inlet-profiles]
```

`precompute` creates immutable, content-addressed workspace assets without launching
the solver. It uses exactly the same provider graph as `run --solve`: grids build
before any Python initial condition or inlet slice that consumes them. Provider
settings, referenced-file checksums, and software identity determine reuse.

The build occurs under `assets/.precompute-*` and publishes to `assets/objects/` only
after the whole requested dependency closure succeeds. The mutable pointer in
`assets/sets/` records which objects match the current case. A normal run reuses those
objects automatically and writes the exact selection to `inputs/assets.lock.yml`.

Providers owned by the C runtime are reported but never approximated by Python. A
selected C-only provider makes precompute fail before publication. See
**@subpage 52_Run_Artifact_Lifecycle_Contract** for the directory and dependency contracts.

@section p05_submit_sec 5b. submit: Execute Existing Staged Artifacts

```bash
./bin/picurv submit [--run-dir <run_dir> | --study-dir <study_dir>] \
  [--stage {all,solve,post-process}] [--force] [--dry-run]
```

Behavior:

- consumes existing `--no-submit` artifacts without regenerating configs or scripts,
- reads `<run.scheduler>/submission.json` to locate staged Slurm scripts or local command tokens,
- executes/submits `solve`, `post-process`, or both,
- wires the post stage dependency automatically for Slurm when `all` is selected,
- executes local staged stages in order when `launch_mode: local`,
- refuses re-submission unless `--force` is explicitly provided.

Examples:

```bash
./bin/picurv submit --run-dir runs/my_case_20260310-120000
./bin/picurv submit --run-dir runs/my_case_20260310-120000 --stage solve
./bin/picurv submit --study-dir studies/my_study_20260310-120000 --dry-run
```

Notes:

- `--run-dir` supports local and Slurm staged runs,
- `--study-dir` remains Slurm-only,
- `--dry-run` prints the exact local command or `sbatch` plan,
- `--force` is the opt-in path for deliberate resubmission.

@section p05_cancel_sec 5c. cancel: Stop A Slurm Run By Run Directory

```bash
./bin/picurv cancel --run-dir <run_dir> [--stage {all,solve,post-process}] [--graceful] [--dry-run]
```

Behavior:

- reads `<run.scheduler>/submission.json` from an existing run directory,
- resolves the recorded Slurm job IDs for `solve` and/or `post-process`,
- runs hard `scancel <job_id>` for the selected stage set by default,
- with `--graceful`, sends `scancel --signal=USR1 --full <job_id>` for solver jobs, reaching the batch process and MPI children so the runtime can write the latest safe off-cadence step at the next checkpoint,
- post-process jobs still use hard cancellation, even when `--graceful` is present,
- avoids manual job-id lookup when the run directory is already known.

Examples:

```bash
./bin/picurv cancel --run-dir runs/my_case_20260310-120000
./bin/picurv cancel --run-dir runs/my_case_20260310-120000 --stage solve
./bin/picurv cancel --run-dir runs/my_case_20260310-120000 --stage solve --graceful
./bin/picurv cancel --run-dir runs/my_case_20260310-120000 --dry-run
```

Notes:

- works only for Slurm-submitted runs that have `<run.scheduler>/submission.json`,
- does not apply to local runs,
- `--graceful` requests solver shutdown and final output; if a job is wedged or not reaching runtime checkpoints, rerun without `--graceful` to hard-cancel it,
- `--dry-run` is useful when you want to confirm the recorded stage/job mapping first.

@section p05_sweep_sec 6. sweep: Parameter Study via Slurm Arrays

```bash
# Launch new study
./bin/picurv sweep \
  --study my_study/study.yml \
  --cluster my_study/cluster.yml [--no-submit]

# Continue a partially-completed study
./bin/picurv sweep --continue --study-dir studies/<study_id> \
  [--cluster cluster_more_time.yml]

# Re-aggregate metrics manually
./bin/picurv sweep --reaggregate --study-dir studies/<study_id>
```

Behavior (new study):
- expands parameter matrix from `study.yml`
- materializes case directories under `studies/<study_id>/cases/`
- generates `solver_array.sbatch`, `post_array.sbatch`, and `metrics_aggregate.sbatch`
- renders `post_array.sbatch` with the same cluster allocation as the solver array
- submits solver → post (`afterok`) → metrics (`afterany`) chain (unless `--no-submit`)

Behavior (`--continue`):
- detects per-case completion status (complete / partial / empty)
- if all cases complete, auto-aggregates metrics and exits
- otherwise prepares continuation for incomplete cases (checkpoint restart via `resolve_restart_source`)
- submits sparse solver array (incomplete cases only) → full post array → metrics aggregation

Behavior (`--reaggregate`):
- re-runs metrics collection and plot generation on existing study outputs

@section p05_validate_sec 7. validate: Config-Only Checks

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

@section p05_command_matrix_sec 8. Full Command and Option Matrix

The table below is **generated by introspecting the assembled `argparse` parser** -
`build_main_parser()` in `picurv_cli/cli.py` together with the registrars delegated from
other modules, such as `add_storage_parser()` in `picurv_cli/storage.py`. Reading
`cli.py` alone understates the command set, which is why the generator exercises the
built parser rather than pattern-matching the source.

It is regenerated by `make docs-cli-reference` and checked by `make audit-cli-reference`,
so a parser change makes it stale and fails CI. The narrative guidance in the sections
above is hand-written and complements it; this section is the exhaustive reference.

@htmlinclude generated/cli_reference.html

`run`:
- stages:
  - `--solve` (requires `--case`, `--solver`, `--monitor`)
  - `--post-process` (requires `--post`; requires `--run-dir` when `--solve` is not selected)
- solver inputs:
  - `--case <path>`
  - `--solver <path>`
  - `--monitor <path>` (logging, diagnostics, profiling, output cadence)
- post inputs:
  - `--post <path>`
  - `--run-dir <path>`
- restart and asset selection:
  - `--restart-from <run-dir|latest>` (alias `--from`; `latest` picks the newest
    compatible workspace run)
  - `--statistics-state {reset,carry}` (branched restart only; entries at
    @ref p52_cap_restart_stats_sec)
  - `--require-precomputed` (refuse to build a missing or stale workspace asset)
  - `--fetch-missing` (try the configured storage profile before rebuilding one)
- launch controls:
  - `-n, --num-procs <int>` (solver and field postprocessor stages)
  - `--cluster <cluster.yml>` (enables Slurm mode)
  - `--scheduler <name>` (must be used with `--cluster`; must match `cluster.yml:scheduler.type`)
  - `--no-submit` (stage files and submission metadata without starting execution)
  - `--dry-run` (no file writes; plan only, including diagnostic artifacts)
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

`precompute`:
- required:
  - `--case <path>`
- optional:
  - `--only <grid,initial-condition,inlet-profiles>` (defaults to all configured providers)

`summarize`:
- required:
  - `--run-dir <path>`
- <run.config>/run views:
  - `--overview`
  - `--case`
  - `--solver`
  - `--monitor`
- selected-step health:
  - mutually exclusive `--step <n>`, `--latest`, or `--max-step`
  - `--snapshot-rows <positive-int>`
- standalone time-history discovery/plotting:
  - mutually exclusive `--list-plot-series` or `--plot <qualified-series>`
  - `--last <positive-int>` (requires `--plot`)
  - `--plot-output <path>` (requires `--plot`)
  - `--linear-y` (requires `--plot`)
  - plot/list modes cannot combine with <run.config>/run views or selected-step health
- output:
  - `--format {text,json}` (`--list-plot-series` supports JSON; `--plot` rejects JSON)

`submit`:
- required:
  - exactly one of `--run-dir <path>` or `--study-dir <path>`
- optional:
  - `--stage {all,solve,post-process}`
  - `--force`
  - `--dry-run`

`cancel`:
- required:
  - `--run-dir <path>`
- optional:
  - `--stage {all,solve,post-process}`
  - `--dry-run`

`sweep`:
- new study mode (default):
  - `--study <study.yml>` (required)
  - `--cluster <cluster.yml>` (required)
  - `--no-submit` (optional)
- continuation mode:
  - `--continue` (required)
  - `--study-dir <path>` (required)
  - `--cluster <cluster.yml>` (optional; overrides original cluster resources)
- reaggregation mode:
  - `--reaggregate` (required)
  - `--study-dir <path>` (required)

`init`:
- positional:
  - `<template_name>`
- optional:
  - `--dest <dir>`
  - `--source-root <repo>`
  - `--pin-binaries` (copy `simulator`/`postprocessor` into the case for version-pinning)

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

For file-backed grid modes, `artifacts` includes the planned staged grid path
(`<run.config>/grid.run`). For `grid.mode: grid_gen`, it also includes the generated
PICGRID path (`<run.config>/grid.generated.picgrid` by default) plus any configured
`stats_file` or `vts_file`. Dry-run still does not run `generators/grid.gen` or
write these files.
For generated prescribed-flow profiles, `artifacts` includes the dimensional
generated `.picslice`, the solver-scale staged `.picslice`, and `profile.info`.

Stage entries under `stages.solve` and `stages.post-process` include:

- `mode` (`local` or `slurm`)
- `num_procs_effective`
- `launch_command` (tokenized command list)
- `launch_command_string` (shell-ready display string)

Additional stage fields:
- `script` (Slurm script path in cluster mode)
- `command` / `command_string` (local staged command tokens/display string)
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
3. keep `monitor.yml` focused on logging, diagnostics, profiling, and I/O cadence,
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

@section p05_binaries_sec 12. Binary Resolution and Rebuild Safety

`picurv` resolves `simulator` and `postprocessor` at launch time using this precedence:

1. **Invocation directory** — if the binary exists as a sibling of the invoked `picurv` script
   (e.g. case-local copies from `--pin-binaries` or `sync-binaries`), it is used first.
2. **Project `bin/` directory** — the default location after `make all`.

`bin/picurv` is a launcher for `picurv_cli/picurv`. This means:

- there is one stable executable entrypoint backed by the `picurv_cli/` package,
- `make conductor` recreates the launcher (idempotent),
- bootstrap-managed installs run with the repo-local `.picurv-venv` Python,
- `etc/picurv.sh` adds `bin/` to PATH so `picurv` works from any directory.

**Rebuilding while jobs are running:**

- Updating `picurv` (the Python script) mid-run is always safe — it is only used to launch jobs,
  not during solver execution.
- Rebuilding `simulator`/`postprocessor` (`make all`) overwrites the binaries in `bin/`.
  If a Slurm job references `bin/simulator` by absolute path and has not yet started, the running
  binary may be replaced before execution begins.
- For a queued production job, use a tagged release checkout or an explicitly pinned
  legacy binary copy. The run manifest records the build identity that staged it;
  `picurv version` checks the currently active identity.

**Recommended workflow for concurrent development and production:**

```bash
picurv init flat_channel --dest production_case --pin-binaries
picurv run --solve --cluster cluster.yml ...   # uses case-local binaries
# safe to rebuild in the repo now — production_case has its own copies
make all
```

@section p05_cap_input_mode_sec 12.1 Workspace Input Import Mode Entries

@htmlinclude generated/capability_inventory_workspace_input_import_mode.html

@subsection p05_cap_input_mode_copy_sub copy

@anchor p05_cap_input_mode_copy

**Identity.** `picurv inputs import <kind> <source> --mode copy`.

**What it does.** Copies the source into its canonical workspace `inputs/` home and
catalogues the source checksum, size, and new relative path.

**When to choose it.** Default to `copy` when portability and independent ownership
matter more than local disk duplication. Unlike links, the workspace remains valid if
the source is moved or deleted.

**Parameters it owns.** Optional `--name` supplies the destination basename.

**Interactions.** Storage includes the copied file when a run locks and materializes
an asset derived from it.

**Diagnostics.** Success prints the catalog ID and relative path; an existing target
or missing source fails before the catalog changes.

**Evidence.** Unit verified — `tests/test_workspace_lifecycle.py` checks the copied
bytes and catalog entry.

**Limitations.** Uses additional disk space equal to the imported file.

@subsection p05_cap_input_mode_reflink_sub reflink

@anchor p05_cap_input_mode_reflink

**Identity.** `picurv inputs import <kind> <source> --mode reflink`.

**What it does.** Requests a copy-on-write clone through `cp --reflink=always`, then
catalogues it like a copy.

**When to choose it.** Prefer it for very large local inputs on a filesystem that
supports reflinks: it begins space-efficiently but remains independently writable,
unlike `hardlink`.

**Parameters it owns.** Optional `--name` supplies the destination basename.

**Interactions.** Failure does not fall back silently; choose `copy` explicitly when
the filesystem reports that reflinks are unsupported.

**Diagnostics.** The import fails with the native copy error if reflink creation is
unavailable and writes no catalog record.

**Evidence.** Unit verified — `tests/test_workspace_lifecycle.py` verifies the mode is
part of the public parser and catalog contract.

**Limitations.** Experimental across cluster filesystems; support depends on the local
`cp` and filesystem.

@subsection p05_cap_input_mode_hardlink_sub hardlink

@anchor p05_cap_input_mode_hardlink

**Identity.** `picurv inputs import <kind> <source> --mode hardlink`.

**What it does.** Creates a second directory entry for the same inode and catalogues
the imported path.

**When to choose it.** Use it only for immutable, same-filesystem inputs when avoiding
a full copy matters and every owner agrees not to modify the bytes.

**Parameters it owns.** Optional `--name` supplies the destination basename.

**Interactions.** Source and destination share mutations. Asset checksums detect a
changed provider input and prevent stale object reuse, but cannot undo the mutation.

**Diagnostics.** Cross-filesystem or permission failures are reported before the
catalog changes.

**Evidence.** Unit verified — `tests/test_workspace_lifecycle.py` verifies the mode is
part of the public parser and catalog contract.

**Limitations.** Experimental because shared-inode ownership is easy to misuse and it
cannot cross filesystems.

@subsection p05_cap_input_mode_reference_sub reference

@anchor p05_cap_input_mode_reference

**Identity.** `picurv inputs import <kind> <source> --mode reference`.

**What it does.** Writes a small `.reference.yml` containing the absolute target,
registration checksum, and byte count; it does not copy the target.

**When to choose it.** Use it for an authoritative shared dataset that should remain
externally owned. Choose `copy` when the workspace or storage archive must be portable.

**Parameters it owns.** Optional `--name` names the reference record.

**Interactions.** Resolution checks the target loudly. Storage records the dependency
but neither archives nor prunes the external file.

**Diagnostics.** Registration prints an external-reference warning; missing targets
fail when registered or consumed.

**Evidence.** Unit verified — `tests/test_workspace_lifecycle.py` checks reference
content and catalog identity.

**Limitations.** Restoring the workspace cannot restore an external target; its owner
must make the same path available or the reference must be replaced.

@section p05_artifacts_sec 13. Generated Runtime Artifacts

Single run (`run`):
- `<run.config>/*.control`, `bcs*.run`, immutable YAML snapshots, plus optional
  `whitelist.run` / `profile.run` sidecars when enabled
- `<run.config.history>/<revision>/` (new active config snapshots for continuation)
- `<run.post_recipes>/<recipe-id>/{post.yml,post.run,state.json}`
- `<run.inputs>/{grid,initial_condition,inlet_profiles,restart}/` and `<run.asset_lock>`
- `<run.checkpoints>/`, `<run.analysis>/`, and `<run.visualization>/`
- `<run.runtime_logs>/` and `<run.scheduler>/` (scripts, stream logs, and submission state)
- `runs/<run_id>/manifest.json`

Sweep (`sweep`):
- `studies/<study_id>/cases/<case_i>/...`
- `studies/<study_id>/scheduler/case_index.tsv`
- `studies/<study_id>/scheduler/solver_array.sbatch`
- `studies/<study_id>/scheduler/post_array.sbatch`
- `studies/<study_id>/scheduler/solver_<array_jobid>_<taskid>.out/.err`, `post_<array_jobid>_<taskid>.out/.err`
- `studies/<study_id>/scheduler/submission.json`
- `studies/<study_id>/output/analysis/metrics_table.csv`
- `studies/<study_id>/output/analysis/plots/*.png` (if plotting enabled and matplotlib available)
- `studies/<study_id>/study_manifest.json`

@section p05_storage_sec 14. `storage`

`storage` is a top-level command registered from `picurv_cli/storage.py`, and is the
one command family this page does not detail. Subcommands: `setup`, `plan`, `protect`,
`offload`, `restore`, `verify`, `list`, `show`, `status`.

Full semantics, retention model, and rclone configuration are documented in
**@subpage 61_Storage_Management_Guide**.

@section p05_next_steps_sec 15. Next Steps

- Config contract: **@subpage 14_Config_Contract**
- User workflows: **@subpage 11_User_How_To_Guides**
- Extensibility: **@subpage 17_Workflow_Extensibility**
- First-run path: **@subpage 02_Tutorial_Programmatic_Grid**
- Grid generator details: **@subpage 48_Grid_Generator_Guide**
- Modular examples and recipes: **@subpage 49_Workflow_Recipes_and_Config_Cookbook**
- Troubleshooting map: **@subpage 39_Common_Fatal_Errors**
- Smoke/CI verification details: **@subpage 40_Testing_and_Quality_Guide**
