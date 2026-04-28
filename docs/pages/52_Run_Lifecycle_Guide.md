@page 52_Run_Lifecycle_Guide Run Lifecycle Guide

@anchor _Run_Lifecycle_Guide

This page explains how a PICurv run moves from a new solve to restart, post-processing reuse, and cluster job generation.
It is the operational view of the run directory lifecycle.

@tableofcontents

@section p52_scope_sec 1. What A Run Lifecycle Means

For PICurv, a "run" is not just a solver launch.
It is the full set of generated artifacts under `runs/<run_id>/`, including:

- normalized runtime config artifacts under `config/`,
- solver and post logs under `logs/`,
- solver outputs and restart files,
- optional scheduler scripts and submission metadata under `scheduler/`.

Key rule:

- every `picurv run --solve ...` creates a fresh run directory,
- `picurv` does not mutate an old run directory in place when you start a new solve,
- restart workflows (`--restart-from`) read from an existing run but still create a new run directory for the restarted continuation,
- continuation workflows (`--continue --run-dir`) resume inside the same run directory.

`run_id` is generated automatically as `<case_basename>_<timestamp>`.

@section p52_newrun_sec 2. Start A New Run

Typical local solve + post:

```bash
./bin/picurv run --solve --post-process -n 4 \
  --case my_case/case.yml \
  --solver my_case/solver.yml \
  --monitor my_case/monitor.yml \
  --post my_case/post.yml
```

Typical cluster solve + post:

```bash
./bin/picurv run --solve --post-process \
  --case my_case/case.yml \
  --solver my_case/solver.yml \
  --monitor my_case/monitor.yml \
  --post my_case/post.yml \
  --cluster my_case/cluster.yml
```

Recommended preflight:

1. `picurv validate ...`
2. `picurv run ... --dry-run`
3. if using Slurm, `picurv run ... --cluster ... --no-submit`

This sequence verifies contract correctness before consuming runtime or queue time.

@section p52_layout_sec 3. Read The Run Directory Correctly

A typical run directory contains:

- `runs/<run_id>/config/`: generated `.control`, BC files, copied YAML inputs, and `post.run`
- `runs/<run_id>/logs/`: solver/postprocessor runtime logs and metrics written by PICurv itself
- `runs/<run_id>/output/`: solver outputs when monitor paths use the default layout
- `runs/<run_id>/scheduler/`: generated Slurm scripts, `submission.json`, and cluster stdout/stderr in cluster mode
- `runs/<run_id>/manifest.json`: top-level run metadata

Practical interpretation:

- if validation succeeds but runtime is wrong, inspect `config/` first,
- if scheduler behavior is wrong, inspect `scheduler/solver.sbatch` or `scheduler/post.sbatch`,
- `scheduler/submission.json` is the source of truth for delayed `submit` and run-directory-based `cancel`,
- if restart/post-only behavior is wrong, confirm the previous run directory contents before changing YAML again.

@section p52_launchers_sec 4. Local, Login-Node, and Batch Launch Resolution

PICurv now separates case physics from site execution policy.

Local multi-rank precedence:

1. `PICURV_MPI_LAUNCHER`
2. `MPI_LAUNCHER`
3. nearest `.picurv-execution.yml`
4. nearest legacy `.picurv-local.yml`
5. built-in `mpiexec`

Cluster batch precedence:

1. `cluster.yml -> execution`
2. nearest `.picurv-execution.yml -> cluster_execution`
3. nearest `.picurv-execution.yml -> default_execution`
4. built-in `srun`

This gives three clean cases:

- workstation users usually need no extra file,
- `picurv init` creates `.picurv-execution.yml` in each new case with inert defaults,
- cluster login-node users can edit `.picurv-execution.yml` when needed,
- batch users can reuse that same file unless `cluster.yml` needs a batch-specific override.

@section p52_restart_sec 5. Restart And Continuation

Restart and continuation use the normal solve workflow.
There is no separate restart command.
The restart source is specified entirely through CLI flags rather than YAML keys.

Three scenarios are supported:

@subsection p52_restart_cfd 5.1 CFD Restart (New Run, Continue Solving)

Use `--restart-from` to create a new run that continues solving from another run's checkpoint data.

```bash
picurv run --solve --restart-from runs/old_run \
  --case case.yml --solver solver.yml --monitor monitor.yml
```

Relevant YAML settings:

- `case.yml`: set `start_step` to the checkpoint step (e.g. 500) and `total_steps` to the desired additional count.
- `solver.yml`: set `eulerian_field_source: "solve"` so the solver advances the Eulerian fields from the restart state.

Operational meaning:

- PICurv copies the checkpoint from `runs/old_run` into the new run's `restart/` directory,
- the solver loads that checkpoint and continues from step 501,
- all new output is written into a fresh `runs/<new_run_id>/` directory.

@subsection p52_restart_particle 5.2 Particle-Tracking Restart (New Run, Pre-Computed Flow)

Use `--restart-from` when the Eulerian flow is already computed and you only need to track particles through it.

```bash
picurv run --solve --restart-from runs/old_flow_run \
  --case case.yml --solver solver.yml --monitor monitor.yml
```

Relevant YAML settings:

- `solver.yml`: set `eulerian_field_source: "load"` so the solver reads pre-computed Eulerian fields instead of advancing them.

Operational meaning:

- PICurv points `restart_dir` directly at the source run's output (no file copy),
- the solver loads the pre-computed flow fields from that source,
- particle tracking proceeds using the loaded fields,
- all new output is written into a fresh `runs/<new_run_id>/` directory.

@subsection p52_continue 5.3 Continue In-Place (Same Run Directory)

Use `--continue --run-dir` to resume a run that was interrupted or stopped early, writing into the same run directory.

```bash
picurv run --solve --continue --run-dir runs/my_run \
  --case case.yml --solver solver.yml --monitor monitor.yml
```

Operational meaning:

- PICurv copies the latest checkpoint from `output/` to `restart/` inside `runs/my_run/`,
- the solver resumes from that checkpoint,
- logs append to the existing log files,
- all output stays within the same `runs/my_run/` directory.

@subsection p52_restart_notes 5.4 Notes

When `start_step > 0`, the initial Eulerian state is always loaded from the restart source regardless of the `eulerian_field_source` setting.
The `eulerian_field_source` value only controls what happens on subsequent steps: `"solve"` advances the fields, `"load"` reads pre-computed fields.

Before launching any restart or continuation:

- verify restart source files exist for the requested step,
- use `start_step` equal to the saved checkpoint step, not the next desired step,
- verify `monitor.yml` directory names match the source run layout.

@section p52_post_sec 6. Postprocess An Existing Run

When solver outputs already exist, reuse the run directory directly:

```bash
./bin/picurv run --post-process \
  --run-dir runs/flat_channel_20260303-120000 \
  --post my_case/alt_analysis.yml
```

Use this when:

- you want a different visualization/statistics recipe,
- solver data are already on disk,
- you do not want to rerun the solver.

PICurv will auto-identify the required case/monitor/control artifacts from `runs/<run_id>/config/`.

Operational patterns for post-only reuse:

- Keep `post.yml` as the full analysis window you want, then use `--continue` to skip steps that were already completed for the same recipe. You do not need to keep editing `start_step` during batch catch-up.

```bash
./bin/picurv run --post-process --continue \
  --run-dir runs/search_robustness_20260322-073415 \
  --post search_robustness_analysis.yml
```

- Live solver example: if `post.yml` requests `0..1000` every `10`, but solver source files currently exist only through step `420`, PICurv launches only `0..420` on the first pass. A later `--continue` run resumes at `430` after those source files appear.
- Interrupted batch example: if `Field_00070.vts` exists but the required MSD CSV still stops at `60`, step `70` is treated as incomplete and the next `--continue` run restarts from `70`.
- Explicit rerun example: if you omit `--continue`, PICurv honors the requested window exactly, rewrites any overlapping VTK files for those steps, and rewrites repeated statistics rows so each step still appears once in the final CSV.
- Changed recipe example: if you point the same `run_dir` at a different `post.yml` recipe, such as adding `Qcrit` or changing the statistics prefix, PICurv starts from that recipe's configured `start_step` instead of inheriting completion from the previous recipe.
- Concurrency rule: PICurv holds a run-directory post lock while the post stage is active. A second post job targeting the same `runs/<run_id>` is refused immediately so two writers cannot race on the same output tree.

@section p52_cluster_sec 7. Batch Job Generation And Reuse

In cluster mode, `picurv` writes scheduler artifacts into the new run directory:

- `scheduler/solver.sbatch`
- `scheduler/solver_<jobid>.out` / `scheduler/solver_<jobid>.err`
- `scheduler/post.sbatch`
- `scheduler/post_<jobid>.out` / `scheduler/post_<jobid>.err`
- `scheduler/submission.json`

Recommended operational pattern:

1. `--dry-run` to confirm launch commands and artifact paths
2. `--no-submit` to inspect generated batch scripts
3. `picurv submit --run-dir runs/<run_id>` only after the scripts look correct
4. `picurv cancel --run-dir runs/<run_id>` when you need to stop a submitted stage without separate job-id bookkeeping

This is especially useful when changing:

- MPI launcher tokens,
- resource counts,
- queue/account settings,
- restart or post-only job behavior.

Operational examples:

```bash
./bin/picurv submit --run-dir runs/<run_id>
./bin/picurv cancel --run-dir runs/<run_id> --stage solve
./bin/picurv cancel --run-dir runs/<run_id> --stage solve --graceful
```

Generated Slurm solver jobs also export runtime walltime metadata into `solver.sbatch`, so the
solver can estimate completed-step cost and request a graceful final write before remaining
walltime gets too tight. If the cluster profile also requests an early signal, PICurv traps
`SIGUSR1`, `SIGTERM`, and `SIGINT`, then uses the same safe-checkpoint final-write path. Use
`signal: "USR1@300"` for `srun`, or `signal: "B:USR1@300"` plus `exec mpirun ...` for direct
`mpirun` batch launches.

For manual cancellation, plain `picurv cancel` is a hard Slurm cancel. Add `--graceful`
when you want the solver to receive `SIGUSR1`, stop at the next safe checkpoint, and write
the latest safe off-cadence output first. Fall back to plain cancel if the job is wedged or
not reaching checkpoints.

@section p52_rules_sec 8. Safe Rules Of Thumb

- Treat `runs/<run_id>/config/` as the ground truth for what the binaries actually consumed.
- Do not hand-edit generated scheduler scripts unless debugging a one-off issue; prefer fixing YAML or `.picurv-execution.yml`.
- Use a fresh restarted run instead of overwriting the previous run directory.
- Use post-only reruns when analysis changes but solver data do not.
- Keep site launcher policy in `.picurv-execution.yml`; keep scheduler policy in `cluster.yml`.
- Keep the shutdown warning window longer than the slowest expected timestep if you rely on the fallback signal path.
- If the runtime walltime guard is too eager or too late for a workload, tune `execution.walltime_guard` in the cluster profile rather than editing generated scripts.

@section p52_refs_sec 9. Related Pages

- **@subpage 06_Simulation_Anatomy**
- **@subpage 36_Cluster_Run_Guide**
- **@subpage 39_Common_Fatal_Errors**
- **@subpage 45_Particle_Initialization_and_Restart**
- **@subpage 49_Workflow_Recipes_and_Config_Cookbook**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Run Lifecycle Guide** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

Treat this page as both a conceptual reference and a runbook. If you are debugging, pair the method/procedure described here with monitor output, generated runtime artifacts under `runs/<run_id>/config`, and the associated solver/post logs so numerical intent and implementation behavior stay aligned.

### What To Extract Before Changing A Case

- Identify which YAML role or runtime stage this page governs.
- List the primary control knobs (paths, stage flags, scheduler settings, restart source, or launcher policy).
- Record expected success indicators (artifact presence, step continuity, job script contents, or stable reused outputs).
- Record failure signals that require rollback or parameter isolation.

### Practical CFD Troubleshooting Pattern

1. Reproduce the workflow on a small case or short time window.
2. Confirm generated `config/` and `scheduler/` artifacts before blaming the solver.
3. Change one lifecycle variable at a time: launcher, resources, restart source, or post recipe.
4. If behavior remains inconsistent, compare against a known-good prior run directory and re-check the generated artifacts.
