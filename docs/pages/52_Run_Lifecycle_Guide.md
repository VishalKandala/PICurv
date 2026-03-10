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
- restart workflows read from an existing run but still create a new run directory for the restarted continuation.

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
- `runs/<run_id>/logs/`: solver/post logs
- `runs/<run_id>/results/`: solver outputs when monitor paths use the default layout
- `runs/<run_id>/scheduler/`: generated Slurm scripts and `submission.json` in cluster mode
- `runs/<run_id>/manifest.json`: top-level run metadata

Practical interpretation:

- if validation succeeds but runtime is wrong, inspect `config/` first,
- if scheduler behavior is wrong, inspect `scheduler/solver.sbatch` or `scheduler/post.sbatch`,
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
- cluster login-node users can add `.picurv-execution.yml`,
- batch users can reuse that same file unless `cluster.yml` needs a batch-specific override.

@section p52_restart_sec 5. Restart From An Existing Run

Restart uses the normal solve workflow.
There is no separate restart command.

Example:

```yaml
run_control:
  start_step: 500
  total_steps: 1000
  restart_from_run_dir: "../runs/flat_channel_20260303-120000"
```

```yaml
operation_mode:
  eulerian_field_source: "load"
```

```yaml
models:
  physics:
    particles:
      restart_mode: "load"   # or "init"
```

Operational meaning:

- previous run ended at step `500`,
- new run loads state from that old run,
- first new advanced step is `501`,
- new run finishes at `1500`,
- restarted continuation is written into a new `runs/<new_run_id>/` directory.

Before launching:

- verify restart source files exist for the requested step,
- use `start_step` equal to the saved step, not the next desired step,
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

@section p52_cluster_sec 7. Batch Job Generation And Reuse

In cluster mode, `picurv` writes scheduler artifacts into the new run directory:

- `scheduler/solver.sbatch`
- `scheduler/post.sbatch`
- `scheduler/submission.json`

Recommended operational pattern:

1. `--dry-run` to confirm launch commands and artifact paths
2. `--no-submit` to inspect generated batch scripts
3. full submit only after the scripts look correct

This is especially useful when changing:

- MPI launcher tokens,
- resource counts,
- queue/account settings,
- restart or post-only job behavior.

@section p52_rules_sec 8. Safe Rules Of Thumb

- Treat `runs/<run_id>/config/` as the ground truth for what the binaries actually consumed.
- Do not hand-edit generated scheduler scripts unless debugging a one-off issue; prefer fixing YAML or `.picurv-execution.yml`.
- Use a fresh restarted run instead of overwriting the previous run directory.
- Use post-only reruns when analysis changes but solver data do not.
- Keep site launcher policy in `.picurv-execution.yml`; keep scheduler policy in `cluster.yml`.

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
