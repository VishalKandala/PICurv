@page 52_Run_Lifecycle_Guide Run Lifecycle Guide

@anchor _Run_Lifecycle_Guide

@pagemeta{How-to, Operators running and restarting jobs, Current behavior}

This page explains how a PICurv run moves from a new solve to restart, post-processing reuse, and cluster job generation.
It is the operational view of the run directory lifecycle.

For verified remote protection, cold offload, and later restoration of these artifacts, see **@subpage 61_Storage_Management_Guide**.

@tableofcontents

@dotfile run_lifecycle.dot "States a run moves through, and the commands that move it"

Every arrow is a command you can issue. `--continue` resumes the same run in place;
`--restart-from` seeds a **new** run from an existing checkpoint. Those are different
operations with different provenance, and choosing between them is the decision this
page exists to support.


@section p52_scope_sec 1. What A Run Lifecycle Means

For PICurv, a "run" is not just a solver launch.
It is the full set of generated artifacts under `runs/<run_id>/`, including:

- normalized runtime config artifacts under `<run.config>/`,
- solver and post logs under `<run.runtime_logs>/`,
- solver outputs and restart files,
- optional scheduler scripts and submission metadata under `<run.scheduler>/`.

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
3. `picurv run ... --no-submit` to stage local commands, or add `--cluster ...` to stage Slurm scripts

This sequence verifies contract correctness before consuming runtime or queue time.

@section p52_layout_sec 3. Read The Run Directory Correctly

A typical run directory contains:

- `<run.config>/`: generated `.control`, BC files, copied YAML inputs,
  and `post.run`
- `<run.runtime_logs>/`: solver/postprocessor runtime logs and metrics written by PICurv itself
- `<run.solver_output>/checkpoints/step_<12 digits>/`: immutable committed
  solver checkpoints when monitor paths use the default layout
- `<run.scheduler>/`: generated Slurm scripts, `submission.json`, and cluster stdout/stderr in cluster mode
- `runs/<run_id>/manifest.json`: top-level run metadata

Practical interpretation:

- if validation succeeds but runtime is wrong, inspect `<run.config>/` first,
- if scheduler behavior is wrong, inspect `<run.scheduler>/solver.sbatch` or `<run.scheduler>/post.sbatch`,
- `<run.scheduler>/submission.json` is the source of truth for delayed `submit` and run-directory-based `cancel`,
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

- PICurv validates the requested committed bundle and copies that whole bundle
  into the new run's `restart/checkpoints/` directory,
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

- `case.yml -> run_control.start_step` must be the saved checkpoint step and must be greater than zero; use a normal run without `--continue` for a fresh start at step zero,
- PICurv validates and reads the requested immutable bundle directly from
  `<run.solver_output>/checkpoints/`; it does not create a second in-place copy,
- the solver resumes from that checkpoint,
- logs append to the existing log files and remain independent of checkpoint payloads,
- all output stays within the same `runs/my_run/` directory.

@subsection p52_restart_notes 5.4 Notes

When `start_step > 0`, the initial Eulerian state is always loaded from the restart source regardless of the `eulerian_field_source` setting.
The `eulerian_field_source` value only controls what happens on subsequent steps: `"solve"` advances the fields, `"load"` reads pre-computed fields.

Before launching any restart or continuation:

- verify a committed bundle exists for the requested step,
- use `start_step` equal to the saved checkpoint step, not the next desired step,
- do not edit bundle contents or generated control files by hand.

`--continue` is for the same physical case. It permits changes to run-control
time settings, solver parameters/type, monitoring, postprocessing, and resource
settings. Other `case.yml` changes are rejected; use `--restart-from` to create
a new branch. `--restart-from` may change case physics, but C still requires the
same grid geometry/layout. A newly enabled optional subsystem such as particles
uses its normal initialization path unless its restart mode explicitly requests
saved state; subsequent checkpoints then include that subsystem.

@section p52_post_sec 6. Postprocess An Existing Run

When solver outputs already exist, reuse the run directory directly:

```bash
./bin/picurv run --post-process \
  --run-dir runs/flat_channel_20260303-120000 \
  --post my_case/alt_analysis.yml
```

Use this when:

- you want a different <run.visualization>/statistics recipe,
- solver data are already on disk,
- you do not want to rerun the solver.

PICurv will auto-identify the required case/monitor/control artifacts from `<run.config>/`.

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

- `<run.scheduler>/solver.sbatch`
- `<run.scheduler>/solver_<jobid>.out` / `<run.scheduler>/solver_<jobid>.err`
- `<run.scheduler>/post.sbatch`
- `<run.scheduler>/post_<jobid>.out` / `<run.scheduler>/post_<jobid>.err`
- `<run.scheduler>/submission.json`

Recommended operational pattern:

1. `--dry-run` to confirm launch commands and artifact paths
2. `--no-submit` to inspect generated local commands or batch scripts
3. `picurv submit --run-dir runs/<run_id>` only after the staged artifacts look correct
4. `picurv cancel --run-dir runs/<run_id>` for Slurm runs when you need to stop a submitted stage without separate job-id bookkeeping

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
when you want the solver process tree to receive `SIGUSR1`, stop at the next safe checkpoint, and write
the latest safe off-cadence output first. Fall back to plain cancel if the job is wedged or
not reaching checkpoints.

@section p52_compat_sec 8. What You Can Change, And How To Continue

Before changing anything about a run you intend to carry forward, decide which of
three operations you want:

- **Continue in place** (`--continue`) — same run, same identity, more timesteps.
- **Restart into a new run** (`--restart-from`) — a new run seeded from an existing
  checkpoint, with its own identity and provenance.
- **Fresh run** — no inherited state.

@subsection p52_compat_axes_sub 8.1 Three Different Questions

The tables below keep three things apart that are easy to conflate:

- **Accepted?** — whether the tooling permits it. This is a fact about the software.
- **Advised?** — whether it is scientifically sensible. The tool permitting something
  does not make it a good idea.
- **Basis** — how the row is known: `enforced` (validation accepts or rejects it),
  `tested` (a test exercises it), `characterized` (measured, not gated), `unknown`
  (not established).

@warning `unknown` means nobody has checked. It is **not** a synonym for "no". Where a
cell says unknown, the tool may well permit the change and produce something
meaningless without complaining.

@subsection p52_compat_inplace_sub 8.2 Continue In Place (`--continue`)

Continuation guards **physical case identity**: `validate_continue_case_identity`
compares a hash of `case.yml` with `run_control` and the particle `restart_mode`
excluded. Everything else in `case.yml` is part of that identity.

| Change | Accepted? | Advised? | Basis |
|---|---|---|---|
| `run_control` (more timesteps) | yes | yes | enforced — excluded from the identity hash |
| Particle `restart_mode` | yes | yes | enforced — excluded from the identity hash |
| Monitor cadence / logging | yes | yes | enforced — `monitor.yml` is outside the identity |
| Post-processing recipe | yes | yes | enforced — `post.yml` is outside the identity |
| Solver tolerances | yes | discouraged | enforced (accepted); effect not characterized |
| Momentum solver selection | yes | discouraged | enforced (accepted); effect unknown |
| Timestep size | yes | discouraged | enforced (accepted); effect not characterized |
| Boundary conditions | **no** | — | enforced — part of the case identity |
| Turbulence model | **no** | — | enforced — part of the case identity |
| Physical properties | **no** | — | enforced — part of the case identity |
| Grid dimensions or geometry | **no** | — | enforced — part of the case identity |
| MPI rank count | unknown | avoid | unknown — no test covers rank change under `--continue` |

@note Anything reachable through `solver.yml`, `monitor.yml`, or `post.yml` is
**accepted** during continuation by design; the guard is on the physical case only.
The "discouraged" rows are a scientific judgement, not a tooling restriction: the
resulting trajectory concatenates two differently-converged segments, which is rarely
what you want in a published result.

@subsection p52_compat_restart_sub 8.3 Restart Into A New Run (`--restart-from`)

A new run carries its own identity, so the case-identity guard does not apply.

| Change | Accepted? | Advised? | Basis |
|---|---|---|---|
| Anything in solver / monitor / post | yes | yes | enforced |
| Boundary conditions | yes | yes, with care | enforced (accepted); physical continuity unknown |
| Turbulence model | yes | yes, with care | enforced (accepted); effect unknown |
| MPI rank count | yes | yes | tested — the field-statistics suite covers a changed rank count across restart |
| Physical properties | yes | treat as a new study | enforced (accepted); the seeded field is no longer consistent with the new physics |
| Grid dimensions or geometry | **no** | — | enforced — the checkpoint cannot be interpreted on a different grid |

@subsection p52_compat_exact_sub 8.4 Restart Is Not Bit-Exact

Even with no changes, a restarted run does not reproduce an uninterrupted one bit for
bit. The offset is a structural floor from boundary re-application at the restart
step, it decays rather than amplifies, and it is **invariant to solver tolerance** —
tightening tolerances does not remove it. Measurements are in
**@subpage 29_Maintenance_Backlog**.

Acceptance criteria for restart equivalence must therefore be set from that measured
floor, not from machine epsilon.

@subsection p52_compat_flux_sub 8.5 Driven-Flux Targets Across Segments

@ref p44_cap_constant_flux "constant_flux" re-reads its target from the bcs file on
every run, so editing it between segments takes effect deliberately.
@ref p44_cap_initial_flux "initial_flux" instead restores its latched target from the
checkpoint, so a resumed run holds the original target rather than re-measuring a
drifted one.

@section p52_rules_sec 9. Safe Rules Of Thumb

- Treat `<run.config>/` as the ground truth for what the binaries actually consumed.
- Do not hand-edit generated scheduler scripts unless debugging a one-off issue; prefer fixing YAML or `.picurv-execution.yml`.
- Use a fresh restarted run instead of overwriting the previous run directory.
- Use post-only reruns when analysis changes but solver data do not.
- Keep site launcher policy in `.picurv-execution.yml`; keep scheduler policy in `cluster.yml`.
- Keep the shutdown warning window longer than the slowest expected timestep if you rely on the fallback signal path.
- If the runtime walltime guard is too eager or too late for a workload, tune `execution.walltime_guard` in the cluster profile rather than editing generated scripts.

@section p52_refs_sec 10. Related Pages

- **@subpage 06_Simulation_Anatomy**
- **@subpage 36_Cluster_Run_Guide**
- **@subpage 39_Common_Fatal_Errors**
- **@subpage 45_Particle_Initialization_and_Restart**
- **@subpage 49_Workflow_Recipes_and_Config_Cookbook**
- **@subpage 58_Field_Statistics** (statistics state inside the committed checkpoint bundle)
