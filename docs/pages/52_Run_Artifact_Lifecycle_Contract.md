@page 52_Run_Artifact_Lifecycle_Contract Run Artifact Lifecycle Contract

@anchor _Run_Artifact_Lifecycle_Contract

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


@section p52_scope_sec 1. Workspace, Run, And Study Identities

`picurv init <template> --dest <workspace>` creates one self-contained campaign
workspace. Editable files, imported inputs, reusable generated assets, standalone
runs, and studies all have fixed homes beneath it:

```text
<workspace>/
  .picurv-workspace.yml
  config/
    case.yml  solver.yml  monitor.yml  post.yml  cluster.yml
    studies/  grids/  initial_conditions/  inlet_profiles/
    variants/
  inputs/
    grids/  initial_conditions/  inlet_profiles/  reference_fields/
  assets/
    objects/{grids,initial_conditions,inlet_profiles}/
    sets/
    catalog.yml
  runs/
  studies/
```

The directories are created at initialization even when empty. Files inside them
are created only when used. This gives every workspace the same navigable shape
without generating unused data.

`case.yml -> title` is the human run label. A fresh run is named
`<sanitized-title>_<timestamp>`; the title, manifest, and asset hashes carry
identity, so users do not need to rename run directories to remember them. The
timestamp has one-second resolution, so two runs of the same case started within the
same second would name the same directory; the second one takes a `-2` suffix rather
than being refused, because creating a run never writes into an existing one. A study
uses its own title the same way and keeps numbered members below `cases/case_####`.
When an umbrella example contains repeated filenames, initialization preserves the
extra configurations below `config/variants/<original-path>` instead of overwriting
them.

@subsection p52_scope_identity_sub 1.1 Identity Comes From The Manifest

A run publishes what it is in its own `<run.manifest>`: `artifact_type`, `run_id`,
`study_id` and `case_id` for a study member, the `paths` map naming every component
directory, and the `lineage` record described in §5.5. A study publishes the same way
in `study_manifest.json`.

Everything downstream reads those records rather than re-deriving them. `picurv
storage` takes a run's identity and its component boundaries from the manifest, so a
run directory that was renamed, copied, or restored under a different name still
reports the identity it was created with, and still classifies its own output
correctly. A committed checkpoint's step likewise comes from the `-checkpoint_step`
its `checkpoint.meta` records, not from the `step_############` directory it sits in,
so a bundle copied without being re-stamped is not mistaken for the step its new name
claims.

A directory carrying no manifest — one staged by an older release — still works: the
fixed topology above is assumed, and the fallback is reported as `identity_source:
directory-name` rather than presented as a recorded answer.

@subsection p52_scope_layout_sub 1.2 The Run Root Is Closed

A run owns exactly five directories directly below `<run.root>`: `<run.config>`,
`<run.inputs>`, `<run.solver_output>`, `<run.runtime_logs>`, and `<run.scheduler>`.
Staging a run, resuming one, and post-processing one all refuse a run root that has
grown a sixth.

The refusal is not tidiness. A directory beside `<run.solver_output>` is routed by
nothing, so the solver never writes there deliberately; it is classified `unclassified`
by storage, so it is archived in every copy of the run and pruned by no policy; and it
is absent from the `paths` map in `<run.manifest>`, so no reader can say what it holds.
Catching it when the run is staged is the last point at which the answer is still
"move it".

Scientific output belongs under `<run.checkpoints>`, `<run.analysis>`, or
`<run.visualization>`. An unexpected *file* at the run root is reported and allowed:
a stray note costs nothing, and refusing to resume a long campaign over one would be a
worse failure than the one being prevented.

@section p52_newrun_sec 2. Editable Configurations And Imported Files

Run commands from the workspace and refer to canonical configuration paths:

```bash
picurv validate \
  --case config/case.yml --solver config/solver.yml --monitor config/monitor.yml

picurv run --solve --post-process -n 4 \
  --case config/case.yml --solver config/solver.yml \
  --monitor config/monitor.yml --post config/post.yml
```

Generator descriptions such as grid `.cfg` files and `expressions.cfg` belong in
`<workspace>/config/grids/`, `<workspace>/config/initial_conditions/`, or
`<workspace>/config/inlet_profiles/`. User-supplied data belong in the matching
`<workspace>/inputs/` directory. Register an external
file explicitly so its checksum and ownership mode are catalogued:

```bash
picurv inputs import grid /archive/mesh.picgrid --mode copy
picurv inputs import initial-condition /archive/ufield.dat --mode reflink
picurv inputs import reference-field /shared/baseline.dat --mode reference
```

`reference` writes a small `.reference.yml`; it never silently copies, archives, or
prunes the external target. Configuration paths resolve from the workspace root.
Missing configured inputs fail loudly before submission.

Recommended preflight is `validate`, then `run --dry-run`, then `run --no-submit`.
The dry run reports the proposed run identity, canonical paths, asset actions, and
blocking provider dependencies without creating the directory.

@section p52_layout_sec 3. Reusable Assets And The Run Snapshot

`picurv precompute --case <workspace>/config/case.yml` resolves a dependency graph for grid,
initial-condition, and inlet-profile providers. File and Python providers build in
an isolated temporary run layout. PICurv publishes their complete output only after
every selected provider succeeds:

```text
assets/objects/<kind>/<content-hash>/
  asset.json          identity, provider settings, checksums, provenance
  validation.json     what was produced: dimensions, bounds, fields, checks
  preview.vts         inspectable geometry, when the grid is small enough
  summary.json        initial-condition statistics, when the provider emits them
  spectrum.csv        initial-condition spectrum, when the provider emits one
  payload/...         the solver-ready files, at their run-relative paths
assets/sets/<case-name>-<path-hash>.yml
```

Precompute exists so a grid, field, or profile can be looked at and corrected before a
solve is committed, so every published object carries inspection material beside its
payload. A preview is skipped, and `validation.json` says so, when the first block
exceeds two million nodes: past that the ASCII geometry costs more than the look is
worth. Inspection material describes the payload and is deliberately excluded from the
object's identity, so changing a preview format does not re-identify every asset.

Generated destinations are PICurv's. `grid.generator.output_file`, `stats_file`, and
`vts_file` are rejected rather than honoured: a configuration file that names its own
output path creates a competing directory outside the store, which is what the fixed
topology exists to prevent.

The payload keeps its run-relative shape rather than a flat set of canonical filenames,
because an inlet-profile asset holds one file per block and face and a flat root cannot
express that. Materialization exposes each payload path into the run at the same
relative location.

The object identity includes normalized provider settings, the case values each build
reads, the identities of the assets it depends on, checksums of referenced
files, and the PICurv build. Changing an equation, grid config, imported field, inlet
parameters, or generator code therefore selects a new object. Unchanged inputs reuse
the existing object. `--only grid,initial-condition` selects a dependency closure;
publication remains all-or-nothing for that invocation.

A provider executed only in C is reported as `runtime-c`. Precompute does not imitate
it in Python. If a Python initial-condition generator needs a grid configured as
`programmatic_c`, precompute and run planning fail and tell the user to choose `file`
or `grid_gen`. Otherwise the simulator announces the runtime provider and generates it
during startup.

Fresh runs reuse valid assets automatically. `--require-precomputed` refuses a missing
or stale object instead of building it, while `--fetch-missing` searches configured
storage archives before rebuilding. Each run receives physical run-local files by
reflink, hardlink, or copy and records the exact object/checksum mapping in
`inputs/assets.lock.yml`.

Every run has this fixed shape:

```text
runs/<title>_<timestamp>/
  manifest.json
  config/
    active.json  case.yml  solver.yml  monitor.yml  cluster.yml
    history/<revision>/
    post-recipes/<recipe-id>/{post.yml,post.run,state.json}
  inputs/{grid,initial_condition,inlet_profiles,restart}/
  output/
    checkpoints/
    analysis/{metrics,statistics,spectra,plots}/
    visualization/
  logs/
  scheduler/
```

No peer `diagnostics/`, `results/`, or arbitrary monitor-selected output root is
created. C runtime metrics use `<run.analysis.metrics>`; post statistics and spectra
use their analysis homes; renderable VTK output uses `<run.visualization>`.

The initial YAML snapshot is immutable evidence. An in-place continuation stores the
new YAML and generated controls under `<run.config.history>/<revision>/` and updates
`<run.config.active>`; it does not erase the original. `manifest.json` records the
workspace identity, active build, canonical paths, stages, locked assets, per-executable
build identity (`binaries`), and each stage's lifecycle state (`components`: e.g.
not_requested, planned, complete, offloaded).

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
  into the new run's `inputs/restart/checkpoints/` directory,
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
- PICurv validates the requested immutable bundle in `<run.checkpoints>` and
  materializes it into `<run.inputs>/restart`, the one restart home every mode uses,
  so the generated control carries a single canonical path. Same-filesystem reflink
  or hardlink materialization makes this metadata-cheap rather than a second copy,
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

@subsection p52_restart_lineage 5.5 A Branch Records What It Branched From

A run created with `--restart-from` writes a `lineage` record into its
`<run.manifest>`:

```json
"lineage": {
  "relationship": "branch",
  "parent_run_id": "channel_20260101-120000",
  "parent_study_id": null,
  "parent_case_id": null,
  "parent_identity_source": "manifest",
  "parent_path": "runs/channel_20260101-120000",
  "checkpoint_step": 4000,
  "statistics_state": "carry",
  "requested_source": "latest",
  "recorded_at": "2026-01-02T09:15:00-06:00"
}
```

Without it a branch is indistinguishable from a fresh run. The bundle copied into the
run's restart input carries the parent's geometry digest and software identity, but names
no run, so the trajectory a result belongs to could not be reconstructed afterwards.

The parent is named from its own manifest, so a parent that was later renamed is still
identified correctly; `parent_identity_source` says whether that lookup succeeded or
fell back to the directory name. `requested_source` records what was asked for, which
is how a `--restart-from latest` is distinguished from the explicit path it resolved
to.

A fresh run records `{"relationship": "root"}` rather than omitting the key, so a
reader never has to distinguish "started from nothing" from "written before lineage
existed". An in-place `--continue` keeps whatever lineage the run was created with:
what a run branched from does not change because it was resumed.

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

PICurv auto-identifies the active case, monitor, and control artifacts from
`<run.config.active>` and its referenced revision. Every normalized recipe gets a stable
ID. Its controls live under `<run.post_recipes>/<recipe-id>/`, its field visualization
under `<run.visualization>/<recipe-id>/`, and its statistics/spectra below the matching
analysis directories. Two changed recipes therefore coexist instead of overwriting one
another.

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
- Concurrency rule: PICurv holds a post lock while the stage is active. A second writer
  targeting the same output lineage is refused so generated controls and result files
  cannot race.

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

@section p52_cap_restart_stats_sec 8.6 Restart Statistics State Entries

@htmlinclude generated/capability_inventory_run_restart_statistics_state.html

@subsection p52_cap_restart_stats_reset_sub reset

@anchor p52_cap_restart_stats_reset

**Identity.** `picurv run --solve --restart-from ... --statistics-state reset`.

**What it does.** Starts field-statistics accumulators empty in the new branched run;
the checkpoint still seeds physical fields.

**When to choose it.** Use `reset` for a changed case, a new measurement window, or
whenever the new run should stand alone statistically. It is safer than `carry` when
solver or monitoring definitions changed.

Neither value is a default. Branching a run with `field_statistics.enabled: true` and
no `--statistics-state` is refused, because guessing is wrong in both directions:
resetting silently reports a shorter average than the user expects, and carrying
averages two trajectories together. With field statistics disabled there is nothing to
decide and the flag stays optional.

**Parameters it owns.** None.

**Interactions.** Applies to a new run created with `--restart-from`; in-place
`--continue` already continues the same run state. It does not change checkpoint
selection or Eulerian restart authority.

**Diagnostics.** The generated control omits `-field_statistics_continue true`, and
the new run's active configuration records `statistics_state: reset`.

**Evidence.** Unit verified — `tests/test_workspace_lifecycle.py` exercises the CLI
choice and generated run control surface.

**Limitations.** It cannot merge old and new statistics later; preserve or postprocess
the old window separately if both are needed.

@subsection p52_cap_restart_stats_carry_sub carry

@anchor p52_cap_restart_stats_carry

**Identity.** `picurv run --solve --restart-from ... --statistics-state carry`.

**What it does.** Requests compatible checkpointed field-statistics accumulator state
for the new run before sampling continues.

**When to choose it.** Use `carry` only when the grid, field-statistics windows,
weighting, fields, and covariance definitions are unchanged and the new run is a
scientific continuation of the same measurement.

**Parameters it owns.** None.

**Interactions.** Emits `-field_statistics_continue true`. Checkpoint compatibility
validation still owns grid layout and required state; changing the statistics recipe
can make the requested state unusable.

**Diagnostics.** The generated control contains the continue flag. A missing or
incompatible checkpoint statistics payload fails during restart setup rather than
silently resetting.

**Evidence.** Unit verified — `tests/test_workspace_lifecycle.py` exercises the CLI
choice and generated run control surface.

**Limitations.** Experimental: numerical equivalence across every change of MPI layout
and optional field combination has not been characterized.

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
