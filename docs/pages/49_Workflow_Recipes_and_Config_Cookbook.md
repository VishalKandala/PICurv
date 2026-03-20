@page 49_Workflow_Recipes_and_Config_Cookbook Workflow Recipes and Config Cookbook

@anchor _Workflow_Recipes_and_Config_Cookbook

This page is the practical companion to the reference docs.
Use it when you already understand the config roles and want proven patterns for combining them cleanly.

@tableofcontents

@section p49_idea_sec 1. Core Workflow Idea

PICurv is designed around modular config roles:

- `case.yml`: physics, grid, BCs, run duration
- `solver.yml`: numerical strategy
- `monitor.yml`: logging, output cadence, directories
- `post.yml`: analysis recipe
- optional `cluster.yml`: scheduler policy
- optional `study.yml`: parameter sweep definition

The intended UX is:

- keep each file focused on one concern,
- reuse stable profiles across many runs,
- only clone and specialize the roles that actually changed.

@section p49_combos_sec 2. Common Mix-and-Match Patterns

Pattern A: same case, different monitoring

- keep one trusted `case.yml`
- pair it with a lightweight `monitor.yml` for fast iteration
- switch to a richer `monitor.yml` for production output

Pattern B: same case, different analysis

- run once with a standard monitor profile
- post-process the same run directory multiple times with different `post.yml` recipes

Pattern C: same solver profile reused across many cases

- keep one validated `solver.yml` for a known momentum strategy
- pair it with multiple physical cases as long as the strategy is appropriate

Pattern D: same cluster profile across many jobs

- reuse one `cluster.yml` for site policy
- vary only case/solver/monitor/post inputs

Pattern E: same case, different local launch environment

- keep case/solver/monitor/post unchanged
- use default `mpiexec` on ordinary workstations
- add a repo-local or case-local `.picurv-execution.yml` when a cluster login node or cluster batch job needs site-specific launcher flags

@section p49_recipes_sec 3. Command Recipes

First-time local workflow:

```bash
./bin/picurv init flat_channel --dest my_case
./bin/picurv validate --case my_case/flat_channel.yml --solver my_case/Imp-MG-Standard.yml --monitor my_case/Standard_Output.yml --post my_case/standard_analysis.yml
./bin/picurv run --solve --post-process -n 4 --case my_case/flat_channel.yml --solver my_case/Imp-MG-Standard.yml --monitor my_case/Standard_Output.yml --post my_case/standard_analysis.yml
```

Solver-only run:

```bash
./bin/picurv run --solve -n 8 --case case.yml --solver solver.yml --monitor monitor.yml
```

Solver-only run on a cluster login node with a persistent site launcher config:

```bash
# Preferred: create one repo-root `.picurv-execution.yml` in your cluster clone.
# `picurv init` seeds new cases from it automatically.
# For an existing case, copy config/runtime/execution.example.yml once and rename it.
./bin/picurv run --solve -n 8 --case case.yml --solver solver.yml --monitor monitor.yml
```

Cluster run that reuses the same site launcher defaults unless `cluster.yml` overrides them:

```bash
# Preferred: create one repo-root `.picurv-execution.yml` in your cluster clone.
# `picurv init` seeds new cases from it automatically.
# For an existing case, copy config/runtime/execution.example.yml once and rename it.
./bin/picurv run --solve --post-process --case case.yml --solver solver.yml --monitor monitor.yml --post post.yml --cluster cluster.yml
```

Post-process an existing run with a different recipe:

```bash
./bin/picurv run --post-process --run-dir runs/<run_id> --post alt_analysis.yml
```

Dry-run planning:

```bash
./bin/picurv run --solve --post-process --case case.yml --solver solver.yml --monitor monitor.yml --post post.yml --dry-run --format json
```

Cluster generation without submit:

```bash
./bin/picurv run --solve --post-process --case case.yml --solver solver.yml --monitor monitor.yml --post post.yml --cluster cluster.yml --no-submit
```

Delayed submit from existing staged artifacts:

```bash
./bin/picurv submit --run-dir runs/<run_id>
```

Cancel a submitted run by directory:

```bash
./bin/picurv cancel --run-dir runs/<run_id> --stage solve
```

Generated Slurm solver jobs already enable the runtime walltime guard by default. Override it only
when needed:

```yaml
execution:
  walltime_guard:
    enabled: true
    warmup_steps: 10
    multiplier: 2.0
    min_seconds: 60
    estimator_alpha: 0.35
```

Keep an early signal as fallback protection for preemption/manual termination:

```yaml
execution:
  extra_sbatch:
    signal: "USR1@300"
```

For direct `mpirun` batch launches, use `signal: "B:USR1@300"` and prefer `exec mpirun ...`.

Run naming note:

- `cluster.yml` does not provide a run-name field.
- `picurv` creates `runs/<case_basename>_<timestamp>/` automatically.
- generated Slurm job names are derived from that run ID.

@section p49_examples_sec 4. Example Config Combinations

Programmatic baseline:

- case: `examples/flat_channel/flat_channel.yml`
- solver: `examples/flat_channel/Imp-MG-Standard.yml`
- monitor: `examples/flat_channel/Standard_Output.yml`
- post: `examples/flat_channel/standard_analysis.yml`

File-grid baseline:

- case: `examples/bent_channel/bent_channel.yml`
- solver: `examples/bent_channel/Imp-MG-Standard.yml`
- monitor: `examples/bent_channel/Standard_Output.yml`
- post: `examples/bent_channel/standard_analysis.yml`

Analytical zero-flow Brownian verification:

- case: `examples/brownian_motion/brownian_motion.yml`
- solver: `examples/brownian_motion/Analytical-Zero.yml`
- monitor: `examples/brownian_motion/Standard_Output.yml`
- post: `examples/brownian_motion/brownian_analysis.yml`

Analytical uniform-flow deterministic drift verification:

- case: `examples/drift_uniform_flow/case.yml`
- solver: `examples/drift_uniform_flow/solver.yml`
- monitor: `examples/drift_uniform_flow/monitor.yml`
- post: `examples/drift_uniform_flow/post.yml`

Analytical zero-flow plus verification-source diffusivity drift verification:

- case: `examples/drift_diffusivity_gradient/case.yml`
- solver: `examples/drift_diffusivity_gradient/solver.yml`
- monitor: `examples/drift_diffusivity_gradient/monitor.yml`
- post: `examples/drift_diffusivity_gradient/post.yml`

These are examples, not fixed bundles.
You can intentionally swap roles when the contract makes sense.

@section p49_drift_verification_sec 5. Drift Verification Workflow

Use the drift examples as separate, end-to-end checks for the two deterministic particle-position terms:

- `examples/drift_uniform_flow/` isolates the carrier-velocity term in `X_new = X_old + (vel + diffusivitygradient) * dt + Brownian`.
- `examples/drift_diffusivity_gradient/` isolates the diffusivity-gradient term while keeping the carrier flow at zero.

The intended policy is:

- prefer ordinary production pathways first,
- use analytical solver modes when the solver already supports them,
- only use `verification.*` overrides when there is no cleaner end-to-end way to prescribe the field needed for verification.

For the current diffusivity-drift case, the curated source override lives behind:

```yaml
verification:
  sources:
    diffusivity:
      mode: "analytical"
      profile: "LINEAR_X"
      gamma0: 1.0e-3
      slope_x: 2.0e-4
```

Contributor note:

- new verification source injections must be implemented in `include/verification_sources.h` and `src/verification_sources.c`
- production call sites should stay thin and delegate into that module
- do not introduce one-off verification flags directly into unrelated solver files unless there is no alternative and the delegation point remains explicit

@section p49_config_patterns_sec 6. Configuration Patterns Worth Reusing

Minimal zero-flow startup:

```yaml
properties:
  initial_conditions:
    mode: "Zero"
```

Poiseuille with a clear scalar input:

```yaml
properties:
  initial_conditions:
    mode: "Poiseuille"
    peak_velocity_physical: 1.25
```

Generated grid workflow:

```yaml
grid:
  mode: grid_gen
  generator:
    config_file: config/grids/coarse_square_tube_curved.cfg
    grid_type: cpipe
    cli_args: ["--ncells-i", "96", "--ncells-j", "96"]
```

Restart with explicit particle re-initialization:

```yaml
run_control:
  start_step: 200

models:
  physics:
    particles:
      count: 10000
      restart_mode: "init"
```

Continue a run that ended at step 500 for 1000 more steps:

```yaml
run_control:
  start_step: 500
  total_steps: 1000
```

```yaml
operation_mode:
  eulerian_field_source: "load"
```

Optional particle choice:

```yaml
models:
  physics:
    particles:
      restart_mode: "load"   # or "init"
```

```bash
./bin/picurv run --solve --post-process \
  --case restart_case/case.yml \
  --solver restart_case/solver.yml \
  --monitor restart_case/monitor.yml \
  --post restart_case/post.yml
```

Expected step range:

- saved state loaded from step `500`
- first new step is `501`
- final step after this run is `1500`

@section p49_choose_sec 7. Which Example To Start From

Choose `flat_channel` when:

- you want the simplest first runnable case,
- you are learning `programmatic_c`,
- you want a clean baseline for solver/monitor/post swaps.

Choose `bent_channel` when:

- you need file-based curvilinear geometry,
- you want to understand `grid.mode: file`,
- you plan to compare against `grid_gen` later.

Choose `brownian_motion` when:

- you want an analytical zero-flow validation case,
- you want a particle-focused test with minimal flow complexity,
- you want a reference for `ZERO_FLOW` analytical mode.

@section p49_next_steps_sec 7. Related Pages

- **@subpage 05_The_Conductor_Script**
- **@subpage 07_Case_Reference**
- **@subpage 14_Config_Contract**
- **@subpage 52_Run_Lifecycle_Guide**
- **@subpage 48_Grid_Generator_Guide**
- **@subpage 45_Particle_Initialization_and_Restart**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Workflow Recipes and Config Cookbook** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
