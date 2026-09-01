@page 06_Simulation_Anatomy Anatomy of a Simulation

@anchor _Simulation_Anatomy

PICurv runs are composed from modular inputs. This lets you swap physics, numerics, monitoring, and post-analysis independently.

@tableofcontents

@dotfile config_pipeline.dot "From authored YAML to generated artifacts to runtime"

Nothing in the C runtime reads your YAML. `picurv validate` is the only translator:
it checks the contracts, applies non-dimensionalization, and writes the generated
control artifacts the runtime actually consumes. That is why a change that does not
pass validation never reaches the solver, and why the generated files under
`<run.config>` are the authoritative record of what a run was actually told to do.


@section p06_roles_sec 1. Logical Inputs

A single-run workflow uses five logical inputs:

1. `case.yml`: physics/domain/grid/BC/run control
2. `solver.yml`: numerical strategy and solver controls
3. `monitor.yml`: logging, profiling, diagnostics, output cadence, directories
4. `post.yml`: post-processing pipeline recipe
5. Launch settings: stage selection + MPI process count (`-n`)

Cluster/sweep extensions add:

6. `cluster.yml`: Slurm resource/scheduler contract (`run --cluster`, `sweep --cluster`)
7. `study.yml`: parameter matrix + metrics/plot contract (`sweep --study`)

You can choose any file names. C binaries do not require fixed YAML names.

@section p06_compose_sec 2. Composition in Practice

Example solve invocation:

```bash
./bin/picurv run --solve -n 16 \
  --case my_cases/channel_case.yml \
  --solver config/solvers/Imp-MG-Standard.yml \
  --monitor config/monitors/Standard_Output.yml
```

Only numerics/monitoring can be swapped without touching case physics.
In practice, users often keep one validated `case.yml`, then compare multiple `solver.yml` or
`monitor.yml` profiles against it, or reuse one `post.yml` across several cases.

@section p06_artifacts_sec 3. YAML -> Runtime Artifacts

`picurv` validates inputs, then generates C-facing artifacts under `<run.config>/`:

- `<run_id>.control`
- `bcs.run` or `bcs_block*.run`
- `whitelist.run`
- `profile.run`
- `post.run`

These are the concrete contract consumed by C-side parsers in `setup.c`, `io.c`, and BC/profile loaders.

@section p06_libraries_sec 4. Reusable Config Libraries

Project-level shared profiles live in:

```text
config/
|- solvers/
|- monitors/
|- postprocessors/
|- schedulers/
`- studies/
```

Use these as reusable baselines and keep case-specific physics in study-local files.

@section p06_next_steps_sec 5. Next Steps

Proceed to **@subpage 07_Case_Reference**.

For full contract and ingestion mapping, also see:
- **@subpage 14_Config_Contract**
- **@subpage 15_Config_Ingestion_Map**
- **@subpage 16_Config_Extension_Playbook**
- **@subpage 52_Run_Lifecycle_Guide**
- **@subpage 49_Workflow_Recipes_and_Config_Cookbook**
- **@subpage 21_Methods_Overview**
