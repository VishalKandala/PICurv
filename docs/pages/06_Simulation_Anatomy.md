@page 06_Simulation_Anatomy Anatomy of a Simulation

PICurv runs are composed from modular inputs. This lets you swap physics, numerics, monitoring, and post-analysis independently.

@tableofcontents

@section roles_sec 1. Logical Inputs

A single-run workflow uses five logical inputs:

1. `case.yml`: physics/domain/grid/BC/run control
2. `solver.yml`: numerical strategy and solver controls
3. `monitor.yml`: logging, profiling, output cadence, directories
4. `post.yml`: post-processing pipeline recipe
5. Launch settings: stage selection + MPI process count (`-n`)

Cluster/sweep extensions add:

6. `cluster.yml`: Slurm resource/scheduler contract (`run --cluster`, `sweep --cluster`)
7. `study.yml`: parameter matrix + metrics/plot contract (`sweep --study`)

You can choose any file names. C binaries do not require fixed YAML names.

@section compose_sec 2. Composition in Practice

Example solve invocation:

```bash
./scripts/pic.flow run --solve -n 16 \
  --case my_cases/channel_case.yml \
  --solver config/solvers/Imp-MG-Standard.yml \
  --monitor config/monitors/Standard_Output.yml
```

Only numerics/monitoring can be swapped without touching case physics.

@section artifacts_sec 3. YAML -> Runtime Artifacts

`pic.flow` validates inputs, then generates C-facing artifacts under `runs/<run_id>/config/`:

- `<run_id>.control`
- `bcs.run` or `bcs_block*.run`
- `whitelist.run`
- `profile.run`
- `post.run`

These are the concrete contract consumed by C-side parsers in `setup.c`, `io.c`, and BC/profile loaders.

@section libraries_sec 4. Reusable Config Libraries

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

@section next_steps_sec 5. Next Steps

Proceed to **@subpage 07_Case_Reference**.

For full contract and ingestion mapping, also see:
- **@subpage 14_Config_Contract**
- **@subpage 15_Config_Ingestion_Map**
- **@subpage 16_Config_Extension_Playbook**
- **@subpage 21_Methods_Overview**
