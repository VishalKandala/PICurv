@mainpage PICurv Solver Documentation

@section intro_sec Introduction

# PICurv: A Hybrid Eulerian-Lagrangian Framework for Scalar Transport

PICurv is a parallel CFD and particle framework for incompressible flow and scalar transport. It combines:

- Eulerian flow solves on curvilinear grids (PETSc `DMDA` based)
- Lagrangian particle transport (`DMSwarm` based)
- A YAML-driven run pipeline through the built conductor `./bin/picurv`

The user workflow is configuration-first: YAML inputs are validated and translated into generated runtime artifacts consumed by the C solver and postprocessor.
Those YAML roles are modular by design: `case.yml`, `solver.yml`, `monitor.yml`, and `post.yml` are meant to be recombined instead of rewritten from scratch for every run.

@section get_started_main Documentation Entry Points (Structural)

Start here first:

- **@subpage 41_Getting_Started_Index**: fast first-run path for installation, first case setup, validation, and first outputs.

Then use the structural map and reference pages:

- **@subpage Documentation_Map**: complete categorized index by workflow, artifact, and runtime layer.
- **@subpage 06_Simulation_Anatomy**: end-to-end run structure and generated artifacts.
- **@subpage 14_Config_Contract**: YAML contracts and Python-to-C handoff.
- **@subpage 49_Workflow_Recipes_and_Config_Cookbook**: practical profile-composition patterns and runnable recipes.
- **@subpage 21_Methods_Overview**: governing methods and model map.
- **@subpage 46_C_Runtime_Execution_Map**: solver execution order and C runtime touchpoints.
- **@subpage 53_Search_Robustness_Metrics_Reference**: authoritative definitions for runtime search observability and paper-grade search signals.
- **Data Structures (`annotated_structured.html`)** and **File List (`files_structured.html`)**:
  organized code-reference views by module and file type.

@section preview_sec Quick Preview

@htmlonly
<div style="text-align:center; margin:1rem 0;">
  <img
    src="assets/curv.gif"
    alt="PICurv simulation preview"
    style="width:100%; max-width:900px; height:auto; display:inline-block;"
  />
</div>

<div style="text-align:center; margin:1rem 0;">
  <img
    src="assets/paraview_flat_channel.png"
    alt="Example visualization output (ParaView)"
    style="width:100%; max-width:900px; height:auto; display:inline-block;"
  />
</div>
@endhtmlonly

@section methodology_sec Core Methodology

PICurv uses a two-way coupled Eulerian-Lagrangian strategy:

1. Eulerian phase: velocity/pressure evolution on a structured curvilinear grid.
2. Lagrangian phase: particle motion and particle-carried quantities.
3. Coupling:
   - Grid -> particle interpolation for advection/source evaluation.
   - Particle -> grid projection for reconstructed Eulerian particle-derived fields.

@section features_sec Key Features

- YAML-driven pipeline (`case.yml`, `solver.yml`, `monitor.yml`, `post.yml`)
- Cluster scheduler integration (`cluster.yml`) for Slurm job generation/submission
- Study/sweep orchestration (`study.yml`) with Slurm arrays, dependency chaining, and metric plots
- Multiple grid ingestion modes: `programmatic_c`, `file`, and `grid_gen`
- Direct and wrapped grid generation through `scripts/grid.gen`
- Named momentum solver strategy with solver-specific option blocks
- Structured post-processing pipelines (Eulerian, Lagrangian, statistics)
- Dedicated config contract and ingestion mapping docs for maintenance and extension

For release notes and recent contract changes, see **@subpage 18_Changelog**.
For solver-method overviews, see **@subpage 21_Methods_Overview**.
