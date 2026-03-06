@page 42_User_Guide_Index User Guide

@anchor _User_Guide_Index

This section is for day-to-day run authoring, execution, and troubleshooting.
If Getting Started proves the toolchain works, this section is where production workflow lives.

@tableofcontents

@section p42_workflow_sec 1. Workflow and Orchestration

- **@subpage 05_The_Conductor_Script**: command model and flag contracts.
- **@subpage 49_Workflow_Recipes_and_Config_Cookbook**: practical run patterns and modular profile combinations.
- **@subpage 06_Simulation_Anatomy**: runtime artifact graph and execution stages.
- **@subpage 36_Cluster_Run_Guide**: Slurm submission path for single runs.
- **@subpage 37_Sweep_Studies_Guide**: parameter studies and array jobs.

@section p42_config_sec 2. Configuration References

- **@subpage 07_Case_Reference**: physics/domain/grid/BC controls.
- **@subpage 08_Solver_Reference**: numerical strategy and solver tuning.
- **@subpage 09_Monitor_Reference**: logging/profiling/output cadence.
- **@subpage 10_Post_Processing_Reference**: analysis tasks and VTK export controls.
- **@subpage 48_Grid_Generator_Guide**: direct `grid.gen` usage and `grid_gen` wrapper semantics.
- **@subpage 44_Boundary_Conditions_Guide**: detailed BC handler options and validation rules.
- **@subpage 45_Particle_Initialization_and_Restart**: particle seeding/restart behavior by mode.

@section p42_practical_sec 3. Practical Recipes and Support

- **@subpage 11_User_How_To_Guides**: goal-driven recipes ("how do I ...").
- **@subpage 12_Capabilities_Summary**: feature matrix and support status.
- **@subpage 39_Common_Fatal_Errors**: known failure signatures and fixes.
- **@subpage 40_Testing_and_Quality_Guide**: smoke/quality checks before pushing changes.

@section p42_support_sec 4. Repository Orientation

- **@subpage 30_Repository_Navigation**: directory map and linked local guides.
- **@subpage 18_Changelog**: recent behavior/contract changes.
- **@subpage Documentation_Map**: categorized index across lifecycle and artifact types.

@section p42_outcomes_sec 5. Competencies You Should Gain

After completing this section, you should be able to:

- create or modify a case safely without breaking YAML contract assumptions,
- run local and cluster workflows with reproducible artifacts,
- execute sweeps and interpret study outputs,
- diagnose and fix common runtime/configuration failures quickly.

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **User Guide** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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

