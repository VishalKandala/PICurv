@page 42_User_Guide_Index User Guide

@anchor _User_Guide_Index

@pagemeta{Hub, Users authoring runs, Routing only - see linked pages for detail}

This section is for day-to-day run authoring, execution, and troubleshooting.
If Getting Started proves the toolchain works, this section is where production workflow lives.

@tableofcontents

@section p42_workflow_sec 1. Workflow and Orchestration

- **@subpage 05_The_Conductor_Script**: command model and flag contracts.
- **@subpage 49_Workflow_Recipes_and_Config_Cookbook**: practical run patterns and modular profile combinations.
- **@subpage 06_Simulation_Anatomy**: runtime artifact graph and execution stages.
- **@subpage 52_Run_Artifact_Lifecycle_Contract**: operational path from new run to restart, post-only reuse, and batch job generation.
- **@subpage 61_Storage_Management_Guide**: recurring run/study archival, offload, catalog, and restore workflows.
- **@subpage 36_Cluster_Run_Guide**: Slurm submission path for single runs.
- **@subpage 37_Sweep_Studies_Guide**: parameter studies and array jobs.

@section p42_config_sec 2. Configuration References

- **@subpage 07_Case_Reference**: physics/domain/grid/BC controls.
- **@subpage 72_LES_Turbulence_Closure**: what the subgrid models do, how to choose their averaging and limiting, and how to read the coefficient diagnostics.
- **@subpage 08_Solver_Reference**: numerical strategy and solver tuning.
- **@subpage 09_Monitor_Reference**: logging/profiling/diagnostics/output cadence.
- **@subpage 10_Post_Processing_Reference**: analysis tasks and VTK export controls.
- **@subpage 48_Grid_Generator_Guide**: direct `grid.gen` usage and `grid_gen` wrapper semantics.
- **@subpage 44_Boundary_Conditions_Guide**: detailed BC handler options and validation rules.
- **@subpage 45_Particle_Initialization_and_Restart**: particle seeding/restart behavior by mode.

@section p42_journeys_sec 3. Task Routes

- **@subpage 65_Example_Catalog**: every shipped example, what it demonstrates, and what evidence it provides.
- **@subpage 70_Case_Design_Guide**: building a case from scratch, in the order the decisions constrain each other.
- **@subpage 67_Troubleshooting**: symptom-driven diagnosis, cheapest checks first.
- **@subpage 66_Evidence_Matrix**: declared evidence sources for every public capability - sources, not verified results.
- **@subpage 68_Glossary**: PICurv terms and renamed spellings.

@section p42_practical_sec 4. Practical Recipes and Support

- **@subpage 11_User_How_To_Guides**: goal-driven recipes ("how do I ...").
- **@subpage 12_Capabilities_Summary**: feature matrix and support status.
- **@subpage 39_Common_Fatal_Errors**: known failure signatures and fixes.
- **@subpage 40_Testing_and_Quality_Guide**: smoke/quality checks before pushing changes.

@section p42_support_sec 5. Repository Orientation

- **@subpage 30_Repository_Navigation**: directory map and linked local guides.
- **@subpage 18_Changelog**: recent behavior/contract changes.
- **@subpage Documentation_Map**: categorized index across lifecycle and artifact types.

@section p42_outcomes_sec 6. Competencies You Should Gain

After completing this section, you should be able to:

- create or modify a case safely without breaking YAML contract assumptions,
- run local and cluster workflows with reproducible artifacts,
- execute sweeps and interpret study outputs,
- diagnose and fix common runtime/configuration failures quickly.
