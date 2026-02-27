@page 42_User_Guide_Index User Guide

This section is for day-to-day run authoring, execution, and troubleshooting.
If Getting Started proves the toolchain works, this section is where production workflow lives.

@tableofcontents

@section workflow_sec 1. Workflow and Orchestration

- **@subpage 05_The_Conductor_Script**: command model and flag contracts.
- **@subpage 06_Simulation_Anatomy**: runtime artifact graph and execution stages.
- **@subpage 36_Cluster_Run_Guide**: Slurm submission path for single runs.
- **@subpage 37_Sweep_Studies_Guide**: parameter studies and array jobs.

@section config_sec 2. Configuration References

- **@subpage 07_Case_Reference**: physics/domain/grid/BC controls.
- **@subpage 08_Solver_Reference**: numerical strategy and solver tuning.
- **@subpage 09_Monitor_Reference**: logging/profiling/output cadence.
- **@subpage 10_Post_Processing_Reference**: analysis tasks and VTK export controls.

@section practical_sec 3. Practical Recipes and Support

- **@subpage 11_User_How_To_Guides**: goal-driven recipes ("how do I ...").
- **@subpage 12_Capabilities_Summary**: feature matrix and support status.
- **@subpage 39_Common_Fatal_Errors**: known failure signatures and fixes.
- **@subpage 40_Testing_and_Quality_Guide**: smoke/quality checks before pushing changes.

@section support_sec 4. Repository Orientation

- **@subpage 30_Repository_Navigation**: directory map and linked local guides.
- **@subpage 18_Changelog**: recent behavior/contract changes.

@section outcomes_sec 5. Competencies You Should Gain

After completing this section, you should be able to:

- create or modify a case safely without breaking YAML contract assumptions,
- run local and cluster workflows with reproducible artifacts,
- execute sweeps and interpret study outputs,
- diagnose and fix common runtime/configuration failures quickly.
