@page Documentation_Map Documentation Map

@anchor _Map

This is the primary structural entry point for PICurv documentation.
Use this page instead of the raw generated page dump.

@tableofcontents

@section p47_lifecycle_axis_sec 1. Simulation Lifecycle

@subsection p47_lifecycle_design_sec 1.1 Design and Setup

- **@subpage 01_Installation**
- **@subpage 02_Tutorial_Programmatic_Grid**
- **@subpage 03_Tutorial_File-Based_Grid**
- **@subpage 07_Case_Reference**
- **@subpage 08_Solver_Reference**
- **@subpage 09_Monitor_Reference**
- **@subpage 10_Post_Processing_Reference**
- **@subpage 48_Grid_Generator_Guide**
- **@subpage 44_Boundary_Conditions_Guide**
- **@subpage 45_Particle_Initialization_and_Restart**

@subsection p47_lifecycle_run_sec 1.2 Run and Monitor

- **@subpage 05_The_Conductor_Script**
- **@subpage 49_Workflow_Recipes_and_Config_Cookbook**
- **@subpage 06_Simulation_Anatomy**
- **@subpage 36_Cluster_Run_Guide**
- **@subpage 37_Sweep_Studies_Guide**
- **@subpage 04_Visualization_Tutorial**

@subsection p47_lifecycle_debug_sec 1.3 Debug and Validate

- **@subpage 39_Common_Fatal_Errors**
- **@subpage 40_Testing_and_Quality_Guide**
- **@subpage 12_Capabilities_Summary**

@section p47_artifact_axis_sec 2. Configuration and Artifacts

@subsection p47_artifact_yaml_sec 2.1 YAML Contracts and Ingestion

- **@subpage 14_Config_Contract**
- **@subpage 15_Config_Ingestion_Map**
- **@subpage 16_Config_Extension_Playbook**
- **@subpage 17_Workflow_Extensibility**
- **@subpage 48_Grid_Generator_Guide**
- **@subpage 49_Workflow_Recipes_and_Config_Cookbook**

@subsection p47_artifact_runtime_sec 2.2 Runtime and Artifacts

- **@subpage 33_Initial_Conditions**
- **@subpage 34_Particle_Model_Overview**
- **@subpage 45_Particle_Initialization_and_Restart**

@section p47_methods_axis_sec 3. Numerical Methods and Models

- **@subpage 21_Methods_Overview**
- **@subpage 22_CURVIB_Method**
- **@subpage 23_Fractional_Step_Method**
- **@subpage 24_Dual_Time_Picard_RK4**
- **@subpage 25_Pressure_Poisson_GMRES_Multigrid**
- **@subpage 26_Walking_Search_Method**
- **@subpage 27_Trilinear_Interpolation_and_Projection**
- **@subpage 28_IEM_and_Statistical_Averaging**
- **@subpage 31_Momentum_Solvers**
- **@subpage 32_Analytical_Solutions**

@section p47_runtime_axis_sec 4. Runtime and Source Architecture

- **@subpage 13_Code_Architecture**
- **@subpage 20_Grid_Cell_Architecture_Guide**
- **@subpage 46_C_Runtime_Execution_Map**
- **@subpage 35_API_Documentation_Status**

@section p47_operations_axis_sec 5. Operations and Quality

- **@subpage 36_Cluster_Run_Guide**
- **@subpage 37_Sweep_Studies_Guide**
- **@subpage 39_Common_Fatal_Errors**
- **@subpage 40_Testing_and_Quality_Guide**
- **@subpage 29_Maintenance_Backlog**

@section p47_repo_axis_sec 6. Repository and Documentation Navigation

- **@subpage 30_Repository_Navigation**
- **@subpage 18_Changelog**

@section p47_ref_axis_sec 7. Generated Reference Views

This section is for maintainers who need code-symbol navigation views instead of prose guides.

- **File List by type (headers/sources/scripts):** `files_structured.html`
- **Data Structures by solver module:** `annotated_structured.html`
- **Raw Doxygen indices (unstructured):** `files.html`, `annotated.html`, `globals.html`

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Documentation Map** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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

