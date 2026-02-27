@page Documentation_Map Documentation Map

This is the primary structural entry point for PICurv documentation.
Use this page instead of the raw generated page dump.

@tableofcontents

@section lifecycle_axis_sec 1. Simulation Lifecycle

@subsection lifecycle_design_sec 1.1 Design and Setup

- **@subpage 01_Installation**
- **@subpage 02_Tutorial_Programmatic_Grid**
- **@subpage 03_Tutorial_File-Based_Grid**
- **@subpage 07_Case_Reference**
- **@subpage 08_Solver_Reference**
- **@subpage 09_Monitor_Reference**
- **@subpage 10_Post_Processing_Reference**
- **@subpage 44_Boundary_Conditions_Guide**
- **@subpage 45_Particle_Initialization_and_Restart**

@subsection lifecycle_run_sec 1.2 Run and Monitor

- **@subpage 05_The_Conductor_Script**
- **@subpage 06_Simulation_Anatomy**
- **@subpage 36_Cluster_Run_Guide**
- **@subpage 37_Sweep_Studies_Guide**
- **@subpage 04_Visualization_Tutorial**

@subsection lifecycle_debug_sec 1.3 Debug and Validate

- **@subpage 39_Common_Fatal_Errors**
- **@subpage 40_Testing_and_Quality_Guide**
- **@subpage 12_Capabilities_Summary**

@section artifact_axis_sec 2. Configuration and Artifacts

@subsection artifact_yaml_sec 2.1 YAML Contracts and Ingestion

- **@subpage 14_Config_Contract**
- **@subpage 15_Config_Ingestion_Map**
- **@subpage 16_Config_Extension_Playbook**
- **@subpage 17_Workflow_Extensibility**

@subsection artifact_runtime_sec 2.2 Runtime and Artifacts

- **@subpage 33_Initial_Conditions**
- **@subpage 34_Particle_Model_Overview**
- **@subpage 45_Particle_Initialization_and_Restart**

@section methods_axis_sec 3. Numerical Methods and Models

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

@section runtime_axis_sec 4. Runtime and Source Architecture

- **@subpage 13_Code_Architecture**
- **@subpage 20_Grid_Cell_Architecture_Guide**
- **@subpage 46_C_Runtime_Execution_Map**
- **@subpage 35_API_Documentation_Status**

@section operations_axis_sec 5. Operations and Quality

- **@subpage 36_Cluster_Run_Guide**
- **@subpage 37_Sweep_Studies_Guide**
- **@subpage 39_Common_Fatal_Errors**
- **@subpage 40_Testing_and_Quality_Guide**
- **@subpage 29_Maintenance_Backlog**

@section repo_axis_sec 6. Repository and Documentation Navigation

- **@subpage 30_Repository_Navigation**
- **@subpage 18_Changelog**

@section ref_axis_sec 7. Generated Reference Views

This section is for maintainers who need code-symbol navigation views instead of prose guides.

- **File List by type (headers/sources/scripts):** `files_structured.html`
- **Data Structures by solver module:** `annotated_structured.html`
- **Raw Doxygen indices (unstructured):** `files.html`, `annotated.html`, `globals.html`
