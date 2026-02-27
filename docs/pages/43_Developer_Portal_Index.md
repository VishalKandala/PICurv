@page 43_Developer_Portal_Index Developer Portal

This section is for maintainers and contributors changing solver behavior, YAML contracts, or workflow logic.
It emphasizes architecture boundaries, method-level reasoning, and safe extension points.

@tableofcontents

@section architecture_sec 1. Architecture and Contracts

- **@subpage 13_Code_Architecture**: module boundaries and runtime ownership model.
- **@subpage 14_Config_Contract**: YAML -> generated artifact -> runtime contract.
- **@subpage 15_Config_Ingestion_Map**: where specific YAML keys enter C/Python paths.
- **@subpage 16_Config_Extension_Playbook**: adding new keys/features safely.
- **@subpage 17_Workflow_Extensibility**: extending orchestration beyond current run/sweep modes.
- **@subpage 46_C_Runtime_Execution_Map**: startup/timestep execution trace across C modules.

@section methods_sec 2. Numerical Methods and Models

- **@subpage 21_Methods_Overview**: entry map across all numerical method pages.
- **@subpage 22_CURVIB_Method** through **@subpage 28_IEM_and_Statistical_Averaging**: flow solve, pressure solve, coupling, and mixing/statistics internals.
- **@subpage 31_Momentum_Solvers**: runtime-dispatched momentum solver implementations.
- **@subpage 32_Analytical_Solutions** and **@subpage 33_Initial_Conditions**: initialization and analytical-mode behavior.
- **@subpage 44_Boundary_Conditions_Guide**: BC handler mapping from YAML to C handler objects.
- **@subpage 45_Particle_Initialization_and_Restart**: seeding and restart/migration behavior.
- **@subpage 34_Particle_Model_Overview**: particle lifecycle and projection coupling.

@section maintenance_sec 3. Documentation and Maintenance

- **@subpage 35_API_Documentation_Status**: current API documentation coverage and gaps.
- **@subpage 20_Grid_Cell_Architecture_Guide**: data-structure contracts in grid/cell storage.
- **@subpage 29_Maintenance_Backlog**: low-priority but known technical-debt items.

@section contribution_flow_sec 4. Suggested Contributor Read Path

1. **@subpage 13_Code_Architecture**
2. **@subpage 14_Config_Contract**
3. **@subpage 15_Config_Ingestion_Map**
4. **@subpage 46_C_Runtime_Execution_Map**
5. **@subpage 21_Methods_Overview**
6. **@subpage 16_Config_Extension_Playbook**

@section developer_outcomes_sec 5. Expected Outcomes

After working through this section, you should be able to:

- trace a new YAML key from schema to runtime consumer,
- identify the right C module for a numerical feature change,
- update docs/tests/validation alongside code changes,
- preserve compatibility and diagnostics while extending behavior.
