@page 43_Developer_Portal_Index Developer Portal

@anchor _Developer_Portal_Index

@pagemeta{Hub, Contributors and maintainers, Routing only - see linked pages for detail}

This section is for maintainers and contributors changing solver behavior, YAML contracts, or workflow logic.
It emphasizes architecture boundaries, method-level reasoning, and safe extension points.

@tableofcontents

@section p43_architecture_sec 1. Architecture and Contracts

- **@subpage 13_Code_Architecture**: module boundaries and runtime ownership model.
- **@subpage 14_Config_Contract**: YAML -> generated artifact -> runtime contract.
- **@subpage 15_Config_Ingestion_Map**: where specific YAML keys enter C/Python paths.
- **@subpage 16_Config_Extension_Playbook**: adding new keys/features safely.
- **@subpage 50_Modular_Selector_Extension_Guide**: selector-by-selector hook points for extending current user-facing options.
- **@subpage 17_Workflow_Extensibility**: extending orchestration beyond current run/sweep modes.
- **@subpage 46_C_Runtime_Execution_Map**: startup/timestep execution trace across C modules.
- **@subpage 56_Field_Identity_and_Layout_Catalog**: typed Eulerian/particle field identities, runtime views, and coupling metadata.
- **@subpage 58_Field_Statistics**: scientific field statistics — accumulation, monitoring, checkpoint state, derived output, and extension points.
- **@subpage 57_Future_Architecture_Specifications**: status and sequencing of proposed future architecture.
- **@subpage 60_Field_Statistics_Planned_Extensions**: spatial targets, further products, and histories, with their dependency order.
- **@subpage 59_Function_Identity_and_Observability_Specification**: deferred, benchmark-gated logging/profiling identity design.

@section p43_contracts_sec 2. Documentation Contracts

- **@subpage 64_Documentation_Extension_Framework**: capability, family, and subsystem contracts, lifecycle gates, and concern modules.
- **@subpage 63_Page_Type_Contract**: which page type you are writing, and what it must refuse to do.
- **@subpage 62_Capability_Status_Vocabulary**: status words, the known-defect disclosure rule, and evidence facets.
- **@subpage 69_Scientific_Writing_Handoff**: what the foundational method pages owe, and the scaffolding already in place.

@section p43_methods_sec 3. Numerical Methods and Models

- **@subpage 21_Methods_Overview**: entry map across all numerical method pages.
- **@subpage 22_CURVIB_Method** through **@subpage 28_IEM_and_Statistical_Averaging**: flow solve, pressure solve, coupling, and mixing/statistics internals.
- **@subpage 31_Momentum_Solvers**: runtime-dispatched momentum solver implementations.
- **@subpage 32_Analytical_Solutions** and **@subpage 33_Initial_Conditions**: initialization and analytical-mode behavior.
- **@subpage 44_Boundary_Conditions_Guide**: BC handler mapping from YAML to C handler objects.
- **@subpage 45_Particle_Initialization_and_Restart**: seeding and restart/migration behavior.
- **@subpage 34_Particle_Model_Overview**: particle lifecycle and projection coupling.

@section p43_maintenance_sec 4. Documentation and Maintenance

- **@subpage 35_API_Documentation_Status**: current API documentation coverage and gaps.
- **@subpage 20_Grid_Cell_Architecture_Guide**: data-structure contracts in grid/cell storage.
- **@subpage 29_Maintenance_Backlog**: low-priority but known technical-debt items.
- **@subpage Documentation_Map**: alternate categorized entry path.

@section p43_contribution_flow_sec 5. Suggested Contributor Read Path

1. **@subpage 13_Code_Architecture**
2. **@subpage 14_Config_Contract**
3. **@subpage 15_Config_Ingestion_Map**
4. **@subpage 16_Config_Extension_Playbook**
5. **@subpage 50_Modular_Selector_Extension_Guide**
6. **@subpage 46_C_Runtime_Execution_Map**
7. **@subpage 56_Field_Identity_and_Layout_Catalog**
8. **@subpage 58_Field_Statistics**
9. **@subpage 21_Methods_Overview**

@section p43_developer_outcomes_sec 6. Expected Outcomes

After working through this section, you should be able to:

- trace a new YAML key from schema to runtime consumer,
- identify the right C module for a numerical feature change,
- update docs/tests/validation alongside code changes,
- preserve diagnostics and contract clarity while extending behavior.
