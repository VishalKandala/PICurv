/**
@page developer_guide Developer Guide

This guide is for developers who wish to understand, modify, or extend the PICurv codebase.

@section arch_overview High-Level Architecture
The solver is built around a central `SimulationContext` (simCtx) object that holds all simulation-wide data and configuration. The program flow follows five main stages as defined in `main()`: Initialize, Configure, Setup, Execute, and Finalize.

@image html high_level_architecture.png "High-Level Simulation Workflow" width=800px

@section program_flow Detailed Program Flow
The `main()` function orchestrates the simulation.
1.  **Initialize**: `PetscInitialize()` sets up the MPI and PETSc environment.
2.  **Configure**: `CreateSimulationContext()` parses all command-line arguments into the `simCtx` struct.
3.  **Setup**: A sequence of setup routines prepares the simulation environment.
    - `SetupGridAndSolvers()`: Builds the multi-grid hierarchy and PETSc solvers.
    - ... (list all setup functions)
4.  **Execute**:
    - `PerformInitialSetup()`: One-time operations before the time loop.
    - `AdvanceSimulation()`: The main time-stepping loop.
5.  **Finalize**: `PetscFinalize()` cleans up resources.

@section data_structures Core Data Structures
- @ref SimCtx : The most important data structure in the solver.
- (Link to other key structs)

@section api_groups API Modules
The codebase is organized into logical modules for clarity.
- @ref setup_routines
- @ref execution_routines
- @ref particle_management
- ... (etc.)

*/
