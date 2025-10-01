/**
@page 10_Capabilities_Summary Capabilities Summary: What You Can Do

This page provides a high-level summary of the features and models currently implemented in PICurv that are fully configurable through the `pic-flow` workflow and `.yml` files. This is the scope of what you can achieve as a user, without needing to modify the C source code.

@tableofcontents

@section grid_cap_sec 1. Grid and Domain Features

- **Grid Generation:**
    - **Programmatic:** Generate stretched Cartesian grids for simple domains. You control the resolution, physical dimensions, and geometric stretching ratio in each direction.
    - **File-Based:** Import complex, single-block or multi-block curvilinear grids from an external file (`.picgrid` format).
- **Domain Topology:**
    - **Multi-Block:** Configure simulations with multiple stitched grid blocks, each with its own boundary conditions.
    - **Periodicity:** Enable periodic boundary conditions in any of the i, j, or k logical directions.

@section physics_cap_sec 2. Physics & Flow Models

- **Flow Regime:**
    - **Incompressible Flow:** The solver is built on the incompressible Navier-Stokes equations.
    - **Laminar or Turbulent:** Simulate both laminar flows and turbulent flows using the implemented turbulence model.
- **Dimensionality:**
    - Configure simulations to run in full **3D** or simplified **2D**.
- **Turbulence Modeling:**
    - **Large Eddy Simulation (LES):** An LES model with a dynamic Smagorinsky sub-grid scale model is available.
- **Lagrangian Particle Transport:**
    - **Particle Seeding:** Add a specified number of massless Lagrangian particles to the domain.
    - **Initialization:** Initialize particles either randomly throughout the entire volume or specifically on the surface of the primary inlet.

@section numerics_cap_sec 3. Numerical Scheme Control

- **Time Integration:**
    - **Explicit (Runge-Kutta):** A fast, explicit time-stepping scheme suitable for flows with a moderate time-step restriction.
    - **Implicit:** A more robust implicit scheme that allows for larger time steps, suitable for stiff or complex problems.
- **Convection Scheme:**
    - **QUICK:** A third-order upwind scheme (default).
    - **Central Differencing:** A second-order central difference scheme.
- **Pressure-Poisson Solver:**
    - **Geometric Multigrid:** A highly configurable geometric multigrid (PCMG) preconditioner is the default. You can control the number of levels, sweeps, and cycle types.
    - **Advanced PETSc Control:** Pass-through any valid PETSc KSP (solver) or PC (preconditioner) option to fully customize the linear algebra strategy.

@section bc_cap_sec 4. Boundary Conditions

The following physical boundary conditions are available through the `boundary_conditions` section of `case.yml`:

-   **Inlets:**
    -   `constant_velocity`: Uniform velocity profile.
    -   `parabolic`: Parabolic (Poiseuille) profile for developed flow.
-   **Outlets:**
    -   `conservation`: A zero-gradient condition that enforces mass conservation.
-   **Walls:**
    -   `noslip`: Standard no-slip wall condition.
-   **Other:**
    -   `symmetry_plane`: A slip-wall condition.
    -   `periodic`: For use with periodic domains.

@section pp_cap_sec 5. Post-Processing Capabilities

The built-in post-processor allows you to perform the following analyses without external tools:

-   **Data Conversion:**
    -   Convert cell-centered data to node-centered data (`CellToNodeAverage`), which is ideal for smooth contour plots in visualization software.
    -   Convert all non-dimensional solver output back to physical, dimensional units.
-   **Derived Quantities (Eulerian):**
    -   Calculate the **Q-Criterion** (`ComputeQCriterion`) to identify vortex structures.
    -   Normalize fields relative to a reference point (`NormalizeRelativeField`).
-   **Derived Quantities (Lagrangian):**
    -   Calculate the **Specific Kinetic Energy** (`specific_ke`) of each particle.

@section next_steps_sec 6. Beyond Configuration: The Developer Portal

If your research requires a feature not listed above—such as a new turbulence model, a custom boundary condition, or a new post-processing kernel—you will need to extend the C source code.

Please proceed to the **@subpage 13_Code_Architecture "Developer Portal"** to learn about the code's structure and how to contribute new features.
