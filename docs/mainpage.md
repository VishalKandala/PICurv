@mainpage PICurv Solver Documentation

@section intro_sec Introduction

# 🌀 PICurv: A Hybrid Eulerian-Lagrangian Framework for Scalar Transport

Welcome to the documentation for **PICurv**, a high-performance, parallel framework for simulating turbulent flows with actively coupled scalar fields using a hybrid Eulerian-Lagrangian method.

PICurv is designed to tackle complex scalar transport problems where traditional grid-based methods suffer from numerical diffusion. It couples a curvilinear immersed boundary (CurvIB) fluid dynamics solver with a parallel Lagrangian particle method. In this framework, the "particles" act as moving computational points that carry and evolve scalar properties. These properties are then projected back onto the Eulerian grid, creating a powerful two-way coupling. This makes PICurv an ideal platform for advanced simulations in turbulent mixing and combustion, such as transported Probability Density Function (t-PDF) or Flamelet/Progress Variable (FPV) models.

@section methodology_sec Core Methodology

PICurv's methodology is a tightly integrated, two-way coupled hybrid Eulerian-Lagrangian scheme. It solves for the fluid mechanics on a stationary grid while tracking the evolution of scalar fields on a vast number of moving Lagrangian markers.

### 1. The Eulerian Phase (The Grid)

The background fluid flow and the grid-based representation of scalar fields are handled in the Eulerian frame.
- **Flow Solver:** The velocity field is solved using the **CurvIB solver**, which employs a pressure-based **fractional-step projection method** for the incompressible Navier-Stokes equations on a structured, body-fitted curvilinear grid.
- **Turbulence Modeling:** Large Eddy Simulation (LES) with a dynamic Smagorinsky sub-grid scale model is used to efficiently capture turbulent flow dynamics.
- **Eulerian Scalar Fields (`phi_field`):** The grid also hosts Eulerian fields for the scalars being tracked (e.g., mixture fraction, temperature). These fields are not directly advected but are instead constructed from the properties of the Lagrangian markers.

### 2. The Lagrangian Phase (The Markers)

The core of the scalar transport model resides in the Lagrangian phase, which consists of millions of computational markers.
- **Computational Markers:** These are abstract, massless markers that act as carriers of scalar information (`phi_particle`). Each marker represents a statistical sample of the flow.
- **Scalar Evolution:** An evolution equation (e.g., a transport or reaction equation) is solved for the scalar properties along each marker's trajectory: `d(phi_particle)/dt = S(phi_particle, ...)`. The source term `S` can depend on local fluid properties interpolated from the grid.
- **Parallel Management:** Marker data and evolution are managed by PETSc's `DMSwarm`, handling memory, data layout, and automatic migration between MPI ranks.

### 3. The Two-Way Coupling Mechanism

The true power of PICurv lies in the continuous, two-way exchange of information between the grid and the markers, using specific numerical schemes.

- **Grid → Marker (Interpolation):**
    1.  The fluid velocity `**u**` is interpolated to each marker's location using a **trilinear interpolation** scheme to advect it accurately.
    2.  Other necessary Eulerian field variables (e.g., velocity gradients, turbulent dissipation) are also interpolated to the markers using the same method to evaluate the source terms in the scalar evolution equation.

- **Marker → Grid (Projection/Deposition):**
    1.  The Eulerian scalar field (`phi_field`) is reconstructed from the Lagrangian marker properties using **cell averaging**. The value of `phi_field` for a given grid cell is computed as the average of the `phi_particle` values from all markers currently residing within that cell's volume.
    2.  This projection step allows the high-resolution, non-diffusive information evolved in the Lagrangian frame to directly influence the Eulerian field, which can then be used for visualization or to affect other physics in the simulation.

@section features_sec Key Features

- **Hybrid Eulerian-Lagrangian Solver:** A modern framework for high-fidelity scalar transport in turbulent flows.
- **Advanced CFD for Transport:** Employs a pressure-based fractional-step method and LES to provide an accurate velocity field for Lagrangian advection.
- **Two-Way Scalar Coupling:** Implements a robust mechanism for projecting Lagrangian scalar information back onto the Eulerian grid.
- **Flexible Scalar Evolution:** The framework is designed to solve arbitrary user-defined transport equations on the Lagrangian markers.
- **Efficient Coupling Kernels:** Employs **trilinear interpolation** for grid-to-marker data transfer and **cell averaging** for marker-to-grid projection, ensuring accurate and robust data exchange.
- **Parallel by Design:** Built on PETSc's `DMDA` and `DMSwarm` to ensure excellent performance and scalability on HPC systems.

@section arch_sec Architecture and Program Flow

The solver follows a clear, five-stage execution path: **Initialize, Configure, Setup, Execute, and Finalize**, orchestrated by `main()`. All simulation data is managed by a central `SimulationContext` object.

@image html high_level_architecture.png "High-Level Simulation Workflow" width=800px

@section nav_links Next Steps

- **@subpage getting_started**: Compile the code and run your first simulation.
- **@subpage user_guide**: Learn how to configure simulations, use command-line options, and interpret results.
- **@subpage developer_guide**: Understand the code's architecture and learn how to extend it with new models.
- **@subpage theory_and_numerics**: Read the detailed description of the governing equations and numerical methods.
