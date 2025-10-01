@mainpage PICurv Solver Documentation

@section intro_sec Introduction

# 🌀 PICurv: A Hybrid Eulerian-Lagrangian Framework for Scalar Transport

Welcome to the documentation for **PICurv**, a high-performance, parallel framework for simulating turbulent flows with actively coupled scalar fields using a hybrid Eulerian-Lagrangian method.

PICurv is designed to tackle complex scalar transport problems where traditional grid-based methods suffer from numerical diffusion. It couples a curvilinear immersed boundary (CurvIB) fluid dynamics solver with a parallel Lagrangian particle method. In this framework, "particles" act as moving computational points that carry and evolve scalar properties. These properties are then projected back onto the Eulerian grid, creating a powerful two-way coupling. This makes PICurv an ideal platform for advanced simulations in turbulent mixing and combustion, such as transported Probability Density Function (t-PDF) or Flamelet/Progress Variable (FPV) models.

---

@section get_started_main To Begin, Choose Your Path:

<div class="main_page_buttons">
    <a href="01_Installation.html" class="main_page_button">
        <h2>Getting Started</h2>
        <p>For new users. Install the code, run your first example simulation, and visualize the results in minutes.</p>
    </a>
    <a href="04_Conductor_Script.html" class="main_page_button">
        <h2>User Guide</h2>
        <p>For power users. Learn how to configure your own simulations using the powerful YAML-based workflow.</p>
    </a>
    <a href="09_Code_Architecture.html" class="main_page_button">
        <h2>Developer Portal</h2>
        <p>For contributors. Understand the C++ architecture and learn how to extend the solver with new models.</p>
    </a>
</div>

---

@section methodology_sec Core Methodology

PICurv's methodology is a tightly integrated, two-way coupled hybrid Eulerian-Lagrangian scheme. It solves for the fluid mechanics on a stationary grid while tracking the evolution of scalar fields on a vast number of moving Lagrangian markers.

@image html high_level_architecture.png "High-Level Simulation Workflow" width=800px

### 1. The Eulerian Phase (The Grid)

The background fluid flow is handled by the **CurvIB solver**, which employs a pressure-based **fractional-step projection method** for the incompressible Navier-Stokes equations. Turbulence is modeled using Large Eddy Simulation (LES). The grid also hosts Eulerian scalar fields (`phi_field`), which are constructed from the properties of the Lagrangian markers.

### 2. The Lagrangian Phase (The Markers)

The core of the scalar transport model resides in the Lagrangian phase, which consists of millions of computational markers. These are abstract, massless markers that act as carriers of scalar information (`phi_particle`). An evolution equation is solved for the scalar properties along each marker's trajectory, with parallel management handled by PETSc's `DMSwarm`.

### 3. The Two-Way Coupling Mechanism

The true power of PICurv lies in the continuous, two-way exchange of information between the grid and the markers.

-   **Grid → Marker (Interpolation):** Fluid velocity and other Eulerian properties are interpolated to each marker's location to advect it and evaluate source terms in its governing equation.
-   **Marker → Grid (Projection):** The Eulerian scalar field is reconstructed from the Lagrangian marker properties using cell averaging. This allows the high-resolution, non-diffusive information from the particles to directly influence the grid-based fields.

@section features_sec Key Features

- **Hybrid Eulerian-Lagrangian Solver:** A modern framework for high-fidelity scalar transport in turbulent flows.
- **YAML-Driven Workflow:** User-friendly configuration managed by a powerful Python conductor script (`pic.flow`).
- **Advanced CFD for Transport:** Employs a pressure-based fractional-step method and LES to provide an accurate velocity field for Lagrangian advection.
- **Two-Way Scalar Coupling:** Implements a robust mechanism for projecting Lagrangian scalar information back onto the Eulerian grid.
- **Flexible Scalar Evolution:** The framework is designed to solve arbitrary user-defined transport equations on the Lagrangian markers.
- **Parallel by Design:** Built on PETSc's `DMDA` and `DMSwarm` to ensure excellent performance and scalability on HPC systems.
