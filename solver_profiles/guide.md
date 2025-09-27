# PIC-Flow Solver Profiles Library

## 1. Purpose of this Directory

This directory contains a library of reusable **Solver Profiles** (`.yml` files). Each file in this directory defines a complete numerical strategy for the C-solver.

The core principle of the PIC-Flow platform is **modularity**. You can combine any `case.yml` with any `solver.yml` at runtime. This allows you to test different numerical methods on the same physical problem without ever changing the case definition.

A solver profile controls:
- The core numerical scheme (e.g., implicit vs. explicit).
- Convergence tolerances for iterative solvers.
- The configuration of the pressure-Poisson solver (including multigrid settings).
- Any advanced PETSc command-line options.

## 2. Included Profiles

- **`Imp-MG-Standard.yml`**: A robust, general-purpose implicit solver that uses a geometric multigrid preconditioner. It is a well-balanced choice for a wide variety of incompressible flow problems and is the recommended starting point for new cases.
- **`master_solver_template.yml`**: A heavily commented reference file that showcases **all possible** solver configuration options. Use this as a guide or a source for copy-pasting advanced sections, but do not use it to run a simulation directly.

*(As you develop more solvers, add them to this list, e.g., `Explicit-RK3.yml`, `Robust-FSI-Solver.yml`, etc.)*

## 3. How to Use

When you launch a simulation, you specify which solver profile to use with the `--solver` flag:

```bash
# Example: Using the standard implicit solver for a specific case
./bin/pic.flow run \
    --case my_studies/my_case.yml \
    --solver solver_profiles/Imp-MG-Standard.yml \
    --monitor ...
```

## 4. Creating a New Solver Profile

You should create a new solver profile whenever you want to test or save a distinct numerical strategy.

**Workflow:**
1.  **Copy an existing profile:** The easiest way to start is to copy a similar profile, like `Imp-MG-Standard.yml`.
2.  **Give it a descriptive name:** The filename should clearly state its purpose (e.g., `Fine-Tolerance-GMRES.yml`, `Quick-Test-Jacobi.yml`).
3.  **Edit the contents:** Modify the settings to implement your new strategy. For example, you might tighten the `tolerances`, change the `ksp_type` in the `petsc_passthrough_options`, or adjust the number of `multigrid` levels.
4.  **Save and use:** Save the new file in this directory. It is now part of your library and can be used with any case via the `--solver` flag.