@page 06_Simulation_Anatomy Anatomy of a Simulation: Cases, Solvers, and Monitors

A key design principle of the PICurv platform is **modularity**. A simulation is not defined by a single, monolithic configuration file. Instead, it is **composed** by mixing and matching different types of configuration files at runtime using the `pic-flow` conductor.

Understanding this concept is the key to becoming a power user. This guide explains the different roles of the configuration files and how they work together.

@tableofcontents

@section roles_sec 1. The Three Roles of Configuration Files

Every simulation run is defined by combining three distinct types of YAML files:

1.  **The Case File (`case.yml`):**
    -   **Role:** Defines the **"what"** and **"where"**.
    -   **Content:** Physics (Reynolds number), geometry (grid dimensions), boundary conditions, and run duration.
    -   **Location:** Typically unique to each study and lives inside a study directory (e.g., `my_flat_channel/case.yml`).

2.  **The Solver Profile (`solver.yml`):**
    -   **Role:** Defines the **"how"**.
    -   **Content:** Numerical strategy, solver types (implicit/explicit), tolerances, and advanced PETSc settings.
    -   **Location:** A library of reusable profiles is located in `config/solvers/`.

3.  **The Monitor Profile (`monitor.yml`):**
    -   **Role:** Defines the **"observation"**.
    -   **Content:** I/O frequency, logging verbosity, and performance profiling settings.
    -   **Location:** A library of reusable profiles is located in `config/monitors/`.

@section workflow_sec 2. The Mix-and-Match Workflow

The `pic-flow run` command requires you to explicitly provide one of each of these files. This allows for powerful combinations:

```bash
./bin/pic-flow run \
    --case    my_studies/turbulent_pipe/case.yml \
    --solver  config/solvers/Quick-Test-Explicit.yml \
    --monitor config/monitors/Debug-Output.yml \
    -n 16 --solve
```

In this example, you are running a `turbulent_pipe` simulation using a fast, low-accuracy solver for a quick test, while getting highly detailed debug output. To run the same case for a final, high-fidelity result, you would simply change the flags:

```bash
./bin/pic-flow run \
    --case    my_studies/turbulent_pipe/case.yml \
    --solver  config/solvers/High-Accuracy-Implicit.yml \
    --monitor config/monitors/Standard-Output.yml \
    -n 128 --solve
```
The physics of the problem (`case.yml`) remain identical, but the numerical method and the level of observation have been completely changed.

@section library_sec 3. The `config/` Library

The `config/` directory at the root of the project is your central library of reusable profiles. It is organized into subdirectories for each type of profile.

```
config/
├── solvers/
│   ├── Imp-MG-Standard.yml
│   └── guide.md
├── monitors/
│   ├── Standard_Output.yml
│   └── guide.md
└── postprocessors/
    ├── standard_analysis.yml
    └── guide.md
```

@subsection best_practices_ssec 3.1. Best Practices

-   **Keep it Centralized:** When you develop a new, useful solver or monitor configuration, save it in the appropriate `config/` subdirectory with a descriptive name.
-   **Don't Modify Templates:** The files in `examples/` are starting points. When you `init` a new study, they are copied to your study directory. You can modify the copies freely, but it's best to pull reusable profiles from the central `config/` library.
-   **Build Your Own Library:** Add your own `.yml` files to the `config/` directories. For example, you could create `config/solvers/CFD-Lab-Cluster-Solver.yml` with settings optimized for your specific HPC environment.

@section next_steps_sec 4. Next Steps

Now that you understand the roles of the different configuration files, you are ready to learn the specific parameters for each.

Proceed to the **@subpage 07_Case_Reference** to dive into the details of the `case.yml` file.
