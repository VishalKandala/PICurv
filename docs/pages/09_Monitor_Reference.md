@page 09_Monitor_Reference Configuration Reference: Monitor Profiles (`monitor.yml`)

For a complete, heavily commented reference file showing every possible option, please see the master template:

@verbinclude master_template/master_monitor.yml

A **Monitor Profile** is a `.yml` file that defines a strategy for observing and recording a simulation run. Changing a monitor profile **does not change the scientific result** of the simulation; it only changes what you see on the console, what performance data is collected, and how often results are saved to disk.

This modularity allows you to run the exact same case with different monitor profiles for different purposes, such as production runs, detailed debugging, or performance analysis.

This document serves as a reference for all available sections and parameters within a `monitor.yml` file.

@tableofcontents

@section io_sec 1. The `io` Section

This section controls the frequency and location of all file-based input and output.

```yaml
io:
  data_output_frequency: 100
  particle_log_interval: 10
  directories:
    output: "results"
    restart: "results"
    log: "logs"
```

| Parameter | Type | Description | C-Solver Flag |
| :--- | :--- | :--- | :--- |
| `data_output_frequency`| Integer | The solver will save all binary output files (Eulerian and Lagrangian) every N time steps. | `-tio` |
| `particle_log_interval`| Integer | The solver will print a detailed table of particle properties to the console every N time steps (only if verbosity is INFO or higher). | `-logfreq` |
| `directories.output` | String | The name of the subdirectory within the run folder where simulation results (`ufield`, `pfield`, etc.) will be saved. | `-output_dir` |
| `directories.restart`| String | The name of the subdirectory where the solver will look for restart files. | `-restart_dir` |
| `directories.log` | String | The name of the subdirectory where text-based log files (e.g., convergence history) will be saved. | `-log_dir` |

@section logging_sec 2. The `logging` Section

This section gives you fine-grained control over the console output for debugging and monitoring.

```yaml
logging:
  verbosity: "INFO"
  enabled_functions:
    - main
    - CreateSimulationContext
```

| Parameter | Type | Description | Environment Variable / C-Solver Flag |
| :--- | :--- | :--- | :--- |
| `verbosity` | String | Sets the maximum level of detail for console messages. Messages with a lower severity are always shown. The levels are, in order of increasing detail: `ERROR`, `WARNING`, `INFO` (default), `DEBUG`, `PROFILE`. | `LOG_LEVEL` |
| `enabled_functions`| String List | A list of C function names. If this list is provided, **only** log messages originating from these specific functions will be printed to the console. This is an extremely powerful tool for targeted debugging. | `-whitelist_config_file` |

@subsection log_levels_ssec 2.1. Log Verbosity Levels

-   **`ERROR`**: Only shows critical, simulation-ending errors.
-   **`WARNING`**: Shows errors and potential issues that do not halt the simulation.
-   **`INFO`**: (Default) Shows high-level progress updates, such as which stage of the solver is running.
-   **`DEBUG`**: Shows highly detailed, low-level output from deep within functions, intended for developers. This can produce a very large amount of text.
-   **`PROFILE`**: Shows messages exclusively related to the performance profiling system.

@section profiling_sec 3. The `profiling` Section

This section controls the built-in performance profiling system, which measures the wall-clock time spent in key C functions.

```yaml
profiling:
  critical_functions:
    - Flow_Solver
    - AdvanceSimulation
```

| Parameter | Type | Description | C-Solver Flag |
| :--- | :--- | :--- | :--- |
| `critical_functions`| String List | A list of C function names that are considered "critical" for performance. The time spent in these functions will be logged at the `INFO` verbosity level at the end of each time step. A full report of all profiled functions is always available at the `PROFILE` verbosity level. | `-profile_config_file` |

@section next_steps_sec 4. Next Steps

You now have a complete understanding of the main configuration files for running a simulation. The final piece of the user workflow is learning how to analyze the data you've generated.

Proceed to the **@subpage 10_Post_Processing_Reference** to learn about all the parameters available in the `post.yml` file for creating analysis recipes.