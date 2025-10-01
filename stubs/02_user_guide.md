@page user_guide User Guide

This guide provides practical information on how to set up the run environment, configure a simulation, and interpret the output from the PICurv solver.

@section run_env_sec 1. Setting Up the Run Environment

The PICurv solver expects a specific directory structure to run properly. It is strongly recommended to perform simulations in a separate directory outside of the source code repository.

A typical run directory should look like this:
```
my_simulation/
├── config/         # Contains all input and configuration files.
├── logs/           # Destination for runtime text logs.
├── results/        # Destination for binary output data.
└── picsolver       # Symbolic link to the solver executable.
```

Here is a step-by-step guide to create this environment:
1.  **Create a Run Directory:**
2.  
    ```bash
    mkdir my_simulation
    cd my_simulation
    ```
3.  **Create Subdirectories:**
    The `logs` and `results` directories must exist before starting a simulation.
    
    ```bash
    mkdir logs results
    ```
    
5.  **Copy a Configuration:**
    The easiest way to start is by using one of the examples from the `test/` directory in the source code.
    
    ```bash
    # Assuming 'my_simulation' and 'picurv' are in the same parent directory
    cp -r ../picurv/test/channel_flow/config ./
    ```
7.  **Link the Executable:**
    Create a symbolic link to the `picsolver` executable that you built.
    ```bash
    # Adjust the relative path to your picurv installation
    ln -s ../picurv/bin/picsolver .
    ```

@section dir_breakdown_sec 2. Directory and File Breakdown

@subsection config_dir_sec The `config/` Directory
This directory is the heart of your simulation setup. It contains all the necessary input files that define the problem. The solver will look for the following files here:

-   `control.dat`: The main control file. It contains essential simulation parameters like time step size, number of steps to run, I/O frequency, and flags for runtime behavior (e.g., whether to generate the grid programmatically).
-   `bcs.dat`: The boundary conditions file. It specifies the type and values of boundary conditions (e.g., inlet, outlet, wall) for all faces of the computational domain.
-   `grid.dat`: The grid file, containing the coordinates of the mesh. **Note:** This file is optional. If the grid is generated programmatically (as specified in `control.dat`), this file is not needed.
-   `whitelist.dat`: An optional file used for debugging. It contains a list of function names, one per line. If this file is present, only log messages originating from these functions will be printed to the console.

@subsection logs_dir_sec The `logs/` Directory
This directory is where the solver writes all runtime text-based output. The solver requires this directory to be present at startup. It will contain files such as:
-   Time histories of solver residuals (e.g., for momentum and pressure-Poisson equations).
-   Continuity error metrics at each time step.
-   Performance and profiling information.
-   General simulation progress and informational messages.

@subsection results_dir_sec The `results/` Directory
This directory stores the binary output data from the simulation, which is used for post-processing and visualization. Each file is typically named with the field and the time step number. Examples include:
-   `ufield_001000.bin`: Eulerian velocity field at step 1000.
-   `pfield_001000.bin`: Pressure field at step 1000.
-   `position_001000.bin`: Positions of all Lagrangian particles at step 1000.
-   `velocity_001000.bin`: Velocities of all Lagrangian particles at step 1000.

@section run_sim_sec 3. Running a Simulation

Once the run environment is set up as described above, you can execute the solver using `mpirun`.

From within your run directory (`my_simulation/` in our example):
```bash
# Run the solver on 4 processor cores
mpirun -np 4 ./picsolver
```
The simulation will start, printing output to the console (as filtered by the whitelist) and writing detailed logs and binary results to their respective directories.

@section output_control_sec 4. Controlling Simulation Output

You have two primary mechanisms for controlling the verbosity of the solver's console output, which is very useful for debugging and monitoring.

@subsection log_levels_sec 4.1. Log Levels
The solver uses a tiered logging system. You can set a maximum log level (typically in `control.dat`) to control what gets printed. The levels are, in order of increasing verbosity:
-   `LOG_ERROR`: Only critical, simulation-ending errors are shown.
-   `LOG_WARNING`: Errors and potential issues are shown.
-   `LOG_INFO`: (Default) Warnings, errors, and general high-level simulation progress.
-   `LOG_DEBUG`: Highly detailed output from deep within functions, intended for developers.
-   `LOG_PROFILE`: Messages exclusively related to performance profiling.

@subsection whitelist_sec 4.2. Function Whitelist
For highly targeted debugging, you can use the `whitelist.dat` file. If this file exists in the `config/` directory, the solver will *only* print console messages from functions whose names are listed in this file. This allows you to isolate the output from a specific part of the code you are working on without being flooded by messages from other functions.

@section next_steps_sec Next Steps
After running a simulation, the next step is to analyze the data. Refer to the documentation for the `postprocessor` utility and any provided visualization scripts to learn how to process the binary files in the `results/` directory.
