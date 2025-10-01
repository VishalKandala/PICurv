@page 02_Tutorial_Programmatic_Grid Tutorial: Your First Simulation (Flat Channel)

This tutorial will guide you through the complete end-to-end workflow for running a simulation with PICurv. We will use the `flat_channel` example, a fundamental validation case that demonstrates laminar, incompressible flow developing in a straight duct.

This case is a perfect starting point because the grid is **generated programmatically** by the solver based on simple parameters, allowing you to quickly get a simulation running without needing external tools.

@tableofcontents

@section init_sec 1. Initializing the Study Directory

The first step is to create a self-contained "study" directory for our simulation. The `pic.flow` conductor script does this for us using the `init` command.

1.  Navigate to the root of your `PICurv` project directory.
2.  Run the following command:

    ```bash
    ./bin/pic-flow init flat_channel --dest my_first_run
    ```

This command finds the `flat_channel` template in the `examples/` directory and copies all its files into a new folder named `my_first_run` in your current directory.

Your new study directory will have this structure:
```
my_first_run/
├── flat_channel.yml
├── Imp-MG-Standard.yml
└── Standard_Output.yml
```
These YAML (`.yml`) files contain all the settings for our simulation.

@section config_sec 2. Understanding the Configuration Files

Let's take a quick look at the three files that were just created.

- **`flat_channel.yml` (The "Case" file):** This is the most important file. It defines the *physics* and *geometry* of the problem. If you open it, you'll see sections for:
    - `properties`: Defines physical constants like fluid density and viscosity, and the reference scales that determine the Reynolds number.
    - `run_control`: Sets the simulation duration (`total_steps`) and timestep (`dt_physical`).
    - `grid`: Crucially, this is set to `mode: programmatic_c` and specifies the grid resolution (`im`, `jm`, `km`).
    - `boundary_conditions`: Defines what happens at the edges of our domain (inlet, outlet, walls).

- **`Imp-MG-Standard.yml` (The "Solver Profile"):** This file defines the *numerical strategy* the C-solver will use. It controls tolerances, the type of linear solver, and multigrid settings. You can think of it as the "engine" of the simulation, which can be swapped out without changing the physics.

- **`Standard_Output.yml` (The "Monitor Profile"):** This file controls the "instrumentation" of the run. It defines how often to save results (`io.data_output_frequency`) and how verbose the console output should be (`logging.verbosity`).

@section run_sec 3. Running the Simulation

With our study directory prepared, we can now execute the simulation. The `pic.flow run` command orchestrates everything.

1.  Make sure you are still in the root of your `PICurv` project directory.
2.  Run the following command. We will run on 4 processor cores (`-n 4`) and tell the script to execute both the solver (`--solve`) and the post-processor (`--post-process`) stages.

    ```bash
    ./bin/pic-flow run \
        --case my_first_run/flat_channel.yml \
        --solver my_first_run/Imp-MG-Standard.yml \
        --monitor my_first_run/Standard_Output.yml \
        --post config/postprocessors/standard_analysis.yml \
        -n 4 --solve --post-process
    ```
    *Note: For the `--post` flag, we are using a standard analysis recipe provided with the code.*

The `pic.flow` script will now:
1.  Create a unique, timestamped run directory inside a new `runs/` folder.
2.  Read all the `.yml` files.
3.  Convert the user-friendly YAML settings into a machine-readable `.control` file.
4.  Launch the C-solver (`picsolver`) using `mpiexec`.
5.  After the solver finishes, it will automatically launch the C post-processor (`postprocessor`).

You will see progress updates printed to your console.

@section results_sec 4. Examining the Output

Once the run is complete, you will have a new directory structure inside the `runs/` folder, something like this:

```
runs/
└── flat_channel_20240401-153000/  (Your case name + timestamp)
    ├── config/                   (A copy of all input files)
    ├── logs/                     (Detailed log files from the C-solver)
    ├── results/                  (Raw, non-dimensional binary output from the solver)
    └── viz/                      (Final, dimensional visualization files from the post-processor)
```

The most important directory for visualization is **`viz/`**. It contains a series of `.vts` (VTK Structured Grid) files, one for each output timestep (e.g., `Field_000100.vts`, `Field_000200.vts`).

@section viz_sec 5. Visualizing the Results

We will use **ParaView**, a powerful, open-source scientific visualization tool, to view the output.

1.  **Install ParaView:** Download and install it from [paraview.org](https://www.paraview.org/download/).
2.  **Open the File Series:**
    -   Launch ParaView.
    -   Go to `File -> Open`.
    -   Navigate into your `runs/flat_channel_.../viz/` directory.
    -   ParaView will automatically group the `Field_....vts` files into a single time series. Select the group (it will look like `Field..vts`) and click **OK**.
    -   In the "Properties" panel on the left, click the blue **Apply** button. Your channel geometry will appear in the 3D view.

3.  **Create a Slice to See the Profile:**
    -   With the `Field..vts` object selected in the "Pipeline Browser", click the **Slice** filter icon in the toolbar (it looks like a plane cutting a cube).
    -   In the Properties panel, change the "Slice Type" to `Plane`. Make sure the "Normal" is pointing in the direction of the flow (e.g., if flow is in Z, set Normal to `0, 0, 1`).
    -   Click **Apply**. You now have a 2D slice through the middle of your channel.

4.  **Color by Velocity:**
    -   In the toolbar at the top, find the dropdown menu that says "Solid Color".
    -   Change it to **`U_nodal`**.
    -   Now find the component dropdown next to it (usually says "Magnitude") and select the component corresponding to the flow direction (e.g., `Z`).
    -   Click the "Rescale to data range" button (a rainbow-colored arrow) to adjust the color map.

You should now see the classic parabolic "Poiseuille" velocity profile, which is dark blue at the walls (zero velocity) and red in the center (maximum velocity). You can use the play buttons at the top to see how this profile develops over time.

@image html paraview_flat_channel.png "Example of a flat channel velocity profile visualized in ParaView." width=700px

@section next_steps_sec 6. Next Steps

You have successfully run your first simulation! You can now experiment by changing parameters in `my_first_run/flat_channel.yml` (like the Reynolds number or grid resolution) and re-running the `pic.flow run` command.

When you are ready, proceed to the next tutorial to learn how to work with more complex, file-based geometries: **@subpage 03_Tutorial_File-Based_Grid**.
