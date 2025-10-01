@page 03_Tutorial_File-Based_Grid Tutorial: Using a File-Based Grid (Bent Channel)

In the previous tutorial, we ran a simulation where the grid was generated automatically. For complex, real-world geometries, you will typically use a dedicated meshing tool to create a grid and provide it to the solver as a file.

This tutorial will guide you through running the `bent_channel` example, which demonstrates how to use a pre-computed curvilinear grid file.

@tableofcontents

@section overview_sec 1. Case Overview: Laminar Flow in a Bend

This template simulates incompressible laminar flow through a 90-degree bent channel. The key feature is that the grid is **read from a file**, demonstrating the workflow for non-trivial geometries.

**Physics:**
- **Flow Type:** Laminar
- **Reynolds Number (Re):** 200
- **Geometry:** 90-Degree Bent Square Duct

@image html bent_channel_grid.png "Visualization of the bent channel grid file in ParaView." width=600px

@section init_sec 2. Initializing the Study

Just like before, we use the `pic.flow init` command to create our study directory.

1.  Navigate to the root of your `PICurv` project directory.
2.  Run the following command:

    ```bash
    ./bin/pic-flow init bent_channel --dest my_bent_channel_run
    ```

This copies the `bent_channel` template into a new `my_bent_channel_run` directory. If you look inside, you will see a new set of files:

```
my_bent_channel_run/
├── bent_channel.yml
├── Imp-MG-Standard.yml
├── Standard_Output.yml
├── standard_analysis.yml
├── bent_channel_coarse.picgrid  <-- The grid coordinate file
├── bent_channel_coarse.vts      <-- A VTK file for visualizing the grid
└── bent_channel_coarse.info     <-- Grid quality metrics (for reference)
```

@section config_sec 3. The File-Based Grid Configuration

Open the `my_bent_channel_run/bent_channel.yml` file and look at the `grid` section. It's different from the previous tutorial:

```yaml
grid:
  mode: file
  source_file: bent_channel_coarse.picgrid
```

- **`mode: file`**: This tells `pic.flow` that the C-solver should not generate its own grid.
- **`source_file: ...`**: This provides the path to the coordinate file (`.picgrid`) that the solver should use.

Notice that the `Imp-MG-Standard.yml` and `Standard_Output.yml` files are identical to the ones from the `flat_channel` example. This highlights a core design principle of PICurv: **physical cases, numerical solvers, and monitoring are modular and interchangeable.**

@section run_sec 4. Running the Simulation

The run command is almost identical to the previous tutorial. We just point the `--case` flag to our new `bent_channel.yml` file.

1.  From the project root directory, execute the `run` command:

    ```bash
    ./bin/pic-flow run \
        --case my_bent_channel_run/bent_channel.yml \
        --solver my_bent_channel_run/Imp-MG-Standard.yml \
        --monitor my_bent_channel_run/Standard_Output.yml \
        --post my_bent_channel_run/standard_analysis.yml \
        -n 4 --solve --post-process
    ```
2.  The solver will start, read the grid file, and run the simulation.

@section viz_sec 5. Visualizing the Results

The output will again be in a timestamped directory inside `runs/`, with the final visualization files in the `viz/` subdirectory.

1.  **Open the File Series in ParaView:**
    -   Go to `File -> Open` and select the `Field..vts` series from the `viz/` directory.
    -   Click **Apply**.

2.  **Visualize the Flow Path:**
    -   The best way to see the flow in a bend is with a **Stream Tracer**.
    -   With `Field..vts` selected, click the Stream Tracer filter icon in the toolbar.
    -   The "Seed" for the tracer is a point or line where streamlines will start. By default, it's a line. Drag the line so it is positioned at the inlet of your channel.
    -   In the Properties panel, increase the number of "Resolution" points on the line to generate more streamlines.
    -   Click **Apply**.
    -   Color the resulting streamlines by the `U_nodal` vector field magnitude.

You will see the flow navigating the bend, with higher velocities on the inner wall and lower velocities on the outer wall, which is the expected physical result.

@image html paraview_bent_channel.png "Example of streamlines in the bent channel, colored by velocity magnitude." width=700px

@section next_steps_sec 6. Next Steps

You have now successfully run a simulation using both a programmatic grid and a file-based grid. You understand the core user workflow of PICurv.Before moving on to create your own cases, it's essential to learn the fundamentals of visualizing and analyzing your results.

Please proceed to the final introductory tutorial: **@subpage 04_Visualization_Tutorial**.

