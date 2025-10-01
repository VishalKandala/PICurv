@page 11_User_How_To_Guides User How-To Guides

This page provides a collection of quick, goal-oriented guides for common simulation and configuration tasks. Use these recipes to solve specific problems without needing to read the full reference manuals.

@tableofcontents

@section setup_guides_sec 1. Simulation Setup

@subsection ht_reynolds_ssec 1.1. How do I change the Reynolds number?

The Reynolds number (`Re`) is calculated automatically by `pic-flow` from the physical properties in your `case.yml` file using the formula `Re = (ρ * U * L) / μ`.

To change `Re`, you can modify any of the four defining parameters in the `properties` section. For example, to double the Reynolds number, you could either double the `density` or halve the `viscosity`.

**`case.yml`:**
```yaml
properties:
  scaling:
    length_ref: 0.1     # L
    velocity_ref: 1.5   # U
  fluid:
    density: 1000.0     # ρ
    viscosity: 0.001    # μ
```

@subsection ht_2d_ssec 1.2. How do I set up a 2D simulation?

Set the `dimensionality` parameter to `"2D"` in your `case.yml`. The C-solver will then ignore computations in the z-direction. For a programmatic grid, you should also set the resolution in the z-direction (`km`) to a small number (e.g., 3) to minimize computational overhead.

**`case.yml`:**
```yaml
models:
  physics:
    dimensionality: "2D"

grid:
  mode: programmatic_c
  programmatic_settings:
    im: [129]
    jm: [65]
    km: [3] # Small number for 2D cases
```

@subsection ht_resolution_ssec 1.3. How do I change the grid resolution?

This depends on your grid `mode` in `case.yml`.

-   **For `programmatic_c` grids:** Change the number of cells in the `im`, `jm`, and `km` lists.
    ```yaml
    grid:
      mode: programmatic_c
      programmatic_settings:
        im: [257] # Doubled resolution
        jm: [129]
        km: [129]
    ```
-   **For `file` grids:** You must generate a new grid file using an external tool and update the `source_file` path.
    ```yaml
    grid:
      mode: file
      source_file: my_grids/new_fine_mesh.picgrid
    ```

@subsection ht_time_ssec 1.4. How do I change the simulation time and output frequency?

Simulation time is controlled in `case.yml`, while output frequency is in `monitor.yml`.

-   **To run for a longer time:** Increase `total_steps` or adjust `dt_physical` in `case.yml`.
    ```yaml
    # In case.yml
    run_control:
      dt_physical: 0.0001
      total_steps: 5000 # Run for 5000 steps instead of 2000
    ```
-   **To save results more frequently:** Decrease `data_output_frequency` in `monitor.yml`.
    ```yaml
    # In monitor.yml
    io:
      data_output_frequency: 50 # Save every 50 steps instead of 100
    ```

@section bc_guides_sec 2. Boundary Conditions

@subsection ht_bc_ssec 2.1. How do I set up a simple boundary condition?

Define it in the `boundary_conditions` list in your `case.yml`. Each entry needs a `face`, `type`, and `handler`. Some handlers require `params`.

**Example: A constant velocity inlet and a no-slip wall.**
```yaml
# In case.yml
boundary_conditions:
  - face: "-Zeta" # K-min face
    type: "INLET"
    handler: "constant_velocity"
    params:
      vx: 0.0
      vy: 0.0
      vz: 1.5
  - face: "-Eta"  # J-min face
    type: "WALL"
    handler: "noslip"
```
Refer to the @subpage 06_Case_Reference "Case Reference" for a full list of handlers.

@subsection ht_periodic_ssec 2.2. How do I set up a periodic domain?

Set the corresponding periodic flag to `true` under `models.domain` in your `case.yml`. For example, for a channel that is periodic in the x-direction:

**`case.yml`:**
```yaml
models:
  domain:
    i_periodic: true
```
You must also ensure that the boundary conditions for the periodic faces (e.g., `-Xi` and `+Xi`) are set to the `INTERFACE` type with the `periodic` handler.

@section run_guides_sec 3. Running & Monitoring

@subsection ht_parallel_ssec 3.1. How do I run in parallel and specify the processor layout?

Use the `-n` flag with `pic-flow run`. To optimize performance, you can optionally suggest a processor decomposition to PETSc in your `case.yml`.

**Run command:**
```bash
# Run on 16 processor cores
./bin/pic-flow run -n 16 --solve ...
```

**Optional `case.yml` for a 4x2x2 layout:**
```yaml
grid:
  mode: programmatic_c
  programmatic_settings:
    da_processors_x: [4]
    da_processors_y: [2]
    da_processors_z: [2]
    # ... other grid settings ...
```

@subsection ht_restart_ssec 3.2. How do I restart a simulation?

Set the `start_step` in `case.yml` to the time step you want to restart *from*. The solver will automatically find the corresponding output files in the `results` directory.

**`case.yml`:**
```yaml
run_control:
  dt_physical: 0.0001
  start_step: 1000 # Start this run from the data saved at step 1000
  total_steps: 2000
```
This configuration would run the simulation from step 1000 to step 3000.

@subsection ht_debug_ssec 3.3. How do I get detailed debug output from a specific part of the code?

Use the `logging` section of your `monitor.yml`. Set the `verbosity` to `DEBUG` and specify the C function(s) you want to inspect in the `enabled_functions` list.

**`monitor.yml` for debugging the projection step:**
```yaml
logging:
  verbosity: "DEBUG"
  enabled_functions:
    - Projection # Only show debug messages from the Projection() C function
    - UpdatePressure
```

@section pp_guides_sec 4. Post-Processing

@subsection ht_pp_ssec 4.1. How do I run post-processing on a finished simulation?

Use the `--post-process` and `--run-dir` flags with `pic.flow run`. You also need to provide a post-processing recipe with `--post`.

```bash
./bin/pic-flow run --post-process \
    --run-dir runs/flat_channel_20240401-153000 \
    --post config/postprocessors/standard_analysis.yml
```

@subsection ht_newfield_ssec 4.2. How do I calculate a new field like Q-Criterion?

Add the task to the `eulerian_pipeline` list in your `post.yml` recipe. Then, add the name of the output field (`Qcrit` in this case) to the `eulerian_fields` list to ensure it gets saved to the VTK file.

**`my_analysis.yml`:**
```yaml
eulerian_pipeline:
  - task: CellToNodeAverage
    input_field: Ucat
    output_field: Ucat_nodal
  - task: ComputeQCriterion # Add this task

io:
  # ...
  eulerian_fields:
    - Ucat_nodal
    - Qcrit # Add the output field here
```

@section next_steps_sec 5. Next Steps

This page serves as a quick reference for common tasks. For a broader overview of all features available through the configuration files, see the summary page.

Proceed to the **@subpage 12_Capabilities_Summary**.

