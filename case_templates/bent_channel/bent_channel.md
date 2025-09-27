# Case Template: Laminar Flow in a Bent Channel

## 1. Description

This template simulates **laminar, incompressible flow** through a 90-degree bent channel with a square cross-section. It is designed to demonstrate the platform's ability to handle complex, curvilinear geometries.

The key feature of this template is that the grid is **read from a file (`grid.picgrid`)**, not generated programmatically. This is the standard workflow for any non-trivial geometry.

**Default Physics:**
- **Flow Type:** Laminar
- **Reynolds Number (Re):** 200
- **Geometry:** 90-Degree Bent Square Duct

## 2. Files in this Template

- **`bent_channel.yml`**: The main configuration file. Note that `grid.mode` is set to **`file`** and `source_file` points to the included grid.
- **`Imp-MG-Standard.yml`**: A general-purpose implicit **solver profile**. This is the same file as in the `flat_channel` template, demonstrating the modularity of solver profiles.
- **`Standard_Output.yml`**: A standard **monitoring profile**. This is also the same as in the `flat_channel` template, showing that monitor profiles are reusable across different cases.

### Included Grid Artifacts
- **`bent_channel_coarse.picgrid`**: The raw grid coordinate file used by the solver.
- **`bent_channel_coarse.vts`**: A VTK file of the grid. **You should open this file in ParaView first** to visualize the mesh and understand its structure before running a simulation.
- **`bent_channel_coarse.info`**: A text file containing quality metrics and statistics about the grid.

## 3. How to Use this Template

### Step 1: Initialize a New Study

From your project root, run the `init` command:
```bash
# Example: Create a new study directory for the bent channel
./bin/pic.flow init bent_channel --dest my_bent_channel_study
```
This will create the folder `my_bent_channel_study/` and copy all files, including all three grid artifacts.

### Step 2: Customize the Simulation

Navigate into `my_bent_channel_study/` and open `case.yml`. You can modify physical properties (like `density` to change the Reynolds number) and `run_control` parameters.

**Important:** The grid is controlled by the included file. To change the geometry, you must:
1.  Generate a new grid using the `grid.gen` utility or another tool.
2.  Place the new `.picgrid`, `.vts`, and `.info` files into your study directory.
3.  Update the `source_file` path in `case.yml` to point to your new grid file.

### Step 3: Run the Simulation

From the project root, launch the simulation:
```bash
# Run the simulation on 4 processor cores
./bin/pic.flow run \
    --case my_bent_channel_study/case.yml \
    --solver my_bent_channel_study/Imp-MG-Standard.yml \
    --monitor my_bent_channel_study/Standard_Output.yml \
    -n 4 --solve
```
Notice you can easily swap `--solver` or `--monitor` to point to a different profile from the central libraries (e.g., `solver_profiles/Robust-Solver.yml`) to test different numerical schemes on this complex geometry.

## 4. Expected Results & Visualization

The output will be a series of `.vts` files in the `runs/` directory.
- Open the files in **ParaView**.
- Use a **Stream Tracer** or **Glyph** filter on the `U_nodal` vector field to visualize the flow path as it navigates the bend. You will see higher velocities on the inner wall of the bend.