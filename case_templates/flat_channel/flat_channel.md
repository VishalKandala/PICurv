# Case Template: Laminar Flow in a Flat Channel

## 1. Description

This template simulates **laminar, incompressible flow** through a straight channel with a square cross-section. The fluid enters with a uniform velocity, and a parabolic "Poiseuille" velocity profile develops as the flow moves down the channel.

This is a fundamental validation case and an excellent starting point for learning the PIC-Flow workflow. The grid for this simulation is **generated programmatically** by the C-solver based on the settings in `case.yml`.

**Default Physics:**
- **Flow Type:** Laminar
- **Reynolds Number (Re):** 200
- **Geometry:** Straight Square Duct

## 2. Files in this Template

- **`flat_channel.yml`**: The main configuration file. Defines the physics, geometry, and run duration. **This is the primary file you will edit.**
- **`Imp-MG-Standard.yml`**: A robust, pre-configured **solver profile** using an implicit method with a multigrid preconditioner. It's a good general-purpose solver.
- **`Standard_Output.yml`**: A standard **monitoring profile** that shows INFO-level progress updates and saves results every 100 timesteps.

## 3. How to Use this Template

### Step 1: Initialize a New Study

Navigate to your main project root directory (the one containing `bin/`) and run the `init` command. It's best practice to give your study a descriptive name using the `--dest` flag.

```bash
# Example: Create a new study directory named "my_Re400_channel_study"
./bin/pic.flow init flat_channel --dest my_Re400_channel_study
```
This will create a new folder `my_Re400_channel_study/` containing copies of all the template files.

### Step 2: Customize the Simulation

Navigate into your new directory (`cd my_Re400_channel_study`) and open `case.yml` in a text editor. The most common parameters you will want to change are:
- **Reynolds number:** Modify `density`, `viscosity`, `length_ref`, or `velocity_ref` under the `properties` section.
- **Run duration:** Modify `total_steps` or `dt_physical` under the `run_control` section.
- **Grid resolution:** Modify `im`, `jm`, or `km` under `grid.programmatic_settings`.

### Step 3: Run the Simulation

Navigate back to the project root (`cd ..`) and launch the simulation using the `run` command, pointing to the files inside your new study directory.

```bash
# Run the simulation on 4 processor cores
./bin/pic.flow run \
    --case my_Re400_channel_study/case.yml \
    --solver my_Re400_channel_study/Imp-MG-Standard.yml \
    --monitor my_Re400_channel_study/Standard_Output.yml \
    -n 4 --solve
```

## 4. Building a New Template from this Example

If you want to create a new template for a different programmatic grid case (e.g., a stretched channel), the easiest way is to:
1.  **Initialize a study from this template:** `pic.flow init flat_channel --dest new_stretched_channel`
2.  **Customize the files:** Go into `new_stretched_channel/` and edit `case.yml` (e.g., change the `rzs` stretching ratio). Rename `Imp-MG-Standard.yml` or `Standard_Output.yml` if you make significant changes to them.
3.  **Add a `README.md`:** Write a new README explaining your changes.
4.  **Copy to `case_templates/`:** Copy the entire `new_stretched_channel` directory into the main `case_templates/` directory to make it a reusable template.

## 5. Expected Results & Visualization

After the run completes, a new directory will be created under `runs/`. Inside its `results/` subdirectory, you will find a series of `.vts` files.
- Open these files in a visualization tool like **ParaView**.
- Create a "Slice" through the center of the channel and color it by `U_nodal` magnitude to see the classic parabolic velocity profile develop.