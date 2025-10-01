@page 07_Case_Reference Configuration Reference: The `case.yml` File

For a complete, heavily commented reference file showing every possible option, please see the master template:

@verbinclude master_template/master_case.yml

The `case.yml` file is the heart of your simulation setup. It defines the fundamental physical properties of the problem, the geometry of the domain, the boundary conditions, and the duration of the run.

This document serves as a comprehensive reference for all available sections and parameters within the `case.yml` file.

@tableofcontents

@section properties_sec 1. The `properties` Section

This section defines the core physical scales and fluid properties. These values are used by the `pic-flow` script to calculate the non-dimensional numbers (like the Reynolds number) that are passed to the C-solver. For more details, see the @subpage 14_Theory_and_Numerics "Theory & Numerics" page on non-dimensionalization.

```yaml
properties:
  scaling:
    length_ref: 0.1
    velocity_ref: 1.5
  fluid:
    density: 1000.0
    viscosity: 0.001
  initial_conditions:
    mode: "Constant"
    u_physical: 1.5
    v_physical: 0.0
    w_physical: 0.0
```

@subsection scaling_ssec 1.1. `scaling`

| Parameter | Type | Description | C-Solver Flag |
| :--- | :--- | :--- | :--- |
| `length_ref` | Real | A characteristic length of the problem (e.g., channel height) in meters. | `-scaling_L_ref` |
| `velocity_ref` | Real | A characteristic velocity (e.g., inlet velocity) in m/s. | `-scaling_U_ref` |

@subsection fluid_ssec 1.2. `fluid`

| Parameter | Type | Description |
| :--- | :--- | :--- |
| `density` | Real | Reference fluid density (`ρ_ref`) in kg/m³. Used to calculate Reynolds number. |
| `viscosity` | Real | Fluid dynamic viscosity (`μ`) in Pa·s. Used to calculate Reynolds number. |

@subsection ic_ssec 1.3. `initial_conditions`

Defines the state of the flow field at t=0.

| Parameter | Type | Description | C-Solver Flag |
| :--- | :--- | :--- | :--- |
| `mode` | String | Sets the initial velocity profile. Valid options are: `"Zero"`, `"Constant"`, `"Poiseuille"`. | `-finit` |
| `u_physical`| Real | The x-component of velocity in m/s, used if `mode` is "Constant". | `-ucont_x` |
| `v_physical`| Real | The y-component of velocity in m/s, used if `mode` is "Constant". | `-ucont_y` |
| `w_physical`| Real | The z-component of velocity in m/s, used if `mode` is "Constant". | `-ucont_z` |


@section run_control_sec 2. The `run_control` Section

This section governs the time-stepping and duration of the simulation.

```yaml
run_control:
  dt_physical: 0.0001
  start_step: 0
  total_steps: 2000
```

| Parameter | Type | Description | C-Solver Flag |
| :--- | :--- | :--- | :--- |
| `dt_physical` | Real | The physical time step size in seconds. `pic-flow` converts this to a non-dimensional `dt*`. | `-dt` |
| `start_step` | Integer | The time step number to start or restart from. A value > 0 indicates a restart. | `-start_step` |
| `total_steps` | Integer | The total number of time steps to execute in this run. | `-totalsteps` |

@section grid_sec 3. The `grid` Section

This section defines the computational mesh. There are two primary modes: `programmatic_c` and `file`.

@subsection programmatic_ssec 3.1. Programmatic Grid (`mode: programmatic_c`)

This mode instructs the C-solver to generate a stretched Cartesian grid internally.

```yaml
grid:
  mode: programmatic_c
  programmatic_settings:
    im: [65]
    jm: [33]
    km: [33]
    xMins: [0.0]
    xMaxs: [1.0]
    # ... and so on for y and z ...
    rxs: [1.0] # Stretching ratio
    da_processors_x: [4] # Optional MPI layout
```

| Parameter | Type | Description | C-Solver Flag |
| :--- | :--- | :--- | :--- |
| `im`, `jm`, `km` | Int Array | Number of **cells** in each direction. Use one entry per block for multi-block. | `-im`, `-jm`, `-km` |
| `xMins`, `xMaxs` | Real Array | Physical domain boundaries for each block in the x-direction. | `-xMins`, `-xMaxs`|
| `yMins`, `yMaxs` | Real Array | Physical domain boundaries for each block in the y-direction. | `-yMins`, `-yMaxs`|
| `zMins`, `zMaxs` | Real Array | Physical domain boundaries for each block in the z-direction. | `-zMins`, `-zMaxs`|
| `rxs`, `rys`, `rzs`| Real Array | Geometric stretching ratio in each direction. `1.0` is uniform. `>1.0` clusters points near the minimum boundary. | `-rxs`, `-rys`, `-rzs`|
| `da_processors_x/y/z`| Int Array | (Optional) Explicitly sets the number of MPI processes to use in each direction for domain decomposition. The product must equal the total number of processes (`-n`). If omitted, PETSc decides automatically. | `-da_processors_x/y/z` |

@subsection file_ssec 3.2. File-Based Grid (`mode: file`)

This mode instructs the solver to use an external grid file, necessary for complex or curvilinear geometries.

```yaml
grid:
  mode: file
  source_file: my_grids/bent_channel.picgrid
```

| Parameter | Type | Description | C-Solver Flag |
| :--- | :--- | :--- | :--- |
| `source_file` | String | Path to the grid coordinate file (e.g., `.picgrid`). The path is relative to the `case.yml` file's location. | `-grid_file` |

@section models_sec 4. The `models` Section

This section is for selecting high-level physics models and domain properties.

```yaml
models:
  domain:
    blocks: 1
    i_periodic: false
  physics:
    dimensionality: "3D"
    turbulence:
      les: false
    particles:
      count: 0
      init_mode: "Surface"
```

@subsection domain_ssec 4.1. `domain`

| Parameter | Type | Description | C-Solver Flag |
| :--- | :--- | :--- | :--- |
| `blocks` | Integer | The number of computational blocks for multi-block simulations. | `-nblk` |
| `i_periodic`| Boolean | If `true`, sets the i-direction (X) as periodic. | `-i_periodic` |
| `j_periodic`| Boolean | If `true`, sets the j-direction (Y) as periodic. | `-j_periodic` |
| `k_periodic`| Boolean | If `true`, sets the k-direction (Z) as periodic. | `-k_periodic` |

@subsection physics_ssec 4.2. `physics`

| Parameter | Type | Description | C-Solver Flag |
| :--- | :--- | :--- | :--- |
| `dimensionality`| String | Set to `"2D"` to run a 2D simulation. Defaults to `"3D"`. | `-TwoD` |
| `turbulence.les`| Boolean | If `true`, enables the Large Eddy Simulation (LES) model. | `-les` |
| `particles.count`| Integer | The total number of Lagrangian particles to simulate. Set to `0` to disable particles. | `-numParticles`|
| `particles.init_mode`| String | How particles are initialized at t=0. `"Surface"` places them on the primary inlet face. `"Volume"` places them randomly throughout the domain. | `-pinit` |
| `particles.restart_mode`| String | For restart runs (`start_step > 0`). `"load"` loads particle positions from files. `"init"` generates a fresh set of particles in the restarted flow field. | `-particle_restart_mode` |

@section bc_sec 5. The `boundary_conditions` Section

This section defines the physical behavior at each of the six faces of the computational domain(s). For multi-block cases, this is a list of lists.

**Syntax:**
```yaml
boundary_conditions:
  # Block 0 Definitions
  - - face: "-Xi"
      type: "INLET"
      handler: "constant_velocity"
      params:
        vx: 1.5
        vy: 0.0
        vz: 0.0
    - face: "+Xi"
      type: "OUTLET"
      handler: "conservation"
    - face: "-Eta"
      type: "WALL"
      handler: "noslip"
    # ... other faces ...
```

**Main Keys:**
| Key | Description |
| :--- | :--- |
| `face` | (Required) Identifies the boundary face. Valid options: `"-Xi"`, `"+Xi"`, `"-Eta"`, `"+Eta"`, `"-Zeta"`, `"+Zeta"`. |
| `type` | (Required) The general physical category of the boundary. See table below. |
| `handler`| (Required) The specific C-function that implements the boundary condition. See table below. |
| `params` | (Optional) A dictionary of key-value pairs providing parameters to the handler. |

**Common `type` and `handler` Combinations:**
| `type` | `handler` | Description & Required `params` |
| :--- | :--- | :--- |
| `INLET` | `constant_velocity` | Specifies a uniform inlet velocity. Requires `vx`, `vy`, `vz` in physical units (m/s). |
| `INLET` | `parabolic` | Specifies a parabolic inlet profile. Requires `u_max` (the centerline velocity). |
| `OUTLET`| `conservation` | A simple zero-gradient outflow condition that enforces global mass conservation. |
| `WALL` | `noslip` | A standard no-slip wall (zero velocity relative to the wall). |
| `SYMMETRY`| `symmetry_plane` | A slip-wall condition for symmetry planes. |
| `NOGRAD`| `allcopy` | A generic zero-gradient condition. |

@section next_steps_sec 6. Next Steps

With a full understanding of the `case.yml` file, you are now equipped to define the physics of your own simulations. Next, you will need to learn how to control the numerical schemes used to solve the problem.

Proceed to the **@subpage 08_Solver_Reference** to learn about all the parameters available in the `solver.yml` file.
