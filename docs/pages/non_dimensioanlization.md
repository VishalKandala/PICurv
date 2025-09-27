# Non-Dimensionalization Strategy in PIC-Flow

## 1. Overview

This document describes the non-dimensionalization strategy employed by the `pic-flow` simulation platform. The entire workflow is built on the **"Dimensional In, Non-Dimensional Core, A-Posteriori Scaling"** model. This approach combines user-friendliness, numerical stability, and scientific rigor.

-   **Dimensional In:** The user defines a simulation case using simple, intuitive configuration files with familiar physical units (e.g., meters, seconds, m/s, Pa·s).
-   **Non-Dimensional Core:** The `pic-flow` Python conductor script reads these physical values and converts them into a consistent set of non-dimensional numbers before launching the solver. The high-performance C solver (`picsolver`) **only** ever operates on these non-dimensional quantities.
-   **A-Posteriori Scaling:** The C solver writes its raw output (grid and field data) in non-dimensional format. The C post-processor (`postprocessor`), also orchestrated by `pic-flow`, reads this non-dimensional data, converts it back to physical units in-memory, and then performs analysis and generates dimensional visualization files.

This strategy ensures that the user experience is intuitive, the core solver is numerically robust, and the results are scientifically general and easy to analyze.

## 2. Characteristic Scales

The conversion between dimensional and non-dimensional quantities is based on a set of user-provided characteristic scales. These are defined in the `[characteristic_scales]` section of the `case.yml` file.

| Parameter | Symbol | Description | Units |
| :--- | :--- | :--- | :--- |
| `length_ref` | `L_ref` | A characteristic length of the problem. | `m` |
| `velocity_ref` | `U_ref` | A characteristic velocity of the problem. | `m/s`|
| `density` | `ρ_ref` | The reference density of the fluid. | `kg/m³`|

From these primary scales, several important secondary scales are derived:

-   **Characteristic Time (`T_ref`):** `T_ref = L_ref / U_ref` [s]
-   **Characteristic Pressure (`P_ref`):** `P_ref = ρ_ref * U_ref²` [Pa]

## 3. The Non-Dimensionalization Process

The `pic-flow` script automatically performs the following conversions. A quantity with an asterisk (`*`) denotes its non-dimensional form.

### 3.1. Primary Simulation Parameters

These are the key dimensionless numbers that govern the physics of the flow.

| Parameter | Symbol | Formula | C Solver Flag |
| :--- | :--- | :--- | :--- |
| **Reynolds Number** | `Re` | `Re = (ρ_ref * U_ref * L_ref) / μ` | `-ren` |
| **Time Step** | `dt*` | `dt* = dt_physical / T_ref` | `-dt` |

The C solver (`picsolver`) receives these values directly in its `control.dat` file and uses them in its governing equations. The `Viscous` kernel, for example, uses `1/Re` as the effective kinematic viscosity.

### 3.2. Field and Coordinate Scaling

All physical fields are converted to their non-dimensional counterparts before the solver begins and are converted back during post-processing.

| Field Description | `UserCtx` / `DMSwarm` Member | Dimensional Units | Scaling Equation | Scaling Factor |
| :--- | :--- | :--- | :--- | :--- |
| **Grid Coordinates** | `DMDA Coordinates` | `L` | `x* = x / L_ref` | `L_ref` |
| **Cartesian Velocity** | `user->Ucat`, particle `velocity` | `L/T` | `u* = u / U_ref` | `U_ref` |
| **Contravariant Volume Flux** | `user->Ucont` | `L³/T` | `Ucont* = Ucont / (U_ref * L_ref²)` | `U_ref * L_ref²` |
| **Pressure** | `user->P` | `M L⁻¹ T⁻²` | `p* = p / P_ref` | `P_ref` |
| **Particle Position** | `DMSwarm` "position" field | `L` | `x_p* = x_p / L_ref` | `L_ref` |

### 3.3. Grid Metrics (Internal Conversion)

The user provides a dimensional grid file. The `pic-flow` script reads this file, scales the coordinates by `1/L_ref`, and writes a new, non-dimensional grid file. The C solver then reads this non-dimensional grid and calculates all metrics. As a result, the metric fields within the solver are inherently non-dimensional. Their scaling factors are listed here for completeness.

| Metric Description | `UserCtx` Member | Dimensional Units | Scaling Factor |
| :--- | :--- | :--- | :--- |
| **Face Area Vectors**| `Csi`, `Eta`, `Zet`, `ICsi`... | `L²` | `L_ref²` |
| **Inverse Jacobian** | `Aj`, `IAj`, `JAj`, `KAj`| `1/L³` | `1 / L_ref³` |
| **Cell/Face Centers** | `Cent`, `Centx`, etc. | `L` | `L_ref` |

## 4. Implementation in the Code

### `pic-flow` Conductor Script

-   **Responsibility:** Performs all calculations listed in sections 3.1 and 3.2.
-   **Action:**
    1.  Reads the dimensional `case.yml`.
    2.  Reads the dimensional grid file, scales it, and writes a temporary non-dimensional grid file.
    3.  Calculates `Re`, `dt*`, etc.
    4.  Generates the `control.dat` file, populating it with all non-dimensional values.
    5.  It also passes the original reference scales (`-scaling_L_ref`, etc.) to the C executables for use in post-processing.

### C Solver (`picsolver`)

-   **Responsibility:** Solve the non-dimensional governing equations.
-   **Action:**
    1.  Reads `control.dat` to get non-dimensional parameters like `-ren` and `-dt`.
    2.  Reads the non-dimensional grid file.
    3.  Computes all grid metrics, which are automatically non-dimensional.
    4.  Runs the simulation entirely in non-dimensional space.
    5.  Writes all output fields (`Ucat`, `P`, `position`, etc.) in their raw, **non-dimensional** format.

### C Post-Processor (`postprocessor`)

-   **Responsibility:** Analyze data and produce human-interpretable, dimensional output.
-   **Action:**
    1.  In `main()`, it calls `CreateSimulationContext`, which reads the `-scaling_*` flags from the configuration provided by `pic-flow`.
    2.  In its main processing loop, for each timestep, it first calls `ReadSimulationFields()` to load the **non-dimensional** data from the solver's output.
    3.  It immediately calls a dedicated function, `DimensionalizeAllLoadedFields()`. This function uses the stored scaling factors to convert all loaded fields (`Ucat`, `P`, coordinates, etc.) into physical, dimensional units **in-memory**.
    4.  All subsequent analysis kernels (`ComputeQCriterion`, etc.) and the final VTK writer (`WriteEulerianFile`) operate on this dimensional data.

This architecture ensures a clean separation of concerns, allowing the user to work in the physical domain while the solver benefits from the stability and generality of a non-dimensional formulation.