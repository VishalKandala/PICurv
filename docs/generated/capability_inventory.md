<!-- GENERATED FILE - do not edit by hand.
     Regenerate with: make docs-inventory
     Source of truth: the Python validation layer named per family below. -->

### Boundary Condition Handlers

Selector: `case.yml -> boundary_conditions[].handler`

| Value | Applies to | Required parameters | Optional parameters |
|---|---|---|---|
| `conservation` | `OUTLET` | none | none |
| `constant_flux` | `PERIODIC` | `target_flux` | `apply_trim`, `enforce_seam_flux` |
| `constant_velocity` | `INLET` | `vx`, `vy`, `vz` | none |
| `geometric` | `PERIODIC` | none | none |
| `initial_flux` | `PERIODIC` | none | `apply_trim`, `enforce_seam_flux` |
| `noslip` | `WALL` | none | none |
| `parabolic` | `INLET` | `v_max` | none |
| `prescribed_flow` | `INLET` | `source` | none |

### Boundary Condition Types

Selector: `case.yml -> boundary_conditions[].type`

| Value | Maps to |
|---|---|
| `inlet` | `INLET` |
| `outlet` | `OUTLET` |
| `periodic` | `PERIODIC` |
| `symmetry` | `SYMMETRY` |
| `wall` | `WALL` |

### Momentum Solvers

Selector: `solver.yml -> strategy.momentum_solver`

| Value | Maps to |
|---|---|
| `Dual Time Picard Jameson RK` | `DUALTIME_PICARD_JAMESON_RK` |
| `Dual Time Picard RK4` | `DUALTIME_PICARD_JAMESON_RK` |
| `Explicit RK4` | `EXPLICIT_RK` |
| `Newton Krylov` | `newton_krylov` |

### LES Subgrid Models

Selector: `case.yml -> models.physics.turbulence.les.model`

| Value | Maps to |
|---|---|
| `constant` | `1` |
| `constant_smagorinsky` | `1` |
| `disabled` | `0` |
| `dynamic` | `2` |
| `dynamic_smagorinsky` | `2` |
| `no_les` | `0` |
| `none` | `0` |
| `off` | `0` |
| `smagorinsky` | `1` |

### Eulerian Field Source

Selector: `solver.yml -> operation_mode.eulerian_source`

| Value | Maps to |
|---|---|
| `analytical` | `analytical` |
| `load` | `load` |
| `solve` | `solve` |

### Field Initialization Modes

Selector: `case.yml -> initial_conditions field mode`

| Value | Maps to |
|---|---|
| `Constant` | `1` |
| `Poiseuille` | `2` |
| `Zero` | `0` |

### Particle Initialization Modes

Selector: `case.yml -> models.physics.particles.init_mode`

| Value | Maps to |
|---|---|
| `PointSource` | `2` |
| `Surface` | `0` |
| `SurfaceEdges` | `3` |
| `Volume` | `1` |

### Grid-to-Particle Interpolation

Selector: `solver.yml -> interpolation.method`

| Value | Maps to |
|---|---|
| `CornerAveraged` | `1` |
| `Trilinear` | `0` |

### Solution Convergence Modes

Selector: `solver.yml -> solution_convergence.mode`

| Value | Maps to |
|---|---|
| `periodic_deterministic` | `PERIODIC_DETERMINISTIC` |
| `statistical_steady` | `STATISTICAL_STEADY` |
| `steady_deterministic` | `STEADY_DETERMINISTIC` |
| `transient` | `TRANSIENT` |

### RANS Models

Selector: `case.yml -> models.physics.turbulence.rans.model`

| Value | Maps to |
|---|---|
| `disabled` | `0` |
| `k_omega` | `1` |
| `komega` | `1` |
| `none` | `0` |
| `off` | `0` |

### LES Test Filters

Selector: `case.yml -> models.physics.turbulence.les.test_filter.kernel`

| Value | Maps to |
|---|---|
| `box` | `0` |
| `simpson_ik` | `1` |
| `volume_weighted_box` | `0` |

### LES Grid Filter Width Models

Selector: `case.yml -> models.physics.turbulence.les.filter_width`

| Value | Maps to |
|---|---|
| `cube_root_volume` | `0` |
| `geometric_mean` | `1` |
| `max_edge` | `2` |

### LES Coefficient Averaging Modes

Selector: `case.yml -> models.physics.turbulence.les.averaging.mode`

| Value | Maps to |
|---|---|
| `global` | `2` |
| `homogeneous` | `1` |
| `local` | `0` |

### LES Coefficient Limiting Modes

Selector: `case.yml -> models.physics.turbulence.les.clipping.mode`

| Value | Maps to |
|---|---|
| `clamp` | `0` |
| `clip_negative` | `1` |
| `none` | `2` |

### Analytical Solution Types

Selector: `solver.yml -> operation_mode.analytical_type`

| Value | Maps to |
|---|---|
| `TGV3D` | `TGV3D` |
| `UNIFORM_FLOW` | `UNIFORM_FLOW` |
| `ZERO_FLOW` | `ZERO_FLOW` |

### Wall Function Models

Selector: `case.yml -> models.physics.turbulence.wall_function.model`

| Value | Maps to |
|---|---|
| `log_law` | `log_law` |
| `loglaw` | `loglaw` |

### Initial Condition Target Field

Selector: `case.yml -> properties.initial_conditions.field`

| Value | Maps to |
|---|---|
| `ucat` | `ucat` |
| `ucont` | `ucont` |

### Particle Statistics Tasks

Selector: `post.yml -> statistics pipeline task`

| Value | Maps to |
|---|---|
| `msd` | `msd` |

### Grid Ingestion Modes

Selector: `case.yml -> grid.mode`

| Value | Maps to |
|---|---|
| `file` | `file` |
| `grid_gen` | `grid_gen` |
| `programmatic_c` | `programmatic_c` |

### Grid Generator Geometries

Selector: `case.yml -> grid.generator.grid_type`

| Value | Maps to |
|---|---|
| `cpipe` | `cpipe` |
| `pipe` | `pipe` |
| `warp` | `warp` |

### Particle Restart Modes

Selector: `case.yml -> models.physics.particles.restart_mode`

| Value | Maps to |
|---|---|
| `init` | `init` |
| `load` | `load` |

### Eulerian Post-Processing Kernels

Selector: `post.yml -> eulerian_pipeline[].task`

| Value | Maps to |
|---|---|
| `nodal_average` | `nodal_average` |
| `normalize_field` | `normalize_field` |
| `q_criterion` | `q_criterion` |

### Lagrangian Post-Processing Kernels

Selector: `post.yml -> lagrangian_pipeline[].task`

| Value | Maps to |
|---|---|
| `specific_ke` | `specific_ke` |

### Spectra Tasks

Selector: `post.yml -> spectra.tasks[].task`

| Value | Maps to |
|---|---|
| `shell_spectrum` | `shell_spectrum` |

### Field Statistics Derived Outputs

Selector: `post.yml -> field_statistics.outputs[]`

| Value | Maps to |
|---|---|
| `flux` | `flux` |
| `mean` | `mean` |
| `reynolds_stress` | `reynolds_stress` |
| `rms` | `rms` |
| `tke` | `tke` |

### Field Statistics Weighting Modes

Selector: `monitor.yml -> field_statistics.windows[].weighting`

| Value | Maps to |
|---|---|
| `physical_time` | `physical_time` |
| `sample` | `sample` |

### Study Types

Selector: `study.yml -> study_type`

| Value | Maps to |
|---|---|
| `grid_independence` | `grid_independence` |
| `sensitivity` | `sensitivity` |
| `timestep_independence` | `timestep_independence` |

### Newton-Krylov Preconditioner Models

Selector: `solver.yml -> momentum_solver.newton_krylov.preconditioner.model`

| Value | Maps to |
|---|---|
| `frozen_momentum_jacobian` | `frozen_momentum_jacobian` |
| `none` | `none` |

### Storage Compression Policies

Selector: `picurv storage --compression / storage.yml -> compression`

| Value | Maps to |
|---|---|
| `auto` | `auto` |
| `balanced` | `balanced` |
| `fast` | `fast` |
| `maximum` | `maximum` |
| `none` | `none` |
