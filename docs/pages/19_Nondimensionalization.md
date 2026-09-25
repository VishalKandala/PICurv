@page 19_Nondimensionalization Units and Non-Dimensionalization

@anchor _Nondimensionalization

PICurv follows one units rule, everywhere:

1. Every value a user writes in `case.yml`, `solver.yml`, or `monitor.yml` is physical,
   in any consistent unit system, the one the reference scales are written in.
2. Each value is converted to solver units exactly once, at the point where the
   reference scales are known.
3. The solver works in non-dimensional form, and checkpoints stay in that form.
4. Output is physical when a post-processing recipe sets
   `global_operations.dimensionalize: true`, and in solver units otherwise.

Choosing `length_ref: 1.0` and `velocity_ref: 1.0` makes physical and solver values
coincide, which is what most test and smoke cases do; the rule still holds, and a case
may change its scales without rewriting any other input.

@tableofcontents

@section p19_refs_sec 1. Reference Scales and Dimensions

From `case.yml`:

- `L = properties.scaling.length_ref`
- `U = properties.scaling.velocity_ref`
- `rho = properties.fluid.density`

A quantity's dimension is recorded as three exponents `(a, b, c)` of those scales; a
physical value is divided by `L^a U^b rho^c` to reach solver units, and a solver value is
multiplied by it to leave them. The launcher (`picurv_cli/core.py`) records input
dimensions with these triples, and the C field catalogs record field dimensions with the
same triples, so input and output conversion share one vocabulary:

| Dimension | Exponents `(L, U, rho)` | Scale | Examples |
|---|---|---|---|
| `LENGTH` | (1, 0, 0) | `L` | grid coordinates, particle position, roughness height |
| `AREA` | (2, 0, 0) | `L^2` | face-area metrics `Csi`, `Eta`, `Zet` |
| `INVERSE_VOLUME` | (-3, 0, 0) | `1/L^3` | Jacobian `Aj` |
| `WAVENUMBER` | (-1, 0, 0) | `1/L` | spectral `k0`, verification `kx` |
| `VELOCITY` | (0, 1, 0) | `U` | `Ucat`, BC velocities, particle velocity |
| `TIME` | (1, -1, 0) | `L/U` | `dt_physical`, statistics window times |
| `VOLUME_FLUX` | (2, 1, 0) | `U L^2` | `Ucont`, `target_flux` |
| `DIFFUSIVITY` | (1, 1, 0) | `U L` | `Nu_t`, `Diffusivity` |
| `PRESSURE` | (0, 2, 1) | `rho U^2` | `P`, `Phi` |
| `INVERSE_TIME_SQUARED` | (-2, 2, 0) | `U^2/L^2` | Q-criterion |
| `DIMENSIONLESS` | (0, 0, 0) | 1 | `Psi`, `CS`, `Nvert`, model constants |

The Reynolds number is `Re = rho U L / mu` and is passed as `-ren`; `T = L/U` and
`P = rho U^2` are the derived time and pressure scales.

@section p19_pipeline_sec 2. Where Conversion Happens

Each physical input is converted at one site, recorded beside its dimension:

| Site | Meaning | Inputs |
|---|---|---|
| `cli` | the launcher divides the value while writing the control file | `dt_physical`, BC velocities and fluxes, `_physical` IC values, point source, roughness height, `uniform_flow`, verification profiles, window times, spectral provider parameters |
| `c` | the runtime divides it, because the input reaches C with no staging step | `programmatic_settings` domain bounds (`ReadGridGenerationInputs()` in `src/io.c`) |
| `staging` | a file payload is divided as it is staged under `<run.inputs>` | `.picgrid` grids, `grid.gen` output, inlet-profile PICSLICEs, a `mode: file` initial condition |
| `provider` | a generator is handed the reference scales and writes solver units | `ic_gen` expressions (`--length-ref`, `--velocity-ref`) |
| `reference` | the input defines a scale | `length_ref`, `velocity_ref`, `density`, `viscosity`, the provenance scales of a file payload |
| `passthrough` | raw solver flags, in solver units by definition | `solver_parameters`, `petsc_passthrough_options` |

The dimension and site of every input are recorded in `INPUT_QUANTITIES` (with
`BC_PARAM_QUANTITIES` and `IC_PARAM_QUANTITIES` for free-form parameter mappings) in
`picurv_cli/core.py`, and every conversion takes its factor from that record.

@section p19_inputs_sec 3. Input Index

Every input that carries a physical dimension. Selectors, counts, steps, grid indices,
tolerances, and dimensionless model constants are recorded as such and are not listed.
`PAYLOAD` marks a file whose dimension follows its content: velocity for `Ucat`, volume
flux for `Ucont`, and physical coordinates and values for `ic_gen` expressions.

| Input | Dimension | Converted by | Scale |
|---|---|---|---|
| `case.yml: boundary_conditions[].params.source` | VELOCITY | staging | U |
| `case.yml: boundary_conditions[].params.target_flux` | VOLUME_FLUX | cli | U L^2 |
| `case.yml: boundary_conditions[].params.v_max` | VELOCITY | cli | U |
| `case.yml: boundary_conditions[].params.vx` | VELOCITY | cli | U |
| `case.yml: boundary_conditions[].params.vy` | VELOCITY | cli | U |
| `case.yml: boundary_conditions[].params.vz` | VELOCITY | cli | U |
| `case.yml: grid.generator` | LENGTH | staging | L |
| `case.yml: grid.programmatic_settings.xMaxs` | LENGTH | c | L |
| `case.yml: grid.programmatic_settings.xMins` | LENGTH | c | L |
| `case.yml: grid.programmatic_settings.yMaxs` | LENGTH | c | L |
| `case.yml: grid.programmatic_settings.yMins` | LENGTH | c | L |
| `case.yml: grid.programmatic_settings.zMaxs` | LENGTH | c | L |
| `case.yml: grid.programmatic_settings.zMins` | LENGTH | c | L |
| `case.yml: grid.source_file` | LENGTH | staging | L |
| `case.yml: initial_conditions.params.bulk_velocity (channel_spectral_velocity)` | VELOCITY | cli | U |
| `case.yml: initial_conditions.params.bulk_velocity (duct_spectral_velocity)` | VELOCITY | cli | U |
| `case.yml: initial_conditions.params.config_file (ic_gen)` | PAYLOAD | provider | by payload field |
| `case.yml: initial_conditions.params.normalization.target (spectral_random_velocity)` | VELOCITY | cli | U |
| `case.yml: initial_conditions.params.peak_velocity_physical (poiseuille)` | VELOCITY | cli | U |
| `case.yml: initial_conditions.params.perturbation_rms (channel_spectral_velocity)` | VELOCITY | cli | U |
| `case.yml: initial_conditions.params.perturbation_rms (duct_spectral_velocity)` | VELOCITY | cli | U |
| `case.yml: initial_conditions.params.random.mean (spectral_random_velocity)` | VELOCITY | cli | U |
| `case.yml: initial_conditions.params.spectrum.k0 (channel_spectral_velocity)` | WAVENUMBER | cli | 1/L |
| `case.yml: initial_conditions.params.spectrum.k0 (duct_spectral_velocity)` | WAVENUMBER | cli | 1/L |
| `case.yml: initial_conditions.params.spectrum.k0 (spectral_random_velocity)` | WAVENUMBER | cli | 1/L |
| `case.yml: initial_conditions.params.spectrum.k_cut (channel_spectral_velocity)` | WAVENUMBER | cli | 1/L |
| `case.yml: initial_conditions.params.spectrum.k_cut (duct_spectral_velocity)` | WAVENUMBER | cli | 1/L |
| `case.yml: initial_conditions.params.spectrum.k_cut (spectral_random_velocity)` | WAVENUMBER | cli | 1/L |
| `case.yml: initial_conditions.params.u_physical (constant)` | VELOCITY | cli | U |
| `case.yml: initial_conditions.params.v_physical (constant)` | VELOCITY | cli | U |
| `case.yml: initial_conditions.params.velocity_physical (constant)` | VELOCITY | cli | U |
| `case.yml: initial_conditions.params.velocity_physical (streamwise_constant)` | VELOCITY | cli | U |
| `case.yml: initial_conditions.params.w_physical (constant)` | VELOCITY | cli | U |
| `case.yml: models.physics.particles.point_source.x` | LENGTH | cli | L |
| `case.yml: models.physics.particles.point_source.y` | LENGTH | cli | L |
| `case.yml: models.physics.particles.point_source.z` | LENGTH | cli | L |
| `case.yml: models.physics.turbulence.wall_function.roughness_height` | LENGTH | cli | L |
| `case.yml: properties.fluid.density` | DENSITY | reference | rho |
| `case.yml: properties.fluid.viscosity` | DYNAMIC_VISCOSITY | reference | rho U L |
| `case.yml: properties.initial_conditions.length_scale` | LENGTH | reference | L |
| `case.yml: properties.initial_conditions.peak_velocity_physical` | VELOCITY | cli | U |
| `case.yml: properties.initial_conditions.source_file` | PAYLOAD | staging | by payload field |
| `case.yml: properties.initial_conditions.u_physical` | VELOCITY | cli | U |
| `case.yml: properties.initial_conditions.v_physical` | VELOCITY | cli | U |
| `case.yml: properties.initial_conditions.velocity_physical` | VELOCITY | cli | U |
| `case.yml: properties.initial_conditions.velocity_scale` | VELOCITY | reference | U |
| `case.yml: properties.initial_conditions.w_physical` | VELOCITY | cli | U |
| `case.yml: properties.scaling.length_ref` | LENGTH | reference | L |
| `case.yml: properties.scaling.velocity_ref` | VELOCITY | reference | U |
| `case.yml: run_control.dt_physical` | TIME | cli | L/U |
| `monitor.yml: field_statistics.windows.[].end_time` | TIME | cli | L/U |
| `monitor.yml: field_statistics.windows.[].start_time` | TIME | cli | L/U |
| `monitor.yml: field_statistics.windows.[].time_cadence` | TIME | cli | L/U |
| `solver.yml: operation_mode.uniform_flow.u` | VELOCITY | cli | U |
| `solver.yml: operation_mode.uniform_flow.v` | VELOCITY | cli | U |
| `solver.yml: operation_mode.uniform_flow.w` | VELOCITY | cli | U |
| `solver.yml: verification.sources.diffusivity.gamma0` | DIFFUSIVITY | cli | U L |
| `solver.yml: verification.sources.diffusivity.slope_x` | VELOCITY | cli | U |
| `solver.yml: verification.sources.scalar.kx` | WAVENUMBER | cli | 1/L |
| `solver.yml: verification.sources.scalar.ky` | WAVENUMBER | cli | 1/L |
| `solver.yml: verification.sources.scalar.kz` | WAVENUMBER | cli | 1/L |
| `solver.yml: verification.sources.scalar.slope_x` | WAVENUMBER | cli | 1/L |

@section p19_fields_sec 4. Field Index

Every catalogued field records the dimension post-processing scales it by, as a
`FIELD_DIM_*` argument of its entry in `src/field_catalog.c` or
`src/particle_field_catalog.c`. `NOT_A_QUANTITY` marks integer bookkeeping, and
`FROM_SOURCE` marks staging storage that holds whatever field was written into it; neither
has a scale of its own, and asking to dimensionalize one is an error. Each field's scale is
resolved from this record by `PicurvFieldReferenceScale()` (`src/io.c`).

| Field | Catalog | Dimension | Scale |
|---|---|---|---|
| `Aj` | Eulerian | `INVERSE_VOLUME` | 1/L^3 |
| `CellScalarAtCorner` | Eulerian | `FROM_SOURCE` | its source field's |
| `CellVectorAtCorner` | Eulerian | `FROM_SOURCE` | its source field's |
| `Cent` | Eulerian | `LENGTH` | L |
| `Centx` | Eulerian | `LENGTH` | L |
| `Centy` | Eulerian | `LENGTH` | L |
| `Centz` | Eulerian | `LENGTH` | L |
| `Coordinates` | Eulerian | `LENGTH` | L |
| `CS` | Eulerian | `DIMENSIONLESS` | 1 |
| `Csi` | Eulerian | `AREA` | L^2 |
| `Diffusivity` | Eulerian | `DIFFUSIVITY` | U L |
| `DiffusivityGradient` | Eulerian | `VELOCITY` | U |
| `Eta` | Eulerian | `AREA` | L^2 |
| `GridSpace` | Eulerian | `LENGTH` | L |
| `IAj` | Eulerian | `INVERSE_VOLUME` | 1/L^3 |
| `ICsi` | Eulerian | `AREA` | L^2 |
| `IEta` | Eulerian | `AREA` | L^2 |
| `IZet` | Eulerian | `AREA` | L^2 |
| `JAj` | Eulerian | `INVERSE_VOLUME` | 1/L^3 |
| `JCsi` | Eulerian | `AREA` | L^2 |
| `JEta` | Eulerian | `AREA` | L^2 |
| `JZet` | Eulerian | `AREA` | L^2 |
| `KAj` | Eulerian | `INVERSE_VOLUME` | 1/L^3 |
| `KCsi` | Eulerian | `AREA` | L^2 |
| `KEta` | Eulerian | `AREA` | L^2 |
| `KZet` | Eulerian | `AREA` | L^2 |
| `Nu_t` | Eulerian | `DIFFUSIVITY` | U L |
| `NuWall` | Eulerian | `DIFFUSIVITY` | U L |
| `Nvert` | Eulerian | `DIMENSIONLESS` | 1 |
| `Nvert_o` | Eulerian | `DIMENSIONLESS` | 1 |
| `P` | Eulerian | `PRESSURE` | rho U^2 |
| `ParticleCount` | Eulerian | `DIMENSIONLESS` | 1 |
| `Phi` | Eulerian | `PRESSURE` | rho U^2 |
| `PostScalar` | Eulerian | `FROM_SOURCE` | its source field's |
| `PostVector` | Eulerian | `FROM_SOURCE` | its source field's |
| `Psi` | Eulerian | `DIMENSIONLESS` | 1 |
| `Qcrit` | Eulerian | `INVERSE_TIME_SQUARED` | U^2/L^2 |
| `Ucat` | Eulerian | `VELOCITY` | U |
| `Ucont` | Eulerian | `VOLUME_FLUX` | U L^2 |
| `Ucont_o` | Eulerian | `VOLUME_FLUX` | U L^2 |
| `Ucont_rm1` | Eulerian | `VOLUME_FLUX` | U L^2 |
| `Utau` | Eulerian | `VELOCITY` | U |
| `Zet` | Eulerian | `AREA` | L^2 |
| `Diffusivity` | particle | `DIFFUSIVITY` | U L |
| `DiffusivityGradient` | particle | `VELOCITY` | U |
| `DMSwarm_CellID` | particle | `NOT_A_QUANTITY` | not scaled |
| `DMSwarm_location_status` | particle | `NOT_A_QUANTITY` | not scaled |
| `DMSwarm_pid` | particle | `NOT_A_QUANTITY` | not scaled |
| `DMSwarm_rank` | particle | `NOT_A_QUANTITY` | not scaled |
| `position` | particle | `LENGTH` | L |
| `Psi` | particle | `DIMENSIONLESS` | 1 |
| `velocity` | particle | `VELOCITY` | U |
| `weight` | particle | `DIMENSIONLESS` | 1 |

@section p19_output_sec 5. Output

`global_operations.dimensionalize` is a setting that reaches every producer, and each
scales its own values once, where they are made:

| Producer | With `dimensionalize: true` |
|---|---|
| Eulerian fields | each field is scaled by its catalog dimension as it is read from a checkpoint (`ReadSimulationFields()`), so a field a step does not reload is never scaled twice |
| Grid coordinates | scaled by `L` once, after the first read |
| Particle fields | every loaded particle field with a fixed dimension, after the statistics pipeline has run |
| Particle MSD (`*_msd.csv`) | computed in solver units against `D = 1/(Re Sc)`; `t` is written times `L/U`, the `MSD_*` columns times `L^2`, radii and centre of mass times `L`; the relative error and shell fractions are ratios and are the same in either system |
| Field-statistics fields | a derived statistic carries its source field's scale raised to the kind's exponent (mean and RMS to the first, Reynolds stress and TKE to the second), and a co-moment the product of its two fields' scales |
| Field-statistics history (`*_statistics_<window>.csv`) | `represented_time` in seconds, and `total_weight` too under `physical_time` weighting |
| Spectra | every spectral quantity by its own combination of `U` and `L`; the `time` column times `L/U` |
| Q-criterion | `U^2/L^2`: the velocity is already physical and the metrics are not, so the kernel divides by `L^2` |
| ParaView collections (`.pvd`) | each frame's time is its checkpoint time times the reference time of the run that wrote it |

Without the setting, every one of these is in solver units, including the collection
times, which then match the frames they index.

Runtime diagnostics written during a solve (`les_coefficient.csv`, `wall_model.csv`,
`interpolation_error.csv`, `scatter_metrics.csv`, `search_metrics.csv`) are in solver
units, with a `physical_time` column beside `time` for plotting in seconds. Checkpoints
are never dimensionalized: they are restart state, and a restart must reproduce the
solver's own units exactly.

@section p19_payloads_sec 6. File Payloads and Providers

- **Grids.** A `.picgrid`, whether supplied or written by `grid.gen`, is physical and is
  divided by `L` when staged; `programmatic_settings` bounds are divided by `L` in C.
- **Inlet profiles.** A `file` or `generated` PICSLICE is a physical velocity. A
  `field_slice` re-samples an earlier run's solver-unit field, and names that run's case
  (`source_case`) or its `velocity_scale`.
- **File initial conditions.** A `mode: file` payload is physical: `Ucat` in velocity,
  `Ucont` in volume flux per face. A field saved by an earlier PICurv run is in that run's
  solver units, so it names the case that wrote it (`source_case`), or gives
  `velocity_scale`, with `length_scale` too for `Ucont`; the staged payload is rescaled
  from those scales to this case's.
- **`ic_gen` expressions** are evaluated at physical coordinates and give physical values.
  The generator receives `--length-ref` and `--velocity-ref` and writes solver units; a
  custom `params.script` receives the same two options.
- **Spectral providers** receive their velocities and wavenumbers already converted, and
  work on the staged non-dimensional grid; their summaries report solver units.

@section p19_transition_sec 7. Inputs Once Read in Solver Units

The point source, the wall-function roughness height, `uniform_flow`, the verification
source profiles, the statistics window times, the `ic_gen` expressions, the spectral
provider parameters, and a `mode: file` initial-condition payload were read in solver
units before this rule covered them. At unit reference scales the two readings are the
same. At any other scale, `picurv validate` names each of these inputs a configuration
sets, since a value written in solver units has to be converted to physical.

@section p19_extend_sec 8. Adding an Input or a Field

- A new configuration key needs an entry in `INPUT_QUANTITIES` (or `BC_PARAM_QUANTITIES`,
  `IC_PARAM_QUANTITIES`) giving its dimension and conversion site;
  `tests/test_units_and_scaling.py` fails until it has one, and a `cli` conversion takes
  its factor from that entry through `to_solver_units()`.
- A new catalogued field needs a `FIELD_DIM_*` argument; an entry without one does not
  compile.
- Both indexes on this page are checked against the code by
  `tests/tooling/audit_units.py`, the checker of the `units.nondimensionalization`
  contract.

@section p19_links_sec 9. References

- Config contract details: **@subpage 14_Config_Contract**
- Field identities and layouts: @ref 56_Field_Identity_and_Layout_Catalog
