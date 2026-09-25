- Every configuration input is now physical and is converted to solver units exactly
  once; post-processing reports every output physically under `dimensionalize`. This
  changes what several inputs mean when the reference scales are not 1 - `picurv
  validate` names them - and fixes outputs that mixed unit systems.
  - Now read in physical units: `particles.point_source`,
    `wall_function.roughness_height`, `solver.yml` `uniform_flow`, the
    `verification.sources` profiles, `field_statistics` window `start_time`,
    `end_time` and `time_cadence`, `ic_gen` expressions (physical coordinates and
    values), the spectral provider velocities and wavenumbers, and a `mode: file`
    initial condition. A file saved by an earlier run declares its scales with
    `source_case`, or `velocity_scale` and `length_scale`.
  - Every input records its dimension and conversion site in `INPUT_QUANTITIES`, and
    every catalogued field a `FIELD_DIM_*` dimension; `units.nondimensionalization` is
    now an enforced contract, checked by `audit_units.py` against page 19's input and
    field indexes.
  - Loaded fields are scaled as they are read, so a field a step does not reload is no
    longer rescaled; `Nu_t`, `Diffusivity`, `Phi` and the other catalogued fields now
    have scales instead of staying non-dimensional with a warning. The
    `DimensionalizeAllLoadedFields` pipeline stage is retired.
  - MSD columns, spectra and field-statistics history times, and ParaView collection
    times are physical under `dimensionalize`; runtime diagnostics CSVs gain a
    `physical_time` column.
