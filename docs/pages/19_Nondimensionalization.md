@page 19_Nondimensionalization Non-Dimensionalization Model

@anchor _Nondimensionalization

PICurv uses a dimensional-input, non-dimensional-core workflow:

1. Users specify physical quantities in YAML.
2. `picurv` converts/control-maps these to non-dimensional solver inputs.
3. Solver evolves non-dimensional fields.
4. Postprocessing can re-dimensionalize outputs for <run.visualization>/analysis.

@tableofcontents

@section p19_refs_sec 1. Reference Scales

From `case.yml`:
- `L_ref = properties.scaling.length_ref`
- `U_ref = properties.scaling.velocity_ref`
- `rho_ref = properties.fluid.density`

Derived scales:
- `T_ref = L_ref / U_ref`
- `P_ref = rho_ref * U_ref^2`

@section p19_primary_sec 2. Primary Converted Quantities

Common quantities prepared before solver launch:

- Reynolds number:
  - `Re = rho_ref * U_ref * L_ref / mu`
  - emitted as `-ren`
- Physical timestep to non-dimensional timestep:
  - `dt* = dt_physical / T_ref`
  - emitted as `-dt`

For file/grid-generated meshes, coordinates are normalized by `L_ref` before C-side ingestion.

@section p19_fields_sec 3. Field/Coordinate Scaling Conventions

Representative scaling conventions used across solver/post:

- Position: `x* = x / L_ref`
- Velocity: `u* = u / U_ref`
- Pressure: `p* = p / P_ref`

Particle and Eulerian data follow the same scale family so interpolation/projection remain consistent.

@section p19_pipeline_sec 4. Where This Happens in the Pipeline

- YAML validation and conversion: `picurv_cli/core.py`
- Runtime option ingestion: `src/setup.c`
- Grid/data read-write scaling hooks: `src/io.c`
- Post re-dimensionalization trigger via `post.yml`:
  - `global_operations.dimensionalize: true`
  - pipeline keyword `DimensionalizeAllLoadedFields`

@subsection p19_derived_sub 4.1 Derived Products Carry Their Own Exponents

`global_operations.dimensionalize` reaches three producers, because a run's derived
output is not made in one place:

| Producer | What it scales |
|---|---|
| Field pipeline (`DimensionalizeAllLoadedFields`) | Loaded fields and grid coordinates, each by its own reference scale |
| Field-statistics accumulator | Derived statistics, by the source field's scale raised to the kind's exponent |
| Spectra generator | Every reported spectral quantity, by its own combination of `U_ref` and `L_ref` |

A single blanket factor cannot serve them, because a derived quantity's dimensions are
not its source field's:

- mean, RMS: source scale to the first power
- Reynolds stress, turbulent kinetic energy: to the second
- co-moment flux: the product of the two fields' scales, which need not be equal
- wavenumber: `1 / L_ref`; energy: `U_ref^2`; spectral second moment: `U_ref^2 / L_ref^2`;
  integral and Taylor scales: `L_ref`
- Parseval residual and shell counts are ratios and counts, and are never scaled

Checkpoints are never dimensionalized: they are restart state, and restart must reproduce
the solver's own units exactly. Only derived output is scaled, so a `.dat` under
`<run.checkpoints>` and a `.vts` under `<run.visualization>` are deliberately in different
unit systems.

The grid is the one quantity converted on the way *in* rather than out, in all three
`grid.mode` values, so the round trip is uniform regardless of which one a case uses:

- `file` and `grid_gen`: a `.picgrid` is supplied dimensional and is divided by
  `L_ref` when the conductor publishes it as an asset, so the solver reads geometry
  already in solver units.
- `programmatic_c`: `xMins`/`xMaxs`/`yMins`/... are likewise supplied dimensional.
  `ReadGridGenerationInputs()` (`src/io.c`) divides each bound by `L_ref` after
  `ParseScalingInformation()` has resolved it, immediately before assigning
  `user->Min_X`/`Max_X`/etc - the same division a `file` grid receives in Python, just
  performed in C because this mode has no Python staging step to perform it in.

Coordinates written by post-processing are multiplied back by `L_ref` unconditionally
(`DimensionalizeField("Coordinates")`), so a case authored with `length_ref: 1.0` is
unaffected either way, and a case using another reference length gets back the same
physical bounds it specified, whichever grid mode produced them.

@section p19_notes_sec 5. Practical Notes

- Choose physically meaningful `L_ref` and `U_ref`; they directly affect `Re` and `dt*`.
- Keep input units consistent across all physical parameters in `case.yml`.
- For reproducibility, keep run-local generated artifacts under `<run.config>/`.

@section p19_links_sec 6. References

- PETSc docs (solver stack context): https://petsc.org/release/docs/manual/
- Config contract details: **@subpage 14_Config_Contract**
