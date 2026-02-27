@page 19_Nondimensionalization Non-Dimensionalization Model

PICurv uses a dimensional-input, non-dimensional-core workflow:

1. Users specify physical quantities in YAML.
2. `pic.flow` converts/control-maps these to non-dimensional solver inputs.
3. Solver evolves non-dimensional fields.
4. Postprocessing can re-dimensionalize outputs for visualization/analysis.

@tableofcontents

@section refs_sec 1. Reference Scales

From `case.yml`:
- `L_ref = properties.scaling.length_ref`
- `U_ref = properties.scaling.velocity_ref`
- `rho_ref = properties.fluid.density`

Derived scales:
- `T_ref = L_ref / U_ref`
- `P_ref = rho_ref * U_ref^2`

@section primary_sec 2. Primary Converted Quantities

Common quantities prepared before solver launch:

- Reynolds number:
  - `Re = rho_ref * U_ref * L_ref / mu`
  - emitted as `-ren`
- Physical timestep to non-dimensional timestep:
  - `dt* = dt_physical / T_ref`
  - emitted as `-dt`

For file/grid-generated meshes, coordinates are normalized by `L_ref` before C-side ingestion.

@section fields_sec 3. Field/Coordinate Scaling Conventions

Representative scaling conventions used across solver/post:

- Position: `x* = x / L_ref`
- Velocity: `u* = u / U_ref`
- Pressure: `p* = p / P_ref`

Particle and Eulerian data follow the same scale family so interpolation/projection remain consistent.

@section pipeline_sec 4. Where This Happens in the Pipeline

- YAML validation and conversion: `scripts/pic.flow`
- Runtime option ingestion: `src/setup.c`
- Grid/data read-write scaling hooks: `src/io.c`
- Post re-dimensionalization trigger via `post.yml`:
  - `global_operations.dimensionalize: true`
  - pipeline keyword `DimensionalizeAllLoadedFields`

@section notes_sec 5. Practical Notes

- Choose physically meaningful `L_ref` and `U_ref`; they directly affect `Re` and `dt*`.
- Keep input units consistent across all physical parameters in `case.yml`.
- For reproducibility, keep run-local generated artifacts under `runs/<run_id>/config/`.

@section links_sec 6. References

- PETSc docs (solver stack context): https://petsc.org/release/docs/manual/
- Config contract details: **@subpage 14_Config_Contract**
