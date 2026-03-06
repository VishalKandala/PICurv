@page 19_Nondimensionalization Non-Dimensionalization Model

@anchor _Nondimensionalization

PICurv uses a dimensional-input, non-dimensional-core workflow:

1. Users specify physical quantities in YAML.
2. `picurv` converts/control-maps these to non-dimensional solver inputs.
3. Solver evolves non-dimensional fields.
4. Postprocessing can re-dimensionalize outputs for visualization/analysis.

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

- YAML validation and conversion: `scripts/picurv`
- Runtime option ingestion: `src/setup.c`
- Grid/data read-write scaling hooks: `src/io.c`
- Post re-dimensionalization trigger via `post.yml`:
  - `global_operations.dimensionalize: true`
  - pipeline keyword `DimensionalizeAllLoadedFields`

@section p19_notes_sec 5. Practical Notes

- Choose physically meaningful `L_ref` and `U_ref`; they directly affect `Re` and `dt*`.
- Keep input units consistent across all physical parameters in `case.yml`.
- For reproducibility, keep run-local generated artifacts under `runs/<run_id>/config/`.

@section p19_links_sec 6. References

- PETSc docs (solver stack context): https://petsc.org/release/docs/manual/
- Config contract details: **@subpage 14_Config_Contract**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Non-Dimensionalization Model** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

Treat this page as both a conceptual reference and a runbook. If you are debugging, pair the method/procedure described here with monitor output, generated runtime artifacts under `runs/<run_id>/config`, and the associated solver/post logs so numerical intent and implementation behavior stay aligned.

### What To Extract Before Changing A Case

- Identify which YAML role or runtime stage this page governs.
- List the primary control knobs (tolerances, cadence, paths, selectors, or mode flags).
- Record expected success indicators (convergence trend, artifact presence, or stable derived metrics).
- Record failure signals that require rollback or parameter isolation.

### Practical CFD Troubleshooting Pattern

1. Reproduce the issue on a tiny case or narrow timestep window.
2. Change one control at a time and keep all other roles/configs fixed.
3. Validate generated artifacts and logs after each change before scaling up.
4. If behavior remains inconsistent, compare against a known-good baseline example and re-check grid/BC consistency.

