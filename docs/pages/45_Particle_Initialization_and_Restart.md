@page 45_Particle_Initialization_and_Restart Particle Initialization and Restart Guide

@anchor _Particle_Initialization_and_Restart

This page documents particle seeding, restart behavior, and early-step migration/settling logic.
It is written for both case authors and contributors working in `src/ParticleSwarm.c` and `src/ParticleMotion.c`.

@tableofcontents

@section p45_contract_sec 1. User Contract in `case.yml`

Particle controls live in:

```yaml
models:
  physics:
    particles:
      count: 50000
      init_mode: "Surface"   # Surface | Volume | PointSource | SurfaceEdges
      restart_mode: "init"   # init | load
      point_source:
        x: 0.5
        y: 0.5
        z: 0.5
```

Mapping to control flags:

- `count` -> `-numParticles`
- `init_mode` -> `-pinit`
- `restart_mode` -> `-particle_restart_mode`
- `point_source` -> `-psrc_x/-psrc_y/-psrc_z` (required when `init_mode` is `PointSource`)

@section p45_mode_values_sec 2. Accepted `init_mode` Values

`picurv` accepts these exact canonical strings:

- `Surface`
- `Volume`
- `PointSource`
- `SurfaceEdges`

Enum mapping in C (`ParticleInitializationType`):

- `0`: `PARTICLE_INIT_SURFACE_RANDOM`
- `1`: `PARTICLE_INIT_VOLUME`
- `2`: `PARTICLE_INIT_POINT_SOURCE`
- `3`: `PARTICLE_INIT_SURFACE_EDGES`

@section p45_mode_behavior_sec 3. Mode Behavior in C

Main setup flow:

1. @ref InitializeParticleSwarm creates the DMSwarm and fields.
2. `AssignInitialPropertiesToSwarm` seeds base particle state.
3. `PerformInitializedParticleSetup` settles/migrates and couples to Eulerian fields.

Mode details:

`Surface` (`0`):

- requires an identified INLET face from BC parsing,
- ranks not servicing inlet place particles at inlet-center fallback (`CMx_c/CMy_c/CMz_c`) before migration,
- `ReinitializeParticlesOnInletSurface` re-spreads particles on inlet partitions after first settlement.

`SurfaceEdges` (`3`):

- deterministic placement by particle ID on inlet-face lattice,
- if deterministic target is non-local after migration, code falls back to random inlet placement.

`Volume` (`1`):

- random logical coordinates inside locally owned cells,
- mapped to physical space via metric interpolation.

`PointSource` (`2`):

- all particles start at fixed coordinates `(psrc_x, psrc_y, psrc_z)`.

@section p45_restart_matrix_sec 4. Restart Behavior Matrix

`InitializeParticleSwarm` behavior depends on `StartStep` and restart mode:

| `StartStep` | `particle_restart_mode` | Behavior |
|---|---|---|
| `0` | any | initialize new population |
| `>0` | `init` | initialize new population in restarted flow |
| `>0` | `load` | load particle fields from restart files |

Operational note:

- For a run completed through step `N`, use `start_step: N`.
- Choose `restart_mode: load` to continue the existing particle swarm.
- Choose `restart_mode: init` to reseed a fresh particle population in the restarted flow field.

For loaded particles, fast migration path:

- @ref MigrateRestartParticlesUsingCellID uses stored cell IDs to migrate directly before full walking-search fallback.

@section p45_settle_sec 5. Early-Step Settlement and Coupling

For initialized particles (`StartStep == 0` path):

1. `LocateAllParticlesInGrid` performs location/migration.
2. surface modes call @ref ReinitializeParticlesOnInletSurface.
3. statuses are reset and location pass is repeated.
4. `InterpolateAllFieldsToSwarm` assigns flow fields at particle positions.
5. optional scatter updates Eulerian particle-derived fields.

For loaded particles:

1. `MigrateRestartParticlesUsingCellID` fast migration,
2. `LocateAllParticlesInGrid` resolves invalid/missing cases,
3. interpolation/scatter synchronize coupling state before stepping.

@section p45_fields_sec 6. Swarm Fields Initialized at Startup

After position/PID/cell placeholders, initialization sets defaults for:

- `velocity` (vector),
- `weight` (vector),
- `Diffusivity` (scalar),
- `DiffusivityGradient` (vector),
- `Psi` (scalar).

Cell IDs start at `-1` until location confirms host cells.

@section p45_diagnostics_sec 7. Diagnostics and Sanity Checks

Check banner/log output for:

- selected particle initialization mode,
- identified inlet face for surface modes,
- migration/lost-particle counters after first steps.

Typical errors:

- no INLET face with surface modes,
- missing `point_source.{x,y,z}` for point source mode,
- restart mode not in `{init, load}`.

@section p45_extension_sec 8. Contributor Extension Points

If adding a new particle initialization mode:

1. add validation and mapping in `scripts/picurv`,
2. extend `ParticleInitializationType` and parser wiring in C,
3. implement placement logic in `InitializeParticleBasicProperties` and any inlet reinit path,
4. update logging string mappings and tests,
5. update this page and related method references.

For the full selector extension checklist, see **@subpage 50_Modular_Selector_Extension_Guide**.

@section p45_refs_sec 9. Related Pages

- **@subpage 33_Initial_Conditions**
- **@subpage 34_Particle_Model_Overview**
- **@subpage 52_Run_Lifecycle_Guide**
- **@subpage 26_Walking_Search_Method**
- **@subpage 27_Trilinear_Interpolation_and_Projection**
- **@subpage 39_Common_Fatal_Errors**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Particle Initialization and Restart Guide** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
