@page 34_Particle_Model_Overview Particle Model and Coupling Overview

@anchor _Particle_Model_Overview

This page documents the particle pipeline exactly as orchestrated in the current solver flow.

@tableofcontents

@section p34_loop_sec 1. Per-Step Particle Pipeline

Typical order inside coupled step:

1. interpolate Eulerian fields to swarm: @ref InterpolateAllFieldsToSwarm
2. advance particle positions: @ref UpdateAllParticlePositions
3. locate/migrate particles: @ref LocateAllParticlesInGrid
4. update particle physics fields: @ref UpdateAllParticleFields
5. scatter particle fields back to Eulerian storage: @ref ScatterAllParticleFieldsToEulerFields

Each stage relies on valid `DMSwarm_CellID`, interpolation weights, and synchronized ghost data.

Step 1 dispatches through @ref InterpolateEulerFieldToSwarm, which selects between two interpolation paths based on `SimCtx.interpolationMethod`:

- **Trilinear** (default): direct interpolation from the 8 nearest cell centers, second-order on curvilinear grids.
- **CornerAveraged** (legacy): two-stage center-to-corner averaging then trilinear from corners.

See **@subpage 27_Trilinear_Interpolation_and_Projection** for method details and configuration.

@section p34_fields_sec 2. Core Particle Fields In Use

Commonly used swarm fields include:

- position (`position`)
- velocity (`velocity`)
- cell ID (`DMSwarm_CellID`)
- location status (`DMSwarm_location_status`)
- diffusivity and scalar (`Diffusivity`, `Psi`)

Status transitions (`NEEDS_LOCATION`, `ACTIVE_AND_LOCATED`, `MIGRATING_OUT`, `LOST`) determine whether particles are walked, migrated, or skipped in settlement passes.

@section p34_physics_sec 3. Current Built-In Physics

Current scalar model path:

- kernel: @ref UpdateParticleField
- batch loop: @ref UpdateFieldForAllParticles
- orchestrator: @ref UpdateAllParticleFields

This presently implements IEM-style relaxation for `Psi`, with model constants sourced from runtime context.

@section p34_statistics_sec 4. Statistics and Diagnostics

Post statistics currently include global kernels such as:

- @ref ComputeParticleMSD

Additional health indicators are available from migration counters and settlement-pass counts stored in `SimCtx` fields updated by location logic.

@section p34_extension_sec 5. Extending To New Closures

Recommended extension pattern:

1. add particle model selector to YAML/flags,
2. branch inside particle physics orchestrator,
3. isolate model-specific kernels from motion/location kernels,
4. keep scatter/interpolation interfaces unchanged,
5. validate with smoke tests and deterministic small cases.

For configuration contract changes, update:

- **@subpage 14_Config_Contract**
- **@subpage 16_Config_Extension_Playbook**
- **@subpage 40_Testing_and_Quality_Guide**

@section p34_refs_sec 6. Related Pages

- **@subpage 33_Initial_Conditions**
- **@subpage 45_Particle_Initialization_and_Restart**
- **@subpage 26_Walking_Search_Method**

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Particle Model and Coupling Overview** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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

