@page 34_Particle_Model_Overview Particle Model and Coupling Overview

This page documents the particle pipeline exactly as orchestrated in the current solver flow.

@tableofcontents

@section loop_sec 1. Per-Step Particle Pipeline

Typical order inside coupled step:

1. interpolate Eulerian fields to swarm: @ref InterpolateAllFieldsToSwarm
2. advance particle positions: @ref UpdateAllParticlePositions
3. locate/migrate particles: @ref LocateAllParticlesInGrid
4. update particle physics fields: @ref UpdateAllParticleFields
5. scatter particle fields back to Eulerian storage: @ref ScatterAllParticleFieldsToEulerFields

Each stage relies on valid `DMSwarm_CellID`, interpolation weights, and synchronized ghost data.

@section fields_sec 2. Core Particle Fields In Use

Commonly used swarm fields include:

- position (`position`)
- velocity (`velocity`)
- cell ID (`DMSwarm_CellID`)
- location status (`DMSwarm_location_status`)
- diffusivity and scalar (`Diffusivity`, `Psi`)

Status transitions (`NEEDS_LOCATION`, `ACTIVE_AND_LOCATED`, `MIGRATING_OUT`, `LOST`) determine whether particles are walked, migrated, or skipped in settlement passes.

@section physics_sec 3. Current Built-In Physics

Current scalar model path:

- kernel: @ref UpdateParticleField
- batch loop: @ref UpdateFieldForAllParticles
- orchestrator: @ref UpdateAllParticleFields

This presently implements IEM-style relaxation for `Psi`, with model constants sourced from runtime context.

@section statistics_sec 4. Statistics and Diagnostics

Post statistics currently include global kernels such as:

- @ref ComputeParticleMSD

Additional health indicators are available from migration counters and settlement-pass counts stored in `SimCtx` fields updated by location logic.

@section extension_sec 5. Extending To New Closures

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

@section refs_sec 6. Related Pages

- **@subpage 33_Initial_Conditions**
- **@subpage 45_Particle_Initialization_and_Restart**
- **@subpage 26_Walking_Search_Method**
