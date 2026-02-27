@page 34_Particle_Model_Overview Particle Model and Coupling Overview

This page maps the particle pipeline: motion, location, coupling, physics update, and statistics.

@tableofcontents

@section pipeline_sec 1. Runtime Pipeline

Typical per-step flow:
1. Interpolate Eulerian fields to particles: @ref InterpolateAllFieldsToSwarm
2. Move particles: @ref UpdateAllParticlePositions
3. Locate/migrate particles: @ref LocateAllParticlesInGrid
4. Update particle physics fields: @ref UpdateAllParticleFields
5. Scatter particle data back to Eulerian fields: @ref ScatterAllParticleFieldsToEulerFields

@section physics_sec 2. Current Particle Physics

Current scalar-mixing update is implemented through:
- @ref UpdateParticleField
- @ref UpdateFieldForAllParticles
- @ref UpdateAllParticleFields

The present built-in model is IEM-style relaxation for `Psi`.

@section stats_sec 3. Statistics Pipeline

Post statistics kernel entry currently includes:
- @ref ComputeParticleMSD

Statistics are emitted as CSV outputs from postprocessing paths.

@section extension_sec 4. Data-Driven Closure Extension Path

For future data-driven closures:
1. Add a model selector in schema + control generation.
2. Add runtime dispatch in particle physics update path.
3. Keep feature extraction and model I/O isolated from motion/location kernels.
4. For tightly coupled inference, treat model invocation as a pluggable kernel behind the same `UpdateFieldForAllParticles` entry.
