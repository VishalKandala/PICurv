@page 33_Initial_Conditions Initial Condition Modes

This page summarizes currently available Eulerian and particle initial-condition paths.

@tableofcontents

@section euler_ic_sec 1. Eulerian Initial Conditions

Main entry point:
- @ref InitializeEulerianState

Interior field initializer:
- @ref SetInitialInteriorField

`-finit` values currently used:
- `0`: zero velocity
- `1`: constant normal velocity
- `2`: Poiseuille normal velocity

Human-readable mapping:
- @ref FieldInitializationToString

@section particle_ic_sec 2. Particle Initialization Modes

Particle swarm entry point:
- @ref InitializeParticleSwarm

Mode enum (`-pinit`) in `ParticleInitializationType`:
- `0`: `PARTICLE_INIT_SURFACE_RANDOM`
- `1`: `PARTICLE_INIT_VOLUME`
- `2`: `PARTICLE_INIT_POINT_SOURCE`
- `3`: `PARTICLE_INIT_SURFACE_EDGES`

Human-readable mapping:
- @ref ParticleInitializationToString

Inlet-surface refresh helper:
- @ref ReinitializeParticlesOnInletSurface

@section restart_sec 3. Restart Coupling

`-particle_restart_mode` currently accepts:
- `load`
- `init`

Final restart coupling is handled in simulation orchestration paths after fluid state initialization.
