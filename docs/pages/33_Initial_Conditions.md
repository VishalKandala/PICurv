@page 33_Initial_Conditions Initial Condition Modes

This page documents the currently implemented Eulerian and particle initialization paths, including their CLI mappings.

@tableofcontents

@section euler_sec 1. Eulerian Initialization

Main initialization entry:

- @ref InitializeEulerianState

Interior field initializer:

- @ref SetInitialInteriorField

Current `-finit` mappings:

- `0`: zero velocity
- `1`: constant normal velocity
- `2`: Poiseuille-like normal profile

Human-readable mapping helper:

- @ref FieldInitializationToString

For the Poiseuille-like profile, code uses index-space normalized coordinates and applies a separable parabolic form in cross-stream directions.

@section euler_formula_sec 2. Contravariant Initialization Note

Initialization is applied to contravariant velocity components, scaled by face-area metrics:

\f[
U_n = v_n\,A_n,
\f]

where face areas are derived from metric vectors (`Csi`, `Eta`, `Zet`) and sign is aligned with inlet face orientation.

This is why equivalent physical inflow speed can map to different `Ucont` magnitudes on stretched/curved meshes.

@section particle_sec 3. Particle Initialization Modes

Particle swarm entry:

- @ref InitializeParticleSwarm

`-pinit` / `ParticleInitializationType` modes:

- `0`: `PARTICLE_INIT_SURFACE_RANDOM`
- `1`: `PARTICLE_INIT_VOLUME`
- `2`: `PARTICLE_INIT_POINT_SOURCE`
- `3`: `PARTICLE_INIT_SURFACE_EDGES`

Human-readable mapping helper:

- @ref ParticleInitializationToString

Inlet refresh helper:

- @ref ReinitializeParticlesOnInletSurface

@section restart_sec 4. Restart Coupling

`-particle_restart_mode` currently accepts:

- `load`
- `init`

For restart files with valid cell IDs, migration can use fast ownership-based path via @ref MigrateRestartParticlesUsingCellID before full walking-search fallback.

@section refs_sec 5. Related Pages

- **@subpage 32_Analytical_Solutions**
- **@subpage 34_Particle_Model_Overview**
- **@subpage 39_Common_Fatal_Errors**
