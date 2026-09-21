@page 34_Particle_Model_Overview Particle Model and Coupling Overview

@anchor _Particle_Model_Overview

This page documents the particle pipeline exactly as orchestrated in the current solver flow.

@tableofcontents

@dotfile particle_step.dot "Order of operations in one coupled particle step"

@note Diagram derived from the implementation in `src/runloop.c`; the advection terms
are those in @ref UpdateParticlePosition. Flagged for owner validation of the
scientific interpretation.


@section p34_loop_sec 1. Per-Step Particle Pipeline

Order inside the coupled step, as implemented in @ref AdvanceSimulation
(`src/runloop.c`):

1. refresh Eulerian transport properties when a turbulence model is active:
   @ref ComputeEulerianDiffusivity and @ref ComputeEulerianDiffusivityGradient
2. advect particle positions using the velocity interpolated at the **previous**
   step: @ref UpdateAllParticlePositions
3. settle particles - locate new host cells and migrate across ranks:
   @ref LocateAllParticlesInGrid
4. remove particles that are now lost or outside the global domain:
   @ref CheckAndRemoveLostParticles
5. interpolate the newly computed Eulerian velocity onto the settled positions:
   @ref InterpolateAllFieldsToSwarm
6. update particle physics fields: @ref UpdateAllParticleFields
7. scatter particle fields back to Eulerian storage:
   @ref CalculateParticleCountPerCell and @ref ScatterAllParticleFieldsToEulerFields

The advect-then-interpolate ordering is deliberate and is the single most
important property of this loop. A particle carries the state interpolated at
\f$t_n\f$ into step 2, and @ref UpdateParticlePosition displaces it by three
contributions - convective velocity, drift along the diffusivity gradient, and a
stochastic Brownian kick:

\f[
\mathbf{X}^{n+1} = \mathbf{X}^{n}
  + \left( \mathbf{u}_p^{\,n} + \nabla D^{\,n} \right) \Delta t
  + \Delta \mathbf{X}_{\mathrm{brownian}},
\f]

where \f$\nabla D\f$ is the interpolated diffusivity gradient carried on the
swarm and \f$\Delta\mathbf{X}_{\mathrm{brownian}}\f$ comes from
@ref CalculateBrownianDisplacement, scaled by the particle's interpolated
diffusivity. Step 5 then supplies the state at \f$t_{n+1}\f$ from the flow field
@ref FlowSolver has just computed, ready for the next step. Interpolation
therefore closes each step rather than opening it.

Particles are not advected without a valid interpolated state: both
@ref PerformInitializedParticleSetup and @ref PerformLoadedParticleSetup settle the
swarm and run one interpolation at \f$t=0\f$ before the step loop begins, so the
first physical advection already has interpolated velocity, diffusivity, and
diffusivity gradient available.

Each stage relies on valid `DMSwarm_CellID`, interpolation weights, and synchronized ghost data.

Step 5 dispatches through @ref InterpolateEulerFieldToSwarm, which selects between two interpolation paths based on `SimCtx.interpolationMethod`:

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

This presently implements IEM-style relaxation for `Psi`, with the mixing constant fixed at
`C_IEM = 2.0` inside @ref UpdateFieldForAllParticles; it is not configurable. It is also
inert in practice: `Psi` starts at zero on every particle and no configuration seeds or
sources a scalar, so the relaxation has nothing to act on (@ref p28_iem_sec).

@section p34_statistics_sec 4. Statistics and Diagnostics

Post statistics currently include global kernels such as:

- @ref ComputeParticleMSD

Additional health indicators are available from migration counters and settlement-pass counts stored in `SimCtx` fields updated by location logic.
For particle-enabled runs, the runtime also writes `<run.analysis.metrics>/search_metrics.csv`
with per-step and cumulative loss plus aggregated search-attempt, traversal,
tie-break, boundary-clamp, bbox-guess, and pass-depth signals. The
`examples/search_robustness/` bundle is the dedicated end-to-end reference for
interpreting that artifact.

@section p34_limits_sec 5. Limitations

- **Particles are not wrapped across periodic boundaries.** Neither the motion update
  nor the walking search handles periodicity, so a particle that crosses a periodic face
  leaves the domain and is removed as lost.
- **A lost particle is removed, not reflected.** Steady loss through a boundary you did
  not intend usually means the carrier velocity or the diffusion step carries particles
  out; the per-step lost count reports it.
- **Diffusion, advection, the gradient drift and interpolation are verified.** A
  10,000-particle cloud in zero flow reproduced the Einstein relation to 0.17% in the MSD
  slope; a cloud in uniform flow drifted at exactly the carrier velocity; paired runs with
  the same seed resolved the drift a linear diffusivity gradient adds to 1.1% of `a t`
  (`diffusivity-gradient-drift-paired-2026-09-18`); and Trilinear interpolation of the
  TGV3D field converged at order 1.97. Only one gradient, and zero or uniform carriers,
  were measured. The examples are listed in @ref p65_verify_sec.
- **Brownian trajectories repeat only for the same rank count.** Placement and Brownian
  draws are seeded from `models.physics.particles.random_seed` (default `12345`) plus the
  MPI rank, so identical inputs on the same number of ranks give identical trajectories,
  while a different rank count gives a different, equally valid, realization. A restart
  does not continue the sequence it stopped in: it draws a fresh sequence determined by
  the seed and the restart step.

@section p34_extension_sec 6. Extending To New Closures

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

@section p34_refs_sec 7. Related Pages

- **@subpage 33_Initial_Conditions**
- **@subpage 45_Particle_Initialization_and_Restart**
- **@subpage 26_Walking_Search_Method**
