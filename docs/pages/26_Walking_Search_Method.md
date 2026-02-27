@page 26_Walking_Search_Method Walking Search for Particle Location

Particle tracking requires robustly determining each particle’s host cell on distributed curvilinear grids. PICurv uses a walking-search based location strategy with migration support.

@tableofcontents

@section purpose_sec 1. Why Walking Search

- Particles move continuously while mesh ownership is discrete by rank/block/cell.
- A local walk from prior host-cell state is more efficient than global brute-force scans.
- The method supports migration decisions when particles leave local ownership.

@section code_sec 2. Code Touchpoints

- Single-particle locator: @ref LocateParticleInGrid
- Multi-particle settlement orchestrator: @ref LocateAllParticlesInGrid
- Migration helpers in motion layer: @ref PerformMigration, @ref PerformSingleParticleMigrationCycle
- Restart-cell migration path: @ref MigrateRestartParticlesUsingCellID

@section workflow_sec 3. Typical Flow

1. Predict new particle position.
2. Attempt local walking search from previous cell.
3. If local ownership fails, determine migration target rank/block.
4. Migrate, then re-locate and refresh interpolation weights.

@section refs_sec 4. References

- PETSc DMSwarm background: https://petsc.org/release/manualpages/DMSWARM/
