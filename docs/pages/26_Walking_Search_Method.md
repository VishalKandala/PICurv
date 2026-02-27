@page 26_Walking_Search_Method Walking Search for Particle Location

PICurv uses a robust locate-and-migrate pipeline so every particle ends each step with a valid owner rank and host cell (or is explicitly marked lost).

@tableofcontents

@section concept_sec 1. Core Search Concept

A particle is tested against candidate cell faces using signed distances.
Search starts from prior cell information when available, then walks to neighboring cells by face-crossing logic until a terminal status is reached.

This avoids expensive global brute-force lookup and is suitable for distributed curvilinear grids.

@section statuses_sec 2. Settlement Status Model

Location/migration logic uses `ParticleLocationStatus`:

- `NEEDS_LOCATION`
- `ACTIVE_AND_LOCATED`
- `MIGRATING_OUT`
- `LOST`
- `UNINITIALIZED`

These statuses drive the iterative migration passes in the orchestrator.

@section implementation_sec 3. Implementation Touchpoints

- single-particle walking search: @ref LocateParticleInGrid
- robust locate-or-migrate decision: `LocateParticleOrFindMigrationTarget` (walking-search module)
- global orchestrator: @ref LocateAllParticlesInGrid
- rank migration communication: @ref PerformMigration
- per-particle migration cycle helper: @ref PerformSingleParticleMigrationCycle
- restart fast-path migration: @ref MigrateRestartParticlesUsingCellID

The current orchestrator includes:

1. PID snapshot,
2. per-particle locate/guess/verify,
3. migration,
4. newcomer flagging,
5. repeat until no global migrations remain (or safety pass cap).

@section restart_sec 4. Restart Path Specifics

For restart-loaded particles with valid `CellID`, PICurv can skip full walking search and directly resolve owner rank via cell ownership map before migration.
This is usually much cheaper than full re-location on first restart step.

@section practical_sec 5. Failure Modes and Signals

Watch for:

- excessive migration pass counts,
- repeated `LOST` particles,
- non-convergence of settlement loop.

These usually indicate domain-decomposition mismatch, bounding-box mismatch, or invalid particle coordinates.

See **@subpage 39_Common_Fatal_Errors** for troubleshooting commands.
