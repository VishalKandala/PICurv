@page 26_Walking_Search_Method Walking Search for Particle Location

@anchor _Walking_Search_Method

PICurv uses a robust locate-and-migrate pipeline so every particle ends each step with a valid owner rank and host cell (or is explicitly marked lost).

@tableofcontents

@section p26_concept_sec 1. Core Search Concept

A particle is tested against candidate cell faces using signed distances.
Search starts from prior cell information when available, then walks to neighboring cells by face-crossing logic until a terminal status is reached.

This avoids expensive global brute-force lookup and is suitable for distributed curvilinear grids.

@section p26_statuses_sec 2. Settlement Status Model

Location/migration logic uses `ParticleLocationStatus`:

- `NEEDS_LOCATION`
- `ACTIVE_AND_LOCATED`
- `MIGRATING_OUT`
- `LOST`
- `UNINITIALIZED`

These statuses drive the iterative migration passes in the orchestrator.

@section p26_implementation_sec 3. Implementation Touchpoints

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

@section p26_restart_sec 4. Restart Path Specifics

For restart-loaded particles with valid `CellID`, PICurv can skip full walking search and directly resolve owner rank via cell ownership map before migration.
This is usually much cheaper than full re-location on first restart step.

@section p26_testing_sec 5. Current test status

Current direct coverage is strongest at the orchestrator level:

- `LocateAllParticlesInGrid` has dedicated direct tests for the prior-cell fast path and the local guess-then-verify path
- restart fast-path migration is covered directly through `MigrateRestartParticlesUsingCellID`

Current direct gaps still targeted for the next remediation are:

- branch-complete unit coverage for `LocateParticleOrFindMigrationTarget`
  - boundary clamp
  - ghost-region handoff
  - stuck/tie-breaker
  - `LOST` and `MIGRATING_OUT`
- direction-complete and failure-path coverage for the `GuessParticleOwnerWithBBox` heuristic
- non-restart MPI migration tests for multi-pass handoff and newcomer flagging

@section p26_practical_sec 6. Failure Modes and Signals

Watch for:

- excessive migration pass counts,
- repeated `LOST` particles,
- non-convergence of settlement loop,
- elevated traversal effort or repeated tie-break/boundary-clamp events in `logs/search_metrics.csv`,
- increasing `lost` or `lost_cumulative` even when migration otherwise appears healthy.

These usually indicate domain-decomposition mismatch, bounding-box mismatch, or invalid particle coordinates.

See **@subpage 39_Common_Fatal_Errors** for troubleshooting commands.

<!-- DOC_EXPANSION_CFD_GUIDANCE -->

## CFD Reader Guidance and Practical Use

This page describes **Walking Search for Particle Location** within the PICurv workflow. For CFD users, the most reliable reading strategy is to map the page content to a concrete run decision: what is configured, what runtime stage it influences, and which diagnostics should confirm expected behavior.

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
