@page 26_Walking_Search_Method Walking Search for Particle Location

@anchor _Walking_Search_Method

PICurv uses a robust locate-and-migrate pipeline so every particle ends each step with a valid owner rank and host cell, or is explicitly marked as having left the domain (`EXITED`) or as a search failure (`LOST`).

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
- `LOST` - the search failed for a particle believed inside the domain
- `UNINITIALIZED`
- `EXITED` - the particle left through a face whose policy is `remove`

These statuses drive the iterative migration passes in the orchestrator. `LOST` and `EXITED`
are final within one settlement: only pass 1 searches them, so a particle beyond the
boundary is never re-walked and re-migrated pass after pass. Cleanup removes both and
counts them separately (@ref p53_definitions_legacy).

@section p26_boundary_sec 2a. Domain Boundaries

When a walk ends outside a boundary cell, the action depends on the crossed faces'
`particle_policy` (@ref p44_ppol_sec):

- any crossed face periodic (`wrap`): the particle is shifted by every crossed seam's
  translation at once and the walk moves to the cell at the other end of those axes, handing
  off to its owner if that cell is remote, so an edge or corner crossing of several seams needs
  one handoff; a wall crossed in the same corner reflects on the next evaluation;
- otherwise, on a pass-1 walk that started from the particle's previous cell, an end more than three
  cell sizes beyond the face is treated as blocked by the domain's shape (a non-convex boundary
  between the start and the target), and the particle is searched again from scratch in the
  same pass;
- any crossed face `remove`: the particle is `EXITED`;
- all crossed faces `reflect`: the particle is mirrored about the most-crossed face plane and
  the cell is evaluated again, up to four reflections per search, after which it is placed
  just inside the boundary cell (counted as a reflection fallback). Only a cell warped so far
  that its centroid is outside the crossed face planes cannot take the fallback, and the
  particle is then `LOST`.

The rank owning a cell is found by a binary search of the decomposition's per-axis range
starts (@ref FindOwnerOfCell), so its cost does not grow with the rank count.

@section p26_implementation_sec 3. Implementation Touchpoints

- single-particle walking search: @ref LocateParticleInGrid
- robust locate-or-migrate decision: `LocateParticleOrFindMigrationTarget` (walking-search module)
- global orchestrator: @ref LocateAllParticlesInGrid
- rank migration communication: @ref PerformMigration (sent particles are removed afterwards in one pass by @ref RemoveFlaggedParticles, keeping the one-at-a-time survivor order; measurement `swarm-one-pass-removal-2026-10-20`)
- owner of a cell: @ref FindOwnerOfCell
- wall reflection: @ref ReflectParticleAcrossCellFace, with the cap fallback @ref PlaceParticleInsideCell
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

Direct coverage:

- `LocateAllParticlesInGrid`: the prior-cell fast path and the local guess-then-verify path
- restart fast-path migration through `MigrateRestartParticlesUsingCellID`
- `LocateParticleOrFindMigrationTarget` boundary branches: reflection about a face
  (`make unit-particles`), an exit classified once across settlement passes on two ranks
  (`make unit-mpi`), and wraps through a face, an x-y corner and an off-axis seam
  (`make unit-periodic`)
- `FindOwnerOfCell` against a scan of every rank's ranges (`make unit-mpi`)

Remaining direct gaps:

- branch coverage of `LocateParticleOrFindMigrationTarget` for the boundary clamp, ghost-region
  handoff, the stuck tie-breaker and the blocked-walk retry
- direction-complete and failure-path coverage for the `GuessParticleOwnerWithBBox` heuristic
- non-restart MPI migration tests for multi-pass handoff and newcomer flagging

@section p26_practical_sec 6. Failure Modes and Signals

Watch for:

- excessive migration pass counts,
- repeated `LOST` particles (`EXITED` through an outlet is outflow, not a failure),
- nonzero `reflection_fallbacks`, or many `blocked_walk_retry_count` events outside a strongly non-convex domain,
- non-convergence of settlement loop,
- elevated traversal effort or repeated tie-break/boundary-clamp events in `<run.analysis.metrics>/search_metrics.csv`,
- increasing `search_failure_fraction`, `search_work_index`, or `re_search_fraction`,
- increasing `lost` or `lost_cumulative` even when migration otherwise appears healthy.

For precise metric definitions and formulas, see **@subpage 53_Search_Robustness_Metrics_Reference**.

These usually indicate domain-decomposition mismatch, bounding-box mismatch, or invalid particle coordinates.

See **@subpage 39_Common_Fatal_Errors** for troubleshooting commands.
