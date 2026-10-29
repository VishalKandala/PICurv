/**
 * @file ParticleMotion.h
 * @brief Header file for Particle Motion and migration related functions.
 *
 * This file contains declarations of functions responsible for moving and migrating particle swarms within a simulation using PETSc's DMSwarm.
 */

 #ifndef PARTICLE_MOTION_H
 #define PARTICLE_MOTION_H

// Include necessary headers
#include <petsc.h>        // PETSc library header
#include <petscdmswarm.h> // PETSc DMSwarm header
#include <stdbool.h>
#include <petscsys.h>     // For PetscRealloc
#include <math.h>
#include "variables.h"       // Common type definitions
#include "logging.h"      // Logging macros and definitions
#include "walkingsearch.h"  // Walking search function for particle migration   

/**
 * @brief Box-Muller transform of two uniforms on (0,1) into two independent N(0,1) values.
 *
 * @param[in]  u1 First uniform; must be strictly positive.
 * @param[in]  u2 Second uniform.
 * @param[out] n1 First Gaussian number.
 * @param[out] n2 Second Gaussian number.
 */
void GenerateGaussianNoise(PetscReal u1, PetscReal u2, PetscReal *n1, PetscReal *n2);

/**
 * @brief Calculates the stochastic displacement vector (Brownian motion) for a single particle.
 *        Equation: dX_stoch = sqrt(2 * Gamma_eff * dt) * N(0,1)
 *
 * The normals come from PicurvParticleDrawHash() in the Brownian key space with the step
 * as salt, so a particle's increment at a step is a function of `particleRandomSeed`, its
 * PID and the step alone: independent of the rank layout, of the particle's position in
 * the local swarm, and of whether the run was restarted.
 *
 * @param[in]  user          Pointer to UserCtx (access to dt, step and the seed).
 * @param[in]  pid           The particle's global ID.
 * @param[in]  diff_eff      The effective diffusivity (Gamma + Gamma_t) at the particle's location.
 * @param[out] displacement  Pointer to a Cmpnts struct to store the resulting (dx, dy, dz).
 *
 * @return PetscErrorCode
 */
PetscErrorCode CalculateBrownianDisplacement(UserCtx *user, PetscInt64 pid, PetscReal diff_eff, Cmpnts *displacement);

/**
 * @brief Updates a particle's position based on its velocity and the timestep dt (stored in user->dt).
 *
 * @param[in]     user     Pointer to your UserCtx (must contain user->dt).
 * @param[in,out] particle Pointer to the particle struct (contains, pos,vel,diffusivity etc).
 *
 * @return PetscErrorCode  Returns 0 on success, or an error code on failure.
 */
PetscErrorCode UpdateParticlePosition(UserCtx *user, Particle *particle);

/**
 * @brief Loops over all local particles in the DMSwarm, updating their positions
 *        based on velocity and the global timestep user->dt.
 * @param[in,out] user    Pointer to UserCtx (must contain dt).
 *
 * @return PetscErrorCode Returns 0 on success, or an error code on failure.
 */
 PetscErrorCode UpdateAllParticlePositions(UserCtx *user);

/**
 * @brief Checks for particles outside the physical domain boundaries and removes them
 *        using DMSwarmRemovePointAtIndex.
 *
 * This function iterates through all particles local to the current MPI rank.
 * It checks if a particle's position (x, y, or z) is outside the specified
 * physical domain boundaries [xMin, xMax], [yMin, yMax], [zMin, zMax].
 *
 * Particles found out of bounds are flagged in a backward scan and removed together
 * by `RemoveFlaggedParticles`.
 *
 * @param user    Pointer to the UserCtx structure.
 * @param[out] removedCountLocal Pointer to store the number of particles removed *on this rank*.
 * @param[out] removedCountGlobal Pointer to store the total number of particles removed *across all ranks*.
 * @param[in]      bboxlist       An array of BoundingBox structures for ALL MPI ranks, indexed 0 to (size-1).
 *                                This array must be up-to-date and available on all ranks.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode CheckAndRemoveOutOfBoundsParticles(UserCtx *user,
                                              PetscInt *removedCountLocal,
					      PetscInt *removedCountGlobal,
					      const BoundingBox *bboxlist);

/**
 * @brief Removes particles that the location algorithm flagged as LOST or EXITED.
 *
 * This function is the designated cleanup utility. It should be called after the
 * `LocateAllParticlesInGrid` orchestrator has run and every particle's status
 * has been definitively determined.
 *
 * It iterates through all locally owned particles and checks their `DMSwarm_location_status`
 * field. Particles whose status is `LOST` (the search failed) or `EXITED` (the particle left
 * through a removing boundary face) are flagged in a backward scan and permanently removed
 * together by `RemoveFlaggedParticles`. The two kinds are counted and summed separately.
 *
 * This approach centralizes the removal logic, making the `DMSwarm_location_status`
 * the single source of truth for a particle's validity, which is more robust than
 * relying on secondary geometric checks (like bounding boxes).
 *
 * @param[in,out]  user              Pointer to the UserCtx structure containing the swarm.
 * @param[out]     removedCountLocal Pointer to store the number of LOST particles removed on this rank.
 * @param[out]     removedCountGlobal Optional, may be NULL: the number of LOST particles removed across all ranks.
 * @param[out]     removedScalarSumGlobal Optional, may be NULL: the sum of `Psi` over every
 *                 LOST particle removed on any rank. With a 0/1 label and micromixing off, it
 *                 counts the removed particles that carried the label.
 * @param[out]     exitedCountGlobal Optional, may be NULL: the number of EXITED particles removed across all ranks.
 * @param[out]     exitedScalarSumGlobal Optional, may be NULL: the sum of `Psi` over every EXITED particle removed.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 */
PetscErrorCode CheckAndRemoveLostParticles(UserCtx *user,
                                           PetscInt *removedCountLocal,
                                           PetscInt *removedCountGlobal,
                                           PetscReal *removedScalarSumGlobal,
                                           PetscInt *exitedCountGlobal,
                                           PetscReal *exitedScalarSumGlobal);

/**
 * @brief Writes migration destinations into the DMSwarm rank field for marked particles.
 *
 * This helper consumes the migration list built with AddToMigrationList() by the settlement
 * loop (LocateAllParticlesInGrid()) or the restart path (MigrateRestartParticlesUsingCellID())
 * and updates each selected particle's destination rank so that a subsequent
 * `PerformMigration` call can transfer ownership correctly.
 *
 * @param[in,out] user           Context containing the swarm and migration rank field.
 * @param[in]     migrationList  Array of migration directives (local index + destination rank).
 * @param[in]     migrationCount Number of valid entries in `migrationList`.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SetMigrationRanks(UserCtx* user, const MigrationInfo *migrationList, PetscInt migrationCount);

/**
 * @brief Performs particle migration based on the pre-populated DMSwarmPICField_rank field.
 *
 * Assumes SetMigrationRanks has already been called to mark particles with their target ranks.
 * Calls DMSwarmMigrate to execute the communication while keeping the sent particles, then
 * removes them in one pass with `RemoveFlaggedParticles`. The resulting local order is the
 * one DMSwarmMigrate's own removal produces: the particles that stayed occupy the first
 * `*n_survivors` entries and the received ones follow them.
 *
 * @param user Pointer to the UserCtx structure containing the swarm.
 * @param[out] n_survivors Number of local particles that stayed, or NULL when not wanted.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode PerformMigration(UserCtx *user, PetscInt *n_survivors);

/**
 * @brief Removes flagged swarm points in one pass, preserving the one-at-a-time removal order.
 *
 * Each flagged slot in `[0, n_scan)` is filled with the last kept point of that range,
 * exactly as `DMSwarmRemovePointAtIndex` would when the points are removed in scan order,
 * so the survivors end in the same order. Points after `n_scan` are kept and moved down
 * behind the survivors. The swarm is resized once, instead of once per removed point,
 * since each resize re-zeroes the whole swarm buffer.
 *
 * @param[in,out] swarm         The DMSwarm. No field may be held (DMSwarmGetField) on entry.
 * @param[in]     n_scan        Number of leading local points the flags cover.
 * @param[in,out] remove        Flags for points `[0, n_scan)`; overwritten as points move.
 * @param[in]     scan_backward PETSC_TRUE for removal in descending index order,
 *                              PETSC_FALSE for ascending order (DMSwarmMigrate's order).
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode RemoveFlaggedParticles(DM swarm, PetscInt n_scan, PetscBool remove[], PetscBool scan_backward);

/**
 * @brief Reorders the local swarm points so that points in the same host cell are contiguous.
 *
 * Points are ordered by their `DMSwarm_CellID` with k slowest and i fastest, stably, so points
 * sharing a cell keep their relative order; points without a valid cell go last. Every field
 * in the particle field catalog (which lists every field the swarm registers, PETSc's pid and
 * rank included) is permuted together, so each point keeps all of its values. Interpolation,
 * location and the other per-point loops then read neighbouring grid data for consecutive
 * points. Per-point results do not depend on the order; sums over points (the Eulerian
 * scatter, and the IEM cell means built from it) are taken in the new order and can differ
 * in rounding.
 *
 * @param[in,out] swarm The DMSwarm. No field may be held (DMSwarmGetField) on entry.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode SortParticlesByCell(DM swarm);

/**
 * @brief Counts particles in each cell of the DMDA 'da' and stores the result in user->ParticleCount.
 *
 * Assumes user->ParticleCount is a pre-allocated global vector associated with user->da
 * and initialized to zero before calling this function (though it resets it internally).
 * Assumes particle 'DMSwarm_CellID' field contains local cell indices.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing da, swarm, and ParticleCount.
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode CalculateParticleCountPerCell(UserCtx *user);

/**
 * @brief Resizes a swarm collectively to a target global particle count.
 *
 * The target is divided by quotient and remainder so every rank receives either
 * `floor(N_target / nranks)` or one additional entry. Resizing establishes the
 * storage layout; callers initialize or overwrite particle fields afterwards.
 *
 * @param[in,out] swarm    Swarm object to resize.
 * @param[in]     N_target Target global particle count.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ResizeSwarmGlobally(DM swarm, PetscInt N_target);

/**
 * @brief Checks particle count in the reference file and resizes the swarm if needed.
 *
 * Reads the specified field file (e.g., position) into a temporary Vec to determine
 * the number of particles (`N_file`) represented in that file for the given timestep.
 * Compares `N_file` with the current swarm size (`N_current`). If they differ,
 * resizes the swarm globally (adds or removes particles) to match `N_file`.
 * The resized population is balanced across the communicator before field input.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing the DMSwarm.
 * @param[in]     ti   Time index for constructing the file name.
 * @param[in]     ext  File extension (e.g., "dat").
 *
 * @return PetscErrorCode 0 on success, non-zero on critical failure.
 */
PetscErrorCode PreCheckAndResizeSwarm(UserCtx *user, PetscInt ti, const char *ext);

/**
 * @brief Creates a sorted snapshot of all Particle IDs (PIDs) from a raw data array.
 * @ingroup ParticleUtils
 *
 * This function is a crucial helper for the migration process. It captures the state of
 * which particles are on the current MPI rank *before* migration occurs by taking a
 * pointer to the swarm's raw PID data array. The resulting sorted array can then be used
 * with an efficient binary search to quickly identify newcomer particles after migration.
 *
 * This function does NOT call DMSwarmGetField/RestoreField. It is the caller's
 * responsibility to acquire the `pid_field` pointer before calling and restore it afterward.
 *
 * @param[in]  pid_field         A read-only pointer to the raw array of PIDs for the local swarm.
 * @param[in]  n_local           The number of particles currently on the local rank.
 * @param[out] pids_snapshot_out A pointer to a `PetscInt64*` array. This function will
 *                               allocate memory for this array, and the caller is
 *                               responsible for freeing it with `PetscFree()` when it
 *                               is no longer needed.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 */
PetscErrorCode GetLocalPIDSnapshot(const PetscInt64 pid_field[], 
                                   PetscInt n_local, 
                                   PetscInt64 **pids_snapshot_out);

/**
 * @brief Safely adds a new migration task to a dynamically sized list.
 *
 * This utility function manages a dynamic array of MigrationInfo structs. It appends
 * a new entry to the list and automatically doubles the array's capacity using
 * `PetscRealloc` if the current capacity is exceeded. This prevents buffer overflows
 * and avoids the need to know the number of migrating particles in advance.
 *
 * @param[in,out] migration_list_p  A pointer to the MigrationInfo array pointer. The function
 *                                  will update this pointer if the array is reallocated.
 * @param[in,out] capacity_p        A pointer to an integer holding the current allocated
 *                                  capacity of the list (in number of elements). This will be
 *                                  updated upon reallocation.
 * @param[in,out] count_p           A pointer to an integer holding the current number of
 *                                  items in the list. This will be incremented by one.
 * @param[in]     particle_local_idx The local index (from 0 to nlocal-1) of the particle
 *                                  that needs to be migrated.
 * @param[in]     destination_rank   The target MPI rank for the particle.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure (e.g., from memory allocation).
 */
PetscErrorCode AddToMigrationList(MigrationInfo **migration_list_p,
                                  PetscInt *capacity_p,
                                  PetscInt *count_p,
                                  PetscInt particle_local_idx,
                                  PetscMPIInt destination_rank);

/**
 * @brief Identifies newly arrived particles after migration and flags them for a location search.
 * @ingroup ParticleMotion
 *
 * This function is a critical component of the iterative migration process managed by
 * the main particle settlement orchestrator (e.g., `SettleParticles`). After a
 * `DMSwarmMigrate` call, each rank's local particle list is a new mix of resident
 * particles and newly received ones. This function's job is to efficiently identify
 * these "newcomers" and set their `DMSwarm_location_status` field to `NEEDS_LOCATION`.
 *
 * This ensures that in the subsequent pass of the migration `do-while` loop, only the
 * newly arrived particles are processed by the expensive location algorithm, preventing
 * redundant work on particles that are already settled on the current rank.
 *
 * PerformMigration() leaves the particles that stayed first and the received ones after
 * them, so the newcomers are exactly the entries from `first_newcomer` on: no PID lookup
 * is needed, whatever order the particles are stored in.
 *
 * @param[in] swarm          The DMSwarm object, which has just completed a migration.
 * @param[in] n_local_before Number of local particles before the migration, for the log.
 * @param[in] first_newcomer Index of the first received particle: the survivor count
 *                           PerformMigration() reported.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 */
PetscErrorCode FlagNewcomersForLocation(DM swarm, PetscInt n_local_before, PetscInt first_newcomer);


/**
 * @brief Fast-path migration for restart particles using preloaded Cell IDs.
 *
 * This function provides an optimized migration path for particles whose Cell ID is
 * already known: particles loaded from restart files, and freshly seeded particles,
 * whose seeding records the cell it placed them in (a point source's cell may belong to
 * another rank). Unlike the standard `LocateAllParticlesInGrid()`, which performs
 * walking searches, this function sends each such particle straight to its cell's owner.
 * Particles without a Cell ID (-1) are left for `LocateAllParticlesInGrid()`.
 *
 * **How It Works:**
 * 1. Iterates through all local particles.
 * 2. For each particle with a valid Cell ID (ci, cj, ck):
 *    - Calls `FindOwnerOfCell(ci, cj, ck)` to determine the correct rank.
 *    - If owner differs from current rank, adds to migration list.
 *    - If owner matches current rank, the existing `ACTIVE_AND_LOCATED` status is preserved.
 * 3. Uses existing `SetMigrationRanks()` and `PerformMigration()` infrastructure.
 * 4. Achieves **single-pass direct migration** (no multi-hop, no walking searches).
 *
 * @param[in,out] user Pointer to UserCtx containing the swarm and RankCellInfoMap.
 *                     The function updates particle status fields and performs migration.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 *
 * @note Testing status:
 *       Direct coverage currently focuses on restart fast-path ownership transfer.
 *       Non-restart multi-pass migration behavior remains part of the next
 *       simulation-core test backlog.
 */
PetscErrorCode MigrateRestartParticlesUsingCellID(UserCtx *user);

/**
 * @brief Orchestrates the complete particle location and migration process for one timestep.
 * @ingroup ParticleLocation
 *
 * This function is the master orchestrator for ensuring every particle is on its correct
 * MPI rank and has a valid host cell index. It is designed to be called once per
 * timestep after particle positions have been updated.
 *
 * The function uses a robust, iterative "Guess and Verify" strategy within a
 * do-while loop to handle complex particle motion across processor boundaries,
 * especially on curvilinear grids.
 *
 * 1.  **State Snapshot:** At the start of each pass, it captures a list of all Particle IDs (PIDs)
 *     on the current rank.
 * 2.  **"Guess" (Heuristic):** For particles that are "lost" (no valid host cell),
 *     it first attempts a fast, bounding-box-based guess to find a potential new owner rank.
 * 3.  **"Verify" (Robust Walk):** For all other particles, or if the guess fails,
 *     it uses a robust cell-walking algorithm (`LocateParticleOrFindMigrationTarget`)
 *     that determines the particle's status: located locally, needs migration, exited the
 *     domain, or lost. A walk from the prior cell that the domain boundary blocked is retried
 *     through the guess path in the same pass.
 *     `LOST` and `EXITED` are final within one call: only pass 1 searches them.
 * 4.  **Migration:** After identifying all migrating particles on a pass, it performs the
 *     MPI communication using the `SetMigrationRanks` and `PerformMigration` helpers.
 * 5.  **Newcomer Flagging:** After migration, it uses the PID snapshot from step 1 to
 *     efficiently identify newly arrived particles and flag them for location on the next pass.
 * 6.  **Iteration:** The process repeats in a `do-while` loop until a pass occurs where
 *     no particles migrate, ensuring the entire swarm is in a stable, consistent state.
 *
 * @param[in,out] user Pointer to the UserCtx, containing the swarm and all necessary
 *                     domain topology information (bboxlist, RankCellInfoMap, etc.).
 * @param[in] bboxlist  An array of BoundingBox structures for ALL MPI ranks, indexed 0 to (size-1).
 *                      This array must be up-to-date and available on all ranks.
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 *
 * @note Testing status:
 *       Direct unit coverage currently pins the prior-cell fast path, the
 *       local guess-then-verify path, and (two ranks) an exit that stays
 *       classified once across a migration pass. Newcomer flagging and several
 *       lost/migration edge cases are still targeted for future bespoke tests.
 */
PetscErrorCode LocateAllParticlesInGrid(UserCtx *user,BoundingBox *bboxlist);

/**
 * @brief Marks all local particles as `NEEDS_LOCATION` for the next settlement pass.
 *
 * This function is designed to be called at the end of a full timestep, after all
 * particle-based calculations are complete. It prepares the swarm for the next
 * timestep by ensuring that after the next position update, every particle will be
 * re-evaluated by the LocateAllParticlesInGrid orchestrator.
 *
 * It iterates through all locally owned particles and sets their
 * `DMSwarm_location_status` field to `NEEDS_LOCATION`.
 *
 * @param[in,out] user Pointer to the UserCtx containing the swarm.
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 */
PetscErrorCode ResetAllParticleStatuses(UserCtx *user);

 #endif // PARTICLE_MOTION_H
