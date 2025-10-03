/**
 * @file ParticleSwarm.h
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
 * @brief Updates a particle's position based on its velocity and the timestep dt (stored in user->dt).
 *
 * @param[in]     user     Pointer to your UserCtx (must contain user->dt).
 * @param[in,out] position Pointer to the particle's current position (Cmpnts).
 * @param[in]     velocity Pointer to the particle's velocity (Cmpnts).
 *
 * @return PetscErrorCode  Returns 0 on success, or an error code on failure.
 */
 PetscErrorCode UpdateParticlePosition(UserCtx *user, Cmpnts *position, const Cmpnts *velocity);

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
 * If a particle is found out of bounds, it is removed using DMSwarmRemovePointAtIndex.
 * NOTE: Removing points changes the indices of subsequent points in the iteration.
 *       Therefore, it's crucial to iterate BACKWARDS or carefully manage indices
 *       after a removal. Iterating backwards is generally safer.
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
 * @brief Removes particles that have been definitively flagged as LOST by the location algorithm.
 *
 * This function is the designated cleanup utility. It should be called after the
 * `LocateAllParticlesInGrid` orchestrator has run and every particle's status
 * has been definitively determined.
 *
 * It iterates through all locally owned particles and checks their `DMSwarm_location_status`
 * field. If a particle's status is `LOST`, it is permanently removed from the simulation
 * using `DMSwarmRemovePointAtIndex`.
 *
 * This approach centralizes the removal logic, making the `DMSwarm_location_status`
 * the single source of truth for a particle's validity, which is more robust than
 * relying on secondary geometric checks (like bounding boxes).
 *
 * @param[in,out]  user              Pointer to the UserCtx structure containing the swarm.
 * @param[out]     removedCountLocal Pointer to store the number of particles removed on this rank.
 * @param[out]     removedCountGlobal Pointer to store the total number of particles removed across all ranks.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 */
PetscErrorCode CheckAndRemoveLostParticles(UserCtx *user,
                                           PetscInt *removedCountLocal,
                                           PetscInt *removedCountGlobal);

/**
 * @brief Defines the basic migration pattern for particles within the swarm.
 *
 * This function establishes the migration pattern that dictates how particles
 * move between different MPI ranks in the simulation. It initializes a migration
 * list where each particle is assigned a target rank based on predefined conditions.
 * The migration pattern can be customized to implement various migration behaviors.
 *
 * @param[in,out] user    Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode DefineBasicMigrationPattern(UserCtx* user);

/**
 * @brief Performs the basic migration of particles based on the defined migration pattern.
 *
 * This function updates the positions of particles within the swarm by migrating them
 * to target MPI ranks as specified in the migration list. It handles the migration process
 * by setting the 'DMSwarm_rank' field for each particle and invokes the DMSwarm migration
 * mechanism to relocate particles across MPI processes. After migration, it cleans up
 * allocated resources and ensures synchronization across all MPI ranks.
 *
 * @param[in,out] user    Pointer to the UserCtx structure containing simulation context.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode PerformBasicMigration(UserCtx* user);

/**
 * @brief Identifies particles leaving the local bounding box and finds their target neighbor rank.
 *
 * Iterates local particles, checks against local bounding box. If outside, checks
 * the pre-computed immediate neighbors (user->neighbors) using the global bboxlist
 * to see if the particle landed in one of them. Populates the migrationList.
 * Does NOT handle particles leaving the global domain (assumes CheckAndRemove was called).
 *
 * @param user           Pointer to the UserCtx (contains local bbox and neighbors).
 * @param bboxlist       Array of BoundingBox structs for all ranks (for checking neighbor boxes).
 * @param migrationList  Pointer to an array of MigrationInfo structs (output, allocated/reallocated by this func).
 * @param migrationCount Pointer to the number of particles marked for migration (output).
 * @param listCapacity   Pointer to the current allocated capacity of migrationList (in/out).
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode IdentifyMigratingParticles(UserCtx *user,
                                        const BoundingBox *bboxlist,
                                        MigrationInfo **migrationList,
                                        PetscInt *migrationCount,
					  PetscInt *listCapacity);

// --- Helper function to set migration rank field ---
// This needs to be called AFTER identifying migrating particles
PetscErrorCode SetMigrationRanks(UserCtx* user, const MigrationInfo *migrationList, PetscInt migrationCount);

/**
 * @brief Performs particle migration based on the pre-populated DMSwarmPICField_rank field.
 *
 * Assumes SetMigrationRanks has already been called to mark particles with their target ranks.
 * Calls DMSwarmMigrate to execute the communication and removal of un-migrated particles.
 *
 * @param user Pointer to the UserCtx structure containing the swarm.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode PerformMigration(UserCtx *user);

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

// --- Helper function to resize swarm globally (add or remove) ---
// This assumes removing excess particles means removing the globally last ones.
PetscErrorCode ResizeSwarmGlobally(DM swarm, PetscInt N_target);

/**
 * @brief Checks particle count in the reference file and resizes the swarm if needed.
 *
 * Reads the specified field file (e.g., position) into a temporary Vec to determine
 * the number of particles (`N_file`) represented in that file for the given timestep.
 * Compares `N_file` with the current swarm size (`N_current`). If they differ,
 * resizes the swarm globally (adds or removes particles) to match `N_file`.
 * Removal assumes excess particles are the globally last ones.
 *
 * @param[in,out] user      Pointer to the UserCtx structure containing the DMSwarm.
 * @param[in]     fieldName Name of the reference field (e.g., "position").
 * @param[in]     ti        Time index for constructing the file name.
 * @param[in]     ext       File extension (e.g., "dat").
 * @param[out]    skipStep  Pointer to boolean flag, set to PETSC_TRUE if the step
 *                          should be skipped (e.g., file not found), PETSC_FALSE otherwise.
 *
 * @return PetscErrorCode 0 on success, non-zero on critical failure.
 *         If the reference file is not found, returns 0 and sets skipStep = PETSC_TRUE.
 */
PetscErrorCode PreCheckAndResizeSwarm(UserCtx *user, PetscInt ti, const char *ext);

/**
 * @brief Performs one full cycle of particle migration: identify, set ranks, and migrate.
 *
 * This function encapsulates the three main steps of migrating particles between MPI ranks:
 * 1. Identify particles on the local rank that need to move based on their current
 *    positions and the domain decomposition (`bboxlist`).
 * 2. Determine the destination rank for each migrating particle.
 * 3. Perform the actual migration using PETSc's `DMSwarmMigrate`.
 * It also calculates and logs the global number of particles migrated.
 *
 * @param user Pointer to the UserCtx structure.
 * @param bboxlist Array of BoundingBox structures defining the spatial domain of each MPI rank.
 * @param migrationList_p Pointer to a pointer for the MigrationInfo array. This array will be
 *                        allocated/reallocated by `IdentifyMigratingParticles` if necessary.
 *                        The caller is responsible for freeing this list eventually.
 * @param migrationCount_p Pointer to store the number of particles identified for migration
 *                         on the local rank. This is reset to 0 after migration for the current cycle.
 * @param migrationListCapacity_p Pointer to store the current capacity of the `migrationList_p` array.
 * @param currentTime Current simulation time (used for logging).
 * @param step Current simulation step number (used for logging).
 * @param migrationCycleName A descriptive name for this migration cycle (e.g., "Preliminary Sort", "Main Loop")
 *                           for logging purposes.
 * @param[out] globalMigrationCount_out Pointer to store the total number of particles migrated
 *                                      across all MPI ranks during this cycle.
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode PerformSingleParticleMigrationCycle(UserCtx *user, const BoundingBox *bboxlist,
                                                   MigrationInfo **migrationList_p, PetscInt *migrationCount_p,
                                                   PetscInt *migrationListCapacity_p,
                                                   PetscReal currentTime, PetscInt step, const char *migrationCycleName,
                                                   PetscInt *globalMigrationCount_out);


/**
 * @brief Re-initializes the positions of particles currently on this rank if this rank owns
 *        part of the designated inlet surface.
 *
 * This function is intended for `user->ParticleInitialization == 0` (Surface Initialization mode)
 * and is typically called after an initial migration step (e.g., in `PerformInitialSetup`).
 * It ensures that all particles that should originate from the inlet surface and are now
 * on the correct MPI rank are properly distributed across that rank's portion of the inlet.
 *
 * @param user Pointer to the UserCtx structure, containing simulation settings and grid information.
 * @param currentTime Current simulation time (used for logging).
 * @param step Current simulation step number (used for logging).
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode ReinitializeParticlesOnInletSurface(UserCtx *user, PetscReal currentTime, PetscInt step);

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
 * The identification is done by comparing the PIDs of particles currently on the rank
 * against a "snapshot" of PIDs taken *before* the migration occurred.
 *
 * @param[in] swarm            The DMSwarm object, which has just completed a migration.
 * @param[in] n_local_before   The number of particles that were on this rank *before* the
 *                             migration was performed.
 * @param[in] pids_before      A pre-sorted array of the PIDs that were on this rank before
 *                             the migration. This is used for fast lookups.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 *
 * @note This function assumes the `pids_before` array is sorted in ascending order to
 *       enable the use of an efficient binary search.
 */
PetscErrorCode FlagNewcomersForLocation(DM swarm,
                                        PetscInt n_local_before,
                                        const PetscInt64 pids_before[]);


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
 *     that determines the particle's status: located locally, needs migration, or is lost.
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
 */
PetscErrorCode LocateAllParticlesInGrid(UserCtx *user,BoundingBox *bboxlist);

/**
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
