// ParticleMotion.c

#include "ParticleMotion.h"

// Define a buffer size for error messages if not already available
#ifndef ERROR_MSG_BUFFER_SIZE
#define ERROR_MSG_BUFFER_SIZE 256 // Or use PETSC_MAX_PATH_LEN if appropriate
#endif

/**
 * @brief Updates a particle's position based on its velocity and the timestep dt (stored in user->dt).
 *
 * @param[in]     user     Pointer to your UserCtx (must contain user->dt).
 * @param[in,out] position Pointer to the particle's current position (Cmpnts).
 * @param[in]     velocity Pointer to the particle's velocity (Cmpnts).
 *
 * @return PetscErrorCode  Returns 0 on success, or an error code on failure.
 */
PetscErrorCode UpdateParticlePosition(UserCtx *user, Cmpnts *position, const Cmpnts *velocity)
{
  PetscFunctionBeginUser; // PETSc macro for error/stack tracing

  /* Update the position with velocity * dt */
  position->x += velocity->x * user->dt;
  position->y += velocity->y * user->dt;
  position->z += velocity->z * user->dt;

  PetscFunctionReturn(0);
}

/**
 * @brief Loops over all local particles in the DMSwarm, updating their positions
 *        based on velocity and the global timestep user->dt. 
 * @param[in,out] user    Pointer to UserCtx (must contain dt).
 *
 * @return PetscErrorCode Returns 0 on success, or an error code on failure.
 */
PetscErrorCode UpdateAllParticlePositions(UserCtx *user)
{
  PetscErrorCode ierr;
  DM swarm = user->swarm;
  PetscInt       nLocal, p;
  PetscReal        *pos = NULL;
  PetscReal        *vel = NULL;
  PetscMPIInt rank;
  Cmpnts temp_pos, temp_vel; 

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscFunctionBeginUser;  // PETSc macro for error/stack tracing

  // 1) Get the number of local particles
  ierr = DMSwarmGetLocalSize(swarm, &nLocal); CHKERRQ(ierr);
  if (nLocal == 0) {
    LOG_ALLOW(LOCAL,LOG_DEBUG,"[Rank %d] No particles to move/transport. \n",rank);
    PetscFunctionReturn(0);   // nothing to do, no fields held 
  }
  // 2) Access the "position" and "velocity" fields
  ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&pos); CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "velocity", NULL, NULL, (void**)&vel); CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL,LOG_DEBUG," [Rank %d] No.of Particles to update: %d.\n",rank,nLocal);

  // 3) Loop over all local particles, updating each position by velocity * dt
  for (p = 0; p < nLocal; p++) {
    // update temporary position struct
    temp_pos.x = pos[3*p];
    temp_pos.y = pos[3*p + 1];
    temp_pos.z = pos[3*p + 2];

    // update temporary velocity struct
    temp_vel.x = vel[3*p];
    temp_vel.y = vel[3*p + 1];
    temp_vel.z = vel[3*p + 2];    
    
    ierr = UpdateParticlePosition(user, &temp_pos, &temp_vel); CHKERRQ(ierr);
    
    // update swarm from temporary position struct
    pos[3*p] = temp_pos.x;
    pos[3*p + 1] = temp_pos.y;
    pos[3*p + 2] = temp_pos.z;

    vel[3*p] = temp_vel.x;
    vel[3*p + 1] = temp_vel.y;
    vel[3*p + 2] = temp_vel.z;
  }

  // 4) Restore the fields
  ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&pos); CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "velocity", NULL, NULL, (void**)&vel); CHKERRQ(ierr);

  LOG_ALLOW(LOCAL,LOG_DEBUG,"Particle moved/transported successfully on Rank %d.\n",rank);
  
  PetscFunctionReturn(0);
}


/**
 * @brief Checks if a particle position is within the bounds of a given bounding box.
 *
 * @param bbox Pointer to the BoundingBox structure.
 * @param pos  Pointer to the particle's position (Cmpnts).
 *
 * @return PetscBool PETSC_TRUE if the particle is inside or on the boundary, PETSC_FALSE otherwise.
 */
static inline PetscBool IsParticleInBox(const BoundingBox *bbox, const Cmpnts *pos) {
    return (pos->x >= bbox->min_coords.x && pos->x <= bbox->max_coords.x &&
            pos->y >= bbox->min_coords.y && pos->y <= bbox->max_coords.y &&
            pos->z >= bbox->min_coords.z && pos->z <= bbox->max_coords.z);
}

/**
 * @brief Checks for particles outside the global physical domain and removes them.
 *
 * This function iterates through all particles local to the current MPI rank. It checks
 * if a particle's position is inside ANY of the bounding boxes that define the subdomains
 * of the MPI ranks. The list of all bounding boxes must be available on all ranks.
 *
 * If a particle is not found within any of the provided bounding boxes, it is considered
 * to be "out of bounds" and is permanently removed from the simulation using
 * DMSwarmRemovePointAtIndex().
 *
 * This function uses a robust "Restore-Remove-Reacquire" pattern for modifying the swarm.
 * When a particle is removed, all pointers to swarm data are restored, the removal operation
 * is performed, and then the pointers are reacquired if the loop needs to continue. This is
 * the safest way to modify the swarm's structure during iteration.
 *
 * @note This function is often redundant if a robust particle location scheme (e.g., one that
 *       sets a `LOST` status) and a corresponding cleanup function (e.g.,
 *       `CheckAndRemoveLostParticles`) are already in use. Using both can introduce
 *       complexity.
 *
 * @warning This function contains a collective MPI operation (MPI_Allreduce). To avoid
 *          deadlocks, it MUST be called by ALL ranks in the communicator, even if they
 *          have no work to do. Do not place this call inside rank-specific conditional
 *          logic in your main loop.
 *
 * @param[in,out]  user             Pointer to the UserCtx structure containing the swarm and MPI info.
 * @param[out]     removedCountLocal Pointer to store the number of particles removed on THIS rank.
 * @param[out]     removedCountGlobal Pointer to store the total number of particles removed ACROSS ALL ranks.
 * @param[in]      bboxlist         An array of BoundingBox structures for ALL MPI ranks, indexed 0 to (size-1).
 *                                  This array must be up-to-date and available on all ranks.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 */
/*
PetscErrorCode CheckAndRemoveOutOfBoundsParticles(UserCtx *user,
                                              PetscInt *removedCountLocal,
                                              PetscInt *removedCountGlobal,
                                              const BoundingBox *bboxlist)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    PetscInt       nLocalInitial;
    PetscReal      *pos_p = NULL;
    PetscInt64     *pid_p = NULL; // For better logging
    PetscInt       local_removed_count = 0;
    PetscMPIInt    global_removed_count_mpi = 0;
    PetscMPIInt    rank, size;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "[Rank %d] Checking for out-of-bounds particles...\n", rank);

    // Initialize output parameters
    *removedCountLocal = 0;
    if (removedCountGlobal) *removedCountGlobal = 0;

    ierr = DMSwarmGetLocalSize(swarm, &nLocalInitial); CHKERRQ(ierr);

    if(nLocalInitial > 0){
    // Get access to swarm fields. Getting PID makes for better log messages.
    ierr = DMSwarmGetField(swarm, "position",    NULL, NULL, (void **)&pos_p); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid_p); CHKERRQ(ierr);

    // --- Iterate BACKWARDS to handle index changes safely during removal ---
    for (PetscInt p = nLocalInitial - 1; p >= 0; p--) {
        PetscBool isInsideAnyBox = PETSC_FALSE;

        // Create a temporary structure for the particle's position to pass to the check function.
        // This assumes a helper function like IsParticleInBox exists. If not, inline the logic.
        Cmpnts current_pos = {pos_p[3*p + 0], pos_p[3*p + 1], pos_p[3*p + 2]};

        // Check if the particle is inside ANY of the rank bounding boxes
        for (PetscMPIInt proc = 0; proc < size; proc++) {
            if (IsParticleInBox(&bboxlist[proc], &current_pos)) {
                isInsideAnyBox = PETSC_TRUE;
                break; // Efficiently exit inner loop once a containing box is found
            }
        }

        if (!isInsideAnyBox) {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Removing out-of-bounds particle [PID %lld] at local index %d. Pos: (%g, %g, %g)\n",
                      rank, (long long)pid_p[p], p, current_pos.x, current_pos.y, current_pos.z);

            // --- Safe Removal Pattern: Restore -> Remove -> Reacquire ---

            // 1. Restore all fields BEFORE modifying the swarm structure.
            ierr = DMSwarmRestoreField(swarm, "position",    NULL, NULL, (void **)&pos_p); CHKERRQ(ierr);
            ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid_p); CHKERRQ(ierr);

            // 2. Remove the particle at the current local index 'p'.
            ierr = DMSwarmRemovePointAtIndex(swarm, p); CHKERRQ(ierr);
            local_removed_count++;

            // 3. After removal, the swarm is modified. We MUST re-acquire pointers if we are to continue.
            PetscInt nLocalCurrent;
            ierr = DMSwarmGetLocalSize(swarm, &nLocalCurrent); CHKERRQ(ierr);
            if (nLocalCurrent > 0) { // If there are still particles left to check
                ierr = DMSwarmGetField(swarm, "position",    NULL, NULL, (void **)&pos_p); CHKERRQ(ierr);
                ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid_p); CHKERRQ(ierr);
            } else {
                // All particles on this rank have been removed. Invalidate pointers and exit the loop.
                pos_p = NULL;
                pid_p = NULL;
		//  break;
            }
        }
      } // End of backwards loop

    // Restore fields if the loop completed without a final removal (i.e., pointers are still valid).
    if (pos_p) {
        ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void **)&pos_p); CHKERRQ(ierr);
     }
    if (pid_p) {
        ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid_p); CHKERRQ(ierr);
     }
    } // End of if (nLocalInitial > 0) 

    // Get the *final* local size for logging
    PetscInt nLocalFinal;
    ierr = DMSwarmGetLocalSize(swarm, &nLocalFinal); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "[Rank %d] Finished removing %d out-of-bounds particles. Final local size: %d.\n", rank, local_removed_count, nLocalFinal);

    // --- Synchronize counts across all ranks ---
    *removedCountLocal = local_removed_count;
    if (removedCountGlobal) {
        ierr = MPI_Allreduce(&local_removed_count, &global_removed_count_mpi, 1, MPI_INT, MPI_SUM, PetscObjectComm((PetscObject)swarm)); CHKERRQ(ierr);
        *removedCountGlobal = global_removed_count_mpi;
        LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "[Rank %d] Removed %d out-of-bounds particles globally.\n", rank, *removedCountGlobal);
    }

    PetscFunctionReturn(0);
}
*/

/**
 * @brief Checks for particles outside the global physical domain and removes them.
 *
 * This function iterates through all particles local to the current MPI rank. It determines
 * if a particle's physical position is outside of ALL subdomains owned by the MPI processes.
 * To perform this check, it requires a list of the bounding boxes for every MPI rank's
 * subdomain.
 *
 * If a particle is not found within any of the subdomains, it is considered "out of bounds"
 * and is permanently removed from the simulation using DMSwarmRemovePointAtIndex().
 *
 * The function is carefully designed to handle modifications to the DMSwarm during iteration
 * safely. It does this by:
 *   1. Iterating BACKWARDS through the local particle list, so removing an element at
 *      index 'p' does not affect the indices of subsequent elements to be checked (p-1, p-2, ...).
 *   2. Using a robust "Restore-Remove-Reacquire" pattern. When a particle is removed,
 *      all pointers to swarm data are first restored, the removal operation is performed,
 *      and then the pointers are re-acquired for the remainder of the loop.
 *
 * @warning This function contains a collective MPI operation (MPI_Allreduce) at the end.
 *          To avoid deadlocks, it MUST be called by ALL ranks in the communicator, even if
 *          a rank has no local particles. Do not place this call inside rank-specific
 *          conditional logic in your main loop.
 *
 * @note This function can be redundant if a robust particle location scheme (e.g., one that
 *       sets a `LOST` status) and a corresponding cleanup function are already in use.
 *
 * @param[in,out]  user             Pointer to the UserCtx structure containing the swarm and MPI info.
 * @param[out]     removedCountLocal Pointer to an integer that will store the number of particles removed on THIS rank.
 * @param[out]     removedCountGlobal Pointer to an integer that will store the total number of particles removed ACROSS ALL ranks.
 * @param[in]      bboxlist         An array of BoundingBox structures for ALL MPI ranks, indexed 0 to (size-1).
 *                                  This array must be up-to-date and available on every rank.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 */
PetscErrorCode CheckAndRemoveOutOfBoundsParticles(UserCtx *user,
                                              PetscInt *removedCountLocal,
                                              PetscInt *removedCountGlobal,
                                              const BoundingBox *bboxlist)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    PetscInt       nLocalInitial;
    PetscReal      *pos_p = NULL;
    PetscInt64     *pid_p = NULL; // For better logging
    PetscInt       local_removed_count = 0;
    PetscMPIInt    global_removed_count_mpi = 0;
    PetscMPIInt    rank, size;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "[Rank %d] Checking for out-of-bounds particles...", rank);

    // Initialize output parameters to ensure clean state
    *removedCountLocal = 0;
    if (removedCountGlobal) *removedCountGlobal = 0;

    ierr = DMSwarmGetLocalSize(swarm, &nLocalInitial); CHKERRQ(ierr);

    // Only proceed if there are particles to check on this rank.
    // All ranks will still participate in the final collective MPI_Allreduce.
    if (nLocalInitial > 0) {
        // Get access to swarm fields once before the loop begins.
        ierr = DMSwarmGetField(swarm, "position",    NULL, NULL, (void **)&pos_p); CHKERRQ(ierr);
        ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid_p); CHKERRQ(ierr);

        // --- Iterate BACKWARDS to handle index changes safely during removal ---
        for (PetscInt p = nLocalInitial - 1; p >= 0; p--) {
            PetscBool isInsideAnyBox = PETSC_FALSE;
            Cmpnts current_pos = {pos_p[3*p + 0], pos_p[3*p + 1], pos_p[3*p + 2]};

            // Check if the particle is inside ANY of the rank bounding boxes
            for (PetscMPIInt proc = 0; proc < size; proc++) {
                if (IsParticleInBox(&bboxlist[proc], &current_pos)) {
                    isInsideAnyBox = PETSC_TRUE;
                    break; // Particle is inside a valid domain, stop checking.
                }
            }

            if (!isInsideAnyBox) {
                LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Removing out-of-bounds particle [PID %lld] at local index %d. Pos: (%g, %g, %g)\n",
                          rank, (long long)pid_p[p], p, current_pos.x, current_pos.y, current_pos.z);

                // --- Safe Removal Pattern: Restore -> Remove -> Reacquire ---
                // This is the fix for the double-restore bug. Pointers are managed carefully
                // within this block and then restored cleanly after the loop.

                // 1. Restore all fields BEFORE modifying the swarm structure. This invalidates pos_p and pid_p.
                ierr = DMSwarmRestoreField(swarm, "position",    NULL, NULL, (void **)&pos_p); CHKERRQ(ierr);
                ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid_p); CHKERRQ(ierr);

                // 2. Remove the particle at the current local index 'p'.
                ierr = DMSwarmRemovePointAtIndex(swarm, p); CHKERRQ(ierr);
                local_removed_count++;

                // 3. After removal, re-acquire pointers ONLY if the loop is not finished.
                PetscInt nLocalCurrent;
                ierr = DMSwarmGetLocalSize(swarm, &nLocalCurrent); CHKERRQ(ierr);

                if (nLocalCurrent > 0 && p > 0) { // Check if there are particles left AND iterations left
                    ierr = DMSwarmGetField(swarm, "position",    NULL, NULL, (void **)&pos_p); CHKERRQ(ierr);
                    ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid_p); CHKERRQ(ierr);
                } else {
                    // All remaining particles were removed OR this was the last particle (p=0).
                    // Invalidate pointers to prevent the final restore call and exit the loop.
                    pos_p = NULL;
                    pid_p = NULL;
                    break;
                }
            }
        } // End of backwards loop

        // At the end, restore any valid pointers. This handles three cases:
        // 1. No particles were removed: restores the original pointers.
        // 2. Particles were removed mid-loop: restores the pointers from the last re-acquisition.
        // 3. All particles were removed: pointers are NULL, so nothing is done.
        if (pos_p) { ierr = DMSwarmRestoreField(swarm, "position",    NULL, NULL, (void **)&pos_p); CHKERRQ(ierr); }
        if (pid_p) { ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid_p); CHKERRQ(ierr); }
    } // End of if (nLocalInitial > 0)

    PetscInt nLocalFinal;
    ierr = DMSwarmGetLocalSize(swarm, &nLocalFinal); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "[Rank %d] Finished removing %d out-of-bounds particles. Final local size: %d.\n", rank, local_removed_count, nLocalFinal);

    // --- Synchronize counts across all ranks ---
    *removedCountLocal = local_removed_count;
    if (removedCountGlobal) {
        ierr = MPI_Allreduce(&local_removed_count, &global_removed_count_mpi, 1, MPI_INT, MPI_SUM, PetscObjectComm((PetscObject)swarm)); CHKERRQ(ierr);
        *removedCountGlobal = global_removed_count_mpi;
        // Use a synchronized log message so only one rank prints the global total.
        LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "[Rank %d] Removed %d out-of-bounds particles globally.\n", rank, *removedCountGlobal);
    }

    PetscFunctionReturn(0);
}

/**
 * @brief Removes particles that have been definitively flagged as LOST by the location algorithm.
 *
 * This function is the designated cleanup utility for particles that have exited the
 * physical domain or could not be located for any other reason. It should be called
 * after a particle location and migration phase is complete.
 *
 * It iterates through all locally owned particles and checks their `DMSwarm_location_status`
 * field. If a particle's status is `LOST`, it is permanently removed from the simulation
 * using `DMSwarmRemovePointAtIndex`.
 *
 * The function uses a robust "Restore-Remove-Reacquire" pattern for modifying the swarm.
 * When a particle is removed, all pointers to swarm data are restored, the removal operation
 * is performed, and then the pointers are reacquired if the loop needs to continue.
 *
 * @warning This function contains a collective MPI operation (MPI_Allreduce). To avoid
 *          deadlocks, it MUST be called by ALL ranks in the communicator. This implementation
 *          is designed to be safe by ensuring all ranks, even those with zero local
 *          particles, participate in the collective call.
 *
 * @param[in,out]  user              Pointer to the UserCtx structure containing the swarm.
 * @param[out]     removedCountLocal Pointer to store the number of particles removed on THIS rank.
 * @param[out]     removedCountGlobal Pointer to store the total number of particles removed ACROSS ALL ranks.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 */
/*
PetscErrorCode CheckAndRemoveLostParticles(UserCtx *user,
                                           PetscInt *removedCountLocal,
                                           PetscInt *removedCountGlobal)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    PetscInt       nLocalInitial;
    PetscInt       *status_p = NULL;
    PetscInt64     *pid_p = NULL; // For better logging
    PetscReal      *pos_p = NULL; // For better logging
    PetscInt       local_removed_count = 0;
    PetscMPIInt    global_removed_count_mpi = 0;
    PetscMPIInt    rank;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Checking for and removing LOST particles...\n", rank);

    // Initialize output parameters
    *removedCountLocal = 0;
    if (removedCountGlobal) *removedCountGlobal = 0;

    ierr = DMSwarmGetLocalSize(swarm, &nLocalInitial); CHKERRQ(ierr);

    // Guard the local processing loop. If a rank has no particles, it will skip
    // this block entirely but will still proceed to the collective MPI_Allreduce call.
    // This is the critical fix to prevent deadlocks.
    if (nLocalInitial > 0) {
        // Get access to swarm fields. We get all three for robust processing and logging.
        ierr = DMSwarmGetField(swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status_p); CHKERRQ(ierr);
        ierr = DMSwarmGetField(swarm, "DMSwarm_pid",             NULL, NULL, (void **)&pid_p);    CHKERRQ(ierr);
        ierr = DMSwarmGetField(swarm, "position",                NULL, NULL, (void **)&pos_p);    CHKERRQ(ierr);

        // --- Iterate BACKWARDS to handle index changes safely during removal ---
        for (PetscInt p = nLocalInitial - 1; p >= 0; p--) {
            if (status_p[p] == LOST) {
	      LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Removing LOST particle [PID %ld] at local index %d. Position: (%.4f, %.4f, %.4f).\n",rank, pid_p[p],p,pos_p[3*p], pos_p[3*p+1], pos_p[3*p+2]);

                // --- Safe Removal Pattern: Restore -> Remove -> Reacquire ---

                // 1. Restore all fields BEFORE modifying the swarm structure.
                ierr = DMSwarmRestoreField(swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status_p); CHKERRQ(ierr);
                ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid",             NULL, NULL, (void **)&pid_p);    CHKERRQ(ierr);
                ierr = DMSwarmRestoreField(swarm, "position",                NULL, NULL, (void **)&pos_p);    CHKERRQ(ierr);

                // 2. Remove the particle at the current local index 'p'.
                ierr = DMSwarmRemovePointAtIndex(swarm, p); CHKERRQ(ierr);
                local_removed_count++;

                // 3. After removal, the swarm is modified. We MUST re-acquire pointers if we are to continue.
                PetscInt nLocalCurrent;
                ierr = DMSwarmGetLocalSize(swarm, &nLocalCurrent); CHKERRQ(ierr);
                if (nLocalCurrent > 0) { // If there are still particles left to check
                    ierr = DMSwarmGetField(swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status_p); CHKERRQ(ierr);
                    ierr = DMSwarmGetField(swarm, "DMSwarm_pid",             NULL, NULL, (void **)&pid_p);    CHKERRQ(ierr);
                    ierr = DMSwarmGetField(swarm, "position",                NULL, NULL, (void **)&pos_p);    CHKERRQ(ierr);
                } else {
                    // All particles on this rank have been removed. Invalidate pointers and exit the loop.
                    status_p = NULL;
                    pid_p = NULL;
                    pos_p = NULL;
                    break;
                }
            }
        } // End of backwards loop

        // Restore fields if the loop completed without a final removal (i.e., pointers are still valid).
        if (status_p) { ierr = DMSwarmRestoreField(swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status_p); CHKERRQ(ierr); }
        if (pid_p)    { ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid",             NULL, NULL, (void **)&pid_p);    CHKERRQ(ierr); }
        if (pos_p)    { ierr = DMSwarmRestoreField(swarm, "position",                NULL, NULL, (void **)&pos_p);    CHKERRQ(ierr); }
    } // End of if (nLocalInitial > 0)

    // Get the *final* local size for logging
    PetscInt nLocalFinal;
    ierr = DMSwarmGetLocalSize(swarm, &nLocalFinal); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Finished removing %d LOST particles. Final local size: %d.\n", rank, local_removed_count, nLocalFinal);

    // --- Synchronize counts across all ranks ---
    // ALL ranks MUST reach this point to participate in the collective call.
    *removedCountLocal = local_removed_count;
    if (removedCountGlobal) {
        ierr = MPI_Allreduce(&local_removed_count, &global_removed_count_mpi, 1, MPI_INT, MPI_SUM, PetscObjectComm((PetscObject)swarm)); CHKERRQ(ierr);
        *removedCountGlobal = global_removed_count_mpi;
        LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "[Rank %d] Removed %d LOST particles globally.\n", rank, *removedCountGlobal);
    }

    PetscFunctionReturn(0);
}
*/

/**
 * @brief Removes particles that have been definitively flagged as LOST by the location algorithm.
 *
 * This function is the designated cleanup utility for particles that have exited the
 * physical domain or could not be located for any other reason. It should be called
 * after a particle location and migration phase is complete.
 *
 * It iterates through all locally owned particles and checks their `DMSwarm_location_status`
 * field. If a particle's status is `LOST`, it is permanently removed from the simulation
 * using `DMSwarmRemovePointAtIndex`.
 *
 * The function is carefully designed to handle modifications to the DMSwarm during iteration
 * safely. It does this by:
 *   1. Iterating BACKWARDS through the local particle list, so removing an element at
 *      index 'p' does not affect the indices of subsequent elements to be checked (p-1, p-2, ...).
 *   2. Using a robust "Restore-Remove-Reacquire" pattern. When a particle is removed,
 *      all pointers to swarm data are first restored, the removal operation is performed,
 *      and then the pointers are re-acquired for the remainder of the loop.
 *
 * @warning This function contains a collective MPI operation (MPI_Allreduce) at the end.
 *          To avoid deadlocks, it MUST be called by ALL ranks in the communicator, even if
 *          a rank has no local particles. Do not place this call inside rank-specific
 *          conditional logic in your main loop.
 *
 * @param[in,out]  user              Pointer to the UserCtx structure containing the swarm.
 * @param[out]     removedCountLocal Pointer to an integer that will store the number of particles removed on THIS rank.
 * @param[out]     removedCountGlobal Pointer to an integer that will store the total number of particles removed ACROSS ALL ranks.
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 */
PetscErrorCode CheckAndRemoveLostParticles(UserCtx *user,
                                           PetscInt *removedCountLocal,
                                           PetscInt *removedCountGlobal)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    PetscInt       nLocalInitial;
    PetscInt       *status_p = NULL;
    PetscInt64     *pid_p = NULL; // For better logging
    PetscReal      *pos_p = NULL; // For better logging
    PetscInt       local_removed_count = 0;
    PetscMPIInt    global_removed_count_mpi = 0;
    PetscMPIInt    rank;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Checking for and removing LOST particles...", rank);

    // Initialize output parameters to ensure clean state
    *removedCountLocal = 0;
    if (removedCountGlobal) *removedCountGlobal = 0;

    ierr = DMSwarmGetLocalSize(swarm, &nLocalInitial); CHKERRQ(ierr);

    // Only proceed if there are particles to check on this rank.
    // All ranks will still participate in the final collective MPI_Allreduce.
    if (nLocalInitial > 0) {
        // Get access to all swarm fields once before the loop begins.
        ierr = DMSwarmGetField(swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status_p); CHKERRQ(ierr);
        ierr = DMSwarmGetField(swarm, "DMSwarm_pid",             NULL, NULL, (void **)&pid_p);    CHKERRQ(ierr);
        ierr = DMSwarmGetField(swarm, "position",                NULL, NULL, (void **)&pos_p);    CHKERRQ(ierr);

        // --- Iterate BACKWARDS to handle index changes safely during removal ---
        for (PetscInt p = nLocalInitial - 1; p >= 0; p--) {
            if (status_p[p] == LOST) {
                LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Removing LOST particle [PID %lld] at local index %d. Position: (%.4f, %.4f, %.4f).\n",
                          rank, (long long)pid_p[p], p, pos_p[3*p], pos_p[3*p+1], pos_p[3*p+2]);

                // --- Safe Removal Pattern: Restore -> Remove -> Reacquire ---
                // This is the fix for the double-restore bug. Pointers are managed carefully
                // within this block and then restored cleanly after the loop.

                // 1. Restore all fields BEFORE modifying the swarm structure. This invalidates all pointers.
                ierr = DMSwarmRestoreField(swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status_p); CHKERRQ(ierr);
                ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid",             NULL, NULL, (void **)&pid_p);    CHKERRQ(ierr);
                ierr = DMSwarmRestoreField(swarm, "position",                NULL, NULL, (void **)&pos_p);    CHKERRQ(ierr);

                // 2. Remove the particle at the current local index 'p'.
                ierr = DMSwarmRemovePointAtIndex(swarm, p); CHKERRQ(ierr);
                local_removed_count++;

                // 3. After removal, re-acquire pointers ONLY if the loop is not finished.
                PetscInt nLocalCurrent;
                ierr = DMSwarmGetLocalSize(swarm, &nLocalCurrent); CHKERRQ(ierr);

                if (nLocalCurrent > 0 && p > 0) { // Check if there are particles left AND iterations left
                    ierr = DMSwarmGetField(swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status_p); CHKERRQ(ierr);
                    ierr = DMSwarmGetField(swarm, "DMSwarm_pid",             NULL, NULL, (void **)&pid_p);    CHKERRQ(ierr);
                    ierr = DMSwarmGetField(swarm, "position",                NULL, NULL, (void **)&pos_p);    CHKERRQ(ierr);
                } else {
                    // All remaining particles were removed OR this was the last particle (p=0).
                    // Invalidate pointers to prevent the final restore call and exit the loop.
                    status_p = NULL;
                    pid_p = NULL;
                    pos_p = NULL;
                    break;
                }
            }
        } // End of backwards loop

        // At the end, restore any valid pointers. This handles three cases:
        // 1. No particles were removed: restores the original pointers.
        // 2. Particles were removed mid-loop: restores the pointers from the last re-acquisition.
        // 3. All particles were removed: pointers are NULL, so nothing is done.
        if (status_p) { ierr = DMSwarmRestoreField(swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status_p); CHKERRQ(ierr); }
        if (pid_p)    { ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid",             NULL, NULL, (void **)&pid_p);    CHKERRQ(ierr); }
        if (pos_p)    { ierr = DMSwarmRestoreField(swarm, "position",                NULL, NULL, (void **)&pos_p);    CHKERRQ(ierr); }
    } // End of if (nLocalInitial > 0)

    PetscInt nLocalFinal;
    ierr = DMSwarmGetLocalSize(swarm, &nLocalFinal); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Finished removing %d LOST particles. Final local size: %d.\n", rank, local_removed_count, nLocalFinal);

    // --- Synchronize counts across all ranks ---
    *removedCountLocal = local_removed_count;
    if (removedCountGlobal) {
        ierr = MPI_Allreduce(&local_removed_count, &global_removed_count_mpi, 1, MPI_INT, MPI_SUM, PetscObjectComm((PetscObject)swarm)); CHKERRQ(ierr);
        *removedCountGlobal = global_removed_count_mpi;
        // Use a synchronized log message so only one rank prints the global total.
        LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "[Rank %d] Removed %d LOST particles globally.\n", rank, *removedCountGlobal);
    }

    PetscFunctionReturn(0);
}


/**
 * @brief Identifies particles that have left the local MPI rank's domain and determines their target rank.
 *
 * This function iterates through all particles currently local to this MPI rank.
 * For each particle, it first checks if the particle is still within the rank's
 * primary bounding box (defined by `user->bbox`).
 *
 * If a particle is outside `user->bbox`:
 * 1. It preferentially checks the 6 immediate Cartesian neighbors (defined in `user->neighbors`).
 *    The bounding boxes for these neighbors are looked up in the global `bboxlist`.
 * 2. If the particle is not found in any of the immediate Cartesian neighbors, a fallback
 *    search is initiated. This fallback search iterates through *all other* MPI ranks
 *    (excluding the current rank and already checked immediate neighbors if optimized)
 *    using the `bboxlist` to find a rank whose bounding box contains the particle.
 *
 * Particles successfully assigned a `targetRank` are added to the `migrationList`.
 * If a particle leaves the local box but is not found in any other rank's bounding box
 * (even after the fallback), a warning is logged, as this particle might be lost from
 * the global computational domain (assuming `bboxlist` covers the entire domain).
 * Such "lost" particles should ideally be handled by a separate global boundary condition
 * check (e.g., `CheckAndRemoveOutOfBoundsParticles`) prior to calling this function.
 *
 * The `migrationList` is dynamically reallocated if its current capacity is exceeded.
 *
 * @param[in]      user           Pointer to the UserCtx structure. It must contain:
 *                                - `swarm`: The DMSwarm object.
 *                                - `bbox`: The BoundingBox of the current MPI rank.
 *                                - `neighbors`: A RankNeighbors struct with the ranks of the 6 Cartesian neighbors.
 *                                - `size`: The total number of MPI ranks in PETSC_COMM_WORLD (implicitly via MPI_Comm_size).
 * @param[in]      bboxlist       An array of BoundingBox structures for ALL MPI ranks, indexed 0 to (size-1).
 *                                This array must be up-to-date and available on all ranks.
 * @param[in,out]  migrationList  Pointer to an array of MigrationInfo structures. This array will be
 *                                populated with particles to be migrated. The function may reallocate
 *                                this array if more space is needed. The caller is responsible for
 *                                eventually freeing this memory if it's not NULL.
 * @param[out]     migrationCount Pointer to a PetscInt that will be set to the number of particles
 *                                identified for migration (i.e., the number of valid entries in `*migrationList`).
 * @param[in,out]  listCapacity   Pointer to a PetscInt representing the current allocated capacity of
 *                                `*migrationList` (in terms of number of MigrationInfo structs).
 *                                This function will update it if `*migrationList` is reallocated.
 *
 * @return PetscErrorCode 0 on success, non-zero on PETSc/MPI errors or if essential fields are missing.
 *
 * @note It is assumed that `user->neighbors` contains valid rank identifiers or `MPI_PROC_NULL`.
 * @note It is assumed that `bboxlist` is correctly populated and broadcast to all ranks before calling this.
 * @note For the fallback search to be effective, `bboxlist` should represent a tiling of the
 *       global domain. Particles not found in any box in `bboxlist` are effectively outside this
 *       tiled global domain.
 */
PetscErrorCode IdentifyMigratingParticles(UserCtx *user,
					  const BoundingBox *bboxlist,
					  MigrationInfo **migrationList,
					  PetscInt *migrationCount,
					  PetscInt *listCapacity)
{
  PetscErrorCode ierr;
  DM             swarm = user->swarm;
  PetscInt       nLocal, p;
  Cmpnts        *pos = NULL;
  PetscMPIInt    rank,size;
  BoundingBox    localBBox = user->bbox;
  RankNeighbors  neighbors = user->neighbors; // Use stored neighbors
  
  //PetscInt       currentMigrationCount = 0;
  //PetscInt       currentListCapacity = *listCapacity;
  //MigrationInfo *currentMigrationList = *migrationList;
  // Add PID pointer if logging PIDs
  PetscInt64    *pids = NULL;
  PetscFunctionBeginUser;

  // --- Input Validation and Initialization ---
  if (!user) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx 'user' is NULL.");
  if (!bboxlist) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Global 'bboxlist' is NULL.");
  if (!migrationList) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "'migrationList' output pointer is NULL.");
  if (!migrationCount) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "'migrationCount' output pointer is NULL.");
  if (!listCapacity) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "'listCapacity' output pointer is NULL.");

  if (!swarm) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "UserCtx->swarm is NULL.");
  
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
  
  ierr = DMSwarmGetLocalSize(swarm, &nLocal); CHKERRQ(ierr);
  if (nLocal == 0) {
    *migrationCount = 0;
    // Ensure output pointers are consistent even if no allocation happened
    // *migrationList = currentMigrationList;
    // *listCapacity = currentListCapacity;
    PetscFunctionReturn(0);
  }
  // Get read-only access to position
  ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void **)&pos); CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void **)&pids); CHKERRQ(ierr); // If logging PIDs
  if (!pos || !pids) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Could not access required DMSwarm fields.");
  /*
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
    "[IdentifyMigratingParticles - Rank %d INCOMING user->neighbors] xm=%d, xp=%d, ym=%d, yp=%d, zm=%d, zp=%d. MPI_PROC_NULL is %d. \n",
    rank, user->neighbors.rank_xm, user->neighbors.rank_xp,
    user->neighbors.rank_ym, user->neighbors.rank_yp,
    user->neighbors.rank_zm, user->neighbors.rank_zp, (int)MPI_PROC_NULL);
  */
  
  *migrationCount = 0;
  //currentMigrationCount = 0; // Reset count for this call
  for (p = 0; p < nLocal; p++) {
    // Check if particle is OUTSIDE the local bounding box
    if (!IsParticleInBox(&localBBox, &pos[p]))
      {
	PetscInt targetRank = MPI_PROC_NULL; // Target rank not yet found

	// Determine likely exit direction(s) to prioritize neighbor check
	PetscBool exit_xm = pos[p].x < localBBox.min_coords.x;
	PetscBool exit_xp = pos[p].x > localBBox.max_coords.x;
	PetscBool exit_ym = pos[p].y < localBBox.min_coords.y;
	PetscBool exit_yp = pos[p].y > localBBox.max_coords.y;
	PetscBool exit_zm = pos[p].z < localBBox.min_coords.z;
	PetscBool exit_zp = pos[p].z > localBBox.max_coords.z;

	// DEBUG ------------
	/*
	  if (rank == 1 && p < 5) { // Log for first few particles on Rank 1
	  LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
	  "[Identify - Rank 1, p=%d] Particle Pos(%.2f,%.2f,%.2f). Exits(xm%d,xp%d,ym%d,yp%d,zm%d,zp%d). Neighbors(xm%d,xp%d,ym%d,yp%d,zm%d,zp%d)",
	  p, pos[p].x, pos[p].y, pos[p].z,
	  exit_xm, exit_xp, exit_ym, exit_yp, exit_zm, exit_zp,
	  neighbors.rank_xm, neighbors.rank_xp, neighbors.rank_ym, neighbors.rank_yp, neighbors.rank_zm, neighbors.rank_zp);
	  }
	*/
	// DEBUG -----------
        
	
	LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Particle %d[PID = %ld] at (%g,%g,%g) left local bbox. Checking neighbors...\n",rank, p,pids[p], pos[p].x, pos[p].y, pos[p].z);

	//1.  Check neighbors preferentially
	if (exit_xm && neighbors.rank_xm != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_xm], &pos[p])) targetRank = neighbors.rank_xm;
	else if (exit_xp && neighbors.rank_xp != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_xp], &pos[p])) targetRank = neighbors.rank_xp;
	else if (exit_ym && neighbors.rank_ym != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_ym], &pos[p])) targetRank = neighbors.rank_ym;
	else if (exit_yp && neighbors.rank_yp != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_yp], &pos[p])) targetRank = neighbors.rank_yp;
	else if (exit_zm && neighbors.rank_zm != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_zm], &pos[p])) targetRank = neighbors.rank_zm;
	else if (exit_zp && neighbors.rank_zp != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors.rank_zp], &pos[p])) targetRank = neighbors.rank_zp;
	// Add checks for edge/corner neighbors if needed and if they were stored

	
	// 2.--- Fallback (if strict Cartesian neighbors aren't enough) ---
	
	if (targetRank == MPI_PROC_NULL){
	  /* ... loop through all ranks in bboxlist ... */
        
	  LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Particle %" PetscInt64_FMT " not in immediate neighbors. Fallback: checking all %d other ranks...",rank, pids[p],size);

	  for(PetscMPIInt r = 0; r < size; r++){

	    if(IsParticleInBox(&bboxlist[r], &pos[p])){

	      targetRank = r;

	      LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d] Particle %ld FOUND in bboxlist[%d] during fallback search. \n",rank, pids[p],r);
	      break;
	    }
	    
	  }
	  
	}
	// 3. ---- Add to migration list if a target rank was found.
	if (targetRank != MPI_PROC_NULL){

	  /* OLD LOGIC
          // Resize list if needed (using PetscRealloc for safety)
	  if (currentMigrationCount >= currentListCapacity) {
	    PetscInt OldCapacity = currentListCapacity;
	    
	      PetscInt newCapacity = (currentListCapacity == 0) ? 16 : currentListCapacity * 2;

	    ierr = PetscRealloc((size_t)newCapacity * sizeof(MigrationInfo),
				migrationList); CHKERRQ(ierr);

	    currentMigrationList = *migrationList;
	    
	    
	    currentListCapacity = newCapacity;

	    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Reallocated migrationList capacity from %d to %d\n", rank,OldCapacity, newCapacity);
	   }
	    
	  // Add to migration list
	  currentMigrationList[currentMigrationCount].local_index = p;
	  currentMigrationList[currentMigrationCount].target_rank = targetRank;
	  // currentMigrationList[currentMigrationCount].pid = pids[p]; // If storing PID
	  currentMigrationCount++;
	  */

	  // NEW LOGIC (TEST)
          
	  ierr = AddToMigrationList(migrationList,    // Pointer to the list pointer
				    listCapacity,     // Pointer to the capacity variable
				    migrationCount,   // Pointer to the count variable
				    p,                // The particle's local index
				    targetRank);      // The destination rank
	  CHKERRQ(ierr);
	  
	  LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Particle %d marked for migration to rank %d.\n", rank, p, targetRank);
	  
	}
	  else {
	  // Particle left local box but was not found in any *checked* neighbor box.
	  // Since CheckAndRemove should have run first, this might indicate a particle
	  // moved more than one cell width into a diagonal neighbor's domain not checked here,
	  // or there's an issue with BBox overlap/gaps.
          LOG_ALLOW(LOCAL, LOG_WARNING, "Rank %d: Particle %d[PID = %ld] at (%g,%g,%g) left local bbox but target neighbor rank not found! (May be lost or need wider neighbor check).\n",
                    rank, p, pids[p],pos[p].x, pos[p].y, pos[p].z);
          // Consider marking for removal here if this case should not happen:
          // Maybe add to a separate 'removalList' or use the marking technique
          // from the CheckAndRemove function if it's readily available.
	  }
      } // end if isOutsideLocal
  } // end particle loop
  ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void **)&pos); CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void **)&pids); CHKERRQ(ierr);
  //  *migrationList = currentMigrationList;
  // *migrationCount = currentMigrationCount;
  //  *listCapacity = currentListCapacity;
  
  //LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Identified %d particles for potential migration.\n", rank, currentMigrationCount);
  LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Identified %d particles for potential migration.\n", rank, *migrationCount);
  PetscFunctionReturn(0);
}


/**
 * @brief Sets the target rank field (DMSwarmPICField_rank) for particles scheduled for migration.
 *
 * @param user           Pointer to UserCtx pa(contains swarm).
 * @param migrationList  Array of MigrationInfo structs containing local indices and target ranks.
 * @param migrationCount Number of particles in the migrationList.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode SetMigrationRanks(UserCtx* user, const MigrationInfo *migrationList, PetscInt migrationCount)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    PetscInt       p_idx;
    PetscInt      *rankField = NULL; // Field storing target rank

    PetscFunctionBeginUser;

    // Ensure the migration rank field exists
    ierr = DMSwarmGetField(swarm, "DMSwarm_rank", NULL, NULL, (void **)&rankField); CHKERRQ(ierr);

    // Set the target rank for migrating particles
    for(p_idx = 0; p_idx < migrationCount; ++p_idx) {
        rankField[migrationList[p_idx].local_index] = migrationList[p_idx].target_rank;
    }

    ierr = DMSwarmRestoreField(swarm, "DMSwarm_rank", NULL, NULL, (void **)&rankField); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

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
PetscErrorCode PerformMigration(UserCtx *user)
{
    PetscErrorCode ierr;
    DM swarm = user->swarm;
    PetscMPIInt rank;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Rank %d: Starting DMSwarmMigrate...\n", rank);

    // Perform the migration - PETSC_TRUE removes particles that fail to land
    // in a valid cell on the target rank (or were marked with an invalid rank).
    ierr = DMSwarmMigrate(swarm, PETSC_TRUE); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Rank %d: Migration complete.\n", rank);
    PetscFunctionReturn(0);
}

//-----------------------------------------------------------------------------
// MODULE (COUNT): Calculates Particle Count Per Cell - REVISED FOR DMDAVecGetArray
//-----------------------------------------------------------------------------
/**
 * @brief Counts particles in each cell of the DMDA 'da' and stores the result in user->ParticleCount.
 * @ingroup scatter_module
 *
 * Zeros the user->ParticleCount vector, then iterates through local particles.
 * Reads the **GLOBAL** cell index (I, J, K) stored in the "DMSwarm_CellID" field.
 * Uses DMDAVecGetArray to access the local portion of the count vector and increments
 * the count at the global index (I, J, K) if it belongs to the local patch (including ghosts).
 *
 * @param[in,out] user Pointer to the UserCtx structure containing da, swarm, and ParticleCount.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.

 */
PetscErrorCode CalculateParticleCountPerCell(UserCtx *user) {
    PetscErrorCode ierr;
    DM             da = user->da;
    DM             swarm = user->swarm;
    Vec            countVec = user->ParticleCount;
    PetscInt       nlocal, p;
    PetscInt       *global_cell_id_arr; // Read GLOBAL cell IDs
    PetscScalar    ***count_arr_3d;     // Use 3D accessor
    PetscInt64       *PID_arr;
    PetscMPIInt    rank;
    char           msg[ERROR_MSG_BUFFER_SIZE];
    PetscInt       particles_counted_locally = 0;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // --- Input Validation ---
    if (!da) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx->da is NULL.");
    if (!swarm) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx->swarm is NULL.");
    if (!countVec) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx->ParticleCount is NULL.");
    // Check DOF of da
    PetscInt count_dof;
    ierr = DMDAGetInfo(da, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &count_dof, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
     if (count_dof != 1) { PetscSNPrintf(msg, sizeof(msg), "countDM must have DOF=1, got %d.", count_dof); SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, msg); }

    // --- Zero the count vector ---
    ierr = VecSet(countVec, 0.0); CHKERRQ(ierr);

    // --- Get Particle Data ---
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "CalculateParticleCountPerCell: Accessing particle data.\n");
    ierr = DMSwarmGetLocalSize(swarm, &nlocal); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm,"DMSwarm_CellID", NULL, NULL, (void **)&global_cell_id_arr); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm,"DMSwarm_pid",NULL,NULL,(void **)&PID_arr);CHKERRQ(ierr);

    // --- Get Grid Vector Array using DMDA accessor ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "CalculateParticleCountPerCell: Accessing ParticleCount vector array (using DMDAVecGetArray).\n");
    ierr = DMDAVecGetArray(da, countVec, &count_arr_3d); CHKERRQ(ierr);

    // Get local owned range for validation/logging if needed, but not for indexing with DMDAVecGetArray
    PetscInt xs, ys, zs, xm, ym, zm;
    ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm); CHKERRQ(ierr);

    // --- Accumulate Counts Locally ---
    LOG_ALLOW(LOCAL, LOG_DEBUG, "CalculateParticleCountPerCell (Rank %d): Processing %d local particles using GLOBAL CellIDs.\n",rank,nlocal);
    for (p = 0; p < nlocal; p++) {
        // Read the GLOBAL indices stored for this particle
        PetscInt64 i = global_cell_id_arr[p * 3 + 0]; // Global i index
        PetscInt64 j = global_cell_id_arr[p * 3 + 1]; // Global j index
        PetscInt64 k = global_cell_id_arr[p * 3 + 2]; // Global k index

        // *** Bounds check is implicitly handled by DMDAVecGetArray for owned+ghost region ***
        // However, accessing outside this region using global indices WILL cause an error.
        // A preliminary check might still be wise if global IDs could be wild.
        // We rely on LocateAllParticles to provide valid global indices [0..IM-1] etc.
    // *** ADD PRINTF HERE ***
	LOG_LOOP_ALLOW(LOCAL,LOG_DEBUG,p,100,"[Rank %d CalcCount] Read CellID for p=%d, PID = %ld: (%ld, %ld, %ld)\n", rank, p,PID_arr[p],i, j, k);
    // *** END PRINTF ***
        // *** Access the local array using GLOBAL indices ***
        // DMDAVecGetArray allows this, mapping global (I,J,K) to the correct
        // location within the local ghosted array segment.
        // This assumes (I,J,K) corresponds to a cell within the local owned+ghost region.
        // If (I,J,K) is for a cell owned by *another* rank (and not in the ghost region
        // of this rank), accessing count_arr_3d[K][J][I] will likely lead to a crash
        // or incorrect results. This highlights why storing LOCAL indices is preferred
        // for parallel runs. But proceeding with GLOBAL as requested:
        if (i >= xs && i < xs + xm  && 
            j >= ys && j < ys + ym  && // Adjust based on actual ghost width
            k >= zs && k < zs + zm  )   // This check prevents definite crashes but doesn't guarantee ownership
        {
             // Increment count at the location corresponding to GLOBAL index (I,J,K)
	  //  LOG_ALLOW(LOCAL, LOG_DEBUG, "CalculateParticleCountPerCell (Rank %d): Particle %d with global CellID (%d, %d, %d) incremented with a particle.\n",rank, p, i, j, k);
             count_arr_3d[k][j][i] += 1.0;
             particles_counted_locally++;
         } else {
              // This particle's global ID is likely outside the range this rank handles (even ghosts)
              LOG_ALLOW(LOCAL, LOG_DEBUG, "CalculateParticleCountPerCell (Rank %d): Skipping particle %ld with global CellID (%ld, %ld, %ld) - likely outside local+ghost range.\n",rank, PID_arr[p] , i, j, k);
	}
    }
    LOG_ALLOW(LOCAL, LOG_DEBUG, "CalculateParticleCountPerCell (Rank %d): Local counting finished. Processed %d particles locally.\n", rank, particles_counted_locally);

    // --- Restore Access ---
    ierr = DMDAVecRestoreArray(da, countVec, &count_arr_3d); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm,"DMSwarm_CellID", NULL, NULL, (void **)&global_cell_id_arr); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm,"DMSwarm_pid",NULL,NULL,(void **)&PID_arr);CHKERRQ(ierr);

    // --- Assemble Global Vector ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "CalculateParticleCountPerCell: Assembling global ParticleCount vector.\n");
    ierr = VecAssemblyBegin(countVec); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(countVec); CHKERRQ(ierr);

    // --- Verification Logging ---
    PetscReal total_counted_particles = 0.0, max_count_in_cell = 0.0;
    ierr = VecSum(countVec, &total_counted_particles); CHKERRQ(ierr);
    PetscInt max_idx_global = -1;
    ierr = VecMax(countVec, &max_idx_global, &max_count_in_cell); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "CalculateParticleCountPerCell: Total counted globally = %.0f, Max count in cell = %.0f\n",
              total_counted_particles, max_count_in_cell);

    // --- ADD THIS DEBUGGING BLOCK ---
    if (max_idx_global >= 0) { // Check if VecMax found a location
         // Need to convert the flat global index back to 3D global index (I, J, K)
         // Get global grid dimensions (Nodes, NOT Cells IM/JM/KM)
         PetscInt M, N, P;
         ierr = DMDAGetInfo(da, NULL, &M, &N, &P, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
         // Note: Assuming DOF=1 for countVec, index mapping uses node dimensions M,N,P from DMDA creation (IM+1, etc)
         // Re-check if your DMDA uses cell counts (IM) or node counts (IM+1) for Vec layout. Let's assume Node counts M,N,P.
         PetscInt Kmax = max_idx_global / (M * N);
         PetscInt Jmax = (max_idx_global % (M * N)) / M;
         PetscInt Imax = max_idx_global % M;
         LOG_ALLOW(GLOBAL, LOG_INFO, "  -> Max count located at global index (I,J,K) = (%d, %d, %d) [Flat index: %d]\n",
                   (int)Imax, (int)Jmax, (int)Kmax, (int)max_idx_global);

        // Also, let's explicitly check the count at (0,0,0)
        PetscScalar count_at_origin = 0.0;
        PetscScalar ***count_arr_for_check;
        ierr = DMDAVecGetArrayRead(da, countVec, &count_arr_for_check); CHKERRQ(ierr);
        // Check bounds before accessing - crucial if using global indices
        PetscInt xs, ys, zs, xm, ym, zm;
        ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm); CHKERRQ(ierr);
        if (0 >= xs && 0 < xs+xm && 0 >= ys && 0 < ys+ym && 0 >= zs && 0 < zs+zm) {
             count_at_origin = count_arr_for_check[0][0][0]; // Access using global index (0,0,0)
        } else {
            // Origin is not on this rank (relevant for parallel, but check anyway)
            count_at_origin = -999.0; // Indicate it wasn't accessible locally
        }
        ierr = DMDAVecRestoreArrayRead(da, countVec, &count_arr_for_check); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_INFO, "  -> Count at global index (0,0,0) = %.1f\n", count_at_origin);

    } else {
         LOG_ALLOW(GLOBAL, LOG_WARNING, "  -> VecMax did not return a location for the maximum value.\n");
    }
    // --- END DEBUGGING BLOCK ---
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "CalculateParticleCountPerCell: Particle counting complete.\n");

    PetscFunctionReturn(0);
}

// NOTE: AccumulateParticleField also needs the same modification if it is to handle global CellIDs
// and use DMDAVecGetArray/DMDAVecGetArrayDOF

// --- Helper function to resize swarm globally (add or remove) ---
// This assumes removing excess particles means removing the globally last ones.
PetscErrorCode ResizeSwarmGlobally(DM swarm, PetscInt N_target)
{
    PetscErrorCode ierr;
    PetscInt       N_current, nlocal_current;
    PetscMPIInt    rank;
    MPI_Comm       comm;

    PetscFunctionBeginUser;
    ierr = PetscObjectGetComm((PetscObject)swarm, &comm); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);
    ierr = DMSwarmGetSize(swarm, &N_current); CHKERRQ(ierr);
    ierr = DMSwarmGetLocalSize(swarm, &nlocal_current); CHKERRQ(ierr);

    PetscInt delta = N_target - N_current;

    if (delta == 0) {
        PetscFunctionReturn(0); // Nothing to do
    }

    if (delta < 0) { // Remove particles
        PetscInt num_to_remove_global = -delta;
        LOG_ALLOW(GLOBAL, LOG_INFO, "ResizeSwarmGlobally: Current size %d > target size %d. Removing %d particles globally.\n", N_current, N_target, num_to_remove_global);

        // --- Strategy: Remove the globally last 'num_to_remove_global' particles ---
        // Each rank needs to determine how many of its *local* particles fall
        // within the range of global indices [N_target, N_current - 1].

        PetscInt rstart = 0;
	PetscInt  rend;
	// Global range owned by this rank [rstart, rend)

	ierr = MPI_Exscan(&nlocal_current, &rstart, 1, MPIU_INT, MPI_SUM, comm); CHKERRMPI(ierr); // Use CHKERRMPI for MPI calls

	rend = rstart + nlocal_current;
	

        // Calculate how many local particles have global indices >= N_target
        PetscInt nlocal_remove_count = 0;
        if (rend > N_target) { // If this rank owns any particles slated for removal
            PetscInt start_remove_local_idx = (N_target > rstart) ? (N_target - rstart) : 0;
            nlocal_remove_count = nlocal_current - start_remove_local_idx;
        }

        if (nlocal_remove_count < 0) nlocal_remove_count = 0; // Sanity check
        if (nlocal_remove_count > nlocal_current) nlocal_remove_count = nlocal_current; // Sanity check

        LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Global range [%d, %d). Target size %d. Need to remove %d local particles (from end).\n", rank, rstart, rend, N_target, nlocal_remove_count);

        // Remove the last 'nlocal_remove_count' particles *locally* by iterating backwards
        PetscInt removal_ops_done = 0;
        for (PetscInt p = nlocal_current - 1; p >= 0 && removal_ops_done < nlocal_remove_count; --p) {
            ierr = DMSwarmRemovePointAtIndex(swarm, p); CHKERRQ(ierr);
            removal_ops_done++;
        }

        if (removal_ops_done != nlocal_remove_count) {
             SETERRQ(comm, PETSC_ERR_PLIB, "Rank %d: Failed to remove the expected number of local particles (%d != %d)", rank, removal_ops_done, nlocal_remove_count);
        }
         LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Removed %d local particles.\n", rank, removal_ops_done);
 
	// Barrier to ensure all removals are done before size check
        ierr = MPI_Barrier(comm); CHKERRMPI(ierr);

    } else { // delta > 0: Add particles
        PetscInt num_to_add_global = delta;
        LOG_ALLOW(GLOBAL, LOG_INFO, "ResizeSwarmGlobally: Current size %d < target size %d. Adding %d particles globally.\n", N_current, N_target, num_to_add_global);
        ierr = DMSwarmAddNPoints(swarm, num_to_add_global); CHKERRQ(ierr);
        // Note: Added particles will have uninitialized field data. Reading will overwrite.
    }

    // Verify final size
    PetscInt N_final;
    ierr = DMSwarmGetSize(swarm, &N_final); CHKERRQ(ierr);
    if (N_final != N_target) {
        SETERRQ(comm, PETSC_ERR_PLIB, "Failed to resize swarm: expected %d particles, got %d", N_target, N_final);
    }
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "ResizeSwarmGlobally: Swarm successfully resized to %d particles.\n", N_final);
    PetscFunctionReturn(0);
}

/**
 * @brief Checks particle count from a saved file and resizes the swarm globally.
 *
 * This function uses a robust parallel pattern: only Rank 0 reads the reference
 * position file to determine the total number of particles saved (`N_file`).
 * This count is then broadcast to all other ranks. Finally, each rank compares
 * N_file with the current swarm size and participates in resizing if necessary.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing the DMSwarm.
 * @param[in]     ti   Time index for constructing the file name.
 * @param[in]     ext  File extension (e.g., "dat").
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode PreCheckAndResizeSwarm(UserCtx *user,
                                      PetscInt ti,
                                      const char *ext)
{
    PetscErrorCode ierr;
    char           filename[PETSC_MAX_PATH_LEN];
    PetscInt       N_file = 0; // The number of particles determined from the file
    PetscInt       N_current = 0;
    MPI_Comm       comm;
    PetscMPIInt    rank;
    const char    *refFieldName = "position";
    const PetscInt bs = 3;
    
    // NOTE: Your filename format has a hardcoded "_0" which is typical for
    // PETSc when writing a parallel object from a single rank.
    // If you ever write in parallel, PETSc might create one file per rank.
    // The current logic assumes a single file written by one process.
    const int      placeholder_int = 0;

    PetscFunctionBeginUser;
    ierr = PetscObjectGetComm((PetscObject)user->swarm, &comm); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);

        // --- Construct filename using the specified format ---
    // results/%s%05<PetscInt_FMT>_%d.%s
    ierr = PetscSNPrintf(filename, sizeof(filename), "results/%s%05" PetscInt_FMT "_%d.%s",
                         refFieldName, ti, placeholder_int, ext); CHKERRQ(ierr);
    // Note: Make sure the "results" directory exists or handle directory creation elsewhere.

    LOG_ALLOW(GLOBAL, LOG_INFO, "PreCheckAndResizeSwarm: Checking particle count for timestep %d using ref file '%s'.\n", ti, filename);

    // --- Rank 0 reads the file to determine the size ---
    if (rank == 0) {
        PetscBool fileExists = PETSC_FALSE;
        ierr = PetscTestFile(filename, 'r', &fileExists); CHKERRQ(ierr);

        if (!fileExists) {
            // Set a special value to indicate file not found, then broadcast it.
            N_file = -1;
            LOG_ALLOW(GLOBAL, LOG_ERROR, "Rank 0: Mandatory reference file '%s' not found for timestep %d.\n", filename, ti);
        } else {
            PetscViewer viewer;
            Vec         tmpVec;
            PetscInt    vecSize;
            
            ierr = VecCreate(PETSC_COMM_SELF, &tmpVec); CHKERRQ(ierr); // Create a SEQUENTIAL vector
            ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF, filename, FILE_MODE_READ, &viewer); CHKERRQ(ierr);
            ierr = VecLoad(tmpVec, viewer); CHKERRQ(ierr);
            ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

            ierr = VecGetSize(tmpVec, &vecSize); CHKERRQ(ierr);
            ierr = VecDestroy(&tmpVec); CHKERRQ(ierr);

            if (vecSize % bs != 0) {
                N_file = -2; // Special error code for bad file format
                LOG_ALLOW(GLOBAL, LOG_ERROR, "Rank 0: Vector size %d from file '%s' is not divisible by block size %d.\n", vecSize, filename, bs);
            } else {
                N_file = vecSize / bs;
                LOG_ALLOW(GLOBAL, LOG_DEBUG, "Rank 0: Found %d particles in file.\n", N_file);
            }
        }
    }

    // --- Broadcast the particle count (or error code) from Rank 0 to all other ranks ---
    ierr = MPI_Bcast(&N_file, 1, MPIU_INT, 0, comm); CHKERRMPI(ierr);

    // --- All ranks check for errors and abort if necessary ---
    if (N_file == -1) {
        SETERRQ(comm, PETSC_ERR_FILE_OPEN, "Mandatory reference file '%s' not found for timestep %d (as determined by Rank 0).", filename, ti);
    }
    if (N_file == -2) {
        SETERRQ(comm, PETSC_ERR_FILE_READ, "Reference file '%s' has incorrect format (as determined by Rank 0).", filename);
    }
    if (N_file < 0) {
         SETERRQ(comm, PETSC_ERR_PLIB, "Received invalid particle count %d from Rank 0.", N_file);
    }


    // --- Now all ranks have the correct N_file, compare and resize if needed ---
    ierr = DMSwarmGetSize(user->swarm, &N_current); CHKERRQ(ierr);

    if (N_file != N_current) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Swarm size %d differs from file size %d. Resizing swarm globally.\n", N_current, N_file);
        ierr = ResizeSwarmGlobally(user->swarm, N_file); CHKERRQ(ierr);
    } else {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Swarm size (%d) already matches file size. No resize needed.\n", N_current);
    }
    
    // Also update the context
    user->NumberofParticles = N_file;

    PetscFunctionReturn(0);
}


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
                                                   PetscInt *globalMigrationCount_out)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank; // MPI rank of the current process

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

        /******************************************************************/
    /*                 START OF NEW TEST LOGIC                        */
    /******************************************************************/

    // --- TEST STEP 1: Take a snapshot of PIDs BEFORE migration ---
    PetscInt nlocal_before;
    ierr = DMSwarmGetLocalSize(user->swarm, &nlocal_before); CHKERRQ(ierr);
    
    PetscInt64 *pids_before_snapshot = NULL;
    // We need to get the PID field to pass to our snapshot function
    PetscInt64 *pid_field_for_snapshot;
    ierr = DMSwarmGetField(user->swarm, "DMSwarm_pid", NULL, NULL, (void**)&pid_field_for_snapshot); CHKERRQ(ierr);
    
    // Call our helper to create the sorted snapshot
    ierr = GetLocalPIDSnapshot(pid_field_for_snapshot, nlocal_before, &pids_before_snapshot); CHKERRQ(ierr);

    // Restore the field immediately
    ierr = DMSwarmRestoreField(user->swarm, "DMSwarm_pid", NULL, NULL, (void**)&pid_field_for_snapshot); CHKERRQ(ierr);
    
    LOG_ALLOW(LOCAL, LOG_INFO, "[TEST HARNESS - Rank %d] Created pre-migration PID snapshot with %d entries.\n", rank, nlocal_before);

    /******************************************************************/
    /*                 END OF NEW TEST LOGIC (PART 1)                 */
    /******************************************************************/

    // Step 1: Identify particles that need to migrate from the current rank
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] %s Migration: Identifying migrating particles.\n", currentTime, step, migrationCycleName);
    ierr = IdentifyMigratingParticles(user, bboxlist, migrationList_p, migrationCount_p, migrationListCapacity_p); CHKERRQ(ierr);

    // Ensure Identification is done for all ranks before sharing.
    ierr = PetscBarrier(NULL);
    
    // Step 2: Get the global count of migrating particles
    PetscInt localMigrationCount = *migrationCount_p; // Use a local variable for MPI_Allreduce
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "[Rank %d] Before MPI_Allreduce. localMigrationCount = %d.\n", rank, localMigrationCount);
    ierr = MPI_Allreduce(&localMigrationCount, globalMigrationCount_out, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "[Rank %d] After MPI_Allreduce. globalMigrationCount_out = %d.\n", rank, *globalMigrationCount_out);
     
    // Step 3: Perform migration if any particles are moving globally
    if (*globalMigrationCount_out > 0) {
        LOG_ALLOW(LOCAL, LOG_INFO, "[T=%.4f, Step=%d] %s Migration: Performing migration (%d particles globally, %d locally from rank %d).\n",
                  currentTime, step, migrationCycleName, *globalMigrationCount_out, localMigrationCount, rank);
        // Set the destination ranks for the locally identified migrating particles
        ierr = SetMigrationRanks(user, *migrationList_p, localMigrationCount); CHKERRQ(ierr);
        // Execute the migration
        ierr = PerformMigration(user); CHKERRQ(ierr);

	// Make sure all ranks finish Migration before proceding.
	ierr = PetscBarrier(NULL);
    } else {
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] %s Migration: No particles identified globally for migration.\n", currentTime, step, migrationCycleName);
    }

        /******************************************************************/
    /*                 START OF NEW TEST LOGIC (PART 2)               */
    /******************************************************************/
    
    // --- TEST STEP 2: Call the function we want to test ---
    // This happens AFTER the migration is complete.
    if (*globalMigrationCount_out > 0) {
        LOG_ALLOW(LOCAL, LOG_INFO, "[TEST HARNESS - Rank %d] Calling FlagNewcomersForLocation to verify identification.\n", rank);
        ierr = FlagNewcomersForLocation(user->swarm, nlocal_before, pids_before_snapshot); CHKERRQ(ierr);
    }

    // --- TEST STEP 3: Cleanup ---
    ierr = PetscFree(pids_before_snapshot); CHKERRQ(ierr);

    /******************************************************************/
    /*                 END OF NEW TEST LOGIC (PART 2)                 */
    /******************************************************************/
    
    // Reset local migration count for the next potential migration cycle by the caller.
    // The migrationList and its capacity persist and are managed by the caller.
    *migrationCount_p = 0;

    PetscFunctionReturn(0);
}

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
PetscErrorCode ReinitializeParticlesOnInletSurface(UserCtx *user, PetscReal currentTime, PetscInt step)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;                        // MPI rank of the current process
    DM             swarm = user->swarm;         // The particle swarm DM
    PetscReal      *positions_field = NULL;     // Pointer to swarm field for physical positions
    PetscReal      *pos_phy_field = NULL;       // Pointer to swarm field for physical positions (backup)
    PetscInt64     *particleIDs = NULL;         // Pointer to swarm field for Particle IDs (for logging)
    PetscInt       *cell_ID_field = NULL;       // Pointer to swarm field for Cell IDs (for resetting after migration)
    const Cmpnts   ***coor_nodes_local_array;   // Read-only access to local node coordinates
    Vec            Coor_local;                  // Local vector for node coordinates
    DMDALocalInfo  info;                        // Local grid information (node-based) from user->da
    PetscInt       xs_gnode_rank, ys_gnode_rank, zs_gnode_rank; // Local starting node indices (incl. ghosts) of rank's DA
    PetscInt       IM_nodes_global, JM_nodes_global, KM_nodes_global; // Global node counts

    PetscRandom    rand_logic_reinit_i, rand_logic_reinit_j, rand_logic_reinit_k; // RNGs for re-placement
    PetscInt       nlocal_current;                // Number of particles currently on this rank
    PetscInt       particles_actually_reinitialized_count = 0; // Counter for logging
    PetscBool      can_this_rank_service_inlet = PETSC_FALSE;  // Flag

    PetscFunctionBeginUser;

    // This function is only relevant for surface initialization mode and if an inlet face is defined.
    if (user->ParticleInitialization != 0 || !user->inletFaceDefined) {
        PetscFunctionReturn(0);
    }

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = DMSwarmGetLocalSize(swarm, &nlocal_current); CHKERRQ(ierr);

    // If no particles on this rank, nothing to do.
    if (nlocal_current == 0) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f, Step=%d] Rank %d has no local particles to re-initialize on inlet.\n", currentTime, step, rank);
        PetscFunctionReturn(0);
    }

    // Get DMDA information for the node-centered coordinate grid (user->da)
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    ierr = DMDAGetInfo(user->da, NULL, &IM_nodes_global, &JM_nodes_global, &KM_nodes_global, NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL); CHKERRQ(ierr);
    ierr = DMDAGetGhostCorners(user->da, &xs_gnode_rank, &ys_gnode_rank, &zs_gnode_rank, NULL, NULL, NULL); CHKERRQ(ierr);

    // Check if this rank is responsible for (part of) the designated inlet surface
    ierr = CanRankServiceInletFace(user, &info, IM_nodes_global, JM_nodes_global, KM_nodes_global, &can_this_rank_service_inlet); CHKERRQ(ierr);

    if (!can_this_rank_service_inlet) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f, Step=%d] Rank %d cannot service inlet face %d. Skipping re-initialization of %d particles.\n", currentTime, step, rank, user->identifiedInletBCFace, nlocal_current);
        PetscFunctionReturn(0);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Rank %d is on inlet face %d. Attempting to re-place %d local particles.\n", currentTime, step, rank, user->identifiedInletBCFace, nlocal_current);

    // Get coordinate array and swarm fields for modification
    ierr = DMGetCoordinatesLocal(user->da, &Coor_local); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, Coor_local, (void*)&coor_nodes_local_array); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "position",       NULL, NULL, (void**)&positions_field); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "pos_phy",        NULL, NULL, (void**)&pos_phy_field);   CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_pid",    NULL, NULL, (void**)&particleIDs);    CHKERRQ(ierr); // For logging
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cell_ID_field);  CHKERRQ(ierr);

    // Initialize fresh RNGs for this re-placement to ensure good distribution
    ierr = InitializeLogicalSpaceRNGs(&rand_logic_reinit_i, &rand_logic_reinit_j, &rand_logic_reinit_k); CHKERRQ(ierr);
    // Optional: Seed RNGs for deterministic behavior if required, e.g., based on rank and step.
    // PetscRandomSetSeed(rand_logic_i, (unsigned long)rank*1000 + step + 100); PetscRandomSeed(rand_logic_i); // Example

    // Loop over all particles currently local to this rank
    for (PetscInt p = 0; p < nlocal_current; p++) {
        PetscInt  ci_metric_lnode, cj_metric_lnode, ck_metric_lnode; // Local node indices (of rank's DA patch) for cell origin
        PetscReal xi_metric_logic, eta_metric_logic, zta_metric_logic; // Intra-cell logical coordinates
        Cmpnts    phys_coords; // To store newly calculated physical coordinates

        // Get random cell on this rank's portion of the inlet and random logical coords within it
        ierr = GetRandomCellAndLogicOnInletFace(user, &info, xs_gnode_rank, ys_gnode_rank, zs_gnode_rank,
                                                IM_nodes_global, JM_nodes_global, KM_nodes_global,
                                                &rand_logic_reinit_i, &rand_logic_reinit_j, &rand_logic_reinit_k,
                                                &ci_metric_lnode, &cj_metric_lnode, &ck_metric_lnode,
                                                &xi_metric_logic, &eta_metric_logic, &zta_metric_logic); CHKERRQ(ierr);
        
        // Convert these logical coordinates to physical coordinates
        ierr = MetricLogicalToPhysical(user, coor_nodes_local_array,
                                       ci_metric_lnode, cj_metric_lnode, ck_metric_lnode,
                                       xi_metric_logic, eta_metric_logic, zta_metric_logic,
                                       &phys_coords); CHKERRQ(ierr);

        // Update the particle's position in the swarm fields
        positions_field[3*p+0] = phys_coords.x; 
        positions_field[3*p+1] = phys_coords.y; 
        positions_field[3*p+2] = phys_coords.z;
        pos_phy_field[3*p+0]   = phys_coords.x; 
        pos_phy_field[3*p+1]   = phys_coords.y; 
        pos_phy_field[3*p+2]   = phys_coords.z;
        particles_actually_reinitialized_count++;

	cell_ID_field[3*p+0] = -1;
	cell_ID_field[3*p+1] = -1;
	cell_ID_field[3*p+2] = -1;

        LOG_LOOP_ALLOW(LOCAL, LOG_DEBUG, p, (nlocal_current > 20 ? nlocal_current/10 : 1), // Sampled logging
            "ReInit - Rank %d: PID %ld (idx %ld) RE-PLACED. CellOriginNode(locDAIdx):(%d,%d,%d). LogicCoords: (%.2e,%.2f,%.2f). PhysCoords: (%.6f,%.6f,%.6f).\n",
            rank, particleIDs[p], (long)p, 
            ci_metric_lnode, cj_metric_lnode, ck_metric_lnode,
            xi_metric_logic, eta_metric_logic, zta_metric_logic, 
            phys_coords.x, phys_coords.y, phys_coords.z);
    }

    // Logging summary of re-initialization
    if (particles_actually_reinitialized_count > 0) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Rank %d (on inlet face %d) successfully re-initialized %d of %d local particles.\n", currentTime, step, rank, user->identifiedInletBCFace, particles_actually_reinitialized_count, nlocal_current);
    } else if (nlocal_current > 0) { // This case should ideally not be hit if can_this_rank_service_inlet was true and particles were present.
         LOG_ALLOW(GLOBAL, LOG_WARNING, "[T=%.4f, Step=%d] Rank %d claimed to service inlet face %d, but re-initialized 0 of %d local particles. This may indicate an issue if particles were expected to be re-placed.\n", currentTime, step, rank, user->identifiedInletBCFace, nlocal_current);
    }

    // Cleanup: Destroy RNGs and restore swarm fields/coordinate array
    ierr = PetscRandomDestroy(&rand_logic_reinit_i); CHKERRQ(ierr);
    ierr = PetscRandomDestroy(&rand_logic_reinit_j); CHKERRQ(ierr);
    ierr = PetscRandomDestroy(&rand_logic_reinit_k); CHKERRQ(ierr);

    ierr = DMSwarmRestoreField(swarm, "position",       NULL, NULL, (void**)&positions_field); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "pos_phy",        NULL, NULL, (void**)&pos_phy_field);   CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid",    NULL, NULL, (void**)&particleIDs);    CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cell_ID_field);  CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, Coor_local, (void*)&coor_nodes_local_array); CHKERRQ(ierr);


    PetscFunctionReturn(0);
}

// Comparison function needed for qsort
static int compare_PetscInt64(const void *a, const void *b) {
    PetscInt64 val_a = *(const PetscInt64*)a;
    PetscInt64 val_b = *(const PetscInt64*)b;
    if (val_a < val_b) return -1;
    if (val_a > val_b) return 1;
    return 0;
}

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
                                   PetscInt64 **pids_snapshot_out)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;

    PetscFunctionBeginUser;

    // --- 1. Input Validation ---
    if (!pids_snapshot_out) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Output pointer pids_snapshot_out is NULL.");
    }
    // If n_local > 0, pid_field must not be NULL.
    if (n_local > 0 && !pid_field) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input pid_field pointer is NULL for n_local > 0.");
    }
    
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d]: Creating PID snapshot for %d local particles.\n", rank, n_local);

    // If there are no local particles, the snapshot is empty (NULL).
    if (n_local == 0) {
        *pids_snapshot_out = NULL;
        PetscFunctionReturn(0);
    }
    
    // --- 2. Allocate Memory for the Snapshot ---
    ierr = PetscMalloc1(n_local, pids_snapshot_out); CHKERRQ(ierr);

    // --- 3. Copy Data ---
    // Perform a fast memory copy from the provided array to our new snapshot array.
    ierr = PetscMemcpy(*pids_snapshot_out, pid_field, n_local * sizeof(PetscInt64)); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d]: Copied %d PIDs.\n", rank, n_local);

    // --- 4. Sort the Snapshot Array ---
    // Sorting enables fast binary search lookups later.
    ierr = PetscSortInt64(n_local, *pids_snapshot_out); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d]: PID snapshot sorted successfully.\n", rank);

    PetscFunctionReturn(0); 
}


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
                                  PetscMPIInt destination_rank)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;

    PetscFunctionBeginUser;

    // --- 1. Input Validation ---
    if (!migration_list_p || !capacity_p || !count_p) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Null pointer provided to AddToMigrationList for list management.");
    }

    // --- 2. Check if the list needs to be resized ---
    if (*count_p >= *capacity_p) {
        PetscInt old_capacity = *capacity_p;
        // Start with a reasonable base capacity, then double for subsequent reallocations.
        PetscInt new_capacity = (old_capacity == 0) ? 16 : old_capacity * 2;

        // Use PetscRealloc for safe memory reallocation.
        // It handles allocating new memory, copying old data, and freeing the old block.
        // The first argument to PetscRealloc is the new size in BYTES.
        ierr = PetscRealloc(new_capacity * sizeof(MigrationInfo), migration_list_p); CHKERRQ(ierr);
        
        *capacity_p = new_capacity; // Update the capacity tracker

        ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
        LOG_ALLOW(LOCAL, LOG_DEBUG, "AddToMigrationList [Rank %d]: Reallocated migrationList capacity from %d to %d.\n",
                  rank, old_capacity, new_capacity);
    }

    // --- 3. Add the new migration data to the list ---
    // Dereference the pointer-to-a-pointer to get the actual array.
    MigrationInfo *list = *migration_list_p;

    list[*count_p].local_index = particle_local_idx;
    list[*count_p].target_rank = destination_rank;

    // --- 4. Increment the count of items in the list ---
    (*count_p)++;
    
    PetscFunctionReturn(0);
}

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
                                        const PetscInt64 pids_before[])
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    PetscInt       n_local_after;
    PetscInt       newcomer_count = 0;
    
    // Pointers to the swarm data fields we will read and modify
    PetscInt64 *pid_field_after    = NULL;
    PetscInt   *status_field_after = NULL;
    PetscInt   *cell_field_after = NULL;

    PetscFunctionBeginUser;

    // --- 1. Input Validation and Basic Setup ---
    if (!swarm) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Input DMSwarm is NULL in FlagNewcomersForLocation.");
    }
    // If n_local_before > 0, the corresponding PID array must not be null.
    if (n_local_before > 0 && !pids_before) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input pids_before array is NULL for n_local_before > 0.");
    }

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    
    // Get the number of particles on this rank *after* the migration.
    ierr = DMSwarmGetLocalSize(swarm, &n_local_after); CHKERRQ(ierr);
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "FlagNewcomersForLocation [Rank %d]: Checking for newcomers. Size before: %d, Size after: %d\n",
              rank, n_local_before, n_local_after);

    // If there are no particles now, there's nothing to do.
    if (n_local_after == 0) {
        PetscFunctionReturn(0);
    }
    
    // --- 2. Access Swarm Data ---
    // Get read-only access to the PIDs and read-write access to the status field.
    ierr = DMSwarmGetField(swarm, "DMSwarm_pid",             NULL, NULL, (void**)&pid_field_after); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_location_status", NULL, NULL, (void**)&status_field_after); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID",          NULL, NULL, (void**)&cell_field_after); CHKERRQ(ierr);
    if (!pid_field_after || !status_field_after) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Failed to get required swarm fields in FlagNewcomersForLocation.");
    }

    // --- 3. Identify and Flag Newcomers ---
    // Loop through all particles currently on this rank.
    for (PetscInt p_idx = 0; p_idx < n_local_after; ++p_idx) {
        PetscInt64 current_pid = pid_field_after[p_idx];
        PetscBool  is_found_in_before_list;

        // Use our custom, efficient helper function for the lookup.
        ierr = BinarySearchInt64(n_local_before, pids_before, current_pid, &is_found_in_before_list); CHKERRQ(ierr);

        // If the PID was NOT found in the "before" list, it must be a newcomer.
        if (!is_found_in_before_list) {
          // Flag it for processing in the next pass of the migration loop.
          status_field_after[p_idx] = NEEDS_LOCATION;
	  // cell_field_after[3*p_idx+0] = -1;
	  // cell_field_after[3*p_idx+1] = -1; 
	  // cell_field_after[3*p_idx+2] = -1; 
            newcomer_count++;
            
            LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d]: Flagged newcomer PID %ld at local index %d as NEEDS_LOCATION.\n",
                      rank, current_pid, p_idx);
        }
    }

    // --- 4. Restore Swarm Fields ---
    // Release the locks on the swarm data arrays.
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid",             NULL, NULL, (void**)&pid_field_after); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_location_status", NULL, NULL, (void**)&status_field_after); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID",          NULL, NULL, (void**)&cell_field_after); CHKERRQ(ierr);

    if (newcomer_count > 0) {
        LOG_ALLOW(LOCAL, LOG_INFO, "[Rank %d]: Identified and flagged %d newcomers.\n", rank, newcomer_count);
    }

    PetscFunctionReturn(0);
}

/**
 * @brief Provides a fast, heuristic-based guess for a particle's owner rank using bounding boxes.
 * @ingroup ParticleLocation
 *
 * This function is part of the "Guess and Verify" strategy, called only for "lost"
 * particles. It attempts to find a candidate owner by checking which rank's bounding box
 * contains the particle's physical position.
 *
 * To optimize the search, it uses the particle's position relative to the local
 * bounding box to intelligently check the most likely neighboring ranks first.
 * For example, if a particle's x-coordinate is less than the local minimum x, it
 * will check the -X neighbor first. If no owner is found in the immediate neighbors,
 * it performs a full search of all other ranks as a fallback.
 *
 * @param[in]  user             Pointer to the UserCtx, which must contain the pre-computed
 *                              `bbox` (local), `neighbors` struct, and the global `bboxlist`.
 * @param[in]  particle         A pointer to the particle whose owner is being guessed.
 * @param[in]  bboxlist       An array of BoundingBox structures for ALL MPI ranks, indexed 0 to (size-1).
 *                                This array must be up-to-date and available on all ranks.
 * @param[out] guess_rank_out   A pointer to a PetscMPIInt. Set to the candidate owner's rank
 *                              if found, otherwise set to -1 (or MPI_PROC_NULL).
 *
 * @return PetscErrorCode 0 on success, or a non-zero PETSc error code on failure.
 */
PetscErrorCode GuessParticleOwnerWithBBox(UserCtx *user, 
                                          const Particle *particle,
					  const BoundingBox *bboxlist,
                                          PetscMPIInt *guess_rank_out)
{
    PetscErrorCode  ierr;
    PetscMPIInt     rank, size;
    const RankNeighbors *neighbors = &user->neighbors; // Use a direct pointer for clarity
    const BoundingBox   *localBBox = &user->bbox;

    PetscFunctionBeginUser;

    // --- 1. Input Validation and Setup ---
    if (!user || !particle || !guess_rank_out || !bboxlist) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Null pointer provided to GuessParticleOwnerWithBBox.");
    }
    if (!localBBox|| !neighbors) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Required user->bboxl or user->neighbors is not initialized.");
    }

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
    
    *guess_rank_out = MPI_PROC_NULL; // Default to "not found"

    LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %ld]: Starting guess for particle at (%.3f, %.3f, %.3f).\n",
              particle->PID, particle->loc.x, particle->loc.y, particle->loc.z);

    // *** THE PRIMARY FIX ***
    // --- Step 0: Check if the particle is inside the CURRENT rank's bounding box FIRST. ---
    // This handles the common case of initial placement where a particle is "lost" but physically local.
    if (IsParticleInBox(localBBox, &particle->loc)) {
      *guess_rank_out = rank;
      LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %ld]: Fast path guess SUCCESS. Particle is within the local (Rank %d) bounding box.\n",
		particle->PID, rank);
      PetscFunctionReturn(0); // Found it, we're done.
    }
    // --- 2. Fast Path: Check Immediate Neighbors Based on Exit Direction ---
    // This is the logic repurposed directly from your IdentifyMigratingParticles.

    // Determine likely exit direction(s) to prioritize neighbor check
    PetscBool exit_xm = particle->loc.x < localBBox->min_coords.x;
    PetscBool exit_xp = particle->loc.x > localBBox->max_coords.x;
    PetscBool exit_ym = particle->loc.y < localBBox->min_coords.y;
    PetscBool exit_yp = particle->loc.y > localBBox->max_coords.y;
    PetscBool exit_zm = particle->loc.z < localBBox->min_coords.z;
    PetscBool exit_zp = particle->loc.z > localBBox->max_coords.z;
    
    if (exit_xm && neighbors->rank_xm != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors->rank_xm], &particle->loc)) {
        *guess_rank_out = neighbors->rank_xm;
    } else if (exit_xp&& neighbors->rank_xp != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors->rank_xp], &particle->loc)) {
        *guess_rank_out = neighbors->rank_xp;
    } else if (exit_ym && neighbors->rank_ym != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors->rank_ym], &particle->loc)) {
        *guess_rank_out = neighbors->rank_ym;
    } else if (exit_yp && neighbors->rank_yp != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors->rank_yp], &particle->loc)) {
        *guess_rank_out = neighbors->rank_yp;
    } else if (exit_zm && neighbors->rank_zm != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors->rank_zm], &particle->loc)) {
        *guess_rank_out = neighbors->rank_zm;
    } else if (exit_zp && neighbors->rank_zp != MPI_PROC_NULL && IsParticleInBox(&bboxlist[neighbors->rank_zp], &particle->loc)) {
        *guess_rank_out = neighbors->rank_zp;
    }
    // Note: This does not handle corner/edge neighbors, which is why the fallback is essential.

    if (*guess_rank_out != MPI_PROC_NULL) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %ld]: Fast path guess SUCCESS. Found in immediate neighbor Rank %d.\n",
                  particle->PID, *guess_rank_out);
        PetscFunctionReturn(0); // Found it, we're done.
    }

    // --- 3. Robust Fallback: Check All Other Ranks ---
    // If we get here, the particle was not in any of the immediate face neighbors' boxes.
    LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %ld]: Not in immediate face neighbors. Starting global fallback search.\n",
              particle->PID);
    
    for (PetscMPIInt r = 0; r < size; ++r) {
        if (r == rank) continue; // Don't check ourselves.

        if (IsParticleInBox(&bboxlist[r], &particle->loc)) {
          PetscBool is_in = PETSC_TRUE;
	  // This detailed, synchronized print will solve the mystery
	  LOG_ALLOW(LOCAL,LOG_DEBUG, "[Guess BBox | Rank %d] Checking PID %lld at (%.4f, %.4f, %.4f) against Rank %d's box: [(%.4f, %.4f, %.4f) to (%.4f, %.4f, %.4f)]. Result: %s\n",
				  (int)rank, (long long)particle->PID,
				  particle->loc.x, particle->loc.y, particle->loc.z,
				  (int)r,
				  bboxlist[r].min_coords.x, bboxlist[r].min_coords.y, bboxlist[r].min_coords.z,
				  bboxlist[r].max_coords.x, bboxlist[r].max_coords.y, bboxlist[r].max_coords.z,
				  is_in ? "INSIDE" : "OUTSIDE");
	  
            *guess_rank_out = r;
            LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %ld]: Fallback search SUCCESS. Found in Rank %d.\n",
                      particle->PID, *guess_rank_out);
            PetscFunctionReturn(0); // Found it, we're done.
        }
    }

    // If the code reaches here, the particle was not found in any rank's bounding box.
    LOG_ALLOW(LOCAL, LOG_WARNING, "[PID %ld]: Guess FAILED. Particle not found in any rank's bounding box.\n",
              particle->PID);
    
    // The guess_rank_out will remain -1, signaling failure to the caller.
    PetscFunctionReturn(0);
}

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
PetscErrorCode LocateAllParticlesInGrid_TEST(UserCtx *user,BoundingBox *bboxlist)
{
    PetscErrorCode ierr;
    PetscInt       passes = 0;
    const PetscInt MAX_MIGRATION_PASSES = 10; // Safety break for runaway loops
    PetscInt       global_migrations_this_pass;
    PetscMPIInt    rank;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_FUNC_TIMER_BEGIN_EVENT(EVENT_GlobalParticleLocation, GLOBAL);
    LOG_ALLOW(GLOBAL, LOG_INFO, "LocateAllParticlesInGrid (Orchestrator) - Beginning particle settlement process.\n");

    // This loop ensures that particles that jump across multiple ranks are
    // handled correctly in successive, iterative handoffs.
    do {
        passes++;
        LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, "[Rank %d] Starting migration pass %d.\n", rank, passes);

        // --- STAGE 1: PER-PASS INITIALIZATION ---
        MigrationInfo  *migrationList = NULL;
        PetscInt       local_migration_count = 0;
        PetscInt       migrationListCapacity = 0;
        PetscInt       nlocal_before;
        PetscInt64     *pids_before_snapshot = NULL;
	PetscInt       local_lost_count = 0;

        ierr = DMSwarmGetLocalSize(user->swarm, &nlocal_before); CHKERRQ(ierr);
        LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d] Pass %d begins with %d local particles.\n", rank, passes, nlocal_before);


        // --- STAGE 2: PRE-MIGRATION SNAPSHOT & MAIN PROCESSING LOOP ---
        if (nlocal_before > 0) {
            // Get pointers to all fields needed for this pass
            PetscReal  *pos_p, *weights_p, *vel_p;
            PetscInt   *cell_p, *status_p;
            PetscInt64 *pid_p;
            ierr = DMSwarmGetField(user->swarm, "position",                NULL, NULL, (void**)&pos_p);    CHKERRQ(ierr);
            ierr = DMSwarmGetField(user->swarm, "velocity",                NULL, NULL, (void**)&vel_p);    CHKERRQ(ierr);
            ierr = DMSwarmGetField(user->swarm, "weight",                  NULL, NULL, (void**)&weights_p);  CHKERRQ(ierr);
            ierr = DMSwarmGetField(user->swarm, "DMSwarm_CellID",          NULL, NULL, (void**)&cell_p);     CHKERRQ(ierr);
            ierr = DMSwarmGetField(user->swarm, "DMSwarm_pid",             NULL, NULL, (void**)&pid_p);      CHKERRQ(ierr);
            ierr = DMSwarmGetField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void**)&status_p);   CHKERRQ(ierr);

            // Create a sorted snapshot of current PIDs to identify newcomers after migration.
            // This helper requires a raw pointer, which we just acquired.
            ierr = GetLocalPIDSnapshot(pid_p, nlocal_before, &pids_before_snapshot); CHKERRQ(ierr);

            for (PetscInt p_idx = 0; p_idx < nlocal_before; p_idx++) {

                // OPTIMIZATION: Skip particles already settled in a previous pass of this do-while loop.

	      LOG_ALLOW(LOCAL,LOG_DEBUG,
			"p_idx=%d, PID=%ld, status=%d, cell=(%d, %d, %d)\n",
			p_idx,
			(long)pid_p[p_idx],
			status_p[p_idx],
			cell_p[3*p_idx],
			cell_p[3*p_idx+1],
			cell_p[3*p_idx+2]);

                if (status_p[p_idx] == ACTIVE_AND_LOCATED) {
		  LOG_ALLOW(LOCAL,LOG_DEBUG," [rank %d][PID %ld] skipped in pass %d as it is already located at (%d,%d,%d).\n",rank,pid_p[p_idx],passes,cell_p[3*p_idx],cell_p[3*p_idx + 1],cell_p[3*p_idx + 2]);
                    continue;
                }

                // UNPACK: Create a temporary C struct for easier processing using our helper.
                Particle current_particle;

		//	LOG_ALLOW(LOCAL,LOG_DEBUG,"about to unpack p_idx=%d (PID=%ld)\n",p_idx, (long)pid_p[p_idx]);
		
                ierr = UnpackSwarmFields(p_idx, pid_p, weights_p, pos_p, cell_p, vel_p, status_p, &current_particle); CHKERRQ(ierr);

		//		LOG_ALLOW(LOCAL,LOG_DEBUG,"unpacked p_idx=%d → cell[0]=%d, status=%d\n",p_idx, current_particle.cell[0], current_particle.location_status);

		ParticleLocationStatus final_status = (ParticleLocationStatus)status_p[p_idx];


		// CASE 1: Particle has a valid prior cell index.
                // It has moved, so we only need to run the robust walk from its last known location.
                if (current_particle.cell[0] >= 0) {
		  LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %ld] has valid prior cell. Strategy: Robust Walk from previous cell.\n", current_particle.PID);
		  ierr = LocateParticleOrFindMigrationTarget_TEST(user, &current_particle, &final_status); CHKERRQ(ierr);
                } 

		/*		
                // --- "GUESS" FAST PATH for lost particles ---
                if (current_particle.cell[0] < 0) {
                    LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %ld] is lost or uninitialzied (cell=%d), attempting fast guess.\n",current_particle.PID, current_particle.cell[0]);
                    ierr = GuessParticleOwnerWithBBox(user, &current_particle, bboxlist, &destination_rank); CHKERRQ(ierr);
                    if (destination_rank != MPI_PROC_NULL && destination_rank != rank) {
                        final_status = MIGRATING_OUT;
                        // The particle struct's destination rank must be updated for consistency
                        current_particle.destination_rank = destination_rank;
                    }
                }

		LOG_ALLOW(LOCAL,LOG_DEBUG,"[PID %ld] Particle status after Initial Guess:%d \n",current_particle.PID,final_status);

                // --- "VERIFY" ROBUST WALK if guess didn't resolve it ---
                if (final_status == NEEDS_LOCATION  || UNINITIALIZED) {
                    LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %ld] Not resolved by guess, starting robust walk.\n", current_particle.PID);
                    // This function will update the particle's status and destination rank internally.
		     ierr = LocateParticleOrFindMigrationTarget_TEST(user, &current_particle, &final_status); CHKERRQ(ierr);
                    destination_rank = current_particle.destination_rank; // Retrieve the result
                }

                // --- PROCESS THE FINAL STATUS AND TAKE ACTION ---
                if (final_status == MIGRATING_OUT) {
                    status_p[p_idx] = MIGRATING_OUT; // Mark for removal by DMSwarm
                    ierr = AddToMigrationList(&migrationList, &migrationListCapacity, &local_migration_count, p_idx, destination_rank); CHKERRQ(ierr);
                    LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %ld] at local index %d marked for migration to rank %d.\n",current_particle.PID, p_idx, destination_rank);
                } else {
                     // Particle's final status is either LOCATED or LOST; update its state in the swarm arrays.
                     current_particle.location_status = final_status;
                     // PACK: Use the helper to write results back to the swarm arrays.
                     ierr = UpdateSwarmFields(p_idx, &current_particle, weights_p, cell_p, status_p); CHKERRQ(ierr);
                }
		*/
                // CASE 2: Particle is "lost" (cell = -1). Strategy: Guess -> Verify.
                else {
		  LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %ld] has invalid cell. Strategy: Guess Owner -> Find Cell.\n",current_particle.PID);
                    
		  PetscMPIInt guessed_owner_rank = MPI_PROC_NULL;
		  ierr = GuessParticleOwnerWithBBox(user, &current_particle, bboxlist, &guessed_owner_rank); CHKERRQ(ierr);

		  // If the guess finds a DIFFERENT rank, we can mark for migration and skip the walk.
		  if (guessed_owner_rank != MPI_PROC_NULL && guessed_owner_rank != rank) {
		    LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %ld] Guess SUCCESS: Found migration target Rank %d. Finalizing.\n", current_particle.PID, guessed_owner_rank);
		    final_status = MIGRATING_OUT;
		    current_particle.destination_rank = guessed_owner_rank;
		  } 
		  else {

		    // This block runs if the guess either failed (rank is NULL) or found the particle is local (rank is self).
		    // In BOTH cases, the situation is unresolved, and we MUST fall back to the robust walk.
		    if (guessed_owner_rank == rank) {
		      LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %ld] Guess determined particle is local. Proceeding to robust walk to find cell.\n", current_particle.PID);
		    } else { // guessed_owner_rank == MPI_PROC_NULL
		      LOG_ALLOW(LOCAL, LOG_WARNING, "[PID %ld] Guess FAILED to find an owner. Proceeding to robust walk for definitive search.\n", current_particle.PID);
		    }
                        
		    ierr = LocateParticleOrFindMigrationTarget_TEST(user, &current_particle, &final_status); CHKERRQ(ierr);
		  }
                }
		
                // --- PROCESS THE FINAL, DEFINITIVE STATUS ---
                current_particle.location_status = final_status;
                ierr = UpdateSwarmFields(p_idx, &current_particle, weights_p, cell_p, status_p); CHKERRQ(ierr);
                
                if (final_status == MIGRATING_OUT) {
		  ierr = AddToMigrationList(&migrationList, &migrationListCapacity, &local_migration_count, p_idx, current_particle.destination_rank); CHKERRQ(ierr);
                } else if (final_status == LOST) {
		  local_lost_count++;
                }
	    		
            } // End of main particle processing loop

            // Restore all the fields acquired for this pass.
            ierr = DMSwarmRestoreField(user->swarm, "position",                NULL, NULL, (void**)&pos_p);    CHKERRQ(ierr);
            ierr = DMSwarmRestoreField(user->swarm, "velocity",                NULL, NULL, (void**)&vel_p);    CHKERRQ(ierr);
            ierr = DMSwarmRestoreField(user->swarm, "weight",                  NULL, NULL, (void**)&weights_p);  CHKERRQ(ierr);
            ierr = DMSwarmRestoreField(user->swarm, "DMSwarm_CellID",          NULL, NULL, (void**)&cell_p);     CHKERRQ(ierr);
            ierr = DMSwarmRestoreField(user->swarm, "DMSwarm_pid",             NULL, NULL, (void**)&pid_p);      CHKERRQ(ierr);
            ierr = DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void**)&status_p);   CHKERRQ(ierr);
        }

        // --- STAGE 3: ACTION & MPI COMMUNICATION ---
        LOG_ALLOW(LOCAL, LOG_INFO, "[Rank %d] Pass %d: Identified %d particles to migrate out.\n", rank, passes, local_migration_count);

	// --- STAGE 3: SYNCHRONIZE AND DECIDE ---
        // FIRST, determine if any rank wants to migrate. This call is safe because
        // all ranks have finished their local work and can participate.
        ierr = MPI_Allreduce(&local_migration_count, &global_migrations_this_pass, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);	

	if(global_migrations_this_pass > 0 ){

	LOG_ALLOW(GLOBAL, LOG_INFO, "Pass %d: Migrating %d particles globally.\n", passes, global_migrations_this_pass);
        
        ierr = SetMigrationRanks(user, migrationList, local_migration_count); CHKERRQ(ierr);
        ierr = PerformMigration(user); CHKERRQ(ierr);

        // --- STAGE 4: POST-MIGRATION RESET ---
        // Identify newly arrived particles and flag them with NEEDS_LOCATION so they are
        // processed in the next pass. This uses the snapshot taken in STAGE 2.
        ierr = FlagNewcomersForLocation(user->swarm, nlocal_before, pids_before_snapshot); CHKERRQ(ierr);
	}
        // --- STAGE 5: LOOP SYNCHRONIZATION AND CLEANUP ---

        ierr = PetscFree(pids_before_snapshot);
        ierr = PetscFree(migrationList);

        LOG_ALLOW(GLOBAL, LOG_INFO, "End of LocateAllParticlesInGrid pass %d. Total particles migrated globally: %d.\n", passes, global_migrations_this_pass);

    } while (global_migrations_this_pass > 0 && passes < MAX_MIGRATION_PASSES);

    // --- FINAL CHECKS ---
    if (passes >= MAX_MIGRATION_PASSES) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_CONV_FAILED, "Particle migration failed to converge after %d passes. Check for particles oscillating between ranks.", MAX_MIGRATION_PASSES);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Particle Location completed in %d passes.\n", passes);
    LOG_FUNC_TIMER_END_EVENT(EVENT_GlobalParticleLocation, GLOBAL);
    PetscFunctionReturn(0);
}

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
PetscErrorCode ResetAllParticleStatuses(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscInt       n_local;
    PetscInt      *status_p;

    PetscFunctionBeginUser;

    ierr = DMSwarmGetLocalSize(user->swarm, &n_local); CHKERRQ(ierr);

    if (n_local > 0) {
        // Get write access to the status field
        ierr = DMSwarmGetField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void**)&status_p); CHKERRQ(ierr);
        
        for (PetscInt p = 0; p < n_local; ++p) {
            // Only reset particles that are considered settled. This is a small optimization
            // to avoid changing the status of a LOST particle, though resetting all would also be fine.
            if (status_p[p] == ACTIVE_AND_LOCATED) {
                status_p[p] = NEEDS_LOCATION;
            }
        }
        
        // Restore the field
        ierr = DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void**)&status_p); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}
