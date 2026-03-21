// ParticleMotion.c

#include "ParticleMotion.h"

// Define a buffer size for error messages if not already available
#ifndef ERROR_MSG_BUFFER_SIZE
#define ERROR_MSG_BUFFER_SIZE 256 // Or use PETSC_MAX_PATH_LEN if appropriate
#endif

#undef __FUNCT__
#define __FUNCT__ "GenerateGaussianNoise"
/**
 * @brief Internal helper implementation: `GenerateGaussianNoise()`.
 * @details Local to this translation unit.
 */
PetscErrorCode GenerateGaussianNoise(PetscRandom rnd, PetscReal *n1, PetscReal *n2)
{
    PetscErrorCode ierr;
    PetscScalar    val1, val2;
    PetscReal      u1, u2;
    PetscReal      magnitude, theta;

    PetscFunctionBeginUser;

    // 1. Get two independent uniform random numbers from the generator
    // PetscRandomGetValue returns a PetscScalar (which might be complex).
    // We take the Real part to ensure this works in both Real and Complex builds.
    ierr = PetscRandomGetValue(rnd, &val1); CHKERRQ(ierr);
    ierr = PetscRandomGetValue(rnd, &val2); CHKERRQ(ierr);

    u1 = PetscRealPart(val1);
    u2 = PetscRealPart(val2);

    // 2. Safety Check: log(0) is undefined (infinity).
    // If the RNG returns exactly 0.0, bump it to a tiny epsilon.
    if (u1 <= 0.0) u1 = 1.0e-14;

    // 3. Box-Muller Transform
    // Formula: R = sqrt(-2 * ln(u1)), Theta = 2 * PI * u2
    magnitude = PetscSqrtReal(-2.0 * PetscLogReal(u1));
    theta     = 2.0 * PETSC_PI * u2; 

    // 4. Calculate independent Normal variables
    *n1 = magnitude * PetscCosReal(theta);
    *n2 = magnitude * PetscSinReal(theta);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CalculateBrownianDisplacement"
/**
 * @brief Internal helper implementation: `CalculateBrownianDisplacement()`.
 * @details Local to this translation unit.
 */
PetscErrorCode CalculateBrownianDisplacement(UserCtx *user, PetscReal diff_eff, Cmpnts *displacement)
{
    PetscErrorCode ierr;
    PetscReal      dt = user->simCtx->dt;
    PetscReal      sigma;
    PetscReal      n_x, n_y, n_z, gaussian_dummy;

    PetscFunctionBeginUser;

    // 1. Initialize output to zero for safety
    displacement->x = 0.0;
    displacement->y = 0.0;
    displacement->z = 0.0;

    // 2. Physical check: Diffusivity cannot be negative. 
    // If 0, there is no Brownian motion.
    if (diff_eff <= 1.0e-12) {
        PetscFunctionReturn(0);
    }

    // 3. Calculate the Scaling Factor (Standard Deviation)
    // Formula: sigma = sqrt(2 * D * dt)
    // Note: dt is inside the root because variance scales linearly with time.
    sigma = PetscSqrtReal(2.0 * diff_eff * dt);

    // 4. Generate 3 Independent Gaussian Random Numbers
    // GenerateGaussianNoise produces 2 numbers at a time. We call it twice.
    
    // Get noise for X and Y
    ierr = GenerateGaussianNoise(user->simCtx->BrownianMotionRNG, &n_x, &n_y); CHKERRQ(ierr);
    
    // Get noise for Z (second sample is intentionally discarded here).
    ierr = GenerateGaussianNoise(user->simCtx->BrownianMotionRNG, &n_z, &gaussian_dummy); CHKERRQ(ierr);

    // 5. Calculate final stochastic displacement
    displacement->x = sigma * n_x;
    displacement->y = sigma * n_y;
    displacement->z = sigma * n_z;

    PetscFunctionReturn(0);
}

#undef __FUNCT__ 
#define __FUNCT__ "UpdateParticlePosition"
/**
 * @brief Internal helper implementation: `UpdateParticlePosition()`.
 * @details Local to this translation unit.
 */
PetscErrorCode UpdateParticlePosition(UserCtx *user, Particle *particle)
{
  PetscFunctionBeginUser; // PETSc macro for error/stack tracing
  PROFILE_FUNCTION_BEGIN;

  PetscErrorCode ierr;
  PetscReal  dt = user->simCtx->dt;
  Cmpnts brownian_disp;

  // 2. Calculate the stochastic kick
  ierr = CalculateBrownianDisplacement(user,particle->diffusivity, &brownian_disp); CHKERRQ(ierr);

  // --- Update Position ---
  // X_new = X_old + ((U_convection + U_diffusivitygradient) * dt) + dX_brownian
  
  particle->loc.x += ((particle->vel.x + particle->diffusivitygradient.x) * dt) + brownian_disp.x;
  particle->loc.y += ((particle->vel.y + particle->diffusivitygradient.y) * dt) + brownian_disp.y;
  particle->loc.z += ((particle->vel.z + particle->diffusivitygradient.z) * dt) + brownian_disp.z;

  PROFILE_FUNCTION_END;
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdateAllParticlePositions"
/**
 * @brief Internal helper implementation: `UpdateAllParticlePositions()`.
 * @details Local to this translation unit.
 */
PetscErrorCode UpdateAllParticlePositions(UserCtx *user)
{
  PetscErrorCode ierr;
  DM swarm = user->swarm;
  PetscInt       nLocal, p;
  PetscReal        *pos = NULL;
  PetscReal        *vel = NULL;
  PetscReal        *diffusivity = NULL;
  Cmpnts           *diffusivitygradient = NULL; 
  PetscReal        *psi = NULL;
  PetscReal        *weights = NULL;
  PetscInt         *cell = NULL;
  PetscInt         *status = NULL;
  PetscInt64       *pid = NULL;
  PetscMPIInt rank;
 
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  PetscFunctionBeginUser;  // PETSc macro for error/stack tracing
  
  PROFILE_FUNCTION_BEGIN;

  // 1) Get the number of local particles
  ierr = DMSwarmGetLocalSize(swarm, &nLocal); CHKERRQ(ierr);
  if (nLocal == 0) {
    LOG_ALLOW(LOCAL,LOG_DEBUG,"[Rank %d] No particles to move/transport. \n",rank);
    PetscFunctionReturn(0);   // nothing to do, no fields held 
  }
  // 2) Access the "position" and "velocity" fields
  ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&pos); CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "velocity", NULL, NULL, (void**)&vel); CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "Diffusivity", NULL, NULL, (void**)&diffusivity); CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "DiffusivityGradient", NULL, NULL, (void**)&diffusivitygradient); CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "Psi", NULL, NULL, (void**)&psi); CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cell); CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "DMSwarm_location_status", NULL, NULL, (void**)&status); CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&pid); CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL,LOG_DEBUG," [Rank %d] No.of Particles to update: %" PetscInt_FMT ".\n",rank,nLocal);

  // 3) Loop over all local particles, updating each position by velocity * dt
  for (p = 0; p < nLocal; p++) {
    // update temporary particle struct
    Particle particle;

    // Unpack: Use the helper to read from swarm arrays into the particle struct
    ierr = UnpackSwarmFields(p, pid, weights, pos, cell, vel, status, diffusivity, diffusivitygradient, psi, &particle); CHKERRQ(ierr);

    // Update position based on velocity and Brownian motion
    ierr = UpdateParticlePosition(user, &particle); CHKERRQ(ierr);
    
    // Update swarm fields
    ierr = UpdateSwarmFields(p, &particle, pos, vel, weights, cell, status, diffusivity, diffusivitygradient, psi); CHKERRQ(ierr);
  }

  // 4) Restore the fields
  ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&pos); CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "velocity", NULL, NULL, (void**)&vel); CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "Diffusivity", NULL, NULL, (void**)&diffusivity); CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "DiffusivityGradient", NULL, NULL, (void**)&diffusivitygradient); CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "Psi", NULL, NULL, (void**)&psi); CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cell); CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "DMSwarm_location_status", NULL, NULL, (void**)&status); CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&pid); CHKERRQ(ierr);


  LOG_ALLOW(LOCAL,LOG_DEBUG,"Particle moved/transported successfully on Rank %d.\n",rank);
  
  PROFILE_FUNCTION_END;

  PetscFunctionReturn(0);
}


/**
 * @brief Internal helper implementation: `IsParticleInBox()`.
 * @details Local to this translation unit.
 */
static inline PetscBool IsParticleInBox(const BoundingBox *bbox, const Cmpnts *pos) {
    return (pos->x >= bbox->min_coords.x && pos->x <= bbox->max_coords.x &&
            pos->y >= bbox->min_coords.y && pos->y <= bbox->max_coords.y &&
            pos->z >= bbox->min_coords.z && pos->z <= bbox->max_coords.z);
}


#undef __FUNCT__
#define __FUNCT__ "CheckAndRemoveOutOfBoundsParticles"

/**
 * @brief Internal helper implementation: `CheckAndRemoveOutOfBoundsParticles()`.
 * @details Local to this translation unit.
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

#undef __FUNCT__
#define __FUNCT__ "CheckAndRemoveLostParticles"
/**
 * @brief Internal helper implementation: `CheckAndRemoveLostParticles()`.
 * @details Local to this translation unit.
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
    PROFILE_FUNCTION_BEGIN;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Checking for and removing LOST particles...\n", rank);

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
    
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "SetMigrationRanks"
/**
 * @brief Internal helper implementation: `SetMigrationRanks()`.
 * @details Local to this translation unit.
 */
PetscErrorCode SetMigrationRanks(UserCtx* user, const MigrationInfo *migrationList, PetscInt migrationCount)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    PetscInt       p_idx;
    PetscInt      *rankField = NULL; // Field storing target rank

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // Ensure the migration rank field exists
    ierr = DMSwarmGetField(swarm, "DMSwarm_rank", NULL, NULL, (void **)&rankField); CHKERRQ(ierr);

    // Set the target rank for migrating particles
    for(p_idx = 0; p_idx < migrationCount; ++p_idx) {
        rankField[migrationList[p_idx].local_index] = migrationList[p_idx].target_rank;
    }

    ierr = DMSwarmRestoreField(swarm, "DMSwarm_rank", NULL, NULL, (void **)&rankField); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PerformMigration"

/**
 * @brief Implementation of \ref PerformMigration().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/ParticleMotion.h`.
 * @see PerformMigration()
 */
PetscErrorCode PerformMigration(UserCtx *user)
{
    PetscErrorCode ierr;
    DM swarm = user->swarm;
    PetscMPIInt rank;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Rank %d: Starting DMSwarmMigrate...\n", rank);

    // Perform the migration - PETSC_TRUE removes particles that fail to land
    // in a valid cell on the target rank (or were marked with an invalid rank).
    ierr = DMSwarmMigrate(swarm, PETSC_TRUE); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Rank %d: Migration complete.\n", rank);
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

//-----------------------------------------------------------------------------
// MODULE (COUNT): Calculates Particle Count Per Cell - REVISED FOR DMDAVecGetArray
//-----------------------------------------------------------------------------

#undef __FUNCT__
#define __FUNCT__ "CalculateParticleCountPerCell"
/**
 * @brief Implementation of \ref CalculateParticleCountPerCell().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/logging.h`.
 * @see CalculateParticleCountPerCell()
 */
PetscErrorCode CalculateParticleCountPerCell(UserCtx *user) {
    PetscErrorCode ierr;
    DM             da = user->da;
    DM             swarm = user->swarm;
    Vec            countVec = user->ParticleCount;
    Vec            localcountVec = user->lParticleCount;
    PetscInt       nlocal, p;
    PetscInt       *global_cell_id_arr; // Read GLOBAL cell IDs
    PetscScalar    ***count_arr_3d;     // Use 3D accessor
    PetscInt64       *PID_arr;
    PetscMPIInt    rank;
    char           msg[ERROR_MSG_BUFFER_SIZE];
    PetscInt       particles_counted_locally = 0;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // --- Input Validation ---
    if (!da) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx->da is NULL.");
    if (!swarm) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx->swarm is NULL.");
    if (!countVec) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx->ParticleCount is NULL.");
    // Check DOF of da
    PetscInt count_dof;
    ierr = DMDAGetInfo(da, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &count_dof, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
    if (count_dof != 1) {
        PetscSNPrintf(msg, sizeof(msg), "countDM must have DOF=1, got %" PetscInt_FMT ".", count_dof);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "%s", msg);
    }

    // --- Zero the local count vector ---
    ierr = VecSet(localcountVec, 0.0); CHKERRQ(ierr);

    // --- Get Particle Data ---
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Accessing particle data.\n");
    ierr = DMSwarmGetLocalSize(swarm, &nlocal); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm,"DMSwarm_CellID", NULL, NULL, (void **)&global_cell_id_arr); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm,"DMSwarm_pid",NULL,NULL,(void **)&PID_arr);CHKERRQ(ierr);

    // --- Get Grid Vector Array using DMDA accessor ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Accessing ParticleCount vector array (using DMDAVecGetArray).\n");
    ierr = DMDAVecGetArray(da, localcountVec, &count_arr_3d); CHKERRQ(ierr);

    // Get local owned + ghosted range for writing into ghost slots.
    PetscInt gxs, gys, gzs, gxm, gym, gzm;
    ierr = DMDAGetGhostCorners(da, &gxs, &gys, &gzs, &gxm, &gym, &gzm); CHKERRQ(ierr);

    // --- Accumulate Counts Locally ---
    LOG_ALLOW(LOCAL, LOG_DEBUG, "CalculateParticleCountPerCell (Rank %d): Processing %" PetscInt_FMT " local particles using GLOBAL CellIDs.\n",rank,nlocal);
    for (p = 0; p < nlocal; p++) {
        // Read the GLOBAL indices stored for this particle
        PetscInt i_geom = global_cell_id_arr[p * 3 + 0]; // Global i index
        PetscInt j_geom = global_cell_id_arr[p * 3 + 1]; // Global j index
        PetscInt k_geom = global_cell_id_arr[p * 3 + 2]; // Global k index

        // Apply the shift to ensure ParticleCount follows the indexing convention for cell-centered data in this codebase.
        PetscInt i = (PetscInt)i_geom + 1; // Shift for cell-centered
        PetscInt j = (PetscInt)j_geom + 1; // Shift for cell-centered
        PetscInt k = (PetscInt)k_geom + 1; // Shift for cell-centered

        // *** Bounds check is implicitly handled by DMDAVecGetArray for owned+ghost region ***
        // However, accessing outside this region using global indices WILL cause an error.
        // A preliminary check might still be wise if global IDs could be wild.
        // We rely on LocateAllParticles to provide valid global indices [0..IM-1] etc.
 
        LOG_LOOP_ALLOW(LOCAL, LOG_VERBOSE, p, 100,
                       "[Rank %d] Read CellID for p=%" PetscInt_FMT ", PID = %" PetscInt64_FMT ": (%" PetscInt_FMT ", %" PetscInt_FMT ", %" PetscInt_FMT ")\n",
                       rank, p, PID_arr[p], i, j, k);
    
        // Check if the global index (i,j,k) falls within the local + ghost range
        if (i >= gxs && i < gxs + gxm  && 
            j >= gys && j < gys + gym  && // Adjust based on actual ghost width
            k >= gzs && k < gzs + gzm  )   // This check prevents definite crashes but doesn't guarantee ownership
        {

             // Increment count at the location corresponding to GLOBAL index (I,J,K)
	  //  LOG_ALLOW(LOCAL, LOG_DEBUG, "CalculateParticleCountPerCell (Rank %d): Particle %d with global CellID (%d, %d, %d) incremented with a particle.\n",rank, p, i, j, k);
             count_arr_3d[k][j][i] += 1.0;
             particles_counted_locally++;
         } else {
              // This particle's global ID is likely outside the range this rank handles (even ghosts)
              // note: this is not necessarily an error if the particle is legitimately outside the local+ghost region
              LOG_ALLOW(LOCAL, LOG_VERBOSE,
                        "(Rank %d): Skipping particle %" PetscInt64_FMT " with global CellID (%" PetscInt_FMT ", %" PetscInt_FMT ", %" PetscInt_FMT ") - likely outside local+ghost range.\n",
                        rank, PID_arr[p], i, j, k);
	}
    }
    LOG_ALLOW(LOCAL, LOG_DEBUG, "(Rank %d): Local counting finished. Processed %" PetscInt_FMT " particles locally.\n", rank, particles_counted_locally);

    // --- Restore Access ---
    ierr = DMDAVecRestoreArray(da, localcountVec, &count_arr_3d); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm,"DMSwarm_CellID", NULL, NULL, (void **)&global_cell_id_arr); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm,"DMSwarm_pid",NULL,NULL,(void **)&PID_arr);CHKERRQ(ierr);

    // --- Assemble Global Vector ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Assembling global ParticleCount vector.\n");
    ierr = VecZeroEntries(countVec); CHKERRQ(ierr); // Ensure global vector is zeroed before accumulation
    ierr = DMLocalToGlobalBegin(da, localcountVec, ADD_VALUES, countVec); CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(da, localcountVec, ADD_VALUES, countVec); CHKERRQ(ierr);
    /* 
     * OPTIONAL: Synchronize Ghosts for Stencil Operations
     * If a future function needs to read ParticleCount from neighbor cells (e.g., density smoothing 
     * or gradient calculations), uncomment the following lines to update the ghost slots 
     * in user->lParticleCount with the final summed values.
     * 
     ierr = UpdateLocalGhosts(user,"ParticleCount"); CHKERRQ(ierr);
     */

    // --- Verification Logging ---
    PetscReal total_counted_particles = 0.0, max_count_in_cell = 0.0;
    ierr = VecSum(countVec, &total_counted_particles); CHKERRQ(ierr);
    PetscInt max_idx_global = -1;
    ierr = VecMax(countVec, &max_idx_global, &max_count_in_cell); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Total counted globally = %.0f, Max count in cell = %.0f\n",
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
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Particle counting complete.\n");


    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "ResizeSwarmGlobally"
/**
 * @brief Implementation of \ref ResizeSwarmGlobally().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/ParticleMotion.h`.
 * @see ResizeSwarmGlobally()
 */

PetscErrorCode ResizeSwarmGlobally(DM swarm, PetscInt N_target)
{
    PetscErrorCode ierr;
    PetscInt       N_current, nlocal_current;
    PetscMPIInt    rank;
    MPI_Comm       comm;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    ierr = PetscObjectGetComm((PetscObject)swarm, &comm); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);
    ierr = DMSwarmGetSize(swarm, &N_current); CHKERRQ(ierr);
    ierr = DMSwarmGetLocalSize(swarm, &nlocal_current); CHKERRQ(ierr);

    PetscInt delta = N_target - N_current;

    if (delta == 0) {
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0); // Nothing to do
    }

    if (delta < 0) { // Remove particles
        PetscInt num_to_remove_global = -delta;
        LOG_ALLOW(GLOBAL, LOG_INFO, "Current size %d > target size %d. Removing %d particles globally.\n", N_current, N_target, num_to_remove_global);

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
        LOG_ALLOW(GLOBAL, LOG_INFO, "Current size %d < target size %d. Adding %d particles globally.\n", N_current, N_target, num_to_add_global);
        ierr = DMSwarmAddNPoints(swarm, num_to_add_global); CHKERRQ(ierr);
        // Note: Added particles will have uninitialized field data. Reading will overwrite.
    }

    // Verify final size
    PetscInt N_final;
    ierr = DMSwarmGetSize(swarm, &N_final); CHKERRQ(ierr);
    if (N_final != N_target) {
        SETERRQ(comm, PETSC_ERR_PLIB, "Failed to resize swarm: expected %d particles, got %d", N_target, N_final);
    }
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Swarm successfully resized to %d particles.\n", N_final);
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PreCheckAndResizeSwarm"
/**
 * @brief Internal helper implementation: `PreCheckAndResizeSwarm()`.
 * @details Local to this translation unit.
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
    SimCtx        *simCtx = user->simCtx;
    char          *source_path;
    
    // NOTE: Your filename format has a hardcoded "_0" which is typical for
    // PETSc when writing a parallel object from a single rank.
    // If you ever write in parallel, PETSc might create one file per rank.
    // The current logic assumes a single file written by one process.
    const int      placeholder_int = 0;

    // Setup the I/O environment


    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    ierr = PetscObjectGetComm((PetscObject)user->swarm, &comm); CHKERRQ(ierr);
    ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);

    // First, determine the top-level source directory based on the execution mode.
    if (simCtx->exec_mode == EXEC_MODE_SOLVER) {
        source_path = simCtx->restart_dir;
    } else if (simCtx->exec_mode == EXEC_MODE_POSTPROCESSOR) {
        source_path = simCtx->pps->source_dir;
    } else {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "Invalid execution mode for reading simulation fields.");
    }

    // Set the current I/O directory context
    ierr = PetscSNPrintf(simCtx->_io_context_buffer, sizeof(simCtx->_io_context_buffer),
                         "%s/%s", source_path, simCtx->particle_subdir); CHKERRQ(ierr);

    simCtx->current_io_directory = simCtx->_io_context_buffer;

        // --- Construct filename using the specified format ---
    ierr = PetscSNPrintf(filename, sizeof(filename), "%s/%s%05" PetscInt_FMT "_%d.%s",
                         simCtx->current_io_directory,refFieldName, ti, placeholder_int, ext); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Checking particle count for timestep %d using ref file '%s'.\n", ti, filename);

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
    user->simCtx->np = N_file;

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ReinitializeParticlesOnInletSurface"
/**
 * @brief Internal helper implementation: `ReinitializeParticlesOnInletSurface()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ReinitializeParticlesOnInletSurface(UserCtx *user, PetscReal currentTime, PetscInt step)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;                        // MPI rank of the current process
    DM             swarm = user->swarm;         // The particle swarm DM
    PetscReal      *positions_field = NULL;     // Pointer to swarm field for physical positions
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

    PROFILE_FUNCTION_BEGIN;

    // This function is only relevant for surface initialization mode and if an inlet face is defined.
    if ((user->simCtx->ParticleInitialization != 0 && user->simCtx->ParticleInitialization !=3) || !user->inletFaceDefined) {
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
    ierr = DMDAGetCorners(user->da, &xs_gnode_rank, &ys_gnode_rank, &zs_gnode_rank, NULL, NULL, NULL); CHKERRQ(ierr);

    // Modification to IM_nodes_global etc. to account for 1-cell halo in each direction.
    IM_nodes_global -= 1; JM_nodes_global -= 1; KM_nodes_global -= 1; 

    const PetscInt IM_cells_global = IM_nodes_global > 0 ? IM_nodes_global - 1 : 0;
    const PetscInt JM_cells_global = JM_nodes_global > 0 ? JM_nodes_global - 1 : 0;
    const PetscInt KM_cells_global = KM_nodes_global > 0 ? KM_nodes_global - 1 : 0;



    // Check if this rank is responsible for (part of) the designated inlet surface
    ierr = CanRankServiceInletFace(user, &info, IM_nodes_global, JM_nodes_global, KM_nodes_global, &can_this_rank_service_inlet); CHKERRQ(ierr);

    // Get coordinate array and swarm fields for modification
    ierr = DMGetCoordinatesLocal(user->da, &Coor_local); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, Coor_local, (void*)&coor_nodes_local_array); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "position",       NULL, NULL, (void**)&positions_field); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_pid",    NULL, NULL, (void**)&particleIDs);    CHKERRQ(ierr); // For logging
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cell_ID_field);  CHKERRQ(ierr);

    if (!can_this_rank_service_inlet) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f, Step=%d] Rank %d cannot service inlet face %s. Skipping re-initialization of %d particles.\n", currentTime, step, rank, BCFaceToString(user->identifiedInletBCFace), nlocal_current);
        
        //  FALLBACK ACTION: Reset position fields to Inlet center for migration and cell ID to -1 for safety.
        LOG_ALLOW(LOCAL, LOG_DEBUG, "[T=%.4f, Step=%d] Rank %d is resetting %d local particles to inlet center (%.6f, %.6f, %.6f) for migration.\n", currentTime, step, rank, nlocal_current, user->simCtx->CMx_c, user->simCtx->CMy_c, user->simCtx->CMz_c);

        for(PetscInt p = 0; p < nlocal_current; p++){
            positions_field[3*p+0] = user->simCtx->CMx_c;
            positions_field[3*p+1] = user->simCtx->CMy_c;
            positions_field[3*p+2] = user->simCtx->CMz_c;

            cell_ID_field[3*p+0] = -1;
            cell_ID_field[3*p+1] = -1;
            cell_ID_field[3*p+2] = -1;
        }
        
        // Cleanup: restore swarm fields/coordinate array
        ierr = DMDAVecRestoreArrayRead(user->fda, Coor_local, (void*)&coor_nodes_local_array); CHKERRQ(ierr);
        ierr = DMSwarmRestoreField(swarm, "position",       NULL, NULL, (void**)&positions_field); CHKERRQ(ierr);
        ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid",    NULL, NULL, (void**)&particleIDs);    CHKERRQ(ierr); // For logging
        ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cell_ID_field);  CHKERRQ(ierr);
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Rank %d is on inlet face %s. Attempting to re-place %d local particles.\n", currentTime, step, rank, BCFaceToString(user->identifiedInletBCFace), nlocal_current);

    // Initialize fresh RNGs for this re-placement to ensure good distribution
    ierr = InitializeLogicalSpaceRNGs(&rand_logic_reinit_i, &rand_logic_reinit_j, &rand_logic_reinit_k); CHKERRQ(ierr);
    // Optional: Seed RNGs for deterministic behavior if required, e.g., based on rank and step.
    // PetscRandomSetSeed(rand_logic_i, (unsigned long)rank*1000 + step + 100); PetscRandomSeed(rand_logic_i); // Example

    // Loop over all particles currently local to this rank
    for (PetscInt p = 0; p < nlocal_current; p++) {
        PetscInt  ci_metric_lnode, cj_metric_lnode, ck_metric_lnode; // Local node indices (of rank's DA patch) for cell origin
        PetscReal xi_metric_logic, eta_metric_logic, zta_metric_logic; // Intra-cell logical coordinates
        Cmpnts    phys_coords = {0.0,0.0,0.0}; // To store newly calculated physical coordinates
        PetscBool particle_was_placed = PETSC_FALSE;

        if(user->simCtx->ParticleInitialization == PARTICLE_INIT_SURFACE_RANDOM){
        // Get random cell on this rank's portion of the inlet and random logical coords within it
            ierr = GetRandomCellAndLogicalCoordsOnInletFace(user, &info, xs_gnode_rank, ys_gnode_rank, zs_gnode_rank,
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
            particle_was_placed = PETSC_TRUE;                            

        }else if(user->simCtx->ParticleInitialization == PARTICLE_INIT_SURFACE_EDGES){
            PetscBool placement_flag = PETSC_FALSE;
            ierr = GetDeterministicFaceGridLocation(user, &info, xs_gnode_rank, ys_gnode_rank, zs_gnode_rank,
                                                IM_cells_global, JM_cells_global, KM_cells_global,
                                                particleIDs[p],
                                                &ci_metric_lnode, &cj_metric_lnode, &ck_metric_lnode,
                                                &xi_metric_logic, &eta_metric_logic, &zta_metric_logic,&placement_flag); CHKERRQ(ierr);
        

            if(placement_flag){
                // Convert these logical coordinates to physical coordinates
                ierr = MetricLogicalToPhysical(user, coor_nodes_local_array,
                                        ci_metric_lnode, cj_metric_lnode, ck_metric_lnode,
                                        xi_metric_logic, eta_metric_logic, zta_metric_logic,
                                        &phys_coords); CHKERRQ(ierr);

                // Update the particle's position in the swarm fields
                positions_field[3*p+0] = phys_coords.x;
                positions_field[3*p+1] = phys_coords.y;
                positions_field[3*p+2] = phys_coords.z;
                particle_was_placed =  PETSC_TRUE;
            } else{
                // Deterministic placement failed (particle migrated to rank where formula says it doesn't belong)
                // Fall back to random placement on this rank's portion of inlet surface
                LOG_ALLOW(GLOBAL, LOG_WARNING, "Rank %d: Particle PID %ld deterministic placement failed (belongs to different rank). Falling back to random placement.\n", rank, particleIDs[p]);

                ierr = GetRandomCellAndLogicalCoordsOnInletFace(user, &info, xs_gnode_rank, ys_gnode_rank, zs_gnode_rank,
                                                        IM_nodes_global, JM_nodes_global, KM_nodes_global,
                                                        &rand_logic_reinit_i, &rand_logic_reinit_j, &rand_logic_reinit_k,
                                                        &ci_metric_lnode, &cj_metric_lnode, &ck_metric_lnode,
                                                        &xi_metric_logic, &eta_metric_logic, &zta_metric_logic); CHKERRQ(ierr);

                // Convert to physical coordinates
                ierr = MetricLogicalToPhysical(user, coor_nodes_local_array,
                                            ci_metric_lnode, cj_metric_lnode, ck_metric_lnode,
                                            xi_metric_logic, eta_metric_logic, zta_metric_logic,
                                            &phys_coords); CHKERRQ(ierr);

                // Update particle position
                positions_field[3*p+0] = phys_coords.x;
                positions_field[3*p+1] = phys_coords.y;
                positions_field[3*p+2] = phys_coords.z;
                particle_was_placed = PETSC_TRUE;
            }                                    

        } else{
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "ReinitializeParticlesOnInletSurface only supports ParticleInitialization modes 0 and 3.");
        }

        if(particle_was_placed){
        particles_actually_reinitialized_count++;

        cell_ID_field[3*p+0] = -1;
        cell_ID_field[3*p+1] = -1;
        cell_ID_field[3*p+2] = -1;

        LOG_LOOP_ALLOW(LOCAL, LOG_VERBOSE, p, (nlocal_current > 20 ? nlocal_current/10 : 1), // Sampled logging
            "Rank %d: PID %ld (idx %ld) RE-PLACED. CellOriginNode(locDAIdx):(%d,%d,%d). LogicCoords: (%.2e,%.2f,%.2f). PhysCoords: (%.6f,%.6f,%.6f).\n",
            rank, particleIDs[p], (long)p, 
            ci_metric_lnode, cj_metric_lnode, ck_metric_lnode,
            xi_metric_logic, eta_metric_logic, zta_metric_logic, 
            phys_coords.x, phys_coords.y, phys_coords.z);
        }
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
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid",    NULL, NULL, (void**)&particleIDs);    CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cell_ID_field);  CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, Coor_local, (void*)&coor_nodes_local_array); CHKERRQ(ierr);


    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GetLocalPIDSnapshot"
/**
 * @brief Internal helper implementation: `GetLocalPIDSnapshot()`.
 * @details Local to this translation unit.
 */
PetscErrorCode GetLocalPIDSnapshot(const PetscInt64 pid_field[], 
                                   PetscInt n_local, 
                                   PetscInt64 **pids_snapshot_out)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

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

        PROFILE_FUNCTION_END;
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


    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "AddToMigrationList"
/**
 * @brief Internal helper implementation: `AddToMigrationList()`.
 * @details Local to this translation unit.
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

    PROFILE_FUNCTION_BEGIN;

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
        LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d]: Reallocated migrationList capacity from %d to %d.\n",
                  rank, old_capacity, new_capacity);
    }

    // --- 3. Add the new migration data to the list ---
    // Dereference the pointer-to-a-pointer to get the actual array.
    MigrationInfo *list = *migration_list_p;

    list[*count_p].local_index = particle_local_idx;
    list[*count_p].target_rank = destination_rank;

    // --- 4. Increment the count of items in the list ---
    (*count_p)++;
    

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "FlagNewComersForLocation"
/**
 * @brief Internal helper implementation: `FlagNewcomersForLocation()`.
 * @details Local to this translation unit.
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

    PROFILE_FUNCTION_BEGIN;

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
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d]: Checking for newcomers. Size before: %d, Size after: %d\n",
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
            
            LOG_ALLOW(LOCAL, LOG_VERBOSE, "[Rank %d]: Flagged newcomer PID %ld at local index %d as NEEDS_LOCATION.\n",
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


    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MigrateRestartParticlesUsingCellID"
/**
 * @brief Internal helper implementation: `MigrateRestartParticlesUsingCellID()`.
 * @details Local to this translation unit.
 */
PetscErrorCode MigrateRestartParticlesUsingCellID(UserCtx *user)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    PetscInt       nlocal;
    PetscInt       *cell_p = NULL;
    PetscInt64     *pid_p = NULL;
    PetscMPIInt    rank;
    
    MigrationInfo  *migrationList = NULL;
    PetscInt       local_migration_count = 0;
    PetscInt       migrationListCapacity = 0;
    PetscInt       global_migration_count = 0;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    ierr = DMSwarmGetLocalSize(swarm, &nlocal); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Checking %d restart particles for direct migration using CellIDs.\n", nlocal);

    if (nlocal > 0) {
        ierr = DMSwarmGetField(swarm, "DMSwarm_CellID",          NULL, NULL, (void**)&cell_p);   CHKERRQ(ierr);
        ierr = DMSwarmGetField(swarm, "DMSwarm_pid",             NULL, NULL, (void**)&pid_p);    CHKERRQ(ierr);
        
        // Note: We do NOT need to modify the status field here. 
        // We trust the loaded status (ACTIVE_AND_LOCATED) is correct for the destination rank.

        for (PetscInt p_idx = 0; p_idx < nlocal; ++p_idx) {
            PetscInt ci = cell_p[3*p_idx + 0];
            PetscInt cj = cell_p[3*p_idx + 1];
            PetscInt ck = cell_p[3*p_idx + 2];

            /* Skip particles with invalid Cell IDs (will be handled by LocateAllParticles) */
            if (ci < 0 || cj < 0 || ck < 0) {
                continue;
            }

            PetscMPIInt owner_rank;
            ierr = FindOwnerOfCell(user, ci, cj, ck, &owner_rank); CHKERRQ(ierr);

            if (owner_rank != -1 && owner_rank != rank) {
                /* Particle belongs to another rank - migrate it */
                ierr = AddToMigrationList(&migrationList, &migrationListCapacity, &local_migration_count,
                                          p_idx, owner_rank); CHKERRQ(ierr);
                
                LOG_ALLOW(LOCAL, LOG_VERBOSE, "[PID %ld] Direct migration: Cell (%d,%d,%d) belongs to Rank %d (Current: %d).\n",
                          (long)pid_p[p_idx], ci, cj, ck, owner_rank, rank);
            }
        }

        ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID",          NULL, NULL, (void**)&cell_p);   CHKERRQ(ierr);
        ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid",             NULL, NULL, (void**)&pid_p);    CHKERRQ(ierr);
    }

    /* Check if any rank needs to migrate particles */
    ierr = MPI_Allreduce(&local_migration_count, &global_migration_count, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);

    if (global_migration_count > 0) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Fast restart migration: Directly migrating %d particles using CellIDs.\n", global_migration_count);
        ierr = SetMigrationRanks(user, migrationList, local_migration_count); CHKERRQ(ierr);
        ierr = PerformMigration(user); CHKERRQ(ierr);
        /* We do NOT flag newcomers here. We trust their loaded status (ACTIVE_AND_LOCATED) */
        /* is valid for their destination rank. */
    } else {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Fast restart migration: All particles are already on correct ranks.\n");
    }

    ierr = PetscFree(migrationList); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GuessParticleOwnerWithBBox"
/**
 * @brief Internal helper implementation: `GuessParticleOwnerWithBBox()`.
 * @details Local to this translation unit.
 * @note Testing status:
 *       The current direct surface reaches this helper through orchestrator
 *       tests, but direction-complete immediate-neighbor coverage and the
 *       explicit "not found in any rank" path are still targeted for future
 *       bespoke tests.
 */
static PetscErrorCode GuessParticleOwnerWithBBox(UserCtx *user,
                                                 const Particle *particle,
                                                 const BoundingBox *bboxlist,
                                                 PetscMPIInt *guess_rank_out)
{
    PetscErrorCode  ierr;
    PetscMPIInt     rank, size;
    const RankNeighbors *neighbors = &user->neighbors; // Use a direct pointer for clarity
    const BoundingBox   *localBBox = &user->bbox;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

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

    // --- Step 0: Check if the particle is inside the CURRENT rank's bounding box FIRST. ---
    // This handles the common case of initial placement where a particle is "lost" but physically local.
    if (IsParticleInBox(localBBox, &particle->loc)) {
      *guess_rank_out = rank;
      LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %ld]: Fast path guess SUCCESS. Particle is within the local (Rank %d) bounding box.\n",
		particle->PID, rank);

      PROFILE_FUNCTION_END;
      PetscFunctionReturn(0); // Found it, we're done.
    }
    // --- 2. Fast Path: Check Immediate Neighbors Based on Exit Direction ---
    
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

        PROFILE_FUNCTION_END;          
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
	  LOG_ALLOW(LOCAL,LOG_VERBOSE, "[Rank %d] Checking PID %lld at (%.4f, %.4f, %.4f) against Rank %d's box: [(%.4f, %.4f, %.4f) to (%.4f, %.4f, %.4f)]. Result: %s\n",
				  (int)rank, (long long)particle->PID,
				  particle->loc.x, particle->loc.y, particle->loc.z,
				  (int)r,
				  bboxlist[r].min_coords.x, bboxlist[r].min_coords.y, bboxlist[r].min_coords.z,
				  bboxlist[r].max_coords.x, bboxlist[r].max_coords.y, bboxlist[r].max_coords.z,
				  is_in ? "INSIDE" : "OUTSIDE");
	  
            *guess_rank_out = r;
            LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %ld]: Fallback search SUCCESS. Found in Rank %d.\n",
                      particle->PID, *guess_rank_out);
            
            PROFILE_FUNCTION_END;
            PetscFunctionReturn(0); // Found it, we're done.
        }
    }

    // If the code reaches here, the particle was not found in any rank's bounding box.
    LOG_ALLOW(LOCAL, LOG_WARNING, "[PID %ld]: Guess FAILED. Particle not found in any rank's bounding box.\n",
              particle->PID);
    
    // The guess_rank_out will remain -1, signaling failure to the caller.
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LocateAllParticlesInGrid"
/**
 * @brief Implementation of \ref LocateAllParticlesInGrid().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/ParticleMotion.h`.
 * @see LocateAllParticlesInGrid()
 */
PetscErrorCode LocateAllParticlesInGrid(UserCtx *user,BoundingBox *bboxlist)
{
    PetscErrorCode ierr;
    PetscInt       passes = 0;
    const PetscInt MAX_MIGRATION_PASSES = 50; // Safety break for runaway loops
    PetscInt       global_migrations_this_pass;
    PetscMPIInt    rank;
    PetscInt total_migrated_this_timestep = 0;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = ResetSearchMetrics(user->simCtx); CHKERRQ(ierr);
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

	      LOG_ALLOW(LOCAL,LOG_VERBOSE,
			"Local Particle idx=%d, PID=%ld, status=%s, cell=(%d, %d, %d)\n",
			p_idx,
			(long)pid_p[p_idx],
			ParticleLocationStatusToString((ParticleLocationStatus)status_p[p_idx]),
			cell_p[3*p_idx],
			cell_p[3*p_idx+1],
			cell_p[3*p_idx+2]);

                if (status_p[p_idx] == ACTIVE_AND_LOCATED) {
		  LOG_ALLOW(LOCAL,LOG_VERBOSE," [rank %d][PID %ld] skipped in pass %d as it is already located at (%d,%d,%d).\n",rank,pid_p[p_idx],passes,cell_p[3*p_idx],cell_p[3*p_idx + 1],cell_p[3*p_idx + 2]);
                    continue;
                }

                // UNPACK: Create a temporary C struct for easier processing using our helper.
                Particle current_particle;

		    //	LOG_ALLOW(LOCAL,LOG_DEBUG,"about to unpack p_idx=%d (PID=%ld)\n",p_idx, (long)pid_p[p_idx]);
		
                ierr = UnpackSwarmFields(p_idx, pid_p, weights_p, pos_p, cell_p, vel_p, status_p,NULL,NULL,NULL,&current_particle); CHKERRQ(ierr);

		    //		LOG_ALLOW(LOCAL,LOG_DEBUG,"unpacked p_idx=%d → cell[0]=%d, status=%s\n",p_idx, current_particle.cell[0], ParticleLocationStatusToString((ParticleLocationStatus)current_particle.location_status));

		    ParticleLocationStatus final_status = (ParticleLocationStatus)status_p[p_idx];


		    // CASE 1: Particle has a valid prior cell index.
                // It has moved, so we only need to run the robust walk from its last known location.
                if (current_particle.cell[0] >= 0) {
		            LOG_ALLOW(LOCAL, LOG_VERBOSE, "[PID %ld] has valid prior cell. Strategy: Robust Walk from previous cell.\n", current_particle.PID);
		            ierr = LocateParticleOrFindMigrationTarget(user, &current_particle, &final_status); CHKERRQ(ierr);
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
		     ierr = LocateParticleOrFindMigrationTarget(user, &current_particle, &final_status); CHKERRQ(ierr);
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
                     ierr = UpdateSwarmFields(p_idx, &current_particle, pos_p, vel_p, weights_p, cell_p, status_p,NULL,NULL,NULL); CHKERRQ(ierr);
                }
		*/
                // CASE 2: Particle is "lost" (cell = -1). Strategy: Guess -> Verify.
                else {
		  LOG_ALLOW(LOCAL, LOG_VERBOSE, "[PID %ld] has invalid cell. Strategy: Guess Owner -> Find Cell.\n",current_particle.PID);
                    
		  PetscMPIInt guessed_owner_rank = MPI_PROC_NULL;
		  ierr = GuessParticleOwnerWithBBox(user, &current_particle, bboxlist, &guessed_owner_rank); CHKERRQ(ierr);

		  // If the guess finds a DIFFERENT rank, we can mark for migration and skip the walk.
		  if (guessed_owner_rank != MPI_PROC_NULL && guessed_owner_rank != rank) {
		    user->simCtx->searchMetrics.bboxGuessSuccessCount++;
		    LOG_ALLOW(LOCAL, LOG_VERBOSE, "[PID %ld] Guess SUCCESS: Found migration target Rank %d. Finalizing.\n", current_particle.PID, guessed_owner_rank);
		    final_status = MIGRATING_OUT;
		    current_particle.destination_rank = guessed_owner_rank;
		  } 
		  else {
		    user->simCtx->searchMetrics.bboxGuessFallbackCount++;

		    // This block runs if the guess either failed (rank is NULL) or found the particle is local (rank is self).
		    // In BOTH cases, the situation is unresolved, and we MUST fall back to the robust walk.
		    if (guessed_owner_rank == rank) {
		      LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %ld] Guess determined particle is local. Proceeding to robust walk to find cell.\n", current_particle.PID);
		    } else { // guessed_owner_rank == MPI_PROC_NULL
		      LOG_ALLOW(LOCAL, LOG_WARNING, "[PID %ld] Guess FAILED to find an owner. Proceeding to robust walk for definitive search.\n", current_particle.PID);
		    }
                        
		    ierr = LocateParticleOrFindMigrationTarget(user, &current_particle, &final_status); CHKERRQ(ierr);
		  }
                }
		
                // --- PROCESS THE FINAL, DEFINITIVE STATUS ---
                current_particle.location_status = final_status;
                ierr = UpdateSwarmFields(p_idx, &current_particle, pos_p, vel_p, weights_p, cell_p, status_p,NULL,NULL,NULL); CHKERRQ(ierr);
                
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

	    total_migrated_this_timestep += global_migrations_this_pass;

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

        LOG_ALLOW(GLOBAL, LOG_INFO, "End of pass %d. Total particles migrated globally: %d.\n", passes, global_migrations_this_pass);

    } while (global_migrations_this_pass > 0 && passes < MAX_MIGRATION_PASSES);

    // --- FINAL CHECKS ---
    if (passes >= MAX_MIGRATION_PASSES) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_CONV_FAILED, "Particle migration failed to converge after %d passes. Check for particles oscillating between ranks.", MAX_MIGRATION_PASSES);
    }

    user->simCtx->particlesMigratedLastStep = total_migrated_this_timestep;
    user->simCtx->migrationPassesLastStep = passes;
    user->simCtx->searchMetrics.maxParticlePassDepth = PetscMax(user->simCtx->searchMetrics.maxParticlePassDepth, passes);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Particle Location completed in %d passes.\n", passes);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ResetAllParticleStatuses"
/**
 * @brief Implementation of \ref ResetAllParticleStatuses().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/ParticleMotion.h`.
 * @see ResetAllParticleStatuses()
 */
PetscErrorCode ResetAllParticleStatuses(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscInt       n_local;
    PetscInt      *status_p;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

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


    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
