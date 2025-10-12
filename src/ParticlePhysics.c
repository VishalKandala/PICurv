// ParticlePhysics.c

#include "ParticlePhysics.h"

// Define a buffer size for error messages if not already available
#ifndef ERROR_MSG_BUFFER_SIZE
#define ERROR_MSG_BUFFER_SIZE 256 // Or use PETSC_MAX_PATH_LEN if appropriate
#endif

#undef __FUNCT__
#define __FUNCT__ "UpdateParticleField"
/**
 * @brief Updates a single particle's field based on its state and the specified field name.
 *
 * This function serves as a switchboard for various particle property calculations.
 * Given a particle's data, it computes a new value for the specified field.
 *
 * Currently supported fields:
 *  - "Psi": Calculates the particle's kinetic energy (0.5 * |v|^2).
 *
 * @param[in]     fieldName   The name of the field to update (e.g., "Psi").
 * @param[in]     t           The current simulation time.
 * @param[in]     pos         The particle's physical position.
 * @param[in]     vel         The particle's velocity.
 * @param[in,out] psi_io      A pointer to the particle's current Psi value (used for accumulation if needed).
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode UpdateParticleField(const char *fieldName,
                                                 PetscReal t,
                                                 Cmpnts pos,
                                                 Cmpnts vel,
                                                 PetscReal *psi_io)
{
    PetscFunctionBeginUser;

    if (strcmp(fieldName, "Psi") == 0) {
        // --- Calculate Kinetic Energy ---
        // Formula: Psi = 0.5 * (vx^2 + vy^2 + vz^2)
        *psi_io = 0.5 * (vel.x * vel.x + vel.y * vel.y + vel.z * vel.z);

    }
    // else if (strcmp(fieldName, "SomeOtherField") == 0) {
    //     // ... logic for another field ...
    // }
    else {
        // If the field name is not recognized, we can choose to do nothing or raise an error.
        // For now, we'll just log a warning and do nothing.
        LOG_ALLOW(GLOBAL, LOG_WARNING, "UpdateParticleField: Field name '%s' not recognized. No update performed.\n", fieldName);
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdateFieldForAllParticles"
/**
 * @brief Loops over all local particles and updates a specified field.
 *
 * This function orchestrates the update of a single particle field across the entire
 * local swarm. It gets access to the necessary particle data arrays and calls the
 * `UpdateParticleField` kernel for each particle.
 *
 * @param[in,out] user      Pointer to the UserCtx containing the swarm and simulation context.
 * @param[in]     fieldName The name of the field to update (e.g., "Psi").
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode UpdateFieldForAllParticles(UserCtx *user, const char *fieldName)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    PetscInt       n_local;
    PetscReal      current_time = user->simCtx->ti;

    // Pointers to swarm data arrays
    const PetscReal *pos_arr;
    const PetscReal *vel_arr;
    PetscReal       *psi_arr; // This will be our target field array

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    ierr = DMSwarmGetLocalSize(swarm, &n_local); CHKERRQ(ierr);
    if (n_local == 0) {
        PetscFunctionReturn(0); // Nothing to do on this rank
    }

    // Get read access to position and velocity
    ierr = DMSwarmGetField(swarm, "position", NULL, NULL, (void**)&pos_arr); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "velocity", NULL, NULL, (void**)&vel_arr); CHKERRQ(ierr);
    
    // Get write access to the target field
    ierr = DMSwarmGetField(swarm, fieldName, NULL, NULL, (void**)&psi_arr); CHKERRQ(ierr);

    // Loop over each local particle
    for (PetscInt p = 0; p < n_local; ++p) {
        // Unpack data for the current particle
        Cmpnts current_pos = {pos_arr[3*p + 0], pos_arr[3*p + 1], pos_arr[3*p + 2]};
        Cmpnts current_vel = {vel_arr[3*p + 0], vel_arr[3*p + 1], vel_arr[3*p + 2]};
        
        // Call the kernel to calculate and update the value in the array directly.
        // We pass the address of the specific particle's entry in the psi_arr.
        ierr = UpdateParticleField(fieldName, current_time, current_pos, current_vel, &psi_arr[p]); CHKERRQ(ierr);
    }

    // Restore all field arrays
    ierr = DMSwarmRestoreField(swarm, "position", NULL, NULL, (void**)&pos_arr); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "velocity", NULL, NULL, (void**)&vel_arr); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, fieldName,  NULL, NULL, (void**)&psi_arr); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Updated particle field '%s' for all local particles.\n", fieldName);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdateAllParticlePhysics"
/**
 * @brief Orchestrates the update of all physical properties for particles.
 *
 * This function serves as the top-level entry point for updating particle-specific
 * physical quantities after their position and the surrounding fluid velocity are known.
 * It calls a sequence of more specific update routines for each property.
 *
 * For example, it can be configured to update:
 *  - Particle kinetic energy (stored in "Psi")
 *  - Particle accumulated shear stress
 *  - Particle temperature
 *  - etc.
 *
 * @param[in,out] user  Pointer to the UserCtx containing the swarm and simulation context.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode UpdateAllParticleFields(UserCtx *user)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Updating all particle physical properties...\n");

    // --- Call updater for 'Psi' field (currently calculating Kinetic Energy) ---
    ierr = UpdateFieldForAllParticles(user, "Psi"); CHKERRQ(ierr);

    // --- FUTURE EXPANSION: Add calls for other fields here ---
    /*
    if (simCtx->enable_temperature_solve) {
        ierr = UpdateFieldForAllParticles(user, "Temperature"); CHKERRQ(ierr);
    }
    if (simCtx->enable_species_transport) {
        ierr = UpdateFieldForAllParticles(user, "SpeciesMassFraction"); CHKERRQ(ierr);
    }
    */

    LOG_ALLOW(GLOBAL, LOG_INFO, "All particle physical properties updated.\n");
    
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}