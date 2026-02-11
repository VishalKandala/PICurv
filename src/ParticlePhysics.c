#include "ParticlePhysics.h"

#ifndef ERROR_MSG_BUFFER_SIZE
#define ERROR_MSG_BUFFER_SIZE 256
#endif

#undef __FUNCT__
#define __FUNCT__ "UpdateParticleField"
/**
 * @brief Updates a single particle's field based on its state and physics model.
 *
 * Implements the IEM (Interaction by Exchange with the Mean) model for scalar mixing.
 * 
 * Physics: dPsi/dt = -Omega * (Psi - <Psi>)
 * Solution: Psi_new = <Psi> + (Psi_old - <Psi>) * exp(-Omega * dt)
 *
 * @param[in]     fieldName   Name of the field (e.g., "Psi").
 * @param[in]     dt          Time step size.
 * @param[in,out] psi_io      Pointer to the particle's scalar value (Psi).
 * @param[in]     diffusivity Particle diffusivity (Gamma + Gamma_t).
 * @param[in]     mean_val    Local Eulerian mean value (<Psi>).
 * @param[in]     cell_vol    Volume of the host cell (1/Jacobian).
 * @param[in]     C_model     Model constant (C_IEM).
 */
PetscErrorCode UpdateParticleField(const char *fieldName,
                                   PetscReal dt,
                                   PetscReal *psi_io,
                                   PetscReal diffusivity,
                                   PetscReal mean_val,
                                   PetscReal cell_vol,
                                   PetscReal C_model)
{
    PetscFunctionBeginUser;

    if (strcmp(fieldName, "Psi") == 0) {
        // --- IEM Mixing Model ---
        
        // 1. Calculate Characteristic Length Scale (Delta)
        // Delta^2 approx Volume^(2/3)
        // Safety: Ensure volume is positive
        if (cell_vol < 1.0e-14) cell_vol = 1.0e-14;
        
        PetscReal delta2 = PetscPowReal(cell_vol, 0.6666667); 

        // 2. Calculate Mixing Frequency (Omega)
        // Omega = C_phi * (Gamma_eff) / Delta^2
        PetscReal omega = C_model * diffusivity / delta2;

        // 3. Analytical Integration
        // Exact solution for linear decay ODE
        PetscReal decay = PetscExpReal(-omega * dt);
        
        // Update Value IN-PLACE
        PetscReal psi_old = *psi_io;
        *psi_io = mean_val + (psi_old - mean_val) * decay;
    }
    else {
        // Logic for other fields can be added here
        // e.g. Temperature might use a different C_model or source term
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdateFieldForAllParticles"
/**
 * @brief Loops over all local particles and updates a specified field.
 *
 * Prepares necessary Eulerian and Lagrangian data structures before looping.
 *
 * @param[in,out] user      Pointer to the UserCtx.
 * @param[in]     fieldName The name of the field to update (e.g., "Psi").
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode UpdateFieldForAllParticles(UserCtx *user, const char *fieldName)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    DM             da = user->da; // Scalar DM for accessing lPsi and lAj
    PetscInt       n_local;
    PetscReal      dt = user->simCtx->dt;
    
    // Model Constant (Default to 2.0 if not set)
    PetscReal      C_IEM = 2.0;//(user->simCtx->C_IEM > 1.0e-6) ? user->simCtx->C_IEM : 2.0;

    // Pointers for Swarm Data
    PetscReal       *psi_arr = NULL;
    PetscReal       *diff_arr = NULL;
    PetscInt        *cell_arr = NULL;

    // Pointers for Eulerian Data (Read-Only)
    PetscReal       ***grid_mean = NULL;
    PetscReal       ***grid_aj = NULL;

    // Flags to track what we accessed (for cleanup)
    PetscBool       accessed_eulerian = PETSC_FALSE;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    ierr = DMSwarmGetLocalSize(swarm, &n_local); CHKERRQ(ierr);
    if (n_local == 0) {
        PetscFunctionReturn(0); 
    }

    // 1. Access Target Field (Write Access)
    ierr = DMSwarmGetField(swarm, fieldName, NULL, NULL, (void**)&psi_arr); CHKERRQ(ierr);

    // 2. Access Auxiliary Data based on Physics Requirements
    //    For "Psi", we need Diffusivity, CellID, Grid Mean, and Grid Metrics.
    if (strcmp(fieldName, "Psi") == 0) {
        
        // --- Lagrangian Data ---
        ierr = DMSwarmGetField(swarm, "Diffusivity", NULL, NULL, (void**)&diff_arr); CHKERRQ(ierr);
        ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cell_arr); CHKERRQ(ierr);

        // --- Eulerian Data ---
        // We need Local vectors (lPsi, lAj) to handle ghost indexing safely
        if (!user->lPsi || !user->lAj) {
             SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "UserCtx lPsi or lAj not initialized.");
        }
        
        ierr = DMDAVecGetArrayRead(da, user->lPsi, &grid_mean); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(da, user->lAj,  &grid_aj);   CHKERRQ(ierr);
        accessed_eulerian = PETSC_TRUE;
    }

    // 3. Loop over particles
    for (PetscInt p = 0; p < n_local; ++p) {
        
        // Default values for generic fields
        PetscReal p_diff = 0.0;
        PetscReal p_mean = 0.0;
        PetscReal p_vol  = 1.0;

        if (strcmp(fieldName, "Psi") == 0) {
            // Unpack indices
            PetscInt i = cell_arr[3*p + 0];
            PetscInt j = cell_arr[3*p + 1];
            PetscInt k = cell_arr[3*p + 2];

            // Fetch Data
            p_diff = diff_arr[p];
            p_mean = grid_mean[k][j][i]; // Direct Lookup
            
            // Calculate Volume = 1.0 / Jacobian
            PetscReal jac = grid_aj[k][j][i];
            p_vol = (jac > 1.0e-14) ? (1.0 / jac) : 1.0e-14;
        }

        // Call Kernel
        ierr = UpdateParticleField(fieldName, dt, &psi_arr[p], p_diff, p_mean, p_vol, C_IEM); 
        CHKERRQ(ierr);
    }

    // 4. Restore Arrays
    ierr = DMSwarmRestoreField(swarm, fieldName, NULL, NULL, (void**)&psi_arr); CHKERRQ(ierr);

    if (strcmp(fieldName, "Psi") == 0) {
        ierr = DMSwarmRestoreField(swarm, "Diffusivity", NULL, NULL, (void**)&diff_arr); CHKERRQ(ierr);
        ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cell_arr); CHKERRQ(ierr);
    }

    if (accessed_eulerian) {
        ierr = DMDAVecRestoreArrayRead(da, user->lPsi, &grid_mean); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(da, user->lAj,  &grid_aj);   CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Updated particle physics for field '%s'.\n", fieldName);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdateAllParticleFields"
PetscErrorCode UpdateAllParticleFields(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Updating all particle physical properties...\n");

    // Update Scalar Mixing (IEM)
    // IMPORTANT: Ensure ScatterParticleFieldToEulerField("Psi") has been called 
    // immediately before this to populate user->lPsi with the current mean.
    ierr = UpdateFieldForAllParticles(user, "Psi"); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "All particle physical properties updated.\n");
    
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}   