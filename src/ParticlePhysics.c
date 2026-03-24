#include "ParticlePhysics.h"
#include "verification_sources.h"

#ifndef ERROR_MSG_BUFFER_SIZE
#define ERROR_MSG_BUFFER_SIZE 256
#endif

#undef __FUNCT__
#define __FUNCT__ "UpdateParticleField"
/**
 * @brief Internal helper implementation: `UpdateParticleField()`.
 * @details Local to this translation unit.
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
        if (cell_vol < 1.0e-14) cell_vol = 1.0e-14;

        PetscReal delta2 = PetscPowReal(cell_vol, 0.6666667);
        PetscReal omega = C_model * diffusivity / delta2;
        PetscReal decay = PetscExpReal(-omega * dt);

        PetscReal psi_old = *psi_io;
        *psi_io = mean_val + (psi_old - mean_val) * decay;
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdateFieldForAllParticles"
/**
 * @brief Internal helper implementation: `UpdateFieldForAllParticles()`.
 * @details Local to this translation unit.
 */
PetscErrorCode UpdateFieldForAllParticles(UserCtx *user, const char *fieldName)
{
    PetscErrorCode ierr;
    DM             swarm = user->swarm;
    DM             da = user->da;
    PetscInt       n_local;
    PetscReal      dt = user->simCtx->dt;
    PetscReal      C_IEM = 2.0;

    PetscReal       *psi_arr = NULL;
    PetscReal       *diff_arr = NULL;
    PetscInt        *cell_arr = NULL;

    PetscReal       ***grid_mean = NULL;
    PetscReal       ***grid_aj = NULL;

    PetscBool       accessed_eulerian = PETSC_FALSE;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    ierr = DMSwarmGetLocalSize(swarm, &n_local); CHKERRQ(ierr);
    if (n_local == 0) {
        PetscFunctionReturn(0);
    }

    ierr = DMSwarmGetField(swarm, fieldName, NULL, NULL, (void**)&psi_arr); CHKERRQ(ierr);

    if (strcmp(fieldName, "Psi") == 0) {
        ierr = DMSwarmGetField(swarm, "Diffusivity", NULL, NULL, (void**)&diff_arr); CHKERRQ(ierr);
        ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cell_arr); CHKERRQ(ierr);

        if (!user->lPsi || !user->lAj) {
             SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "UserCtx lPsi or lAj not initialized.");
        }

        ierr = DMDAVecGetArrayRead(da, user->lPsi, &grid_mean); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(da, user->lAj,  &grid_aj);   CHKERRQ(ierr);
        accessed_eulerian = PETSC_TRUE;
    }

    for (PetscInt p = 0; p < n_local; ++p) {
        PetscReal p_diff = 0.0;
        PetscReal p_mean = 0.0;
        PetscReal p_vol  = 1.0;

        if (strcmp(fieldName, "Psi") == 0) {
            PetscInt i = cell_arr[3*p + 0];
            PetscInt j = cell_arr[3*p + 1];
            PetscInt k = cell_arr[3*p + 2];

            p_diff = diff_arr[p];
            p_mean = grid_mean[k][j][i];

            PetscReal jac = grid_aj[k][j][i];
            p_vol = (jac > 1.0e-14) ? (1.0 / jac) : 1.0e-14;
        }

        ierr = UpdateParticleField(fieldName, dt, &psi_arr[p], p_diff, p_mean, p_vol, C_IEM);
        CHKERRQ(ierr);
    }

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
/**
 * @brief Implementation of \ref UpdateAllParticleFields().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/ParticlePhysics.h`.
 * @see UpdateAllParticleFields()
 */
PetscErrorCode UpdateAllParticleFields(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Updating all particle physical properties...\n");

    if (VerificationScalarOverrideActive(user->simCtx)) {
        ierr = ApplyVerificationScalarOverrideToParticles(user); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_INFO, "Verification scalar override active; skipped model-driven Psi update.\n");
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    ierr = UpdateFieldForAllParticles(user, "Psi"); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "All particle physical properties updated.\n");

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
