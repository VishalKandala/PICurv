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
PetscErrorCode UpdateParticleField(ParticleFieldId field_id,
                                   PetscReal dt,
                                   PetscReal *psi_io,
                                   PetscReal diffusivity,
                                   PetscReal mean_val,
                                   PetscReal cell_vol,
                                   PetscReal C_model)
{
    PetscFunctionBeginUser;

    if (field_id == PARTICLE_FIELD_ID_PSI) {
        // Guard the LES mixing time scale against degenerate or cut-cell volumes.
        if (cell_vol < 1.0e-14) cell_vol = 1.0e-14;

        // The IEM model relaxes Psi exponentially toward the cell mean over dt.
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
PetscErrorCode UpdateFieldForAllParticles(UserCtx *user, ParticleFieldId field_id)
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
    const ParticleFieldDescriptor *descriptor = NULL;
    const char      *fieldName = NULL;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    ierr = ParticleFieldGetDescriptor(field_id, &descriptor); CHKERRQ(ierr);
    PetscCheck((descriptor->capabilities & PARTICLE_FIELD_CAPABILITY_MODEL_UPDATE) != 0,
               PETSC_COMM_SELF, PETSC_ERR_SUP,
               "Particle field '%s' has no registered model-update kernel.",
               descriptor->canonical_name);
    fieldName = descriptor->canonical_name;

    ierr = DMSwarmGetLocalSize(swarm, &n_local); CHKERRQ(ierr);
    if (n_local == 0) {
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    ierr = DMSwarmGetField(swarm, fieldName, NULL, NULL, (void**)&psi_arr); CHKERRQ(ierr);

    if (field_id == PARTICLE_FIELD_ID_PSI) {
        ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_DIFFUSIVITY), NULL, NULL, (void**)&diff_arr); CHKERRQ(ierr);
        ierr = DMSwarmGetField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_CELL_ID), NULL, NULL, (void**)&cell_arr); CHKERRQ(ierr);

        // Psi relaxation requires ghosted Eulerian mean and Jacobian fields.
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

        if (field_id == PARTICLE_FIELD_ID_PSI) {
            PetscInt i = cell_arr[3*p + 0];
            PetscInt j = cell_arr[3*p + 1];
            PetscInt k = cell_arr[3*p + 2];

            p_diff = diff_arr[p];
            p_mean = grid_mean[k][j][i];

            // Aj is the reciprocal cell volume in the curvilinear-grid representation.
            PetscReal jac = grid_aj[k][j][i];
            p_vol = (jac > 1.0e-14) ? (1.0 / jac) : 1.0e-14;
        }

        ierr = UpdateParticleField(field_id, dt, &psi_arr[p], p_diff, p_mean, p_vol, C_IEM);
        CHKERRQ(ierr);
    }

    ierr = DMSwarmRestoreField(swarm, fieldName, NULL, NULL, (void**)&psi_arr); CHKERRQ(ierr);

    if (field_id == PARTICLE_FIELD_ID_PSI) {
        ierr = DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_DIFFUSIVITY), NULL, NULL, (void**)&diff_arr); CHKERRQ(ierr);
        ierr = DMSwarmRestoreField(swarm, ParticleFieldName(PARTICLE_FIELD_ID_CELL_ID), NULL, NULL, (void**)&cell_arr); CHKERRQ(ierr);
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
        // Verification profiles define Psi exactly, so bypass the model-driven update.
        ierr = ApplyVerificationScalarOverrideToParticles(user); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL, LOG_INFO, "Verification scalar override active; skipped model-driven Psi update.\n");
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    ierr = UpdateFieldForAllParticles(user, PARTICLE_FIELD_ID_PSI); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "All particle physical properties updated.\n");

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
