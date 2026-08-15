#include "verification_sources.h"

#include "AnalyticalSolutions.h"
#include "logging.h"
#include "setup.h"

/**
 * @brief Reports whether the verification-only diffusivity override is enabled.
 * @details See `include/verification_sources.h` for the public parameter and return contract.
 */
PetscBool VerificationDiffusivityOverrideActive(const SimCtx *simCtx)
{
    if (!simCtx) return PETSC_FALSE;
    return simCtx->verificationDiffusivity.enabled;
}

/**
 * @brief Reports whether the verification-only scalar override is enabled.
 * @details See `include/verification_sources.h` for the public parameter and return contract.
 */
PetscBool VerificationScalarOverrideActive(const SimCtx *simCtx)
{
    if (!simCtx) return PETSC_FALSE;
    return simCtx->verificationScalar.enabled;
}

/**
 * @brief Fills the Eulerian diffusivity field from the configured verification profile.
 * @details See `include/verification_sources.h` for the public parameter and return contract.
 */
PetscErrorCode ApplyVerificationDiffusivityOverride(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx        *simCtx = user->simCtx;
    PetscReal    ***diff_arr = NULL;
    const Cmpnts ***cent = NULL;
    DMDALocalInfo  info = user->info;
    PetscReal      min_gamma = PETSC_MAX_REAL;

    PetscFunctionBeginUser;

    // Leave the production diffusivity untouched unless this verification mode is selected.
    if (!VerificationDiffusivityOverrideActive(simCtx)) PetscFunctionReturn(0);

    ierr = DMDAVecGetArray(user->da, user->Diffusivity, &diff_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->Cent, &cent); CHKERRQ(ierr);

    for (PetscInt k = info.zs; k < info.zs + info.zm; ++k) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; ++j) {
            for (PetscInt i = info.xs; i < info.xs + info.xm; ++i) {
                // The manufactured profile is affine in the physical x coordinate.
                const PetscReal gamma =
                  simCtx->verificationDiffusivity.gamma0 +
                  simCtx->verificationDiffusivity.slope_x * cent[k][j][i].x;
                diff_arr[k][j][i] = gamma;
                min_gamma = PetscMin(min_gamma, gamma);
            }
        }
    }

    ierr = DMDAVecRestoreArrayRead(user->fda, user->Cent, &cent); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, user->Diffusivity, &diff_arr); CHKERRQ(ierr);

    {
        PetscReal global_min_gamma = 0.0;
        ierr = MPI_Allreduce(&min_gamma, &global_min_gamma, 1, MPIU_REAL, MPI_MIN, PETSC_COMM_WORLD); CHKERRMPI(ierr);
        // Diffusivity must remain positive on every rank for the PDE operator to be valid.
        if (global_min_gamma <= 0.0) {
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                    "verification diffusivity override produced non-positive Gamma values (min=%g).",
                    (double)global_min_gamma);
        }
    }

    {
        const FieldId periodic_fields[] = {FIELD_ID_DIFFUSIVITY};
        // Make the analytic field coherent across periodic seams before ghost updates.
        ierr = SynchronizePeriodicCellFields(user, 1, periodic_fields); CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, FIELD_ID_DIFFUSIVITY); CHKERRQ(ierr);
    }
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "Applied verification diffusivity override profile '%s' (gamma0=%g, slope_x=%g).\n",
              simCtx->verificationDiffusivity.profile,
              (double)simCtx->verificationDiffusivity.gamma0,
              (double)simCtx->verificationDiffusivity.slope_x);

    PetscFunctionReturn(0);
}

/**
 * @brief Fills the particle `Psi` field from the configured verification scalar profile.
 * @details See `include/verification_sources.h` for the public parameter and return contract.
 */
PetscErrorCode ApplyVerificationScalarOverrideToParticles(UserCtx *user)
{
    PetscFunctionBeginUser;
    if (!user) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx cannot be NULL.");
    // The normal particle model owns Psi unless the verification override is enabled.
    if (!VerificationScalarOverrideActive(user->simCtx)) PetscFunctionReturn(0);

    PetscCall(SetAnalyticalScalarFieldOnParticles(user, PARTICLE_FIELD_ID_PSI));
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "Applied verification scalar override profile '%s' to particle field 'Psi'.\n",
              user->simCtx->verificationScalar.profile);
    PetscFunctionReturn(0);
}
