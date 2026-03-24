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

    if (!VerificationDiffusivityOverrideActive(simCtx)) PetscFunctionReturn(0);

    ierr = DMDAVecGetArray(user->da, user->Diffusivity, &diff_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->Cent, &cent); CHKERRQ(ierr);

    for (PetscInt k = info.zs; k < info.zs + info.zm; ++k) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; ++j) {
            for (PetscInt i = info.xs; i < info.xs + info.xm; ++i) {
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
        if (global_min_gamma <= 0.0) {
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                    "verification diffusivity override produced non-positive Gamma values (min=%g).",
                    (double)global_min_gamma);
        }
    }

    ierr = UpdateLocalGhosts(user, "Diffusivity"); CHKERRQ(ierr);
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
    if (!VerificationScalarOverrideActive(user->simCtx)) PetscFunctionReturn(0);

    PetscCall(SetAnalyticalScalarFieldOnParticles(user, "Psi"));
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "Applied verification scalar override profile '%s' to particle field 'Psi'.\n",
              user->simCtx->verificationScalar.profile);
    PetscFunctionReturn(0);
}
