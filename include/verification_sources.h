#ifndef VERIFICATION_SOURCES_H
#define VERIFICATION_SOURCES_H

#include "variables.h"

/**
 * @brief Reports whether a verification-only diffusivity override is active.
 *
 * @param simCtx Simulation context.
 * @return PETSC_TRUE when a verification diffusivity source override is enabled.
 */
PetscBool VerificationDiffusivityOverrideActive(const SimCtx *simCtx);

/**
 * @brief Populates the Eulerian diffusivity field from a verification-only source override.
 *
 * @param[in,out] user Block-local context whose `Diffusivity` vector will be filled.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ApplyVerificationDiffusivityOverride(UserCtx *user);

#endif /* VERIFICATION_SOURCES_H */
