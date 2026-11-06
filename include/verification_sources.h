#ifndef VERIFICATION_SOURCES_H
#define VERIFICATION_SOURCES_H

#include "variables.h"
#include "particle_field_catalog.h"

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

/**
 * @brief Reports whether a verification-only scalar override is active.
 *
 * @param simCtx Simulation context.
 * @return PETSC_TRUE when a verification scalar source override is enabled.
 */
PetscBool VerificationScalarOverrideActive(const SimCtx *simCtx);

/**
 * @brief Populates the particle `Psi` field from a verification-only source override.
 *
 * @param[in,out] user Block-local context whose particle scalar field will be filled.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ApplyVerificationScalarOverrideToParticles(UserCtx *user);

/**
 * @brief Compiles a verification field's expression, in physical x, y, z and t, and enables the field.
 * @param[in]     text  Expression text.
 * @param[in,out] field Field to set; a previous expression is released.
 * @return PETSc error code.
 */
PetscErrorCode VerificationFieldCompile(const char *text, VerificationFieldConfig *field);

/**
 * @brief Reads and compiles the configured verification fields (`-verification_*_expr`).
 * @param simCtx Simulation context.
 * @return PETSc error code; an override outside an analytical run is refused.
 */
PetscErrorCode VerificationFieldsCreate(SimCtx *simCtx);

/**
 * @brief Releases the compiled verification fields.
 * @param simCtx Simulation context.
 * @return PETSc error code.
 */
PetscErrorCode VerificationFieldsDestroy(SimCtx *simCtx);

/**
 * @brief Evaluates a verification field at a point.
 * @param[in]  simCtx      Simulation context, for the scales.
 * @param[in]  field       Enabled verification field.
 * @param[in]  point       Solver-unit position.
 * @param[in]  solver_time Solver-unit time.
 * @param[out] value       The expression's physical value.
 * @return PETSc error code.
 */
PetscErrorCode VerificationFieldEvaluate(const SimCtx *simCtx, const VerificationFieldConfig *field, Cmpnts point,
                                         PetscReal solver_time, PetscReal *value);

/**
 * @brief Writes the configured verification scalar profile onto a particle swarm scalar field.
 *
 * @param[in,out] user             Block-local context providing particle positions.
 * @param[in]     particle_field_id Typed identity of the scalar swarm field to overwrite.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SetVerificationScalarOnParticles(UserCtx *user, ParticleFieldId particle_field_id);

/**
 * @brief Writes the configured verification scalar profile at physical cell centers into a scalar Vec.
 *
 * @param[in]     user      Block-local context providing the cell-center coordinates.
 * @param[in,out] targetVec Scalar Vec on `user->da` storage to fill with analytical reference values.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SetVerificationScalarAtCellCenters(UserCtx *user, Vec targetVec);

#endif /* VERIFICATION_SOURCES_H */
