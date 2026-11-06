#ifndef ANALYTICALSOLUTIONS_H
#define ANALYTICALSOLUTIONS_H

#include "variables.h"

/**
 * @brief Sets the Eulerian fields from the case's source (`-ic_source_*`) at the step's time.
 * @details Sets the fluxes on every face, `Ubcs` on every boundary face, and `Ucat` and `P` at the
 *          cell centres. A source that does not read `t` is evaluated once. Collective.
 * @param simCtx The main simulation context, containing configuration and state.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode AnalyticalSolutionEngine(SimCtx *simCtx);

/**
 * @brief Release what AnalyticalSolutionEngine() keeps from step to step, if anything.
 * @param simCtx Simulation context.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode AnalyticalFieldDestroy(SimCtx *simCtx);

/**
 * @brief Reports whether the analytical field can be compared with interpolated particle velocities.
 * @details It must be evaluable at any point, a Cartesian source reading only x, y, z and t, and
 *          must vary in space: interpolating a uniform field has no error to measure.
 * @param simCtx Simulation context, after AnalyticalSolutionEngine().
 * @return PETSC_TRUE if interpolation error is meaningful, PETSC_FALSE otherwise.
 */
PetscBool AnalyticalFieldSupportsInterpolationError(const SimCtx *simCtx);

/**
 * @brief Evaluates the analytical velocity at particle positions.
 * @param tempVec The PETSc Vec containing particle positions, replaced by the velocities.
 * @param simCtx The simulation context.
 * @return PetscErrorCode Returns 0 on success.
 */
PetscErrorCode SetAnalyticalSolutionForParticles(Vec tempVec, SimCtx *simCtx);

#endif // ANALYTICALSOLUTIONS_H
