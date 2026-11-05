#ifndef ANALYTICALSOLUTIONS_H
#define ANALYTICALSOLUTIONS_H

#include <petscpf.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <petsctime.h>
#include <petscdmcomposite.h>

// Include additional headers
#include "variables.h"         // Shared type definitions
#include "ParticleSwarm.h"  // Particle swarm functions
#include "walkingsearch.h"  // Particle location functions
#include "grid.h"           // Grid functions
#include "logging.h"        // Logging macros
#include "io.h"             // Data Input and Output functions
#include "setup.h"          // Setup Module required for array allocation and deallocation
#include "interpolation.h"  // All the different interpolation routines required 

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

/**
 * @brief Evaluates the configured verification scalar profile at one physical point.
 *
 * @param[in]  simCtx Simulation context providing the scalar verification profile.
 * @param[in]  x      Physical x coordinate.
 * @param[in]  y      Physical y coordinate.
 * @param[in]  z      Physical z coordinate.
 * @param[in]  t      Simulation time.
 * @param[out] value  Evaluated scalar value.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode EvaluateAnalyticalScalarProfile(const SimCtx *simCtx,
                                               PetscReal x,
                                               PetscReal y,
                                               PetscReal z,
                                               PetscReal t,
                                               PetscReal *value);

/**
 * @brief Writes the configured verification scalar profile onto a particle swarm scalar field.
 *
 * @param[in,out] user             Block-local context providing particle positions.
 * @param[in]     particle_field_id Typed identity of the scalar swarm field to overwrite.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SetAnalyticalScalarFieldOnParticles(UserCtx *user, ParticleFieldId particle_field_id);

/**
 * @brief Writes the configured verification scalar profile at physical cell centers into a scalar Vec.
 *
 * @param[in]     user      Block-local context providing the cell-center coordinates.
 * @param[in,out] targetVec Scalar Vec on `user->da` storage to fill with analytical reference values.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SetAnalyticalScalarFieldAtCellCenters(UserCtx *user, Vec targetVec);

#endif // ANALYTICALSOLUTIONS_H
