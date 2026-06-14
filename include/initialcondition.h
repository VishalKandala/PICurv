#ifndef INITIALCONDITION_H
#define INITIALCONDITION_H

#include <petscpf.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <petsctime.h>
#include <petscsys.h>
#include <petscdmcomposite.h>
#include <petscsystypes.h>
// Include additional headers
#include "variables.h"         // Shared type definitions
#include "ParticleSwarm.h"  // Particle swarm functions
#include "walkingsearch.h"  // Particle location functions
#include "grid.h"           // Grid functions
#include "logging.h"        // Logging macros
#include "io.h"             // Data Input and Output functions
#include "interpolation.h"  // Interpolation routines
#include "AnalyticalSolutions.h" // Analytical Solution for testing
#include "ParticleMotion.h" // Functions related to motion of particles
#include "Boundaries.h"     //  Functions related to Boundary condition
#include "runloop.h"

/**
 * @brief Sets the initial values for the INTERIOR of a specified Eulerian field.
 *
 * This function initializes the interior nodes of `Ucont` based on the mode selected
 * by `simCtx->initialConditionMode`.
 *
 * Supported profiles for "Ucont":
 *  - IC_MODE_ZERO: All interior contravariant components are set to zero.
 *  - IC_MODE_CONSTANT_CARTESIAN: `UniformCart2Contra` dots the Cartesian
 *       vector `(InitialConstantContra.x/y/z)` with the local metric vectors to fill
 *       all three contravariant components correctly across the entire interior.
 *  - IC_MODE_CONSTANT_STREAMWISE: Sets only the
 *       contravariant component along the streamwise axis (from `flowDirection` or the
 *       identified INLET face) proportional to `icVelocityPhysical * |A_n|`.
 *  - IC_MODE_POISEUILLE: Separable parabolic profile in the two cross-stream index directions;
 *       centerline speed is `icVelocityPhysical`; streamwise axis from `flowDirection`
 *       or the identified INLET face.
 *
 * @param user      The main UserCtx struct, containing all simulation data and configuration.
 * @param fieldName A string ("Ucont" or "P") identifying which field to initialize.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SetInitialInteriorField(UserCtx *user, const char *fieldName);

/**
 * @brief Populate Ucont for one fresh-start block from the configured IC mode.
 *
 * Built-in modes generate Ucont directly. File mode loads either Ucat or Ucont
 * from the staged IC directory and converts Ucat when necessary.
 * @param[in,out] user Block context whose velocity field is populated.
 * @return PETSc error code.
 */
PetscErrorCode PopulateInitialUcont(UserCtx *user);

/**
 * @brief High-level orchestrator to set the complete initial state of the Eulerian solver.
 *
 * This function is called once from main() before the time loop begins. It inspects
 * the simulation context to determine whether to perform a fresh start (t=0) or
 * restart from saved files. It then delegates to the appropriate helper function.
 * Finally, it initializes the solver's history vectors (Ucont_o, P_o, etc.)
 * to ensure the first time step has the necessary data.
 *
 * @param simCtx Simulation context controlling the operation.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InitializeEulerianState(SimCtx *simCtx);

#endif // INITIALCONDITION_H
