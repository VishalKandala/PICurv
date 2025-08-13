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
#include "AnalyticalSolution.h" // Analytical Solution for testing
#include "ParticleMotion.h" // Functions related to motion of particles
#include "Boundaries.h"     //  Functions related to Boundary condition
#include "simulation.h"

/**
 * @brief Sets the initial values for the INTERIOR of a specified Eulerian field.
 *
 * This function initializes the interior nodes of `Ucont` based on a profile selected
 * by `user->FieldInitialization`. It explicitly skips any node that lies on a global
 * boundary, as those values are set by the Boundary System's `Initialize` methods.
 *
 * The initialization is directional, aligned with the primary INLET face that was
 * identified by the parser. This ensures the initial flow is physically meaningful.
 *
 * Supported `user->FieldInitialization` profiles for "Ucont":
 *  - 0: Zero Velocity. All interior components of Ucont are set to 0.
 *  - 1: Constant Normal Velocity. The contravariant velocity component normal to the
 *       inlet direction is set such that the physical velocity normal to those grid
 *       planes is a constant `uin`. Other contravariant components are zero.
 *  - 2: Poiseuille Normal Velocity. The contravariant component normal to the
 *       inlet direction is set with a parabolic profile.
 *
 * @param user      The main UserCtx struct, containing all simulation data and configuration.
 * @param fieldName A string ("Ucont" or "P") identifying which field to initialize.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SetInitialInteriorField(UserCtx *user, const char *fieldName);

/**
 * @brief High-level orchestrator to set the complete initial state of the Eulerian solver.
 *
 * This function is called once from main() before the time loop begins. It inspects
 * the simulation context to determine whether to perform a fresh start (t=0) or
 * restart from saved files. It then delegates to the appropriate helper function.
 * Finally, it initializes the solver's history vectors (Ucont_o, P_o, etc.)
 * to ensure the first time step has the necessary data.
 */
PetscErrorCode InitializeEulerianState(SimCtx *simCtx);

#endif // INITIALCONDITION_H
