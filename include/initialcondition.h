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
 * @brief Populate Ucont for one fresh-start block from the initial source (`-ic_source_*`).
 * @details A staged field, when there is one, is read first: Ucat, or Ucont. A Cartesian source
 *          is then evaluated at every owned cell centre into Ucat, and a curvilinear one sets
 *          each owned face's flux directly: its component along the grid line, evaluated at the
 *          face centre, times the face area. With `-ic_source_field_weight w` the staged field
 *          is scaled by w. With a bulk velocity, the source is scaled so that the volume mean of
 *          the whole field along the source's own direction equals it. Collective.
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
