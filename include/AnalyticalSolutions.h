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
 * @brief Sets the grid domain and resolution for analytical solution cases.
 *
 * @details This function is called when `eulerianSource` is "analytical". It is responsible for
 * automatically configuring the grid based on the chosen `AnalyticalSolutionType`.
 *
 * @par TGV3D Multi-Block Decomposition
 * If the analytical solution is "TGV3D", this function automatically decomposes the
 * required `[0, 2*PI]` physical domain among the available blocks.
 * - **Single Block (`nblk=1`):** The single block is assigned the full `[0, 2*PI]` domain.
 * - **Multiple Blocks (`nblk>1`):** It requires that the number of blocks be a **perfect square**
 *   (e.g., 4, 9, 16). It then arranges the blocks in a `sqrt(nblk)` by `sqrt(nblk)` grid in the
 *   X-Y plane, partitioning the `[0, 2*PI]` domain in X and Y accordingly. The Z domain for all
 *   blocks remains `[0, 2*PI]`. If `nblk` is not a perfect square, the simulation is aborted
 *   with an error.
 *
 * After setting the domain bounds, it proceeds to read the grid resolution options
 * (`-im`, `-jm`, `-km`) from the command line for the specific block.
 *
 * @param user Pointer to the `UserCtx` for a specific block. The function will
 *             populate the geometric fields (`IM`, `JM`, `KM`, `Min_X`, `Max_X`, etc.)
 *             within this struct.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode SetAnalyticalGridInfo(UserCtx *user);

#undef __FUNCT__
#define __FUNCT__ "AnalyticalSolutionEngine"
/**
 * @brief Dispatches to the appropriate analytical solution function based on simulation settings.
 *
 * This function acts as a router. It reads the `AnalyticalSolutionType` from the
 * simulation context and calls the corresponding private implementation function
 * (e.g., for Taylor-Green Vortex, lid-driven cavity, etc.). This design keeps
 * the main simulation code clean and makes it easy to add new analytical test cases.
 *
 * @param simCtx The main simulation context, containing configuration and state.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode AnalyticalSolutionEngine(SimCtx *simCtx);

/**
@brief Applies the analytical solution to particle velocity vector.

@details Dispatcher function that calls the appropriate analytical solution based on 
         simCtx->AnalyticalSolutionType. Supports multiple solution types.

@param tempVec The PETSc Vec containing particle positions which will be used to store velocities.
@param simCtx The simulation context.
@return PetscErrorCode Returns 0 on success.
*/
PetscErrorCode SetAnalyticalSolutionForParticles(Vec tempVec, SimCtx *simCtx);

#endif // ANALYTICALSOLUTIONS_H