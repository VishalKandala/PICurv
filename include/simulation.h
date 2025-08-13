#ifndef SIMULATION_H
#define SIMULATION_H

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
#include "initialcondition.h" // Analytical Solution for testing
#include "ParticleMotion.h" // Functions related to motion of particles
#include "Boundaries.h"     //  Functions related to Boundary conditions
#include "setup.h"          // Functions  related to setup
#include "solvers.h"


/**
 * @brief Copies the current time step's solution fields into history vectors
 *        (e.g., U(t_n) -> U_o, U_o -> U_rm1) for the next time step's calculations.
 *
 * This function is critical for multi-step time integration schemes (like BDF2)
 * used by the legacy solver. It must be called at the end of every time step,
 * after the new solution has been fully computed.
 *
 * The order of operations is important to avoid overwriting data prematurely.
 *
 * @param user The UserCtx for a single block. The function modifies the history
 *             vectors (Ucont_o, Ucont_rm1, etc.) within this context.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode UpdateSolverHistoryVectors(UserCtx *user);

/**
 * @brief Initializes or updates the complete, consistent state of all Eulerian fields for a given timestep.
 *
 * This function is a high-level wrapper that orchestrates the entire process of preparing
 * the fluid fields for a single time step. It follows the standard procedure for a
 * curvilinear solver: first resolving contravariant velocities (`Ucont`) and then
 * converting them to Cartesian (`Ucat`).
 *
 * Its sequential operations are:
 * 1.  Update the INTERIOR of the domain:
 *     - For the initial step, it calls `SetInitialInteriorField` to generate values.
 *     - For subsequent steps, it calls the main fluid solver.
 *     - If restarting from a file, it reads the data, overwriting the whole field.
 *
 * 2.  Apply Boundary Conditions:
 *     - It then calls the modular `BoundarySystem_ExecuteStep` to enforce all configured
 *       boundary conditions on the domain edges.
 *
 * 3.  Convert to Cartesian and Finalize:
 *     - It calls `Contra2Cart` to compute `Ucat` from `Ucont`.
 *     - It calls `UpdateLocalGhosts` to ensure all parallel data is synchronized.
 *
 * @param user        Pointer to the UserCtx structure, containing all simulation data.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SetEulerianFields(UserCtx *user);

/**
 * @brief Executes the main time-marching loop for the particle simulation. [TEST VERSION]
 *
 * This version uses the new, integrated `LocateAllParticlesInGrid_TEST` orchestrator
 * and the `ResetAllParticleStatuses` helper for a clean, robust, and understandable workflow.
 *
 * For each timestep, it performs:
 *  1. Sets the background fluid velocity field (Ucat) for the current step.
 *  2. Updates particle positions using velocity from the *previous* step's interpolation.
 *  3. Removes any particles that have left the global domain.
 *  4. A single call to `LocateAllParticlesInGrid_TEST`, which handles all
 *     particle location and migration until the swarm is fully settled.
 *  5. Interpolates the current fluid velocity to the newly settled particle locations.
 *  6. Scatters particle data back to Eulerian fields.
 *  7. Outputs data at specified intervals.
 *
 * @param user       Pointer to the UserCtx structure..
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode AdvanceSimulation(SimCtx *simCtx);

/**
 * @brief Performs the complete initial setup for the particle simulation at time t=0. [TEST VERSION]
 *
 * This version uses the new, integrated `LocateAllParticlesInGrid_TEST` orchestrator,
 * which handles both location and migration in a single, robust, iterative process.
 *
 * Its sequential operations are:
 * 1. A single, comprehensive call to `LocateAllParticlesInGrid_TEST` to sort all particles
 *    to their correct owner ranks and find their initial host cells.
 * 2. If `user->ParticleInitialization == 0` (Surface Init), it re-initializes particles on the
 *    designated inlet surface, now that they are on the correct MPI ranks.
 * 3. A second call to `LocateAllParticlesInGrid_TEST` is needed after re-initialization to
 *    find the new, correct host cells for the surface-placed particles.
 * 4. Interpolates initial Eulerian fields to the settled particles.
 * 5. Scatters particle data to Eulerian fields (if applicable).
 * 6. Outputs initial data if requested.
 *
 * @param user Pointer to the UserCtx structure.
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode PerformInitialSetup(UserCtx *user, BoundingBox *bboxlist);

#endif  // SIMULATION_H
