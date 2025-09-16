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
 * @brief Finalizes the simulation setup at t=0, ensuring a consistent state before time marching.
 *
 * This function is called from main() after the initial Eulerian and Lagrangian states have been
 * created but before the main time loop begins. Its responsibilities are:
 *
 * 1.  Settling the particle swarm: Migrates particles to their correct owner ranks and finds their
 *     initial host cells. This includes handling special surface initializations.
 * 2.  Coupling the fields: Interpolates the initial Eulerian fields to the settled particle locations.
 * 3.  Preparing for the first step: Scatters particle data back to the grid.
 * 4.  Writing the initial output for step 0.
 *
 * @param simCtx Pointer to the main simulation context structure.
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode PerformInitialSetup(SimCtx *simCtx);

/**
 * @brief Performs post-load/post-init consistency checks for a restarted simulation.
 *
 * This function is called from main() ONLY when a restart is being performed
 * (i.e., StartStep > 0). It inspects the particle restart mode to determine the
 * correct finalization procedure for the Lagrangian swarm.
 *
 * - If particles were loaded from a file (`mode == "load"`), it verifies their
 *   locations within the grid to establish necessary runtime links.
 * - If new particles were initialized into the restarted flow (`mode == "init"`),
 *   it runs the full `PerformInitialSetup` sequence to migrate, locate, and
 *   couple the new particles with the existing fluid state.
 *
 * @param simCtx The main simulation context.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode FinalizeRestartState(SimCtx *simCtx);

#endif  // SIMULATION_H
