#ifndef RUNLOOP_H
#define RUNLOOP_H

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
#include "AnalyticalSolutions.h"
#include "ParticleMotion.h" // Functions related to motion of particles
#include "ParticlePhysics.h" // Functions related to particle scalar/state updates
#include "Boundaries.h"     //  Functions related to Boundary conditions
#include "setup.h"          // Functions  related to setup
#include "solvers.h"

/**
 * @brief Installs lightweight signal handlers for graceful shutdown requests.
 *
 * The handlers only record that a shutdown signal was received. The actual
 * output flush and exit path happens later at safe checkpoints in the run loop.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InitializeRuntimeSignalHandlers(void);

/**
 * @brief Update an EWMA estimate for timestep wall-clock duration.
 *
 * @param[in] has_previous          Whether a previous EWMA estimate exists.
 * @param[in] previous_ewma_seconds Prior EWMA estimate in seconds.
 * @param[in] latest_step_seconds   Latest completed timestep duration in seconds.
 * @param[in] alpha                 EWMA weighting factor in `(0, 1]`.
 * @return PetscReal Updated EWMA estimate in seconds.
 */
PetscReal RuntimeWalltimeGuardUpdateEWMA(PetscBool has_previous, PetscReal previous_ewma_seconds, PetscReal latest_step_seconds, PetscReal alpha);

/**
 * @brief Return the conservative timestep estimate used by the walltime guard.
 *
 * @param[in] warmup_average_seconds Average duration across warmup steps.
 * @param[in] ewma_seconds           Current EWMA duration estimate.
 * @param[in] latest_step_seconds    Most recent completed timestep duration.
 * @return PetscReal Conservative timestep estimate in seconds.
 */
PetscReal RuntimeWalltimeGuardConservativeEstimate(PetscReal warmup_average_seconds, PetscReal ewma_seconds, PetscReal latest_step_seconds);

/**
 * @brief Compute the required shutdown headroom from timestep estimate and floor.
 *
 * @param[in] min_seconds                   Absolute minimum shutdown headroom.
 * @param[in] multiplier                    Safety multiplier applied to the timestep estimate.
 * @param[in] conservative_estimate_seconds Conservative timestep estimate in seconds.
 * @return PetscReal Required headroom in seconds.
 */
PetscReal RuntimeWalltimeGuardRequiredHeadroom(PetscReal min_seconds, PetscReal multiplier, PetscReal conservative_estimate_seconds);

/**
 * @brief Decide whether the runtime walltime guard should stop before another step.
 *
 * @param[in]  completed_steps                Number of completed timesteps observed so far.
 * @param[in]  warmup_steps                   Minimum completed timesteps required before guarding.
 * @param[in]  remaining_seconds              Remaining walltime in seconds.
 * @param[in]  min_seconds                    Absolute minimum shutdown headroom.
 * @param[in]  multiplier                     Safety multiplier applied to timestep estimate.
 * @param[in]  warmup_average_seconds         Average duration across warmup steps.
 * @param[in]  ewma_seconds                   Current EWMA duration estimate.
 * @param[in]  latest_step_seconds            Latest completed timestep duration.
 * @param[out] required_headroom_seconds_out  Computed required headroom in seconds.
 * @return PetscBool `PETSC_TRUE` when shutdown should be requested.
 */
PetscBool RuntimeWalltimeGuardShouldTrigger(PetscInt completed_steps, PetscInt warmup_steps, PetscReal remaining_seconds, PetscReal min_seconds, PetscReal multiplier, PetscReal warmup_average_seconds, PetscReal ewma_seconds, PetscReal latest_step_seconds, PetscReal *required_headroom_seconds_out);
 
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
 * @brief Executes the main time-marching loop for the particle simulation.
 *
 * This version uses the new, integrated `LocateAllParticlesInGrid` orchestrator
 * and the `ResetAllParticleStatuses` helper for a clean, robust, and understandable workflow.
 *
 * For each timestep, it performs:
 *  1. Sets the background fluid velocity field (Ucat) for the current step.
 *  2. Updates particle positions using velocity from the *previous* step's interpolation.
 *  3. Removes any particles that have left the global domain.
 *  4. A single call to `LocateAllParticlesInGrid', which handles all
 *     particle location and migration until the swarm is fully settled.
 *  5. Interpolates the current fluid velocity to the newly settled particle locations.
 *  6. Scatters particle data back to Eulerian fields.
 *  7. Outputs data at specified intervals.
 *
 * @param simCtx     Pointer to the master simulation context.
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
PetscErrorCode PerformInitializedParticleSetup(SimCtx *simCtx);

/**
 * @brief Finalizes the simulation state after particle and fluid data have been loaded from a restart.
 *
 * This helper function performs the critical sequence of operations required to ensure
 * the loaded Lagrangian and Eulerian states are fully consistent and the solver is
 * ready to proceed. This includes:
 * 1. Verifying particle locations in the grid and building runtime links.
 * 2. Synchronizing particle velocity with the authoritative grid velocity via interpolation.
 * 3. Scattering particle source terms (e.g., volume fraction) back to the grid.
 * 4. Updating the solver's history vectors with the final, fully-coupled state.
 * 5. Writing the complete, consistent state to output files for the restart step.
 *
 * @param simCtx The main simulation context.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode PerformLoadedParticleSetup(SimCtx *simCtx);

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

#endif  // RUNLOOP_H
