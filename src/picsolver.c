/**
 * @file main.c
 * @brief Main driver for the unified CFD and Lagrangian particle solver.
 *
 * This version focuses on establishing the core architecture for a particle-coupled
 * flow simulation. It initializes the grid, solver objects, and a DMSwarm of
 * Lagrangian particles. Hooks for the Immersed Boundary Method (IBM) are present
 * but are currently inactive to allow for iterative development.
 *
 * ARCHITECTURE:
 * 1. A central `SimulationContext` (simCtx) holds all simulation-wide configuration.
 * 2. `CreateSimulationContext()` handles all command-line parsing.
 * 3. The `main()` function provides a high-level workflow for the simulation.
 */

#include "simulation.h" // The new, unified header file

/*================================================================================*
 *                           MAIN FUNCTION                                        *
 *================================================================================*/

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    PetscErrorCode    ierr;
    SimCtx *simCtx = NULL; // The single, top-level context object

    // === I. INITIALIZE =======================================================
    ierr = PetscInitialize(&argc, &argv, (char *)0, "PIC-Solver"); CHKERRQ(ierr);

    // === II. CONFIGURE =======================================================
    // Create and populate the entire simulation configuration from command line.
    ierr = CreateSimulationContext(argc, argv, &simCtx); CHKERRQ(ierr);

    // === III. SETUP ==========================================================
    // Build the simulation environment step-by-step.

    // 1. Setup the multi-grid hierarchy, DMs, and PETSc vectors/solvers.
    ierr = SetupGridAndSolvers(simCtx); CHKERRQ(ierr);

    // 2. Setup the boundary condition handlers for the flow solver.
    ierr = SetupBoundaryConditions(simCtx); CHKERRQ(ierr);

    // 3. Setup domain decomposition info needed for particle migration.
    ierr = SetupDomainRankInfo(simCtx); CHKERRQ(ierr);

    // 4. Setup Initial State of Eulerian Fields
    ierr = InitializeEulerianState(simCtx); CHKERRQ(ierr);

    // 5. If requested, initialize the Lagrangian Particle Swarm.
    if (simCtx->np > 0) {
      ierr = InitializeParticleSwarm(simCtx); CHKERRQ(ierr);
     }

    
    
    /*
     * --- IBM/FSI FEATURE HOOK (Currently Inactive) ---
     * The logic for initializing immersed bodies would be enabled here.
     * The `-immersed` command-line flag would control this block.
     *
    if (simCtx->immersed) {
        ierr = InitializeImmersedBoundary(simCtx); CHKERRQ(ierr);
    }
    */

    // === IV. EXECUTE =========================================================
    // Display a summary banner and run the main time-stepping loop.

    // ierr = DisplayBanner(simCtx); CHKERRQ(ierr);

    if (!simCtx->OnlySetup) {
        // This function will contain the main time loop, orchestrating
        // the flow solve and particle advection steps.
         ierr = AdvanceSimulation(simCtx); CHKERRQ(ierr);
    }// else {
      //   PetscPrintf(simCtx->rank == 0 ? PETSC_COMM_SELF : PETSC_COMM_NULL,
      //              "INFO: SETUP ONLY MODE enabled. Skipping time loop.\n");
// }

    // === V. FINALIZE =========================================================
    // Cleanly destroy all objects and free all memory.

    // ierr = FinalizeSimulation(simCtx); CHKERRQ(ierr);
    ierr = PetscFinalize();

    return ierr;
}
