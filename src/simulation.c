/**
 * @file simulation.c  // code for simulation loop 
 * @brief Test program for DMSwarm interpolation using the fdf-curvIB method.
 * Provides the setup to start any simulation with DMSwarm and DMDAs.
 **/

#include "simulation.h"

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
PetscErrorCode UpdateSolverHistoryVectors(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx         *simCtx = user->simCtx; // Access global settings if needed

    PetscFunctionBeginUser;
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d, Block %d: Updating solver history vectors.\n",
              simCtx->rank, user->_this);

    // --- Primary Contravariant Velocity History ---
    // The order is critical here.
    // 1. First, move the n-1 state (Ucont_o) to the n-2 slot (Ucont_rm1).
    ierr = VecCopy(user->Ucont_o, user->Ucont_rm1); CHKERRQ(ierr);
    // 2. Then, move the new n state (Ucont) to the n-1 slot (Ucont_o).
    ierr = VecCopy(user->Ucont, user->Ucont_o); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL,LOG_DEBUG, "Rank %d, Block %d, Ucont history updated.\n",simCtx->rank,user->_this); 
    
    // --- Update History for Other Fields ---
    // These are typically only needed at the n-1 state.
    ierr = VecCopy(user->Ucat, user->Ucat_o); CHKERRQ(ierr);
    ierr = VecCopy(user->P, user->P_o); CHKERRQ(ierr);
    LOG_ALLOW(LOCAL,LOG_DEBUG, "Rank %d, Block %d, Ucat & P  history updated.\n",simCtx->rank,user->_this); 
    
    if (simCtx->immersed) {
        ierr = VecCopy(user->Nvert, user->Nvert_o); CHKERRQ(ierr);
    }

    // --- Update History for Turbulence Models (if active) ---
    if (simCtx->rans) {
       ierr = VecCopy(user->K_Omega, user->K_Omega_o); CHKERRQ(ierr);
    }
    
    // --- Synchronize Local Ghost Regions for the new history vectors ---
    // This is essential so that stencils in the next time step's calculations
    // have correct values from neighboring processes.
    ierr = DMGlobalToLocalBegin(user->fda, user->Ucont_o, INSERT_VALUES, user->lUcont_o); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->fda, user->Ucont_o, INSERT_VALUES, user->lUcont_o); CHKERRQ(ierr);

    ierr = DMGlobalToLocalBegin(user->fda, user->Ucont_rm1, INSERT_VALUES, user->lUcont_rm1); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->fda, user->Ucont_rm1, INSERT_VALUES, user->lUcont_rm1); CHKERRQ(ierr);
    
    if (simCtx->immersed) {
        ierr = DMGlobalToLocalBegin(user->da, user->Nvert_o, INSERT_VALUES, user->lNvert_o); CHKERRQ(ierr);
        ierr = DMGlobalToLocalEnd(user->da, user->Nvert_o, INSERT_VALUES, user->lNvert_o); CHKERRQ(ierr);
    }
    
    if (simCtx->rans) {
       ierr = DMGlobalToLocalBegin(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o); CHKERRQ(ierr);
       ierr = DMGlobalToLocalEnd(user->fda2, user->K_Omega_o, INSERT_VALUES, user->lK_Omega_o); CHKERRQ(ierr);
    }
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PerformInitialSetup"
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
PetscErrorCode PerformInitializedParticleSetup(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    // --- Get pointers from SimCtx instead of passing them as arguments ---
    UserCtx     *user     = simCtx->usermg.mgctx[simCtx->usermg.mglevels-1].user;
    BoundingBox *bboxlist = simCtx->bboxlist;

    PetscFunctionBeginUser;
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Performing initial particle setup procedures.\n", simCtx->ti, simCtx->step);

    // --- 1. Initial Particle Settlement (Location and Migration) ---
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Initial Settlement: Locating and migrating all particles...\n", simCtx->ti, simCtx->step);
    ierr = LocateAllParticlesInGrid(user, bboxlist); CHKERRQ(ierr);

    if(get_log_level() == LOG_DEBUG && is_function_allowed(__FUNCT__)==true){
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Particle field states after Initial settlement...\n", simCtx->ti, simCtx->step);
        ierr = LOG_PARTICLE_FIELDS(user,simCtx->LoggingFrequency); CHKERRQ(ierr);
    }

    // --- 2. Re-initialize Particles on Inlet Surface (if applicable) ---
    if ((simCtx->ParticleInitialization == 0 || simCtx->ParticleInitialization == 3) && user->inletFaceDefined) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Re-initializing particles on inlet surface...\n", simCtx->ti, simCtx->step);
        ierr = ReinitializeParticlesOnInletSurface(user, simCtx->ti, simCtx->step); CHKERRQ(ierr);

        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Resetting statuses for post-reinitialization settlement.\n", simCtx->ti, simCtx->step);
        ierr = ResetAllParticleStatuses(user); CHKERRQ(ierr);
    
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Post-Reinitialization Settlement...\n", simCtx->ti, simCtx->step);
        ierr = LocateAllParticlesInGrid(user, bboxlist); CHKERRQ(ierr);

    }
    
    // --- 3. Finalize State for t=0 ---
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Interpolating initial fields to settled particles.\n", simCtx->ti, simCtx->step);
    ierr = InterpolateAllFieldsToSwarm(user); CHKERRQ(ierr);
    //ierr = ScatterAllParticleFieldsToEulerFields(user); CHKERRQ(ierr);

    // --- 4. Initial History and Output ---
    // Update solver history vectors with the t=0 state before the first real step
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        ierr = UpdateSolverHistoryVectors(&user[bi]); CHKERRQ(ierr);
    }

    if (simCtx->OutputFreq > 0 || (simCtx->StepsToRun == 0 && simCtx->StartStep == 0)) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Writing initial simulation data.\n", simCtx->ti, simCtx->step);
        
        // --- Particle Output (assumes functions operate on the master user context) ---
        if(get_log_level() == LOG_DEBUG){
            LOG(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Initial Particle field.\n", simCtx->ti, simCtx->step);
            ierr = LOG_PARTICLE_FIELDS(user,simCtx->LoggingFrequency); CHKERRQ(ierr);
        }
        ierr = WriteAllSwarmFields(user); CHKERRQ(ierr);
        
        // --- Eulerian Field Output (MUST loop over all blocks) --- // <<< CHANGED/FIXED
        for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
            ierr = WriteSimulationFields(&user[bi]); CHKERRQ(ierr);
        }
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "--- Initial setup complete. Ready for time marching. ---\n");
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PerformLoadedParticleSetup"
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
PetscErrorCode PerformLoadedParticleSetup(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    UserCtx     *user     = simCtx->usermg.mgctx[simCtx->usermg.mglevels-1].user;

    // 0. Reset all particle statuses to ensure proper location checks.
    ierr = ResetAllParticleStatuses(user); CHKERRQ(ierr);

    // 1. Rebuild grid-to-particle links based on loaded coordinates.
    ierr = LocateAllParticlesInGrid(user, simCtx->bboxlist); CHKERRQ(ierr);

    if(get_log_level() == LOG_DEBUG){
        LOG(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Particle field states after locating loaded particles...\n", simCtx->ti, simCtx->step);
        ierr = LOG_PARTICLE_FIELDS(user,simCtx->LoggingFrequency); CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Interpolating initial fields to settled particles.\n", simCtx->ti, simCtx->step);

    // 2. Ensure particles have velocity from the authoritative loaded grid for consistency.
    ierr = InterpolateAllFieldsToSwarm(user); CHKERRQ(ierr);

    // 3. Update Eulerian source terms from the loaded particle data.
    //ierr = ScatterAllParticleFieldsToEulerFields(user); CHKERRQ(ierr);

    // --- 4. Initial History and Output ---
    // Update solver history vectors with the t=0 state before the first real step
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        ierr = UpdateSolverHistoryVectors(&user[bi]); CHKERRQ(ierr);
    }

    if (simCtx->OutputFreq > 0 || (simCtx->StepsToRun == 0 && simCtx->StartStep == 0)) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Writing initial simulation data.\n", simCtx->ti, simCtx->step);
        
        // --- Particle Output (assumes functions operate on the master user context) ---
        if(get_log_level() >=LOG_INFO) ierr = LOG_PARTICLE_FIELDS(user, simCtx->LoggingFrequency); CHKERRQ(ierr);
        ierr = WriteAllSwarmFields(user); CHKERRQ(ierr);
        
        // --- Eulerian Field Output (MUST loop over all blocks) --- // <<< CHANGED/FIXED
        for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
            ierr = WriteSimulationFields(&user[bi]); CHKERRQ(ierr);
        }
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "--- Initial setup complete. Ready for time marching. ---\n");
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FinalizeRestartState"
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
PetscErrorCode FinalizeRestartState(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserCtx       *user = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;

    PetscFunctionBeginUser;

    LOG_ALLOW(GLOBAL, LOG_INFO, "--- Finalizing RESTART from state (step=%d, t=%.4f) ---\n", simCtx->StartStep, simCtx->ti);

    // This function only needs to handle the particle finalization logic.
    // The Eulerian state is assumed to be fully loaded and consistent at this point.
    if (simCtx->np > 0) {

        // Use the particle restart mode to decide the workflow.
        if (strcmp(simCtx->particleRestartMode, "load") == 0) {
            // PARTICLES WERE LOADED: The state is complete, but we must verify
            // the loaded CellIDs and build the in-memory grid-to-particle links.
            LOG_ALLOW(GLOBAL, LOG_INFO, "Particle Mode 'load': Verifying particle locations and building grid links...\n");
            ierr = PerformLoadedParticleSetup(simCtx); CHKERRQ(ierr);

        } else { // Mode must be "init"
            // PARTICLES WERE RE-INITIALIZED: They need to be fully settled and coupled
            // to the surrounding (restarted) fluid state.
            LOG_ALLOW(GLOBAL, LOG_INFO, "Particle Mode 'init': Running full initial setup for new particles in restarted flow.\n");
            ierr = PerformInitializedParticleSetup(simCtx); CHKERRQ(ierr);
        }
    } else {
        LOG_ALLOW(GLOBAL, LOG_INFO, "No particles in simulation, restart finalization is complete.\n");

        // Write the initial eulerian fields (this is done in PerformInitialSetup if particles exist.)
        for(PetscInt bi = 0; bi < simCtx->block_number; bi ++){
            ierr = WriteSimulationFields(&user[bi]); CHKERRQ(ierr);
        }
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "--- Restart state successfully finalized. --\n");

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "AdvanceSimulation"
/**
 * @brief Executes the main time-marching loop for the coupled Euler-Lagrange simulation.
 *
 * This function orchestrates the advancement of the simulation from the configured
 * StartStep to the final step. It does NOT perform the initial t=0 setup, as that
 * is handled by InitializeEulerianState and PerformInitialSetup in main().
 *
 * For each timestep, it performs the following sequence:
 *  1.  **Pre-Solver Actions:** Updates time-dependent boundary conditions (e.g., fluxin)
 *      and resets particle states for the new step.
 *  2.  **Eulerian Solve:** Calls the refactored legacy FlowSolver to advance the
 *      entire fluid field by one time step.
 *  3.  **Lagrangian Update:** Executes the full particle workflow: advection based
 *      on the previous step's velocity, followed by settling (location/migration)
 *      in the new grid, and finally interpolation of the new fluid velocity.
 *  4.  **Two-Way Coupling:** Scatters particle data back to the grid to act as source
 *      terms for the subsequent time step.
 *  5.  **History & I/O:** Updates the solver's history vectors and writes output files
 *      at the specified frequency.
 *
 * @param user       Array of UserCtx structures for the finest grid level.
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode AdvanceSimulation(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    // Get the master context from the first block. All blocks share it.
    UserCtx *user = simCtx->usermg.mgctx[simCtx->usermg.mglevels-1].user;

    // Retrieve control parameters from SimCtx for clarity
    const PetscInt  StartStep = simCtx->StartStep;
    const PetscInt  StepsToRun = simCtx->StepsToRun;
    const PetscInt  OutputFreq = simCtx->OutputFreq;
    const PetscReal dt = simCtx->dt;
    
    // Variables for particle removal statistics
    PetscInt removed_local_ob, removed_global_ob;
    PetscInt removed_local_lost, removed_global_lost;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting main time-marching loop: %d steps from step %d (t=%.4f), dt=%.4f\n",
              StepsToRun, StartStep, simCtx->StartTime, dt);

    // --- Main Time-Marching Loop ---
    for (PetscInt step = StartStep; step < StartStep + StepsToRun; step++) {
        
        // =================================================================
        //     1. PRE-STEP SETUP
        // =================================================================
        
        // Update simulation time and step counters in the master context
        simCtx->step = step + 1;
        simCtx->ti   += simCtx->dt; //simCtx->StartTime + step * simCtx->dt;

        LOG_ALLOW(GLOBAL, LOG_INFO, "--- Advancing Step %d (To t=%.4f) ---\n", simCtx->step, simCtx->ti);
        

        // For particles, reset their status to prepare for the new advection/location cycle
        if (simCtx->np > 0) {
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "Resetting all particle statuses to NEEDS_LOCATION.\n");
            ierr = ResetAllParticleStatuses(user); CHKERRQ(ierr);
        }

        // =================================================================
        //     2. EULERIAN SOLVER STEP
        // =================================================================
        if(get_log_level() == LOG_VERBOSE && is_function_allowed(__FUNCT__)==true){
            ierr = LOG_FIELD_ANATOMY(&user[0],"Coordinates","PreFlowSolver"); CHKERRQ(ierr);
            ierr = LOG_FIELD_ANATOMY(&user[0],"Csi","PreFlowSolver"); CHKERRQ(ierr);
            ierr = LOG_FIELD_ANATOMY(&user[0],"Eta","PreFlowSolver"); CHKERRQ(ierr);
            ierr = LOG_FIELD_ANATOMY(&user[0],"Zet","PreFlowSolver"); CHKERRQ(ierr);
        }
        LOG_ALLOW(GLOBAL, LOG_INFO, "Updating Eulerian Field ...\n");
        if(strcmp(simCtx->eulerianSource,"load")==0){
            //LOAD mode: Read pre-computed fields for the current step.
            LOG_ALLOW(GLOBAL,LOG_INFO,"Eulerian Source 'load': Reading fields (t=%.4f,step=%d)...\n",simCtx->ti,simCtx->step);
            for(PetscInt bi = 0; bi < simCtx->block_number;bi++){
                ierr = ReadSimulationFields(&user[bi],simCtx->step); CHKERRQ(ierr);
            }
        }else if(strcmp(simCtx->eulerianSource,"solve")==0){
            // SOLVE mode:Call the refactored, high-level legacy solver. This single function
            // advances the entire multi-block fluid field from t_n to t_{n+1}.
            LOG_ALLOW(GLOBAL,LOG_INFO,"Eulerian Source 'solve'. Updating Eulerian field via Solver...\n");
            ierr = FlowSolver(simCtx); CHKERRQ(ierr);
        }
        LOG_ALLOW(GLOBAL, LOG_INFO, "Eulerian Field Updated ...\n");
        if(get_log_level() == LOG_VERBOSE && is_function_allowed(__FUNCT__)==true){
            LOG_ALLOW(GLOBAL, LOG_VERBOSE, "Post FlowSolver field states:\n");
            ierr = LOG_FIELD_ANATOMY(&user[0],"Ucat","PostFlowSolver"); CHKERRQ(ierr);
            ierr = LOG_FIELD_ANATOMY(&user[0],"P","PostFlowSolver"); CHKERRQ(ierr);
            ierr = LOG_FIELD_ANATOMY(&user[0],"Ucont","PostFlowSolver"); CHKERRQ(ierr);
        }


        // =================================================================
        //     3. LAGRANGIAN PARTICLE STEP
        // =================================================================
        
        if (simCtx->np > 0) {
            LOG_ALLOW(GLOBAL, LOG_INFO, "Updating Lagrangian particle system...\n");

            // a. Advect particles using the velocity interpolated from the *previous* step.
            //    P(t_{n+1}) = P(t_n) + V_p(t_n) * dt
            ierr = UpdateAllParticlePositions(user); CHKERRQ(ierr);

            // b. Settle all particles: find their new host cells and migrate them across ranks.
            ierr = LocateAllParticlesInGrid(user, simCtx->bboxlist); CHKERRQ(ierr);
 
            // c. Remove any particles that are now lost or out of the global domain.
            ierr = CheckAndRemoveLostParticles(user, &removed_local_lost, &removed_global_lost); CHKERRQ(ierr);
            ierr = CheckAndRemoveOutOfBoundsParticles(user, &removed_local_ob, &removed_global_ob, simCtx->bboxlist); CHKERRQ(ierr);
            if (removed_global_lost + removed_global_ob > 0) {
                LOG_ALLOW(GLOBAL, LOG_INFO, "Removed %d particles globally this step.\n", removed_global_lost + removed_global_ob);
            }

            // d. Interpolate the NEW fluid velocity (just computed by FlowSolver) onto the
            //    particles' new positions. This gives them V_p(t_{n+1}) for the *next* advection step.
            ierr = InterpolateAllFieldsToSwarm(user); CHKERRQ(ierr);

            // e. Update the Particle Fields (e.g., temperature, concentration) if applicable.
            //    This can be extended to include reactions, growth, etc.
            ierr = UpdateAllParticleFields(user); CHKERRQ(ierr);
            
            // f. (For Two-Way Coupling) Scatter particle data back to the grid to act as a source term.
           // ierr = CalculateParticleCountPerCell(user); CHKERRQ(ierr);
           // ierr = ScatterAllParticleFieldsToEulerFields(user); CHKERRQ(ierr);

            if(get_log_level() == LOG_VERBOSE && is_function_allowed(__FUNCT__)==true){
                LOG_ALLOW(GLOBAL, LOG_VERBOSE, "Post Lagrangian update field states:\n");
                ierr = LOG_FIELD_MIN_MAX(&user[0],"Psi"); CHKERRQ(ierr);
            }
        }

        // =================================================================
        //     4. UPDATE HISTORY & I/O
        // =================================================================
        
        // Copy the newly computed fields (Ucont, P, etc.) to the history vectors
        // (_o, _rm1) to prepare for the next time step.
        for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
            ierr = UpdateSolverHistoryVectors(&user[bi]); CHKERRQ(ierr);
        }
        
        //ierr = LOG_UCAT_ANATOMY(&user[0],"Final"); CHKERRQ(ierr);
        
        // Handle periodic file output
        if (OutputFreq > 0 && (step + 1) % OutputFreq == 0) {
            LOG_ALLOW(GLOBAL, LOG_INFO, "Writing output for step %d.\n",simCtx->step);
            for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
                ierr = WriteSimulationFields(&user[bi]); CHKERRQ(ierr);
            }
            if (simCtx->np > 0) {
                ierr = WriteAllSwarmFields(user); CHKERRQ(ierr);
                if(get_log_level() >= LOG_INFO){
                LOG(GLOBAL, LOG_INFO, "Particle states at step %d:\n", simCtx->step);
                LOG_PARTICLE_FIELDS(user,simCtx->LoggingFrequency);
                }
            }
        }

        ProfilingLogTimestepSummary(simCtx->step);

        // Update Progress Bar
        if(simCtx->rank == 0) {
            PrintProgressBar(step,StartStep,StepsToRun,simCtx->ti);
            if(get_log_level()>=LOG_WARNING) PetscPrintf(PETSC_COMM_SELF,"\n");
        }
    } // --- End of Time-Marching Loop ---

    // After the loop, print the 100% complete bar on rank 0 and add a newline
    // to ensure subsequent terminal output starts on a fresh line.
    if (simCtx->rank == 0 && StepsToRun > 0) {
        PrintProgressBar(StartStep + StepsToRun - 1, StartStep, StepsToRun, simCtx->ti);
        PetscPrintf(PETSC_COMM_SELF, "\n");
        fflush(stdout);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Time marching completed. Final time t=%.4f.\n", simCtx->ti);
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
