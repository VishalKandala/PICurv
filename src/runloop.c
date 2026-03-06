/**
 * @file runloop.c
 * @brief Test program for DMSwarm interpolation using the fdf-curvIB method.
 * Provides the setup to start any simulation with DMSwarm and DMDAs.
 **/

#include "runloop.h"

/**
 * @brief Internal helper implementation: `UpdateSolverHistoryVectors()`.
 * @details Local to this translation unit.
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
 * @brief Internal helper implementation: `PerformInitializedParticleSetup()`.
 * @details Local to this translation unit.
 */
PetscErrorCode PerformInitializedParticleSetup(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    // --- Get pointers from SimCtx instead of passing them as arguments ---
    UserCtx     *user     = simCtx->usermg.mgctx[simCtx->usermg.mglevels-1].user;
    BoundingBox *bboxlist = simCtx->bboxlist;

    PetscFunctionBeginUser;
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Performing initial particle setup procedures.\n", simCtx->ti, simCtx->step);

    // --- 0. Loop over all blocks to compute Eulerian diffusivity.
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        ierr = ComputeEulerianDiffusivity(&user[bi]); CHKERRQ(ierr);
        ierr = ComputeEulerianDiffusivityGradient(&user[bi]); CHKERRQ(ierr);
    }

    // --- 1. Initial Particle Settlement (Location and Migration) ---
    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Initial Settlement: Locating and migrating all particles...\n", simCtx->ti, simCtx->step);
    ierr = LocateAllParticlesInGrid(user, bboxlist); CHKERRQ(ierr);

    if(get_log_level() >= LOG_DEBUG && is_function_allowed(__FUNCT__)==true){
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Particle field states after Initial settlement...\n", simCtx->ti, simCtx->step);
        ierr = LOG_PARTICLE_FIELDS(user,simCtx->LoggingFrequency); CHKERRQ(ierr);
    }

    // --- 2. Re-initialize Particles on Inlet Surface (if applicable) ---
    if ((simCtx->ParticleInitialization == PARTICLE_INIT_SURFACE_RANDOM ||
         simCtx->ParticleInitialization == PARTICLE_INIT_SURFACE_EDGES) && user->inletFaceDefined) {
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

    if (simCtx->tiout > 0 || (simCtx->StepsToRun == 0 && simCtx->StartStep == 0)) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Writing initial simulation data.\n", simCtx->ti, simCtx->step);
        ierr = WriteAllSwarmFields(user); CHKERRQ(ierr);
        
        // --- Eulerian Field Output (MUST loop over all blocks) --- // <<< CHANGED/FIXED
        for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
            ierr = WriteSimulationFields(&user[bi]); CHKERRQ(ierr);
        }
    }

    if (IsParticleConsoleSnapshotEnabled(simCtx)) {
        ierr = EmitParticleConsoleSnapshot(user, simCtx, simCtx->step); CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "--- Initial setup complete. Ready for time marching. ---\n");
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PerformLoadedParticleSetup"
/**
 * @brief Internal helper implementation: `PerformLoadedParticleSetup()`.
 * @details Local to this translation unit.
 */
PetscErrorCode PerformLoadedParticleSetup(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    UserCtx     *user     = simCtx->usermg.mgctx[simCtx->usermg.mglevels-1].user;

    // --- 0. Re-compute Eulerian Diffusivity from loaded fields.
    LOG_ALLOW(GLOBAL, LOG_INFO, "Re-computing Eulerian Diffusivity from loaded fields...\n");
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        ierr = ComputeEulerianDiffusivity(&user[bi]); CHKERRQ(ierr);
        ierr = ComputeEulerianDiffusivityGradient(&user[bi]); CHKERRQ(ierr);
    }

    // 0.1 This moves particles to their correct ranks immediately using the loaded Cell ID.
    LOG_ALLOW(GLOBAL, LOG_INFO, "Performing fast restart migration using preloaded Cell IDs...\n");
    ierr = MigrateRestartParticlesUsingCellID(user); CHKERRQ(ierr);

    // 1. To catch any edge cases (particles with invalid CellIDs or newcomers).
    // Because we kept the statuses, this function will now SKIP all the particles 
    // that are already on the correct rank,
    ierr = LocateAllParticlesInGrid(user, simCtx->bboxlist); CHKERRQ(ierr);

    if(get_log_level() == LOG_DEBUG){
        LOG(GLOBAL, LOG_DEBUG, "[T=%.4f, Step=%d] Particle field states after locating loaded particles...\n", simCtx->ti, simCtx->step);
        ierr = LOG_PARTICLE_FIELDS(user,simCtx->LoggingFrequency); CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Interpolating initial fields to settled particles.\n", simCtx->ti, simCtx->step);

    // 2. Ensure particles have velocity from the authoritative loaded grid for consistency.
    ierr = InterpolateAllFieldsToSwarm(user); CHKERRQ(ierr);

    // 3. Update Eulerian source terms from the loaded particle data.
    ierr = ScatterAllParticleFieldsToEulerFields(user); CHKERRQ(ierr);

    // --- 4. Initial History and Output ---
    // Update solver history vectors with the t=0 state before the first real step
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        ierr = UpdateSolverHistoryVectors(&user[bi]); CHKERRQ(ierr);
    }

    if (simCtx->tiout > 0 || (simCtx->StepsToRun == 0 && simCtx->StartStep == 0)) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "[T=%.4f, Step=%d] Writing initial simulation data.\n", simCtx->ti, simCtx->step);
        ierr = WriteAllSwarmFields(user); CHKERRQ(ierr);
        
        // --- Eulerian Field Output (MUST loop over all blocks) --- // <<< CHANGED/FIXED
        for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
            ierr = WriteSimulationFields(&user[bi]); CHKERRQ(ierr);
        }
    }

    if (IsParticleConsoleSnapshotEnabled(simCtx)) {
        ierr = EmitParticleConsoleSnapshot(user, simCtx, simCtx->step); CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "--- Initial setup complete. Ready for time marching. ---\n");
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FinalizeRestartState"
/**
 * @brief Internal helper implementation: `FinalizeRestartState()`.
 * @details Local to this translation unit.
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
 * @brief Internal helper implementation: `AdvanceSimulation()`.
 * @details Local to this translation unit.
 */
PetscErrorCode AdvanceSimulation(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    // Get the master context from the first block. All blocks share it.
    UserCtx *user = simCtx->usermg.mgctx[simCtx->usermg.mglevels-1].user;

    // Retrieve control parameters from SimCtx for clarity
    const PetscInt  StartStep = simCtx->StartStep;
    const PetscInt  StepsToRun = simCtx->StepsToRun;
    const PetscReal dt = simCtx->dt;
    
    // Variables for particle removal statistics
    //PetscInt removed_local_ob, removed_global_ob;
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
            ierr = LOG_FIELD_ANATOMY(&user[0],"Center-Coordinates","PreFlowSolver"); CHKERRQ(ierr);
            ierr = LOG_FIELD_ANATOMY(&user[0],"X-Face-Centers","PreFlowSolver"); CHKERRQ(ierr);
            ierr = LOG_FIELD_ANATOMY(&user[0],"Y-Face-Centers","PreFlowSolver"); CHKERRQ(ierr);
            ierr = LOG_FIELD_ANATOMY(&user[0],"Z-Face-Centers","PreFlowSolver"); CHKERRQ(ierr);
            ierr = LOG_FIELD_ANATOMY(&user[0],"Ucat","PreFlowSolver"); CHKERRQ(ierr);
        }
        LOG_ALLOW(GLOBAL, LOG_INFO, "Updating Eulerian Field ...\n");
        if(strcmp(simCtx->eulerianSource,"load")==0){
            //LOAD mode: Read pre-computed fields for the current step.
            LOG_ALLOW(GLOBAL,LOG_INFO,"Eulerian Source 'load': Reading fields (t=%.4f,step=%d)...\n",simCtx->ti,simCtx->step);
            for(PetscInt bi = 0; bi < simCtx->block_number;bi++){
                ierr = ReadSimulationFields(&user[bi],simCtx->step); CHKERRQ(ierr);
            }
        }else if(strcmp(simCtx->eulerianSource,"analytical")==0){
            // ANALYTICAL mode:Call the Analytical Solution Prescription Engine to enable a variety of analytical functions
            LOG_ALLOW(GLOBAL,LOG_INFO,"Eulerian Source 'analytical'. Updating Eulerian field via the Analytical Solution Engine ...\n");
            ierr = AnalyticalSolutionEngine(simCtx); CHKERRQ(ierr);
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

            // a. Update Eulerian Transport Properties:
            // Optimization: Only recalculate if turbulence is active (Nu_t changes).
            // For Laminar flow, the value calculated at Setup is constant.
            if (simCtx->les || simCtx->rans) {
                for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
                    ierr = ComputeEulerianDiffusivity(&user[bi]); CHKERRQ(ierr);
                    ierr = ComputeEulerianDiffusivityGradient(&user[bi]); CHKERRQ(ierr);
                }
            }

            // a.1 (Optional) Log Eulerian Diffusivity min/max and anatomy for debugging.
            if(get_log_level() == LOG_VERBOSE && is_function_allowed(__FUNCT__)==true){
                LOG_ALLOW(GLOBAL, LOG_VERBOSE, "Updated Diffusivity Min/Max:\n");
                ierr = LOG_FIELD_MIN_MAX(&user[0],"Diffusivity"); CHKERRQ(ierr);
                ierr = LOG_FIELD_MIN_MAX(&user[0],"DiffusivityGradient"); CHKERRQ(ierr);
                //LOG_ALLOW(GLOBAL, LOG_VERBOSE, "Updated Diffusivity Anatomy:\n");
                ierr = LOG_FIELD_ANATOMY(&user[0],"Diffusivity","PostDiffusivityUpdate"); CHKERRQ(ierr);
            }
            // b. Advect particles using the velocity interpolated from the *previous* step.
            //    P(t_{n+1}) = P(t_n) + V_p(t_n) * dt
            ierr = UpdateAllParticlePositions(user); CHKERRQ(ierr);

            // c. Settle all particles: find their new host cells and migrate them across ranks.
            ierr = LocateAllParticlesInGrid(user, simCtx->bboxlist); CHKERRQ(ierr);
 
            // d. Remove any particles that are now lost or out of the global domain.
            ierr = CheckAndRemoveLostParticles(user, &removed_local_lost, &removed_global_lost); CHKERRQ(ierr);
            //ierr = CheckAndRemoveOutOfBoundsParticles(user, &removed_local_ob, &removed_global_ob, simCtx->bboxlist); CHKERRQ(ierr);
            if (removed_global_lost> 0) { // if(removed_global_lost + removed_global_ob > 0){
                LOG_ALLOW(GLOBAL, LOG_INFO, "Removed %d particles globally this step.\n", removed_global_lost); // removed_global_lost + removed_global_ob;
                simCtx->particlesLostLastStep = removed_global_lost;
            }

            // e. Interpolate the NEW fluid velocity (just computed by FlowSolver) onto the
            //    particles' new positions. This gives them V_p(t_{n+1}) for the *next* advection step.
            ierr = InterpolateAllFieldsToSwarm(user); CHKERRQ(ierr);

            // f. Update the Particle Fields (e.g., temperature, concentration) if applicable.
            //    This can be extended to include reactions, growth, etc.
            ierr = UpdateAllParticleFields(user); CHKERRQ(ierr);
            
            // g. (For Two-Way Coupling) Scatter particle data back to the grid to act as a source term.
            ierr = CalculateParticleCountPerCell(user); CHKERRQ(ierr);
            ierr = ScatterAllParticleFieldsToEulerFields(user); CHKERRQ(ierr);
            
            // h. (Optional) Calculate advanced particle metrics for logging/debugging.
            ierr = CalculateAdvancedParticleMetrics(user); CHKERRQ(ierr);

            ierr = LOG_PARTICLE_METRICS(user, "Timestep Metrics"); CHKERRQ(ierr);


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
        if (ShouldWriteDataOutput(simCtx, simCtx->step)) {
            LOG_ALLOW(GLOBAL, LOG_INFO, "Writing output for step %d.\n",simCtx->step);
            for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
                ierr = WriteSimulationFields(&user[bi]); CHKERRQ(ierr);
            }
            if (simCtx->np > 0) {
                ierr = WriteAllSwarmFields(user); CHKERRQ(ierr);
                if (get_log_level() >= LOG_INFO && strcmp(simCtx->eulerianSource,"analytical") == 0) {
                    LOG_INTERPOLATION_ERROR(user);
                }
            }
        }

        if (ShouldEmitPeriodicParticleConsoleSnapshot(simCtx, simCtx->step)) {
            ierr = EmitParticleConsoleSnapshot(user, simCtx, simCtx->step); CHKERRQ(ierr);
        }

        ProfilingLogTimestepSummary(simCtx, simCtx->step);

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
