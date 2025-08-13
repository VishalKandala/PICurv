#include "solvers.h" // The new header we will create

#undef __FUNCT__
#define __FUNCT__ "Flow_Solver"
/**
 * @brief Orchestrates a single time step of the Eulerian fluid solver.
 *
 * This function is the refactored, high-level entry point for advancing the
 * fluid state from time t_n to t_{n+1}. It is a direct adaptation of the
 * legacy Flow_Solver, but it now takes the master SimCtx as its primary
 * argument to eliminate dependencies on global variables.
 *
 * The sequence of operations is:
 * 1.  (Optional) Update turbulence models (RANS/LES) to compute eddy viscosity.
 * 2.  Call the core momentum solver (explicit Runge-Kutta or an implicit scheme)
 *     to get an intermediate velocity field.
 * 3.  Call the pressure-Poisson solver to compute the pressure correction.
 * 4.  Call the projection step to correct the velocity field, ensuring it is
 *     divergence-free.
 * 5.  Perform final state updates, diagnostics, and I/O.
 *
 * @param simCtx The master simulation context, containing all solver settings,
 *               multigrid structures, and data vectors.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode Flow_Solver(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserMG         *usermg = &simCtx->usermg;
    PetscInt       level = usermg->mglevels - 1;
    UserCtx        *user = usermg->mgctx[level].user;
    PetscReal      tm_s, tm_e, tp_e, tpr_e; // Timers for profiling

    PetscFunctionBeginUser;
    LOG_ALLOW(GLOBAL, LOG_INFO, "[Step %d] Entering Flow_Solver orchestrator...\n", simCtx->step);

    /*
    // ========================================================================
    //   SECTION: O-Grid Specific Force Calculations (Legacy Feature)
    // ========================================================================
    // This was a specialized calculation for non-immersed O-grid cases.
    if (simCtx->Ogrid && !simCtx->immersed) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Calculating O-grid forces...\n");
        // Calc_forces_Ogrid(&user[0], simCtx->step, 0);
    }
    */


    /*
    // ========================================================================
    //   SECTION: Turbulence Models (RANS/LES)
    // ========================================================================
    // These models compute the turbulent eddy viscosity (Nu_t) which is then
    // used by the momentum solver in the diffusion term.
    
    if (simCtx->rans) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Updating RANS (k-omega) model...\n");
        for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
            K_Omega_Set_Constant(&user[bi]);
            if (simCtx->step == simCtx->StartStep) {
                LOG_ALLOW(LOCAL, LOG_DEBUG, "  Initializing K-Omega field for block %d.\n", bi);
                K_Omega_IC(&user[bi]);
            }
            // In a full implementation, the K-Omega transport equations would be solved here.
            // Solve_K_Omega(&user[bi]);
        }
    }

    if (simCtx->les) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Updating LES (Smagorinsky) model...\n");
        for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
            // LES models require Cartesian velocity to compute strain rates
            ierr = Contra2Cart(&user[bi]); CHKERRQ(ierr);
            if (simCtx->step % simCtx->dynamic_freq == 0) {
                Compute_Smagorinsky_Constant_1(&user[bi], user[bi].lUcont, user[bi].lUcat);
            }
            Compute_eddy_viscosity_LES(&user[bi]);
        }
    }
    */


    // ========================================================================
    //   SECTION: Momentum Equation Solver
    // ========================================================================
    // This is the core of the time step. It computes an intermediate velocity
    // field by solving the momentum equations.
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Beginning momentum solve (Implicit Flag = %d)...\n", simCtx->implicit);
    ierr = PetscTime(&tm_s); CHKERRQ(ierr);
    
    if (simCtx->implicit > 0) {
        // --- Implicit Path ---
        LOG_ALLOW(GLOBAL, LOG_INFO, "Executing implicit momentum solver...\n");
        // Since IBM is disabled, we pass NULL for ibm and fsi arguments.
        ierr = ImpRK(user, NULL, NULL); CHKERRQ(ierr);
    } else {
        // --- Explicit Path (Default) ---
        LOG_ALLOW(GLOBAL, LOG_INFO, "Executing Iterative implicit momentum solver (Implicit Runge-Kutta)...\n");
        // Since IBM is disabled, we pass NULL for ibm and fsi arguments.
	ierr = RungeKutta(user, NULL, NULL); CHKERRQ(ierr);
    }
    
    ierr = PetscTime(&tm_e); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Momentum solve completed in %.4f seconds.\n", tm_e - tm_s);

    // ========================================================================
    //   !!! DIAGNOSTIC CHECKPOINT 1: AFTER ImpRK !!!
    // ========================================================================
    PetscReal norm_global, norm_local;
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
      VecNorm(user[bi].Ucont, NORM_2, &norm_global);
      VecNorm(user[bi].lUcont, NORM_2, &norm_local);
      LOG_ALLOW(GLOBAL, LOG_INFO, "[DIAGNOSTIC] After ImpRK (Block %d): Global Ucont Norm = %e, Local lUcont Norm = %e\n", bi, norm_global, norm_local);
    }
    // ========================================================================

    // ========================================================================
//   !!! THE FIX: EXPLICIT SYNCHRONIZATION !!!
// ========================================================================
LOG_ALLOW(GLOBAL, LOG_INFO, "SYNC: Forcing update of local Ucont from global Ucont before Poisson solve.\n");
for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
    ierr = DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); CHKERRQ(ierr);
}
// ========================================================================

// ========================================================================
//   !!! DIAGNOSTIC CHECKPOINT 2: AFTER SYNC !!!
// ========================================================================
for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
    VecNorm(user[bi].Ucont, NORM_2, &norm_global);
    VecNorm(user[bi].lUcont, NORM_2, &norm_local);
    LOG_ALLOW(GLOBAL, LOG_INFO, "[DIAGNOSTIC] After Sync (Block %d): Global Ucont Norm = %e, Local lUcont Norm = %e\n", bi, norm_global, norm_local);
}
// ========================================================================
    
// ========================================================================
//   SECTION: Pressure-Poisson Solver
// ========================================================================
// This step enforces the continuity equation (incompressibility) by solving
// for a pressure correction field.
    
 LOG_ALLOW(GLOBAL, LOG_INFO, "Beginning pressure-Poisson solve (Poisson Flag = %d)...\n", simCtx->poisson);
    
 if (simCtx->poisson == 0) {
   ierr = PoissonSolver_MG(usermg); CHKERRQ(ierr);
 } else {
   // Logic for other Poisson solvers (e.g., Hypre) would go here.
   LOG_ALLOW(GLOBAL, LOG_WARNING, "Poisson solver type %d is not currently enabled.\n", simCtx->poisson);
 }
    
 ierr = PetscTime(&tp_e); CHKERRQ(ierr);
 LOG_ALLOW(GLOBAL, LOG_DEBUG, "Pressure-Poisson solve completed in %.4f seconds.\n", tp_e - tm_e);

    // ========================================================================
    //   SECTION: Velocity Correction (Projection)
    // ========================================================================
    // The pressure correction is used to update the pressure field and project
    // the intermediate velocity onto a divergence-free space.
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Applying velocity correction/projection step...\n");
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        ierr = UpdatePressure(&user[bi]); CHKERRQ(ierr);
        ierr = Projection(&user[bi]); CHKERRQ(ierr);
        
        // Ensure local ghost cells for the final pressure field are correct
        ierr = DMGlobalToLocalBegin(user[bi].da, user[bi].P, INSERT_VALUES, user[bi].lP); CHKERRQ(ierr);
        ierr = DMGlobalToLocalEnd(user[bi].da, user[bi].P, INSERT_VALUES, user[bi].lP); CHKERRQ(ierr);
    }
    
    ierr = PetscTime(&tpr_e); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Velocity correction completed in %.4f seconds.\n", tpr_e - tp_e);


    // ========================================================================
    //   SECTION: Final Diagnostics and I/O
    // ========================================================================
    
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        LOG_ALLOW(LOCAL, LOG_INFO, "Finalizing state for block %d...\n", bi);
        
        // --- Perform Divergence Check ---
        // This is a diagnostic to verify the quality of the velocity correction.
	//  ierr = Divergence(&user[bi]); CHKERRQ(ierr);
        
        /*
        // --- Immersed Boundary Interpolation (Post-Correction) ---
        // This step would update the velocity values AT the IB nodes to match the
        // newly corrected fluid field. Important for the next time step.
        if (simCtx->immersed) {
            for (PetscInt ibi = 0; ibi < simCtx->NumberOfBodies; ibi++) {
                ibm_interpolation_advanced(&user[bi], &simCtx->ibm[ibi], ibi, 1);
            }
        }
        */

        // --- Averaging and Statistics (if enabled) ---
        /*
        if (simCtx->averaging) {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Performing statistical averaging for block %d.\n", bi);
            Do_averaging(&user[bi]);
        }
        */

        // --- Per-Step I/O ---
        // The I/O is now handled in the main AdvanceSimulation loop to give
        // the user more control. The legacy call is commented out here.
        // if (simCtx->step % simCtx->tiout == 0) {
        //     Ucont_P_Binary_Output(&user[bi]);
        // }
    }
    
    // --- Profiling Output ---
    // (This can be kept as is, but it's good practice to wrap it in a LOG_PROFILE check)
    LOG_PROFILE_MSG(GLOBAL, "TIMING [Step %d]: Momentum=%.4fs, Poisson=%.4fs, Projection=%.4fs, TOTAL=%.4fs\n",
                    simCtx->step, tm_e - tm_s, tp_e - tm_e, tpr_e - tp_e, tpr_e - tm_s);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Flow_Solver orchestrator finished for step %d.\n", simCtx->step);
    PetscFunctionReturn(0);
}
