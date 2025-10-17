#include "solvers.h" // The new header we will create

#undef __FUNCT__
#define __FUNCT__ "FlowSolver"
/**
 * @brief Orchestrates a single time step of the Eulerian fluid solver.
 *
 * This function is the refactored, high-level entry point for advancing the
 * fluid state from time t_n to t_{n+1}. It is a direct adaptation of the
 * legacy FlowSolver, but it now takes the master SimCtx as its primary
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
PetscErrorCode FlowSolver(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserMG         *usermg = &simCtx->usermg;
    PetscInt       level = usermg->mglevels - 1;
    UserCtx        *user = usermg->mgctx[level].user;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(GLOBAL, LOG_INFO, "[Step %d] Entering orchestrator...\n", simCtx->step);

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


    
    // ========================================================================
    //   SECTION: Turbulence Models (RANS/LES)
    // ========================================================================
    // These models compute the turbulent eddy viscosity (Nu_t) which is then
    // used by the momentum solver in the diffusion term.

    /*
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
    */

    if (simCtx->les) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Updating LES (Smagorinsky) model...\n");
        for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
            // LES models require Cartesian velocity to compute strain rates
            UpdateLocalGhosts(&user[bi], "Ucont");
            ierr = Contra2Cart(&user[bi]); CHKERRQ(ierr);
            UpdateLocalGhosts(&user[bi], "Ucat");
            if (simCtx->step % simCtx->dynamic_freq == 0) {
                LOG_ALLOW(LOCAL, LOG_DEBUG, "  Computing Smagorinsky constant for block %d.\n", bi);
                ComputeSmagorinskyConstant(&user[bi]);
            }
          //  LOG_ALLOW(LOCAL, LOG_DEBUG, "  Computing eddy viscosity for block %d.\n", bi);
          ComputeEddyViscosityLES(&user[bi]);
        }
    }
    

    // ========================================================================
    //   SECTION: Momentum Equation Solver
    // ========================================================================
    // This is the core of the time step. It computes an intermediate velocity
    // field by solving the momentum equations.
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Beginning momentum solve (Implicit Flag = %d)...\n", simCtx->implicit);
    
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
    
    // ========================================================================
    //   SECTION: Velocity Correction (Projection)
    // ========================================================================
    // The pressure correction is used to update the pressure field and project
    // the intermediate velocity onto a divergence-free space.
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Applying velocity correction/projection step...\n");
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        ierr = UpdatePressure(&user[bi]); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_INFO," Pressure Updated for Block %d.\n",bi);
	
        ierr = Projection(&user[bi]); CHKERRQ(ierr);
	
        LOG_ALLOW(GLOBAL,LOG_INFO," Velocity corrected for Block %d.\n",bi);

	// Ensure local ghost cells for the final pressure field are correct
	ierr = UpdateLocalGhosts(&user[bi],"P");
    }
    
    // ========================================================================
    
    // ========================================================================
    //   SECTION: Final Diagnostics and I/O
    // ========================================================================
    
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Finalizing state & Diagnostics for block %d...\n", bi);
        
        // --- Perform Divergence Check ---
        // This is a diagnostic to verify the quality of the velocity correction.
	ierr = ComputeDivergence(&user[bi]); CHKERRQ(ierr);

	// -- Log Continuity metrics ----
	ierr = LOG_CONTINUITY_METRICS(&user[bi]);
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

        // }
    }
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "orchestrator finished for step %d.\n", simCtx->step);
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
