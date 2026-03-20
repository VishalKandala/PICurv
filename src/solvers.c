#include "solvers.h" // The new header we will create

#undef __FUNCT__
#define __FUNCT__ "FlowSolver"
/**
 * @brief Implementation of \ref FlowSolver().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/solvers.h`.
 * @see FlowSolver()
 */
PetscErrorCode FlowSolver(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserMG         *usermg = NULL;
    PetscInt       level;
    UserCtx        *user = NULL;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    if (simCtx->mom_solver_type != MOMENTUM_SOLVER_DUALTIME_PICARD_RK4 &&
        simCtx->mom_solver_type != MOMENTUM_SOLVER_EXPLICIT_RK) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "Unknown momentum solver type %d. Supported values are EXPLICIT_RK and DUALTIME_PICARD_RK4.",
                simCtx->mom_solver_type);
    }

    usermg = &simCtx->usermg;
    level  = usermg->mglevels - 1;
    user   = usermg->mgctx[level].user;
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
            ierr = UpdateLocalGhosts(&user[bi], "Ucont");
            ierr = Contra2Cart(&user[bi]); CHKERRQ(ierr);
            ierr = UpdateLocalGhosts(&user[bi], "Ucat");
            if(simCtx->les == CONSTANT_SMAGORINSKY) {
                LOG_ALLOW(LOCAL, LOG_INFO, "  Using constant Smagorinsky model for block %d.\n", bi);
                // Constant Smagorinsky does not require dynamic computation
                if(simCtx->step == simCtx->StartStep + 1){ // since step is updated before flowsolver call.
                    ierr = ComputeSmagorinskyConstant(&user[bi]); CHKERRQ(ierr);
                }    
            } else if(simCtx->les == DYNAMIC_SMAGORINSKY) {
            if (simCtx->step % simCtx->dynamic_freq == 0) {
                LOG_ALLOW(LOCAL, LOG_DEBUG, "  Computing Smagorinsky constant for block %d.\n", bi);
                ierr = ComputeSmagorinskyConstant(&user[bi]);
                }
            }
          //  LOG_ALLOW(LOCAL, LOG_DEBUG, "  Computing eddy viscosity for block %d.\n", bi);
          ierr = ComputeEddyViscosityLES(&user[bi]);
        }
    }
    

    // ========================================================================
    //   SECTION: Momentum Equation Solver
    // ========================================================================
    // This is the core of the time step. It computes an intermediate velocity
    // field by solving the momentum equations.
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Beginning momentum step solve (Solver = %s)...\n", MomentumSolverTypeToString(simCtx->mom_solver_type));
    
    // Since IBM is disabled, we pass NULL for ibm and fsi arguments.
    // ierr = ImpRK(user, NULL, NULL); CHKERRQ(ierr);
    // Add new momentum solver types here only after wiring the enum, parser, docs, and tests.
    if(simCtx->mom_solver_type == MOMENTUM_SOLVER_DUALTIME_PICARD_RK4) {
        ierr = MomentumSolver_DualTime_Picard_RK4(user,NULL,NULL); CHKERRQ(ierr);
    } else if(simCtx->mom_solver_type == MOMENTUM_SOLVER_EXPLICIT_RK) {   
    // Since IBM is disabled, we pass NULL for ibm and fsi arguments.
    ierr = MomentumSolver_Explicit_RungeKutta4(user, NULL, NULL); CHKERRQ(ierr);
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
   SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
           "Unsupported Poisson solver type %d. The current runtime supports only the multigrid path (poisson = 0).",
           simCtx->poisson);
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
