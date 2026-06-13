#include "momentumsolvers.h"


#undef __FUNCT__
#define __FUNCT__ "ComputeTotalResidual"
/**
 * @brief Internal helper implementation: `ComputeTotalResidual()`.
 * @details Local to this translation unit.
 */
static PetscErrorCode ComputeTotalResidual(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx *simCtx = user->simCtx;
    
    // Extract Time Parameters from Context
    const PetscInt  ti      = simCtx->step;
    const PetscInt  tistart = simCtx->StartStep;
    const PetscReal dt      = simCtx->dt;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    // 1. Calculate Spatial Terms (stored in user->Rhs)
    //    Rhs = -Div(Flux) + Viscous + Source
    ierr = ComputeRHS(user, user->Rhs); CHKERRQ(ierr);

    // 2. Add Physical Time Derivative Terms (BDF Discretization)
    //    The equation solved is: dU/dtau = RHS_Spatial + RHS_Temporal
    
    if (COEF_TIME_ACCURACY > 1.1 && ti != tistart && ti != 1) {
        // --- BDF2 (Second Order Backward Difference) ---
        // Formula: (1.5*U^{n} - 2.0*U^{n-1} + 0.5*U^{n-2}) / dt  = RHS_Spatial(U^{n})
        // Moved to : -1.5/dt * U^{n} + 2.0/dt * U^{n-1} - 0.5/dt * U^{n-2} + RHS_Spatial(U^{n})
        ierr = VecAXPY(user->Rhs, -COEF_TIME_ACCURACY/dt, user->Ucont);
        ierr = VecAXPY(user->Rhs, +2.0/dt, user->Ucont_o);
        ierr = VecAXPY(user->Rhs, -0.5/dt, user->Ucont_rm1);
    } else {
        // --- BDF1 (First Order / Euler Implicit) ---
        // Formula: (U^{n} - U^{n-1}) / dt = RHS_Spatial(U^{n})
        // Moved to RHS: -1.0/dt * U^{n} + 1.0/dt * U^{n-1} + RHS_Spatial(U^{n})
        ierr = VecAXPY(user->Rhs, -1.0/dt, user->Ucont);
        ierr = VecAXPY(user->Rhs, +1.0/dt, user->Ucont_o);
    }

    // 3. Enforce Boundary Conditions on the Residual
    ierr = EnforceRHSBoundaryConditions(user); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "Momentum_Solver_Explicit_RungeKutta4"
/**
 * @brief Internal helper implementation: `MomentumSolver_Explicit_RungeKutta4()`.
 * @details Local to this translation unit.
 */
PetscErrorCode MomentumSolver_Explicit_RungeKutta4(UserCtx *user, IBMNodes *ibm, FSInfo *fsi)
{
    PetscErrorCode ierr;
    (void)fsi;
    // --- Context Acquisition ---
    // Get the master simulation context from the first block's UserCtx.
    // This is the bridge to access all former global variables.
    SimCtx *simCtx = user[0].simCtx;
    PetscReal dt = simCtx->dt;
    //PetscReal st = simCtx->st;
    PetscInt  istage;
    PetscReal alfa[] = {0.25, 1.0/3.0, 0.5, 1.0};

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Executing explicit momentum solver (Runge-Kutta) for %d block(s).\n",simCtx->block_number);

    // --- 1. Pre-Loop Initialization (Legacy Logic) ---
    // This block prepares boundary conditions and allocates the RHS vector for all blocks
    // before the main RK loop begins. This logic is preserved from the original code.
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Preparing all blocks for Runge-Kutta solve...\n");
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {

        ierr = ApplyBoundaryConditions(&user[bi]); CHKERRQ(ierr);

        /*
        // Immersed boundary interpolation (if enabled)
        if (simCtx->immersed) {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "  Performing pre-RK IBM interpolation for block %d.\n", bi);
            for (PetscInt ibi = 0; ibi < simCtx->NumberOfBodies; ibi++) {
                // The 'ibm' and 'fsi' pointers are passed directly from FlowSolver
                ierr = ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi, 1); CHKERRQ(ierr);
            }
        }
        */
        
        // Allocate the persistent RHS vector for this block's context
        ierr = VecDuplicate(user[bi].Ucont, &user[bi].Rhs); CHKERRQ(ierr);
    }
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Pre-loop initialization complete.\n");

    // --- 2. Main Runge-Kutta Loop ---
    // The legacy code had an outer `pseudot` loop that only ran once. We preserve it.
    for (PetscInt pseudot = 0; pseudot < 1; pseudot++) {
        // Loop over each block to perform the RK stages
        for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
            for (istage = 0; istage < 4; istage++) {
                LOG_ALLOW(LOCAL, LOG_DEBUG, "  Block %d, RK Stage %d (alpha=%.4f)...\n", bi, istage, alfa[istage]);

                // a. Calculate the Right-Hand Side (RHS) of the momentum equation.
                ierr = ComputeRHS(&user[bi], user[bi].Rhs); CHKERRQ(ierr);

                // b. Advance Ucont to the next intermediate stage using the RK coefficient.
                //    Ucont_new = Ucont_old + alpha * dt * RHS
                ierr = VecWAXPY(user[bi].Ucont, alfa[istage] * dt, user[bi].Rhs, user[bi].Ucont_o); CHKERRQ(ierr);
 
                // c. Synchronize periodic endpoints and local ghosts for the new intermediate Ucont.
                const char *staggered_fields[] = {"Ucont"};
                ierr = SynchronizePeriodicStaggeredFields(&user[bi], 1, staggered_fields); CHKERRQ(ierr);

                // d. Re-apply boundary conditions for the new intermediate velocity.
                //    This is crucial for the stability and accuracy of the multi-stage scheme.
                ierr = ApplyBoundaryConditions(&user[bi]); CHKERRQ(ierr);
                
            } // End of RK stages for one block

            /*
            // Final IBM Interpolation for the block (if enabled)
            if (simCtx->immersed) {
                LOG_ALLOW(LOCAL, LOG_DEBUG, "  Performing post-RK IBM interpolation for block %d.\n", bi);
                for (PetscInt ibi = 0; ibi < simCtx->NumberOfBodies; ibi++) {
                    ierr = ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi, 1); CHKERRQ(ierr);
                }
            }
            */

        } // End loop over blocks

        // --- 3. Inter-Block Communication (Legacy Logic) ---
        // This is called after all blocks have completed their RK stages.
        if (simCtx->block_number > 1) {
	  //    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Updating multi-block interfaces after RK stages.\n");
	    //   ierr = Block_Interface_U(user); CHKERRQ(ierr);
        }

    } // End of pseudo-time loop

    // --- 4. Cleanup ---
    // Destroy the RHS vectors that were created at the start of this function.
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        ierr = VecDestroy(&user[bi].Rhs); CHKERRQ(ierr);
    }
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Runge-Kutta solve completed for all blocks.\n");

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MomentumSolver_DualTime_Picard_JamesonRK"
/**
 * @brief Internal helper implementation: `MomentumSolver_DualTime_Picard_JamesonRK()`.
 * @details Local to this translation unit.
 */
PetscErrorCode MomentumSolver_DualTime_Picard_JamesonRK(UserCtx *user, IBMNodes *ibm, FSInfo *fsi)
{
    (void)ibm;
    (void)fsi;
    // --- CONTEXT ACQUISITION BLOCK ---
    SimCtx *simCtx = user[0].simCtx;
    const PetscInt block_number = simCtx->block_number;
    
    // Legacy Names (Physics/Time) - Kept for brevity in formulas
    const PetscInt ti = simCtx->step;
    const PetscReal dt = simCtx->dt;       // Physical Time Step
    const PetscReal cfl = simCtx->pseudo_cfl;     // Initial CFL Scaling
    const PetscReal alfa[] = {0.25, 1.0/3.0, 0.5, 1.0}; // Jameson RK smoothing coefficients
    Vec *pRhs; // Per-block backup of Rhs at last accepted pseudo-state.

    // State flags
    PetscBool force_restart = PETSC_FALSE;

    // Renamed Solver Parameters
    const PetscInt  max_pseudo_steps = simCtx->mom_max_pseudo_steps; // Max Dual-Time Iterations
    const PetscReal tol_abs_delta    = simCtx->mom_atol; // Stop if |dU| < tol
    const PetscReal tol_rtol_delta   = simCtx->mom_rtol; // Stop if |dU|/|dU0| < tol
    // --- END CONTEXT ACQUISITION BLOCK ---

    PetscErrorCode ierr;
    PetscMPIInt    rank;
    PetscInt       istage, pseudo_iter;
    PetscInt       accepted_iter = 0, rejected_iter = 0, recovery_streak = 0;
    PetscReal      ts, te, cput;

    // --- Global Convergence Metrics ---
    PetscReal global_norm_delta = 10.0; 
    PetscReal global_rel_delta  = 1.0;  
    PetscReal global_norm_resid = 1.0;
    PetscReal global_rel_resid  = 1.0;  

    // --- Local Metric Arrays (Per Block) ---
    PetscReal *delta_sol_norm_init, *delta_sol_norm_prev, *delta_sol_norm_curr, *delta_sol_rel_curr;
    PetscReal *resid_norm_init,     *resid_norm_prev,     *resid_norm_curr,     *resid_rel_curr;
    PetscReal *pseudo_dt_scaling;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Executing Dual-Time Momentum Solver for %d block(s).\n", block_number);

    // --- Allocate metric arrays ---
    ierr = PetscMalloc2(block_number, &delta_sol_norm_init, block_number, &delta_sol_norm_prev); CHKERRQ(ierr);
    ierr = PetscMalloc2(block_number, &delta_sol_norm_curr, block_number, &delta_sol_rel_curr); CHKERRQ(ierr);
    ierr = PetscMalloc2(block_number, &resid_norm_init,     block_number, &resid_norm_prev); CHKERRQ(ierr);
    ierr = PetscMalloc2(block_number, &resid_norm_curr,     block_number, &pseudo_dt_scaling); CHKERRQ(ierr);
    ierr = PetscMalloc1(block_number, &resid_rel_curr); CHKERRQ(ierr);
    ierr = PetscMalloc1(block_number, &pRhs); CHKERRQ(ierr);

    ierr = PetscTime(&ts); CHKERRQ(ierr);

    // --- 1. Pre-Loop Initialization ---
    
    if (block_number > 1) {
      // ierr = Block_Interface_U(user); CHKERRQ(ierr);
    }

    for (PetscInt bi = 0; bi < block_number; bi++) {
        ierr = ApplyBoundaryConditions(&user[bi]); CHKERRQ(ierr);
        
        // Immersed boundary interpolation (if enabled)

        // Allocate workspace
        ierr = VecDuplicate(user[bi].Ucont, &user[bi].Rhs); CHKERRQ(ierr);
        ierr = VecDuplicate(user[bi].Rhs, &pRhs[bi]); CHKERRQ(ierr);
        ierr = VecDuplicate(user[bi].Ucont, &user[bi].dUcont); CHKERRQ(ierr);
        ierr = VecDuplicate(user[bi].Ucont, &user[bi].pUcont); CHKERRQ(ierr);
        
        // Initialize Backup (pUcont) with current state
        ierr = VecCopy(user[bi].Ucont, user[bi].pUcont); CHKERRQ(ierr);

        // --- Calculate Initial Total Residual (Spatial + Temporal) ---
        ierr = ComputeTotalResidual(&user[bi]); CHKERRQ(ierr);

        // Backup Rhs with current state
        ierr = VecCopy(user[bi].Rhs, pRhs[bi]); CHKERRQ(ierr);

        // Compute Initial Norms
        ierr = VecNorm(user[bi].Rhs, NORM_INFINITY, &resid_norm_init[bi]); CHKERRQ(ierr);
        
        // Initialize history for backtracking logic
        resid_norm_prev[bi]     = resid_norm_init[bi]; // Meaningful ratio on first iteration
        delta_sol_norm_prev[bi] = 1000.0;
        pseudo_dt_scaling[bi] = cfl; // Initialize adaptive scalar

	    LOG_ALLOW(GLOBAL,LOG_INFO," Block %d | Max RHS = %.6f | Pseudo-CFL = %.4f .\n", bi, resid_norm_init[bi], cfl);
    }
    
    // --- 2. Main Pseudo-Time Iteration Loop ---
    pseudo_iter = 0;
    PetscBool residual_convergence_enabled =
        (PetscBool)(simCtx->mom_resid_atol > 0.0 || simCtx->mom_resid_rtol > 0.0);
    PetscBool converged = PETSC_FALSE;
    PetscBool last_trial_nonfinite = PETSC_FALSE;
    while (!converged && pseudo_iter < max_pseudo_steps)
    {
        pseudo_iter++;
        force_restart = PETSC_FALSE;

        // Every attempt starts from the last globally accepted state.
        for (PetscInt bi = 0; bi < block_number; bi++) {
            ierr = VecCopy(user[bi].pUcont, user[bi].Ucont); CHKERRQ(ierr);
            ierr = VecCopy(pRhs[bi], user[bi].Rhs); CHKERRQ(ierr);
            const char *staggered_fields[] = {"Ucont"};
            ierr = SynchronizePeriodicStaggeredFields(&user[bi], 1, staggered_fields); CHKERRQ(ierr);
            ierr = ApplyBoundaryConditions(&user[bi]); CHKERRQ(ierr);
        }
	
        for (PetscInt bi = 0; bi < block_number; bi++) {
            
            // === 4-Stage Jameson RK Smoothing Loop ===
            for (istage = 0; istage < 4; istage++) {
                
                LOG_ALLOW(GLOBAL,LOG_TRACE," Pseudo-Iter: %d | RK-Stage: %d\n", pseudo_iter, istage);
                
                // RK Update: U_new = U_old + (Scaler * Alpha * dt) * Total_Residual
                ierr = VecWAXPY(user[bi].Ucont, 
                                pseudo_dt_scaling[bi] * alfa[istage] * dt, 
                                user[bi].Rhs, 
                                user[bi].pUcont); CHKERRQ(ierr);

                // Sync Ghosts & Re-apply BCs for intermediate stage
                const char *staggered_fields[] = {"Ucont"};
                ierr = SynchronizePeriodicStaggeredFields(&user[bi], 1, staggered_fields); CHKERRQ(ierr);
                
                ierr = ApplyBoundaryConditions(&user[bi]); CHKERRQ(ierr);

                // --- Re-calculate Total Residual for next stage ---
                ierr = ComputeTotalResidual(&user[bi]); CHKERRQ(ierr);
      
            } // End RK Stages

            // Immersed boundary interpolation (if enabled)

            // === Convergence Metrics Calculation ===
            
            // Calculate dU = U_current - U_backup
            ierr = VecWAXPY(user[bi].dUcont, -1.0, user[bi].pUcont, user[bi].Ucont); CHKERRQ(ierr);

            // Compute Infinity Norms
            ierr = VecNorm(user[bi].dUcont, NORM_INFINITY, &delta_sol_norm_curr[bi]); CHKERRQ(ierr);
            ierr = VecNorm(user[bi].Rhs, NORM_INFINITY, &resid_norm_curr[bi]); CHKERRQ(ierr);
	    
            // Normalize relative metrics
            if (pseudo_iter == 1) {
                delta_sol_norm_init[bi] = delta_sol_norm_curr[bi];
                delta_sol_rel_curr[bi]  = 1.0;
                resid_rel_curr[bi]      = 1.0;
                // resid_norm_init[bi] set correctly before the loop — do not overwrite
            } else {
                if (delta_sol_norm_init[bi] > 1.0e-10) 
                    delta_sol_rel_curr[bi] = delta_sol_norm_curr[bi] / delta_sol_norm_init[bi];
                else 
                    delta_sol_rel_curr[bi] = 0.0;

                if(resid_norm_init[bi] > 1.0e-10) 
                    resid_rel_curr[bi] = resid_norm_curr[bi] / resid_norm_init[bi];
                else 
                    resid_rel_curr[bi] = 0.0;
	        }
	      
            // --- File Logging ---
            if (!rank) {
                FILE *f;
                char filen[PETSC_MAX_PATH_LEN + 128];
                ierr = PetscSNPrintf(filen, sizeof(filen), "%s/Momentum_Solver_Convergence_History_Block_%1d.log", simCtx->log_dir, bi); CHKERRQ(ierr);
                if(simCtx->step == simCtx->StartStep + 1 && pseudo_iter == 1 && !simCtx->continueMode) f = fopen(filen, "w");
                else f = fopen(filen, "a");
                if (simCtx->continueMode && simCtx->step == simCtx->StartStep + 1 && pseudo_iter == 1) {
                    PetscFPrintf(PETSC_COMM_WORLD, f, "# Continuation from step %" PetscInt_FMT "\n", simCtx->StartStep);
                }
                PetscFPrintf(PETSC_COMM_WORLD, f, "Step: %d | PseudoIter(k): %d| | Pseudo-cfl: %.4f |dUk|: %le | |dUk|/|dUprev|: %le | |Rk|: %le | |Rk|/|Rprev|: %le \n",
                             (int)ti, (int)pseudo_iter, pseudo_dt_scaling[bi], delta_sol_norm_curr[bi], delta_sol_rel_curr[bi], resid_norm_curr[bi],resid_rel_curr[bi]);
                fclose(f);
            }
        } // End loop over blocks

        // --- Update Global Convergence Criteria ---
        global_norm_delta = -1.0e20; 
        global_rel_delta  = -1.0e20; 
        global_norm_resid = -1.0e20;
        global_rel_resid  = -1.0e20;
        
        for (PetscInt bi = 0; bi < block_number; bi++) {
            global_norm_delta = PetscMax(delta_sol_norm_curr[bi], global_norm_delta);
            global_rel_delta  = PetscMax(delta_sol_rel_curr[bi],  global_rel_delta);
            global_norm_resid = PetscMax(resid_norm_curr[bi],      global_norm_resid);
            global_rel_resid  = PetscMax(resid_rel_curr[bi],      global_rel_resid);
        }
        ierr = PetscTime(&te); CHKERRQ(ierr);
        cput = te - ts;
        LOG_ALLOW(GLOBAL, LOG_INFO, "  Pseudo-Iter(k) %d: |dUk|=%e, |dUk|/|dU0| = %e, |Rk|/|R0| = %e, CPU=%.2fs\n",
              pseudo_iter, global_norm_delta, global_rel_delta, global_rel_resid, cput);

        // === Adaptive Pseudo-CFL Trial Acceptance and Rollback ===
        const PetscReal resid_floor = 1.0e-30;
        PetscReal global_trial_ratio = 0.0;
        PetscBool global_nonfinite = PETSC_FALSE;
        for (PetscInt bi = 0; bi < block_number; bi++) {
            PetscReal ratio;
            if (resid_norm_prev[bi] > resid_floor) {
                ratio = resid_norm_curr[bi] / resid_norm_prev[bi];
            } else if (resid_norm_curr[bi] <= resid_floor) {
                ratio = 0.0;
            } else {
                ratio = PETSC_MAX_REAL;
            }
            global_trial_ratio = PetscMax(global_trial_ratio, ratio);
            global_nonfinite = (PetscBool)(global_nonfinite ||
                PetscIsInfOrNanReal(delta_sol_norm_curr[bi]) ||
                PetscIsInfOrNanReal(resid_norm_curr[bi]) ||
                PetscIsInfOrNanReal(ratio));
        }

        PetscReal reduced_trial_ratio;
        PetscMPIInt local_nonfinite = global_nonfinite ? 1 : 0;
        PetscMPIInt reduced_nonfinite = 0;
        ierr = MPI_Allreduce(&global_trial_ratio, &reduced_trial_ratio, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
        ierr = MPI_Allreduce(&local_nonfinite, &reduced_nonfinite, 1, MPI_INT, MPI_LOR, PETSC_COMM_WORLD); CHKERRQ(ierr);
        global_trial_ratio = reduced_trial_ratio;
        global_nonfinite = reduced_nonfinite ? PETSC_TRUE : PETSC_FALSE;
        last_trial_nonfinite = global_nonfinite;
        force_restart = (PetscBool)(global_nonfinite ||
            global_trial_ratio > simCtx->mom_dt_jameson_residual_norm_noise_allowance_factor);

        if (force_restart) {
            PetscReal old_cfl = pseudo_dt_scaling[0];
            PetscReal next_cfl = PetscMax(simCtx->min_pseudo_cfl,
                old_cfl * simCtx->pseudo_cfl_reduction_factor);
            rejected_iter++;
            recovery_streak = 0;
            for (PetscInt bi = 0; bi < block_number; bi++) {
                ierr = VecCopy(user[bi].pUcont, user[bi].Ucont); CHKERRQ(ierr);
                ierr = VecCopy(pRhs[bi], user[bi].Rhs); CHKERRQ(ierr);
                const char *staggered_fields[] = {"Ucont"};
                ierr = SynchronizePeriodicStaggeredFields(&user[bi], 1, staggered_fields); CHKERRQ(ierr);
                ierr = ApplyBoundaryConditions(&user[bi]); CHKERRQ(ierr);
                pseudo_dt_scaling[bi] = next_cfl;
            }
            LOG_ALLOW(GLOBAL, LOG_DEBUG,
                      "  Momentum trial %d rejected (ratio=%e, nonfinite=%d); CFL %.6f -> %.6f.\n",
                      pseudo_iter, global_trial_ratio, (int)global_nonfinite, old_cfl, next_cfl);
            if (global_nonfinite && old_cfl <= simCtx->min_pseudo_cfl) {
                SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_CONV_FAILED,
                        "Momentum solver produced a non-finite trial at minimum pseudo-CFL.");
            }
            continue;
        }

        accepted_iter++;
        last_trial_nonfinite = PETSC_FALSE;
        for (PetscInt bi = 0; bi < block_number; bi++) {
            ierr = VecCopy(user[bi].Ucont, user[bi].pUcont); CHKERRQ(ierr);
            ierr = VecCopy(user[bi].Rhs, pRhs[bi]); CHKERRQ(ierr);
            resid_norm_prev[bi] = resid_norm_curr[bi];
            delta_sol_norm_prev[bi] = delta_sol_norm_curr[bi];
        }

        PetscReal next_cfl = pseudo_dt_scaling[0];
        if (global_trial_ratio < 0.90) {
            next_cfl *= simCtx->pseudo_cfl_growth_factor;
            recovery_streak = 0;
        } else if (global_trial_ratio <= 1.0) {
            recovery_streak++;
            if (recovery_streak >= 3) {
                next_cfl *= simCtx->pseudo_cfl_growth_factor;
                recovery_streak = 0;
            }
        } else {
            next_cfl *= simCtx->pseudo_cfl_reduction_factor;
            recovery_streak = 0;
        }
        next_cfl = PetscMin(next_cfl, simCtx->max_pseudo_cfl);
        next_cfl = PetscMax(next_cfl, simCtx->min_pseudo_cfl);
        for (PetscInt bi = 0; bi < block_number; bi++) pseudo_dt_scaling[bi] = next_cfl;

        PetscBool update_pass;
        PetscBool residual_pass = PETSC_TRUE;
        if (residual_convergence_enabled) {
            update_pass = (PetscBool)(global_norm_delta <= tol_abs_delta || global_rel_delta <= tol_rtol_delta);
            residual_pass = (PetscBool)(
                (simCtx->mom_resid_atol > 0.0 && global_norm_resid <= simCtx->mom_resid_atol) ||
                (simCtx->mom_resid_rtol > 0.0 && global_rel_resid <= simCtx->mom_resid_rtol));
        } else {
            update_pass = (PetscBool)(global_norm_delta <= tol_abs_delta && global_rel_delta <= tol_rtol_delta);
        }
        converged = (PetscBool)(update_pass && residual_pass);

        if (block_number > 1) {
             // ierr = Block_Interface_U(user); CHKERRQ(ierr);
        }
    } // End while loop

    if (last_trial_nonfinite) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_CONV_FAILED,
                "Momentum solver exhausted its attempt limit while recovering from a non-finite trial.");
    }

    // Carry the controller-selected next CFL directly into the next physical step.
    PetscReal next_start_cfl = pseudo_dt_scaling[0];

    // 5. Save for next Physical Time Step
    if (!rank) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "  CFL Adaptation: Step %d finished at CFL %.3f. Next step will start at CFL %.3f.\n", 
                  simCtx->step, pseudo_dt_scaling[0], next_start_cfl);
    }
    
    // Update the context so the next call to MomentumSolver reads this new value
    simCtx->pseudo_cfl = next_start_cfl;    

    // Ensure every nonfatal exit exposes the last accepted finite state.
    for (PetscInt bi = 0; bi < block_number; bi++) {
        ierr = VecCopy(user[bi].pUcont, user[bi].Ucont); CHKERRQ(ierr);
        ierr = VecCopy(pRhs[bi], user[bi].Rhs); CHKERRQ(ierr);
        const char *staggered_fields[] = {"Ucont"};
        ierr = SynchronizePeriodicStaggeredFields(&user[bi], 1, staggered_fields); CHKERRQ(ierr);
        ierr = ApplyBoundaryConditions(&user[bi]); CHKERRQ(ierr);
    }

    if (!converged) {
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "Momentum solver reached %d attempts without convergence; continuing from the last accepted finite state.\n",
                  pseudo_iter);
    }
    if (accepted_iter == 0) {
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "Momentum solver accepted no pseudo-time trials; retaining the physical-step entry state.\n");
    }
    LOG_ALLOW(GLOBAL, LOG_INFO,
              "Momentum solver finished after %d attempts (%d accepted, %d rejected), next CFL=%.6f.\n",
              pseudo_iter, accepted_iter, rejected_iter, next_start_cfl);

    // --- Final Cleanup ---
    for (PetscInt bi = 0; bi < block_number; bi++) {
        ierr = VecDestroy(&user[bi].Rhs); CHKERRQ(ierr);
        ierr = VecDestroy(&user[bi].dUcont); CHKERRQ(ierr);
        ierr = VecDestroy(&user[bi].pUcont); CHKERRQ(ierr);
        ierr = VecDestroy(&pRhs[bi]); CHKERRQ(ierr);
    }
    ierr = PetscFree(pRhs); CHKERRQ(ierr);

    ierr = PetscFree2(delta_sol_norm_init, delta_sol_norm_prev);CHKERRQ(ierr);
    ierr = PetscFree2(delta_sol_norm_curr, delta_sol_rel_curr);CHKERRQ(ierr);
    ierr = PetscFree2(resid_norm_init,     resid_norm_prev);CHKERRQ(ierr);
    ierr = PetscFree2(resid_norm_curr,     pseudo_dt_scaling); CHKERRQ(ierr);
    ierr = PetscFree(resid_rel_curr); CHKERRQ(ierr);
    
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
