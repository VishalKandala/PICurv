#include "momentumsolvers.h"


#undef __FUNCT__
#define __FUNCT__ "ComputeTotalResidual"
/**
 * @brief Computes the Total Residual for the Dual-Time stepping method.
 *
 * @details 
 * Calculates R_total(U) = R_spatial(U) + R_temporal(U).
 * 1. Calls ComputeRHS() to get spatial fluxes (Convection, Diffusion, Pressure).
 * 2. Adds the Physical Time Derivative terms (BDF1 or BDF2) to the RHS vector.
 *    Time parameters (ti, dt) are accessed directly from user->simCtx.
 * 3. Enforces RHS boundary conditions.
 *
 * @param user Pointer to the User context for the specific block.
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
        // Formula: (1.5*U^{n} - 2.0*U^{n-1} + 0.5*U^{n-2}) / dt
        // Moved to RHS: -1.5/dt * U^{n} + ...
        ierr = VecAXPY(user->Rhs, -COEF_TIME_ACCURACY/dt, user->Ucont);
        ierr = VecAXPY(user->Rhs, +2.0/dt, user->Ucont_o);
        ierr = VecAXPY(user->Rhs, -0.5/dt, user->Ucont_rm1);
    } else {
        // --- BDF1 (First Order / Euler Implicit) ---
        // Formula: (U^{n} - U^{n-1}) / dt
        // Moved to RHS: -1.0/dt * U^{n} + ...
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
 * @brief Advances the momentum equations for ALL blocks by one time step
 *        using an explicit 4-stage, 4th-order Runge-Kutta scheme.
 *
 * This function computes an intermediate, non-divergence-free contravariant
 * velocity field (Ucont) at time t_{n+1} for all computational blocks.
 *
 * This is a minimally-edited version of the legacy solver. It retains its
 * internal loop over all blocks and is intended to be called once per time step
 * from the main FlowSolver orchestrator. All former global variables are now
 * accessed via the SimCtx passed in through the first block's UserCtx.
 * 
 * @note This does not employ a Backward Euler implicit treatment of the time derivative.
 * The physical timestep itself is advanced explicitly. The dual-time stepping is not used here.
 *
 * @param user The array of UserCtx structs for all blocks.
 * @param ibm  (Optional) Pointer to the full array of IBM data structures. Pass NULL if disabled.
 * @param fsi  (Optional) Pointer to the full array of FSI data structures. Pass NULL if disabled.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode MomentumSolver_Explicit_RungeKutta4(UserCtx *user, IBMNodes *ibm, FSInfo *fsi)
{
    PetscErrorCode ierr;
    // --- Context Acquisition ---
    // Get the master simulation context from the first block's UserCtx.
    // This is the bridge to access all former global variables.
    SimCtx *simCtx = user[0].simCtx;
    PetscReal dt = simCtx->dt;
    PetscReal st = simCtx->st;
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
                ierr = VecWAXPY(user[bi].Ucont, alfa[istage] * dt * st, user[bi].Rhs, user[bi].Ucont_o); CHKERRQ(ierr);
 
                // c. Update local ghost values for the new intermediate Ucont.
                ierr = DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); CHKERRQ(ierr);
                ierr = DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); CHKERRQ(ierr);

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
#define __FUNCT__ "ImpRK"
/**
 * @brief Advances the momentum equations using an iterative explicit Runge-Kutta scheme.
 *
 * This function uses a pseudo-time-stepping approach. It repeatedly applies an
 * explicit RK scheme within a while loop until the change between iterations
 * converges below specified tolerances (imp_atol, imp_rtol) or a maximum
 * number of iterations (imp_MAX_IT) is reached.
 *
 * This refactored version keeps the original internal loop over all blocks. All
 * former global variables are accessed via the SimCtx. It directly calculates the
 * full RHS (spatial + temporal terms) instead of calling a separate CalcRHS function.
 *
 * @param user The array of UserCtx structs for all blocks.
 * @param ibm  (Optional) Pointer to the full array of IBM data structures. Pass NULL if disabled.
 * @param fsi  (Optional) Pointer to the full array of FSI data structures. Pass NULL if disabled.
 * @return PetscErrorCode 0 on success.
 *
PetscErrorCode ImpRK(UserCtx *user, IBMNodes *ibm, FSInfo *fsi)
{
    // --- CONTEXT ACQUISITION BLOCK ---
    SimCtx *simCtx = user[0].simCtx;
    const PetscInt block_number = simCtx->block_number;
    const PetscInt immersed = simCtx->immersed;
    const PetscInt ti = simCtx->step;
    const PetscInt tistart = simCtx->StartStep;
    const PetscInt imp_MAX_IT = simCtx->imp_MAX_IT;
    const PetscReal imp_atol = simCtx->imp_atol;
    const PetscReal imp_rtol = simCtx->imp_rtol;
    const PetscReal dt = simCtx->dt;
    const PetscReal st = simCtx->st;
    const PetscReal cfl = simCtx->cfl;
    // --- END CONTEXT ACQUISITION BLOCK ---

    PetscErrorCode ierr;
    PetscMPIInt    rank;
    PetscInt       istage, pseudot, ibi;
    PetscReal      alfa[] = {0.25, 1.0/3.0, 0.5, 1.0};
    PetscReal      ts, te, cput;

    // --- Convergence Metrics ---
    PetscReal normdU = 10.0, reldU = 1.0, relF = 1.0;
    PetscReal *normdU0_bk, *normdU1_bk, *normdU_bk, *reldU_bk;
    PetscReal *normF0_bk, *normF1_bk, *normF_bk, *relF_bk, *lambda;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Executing implicit-like RK solver (ImpRK) for %d block(s).\n", block_number);

    // --- Allocate metric arrays ---
    ierr = PetscMalloc2(block_number, &normdU0_bk, block_number, &normdU1_bk); CHKERRQ(ierr);
    ierr = PetscMalloc2(block_number, &normdU_bk, block_number, &reldU_bk); CHKERRQ(ierr);
    ierr = PetscMalloc2(block_number, &normF0_bk,  block_number, &normF1_bk); CHKERRQ(ierr);
    ierr = PetscMalloc2(block_number, &normF_bk, block_number, &lambda); CHKERRQ(ierr);
    ierr = PetscMalloc1(block_number, &relF_bk); CHKERRQ(ierr);

    ierr = PetscTime(&ts); CHKERRQ(ierr);

    // --- 1. Pre-Loop Initialization (Legacy Logic) ---
    if (block_number > 1) {
      //       ierr = Block_Interface_U(user); CHKERRQ(ierr);
    }

    for (PetscInt bi = 0; bi < block_number; bi++) {
        // Prepare BCs and workspace vectors

        ierr = ApplyBoundaryConditions(&user[bi]); CHKERRQ(ierr);
        
        
        //if (immersed) {
        //    for (ibi = 0; ibi < simCtx->NumberOfBodies; ibi++) {
        //        ierr = ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi, 1); CHKERRQ(ierr);
        //    }
        //}
        
        
        // Allocate workspace vectors
        ierr = VecDuplicate(user[bi].Ucont, &user[bi].Rhs); CHKERRQ(ierr);
        ierr = VecDuplicate(user[bi].Ucont, &user[bi].dUcont); CHKERRQ(ierr);
        ierr = VecDuplicate(user[bi].Ucont, &user[bi].pUcont); CHKERRQ(ierr);
        // Initialize pUcont with current Ucont
        ierr = VecCopy(user[bi].Ucont, user[bi].pUcont); CHKERRQ(ierr);


	LOG_ALLOW(GLOBAL,LOG_DEBUG," BCs and workspace vectors prepared for Initial RHS calculation.\n");
	
        // --- Calculate INITIAL RHS for convergence check ---
        ierr = ComputeRHS(&user[bi], user[bi].Rhs); CHKERRQ(ierr);

	LOG_ALLOW(GLOBAL,LOG_DEBUG, " Initial RHS calculated for convergence check .\n");
	
        // Add time derivative part
        if (COEF_TIME_ACCURACY > 1.1 && ti != tistart && ti != 1) {
            // Second Order Backward Difference (BDF2)
            ierr = VecAXPY(user[bi].Rhs, -COEF_TIME_ACCURACY/dt, user[bi].Ucont); CHKERRQ(ierr);
            ierr = VecAXPY(user[bi].Rhs, +2.0/dt, user[bi].Ucont_o); CHKERRQ(ierr);
            ierr = VecAXPY(user[bi].Rhs, -0.5/dt, user[bi].Ucont_rm1); CHKERRQ(ierr);
        } else {
            // For Initial Step: First Order Backward Difference (BDF1)
            ierr = VecAXPY(user[bi].Rhs, -1.0/dt, user[bi].Ucont); CHKERRQ(ierr);
            ierr = VecAXPY(user[bi].Rhs, +1.0/dt, user[bi].Ucont_o); CHKERRQ(ierr);
        }


    // Enforce RHS Boundary Condition (Set to Zero)
    ierr = EnforceRHSBoundaryConditions(&user[bi]); CHKERRQ(ierr);    

	LOG_ALLOW(GLOBAL,LOG_DEBUG," Time derivative part added to RHS.\n");
        ierr = VecNorm(user[bi].Rhs, NORM_INFINITY, &normF0_bk[bi]); CHKERRQ(ierr);
        normF1_bk[bi] = 1000.0; // Initialize for line search logic
        normdU1_bk[bi] = 1000.0;
        lambda[bi] = cfl; // Initialize adaptive step size

	LOG_ALLOW(GLOBAL,LOG_INFO," Block %d | Max RHS = %.6f | CFL = %.4f .\n", bi,normF0_bk[bi],cfl);
    }
    
    // --- 2. Main Pseudo-Time Iteration Loop ---
    pseudot = 0;
    while (((normdU > imp_atol && reldU > imp_rtol) || pseudot < 1) && pseudot < imp_MAX_IT) {
        pseudot++;
	
        for (PetscInt bi = 0; bi < block_number; bi++) {
            for (istage = 0; istage < 4; istage++) {
                // Advance in time using RK scheme with adaptive step size `lambda`
	            LOG_ALLOW(GLOBAL,LOG_TRACE," Pseudo-Timestep Solver | RK-Stage : %d | Pseudo-Timestep :%d .\n",istage,pseudot);
                ierr = VecWAXPY(user[bi].Ucont, lambda[bi] * alfa[istage] * dt * st, user[bi].Rhs, user[bi].pUcont); CHKERRQ(ierr);

		        LOG_ALLOW(GLOBAL,LOG_VERBOSE, " Ucont updated as Ucont = [lambda(%.4f)]x[alfa(%.4f)]x[dt(%.4f)]x[st(%.4f)]xRhs + pUcont.\n",lambda[bi],alfa[istage],dt,st);
		
                // Sync ghosts and re-apply BCs for the intermediate stage
                ierr = DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); CHKERRQ(ierr);
                ierr = DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); CHKERRQ(ierr);

                ierr = ApplyBoundaryConditions(&user[bi]); CHKERRQ(ierr);

		        LOG_ALLOW(GLOBAL,LOG_TRACE, " Ghosts synced and BCs reapplied.\n");

                // --- Re-calculate the full RHS for the next RK stage ---
                ierr = ComputeRHS(&user[bi], user[bi].Rhs); CHKERRQ(ierr);

		        LOG_ALLOW(GLOBAL,LOG_TRACE, " RHS calculated for the next RK stage.\n");
		
                if (COEF_TIME_ACCURACY > 1.1 && ti != tistart && ti != 1) {
                    ierr = VecAXPY(user[bi].Rhs, -COEF_TIME_ACCURACY/dt, user[bi].Ucont); CHKERRQ(ierr);
                    ierr = VecAXPY(user[bi].Rhs, +2.0/dt, user[bi].Ucont_o); CHKERRQ(ierr);
                    ierr = VecAXPY(user[bi].Rhs, -0.5/dt, user[bi].Ucont_rm1); CHKERRQ(ierr);
                } else {
                    ierr = VecAXPY(user[bi].Rhs, -1.0/dt, user[bi].Ucont); CHKERRQ(ierr);
                    ierr = VecAXPY(user[bi].Rhs, +1.0/dt, user[bi].Ucont_o); CHKERRQ(ierr);
                }


                // Enforce RHS Boundary Condition (Set to Zero)
                ierr = EnforceRHSBoundaryConditions(&user[bi]); CHKERRQ(ierr);        

		        if(istage <3) LOG_ALLOW(GLOBAL,LOG_TRACE," Time derivative part updated for next RK stage .\n");
		        else LOG_ALLOW(GLOBAL,LOG_DEBUG," All RK stages complete.\n");
            } // End RK stages

            
            //if (immersed) {
            //    for (ibi = 0; ibi < simCtx->NumberOfBodies; ibi++) {
            //        ierr = ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi, 1); CHKERRQ(ierr);
            //    }
            //}
            

            // --- Calculate Convergence Metrics for this block ---
	    
            ierr = VecWAXPY(user[bi].dUcont, -1.0, user[bi].Ucont, user[bi].pUcont); CHKERRQ(ierr);

	    LOG_ALLOW(GLOBAL,LOG_TRACE," Block %d | dU calculated .\n",bi);
            
            normdU1_bk[bi] = normdU_bk[bi];
	    
            ierr = VecNorm(user[bi].dUcont, NORM_INFINITY, &normdU_bk[bi]); CHKERRQ(ierr);
            ierr = VecNorm(user[bi].Rhs, NORM_INFINITY, &normF_bk[bi]); CHKERRQ(ierr);
	    
	    LOG_ALLOW(GLOBAL,LOG_TRACE," Block %d | |dU|  =  %le | |U| = %le.\n",bi,normdU_bk[bi],normF_bk[bi]);
	    
            if (pseudot == 1) {
                normdU0_bk[bi] = normdU_bk[bi];
                reldU_bk[bi] = 1.0;
	        normF0_bk[bi] = normF_bk[bi];
		relF_bk[bi] = 1.0;
            } else { //pseudot > 1
	      if (normdU0_bk[bi] > 1.0e-10) {
                reldU_bk[bi] = normdU_bk[bi] / normdU0_bk[bi];
	      }
              else {
                reldU_bk[bi] = 0.0;
	      }
	      if(normF0_bk[bi] > 1.0e-10) {
		relF_bk[bi] = normF_bk[bi] / normF0_bk[bi];
	      }else {
              relF_bk[bi] = 0.0;
	      }
	    }
	      

	    LOG_ALLOW(GLOBAL,LOG_TRACE, " Block %d | |dU|/|dU_0| = %le. \n",bi,reldU_bk[bi]);

	    LOG_ALLOW(GLOBAL,LOG_TRACE, " Block %d | |R|/|R_0| = %le. \n",bi,relF_bk[bi]);
	    
            // file logging
            if (!rank) {
                FILE *f;
                char filen[80];
                sprintf(filen, "%s/Momentum_Solver_Convergence_History_Block_%1d.log", simCtx->log_dir,bi);
                if(simCtx->step == simCtx->StartStep + 1 && pseudot ==1) {
                    f = fopen(filen, "w");
                } else {
                f = fopen(filen, "a");
                }
                PetscFPrintf(PETSC_COMM_WORLD, f, "Block %d | Step: %d | Pseudo-Timestep %d | |dU|  %le |  |dU|/|dU_0| %le | |U| = %le \n",bi,(int)ti, (int)pseudot, normdU_bk[bi], reldU_bk[bi], normF_bk[bi]);
                fclose(f);
            }

	    LOG_ALLOW(GLOBAL,LOG_DEBUG, " Block %d | Step: %d | Pseudo-Timestep %d | |dU|  %le |  |dU|/|dU_0| %le | |U| = %le \n",bi,ti,pseudot,normdU_bk[bi], reldU_bk[bi], normF_bk[bi]);
        } // End loop over blocks

        // --- Update Global Convergence Criteria and Perform Line Search ---
        normdU = -1.0e20; reldU = -1.0e20; relF = -1.0e20;
        for (PetscInt bi = 0; bi < block_number; bi++) {
            normdU = PetscMax(normdU_bk[bi], normdU);
            reldU = PetscMax(reldU_bk[bi], reldU);
            relF = PetscMax(relF_bk[bi], relF);
        }

        ierr = PetscTime(&te); CHKERRQ(ierr);
        cput = te - ts;
        LOG_ALLOW(GLOBAL, LOG_INFO, "  Iter(Pseudo-Timestep) %d: |dU|=%e, |dU|/|dU_0|=%e, |U|=%e, CPU=%.2fs\n",
                  pseudot, normdU, reldU, relF, cput);

        for (PetscInt bi = 0; bi < block_number; bi++) {
            // Adaptive step size logic (line search)
            if (pseudot > 1 && normF_bk[bi] > normF1_bk[bi] && normdU_bk[bi] > normdU1_bk[bi] && lambda[bi] > 0.01) {
                lambda[bi] *= 0.5;
                pseudot--; // Retry this iteration with a smaller step
                LOG_ALLOW(LOCAL, LOG_DEBUG, "  Block %d: Line search failed. Reducing lambda to %f and retrying iter.\n", bi, lambda[bi]);
                ierr = VecCopy(user[bi].pUcont, user[bi].Ucont); CHKERRQ(ierr);
                ierr = DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); CHKERRQ(ierr);
                ierr = DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); CHKERRQ(ierr);
            } else {
                ierr = VecCopy(user[bi].Ucont, user[bi].pUcont); CHKERRQ(ierr);
                normF1_bk[bi] = normF_bk[bi];
                normdU1_bk[bi] = normdU_bk[bi];
            }
        }
        
        if (block_number > 1) {
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "  Updating multi-block interfaces after iteration %d.\n", pseudot);
	    //  ierr = Block_Interface_U(user); CHKERRQ(ierr);
        }
    } // End while loop

    // --- Final Cleanup ---
    
    //if (block_number > 1) {
    //    Block_Interface_Correction(user);
    //    if (simCtx->blank) {
    //        Block_Blank_Correction_adv(&user[0], 0);
    //    }
    //}
    
    for (PetscInt bi = 0; bi < block_number; bi++) {
        ierr = VecDestroy(&user[bi].Rhs); CHKERRQ(ierr);
        ierr = VecDestroy(&user[bi].dUcont); CHKERRQ(ierr);
        ierr = VecDestroy(&user[bi].pUcont); CHKERRQ(ierr);
    }
    
    ierr = PetscFree2(normdU0_bk, normdU1_bk);CHKERRQ(ierr);
    ierr = PetscFree2(normdU_bk, reldU_bk);CHKERRQ(ierr);
    ierr = PetscFree2(normF0_bk, normF1_bk);CHKERRQ(ierr);
    ierr = PetscFree2(normF_bk, lambda); CHKERRQ(ierr);
    ierr = PetscFree(relF_bk); CHKERRQ(ierr);
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "ImpRK solve completed after %d iterations.\n", pseudot);
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
*/

#undef __FUNCT__
#define __FUNCT__ "MomentumSolver_DualTime_Picard_RK4"
/**
 * @brief Solves the Momentum Equations using Dual-Time Stepping with a Fixed-Point RK4 Smoother.
 *
 * =================================================================================================
 * GLOSSARY & THEORETICAL BASIS
 * =================================================================================================
 * 
 * 1. METHODOLOGY: Dual-Time Stepping (Pseudo-Time Integration)
 *    We aim to solve the implicit BDF equation:  R_spatial(U) + dU/dt_physical = 0.
 *    We do this by introducing a fictitious "Pseudo-Time" (tau) and iterating to steady state:
 *       dU/d(tau) = - [ R_spatial(U) + BDF_Terms(U) ]
 *    When dU/d(tau) -> 0, the physical time step is satisfied.
 *
 * 2. ALGORITHM: Fixed-Point Iteration with Explicit Runge-Kutta
 *    This is technically a Fixed-Point iteration on the operator:
 *       U_new = U_old + pseudo_dt_scaling * dt_pseudo * Total_Residual(U_old)
 *    We use a 4-Stage Explicit RK scheme (Jameson-Schmidt-Turkel coeffs) to smooth errors.
 *
 * 3. STABILITY: Backtracking Line Search
 *    If a pseudo-time step causes the Residual or Solution Error to GROW (Divergence),
 *    the solver "Backtracks": it restores the previous solution, cuts the pseudo-time step
 *    scaling factor (pseudo_dt_scaling) in half, and retries the iteration.
 *
 * =================================================================================================
 * VARIABLE MAPPING
 * =================================================================================================
 * -- Physics Variables (Legacy Names Kept) --
 * ti   : Physical Time Step Index.
 * dt   : Physical Time Step size (Delta t).
 * st   : Pseudo-Time Step size (Delta tau).
 * alfa : Runge-Kutta stage coefficients {1/4, 1/3, 1/2, 1}.
 *
 * -- Convergence & Solver Control (Renamed) --
 * pseudo_iter       : Counter for the inner dual-time loop.
 * pseudo_dt_scaling : Adaptive scalar for the pseudo-time step (formerly lambda).
 * delta_sol_norm    : The L_inf norm of the change in solution (dU).
 * resid_norm        : The L_inf norm of the Total Residual (RHS).
 * =================================================================================================
 */
PetscErrorCode MomentumSolver_DualTime_Picard_RK4(UserCtx *user, IBMNodes *ibm, FSInfo *fsi)
{
    // --- CONTEXT ACQUISITION BLOCK ---
    SimCtx *simCtx = user[0].simCtx;
    const PetscInt block_number = simCtx->block_number;
    
    // Legacy Names (Physics/Time) - Kept for brevity in formulas
    const PetscInt ti = simCtx->step;
    const PetscReal dt = simCtx->dt;       // Physical Time Step
    const PetscReal st = simCtx->st;       // Pseudo Time Step
    const PetscReal cfl = simCtx->cfl;     // Initial CFL Scaling
    const PetscReal alfa[] = {0.25, 1.0/3.0, 0.5, 1.0}; // RK4 Coefficients

    // Renamed Solver Parameters
    const PetscInt  max_pseudo_steps = simCtx->mom_max_pseudo_steps; // Max Dual-Time Iterations
    const PetscReal tol_abs_delta    = simCtx->mom_atol; // Stop if |dU| < tol
    const PetscReal tol_rtol_delta   = simCtx->mom_rtol; // Stop if |dU|/|dU0| < tol
    // --- END CONTEXT ACQUISITION BLOCK ---

    PetscErrorCode ierr;
    PetscMPIInt    rank;
    PetscInt       istage, pseudo_iter, ibi;
    PetscReal      ts, te, cput;

    // --- Global Convergence Metrics ---
    PetscReal global_norm_delta = 10.0; 
    PetscReal global_rel_delta  = 1.0;  
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

    ierr = PetscTime(&ts); CHKERRQ(ierr);

    // --- 1. Pre-Loop Initialization ---
    
    if (block_number > 1) {
      // ierr = Block_Interface_U(user); CHKERRQ(ierr);
    }

    for (PetscInt bi = 0; bi < block_number; bi++) {
        ierr = ApplyBoundaryConditions(&user[bi]); CHKERRQ(ierr);
        
        // Allocate workspace
        ierr = VecDuplicate(user[bi].Ucont, &user[bi].Rhs); CHKERRQ(ierr);
        ierr = VecDuplicate(user[bi].Ucont, &user[bi].dUcont); CHKERRQ(ierr);
        ierr = VecDuplicate(user[bi].Ucont, &user[bi].pUcont); CHKERRQ(ierr);
        
        // Initialize Backup (pUcont) with current state
        ierr = VecCopy(user[bi].Ucont, user[bi].pUcont); CHKERRQ(ierr);

        // --- Calculate Initial Total Residual (Spatial + Temporal) ---
        ierr = ComputeTotalResidual(&user[bi]); CHKERRQ(ierr);

        // Compute Initial Norms
        ierr = VecNorm(user[bi].Rhs, NORM_INFINITY, &resid_norm_init[bi]); CHKERRQ(ierr);
        
        // Initialize history for backtracking logic
        resid_norm_prev[bi] = 1000.0; 
        delta_sol_norm_prev[bi] = 1000.0;
        pseudo_dt_scaling[bi] = cfl; // Initialize adaptive scalar

	    LOG_ALLOW(GLOBAL,LOG_INFO," Block %d | Max RHS = %.6f | CFL = %.4f .\n", bi, resid_norm_init[bi], cfl);
    }
    
    // --- 2. Main Pseudo-Time Iteration Loop ---
    pseudo_iter = 0;
    while (((global_norm_delta > tol_abs_delta && global_rel_delta > tol_rtol_delta) || pseudo_iter < 1) 
            && pseudo_iter < max_pseudo_steps) 
    {
        pseudo_iter++;
	
        for (PetscInt bi = 0; bi < block_number; bi++) {
            
            // === 4-Stage Explicit Runge-Kutta Loop ===
            for (istage = 0; istage < 4; istage++) {
                
                LOG_ALLOW(GLOBAL,LOG_TRACE," Pseudo-Iter: %d | RK-Stage: %d\n", pseudo_iter, istage);
                
                // RK Update: U_new = U_old + (Scaler * Alpha * dt * st) * Total_Residual
                ierr = VecWAXPY(user[bi].Ucont, 
                                pseudo_dt_scaling[bi] * alfa[istage] * dt * st, 
                                user[bi].Rhs, 
                                user[bi].pUcont); CHKERRQ(ierr);

                // Sync Ghosts & Re-apply BCs for intermediate stage
                ierr = DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); CHKERRQ(ierr);
                ierr = DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); CHKERRQ(ierr);
                ierr = ApplyBoundaryConditions(&user[bi]); CHKERRQ(ierr);

                // --- Re-calculate Total Residual for next stage ---
                ierr = ComputeTotalResidual(&user[bi]); CHKERRQ(ierr);
      
            } // End RK Stages

            // === Convergence Metrics Calculation ===
            
            // Calculate dU = U_current - U_backup
            ierr = VecWAXPY(user[bi].dUcont, -1.0, user[bi].Ucont, user[bi].pUcont); CHKERRQ(ierr);

            // Shift history
            delta_sol_norm_prev[bi] = delta_sol_norm_curr[bi];
	    
            // Compute Infinity Norms
            ierr = VecNorm(user[bi].dUcont, NORM_INFINITY, &delta_sol_norm_curr[bi]); CHKERRQ(ierr);
            ierr = VecNorm(user[bi].Rhs, NORM_INFINITY, &resid_norm_curr[bi]); CHKERRQ(ierr);
	    
            // Normalize relative metrics
            if (pseudo_iter == 1) {
                delta_sol_norm_init[bi] = delta_sol_norm_curr[bi];
                delta_sol_rel_curr[bi] = 1.0;
	            resid_norm_init[bi] = resid_norm_curr[bi];
		        resid_rel_curr[bi] = 1.0;
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
                char filen[128];
                sprintf(filen, "%s/Momentum_Solver_Convergence_History_Block_%1d.log", simCtx->log_dir, bi);
                if(simCtx->step == simCtx->StartStep + 1 && pseudo_iter == 1) f = fopen(filen, "w");
                else f = fopen(filen, "a");
                
                PetscFPrintf(PETSC_COMM_WORLD, f, "Step: %d | PseudoIter: %d | |dU|: %le | Rel_dU: %le | |Resid|: %le \n",
                             (int)ti, (int)pseudo_iter, delta_sol_norm_curr[bi], delta_sol_rel_curr[bi], resid_norm_curr[bi]);
                fclose(f);
            }
        } // End loop over blocks

        // --- Update Global Convergence Criteria ---
        global_norm_delta = -1.0e20; 
        global_rel_delta  = -1.0e20; 
        global_rel_resid  = -1.0e20;
        
        for (PetscInt bi = 0; bi < block_number; bi++) {
            global_norm_delta = PetscMax(delta_sol_norm_curr[bi], global_norm_delta);
            global_rel_delta  = PetscMax(delta_sol_rel_curr[bi],  global_rel_delta);
            global_rel_resid  = PetscMax(resid_rel_curr[bi],      global_rel_resid);
        }

        ierr = PetscTime(&te); CHKERRQ(ierr);
        cput = te - ts;
        LOG_ALLOW(GLOBAL, LOG_INFO, "  Pseudo-Iter %d: |dU|=%e, Rel_dU=%e, |Resid|=%e, CPU=%.2fs\n",
                  pseudo_iter, global_norm_delta, global_rel_delta, global_rel_resid, cput);

        // === Backtracking Line Search ===
        for (PetscInt bi = 0; bi < block_number; bi++) {
            
            // TRIGGER: Divergence Detected
            if (pseudo_iter > 1 && 
                resid_norm_curr[bi] > resid_norm_prev[bi] && 
                delta_sol_norm_curr[bi] > delta_sol_norm_prev[bi] && 
                pseudo_dt_scaling[bi] > 0.01) 
            {
                // 1. Reduce Step Size
                pseudo_dt_scaling[bi] *= 0.5;
                // 2. Decrement iterator to retry
                pseudo_iter--; 
                
                LOG_ALLOW(LOCAL, LOG_DEBUG, "  Block %d: Backtracking triggered. Reducing scaling to %f.\n", bi, pseudo_dt_scaling[bi]);
                
                // 3. RESTORE from Backup
                ierr = VecCopy(user[bi].pUcont, user[bi].Ucont); CHKERRQ(ierr);
                
                // 4. SYNC Ghost Cells
                ierr = DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); CHKERRQ(ierr);
                ierr = DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); CHKERRQ(ierr);
            } 
            else {
                // SUCCESS: Step accepted. Update Backup and History.
                ierr = VecCopy(user[bi].Ucont, user[bi].pUcont); CHKERRQ(ierr);
                resid_norm_prev[bi]     = resid_norm_curr[bi];
                delta_sol_norm_prev[bi] = delta_sol_norm_curr[bi];
            }
        }
        
        if (block_number > 1) {
             // ierr = Block_Interface_U(user); CHKERRQ(ierr);
        }
    } // End while loop

    // --- Final Cleanup ---
    for (PetscInt bi = 0; bi < block_number; bi++) {
        ierr = VecDestroy(&user[bi].Rhs); CHKERRQ(ierr);
        ierr = VecDestroy(&user[bi].dUcont); CHKERRQ(ierr);
        ierr = VecDestroy(&user[bi].pUcont); CHKERRQ(ierr);
    }
    
    ierr = PetscFree2(delta_sol_norm_init, delta_sol_norm_prev);CHKERRQ(ierr);
    ierr = PetscFree2(delta_sol_norm_curr, delta_sol_rel_curr);CHKERRQ(ierr);
    ierr = PetscFree2(resid_norm_init,     resid_norm_prev);CHKERRQ(ierr);
    ierr = PetscFree2(resid_norm_curr,     pseudo_dt_scaling); CHKERRQ(ierr);
    ierr = PetscFree(resid_rel_curr); CHKERRQ(ierr);
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Momentum Solver converged/completed after %d pseudo-iterations.\n", pseudo_iter);
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}