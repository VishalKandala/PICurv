#include "implicitsolvers.h"

#undef __FUNCT__
#define __FUNCT__ "RungeKutta"
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
 * @param user The array of UserCtx structs for all blocks.
 * @param ibm  (Optional) Pointer to the full array of IBM data structures. Pass NULL if disabled.
 * @param fsi  (Optional) Pointer to the full array of FSI data structures. Pass NULL if disabled.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode RungeKutta(UserCtx *user, IBMNodes *ibm, FSInfo *fsi)
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
        ierr = InflowFlux(&user[bi]); CHKERRQ(ierr);
        ierr = OutflowFlux(&user[bi]); CHKERRQ(ierr);
        ierr = FormBCS(&user[bi]); CHKERRQ(ierr);

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
                ierr = InflowFlux(&user[bi]); CHKERRQ(ierr);
                ierr = OutflowFlux(&user[bi]); CHKERRQ(ierr);
                ierr = FormBCS(&user[bi]); CHKERRQ(ierr);
                
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
 */
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
        ierr = InflowFlux(&user[bi]); CHKERRQ(ierr);
        ierr = OutflowFlux(&user[bi]); CHKERRQ(ierr);
        ierr = FormBCS(&user[bi]); CHKERRQ(ierr);
        
        /*
        if (immersed) {
            for (ibi = 0; ibi < simCtx->NumberOfBodies; ibi++) {
                ierr = ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi, 1); CHKERRQ(ierr);
            }
        }
        */
        
        ierr = VecDuplicate(user[bi].Ucont, &user[bi].Rhs); CHKERRQ(ierr);
        ierr = VecDuplicate(user[bi].Ucont, &user[bi].dUcont); CHKERRQ(ierr);
        ierr = VecDuplicate(user[bi].Ucont, &user[bi].pUcont); CHKERRQ(ierr);
        ierr = VecCopy(user[bi].Ucont, user[bi].pUcont); CHKERRQ(ierr);


	LOG_ALLOW(GLOBAL,LOG_DEBUG," BCs and workspace vectors prepared for Initial RHS calculation.\n");
	
        // --- Calculate INITIAL RHS for convergence check ---
        ierr = ComputeRHS(&user[bi], user[bi].Rhs); CHKERRQ(ierr);

	LOG_ALLOW(GLOBAL,LOG_DEBUG, " Initial RHS calculated for convergence check .\n");
	
        // Add time derivative part
        if (COEF_TIME_ACCURACY > 1.1 && ti != tistart && ti != 1) {
            ierr = VecAXPY(user[bi].Rhs, -COEF_TIME_ACCURACY/dt, user[bi].Ucont); CHKERRQ(ierr);
            ierr = VecAXPY(user[bi].Rhs, +2.0/dt, user[bi].Ucont_o); CHKERRQ(ierr);
            ierr = VecAXPY(user[bi].Rhs, -0.5/dt, user[bi].Ucont_rm1); CHKERRQ(ierr);
        } else {
            ierr = VecAXPY(user[bi].Rhs, -1.0/dt, user[bi].Ucont); CHKERRQ(ierr);
            ierr = VecAXPY(user[bi].Rhs, +1.0/dt, user[bi].Ucont_o); CHKERRQ(ierr);
        }

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
	      LOG_ALLOW(GLOBAL,LOG_INFO," Pseudo-Timestep Solver | RK-Stage : %d | Pseudo-Timestep :%d .\n",istage,pseudot);
                ierr = VecWAXPY(user[bi].Ucont, lambda[bi] * alfa[istage] * dt * st, user[bi].Rhs, user[bi].pUcont); CHKERRQ(ierr);

		LOG_ALLOW(GLOBAL,LOG_INFO, " Ucont updated as Ucont = [lambda(%.4f)]x[alfa(%.4f)]x[dt(%.4f)]x[st(%.4f)]xRhs + pUcont.\n",lambda[bi],alfa[istage],dt,st);
		
                // Sync ghosts and re-apply BCs for the intermediate stage
                ierr = DMGlobalToLocalBegin(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); CHKERRQ(ierr);
                ierr = DMGlobalToLocalEnd(user[bi].fda, user[bi].Ucont, INSERT_VALUES, user[bi].lUcont); CHKERRQ(ierr);
                ierr = InflowFlux(&user[bi]); CHKERRQ(ierr);
                ierr = OutflowFlux(&user[bi]); CHKERRQ(ierr);
                ierr = FormBCS(&user[bi]); CHKERRQ(ierr);

		LOG_ALLOW(GLOBAL,LOG_INFO, " Ghosts synced and BCs reapplied.\n");

                // --- Re-calculate the full RHS for the next RK stage ---
                ierr = ComputeRHS(&user[bi], user[bi].Rhs); CHKERRQ(ierr);

		LOG_ALLOW(GLOBAL,LOG_INFO, " RHS calculated for the next RK stage.\n");
		
                if (COEF_TIME_ACCURACY > 1.1 && ti != tistart && ti != 1) {
                    ierr = VecAXPY(user[bi].Rhs, -COEF_TIME_ACCURACY/dt, user[bi].Ucont); CHKERRQ(ierr);
                    ierr = VecAXPY(user[bi].Rhs, +2.0/dt, user[bi].Ucont_o); CHKERRQ(ierr);
                    ierr = VecAXPY(user[bi].Rhs, -0.5/dt, user[bi].Ucont_rm1); CHKERRQ(ierr);
                } else {
                    ierr = VecAXPY(user[bi].Rhs, -1.0/dt, user[bi].Ucont); CHKERRQ(ierr);
                    ierr = VecAXPY(user[bi].Rhs, +1.0/dt, user[bi].Ucont_o); CHKERRQ(ierr);
                }

		if(istage <3) LOG_ALLOW(GLOBAL,LOG_INFO," Time derivative part updated for next RK stage .\n");
		else LOG_ALLOW(GLOBAL,LOG_INFO," All RK stages complete.\n");
            } // End RK stages

            /*
            if (immersed) {
                for (ibi = 0; ibi < simCtx->NumberOfBodies; ibi++) {
                    ierr = ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi, 1); CHKERRQ(ierr);
                }
            }
            */

            // --- Calculate Convergence Metrics for this block ---
	    
            ierr = VecWAXPY(user[bi].dUcont, -1.0, user[bi].Ucont, user[bi].pUcont); CHKERRQ(ierr);

	    LOG_ALLOW(GLOBAL,LOG_DEBUG," Block %d | dU calculated .\n",bi);
            
            normdU1_bk[bi] = normdU_bk[bi];
	    
            ierr = VecNorm(user[bi].dUcont, NORM_INFINITY, &normdU_bk[bi]); CHKERRQ(ierr);
            ierr = VecNorm(user[bi].Rhs, NORM_INFINITY, &normF_bk[bi]); CHKERRQ(ierr);
	    
	    LOG_ALLOW(GLOBAL,LOG_DEBUG," Block %d | |dU|  =  %le | |U| = %le.\n",bi,normdU_bk[bi],normF_bk[bi]);
	    
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
	      

	    LOG_ALLOW(GLOBAL,LOG_DEBUG, " Block %d | |dU|/|dU_0| = %le. \n",bi,reldU_bk[bi]);

	    LOG_ALLOW(GLOBAL,LOG_DEBUG, " Block %d | |R|/|R_0| = %le. \n",bi,relF_bk[bi]);
	    
            // file logging
            if (!rank) {
                FILE *f;
                char filen[80];
                sprintf(filen, "%s/Momentum_Solver_Convergence_History_Block_%1d.log", simCtx->log_dir,bi);
                f = fopen(filen, "a");
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
    /*
    if (block_number > 1) {
        Block_Interface_Correction(user);
        if (simCtx->blank) {
            Block_Blank_Correction_adv(&user[0], 0);
        }
    }
    */
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
