/**
 * @file postprocessor.c
 * @brief Phase 2 implementation of the post-processing tool.
 *
 * This phase introduces a dedicated configuration system and performs a
 * single-step data load to verify the I/O and data structures.
 */

#include "postprocessor.h" // Use our new header

/**
 * @brief Parses the processing pipeline string and executes the requested kernels.
 * @param user The UserCtx containing the data to be transformed.
 * @param config The PostProcessConfig containing the pipeline string.
 * @return PetscErrorCode
 */
PetscErrorCode RunProcessingPipeline(UserCtx* user, PostProcessParams* pps)
{
    PetscErrorCode ierr;
    char *pipeline_copy, *token;

    PetscFunctionBeginUser;
    LOG_ALLOW(GLOBAL, LOG_INFO, "--- Starting Data Transformation Pipeline ---\n");

    // Make a writable copy of the pipeline string for strtok
    ierr = PetscStrallocpy(pps->process_pipeline, &pipeline_copy); CHKERRQ(ierr);

    token = strtok(pipeline_copy, ",");
    while (token) {
        TrimWhitespace(token); // Use the helper from postprocess_config.c
        LOG_ALLOW(GLOBAL, LOG_INFO, "Executing pipeline step: %s\n", token);

        if (strcmp(token, "ComputeQCriterion") == 0) {
            ierr = ComputeQCriterion(user); CHKERRQ(ierr);
        }
        // Add more kernels here
        // A more advanced parser could handle arguments like NodalAverage(P->P_nodal)
        else if (strcmp(token, "NodalAverage_P") == 0) { // Simple name for now
             ierr = ComputeNodalAverage(user, "P", "P_nodal"); CHKERRQ(ierr);
        }
        else {
             LOG_ALLOW(GLOBAL, LOG_WARNING, "Unknown processing kernel '%s' requested. Skipping.\n", token);
        }

        token = strtok(NULL, ",");
    }

    ierr = PetscFree(pipeline_copy); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "--- Data Transformation Pipeline Complete ---\n");
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    PetscErrorCode    ierr;
    SimCtx            *simCtx = NULL;

    // === I. INITIALIZE PETSC & MPI ===========================================
    ierr = PetscInitialize(&argc, &argv, (char *)0, "Unified Post-Processing Tool"); CHKERRQ(ierr);

    // === II. CONFIGURE SIMULATION & POST-PROCESSING CONTEXTS =================
    ierr = CreateSimulationContext(argc, argv, &simCtx); CHKERRQ(ierr);
    // === IIA. SET EXECUTION MODE (SOLVER vs POST-PROCESSOR) =====
    simCtx->exec_mode = EXEC_MODE_POSTPROCESSOR;
    // === III. SETUP GRID & DATA STRUCTURES ===================================
    ierr = SetupGridAndSolvers(simCtx); CHKERRQ(ierr);
    
    // Get the finest-level user context, as this is where we'll load data
    UserCtx *user = simCtx->usermg.mgctx[simCtx->usermg.mglevels-1].user;

// === IV. PHASE 3: SINGLE-STEP DATA PROCESSING ===========================
    LOG_ALLOW(GLOBAL, LOG_INFO, "=============================================================\n");
    LOG_ALLOW(GLOBAL, LOG_INFO, "PHASE 3: Loading and processing data for t=%d\n", simCtx->pps->startTime);

    PetscInt ti = simCtx->pps->startTime;

    // 1. Load Data
    ierr = ReadSimulationFields(user, ti); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "P"); CHKERRQ(ierr);
    
    // 2. Transform Data
    ierr = RunProcessingPipeline(user, simCtx->pps); CHKERRQ(ierr);

    // 3. Verification
    PetscReal Umax,Qcritmax;
    LOG_ALLOW(GLOBAL, LOG_INFO, "-> Verification: Max Q-Criterion vector \n");
    ierr = VecNorm(user->Qcrit, NORM_INFINITY, &Qcritmax); CHKERRQ(ierr);
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "-> Verification: Max Nodal Pressure vector \n");
    ierr = VecNorm(user->P_nodal, NORM_INFINITY, &Umax); CHKERRQ(ierr);    

    LOG_ALLOW(GLOBAL, LOG_INFO, "   Max |Qcrit| = %g \n", Qcritmax);
    LOG_ALLOW(GLOBAL, LOG_INFO, "   Max |P_nodal| = %g \n", Umax);
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "PHASE 3 VERIFICATION COMPLETE\n");    

    // === V. FINALIZE =========================================================
   // ierr = FinalizeSimulation(simCtx); CHKERRQ(ierr);
    ierr = PetscFinalize();
    return ierr;
}