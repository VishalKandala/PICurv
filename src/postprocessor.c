/**
 * @file postprocess.c
 * @brief Phase 2 implementation of the post-processing tool.
 *
 * This phase introduces a dedicated configuration system and performs a
 * single-step data load to verify the I/O and data structures.
 */

#include "postprocess.h" // Use our new header

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    PetscErrorCode    ierr;
    SimCtx            *simCtx = NULL;
    PostProcessConfig config;
    char              ppConfigFile[PETSC_MAX_PATH_LEN] = "postprocess.cfg";

    // === I. INITIALIZE PETSC & MPI ===========================================
    ierr = PetscInitialize(&argc, &argv, (char *)0, "Unified Post-Processing Tool"); CHKERRQ(ierr);

    // === II. CONFIGURE SIMULATION & POST-PROCESSING CONTEXTS =================
    ierr = CreateSimulationContext(argc, argv, &simCtx); CHKERRQ(ierr);

    // Get the post-processing config file path from command line, default to "postprocess.cfg"
    ierr = PetscOptionsGetString(NULL, NULL, "-pp_config", ppConfigFile, sizeof(ppConfigFile), NULL); CHKERRQ(ierr);
    
    // Parse the dedicated post-processing configuration
    ierr = ParsePostProcessConfig(ppConfigFile, &config); CHKERRQ(ierr);

    // === III. SETUP GRID & DATA STRUCTURES ===================================
    ierr = SetupGridAndSolvers(simCtx); CHKERRQ(ierr);
    
    // Get the finest-level user context, as this is where we'll load data
    UserCtx *user = simCtx->usermg.mgctx[simCtx->usermg.mglevels-1].user;

    // === IV. PHASE 2: SINGLE-STEP DATA LOADING ===============================
    LOG_ALLOW(GLOBAL, LOG_INFO, "=============================================================\n");
    LOG_ALLOW(GLOBAL, LOG_INFO, "PHASE 2: Attempting to load data for single time step t=%d\n", config.startTime);

    PetscInt ti = config.startTime;

    // 1. Load all raw Eulerian data from files into the pre-allocated Vecs
    ierr = ReadSimulationFields(user, ti); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "-> ReadSimulationFields completed for t=%d.\n", ti);

    // 2. Distribute ghost cell information for the loaded fields. This is critical
    //    for any subsequent processing that requires stencil operations.
    ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "P"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "Nvert"); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "-> UpdateLocalGhosts completed for primary fields.\n", ti);

    // 3. (Verification) Print a vector to the screen to prove it was loaded
    LOG_ALLOW(GLOBAL, LOG_INFO, "-> Verification: Printing Pressure vector (P) to screen...\n");
    ierr = VecView(user->P, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "PHASE 2 VERIFICATION COMPLETE\n");
    LOG_ALLOW(GLOBAL, LOG_INFO, "=============================================================\n");

    // === V. FINALIZE =========================================================
    ierr = FinalizeSimulation(simCtx); CHKERRQ(ierr);
    ierr = PetscFinalize();
    return ierr;
}