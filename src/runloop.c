/**
 * @file runloop.c
 * @brief Test program for DMSwarm interpolation using the fdf-curvIB method.
 * Provides the setup to start any simulation with DMSwarm and DMDAs.
 **/

#include <ctype.h>
#include <errno.h>
#include <signal.h>
#include "runloop.h"

static volatile sig_atomic_t g_runtime_shutdown_signal = 0;
static PetscBool             g_runtime_shutdown_auto_requested = PETSC_FALSE;

/**
 * @brief Internal helper implementation: `RuntimeShutdownSignalHandler()`.
 * @details Local to this translation unit.
 */
static void RuntimeShutdownSignalHandler(int signum)
{
    if (g_runtime_shutdown_signal == 0) {
        g_runtime_shutdown_signal = signum;
    }
}

/**
 * @brief Internal helper implementation: `RuntimeRequestAutoWalltimeGuard()`.
 * @details Local to this translation unit.
 */
static void RuntimeRequestAutoWalltimeGuard(void)
{
    if (g_runtime_shutdown_signal == 0) {
        g_runtime_shutdown_auto_requested = PETSC_TRUE;
    }
}

/**
 * @brief Internal helper implementation: `RuntimeShutdownRequested()`.
 * @details Local to this translation unit.
 */
static PetscBool RuntimeShutdownRequested(void)
{
    return (PetscBool)(g_runtime_shutdown_signal != 0 || g_runtime_shutdown_auto_requested);
}

/**
 * @brief Internal helper implementation: `RuntimeShutdownSignalName()`.
 * @details Local to this translation unit.
 */
static const char *RuntimeShutdownSignalName(PetscInt signum)
{
    switch (signum) {
#ifdef SIGTERM
        case SIGTERM: return "SIGTERM";
#endif
#ifdef SIGUSR1
        case SIGUSR1: return "SIGUSR1";
#endif
#ifdef SIGINT
        case SIGINT: return "SIGINT";
#endif
        default: return "UNKNOWN";
    }
}

/**
 * @brief Internal helper implementation: `RuntimeShutdownReasonName()`.
 * @details Local to this translation unit.
 */
static const char *RuntimeShutdownReasonName(void)
{
    if (g_runtime_shutdown_signal != 0) {
        return RuntimeShutdownSignalName((PetscInt)g_runtime_shutdown_signal);
    }
    if (g_runtime_shutdown_auto_requested) {
        return "AUTO_WALLTIME_GUARD";
    }
    return "NONE";
}

/**
 * @brief Internal helper implementation: `RegisterRuntimeSignalHandler()`.
 * @details Local to this translation unit.
 */
static PetscErrorCode RegisterRuntimeSignalHandler(int signum)
{
    struct sigaction action;

    PetscFunctionBeginUser;
    memset(&action, 0, sizeof(action));
    action.sa_handler = RuntimeShutdownSignalHandler;

    if (sigemptyset(&action.sa_mask) != 0) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SYS, "sigemptyset failed for signal %d.", signum);
    }

    if (sigaction(signum, &action, NULL) != 0) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SYS, "sigaction failed for signal %d.", signum);
    }

    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref InitializeRuntimeSignalHandlers().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the matching public header declaration.
 * @see InitializeRuntimeSignalHandlers()
 */
PetscErrorCode InitializeRuntimeSignalHandlers(void)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    g_runtime_shutdown_signal = 0;
    g_runtime_shutdown_auto_requested = PETSC_FALSE;

#ifdef SIGTERM
    ierr = RegisterRuntimeSignalHandler(SIGTERM); CHKERRQ(ierr);
#endif
#ifdef SIGUSR1
    ierr = RegisterRuntimeSignalHandler(SIGUSR1); CHKERRQ(ierr);
#endif
#ifdef SIGINT
    ierr = RegisterRuntimeSignalHandler(SIGINT); CHKERRQ(ierr);
#endif

    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref RuntimeWalltimeGuardParsePositiveSeconds().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the matching public header declaration.
 * @see RuntimeWalltimeGuardParsePositiveSeconds()
 */
PetscBool RuntimeWalltimeGuardParsePositiveSeconds(const char *text, PetscReal *seconds_out)
{
    char   *endptr = NULL;
    double parsed_value;

    if (seconds_out) *seconds_out = 0.0;
    if (!text || text[0] == '\0') return PETSC_FALSE;

    errno        = 0;
    parsed_value = strtod(text, &endptr);
    if (endptr == text || errno == ERANGE || !isfinite(parsed_value) || parsed_value <= 0.0) {
        return PETSC_FALSE;
    }

    while (*endptr != '\0' && isspace((unsigned char)*endptr)) {
        endptr++;
    }
    if (*endptr != '\0') return PETSC_FALSE;

    if (seconds_out) *seconds_out = (PetscReal)parsed_value;
    return PETSC_TRUE;
}

/**
 * @brief Implementation of \ref RuntimeWalltimeGuardUpdateEWMA().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the matching public header declaration.
 * @see RuntimeWalltimeGuardUpdateEWMA()
 */
PetscReal RuntimeWalltimeGuardUpdateEWMA(PetscBool has_previous, PetscReal previous_ewma_seconds, PetscReal latest_step_seconds, PetscReal alpha)
{
    if (!has_previous) return latest_step_seconds;
    return alpha * latest_step_seconds + (1.0 - alpha) * previous_ewma_seconds;
}

/**
 * @brief Implementation of \ref RuntimeWalltimeGuardConservativeEstimate().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the matching public header declaration.
 * @see RuntimeWalltimeGuardConservativeEstimate()
 */
PetscReal RuntimeWalltimeGuardConservativeEstimate(PetscReal warmup_average_seconds, PetscReal ewma_seconds, PetscReal latest_step_seconds)
{
    return PetscMax(warmup_average_seconds, PetscMax(ewma_seconds, latest_step_seconds));
}

/**
 * @brief Implementation of \ref RuntimeWalltimeGuardRequiredHeadroom().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the matching public header declaration.
 * @see RuntimeWalltimeGuardRequiredHeadroom()
 */
PetscReal RuntimeWalltimeGuardRequiredHeadroom(PetscReal min_seconds, PetscReal multiplier, PetscReal conservative_estimate_seconds)
{
    return PetscMax(min_seconds, multiplier * conservative_estimate_seconds);
}

/**
 * @brief Implementation of \ref RuntimeWalltimeGuardShouldTrigger().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the matching public header declaration.
 * @see RuntimeWalltimeGuardShouldTrigger()
 */
PetscBool RuntimeWalltimeGuardShouldTrigger(PetscInt completed_steps, PetscInt warmup_steps, PetscReal remaining_seconds, PetscReal min_seconds, PetscReal multiplier, PetscReal warmup_average_seconds, PetscReal ewma_seconds, PetscReal latest_step_seconds, PetscReal *required_headroom_seconds_out)
{
    PetscReal conservative_estimate = 0.0;
    PetscReal required_headroom     = 0.0;

    if (required_headroom_seconds_out) *required_headroom_seconds_out = 0.0;
    if (completed_steps < warmup_steps) return PETSC_FALSE;

    conservative_estimate = RuntimeWalltimeGuardConservativeEstimate(warmup_average_seconds, ewma_seconds, latest_step_seconds);
    required_headroom     = RuntimeWalltimeGuardRequiredHeadroom(min_seconds, multiplier, conservative_estimate);
    if (required_headroom_seconds_out) *required_headroom_seconds_out = required_headroom;
    return (PetscBool)(remaining_seconds <= required_headroom);
}

/**
 * @brief Internal helper implementation: `RuntimeWalltimeGuardRemainingSeconds()`.
 * @details Local to this translation unit.
 */
static PetscReal RuntimeWalltimeGuardRemainingSeconds(const SimCtx *simCtx)
{
    time_t now = time(NULL);

    if (now == (time_t)-1) return -1.0;
    return simCtx->walltimeGuardLimitSeconds - ((PetscReal)now - simCtx->walltimeGuardJobStartEpochSeconds);
}

/**
 * @brief Internal helper implementation: `UpdateRuntimeWalltimeGuardEstimator()`.
 * @details Local to this translation unit.
 */
static void UpdateRuntimeWalltimeGuardEstimator(SimCtx *simCtx, PetscReal completed_step_seconds)
{
    simCtx->walltimeGuardCompletedSteps++;
    simCtx->walltimeGuardLatestStepSeconds = completed_step_seconds;

    if (simCtx->walltimeGuardCompletedSteps <= simCtx->walltimeGuardWarmupSteps) {
        simCtx->walltimeGuardWarmupTotalSeconds += completed_step_seconds;
        if (simCtx->walltimeGuardCompletedSteps == simCtx->walltimeGuardWarmupSteps) {
            simCtx->walltimeGuardWarmupAverageSeconds = simCtx->walltimeGuardWarmupTotalSeconds / (PetscReal)simCtx->walltimeGuardWarmupSteps;
            simCtx->walltimeGuardEWMASeconds          = simCtx->walltimeGuardWarmupAverageSeconds;
            simCtx->walltimeGuardHasEWMA             = PETSC_TRUE;
        }
        return;
    }

    simCtx->walltimeGuardEWMASeconds = RuntimeWalltimeGuardUpdateEWMA(
        simCtx->walltimeGuardHasEWMA,
        simCtx->walltimeGuardEWMASeconds,
        completed_step_seconds,
        simCtx->walltimeGuardEstimatorAlpha
    );
    simCtx->walltimeGuardHasEWMA = PETSC_TRUE;
}

/**
 * @brief Internal helper implementation: `MaybeRequestRuntimeWalltimeGuardShutdown()`.
 * @details Local to this translation unit.
 */
static PetscErrorCode MaybeRequestRuntimeWalltimeGuardShutdown(SimCtx *simCtx, const char *checkpoint_name)
{
    PetscReal remaining_seconds = 0.0;
    PetscReal required_headroom = 0.0;

    PetscFunctionBeginUser;
    if (!simCtx->walltimeGuardActive || RuntimeShutdownRequested()) PetscFunctionReturn(0);

    remaining_seconds = RuntimeWalltimeGuardRemainingSeconds(simCtx);
    if (!RuntimeWalltimeGuardShouldTrigger(
            simCtx->walltimeGuardCompletedSteps,
            simCtx->walltimeGuardWarmupSteps,
            remaining_seconds,
            simCtx->walltimeGuardMinSeconds,
            simCtx->walltimeGuardMultiplier,
            simCtx->walltimeGuardWarmupAverageSeconds,
            simCtx->walltimeGuardEWMASeconds,
            simCtx->walltimeGuardLatestStepSeconds,
            &required_headroom)) {
        PetscFunctionReturn(0);
    }

    RuntimeRequestAutoWalltimeGuard();
    LOG_ALLOW(
        GLOBAL,
        LOG_WARNING,
        "[T=%.4f, Step=%d] AUTO_WALLTIME_GUARD requested at %s: remaining walltime %.1f s <= required headroom %.1f s "
        "(warmup avg %.1f s, ewma %.1f s, latest %.1f s).\n",
        simCtx->ti,
        simCtx->step,
        checkpoint_name,
        (double)remaining_seconds,
        (double)required_headroom,
        (double)simCtx->walltimeGuardWarmupAverageSeconds,
        (double)simCtx->walltimeGuardEWMASeconds,
        (double)simCtx->walltimeGuardLatestStepSeconds
    );

    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `WriteForcedTerminationOutput()`.
 * @details Local to this translation unit.
 */
static PetscErrorCode WriteForcedTerminationOutput(SimCtx *simCtx, UserCtx *user, const char *phase)
{
    PetscErrorCode ierr;

    PetscFunctionBeginUser;

    LOG(GLOBAL, LOG_WARNING,
        "[T=%.4f, Step=%d] Shutdown requested by %s during %s. Writing final output outside the normal cadence before exiting.\n",
        simCtx->ti, simCtx->step, RuntimeShutdownReasonName(), phase);

    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        ierr = WriteSimulationFields(&user[bi]); CHKERRQ(ierr);
    }

    if (simCtx->np > 0) {
        ierr = WriteAllSwarmFields(user); CHKERRQ(ierr);
    }

    if (IsParticleConsoleSnapshotEnabled(simCtx)) {
        ierr = EmitParticleConsoleSnapshot(user, simCtx, simCtx->step); CHKERRQ(ierr);
    }

    ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRMPI(ierr);
    fflush(stdout);

    PetscFunctionReturn(0);
}

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
    PetscReal      step_start_seconds = 0.0;
    PetscReal      step_elapsed_local = 0.0;
    PetscReal      step_elapsed_max = 0.0;
    // Get the master context from the first block. All blocks share it.
    UserCtx *user = simCtx->usermg.mgctx[simCtx->usermg.mglevels-1].user;

    // Retrieve control parameters from SimCtx for clarity
    const PetscInt  StartStep = simCtx->StartStep;
    const PetscInt  StepsToRun = simCtx->StepsToRun;
    const PetscReal dt = simCtx->dt;
    
    // Variables for particle removal statistics
    //PetscInt removed_local_ob, removed_global_ob;
    PetscInt removed_local_lost, removed_global_lost;
    PetscBool terminated_early = PETSC_FALSE;
    PetscInt  last_completed_loop_index = StartStep - 1;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting main time-marching loop: %d steps from step %d (t=%.4f), dt=%.4f\n",
              StepsToRun, StartStep, simCtx->StartTime, dt);

    // --- Main Time-Marching Loop ---
    for (PetscInt step = StartStep; step < StartStep + StepsToRun; step++) {
        ierr = MaybeRequestRuntimeWalltimeGuardShutdown(simCtx, "pre-step checkpoint"); CHKERRQ(ierr);
        if (RuntimeShutdownRequested()) {
            ierr = WriteForcedTerminationOutput(simCtx, user, "pre-step checkpoint"); CHKERRQ(ierr);
            terminated_early = PETSC_TRUE;
            break;
        }

        step_start_seconds = MPI_Wtime();
        
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

        last_completed_loop_index = step;

        step_elapsed_local = MPI_Wtime() - step_start_seconds;
        step_elapsed_max   = step_elapsed_local;
        ierr               = MPI_Allreduce(&step_elapsed_local, &step_elapsed_max, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD); CHKERRMPI(ierr);
        if (simCtx->walltimeGuardActive) {
            UpdateRuntimeWalltimeGuardEstimator(simCtx, step_elapsed_max);
            ierr = MaybeRequestRuntimeWalltimeGuardShutdown(simCtx, "post-step checkpoint"); CHKERRQ(ierr);
        }

        if (RuntimeShutdownRequested()) {
            ierr = WriteForcedTerminationOutput(simCtx, user, "post-step checkpoint"); CHKERRQ(ierr);
            terminated_early = PETSC_TRUE;
            break;
        }
    } // --- End of Time-Marching Loop ---

    // After the loop, print the final progress state on rank 0 and add a newline
    // to ensure subsequent terminal output starts on a fresh line.
    if (simCtx->rank == 0 && StepsToRun > 0) {
        if (!terminated_early && last_completed_loop_index >= StartStep) {
            PrintProgressBar(StartStep + StepsToRun - 1, StartStep, StepsToRun, simCtx->ti);
        } else if (terminated_early && last_completed_loop_index >= StartStep) {
            PrintProgressBar(last_completed_loop_index, StartStep, StepsToRun, simCtx->ti);
        }
        PetscPrintf(PETSC_COMM_SELF, "\n");
        fflush(stdout);
    }

    if (terminated_early) {
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "Time marching stopped early after %s. Final retained state is step %d at t=%.4f.\n",
                  RuntimeShutdownReasonName(), simCtx->step, simCtx->ti);
    } else {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Time marching completed. Final time t=%.4f.\n", simCtx->ti);
    }
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
