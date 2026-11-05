/**
 * @file initialcondition.c  //  Setup the Initial conditions for different cases. 
 * @brief Test program for DMSwarm interpolation using the fdf-curvIB method.
 **/

 #include "initialcondition.h"
 #include "expression.h"
 #include "Metric.h"
 #include "poisson.h"

/**
 * @brief Derive Ucont from a Ucat whose interior cells are set: sync the periodic and ghost
 *        cells, then interpolate to the faces.
 */
static PetscErrorCode UcontFromInteriorUcat(UserCtx *user)
{
    const FieldId cell_fields[] = {FIELD_ID_UCAT};

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCall(SynchronizePeriodicCellFields(user, 1, cell_fields));
    PetscCall(UpdateLocalGhosts(user, FIELD_ID_UCAT));
    PetscCall(Cart2Contra(user));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Read one staged initial field (`ufield00000_0` or `vfield00000_0`) from the
 *        initial-condition directory into @p target.
 */
static PetscErrorCode ReadStagedInitialField(UserCtx *user, const char *name, Vec target)
{
    SimCtx *simCtx = user->simCtx;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(PetscStrncpy(simCtx->_io_context_buffer, simCtx->initialConditionDirectory,
                           sizeof(simCtx->_io_context_buffer)));
    simCtx->current_io_directory = simCtx->_io_context_buffer;
    PetscCall(ReadFieldData(user, name, target, "dat"));
    simCtx->current_io_directory = NULL;
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PopulateInitialUcont().
 * @see PopulateInitialUcont()
 */
PetscErrorCode PopulateInitialUcont(UserCtx *user)
{
    SimCtx *simCtx = user->simCtx;
    const PetscBool staged_flux = (PetscBool)(simCtx->initialConditionStaged && simCtx->initialConditionField == IC_FIELD_UCONT);
    PicurvVelocitySource source;
    PetscReal physical_time, weight = 1.0, factor = 1.0;
    Vec fraction, length;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCall(PicurvVelocitySourceCreate("ic_source", PETSC_FALSE, simCtx->initialConditionStaged, &source));
    PetscCall(PicurvPhysicalTime(simCtx, simCtx->ti, &physical_time));
    PetscCall(VecSet(user->Ucat, 0.0));
    PetscCall(VecSet(user->Ucont, 0.0));
    if (simCtx->initialConditionStaged) {
        PetscCall(PetscOptionsGetReal(NULL, NULL, "-ic_source_field_weight", &weight, NULL));
        PetscCall(ReadStagedInitialField(user, staged_flux ? "vfield00000_0" : "ufield00000_0", staged_flux ? user->Ucont : user->Ucat));
        if (weight != 1.0) PetscCall(VecScale(staged_flux ? user->Ucont : user->Ucat, weight));
        if (!staged_flux) PetscCall(UcontFromInteriorUcat(user));
    }
    if (source.count > 0) {
        PetscCall(DMCreateGlobalVector(user->fda, &fraction));
        PetscCall(DMCreateGlobalVector(user->fda, &length));
        PetscCall(ComputeGridLineFractions(user, fraction, length));
        if (source.bulk) PetscCall(PicurvSourceCellBulkFactor(user, &source, fraction, length, physical_time, &factor));
        PetscCall(PicurvSourceAddFluxes(user, &source, fraction, length, physical_time, factor, PETSC_FALSE));
        PetscCall(VecDestroy(&fraction));
        PetscCall(VecDestroy(&length));
    }
    PetscCall(PicurvVelocitySourceDestroy(&source));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Complete one initial-condition block after its interior values are assigned.
 */
static PetscErrorCode FinalizeBlockState(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    // This sequence ensures a fully consistent state for a single block.
    ierr = ApplyBoundaryConditions(user); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_TRACE," Boundary condition applied.\n");
    // 2. Sync contravariant velocity field.
    const FieldId staggered_fields[] = {FIELD_ID_UCONT};
    ierr = SynchronizePeriodicStaggeredFields(user, 1, staggered_fields); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_TRACE," Ucont field ghosts updated.\n");
    
    // 3. Convert to Cartesian velocity.
    ierr = Contra2Cart(user); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_TRACE," Converted Ucont to Ucat.\n");

    // 4. Finalize periodic endpoint values, then refresh local Cartesian velocity.
    {
        const FieldId cell_fields[] = {FIELD_ID_UCAT};
        ierr = SynchronizePeriodicCellFields(user, 1, cell_fields); CHKERRQ(ierr);
    }
    ierr = UpdateLocalGhosts(user, FIELD_ID_UCAT); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_TRACE," Ucat field ghosts updated.\n"); 

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}


/**
 * @brief Make the fresh-start velocity divergence-free with the time step's own projection:
 *        solve the pressure Poisson equation for the initial Ucont, whose boundary fluxes are
 *        already set, and correct the interior faces. As after a step's projection, the
 *        boundary handlers are not reapplied, so the field stays divergence-free. The
 *        potential is discarded, so pressure keeps its initial value. Collective.
 */
static PetscErrorCode ProjectInitialVelocity(SimCtx *simCtx)
{
    UserCtx *user_finest = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;
    const FieldId potential[] = {FIELD_ID_PHI};

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCall(PoissonSolver_Multigrid(&simCtx->usermg));
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        /* UpdatePressure's refresh of Phi, without adding it to the pressure. */
        PetscCall(SynchronizePeriodicCellFields(&user_finest[bi], 1, potential));
        PetscCall(UpdateLocalGhosts(&user_finest[bi], FIELD_ID_PHI));
        PetscCall(ProjectVelocity(&user_finest[bi]));
        PetscCall(VecSet(user_finest[bi].Phi, 0.0));
        PetscCall(UpdateLocalGhosts(&user_finest[bi], FIELD_ID_PHI));
    }
    LOG(GLOBAL, LOG_INFO, "Initial velocity projected to be divergence-free (Poisson source imbalance %g).\n",
              (double)simCtx->poissonSourceImbalance);
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Initialize Eulerian fields for a new simulation without restart data.
 */
static PetscErrorCode SetInitialFluidState_FreshStart(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserCtx *user_finest = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;
    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d, Block %" PetscInt_FMT ": Setting t=0 state.\n", simCtx->rank, bi);

        // 1. Set an initial guess for the INTERIOR of the domain.
        //    Replaces the legacy `if(InitialGuessOne)` block.

	LOG_ALLOW(GLOBAL,LOG_TRACE," Initializing Interior Ucont field.\n");
        ierr = PopulateInitialUcont(&user_finest[bi]); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_TRACE," Interior Ucont field initialized.\n");

    // 2. Apply all boundary conditions, convert to Cartesian, and sync ghosts.
    LOG_ALLOW(GLOBAL,LOG_TRACE," Boundary condition application and state finalization initiated.\n");
    ierr = FinalizeBlockState(&user_finest[bi]); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_TRACE," Boundary condition application and state finalization complete.\n");
    }

    // If using multiple grid blocks, handle the interface conditions between them.
    if (simCtx->block_number > 1) {
      //   LOG_ALLOW(GLOBAL, LOG_INFO, "Updating multi-block interfaces for t=0.\n");
	//	ierr = Block_Interface_U(user_finest); CHKERRQ(ierr);
        // After interface update, ghost regions might be stale. Refresh them.
        for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
             const FieldId staggered_fields[] = {FIELD_ID_UCONT};
             ierr = SynchronizePeriodicStaggeredFields(&user_finest[bi], 1, staggered_fields); CHKERRQ(ierr);
             ierr = UpdateLocalGhosts(&user_finest[bi], FIELD_ID_UCAT); CHKERRQ(ierr);
        }
    }

    {
        PetscBool project = PETSC_FALSE;
        ierr = PetscOptionsGetBool(NULL, NULL, "-ic_source_project", &project, NULL); CHKERRQ(ierr);
        if (project) { ierr = ProjectInitialVelocity(simCtx); CHKERRQ(ierr); }
    }

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

/**
 * @brief Restore Eulerian fields from checkpoint files for a restart simulation.
 */
static PetscErrorCode SetInitialFluidState_Load(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserCtx *user_finest = simCtx->usermg.mgctx[simCtx->mglevels - 1].user;
    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d, Block %" PetscInt_FMT ": Reading restart files for step %" PetscInt_FMT ".\n",
                  simCtx->rank, bi, simCtx->StartStep);

        // ReadSimulationFields handles all file I/O for one block.
        ierr = ReadSimulationFields(&user_finest[bi], simCtx->StartStep); CHKERRQ(ierr);
        
        // Apply Boundary Conditions on Read fields
        ierr = ApplyBoundaryConditions(&user_finest[bi]); CHKERRQ(ierr);
        // After reading from a file, the local ghost regions MUST be updated
        // to ensure consistency across process boundaries for the first time step.
        //ierr = UpdateLocalGhosts(&user_finest[bi], FIELD_ID_UCAT); CHKERRQ(ierr);
        //ierr = UpdateLocalGhosts(&user_finest[bi], FIELD_ID_P); CHKERRQ(ierr);
        // ... add ghost updates for any other fields read from file ...
    }
    
    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `InitializeEulerianState()`.
 * @details Local to this translation unit.
 */
PetscErrorCode InitializeEulerianState(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserCtx *user_finest = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_INFO, "--- Initializing Eulerian State ---\n");

    if (simCtx->StartStep > 0) {
        if(strcmp(simCtx->eulerianSource,"analytical")==0){
            LOG_ALLOW(GLOBAL,LOG_INFO,"Initializing analytical field from its source (t=%.4f, step=%" PetscInt_FMT ").\n",simCtx->StartTime,simCtx->StartStep);
            ierr = AnalyticalSolutionEngine(simCtx);
        }
        else{
            LOG_ALLOW(GLOBAL, LOG_INFO, "Starting from RESTART files (t=%.4f, step=%" PetscInt_FMT ").\n",
                    simCtx->StartTime, simCtx->StartStep);
            ierr = SetInitialFluidState_Load(simCtx); CHKERRQ(ierr);
        }
        /* Statistics resume from the same bundle the flow state came from, and are
         * restored regardless of the Eulerian source: an analytical restart still
         * continues a window that was accumulating before it. */
        ierr = RestoreFieldStatisticsState(simCtx, simCtx->StartStep); CHKERRQ(ierr);
    } else { // StartStep = 0
        LOG_ALLOW(GLOBAL, LOG_INFO, "Performing a FRESH START (t=0, step=0).\n");
        if(strcmp(simCtx->eulerianSource,"solve")==0){
            ierr = SetInitialFluidState_FreshStart(simCtx); CHKERRQ(ierr);
        }else if(strcmp(simCtx->eulerianSource,"load")==0){
            LOG_ALLOW(GLOBAL,LOG_INFO,"FRESH START in LOAD mode. Reading files (t=%.4f,step=%" PetscInt_FMT ").\n",
                      simCtx->StartTime,simCtx->StartStep);
            ierr=SetInitialFluidState_Load(simCtx);CHKERRQ(ierr);
        }else if(strcmp(simCtx->eulerianSource,"analytical")==0){
            LOG_ALLOW(GLOBAL,LOG_INFO,"FRESH START in ANALYTICAL mode (t=%.4f,step=%" PetscInt_FMT ").\n",
                      simCtx->StartTime,simCtx->StartStep);
            ierr=AnalyticalSolutionEngine(simCtx);CHKERRQ(ierr);
        }
    }

    // This crucial step, taken from the end of the legacy setup, ensures
    // that the history vectors (Ucont_o, Ucont_rm1, etc.) are correctly
    // populated before the first call to the time-stepping loop.
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        ierr = UpdateSolverHistoryVectors(&user_finest[bi],
                                           (PetscBool)(simCtx->StartStep > 0 && simCtx->restartHistoryAvailable)); CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "--- Eulerian State Initialized and History Vectors Populated ---\n");

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
