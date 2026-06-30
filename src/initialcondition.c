/**
 * @file initialcondition.c  //  Setup the Initial conditions for different cases. 
 * @brief Test program for DMSwarm interpolation using the fdf-curvIB method.
 **/

 #include "initialcondition.h"

#undef __FUNCT__
#define __FUNCT__ "SetInitialInteriorField"
/**
 * @brief Internal helper implementation: `SetInitialInteriorField()`.
 * @details Local to this translation unit.
 */
PetscErrorCode SetInitialInteriorField(UserCtx *user, const char *fieldName)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    SimCtx *simCtx = user->simCtx;
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Setting initial INTERIOR field for '%s' with mode %d.\n", fieldName, simCtx->initialConditionMode);

    // This function currently only implements logic for Ucont.
    if (strcmp(fieldName, "Ucont") != 0) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Skipping SetInitialInteriorField for non-Ucont field '%s'.\n", fieldName);

        PROFILE_FUNCTION_END;

        PetscFunctionReturn(0);
    }

    // --- 1. Get DMDA info and grid dimensions ---
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(user->fda, &info); CHKERRQ(ierr);

    const PetscInt im_phys = info.mx - 1;
    const PetscInt jm_phys = info.my - 1;
    const PetscInt km_phys = info.mz - 1;

    const PetscReal u_cart = simCtx->InitialConstantContra.x;
    const PetscReal v_cart = simCtx->InitialConstantContra.y;
    const PetscReal w_cart = simCtx->InitialConstantContra.z;

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "IC cartesian=(%.3f,%.3f,%.3f) ic_velocity_physical=%.3f mode=%d\n",
              (double)u_cart, (double)v_cart, (double)w_cart,
              (double)simCtx->icVelocityPhysical, (int)simCtx->initialConditionMode);

    // --- 2. Early dispatch: cartesian Constant delegates to the uniform converter ---
    if (simCtx->initialConditionMode == IC_MODE_CONSTANT_CARTESIAN) {
        ierr = UniformCart2Contra(user, u_cart, v_cart, w_cart); CHKERRQ(ierr);
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    // --- 3. Resolve flow direction for streamwise Constant and Poiseuille ---
    const PetscBool needs_flow_dir = (PetscBool)(
        simCtx->initialConditionMode == IC_MODE_POISEUILLE ||
        simCtx->initialConditionMode == IC_MODE_CONSTANT_STREAMWISE);
    FlowDirection fd          = FLOW_DIR_UNSET;
    PetscInt      flow_axis   = 0;
    PetscReal     flow_dir_sign = 1.0;

    if (needs_flow_dir) {
        if (user->inletFaceDefined)
            fd = (FlowDirection)user->identifiedInletBCFace;
        else if (simCtx->flowDirection != FLOW_DIR_UNSET)
            fd = simCtx->flowDirection;
        else
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER,
                    "Streamwise Constant and Poiseuille IC modes require either an INLET face or -flow_direction.");
        flow_axis     = (PetscInt)fd / 2;
        flow_dir_sign = ((PetscInt)fd % 2 == 0) ? 1.0 : -1.0;
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "IC flow_direction=%d (axis=%d sign=%.1f)\n",
                  (int)fd, (int)flow_axis, (double)flow_dir_sign);
    }

    // --- 4. Open arrays for non-cartesian modes ---
    Cmpnts ***csi_arr, ***eta_arr, ***zet_arr;
    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi, &csi_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta, &eta_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet, &zet_arr); CHKERRQ(ierr);

    Cmpnts ***ucont_arr;
    ierr = DMDAVecGetArray(user->fda, user->Ucont, &ucont_arr); CHKERRQ(ierr);

    PetscInt i, j, k;
    const PetscInt xs = info.xs, xe = info.xs + info.xm;
    const PetscInt ys = info.ys, ye = info.ys + info.ym;
    const PetscInt zs = info.zs, ze = info.zs + info.zm;
        
    for (k = zs; k < ze; k++) {
        for (j = ys; j < ye; j++) {
            for (i = xs; i < xe; i++) {

                // Check to ensure we only set initial conditions for PHYSICAL cells, not ghost cells.
                // Ghost cells (at indices 0 and n) will be set later by ApplyBoundaryConditions.
                //
                // Grid structure: For n physical grid points, DMDA has size n+1
                //   - im_phys = mx - 1 = n (number of coordinate points, also equals number of cells + 1)
                //   - Physical cell indices: [1, im_phys-1] = [1, n-1] (gives n-1 physical cells)
                //   - Ghost cells at boundaries: index 0 and index im_phys (= n)
                //
                // Example: n=25 physical points → im_phys=25
                //   - Physical cells: indices 1..24 (24 cells)
                //   - Ghost cells: indices 0 and 25
                const PetscBool is_interior = (i > 0 && i < im_phys &&
                                               j > 0 && j < jm_phys &&
                                               k > 0 && k < km_phys);

                if (is_interior) {
                    Cmpnts ucont_val = {0.0, 0.0, 0.0}; // Default to zero velocity
                    PetscReal normal_velocity_mag = 0.0;

                    switch (simCtx->initialConditionMode) {
                        case IC_MODE_ZERO:
                            break;
                        case IC_MODE_CONSTANT_STREAMWISE:
                            normal_velocity_mag = simCtx->icVelocityPhysical;
                            break;
                        case IC_MODE_POISEUILLE:
                            {
                                PetscInt cs1, cs2, n1, n2;
                                if (flow_axis == 0)      { cs1 = j; cs2 = k; n1 = jm_phys; n2 = km_phys; }
                                else if (flow_axis == 1) { cs1 = i; cs2 = k; n1 = im_phys; n2 = km_phys; }
                                else                     { cs1 = i; cs2 = j; n1 = im_phys; n2 = jm_phys; }
                                const PetscReal w1 = (PetscReal)(n1 - 2);
                                const PetscReal w2 = (PetscReal)(n2 - 2);
                                const PetscReal n1_norm = (cs1 - (1.0 + w1 / 2.0)) / (w1 / 2.0);
                                const PetscReal n2_norm = (cs2 - (1.0 + w2 / 2.0)) / (w2 / 2.0);
                                normal_velocity_mag = simCtx->icVelocityPhysical *
                                                      (1.0 - n1_norm * n1_norm) * (1.0 - n2_norm * n2_norm);
                                if (normal_velocity_mag < 0.0) normal_velocity_mag = 0.0;
                            }
                            break;
                        default:
                            LOG_ALLOW(LOCAL, LOG_WARNING, "Unrecognized initial-condition mode %d. Defaulting to zero.\n", simCtx->initialConditionMode);
                            break;
                    }

                    // Step B: apply flow direction and set the single contravariant flux component.
                    if (normal_velocity_mag != 0.0) {
                        const PetscReal signed_vel = normal_velocity_mag * flow_dir_sign * user->GridOrientation;
                        if (flow_axis == 0) {
                            const PetscReal area = sqrt(csi_arr[k][j][i].x * csi_arr[k][j][i].x +
                                                        csi_arr[k][j][i].y * csi_arr[k][j][i].y +
                                                        csi_arr[k][j][i].z * csi_arr[k][j][i].z);
                            ucont_val.x = signed_vel * area;
                        } else if (flow_axis == 1) {
                            const PetscReal area = sqrt(eta_arr[k][j][i].x * eta_arr[k][j][i].x +
                                                        eta_arr[k][j][i].y * eta_arr[k][j][i].y +
                                                        eta_arr[k][j][i].z * eta_arr[k][j][i].z);
                            ucont_val.y = signed_vel * area;
                        } else {
                            const PetscReal area = sqrt(zet_arr[k][j][i].x * zet_arr[k][j][i].x +
                                                        zet_arr[k][j][i].y * zet_arr[k][j][i].y +
                                                        zet_arr[k][j][i].z * zet_arr[k][j][i].z);
                            ucont_val.z = signed_vel * area;
                        }
                    }
                    ucont_arr[k][j][i] = ucont_val;
                } // end if(is_interior)
            }
        }
    }
    ierr = DMDAVecRestoreArray(user->fda, user->Ucont, &ucont_arr); CHKERRQ(ierr);

    // --- 5. Restore arrays ---
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi, &csi_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta, &eta_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet, &zet_arr); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LoadInitialUcont"
/**
 * @brief Load a staged file IC and return with Ucont populated.
 */
static PetscErrorCode LoadInitialUcont(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx *simCtx = user->simCtx;

    PetscFunctionBeginUser;
    ierr = PetscStrncpy(simCtx->_io_context_buffer, simCtx->initialConditionDirectory,
                        sizeof(simCtx->_io_context_buffer)); CHKERRQ(ierr);
    simCtx->current_io_directory = simCtx->_io_context_buffer;

    if (simCtx->initialConditionField == IC_FIELD_UCAT) {
        ierr = ReadFieldData(user, "ufield", user->Ucat, 0, "dat"); CHKERRQ(ierr);
        {
            const char *cell_fields[] = {"Ucat"};
            ierr = SynchronizePeriodicCellFields(user, 1, cell_fields); CHKERRQ(ierr);
        }
        ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);
        ierr = Cart2Contra(user); CHKERRQ(ierr);
    } else if (simCtx->initialConditionField == IC_FIELD_UCONT) {
        ierr = ReadFieldData(user, "vfield", user->Ucont, 0, "dat"); CHKERRQ(ierr);
    } else {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
                "Unsupported file initial-condition field selector %d.",
                simCtx->initialConditionField);
    }

    simCtx->current_io_directory = NULL;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PopulateInitialUcont"
/**
 * @brief Dispatch one fresh-start IC and return with Ucont populated.
 */
PetscErrorCode PopulateInitialUcont(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx *simCtx = user->simCtx;

    PetscFunctionBeginUser;
    if (simCtx->initialConditionMode == IC_MODE_FILE) {
        ierr = LoadInitialUcont(user); CHKERRQ(ierr);
    } else {
        ierr = SetInitialInteriorField(user, "Ucont"); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FinalizeBlockState"
/**
 * @brief Internal helper implementation: `FinalizeBlockState()`.
 * @details Local to this translation unit.
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
    const char *staggered_fields[] = {"Ucont"};
    ierr = SynchronizePeriodicStaggeredFields(user, 1, staggered_fields); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_TRACE," Ucont field ghosts updated.\n");
    
    // 3. Convert to Cartesian velocity.
    ierr = Contra2Cart(user); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_TRACE," Converted Ucont to Ucat.\n");

    // 4. Finalize periodic endpoint values, then refresh local Cartesian velocity.
    {
        const char *cell_fields[] = {"Ucat"};
        ierr = SynchronizePeriodicCellFields(user, 1, cell_fields); CHKERRQ(ierr);
    }
    ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_TRACE," Ucat field ghosts updated.\n"); 

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "SetInitialFluidState_FreshStart"
/**
 * @brief Internal helper implementation: `SetInitialFluidState_FreshStart()`.
 * @details Local to this translation unit.
 */
static PetscErrorCode SetInitialFluidState_FreshStart(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserCtx *user_finest = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;
    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d, Block %d: Setting t=0 state.\n", simCtx->rank, bi);

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
             const char *staggered_fields[] = {"Ucont"};
             ierr = SynchronizePeriodicStaggeredFields(&user_finest[bi], 1, staggered_fields); CHKERRQ(ierr);
             ierr = UpdateLocalGhosts(&user_finest[bi], "Ucat"); CHKERRQ(ierr);
        }
    }

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetInitialFluidState_Load"
/**
 * @brief Internal helper implementation: `SetInitialFluidState_Load()`.
 * @details Local to this translation unit.
 */
static PetscErrorCode SetInitialFluidState_Load(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserCtx *user_finest = simCtx->usermg.mgctx[simCtx->mglevels - 1].user;
    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d, Block %d: Reading restart files for step %d.\n",
                  simCtx->rank, bi, simCtx->StartStep);

        // ReadSimulationFields handles all file I/O for one block.
        ierr = ReadSimulationFields(&user_finest[bi], simCtx->StartStep); CHKERRQ(ierr);
        
        // Apply Boundary Conditions on Read fields
        ierr = ApplyBoundaryConditions(&user_finest[bi]); CHKERRQ(ierr);
        // After reading from a file, the local ghost regions MUST be updated
        // to ensure consistency across process boundaries for the first time step.
        //ierr = UpdateLocalGhosts(&user_finest[bi], "Ucat"); CHKERRQ(ierr);
        //ierr = UpdateLocalGhosts(&user_finest[bi], "P"); CHKERRQ(ierr);
        // ... add ghost updates for any other fields read from file ...
    }
    
    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "InitializeEulerianState"
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
            LOG_ALLOW(GLOBAL,LOG_INFO,"Initializing Analytical Solution type: %s (t=%.4f, step=%d).\n",simCtx->AnalyticalSolutionType,simCtx->StartTime,simCtx->StartStep);
            ierr = AnalyticalSolutionEngine(simCtx);
        }
        else{
            LOG_ALLOW(GLOBAL, LOG_INFO, "Starting from RESTART files (t=%.4f, step=%d).\n",
                    simCtx->StartTime, simCtx->StartStep);
            ierr = SetInitialFluidState_Load(simCtx); CHKERRQ(ierr);
        }
    } else { // StartStep = 0
        LOG_ALLOW(GLOBAL, LOG_INFO, "Performing a FRESH START (t=0, step=0).\n");
        if(strcmp(simCtx->eulerianSource,"solve")==0){
            ierr = SetInitialFluidState_FreshStart(simCtx); CHKERRQ(ierr);
        }else if(strcmp(simCtx->eulerianSource,"load")==0){
            LOG_ALLOW(GLOBAL,LOG_INFO,"FRESH START in LOAD mode. Reading files (t=%.4f,step=%d).\n",
                      simCtx->StartTime,simCtx->StartStep);
            ierr=SetInitialFluidState_Load(simCtx);CHKERRQ(ierr);
        }else if(strcmp(simCtx->eulerianSource,"analytical")==0){
            LOG_ALLOW(GLOBAL,LOG_INFO,"FRESH START in ANALYTICAL mode. Initializing Analytical Solution type: %s (t=%.4f,step=%d).\n",
                      simCtx->AnalyticalSolutionType,simCtx->StartTime,simCtx->StartStep);
            ierr=AnalyticalSolutionEngine(simCtx);CHKERRQ(ierr);
        }
    }

    // This crucial step, taken from the end of the legacy setup, ensures
    // that the history vectors (Ucont_o, Ucont_rm1, etc.) are correctly
    // populated before the first call to the time-stepping loop.
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        ierr = UpdateSolverHistoryVectors(&user_finest[bi]); CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "--- Eulerian State Initialized and History Vectors Populated ---\n");

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
