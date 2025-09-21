/**
 * @file initialcondition.c  //  Setup the Initial conditions for different cases. 
 * @brief Test program for DMSwarm interpolation using the fdf-curvIB method.
 **/

 #include "initialcondition.h"

/**
 * @brief Sets the initial values for the INTERIOR of a specified Eulerian field.
 *
 * This function initializes the interior nodes of `Ucont` based on a profile selected
 * by `user->FieldInitialization`. It explicitly skips any node that lies on a global
 * boundary, as those values are set by the Boundary System's `Initialize` methods.
 *
 * The initialization is directional, aligned with the primary INLET face that was
 * identified by the parser. This ensures the initial flow is physically meaningful.
 *
 * Supported `user->FieldInitialization` profiles for "Ucont":
 *  - 0: Zero Velocity. All interior components of Ucont are set to 0.
 *  - 1: Constant Normal Velocity. The contravariant velocity component normal to the
 *       inlet direction is set such that the physical velocity normal to those grid
 *       planes is a constant `uin`. Other contravariant components are zero.
 *  - 2: Poiseuille Normal Velocity. The contravariant component normal to the
 *       inlet direction is set with a parabolic profile.
 *
 * @param user      The main UserCtx struct, containing all simulation data and configuration.
 * @param fieldName A string ("Ucont" or "P") identifying which field to initialize.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode SetInitialInteriorField(UserCtx *user, const char *fieldName)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    SimCtx *simCtx = user->simCtx;
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Setting initial INTERIOR field for '%s' with profile %d.\n", fieldName, simCtx->FieldInitialization);

    // This function currently only implements logic for Ucont.
    if (strcmp(fieldName, "Ucont") != 0) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Skipping SetInitialInteriorField for non-Ucont field '%s'.\n", fieldName);
        PetscFunctionReturn(0);
    }

    // --- 1. Get references to required data and PETSc arrays ---
    DM            fieldDM = user->fda;
    Vec           fieldVec = user->Ucont;
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(fieldDM, &info); CHKERRQ(ierr);

    Vec      localCoor;
    Cmpnts ***coor_arr;
    ierr = DMGetCoordinatesLocal(user->da, &localCoor); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, localCoor, &coor_arr); CHKERRQ(ierr);
    
    Cmpnts ***csi_arr, ***eta_arr, ***zet_arr;
    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi, &csi_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta, &eta_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet, &zet_arr); CHKERRQ(ierr);
   
    // --- 2. Compute Cell-Center Coordinates (only if needed by the selected profile) ---
    Cmpnts ***cent_coor = NULL;
    PetscInt xs_cell=0, xm_cell=0, ys_cell=0, ym_cell=0, zs_cell=0, zm_cell=0;
    
    if (simCtx->FieldInitialization == 2) { // Profile 2 (Poiseuille) requires cell centers.
      /*
        ierr = GetOwnedCellRange(&info, 0, &xs_cell, &xm_cell); CHKERRQ(ierr);
        ierr = GetOwnedCellRange(&info, 1, &ys_cell, &ym_cell); CHKERRQ(ierr);
        ierr = GetOwnedCellRange(&info, 2, &zs_cell, &zm_cell); CHKERRQ(ierr);
      */
      ierr = GetGhostedCellRange(&info,&user->neighbors,0, &xs_cell, &xm_cell); CHKERRQ(ierr);
      ierr = GetGhostedCellRange(&info,&user->neighbors,1, &ys_cell, &ym_cell); CHKERRQ(ierr);
      ierr = GetGhostedCellRange(&info,&user->neighbors,2, &zs_cell, &zm_cell); CHKERRQ(ierr);
      
        if (xm_cell > 0 && ym_cell > 0 && zm_cell > 0) {
            ierr = Allocate3DArray(&cent_coor, zm_cell, ym_cell, xm_cell); CHKERRQ(ierr);
            ierr = InterpolateFieldFromCornerToCenter(coor_arr, cent_coor, user); CHKERRQ(ierr);
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Computed temporary cell-center coordinates for Poiseuille profile.\n");
        }
    }
    
    // --- 3. Loop Over Owned Nodes and Apply Initial Condition to Interior ---
    Cmpnts ***ucont_arr;
    ierr = DMDAVecGetArray(fieldDM, fieldVec, &ucont_arr); CHKERRQ(ierr);
    
    PetscInt i, j, k;
    const PetscInt mx = info.mx, my = info.my, mz = info.mz; // Global node dimensions
    const PetscInt xs = info.xs, xe = info.xs + info.xm;
    const PetscInt ys = info.ys, ye = info.ys + info.ym;
    const PetscInt zs = info.zs, ze = info.zs + info.zm;
    
    const Cmpnts uin = simCtx->InitialConstantContra; // Max/average velocity from user options
    // Flow into a negative face (e.g., -Zeta at k=0) is in the positive physical direction (+z).
    // Flow into a positive face (e.g., +Zeta at k=mz-1) is in the negative physical direction (-z).

    LOG_ALLOW(GLOBAL,LOG_DEBUG,"InitialConstantContra = (%.3f, %.3f, %.3f)\n",(double)uin.x, (double)uin.y, (double)uin.z);
    
    const PetscReal flow_direction_sign = (user->identifiedInletBCFace % 2 == 0) ? 1.0 : -1.0;
        
    for (k = zs; k < ze; k++) {
        for (j = ys; j < ye; j++) {
            for (i = xs; i < xe; i++) {
                
                // The crucial check to ensure we only modify interior nodes.
                const PetscBool is_interior = (i > 0 && i < mx - 1 &&
                                               j > 0 && j < my - 1 &&
                                               k > 0 && k < mz - 1);

                if (is_interior) {
                    Cmpnts ucont_val = {0.0, 0.0, 0.0}; // Default to zero velocity
                    PetscReal normal_velocity_mag = 0.0;

                    // Step A: Determine the magnitude of the desired physical normal velocity.
                    switch (simCtx->FieldInitialization) {
                        case 0: // Zero initial velocity
                            normal_velocity_mag = 0.0;
                            break;
 		        case 1: /* Constant Normal Velocity */
		            if (user->identifiedInletBCFace == BC_FACE_NEG_X ||
			    user->identifiedInletBCFace == BC_FACE_POS_X) {
			    normal_velocity_mag = uin.x;
		            } else if (user->identifiedInletBCFace == BC_FACE_NEG_Y || user->identifiedInletBCFace == BC_FACE_POS_Y) {
			   normal_velocity_mag = uin.y;
		            } else {
			   normal_velocity_mag = uin.z;
		            }
		      break; 
                        case 2: // Poiseuille Normal Velocity
                            {
			                                      // This profile assumes flow is aligned with the k-index direction.
                                // It uses grid indices (i,j) to define the cross-section, which works for bent geometries.
                                PetscReal u0 = 0.0;
                                if (user->identifiedInletBCFace <= BC_FACE_POS_X) u0 = uin.x;
                                else if (user->identifiedInletBCFace <= BC_FACE_POS_Y) u0 = uin.y;
                                else u0 = uin.z; // Assumes Z-like inlet direction

                                // Define channel geometry in "index space" based on global grid dimensions
                                // We subtract 2.0 because the interior runs from index 1 to mx-2 (or my-2).
                                const PetscReal i_width  = (PetscReal)(mx - 2);
                                const PetscReal j_width  = (PetscReal)(my - 2);
                                const PetscReal i_center = 1.0 + i_width / 2.0;
                                const PetscReal j_center = 1.0 + j_width / 2.0;

                                // Create normalized coordinates for the current point (i,j), ranging from -1 to 1
                                const PetscReal i_norm = (i - i_center) / (i_width / 2.0);
                                const PetscReal j_norm = (j - j_center) / (j_width / 2.0);

                                // Apply the parabolic profile for a rectangular/square channel
                                // V(i,j) = V_max * (1 - i_norm^2) * (1 - j_norm^2)
                                const PetscReal profile_i = 1.0 - i_norm * i_norm;
                                const PetscReal profile_j = 1.0 - j_norm * j_norm;
                                normal_velocity_mag = u0 * profile_i * profile_j;

                                // Clamp to zero for any points outside the channel (or due to minor float errors)
                                if (normal_velocity_mag < 0.0) {
                                    normal_velocity_mag = 0.0;
                                }
                            }
                            break;
			    /*
                                PetscReal r_sq = 0.0;
                                const PetscInt i_local = i - xs_cell, j_local = j - ys_cell, k_local = k - zs_cell;
                                if (cent_coor && i_local >= 0 && i_local < xm_cell && j_local >= 0 && j_local < ym_cell && k_local >= 0 && k_local < zm_cell) {
                                    const Cmpnts* center = &cent_coor[k_local][j_local][i_local];
                                    if (user->identifiedInletBCFace <= BC_FACE_POS_X) r_sq = center->y * center->y + center->z * center->z;
                                    else if (user->identifiedInletBCFace <= BC_FACE_POS_Y) r_sq = center->x * center->x + center->z * center->z;
                                    else r_sq = center->x * center->x + center->y * center->y;
				    // pick the correct contravariant component for center‐line speed 
				    PetscReal u0;
				    if (user->identifiedInletBCFace == BC_FACE_NEG_X ||
					user->identifiedInletBCFace == BC_FACE_POS_X) {
				      u0 = uin.x;
				    } else if (user->identifiedInletBCFace == BC_FACE_NEG_Y ||
					       user->identifiedInletBCFace == BC_FACE_POS_Y) {
				      u0 = uin.y;
				    } else {
				      u0 = uin.z;
				    }
				    // now form the parabolic profile as before 
                                    normal_velocity_mag = 2.0 * u0 * (1.0 - 4.0 * r_sq);
                                }
                            }
                            break;
			    */
                        default:
                            LOG_ALLOW(LOCAL, LOG_WARNING, "Unrecognized FieldInitialization profile %d. Defaulting to zero.\n", simCtx->FieldInitialization);
                            normal_velocity_mag = 0.0;
                            break;
                    }

                    // Step B: Apply direction sign and set the correct contravariant component.
                    // The contravariant component U^n = v_n * Area_n, where v_n is the physical normal velocity.
                    if (normal_velocity_mag != 0.0) {
		       const PetscReal signed_normal_vel = normal_velocity_mag * flow_direction_sign*user->GridOrientation;
                        
                        if (user->identifiedInletBCFace == BC_FACE_NEG_X || user->identifiedInletBCFace == BC_FACE_POS_X) {
                            const PetscReal area_i = sqrt(csi_arr[k][j][i].x * csi_arr[k][j][i].x + csi_arr[k][j][i].y * csi_arr[k][j][i].y + csi_arr[k][j][i].z * csi_arr[k][j][i].z);
			    
                            ucont_val.x = signed_normal_vel * area_i;

			    LOG_LOOP_ALLOW(GLOBAL,LOG_DEBUG,k,50," ucont_val.x = %.6f (signed_normal_vel=%.3f × area=%.4f)\n",ucont_val.x, signed_normal_vel, area_i);
                        } 
                        else if (user->identifiedInletBCFace == BC_FACE_NEG_Y || user->identifiedInletBCFace == BC_FACE_POS_Y) {
                            const PetscReal area_j = sqrt(eta_arr[k][j][i].x * eta_arr[k][j][i].x + eta_arr[k][j][i].y * eta_arr[k][j][i].y + eta_arr[k][j][i].z * eta_arr[k][j][i].z);
			     
                            ucont_val.y = signed_normal_vel * area_j;

			    LOG_LOOP_ALLOW(GLOBAL,LOG_DEBUG,k,50," ucont_val.y = %.6f (signed_normal_vel=%.3f × area=%.4f)\n",ucont_val.y, signed_normal_vel, area_j);
                        } 
                        else { // Z-inlet
                            const PetscReal area_k = sqrt(zet_arr[k][j][i].x * zet_arr[k][j][i].x + zet_arr[k][j][i].y * zet_arr[k][j][i].y + zet_arr[k][j][i].z * zet_arr[k][j][i].z);

                            ucont_val.z = signed_normal_vel * area_k;

			    LOG_LOOP_ALLOW(GLOBAL,LOG_DEBUG,k,50," i,j,k,ucont_val.z = %d, %d, %d, %.6f (signed_normal_vel=%.3f × area=%.4f)\n",i,j,k,ucont_val.z, signed_normal_vel, area_k);
                        }
                    }
                    ucont_arr[k][j][i] = ucont_val;
                } // end if(is_interior)
            }
        }
    }
    ierr = DMDAVecRestoreArray(fieldDM, fieldVec, &ucont_arr); CHKERRQ(ierr);

    // --- 5. Cleanup: Restore arrays and free temporary memory ---
    ierr = DMDAVecRestoreArrayRead(user->fda, localCoor, &coor_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi, &csi_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta, &eta_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet, &zet_arr); CHKERRQ(ierr);
    
    if (cent_coor) {
        ierr = Deallocate3DArray(cent_coor, zm_cell, ym_cell); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

/**
 * @brief (HELPER) Finalizes a block's state by applying BCs, converting to Cartesian,
 * and synchronizing ghost regions.
 */
static PetscErrorCode FinalizeBlockState(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    // This sequence ensures a fully consistent state for a single block.
    // 1. Apply BCs using the information from InflowFlux/OutflowFlux.
    ierr = FormBCS(user); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_DEBUG," Boundary condition applied.\n");
    // 2. Sync contravariant velocity field.
    ierr = UpdateLocalGhosts(user, "Ucont"); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_DEBUG," Ucont field ghosts updated.\n");
    
    // 3. Convert to Cartesian velocity.
    ierr = Contra2Cart(user); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_DEBUG," Converted Ucont to Ucat.\n");

    // 4. Sync the new Cartesian velocity field.
    ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_DEBUG," Ucat field ghosts updated.\n"); 

    PetscFunctionReturn(0);
}


/**
 * @brief (HELPER) Sets the t=0 fluid state for all blocks.
 *
 * This function replicates the initialization sequence for a fresh start from
 * the legacy code. It sets an initial guess for the interior, establishes the
 * boundary conditions, and ensures the final state is consistent across MPI ranks.
 */
static PetscErrorCode SetInitialFluidState_FreshStart(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserCtx *user_finest = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;
    PetscFunctionBeginUser;

    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d, Block %d: Setting t=0 state.\n", simCtx->rank, bi);

        // 1. Set an initial guess for the INTERIOR of the domain.
        //    Replaces the legacy `if(InitialGuessOne)` block.

	LOG_ALLOW(GLOBAL,LOG_DEBUG," Initializing Interior Ucont field.\n");
        ierr = SetInitialInteriorField(&user_finest[bi], "Ucont"); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_DEBUG," Interior Ucont field initialized.\n");
	
        // 2. Establish the initial boundary flux values and set inlet profiles.
        //    This sequence is critical and comes directly from the legacy setup.
	LOG_ALLOW(GLOBAL,LOG_DEBUG," Boundary condition Preparation steps initiated.\n");
        ierr = InflowFlux(&user_finest[bi]); CHKERRQ(ierr);
        ierr = OutflowFlux(&user_finest[bi]); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_DEBUG," Boundary condition Preparation steps completed.\n");
        // 3. Apply all boundary conditions, convert to Cartesian, and sync ghosts.
        LOG_ALLOW(GLOBAL,LOG_DEBUG," Boundary condition application and state finalization initiated.\n");
        ierr = FinalizeBlockState(&user_finest[bi]); CHKERRQ(ierr);
	LOG_ALLOW(GLOBAL,LOG_DEBUG," Boundary condition application and state finalization complete.\n");
    }

    // If using multiple grid blocks, handle the interface conditions between them.
    if (simCtx->block_number > 1) {
      //   LOG_ALLOW(GLOBAL, LOG_INFO, "Updating multi-block interfaces for t=0.\n");
	//	ierr = Block_Interface_U(user_finest); CHKERRQ(ierr);
        // After interface update, ghost regions might be stale. Refresh them.
        for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
             ierr = UpdateLocalGhosts(&user_finest[bi], "Ucont"); CHKERRQ(ierr);
             ierr = UpdateLocalGhosts(&user_finest[bi], "Ucat"); CHKERRQ(ierr);
        }
    }

    PetscFunctionReturn(0);
}

/**
 * @brief (HELPER) Reads fluid state for all blocks from restart files.
 */
static PetscErrorCode SetInitialFluidState_Load(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserCtx *user_finest = simCtx->usermg.mgctx[simCtx->mglevels - 1].user;
    PetscFunctionBeginUser;

    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d, Block %d: Reading restart files for step %d.\n",
                  simCtx->rank, bi, simCtx->StartStep);

        // ReadSimulationFields handles all file I/O for one block.
        ierr = ReadSimulationFields(&user_finest[bi], simCtx->StartStep); CHKERRQ(ierr);
        
        // After reading from a file, the local ghost regions MUST be updated
        // to ensure consistency across process boundaries for the first time step.
        ierr = UpdateLocalGhosts(&user_finest[bi], "Ucont"); CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(&user_finest[bi], "Ucat"); CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(&user_finest[bi], "P"); CHKERRQ(ierr);
        // ... add ghost updates for any other fields read from file ...
    }
    
    PetscFunctionReturn(0);
}

/**
 * @brief High-level orchestrator to set the complete initial state of the Eulerian solver.
 *
 * This function is called once from main() before the time loop begins. It inspects
 * the simulation context to determine whether to perform a fresh start (t=0) or
 * restart from saved files. It then delegates to the appropriate helper function.
 * Finally, it initializes the solver's history vectors (Ucont_o, P_o, etc.)
 * to ensure the first time step has the necessary data.
 */
PetscErrorCode InitializeEulerianState(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserCtx *user_finest = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;

    PetscFunctionBeginUser;

    LOG_ALLOW(GLOBAL, LOG_INFO, "--- Initializing Eulerian State ---\n");

    if (simCtx->StartStep > 0) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Starting from RESTART files (t=%.4f, step=%d).\n",
                  simCtx->StartTime, simCtx->StartStep);
        ierr = SetInitialFluidState_Load(simCtx); CHKERRQ(ierr);
    } else {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Performing a FRESH START (t=0, step=0).\n");
        if(strcmp(simCtx->eulerianSource,"solve")==0){
            ierr = SetInitialFluidState_FreshStart(simCtx); CHKERRQ(ierr);
        }else if(strcmp(simCtx->eulerianSource,"load")==0){
            
            LOG_ALLOW(GLOBAL,LOG_INFO,"FRESH START in LOAD mode. Reading files (t=%.4f,step=%d).\n",
                      simCtx->StartTime,simCtx->StartStep);

            ierr=SetInitialFluidState_Load(simCtx);CHKERRQ(ierr);
        }
    }

    // This crucial step, taken from the end of the legacy setup, ensures
    // that the history vectors (Ucont_o, Ucont_rm1, etc.) are correctly
    // populated before the first call to the time-stepping loop.
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        ierr = UpdateSolverHistoryVectors(&user_finest[bi]); CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "--- Eulerian State Initialized and History Vectors Populated ---\n");
    PetscFunctionReturn(0);
}
