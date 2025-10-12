#include "postprocessing_kernels.h"

// =========== Dimensionalization Kernels ========================
#undef __FUNCT__
#define __FUNCT__ "DimensionalizeField"
/**
 * @brief Scales a specified field from non-dimensional to dimensional units in-place.
 *
 * This function acts as a dispatcher. It takes the string name of a field,
 * identifies the corresponding PETSc Vec object and the correct physical
 * scaling factor (e.g., U_ref for velocity, P_ref for pressure), and then
 * performs an in-place VecScale operation. It correctly handles the different
 * physical dimensions of Cartesian velocity vs. contravariant volume flux.
 *
 * @param[in,out] user        The UserCtx containing the PETSc Vecs to be modified.
 * @param[in]     field_name  The case-insensitive string name of the field to dimensionalize
 *                            (e.g., "Ucat", "P", "Ucont", "Coordinates", "ParticlePosition", "ParticleVelocity").
 * @return PetscErrorCode
 */
PetscErrorCode DimensionalizeField(UserCtx *user, const char *field_name)
{
    PetscErrorCode ierr;
    SimCtx         *simCtx = user->simCtx;
    Vec            target_vec = NULL;
    PetscReal      scale_factor = 1.0;
    char           field_type[64] = "Unknown";
    PetscBool      is_swarm_field = PETSC_FALSE; // Flag for special swarm handling
    const char     *swarm_field_name = NULL;    // Name of the field within the swarm

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    if (!user) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx is NULL.");
    if (!field_name) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "field_name is NULL.");

    // --- 1. Identify the target Vec and the correct scaling factor ---
    if (strcasecmp(field_name, "Ucat") == 0) {
        target_vec = user->Ucat;
        scale_factor = simCtx->scaling.U_ref;
        strcpy(field_type, "Cartesian Velocity (L/T)");
    } else if (strcasecmp(field_name, "Ucont") == 0) {
        target_vec = user->Ucont;
        scale_factor = simCtx->scaling.U_ref * simCtx->scaling.L_ref * simCtx->scaling.L_ref;
        strcpy(field_type, "Contravariant Volume Flux (L^3/T)");
    } else if (strcasecmp(field_name, "P") == 0) {
        target_vec = user->P;
        scale_factor = simCtx->scaling.P_ref;
        strcpy(field_type, "Pressure (M L^-1 T^-2)");
    } else if (strcasecmp(field_name, "Coordinates") == 0) {
        ierr = DMGetCoordinates(user->da, &target_vec); CHKERRQ(ierr);
        scale_factor = simCtx->scaling.L_ref;
        strcpy(field_type, "Grid Coordinates (L)");
    } else if (strcasecmp(field_name, "ParticlePosition") == 0) {
        is_swarm_field = PETSC_TRUE;
        swarm_field_name = "position";
        scale_factor = simCtx->scaling.L_ref;
        strcpy(field_type, "Particle Position (L)");
    } else if (strcasecmp(field_name, "ParticleVelocity") == 0) {
        is_swarm_field = PETSC_TRUE;
        swarm_field_name = "velocity";
        scale_factor = simCtx->scaling.U_ref;
        strcpy(field_type, "Particle Velocity (L/T)");
    } else {
        LOG(GLOBAL, LOG_WARNING, "DimensionalizeField: Unknown or unhandled field_name '%s'. Field will not be scaled.\n", field_name);
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }
    
    // --- 2. Check for trivial scaling ---
    if (PetscAbsReal(scale_factor - 1.0) < PETSC_MACHINE_EPSILON) {
        LOG(GLOBAL, LOG_DEBUG, "DimensionalizeField: Scaling factor for '%s' is 1.0. Skipping operation.\n", field_name);
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    // --- 3. Perform the in-place scaling operation ---
    LOG(GLOBAL, LOG_INFO, "Scaling '%s' field (%s) by factor %.4e.\n", field_name, field_type, scale_factor);

    if (is_swarm_field) {
        // Special handling for DMSwarm fields
        ierr = DMSwarmCreateGlobalVectorFromField(user->swarm, swarm_field_name, &target_vec); CHKERRQ(ierr);
        ierr = VecScale(target_vec, scale_factor); CHKERRQ(ierr);
        ierr = DMSwarmDestroyGlobalVectorFromField(user->swarm, swarm_field_name, &target_vec); CHKERRQ(ierr);
    } else {
        // Standard handling for PETSc Vecs
        if (target_vec) {
            ierr = VecScale(target_vec, scale_factor); CHKERRQ(ierr);
        } else {
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "Target vector for field '%s' was not found or is NULL.", field_name);
        }
    }
    
    // --- 4. Post-scaling updates for special cases ---
    if (strcasecmp(field_name, "Coordinates") == 0) {
        Vec l_coords;
        ierr = DMGetCoordinates(user->da, &target_vec); CHKERRQ(ierr); // Re-get the global vector
        ierr = DMGetCoordinatesLocal(user->da, &l_coords); CHKERRQ(ierr);
        ierr = DMGlobalToLocalBegin(user->fda, target_vec, INSERT_VALUES, l_coords); CHKERRQ(ierr);
        ierr = DMGlobalToLocalEnd(user->fda, target_vec, INSERT_VALUES, l_coords); CHKERRQ(ierr);
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DimensionalizeAllLoadedFields"
/**
 * @brief Orchestrates the dimensionalization of all relevant fields loaded from a file.
 *
 * This function is intended to be called in the post-processor immediately after
 * all solver output has been read into memory. It calls DimensionalizeField() for each of the core
 * physical quantities to convert the entire loaded state from non-dimensional to
 * dimensional units, preparing it for analysis and visualization.
 *
 * @param[in,out] user The UserCtx containing all the fields to be dimensionalized.
 * @return PetscErrorCode
 */
PetscErrorCode DimensionalizeAllLoadedFields(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx         *simCtx = user->simCtx;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    LOG(GLOBAL, LOG_INFO, "--- Converting all loaded fields to dimensional units ---\n");

    // Scale the grid itself first
    ierr = DimensionalizeField(user, "Coordinates"); CHKERRQ(ierr);

    // Scale primary fluid fields
    ierr = DimensionalizeField(user, "Ucat"); CHKERRQ(ierr);
    ierr = DimensionalizeField(user, "Ucont"); CHKERRQ(ierr);
    ierr = DimensionalizeField(user, "P"); CHKERRQ(ierr);

    // If particles are present, scale their fields
    if (simCtx->np > 0 && user->swarm) {
        ierr = DimensionalizeField(user, "ParticlePosition"); CHKERRQ(ierr);
        ierr = DimensionalizeField(user, "ParticleVelocity"); CHKERRQ(ierr);
    }

    LOG(GLOBAL, LOG_INFO, "--- Field dimensionalization complete ---\n");

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

//============ Post-Processing Kernels ===========================

#undef __FUNCT__
#define __FUNCT__ "ComputeNodalAverage"
/**
 * @brief Computes node-centered data by averaging 8 surrounding cell-centered values,
 *        exactly replicating the legacy code's indexing and boundary handling.
 *
 * This kernel uses a stencil that averages the 8 cells from the corner (i,j,k) to
 * (i+1, j+1, k+1) and stores the result at the output node (i,j,k). This matches
 * the legacy code's behavior. It operates on the full range of output nodes necessary
 * for the subsampled grid, preventing zero-padding at the boundaries.
 *
 * @param user The UserCtx, providing access to DMs and Vecs.
 * @param in_field_name The string name of the source Vec (e.g., "P", "Ucat").
 * @param out_field_name The string name of the destination Vec (e.g., "P_nodal").
 * @return PetscErrorCode
 */
PetscErrorCode ComputeNodalAverage(UserCtx* user, const char* in_field_name, const char* out_field_name)
{
    PetscErrorCode ierr;
    Vec            in_vec_local = NULL, out_vec_global = NULL;
    DM             dm_in = NULL, dm_out = NULL;
    PetscInt       dof = 0;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(GLOBAL, LOG_INFO, "-> KERNEL: Running ComputeNodalAverage on '%s' -> '%s'.\n", in_field_name, out_field_name);

    // --- 1. Map string names to PETSc objects ---
    if (strcasecmp(in_field_name, "P") == 0)             { in_vec_local = user->lP;         dm_in = user->da;   dof = 1; }
    else if (strcasecmp(in_field_name, "Ucat") == 0)    { in_vec_local = user->lUcat;      dm_in = user->fda;  dof = 3; }
    else if (strcasecmp(in_field_name, "Psi") == 0)     { in_vec_local = user->lPsi;       dm_in = user->da;   dof = 1; }
    // ... (add other fields as needed) ...
    else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unknown input field name for nodal averaging: %s", in_field_name);

    if (strcasecmp(out_field_name, "P_nodal") == 0)      { out_vec_global = user->P_nodal;    dm_out = user->da; }
    else if (strcasecmp(out_field_name, "Ucat_nodal") == 0) { out_vec_global = user->Ucat_nodal; dm_out = user->fda; }
    else if (strcasecmp(out_field_name, "Psi_nodal") == 0)   { out_vec_global = user->Psi_nodal;  dm_out = user->da; }
    // ... (add other fields as needed) ...
    else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unknown output field name for nodal averaging: %s", out_field_name);

    // --- 2. Ensure Input Data Ghosts are Up-to-Date ---
    ierr = UpdateLocalGhosts(user, in_field_name); CHKERRQ(ierr);

    // --- 3. Get DMDA info and array pointers ---
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(dm_out, &info); CHKERRQ(ierr);

    if (dof == 1) { // --- Scalar Field Averaging ---
        const PetscReal ***l_in_arr;
        PetscReal       ***g_out_arr;
        ierr = DMDAVecGetArrayRead(dm_in,in_vec_local, (void*)&l_in_arr); CHKERRQ(ierr);
        ierr = DMDAVecGetArray(dm_out,out_vec_global, (void*)&g_out_arr); CHKERRQ(ierr);

        // Loop over the output NODE locations. The loop bounds match the required
        // size of the final subsampled grid.
        for (PetscInt k = info.zs; k < info.zs + info.zm - 1; k++) {
            for (PetscInt j = info.ys; j < info.ys + info.ym - 1; j++) {
                for (PetscInt i = info.xs; i < info.xs + info.xm - 1; i++) {
                    g_out_arr[k][j][i] = 0.125 * (l_in_arr[k][j][i]     + l_in_arr[k][j][i+1] +
                                                  l_in_arr[k][j+1][i]   + l_in_arr[k][j+1][i+1] +
                                                  l_in_arr[k+1][j][i]   + l_in_arr[k+1][j][i+1] +
                                                  l_in_arr[k+1][j+1][i] + l_in_arr[k+1][j+1][i+1]);
                }
            }
        }
        ierr = DMDAVecRestoreArrayRead(dm_in,in_vec_local, (void*)&l_in_arr); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(dm_out,out_vec_global, (void*)&g_out_arr); CHKERRQ(ierr);

    } else if (dof == 3) { // --- Vector Field Averaging ---
        const Cmpnts ***l_in_arr;
        Cmpnts       ***g_out_arr;
        ierr = DMDAVecGetArrayRead(dm_in,in_vec_local, (void*)&l_in_arr); CHKERRQ(ierr);
        ierr = DMDAVecGetArray(dm_out,out_vec_global, (void*)&g_out_arr); CHKERRQ(ierr);

        for (PetscInt k = info.zs; k < info.zs + info.zm - 1; k++) {
            for (PetscInt j = info.ys; j < info.ys + info.ym - 1; j++) {
                for (PetscInt i = info.xs; i < info.xs + info.xm - 1; i++) {
                    g_out_arr[k][j][i].x = 0.125 * (l_in_arr[k][j][i].x + l_in_arr[k][j][i+1].x +
                                                    l_in_arr[k][j+1][i].x + l_in_arr[k][j+1][i+1].x +
                                                    l_in_arr[k+1][j][i].x + l_in_arr[k+1][j][i+1].x +
                                                    l_in_arr[k+1][j+1][i].x + l_in_arr[k+1][j+1][i+1].x);

                    g_out_arr[k][j][i].y = 0.125 * (l_in_arr[k][j][i].y + l_in_arr[k][j][i+1].y +
                                                    l_in_arr[k][j+1][i].y + l_in_arr[k][j+1][i+1].y +
                                                    l_in_arr[k+1][j][i].y + l_in_arr[k+1][j][i+1].y +
                                                    l_in_arr[k+1][j+1][i].y + l_in_arr[k+1][j+1][i+1].y);

                    g_out_arr[k][j][i].z = 0.125 * (l_in_arr[k][j][i].z + l_in_arr[k][j][i+1].z +
                                                    l_in_arr[k][j+1][i].z + l_in_arr[k][j+1][i+1].z +
                                                    l_in_arr[k+1][j][i].z + l_in_arr[k+1][j][i+1].z +
                                                    l_in_arr[k+1][j+1][i].z + l_in_arr[k+1][j+1][i+1].z);
                }
            }
        }
        ierr = DMDAVecRestoreArrayRead(dm_in,in_vec_local, (void*)&l_in_arr); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(dm_out,out_vec_global, (void*)&g_out_arr); CHKERRQ(ierr);
    }
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ComputeQCriterion"
/**
 * @brief Computes the Q-Criterion, a scalar value identifying vortex cores.
 *
 * This function is self-contained. It ensures all its required input data
 * (Ucat and grid metrics) have up-to-date ghost values before proceeding with
 * the calculation. The result is stored in the global vector user->Qcrit.
 *
 * @param user The UserCtx, providing access to all necessary data.
 * @return PetscErrorCode
 */
PetscErrorCode ComputeQCriterion(UserCtx* user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    const Cmpnts   ***lucat, ***lcsi, ***leta, ***lzet;
    const PetscReal***laj, ***lnvert;
    PetscReal      ***gq;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(GLOBAL, LOG_INFO, "-> KERNEL: Running ComputeQCriterion.\n");

    // --- 1. Ensure all required ghost values are up-to-date ---
    ierr = UpdateLocalGhosts(user, "Ucat");  CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "Csi");   CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "Eta");   CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "Zet");   CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "Aj");    CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "Nvert"); CHKERRQ(ierr);

    // --- 2. Get DMDA info and array pointers ---
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    
    ierr = DMDAVecGetArrayRead(user->fda, user->lUcat,  (void*)&lucat);   CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi,   (void*)&lcsi);    CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta,   (void*)&leta);    CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet,   (void*)&lzet);    CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da,  user->lAj,    (void*)&laj);     CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da,  user->lNvert, (void*)&lnvert);  CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da,  user->Qcrit, (void*)&gq);       CHKERRQ(ierr);

    // --- 3. Define Loop Bounds for INTERIOR Cells ---
    PetscInt i_start = (info.xs == 0) ? 1 : info.xs;
    PetscInt i_end   = (info.xs + info.xm == info.mx) ? info.mx - 1 : info.xs + info.xm;
    PetscInt j_start = (info.ys == 0) ? 1 : info.ys;
    PetscInt j_end   = (info.ys + info.ym == info.my) ? info.my - 1 : info.ys + info.ym;
    PetscInt k_start = (info.zs == 0) ? 1 : info.zs;
    PetscInt k_end   = (info.zs + info.zm == info.mz) ? info.mz - 1 : info.zs + info.zm;

    // --- 4. Main Computation Loop ---
    for (PetscInt k = k_start; k < k_end; k++) {
        for (PetscInt j = j_start; j < j_end; j++) {
            for (PetscInt i = i_start; i < i_end; i++) {
                
                // Calculate velocity derivatives in computational space (central differences)
                PetscReal uc = 0.5 * (lucat[k][j][i+1].x - lucat[k][j][i-1].x);
                PetscReal vc = 0.5 * (lucat[k][j][i+1].y - lucat[k][j][i-1].y);
                PetscReal wc = 0.5 * (lucat[k][j][i+1].z - lucat[k][j][i-1].z);

                PetscReal ue = 0.5 * (lucat[k][j+1][i].x - lucat[k][j-1][i].x);
                PetscReal ve = 0.5 * (lucat[k][j+1][i].y - lucat[k][j-1][i].y);
                PetscReal we = 0.5 * (lucat[k][j+1][i].z - lucat[k][j-1][i].z);

                PetscReal uz = 0.5 * (lucat[k+1][j][i].x - lucat[k-1][j][i].x);
                PetscReal vz = 0.5 * (lucat[k+1][j][i].y - lucat[k-1][j][i].y);
                PetscReal wz = 0.5 * (lucat[k+1][j][i].z - lucat[k-1][j][i].z);

                // Average metrics to the cell center
                PetscReal csi1 = 0.5 * (lcsi[k][j][i].x + lcsi[k][j][i-1].x) * laj[k][j][i];
                PetscReal csi2 = 0.5 * (lcsi[k][j][i].y + lcsi[k][j][i-1].y) * laj[k][j][i];
                PetscReal csi3 = 0.5 * (lcsi[k][j][i].z + lcsi[k][j][i-1].z) * laj[k][j][i];
                
                PetscReal eta1 = 0.5 * (leta[k][j][i].x + leta[k][j-1][i].x) * laj[k][j][i];
                PetscReal eta2 = 0.5 * (leta[k][j][i].y + leta[k][j-1][i].y) * laj[k][j][i];
                PetscReal eta3 = 0.5 * (leta[k][j][i].z + leta[k][j-1][i].z) * laj[k][j][i];

                PetscReal zet1 = 0.5 * (lzet[k][j][i].x + lzet[k-1][j][i].x) * laj[k][j][i];
                PetscReal zet2 = 0.5 * (lzet[k][j][i].y + lzet[k-1][j][i].y) * laj[k][j][i];
                PetscReal zet3 = 0.5 * (lzet[k][j][i].z + lzet[k-1][j][i].z) * laj[k][j][i];

                // Calculate velocity gradient tensor components d_ij = du_i/dx_j
                PetscReal d11 = uc * csi1 + ue * eta1 + uz * zet1;
                PetscReal d12 = uc * csi2 + ue * eta2 + uz * zet2;
                PetscReal d13 = uc * csi3 + ue * eta3 + uz * zet3;

                PetscReal d21 = vc * csi1 + ve * eta1 + vz * zet1;
                PetscReal d22 = vc * csi2 + ve * eta2 + vz * zet2;
                PetscReal d23 = vc * csi3 + ve * eta3 + vz * zet3;
                
                PetscReal d31 = wc * csi1 + we * eta1 + wz * zet1;
                PetscReal d32 = wc * csi2 + we * eta2 + wz * zet2;
                PetscReal d33 = wc * csi3 + we * eta3 + wz * zet3;

                // Strain-Rate Tensor S_ij = 0.5 * (d_ij + d_ji)
                PetscReal s11 = d11;
                PetscReal s12 = 0.5 * (d12 + d21);
                PetscReal s13 = 0.5 * (d13 + d31);
                PetscReal s22 = d22;
                PetscReal s23 = 0.5 * (d23 + d32);
                PetscReal s33 = d33;

                // Vorticity Tensor Omega_ij = 0.5 * (d_ij - d_ji)
                PetscReal w12 = 0.5 * (d12 - d21);
                PetscReal w13 = 0.5 * (d13 - d31);
                PetscReal w23 = 0.5 * (d23 - d32);

                // Squared norms of the tensors
                PetscReal s_norm_sq = s11*s11 + s22*s22 + s33*s33 + 2.0*(s12*s12 + s13*s13 + s23*s23);
                PetscReal w_norm_sq = 2.0 * (w12*w12 + w13*w13 + w23*w23);
                
                gq[k][j][i] = 0.5 * (w_norm_sq - s_norm_sq);

                if (lnvert[k][j][i] > 0.1) {
                    gq[k][j][i] = 0.0;
                }
            }
        }
    }

    // --- 5. Restore arrays ---
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lUcat,  (void*)&lucat);   CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi,   (void*)&lcsi);    CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta,   (void*)&leta);    CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet,   (void*)&lzet);    CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da,  user->lAj,    (void*)&laj);     CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da,  user->lNvert, (void*)&lnvert);  CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da,  user->Qcrit, (void*)&gq);       CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "NormalizeRelativeField"
/**
 * @brief Normalizes a relative field by subtracting a reference value.
 *
 * This kernel finds the relative field at a specific grid point (i,j,k) and subtracts
 * this value from the entire field. The reference point is configurable via
 * command-line options (-ip, -jp, -kp). The operation is performed in-place
 * on the provided relative field vector.
 *
 * @param user The UserCtx, providing access to DMs and Vecs.
 * @param relative_field_name The string name of the relative field Vec to normalize (e.g., "P").
 * @return PetscErrorCode
 */
PetscErrorCode NormalizeRelativeField(UserCtx* user, const char* relative_field_name)
{
    PetscErrorCode ierr;
    Vec            P_vec = NULL;
    PetscMPIInt    rank;
    PetscInt       ip=1, jp=1, kp=1; // Default reference point
    PetscReal      p_ref = 0.0;
    PetscInt       ref_point_global_idx[1];
    PetscScalar    ref_value_local[1];
    IS             is_from, is_to;
    VecScatter     scatter_ctx;
    Vec            ref_vec_seq;
    PostProcessParams *pps = user->simCtx->pps;

    // Fetch referenc point from pps.
    ip = pps->reference[0];
    jp = pps->reference[1];
    kp = pps->reference[2];

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(GLOBAL, LOG_INFO, "-> KERNEL: Running NormalizeRelativeField on '%s'.\n", relative_field_name);

    // --- 1. Map string argument to the PETSc Vec ---
    if (strcasecmp(relative_field_name, "P") == 0) {
        P_vec = user->P;
    } else {
        SETERRQ(PETSC_COMM_SELF, 1, "NormalizeRelativeField only supports the primary 'P' field , not '%s' currently.", relative_field_name);
    }

    // --- 2. Get reference point from options and calculate its global DA index ---
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Convert the (i,j,k) logical grid coordinates to the global 1D index used by the DMDA vector
    ref_point_global_idx[0] = kp * (user->IM * user->JM) + jp * user->IM + ip;

    // --- 3. Robustly Scatter the single reference value to rank 0 ---
    // Create an Index Set (IS) for the source point (from the global vector)
    ierr = ISCreateGeneral(PETSC_COMM_WORLD, 1, ref_point_global_idx, PETSC_COPY_VALUES, &is_from); CHKERRQ(ierr);

    // Create a sequential vector on rank 0 to hold the result
    ierr = VecCreateSeq(PETSC_COMM_SELF, 1, &ref_vec_seq); CHKERRQ(ierr);
    
    // Create an Index Set for the destination point (index 0 of the new sequential vector)
    PetscInt dest_idx[1] = {0};
    ierr = ISCreateGeneral(PETSC_COMM_SELF, 1, dest_idx, PETSC_COPY_VALUES, &is_to); CHKERRQ(ierr);

    // Create the scatter context and perform the scatter
    ierr = VecScatterCreate(P_vec, is_from, ref_vec_seq, is_to, &scatter_ctx); CHKERRQ(ierr);
    ierr = VecScatterBegin(scatter_ctx, P_vec, ref_vec_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
    ierr = VecScatterEnd(scatter_ctx, P_vec, ref_vec_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);

    // On rank 0, get the value. On other ranks, this will do nothing.
    if (rank == 0) {
        ierr = VecGetValues(ref_vec_seq, 1, dest_idx, ref_value_local); CHKERRQ(ierr);
        p_ref = ref_value_local[0];
        LOG_ALLOW(LOCAL, LOG_DEBUG, "%s reference point (%" PetscInt_FMT ", %" PetscInt_FMT ", %" PetscInt_FMT ") has value %g.\n", relative_field_name, jp, kp, ip, p_ref);
    }
    
    // --- 4. Broadcast the reference value from rank 0 to all other processes ---
    ierr = MPI_Bcast(&p_ref, 1, MPIU_REAL, 0, PETSC_COMM_WORLD); CHKERRQ(ierr);

    // --- 5. Perform the normalization (in-place shift) on the full distributed vector ---
    ierr = VecShift(P_vec, -p_ref); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "%s field normalized by subtracting %g.\n", relative_field_name, p_ref);
    
    // --- 6. Cleanup ---
    ierr = ISDestroy(&is_from); CHKERRQ(ierr);
    ierr = ISDestroy(&is_to); CHKERRQ(ierr);
    ierr = VecScatterDestroy(&scatter_ctx); CHKERRQ(ierr);
    ierr = VecDestroy(&ref_vec_seq); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

// ===========================================================================
// Particle Post-Processing Kernels
// ===========================================================================
#undef __FUNCT__
#define __FUNCT__ "ComputeSpecificKE"
/**
 * @brief Computes the specific kinetic energy (KE per unit mass) for each particle.
 *
 * This kernel calculates SKE = 0.5 * |velocity|^2. It requires that the
 * velocity field exists and will populate the specific kinetic energy field.
 * The output field must be registered before this kernel is called.
 *
 * @param user           The UserCtx containing the DMSwarm.
 * @param velocity_field The name of the input vector field for particle velocity.
 * @param ske_field      The name of the output scalar field to store specific KE.
 * @return PetscErrorCode
 */
PetscErrorCode ComputeSpecificKE(UserCtx* user, const char* velocity_field, const char* ske_field)
{
    PetscErrorCode ierr;
    PetscInt n_local;
    const PetscScalar (*vel_arr)[3]; // Access velocity as array of 3-component vectors
    PetscScalar *ske_arr;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(GLOBAL, LOG_INFO, "-> KERNEL: Running ComputeSpecificKE ('%s' -> '%s').\n", velocity_field, ske_field);

    // Get local data arrays from the DMSwarm
    ierr = DMSwarmGetLocalSize(user->swarm, &n_local); CHKERRQ(ierr);
    if (n_local == 0) PetscFunctionReturn(0);

    // Get read-only access to velocity and write access to the output field
    ierr = DMSwarmGetField(user->swarm, velocity_field, NULL, NULL, (const void**)&vel_arr); CHKERRQ(ierr);
    ierr = DMSwarmGetField(user->post_swarm, ske_field, NULL, NULL, (void**)&ske_arr); CHKERRQ(ierr);

    // Main computation loop
    for (PetscInt p = 0; p < n_local; p++) {
        const PetscScalar u = vel_arr[p][0];
        const PetscScalar v = vel_arr[p][1];
        const PetscScalar w = vel_arr[p][2];
        const PetscScalar vel_sq = u*u + v*v + w*w;
        ske_arr[p] = 0.5 * vel_sq;
    }

    // Restore arrays
    ierr = DMSwarmRestoreField(user->swarm, velocity_field, NULL, NULL, (const void**)&vel_arr); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(user->post_swarm, ske_field, NULL, NULL, (void**)&ske_arr); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

