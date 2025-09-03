#include "postprocessing_kernels.h"

/**
 * @brief Computes node-centered data by averaging 8 surrounding cell-centered values.
 *
 * This kernel is a stateless transformation. It first ensures its input data has
 * up-to-date ghost values by calling UpdateLocalGhosts, then performs the averaging
 * stencil, writing the result to a global output vector.
 *
 * @param user The UserCtx, providing access to DMs and Vecs.
 * @param in_field_name The string name of the source Vec (e.g., "P", "Ucat").
 * @param out_field_name The string name of the destination Vec (e.g., "P_nodal").
 * @return PetscErrorCode
 */
PetscErrorCode ComputeNodalAverage(UserCtx* user, const char* in_field_name, const char* out_field_name)
{
    PetscErrorCode ierr;
    Vec            in_vec_global = NULL, in_vec_local = NULL, out_vec_global = NULL;
    DM             dm = NULL; // The DM associated with the data layout
    PetscInt       dof = 0;

    PetscFunctionBeginUser;
    LOG_ALLOW(GLOBAL, LOG_INFO, "-> KERNEL: Running ComputeNodalAverage on '%s' -> '%s'.\n", in_field_name, out_field_name);

    // --- 1. Map string names to PETSc objects ---
    if (strcmp(in_field_name, "P") == 0)             { in_vec_global = user->P;         in_vec_local = user->lP;         dm = user->da;   dof = 1; }
    else if (strcmp(in_field_name, "Ucat") == 0)    { in_vec_global = user->Ucat;      in_vec_local = user->lUcat;      dm = user->fda;  dof = 3; }
    else if (strcmp(in_field_name, "Qcrit") == 0) { in_vec_global = user->Qcrit;     in_vec_local = user->lP; /* Re-use lP for scalar data */ dm = user->da; dof = 1; }
    else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unknown input field name for nodal averaging: %s", in_field_name);

    if (strcmp(out_field_name, "P_nodal") == 0)      { out_vec_global = user->P_nodal; }
    else if (strcmp(out_field_name, "Ucat_nodal") == 0) { out_vec_global = user->Ucat_nodal; }
    else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unknown output field name for nodal averaging: %s", out_field_name);

    // --- 2. Ensure Input Data Ghosts are Up-to-Date ---
    ierr = UpdateLocalGhosts(user, in_field_name); CHKERRQ(ierr);

    // --- 3. Get DMDA info and array pointers ---
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(dm, &info); CHKERRQ(ierr);

    if (dof == 1) { // --- Scalar Field Averaging ---
        const PetscReal ***l_in_arr;
        PetscReal       ***g_out_arr;
        ierr = DMDAVecGetArrayRead(dm,in_vec_local, (void*)&l_in_arr); CHKERRQ(ierr);
        ierr = DMDAVecGetArray(dm,out_vec_global, (void*)&g_out_arr); CHKERRQ(ierr);

        for (PetscInt k = info.zs; k < info.zs + info.zm; k++) {
            if (k == 0 || k == info.mz - 1) continue;
            for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
                if (j == 0 || j == info.my - 1) continue;
                for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
                    if (i == 0 || i == info.mx - 1) continue;
                    
                    g_out_arr[k][j][i] = 0.125 * (l_in_arr[k-1][j-1][i-1] + l_in_arr[k-1][j-1][i] +
                                                  l_in_arr[k-1][j][i-1]   + l_in_arr[k-1][j][i]   +
                                                  l_in_arr[k][j-1][i-1]   + l_in_arr[k][j-1][i]   +
                                                  l_in_arr[k][j][i-1]     + l_in_arr[k][j][i]);
                }
            }
        }
        ierr = DMDAVecRestoreArrayRead(dm,in_vec_local, (void*)&l_in_arr); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(dm,out_vec_global, (void*)&g_out_arr); CHKERRQ(ierr);

    } else if (dof == 3) { // --- Vector Field Averaging ---
        const Cmpnts ***l_in_arr;
        Cmpnts       ***g_out_arr;
        ierr = DMDAVecGetArrayRead(dm,in_vec_local, (void*)&l_in_arr); CHKERRQ(ierr);
        ierr = DMDAVecGetArray(dm,out_vec_global, (void*)&g_out_arr); CHKERRQ(ierr);

        for (PetscInt k = info.zs; k < info.zs + info.zm; k++) {
            if (k == 0 || k == info.mz - 1) continue;
            for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
                if (j == 0 || j == info.my - 1) continue;
                for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
                    if (i == 0 || i == info.mx - 1) continue;

                    g_out_arr[k][j][i].x = 0.125 * (l_in_arr[k-1][j-1][i-1].x + l_in_arr[k-1][j-1][i].x +
                                                    l_in_arr[k-1][j][i-1].x   + l_in_arr[k-1][j][i].x   +
                                                    l_in_arr[k][j-1][i-1].x   + l_in_arr[k][j-1][i].x   +
                                                    l_in_arr[k][j][i-1].x     + l_in_arr[k][j][i].x);
                    g_out_arr[k][j][i].y = 0.125 * (l_in_arr[k-1][j-1][i-1].y + l_in_arr[k-1][j-1][i].y +
                                                    l_in_arr[k-1][j][i-1].y   + l_in_arr[k-1][j][i].y   +
                                                    l_in_arr[k][j-1][i-1].y   + l_in_arr[k][j-1][i].y   +
                                                    l_in_arr[k][j][i-1].y     + l_in_arr[k][j][i].y);
                    g_out_arr[k][j][i].z = 0.125 * (l_in_arr[k-1][j-1][i-1].z + l_in_arr[k-1][j-1][i].z +
                                                    l_in_arr[k-1][j][i-1].z   + l_in_arr[k-1][j][i].z   +
                                                    l_in_arr[k][j-1][i-1].z   + l_in_arr[k][j-1][i].z   +
                                                    l_in_arr[k][j][i-1].z     + l_in_arr[k][j][i].z);
                }
            }
        }
        ierr = DMDAVecRestoreArrayRead(dm,in_vec_local, (void*)&l_in_arr); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(dm,out_vec_global, (void*)&g_out_arr); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

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

    PetscFunctionReturn(0);
}