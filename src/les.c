/**
 * @file les.c
 * @brief Implements the Large Eddy Simulation (LES) turbulence models.
 *
 * This file contains the core logic for the dynamic Smagorinsky LES model.
 * It is responsible for two primary tasks, orchestrated by the main FlowSolver:
 * 1.  `ComputeSmagorinskyConstant`: Dynamically calculates the Smagorinsky coefficient (Cs)
 *     at each grid point by applying a test filter to the resolved flow field. This
 *     is computationally intensive and is typically not performed every time step.
 * 2.  `ComputeEddyViscosityLES`: Uses the calculated Cs to compute the turbulent
 *     eddy viscosity (nu_t), which is then added to the molecular viscosity in the
 *     momentum equations to model the dissipative effects of sub-grid scale turbulence.
 */

#include "les.h"

// A small constant to prevent division by zero in sensitive calculations.
const double LES_EPSILON = 1.0e-12;


#undef __FUNCT__
#define __FUNCT__ "ComputeSmagorinskyConstant"
PetscErrorCode ComputeSmagorinskyConstant(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx *simCtx = user->simCtx; // Get the global context for simulation parameters

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // --- Initial Condition / Simple Model Handling ---

    // For the first couple of steps of a simulation from t=0, the flow is not developed
    // enough for the dynamic procedure to be stable. Set Cs to zero.
    if (simCtx->step < 2 && simCtx->StartStep == 0) {
        ierr = VecSet(user->CS, 0.0); CHKERRQ(ierr);
        LOG_ALLOW(GLOBAL,LOG_DEBUG,"Setting Smagorinsky coefficient Cs=0.0 for initial steps (step=%d)\n", simCtx->step);
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    // If the user requests the non-dynamic, constant-coefficient Smagorinsky model (les=1),
    // set a constant value and exit.
    if (simCtx->les == 1) {
        LOG_ALLOW(GLOBAL,LOG_INFO,"Using constant-coefficient Smagorinsky model with Cs=%.4f \n", simCtx->Const_CS);
        ierr = VecSet(user->CS, simCtx->Const_CS); CHKERRQ(ierr); // A typical constant value
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    // --- Setup for Dynamic Procedure ---

	DM		        da = user->da, fda = user->fda;
	DMDALocalInfo	info;
	PetscInt	    i, j, k, p, q, r;
	PetscInt	    xs, xe, ys, ye, zs, ze; // Local grid corner indices
	PetscInt	    mx, my, mz;             // Global grid dimensions
	PetscInt	    lxs, lxe, lys, lye, lzs, lze; // Local interior loop bounds

    // Temporary PETSc vectors required for the calculation. These are allocated
    // here and destroyed at the end to keep the UserCtx clean.
    Vec lSx, lSy, lSz, lS_abs;
    Vec lLM, lMM; // Leonard & cross-term stress tensors

    // Get grid dimensions and local bounds
	ierr = DMDAGetLocalInfo(da, &info); CHKERRQ(ierr);
	mx = info.mx; my = info.my; mz = info.mz;
	xs = info.xs; xe = xs + info.xm;
	ys = info.ys; ye = ys + info.ym;
	zs = info.zs; ze = zs + info.zm;

    // Define loop bounds to exclude physical boundaries where stencils are invalid.
	lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;
	if (xs==0) lxs = xs+1;
	if (ys==0) lys = ys+1;
	if (zs==0) lzs = zs+1;
	if (xe==mx) lxe = xe-1;
	if (ye==my) lye = ye-1;
	if (ze==mz) lze = ze-1;

    // Allocate temporary vectors
    ierr = VecDuplicate(user->lUcont, &lSx); CHKERRQ(ierr);
    ierr = VecDuplicate(user->lUcont, &lSy); CHKERRQ(ierr);
    ierr = VecDuplicate(user->lUcont, &lSz); CHKERRQ(ierr);
    ierr = VecDuplicate(user->lP, &lS_abs); CHKERRQ(ierr);
    ierr = VecDuplicate(user->lP, &lLM); CHKERRQ(ierr);
    ierr = VecDuplicate(user->lP, &lMM); CHKERRQ(ierr);
    ierr = VecSet(lLM, 0.0); CHKERRQ(ierr);
    ierr = VecSet(lMM, 0.0); CHKERRQ(ierr);

    // Get read/write access to PETSc vector data as multidimensional arrays
    Cmpnts	    ***ucat, ***csi, ***eta, ***zet;
    PetscReal	***nvert, ***Cs_arr, ***aj, ***S_abs_arr, ***LM_arr, ***MM_arr;
    Cmpnts      ***Sx_arr, ***Sy_arr, ***Sz_arr;

    ierr = DMDAVecGetArray(fda, user->lUcat, &ucat); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(fda, user->lCsi, (Cmpnts***)&csi); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(fda, user->lEta, (Cmpnts***)&eta); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(fda, user->lZet, (Cmpnts***)&zet); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da, user->lNvert, (PetscReal***)&nvert); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da, user->lAj, (PetscReal***)&aj); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, user->CS, &Cs_arr); CHKERRQ(ierr);

    ierr = DMDAVecGetArray(fda, lSx, &Sx_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(fda, lSy, &Sy_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(fda, lSz, &Sz_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, lS_abs, &S_abs_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, lLM, &LM_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, lMM, &MM_arr); CHKERRQ(ierr);

    // --- 1. Compute and store the strain rate tensor |S| and velocity gradients at all points ---
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if( nvert[k][j][i] > 1.1) continue; // Skip points inside solid bodies

        Cmpnts dudx, dvdx, dwdx;
        ierr = ComputeVectorFieldDerivatives(user, i, j, k, ucat, &dudx, &dvdx, &dwdx); CHKERRQ(ierr);

		double Sxx = dudx.x;
		double Sxy = 0.5 * (dudx.y + dvdx.x);
		double Sxz = 0.5 * (dudx.z + dwdx.x);
		double Syy = dvdx.y;
		double Syz = 0.5 * (dvdx.z + dwdx.y);
		double Szz = dwdx.z;
		double Syx = Sxy, Szx = Sxz, Szy = Syz; // Symmetry

		S_abs_arr[k][j][i] = sqrt( 2.0 * (Sxx*Sxx + Sxy*Sxy + Sxz*Sxz +
                                         Syx*Syx + Syy*Syy + Syz*Syz +
                                         Szx*Szx + Szy*Szy + Szz*Szz) );

        // Store the full velocity gradient tensor for use in the next step
        Sx_arr[k][j][i] = dudx; // Contains {du/dx, du/dy, du/dz}
        Sy_arr[k][j][i] = dvdx; // Contains {dv/dx, dv/dy, dv/dz}
        Sz_arr[k][j][i] = dwdx; // Contains {dw/dx, dw/dy, dw/dz}
	}

    // --- 2. Main loop to compute the dynamic coefficient ---
 	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if(nvert[k][j][i] > 1.1) {
			LM_arr[k][j][i] = 0.0;
            MM_arr[k][j][i] = 0.0;
			continue;
		}

        // --- 2a. Gather data from the 3x3x3 stencil around the current point (i,j,k) ---
		double u[3][3][3], v[3][3][3], w[3][3][3];
		double S_abs_stencil[3][3][3];
		double S11[3][3][3], S12[3][3][3], S13[3][3][3];
        double S22[3][3][3], S23[3][3][3], S33[3][3][3];
		double weights[3][3][3];

		for(r=-1; r<=1; r++)
		for(q=-1; q<=1; q++)
		for(p=-1; p<=1; p++) {
			int R=r+1, Q=q+1, P=p+1; // Stencil array indices (0-2)
			int KK=k+r, JJ=j+q, II=i+p; // Global grid indices

			u[R][Q][P] = ucat[KK][JJ][II].x;
			v[R][Q][P] = ucat[KK][JJ][II].y;
			w[R][Q][P] = ucat[KK][JJ][II].z;

            // Strain rate tensor components (Sij = 0.5 * (dui/dxj + duj/dxi))
            S11[R][Q][P] = Sx_arr[KK][JJ][II].x;                      // du/dx
            S12[R][Q][P] = 0.5 * (Sx_arr[KK][JJ][II].y + Sy_arr[KK][JJ][II].x); // 0.5 * (du/dy + dv/dx)
            S13[R][Q][P] = 0.5 * (Sx_arr[KK][JJ][II].z + Sz_arr[KK][JJ][II].x); // 0.5 * (du/dz + dw/dx)
            S22[R][Q][P] = Sy_arr[KK][JJ][II].y;                      // dv/dy
            S23[R][Q][P] = 0.5 * (Sy_arr[KK][JJ][II].z + Sz_arr[KK][JJ][II].y); // 0.5 * (dv/dz + dw/dy)
            S33[R][Q][P] = Sz_arr[KK][JJ][II].z;                      // dw/dz

			S_abs_stencil[R][Q][P] = S_abs_arr[KK][JJ][II];
			weights[R][Q][P] = aj[KK][JJ][II]; // Weight is Jacobian (1/volume)
		}

        // --- 2b. Apply the test filter to all required quantities ---
        double u_filt  = ApplyLESTestFilter(simCtx, u, weights);
        double v_filt  = ApplyLESTestFilter(simCtx, v, weights);
        double w_filt  = ApplyLESTestFilter(simCtx, w, weights);

        double S11_filt = ApplyLESTestFilter(simCtx, S11, weights);
        double S12_filt = ApplyLESTestFilter(simCtx, S12, weights);
        double S13_filt = ApplyLESTestFilter(simCtx, S13, weights);
        double S22_filt = ApplyLESTestFilter(simCtx, S22, weights);
        double S23_filt = ApplyLESTestFilter(simCtx, S23, weights);
        double S33_filt = ApplyLESTestFilter(simCtx, S33, weights);

        double S_abs_filt = ApplyLESTestFilter(simCtx, S_abs_stencil, weights);

        // Filter quadratic terms: <u*u>, <u*v>, etc.
        double uu[3][3][3], uv[3][3][3], uw[3][3][3], vv[3][3][3], vw[3][3][3], ww[3][3][3];
        for(r=0; r<3; r++) for(q=0; q<3; q++) for(p=0; p<3; p++) {
            uu[r][q][p] = u[r][q][p] * u[r][q][p];
            uv[r][q][p] = u[r][q][p] * v[r][q][p];
            uw[r][q][p] = u[r][q][p] * w[r][q][p];
            vv[r][q][p] = v[r][q][p] * v[r][q][p];
            vw[r][q][p] = v[r][q][p] * w[r][q][p];
            ww[r][q][p] = w[r][q][p] * w[r][q][p];
        }
        double uu_filt = ApplyLESTestFilter(simCtx, uu, weights);
        double uv_filt = ApplyLESTestFilter(simCtx, uv, weights);
        double uw_filt = ApplyLESTestFilter(simCtx, uw, weights);
        double vv_filt = ApplyLESTestFilter(simCtx, vv, weights);
        double vw_filt = ApplyLESTestFilter(simCtx, vw, weights);
        double ww_filt = ApplyLESTestFilter(simCtx, ww, weights);

        // --- 2c. Compute the Leonard stress tensor (L_ij = <u_i*u_j> - <u_i>*<u_j>) ---
		double L11 = uu_filt - u_filt * u_filt;
		double L12 = uv_filt - u_filt * v_filt;
		double L13 = uw_filt - u_filt * w_filt;
		double L22 = vv_filt - v_filt * v_filt;
		double L23 = vw_filt - v_filt * w_filt;
		double L33 = ww_filt - w_filt * w_filt;

        // --- 2d. Compute the model tensor M_ij ---
        double grid_filter_width, test_filter_width;
        ierr = ComputeCellCharacteristicLengthScale(aj[k][j][i], csi[k][j][i], eta[k][j][i], zet[k][j][i], &grid_filter_width, &test_filter_width, &test_filter_width); CHKERRQ(ierr); // Simplified for now
        grid_filter_width = pow(1.0 / aj[k][j][i], 1.0/3.0);
        test_filter_width = 2.0 * grid_filter_width; // Standard box filter definition

        double alpha = pow(test_filter_width / grid_filter_width, 2.0);

        double M11 = -2.0 * grid_filter_width * grid_filter_width * (alpha * S_abs_filt * S11_filt - S_abs_filt * S11_filt);
        double M12 = -2.0 * grid_filter_width * grid_filter_width * (alpha * S_abs_filt * S12_filt - S_abs_filt * S12_filt);
        double M13 = -2.0 * grid_filter_width * grid_filter_width * (alpha * S_abs_filt * S13_filt - S_abs_filt * S13_filt);
        double M22 = -2.0 * grid_filter_width * grid_filter_width * (alpha * S_abs_filt * S22_filt - S_abs_filt * S22_filt);
        double M23 = -2.0 * grid_filter_width * grid_filter_width * (alpha * S_abs_filt * S23_filt - S_abs_filt * S23_filt);
        double M33 = -2.0 * grid_filter_width * grid_filter_width * (alpha * S_abs_filt * S33_filt - S_abs_filt * S33_filt);

        // --- 2e. Contract tensors to find Cs^2 (L_ij * M_ij / (M_kl * M_kl)) ---
		LM_arr[k][j][i] = L11*M11 + 2*L12*M12 + 2*L13*M13 + L22*M22 + 2*L23*M23 + L33*M33;
		MM_arr[k][j][i] = M11*M11 + 2*M12*M12 + 2*M13*M13 + M22*M22 + 2*M23*M23 + M33*M33;
    }

    // --- 3. Averaging and finalization ---

    // At this point, LM_arr and MM_arr contain raw, unaveraged values.
    // To stabilize the solution, these are often averaged over homogeneous directions or locally.
    // The logic for homogeneous averaging would go here if simCtx->i_homo_filter is true.
    // For this general version, we proceed with the local (unaveraged) values.

	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if(nvert[k][j][i] > 1.1) {
			Cs_arr[k][j][i] = 0.0;
			continue;
		}

        double C_sq = 0.0;
        if (MM_arr[k][j][i] > LES_EPSILON) {
            C_sq = LM_arr[k][j][i] / MM_arr[k][j][i];
        }

        // Final clipping and assignment: Cs must be positive and not excessively large.
        double Cs_val = (C_sq > 0.0) ? sqrt(C_sq) : 0.0;
		Cs_arr[k][j][i] = PetscMin(PetscMax(Cs_val, 0.0), simCtx->max_cs);
	}

    // --- 4. Cleanup ---

    // Restore all PETSc arrays
    ierr = DMDAVecRestoreArray(fda, user->lUcat, &ucat); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(fda, user->lCsi, (Cmpnts***)&csi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(fda, user->lEta, (Cmpnts***)&eta); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(fda, user->lZet, (Cmpnts***)&zet); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da, user->lNvert, (PetscReal***)&nvert); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(da, user->lAj, (PetscReal***)&aj); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da, user->CS, &Cs_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(fda, lSx, &Sx_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(fda, lSy, &Sy_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(fda, lSz, &Sz_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da, lS_abs, &S_abs_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da, lLM, &LM_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da, lMM, &MM_arr); CHKERRQ(ierr);

    // Destroy temporary vectors to prevent memory leaks
	ierr = VecDestroy(&lSx); CHKERRQ(ierr);
	ierr = VecDestroy(&lSy); CHKERRQ(ierr);
	ierr = VecDestroy(&lSz); CHKERRQ(ierr);
	ierr = VecDestroy(&lS_abs); CHKERRQ(ierr);
	ierr = VecDestroy(&lLM); CHKERRQ(ierr);
	ierr = VecDestroy(&lMM); CHKERRQ(ierr);

    // Communicate ghost point data for the newly computed Cs field
	ierr = DMGlobalToLocalBegin(da, user->CS, INSERT_VALUES, user->lCs); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da, user->CS, INSERT_VALUES, user->lCs); CHKERRQ(ierr);

    PetscReal max_norm;
    ierr = VecMax(user->CS, NULL, &max_norm); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "  Max dynamic Smagorinsky constant (Cs) computed: %e (clip value: %e)\n", max_norm, simCtx->max_cs);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ComputeEddyViscosityLES"
PetscErrorCode ComputeEddyViscosityLES(UserCtx *user)
{
    PetscErrorCode ierr;
	DM		        da = user->da, fda = user->fda;
	DMDALocalInfo	info;
	PetscInt	    i, j, k;
	PetscInt	    xs, xe, ys, ye, zs, ze;
    PetscInt        lxs, lxe, lys, lye, lzs, lze; // Local interior loop bounds
	PetscInt	    mx, my, mz;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // Get grid dimensions and local bounds
	ierr = DMDAGetLocalInfo(da, &info); CHKERRQ(ierr);
	mx = info.mx; my = info.my; mz = info.mz;
	xs = info.xs; xe = xs + info.xm;
	ys = info.ys; ye = ys + info.ym;
	zs = info.zs; ze = zs + info.zm;

    // Define loop bounds to exclude physical boundaries where stencils are invalid.
    lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;
    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;
    if (xe==mx) lxe = xe-1;
    if (ye==my) lye = ye-1;
    if (ze==mz) lze = ze-1;



    // Get read/write access to PETSc data arrays
    Cmpnts	    ***ucat;
    PetscReal	***Cs_arr, ***nu_t_arr, ***nvert, ***aj;

	ierr = DMDAVecGetArrayRead(fda, user->lUcat, (Cmpnts***)&ucat); CHKERRQ(ierr);
	ierr = DMDAVecGetArrayRead(da, user->lCs, (PetscReal***)&Cs_arr); CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, user->Nu_t, &nu_t_arr); CHKERRQ(ierr);
	ierr = DMDAVecGetArrayRead(da, user->lNvert, (PetscReal***)&nvert); CHKERRQ(ierr);
	ierr = DMDAVecGetArrayRead(da, user->lAj, (PetscReal***)&aj); CHKERRQ(ierr);

    // Loop over the interior of the local domain
	for (k=lzs; k<lze; k++)
	for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
		if(nvert[k][j][i] > 0.1) {
			nu_t_arr[k][j][i] = 0.0;
			continue;
		}

        LOG_ALLOW(GLOBAL, LOG_VERBOSE, "  Computing eddy viscosity at point (%d,%d,%d)\n", i, j, k);
        // 1. Compute the local strain rate magnitude |S|
        Cmpnts dudx, dvdx, dwdx;
        ierr = ComputeVectorFieldDerivatives(user, i, j, k, ucat, &dudx, &dvdx, &dwdx); CHKERRQ(ierr);

		double Sxx = dudx.x;
		double Sxy = 0.5 * (dudx.y + dvdx.x);
		double Sxz = 0.5 * (dudx.z + dwdx.x);
		double Syy = dvdx.y;
		double Syz = 0.5 * (dvdx.z + dwdx.y);
		double Szz = dwdx.z;
		double strain_rate_mag = sqrt( 2.0 * (Sxx*Sxx + Sxy*Sxy + Sxz*Sxz + Sxy*Sxy + Syy*Syy + Syz*Syz + Sxz*Sxz + Syz*Syz + Szz*Szz) );

        // 2. Determine the grid filter width, Delta = (cell volume)^(1/3)
		double filter_width  = pow( 1.0/aj[k][j][i], 1.0/3.0 );

        // 3. Compute eddy viscosity: nu_t = (Cs * Delta)^2 * |S|
		nu_t_arr[k][j][i] = pow(Cs_arr[k][j][i] * filter_width, 2.0) * strain_rate_mag;

        LOG_ALLOW(GLOBAL, LOG_VERBOSE, "    Cs=%.4e, Delta=%.4e, |S|=%.4e => nu_t=%.4e\n",
                    Cs_arr[k][j][i], filter_width, strain_rate_mag, nu_t_arr[k][j][i]);
	}

    // Restore PETSc data arrays
	ierr = DMDAVecRestoreArrayRead(fda, user->lUcat, (Cmpnts***)&ucat); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayRead(da, user->lCs, (PetscReal***)&Cs_arr); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, user->Nu_t, &nu_t_arr); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayRead(da, user->lNvert, (PetscReal***)&nvert); CHKERRQ(ierr);
	ierr = DMDAVecRestoreArrayRead(da, user->lAj, (PetscReal***)&aj); CHKERRQ(ierr);

    // Update ghost points for the newly computed eddy viscosity
	ierr = DMGlobalToLocalBegin(da, user->Nu_t, INSERT_VALUES, user->lNu_t); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da, user->Nu_t, INSERT_VALUES, user->lNu_t); CHKERRQ(ierr);

    PetscReal max_norm;
    ierr = VecMax(user->Nu_t, NULL, &max_norm); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "  Max eddy viscosity (Nu_t) computed: %e\n", max_norm);

    PROFILE_FUNCTION_END;
	PetscFunctionReturn(0);
}