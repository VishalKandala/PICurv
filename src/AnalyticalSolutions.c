/**
 * @file AnalyticalSolutions.c
 * @brief Implements the analytical solution engine for initializing or driving the simulation.
 *
 * @details This file provides a modular and extensible framework for applying analytical solutions
 * to the Eulerian fields. The primary entry point is `AnalyticalSolutionEngine`, which acts
 * as a dispatcher based on user configuration.
 *
 * --- DESIGN PHILOSOPHY ---
 * 1.  **Non-Dimensional Core:** All calculations within this engine are performed in
 *     **non-dimensional units** (e.g., reference velocity U_ref=1.0, reference length L_ref=1.0).
 *     This is critical for consistency with the core numerical solver, which also operates on
 *     non-dimensional equations. The `simCtx->scaling` parameters are intentionally NOT used here;
 *     they are reserved for dimensionalization during I/O and post-processing only.
 *
 * 2.  **Separation of Concerns:** The role of this engine is to set the **physical state** of the
 *     fluid at a given time `t`. This involves:
 *     - Setting the values of fields (`Ucat`, `P`) in the physical interior of the domain.
 *     - Declaring the physical values on the boundaries by populating the boundary condition
 *       vector (`user->Bcs.Ubcs`).
 *     It does NOT implement the numerical scheme for ghost cells. Instead, after setting the
 *     physical state, it relies on the solver's standard utility functions (`UpdateDummyCells`,
 *     `UpdateCornerNodes`) to correctly populate all ghost cell layers.
 *
 * 3.  **Extensibility:** The dispatcher design makes it straightforward to add new analytical
 *     solutions. A developer only needs to add a new `else if` condition and a corresponding
 *     `static` implementation function, without modifying any other part of the solver.
 */

#include "AnalyticalSolutions.h"     // The header that declares this file's functions
#include <petscmath.h>               // Provides PETSc-compatible math functions like sin, cos, exp

// Forward-declare the private function for the TGV3D case.
// This function is not visible outside this file, enforcing modularity.
static PetscErrorCode SetAnalyticalSolution_TGV3D(SimCtx *simCtx);

#undef __FUNCT__
#define __FUNCT__ "SetAnalyticalGridInfo"
/**
 * @brief Sets the grid domain and resolution for analytical solution cases.
 *
 * @details This function is called when `eulerianSource` is "analytical". It is responsible for
 * automatically configuring the grid based on the chosen `AnalyticalSolutionType`.
 *
 * @par TGV3D Multi-Block Decomposition
 * If the analytical solution is "TGV3D", this function automatically decomposes the
 * required `[0, 2*PI]` physical domain among the available blocks.
 * - **Single Block (`nblk=1`):** The single block is assigned the full `[0, 2*PI]` domain.
 * - **Multiple Blocks (`nblk>1`):** It requires that the number of blocks be a **perfect square**
 *   (e.g., 4, 9, 16). It then arranges the blocks in a `sqrt(nblk)` by `sqrt(nblk)` grid in the
 *   X-Y plane, partitioning the `[0, 2*PI]` domain in X and Y accordingly. The Z domain for all
 *   blocks remains `[0, 2*PI]`. If `nblk` is not a perfect square, the simulation is aborted
 *   with an error.
 *
 * After setting the domain bounds, it proceeds to read the grid resolution options
 * (`-im`, `-jm`, `-km`) from the command line for the specific block.
 *
 * @param user Pointer to the `UserCtx` for a specific block. The function will
 *             populate the geometric fields (`IM`, `JM`, `KM`, `Min_X`, `Max_X`, etc.)
 *             within this struct.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode SetAnalyticalGridInfo(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx         *simCtx = user->simCtx;
    PetscInt       nblk = simCtx->block_number;
    PetscInt       block_index = user->_this;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    if (strcmp(simCtx->AnalyticalSolutionType, "TGV3D") == 0) {
        LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "Rank %d: Configuring grid for TGV3D analytical solution, block %d.\n", simCtx->rank, block_index);

        if (nblk == 1) {
            // --- Single Block Case ---
            if (block_index == 0) {
                LOG_ALLOW(GLOBAL, LOG_INFO, "Single block detected. Setting domain to [0, 2*PI].\n");
            }
            user->Min_X = 0.0; user->Max_X = 2.0 * PETSC_PI;
            user->Min_Y = 0.0; user->Max_Y = 2.0 * PETSC_PI;
            user->Min_Z = 0.0; user->Max_Z = 1.5; //2.0 * PETSC_PI;

        } else { // --- Multi-Block Case ---
            PetscReal s = sqrt((PetscReal)nblk);
            
            // Validate that nblk is a perfect square.
            if (fabs(s - floor(s)) > 1e-9) {
                SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
                         "\n\n*** CONFIGURATION ERROR FOR TGV3D ***\n"
                         "For multi-block TGV3D cases, the number of blocks must be a perfect square (e.g., 4, 9, 16).\n"
                         "You have specified %d blocks. Please adjust `-block_number`.\n", nblk);
            }
            PetscInt blocks_per_dim = (PetscInt)s;

            if (block_index == 0) {
                 LOG_ALLOW(GLOBAL, LOG_INFO, "%d blocks detected. Decomposing domain into a %d x %d grid in the X-Y plane.\n", nblk, blocks_per_dim, blocks_per_dim);
            }

            // Determine the (row, col) position of this block in the 2D decomposition
            PetscInt row = block_index / blocks_per_dim;
            PetscInt col = block_index % blocks_per_dim;

            // Calculate the width/height of each sub-domain
            PetscReal block_width  = (2.0 * PETSC_PI) / (PetscReal)blocks_per_dim;
            PetscReal block_height = (2.0 * PETSC_PI) / (PetscReal)blocks_per_dim;
            
            // Assign this block its specific sub-domain
            user->Min_X = col * block_width;
            user->Max_X = (col + 1) * block_width;
            user->Min_Y = row * block_height;
            user->Max_Y = (row + 1) * block_height;
            user->Min_Z = 0.0;
            user->Max_Z = 2.0 * PETSC_PI; // Z-domain is not decomposed
        }
    }
    /*
     * --- EXTENSIBILITY HOOK ---
     * To add another analytical case with special grid requirements:
     *
     * else if (strcmp(simCtx->AnalyticalSolutionType, "ChannelFlow") == 0) {
     *     // ... implement logic to set domain for ChannelFlow case ...
     * }
     */
    else {
        // Fallback for other analytical solutions: require manual grid specification.
        LOG(GLOBAL, LOG_INFO, "Analytical solution '%s' detected. Using user-provided grid generation inputs.\n", simCtx->AnalyticalSolutionType);
        ierr = ReadGridGenerationInputs(user); CHKERRQ(ierr);
        PetscFunctionReturn(0); // Return early as ReadGridGenerationInputs already handled everything.
    }
    
    // --- For TGV3D, we now read the grid resolution, as this is independent of the domain ---
    PetscInt  *IMs, *JMs, *KMs;
    PetscBool found;
    PetscInt  count;
    
    ierr = PetscMalloc3(nblk, &IMs, nblk, &JMs, nblk, &KMs); CHKERRQ(ierr);
    
    // Set defaults first
    for(PetscInt i=0; i<nblk; ++i) { IMs[i] = 10; JMs[i] = 10; KMs[i] = 10; }
    
    count = nblk; ierr = PetscOptionsGetIntArray(NULL, NULL, "-im", IMs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetIntArray(NULL, NULL, "-jm", JMs, &count, &found); CHKERRQ(ierr);
    count = nblk; ierr = PetscOptionsGetIntArray(NULL, NULL, "-km", KMs, &count, &found); CHKERRQ(ierr);
    
    user->IM = IMs[block_index];
    user->JM = JMs[block_index];
    user->KM = KMs[block_index];

    // We can also read stretching ratios, as they are independent of the domain size
    // For simplicity, we assume uniform grid unless specified.
    user->rx = 1.0; user->ry = 1.0; user->rz = 1.0;
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Block %d grid resolution set: IM=%d, JM=%d, KM=%d\n",
              simCtx->rank, block_index, user->IM, user->JM, user->KM);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Block %d final bounds: X=[%.4f, %.4f], Y=[%.4f, %.4f], Z=[%.4f, %.4f]\n",
              simCtx->rank, block_index, user->Min_X, user->Max_X, user->Min_Y, user->Max_Y, user->Min_Z, user->Max_Z);

    ierr = PetscFree3(IMs, JMs, KMs); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


/*================================================================================*
 *                       PUBLIC ANALYTICAL SOLUTION ENGINE                        *
 *================================================================================*/

#undef __FUNCT__
#define __FUNCT__ "AnalyticalSolutionEngine"
/**
 * @brief Dispatches to the appropriate analytical solution function based on simulation settings.
 *
 * This function acts as a router. It reads the `AnalyticalSolutionType` string from the
 * simulation context and calls the corresponding private implementation function
 * (e.g., for Taylor-Green Vortex). This design keeps the main simulation code clean
 * and makes it easy to add new analytical test cases.
 *
 * @param[in] simCtx The main simulation context, containing all configuration and state.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode AnalyticalSolutionEngine(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // -- Before any operation, here is defensive test to ensure that the Corner->Center Interpolation method works
    //ierr = TestCornerToCenterInterpolation(&(simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1]->user[0])); 

    // --- Dispatch based on the string provided by the user ---
    if (strcmp(simCtx->AnalyticalSolutionType, "TGV3D") == 0) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Applying Analytical Solution: 3D Taylor-Green Vortex (TGV3D).\n");
        ierr = SetAnalyticalSolution_TGV3D(simCtx); CHKERRQ(ierr);
    }
    /*
     * --- EXTENSIBILITY HOOK ---
     * To add a new analytical solution (e.g., "ChannelFlow"):
     * 1. Add an `else if` block here:
     *
     *    else if (strcmp(simCtx->AnalyticalSolutionType, "ChannelFlow") == 0) {
     *        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Applying Analytical Solution: Channel Flow.\n");
     *        ierr = SetAnalyticalSolution_ChannelFlow(simCtx); CHKERRQ(ierr);
     *    }
     *
     * 2. Implement the static function `SetAnalyticalSolution_ChannelFlow(SimCtx *simCtx)`
     *    below, following the TGV3D pattern.
     */
    else {
        // If the type is unknown, raise a fatal error to prevent silent failures.
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_ARG_UNKNOWN_TYPE, "Unknown AnalyticalSolutionType specified: '%s'", simCtx->AnalyticalSolutionType);
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


/*================================================================================*
 *                       PRIVATE IMPLEMENTATIONS                                  *
 *================================================================================*/

#undef __FUNCT__
#define __FUNCT__ "SetAnalyticalSolution_TGV3D"
/**
 * @brief Sets the non-dimensional velocity and pressure fields to the 3D Taylor-Green Vortex solution.
 *
 * @details This function implements the classical decaying Taylor-Green vortex in non-dimensional
 * form, consistent with the solver's internal state. It is fully compatible with the
 * solver's curvilinear, parallel, cell-centered architecture.
 *
 * @par WORKFLOW
 * 1.  Defines non-dimensional parameters for the TGV solution (V0=1.0, rho=1.0).
 * 2.  Defines correct loop bounds (`lxs`, `lxe`, etc.) to iterate over ONLY the physical
 *     interior cells owned by the current MPI rank.
 * 3.  Sets the **interior** cell-centered velocity (`Ucat`) and pressure (`P`) by evaluating
 *     the analytical formula at the cell-center coordinates (read from the local `user->Cent` vector).
 * 4.  Sets the **boundary condition** vector (`user->Bcs.Ubcs`) by evaluating the analytical
 *     velocity at the true physical face-center locations (read from the local `user->Centx`,
 *     `user->Centy`, and `user->Centz` vectors).
 * 5.  Manually sets the **pressure ghost cells** on the faces (excluding edges/corners) to
 *     enforce a zero-gradient (Neumann) boundary condition, which is standard for pressure.
 * 6.  Calls the solver's utility functions `UpdateDummyCells` (for `Ucat`) and `UpdateCornerNodes`
 *     (for both `Ucat` and `P`) to finalize all ghost cell values, ensuring a fully consistent field.
 *
 * @param[in] simCtx The main simulation context.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
static PetscErrorCode SetAnalyticalSolution_TGV3D(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserCtx       *user_finest = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;
    
    // --- NON-DIMENSIONAL TGV Parameters ---
    const PetscReal V0  = 1.0; // Non-dimensional reference velocity.
    const PetscReal rho = 1.0; // Non-dimensional reference density.
    const PetscReal p0  = 0.0; // Non-dimensional reference pressure.
    
    // Kinematic viscosity is derived from the non-dimensional Reynolds number.
    const PetscReal nu  = (simCtx->ren > 0) ? (1.0 / simCtx->ren) : 0.0; 

    const PetscReal k   = 1.0; // Wavenumber, assumes a non-dimensional [0, 2*pi] domain.
    const PetscReal t   = simCtx->ti;

    LOG_ALLOW(GLOBAL,LOG_TRACE,"TGV Setup: t = %.4f, V0* = %.4f, rho* = %.4f, k = %.4f, p0* = %4.f, nu = %.6f.\n",simCtx->ti,V0,rho,k,p0,nu);

    const PetscReal vel_decay = exp(-2.0 * nu * k * k * t);
    const PetscReal prs_decay = exp(-4.0 * nu * k * k * t);

    PetscFunctionBeginUser;

    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        UserCtx* user = &user_finest[bi];
        DMDALocalInfo info = user->info;
        PetscInt      xs = info.xs, xe = info.xs + info.xm;
        PetscInt      ys = info.ys, ye = info.ys + info.ym;
        PetscInt      zs = info.zs, ze = info.zs + info.zm;
        PetscInt      mx = info.mx, my = info.my, mz = info.mz;

        Cmpnts      ***ucat, ***ubcs;
        const Cmpnts ***cent, ***cent_x, ***cent_y, ***cent_z;
        PetscReal   ***p;
        
        // Define loop bounds for physical interior cells owned by this rank.
        PetscInt lxs = (xs == 0) ? xs + 1 : xs, lxe = (xe == mx) ? xe - 1 : xe;
        PetscInt lys = (ys == 0) ? ys + 1 : ys, lye = (ye == my) ? ye - 1 : ye;
        PetscInt lzs = (zs == 0) ? zs + 1 : zs, lze = (ze == mz) ? ze - 1 : ze;

        // --- Get Arrays ---
        ierr = DMDAVecGetArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);
        ierr = DMDAVecGetArray(user->da, user->P, &p); CHKERRQ(ierr);
        ierr = DMDAVecGetArray(user->fda, user->Bcs.Ubcs, &ubcs); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(user->fda, user->Cent, &cent); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(user->fda, user->Centx, &cent_x); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(user->fda, user->Centy, &cent_y); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(user->fda, user->Centz, &cent_z); CHKERRQ(ierr);

        // --- Set INTERIOR cell-centered velocity (Ucat) ---
        for (PetscInt k_cell = lzs; k_cell < lze; k_cell++) {
            for (PetscInt j_cell = lys; j_cell < lye; j_cell++) {
                for (PetscInt i_cell = lxs; i_cell < lxe; i_cell++) {
                    const PetscReal cx = cent[k_cell][j_cell][i_cell].x, cy = cent[k_cell][j_cell][i_cell].y, cz = cent[k_cell][j_cell][i_cell].z;
                    ucat[k_cell][j_cell][i_cell].x =  V0 * sin(k*cx) * cos(k*cy) * cos(k*cz) * vel_decay;
                    ucat[k_cell][j_cell][i_cell].y = -V0 * cos(k*cx) * sin(k*cy) * cos(k*cz) * vel_decay;
                    ucat[k_cell][j_cell][i_cell].z =  0.0;
                }
            }
        }


        // --- Set INTERIOR cell-centered pressure (P) ---
        for (PetscInt k_cell = lzs; k_cell < lze; k_cell++) {
            for (PetscInt j_cell = lys; j_cell < lye; j_cell++) {
                for (PetscInt i_cell = lxs; i_cell < lxe; i_cell++) {
                    const PetscReal cx = cent[k_cell][j_cell][i_cell].x, cy = cent[k_cell][j_cell][i_cell].y;
                    p[k_cell][j_cell][i_cell] = p0 + (rho * V0 * V0 / 4.0) * (cos(2*k*cx) + cos(2*k*cy)) * prs_decay;
                }
            }
        }
        

        // --- Set BOUNDARY condition vector for velocity (Ubcs) ---
        if (xs == 0) for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) {
            const PetscReal fcx=cent_x[k][j][xs].x, fcy=cent_x[k][j][xs].y, fcz=cent_x[k][j][xs].z;
            ubcs[k][j][xs].x = V0*sin(k*fcx)*cos(k*fcy)*cos(k*fcz)*vel_decay; ubcs[k][j][xs].y = -V0*cos(k*fcx)*sin(k*fcy)*cos(k*fcz)*vel_decay; ubcs[k][j][xs].z = 0.0;
        }
        if (xe == mx) for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) {
            const PetscReal fcx=cent_x[k][j][xe-1].x, fcy=cent_x[k][j][xe-1].y, fcz=cent_x[k][j][xe-1].z;
            ubcs[k][j][xe-1].x = V0*sin(k*fcx)*cos(k*fcy)*cos(k*fcz)*vel_decay; ubcs[k][j][xe-1].y = -V0*cos(k*fcx)*sin(k*fcy)*cos(k*fcz)*vel_decay; ubcs[k][j][xe-1].z = 0.0;
        }
        if (ys == 0) for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) {
            const PetscReal fcx=cent_y[k][ys][i].x, fcy=cent_y[k][ys][i].y, fcz=cent_y[k][ys][i].z;
            ubcs[k][ys][i].x = V0*sin(k*fcx)*cos(k*fcy)*cos(k*fcz)*vel_decay; ubcs[k][ys][i].y = -V0*cos(k*fcx)*sin(k*fcy)*cos(k*fcz)*vel_decay; ubcs[k][ys][i].z = 0.0;
        }
        if (ye == my) for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) {
            const PetscReal fcx=cent_y[k][ye-1][i].x, fcy=cent_y[k][ye-1][i].y, fcz=cent_y[k][ye-1][i].z;
            ubcs[k][ye-1][i].x = V0*sin(k*fcx)*cos(k*fcy)*cos(k*fcz)*vel_decay; ubcs[k][ye-1][i].y = -V0*cos(k*fcx)*sin(k*fcy)*cos(k*fcz)*vel_decay; ubcs[k][ye-1][i].z = 0.0;
        }
        if (zs == 0) for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) {
            const PetscReal fcx=cent_z[zs][j][i].x, fcy=cent_z[zs][j][i].y, fcz=cent_z[zs][j][i].z;
            ubcs[zs][j][i].x = V0*sin(k*fcx)*cos(k*fcy)*cos(k*fcz)*vel_decay; ubcs[zs][j][i].y = -V0*cos(k*fcx)*sin(k*fcy)*cos(k*fcz)*vel_decay; ubcs[zs][j][i].z = 0.0;
        }
        if (ze == mz) for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) {
            const PetscReal fcx=cent_z[ze-1][j][i].x, fcy=cent_z[ze-1][j][i].y, fcz=cent_z[ze-1][j][i].z;
            ubcs[ze-1][j][i].x = V0*sin(k*fcx)*cos(k*fcy)*cos(k*fcz)*vel_decay; ubcs[ze-1][j][i].y = -V0*cos(k*fcx)*sin(k*fcy)*cos(k*fcz)*vel_decay; ubcs[ze-1][j][i].z = 0.0;
        }

        // --- Set PRESSURE GHOST CELLS (Neumann BC: P_ghost = P_interior) ---
        if (xs == 0) for (PetscInt k=lzs; k<lze; k++) for (PetscInt j=lys; j<lye; j++) p[k][j][xs] = p[k][j][xs+1];
        if (xe == mx) for (PetscInt k=lzs; k<lze; k++) for (PetscInt j=lys; j<lye; j++) p[k][j][xe-1] = p[k][j][xe-2];
        
        if (ys == 0) for (PetscInt k=lzs; k<lze; k++) for (PetscInt i=lxs; i<lxe; i++) p[k][ys][i] = p[k][ys+1][i];
        if (ye == my) for (PetscInt k=lzs; k<lze; k++) for (PetscInt i=lxs; i<lxe; i++) p[k][ye-1][i] = p[k][ye-2][i];

        if (zs == 0) for (PetscInt j=lys; j<lye; j++) for (PetscInt i=lxs; i<lxe; i++) p[zs][j][i] = p[zs+1][j][i];
        if (ze == mz) for (PetscInt j=lys; j<lye; j++) for (PetscInt i=lxs; i<lxe; i++) p[ze-1][j][i] = p[ze-2][j][i];

        // --- Restore all arrays ---
        ierr = DMDAVecRestoreArray(user->fda, user->Ucat, &ucat); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(user->da, user->P, &p); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(user->fda, user->Bcs.Ubcs, &ubcs); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(user->fda, user->Cent, &cent); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(user->fda, user->Centx, &cent_x); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(user->fda, user->Centy, &cent_y); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(user->fda, user->Centz, &cent_z); CHKERRQ(ierr);

        // Pre-Dummy cell update synchronization.
        ierr = UpdateLocalGhosts(user,"Ucat");
        ierr = UpdateLocalGhosts(user,"P");

        // --- Finalize all ghost cell values ---
        ierr = UpdateDummyCells(user); CHKERRQ(ierr);
        ierr = UpdateCornerNodes(user); CHKERRQ(ierr);

        // Final Synchronization.
        ierr = UpdateLocalGhosts(user,"Ucat");
        ierr = UpdateLocalGhosts(user,"P");

    }

    PetscFunctionReturn(0);
}