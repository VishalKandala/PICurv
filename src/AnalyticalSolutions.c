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
 *     - For most types (TGV3D, ZERO_FLOW): setting `Ucat` and `P` directly.
 *     - For curvilinear-aware types (UNIFORM_FLOW): setting `Ucont` via metric
 *       dot products and deriving `Ucat` through `Contra2Cart`.
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
static PetscErrorCode SetAnalyticalSolution_ZeroFlow(SimCtx *simCtx);
static PetscErrorCode SetAnalyticalSolution_UniformFlow(SimCtx *simCtx);
static PetscErrorCode SetAnalyticalSolutionForParticles_TGV3D(Vec tempVec, SimCtx *simCtx);
static PetscErrorCode SetAnalyticalSolutionForParticles_UniformFlow(Vec tempVec, SimCtx *simCtx);
static PetscErrorCode EvaluateConfiguredScalarProfile(const VerificationScalarConfig *cfg,
                                                      PetscReal x,
                                                      PetscReal y,
                                                      PetscReal z,
                                                      PetscReal *value);

/**
 * @brief Implementation of \ref AnalyticalTypeRequiresCustomGeometry().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/AnalyticalSolutions.h`.
 * @see AnalyticalTypeRequiresCustomGeometry()
 */

PetscBool AnalyticalTypeRequiresCustomGeometry(const char *analytical_type)
{
    if (!analytical_type) return PETSC_FALSE;
    return (strcmp(analytical_type, "TGV3D") == 0) ? PETSC_TRUE : PETSC_FALSE;
}

/**
 * @brief Implementation of \ref AnalyticalTypeSupportsInterpolationError().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/AnalyticalSolutions.h`.
 * @see AnalyticalTypeSupportsInterpolationError()
 */
PetscBool AnalyticalTypeSupportsInterpolationError(const char *analytical_type)
{
    if (!analytical_type) return PETSC_FALSE;
    if (strcmp(analytical_type, "ZERO_FLOW") == 0 || strcmp(analytical_type, "UNIFORM_FLOW") == 0) return PETSC_FALSE;
    return PETSC_TRUE;
}

#undef __FUNCT__
#define __FUNCT__ "SetAnalyticalGridInfo"
/**
 * @brief Internal helper implementation: `SetAnalyticalGridInfo()`.
 * @details Local to this translation unit.
 */
PetscErrorCode SetAnalyticalGridInfo(UserCtx *user)
{
    SimCtx         *simCtx = user->simCtx;
    PetscInt       nblk = simCtx->block_number;
    PetscInt       block_index = user->_this;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    if (!AnalyticalTypeRequiresCustomGeometry(simCtx->AnalyticalSolutionType)) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE,
                 "SetAnalyticalGridInfo called for analytical type '%s' that does not require custom geometry.",
                 simCtx->AnalyticalSolutionType);
    }
    if (user->IM <= 0 || user->JM <= 0 || user->KM <= 0) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE,
                "Analytical grid resolution is not initialized. Ensure IM/JM/KM are preloaded before SetAnalyticalGridInfo.");
    }

    if (strcmp(simCtx->AnalyticalSolutionType, "TGV3D") == 0) {
        LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "Rank %d: Configuring grid for TGV3D analytical solution, block %d.\n", simCtx->rank, block_index);

        if (nblk == 1) {
            // --- Single Block Case ---
            if (block_index == 0) {
                LOG_ALLOW(GLOBAL, LOG_INFO, "Single block detected. Setting domain to [0, 2*PI].\n");
            }
            user->Min_X = 0.0; user->Max_X = 2.0 * PETSC_PI;
            user->Min_Y = 0.0; user->Max_Y = 2.0 * PETSC_PI;
            user->Min_Z = 0.0; user->Max_Z = 0.2 * PETSC_PI; //2.0 * PETSC_PI;

        } else { // --- Multi-Block Case ---
            PetscReal s = sqrt((PetscReal)nblk);
            
            // Validate that nblk is a perfect square.
            if (fabs(s - floor(s)) > 1e-9) {
                SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
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
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_UNKNOWN_TYPE,
                 "Analytical type '%s' has no custom geometry implementation.",
                 simCtx->AnalyticalSolutionType);
    }

    // We can also read stretching ratios, as they are independent of the domain size
    // For simplicity, we assume uniform grid unless specified.
    user->rx = 1.0; user->ry = 1.0; user->rz = 1.0;
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Block %d grid resolution set: IM=%d, JM=%d, KM=%d\n",
              simCtx->rank, block_index, user->IM, user->JM, user->KM);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Block %d final bounds: X=[%.4f, %.4f], Y=[%.4f, %.4f], Z=[%.4f, %.4f]\n",
              simCtx->rank, block_index, user->Min_X, user->Max_X, user->Min_Y, user->Max_Y, user->Min_Z, user->Max_Z);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


/*================================================================================*
 *                       PUBLIC ANALYTICAL SOLUTION ENGINE                        *
 *================================================================================*/

#undef __FUNCT__
#define __FUNCT__ "AnalyticalSolutionEngine"
/**
 * @brief Implementation of \ref AnalyticalSolutionEngine().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/AnalyticalSolutions.h`.
 * @see AnalyticalSolutionEngine()
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
    else if (strcmp(simCtx->AnalyticalSolutionType, "ZERO_FLOW") == 0) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Applying Analytical Solution: Zero Background Flow (ZERO_FLOW).\n");
        ierr = SetAnalyticalSolution_ZeroFlow(simCtx); CHKERRQ(ierr);
    }
    else if (strcmp(simCtx->AnalyticalSolutionType, "UNIFORM_FLOW") == 0) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Applying Analytical Solution: Uniform Background Flow (UNIFORM_FLOW).\n");
        ierr = SetAnalyticalSolution_UniformFlow(simCtx); CHKERRQ(ierr);
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
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_UNKNOWN_TYPE, "Unknown AnalyticalSolutionType specified: '%s'", simCtx->AnalyticalSolutionType);
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
 * @brief Internal helper implementation: `SetAnalyticalSolution_TGV3D()`.
 * @details Local to this translation unit.
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

#undef __FUNCT__
#define __FUNCT__ "SetAnalyticalSolution_ZeroFlow"
/**
 * @brief Internal helper implementation: `SetAnalyticalSolution_ZeroFlow()`.
 * @details Local to this translation unit.
 */
static PetscErrorCode SetAnalyticalSolution_ZeroFlow(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserCtx *user_finest = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;

    PetscFunctionBeginUser;

    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        UserCtx *user = &user_finest[bi];

        ierr = VecZeroEntries(user->Ucat);     CHKERRQ(ierr);
        ierr = VecZeroEntries(user->P);        CHKERRQ(ierr);
        ierr = VecZeroEntries(user->Bcs.Ubcs); CHKERRQ(ierr);

        // Ghost-cell finalization — identical sequence to TGV3D
        ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, "P");    CHKERRQ(ierr);
        ierr = UpdateDummyCells(user);          CHKERRQ(ierr);
        ierr = UpdateCornerNodes(user);         CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, "P");    CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetAnalyticalSolution_UniformFlow"
/**
 * @brief Internal helper implementation: `SetAnalyticalSolution_UniformFlow()`.
 * @details Local to this translation unit.
 */
static PetscErrorCode SetAnalyticalSolution_UniformFlow(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserCtx *user_finest = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;
    const Cmpnts uniform_velocity = simCtx->AnalyticalUniformVelocity;
    const PetscReal u = uniform_velocity.x;
    const PetscReal v = uniform_velocity.y;
    const PetscReal w = uniform_velocity.z;

    PetscFunctionBeginUser;

    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        UserCtx *user = &user_finest[bi];
        DMDALocalInfo info = user->info;
        PetscInt      xs = info.xs, xe = info.xs + info.xm;
        PetscInt      ys = info.ys, ye = info.ys + info.ym;
        PetscInt      zs = info.zs, ze = info.zs + info.zm;
        PetscInt      mx = info.mx, my = info.my, mz = info.mz;

        // --- Step 1: Set Ucont (contravariant flux) using metric dot products ---
        // U^i = g_i · u_phys, where g_i is the face area vector (Csi, Eta, Zet).
        // This is the curvilinear form: the volume flux through each face.
        // We loop over ALL owned nodes (including boundary faces) because Contra2Cart
        // averages adjacent face values — boundary faces must carry valid fluxes.
        Cmpnts ***ucont, ***ubcs;
        const Cmpnts ***csi, ***eta, ***zet;

        ierr = DMDAVecGetArray(user->fda, user->Ucont, &ucont);         CHKERRQ(ierr);
        ierr = DMDAVecGetArray(user->fda, user->Bcs.Ubcs, &ubcs);      CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(user->fda, user->lCsi, &csi);       CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(user->fda, user->lEta, &eta);       CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(user->fda, user->lZet, &zet);       CHKERRQ(ierr);

        for (PetscInt k = zs; k < ze; k++) {
            for (PetscInt j = ys; j < ye; j++) {
                for (PetscInt i = xs; i < xe; i++) {
                    ucont[k][j][i].x = csi[k][j][i].x * u + csi[k][j][i].y * v + csi[k][j][i].z * w;
                    ucont[k][j][i].y = eta[k][j][i].x * u + eta[k][j][i].y * v + eta[k][j][i].z * w;
                    ucont[k][j][i].z = zet[k][j][i].x * u + zet[k][j][i].y * v + zet[k][j][i].z * w;
                }
            }
        }

        // --- Step 2: Set Ubcs at boundaries (physical Cartesian velocity) ---
        if (xs == 0) for (PetscInt k = zs; k < ze; k++) for (PetscInt j = ys; j < ye; j++) ubcs[k][j][xs] = uniform_velocity;
        if (xe == mx) for (PetscInt k = zs; k < ze; k++) for (PetscInt j = ys; j < ye; j++) ubcs[k][j][xe - 1] = uniform_velocity;
        if (ys == 0) for (PetscInt k = zs; k < ze; k++) for (PetscInt i = xs; i < xe; i++) ubcs[k][ys][i] = uniform_velocity;
        if (ye == my) for (PetscInt k = zs; k < ze; k++) for (PetscInt i = xs; i < xe; i++) ubcs[k][ye - 1][i] = uniform_velocity;
        if (zs == 0) for (PetscInt j = ys; j < ye; j++) for (PetscInt i = xs; i < xe; i++) ubcs[zs][j][i] = uniform_velocity;
        if (ze == mz) for (PetscInt j = ys; j < ye; j++) for (PetscInt i = xs; i < xe; i++) ubcs[ze - 1][j][i] = uniform_velocity;

        ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet, &zet);   CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta, &eta);   CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi, &csi);   CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(user->fda, user->Bcs.Ubcs, &ubcs);  CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(user->fda, user->Ucont, &ucont);    CHKERRQ(ierr);

        // --- Step 3: Zero pressure ---
        ierr = VecZeroEntries(user->P); CHKERRQ(ierr);

        // --- Step 4: Finalize state — derive Ucat from Ucont via metric inversion ---
        ierr = UpdateLocalGhosts(user, "Ucont"); CHKERRQ(ierr);
        ierr = Contra2Cart(user);                CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, "Ucat");  CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, "P");     CHKERRQ(ierr);
        ierr = UpdateDummyCells(user);           CHKERRQ(ierr);
        ierr = UpdateCornerNodes(user);          CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, "Ucat");  CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, "P");     CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetAnalyticalSolutionForParticles_TGV3D"
/**
 * @brief Internal helper implementation: `SetAnalyticalSolutionForParticles_TGV3D()`.
 * @details Local to this translation unit.
 */
static PetscErrorCode SetAnalyticalSolutionForParticles_TGV3D(Vec tempVec, SimCtx *simCtx)
{
    PetscErrorCode ierr;
    PetscInt nLocal;
    PetscReal *data;
    
    PetscFunctionBeginUser;

    // TGV3D parameters (matching your Eulerian implementation)
    const PetscReal V0  = 1.0;
    const PetscReal k   = 1.0;
    const PetscReal nu  = (simCtx->ren > 0) ? (1.0 / simCtx->ren) : 0.0;
    const PetscReal t   = simCtx->ti;
    const PetscReal vel_decay = exp(-2.0 * nu * k * k * t);
    
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "TGV3D Particles: t=%.4f, V0=%.4f, k=%.4f, nu=%.6f\n", t, V0, k, nu);
    
    ierr = VecGetLocalSize(tempVec, &nLocal); CHKERRQ(ierr);
    ierr = VecGetArray(tempVec, &data); CHKERRQ(ierr);
    
    // Process particles: data is interleaved [x0,y0,z0, x1,y1,z1, ...]
    for (PetscInt i = 0; i < nLocal; i += 3) {
        const PetscReal x = data[i];
        const PetscReal y = data[i+1];
        const PetscReal z = data[i+2];
        
        // TGV3D velocity field
        data[i]   =  V0 * sin(k*x) * cos(k*y) * cos(k*z) * vel_decay;  // u
        data[i+1] = -V0 * cos(k*x) * sin(k*y) * cos(k*z) * vel_decay;  // v
        data[i+2] =  0.0;                                                // w
    }
    
    ierr = VecRestoreArray(tempVec, &data); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetAnalyticalSolutionForParticles_UniformFlow"
/**
 * @brief Internal helper implementation: `SetAnalyticalSolutionForParticles_UniformFlow()`.
 * @details Local to this translation unit.
 */
static PetscErrorCode SetAnalyticalSolutionForParticles_UniformFlow(Vec tempVec, SimCtx *simCtx)
{
    PetscErrorCode ierr;
    PetscInt nLocal;
    PetscReal *data;
    const Cmpnts uniform_velocity = simCtx->AnalyticalUniformVelocity;

    PetscFunctionBeginUser;

    ierr = VecGetLocalSize(tempVec, &nLocal); CHKERRQ(ierr);
    ierr = VecGetArray(tempVec, &data); CHKERRQ(ierr);
    for (PetscInt i = 0; i < nLocal; i += 3) {
        data[i] = uniform_velocity.x;
        data[i + 1] = uniform_velocity.y;
        data[i + 2] = uniform_velocity.z;
    }
    ierr = VecRestoreArray(tempVec, &data); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetAnalyticalSolutionForParticles"
/**
 * @brief Implementation of \ref SetAnalyticalSolutionForParticles().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/AnalyticalSolutions.h`.
 * @see SetAnalyticalSolutionForParticles()
 */
PetscErrorCode SetAnalyticalSolutionForParticles(Vec tempVec, SimCtx *simCtx)
{
    PetscErrorCode ierr;
    const char *analytical_type = simCtx->AnalyticalSolutionType[0] ? simCtx->AnalyticalSolutionType : "default";

    PetscFunctionBeginUser;
    
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Type: %s\n", analytical_type);
    
    // Check for specific analytical solution types
    if (strcmp(simCtx->AnalyticalSolutionType, "TGV3D") == 0) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Using TGV3D solution.\n");
        ierr = SetAnalyticalSolutionForParticles_TGV3D(tempVec, simCtx); CHKERRQ(ierr);
        return 0;
    }
    if (strcmp(simCtx->AnalyticalSolutionType, "UNIFORM_FLOW") == 0) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Using UNIFORM_FLOW solution.\n");
        ierr = SetAnalyticalSolutionForParticles_UniformFlow(tempVec, simCtx); CHKERRQ(ierr);
        return 0;
    }
    
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper that evaluates the configured scalar verification profile.
 * @details Local to this translation unit.
 */
static PetscErrorCode EvaluateConfiguredScalarProfile(const VerificationScalarConfig *cfg,
                                                      PetscReal x,
                                                      PetscReal y,
                                                      PetscReal z,
                                                      PetscReal *value)
{
    PetscFunctionBeginUser;
    if (!cfg) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "VerificationScalarConfig cannot be NULL.");
    if (!value) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Scalar output pointer cannot be NULL.");

    if (strcmp(cfg->profile, "CONSTANT") == 0) {
        *value = cfg->value;
    } else if (strcmp(cfg->profile, "LINEAR_X") == 0) {
        *value = cfg->phi0 + cfg->slope_x * x;
    } else if (strcmp(cfg->profile, "SIN_PRODUCT") == 0) {
        *value = cfg->amplitude *
                 PetscSinReal(cfg->kx * x) *
                 PetscSinReal(cfg->ky * y) *
                 PetscSinReal(cfg->kz * z);
    } else {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
                "Unsupported verification scalar profile '%s'.", cfg->profile);
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EvaluateAnalyticalScalarProfile"
/**
 * @brief Implementation of \ref EvaluateAnalyticalScalarProfile().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/AnalyticalSolutions.h`.
 * @see EvaluateAnalyticalScalarProfile()
 */
PetscErrorCode EvaluateAnalyticalScalarProfile(const SimCtx *simCtx,
                                               PetscReal x,
                                               PetscReal y,
                                               PetscReal z,
                                               PetscReal t,
                                               PetscReal *value)
{
    PetscFunctionBeginUser;
    (void)t;
    if (!simCtx) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "SimCtx cannot be NULL.");
    if (!simCtx->verificationScalar.enabled) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
                "EvaluateAnalyticalScalarProfile requires verification scalar mode to be enabled.");
    }
    PetscCall(EvaluateConfiguredScalarProfile(&simCtx->verificationScalar, x, y, z, value));
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetAnalyticalScalarFieldOnParticles"
/**
 * @brief Implementation of \ref SetAnalyticalScalarFieldOnParticles().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/AnalyticalSolutions.h`.
 * @see SetAnalyticalScalarFieldOnParticles()
 */
PetscErrorCode SetAnalyticalScalarFieldOnParticles(UserCtx *user, const char *swarm_field_name)
{
    PetscErrorCode ierr;
    PetscInt       nlocal = 0;
    PetscReal     *positions = NULL;
    PetscReal     *scalar_values = NULL;

    PetscFunctionBeginUser;
    if (!user) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx cannot be NULL.");
    if (!user->swarm) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "UserCtx->swarm is NULL.");
    if (!swarm_field_name) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "swarm_field_name cannot be NULL.");

    ierr = DMSwarmGetLocalSize(user->swarm, &nlocal); CHKERRQ(ierr);
    if (nlocal == 0) PetscFunctionReturn(0);

    ierr = DMSwarmGetField(user->swarm, "position", NULL, NULL, (void **)&positions); CHKERRQ(ierr);
    ierr = DMSwarmGetField(user->swarm, swarm_field_name, NULL, NULL, (void **)&scalar_values); CHKERRQ(ierr);

    for (PetscInt p = 0; p < nlocal; ++p) {
        PetscReal value = 0.0;
        ierr = EvaluateAnalyticalScalarProfile(user->simCtx,
                                               positions[3 * p + 0],
                                               positions[3 * p + 1],
                                               positions[3 * p + 2],
                                               user->simCtx->ti,
                                               &value); CHKERRQ(ierr);
        scalar_values[p] = value;
    }

    ierr = DMSwarmRestoreField(user->swarm, swarm_field_name, NULL, NULL, (void **)&scalar_values); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void **)&positions); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SetAnalyticalScalarFieldAtCellCenters"
/**
 * @brief Implementation of \ref SetAnalyticalScalarFieldAtCellCenters().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/AnalyticalSolutions.h`.
 * @see SetAnalyticalScalarFieldAtCellCenters()
 */
PetscErrorCode SetAnalyticalScalarFieldAtCellCenters(UserCtx *user, Vec targetVec)
{
    PetscErrorCode  ierr;
    PetscReal     ***target = NULL;
    const Cmpnts  ***cent = NULL;
    DMDALocalInfo   info;
    PetscInt        xs, xe, ys, ye, zs, ze, mx, my, mz;
    PetscInt        lxs, lxe, lys, lye, lzs, lze;

    PetscFunctionBeginUser;
    if (!user) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx cannot be NULL.");
    if (!targetVec) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "targetVec cannot be NULL.");

    info = user->info;
    xs = info.xs; xe = info.xs + info.xm;
    ys = info.ys; ye = info.ys + info.ym;
    zs = info.zs; ze = info.zs + info.zm;
    mx = info.mx; my = info.my; mz = info.mz;
    lxs = (xs == 0) ? xs + 1 : xs; lxe = (xe == mx) ? xe - 1 : xe;
    lys = (ys == 0) ? ys + 1 : ys; lye = (ye == my) ? ye - 1 : ye;
    lzs = (zs == 0) ? zs + 1 : zs; lze = (ze == mz) ? ze - 1 : ze;

    ierr = VecSet(targetVec, 0.0); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, targetVec, &target); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->Cent, &cent); CHKERRQ(ierr);

    for (PetscInt k = lzs; k < lze; ++k) {
        for (PetscInt j = lys; j < lye; ++j) {
            for (PetscInt i = lxs; i < lxe; ++i) {
                PetscReal value = 0.0;
                ierr = EvaluateAnalyticalScalarProfile(user->simCtx,
                                                       cent[k][j][i].x,
                                                       cent[k][j][i].y,
                                                       cent[k][j][i].z,
                                                       user->simCtx->ti,
                                                       &value); CHKERRQ(ierr);
                target[k][j][i] = value;
            }
        }
    }

    ierr = DMDAVecRestoreArrayRead(user->fda, user->Cent, &cent); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, targetVec, &target); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}
