#include "momentumsolvers.h"


/*================================================================================*
 *               SHARED PHYSICAL-TIME (BDF) COEFFICIENT PLUMBING                  *
 *================================================================================*/

#undef __FUNCT__
#define __FUNCT__ "MomentumUsesBDF2"
/**
 * @brief Single source of truth for the BDF1/BDF2 selection. See header.
 *
 * This reproduces exactly the predicate historically inlined in
 * ComputeTotalResidual(): BDF1 on the cold-start first step (ti == 1) and on the
 * first step after a restart (ti == tistart); BDF2 thereafter (when the second-
 * order time accuracy is enabled). Restart-history validity is a separate concern
 * and is intentionally not addressed here.
 */
PetscBool MomentumUsesBDF2(SimCtx *simCtx)
{
    const PetscInt ti      = simCtx->step;
    const PetscInt tistart = simCtx->StartStep;
    return (PetscBool)(COEF_TIME_ACCURACY > 1.1 && ti != tistart && ti != 1);
}

#undef __FUNCT__
#define __FUNCT__ "MomentumBDFCoefficient"
/**
 * @brief Returns a0 in {1.0, 1.5} for the current physical step. See header.
 */
PetscReal MomentumBDFCoefficient(SimCtx *simCtx)
{
    return MomentumUsesBDF2(simCtx) ? COEF_TIME_ACCURACY : 1.0;
}


#undef __FUNCT__
#define __FUNCT__ "ComputeTotalResidual"
/**
 * @brief Shared implementation of `ComputeTotalResidual()`.
 */
PetscErrorCode ComputeTotalResidual(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx *simCtx = user->simCtx;
    
    // Extract Time Parameters from Context.
    // BDF step selection now lives in MomentumUsesBDF2() (shared with the stability estimate).
    const PetscReal dt      = simCtx->dt;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    // 1. Calculate Spatial Terms (stored in user->Rhs)
    //    Rhs = -Div(Flux) + Viscous + Source
    ierr = ComputeRHS(user, user->Rhs); CHKERRQ(ierr);

    // 2. Add Physical Time Derivative Terms (BDF Discretization)
    //    The equation solved is: dU/dtau = RHS_Spatial + RHS_Temporal
    
    /* BDF order + leading coefficient from the shared helpers (single source of truth,
       shared with the momentum stability estimate). a0 = 1.5 (==COEF_TIME_ACCURACY) for
       BDF2, 1.0 for BDF1; using a0 for the current-state term keeps the residual
       numerically identical to the historical inlined coefficients. */
    const PetscBool use_bdf2 = MomentumUsesBDF2(simCtx);
    const PetscReal a0       = MomentumBDFCoefficient(simCtx);
    if (use_bdf2) {
        // --- BDF2 (Second Order Backward Difference) ---
        // (a0*U^{n} - 2.0*U^{n-1} + 0.5*U^{n-2}) / dt = RHS_Spatial(U^{n})
        ierr = VecAXPY(user->Rhs, -a0/dt,  user->Ucont);     CHKERRQ(ierr);
        ierr = VecAXPY(user->Rhs, +2.0/dt, user->Ucont_o);   CHKERRQ(ierr);
        ierr = VecAXPY(user->Rhs, -0.5/dt, user->Ucont_rm1); CHKERRQ(ierr);
    } else {
        // --- BDF1 (First Order / Euler Implicit) ---
        // (a0*U^{n} - U^{n-1}) / dt = RHS_Spatial(U^{n}), with a0 = 1.0
        ierr = VecAXPY(user->Rhs, -a0/dt,  user->Ucont);     CHKERRQ(ierr);
        ierr = VecAXPY(user->Rhs, +1.0/dt, user->Ucont_o);   CHKERRQ(ierr);
    }

    // 3. Enforce Boundary Conditions on the Residual
    ierr = EnforceRHSBoundaryConditions(user); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "Momentum_Solver_Explicit_RungeKutta4"
/**
 * @brief Internal helper implementation: `MomentumSolver_Explicit_RungeKutta4()`.
 * @details Local to this translation unit.
 */
PetscErrorCode MomentumSolver_Explicit_RungeKutta4(UserCtx *user, IBMNodes *ibm, FSInfo *fsi)
{
    PetscErrorCode ierr;
    (void)fsi;
    // --- Context Acquisition ---
    // Get the master simulation context from the first block's UserCtx.
    // This is the bridge to access all former global variables.
    SimCtx *simCtx = user[0].simCtx;
    PetscReal dt = simCtx->dt;
    //PetscReal st = simCtx->st;
    PetscInt  istage;
    PetscReal alfa[] = {0.25, 1.0/3.0, 0.5, 1.0};

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Executing explicit momentum solver (Runge-Kutta) for %d block(s).\n",simCtx->block_number);

    // --- 1. Pre-Loop Initialization (Legacy Logic) ---
    // This block prepares boundary conditions and allocates the RHS vector for all blocks
    // before the main RK loop begins. This logic is preserved from the original code.
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Preparing all blocks for Runge-Kutta solve...\n");
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {

        ierr = ApplyBoundaryConditions(&user[bi]); CHKERRQ(ierr);

        /*
        // Immersed boundary interpolation (if enabled)
        if (simCtx->immersed) {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "  Performing pre-RK IBM interpolation for block %d.\n", bi);
            for (PetscInt ibi = 0; ibi < simCtx->NumberOfBodies; ibi++) {
                // The 'ibm' and 'fsi' pointers are passed directly from FlowSolver
                ierr = ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi, 1); CHKERRQ(ierr);
            }
        }
        */
        
        // Allocate the persistent RHS vector for this block's context
        ierr = VecDuplicate(user[bi].Ucont, &user[bi].Rhs); CHKERRQ(ierr);
    }
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Pre-loop initialization complete.\n");

    // --- 2. Main Runge-Kutta Loop ---
    // The legacy code had an outer `pseudot` loop that only ran once. We preserve it.
    for (PetscInt pseudot = 0; pseudot < 1; pseudot++) {
        // Loop over each block to perform the RK stages
        for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
            for (istage = 0; istage < 4; istage++) {
                LOG_ALLOW(LOCAL, LOG_DEBUG, "  Block %d, RK Stage %d (alpha=%.4f)...\n", bi, istage, alfa[istage]);

                // a. Calculate the Right-Hand Side (RHS) of the momentum equation.
                ierr = ComputeRHS(&user[bi], user[bi].Rhs); CHKERRQ(ierr);

                // b. Advance Ucont to the next intermediate stage using the RK coefficient.
                //    Ucont_new = Ucont_old + alpha * dt * RHS
                ierr = VecWAXPY(user[bi].Ucont, alfa[istage] * dt, user[bi].Rhs, user[bi].Ucont_o); CHKERRQ(ierr);
 
                // c. Synchronize periodic endpoints and local ghosts for the new intermediate Ucont.
                const char *staggered_fields[] = {"Ucont"};
                ierr = SynchronizePeriodicStaggeredFields(&user[bi], 1, staggered_fields); CHKERRQ(ierr);

                // d. Re-apply boundary conditions for the new intermediate velocity.
                //    This is crucial for the stability and accuracy of the multi-stage scheme.
                ierr = ApplyBoundaryConditions(&user[bi]); CHKERRQ(ierr);
                
            } // End of RK stages for one block

            /*
            // Final IBM Interpolation for the block (if enabled)
            if (simCtx->immersed) {
                LOG_ALLOW(LOCAL, LOG_DEBUG, "  Performing post-RK IBM interpolation for block %d.\n", bi);
                for (PetscInt ibi = 0; ibi < simCtx->NumberOfBodies; ibi++) {
                    ierr = ibm_interpolation_advanced(&user[bi], &ibm[ibi], ibi, 1); CHKERRQ(ierr);
                }
            }
            */

        } // End loop over blocks

        // --- 3. Inter-Block Communication (Legacy Logic) ---
        // This is called after all blocks have completed their RK stages.
        if (simCtx->block_number > 1) {
	  //    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Updating multi-block interfaces after RK stages.\n");
	    //   ierr = Block_Interface_U(user); CHKERRQ(ierr);
        }

    } // End of pseudo-time loop

    // --- 4. Cleanup ---
    // Destroy the RHS vectors that were created at the start of this function.
    for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
        ierr = VecDestroy(&user[bi].Rhs); CHKERRQ(ierr);
    }
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Runge-Kutta solve completed for all blocks.\n");

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeGlobalSpectralRadiusEstimate"
/**
 * @brief Compute a conservative global pseudo-time spectral radius estimate.
 *
 * The spectral radius lambda_max (units: 1/s) governs the pseudo-time stability limit
 * for the explicit Jameson 4-stage RK smoother. The pseudo-time step is then:
 *
 *     dtau = pseudo_cfl / lambda_max
 *
 * making pseudo_cfl a true dimensionless Courant number independent of the physical
 * timestep dt. Note: 2.83 is only the 4-stage RK *imaginary-axis* scalar limit; it is NOT
 * a generally-stable bound for the actual nonlinear/non-normal operator (see A4 validation).
 *
 * Per owned cell (k,j,i) the convective spectral radius contribution is:
 *
 *     lambda_cell = (|U_xi_flux| + |U_eta_flux| + |U_zeta_flux|) x (1/V)
 *
 * where:
 *   - U_xi_flux   = ucont[k][j][i].x : volumetric flux through the xi-face   [m3/s]
 *   - U_eta_flux  = ucont[k][j][i].y : volumetric flux through the eta-face   [m3/s]
 *   - U_zeta_flux = ucont[k][j][i].z : volumetric flux through the zeta-face  [m3/s]
 *   - 1/V = lAj[k][j][i]             : inverse cell volume (cell Jacobian)    [1/m3]
 *
 * Product units: [m3/s * 1/m3] = [1/s] = spectral radius.
 *
 * Only the 'positive' face at index (k,j,i) is used per direction (no neighbor
 * access), so no ghost-cell scatter is needed. lUcont is valid after the pre-loop
 * SynchronizePeriodicStaggeredFields; lAj is set at grid initialisation and never
 * changes. No DMGlobalToLocal calls are added.
 *
 * A BDF2 lower bound COEF_TIME_ACCURACY/dt is applied after the MPI reduction so
 * that zero-flow startup (all ucont == 0) gives dtau = pseudo_cfl * dt / 1.5 rather
 * than infinity.
 *
 * @param[in]  user         Array of UserCtx (one per block).
 * @param[in]  block_number Number of blocks.
 * @param[in]  dt           Physical time step (for the BDF2 lower bound).
 * @param[out] lambda_max_out  Global spectral radius [1/s].
 */
static PetscErrorCode ComputeGlobalSpectralRadiusEstimate(
    UserCtx *user, PetscInt block_number, PetscReal dt, PetscReal *lambda_max_out)
{
    PetscErrorCode  ierr;
    PetscReal       local_max = 0.0;

    PetscFunctionBeginUser;

    for (PetscInt bi = 0; bi < block_number; bi++) {
        DMDALocalInfo info = user[bi].info; /* pre-stored local ownership info */
        Cmpnts      ***ucont;
        PetscReal   ***aj;

        /* Read-only access: lUcont (valid after SynchronizePeriodicStaggeredFields),
           lAj (grid metric, unchanged after initialisation). */
        ierr = DMDAVecGetArrayRead(user[bi].fda, user[bi].lUcont, &ucont); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(user[bi].da,  user[bi].lAj,    &aj);    CHKERRQ(ierr);

        /* Loop over owned cells only — no ghost-cell indices accessed. */
        for (PetscInt k = info.zs; k < info.zs + info.zm; k++) {
            for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
                for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
                    /* Sum absolute face-flux magnitudes (one face per coordinate direction).
                     * Using only the 'positive' face (index i/j/k, not i-1/j-1/k-1) avoids
                     * neighbor access. For smooth velocity fields the error is < 2x and the
                     * adaptive CFL controller handles any residual conservatism. */
                    PetscReal flux_sum = PetscAbsReal(ucont[k][j][i].x)
                                      + PetscAbsReal(ucont[k][j][i].y)
                                      + PetscAbsReal(ucont[k][j][i].z);

                    /* aj = 1/J = 1/cell_volume; lambda = flux/volume [1/s] */
                    PetscReal lambda = flux_sum * aj[k][j][i];
                    local_max = PetscMax(local_max, lambda);
                }
            }
        }

        ierr = DMDAVecRestoreArrayRead(user[bi].fda, user[bi].lUcont, &ucont); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(user[bi].da,  user[bi].lAj,    &aj);    CHKERRQ(ierr);
    }

    /* Reduce across all MPI ranks to get the global maximum. */
    ierr = MPI_Allreduce(&local_max, lambda_max_out, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD);
    CHKERRQ(ierr);

    /* BDF2 lower bound: if the velocity field is zero (startup, stagnant IC) lambda_max
       would be 0, giving dtau = inf. Clamp to 1.5/dt so dtau = pseudo_cfl * dt/1.5 in
       that degenerate case (conservative but finite). */
    *lambda_max_out = PetscMax(*lambda_max_out, COEF_TIME_ACCURACY / dt);

    PetscFunctionReturn(0);
}


/* Absolute metric-row sum for one face: |N.N| + |A.N| + |B.N|, where N is the
   face's own normal-metric vector and A,B are the other two face-metric vectors.
   Mirrors the g11/g21/g31 rows assembled per face in Viscous(). */
static inline PetscReal MomFaceGabs(Cmpnts N, Cmpnts A, Cmpnts B)
{
    const PetscReal NN = N.x*N.x + N.y*N.y + N.z*N.z;
    const PetscReal AN = A.x*N.x + A.y*N.y + A.z*N.z;
    const PetscReal BN = B.x*N.x + B.y*N.y + B.z*N.z;
    return PetscAbsReal(NN) + PetscAbsReal(AN) + PetscAbsReal(BN);
}

/* nvert bands matching the physical operator:
 *   SKIP    : cell excluded from the estimate (solid/blanked), nvert > 0.1.
 *   VISC1S  : Viscous() one-sided cross-derivative trigger band, (0.5, 7).
 *   QUICKBLK: Convection() QUICK-modifier band, [0.1, 7] (mirrors NOT(<0.1 || >innerblank)). */
#define MOM_SKIP_SOLID(v)    ((v) > 0.1)
#define MOM_VISC_ONESIDED(v) ((v) > 0.5 && (v) < 7.0)
#define MOM_QUICK_BLOCKS(v)  ((v) >= 0.1 && (v) <= 7.0)

/* Read-only: does ANY of the six viscous faces entering cell (k,j,i) use an
 * nvert-selected one-sided cross-derivative? Conservative solid-band check over
 * the full 3x3x3 neighborhood that the per-face Viscous() predicates can reach.
 * Indices i +/- 1, j +/- 1, k +/- 1 are always valid on the interior loop range. */
static PetscBool MomCellUsesOneSidedViscousStencil(PetscReal ***nvert, PetscInt k, PetscInt j, PetscInt i)
{
    for (PetscInt dk = -1; dk <= 1; dk++)
        for (PetscInt dj = -1; dj <= 1; dj++)
            for (PetscInt di = -1; di <= 1; di++) {
                if (!dk && !dj && !di) continue;
                if (MOM_VISC_ONESIDED(nvert[k+dk][j+dj][i+di])) return PETSC_TRUE;
            }
    return PETSC_FALSE;
}

/* Read-only: is cell (k,j,i) IB-adjacent (any solid in its 3x3x3 neighborhood)?
 * Uses the SKIP band so it is consistent with the one-sided viscous trigger. */
static PetscBool MomCellHasSolidNeighbor(PetscReal ***nvert, PetscInt k, PetscInt j, PetscInt i)
{
    for (PetscInt dk = -1; dk <= 1; dk++)
        for (PetscInt dj = -1; dj <= 1; dj++)
            for (PetscInt di = -1; di <= 1; di++) {
                if (!dk && !dj && !di) continue;
                if (MOM_SKIP_SOLID(nvert[k+dk][j+dj][i+di])) return PETSC_TRUE;
            }
    return PETSC_FALSE;
}

/* Read-only: is the QUICK reconstruction in one computational direction modified for
 * cell (k,j,i)? The cell uses BOTH faces in that direction (index f and f-1), whose
 * QUICK branches are gated by a non-periodic domain edge or by nvert blockers within
 * the stencil reach (including the second neighbor f-2). The +/-2 reads are guarded by
 * the ghost range [g0, g1) so they are safe for both periodic and non-periodic DMs. */
static PetscBool MomQuickDirModified(PetscReal ***nvert, PetscInt a, PetscInt b, PetscInt c,
                                     PetscInt idx, PetscInt m, PetscBool np0, PetscBool np1,
                                     PetscInt g0, PetscInt g1, char dir)
{
    /* (a,b,c) are the fixed indices; idx is the moving index in direction `dir`.
       The two contributing faces (idx, idx-1) have QUICK stencils spanning idx-2..idx+2. */
    if (np0 && idx <= 1)   return PETSC_TRUE;     /* faces idx,idx-1 reach the negative edge */
    if (np1 && idx >= m-2) return PETSC_TRUE;     /* ... or the positive edge */
    for (PetscInt d = -2; d <= 2; d++) {
        if (d == 0) continue;
        const PetscInt p = idx + d;
        /* Required stencil information unavailable in the local ghost range: do NOT assume
           fluid -> conservatively classify the direction as modified (use 2.5). */
        if (p < g0 || p >= g1) return PETSC_TRUE;
        PetscReal v;
        if      (dir == 'i') v = nvert[a][b][p];
        else if (dir == 'j') v = nvert[a][p][c];
        else                 v = nvert[p][b][c];
        if (MOM_QUICK_BLOCKS(v)) return PETSC_TRUE;
    }
    return PETSC_FALSE;
}

/* Active-row mask (3 bits: xi=1, eta=2, zeta=4) mirroring the final RHS masking.
 * A solid cell disables all rows; a positive solid neighbour or a positive non-periodic
 * physical face disables the corresponding normal (staggered) row; TwoD disables the
 * homogeneous direction. Returns 0 only when the location carries no active unknown. */
PetscInt MomCellActiveRows(PetscReal ***nvert, PetscInt k, PetscInt j, PetscInt i,
                           PetscInt mx, PetscInt my, PetscInt mz,
                           PetscBool np_x1, PetscBool np_y1, PetscBool np_z1, PetscInt twoD)
{
    if (MOM_SKIP_SOLID(nvert[k][j][i])) return 0;     /* solid cell: all rows inactive */
    PetscInt bits = 0x7;
    if (MOM_SKIP_SOLID(nvert[k][j][i+1])) bits &= ~0x1;  /* positive xi neighbour solid   */
    if (MOM_SKIP_SOLID(nvert[k][j+1][i])) bits &= ~0x2;  /* positive eta neighbour solid  */
    if (MOM_SKIP_SOLID(nvert[k+1][j][i])) bits &= ~0x4;  /* positive zeta neighbour solid */
    if (np_x1 && i == mx-2) bits &= ~0x1;                /* positive non-periodic xi face  */
    if (np_y1 && j == my-2) bits &= ~0x2;
    if (np_z1 && k == mz-2) bits &= ~0x4;
    if      (twoD == 1) bits &= ~0x1;                    /* TwoD homogeneous direction: xi  */
    else if (twoD == 2) bits &= ~0x2;                    /*                             eta  */
    else if (twoD == 3) bits &= ~0x4;                    /*                             zeta */
    return bits;
}

#undef __FUNCT__
#define __FUNCT__ "ComputeMomentumStabilityEstimate"
/**
 * @brief Practical conservative momentum pseudo-time stability estimate (shadow). See header.
 *
 * Operator-scaled estimate lambda = max_cell (a0/dt + lambda_c + lambda_nu) over active,
 * non-solid interior cells, blocks and MPI ranks (lambda_c already carries the per-direction
 * QUICK scheme factors -- it is NOT multiplied by f_c again here). Read-only; no halo exchange,
 * 5 global scalar collectives. NOT a proven spectral radius. Drives nothing in shadow mode.
 */
PetscErrorCode ComputeMomentumStabilityEstimate(UserCtx *user, PetscInt block_number,
                                                PetscReal dt, MomStabCandidate candidate,
                                                MomStabilityReport *rep)
{
    PetscErrorCode ierr;
    PetscMPIInt rank;
    PetscFunctionBeginUser;

    /* ---- input validation (no silent fallback) ---- */
    PetscCheck(user != NULL, PETSC_COMM_WORLD, PETSC_ERR_ARG_NULL, "user is NULL");
    PetscCheck(rep  != NULL, PETSC_COMM_WORLD, PETSC_ERR_ARG_NULL, "rep is NULL");
    PetscCheck(block_number > 0, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "block_number must be > 0");
    PetscCheck(candidate >= MOM_STAB_CAND_B && candidate <= MOM_STAB_CAND_D,
               PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "unsupported stability candidate enum");
    PetscCheck(PetscIsNormalReal(dt) && dt > 0.0, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
               "dt must be finite and positive");

    SimCtx         *simCtx   = user[0].simCtx;
    PetscCheck(simCtx != NULL, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "user[0].simCtx is NULL");
    const PetscReal a0       = MomentumBDFCoefficient(simCtx);
    const PetscBool centered = (PetscBool)(simCtx->les || simCtx->central);
    const PetscBool has_nut  = (PetscBool)(simCtx->les || simCtx->rans);
    const PetscBool inviscid = (PetscBool)simCtx->invicid;
    PetscCheck(inviscid || (PetscIsNormalReal(simCtx->ren) && simCtx->ren > 0.0),
               PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE, "Reynolds number must be finite and positive when viscous");
    const PetscReal lambda_t = a0 / dt;
    const PetscReal nu_mol   = inviscid ? 0.0 : 1.0 / simCtx->ren;
    const PetscInt  twoD     = simCtx->TwoD;

    /* Shadow-mode completeness flag. The estimate does NOT cover the Clark nonlinear stress
       Jacobian, nor RANS eddy-viscosity sign behaviour (not verified sign-definite).
       Body forces in the supported configs read a per-timestep-frozen scalar
       (simCtx->bulkVelocityCorrection) inside ComputeRHS, so within a pseudo-solve they are a
       constant forcing with ZERO velocity Jacobian (consistent with the frozen-pressure
       treatment) -- they do not make the estimate incomplete. */
    rep->estimate_incomplete = (PetscBool)(simCtx->clark || simCtx->rans);

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    /* Per-candidate local maxima; selected candidate also tracks the controlling cell. */
    PetscReal locmaxB = 0.0, locmaxC = 0.0, locmaxD = 0.0, sel_max = 0.0;
    PetscReal sel_lc = 0.0, sel_lv = 0.0;
    PetscInt  sel_ci = -1, sel_cj = -1, sel_ck = -1, sel_blk = -1, sel_class = 0, sel_os = 0;
    PetscInt  local_active = 0;

    for (PetscInt bi = 0; bi < block_number; bi++) {
        UserCtx *u = &user[bi];
        DMDALocalInfo info = u->info;
        const PetscInt mx = info.mx, my = info.my, mz = info.mz;
        const PetscInt xs = info.xs, xe = xs + info.xm;
        const PetscInt ys = info.ys, ye = ys + info.ym;
        const PetscInt zs = info.zs, ze = zs + info.zm;
        /* Interior range == where ComputeRHS evaluates the residual; boundary faces
           (0 and m-1) are RHS-masked, so excluding them implements active-row masking
           for the normal boundary rows. The last physical face (m-2) stays IN with the
           full cell estimate (its tangential rows remain active). */
        const PetscInt lxs = (xs==0)?xs+1:xs, lxe = (xe==mx)?xe-1:xe;
        const PetscInt lys = (ys==0)?ys+1:ys, lye = (ye==my)?ye-1:ye;
        const PetscInt lzs = (zs==0)?zs+1:zs, lze = (ze==mz)?ze-1:ze;

        /* Non-periodic boundary flags (QUICK stencil is degraded near these). */
        const PetscBool npx0 = (PetscBool)(u->boundary_faces[BC_FACE_NEG_X].mathematical_type != PERIODIC);
        const PetscBool npx1 = (PetscBool)(u->boundary_faces[BC_FACE_POS_X].mathematical_type != PERIODIC);
        const PetscBool npy0 = (PetscBool)(u->boundary_faces[BC_FACE_NEG_Y].mathematical_type != PERIODIC);
        const PetscBool npy1 = (PetscBool)(u->boundary_faces[BC_FACE_POS_Y].mathematical_type != PERIODIC);
        const PetscBool npz0 = (PetscBool)(u->boundary_faces[BC_FACE_NEG_Z].mathematical_type != PERIODIC);
        const PetscBool npz1 = (PetscBool)(u->boundary_faces[BC_FACE_POS_Z].mathematical_type != PERIODIC);
        const PetscBool wxn = (PetscBool)(u->boundary_faces[BC_FACE_NEG_X].mathematical_type == WALL);
        const PetscBool wxp = (PetscBool)(u->boundary_faces[BC_FACE_POS_X].mathematical_type == WALL);
        const PetscBool wyn = (PetscBool)(u->boundary_faces[BC_FACE_NEG_Y].mathematical_type == WALL);
        const PetscBool wyp = (PetscBool)(u->boundary_faces[BC_FACE_POS_Y].mathematical_type == WALL);
        const PetscBool wzn = (PetscBool)(u->boundary_faces[BC_FACE_NEG_Z].mathematical_type == WALL);
        const PetscBool wzp = (PetscBool)(u->boundary_faces[BC_FACE_POS_Z].mathematical_type == WALL);

        Cmpnts ***ucont, ***ucat = NULL;
        Cmpnts ***icsi, ***ieta, ***izet, ***jcsi, ***jeta, ***jzet, ***kcsi, ***keta, ***kzet;
        Cmpnts ***csi = NULL, ***eta = NULL, ***zet = NULL;
        PetscReal ***aj, ***iaj, ***jaj, ***kaj, ***nvert, ***nut = NULL;

        ierr = DMDAVecGetArrayRead(u->fda, u->lUcont, &ucont); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(u->da,  u->lAj,   &aj);     CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(u->da,  u->lIAj,  &iaj);    CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(u->da,  u->lJAj,  &jaj);    CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(u->da,  u->lKAj,  &kaj);    CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(u->fda, u->lICsi, &icsi);   CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(u->fda, u->lIEta, &ieta);   CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(u->fda, u->lIZet, &izet);   CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(u->fda, u->lJCsi, &jcsi);   CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(u->fda, u->lJEta, &jeta);   CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(u->fda, u->lJZet, &jzet);   CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(u->fda, u->lKCsi, &kcsi);   CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(u->fda, u->lKEta, &keta);   CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(u->fda, u->lKZet, &kzet);   CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(u->da,  u->lNvert, &nvert); CHKERRQ(ierr);
        if (has_nut) { ierr = DMDAVecGetArrayRead(u->da, u->lNu_t, &nut); CHKERRQ(ierr); }
        /* Cartesian velocity + cell metrics: always read so candidate D's gradient term
           (and hence rep->lambda_D) is meaningful regardless of the selected candidate.
           Precondition: lUcat fresh (Contra2Cart has run in the preceding residual eval). */
        ierr = DMDAVecGetArrayRead(u->fda, u->lUcat, &ucat); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(u->fda, u->lCsi,  &csi);  CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(u->fda, u->lEta,  &eta);  CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(u->fda, u->lZet,  &zet);  CHKERRQ(ierr);

        /* Error flag for the cell loop: on a bad metric/estimate we record the location and
           jump to block_cleanup so EVERY acquired array is restored before the error returns. */
        PetscErrorCode cell_err = 0;
        PetscInt       be_i = -1, be_j = -1, be_k = -1;
        const char    *be_what = NULL;

        for (PetscInt k = lzs; k < lze; k++) {
            for (PetscInt j = lys; j < lye; j++) {
                for (PetscInt i = lxs; i < lxe; i++) {
                    /* ---- active-row mask: skip only locations with no active unknown ---- */
                    const PetscInt rows = MomCellActiveRows(nvert, k, j, i, mx, my, mz,
                                                            npx1, npy1, npz1, twoD);
                    if (rows == 0) continue;
                    local_active++;
                    const PetscReal Ajc = aj[k][j][i];
                    if (!(PetscIsNormalReal(Ajc) && Ajc > 0.0)) {
                        cell_err = PETSC_ERR_FP; be_i = i; be_j = j; be_k = k;
                        be_what = "non-finite/non-positive cell inverse-Jacobian";
                        goto block_cleanup;
                    }

                    /* ---- directional QUICK modification (per direction; both faces; 2nd neighbour) ---- */
                    const PetscBool mod_x = MomQuickDirModified(nvert, k, j, i, i, mx, npx0, npx1,
                                                               info.gxs, info.gxs+info.gxm, 'i');
                    const PetscBool mod_y = MomQuickDirModified(nvert, k, j, i, j, my, npy0, npy1,
                                                               info.gys, info.gys+info.gym, 'j');
                    const PetscBool mod_z = MomQuickDirModified(nvert, k, j, i, k, mz, npz0, npz1,
                                                               info.gzs, info.gzs+info.gzm, 'k');

                    /* ---- classify cell: interior / physical-boundary / IB-adjacent ---- */
                    const PetscBool bnd = (PetscBool)((npx0 && i<=1) || (npx1 && i>=mx-2) ||
                                                      (npy0 && j<=1) || (npy1 && j>=my-2) ||
                                                      (npz0 && k<=1) || (npz1 && k>=mz-2));
                    const PetscBool ib  = MomCellHasSolidNeighbor(nvert, k, j, i);
                    const PetscInt cls = ib ? 2 : (bnd ? 1 : 0);

                    /* ---- convective: six contravariant face fluxes, DIRECTIONAL QUICK factors ---- */
                    const PetscReal Uxp = ucont[k][j][i].x,   Uxm = ucont[k][j][i-1].x;
                    const PetscReal Uyp = ucont[k][j][i].y,   Uym = ucont[k][j-1][i].y;
                    const PetscReal Uzp = ucont[k][j][i].z,   Uzm = ucont[k-1][j][i].z;
                    const PetscReal divU = PetscAbsReal((Uxp-Uxm)+(Uyp-Uym)+(Uzp-Uzm));
                    const PetscReal fx = centered ? 1.0 : (mod_x ? 2.5 : (4.0/3.0));
                    const PetscReal fy = centered ? 1.0 : (mod_y ? 2.5 : (4.0/3.0));
                    const PetscReal fz = centered ? 1.0 : (mod_z ? 2.5 : (4.0/3.0));

                    /* per-direction transport scale; C's divergence term stays unscaled by f_d. */
                    const PetscReal lcB = 0.5 * Ajc * (
                          fx * (PetscAbsReal(Uxp)+PetscAbsReal(Uxm))
                        + fy * (PetscAbsReal(Uyp)+PetscAbsReal(Uym))
                        + fz * (PetscAbsReal(Uzp)+PetscAbsReal(Uzm)) );
                    const PetscReal lcC = lcB + 0.5 * Ajc * divU;
                    /* lambda_grad_u = max_i sum_j |du_i/dx_j| from fresh Cartesian velocity
                       (the nonlinear zero-order Jacobian term delta_u . grad u).
                       Physical gradient: d/dx_dir = Ajc*(Csi.dir d/dxi + Eta.dir d/deta + Zet.dir d/dzeta). */
                    const Cmpnts duc = { 0.5*(ucat[k][j][i+1].x-ucat[k][j][i-1].x),
                                         0.5*(ucat[k][j][i+1].y-ucat[k][j][i-1].y),
                                         0.5*(ucat[k][j][i+1].z-ucat[k][j][i-1].z) };
                    const Cmpnts due = { 0.5*(ucat[k][j+1][i].x-ucat[k][j-1][i].x),
                                         0.5*(ucat[k][j+1][i].y-ucat[k][j-1][i].y),
                                         0.5*(ucat[k][j+1][i].z-ucat[k][j-1][i].z) };
                    const Cmpnts duz = { 0.5*(ucat[k+1][j][i].x-ucat[k-1][j][i].x),
                                         0.5*(ucat[k+1][j][i].y-ucat[k-1][j][i].y),
                                         0.5*(ucat[k+1][j][i].z-ucat[k-1][j][i].z) };
                    const Cmpnts C = csi[k][j][i], E = eta[k][j][i], Z = zet[k][j][i];
                    #define MOM_ROW(cmp) ( \
                        PetscAbsReal(Ajc*(C.x*duc.cmp + E.x*due.cmp + Z.x*duz.cmp)) + \
                        PetscAbsReal(Ajc*(C.y*duc.cmp + E.y*due.cmp + Z.y*duz.cmp)) + \
                        PetscAbsReal(Ajc*(C.z*duc.cmp + E.z*due.cmp + Z.z*duz.cmp)) )
                    const PetscReal lgrad = PetscMax(MOM_ROW(x), PetscMax(MOM_ROW(y), MOM_ROW(z)));
                    #undef MOM_ROW
                    const PetscReal lcD = lcC + lgrad;

                    /* ---- viscous: six faces, face Jacobians, full metric rows ---- */
                    PetscReal lv = 0.0;
                    if (!inviscid) {
                        /* xi+ (i), xi- (i-1) */
                        for (PetscInt s = 0; s < 2; s++) {
                            const PetscInt fi = (s==0) ? i : i-1;
                            if (!(PetscIsNormalReal(iaj[k][j][fi]) && iaj[k][j][fi] > 0.0)) {
                                cell_err = PETSC_ERR_FP; be_i=i; be_j=j; be_k=k;
                                be_what = "non-finite/non-positive xi-face inverse-Jacobian"; goto block_cleanup; }
                            PetscReal nuf = nu_mol;
                            if (has_nut) {
                                PetscReal nt = 0.5*(nut[k][j][fi] + nut[k][j][fi+1]);
                                if ((wxn && fi==0) || (wxp && fi==mx-2)) nt = 0.0;
                                nuf += nt;
                            }
                            lv += PetscAbsReal(nuf) * Ajc * iaj[k][j][fi]
                                * MomFaceGabs(icsi[k][j][fi], ieta[k][j][fi], izet[k][j][fi]);
                        }
                        /* eta+ (j), eta- (j-1) */
                        for (PetscInt s = 0; s < 2; s++) {
                            const PetscInt fj = (s==0) ? j : j-1;
                            if (!(PetscIsNormalReal(jaj[k][fj][i]) && jaj[k][fj][i] > 0.0)) {
                                cell_err = PETSC_ERR_FP; be_i=i; be_j=j; be_k=k;
                                be_what = "non-finite/non-positive eta-face inverse-Jacobian"; goto block_cleanup; }
                            PetscReal nuf = nu_mol;
                            if (has_nut) {
                                PetscReal nt = 0.5*(nut[k][fj][i] + nut[k][fj+1][i]);
                                if ((wyn && fj==0) || (wyp && fj==my-2)) nt = 0.0;
                                nuf += nt;
                            }
                            /* eta-face normal is the eta metric (jeta) -> pass as N. */
                            lv += PetscAbsReal(nuf) * Ajc * jaj[k][fj][i]
                                * MomFaceGabs(jeta[k][fj][i], jcsi[k][fj][i], jzet[k][fj][i]);
                        }
                        /* zeta+ (k), zeta- (k-1) */
                        for (PetscInt s = 0; s < 2; s++) {
                            const PetscInt fk = (s==0) ? k : k-1;
                            if (!(PetscIsNormalReal(kaj[fk][j][i]) && kaj[fk][j][i] > 0.0)) {
                                cell_err = PETSC_ERR_FP; be_i=i; be_j=j; be_k=k;
                                be_what = "non-finite/non-positive zeta-face inverse-Jacobian"; goto block_cleanup; }
                            PetscReal nuf = nu_mol;
                            if (has_nut) {
                                PetscReal nt = 0.5*(nut[fk][j][i] + nut[fk+1][j][i]);
                                if ((wzn && fk==0) || (wzp && fk==mz-2)) nt = 0.0;
                                nuf += nt;
                            }
                            /* zeta-face normal is the zeta metric (kzet) -> pass as N. */
                            lv += PetscAbsReal(nuf) * Ajc * kaj[fk][j][i]
                                * MomFaceGabs(kzet[fk][j][i], kcsi[fk][j][i], keta[fk][j][i]);
                        }
                        lv *= 4.0;   /* 2 (two-face row-sum) x 2 (longitudinal full-stress) */
                    }

                    /* one-sided viscous cross-derivative branch near a solid: conservative x2,
                       applied ONCE per cell, after a complete 3x3x3 solid-band check. */
                    PetscInt one_sided = 0;
                    if (!inviscid && MomCellUsesOneSidedViscousStencil(nvert, k, j, i)) {
                        one_sided = 1;
                        lv *= 2.0;
                    }

                    /* ---- candidate totals and running maxima ---- */
                    const PetscReal lB = lambda_t + lcB + lv;
                    const PetscReal lC = lambda_t + lcC + lv;
                    const PetscReal lD = lambda_t + lcD + lv;
                    if (PetscIsInfOrNanReal(lD)) {
                        cell_err = PETSC_ERR_FP; be_i=i; be_j=j; be_k=k;
                        be_what = "non-finite stability estimate"; goto block_cleanup;
                    }
                    locmaxB = PetscMax(locmaxB, lB);
                    locmaxC = PetscMax(locmaxC, lC);
                    locmaxD = PetscMax(locmaxD, lD);

                    const PetscReal lsel  = (candidate==MOM_STAB_CAND_B)?lB:(candidate==MOM_STAB_CAND_C)?lC:lD;
                    const PetscReal lcsel = (candidate==MOM_STAB_CAND_B)?lcB:(candidate==MOM_STAB_CAND_C)?lcC:lcD;
                    if (lsel > sel_max) {
                        sel_max = lsel; sel_lc = lcsel; sel_lv = lv;
                        sel_ci = i; sel_cj = j; sel_ck = k; sel_blk = bi;
                        sel_class = cls; sel_os = one_sided;
                    }
                }
            }
        }

        /* Single cleanup point: reached on normal completion AND on any cell-loop error.
           Every successful GetArrayRead above has exactly one matching restore here. */
        block_cleanup:
        ierr = DMDAVecRestoreArrayRead(u->fda, u->lUcont, &ucont); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(u->da,  u->lAj,   &aj);     CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(u->da,  u->lIAj,  &iaj);    CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(u->da,  u->lJAj,  &jaj);    CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(u->da,  u->lKAj,  &kaj);    CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(u->fda, u->lICsi, &icsi);   CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(u->fda, u->lIEta, &ieta);   CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(u->fda, u->lIZet, &izet);   CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(u->fda, u->lJCsi, &jcsi);   CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(u->fda, u->lJEta, &jeta);   CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(u->fda, u->lJZet, &jzet);   CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(u->fda, u->lKCsi, &kcsi);   CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(u->fda, u->lKEta, &keta);   CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(u->fda, u->lKZet, &kzet);   CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(u->da,  u->lNvert, &nvert); CHKERRQ(ierr);
        if (has_nut) { ierr = DMDAVecRestoreArrayRead(u->da, u->lNu_t, &nut); CHKERRQ(ierr); }
        ierr = DMDAVecRestoreArrayRead(u->fda, u->lUcat, &ucat); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(u->fda, u->lCsi,  &csi);  CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(u->fda, u->lEta,  &eta);  CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(u->fda, u->lZet,  &zet);  CHKERRQ(ierr);

        /* Now that this block's arrays are restored, surface any cell-loop error. */
        PetscCheck(cell_err == 0, PETSC_COMM_SELF, cell_err, "%s at (%" PetscInt_FMT
                   ",%" PetscInt_FMT ",%" PetscInt_FMT ")", be_what ? be_what : "error", be_i, be_j, be_k);
    }

    /* ---- global reductions (no ghost/halo exchange; scalar collectives only) ----
       5 collectives total: one 3-element MAX (B/C/D), one SUM (active-cell count),
       one MIN (portable owner selection), and two broadcasts (the controlling-cell
       real breakdown and integer indices). */
    PetscReal loc3[3] = { locmaxB, locmaxC, locmaxD }, glo3[3];
    ierr = MPI_Allreduce(loc3, glo3, 3, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);

    PetscInt global_active = 0;
    ierr = MPI_Allreduce(&local_active, &global_active, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);

    /* The selected candidate's global max equals glo3[candidate] (no extra reduction).
       Portable owner pick: lowest rank whose local sel_max equals the global max. */
    const PetscReal sel_global = glo3[(int)candidate];
    PetscMPIInt nproc, claim, owner;
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &nproc); CHKERRQ(ierr);
    claim = (sel_max == sel_global) ? rank : nproc;
    ierr = MPI_Allreduce(&claim, &owner, 1, MPI_INT, MPI_MIN, PETSC_COMM_WORLD); CHKERRQ(ierr);
    if (owner == nproc) owner = 0;   /* no active cell anywhere: deterministic fallback owner */

    PetscReal dbuf[2] = { sel_lc, sel_lv };
    PetscInt  ibuf[6] = { sel_ci, sel_cj, sel_ck, sel_blk, sel_class, sel_os };
    ierr = MPI_Bcast(dbuf, 2, MPIU_REAL, owner, PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Bcast(ibuf, 6, MPIU_INT,  owner, PETSC_COMM_WORLD); CHKERRQ(ierr);

    rep->lambda_t = lambda_t;
    rep->lambda_c = dbuf[0];
    rep->lambda_v = dbuf[1];
    rep->lambda   = sel_global;
    rep->lambda_B = glo3[0]; rep->lambda_C = glo3[1]; rep->lambda_D = glo3[2];
    rep->ci = ibuf[0]; rep->cj = ibuf[1]; rep->ck = ibuf[2];
    rep->cblock = ibuf[3]; rep->cclass = ibuf[4]; rep->one_sided = ibuf[5];
    rep->active_cells = global_active;

    /* Only a genuinely empty active set (all cells masked) may fall back to the temporal
       term alone. A non-finite or non-positive estimate with active cells is a hard error. */
    if (global_active == 0) {
        rep->lambda = lambda_t;
    } else {
        PetscCheck(PetscIsNormalReal(rep->lambda) && rep->lambda > 0.0, PETSC_COMM_WORLD, PETSC_ERR_FP,
                   "momentum stability estimate is non-finite or non-positive with %" PetscInt_FMT
                   " active cells", global_active);
    }

    /* Dominant limiter at the controlling cell. */
    if (rep->lambda_t >= rep->lambda_c && rep->lambda_t >= rep->lambda_v) rep->limiter = MOM_STAB_LIMITER_TIME;
    else if (rep->lambda_c >= rep->lambda_v)                              rep->limiter = MOM_STAB_LIMITER_CONVECTION;
    else                                                                  rep->limiter = MOM_STAB_LIMITER_VISCOSITY;

    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "MomentumSolver_DualTime_Picard_JamesonRK"
/**
 * @brief Internal helper implementation: `MomentumSolver_DualTime_Picard_JamesonRK()`.
 * @details Local to this translation unit.
 */
PetscErrorCode MomentumSolver_DualTime_Picard_JamesonRK(UserCtx *user, IBMNodes *ibm, FSInfo *fsi)
{
    (void)ibm;
    (void)fsi;
    // --- CONTEXT ACQUISITION BLOCK ---
    SimCtx *simCtx = user[0].simCtx;
    const PetscInt block_number = simCtx->block_number;
    
    // Legacy Names (Physics/Time) - Kept for brevity in formulas
    const PetscInt ti = simCtx->step;
    const PetscReal dt = simCtx->dt;       // Physical Time Step
    /* cfl: initial pseudo-CFL Courant number (dimensionless, flow-independent since Phase 3).
     * Actual dtau = cfl / lambda_max, computed after the spectral radius estimate below. */
    const PetscReal cfl = simCtx->pseudo_cfl;
    const PetscReal alfa[] = {0.25, 1.0/3.0, 0.5, 1.0}; // Jameson RK smoothing coefficients
    Vec *pRhs; // Per-block backup of Rhs at last accepted pseudo-state.

    // State flags
    PetscBool force_restart = PETSC_FALSE;

    // Renamed Solver Parameters
    const PetscInt  max_pseudo_steps = simCtx->mom_max_pseudo_steps; // Max Dual-Time Iterations
    const PetscReal tol_abs_delta    = simCtx->mom_atol; // Stop if |dU| < tol
    const PetscReal tol_rtol_delta   = simCtx->mom_rtol; // Stop if |dU|/|dU0| < tol
    // --- END CONTEXT ACQUISITION BLOCK ---

    PetscErrorCode ierr;
    PetscMPIInt    rank;
    PetscInt       istage, pseudo_iter;
    PetscInt       accepted_iter = 0, rejected_iter = 0, recovery_streak = 0;
    PetscReal      ts, te, cput;

    // --- Global Convergence Metrics ---
    PetscReal global_norm_delta = 10.0;
    PetscReal global_rel_delta  = 1.0;
    PetscReal global_norm_resid = 1.0;
    PetscReal global_rel_resid  = 1.0;
    PetscReal smoothed_trial_ratio = 1.0; /* EMA-smoothed ratio; initialized neutral */

    // --- Spectral-radius-based pseudo-time step ---
    /* lambda_max: global maximum spectral radius [1/s], computed once per physical timestep.
     * pseudo_dtau[bi]: physical pseudo-time step dtau = pseudo_cfl / lambda_max [time].
     *   Unlike the old scheme (dtau = pseudo_cfl * dt), this is independent of dt and makes
     *   pseudo_cfl a true Courant number. (2.83 is only the scalar imaginary-axis RK limit,
     *   not a generally-stable value for the actual operator -- to be characterized in A4.)
     * dtau_min / dtau_max: per-timestep bounds derived from simCtx->min/max_pseudo_cfl. */
    PetscReal  lambda_max = 0.0;
    PetscReal  dtau_min, dtau_max;

    // --- Local Metric Arrays (Per Block) ---
    PetscReal *delta_sol_norm_init, *delta_sol_norm_prev, *delta_sol_norm_curr, *delta_sol_rel_curr;
    PetscReal *resid_norm_init,     *resid_norm_prev,     *resid_norm_curr,     *resid_rel_curr;
    PetscReal *pseudo_dtau;          /* per-block pseudo-time step dtau [physical time, NOT dtau/dt] */
    PetscReal *trial_ratio_log;        /* per-block step-to-step residual ratio (for diagnostics) */
    PetscReal  last_accepted_resid;    /* last accepted |R|, for final summary */

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Executing Dual-Time Momentum Solver for %d block(s).\n", block_number);

    // --- Allocate metric arrays ---
    ierr = PetscMalloc2(block_number, &delta_sol_norm_init, block_number, &delta_sol_norm_prev); CHKERRQ(ierr);
    ierr = PetscMalloc2(block_number, &delta_sol_norm_curr, block_number, &delta_sol_rel_curr); CHKERRQ(ierr);
    ierr = PetscMalloc2(block_number, &resid_norm_init,     block_number, &resid_norm_prev); CHKERRQ(ierr);
    ierr = PetscMalloc2(block_number, &resid_norm_curr,     block_number, &pseudo_dtau); CHKERRQ(ierr);
    ierr = PetscMalloc1(block_number, &resid_rel_curr); CHKERRQ(ierr);
    ierr = PetscMalloc1(block_number, &trial_ratio_log); CHKERRQ(ierr);
    ierr = PetscMalloc1(block_number, &pRhs); CHKERRQ(ierr);
    last_accepted_resid = 0.0;

    ierr = PetscTime(&ts); CHKERRQ(ierr);

    // --- 1. Pre-Loop Initialization ---
    
    if (block_number > 1) {
      // ierr = Block_Interface_U(user); CHKERRQ(ierr);
    }

    for (PetscInt bi = 0; bi < block_number; bi++) {
        ierr = ApplyBoundaryConditions(&user[bi]); CHKERRQ(ierr);
        
        // Immersed boundary interpolation (if enabled)

        // Allocate workspace
        ierr = VecDuplicate(user[bi].Ucont, &user[bi].Rhs); CHKERRQ(ierr);
        ierr = VecDuplicate(user[bi].Rhs, &pRhs[bi]); CHKERRQ(ierr);
        ierr = VecDuplicate(user[bi].Ucont, &user[bi].dUcont); CHKERRQ(ierr);
        ierr = VecDuplicate(user[bi].Ucont, &user[bi].pUcont); CHKERRQ(ierr);
        
        // Initialize Backup (pUcont) with current state
        ierr = VecCopy(user[bi].Ucont, user[bi].pUcont); CHKERRQ(ierr);

        // --- Calculate Initial Total Residual (Spatial + Temporal) ---
        ierr = ComputeTotalResidual(&user[bi]); CHKERRQ(ierr);

        // Backup Rhs with current state
        ierr = VecCopy(user[bi].Rhs, pRhs[bi]); CHKERRQ(ierr);

        // Compute Initial Norms
        ierr = VecNorm(user[bi].Rhs, NORM_INFINITY, &resid_norm_init[bi]); CHKERRQ(ierr);
        
        // Initialize history for backtracking logic
        resid_norm_prev[bi]     = resid_norm_init[bi]; /* meaningful ratio on first iteration */
        delta_sol_norm_prev[bi] = 1000.0;
        /* pseudo_dtau[bi] is set after ComputeGlobalSpectralRadiusEstimate below. */

	    LOG_ALLOW(GLOBAL,LOG_INFO," Block %d | Max RHS = %.6f | initial pseudo-CFL = %.4f .\n", bi, resid_norm_init[bi], cfl);
    }

    /* --- Spectral Radius Estimate (computed once per physical timestep) ---
     * lambda_max [1/s] is the global maximum of (face_flux_sum * inverse_cell_volume)
     * over all owned cells across all MPI ranks. See ComputeGlobalSpectralRadiusEstimate().
     *
     * dtau = pseudo_cfl / lambda_max   (pseudo_cfl is now a true Courant number)
     * dtau_min/max are derived from the user-specified min/max_pseudo_cfl bounds.
     *
     * The BDF2 lower bound inside the helper ensures lambda_max >= 1.5/dt, so dtau
     * remains finite even when the velocity field is zero (startup / stagnant ICs). */
    ierr = ComputeGlobalSpectralRadiusEstimate(user, block_number, dt, &lambda_max); CHKERRQ(ierr);
    dtau_min = simCtx->min_pseudo_cfl / lambda_max;
    dtau_max = simCtx->max_pseudo_cfl / lambda_max;

    /* Initialise per-block pseudo-time step from the warm-start CFL (carried from last timestep). */
    for (PetscInt bi = 0; bi < block_number; bi++) {
        pseudo_dtau[bi] = cfl / lambda_max;   /* dtau [physical time], NOT a fraction of dt */
    }

    /* --- SHADOW MODE: new operator-scaled stability estimate (diagnostic only) ---
     * Computes the Workstream-A estimate alongside the legacy one and logs a compact
     * comparison. The legacy estimate above continues to drive production dtau; the
     * shadow estimate changes nothing. Enable with -mom_stability_shadow. */
    {
        PetscBool shadow = PETSC_FALSE;
        ierr = PetscOptionsGetBool(NULL, NULL, "-mom_stability_shadow", &shadow, NULL); CHKERRQ(ierr);
        if (shadow) {
            MomStabilityReport rep;
            ierr = ComputeMomentumStabilityEstimate(user, block_number, dt, MOM_STAB_CAND_C, &rep); CHKERRQ(ierr);
            const char *lim = (rep.limiter==MOM_STAB_LIMITER_TIME) ? "time"
                            : (rep.limiter==MOM_STAB_LIMITER_CONVECTION) ? "convection" : "viscosity";
            /* Detailed shadow comparison at DEBUG only (gated by level + function allowlist);
               never printed unconditionally. Enable in tests via the logging controls. */
            LOG_ALLOW(GLOBAL, LOG_DEBUG,
                "Momentum scale [shadow]: legacy_dtau=%.4e new_dtau=%.4e ratio=%.4f limiter=%s | "
                "lambda_legacy=%.4e lambda_new=%.4e (B=%.4e C=%.4e D=%.4e) | "
                "lt=%.4e lc=%.4e lv=%.4e | cell=(%d,%d,%d) blk=%d class=%d onesided=%d\n",
                cfl / lambda_max, cfl / rep.lambda, rep.lambda / lambda_max, lim,
                lambda_max, rep.lambda, rep.lambda_B, rep.lambda_C, rep.lambda_D,
                rep.lambda_t, rep.lambda_c, rep.lambda_v,
                (int)rep.ci, (int)rep.cj, (int)rep.ck, (int)rep.cblock, (int)rep.cclass, (int)rep.one_sided);
        }
    }

    LOG_ALLOW(GLOBAL, LOG_INFO,
        "Dual-time solver: lambda_max=%.4e [1/s]  dtau_init=%.4e  dtau range [%.4e, %.4e]  "
        "CFL range [%.4f, %.4f]  rejection_threshold=%.3f (EMA alpha=%.2f)  "
        "growth=%.3f  reduction=%.3f  max_accepted=%d.\n",
        lambda_max, cfl / lambda_max, dtau_min, dtau_max,
        simCtx->min_pseudo_cfl, simCtx->max_pseudo_cfl,
        simCtx->mom_dt_jameson_residual_norm_noise_allowance_factor,
        simCtx->mom_ratio_ema_alpha,
        simCtx->pseudo_cfl_growth_factor, simCtx->pseudo_cfl_reduction_factor,
        max_pseudo_steps);

    // --- 2. Main Pseudo-Time Iteration Loop ---
    pseudo_iter = 0;
    PetscBool residual_convergence_enabled =
        (PetscBool)(simCtx->mom_resid_atol > 0.0 || simCtx->mom_resid_rtol > 0.0);
    PetscBool converged = PETSC_FALSE;
    PetscBool last_trial_nonfinite = PETSC_FALSE;
    /* Rejected iterations no longer consume the accepted-iteration budget.
       A hard cap of 3x max_pseudo_steps guards against infinite rejection loops. */
    const PetscInt max_total_attempts = max_pseudo_steps * 3;
    while (!converged && accepted_iter < max_pseudo_steps && pseudo_iter < max_total_attempts)
    {
        pseudo_iter++;
        force_restart = PETSC_FALSE;

        // Every attempt starts from the last globally accepted state.
        for (PetscInt bi = 0; bi < block_number; bi++) {
            ierr = VecCopy(user[bi].pUcont, user[bi].Ucont); CHKERRQ(ierr);
            ierr = VecCopy(pRhs[bi], user[bi].Rhs); CHKERRQ(ierr);
            const char *staggered_fields[] = {"Ucont"};
            ierr = SynchronizePeriodicStaggeredFields(&user[bi], 1, staggered_fields); CHKERRQ(ierr);
            ierr = ApplyBoundaryConditions(&user[bi]); CHKERRQ(ierr);
        }
	
        for (PetscInt bi = 0; bi < block_number; bi++) {
            
            // === 4-Stage Jameson RK Smoothing Loop ===
            for (istage = 0; istage < 4; istage++) {
                
                LOG_ALLOW(GLOBAL,LOG_TRACE," Pseudo-Iter: %d | RK-Stage: %d\n", pseudo_iter, istage);
                
                /* RK Update: U_new = U_old + (dtau * alpha) * Residual
                 * pseudo_dtau[bi] IS dtau [physical time] — no further multiplication by dt.
                 * (In the old scheme dtau = pseudo_cfl * dt; here dtau = pseudo_cfl / lambda_max.) */
                ierr = VecWAXPY(user[bi].Ucont,
                                pseudo_dtau[bi] * alfa[istage],
                                user[bi].Rhs,
                                user[bi].pUcont); CHKERRQ(ierr);

                // Sync Ghosts & Re-apply BCs for intermediate stage
                const char *staggered_fields[] = {"Ucont"};
                ierr = SynchronizePeriodicStaggeredFields(&user[bi], 1, staggered_fields); CHKERRQ(ierr);
                
                ierr = ApplyBoundaryConditions(&user[bi]); CHKERRQ(ierr);

                // --- Re-calculate Total Residual for next stage ---
                ierr = ComputeTotalResidual(&user[bi]); CHKERRQ(ierr);
      
            } // End RK Stages

            // Immersed boundary interpolation (if enabled)

            // === Convergence Metrics Calculation ===
            
            // Calculate dU = U_current - U_backup
            ierr = VecWAXPY(user[bi].dUcont, -1.0, user[bi].pUcont, user[bi].Ucont); CHKERRQ(ierr);

            // Compute Infinity Norms
            ierr = VecNorm(user[bi].dUcont, NORM_INFINITY, &delta_sol_norm_curr[bi]); CHKERRQ(ierr);
            ierr = VecNorm(user[bi].Rhs, NORM_INFINITY, &resid_norm_curr[bi]); CHKERRQ(ierr);
	    
            // Normalize relative metrics
            if (pseudo_iter == 1) {
                delta_sol_norm_init[bi] = delta_sol_norm_curr[bi];
                delta_sol_rel_curr[bi]  = 1.0;
                resid_rel_curr[bi]      = 1.0;
                // resid_norm_init[bi] set correctly before the loop — do not overwrite
            } else {
                if (delta_sol_norm_init[bi] > 1.0e-10) 
                    delta_sol_rel_curr[bi] = delta_sol_norm_curr[bi] / delta_sol_norm_init[bi];
                else 
                    delta_sol_rel_curr[bi] = 0.0;

                if(resid_norm_init[bi] > 1.0e-10) 
                    resid_rel_curr[bi] = resid_norm_curr[bi] / resid_norm_init[bi];
                else 
                    resid_rel_curr[bi] = 0.0;
	        }
	      
            /* Pre-compute per-block step-to-step trial ratio for later file logging. */
            {
                const PetscReal resid_floor_log = 1.0e-30;
                if (resid_norm_prev[bi] > resid_floor_log)
                    trial_ratio_log[bi] = resid_norm_curr[bi] / resid_norm_prev[bi];
                else if (resid_norm_curr[bi] <= resid_floor_log)
                    trial_ratio_log[bi] = 0.0;
                else
                    trial_ratio_log[bi] = PETSC_MAX_REAL;
            }
        } // End loop over blocks

        // --- Update Global Convergence Criteria ---
        global_norm_delta = -1.0e20; 
        global_rel_delta  = -1.0e20; 
        global_norm_resid = -1.0e20;
        global_rel_resid  = -1.0e20;
        
        for (PetscInt bi = 0; bi < block_number; bi++) {
            global_norm_delta = PetscMax(delta_sol_norm_curr[bi], global_norm_delta);
            global_rel_delta  = PetscMax(delta_sol_rel_curr[bi],  global_rel_delta);
            global_norm_resid = PetscMax(resid_norm_curr[bi],      global_norm_resid);
            global_rel_resid  = PetscMax(resid_rel_curr[bi],      global_rel_resid);
        }
        ierr = PetscTime(&te); CHKERRQ(ierr);
        cput = te - ts;
        LOG_ALLOW(GLOBAL, LOG_INFO, "  Pseudo-Iter(k) %d: |dUk|=%e, |dUk|/|dU0| = %e, |Rk|/|R0| = %e, CPU=%.2fs\n",
              pseudo_iter, global_norm_delta, global_rel_delta, global_rel_resid, cput);

        // === Adaptive Pseudo-CFL Trial Acceptance and Rollback ===
        const PetscReal resid_floor = 1.0e-30;
        PetscReal global_trial_ratio = 0.0;
        PetscBool global_nonfinite = PETSC_FALSE;
        for (PetscInt bi = 0; bi < block_number; bi++) {
            PetscReal ratio;
            if (resid_norm_prev[bi] > resid_floor) {
                ratio = resid_norm_curr[bi] / resid_norm_prev[bi];
            } else if (resid_norm_curr[bi] <= resid_floor) {
                ratio = 0.0;
            } else {
                ratio = PETSC_MAX_REAL;
            }
            global_trial_ratio = PetscMax(global_trial_ratio, ratio);
            global_nonfinite = (PetscBool)(global_nonfinite ||
                PetscIsInfOrNanReal(delta_sol_norm_curr[bi]) ||
                PetscIsInfOrNanReal(resid_norm_curr[bi]) ||
                PetscIsInfOrNanReal(ratio));
        }

        PetscReal reduced_trial_ratio;
        PetscMPIInt local_nonfinite = global_nonfinite ? 1 : 0;
        PetscMPIInt reduced_nonfinite = 0;
        ierr = MPI_Allreduce(&global_trial_ratio, &reduced_trial_ratio, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD); CHKERRQ(ierr);
        ierr = MPI_Allreduce(&local_nonfinite, &reduced_nonfinite, 1, MPI_INT, MPI_LOR, PETSC_COMM_WORLD); CHKERRQ(ierr);
        global_trial_ratio = reduced_trial_ratio;
        global_nonfinite = reduced_nonfinite ? PETSC_TRUE : PETSC_FALSE;
        last_trial_nonfinite = global_nonfinite;

        /* EMA-smooth the trial ratio to tolerate occasional non-monotonic residual bumps.
           Non-finite trials bypass smoothing and always trigger rollback. */
        if (!global_nonfinite) {
            smoothed_trial_ratio = simCtx->mom_ratio_ema_alpha * global_trial_ratio
                                 + (1.0 - simCtx->mom_ratio_ema_alpha) * smoothed_trial_ratio;
        }
        LOG_ALLOW(GLOBAL, LOG_DEBUG,
                  "    [k=%d] raw_ratio=%.4e smoothed_ratio=%.4e (threshold=%.3f) | "
                  "|R_prev|=%.6e | |R_curr|=%.6e | CFL=%.6f\n",
                  pseudo_iter, global_trial_ratio, smoothed_trial_ratio,
                  simCtx->mom_dt_jameson_residual_norm_noise_allowance_factor,
                  resid_norm_prev[0], resid_norm_curr[0], pseudo_dtau[0]);

        if (simCtx->no_pseudo_cfl_backtrack) {
            /* Diagnostic mode: commit every finite trial; only non-finite triggers rollback. */
            force_restart = global_nonfinite;
        } else {
            force_restart = (PetscBool)(global_nonfinite ||
                smoothed_trial_ratio > simCtx->mom_dt_jameson_residual_norm_noise_allowance_factor);
        }

        if (force_restart) {
            /* old_dtau: dtau used for this (rejected) trial [physical time].
             * next_dtau: reduced dtau for the retry, clamped to dtau_min.
             * Both are physical-time values; the effective CFL = dtau * lambda_max. */
            PetscReal old_dtau  = pseudo_dtau[0];
            PetscReal next_dtau = PetscMax(dtau_min, old_dtau * simCtx->pseudo_cfl_reduction_factor);
            rejected_iter++;
            recovery_streak = 0;
            for (PetscInt bi = 0; bi < block_number; bi++) {
                ierr = VecCopy(user[bi].pUcont, user[bi].Ucont); CHKERRQ(ierr);
                ierr = VecCopy(pRhs[bi], user[bi].Rhs); CHKERRQ(ierr);
                const char *staggered_fields[] = {"Ucont"};
                ierr = SynchronizePeriodicStaggeredFields(&user[bi], 1, staggered_fields); CHKERRQ(ierr);
                ierr = ApplyBoundaryConditions(&user[bi]); CHKERRQ(ierr);
                pseudo_dtau[bi] = next_dtau;
            }
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "  Trial %d REJECTED (raw_ratio=%.4e, smoothed=%.4e, nonfinite=%d); "
                      "dtau %.4e -> %.4e (cfl_eff %.4f -> %.4f)%s\n",
                      pseudo_iter, global_trial_ratio, smoothed_trial_ratio, (int)global_nonfinite,
                      old_dtau, next_dtau, old_dtau * lambda_max, next_dtau * lambda_max,
                      (old_dtau == next_dtau) ? " [AT FLOOR — no dtau reduction]" : "");
            /* Log rejected trial to file before rolling back. */
            if (!rank) {
                for (PetscInt bi = 0; bi < block_number; bi++) {
                    FILE *f;
                    char filen[PETSC_MAX_PATH_LEN + 128];
                    ierr = PetscSNPrintf(filen, sizeof(filen),
                        "%s/Momentum_Solver_Convergence_History_Block_%1d.log",
                        simCtx->log_dir, bi); CHKERRQ(ierr);
                    if (simCtx->step == simCtx->StartStep + 1 && pseudo_iter == 1 && !simCtx->continueMode)
                        f = fopen(filen, "w");
                    else
                        f = fopen(filen, "a");
                    if (simCtx->continueMode && simCtx->step == simCtx->StartStep + 1 && pseudo_iter == 1)
                        PetscFPrintf(PETSC_COMM_WORLD, f, "# Continuation from step %" PetscInt_FMT "\n", simCtx->StartStep);
                    /* Log format: dtau [physical time] + cfl_eff [dimensionless Courant number = dtau*lambda_max] */
                    PetscFPrintf(PETSC_COMM_WORLD, f,
                        "Step: %d | PseudoIter(k): %d | dtau: %.6e | cfl_eff: %.4f | |dUk|: %le | "
                        "|dUk|/|dU0|: %le | |Rk|: %le | |Rk|/|R0|: %le | "
                        "trial_ratio: %le | smoothed_ratio: %le | status: rejected | "
                        "dtau_after: %.6e | cfl_eff_after: %.4f\n",
                        (int)ti, (int)pseudo_iter, old_dtau, old_dtau * lambda_max,
                        delta_sol_norm_curr[bi], delta_sol_rel_curr[bi],
                        resid_norm_curr[bi], resid_rel_curr[bi],
                        trial_ratio_log[bi], smoothed_trial_ratio,
                        next_dtau, next_dtau * lambda_max);
                    fclose(f);
                }
            }
            if (old_dtau <= dtau_min) {
                if (global_nonfinite) {
                    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_CONV_FAILED,
                            "Momentum solver produced a non-finite trial at minimum pseudo-CFL.");
                }
                /* Finite rejection at the dtau floor: reducing further is impossible.
                   Break instead of retrying bit-identically until max_pseudo_steps. */
                LOG_ALLOW(GLOBAL, LOG_WARNING,
                          "  Trial %d REJECTED (ratio=%.4e) at minimum dtau (%.4e, cfl_eff=%.4f) "
                          "with no further reduction possible; breaking retry loop.\n",
                          pseudo_iter, global_trial_ratio, old_dtau, old_dtau * lambda_max);
                break;
            }
            continue;
        }

        accepted_iter++;
        last_trial_nonfinite = PETSC_FALSE;
        last_accepted_resid  = resid_norm_curr[0]; /* track for final summary */
        for (PetscInt bi = 0; bi < block_number; bi++) {
            ierr = VecCopy(user[bi].Ucont, user[bi].pUcont); CHKERRQ(ierr);
            ierr = VecCopy(user[bi].Rhs, pRhs[bi]); CHKERRQ(ierr);
            resid_norm_prev[bi] = resid_norm_curr[bi];
            delta_sol_norm_prev[bi] = delta_sol_norm_curr[bi];
        }

        /* old_dtau: dtau used for this (accepted) trial.
         * next_dtau: dtau for the next trial, grown/reduced based on convergence rate.
         * Bounds (dtau_min, dtau_max) are derived from simCtx->min/max_pseudo_cfl / lambda_max. */
        PetscReal old_dtau  = pseudo_dtau[0]; /* save before any update */
        PetscReal next_dtau = old_dtau;
        if (global_trial_ratio < 0.90) {
            /* Residual decreasing fast: grow dtau to accelerate convergence. */
            next_dtau *= simCtx->pseudo_cfl_growth_factor;
            recovery_streak = 0;
        } else if (global_trial_ratio <= 1.0) {
            /* Residual barely decreasing: grow only after 3 clean consecutive trials. */
            recovery_streak++;
            if (recovery_streak >= 3) {
                next_dtau *= simCtx->pseudo_cfl_growth_factor;
                recovery_streak = 0;
            }
        } else {
            /* Accepted but residual slightly grew (within EMA noise allowance): reduce dtau. */
            next_dtau *= simCtx->pseudo_cfl_reduction_factor;
            recovery_streak = 0;
        }
        next_dtau = PetscMin(next_dtau, dtau_max);
        next_dtau = PetscMax(next_dtau, dtau_min);
        for (PetscInt bi = 0; bi < block_number; bi++) pseudo_dtau[bi] = next_dtau;

        LOG_ALLOW(GLOBAL, LOG_DEBUG,
                  "  Trial %d ACCEPTED (raw_ratio=%.4e, smoothed=%.4e); |dU|=%.6e | "
                  "dtau %.4e -> %.4e (cfl_eff %.4f -> %.4f)\n",
                  pseudo_iter, global_trial_ratio, smoothed_trial_ratio, global_norm_delta,
                  old_dtau, next_dtau, old_dtau * lambda_max, next_dtau * lambda_max);

        /* --- Post-decision file logging (one row per accepted trial) --- */
        if (!rank) {
            for (PetscInt bi = 0; bi < block_number; bi++) {
                FILE *f;
                char filen[PETSC_MAX_PATH_LEN + 128];
                ierr = PetscSNPrintf(filen, sizeof(filen),
                    "%s/Momentum_Solver_Convergence_History_Block_%1d.log",
                    simCtx->log_dir, bi); CHKERRQ(ierr);
                if (simCtx->step == simCtx->StartStep + 1 && pseudo_iter == 1 && !simCtx->continueMode)
                    f = fopen(filen, "w");
                else
                    f = fopen(filen, "a");
                if (simCtx->continueMode && simCtx->step == simCtx->StartStep + 1 && pseudo_iter == 1)
                    PetscFPrintf(PETSC_COMM_WORLD, f, "# Continuation from step %" PetscInt_FMT "\n", simCtx->StartStep);
                /* Log format: dtau [physical time] + cfl_eff [dimensionless Courant number = dtau*lambda_max] */
                PetscFPrintf(PETSC_COMM_WORLD, f,
                    "Step: %d | PseudoIter(k): %d | dtau: %.6e | cfl_eff: %.4f | |dUk|: %le | "
                    "|dUk|/|dU0|: %le | |Rk|: %le | |Rk|/|R0|: %le | "
                    "trial_ratio: %le | smoothed_ratio: %le | status: accepted | "
                    "dtau_after: %.6e | cfl_eff_after: %.4f\n",
                    (int)ti, (int)pseudo_iter, old_dtau, old_dtau * lambda_max,
                    delta_sol_norm_curr[bi], delta_sol_rel_curr[bi],
                    resid_norm_curr[bi], resid_rel_curr[bi],
                    trial_ratio_log[bi], smoothed_trial_ratio,
                    next_dtau, next_dtau * lambda_max);
                fclose(f);
            }
        }

        PetscBool update_pass;
        PetscBool residual_pass = PETSC_TRUE;
        if (residual_convergence_enabled) {
            update_pass = (PetscBool)(global_norm_delta <= tol_abs_delta && global_rel_delta <= tol_rtol_delta);
            residual_pass = (PetscBool)(
                (simCtx->mom_resid_atol > 0.0 && global_norm_resid <= simCtx->mom_resid_atol) ||
                (simCtx->mom_resid_rtol > 0.0 && global_rel_resid <= simCtx->mom_resid_rtol));
        } else {
            update_pass = (PetscBool)(global_norm_delta <= tol_abs_delta && global_rel_delta <= tol_rtol_delta);
        }
        converged = (PetscBool)(update_pass && residual_pass);

        if (block_number > 1) {
             // ierr = Block_Interface_U(user); CHKERRQ(ierr);
        }
    } // End while loop

    if (last_trial_nonfinite) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_CONV_FAILED,
                "Momentum solver exhausted its attempt limit while recovering from a non-finite trial.");
    }

    /* Convert the final dtau back to a CFL number for the warm-start of the next physical timestep.
     * simCtx->pseudo_cfl stores a dimensionless Courant number (dtau * lambda_max); the next
     * timestep's lambda_max may differ as the flow evolves, so storing the CFL number (not dtau)
     * ensures the warm-start adapts correctly to the new spectral radius. */
    PetscReal next_dtau_start  = pseudo_dtau[0];
    PetscReal next_cfl_warmstart = next_dtau_start * lambda_max;  /* CFL = dtau * lambda_max */

    if (!rank) {
        LOG_ALLOW(GLOBAL, LOG_INFO,
                  "  Step %d finished: dtau=%.4e  cfl_eff=%.4f  lambda_max=%.4e [1/s]. "
                  "Next step warm-starts at cfl=%.4f.\n",
                  simCtx->step, next_dtau_start, next_cfl_warmstart, lambda_max, next_cfl_warmstart);
    }

    /* Store the CFL number (not dtau) so it remains meaningful across physical timesteps. */
    simCtx->pseudo_cfl        = next_cfl_warmstart;
    simCtx->mom_last_lambda_max = lambda_max;    

    /* Expose the last accepted finite state. On converged or floor-break exits pUcont==Ucont
       already (VecCopy would be an identity), but ghosts may have gone stale; always re-sync
       and re-apply BCs. Only restore the solution vector itself when nothing was ever accepted. */
    if (accepted_iter == 0) {
        for (PetscInt bi = 0; bi < block_number; bi++) {
            ierr = VecCopy(user[bi].pUcont, user[bi].Ucont); CHKERRQ(ierr);
            ierr = VecCopy(pRhs[bi], user[bi].Rhs); CHKERRQ(ierr);
        }
    }
    for (PetscInt bi = 0; bi < block_number; bi++) {
        const char *staggered_fields[] = {"Ucont"};
        ierr = SynchronizePeriodicStaggeredFields(&user[bi], 1, staggered_fields); CHKERRQ(ierr);
        ierr = ApplyBoundaryConditions(&user[bi]); CHKERRQ(ierr);
    }

    simCtx->mom_last_converged = converged;
    if (!converged) {
        PetscPrintf(PETSC_COMM_WORLD,
            "[WARNING] Momentum solver step %d: reached %d total attempts (%d accepted, %d rejected) "
            "without convergence; continuing from last accepted finite state.\n",
            (int)ti, pseudo_iter, accepted_iter, rejected_iter);
    }
    if (accepted_iter == 0) {
        PetscPrintf(PETSC_COMM_WORLD,
            "[WARNING] Momentum solver step %d: no pseudo-time trials were accepted; "
            "retaining physical-step entry state.\n", (int)ti);
    }
    LOG_ALLOW(GLOBAL, LOG_INFO,
              "Momentum solver finished: %d attempts (%d accepted, %d rejected) of %d max accepted / %d hard cap, "
              "converged=%s, last accepted |R|=%.6e, last accepted |dU|=%.6e, "
              "next_dtau=%.4e (cfl=%.4f).\n",
              pseudo_iter, accepted_iter, rejected_iter, max_pseudo_steps, max_total_attempts,
              converged ? "yes" : "no",
              (accepted_iter > 0) ? last_accepted_resid : resid_norm_init[0],
              (accepted_iter > 0) ? delta_sol_norm_prev[0] : 0.0,
              next_dtau_start, next_cfl_warmstart);

    // --- Final Cleanup ---
    for (PetscInt bi = 0; bi < block_number; bi++) {
        ierr = VecDestroy(&user[bi].Rhs); CHKERRQ(ierr);
        ierr = VecDestroy(&user[bi].dUcont); CHKERRQ(ierr);
        ierr = VecDestroy(&user[bi].pUcont); CHKERRQ(ierr);
        ierr = VecDestroy(&pRhs[bi]); CHKERRQ(ierr);
    }
    ierr = PetscFree(pRhs); CHKERRQ(ierr);

    ierr = PetscFree2(delta_sol_norm_init, delta_sol_norm_prev);CHKERRQ(ierr);
    ierr = PetscFree2(delta_sol_norm_curr, delta_sol_rel_curr);CHKERRQ(ierr);
    ierr = PetscFree2(resid_norm_init,     resid_norm_prev);CHKERRQ(ierr);
    ierr = PetscFree2(resid_norm_curr,     pseudo_dtau); CHKERRQ(ierr);
    ierr = PetscFree(resid_rel_curr); CHKERRQ(ierr);
    ierr = PetscFree(trial_ratio_log); CHKERRQ(ierr);
    
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
