/**
 * @file test_momentum_newton_boundary_fixedpoint.c
 * @brief Opt-in integration regression for the Newton--Krylov deterministic
 *        Cartesian-seed correction on the production-sized straight duct.
 *
 * The diagnostic phase established that the production Newton--Krylov residual
 * must reconstruct Ucat/lUcat from the current SNES trial vector X before the
 * first boundary sweep (see src/momentum_newton_krylov.c). That correction is
 * now in production. This file locks the resulting guarantees in place using the
 * *actual production* residual callback (included via a privately renamed copy of
 * the production translation unit) -- it deliberately does NOT keep a duplicate
 * test-local seeded callback, so it exercises exactly what ships:
 *
 *   1. residual purity at the reproduced physical step-2 state (immediate and
 *      after real PETSc MFFD products);
 *   2. a complete production step-2 solve from the true projected step-1 state
 *      using PETSc's default classical Gram--Schmidt;
 *   3. boundary-map compatibility: the converged three-pass solution also zeros
 *      the 24-pass outlet residual;
 *   4. clean pressure projection of the converged momentum state.
 *
 * The heavy fixture keeps this target out of `make check`; the cheap
 * seed-removal detector lives in unit-newton-krylov (default suite). Central
 * checks run on one and four MPI ranks.
 */

#include "test_support.h"
#include "initialcondition.h"
#include "runloop.h"
#include "solvers.h"

#define MomentumSolver_NewtonKrylov MomentumSolver_NewtonKrylov_BoundaryFixedpointPrivateCopy
#include "../../src/momentum_newton_krylov.c"
#undef MomentumSolver_NewtonKrylov

/* ======================================================================== */
/*  Persistent-state snapshot machinery (solve-entry restore).              */
/* ======================================================================== */

#define FP_FIELD_COUNT 16
#define FP_SCALAR_COUNT 10

typedef struct {
    Vec value[FP_FIELD_COUNT];
    PetscReal scalar[FP_SCALAR_COUNT];
} FpSnapshot;

/** @brief Collects the persistent vectors observed by the residual/boundary path. */
static void FpGetVectors(UserCtx *user, Vec field[FP_FIELD_COUNT])
{
    field[0] = user->Ucont; field[1] = user->lUcont;
    field[2] = user->Ucat; field[3] = user->lUcat;
    field[4] = user->P; field[5] = user->lP;
    field[6] = user->Bcs.Ubcs; field[7] = user->Bcs.Uch;
    field[8] = user->Rhs;
    field[9] = user->Ucont_o; field[10] = user->lUcont_o;
    field[11] = user->Ucont_rm1; field[12] = user->lUcont_rm1;
    field[13] = user->Ucat_o; field[14] = user->P_o;
    field[15] = user->CellFieldAtCorner;
}

/** @brief Collects the persistent scalar diagnostics maintained by boundary processing. */
static void FpGetScalars(UserCtx *user, PetscReal scalar[FP_SCALAR_COUNT])
{
    SimCtx *s = user->simCtx;
    scalar[0] = s->FluxInSum; scalar[1] = s->FluxOutSum;
    scalar[2] = s->FarFluxInSum; scalar[3] = s->FarFluxOutSum;
    scalar[4] = s->Fluxsum; scalar[5] = s->AreaInSum;
    scalar[6] = s->AreaOutSum; scalar[7] = s->bulkVelocityCorrection;
    scalar[8] = s->boundaryVelocityCorrection; scalar[9] = user->FluxIntpSum;
}

/** @brief Deep-copies the audited persistent state into a test-owned snapshot. */
static PetscErrorCode FpCapture(UserCtx *user, FpSnapshot *snap)
{
    Vec field[FP_FIELD_COUNT];
    PetscFunctionBeginUser;
    memset(snap, 0, sizeof(*snap));
    FpGetVectors(user, field);
    for (PetscInt i = 0; i < FP_FIELD_COUNT; ++i) {
        if (field[i]) {
            PetscCall(VecDuplicate(field[i], &snap->value[i]));
            PetscCall(VecCopy(field[i], snap->value[i]));
        }
    }
    FpGetScalars(user, snap->scalar);
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Restores a previously captured persistent-state snapshot. */
static PetscErrorCode FpRestore(UserCtx *user, const FpSnapshot *snap)
{
    Vec field[FP_FIELD_COUNT];
    SimCtx *s = user->simCtx;
    PetscFunctionBeginUser;
    FpGetVectors(user, field);
    for (PetscInt i = 0; i < FP_FIELD_COUNT; ++i)
        if (field[i] && snap->value[i]) PetscCall(VecCopy(snap->value[i], field[i]));
    s->FluxInSum = snap->scalar[0]; s->FluxOutSum = snap->scalar[1];
    s->FarFluxInSum = snap->scalar[2]; s->FarFluxOutSum = snap->scalar[3];
    s->Fluxsum = snap->scalar[4]; s->AreaInSum = snap->scalar[5];
    s->AreaOutSum = snap->scalar[6]; s->bulkVelocityCorrection = snap->scalar[7];
    s->boundaryVelocityCorrection = snap->scalar[8]; user->FluxIntpSum = snap->scalar[9];
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Releases all vectors owned by one persistent-state snapshot. */
static PetscErrorCode FpDestroy(FpSnapshot *snap)
{
    PetscFunctionBeginUser;
    for (PetscInt i = 0; i < FP_FIELD_COUNT; ++i) PetscCall(VecDestroy(&snap->value[i]));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/* ======================================================================== */
/*  Collective norm helpers.                                                */
/* ======================================================================== */

/** @brief Computes collective global difference norms for global or rank-local vectors. */
static PetscErrorCode CollectiveDifference(Vec a, Vec b, PetscReal *norm_inf, PetscReal *norm_2)
{
    Vec delta = NULL;
    MPI_Comm comm;
    PetscMPIInt comm_size, world_size;
    PetscReal local_inf, local_2;
    PetscFunctionBeginUser;
    PetscCall(VecDuplicate(a, &delta));
    PetscCall(VecWAXPY(delta, -1.0, a, b));
    PetscCall(VecNorm(delta, NORM_INFINITY, &local_inf));
    PetscCall(VecNorm(delta, NORM_2, &local_2));
    comm = PetscObjectComm((PetscObject)a);
    PetscCallMPI(MPI_Comm_size(comm, &comm_size));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &world_size));
    if (comm_size == 1 && world_size > 1) {
        PetscReal square = local_2 * local_2;
        PetscCallMPI(MPI_Allreduce(&local_inf, norm_inf, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD));
        PetscCallMPI(MPI_Allreduce(&square, norm_2, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD));
        *norm_2 = PetscSqrtReal(*norm_2);
    } else {
        *norm_inf = local_inf;
        *norm_2 = local_2;
    }
    PetscCall(VecDestroy(&delta));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Finds the MPI rank containing the largest local entry difference. */
static PetscErrorCode DifferenceMaxRank(Vec a, Vec b, PetscMPIInt *max_rank)
{
    const PetscScalar *aa = NULL, *bb = NULL;
    PetscInt nlocal;
    PetscReal local_max = 0.0, *all_max = NULL;
    PetscMPIInt size;
    PetscFunctionBeginUser;
    PetscCall(VecGetLocalSize(a, &nlocal));
    PetscCall(VecGetArrayRead(a, &aa)); PetscCall(VecGetArrayRead(b, &bb));
    for (PetscInt i = 0; i < nlocal; ++i) local_max = PetscMax(local_max, PetscAbsScalar(bb[i] - aa[i]));
    PetscCall(VecRestoreArrayRead(b, &bb)); PetscCall(VecRestoreArrayRead(a, &aa));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    PetscCall(PetscMalloc1(size, &all_max));
    PetscCallMPI(MPI_Allgather(&local_max, 1, MPIU_REAL, all_max, 1, MPIU_REAL, PETSC_COMM_WORLD));
    *max_rank = 0;
    for (PetscMPIInt r = 1; r < size; ++r) if (all_max[r] > all_max[*max_rank]) *max_rank = r;
    PetscCall(PetscFree(all_max));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/* ======================================================================== */
/*  Deterministic Cartesian seed and an N-pass residual for compatibility.  */
/* ======================================================================== */

/**
 * @brief Reconstructs Ucat/lUcat deterministically from X, mirroring the seed
 *        now performed by the production callback.
 *
 * Used only by the N-pass compatibility residual (Experiment 3), which needs a
 * boundary-pass count other than the production default of one call (three
 * passes). For a single call it reproduces the production callback exactly.
 */
static PetscErrorCode SeedCartesianFromX(UserCtx *user, Vec X)
{
    const char *staggered[] = {"Ucont"};
    const char *cell[] = {"Ucat"};
    PetscFunctionBeginUser;
    PetscCall(VecCopy(X, user->Ucont));
    PetscCall(SynchronizePeriodicStaggeredFields(user, 1, staggered));
    PetscCall(Contra2Cart(user));
    PetscCall(SynchronizePeriodicCellFields(user, 1, cell));
    PetscCall(UpdateLocalGhosts(user, "Ucat"));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/**
 * @brief Deterministic residual at X using an explicit number of boundary-condition
 *        calls (each three internal passes), for the boundary-map compatibility test.
 */
static PetscErrorCode ResidualNCalls(MomentumNewtonKrylovContext *ctx, Vec X,
                                     PetscInt boundary_calls, Vec F)
{
    UserCtx *user = ctx->user;
    PetscFunctionBeginUser;
    PetscCall(SeedCartesianFromX(user, X));
    for (PetscInt i = 0; i < boundary_calls; ++i) PetscCall(ApplyBoundaryConditions(user));
    PetscCall(ComputeTotalResidual(user));
    PetscCall(VecCopy(user->Rhs, F));
    PetscCall(VecScale(F, -1.0));
    PetscCall(MomentumNewtonKrylov_ApplyConstraints(ctx, X, F));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/* ======================================================================== */
/*  Newton solve configuration and monitor.                                 */
/* ======================================================================== */

/** @brief Installs the production Newton, line-search, and GMRES options. */
static PetscErrorCode ConfigureNewtonOptions(void)
{
    PetscFunctionBeginUser;
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_snes_type", "newtonls"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_snes_atol", "1e-10"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_snes_rtol", "1e-8"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_snes_stol", "1e-12"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_snes_max_it", "25"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_snes_linesearch_type", "bt"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_ksp_type", "gmres"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_ksp_atol", "1e-10"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_ksp_rtol", "1e-6"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_ksp_max_it", "400"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_ksp_gmres_restart", "80"));
    PetscCall(PetscOptionsSetValue(NULL, "-mom_nk_pc_type", "none"));
    PetscFunctionReturn(PETSC_SUCCESS);
}

typedef struct { const char *label; PetscReal *initial_norm; } SolveMonitorCtx;

/** @brief Records per-Newton nonlinear norm, KSP iterations/reason, and accepted lambda. */
static PetscErrorCode SolveMonitor(SNES snes, PetscInt it, PetscReal norm, void *vctx)
{
    SolveMonitorCtx *m = (SolveMonitorCtx *)vctx;
    SNESLineSearch ls = NULL;
    KSP ksp = NULL;
    KSPConvergedReason kr = KSP_CONVERGED_ITERATING;
    PetscReal lambda = -1.0;
    PetscInt cum_ksp = 0, step_ksp = 0;
    PetscFunctionBeginUser;
    if (it == 0 && m->initial_norm) *m->initial_norm = norm;
    PetscCall(SNESGetLinearSolveIterations(snes, &cum_ksp));
    PetscCall(SNESGetKSP(snes, &ksp));
    PetscCall(KSPGetConvergedReason(ksp, &kr));
    PetscCall(KSPGetIterationNumber(ksp, &step_ksp));
    PetscCall(SNESGetLineSearch(snes, &ls));
    if (it > 0) PetscCall(SNESLineSearchGetLambda(ls, &lambda));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "FULLSOLVE %s newton=%d nonlinear_norm=% .16e last_ksp_its=%d cumulative_ksp_its=%d "
        "last_ksp_reason=%d accepted_lambda=% .6e\n",
        m->label, (int)it, (double)norm, (int)step_ksp, (int)cum_ksp, (int)kr, (double)lambda));
    PetscFunctionReturn(PETSC_SUCCESS);
}

typedef struct {
    SNESConvergedReason reason;
    PetscInt  its, fevals, ksp_total;
    PetscReal initial_norm, final_norm;
    PetscBool converged;
    Vec       Xfinal;
} SolveResult;

/**
 * @brief Creates a test-owned matrix-free SNES driven by the *production* residual callback.
 * @param use_mgs When false (production default), leaves GMRES orthogonalization at
 *        PETSc's default classical Gram--Schmidt -- the path the production solver uses.
 */
static PetscErrorCode CreateProductionSNES(UserCtx *user, MomentumNewtonKrylovContext *pctx,
                                           PetscBool use_mgs, SNES *snes, Mat *J, KSP *ksp)
{
    PC pc = NULL;
    PetscFunctionBeginUser;
    PetscCall(SNESCreate(PETSC_COMM_WORLD, snes));
    PetscCall(SNESSetOptionsPrefix(*snes, "mom_nk_"));
    PetscCall(SNESSetType(*snes, SNESNEWTONLS));
    PetscCall(SNESSetDM(*snes, user->fda));
    PetscCall(SNESSetFunction(*snes, NULL, MomentumNewtonKrylov_FormResidual, pctx));
    PetscCall(MatCreateSNESMF(*snes, J));
    PetscCall(SNESSetJacobian(*snes, *J, *J, MatMFFDComputeJacobian, NULL));
    PetscCall(SNESGetKSP(*snes, ksp));
    PetscCall(KSPSetType(*ksp, KSPGMRES));
    if (use_mgs) PetscCall(KSPGMRESSetOrthogonalization(*ksp, KSPGMRESModifiedGramSchmidtOrthogonalization));
    PetscCall(KSPGMRESSetRestart(*ksp, 80));
    PetscCall(KSPGetPC(*ksp, &pc));
    PetscCall(PCSetType(pc, PCNONE));
    PetscCall(SNESSetFromOptions(*snes));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Runs one complete production-callback step-2 solve from the restored entry state. */
static PetscErrorCode RunProductionFullSolve(UserCtx *user, const FpSnapshot *entry, Vec X0,
                                             PetscBool use_mgs, const char *label, SolveResult *res)
{
    MomentumNewtonKrylovContext pctx = {user, NULL, PETSC_FALSE, 0.0};
    SNES snes = NULL; Mat J = NULL; KSP ksp = NULL; Vec solution = NULL;
    SolveMonitorCtx m;
    PetscFunctionBeginUser;
    memset(res, 0, sizeof(*res));
    res->initial_norm = -1.0;
    PetscCall(FpRestore(user, entry));
    PetscCall(VecDuplicate(X0, &solution)); PetscCall(VecCopy(X0, solution));
    PetscCall(CreateProductionSNES(user, &pctx, use_mgs, &snes, &J, &ksp));
    m.label = label; m.initial_norm = &res->initial_norm;
    PetscCall(SNESMonitorSet(snes, SolveMonitor, &m, NULL));
    PetscCall(SNESSolve(snes, NULL, solution));
    PetscCall(SNESGetConvergedReason(snes, &res->reason));
    PetscCall(SNESGetIterationNumber(snes, &res->its));
    PetscCall(SNESGetNumberFunctionEvals(snes, &res->fevals));
    PetscCall(SNESGetLinearSolveIterations(snes, &res->ksp_total));
    PetscCall(SNESGetFunctionNorm(snes, &res->final_norm));
    res->converged = (PetscBool)(res->reason > 0);
    PetscCall(VecDuplicate(solution, &res->Xfinal)); PetscCall(VecCopy(solution, res->Xfinal));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "FULLSOLVE_SUMMARY %s reason=%s (%d) newton=%d fevals=%d ksp_total=%d orthogonalization=%s "
        "initial_norm=% .16e final_norm=% .16e -> %s\n",
        label, SNESConvergedReasons[res->reason], (int)res->reason, (int)res->its, (int)res->fevals,
        (int)res->ksp_total, use_mgs ? "modified-GS" : "classical-GS(default)",
        (double)res->initial_norm, (double)res->final_norm, res->converged ? "CONVERGED" : "FAILED"));
    PetscCall(VecDestroy(&solution)); PetscCall(MatDestroy(&J)); PetscCall(SNESDestroy(&snes));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Stops SNES immediately after its first accepted Newton correction. */
static PetscErrorCode StopAfterOne(SNES snes, PetscInt it, PetscReal xnorm, PetscReal gnorm,
                                   PetscReal fnorm, SNESConvergedReason *reason, void *vctx)
{
    PetscReal *norms = (PetscReal *)vctx;
    PetscFunctionBeginUser;
    (void)snes; (void)xnorm; (void)gnorm;
    if (it == 0) norms[0] = fnorm;
    if (it == 1) norms[1] = fnorm;
    *reason = it >= 1 ? SNES_CONVERGED_ITS : SNES_CONVERGED_ITERATING;
    PetscFunctionReturn(PETSC_SUCCESS);
}

/* ======================================================================== */
/*  The regression itself.                                                  */
/* ======================================================================== */

/**
 * @brief Builds the production-sized duct, advances the true physical step 1, and
 *        exercises purity, full-solve convergence, boundary-map compatibility, and
 *        projection using the production Newton--Krylov callback and solver.
 */
static PetscErrorCode TestNewtonKrylovProductionRegression(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN] = "";
    FpSnapshot entry = {0};
    MomentumNewtonKrylovContext pctx;
    Vec X0 = NULL, X1 = NULL, Xfinal = NULL;
    Vec F0 = NULL, F1 = NULL, Fafter = NULL, v = NULL, jv = NULL;
    Vec F3 = NULL, F24 = NULL;
    SolveResult res_cgs = {0};
    PetscReal norms[2] = {0.0, 0.0};
    PetscReal imm_inf = 0.0, imm_2 = 0.0, mffd_inf = 0.0, mffd_2 = 0.0;
    PetscReal f3_2 = 0.0, f24_2 = 0.0, f24_inf = 0.0, f24m3_inf = 0.0, f24m3_2 = 0.0;
    PetscReal max_div = -1.0;
    const char *fields[] = {"Ucont"};
    const char *stag[] = {"Ucont"}; const char *cell[] = {"Ucat"};
    const PetscReal purity_inf_tol = 1e-13, purity_2_tol = 1e-12;
    const PetscReal compat_tol = 1e-8, div_tol = 1e-6;
    PetscMPIInt world, mffd_rank = 0;
    PetscErrorCode prod_ierr;

    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &world));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "\n########## NEWTON-KRYLOV PRODUCTION REGRESSION (ranks=%d) ##########\n", (int)world));
    PetscCall(PicurvBuildMomentumPurityRuntimeContext(NULL, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(InitializeEulerianState(simCtx));
    PetscCall(ConfigureNewtonOptions());

    /* Complete physical step 1 (momentum + Poisson + projection + finalization)
     * through the real run path, advance histories, then enter step 2. */
    simCtx->step = 1; simCtx->ti = simCtx->dt;
    PetscCall(FlowSolver(simCtx));
    PetscCall(UpdateSolverHistoryVectors(user));
    simCtx->step = 2; simCtx->ti = 2.0 * simCtx->dt;

    PetscCall(VecDuplicate(user->Ucont, &user->Rhs));
    pctx.user = user; pctx.history_file = NULL; pctx.have_initial_norm = PETSC_FALSE; pctx.initial_norm = 0.0;

    /* Canonicalize the entry exactly as production MomentumSolver_NewtonKrylov does. */
    PetscCall(SynchronizePeriodicStaggeredFields(user, 1, fields));
    PetscCall(ApplyBoundaryConditions(user));
    PetscCall(FpCapture(user, &entry));
    PetscCall(VecDuplicate(user->Ucont, &X0)); PetscCall(VecCopy(user->Ucont, X0));

    /* -------- 4.1 residual purity at the step-2 / Newton-1 state -------- */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n===== 4.1 production-callback residual purity =====\n"));
    {   /* Reach the accepted Newton-1 state X1 with the production callback + default CGS. */
        SNES snes = NULL; Mat J = NULL; KSP ksp = NULL;
        PetscCall(VecDuplicate(X0, &X1)); PetscCall(VecCopy(X0, X1));
        PetscCall(CreateProductionSNES(user, &pctx, PETSC_FALSE, &snes, &J, &ksp));
        PetscCall(SNESSetConvergenceTest(snes, StopAfterOne, norms, NULL));
        PetscCall(SNESSolve(snes, NULL, X1));
        PetscCall(MatDestroy(&J)); PetscCall(SNESDestroy(&snes));
    }
    PetscCall(VecDuplicate(X1, &F0)); PetscCall(VecDuplicate(X1, &F1));
    PetscCall(VecDuplicate(X1, &Fafter)); PetscCall(VecDuplicate(X1, &v)); PetscCall(VecDuplicate(X1, &jv));

    /* Immediate repeatability with hidden state poisoned between the two calls. */
    PetscCall(MomentumNewtonKrylov_FormResidual(NULL, X1, F0, &pctx));
    simCtx->FluxInSum = 1234.0; simCtx->FluxOutSum = -4321.0;
    simCtx->FarFluxInSum = 77.0; simCtx->FarFluxOutSum = -88.0;
    PetscCall(VecSet(user->Ucat, 7.0)); PetscCall(VecSet(user->lUcat, 7.0));
    PetscCall(MomentumNewtonKrylov_FormResidual(NULL, X1, F1, &pctx));
    PetscCall(CollectiveDifference(F0, F1, &imm_inf, &imm_2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "PURITY immediate inf=% .3e l2=% .3e (tol inf=%.0e l2=%.0e)\n",
        (double)imm_inf, (double)imm_2, (double)purity_inf_tol, (double)purity_2_tol));

    /* Post-MFFD repeatability: real PETSc MFFD products, then reevaluate. */
    {
        SNES snes = NULL; Mat J = NULL; KSP ksp = NULL;
        PetscInt lo, hi; PetscScalar *a = NULL; PetscReal vn;
        PetscCall(VecGetOwnershipRange(v, &lo, &hi));
        PetscCall(VecGetArray(v, &a));
        for (PetscInt i = lo; i < hi; ++i) a[i - lo] = PetscSinReal(0.17 * (PetscReal)(i + 1));
        PetscCall(VecRestoreArray(v, &a));
        PetscCall(VecNorm(v, NORM_2, &vn)); PetscCall(VecScale(v, 1.0 / vn));
        PetscCall(CreateProductionSNES(user, &pctx, PETSC_FALSE, &snes, &J, &ksp));
        PetscCall(MatMFFDSetBase(J, X1, F0));
        for (PetscInt i = 0; i < 40; ++i) PetscCall(MatMult(J, v, jv));
        PetscCall(MomentumNewtonKrylov_FormResidual(NULL, X1, Fafter, &pctx));
        PetscCall(CollectiveDifference(F0, Fafter, &mffd_inf, &mffd_2));
        PetscCall(DifferenceMaxRank(F0, Fafter, &mffd_rank));
        PetscCall(MatDestroy(&J)); PetscCall(SNESDestroy(&snes));
    }
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "PURITY post_mffd inf=% .3e l2=% .3e max_rank=%d\n",
        (double)mffd_inf, (double)mffd_2, (int)mffd_rank));

    /* -------- 4.2 full step-2 solve with default classical Gram-Schmidt -------- */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n===== 4.2 full production step-2 solve (default classical GS) =====\n"));
    PetscCall(RunProductionFullSolve(user, &entry, X0, PETSC_FALSE, "CGS-default", &res_cgs));

    /* Also drive the ACTUAL production solver wrapper (default GS) and take X_final
     * from it. The wrapper requires user->Rhs unallocated on entry. */
    PetscCall(FpRestore(user, &entry));
    PetscCall(VecDestroy(&user->Rhs));
    prod_ierr = MomentumSolver_NewtonKrylov_BoundaryFixedpointPrivateCopy(user, NULL, NULL);
    PetscCall(VecDuplicate(user->Ucont, &user->Rhs));
    PetscCall(VecDuplicate(user->Ucont, &Xfinal)); PetscCall(VecCopy(user->Ucont, Xfinal));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "PROD_SOLVER ierr=%d mom_last_converged=%d\n",
        (int)prod_ierr, (int)simCtx->mom_last_converged));
    {
        PetscReal dinf, d2; PetscMPIInt r;
        PetscCall(CollectiveDifference(res_cgs.Xfinal, Xfinal, &dinf, &d2));
        PetscCall(DifferenceMaxRank(res_cgs.Xfinal, Xfinal, &r));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,
            "PROD_SOLVER test_snes_vs_wrapper inf=% .3e l2=% .3e max_rank=%d\n",
            (double)dinf, (double)d2, (int)r));
    }

    /* -------- 4.3 boundary-map compatibility at the converged solution -------- */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n===== 4.3 boundary-map compatibility (3 vs 24 passes) =====\n"));
    PetscCall(VecDuplicate(Xfinal, &F3)); PetscCall(VecDuplicate(Xfinal, &F24));
    PetscCall(FpRestore(user, &entry));
    PetscCall(ResidualNCalls(&pctx, Xfinal, 1, F3));    /* 3 passes == production callback */
    PetscCall(VecNorm(F3, NORM_2, &f3_2));
    PetscCall(FpRestore(user, &entry));
    PetscCall(ResidualNCalls(&pctx, Xfinal, 8, F24));   /* 24 passes */
    PetscCall(VecNorm(F24, NORM_2, &f24_2));
    PetscCall(VecNorm(F24, NORM_INFINITY, &f24_inf));
    PetscCall(CollectiveDifference(F3, F24, &f24m3_inf, &f24m3_2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "COMPAT F3_l2=% .6e F24_l2=% .6e F24_inf=% .6e F24_minus_F3_l2=% .6e (tol %.0e)\n",
        (double)f3_2, (double)f24_2, (double)f24_inf, (double)f24m3_2, (double)compat_tol));

    /* -------- 4.4 projection compatibility -------- */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n===== 4.4 projection compatibility =====\n"));
    PetscCall(FpRestore(user, &entry));
    PetscCall(VecCopy(Xfinal, user->Ucont));
    PetscCall(SynchronizePeriodicStaggeredFields(user, 1, stag));
    PetscCall(Contra2Cart(user));
    PetscCall(SynchronizePeriodicCellFields(user, 1, cell));
    PetscCall(UpdateLocalGhosts(user, "Ucat"));
    PetscCall(ApplyBoundaryConditions(user));
    PetscCall(PoissonSolver_MG(&simCtx->usermg));
    PetscCall(UpdatePressure(user));
    PetscCall(Projection(user));
    PetscCall(UpdateLocalGhosts(user, "P"));
    PetscCall(ComputeDivergence(user));
    max_div = simCtx->MaxDiv;
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "PROJECTION max_divergence=% .6e net_divergence_sum=% .6e (tol %.0e)\n",
        (double)max_div, (double)simCtx->summationRHS, (double)div_tol));

    /* -------- consolidated verdicts (printed before assertions) -------- */
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\nREGRESSION_VERDICTS (ranks=%d)\n", (int)world));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  PRODUCTION_PURITY: %s\n",
        (imm_inf <= purity_inf_tol && imm_2 <= purity_2_tol &&
         mffd_inf <= purity_inf_tol && mffd_2 <= purity_2_tol) ? "PASS" : "FAIL"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  FULL_STEP2_DEFAULT_GS: %s\n",
        res_cgs.converged ? "PASS" : "FAIL"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  PRODUCTION_SOLVER_WRAPPER: %s\n",
        (prod_ierr == 0 && simCtx->mom_last_converged) ? "PASS" : "FAIL"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  THREE_PASS_24_PASS_COMPATIBLE: %s\n",
        (f24_2 <= compat_tol) ? "PASS" : "FAIL"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  PROJECTION_DIVERGENCE_FREE: %s\n",
        (max_div <= div_tol) ? "PASS" : "FAIL"));

    /* -------- assertions -------- */
    PetscCheck(imm_inf <= purity_inf_tol && imm_2 <= purity_2_tol, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
        "production residual not immediately repeatable (inf=%g l2=%g)", (double)imm_inf, (double)imm_2);
    PetscCheck(mffd_inf <= purity_inf_tol && mffd_2 <= purity_2_tol, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
        "production residual changed after MFFD products (inf=%g l2=%g)", (double)mffd_inf, (double)mffd_2);
    PetscCheck(res_cgs.converged, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
        "full step-2 solve (default GS) did not converge (reason=%d)", (int)res_cgs.reason);
    PetscCheck(res_cgs.reason != SNES_DIVERGED_LINE_SEARCH && res_cgs.reason != SNES_DIVERGED_LINEAR_SOLVE,
        PETSC_COMM_WORLD, PETSC_ERR_PLIB, "full step-2 solve diverged via line search or linear solve (reason=%d)",
        (int)res_cgs.reason);
    PetscCheck(prod_ierr == 0 && simCtx->mom_last_converged, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
        "production Newton-Krylov solver wrapper did not converge (ierr=%d)", (int)prod_ierr);
    PetscCheck(f24_2 <= compat_tol, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
        "converged three-pass solution does not satisfy the 24-pass residual (||F24||2=%g)", (double)f24_2);
    PetscCheck(max_div <= div_tol, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
        "projection of the converged momentum state is not divergence free (max_div=%g)", (double)max_div);

    PetscCall(VecDestroy(&F24)); PetscCall(VecDestroy(&F3));
    PetscCall(VecDestroy(&Xfinal)); PetscCall(VecDestroy(&res_cgs.Xfinal));
    PetscCall(VecDestroy(&jv)); PetscCall(VecDestroy(&v));
    PetscCall(VecDestroy(&Fafter)); PetscCall(VecDestroy(&F1)); PetscCall(VecDestroy(&F0));
    PetscCall(VecDestroy(&X1)); PetscCall(VecDestroy(&X0));
    PetscCall(FpDestroy(&entry));
    PetscCall(VecDestroy(&user->Rhs));
    PetscCall(PicurvDestroyRuntimeContext(&simCtx));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscFunctionReturn(PETSC_SUCCESS);
}

/** @brief Runs the opt-in Newton--Krylov production regression executable. */
int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"newton-krylov-production-regression", TestNewtonKrylovProductionRegression}
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv Newton-Krylov production regression");
    if (ierr) return (int)ierr;
    ierr = PicurvRunTests("unit-momentum-newton-boundary-fixedpoint", cases, 1);
    if (PetscFinalize()) return 1;
    return (int)ierr;
}
