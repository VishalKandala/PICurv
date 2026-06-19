/**
 * @file test_momentum_convective_candidates.c
 * @brief A4a focused convective-candidate study (states A-C only).
 *
 * Builds the finite-difference Jacobian of the ACTUAL production convection-only residual
 * (ComputeRHS with inviscid, P=0, centered, periodic, single block) on a tiny periodic
 * Cartesian grid, then compares the B/C/D estimator candidates against rho(J), sigma_max(J),
 * the exact 4-stage RK matrix polynomial P(z)=1+z+z^2/2+z^3/6+z^4/24, and a direct anchored
 * 4-stage perturbation cross-check. Shadow-only: changes no production default.
 *
 * States: A uniform divergence-free; B nonzero discrete divergence; C divergence-free shear.
 */

#include "test_support.h"
#include "momentumsolvers.h"
#include "rhs.h"
#include "setup.h"
#include "Boundaries.h"
#include <petscblaslapack.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

/* ----------------------------------------------------------------------------------- *
 *  Periodic independent staggered-DOF map.                                             *
 * ----------------------------------------------------------------------------------- */
typedef struct { PetscInt n, expected_n; PetscInt *comp, *ci, *cj, *ck; } DofMap;
typedef struct {
    PetscReal declared[3];
    PetscReal ucat_global[3];
    PetscReal ucat_ghost[3];
    PetscReal ucont_global[3];
} SeamDiagnostics;
typedef struct { PetscReal n2, ninf, checksum; } GlobalVecStats;

static char g_ref_path[PETSC_MAX_PATH_LEN] = "";
static char g_ref_token[128] = "";
static PetscBool g_ref_path_set = PETSC_FALSE;
static PetscBool g_ref_token_set = PETSC_FALSE;

typedef enum {
    STABLE_CFL_NONE,
    STABLE_CFL_FINITE,
    STABLE_CFL_EXCEEDS_SCAN
} StableCFLStatus;

typedef struct {
    StableCFLStatus status;
    PetscReal cfl;
} StableCFLResult;

/**
 * @brief Returns the number of independent periodic representatives in one direction.
 */
static inline PetscInt PeriodicRepCount(PetscInt npts) { return npts - 2; }

/**
 * @brief Counts all independent component-staggered representatives used by ComputeRHS.
 */
static inline PetscInt DofMapExpectedCount(DMDALocalInfo info)
{
    return 3 * PeriodicRepCount(info.mx) * PeriodicRepCount(info.my) * PeriodicRepCount(info.mz);
}

/**
 * @brief Builds the serial periodic independent face-DOF map used by dense Jacobians.
 *
 * Production periodic synchronization copies global plane 0 from mx-2 and plane mx-1
 * from 1 (and analogously in y/z), so representatives 1..m-2 contain each independent
 * component-staggered Ucont face DOF exactly once. Perturbations only touch these reps;
 * EvalConvResidual() then calls SynchronizePeriodicStaggeredFields() to update duplicates.
 */
static PetscErrorCode DofMapBuild(UserCtx *user, DofMap *map)
{
    DMDALocalInfo info = user->info;
    const PetscInt mx = info.mx, my = info.my, mz = info.mz;
    const PetscInt lxs = 1, lxe = mx-1, lys = 1, lye = my-1, lzs = 1, lze = mz-1;
    PetscInt cnt = 0;
    PetscFunctionBeginUser;
    map->expected_n = DofMapExpectedCount(info);
    map->n = (lxe-lxs)*(lye-lys)*(lze-lzs)*3;
    PetscCall(PetscMalloc4(map->n, &map->comp, map->n, &map->ci, map->n, &map->cj, map->n, &map->ck));
    for (PetscInt k = lzs; k < lze; k++)
        for (PetscInt j = lys; j < lye; j++)
            for (PetscInt i = lxs; i < lxe; i++)
                for (PetscInt c = 0; c < 3; c++) {
                    map->comp[cnt] = c; map->ci[cnt] = i; map->cj[cnt] = j; map->ck[cnt] = k; cnt++;
                }
    PetscCheck(map->n == map->expected_n, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ,
               "periodic independent DOF count mismatch: got %" PetscInt_FMT ", expected %" PetscInt_FMT,
               map->n, map->expected_n);
    PetscFunctionReturn(0);
}

/**
 * @brief Builds the rank-owned periodic independent face-DOF map for MPI checks.
 */
static PetscErrorCode DofMapBuildOwned(UserCtx *user, DofMap *map)
{
    DMDALocalInfo info = user->info;
    const PetscInt mx = info.mx, my = info.my, mz = info.mz;
    const PetscInt xs = info.xs, xe = info.xs + info.xm;
    const PetscInt ys = info.ys, ye = info.ys + info.ym;
    const PetscInt zs = info.zs, ze = info.zs + info.zm;
    const PetscInt lxs = PetscMax(xs, 1), lxe = PetscMin(xe, mx-1);
    const PetscInt lys = PetscMax(ys, 1), lye = PetscMin(ye, my-1);
    const PetscInt lzs = PetscMax(zs, 1), lze = PetscMin(ze, mz-1);
    PetscInt cnt = 0;
    PetscFunctionBeginUser;
    map->expected_n = DofMapExpectedCount(info);
    map->n = PetscMax(0,lxe-lxs)*PetscMax(0,lye-lys)*PetscMax(0,lze-lzs)*3;
    PetscCall(PetscMalloc4(map->n, &map->comp, map->n, &map->ci, map->n, &map->cj, map->n, &map->ck));
    for (PetscInt k = lzs; k < lze; k++)
        for (PetscInt j = lys; j < lye; j++)
            for (PetscInt i = lxs; i < lxe; i++)
                for (PetscInt c = 0; c < 3; c++) {
                    map->comp[cnt] = c; map->ci[cnt] = i; map->cj[cnt] = j; map->ck[cnt] = k; cnt++;
                }
    PetscFunctionReturn(0);
}

/**
 * @brief Releases storage owned by an active-DOF map.
 */
static PetscErrorCode DofMapDestroy(DofMap *map)
{ PetscFunctionBeginUser; PetscCall(PetscFree4(map->comp, map->ci, map->cj, map->ck)); PetscFunctionReturn(0); }

/* component accessor: Cmpnts is 3 contiguous PetscReal (x,y,z) in the real build. */
static inline PetscReal CmpGet(Cmpnts c, PetscInt comp) { const PetscReal *p = (const PetscReal*)&c; return p[comp]; }

/* ----------------------------------------------------------------------------------- *
 *  Deterministic convection-only residual wrapper using the real production path.      *
 * ----------------------------------------------------------------------------------- */
static PetscErrorCode EvalConvResidual(UserCtx *user, Vec Ucont_in, Vec Rhs, const DofMap *map, PetscReal *Ract)
{
    Cmpnts ***r;
    PetscFunctionBeginUser;
    PetscCall(VecCopy(Ucont_in, user->Ucont));
    {
        const char *fld[] = {"Ucont"};
        PetscCall(SynchronizePeriodicStaggeredFields(user, 1, fld));   /* global->local + periodic */
    }
    PetscCall(ComputeRHS(user, Rhs));                                  /* Contra2Cart + Convection + mapping */
    PetscCall(DMDAVecGetArrayRead(user->fda, Rhs, &r));
    for (PetscInt m = 0; m < map->n; m++)
        Ract[m] = CmpGet(r[map->ck[m]][map->cj[m]][map->ci[m]], map->comp[m]);
    PetscCall(DMDAVecRestoreArrayRead(user->fda, Rhs, &r));
    PetscFunctionReturn(0);
}

/* perturb one active global Ucont DOF (then caller re-syncs through EvalConvResidual). */
static PetscErrorCode PerturbDof(UserCtx *user, Vec Ucont, const DofMap *map, PetscInt m, PetscReal delta)
{
    Cmpnts ***a;
    PetscFunctionBeginUser;
    PetscCall(DMDAVecGetArray(user->fda, Ucont, &a));
    { PetscReal *p = (PetscReal*)&a[map->ck[m]][map->cj[m]][map->ci[m]]; p[map->comp[m]] += delta; }
    PetscCall(DMDAVecRestoreArray(user->fda, Ucont, &a));
    PetscFunctionReturn(0);
}

/**
 * @brief Reads one active contravariant component from a global vector.
 */
static PetscErrorCode GetDof(UserCtx *user, Vec Ucont, const DofMap *map, PetscInt m, PetscReal *val)
{
    Cmpnts ***a;
    PetscFunctionBeginUser;
    PetscCall(DMDAVecGetArray(user->fda, Ucont, &a));
    { const PetscReal *p = (const PetscReal*)&a[map->ck[m]][map->cj[m]][map->ci[m]]; *val = p[map->comp[m]]; }
    PetscCall(DMDAVecRestoreArray(user->fda, Ucont, &a));
    PetscFunctionReturn(0);
}

/* ----------------------------------------------------------------------------------- *
 *  Base-state construction: physical Cartesian velocity -> production contravariant.   *
 * ----------------------------------------------------------------------------------- */
typedef enum { STATE_A, STATE_B, STATE_C } CandState;

/**
 * @brief Returns the cell-centered periodic angle using duplicated endpoint planes.
 */
static inline PetscReal PeriodicCellAngle(PetscInt idx, PetscInt npts)
{
    const PetscInt nuniq = npts - 1;
    const PetscInt ip = (idx == npts - 1) ? 0 : idx;
    return 2.0*PETSC_PI*((PetscReal)ip)/((PetscReal)nuniq);
}

/**
 * @brief Returns a face-representative periodic angle for component-staggered Ucont.
 */
static inline PetscReal PeriodicFaceAngle(PetscInt idx, PetscInt npts)
{
    const PetscInt nuniq = PeriodicRepCount(npts);
    PetscInt ip = idx - 1;
    if (idx == 0) ip = nuniq - 1;
    else if (idx == npts - 1) ip = 0;
    return 2.0*PETSC_PI*((PetscReal)ip)/((PetscReal)nuniq);
}

/**
 * @brief Evaluates one of the three analytic Cartesian candidate states.
 */
static inline Cmpnts AnalyticVelocity(CandState st, PetscInt i, PetscInt j, PetscInt k,
                                      PetscInt mx, PetscInt my, PetscInt mz)
{
    Cmpnts v;
    const PetscReal x = PeriodicCellAngle(i, mx);
    const PetscReal y = PeriodicCellAngle(j, my);
    (void)k; (void)mz;
    if (st == STATE_A) {
        v.x = 0.7; v.y = -0.4; v.z = 0.0;
    } else if (st == STATE_B) {
        (void)x; v.x = 0.0; v.y = 0.0; v.z = 0.0;
    } else {
        v.x = 1.0 + 0.5*PetscSinReal(y); v.y = 0.0; v.z = 0.0;
    }
    return v;
}

/**
 * @brief Evaluates the declared direct component-staggered State B Ucont field.
 */
static inline Cmpnts DirectUcontVelocity(CandState st, PetscInt i, PetscInt j, PetscInt k,
                                         PetscInt mx, PetscInt my, PetscInt mz)
{
    Cmpnts v = {0.0, 0.0, 0.0};
    (void)j; (void)k; (void)my; (void)mz;
    if (st == STATE_B) v.x = PetscSinReal(PeriodicFaceAngle(i, mx));
    return v;
}

/**
 * @brief Computes the componentwise infinity norm of the difference between two vectors.
 */
static inline PetscReal CmpDiffInf(Cmpnts a, Cmpnts b)
{
    return PetscMax(PetscAbsReal(a.x-b.x), PetscMax(PetscAbsReal(a.y-b.y), PetscAbsReal(a.z-b.z)));
}

/**
 * @brief Computes analytic periodic seam mismatches for each coordinate direction.
 */
static PetscErrorCode ComputeDeclaredSeamMismatch(CandState st, DMDALocalInfo info, PetscReal seam[3])
{
    PetscReal loc[3] = {0.0, 0.0, 0.0}, glo[3];
    const PetscInt mx = info.mx, my = info.my, mz = info.mz;
    PetscFunctionBeginUser;
    if (st == STATE_B) {
        for (PetscInt k = 1; k < mz-1; k++)
            for (PetscInt j = 1; j < my-1; j++) {
                loc[0] = PetscMax(loc[0], CmpDiffInf(DirectUcontVelocity(st, 0, j, k, mx, my, mz),
                                                     DirectUcontVelocity(st, mx-2, j, k, mx, my, mz)));
                loc[0] = PetscMax(loc[0], CmpDiffInf(DirectUcontVelocity(st, mx-1, j, k, mx, my, mz),
                                                     DirectUcontVelocity(st, 1, j, k, mx, my, mz)));
            }
        for (PetscInt k = 1; k < mz-1; k++)
            for (PetscInt i = 1; i < mx-1; i++) {
                loc[1] = PetscMax(loc[1], CmpDiffInf(DirectUcontVelocity(st, i, 0, k, mx, my, mz),
                                                     DirectUcontVelocity(st, i, my-2, k, mx, my, mz)));
                loc[1] = PetscMax(loc[1], CmpDiffInf(DirectUcontVelocity(st, i, my-1, k, mx, my, mz),
                                                     DirectUcontVelocity(st, i, 1, k, mx, my, mz)));
            }
        for (PetscInt j = 1; j < my-1; j++)
            for (PetscInt i = 1; i < mx-1; i++) {
                loc[2] = PetscMax(loc[2], CmpDiffInf(DirectUcontVelocity(st, i, j, 0, mx, my, mz),
                                                     DirectUcontVelocity(st, i, j, mz-2, mx, my, mz)));
                loc[2] = PetscMax(loc[2], CmpDiffInf(DirectUcontVelocity(st, i, j, mz-1, mx, my, mz),
                                                     DirectUcontVelocity(st, i, j, 1, mx, my, mz)));
            }
    } else {
        for (PetscInt k = 0; k < mz; k++)
            for (PetscInt j = 0; j < my; j++)
                loc[0] = PetscMax(loc[0], CmpDiffInf(AnalyticVelocity(st, 0, j, k, mx, my, mz),
                                                     AnalyticVelocity(st, mx-1, j, k, mx, my, mz)));
        for (PetscInt k = 0; k < mz; k++)
            for (PetscInt i = 0; i < mx; i++)
                loc[1] = PetscMax(loc[1], CmpDiffInf(AnalyticVelocity(st, i, 0, k, mx, my, mz),
                                                     AnalyticVelocity(st, i, my-1, k, mx, my, mz)));
        for (PetscInt j = 0; j < my; j++)
            for (PetscInt i = 0; i < mx; i++)
                loc[2] = PetscMax(loc[2], CmpDiffInf(AnalyticVelocity(st, i, j, 0, mx, my, mz),
                                                     AnalyticVelocity(st, i, j, mz-1, mx, my, mz)));
    }
    PetscCallMPI(MPI_Allreduce(loc, glo, 3, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD));
    seam[0] = glo[0]; seam[1] = glo[1]; seam[2] = glo[2];
    PetscFunctionReturn(0);
}

/**
 * @brief Reports whether a global index is present in a rank's local ghosted range.
 */
static inline PetscBool InGhostRange(PetscInt idx, PetscInt lo, PetscInt n)
{ return (PetscBool)(idx >= lo && idx < lo + n); }

/**
 * @brief Computes duplicate-plane mismatch in a local vector view.
 */
static PetscErrorCode ComputeLocalDuplicateMismatch(UserCtx *user, Vec local, PetscReal seam[3])
{
    DMDALocalInfo info = user->info;
    Cmpnts ***a;
    PetscReal loc[3] = {0.0, 0.0, 0.0}, glo[3];
    const PetscInt mx = info.mx, my = info.my, mz = info.mz;
    PetscFunctionBeginUser;
    PetscCall(DMDAVecGetArrayRead(user->fda, local, &a));
#define HAVE_I(ii) InGhostRange((ii), info.gxs, info.gxm)
#define HAVE_J(jj) InGhostRange((jj), info.gys, info.gym)
#define HAVE_K(kk) InGhostRange((kk), info.gzs, info.gzm)
    if (HAVE_I(0) && HAVE_I(mx-2))
        for (PetscInt k = 1; k < mz-1; k++) if (HAVE_K(k))
            for (PetscInt j = 1; j < my-1; j++) if (HAVE_J(j))
                loc[0] = PetscMax(loc[0], CmpDiffInf(a[k][j][0], a[k][j][mx-2]));
    if (HAVE_I(mx-1) && HAVE_I(1))
        for (PetscInt k = 1; k < mz-1; k++) if (HAVE_K(k))
            for (PetscInt j = 1; j < my-1; j++) if (HAVE_J(j))
                loc[0] = PetscMax(loc[0], CmpDiffInf(a[k][j][mx-1], a[k][j][1]));
    if (HAVE_J(0) && HAVE_J(my-2))
        for (PetscInt k = 1; k < mz-1; k++) if (HAVE_K(k))
            for (PetscInt i = 1; i < mx-1; i++) if (HAVE_I(i))
                loc[1] = PetscMax(loc[1], CmpDiffInf(a[k][0][i], a[k][my-2][i]));
    if (HAVE_J(my-1) && HAVE_J(1))
        for (PetscInt k = 1; k < mz-1; k++) if (HAVE_K(k))
            for (PetscInt i = 1; i < mx-1; i++) if (HAVE_I(i))
                loc[1] = PetscMax(loc[1], CmpDiffInf(a[k][my-1][i], a[k][1][i]));
    if (HAVE_K(0) && HAVE_K(mz-2))
        for (PetscInt j = 1; j < my-1; j++) if (HAVE_J(j))
            for (PetscInt i = 1; i < mx-1; i++) if (HAVE_I(i))
                loc[2] = PetscMax(loc[2], CmpDiffInf(a[0][j][i], a[mz-2][j][i]));
    if (HAVE_K(mz-1) && HAVE_K(1))
        for (PetscInt j = 1; j < my-1; j++) if (HAVE_J(j))
            for (PetscInt i = 1; i < mx-1; i++) if (HAVE_I(i))
                loc[2] = PetscMax(loc[2], CmpDiffInf(a[mz-1][j][i], a[1][j][i]));
#undef HAVE_I
#undef HAVE_J
#undef HAVE_K
    PetscCall(DMDAVecRestoreArrayRead(user->fda, local, &a));
    PetscCallMPI(MPI_Allreduce(loc, glo, 3, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD));
    seam[0] = glo[0]; seam[1] = glo[1]; seam[2] = glo[2];
    PetscFunctionReturn(0);
}

/**
 * @brief Computes outer periodic ghost mismatch for local Ucat.
 */
static PetscErrorCode ComputeLocalOuterGhostMismatch(UserCtx *user, Vec local, PetscReal seam[3])
{
    DMDALocalInfo info = user->info;
    Cmpnts ***a;
    PetscReal loc[3] = {0.0, 0.0, 0.0}, glo[3];
    const PetscInt mx = info.mx, my = info.my, mz = info.mz;
    PetscFunctionBeginUser;
    PetscCall(DMDAVecGetArrayRead(user->fda, local, &a));
#define HAVE_I(ii) InGhostRange((ii), info.gxs, info.gxm)
#define HAVE_J(jj) InGhostRange((jj), info.gys, info.gym)
#define HAVE_K(kk) InGhostRange((kk), info.gzs, info.gzm)
    if (HAVE_I(-1) && HAVE_I(1))
        for (PetscInt k = 1; k < mz-1; k++) if (HAVE_K(k))
            for (PetscInt j = 1; j < my-1; j++) if (HAVE_J(j))
                loc[0] = PetscMax(loc[0], CmpDiffInf(a[k][j][-1], a[k][j][1]));
    if (HAVE_I(mx) && HAVE_I(mx-2))
        for (PetscInt k = 1; k < mz-1; k++) if (HAVE_K(k))
            for (PetscInt j = 1; j < my-1; j++) if (HAVE_J(j))
                loc[0] = PetscMax(loc[0], CmpDiffInf(a[k][j][mx], a[k][j][mx-2]));
    if (HAVE_J(-1) && HAVE_J(1))
        for (PetscInt k = 1; k < mz-1; k++) if (HAVE_K(k))
            for (PetscInt i = 1; i < mx-1; i++) if (HAVE_I(i))
                loc[1] = PetscMax(loc[1], CmpDiffInf(a[k][-1][i], a[k][1][i]));
    if (HAVE_J(my) && HAVE_J(my-2))
        for (PetscInt k = 1; k < mz-1; k++) if (HAVE_K(k))
            for (PetscInt i = 1; i < mx-1; i++) if (HAVE_I(i))
                loc[1] = PetscMax(loc[1], CmpDiffInf(a[k][my][i], a[k][my-2][i]));
    if (HAVE_K(-1) && HAVE_K(1))
        for (PetscInt j = 1; j < my-1; j++) if (HAVE_J(j))
            for (PetscInt i = 1; i < mx-1; i++) if (HAVE_I(i))
                loc[2] = PetscMax(loc[2], CmpDiffInf(a[-1][j][i], a[1][j][i]));
    if (HAVE_K(mz) && HAVE_K(mz-2))
        for (PetscInt j = 1; j < my-1; j++) if (HAVE_J(j))
            for (PetscInt i = 1; i < mx-1; i++) if (HAVE_I(i))
                loc[2] = PetscMax(loc[2], CmpDiffInf(a[mz][j][i], a[mz-2][j][i]));
#undef HAVE_I
#undef HAVE_J
#undef HAVE_K
    PetscCall(DMDAVecRestoreArrayRead(user->fda, local, &a));
    PetscCallMPI(MPI_Allreduce(loc, glo, 3, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD));
    seam[0] = glo[0]; seam[1] = glo[1]; seam[2] = glo[2];
    PetscFunctionReturn(0);
}

/**
 * @brief Configures the minimal context for centered inviscid periodic convection tests.
 */
static PetscErrorCode ConfigureCandidateFixture(SimCtx *simCtx, UserCtx *user)
{
    PetscFunctionBeginUser;
    for (int f = 0; f < 6; f++) user->boundary_faces[f].mathematical_type = PERIODIC;
    simCtx->dt = 0.1; simCtx->step = 5; simCtx->StartStep = 0;     /* BDF2 -> a0=1.5 */
    simCtx->ren = 1.0e6; simCtx->invicid = 1; simCtx->les = 0; simCtx->rans = 0;
    simCtx->central = 1; simCtx->clark = 0; simCtx->TwoD = 0; simCtx->block_number = 1;
    simCtx->bulkVelocityCorrection = 0.0; simCtx->moveframe = 0; simCtx->rotateframe = 0;
    if (!user->lNu_t) PetscCall(DMCreateLocalVector(user->da, &user->lNu_t));
    PetscCall(VecSet(user->P, 0.0)); PetscCall(UpdateLocalGhosts(user, "P"));
    PetscCall(VecSet(user->lNvert, 0.0)); PetscCall(VecSet(user->Nvert, 0.0));
    PetscFunctionReturn(0);
}

/* set global Ucat to the analytic velocity for the state at cell-centre (i,j,k). */
static PetscErrorCode SetUcatField(UserCtx *user, CandState st)
{
    Cmpnts ***u;
    DMDALocalInfo info = user->info;
    const PetscInt mx = info.mx, my = info.my, mz = info.mz;
    PetscFunctionBeginUser;
    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &u));
    /* only owned region for a GLOBAL vec */
    for (PetscInt k = info.zs; k < info.zs+info.zm; k++)
        for (PetscInt j = info.ys; j < info.ys+info.ym; j++)
            for (PetscInt i = info.xs; i < info.xs+info.xm; i++) {
                u[k][j][i] = AnalyticVelocity(st, i, j, k, mx, my, mz);
            }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &u));
    PetscFunctionReturn(0);
}

/**
 * @brief Sets the direct State B component-staggered Ucont field on owned entries.
 */
static PetscErrorCode SetDirectUcontField(UserCtx *user, CandState st)
{
    Cmpnts ***u;
    DMDALocalInfo info = user->info;
    const PetscInt mx = info.mx, my = info.my, mz = info.mz;
    PetscFunctionBeginUser;
    PetscCall(VecSet(user->Ucont, 0.0));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucont, &u));
    for (PetscInt k = info.zs; k < info.zs+info.zm; k++)
        for (PetscInt j = info.ys; j < info.ys+info.ym; j++)
            for (PetscInt i = info.xs; i < info.xs+info.xm; i++)
                u[k][j][i] = DirectUcontVelocity(st, i, j, k, mx, my, mz);
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucont, &u));
    PetscFunctionReturn(0);
}

/**
 * @brief Computes max discrete divergence of the current local Ucont field.
 */
static PetscErrorCode ComputeMaxDiscreteDivergence(UserCtx *user, PetscReal *maxdiv)
{
    DMDALocalInfo info = user->info;
    const PetscInt mx = info.mx, my = info.my, mz = info.mz;
    Cmpnts ***uc;
    PetscReal dv = 0.0, dv_global;
    PetscFunctionBeginUser;
    PetscCall(UpdateLocalGhosts(user, "Ucont"));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lUcont, &uc));
    for (PetscInt k = info.zs; k < info.zs+info.zm; k++)
        for (PetscInt j = info.ys; j < info.ys+info.ym; j++)
            for (PetscInt i = info.xs; i < info.xs+info.xm; i++) {
                if (i<1||i>mx-2||j<1||j>my-2||k<1||k>mz-2) continue;
                const PetscReal d = (uc[k][j][i].x - uc[k][j][i-1].x)
                                  + (uc[k][j][i].y - uc[k][j-1][i].y)
                                  + (uc[k][j][i].z - uc[k-1][j][i].z);
                dv = PetscMax(dv, PetscAbsReal(d));
            }
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lUcont, &uc));
    PetscCallMPI(MPI_Allreduce(&dv, &dv_global, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD));
    *maxdiv = dv_global;
    PetscFunctionReturn(0);
}

/* Build the base state. States A/C are declared in Cartesian space and converted through
   Cart2Contra; State B is declared directly in synchronized component-staggered Ucont space. */
static PetscErrorCode BuildBaseState(UserCtx *user, CandState st, Vec Ubase,
                                     PetscReal *repeat_inf, PetscReal *maxdiv,
                                     SeamDiagnostics *seam)
{
    DMDALocalInfo info = user->info;
    const PetscInt mx = info.mx, my = info.my, mz = info.mz;
    PetscReal err = 0.0, err_global;
    Vec target;
    Cmpnts ***ur, ***ut;
    PetscFunctionBeginUser;

    PetscCall(VecDuplicate(user->Ucat, &target));

    if (seam) PetscCall(ComputeDeclaredSeamMismatch(st, info, seam->declared));

    if (st == STATE_B) {
        const char *ufld[] = {"Ucont"};
        const char *cfld[] = {"Ucat"};
        PetscCall(SetDirectUcontField(user, st));
        PetscCall(SynchronizePeriodicStaggeredFields(user, 1, ufld));
        PetscCall(VecCopy(user->Ucont, Ubase));
        if (seam) PetscCall(ComputeLocalDuplicateMismatch(user, user->lUcont, seam->ucont_global));
        PetscCall(Contra2Cart(user));
        PetscCall(SynchronizePeriodicCellFields(user, 1, cfld));
        PetscCall(VecCopy(user->Ucat, target));          /* saved recovered Cartesian state */
        PetscCall(Contra2Cart(user));
        PetscCall(SynchronizePeriodicCellFields(user, 1, cfld));
    } else {
        const char *cfld[] = {"Ucat"};
        const char *ufld[] = {"Ucont"};
        PetscCall(SetUcatField(user, st));
        PetscCall(SynchronizePeriodicCellFields(user, 1, cfld));
        PetscCall(VecCopy(user->Ucat, target));          /* saved declared Cartesian state */
        PetscCall(Cart2Contra(user));                    /* global Ucont from lUcat + metrics */
        PetscCall(SynchronizePeriodicStaggeredFields(user, 1, ufld));
        PetscCall(VecCopy(user->Ucont, Ubase));
        if (seam) PetscCall(ComputeLocalDuplicateMismatch(user, user->lUcont, seam->ucont_global));
        PetscCall(Contra2Cart(user));                    /* recovered Cartesian from Ucont */
        PetscCall(SynchronizePeriodicCellFields(user, 1, cfld));
    }

    if (seam) {
        PetscCall(UpdateLocalGhosts(user, "Ucat"));
        PetscCall(ComputeLocalDuplicateMismatch(user, user->lUcat, seam->ucat_global));
        PetscCall(ComputeLocalOuterGhostMismatch(user, user->lUcat, seam->ucat_ghost));
    }

    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucat, &ur));
    PetscCall(DMDAVecGetArrayRead(user->fda, target, &ut));
    for (PetscInt k = info.zs; k < info.zs+info.zm; k++)
        for (PetscInt j = info.ys; j < info.ys+info.ym; j++)
            for (PetscInt i = info.xs; i < info.xs+info.xm; i++) {
                if (i<1||i>mx-2||j<1||j>my-2||k<1||k>mz-2) continue;   /* interior cells only */
                err = PetscMax(err, PetscAbsReal(ur[k][j][i].x - ut[k][j][i].x));
                err = PetscMax(err, PetscAbsReal(ur[k][j][i].y - ut[k][j][i].y));
                err = PetscMax(err, PetscAbsReal(ur[k][j][i].z - ut[k][j][i].z));
            }
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucat, &ur));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, target, &ut));
    PetscCall(VecDestroy(&target));
    PetscCall(UpdateLocalGhosts(user, "Ucat"));         /* lUcat now consistent with base Ucont */

    PetscCallMPI(MPI_Allreduce(&err, &err_global, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD));
    PetscCall(ComputeMaxDiscreteDivergence(user, maxdiv));
    *repeat_inf = err_global;
    PetscFunctionReturn(0);
}

/**
 * @brief Computes the global maximum Cartesian velocity-gradient row-sum used by Candidate D.
 */
static PetscErrorCode ComputeMaxGradientContribution(UserCtx *user, PetscReal *gradmax)
{
    DMDALocalInfo info = user->info;
    Cmpnts ***ucat, ***csi, ***eta, ***zet;
    PetscReal ***aj;
    PetscReal loc = 0.0, glo;
    const PetscInt mx = info.mx, my = info.my, mz = info.mz;
    PetscFunctionBeginUser;
    PetscCall(UpdateLocalGhosts(user, "Ucat"));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lUcat, &ucat));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCsi,  &csi));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lEta,  &eta));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lZet,  &zet));
    PetscCall(DMDAVecGetArrayRead(user->da,  user->lAj,   &aj));
    for (PetscInt k = info.zs; k < info.zs+info.zm; k++)
        for (PetscInt j = info.ys; j < info.ys+info.ym; j++)
            for (PetscInt i = info.xs; i < info.xs+info.xm; i++) {
                if (i<1||i>mx-2||j<1||j>my-2||k<1||k>mz-2) continue;
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
                const PetscReal Ajc = aj[k][j][i];
#define ROWSUM(cmp) ( \
    PetscAbsReal(Ajc*(C.x*duc.cmp + E.x*due.cmp + Z.x*duz.cmp)) + \
    PetscAbsReal(Ajc*(C.y*duc.cmp + E.y*due.cmp + Z.y*duz.cmp)) + \
    PetscAbsReal(Ajc*(C.z*duc.cmp + E.z*due.cmp + Z.z*duz.cmp)) )
                loc = PetscMax(loc, PetscMax(ROWSUM(x), PetscMax(ROWSUM(y), ROWSUM(z))));
#undef ROWSUM
            }
    PetscCall(DMDAVecRestoreArrayRead(user->da,  user->lAj,   &aj));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lZet,  &zet));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lEta,  &eta));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCsi,  &csi));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lUcat, &ucat));
    PetscCallMPI(MPI_Allreduce(&loc, &glo, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD));
    *gradmax = glo;
    PetscFunctionReturn(0);
}

/* ----------------------------------------------------------------------------------- *
 *  Dense linear algebra on the (column-major) active Jacobian via PETSc LAPACK.        *
 * ----------------------------------------------------------------------------------- */
/* rho(J): max |eigenvalue|. A is column-major n x n and is COPIED (dgeev overwrites). */
static PetscErrorCode DenseSpectralRadius(const PetscReal *A, PetscInt n, PetscReal *rho,
                                          PetscReal *maxRealPart)
{
    PetscBLASInt N, lda, lwork, info;
    PetscReal *Acopy, *wr, *wi, *work, dummy = 0.0;
    PetscFunctionBeginUser;
    PetscCall(PetscBLASIntCast(n, &N)); lda = N; lwork = 8*N;
    PetscCall(PetscMalloc4(n*n, &Acopy, n, &wr, n, &wi, (size_t)lwork, &work));
    for (PetscInt t = 0; t < n*n; t++) Acopy[t] = A[t];
    {
        char nochar = 'N';
        LAPACKgeev_(&nochar, &nochar, &N, Acopy, &lda, wr, wi, &dummy, &lda, &dummy, &lda, work, &lwork, &info);
    }
    PetscCheck(info == 0, PETSC_COMM_SELF, PETSC_ERR_LIB, "LAPACK dgeev failed: info=%d", (int)info);
    *rho = 0.0; *maxRealPart = -PETSC_MAX_REAL;
    for (PetscInt t = 0; t < n; t++) {
        *rho = PetscMax(*rho, PetscSqrtReal(wr[t]*wr[t] + wi[t]*wi[t]));
        *maxRealPart = PetscMax(*maxRealPart, wr[t]);
    }
    PetscCall(PetscFree4(Acopy, wr, wi, work));
    PetscFunctionReturn(0);
}

/**
 * @brief Computes the spectral radius of the RK polynomial by applying it to eig(J).
 */
static PetscErrorCode DenseRKPolynomialSpectralRadius(const PetscReal *J, PetscInt n,
                                                      PetscReal dtau, PetscReal *rho)
{
    PetscBLASInt N, lda, lwork, info;
    PetscReal *Jcopy, *wr, *wi, *work, dummy = 0.0;
    PetscFunctionBeginUser;
    PetscCall(PetscBLASIntCast(n, &N)); lda = N; lwork = 8*N;
    PetscCall(PetscMalloc4(n*n, &Jcopy, n, &wr, n, &wi, (size_t)lwork, &work));
    for (PetscInt t = 0; t < n*n; t++) Jcopy[t] = J[t];
    {
        char nochar = 'N';
        LAPACKgeev_(&nochar, &nochar, &N, Jcopy, &lda, wr, wi, &dummy, &lda, &dummy, &lda, work, &lwork, &info);
    }
    PetscCheck(info == 0, PETSC_COMM_SELF, PETSC_ERR_LIB, "LAPACK dgeev failed: info=%d", (int)info);
    *rho = 0.0;
    for (PetscInt t = 0; t < n; t++) {
        const PetscReal zr = dtau*wr[t], zi = dtau*wi[t];
        const PetscReal z2r = zr*zr - zi*zi, z2i = 2.0*zr*zi;
        const PetscReal z3r = z2r*zr - z2i*zi, z3i = z2r*zi + z2i*zr;
        const PetscReal z4r = z3r*zr - z3i*zi, z4i = z3r*zi + z3i*zr;
        const PetscReal pr = 1.0 + zr + 0.5*z2r + z3r/6.0 + z4r/24.0;
        const PetscReal pi =       zi + 0.5*z2i + z3i/6.0 + z4i/24.0;
        *rho = PetscMax(*rho, PetscSqrtReal(pr*pr + pi*pi));
    }
    PetscCall(PetscFree4(Jcopy, wr, wi, work));
    PetscFunctionReturn(0);
}

/* sigma_max(J) = ||J||_2 (and optionally the dominant right singular vector v1). */
static PetscErrorCode DenseSigmaMax(const PetscReal *A, PetscInt n, PetscReal *smax, PetscReal *v1)
{
    PetscBLASInt N, lda, lwork, info;
    PetscReal *Acopy, *S, *VT, *work, ufake = 0.0;
    PetscFunctionBeginUser;
    PetscCall(PetscBLASIntCast(n, &N)); lda = N; lwork = 8*N + 4*N;
    PetscCall(PetscMalloc4(n*n, &Acopy, n, &S, n*n, &VT, (size_t)lwork, &work));
    for (PetscInt t = 0; t < n*n; t++) Acopy[t] = A[t];
    {
        char jobu = 'N', jobvt = v1 ? 'S' : 'N';
        LAPACKgesvd_(&jobu, &jobvt, &N, &N, Acopy, &lda, S, &ufake, &lda, VT, &lda, work, &lwork, &info);
    }
    PetscCheck(info == 0, PETSC_COMM_SELF, PETSC_ERR_LIB, "LAPACK dgesvd failed: info=%d", (int)info);
    *smax = S[0];
    if (v1) { for (PetscInt r = 0; r < n; r++) v1[r] = VT[0 + r*n]; }  /* first row of VT = v1^T */
    PetscCall(PetscFree4(Acopy, S, VT, work));
    PetscFunctionReturn(0);
}

/* Frobenius non-normality ||J^T J - J J^T||_F / max(||J||_F^2, eps). */
static PetscReal DenseNonNormality(const PetscReal *A, PetscInt n)
{
    PetscReal fro2 = 0.0, comm = 0.0;
    for (PetscInt t = 0; t < n*n; t++) fro2 += A[t]*A[t];
    for (PetscInt p = 0; p < n; p++)
        for (PetscInt q = 0; q < n; q++) {
            PetscReal ata = 0.0, aat = 0.0;
            for (PetscInt r = 0; r < n; r++) { ata += A[p + r*n]*A[q + r*n]; aat += A[r + p*n]*A[r + q*n]; }
            const PetscReal d = ata - aat; comm += d*d;
        }
    return PetscSqrtReal(comm) / PetscMax(fro2, PETSC_MACHINE_EPSILON);
}

/**
 * @brief Computes the normalized Frobenius defect from skew symmetry.
 */
static PetscReal DenseSkewnessDefect(const PetscReal *A, PetscInt n)
{
    PetscReal fro2 = 0.0, sym2 = 0.0;
    for (PetscInt i = 0; i < n*n; i++) fro2 += A[i]*A[i];
    for (PetscInt c = 0; c < n; c++)
        for (PetscInt r = 0; r < n; r++) {
            const PetscReal s = A[r + c*n] + A[c + r*n];
            sym2 += s*s;
        }
    return PetscSqrtReal(sym2) / PetscMax(PetscSqrtReal(fro2), PETSC_MACHINE_EPSILON);
}

/**
 * @brief Copies a dense matrix and adds a scalar shift to its diagonal.
 */
static PetscErrorCode DenseShiftIdentity(const PetscReal *A, PetscInt n, PetscReal shift, PetscReal *B)
{
    PetscFunctionBeginUser;
    PetscCall(PetscArraycpy(B, A, (size_t)n*n));
    for (PetscInt d = 0; d < n; d++) B[d + d*n] += shift;
    PetscFunctionReturn(0);
}

/* B = alpha*A (column-major) ; C = A*Bm ; returns into out (n x n). */
static void MatMul(const PetscReal *A, const PetscReal *B, PetscReal *out, PetscInt n)
{
    for (PetscInt c = 0; c < n; c++)
        for (PetscInt r = 0; r < n; r++) {
            PetscReal s = 0.0;
            for (PetscInt t = 0; t < n; t++) s += A[r + t*n]*B[t + c*n];
            out[r + c*n] = s;
        }
}

/* P(M)=I+M+M^2/2+M^3/6+M^4/24 via Horner: I+M(I+M/2(I+M/3(I+M/4))). M=dtau*J. */
static PetscErrorCode RKPolynomial(const PetscReal *J, PetscReal dtau, PetscInt n, PetscReal *P)
{
    PetscReal *M, *T, *T2;
    PetscFunctionBeginUser;
    PetscCall(PetscMalloc3(n*n, &M, n*n, &T, n*n, &T2));
    for (PetscInt t = 0; t < n*n; t++) M[t] = dtau*J[t];
    /* start S = I + M/4 */
    for (PetscInt t = 0; t < n*n; t++) T[t] = M[t]/4.0;
    for (PetscInt d = 0; d < n; d++) T[d + d*n] += 1.0;
    const PetscReal coef[3] = {3.0, 2.0, 1.0};   /* divide by 3, then 2, then 1 */
    for (int s = 0; s < 3; s++) {
        MatMul(M, T, T2, n);                      /* T2 = M*S */
        for (PetscInt t = 0; t < n*n; t++) T[t] = T2[t]/coef[s];
        for (PetscInt d = 0; d < n; d++) T[d + d*n] += 1.0;   /* S = I + M/coef * S */
    }
    for (PetscInt t = 0; t < n*n; t++) P[t] = T[t];
    PetscCall(PetscFree3(M, T, T2));
    PetscFunctionReturn(0);
}

typedef enum { METRIC_RHO, METRIC_NORM } PMetric;

/**
 * @brief Evaluates either spectral-radius or 2-norm amplification for one pseudo-time step.
 */
static PetscErrorCode AmplificationMetric(const PetscReal *J, PetscInt n, PetscReal dtau,
                                          PMetric which, PetscReal *metric)
{
    PetscReal *Pm;
    PetscFunctionBeginUser;
    if (which == METRIC_RHO) {
        PetscCall(DenseRKPolynomialSpectralRadius(J, n, dtau, metric));
        PetscFunctionReturn(0);
    }
    PetscCall(PetscMalloc1(n*n, &Pm));
    PetscCall(RKPolynomial(J, dtau, n, Pm));
    PetscCall(DenseSigmaMax(Pm, n, metric, NULL));
    PetscCall(PetscFree(Pm));
    PetscFunctionReturn(0);
}

/* Stable interval connected to CFL=0; tolerance=1e-8, initial probe=1e-8, scan step=0.01, max=4.0. */
static PetscErrorCode StableCFL(const PetscReal *J, PetscInt n, PetscReal lam, PMetric which,
                                StableCFLResult *result)
{
    const PetscReal tol = 1e-8, probe = 1e-8, scan_step = 0.01, hi = 4.0;
    const PetscReal probe_tol = 1e-12, min_positive_cfl = 1e-6;
    PetscReal met;
    PetscFunctionBeginUser;
    result->status = STABLE_CFL_NONE;
    result->cfl = 0.0;
    if (which == METRIC_RHO) {
        PetscReal rhoJ, maxreJ;
        PetscCall(DenseSpectralRadius(J, n, &rhoJ, &maxreJ));
        if (maxreJ > 1e-8) PetscFunctionReturn(0);
    }
    PetscCall(AmplificationMetric(J, n, probe/lam, which, &met));
    if (met > 1.0 + probe_tol) PetscFunctionReturn(0);

    PetscReal stable = probe, cross_b = -1.0;
    for (PetscReal cfl = scan_step; cfl <= hi + 1e-12; cfl += scan_step) {
        PetscCall(AmplificationMetric(J, n, cfl/lam, which, &met));
        if (met > 1.0 + tol) { cross_b = cfl; break; }
        stable = cfl;
    }
    if (cross_b < 0.0) {
        result->status = STABLE_CFL_EXCEEDS_SCAN;
        result->cfl = hi;
        PetscFunctionReturn(0);
    }
    PetscReal cross_a = stable;
    for (int it = 0; it < 40; it++) {
        const PetscReal mid = 0.5*(cross_a + cross_b);
        PetscCall(AmplificationMetric(J, n, mid/lam, which, &met));
        if (met > 1.0 + tol) cross_b = mid; else cross_a = mid;
    }
    if (cross_a < min_positive_cfl) PetscFunctionReturn(0);
    result->status = STABLE_CFL_FINITE;
    result->cfl = cross_a;
    PetscFunctionReturn(0);
}

/**
 * @brief Returns human-readable text for a stable-CFL search result.
 */
static const char *StableCFLStatusText(StableCFLResult r)
{
    if (r.status == STABLE_CFL_NONE) return "no positive stable interval connected to zero";
    if (r.status == STABLE_CFL_EXCEEDS_SCAN) return "stable through CFL >= 4.0";
    return "stable CFL";
}

/**
 * @brief Prints one candidate's eigenvalue and norm stable-CFL statuses.
 */
static PetscErrorCode PrintStableCFLLine(const char *candidate,
                                         StableCFLResult eig,
                                         StableCFLResult norm)
{
    PetscFunctionBeginUser;
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "    cand %s: eig %s", candidate, StableCFLStatusText(eig)));
    if (eig.status == STABLE_CFL_FINITE) PetscCall(PetscPrintf(PETSC_COMM_WORLD, " = %.4f", (double)eig.cfl));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, " | norm %s", StableCFLStatusText(norm)));
    if (norm.status == STABLE_CFL_FINITE) PetscCall(PetscPrintf(PETSC_COMM_WORLD, " = %.4f", (double)norm.cfl));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "\n"));
    PetscFunctionReturn(0);
}

/* P(M) applied to a vector: out = P(dtau*J) * x  (dense). */
static void ApplyP(const PetscReal *P, const PetscReal *x, PetscReal *out, PetscInt n)
{
    for (PetscInt r = 0; r < n; r++) { PetscReal s = 0.0; for (PetscInt c = 0; c < n; c++) s += P[r + c*n]*x[c]; out[r] = s; }
}

/**
 * @brief Prints frozen RK amplification tables for the supplied operator and candidates.
 */
static PetscErrorCode PrintFrozenAmplificationTable(const char *title, const PetscReal *J,
                                                    PetscInt n, const PetscReal lams[3],
                                                    const char *cn[3])
{
    const PetscReal cfls[5] = {0.25, 0.50, 1.00, 1.50, 2.00};
    PetscReal *Pm;
    PetscFunctionBeginUser;
    PetscCall(PetscMalloc1(n*n, &Pm));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "  --- %s ---\n"
        "    candidate  CFL    dtau          rho(P)        sigma_max(P)\n", title));
    for (int c = 0; c < 3; c++) {
        if (!(lams[c] > 0.0)) continue;
        for (int q = 0; q < 5; q++) {
            const PetscReal dtau = cfls[q]/lams[c];
            PetscReal rhoP, smaxP;
            PetscCall(RKPolynomial(J, dtau, n, Pm));
            PetscCall(DenseRKPolynomialSpectralRadius(J, n, dtau, &rhoP));
            PetscCall(DenseSigmaMax(Pm, n, &smaxP, NULL));
            PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                "    %-9s %.2f   %.6e  %.6e  %.6e\n",
                cn[c], (double)cfls[q], (double)dtau, (double)rhoP, (double)smaxP));
        }
    }
    PetscCall(PetscFree(Pm));
    PetscFunctionReturn(0);
}

/**
 * @brief Computes the Euclidean norm of a dense vector.
 */
static PetscReal VecNorm2Array(const PetscReal *x, PetscInt n)
{
    PetscReal s = 0.0;
    for (PetscInt i = 0; i < n; i++) s += x[i]*x[i];
    return PetscSqrtReal(s);
}

/**
 * @brief Computes the infinity norm of a dense vector.
 */
static PetscReal VecNormInfArray(const PetscReal *x, PetscInt n)
{
    PetscReal s = 0.0;
    for (PetscInt i = 0; i < n; i++) s = PetscMax(s, PetscAbsReal(x[i]));
    return s;
}

/**
 * @brief Computes the Frobenius norm of a dense column-major matrix.
 */
static PetscReal DenseFrobenius(const PetscReal *A, PetscInt n)
{
    PetscReal s = 0.0;
    for (PetscInt i = 0; i < n*n; i++) s += A[i]*A[i];
    return PetscSqrtReal(s);
}

/**
 * @brief Computes a normalized Frobenius difference between two dense matrices.
 */
static PetscReal DenseRelativeDiff(const PetscReal *A, const PetscReal *B, PetscInt n, PetscReal denom_ref)
{
    PetscReal s = 0.0;
    for (PetscInt i = 0; i < n*n; i++) {
        const PetscReal d = A[i] - B[i];
        s += d*d;
    }
    return PetscSqrtReal(s) / PetscMax(1.0, denom_ref);
}

/**
 * @brief Adds a dense active-space vector into a global contravariant vector.
 */
static PetscErrorCode AddActiveVector(UserCtx *user, Vec U, const DofMap *map,
                                      const PetscReal *x, PetscReal scale)
{
    Cmpnts ***a;
    PetscFunctionBeginUser;
    PetscCall(DMDAVecGetArray(user->fda, U, &a));
    for (PetscInt m = 0; m < map->n; m++) {
        PetscReal *p = (PetscReal*)&a[map->ck[m]][map->cj[m]][map->ci[m]];
        p[map->comp[m]] += scale*x[m];
    }
    PetscCall(DMDAVecRestoreArray(user->fda, U, &a));
    PetscFunctionReturn(0);
}

/**
 * @brief Extracts active-space entries from a global contravariant vector.
 */
static PetscErrorCode ExtractActiveVector(UserCtx *user, Vec U, const DofMap *map, PetscReal *x)
{
    Cmpnts ***a;
    PetscFunctionBeginUser;
    PetscCall(DMDAVecGetArrayRead(user->fda, U, &a));
    for (PetscInt m = 0; m < map->n; m++) {
        const PetscReal *p = (const PetscReal*)&a[map->ck[m]][map->cj[m]][map->ci[m]];
        x[m] = p[map->comp[m]];
    }
    PetscCall(DMDAVecRestoreArrayRead(user->fda, U, &a));
    PetscFunctionReturn(0);
}

/**
 * @brief Returns a deterministic checksum weight for an active DOF.
 */
static PetscReal DofWeight(const DofMap *map, PetscInt m)
{
    return 1.0 + 0.013*(PetscReal)(map->comp[m]+1)
               + 0.017*(PetscReal)map->ci[m]
               + 0.019*(PetscReal)map->cj[m]
               + 0.023*(PetscReal)map->ck[m];
}

/**
 * @brief Computes global active-vector norms and checksum.
 */
static PetscErrorCode ActiveStats(const DofMap *map, const PetscReal *x, GlobalVecStats *stats)
{
    PetscReal loc2 = 0.0, locinf = 0.0, locsum = 0.0;
    PetscReal glo2, gloinf, glosum;
    PetscFunctionBeginUser;
    for (PetscInt m = 0; m < map->n; m++) {
        loc2 += x[m]*x[m];
        locinf = PetscMax(locinf, PetscAbsReal(x[m]));
        locsum += DofWeight(map, m)*x[m];
    }
    PetscCallMPI(MPI_Allreduce(&loc2, &glo2, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Allreduce(&locinf, &gloinf, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Allreduce(&locsum, &glosum, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD));
    stats->n2 = PetscSqrtReal(glo2); stats->ninf = gloinf; stats->checksum = glosum;
    PetscFunctionReturn(0);
}

/**
 * @brief Fills a globally normalized deterministic active-space perturbation direction.
 */
static PetscErrorCode FillDeterministicDirection(UserCtx *user, Vec V, const DofMap *map, PetscReal *x)
{
    PetscReal loc2 = 0.0, glo2;
    PetscFunctionBeginUser;
    PetscCall(VecSet(V, 0.0));
    for (PetscInt m = 0; m < map->n; m++) {
        x[m] = PetscSinReal(0.37*(PetscReal)(map->ci[m]+1)
                          + 0.51*(PetscReal)(map->cj[m]+1)
                          + 0.73*(PetscReal)(map->ck[m]+1)
                          + 0.29*(PetscReal)(map->comp[m]+1));
        loc2 += x[m]*x[m];
    }
    PetscCallMPI(MPI_Allreduce(&loc2, &glo2, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD));
    const PetscReal invn = 1.0/PetscSqrtReal(glo2);
    for (PetscInt m = 0; m < map->n; m++) x[m] *= invn;
    PetscCall(AddActiveVector(user, V, map, x, 1.0));
    PetscFunctionReturn(0);
}

/* ----------------------------------------------------------------------------------- *
 *  Direct anchored 4-stage recurrence on the real residual (no controller).            *
 * ----------------------------------------------------------------------------------- */
/* Phi(Ufull) -> result stored in Uout (global), using EvalConvResidual at each stage. */
static PetscErrorCode FourStage(UserCtx *user, Vec U0full, PetscReal dtau, Vec Rhs, const DofMap *map,
                                PetscReal *Rscratch, Vec Uwork, Vec Uout)
{
    const PetscReal alfa[4] = {0.25, 1.0/3.0, 0.5, 1.0};
    PetscFunctionBeginUser;
    PetscCall(VecCopy(U0full, Uwork));                 /* stage state */
    for (int s = 0; s < 4; s++) {
        PetscCall(EvalConvResidual(user, Uwork, Rhs, map, Rscratch));   /* R(U^{(s)}) into active rows */
        PetscCall(VecCopy(U0full, Uout));               /* U^{(s+1)} = U0 + alfa*dtau*R */
        Cmpnts ***a;
        PetscCall(DMDAVecGetArray(user->fda, Uout, &a));
        for (PetscInt m = 0; m < map->n; m++) {
            PetscReal *p = (PetscReal*)&a[map->ck[m]][map->cj[m]][map->ci[m]];
            p[map->comp[m]] += alfa[s]*dtau*Rscratch[m];
        }
        PetscCall(DMDAVecRestoreArray(user->fda, Uout, &a));
    PetscCall(VecCopy(Uout, Uwork));
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Forms one anchored RK stage state from the base state and active residual.
 */
static PetscErrorCode SetAnchoredStage(UserCtx *user, Vec U0full, PetscReal scale,
                                       const DofMap *map, const PetscReal *Ract, Vec Ustage)
{
    PetscFunctionBeginUser;
    PetscCall(VecCopy(U0full, Ustage));
    PetscCall(AddActiveVector(user, Ustage, map, Ract, scale));
    PetscFunctionReturn(0);
}

/**
 * @brief Builds the first three anchored RK stage states for a base vector.
 */
static PetscErrorCode BuildStageStates(UserCtx *user, Vec U0full, PetscReal dtau, Vec Rhs,
                                       const DofMap *map, PetscReal *Rscratch,
                                       Vec Y1, Vec Y2, Vec Y3)
{
    const PetscReal alfa[3] = {0.25, 1.0/3.0, 0.5};
    PetscFunctionBeginUser;
    PetscCall(EvalConvResidual(user, U0full, Rhs, map, Rscratch));
    PetscCall(SetAnchoredStage(user, U0full, alfa[0]*dtau, map, Rscratch, Y1));
    PetscCall(EvalConvResidual(user, Y1, Rhs, map, Rscratch));
    PetscCall(SetAnchoredStage(user, U0full, alfa[1]*dtau, map, Rscratch, Y2));
    PetscCall(EvalConvResidual(user, Y2, Rhs, map, Rscratch));
    PetscCall(SetAnchoredStage(user, U0full, alfa[2]*dtau, map, Rscratch, Y3));
    PetscFunctionReturn(0);
}

/**
 * @brief Builds a centered finite-difference Jacobian of the production convective residual.
 */
static PetscErrorCode BuildFDJacobian(UserCtx *user, Vec Ucenter, PetscReal epsrel, Vec Rhs,
                                      const DofMap *map, PetscReal *Rp, PetscReal *Rm,
                                      Vec Uwork, PetscReal *J)
{
    PetscFunctionBeginUser;
    for (PetscInt col = 0; col < map->n; col++) {
        PetscReal u0;
        PetscCall(GetDof(user, Ucenter, map, col, &u0));
        const PetscReal eps = epsrel*PetscMax(1.0, PetscAbsReal(u0));
        PetscCall(VecCopy(Ucenter, Uwork));
        PetscCall(PerturbDof(user, Uwork, map, col, +eps));
        PetscCall(EvalConvResidual(user, Uwork, Rhs, map, Rp));
        PetscCall(VecCopy(Ucenter, Uwork));
        PetscCall(PerturbDof(user, Uwork, map, col, -eps));
        PetscCall(EvalConvResidual(user, Uwork, Rhs, map, Rm));
        for (PetscInt row = 0; row < map->n; row++) J[row + col*map->n] = (Rp[row]-Rm[row])/(2.0*eps);
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Builds the exact four-stage tangent from stage-dependent Jacobians.
 */
static PetscErrorCode BuildStageTangent(const PetscReal *J0, const PetscReal *J1,
                                        const PetscReal *J2, const PetscReal *J3,
                                        PetscReal dtau, PetscInt n, PetscReal *T4)
{
    const PetscReal alfa[4] = {0.25, 1.0/3.0, 0.5, 1.0};
    PetscReal *T1, *T2, *T3, *Tmp;
    PetscFunctionBeginUser;
    PetscCall(PetscMalloc4(n*n, &T1, n*n, &T2, n*n, &T3, n*n, &Tmp));

    for (PetscInt i = 0; i < n*n; i++) T1[i] = alfa[0]*dtau*J0[i];
    for (PetscInt d = 0; d < n; d++) T1[d + d*n] += 1.0;

    MatMul(J1, T1, Tmp, n);
    for (PetscInt i = 0; i < n*n; i++) T2[i] = alfa[1]*dtau*Tmp[i];
    for (PetscInt d = 0; d < n; d++) T2[d + d*n] += 1.0;

    MatMul(J2, T2, Tmp, n);
    for (PetscInt i = 0; i < n*n; i++) T3[i] = alfa[2]*dtau*Tmp[i];
    for (PetscInt d = 0; d < n; d++) T3[d + d*n] += 1.0;

    MatMul(J3, T3, Tmp, n);
    for (PetscInt i = 0; i < n*n; i++) T4[i] = alfa[3]*dtau*Tmp[i];
    for (PetscInt d = 0; d < n; d++) T4[d + d*n] += 1.0;

    PetscCall(PetscFree4(T1, T2, T3, Tmp));
    PetscFunctionReturn(0);
}

/**
 * @brief Builds a finite-difference Jacobian of the complete nonlinear four-stage map.
 */
static PetscErrorCode BuildPhiJacobian(UserCtx *user, Vec Ucenter, PetscReal dtau, Vec Rhs,
                                       const DofMap *map, PetscReal *Rscratch,
                                       Vec Upert, Vec Ustage, Vec PhiP, Vec PhiM,
                                       PetscReal *xp, PetscReal *xm, PetscReal epsrel,
                                       PetscReal *JPhi)
{
    PetscFunctionBeginUser;
    for (PetscInt col = 0; col < map->n; col++) {
        PetscReal u0;
        PetscCall(GetDof(user, Ucenter, map, col, &u0));
        const PetscReal eps = epsrel*PetscMax(1.0, PetscAbsReal(u0));
        PetscCall(VecCopy(Ucenter, Upert));
        PetscCall(PerturbDof(user, Upert, map, col, +eps));
        PetscCall(FourStage(user, Upert, dtau, Rhs, map, Rscratch, Ustage, PhiP));
        PetscCall(VecCopy(Ucenter, Upert));
        PetscCall(PerturbDof(user, Upert, map, col, -eps));
        PetscCall(FourStage(user, Upert, dtau, Rhs, map, Rscratch, Ustage, PhiM));
        PetscCall(ExtractActiveVector(user, PhiP, map, xp));
        PetscCall(ExtractActiveVector(user, PhiM, map, xm));
        for (PetscInt row = 0; row < map->n; row++) JPhi[row + col*map->n] = (xp[row]-xm[row])/(2.0*eps);
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Runs stage-dependent RK tangent and direct nonlinear perturbation diagnostics.
 */
static PetscErrorCode RunRKTangentDiagnostics(UserCtx *user, CandState st, Vec Ubase, Vec Rhs,
                                              const DofMap *map, PetscReal epsrel,
                                              const PetscReal lams[3], const char *cn[3],
                                              const PetscReal *J0)
{
    const PetscReal cflsB[4] = {0.1, 0.25, 0.5, 1.0};
    const PetscReal cflsOther[1] = {0.5};
    const PetscReal *cfls = (st == STATE_B) ? cflsB : cflsOther;
    const PetscInt ncfl = (st == STATE_B) ? 4 : 1;
    const PetscReal amps[3] = {1e-4, 1e-5, 1e-6};
    Vec Y1, Y2, Y3, Upert, Ustage, PhiP, PhiM, Phi0;
    PetscReal *Rtmp, *J1, *J2, *J3, *T4, *Pm, *JPhi, *v1, *xrand;
    PetscReal *meas, *phi0, *phip, *xscaled, *predP, *predT;
    PetscFunctionBeginUser;

    PetscCall(VecDuplicate(Ubase, &Y1));
    PetscCall(VecDuplicate(Ubase, &Y2));
    PetscCall(VecDuplicate(Ubase, &Y3));
    PetscCall(VecDuplicate(Ubase, &Upert));
    PetscCall(VecDuplicate(Ubase, &Ustage));
    PetscCall(VecDuplicate(Ubase, &PhiP));
    PetscCall(VecDuplicate(Ubase, &PhiM));
    PetscCall(VecDuplicate(Ubase, &Phi0));
    PetscCall(PetscMalloc5(map->n, &Rtmp, (size_t)map->n*map->n, &J1,
                           (size_t)map->n*map->n, &J2, (size_t)map->n*map->n, &J3,
                           (size_t)map->n*map->n, &T4));
    PetscCall(PetscMalloc5((size_t)map->n*map->n, &Pm, (size_t)map->n*map->n, &JPhi,
                           map->n, &v1, map->n, &xrand, map->n, &meas));
    PetscCall(PetscMalloc5(map->n, &phi0, map->n, &phip, map->n, &xscaled,
                           map->n, &predP, map->n, &predT));

    PetscReal nrm = 0.0;
    for (PetscInt m = 0; m < map->n; m++) {
        xrand[m] = PetscSinReal((PetscReal)(m+1)*1.2345);
        nrm += xrand[m]*xrand[m];
    }
    nrm = PetscSqrtReal(nrm);
    for (PetscInt m = 0; m < map->n; m++) xrand[m] /= nrm;

    for (int cand = 0; cand < 3; cand++) {
      if (!(lams[cand] > 0.0)) continue;
      for (PetscInt icfl = 0; icfl < ncfl; icfl++) {
        const PetscReal cfl = cfls[icfl], dtau = cfl/lams[cand];
        PetscCall(BuildStageStates(user, Ubase, dtau, Rhs, map, Rtmp, Y1, Y2, Y3));
        PetscCall(BuildFDJacobian(user, Y1, epsrel, Rhs, map, phi0, phip, Upert, J1));
        PetscCall(BuildFDJacobian(user, Y2, epsrel, Rhs, map, phi0, phip, Upert, J2));
        PetscCall(BuildFDJacobian(user, Y3, epsrel, Rhs, map, phi0, phip, Upert, J3));
        PetscCall(BuildStageTangent(J0, J1, J2, J3, dtau, map->n, T4));
        PetscCall(RKPolynomial(J0, dtau, map->n, Pm));
        PetscCall(BuildPhiJacobian(user, Ubase, dtau, Rhs, map, Rtmp,
                                   Upert, Ustage, PhiP, PhiM, phip, phi0, epsrel, JPhi));

        const PetscReal froT = DenseFrobenius(T4, map->n);
        const PetscReal froPhi = DenseFrobenius(JPhi, map->n);
        const PetscReal relTP = DenseRelativeDiff(T4, Pm, map->n, froT);
        const PetscReal relPhiT = DenseRelativeDiff(JPhi, T4, map->n, froPhi);
        const PetscReal relPhiP = DenseRelativeDiff(JPhi, Pm, map->n, froPhi);
        PetscReal smaxT, rhoT, dummyT;
        PetscCall(DenseSpectralRadius(T4, map->n, &rhoT, &dummyT));
        PetscCall(DenseSigmaMax(T4, map->n, &smaxT, v1));

        PetscCall(PetscPrintf(PETSC_COMM_WORLD,
            "  --- convection-only stage-dependent RK tangent (cand %s, CFL=%.2f, dtau=%.6e) ---\n"
            "    rho(T4)=%.6e  sigma_max(T4)=%.6e\n"
            "    ||T4-P(hJ0)||F/max(1,||T4||F)      = %.3e\n"
            "    ||J_Phi-T4||F/max(1,||J_Phi||F)    = %.3e\n"
            "    ||J_Phi-P(hJ0)||F/max(1,||J_Phi||F)= %.3e\n",
            cn[cand], (double)cfl, (double)dtau, (double)rhoT, (double)smaxT,
            (double)relTP, (double)relPhiT, (double)relPhiP));

        PetscCall(FourStage(user, Ubase, dtau, Rhs, map, Rtmp, Ustage, Phi0));
        PetscCall(ExtractActiveVector(user, Phi0, map, phi0));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD,
            "    direction  amp    amp(meas)  amp(frozen)  amp(stage)  err(frozen)  err(stage)\n"));
        for (int dir = 0; dir < 2; dir++) {
            const PetscReal *xv = (dir == 0) ? xrand : v1;
            const char *dname = (dir == 0) ? "random" : "v1(T4)";
            for (int a = 0; a < 3; a++) {
                const PetscReal amp = amps[a];
                PetscCall(VecCopy(Ubase, Upert));
                PetscCall(AddActiveVector(user, Upert, map, xv, amp));
                PetscCall(FourStage(user, Upert, dtau, Rhs, map, Rtmp, Ustage, PhiP));
                PetscCall(ExtractActiveVector(user, PhiP, map, phip));
                for (PetscInt m = 0; m < map->n; m++) meas[m] = phip[m] - phi0[m];
                for (PetscInt m = 0; m < map->n; m++) xscaled[m] = amp*xv[m];
                ApplyP(Pm, xscaled, predP, map->n);
                ApplyP(T4, xscaled, predT, map->n);

                PetscReal eP = 0.0, eT = 0.0;
                for (PetscInt m = 0; m < map->n; m++) {
                    const PetscReal dP = meas[m] - predP[m];
                    const PetscReal dT = meas[m] - predT[m];
                    eP += dP*dP; eT += dT*dT;
                }
                const PetscReal nm = VecNorm2Array(meas, map->n);
                const PetscReal nP = VecNorm2Array(predP, map->n);
                const PetscReal nT = VecNorm2Array(predT, map->n);
                PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                    "    %-8s %.0e  %.6e  %.6e  %.6e  %.3e  %.3e\n",
                    dname, (double)amp, (double)(nm/amp), (double)(nP/amp), (double)(nT/amp),
                    (double)(PetscSqrtReal(eP)/PetscMax(PETSC_MACHINE_EPSILON, nP)),
                    (double)(PetscSqrtReal(eT)/PetscMax(PETSC_MACHINE_EPSILON, nT))));
            }
        }

        PetscCall(PicurvAssertBool((PetscBool)(relPhiT < 1e-6),
                                   "stage tangent matches direct finite-difference RK map"));
        if (st == STATE_B) {
            PetscCall(PicurvAssertBool((PetscBool)(relTP > 1e-5),
                                       "B: non-steady base makes frozen-Jacobian RK map measurably different"));
            PetscCall(PicurvAssertBool((PetscBool)(relPhiP > 1e-5),
                                       "B: direct RK map confirms frozen-Jacobian error"));
        } else {
            PetscCall(PicurvAssertBool((PetscBool)(relTP < 1e-8),
                                       "steady base reduces stage tangent to frozen RK polynomial"));
        }
      }
    }

    PetscCall(PetscFree5(Rtmp, J1, J2, J3, T4));
    PetscCall(PetscFree5(Pm, JPhi, v1, xrand, meas));
    PetscCall(PetscFree5(phi0, phip, xscaled, predP, predT));
    PetscCall(VecDestroy(&Y1)); PetscCall(VecDestroy(&Y2)); PetscCall(VecDestroy(&Y3));
    PetscCall(VecDestroy(&Upert)); PetscCall(VecDestroy(&Ustage));
    PetscCall(VecDestroy(&PhiP)); PetscCall(VecDestroy(&PhiM)); PetscCall(VecDestroy(&Phi0));
    PetscFunctionReturn(0);
}

/* ----------------------------------------------------------------------------------- *
 *  The A4a study for one state.                                                        *
 * ----------------------------------------------------------------------------------- */
static PetscErrorCode RunState(CandState st, const char *name)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL; DofMap map;
    Vec Ubase, Rhs, Uwork, Uout;
    PetscReal *Rref, *Rrep, *Jbest;
    PetscReal repeat_err, maxdiv, det_err;
    SeamDiagnostics seam;
    MomStabilityReport rep;
    const PetscInt N = (st == STATE_C) ? 5 : 4;   /* C uses one extra point for canonical shear. */
    PetscFunctionBeginUser;

    /* periodic Cartesian fixture, inviscid + centered + P=0, no LES/RANS/Clark/IB/body force. */
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, N, N, N, PETSC_TRUE, PETSC_TRUE, PETSC_TRUE));
    PetscCall(ConfigureCandidateFixture(simCtx, user));

    PetscCall(DofMapBuild(user, &map));
    PetscCall(VecDuplicate(user->Ucont, &Ubase));
    PetscCall(VecDuplicate(user->Ucont, &Rhs));
    PetscCall(VecDuplicate(user->Ucont, &Uwork));
    PetscCall(VecDuplicate(user->Ucont, &Uout));
    PetscCall(PetscMalloc3(map.n, &Rref, map.n, &Rrep, (size_t)map.n*map.n, &Jbest));

    PetscCall(BuildBaseState(user, st, Ubase, &repeat_err, &maxdiv, &seam));

    /* residual determinism: evaluate the base residual twice. */
    PetscCall(EvalConvResidual(user, Ubase, Rhs, &map, Rref));
    PetscCall(EvalConvResidual(user, Ubase, Rhs, &map, Rrep));
    det_err = 0.0; for (PetscInt m = 0; m < map.n; m++) det_err = PetscMax(det_err, PetscAbsReal(Rref[m]-Rrep[m]));
    const PetscReal R0_2 = VecNorm2Array(Rref, map.n);
    const PetscReal R0_inf = VecNormInfArray(Rref, map.n);

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "\n================ STATE %s ================\n"
        "  grid: DMDA %d^3 (periodic) | independent face DOFs: actual=%d expected=%d\n"
        "  declared endpoint mismatch       (x,y,z) = %.3e %.3e %.3e\n"
        "  actual Ucat duplicate mismatch   (x,y,z) = %.3e %.3e %.3e\n"
        "  local lUcat ghost mismatch       (x,y,z) = %.3e %.3e %.3e\n"
        "  actual Ucont duplicate mismatch  (x,y,z) = %.3e %.3e %.3e\n"
        "  ||Ucat_repeat - Ucat_reference||inf = %.3e\n"
        "  max|div_h Ucont|                    = %.3e\n"
        "  ||R(U0)||2 = %.6e  ||R(U0)||inf = %.6e\n"
        "  residual determinism ||R_rep-R_ref||inf = %.3e\n",
        name, (int)(N+1), (int)map.n,
        (int)map.expected_n,
        (double)seam.declared[0], (double)seam.declared[1], (double)seam.declared[2],
        (double)seam.ucat_global[0], (double)seam.ucat_global[1], (double)seam.ucat_global[2],
        (double)seam.ucat_ghost[0], (double)seam.ucat_ghost[1], (double)seam.ucat_ghost[2],
        (double)seam.ucont_global[0], (double)seam.ucont_global[1], (double)seam.ucont_global[2],
        (double)repeat_err, (double)maxdiv,
        (double)R0_2, (double)R0_inf, (double)det_err));

    /* ---- epsilon-convergence study: build J at several eps, compare to next-finer ---- */
    const PetscReal epsrel[5] = {1e-4, 1e-5, 1e-6, 1e-7, 1e-8};
    PetscReal *Jprev, *Jcur;
    PetscCall(PetscMalloc2((size_t)map.n*map.n, &Jprev, (size_t)map.n*map.n, &Jcur));
    PetscReal best_rel = PETSC_MAX_REAL; PetscInt best_e = 2;
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  epsilon-convergence (||J_e - J_e/10||_F / ||J||_F):\n"));
    for (int e = 0; e < 5; e++) {
        /* build column-major J at epsrel[e] */
        for (PetscInt col = 0; col < map.n; col++) {
            PetscReal u0; PetscCall(GetDof(user, Ubase, &map, col, &u0));
            const PetscReal eps = epsrel[e]*PetscMax(1.0, PetscAbsReal(u0));
            PetscReal *Rp = Rref, *Rm = Rrep;   /* reuse scratch */
            PetscCall(VecCopy(Ubase, Uwork));
            PetscCall(PerturbDof(user, Uwork, &map, col, +eps));
            PetscCall(EvalConvResidual(user, Uwork, Rhs, &map, Rp));
            PetscCall(VecCopy(Ubase, Uwork));
            PetscCall(PerturbDof(user, Uwork, &map, col, -eps));
            PetscCall(EvalConvResidual(user, Uwork, Rhs, &map, Rm));
            for (PetscInt row = 0; row < map.n; row++) Jcur[row + col*map.n] = (Rp[row]-Rm[row])/(2.0*eps);
        }
        if (e > 0) {
            PetscReal num = 0.0, den = 0.0;
            for (PetscInt t = 0; t < map.n*map.n; t++) { const PetscReal d = Jprev[t]-Jcur[t]; num += d*d; den += Jcur[t]*Jcur[t]; }
            const PetscReal rel = PetscSqrtReal(num)/PetscMax(1.0, PetscSqrtReal(den));
            PetscCall(PetscPrintf(PETSC_COMM_WORLD, "    eps=%.0e -> %.0e : rel=%.3e\n",
                                  (double)epsrel[e-1], (double)epsrel[e], (double)rel));
            if (rel < best_rel) { best_rel = rel; best_e = e; }
        }
        if (st == STATE_A) {
            PetscReal erho, emaxre, esmax;
            PetscCall(DenseSpectralRadius(Jcur, map.n, &erho, &emaxre));
            PetscCall(DenseSigmaMax(Jcur, map.n, &esmax, NULL));
            PetscCall(PetscPrintf(PETSC_COMM_WORLD,
                "    eps=%.0e metrics: rho=%.6e sigma=%.6e max_real=%.6e fro=%.6e skew=%.3e\n",
                (double)epsrel[e], (double)erho, (double)esmax, (double)emaxre,
                (double)DenseFrobenius(Jcur, map.n), (double)DenseSkewnessDefect(Jcur, map.n)));
        }
        PetscCall(PetscArraycpy(Jprev, Jcur, (size_t)map.n*map.n));   /* prev <- current (no pointer swap) */
    }
    /* rebuild J at the plateau epsilon into Jbest */
    {
        const PetscReal er = epsrel[best_e];
        for (PetscInt col = 0; col < map.n; col++) {
            PetscReal u0; PetscCall(GetDof(user, Ubase, &map, col, &u0));
            const PetscReal eps = er*PetscMax(1.0, PetscAbsReal(u0));
            PetscCall(VecCopy(Ubase, Uwork)); PetscCall(PerturbDof(user, Uwork, &map, col, +eps));
            PetscCall(EvalConvResidual(user, Uwork, Rhs, &map, Rref));
            PetscCall(VecCopy(Ubase, Uwork)); PetscCall(PerturbDof(user, Uwork, &map, col, -eps));
            PetscCall(EvalConvResidual(user, Uwork, Rhs, &map, Rrep));
            for (PetscInt row = 0; row < map.n; row++) Jbest[row + col*map.n] = (Rref[row]-Rrep[row])/(2.0*eps);
        }
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  selected plateau eps = %.0e\n", (double)er));
    }
    PetscCall(PetscFree2(Jprev, Jcur));

    /* The FD probes leave the production vectors at the final perturbed evaluation;
       restore the base state before every downstream diagnostic. */
    PetscCall(VecCopy(Ubase, user->Ucont));
    {
        const char *fld[] = {"Ucont"};
        PetscCall(SynchronizePeriodicStaggeredFields(user, 1, fld));
    }

    /* base state must be restored after the Jacobian build. */
    PetscBool eqbase; { Vec chk; PetscCall(VecDuplicate(Ubase,&chk)); PetscCall(VecCopy(user->Ucont,chk));
        PetscCall(VecAXPY(chk, -1.0, Ubase)); PetscReal nb; PetscCall(VecNorm(chk, NORM_INFINITY, &nb));
        eqbase = (PetscBool)(nb < 1e-12); PetscCall(VecDestroy(&chk));
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  base Ucont restored after Jacobian: %s\n", eqbase?"yes":"NO")); }

    /* The FD operator acts on Ucont; J is on the convective residual (sign per ComputeRHS). */
    PetscReal rho, maxre, smax;
    PetscCall(DenseSpectralRadius(Jbest, map.n, &rho, &maxre));
    PetscCall(DenseSigmaMax(Jbest, map.n, &smax, NULL));
    const PetscReal eta = DenseNonNormality(Jbest, map.n);

    /* ---- candidate estimates (convection-only: lambda_cX = lambda_X - lambda_t) ---- */
    PetscCall(EvalConvResidual(user, Ubase, Rhs, &map, Rref));   /* restore lUcat consistent w/ base */
    PetscCall(UpdateLocalGhosts(user, "Ucont"));
    PetscCall(Contra2Cart(user)); PetscCall(UpdateLocalGhosts(user, "Ucat"));
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_C, &rep));
    const PetscReal lcB = rep.lambda_B - rep.lambda_t;
    const PetscReal lcC = rep.lambda_C - rep.lambda_t;
    const PetscReal lcD = rep.lambda_D - rep.lambda_t;
    PetscReal gradmax = 0.0;
    PetscCall(ComputeMaxGradientContribution(user, &gradmax));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "  --- convective Jacobian spectrum ---\n"
        "    rho(J)=%.6e  sigma_max(J)=%.6e  nonnormality=%.3e  max_real_eig=%.3e\n"
        "  --- candidate convective estimates ---\n"
        "    lambda_cB=%.6e  lambda_cC=%.6e  lambda_cD=%.6e\n"
        "    max local |grad u| row-sum contribution = %.6e\n"
        "    rB/rho=%.3f rC/rho=%.3f rD/rho=%.3f | rB/smax=%.3f rC/smax=%.3f rD/smax=%.3f\n",
        (double)rho, (double)smax, (double)eta, (double)maxre,
        (double)lcB, (double)lcC, (double)lcD,
        (double)gradmax,
        (double)(lcB/rho), (double)(lcC/rho), (double)(lcD/rho),
        (double)(lcB/smax), (double)(lcC/smax), (double)(lcD/smax)));

    /* ---- RK pseudo-CFL stability per candidate (convective Jacobian) ---- */
    const PetscReal lams[3] = {lcB, lcC, lcD}; const char *cn[3] = {"B","C","D"};
    PetscReal *Jpseudo;
    PetscCall(PetscMalloc1((size_t)map.n*map.n, &Jpseudo));
    PetscCall(DenseShiftIdentity(Jbest, map.n, -rep.lambda_t, Jpseudo));
    const PetscReal lams_full[3] = {rep.lambda_B, rep.lambda_C, rep.lambda_D};
    PetscCall(PrintFrozenAmplificationTable("CONVECTION-ONLY FROZEN JACOBIAN", Jbest, map.n, lams, cn));
    PetscCall(PrintFrozenAmplificationTable("COMPLETE PSEUDO-TIME FROZEN OPERATOR", Jpseudo, map.n, lams_full, cn));
    PetscCall(PetscFree(Jpseudo));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "  --- RK stable pseudo-CFL (convective J) ---\n"));
    for (int c = 0; c < 3; c++) {
        if (!(lams[c] > 0.0)) continue;
        StableCFLResult cfl_rho, cfl_nrm;
        PetscCall(StableCFL(Jbest, map.n, lams[c], METRIC_RHO,  &cfl_rho));
        PetscCall(StableCFL(Jbest, map.n, lams[c], METRIC_NORM, &cfl_nrm));
        PetscCall(PrintStableCFLLine(cn[c], cfl_rho, cfl_nrm));
    }

    /* ---- exact stage-dependent tangent map and direct nonlinear RK cross-check ---- */
    if (lcC > 0.0) PetscCall(RunRKTangentDiagnostics(user, st, Ubase, Rhs, &map,
                                                     epsrel[best_e], lams, cn, Jbest));

    /* ---- robust automated assertions per state ---- */
    PetscCall(PicurvAssertBool(eqbase, "base Ucont restored after Jacobian build"));
    PetscCall(PicurvAssertBool((PetscBool)(det_err < 1e-10), "residual evaluation deterministic"));
    PetscCall(PicurvAssertBool((PetscBool)(rho > 0.0 && PetscIsNormalReal(rho)), "finite rho(J)"));
    PetscCall(PicurvAssertBool((PetscBool)(smax > 0.0 && PetscIsNormalReal(smax)), "finite sigma_max(J)"));
    if (st == STATE_A) {
        for (int d = 0; d < 3; d++) {
            PetscCall(PicurvAssertRealNear(seam.declared[d], 0.0, 1e-12, "A: declared periodic seam"));
            PetscCall(PicurvAssertRealNear(seam.ucat_global[d], 0.0, 1e-12, "A: Ucat duplicate seam"));
            PetscCall(PicurvAssertRealNear(seam.ucat_ghost[d], 0.0, 1e-12, "A: lUcat ghost seam"));
            PetscCall(PicurvAssertRealNear(seam.ucont_global[d], 0.0, 1e-12, "A: Ucont duplicate seam"));
        }
        PetscCall(PicurvAssertRealNear(repeat_err, 0.0, 1e-9, "A: recovered Ucat repeat near roundoff"));
        PetscCall(PicurvAssertRealNear(R0_2, 0.0, 1e-12, "A: base residual 2-norm near roundoff"));
        PetscCall(PicurvAssertRealNear(R0_inf, 0.0, 1e-12, "A: base residual inf-norm near roundoff"));
        PetscCall(PicurvAssertRealNear(lcB, lcC, 1e-9, "A: B == C (no divergence)"));
        PetscCall(PicurvAssertRealNear(lcC, lcD, 1e-9, "A: C == D (no shear)"));
    } else if (st == STATE_B) {
        for (int d = 0; d < 3; d++) {
            PetscCall(PicurvAssertRealNear(seam.declared[d], 0.0, 1e-12, "B: declared periodic seam"));
            PetscCall(PicurvAssertRealNear(seam.ucat_global[d], 0.0, 1e-12, "B: Ucat duplicate seam"));
            PetscCall(PicurvAssertRealNear(seam.ucat_ghost[d], 0.0, 1e-12, "B: lUcat ghost seam"));
            PetscCall(PicurvAssertRealNear(seam.ucont_global[d], 0.0, 1e-12, "B: Ucont duplicate seam"));
        }
        PetscCall(PicurvAssertRealNear(repeat_err, 0.0, 1e-9, "B: recovered Ucat repeat near roundoff"));
        PetscCall(PicurvAssertBool((PetscBool)(maxdiv > 1e-3), "B: nonzero discrete divergence"));
        PetscCall(PicurvAssertBool((PetscBool)(R0_2 > 1e-3), "B: materially nonzero base residual 2-norm"));
        PetscCall(PicurvAssertBool((PetscBool)(R0_inf > 1e-3), "B: materially nonzero base residual inf-norm"));
        PetscCall(PicurvAssertBool((PetscBool)(lcC > lcB + 1e-9), "B: C > B"));
    } else {
        for (int d = 0; d < 3; d++) {
            PetscCall(PicurvAssertRealNear(seam.declared[d], 0.0, 1e-12, "C: declared periodic seam"));
            PetscCall(PicurvAssertRealNear(seam.ucat_global[d], 0.0, 1e-12, "C: Ucat duplicate seam"));
            PetscCall(PicurvAssertRealNear(seam.ucat_ghost[d], 0.0, 1e-12, "C: lUcat ghost seam"));
            PetscCall(PicurvAssertRealNear(seam.ucont_global[d], 0.0, 1e-12, "C: Ucont duplicate seam"));
        }
        PetscCall(PicurvAssertRealNear(repeat_err, 0.0, 1e-9, "C: recovered Ucat repeat near roundoff"));
        PetscCall(PicurvAssertBool((PetscBool)(maxdiv < 1e-6), "C: divergence near zero"));
        PetscCall(PicurvAssertRealNear(R0_2, 0.0, 1e-12, "C: base residual 2-norm near roundoff"));
        PetscCall(PicurvAssertRealNear(R0_inf, 0.0, 1e-12, "C: base residual inf-norm near roundoff"));
        PetscCall(PicurvAssertBool((PetscBool)(PetscAbsReal(lcC-lcB) < 1e-6), "C: C ~= B"));
        PetscCall(PicurvAssertBool((PetscBool)(gradmax > 1e-6), "C: canonical shear has nonzero local gradient contribution"));
    }

    PetscCall(PetscFree3(Rref, Rrep, Jbest));
    PetscCall(VecDestroy(&Ubase)); PetscCall(VecDestroy(&Rhs));
    PetscCall(VecDestroy(&Uwork)); PetscCall(VecDestroy(&Uout));
    PetscCall(DofMapDestroy(&map));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Runs the State A candidate harness.
 */
static PetscErrorCode TestStateA(void) { return RunState(STATE_A, "A (uniform div-free)"); }

/**
 * @brief Runs the State B candidate harness.
 */
static PetscErrorCode TestStateB(void) { return RunState(STATE_B, "B (nonzero divergence)"); }

/**
 * @brief Runs the State C candidate harness.
 */
static PetscErrorCode TestStateC(void) { return RunState(STATE_C, "C (div-free shear)"); }

/**
 * @brief Runs one State A grid-size audit case.
 */
static PetscErrorCode RunStateAGridAuditOne(PetscInt N)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL; DofMap map;
    Vec Ubase, Rhs, Uwork;
    PetscReal *Rp, *Rm, *J;
    PetscReal repeat_err, maxdiv;
    SeamDiagnostics seam;
    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, N, N, N, PETSC_TRUE, PETSC_TRUE, PETSC_TRUE));
    PetscCall(ConfigureCandidateFixture(simCtx, user));
    PetscCall(DofMapBuild(user, &map));
    PetscCall(VecDuplicate(user->Ucont, &Ubase));
    PetscCall(VecDuplicate(user->Ucont, &Rhs));
    PetscCall(VecDuplicate(user->Ucont, &Uwork));
    PetscCall(PetscMalloc3(map.n, &Rp, map.n, &Rm, (size_t)map.n*map.n, &J));
    PetscCall(BuildBaseState(user, STATE_A, Ubase, &repeat_err, &maxdiv, &seam));
    PetscCall(BuildFDJacobian(user, Ubase, 1e-5, Rhs, &map, Rp, Rm, Uwork, J));
    PetscReal rho, maxre, smax;
    PetscCall(DenseSpectralRadius(J, map.n, &rho, &maxre));
    PetscCall(DenseSigmaMax(J, map.n, &smax, NULL));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "  State A grid audit: DMDA %d^3 independent=%d expected=%d max_real=%.6e rho=%.6e sigma=%.6e skew=%.3e nonnormality=%.3e repeat=%.3e\n",
        (int)(N+1), (int)map.n, (int)map.expected_n, (double)maxre, (double)rho, (double)smax,
        (double)DenseSkewnessDefect(J, map.n), (double)DenseNonNormality(J, map.n), (double)repeat_err));
    PetscCall(PetscFree3(Rp, Rm, J));
    PetscCall(VecDestroy(&Ubase)); PetscCall(VecDestroy(&Rhs)); PetscCall(VecDestroy(&Uwork));
    PetscCall(DofMapDestroy(&map));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Runs the State A grid-dependence and active-space audit.
 */
static PetscErrorCode TestStateAGridAudit(void)
{
    PetscFunctionBeginUser;
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "\n================ STATE A GRID / ACTIVE-SPACE AUDIT ================\n"
        "  residual sign convention: ComputeRHS returns the production convection residual used by pseudo-time updates.\n"
        "  component-staggered periodic duplicate planes: 0<-m-2 and m-1<-1 in each direction.\n"
        "  independent map: all three Ucont components use representatives i,j,k=1..m-2; count = 3*(m-2)^3.\n"));
    PetscCall(RunStateAGridAuditOne(4));
    PetscCall(RunStateAGridAuditOne(5));
    PetscFunctionReturn(0);
}

/**
 * @brief Writes the one-rank State A matrix-free decomposition reference.
 */
static PetscErrorCode WriteStateADecompReference(GlobalVecStats r0, GlobalVecStats jv,
                                                 GlobalVecStats phi, const MomStabilityReport *rep)
{
    PetscMPIInt rank;
    PetscFunctionBeginUser;
    if (!g_ref_path_set) PetscFunctionReturn(0);
    PetscCheck(g_ref_token_set, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE,
               "-candidate_ref_token is required when -candidate_ref_path is set");
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    if (rank == 0) {
        FILE *fp = fopen(g_ref_path, "w");
        PetscCheck(fp != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "could not write State A MPI reference");
        fprintf(fp, "PICURV_CANDIDATE_STATEA_REF_V2 %s\n", g_ref_token);
        fprintf(fp, "%.17e %.17e %.17e\n", (double)r0.n2, (double)r0.ninf, (double)r0.checksum);
        fprintf(fp, "%.17e %.17e %.17e\n", (double)jv.n2, (double)jv.ninf, (double)jv.checksum);
        fprintf(fp, "%.17e %.17e %.17e\n", (double)phi.n2, (double)phi.ninf, (double)phi.checksum);
        fprintf(fp, "%.17e %.17e %.17e %d %d %d %d\n",
                (double)(rep->lambda_B - rep->lambda_t),
                (double)(rep->lambda_C - rep->lambda_t),
                (double)(rep->lambda_D - rep->lambda_t),
                (int)rep->active_cells, (int)rep->cblock, (int)rep->ci, (int)rep->cj);
        fclose(fp);
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Compares a distributed State A matrix-free check against the one-rank reference.
 */
static PetscErrorCode ReadAndCompareStateADecompReference(GlobalVecStats r0, GlobalVecStats jv,
                                                          GlobalVecStats phi, const MomStabilityReport *rep)
{
    PetscMPIInt rank;
    PetscReal vals[16] = {0.0};
    PetscFunctionBeginUser;
    PetscCheck(g_ref_path_set && g_ref_token_set, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE,
               "two-rank State A decomp check requires matching -candidate_ref_path and -candidate_ref_token from a preceding one-rank run");
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    if (rank == 0) {
        char magic[64], token[128];
        FILE *fp = fopen(g_ref_path, "r");
        PetscCheck(fp != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
                   "State A one-rank decomp reference missing; run the Makefile target so -n 1 precedes -n 2");
        PetscCheck(fscanf(fp, "%63s %127s", magic, token) == 2, PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "bad reference header");
        PetscCheck(strcmp(magic, "PICURV_CANDIDATE_STATEA_REF_V2") == 0, PETSC_COMM_SELF, PETSC_ERR_FILE_READ,
                   "bad reference magic");
        PetscCheck(strcmp(token, g_ref_token) == 0, PETSC_COMM_SELF, PETSC_ERR_FILE_READ,
                   "State A one-rank decomp reference token mismatch");
        int active, cblock, ci, cj;
        PetscCheck(fscanf(fp, "%le %le %le", &vals[0], &vals[1], &vals[2]) == 3, PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "bad reference R0");
        PetscCheck(fscanf(fp, "%le %le %le", &vals[3], &vals[4], &vals[5]) == 3, PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "bad reference Jv");
        PetscCheck(fscanf(fp, "%le %le %le", &vals[6], &vals[7], &vals[8]) == 3, PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "bad reference Phi");
        PetscCheck(fscanf(fp, "%le %le %le %d %d %d %d", &vals[9], &vals[10], &vals[11],
                          &active, &cblock, &ci, &cj) == 7, PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "bad reference estimator");
        vals[12] = (PetscReal)active; vals[13] = (PetscReal)cblock; vals[14] = (PetscReal)ci; vals[15] = (PetscReal)cj;
        fclose(fp);
        PetscCheck(remove(g_ref_path) == 0, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
                   "could not remove consumed State A MPI reference '%s'", g_ref_path);
    }
    PetscCallMPI(MPI_Bcast(vals, 16, MPIU_REAL, 0, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertRealNear(r0.n2, vals[0], 1e-10, "MPI State A R0 L2"));
    PetscCall(PicurvAssertRealNear(r0.ninf, vals[1], 1e-10, "MPI State A R0 Linf"));
    PetscCall(PicurvAssertRealNear(r0.checksum, vals[2], 1e-10, "MPI State A R0 checksum"));
    PetscCall(PicurvAssertRealNear(jv.n2, vals[3], 1e-8, "MPI State A Jv L2"));
    PetscCall(PicurvAssertRealNear(jv.ninf, vals[4], 1e-8, "MPI State A Jv Linf"));
    PetscCall(PicurvAssertRealNear(jv.checksum, vals[5], 1e-8, "MPI State A Jv checksum"));
    PetscCall(PicurvAssertRealNear(phi.n2, vals[6], 1e-8, "MPI State A DPhi v L2"));
    PetscCall(PicurvAssertRealNear(phi.ninf, vals[7], 1e-8, "MPI State A DPhi v Linf"));
    PetscCall(PicurvAssertRealNear(phi.checksum, vals[8], 1e-8, "MPI State A DPhi v checksum"));
    PetscCall(PicurvAssertRealNear(rep->lambda_B - rep->lambda_t, vals[9], 1e-10, "MPI State A lambda_cB"));
    PetscCall(PicurvAssertRealNear(rep->lambda_C - rep->lambda_t, vals[10], 1e-10, "MPI State A lambda_cC"));
    PetscCall(PicurvAssertRealNear(rep->lambda_D - rep->lambda_t, vals[11], 1e-10, "MPI State A lambda_cD"));
    PetscCall(PicurvAssertBool((PetscBool)(rep->active_cells == (PetscInt)vals[12]), "MPI State A active cell count"));
    PetscCall(PicurvAssertBool((PetscBool)(rep->cblock == (PetscInt)vals[13]), "MPI State A controlling block"));
    PetscFunctionReturn(0);
}

/**
 * @brief Runs State A matrix-free residual, Jv, and four-stage MPI decomposition checks.
 */
static PetscErrorCode StateADecompDirectionalCheck(UserCtx *user, Vec Ubase, Vec Rhs,
                                                   PetscReal lcC, const MomStabilityReport *rep)
{
    DofMap map;
    Vec V, Up, Um, PhiP, PhiM, Ustage;
    PetscReal *v, *R0, *Rp, *Rm, *ActP, *ActM, *Out;
    GlobalVecStats sR0, sJv, sPhi;
    const PetscReal eps = 1e-6, dtau = 0.5/lcC;
    PetscMPIInt size;
    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    PetscCall(DofMapBuildOwned(user, &map));
    PetscCall(VecDuplicate(Ubase, &V)); PetscCall(VecDuplicate(Ubase, &Up));
    PetscCall(VecDuplicate(Ubase, &Um)); PetscCall(VecDuplicate(Ubase, &PhiP));
    PetscCall(VecDuplicate(Ubase, &PhiM)); PetscCall(VecDuplicate(Ubase, &Ustage));
    PetscCall(PetscMalloc6(map.n, &v, map.n, &R0, map.n, &Rp, map.n, &Rm, map.n, &ActP, map.n, &ActM));
    PetscCall(PetscMalloc1(map.n, &Out));
    PetscCall(FillDeterministicDirection(user, V, &map, v));

    PetscCall(EvalConvResidual(user, Ubase, Rhs, &map, R0));
    PetscCall(ActiveStats(&map, R0, &sR0));
    PetscCall(VecCopy(Ubase, Up)); PetscCall(VecAXPY(Up, eps, V));
    PetscCall(VecCopy(Ubase, Um)); PetscCall(VecAXPY(Um, -eps, V));
    PetscCall(EvalConvResidual(user, Up, Rhs, &map, Rp));
    PetscCall(EvalConvResidual(user, Um, Rhs, &map, Rm));
    for (PetscInt m = 0; m < map.n; m++) Out[m] = (Rp[m]-Rm[m])/(2.0*eps);
    PetscCall(ActiveStats(&map, Out, &sJv));

    PetscCall(FourStage(user, Up, dtau, Rhs, &map, Rp, Ustage, PhiP));
    PetscCall(FourStage(user, Um, dtau, Rhs, &map, Rm, Ustage, PhiM));
    PetscCall(ExtractActiveVector(user, PhiP, &map, ActP));
    PetscCall(ExtractActiveVector(user, PhiM, &map, ActM));
    for (PetscInt m = 0; m < map.n; m++) Out[m] = (ActP[m]-ActM[m])/(2.0*eps);
    PetscCall(ActiveStats(&map, Out, &sPhi));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "  State A matrix-free MPI check: R0 L2=%.6e Linf=%.6e checksum=%.6e\n"
        "                                Jv L2=%.6e Linf=%.6e checksum=%.6e\n"
        "                           DPhi v L2=%.6e Linf=%.6e checksum=%.6e\n",
        (double)sR0.n2, (double)sR0.ninf, (double)sR0.checksum,
        (double)sJv.n2, (double)sJv.ninf, (double)sJv.checksum,
        (double)sPhi.n2, (double)sPhi.ninf, (double)sPhi.checksum));
    if (size == 1) PetscCall(WriteStateADecompReference(sR0, sJv, sPhi, rep));
    else           PetscCall(ReadAndCompareStateADecompReference(sR0, sJv, sPhi, rep));

    PetscCall(PetscFree6(v, R0, Rp, Rm, ActP, ActM)); PetscCall(PetscFree(Out));
    PetscCall(VecDestroy(&V)); PetscCall(VecDestroy(&Up)); PetscCall(VecDestroy(&Um));
    PetscCall(VecDestroy(&PhiP)); PetscCall(VecDestroy(&PhiM)); PetscCall(VecDestroy(&Ustage));
    PetscCall(DofMapDestroy(&map));
    PetscFunctionReturn(0);
}

/**
 * @brief Runs one decomp baseline and compares scalar estimates with regenerated references.
 */
static PetscErrorCode RunDecompBaseline(CandState st, const char *name,
                                        PetscReal expB, PetscReal expC, PetscReal expD)
{
    SimCtx *simCtx = NULL; UserCtx *user = NULL;
    Vec Ubase, Rhs;
    PetscReal repeat_err, maxdiv;
    SeamDiagnostics seam;
    MomStabilityReport rep;
    const PetscInt N = 8;
    PetscFunctionBeginUser;

    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, N, N, N, PETSC_TRUE, PETSC_TRUE, PETSC_TRUE));
    PetscCall(ConfigureCandidateFixture(simCtx, user));
    PetscCall(VecDuplicate(user->Ucont, &Ubase));
    PetscCall(VecDuplicate(user->Ucont, &Rhs));
    PetscCall(BuildBaseState(user, st, Ubase, &repeat_err, &maxdiv, &seam));
    PetscCall(ComputeMomentumStabilityEstimate(user, 1, simCtx->dt, MOM_STAB_CAND_C, &rep));

    const PetscReal lcB = rep.lambda_B - rep.lambda_t;
    const PetscReal lcC = rep.lambda_C - rep.lambda_t;
    const PetscReal lcD = rep.lambda_D - rep.lambda_t;
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "\n================ DECOMP BASELINE %s ================\n"
        "  active cells: %d | repeat=%.3e | max|div_h Ucont|=%.3e\n"
        "  lambda_cB=%.6e  lambda_cC=%.6e  lambda_cD=%.6e\n",
        name, (int)rep.active_cells, (double)repeat_err, (double)maxdiv,
        (double)lcB, (double)lcC, (double)lcD));

    PetscCall(PicurvAssertBool((PetscBool)(rep.active_cells == 343), "decomp: active-cell count invariant"));
    PetscCall(PicurvAssertRealNear(expB, lcB, 5e-7, "decomp: lambda_cB matches one-rank baseline"));
    PetscCall(PicurvAssertRealNear(expC, lcC, 5e-7, "decomp: lambda_cC matches one-rank baseline"));
    PetscCall(PicurvAssertRealNear(expD, lcD, 5e-7, "decomp: lambda_cD matches one-rank baseline"));
    if (st == STATE_A) {
        PetscCall(StateADecompDirectionalCheck(user, Ubase, Rhs, lcC, &rep));
        PetscCall(PicurvAssertRealNear(maxdiv, 0.0, 1e-12, "decomp A: divergence-free"));
        PetscCall(PicurvAssertRealNear(lcB, lcC, 1e-12, "decomp A: B == C"));
        PetscCall(PicurvAssertRealNear(lcC, lcD, 1e-12, "decomp A: C == D"));
    } else if (st == STATE_B) {
        PetscCall(PicurvAssertBool((PetscBool)(maxdiv > 1e-3), "decomp B: nonzero divergence"));
        PetscCall(PicurvAssertBool((PetscBool)(lcC > lcB), "decomp B: C > B"));
    } else {
        PetscCall(PicurvAssertRealNear(maxdiv, 0.0, 1e-12, "decomp C: divergence-free"));
        PetscCall(PicurvAssertRealNear(lcB, lcC, 1e-12, "decomp C: B == C"));
    }

    PetscCall(VecDestroy(&Ubase));
    PetscCall(VecDestroy(&Rhs));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Runs the State A decomposition baseline.
 */
static PetscErrorCode TestDecompA(void) { return RunDecompBaseline(STATE_A, "A (uniform div-free)", 1.1000000, 1.1000000, 1.1000000); }

/**
 * @brief Runs the State B decomposition baseline.
 */
static PetscErrorCode TestDecompB(void) { return RunDecompBaseline(STATE_B, "B (nonzero divergence)", 0.878379697325, 0.974927912182, 1.572173303885); }

/**
 * @brief Runs the State C decomposition baseline.
 */
static PetscErrorCode TestDecompC(void) { return RunDecompBaseline(STATE_C, "C (div-free shear)", 1.5000000, 1.5000000, 1.780330085890); }

/**
 * @brief Reads optional paired-run MPI reference path and token.
 */
static PetscErrorCode ConfigureMPIReferenceOptions(void)
{
    PetscFunctionBeginUser;
    PetscCall(PetscOptionsGetString(NULL, NULL, "-candidate_ref_path",
                                    g_ref_path, sizeof(g_ref_path), &g_ref_path_set));
    PetscCall(PetscOptionsGetString(NULL, NULL, "-candidate_ref_token",
                                    g_ref_token, sizeof(g_ref_token), &g_ref_token_set));
    if (g_ref_path_set || g_ref_token_set) {
        PetscCheck(g_ref_path_set && g_ref_token_set, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE,
                   "-candidate_ref_path and -candidate_ref_token must be supplied together");
    }
    PetscFunctionReturn(0);
}

/**
 * @brief PETSc entry point for the focused convective-candidate harness.
 */
int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    PetscMPIInt size;
    const PicurvTestCase cases[] = {
        {"candidate-state-A-uniform", TestStateA},
        {"candidate-state-B-divergence", TestStateB},
        {"candidate-state-C-shear", TestStateC},
        {"candidate-state-A-grid-audit", TestStateAGridAudit},
    };
    const PicurvTestCase decomp_cases[] = {
        {"candidate-decomp-A-uniform", TestDecompA},
        {"candidate-decomp-B-divergence", TestDecompB},
        {"candidate-decomp-C-shear", TestDecompC},
    };
    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv A4a convective-candidate study");
    if (ierr) return (int)ierr;
    ierr = ConfigureMPIReferenceOptions(); if (ierr) { PetscFinalize(); return (int)ierr; }
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); if (ierr) { PetscFinalize(); return (int)ierr; }
    if (size == 1) {
        ierr = PicurvRunTests("unit-momentum-candidates", cases, sizeof(cases)/sizeof(cases[0]));
        if (!ierr) ierr = PicurvRunTests("unit-momentum-candidates-decomp", decomp_cases, sizeof(decomp_cases)/sizeof(decomp_cases[0]));
    } else {
        ierr = PicurvRunTests("unit-momentum-candidates-decomp", decomp_cases, sizeof(decomp_cases)/sizeof(decomp_cases[0]));
    }
    if (ierr) { PetscFinalize(); return (int)ierr; }
    return (int)PetscFinalize();
}
