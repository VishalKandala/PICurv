/**
 * @file test_statistics_accumulator.c
 * @brief C unit tests for per-window accumulator storage and pointwise application.
 *
 * The field series reuses the values verified analytically by the moment-kernel
 * suite, so a failure here isolates the field-to-accumulator wiring rather than
 * the numerics: three samples of (1,2,3), (3,6,5), (5,4,7) have means (3,4,5) and
 * symmetric centered products (8,4,8,8,4,8) in (xx,xy,xz,yy,yz,zz) order.
 */

#include "test_support.h"

#include "statistics_accumulator.h"
#include "statistics_moments.h"
#include "field_catalog.h"
#include "statistics_window.h"
#include "statistics_target.h"

#define ACC_N 6                 /* fixture size; cell domain is 5x5x5 */
#define ACC_CELLS 125

/** @brief Builds a definition requesting Ucat and P with second moments. */
static PicurvWindowDefinition AccDefinition(PetscBool want_second)
{
    PicurvWindowDefinition d;
    memset(&d, 0, sizeof(d));
    strncpy(d.name, "acc", PICURV_WINDOW_NAME_LENGTH - 1);
    d.weighting = PICURV_WEIGHTING_SAMPLE;
    d.cadence_kind = PICURV_CADENCE_STEP;
    d.step_cadence = 1;
    d.field_count = 2;
    d.fields[0].field_id = FIELD_ID_UCAT; d.fields[0].want_second = want_second;
    d.fields[1].field_id = FIELD_ID_P;    d.fields[1].want_second = want_second;
    return d;
}

/** @brief Sets Ucat to a uniform vector and P to a uniform scalar. */
static PetscErrorCode SetUniform(UserCtx *user, PetscReal x, PetscReal y, PetscReal z, PetscReal p)
{
    Cmpnts ***ucat = NULL;
    PetscReal ***pp = NULL;
    const DMDALocalInfo info = user->info;

    PetscFunctionBeginUser;
    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
    PetscCall(DMDAVecGetArray(user->da, user->P, &pp));
    for (PetscInt k = info.zs; k < info.zs + info.zm; ++k)
        for (PetscInt j = info.ys; j < info.ys + info.ym; ++j)
            for (PetscInt i = info.xs; i < info.xs + info.xm; ++i) {
                ucat[k][j][i].x = x; ucat[k][j][i].y = y; ucat[k][j][i].z = z;
                pp[k][j][i] = p;
            }
    PetscCall(DMDAVecRestoreArray(user->da, user->P, &pp));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));
    PetscFunctionReturn(0);
}

/** @brief Reads one interior point of a scalar accumulator vector. */
static PetscErrorCode ReadScalarAt(UserCtx *user, Vec v, PetscInt i, PetscInt j, PetscInt k, PetscReal *out)
{
    PetscReal ***a = NULL;
    PetscFunctionBeginUser;
    PetscCall(DMDAVecGetArrayRead(user->da, v, &a));
    *out = a[k][j][i];
    PetscCall(DMDAVecRestoreArrayRead(user->da, v, &a));
    PetscFunctionReturn(0);
}

/** @brief Reads one component of a multi-component accumulator vector at an interior point. */
static PetscErrorCode ReadComponentAt(UserCtx *user, Vec v, PetscInt components, PetscInt component,
                                      PetscInt i, PetscInt j, PetscInt k, PetscReal *out)
{
    PetscScalar ****a = NULL;
    DM dm = NULL;
    PetscFunctionBeginUser;
    PetscCall(PicurvStatisticsComponentDM(user, components, &dm));
    PetscCall(DMDAVecGetArrayDOFRead(dm, v, &a));
    *out = a[k][j][i][component];
    PetscCall(DMDAVecRestoreArrayDOFRead(dm, v, &a));
    PetscFunctionReturn(0);
}

/** @brief Component counts follow the documented product shapes. */
static PetscErrorCode TestComponentCounts(void)
{
    PetscInt n = 0;
    PetscErrorCode bad = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvProductComponentCount(1, &n));
    PetscCall(PicurvAssertIntEqual(1, n, "a scalar self-product has one component"));
    PetscCall(PicurvProductComponentCount(3, &n));
    PetscCall(PicurvAssertIntEqual(6, n, "a three-vector self-product has six symmetric components"));
    PetscCall(PicurvCovarianceComponentCount(3, 1, &n));
    PetscCall(PicurvAssertIntEqual(3, n, "a vector-scalar covariance has three components"));
    PetscCall(PicurvCovarianceComponentCount(1, 1, &n));
    PetscCall(PicurvAssertIntEqual(1, n, "a scalar-scalar covariance has one component"));

    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    bad = PicurvCovarianceComponentCount(3, 3, &n);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_OUTOFRANGE, bad,
                                   "vector-vector cross products are an explicit non-goal"));
    PetscFunctionReturn(0);
}

/** @brief Storage allocates one mean per field and the right number of product components. */
static PetscErrorCode TestStorageShape(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindowStorage storage;
    PicurvWindowDefinition d = AccDefinition(PETSC_TRUE);
    PetscInt bs = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, ACC_N, ACC_N, ACC_N));
    PetscCall(PicurvWindowStorageCreate(user, &d, &storage));

    PetscCall(PicurvAssertIntEqual(2, storage.field_count, "two fields requested"));

    /* Each product is one vector carrying all of its components, so the symmetric
     * second-order tensor stays a single object rather than six scalars. */
    PetscCall(VecGetBlockSize(storage.m2[0], &bs));
    PetscCall(PicurvAssertIntEqual(6, bs, "the Ucat product is one six-component tensor"));
    PetscCall(VecGetBlockSize(storage.m2[1], &bs));
    PetscCall(PicurvAssertIntEqual(1, bs, "the P product is a single scalar"));

    /* Means inherit their source field's layout. */
    PetscCall(VecGetBlockSize(storage.mean[0], &bs));
    PetscCall(PicurvAssertIntEqual(3, bs, "the Ucat mean is a three-vector"));
    PetscCall(VecGetBlockSize(storage.mean[1], &bs));
    PetscCall(PicurvAssertIntEqual(1, bs, "the P mean is a scalar"));

    PetscCall(PicurvWindowStorageDestroy(&storage));
    PetscCall(PicurvAssertIntEqual(0, storage.field_count, "destroy zeroes the storage"));
    PetscCall(PicurvAssertBool((PetscBool)(storage.mean == NULL), "destroy releases the mean array"));
    /* Destroying zeroed storage is safe. */
    PetscCall(PicurvWindowStorageDestroy(&storage));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief Three accepted states reproduce the analytically known moments at every point. */
static PetscErrorCode TestKnownAccumulationAcrossField(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindowStorage storage;
    PicurvWindowDefinition d = AccDefinition(PETSC_TRUE);
    const PetscReal series[3][3] = {{1.0, 2.0, 3.0}, {3.0, 6.0, 5.0}, {5.0, 4.0, 7.0}};
    const PetscReal scalars[3] = {1.0, 2.0, 6.0};
    const PetscReal expected_product[6] = {8.0, 4.0, 8.0, 8.0, 4.0, 8.0};
    const char *labels[6] = {"xx", "xy", "xz", "yy", "yz", "zz"};
    Cmpnts ***mean_vec = NULL;
    PetscReal value = 0.0;
    char context[128];

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, ACC_N, ACC_N, ACC_N));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(PicurvWindowStorageCreate(user, &d, &storage));

    for (PetscInt s = 0; s < 3; ++s) {
        PetscCall(SetUniform(user, series[s][0], series[s][1], series[s][2], scalars[s]));
        PetscCall(PicurvWindowAccumulate(user, &d, &storage, 1.0));
    }

    /* Occupancy: every fluid point saw all three states. */
    PetscCall(ReadScalarAt(user, storage.count, 2, 2, 2, &value));
    PetscCall(PicurvAssertRealNear(3.0, value, 1.0e-12, "per-point sample count"));
    PetscCall(ReadScalarAt(user, storage.weight, 2, 2, 2, &value));
    PetscCall(PicurvAssertRealNear(3.0, value, 1.0e-12, "per-point total weight"));
    PetscCall(ReadScalarAt(user, storage.weight_sq, 2, 2, 2, &value));
    PetscCall(PicurvAssertRealNear(3.0, value, 1.0e-12, "per-point squared-weight sum"));

    /* Vector mean. */
    PetscCall(DMDAVecGetArrayRead(user->fda, storage.mean[0], &mean_vec));
    PetscCall(PicurvAssertRealNear(3.0, mean_vec[2][2][2].x, 1.0e-12, "Ucat mean x"));
    PetscCall(PicurvAssertRealNear(4.0, mean_vec[2][2][2].y, 1.0e-12, "Ucat mean y"));
    PetscCall(PicurvAssertRealNear(5.0, mean_vec[2][2][2].z, 1.0e-12, "Ucat mean z"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, storage.mean[0], &mean_vec));

    /* All six symmetric centered products, in fixed catalog order. */
    for (PetscInt c = 0; c < 6; ++c) {
        PetscCall(ReadComponentAt(user, storage.m2[0], 6, c, 2, 2, 2, &value));
        PetscCall(PetscSNPrintf(context, sizeof(context), "Ucat centered product %s", labels[c]));
        PetscCall(PicurvAssertRealNear(expected_product[c], value, 1.0e-11, context));
    }

    /* Scalar mean and second moment. */
    PetscCall(ReadScalarAt(user, storage.mean[1], 2, 2, 2, &value));
    PetscCall(PicurvAssertRealNear(3.0, value, 1.0e-12, "P mean"));
    PetscCall(ReadScalarAt(user, storage.m2[1], 2, 2, 2, &value));
    PetscCall(PicurvAssertRealNear(14.0, value, 1.0e-11, "P centered second moment"));

    PetscCall(PicurvWindowStorageDestroy(&storage));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief The valid-fraction range reports mask coverage, including a never-valid point.
 *
 * This is the mask-health indicator the console snapshot prints, and it is the only
 * signal that distinguishes a mean built from every state from one built from none.
 * A blanked point keeps a zero mean that looks like a legitimate value, so the range
 * has to surface it.
 */
static PetscErrorCode TestValidFractionRange(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindowStorage storage;
    PicurvWindowDefinition d = AccDefinition(PETSC_FALSE);
    PetscReal ***nvert = NULL;
    PetscReal lowest = 0.0, highest = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, ACC_N, ACC_N, ACC_N));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(PicurvWindowStorageCreate(user, &d, &storage));
    PetscCall(SetUniform(user, 1.0, 2.0, 3.0, 7.0));

    /* Before any sample the range is the documented degenerate, not a division. */
    PetscCall(PicurvWindowValidFractionRange(user, &d, &storage, 0, &lowest, &highest));
    PetscCall(PicurvAssertRealNear(1.0, lowest, 1.0e-12, "an unsampled window reports full coverage"));
    PetscCall(PicurvAssertRealNear(0.0, highest, 1.0e-12, "an unsampled window reports no maximum"));

    /* An unobstructed domain: every targeted point saw every state. */
    PetscCall(PicurvWindowAccumulate(user, &d, &storage, 1.0));
    PetscCall(PicurvWindowAccumulate(user, &d, &storage, 1.0));
    PetscCall(PicurvWindowValidFractionRange(user, &d, &storage, 2, &lowest, &highest));
    PetscCall(PicurvAssertRealNear(1.0, lowest, 1.0e-12, "a clear domain is fully covered"));
    PetscCall(PicurvAssertRealNear(1.0, highest, 1.0e-12, "a clear domain has no over-counted point"));

    /* Blank one interior cell and take a third state: that point now trails. */
    PetscCall(DMDAVecGetArray(user->da, user->Nvert, &nvert));
    nvert[3][3][3] = 1.0;
    PetscCall(DMDAVecRestoreArray(user->da, user->Nvert, &nvert));
    PetscCall(PicurvWindowAccumulate(user, &d, &storage, 1.0));
    PetscCall(PicurvWindowValidFractionRange(user, &d, &storage, 3, &lowest, &highest));
    PetscCall(PicurvAssertRealNear(2.0 / 3.0, lowest, 1.0e-12,
                                   "a point blanked for one state trails the rest"));
    PetscCall(PicurvAssertRealNear(1.0, highest, 1.0e-12, "unobstructed points stay fully covered"));

    PetscCall(PicurvWindowStorageDestroy(&storage));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief A point blanked from the start is reported as never valid.
 */
static PetscErrorCode TestValidFractionDetectsNeverValidPoint(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindowStorage storage;
    PicurvWindowDefinition d = AccDefinition(PETSC_FALSE);
    PetscReal ***nvert = NULL;
    PetscReal lowest = 1.0, highest = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, ACC_N, ACC_N, ACC_N));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(DMDAVecGetArray(user->da, user->Nvert, &nvert));
    nvert[3][3][3] = 1.0;
    PetscCall(DMDAVecRestoreArray(user->da, user->Nvert, &nvert));

    PetscCall(PicurvWindowStorageCreate(user, &d, &storage));
    PetscCall(SetUniform(user, 1.0, 2.0, 3.0, 7.0));
    PetscCall(PicurvWindowAccumulate(user, &d, &storage, 1.0));
    PetscCall(PicurvWindowAccumulate(user, &d, &storage, 1.0));

    PetscCall(PicurvWindowValidFractionRange(user, &d, &storage, 2, &lowest, &highest));
    PetscCall(PicurvAssertRealNear(0.0, lowest, 1.0e-12,
                                   "a permanently blanked point reports zero coverage"));
    PetscCall(PicurvAssertRealNear(1.0, highest, 1.0e-12, "the rest of the domain is unaffected"));

    PetscCall(PicurvWindowStorageDestroy(&storage));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief A vector-scalar covariance reproduces its analytically known components. */
static PetscErrorCode TestCovarianceAccumulation(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindowStorage storage;
    PicurvWindowDefinition d = AccDefinition(PETSC_TRUE);
    const PetscReal series[3][3] = {{1.0, 2.0, 3.0}, {3.0, 6.0, 5.0}, {5.0, 4.0, 7.0}};
    const PetscReal scalars[3] = {1.0, 2.0, 6.0};
    /* Means are (3,4,5) and 3, so the centered sums are
     * x: (-2)(-2)+(0)(-1)+(2)(3) = 10, y: (-2)(-2)+(2)(-1)+(0)(3) = 2,
     * z: (-2)(-2)+(0)(-1)+(2)(3) = 10. */
    const PetscReal expected[3] = {10.0, 2.0, 10.0};
    const char *labels[3] = {"x", "y", "z"};
    PetscReal value = 0.0;
    char context[128];

    PetscFunctionBeginUser;
    d.covariance_count = 1;
    d.covariances[0].first = FIELD_ID_UCAT;
    d.covariances[0].second = FIELD_ID_P;

    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, ACC_N, ACC_N, ACC_N));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(PicurvWindowStorageCreate(user, &d, &storage));
    {
        PetscInt bs = 0;
        PetscCall(VecGetBlockSize(storage.cm[0], &bs));
        PetscCall(PicurvAssertIntEqual(3, bs, "a vector-scalar pair is one three-component vector"));
    }

    for (PetscInt s = 0; s < 3; ++s) {
        PetscCall(SetUniform(user, series[s][0], series[s][1], series[s][2], scalars[s]));
        PetscCall(PicurvWindowAccumulate(user, &d, &storage, 1.0));
    }

    for (PetscInt c = 0; c < 3; ++c) {
        PetscCall(ReadComponentAt(user, storage.cm[0], 3, c, 2, 2, 2, &value));
        PetscCall(PetscSNPrintf(context, sizeof(context), "Ucat-P co-moment %s", labels[c]));
        PetscCall(PicurvAssertRealNear(expected[c], value, 1.0e-11, context));
    }

    /* The co-moment pass must not disturb the moments computed alongside it. */
    PetscCall(ReadScalarAt(user, storage.mean[1], 2, 2, 2, &value));
    PetscCall(PicurvAssertRealNear(3.0, value, 1.0e-12, "P mean is unaffected by the covariance pass"));
    PetscCall(ReadComponentAt(user, storage.m2[0], 6, 0, 2, 2, 2, &value));
    PetscCall(PicurvAssertRealNear(8.0, value, 1.0e-11, "Ucat xx product is unaffected by the covariance pass"));

    PetscCall(PicurvWindowStorageDestroy(&storage));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief A field's covariance with itself reduces exactly to its own second moment. */
static PetscErrorCode TestScalarSelfCovarianceMatchesSecondMoment(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindowStorage storage;
    PicurvWindowDefinition d = AccDefinition(PETSC_TRUE);
    const PetscReal scalars[3] = {1.0, 2.0, 6.0};
    PetscReal co_moment = 0.0;
    PetscReal second_moment = 0.0;

    PetscFunctionBeginUser;
    /* Pairing P with itself makes the co-moment pass and the moment pass compute the
     * same quantity by different routes, so any drift between them shows up here
     * rather than only in a cross-field result no closed form covers. */
    d.covariance_count = 1;
    d.covariances[0].first = FIELD_ID_P;
    d.covariances[0].second = FIELD_ID_P;

    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, ACC_N, ACC_N, ACC_N));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(PicurvWindowStorageCreate(user, &d, &storage));

    for (PetscInt s = 0; s < 3; ++s) {
        PetscCall(SetUniform(user, 0.0, 0.0, 0.0, scalars[s]));
        PetscCall(PicurvWindowAccumulate(user, &d, &storage, 1.0));
    }

    PetscCall(ReadScalarAt(user, storage.cm[0], 2, 2, 2, &co_moment));
    PetscCall(ReadScalarAt(user, storage.m2[1], 2, 2, 2, &second_moment));
    PetscCall(PicurvAssertRealNear(14.0, co_moment, 1.0e-11, "P self-covariance equals its centered sum"));
    PetscCall(PicurvAssertBool((PetscBool)(co_moment == second_moment),
                               "self-covariance and second moment agree bit for bit"));

    PetscCall(PicurvWindowStorageDestroy(&storage));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief A covariance member missing from the field list is rejected, not silently skipped. */
static PetscErrorCode TestCovarianceRequiresFieldMembership(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindowStorage storage;
    PicurvWindowDefinition d = AccDefinition(PETSC_FALSE);
    PetscErrorCode bad = 0;

    PetscFunctionBeginUser;
    /* Nvert is catalogued but is not in the window's field list, so no running mean
     * exists for it and the co-moment update has nothing to center against. */
    d.covariance_count = 1;
    d.covariances[0].first = FIELD_ID_NVERT;
    d.covariances[0].second = FIELD_ID_P;

    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, ACC_N, ACC_N, ACC_N));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(PicurvWindowStorageCreate(user, &d, &storage));
    PetscCall(SetUniform(user, 1.0, 2.0, 3.0, 4.0));

    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    bad = PicurvWindowAccumulate(user, &d, &storage, 1.0);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_WRONGSTATE, bad,
                                   "a covariance member outside the field list is rejected"));

    PetscCall(PicurvWindowStorageDestroy(&storage));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Derived quantities reproduce the analytic values the accumulated state implies.
 *
 * The series is the one the moment suite verifies by hand, so a failure here isolates
 * the normalization and component selection rather than the accumulation: with three
 * unit-weight samples, `R_ij = C_ij / 3`, and the diagonal of the Ucat tensor is 8/3
 * in each direction, giving TKE = 4 and an RMS of sqrt(8/3) per component.
 */
static PetscErrorCode TestDerivedQuantities(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindowStorage storage;
    PicurvWindowDefinition d = AccDefinition(PETSC_TRUE);
    const PetscReal series[3][3] = {{1.0, 2.0, 3.0}, {3.0, 6.0, 5.0}, {5.0, 4.0, 7.0}};
    const PetscReal scalars[3] = {1.0, 2.0, 6.0};
    Vec scalar_target = NULL, vector_target = NULL;
    PetscInt count = 0;
    PetscReal value = 0.0;

    PetscFunctionBeginUser;
    d.covariance_count = 1;
    d.covariances[0].first = FIELD_ID_UCAT;
    d.covariances[0].second = FIELD_ID_P;

    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, ACC_N, ACC_N, ACC_N));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(PicurvWindowStorageCreate(user, &d, &storage));
    for (PetscInt s = 0; s < 3; ++s) {
        PetscCall(SetUniform(user, series[s][0], series[s][1], series[s][2], scalars[s]));
        PetscCall(PicurvWindowAccumulate(user, &d, &storage, 1.0));
    }
    PetscCall(DMCreateGlobalVector(user->da, &scalar_target));
    PetscCall(DMCreateGlobalVector(user->fda, &vector_target));

    /* Requesting every output resolves each kind against what the window carries. */
    PetscCall(PicurvWindowDerivedCount(&d, &storage, "mean,reynolds_stress,rms,tke,flux", &count));
    /* 2 means + (6 Ucat + 1 P) stresses + (3 Ucat + 1 P) RMS + 1 TKE + 1 flux. */
    PetscCall(PicurvAssertIntEqual(15, count, "every requested output is enumerated"));

    for (PetscInt index = 0; index < count; ++index) {
        PicurvDerivedField field;
        PetscBool is_tke = PETSC_FALSE, is_rxx = PETSC_FALSE, is_rmsx = PETSC_FALSE;
        PetscBool is_flux = PETSC_FALSE, is_pvar = PETSC_FALSE;

        PetscCall(PicurvWindowDerive(user, &d, &storage, "mean,reynolds_stress,rms,tke,flux",
                                     index, scalar_target, vector_target, &field));
        PetscCall(PetscStrcmp(field.name, "acc_Ucat_tke", &is_tke));
        PetscCall(PetscStrcmp(field.name, "acc_Ucat_R_xx", &is_rxx));
        PetscCall(PetscStrcmp(field.name, "acc_Ucat_rmsx", &is_rmsx));
        PetscCall(PetscStrcmp(field.name, "acc_Ucat_P_flux", &is_flux));
        PetscCall(PetscStrcmp(field.name, "acc_P_variance", &is_pvar));

        if (is_rxx) {
            PetscCall(ReadScalarAt(user, scalar_target, 2, 2, 2, &value));
            PetscCall(PicurvAssertRealNear(8.0 / 3.0, value, 1.0e-11, "R_xx is C_xx over the weight"));
        } else if (is_rmsx) {
            PetscCall(ReadScalarAt(user, scalar_target, 2, 2, 2, &value));
            PetscCall(PicurvAssertRealNear(PetscSqrtReal(8.0 / 3.0), value, 1.0e-11,
                                           "RMS is the root of the normal stress"));
        } else if (is_tke) {
            PetscCall(ReadScalarAt(user, scalar_target, 2, 2, 2, &value));
            PetscCall(PicurvAssertRealNear(4.0, value, 1.0e-11, "TKE is half the trace"));
        } else if (is_pvar) {
            PetscCall(ReadScalarAt(user, scalar_target, 2, 2, 2, &value));
            PetscCall(PicurvAssertRealNear(14.0 / 3.0, value, 1.0e-11, "a scalar product normalizes too"));
        } else if (is_flux) {
            Cmpnts ***flux = NULL;

            PetscCall(DMDAVecGetArrayRead(user->fda, vector_target, &flux));
            PetscCall(PicurvAssertRealNear(10.0 / 3.0, flux[2][2][2].x, 1.0e-11, "flux x"));
            PetscCall(PicurvAssertRealNear(2.0 / 3.0, flux[2][2][2].y, 1.0e-11, "flux y"));
            PetscCall(PicurvAssertRealNear(10.0 / 3.0, flux[2][2][2].z, 1.0e-11, "flux z"));
            PetscCall(DMDAVecRestoreArrayRead(user->fda, vector_target, &flux));
        }
    }

    /* A narrower request produces only what it asked for. */
    PetscCall(PicurvWindowDerivedCount(&d, &storage, "tke", &count));
    PetscCall(PicurvAssertIntEqual(1, count, "a single output kind enumerates once"));

    /* An unsampled point is left at zero rather than divided by a zero weight. */
    {
        PicurvDerivedField field;
        PetscReal ***nvert = NULL;

        PetscCall(DMDAVecGetArray(user->da, user->Nvert, &nvert));
        nvert[3][3][3] = 1.0;
        PetscCall(DMDAVecRestoreArray(user->da, user->Nvert, &nvert));
        PetscCall(PicurvWindowStorageDestroy(&storage));
        PetscCall(PicurvWindowStorageCreate(user, &d, &storage));
        PetscCall(SetUniform(user, 1.0, 2.0, 3.0, 4.0));
        PetscCall(PicurvWindowAccumulate(user, &d, &storage, 1.0));
        PetscCall(PicurvWindowDerive(user, &d, &storage, "tke", 0, scalar_target, vector_target, &field));
        PetscCall(ReadScalarAt(user, scalar_target, 3, 3, 3, &value));
        PetscCall(PicurvAssertRealNear(0.0, value, 1.0e-12,
                                       "a never-sampled point derives to zero, not a division"));
    }

    /* Dimensionalization raises the field's own scale to the kind's own power: a mean
     * is linear in it and a covariance is quadratic. One blanket velocity factor would
     * be right for two of the five kinds and wrong for the other three, which is why
     * derived statistics were left non-dimensional before this existed. */
    {
        PicurvDerivedField field;
        PetscReal nondimensional_mean = 0.0, nondimensional_tke = 0.0;
        PetscReal dimensional_mean = 0.0, dimensional_tke = 0.0;
        const PetscReal velocity = 3.0;

        PetscCall(PicurvWindowStorageDestroy(&storage));
        PetscCall(PicurvWindowStorageCreate(user, &d, &storage));
        PetscCall(SetUniform(user, 1.0, 2.0, 3.0, 4.0));
        PetscCall(PicurvWindowAccumulate(user, &d, &storage, 1.0));
        PetscCall(SetUniform(user, 3.0, 4.0, 5.0, 6.0));
        PetscCall(PicurvWindowAccumulate(user, &d, &storage, 1.0));

        /* The minimal fixture carries no post-processing settings; the harness teardown
         * already frees them, so allocating here is symmetric. */
        if (!user->simCtx->pps) PetscCall(PetscCalloc1(1, &user->simCtx->pps));
        user->simCtx->pps->dimensionalize = PETSC_FALSE;
        PetscCall(PicurvWindowDerive(user, &d, &storage, "mean", 0, scalar_target, vector_target, &field));
        PetscCall(ReadComponentAt(user, vector_target, 3, 0, 2, 2, 2, &nondimensional_mean));
        PetscCall(PicurvWindowDerive(user, &d, &storage, "tke", 0, scalar_target, vector_target, &field));
        PetscCall(ReadScalarAt(user, scalar_target, 2, 2, 2, &nondimensional_tke));

        user->simCtx->scaling.U_ref = velocity;
        user->simCtx->pps->dimensionalize = PETSC_TRUE;
        PetscCall(PicurvWindowDerive(user, &d, &storage, "mean", 0, scalar_target, vector_target, &field));
        PetscCall(ReadComponentAt(user, vector_target, 3, 0, 2, 2, 2, &dimensional_mean));
        PetscCall(PicurvWindowDerive(user, &d, &storage, "tke", 0, scalar_target, vector_target, &field));
        PetscCall(ReadScalarAt(user, scalar_target, 2, 2, 2, &dimensional_tke));

        PetscCall(PicurvAssertRealNear(velocity * nondimensional_mean, dimensional_mean, 1.0e-11,
                                       "a mean carries the field's scale once"));
        PetscCall(PicurvAssertRealNear(velocity * velocity * nondimensional_tke, dimensional_tke,
                                       1.0e-11, "a turbulent kinetic energy carries it squared"));
        user->simCtx->pps->dimensionalize = PETSC_FALSE;
        user->simCtx->scaling.U_ref = 1.0;
    }

    /* An unknown output kind is refused rather than silently ignored. */
    {
        PetscErrorCode bad = 0;

        PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
        bad = PicurvWindowDerivedCount(&d, &storage, "skewness", &count);
        PetscCall(PetscPopErrorHandler());
        PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_WRONG, bad, "an unknown output kind is refused"));
    }

    PetscCall(VecDestroy(&vector_target));
    PetscCall(VecDestroy(&scalar_target));
    PetscCall(PicurvWindowStorageDestroy(&storage));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief The spatial mean divides by sampled points, not by the whole vector.
 *
 * A derived field is zero outside the target domain and wherever the mask never
 * admitted a point. Those zeros are absences, and averaging over them silently
 * scales the answer down by the fraction of the vector the window never covered —
 * a wrong number that still looks plausible.
 */
static PetscErrorCode TestSpatialMeanExcludesUnsampledPoints(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindowStorage storage;
    PicurvWindowDefinition d = AccDefinition(PETSC_FALSE);
    SpatialTargetPlan plan;
    PetscReal ***nvert = NULL;
    PetscReal mean = 0.0;
    PetscReal whole_vector_mean = 0.0;
    PetscReal sum = 0.0;
    PetscInt vector_size = 0;
    PetscInt targeted = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, ACC_N, ACC_N, ACC_N));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(PicurvWindowStorageCreate(user, &d, &storage));
    PetscCall(SetUniform(user, 1.0, 2.0, 3.0, 4.0));
    PetscCall(PicurvWindowAccumulate(user, &d, &storage, 1.0));

    /* A field holding exactly one over every sampled point: its masked mean must be
     * one, whatever fraction of the vector the target domain occupies. */
    PetscCall(VecZeroEntries(user->PostScalar));
    PetscCall(SpatialTargetPlanCreate(user, FIELD_ID_P, PICURV_STATISTICS_MASK_FLUID, &plan));
    PetscCall(SpatialTargetPlanGlobalPointCount(&plan, PETSC_COMM_WORLD, &targeted));
    {
        PetscReal ***values = NULL;

        PetscCall(DMDAVecGetArray(user->da, user->PostScalar, &values));
        for (PetscInt k = plan.start[2]; k < plan.end[2]; ++k)
            for (PetscInt j = plan.start[1]; j < plan.end[1]; ++j)
                for (PetscInt i = plan.start[0]; i < plan.end[0]; ++i) values[k][j][i] = 1.0;
        PetscCall(DMDAVecRestoreArray(user->da, user->PostScalar, &values));
    }

    PetscCall(PicurvWindowSpatialMean(user, &d, &storage, user->PostScalar, &mean));
    PetscCall(PicurvAssertRealNear(1.0, mean, 1.0e-12,
                                   "the mean over sampled points ignores untargeted entries"));

    /* The whole-vector average is genuinely different, which is what makes the
     * distinction worth asserting rather than assuming. */
    PetscCall(VecSum(user->PostScalar, &sum));
    PetscCall(VecGetSize(user->PostScalar, &vector_size));
    whole_vector_mean = sum / (PetscReal)vector_size;
    PetscCall(PicurvAssertBool((PetscBool)(vector_size > targeted),
                               "the vector is larger than the targeted domain"));
    PetscCall(PicurvAssertBool((PetscBool)(whole_vector_mean < 0.9),
                               "a whole-vector average would understate the result"));

    /* Blanking a point removes it from both the sum and the count, so a uniform
     * field still averages to its own value rather than being diluted. */
    PetscCall(DMDAVecGetArray(user->da, user->Nvert, &nvert));
    nvert[3][3][3] = 1.0;
    PetscCall(DMDAVecRestoreArray(user->da, user->Nvert, &nvert));
    PetscCall(PicurvWindowStorageDestroy(&storage));
    PetscCall(PicurvWindowStorageCreate(user, &d, &storage));
    PetscCall(PicurvWindowAccumulate(user, &d, &storage, 1.0));
    PetscCall(PicurvWindowSpatialMean(user, &d, &storage, user->PostScalar, &mean));
    PetscCall(PicurvAssertRealNear(1.0, mean, 1.0e-12,
                                   "a never-sampled point is excluded from the average"));

    PetscCall(PicurvWindowStorageDestroy(&storage));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief Payload enumeration covers every vector exactly once with stable names. */
static PetscErrorCode TestPayloadEnumeration(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindowStorage storage;
    PicurvWindowDefinition d = AccDefinition(PETSC_TRUE);
    /* 3 occupancy + 2 means + 2 products + 1 co-moment. Each product and co-moment
     * is one payload carrying all of its components. */
    const char *expected[8] = {
        "count", "weight", "weight_sq",
        "Ucat_mean", "P_mean",
        "Ucat_m2", "P_m2",
        "Ucat_P_cm"
    };
    const PetscInt expected_components[8] = {1, 1, 1, 3, 1, 6, 1, 3};
    PetscInt count = 0;
    PetscErrorCode bad = 0;

    PetscFunctionBeginUser;
    d.covariance_count = 1;
    d.covariances[0].first = FIELD_ID_UCAT;
    d.covariances[0].second = FIELD_ID_P;

    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, ACC_N, ACC_N, ACC_N));
    PetscCall(PicurvWindowStorageCreate(user, &d, &storage));
    PetscCall(PicurvWindowStoragePayloadCount(&storage, &count));
    PetscCall(PicurvAssertIntEqual(8, count, "every accumulator vector is enumerated once"));

    for (PetscInt index = 0; index < count; ++index) {
        PicurvStatisticsPayload payload;
        PetscBool matches = PETSC_FALSE;
        PetscInt block_size = 0;
        char context[160];

        PetscCall(PicurvWindowStoragePayload(user, &d, &storage, index, &payload));
        PetscCall(PetscStrcmp(payload.name, expected[index], &matches));
        PetscCall(PetscSNPrintf(context, sizeof(context), "payload %d is named '%s', got '%s'",
                                (int)index, expected[index], payload.name));
        PetscCall(PicurvAssertBool(matches, context));
        PetscCall(PicurvAssertBool((PetscBool)(payload.vec != NULL), "an enumerated payload has a vector"));
        PetscCall(PicurvAssertBool((PetscBool)(payload.role != NULL), "an enumerated payload has a role"));

        /* The declared component count must match the vector actually enumerated,
         * because the manifest inventory records it and the reader trusts it. */
        PetscCall(VecGetBlockSize(payload.vec, &block_size));
        PetscCall(PetscSNPrintf(context, sizeof(context), "payload '%s' declares %d components",
                                payload.name, (int)expected_components[index]));
        PetscCall(PicurvAssertIntEqual(expected_components[index], payload.components, context));
        PetscCall(PicurvAssertIntEqual(payload.components, block_size,
                                       "the declared component count matches the vector"));
    }

    /* A field without a requested second moment contributes no product payload. */
    {
        PicurvWindowDefinition partial = AccDefinition(PETSC_FALSE);
        PicurvWindowStorage lean;

        partial.fields[0].want_second = PETSC_TRUE;   /* Ucat keeps its product, P does not. */
        PetscCall(PicurvWindowStorageCreate(user, &partial, &lean));
        PetscCall(PicurvWindowStoragePayloadCount(&lean, &count));
        PetscCall(PicurvAssertIntEqual(6, count, "a field without a second moment adds no payload"));
        PetscCall(PicurvWindowStorageDestroy(&lean));
    }

    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    bad = PicurvWindowStoragePayload(user, &d, &storage, count, NULL);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertBool((PetscBool)(bad != 0), "an out-of-range payload index is rejected"));

    PetscCall(PicurvWindowStorageDestroy(&storage));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief Unequal weights reproduce the weighted result the kernels define. */
static PetscErrorCode TestWeightedAccumulation(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindowStorage storage;
    PicurvWindowDefinition d = AccDefinition(PETSC_TRUE);
    PetscReal value = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, ACC_N, ACC_N, ACC_N));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(PicurvWindowStorageCreate(user, &d, &storage));

    /* P = 1 with weight 1, then P = 3 with weight 3: mean 2.5, M2 = 3. */
    PetscCall(SetUniform(user, 0.0, 0.0, 0.0, 1.0));
    PetscCall(PicurvWindowAccumulate(user, &d, &storage, 1.0));
    PetscCall(SetUniform(user, 0.0, 0.0, 0.0, 3.0));
    PetscCall(PicurvWindowAccumulate(user, &d, &storage, 3.0));

    PetscCall(ReadScalarAt(user, storage.weight, 2, 2, 2, &value));
    PetscCall(PicurvAssertRealNear(4.0, value, 1.0e-12, "weighted total weight"));
    PetscCall(ReadScalarAt(user, storage.weight_sq, 2, 2, 2, &value));
    PetscCall(PicurvAssertRealNear(10.0, value, 1.0e-12, "weighted squared-weight sum"));
    PetscCall(ReadScalarAt(user, storage.mean[1], 2, 2, 2, &value));
    PetscCall(PicurvAssertRealNear(2.5, value, 1.0e-12, "weighted mean"));
    PetscCall(ReadScalarAt(user, storage.m2[1], 2, 2, 2, &value));
    PetscCall(PicurvAssertRealNear(3.0, value, 1.0e-11, "weighted centered second moment"));

    PetscCall(PicurvWindowStorageDestroy(&storage));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief Masked points accumulate nothing, and remain distinguishable from unsampled ones. */
static PetscErrorCode TestMaskedPointsAreExcluded(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindowStorage storage;
    PicurvWindowDefinition d = AccDefinition(PETSC_FALSE);
    PetscReal ***nvert = NULL;
    PetscReal value = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, ACC_N, ACC_N, ACC_N));
    PetscCall(VecSet(user->Nvert, 0.0));
    /* Blank a single interior cell. */
    PetscCall(DMDAVecGetArray(user->da, user->Nvert, &nvert));
    nvert[3][3][3] = 1.0;
    PetscCall(DMDAVecRestoreArray(user->da, user->Nvert, &nvert));

    PetscCall(PicurvWindowStorageCreate(user, &d, &storage));
    PetscCall(SetUniform(user, 1.0, 2.0, 3.0, 7.0));
    PetscCall(PicurvWindowAccumulate(user, &d, &storage, 1.0));
    PetscCall(PicurvWindowAccumulate(user, &d, &storage, 1.0));

    PetscCall(ReadScalarAt(user, storage.count, 2, 2, 2, &value));
    PetscCall(PicurvAssertRealNear(2.0, value, 1.0e-12, "a fluid point sees every state"));
    PetscCall(ReadScalarAt(user, storage.count, 3, 3, 3, &value));
    PetscCall(PicurvAssertRealNear(0.0, value, 1.0e-12, "a blanked point accumulates no sample"));
    PetscCall(ReadScalarAt(user, storage.weight, 3, 3, 3, &value));
    PetscCall(PicurvAssertRealNear(0.0, value, 1.0e-12, "a blanked point accumulates no weight"));
    PetscCall(ReadScalarAt(user, storage.mean[1], 3, 3, 3, &value));
    PetscCall(PicurvAssertRealNear(0.0, value, 1.0e-12, "a blanked point keeps an untouched mean"));

    PetscCall(PicurvWindowStorageDestroy(&storage));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief A constant field yields exactly zero product at every point. */
static PetscErrorCode TestConstantFieldZeroProduct(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindowStorage storage;
    PicurvWindowDefinition d = AccDefinition(PETSC_TRUE);
    PetscReal value = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, ACC_N, ACC_N, ACC_N));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(PicurvWindowStorageCreate(user, &d, &storage));

    for (PetscInt s = 0; s < 5; ++s) {
        PetscCall(SetUniform(user, 2.5, -1.5, 4.0, 9.0));
        PetscCall(PicurvWindowAccumulate(user, &d, &storage, 1.0 + 0.5 * (PetscReal)s));
    }
    for (PetscInt c = 0; c < 6; ++c) {
        PetscCall(ReadComponentAt(user, storage.m2[0], 6, c, 2, 2, 2, &value));
        PetscCall(PicurvAssertBool((PetscBool)(value == 0.0),
                                   "a constant field gives bitwise-zero products"));
    }
    PetscCall(ReadScalarAt(user, storage.mean[1], 2, 2, 2, &value));
    PetscCall(PicurvAssertRealNear(9.0, value, 1.0e-12, "a constant field keeps its value as the mean"));

    PetscCall(PicurvWindowStorageDestroy(&storage));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief The runloop entry point applies exactly the states the schedule accepts.
 *
 * Drives FieldStatisticsUpdateWindows directly, so the scheduling decision and the
 * field accumulation are exercised through the same path the solver uses rather
 * than being tested only in isolation.
 */
static PetscErrorCode TestRunloopDriverAppliesScheduledStates(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindow window;
    PicurvWindowDefinition d = AccDefinition(PETSC_TRUE);
    PicurvWindowStorage storage;
    PetscReal value = 0.0;
    const PetscReal dt = 0.5;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, ACC_N, ACC_N, ACC_N));
    PetscCall(VecSet(user->Nvert, 0.0));

    /* Stride two, physical-time weighting, starting at the origin. */
    d.weighting = PICURV_WEIGHTING_PHYSICAL_TIME;
    d.step_cadence = 2;
    PetscCall(PicurvWindowInit(&window, &d));
    PetscCall(PicurvWindowStorageCreate(user, &d, &storage));

    simCtx->fieldStatisticsEnabled = PETSC_TRUE;
    simCtx->fieldStatisticsWindowCount = 1;
    simCtx->fieldStatisticsWindows = &window;
    user->fieldStatisticsStorage = &storage;

    /* Steps 0,2,4 are due. Step 0 anchors; steps 2 and 4 each represent 1.0 of time. */
    for (PetscInt step = 0; step <= 4; ++step) {
        PetscCall(SetUniform(user, 1.0, 2.0, 3.0, (PetscReal)(step + 1)));
        PetscCall(FieldStatisticsUpdateWindows(simCtx, step, (PetscReal)step * dt));
    }

    PetscCall(PicurvAssertIntEqual(2, window.sample_count, "stride two over five steps accepts two samples"));
    PetscCall(PicurvAssertRealNear(2.0, window.total_weight, 1.0e-12, "accepted weights span the elapsed time"));

    /* Only the due states reached the field accumulators. */
    PetscCall(ReadScalarAt(user, storage.count, 2, 2, 2, &value));
    PetscCall(PicurvAssertRealNear(2.0, value, 1.0e-12, "per-point count matches the accepted samples"));
    PetscCall(ReadScalarAt(user, storage.weight, 2, 2, 2, &value));
    PetscCall(PicurvAssertRealNear(2.0, value, 1.0e-12, "per-point weight matches the accepted weights"));

    /* P carried 3 at step 2 and 5 at step 4, each with weight 1: mean 4, M2 = 2. */
    PetscCall(ReadScalarAt(user, storage.mean[1], 2, 2, 2, &value));
    PetscCall(PicurvAssertRealNear(4.0, value, 1.0e-12, "only the scheduled samples contribute to the mean"));
    PetscCall(ReadScalarAt(user, storage.m2[1], 2, 2, 2, &value));
    PetscCall(PicurvAssertRealNear(2.0, value, 1.0e-11, "second moment reflects only the scheduled samples"));

    /* A disabled subsystem is inert on the same path. */
    simCtx->fieldStatisticsEnabled = PETSC_FALSE;
    PetscCall(FieldStatisticsUpdateWindows(simCtx, 6, 3.0));
    PetscCall(PicurvAssertIntEqual(2, window.sample_count, "a disabled subsystem accepts nothing"));

    simCtx->fieldStatisticsWindows = NULL;
    simCtx->fieldStatisticsWindowCount = 0;
    user->fieldStatisticsStorage = NULL;
    PetscCall(PicurvWindowStorageDestroy(&storage));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Entry point for the accumulator suite.
 */
int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"component-counts", TestComponentCounts},
        {"storage-shape", TestStorageShape},
        {"known-accumulation-across-field", TestKnownAccumulationAcrossField},
        {"weighted-accumulation", TestWeightedAccumulation},
        {"masked-points-excluded", TestMaskedPointsAreExcluded},
        {"constant-field-zero-product", TestConstantFieldZeroProduct},
        {"valid-fraction-range", TestValidFractionRange},
        {"valid-fraction-detects-never-valid-point", TestValidFractionDetectsNeverValidPoint},
        {"covariance-accumulation", TestCovarianceAccumulation},
        {"scalar-self-covariance-matches-second-moment", TestScalarSelfCovarianceMatchesSecondMoment},
        {"covariance-requires-field-membership", TestCovarianceRequiresFieldMembership},
        {"payload-enumeration", TestPayloadEnumeration},
        {"derived-quantities", TestDerivedQuantities},
        {"spatial-mean-excludes-unsampled-points", TestSpatialMeanExcludesUnsampledPoints},
        {"runloop-driver-applies-scheduled-states", TestRunloopDriverAppliesScheduledStates},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv statistics accumulator tests");
    if (ierr) return (int)ierr;
    ierr = PicurvRunTests("unit-statistics-accumulator", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) { PetscFinalize(); return (int)ierr; }
    ierr = PetscFinalize();
    return (int)ierr;
}
