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
    /* Ucat contributes six product components and P one. */
    PetscCall(PicurvAssertIntEqual(0, storage.m2_offset[0], "first field starts at offset zero"));
    PetscCall(PicurvAssertIntEqual(6, storage.m2_offset[1], "Ucat occupies six product components"));
    PetscCall(PicurvAssertIntEqual(7, storage.m2_offset[2], "P adds one further component"));

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
        PetscCall(ReadScalarAt(user, storage.m2[storage.m2_offset[0] + c], 2, 2, 2, &value));
        PetscCall(PetscSNPrintf(context, sizeof(context), "Ucat centered product %s", labels[c]));
        PetscCall(PicurvAssertRealNear(expected_product[c], value, 1.0e-11, context));
    }

    /* Scalar mean and second moment. */
    PetscCall(ReadScalarAt(user, storage.mean[1], 2, 2, 2, &value));
    PetscCall(PicurvAssertRealNear(3.0, value, 1.0e-12, "P mean"));
    PetscCall(ReadScalarAt(user, storage.m2[storage.m2_offset[1]], 2, 2, 2, &value));
    PetscCall(PicurvAssertRealNear(14.0, value, 1.0e-11, "P centered second moment"));

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
    PetscCall(ReadScalarAt(user, storage.m2[storage.m2_offset[1]], 2, 2, 2, &value));
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
        PetscCall(ReadScalarAt(user, storage.m2[storage.m2_offset[0] + c], 2, 2, 2, &value));
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
    PetscCall(ReadScalarAt(user, storage.m2[storage.m2_offset[1]], 2, 2, 2, &value));
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
        {"runloop-driver-applies-scheduled-states", TestRunloopDriverAppliesScheduledStates},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv statistics accumulator tests");
    if (ierr) return (int)ierr;
    ierr = PicurvRunTests("unit-statistics-accumulator", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) { PetscFinalize(); return (int)ierr; }
    ierr = PetscFinalize();
    return (int)ierr;
}
