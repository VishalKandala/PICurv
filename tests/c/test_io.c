/**
 * @file test_io.c
 * @brief C unit tests for I/O helpers, parsers, and startup-banner output.
 */

#include "test_support.h"

#include "checksum.h"
#include "io.h"
#include "statistics_accumulator.h"
#include "statistics_window.h"
#include "field_catalog.h"

#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
/**
 * @brief Tests cadence-based Eulerian output triggering.
 */

static PetscErrorCode TestShouldWriteDataOutput(void)
{
    SimCtx simCtx;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&simCtx, sizeof(simCtx)));
    simCtx.tiout = 5;

    PetscCall(PicurvAssertBool((PetscBool)!ShouldWriteDataOutput(NULL, 5), "NULL SimCtx should never request output"));
    PetscCall(PicurvAssertBool((PetscBool)!ShouldWriteDataOutput(&simCtx, 4), "non-cadence step should not trigger output"));
    PetscCall(PicurvAssertBool(ShouldWriteDataOutput(&simCtx, 10), "cadence-aligned step should trigger output"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests filesystem existence checks for files and directories.
 */

static PetscErrorCode TestVerifyPathExistence(void)
{
    char tmpdir[PETSC_MAX_PATH_LEN];
    char filepath[PETSC_MAX_PATH_LEN];
    FILE *file = NULL;
    PetscBool exists = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(filepath, sizeof(filepath), "%s/sample.txt", tmpdir));

    file = fopen(filepath, "w");
    PetscCheck(file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to create temp file '%s'.", filepath);
    fputs("picurv\n", file);
    fclose(file);

    PetscCall(VerifyPathExistence(tmpdir, PETSC_TRUE, PETSC_FALSE, "temp directory", &exists));
    PetscCall(PicurvAssertBool(exists, "VerifyPathExistence should find the temp directory"));

    PetscCall(VerifyPathExistence(filepath, PETSC_FALSE, PETSC_FALSE, "temp file", &exists));
    PetscCall(PicurvAssertBool(exists, "VerifyPathExistence should find the temp file"));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests writing and reloading core Eulerian field vectors.
 */

static PetscErrorCode TestWriteAndReadSimulationFields(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char path[PETSC_MAX_PATH_LEN];
    PetscBool exists = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    tmpdir[0] = '\0';
    if (simCtx->rank == 0) PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCallMPI(MPI_Bcast(tmpdir, sizeof(tmpdir), MPI_CHAR, 0, PETSC_COMM_WORLD));

    PetscCall(PetscStrncpy(simCtx->output_dir, tmpdir, sizeof(simCtx->output_dir)));
    PetscCall(PetscStrncpy(simCtx->restart_dir, tmpdir, sizeof(simCtx->restart_dir)));
    PetscCall(VecSet(user->P, 4.5));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(VecSet(user->Ucat, 2.0));
    PetscCall(VecSet(user->Ucont, 3.0));
    PetscCall(VecSet(user->Ucont_rm1, 1.25));
    simCtx->ti = 0.35;
    PetscCall(PicurvPopulateIdentityMetrics(user));

    PetscCall(WriteCheckpointBundle(simCtx, "test"));
    PetscCall(PetscSNPrintf(path, sizeof(path),
                            "%s/checkpoints/step_000000000001/checkpoint.meta", tmpdir));
    PetscCall(PetscTestFile(path, 'r', &exists));
    PetscCall(PicurvAssertBool(exists, "checkpoint coordinator should write checkpoint.meta"));
    PetscCall(PetscSNPrintf(path, sizeof(path),
                            "%s/checkpoints/step_000000000001/COMMITTED", tmpdir));
    PetscCall(PetscTestFile(path, 'r', &exists));
    PetscCall(PicurvAssertBool(exists, "checkpoint coordinator should write COMMITTED last"));
    PetscCall(PetscSNPrintf(path, sizeof(path),
                            "%s/checkpoints/step_000000000001/eulerian/block_0000/Ucat.dat", tmpdir));
    PetscCall(PetscTestFile(path, 'r', &exists));
    PetscCall(PicurvAssertBool(exists, "checkpoint should use catalogued canonical field names"));
    PetscCall(VecZeroEntries(user->P));
    PetscCall(VecZeroEntries(user->Ucat));
    PetscCall(VecZeroEntries(user->Ucont));
    PetscCall(VecZeroEntries(user->Ucont_rm1));
    simCtx->ti = 0.0;

    PetscCall(ReadSimulationFields(user, simCtx->step));
    PetscCall(PicurvAssertVecConstant(user->P, 4.5, 1.0e-12, "ReadSimulationFields should restore P"));
    PetscCall(PicurvAssertVecConstant(user->Ucat, 2.0, 1.0e-12, "ReadSimulationFields should restore Ucat"));
    PetscCall(PicurvAssertVecConstant(user->Ucont, 3.0, 1.0e-12, "ReadSimulationFields should restore Ucont"));
    PetscCall(PicurvAssertVecConstant(user->Ucont_rm1, 1.25, 1.0e-12,
                                      "ReadSimulationFields should restore BDF2 history"));
    PetscCall(PicurvAssertRealNear(0.35, simCtx->ti, 1.0e-12,
                                   "checkpoint physical time should be authoritative"));

    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    if (simCtx->rank == 0) PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief Builds the statistics fixture: one Ucat/P window with a second moment. */
static PicurvWindowDefinition StatisticsFixtureDefinition(void)
{
    PicurvWindowDefinition d;

    memset(&d, 0, sizeof(d));
    strncpy(d.name, "production", PICURV_WINDOW_NAME_LENGTH - 1);
    d.weighting = PICURV_WEIGHTING_SAMPLE;
    d.cadence_kind = PICURV_CADENCE_STEP;
    d.step_cadence = 1;
    d.field_count = 2;
    d.fields[0].field_id = FIELD_ID_UCAT; d.fields[0].want_second = PETSC_TRUE;
    d.fields[1].field_id = FIELD_ID_P;    d.fields[1].want_second = PETSC_TRUE;
    return d;
}

/** @brief Attaches one accumulating window to a fixture context and primes it with samples. */
static PetscErrorCode AttachStatisticsFixture(SimCtx *simCtx, UserCtx *user,
                                              PicurvWindow *window, PicurvWindowStorage *storage,
                                              const PicurvWindowDefinition *definition)
{
    PetscFunctionBeginUser;
    PetscCall(PicurvWindowInit(window, definition));
    PetscCall(PicurvWindowStorageCreate(user, definition, storage));
    simCtx->fieldStatisticsEnabled = PETSC_TRUE;
    simCtx->fieldStatisticsWindowCount = 1;
    simCtx->fieldStatisticsWindows = window;
    user->fieldStatisticsStorage = storage;
    PetscFunctionReturn(0);
}

/** @brief Detaches the statistics fixture without disturbing the shared teardown. */
static PetscErrorCode DetachStatisticsFixture(SimCtx *simCtx, UserCtx *user,
                                              PicurvWindowStorage *storage)
{
    PetscFunctionBeginUser;
    simCtx->fieldStatisticsEnabled = PETSC_FALSE;
    simCtx->fieldStatisticsWindowCount = 0;
    simCtx->fieldStatisticsWindows = NULL;
    user->fieldStatisticsStorage = NULL;
    PetscCall(PicurvWindowStorageDestroy(storage));
    PetscFunctionReturn(0);
}

/**
 * @brief Rewrites the bundle holding a window that closed well before the checkpoint time.
 *
 * Reproduces the shape a real run reaches when a bounded window finished and the
 * simulation kept going: the saved state is complete, and its end sits behind the
 * checkpoint by far more than one step.
 */
static PetscErrorCode WriteCompletedStatisticsWindow(SimCtx *simCtx, UserCtx *user,
                                                     PicurvWindow *window,
                                                     PicurvWindowStorage *storage,
                                                     const PicurvWindowDefinition *definition)
{
    PicurvWindowDefinition bounded = *definition;
    PetscBool accepted = PETSC_FALSE;
    PetscReal weight = 0.0;
    char step_directory[PETSC_MAX_PATH_LEN];

    PetscFunctionBeginUser;
    bounded.bounded = PETSC_TRUE;
    bounded.end_time = 1.0;
    PetscCall(PicurvWindowInit(window, &bounded));
    PetscCall(PicurvWindowOfferState(window, 0, 0.0, &accepted, &weight));
    PetscCall(PicurvWindowOfferState(window, 1, 1.0, &accepted, &weight));
    if (accepted) PetscCall(PicurvWindowAccumulate(user, &bounded, storage, weight));
    PetscCall(PicurvAssertBool((PetscBool)(window->state == PICURV_WINDOW_COMPLETE),
                               "the fixture window reaches its bounded end"));

    /* The run continued well past the window's end before this checkpoint. */
    simCtx->ti = 5.0;
    PetscCall(PetscSNPrintf(step_directory, sizeof(step_directory),
                            "%s/checkpoints/step_000000000001", simCtx->output_dir));
    if (simCtx->rank == 0) PetscCall(PetscRMTree(step_directory));
    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    PetscCall(WriteCheckpointBundle(simCtx, "test"));
    PetscFunctionReturn(0);
}

/**
 * @brief Verifies accumulated statistics survive a checkpoint round trip unchanged.
 *
 * Both halves of the state matter and are checked separately: the per-point
 * accumulator vectors, and the window's scalar bookkeeping. Restoring only the
 * vectors would silently resume with a broken schedule, and restoring only the
 * scalars would report a sample count the fields do not support.
 */
static PetscErrorCode TestCheckpointStatisticsRoundTrip(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindow window;
    PicurvWindowStorage storage;
    const PicurvWindowDefinition definition = StatisticsFixtureDefinition();
    char tmpdir[PETSC_MAX_PATH_LEN];
    char path[PETSC_MAX_PATH_LEN];
    PetscBool exists = PETSC_FALSE;
    PetscInt payload_count = 0;
    Vec *reference = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    tmpdir[0] = '\0';
    if (simCtx->rank == 0) PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCallMPI(MPI_Bcast(tmpdir, sizeof(tmpdir), MPI_CHAR, 0, PETSC_COMM_WORLD));
    PetscCall(PetscStrncpy(simCtx->output_dir, tmpdir, sizeof(simCtx->output_dir)));
    PetscCall(PetscStrncpy(simCtx->restart_dir, tmpdir, sizeof(simCtx->restart_dir)));
    PetscCall(PicurvPopulateIdentityMetrics(user));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(AttachStatisticsFixture(simCtx, user, &window, &storage, &definition));

    /* Three accepted states at unit weight: P takes 1, 2 and 6, so its per-point
     * mean is 3 and its centered second moment is 14. */
    {
        const PetscReal scalars[3] = {1.0, 2.0, 6.0};

        for (PetscInt s = 0; s < 3; ++s) {
            PetscBool accepted = PETSC_FALSE;
            PetscReal weight = 0.0;

            PetscCall(VecSet(user->P, scalars[s]));
            PetscCall(VecSet(user->Ucat, scalars[s]));
            PetscCall(PicurvWindowOfferState(&window, s, (PetscReal)s, &accepted, &weight));
            if (accepted) PetscCall(PicurvWindowAccumulate(user, &definition, &storage, weight));
        }
    }
    PetscCall(PicurvAssertIntEqual(2, window.sample_count, "the anchoring state is not itself a sample"));

    PetscCall(WriteCheckpointBundle(simCtx, "test"));
    PetscCall(PetscSNPrintf(path, sizeof(path),
                            "%s/checkpoints/step_000000000001/statistics/window_0000/block_0000/P_mean.dat",
                            tmpdir));
    PetscCall(PetscTestFile(path, 'r', &exists));
    PetscCall(PicurvAssertBool(exists, "statistics payloads use their enumerated names"));
    PetscCall(PetscSNPrintf(path, sizeof(path),
                            "%s/checkpoints/step_000000000001/statistics/window_0000/block_0000/Ucat_m2.dat",
                            tmpdir));
    PetscCall(PetscTestFile(path, 'r', &exists));
    PetscCall(PicurvAssertBool(exists, "the symmetric product is one payload, not six"));

    /* Keep a reference copy of every payload, then discard the in-memory state so a
     * restore that did nothing cannot pass. */
    PetscCall(PicurvWindowStoragePayloadCount(&storage, &payload_count));
    PetscCall(PetscCalloc1((size_t)payload_count, &reference));
    for (PetscInt index = 0; index < payload_count; ++index) {
        PicurvStatisticsPayload payload;

        PetscCall(PicurvWindowStoragePayload(user, &definition, &storage, index, &payload));
        PetscCall(VecDuplicate(payload.vec, &reference[index]));
        PetscCall(VecCopy(payload.vec, reference[index]));
        PetscCall(VecZeroEntries(payload.vec));
    }
    PetscCall(PicurvWindowInit(&window, &definition));
    PetscCall(PicurvAssertIntEqual(0, window.sample_count, "the fixture is reset before restoring"));

    simCtx->fieldStatisticsContinue = PETSC_TRUE;
    PetscCall(RestoreFieldStatisticsState(simCtx, simCtx->step));

    PetscCall(PicurvAssertIntEqual(2, window.sample_count, "the sample count is restored"));
    PetscCall(PicurvAssertRealNear(2.0, window.total_weight, 1.0e-12, "the total weight is restored"));
    PetscCall(PicurvAssertRealNear(2.0, window.represented_time, 1.0e-12,
                                   "the represented time is restored"));
    PetscCall(PicurvAssertRealNear(2.0, window.last_accepted_time, 1.0e-12,
                                   "the quadrature origin is restored"));
    PetscCall(PicurvAssertIntEqual(0, window.activation_step, "the schedule anchor is restored"));
    PetscCall(PicurvAssertIntEqual(2, window.last_event_step,
                                   "the duplicate-event guard survives the restart"));
    PetscCall(PicurvAssertIntEqual(1, window.restart_count, "the restart lineage advances on resume"));
    PetscCall(PicurvAssertBool((PetscBool)(window.state == PICURV_WINDOW_ACTIVE),
                               "an open window resumes active"));

    /* Every payload must come back bit for bit. The binary format stores IEEE
     * doubles and the natural ordering is decomposition independent, so anything
     * short of an exact match is a defect rather than rounding. */
    for (PetscInt index = 0; index < payload_count; ++index) {
        PicurvStatisticsPayload payload;
        PetscReal difference = 0.0;
        char context[192];

        PetscCall(PicurvWindowStoragePayload(user, &definition, &storage, index, &payload));
        PetscCall(VecAXPY(reference[index], -1.0, payload.vec));
        PetscCall(VecNorm(reference[index], NORM_INFINITY, &difference));
        PetscCall(PetscSNPrintf(context, sizeof(context),
                                "payload '%s' is restored bit for bit", payload.name));
        PetscCall(PicurvAssertBool((PetscBool)(difference == 0.0), context));
        PetscCall(VecDestroy(&reference[index]));
    }
    PetscCall(PetscFree(reference));

    /* Spot-check one interior point against the analytically known value, so the
     * comparison above cannot pass by restoring two identically wrong vectors. */
    {
        PetscReal ***mean = NULL;

        PetscCall(DMDAVecGetArrayRead(user->da, storage.mean[1], &mean));
        if (user->info.xs <= 2 && 2 < user->info.xs + user->info.xm &&
            user->info.ys <= 2 && 2 < user->info.ys + user->info.ym &&
            user->info.zs <= 2 && 2 < user->info.zs + user->info.zm) {
            PetscCall(PicurvAssertRealNear(4.0, mean[2][2][2], 1.0e-12,
                                           "the restored scalar mean holds the accumulated value"));
        }
        PetscCall(DMDAVecRestoreArrayRead(user->da, storage.mean[1], &mean));
    }

    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    if (simCtx->rank == 0) PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    PetscCall(DetachStatisticsFixture(simCtx, user, &storage));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Verifies continuation refuses every state that would corrupt an average.
 *
 * A silent reset is the failure mode this guards: resuming from zero, or merging
 * samples taken under a different definition, both produce a window whose reported
 * sample count no longer describes the numbers it carries.
 */
static PetscErrorCode TestCheckpointStatisticsContinuationGuards(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindow window;
    PicurvWindowStorage storage;
    const PicurvWindowDefinition definition = StatisticsFixtureDefinition();
    PicurvWindowDefinition changed;
    char tmpdir[PETSC_MAX_PATH_LEN];
    PetscErrorCode bad = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    tmpdir[0] = '\0';
    if (simCtx->rank == 0) PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCallMPI(MPI_Bcast(tmpdir, sizeof(tmpdir), MPI_CHAR, 0, PETSC_COMM_WORLD));
    PetscCall(PetscStrncpy(simCtx->output_dir, tmpdir, sizeof(simCtx->output_dir)));
    PetscCall(PetscStrncpy(simCtx->restart_dir, tmpdir, sizeof(simCtx->restart_dir)));
    PetscCall(PicurvPopulateIdentityMetrics(user));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(AttachStatisticsFixture(simCtx, user, &window, &storage, &definition));
    {
        PetscBool accepted = PETSC_FALSE;
        PetscReal weight = 0.0;

        PetscCall(PicurvWindowOfferState(&window, 0, 0.0, &accepted, &weight));
        PetscCall(PicurvWindowOfferState(&window, 1, 1.0, &accepted, &weight));
        if (accepted) PetscCall(PicurvWindowAccumulate(user, &definition, &storage, weight));
    }
    PetscCall(WriteCheckpointBundle(simCtx, "test"));

    /* A hashed property changed: the saved samples describe a different average. */
    changed = definition;
    changed.weighting = PICURV_WEIGHTING_PHYSICAL_TIME;
    PetscCall(PicurvWindowInit(&window, &changed));
    simCtx->fieldStatisticsContinue = PETSC_TRUE;
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    bad = RestoreFieldStatisticsState(simCtx, simCtx->step);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_INCOMP, bad,
                                   "a changed definition is refused rather than merged"));

    /* A window renamed under the same index is a different window. */
    changed = definition;
    strncpy(changed.name, "other", PICURV_WINDOW_NAME_LENGTH - 1);
    PetscCall(PicurvWindowInit(&window, &changed));
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    bad = RestoreFieldStatisticsState(simCtx, simCtx->step);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_INCOMP, bad, "a renamed window is refused"));

    /* Shortening a window would discard represented time it already claims. */
    changed = definition;
    changed.bounded = PETSC_TRUE;
    changed.end_time = 0.5;
    PetscCall(PicurvWindowInit(&window, &changed));
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    bad = RestoreFieldStatisticsState(simCtx, simCtx->step);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_INCOMP, bad, "a shortened window is refused"));

    /* Extending it forward is permitted, because end_time is outside the hash. The
     * saved window is still open here, so there is no former end to leave a gap
     * after and the extension is accepted. */
    changed = definition;
    changed.bounded = PETSC_TRUE;
    changed.end_time = 40.0;
    PetscCall(PicurvWindowInit(&window, &changed));
    PetscCall(RestoreFieldStatisticsState(simCtx, simCtx->step));
    PetscCall(PicurvAssertIntEqual(1, window.sample_count, "an extended window keeps its samples"));
    PetscCall(PicurvAssertBool((PetscBool)(window.state == PICURV_WINDOW_ACTIVE),
                               "an extended window resumes active"));

    /* Reopening a window that already closed, across time it never sampled, would
     * weight that gap into the first new interval. */
    simCtx->dt = 0.1;
    PetscCall(WriteCompletedStatisticsWindow(simCtx, user, &window, &storage, &definition));
    changed = definition;
    changed.bounded = PETSC_TRUE;
    changed.end_time = 40.0;
    PetscCall(PicurvWindowInit(&window, &changed));
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    bad = RestoreFieldStatisticsState(simCtx, simCtx->step);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_INCOMP, bad,
                                   "reopening a closed window across an unsampled gap is refused"));

    /* Without an explicit request, a restart starts from zero and says so. */
    PetscCall(PicurvWindowInit(&window, &definition));
    simCtx->fieldStatisticsContinue = PETSC_FALSE;
    PetscCall(RestoreFieldStatisticsState(simCtx, simCtx->step));
    PetscCall(PicurvAssertIntEqual(0, window.sample_count,
                                   "statistics start fresh when continuation is not requested"));

    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    if (simCtx->rank == 0) PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    PetscCall(DetachStatisticsFixture(simCtx, user, &storage));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief A damaged statistics payload must fail bundle validation, not load silently.
 *
 * Statistics payloads are validated because they enter the manifest inventory the
 * existing validator already walks. That is an easy property to believe and an easy
 * one to lose, so it is checked by damaging a payload rather than by inspection.
 */
static PetscErrorCode TestCheckpointStatisticsPayloadIsValidated(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindow window;
    PicurvWindowStorage storage;
    const PicurvWindowDefinition definition = StatisticsFixtureDefinition();
    char tmpdir[PETSC_MAX_PATH_LEN];
    char payload_path[PETSC_MAX_PATH_LEN];
    PetscErrorCode bad = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    tmpdir[0] = '\0';
    if (simCtx->rank == 0) PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCallMPI(MPI_Bcast(tmpdir, sizeof(tmpdir), MPI_CHAR, 0, PETSC_COMM_WORLD));
    PetscCall(PetscStrncpy(simCtx->output_dir, tmpdir, sizeof(simCtx->output_dir)));
    PetscCall(PetscStrncpy(simCtx->restart_dir, tmpdir, sizeof(simCtx->restart_dir)));
    PetscCall(PicurvPopulateIdentityMetrics(user));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(AttachStatisticsFixture(simCtx, user, &window, &storage, &definition));
    {
        PetscBool accepted = PETSC_FALSE;
        PetscReal weight = 0.0;

        PetscCall(PicurvWindowOfferState(&window, 0, 0.0, &accepted, &weight));
        PetscCall(PicurvWindowOfferState(&window, 1, 1.0, &accepted, &weight));
        if (accepted) PetscCall(PicurvWindowAccumulate(user, &definition, &storage, weight));
    }
    PetscCall(WriteCheckpointBundle(simCtx, "test"));

    /* Truncating one payload leaves the manifest and its commit marker intact, so
     * only the inventory's recorded byte size can catch it. */
    PetscCall(PetscSNPrintf(payload_path, sizeof(payload_path),
                            "%s/checkpoints/step_000000000001/statistics/window_0000/block_0000/P_mean.dat",
                            tmpdir));
    if (simCtx->rank == 0) {
        FILE *damaged = fopen(payload_path, "w");

        PetscCheck(damaged != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
                   "Unable to truncate '%s'.", payload_path);
        PetscCheck(fclose(damaged) == 0, PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
                   "Unable to close '%s'.", payload_path);
    }
    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));

    simCtx->fieldStatisticsContinue = PETSC_TRUE;
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    bad = RestoreFieldStatisticsState(simCtx, simCtx->step);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_FILE_UNEXPECTED, bad,
                                   "a truncated statistics payload fails bundle validation"));

    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    if (simCtx->rank == 0) PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    PetscCall(DetachStatisticsFixture(simCtx, user, &storage));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief A run without statistics writes no statistics subtree and refuses to fake one. */
static PetscErrorCode TestCheckpointStatisticsAbsentWhenDisabled(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindow window;
    PicurvWindowStorage storage;
    const PicurvWindowDefinition definition = StatisticsFixtureDefinition();
    char tmpdir[PETSC_MAX_PATH_LEN];
    char path[PETSC_MAX_PATH_LEN];
    PetscBool exists = PETSC_TRUE;
    PetscErrorCode bad = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    tmpdir[0] = '\0';
    if (simCtx->rank == 0) PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCallMPI(MPI_Bcast(tmpdir, sizeof(tmpdir), MPI_CHAR, 0, PETSC_COMM_WORLD));
    PetscCall(PetscStrncpy(simCtx->output_dir, tmpdir, sizeof(simCtx->output_dir)));
    PetscCall(PetscStrncpy(simCtx->restart_dir, tmpdir, sizeof(simCtx->restart_dir)));
    PetscCall(PicurvPopulateIdentityMetrics(user));

    PetscCall(WriteCheckpointBundle(simCtx, "test"));
    PetscCall(PetscSNPrintf(path, sizeof(path), "%s/checkpoints/step_000000000001/statistics", tmpdir));
    PetscCall(PetscTestDirectory(path, 'r', &exists));
    PetscCall(PicurvAssertBool((PetscBool)!exists,
                               "a run without statistics writes no statistics subtree"));

    /* Continuing from that bundle must fail loudly rather than resume from zero. */
    PetscCall(AttachStatisticsFixture(simCtx, user, &window, &storage, &definition));
    simCtx->fieldStatisticsContinue = PETSC_TRUE;
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    bad = RestoreFieldStatisticsState(simCtx, simCtx->step);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_INCOMP, bad,
                                   "missing statistics state is fatal, never silently zeroed"));

    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    if (simCtx->rank == 0) PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    PetscCall(DetachStatisticsFixture(simCtx, user, &storage));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Verifies a committed checkpoint step is rewritten neither silently nor inconsistently.
 *
 * `WriteCheckpointBundle` revalidates an already-committed step instead of rewriting it,
 * so a repeated call at the same step must leave the payloads byte-identical even when the
 * in-memory state has since diverged. The same guard must reject a repeat call whose
 * physical time disagrees with the committed bundle, because that means the step number no
 * longer identifies the same solver state.
 */
static PetscErrorCode TestCheckpointSameStepRewriteIsRejected(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char payload_path[PETSC_MAX_PATH_LEN];
    struct stat first_stat;
    struct stat second_stat;
    PetscErrorCode mismatch_ierr = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    tmpdir[0] = '\0';
    if (simCtx->rank == 0) PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCallMPI(MPI_Bcast(tmpdir, sizeof(tmpdir), MPI_CHAR, 0, PETSC_COMM_WORLD));

    PetscCall(PetscStrncpy(simCtx->output_dir, tmpdir, sizeof(simCtx->output_dir)));
    PetscCall(PetscStrncpy(simCtx->restart_dir, tmpdir, sizeof(simCtx->restart_dir)));
    PetscCall(VecSet(user->P, 4.5));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(VecSet(user->Ucat, 2.0));
    PetscCall(VecSet(user->Ucont, 3.0));
    PetscCall(VecSet(user->Ucont_rm1, 1.25));
    simCtx->ti = 0.35;
    PetscCall(PicurvPopulateIdentityMetrics(user));

    PetscCall(WriteCheckpointBundle(simCtx, "cadence"));
    PetscCall(PetscSNPrintf(payload_path, sizeof(payload_path),
                            "%s/checkpoints/step_000000000001/eulerian/block_0000/Ucat.dat", tmpdir));
    PetscCall(PicurvAssertBool((PetscBool)(stat(payload_path, &first_stat) == 0),
                               "first checkpoint write should produce a Ucat payload"));

    /* Diverge the in-memory state so a rewrite would be observable in the payload bytes. */
    PetscCall(VecSet(user->Ucat, 99.0));
    PetscCall(WriteCheckpointBundle(simCtx, "cadence"));

    PetscCall(PicurvAssertBool((PetscBool)(stat(payload_path, &second_stat) == 0),
                               "repeated checkpoint write should leave the Ucat payload in place"));
    /* Nanosecond resolution: a same-second rewrite would still move this. */
    PetscCall(PicurvAssertBool((PetscBool)(first_stat.st_mtim.tv_sec == second_stat.st_mtim.tv_sec &&
                                           first_stat.st_mtim.tv_nsec == second_stat.st_mtim.tv_nsec),
                               "repeated checkpoint write at a committed step should not touch payload mtime"));
    PetscCall(PicurvAssertBool((PetscBool)(first_stat.st_size == second_stat.st_size),
                               "repeated checkpoint write at a committed step should not resize payloads"));

    PetscCall(VecZeroEntries(user->Ucat));
    PetscCall(ReadSimulationFields(user, simCtx->step));
    PetscCall(PicurvAssertVecConstant(user->Ucat, 2.0, 1.0e-12,
                                      "repeated checkpoint write must not overwrite committed payload contents"));

    /* A same-step call whose physical time disagrees is a state-identity error, not a no-op. */
    simCtx->ti = 0.75;
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    mismatch_ierr = WriteCheckpointBundle(simCtx, "cadence");
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_FILE_UNEXPECTED, mismatch_ierr,
                                   "same-step checkpoint write with a different physical time should fail"));

    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    if (simCtx->rank == 0) PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief Verifies the dependency-free SHA-256 implementation against a standard vector. */
static PetscErrorCode TestCheckpointSHA256KnownVector(void)
{
    PicurvSHA256Context context;
    char digest[65];

    PetscFunctionBeginUser;
    PicurvSHA256Init(&context);
    PicurvSHA256Update(&context, "abc", 3);
    PicurvSHA256FinalHex(&context, digest);
    PetscCall(PicurvAssertBool(
        (PetscBool)!strcmp(digest, "ba7816bf8f01cfea414140de5dae2223"
                                   "b00361a396177a9cb410ff61f20015ad"),
        "SHA-256 implementation should match the standard abc test vector"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests parsing of post-processing control settings from a file.
 */

static PetscErrorCode TestParsePostProcessingSettings(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char cfg_path[PETSC_MAX_PATH_LEN];
    FILE *file = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PetscCalloc1(1, &simCtx->pps));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(cfg_path, sizeof(cfg_path), "%s/post.run", tmpdir));
    PetscCall(PetscStrncpy(simCtx->PostprocessingControlFile, cfg_path, sizeof(simCtx->PostprocessingControlFile)));

    file = fopen(cfg_path, "w");
    PetscCheck(file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to create temp config file '%s'.", cfg_path);
    fputs("startTime = 2\n", file);
    fputs("endTime = 6\n", file);
    fputs("timeStep = 2\n", file);
    fputs("output_particles = true\n", file);
    fputs("output_prefix = SmokeField\n", file);
    fclose(file);

    PetscCall(ParsePostProcessingSettings(simCtx));
    PetscCall(PicurvAssertIntEqual(2, simCtx->pps->startTime, "ParsePostProcessingSettings should parse startTime"));
    PetscCall(PicurvAssertIntEqual(6, simCtx->pps->endTime, "ParsePostProcessingSettings should parse endTime"));
    PetscCall(PicurvAssertIntEqual(2, simCtx->pps->timeStep, "ParsePostProcessingSettings should parse timeStep"));
    PetscCall(PicurvAssertBool(simCtx->pps->outputParticles, "ParsePostProcessingSettings should parse output_particles"));
    PetscCall(PicurvAssertBool((PetscBool)(strcmp(simCtx->pps->output_prefix, "SmokeField") == 0),
                               "ParsePostProcessingSettings should parse output_prefix"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests trimming of leading and trailing whitespace.
 */

static PetscErrorCode TestTrimWhitespace(void)
{
    char value_a[] = "   inlet_value   ";
    char value_b[] = "   ";

    PetscFunctionBeginUser;
    TrimWhitespace(value_a);
    PetscCall(PicurvAssertBool((PetscBool)(strcmp(value_a, "inlet_value") == 0),
                               "TrimWhitespace should remove leading and trailing whitespace"));

    TrimWhitespace(value_b);
    PetscCall(PicurvAssertBool((PetscBool)(strcmp(value_b, "") == 0),
                               "TrimWhitespace should reduce all-whitespace strings to empty"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests boundary-condition string parsers for face, type, and handler names.
 */

static PetscErrorCode TestBoundaryConditionStringParsers(void)
{
    BCFace face = BC_FACE_NEG_X;
    BCType type = WALL;
    BCHandlerType handler = BC_HANDLER_WALL_NOSLIP;

    PetscFunctionBeginUser;
    PetscCall(StringToBCFace("+Zeta", &face));
    PetscCall(PicurvAssertIntEqual(BC_FACE_POS_Z, face, "StringToBCFace should parse +Zeta"));

    PetscCall(StringToBCType("periodic", &type));
    PetscCall(PicurvAssertIntEqual(PERIODIC, type, "StringToBCType should parse PERIODIC case-insensitively"));

    PetscCall(StringToBCHandlerType("constant_flux", &handler));
    PetscCall(PicurvAssertIntEqual(BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX, handler,
                                   "StringToBCHandlerType should parse constant_flux"));
    PetscCall(StringToBCHandlerType("prescribed_flow", &handler));
    PetscCall(PicurvAssertIntEqual(BC_HANDLER_INLET_PROFILE_FROM_FILE, handler,
                                   "StringToBCHandlerType should parse prescribed_flow"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests validation of boundary-type and handler compatibility.
 */

static PetscErrorCode TestValidateBCHandlerForBCType(void)
{
    PetscFunctionBeginUser;
    PetscCall(PicurvAssertBool((PetscBool)(ValidateBCHandlerForBCType(WALL, BC_HANDLER_WALL_NOSLIP) == 0),
                               "WALL + noslip should be a valid combination"));
    PetscCall(PicurvAssertBool((PetscBool)(ValidateBCHandlerForBCType(PERIODIC, BC_HANDLER_PERIODIC_GEOMETRIC) == 0),
                               "PERIODIC + geometric should be a valid combination"));
    PetscCall(PicurvAssertBool((PetscBool)(ValidateBCHandlerForBCType(INLET, BC_HANDLER_INLET_PROFILE_FROM_FILE) == 0),
                               "INLET + prescribed_flow should be a valid combination"));
    PetscCall(PicurvAssertBool((PetscBool)(ValidateBCHandlerForBCType(INLET, BC_HANDLER_WALL_NOSLIP) != 0),
                               "INLET + noslip should be rejected"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests scaling-reference parsing and derived pressure scaling.
 */

static PetscErrorCode TestParseScalingInformation(void)
{
    SimCtx simCtx;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&simCtx, sizeof(simCtx)));

    PetscCall(PetscOptionsClearValue(NULL, "-scaling_L_ref"));
    PetscCall(PetscOptionsClearValue(NULL, "-scaling_U_ref"));
    PetscCall(PetscOptionsClearValue(NULL, "-scaling_rho_ref"));

    PetscCall(ParseScalingInformation(&simCtx));
    PetscCall(PicurvAssertRealNear(1.0, simCtx.scaling.L_ref, 1.0e-12, "Default scaling_L_ref should be 1.0"));
    PetscCall(PicurvAssertRealNear(1.0, simCtx.scaling.U_ref, 1.0e-12, "Default scaling_U_ref should be 1.0"));
    PetscCall(PicurvAssertRealNear(1.0, simCtx.scaling.rho_ref, 1.0e-12, "Default scaling_rho_ref should be 1.0"));
    PetscCall(PicurvAssertRealNear(1.0, simCtx.scaling.P_ref, 1.0e-12, "Default scaling_P_ref should be 1.0"));

    PetscCall(PetscOptionsSetValue(NULL, "-scaling_L_ref", "2.5"));
    PetscCall(PetscOptionsSetValue(NULL, "-scaling_U_ref", "4.0"));
    PetscCall(PetscOptionsSetValue(NULL, "-scaling_rho_ref", "1.2"));

    PetscCall(ParseScalingInformation(&simCtx));
    PetscCall(PicurvAssertRealNear(2.5, simCtx.scaling.L_ref, 1.0e-12, "scaling_L_ref should honor options"));
    PetscCall(PicurvAssertRealNear(4.0, simCtx.scaling.U_ref, 1.0e-12, "scaling_U_ref should honor options"));
    PetscCall(PicurvAssertRealNear(1.2, simCtx.scaling.rho_ref, 1.0e-12, "scaling_rho_ref should honor options"));
    PetscCall(PicurvAssertRealNear(19.2, simCtx.scaling.P_ref, 1.0e-12, "scaling_P_ref should be rho_ref*U_ref^2"));

    PetscCall(PetscOptionsClearValue(NULL, "-scaling_L_ref"));
    PetscCall(PetscOptionsClearValue(NULL, "-scaling_U_ref"));
    PetscCall(PetscOptionsClearValue(NULL, "-scaling_rho_ref"));
    PetscFunctionReturn(0);
}
/**
 * @brief Captures the startup banner into a temporary file-backed buffer.
 */
static PetscErrorCode CaptureBannerOutput(SimCtx *simCtx, char *captured, size_t captured_len)
{
    char tmpdir[PETSC_MAX_PATH_LEN];
    char capture_path[PETSC_MAX_PATH_LEN];
    FILE *capture_file = NULL;
    int saved_stdout = -1;
    int capture_fd = -1;
    size_t bytes_read = 0;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    PetscCheck(simCtx != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "SimCtx cannot be NULL.");
    PetscCheck(captured != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Capture buffer cannot be NULL.");
    PetscCheck(captured_len > 0, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "Capture buffer must be non-empty.");

    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(capture_path, sizeof(capture_path), "%s/banner.log", tmpdir));

    fflush(stdout);
    saved_stdout = dup(STDOUT_FILENO);
    PetscCheck(saved_stdout >= 0, PETSC_COMM_SELF, PETSC_ERR_SYS, "dup(STDOUT_FILENO) failed.");
    capture_fd = open(capture_path, O_CREAT | O_TRUNC | O_WRONLY, 0600);
    PetscCheck(capture_fd >= 0, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to open capture file '%s'.", capture_path);
    PetscCheck(dup2(capture_fd, STDOUT_FILENO) >= 0, PETSC_COMM_SELF, PETSC_ERR_SYS, "dup2() failed while redirecting stdout.");
    close(capture_fd);
    capture_fd = -1;

    ierr = DisplayBanner(simCtx);
    fflush(stdout);
    PetscCheck(dup2(saved_stdout, STDOUT_FILENO) >= 0, PETSC_COMM_SELF, PETSC_ERR_SYS, "dup2() failed while restoring stdout.");
    close(saved_stdout);
    saved_stdout = -1;
    PetscCall(ierr);

    capture_file = fopen(capture_path, "r");
    PetscCheck(capture_file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to read capture file '%s'.", capture_path);
    bytes_read = fread(captured, 1, captured_len - 1, capture_file);
    captured[bytes_read] = '\0';
    fclose(capture_file);
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCallMPI(MPI_Bcast(captured, (PetscMPIInt)captured_len, MPI_CHAR, 0, PETSC_COMM_WORLD));
    PetscFunctionReturn(0);
}
/**
 * @brief Asserts that captured banner output contains one expected substring.
 */
static PetscErrorCode AssertCapturedContains(const char *captured, const char *needle, const char *message)
{
    PetscFunctionBeginUser;
    PetscCall(PicurvAssertBool((PetscBool)(strstr(captured, needle) != NULL), message));
    PetscFunctionReturn(0);
}
/**
 * @brief Asserts that captured banner output omits one forbidden substring.
 */
static PetscErrorCode AssertCapturedOmits(const char *captured, const char *needle, const char *message)
{
    PetscFunctionBeginUser;
    PetscCall(PicurvAssertBool((PetscBool)(strstr(captured, needle) == NULL), message));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that the startup banner reports statistics monitoring in every state.
 *
 * The banner is the one place a log records whether monitoring was active, so all three
 * states must be distinguishable after the fact: a live cadence, a subsystem that is
 * accumulating with the console silenced, and a run that configured no window at all.
 * The banner reads only the window array and the console cadence, so this fixture sets
 * those directly rather than allocating accumulator storage it would never touch.
 */

static PetscErrorCode TestDisplayBannerReportsStatisticsCadence(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindow window;
    PicurvWindowDefinition definition = StatisticsFixtureDefinition();
    char captured[16384];

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvWindowInit(&window, &definition));
    simCtx->OnlySetup = PETSC_FALSE;
    simCtx->StepsToRun = 5;
    PetscCall(PetscStrncpy(simCtx->eulerianSource, "solve", sizeof(simCtx->eulerianSource)));

    /* A run with statistics accumulating and the console cadence live. */
    simCtx->fieldStatisticsEnabled = PETSC_TRUE;
    simCtx->fieldStatisticsWindowCount = 1;
    simCtx->fieldStatisticsWindows = &window;
    simCtx->statisticsConsoleOutputFreq = 5;
    PetscCall(CaptureBannerOutput(simCtx, captured, sizeof(captured)));
    PetscCall(AssertCapturedContains(captured, "Statistics Console Cadence : every 5 step(s), 1 window(s)",
                                     "DisplayBanner should report the live statistics console cadence"));

    /* Accumulating, but with console reporting switched off. */
    simCtx->statisticsConsoleOutputFreq = 0;
    PetscCall(CaptureBannerOutput(simCtx, captured, sizeof(captured)));
    PetscCall(AssertCapturedContains(captured, "Statistics Console Cadence : DISABLED (1 window(s) accumulating)",
                                     "DisplayBanner should distinguish a silenced console from an inactive subsystem"));

    /* No window configured: the row still appears, so its absence is never ambiguous. */
    simCtx->fieldStatisticsEnabled = PETSC_FALSE;
    simCtx->fieldStatisticsWindowCount = 0;
    simCtx->fieldStatisticsWindows = NULL;
    PetscCall(CaptureBannerOutput(simCtx, captured, sizeof(captured)));
    PetscCall(AssertCapturedContains(captured, "Statistics Console Cadence : DISABLED (no window configured)",
                                     "DisplayBanner should record that no statistics window was configured"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests conditional startup-banner fields across particle and analytical cases.
 */

static PetscErrorCode TestDisplayBannerTracksConditionalStartupFields(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char captured[16384];

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    simCtx->OnlySetup = PETSC_FALSE;
    simCtx->StepsToRun = 5;
    simCtx->immersed = PETSC_FALSE;
    PetscCall(PetscStrncpy(simCtx->eulerianSource, "solve", sizeof(simCtx->eulerianSource)));
    simCtx->particleConsoleOutputFreq = 7;
    simCtx->LoggingFrequency = 4;
    simCtx->ParticleInitialization = PARTICLE_INIT_POINT_SOURCE;
    simCtx->mom_solver_type = MOMENTUM_SOLVER_DUALTIME_PICARD_JAMESON_RK;

    PetscCall(CaptureBannerOutput(simCtx, captured, sizeof(captured)));
    PetscCall(AssertCapturedContains(captured, "Run Mode                   : Full Simulation",
                                     "DisplayBanner should include the run mode"));
    PetscCall(AssertCapturedContains(captured, "Field/Restart Cadence      : every 2 step(s)",
                                     "DisplayBanner should include field/restart cadence"));
    PetscCall(AssertCapturedContains(captured, "Immersed Boundary          : DISABLED",
                                     "DisplayBanner should include immersed-boundary state"));
    PetscCall(AssertCapturedContains(captured, "Periodic Axes (BC-derived)  : I=NO, J=NO, K=NO",
                                     "DisplayBanner should include BC-derived periodic axes"));
    PetscCall(AssertCapturedContains(captured, "Number of Particles         : 0",
                                     "DisplayBanner should include the total particle count"));
    PetscCall(AssertCapturedContains(captured, "Eulerian State Source       : initial condition (Zero)",
                                     "DisplayBanner should identify a fresh solve IC as the Eulerian source"));
    PetscCall(AssertCapturedOmits(captured, "Particle Console Cadence",
                                  "DisplayBanner should omit particle console cadence when no particles are configured"));
    PetscCall(AssertCapturedOmits(captured, "Particle Log Row Sampling",
                                  "DisplayBanner should omit particle row sampling when no particles are configured"));
    PetscCall(AssertCapturedOmits(captured, "Particle Restart Mode",
                                  "DisplayBanner should omit particle restart mode when no particles are configured"));
    PetscCall(AssertCapturedOmits(captured, "Particle Initialization Mode",
                                  "DisplayBanner should omit particle initialization mode when no particles are configured"));
    PetscCall(AssertCapturedOmits(captured, "Interpolation Method",
                                  "DisplayBanner should omit interpolation method when no particles are configured"));
    PetscCall(AssertCapturedContains(captured, "Initial Pseudo-CFL (Courant)",
                                     "DisplayBanner should report pseudo-CFL for the dual-time momentum solver"));
    PetscCall(AssertCapturedContains(captured, "Pseudo-CFL Adaptation",
                                     "DisplayBanner should report the active dual-time controller"));
    PetscCall(AssertCapturedContains(captured, "Console Log Level",
                                     "DisplayBanner should report the effective console log level"));
    PetscCall(AssertCapturedContains(captured, "Profiling Timestep Output",
                                     "DisplayBanner should report profiling mode"));
    PetscCall(AssertCapturedContains(captured, "Runtime Memory Log",
                                     "DisplayBanner should report runtime-memory logging state"));

    simCtx->mom_solver_type = MOMENTUM_SOLVER_NEWTON_KRYLOV;
    PetscCall(CaptureBannerOutput(simCtx, captured, sizeof(captured)));
    PetscCall(AssertCapturedContains(captured, "Momentum Equation Solver    : Newton Krylov",
                                     "DisplayBanner should identify the selected Newton-Krylov solver"));
    PetscCall(AssertCapturedOmits(captured, "Initial Pseudo-CFL (Courant)",
                                  "DisplayBanner must not report pseudo-CFL for Newton-Krylov"));
    PetscCall(AssertCapturedContains(captured, "Newton-Krylov PETSc Controls",
                                     "DisplayBanner should report the active Newton-Krylov control family"));
    simCtx->mom_solver_type = MOMENTUM_SOLVER_DUALTIME_PICARD_JAMESON_RK;

    simCtx->solutionConvergenceMode = SOLUTION_CONVERGENCE_PERIODIC_DETERMINISTIC;
    simCtx->solutionConvergencePeriodSteps = 11;
    PetscCall(CaptureBannerOutput(simCtx, captured, sizeof(captured)));
    PetscCall(AssertCapturedContains(captured, "Solution Convergence Mode   : PERIODIC_DETERMINISTIC",
                                     "DisplayBanner should identify the active periodic convergence mode"));
    PetscCall(AssertCapturedContains(captured, "Convergence Period          : 11 step(s)",
                                     "DisplayBanner should report the active periodic convergence period"));
    simCtx->solutionConvergenceMode = SOLUTION_CONVERGENCE_STEADY_DETERMINISTIC;
    PetscCall(CaptureBannerOutput(simCtx, captured, sizeof(captured)));
    PetscCall(AssertCapturedOmits(captured, "Convergence Period",
                                  "DisplayBanner should omit periodic-only convergence fields when inactive"));

    simCtx->StartStep = 3;
    simCtx->np = 8;
    simCtx->i_periodic = 1;
    user->periodic_translation_valid[0] = PETSC_TRUE;
    user->periodic_translation[0] = (Cmpnts){1.0, 0.0, 0.0};
    simCtx->particleConsoleOutputFreq = 0;
    simCtx->LoggingFrequency = 4;
    PetscCall(PetscStrncpy(simCtx->particleRestartMode, "load", sizeof(simCtx->particleRestartMode)));
    PetscCall(CaptureBannerOutput(simCtx, captured, sizeof(captured)));
    PetscCall(AssertCapturedContains(captured, "Eulerian State Source       : restart step 3",
                                     "DisplayBanner should identify restart authority"));
    PetscCall(AssertCapturedOmits(captured, "initial condition (",
                                  "DisplayBanner should not present an IC as active during restart"));
    PetscCall(AssertCapturedContains(captured, "Number of Particles         : 8",
                                     "DisplayBanner should include the active particle count"));
    PetscCall(AssertCapturedContains(captured, "Periodic I Translation",
                                     "DisplayBanner should include validated periodic translation"));
    PetscCall(AssertCapturedContains(captured, "Particle Periodic Wrapping  : UNSUPPORTED",
                                     "DisplayBanner should distinguish Eulerian periodicity from particle wrapping"));
    PetscCall(AssertCapturedContains(captured, "Particle Console Cadence   : DISABLED",
                                     "DisplayBanner should show disabled particle console cadence when particles are configured"));
    PetscCall(AssertCapturedContains(captured, "Particle Log Row Sampling  : every 4 particle(s)",
                                     "DisplayBanner should include particle row sampling when particles are configured"));
    PetscCall(AssertCapturedContains(captured, "Particle Restart Mode      : load",
                                     "DisplayBanner should include particle restart mode for restarted particle runs"));
    PetscCall(AssertCapturedContains(captured, "Particle Initialization Mode: Point Source",
                                     "DisplayBanner should include particle initialization mode when particles are configured"));
    PetscCall(AssertCapturedContains(captured, "Interpolation Method       : Trilinear (direct cell-center)",
                                     "DisplayBanner should include default interpolation method when particles are configured"));
    PetscCall(AssertCapturedOmits(captured, "Particles Initialized At",
                                  "DisplayBanner should omit inlet-face placement details for point-source particle initialization"));

    simCtx->StartStep = 0;
    PetscCall(PetscStrncpy(simCtx->eulerianSource, "load", sizeof(simCtx->eulerianSource)));
    PetscCall(CaptureBannerOutput(simCtx, captured, sizeof(captured)));
    PetscCall(AssertCapturedContains(captured, "Eulerian State Source       : load",
                                     "DisplayBanner should identify load authority"));
    PetscCall(AssertCapturedOmits(captured, "initial condition (",
                                  "DisplayBanner should not present an IC as active in load mode"));

    simCtx->StartStep = 0;
    simCtx->particleConsoleOutputFreq = 6;
    simCtx->ParticleInitialization = PARTICLE_INIT_SURFACE_RANDOM;
    user->inletFaceDefined = PETSC_TRUE;
    user->identifiedInletBCFace = (BCFace)0;
    PetscCall(PetscStrncpy(simCtx->eulerianSource, "analytical", sizeof(simCtx->eulerianSource)));
    PetscCall(PetscStrncpy(simCtx->AnalyticalSolutionType, "ZERO_FLOW", sizeof(simCtx->AnalyticalSolutionType)));
    PetscCall(CaptureBannerOutput(simCtx, captured, sizeof(captured)));
    PetscCall(AssertCapturedContains(captured, "Analytical Solution Type : ZERO_FLOW",
                                     "DisplayBanner should include the analytical solution type for analytical runs"));
    PetscCall(AssertCapturedContains(captured, "Particle Console Cadence   : every 6 step(s)",
                                     "DisplayBanner should include active particle console cadence when particles are configured"));
    PetscCall(AssertCapturedContains(captured, "Particle Initialization Mode: Surface: Random",
                                     "DisplayBanner should include particle initialization mode for analytical particle runs"));
    PetscCall(AssertCapturedContains(captured, "Particles Initialized At",
                                     "DisplayBanner should include inlet-face placement details for surface particle initialization"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Runs the unit-io PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"should-write-data-output", TestShouldWriteDataOutput},
        {"verify-path-existence", TestVerifyPathExistence},
        {"write-and-read-simulation-fields", TestWriteAndReadSimulationFields},
        {"checkpoint-same-step-rewrite-rejected", TestCheckpointSameStepRewriteIsRejected},
        {"display-banner-reports-statistics-cadence", TestDisplayBannerReportsStatisticsCadence},
        {"checkpoint-statistics-round-trip", TestCheckpointStatisticsRoundTrip},
        {"checkpoint-statistics-continuation-guards", TestCheckpointStatisticsContinuationGuards},
        {"checkpoint-statistics-payload-is-validated", TestCheckpointStatisticsPayloadIsValidated},
        {"checkpoint-statistics-absent-when-disabled", TestCheckpointStatisticsAbsentWhenDisabled},
        {"checkpoint-sha256-known-vector", TestCheckpointSHA256KnownVector},
        {"parse-post-processing-settings", TestParsePostProcessingSettings},
        {"trim-whitespace", TestTrimWhitespace},
        {"bc-string-parsers", TestBoundaryConditionStringParsers},
        {"validate-bc-handler-for-type", TestValidateBCHandlerForBCType},
        {"parse-scaling-information", TestParseScalingInformation},
        {"display-banner-startup-summary", TestDisplayBannerTracksConditionalStartupFields},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv I/O tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-io", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
