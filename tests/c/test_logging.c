/**
 * @file test_logging.c
 * @brief C unit tests for runtime log-level, allow-list, conversion, and profiling helpers.
 */

#include "test_support.h"

#include "logging.h"

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
/**
 * @brief Asserts that one text file contains a required substring.
 */

static PetscErrorCode AssertFileContains(const char *path, const char *needle, const char *context)
{
    FILE *fp = NULL;
    long file_size = 0;
    char *buffer = NULL;

    PetscFunctionBeginUser;
    fp = fopen(path, "rb");
    if (!fp) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to open '%s' for assertion.", path);
    }
    if (fseek(fp, 0, SEEK_END) != 0) {
        fclose(fp);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to seek '%s'.", path);
    }
    file_size = ftell(fp);
    if (file_size < 0) {
        fclose(fp);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to measure '%s'.", path);
    }
    if (fseek(fp, 0, SEEK_SET) != 0) {
        fclose(fp);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to rewind '%s'.", path);
    }

    PetscCall(PetscMalloc1((size_t)file_size + 1, &buffer));
    if (file_size > 0) {
        size_t bytes_read = fread(buffer, 1, (size_t)file_size, fp);
        if (bytes_read != (size_t)file_size) {
            fclose(fp);
            PetscCall(PetscFree(buffer));
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to read '%s'.", path);
        }
    }
    buffer[file_size] = '\0';
    fclose(fp);

    PetscCall(PicurvAssertBool((PetscBool)(strstr(buffer, needle) != NULL), context));
    PetscCall(PetscFree(buffer));
    PetscFunctionReturn(0);
}
/**
 * @brief Asserts that one text file does not contain an excluded substring.
 */

static PetscErrorCode AssertFileNotContains(const char *path, const char *needle, const char *context)
{
    FILE *fp = NULL;
    long file_size = 0;
    char *buffer = NULL;

    PetscFunctionBeginUser;
    fp = fopen(path, "rb");
    if (!fp) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to open '%s' for assertion.", path);
    }
    if (fseek(fp, 0, SEEK_END) != 0) {
        fclose(fp);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to seek '%s'.", path);
    }
    file_size = ftell(fp);
    if (file_size < 0) {
        fclose(fp);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to measure '%s'.", path);
    }
    if (fseek(fp, 0, SEEK_SET) != 0) {
        fclose(fp);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to rewind '%s'.", path);
    }

    PetscCall(PetscMalloc1((size_t)file_size + 1, &buffer));
    if (file_size > 0) {
        size_t bytes_read = fread(buffer, 1, (size_t)file_size, fp);
        if (bytes_read != (size_t)file_size) {
            fclose(fp);
            PetscCall(PetscFree(buffer));
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to read '%s'.", path);
        }
    }
    buffer[file_size] = '\0';
    fclose(fp);

    PetscCall(PicurvAssertBool((PetscBool)(strstr(buffer, needle) == NULL), context));
    PetscCall(PetscFree(buffer));
    PetscFunctionReturn(0);
}

typedef PetscErrorCode (*CapturedLoggingFn)(UserCtx *user, SimCtx *simCtx, void *ctx);

typedef struct AnatomyCaptureCtx {
    const char *field_name;
    const char *stage_name;
} AnatomyCaptureCtx;

/**
 * @brief Captures stdout emitted by one logging helper into a temporary file-backed buffer.
 */
static PetscErrorCode CaptureLoggingOutput(UserCtx *user,
                                           SimCtx *simCtx,
                                           CapturedLoggingFn fn,
                                           void *ctx,
                                           char *captured,
                                           size_t captured_len)
{
    char tmpdir[PETSC_MAX_PATH_LEN];
    char capture_path[PETSC_MAX_PATH_LEN];
    FILE *capture_file = NULL;
    int saved_stdout = -1;
    int capture_fd = -1;
    size_t bytes_read = 0;
    PetscErrorCode ierr;

    PetscFunctionBeginUser;
    PetscCheck(fn != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Capture callback cannot be NULL.");
    PetscCheck(captured != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Capture buffer cannot be NULL.");
    PetscCheck(captured_len > 0, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, "Capture buffer must be non-empty.");

    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(capture_path, sizeof(capture_path), "%s/logging.out", tmpdir));

    fflush(stdout);
    saved_stdout = dup(STDOUT_FILENO);
    PetscCheck(saved_stdout >= 0, PETSC_COMM_SELF, PETSC_ERR_SYS, "dup(STDOUT_FILENO) failed.");
    capture_fd = open(capture_path, O_CREAT | O_TRUNC | O_WRONLY, 0600);
    PetscCheck(capture_fd >= 0, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to open capture file '%s'.", capture_path);
    PetscCheck(dup2(capture_fd, STDOUT_FILENO) >= 0, PETSC_COMM_SELF, PETSC_ERR_SYS, "dup2() failed while redirecting stdout.");
    close(capture_fd);
    capture_fd = -1;

    ierr = fn(user, simCtx, ctx);
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
    PetscFunctionReturn(0);
}

/**
 * @brief Adapts `LOG_PARTICLE_FIELDS()` to the generic stdout-capture callback shape.
 */
static PetscErrorCode InvokeParticleFieldLog(UserCtx *user, SimCtx *simCtx, void *ctx)
{
    PetscInt print_interval = *((PetscInt *)ctx);

    PetscFunctionBeginUser;
    (void)simCtx;
    PetscCall(LOG_PARTICLE_FIELDS(user, print_interval));
    PetscFunctionReturn(0);
}

/**
 * @brief Adapts `EmitParticleConsoleSnapshot()` to the generic stdout-capture callback shape.
 */
static PetscErrorCode InvokeParticleConsoleSnapshot(UserCtx *user, SimCtx *simCtx, void *ctx)
{
    PetscInt step = *((PetscInt *)ctx);

    PetscFunctionBeginUser;
    PetscCall(EmitParticleConsoleSnapshot(user, simCtx, step));
    PetscFunctionReturn(0);
}

/**
 * @brief Adapts `LOG_FIELD_ANATOMY()` to the generic stdout-capture callback shape.
 */
static PetscErrorCode InvokeFieldAnatomyLog(UserCtx *user, SimCtx *simCtx, void *ctx)
{
    const AnatomyCaptureCtx *anatomy_ctx = (const AnatomyCaptureCtx *)ctx;

    PetscFunctionBeginUser;
    (void)simCtx;
    PetscCall(LOG_FIELD_ANATOMY(user, anatomy_ctx->field_name, anatomy_ctx->stage_name));
    PetscFunctionReturn(0);
}

/**
 * @brief Creates a small particle-bearing runtime fixture used by logging tests.
 */
static PetscErrorCode SeedLoggingParticleFixture(SimCtx **simCtx_out, UserCtx **user_out)
{
    PetscReal *positions = NULL;
    PetscReal *velocities = NULL;
    PetscReal *weights = NULL;
    PetscInt *cell_ids = NULL;
    PetscInt *status = NULL;
    PetscReal *psi = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(simCtx_out, user_out, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(*user_out, 2, "ske"));
    (*simCtx_out)->np = 2;
    (*simCtx_out)->particleConsoleOutputFreq = 2;
    (*simCtx_out)->LoggingFrequency = 1;

    PetscCall(DMSwarmGetField((*user_out)->swarm, "position", NULL, NULL, (void **)&positions));
    PetscCall(DMSwarmGetField((*user_out)->swarm, "velocity", NULL, NULL, (void **)&velocities));
    PetscCall(DMSwarmGetField((*user_out)->swarm, "weight", NULL, NULL, (void **)&weights));
    PetscCall(DMSwarmGetField((*user_out)->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    PetscCall(DMSwarmGetField((*user_out)->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(DMSwarmGetField((*user_out)->swarm, "Psi", NULL, NULL, (void **)&psi));

    positions[0] = 0.25; positions[1] = 0.50; positions[2] = 0.75;
    positions[3] = 0.50; positions[4] = 0.50; positions[5] = 0.50;
    velocities[0] = 1.0; velocities[1] = 2.0; velocities[2] = 3.0;
    velocities[3] = 4.0; velocities[4] = 5.0; velocities[5] = 6.0;
    weights[0] = 0.2; weights[1] = 0.3; weights[2] = 0.4;
    weights[3] = 0.5; weights[4] = 0.5; weights[5] = 0.5;
    cell_ids[0] = 0; cell_ids[1] = 0; cell_ids[2] = 0;
    cell_ids[3] = 1; cell_ids[4] = 1; cell_ids[5] = 1;
    status[0] = ACTIVE_AND_LOCATED;
    status[1] = ACTIVE_AND_LOCATED;
    psi[0] = 1.0;
    psi[1] = 3.0;

    PetscCall(DMSwarmRestoreField((*user_out)->swarm, "Psi", NULL, NULL, (void **)&psi));
    PetscCall(DMSwarmRestoreField((*user_out)->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(DMSwarmRestoreField((*user_out)->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    PetscCall(DMSwarmRestoreField((*user_out)->swarm, "weight", NULL, NULL, (void **)&weights));
    PetscCall(DMSwarmRestoreField((*user_out)->swarm, "velocity", NULL, NULL, (void **)&velocities));
    PetscCall(DMSwarmRestoreField((*user_out)->swarm, "position", NULL, NULL, (void **)&positions));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests string-conversion helpers for configured enums and unknown values.
 */

static PetscErrorCode TestStringConversionHelpers(void)
{
    PetscFunctionBeginUser;
    PetscCall(PicurvAssertBool(strcmp(BCFaceToString(BC_FACE_NEG_X), "-Xi (I-Min)") == 0,
                               "BCFaceToString should report the negative-x face"));
    PetscCall(PicurvAssertBool(strcmp(FieldInitializationToString(0), "Zero") == 0,
                               "FieldInitializationToString should report the zero mode"));
    PetscCall(PicurvAssertBool(strcmp(FieldInitializationToString(99), "Unknown Field Initialization") == 0,
                               "FieldInitializationToString should reject unknown selectors"));
    PetscCall(PicurvAssertBool(strcmp(ParticleInitializationToString(PARTICLE_INIT_VOLUME), "Volume") == 0,
                               "ParticleInitializationToString should report the volume mode"));
    PetscCall(PicurvAssertBool(strcmp(LESModelToString(CONSTANT_SMAGORINSKY), "Constant Smagorinsky") == 0,
                               "LESModelToString should report the constant model"));
    PetscCall(PicurvAssertBool(strcmp(MomentumSolverTypeToString(MOMENTUM_SOLVER_EXPLICIT_RK), "Explicit 4 stage Runge-Kutta ") == 0,
                               "MomentumSolverTypeToString should report the explicit solver"));
    PetscCall(PicurvAssertBool(strcmp(BCTypeToString(PERIODIC), "PERIODIC") == 0,
                               "BCTypeToString should report periodic boundaries"));
    PetscCall(PicurvAssertBool(strcmp(BCHandlerTypeToString(BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX), "constant flux") == 0,
                               "BCHandlerTypeToString should report the driven periodic handler"));
    PetscCall(PicurvAssertBool(strcmp(ParticleLocationStatusToString(LOST), "LOST") == 0,
                               "ParticleLocationStatusToString should report LOST state"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that log level selection honors the environment variable.
 */

static PetscErrorCode TestGetLogLevelFromEnvironment(void)
{
    PetscFunctionBeginUser;
    PetscCall(PicurvAssertIntEqual(LOG_INFO, get_log_level(),
                                   "get_log_level should honor LOG_LEVEL=INFO in this test binary"));
    PetscCall(print_log_level());
    PetscFunctionReturn(0);
}
/**
 * @brief Tests the function allow-list filter used by the logging layer.
 */

static PetscErrorCode TestAllowedFunctionsFilter(void)
{
    const char *allow_list[] = {"ComputeSpecificKE", "WriteEulerianFile"};

    PetscFunctionBeginUser;
    set_allowed_functions(allow_list, 2);
    PetscCall(PicurvAssertBool(is_function_allowed("ComputeSpecificKE"),
                               "Allowed list should include ComputeSpecificKE"));
    PetscCall(PicurvAssertBool((PetscBool)!is_function_allowed("UnlistedFunction"),
                               "Allowed list should exclude unknown function names"));

    set_allowed_functions(NULL, 0);
    PetscCall(PicurvAssertBool(is_function_allowed("AnyFunction"),
                               "Empty allow-list should permit all functions"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests periodic particle console snapshot enablement and cadence.
 */

static PetscErrorCode TestParticleConsoleSnapshotCadence(void)
{
    SimCtx simCtx;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&simCtx, sizeof(simCtx)));
    simCtx.np = 32;
    simCtx.particleConsoleOutputFreq = 4;

    PetscCall(PicurvAssertBool(IsParticleConsoleSnapshotEnabled(&simCtx),
                               "Particle snapshot contract should be enabled when particles and cadence are configured"));
    PetscCall(PicurvAssertBool(ShouldEmitPeriodicParticleConsoleSnapshot(&simCtx, 8),
                               "Snapshot should emit on cadence-aligned completed steps"));
    PetscCall(PicurvAssertBool((PetscBool)!ShouldEmitPeriodicParticleConsoleSnapshot(&simCtx, 7),
                               "Snapshot should not emit off-cadence"));

    simCtx.particleConsoleOutputFreq = 0;
    PetscCall(PicurvAssertBool((PetscBool)!IsParticleConsoleSnapshotEnabled(&simCtx),
                               "Zero cadence should disable periodic particle snapshots"));
    PetscCall(PicurvAssertBool((PetscBool)!ShouldEmitPeriodicParticleConsoleSnapshot(NULL, 4),
                               "NULL SimCtx should never emit periodic snapshots"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests logging-side file parsing, helper formatting, and progress utilities.
 */

static PetscErrorCode TestLoggingFileParsingAndFormattingHelpers(void)
{
    char tmpdir[PETSC_MAX_PATH_LEN];
    char allow_path[PETSC_MAX_PATH_LEN];
    char dual_log_path[PETSC_MAX_PATH_LEN];
    FILE *file = NULL;
    char **funcs = NULL;
    PetscInt nfuncs = 0;
    Cell cell;
    PetscReal distances[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    DualMonitorCtx *monctx = NULL;
    void *ctx = NULL;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&cell, sizeof(cell)));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(allow_path, sizeof(allow_path), "%s/allowed_functions.txt", tmpdir));
    PetscCall(PetscSNPrintf(dual_log_path, sizeof(dual_log_path), "%s/dual-monitor.log", tmpdir));

    file = fopen(allow_path, "w");
    PetscCheck(file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to create allow-list file '%s'.", allow_path);
    fputs("  ComputeSpecificKE  \n", file);
    fputs("# comment-only line\n", file);
    for (PetscInt i = 0; i < 17; ++i) {
        fprintf(file, "Helper_%02d   # trailing comment\n", (int)i);
    }
    fclose(file);
    file = NULL;

    PetscCall(LoadAllowedFunctionsFromFile(allow_path, &funcs, &nfuncs));
    PetscCall(PicurvAssertIntEqual(18, nfuncs, "LoadAllowedFunctionsFromFile should trim comments and keep all identifiers"));
    PetscCall(PicurvAssertBool((PetscBool)(strcmp(funcs[0], "ComputeSpecificKE") == 0),
                               "LoadAllowedFunctionsFromFile should trim leading and trailing whitespace"));
    PetscCall(PicurvAssertBool((PetscBool)(strcmp(funcs[17], "Helper_16") == 0),
                               "LoadAllowedFunctionsFromFile should grow past the initial pointer capacity"));
    PetscCall(FreeAllowedFunctions(funcs, nfuncs));

    for (PetscInt i = 0; i < 8; ++i) {
        cell.vertices[i].x = (PetscReal)i;
        cell.vertices[i].y = (PetscReal)(i + 1);
        cell.vertices[i].z = (PetscReal)(i + 2);
    }
    PetscCall(LOG_CELL_VERTICES(&cell, 0));
    PetscCall(LOG_FACE_DISTANCES(distances));

    PetscCall(PetscCalloc1(1, &monctx));
    monctx->file_handle = fopen(dual_log_path, "w");
    PetscCheck(monctx->file_handle != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to create dual-monitor log '%s'.", dual_log_path);
    ctx = monctx;
    PetscCall(DualMonitorDestroy(&ctx));
    PetscCall(PicurvAssertBool((PetscBool)(ctx == NULL),
                               "DualMonitorDestroy should clear the caller-owned context pointer"));

    PrintProgressBar(0, 0, 4, 0.10);
    PrintProgressBar(3, 0, 4, 0.40);
    PrintProgressBar(0, 0, 0, 0.00);
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "\n"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests continuity, min/max, and anatomy logging helpers on minimal runtime fixtures.
 */

static PetscErrorCode TestLoggingContinuityAndFieldDiagnostics(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char continuity_path[PETSC_MAX_PATH_LEN];
    PetscReal ***p = NULL;
    Cmpnts ***ucat = NULL;
    Cmpnts ***ucont = NULL;
    PetscErrorCode ierr_minmax = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscStrncpy(simCtx->log_dir, tmpdir, sizeof(simCtx->log_dir)));

    simCtx->StartStep = 0;
    simCtx->step = 1;
    simCtx->MaxDiv = 1.25;
    simCtx->MaxDivx = 1;
    simCtx->MaxDivy = 2;
    simCtx->MaxDivz = 3;
    simCtx->MaxDivFlatArg = 17;
    simCtx->summationRHS = 8.5;
    simCtx->FluxInSum = 5.0;
    simCtx->FluxOutSum = 3.25;
    PetscCall(LOG_CONTINUITY_METRICS(user));

    simCtx->step = 2;
    simCtx->MaxDiv = 0.75;
    simCtx->summationRHS = 4.5;
    simCtx->FluxInSum = 2.5;
    simCtx->FluxOutSum = 1.0;
    PetscCall(LOG_CONTINUITY_METRICS(user));

    PetscCall(PetscSNPrintf(continuity_path, sizeof(continuity_path), "%s/Continuity_Metrics.log", simCtx->log_dir));
    PetscCall(PicurvAssertFileExists(continuity_path, "continuity metrics log should be written"));
    PetscCall(AssertFileContains(continuity_path, "Timestep", "continuity metrics log should include the header"));
    PetscCall(AssertFileContains(continuity_path, "([3][2][1] = 17)", "continuity metrics log should include the divergence location"));
    PetscCall(AssertFileContains(continuity_path, "2          | 0", "continuity metrics log should append later timesteps"));

    PetscCall(DMDAVecGetArray(user->da, user->P, &p));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucont, &ucont));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
                p[k][j][i] = (PetscReal)(i + j + k);
                ucat[k][j][i].x = (PetscReal)i;
                ucat[k][j][i].y = (PetscReal)(-j);
                ucat[k][j][i].z = (PetscReal)(2 * k);
                ucont[k][j][i].x = (PetscReal)(10 + i);
                ucont[k][j][i].y = (PetscReal)(20 + j);
                ucont[k][j][i].z = (PetscReal)(30 + k);
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucont, &ucont));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));
    PetscCall(DMDAVecRestoreArray(user->da, user->P, &p));

    PetscCall(DMGlobalToLocalBegin(user->da, user->P, INSERT_VALUES, user->lP));
    PetscCall(DMGlobalToLocalEnd(user->da, user->P, INSERT_VALUES, user->lP));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat));
    PetscCall(DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));
    PetscCall(DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont));

    PetscCall(LOG_FIELD_MIN_MAX(user, "P"));
    PetscCall(LOG_FIELD_MIN_MAX(user, "Ucat"));
    PetscCall(LOG_FIELD_MIN_MAX(user, "Coordinates"));
    PetscCall(LOG_FIELD_MIN_MAX(user, "Ucont"));

    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    ierr_minmax = LOG_FIELD_MIN_MAX(user, "NotARealField");
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertBool((PetscBool)(ierr_minmax != 0),
                               "LOG_FIELD_MIN_MAX should reject unknown field names"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests interpolation-error logging against an analytically matched particle field.
 */

static PetscErrorCode TestInterpolationErrorLogging(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal (*pos_arr)[3] = NULL;
    PetscReal (*vel_arr)[3] = NULL;
    Vec position_vec = NULL;
    Vec analytical_vec = NULL;
    const PetscScalar *analytical_arr = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 2, "ske"));
    PetscCall(PetscStrncpy(simCtx->AnalyticalSolutionType, "TGV3D", sizeof(simCtx->AnalyticalSolutionType)));
    simCtx->ren = 1.0;
    simCtx->ti = 0.0;

    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));
    pos_arr[0][0] = 0.5 * PETSC_PI; pos_arr[0][1] = 0.0;          pos_arr[0][2] = 0.0;
    pos_arr[1][0] = 0.0;           pos_arr[1][1] = 0.5 * PETSC_PI; pos_arr[1][2] = 0.0;
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));

    PetscCall(DMSwarmCreateGlobalVectorFromField(user->swarm, "position", &position_vec));
    PetscCall(VecDuplicate(position_vec, &analytical_vec));
    PetscCall(VecCopy(position_vec, analytical_vec));
    PetscCall(SetAnalyticalSolutionForParticles(analytical_vec, simCtx));

    PetscCall(DMSwarmGetField(user->swarm, "velocity", NULL, NULL, (void *)&vel_arr));
    PetscCall(VecGetArrayRead(analytical_vec, &analytical_arr));
    for (PetscInt particle = 0; particle < 2; ++particle) {
        vel_arr[particle][0] = PetscRealPart(analytical_arr[3 * particle + 0]);
        vel_arr[particle][1] = PetscRealPart(analytical_arr[3 * particle + 1]);
        vel_arr[particle][2] = PetscRealPart(analytical_arr[3 * particle + 2]);
    }
    PetscCall(VecRestoreArrayRead(analytical_vec, &analytical_arr));
    PetscCall(DMSwarmRestoreField(user->swarm, "velocity", NULL, NULL, (void *)&vel_arr));

    PetscCall(VecDestroy(&analytical_vec));
    PetscCall(DMSwarmDestroyGlobalVectorFromField(user->swarm, "position", &position_vec));
    PetscCall(LOG_INTERPOLATION_ERROR(user));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests stdout particle-table logging on a production-like swarm fixture.
 */

static PetscErrorCode TestParticleFieldTableLogging(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscInt print_interval = 1;
    char captured[8192];

    PetscFunctionBeginUser;
    PetscCall(SeedLoggingParticleFixture(&simCtx, &user));
    PetscCall(CaptureLoggingOutput(user, simCtx, InvokeParticleFieldLog, &print_interval, captured, sizeof(captured)));
    PetscCall(PicurvAssertBool((PetscBool)(strstr(captured, "Position (x,y,z)") != NULL),
                               "LOG_PARTICLE_FIELDS should print the particle table header"));
    PetscCall(PicurvAssertBool((PetscBool)(strstr(captured, "Weights (a1,a2,a3)") != NULL),
                               "LOG_PARTICLE_FIELDS should print the weight-column header"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests console snapshot logging against the public periodic-snapshot helper.
 */

static PetscErrorCode TestParticleConsoleSnapshotLogging(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscInt step = 4;
    char captured[8192];

    PetscFunctionBeginUser;
    PetscCall(SeedLoggingParticleFixture(&simCtx, &user));
    PetscCall(CaptureLoggingOutput(user, simCtx, InvokeParticleConsoleSnapshot, &step, captured, sizeof(captured)));
    PetscCall(PicurvAssertBool((PetscBool)(strstr(captured, "Particle states at step 4") != NULL),
                               "EmitParticleConsoleSnapshot should print the step banner"));
    PetscCall(PicurvAssertBool((PetscBool)(strstr(captured, "Position (x,y,z)") != NULL),
                               "EmitParticleConsoleSnapshot should reuse the particle table output"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests file-backed particle metrics logging after derived metrics are computed.
 */

static PetscErrorCode TestParticleMetricsLogging(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char metrics_path[PETSC_MAX_PATH_LEN];

    PetscFunctionBeginUser;
    PetscCall(SeedLoggingParticleFixture(&simCtx, &user));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscStrncpy(simCtx->log_dir, tmpdir, sizeof(simCtx->log_dir)));
    simCtx->StartStep = 0;
    simCtx->step = 1;
    simCtx->particlesLostLastStep = 1;
    simCtx->particlesMigratedLastStep = 2;
    simCtx->migrationPassesLastStep = 3;

    PetscCall(CalculateParticleCountPerCell(user));
    PetscCall(CalculateAdvancedParticleMetrics(user));
    PetscCall(LOG_PARTICLE_METRICS(user, "Timestep Metrics"));

    PetscCall(PetscSNPrintf(metrics_path, sizeof(metrics_path), "%s/Particle_Metrics.log", simCtx->log_dir));
    PetscCall(PicurvAssertFileExists(metrics_path, "LOG_PARTICLE_METRICS should write Particle_Metrics.log"));
    PetscCall(AssertFileContains(metrics_path, "Timestep Metrics", "Particle metrics log should include the caller-provided stage label"));
    PetscCall(AssertFileContains(metrics_path, "Occupied Cells", "Particle metrics log should include the metrics table header"));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests stdout field-anatomy logging on the corrected production-like DM fixture.
 */

static PetscErrorCode TestFieldAnatomyLogging(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char captured[8192];
    AnatomyCaptureCtx anatomy_ctx = {"P", "unit-test"};

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(VecSet(user->P, 7.0));
    PetscCall(DMGlobalToLocalBegin(user->da, user->P, INSERT_VALUES, user->lP));
    PetscCall(DMGlobalToLocalEnd(user->da, user->P, INSERT_VALUES, user->lP));
    PetscCall(CaptureLoggingOutput(user, simCtx, InvokeFieldAnatomyLog, &anatomy_ctx, captured, sizeof(captured)));
    PetscCall(PicurvAssertBool((PetscBool)(strstr(captured, "Field Anatomy Log: [P]") != NULL),
                               "LOG_FIELD_ANATOMY should print the requested field name"));
    PetscCall(PicurvAssertBool((PetscBool)(strstr(captured, "Layout: [Cell-Centered]") != NULL),
                               "LOG_FIELD_ANATOMY should report the inferred data layout"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests profiling helper lifecycle logging for timestep and final-summary outputs.
 */

static PetscErrorCode TestProfilingLifecycleHelpers(void)
{
    SimCtx simCtx;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char timestep_path[PETSC_MAX_PATH_LEN];
    char summary_path[PETSC_MAX_PATH_LEN];
    static char selected_name[] = "FlowSolver";
    char *selected_funcs[] = {selected_name};

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&simCtx, sizeof(simCtx)));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscStrncpy(simCtx.log_dir, tmpdir, sizeof(simCtx.log_dir)));
    PetscCall(PetscStrncpy(simCtx.profilingTimestepMode, "selected", sizeof(simCtx.profilingTimestepMode)));
    PetscCall(PetscStrncpy(simCtx.profilingTimestepFile, "Profiling_Timestep_Summary.csv", sizeof(simCtx.profilingTimestepFile)));
    simCtx.rank = 0;
    simCtx.exec_mode = EXEC_MODE_SOLVER;
    simCtx.StartStep = 0;
    simCtx.nProfilingSelectedFuncs = 1;
    simCtx.profilingSelectedFuncs = selected_funcs;
    simCtx.profilingFinalSummary = PETSC_TRUE;

    PetscCall(ProfilingInitialize(&simCtx));

    _ProfilingStart("FlowSolver");
    _ProfilingEnd("FlowSolver");
    _ProfilingStart("UnselectedHelper");
    _ProfilingEnd("UnselectedHelper");
    PetscCall(ProfilingLogTimestepSummary(&simCtx, 1));

    _ProfilingStart("FlowSolver");
    _ProfilingEnd("FlowSolver");
    PetscCall(ProfilingResetTimestepCounters());
    PetscCall(ProfilingLogTimestepSummary(&simCtx, 2));
    PetscCall(ProfilingFinalize(&simCtx));

    PetscCall(PetscSNPrintf(timestep_path, sizeof(timestep_path), "%s/%s", simCtx.log_dir, simCtx.profilingTimestepFile));
    PetscCall(PetscSNPrintf(summary_path, sizeof(summary_path), "%s/ProfilingSummary_Solver.log", simCtx.log_dir));
    PetscCall(PicurvAssertFileExists(timestep_path, "profiling timestep summary should be written"));
    PetscCall(PicurvAssertFileExists(summary_path, "profiling final summary should be written"));
    PetscCall(AssertFileContains(timestep_path, "step,function,calls,step_time_s",
                                 "profiling timestep summary should contain the CSV header"));
    PetscCall(AssertFileContains(timestep_path, "1,FlowSolver,1,",
                                 "profiling timestep summary should log selected functions"));
    PetscCall(AssertFileNotContains(timestep_path, "UnselectedHelper",
                                    "profiling timestep summary should omit unselected functions in selected mode"));
    PetscCall(AssertFileContains(summary_path, "FINAL PROFILING SUMMARY",
                                 "profiling final summary should include its table banner"));
    PetscCall(AssertFileContains(summary_path, "FlowSolver",
                                 "profiling final summary should include selected functions"));
    PetscCall(AssertFileContains(summary_path, "UnselectedHelper",
                                 "profiling final summary should include total-time entries for unselected functions"));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscFunctionReturn(0);
}
/**
 * @brief Runs the unit-logging PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"string-conversion-helpers", TestStringConversionHelpers},
        {"get-log-level-from-environment", TestGetLogLevelFromEnvironment},
        {"allowed-functions-filter", TestAllowedFunctionsFilter},
        {"particle-console-snapshot-cadence", TestParticleConsoleSnapshotCadence},
        {"logging-file-parsing-and-formatting-helpers", TestLoggingFileParsingAndFormattingHelpers},
        {"logging-continuity-and-field-diagnostics", TestLoggingContinuityAndFieldDiagnostics},
        {"interpolation-error-logging", TestInterpolationErrorLogging},
        {"particle-field-table-logging", TestParticleFieldTableLogging},
        {"particle-console-snapshot-logging", TestParticleConsoleSnapshotLogging},
        {"particle-metrics-logging", TestParticleMetricsLogging},
        {"field-anatomy-logging", TestFieldAnatomyLogging},
        {"profiling-lifecycle-helpers", TestProfilingLifecycleHelpers},
    };

    (void)setenv("LOG_LEVEL", "INFO", 1);

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv logging tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-logging", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    set_allowed_functions(NULL, 0);

    ierr = PetscFinalize();
    return (int)ierr;
}
