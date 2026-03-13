/**
 * @file test_setup_lifecycle.c
 * @brief C unit tests for setup, initialization, and cleanup lifecycle entry points.
 */

#include "test_support.h"

#include "ParticleSwarm.h"
#include "initialcondition.h"
#include "runloop.h"
#include "setup.h"

#include <stdio.h>
#include <string.h>

/**
 * @brief Asserts that a directory path exists and is readable.
 */
static PetscErrorCode AssertDirectoryExists(const char *path, const char *context)
{
    PetscBool exists = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(PetscTestDirectory(path, 'r', &exists));
    PetscCall(PicurvAssertBool(exists, context));
    PetscFunctionReturn(0);
}

/**
 * @brief Writes one small temporary text file used by the partial-lifecycle context-only fixture.
 */
static PetscErrorCode WriteContextOnlyFile(const char *path, const char *contents)
{
    FILE *file = NULL;

    PetscFunctionBeginUser;
    file = fopen(path, "w");
    PetscCheck(file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to open '%s' for writing.", path);
    fputs(contents, file);
    fclose(file);
    PetscFunctionReturn(0);
}

/**
 * @brief Creates the control-file bundle needed for the context-only cleanup test.
 */
static PetscErrorCode PrepareContextOnlyConfig(char *tmpdir,
                                               size_t tmpdir_len,
                                               char *control_path,
                                               size_t control_path_len)
{
    char bcs_path[PETSC_MAX_PATH_LEN];
    char post_path[PETSC_MAX_PATH_LEN];
    char output_dir[PETSC_MAX_PATH_LEN];
    char log_dir[PETSC_MAX_PATH_LEN];
    char control_buffer[8192];

    PetscFunctionBeginUser;
    PetscCall(PicurvMakeTempDir(tmpdir, tmpdir_len));
    PetscCall(PetscSNPrintf(bcs_path, sizeof(bcs_path), "%s/bcs.run", tmpdir));
    PetscCall(PetscSNPrintf(post_path, sizeof(post_path), "%s/post.run", tmpdir));
    PetscCall(PetscSNPrintf(output_dir, sizeof(output_dir), "%s/results", tmpdir));
    PetscCall(PetscSNPrintf(log_dir, sizeof(log_dir), "%s/logs", tmpdir));
    PetscCall(PetscSNPrintf(control_path, control_path_len, "%s/test.control", tmpdir));

    PetscCall(WriteContextOnlyFile(
        bcs_path,
        "-Xi WALL noslip\n"
        "+Xi WALL noslip\n"
        "-Eta WALL noslip\n"
        "+Eta WALL noslip\n"
        "-Zeta INLET constant_velocity vx=0.0 vy=0.0 vz=1.5\n"
        "+Zeta OUTLET conservation\n"));
    PetscCall(WriteContextOnlyFile(
        post_path,
        "startTime = 0\n"
        "endTime = 1\n"
        "timeStep = 1\n"
        "output_particles = false\n"));
    PetscCall(PetscSNPrintf(
        control_buffer,
        sizeof(control_buffer),
        "-start_step 0\n"
        "-totalsteps 2\n"
        "-ren 100.0\n"
        "-dt 0.001\n"
        "-finit 1\n"
        "-ucont_x 0.0\n"
        "-ucont_y 0.0\n"
        "-ucont_z 1.5\n"
        "-bcs_files %s\n"
        "-profiling_timestep_mode off\n"
        "-profiling_final_summary true\n"
        "-postprocessing_config_file %s\n"
        "-grid\n"
        "-im 6\n"
        "-jm 6\n"
        "-km 6\n"
        "-xMins 0.0\n"
        "-xMaxs 1.0\n"
        "-yMins 0.0\n"
        "-yMaxs 1.0\n"
        "-zMins 0.0\n"
        "-zMaxs 1.0\n"
        "-rxs 1.0\n"
        "-rys 1.0\n"
        "-rzs 1.0\n"
        "-cgrids 0\n"
        "-nblk 1\n"
        "-euler_field_source solve\n"
        "-mom_solver_type EXPLICIT_RK\n"
        "-mg_level 1\n"
        "-poisson 0\n"
        "-tio 0\n"
        "-particle_console_output_freq 0\n"
        "-logfreq 1\n"
        "-output_dir %s\n"
        "-restart_dir %s\n"
        "-log_dir %s\n"
        "-numParticles 0\n"
        "-pinit 2\n",
        bcs_path,
        post_path,
        output_dir,
        output_dir,
        log_dir));
    PetscCall(WriteContextOnlyFile(control_path, control_buffer));
    PetscFunctionReturn(0);
}

/**
 * @brief Finalizes and frees one lifecycle test context, then clears any PETSc options used to build it.
 */
static PetscErrorCode FreeLifecycleContext(SimCtx **simCtx_ptr)
{
    PetscFunctionBeginUser;
    PetscCall(PicurvDestroyRuntimeContext(simCtx_ptr));
    PetscFunctionReturn(0);
}

/**
 * @brief Builds a full setup fixture through environment, grid, BC, and domain-rank initialization.
 */
static PetscErrorCode BuildLifecycleContext(PetscBool enable_particles, SimCtx **simCtx_out, char *tmpdir, size_t tmpdir_len)
{
    PetscFunctionBeginUser;
    PetscCall(PicurvBuildTinyRuntimeContext(NULL, enable_particles, simCtx_out, NULL, tmpdir, tmpdir_len));
    PetscFunctionReturn(0);
}

/**
 * @brief Builds only the top-level simulation context used by partial-initialization cleanup tests.
 */
static PetscErrorCode BuildContextOnly(SimCtx **simCtx_out, char *tmpdir, size_t tmpdir_len)
{
    char control_path[PETSC_MAX_PATH_LEN];
    SimCtx *simCtx = NULL;

    PetscFunctionBeginUser;
    PetscCall(PetscOptionsClear(NULL));
    PetscCall(PrepareContextOnlyConfig(tmpdir, tmpdir_len, control_path, sizeof(control_path)));
    PetscCall(PetscOptionsSetValue(NULL, "-control_file", control_path));
    PetscCall(CreateSimulationContext(0, NULL, &simCtx));
    simCtx->exec_mode = EXEC_MODE_SOLVER;
    *simCtx_out = simCtx;
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that the shared richer runtime fixture mirrors normalized production setup contracts.
 */
static PetscErrorCode TestSharedRuntimeFixtureContracts(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];

    PetscFunctionBeginUser;
    PetscCall(PicurvBuildTinyRuntimeContext(NULL, PETSC_FALSE, &simCtx, &user, tmpdir, sizeof(tmpdir)));

    PetscCall(PicurvAssertBool((PetscBool)(user->bbox.min_coords.x >= -2.0e-6), "runtime fixture bbox xmin should stay inside normalized domain tolerance"));
    PetscCall(PicurvAssertBool((PetscBool)(user->bbox.max_coords.x <= 1.0 + 2.0e-6), "runtime fixture bbox xmax should stay inside normalized domain tolerance"));
    PetscCall(PicurvAssertBool((PetscBool)(user->bbox.max_coords.y <= 1.0 + 2.0e-6), "runtime fixture bbox ymax should stay inside normalized domain tolerance"));
    PetscCall(PicurvAssertBool((PetscBool)(user->bbox.max_coords.z <= 1.0 + 2.0e-6), "runtime fixture bbox zmax should stay inside normalized domain tolerance"));
    PetscCall(PicurvAssertBool((PetscBool)(simCtx->bboxlist != NULL), "runtime fixture should gather bboxlist through SetupDomainRankInfo"));
    PetscCall(PicurvAssertBool((PetscBool)(user->RankCellInfoMap != NULL), "runtime fixture should gather rank-cell ownership metadata"));
    PetscCall(PicurvAssertBool((PetscBool)(simCtx->bboxlist[simCtx->rank].max_coords.x >= user->bbox.max_coords.x - 1.0e-10),
                               "runtime fixture rank bbox entry should include the local bbox extent"));

    PetscCall(FreeLifecycleContext(&simCtx));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests the core setup lifecycle through environment, grid, BC, rank-info, and Eulerian-state initialization.
 */
static PetscErrorCode TestSetupLifecycleCoreSolverSetup(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char results_dir[PETSC_MAX_PATH_LEN];
    char logs_dir[PETSC_MAX_PATH_LEN];
    Cmpnts ***ucont = NULL;

    PetscFunctionBeginUser;
    PetscCall(BuildLifecycleContext(PETSC_FALSE, &simCtx, tmpdir, sizeof(tmpdir)));
    user = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;

    PetscCall(PetscSNPrintf(results_dir, sizeof(results_dir), "%s/results", tmpdir));
    PetscCall(PetscSNPrintf(logs_dir, sizeof(logs_dir), "%s/logs", tmpdir));
    PetscCall(PicurvAssertBool((PetscBool)(user->da != NULL), "SetupGridAndSolvers should allocate da"));
    PetscCall(PicurvAssertBool((PetscBool)(user->fda != NULL), "SetupGridAndSolvers should allocate coordinate DM"));
    PetscCall(PicurvAssertIntEqual(user->IM + 1, user->info.mx, "DM node count should match IM+1"));
    PetscCall(PicurvAssertBool(user->inletFaceDefined, "SetupBoundaryConditions should identify the inlet face"));
    PetscCall(PicurvAssertIntEqual(BC_FACE_NEG_Z, user->identifiedInletBCFace, "NEG_Z should be the configured inlet face"));
    PetscCall(PicurvAssertBool((PetscBool)(simCtx->bboxlist != NULL), "SetupDomainRankInfo should build bboxlist"));
    PetscCall(PicurvAssertBool((PetscBool)(user->RankCellInfoMap != NULL), "SetupDomainRankInfo should build rank-cell decomposition map"));
    PetscCall(AssertDirectoryExists(results_dir, "SetupSimulationEnvironment should create the output directory"));
    PetscCall(AssertDirectoryExists(logs_dir, "SetupSimulationEnvironment should create the log directory"));

    PetscCall(InitializeEulerianState(simCtx));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, &ucont));
    PetscCall(PicurvAssertBool((PetscBool)(ucont[1][1][1].z > 0.0),
                               "InitializeEulerianState should seed a positive inlet-aligned interior field"));
    PetscCall(PicurvAssertRealNear(ucont[1][1][1].z, ucont[2][2][2].z, 1.0e-12,
                                   "InitializeEulerianState should initialize a spatially consistent interior field"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, &ucont));

    PetscCall(FreeLifecycleContext(&simCtx));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests particle-swarm initialization and deterministic settlement on a tiny fully initialized case.
 */
static PetscErrorCode TestSetupLifecycleParticleInitialization(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    PetscInt nlocal = 0;
    PetscInt *cell_ids = NULL;
    PetscInt *status = NULL;

    PetscFunctionBeginUser;
    PetscCall(BuildLifecycleContext(PETSC_TRUE, &simCtx, tmpdir, sizeof(tmpdir)));
    PetscCall(InitializeEulerianState(simCtx));
    PetscCall(InitializeParticleSwarm(simCtx));
    PetscCall(PerformInitializedParticleSetup(simCtx));

    user = simCtx->usermg.mgctx[simCtx->usermg.mglevels - 1].user;
    PetscCall(PicurvAssertBool((PetscBool)(user->swarm != NULL), "InitializeParticleSwarm should allocate the solver swarm"));
    PetscCall(DMSwarmGetLocalSize(user->swarm, &nlocal));
    PetscCall(PicurvAssertIntEqual(8, nlocal, "single-rank lifecycle particle setup should own all seeded particles"));

    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    for (PetscInt p = 0; p < nlocal; ++p) {
        PetscCall(PicurvAssertBool((PetscBool)(cell_ids[3 * p + 0] >= 0), "settled particles should have a valid i cell id"));
        PetscCall(PicurvAssertBool((PetscBool)(cell_ids[3 * p + 1] >= 0), "settled particles should have a valid j cell id"));
        PetscCall(PicurvAssertBool((PetscBool)(cell_ids[3 * p + 2] >= 0), "settled particles should have a valid k cell id"));
        PetscCall(PicurvAssertIntEqual(ACTIVE_AND_LOCATED, status[p], "settled particles should be marked ACTIVE_AND_LOCATED"));
    }
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));

    PetscCall(FreeLifecycleContext(&simCtx));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests standalone RNG initialization helpers and minimal-context cleanup.
 */
static PetscErrorCode TestSetupLifecycleRandomGeneratorsAndCleanup(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscRandom randx = NULL, randy = NULL, randz = NULL;
    PetscRandom rand_i = NULL, rand_j = NULL, rand_k = NULL;
    PetscScalar sample_x = 0.0, sample_i = 0.0;
    PetscReal seconds = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    user->bbox.min_coords.x = 0.0;
    user->bbox.min_coords.y = 0.0;
    user->bbox.min_coords.z = 0.0;
    user->bbox.max_coords.x = 1.0;
    user->bbox.max_coords.y = 2.0;
    user->bbox.max_coords.z = 3.0;

    PetscCall(InitializeRandomGenerators(user, &randx, &randy, &randz));
    PetscCall(InitializeLogicalSpaceRNGs(&rand_i, &rand_j, &rand_k));
    PetscCall(InitializeBrownianRNG(simCtx));
    PetscCall(PicurvAssertBool((PetscBool)(simCtx->BrownianMotionRNG != NULL), "InitializeBrownianRNG should allocate the Brownian RNG"));
    PetscCall(PicurvAssertBool(RuntimeWalltimeGuardParsePositiveSeconds("12.5", &seconds), "RuntimeWalltimeGuardParsePositiveSeconds should parse positive numeric strings"));
    PetscCall(PicurvAssertRealNear(12.5, seconds, 1.0e-12, "RuntimeWalltimeGuardParsePositiveSeconds parsed value"));

    PetscCall(PetscRandomGetValue(randx, &sample_x));
    PetscCall(PetscRandomGetValue(rand_i, &sample_i));
    PetscCall(PicurvAssertBool((PetscBool)(PetscRealPart(sample_x) >= 0.0 && PetscRealPart(sample_x) <= 1.0),
                               "InitializeRandomGenerators should honor the configured bbox interval"));
    PetscCall(PicurvAssertBool((PetscBool)(PetscRealPart(sample_i) >= 0.0 && PetscRealPart(sample_i) <= 1.0),
                               "InitializeLogicalSpaceRNGs should generate logical coordinates in [0,1]"));

    PetscCall(PetscRandomDestroy(&randx));
    PetscCall(PetscRandomDestroy(&randy));
    PetscCall(PetscRandomDestroy(&randz));
    PetscCall(PetscRandomDestroy(&rand_i));
    PetscCall(PetscRandomDestroy(&rand_j));
    PetscCall(PetscRandomDestroy(&rand_k));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests cleanup after partial and fuller setup states without requiring unsupported double-finalization behavior.
 */
static PetscErrorCode TestSetupLifecycleCleanupAcrossInitializationStates(void)
{
    SimCtx *context_only = NULL;
    SimCtx *grid_only = NULL;
    char context_tmpdir[PETSC_MAX_PATH_LEN];
    char grid_tmpdir[PETSC_MAX_PATH_LEN];

    PetscFunctionBeginUser;
    PetscCall(BuildContextOnly(&context_only, context_tmpdir, sizeof(context_tmpdir)));
    PetscCall(PicurvAssertBool((PetscBool)(context_only != NULL), "CreateSimulationContext should allocate the top-level SimCtx"));
    PetscCall(PicurvAssertIntEqual(1, context_only->block_number, "CreateSimulationContext should parse the configured block count"));
    PetscCall(PetscOptionsClear(NULL));
    PetscCall(FreeLifecycleContext(&context_only));

    PetscCall(BuildContextOnly(&grid_only, grid_tmpdir, sizeof(grid_tmpdir)));
    PetscCall(SetupSimulationEnvironment(grid_only));
    PetscCall(SetupGridAndSolvers(grid_only));
    PetscCall(PicurvAssertBool((PetscBool)(grid_only->usermg.mgctx[grid_only->usermg.mglevels - 1].user->Ucont != NULL),
                               "SetupGridAndSolvers should allocate baseline Eulerian vectors"));
    PetscCall(FreeLifecycleContext(&grid_only));
    PetscFunctionReturn(0);
}

/**
 * @brief Runs the unit-setup PETSc test binary.
 */
int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"setup-lifecycle-core-solver-setup", TestSetupLifecycleCoreSolverSetup},
        {"setup-lifecycle-particle-initialization", TestSetupLifecycleParticleInitialization},
        {"setup-lifecycle-random-generators-and-cleanup", TestSetupLifecycleRandomGeneratorsAndCleanup},
        {"setup-lifecycle-cleanup-across-initialization-states", TestSetupLifecycleCleanupAcrossInitializationStates},
        {"shared-runtime-fixture-contracts", TestSharedRuntimeFixtureContracts},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv setup lifecycle tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-setup", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
