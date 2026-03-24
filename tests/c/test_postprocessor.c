/**
 * @file test_postprocessor.c
 * @brief C unit tests for postprocessor orchestration and output paths.
 */

#include "test_support.h"

#include "postprocessor.h"
#include "setup.h"

#include <stdio.h>
#include <string.h>

/**
 * @brief Assert that a file path does not exist.
 * @param[in] path File path to check.
 * @param[in] context Failure context message.
 * @return Petsc error code.
 */
static PetscErrorCode PicurvAssertFileMissing(const char *path, const char *context)
{
    PetscBool exists = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(PetscTestFile(path, 'r', &exists));
    if (exists) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "[FAIL] %s | unexpected file: %s\n", context, path));
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Unexpected file is present.");
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Count file lines that begin with a given prefix.
 * @param[in] path File path to scan.
 * @param[in] prefix Prefix to match at the beginning of each line.
 * @param[out] count_out Count of matching lines.
 * @return Petsc error code.
 */
static PetscErrorCode CountLinesWithPrefix(const char *path, const char *prefix, PetscInt *count_out)
{
    FILE *file = NULL;
    char line[4096];
    PetscInt count = 0;

    PetscFunctionBeginUser;
    PetscCheck(count_out != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "count_out must not be NULL.");
    file = fopen(path, "r");
    PetscCheck(file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to open %s.", path);
    while (fgets(line, sizeof(line), file)) {
        if (strncmp(line, prefix, strlen(prefix)) == 0) count++;
    }
    fclose(file);
    *count_out = count;
    PetscFunctionReturn(0);
}
/**
 * @brief Tests post-processing swarm setup and pipeline-field registration.
 */

static PetscErrorCode TestSetupPostProcessSwarmRegistersPipelineFields(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PostProcessParams pps;
    PetscInt bs = 0;
    void *field_ptr = NULL;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&pps, sizeof(pps)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PetscStrncpy(pps.particle_pipeline, "ComputeSpecificKE:velocity>ske", sizeof(pps.particle_pipeline)));

    PetscCall(SetupPostProcessSwarm(user, &pps));
    PetscCall(PicurvAssertBool((PetscBool)(user->post_swarm != NULL),
                               "SetupPostProcessSwarm should allocate a dedicated post swarm"));

    PetscCall(DMSwarmSetLocalSizes(user->post_swarm, 0, 0));
    PetscCall(DMSwarmGetField(user->post_swarm, "ske", &bs, NULL, &field_ptr));
    PetscCall(PicurvAssertIntEqual(1, bs, "SetupPostProcessSwarm should register scalar output field 'ske'"));
    PetscCall(DMSwarmRestoreField(user->post_swarm, "ske", &bs, NULL, &field_ptr));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests Eulerian post-processing pipeline dispatch for configured kernels.
 */

static PetscErrorCode TestEulerianDataProcessingPipelineRunsConfiguredKernels(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PostProcessParams pps;
    const PetscScalar ***p_nodal_arr = NULL;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&pps, sizeof(pps)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PetscStrncpy(pps.process_pipeline, "UnknownStep;CellToNodeAverage:P>P_nodal", sizeof(pps.process_pipeline)));

    PetscCall(VecSet(user->P, 9.0));
    PetscCall(VecSet(user->P_nodal, -1.0));
    PetscCall(EulerianDataProcessingPipeline(user, &pps));

    PetscCall(DMDAVecGetArrayRead(user->da, user->P_nodal, (void *)&p_nodal_arr));
    PetscCall(PicurvAssertRealNear(9.0, PetscRealPart(p_nodal_arr[0][0][0]), 1.0e-12,
                                   "Eulerian pipeline should execute CellToNodeAverage for P->P_nodal"));
    PetscCall(PicurvAssertRealNear(-1.0, PetscRealPart(p_nodal_arr[user->KM][user->JM][user->IM]), 1.0e-12,
                                   "Eulerian pipeline should leave the extra non-physical boundary node unchanged"));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->P_nodal, (void *)&p_nodal_arr));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests particle post-processing pipeline dispatch for specific-KE output.
 */

static PetscErrorCode TestParticleDataProcessingPipelineComputesSpecificKE(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PostProcessParams pps;
    PetscReal (*vel_arr)[3] = NULL;
    PetscScalar *ske_arr = NULL;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&pps, sizeof(pps)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 2, "ske"));
    PetscCall(PetscStrncpy(pps.particle_pipeline, "ComputeSpecificKE:velocity>ske", sizeof(pps.particle_pipeline)));

    PetscCall(DMSwarmGetField(user->swarm, "velocity", NULL, NULL, (void *)&vel_arr));
    vel_arr[0][0] = 1.0;
    vel_arr[0][1] = 2.0;
    vel_arr[0][2] = 2.0;
    vel_arr[1][0] = 0.0;
    vel_arr[1][1] = 3.0;
    vel_arr[1][2] = 4.0;
    PetscCall(DMSwarmRestoreField(user->swarm, "velocity", NULL, NULL, (void *)&vel_arr));

    PetscCall(ParticleDataProcessingPipeline(user, &pps));
    PetscCall(DMSwarmGetField(user->post_swarm, "ske", NULL, NULL, (void *)&ske_arr));
    PetscCall(PicurvAssertRealNear(4.5, PetscRealPart(ske_arr[0]), 1.0e-12,
                                   "Particle pipeline should write first specific kinetic energy value"));
    PetscCall(PicurvAssertRealNear(12.5, PetscRealPart(ske_arr[1]), 1.0e-12,
                                   "Particle pipeline should write second specific kinetic energy value"));
    PetscCall(DMSwarmRestoreField(user->post_swarm, "ske", NULL, NULL, (void *)&ske_arr));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests global statistics pipeline generation of MSD CSV output.
 */

static PetscErrorCode TestGlobalStatisticsPipelineWritesMSDCSV(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PostProcessParams pps;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char csv_prefix[PETSC_MAX_PATH_LEN];
    char csv_path[PETSC_MAX_PATH_LEN];
    PetscReal (*pos_arr)[3] = NULL;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&pps, sizeof(pps)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 2, "ske"));

    simCtx->dt = 1.0;
    simCtx->ren = 1.0;
    simCtx->schmidt_number = 1.0;
    simCtx->psrc_x = 0.0;
    simCtx->psrc_y = 0.0;
    simCtx->psrc_z = 0.0;

    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(csv_prefix, sizeof(csv_prefix), "%s/stats", tmpdir));
    PetscCall(PetscSNPrintf(csv_path, sizeof(csv_path), "%s_msd.csv", csv_prefix));

    PetscCall(PetscStrncpy(pps.statistics_pipeline, "ComputeMSD", sizeof(pps.statistics_pipeline)));
    PetscCall(PetscStrncpy(pps.statistics_output_prefix, csv_prefix, sizeof(pps.statistics_output_prefix)));

    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));
    pos_arr[0][0] = 1.0;
    pos_arr[0][1] = 0.0;
    pos_arr[0][2] = 0.0;
    pos_arr[1][0] = -1.0;
    pos_arr[1][1] = 0.0;
    pos_arr[1][2] = 0.0;
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));

    PetscCall(GlobalStatisticsPipeline(user, &pps, 1));
    PetscCall(PicurvAssertFileExists(csv_path, "GlobalStatisticsPipeline should dispatch ComputeMSD and emit CSV output"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that reprocessing the same MSD step rewrites the CSV row instead of duplicating it.
 */
static PetscErrorCode TestGlobalStatisticsPipelineRewritesExistingMSDStep(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PostProcessParams pps;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char csv_prefix[PETSC_MAX_PATH_LEN];
    char csv_path[PETSC_MAX_PATH_LEN];
    char csv_tmp_path[PETSC_MAX_PATH_LEN];
    PetscReal (*pos_arr)[3] = NULL;
    PetscInt row_count = 0;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&pps, sizeof(pps)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 2, "ske"));

    simCtx->dt = 1.0;
    simCtx->ren = 1.0;
    simCtx->schmidt_number = 1.0;
    simCtx->psrc_x = 0.0;
    simCtx->psrc_y = 0.0;
    simCtx->psrc_z = 0.0;

    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(csv_prefix, sizeof(csv_prefix), "%s/stats", tmpdir));
    PetscCall(PetscSNPrintf(csv_path, sizeof(csv_path), "%s_msd.csv", csv_prefix));
    PetscCall(PetscSNPrintf(csv_tmp_path, sizeof(csv_tmp_path), "%s.tmp", csv_path));

    PetscCall(PetscStrncpy(pps.statistics_pipeline, "ComputeMSD", sizeof(pps.statistics_pipeline)));
    PetscCall(PetscStrncpy(pps.statistics_output_prefix, csv_prefix, sizeof(pps.statistics_output_prefix)));

    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));
    pos_arr[0][0] = 1.0;
    pos_arr[0][1] = 0.0;
    pos_arr[0][2] = 0.0;
    pos_arr[1][0] = -1.0;
    pos_arr[1][1] = 0.0;
    pos_arr[1][2] = 0.0;
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));
    PetscCall(GlobalStatisticsPipeline(user, &pps, 1));

    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));
    pos_arr[0][0] = 2.0;
    pos_arr[0][1] = 0.0;
    pos_arr[0][2] = 0.0;
    pos_arr[1][0] = -2.0;
    pos_arr[1][1] = 0.0;
    pos_arr[1][2] = 0.0;
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));
    PetscCall(GlobalStatisticsPipeline(user, &pps, 1));

    PetscCall(CountLinesWithPrefix(csv_path, "1,", &row_count));
    PetscCall(PicurvAssertIntEqual(1, row_count,
                                   "Reprocessing the same MSD step should keep one CSV row for that step"));
    PetscCall(PicurvAssertFileMissing(csv_tmp_path,
                                      "Reprocessing the same MSD step should not leave a temporary CSV artifact behind"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that postprocessor setup creates nested statistics directories without clobbering existing contents.
 */
static PetscErrorCode TestSetupSimulationEnvironmentCreatesStatisticsDirectoryWithoutClobbering(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PostProcessParams *pps = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char post_path[PETSC_MAX_PATH_LEN];
    char source_dir[PETSC_MAX_PATH_LEN];
    char stats_dir[PETSC_MAX_PATH_LEN];
    char stats_prefix[PETSC_MAX_PATH_LEN];
    char sentinel_path[PETSC_MAX_PATH_LEN];
    PetscBool exists = PETSC_FALSE;
    FILE *file = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PetscNew(&simCtx->pps));
    pps = simCtx->pps;
    PetscCall(PetscMemzero(pps, sizeof(*pps)));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));

    PetscCall(PetscSNPrintf(post_path, sizeof(post_path), "%s/post.run", tmpdir));
    file = fopen(post_path, "w");
    PetscCheck(file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to create post config %s.", post_path);
    fputs("startTime = 0\nendTime = 0\ntimeStep = 1\n", file);
    fclose(file);
    file = NULL;

    PetscCall(PetscSNPrintf(source_dir, sizeof(source_dir), "%s/output", tmpdir));
    PetscCall(PicurvEnsureDir(source_dir));
    PetscCall(PetscSNPrintf(stats_dir, sizeof(stats_dir), "%s/output/statistics", tmpdir));
    PetscCall(PetscSNPrintf(stats_prefix, sizeof(stats_prefix), "%s/BrownianStats", stats_dir));
    PetscCall(PetscSNPrintf(sentinel_path, sizeof(sentinel_path), "%s/sentinel.keep", stats_dir));

    simCtx->exec_mode = EXEC_MODE_POSTPROCESSOR;
    simCtx->generate_grid = PETSC_TRUE;
    PetscCall(PetscStrncpy(simCtx->PostprocessingControlFile, post_path, sizeof(simCtx->PostprocessingControlFile)));
    PetscCall(PetscStrncpy(pps->source_dir, source_dir, sizeof(pps->source_dir)));
    PetscCall(PetscSNPrintf(pps->output_prefix, sizeof(pps->output_prefix), "%s/viz/Field", tmpdir));
    PetscCall(PetscStrncpy(pps->statistics_pipeline, "ComputeMSD", sizeof(pps->statistics_pipeline)));
    PetscCall(PetscStrncpy(pps->statistics_output_prefix, stats_prefix, sizeof(pps->statistics_output_prefix)));

    PetscCall(SetupSimulationEnvironment(simCtx));
    PetscCall(PetscTestDirectory(stats_dir, 'r', &exists));
    PetscCall(PicurvAssertBool(exists,
                               "SetupSimulationEnvironment should create the statistics directory when it is missing"));

    file = fopen(sentinel_path, "w");
    PetscCheck(file != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Failed to create statistics sentinel %s.", sentinel_path);
    fputs("keep\n", file);
    fclose(file);
    file = NULL;

    PetscCall(SetupSimulationEnvironment(simCtx));
    PetscCall(PicurvAssertFileExists(sentinel_path,
                                     "SetupSimulationEnvironment should preserve existing statistics directory contents on repeat runs"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests VTS emission for Eulerian post-processing output.
 */

static PetscErrorCode TestWriteEulerianFileWritesVTS(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PostProcessParams pps;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char vtk_path[PETSC_MAX_PATH_LEN];

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&pps, sizeof(pps)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(vtk_path, sizeof(vtk_path), "%s/field_00003.vts", tmpdir));

    PetscCall(PetscStrncpy(pps.output_fields_instantaneous, "P_nodal", sizeof(pps.output_fields_instantaneous)));
    PetscCall(PetscSNPrintf(pps.output_prefix, sizeof(pps.output_prefix), "%s/field", tmpdir));
    PetscCall(VecSet(user->P_nodal, 7.0));

    PetscCall(WriteEulerianFile(user, &pps, 3));
    PetscCall(PicurvAssertFileExists(vtk_path, "WriteEulerianFile should emit a .vts file for requested output fields"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests VTP emission for particle post-processing output.
 */

static PetscErrorCode TestWriteParticleFileWritesVTP(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PostProcessParams pps;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char vtk_path[PETSC_MAX_PATH_LEN];
    PetscReal (*pos_arr)[3] = NULL;
    PetscReal (*vel_arr)[3] = NULL;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&pps, sizeof(pps)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 3, "ske"));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(vtk_path, sizeof(vtk_path), "%s/particles_00004.vtp", tmpdir));

    pps.outputParticles = PETSC_TRUE;
    pps.particle_output_freq = 1;
    PetscCall(PetscStrncpy(pps.particle_fields, "position,velocity", sizeof(pps.particle_fields)));
    PetscCall(PetscSNPrintf(pps.particle_output_prefix, sizeof(pps.particle_output_prefix), "%s/particles", tmpdir));

    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));
    PetscCall(DMSwarmGetField(user->swarm, "velocity", NULL, NULL, (void *)&vel_arr));
    for (PetscInt p = 0; p < 3; ++p) {
        pos_arr[p][0] = (PetscReal)p;
        pos_arr[p][1] = (PetscReal)(p + 1);
        pos_arr[p][2] = (PetscReal)(p + 2);
        vel_arr[p][0] = 2.0 * (PetscReal)p;
        vel_arr[p][1] = 3.0 * (PetscReal)p;
        vel_arr[p][2] = 4.0 * (PetscReal)p;
    }
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));
    PetscCall(DMSwarmRestoreField(user->swarm, "velocity", NULL, NULL, (void *)&vel_arr));

    PetscCall(WriteParticleFile(user, &pps, 4));
    PetscCall(PicurvAssertFileExists(vtk_path, "WriteParticleFile should emit a .vtp file for requested particle fields"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests repeated same-step VTS emission leaves only the final artifact.
 */
static PetscErrorCode TestWriteEulerianFileRewritesSameStepCleanly(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PostProcessParams pps;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char vtk_path[PETSC_MAX_PATH_LEN];
    char vtk_tmp_path[PETSC_MAX_PATH_LEN];

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&pps, sizeof(pps)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(vtk_path, sizeof(vtk_path), "%s/field_00003.vts", tmpdir));
    PetscCall(PetscSNPrintf(vtk_tmp_path, sizeof(vtk_tmp_path), "%s.tmp", vtk_path));

    PetscCall(PetscStrncpy(pps.output_fields_instantaneous, "P_nodal", sizeof(pps.output_fields_instantaneous)));
    PetscCall(PetscSNPrintf(pps.output_prefix, sizeof(pps.output_prefix), "%s/field", tmpdir));
    PetscCall(VecSet(user->P_nodal, 7.0));

    PetscCall(WriteEulerianFile(user, &pps, 3));
    PetscCall(VecSet(user->P_nodal, 9.0));
    PetscCall(WriteEulerianFile(user, &pps, 3));

    PetscCall(PicurvAssertFileExists(vtk_path,
                                     "Rewriting the same Eulerian step should leave the final .vts file present"));
    PetscCall(PicurvAssertFileMissing(vtk_tmp_path,
                                      "Rewriting the same Eulerian step should not leave a temporary .vts file behind"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests repeated same-step VTP emission leaves only the final artifact.
 */

/**
 * @brief Tests repeated same-step VTP emission leaves only the final artifact.
 */
static PetscErrorCode TestWriteParticleFileRewritesSameStepCleanly(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PostProcessParams pps;
    char tmpdir[PETSC_MAX_PATH_LEN];
    char vtk_path[PETSC_MAX_PATH_LEN];
    char vtk_tmp_path[PETSC_MAX_PATH_LEN];
    PetscReal (*pos_arr)[3] = NULL;
    PetscReal (*vel_arr)[3] = NULL;

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&pps, sizeof(pps)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 3, "ske"));
    PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCall(PetscSNPrintf(vtk_path, sizeof(vtk_path), "%s/particles_00004.vtp", tmpdir));
    PetscCall(PetscSNPrintf(vtk_tmp_path, sizeof(vtk_tmp_path), "%s.tmp", vtk_path));

    pps.outputParticles = PETSC_TRUE;
    pps.particle_output_freq = 1;
    PetscCall(PetscStrncpy(pps.particle_fields, "position,velocity", sizeof(pps.particle_fields)));
    PetscCall(PetscSNPrintf(pps.particle_output_prefix, sizeof(pps.particle_output_prefix), "%s/particles", tmpdir));

    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));
    PetscCall(DMSwarmGetField(user->swarm, "velocity", NULL, NULL, (void *)&vel_arr));
    for (PetscInt p = 0; p < 3; ++p) {
        pos_arr[p][0] = (PetscReal)p;
        pos_arr[p][1] = (PetscReal)(p + 1);
        pos_arr[p][2] = (PetscReal)(p + 2);
        vel_arr[p][0] = 2.0 * (PetscReal)p;
        vel_arr[p][1] = 3.0 * (PetscReal)p;
        vel_arr[p][2] = 4.0 * (PetscReal)p;
    }
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void *)&pos_arr));
    PetscCall(DMSwarmRestoreField(user->swarm, "velocity", NULL, NULL, (void *)&vel_arr));

    PetscCall(WriteParticleFile(user, &pps, 4));
    PetscCall(WriteParticleFile(user, &pps, 4));

    PetscCall(PicurvAssertFileExists(vtk_path,
                                     "Rewriting the same particle step should leave the final .vtp file present"));
    PetscCall(PicurvAssertFileMissing(vtk_tmp_path,
                                      "Rewriting the same particle step should not leave a temporary .vtp file behind"));

    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Runs the unit-postprocessor PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"setup-postprocess-swarm-registers-pipeline-fields", TestSetupPostProcessSwarmRegistersPipelineFields},
        {"eulerian-data-processing-pipeline-runs-configured-kernels", TestEulerianDataProcessingPipelineRunsConfiguredKernels},
        {"particle-data-processing-pipeline-computes-specific-ke", TestParticleDataProcessingPipelineComputesSpecificKE},
        {"global-statistics-pipeline-writes-msd-csv", TestGlobalStatisticsPipelineWritesMSDCSV},
        {"global-statistics-pipeline-rewrites-existing-msd-step", TestGlobalStatisticsPipelineRewritesExistingMSDStep},
        {"setup-simulation-environment-creates-statistics-directory-without-clobbering", TestSetupSimulationEnvironmentCreatesStatisticsDirectoryWithoutClobbering},
        {"write-eulerian-file-writes-vts", TestWriteEulerianFileWritesVTS},
        {"write-particle-file-writes-vtp", TestWriteParticleFileWritesVTP},
        {"write-eulerian-file-rewrites-same-step-cleanly", TestWriteEulerianFileRewritesSameStepCleanly},
        {"write-particle-file-rewrites-same-step-cleanly", TestWriteParticleFileRewritesSameStepCleanly},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv postprocessor tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-postprocessor", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
