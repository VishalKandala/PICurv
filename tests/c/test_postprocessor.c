/**
 * @file test_postprocessor.c
 * @brief C unit tests for postprocessor orchestration and output paths.
 */

#include "test_support.h"

#include "postprocessor.h"

#include <stdio.h>
#include <string.h>
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
        {"write-eulerian-file-writes-vts", TestWriteEulerianFileWritesVTS},
        {"write-particle-file-writes-vtp", TestWriteParticleFileWritesVTP},
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
