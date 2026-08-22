/**
 * @file test_post_compute_mpi.c
 * @brief Decomposition-independent regression tests for post-processing compute.
 */

#include "test_support.h"

#include "postprocessor.h"
#include "statistics_accumulator.h"
#include "statistics_window.h"

static const PetscReal kSentinel = -9876.5;

/** @brief Analytic scalar field used to expose misplaced or missing grid points. */
static PetscReal ScalarValue(PetscInt i, PetscInt j, PetscInt k)
{
    return (PetscReal)(i + 10 * j + 100 * k);
}

/** @brief Requires every rank's maximum error to be within tolerance. */
static PetscErrorCode AssertGlobalError(PetscReal local_error, PetscReal tolerance,
                                        const char *context)
{
    PetscReal global_error = 0.0;

    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Allreduce(&local_error, &global_error, 1, MPIU_REAL, MPI_MAX,
                               PETSC_COMM_WORLD));
    PetscCall(PicurvAssertBool((PetscBool)(global_error <= tolerance), context));
    PetscFunctionReturn(0);
}

/** @brief Exercises every Eulerian transformation currently dispatched by the pipeline. */
static PetscErrorCode TestEulerianPipelineDecompositionIndependent(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PostProcessParams *pps = NULL;
    PetscReal ***pressure = NULL;
    Cmpnts ***velocity = NULL;
    const PetscReal ***pressure_result = NULL;
    const PetscReal ***pressure_nodal = NULL;
    const PetscReal ***qcrit = NULL;
    const Cmpnts ***velocity_nodal = NULL;
    PetscReal pressure_error = 0.0;
    PetscReal pressure_nodal_error = 0.0;
    PetscReal velocity_nodal_error = 0.0;
    PetscReal qcrit_error = 0.0;
    PetscInt local_physical_nodes = 0;
    PetscInt global_physical_nodes = 0;
    PetscInt local_q_points = 0;
    PetscInt global_q_points = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 8, 6, 4));
    PetscCall(PetscCalloc1(1, &simCtx->pps));
    pps = simCtx->pps;
    pps->reference[0] = 1;
    pps->reference[1] = 1;
    pps->reference[2] = 1;
    simCtx->scaling.L_ref = 1.0;
    simCtx->scaling.U_ref = 1.0;
    simCtx->scaling.P_ref = 2.0;

    PetscCall(DMDAVecGetArray(user->da, user->P, &pressure));
    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &velocity));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
                pressure[k][j][i] = ScalarValue(i, j, k);
                velocity[k][j][i] = (Cmpnts){1.0, 2.0, 3.0};
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &velocity));
    PetscCall(DMDAVecRestoreArray(user->da, user->P, &pressure));
    PetscCall(VecSet(user->P_nodal, kSentinel));
    PetscCall(VecSet(user->Ucat_nodal, kSentinel));
    PetscCall(VecSet(user->Qcrit, kSentinel));

    PetscCall(PetscStrncpy(
        pps->process_pipeline,
        "NormalizeRelativeField:P;DimensionalizeAllLoadedFields;"
        "CellToNodeAverage:P>P_nodal;CellToNodeAverage:Ucat>Ucat_nodal;"
        "ComputeQCriterion",
        sizeof(pps->process_pipeline)));
    PetscCall(EulerianDataProcessingPipeline(user, pps));

    PetscCall(DMDAVecGetArrayRead(user->da, user->P, &pressure_result));
    PetscCall(DMDAVecGetArrayRead(user->da, user->P_nodal, &pressure_nodal));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucat_nodal, &velocity_nodal));
    PetscCall(DMDAVecGetArrayRead(user->da, user->Qcrit, &qcrit));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j) {
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i) {
                const PetscBool physical_node =
                    (PetscBool)(i < user->info.mx - 1 && j < user->info.my - 1 &&
                                k < user->info.mz - 1);
                const PetscBool q_point =
                    (PetscBool)(i >= 1 && i < user->info.mx - 1 &&
                                j >= 1 && j < user->info.my - 1 &&
                                k >= 1 && k < user->info.mz - 1);
                const PetscReal expected_pressure = 2.0 * (ScalarValue(i, j, k) - 111.0);

                pressure_error = PetscMax(
                    pressure_error,
                    PetscAbsReal(pressure_result[k][j][i] - expected_pressure));
                if (physical_node) {
                    const PetscReal expected_nodal =
                        2.0 * (ScalarValue(i, j, k) + 55.5 - 111.0);

                    pressure_nodal_error = PetscMax(
                        pressure_nodal_error,
                        PetscAbsReal(pressure_nodal[k][j][i] - expected_nodal));
                    velocity_nodal_error = PetscMax(
                        velocity_nodal_error,
                        PetscAbsReal(velocity_nodal[k][j][i].x - 1.0));
                    velocity_nodal_error = PetscMax(
                        velocity_nodal_error,
                        PetscAbsReal(velocity_nodal[k][j][i].y - 2.0));
                    velocity_nodal_error = PetscMax(
                        velocity_nodal_error,
                        PetscAbsReal(velocity_nodal[k][j][i].z - 3.0));
                    ++local_physical_nodes;
                } else {
                    pressure_nodal_error = PetscMax(
                        pressure_nodal_error,
                        PetscAbsReal(pressure_nodal[k][j][i] - kSentinel));
                    velocity_nodal_error = PetscMax(
                        velocity_nodal_error,
                        PetscAbsReal(velocity_nodal[k][j][i].x - kSentinel));
                    velocity_nodal_error = PetscMax(
                        velocity_nodal_error,
                        PetscAbsReal(velocity_nodal[k][j][i].y - kSentinel));
                    velocity_nodal_error = PetscMax(
                        velocity_nodal_error,
                        PetscAbsReal(velocity_nodal[k][j][i].z - kSentinel));
                }
                qcrit_error = PetscMax(
                    qcrit_error,
                    PetscAbsReal(qcrit[k][j][i] - (q_point ? 0.0 : kSentinel)));
                if (q_point) ++local_q_points;
            }
        }
    }
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->Qcrit, &qcrit));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucat_nodal, &velocity_nodal));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->P_nodal, &pressure_nodal));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->P, &pressure_result));

    PetscCallMPI(MPI_Allreduce(&local_physical_nodes, &global_physical_nodes, 1,
                               MPIU_INT, MPI_SUM, PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Allreduce(&local_q_points, &global_q_points, 1,
                               MPIU_INT, MPI_SUM, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertIntEqual(user->IM * user->JM * user->KM,
                                   global_physical_nodes,
                                   "the Eulerian pipeline must visit every physical node exactly once"));
    PetscCall(PicurvAssertIntEqual((user->IM - 1) * (user->JM - 1) * (user->KM - 1),
                                   global_q_points,
                                   "Q-criterion must visit every interior cell exactly once"));
    PetscCall(AssertGlobalError(pressure_error, 1.0e-12,
                                "pressure normalization and dimensionalization are decomposition independent"));
    PetscCall(AssertGlobalError(pressure_nodal_error, 1.0e-12,
                                "scalar nodal averaging covers rank interfaces and physical boundaries"));
    PetscCall(AssertGlobalError(velocity_nodal_error, 1.0e-12,
                                "vector nodal averaging covers rank interfaces and physical boundaries"));
    PetscCall(AssertGlobalError(qcrit_error, 1.0e-12,
                                "Q-criterion covers the complete distributed interior"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief Verifies the production-created derived swarm mirrors source ownership. */
static PetscErrorCode TestParticlePipelineMatchesSourceOwnership(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PostProcessParams pps;
    PetscMPIInt rank = 0, size = 1;
    PetscInt source_local = 0;
    PetscInt source_global = 0;
    PetscInt post_local = 0;
    PetscInt post_global = 0;
    PetscInt offset = 0;
    PetscScalar (*velocity)[3] = NULL;
    const PetscScalar *ske = NULL;
    PetscReal local_error = 0.0;
    const PetscInt total_particles = 7;

    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    source_local = total_particles / size + (rank < total_particles % size ? 1 : 0);
    PetscCall(PetscMemzero(&pps, sizeof(pps)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, source_local, "unused"));
    PetscCall(DMDestroy(&user->post_swarm));
    PetscCall(PetscStrncpy(pps.particle_pipeline, "ComputeSpecificKE:velocity>ske",
                           sizeof(pps.particle_pipeline)));
    PetscCall(SetupPostProcessSwarm(user, &pps));

    if (size > 1) {
        PetscCallMPI(MPI_Exscan(&source_local, &offset, 1, MPIU_INT, MPI_SUM,
                                PETSC_COMM_WORLD));
        if (rank == 0) offset = 0;
    }
    PetscCall(DMSwarmGetField(user->swarm, "velocity", NULL, NULL, (void **)&velocity));
    for (PetscInt p = 0; p < source_local; ++p) {
        const PetscReal tag = (PetscReal)(offset + p + 1);

        velocity[p][0] = tag;
        velocity[p][1] = 2.0;
        velocity[p][2] = -1.0;
    }
    PetscCall(DMSwarmRestoreField(user->swarm, "velocity", NULL, NULL,
                                  (void **)&velocity));

    PetscCall(ParticleDataProcessingPipeline(user, &pps));
    PetscCall(DMSwarmGetLocalSize(user->swarm, &source_local));
    PetscCall(DMSwarmGetSize(user->swarm, &source_global));
    PetscCall(DMSwarmGetLocalSize(user->post_swarm, &post_local));
    PetscCall(DMSwarmGetSize(user->post_swarm, &post_global));
    PetscCall(PicurvAssertIntEqual(source_local, post_local,
                                   "derived particle fields must share source local ownership"));
    PetscCall(PicurvAssertIntEqual(source_global, post_global,
                                   "derived and source swarms must have the same global size"));

    PetscCall(DMSwarmGetField(user->post_swarm, "ske", NULL, NULL, (void **)&ske));
    for (PetscInt p = 0; p < source_local; ++p) {
        const PetscReal tag = (PetscReal)(offset + p + 1);
        const PetscReal expected = 0.5 * (tag * tag + 5.0);

        local_error = PetscMax(local_error,
                               PetscAbsReal(PetscRealPart(ske[p]) - expected));
    }
    PetscCall(DMSwarmRestoreField(user->post_swarm, "ske", NULL, NULL, (void **)&ske));
    PetscCall(AssertGlobalError(local_error, 1.0e-12,
                                "specific kinetic energy must match every source particle"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief Checks that global resizing produces one balanced global population. */
static PetscErrorCode TestResizeSwarmGloballyBalanced(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscMPIInt rank = 0, size = 1;
    const PetscInt targets[] = {7, 3, 9, 0};

    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 0, "unused"));

    for (size_t target_index = 0;
         target_index < sizeof(targets) / sizeof(targets[0]); ++target_index) {
        PetscInt local = 0;
        PetscInt global = 0;
        const PetscInt target = targets[target_index];
        const PetscInt expected_local =
            target / size + (rank < target % size ? 1 : 0);

        PetscCall(ResizeSwarmGlobally(user->swarm, target));
        PetscCall(DMSwarmGetLocalSize(user->swarm, &local));
        PetscCall(DMSwarmGetSize(user->swarm, &global));
        PetscCall(PicurvAssertIntEqual(expected_local, local,
                                       "global swarm resize must use balanced local sizes"));
        PetscCall(PicurvAssertIntEqual(target, global,
                                       "global swarm resize must conserve the requested total"));
    }

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief Exercises accumulator-to-derived-nodal field statistics across ownership boundaries. */
static PetscErrorCode TestFieldStatisticsDerivedNodalMultiRank(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindowDefinition definition;
    PicurvWindow window;
    PicurvWindowStorage storage;
    PetscReal ***pressure = NULL;
    const PetscReal ***nodal = NULL;
    Vec result = NULL;
    PetscInt components = 0;
    PetscInt local_checked = 0;
    PetscInt global_checked = 0;
    PetscReal local_error = 0.0;
    char name[PICURV_STATISTICS_PAYLOAD_NAME_LENGTH];

    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(&definition, sizeof(definition)));
    PetscCall(PetscStrncpy(definition.name, "analytic", sizeof(definition.name)));
    definition.weighting = PICURV_WEIGHTING_SAMPLE;
    definition.cadence_kind = PICURV_CADENCE_STEP;
    definition.step_cadence = 1;
    definition.field_count = 1;
    definition.fields[0].field_id = FIELD_ID_P;

    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 8, 6, 4));
    PetscCall(VecSet(user->Nvert, 0.0));
    PetscCall(PicurvWindowInit(&window, &definition));
    PetscCall(PicurvWindowStorageCreate(user, &definition, &storage));
    simCtx->fieldStatisticsEnabled = PETSC_TRUE;
    simCtx->fieldStatisticsWindowCount = 1;
    simCtx->fieldStatisticsWindows = &window;
    user->fieldStatisticsStorage = &storage;

    PetscCall(DMDAVecGetArray(user->da, user->P, &pressure));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; ++k)
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; ++j)
            for (PetscInt i = user->info.xs; i < user->info.xs + user->info.xm; ++i)
                pressure[k][j][i] = ScalarValue(i, j, k);
    PetscCall(DMDAVecRestoreArray(user->da, user->P, &pressure));
    PetscCall(PicurvWindowAccumulate(user, &definition, &storage, 1.0));
    PetscCall(ComputeWindowStatisticNodal(user, 0, "mean", 0, name, sizeof(name),
                                          &result, &components));
    PetscCall(PicurvAssertIntEqual(1, components,
                                   "the pressure mean must derive as a scalar field"));
    PetscCall(PicurvAssertBool((PetscBool)(result == user->PostScalarNodal),
                               "the pressure mean must use the scalar nodal staging vector"));

    PetscCall(DMDAVecGetArrayRead(user->da, result, &nodal));
    for (PetscInt k = PetscMax(user->info.zs, 1);
         k < PetscMin(user->info.zs + user->info.zm, user->info.mz - 2); ++k) {
        for (PetscInt j = PetscMax(user->info.ys, 1);
             j < PetscMin(user->info.ys + user->info.ym, user->info.my - 2); ++j) {
            for (PetscInt i = PetscMax(user->info.xs, 1);
                 i < PetscMin(user->info.xs + user->info.xm, user->info.mx - 2); ++i) {
                local_error = PetscMax(
                    local_error,
                    PetscAbsReal(nodal[k][j][i] - (ScalarValue(i, j, k) + 55.5)));
                ++local_checked;
            }
        }
    }
    PetscCall(DMDAVecRestoreArrayRead(user->da, result, &nodal));
    PetscCallMPI(MPI_Allreduce(&local_checked, &global_checked, 1, MPIU_INT, MPI_SUM,
                               PETSC_COMM_WORLD));
    PetscCall(PicurvAssertIntEqual((user->IM - 2) * (user->JM - 2) * (user->KM - 2),
                                   global_checked,
                                   "the derived mean must cover every fully resolved nodal point"));
    PetscCall(AssertGlobalError(local_error, 1.0e-12,
                                "field-statistics nodal derivation must cross rank interfaces"));

    simCtx->fieldStatisticsEnabled = PETSC_FALSE;
    simCtx->fieldStatisticsWindowCount = 0;
    simCtx->fieldStatisticsWindows = NULL;
    user->fieldStatisticsStorage = NULL;
    PetscCall(PicurvWindowStorageDestroy(&storage));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief Runs the serial/MPI post-processing compute regression suite. */
int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"eulerian-pipeline-decomposition-independent", TestEulerianPipelineDecompositionIndependent},
        {"particle-pipeline-matches-source-ownership", TestParticlePipelineMatchesSourceOwnership},
        {"resize-swarm-globally-balanced", TestResizeSwarmGloballyBalanced},
        {"field-statistics-derived-nodal-multi-rank", TestFieldStatisticsDerivedNodalMultiRank},
    };

    ierr = PetscInitialize(&argc, &argv, NULL,
                           "PICurv serial/MPI post-processing compute tests");
    if (ierr) return (int)ierr;
    ierr = PicurvRunTests("unit-post-compute-mpi", cases,
                          sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }
    ierr = PetscFinalize();
    return (int)ierr;
}
