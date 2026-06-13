/**
 * @file test_mpi_kernels.c
 * @brief C unit tests for MPI collective helper routines.
 */

#include "test_support.h"

#include "ParticleMotion.h"
#include "ParticleSwarm.h"
#include "Boundaries.h"
#include "grid.h"

#include <stdlib.h>

/**
 * @brief Marks both x boundary faces as geometric periodic boundaries.
 */
static void MarkXPeriodic(UserCtx *user)
{
    user->boundary_faces[BC_FACE_NEG_X].face_id = BC_FACE_NEG_X;
    user->boundary_faces[BC_FACE_NEG_X].mathematical_type = PERIODIC;
    user->boundary_faces[BC_FACE_NEG_X].handler_type = BC_HANDLER_PERIODIC_GEOMETRIC;
    user->boundary_faces[BC_FACE_POS_X].face_id = BC_FACE_POS_X;
    user->boundary_faces[BC_FACE_POS_X].mathematical_type = PERIODIC;
    user->boundary_faces[BC_FACE_POS_X].handler_type = BC_HANDLER_PERIODIC_GEOMETRIC;
}

/**
 * @brief Returns a deterministic scalar value for periodic-transfer assertions.
 */
static PetscReal PeriodicScalarValue(PetscInt i, PetscInt j, PetscInt k, PetscReal offset)
{
    return offset + (PetscReal)(i + 10 * j + 100 * k);
}
/**
 * @brief Tests collective particle distribution consistency across MPI ranks.
 */

static PetscErrorCode TestDistributeParticlesCollectiveConsistency(void)
{
    PetscMPIInt rank = 0, size = 1;
    PetscInt local_particles = 0;
    PetscInt remainder = 0;
    PetscInt global_particles = 0;
    PetscInt remainder_min = 0, remainder_max = 0;
    const PetscInt total_particles = 137;

    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    PetscCall(PicurvAssertBool((PetscBool)(size >= 2), "unit-mpi requires at least two MPI ranks"));

    PetscCall(DistributeParticles(total_particles, rank, size, &local_particles, &remainder));

    PetscCall(PicurvAssertIntEqual(
        total_particles / size + (((PetscInt)rank < remainder) ? 1 : 0),
        local_particles,
        "local particle share should match quotient+remainder policy"));

    PetscCallMPI(MPI_Allreduce(&local_particles, &global_particles, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertIntEqual(total_particles, global_particles, "distributed particle count must conserve total particles"));

    PetscCallMPI(MPI_Allreduce(&remainder, &remainder_min, 1, MPIU_INT, MPI_MIN, PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Allreduce(&remainder, &remainder_max, 1, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertIntEqual(remainder_min, remainder_max, "all ranks should report the same remainder"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests multi-rank bounding-box gather and broadcast helpers.
 */

static PetscErrorCode TestBoundingBoxCollectivesMultiRank(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    BoundingBox local_bbox;
    BoundingBox *boxes = NULL;
    PetscMPIInt rank = 0, size = 1;
    PetscReal global_min_x = 0.0;
    PetscReal global_max_x = 0.0;
    PetscReal expected_global_max_x = 0.0;

    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    PetscCall(PicurvAssertBool((PetscBool)(size >= 2), "unit-mpi requires at least two MPI ranks"));

    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 8, 6, 4));
    expected_global_max_x = ((PetscReal)(user->IM - 1) / (PetscReal)user->IM) + 1.0e-6;
    PetscCall(ComputeLocalBoundingBox(user, &local_bbox));
    PetscCall(GatherAllBoundingBoxes(user, &boxes));
    PetscCall(BroadcastAllBoundingBoxes(user, &boxes));
    PetscCall(PicurvAssertBool((PetscBool)(boxes != NULL), "all ranks should hold the gathered bounding-box table"));

    for (PetscMPIInt r = 0; r < size; ++r) {
        PetscCall(PicurvAssertBool((PetscBool)(boxes[r].min_coords.x <= boxes[r].max_coords.x), "bbox x-range should be valid"));
        PetscCall(PicurvAssertBool((PetscBool)(boxes[r].min_coords.y <= boxes[r].max_coords.y), "bbox y-range should be valid"));
        PetscCall(PicurvAssertBool((PetscBool)(boxes[r].min_coords.z <= boxes[r].max_coords.z), "bbox z-range should be valid"));
    }

    PetscCall(PicurvAssertRealNear(local_bbox.min_coords.x, boxes[rank].min_coords.x, 1.0e-10, "local bbox min x preserved"));
    PetscCall(PicurvAssertRealNear(local_bbox.max_coords.x, boxes[rank].max_coords.x, 1.0e-10, "local bbox max x preserved"));

    PetscCallMPI(MPI_Allreduce(&local_bbox.min_coords.x, &global_min_x, 1, MPIU_REAL, MPI_MIN, PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Allreduce(&local_bbox.max_coords.x, &global_max_x, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertBool((PetscBool)(global_min_x <= 0.0), "global min x should include domain start"));
    PetscCall(PicurvAssertRealNear(expected_global_max_x, global_max_x, 1.0e-10, "global max x should match the normalized physical-node domain end"));

    free(boxes);
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests distributed synchronization of persistent cell-centered scalar fields.
 */

static PetscErrorCode TestPeriodicCellFieldSynchronizationMultiRank(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal ***p = NULL;
    PetscReal ***cs = NULL;
    PetscReal ***diffusivity = NULL;
    const char *fields[] = {"P", "CS", "Diffusivity"};
    PetscInt xs, xe, ys, ye, zs, ze, mx;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 8, 4, 4, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
    MarkXPeriodic(user);
    PetscCall(PicurvAssertBool((PetscBool)(user->info.xm < user->info.mx),
                               "periodic cell synchronization test requires the x axis to be partitioned"));
    PetscCall(DMCreateGlobalVector(user->da, &user->CS));
    PetscCall(DMCreateLocalVector(user->da, &user->lCs));

    xs = user->info.xs;
    xe = xs + user->info.xm;
    ys = user->info.ys;
    ye = ys + user->info.ym;
    zs = user->info.zs;
    ze = zs + user->info.zm;
    mx = user->info.mx;

    PetscCall(DMDAVecGetArray(user->da, user->P, &p));
    PetscCall(DMDAVecGetArray(user->da, user->CS, &cs));
    PetscCall(DMDAVecGetArray(user->da, user->Diffusivity, &diffusivity));
    for (PetscInt k = zs; k < ze; ++k) {
        for (PetscInt j = ys; j < ye; ++j) {
            for (PetscInt i = xs; i < xe; ++i) {
                p[k][j][i] = PeriodicScalarValue(i, j, k, 0.0);
                cs[k][j][i] = PeriodicScalarValue(i, j, k, 1000.0);
                diffusivity[k][j][i] = PeriodicScalarValue(i, j, k, 2000.0);
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->da, user->P, &p));
    PetscCall(DMDAVecRestoreArray(user->da, user->CS, &cs));
    PetscCall(DMDAVecRestoreArray(user->da, user->Diffusivity, &diffusivity));

    PetscCall(SynchronizePeriodicCellFields(user, 3, fields));

    PetscCall(DMDAVecGetArrayRead(user->da, user->P, &p));
    PetscCall(DMDAVecGetArrayRead(user->da, user->CS, &cs));
    PetscCall(DMDAVecGetArrayRead(user->da, user->Diffusivity, &diffusivity));
    if (xs == 0) {
        for (PetscInt k = zs; k < ze; ++k) {
            for (PetscInt j = ys; j < ye; ++j) {
                PetscCall(PicurvAssertRealNear(PeriodicScalarValue(mx - 2, j, k, 0.0), p[k][j][0], 1.0e-12,
                                               "distributed P negative endpoint should copy the opposite interior"));
                PetscCall(PicurvAssertRealNear(PeriodicScalarValue(mx - 2, j, k, 1000.0), cs[k][j][0], 1.0e-12,
                                               "distributed CS negative endpoint should copy the opposite interior"));
                PetscCall(PicurvAssertRealNear(PeriodicScalarValue(mx - 2, j, k, 2000.0), diffusivity[k][j][0], 1.0e-12,
                                               "distributed Diffusivity negative endpoint should copy the opposite interior"));
            }
        }
    }
    if (xe == mx) {
        for (PetscInt k = zs; k < ze; ++k) {
            for (PetscInt j = ys; j < ye; ++j) {
                PetscCall(PicurvAssertRealNear(PeriodicScalarValue(1, j, k, 0.0), p[k][j][mx - 1], 1.0e-12,
                                               "distributed P positive endpoint should copy the opposite interior"));
                PetscCall(PicurvAssertRealNear(PeriodicScalarValue(1, j, k, 1000.0), cs[k][j][mx - 1], 1.0e-12,
                                               "distributed CS positive endpoint should copy the opposite interior"));
                PetscCall(PicurvAssertRealNear(PeriodicScalarValue(1, j, k, 2000.0), diffusivity[k][j][mx - 1], 1.0e-12,
                                               "distributed Diffusivity positive endpoint should copy the opposite interior"));
            }
        }
    }
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->P, &p));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->CS, &cs));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->Diffusivity, &diffusivity));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests periodic seam translation discovery when seam nodes are on different ranks.
 */
static PetscErrorCode TestPeriodicGeometryValidationMultiRank(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal expected_translation = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 8, 4, 4, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
    MarkXPeriodic(user);
    PetscCall(PicurvAssertBool((PetscBool)(user->info.xm < user->info.mx),
                               "periodic geometry validation test requires the x axis to be partitioned"));

    expected_translation = (PetscReal)(user->info.mx - 2) / (PetscReal)user->info.mx;
    PetscCall(ValidatePeriodicGeometry(user));

    PetscCall(PicurvAssertBool(user->periodic_translation_valid[0],
                               "distributed X-periodic translation should be marked valid"));
    PetscCall(PicurvAssertRealNear(expected_translation, user->periodic_translation[0].x, 1.0e-12,
                                   "distributed X-periodic translation"));
    PetscCall(PicurvAssertRealNear(0.0, user->periodic_translation[0].y, 1.0e-12,
                                   "distributed X-periodic translation y component"));
    PetscCall(PicurvAssertRealNear(0.0, user->periodic_translation[0].z, 1.0e-12,
                                   "distributed X-periodic translation z component"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests translated face-center ghosts when the periodic axis is partitioned.
 */
static PetscErrorCode TestPeriodicFaceCenterCoordinatesMultiRank(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Cmpnts ***centx = NULL;
    const Cmpnts ***lcentx = NULL;
    const char *fields[] = {"Centx"};
    PetscReal translation, spacing;
    PetscInt xs, xe, mx;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 8, 4, 4, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
    MarkXPeriodic(user);
    PetscCall(ValidatePeriodicGeometry(user));
    PetscCall(PicurvAssertBool((PetscBool)(user->info.xm < user->info.mx),
                               "face-center coordinate test requires the x axis to be partitioned"));
    translation = user->periodic_translation[0].x;
    spacing = translation / (PetscReal)(user->info.mx - 2);
    xs = user->info.xs;
    xe = xs + user->info.xm;
    mx = user->info.mx;

    PetscCall(DMDAVecGetArray(user->fda, user->Centx, &centx));
    for (PetscInt k = user->info.zs; k < user->info.zs + user->info.zm; k++) {
        for (PetscInt j = user->info.ys; j < user->info.ys + user->info.ym; j++) {
            for (PetscInt i = xs; i < xe; i++) {
                centx[k][j][i] = (Cmpnts){spacing * i, 2.0 + j, 3.0 + k};
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Centx, &centx));

    PetscCall(SynchronizePeriodicFaceFields(user, 'i', 1, fields));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCentx, &lcentx));
    if (xs == 0) {
        PetscCall(PicurvAssertRealNear(-spacing, lcentx[2][2][-1].x, 1.0e-12,
                                       "distributed translated Centx negative adjacent ghost"));
        PetscCall(PicurvAssertRealNear(0.0, lcentx[2][2][0].x, 1.0e-12,
                                       "distributed translated Centx negative endpoint"));
    }
    if (xe == mx) {
        PetscCall(PicurvAssertRealNear(translation + spacing, lcentx[2][2][mx - 1].x, 1.0e-12,
                                       "distributed translated Centx positive endpoint"));
        PetscCall(PicurvAssertRealNear(translation + 2.0 * spacing, lcentx[2][2][mx].x, 1.0e-12,
                                       "distributed translated Centx positive adjacent ghost"));
    }
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCentx, &lcentx));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests persistent I-face synchronization across an MPI-partitioned seam.
 */
static PetscErrorCode TestPeriodicFaceFieldSynchronizationMultiRank(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscReal ***iaj = NULL;
    const char *fields[] = {"IAj"};
    PetscInt xs, xe, ys, ye, zs, ze, mx;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 8, 4, 4, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
    MarkXPeriodic(user);
    PetscCall(PicurvAssertBool((PetscBool)(user->info.xm < user->info.mx),
                               "periodic face synchronization test requires the x axis to be partitioned"));

    xs = user->info.xs; xe = xs + user->info.xm;
    ys = user->info.ys; ye = ys + user->info.ym;
    zs = user->info.zs; ze = zs + user->info.zm;
    mx = user->info.mx;

    PetscCall(DMDAVecGetArray(user->da, user->IAj, &iaj));
    for (PetscInt k = zs; k < ze; k++) for (PetscInt j = ys; j < ye; j++) for (PetscInt i = xs; i < xe; i++) {
        iaj[k][j][i] = PeriodicScalarValue(i, j, k, 3000.0);
    }
    PetscCall(DMDAVecRestoreArray(user->da, user->IAj, &iaj));

    PetscCall(SynchronizePeriodicFaceFields(user, 'i', 1, fields));

    PetscCall(DMDAVecGetArrayRead(user->da, user->IAj, &iaj));
    if (xs == 0) {
        for (PetscInt k = zs; k < ze; k++) for (PetscInt j = ys; j < ye; j++) {
            PetscCall(PicurvAssertRealNear(PeriodicScalarValue(mx - 2, j, k, 3000.0), iaj[k][j][0], 1.0e-12,
                                           "distributed I-face negative seam should copy the opposite physical seam"));
        }
    }
    if (xe == mx) {
        for (PetscInt k = zs; k < ze; k++) for (PetscInt j = ys; j < ye; j++) {
            PetscCall(PicurvAssertRealNear(PeriodicScalarValue(1, j, k, 3000.0), iaj[k][j][mx - 1], 1.0e-12,
                                           "distributed I-face positive dummy should copy the leading physical face"));
        }
    }
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->IAj, &iaj));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests persistent Ucont synchronization across an MPI-partitioned seam.
 */
static PetscErrorCode TestPeriodicStaggeredFieldSynchronizationMultiRank(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Cmpnts ***ucont = NULL;
    const char *fields[] = {"Ucont"};
    PetscInt xs, xe, ys, ye, zs, ze, mx;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 8, 4, 4, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
    MarkXPeriodic(user);
    PetscCall(PicurvAssertBool((PetscBool)(user->info.xm < user->info.mx),
                               "periodic staggered synchronization test requires the x axis to be partitioned"));

    xs = user->info.xs; xe = xs + user->info.xm;
    ys = user->info.ys; ye = ys + user->info.ym;
    zs = user->info.zs; ze = zs + user->info.zm;
    mx = user->info.mx;

    PetscCall(DMDAVecGetArray(user->fda, user->Ucont, &ucont));
    for (PetscInt k = zs; k < ze; k++) for (PetscInt j = ys; j < ye; j++) for (PetscInt i = xs; i < xe; i++) {
        PetscReal value = PeriodicScalarValue(i, j, k, 0.0);
        ucont[k][j][i] = (Cmpnts){value + 1000.0, value + 2000.0, value + 3000.0};
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucont, &ucont));

    PetscCall(SynchronizePeriodicStaggeredFields(user, 1, fields));

    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, &ucont));
    if (xs == 0) {
        for (PetscInt k = zs; k < ze; k++) for (PetscInt j = ys; j < ye; j++) {
            PetscReal expected = PeriodicScalarValue(mx - 2, j, k, 0.0);
            PetscCall(PicurvAssertRealNear(expected + 1000.0, ucont[k][j][0].x, 1.0e-12,
                                           "distributed Ucont.x negative seam"));
            PetscCall(PicurvAssertRealNear(expected + 2000.0, ucont[k][j][0].y, 1.0e-12,
                                           "distributed Ucont.y negative X endpoint"));
        }
    }
    if (xe == mx) {
        for (PetscInt k = zs; k < ze; k++) for (PetscInt j = ys; j < ye; j++) {
            PetscReal expected = PeriodicScalarValue(1, j, k, 0.0);
            PetscCall(PicurvAssertRealNear(expected + 2000.0, ucont[k][j][mx - 1].y, 1.0e-12,
                                           "distributed Ucont.y positive X dummy"));
            PetscCall(PicurvAssertRealNear(expected + 3000.0, ucont[k][j][mx - 1].z, 1.0e-12,
                                           "distributed Ucont.z positive X dummy"));
        }
    }
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, &ucont));

    PetscCall(DMDAVecGetArrayRead(user->fda, user->lUcont, &ucont));
    if (xs == 0) {
        for (PetscInt k = zs; k < ze; k++) for (PetscInt j = ys; j < ye; j++) {
            PetscCall(PicurvAssertRealNear(PeriodicScalarValue(mx - 3, j, k, 1000.0), ucont[k][j][-1].x, 1.0e-12,
                                           "distributed Ucont.x negative adjacent ghost repair"));
            PetscCall(PicurvAssertRealNear(PeriodicScalarValue(1, j, k, 2000.0), ucont[k][j][-1].y, 1.0e-12,
                                           "distributed Ucont.y tangential X wraparound"));
        }
    }
    if (xe == mx) {
        for (PetscInt k = zs; k < ze; k++) for (PetscInt j = ys; j < ye; j++) {
            PetscCall(PicurvAssertRealNear(PeriodicScalarValue(2, j, k, 1000.0), ucont[k][j][mx].x, 1.0e-12,
                                           "distributed Ucont.x positive adjacent ghost repair"));
            PetscCall(PicurvAssertRealNear(PeriodicScalarValue(mx - 2, j, k, 2000.0), ucont[k][j][mx].y, 1.0e-12,
                                           "distributed Ucont.y positive tangential X wraparound"));
        }
    }
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lUcont, &ucont));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests restart fast-path migration using preloaded cell ownership metadata.
 */

static PetscErrorCode TestRestartCellIdMigrationMovesParticleToOwner(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    RankCellInfo my_cell_info;
    PetscMPIInt rank = 0, size = 1;
    PetscInt *cell_ids = NULL;
    PetscReal *positions = NULL;
    PetscInt nlocal = 0;

    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    PetscCall(PicurvAssertIntEqual(2, size, "restart migration unit test expects exactly two MPI ranks"));

    PetscCall(PetscMemzero(&my_cell_info, sizeof(my_cell_info)));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 8, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 1, "ske"));
    PetscCall(GetOwnedCellRange(&user->info, 0, &my_cell_info.xs_cell, &my_cell_info.xm_cell));
    PetscCall(GetOwnedCellRange(&user->info, 1, &my_cell_info.ys_cell, &my_cell_info.ym_cell));
    PetscCall(GetOwnedCellRange(&user->info, 2, &my_cell_info.zs_cell, &my_cell_info.zm_cell));
    PetscCall(PetscMalloc1(size, &user->RankCellInfoMap));
    PetscCallMPI(MPI_Allgather(&my_cell_info, sizeof(RankCellInfo), MPI_BYTE,
                               user->RankCellInfoMap, sizeof(RankCellInfo), MPI_BYTE,
                               PETSC_COMM_WORLD));

    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void **)&positions));
    if (rank == 0) {
        cell_ids[0] = user->RankCellInfoMap[1].xs_cell;
        cell_ids[1] = user->RankCellInfoMap[1].ys_cell;
        cell_ids[2] = user->RankCellInfoMap[1].zs_cell;
        positions[0] = 0.875;
        positions[1] = 0.5;
        positions[2] = 0.5;
    } else {
        cell_ids[0] = user->RankCellInfoMap[1].xs_cell;
        cell_ids[1] = user->RankCellInfoMap[1].ys_cell;
        cell_ids[2] = user->RankCellInfoMap[1].zs_cell;
        positions[0] = 0.95;
        positions[1] = 0.5;
        positions[2] = 0.5;
    }
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void **)&positions));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));

    PetscCall(MigrateRestartParticlesUsingCellID(user));
    PetscCall(DMSwarmGetLocalSize(user->swarm, &nlocal));
    PetscCall(PicurvAssertIntEqual((rank == 0) ? 0 : 2, nlocal,
                                   "restart migration should move the foreign particle onto the owning rank"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Runs the unit-mpi PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"distribute-particles-collective-consistency", TestDistributeParticlesCollectiveConsistency},
        {"bounding-box-collectives-multi-rank", TestBoundingBoxCollectivesMultiRank},
        {"periodic-cell-field-synchronization-multi-rank", TestPeriodicCellFieldSynchronizationMultiRank},
        {"periodic-geometry-validation-multi-rank", TestPeriodicGeometryValidationMultiRank},
        {"periodic-face-center-coordinates-multi-rank", TestPeriodicFaceCenterCoordinatesMultiRank},
        {"periodic-face-field-synchronization-multi-rank", TestPeriodicFaceFieldSynchronizationMultiRank},
        {"periodic-staggered-field-synchronization-multi-rank", TestPeriodicStaggeredFieldSynchronizationMultiRank},
        {"restart-cellid-migration-moves-particle-to-owner", TestRestartCellIdMigrationMovesParticleToOwner},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv MPI-focused runtime tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-mpi", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
