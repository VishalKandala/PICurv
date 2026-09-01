/**
 * @file test_mpi_kernels.c
 * @brief C unit tests for MPI collective helper routines.
 */

#include "test_support.h"

#include "ParticleMotion.h"
#include "ParticleSwarm.h"
#include "Boundaries.h"
#include "grid.h"
#include "statistics_target.h"
#include "statistics_accumulator.h"
#include "statistics_window.h"
#include "io.h"
#include "field_catalog.h"
#include "les.h"

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
    Cmpnts ***ucat = NULL;
    const FieldId fields[] = {FIELD_ID_P, FIELD_ID_CS, FIELD_ID_DIFFUSIVITY, FIELD_ID_UCAT};
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
    PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
    for (PetscInt k = zs; k < ze; ++k) {
        for (PetscInt j = ys; j < ye; ++j) {
            for (PetscInt i = xs; i < xe; ++i) {
                p[k][j][i] = PeriodicScalarValue(i, j, k, 0.0);
                cs[k][j][i] = PeriodicScalarValue(i, j, k, 1000.0);
                diffusivity[k][j][i] = PeriodicScalarValue(i, j, k, 2000.0);
                ucat[k][j][i] = (Cmpnts){PeriodicScalarValue(i,j,k,3000.0),0.0,0.0};
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->da, user->P, &p));
    PetscCall(DMDAVecRestoreArray(user->da, user->CS, &cs));
    PetscCall(DMDAVecRestoreArray(user->da, user->Diffusivity, &diffusivity));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));

    PetscCall(SynchronizePeriodicCellFields(user, 4, fields));

    PetscCall(DMDAVecGetArrayRead(user->da, user->P, &p));
    PetscCall(DMDAVecGetArrayRead(user->da, user->CS, &cs));
    PetscCall(DMDAVecGetArrayRead(user->da, user->Diffusivity, &diffusivity));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucat, &ucat));
    if (xs == 0) {
        for (PetscInt k = zs; k < ze; ++k) {
            for (PetscInt j = ys; j < ye; ++j) {
                PetscCall(PicurvAssertRealNear(PeriodicScalarValue(mx - 2, j, k, 0.0), p[k][j][0], 1.0e-12,
                                               "distributed P negative endpoint should copy the opposite interior"));
                PetscCall(PicurvAssertRealNear(PeriodicScalarValue(mx - 2, j, k, 1000.0), cs[k][j][0], 1.0e-12,
                                               "distributed CS negative endpoint should copy the opposite interior"));
                PetscCall(PicurvAssertRealNear(PeriodicScalarValue(mx - 2, j, k, 2000.0), diffusivity[k][j][0], 1.0e-12,
                                               "distributed Diffusivity negative endpoint should copy the opposite interior"));
                PetscCall(PicurvAssertRealNear(PeriodicScalarValue(mx-2,j,k,3000.0),ucat[k][j][0].x,1e-12,
                                               "distributed Ucat negative endpoint should copy the opposite interior"));
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
                PetscCall(PicurvAssertRealNear(PeriodicScalarValue(1,j,k,3000.0),ucat[k][j][mx-1].x,1e-12,
                                               "distributed Ucat positive endpoint should copy the leading interior"));
            }
        }
    }
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->P, &p));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->CS, &cs));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->Diffusivity, &diffusivity));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucat, &ucat));
    PetscCall(DMDAVecGetArrayRead(user->fda,user->lUcat,&ucat));
    if(xs==0) PetscCall(PicurvAssertRealNear(PeriodicScalarValue(1,2,2,3000.0),ucat[2][2][-1].x,1e-12,
                                             "distributed local Ucat negative ghost follows corrected global endpoint"));
    if(xe==mx) PetscCall(PicurvAssertRealNear(PeriodicScalarValue(mx-2,2,2,3000.0),ucat[2][2][mx].x,1e-12,
                                              "distributed local Ucat positive ghost follows corrected global endpoint"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda,user->lUcat,&ucat));

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
    const FieldId fields[] = {FIELD_ID_CENTX};
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
    const FieldId fields[] = {FIELD_ID_IAJ};
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
    const FieldId fields[] = {FIELD_ID_UCONT};
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
 * @brief Spatial target domains must not depend on how ranks divide the block.
 *
 * The global point count is a function of layout, global dimensions, and
 * periodicity alone. Summing the per-rank counts must therefore reproduce the
 * same analytic totals a single rank computes, with no point counted twice at a
 * decomposition seam and none dropped.
 */
static PetscErrorCode TestSpatialTargetIsDecompositionIndependent(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    SpatialTargetPlan plan;
    PetscInt local = 0, global = 0, rank_count = 0;
    PetscMPIInt size = 0;
    const PetscInt n = 6;          /* DMDA is n + 1 = 7 per dimension */
    const PetscInt cell = 5;       /* cell-like span [1, 6) */
    const PetscInt node = 6;       /* node-like nonperiodic span [0, 6) */

    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, n, n, n));

    PetscCall(SpatialTargetPlanCreate(user, FIELD_ID_UCAT, PICURV_STATISTICS_MASK_FLUID, &plan));
    PetscCall(SpatialTargetPlanLocalPointCount(&plan, &local));
    PetscCall(SpatialTargetPlanGlobalPointCount(&plan, PETSC_COMM_WORLD, &global));
    PetscCall(PicurvAssertIntEqual(cell * cell * cell, global,
                                   "cell-centered global count must be decomposition independent"));

    /* Every rank must contribute a non-negative count, and the ranks together
     * must contribute the whole domain exactly once. */
    PetscCall(PicurvAssertBool((PetscBool)(local >= 0), "local point count must be non-negative"));
    PetscCallMPI(MPI_Allreduce(&local, &rank_count, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertIntEqual(global, rank_count,
                                   "summed local counts must equal the reduced global count"));

    PetscCall(SpatialTargetPlanCreate(user, FIELD_ID_CSI, PICURV_STATISTICS_MASK_FLUID, &plan));
    PetscCall(SpatialTargetPlanGlobalPointCount(&plan, PETSC_COMM_WORLD, &global));
    PetscCall(PicurvAssertIntEqual(node * cell * cell, global,
                                   "I-face global count must be decomposition independent"));

    PetscCall(SpatialTargetPlanCreate(user, FIELD_ID_COORDINATES, PICURV_STATISTICS_MASK_FLUID, &plan));
    PetscCall(SpatialTargetPlanGlobalPointCount(&plan, PETSC_COMM_WORLD, &global));
    PetscCall(PicurvAssertIntEqual(node * node * node, global,
                                   "node-centered global count must be decomposition independent"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));

    /* Fully periodic: duplicate planes must be dropped once globally, not once per rank. */
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, n, n, n,
                                                         PETSC_TRUE, PETSC_TRUE, PETSC_TRUE));
    PetscCall(SpatialTargetPlanCreate(user, FIELD_ID_COORDINATES, PICURV_STATISTICS_MASK_FLUID, &plan));
    PetscCall(SpatialTargetPlanGlobalPointCount(&plan, PETSC_COMM_WORLD, &global));
    PetscCall(PicurvAssertIntEqual(cell * cell * cell, global,
                                   "fully periodic node-centered count must drop duplicates exactly once"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Statistics payloads must survive a checkpoint round trip under decomposition.
 *
 * The write path gathers each accumulator into natural ordering before it reaches
 * disk, and the read path scatters it back. That pairing is where a decomposition
 * bug would live, so this drives the whole cycle on more than one rank and requires
 * an exact match: the format stores IEEE doubles, and natural ordering is defined
 * to be independent of how ranks divide the block.
 */
static PetscErrorCode TestStatisticsCheckpointRoundTripMultiRank(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PicurvWindow window;
    PicurvWindowStorage storage;
    PicurvWindowDefinition definition;
    Vec *reference = NULL;
    char tmpdir[PETSC_MAX_PATH_LEN];
    PetscInt payload_count = 0;
    PetscMPIInt size = 0;

    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    PetscCall(PicurvAssertBool((PetscBool)(size >= 2), "unit-mpi requires at least two MPI ranks"));

    memset(&definition, 0, sizeof(definition));
    strncpy(definition.name, "production", PICURV_WINDOW_NAME_LENGTH - 1);
    definition.weighting = PICURV_WEIGHTING_PHYSICAL_TIME;
    definition.cadence_kind = PICURV_CADENCE_STEP;
    definition.step_cadence = 1;
    definition.field_count = 2;
    definition.fields[0].field_id = FIELD_ID_UCAT; definition.fields[0].want_second = PETSC_TRUE;
    definition.fields[1].field_id = FIELD_ID_P;    definition.fields[1].want_second = PETSC_TRUE;
    definition.covariance_count = 1;
    definition.covariances[0].first = FIELD_ID_UCAT;
    definition.covariances[0].second = FIELD_ID_P;

    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    tmpdir[0] = '\0';
    if (simCtx->rank == 0) PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCallMPI(MPI_Bcast(tmpdir, sizeof(tmpdir), MPI_CHAR, 0, PETSC_COMM_WORLD));
    PetscCall(PetscStrncpy(simCtx->output_dir, tmpdir, sizeof(simCtx->output_dir)));
    PetscCall(PetscStrncpy(simCtx->restart_dir, tmpdir, sizeof(simCtx->restart_dir)));
    PetscCall(PicurvPopulateIdentityMetrics(user));
    PetscCall(VecSet(user->Nvert, 0.0));

    PetscCall(PicurvWindowInit(&window, &definition));
    PetscCall(PicurvWindowStorageCreate(user, &definition, &storage));
    simCtx->fieldStatisticsEnabled = PETSC_TRUE;
    simCtx->fieldStatisticsWindowCount = 1;
    simCtx->fieldStatisticsWindows = &window;
    user->fieldStatisticsStorage = &storage;

    /* A spatially varying field makes the assertion meaningful: a payload gathered
     * or scattered through the wrong owner comes back permuted, which a uniform
     * field could not reveal. */
    for (PetscInt s = 0; s < 4; ++s) {
        PetscBool accepted = PETSC_FALSE;
        PetscReal weight = 0.0;
        PetscReal ***p = NULL;
        Cmpnts ***ucat = NULL;
        const DMDALocalInfo info = user->info;

        PetscCall(DMDAVecGetArray(user->da, user->P, &p));
        PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
        for (PetscInt k = info.zs; k < info.zs + info.zm; ++k)
            for (PetscInt j = info.ys; j < info.ys + info.ym; ++j)
                for (PetscInt i = info.xs; i < info.xs + info.xm; ++i) {
                    const PetscReal tag = (PetscReal)(100 * k + 10 * j + i);

                    p[k][j][i] = tag + 0.25 * (PetscReal)s;
                    ucat[k][j][i].x = tag;
                    ucat[k][j][i].y = tag + 1.0 + (PetscReal)s;
                    ucat[k][j][i].z = tag - 1.0;
                }
        PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));
        PetscCall(DMDAVecRestoreArray(user->da, user->P, &p));

        PetscCall(PicurvWindowOfferState(&window, s, 0.5 * (PetscReal)s, &accepted, &weight));
        if (accepted) PetscCall(PicurvWindowAccumulate(user, &definition, &storage, weight));
    }
    PetscCall(PicurvAssertIntEqual(3, window.sample_count, "three intervals are represented"));

    PetscCall(WriteCheckpointBundle(simCtx, "test"));

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

    simCtx->fieldStatisticsContinue = PETSC_TRUE;
    PetscCall(RestoreFieldStatisticsState(simCtx, simCtx->step));
    PetscCall(PicurvAssertIntEqual(3, window.sample_count, "window bookkeeping survives on every rank"));
    PetscCall(PicurvAssertRealNear(1.5, window.total_weight, 1.0e-12,
                                   "physical-time weights survive on every rank"));

    for (PetscInt index = 0; index < payload_count; ++index) {
        PicurvStatisticsPayload payload;
        PetscReal difference = 0.0;
        char context[192];

        PetscCall(PicurvWindowStoragePayload(user, &definition, &storage, index, &payload));
        PetscCall(VecAXPY(reference[index], -1.0, payload.vec));
        PetscCall(VecNorm(reference[index], NORM_INFINITY, &difference));
        PetscCall(PetscSNPrintf(context, sizeof(context),
                                "payload '%s' round trips exactly across %d ranks",
                                payload.name, (int)size));
        PetscCall(PicurvAssertBool((PetscBool)(difference == 0.0), context));
        PetscCall(VecDestroy(&reference[index]));
    }
    PetscCall(PetscFree(reference));

    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    if (simCtx->rank == 0) PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));
    simCtx->fieldStatisticsEnabled = PETSC_FALSE;
    simCtx->fieldStatisticsWindowCount = 0;
    simCtx->fieldStatisticsWindows = NULL;
    user->fieldStatisticsStorage = NULL;
    PetscCall(PicurvWindowStorageDestroy(&storage));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that LES coefficient averaging does not depend on the decomposition.
 *
 * The averaged coefficient is a physical property of the block, so splitting the same
 * field across more ranks must not change it. This case catches a wrong communicator,
 * a halo entry counted as if it were owned, or a layout boundary plane included in one
 * decomposition but not another. The reference is computed from the analytic
 * definition over the interior index box rather than from a single-rank run, so the
 * test does not merely compare one decomposition against itself.
 */
static PetscErrorCode TestLESAveragingIsDecompositionIndependent(void)
{
    SimCtx      *simCtx = NULL;
    UserCtx     *user = NULL;
    Vec          numerator, denominator, ratio;
    SpatialTargetPlan plan;
    PetscReal ***num = NULL, ***den = NULL, ***out = NULL;
    DMDALocalInfo info;
    PetscInt      xs, xe, ys, ye, zs, ze;
    PetscReal     sampled = 0.0;
    PetscReal     numerator_total = 0.0, denominator_total = 0.0;
    const PetscBool all[3] = {PETSC_TRUE, PETSC_TRUE, PETSC_TRUE};

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 8, 8, 8));
    PetscCall(DMDAGetLocalInfo(user->da, &info));
    PetscCall(VecSet(user->lNvert, 0.0));
    PetscCall(VecSet(user->lAj, 1.0));

    xs = (info.xs == 0) ? 1 : info.xs;
    ys = (info.ys == 0) ? 1 : info.ys;
    zs = (info.zs == 0) ? 1 : info.zs;
    xe = (info.xs + info.xm == info.mx) ? info.mx - 1 : info.xs + info.xm;
    ye = (info.ys + info.ym == info.my) ? info.my - 1 : info.ys + info.ym;
    ze = (info.zs + info.zm == info.mz) ? info.mz - 1 : info.zs + info.zm;

    PetscCall(VecDuplicate(user->lNvert, &numerator));
    PetscCall(VecDuplicate(user->lNvert, &denominator));
    PetscCall(VecDuplicate(user->lNvert, &ratio));
    PetscCall(VecSet(numerator, 0.0));
    PetscCall(VecSet(denominator, 0.0));
    PetscCall(VecSet(ratio, 0.0));

    /* Each rank writes only the cells it owns, from a closed-form function of the
       global indices, so the distributed sums are exactly the serial ones below. */
    PetscCall(DMDAVecGetArray(user->da, numerator, &num));
    PetscCall(DMDAVecGetArray(user->da, denominator, &den));
    for (PetscInt k = zs; k < ze; ++k)
    for (PetscInt j = ys; j < ye; ++j)
    for (PetscInt i = xs; i < xe; ++i) {
        num[k][j][i] = (PetscReal)(i * j + k);
        den[k][j][i] = (PetscReal)(i + 1);
    }
    PetscCall(DMDAVecRestoreArray(user->da, numerator, &num));
    PetscCall(DMDAVecRestoreArray(user->da, denominator, &den));

    PetscCall(SpatialTargetPlanCreate(user, FIELD_ID_CS,
                                      PICURV_STATISTICS_MASK_FLUID, &plan));
    PetscCall(PicurvSpatialRatioAverage(user, &plan, numerator, denominator, NULL, all,
                                        PETSC_COMM_WORLD, ratio, NULL));

    for (PetscInt k = 1; k < info.mz - 1; ++k)
    for (PetscInt j = 1; j < info.my - 1; ++j)
    for (PetscInt i = 1; i < info.mx - 1; ++i) {
        numerator_total   += (PetscReal)(i * j + k);
        denominator_total += (PetscReal)(i + 1);
    }

    PetscCall(DMDAVecGetArray(user->da, ratio, &out));
    sampled = out[zs][ys][xs];
    PetscCall(DMDAVecRestoreArray(user->da, ratio, &out));

    PetscCall(PicurvAssertRealNear(numerator_total / denominator_total, sampled, 1.0e-10,
                                   "the globally averaged ratio must match the serial sums "
                                   "however the block is decomposed"));

    PetscCall(VecDestroy(&numerator));
    PetscCall(VecDestroy(&denominator));
    PetscCall(VecDestroy(&ratio));
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
        {"spatial-target-decomposition-independent", TestSpatialTargetIsDecompositionIndependent},
        {"statistics-checkpoint-round-trip-multi-rank", TestStatisticsCheckpointRoundTripMultiRank},
        {"les-averaging-decomposition-independent", TestLESAveragingIsDecompositionIndependent},
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
