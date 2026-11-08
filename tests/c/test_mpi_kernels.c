/**
 * @file test_mpi_kernels.c
 * @brief C unit tests for MPI collective helper routines.
 */

#include "test_support.h"

#include "ParticleMotion.h"
#include "ParticleSwarm.h"
#include "ParticlePhysics.h"
#include "runloop.h"
#include "expression.h"
#include "initialcondition.h"
#include "Boundaries.h"
#include "grid.h"
#include "Metric.h"
#include "walkingsearch.h"
#include "setup.h"
#include "statistics_target.h"
#include "statistics_accumulator.h"
#include "statistics_window.h"
#include "io.h"
#include "field_catalog.h"
#include "les.h"
#include "logging.h"

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
/** @brief Width of cell n (1-based interior index) along an axis in the line-fraction test. */
static PetscReal LineFractionTestWidth(PetscInt axis, PetscInt n)
{
    if (axis == 0) return PetscPowReal(1.3, (PetscReal)n);   /* geometric stretching */
    if (axis == 1) return 0.5;                              /* uniform */
    return 1.0 + 0.1 * (PetscReal)n;                        /* linear growth */
}

/**
 * @brief Grid-line fractions follow arc length on a stretched orthogonal grid, and the
 *        lines that cross ranks are summed in order.
 */
static PetscErrorCode TestGridLineFractionsFollowArcLength(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Vec fraction, length;
    Cmpnts ***csi, ***eta, ***zet;
    const Cmpnts ***frac, ***len;
    PetscReal ***aj;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 9, 7, 8));
    const DMDALocalInfo info = user->info;
    const PetscInt mx = info.mx, my = info.my, mz = info.mz, n_nodes[3] = {mx, my, mz};

    /* An orthogonal box cell of widths (dx, dy, dz) has face-area vectors (dy dz, 0, 0),
       (0, dx dz, 0), (0, 0, dx dy) and inverse volume 1/(dx dy dz). */
    PetscCall(DMDAVecGetArray(user->fda, user->Csi, &csi));
    PetscCall(DMDAVecGetArray(user->fda, user->Eta, &eta));
    PetscCall(DMDAVecGetArray(user->fda, user->Zet, &zet));
    PetscCall(DMDAVecGetArray(user->da, user->Aj, &aj));
    for (PetscInt k = info.zs; k < info.zs + info.zm; ++k)
    for (PetscInt j = info.ys; j < info.ys + info.ym; ++j)
    for (PetscInt i = info.xs; i < info.xs + info.xm; ++i) {
        const PetscReal dx = LineFractionTestWidth(0, i), dy = LineFractionTestWidth(1, j), dz = LineFractionTestWidth(2, k);
        csi[k][j][i] = (Cmpnts){dy * dz, 0.0, 0.0};
        eta[k][j][i] = (Cmpnts){0.0, dx * dz, 0.0};
        zet[k][j][i] = (Cmpnts){0.0, 0.0, dx * dy};
        aj[k][j][i] = 1.0 / (dx * dy * dz);
    }
    PetscCall(DMDAVecRestoreArray(user->da, user->Aj, &aj));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Zet, &zet));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Eta, &eta));
    PetscCall(DMDAVecRestoreArray(user->fda, user->Csi, &csi));
    PetscCall(DMGlobalToLocal(user->fda, user->Csi, INSERT_VALUES, user->lCsi));
    PetscCall(DMGlobalToLocal(user->fda, user->Eta, INSERT_VALUES, user->lEta));
    PetscCall(DMGlobalToLocal(user->fda, user->Zet, INSERT_VALUES, user->lZet));
    PetscCall(DMGlobalToLocal(user->da, user->Aj, INSERT_VALUES, user->lAj));

    PetscCall(DMCreateGlobalVector(user->fda, &fraction));
    PetscCall(DMCreateGlobalVector(user->fda, &length));
    PetscCall(ComputeGridLineFractions(user, fraction, length));
    PetscCall(DMDAVecGetArrayRead(user->fda, fraction, &frac));
    PetscCall(DMDAVecGetArrayRead(user->fda, length, &len));
    for (PetscInt k = PetscMax(1, info.zs); k < PetscMin(mz - 1, info.zs + info.zm); ++k)
    for (PetscInt j = PetscMax(1, info.ys); j < PetscMin(my - 1, info.ys + info.ym); ++j)
    for (PetscInt i = PetscMax(1, info.xs); i < PetscMin(mx - 1, info.xs + info.xm); ++i) {
        const PetscInt idx[3] = {i, j, k};
        const PetscReal got_f[3] = {frac[k][j][i].x, frac[k][j][i].y, frac[k][j][i].z};
        const PetscReal got_l[3] = {len[k][j][i].x, len[k][j][i].y, len[k][j][i].z};
        for (PetscInt axis = 0; axis < 3; ++axis) {
            PetscReal before = 0.0, total = 0.0;
            for (PetscInt n = 1; n < n_nodes[axis] - 1; ++n) {
                if (n < idx[axis]) before += LineFractionTestWidth(axis, n);
                total += LineFractionTestWidth(axis, n);
            }
            const PetscReal expected = (before + 0.5 * LineFractionTestWidth(axis, idx[axis])) / total;
            PetscCall(PicurvAssertRealNear(expected, got_f[axis], 1.0e-12, "arc-length fraction"));
            PetscCall(PicurvAssertRealNear(total, got_l[axis], 1.0e-12, "grid-line length"));
        }
    }
    PetscCall(DMDAVecRestoreArrayRead(user->fda, length, &len));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, fraction, &frac));
    PetscCall(VecDestroy(&fraction));
    PetscCall(VecDestroy(&length));
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
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 8, 6, 6, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
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
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 8, 6, 6, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
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
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 8, 6, 6, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
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
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 8, 6, 6, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
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
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 8, 6, 6, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
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
 * @brief A profiled function that takes a measurable time; the profiling macros name it by __func__.
 */
static PetscErrorCode ProfiledOnRankOneOnly(void)
{
    PetscLogDouble start = 0.0, now = 0.0;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    /* Spin until the clock advances: a zero elapsed time ties both ranks and the maximum falls to rank 0. */
    PetscCall(PetscTime(&start));
    do PetscCall(PetscTime(&now)); while (now <= start);
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief The final profiling summary reduces each function's total time over ranks.
 *
 * A function profiled on rank 1 only must reach rank 0's summary with its maximum on rank 1,
 * its minimum (zero) on rank 0 and one call in total.
 */
static PetscErrorCode TestProfilingSummaryReducesOverRanks(void)
{
    SimCtx simCtx;
    PetscMPIInt rank = 0, size = 1;
    char tmpdir[PETSC_MAX_PATH_LEN], summary_path[PETSC_MAX_PATH_LEN + 64], line[1024];

    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    PetscCall(PicurvAssertBool((PetscBool)(size >= 2), "the profiling summary test needs a second rank"));
    PetscCall(PetscMemzero(&simCtx, sizeof(simCtx)));
    if (rank == 0) PetscCall(PicurvMakeTempDir(tmpdir, sizeof(tmpdir)));
    PetscCallMPI(MPI_Bcast(tmpdir, PETSC_MAX_PATH_LEN, MPI_CHAR, 0, PETSC_COMM_WORLD));
    PetscCall(PetscStrncpy(simCtx.log_dir, tmpdir, sizeof(simCtx.log_dir)));
    simCtx.rank = rank;
    simCtx.exec_mode = EXEC_MODE_SOLVER;
    simCtx.profilingFinalSummary = PETSC_TRUE;

    if (rank == 1) PetscCall(ProfiledOnRankOneOnly());
    PetscCall(ProfilingFinalize(&simCtx));

    if (rank == 0) {
        PetscBool found = PETSC_FALSE;
        FILE *f;
        PetscCall(PetscSNPrintf(summary_path, sizeof(summary_path), "%s/ProfilingSummary_Solver.csv", tmpdir));
        f = fopen(summary_path, "r");
        PetscCheck(f, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Cannot open %s", summary_path);
        while (fgets(line, sizeof(line), f)) {
            char name[256];
            double max_time, mean_time, min_time;
            int max_rank, min_rank;
            long long calls;
            double avg_ms, min_calls, max_calls;
            int min_calls_rank, max_calls_rank;
            if (sscanf(line, "picurv,%255[^,],%lf,%d,%lf,%lf,%d,%lld,%lf,%lf,%d,%lf,%d", name, &max_time,
                       &max_rank, &mean_time, &min_time, &min_rank, &calls, &avg_ms, &min_calls, &min_calls_rank,
                       &max_calls, &max_calls_rank) == 12 && strcmp(name, "ProfiledOnRankOneOnly") == 0) {
                found = PETSC_TRUE;
                PetscCall(PicurvAssertIntEqual(1, max_rank, "the maximum should come from the rank that profiled it"));
                PetscCall(PicurvAssertIntEqual(0, min_rank, "the minimum should come from the rank that did not"));
                PetscCall(PicurvAssertRealNear(0.0, min_time, 0.0, "a rank that never called it contributes zero time"));
                PetscCall(PicurvAssertIntEqual(1, (PetscInt)calls, "calls should be summed over ranks"));
                PetscCall(PicurvAssertIntEqual(0, (PetscInt)min_calls, "the fewest calls per rank should be the uncalled rank's zero"));
                PetscCall(PicurvAssertIntEqual(0, min_calls_rank, "the fewest calls should come from rank 0"));
                PetscCall(PicurvAssertIntEqual(1, (PetscInt)max_calls, "the most calls per rank should be rank 1's one"));
                PetscCall(PicurvAssertIntEqual(1, max_calls_rank, "the most calls should come from rank 1"));
            }
        }
        fclose(f);
        PetscCall(PicurvAssertBool(found, "the summary should list a function profiled on another rank"));
        PetscCall(PicurvRemoveTempDir(tmpdir));
    }
    PetscFunctionReturn(0);
}
/**
 * @brief The per-axis owner lookup agrees with a scan of every rank's cell ranges.
 *
 * Checks every cell of the grid, plus cells just outside each face, which have no owner.
 */
static PetscErrorCode TestFindOwnerOfCellMatchesRankRanges(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscMPIInt size = 1;

    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 8, 6, 4));
    PetscCall(SetupDomainCellDecompositionMap(user));

    for (PetscInt k = -1; k <= user->KM - 1; ++k) {
        for (PetscInt j = -1; j <= user->JM - 1; ++j) {
            for (PetscInt i = -1; i <= user->IM - 1; ++i) {
                PetscMPIInt owner, expected = -1;
                for (PetscMPIInt r = 0; r < size; ++r) {
                    const RankCellInfo *info = &user->RankCellInfoMap[r];
                    if (i >= info->xs_cell && i < info->xs_cell + info->xm_cell &&
                        j >= info->ys_cell && j < info->ys_cell + info->ym_cell &&
                        k >= info->zs_cell && k < info->zs_cell + info->zm_cell) expected = r;
                }
                PetscCall(FindOwnerOfCell(user, i, j, k, &owner));
                PetscCall(PicurvAssertIntEqual(expected, owner, "owner lookup should match the rank cell ranges"));
            }
        }
    }

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
    PetscMPIInt rank = 0, size = 1;
    PetscInt *cell_ids = NULL;
    PetscReal *positions = NULL;
    PetscInt nlocal = 0;

    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    /* Two ranks only: it places particles on rank 0 for rank 1 to own. */
    if (size != 2) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "     [SKIP] needs exactly two ranks\n"));
        PetscFunctionReturn(0);
    }

    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 8, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, 1, "ske"));
    PetscCall(SetupDomainCellDecompositionMap(user));

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
 * @brief A particle that left the domain is classified once and never re-searched in later passes.
 *
 * Rank 0 holds one particle beyond the -Xi face and one that has moved into rank 1, so the
 * settlement needs a second pass. The leaving particle must end EXITED after pass 1 and be
 * skipped afterwards: re-searching it is what used to count it as lost once per pass and, on
 * curved boundaries, migrate it between ranks until the pass limit aborted the run.
 */
static PetscErrorCode TestSettlementClassifiesExitOnceAcrossPasses(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscMPIInt rank = 0, size = 1;
    PetscInt *cell_ids = NULL, *status = NULL;
    PetscReal *positions = NULL;
    PetscInt nlocal = 0;
    PetscInt64 local_counts[2], global_counts[2];

    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    /* Two ranks only: it moves a particle from rank 0 to rank 1. */
    if (size != 2) {
        PetscCall(PetscPrintf(PETSC_COMM_WORLD, "     [SKIP] needs exactly two ranks\n"));
        PetscFunctionReturn(0);
    }

    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 8, 4, 4));
    PetscCall(PicurvCreateSwarmPair(user, (rank == 0) ? 2 : 0, "ske"));
    PetscCall(SetupDomainCellDecompositionMap(user));

    if (rank == 0) {
        PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
        PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void **)&positions));
        PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
        /* Particle 0 left through -Xi from boundary cell (0,1,1), well inside the blocked-walk guard. */
        cell_ids[0] = 0; cell_ids[1] = 1; cell_ids[2] = 1;
        positions[0] = -0.05; positions[1] = 0.5; positions[2] = 0.5;
        /* Particle 1 moved from rank 0's last owned cell into rank 1's subdomain. */
        cell_ids[3] = user->RankCellInfoMap[rank].xs_cell + user->RankCellInfoMap[rank].xm_cell - 1; cell_ids[4] = 1; cell_ids[5] = 1;
        positions[3] = 0.875; positions[4] = 0.5; positions[5] = 0.5;
        status[0] = NEEDS_LOCATION;
        status[1] = NEEDS_LOCATION;
        PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
        PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void **)&positions));
        PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_ids));
    }

    PetscCall(LocateAllParticlesInGrid(user, simCtx->bboxlist));

    PetscCall(PicurvAssertIntEqual(2, simCtx->migrationPassesLastStep,
                                   "one migration and one confirming pass should settle both particles"));
    PetscCall(DMSwarmGetLocalSize(user->swarm, &nlocal));
    PetscCall(PicurvAssertIntEqual(1, nlocal, "each rank should hold one particle after settlement"));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
    PetscCall(PicurvAssertIntEqual((rank == 0) ? EXITED : ACTIVE_AND_LOCATED, status[0],
                                   "the leaving particle should be EXITED and the migrated one located"));
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));

    local_counts[0] = simCtx->searchMetrics.searchExitedCount;
    local_counts[1] = simCtx->searchMetrics.searchLostCount;
    PetscCallMPI(MPI_Allreduce(local_counts, global_counts, 2, MPIU_INT64, MPI_SUM, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertIntEqual(1, (PetscInt)global_counts[0], "the exit should be counted once, not once per pass"));
    PetscCall(PicurvAssertIntEqual(0, (PetscInt)global_counts[1], "no search failed"));

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
    /* Psi is particle-sampled, so the window also carries the occupied sample set and
     * Ucat's occupied-state mean for the Ucat-Psi pair; every payload must round-trip. */
    definition.field_count = 3;
    definition.fields[0].field_id = FIELD_ID_UCAT; definition.fields[0].want_second = PETSC_TRUE;
    definition.fields[1].field_id = FIELD_ID_P;    definition.fields[1].want_second = PETSC_TRUE;
    definition.fields[2].field_id = FIELD_ID_PSI;  definition.fields[2].want_second = PETSC_TRUE;
    definition.covariance_count = 2;
    definition.covariances[0].first = FIELD_ID_UCAT;
    definition.covariances[0].second = FIELD_ID_P;
    definition.covariances[1].first = FIELD_ID_UCAT;
    definition.covariances[1].second = FIELD_ID_PSI;

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
        PetscReal ***p = NULL, ***psi = NULL, ***particles = NULL;
        Cmpnts ***ucat = NULL;
        const DMDALocalInfo info = user->info;

        PetscCall(DMDAVecGetArray(user->da, user->P, &p));
        PetscCall(DMDAVecGetArray(user->da, user->Psi, &psi));
        PetscCall(DMDAVecGetArray(user->da, user->ParticleCount, &particles));
        PetscCall(DMDAVecGetArray(user->fda, user->Ucat, &ucat));
        for (PetscInt k = info.zs; k < info.zs + info.zm; ++k)
            for (PetscInt j = info.ys; j < info.ys + info.ym; ++j)
                for (PetscInt i = info.xs; i < info.xs + info.xm; ++i) {
                    const PetscReal tag = (PetscReal)(100 * k + 10 * j + i);

                    p[k][j][i] = tag + 0.25 * (PetscReal)s;
                    /* Occupancy varies with both the point and the state. */
                    particles[k][j][i] = ((i + j + k + s) % 3) ? 2.0 : 0.0;
                    psi[k][j][i] = 0.01 * tag + (PetscReal)s;
                    ucat[k][j][i].x = tag;
                    ucat[k][j][i].y = tag + 1.0 + (PetscReal)s;
                    ucat[k][j][i].z = tag - 1.0;
                }
        PetscCall(DMDAVecRestoreArray(user->fda, user->Ucat, &ucat));
        PetscCall(DMDAVecRestoreArray(user->da, user->ParticleCount, &particles));
        PetscCall(DMDAVecRestoreArray(user->da, user->Psi, &psi));
        PetscCall(DMDAVecRestoreArray(user->da, user->P, &p));

        PetscCall(PicurvWindowOfferState(&window, s, 0.5 * (PetscReal)s, &accepted, &weight));
        if (accepted) PetscCall(PicurvWindowAccumulate(user, &definition, &storage, weight));
    }
    PetscCall(PicurvAssertIntEqual(3, window.sample_count, "three intervals are represented"));

    PetscCall(WriteCheckpointBundle(simCtx, "test"));

    PetscCall(PicurvWindowStoragePayloadCount(&storage, &payload_count));
    /* Six occupancy vectors, three means, Ucat's occupied-state mean, three second
     * moments and two co-moments. */
    PetscCall(PicurvAssertIntEqual(15, payload_count, "the occupied sample set adds its own payloads"));
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
 * @brief Tests that injection follows the flux, is keyed on the grid, and carries no state.
 * @details An empty domain, with a different flux through each inlet face cell and outflow through
 *          every third. The global count and identities are computed from the formula over all
 *          face cells, so they hold whatever the decomposition. Each particle carries the value
 *          at the face, sits within u dt of it, awaits location, and a second call for the same
 *          step reproduces the first. With a diffusivity, each particle moves further by exactly
 *          its own keyed Brownian displacement for the step.
 */
static PetscErrorCode TestInjectionFollowsFluxAndIsKeyedOnTheGrid(void)
{
    const char *options = "-particle_inlet_count 1\n"
                          "-particle_inlet_0_face 4\n"
                          "-particle_inlet_0_particles_per_cell \"16.0 + 0*t\"\n"
                          "-particle_inlet_0_fields_count 1\n"
                          "-particle_inlet_0_fields_0_name Psi\n"
                          "-particle_inlet_0_fields_0_expr_0 \"z + 1\"\n";
    SimCtx    *simCtx = NULL;
    UserCtx   *user = NULL;
    char       tmpdir[PETSC_MAX_PATH_LEN];
    Cmpnts  ***ucont = NULL;
    PetscReal  positions_first = 0.0, positions_second = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvBuildTinyRuntimeContextWithOptions(NULL, PETSC_FALSE, options, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(PicurvAssertBool((PetscBool)(simCtx->np == 0 && ParticlesEnabled(simCtx)), "an inlet enables particles in an empty domain"));
    PetscCall(InitializeEulerianState(simCtx));
    PetscCall(InitializeParticleSwarm(simCtx));

    const PetscInt   nx = user->info.mx - 2, ny = user->info.my - 2, nz = user->info.mz - 2;
    const PetscInt64 faces = 2 * ((PetscInt64)nx * ny + (PetscInt64)ny * nz + (PetscInt64)nz * nx);
    const PetscInt64 face_start = 2 * ((PetscInt64)ny * nz + (PetscInt64)nz * nx);
    const PetscReal  L = simCtx->scaling.L_ref;
    simCtx->step = 3;
    simCtx->dt = 0.1;

    /* A flux per face cell that the decomposition cannot change: outflow through every third. */
#define TEST_INLET_FLUX(i_, j_) (((i_) % 3 == 0) ? -0.01 : 0.002 * ((i_) + 2 * (j_)))
    PetscCall(VecSet(user->lUcont, 0.0));
    PetscCall(DMDAVecGetArray(user->fda, user->lUcont, &ucont));
    if (user->info.gzs <= 0) {
        for (PetscInt j = user->info.gys; j < user->info.gys + user->info.gym; ++j)
        for (PetscInt i = user->info.gxs; i < user->info.gxs + user->info.gxm; ++i) {
            ucont[0][j][i].z = TEST_INLET_FLUX(i, j);
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, user->lUcont, &ucont));

    PetscInt64 expected_count = 0, expected_pid_sum = 0;
    for (PetscInt j = 1; j <= ny; ++j)
    for (PetscInt i = 1; i <= nx; ++i) {
        const PetscInt64 id = simCtx->step * faces + face_start + (i - 1) + (PetscInt64)nx * (j - 1);
        const PetscReal  mean = PetscMax(TEST_INLET_FLUX(i, j), 0.0) * 16.0 * (PetscReal)(nx * ny * nz) * simCtx->dt; /* a unit box */
        const PetscReal  whole = PetscFloorReal(mean);
        const uint64_t   hash = PicurvDrawHash((PetscInt64)simCtx->particleRandomSeed, PICURV_DRAW_SPACE_INLET, 0, id, 0);
        const PetscInt   n = (PetscInt)whole + ((PicurvUnitFromHash(hash) < mean - whole) ? 1 : 0);
        for (PetscInt k = 0; k < n; ++k) expected_pid_sum += ((PetscInt64)1 << 40) + id * 1024 + k;
        expected_count += n;
    }
    PetscCall(PicurvAssertBool((PetscBool)(expected_count > 0), "the fixture must inject something"));

    for (PetscInt pass = 0; pass < 2; ++pass) {
        PetscInt          nlocal = 0, nglobal = 0;
        PetscInt64        pid_sum = 0;
        PetscReal         position_sum = 0.0;
        const PetscReal  *position = NULL, *psi = NULL;
        const PetscInt64 *pid = NULL;
        const PetscInt   *status = NULL;

        PetscCall(DMSwarmSetLocalSizes(user->swarm, 0, -1));
        PetscCall(InjectParticlesAtInlets(user));
        PetscCall(DMSwarmGetLocalSize(user->swarm, &nlocal));
        PetscCall(DMSwarmGetSize(user->swarm, &nglobal));
        PetscCall(PicurvAssertIntEqual((PetscInt)expected_count, nglobal, "the count follows the flux through every face cell"));
        PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void **)&position));
        PetscCall(DMSwarmGetField(user->swarm, "Psi", NULL, NULL, (void **)&psi));
        PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid));
        PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
        for (PetscInt p = 0; p < nlocal; ++p) {
            PetscCall(PicurvAssertBool((PetscBool)(status[p] == NEEDS_LOCATION), "a new particle awaits location"));
            PetscCall(PicurvAssertBool((PetscBool)(position[3 * p + 2] > user->Min_Z), "a new particle is inside the inlet"));
            /* The value was set on the face, before the particle moved in. */
            PetscCall(PicurvAssertRealNear(user->Min_Z * L + 1.0, psi[p], 1.0e-3, "a new particle carries the value at the face"));
            pid_sum += pid[p];
            position_sum += position[3 * p] + 2.0 * position[3 * p + 1] + 3.0 * position[3 * p + 2];
        }
        PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_location_status", NULL, NULL, (void **)&status));
        PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid));
        PetscCall(DMSwarmRestoreField(user->swarm, "Psi", NULL, NULL, (void **)&psi));
        PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void **)&position));
        PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, &pid_sum, 1, MPIU_INT64, MPI_SUM, PETSC_COMM_WORLD));
        PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, &position_sum, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD));
        PetscCall(PicurvAssertBool((PetscBool)(pid_sum == expected_pid_sum), "identities come from the grid position, the step and the count"));
        if (pass == 0) positions_first = position_sum; else positions_second = position_sum;
    }
    PetscCall(PicurvAssertRealNear(positions_first, positions_second, 1.0e-9, "injection carries no state from one call to the next"));
#undef TEST_INLET_FLUX

    {
        const PetscReal   diffusivity = 0.5;
        PetscInt          nlocal = 0;
        PetscReal        *without_kick = NULL;
        const PetscReal  *position = NULL;
        const PetscInt64 *pid = NULL;

        PetscCall(DMSwarmGetLocalSize(user->swarm, &nlocal));
        PetscCall(PetscMalloc1(3 * nlocal, &without_kick));
        PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void **)&position));
        PetscCall(PetscArraycpy(without_kick, position, 3 * nlocal));
        PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void **)&position));

        PetscCall(VecSet(user->lDiffusivity, diffusivity));
        PetscCall(DMSwarmSetLocalSizes(user->swarm, 0, -1));
        PetscCall(InjectParticlesAtInlets(user));
        PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void **)&position));
        PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid));
        for (PetscInt p = 0; p < nlocal; ++p) {
            Cmpnts kick;

            PetscCall(CalculateBrownianDisplacement(user, pid[p], diffusivity, &kick));
            PetscCall(PicurvAssertRealNear(without_kick[3 * p + 0] + kick.x, position[3 * p + 0], 1.0e-12, "a new particle takes the step's Brownian displacement in x"));
            PetscCall(PicurvAssertRealNear(without_kick[3 * p + 1] + kick.y, position[3 * p + 1], 1.0e-12, "a new particle takes the step's Brownian displacement in y"));
            PetscCall(PicurvAssertRealNear(without_kick[3 * p + 2] + kick.z, position[3 * p + 2], 1.0e-12, "a new particle takes the step's Brownian displacement in z"));
        }
        PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_pid", NULL, NULL, (void **)&pid));
        PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void **)&position));
        PetscCall(PetscFree(without_kick));
    }

    PetscCall(PicurvDestroyRuntimeContext(&simCtx));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscFunctionReturn(0);
}

/**
 * @brief Checks that IEM's linear reconstruction holds across subdomain seams.
 * @details Particles fill the box and carry Psi = 2 (x + y + z), so every cell's reconstruction
 *          reproduces the line up to the noise of its gradient. A seam that hid a neighbour's
 *          particles from the gradient would leave the cells beside it relaxing toward their
 *          stepped mean, as far from the line as the stepped mean leaves every cell. Each cell's
 *          particles are on its owner rank, so the stepped reference is computed there.
 */
static PetscErrorCode TestIEMReconstructionAcrossSeams(void)
{
    const char *options = "-numParticles 5000\n-pinit 1\n"
                          "-particle_fields_count 1\n"
                          "-particle_fields_0_name Psi\n"
                          "-particle_fields_0_expr_0 \"2*(x + y + z)\"\n";
    SimCtx     *simCtx = NULL;
    UserCtx    *user = NULL;
    char        tmpdir[PETSC_MAX_PATH_LEN];
    PetscInt    nlocal = 0, *cell = NULL, ncells, cx_n[3];
    PetscReal  *psi = NULL, *diffusivity = NULL, *pos = NULL, *cx = NULL, *cn = NULL;
    PetscReal   sums[4] = {0.0, 0.0, 0.0, 0.0}; /* sum before, sum after, squared departure from the line, the same for the stepped mean */
    PetscReal   decay;

    PetscFunctionBeginUser;
    PetscCall(PicurvBuildTinyRuntimeContextWithOptions(NULL, PETSC_TRUE, options, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(InitializeEulerianState(simCtx));
    PetscCall(InitializeParticleSwarm(simCtx));
    PetscCall(PerformInitializedParticleSetup(simCtx));
    simCtx->dt = 0.5;
    cx_n[0] = user->info.mx - 2; cx_n[1] = user->info.my - 2; cx_n[2] = user->info.mz - 2;
    ncells = cx_n[0] * cx_n[1] * cx_n[2];
    decay = PetscExpReal(-simCtx->iem_constant * 0.2 * 0.5 / PetscPowReal(1.0 / (PetscReal)ncells, 0.6666667));

    PetscCall(DMSwarmGetLocalSize(user->swarm, &nlocal));
    PetscCall(DMSwarmGetField(user->swarm, "Psi", NULL, NULL, (void **)&psi));
    PetscCall(DMSwarmGetField(user->swarm, "Diffusivity", NULL, NULL, (void **)&diffusivity));
    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void **)&pos));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell));
    PetscCall(PetscCalloc2(3 * ncells, &cx, ncells, &cn));
    for (PetscInt p = 0; p < nlocal; ++p) {
        const PetscInt c = cell[3*p] + cx_n[0] * (cell[3*p+1] + cx_n[1] * cell[3*p+2]);
        diffusivity[p] = 0.2;
        sums[0] += psi[p];
        for (PetscInt d = 0; d < 3; ++d) cx[3*c + d] += pos[3*p + d];
        cn[c] += 1.0;
    }
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell));
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void **)&pos));
    PetscCall(DMSwarmRestoreField(user->swarm, "Diffusivity", NULL, NULL, (void **)&diffusivity));
    PetscCall(DMSwarmRestoreField(user->swarm, "Psi", NULL, NULL, (void **)&psi));

    PetscCall(UpdateAllParticleFields(user));

    PetscCall(DMSwarmGetField(user->swarm, "Psi", NULL, NULL, (void **)&psi));
    PetscCall(DMSwarmGetField(user->swarm, "position", NULL, NULL, (void **)&pos));
    PetscCall(DMSwarmGetField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell));
    for (PetscInt p = 0; p < nlocal; ++p) {
        const PetscInt  c = cell[3*p] + cx_n[0] * (cell[3*p+1] + cx_n[1] * cell[3*p+2]);
        const PetscReal line = 2.0 * (pos[3*p] + pos[3*p+1] + pos[3*p+2]);
        const PetscReal mean = 2.0 * (cx[3*c] + cx[3*c+1] + cx[3*c+2]) / cn[c];
        const PetscReal stepped = mean + (line - mean) * decay;
        sums[1] += psi[p];
        sums[2] += (psi[p] - line) * (psi[p] - line);
        sums[3] += (stepped - line) * (stepped - line);
    }
    PetscCall(DMSwarmRestoreField(user->swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell));
    PetscCall(DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void **)&pos));
    PetscCall(DMSwarmRestoreField(user->swarm, "Psi", NULL, NULL, (void **)&psi));
    PetscCall(PetscFree2(cx, cn));
    PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, sums, 4, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD));

    PetscCall(PicurvAssertRealNear(sums[0], sums[1], 1.0e-9, "the reconstruction keeps every cell's scalar sum"));
    PetscCall(PicurvAssertBool((PetscBool)(sums[3] > 0.0 && sums[2] < 0.05 * sums[3]),
                               "relaxing toward the reconstruction must keep a linear field close to its line in every subdomain"));

    PetscCall(PicurvDestroyRuntimeContext(&simCtx));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscFunctionReturn(0);
}

/**
 * @brief Velocity sources draw by global index, whatever the decomposition.
 * @details Every value is compared with a draw this test makes from the cell's or face cell's
 *          global indices alone, so a key that depended on the rank layout would fail at any
 *          rank count but one. A cell reads one draw per component; a face cell a new one each
 *          step, and the same one when a step is revisited.
 */
static PetscErrorCode TestSourceDrawsAreKeyedOnTheGrid(void)
{
    const char *options = "-source_random_seed 77\n"
                          "-ic_source_expr_0 \"uniform()\"\n"
                          "-ic_source_expr_1 \"uniform()\"\n"
                          "-ic_source_expr_2 \"1 + uniform(3)\"\n"
                          "-face_source_0_expr_0 \"1 + 0.1*normal()\"\n";
    static const char *const kNames[] = {"t"};
    SimCtx           *simCtx = NULL;
    UserCtx          *user = NULL;
    char              tmpdir[PETSC_MAX_PATH_LEN];
    PicurvExpression *cell[3] = {NULL, NULL, NULL}, *face = NULL;
    PetscReal         worst = 0.0, changed = 0.0;
    PetscInt          compared[2] = {0, 0};

    PetscFunctionBeginUser;
    PetscCall(PicurvBuildTinyRuntimeContextWithOptions(NULL, PETSC_FALSE, options, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(InitializeEulerianState(simCtx));
    PetscCall(PicurvExpressionCompile("uniform()", kNames, 1, PETSC_TRUE, &cell[0]));
    PetscCall(PicurvExpressionCompile("uniform()", kNames, 1, PETSC_TRUE, &cell[1]));
    PetscCall(PicurvExpressionCompile("1 + uniform(3)", kNames, 1, PETSC_TRUE, &cell[2]));
    PetscCall(PicurvExpressionCompile("1 + 0.1*normal()", kNames, 1, PETSC_TRUE, &face));

    const DMDALocalInfo *info = &user->info;
    const PetscInt   nodes[3] = {info->mx, info->my, info->mz};
    const PetscReal  U_ref = simCtx->scaling.U_ref;
    const Cmpnts  ***ucont = NULL, ***csi = NULL, ***eta = NULL, ***zet = NULL;

    /* The initial field: each interior face is the draw of the cell behind it. */
    PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, (void *)&ucont));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCsi, (void *)&csi));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lEta, (void *)&eta));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lZet, (void *)&zet));
    for (PetscInt k = PetscMax(info->zs, 1); k < PetscMin(info->zs + info->zm, nodes[2] - 2); ++k)
    for (PetscInt j = PetscMax(info->ys, 1); j < PetscMin(info->ys + info->ym, nodes[1] - 2); ++j)
    for (PetscInt i = PetscMax(info->xs, 1); i < PetscMin(info->xs + info->xm, nodes[0] - 2); ++i) {
        PicurvExpressionDrawKey key = {77, i + (PetscInt64)nodes[0] * (j + (PetscInt64)nodes[1] * k),
                                       PICURV_DRAW_SALT_INITIAL_FIELD, 0, PETSC_TRUE};
        Cmpnts v;
        key.field = 0; PetscCall(PicurvExpressionEvaluate(cell[0], NULL, &key, &v.x));
        key.field = 1; PetscCall(PicurvExpressionEvaluate(cell[1], NULL, &key, &v.y));
        key.field = 2; PetscCall(PicurvExpressionEvaluate(cell[2], NULL, &key, &v.z));
        const Cmpnts m[3] = {csi[k][j][i], eta[k][j][i], zet[k][j][i]};
        const PetscReal found[3] = {ucont[k][j][i].x, ucont[k][j][i].y, ucont[k][j][i].z};
        for (PetscInt a = 0; a < 3; ++a) {
            const PetscReal expected = (v.x * m[a].x + v.y * m[a].y + v.z * m[a].z) / U_ref;
            worst = PetscMax(worst, PetscAbsReal(found[a] - expected));
        }
        compared[0]++;
    }
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, (void *)&ucont));

    /* The inlet: a draw per face cell and step. */
    PetscReal *first = NULL;
    const PetscInt owned[2] = {PetscMin(info->xs + info->xm, nodes[0] - 1) - PetscMax(info->xs, 1),
                               PetscMin(info->ys + info->ym, nodes[1] - 1) - PetscMax(info->ys, 1)};
    PetscCall(PetscCalloc1(owned[0] * owned[1], &first));
    for (PetscInt pass = 0; pass < 3; ++pass) {
        const PetscInt64 step = (pass == 1) ? 4 : 3;
        simCtx->step = step;
        PetscCall(BoundarySystem_ExecuteStep(user));
        if (info->zs != 0) continue;
        PetscCall(DMDAVecGetArrayRead(user->fda, user->Ucont, (void *)&ucont));
        for (PetscInt j = 1, n = 0; j < nodes[1] - 1; ++j)
        for (PetscInt i = 1; i < nodes[0] - 1; ++i) {
            if (j < info->ys || j >= info->ys + info->ym || i < info->xs || i >= info->xs + info->xm) continue;
            const PicurvExpressionDrawKey key = {77, PicurvFaceCellId(nodes, BC_FACE_NEG_Z, i - 1, j - 1), step, 0, PETSC_TRUE};
            const Cmpnts m = zet[0][j][i];
            PetscReal speed;
            PetscCall(PicurvExpressionEvaluate(face, NULL, &key, &speed));
            worst = PetscMax(worst, PetscAbsReal(ucont[0][j][i].z - speed / U_ref * PetscSqrtReal(m.x * m.x + m.y * m.y + m.z * m.z)));
            if (pass == 0) first[n] = ucont[0][j][i].z;
            if (pass == 1) changed += PetscAbsReal(ucont[0][j][i].z - first[n]);
            if (pass == 2) worst = PetscMax(worst, PetscAbsReal(ucont[0][j][i].z - first[n]));
            compared[1]++;
            ++n;
        }
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Ucont, (void *)&ucont));
    }
    PetscCall(PetscFree(first));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lZet, (void *)&zet));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lEta, (void *)&eta));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCsi, (void *)&csi));
    PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, &worst, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, &changed, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, compared, 2, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertIntEqual((nodes[0] - 3) * (nodes[1] - 3) * (nodes[2] - 3), compared[0], "every interior cell was compared"));
    PetscCall(PicurvAssertIntEqual(3 * (nodes[0] - 2) * (nodes[1] - 2), compared[1], "every inlet cell was compared on each pass"));
    PetscCall(PicurvAssertRealNear(0.0, worst, 1.0e-12, "a value is the draw its global index keys"));
    PetscCall(PicurvAssertBool((PetscBool)(changed > 0.0), "a later step draws anew"));
    for (PetscInt a = 0; a < 3; ++a) PetscCall(PicurvExpressionDestroy(&cell[a]));
    PetscCall(PicurvExpressionDestroy(&face));
    PetscCall(PicurvDestroyRuntimeContext(&simCtx));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscFunctionReturn(0);
}

/**
 * @brief A filtered inflow draw at face cell (g0, g1) of -Zeta: the separable kernel applied to
 *        the keyed normals on a lattice padded by the half-widths.
 */
static PetscReal InflowOracle(PetscInt64 seed, PetscInt c, PetscInt g0, PetscInt g1, PetscInt64 step,
                              const PetscInt cells[2], const PetscInt half[2], PetscReal *const kernel[2])
{
    PetscReal sum = 0.0;
    for (PetscInt k1 = -half[1]; k1 <= half[1]; ++k1)
    for (PetscInt k0 = -half[0]; k0 <= half[0]; ++k0) {
        const PetscInt64 id = ((PetscInt64)BC_FACE_NEG_Z * (cells[1] + 2 * half[1]) + (g1 + k1 + half[1])) * (cells[0] + 2 * half[0])
                              + (g0 + k0 + half[0]);
        sum += kernel[0][k0 + half[0]] * kernel[1][k1 + half[1]]
               * PicurvNormalFromHash(PicurvDrawHash(seed, PICURV_DRAW_SPACE_INFLOW, (uint64_t)c, id, step));
    }
    return sum;
}

/**
 * @brief Tests that an inlet's synthetic fluctuation is the digital filter of draws keyed on the
 *        grid, scaled by the Cholesky factor of the Reynolds stresses, on any decomposition.
 * @details The oracle rebuilds the kernels from the grid spacing, filters the keyed normals at
 *          each owned inlet cell, and checks the state in the dummy plane and the boundary
 *          velocity at the first step and after one time-filter update. Lengths are chosen so
 *          the halo passes the face's closed edges.
 */
static PetscErrorCode TestInletFluctuationIsTheFilteredKeyedDraw(void)
{
    const char *options = "-source_random_seed 91\n"
                          "-inflow_fluctuations true\n"
                          "-face_source_0_fluct true\n"
                          "-face_source_0_fluct_stress_0_frame cartesian\n"
                          "-face_source_0_fluct_stress_0_expr_0 0.04\n"
                          "-face_source_0_fluct_stress_0_expr_1 0.01\n"
                          "-face_source_0_fluct_stress_0_expr_2 0.09\n"
                          "-face_source_0_fluct_stress_1_frame cartesian\n"
                          "-face_source_0_fluct_stress_1_expr_0 0.01\n"
                          "-face_source_0_fluct_stress_1_expr_1 0.0\n"
                          "-face_source_0_fluct_stress_1_expr_2 0.005\n"
                          "-face_source_0_fluct_length 0.25,0.13\n"
                          "-face_source_0_fluct_time 0.05\n";
    const PetscReal length[2] = {0.25, 0.13}, T = 0.05;
    SimCtx         *simCtx = NULL;
    UserCtx        *user = NULL;
    char            tmpdir[PETSC_MAX_PATH_LEN];
    PetscReal       worst = 0.0, moved = 0.0;
    PetscInt        compared = 0;
    PetscReal      *kernel[2] = {NULL, NULL};
    PetscInt        half[2];

    PetscFunctionBeginUser;
    PetscCall(PicurvBuildTinyRuntimeContextWithOptions(NULL, PETSC_FALSE, options, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(InitializeEulerianState(simCtx));

    const DMDALocalInfo *info = &user->info;
    const PetscInt   nodes[3] = {info->mx, info->my, info->mz};
    const PetscInt   cells[2] = {nodes[0] - 2, nodes[1] - 2};
    const PetscReal  U_ref = simCtx->scaling.U_ref, U2 = U_ref * U_ref;
    const PetscReal  decay = PetscExpReal(-PETSC_PI * simCtx->dt / (2.0 * T)), gain = PetscSqrtReal(1.0 - decay * decay);
    /* A = chol(R / U_ref^2) for R = [[0.04, 0.01, 0], [0.01, 0.01, 0.005], [0, 0.005, 0.09]]. */
    PetscReal        A[3][3] = {{0.0}};
    A[0][0] = 0.2 / U_ref;
    A[1][0] = 0.01 / U2 / A[0][0];
    A[1][1] = PetscSqrtReal(0.01 / U2 - A[1][0] * A[1][0]);
    A[2][1] = 0.005 / U2 / A[1][1];
    A[2][2] = PetscSqrtReal(0.09 / U2 - A[2][1] * A[2][1]);

    /* The unit box: h = 1 / cells, n = L / h, weights exp(-pi |k| / n) to a unit sum of squares. */
    for (PetscInt d = 0; d < 2; ++d) {
        const PetscReal n = length[d] * (PetscReal)cells[d];
        PetscReal       norm = 0.0;
        half[d] = PetscMax(1, (PetscInt)PetscCeilReal(2.0 * n));
        PetscCall(PetscMalloc1(2 * half[d] + 1, &kernel[d]));
        for (PetscInt k = -half[d]; k <= half[d]; ++k) {
            kernel[d][k + half[d]] = PetscExpReal(-PETSC_PI * PetscAbsReal((PetscReal)k) / n);
            norm += kernel[d][k + half[d]] * kernel[d][k + half[d]];
        }
        for (PetscInt k = 0; k <= 2 * half[d]; ++k) kernel[d][k] /= PetscSqrtReal(norm);
    }
    PetscCall(PicurvAssertBool((PetscBool)(half[0] > cells[0] / 2), "the halo passes a closed edge of the face"));

    PetscReal *previous = NULL;
    PetscCall(PetscCalloc1(3 * info->xm * info->ym, &previous));
    for (PetscInt pass = 0; pass < 2; ++pass) {
        if (pass == 1) {
            simCtx->step = 1;
            PetscCall(BoundarySystem_ExecuteStep(user));
        }
        if (info->zs != 0) continue;
        const Cmpnts ***state = NULL, ***ubcs = NULL;
        PetscCall(DMDAVecGetArrayRead(user->fda, user->InflowFluctuation, (void *)&state));
        PetscCall(DMDAVecGetArrayRead(user->fda, user->Bcs.Ubcs, (void *)&ubcs));
        for (PetscInt j = PetscMax(info->ys, 1); j < PetscMin(info->ys + info->ym, nodes[1] - 1); ++j)
        for (PetscInt i = PetscMax(info->xs, 1); i < PetscMin(info->xs + info->xm, nodes[0] - 1); ++i) {
            PetscReal *psi = &previous[3 * ((i - info->xs) + info->xm * (j - info->ys))];
            for (PetscInt c = 0; c < 3; ++c) {
                const PetscReal draw = InflowOracle(91, c, i - 1, j - 1, pass, cells, half, kernel);
                psi[c] = (pass == 0) ? draw : decay * psi[c] + gain * draw;
            }
            const PetscReal found[3] = {state[0][j][i].x, state[0][j][i].y, state[0][j][i].z};
            const PetscReal u[3] = {ubcs[0][j][i].x, ubcs[0][j][i].y, ubcs[0][j][i].z - 1.5 / U_ref};
            for (PetscInt c = 0; c < 3; ++c) {
                const PetscReal expected = A[c][0] * psi[0] + A[c][1] * psi[1] + A[c][2] * psi[2];
                worst = PetscMax(worst, PetscAbsReal(found[c] - psi[c]));
                worst = PetscMax(worst, PetscAbsReal(u[c] - expected));
                if (pass == 1) moved += PetscAbsReal(u[c]);
            }
            compared++;
        }
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->Bcs.Ubcs, (void *)&ubcs));
        PetscCall(DMDAVecRestoreArrayRead(user->fda, user->InflowFluctuation, (void *)&state));
    }
    PetscCall(PetscFree(previous));
    PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, &worst, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, &moved, 1, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, &compared, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertIntEqual(2 * cells[0] * cells[1], compared, "every inlet cell was compared on each pass"));
    PetscCall(PicurvAssertRealNear(0.0, worst, 1.0e-12, "the state and inlet velocity are the filtered keyed draw"));
    PetscCall(PicurvAssertBool((PetscBool)(moved > 0.0), "the inlet carries a fluctuation"));
    for (PetscInt d = 0; d < 2; ++d) PetscCall(PetscFree(kernel[d]));
    PetscCall(PicurvDestroyRuntimeContext(&simCtx));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that a face that is not an INLET refuses to inject.
 */
static PetscErrorCode TestInjectionRefusesAFaceThatIsNotAnInlet(void)
{
    const char    *options = "-particle_inlet_count 1\n-particle_inlet_0_face 0\n-particle_inlet_0_particles_per_cell \"1.0\"\n";
    SimCtx        *simCtx = NULL;
    UserCtx       *user = NULL;
    char           tmpdir[PETSC_MAX_PATH_LEN];
    PetscErrorCode refused;

    PetscFunctionBeginUser;
    PetscCall(PicurvBuildTinyRuntimeContextWithOptions(NULL, PETSC_FALSE, options, &simCtx, &user, tmpdir, sizeof(tmpdir)));
    PetscCall(InitializeEulerianState(simCtx));
    PetscCall(InitializeParticleSwarm(simCtx));
    simCtx->step = 1;
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    refused = InjectParticlesAtInlets(user);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertBool((PetscBool)(refused != 0), "a wall does not inject particles"));
    PetscCall(PicurvDestroyRuntimeContext(&simCtx));
    PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that, across a periodic seam and rank boundaries, the face planes the search and
 *        interpolation read after a ghost update are each cell's own planes.
 * @details Every grid cell made of physical nodes and every dual cell the interpolation can form
 *          must give, to rounding, the distances CalculateDistancesToCellFaces() computes from that
 *          cell's own vertices, whatever the fields' sync classes do to the ghosts beyond.
 */
static PetscErrorCode TestFacePlanesSurvivePeriodicGhostUpdate(void)
{
    SimCtx         *simCtx = NULL;
    UserCtx        *user = NULL;
    FacePlaneArrays arrays;
    CellFacePlanes  shared;
    Cell            cell;
    PetscReal       d_shared[NUM_FACES], d_cell[NUM_FACES];
    const Cmpnts ***cent = NULL;
    const FieldId   fields[] = {FIELD_ID_PLANE_NORMAL_I, FIELD_ID_PLANE_NORMAL_J, FIELD_ID_PLANE_NORMAL_K,
                                FIELD_ID_PLANE_CENTROID_I, FIELD_ID_PLANE_CENTROID_J, FIELD_ID_PLANE_CENTROID_K,
                                FIELD_ID_SEARCH_CELL_SIZE, FIELD_ID_DUAL_PLANE_NORMAL_I, FIELD_ID_DUAL_PLANE_NORMAL_J,
                                FIELD_ID_DUAL_PLANE_NORMAL_K, FIELD_ID_DUAL_PLANE_CENTROID_I,
                                FIELD_ID_DUAL_PLANE_CENTROID_J, FIELD_ID_DUAL_PLANE_CENTROID_K};
    PetscInt        checked = 0, checked_global = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user, 8, 6, 6, PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
    MarkXPeriodic(user);
    for (size_t f = 0; f < sizeof(fields) / sizeof(fields[0]); ++f) PetscCall(UpdateLocalGhosts(user, fields[f]));

    const DMDALocalInfo *info = &user->info;
    const PetscInt gs[3] = {info->gxs, info->gys, info->gzs};
    const PetscInt ge[3] = {info->gxs + info->gxm, info->gys + info->gym, info->gzs + info->gzm};
    const PetscInt mx[3] = {info->mx, info->my, info->mz};
    const Cmpnts   p = {0.43, 0.38, 0.61};

    /* Grid cells: origins whose nodes are physical (0..mx-2) and local. */
    PetscCall(MetricGetFacePlaneArrays(user, PETSC_FALSE, &arrays));
    for (PetscInt k = PetscMax(gs[2], 0); k + 1 <= PetscMin(ge[2] - 1, mx[2] - 2); ++k)
        for (PetscInt j = PetscMax(gs[1], 0); j + 1 <= PetscMin(ge[1] - 1, mx[1] - 2); ++j)
            for (PetscInt i = PetscMax(gs[0], 0); i + 1 <= PetscMin(ge[0] - 1, mx[0] - 2); ++i) {
                MetricGatherCellFacePlanes(&arrays, info, i, j, k, &shared);
                PetscCall(RetrieveCurrentCell(user, i, j, k, &cell));
                PetscCall(CalculateDistancesToCellFaces(p, &cell, d_cell, 0.0, NULL));
                PetscCall(CalculateDistancesToCellFacePlanes(p, &shared, d_shared, 0.0, NULL));
                for (PetscInt f = 0; f < NUM_FACES; ++f)
                    PetscCall(PicurvAssertRealNear(d_cell[f], d_shared[f], 1.0e-14, "grid cell plane after ghost update"));
                checked++;
            }
    PetscCall(MetricRestoreFacePlaneArrays(user, &arrays));

    /* Dual cells: centres 0..mx-1 on the periodic x axis, 1..mx-2 across the walls. */
    PetscCall(MetricGetFacePlaneArrays(user, PETSC_TRUE, &arrays));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCent, (void *)&cent));
    for (PetscInt k = PetscMax(gs[2], 1); k + 1 <= PetscMin(ge[2] - 1, mx[2] - 2); ++k)
        for (PetscInt j = PetscMax(gs[1], 1); j + 1 <= PetscMin(ge[1] - 1, mx[1] - 2); ++j)
            for (PetscInt i = PetscMax(gs[0], 0); i + 1 <= PetscMin(ge[0] - 1, mx[0] - 1); ++i) {
                cell.vertices[0] = cent[k][j][i];         cell.vertices[1] = cent[k][j][i + 1];
                cell.vertices[2] = cent[k][j + 1][i + 1]; cell.vertices[3] = cent[k][j + 1][i];
                cell.vertices[4] = cent[k + 1][j + 1][i]; cell.vertices[5] = cent[k + 1][j + 1][i + 1];
                cell.vertices[6] = cent[k + 1][j][i + 1]; cell.vertices[7] = cent[k + 1][j][i];
                MetricGatherCellFacePlanes(&arrays, info, i, j, k, &shared);
                PetscCall(CalculateDistancesToCellFaces(p, &cell, d_cell, 0.0, NULL));
                PetscCall(CalculateDistancesToCellFacePlanes(p, &shared, d_shared, 0.0, NULL));
                for (PetscInt f = 0; f < NUM_FACES; ++f)
                    PetscCall(PicurvAssertRealNear(d_cell[f], d_shared[f], 1.0e-14, "dual cell plane after ghost update"));
                checked++;
            }
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCent, (void *)&cent));
    PetscCall(MetricRestoreFacePlaneArrays(user, &arrays));

    PetscCallMPI(MPI_Allreduce(&checked, &checked_global, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD));
    PetscCall(PicurvAssertBool((PetscBool)(checked_global > 0), "some cells were checked"));
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
        {"face-planes-survive-periodic-ghost-update", TestFacePlanesSurvivePeriodicGhostUpdate},
        {"spatial-target-decomposition-independent", TestSpatialTargetIsDecompositionIndependent},
        {"statistics-checkpoint-round-trip-multi-rank", TestStatisticsCheckpointRoundTripMultiRank},
        {"les-averaging-decomposition-independent", TestLESAveragingIsDecompositionIndependent},
        {"find-owner-of-cell-matches-rank-ranges", TestFindOwnerOfCellMatchesRankRanges},
        {"profiling-summary-reduces-over-ranks", TestProfilingSummaryReducesOverRanks},
        {"restart-cellid-migration-moves-particle-to-owner", TestRestartCellIdMigrationMovesParticleToOwner},
        {"settlement-classifies-exit-once-across-passes", TestSettlementClassifiesExitOnceAcrossPasses},
        {"grid-line-fractions-follow-arc-length", TestGridLineFractionsFollowArcLength},
        {"injection-follows-flux-and-is-keyed-on-the-grid", TestInjectionFollowsFluxAndIsKeyedOnTheGrid},
        {"injection-refuses-a-face-that-is-not-an-inlet", TestInjectionRefusesAFaceThatIsNotAnInlet},
        {"source-draws-are-keyed-on-the-grid", TestSourceDrawsAreKeyedOnTheGrid},
        {"inlet-fluctuation-is-the-filtered-keyed-draw", TestInletFluctuationIsTheFilteredKeyedDraw},
        {"iem-reconstruction-across-seams", TestIEMReconstructionAcrossSeams},
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
