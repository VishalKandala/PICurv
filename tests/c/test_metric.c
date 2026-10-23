/**
 * @file test_metric.c
 * @brief C unit tests for metric tensors, face geometry, and length-scale helpers.
 */

#include "test_support.h"

#include "Metric.h"
/**
 * @brief Tests inversion of a diagonal covariant metric tensor.
 */

static PetscErrorCode TestInvertCovariantMetricTensorDiagonal(void)
{
    double covariant[3][3] = {
        {4.0, 0.0, 0.0},
        {0.0, 9.0, 0.0},
        {0.0, 0.0, 16.0},
    };
    double contravariant[3][3] = {{0.0}};

    PetscFunctionBeginUser;
    PetscCall(InvertCovariantMetricTensor(covariant, contravariant));
    PetscCall(PicurvAssertRealNear(0.25, contravariant[0][0], 1.0e-12, "inverse(0,0)"));
    PetscCall(PicurvAssertRealNear(1.0 / 9.0, contravariant[1][1], 1.0e-12, "inverse(1,1)"));
    PetscCall(PicurvAssertRealNear(0.0625, contravariant[2][2], 1.0e-12, "inverse(2,2)"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests contravariant velocity computation for identity metrics.
 */

static PetscErrorCode TestMetricVelocityContravariantIdentity(void)
{
    const PetscReal J[3][3] = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
    };
    const PetscReal u[3] = {1.5, -2.0, 4.0};
    PetscReal uc[3] = {0.0, 0.0, 0.0};

    PetscFunctionBeginUser;
    PetscCall(MetricVelocityContravariant(J, 1.0, u, uc));
    PetscCall(PicurvAssertRealNear(1.5, uc[0], 1.0e-12, "identity map uc[0]"));
    PetscCall(PicurvAssertRealNear(-2.0, uc[1], 1.0e-12, "identity map uc[1]"));
    PetscCall(PicurvAssertRealNear(4.0, uc[2], 1.0e-12, "identity map uc[2]"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests contravariant velocity computation for scaled axes.
 */

static PetscErrorCode TestMetricVelocityContravariantScaledAxes(void)
{
    const PetscReal J[3][3] = {
        {2.0, 0.0, 0.0},
        {0.0, 3.0, 0.0},
        {0.0, 0.0, 4.0},
    };
    const PetscReal u[3] = {2.0, 3.0, 4.0};
    PetscReal uc[3] = {0.0, 0.0, 0.0};

    PetscFunctionBeginUser;
    PetscCall(MetricVelocityContravariant(J, 24.0, u, uc));
    PetscCall(PicurvAssertRealNear(1.0, uc[0], 1.0e-12, "scaled map uc[0]"));
    PetscCall(PicurvAssertRealNear(1.0, uc[1], 1.0e-12, "scaled map uc[1]"));
    PetscCall(PicurvAssertRealNear(1.0, uc[2], 1.0e-12, "scaled map uc[2]"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests face-normal and face-area calculation on an axis-aligned cell.
 */

static PetscErrorCode TestCalculateFaceNormalAndAreaAxisAligned(void)
{
    const Cmpnts csi = {2.0, 0.0, 0.0};
    const Cmpnts eta = {0.0, 3.0, 0.0};
    const Cmpnts zet = {0.0, 0.0, 4.0};
    double ni[3] = {0.0, 0.0, 0.0};
    double nj[3] = {0.0, 0.0, 0.0};
    double nk[3] = {0.0, 0.0, 0.0};
    double Ai = 0.0, Aj = 0.0, Ak = 0.0;

    PetscFunctionBeginUser;
    PetscCall(CalculateFaceNormalAndArea(csi, eta, zet, ni, nj, nk, &Ai, &Aj, &Ak));
    PetscCall(PicurvAssertRealNear(2.0, Ai, 1.0e-12, "Ai"));
    PetscCall(PicurvAssertRealNear(3.0, Aj, 1.0e-12, "Aj"));
    PetscCall(PicurvAssertRealNear(4.0, Ak, 1.0e-12, "Ak"));
    PetscCall(PicurvAssertRealNear(1.0, ni[0], 1.0e-12, "ni.x"));
    PetscCall(PicurvAssertRealNear(0.0, ni[1], 1.0e-12, "ni.y"));
    PetscCall(PicurvAssertRealNear(0.0, ni[2], 1.0e-12, "ni.z"));
    PetscCall(PicurvAssertRealNear(0.0, nj[0], 1.0e-12, "nj.x"));
    PetscCall(PicurvAssertRealNear(1.0, nj[1], 1.0e-12, "nj.y"));
    PetscCall(PicurvAssertRealNear(0.0, nj[2], 1.0e-12, "nj.z"));
    PetscCall(PicurvAssertRealNear(0.0, nk[0], 1.0e-12, "nk.x"));
    PetscCall(PicurvAssertRealNear(0.0, nk[1], 1.0e-12, "nk.y"));
    PetscCall(PicurvAssertRealNear(1.0, nk[2], 1.0e-12, "nk.z"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests characteristic length-scale computation on an axis-aligned cell.
 */

static PetscErrorCode TestComputeCellCharacteristicLengthScaleAxisAligned(void)
{
    const Cmpnts csi = {2.0, 0.0, 0.0};
    const Cmpnts eta = {0.0, 3.0, 0.0};
    const Cmpnts zet = {0.0, 0.0, 4.0};
    double dx = 0.0, dy = 0.0, dz = 0.0;

    PetscFunctionBeginUser;
    PetscCall(ComputeCellCharacteristicLengthScale(1.0, csi, eta, zet, &dx, &dy, &dz));
    PetscCall(PicurvAssertRealNear(0.5, dx, 1.0e-12, "dx"));
    PetscCall(PicurvAssertRealNear(1.0 / 3.0, dy, 1.0e-12, "dy"));
    PetscCall(PicurvAssertRealNear(0.25, dz, 1.0e-12, "dz"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests the length scales of a wall-resolved cell, whose metrics are tiny.
 *
 * @details The determinant of the face-area vectors is the square of the cell volume for
 *          an orthogonal cell, so it falls as the sixth power of the spacing. The cell
 *          here is one from the Humphrey wall-resolved grid: 4.73e-4 x 4.73e-4 x 3.05e-2,
 *          whose determinant is about 4.7e-17. An absolute singularity floor rejected it
 *          and left the length scales uninitialized, which reached the solver as a NaN
 *          eddy viscosity rather than as an error.
 */
static PetscErrorCode TestComputeCellCharacteristicLengthScaleWallResolvedCell(void)
{
    const double wall_normal = 4.73e-4, streamwise = 3.049381e-2;
    const Cmpnts csi = {wall_normal * streamwise, 0.0, 0.0};
    const Cmpnts eta = {0.0, wall_normal * streamwise, 0.0};
    const Cmpnts zet = {0.0, 0.0, wall_normal * wall_normal};
    const double volume = wall_normal * wall_normal * streamwise;
    double dx = 0.0, dy = 0.0, dz = 0.0;

    PetscFunctionBeginUser;
    PetscCall(ComputeCellCharacteristicLengthScale(1.0 / volume, csi, eta, zet, &dx, &dy, &dz));
    PetscCall(PicurvAssertRealNear(wall_normal, dx, 1.0e-15, "dx of a wall-resolved cell"));
    PetscCall(PicurvAssertRealNear(wall_normal, dy, 1.0e-15, "dy of a wall-resolved cell"));
    PetscCall(PicurvAssertRealNear(streamwise, dz, 1.0e-15, "dz of a wall-resolved cell"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that MetricJacobian returns the physical Jacobian of a sheared, scaled cell, whose determinant is its volume.
 */
static PetscErrorCode TestMetricJacobianIsCellVolumeForParallelepiped(void)
{
    const PetscReal a = 0.5, b = 2.0, c = 3.0, shear = 0.7;
    Cmpnts          nodes[2][2][2], *rows[2][2], **planes[2];
    PetscReal       J[3][3], detJ;

    PetscFunctionBeginUser;
    for (PetscInt k = 0; k < 2; ++k) {
        for (PetscInt j = 0; j < 2; ++j) {
            for (PetscInt i = 0; i < 2; ++i) {
                nodes[k][j][i].x = a * i + shear * k;
                nodes[k][j][i].y = b * j;
                nodes[k][j][i].z = c * k;
            }
            rows[k][j] = nodes[k][j];
        }
        planes[k] = rows[k];
    }
    PetscCall(MetricJacobian(NULL, (const Cmpnts ***)planes, 0, 0, 0, 0.3, 0.6, 0.9, J, &detJ));
    PetscCall(PicurvAssertRealNear(a * b * c, detJ, 1.0e-14, "detJ of a parallelepiped is its volume"));
    PetscCall(PicurvAssertRealNear(a, J[0][0], 1.0e-14, "dx/dxi is the i edge length"));
    PetscCall(PicurvAssertRealNear(shear, J[0][2], 1.0e-14, "dx/dzta is the shear of the k edge"));
    PetscCall(PicurvAssertRealNear(c, J[2][2], 1.0e-14, "dz/dzta is the k edge length"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that a degenerate cell is still refused, and that the refusal propagates.
 */

static PetscErrorCode TestDegenerateCellIsRefused(void)
{
    const Cmpnts csi = {2.0, 0.0, 0.0};
    const Cmpnts eta = {4.0, 0.0, 0.0};   /* parallel to csi: the rows are coplanar */
    const Cmpnts zet = {0.0, 0.0, 3.0};
    double dx = 0.0, dy = 0.0, dz = 0.0;
    PetscErrorCode refused = 0;

    PetscFunctionBeginUser;
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    refused = ComputeCellCharacteristicLengthScale(1.0, csi, eta, zet, &dx, &dy, &dz);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_MAT_LU_ZRPVT, refused,
                                   "a coplanar cell is refused, and the caller is told"));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that the shared per-face planes reproduce each cell's own planes and are observable.
 * @details On the fixture's uniform grid every face is planar, so a grid cell's gathered planes
 *          and a dual cell's must give the per-cell distances CalculateDistancesToCellFaces()
 *          computes from the same vertices, to rounding; the catalog's diagnostics must accept them.
 */
static PetscErrorCode TestFacePlanesMatchPerCellPlanes(void)
{
    SimCtx         *simCtx = NULL;
    UserCtx        *user = NULL;
    FacePlaneArrays arrays;
    CellFacePlanes  shared;
    Cell            cell;
    PetscReal       d_shared[NUM_FACES], d_cell[NUM_FACES], size;
    const Cmpnts    p = {0.31, 0.47, 0.52};
    const Cmpnts ***cent = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 6, 6, 6));
    PetscCall(MetricComputeFacePlanes(user));

    PetscCall(MetricGetFacePlaneArrays(user, PETSC_FALSE, &arrays));
    for (PetscInt c = 0; c < 3; ++c) {
        MetricGatherCellFacePlanes(&arrays, &user->info, c, c + 1, c, &shared);
        PetscCall(RetrieveCurrentCell(user, c, c + 1, c, &cell));
        PetscCall(CalculateDistancesToCellFaces(p, &cell, d_cell, 0.0, NULL));
        PetscCall(CalculateDistancesToCellFacePlanes(p, &shared, d_shared, 0.0, NULL));
        PetscCall(GetCellCharacteristicSize(&cell, &size));
        PetscCall(PicurvAssertIntEqual(0, shared.degenerate, "grid cell faces all have planes"));
        PetscCall(PicurvAssertRealNear(size, shared.size, 1.0e-15, "grid cell size"));
        for (PetscInt f = 0; f < NUM_FACES; ++f) PetscCall(PicurvAssertRealNear(d_cell[f], d_shared[f], 1.0e-14, "grid cell face distance"));
    }
    PetscCall(MetricRestoreFacePlaneArrays(user, &arrays));

    PetscCall(MetricGetFacePlaneArrays(user, PETSC_TRUE, &arrays));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCent, (void *)&cent));
    for (PetscInt c = 1; c < 4; ++c) {
        const PetscInt i = c, j = 2, k = 4 - c;
        cell.vertices[0] = cent[k][j][i];         cell.vertices[1] = cent[k][j][i + 1];
        cell.vertices[2] = cent[k][j + 1][i + 1]; cell.vertices[3] = cent[k][j + 1][i];
        cell.vertices[4] = cent[k + 1][j + 1][i]; cell.vertices[5] = cent[k + 1][j + 1][i + 1];
        cell.vertices[6] = cent[k + 1][j][i + 1]; cell.vertices[7] = cent[k + 1][j][i];
        MetricGatherCellFacePlanes(&arrays, &user->info, i, j, k, &shared);
        PetscCall(CalculateDistancesToCellFaces(p, &cell, d_cell, 0.0, NULL));
        PetscCall(CalculateDistancesToCellFacePlanes(p, &shared, d_shared, 0.0, NULL));
        for (PetscInt f = 0; f < NUM_FACES; ++f) PetscCall(PicurvAssertRealNear(d_cell[f], d_shared[f], 1.0e-14, "dual cell face distance"));
    }
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCent, (void *)&cent));
    PetscCall(MetricRestoreFacePlaneArrays(user, &arrays));

    PetscCall(LOG_FIELD_MIN_MAX(user, FIELD_ID_PLANE_NORMAL_I));
    PetscCall(LOG_FIELD_ANATOMY(user, FIELD_ID_DUAL_PLANE_CENTROID_K, "face-plane test"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Runs the unit-metric PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"invert-covariant-metric-diagonal", TestInvertCovariantMetricTensorDiagonal},
        {"metric-velocity-contravariant-identity", TestMetricVelocityContravariantIdentity},
        {"metric-velocity-contravariant-scaled-axes", TestMetricVelocityContravariantScaledAxes},
        {"face-normal-and-area-axis-aligned", TestCalculateFaceNormalAndAreaAxisAligned},
        {"characteristic-length-scale-axis-aligned", TestComputeCellCharacteristicLengthScaleAxisAligned},
        {"characteristic-length-scale-wall-resolved-cell", TestComputeCellCharacteristicLengthScaleWallResolvedCell},
        {"degenerate-cell-is-refused", TestDegenerateCellIsRefused},
        {"metric-jacobian-is-cell-volume-for-parallelepiped", TestMetricJacobianIsCellVolumeForParallelepiped},
        {"face-planes-match-per-cell-planes", TestFacePlanesMatchPerCellPlanes},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv metric tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-metric", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
