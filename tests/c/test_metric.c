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
