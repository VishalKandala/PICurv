/**
 * @file test_geometry.c
 * @brief C unit tests for interpolation and signed-distance geometry helpers.
 */

#include "test_support.h"

#include "interpolation.h"
#include "walkingsearch.h"
/**
 * @brief Tests scalar interpolation on a constant field.
 */

static PetscErrorCode TestScalarInterpolationConstantField(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Vec input = NULL, output = NULL;
    PetscReal ***in_arr = NULL, ***out_arr = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(VecDuplicate(user->P, &input));
    PetscCall(VecDuplicate(user->P, &output));
    PetscCall(VecSet(input, 8.0));
    PetscCall(VecSet(output, 0.0));

    PetscCall(DMDAVecGetArrayRead(user->da, input, &in_arr));
    PetscCall(DMDAVecGetArray(user->da, output, &out_arr));
    PetscCall(InterpolateFieldFromCornerToCenter_Scalar(in_arr, out_arr, user));
    PetscCall(PicurvAssertRealNear(8.0, out_arr[1][1][1], 1.0e-12, "scalar corner->center preserves constants"));
    PetscCall(DMDAVecRestoreArrayRead(user->da, input, &in_arr));
    PetscCall(DMDAVecRestoreArray(user->da, output, &out_arr));

    PetscCall(VecDestroy(&input));
    PetscCall(VecDestroy(&output));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests vector interpolation on a constant field.
 */

static PetscErrorCode TestVectorInterpolationConstantField(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    Vec input = NULL, output = NULL;
    Cmpnts ***in_arr = NULL, ***out_arr = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(VecDuplicate(user->Ucat, &input));
    PetscCall(VecDuplicate(user->Ucat, &output));
    PetscCall(VecSet(input, 3.5));
    PetscCall(VecSet(output, 0.0));

    PetscCall(DMDAVecGetArrayRead(user->fda, input, &in_arr));
    PetscCall(DMDAVecGetArray(user->fda, output, &out_arr));
    PetscCall(InterpolateFieldFromCornerToCenter_Vector(in_arr, out_arr, user));
    PetscCall(PicurvAssertRealNear(3.5, out_arr[1][1][1].x, 1.0e-12, "vector corner->center x preserves constants"));
    PetscCall(PicurvAssertRealNear(3.5, out_arr[1][1][1].y, 1.0e-12, "vector corner->center y preserves constants"));
    PetscCall(PicurvAssertRealNear(3.5, out_arr[1][1][1].z, 1.0e-12, "vector corner->center z preserves constants"));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, input, &in_arr));
    PetscCall(DMDAVecRestoreArray(user->fda, output, &out_arr));

    PetscCall(VecDestroy(&input));
    PetscCall(VecDestroy(&output));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests signed-distance computation and inside/outside classification.
 */

static PetscErrorCode TestSignedDistanceAndClassification(void)
{
    /* The unit cube in walking-search vertex order; its z = 0 face is BACK. */
    const PetscReal corner[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                    {0, 1, 1}, {1, 1, 1}, {1, 0, 1}, {0, 0, 1}};
    const Cmpnts inside_point = {0.5, 0.5, 0.25};
    const Cmpnts outside_point = {0.5, 0.5, -1.0};
    Cell           cell;
    CellFacePlanes planes;
    PetscReal d_inside[NUM_FACES], d_outside[NUM_FACES];
    PetscReal distances[NUM_FACES] = {1.0, 1.0, 0.0, 1.0, 1.0, 1.0};
    PetscInt result = -99;

    PetscFunctionBeginUser;
    for (int v = 0; v < 8; ++v) {
        cell.vertices[v].x = corner[v][0]; cell.vertices[v].y = corner[v][1]; cell.vertices[v].z = corner[v][2];
    }
    MetricBuildCellFacePlanes(&cell, &planes);
    PetscCall(CalculateDistancesToCellFacePlanes(inside_point, &planes, d_inside, 1.0e-12, NULL));
    PetscCall(CalculateDistancesToCellFacePlanes(outside_point, &planes, d_outside, 1.0e-12, NULL));
    PetscCall(PicurvAssertRealNear(0.25, d_inside[BACK], 1.0e-14, "inside point should yield its positive distance to the face"));
    PetscCall(PicurvAssertRealNear(-1.0, d_outside[BACK], 1.0e-14, "outside point should yield its negative distance to the face"));

    PetscCall(DeterminePointPosition(distances, &result));
    PetscCall(PicurvAssertIntEqual(1, result, "one zero face distance should classify as on-face"));

    distances[0] = -0.1;
    PetscCall(DeterminePointPosition(distances, &result));
    PetscCall(PicurvAssertIntEqual(-1, result, "negative face distance should classify as outside"));
    PetscFunctionReturn(0);
}
/**
 * @brief Runs the unit-geometry PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"scalar-interpolation-constant", TestScalarInterpolationConstantField},
        {"vector-interpolation-constant", TestVectorInterpolationConstantField},
        {"signed-distance-and-classification", TestSignedDistanceAndClassification},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv geometry tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-geometry", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
