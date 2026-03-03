#include "test_support.h"

#include "interpolation.h"
#include "walkingsearch.h"

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

static PetscErrorCode TestSignedDistanceAndClassification(void)
{
    const Cmpnts v1 = {0.0, 0.0, 0.0};
    const Cmpnts v2 = {1.0, 0.0, 0.0};
    const Cmpnts v3 = {1.0, 1.0, 0.0};
    const Cmpnts v4 = {0.0, 1.0, 0.0};
    const Cmpnts centroid = {0.5, 0.5, 0.5};
    const Cmpnts inside_point = {0.5, 0.5, 0.25};
    const Cmpnts outside_point = {0.5, 0.5, -1.0};
    PetscReal d_inside = 0.0, d_outside = 0.0;
    PetscReal distances[NUM_FACES] = {1.0, 1.0, 0.0, 1.0, 1.0, 1.0};
    PetscInt result = -99;

    PetscFunctionBeginUser;
    PetscCall(ComputeSignedDistanceToPlane(v1, v2, v3, v4, centroid, inside_point, &d_inside, 1.0e-12));
    PetscCall(ComputeSignedDistanceToPlane(v1, v2, v3, v4, centroid, outside_point, &d_outside, 1.0e-12));
    PetscCall(PicurvAssertBool((PetscBool)(d_inside > 0.0), "inside point should yield positive signed distance"));
    PetscCall(PicurvAssertBool((PetscBool)(d_outside < 0.0), "outside point should yield negative signed distance"));

    PetscCall(DeterminePointPosition(distances, &result));
    PetscCall(PicurvAssertIntEqual(1, result, "one zero face distance should classify as on-face"));

    distances[0] = -0.1;
    PetscCall(DeterminePointPosition(distances, &result));
    PetscCall(PicurvAssertIntEqual(-1, result, "negative face distance should classify as outside"));
    PetscFunctionReturn(0);
}

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
