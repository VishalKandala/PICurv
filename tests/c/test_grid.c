/**
 * @file test_grid.c
 * @brief C unit tests for grid bounding-box helper routines.
 */

#include "test_support.h"

#include "grid.h"
#include "io.h"

#include <stdlib.h>
/**
 * @brief Tests local bounding-box construction on a uniform grid.
 */

static PetscErrorCode TestComputeLocalBoundingBoxUniformGrid(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    BoundingBox bbox;
    PetscReal expected_max = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(ComputeLocalBoundingBox(user, &bbox));
    expected_max = ((PetscReal)(user->IM - 1) / (PetscReal)user->IM) + 1.0e-6;

    PetscCall(PicurvAssertRealNear(-1.0e-6, bbox.min_coords.x, 1.0e-10, "bbox min x"));
    PetscCall(PicurvAssertRealNear(-1.0e-6, bbox.min_coords.y, 1.0e-10, "bbox min y"));
    PetscCall(PicurvAssertRealNear(-1.0e-6, bbox.min_coords.z, 1.0e-10, "bbox min z"));
    PetscCall(PicurvAssertRealNear(expected_max, bbox.max_coords.x, 1.0e-10, "bbox max x"));
    PetscCall(PicurvAssertRealNear(expected_max, bbox.max_coords.y, 1.0e-10, "bbox max y"));
    PetscCall(PicurvAssertRealNear(expected_max, bbox.max_coords.z, 1.0e-10, "bbox max z"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests collective gather and broadcast of rank-local bounding boxes.
 */

static PetscErrorCode TestGatherAndBroadcastBoundingBoxes(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    BoundingBox local_bbox;
    BoundingBox *boxes = NULL;
    PetscMPIInt rank = 0, size = 1;

    PetscFunctionBeginUser;
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));

    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(ComputeLocalBoundingBox(user, &local_bbox));
    PetscCall(GatherAllBoundingBoxes(user, &boxes));

    if (rank == 0) {
        PetscCall(PicurvAssertBool((PetscBool)(boxes != NULL), "root must receive gathered bounding boxes"));
    } else {
        PetscCall(PicurvAssertBool((PetscBool)(boxes == NULL), "non-root should have NULL pre-broadcast bbox pointer"));
    }

    PetscCall(BroadcastAllBoundingBoxes(user, &boxes));
    PetscCall(PicurvAssertBool((PetscBool)(boxes != NULL), "all ranks must have bbox array after broadcast"));

    for (PetscMPIInt r = 0; r < size; ++r) {
        PetscCall(PicurvAssertBool((PetscBool)(boxes[r].min_coords.x <= boxes[r].max_coords.x), "bbox x range is valid"));
        PetscCall(PicurvAssertBool((PetscBool)(boxes[r].min_coords.y <= boxes[r].max_coords.y), "bbox y range is valid"));
        PetscCall(PicurvAssertBool((PetscBool)(boxes[r].min_coords.z <= boxes[r].max_coords.z), "bbox z range is valid"));
    }

    PetscCall(PicurvAssertRealNear(local_bbox.min_coords.x, boxes[rank].min_coords.x, 1.0e-10, "rank-local bbox min x"));
    PetscCall(PicurvAssertRealNear(local_bbox.max_coords.x, boxes[rank].max_coords.x, 1.0e-10, "rank-local bbox max x"));

    free(boxes);
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}
/**
 * @brief Tests that programmatic domain bounds are non-dimensionalized by L_ref.
 *
 * A file-backed grid is divided by the reference length when the conductor publishes
 * it as an asset; a programmatic grid has no such staging step, so ReadGridGenerationInputs()
 * must perform the same division so all three grid modes take dimensional user input
 * and the postprocessor's unconditional `Coordinates * L_ref` round-trips correctly.
 */
static PetscErrorCode TestReadGridGenerationInputsNonDimensionalizesBounds(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));

    PetscCall(PetscOptionsSetValue(NULL, "-im", "4"));
    PetscCall(PetscOptionsSetValue(NULL, "-jm", "4"));
    PetscCall(PetscOptionsSetValue(NULL, "-km", "4"));
    PetscCall(PetscOptionsSetValue(NULL, "-xMins", "0.0"));
    PetscCall(PetscOptionsSetValue(NULL, "-xMaxs", "0.2"));
    PetscCall(PetscOptionsSetValue(NULL, "-yMins", "0.0"));
    PetscCall(PetscOptionsSetValue(NULL, "-yMaxs", "0.2"));
    PetscCall(PetscOptionsSetValue(NULL, "-zMins", "0.0"));
    PetscCall(PetscOptionsSetValue(NULL, "-zMaxs", "0.4"));

    /* Authored in physical units (metres); L_ref converts to the solver units the
     * simulator evolves in. 0.2 m at L_ref=0.05 is a domain of 4 solver units, the same
     * shape a nondimensional-input user would have written by hand. */
    simCtx->scaling.L_ref = 0.05;
    PetscCall(ReadGridGenerationInputs(user));

    PetscCall(PicurvAssertRealNear(0.0, user->Min_X, 1.0e-12, "xMin divides by L_ref"));
    PetscCall(PicurvAssertRealNear(4.0, user->Max_X, 1.0e-12, "0.2 / 0.05 = 4.0 solver units"));
    PetscCall(PicurvAssertRealNear(4.0, user->Max_Y, 1.0e-12, "yMax non-dimensionalized"));
    PetscCall(PicurvAssertRealNear(8.0, user->Max_Z, 1.0e-12, "0.4 / 0.05 = 8.0 solver units"));

    /* L_ref = 1.0 is a no-op, so an author already working in solver units is
     * unaffected - this is the common case and must not silently change behaviour. */
    simCtx->scaling.L_ref = 1.0;
    PetscCall(ReadGridGenerationInputs(user));
    PetscCall(PicurvAssertRealNear(0.2, user->Max_X, 1.0e-12, "L_ref=1.0 leaves bounds unchanged"));

    /* A non-positive reference length is refused rather than dividing by zero or
     * flipping the domain's sign silently. */
    {
        PetscErrorCode bad = 0;

        simCtx->scaling.L_ref = 0.0;
        PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
        bad = ReadGridGenerationInputs(user);
        PetscCall(PetscPopErrorHandler());
        PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_OUTOFRANGE, bad,
                                       "L_ref <= 0 is refused rather than dividing by zero"));
    }

    PetscCall(PetscOptionsClearValue(NULL, "-im"));
    PetscCall(PetscOptionsClearValue(NULL, "-jm"));
    PetscCall(PetscOptionsClearValue(NULL, "-km"));
    PetscCall(PetscOptionsClearValue(NULL, "-xMins"));
    PetscCall(PetscOptionsClearValue(NULL, "-xMaxs"));
    PetscCall(PetscOptionsClearValue(NULL, "-yMins"));
    PetscCall(PetscOptionsClearValue(NULL, "-yMaxs"));
    PetscCall(PetscOptionsClearValue(NULL, "-zMins"));
    PetscCall(PetscOptionsClearValue(NULL, "-zMaxs"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Runs the unit-grid PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"compute-local-bbox-uniform-grid", TestComputeLocalBoundingBoxUniformGrid},
        {"gather-and-broadcast-bboxes", TestGatherAndBroadcastBoundingBoxes},
        {"read-grid-generation-inputs-nondimensionalizes-bounds",
         TestReadGridGenerationInputsNonDimensionalizesBounds},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv grid tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-grid", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
