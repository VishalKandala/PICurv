/**
 * @file test_grid.c
 * @brief C unit tests for grid bounding-box helper routines.
 */

#include "test_support.h"

#include "grid.h"

#include <stdlib.h>
/**
 * @brief Tests local bounding-box construction on a uniform grid.
 */

static PetscErrorCode TestComputeLocalBoundingBoxUniformGrid(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    BoundingBox bbox;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, 4, 4, 4));
    PetscCall(ComputeLocalBoundingBox(user, &bbox));

    PetscCall(PicurvAssertRealNear(-1.0e-6, bbox.min_coords.x, 1.0e-10, "bbox min x"));
    PetscCall(PicurvAssertRealNear(-1.0e-6, bbox.min_coords.y, 1.0e-10, "bbox min y"));
    PetscCall(PicurvAssertRealNear(-1.0e-6, bbox.min_coords.z, 1.0e-10, "bbox min z"));
    PetscCall(PicurvAssertRealNear(3.0 + 1.0e-6, bbox.max_coords.x, 1.0e-10, "bbox max x"));
    PetscCall(PicurvAssertRealNear(3.0 + 1.0e-6, bbox.max_coords.y, 1.0e-10, "bbox max y"));
    PetscCall(PicurvAssertRealNear(3.0 + 1.0e-6, bbox.max_coords.z, 1.0e-10, "bbox max z"));

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
 * @brief Runs the unit-grid PETSc test binary.
 */

int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"compute-local-bbox-uniform-grid", TestComputeLocalBoundingBoxUniformGrid},
        {"gather-and-broadcast-bboxes", TestGatherAndBroadcastBoundingBoxes},
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
