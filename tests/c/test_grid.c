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
 * @brief Tests that the coarsest split puts each refined finest-level boundary nearest its share.
 * @details 128 cells on 4 ranks at 3 levels: PETSc's even coarsest split (9,9,8,8 of 34
 *          points) refines to finest owned cells 33/36/32/27; 9,8,8,9 refines to 33/32/32/31,
 *          the best refinement allows. Unrefined, the last rank takes both spare points of 130,
 *          32/32/32/34, so each owns 32 cells, which PETSc's 33/33/32/32 does not.
 */
static PetscErrorCode TestBalanceCoarsestOwnershipBalancesFinestLevel(void)
{
    const PetscInt refined[4] = {9, 8, 8, 9}, unrefined[4] = {25, 25, 25, 24}, one_level[4] = {32, 32, 32, 34};
    PetscInt       l[4];

    PetscFunctionBeginUser;
    PetscCall(BalanceCoarsestOwnership(34, 128, 2, 4, l));
    for (PetscInt r = 0; r < 4; ++r) PetscCall(PicurvAssertIntEqual(refined[r], l[r], "3-level coarsest split"));
    PetscCall(BalanceCoarsestOwnership(99, 97, 0, 4, l));
    for (PetscInt r = 0; r < 4; ++r) PetscCall(PicurvAssertIntEqual(unrefined[r], l[r], "unrefined axis keeps PETSc's split"));
    PetscCall(BalanceCoarsestOwnership(130, 128, 0, 4, l));
    for (PetscInt r = 0; r < 4; ++r) PetscCall(PicurvAssertIntEqual(one_level[r], l[r], "one level gives the last rank both spare points"));
    PetscFunctionReturn(0);
}

/**
 * @brief Tests that a grid too large for 32-bit PetscInt is refused before its DMs are created.
 * @details 1300^3 cells give 1301^3 nodes, whose 3-component fields need about 6.6e9 entries.
 *          A 64-bit-index build indexes them, so the test applies to 32-bit builds only.
 */
static PetscErrorCode TestOversizedGridIsRefusedOn32BitIndices(void)
{
    SimCtx        *simCtx = NULL;
    char           tmpdir[PETSC_MAX_PATH_LEN] = "";
    PetscErrorCode setup_ierr;

    PetscFunctionBeginUser;
    if (sizeof(PetscInt) > 4) PetscFunctionReturn(0);
    PetscCall(PetscPushErrorHandler(PetscReturnErrorHandler, NULL));
    PetscCall(PicurvCheckPetscIntRange(PETSC_COMM_WORLD, (PetscInt64)PETSC_MAX_INT, "an object at the limit"));
    setup_ierr = PicurvCheckPetscIntRange(PETSC_COMM_WORLD, (PetscInt64)PETSC_MAX_INT + 1, "an object past the limit");
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_OUTOFRANGE, setup_ierr, "one entry past PETSC_MAX_INT must be refused"));

    PetscCall(PetscPushErrorHandler(PetscReturnErrorHandler, NULL));
    setup_ierr = PicurvBuildTinyRuntimeContextWithOptions(NULL, PETSC_FALSE, "-im 1300\n-jm 1300\n-km 1300\n-mg_level 1\n",
                                                          &simCtx, NULL, tmpdir, sizeof(tmpdir));
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_OUTOFRANGE, setup_ierr,
                                   "grid setup must refuse a grid whose fields outgrow 32-bit PetscInt"));
    if (simCtx) PetscCall(PicurvDestroyRuntimeContext(&simCtx));
    if (tmpdir[0]) PetscCall(PicurvRemoveTempDir(tmpdir));
    PetscCall(PetscOptionsClear(NULL));
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
        {"balance-coarsest-ownership-balances-finest-level", TestBalanceCoarsestOwnershipBalancesFinestLevel},
        {"oversized-grid-is-refused-on-32-bit-indices", TestOversizedGridIsRefusedOn32BitIndices},
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
