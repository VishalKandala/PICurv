/**
 * @file test_statistics_target.c
 * @brief C unit tests for pointwise spatial target resolution.
 *
 * Covers the Stage 2 acceptance items in
 * @ref 60_Field_Statistics_Phase2_Implementation_Specification section 13: cell,
 * node, and I/J/K-face layouts, and nonperiodic, mixed, and fully periodic
 * domains with no duplicate planes.
 *
 * The fixture builds a DMDA of size `n + 1` per dimension, so with `n = 6` each
 * dimension has seven slots: one extra non-physical high-side slot, and under the
 * solver's shifted convention a boundary/dummy slot at index zero.
 */

#include "test_support.h"

#include "statistics_target.h"

#define TARGET_TEST_N 6            /* fixture size; DMDA is TARGET_TEST_N + 1 per dimension */
#define TARGET_CELL_EXTENT 5       /* cell-like span: indices [1, 6) */
#define TARGET_NODE_EXTENT 6       /* node-like nonperiodic span: indices [0, 6) */

/** @brief Builds a plan and returns its global point count. */
static PetscErrorCode PlanCount(UserCtx *user, FieldId field_id, PetscInt *count)
{
    SpatialTargetPlan plan;

    PetscFunctionBeginUser;
    PetscCall(SpatialTargetPlanCreate(user, field_id, PICURV_STATISTICS_MASK_FLUID, &plan));
    PetscCall(SpatialTargetPlanGlobalPointCount(&plan, PETSC_COMM_WORLD, count));
    PetscFunctionReturn(0);
}

/** @brief Cell-centered fields must exclude the dummy slot and the extra high slot. */
static PetscErrorCode TestCellCenteredExcludesBoundaryAndDummy(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    SpatialTargetPlan plan;
    PetscInt local = 0, global = 0;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, TARGET_TEST_N, TARGET_TEST_N, TARGET_TEST_N));

    PetscCall(SpatialTargetPlanCreate(user, FIELD_ID_P, PICURV_STATISTICS_MASK_FLUID, &plan));
    PetscCall(PicurvAssertIntEqual(PICURV_TARGET_POINTWISE, plan.kind, "plan kind should be pointwise"));
    PetscCall(PicurvAssertIntEqual(FIELD_LAYOUT_CELL_CENTERED, plan.descriptor->layout,
                                   "P should resolve as cell-centered"));
    /* On one rank the owned range spans the block, so the plan is the layout span. */
    if (simCtx->size == 1) {
        for (PetscInt dim = 0; dim < 3; ++dim) {
            PetscCall(PicurvAssertIntEqual(1, plan.start[dim],
                                           "cell-centered span must skip the index-zero dummy slot"));
            PetscCall(PicurvAssertIntEqual(TARGET_TEST_N, plan.end[dim],
                                           "cell-centered span must skip the extra high-side slot"));
        }
    }
    PetscCall(SpatialTargetPlanLocalPointCount(&plan, &local));
    PetscCall(SpatialTargetPlanGlobalPointCount(&plan, PETSC_COMM_WORLD, &global));
    PetscCall(PicurvAssertIntEqual(TARGET_CELL_EXTENT * TARGET_CELL_EXTENT * TARGET_CELL_EXTENT, global,
                                   "cell-centered global point count"));
    PetscCall(PicurvAssertBool((PetscBool)(local <= global), "local count cannot exceed global count"));

    /* Vector and scalar cell-centered fields share one domain. */
    PetscCall(PlanCount(user, FIELD_ID_UCAT, &global));
    PetscCall(PicurvAssertIntEqual(TARGET_CELL_EXTENT * TARGET_CELL_EXTENT * TARGET_CELL_EXTENT, global,
                                   "Ucat shares the cell-centered domain regardless of dof"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief Each face family is node-like in exactly its own direction. */
static PetscErrorCode TestFaceLayoutsAreNodeLikeInOwnDirection(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscInt global = 0;
    const PetscInt cell = TARGET_CELL_EXTENT;
    const PetscInt node = TARGET_NODE_EXTENT;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, TARGET_TEST_N, TARGET_TEST_N, TARGET_TEST_N));

    PetscCall(PlanCount(user, FIELD_ID_CSI, &global));
    PetscCall(PicurvAssertIntEqual(node * cell * cell, global, "I-face domain is node-like in i only"));
    PetscCall(PlanCount(user, FIELD_ID_ETA, &global));
    PetscCall(PicurvAssertIntEqual(cell * node * cell, global, "J-face domain is node-like in j only"));
    PetscCall(PlanCount(user, FIELD_ID_ZET, &global));
    PetscCall(PicurvAssertIntEqual(cell * cell * node, global, "K-face domain is node-like in k only"));
    PetscCall(PlanCount(user, FIELD_ID_COORDINATES, &global));
    PetscCall(PicurvAssertIntEqual(node * node * node, global, "node-centered domain is node-like in all directions"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Periodic directions must drop the wrapped duplicate plane.
 *
 * The periodic repair algorithms write index zero and the final slot from the
 * opposite side, so both are dependent duplicates. A node-like direction
 * therefore loses one entry under periodicity, while a cell-like direction is
 * unchanged because its duplicates were already outside the span.
 */
static PetscErrorCode TestPeriodicDirectionsDropDuplicatePlane(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    PetscInt global = 0;
    const PetscInt cell = TARGET_CELL_EXTENT;
    const PetscInt node = TARGET_NODE_EXTENT;

    PetscFunctionBeginUser;
    /* Mixed: periodic in i only. */
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user,
                                                         TARGET_TEST_N, TARGET_TEST_N, TARGET_TEST_N,
                                                         PETSC_TRUE, PETSC_FALSE, PETSC_FALSE));
    PetscCall(PlanCount(user, FIELD_ID_P, &global));
    PetscCall(PicurvAssertIntEqual(cell * cell * cell, global,
                                   "cell-centered count is unchanged by periodicity"));
    PetscCall(PlanCount(user, FIELD_ID_CSI, &global));
    PetscCall(PicurvAssertIntEqual(cell * cell * cell, global,
                                   "I-face loses its duplicate plane when i is periodic"));
    PetscCall(PlanCount(user, FIELD_ID_ETA, &global));
    PetscCall(PicurvAssertIntEqual(cell * node * cell, global,
                                   "J-face keeps its node-like span when only i is periodic"));
    PetscCall(PlanCount(user, FIELD_ID_COORDINATES, &global));
    PetscCall(PicurvAssertIntEqual(cell * node * node, global,
                                   "node-centered loses only the periodic direction"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));

    /* Fully periodic: every node-like direction collapses to the cell extent. */
    PetscCall(PicurvCreateMinimalContextsWithPeriodicity(&simCtx, &user,
                                                         TARGET_TEST_N, TARGET_TEST_N, TARGET_TEST_N,
                                                         PETSC_TRUE, PETSC_TRUE, PETSC_TRUE));
    PetscCall(PlanCount(user, FIELD_ID_COORDINATES, &global));
    PetscCall(PicurvAssertIntEqual(cell * cell * cell, global,
                                   "fully periodic node-centered domain has no duplicate planes"));
    PetscCall(PlanCount(user, FIELD_ID_CSI, &global));
    PetscCall(PicurvAssertIntEqual(cell * cell * cell, global,
                                   "fully periodic I-face domain has no duplicate planes"));
    PetscCall(PlanCount(user, FIELD_ID_P, &global));
    PetscCall(PicurvAssertIntEqual(cell * cell * cell, global,
                                   "fully periodic cell-centered domain is unchanged"));
    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief Component-staggered fields cannot share one pointwise domain and are rejected. */
static PetscErrorCode TestComponentStaggeredRejected(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    SpatialTargetPlan plan;
    PetscErrorCode staggered_ierr = 0;
    PetscErrorCode classify_ierr = 0;
    PetscBool node_like = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, TARGET_TEST_N, TARGET_TEST_N, TARGET_TEST_N));

    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    staggered_ierr = SpatialTargetPlanCreate(user, FIELD_ID_UCONT, PICURV_STATISTICS_MASK_FLUID, &plan);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_WRONGSTATE, staggered_ierr,
                                   "Ucont is component-staggered and must not resolve a pointwise domain"));

    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    classify_ierr = PicurvLayoutDimensionIsNodeLike(FIELD_LAYOUT_COMPONENT_STAGGERED, 0, &node_like);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_WRONGSTATE, classify_ierr,
                                   "component-staggered layout has no single per-dimension classification"));

    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    classify_ierr = PicurvLayoutDimensionIsNodeLike(FIELD_LAYOUT_CELL_CENTERED, 3, &node_like);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_OUTOFRANGE, classify_ierr,
                                   "dimension index outside 0..2 must be rejected"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/** @brief The fluid mask admits fluid cells and rejects blanked ones at the documented threshold. */
static PetscErrorCode TestFluidMaskThreshold(void)
{
    SimCtx *simCtx = NULL;
    UserCtx *user = NULL;
    SpatialTargetPlan plan;

    PetscFunctionBeginUser;
    PetscCall(PicurvCreateMinimalContexts(&simCtx, &user, TARGET_TEST_N, TARGET_TEST_N, TARGET_TEST_N));
    PetscCall(SpatialTargetPlanCreate(user, FIELD_ID_UCAT, PICURV_STATISTICS_MASK_FLUID, &plan));

    PetscCall(PicurvAssertBool(SpatialTargetPlanMaskAllows(&plan, 0.0),
                               "open fluid must pass the fluid mask"));
    PetscCall(PicurvAssertBool(SpatialTargetPlanMaskAllows(&plan, 0.05),
                               "values below the threshold must pass the fluid mask"));
    PetscCall(PicurvAssertBool((PetscBool)!SpatialTargetPlanMaskAllows(&plan, 1.0),
                               "solid cells must fail the fluid mask"));
    PetscCall(PicurvAssertBool((PetscBool)!SpatialTargetPlanMaskAllows(&plan, PICURV_STATISTICS_FLUID_THRESHOLD),
                               "the threshold itself is excluded, matching the existing Nvert < 0.1 rule"));
    PetscCall(PicurvAssertBool((PetscBool)!SpatialTargetPlanMaskAllows(NULL, 0.0),
                               "a null plan must not admit points"));

    PetscCall(PicurvDestroyMinimalContexts(&simCtx, &user));
    PetscFunctionReturn(0);
}

/**
 * @brief Entry point for the spatial target suite.
 */
int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"cell-centered-excludes-boundary-and-dummy", TestCellCenteredExcludesBoundaryAndDummy},
        {"face-layouts-node-like-in-own-direction", TestFaceLayoutsAreNodeLikeInOwnDirection},
        {"periodic-directions-drop-duplicate-plane", TestPeriodicDirectionsDropDuplicatePlane},
        {"component-staggered-rejected", TestComponentStaggeredRejected},
        {"fluid-mask-threshold", TestFluidMaskThreshold},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv spatial target tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-statistics-target", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
