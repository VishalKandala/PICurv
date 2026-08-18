/**
 * @file statistics_target.c
 * @brief Spatial target resolution for the field-statistics pipeline.
 *
 * Full API contract is documented with the declarations in
 * `include/statistics_target.h`.
 */

#include "statistics_target.h"

/**
 * @brief Implementation of \ref PicurvLayoutDimensionIsNodeLike().
 * @see PicurvLayoutDimensionIsNodeLike()
 */
PetscErrorCode PicurvLayoutDimensionIsNodeLike(FieldLayout layout, PetscInt dim, PetscBool *node_like)
{
    PetscFunctionBeginUser;
    PetscCheck(node_like != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Classification output is required.");
    PetscCheck(dim >= 0 && dim < 3, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Dimension must be 0, 1, or 2; got %" PetscInt_FMT ".", dim);

    switch (layout) {
    case FIELD_LAYOUT_NODE_CENTERED:
        *node_like = PETSC_TRUE;
        break;
    case FIELD_LAYOUT_CELL_CENTERED:
        *node_like = PETSC_FALSE;
        break;
    case FIELD_LAYOUT_I_FACE:
        *node_like = (PetscBool)(dim == 0);
        break;
    case FIELD_LAYOUT_J_FACE:
        *node_like = (PetscBool)(dim == 1);
        break;
    case FIELD_LAYOUT_K_FACE:
        *node_like = (PetscBool)(dim == 2);
        break;
    case FIELD_LAYOUT_COMPONENT_STAGGERED:
        /* x, y, and z live on I-, J-, and K-faces respectively, so no single
         * classification describes the packed vector. */
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
                "Component-staggered layout has no single per-dimension classification.");
    default:
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
                "Unknown field layout %d.", (int)layout);
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper: resolves the layout-valid index span for one dimension.
 * @details Local to this translation unit.
 *
 * The DMDA carries one extra high-side slot beyond the physical grid, and the
 * solver's shifted convention places boundary/dummy values at index zero. Under
 * periodicity the repair algorithms in `Boundaries.c` write index `0` and index
 * `size-1` from the opposite side, so both are dependent duplicates and the
 * independent span starts at one in every layout.
 */
static void ResolveLayoutSpan(PetscBool node_like, PetscBool periodic, PetscInt size,
                              PetscInt *lo, PetscInt *hi_exclusive)
{
    /* Node-like and non-periodic is the only case whose first physical entry
     * sits at index zero; everywhere else index zero is a boundary, dummy, or
     * periodic duplicate. */
    *lo = (node_like && !periodic) ? 0 : 1;
    /* The final slot is the DMDA's extra non-physical entry under every layout. */
    *hi_exclusive = size - 1;
}

/**
 * @brief Implementation of \ref SpatialTargetPlanCreate().
 * @see SpatialTargetPlanCreate()
 */
PetscErrorCode SpatialTargetPlanCreate(UserCtx *user, FieldId field_id,
                                       PicurvStatisticsMask mask, SpatialTargetPlan *plan)
{
    const FieldDescriptor *descriptor = NULL;
    const DMDALocalInfo   *info = NULL;
    SimCtx                *simCtx = NULL;
    PetscInt owned_start[3];
    PetscInt owned_end[3];
    PetscInt global_size[3];
    PetscBool periodic[3];

    PetscFunctionBeginUser;
    PetscCheck(user != NULL && plan != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Block context and plan output are required.");
    PetscCheck(mask == PICURV_STATISTICS_MASK_FLUID, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Only the fluid mask is implemented.");
    simCtx = user->simCtx;
    PetscCheck(simCtx != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "Block context must carry a simulation context.");

    PetscCall(FieldGetDescriptor(field_id, &descriptor));
    PetscCheck(descriptor->layout != FIELD_LAYOUT_COMPONENT_STAGGERED,
               PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "Field '%s' is component-staggered; its components live on different face "
               "families and cannot share one pointwise target domain.",
               descriptor->canonical_name);

    info = &user->info;
    owned_start[0] = info->xs; owned_end[0] = info->xs + info->xm; global_size[0] = info->mx;
    owned_start[1] = info->ys; owned_end[1] = info->ys + info->ym; global_size[1] = info->my;
    owned_start[2] = info->zs; owned_end[2] = info->zs + info->zm; global_size[2] = info->mz;
    /* The simCtx flags are authoritative here because they are what selects
     * DM_BOUNDARY_PERIODIC when the DMDA is built, and therefore what decides
     * whether the end planes are wrapped duplicates. */
    periodic[0] = (PetscBool)(simCtx->i_periodic != 0);
    periodic[1] = (PetscBool)(simCtx->j_periodic != 0);
    periodic[2] = (PetscBool)(simCtx->k_periodic != 0);

    plan->kind = PICURV_TARGET_POINTWISE;
    plan->mask = mask;
    plan->descriptor = descriptor;

    for (PetscInt dim = 0; dim < 3; ++dim) {
        PetscBool node_like = PETSC_FALSE;
        PetscInt  layout_lo = 0;
        PetscInt  layout_hi = 0;

        PetscCall(PicurvLayoutDimensionIsNodeLike(descriptor->layout, dim, &node_like));
        ResolveLayoutSpan(node_like, periodic[dim], global_size[dim], &layout_lo, &layout_hi);

        /* Intersect the layout-valid span with this rank's owned range: the
         * first exclusion removes solver-layout indices, the second removes
         * PETSc halo storage. */
        plan->start[dim] = PetscMax(owned_start[dim], layout_lo);
        plan->end[dim]   = PetscMin(owned_end[dim], layout_hi);
        if (plan->end[dim] < plan->start[dim]) plan->end[dim] = plan->start[dim];
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref SpatialTargetPlanLocalPointCount().
 * @see SpatialTargetPlanLocalPointCount()
 */
PetscErrorCode SpatialTargetPlanLocalPointCount(const SpatialTargetPlan *plan, PetscInt *count)
{
    PetscInt total = 1;

    PetscFunctionBeginUser;
    PetscCheck(plan != NULL && count != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Plan and count output are required.");
    for (PetscInt dim = 0; dim < 3; ++dim) {
        const PetscInt extent = plan->end[dim] - plan->start[dim];
        if (extent <= 0) { *count = 0; PetscFunctionReturn(0); }
        total *= extent;
    }
    *count = total;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref SpatialTargetPlanGlobalPointCount().
 * @see SpatialTargetPlanGlobalPointCount()
 */
PetscErrorCode SpatialTargetPlanGlobalPointCount(const SpatialTargetPlan *plan, MPI_Comm comm, PetscInt *count)
{
    PetscInt local = 0;

    PetscFunctionBeginUser;
    PetscCheck(plan != NULL && count != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Plan and count output are required.");
    PetscCall(SpatialTargetPlanLocalPointCount(plan, &local));
    PetscCallMPI(MPI_Allreduce(&local, count, 1, MPIU_INT, MPI_SUM, comm));
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref SpatialTargetPlanMaskAllows().
 * @see SpatialTargetPlanMaskAllows()
 */
PetscBool SpatialTargetPlanMaskAllows(const SpatialTargetPlan *plan, PetscReal nvert_value)
{
    if (plan == NULL) return PETSC_FALSE;
    return (PetscBool)(nvert_value < PICURV_STATISTICS_FLUID_THRESHOLD);
}
