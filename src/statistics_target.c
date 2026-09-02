/**
 * @file statistics_target.c
 * @brief Spatial target resolution for the field-statistics pipeline.
 *
 * Full API contract is documented with the declarations in
 * `include/statistics_target.h`.
 */

#include "statistics_target.h"
#include "logging.h"

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

#undef __FUNCT__
#define __FUNCT__ "SpatialTargetPlanCreate"
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
    PROFILE_FUNCTION_BEGIN;
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
    PROFILE_FUNCTION_END;
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

#undef __FUNCT__
#define __FUNCT__ "SpatialTargetPlanGlobalPointCount"
/**
 * @brief Implementation of \ref SpatialTargetPlanGlobalPointCount().
 * @see SpatialTargetPlanGlobalPointCount()
 */
PetscErrorCode SpatialTargetPlanGlobalPointCount(const SpatialTargetPlan *plan, MPI_Comm comm, PetscInt *count)
{
    PetscInt local = 0;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCheck(plan != NULL && count != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Plan and count output are required.");
    PetscCall(SpatialTargetPlanLocalPointCount(plan, &local));
    PetscCallMPI(MPI_Allreduce(&local, count, 1, MPIU_INT, MPI_SUM, comm));
    PROFILE_FUNCTION_END;
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

/**
 * @brief Flattens the retained global indices of one point into a reduction slot.
 *
 * A direction being averaged over collapses to a single entry, so the buffer is
 * indexed only by the directions left out.
 */
static PetscInt SpatialAverageSlot(const PetscBool average_direction[3],
                                   const PetscInt retained_extent[3],
                                   PetscInt i, PetscInt j, PetscInt k)
{
    return (average_direction[0] ? 0 : i) +
           retained_extent[0] * ((average_direction[1] ? 0 : j) +
           retained_extent[1] *  (average_direction[2] ? 0 : k));
}

/**
 * @brief Implementation of \ref PicurvSpatialRatioAverage().
 * @see PicurvSpatialRatioAverage()
 */
PetscErrorCode PicurvSpatialRatioAverage(UserCtx *user, const SpatialTargetPlan *plan,
                                         Vec numerator, Vec denominator, Vec inclusion,
                                         const PetscBool average_direction[3], MPI_Comm comm,
                                         Vec ratio, PetscReal *scalar)
{
    DM             da = NULL;
    DMDALocalInfo  info;
    PetscInt       retained_extent[3];
    PetscInt       buffer_size = 1;
    PetscInt       averaged_count = 0;
    PetscReal   ***num = NULL, ***den = NULL, ***inc = NULL, ***out = NULL, ***nvert = NULL;
    PetscReal     *num_sum = NULL, *den_sum = NULL;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCheck(user != NULL && plan != NULL && numerator != NULL && average_direction != NULL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Context, plan, numerator, and direction selector are required.");

    da = user->da;
    PetscCall(DMDAGetLocalInfo(da, &info));

    retained_extent[0] = average_direction[0] ? 1 : info.mx;
    retained_extent[1] = average_direction[1] ? 1 : info.my;
    retained_extent[2] = average_direction[2] ? 1 : info.mz;
    for (PetscInt axis = 0; axis < 3; axis++) {
        if (average_direction[axis]) averaged_count++;
        buffer_size *= retained_extent[axis];
    }
    PetscCheck(scalar == NULL || averaged_count == 3, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "A single averaged value exists only when every direction is averaged over; "
               "%" PetscInt_FMT " of 3 were selected.", averaged_count);

    PetscCall(DMDAVecGetArrayRead(da, numerator, &num));
    if (denominator) PetscCall(DMDAVecGetArrayRead(da, denominator, &den));
    if (inclusion) PetscCall(DMDAVecGetArrayRead(da, inclusion, &inc));
    if (ratio) PetscCall(DMDAVecGetArray(da, ratio, &out));
    PetscCall(DMDAVecGetArrayRead(da, user->lNvert, &nvert));

    if (averaged_count == 0) {
        /* Pointwise: the empty averaging set, which is a legitimate request rather
         * than a degenerate one. It is what a local formulation asks for. The masks
         * mean the same thing here as they do for an averaged set - a point they
         * exclude contributes nothing and receives nothing - so that a caller can
         * switch between the two without its exclusions quietly changing meaning. */
        for (PetscInt k = plan->start[2]; k < plan->end[2]; k++)
        for (PetscInt j = plan->start[1]; j < plan->end[1]; j++)
        for (PetscInt i = plan->start[0]; i < plan->end[0]; i++) {
            const PetscReal divisor = denominator ? den[k][j][i] : 1.0;

            if (!SpatialTargetPlanMaskAllows(plan, nvert[k][j][i]) ||
                (inc && inc[k][j][i] <= 0.0)) {
                out[k][j][i] = 0.0;
                continue;
            }
            out[k][j][i] = (PetscAbsReal(divisor) > 0.0) ? (num[k][j][i] / divisor) : 0.0;
        }
    } else {
        PetscCheck(buffer_size <= PICURV_SPATIAL_AVERAGE_MAX_BUFFER, PETSC_COMM_SELF,
                   PETSC_ERR_ARG_OUTOFRANGE,
                   "Averaging over %" PetscInt_FMT " direction(s) would need a reduction buffer "
                   "of %" PetscInt_FMT " entries, above the %d entry limit. Average over more "
                   "directions.", averaged_count, buffer_size, PICURV_SPATIAL_AVERAGE_MAX_BUFFER);

        PetscCall(PetscCalloc2(buffer_size, &num_sum, buffer_size, &den_sum));

        for (PetscInt k = plan->start[2]; k < plan->end[2]; k++)
        for (PetscInt j = plan->start[1]; j < plan->end[1]; j++)
        for (PetscInt i = plan->start[0]; i < plan->end[0]; i++) {
            PetscInt slot;

            if (!SpatialTargetPlanMaskAllows(plan, nvert[k][j][i])) continue;
            /* A point the caller marks unmeasured holds a zero that means absence, not
             * a measurement; counting it would scale the answer toward zero. */
            if (inc && inc[k][j][i] <= 0.0) continue;

            slot = SpatialAverageSlot(average_direction, retained_extent, i, j, k);
            num_sum[slot] += num[k][j][i];
            den_sum[slot] += denominator ? den[k][j][i] : 1.0;
        }

        PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, num_sum, (PetscMPIInt)buffer_size,
                                   MPIU_REAL, MPI_SUM, comm));
        PetscCallMPI(MPI_Allreduce(MPI_IN_PLACE, den_sum, (PetscMPIInt)buffer_size,
                                   MPIU_REAL, MPI_SUM, comm));

        if (out) {
            for (PetscInt k = plan->start[2]; k < plan->end[2]; k++)
            for (PetscInt j = plan->start[1]; j < plan->end[1]; j++)
            for (PetscInt i = plan->start[0]; i < plan->end[0]; i++) {
                const PetscInt slot = SpatialAverageSlot(average_direction, retained_extent, i, j, k);

                out[k][j][i] = (PetscAbsReal(den_sum[slot]) > 0.0)
                                   ? (num_sum[slot] / den_sum[slot]) : 0.0;
            }
        }
        if (scalar) {
            *scalar = (PetscAbsReal(den_sum[0]) > 0.0) ? (num_sum[0] / den_sum[0]) : 0.0;
        }

        PetscCall(PetscFree2(num_sum, den_sum));
    }

    PetscCall(DMDAVecRestoreArrayRead(da, numerator, &num));
    if (denominator) PetscCall(DMDAVecRestoreArrayRead(da, denominator, &den));
    if (inclusion) PetscCall(DMDAVecRestoreArrayRead(da, inclusion, &inc));
    if (ratio) PetscCall(DMDAVecRestoreArray(da, ratio, &out));
    PetscCall(DMDAVecRestoreArrayRead(da, user->lNvert, &nvert));

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
