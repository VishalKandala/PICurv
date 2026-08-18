/**
 * @file statistics_target.h
 * @brief Spatial target resolution for the field-statistics pipeline.
 *
 * Implements the pointwise identity mapping described in @ref p58_target_sec. The
 * abstraction exists with only that one mapping so later spatial bins, profiles,
 * regions, and surfaces extend it rather than retrofit it; see
 * @ref 60_Field_Statistics_Planned_Extensions.
 *
 * The central responsibility is producing an iteration domain that excludes two
 * distinct categories which @ref 56_Field_Identity_and_Layout_Catalog section 4
 * warns must not be conflated:
 *
 * 1. PETSc/MPI halo entries created by DMDA decomposition, and
 * 2. solver-layout boundary, dummy, or duplicate-periodic indices, which exist
 *    even on a single MPI rank.
 *
 * Nothing here resolves field identity or layout metadata itself; that comes
 * from the typed catalog through `FieldGetDescriptor`.
 */

#ifndef PICURV_STATISTICS_TARGET_H
#define PICURV_STATISTICS_TARGET_H

#include "variables.h"
#include "field_catalog.h"

/** @brief Threshold below which a cell counts as fluid for the default mask. */
#define PICURV_STATISTICS_FLUID_THRESHOLD 0.1

/** @brief Spatial mapping kind. Only the pointwise identity is implemented. */
typedef enum {
    PICURV_TARGET_POINTWISE = 0
} PicurvTargetKind;

/** @brief Point-eligibility mask. Only the fluid mask is implemented. */
typedef enum {
    PICURV_STATISTICS_MASK_FLUID = 0
} PicurvStatisticsMask;

/**
 * @brief Resolved iteration domain for one field on one block.
 *
 * Bounds are half-open `[start, end)` in DMDA index order `(i, j, k)` and are
 * already intersected with this rank's owned range, so iterating them touches
 * neither halo storage nor layout boundary/dummy/duplicate-periodic indices.
 * An empty domain on a rank is represented by `end[d] <= start[d]`.
 */
typedef struct {
    PicurvTargetKind        kind;        /**< Spatial mapping; always pointwise. */
    PicurvStatisticsMask    mask;        /**< Point-eligibility mask. */
    const FieldDescriptor  *descriptor;  /**< Catalog metadata for the targeted field. */
    PetscInt                start[3];    /**< Inclusive start per dimension (i, j, k). */
    PetscInt                end[3];      /**< Exclusive end per dimension (i, j, k). */
} SpatialTargetPlan;

/**
 * @brief Reports whether a layout is node-like in one dimension.
 *
 * Node-like dimensions carry one more physical entry than cell-like dimensions,
 * because they sit on grid lines rather than between them. `I_FACE` is node-like
 * in i and cell-like in j and k, and correspondingly for the other face
 * families.
 *
 * @param[in]  layout    Catalog layout.
 * @param[in]  dim       Dimension index: 0 for i, 1 for j, 2 for k.
 * @param[out] node_like Resolved classification.
 * @return Zero on success, or `PETSC_ERR_ARG_OUTOFRANGE` for an invalid
 *         dimension, or `PETSC_ERR_ARG_WRONGSTATE` for a component-staggered
 *         layout, whose components do not share one classification.
 */
PetscErrorCode PicurvLayoutDimensionIsNodeLike(FieldLayout layout, PetscInt dim, PetscBool *node_like);

/**
 * @brief Resolves the iteration domain for one field on one block.
 *
 * Rejects component-staggered fields, whose x, y, and z components live on
 * different face families and therefore cannot share a single pointwise domain.
 * Configured statistics request only cell-centered fields, but the plan resolves
 * node and face layouts correctly so later phases inherit a verified contract.
 *
 * @param[in]  user     Block context supplying the DMDA layout and periodicity.
 * @param[in]  field_id Catalogued field to target.
 * @param[in]  mask     Point-eligibility mask.
 * @param[out] plan     Resolved plan.
 * @return Zero on success, or a PETSc error for a null argument, an unknown
 *         field, or an unsupported layout.
 */
PetscErrorCode SpatialTargetPlanCreate(UserCtx *user, FieldId field_id,
                                       PicurvStatisticsMask mask, SpatialTargetPlan *plan);

/**
 * @brief Counts the points this rank contributes.
 * @param[in]  plan  Resolved plan.
 * @param[out] count Local point count; zero for an empty domain.
 * @return Zero on success, or `PETSC_ERR_ARG_NULL` for a null argument.
 */
PetscErrorCode SpatialTargetPlanLocalPointCount(const SpatialTargetPlan *plan, PetscInt *count);

/**
 * @brief Counts the points contributed across a communicator.
 *
 * The result is decomposition independent: it depends only on the layout,
 * global dimensions, and periodicity, never on how ranks divide the domain.
 *
 * @param[in]  plan  Resolved plan.
 * @param[in]  comm  Communicator to reduce over.
 * @param[out] count Global point count.
 * @return Zero on success, or `PETSC_ERR_ARG_NULL` for a null argument.
 */
PetscErrorCode SpatialTargetPlanGlobalPointCount(const SpatialTargetPlan *plan, MPI_Comm comm, PetscInt *count);

/**
 * @brief Reports whether a point passes the plan's mask.
 *
 * Because `Nvert` changes when immersed bodies move, the mask is treated as
 * potentially moving: callers accumulate a per-point valid count and weight
 * rather than assuming a point contributes to every accepted state.
 *
 * @param[in] plan        Resolved plan.
 * @param[in] nvert_value Node-blanking value at the point.
 * @return `PETSC_TRUE` when the point is eligible.
 */
PetscBool SpatialTargetPlanMaskAllows(const SpatialTargetPlan *plan, PetscReal nvert_value);

#endif /* PICURV_STATISTICS_TARGET_H */
