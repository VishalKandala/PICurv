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

/** @brief Largest reduction buffer a directional average will allocate, in entries. */
#define PICURV_SPATIAL_AVERAGE_MAX_BUFFER 4000000

/**
 * @brief Averages two fields over a target domain and divides the results.
 *
 * Computes `ratio = <numerator> / <denominator>`, where each average is taken
 * separately over the points the plan admits, and then divides. The order matters:
 * the mean of the pointwise quotients is a different number from the quotient of the
 * means, and several callers need the second one. The LES dynamic procedure is one,
 * because Lilly's least-squares closure is defined that way.
 *
 * `average_direction` selects which logical directions collapse into the average.
 * Selecting none makes the operation pointwise, which lets a caller express a local
 * and an averaged formulation through one code path rather than two. Selecting all
 * three yields one number for the whole domain. Selecting a subset retains a profile
 * along the directions left out — averaging over xi and zeta on a plane channel, for
 * example, leaves a wall-normal profile.
 *
 * Weighting is the caller's, not this function's: scale @p numerator and
 * @p denominator by whatever weight the average should carry before calling. A
 * caller wanting a volume-weighted average multiplies both by the cell volume; a
 * caller wanting each admitted point to count once passes them unscaled.
 *
 * The reduction runs over @p comm, so a result is independent of how the domain is
 * decomposed across ranks. The plan already excludes PETSc halo storage and the
 * solver's boundary, dummy, and duplicate-periodic indices, so no physical point is
 * counted twice.
 *
 * @param[in]  user        Block context supplying the DMDA and the solid mask.
 * @param[in]  plan        Resolved iteration domain, from `SpatialTargetPlanCreate()`.
 * @param[in]  numerator   Ghosted local field to average in the numerator.
 * @param[in]  denominator Ghosted local field to average in the denominator, or NULL
 *                         to divide by the number of admitted points.
 * @param[in]  inclusion   Optional second mask: points where this field is not
 *                         positive are skipped. Use it to exclude points that hold a
 *                         zero meaning "never measured" rather than a measurement.
 *                         NULL applies the plan's mask alone.
 * @param[in]  average_direction Directions to average over, in `(xi, eta, zeta)` order.
 * @param[in]  comm        Communicator to reduce over.
 * @param[out] ratio       Ghosted local field receiving the quotient, or NULL when only
 *                         the scalar is wanted. Zero wherever the averaged denominator
 *                         underflows.
 * @param[out] scalar      Receives the single averaged value, or NULL. Valid only when
 *                         every direction is averaged, since otherwise there is no
 *                         single value to report.
 * @return Zero on success, `PETSC_ERR_ARG_OUTOFRANGE` when the retained directions
 *         would need an unreasonably large reduction buffer, or
 *         `PETSC_ERR_ARG_WRONGSTATE` when a scalar is requested from a result that
 *         still varies in space.
 */
PetscErrorCode PicurvSpatialRatioAverage(UserCtx *user, const SpatialTargetPlan *plan,
                                         Vec numerator, Vec denominator, Vec inclusion,
                                         const PetscBool average_direction[3], MPI_Comm comm,
                                         Vec ratio, PetscReal *scalar);

#endif /* PICURV_STATISTICS_TARGET_H */
