/**
 * @file statistics_accumulator.h
 * @brief Per-window PETSc accumulator storage and pointwise application.
 *
 * Holds the independent state each named window owns, and applies one accepted
 * completed state to it, per
 * @ref 60_Field_Statistics_Phase2_Implementation_Specification sections 5 and 13.
 *
 * Storage is allocated once by the vector factory and released once at teardown.
 * Application is strictly pointwise: it reads a source field value at an owned
 * point and updates that point's accumulator, so it performs no halo exchange and
 * allocates nothing.
 */

#ifndef PICURV_STATISTICS_ACCUMULATOR_H
#define PICURV_STATISTICS_ACCUMULATOR_H

#include "variables.h"
#include "statistics_window.h"

/**
 * @brief Independent accumulator state for one window on one block.
 *
 * Per-point occupancy is tracked separately from the field moments because the
 * fluid mask can move: a point contributes only to the states in which it was
 * valid, so its own count and weight are what normalize its moments.
 *
 * Product components are stored one degree of freedom per vector, which keeps
 * every accumulator on a DM the factory already built and maps one-to-one onto
 * the checkpoint payload inventory.
 */
typedef struct PicurvWindowStorage {
    PetscInt  field_count;       /**< Fields accumulated. */
    PetscInt  covariance_count;  /**< Covariance pairs accumulated. */
    Vec       count;             /**< Per-point accepted sample count. */
    Vec       weight;            /**< Per-point valid weight. */
    Vec       weight_sq;         /**< Per-point squared-weight sum. */
    Vec      *mean;              /**< One per field, matching that field's layout. */
    PetscInt *m2_offset;         /**< field_count+1 offsets into @c m2. */
    Vec      *m2;                /**< Flattened centered second-moment components. */
    PetscInt *cm_offset;         /**< covariance_count+1 offsets into @c cm. */
    Vec      *cm;                /**< Flattened centered co-moment components. */
} PicurvWindowStorage;

/**
 * @brief Reports how many symmetric product components a field's second moment needs.
 * @param[in]  dof   Degree of freedom of the field.
 * @param[out] count Component count: one for a scalar, six for a three-vector.
 * @return Zero on success, or `PETSC_ERR_ARG_OUTOFRANGE` for an unsupported dof.
 */
PetscErrorCode PicurvProductComponentCount(PetscInt dof, PetscInt *count);

/**
 * @brief Reports how many components a covariance between two fields needs.
 * @param[in]  dof_a Degree of freedom of the first member.
 * @param[in]  dof_b Degree of freedom of the second member.
 * @param[out] count Component count.
 * @return Zero on success, or `PETSC_ERR_ARG_OUTOFRANGE` for an unsupported pairing.
 */
PetscErrorCode PicurvCovarianceComponentCount(PetscInt dof_a, PetscInt dof_b, PetscInt *count);

/**
 * @brief Allocates the accumulator state one window owns on one block.
 *
 * Every vector is duplicated from one the factory already built, so no new DM or
 * layout decision is introduced.
 *
 * @param[in]  user       Block context supplying the source fields.
 * @param[in]  definition Window definition naming the requested fields and pairs.
 * @param[out] storage    Storage to populate; zeroed on entry.
 * @return Zero on success, or a PETSc error for an unknown field or unsupported layout.
 */
PetscErrorCode PicurvWindowStorageCreate(UserCtx *user, const PicurvWindowDefinition *definition,
                                         PicurvWindowStorage *storage);

/**
 * @brief Releases accumulator state previously created for one window.
 * @param[in,out] storage Storage to release; safe to call on zeroed storage.
 * @return Zero on success.
 */
PetscErrorCode PicurvWindowStorageDestroy(PicurvWindowStorage *storage);

/**
 * @brief Applies one accepted completed state to a window's accumulators.
 *
 * Iterates the pointwise target domain, skipping points the mask rejects, and
 * updates each point's occupancy and every requested moment and co-moment through
 * the centered kernels.
 *
 * @param[in]     user       Block context supplying the source fields.
 * @param[in]     definition Window definition.
 * @param[in,out] storage    Accumulator state to update.
 * @param[in]     weight     Weight the window assigned to this state; must be positive.
 * @return Zero on success, or a PETSc error.
 */
PetscErrorCode PicurvWindowAccumulate(UserCtx *user, const PicurvWindowDefinition *definition,
                                      PicurvWindowStorage *storage, PetscReal weight);

#endif /* PICURV_STATISTICS_ACCUMULATOR_H */
