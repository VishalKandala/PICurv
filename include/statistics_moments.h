/**
 * @file statistics_moments.h
 * @brief Weighted centered-moment kernels for the field-statistics pipeline.
 *
 * These are pure numerical kernels: no configuration, no PETSc vectors, and no
 * knowledge of windows, fields, or scheduling. They implement the accumulator
 * contract fixed in
 * @ref 60_Field_Statistics_Phase2_Implementation_Specification section 5, which
 * in turn follows @ref 58_Turbulence_Statistics_Pipeline_Specification section 7.
 *
 * Centered state is used rather than raw sums because it retains everything the
 * first and second moments need while staying numerically stable for
 * high-mean/low-fluctuation signals. Raw sums are recoverable as
 * `S_x = W*mu` and `Q_xx = M2 + W*mu^2`.
 */

#ifndef PICURV_STATISTICS_MOMENTS_H
#define PICURV_STATISTICS_MOMENTS_H

#include <petscsys.h>

/**
 * @brief Weighted centered state for one scalar quantity at one point.
 *
 * `weight_sq` is retained because physical-time weighting produces unequal
 * weights, so an effective sample size cannot be inferred from `count` alone.
 */
typedef struct {
    PetscReal count;      /**< Number of accepted samples. */
    PetscReal weight;     /**< Total weight W. */
    PetscReal weight_sq;  /**< Sum of squared weights W2. */
    PetscReal mean;       /**< Weighted mean mu. */
    PetscReal m2;         /**< Centered second-moment sum M2. */
} PicurvMomentState;

/**
 * @brief Weighted centered co-moment state for one ordered pair of quantities.
 *
 * The two means are tracked inside the pair rather than referenced from two
 * `PicurvMomentState` values, so a co-moment update is self-contained and cannot
 * be corrupted by update ordering between the members.
 */
typedef struct {
    PetscReal count;      /**< Number of accepted samples. */
    PetscReal weight;     /**< Total weight W. */
    PetscReal weight_sq;  /**< Sum of squared weights W2. */
    PetscReal mean_x;     /**< Weighted mean of the first member. */
    PetscReal mean_y;     /**< Weighted mean of the second member. */
    PetscReal cm;         /**< Centered co-moment sum C. */
} PicurvCoMomentState;

/**
 * @brief Resets a scalar moment accumulator to the empty state.
 * @param[out] state Accumulator to clear; ignored when NULL.
 */
void PicurvMomentStateReset(PicurvMomentState *state);

/**
 * @brief Resets a co-moment accumulator to the empty state.
 * @param[out] state Accumulator to clear; ignored when NULL.
 */
void PicurvCoMomentStateReset(PicurvCoMomentState *state);

/**
 * @brief Applies one weighted sample to a scalar moment accumulator.
 *
 * Implements the stable weighted update
 * `W' = W+w`, `d = x-mu`, `mu' = mu + (w/W')d`, `M2' = M2 + w*d*(x-mu')`.
 *
 * @param[in,out] state  Accumulator to update.
 * @param[in]     value  Sample value x.
 * @param[in]     weight Sample weight w; must be strictly positive.
 * @return Zero on success, or `PETSC_ERR_ARG_OUTOFRANGE` for a non-positive weight.
 */
PetscErrorCode PicurvMomentStateUpdate(PicurvMomentState *state, PetscReal value, PetscReal weight);

/**
 * @brief Applies one weighted paired sample to a co-moment accumulator.
 *
 * Uses the pre-update mean of x and the post-update mean of y, which is the form
 * that reduces exactly to the `M2` update when x and y are the same signal.
 *
 * @param[in,out] state   Accumulator to update.
 * @param[in]     value_x Sample value of the first member.
 * @param[in]     value_y Sample value of the second member.
 * @param[in]     weight  Sample weight w; must be strictly positive.
 * @return Zero on success, or `PETSC_ERR_ARG_OUTOFRANGE` for a non-positive weight.
 */
PetscErrorCode PicurvCoMomentStateUpdate(PicurvCoMomentState *state,
                                         PetscReal value_x, PetscReal value_y, PetscReal weight);

/**
 * @brief Merges two independently accumulated scalar states.
 *
 * Implements the stable weighted parallel combination
 * `d = mu_b - mu_a`, `W = W_a + W_b`, `mu = mu_a + (W_b/W)d`,
 * `M2 = M2_a + M2_b + (W_a*W_b/W)d^2`.
 * Merging an empty state is a no-op, so partitions that received no samples are
 * safe to combine.
 *
 * @param[out] result Merged state; may alias @p a or @p b.
 * @param[in]  a      First partition.
 * @param[in]  b      Second partition.
 * @return Zero on success, or `PETSC_ERR_ARG_NULL` when any argument is NULL.
 */
PetscErrorCode PicurvMomentStateMerge(PicurvMomentState *result,
                                      const PicurvMomentState *a, const PicurvMomentState *b);

/**
 * @brief Merges two independently accumulated co-moment states.
 * @param[out] result Merged state; may alias @p a or @p b.
 * @param[in]  a      First partition.
 * @param[in]  b      Second partition.
 * @return Zero on success, or `PETSC_ERR_ARG_NULL` when any argument is NULL.
 */
PetscErrorCode PicurvCoMomentStateMerge(PicurvCoMomentState *result,
                                        const PicurvCoMomentState *a, const PicurvCoMomentState *b);

/**
 * @brief Returns the weighted variance `M2/W`, or zero when no weight accumulated.
 *
 * This is the population (weight-normalized) variance the postprocessing contract
 * uses for `R_ii`. It is never negative for a state built only through
 * `PicurvMomentStateUpdate`, but callers taking a square root must still clamp,
 * because cancellation can drive `M2` slightly negative for degenerate inputs.
 *
 * @param[in] state Accumulator to read.
 * @return Weighted variance, or zero for a NULL or unsampled accumulator.
 */
PetscReal PicurvMomentStateVariance(const PicurvMomentState *state);

/**
 * @brief Returns the weighted covariance `C/W`, or zero when no weight accumulated.
 * @param[in] state Accumulator to read.
 * @return Weighted covariance, or zero for a NULL or unsampled accumulator.
 */
PetscReal PicurvCoMomentStateCovariance(const PicurvCoMomentState *state);

/**
 * @brief Returns Kish effective sample size `W^2/W2`, or zero when no weight accumulated.
 *
 * Equals the sample count exactly under equal weights, and degrades toward one as
 * the weight distribution becomes dominated by a single sample.
 *
 * @param[in] state Accumulator to read.
 * @return Effective sample size, or zero for a NULL or unsampled accumulator.
 */
PetscReal PicurvMomentStateEffectiveCount(const PicurvMomentState *state);

#endif /* PICURV_STATISTICS_MOMENTS_H */
