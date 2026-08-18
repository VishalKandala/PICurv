/**
 * @file statistics_moments.c
 * @brief Weighted centered-moment kernels for the field-statistics pipeline.
 *
 * Full API contract is documented with the declarations in
 * `include/statistics_moments.h`.
 */

#include "statistics_moments.h"

/**
 * @brief Implementation of \ref PicurvMomentStateReset().
 * @see PicurvMomentStateReset()
 */
void PicurvMomentStateReset(PicurvMomentState *state)
{
    if (state == NULL) return;
    state->count = 0.0;
    state->weight = 0.0;
    state->weight_sq = 0.0;
    state->mean = 0.0;
    state->m2 = 0.0;
}

/**
 * @brief Implementation of \ref PicurvCoMomentStateReset().
 * @see PicurvCoMomentStateReset()
 */
void PicurvCoMomentStateReset(PicurvCoMomentState *state)
{
    if (state == NULL) return;
    state->count = 0.0;
    state->weight = 0.0;
    state->weight_sq = 0.0;
    state->mean_x = 0.0;
    state->mean_y = 0.0;
    state->cm = 0.0;
}

/**
 * @brief Implementation of \ref PicurvMomentStateUpdate().
 * @see PicurvMomentStateUpdate()
 */
PetscErrorCode PicurvMomentStateUpdate(PicurvMomentState *state, PetscReal value, PetscReal weight)
{
    PetscReal new_weight = 0.0;
    PetscReal delta = 0.0;

    PetscFunctionBeginUser;
    PetscCheck(state != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Moment accumulator is required.");
    PetscCheck(weight > 0.0, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Sample weight must be strictly positive, got %g.", (double)weight);

    new_weight = state->weight + weight;
    delta = value - state->mean;
    state->mean += (weight / new_weight) * delta;
    /* Uses the post-update mean deliberately: this is the pairing that keeps the
     * centered sum stable when the mean is large relative to the fluctuation. */
    state->m2 += weight * delta * (value - state->mean);
    state->weight = new_weight;
    state->weight_sq += weight * weight;
    state->count += 1.0;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvCoMomentStateUpdate().
 * @see PicurvCoMomentStateUpdate()
 */
PetscErrorCode PicurvCoMomentStateUpdate(PicurvCoMomentState *state,
                                         PetscReal value_x, PetscReal value_y, PetscReal weight)
{
    PetscReal new_weight = 0.0;
    PetscReal delta_x = 0.0;

    PetscFunctionBeginUser;
    PetscCheck(state != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Co-moment accumulator is required.");
    PetscCheck(weight > 0.0, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Sample weight must be strictly positive, got %g.", (double)weight);

    new_weight = state->weight + weight;
    delta_x = value_x - state->mean_x;
    state->mean_x += (weight / new_weight) * delta_x;
    state->mean_y += (weight / new_weight) * (value_y - state->mean_y);
    /* Pre-update mean of x against post-update mean of y. With x == y this
     * reduces exactly to the M2 update, which the unit tests assert. */
    state->cm += weight * delta_x * (value_y - state->mean_y);
    state->weight = new_weight;
    state->weight_sq += weight * weight;
    state->count += 1.0;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvMomentStateMerge().
 * @see PicurvMomentStateMerge()
 */
PetscErrorCode PicurvMomentStateMerge(PicurvMomentState *result,
                                      const PicurvMomentState *a, const PicurvMomentState *b)
{
    PetscReal total_weight = 0.0;
    PetscReal delta = 0.0;
    PicurvMomentState merged;

    PetscFunctionBeginUser;
    PetscCheck(result != NULL && a != NULL && b != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Merge requires a destination and two source accumulators.");

    if (b->weight == 0.0) { *result = *a; PetscFunctionReturn(0); }
    if (a->weight == 0.0) { *result = *b; PetscFunctionReturn(0); }

    total_weight = a->weight + b->weight;
    delta = b->mean - a->mean;
    merged.count = a->count + b->count;
    merged.weight = total_weight;
    merged.weight_sq = a->weight_sq + b->weight_sq;
    merged.mean = a->mean + (b->weight / total_weight) * delta;
    merged.m2 = a->m2 + b->m2 + (a->weight * b->weight / total_weight) * delta * delta;
    *result = merged;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvCoMomentStateMerge().
 * @see PicurvCoMomentStateMerge()
 */
PetscErrorCode PicurvCoMomentStateMerge(PicurvCoMomentState *result,
                                        const PicurvCoMomentState *a, const PicurvCoMomentState *b)
{
    PetscReal total_weight = 0.0;
    PetscReal delta_x = 0.0;
    PetscReal delta_y = 0.0;
    PicurvCoMomentState merged;

    PetscFunctionBeginUser;
    PetscCheck(result != NULL && a != NULL && b != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Merge requires a destination and two source accumulators.");

    if (b->weight == 0.0) { *result = *a; PetscFunctionReturn(0); }
    if (a->weight == 0.0) { *result = *b; PetscFunctionReturn(0); }

    total_weight = a->weight + b->weight;
    delta_x = b->mean_x - a->mean_x;
    delta_y = b->mean_y - a->mean_y;
    merged.count = a->count + b->count;
    merged.weight = total_weight;
    merged.weight_sq = a->weight_sq + b->weight_sq;
    merged.mean_x = a->mean_x + (b->weight / total_weight) * delta_x;
    merged.mean_y = a->mean_y + (b->weight / total_weight) * delta_y;
    merged.cm = a->cm + b->cm + (a->weight * b->weight / total_weight) * delta_x * delta_y;
    *result = merged;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PicurvMomentStateVariance().
 * @see PicurvMomentStateVariance()
 */
PetscReal PicurvMomentStateVariance(const PicurvMomentState *state)
{
    if (state == NULL || state->weight == 0.0) return 0.0;
    return state->m2 / state->weight;
}

/**
 * @brief Implementation of \ref PicurvCoMomentStateCovariance().
 * @see PicurvCoMomentStateCovariance()
 */
PetscReal PicurvCoMomentStateCovariance(const PicurvCoMomentState *state)
{
    if (state == NULL || state->weight == 0.0) return 0.0;
    return state->cm / state->weight;
}

/**
 * @brief Implementation of \ref PicurvMomentStateEffectiveCount().
 * @see PicurvMomentStateEffectiveCount()
 */
PetscReal PicurvMomentStateEffectiveCount(const PicurvMomentState *state)
{
    if (state == NULL || state->weight_sq == 0.0) return 0.0;
    return (state->weight * state->weight) / state->weight_sq;
}
