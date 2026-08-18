/**
 * @file test_statistics_moments.c
 * @brief C unit tests for the weighted centered-moment kernels.
 *
 * Covers the Stage 1 acceptance items in
 * @ref 60_Field_Statistics_Phase2_Implementation_Specification section 13:
 * constant fields yielding exactly zero covariance, known two- and three-sample
 * moments across all six symmetric velocity components, high-mean/low-fluctuation
 * precision, and merge-equals-sequential.
 */

#include "test_support.h"

#include "statistics_moments.h"

/** @brief A constant signal must produce exactly zero variance, not merely small. */
static PetscErrorCode TestConstantFieldHasExactlyZeroVariance(void)
{
    PicurvMomentState scalar;
    PicurvCoMomentState pair;

    PetscFunctionBeginUser;
    PicurvMomentStateReset(&scalar);
    PicurvCoMomentStateReset(&pair);

    /* Unequal weights so this also covers the physical-time weighting path. */
    for (PetscInt i = 0; i < 8; ++i) {
        PetscCall(PicurvMomentStateUpdate(&scalar, 7.25, 1.0 + 0.5 * (PetscReal)i));
        PetscCall(PicurvCoMomentStateUpdate(&pair, 7.25, -3.5, 1.0 + 0.5 * (PetscReal)i));
    }

    PetscCall(PicurvAssertBool((PetscBool)(scalar.m2 == 0.0),
                               "constant signal must give bitwise-zero centered second moment"));
    PetscCall(PicurvAssertBool((PetscBool)(pair.cm == 0.0),
                               "constant signal pair must give bitwise-zero co-moment"));
    PetscCall(PicurvAssertRealNear(7.25, scalar.mean, 1.0e-14, "constant signal mean"));
    PetscCall(PicurvAssertBool((PetscBool)(PicurvMomentStateVariance(&scalar) == 0.0),
                               "constant signal variance must be exactly zero"));
    PetscFunctionReturn(0);
}

/** @brief Known scalar moments for equal and unequal weights. */
static PetscErrorCode TestKnownScalarMoments(void)
{
    PicurvMomentState equal;
    PicurvMomentState weighted;
    const PetscReal samples[3] = {1.0, 2.0, 6.0};

    PetscFunctionBeginUser;
    PicurvMomentStateReset(&equal);
    for (PetscInt i = 0; i < 3; ++i) PetscCall(PicurvMomentStateUpdate(&equal, samples[i], 1.0));
    /* mean 3, M2 = 4 + 1 + 9 = 14, variance 14/3 */
    PetscCall(PicurvAssertRealNear(3.0, equal.mean, 1.0e-14, "three-sample mean"));
    PetscCall(PicurvAssertRealNear(14.0, equal.m2, 1.0e-13, "three-sample centered second moment"));
    PetscCall(PicurvAssertRealNear(14.0 / 3.0, PicurvMomentStateVariance(&equal), 1.0e-14,
                                   "three-sample variance"));
    PetscCall(PicurvAssertRealNear(3.0, PicurvMomentStateEffectiveCount(&equal), 1.0e-14,
                                   "equal weights make effective count equal sample count"));

    /* Unequal weights: x=1 (w=1), x=3 (w=3) -> W=4, mean 2.5, M2 = 3, variance 0.75 */
    PicurvMomentStateReset(&weighted);
    PetscCall(PicurvMomentStateUpdate(&weighted, 1.0, 1.0));
    PetscCall(PicurvMomentStateUpdate(&weighted, 3.0, 3.0));
    PetscCall(PicurvAssertRealNear(4.0, weighted.weight, 1.0e-14, "unequal-weight total weight"));
    PetscCall(PicurvAssertRealNear(2.5, weighted.mean, 1.0e-14, "unequal-weight mean"));
    PetscCall(PicurvAssertRealNear(3.0, weighted.m2, 1.0e-13, "unequal-weight centered second moment"));
    PetscCall(PicurvAssertRealNear(0.75, PicurvMomentStateVariance(&weighted), 1.0e-14,
                                   "unequal-weight variance"));
    /* Kish effective count: W^2/W2 = 16/10 */
    PetscCall(PicurvAssertRealNear(1.6, PicurvMomentStateEffectiveCount(&weighted), 1.0e-14,
                                   "unequal weights reduce effective count"));
    PetscFunctionReturn(0);
}

/** @brief Known two-sample covariance through the co-moment update. */
static PetscErrorCode TestKnownTwoSampleCovariance(void)
{
    PicurvCoMomentState pair;

    PetscFunctionBeginUser;
    PicurvCoMomentStateReset(&pair);
    PetscCall(PicurvCoMomentStateUpdate(&pair, 1.0, 2.0, 1.0));
    PetscCall(PicurvCoMomentStateUpdate(&pair, 3.0, 6.0, 1.0));
    /* means (2,4); deviations (-1,-2) and (1,2); C = 2 + 2 = 4; covariance 2 */
    PetscCall(PicurvAssertRealNear(2.0, pair.mean_x, 1.0e-14, "two-sample co-moment mean x"));
    PetscCall(PicurvAssertRealNear(4.0, pair.mean_y, 1.0e-14, "two-sample co-moment mean y"));
    PetscCall(PicurvAssertRealNear(4.0, pair.cm, 1.0e-13, "two-sample centered co-moment"));
    PetscCall(PicurvAssertRealNear(2.0, PicurvCoMomentStateCovariance(&pair), 1.0e-14,
                                   "two-sample covariance"));
    PetscFunctionReturn(0);
}

/**
 * @brief All six symmetric components of a three-sample vector self-product.
 *
 * Component order is the fixed upper-triangular row-major order required by the
 * accumulator contract: (xx, xy, xz, yy, yz, zz).
 */
static PetscErrorCode TestSixSymmetricVelocityComponents(void)
{
    const PetscReal series[3][3] = {{1.0, 2.0, 3.0}, {3.0, 6.0, 5.0}, {5.0, 4.0, 7.0}};
    const PetscInt first[6]  = {0, 0, 0, 1, 1, 2};
    const PetscInt second[6] = {0, 1, 2, 1, 2, 2};
    const PetscReal expected_cm[6] = {8.0, 4.0, 8.0, 8.0, 4.0, 8.0};
    const char *labels[6] = {"xx", "xy", "xz", "yy", "yz", "zz"};
    PicurvCoMomentState products[6];
    PetscReal tke = 0.0;
    char context[128];

    PetscFunctionBeginUser;
    for (PetscInt c = 0; c < 6; ++c) PicurvCoMomentStateReset(&products[c]);

    for (PetscInt s = 0; s < 3; ++s) {
        for (PetscInt c = 0; c < 6; ++c) {
            PetscCall(PicurvCoMomentStateUpdate(&products[c],
                                                series[s][first[c]], series[s][second[c]], 1.0));
        }
    }

    /* means (3,4,5); C_xx=C_xz=C_yy=C_zz=8, C_xy=C_yz=4 */
    for (PetscInt c = 0; c < 6; ++c) {
        PetscCall(PetscSNPrintf(context, sizeof(context), "centered co-moment component %s", labels[c]));
        PetscCall(PicurvAssertRealNear(expected_cm[c], products[c].cm, 1.0e-13, context));
    }

    /* TKE = 0.5 * (R_xx + R_yy + R_zz) with R_ii = C_ii / W = 8/3 each */
    tke = 0.5 * (PicurvCoMomentStateCovariance(&products[0]) +
                 PicurvCoMomentStateCovariance(&products[3]) +
                 PicurvCoMomentStateCovariance(&products[5]));
    PetscCall(PicurvAssertRealNear(4.0, tke, 1.0e-13, "turbulent kinetic energy from the trace"));
    PetscFunctionReturn(0);
}

/** @brief A co-moment of a signal with itself must reproduce the scalar second moment bitwise. */
static PetscErrorCode TestCoMomentOfSelfMatchesSecondMoment(void)
{
    PicurvMomentState scalar;
    PicurvCoMomentState self;
    const PetscReal samples[5] = {2.5, -1.75, 9.0, 0.25, 4.5};
    const PetscReal weights[5] = {1.0, 0.5, 2.25, 3.0, 0.125};

    PetscFunctionBeginUser;
    PicurvMomentStateReset(&scalar);
    PicurvCoMomentStateReset(&self);
    for (PetscInt i = 0; i < 5; ++i) {
        PetscCall(PicurvMomentStateUpdate(&scalar, samples[i], weights[i]));
        PetscCall(PicurvCoMomentStateUpdate(&self, samples[i], samples[i], weights[i]));
    }
    PetscCall(PicurvAssertBool((PetscBool)(self.cm == scalar.m2),
                               "co-moment of a signal with itself must equal its centered second moment"));
    PetscCall(PicurvAssertBool((PetscBool)(self.mean_x == scalar.mean && self.mean_y == scalar.mean),
                               "co-moment self-pair means must equal the scalar mean"));
    PetscFunctionReturn(0);
}

/**
 * @brief High mean with small fluctuation, where a naive sum-of-squares cancels catastrophically.
 *
 * Samples 1e8-1 and 1e8+1 have variance exactly 1. Accumulating raw second sums
 * would compute a difference of two values near 2e16, whose spacing exceeds the
 * answer, so the centered update is what makes this recoverable.
 */
static PetscErrorCode TestHighMeanLowFluctuationPrecision(void)
{
    PicurvMomentState state;
    PicurvCoMomentState pair;

    PetscFunctionBeginUser;
    PicurvMomentStateReset(&state);
    PicurvCoMomentStateReset(&pair);
    PetscCall(PicurvMomentStateUpdate(&state, 1.0e8 - 1.0, 1.0));
    PetscCall(PicurvMomentStateUpdate(&state, 1.0e8 + 1.0, 1.0));
    PetscCall(PicurvCoMomentStateUpdate(&pair, 1.0e8 - 1.0, 1.0e8 + 1.0, 1.0));
    PetscCall(PicurvCoMomentStateUpdate(&pair, 1.0e8 + 1.0, 1.0e8 - 1.0, 1.0));

    PetscCall(PicurvAssertRealNear(1.0e8, state.mean, 1.0e-6, "high-mean signal mean"));
    PetscCall(PicurvAssertRealNear(1.0, PicurvMomentStateVariance(&state), 1.0e-9,
                                   "high-mean low-fluctuation variance must not cancel"));
    /* Anti-correlated pair of the same fluctuation: covariance is exactly -1. */
    PetscCall(PicurvAssertRealNear(-1.0, PicurvCoMomentStateCovariance(&pair), 1.0e-9,
                                   "high-mean anti-correlated covariance must not cancel"));
    PetscFunctionReturn(0);
}

/** @brief Merging two partitions must reproduce a single sequential accumulation. */
static PetscErrorCode TestMergeEqualsSequential(void)
{
    PicurvMomentState sequential, part_a, part_b, merged;
    PicurvCoMomentState co_sequential, co_a, co_b, co_merged;
    PicurvMomentState empty;
    const PetscReal samples[8] = {1.5, -2.0, 7.25, 0.5, 3.0, -4.5, 6.0, 2.25};
    const PetscReal weights[8] = {1.0, 2.0, 0.5, 1.25, 3.0, 0.75, 1.0, 2.5};

    PetscFunctionBeginUser;
    PicurvMomentStateReset(&sequential);
    PicurvMomentStateReset(&part_a);
    PicurvMomentStateReset(&part_b);
    PicurvMomentStateReset(&empty);
    PicurvCoMomentStateReset(&co_sequential);
    PicurvCoMomentStateReset(&co_a);
    PicurvCoMomentStateReset(&co_b);

    for (PetscInt i = 0; i < 8; ++i) {
        PetscCall(PicurvMomentStateUpdate(&sequential, samples[i], weights[i]));
        PetscCall(PicurvCoMomentStateUpdate(&co_sequential, samples[i], 2.0 * samples[i] + 1.0, weights[i]));
        if (i < 3) {
            PetscCall(PicurvMomentStateUpdate(&part_a, samples[i], weights[i]));
            PetscCall(PicurvCoMomentStateUpdate(&co_a, samples[i], 2.0 * samples[i] + 1.0, weights[i]));
        } else {
            PetscCall(PicurvMomentStateUpdate(&part_b, samples[i], weights[i]));
            PetscCall(PicurvCoMomentStateUpdate(&co_b, samples[i], 2.0 * samples[i] + 1.0, weights[i]));
        }
    }

    PetscCall(PicurvMomentStateMerge(&merged, &part_a, &part_b));
    PetscCall(PicurvAssertRealNear(sequential.weight, merged.weight, 1.0e-14, "merged total weight"));
    PetscCall(PicurvAssertRealNear(sequential.count, merged.count, 1.0e-14, "merged sample count"));
    PetscCall(PicurvAssertRealNear(sequential.mean, merged.mean, 1.0e-13, "merged mean"));
    PetscCall(PicurvAssertRealNear(sequential.m2, merged.m2, 1.0e-11, "merged centered second moment"));

    PetscCall(PicurvCoMomentStateMerge(&co_merged, &co_a, &co_b));
    PetscCall(PicurvAssertRealNear(co_sequential.mean_x, co_merged.mean_x, 1.0e-13, "merged co-moment mean x"));
    PetscCall(PicurvAssertRealNear(co_sequential.mean_y, co_merged.mean_y, 1.0e-13, "merged co-moment mean y"));
    PetscCall(PicurvAssertRealNear(co_sequential.cm, co_merged.cm, 1.0e-11, "merged centered co-moment"));

    /* Merging an unsampled partition must be a no-op in both directions. */
    PetscCall(PicurvMomentStateMerge(&merged, &sequential, &empty));
    PetscCall(PicurvAssertBool((PetscBool)(merged.m2 == sequential.m2 && merged.mean == sequential.mean),
                               "merging an empty partition on the right must not change the state"));
    PetscCall(PicurvMomentStateMerge(&merged, &empty, &sequential));
    PetscCall(PicurvAssertBool((PetscBool)(merged.m2 == sequential.m2 && merged.mean == sequential.mean),
                               "merging an empty partition on the left must not change the state"));
    PetscFunctionReturn(0);
}

/** @brief Non-positive weights are rejected rather than silently corrupting the accumulator. */
static PetscErrorCode TestNonPositiveWeightRejected(void)
{
    PicurvMomentState state;
    PicurvCoMomentState pair;
    PetscErrorCode zero_ierr = 0, negative_ierr = 0, co_zero_ierr = 0;

    PetscFunctionBeginUser;
    PicurvMomentStateReset(&state);
    PicurvCoMomentStateReset(&pair);

    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    zero_ierr = PicurvMomentStateUpdate(&state, 1.0, 0.0);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_OUTOFRANGE, zero_ierr,
                                   "zero sample weight should be rejected"));

    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    negative_ierr = PicurvMomentStateUpdate(&state, 1.0, -2.0);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_OUTOFRANGE, negative_ierr,
                                   "negative sample weight should be rejected"));

    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    co_zero_ierr = PicurvCoMomentStateUpdate(&pair, 1.0, 2.0, 0.0);
    PetscCall(PetscPopErrorHandler());
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_OUTOFRANGE, co_zero_ierr,
                                   "zero co-moment weight should be rejected"));

    PetscCall(PicurvAssertBool((PetscBool)(state.count == 0.0 && state.weight == 0.0),
                               "a rejected sample must leave the accumulator untouched"));
    PetscFunctionReturn(0);
}

/**
 * @brief Entry point for the centered-moment kernel suite.
 */
int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"constant-field-zero-variance", TestConstantFieldHasExactlyZeroVariance},
        {"known-scalar-moments", TestKnownScalarMoments},
        {"known-two-sample-covariance", TestKnownTwoSampleCovariance},
        {"six-symmetric-velocity-components", TestSixSymmetricVelocityComponents},
        {"co-moment-of-self-matches-second-moment", TestCoMomentOfSelfMatchesSecondMoment},
        {"high-mean-low-fluctuation-precision", TestHighMeanLowFluctuationPrecision},
        {"merge-equals-sequential", TestMergeEqualsSequential},
        {"non-positive-weight-rejected", TestNonPositiveWeightRejected},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv centered-moment kernel tests");
    if (ierr) {
        return (int)ierr;
    }

    ierr = PicurvRunTests("unit-statistics", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) {
        PetscFinalize();
        return (int)ierr;
    }

    ierr = PetscFinalize();
    return (int)ierr;
}
