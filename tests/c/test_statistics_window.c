/**
 * @file test_statistics_window.c
 * @brief C unit tests for window lifecycle, scheduling, and weighting.
 *
 * Covers the Stage 3 acceptance items in
 * @ref 60_Field_Statistics_Phase2_Implementation_Specification section 14:
 * cadence stride, start and end clipping, variable timestep weighting, duplicate
 * event rejection, off-schedule no-ops, and the property that sample and
 * physical-time weighting agree on a constant-timestep run.
 */

#include "test_support.h"

#include "statistics_window.h"
#include "field_catalog.h"

/** @brief Builds a step-cadence definition. */
static PicurvWindowDefinition StepWindow(const char *name, PetscReal start, PetscReal end,
                                         PetscBool bounded, PicurvWeighting weighting,
                                         PetscInt cadence)
{
    PicurvWindowDefinition d;
    memset(&d, 0, sizeof(d));
    strncpy(d.name, name, PICURV_WINDOW_NAME_LENGTH - 1);
    d.start_time = start; d.end_time = end; d.bounded = bounded;
    d.weighting = weighting; d.cadence_kind = PICURV_CADENCE_STEP; d.step_cadence = cadence;
    return d;
}

/** @brief Drives a uniform-dt sequence and returns the accumulated weight and count. */
static PetscErrorCode RunUniform(PicurvWindow *w, PetscInt steps, PetscReal dt,
                                 PetscReal *total_weight, PetscInt *count)
{
    PetscFunctionBeginUser;
    for (PetscInt s = 0; s <= steps; ++s) {
        PetscBool accepted = PETSC_FALSE;
        PetscReal weight = 0.0;
        PetscCall(PicurvWindowOfferState(w, s, (PetscReal)s * dt, &accepted, &weight));
    }
    *total_weight = w->total_weight;
    *count = w->sample_count;
    PetscFunctionReturn(0);
}

/**
 * @brief The two weightings must agree on a constant-timestep run.
 *
 * This is the property that fixed the initial-state rule: counting the state at
 * the window origin under sample weighting alone would make the two disagree by
 * O(1/N) for no physical reason.
 */
static PetscErrorCode TestWeightingsAgreeAtConstantTimestep(void)
{
    PicurvWindow sample_w, time_w;
    PetscReal sample_weight = 0.0, time_weight = 0.0;
    PetscInt sample_count = 0, time_count = 0;
    const PetscReal dt = 0.25;
    const PetscInt steps = 20;

    PetscFunctionBeginUser;
    {
        PicurvWindowDefinition d = StepWindow("s", 0.0, 0.0, PETSC_FALSE, PICURV_WEIGHTING_SAMPLE, 1);
        PetscCall(PicurvWindowInit(&sample_w, &d));
    }
    {
        PicurvWindowDefinition d = StepWindow("t", 0.0, 0.0, PETSC_FALSE, PICURV_WEIGHTING_PHYSICAL_TIME, 1);
        PetscCall(PicurvWindowInit(&time_w, &d));
    }
    PetscCall(RunUniform(&sample_w, steps, dt, &sample_weight, &sample_count));
    PetscCall(RunUniform(&time_w, steps, dt, &time_weight, &time_count));

    /* The state at t=0 anchors the origin and is not a sample under either rule. */
    PetscCall(PicurvAssertIntEqual(steps, sample_count, "sample weighting counts one per interval"));
    PetscCall(PicurvAssertIntEqual(steps, time_count, "time weighting counts one per interval"));
    PetscCall(PicurvAssertIntEqual(sample_count, time_count,
                                   "both weightings must accept exactly the same states"));
    PetscCall(PicurvAssertRealNear((PetscReal)steps, sample_weight, 1.0e-12, "sample total weight"));
    PetscCall(PicurvAssertRealNear((PetscReal)steps * dt, time_weight, 1.0e-12, "time total weight"));
    /* A mean is weight-normalized, so equal counts with proportional weights means
     * the two weightings produce the identical mean for a uniform-dt run. */
    PetscCall(PicurvAssertRealNear(dt, time_weight / sample_weight, 1.0e-12,
                                   "weights differ only by the constant timestep"));
    PetscCall(PicurvAssertRealNear((PetscReal)steps * dt, time_w.represented_time, 1.0e-12,
                                   "represented time spans the whole run"));
    PetscFunctionReturn(0);
}

/** @brief A state at the window origin anchors without becoming a sample. */
static PetscErrorCode TestOriginStateAnchorsWithoutSampling(void)
{
    PicurvWindow w;
    PicurvWindowDefinition d = StepWindow("anchor", 5.0, 0.0, PETSC_FALSE,
                                          PICURV_WEIGHTING_PHYSICAL_TIME, 1);
    PetscBool accepted = PETSC_FALSE;
    PetscReal weight = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvWindowInit(&w, &d));
    PetscCall(PicurvAssertIntEqual(PICURV_WINDOW_PENDING, w.state, "window starts pending"));

    /* Before the start: no effect at all. */
    PetscCall(PicurvWindowOfferState(&w, 10, 4.5, &accepted, &weight));
    PetscCall(PicurvAssertBool((PetscBool)!accepted, "states before the start are ignored"));
    PetscCall(PicurvAssertIntEqual(PICURV_WINDOW_PENDING, w.state, "an early state does not activate"));

    /* Exactly at the start: activates and anchors, but is not a sample. */
    PetscCall(PicurvWindowOfferState(&w, 11, 5.0, &accepted, &weight));
    PetscCall(PicurvAssertBool((PetscBool)!accepted, "the origin state is not a sample"));
    PetscCall(PicurvAssertIntEqual(PICURV_WINDOW_ACTIVE, w.state, "the origin state activates the window"));
    PetscCall(PicurvAssertIntEqual(0, w.sample_count, "anchoring records no sample"));

    /* Next state carries the interval back to the origin. */
    PetscCall(PicurvWindowOfferState(&w, 12, 5.5, &accepted, &weight));
    PetscCall(PicurvAssertBool(accepted, "the state after the origin is a sample"));
    PetscCall(PicurvAssertRealNear(0.5, weight, 1.0e-12, "first sample carries the interval from the origin"));
    PetscFunctionReturn(0);
}

/** @brief A window first seen after its requested start moves its origin forward. */
static PetscErrorCode TestLateFirstObservationMovesOrigin(void)
{
    PicurvWindow w;
    PicurvWindowDefinition d = StepWindow("late", 5.0, 0.0, PETSC_FALSE,
                                          PICURV_WEIGHTING_PHYSICAL_TIME, 1);
    PetscBool accepted = PETSC_FALSE;
    PetscReal weight = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvWindowInit(&w, &d));
    /* Resumed at t=8 with no earlier observation: the window must not claim [5,8]. */
    PetscCall(PicurvWindowOfferState(&w, 40, 8.0, &accepted, &weight));
    PetscCall(PicurvAssertBool((PetscBool)!accepted, "the first observed state anchors the moved origin"));
    PetscCall(PicurvAssertRealNear(8.0, w.effective_start, 1.0e-12,
                                   "effective start moves to the first observed time"));
    PetscCall(PicurvWindowOfferState(&w, 41, 8.5, &accepted, &weight));
    PetscCall(PicurvAssertBool(accepted, "the next state is a sample"));
    PetscCall(PicurvAssertRealNear(0.5, weight, 1.0e-12, "weight is measured from the moved origin"));
    PetscCall(PicurvAssertRealNear(0.5, w.represented_time, 1.0e-12,
                                   "a window never claims time it did not observe"));
    PetscFunctionReturn(0);
}

/** @brief Variable timestep weighting follows the actual elapsed intervals. */
static PetscErrorCode TestVariableTimestepWeighting(void)
{
    PicurvWindow w;
    PicurvWindowDefinition d = StepWindow("vardt", 0.0, 0.0, PETSC_FALSE,
                                          PICURV_WEIGHTING_PHYSICAL_TIME, 1);
    const PetscReal times[4] = {0.0, 0.5, 2.0, 2.25};
    const PetscReal expected[4] = {0.0, 0.5, 1.5, 0.25};
    PetscBool accepted = PETSC_FALSE;
    PetscReal weight = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvWindowInit(&w, &d));
    for (PetscInt i = 0; i < 4; ++i) {
        PetscCall(PicurvWindowOfferState(&w, i, times[i], &accepted, &weight));
        PetscCall(PicurvAssertRealNear(expected[i], weight, 1.0e-12,
                                       "variable-dt weight equals the elapsed interval"));
    }
    PetscCall(PicurvAssertIntEqual(3, w.sample_count, "three intervals become samples"));
    PetscCall(PicurvAssertRealNear(2.25, w.total_weight, 1.0e-12, "weights sum to the span"));
    PetscFunctionReturn(0);
}

/** @brief Stride skips states, and the accepted weight still spans the whole gap. */
static PetscErrorCode TestStrideAndOffScheduleNoOp(void)
{
    PicurvWindow w;
    PicurvWindowDefinition d = StepWindow("stride", 0.0, 0.0, PETSC_FALSE,
                                          PICURV_WEIGHTING_PHYSICAL_TIME, 3);
    PetscBool accepted = PETSC_FALSE;
    PetscReal weight = 0.0;
    const PetscReal dt = 0.1;

    PetscFunctionBeginUser;
    PetscCall(PicurvWindowInit(&w, &d));
    for (PetscInt s = 0; s <= 6; ++s) {
        PetscReal before_weight = w.total_weight;
        PetscInt before_count = w.sample_count;
        PetscCall(PicurvWindowOfferState(&w, s, (PetscReal)s * dt, &accepted, &weight));
        if (s % 3 != 0) {
            PetscCall(PicurvAssertBool((PetscBool)!accepted, "off-schedule states are not sampled"));
            PetscCall(PicurvAssertRealNear(before_weight, w.total_weight, 1.0e-15,
                                           "an off-schedule state changes no scientific state"));
            PetscCall(PicurvAssertIntEqual(before_count, w.sample_count,
                                           "an off-schedule state records no sample"));
        }
    }
    /* Steps 3 and 6 are sampled; each represents three timesteps. */
    PetscCall(PicurvAssertIntEqual(2, w.sample_count, "stride 3 over 6 steps yields two samples"));
    PetscCall(PicurvAssertRealNear(0.6, w.total_weight, 1.0e-12,
                                   "strided weights still cover the full elapsed span"));
    PetscFunctionReturn(0);
}

/** @brief A bounded window clips its final interval and then accepts nothing. */
static PetscErrorCode TestEndClippingAndCompletion(void)
{
    PicurvWindow w;
    PicurvWindowDefinition d = StepWindow("bounded", 0.0, 1.0, PETSC_TRUE,
                                          PICURV_WEIGHTING_PHYSICAL_TIME, 1);
    PetscBool accepted = PETSC_FALSE;
    PetscReal weight = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvWindowInit(&w, &d));
    PetscCall(PicurvWindowOfferState(&w, 0, 0.0, &accepted, &weight));   /* anchor */
    PetscCall(PicurvWindowOfferState(&w, 1, 0.4, &accepted, &weight));
    PetscCall(PicurvAssertRealNear(0.4, weight, 1.0e-12, "interior interval"));
    PetscCall(PicurvWindowOfferState(&w, 2, 0.8, &accepted, &weight));
    PetscCall(PicurvAssertRealNear(0.4, weight, 1.0e-12, "interior interval"));

    /* Overshoots the bound: the final weight clips to the requested end. */
    PetscCall(PicurvWindowOfferState(&w, 3, 1.3, &accepted, &weight));
    PetscCall(PicurvAssertBool(accepted, "the overshooting state still contributes its clipped interval"));
    PetscCall(PicurvAssertRealNear(0.2, weight, 1.0e-12, "final interval clips to the requested end"));
    PetscCall(PicurvAssertIntEqual(PICURV_WINDOW_COMPLETE, w.state, "reaching the bound completes the window"));
    PetscCall(PicurvAssertRealNear(1.0, w.represented_time, 1.0e-12,
                                   "represented time equals the requested span exactly"));
    PetscCall(PicurvAssertRealNear(1.0, PicurvWindowProgress(&w), 1.0e-12, "a complete window reports full progress"));

    /* Nothing further is accepted. */
    PetscCall(PicurvWindowOfferState(&w, 4, 1.5, &accepted, &weight));
    PetscCall(PicurvAssertBool((PetscBool)!accepted, "a complete window accepts nothing further"));
    PetscCall(PicurvAssertIntEqual(3, w.sample_count, "completion does not add samples"));
    PetscFunctionReturn(0);
}

/** @brief The same completed step offered twice is counted once. */
static PetscErrorCode TestDuplicateEventRejected(void)
{
    PicurvWindow w;
    PicurvWindowDefinition d = StepWindow("dup", 0.0, 0.0, PETSC_FALSE,
                                          PICURV_WEIGHTING_PHYSICAL_TIME, 1);
    PetscBool accepted = PETSC_FALSE;
    PetscReal weight = 0.0;

    PetscFunctionBeginUser;
    PetscCall(PicurvWindowInit(&w, &d));
    PetscCall(PicurvWindowOfferState(&w, 0, 0.0, &accepted, &weight));
    PetscCall(PicurvWindowOfferState(&w, 1, 0.5, &accepted, &weight));
    PetscCall(PicurvAssertBool(accepted, "first offer of a step is accepted"));

    PetscCall(PicurvWindowOfferState(&w, 1, 0.5, &accepted, &weight));
    PetscCall(PicurvAssertBool((PetscBool)!accepted, "the same step must not be counted twice"));
    PetscCall(PicurvAssertIntEqual(1, w.sample_count, "a duplicate offer records no extra sample"));
    PetscCall(PicurvAssertRealNear(0.5, w.total_weight, 1.0e-12, "a duplicate offer adds no weight"));
    PetscFunctionReturn(0);
}

/**
 * @brief A step overshooting several time targets is accepted once, losing no time.
 *
 * This is the self-correcting property of right-rectangle weighting: the accepted
 * weight is the actual elapsed interval, not the nominal cadence, so skipped
 * targets neither drop nor double-count represented time.
 */
static PetscErrorCode TestTimeCadenceOvershootAcceptedOnce(void)
{
    PicurvWindow w;
    PicurvWindowDefinition d;
    PetscBool accepted = PETSC_FALSE;
    PetscReal weight = 0.0;

    PetscFunctionBeginUser;
    memset(&d, 0, sizeof(d));
    strncpy(d.name, "tcad", PICURV_WINDOW_NAME_LENGTH - 1);
    d.start_time = 0.0; d.bounded = PETSC_FALSE;
    d.weighting = PICURV_WEIGHTING_PHYSICAL_TIME;
    d.cadence_kind = PICURV_CADENCE_TIME; d.time_cadence = 0.1;
    PetscCall(PicurvWindowInit(&w, &d));

    PetscCall(PicurvWindowOfferState(&w, 0, 0.0, &accepted, &weight));  /* anchor */
    PetscCall(PicurvAssertBool((PetscBool)!accepted, "origin anchors under time cadence too"));

    /* One large step jumps past targets 0.1, 0.2, 0.3, 0.4 and lands on 0.45. */
    PetscCall(PicurvWindowOfferState(&w, 1, 0.45, &accepted, &weight));
    PetscCall(PicurvAssertBool(accepted, "the overshooting state is accepted"));
    PetscCall(PicurvAssertRealNear(0.45, weight, 1.0e-12,
                                   "weight is the actual elapsed interval, not the nominal cadence"));
    PetscCall(PicurvAssertIntEqual(1, w.sample_count, "several skipped targets yield exactly one sample"));

    /* The schedule resumes on the absolute grid rather than drifting from 0.45. */
    PetscCall(PicurvWindowOfferState(&w, 2, 0.46, &accepted, &weight));
    PetscCall(PicurvAssertBool((PetscBool)!accepted, "the next target has not been reached yet"));
    PetscCall(PicurvWindowOfferState(&w, 3, 0.51, &accepted, &weight));
    PetscCall(PicurvAssertBool(accepted, "the state reaching the 0.5 target is accepted"));
    PetscCall(PicurvAssertRealNear(0.06, weight, 1.0e-12, "weight resumes from the last accepted state"));
    PetscCall(PicurvAssertRealNear(0.51, w.represented_time, 1.0e-12,
                                   "no represented time is lost or double counted across the overshoot"));
    PetscFunctionReturn(0);
}

/** @brief Invalid definitions are rejected at initialization. */
static PetscErrorCode TestInvalidDefinitionsRejected(void)
{
    PicurvWindow w;
    PetscErrorCode ierr_cadence = 0, ierr_name = 0, ierr_span = 0;

    PetscFunctionBeginUser;
    {
        PicurvWindowDefinition d = StepWindow("zero", 0.0, 0.0, PETSC_FALSE, PICURV_WEIGHTING_SAMPLE, 0);
        PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
        ierr_cadence = PicurvWindowInit(&w, &d);
        PetscCall(PetscPopErrorHandler());
    }
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_OUTOFRANGE, ierr_cadence, "non-positive stride is rejected"));
    {
        PicurvWindowDefinition d = StepWindow("", 0.0, 0.0, PETSC_FALSE, PICURV_WEIGHTING_SAMPLE, 1);
        PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
        ierr_name = PicurvWindowInit(&w, &d);
        PetscCall(PetscPopErrorHandler());
    }
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_OUTOFRANGE, ierr_name, "an empty window name is rejected"));
    {
        PicurvWindowDefinition d = StepWindow("bad", 5.0, 5.0, PETSC_TRUE, PICURV_WEIGHTING_SAMPLE, 1);
        PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
        ierr_span = PicurvWindowInit(&w, &d);
        PetscCall(PetscPopErrorHandler());
    }
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_OUTOFRANGE, ierr_span, "end must exceed start"));
    PetscFunctionReturn(0);
}

/** @brief Hashes a definition and reports which property group changed against a baseline. */
static PetscErrorCode HashAndCompare(const PicurvWindowDefinition *baseline,
                                     const PicurvWindowDefinition *variant,
                                     PetscBool *same, PetscInt *first_difference)
{
    char base_digest[65], variant_digest[65];
    char base_groups[PICURV_WINDOW_HASH_GROUP_COUNT][PICURV_WINDOW_HASH_GROUP_LENGTH];
    char variant_groups[PICURV_WINDOW_HASH_GROUP_COUNT][PICURV_WINDOW_HASH_GROUP_LENGTH];

    PetscFunctionBeginUser;
    PetscCall(PicurvWindowComputeHash(baseline, base_digest, base_groups));
    PetscCall(PicurvWindowComputeHash(variant, variant_digest, variant_groups));
    PetscCall(PetscStrcmp(base_digest, variant_digest, same));
    *first_difference = -1;
    for (PetscInt group = 0; group < PICURV_WINDOW_HASH_GROUP_COUNT; ++group) {
        PetscBool group_same = PETSC_FALSE;

        PetscCall(PetscStrcmp(base_groups[group], variant_groups[group], &group_same));
        if (!group_same) { *first_difference = group; break; }
    }
    PetscFunctionReturn(0);
}

/** @brief Builds the hash fixture: a Ucat/P window with a second moment and a covariance. */
static PicurvWindowDefinition HashWindow(void)
{
    PicurvWindowDefinition d = StepWindow("production", 10.0, 20.0, PETSC_TRUE,
                                          PICURV_WEIGHTING_PHYSICAL_TIME, 5);
    d.field_count = 2;
    d.fields[0].field_id = FIELD_ID_UCAT; d.fields[0].want_second = PETSC_TRUE;
    d.fields[1].field_id = FIELD_ID_P;    d.fields[1].want_second = PETSC_FALSE;
    d.covariance_count = 1;
    d.covariances[0].first = FIELD_ID_UCAT;
    d.covariances[0].second = FIELD_ID_P;
    return d;
}

/** @brief The hash is stable, and excludes exactly the properties the spec excludes. */
static PetscErrorCode TestHashExclusionsAndStability(void)
{
    const PicurvWindowDefinition baseline = HashWindow();
    PicurvWindowDefinition variant;
    PetscBool same = PETSC_FALSE;
    PetscInt group = -1;

    PetscFunctionBeginUser;
    /* Stability: hashing the same definition twice gives the same digest. */
    variant = baseline;
    PetscCall(HashAndCompare(&baseline, &variant, &same, &group));
    PetscCall(PicurvAssertBool(same, "hashing an unchanged definition is stable"));

    /* end_time is excluded so a bounded window can be extended forward. */
    variant = baseline;
    variant.end_time = 40.0;
    PetscCall(HashAndCompare(&baseline, &variant, &same, &group));
    PetscCall(PicurvAssertBool(same, "extending end_time does not change the identity"));

    /* Listing order is excluded so a reordered configuration still continues. */
    variant = baseline;
    variant.fields[0] = baseline.fields[1];
    variant.fields[1] = baseline.fields[0];
    PetscCall(HashAndCompare(&baseline, &variant, &same, &group));
    PetscCall(PicurvAssertBool(same, "reordering the field list does not change the identity"));

    /* A covariance pair is symmetric, so its member order is excluded too. */
    variant = baseline;
    variant.covariances[0].first = FIELD_ID_P;
    variant.covariances[0].second = FIELD_ID_UCAT;
    PetscCall(HashAndCompare(&baseline, &variant, &same, &group));
    PetscCall(PicurvAssertBool(same, "swapping covariance members does not change the identity"));
    PetscFunctionReturn(0);
}

/** @brief Every hashed property changes the digest and is named by its group digest. */
static PetscErrorCode TestHashDetectsEachProperty(void)
{
    const PicurvWindowDefinition baseline = HashWindow();
    PicurvWindowDefinition variant;
    PetscBool same = PETSC_FALSE;
    PetscInt group = -1;
    PetscBool named = PETSC_FALSE;

    PetscFunctionBeginUser;
    variant = baseline;
    strncpy(variant.name, "other", PICURV_WINDOW_NAME_LENGTH - 1);
    PetscCall(HashAndCompare(&baseline, &variant, &same, &group));
    PetscCall(PicurvAssertBool((PetscBool)!same, "a renamed window is a different window"));
    PetscCall(PicurvAssertIntEqual(0, group, "the name group reports the difference"));

    variant = baseline;
    variant.start_time = 11.0;
    PetscCall(HashAndCompare(&baseline, &variant, &same, &group));
    PetscCall(PicurvAssertBool((PetscBool)!same, "a moved start is a different window"));
    PetscCall(PicurvAssertIntEqual(1, group, "the start_time group reports the difference"));

    variant = baseline;
    variant.weighting = PICURV_WEIGHTING_SAMPLE;
    PetscCall(HashAndCompare(&baseline, &variant, &same, &group));
    PetscCall(PicurvAssertBool((PetscBool)!same, "a changed weighting is a different window"));
    PetscCall(PicurvAssertIntEqual(2, group, "the weighting group reports the difference"));

    variant = baseline;
    variant.step_cadence = 7;
    PetscCall(HashAndCompare(&baseline, &variant, &same, &group));
    PetscCall(PicurvAssertBool((PetscBool)!same, "a changed cadence is a different window"));
    PetscCall(PicurvAssertIntEqual(3, group, "the cadence group reports the difference"));

    /* Dropping a requested moment changes what the saved state means. */
    variant = baseline;
    variant.fields[0].want_second = PETSC_FALSE;
    PetscCall(HashAndCompare(&baseline, &variant, &same, &group));
    PetscCall(PicurvAssertBool((PetscBool)!same, "a changed moment set is a different window"));
    PetscCall(PicurvAssertIntEqual(4, group, "the fields group reports the difference"));

    variant = baseline;
    variant.covariance_count = 0;
    PetscCall(HashAndCompare(&baseline, &variant, &same, &group));
    PetscCall(PicurvAssertBool((PetscBool)!same, "a dropped covariance is a different window"));
    PetscCall(PicurvAssertIntEqual(5, group, "the covariances group reports the difference"));

    /* Group names exist for every group and are bounded. */
    for (PetscInt g = 0; g < PICURV_WINDOW_HASH_GROUP_COUNT; ++g) {
        PetscCall(PetscStrcmp(PicurvWindowHashGroupName(g), "unknown", &named));
        PetscCall(PicurvAssertBool((PetscBool)!named, "every hashed group has a stable name"));
    }
    PetscCall(PetscStrcmp(PicurvWindowHashGroupName(PICURV_WINDOW_HASH_GROUP_COUNT), "unknown", &named));
    PetscCall(PicurvAssertBool(named, "an out-of-range group index is named unknown"));
    PetscFunctionReturn(0);
}

/** @brief Serializes a definition's group digests the way a checkpoint records them. */
static PetscErrorCode SerializeGroupDigests(const PicurvWindowDefinition *definition,
                                            char *out, size_t out_size)
{
    char digest[65];
    char groups[PICURV_WINDOW_HASH_GROUP_COUNT][PICURV_WINDOW_HASH_GROUP_LENGTH];
    size_t used = 0;

    PetscFunctionBeginUser;
    out[0] = '\0';
    PetscCall(PicurvWindowComputeHash(definition, digest, groups));
    for (PetscInt group = 0; group < PICURV_WINDOW_HASH_GROUP_COUNT; ++group) {
        PetscCall(PetscStrlen(out, &used));
        PetscCall(PetscSNPrintf(out + used, out_size - used, "%s%s", group ? "," : "", groups[group]));
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Saved group digests must localize a change to the property that caused it.
 *
 * This is the only path that turns a refused continuation into an actionable
 * message, and it is reachable at runtime only through a fatal error, so it is
 * exercised directly here rather than through that path.
 */
static PetscErrorCode TestFirstHashDifferenceLocalization(void)
{
    const PicurvWindowDefinition baseline = HashWindow();
    PicurvWindowDefinition variant;
    char saved[PICURV_WINDOW_HASH_GROUP_COUNT * (PICURV_WINDOW_HASH_GROUP_LENGTH + 4)];
    PetscInt group = 0;

    PetscFunctionBeginUser;
    PetscCall(SerializeGroupDigests(&baseline, saved, sizeof(saved)));

    /* An unchanged definition reports no differing group. */
    PetscCall(PicurvWindowFirstHashDifference(&baseline, saved, &group));
    PetscCall(PicurvAssertIntEqual(-1, group, "an unchanged definition localizes to no group"));

    /* Each changed property is localized to its own group. */
    variant = baseline;
    variant.weighting = PICURV_WEIGHTING_SAMPLE;
    PetscCall(PicurvWindowFirstHashDifference(&variant, saved, &group));
    PetscCall(PicurvAssertIntEqual(2, group, "a changed weighting localizes to the weighting group"));

    variant = baseline;
    variant.covariance_count = 0;
    PetscCall(PicurvWindowFirstHashDifference(&variant, saved, &group));
    PetscCall(PicurvAssertIntEqual(5, group, "a dropped covariance localizes to the covariances group"));

    /* The earliest differing group wins when several changed at once. */
    variant = baseline;
    variant.start_time = 99.0;
    variant.step_cadence = 3;
    PetscCall(PicurvWindowFirstHashDifference(&variant, saved, &group));
    PetscCall(PicurvAssertIntEqual(1, group, "the first differing group is the one reported"));

    /* Degenerate inputs report no localization rather than a wrong property. */
    PetscCall(PicurvWindowFirstHashDifference(&baseline, NULL, &group));
    PetscCall(PicurvAssertIntEqual(-1, group, "absent group digests localize to nothing"));
    PetscCall(PicurvWindowFirstHashDifference(&baseline, "", &group));
    PetscCall(PicurvAssertIntEqual(-1, group, "an empty digest list localizes to nothing"));

    variant = baseline;
    variant.covariance_count = 0;
    PetscCall(SerializeGroupDigests(&baseline, saved, sizeof(saved)));
    saved[40] = '\0';   /* truncate mid-digest, so the last segment is short */
    PetscCall(PicurvWindowFirstHashDifference(&variant, saved, &group));
    PetscCall(PicurvAssertIntEqual(-1, group,
                                   "a truncated digest list reports no property rather than a wrong one"));
    PetscFunctionReturn(0);
}

/**
 * @brief Entry point for the window lifecycle suite.
 */
int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"weightings-agree-at-constant-timestep", TestWeightingsAgreeAtConstantTimestep},
        {"origin-state-anchors-without-sampling", TestOriginStateAnchorsWithoutSampling},
        {"late-first-observation-moves-origin", TestLateFirstObservationMovesOrigin},
        {"variable-timestep-weighting", TestVariableTimestepWeighting},
        {"stride-and-off-schedule-no-op", TestStrideAndOffScheduleNoOp},
        {"end-clipping-and-completion", TestEndClippingAndCompletion},
        {"duplicate-event-rejected", TestDuplicateEventRejected},
        {"time-cadence-overshoot-accepted-once", TestTimeCadenceOvershootAcceptedOnce},
        {"invalid-definitions-rejected", TestInvalidDefinitionsRejected},
        {"hash-exclusions-and-stability", TestHashExclusionsAndStability},
        {"hash-detects-each-property", TestHashDetectsEachProperty},
        {"first-hash-difference-localization", TestFirstHashDifferenceLocalization},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv statistics window tests");
    if (ierr) return (int)ierr;
    ierr = PicurvRunTests("unit-statistics-window", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) { PetscFinalize(); return (int)ierr; }
    ierr = PetscFinalize();
    return (int)ierr;
}
