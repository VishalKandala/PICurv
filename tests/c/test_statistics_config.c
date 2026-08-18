/**
 * @file test_statistics_config.c
 * @brief C unit tests for field-statistics control ingress.
 *
 * These exercise the parse site directly against the PETSc options database, which
 * is the same path a generated `control` file takes. Python rejects malformed
 * configuration earlier and with better messages, so the cases here are the ones a
 * hand-written or hand-edited control file can still reach.
 */

#include "test_support.h"

#include "statistics_config.h"
#include "statistics_window.h"
#include "field_catalog.h"

/** @brief Clears any statistics options left by a previous case. */
static PetscErrorCode ClearStatisticsOptions(void)
{
    PetscFunctionBeginUser;
    PetscCall(PetscOptionsClear(NULL));
    PetscFunctionReturn(0);
}

/** @brief Installs one option as the control file would. */
static PetscErrorCode SetOption(const char *name, const char *value)
{
    PetscFunctionBeginUser;
    PetscCall(PetscOptionsSetValue(NULL, name, value));
    PetscFunctionReturn(0);
}

/** @brief Installs a complete, valid single-window configuration. */
static PetscErrorCode SetValidWindow(void)
{
    PetscFunctionBeginUser;
    PetscCall(ClearStatisticsOptions());
    PetscCall(SetOption("-field_statistics_enabled", "true"));
    PetscCall(SetOption("-field_statistics_window_count", "1"));
    PetscCall(SetOption("-field_statistics_window_0_name", "production"));
    PetscCall(SetOption("-field_statistics_window_0_start_time", "50.0"));
    PetscCall(SetOption("-field_statistics_window_0_end_time", "250.0"));
    PetscCall(SetOption("-field_statistics_window_0_weighting", "physical_time"));
    PetscCall(SetOption("-field_statistics_window_0_step_cadence", "5"));
    PetscCall(SetOption("-field_statistics_window_0_field_count", "2"));
    PetscCall(SetOption("-field_statistics_window_0_field_0_name", "Ucat"));
    PetscCall(SetOption("-field_statistics_window_0_field_0_moments", "first,second"));
    PetscCall(SetOption("-field_statistics_window_0_field_1_name", "P"));
    PetscCall(SetOption("-field_statistics_window_0_field_1_moments", "first"));
    PetscCall(SetOption("-field_statistics_window_0_covariance_count", "1"));
    PetscCall(SetOption("-field_statistics_window_0_covariance_0", "Ucat,P"));
    PetscFunctionReturn(0);
}

/** @brief Runs the parse site against a zeroed context and reports the error code. */
static PetscErrorCode ParseIntoFreshContext(SimCtx *simCtx, PetscErrorCode *result)
{
    PetscFunctionBeginUser;
    PetscCall(PetscMemzero(simCtx, sizeof(*simCtx)));
    PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
    *result = ParseFieldStatisticsConfig(simCtx);
    PetscCall(PetscPopErrorHandler());
    PetscFunctionReturn(0);
}

/** @brief A complete configuration resolves into the definition the pipeline consumes. */
static PetscErrorCode TestResolvesCompleteWindow(void)
{
    SimCtx simCtx;
    PetscBool name_matches = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(SetValidWindow());
    PetscCall(PetscMemzero(&simCtx, sizeof(simCtx)));
    PetscCall(ParseFieldStatisticsConfig(&simCtx));

    PetscCall(PicurvAssertBool(simCtx.fieldStatisticsEnabled, "the enabled flag is read"));
    PetscCall(PicurvAssertIntEqual(1, simCtx.fieldStatisticsWindowCount, "one window is resolved"));
    {
        const PicurvWindowDefinition *d = &simCtx.fieldStatisticsWindows[0].definition;

        PetscCall(PetscStrcmp(d->name, "production", &name_matches));
        PetscCall(PicurvAssertBool(name_matches, "the window name is read"));
        PetscCall(PicurvAssertRealNear(50.0, d->start_time, 1.0e-12, "the start time is read"));
        PetscCall(PicurvAssertBool(d->bounded, "a present end time bounds the window"));
        PetscCall(PicurvAssertRealNear(250.0, d->end_time, 1.0e-12, "the end time is read"));
        PetscCall(PicurvAssertBool((PetscBool)(d->weighting == PICURV_WEIGHTING_PHYSICAL_TIME),
                                   "the weighting is read"));
        PetscCall(PicurvAssertBool((PetscBool)(d->cadence_kind == PICURV_CADENCE_STEP),
                                   "a step cadence selects the step schedule"));
        PetscCall(PicurvAssertIntEqual(5, d->step_cadence, "the step cadence is read"));
        PetscCall(PicurvAssertIntEqual(2, d->field_count, "both fields are resolved"));
        PetscCall(PicurvAssertIntEqual(FIELD_ID_UCAT, d->fields[0].field_id,
                                       "field names resolve through the typed catalog"));
        PetscCall(PicurvAssertBool(d->fields[0].want_second, "a requested second moment is recorded"));
        PetscCall(PicurvAssertIntEqual(FIELD_ID_P, d->fields[1].field_id, "the second field resolves"));
        PetscCall(PicurvAssertBool((PetscBool)!d->fields[1].want_second,
                                   "a field asking only for the first moment keeps no product"));
        PetscCall(PicurvAssertIntEqual(1, d->covariance_count, "the covariance is resolved"));
        PetscCall(PicurvAssertIntEqual(FIELD_ID_UCAT, d->covariances[0].first, "the first member resolves"));
        PetscCall(PicurvAssertIntEqual(FIELD_ID_P, d->covariances[0].second, "the second member resolves"));
    }
    /* The window starts pending: resolution configures, it does not activate. */
    PetscCall(PicurvAssertBool((PetscBool)(simCtx.fieldStatisticsWindows[0].state == PICURV_WINDOW_PENDING),
                               "a resolved window begins pending"));

    PetscCall(DestroyFieldStatisticsConfig(&simCtx));
    PetscCall(PicurvAssertIntEqual(0, simCtx.fieldStatisticsWindowCount, "teardown clears the count"));
    PetscCall(ClearStatisticsOptions());
    PetscFunctionReturn(0);
}

/** @brief An omitted end time leaves the window open rather than defaulting one. */
static PetscErrorCode TestOmittedEndTimeLeavesWindowOpen(void)
{
    SimCtx simCtx;

    PetscFunctionBeginUser;
    PetscCall(SetValidWindow());
    PetscCall(PetscOptionsClearValue(NULL, "-field_statistics_window_0_end_time"));
    PetscCall(PetscMemzero(&simCtx, sizeof(simCtx)));
    PetscCall(ParseFieldStatisticsConfig(&simCtx));

    PetscCall(PicurvAssertBool((PetscBool)!simCtx.fieldStatisticsWindows[0].definition.bounded,
                               "an absent end time leaves the window open ended"));
    PetscCall(DestroyFieldStatisticsConfig(&simCtx));
    PetscCall(ClearStatisticsOptions());
    PetscFunctionReturn(0);
}

/** @brief A disabled subsystem resolves nothing and allocates nothing. */
static PetscErrorCode TestDisabledResolvesNothing(void)
{
    SimCtx simCtx;

    PetscFunctionBeginUser;
    PetscCall(SetValidWindow());
    PetscCall(SetOption("-field_statistics_enabled", "false"));
    PetscCall(PetscMemzero(&simCtx, sizeof(simCtx)));
    PetscCall(ParseFieldStatisticsConfig(&simCtx));

    PetscCall(PicurvAssertIntEqual(0, simCtx.fieldStatisticsWindowCount,
                                   "a disabled subsystem resolves no window"));
    PetscCall(PicurvAssertBool((PetscBool)(simCtx.fieldStatisticsWindows == NULL),
                               "a disabled subsystem allocates nothing"));

    /* Absent entirely is the same as disabled, and must not error. */
    PetscCall(ClearStatisticsOptions());
    PetscCall(PetscMemzero(&simCtx, sizeof(simCtx)));
    PetscCall(ParseFieldStatisticsConfig(&simCtx));
    PetscCall(PicurvAssertIntEqual(0, simCtx.fieldStatisticsWindowCount,
                                   "an absent configuration resolves no window"));
    PetscCall(ClearStatisticsOptions());
    PetscFunctionReturn(0);
}

/** @brief The console cadence rides the same chain and is independent of the windows. */
static PetscErrorCode TestConsoleCadenceIsIndependent(void)
{
    SimCtx simCtx;

    PetscFunctionBeginUser;
    PetscCall(ClearStatisticsOptions());
    PetscCall(SetOption("-statistics_console_output_freq", "25"));
    PetscCall(PetscMemzero(&simCtx, sizeof(simCtx)));
    PetscCall(ParseFieldStatisticsConfig(&simCtx));
    /* Read even when no window exists, because it is reporting configuration rather
     * than science and the banner reports it either way. */
    PetscCall(PicurvAssertIntEqual(25, simCtx.statisticsConsoleOutputFreq,
                                   "the console cadence is read without any window"));
    PetscCall(ClearStatisticsOptions());
    PetscFunctionReturn(0);
}

/** @brief Malformed configuration is rejected rather than silently partially applied. */
static PetscErrorCode TestMalformedConfigurationRejected(void)
{
    SimCtx simCtx;
    PetscErrorCode bad = 0;

    PetscFunctionBeginUser;

    /* Enabled with no window at all. */
    PetscCall(ClearStatisticsOptions());
    PetscCall(SetOption("-field_statistics_enabled", "true"));
    PetscCall(ParseIntoFreshContext(&simCtx, &bad));
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_WRONG, bad, "enabled with no window is rejected"));

    /* Both cadences: the schedule would be ambiguous. */
    PetscCall(SetValidWindow());
    PetscCall(SetOption("-field_statistics_window_0_time_cadence", "0.5"));
    PetscCall(ParseIntoFreshContext(&simCtx, &bad));
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_WRONG, bad, "two cadences are rejected"));

    /* Neither cadence. */
    PetscCall(SetValidWindow());
    PetscCall(PetscOptionsClearValue(NULL, "-field_statistics_window_0_step_cadence"));
    PetscCall(ParseIntoFreshContext(&simCtx, &bad));
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_WRONG, bad, "no cadence is rejected"));

    /* An unknown weighting. */
    PetscCall(SetValidWindow());
    PetscCall(SetOption("-field_statistics_window_0_weighting", "inverse_variance"));
    PetscCall(ParseIntoFreshContext(&simCtx, &bad));
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_WRONG, bad, "an unknown weighting is rejected"));

    /* An unknown moment. */
    PetscCall(SetValidWindow());
    PetscCall(SetOption("-field_statistics_window_0_field_0_moments", "first,third"));
    PetscCall(ParseIntoFreshContext(&simCtx, &bad));
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_WRONG, bad, "an unknown moment is rejected"));

    /* A covariance member missing from the field list has no mean to center against. */
    PetscCall(SetValidWindow());
    PetscCall(SetOption("-field_statistics_window_0_covariance_0", "Ucat,Nvert"));
    PetscCall(ParseIntoFreshContext(&simCtx, &bad));
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_WRONG, bad,
                                   "a covariance member outside the field list is rejected"));

    /* A field paired with itself is the self-product, requested another way. */
    PetscCall(SetValidWindow());
    PetscCall(SetOption("-field_statistics_window_0_covariance_0", "Ucat,Ucat"));
    PetscCall(ParseIntoFreshContext(&simCtx, &bad));
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_WRONG, bad, "a self-pairing covariance is rejected"));

    /* Two windows sharing a name cannot be told apart in a checkpoint. */
    PetscCall(SetValidWindow());
    PetscCall(SetOption("-field_statistics_window_count", "2"));
    PetscCall(SetOption("-field_statistics_window_1_name", "production"));
    PetscCall(SetOption("-field_statistics_window_1_start_time", "0.0"));
    PetscCall(SetOption("-field_statistics_window_1_weighting", "sample"));
    PetscCall(SetOption("-field_statistics_window_1_step_cadence", "1"));
    PetscCall(SetOption("-field_statistics_window_1_field_count", "1"));
    PetscCall(SetOption("-field_statistics_window_1_field_0_name", "P"));
    PetscCall(SetOption("-field_statistics_window_1_field_0_moments", "first"));
    PetscCall(ParseIntoFreshContext(&simCtx, &bad));
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_WRONG, bad, "a duplicate window name is rejected"));

    /* A bounded window that ends before it starts. */
    PetscCall(SetValidWindow());
    PetscCall(SetOption("-field_statistics_window_0_end_time", "10.0"));
    PetscCall(ParseIntoFreshContext(&simCtx, &bad));
    PetscCall(PicurvAssertIntEqual(PETSC_ERR_ARG_OUTOFRANGE, bad,
                                   "an end before the start is rejected by window validation"));

    PetscCall(ClearStatisticsOptions());
    PetscFunctionReturn(0);
}

/**
 * @brief Entry point for the statistics configuration suite.
 */
int main(int argc, char **argv)
{
    PetscErrorCode ierr;
    const PicurvTestCase cases[] = {
        {"resolves-complete-window", TestResolvesCompleteWindow},
        {"omitted-end-time-leaves-window-open", TestOmittedEndTimeLeavesWindowOpen},
        {"disabled-resolves-nothing", TestDisabledResolvesNothing},
        {"console-cadence-is-independent", TestConsoleCadenceIsIndependent},
        {"malformed-configuration-rejected", TestMalformedConfigurationRejected},
    };

    ierr = PetscInitialize(&argc, &argv, NULL, "PICurv statistics configuration tests");
    if (ierr) return (int)ierr;
    ierr = PicurvRunTests("unit-statistics-config", cases, sizeof(cases) / sizeof(cases[0]));
    if (ierr) { PetscFinalize(); return (int)ierr; }
    ierr = PetscFinalize();
    return (int)ierr;
}
