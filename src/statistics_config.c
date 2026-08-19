/**
 * @file statistics_config.c
 * @brief Control ingress for the field-statistics pipeline.
 *
 * Full API contract is documented with the declarations in
 * `include/statistics_config.h`.
 */

#include "statistics_config.h"
#include "statistics_window.h"
#include "field_catalog.h"
#include "logging.h"
#include "io.h"

/** @brief Longest constructed option name this module builds. */
#define STATISTICS_OPTION_NAME_LENGTH 128

/** @brief Longest option value this module reads, sized for a moment or pair list. */
#define STATISTICS_OPTION_VALUE_LENGTH 256

/**
 * @brief Internal helper: resolves one field name against the typed catalog.
 * @details Local to this translation unit. The error names the window so a
 *          hand-written control file reports the same context Python would.
 */
static PetscErrorCode ResolveStatisticsField(const char *window_name, const char *field_name,
                                             PetscInt *field_id)
{
    FieldId resolved = FIELD_ID_UCAT;

    PetscFunctionBeginUser;
    PetscCheck(field_name != NULL && field_name[0] != '\0', PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
               "Statistics window '%s' names an empty field.", window_name);
    PetscCall(FieldIdFromName(field_name, &resolved));
    *field_id = (PetscInt)resolved;
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper: splits a comma-separated value into trimmed tokens.
 * @details Local to this translation unit. Rewrites @p value in place, so callers
 *          pass a mutable buffer they own.
 */
static PetscErrorCode SplitStatisticsList(char *value, char **tokens, PetscInt capacity,
                                          PetscInt *count)
{
    char *cursor = value;

    PetscFunctionBeginUser;
    *count = 0;
    while (cursor && *cursor) {
        char *comma = strchr(cursor, ',');

        if (comma) *comma = '\0';
        TrimWhitespace(cursor);
        if (cursor[0] != '\0') {
            PetscCheck(*count < capacity, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                       "A statistics list holds at most %" PetscInt_FMT " entries.", capacity);
            tokens[(*count)++] = cursor;
        }
        cursor = comma ? comma + 1 : NULL;
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper: reads one window's field list into its definition.
 * @details Local to this translation unit.
 */
static PetscErrorCode ParseStatisticsFields(PetscInt window, PicurvWindowDefinition *definition)
{
    char name[STATISTICS_OPTION_NAME_LENGTH];
    PetscBool found = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(PetscSNPrintf(name, sizeof(name),
                            "-field_statistics_window_%" PetscInt_FMT "_field_count", window));
    PetscCall(PetscOptionsGetInt(NULL, NULL, name, &definition->field_count, &found));
    PetscCheck(found && definition->field_count > 0, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
               "Statistics window '%s' requests no field.", definition->name);
    PetscCheck(definition->field_count <= PICURV_WINDOW_MAX_REQUESTS, PETSC_COMM_WORLD,
               PETSC_ERR_ARG_OUTOFRANGE,
               "Statistics window '%s' requests %" PetscInt_FMT " fields; at most %d are supported.",
               definition->name, definition->field_count, PICURV_WINDOW_MAX_REQUESTS);

    for (PetscInt field_index = 0; field_index < definition->field_count; ++field_index) {
        char field_name[STATISTICS_OPTION_VALUE_LENGTH];
        char moments[STATISTICS_OPTION_VALUE_LENGTH];
        char *tokens[4];
        PetscInt token_count = 0;

        PetscCall(PetscSNPrintf(name, sizeof(name),
                                "-field_statistics_window_%" PetscInt_FMT "_field_%" PetscInt_FMT "_name",
                                window, field_index));
        PetscCall(PetscOptionsGetString(NULL, NULL, name, field_name, sizeof(field_name), &found));
        PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                   "Statistics window '%s' field %" PetscInt_FMT " has no name.",
                   definition->name, field_index);
        PetscCall(ResolveStatisticsField(definition->name, field_name,
                                         &definition->fields[field_index].field_id));

        PetscCall(PetscSNPrintf(name, sizeof(name),
                                "-field_statistics_window_%" PetscInt_FMT "_field_%" PetscInt_FMT "_moments",
                                window, field_index));
        PetscCall(PetscOptionsGetString(NULL, NULL, name, moments, sizeof(moments), &found));
        PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                   "Statistics window '%s' field '%s' requests no moment.", definition->name, field_name);
        PetscCall(SplitStatisticsList(moments, tokens, 4, &token_count));
        PetscCheck(token_count > 0, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                   "Statistics window '%s' field '%s' requests no moment.", definition->name, field_name);

        definition->fields[field_index].want_second = PETSC_FALSE;
        for (PetscInt t = 0; t < token_count; ++t) {
            if (!strcmp(tokens[t], "first")) continue;
            if (!strcmp(tokens[t], "second")) {
                definition->fields[field_index].want_second = PETSC_TRUE;
                continue;
            }
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                    "Statistics window '%s' field '%s' requests unknown moment '%s'; "
                    "the accumulable moments are 'first' and 'second'.",
                    definition->name, field_name, tokens[t]);
        }
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper: reads one window's covariance list into its definition.
 * @details Local to this translation unit.
 */
static PetscErrorCode ParseStatisticsCovariances(PetscInt window, PicurvWindowDefinition *definition)
{
    char name[STATISTICS_OPTION_NAME_LENGTH];
    PetscBool found = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(PetscSNPrintf(name, sizeof(name),
                            "-field_statistics_window_%" PetscInt_FMT "_covariance_count", window));
    PetscCall(PetscOptionsGetInt(NULL, NULL, name, &definition->covariance_count, &found));
    if (!found) definition->covariance_count = 0;
    PetscCheck(definition->covariance_count >= 0 &&
               definition->covariance_count <= PICURV_WINDOW_MAX_REQUESTS,
               PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
               "Statistics window '%s' requests %" PetscInt_FMT " covariances; at most %d are supported.",
               definition->name, definition->covariance_count, PICURV_WINDOW_MAX_REQUESTS);

    for (PetscInt pair_index = 0; pair_index < definition->covariance_count; ++pair_index) {
        char pair[STATISTICS_OPTION_VALUE_LENGTH];
        char *tokens[4];
        PetscInt token_count = 0;
        PetscBool first_present = PETSC_FALSE;
        PetscBool second_present = PETSC_FALSE;

        PetscCall(PetscSNPrintf(name, sizeof(name),
                                "-field_statistics_window_%" PetscInt_FMT "_covariance_%" PetscInt_FMT,
                                window, pair_index));
        PetscCall(PetscOptionsGetString(NULL, NULL, name, pair, sizeof(pair), &found));
        PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                   "Statistics window '%s' covariance %" PetscInt_FMT " has no members.",
                   definition->name, pair_index);
        PetscCall(SplitStatisticsList(pair, tokens, 4, &token_count));
        PetscCheck(token_count == 2, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                   "Statistics window '%s' covariance %" PetscInt_FMT " needs exactly two members.",
                   definition->name, pair_index);
        PetscCall(ResolveStatisticsField(definition->name, tokens[0], &definition->covariances[pair_index].first));
        PetscCall(ResolveStatisticsField(definition->name, tokens[1], &definition->covariances[pair_index].second));
        PetscCheck(definition->covariances[pair_index].first != definition->covariances[pair_index].second,
                   PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                   "Statistics window '%s' pairs '%s' with itself; that is the self-product already "
                   "requested through moments: [second].", definition->name, tokens[0]);

        /* Both members must carry a running mean, because a co-moment centers
         * against them. This is the C-side guard behind the same Python rule. */
        for (PetscInt field_index = 0; field_index < definition->field_count; ++field_index) {
            if (definition->fields[field_index].field_id == definition->covariances[pair_index].first) first_present = PETSC_TRUE;
            if (definition->fields[field_index].field_id == definition->covariances[pair_index].second) second_present = PETSC_TRUE;
        }
        PetscCheck(first_present && second_present, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                   "Statistics window '%s' pairs '%s' with '%s', but both must also appear in its "
                   "field list so their means exist to center against.",
                   definition->name, tokens[0], tokens[1]);
    }
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper: reads one window's scheduling and weighting properties.
 * @details Local to this translation unit.
 */
static PetscErrorCode ParseStatisticsWindowSchedule(PetscInt window, PicurvWindowDefinition *definition)
{
    char name[STATISTICS_OPTION_NAME_LENGTH];
    char weighting[STATISTICS_OPTION_VALUE_LENGTH];
    PetscBool found = PETSC_FALSE;
    PetscBool step_found = PETSC_FALSE;
    PetscBool time_found = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(PetscSNPrintf(name, sizeof(name),
                            "-field_statistics_window_%" PetscInt_FMT "_start_time", window));
    PetscCall(PetscOptionsGetReal(NULL, NULL, name, &definition->start_time, &found));
    PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
               "Statistics window '%s' has no start time.", definition->name);

    /* An absent end time is what makes a window open ended, so its absence is
     * meaningful rather than a defaulted value. */
    PetscCall(PetscSNPrintf(name, sizeof(name),
                            "-field_statistics_window_%" PetscInt_FMT "_end_time", window));
    PetscCall(PetscOptionsGetReal(NULL, NULL, name, &definition->end_time, &definition->bounded));

    PetscCall(PetscSNPrintf(name, sizeof(name),
                            "-field_statistics_window_%" PetscInt_FMT "_weighting", window));
    PetscCall(PetscOptionsGetString(NULL, NULL, name, weighting, sizeof(weighting), &found));
    PetscCheck(found, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
               "Statistics window '%s' has no weighting.", definition->name);
    if (!strcmp(weighting, "sample")) {
        definition->weighting = PICURV_WEIGHTING_SAMPLE;
    } else if (!strcmp(weighting, "physical_time")) {
        definition->weighting = PICURV_WEIGHTING_PHYSICAL_TIME;
    } else {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                "Statistics window '%s' has weighting '%s'; only 'sample' and 'physical_time' exist.",
                definition->name, weighting);
    }

    PetscCall(PetscSNPrintf(name, sizeof(name),
                            "-field_statistics_window_%" PetscInt_FMT "_step_cadence", window));
    PetscCall(PetscOptionsGetInt(NULL, NULL, name, &definition->step_cadence, &step_found));
    PetscCall(PetscSNPrintf(name, sizeof(name),
                            "-field_statistics_window_%" PetscInt_FMT "_time_cadence", window));
    PetscCall(PetscOptionsGetReal(NULL, NULL, name, &definition->time_cadence, &time_found));
    PetscCheck(step_found != time_found, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
               "Statistics window '%s' must select exactly one of step_cadence and time_cadence.",
               definition->name);
    definition->cadence_kind = step_found ? PICURV_CADENCE_STEP : PICURV_CADENCE_TIME;
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper: reports one resolved window definition to the log.
 * @details Local to this translation unit. Runs once per window at startup, so a
 *          run's log records exactly what it accumulated without having to be
 *          cross-referenced against the configuration that produced it.
 */
static PetscErrorCode LogResolvedWindow(const PicurvWindowDefinition *definition)
{
    char fields[STATISTICS_OPTION_VALUE_LENGTH] = "";
    size_t used = 0;

    PetscFunctionBeginUser;
    for (PetscInt field_index = 0; field_index < definition->field_count; ++field_index) {
        PetscCall(PetscStrlen(fields, &used));
        PetscCall(PetscSNPrintf(fields + used, sizeof(fields) - used, "%s%s%s",
                                field_index ? " " : "",
                                FieldCanonicalName((FieldId)definition->fields[field_index].field_id),
                                definition->fields[field_index].want_second ? "(1,2)" : "(1)"));
    }
    if (definition->cadence_kind == PICURV_CADENCE_STEP) {
        LOG_ALLOW(GLOBAL, LOG_INFO,
                  "Statistics window '%s': t=[%.6g,%s] weighting=%s every %d step(s), fields: %s, "
                  "%d covariance(s).\n",
                  definition->name, (double)definition->start_time,
                  definition->bounded ? "bounded" : "open",
                  definition->weighting == PICURV_WEIGHTING_SAMPLE ? "sample" : "physical_time",
                  (int)definition->step_cadence, fields, (int)definition->covariance_count);
    } else {
        LOG_ALLOW(GLOBAL, LOG_INFO,
                  "Statistics window '%s': t=[%.6g,%s] weighting=%s every %.6g of time, fields: %s, "
                  "%d covariance(s).\n",
                  definition->name, (double)definition->start_time,
                  definition->bounded ? "bounded" : "open",
                  definition->weighting == PICURV_WEIGHTING_SAMPLE ? "sample" : "physical_time",
                  (double)definition->time_cadence, fields, (int)definition->covariance_count);
    }
    if (definition->bounded) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Statistics window '%s' ends at t=%.6g.\n",
                  definition->name, (double)definition->end_time);
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ParseFieldStatisticsConfig"
/**
 * @brief Implementation of \ref ParseFieldStatisticsConfig().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/statistics_config.h`.
 * @see ParseFieldStatisticsConfig()
 */
PetscErrorCode ParseFieldStatisticsConfig(SimCtx *simCtx)
{
    PetscInt window_count = 0;
    PetscBool found = PETSC_FALSE;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCheck(simCtx != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "SimCtx cannot be NULL.");

    PetscCall(PetscOptionsGetInt(NULL, NULL, "-statistics_console_output_freq",
                                 &simCtx->statisticsConsoleOutputFreq, NULL));
    PetscCall(PetscOptionsGetBool(NULL, NULL, "-field_statistics_enabled",
                                  &simCtx->fieldStatisticsEnabled, NULL));
    if (!simCtx->fieldStatisticsEnabled) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Field statistics are disabled.\n");
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    PetscCall(PetscOptionsGetInt(NULL, NULL, "-field_statistics_window_count", &window_count, &found));
    PetscCheck(found && window_count > 0, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
               "Field statistics are enabled but no window is configured.");
    PetscCall(PetscCalloc1((size_t)window_count, &simCtx->fieldStatisticsWindows));
    simCtx->fieldStatisticsWindowCount = window_count;

    for (PetscInt window_index = 0; window_index < window_count; ++window_index) {
        PicurvWindowDefinition definition;
        char option_name[STATISTICS_OPTION_NAME_LENGTH];

        PetscCall(PetscMemzero(&definition, sizeof(definition)));
        PetscCall(PetscSNPrintf(option_name, sizeof(option_name),
                                "-field_statistics_window_%" PetscInt_FMT "_name", window_index));
        PetscCall(PetscOptionsGetString(NULL, NULL, option_name, definition.name,
                                        sizeof(definition.name), &found));
        PetscCheck(found && definition.name[0] != '\0', PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                   "Statistics window %" PetscInt_FMT " has no name.", window_index);

        /* A duplicate name would make two windows indistinguishable in the
         * checkpoint, where window state is matched by name on restart. */
        for (PetscInt prior = 0; prior < window_index; ++prior) {
            PetscBool clash = PETSC_FALSE;

            PetscCall(PetscStrcmp(simCtx->fieldStatisticsWindows[prior].definition.name,
                                  definition.name, &clash));
            PetscCheck(!clash, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,
                       "Statistics window name '%s' is used more than once; names identify saved "
                       "state across a restart and must be unique.", definition.name);
        }

        PetscCall(ParseStatisticsWindowSchedule(window_index, &definition));
        PetscCall(ParseStatisticsFields(window_index, &definition));
        PetscCall(ParseStatisticsCovariances(window_index, &definition));
        PetscCall(PicurvWindowInit(&simCtx->fieldStatisticsWindows[window_index], &definition));
        PetscCall(LogResolvedWindow(&simCtx->fieldStatisticsWindows[window_index].definition));
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Field statistics resolved %d window(s).\n", (int)window_count);
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref DestroyFieldStatisticsConfig().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/statistics_config.h`.
 * @see DestroyFieldStatisticsConfig()
 */
PetscErrorCode DestroyFieldStatisticsConfig(SimCtx *simCtx)
{
    PetscFunctionBeginUser;
    if (!simCtx || !simCtx->fieldStatisticsWindows) PetscFunctionReturn(0);
    PetscCall(PetscFree(simCtx->fieldStatisticsWindows));
    simCtx->fieldStatisticsWindows = NULL;
    simCtx->fieldStatisticsWindowCount = 0;
    PetscFunctionReturn(0);
}
