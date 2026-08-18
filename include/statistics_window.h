/**
 * @file statistics_window.h
 * @brief Window lifecycle, scheduling, and weighting for the field-statistics pipeline.
 *
 * Implements the window semantics fixed in
 * @ref 60_Field_Statistics_Phase2_Implementation_Specification sections 2 and 4:
 * right-rectangle weighting, final-interval clipping, and the rule that a state
 * representing a zero-length interval is not a sample.
 *
 * This module decides **whether** a completed state is accepted and **what weight**
 * it carries. It holds no PETSc objects and performs no field accumulation; the
 * caller applies the returned weight through the moment kernels.
 */

#ifndef PICURV_STATISTICS_WINDOW_H
#define PICURV_STATISTICS_WINDOW_H

#include <petscsys.h>

struct SimCtx;

/** @brief Maximum stored length of a window name, including the terminator. */
#define PICURV_WINDOW_NAME_LENGTH 64

/** @brief Lifecycle state of one window. */
typedef enum {
    PICURV_WINDOW_PENDING = 0,  /**< Requested start not yet reached. */
    PICURV_WINDOW_ACTIVE,       /**< Accepting due states. */
    PICURV_WINDOW_COMPLETE      /**< Bounded end reached; accepts nothing further. */
} PicurvWindowState;

/** @brief How an accepted state's weight is determined. */
typedef enum {
    PICURV_WEIGHTING_SAMPLE = 0,    /**< Equal weight per accepted state. */
    PICURV_WEIGHTING_PHYSICAL_TIME  /**< Weight is the represented interval. */
} PicurvWeighting;

/** @brief Which schedule selects due states. Exactly one is used. */
typedef enum {
    PICURV_CADENCE_STEP = 0,  /**< Every n completed steps from activation. */
    PICURV_CADENCE_TIME       /**< First state at or past each nominal time target. */
} PicurvCadenceKind;

/** @brief Maximum fields or covariance pairs one window may request. */
#define PICURV_WINDOW_MAX_REQUESTS 16

/** @brief One field a window accumulates. The first moment is always kept. */
typedef struct {
    PetscInt  field_id;     /**< Catalogued Eulerian field identity. */
    PetscBool want_second;  /**< Also keep the centered second moment. */
} PicurvWindowFieldRequest;

/** @brief One cross-field covariance a window accumulates. */
typedef struct {
    PetscInt first;   /**< First member; must also appear in the field list. */
    PetscInt second;  /**< Second member; must also appear in the field list. */
} PicurvWindowCovarianceRequest;

/** @brief The scientifically immutable definition of one window. */
typedef struct {
    char              name[PICURV_WINDOW_NAME_LENGTH];
    PetscReal         start_time;    /**< Requested start. */
    PetscReal         end_time;      /**< Requested end; ignored when @c bounded is false. */
    PetscBool         bounded;       /**< False for an open-ended window. */
    PicurvWeighting   weighting;
    PicurvCadenceKind cadence_kind;
    PetscInt          step_cadence;  /**< Used when cadence_kind is step; must be positive. */
    PetscReal         time_cadence;  /**< Used when cadence_kind is time; must be positive. */
    PetscInt                      field_count;
    PicurvWindowFieldRequest      fields[PICURV_WINDOW_MAX_REQUESTS];
    PetscInt                      covariance_count;
    PicurvWindowCovarianceRequest covariances[PICURV_WINDOW_MAX_REQUESTS];
} PicurvWindowDefinition;

/** @brief Runtime state of one window. */
typedef struct PicurvWindow {
    PicurvWindowDefinition definition;
    PicurvWindowState      state;
    PetscReal              effective_start;      /**< Origin of the first represented interval. */
    PetscReal              effective_end;        /**< End of the last represented interval. */
    PetscReal              last_accepted_time;   /**< Right edge of the last represented interval. */
    PetscInt               sample_count;
    PetscReal              total_weight;
    PetscReal              represented_time;     /**< Physical time the window covers. */
    PetscInt               activation_step;      /**< Step at which the window became active. */
    PetscInt               last_event_step;      /**< Guards against a step being offered twice. */
    PetscInt               next_time_target;     /**< k in effective_start + k*time_cadence. */
    PetscInt               restart_count;        /**< Restart segments this state descends from. */
} PicurvWindow;

/**
 * @brief Validates a definition and initializes a window to the pending state.
 * @param[out] window     Window to initialize.
 * @param[in]  definition Requested definition; copied into the window.
 * @return Zero on success, or `PETSC_ERR_ARG_OUTOFRANGE` for a non-positive cadence,
 *         an empty name, or a bounded window whose end does not exceed its start.
 */
PetscErrorCode PicurvWindowInit(PicurvWindow *window, const PicurvWindowDefinition *definition);

/**
 * @brief Offers one completed state to a window and reports the decision.
 *
 * Applies the interval convention in full: the state carries the interval ending
 * at it, measured from the previous accepted state or from the effective start;
 * a zero-length interval is not a sample; and a bounded window clips its final
 * interval to the requested end and then completes.
 *
 * When @p accepted is returned true the window's bookkeeping has already been
 * advanced, and the caller applies @p weight through the moment kernels. When it
 * is false the window is scientifically unchanged.
 *
 * Offering the same step twice is rejected, so a completed state cannot be
 * counted more than once.
 *
 * @param[in,out] window   Window to offer the state to.
 * @param[in]     step     Completed step number.
 * @param[in]     time     Physical time of the completed state.
 * @param[out]    accepted Whether the state became a sample.
 * @param[out]    weight   Weight to apply; zero when not accepted.
 * @return Zero on success, or a PETSc error for a null argument.
 */
PetscErrorCode PicurvWindowOfferState(PicurvWindow *window, PetscInt step, PetscReal time,
                                      PetscBool *accepted, PetscReal *weight);

/** @brief Number of independently hashed property groups in a window definition. */
#define PICURV_WINDOW_HASH_GROUP_COUNT 8

/** @brief Stored length of one truncated group digest, including the terminator. */
#define PICURV_WINDOW_HASH_GROUP_LENGTH 17

/**
 * @brief Computes the resolved identity hash of one window definition.
 *
 * Hashes the canonical serialization defined in
 * @ref 60_Field_Statistics_Phase2_Implementation_Specification section 7, in that
 * fixed order, so a saved window can be matched against a resolved one without
 * storing the definition itself.
 *
 * `end_time` and the enabled flag are deliberately excluded, which is what lets a
 * bounded window be extended forward and lets statistics be switched off and on
 * without invalidating saved state.
 *
 * Field and covariance entries are serialized in catalog order rather than the
 * order the user listed them, so a reordered but otherwise identical configuration
 * continues rather than being rejected.
 *
 * Each property group is additionally hashed on its own. A restart that finds a
 * mismatched full digest compares the group digests to name the first differing
 * property, which a single digest could not do.
 *
 * @param[in]  definition       Window definition to hash.
 * @param[out] digest_hex       Full 64-character digest plus terminator.
 * @param[out] group_digest_hex Optional per-group truncated digests; pass NULL to skip.
 * @return Zero on success, or a PETSc error for a null argument or unknown field.
 */
PetscErrorCode PicurvWindowComputeHash(const PicurvWindowDefinition *definition,
                                       char digest_hex[65],
                                       char group_digest_hex[][PICURV_WINDOW_HASH_GROUP_LENGTH]);

/**
 * @brief Reports which hashed property group first differs from saved group digests.
 *
 * A checkpoint stores the group digests but never the definition itself, so this is
 * what turns "two hashes differ" into a message naming the property that changed.
 *
 * @param[in]  definition          Resolved definition to compare against.
 * @param[in]  saved_group_digests Comma-separated group digests from a checkpoint.
 * @param[out] group               First differing group index, or -1 when the saved
 *                                 digests match or are too malformed to compare.
 * @return Zero on success, or a PETSc error for a null argument or unknown field.
 */
PetscErrorCode PicurvWindowFirstHashDifference(const PicurvWindowDefinition *definition,
                                               const char *saved_group_digests,
                                               PetscInt *group);

/**
 * @brief Returns the stable name of one hashed property group.
 * @param[in] group Group index in `[0, PICURV_WINDOW_HASH_GROUP_COUNT)`.
 * @return Static string; `"unknown"` for an out-of-range index, never NULL.
 */
const char *PicurvWindowHashGroupName(PetscInt group);

/**
 * @brief Reports the fraction of a bounded window's span that has been represented.
 * @param[in] window Window to query.
 * @return Value in [0,1] for a bounded window, or zero for an open one.
 */
PetscReal PicurvWindowProgress(const PicurvWindow *window);

/**
 * @brief Returns a stable human-readable name for a window state.
 * @param[in] state Window lifecycle state.
 * @return Static string; never NULL.
 */
const char *PicurvWindowStateName(PicurvWindowState state);

/**
 * @brief Reports whether this run has live field-statistics state.
 *
 * The subsystem is active only when it is enabled, at least one window is
 * configured, and the window array exists. Every caller that touches window or
 * accumulator state asks this rather than restating the condition, so the three
 * parts cannot drift apart between the runloop, the checkpoint writer, and the
 * console monitor.
 *
 * @param[in] simCtx Simulation context; may be NULL.
 * @return `PETSC_TRUE` when window state exists and may be touched.
 */
PetscBool FieldStatisticsIsActive(const struct SimCtx *simCtx);

/**
 * @brief Offers one completed state to every configured window.
 *
 * Called once per completed step from the runloop. Each due window advances its
 * own bookkeeping independently; windows share the source state but never share
 * accumulator state. Does nothing when field statistics are disabled or no
 * window is configured, which is the case until configuration ingress exists.
 *
 * @param[in,out] simCtx Simulation context carrying the window array.
 * @param[in]     step   Completed step number.
 * @param[in]     time   Physical time of the completed state.
 * @return Zero on success, or a PETSc error propagated from a window update.
 */
PetscErrorCode FieldStatisticsUpdateWindows(struct SimCtx *simCtx, PetscInt step, PetscReal time);

#endif /* PICURV_STATISTICS_WINDOW_H */
