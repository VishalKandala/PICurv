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
 * Every product is one vector carrying all of its components, not one vector per
 * component. A symmetric second-order tensor is a single object: splitting it would
 * cost six memory streams in the per-step accumulation loop instead of one cache
 * line, and six collective gathers per checkpoint instead of one. Component counts
 * that neither `da` nor `fda` provides are carried by a DM mirroring the block
 * decomposition at that degree of freedom.
 */
typedef struct PicurvWindowStorage {
    PetscInt  field_count;       /**< Fields accumulated. */
    PetscInt  covariance_count;  /**< Covariance pairs accumulated. */
    Vec       count;             /**< Per-point accepted sample count. */
    Vec       weight;            /**< Per-point valid weight. */
    Vec       weight_sq;         /**< Per-point squared-weight sum. */
    Vec      *mean;              /**< One per field, matching that field's layout. */
    Vec      *m2;                /**< One per field; NULL when no second moment was requested. */
    Vec      *cm;                /**< One per covariance pair. */
} PicurvWindowStorage;

/** @brief Maximum stored length of a payload name, including the terminator. */
#define PICURV_STATISTICS_PAYLOAD_NAME_LENGTH 96

/**
 * @brief One checkpointable accumulator vector, resolved by enumeration index.
 *
 * The enumeration order is the persistence contract: the manifest inventory, the
 * checkpoint writer, and the restart reader all walk it identically, so a payload
 * lands in the vector it came from without a separate lookup table.
 */
typedef struct {
    char        name[PICURV_STATISTICS_PAYLOAD_NAME_LENGTH]; /**< File basename, no extension. */
    Vec         vec;         /**< Borrowed accumulator vector; never owned by the caller. */
    PetscInt    components;  /**< Degrees of freedom the vector carries. */
    const char *role;        /**< Inventory role: occupancy, mean, second_moment, co_moment. */
    const char *layout;      /**< Catalog layout name for the inventory entry. */
} PicurvStatisticsPayload;

/**
 * @brief Reports how many checkpointable vectors one window's storage holds.
 * @param[in]  storage Storage to measure.
 * @param[out] count   Payload count, including the three occupancy vectors.
 * @return Zero on success, or `PETSC_ERR_ARG_NULL` for a null argument.
 */
PetscErrorCode PicurvWindowStoragePayloadCount(const PicurvWindowStorage *storage, PetscInt *count);

/**
 * @brief Resolves one enumerated payload of a window's storage.
 *
 * Names are derived from catalogued field names and fixed component suffixes, so
 * they are stable across runs, rank counts, and configuration reorderings.
 *
 * @param[in]  user       Block context the storage belongs to.
 * @param[in]  definition Window definition naming the accumulated fields and pairs.
 * @param[in]  storage    Storage to enumerate.
 * @param[in]  index      Payload index in `[0, count)`.
 * @param[out] payload    Resolved payload; the vector is borrowed, not duplicated.
 * @return Zero on success, or `PETSC_ERR_ARG_OUTOFRANGE` for an index outside the range.
 */
PetscErrorCode PicurvWindowStoragePayload(UserCtx *user, const PicurvWindowDefinition *definition,
                                          const PicurvWindowStorage *storage, PetscInt index,
                                          PicurvStatisticsPayload *payload);

/**
 * @brief Resolves the DM carrying a given number of accumulator components.
 *
 * Component counts of one and three reuse the DMs the block already owns; six is
 * carried by the symmetric-tensor DM created alongside them. Every one of these
 * mirrors the block decomposition exactly, so a pointwise loop can read a source
 * field and write an accumulator at the same index.
 *
 * @param[in]  user       Block context owning the DMs.
 * @param[in]  components Component count to place.
 * @param[out] dm         Resolved DM; borrowed, never destroyed by the caller.
 * @return Zero on success, or `PETSC_ERR_ARG_OUTOFRANGE` for an unsupported count.
 */
PetscErrorCode PicurvStatisticsComponentDM(UserCtx *user, PetscInt components, DM *dm);

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

/** @brief Tolerance within which a negative variance is treated as floating-point noise. */
#define PICURV_STATISTICS_VARIANCE_FLOOR 1.0e-12

/** @brief One derived output field, resolved by enumeration index. */
typedef struct {
    char     name[PICURV_STATISTICS_PAYLOAD_NAME_LENGTH]; /**< Output field name, window qualified. */
    PetscInt components;                                   /**< One or three. */
} PicurvDerivedField;

/**
 * @brief Reports how many derived fields a requested output set produces.
 * @param[in]  definition Window definition naming the accumulated fields and pairs.
 * @param[in]  storage    Accumulator state the outputs are derived from.
 * @param[in]  outputs    Comma-separated output kinds: mean, reynolds_stress, rms, tke, flux.
 * @param[out] count      Number of derived fields.
 * @return Zero on success, or `PETSC_ERR_ARG_WRONG` for an unknown output kind.
 */
PetscErrorCode PicurvWindowDerivedCount(const PicurvWindowDefinition *definition,
                                        const PicurvWindowStorage *storage,
                                        const char *outputs, PetscInt *count);

/**
 * @brief Derives one output field from centered accumulator state.
 *
 * Normalizes in exactly one place: `R_ij = C_ij / W`, `RMS_i = sqrt(R_ii)`,
 * `k = (R_xx + R_yy + R_zz) / 2`, and a flux is the co-moment over the same weight.
 * Each uses the moment kernels rather than repeating the division, so the online and
 * offline halves of the pipeline cannot disagree about what centered state means.
 *
 * A variance that comes out slightly negative through floating-point cancellation is
 * clamped only where a square root would otherwise fail, and only within
 * `PICURV_STATISTICS_VARIANCE_FLOOR`. Stored state is never modified.
 *
 * Points the window never sampled are left at zero rather than divided by a zero
 * weight; the valid-fraction range reports how much of the domain that covers.
 *
 * @param[in]  user           Block context supplying the target domain.
 * @param[in]  definition     Window definition.
 * @param[in]  storage        Accumulator state to read.
 * @param[in]  outputs        Comma-separated output kinds.
 * @param[in]  index          Derived field index in `[0, count)`.
 * @param[out] scalar_target  Scalar destination, used when the field has one component.
 * @param[out] vector_target  Vector destination, used when the field has three.
 * @param[out] field          Resolved name and component count of the derived field.
 * @return Zero on success, or a PETSc error.
 */
PetscErrorCode PicurvWindowDerive(UserCtx *user, const PicurvWindowDefinition *definition,
                                  const PicurvWindowStorage *storage, const char *outputs,
                                  PetscInt index, Vec scalar_target, Vec vector_target,
                                  PicurvDerivedField *field);

/**
 * @brief Reports the spatial mean of a derived field over the points a window sampled.
 *
 * The average is taken over targeted points that actually accumulated weight, not
 * over the whole vector. A derived field is zero everywhere outside the target
 * domain and at any point the mask never admitted, and those zeros are absences
 * rather than measurements: including them would scale the answer down by the
 * fraction of the vector the window never covered.
 *
 * Performs a collective reduction, so callers use it for reporting rather than per
 * point.
 *
 * @param[in]  user       Block context supplying the target domain.
 * @param[in]  definition Window definition naming the accumulated fields.
 * @param[in]  storage    Accumulator state supplying per-point occupancy.
 * @param[in]  field      Derived field to average, on the cell-centred scalar DM.
 * @param[out] mean       Spatial mean; zero when the window sampled no point.
 * @return Zero on success, or a PETSc error.
 */
PetscErrorCode PicurvWindowSpatialMean(UserCtx *user, const PicurvWindowDefinition *definition,
                                       const PicurvWindowStorage *storage, Vec field,
                                       PetscReal *mean);

/**
 * @brief Reports the range of per-point valid fraction across a window's domain.
 *
 * A point contributes only to the states in which the mask accepted it, so with a
 * moving immersed body different points carry different sample counts. The ratio of
 * a point's own count to the window's accepted-sample count is its valid fraction,
 * and the range of that ratio is the window's mask-health indicator: a minimum of
 * one means every point saw every state, and a minimum of zero means some point
 * contributed nothing at all.
 *
 * Performs a collective reduction, so callers gate it on an already-active reporting
 * path rather than calling it every step.
 *
 * @param[in]  user         Block context supplying the target domain and mask.
 * @param[in]  definition   Window definition naming the accumulated fields.
 * @param[in]  storage      Accumulator state to inspect.
 * @param[in]  sample_count Accepted states the window has recorded.
 * @param[out] minimum      Smallest valid fraction; one when no point is targeted.
 * @param[out] maximum      Largest valid fraction; zero when no point is targeted.
 * @return Zero on success, or a PETSc error.
 */
PetscErrorCode PicurvWindowValidFractionRange(UserCtx *user, const PicurvWindowDefinition *definition,
                                              const PicurvWindowStorage *storage, PetscInt sample_count,
                                              PetscReal *minimum, PetscReal *maximum);

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
