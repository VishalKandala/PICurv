/**
 * @file statistics_config.h
 * @brief Control ingress for the field-statistics pipeline.
 *
 * Resolves the generated master `control` into the window definitions the rest of
 * the pipeline consumes. The user-facing keys these options come from are
 * documented in @ref p09_field_statistics_sec.
 *
 * A window list is variable arity, so its option names are constructed rather than
 * literal. Every constructed name belongs to a family declared in the ingress audit
 * manifest, which is what keeps the audit's guarantee intact rather than opening a
 * hole in it.
 *
 * This module owns only configuration. It creates no PETSc vectors: the accumulator
 * factory runs later and sizes itself from the definitions resolved here.
 */

#ifndef PICURV_STATISTICS_CONFIG_H
#define PICURV_STATISTICS_CONFIG_H

#include "variables.h"

/**
 * @brief Resolves field-statistics configuration from the control file.
 *
 * Reads the option families, resolves every field name through the typed catalog,
 * validates each definition through `PicurvWindowInit`, and stores the resulting
 * windows on the context. Each resolved definition is reported at `LOG_INFO` through
 * `LOG_ALLOW`, so it carries the usual dual gate and is surfaced by allow-listing the
 * emitting function. Whether statistics were active at all is recorded
 * unconditionally by the startup banner.
 *
 * Leaves the context untouched when statistics are disabled or no window is
 * configured, so a run that does not ask for statistics allocates nothing.
 *
 * Must run before `CreateAndInitializeAllVectors`, which sizes the per-window
 * accumulators from the window count resolved here.
 *
 * @param[in,out] simCtx Simulation context to populate.
 * @return Zero on success, or a PETSc error for a malformed or inconsistent
 *         configuration. Python validation rejects these earlier with better
 *         messages; these checks guard a hand-written control file.
 */
PetscErrorCode ParseFieldStatisticsConfig(SimCtx *simCtx);

/**
 * @brief Releases the window definitions resolved by ParseFieldStatisticsConfig().
 * @param[in,out] simCtx Simulation context to clear; safe when nothing was resolved.
 * @return Zero on success.
 */
PetscErrorCode DestroyFieldStatisticsConfig(SimCtx *simCtx);

#endif /* PICURV_STATISTICS_CONFIG_H */
