#ifndef PARTICLE_STATISTICS_H
#define PARTICLE_STATISTICS_H

#include "variables.h"
#include "logging.h"

/**
 * @file particle_statistics.h
 * @brief Global statistics kernels for the Statistics Pipeline.
 *
 * These kernels compute global aggregate quantities via MPI reduction — they do NOT
 * produce per-particle VTK data.  Each kernel appends one CSV row per call and logs
 * a one-line summary via LOG_INFO.
 *
 * Output file convention:  {stats_prefix}_{kernel_name}.csv
 * All MPI_Allreduce operations are internal; file I/O is rank-0-only.
 *
 * To add a new statistic:
 *   1. Implement PetscErrorCode ComputeXxx(UserCtx*, const char*, PetscInt) in
 *      src/particle_statistics.c
 *   2. Declare it below.
 *   3. Add one else-if in GlobalStatisticsPipeline() in src/postprocessor.c
 *   No other files need to change.
 */

/**
 * @brief Computes the mean-squared displacement (MSD) of a particle cloud.
 *
 * Reference point r0 = (simCtx->psrc_x, psrc_y, psrc_z).
 * Diffusivity D = 1 / (Re * Sc), time t = ti * dt.
 * Computes MSD_x, MSD_y, MSD_z (isotropy check), MSD_total, r_rms,
 * centre-of-mass drift, and fractions inside 1σ/2σ/3σ theoretical shells.
 * Appends one row to {stats_prefix}_msd.csv; writes header on first call.
 *
 * @param user         The UserCtx containing the DMSwarm (user->swarm).
 * @param stats_prefix Base filename prefix (e.g. "brownian_stats").
 * @param ti           Current time-step index.
 * @return PetscErrorCode
 */
PetscErrorCode ComputeParticleMSD(UserCtx *user, const char *stats_prefix, PetscInt ti);

/* -----------------------------------------------------------------------
 * Future kernels — add declarations here:
 *
 * PetscErrorCode ComputeVelocityPDF(UserCtx *user, const char *stats_prefix, PetscInt ti);
 * PetscErrorCode ComputeConcentrationField(UserCtx *user, const char *stats_prefix, PetscInt ti);
 * PetscErrorCode ComputeResidenceTimeDistribution(UserCtx *user, const char *stats_prefix, PetscInt ti);
 * ----------------------------------------------------------------------- */

#endif /* PARTICLE_STATISTICS_H */
