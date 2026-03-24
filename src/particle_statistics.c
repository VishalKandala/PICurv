/**
 * @file particle_statistics.c
 * @brief Global statistical reduction kernels for the Statistics Pipeline.
 *
 * Each kernel in this file:
 *  - Performs an MPI_Allreduce over all particles across all ranks.
 *  - Maintains one logical CSV row per processed step in its output file.
 *  - Logs a one-line summary via LOG_INFO on rank 0.
 *  - Produces O(1) output regardless of particle count, making these kernels
 *    the primary quantitative diagnostic path for large (100M+) particle runs.
 *
 * To add a new statistic:
 *   1. Implement PetscErrorCode ComputeXxx(UserCtx*, const char*, PetscInt) here.
 *   2. Declare it in include/particle_statistics.h.
 *   3. Add one else-if in GlobalStatisticsPipeline() in src/postprocessor.c.
 */

#include "particle_statistics.h"
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

static const char *PARTICLE_MSD_CSV_HEADER =
    "step,t,N,MSD_x,MSD_y,MSD_z,MSD_total,"
    "r_rms_meas,r_rms_theory,rel_err_pct,"
    "com_x,com_y,com_z,"
    "frac_1sigma_pct,frac_2sigma_pct,frac_3sigma_pct\n";

/**
 * @brief Return whether a CSV row belongs to the requested timestep.
 * @param[in] line CSV row text to inspect.
 * @param[in] ti Requested timestep id.
 * @return PETSC_TRUE when the row begins with the requested step id.
 */
static PetscBool MSDCSVLineMatchesStep(const char *line, PetscInt ti)
{
    char *endptr = NULL;
    long parsed_step;

    if (!line) return PETSC_FALSE;
    parsed_step = strtol(line, &endptr, 10);
    if (endptr == line) return PETSC_FALSE;
    while (*endptr == ' ' || *endptr == '\t') endptr++;
    if (*endptr != ',') return PETSC_FALSE;
    return (PetscBool)(parsed_step == (long)ti ? PETSC_TRUE : PETSC_FALSE);
}

/**
 * @brief Rewrite the MSD CSV so the requested timestep appears exactly once.
 * @param[in] csv_path Target CSV path.
 * @param[in] ti Timestep being written.
 * @param[in] row_line Fully formatted CSV row for the timestep.
 * @return Petsc error code.
 */
static PetscErrorCode RewriteParticleMSDCSV(const char *csv_path, PetscInt ti, const char *row_line)
{
    char tmp_path[PETSC_MAX_PATH_LEN];
    FILE *src = NULL;
    FILE *dst = NULL;
    char line[4096];
    PetscBool wrote_header = PETSC_FALSE;
    PetscBool inserted_row = PETSC_FALSE;

    PetscFunctionBeginUser;
    PetscCall(PetscSNPrintf(tmp_path, sizeof(tmp_path), "%s.tmp", csv_path));

    dst = fopen(tmp_path, "w");
    if (!dst) {
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "ComputeParticleMSD: could not open '%s' for writing.\n", tmp_path);
        PetscFunctionReturn(0);
    }

    src = fopen(csv_path, "r");
    if (src) {
        while (fgets(line, sizeof(line), src)) {
            if (!wrote_header) {
                if (strncmp(line, "step,", 5) == 0) {
                    fputs(line, dst);
                    wrote_header = PETSC_TRUE;
                    continue;
                }
                fputs(PARTICLE_MSD_CSV_HEADER, dst);
                wrote_header = PETSC_TRUE;
            }

            if (MSDCSVLineMatchesStep(line, ti)) {
                if (!inserted_row) {
                    fputs(row_line, dst);
                    inserted_row = PETSC_TRUE;
                }
                continue;
            }
            fputs(line, dst);
        }
        fclose(src);
        src = NULL;
    }

    if (!wrote_header) {
        fputs(PARTICLE_MSD_CSV_HEADER, dst);
    }
    if (!inserted_row) {
        fputs(row_line, dst);
    }

    if (fclose(dst) != 0) {
        remove(tmp_path);
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "ComputeParticleMSD: could not finalize '%s'.\n", tmp_path);
        PetscFunctionReturn(0);
    }

    if (rename(tmp_path, csv_path) != 0) {
        remove(tmp_path);
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "ComputeParticleMSD: could not replace '%s'.\n", csv_path);
    }
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeParticleMSD"
/**
 * @brief Internal helper implementation: `ComputeParticleMSD()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ComputeParticleMSD(UserCtx *user, const char *stats_prefix, PetscInt ti)
{
    PetscErrorCode ierr;
    SimCtx        *simCtx = user->simCtx;
    PetscMPIInt    rank;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    /* ------------------------------------------------------------------ *
     * Physics parameters                                                   *
     * ------------------------------------------------------------------ */
    const PetscReal Re  = (simCtx->ren > 0.0) ? simCtx->ren : 1.0;
    const PetscReal Sc  = (simCtx->schmidt_number > 0.0) ? simCtx->schmidt_number : 1.0;
    const PetscReal D   = 1.0 / (Re * Sc);
    const PetscReal t   = (PetscReal)ti * simCtx->dt;
    const PetscReal x0  = simCtx->psrc_x;
    const PetscReal y0  = simCtx->psrc_y;
    const PetscReal z0  = simCtx->psrc_z;

    /* ------------------------------------------------------------------ *
     * Pass 1: local accumulation                                           *
     * ------------------------------------------------------------------ */
    PetscInt       n_local;
    const PetscReal (*pos_arr)[3];

    ierr = DMSwarmGetLocalSize(user->swarm, &n_local); CHKERRQ(ierr);
    ierr = DMSwarmGetField(user->swarm, "position", NULL, NULL, (void**)&pos_arr); CHKERRQ(ierr);

    PetscReal local_sq_x = 0.0, local_sq_y = 0.0, local_sq_z = 0.0;
    PetscReal local_sx   = 0.0, local_sy   = 0.0, local_sz   = 0.0;

    for (PetscInt p = 0; p < n_local; p++) {
        const PetscReal dx = pos_arr[p][0] - x0;
        const PetscReal dy = pos_arr[p][1] - y0;
        const PetscReal dz = pos_arr[p][2] - z0;
        local_sq_x += dx * dx;
        local_sq_y += dy * dy;
        local_sq_z += dz * dz;
        local_sx   += dx;
        local_sy   += dy;
        local_sz   += dz;
    }

    ierr = DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void**)&pos_arr); CHKERRQ(ierr);

    /* ------------------------------------------------------------------ *
     * MPI reduction — 7 doubles + count                                   *
     * ------------------------------------------------------------------ */
    PetscReal local_buf[7] = { local_sq_x, local_sq_y, local_sq_z,
                                local_sx,   local_sy,   local_sz,
                                (PetscReal)n_local };
    PetscReal global_buf[7] = {0};
    ierr = MPI_Allreduce(local_buf, global_buf, 7, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRMPI(ierr);

    const PetscReal N_total  = global_buf[6];
    if (N_total < 1.0) PetscFunctionReturn(0); /* no particles */

    const PetscReal MSD_x    = global_buf[0] / N_total;
    const PetscReal MSD_y    = global_buf[1] / N_total;
    const PetscReal MSD_z    = global_buf[2] / N_total;
    const PetscReal com_x    = global_buf[3] / N_total;
    const PetscReal com_y    = global_buf[4] / N_total;
    const PetscReal com_z    = global_buf[5] / N_total;

    const PetscReal MSD_total    = MSD_x + MSD_y + MSD_z;
    const PetscReal r_rms_meas   = PetscSqrtReal(MSD_total);
    const PetscReal r_theory     = (t > 1e-300) ? PetscSqrtReal(6.0 * D * t) : 0.0;
    const PetscReal rel_err_pct  = (r_theory > 1e-12)
                                   ? PetscAbsReal(r_rms_meas - r_theory) / r_theory * 100.0
                                   : 0.0;

    /* ------------------------------------------------------------------ *
     * Pass 2: fraction inside σ-shells                                    *
     * ------------------------------------------------------------------ */
    ierr = DMSwarmGetField(user->swarm, "position", NULL, NULL, (void**)&pos_arr); CHKERRQ(ierr);

    PetscInt local_n1 = 0, local_n2 = 0, local_n3 = 0;
    const PetscReal r1 = r_theory;
    const PetscReal r2 = 2.0 * r_theory;
    const PetscReal r3 = 3.0 * r_theory;

    for (PetscInt p = 0; p < n_local; p++) {
        const PetscReal dx = pos_arr[p][0] - x0;
        const PetscReal dy = pos_arr[p][1] - y0;
        const PetscReal dz = pos_arr[p][2] - z0;
        const PetscReal r  = PetscSqrtReal(dx*dx + dy*dy + dz*dz);
        if (r < r1) local_n1++;
        if (r < r2) local_n2++;
        if (r < r3) local_n3++;
    }

    ierr = DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (void**)&pos_arr); CHKERRQ(ierr);

    PetscReal local_counts[3] = { (PetscReal)local_n1,
                                   (PetscReal)local_n2,
                                   (PetscReal)local_n3 };
    PetscReal global_counts[3] = {0};
    ierr = MPI_Allreduce(local_counts, global_counts, 3, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRMPI(ierr);

    const PetscReal frac_1s = global_counts[0] / N_total * 100.0;
    const PetscReal frac_2s = global_counts[1] / N_total * 100.0;
    const PetscReal frac_3s = global_counts[2] / N_total * 100.0;

    /* ------------------------------------------------------------------ *
     * Output: CSV + LOG_INFO (rank 0 only)                                *
     * ------------------------------------------------------------------ */
    if (rank == 0) {
        char csv_path[PETSC_MAX_PATH_LEN];
        char row_line[1024];
        PetscSNPrintf(csv_path, sizeof(csv_path), "%s_msd.csv", stats_prefix);
        PetscSNPrintf(
            row_line,
            sizeof(row_line),
            "%d,%.6e,%.0f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.4f,%.6e,%.6e,%.6e,%.2f,%.2f,%.2f\n",
            (int)ti, t, N_total,
            MSD_x, MSD_y, MSD_z, MSD_total,
            r_rms_meas, r_theory, rel_err_pct,
            com_x, com_y, com_z,
            frac_1s, frac_2s, frac_3s
        );
        PetscCall(RewriteParticleMSDCSV(csv_path, ti, row_line));

        LOG_ALLOW(GLOBAL, LOG_INFO,
                  "[MSD ti=%d t=%.4f] total=%.4e theory=%.4e err=%.2f%% | "
                  "fracs: %.1f/%.1f/%.1f%% | COM: (%.2e,%.2e,%.2e)\n",
                  (int)ti, t, MSD_total, r_theory * r_theory, rel_err_pct,
                  frac_1s, frac_2s, frac_3s,
                  com_x, com_y, com_z);
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
