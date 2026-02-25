/**
 * @file particle_statistics.c
 * @brief Global statistical reduction kernels for the Statistics Pipeline.
 *
 * Each kernel in this file:
 *  - Performs an MPI_Allreduce over all particles across all ranks.
 *  - Appends one CSV row per call to its own output file.
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

#undef __FUNCT__
#define __FUNCT__ "ComputeParticleMSD"
/**
 * @brief Computes the mean-squared displacement (MSD) of the particle cloud.
 *
 * Two-pass algorithm:
 *   Pass 1 — accumulate per-component squared displacements and centre-of-mass sums.
 *   MPI reduction — allreduce 7 values (6 sums + count) to get global MSD and COM.
 *   Compute scalars — MSD_total, r_rms_meas, r_rms_theory, rel_err_pct.
 *   Pass 2 — count particles inside 1σ/2σ/3σ theoretical shells; allreduce counts.
 *   Output (rank 0) — append one row to {stats_prefix}_msd.csv; LOG_INFO summary.
 *
 * Physics:  D = 1/(Re·Sc),  t = ti·dt,  r_theory = sqrt(6·D·t).
 * Expected fractions inside r_theory / 2·r_theory / 3·r_theory ≈ 61%/97%/99.9%
 * for a 3D Maxwell-Boltzmann distribution.
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
    ierr = DMSwarmGetField(user->swarm, "position", NULL, NULL, (const void**)&pos_arr); CHKERRQ(ierr);

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

    ierr = DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (const void**)&pos_arr); CHKERRQ(ierr);

    /* ------------------------------------------------------------------ *
     * MPI reduction — 7 doubles + count                                   *
     * ------------------------------------------------------------------ */
    PetscReal local_buf[7] = { local_sq_x, local_sq_y, local_sq_z,
                                local_sx,   local_sy,   local_sz,
                                (PetscReal)n_local };
    PetscReal global_buf[7] = {0};
    MPI_Allreduce(local_buf, global_buf, 7, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

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
    ierr = DMSwarmGetField(user->swarm, "position", NULL, NULL, (const void**)&pos_arr); CHKERRQ(ierr);

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

    ierr = DMSwarmRestoreField(user->swarm, "position", NULL, NULL, (const void**)&pos_arr); CHKERRQ(ierr);

    PetscReal local_counts[3] = { (PetscReal)local_n1,
                                   (PetscReal)local_n2,
                                   (PetscReal)local_n3 };
    PetscReal global_counts[3] = {0};
    MPI_Allreduce(local_counts, global_counts, 3, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

    const PetscReal frac_1s = global_counts[0] / N_total * 100.0;
    const PetscReal frac_2s = global_counts[1] / N_total * 100.0;
    const PetscReal frac_3s = global_counts[2] / N_total * 100.0;

    /* ------------------------------------------------------------------ *
     * Output: CSV + LOG_INFO (rank 0 only)                                *
     * ------------------------------------------------------------------ */
    if (rank == 0) {
        char csv_path[PETSC_MAX_PATH_LEN];
        PetscSNPrintf(csv_path, sizeof(csv_path), "%s_msd.csv", stats_prefix);

        FILE *f = fopen(csv_path, "a");
        if (!f) {
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "ComputeParticleMSD: could not open '%s' for writing.\n", csv_path);
        } else {
            /* Write header if file is empty */
            if (ftell(f) == 0) {
                fprintf(f, "step,t,N,MSD_x,MSD_y,MSD_z,MSD_total,"
                        "r_rms_meas,r_rms_theory,rel_err_pct,"
                        "com_x,com_y,com_z,"
                        "frac_1sigma_pct,frac_2sigma_pct,frac_3sigma_pct\n");
            }
            fprintf(f, "%d,%.6e,%.0f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.4f,%.6e,%.6e,%.6e,%.2f,%.2f,%.2f\n",
                    (int)ti, t, N_total,
                    MSD_x, MSD_y, MSD_z, MSD_total,
                    r_rms_meas, r_theory, rel_err_pct,
                    com_x, com_y, com_z,
                    frac_1s, frac_2s, frac_3s);
            fclose(f);
        }

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
