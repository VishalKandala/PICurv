#include "postprocessing_kernels.h"
#include "statistics_accumulator.h"
#include "statistics_window.h"
#include "particle_field_catalog.h"

// =========== Dimensionalization Kernels ========================
#undef __FUNCT__
#define __FUNCT__ "DimensionalizeField"
/**
 * @brief Implementation of \ref DimensionalizeField().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/postprocessing_kernels.h`.
 * @see DimensionalizeField()
 */
PetscErrorCode DimensionalizeField(UserCtx *user, const char *field_name)
{
    PetscErrorCode ierr;
    SimCtx         *simCtx = user->simCtx;
    Vec            target_vec = NULL;
    PetscReal      scale_factor = 1.0;
    char           field_type[64] = "Unknown";
    PetscBool      is_swarm_field = PETSC_FALSE; // Flag for special swarm handling
    const char     *swarm_field_name = NULL;    // Name of the field within the swarm

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    if (!user) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx is NULL.");
    if (!field_name) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "field_name is NULL.");

    // --- 1. Identify the target Vec and the correct scaling factor ---
    if (strcasecmp(field_name, "Ucat") == 0) {
        target_vec = user->Ucat;
        scale_factor = simCtx->scaling.U_ref;
        strcpy(field_type, "Cartesian Velocity (L/T)");
    } else if (strcasecmp(field_name, "Ucont") == 0) {
        target_vec = user->Ucont;
        scale_factor = simCtx->scaling.U_ref * simCtx->scaling.L_ref * simCtx->scaling.L_ref;
        strcpy(field_type, "Contravariant Volume Flux (L^3/T)");
    } else if (strcasecmp(field_name, "P") == 0) {
        target_vec = user->P;
        scale_factor = simCtx->scaling.P_ref;
        strcpy(field_type, "Pressure (M L^-1 T^-2)");
    } else if (strcasecmp(field_name, "Coordinates") == 0) {
        ierr = DMGetCoordinates(user->da, &target_vec); CHKERRQ(ierr);
        scale_factor = simCtx->scaling.L_ref;
        strcpy(field_type, "Grid Coordinates (L)");
    } else if (strcasecmp(field_name, "ParticlePosition") == 0) {
        is_swarm_field = PETSC_TRUE;
        swarm_field_name = "position";
        scale_factor = simCtx->scaling.L_ref;
        strcpy(field_type, "Particle Position (L)");
    } else if (strcasecmp(field_name, "ParticleVelocity") == 0) {
        is_swarm_field = PETSC_TRUE;
        swarm_field_name = "velocity";
        scale_factor = simCtx->scaling.U_ref;
        strcpy(field_type, "Particle Velocity (L/T)");
    } else {
        LOG(GLOBAL, LOG_WARNING, "DimensionalizeField: Unknown or unhandled field_name '%s'. Field will not be scaled.\n", field_name);
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }
    
    // --- 2. Check for trivial scaling ---
    if (PetscAbsReal(scale_factor - 1.0) < PETSC_MACHINE_EPSILON) {
        LOG(GLOBAL, LOG_DEBUG, "DimensionalizeField: Scaling factor for '%s' is 1.0. Skipping operation.\n", field_name);
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    // --- 3. Perform the in-place scaling operation ---
    LOG(GLOBAL, LOG_INFO, "Scaling '%s' field (%s) by factor %.4e.\n", field_name, field_type, scale_factor);

    if (is_swarm_field) {
        // Special handling for DMSwarm fields
        ierr = DMSwarmCreateGlobalVectorFromField(user->swarm, swarm_field_name, &target_vec); CHKERRQ(ierr);
        ierr = VecScale(target_vec, scale_factor); CHKERRQ(ierr);
        ierr = DMSwarmDestroyGlobalVectorFromField(user->swarm, swarm_field_name, &target_vec); CHKERRQ(ierr);
    } else {
        // Standard handling for PETSc Vecs
        if (target_vec) {
            ierr = VecScale(target_vec, scale_factor); CHKERRQ(ierr);
        } else {
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "Target vector for field '%s' was not found or is NULL.", field_name);
        }
    }
    
    // --- 4. Post-scaling updates for special cases ---
    if (strcasecmp(field_name, "Coordinates") == 0) {
        ierr = UpdateLocalGhosts(user, FIELD_ID_COORDINATES); CHKERRQ(ierr);
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "DimensionalizeAllLoadedFields"
/**
 * @brief Internal helper implementation: `DimensionalizeAllLoadedFields()`.
 * @details Local to this translation unit.
 */
PetscErrorCode DimensionalizeAllLoadedFields(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx         *simCtx = user->simCtx;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    LOG(GLOBAL, LOG_INFO, "--- Converting all loaded fields to dimensional units ---\n");

    // Scale the grid itself first
    ierr = DimensionalizeField(user, "Coordinates"); CHKERRQ(ierr);

    // Scale primary fluid fields
    ierr = DimensionalizeField(user, "Ucat"); CHKERRQ(ierr);
    ierr = DimensionalizeField(user, "Ucont"); CHKERRQ(ierr);
    ierr = DimensionalizeField(user, "P"); CHKERRQ(ierr);

    // If particles are present, scale their fields
    if (simCtx->np > 0 && user->swarm) {
        ierr = DimensionalizeField(user, "ParticlePosition"); CHKERRQ(ierr);
        ierr = DimensionalizeField(user, "ParticleVelocity"); CHKERRQ(ierr);
    }

    LOG(GLOBAL, LOG_INFO, "--- Field dimensionalization complete ---\n");

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

//============ Post-Processing Kernels ===========================

#undef __FUNCT__
#define __FUNCT__ "ExtendToLayoutBoundary"
/**
 * @brief Implementation of \ref ExtendToLayoutBoundary().
 * @details Full API contract is documented with the header declaration in
 *          `include/postprocessing_kernels.h`.
 * @see ExtendToLayoutBoundary()
 */
PetscErrorCode ExtendToLayoutBoundary(UserCtx *user, Vec global, PetscInt components)
{
    SimCtx        *simCtx = NULL;
    DM             dm = NULL;
    Vec            local = NULL;
    PetscInt       periodic[3];

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCheck(user != NULL && global != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Context and field are required.");
    PetscCheck(components == 1 || components == 3, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "ExtendToLayoutBoundary handles 1- or 3-component fields; got %d.", (int)components);
    simCtx = user->simCtx;
    dm = (components == 1) ? user->da : user->fda;

    periodic[0] = simCtx->i_periodic;
    periodic[1] = simCtx->j_periodic;
    periodic[2] = simCtx->k_periodic;
    if (!periodic[0] && !periodic[1] && !periodic[2]) {
        /* Nothing is defined on a non-periodic layout boundary, so nothing is written.
         * The boundary node keeps the value the producing kernel left there, and the
         * spatial reductions that matter scientifically never read it. */
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    PetscCall(DMGetLocalVector(dm, &local));
    /* One pass per periodic direction, each preceded by its own scatter. The order
     * matters: after the i pass the i layout boundary carries data, so the j pass
     * reading across it fills edges correctly, and the k pass then fills corners. */
    for (PetscInt dir = 0; dir < 3; ++dir) {
        DMDALocalInfo info;
        const PetscReal ****source = NULL;
        PetscReal      ****target = NULL;
        PetscInt        extent = 0;

        if (!periodic[dir]) continue;

        PetscCall(DMGlobalToLocalBegin(dm, global, INSERT_VALUES, local));
        PetscCall(DMGlobalToLocalEnd(dm, global, INSERT_VALUES, local));
        PetscCall(DMDAGetLocalInfo(dm, &info));
        extent = (dir == 0) ? info.mx : ((dir == 1) ? info.my : info.mz);

        PetscCall(DMDAVecGetArrayDOFRead(dm, local, &source));
        PetscCall(DMDAVecGetArrayDOF(dm, global, &target));
        for (PetscInt k = info.zs; k < info.zs + info.zm; ++k) {
            for (PetscInt j = info.ys; j < info.ys + info.ym; ++j) {
                for (PetscInt i = info.xs; i < info.xs + info.xm; ++i) {
                    const PetscInt index = (dir == 0) ? i : ((dir == 1) ? j : k);
                    PetscInt       ss = i, sj = j, sk = k;

                    /* The layout wraps one plane inside the DMDA extent: the low dummy
                     * plane repeats the last physical plane, and the high dummy plane
                     * repeats the first. PETSc's periodic ghosts supply both sources
                     * even when they live on another rank. */
                    if (index == 0) {
                        if (dir == 0) ss = extent - 2;
                        else if (dir == 1) sj = extent - 2;
                        else sk = extent - 2;
                    } else if (index == extent - 1) {
                        if (dir == 0) ss = extent + 1;
                        else if (dir == 1) sj = extent + 1;
                        else sk = extent + 1;
                    } else {
                        continue;
                    }
                    for (PetscInt c = 0; c < components; ++c) {
                        target[k][j][i][c] = source[sk][sj][ss][c];
                    }
                }
            }
        }
        PetscCall(DMDAVecRestoreArrayDOF(dm, global, &target));
        PetscCall(DMDAVecRestoreArrayDOFRead(dm, local, &source));
    }
    PetscCall(DMRestoreLocalVector(dm, &local));

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "-> KERNEL: Extended %d-component field across the periodic layout boundary "
              "(i=%d, j=%d, k=%d).\n", (int)components,
              (int)periodic[0], (int)periodic[1], (int)periodic[2]);
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeNodalAverage"
/**
 * @brief Implementation of \ref ComputeNodalAverage().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/postprocessing_kernels.h`.
 * @see ComputeNodalAverage()
 */
PetscErrorCode ComputeNodalAverage(UserCtx* user, const char* in_field_name, const char* out_field_name)
{
    PetscErrorCode ierr;
    FieldId        in_field_id;
    Vec            in_vec_local = NULL, out_vec_global = NULL;
    DM             dm_in = NULL, dm_out = NULL;
    PetscInt       dof = 0;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(GLOBAL, LOG_INFO, "-> KERNEL: Running ComputeNodalAverage on '%s' -> '%s'.\n", in_field_name, out_field_name);

    // --- 1. Map string names to PETSc objects ---
    if (strcasecmp(in_field_name, "P") == 0)             { in_vec_local = user->lP;         dm_in = user->da;   dof = 1; }
    else if (strcasecmp(in_field_name, "Ucat") == 0)    { in_vec_local = user->lUcat;      dm_in = user->fda;  dof = 3; }
    else if (strcasecmp(in_field_name, "Psi") == 0)     { in_vec_local = user->lPsi;       dm_in = user->da;   dof = 1; }
    /* The staging pair carries derived statistics, which are config-counted and so
     * cannot be named by a compile-time member of their own. */
    else if (strcasecmp(in_field_name, "PostScalar") == 0) { in_vec_local = user->lPostScalar; dm_in = user->da;  dof = 1; }
    else if (strcasecmp(in_field_name, "PostVector") == 0) { in_vec_local = user->lPostVector; dm_in = user->fda; dof = 3; }
    // ... (add other fields as needed) ...
    else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unknown input field name for nodal averaging: %s", in_field_name);

    if (strcasecmp(out_field_name, "P_nodal") == 0)      { out_vec_global = user->P_nodal;    dm_out = user->da; }
    else if (strcasecmp(out_field_name, "Ucat_nodal") == 0) { out_vec_global = user->Ucat_nodal; dm_out = user->fda; }
    else if (strcasecmp(out_field_name, "Psi_nodal") == 0)   { out_vec_global = user->Psi_nodal;  dm_out = user->da; }
    else if (strcasecmp(out_field_name, "PostScalarNodal") == 0) { out_vec_global = user->PostScalarNodal; dm_out = user->da; }
    else if (strcasecmp(out_field_name, "PostVectorNodal") == 0) { out_vec_global = user->PostVectorNodal; dm_out = user->fda; }
    // ... (add other fields as needed) ...
    else SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Unknown output field name for nodal averaging: %s", out_field_name);

    // --- 2. Ensure Input Data Ghosts are Up-to-Date ---
    ierr = FieldIdFromName(in_field_name, &in_field_id); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, in_field_id); CHKERRQ(ierr);

    // --- 3. Get DMDA info and array pointers ---
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(dm_out, &info); CHKERRQ(ierr);
    /* Every owned output point is valid except the global high layout plane.
     * A rank interface is not a boundary: its +1 source value is in the halo. */
    const PetscInt i_end = PetscMin(info.xs + info.xm, info.mx - 1);
    const PetscInt j_end = PetscMin(info.ys + info.ym, info.my - 1);
    const PetscInt k_end = PetscMin(info.zs + info.zm, info.mz - 1);

    if (dof == 1) { // --- Scalar Field Averaging ---
        const PetscReal ***l_in_arr;
        PetscReal       ***g_out_arr;
        ierr = DMDAVecGetArrayRead(dm_in,in_vec_local, (void*)&l_in_arr); CHKERRQ(ierr);
        ierr = DMDAVecGetArray(dm_out,out_vec_global, (void*)&g_out_arr); CHKERRQ(ierr);

        // Loop over the output NODE locations. The loop bounds match the required
        // size of the final subsampled grid.
        for (PetscInt k = info.zs; k < k_end; k++) {
            for (PetscInt j = info.ys; j < j_end; j++) {
                for (PetscInt i = info.xs; i < i_end; i++) {
                    g_out_arr[k][j][i] = 0.125 * (l_in_arr[k][j][i]     + l_in_arr[k][j][i+1] +
                                                  l_in_arr[k][j+1][i]   + l_in_arr[k][j+1][i+1] +
                                                  l_in_arr[k+1][j][i]   + l_in_arr[k+1][j][i+1] +
                                                  l_in_arr[k+1][j+1][i] + l_in_arr[k+1][j+1][i+1]);
                }
            }
        }
        ierr = DMDAVecRestoreArrayRead(dm_in,in_vec_local, (void*)&l_in_arr); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(dm_out,out_vec_global, (void*)&g_out_arr); CHKERRQ(ierr);

    } else if (dof == 3) { // --- Vector Field Averaging ---
        const Cmpnts ***l_in_arr;
        Cmpnts       ***g_out_arr;
        ierr = DMDAVecGetArrayRead(dm_in,in_vec_local, (void*)&l_in_arr); CHKERRQ(ierr);
        ierr = DMDAVecGetArray(dm_out,out_vec_global, (void*)&g_out_arr); CHKERRQ(ierr);

        for (PetscInt k = info.zs; k < k_end; k++) {
            for (PetscInt j = info.ys; j < j_end; j++) {
                for (PetscInt i = info.xs; i < i_end; i++) {
                    g_out_arr[k][j][i].x = 0.125 * (l_in_arr[k][j][i].x + l_in_arr[k][j][i+1].x +
                                                    l_in_arr[k][j+1][i].x + l_in_arr[k][j+1][i+1].x +
                                                    l_in_arr[k+1][j][i].x + l_in_arr[k+1][j][i+1].x +
                                                    l_in_arr[k+1][j+1][i].x + l_in_arr[k+1][j+1][i+1].x);

                    g_out_arr[k][j][i].y = 0.125 * (l_in_arr[k][j][i].y + l_in_arr[k][j][i+1].y +
                                                    l_in_arr[k][j+1][i].y + l_in_arr[k][j+1][i+1].y +
                                                    l_in_arr[k+1][j][i].y + l_in_arr[k+1][j][i+1].y +
                                                    l_in_arr[k+1][j+1][i].y + l_in_arr[k+1][j+1][i+1].y);

                    g_out_arr[k][j][i].z = 0.125 * (l_in_arr[k][j][i].z + l_in_arr[k][j][i+1].z +
                                                    l_in_arr[k][j+1][i].z + l_in_arr[k][j+1][i+1].z +
                                                    l_in_arr[k+1][j][i].z + l_in_arr[k+1][j][i+1].z +
                                                    l_in_arr[k+1][j+1][i].z + l_in_arr[k+1][j+1][i+1].z);
                }
            }
        }
        ierr = DMDAVecRestoreArrayRead(dm_in,in_vec_local, (void*)&l_in_arr); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(dm_out,out_vec_global, (void*)&g_out_arr); CHKERRQ(ierr);
    }
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ComputeWindowStatisticNodal"
/**
 * @brief Implementation of \ref ComputeWindowStatisticNodal().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/postprocessing_kernels.h`.
 * @see ComputeWindowStatisticNodal()
 */
PetscErrorCode ComputeWindowStatisticNodal(UserCtx *user, PetscInt window_index,
                                           const char *outputs, PetscInt output_index,
                                           char *out_name, size_t name_size,
                                           Vec *out_vec, PetscInt *out_components)
{
    PetscErrorCode ierr;
    SimCtx *simCtx = NULL;
    PicurvDerivedField derived;
    const PicurvWindowDefinition *definition = NULL;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCheck(user != NULL && out_name != NULL && out_vec != NULL && out_components != NULL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Context and outputs are required.");
    simCtx = user->simCtx;
    PetscCheck(FieldStatisticsIsActive(simCtx) && user->fieldStatisticsStorage != NULL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "No accumulated window state exists to derive.");
    definition = &simCtx->fieldStatisticsWindows[window_index].definition;

    /* Derive into the cell-centred staging pair, then reach the nodal path by
     * catalogued name. A per-window accumulator has no compile-time offset, so
     * staging is what lets the shared ghost and nodal kernels address it at all. */
    ierr = PicurvWindowDerive(user, definition, &user->fieldStatisticsStorage[window_index],
                              outputs, output_index, user->PostScalar, user->PostVector,
                              &derived); CHKERRQ(ierr);
    /* The derivation writes the physical interior only, so the layout boundary would
     * otherwise reach the nodal average as structural zeros. Extend it before the
     * ghost update, which then carries the written values outward. */
    if (derived.components == 1) {
        ierr = ExtendToLayoutBoundary(user, user->PostScalar, 1); CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, FIELD_ID_POST_SCALAR); CHKERRQ(ierr);
        ierr = ComputeNodalAverage(user, "PostScalar", "PostScalarNodal"); CHKERRQ(ierr);
        *out_vec = user->PostScalarNodal;
    } else {
        ierr = ExtendToLayoutBoundary(user, user->PostVector, 3); CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, FIELD_ID_POST_VECTOR); CHKERRQ(ierr);
        ierr = ComputeNodalAverage(user, "PostVector", "PostVectorNodal"); CHKERRQ(ierr);
        *out_vec = user->PostVectorNodal;
    }
    *out_components = derived.components;
    ierr = PetscStrncpy(out_name, derived.name, name_size); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "-> KERNEL: Derived '%s' (%d component(s)).\n",
              out_name, (int)derived.components);
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeWindowStatisticsSummary"
/**
 * @brief Implementation of \ref ComputeWindowStatisticsSummary().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/postprocessing_kernels.h`.
 * @see ComputeWindowStatisticsSummary()
 */
PetscErrorCode ComputeWindowStatisticsSummary(UserCtx *user, PetscInt window_index,
                                              const char *output_prefix, PetscInt ti)
{
    PetscErrorCode ierr;
    SimCtx *simCtx = NULL;
    const PicurvWindow *window = NULL;
    const PicurvWindowStorage *storage = NULL;
    PetscReal lowest = 1.0, highest = 0.0;
    PetscReal mean_tke = 0.0;
    PetscBool has_tke = PETSC_FALSE;
    PetscInt derived_count = 0;
    char path[PETSC_MAX_PATH_LEN];

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCheck(user != NULL && output_prefix != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Context and output prefix are required.");
    simCtx = user->simCtx;
    PetscCheck(FieldStatisticsIsActive(simCtx) && user->fieldStatisticsStorage != NULL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "No accumulated window state exists to summarize.");
    window = &simCtx->fieldStatisticsWindows[window_index];
    storage = &user->fieldStatisticsStorage[window_index];

    ierr = PicurvWindowValidFractionRange(user, &window->definition, storage,
                                          window->sample_count, &lowest, &highest); CHKERRQ(ierr);

    /* A single domain number is what makes the row plottable against step. It is
     * produced by the same derivation the field output uses, so the scalar and the
     * field can never disagree about what the window holds. */
    ierr = PicurvWindowDerivedCount(&window->definition, storage, "tke", &derived_count); CHKERRQ(ierr);
    if (derived_count > 0) {
        PicurvDerivedField derived;

        ierr = PicurvWindowDerive(user, &window->definition, storage, "tke", 0,
                                  user->PostScalar, user->PostVector, &derived); CHKERRQ(ierr);
        /* Averaged over the fluid cells the window actually sampled. A whole-vector
         * mean would divide by boundary and dummy entries the derivation never
         * writes, scaling the answer down by the fraction it never covered. */
        ierr = PicurvWindowSpatialMean(user, &window->definition, storage,
                                       user->PostScalar, &mean_tke); CHKERRQ(ierr);
        has_tke = PETSC_TRUE;
    }

    if (simCtx->rank == 0) {
        FILE *csv = NULL;
        PetscBool exists = PETSC_FALSE;

        ierr = PetscSNPrintf(path, sizeof(path), "%s_statistics_%s.csv",
                             output_prefix, window->definition.name); CHKERRQ(ierr);
        ierr = PetscTestFile(path, 'r', &exists); CHKERRQ(ierr);
        csv = fopen(path, exists ? "a" : "w");
        PetscCheck(csv != NULL, PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN,
                   "Unable to open statistics summary '%s'.", path);
        if (!exists) {
            fprintf(csv, "step,state,samples,total_weight,represented_time,"
                         "valid_fraction_min,valid_fraction_max,mean_tke\n");
        }
        fprintf(csv, "%" PetscInt_FMT ",%s,%d,%.10e,%.10e,%.6f,%.6f,",
                ti, PicurvWindowStateName(window->state), window->sample_count,
                (double)window->total_weight, (double)window->represented_time,
                (double)lowest, (double)highest);
        if (has_tke) fprintf(csv, "%.10e\n", (double)mean_tke);
        else fprintf(csv, "\n");
        PetscCheck(fclose(csv) == 0, PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
                   "Unable to close statistics summary '%s'.", path);
    }
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "-> KERNEL: Summarized window '%s' at step %" PetscInt_FMT ".\n",
              window->definition.name, ti);
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeQCriterion"
/**
 * @brief Implementation of \ref ComputeQCriterion().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/postprocessing_kernels.h`.
 * @see ComputeQCriterion()
 */
PetscErrorCode ComputeQCriterion(UserCtx* user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    const Cmpnts   ***lucat, ***lcsi, ***leta, ***lzet;
    const PetscReal***laj, ***lnvert;
    PetscReal      ***gq;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(GLOBAL, LOG_INFO, "-> KERNEL: Running ComputeQCriterion.\n");

    // --- 1. Ensure all required ghost values are up-to-date ---
    ierr = UpdateLocalGhosts(user, FIELD_ID_UCAT);  CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_CSI);   CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_ETA);   CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_ZET);   CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_AJ);    CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_NVERT); CHKERRQ(ierr);

    // --- 2. Get DMDA info and array pointers ---
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    
    ierr = DMDAVecGetArrayRead(user->fda, user->lUcat,  (void*)&lucat);   CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi,   (void*)&lcsi);    CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta,   (void*)&leta);    CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet,   (void*)&lzet);    CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da,  user->lAj,    (void*)&laj);     CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da,  user->lNvert, (void*)&lnvert);  CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da,  user->Qcrit, (void*)&gq);       CHKERRQ(ierr);

    // --- 3. Define Loop Bounds for INTERIOR Cells ---
    PetscInt i_start = (info.xs == 0) ? 1 : info.xs;
    PetscInt i_end   = (info.xs + info.xm == info.mx) ? info.mx - 1 : info.xs + info.xm;
    PetscInt j_start = (info.ys == 0) ? 1 : info.ys;
    PetscInt j_end   = (info.ys + info.ym == info.my) ? info.my - 1 : info.ys + info.ym;
    PetscInt k_start = (info.zs == 0) ? 1 : info.zs;
    PetscInt k_end   = (info.zs + info.zm == info.mz) ? info.mz - 1 : info.zs + info.zm;

    // --- 4. Main Computation Loop ---
    for (PetscInt k = k_start; k < k_end; k++) {
        for (PetscInt j = j_start; j < j_end; j++) {
            for (PetscInt i = i_start; i < i_end; i++) {
                
                // Calculate velocity derivatives in computational space (central differences)
                PetscReal uc = 0.5 * (lucat[k][j][i+1].x - lucat[k][j][i-1].x);
                PetscReal vc = 0.5 * (lucat[k][j][i+1].y - lucat[k][j][i-1].y);
                PetscReal wc = 0.5 * (lucat[k][j][i+1].z - lucat[k][j][i-1].z);

                PetscReal ue = 0.5 * (lucat[k][j+1][i].x - lucat[k][j-1][i].x);
                PetscReal ve = 0.5 * (lucat[k][j+1][i].y - lucat[k][j-1][i].y);
                PetscReal we = 0.5 * (lucat[k][j+1][i].z - lucat[k][j-1][i].z);

                PetscReal uz = 0.5 * (lucat[k+1][j][i].x - lucat[k-1][j][i].x);
                PetscReal vz = 0.5 * (lucat[k+1][j][i].y - lucat[k-1][j][i].y);
                PetscReal wz = 0.5 * (lucat[k+1][j][i].z - lucat[k-1][j][i].z);

                // Average metrics to the cell center
                PetscReal csi1 = 0.5 * (lcsi[k][j][i].x + lcsi[k][j][i-1].x) * laj[k][j][i];
                PetscReal csi2 = 0.5 * (lcsi[k][j][i].y + lcsi[k][j][i-1].y) * laj[k][j][i];
                PetscReal csi3 = 0.5 * (lcsi[k][j][i].z + lcsi[k][j][i-1].z) * laj[k][j][i];
                
                PetscReal eta1 = 0.5 * (leta[k][j][i].x + leta[k][j-1][i].x) * laj[k][j][i];
                PetscReal eta2 = 0.5 * (leta[k][j][i].y + leta[k][j-1][i].y) * laj[k][j][i];
                PetscReal eta3 = 0.5 * (leta[k][j][i].z + leta[k][j-1][i].z) * laj[k][j][i];

                PetscReal zet1 = 0.5 * (lzet[k][j][i].x + lzet[k-1][j][i].x) * laj[k][j][i];
                PetscReal zet2 = 0.5 * (lzet[k][j][i].y + lzet[k-1][j][i].y) * laj[k][j][i];
                PetscReal zet3 = 0.5 * (lzet[k][j][i].z + lzet[k-1][j][i].z) * laj[k][j][i];

                // Calculate velocity gradient tensor components d_ij = du_i/dx_j
                PetscReal d11 = uc * csi1 + ue * eta1 + uz * zet1;
                PetscReal d12 = uc * csi2 + ue * eta2 + uz * zet2;
                PetscReal d13 = uc * csi3 + ue * eta3 + uz * zet3;

                PetscReal d21 = vc * csi1 + ve * eta1 + vz * zet1;
                PetscReal d22 = vc * csi2 + ve * eta2 + vz * zet2;
                PetscReal d23 = vc * csi3 + ve * eta3 + vz * zet3;
                
                PetscReal d31 = wc * csi1 + we * eta1 + wz * zet1;
                PetscReal d32 = wc * csi2 + we * eta2 + wz * zet2;
                PetscReal d33 = wc * csi3 + we * eta3 + wz * zet3;

                // Strain-Rate Tensor S_ij = 0.5 * (d_ij + d_ji)
                PetscReal s11 = d11;
                PetscReal s12 = 0.5 * (d12 + d21);
                PetscReal s13 = 0.5 * (d13 + d31);
                PetscReal s22 = d22;
                PetscReal s23 = 0.5 * (d23 + d32);
                PetscReal s33 = d33;

                // Vorticity Tensor Omega_ij = 0.5 * (d_ij - d_ji)
                PetscReal w12 = 0.5 * (d12 - d21);
                PetscReal w13 = 0.5 * (d13 - d31);
                PetscReal w23 = 0.5 * (d23 - d32);

                // Squared norms of the tensors
                PetscReal s_norm_sq = s11*s11 + s22*s22 + s33*s33 + 2.0*(s12*s12 + s13*s13 + s23*s23);
                PetscReal w_norm_sq = 2.0 * (w12*w12 + w13*w13 + w23*w23);
                
                gq[k][j][i] = 0.5 * (w_norm_sq - s_norm_sq);

                if (lnvert[k][j][i] > 0.1) {
                    gq[k][j][i] = 0.0;
                }
            }
        }
    }

    // --- 5. Restore arrays ---
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lUcat,  (void*)&lucat);   CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi,   (void*)&lcsi);    CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta,   (void*)&leta);    CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet,   (void*)&lzet);    CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da,  user->lAj,    (void*)&laj);     CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da,  user->lNvert, (void*)&lnvert);  CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da,  user->Qcrit, (void*)&gq);       CHKERRQ(ierr);

    /* The loop above skips the layout boundary, so extend it before anything reads
     * Qcrit with a stencil. Q is not a moment: it does not vanish at a wall, so only
     * the periodic case is defined and only that case is written. */
    ierr = ExtendToLayoutBoundary(user, user->Qcrit, 1); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "NormalizeRelativeField"
/**
 * @brief Implementation of \ref NormalizeRelativeField().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/postprocessing_kernels.h`.
 * @see NormalizeRelativeField()
 */
PetscErrorCode NormalizeRelativeField(UserCtx* user, const char* relative_field_name)
{
    PetscErrorCode ierr;
    Vec            P_vec = NULL;
    DMDALocalInfo  info;
    PetscInt       ip=1, jp=1, kp=1; // Default reference point
    PetscReal      p_ref = 0.0;
    PetscReal      p_ref_local = 0.0;
    PetscInt       found_local = 0, found_global = 0;
    PostProcessParams *pps = user->simCtx->pps;

    // Fetch the logical reference point from pps.
    ip = pps->reference[0];
    jp = pps->reference[1];
    kp = pps->reference[2];

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(GLOBAL, LOG_INFO, "-> KERNEL: Running NormalizeRelativeField on '%s'.\n", relative_field_name);

    // --- 1. Map string argument to the PETSc Vec ---
    if (strcasecmp(relative_field_name, "P") == 0) {
        P_vec = user->P;
    } else {
        SETERRQ(PETSC_COMM_SELF, 1, "NormalizeRelativeField only supports the primary 'P' field , not '%s' currently.", relative_field_name);
    }

    // --- 2. Read the logical reference point from whichever rank owns it ---
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    PetscCheck(ip >= 0 && ip < info.mx && jp >= 0 && jp < info.my &&
               kp >= 0 && kp < info.mz,
               PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
               "Reference point (%" PetscInt_FMT ", %" PetscInt_FMT ", %" PetscInt_FMT
               ") lies outside the %" PetscInt_FMT "x%" PetscInt_FMT "x%" PetscInt_FMT
               " pressure layout.", ip, jp, kp, info.mx, info.my, info.mz);
    if (ip >= info.xs && ip < info.xs + info.xm &&
        jp >= info.ys && jp < info.ys + info.ym &&
        kp >= info.zs && kp < info.zs + info.zm) {
        const PetscReal ***pressure = NULL;

        ierr = DMDAVecGetArrayRead(user->da, P_vec, &pressure); CHKERRQ(ierr);
        p_ref_local = pressure[kp][jp][ip];
        ierr = DMDAVecRestoreArrayRead(user->da, P_vec, &pressure); CHKERRQ(ierr);
        found_local = 1;
    }
    ierr = MPI_Allreduce(&p_ref_local, &p_ref, 1, MPIU_REAL, MPI_SUM,
                         PETSC_COMM_WORLD); CHKERRQ(ierr);
    ierr = MPI_Allreduce(&found_local, &found_global, 1, MPIU_INT, MPI_SUM,
                         PETSC_COMM_WORLD); CHKERRQ(ierr);
    PetscCheck(found_global == 1, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
               "Reference pressure point must have exactly one owner; found %" PetscInt_FMT ".",
               found_global);
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "%s reference point (%" PetscInt_FMT ", %" PetscInt_FMT ", %" PetscInt_FMT
              ") has value %g.\n", relative_field_name, ip, jp, kp, (double)p_ref);

    // --- 3. Perform the normalization (in-place shift) on the full distributed vector ---
    ierr = VecShift(P_vec, -p_ref); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "%s field normalized by subtracting %g.\n", relative_field_name, p_ref);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

// ===========================================================================
// Particle Post-Processing Kernels
// ===========================================================================
#undef __FUNCT__
#define __FUNCT__ "ComputeSpecificKE"
/**
 * @brief Internal helper implementation: `ComputeSpecificKE()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ComputeSpecificKE(UserCtx* user, const char* velocity_field, const char* ske_field)
{
    PetscErrorCode ierr;
    PetscInt n_local;
    const PetscScalar (*vel_arr)[3]; // Access velocity as array of 3-component vectors
    PetscScalar *ske_arr;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(GLOBAL, LOG_INFO, "-> KERNEL: Running ComputeSpecificKE ('%s' -> '%s').\n", velocity_field, ske_field);

    // Get local data arrays from the DMSwarm
    ierr = DMSwarmGetLocalSize(user->swarm, &n_local); CHKERRQ(ierr);
    if (n_local == 0) { PROFILE_FUNCTION_END; PetscFunctionReturn(0); }

    // Get read-only access to velocity and write access to the output field
    ierr = DMSwarmGetField(user->swarm, velocity_field, NULL, NULL, (void**)&vel_arr); CHKERRQ(ierr);
    ierr = DMSwarmGetField(user->post_swarm, ske_field, NULL, NULL, (void**)&ske_arr); CHKERRQ(ierr);

    // Main computation loop
    for (PetscInt p = 0; p < n_local; p++) {
        const PetscScalar u = vel_arr[p][0];
        const PetscScalar v = vel_arr[p][1];
        const PetscScalar w = vel_arr[p][2];
        const PetscScalar vel_sq = u*u + v*v + w*w;
        ske_arr[p] = 0.5 * vel_sq;
    }

    // Restore arrays
    ierr = DMSwarmRestoreField(user->swarm, velocity_field, NULL, NULL, (void**)&vel_arr); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(user->post_swarm, ske_field, NULL, NULL, (void**)&ske_arr); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeDisplacement"
/**
 * @brief Internal helper implementation: `ComputeDisplacement()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ComputeDisplacement(UserCtx *user, const char *disp_field)
{
    PetscErrorCode  ierr;
    PetscInt        n_local;
    const PetscReal (*pos_arr)[3];
    PetscScalar    *disp_out;
    SimCtx         *simCtx = user->simCtx;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(GLOBAL, LOG_INFO, "-> KERNEL: Running ComputeDisplacement (-> '%s').\n", disp_field);

    ierr = DMSwarmGetLocalSize(user->swarm, &n_local); CHKERRQ(ierr);
    if (n_local == 0) { PROFILE_FUNCTION_END; PetscFunctionReturn(0); }

    const PetscReal x0 = simCtx->psrc_x;
    const PetscReal y0 = simCtx->psrc_y;
    const PetscReal z0 = simCtx->psrc_z;

    ierr = DMSwarmGetField(user->swarm,      ParticleFieldName(PARTICLE_FIELD_ID_POSITION), NULL, NULL, (void**)&pos_arr); CHKERRQ(ierr);
    ierr = DMSwarmGetField(user->post_swarm, disp_field, NULL, NULL, (void**)&disp_out);       CHKERRQ(ierr);

    for (PetscInt p = 0; p < n_local; p++) {
        const PetscReal dx = pos_arr[p][0] - x0;
        const PetscReal dy = pos_arr[p][1] - y0;
        const PetscReal dz = pos_arr[p][2] - z0;
        disp_out[p] = PetscSqrtReal(dx*dx + dy*dy + dz*dz);
    }

    ierr = DMSwarmRestoreField(user->swarm,      ParticleFieldName(PARTICLE_FIELD_ID_POSITION), NULL, NULL, (void**)&pos_arr); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(user->post_swarm, disp_field, NULL, NULL, (void**)&disp_out);       CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
