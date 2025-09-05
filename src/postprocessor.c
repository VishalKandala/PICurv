/**
 * @file postprocessor.c
 * @brief Phase 2 implementation of the post-processing tool.
 *
 * This phase introduces a dedicated configuration system and performs a
 * single-step data load to verify the I/O and data structures.
 */

#include "postprocessor.h" // Use our new header

/**
 * @brief Parses the processing pipeline string from the config and executes the requested kernels in sequence.
 *
 * This function uses a general-purpose parser to handle a syntax of the form:
 * "Keyword1:in1>out1; Keyword2:in1,in2>out2; Keyword3:arg1;"
 *
 * It tokenizes the pipeline string and dispatches to the appropriate kernel function
 * from processing_kernels.c with the specified field name arguments.
 *
 * @param user The UserCtx containing the data to be transformed.
 * @param pps  The PostProcessParams struct containing the pipeline string.
 * @return PetscErrorCode
 */
PetscErrorCode RunProcessingPipeline(UserCtx* user, PostProcessParams* pps)
{
    PetscErrorCode ierr;
    char *pipeline_copy, *step_token, *step_saveptr;

    PetscFunctionBeginUser;
    LOG_ALLOW(GLOBAL, LOG_INFO, "--- Starting Data Transformation Pipeline ---\n");
    
    // Do nothing if the pipeline string is empty
    if (pps->process_pipeline[0] == '\0') {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Processing pipeline is empty. No transformations will be run.\n");
        PetscFunctionReturn(0);
    }
    
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Pipeline string: [%s]\n", pps->process_pipeline);
    
    // Make a writable copy for strtok_r, as it modifies the string
    ierr = PetscStrallocpy(pps->process_pipeline, &pipeline_copy); CHKERRQ(ierr);

    // --- Outer Loop: Tokenize by Semicolon (;) to get each processing step ---
    step_token = strtok_r(pipeline_copy, ";", &step_saveptr);
    while (step_token) {
        TrimWhitespace(step_token);
        if (strlen(step_token) == 0) {
            step_token = strtok_r(NULL, ";", &step_saveptr);
            continue;
        }

        char *keyword = strtok(step_token, ":");
        char *args_str = strtok(NULL, ""); // Get the rest of the string as arguments

        if (!keyword) { // Should not happen with TrimWhitespace, but is a safe check
            step_token = strtok_r(NULL, ";", &step_saveptr);
            continue;
        }
        
        TrimWhitespace(keyword);
        if (args_str) TrimWhitespace(args_str);

        LOG_ALLOW(GLOBAL, LOG_INFO, "Executing Transformation: '%s' on args: '%s'\n", keyword, args_str ? args_str : "None");

        // --- DISPATCHER: Route to the correct kernel based on the keyword ---
        if (strcasecmp(keyword, "CellToNodeAverage") == 0) {
            if (!args_str) SETERRQ(PETSC_COMM_SELF, 1, "CellToNodeAverage requires arguments in 'in_field>out_field' format.");
            char *in_field = strtok(args_str, ">");
            char *out_field = strtok(NULL, ">");
            if (!in_field || !out_field) SETERRQ(PETSC_COMM_SELF, 1, "CellToNodeAverage requires 'in>out' syntax (e.g., P>P_nodal).");
            TrimWhitespace(in_field); TrimWhitespace(out_field);
            ierr = ComputeNodalAverage(user, in_field, out_field); CHKERRQ(ierr);
        }
        else if (strcasecmp(keyword, "ComputeQCriterion") == 0) {
            ierr = ComputeQCriterion(user); CHKERRQ(ierr);
        }
        else if (strcasecmp(keyword, "NormalizeRelativeField") == 0) {
            if (!args_str) SETERRQ(PETSC_COMM_SELF, 1, "NormalizePressure requires the pressure field name (e.g., 'P') as an argument.");
            ierr = NormalizeRelativeField(user, args_str); CHKERRQ(ierr);
        }
        // *** Add new kernels here in the future using 'else if' ***
        // else if (strcasecmp(keyword, "ComputeVorticity") == 0) { ... }
        else {
            LOG_ALLOW(GLOBAL, LOG_WARNING, "Unknown transformation keyword '%s'. Skipping.\n", keyword);
        }

        step_token = strtok_r(NULL, ";", &step_saveptr);
    }

    ierr = PetscFree(pipeline_copy); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "--- Data Transformation Pipeline Complete ---\n");
    PetscFunctionReturn(0);
}

PetscErrorCode WriteCombinedVTKFile(UserCtx* user, PostProcessParams* pps, PetscInt ti)
{
    PetscErrorCode ierr;
    VTKMetaData    meta;
    char           filename[MAX_FILENAME_LENGTH];

    PetscFunctionBeginUser;

    if (pps->output_fields_instantaneous[0] == '\0') {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "No instantaneous fields requested for output at ti=%" PetscInt_FMT ". Skipping.\n", ti);
        PetscFunctionReturn(0);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "--- Starting VTK File Writing for ti = %" PetscInt_FMT " ---\n", ti);

    /* 1) Metadata init */
    ierr = PetscMemzero(&meta, sizeof(VTKMetaData)); CHKERRQ(ierr);
    meta.fileType = VTK_STRUCTURED;

    /* 2) Coordinates (subsampled interior) */
    ierr = PrepareOutputCoordinates(user, &meta.coords, &meta.mx, &meta.my, &meta.mz, &meta.npoints); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Using coords linearization order: fast=i mid=j slow=k  (sizes: %" PetscInt_FMT " x %" PetscInt_FMT " x %" PetscInt_FMT ")\n",
              meta.mx, meta.my, meta.mz);

    /* 3) Fields (rank 0) */
    if (user->simCtx->rank == 0) {
        char *fields_copy, *field_name;
        ierr = PetscStrallocpy(pps->output_fields_instantaneous, &fields_copy); CHKERRQ(ierr);

        field_name = strtok(fields_copy, ",");
        while (field_name) {
            TrimWhitespace(field_name);
            if (!*field_name) { field_name = strtok(NULL, ","); continue; }

            LOG_ALLOW(LOCAL, LOG_DEBUG, "Preparing field '%s' for output.\n", field_name);

            Vec       field_vec = NULL;
            PetscInt  num_components = 0;

            if (!strcasecmp(field_name, "P_nodal")) {
                field_vec = user->P_nodal; num_components = 1;
            } else if (!strcasecmp(field_name, "Ucat_nodal")) {
                field_vec = user->Ucat_nodal; num_components = 3;
            } else if (!strcasecmp(field_name, "Qcrit")) {
                field_vec = user->Qcrit;    num_components = 1;
            } else {
                LOG_ALLOW(LOCAL, LOG_WARNING, "Field '%s' not recognized. Skipping.\n", field_name);
                field_name = strtok(NULL, ",");
                continue;
            }

            if (meta.num_point_data_fields >= MAX_POINT_DATA_FIELDS) {
                LOG_ALLOW(LOCAL, LOG_WARNING, "MAX_POINT_DATA_FIELDS reached. Cannot add '%s'.\n", field_name);
                field_name = strtok(NULL, ",");
                continue;
            }

            // --- Add field to metadata ---
            VTKFieldInfo* current_field = &meta.point_data_fields[meta.num_point_data_fields];
            strncpy(current_field->name, field_name, MAX_VTK_FIELD_NAME_LENGTH-1);
            current_field->name[MAX_VTK_FIELD_NAME_LENGTH-1] = '\0';
            current_field->num_components = num_components;

            /* Build interior AoS from NATURAL gathered Vec */
            ierr = PrepareOutputFieldData(user, field_vec, num_components, &current_field->data); CHKERRQ(ierr);
            
            /*
            // *** DEBUG: Dump Ucat_nodal details and add scalar companions Ux, Uy, Uz for easier visualization ***
            
            // If this is Ucat_nodal, dump a few tuples and add scalar companions 
            if (!strcasecmp(field_name, "Ucat_nodal")) {
                const PetscInt npts = meta.npoints;
                const PetscScalar *a = (const PetscScalar*)current_field->data;

                LOG_ALLOW(GLOBAL, LOG_INFO, "DBG Ucat_nodal: ptr=%p  npoints=%" PetscInt_FMT "  num_components=%" PetscInt_FMT "\n",
                          (void*)a, npts, current_field->num_components);

                if (a && current_field->num_components == 3 && npts > 0) {
                    const PetscInt nshow = (npts < 5) ? npts : 5;
                    LOG_ALLOW(GLOBAL, LOG_INFO, "DBG Ucat_nodal: showing first %d of %" PetscInt_FMT " tuples (AoS x,y,z):\n",
                              (int)nshow, npts);
                    for (PetscInt t = 0; t < nshow; ++t) {
                        LOG_ALLOW(GLOBAL, LOG_INFO, "  Ucat_nodal[%3" PetscInt_FMT "] = (%g, %g, %g)\n",
                                  t, (double)a[3*t+0], (double)a[3*t+1], (double)a[3*t+2]);
                    }
                    if (npts > 10) {
                        PetscInt mid = npts / 2;
                        LOG_ALLOW(GLOBAL, LOG_INFO, "  Ucat_nodal[mid=%" PetscInt_FMT "] = (%g, %g, %g)\n",
                                  mid, (double)a[3*mid+0], (double)a[3*mid+1], (double)a[3*mid+2]);
                    }

                    // Add scalar companions from the AoS we just created 
                    PetscScalar *Ux=NULL,*Uy=NULL,*Uz=NULL;
                    ierr = PetscMalloc1(meta.npoints, &Ux); CHKERRQ(ierr);
                    ierr = PetscMalloc1(meta.npoints, &Uy); CHKERRQ(ierr);
                    ierr = PetscMalloc1(meta.npoints, &Uz); CHKERRQ(ierr);
                    for (PetscInt i = 0; i < meta.npoints; ++i) {
                        Ux[i] = a[3*i+0]; Uy[i] = a[3*i+1]; Uz[i] = a[3*i+2];
                    }

                    if (meta.num_point_data_fields + 3 <= MAX_POINT_DATA_FIELDS) {
                        VTKFieldInfo *fx = &meta.point_data_fields[++meta.num_point_data_fields];
                        strncpy(fx->name, "Ux_debug", MAX_VTK_FIELD_NAME_LENGTH-1);
                        fx->name[MAX_VTK_FIELD_NAME_LENGTH-1] = '\0';
                        fx->num_components = 1; fx->data = Ux;

                        VTKFieldInfo *fy = &meta.point_data_fields[++meta.num_point_data_fields];
                        strncpy(fy->name, "Uy_debug", MAX_VTK_FIELD_NAME_LENGTH-1);
                        fy->name[MAX_VTK_FIELD_NAME_LENGTH-1] = '\0';
                        fy->num_components = 1; fy->data = Uy;

                        VTKFieldInfo *fz = &meta.point_data_fields[++meta.num_point_data_fields];
                        strncpy(fz->name, "Uz_debug", MAX_VTK_FIELD_NAME_LENGTH-1);
                        fz->name[MAX_VTK_FIELD_NAME_LENGTH-1] = '\0';
                        fz->num_components = 1; fz->data = Uz;

                        LOG_ALLOW(GLOBAL, LOG_INFO, "DBG: Added scalar companions Ux_debug, Uy_debug, Uz_debug.\n");
                    } else {
                        LOG_ALLOW(GLOBAL, LOG_WARNING, "DBG: Not enough slots to add Ux/Uy/Uz debug fields.\n");
                        PetscFree(Ux); PetscFree(Uy); PetscFree(Uz);
                    }

                    // Mid-plane CSV + AoS vs NATURAL compare (component X) 
                    
                    // Gather NATURAL again (small cost, but isolated and clear) 
                    PetscInt   Ng = 0;
                    double    *nat_d = NULL;
                    DM         dmU = NULL;
                    DMDALocalInfo infU;
                    ierr = VecGetDM(field_vec, &dmU); CHKERRQ(ierr);
                    ierr = DMDAGetLocalInfo(dmU, &infU); CHKERRQ(ierr);

                    const PetscInt M=infU.mx, N=infU.my, P=infU.mz;
                    const PetscInt mx = meta.mx, my = meta.my, mz = meta.mz;
                    const PetscInt iInnerMid = mx/2; // interior index [0..mx-1] 
                    const PetscInt iGlob     = iInnerMid;

                    ierr = VecToArrayOnRank0(field_vec, &Ng, &nat_d); CHKERRQ(ierr);

                    if (nat_d) {
                        const PetscScalar *nar = (const PetscScalar*)nat_d;
                        const char *base = pps->output_prefix;
                        char fn[512], fnc[512];
                        snprintf(fn,  sizeof(fn),  "%s_%05" PetscInt_FMT "_iMid.csv",       base, ti);
                        snprintf(fnc, sizeof(fnc), "%s_%05" PetscInt_FMT "_iMid_compare.csv", base, ti);

                        FILE *fp  = fopen(fn,  "w");
                        FILE *fpc = fopen(fnc, "w");
                        if (fp)  fprintf(fp,  "jInner,kInner,Ux,Uy,Uz\n");
                        if (fpc) fprintf(fpc, "jInner,kInner,Ux_AoS,Ux_NAT,abs_diff\n");

                        double maxAbsDiff = 0.0, sumAbs = 0.0;
                        PetscInt count = 0;

                        for (PetscInt kInner = 0; kInner < mz; ++kInner) {
                            const PetscInt k = kInner;
                            for (PetscInt jInner = 0; jInner < my; ++jInner) {
                                const PetscInt j = jInner;

                                // AoS tuple index 
                                const PetscInt t = iInnerMid + mx * (jInner + my * kInner);
                                const PetscScalar ux = a[3*t+0], uy = a[3*t+1], uz = a[3*t+2];

                                if (fp) fprintf(fp, "%d,%d,%.15e,%.15e,%.15e\n",
                                                (int)jInner,(int)kInner,(double)ux,(double)uy,(double)uz);

                                // NATURAL base for (iGlob,j,k) 
                                const PetscInt baseNat = 3 * (((k)*N + j)*M + iGlob);
                                const PetscScalar uxN = nar[baseNat + 0];

                                const double diff = fabs((double)ux - (double)uxN);
                                if (diff > maxAbsDiff) maxAbsDiff = diff;
                                sumAbs += diff; ++count;

                                if (fpc) fprintf(fpc, "%d,%d,%.15e,%.15e,%.15e\n",
                                                    (int)jInner,(int)kInner,(double)ux,(double)uxN,diff);
                            } // for jInner
                        } // for kInner
                        if (fp)  fclose(fp);
                        if (fpc) fclose(fpc);

                        if (count > 0) {
                            const double meanAbs = sumAbs / (double)count;
                            LOG_ALLOW(GLOBAL, LOG_INFO,
                                        "PETSc-Vec vs AoS (i-mid, Ux): max|Δ|=%.6e, mean|Δ|=%.6e  -> CSV: %s\n",
                                        maxAbsDiff, meanAbs, fnc);
                            LOG_ALLOW(GLOBAL, LOG_INFO, "Wrote i-mid plane CSV: %s\n", fn);
                        } // if count>0

                        ierr = PetscFree(nat_d); CHKERRQ(ierr);
                        
                        
                    } // if nat_d

                } // if a && num_components==3 && npts>0
            } // if Ucat_nodal        
            // --- END DEBUG BLOCK (Ucat_nodal) ---
            */
            
            meta.num_point_data_fields++; /* count the main field we just added */
            field_name = strtok(NULL, ",");
        }

        ierr = PetscFree(fields_copy); CHKERRQ(ierr);

        // --- DEBUG: Add sanity fields i_idx, j_idx, k_idx and x_pos, y_pos, z_pos ---
        // These are the logical indices and physical coordinates of each point in the subsampled grid.
        // They can be used to verify the grid structure and orientation in visualization tools.
        // They are added as scalar fields with names "i_idx", "j_idx", "k_idx" and "x_pos", "y_pos", "z_pos".
        // Note: these are only added if there is room in the MAX_POINT_DATA_FIELDS limit.
        // They are allocated and owned here, and will be freed below.
        // They are in the same linearization order as meta.coords (AoS x,y,z by point).
        /*
        // Append sanity fields i/j/k indices and coordinates 

        // Build i/j/k and x/y/z (length = npoints); these match the same linearization as coords 
        const PetscInt n = meta.npoints;
        PetscScalar *i_idx=NULL,*j_idx=NULL,*k_idx=NULL,*x_pos=NULL,*y_pos=NULL,*z_pos=NULL;

        ierr = PetscMalloc1(n, &i_idx); CHKERRQ(ierr);
        ierr = PetscMalloc1(n, &j_idx); CHKERRQ(ierr);
        ierr = PetscMalloc1(n, &k_idx); CHKERRQ(ierr);
        ierr = PetscMalloc1(n, &x_pos); CHKERRQ(ierr);
        ierr = PetscMalloc1(n, &y_pos); CHKERRQ(ierr);
        ierr = PetscMalloc1(n, &z_pos); CHKERRQ(ierr);

        // coords is length 3*n, AoS: (x,y,z) by point 
        const PetscScalar *c = (const PetscScalar*)meta.coords;
        for (PetscInt k = 0; k < meta.mz; ++k) {
            for (PetscInt j = 0; j < meta.my; ++j) {
                for (PetscInt i = 0; i < meta.mx; ++i) {
                    const PetscInt t = i + meta.mx * (j + meta.my * k);
                    i_idx[t] = (PetscScalar)i;
                    j_idx[t] = (PetscScalar)j;
                    k_idx[t] = (PetscScalar)k;
                    x_pos[t] = c[3*t+0];
                    y_pos[t] = c[3*t+1];
                    z_pos[t] = c[3*t+2];
                }
            }
        }

        const char *nf[6] = {"i_idx","j_idx","k_idx","x_pos","y_pos","z_pos"};
        PetscScalar *arrs[6] = {i_idx,j_idx,k_idx,x_pos,y_pos,z_pos};
        for (int s=0; s<6; ++s) {
            if (meta.num_point_data_fields < MAX_POINT_DATA_FIELDS) {
                VTKFieldInfo *f = &meta.point_data_fields[meta.num_point_data_fields++];
                strncpy(f->name, nf[s], MAX_VTK_FIELD_NAME_LENGTH-1);
                f->name[MAX_VTK_FIELD_NAME_LENGTH-1] = '\0';
                f->num_components = 1;
                f->data = arrs[s];
            } else {
                LOG_ALLOW(GLOBAL, LOG_WARNING, "Sanity field '%s' dropped: MAX_POINT_DATA_FIELDS reached.\n", nf[s]);
                PetscFree(arrs[s]);
            }
        }
        LOG_ALLOW(GLOBAL, LOG_INFO, "DBG: Added sanity fields i_idx/j_idx/k_idx and x_pos/y_pos/z_pos.\n");
        */
        // --- END DEBUG BLOCK (sanity fields) ---

        /* Field summary */
        LOG_ALLOW(GLOBAL, LOG_INFO, "PointData fields to write: %d\n", (int)meta.num_point_data_fields);
        for (PetscInt ii=0; ii<meta.num_point_data_fields; ++ii) {
            LOG_ALLOW(GLOBAL, LOG_INFO, "  # %2" PetscInt_FMT "  Field Name = %s  Components = %d\n",
                      ii, meta.point_data_fields[ii].name, (int)meta.point_data_fields[ii].num_components);
        }
    } // if rank 0

    /* 4) Write the VTS */
    snprintf(filename, sizeof(filename), "%s_%05" PetscInt_FMT ".vts", pps->output_prefix, ti);
    ierr = CreateVTKFileFromMetadata(filename, &meta, PETSC_COMM_WORLD); CHKERRQ(ierr);

    /* 5) Cleanup (rank 0) */
    if (user->simCtx->rank == 0) {
        ierr = PetscFree(meta.coords); CHKERRQ(ierr);
        for (PetscInt i = 0; i < meta.num_point_data_fields; i++) {
            ierr = PetscFree(meta.point_data_fields[i].data); CHKERRQ(ierr);
        }
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "--- VTK File Writing for ti = %" PetscInt_FMT " Complete ---\n", ti);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
    PetscErrorCode    ierr;
    SimCtx            *simCtx = NULL;

    // === I. INITIALIZE PETSC & MPI ===========================================
    ierr = PetscInitialize(&argc, &argv, (char *)0, "Unified Post-Processing Tool"); CHKERRQ(ierr);

    // === II. CONFIGURE SIMULATION & POST-PROCESSING CONTEXTS =================
    ierr = CreateSimulationContext(argc, argv, &simCtx); CHKERRQ(ierr);
    // === IIA. SET EXECUTION MODE (SOLVER vs POST-PROCESSOR) =====
    simCtx->exec_mode = EXEC_MODE_POSTPROCESSOR;
    // === III. SETUP GRID & DATA STRUCTURES ===================================
    ierr = SetupGridAndSolvers(simCtx); CHKERRQ(ierr);
    
    // Get the finest-level user context, as this is where we'll load data
    UserCtx *user = simCtx->usermg.mgctx[simCtx->usermg.mglevels-1].user;
    PostProcessParams *pps = simCtx->pps;

    LOG_ALLOW(GLOBAL, LOG_INFO, "=============================================================\n");
    LOG_ALLOW(GLOBAL, LOG_INFO, "PHASE 4: Starting full processing and writing loop...\n");

    // --- FINAL MAIN LOOP ---
    for (PetscInt ti = pps->startTime; ti <= pps->endTime; ti += pps->timeStep) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "--- Processing Time Step %" PetscInt_FMT " ---\n", ti);

        // 1. Load Data (UpdateLocalGhosts is called inside the kernels)
        ierr = ReadSimulationFields(user, ti); CHKERRQ(ierr);
        
        // 2. Transform Data
        ierr = RunProcessingPipeline(user, pps); CHKERRQ(ierr);

        // 3. Write Output
        ierr = WriteCombinedVTKFile(user, pps, ti); CHKERRQ(ierr);
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "=============================================================\n");
    LOG_ALLOW(GLOBAL, LOG_INFO, "PHASE 4 COMPLETE: Post-processing finished successfully.\n");
    



    // === V. FINALIZE =========================================================
   // ierr = FinalizeSimulation(simCtx); CHKERRQ(ierr);
    ierr = PetscFinalize();
    return ierr;
}