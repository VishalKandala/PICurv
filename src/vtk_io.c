#include "vtk_io.h"

//================================================================================
//                 STATIC (PRIVATE) HELPER FUNCTION PROTOTYPES
//================================================================================

static PetscErrorCode WriteVTKAppendedBlock(FILE *fp, const void *data, PetscInt num_elements, size_t element_size);
static PetscErrorCode WriteVTSXMLHeader(FILE *fp, const VTKMetaData *meta, PetscInt *boffset);
static PetscErrorCode WriteVTPXMLHeader(FILE *fp, const VTKMetaData *meta, PetscInt *boffset);
static PetscErrorCode WriteVTKFileHeader(FILE *fp, const VTKMetaData *meta, PetscInt *boffset);
static PetscErrorCode WriteVTKFileFooter(FILE *fp, const VTKMetaData *meta);

//================================================================================
//                 IMPLEMENTATION OF PRIVATE HELPERS
//================================================================================

static PetscErrorCode WriteVTKAppendedBlock(FILE *fp, const void *data, PetscInt num_elements, size_t element_size) {
    uint32_t block_size = num_elements * (uint32_t)element_size;
    if (fwrite(&block_size, sizeof(uint32_t), 1, fp) != 1) return PETSC_ERR_FILE_WRITE;
    if (fwrite(data, element_size, num_elements, fp) != (size_t)num_elements) return PETSC_ERR_FILE_WRITE;
    return 0;
}

static PetscErrorCode WriteVTSXMLHeader(FILE *fp, const VTKMetaData *meta, PetscInt *boffset)
{
    const char *byte_order = "LittleEndian";
    const char *precision  = "Float64";

    const PetscInt header = (PetscInt)sizeof(uint32_t);

    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"%s\">\n", byte_order);
    fprintf(fp, "  <StructuredGrid WholeExtent=\"0 %" PetscInt_FMT " 0 %" PetscInt_FMT " 0 %" PetscInt_FMT "\">\n", meta->mx - 1, meta->my - 1, meta->mz - 1);
    fprintf(fp, "    <Piece Extent=\"0 %" PetscInt_FMT " 0 %" PetscInt_FMT " 0 %" PetscInt_FMT "\">\n", meta->mx - 1, meta->my - 1, meta->mz - 1);

    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%" PetscInt_FMT "\" />\n",
            precision, *boffset);
    *boffset += header + 3 * meta->npoints * (PetscInt)sizeof(PetscScalar);
    fprintf(fp, "      </Points>\n");

    if (meta->num_point_data_fields > 0) {
        fprintf(fp, "      <PointData>\n");
        for (PetscInt i = 0; i < meta->num_point_data_fields; i++) {
            const VTKFieldInfo* field = &meta->point_data_fields[i];
            fprintf(fp, "        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"%" PetscInt_FMT "\" format=\"appended\" offset=\"%" PetscInt_FMT "\" />\n",
                    precision, field->name, field->num_components, *boffset);
            *boffset += header + field->num_components * meta->npoints * (PetscInt)sizeof(PetscScalar);
        }
        fprintf(fp, "      </PointData>\n");
    }

    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </StructuredGrid>\n");
    fprintf(fp, "  <AppendedData encoding=\"raw\">\n_");
    
    return 0;
}

static PetscErrorCode WriteVTPXMLHeader(FILE *fp, const VTKMetaData *meta, PetscInt *boffset)
{
    const char *byte_order = "LittleEndian";
    const char *precision     = "Float64";
    const char *int_type_str  = (sizeof(PetscInt) == 8) ? "Int64" : "Int32";
    const PetscInt header = (PetscInt)sizeof(uint32_t);

    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"%s\">\n", byte_order);
    fprintf(fp, "  <PolyData>\n");
    fprintf(fp, "    <Piece NumberOfPoints=\"%" PetscInt_FMT "\" NumberOfVerts=\"%" PetscInt_FMT "\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n",
            meta->npoints, meta->npoints);

    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%" PetscInt_FMT "\" />\n",
            precision, *boffset);
    *boffset += header + 3 * meta->npoints * (PetscInt)sizeof(PetscScalar);
    fprintf(fp, "      </Points>\n");
    
    if (meta->num_point_data_fields > 0) {
        fprintf(fp, "      <PointData>\n");
        for (PetscInt i = 0; i < meta->num_point_data_fields; i++) {
            const VTKFieldInfo* field = &meta->point_data_fields[i];
            fprintf(fp, "        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"%" PetscInt_FMT "\" format=\"appended\" offset=\"%" PetscInt_FMT "\" />\n",
                    precision, field->name, field->num_components, *boffset);
            *boffset += header + field->num_components * meta->npoints * (PetscInt)sizeof(PetscScalar);
        }
        fprintf(fp, "      </PointData>\n");
    }

    fprintf(fp, "      <Verts>\n");
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"connectivity\" format=\"appended\" offset=\"%" PetscInt_FMT "\" />\n",
            int_type_str, *boffset);
    *boffset += (uint32_t)sizeof(uint32_t) + meta->npoints * (PetscInt)sizeof(PetscInt);
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"offsets\" format=\"appended\" offset=\"%" PetscInt_FMT "\" />\n",
            int_type_str, *boffset);
    *boffset += (uint32_t)sizeof(uint32_t) + meta->npoints * (PetscInt)sizeof(PetscInt);
    fprintf(fp, "      </Verts>\n");

    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </PolyData>\n");
    fprintf(fp, "  <AppendedData encoding=\"raw\">\n_");

    return 0;
}

static PetscErrorCode WriteVTKFileHeader(FILE *fp, const VTKMetaData *meta, PetscInt *boffset)
{
    if (meta->fileType == VTK_STRUCTURED) {
        return WriteVTSXMLHeader(fp, meta, boffset);
    } else if (meta->fileType == VTK_POLYDATA) {
        return WriteVTPXMLHeader(fp, meta, boffset);
    }
    return PETSC_ERR_ARG_WRONG;
}

static PetscErrorCode WriteVTKFileFooter(FILE *fp, const VTKMetaData *meta)
{
    fprintf(fp, "\n  </AppendedData>\n");
    fprintf(fp, "</VTKFile>\n");
    return 0;
}


//================================================================================
//                 PUBLIC VTK WRITER FUNCTION
//================================================================================

PetscErrorCode CreateVTKFileFromMetadata(const char *filename, const VTKMetaData *meta, MPI_Comm comm)
{
    PetscMPIInt rank;
    MPI_Comm_rank(comm, &rank);
    PetscErrorCode ierr = 0;

    if (!rank) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Rank 0 writing combined VTK file '%s'.\n", filename);
        FILE *fp = fopen(filename, "wb");
        if (!fp) {
            LOG_ALLOW(GLOBAL, LOG_ERROR, "fopen failed for %s.\n", filename);
            return PETSC_ERR_FILE_OPEN;
        }

        PetscInt boffset = 0;
        
        ierr = WriteVTKFileHeader(fp, meta, &boffset);
        if(ierr) { fclose(fp); return ierr; }

        if (meta->coords) {
            ierr = WriteVTKAppendedBlock(fp, meta->coords, 3 * meta->npoints, sizeof(PetscScalar));
            if(ierr) { fclose(fp); return ierr; }
        }
        /*
        // ======== DEBUG: Dump first few values of Ucat_nodal if present ========
        if (!rank) {
            for (PetscInt i = 0; i < meta->num_point_data_fields; i++) {
                const VTKFieldInfo* f = &meta->point_data_fields[i];
                if (strcasecmp(f->name, "Ucat_nodal") == 0) {
                    const PetscScalar *a = f->data;
                    const PetscInt npts = meta->npoints;
                    const PetscInt nshow = (npts < 5) ? npts : 5;

                    LOG_ALLOW(GLOBAL, LOG_INFO,
                            "DBG (writer) Ucat_nodal: first %d of %" PetscInt_FMT " tuples:\n",
                            (int)nshow, npts);
                    for (PetscInt t = 0; t < nshow; ++t) {
                        LOG_ALLOW(GLOBAL, LOG_INFO,
                                "  Ucat_nodal[%3" PetscInt_FMT "] = (%g, %g, %g)\n",
                                t, (double)a[3*t+0], (double)a[3*t+1], (double)a[3*t+2]);
                    }
                }
            }
        }
        */
        // ======== END DEBUG ========
               

        for (PetscInt i = 0; i < meta->num_point_data_fields; i++) {
            const VTKFieldInfo* field = &meta->point_data_fields[i];
            if (field->data) {
                ierr = WriteVTKAppendedBlock(fp, field->data, field->num_components * meta->npoints, sizeof(PetscScalar));
                if(ierr) { fclose(fp); return ierr; }
            }
        }
        if (meta->fileType == VTK_POLYDATA) {
            if (meta->connectivity) {
                ierr = WriteVTKAppendedBlock(fp, meta->connectivity, meta->npoints, sizeof(PetscInt));
                if(ierr) { fclose(fp); return ierr; }
            }
            if (meta->offsets) {
                ierr = WriteVTKAppendedBlock(fp, meta->offsets, meta->npoints, sizeof(PetscInt));
                if(ierr) { fclose(fp); return ierr; }
            }
        }
        ierr = WriteVTKFileFooter(fp, meta);
        if(ierr) { fclose(fp); return ierr; }
        
        fclose(fp);
        LOG_ALLOW(GLOBAL, LOG_INFO, "Rank 0 finished writing VTK file '%s'.\n", filename);
    }
    return 0;
}

//================================================================================
//                 PUBLIC COORDINATE PREPARATION FUNCTION
//================================================================================

/**
 * @brief Creates a C array of coordinates corresponding to a subsampled (legacy-style) grid.
 *
 * This function gathers the full, distributed grid coordinates onto rank 0. On rank 0,
 * it then allocates a new, smaller C array and copies only the coordinates for the
 * nodes within the range [0..IM-2, 0..JM-2, 0..KM-2]. This produces a contiguous
 * array of points for a grid of size (IM-1)x(JM-1)x(KM-1), matching the legacy output.
 * The output arrays are only allocated and valid on rank 0.
 *
 * @param[in]  user        The UserCtx containing the grid information (DM, IM/JM/KM).
 * @param[out] out_coords  On rank 0, a pointer to the newly allocated C array for coordinate data. NULL on other ranks.
 * @param[out] out_nx      The number of points in the x-dimension for the new grid (IM-1).
 * @param[out] out_ny      The number of points in the y-dimension for the new grid (JM-1).
 * @param[out] out_nz      The number of points in the z-dimension for the new grid (KM-1).
 * @param[out] out_npoints The total number of points in the new grid.
 * @return PetscErrorCode
 */
PetscErrorCode PrepareOutputCoordinates(UserCtx* user, PetscScalar** out_coords, PetscInt* out_nx, PetscInt* out_ny, PetscInt* out_nz, PetscInt* out_npoints)
{
    PetscErrorCode ierr;
    Vec            coords_global;
    PetscInt       IM, JM, KM;
    PetscInt       N_full_coords;
    PetscScalar    *full_coords_arr = NULL;

    PetscFunctionBeginUser;
    ierr = DMDAGetInfo(user->da, NULL, &IM, &JM, &KM, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);

    // Set the dimensions of the output grid
    *out_nx = IM - 1;
    *out_ny = JM - 1;
    *out_nz = KM - 1;
    *out_npoints = (*out_nx) * (*out_ny) * (*out_nz);

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Preparing subsampled coordinates for a %" PetscInt_FMT "x%" PetscInt_FMT "x%" PetscInt_FMT " grid.\n", *out_nx, *out_ny, *out_nz);

    // --- Step 1: Gather the full coordinate vector to rank 0 ---
    ierr = DMGetCoordinates(user->da, &coords_global); CHKERRQ(ierr);
    // Reuse of your existing, proven utility function's logic
    ierr = VecToArrayOnRank0(coords_global, &N_full_coords, &full_coords_arr); CHKERRQ(ierr);

    // --- Step 2: On rank 0, subsample the gathered array ---
    if (user->simCtx->rank == 0) {
        // We expect N_full_coords to be IM * JM * KM * 3
        if (N_full_coords != IM * JM * KM * 3) {
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_UNEXPECTED, "Gathered coordinate array has wrong size. Expected %" PetscInt_FMT ", got %" PetscInt_FMT, IM * JM * KM * 3, N_full_coords);
        }

        // Allocate the smaller output C array
        ierr = PetscMalloc1(3 * (*out_npoints), out_coords); CHKERRQ(ierr);

        // Loop over the smaller grid dimensions and copy data
        PetscInt p_out = 0; // Index for the small output array
        for (PetscInt k = 0; k < *out_nz; k++) {
            for (PetscInt j = 0; j < *out_ny; j++) {
                for (PetscInt i = 0; i < *out_nx; i++) {
                    // Calculate the index in the full, 1D source array
                    PetscInt p_in = 3 * (k * (JM * IM) + j * IM + i);
                    
                    (*out_coords)[p_out++] = full_coords_arr[p_in + 0]; // x
                    (*out_coords)[p_out++] = full_coords_arr[p_in + 1]; // y
                    (*out_coords)[p_out++] = full_coords_arr[p_in + 2]; // z
                }
            }
        }
        // Free the temporary full array that was allocated by VecToArrayOnRank0
        ierr = PetscFree(full_coords_arr); CHKERRQ(ierr);
    } else {
        // On other ranks, the output pointer is NULL
        *out_coords = NULL;
    }

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Subsampled coordinates prepared on rank 0 with total %" PetscInt_FMT " points.\n", *out_npoints);

    PetscFunctionReturn(0);
}

//================================================================================
//               PUBLIC FIELD PREPARATION FUNCTION
//================================================================================

PetscErrorCode PrepareOutputEulerianFieldData(UserCtx *user, Vec field_vec, PetscInt num_components, PetscScalar** out_data)
{
    PetscErrorCode ierr;
    PetscInt       IM, JM, KM, nx, ny, nz, npoints;
    PetscInt       N_full_field = 0;
    PetscScalar   *full_field_arr = NULL;

    PetscFunctionBeginUser;
    ierr = DMDAGetInfo(user->da, NULL, &IM, &JM, &KM, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
    nx = IM - 1; ny = JM - 1; nz = KM - 1;
    npoints = nx * ny * nz;

    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "Preparing subsampled field data (dof=%" PetscInt_FMT ") for a %" PetscInt_FMT "x%" PetscInt_FMT "x%" PetscInt_FMT " grid.\n",
              num_components, nx, ny, nz);

    // --- Step 1: Gather the full field vector to rank 0 ---
    ierr = VecToArrayOnRank0(field_vec, &N_full_field, &full_field_arr); CHKERRQ(ierr);

    // --- Step 2: On rank 0, subsample the gathered array ---
    if (user->simCtx->rank == 0) {
        // Sanity check the size of the gathered array
        if (N_full_field != IM * JM * KM * num_components) {
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FILE_UNEXPECTED,
                    "Gathered field array has wrong size. Expected %" PetscInt_FMT ", got %" PetscInt_FMT,
                    IM * JM * KM * num_components, N_full_field);
        }

        // Allocate the smaller output C array (AoS layout expected by your writer)
        ierr = PetscMalloc1(num_components * npoints, out_data); CHKERRQ(ierr);

        /*
        // ======== Layout diagnostics (rank 0 only, read-only) ========
        if (num_components > 1 && full_field_arr) {
            const PetscInt n_full = IM * JM * KM;

            // Choose a safe interior probe point and its +i neighbor
            const PetscInt ic = (IM > 4) ? IM/2 : 1;
            const PetscInt jc = (JM > 4) ? JM/2 : 1;
            const PetscInt kc = (KM > 4) ? KM/2 : 1;
            const PetscInt ic1 = (ic + 1 < IM-1) ? (ic + 1) : ic; // stay interior

            const PetscInt t_full_c  = (kc * JM + jc) * IM + ic;
            const PetscInt t_full_c1 = (kc * JM + jc) * IM + ic1;

            // If source were AoS:  full[num_components * t + comp]
            PetscScalar aos_c0   = full_field_arr[num_components * t_full_c  + 0];
            PetscScalar aos_c1   = full_field_arr[num_components * t_full_c1 + 0];
            PetscScalar aos_c0_y = (num_components > 1) ? full_field_arr[num_components * t_full_c + 1] : 0.0;
            PetscScalar aos_c0_z = (num_components > 2) ? full_field_arr[num_components * t_full_c + 2] : 0.0;

            // If source were SoA:  full[comp * n_full + t]
            PetscScalar soa_c0   = full_field_arr[0 * n_full + t_full_c];
            PetscScalar soa_c1   = full_field_arr[0 * n_full + t_full_c1];
            PetscScalar soa_c0_y = (num_components > 1) ? full_field_arr[1 * n_full + t_full_c] : 0.0;
            PetscScalar soa_c0_z = (num_components > 2) ? full_field_arr[2 * n_full + t_full_c] : 0.0;

            double dAOS = fabs((double)(aos_c1 - aos_c0)); // spatial delta along +i
            double dSOA = fabs((double)(soa_c1 - soa_c0));

            const char *verdict = (dSOA < dAOS) ? "LIKELY SoA (component-major)" : "LIKELY AoS (interleaved)";

            LOG_ALLOW(GLOBAL, LOG_INFO,
                "Layout probe @ (i=%" PetscInt_FMT ", j=%" PetscInt_FMT ", k=%" PetscInt_FMT ") vs +i neighbor:\n"
                "  AoS-candidate:  Ux(center)=%.6e  Ux(+i)=%.6e  |Δ|=%.6e ; Uy(center)=%.6e  Uz(center)=%.6e\n"
                "  SoA-candidate:  Ux(center)=%.6e  Ux(+i)=%.6e  |Δ|=%.6e ; Uy(center)=%.6e  Uz(center)=%.6e\n"
                "  Verdict: %s\n",
                ic, jc, kc,
                (double)aos_c0,  (double)aos_c1,  dAOS, (double)aos_c0_y, (double)aos_c0_z,
                (double)soa_c0,  (double)soa_c1,  dSOA, (double)soa_c0_y, (double)soa_c0_z,
                verdict
            );
        } else if (num_components == 1) {
            LOG_ALLOW(GLOBAL, LOG_INFO, "Layout probe: scalar field (dof=1) — no component interleave to detect.\n");
        }
        // ======== END diagnostics ========
        */

        // --- PACK: assumes source is already AoS ---
        {
            PetscInt p_out = 0;
            for (PetscInt k = 0; k < nz; k++) {
                for (PetscInt j = 0; j < ny; j++) {
                    for (PetscInt i = 0; i < nx; i++) {
                        const PetscInt p_in_start = num_components * (k * (JM * IM) + j * IM + i);
                        for (PetscInt c = 0; c < num_components; c++) {
                            (*out_data)[p_out++] = full_field_arr[p_in_start + c];
                        }
                    }
                }
            }
        }

        // Free temporary gathered array
        ierr = PetscFree(full_field_arr); CHKERRQ(ierr);
    } else {
        // Other ranks: no allocation
        *out_data = NULL;
    }

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Subsampled field data prepared on rank 0 with total %" PetscInt_FMT " points.\n", npoints);

    PetscFunctionReturn(0);
}

/**
 * @brief Gathers, merges, subsamples, and prepares particle data for VTK output.
 *
 * This function orchestrates the preparation of particle data for writing. It is a
 * COLLECTIVE operation that all MPI ranks must enter, though the primary work of
 * data aggregation and preparation is performed on rank 0.
 *
 * The function follows a "merge-and-prepare" strategy:
 * 1.  It treats `user->swarm` as the read-only source of primary data (e.g., positions, velocity)
 *     that was loaded from disk.
 * 2.  It treats `user->pp_swarm` as the source of derived data (e.g., specific_ke) that
 *     was computed by the post-processing pipeline.
 * 3.  **On Rank 0 only**:
 *     a. It gathers the full particle coordinates from the source swarm. This determines
 *        the total number of particles before any subsampling.
 *     b. It parses the `pps->particle_fields` list. For each requested field, it
 *        determines whether to source it from `user->swarm` or `user->pp_swarm`.
 *     c. It gathers the full data for each requested field from its appropriate source.
 *     d. It performs strided subsampling on the coordinates and all gathered fields based on
 *        `pps->particle_output_freq`.
 *     e. It populates the VTKMetaData struct with the final, smaller, subsampled data arrays,
 *        making it ready for the file writer.
 *
 * @param[in]  user      The UserCtx, containing both user->swarm and user->pp_swarm.
 * @param[in]  pps       The PostProcessParams struct for configuration (field list, frequency).
 * @param[out] meta      A pointer to the VTKMetaData struct to be populated (on rank 0).
 * @param[out] p_n_total A pointer to store the total number of particles before subsampling (on rank 0).
 * @return PetscErrorCode
 */
PetscErrorCode PrepareOutputParticleData(UserCtx* user, PostProcessParams* pps, VTKMetaData* meta, PetscInt* p_n_total)
{
    PetscErrorCode ierr;
    PetscInt n_total_particles, n_components;

    PetscFunctionBeginUser;

    // Initialize output parameters on all ranks
    *p_n_total = 0;
    ierr = PetscMemzero(meta, sizeof(VTKMetaData)); CHKERRQ(ierr);

    // --- The entire preparation process is a series of collective gathers followed by rank-0 processing ---
    
    // --- Step 1: Gather Full Coordinates from SOURCE Swarm (Collective) ---
    // This establishes the ground truth for particle positions and total count.
    PetscScalar *full_coords_arr = NULL;
    ierr = SwarmFieldToArrayOnRank0(user->swarm, "position", &n_total_particles, &n_components, &full_coords_arr); CHKERRQ(ierr);
    
    // --- Step 2: Prepare and Subsample Data (Rank 0 Only) ---
    if (user->simCtx->rank == 0) {
        *p_n_total = n_total_particles; // Report original count back to the caller.

        if (n_total_particles == 0) {
            ierr = PetscFree(full_coords_arr); CHKERRQ(ierr); // Cleanup the gathered coords array
            PetscFunctionReturn(0); // Nothing to prepare.
        }
        if (n_components != 3) {
            SETERRQ(PETSC_COMM_SELF, 1, "Coordinate field position must have 3 components, but has %" PetscInt_FMT, n_components);
        }

        // --- Subsampling Calculation ---
        PetscInt stride = pps->particle_output_freq > 0 ? pps->particle_output_freq : 1;
        meta->npoints = (n_total_particles > 0) ? (n_total_particles - 1) / stride + 1 : 0;
        
        LOG_ALLOW(LOCAL, LOG_DEBUG, "Subsampling %" PetscInt_FMT " total particles with stride %" PetscInt_FMT " -> %" PetscInt_FMT " output particles.\n",
                  n_total_particles, stride, meta->npoints);

        // --- Prepare Final Coordinates Array ---
        ierr = PetscMalloc1(3 * meta->npoints, &meta->coords); CHKERRQ(ierr);
        for (PetscInt i = 0; i < meta->npoints; i++) {
            PetscInt source_idx = i * stride;
            for (int d = 0; d < 3; d++) meta->coords[3 * i + d] = full_coords_arr[3 * source_idx + d];
        }
        ierr = PetscFree(full_coords_arr); // Free the temporary full array immediately after use.


        // --- Prepare Final Field Arrays (Gather -> Subsample -> Store) ---
        char *fields_copy, *field_name;
        ierr = PetscStrallocpy(pps->particle_fields, &fields_copy); CHKERRQ(ierr);
        field_name = strtok(fields_copy, ",");
        while (field_name && meta->num_point_data_fields < MAX_POINT_DATA_FIELDS) {
            TrimWhitespace(field_name);
            if (!*field_name || strcasecmp(field_name, "position") == 0) {
                field_name = strtok(NULL, ","); continue;
            }

            // A. Determine which swarm is the source for this field.
            DM swarm_to_use = user->swarm; // Default to the main solver swarm.
            const char* internal_name = field_name; // Default internal name to user-facing name.

            // Use the "try-get" pattern to see if the field exists on the post_swarm.
            PetscErrorCode check_ierr;
            ierr = PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL); CHKERRQ(ierr);
            check_ierr = DMSwarmGetField(user->post_swarm, field_name, NULL, NULL, NULL);
            ierr = PetscPopErrorHandler(); CHKERRQ(ierr);
            if (!check_ierr) { // Success! The field exists on post_swarm.
                swarm_to_use = user->post_swarm;
                ierr = DMSwarmRestoreField(user->post_swarm, field_name, NULL, NULL, NULL); CHKERRQ(ierr);
                LOG_ALLOW(LOCAL, LOG_DEBUG, "Field '%s' will be sourced from the post-processing swarm.\n", field_name);
            } else {
                LOG_ALLOW(LOCAL, LOG_DEBUG, "Field '%s' will be sourced from the main solver swarm.\n", field_name);
                // Map friendly names to internal names for standard solver fields
                if (strcasecmp(field_name, "pid") == 0) internal_name = "DMSwarm_pid";
                else if (strcasecmp(field_name, "CellID") == 0) internal_name = "DMSwarm_CellID";
                else if (strcasecmp(field_name, "Migration Status") == 0) internal_name = "DMSwarm_location_status";
                // Add other mappings if necessary...
            }
            
            // B. Gather the full data for this field from its determined source.
            PetscScalar *full_field_arr = NULL;
            PetscInt field_total_particles, field_num_components;
            ierr = SwarmFieldToArrayOnRank0(swarm_to_use, internal_name, &field_total_particles, &field_num_components, &full_field_arr); CHKERRQ(ierr);

            if (field_total_particles != n_total_particles) {
                LOG_ALLOW(LOCAL, LOG_WARNING, "Field '%s' has %" PetscInt_FMT " particles, but expected %" PetscInt_FMT ". Skipping.\n", field_name, field_total_particles, n_total_particles);
                ierr = PetscFree(full_field_arr); CHKERRQ(ierr);
                field_name = strtok(NULL, ","); continue;
            }
            
            // C. Allocate final array and copy subsampled data into the meta struct.
            VTKFieldInfo* current_field = &meta->point_data_fields[meta->num_point_data_fields];
            strncpy(current_field->name, field_name, MAX_VTK_FIELD_NAME_LENGTH - 1);
            current_field->name[MAX_VTK_FIELD_NAME_LENGTH - 1] = '\0';
            current_field->num_components = field_num_components;
            ierr = PetscMalloc1(current_field->num_components * meta->npoints, &current_field->data); CHKERRQ(ierr);

            for (PetscInt i = 0; i < meta->npoints; i++) {
                PetscInt source_idx = i * stride;
                for (PetscInt c = 0; c < current_field->num_components; c++) {
                    current_field->data[current_field->num_components * i + c] = full_field_arr[current_field->num_components * source_idx + c];
                }
            }
            ierr = PetscFree(full_field_arr); // Free the temporary full array.
            
            meta->num_point_data_fields++;
            field_name = strtok(NULL, ",");
        }
        ierr = PetscFree(fields_copy); CHKERRQ(ierr);

        // --- Finalize VTK MetaData for PolyData ---
        if (meta->npoints > 0) {
            meta->fileType = VTK_POLYDATA;
            ierr = PetscMalloc1(meta->npoints, &meta->connectivity); CHKERRQ(ierr);
            ierr = PetscMalloc1(meta->npoints, &meta->offsets); CHKERRQ(ierr);
            for (PetscInt i = 0; i < meta->npoints; i++) {
                meta->connectivity[i] = i;
                meta->offsets[i] = i + 1;
            }
        }
    } // End of rank 0 block

    // A barrier is a good practice to ensure rank 0 finishes preparing before other ranks
    // potentially proceed to the next timestep where new collective calls might occur.
    ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}