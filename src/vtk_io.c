#include "vtk_io.h"
#include "particle_field_catalog.h"

//================================================================================
//                 STATIC (PRIVATE) HELPER FUNCTION PROTOTYPES
//================================================================================

static PetscErrorCode WriteVTSXMLHeader(FILE *fp, const VTKMetaData *meta, PetscInt *boffset);
static PetscErrorCode WriteVTPXMLHeader(FILE *fp, const VTKMetaData *meta, PetscInt *boffset);
static PetscErrorCode WriteVTKFileHeader(FILE *fp, const VTKMetaData *meta, PetscInt *boffset);
static PetscErrorCode CopyOwnedOutputBox(Vec vec, const PetscInt extent[3], PetscInt start[3], PetscInt size[3],
                                         PetscScalar **out);

//================================================================================
//                 IMPLEMENTATION OF PRIVATE HELPERS
//================================================================================

/**
 * @brief Write the XML header for a VTK structured-grid (`.vts`) file.
 */
static PetscErrorCode WriteVTSXMLHeader(FILE *fp, const VTKMetaData *meta, PetscInt *boffset)
{
    PROFILE_HELPER_BEGIN;
    const char *byte_order = "LittleEndian";
    const char *precision  = "Float64";

    const PetscInt header = (PetscInt)sizeof(uint64_t);

    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"%s\" header_type=\"UInt64\">\n", byte_order);
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
    
    PROFILE_HELPER_END;
    return 0;
}

/**
 * @brief Write the XML header for a VTK polydata (`.vtp`) file.
 */
static PetscErrorCode WriteVTPXMLHeader(FILE *fp, const VTKMetaData *meta, PetscInt *boffset)
{
    PROFILE_HELPER_BEGIN;
    const char *byte_order = "LittleEndian";
    const char *precision     = "Float64";
    const char *int_type_str  = (sizeof(PetscInt) == 8) ? "Int64" : "Int32";
    const PetscInt header = (PetscInt)sizeof(uint64_t);

    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"%s\" header_type=\"UInt64\">\n", byte_order);
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
    *boffset += header + meta->npoints * (PetscInt)sizeof(PetscInt);
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"offsets\" format=\"appended\" offset=\"%" PetscInt_FMT "\" />\n",
            int_type_str, *boffset);
    *boffset += header + meta->npoints * (PetscInt)sizeof(PetscInt);
    fprintf(fp, "      </Verts>\n");

    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </PolyData>\n");
    fprintf(fp, "  <AppendedData encoding=\"raw\">\n_");

    PROFILE_HELPER_END;
    return 0;
}

/**
 * @brief Open a VTK XML document and write its file-level header.
 */
static PetscErrorCode WriteVTKFileHeader(FILE *fp, const VTKMetaData *meta, PetscInt *boffset)
{
    PROFILE_FUNCTION_BEGIN;
    if (meta->fileType == VTK_STRUCTURED) {
        PROFILE_FUNCTION_END;
        return WriteVTSXMLHeader(fp, meta, boffset);
    } else if (meta->fileType == VTK_POLYDATA) {
        PROFILE_FUNCTION_END;
        return WriteVTPXMLHeader(fp, meta, boffset);
    }
    PROFILE_FUNCTION_END;
    return PETSC_ERR_ARG_WRONG;
}

/**
 * @brief Copy this rank's part of the output grid out of a DMDA global vector.
 * @details The output grid is the first `extent` nodes on each axis (the last DA node is a
 *          dummy). This rank's part is its owned box clipped to that grid, copied point by
 *          point (i fastest) with the vector's components interleaved.
 */
static PetscErrorCode CopyOwnedOutputBox(Vec vec, const PetscInt extent[3], PetscInt start[3], PetscInt size[3],
                                         PetscScalar **out)
{
    DM                 dm = NULL;
    DMDALocalInfo      info;
    const PetscScalar *values = NULL;
    PetscInt           p = 0;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCall(VecGetDM(vec, &dm));
    PetscCheck(dm != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "VTK output needs a vector attached to a DMDA.");
    PetscCall(DMDAGetLocalInfo(dm, &info));
    const PetscInt owned_start[3] = {info.xs, info.ys, info.zs}, owned[3] = {info.xm, info.ym, info.zm};
    for (int d = 0; d < 3; ++d) {
        start[d] = owned_start[d];
        size[d]  = PetscMax(0, PetscMin(owned_start[d] + owned[d], extent[d]) - owned_start[d]);
    }
    PetscCall(PetscMalloc1(info.dof * size[0] * size[1] * size[2], out));
    PetscCall(VecGetArrayRead(vec, &values));
    for (PetscInt k = 0; k < size[2]; ++k)
        for (PetscInt j = 0; j < size[1]; ++j)
            for (PetscInt i = 0; i < size[0]; ++i)
                for (PetscInt c = 0; c < info.dof; ++c)
                    (*out)[p++] = values[((k * owned[1] + j) * owned[0] + i) * info.dof + c];
    PetscCall(VecRestoreArrayRead(vec, &values));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}



/**
 * @brief Implementation of \ref CreateVTKFileFromMetadata().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/vtk_io.h`.
 * @see CreateVTKFileFromMetadata()
 */

PetscErrorCode CreateVTKFileFromMetadata(const char *filename, const VTKMetaData *meta, MPI_Comm comm)
{
    static const char  footer[] = "\n  </AppendedData>\n</VTKFile>\n";
    const PetscBool    structured = (PetscBool)(meta->fileType == VTK_STRUCTURED);
    char               tmp_filename[PETSC_MAX_PATH_LEN];
    PetscMPIInt        rank;
    long long          header_bytes = 0; /* -1: rank 0 could not write the header */
    MPI_File           fh;
    MPI_Offset         position;
    PetscInt           blocks = 0, local_points = 1, boffset = 0;
    struct VTKBlock { PetscInt components; MPI_Datatype type; size_t size; const void *data; } block[MAX_POINT_DATA_FIELDS + 3];

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCheck(structured || meta->fileType == VTK_POLYDATA, comm, PETSC_ERR_ARG_WRONG, "Unknown VTK file type.");
    PetscCallMPI(MPI_Comm_rank(comm, &rank));
    PetscCall(PetscSNPrintf(tmp_filename, sizeof(tmp_filename), "%s.tmp", filename));

    /* The appended blocks, in the order the header declares them. */
    block[blocks++] = (struct VTKBlock){3, MPIU_SCALAR, sizeof(PetscScalar), meta->coords};
    for (PetscInt f = 0; f < meta->num_point_data_fields; ++f)
        block[blocks++] = (struct VTKBlock){meta->point_data_fields[f].num_components, MPIU_SCALAR, sizeof(PetscScalar),
                                             meta->point_data_fields[f].data};
    if (!structured) {
        block[blocks++] = (struct VTKBlock){1, MPIU_INT, sizeof(PetscInt), meta->connectivity};
        block[blocks++] = (struct VTKBlock){1, MPIU_INT, sizeof(PetscInt), meta->offsets};
    }
    for (int d = 0; d < (structured ? 3 : 1); ++d) local_points *= meta->local_size[d];

    if (rank == 0) {
        FILE *fp = fopen(tmp_filename, "wb");
        if (!fp || WriteVTKFileHeader(fp, meta, &boffset)) header_bytes = -1;
        else header_bytes = (long long)ftell(fp);
        if (fp && fclose(fp)) header_bytes = -1;
    }
    PetscCallMPI(MPI_Bcast(&header_bytes, 1, MPI_LONG_LONG, 0, comm));
    PetscCheck(header_bytes >= 0, comm, PETSC_ERR_FILE_OPEN, "Could not write the VTK header of '%s'.", tmp_filename);

    /* Rank 0 writes each block's byte count; every rank writes its own points into the block. */
    PetscCallMPI(MPI_File_open(comm, tmp_filename, MPI_MODE_WRONLY, MPI_INFO_NULL, &fh));
    position = (MPI_Offset)header_bytes;
    for (PetscInt b = 0; b < blocks; ++b) {
        const uint64_t bytes = (uint64_t)block[b].components * (uint64_t)meta->npoints * block[b].size;
        MPI_Datatype   point, view;
        PetscMPIInt    count;

        PetscCallMPI(MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL));
        if (rank == 0) PetscCallMPI(MPI_File_write_at(fh, position, &bytes, (int)sizeof(bytes), MPI_BYTE, MPI_STATUS_IGNORE));
        position += (MPI_Offset)sizeof(bytes);

        PetscCallMPI(MPI_Type_contiguous((PetscMPIInt)block[b].components, block[b].type, &point));
        PetscCallMPI(MPI_Type_commit(&point));
        view = point;
        if (structured && local_points > 0) {
            const PetscMPIInt sizes[3]    = {(PetscMPIInt)meta->mz, (PetscMPIInt)meta->my, (PetscMPIInt)meta->mx};
            const PetscMPIInt subsizes[3] = {(PetscMPIInt)meta->local_size[2], (PetscMPIInt)meta->local_size[1], (PetscMPIInt)meta->local_size[0]};
            const PetscMPIInt starts[3]   = {(PetscMPIInt)meta->local_start[2], (PetscMPIInt)meta->local_start[1], (PetscMPIInt)meta->local_start[0]};
            PetscCallMPI(MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C, point, &view));
            PetscCallMPI(MPI_Type_commit(&view));
        }
        PetscCallMPI(MPI_File_set_view(fh, position, point, view, "native", MPI_INFO_NULL));
        PetscCall(PetscMPIIntCast(local_points, &count));
        if (structured) PetscCallMPI(MPI_File_write_all(fh, block[b].data, count, point, MPI_STATUS_IGNORE));
        else PetscCallMPI(MPI_File_write_at_all(fh, (MPI_Offset)meta->local_start[0], block[b].data, count, point, MPI_STATUS_IGNORE));
        if (view != point) PetscCallMPI(MPI_Type_free(&view));
        PetscCallMPI(MPI_Type_free(&point));
        position += (MPI_Offset)bytes;
    }
    PetscCallMPI(MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL));
    if (rank == 0) PetscCallMPI(MPI_File_write_at(fh, position, footer, (int)strlen(footer), MPI_BYTE, MPI_STATUS_IGNORE));
    PetscCallMPI(MPI_File_close(&fh));
    if (rank == 0) PetscCheck(rename(tmp_filename, filename) == 0, PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
                              "Could not move '%s' to '%s'.", tmp_filename, filename);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Wrote VTK file '%s'.\n", filename);
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

//================================================================================
//                 PUBLIC COORDINATE PREPARATION FUNCTION
//================================================================================

/**
 * @brief Internal helper implementation: `PrepareOutputCoordinates()`.
 * @details Local to this translation unit.
 */
PetscErrorCode PrepareOutputCoordinates(UserCtx* user, PetscScalar** out_coords, PetscInt* out_nx, PetscInt* out_ny, PetscInt* out_nz,
                                        PetscInt* out_npoints, PetscInt local_start[3], PetscInt local_size[3])
{
    Vec      coords_global;
    PetscInt IM, JM, KM;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    IM = user->info.mx; JM = user->info.my; KM = user->info.mz;
    *out_nx = IM - 1;
    *out_ny = JM - 1;
    *out_nz = KM - 1;
    *out_npoints = (*out_nx) * (*out_ny) * (*out_nz);
    const PetscInt extent[3] = {*out_nx, *out_ny, *out_nz};
    coords_global = user->Coor;
    PetscCall(CopyOwnedOutputBox(coords_global, extent, local_start, local_size, out_coords));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


/**
 * @brief Implementation of \ref BeginStructuredVTKOutput().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/vtk_io.h`.
 * @see BeginStructuredVTKOutput()
 */
PetscErrorCode BeginStructuredVTKOutput(UserCtx *user, VTKMetaData *meta)
{
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCheck(user != NULL && meta != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Context and metadata are required.");
    PetscCall(PetscMemzero(meta, sizeof(*meta)));
    meta->fileType = VTK_STRUCTURED;
    PetscCall(PrepareOutputCoordinates(user, &meta->coords, &meta->mx, &meta->my, &meta->mz,
                                       &meta->npoints, meta->local_start, meta->local_size));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref AppendStructuredVTKField().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/vtk_io.h`.
 * @see AppendStructuredVTKField()
 */
PetscErrorCode AppendStructuredVTKField(UserCtx *user, VTKMetaData *meta, const char *name,
                                        Vec field_vec, PetscInt components)
{
    VTKFieldInfo *entry = NULL;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCheck(user != NULL && meta != NULL && name != NULL && field_vec != NULL,
               PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Context, metadata, name, and vector are required.");
    if (meta->num_point_data_fields >= MAX_POINT_DATA_FIELDS) {
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "Reached the %d point-data field limit; '%s' is not written.\n",
                  MAX_POINT_DATA_FIELDS, name);
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }
    entry = &meta->point_data_fields[meta->num_point_data_fields];
    PetscCall(PetscStrncpy(entry->name, name, MAX_VTK_FIELD_NAME_LENGTH));
    entry->num_components = components;
    PetscCall(PrepareOutputEulerianFieldData(user, field_vec, components, &entry->data));
    meta->num_point_data_fields++;
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref FinishStructuredVTKOutput().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/vtk_io.h`.
 * @see FinishStructuredVTKOutput()
 */
PetscErrorCode FinishStructuredVTKOutput(VTKMetaData *meta, const char *filename)
{
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    PetscCheck(meta != NULL && filename != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Metadata and filename are required.");
    PetscCall(CreateVTKFileFromMetadata(filename, meta, PETSC_COMM_WORLD));
    for (PetscInt index = 0; index < meta->num_point_data_fields; ++index) {
        if (meta->point_data_fields[index].data) {
            PetscCall(PetscFree(meta->point_data_fields[index].data));
        }
    }
    if (meta->coords) PetscCall(PetscFree(meta->coords));
    meta->num_point_data_fields = 0;
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref PrepareOutputEulerianFieldData().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/vtk_io.h`.
 * @see PrepareOutputEulerianFieldData()
 */

PetscErrorCode PrepareOutputEulerianFieldData(UserCtx *user, Vec field_vec, PetscInt num_components, PetscScalar** out_data)
{
    PetscInt IM, JM, KM, dof, start[3], size[3];

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    IM = user->info.mx; JM = user->info.my; KM = user->info.mz;
    PetscCall(VecGetBlockSize(field_vec, &dof));
    PetscCheck(dof == num_components, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ,
               "Field has %" PetscInt_FMT " components per node, not the %" PetscInt_FMT " requested.", dof, num_components);
    const PetscInt extent[3] = {IM - 1, JM - 1, KM - 1};
    PetscCall(CopyOwnedOutputBox(field_vec, extent, start, size, out_data));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


/**
 * @brief Internal helper implementation: `PrepareOutputParticleData()`.
 * @details Local to this translation unit.
 */
PetscErrorCode PrepareOutputParticleData(UserCtx* user, PostProcessParams* pps, VTKMetaData* meta, PetscInt* p_n_total)
{
    const PetscInt stride = pps->particle_output_freq > 0 ? pps->particle_output_freq : 1;
    PetscInt       nlocal = 0, first = 0, before = 0, kept = 0, bs = 0;
    PetscDataType  type;
    PetscReal     *position = NULL;
    char          *fields_copy = NULL, *field_name = NULL;

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;

    /* Particles are numbered in rank order; output keeps every stride-th one, so this rank
       keeps its particles from local index `first` on, written from output point `before`. */
    PetscCall(DMSwarmGetLocalSize(user->swarm, &nlocal));
    PetscCallMPI(MPI_Exscan(&nlocal, &before, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD));
    if (user->simCtx->rank == 0) before = 0;
    PetscCall(DMSwarmGetSize(user->swarm, p_n_total));
    if (*p_n_total == 0) { PROFILE_HELPER_END; PetscFunctionReturn(0); }
    first = (stride - before % stride) % stride;
    kept  = first < nlocal ? (nlocal - 1 - first) / stride + 1 : 0;

    meta->fileType       = VTK_POLYDATA;
    meta->npoints        = (*p_n_total - 1) / stride + 1;
    meta->local_start[0] = (before + first) / stride;
    meta->local_size[0]  = kept;
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Subsampling %" PetscInt_FMT " particles with stride %" PetscInt_FMT " -> %" PetscInt_FMT " output particles.\n",
              *p_n_total, stride, meta->npoints);

    PetscCall(DMSwarmGetField(user->swarm, ParticleFieldName(PARTICLE_FIELD_ID_POSITION), &bs, &type, (void **)&position));
    PetscCheck(type == PETSC_REAL && bs == 3, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
               "Particle position field must be 3 PETSC_REAL components, not %" PetscInt_FMT " of %s.", bs, PetscDataTypes[type]);
    PetscCall(PetscMalloc1(3 * kept, &meta->coords));
    for (PetscInt i = 0; i < kept; ++i)
        for (int d = 0; d < 3; ++d) meta->coords[3 * i + d] = position[3 * (first + i * stride) + d];
    PetscCall(DMSwarmRestoreField(user->swarm, ParticleFieldName(PARTICLE_FIELD_ID_POSITION), &bs, &type, (void **)&position));

    PetscCall(PetscStrallocpy(pps->particle_fields, &fields_copy));
    for (field_name = strtok(fields_copy, ","); field_name && meta->num_point_data_fields < MAX_POINT_DATA_FIELDS;
         field_name = strtok(NULL, ",")) {
        DM             swarm_to_use = user->swarm;
        const char    *internal_name = field_name;
        PetscErrorCode check_ierr;
        PetscInt       field_local = 0, field_total = 0;
        void          *source = NULL;

        TrimWhitespace(field_name);
        if (!*field_name || strcasecmp(field_name, "position") == 0) continue;

        /* A field registered on the post-processing swarm comes from there, otherwise from the solver swarm. */
        PetscCall(PetscPushErrorHandler(PetscIgnoreErrorHandler, NULL));
        check_ierr = DMSwarmGetField(user->post_swarm, field_name, NULL, NULL, NULL);
        PetscCall(PetscPopErrorHandler());
        if (!check_ierr) {
            swarm_to_use = user->post_swarm;
            PetscCall(DMSwarmRestoreField(user->post_swarm, field_name, NULL, NULL, NULL));
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "Field '%s' will be sourced from the post-processing swarm.\n", field_name);
        } else {
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "Field '%s' will be sourced from the main solver swarm.\n", field_name);
            if (strcasecmp(field_name, "pid") == 0) internal_name = "DMSwarm_pid";
            else if (strcasecmp(field_name, "CellID") == 0) internal_name = "DMSwarm_CellID";
            else if (strcasecmp(field_name, "Migration Status") == 0) internal_name = "DMSwarm_location_status";
        }
        PetscCall(DMSwarmGetLocalSize(swarm_to_use, &field_local));
        PetscCall(DMSwarmGetSize(swarm_to_use, &field_total));
        if (field_total != *p_n_total) {
            LOG_ALLOW(GLOBAL, LOG_WARNING, "Field '%s' has %" PetscInt_FMT " particles, but expected %" PetscInt_FMT ". Skipping.\n",
                      field_name, field_total, *p_n_total);
            continue;
        }
        PetscCheck(field_local == nlocal, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ,
                   "Field '%s' holds %" PetscInt_FMT " particles on this rank, the solver swarm %" PetscInt_FMT ".",
                   field_name, field_local, nlocal);

        VTKFieldInfo *current = &meta->point_data_fields[meta->num_point_data_fields];
        PetscCall(DMSwarmGetField(swarm_to_use, internal_name, &bs, &type, &source));
        PetscCall(PetscStrncpy(current->name, field_name, MAX_VTK_FIELD_NAME_LENGTH));
        current->num_components = bs;
        PetscCall(PetscMalloc1(bs * kept, &current->data));
        for (PetscInt i = 0; i < kept; ++i) {
            const PetscInt from = first + i * stride;
            for (PetscInt c = 0; c < bs; ++c) {
                PetscScalar *to = &current->data[bs * i + c];
                if (type == PETSC_INT64) *to = (PetscScalar)((PetscInt64 *)source)[bs * from + c];
                else if (type == PETSC_INT) *to = (PetscScalar)((PetscInt *)source)[bs * from + c];
                else if (type == PETSC_REAL) *to = ((PetscReal *)source)[bs * from + c];
#if defined(PETSC_USE_COMPLEX)
                else if (type == PETSC_SCALAR) *to = ((PetscScalar *)source)[bs * from + c];
#endif
                else SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "VTK particle output does not support field '%s' with PETSc type %s.",
                             field_name, PetscDataTypes[type]);
            }
        }
        PetscCall(DMSwarmRestoreField(swarm_to_use, internal_name, &bs, &type, &source));
        meta->num_point_data_fields++;
    }
    PetscCall(PetscFree(fields_copy));

    PetscCall(PetscMalloc1(kept, &meta->connectivity));
    PetscCall(PetscMalloc1(kept, &meta->offsets));
    for (PetscInt i = 0; i < kept; ++i) {
        meta->connectivity[i] = meta->local_start[0] + i;
        meta->offsets[i]      = meta->local_start[0] + i + 1;
    }
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}
