#include "vtk_io.h"
//================================================================================
//                 IMPLEMENTATION OF PRIVATE HELPERS
//================================================================================

/**
 * @brief Writes a data block in appended format for a VTK file, including a 4-byte size prefix.
 *
 * This function writes the size of the data block (in bytes) as a 4-byte integer, 
 * followed immediately by the raw data bytes. It returns 0 on success, or 1 on error.
 *
 * @param fp            File pointer to write to (must be open for binary write).
 * @param data          Pointer to the data to be written.
 * @param num_elements  Number of elements in the data array.
 * @param element_size  Size (in bytes) of each element.
 *
 * @return PetscInt  Returns 0 if successful, non-zero otherwise.
 */
static PetscErrorCode WriteVTKAppendedBlock(FILE *fp, const void *data, PetscInt num_elements, size_t element_size) {

    // Log the function call with parameters
  // LOG_ALLOW_SYNC(LOCAL,LOG_INFO,"WriteVTKAppendedBlock - Called with %d elements, each of size %zu bytes.\n",num_elements, element_size);

    // Calculate the block size
    PetscInt block_size = num_elements * (PetscInt)element_size;
// LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "WriteVTKAppendedBlock - Calculated block size: %d bytes.\n", block_size);

    // Write the block size as a 4-byte integer
    if (fwrite(&block_size, sizeof(PetscInt), 1, fp) != 1) {
        fprintf(stderr, "[ERROR] Failed to write block size.\n");
	//    LOG_ALLOW_SYNC(LOCAL, LOG_ERROR, "WriteVTKAppendedBlock - Error writing block size.\n");
        return 1;
    }

    // Write the actual data
    if (fwrite(data, element_size, num_elements, fp) != (size_t)num_elements) {
        fprintf(stderr, "[ERROR] Failed to write data block.\n");
	//   LOG_ALLOW_SYNC(LOCAL, LOG_ERROR, "WriteVTKAppendedBlock - Error writing data block.\n");
        return 1;
    }

    // Log success
// LOG_ALLOW_SYNC(LOCAL, LOG_INFO, "WriteVTKAppendedBlock - Successfully wrote block of %d bytes.\n", block_size);

    return 0; // Success
}

/**
 * @brief Writes the XML header for a .vts (structured grid) file, supporting multiple data arrays.
 *
 * @param[in]  fp         File pointer for writing.
 * @param[in]  meta       The metadata struct containing grid dimensions and field information.
 * @param[in]  boffset    The starting byte offset for the first appended data block.
 * @param[out] boffsetOut The final byte offset after accounting for all data blocks described in the header.
 * @return 0 on success.
 */
static PetscInt WriteVTSXMLHeader(FILE *fp, const VTKMetaData *meta, PetscInt boffset, PetscInt *boffsetOut)
{
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Writing multi-field .vts header for grid %D x %D x %D.\n", meta->mx, meta->my, meta->mz);

    const char *byte_order = "LittleEndian";
    const char *precision  = "Float64"; // Assumes double precision for PetscScalar

    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"%s\">\n", byte_order);
    fprintf(fp, "  <StructuredGrid WholeExtent=\"0 %D 0 %D 0 %D\">\n", meta->mx - 1, meta->my - 1, meta->mz - 1);
    fprintf(fp, "    <Piece Extent=\"0 %D 0 %D 0 %D\">\n", meta->mx - 1, meta->my - 1, meta->mz - 1);

    // --- Points Section ---
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%D\" />\n",
            precision, boffset);
    boffset += (PetscInt)sizeof(PetscInt) + 3 * meta->npoints * (PetscInt)sizeof(PetscScalar);
    fprintf(fp, "      </Points>\n");

    // --- PointData Section (with a loop for multiple fields) ---
    if (meta->num_point_data_fields > 0) {
        fprintf(fp, "      <PointData>\n");
        for (PetscInt i = 0; i < meta->num_point_data_fields; i++) {
            const VTKFieldInfo* field = &meta->point_data_fields[i];
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "  -> Adding header for field '%s' (offset %D)\n", field->name, boffset);
            fprintf(fp, "        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"%D\" format=\"appended\" offset=\"%D\" />\n",
                    precision, field->name, field->num_components, boffset);
            boffset += (PetscInt)sizeof(PetscInt) + field->num_components * meta->npoints * (PetscInt)sizeof(PetscScalar);
        }
        fprintf(fp, "      </PointData>\n");
    }

    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </StructuredGrid>\n");
    fprintf(fp, "  <AppendedData encoding=\"raw\">\n_"); // Underscore is the required separator

    *boffsetOut = boffset;
    return 0;
}

/**
 * @brief Finalizes the .vts file by closing the appended data section and the VTKFile XML tag.
 *
 * This function writes the closing XML tags for the appended data section and the 
 * overall VTK file, ensuring a well-formed .vts file structure.
 *
 * @param[in] fp  File pointer (already open) for writing the .vts footer.
 *
 * @return PetscInt  Returns 0 on success.
 */
static PetscInt WriteVTSXMLFooter(FILE *fp)
{
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTSXMLFooter - Closing .vts XML.\n");

    fprintf(fp, "\n </AppendedData>\n");
    fprintf(fp, "</VTKFile>\n");

    // Log completion
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTSXMLFooter - Successfully wrote .vts footer.\n");

    return 0;
}

/**
 * @brief Writes the XML header for a .vtp (polydata) file, supporting multiple data arrays.
 *
 * @param[in]  fp         File pointer for writing.
 * @param[in]  meta       The metadata struct containing point count and field information.
 * @param[in]  boffset    The starting byte offset for the first appended data block.
 * @param[out] boffsetOut The final byte offset after accounting for all data blocks described in the header.
 * @return 0 on success.
 */
static PetscInt WriteVTPXMLHeader(FILE *fp, const VTKMetaData *meta, PetscInt boffset, PetscInt *boffsetOut)
{
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Writing multi-field .vtp header for %D points.\n", meta->npoints);

    const char *byte_order = "LittleEndian";
    const char *precision     = "Float64";
    const char *int_type_str  = (sizeof(PetscInt) == 8) ? "Int64" : "Int32";

    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"%s\">\n", byte_order);
    fprintf(fp, "  <PolyData>\n");
    fprintf(fp, "    <Piece NumberOfPoints=\"%D\" NumberOfVerts=\"%D\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n",
            meta->npoints, meta->npoints);

    // --- Points Section ---
    fprintf(fp, "      <Points>\n");
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%D\" />\n",
            precision, boffset);
    boffset += (PetscInt)sizeof(PetscInt) + 3 * meta->npoints * (PetscInt)sizeof(PetscScalar);
    fprintf(fp, "      </Points>\n");
    
    // --- PointData Section (with a loop for multiple fields) ---
    if (meta->num_point_data_fields > 0) {
        fprintf(fp, "      <PointData>\n");
        for (PetscInt i = 0; i < meta->num_point_data_fields; i++) {
            const VTKFieldInfo* field = &meta->point_data_fields[i];
            LOG_ALLOW(GLOBAL, LOG_DEBUG, "  -> Adding header for field '%s' (offset %D)\n", field->name, boffset);
            fprintf(fp, "        <DataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"%D\" format=\"appended\" offset=\"%D\" />\n",
                    precision, field->name, field->num_components, boffset);
            boffset += (PetscInt)sizeof(PetscInt) + field->num_components * meta->npoints * (PetscInt)sizeof(PetscScalar);
        }
        fprintf(fp, "      </PointData>\n");
    }

    // --- Verts Section (for simple point cloud) ---
    fprintf(fp, "      <Verts>\n");
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"connectivity\" format=\"appended\" offset=\"%D\" />\n",
            int_type_str, boffset);
    boffset += (PetscInt)sizeof(PetscInt) + meta->npoints * (PetscInt)sizeof(PetscInt);
    fprintf(fp, "        <DataArray type=\"%s\" Name=\"offsets\" format=\"appended\" offset=\"%D\" />\n",
            int_type_str, boffset);
    boffset += (PetscInt)sizeof(PetscInt) + meta->npoints * (PetscInt)sizeof(PetscInt);
    fprintf(fp, "      </Verts>\n");

    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "  </PolyData>\n");
    fprintf(fp, "  <AppendedData encoding=\"raw\">\n_");

    *boffsetOut = boffset;
    return 0;
}

/**
 * @brief Finalizes the .vtp file by closing the appended data section and the VTKFile tag.
 *
 * This function writes the closing XML tags for the appended data section 
 * and completes the overall VTK structure, ensuring a well-formed .vtp file.
 *
 * @param[in] fp  File pointer to the .vtp file (open for writing).
 *
 * @return PetscInt  Returns 0 on success.
 */
static PetscInt WriteVTPXMLFooter(FILE *fp)
{
    // Log the function entry
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteVTPXMLFooter - Closing VTP XML.\n");

    // no stray prPetscInts or data -> directly close
    fprintf(fp, "\n  </AppendedData>\n");
    fprintf(fp, "</VTKFile>\n");

    // Log successful completion
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteVTPXMLFooter - Completed writing VTP XML footer.\n");

    return 0;
}

/**
 * @brief Writes the initial XML header for a VTK file based on the provided metadata.
 *
 * This function writes the XML header for VTK output. It now supports both polydata
 * (VTK_POLYDATA) and structured grid (VTK_STRUCTURED) file types. For each type, the
 * function determines whether to write a scalar or vector field header by inspecting
 * the metadata. For polydata, it calls \c WriteVTPXMLHeader and for structured grids,
 * it calls \c WriteVTSXMLHeader.
 *
 * @param[in]  fp          File pointer (open for writing) to which the header will be written.
 * @param[in]  meta        Pointer to a \c VTKMetaData structure containing fileType, number of points,
 *                         grid dimensions, field names, and other necessary information.
 * @param[in]  boffset     Current byte offset in the appended data section.
 * @param[out] boffsetOut  Updated byte offset after writing any header-related data.
 *
 * @return PetscInt Returns 0 on success, or -1 if the fileType is not handled.
 *
 * @note The function uses PETSc-style logging (LOG_ALLOW) to trace its execution.
 *       For structured grid output, it assumes that grid dimensions (mx, my, mz) and
 *       the number of nodes (nnodes) have been set in the metadata.
 */
static PetscInt WriteVTKFileHeader(FILE *fp, const VTKMetaData *meta, PetscInt boffset, PetscInt *boffsetOut)
{
    const char *fieldName = NULL;
    PetscInt numComponents = 1;

    /*-----------------------------------------------------------------------
     * Log the entry into this function with the key metadata parameters.
     *-----------------------------------------------------------------------*/
    LOG_ALLOW(GLOBAL, LOG_DEBUG,
              "WriteVTKFileHeader - Entered: fileType=%d, npoints=%d, scalarField=%s, vectorField=%s, "
              "numVectorFields=%d, mx=%d, my=%d, mz=%d, nnodes=%d.\n",
              meta->fileType, meta->npoints, meta->scalarFieldName, meta->vectorFieldName,
              meta->numVectorFields, meta->mx, meta->my, meta->mz, meta->nnodes);

    if (meta->fileType == VTK_POLYDATA) {
        /*===================================================================
         * Polydata Branch: VTK_POLYDATA
         *===================================================================*/
        LOG_ALLOW(GLOBAL, LOG_INFO, "WriteVTKFileHeader - Processing VTK_POLYDATA output.\n");

        /* Determine the field type for polydata:
         * If a vector field is available, use it (with 3 components);
         * otherwise, use the scalar field (with 1 component).
         */
        if (meta->numVectorFields > 0 && meta->vectorFieldName) {
            fieldName = meta->vectorFieldName;
            numComponents = 3;
            LOG_ALLOW(GLOBAL, LOG_INFO,
                      "WriteVTKFileHeader - Using vector field '%s' with 3 components for polydata.\n",
                      fieldName);
        } else {
            fieldName = meta->scalarFieldName;
            numComponents = 1;
            LOG_ALLOW(GLOBAL, LOG_INFO,
                      "WriteVTKFileHeader - Using scalar field '%s' with 1 component for polydata.\n",
                      fieldName);
        }

        /* Log the call to the polydata header-writing function */
        LOG_ALLOW(GLOBAL, LOG_DEBUG,
                  "WriteVTKFileHeader - Calling WriteVTPXMLHeader with boffset=%d.\n", boffset);
        return WriteVTPXMLHeader(fp,
                                 meta,
                                 boffset,
                                 boffsetOut);
    } else if (meta->fileType == VTK_STRUCTURED) {
        /*===================================================================
         * Structured Grid Branch: VTK_STRUCTURED
         *===================================================================*/
        LOG_ALLOW(GLOBAL, LOG_INFO, "WriteVTKFileHeader - Processing VTK_STRUCTURED output.\n");

        /* Determine the field type for structured grid:
         * If a vector field is available, use it (with 3 components);
         * otherwise, use the scalar field (with 1 component).
         */
        if (meta->numVectorFields > 0 && meta->vectorFieldName) {
            fieldName = meta->vectorFieldName;
            numComponents = 3;
            LOG_ALLOW(GLOBAL, LOG_INFO,
                      "WriteVTKFileHeader - Using vector field '%s' with 3 components for structured grid.\n",
                      fieldName);
        } else {
            fieldName = meta->scalarFieldName;
            numComponents = 1;
            LOG_ALLOW(GLOBAL, LOG_INFO,
                      "WriteVTKFileHeader - Using scalar field '%s' with 1 component for structured grid.\n",
                      fieldName);
        }

        /* Log the grid dimensions and number of nodes, then call the structured grid header-writing function */
        LOG_ALLOW(GLOBAL, LOG_DEBUG,
                  "WriteVTKFileHeader - Calling WriteVTSXMLHeader with boffset=%d, grid dimensions (%d, %d, %d) and nnodes=%d.\n",
                  boffset, meta->mx, meta->my, meta->mz, meta->nnodes);
        return WriteVTSXMLHeader(fp,
                                 meta,
                                 boffset,
                                 boffsetOut);
    } else {
        /*-------------------------------------------------------------------
         * Unsupported File Type: Log a warning and return an error.
         *-------------------------------------------------------------------*/
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "WriteVTKFileHeader - Unsupported file type %d encountered.\n", meta->fileType);
        return -1;
    }
}

/**
 * @brief Writes the XML footer for a VTK file based on the provided metadata.
 *
 * This function currently handles only \c VTK_POLYDATA. If the file type is 
 * \c VTK_POLYDATA, it delegates to \c WriteVTPXMLFooter. Otherwise, 
 * it logs a warning and returns -1 (indicating that other file types are not yet implemented).
 *
 * @param[in] fp    File pointer (open for writing) to which the footer will be written.
 * @param[in] meta  Pointer to a \c VTKMetaData structure containing the fileType 
 *                  (e.g., \c VTK_POLYDATA or \c VTK_STRUCTURED).
 *
 * @return PetscInt  Returns 0 or the value from \c WriteVTPXMLFooter on success, 
 *              or -1 if the file type is not handled.
 */
static PetscInt WriteVTKFileFooter(FILE *fp, const VTKMetaData *meta)
{
    // Log the entry into this function
    LOG_ALLOW(GLOBAL, LOG_DEBUG, 
              "WriteVTKFileFooter - Called with fileType=%d.\n", meta->fileType);

    if (meta->fileType == VTK_POLYDATA) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteVTKFileFooter - Writing VTP XML footer.\n");
        return WriteVTPXMLFooter(fp);
    } else {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "WriteVTKFileFooter - Writing VTS XML footer.\n");
	return WriteVTSXMLFooter(fp);
    }
}

//================================================================================
//                 PUBLIC VTK WRITER FUNCTIONS
//================================================================================

/**
 * @brief Creates and writes a VTK file (either .vts or .vtp) based on the provided metadata.
 *
 * This function orchestrates the writing process. It calls helper functions to write
 * the XML header, then writes all the binary data blocks in the correct order, and
 * finally writes the XML footer. The I/O occurs only on rank 0.
 *
 * @param[in]  filename  The output file name (e.g., "output.vtp" or "output.vts").
 * @param[in]  meta      Pointer to a VTKMetaData structure containing all necessary fields.
 * @param[in]  comm      The MPI communicator used for parallel execution.
 * @return 0 on success, or 1 on file-opening failure.
 */
PetscInt CreateVTKFileFromMetadata(const char *filename, const VTKMetaData *meta, MPI_Comm comm)
{
    PetscMPIInt rank;
    MPI_Comm_rank(comm, &rank);

    // Only rank 0 performs file I/O.
    if (!rank) {
        LOG_ALLOW(GLOBAL, LOG_INFO, "Rank 0 writing combined VTK file '%s'.\n", filename);

        FILE *fp = fopen(filename, "wb");
        if (!fp) {
            LOG_ALLOW(GLOBAL, LOG_ERROR, "fopen failed for %s.\n", filename);
            return 1; // Return an error code
        }

        PetscInt boffset = 0;
        
        // --- 1. Write the appropriate XML header ---
        WriteVTKFileHeader(fp, meta, boffset, &boffset);

        // --- 2. Write all appended data blocks in the EXACT order they were declared in the header ---
        
        // (a) Coordinates
        if (meta->coords) {
            WriteVTKAppendedBlock(fp, meta->coords, 3 * meta->npoints, sizeof(PetscScalar));
        }

        // (b) All Point Data fields
        for (PetscInt i = 0; i < meta->num_point_data_fields; i++) {
            const VTKFieldInfo* field = &meta->point_data_fields[i];
            if (field->data) {
                WriteVTKAppendedBlock(fp, field->data, field->num_components * meta->npoints, sizeof(PetscScalar));
            }
        }

        // (c) Connectivity and Offsets (only for PolyData)
        if (meta->fileType == VTK_POLYDATA) {
            if (meta->connectivity) {
                WriteVTKAppendedBlock(fp, meta->connectivity, meta->npoints, sizeof(PetscInt));
            }
            if (meta->offsets) {
                WriteVTKAppendedBlock(fp, meta->offsets, meta->npoints, sizeof(PetscInt));
            }
        }

        // --- 3. Write the XML footer ---
        WriteVTKFileFooter(fp, meta);
        
        fclose(fp);
        LOG_ALLOW(GLOBAL, LOG_INFO, "Successfully wrote and closed file '%s'.\n", filename);
    }

    return 0; // Success
}
