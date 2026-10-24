#ifndef VTK_IO_H
#define VTK_IO_H

#include "variables.h"
#include "logging.h"
#include "io.h"

// --- Public Function Prototypes ---

/**
 * @brief Writes one VTK file (.vts or .vtp) collectively from every rank's part of it.
 * @details Rank 0 writes the XML header, each block's byte count and the footer; every rank
 *          writes its own points into each appended block through MPI-IO, so no rank holds
 *          more than its own part. The file is the same whatever the rank count. Block byte
 *          counts are UInt64 (`header_type="UInt64"`), so no block size limit applies.
 * @param[in]  filename  The output file name.
 * @param[in]  meta      This rank's part of the file (see VTKMetaData).
 * @param[in]  comm      The communicator all of whose ranks call this.
 * @return 0 on success, or a PETSc error.
 */
PetscErrorCode CreateVTKFileFromMetadata(const char *filename, const VTKMetaData *meta, MPI_Comm comm);

/**
 * @brief Copies this rank's part of the output grid's coordinates.
 *
 * The output grid is the nodes [0..IM-2, 0..JM-2, 0..KM-2], a grid of (IM-1)x(JM-1)x(KM-1)
 * points. This rank's part is its owned nodes within it, i fastest.
 *
 * @param[in]  user        The UserCtx containing the grid information (DM, IM/JM/KM).
 * @param[out] out_coords  Newly allocated coordinates of this rank's points (x, y, z each).
 * @param[out] out_nx      The number of points in the x-dimension for the new grid (IM-1).
 * @param[out] out_ny      The number of points in the y-dimension for the new grid (JM-1).
 * @param[out] out_nz      The number of points in the z-dimension for the new grid (KM-1).
 * @param[out] out_npoints The total number of points in the new grid.
 * @param[out] local_start First output point of this rank's box on each axis.
 * @param[out] local_size  Points in this rank's box on each axis (zero for a rank holding only dummy nodes).
 * @return PetscErrorCode
 */
PetscErrorCode PrepareOutputCoordinates(UserCtx* user, PetscScalar** out_coords, PetscInt* out_nx, PetscInt* out_ny, PetscInt* out_nz,
                                        PetscInt* out_npoints, PetscInt local_start[3], PetscInt local_size[3]);

/**
 * @brief Copies this rank's part of a nodal field on the output grid.
 *
 * The points and their order are those of PrepareOutputCoordinates() on the same rank, with
 * the field's components interleaved.
 *
 * @param[in]  user           The UserCtx for grid information.
 * @param[in]  field_vec      The full-sized PETSc vector containing the field data (e.g., user->P_nodal).
 * @param[in]  num_components The number of components for this field (1 for scalar, 3 for vector).
 * @param[out] out_data       Newly allocated values of this rank's points.
 * @return PetscErrorCode
 */
PetscErrorCode PrepareOutputEulerianFieldData(UserCtx* user, Vec field_vec, PetscInt num_components, PetscScalar** out_data);

/**
 * @brief Begins one structured VTK file: clears the metadata and builds its coordinates.
 *
 * Every Eulerian VTK file is assembled the same way — coordinates, then a set of
 * point-data fields, then the write. This trio exists so the two producers share that
 * shape instead of restating it, and so the buffers each producer allocates are
 * released in one place rather than in neither.
 *
 * @param[in]  user Block context supplying the grid.
 * @param[out] meta Metadata to initialize.
 * @return Zero on success, or a PETSc error.
 */
PetscErrorCode BeginStructuredVTKOutput(UserCtx *user, VTKMetaData *meta);

/**
 * @brief Adds one point-data field to a structured VTK file being assembled.
 *
 * Skips the field with a warning when the file is already at its field limit, so a
 * long output list degrades rather than corrupting the metadata.
 *
 * @param[in]     user       Block context supplying the grid.
 * @param[in,out] meta       Metadata to append to.
 * @param[in]     name       Field name as it will appear in the file.
 * @param[in]     field_vec  Nodal vector holding the values.
 * @param[in]     components One for a scalar, three for a vector.
 * @return Zero on success, or a PETSc error.
 */
PetscErrorCode AppendStructuredVTKField(UserCtx *user, VTKMetaData *meta, const char *name,
                                        Vec field_vec, PetscInt components);

/**
 * @brief Writes an assembled structured VTK file and releases its buffers.
 * @param[in,out] meta     Metadata to write and then release.
 * @param[in]     filename Destination path.
 * @return Zero on success, or a PETSc error.
 */
PetscErrorCode FinishStructuredVTKOutput(VTKMetaData *meta, const char *filename);

/**
 * @brief Prepares this rank's part of the subsampled particle VTK output.
 *
 * Collective: all ranks must enter it. Particles are numbered in rank order and every
 * pps->particle_output_freq-th one is kept; each rank copies only its own kept particles
 * (positions, the requested fields, connectivity and offsets), with no gather.
 *
 * @param[in]  user      The UserCtx containing the DMSwarm.
 * @param[in]  pps       The PostProcessParams struct for configuration.
 * @param[out] meta      Zeroed metadata to populate with this rank's part (see VTKMetaData).
 * @param[out] p_n_total The total number of particles before subsampling, on every rank.
 * @return PetscErrorCode
 */
PetscErrorCode PrepareOutputParticleData(UserCtx* user, PostProcessParams* pps, VTKMetaData* meta, PetscInt* p_n_total);

#endif // VTK_IO_H

