#ifndef VTK_IO_H
#define VTK_IO_H

#include "variables.h"
#include "logging.h"
#include "io.h"

// --- Public Function Prototypes ---

/**
 * @brief Creates and writes a VTK file (either .vts or .vtp) from a populated metadata struct.
 * @param[in]  filename  The output file name.
 * @param[in]  meta      Pointer to a VTKMetaData structure containing all necessary fields.
 * @param[in]  comm      The MPI communicator.
 * @return 0 on success, non-zero on failure.
 */
PetscInt CreateVTKFileFromMetadata(const char *filename, const VTKMetaData *meta, MPI_Comm comm);

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
PetscErrorCode PrepareOutputCoordinates(UserCtx* user, PetscScalar** out_coords, PetscInt* out_nx, PetscInt* out_ny, PetscInt* out_nz, PetscInt* out_npoints);

/**
 * @brief Creates a C array of field data corresponding to a subsampled (legacy-style) grid.
 *
 * This function gathers a full, distributed PETSc vector to rank 0. On rank 0,
 * it then allocates a new, smaller C array and copies only the data components
 * for nodes within the range [0..IM-2, 0..JM-2, 0..KM-2]. This produces a contiguous
 * data array that perfectly matches the point ordering of the subsampled coordinates.
 * The output array is only allocated and valid on rank 0.
 *
 * @param[in]  user           The UserCtx for grid information.
 * @param[in]  field_vec      The full-sized PETSc vector containing the field data (e.g., user->P_nodal).
 * @param[in]  num_components The number of components for this field (1 for scalar, 3 for vector).
 * @param[out] out_data       On rank 0, a pointer to the newly allocated C array for the field data. NULL on other ranks.
 * @return PetscErrorCode
 */
PetscErrorCode PrepareOutputEulerianFieldData(UserCtx* user, Vec field_vec, PetscInt num_components, PetscScalar** out_data);
#endif // VTK_IO_H

