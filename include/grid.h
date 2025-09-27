// in include/grid.h

#ifndef GRID_H
#define GRID_H

#include "variables.h" // Essential for SimulationContext and UserCtx definitions
#include "logging.h"  // Logging macros and definitions
#include "io.h"
#include "setup.h"    // For SetDMDAProcLayout
/**
 * @file grid.h
 * @brief Public interface for grid, solver, and metric setup routines.
 */

/**
 * @brief Orchestrates the parsing and setting of grid dimensions for all blocks.
 *
 * This function serves as the high-level entry point for defining the geometric
 * properties of each grid block in the simulation. It iterates through every
 * block defined by `simCtx->block_number`.
 *
 * For each block, it performs two key actions:
 * 1.  It explicitly sets the block's index (`_this`) in the corresponding `UserCtx`
 *     struct for the finest multigrid level. This makes the context "self-aware".
 * 2.  It calls a helper function (`ParseAndSetGridInputs`) to handle the detailed
 *     work of parsing options or files to populate the rest of the geometric
 *     properties for that specific block (e.g., `IM`, `Min_X`, `rx`).
 *
 * @param simCtx The master SimCtx, which contains the number of blocks and the
 *               UserCtx hierarchy to be configured.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode DefineAllGridDimensions(SimCtx *simCtx);

/**
 * @brief Orchestrates the creation of DMDA objects for every block and multigrid level.
 *
 * This function systematically builds the entire DMDA hierarchy. It first
 * calculates the dimensions (IM, JM, KM) for all coarse grids based on the
 * finest grid's dimensions and the semi-coarsening flags. It then iterates
 * from the coarsest to the finest level, calling a powerful helper function
 * (`InitializeSingleGridDM`) to create the DMs for each block, ensuring that
 * finer grids are properly aligned with their coarser parents for multigrid efficiency.
 *
 * @param simCtx The master SimCtx, containing the configured UserCtx hierarchy.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode InitializeAllGridDMs(SimCtx *simCtx);

/**
 * @brief Orchestrates the assignment of physical coordinates to all DMDA objects.
 *
 * This function manages the entire process of populating the coordinate vectors
 * for every DMDA across all multigrid levels and blocks. It follows a two-part
 * strategy that is essential for multigrid methods:
 *
 * 1.  **Populate Finest Level:** It first loops through each block and calls a
 *     helper (`SetFinestLevelCoordinates`) to set the physical coordinates for
 *     the highest-resolution grid (the finest multigrid level).
 * 2.  **Restrict to Coarser Levels:** It then iterates downwards from the finest
 *     level, calling a helper (`RestrictCoordinates`) to copy the coordinate
 *     values from the fine grid nodes to their corresponding parent nodes on the
 *     coarser grids. This ensures all levels represent the exact same geometry.
 *
 * @param simCtx The master SimCtx, containing the configured UserCtx hierarchy.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode AssignAllGridCoordinates(SimCtx *simCtx);


/**
 * @brief Computes the local bounding box of the grid on the current process.
 *
 * This function calculates the minimum and maximum coordinates of the local grid points owned
 * by the current MPI process and stores the computed bounding box in the provided structure.
 *
 * @param[in]  user      Pointer to the user-defined context containing grid information.
 * @param[out] localBBox Pointer to the BoundingBox structure to store the computed bounding box.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode ComputeLocalBoundingBox(UserCtx *user, BoundingBox *localBBox);

/**
 * @brief Gathers local bounding boxes from all MPI processes to rank 0.
 *
 * This function computes the local bounding box on each process, then collects all local
 * bounding boxes on the root process (rank 0) using MPI. The result is stored in an array
 * of BoundingBox structures on rank 0.
 *
 * @param[in]  user       Pointer to the user-defined context containing grid information.
 * @param[out] allBBoxes  Pointer to a pointer where the array of gathered bounding boxes
 *                        will be stored on rank 0. The caller on rank 0 must free this array.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode GatherAllBoundingBoxes(UserCtx *user, BoundingBox **allBBoxes);

/**
 * @brief Broadcasts the bounding box information collected on rank 0 to all other ranks.
 *
 * This function assumes that `GatherAllBoundingBoxes()` was previously called, so `bboxlist`
 * is allocated and populated on rank 0. All other ranks will allocate memory for `bboxlist`,
 * and this function will use MPI_Bcast to distribute the bounding box data to them.
 *
 * @param[in]     user      Pointer to the UserCtx structure. (Currently unused in this function, but kept for consistency.)
 * @param[in,out] bboxlist  Pointer to the array of BoundingBoxes. On rank 0, this should point to
 *                          a valid array of size 'size' (where size is the number of MPI ranks).
 *                          On non-root ranks, this function will allocate memory for `bboxlist`.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on MPI or PETSc-related errors.
 */
PetscErrorCode BroadcastAllBoundingBoxes(UserCtx *user, BoundingBox **bboxlist);

/**
 * @brief Calculates the geometric center of the primary inlet face.
 *
 * This function identifies the first face designated as an INLET in the boundary
 * condition configuration. It then iterates over all grid nodes on that physical
 * face across all MPI processes, calculates the average of their coordinates,
 * and stores the result in the user's SimCtx (CMx_c, CMy_c, CMz_c).
 *
 * This provides an automatic, robust way to determine the center for profiles
 * like parabolic flow, removing the need for manual user input.
 *
 * @param user The main UserCtx struct, containing BC config and the grid coordinate vector.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode CalculateInletCenter(UserCtx *user);

#endif // GRID_H
