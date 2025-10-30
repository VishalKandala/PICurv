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
 * @brief Calculates the center and area of the primary INLET face.
 *
 * This function identifies the primary INLET face from the boundary face
 * configurations, computes its geometric center and total area using a
 * generic utility function, and stores these results in the simulation context.
 *
 * @param user Pointer to the UserCtx containing boundary face information.
 * @return PetscErrorCode
 */
PetscErrorCode CalculateInletProperties(UserCtx *user);

/**
 * @brief Calculates the center and area of the primary OUTLET face.
 *
 * This function identifies the primary OUTLET face from the boundary face
 * configurations, computes its geometric center and total area using a
 * generic utility function, and stores these results in the simulation context.
 * @param user Pointer to the UserCtx containing boundary face information.
 * @return PetscErrorCode
 */
PetscErrorCode CalculateOutletProperties(UserCtx *user);

/**
 * @brief Calculates the geometric center and total area of a specified boundary face.
 *
 * This function computes two key properties of a boundary face in the computational domain:
 * 1. **Geometric Center**: The average (x,y,z) position of all physical nodes on the face
 * 2. **Total Area**: The sum of face area vector magnitudes from all non-solid cells adjacent to the face
 *
 * @section architecture Indexing Architecture
 * 
 * The solver uses different indexing conventions for different field types:
 * 
 * **Node-Centered Fields (Coordinates):**
 * - Direct indexing: Node n stored at coor[n]
 * - For mx=26: Physical nodes [0-24], Dummy at [25]
 * - For mz=98: Physical nodes [0-96], Dummy at [97]
 * 
 * **Face-Centered Fields (Metrics: csi, eta, zet):**
 * - Direct indexing: Face n stored at csi/eta/zet[n]
 * - For mx=26: Physical faces [0-24], Dummy at [25]
 * - For mz=98: Physical faces [0-96], Dummy at [97]
 * - Face at index k bounds cells k-1 and k
 * 
 * **Cell-Centered Fields (nvert):**
 * - Shifted indexing: Physical cell c stored at nvert[c+1]
 * - For mx=26 (25 cells): Cell 0→nvert[1], Cell 23→nvert[24]
 * - For mz=98 (96 cells): Cell 0→nvert[1], Cell 95→nvert[96]
 * - nvert[0] and nvert[mx-1] are ghost values
 *
 * @section face_geometry Face-to-Index Mapping
 * 
 * Example for a domain with mx=26, my=26, mz=98:
 * 
 * | Face ID       | Node Index | Face Metric      | Adjacent Cell (shifted) | Physical Extent |
 * |---------------|------------|------------------|-------------------------|-----------------|
 * | BC_FACE_NEG_X | i=0        | csi[k][j][0]     | nvert[k][j][1] (Cell 0) | j∈[0,24], k∈[0,96] |
 * | BC_FACE_POS_X | i=24       | csi[k][j][24]    | nvert[k][j][24] (Cell 23)| j∈[0,24], k∈[0,96] |
 * | BC_FACE_NEG_Y | j=0        | eta[k][0][i]     | nvert[k][1][i] (Cell 0) | i∈[0,24], k∈[0,96] |
 * | BC_FACE_POS_Y | j=24       | eta[k][24][i]    | nvert[k][24][i] (Cell 23)| i∈[0,24], k∈[0,96] |
 * | BC_FACE_NEG_Z | k=0        | zet[0][j][i]     | nvert[1][j][i] (Cell 0) | i∈[0,24], j∈[0,24] |
 * | BC_FACE_POS_Z | k=96       | zet[96][j][i]    | nvert[96][j][i] (Cell 95)| i∈[0,24], j∈[0,24] |
 *
 * @section algorithm Algorithm
 * 
 * The function performs two separate computations with different loop bounds:
 * 
 * **1. Center Calculation (uses ALL physical nodes):**
 * - Loop over all physical nodes on the face (excluding dummy indices)
 * - Accumulate coordinate sums: Σx, Σy, Σz
 * - Count number of nodes
 * - Average: center = (Σx/n, Σy/n, Σz/n)
 * 
 * **2. Area Calculation (uses INTERIOR cells only):**
 * - Loop over interior cell range to avoid accessing ghost values in nvert
 * - For each face adjacent to a fluid cell (nvert < 0.1):
 *   - Compute area magnitude: |csi/eta/zet| = √(x² + y² + z²)
 *   - Accumulate to total area
 * 
 * @section loop_bounds Loop Bound Details
 * 
 * **Why different bounds for center vs. area?**
 * 
 * For BC_FACE_NEG_X at i=0 with my=26, mz=98:
 * 
 * *Center calculation (coordinates):*
 * - j ∈ [ys, j_max): Includes j=[0,24] (25 nodes), excludes dummy at j=25
 * - k ∈ [zs, k_max): Includes k=[0,96] (97 nodes), excludes dummy at k=97
 * - Total: 25 × 97 = 2,425 nodes
 * 
 * *Area calculation (nvert checks):*
 * - j ∈ [lys, lye): j=[1,24] (24 values), excludes boundaries
 * - k ∈ [lzs, lze): k=[1,96] (96 values), excludes boundaries
 * - Why restricted?
 *   - At j=0: nvert[k][0][1] is ghost (no cell at j=-1)
 *   - At j=25: nvert[k][25][1] is ghost (no cell at j=24, index 25 is dummy)
 *   - At k=0: nvert[0][j][1] is ghost (no cell at k=-1)
 *   - At k=97: nvert[97][j][1] is ghost (no cell at k=96, index 97 is dummy)
 * - Total: 24 × 96 = 2,304 interior cells adjacent to face
 *
 * @section area_formulas Area Calculation Formulas
 * 
 * Face area contributions are computed from metric tensor magnitudes:
 * - **i-faces (±Xi)**: Area = |csi| = √(csi_x² + csi_y² + csi_z²)
 * - **j-faces (±Eta)**: Area = |eta| = √(eta_x² + eta_y² + eta_z²)
 * - **k-faces (±Zeta)**: Area = |zet| = √(zet_x² + zet_y² + zet_z²)
 *
 * @param[in]  user        Pointer to UserCtx containing grid info, DMs, and field vectors
 * @param[in]  face_id     Enum identifying which boundary face to analyze (BC_FACE_NEG_X, etc.)
 * @param[out] face_center Pointer to Cmpnts structure to store computed geometric center (x,y,z)
 * @param[out] face_area   Pointer to PetscReal to store computed total face area
 *
 * @return PetscErrorCode Returns 0 on success, non-zero PETSc error code on failure
 *
 * @note This function uses MPI_Allreduce, so it must be called collectively by all ranks
 * @note Only ranks that own the specified boundary face contribute to the calculation
 * @note Center calculation includes ALL physical nodes on the face
 * @note Area calculation ONLY includes faces adjacent to fluid cells (nvert < 0.1)
 * @note Dummy/unused indices (e.g., k=97, j=25 for standard test case) are excluded
 *
 * @warning Assumes grid and field arrays have been properly initialized
 * @warning Incorrect face_id values will result in zero contribution from all ranks
 *
 * @see CanRankServiceFace() for determining rank ownership of boundary faces
 * @see BCFace enum for valid face_id values
 * @see LOG_FIELD_ANATOMY() for debugging field indexing
 *
 * @par Example Usage:
 * @code
 * Cmpnts inlet_center;
 * PetscReal inlet_area;
 * ierr = CalculateFaceCenterAndArea(user, BC_FACE_NEG_Z, &inlet_center, &inlet_area);
 * PetscPrintf(PETSC_COMM_WORLD, "Inlet center: (%.4f, %.4f, %.4f), Area: %.6f\n",
 *             inlet_center.x, inlet_center.y, inlet_center.z, inlet_area);
 * @endcode
 */
PetscErrorCode CalculateFaceCenterAndArea(UserCtx *user, BCFace face_id, 
                                          Cmpnts *face_center, PetscReal *face_area);
                                          
#endif // GRID_H
