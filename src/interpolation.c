/**
 * @file Interpolatetion.c
 * @brief Main program for DMSwarm Interpolatetion using the fdf-curvIB method.
 *
 * Provides routines for Interpolatetion between corner-based and center-based
 * fields in the cell-centered DM (fda), plus partial usage examples for
 * DMSwarm-based field sampling.
 */

#include "interpolation.h"

// Number of weights used in certain trilinear Interpolatetion examples
#define NUM_WEIGHTS 8
// Define a buffer size for error messages if not already available
#ifndef ERROR_MSG_BUFFER_SIZE
#define ERROR_MSG_BUFFER_SIZE 256 // Or use PETSC_MAX_PATH_LEN if appropriate
#endif


/**
 * @brief Interpolates a vector field from corner nodes to cell centers.
 *
 * This function estimates the value of a vector field at the center of each
 * computational cell (whose origin node is owned by this rank) by averaging
 * the values from the 8 corner nodes defining that cell.
 *
 * @param[in]  field_arr       Input: 3D array view (ghosted) of vector data at nodes,
 *                             obtained from a Vec associated with user->fda. Accessed via GLOBAL node indices.
 *                             Ghost values must be up-to-date.
 * @param[out] centfield_arr   Output: 3D local array (0-indexed) where interpolated cell center values
 *                             are stored. Sized [zm_cell_local][ym_cell_local][xm_cell_local]
 *                             by the caller. centfield_arr[k_loc][j_loc][i_loc] stores the value for
 *                             the cell corresponding to the (i_loc, j_loc, k_loc)-th owned cell origin.
 * @param[in]  user            User context containing DMDA information (fda, IM, JM, KM).
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateFieldFromCornerToCenter_Vector(
    Cmpnts ***field_arr,     /* Input: Ghosted local array view from user->fda (global node indices) */
    Cmpnts ***centfield_arr, /* Output: Local 0-indexed array */
    UserCtx *user)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
                   "Rank %d starting InterpolateFieldFromCornerToCenter_Vector.\n", rank);

    // Get local info based on the NODE-based DMDA (user->fda)
    // This DM defines the ownership of the nodes from which we interpolate,
    // and thus defines the cells whose centers we will compute.
    DMDALocalInfo info_nodes_fda;
    ierr = DMDAGetLocalInfo(user->fda, &info_nodes_fda); CHKERRQ(ierr);

    // Determine owned CELL ranges based on owning their origin NODES (from user->fda)
    PetscInt xs_cell_global_i, xm_cell_local_i; // Global start cell index, local count of owned cells
    PetscInt ys_cell_global_j, ym_cell_local_j;
    PetscInt zs_cell_global_k, zm_cell_local_k;

    
    ierr = GetOwnedCellRange(&info_nodes_fda, 0, &xs_cell_global_i, &xm_cell_local_i); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&info_nodes_fda, 1, &ys_cell_global_j, &ym_cell_local_j); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&info_nodes_fda, 2, &zs_cell_global_k, &zm_cell_local_k); CHKERRQ(ierr);

    // Exclusive end global cell indices
    PetscInt xe_cell_global_i_excl = xs_cell_global_i + xm_cell_local_i;
    PetscInt ye_cell_global_j_excl = ys_cell_global_j + ym_cell_local_j;
    PetscInt ze_cell_global_k_excl = zs_cell_global_k + zm_cell_local_k;

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d : Processing Owned Cells (global indices) k=%d..%d (count %d), j=%d..%d (count %d), i=%d..%d (count %d)\n",
              rank,
              zs_cell_global_k, ze_cell_global_k_excl -1, zm_cell_local_k,
              ys_cell_global_j, ye_cell_global_j_excl -1, ym_cell_local_j,
              xs_cell_global_i, xe_cell_global_i_excl -1, xm_cell_local_i);

    // Loop over the GLOBAL indices of the CELLS whose origin nodes are owned by this processor
    // (according to user->fda partitioning)
    for (PetscInt k_glob_cell = zs_cell_global_k; k_glob_cell < ze_cell_global_k_excl; k_glob_cell++) {
        PetscInt k_local = k_glob_cell - zs_cell_global_k; // 0-based local index for centfield_arr

        for (PetscInt j_glob_cell = ys_cell_global_j; j_glob_cell < ye_cell_global_j_excl; j_glob_cell++) {
            PetscInt j_local = j_glob_cell - ys_cell_global_j; // 0-based local index

            for (PetscInt i_glob_cell = xs_cell_global_i; i_glob_cell < xe_cell_global_i_excl; i_glob_cell++) {
                PetscInt i_local = i_glob_cell - xs_cell_global_i; // 0-based local index

                Cmpnts sum = {0.0, 0.0, 0.0};
                // Count is always 8 for a hexahedral cell
                // PetscInt count = 0; // Not strictly needed if we assume 8 corners

                // Loop over the 8 corner NODES that define cell C(i_glob_cell, j_glob_cell, k_glob_cell)
                for (PetscInt dk = 0; dk < 2; dk++) {
                    for (PetscInt dj = 0; dj < 2; dj++) {
                        for (PetscInt di = 0; di < 2; di++) {
                            PetscInt ni_glob = i_glob_cell + di; // Global NODE index i or i+1
                            PetscInt nj_glob = j_glob_cell + dj; // Global NODE index j or j+1
                            PetscInt nk_glob = k_glob_cell + dk; // Global NODE index k or k+1

                            // Access the input 'field_arr' (node values) using GLOBAL node indices.
                            // Ghosts in field_arr must be valid.
                            // DMDAVecGetArrayRead handles mapping global indices to local memory.
                            sum.x += field_arr[nk_glob][nj_glob][ni_glob].x;
                            sum.y += field_arr[nk_glob][nj_glob][ni_glob].y;
                            sum.z += field_arr[nk_glob][nj_glob][ni_glob].z;
                            // count++;
                        }
                    }
                }

                // Calculate average value representing the cell center
                Cmpnts center_value;
                center_value.x = sum.x / 8.0;
                center_value.y = sum.y / 8.0;
                center_value.z = sum.z / 8.0;

                // Store result into the local, 0-indexed output array 'centfield_arr'
                // Defensive check (should not be strictly necessary if caller allocated centfield_arr
                // with dimensions xm_cell_local_i, ym_cell_local_j, zm_cell_local_k)
                if (i_local < 0 || i_local >= xm_cell_local_i ||
                    j_local < 0 || j_local >= ym_cell_local_j ||
                    k_local < 0 || k_local >= zm_cell_local_k) {
                    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Local index out of bounds for centfield_arr write!");
                }
                centfield_arr[k_local][j_local][i_local] = center_value;

            } // End loop i_glob_cell
        } // End loop j_glob_cell
    } // End loop k_glob_cell

    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
                   "Rank %d completed InterpolateFieldFromCornerToCenter_Vector.\n", rank);
    return 0;
}

// -----------------------------------------------------------------------------
// Interpolation: Corner -> Center (scalar)
// -----------------------------------------------------------------------------

/**
 * @brief Interpolates a scalar field from corner nodes to cell centers.
 *
 * This function estimates the value of a scalar field at the center of each
 * computational cell (whose origin node is owned by this rank) by averaging
 * the values from the 8 corner nodes defining that cell.
 *
 * @param[in]  field_arr       Input: 3D array view (ghosted) of scalar data at nodes,
 *                             obtained from a Vec associated with user->da. Accessed via GLOBAL node indices.
 *                             Ghost values must be up-to-date.
 * @param[out] centfield_arr   Output: 3D local array (0-indexed) where interpolated cell center values
 *                             are stored. Sized [zm_cell_local][ym_cell_local][xm_cell_local]
 *                             by the caller. centfield_arr[k_loc][j_loc][i_loc] stores the value for
 *                             the cell corresponding to the (i_loc, j_loc, k_loc)-th owned cell origin.
 * @param[in]  user            User context containing DMDA information (da, IM, JM, KM).
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateFieldFromCornerToCenter_Scalar(
    PetscReal ***field_arr,     /* Input: Ghosted local array view from user->da (global node indices) */
    PetscReal ***centfield_arr, /* Output: Local 0-indexed array for cell centers */
    UserCtx *user)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
                   "Rank %d starting InterpolateFieldFromCornerToCenter_Scalar.\n", rank);

    // Get local info based on user->da (which is node-based but defines cell origins)
    DMDALocalInfo info_da_nodes;
    ierr = DMDAGetLocalInfo(user->da, &info_da_nodes); CHKERRQ(ierr);

    // Determine owned CELL ranges based on owning their origin NODES (from user->da)
    PetscInt xs_cell_global_i, xm_cell_local_i; // Global start cell index, local count of owned cells
    PetscInt ys_cell_global_j, ym_cell_local_j;
    PetscInt zs_cell_global_k, zm_cell_local_k;

    ierr = GetOwnedCellRange(&info_da_nodes, 0, &xs_cell_global_i, &xm_cell_local_i); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&info_da_nodes, 1, &ys_cell_global_j, &ym_cell_local_j); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&info_da_nodes, 2, &zs_cell_global_k, &zm_cell_local_k); CHKERRQ(ierr);

    // Exclusive end global cell indices
    PetscInt xe_cell_global_i_excl = xs_cell_global_i + xm_cell_local_i;
    PetscInt ye_cell_global_j_excl = ys_cell_global_j + ym_cell_local_j;
    PetscInt ze_cell_global_k_excl = zs_cell_global_k + zm_cell_local_k;

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d (IFCTCS): Processing Owned Cells (global indices) k=%d..%d (count %d), j=%d..%d (count %d), i=%d..%d (count %d)\n",
              rank,
              zs_cell_global_k, ze_cell_global_k_excl -1, zm_cell_local_k,
              ys_cell_global_j, ye_cell_global_j_excl -1, ym_cell_local_j,
              xs_cell_global_i, xe_cell_global_i_excl -1, xm_cell_local_i);

    // Loop over the GLOBAL indices of the CELLS whose origin nodes are owned by this processor
    // (according to user->da partitioning)
    for (PetscInt k_glob_cell = zs_cell_global_k; k_glob_cell < ze_cell_global_k_excl; k_glob_cell++) {
        PetscInt k_local = k_glob_cell - zs_cell_global_k; // 0-based local index for centfield_arr

        for (PetscInt j_glob_cell = ys_cell_global_j; j_glob_cell < ye_cell_global_j_excl; j_glob_cell++) {
            PetscInt j_local = j_glob_cell - ys_cell_global_j; // 0-based local index

            for (PetscInt i_glob_cell = xs_cell_global_i; i_glob_cell < xe_cell_global_i_excl; i_glob_cell++) {
                PetscInt i_local = i_glob_cell - xs_cell_global_i; // 0-based local index

                PetscReal sum = 0.0;
                // Count is always 8 for a hexahedral cell if all nodes are accessible
                // PetscInt count = 0;

                // Loop over the 8 NODE indices that define cell C(i_glob_cell, j_glob_cell, k_glob_cell)
                for (PetscInt dk = 0; dk < 2; dk++) {
                    for (PetscInt dj = 0; dj < 2; dj++) {
                        for (PetscInt di = 0; di < 2; di++) {
                            PetscInt ni_glob = i_glob_cell + di; // Global NODE index for corner
                            PetscInt nj_glob = j_glob_cell + dj; // Global NODE index for corner
                            PetscInt nk_glob = k_glob_cell + dk; // Global NODE index for corner

                            // Access the input 'field_arr' (node values) using GLOBAL node indices.
                            sum += field_arr[nk_glob][nj_glob][ni_glob];
                            // count++;
                        }
                    }
                }

                PetscReal center_value = sum / 8.0;

                // Store the calculated cell center value into the output array 'centfield_arr'
                // using the LOCAL 0-based cell index.
                // Defensive check (should not be strictly necessary if caller allocated centfield_arr correctly)
                if (i_local < 0 || i_local >= xm_cell_local_i ||
                    j_local < 0 || j_local >= ym_cell_local_j ||
                    k_local < 0 || k_local >= zm_cell_local_k) {
                    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Local index out of bounds for centfield_arr write in Scalar version!");
                }
                centfield_arr[k_local][j_local][i_local] = center_value;

            } // End loop i_glob_cell
        } // End loop j_glob_cell
    } // End loop k_glob_cell

    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
                   "Rank %d completed InterpolateFieldFromCornerToCenter_Scalar.\n", rank);
    return 0;
}

// -----------------------------------------------------------------------------
// Interpolation: Center -> Corner (scalar)
// -----------------------------------------------------------------------------
/**
 * @brief Interpolates a scalar field from cell centers to corner nodes.
 *
 * This function estimates the value of a scalar field at each grid node by averaging
 * the values from the cell centers of the cells surrounding that node (up to 8).
 * It handles physical boundaries by averaging only the available adjacent cells.
 *
 * Assumes input `centfield_arr` is from a ghosted local vector associated with `user->da` (DOF=1, s=2)
 * and output `field_arr` is from a ghosted local vector also associated with `user->da` (DOF=1, s=2).
 * Input array uses GLOBAL cell indices, output array uses GLOBAL node indices.
 *
 * @param[in]  centfield_arr  Input: 3D array (ghosted) of scalar data at cell centers,
 *                            accessed via GLOBAL cell indices (k=0..KM-1, j=0..JM-1, i=0..IM-1).
 * @param[out] field_arr      Output: 3D array (ghosted) where interpolated node values are stored,
 *                            accessed via GLOBAL node indices (k=0..KM, j=0..JM, i=0..IM).
 * @param[in]  user           User context containing DMDA information (da).
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateFieldFromCenterToCorner_Scalar(
    PetscReal ***centfield_arr, /* Input: Ghosted local array based on da (cell indices) */
    PetscReal ***field_arr,     /* Output: Ghosted local array based on da (node indices) */
    UserCtx *user)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
      "Rank %d starting interpolation.\n", rank);

    // Get local info based on the DMDA (da). This info primarily describes owned nodes.
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);

    // Define the range of NODES owned by this processor using GLOBAL indices
    PetscInt xs_node = info.xs;
    PetscInt xm_node = info.xm;
    PetscInt xe_node = xs_node + xm_node;
    PetscInt ys_node = info.ys;
    PetscInt ym_node = info.ym;
    PetscInt ye_node = ys_node + ym_node;
    PetscInt zs_node = info.zs;
    PetscInt zm_node = info.zm;
    PetscInt ze_node = zs_node + zm_node;

    // Get Global dimensions (number of cells IM, JM, KM)
    PetscInt IM = info.mx - 1;
    PetscInt JM = info.my - 1;
    PetscInt KM = info.mz - 1;


    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Rank %d: Interpolating for owned NODES k=%d..%d, j=%d..%d, i=%d..%d\n",
              rank, zs_node, ze_node, ys_node, ye_node, xs_node, xe_node);

    // Loop over the GLOBAL indices of the NODES owned by this processor
    for (PetscInt k = zs_node; k < ze_node; k++) { // Global NODE index k (0..KM)
        for (PetscInt j = ys_node; j < ye_node; j++) { // Global NODE index j (0..JM)
            for (PetscInt i = xs_node; i < xe_node; i++) { // Global NODE index i (0..IM)

                PetscReal sum = 0.0;
                PetscInt count = 0;

                // Loop over the potential 8 CELL indices surrounding node (i,j,k)
                // Cell indices are (i & i-1), (j & j-1), (k & k-1)
                for (PetscInt dk = -1; dk <= 0; dk++) { // Relative cell k-index offset
                    for (PetscInt dj = -1; dj <= 0; dj++) { // Relative cell j-index offset
                        for (PetscInt di = -1; di <= 0; di++) { // Relative cell i-index offset

                            PetscInt ci = i + di; // Global CELL index of neighbor
                            PetscInt cj = j + dj; // Global CELL index of neighbor
                            PetscInt ck = k + dk; // Global CELL index of neighbor

                            // Check if this CELL index is within the valid global cell range (0..IM-1, etc.)
                            if (ci >= 0 && ci <= IM-1 &&
                                cj >= 0 && cj <= JM-1 &&
                                ck >= 0 && ck <= KM-1)
                            {
                                // Access the input 'centfield_arr' using GLOBAL cell indices.
                                // Relies on centfield_arr being from a ghosted local vector.
                                sum += centfield_arr[ck][cj][ci];
                                count++;
                                // Optional Debug Log
                                // LOG_LOOP_ALLOW(GLOBAL, LOG_DEBUG, i*j*k, 10000,
                                // " Rank %d | Node(i,j,k)=%d,%d,%d | Using Cell(ci,cj,ck)=%d,%d,%d | Value=%.3f | Count=%d\n",
                                // rank, i, j, k, ci, cj, ck, centfield_arr[ck][cj][ci], count);
                            }
                        } // end di loop
                    } // end dj loop
                } // end dk loop

                // Calculate average and store in the output 'field_arr' at the NODE index (i,j,k)
                if (count > 0) {
                    // Store the result using the GLOBAL node index [k][j][i]
                    field_arr[k][j][i] = sum / (PetscReal)count;
                } else {
                    // This indicates an issue - a node should always be adjacent to at least one cell
                    LOG_ALLOW(GLOBAL, LOG_ERROR,
                              "Rank %d: Node (i=%d,j=%d,k=%d) had count=0 surrounding cells! Check logic/ghosting.\n", rank, i, j, k);
                    // Assign a default value or handle error
                    field_arr[k][j][i] = 0.0; // Defaulting to zero might hide issues
                }
            } // End loop i (nodes)
        } // End loop j (nodes)
    } // End loop k (nodes)

    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
              "Rank %d completed interpolation.\n", rank);
    return 0;
}


// -----------------------------------------------------------------------------
// Interpolation: Center -> Corner (vector)
// -----------------------------------------------------------------------------

/**
 * @brief Interpolates a vector field from cell centers to corner nodes.
 *
 * This function estimates the value of a vector field at each grid node by averaging
 * the vector values from the cell centers of the cells surrounding that node (up to 8).
 * It handles physical boundaries by averaging only the available adjacent cells.
 *
 * Assumes input `centfield_arr` is from a ghosted local vector (e.g., representing ucat,
 * stored using node-indexing convention) and output `field_arr` is a ghosted local
 * vector associated with `user->fda` (DOF=3, s=2), accessed using global node indices.
 *
 * @param[in]  centfield_arr  Input: 3D array (ghosted) of vector data conceptually at cell centers,
 *                            accessed via GLOBAL indices respecting the storage convention
 *                            (e.g., `ucat[k][j][i]` uses node index `i` but represents cell `C(i,j,k)` for interior).
 * @param[out] field_arr      Output: 3D array (ghosted) where interpolated node values are stored,
 *                            accessed via GLOBAL node indices (k=0..KM, j=0..JM, i=0..IM).
 * @param[in]  user           User context containing DMDA information (da and fda).
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateFieldFromCenterToCorner_Vector( // Or original name if you prefer
    Cmpnts ***centfield_arr, /* Input: Ghosted local array (fda-based, cell data at node indices) */
    Cmpnts ***field_arr,     /* Output: local array (fda-based, true node data) */
    UserCtx *user)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    DMDALocalInfo info_nodes;
    ierr = DMDAGetLocalInfo(user->fda, &info_nodes); CHKERRQ(ierr);

    // Node ownership range (GLOBAL indices)
    PetscInt xs_node = info_nodes.xs, xm_node = info_nodes.xm, xe_node = xs_node + xm_node;
    PetscInt ys_node = info_nodes.ys, ym_node = info_nodes.ym, ye_node = ys_node + ym_node;
    PetscInt zs_node = info_nodes.zs, zm_node = info_nodes.zm, ze_node = zs_node + zm_node;

    // Global grid dimensions (NODES) - Used for global cell index check
    PetscInt MX_node = info_nodes.mx;
    PetscInt MY_node = info_nodes.my;
    PetscInt MZ_node = info_nodes.mz;
    PetscInt IM = MX_node - 1; // Max cell index i
    PetscInt JM = MY_node - 1; // Max cell index j
    PetscInt KM = MZ_node - 1; // Max cell index k


    // Valid range for accessing the INPUT ghosted array (using NODE indices)
    PetscInt gxs_node = info_nodes.gxs, gxm_node = info_nodes.gxm, gxe_node = gxs_node + gxm_node;
    PetscInt gys_node = info_nodes.gys, gym_node = info_nodes.gym, gye_node = gys_node + gym_node;
    PetscInt gzs_node = info_nodes.gzs, gzm_node = info_nodes.gzm, gze_node = gzs_node + gzm_node;

    // Log only if this function is allowed by the list set in main()
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Rank %d: Interpolating for owned NODES k=%d..%d, j=%d..%d, i=%d..%d\n",
              rank, zs_node, ze_node-1, ys_node, ye_node-1, xs_node, xe_node-1);

    LOG_ALLOW(GLOBAL,LOG_DEBUG, "[%d] Dumping DM state (user->fda) BEFORE GetArray:\n", rank);

    
    if(get_log_level() == LOG_DEBUG && is_function_allowed(__func__)){
    DMView(user->fda, PETSC_VIEWER_STDOUT_SELF);
    // Inside InterpolateFieldFromCenterToCorner_Vector, before the loops:
    
    PetscMPIInt    rank_check;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank_check);
    if (rank_check == 1) { // Only for Rank 1
      LOG_ALLOW(LOCAL,LOG_DEBUG, "[1] Attempting test read of OWNED INTERIOR centfield_arr[3][1][1]\n");
      Cmpnts test_val_owned_interior = centfield_arr[7][1][1];
      LOG_ALLOW(LOCAL,LOG_DEBUG, "[1] SUCCESS reading owned interior: x=%f\n", test_val_owned_interior.x);
      
      LOG_ALLOW(LOCAL,LOG_DEBUG, "[1] Attempting test read of OWNED BOUNDARY centfield_arr[3][0][0]\n");
      Cmpnts test_val_owned_boundary = centfield_arr[7][0][0];
      LOG_ALLOW(LOCAL,LOG_DEBUG, "[1] SUCCESS reading owned boundary: x=%f\n", test_val_owned_boundary.x);

      LOG_ALLOW(LOCAL,LOG_DEBUG, "[1] Attempting test read of GHOST centfield_arr[2][0][0]\n");
      Cmpnts test_val_ghost = centfield_arr[7][0][0]; // This is the line that likely crashes
      LOG_ALLOW(LOCAL,LOG_DEBUG, "[1] SUCCESS reading ghost: x=%f\n", test_val_ghost.x);
      
      }
    }

    ierr =  PetscBarrier(NULL);    
    
    // Proceed with the original loops...
    // Loop over the GLOBAL indices of the NODES owned by this processor (k, j, i)
    for (PetscInt k = zs_node; k < ze_node; k++) {
        for (PetscInt j = ys_node; j < ye_node; j++) {
            for (PetscInt i = xs_node; i < xe_node; i++) {
	      /*
              if(get_log_level() == LOG_DEBUG && is_function_allowed(__func__)){
	      // --- TEMPORARY TEST --- WORKS ONLY WHEN IM,JM,KM=5 !!!
	      if (rank == 1 && k == 3 && j == 0 && i == 0) {
		LOG_ALLOW(GLOBAL,LOG_DEBUG, "[1] Test read inside loop for (3,0,0): Accessing centfield_arr[2][0][0]\n");
		Cmpnts test_val_loop = centfield_arr[2][0][0]; // The crashing access
		LOG_ALLOW(GLOBAL,LOG_DEBUG, "[1] Test read inside loop SUCCESS: x=%f\n", test_val_loop.x);
	        }
	      // --- END TEMPORARY TEST ---
	      }
	      */
	      
                Cmpnts sum = {0.0, 0.0, 0.0};
                PetscInt count = 0;
                PetscBool attempted_read = PETSC_FALSE; // Flag to track if read was tried

                // Loop over the 8 potential cells surrounding node N(k,j,i)
                for (PetscInt dk_offset = -1; dk_offset <= 0; dk_offset++) {
                    for (PetscInt dj_offset = -1; dj_offset <= 0; dj_offset++) {
                        for (PetscInt di_offset = -1; di_offset <= 0; di_offset++) {

                            // Calculate the NODE index where the relevant cell's data is stored
                            PetscInt node_idx_k = k + dk_offset;
                            PetscInt node_idx_j = j + dj_offset;
                            PetscInt node_idx_i = i + di_offset;

			    
                            // Check if this NODE index corresponds to a valid GLOBAL cell index
                            PetscInt cell_idx_k = node_idx_k; // Cell index is same as node index for storage
                            PetscInt cell_idx_j = node_idx_j;
                            PetscInt cell_idx_i = node_idx_i;
			    

			    
                            if (cell_idx_i >= 0 && cell_idx_i <= IM - 1 && // Cell index check
                                cell_idx_j >= 0 && cell_idx_j <= JM - 1 &&
                                cell_idx_k >= 0 && cell_idx_k <= KM - 1)
                            {
                            
			      /*
			    // Check if this NODE index is valid within the GLOBAL node domain (0..Mx-1)
                            // This implicitly checks if the corresponding cell is valid (0..Mx-2)
                            if (node_idx_i >= 0 && node_idx_i < MX_node -1 &&
                                node_idx_j >= 0 && node_idx_j < MY_node -1 &&
                                node_idx_k >= 0 && node_idx_k < MZ_node -1)
                            {
			    */
			      
			    // Check if the NODE index is within the accessible LOCAL+GHOST range of centfield_arr
                                if (node_idx_i >= gxs_node && node_idx_i < gxe_node &&
                                    node_idx_j >= gys_node && node_idx_j < gye_node &&
                                    node_idx_k >= gzs_node && node_idx_k < gze_node)
                                {
                                    // Log attempt just before read
				  LOG_LOOP_ALLOW_EXACT(LOCAL, LOG_DEBUG,k,6,"PRE-READ: Rank %d targeting Node(k,j,i)=%d,%d,%d. Reading input centfield_arr[%d][%d][%d] (for cell C(%d,%d,%d))\n",
					      rank,
					      k,j,i,
					      node_idx_k,node_idx_j,node_idx_i,
					      cell_idx_k, cell_idx_j, cell_idx_i);

				    attempted_read = PETSC_TRUE; // Mark that we are attempting a read

                                    // ---> READ <---
                                    Cmpnts cell_val = centfield_arr[node_idx_k][node_idx_j][node_idx_i];

                                    // Log success immediately after read
                                    LOG_LOOP_ALLOW_EXACT(LOCAL, LOG_DEBUG,k,6,"POST-READ: Rank %d successful read from [%d][%d][%d] -> (%.2f, %.2f, %.2f)\n",
					      rank,
					      node_idx_k,node_idx_j,node_idx_i,
					      cell_val.x, cell_val.y, cell_val.z);

                                    sum.x += cell_val.x;
                                    sum.y += cell_val.y;
                                    sum.z += cell_val.z;
                                    count++;

                                } else {
                                     LOG_ALLOW(GLOBAL, LOG_WARNING, /* ... Ghost range warning ... */);
                                }
                            } // end global cell check
                        } // end di_offset
                    } // end dj_offset
                } // end dk_offset

		// ---> Convert GLOBAL node indices (k,j,i) to LOCAL indices for field_arr <---
		// field_arr is dimensioned nx * ny * nz (local node dimensions)
		// Global indices (k,j,i) range from (zs,ys,xs) to (ze-1, ye-1, xe-1)
		// Local indices range from 0 to (zm-1, ym-1, xm-1)
		PetscInt k_local = k - zs_node; // Offset by the starting global index
		PetscInt j_local = j - ys_node;
                PetscInt i_local = i - xs_node;
                // Calculate average and write to output node (k,j,i)
                if (count > 0) {
		  LOG_LOOP_ALLOW_EXACT(LOCAL, LOG_DEBUG,k,6,"PRE-WRITE: Rank %d targeting Node(k,j,i)=%d,%d,%d. Writing avg (count=%d)\n", rank,k,j,i, count);

                     // --- Defensive check (optional but recommended) ---
                     if (k_local < 0 || k_local >= info_nodes.zm ||
                         j_local < 0 || j_local >= info_nodes.ym ||
                         i_local < 0 || i_local >= info_nodes.xm) {
		       SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Calculated local write index out of bounds!");
                     }
                     // --- End check ---

                     // ---> Write using LOCAL indices <---
                     field_arr[k_local][j_local][i_local].x = sum.x / (PetscReal)count;
                     field_arr[k_local][j_local][i_local].y = sum.y / (PetscReal)count;
                     field_arr[k_local][j_local][i_local].z = sum.z / (PetscReal)count;

                     LOG_LOOP_ALLOW_EXACT(LOCAL, LOG_DEBUG,k,6,"POST-WRITE: Rank %d successful write to field_arr[%d][%d][%d] (local)\n", rank, k_local,j_local,i_local); // Log local indices

		     

                } else {
                     LOG_ALLOW(GLOBAL, LOG_WARNING, // Use WARNING or ERROR
                               "Rank %d: Node (i=%d,j=%d,k=%d) had count=0 surrounding valid cells! Check logic/ghosting. Writing zero.\n", rank, i, j, k);
                     field_arr[k_local][j_local][i_local] = (Cmpnts){0.0, 0.0, 0.0};
                }
                 // Add a log entry even if count=0 or no read was attempted for this node
                 if (!attempted_read && count==0) {
                     LOG_ALLOW(LOCAL, LOG_DEBUG,"NODE-COMPLETE: Rank %d Node(k,j,i)=%d,%d,%d finished loops, no valid cells/reads attempted.\n", rank, k, j, i);
                 } else if (count == 0 && attempted_read) {
                     // This case should ideally not happen if the ghost region check is correct
                      LOG_ALLOW(LOCAL, LOG_ERROR,"NODE-COMPLETE: Rank %d Node(k,j,i)=%d,%d,%d finished loops, attempted reads but count=0!\n", rank, k, j, i);
                 } else {
                     // This is the normal completion case after writing
		   LOG_LOOP_ALLOW_EXACT(LOCAL, LOG_DEBUG,k,6,"NODE-COMPLETE: Rank %d Node(k,j,i)=%d,%d,%d finished loops and write.\n", rank, k, j, i);
                 }

            } // End loop node i
        } // End loop node j
    } // End loop node k

    LOG_ALLOW_SYNC(GLOBAL, LOG_INFO, // Use INFO for completion
              "Rank %d completed interpolation function.\n", rank); // Changed message slightly
    return 0;
}

/**
 * @brief Interpolates a vector field from cell centers to corner nodes.
 *
 * This version is adapted to write directly into a ghosted local array obtained from DMDAVecGetArray(),
 * which allows using GLOBAL indices for writing to the OWNED portion of the array.
 *
 * @param[in]  centfield_arr  Input: 3D array (ghosted) of cell-centered data, accessed via GLOBAL indices.
 * @param[out] corner_arr     Output: 3D array (ghosted) where interpolated node values are stored,
 *                            also accessed via GLOBAL indices for the owned part.
 * @param[in]  user           User context containing DMDA information.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateFieldFromCenterToCorner_Vector_Petsc(
    Cmpnts ***centfield_arr, /* Input: Ghosted local array from Vec (read) */
    Cmpnts ***corner_arr,    /* Output: Ghosted local array from Vec (write) */
    UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info_nodes;
    PetscMPIInt    rank;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
    ierr = DMDAGetLocalInfo(user->fda, &info_nodes); CHKERRQ(ierr);

    // Node ownership range (GLOBAL indices)
    PetscInt xs_node = info_nodes.xs, xm_node = info_nodes.xm, xe_node = xs_node + xm_node;
    PetscInt ys_node = info_nodes.ys, ym_node = info_nodes.ym, ye_node = ys_node + ym_node;
    PetscInt zs_node = info_nodes.zs, zm_node = info_nodes.zm, ze_node = zs_node + zm_node;

    // Global grid dimensions (used for valid cell check)
    PetscInt IM = info_nodes.mx; // Total nodes in i-direction
    PetscInt JM = info_nodes.my; // Total nodes in j-direction
    PetscInt KM = info_nodes.mz; // Total nodes in k-direction

    // Ghosted array starting indices (from the same info struct)
    PetscInt gxs = info_nodes.gxs;
    PetscInt gys = info_nodes.gys;
    PetscInt gzs = info_nodes.gzs;

    // Loop over the GLOBAL indices of the NODES owned by this processor
    for (PetscInt k = zs_node; k < ze_node; k++) {
        for (PetscInt j = ys_node; j < ye_node; j++) {
            for (PetscInt i = xs_node; i < xe_node; i++) {
                Cmpnts sum = {0.0, 0.0, 0.0};
                PetscInt count = 0;
		
                // Loop over the 8 potential cells surrounding node N(k,j,i)
                for (PetscInt dk_offset = -1; dk_offset <= 0; dk_offset++) {
                    for (PetscInt dj_offset = -1; dj_offset <= 0; dj_offset++) {
                        for (PetscInt di_offset = -1; di_offset <= 0; di_offset++) {

			    // These are still GLOBAL cell indices
                            PetscInt global_cell_k = k + dk_offset;
                            PetscInt global_cell_j = j + dj_offset;
                            PetscInt global_cell_i = i + di_offset;

                            // Check if this corresponds to a valid GLOBAL cell index
                            if (global_cell_i >= 0 && global_cell_i < IM - 1 &&
                                global_cell_j >= 0 && global_cell_j < JM - 1 &&
                                global_cell_k >= 0 && global_cell_k < KM - 1)
                            {
			      
                                // Now, read from the local array using the local index
                                Cmpnts cell_val = centfield_arr[global_cell_k][global_cell_j][global_cell_i];

				LOG_LOOP_ALLOW_EXACT(LOCAL, LOG_DEBUG,k,6,"[Rank %d] successful read from [%d][%d][%d] -> (%.2f, %.2f, %.2f)\n",
						     rank,global_cell_k,global_cell_j,global_cell_i,cell_val.x, cell_val.y, cell_val.z);
				
                                sum.x += cell_val.x;
                                sum.y += cell_val.y;
                                sum.z += cell_val.z;
                                count++;
                            }
                        }
                    }
                }

                // --- MODIFICATION: Write using GLOBAL indices ---
                // The calculation of i_local, j_local, k_local is removed.
                // We write directly into the array using the global loop indices.
                if (count > 0) {
                    corner_arr[k][j][i].x = sum.x / (PetscReal)count;
                    corner_arr[k][j][i].y = sum.y / (PetscReal)count;
                    corner_arr[k][j][i].z = sum.z / (PetscReal)count;
                } else {
                    // This case should ideally not happen for a valid owned node, but as a failsafe:
                    corner_arr[k][j][i] = (Cmpnts){0.0, 0.0, 0.0};
                }

		LOG_LOOP_ALLOW_EXACT(LOCAL, LOG_DEBUG,k,6,"[Rank %d] Node(k,j,i)=%d,%d,%d finished loops and write.\n", rank, k, j, i);
            }
        }
    }
    return 0;
}

/**
 * @brief Interpolates a scalar field from cell centers to corner nodes.
 *
 * This version is adapted to write directly into a ghosted local array obtained from DMDAVecGetArray(),
 * which allows using GLOBAL indices for writing to the OWNED portion of the array.
 *
 * @param[in]  centfield_arr  Input: 3D array (ghosted) of scalar data at cell centers, accessed via GLOBAL indices.
 * @param[out] corner_arr     Output: 3D array (ghosted) where interpolated node values are stored,
 *                            also accessed via GLOBAL indices for the owned part.
 * @param[in]  user           User context containing DMDA information.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateFieldFromCenterToCorner_Scalar_Petsc(
    PetscReal ***centfield_arr, /* Input: Ghosted local array from Vec (read) */
    PetscReal ***corner_arr,    /* Output: Ghosted local array from Vec (write) */
    UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info_nodes;
    ierr = DMDAGetLocalInfo(user->fda, &info_nodes); CHKERRQ(ierr);

    // Node ownership range (GLOBAL indices)
    PetscInt xs_node = info_nodes.xs, xe_node = info_nodes.xs + info_nodes.xm;
    PetscInt ys_node = info_nodes.ys, ye_node = info_nodes.ys + info_nodes.ym;
    PetscInt zs_node = info_nodes.zs, ze_node = info_nodes.zs + info_nodes.zm;

    // Global grid dimensions (number of nodes)
    PetscInt IM = info_nodes.mx;
    PetscInt JM = info_nodes.my;
    PetscInt KM = info_nodes.mz;

    // Loop over the GLOBAL indices of the NODES owned by this processor
    for (PetscInt k = zs_node; k < ze_node; k++) {
        for (PetscInt j = ys_node; j < ye_node; j++) {
            for (PetscInt i = xs_node; i < xe_node; i++) {
                PetscReal sum = 0.0;
                PetscInt count = 0;

                // Loop over the 8 potential cells surrounding node N(k,j,i)
                for (PetscInt dk_offset = -1; dk_offset <= 0; dk_offset++) {
                    for (PetscInt dj_offset = -1; dj_offset <= 0; dj_offset++) {
                        for (PetscInt di_offset = -1; di_offset <= 0; di_offset++) {
                            // This is the global index of the cell-centered data point
                            PetscInt cell_idx_k = k + dk_offset;
                            PetscInt cell_idx_j = j + dj_offset;
                            PetscInt cell_idx_i = i + di_offset;

                            // Check if this corresponds to a valid GLOBAL cell index.
                            // Cell indices run from 0 to NumberOfNodes-2.
                            if (cell_idx_i >= 0 && cell_idx_i < IM - 1 &&
                                cell_idx_j >= 0 && cell_idx_j < JM - 1 &&
                                cell_idx_k >= 0 && cell_idx_k < KM - 1)
                            {
                                sum += centfield_arr[cell_idx_k][cell_idx_j][cell_idx_i];
                                count++;
                            }
                        }
                    }
                }

                // Calculate average and write to the output array using GLOBAL indices
                if (count > 0) {
                    corner_arr[k][j][i] = sum / (PetscReal)count;
                } else {
                    corner_arr[k][j][i] = 0.0; // Failsafe
                }
            }
        }
    }
    return 0;
}

/**
 * @brief Retrieves the scalar value at the cell (iCell, jCell, kCell).
 *
 * This function does a simple piecewise (zeroth‐order) Interpolatetion
 * by returning fieldScal[kCell][jCell][iCell], ignoring any fractional coordinates.
 *
 * @param[in]  fieldName  A string identifying the field (e.g. "temperature").
 * @param[in]  fieldScal  3D array of PetscReal (indexed as [k][j][i]).
 * @param[in]  iCell, jCell, kCell  Integral cell indices to sample.
 * @param[out] val        Pointer to a PetscReal that receives the sampled scalar.
 *
 * @return PetscErrorCode  0 on success
 */
PetscErrorCode PieceWiseLinearInterpolation_Scalar(
    const char   *fieldName,
    PetscReal  ***fieldScal,
    PetscInt      iCell,
    PetscInt      jCell,
    PetscInt      kCell,
    PetscReal    *val)
{
  PetscFunctionBegin;
  *val = fieldScal[kCell][jCell][iCell];
  
  // Optional logging
  LOG_ALLOW(LOCAL, LOG_DEBUG,
      "PieceWiseLinearInterpolation_Scalar: Field '%s' at (i=%d, j=%d, k=%d) => val=%.6f\n",
      fieldName, iCell, jCell, kCell, *val);

  PetscFunctionReturn(0);
}


/**
 * @brief Retrieves the vector (Cmpnts) at the cell (iCell, jCell, kCell).
 *
 * This function simply sets:
 *   vec->x = fieldVec[kCell][jCell][iCell].x
 *   vec->y = fieldVec[kCell][jCell][iCell].y
 *   vec->z = fieldVec[kCell][jCell][iCell].z
 * effectively a nearest‐cell or piecewise approach.
 *
 * @param[in]  fieldName  String identifying the field (e.g. "velocity").
 * @param[in]  fieldVec   3D array of Cmpnts (indexed as [k][j][i]).
 * @param[in]  iCell, jCell, kCell  Integral cell indices to sample.
 * @param[out] vec        Pointer to a Cmpnts struct that receives the sampled vector.
 *
 * @return PetscErrorCode  0 on success
 */
PetscErrorCode PieceWiseLinearInterpolation_Vector(
    const char   *fieldName,
    Cmpnts     ***fieldVec,
    PetscInt      iCell,
    PetscInt      jCell,
    PetscInt      kCell,
    Cmpnts       *vec)
{
  PetscFunctionBegin;
  vec->x = fieldVec[kCell][jCell][iCell].x;
  vec->y = fieldVec[kCell][jCell][iCell].y;
  vec->z = fieldVec[kCell][jCell][iCell].z;

  // Optional logging
  LOG_ALLOW(LOCAL, LOG_DEBUG,
      "PieceWiseLinearInterpolation_Vector: Field '%s' at (i=%d, j=%d, k=%d) => (x=%.6f, y=%.6f, z=%.6f)\n",
      fieldName, iCell, jCell, kCell, vec->x, vec->y, vec->z);

  PetscFunctionReturn(0);
}

/**
 * @brief Computes the trilinear Interpolatetion weights from the Interpolatetion coefficients.
 *
 * This function computes the weights for trilinear Interpolatetion at the eight corners of a cell
 * using the Interpolatetion coefficients provided along the x, y, and z directions.
 *
 * @param[in]  a1 Interpolation coefficient along the x-direction (normalized coordinate within the cell).
 * @param[in]  a2 Interpolation coefficient along the y-direction (normalized coordinate within the cell).
 * @param[in]  a3 Interpolation coefficient along the z-direction (normalized coordinate within the cell).
 * @param[out] w  Array of 8 weights, each corresponding to one corner of the cell.
 *
 * @note
 * - The coefficients `a1`, `a2`, and `a3` should be in the range [0, 1].
 * - The order of weights corresponds to the eight corners of a hexahedral cell.
 */
static inline void ComputeTrilinearWeights(PetscReal a1, PetscReal a2, PetscReal a3, PetscReal *w) {
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Computing weights for a1=%f, a2=%f, a3=%f.\n", a1, a2, a3);

    // Ensure a1, a2, a3 are within [0,1]
    a1 = PetscMax(0.0, PetscMin(1.0, a1));
    a2 = PetscMax(0.0, PetscMin(1.0, a2));
    a3 = PetscMax(0.0, PetscMin(1.0, a3));

    const PetscReal oa1 = 1.0 - a1;
    const PetscReal oa2 = 1.0 - a2;
    const PetscReal oa3 = 1.0 - a3;

    w[0] = oa1 * oa2 * oa3;  /* cornerOffsets[0] => (0,0,0) */
    w[1] = a1  * oa2 * oa3;  /* cornerOffsets[1] => (1,0,0) */
    w[2] = oa1 * a2  * oa3;  /* cornerOffsets[2] => (0,1,0) */
    w[3] = a1  * a2  * oa3;  /* cornerOffsets[3] => (1,1,0) */
    w[4] = oa1 * oa2 * a3;   /* cornerOffsets[4] => (0,0,1) */
    w[5] = a1  * oa2 * a3;   /* cornerOffsets[5] => (1,0,1) */
    w[6] = oa1 * a2  * a3;   /* cornerOffsets[6] => (0,1,1) */
    w[7] = a1  * a2  * a3;   /* cornerOffsets[7] => (1,1,1) */

    // Log the computed weights for debugging
    LOG_ALLOW(LOCAL,LOG_DEBUG, "Weights computed - "
        "w0=%f, w1=%f, w2=%f, w3=%f, w4=%f, w5=%f, w6=%f, w7=%f. \n",
        w[0], w[1], w[2], w[3], w[4], w[5], w[6], w[7]);
}

/**
 * @brief Computes the trilinear Interpolateted scalar at a given point.
 *
 * @param[in]  fieldName A string representing the name of the scalar field (e.g., "temperature").
 * @param[in]  fieldScal 3D array of the field from a DMDA (indexed as [k][j][i]),
 *                      each cell a PetscReal.
 * @param[in]  i, j, k   Integral cell indices (the "lower" corner in each dimension).
 * @param[in]  a1, a2, a3 Normalized coordinates within the cell ([0,1] range).
 * @param[out] val       Pointer to a PetscReal that will store the Interpolateted scalar.
 *
 * This function uses the standard 8-corner trilinear formula via `ComputeTrilinearWeights()`.
 * If a different scheme is desired, implement a new function with the same PetscInterface.
 */
PetscErrorCode TrilinearInterpolation_Scalar(
    const char   *fieldName,
    PetscReal  ***fieldScal,
    PetscInt      i,
    PetscInt      j,
    PetscInt      k,
    PetscReal     a1,
    PetscReal     a2,
    PetscReal     a3,
    PetscReal    *val)
{
    PetscFunctionBegin; // PETSc macro for error/stack tracing

    // Compute the 8 corner weights
    PetscReal wcorner[8];
    ComputeTrilinearWeights(a1, a2, a3, wcorner);

    // Offsets for cell corners
    PetscInt i1 = i + 1;
    PetscInt j1 = j + 1;
    PetscInt k1 = k + 1;

    // Initialize the output scalar
    PetscReal sum = 0.0;

    // Corner 0 => (i, j, k)
    sum += wcorner[0] * fieldScal[k ][j ][i ];
    // Corner 1 => (i+1, j, k)
    sum += wcorner[1] * fieldScal[k ][j ][i1];
    // Corner 2 => (i, j+1, k)
    sum += wcorner[2] * fieldScal[k ][j1][i ];
    // Corner 3 => (i+1, j+1, k)
    sum += wcorner[3] * fieldScal[k ][j1][i1];
    // Corner 4 => (i, j, k+1)
    sum += wcorner[4] * fieldScal[k1][j ][i ];
    // Corner 5 => (i+1, j, k+1)
    sum += wcorner[5] * fieldScal[k1][j ][i1];
    // Corner 6 => (i, j+1, k+1)
    sum += wcorner[6] * fieldScal[k1][j1][i ];
    // Corner 7 => (i+1, j+1, k+1)
    sum += wcorner[7] * fieldScal[k1][j1][i1];

    *val = sum;

    // Logging (optional)
    LOG_ALLOW(LOCAL, LOG_DEBUG,
        "TrilinearInterpolation_Scalar: Field '%s' at (i=%d, j=%d, k=%d), "
        "a1=%.6f, a2=%.6f, a3=%.6f -> val=%.6f.\n",
        fieldName, i, j, k, a1, a2, a3, *val);

    // LOG_ALLOW_SYNC(GLOBAL, LOG_INFO,
    //    "TrilinearInterpolation_Scalar: Completed Interpolatetion for field '%s' across local cells.\n",
	  //    fieldName);

    PetscFunctionReturn(0);
}


/**
 * @brief Computes the trilinear Interpolateted vector (e.g., velocity) at a given point.
 *
 * @param[in]  fieldName  A string representing the name of the vector field (e.g., "velocity").
 * @param[in]  fieldVec   3D array of the field from a DMDA (indexed as [k][j][i]),
 *                        each cell of type Cmpnts.
 * @param[in]  i, j, k    Integral cell indices (the "lower" corner in each dimension).
 * @param[in]  a1, a2, a3 Normalized coordinates within the cell ([0,1] range).
 * @param[out] vec        Pointer to a Cmpnts struct that will store the Interpolateted vector (x, y, z).
 *
 * This function uses the standard 8-corner trilinear formula via `ComputeTrilinearWeights()`.
 * If a different scheme is desired, implement a new function with the same PetscInterface.
 */
/**
 * @brief Computes the trilinear Interpolateted vector at a given point, with partial weighting if corners go out of range.
 *
 * If any of the 8 corners (i or i+1, j or j+1, k or k+1) is out of [0..mx), [0..my), [0..mz),
 * that corner is skipped and the corresponding weight is omitted. The total is normalized
 * by the sum of used weights, so we get a partial Interpolatetion near boundaries.
 */
PetscErrorCode TrilinearInterpolation_Vector(
    const char   *fieldName,
    Cmpnts     ***fieldVec,  /* 3D array [k][j][i], dimension [mz][my][mx] */
    PetscInt      i, // local cell index i
    PetscInt      j, // local cell index j
    PetscInt      k, // local cell index k
    PetscReal     a1,
    PetscReal     a2,
    PetscReal     a3,
    Cmpnts       *vec)
{
  PetscFunctionBegin; // PETSc macro for error/stack tracing

  // Compute the 8 corner weights
  PetscErrorCode ierr;
  PetscReal wcorner[8];
  PetscMPIInt rank;

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  
  LOG_ALLOW(LOCAL,LOG_DEBUG,"[Rank %d] Computing Trilinear Weights.\n",rank);
  ComputeTrilinearWeights(a1, a2, a3, wcorner);
  
  LOG_ALLOW(LOCAL,LOG_DEBUG,"[Rank %d]  Trilinear weights computed for local cell %d,%d,%d.\n",rank,i,j,k);

  // For partial Interpolatetion, we'll keep track of how many corners are valid
  // and how much sum of weights is used. Then we do a final normalization.
  PetscReal sumW = 0.0;
  Cmpnts    accum = {0.0, 0.0, 0.0};

  // The eight corner indices, with their weights:
  // corners: (i,j,k), (i+1,j,k), (i,j+1,k), (i+1,j+1,k), etc.
  // We store them in an array to iterate cleanly.
  const PetscInt cornerOffsets[8][3] = {
    {0, 0, 0},
    {1, 0, 0},
    {0, 1, 0},
    {1, 1, 0},
    {0, 0, 1},
    {1, 0, 1},
    {0, 1, 1},
    {1, 1, 1}
  };

  // Weighted partial sum
  for (PetscInt c = 0; c < 8; c++) {
    const PetscInt di = cornerOffsets[c][0];
    const PetscInt dj = cornerOffsets[c][1];
    const PetscInt dk = cornerOffsets[c][2];
    PetscInt iC = i + di;
    PetscInt jC = j + dj;
    PetscInt kC = k + dk;

    /*
    // skip if out of domain
    // (Assuming you know global domain is [0..mx), [0..my), [0..mz).)
    if (iC < 0 || iC >= (PetscInt)userGlobalMx ||
        jC < 0 || jC >= (PetscInt)userGlobalMy ||
        kC < 0 || kC >= (PetscInt)userGlobalMz)
    {
      // skip this corner
      continue;
    }

    */

    LOG_ALLOW(LOCAL,LOG_DEBUG,"[Rank %d] %s[%d][%d][%d] = (%.4f,%.4f,%.4f).\n",rank,fieldName,kC,jC,iC,fieldVec[kC][jC][iC].x,fieldVec[kC][jC][iC].y,fieldVec[kC][jC][iC].z);
    
    // Otherwise, accumulate
    accum.x += wcorner[c] * fieldVec[kC][jC][iC].x;
    accum.y += wcorner[c] * fieldVec[kC][jC][iC].y;
    accum.z += wcorner[c] * fieldVec[kC][jC][iC].z;
    sumW += wcorner[c];
  }

  // If sumW=0 => out-of-range or zero weighting => set (0,0,0)
  if (sumW > 1.0e-14) {
    vec->x = accum.x / sumW;
    vec->y = accum.y / sumW;
    vec->z = accum.z / sumW;
  } else {
    vec->x = 0.0;  vec->y = 0.0;  vec->z = 0.0;
  }

  PetscFunctionReturn(0);
}


/**
 * @brief Interpolates a single field (scalar or vector) for one particle, storing the result in a swarm array.
 *
 * This helper function is used inside a larger loop over local particles.
 *
 * @param[in]  fieldName   A string identifying the field (e.g. "velocity", "temperature", etc.).
 * @param[in]  fieldPtr    A Pointer to the local DMDA array for the field:
 *                         - (PetscReal ***) if scalar (blockSize = 1)
 *                         - (Cmpnts    ***) if vector (blockSize = 3)
 * @param[in]  iCell, jCell, kCell  The cell indices for this particle (already clamped if needed).
 * @param[in]  a1, a2, a3   Interpolation coefficients in [0..1].
 *                          If using PiecewiseLinearInterpolation, these are currently ignored.
 * @param[out] swarmOut     A Pointer to the DMSwarm output array:
 *                            - (PetscReal*) if scalar dof=1
 *                            - (PetscReal*) if vector dof=3 (storing x,y,z in consecutive reals)
 *                            - or (Cmpnts*) if you store the result as a Cmpnts struct
 * @param[in]  p            The local particle index (used to compute the correct offset in swarmOut).
 * @param[in]  blockSize    The number of degrees of freedom (1 => scalar, 3 => vector).
 *
 * This routine demonstrates a switch between:
 *  - PiecewiseLinearInterpolation (zeroth order / nearest cell)
 *  - TrilinearInterpolation (8-corner weighted).
 *
 * By default, PiecewiseLinearInterpolation is active, while the TrilinearInterpolation calls are commented out.
 * To switch to trilinear, simply comment/uncomment appropriately.
 *
 * Logging calls are provided (LOG_ALLOW) that you can adapt to your existing logging system.
 *
 * @return PetscErrorCode 0 on success, non-zero on error.
 */
static inline PetscErrorCode InterpolateEulerFieldToSwarmForParticle(
    const char  *fieldName,
    void        *fieldPtr,   /* typed Pointer => either (PetscReal***) or (Cmpnts***) */
    PetscInt     iCell,
    PetscInt     jCell,
    PetscInt     kCell,
    PetscReal    a1,
    PetscReal    a2,
    PetscReal    a3,
    void        *swarmOut,   /* typed Pointer => (PetscReal*) or (Cmpnts*) or dof=3 array */
    PetscInt     p,          /* particle index */
    PetscInt     blockSize)  /* dof=1 => scalar, dof=3 => vector */
{
  PetscErrorCode ierr;
  PetscFunctionBegin;

  // Optional logging at start
  LOG_ALLOW(LOCAL, LOG_DEBUG,
    "InterpolateEulerFieldToSwarmForParticle: field='%s', blockSize=%d, "
    "cell IDs=(%d,%d,%d), weights=(%.4f,%.4f,%.4f)\n",
    fieldName, blockSize, iCell, jCell, kCell, a1, a2, a3);

  /*
     If blockSize=1, we PetscInterpret the fieldPtr as a 3D array of PetscReal (scalar).
     If blockSize=3, we PetscInterpret the fieldPtr as a 3D array of Cmpnts (vector).
   */
  if (blockSize == 1) {
    /* Scalar field: Cast fieldPtr to (PetscReal ***). */
    PetscReal ***fieldScal = (PetscReal ***) fieldPtr;
    PetscReal val;

    // Currently using trilinear.
    ierr = TrilinearInterpolation(fieldName, fieldScal,
                                   iCell, jCell, kCell,
                                   a1, a2, a3,
                                   &val);
     CHKERRQ(ierr);

    // Alternative (commented) call to PiecewiseLinearInterpolation (zeroth order) :
    //   PetscErrorCode ierr = PiecewiseLinearInterpolation(fieldName,
    //                                                   fieldScal,
    //                                                   iCell, jCell, kCell,
    //                                                   &val);
    // CHKERRQ(ierr);

    // Write the scalar result to the swarm output at index [p].
    ((PetscReal*)swarmOut)[p] = val;

    LOG_ALLOW(LOCAL, LOG_DEBUG,
      "InterpolateEulerFieldToSwarmForParticle [Scalar]: field='%s', result=%.6f "
      "stored at swarmOut index p=%d.\n", fieldName, val, (PetscInt)p);
  }
  else if (blockSize == 3) {
    /* Vector field: Cast fieldPtr to (Cmpnts ***). */
    Cmpnts ***fieldVec = (Cmpnts ***) fieldPtr;
    Cmpnts vec;


    
    
    // Piecewise Interpolatetion (zeroth order).
    //  PetscErrorCode ierr = PieceWiseLinearInterpolation(fieldName,
    //                                                   fieldVec,
    //                                                   iCell, jCell, kCell,
    //                                                   &vec);
    // CHKERRQ(ierr);

    // Alternative (commented) call to trilinear:
     ierr = TrilinearInterpolation(fieldName, fieldVec,
                                   iCell, jCell, kCell,
                                   a1, a2, a3,
                                   &vec);
     CHKERRQ(ierr);

    // If swarmOut is an array of 3 reals per particle:
    ((PetscReal*)swarmOut)[3*p + 0] = vec.x;
    ((PetscReal*)swarmOut)[3*p + 1] = vec.y;
    ((PetscReal*)swarmOut)[3*p + 2] = vec.z;

    LOG_ALLOW(LOCAL, LOG_DEBUG,
      "InterpolateEulerFieldToSwarmForParticle [Vector]: field='%s', result=(%.6f,%.6f,%.6f) "
      "stored at swarmOut[3p..3p+2], p=%d.\n",
      fieldName, vec.x, vec.y, vec.z, (PetscInt)p);

    /*
       If you store the vector result as a Cmpnts in the swarm, do instead:
         ((Cmpnts*)swarmOut)[p] = vec;
       but ensure your DMSwarm field is sized for a Cmpnts struct.
    */
  }
  else {
    /* If blockSize isn't 1 or 3, we raise an error. */
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
      "InterpolateEulerFieldToSwarmForParticle: only blockSize=1 or 3 supported, got %d.",
      (PetscInt)blockSize);
  }

  PetscFunctionReturn(0);
}

/**
 * @brief Interpolates a cell-centered field (scalar or vector) onto DMSwarm particles,
 *        using a robust, PETSc-idiomatic two-stage process.
 *
 * This function first converts the cell-centered input data to corner-node data,
 * storing this intermediate result in a PETSc Vec to correctly handle the communication
 * of ghost-point information across parallel ranks. It then performs a final trilinear
 * interpolation from the ghosted corner data to each particle's location.
 *
 * Workflow:
 *   1.  Create temporary PETSc Vecs (`cornerGlobal`, `cornerLocal`) to manage the
 *       intermediate corner-node data, using the existing nodal DMDA (`user->fda`).
 *   2.  Call a dispatch macro that uses the runtime block size (`bs`) to select the
 *       correct underlying center-to-corner function, writing results into `cornerGlobal`.
 *   3.  Perform a ghost-point exchange (`DMGlobalToLocal`) to transfer the boundary
 *       data from `cornerGlobal` into the ghost regions of `cornerLocal`.
 *   4.  Loop over all local particles. For each particle:
 *       a. Convert its global cell index to a local index relative to the ghosted array.
 *       b. Check if the particle's interpolation stencil is fully contained within
 *          the owned+ghost region. If not, log a warning and set the result to zero.
 *       c. Perform the final trilinear interpolation using the ghosted `cornerLocal` data.
 *   5.  Restore all PETSc objects to prevent memory leaks.
 *
 * @param[in]  user                     User context with DMDA, DMSwarm, etc.
 * @param[in]  fieldGlobal_cellCentered Vec (from fda) with cell-centered data (e.g., Ucat).
 * @param[in]  fieldName                Human-readable field name for logging (e.g., "Ucat").
 * @param[in]  swarmOutFieldName        Name of the DMSwarm field where interpolation results go.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateEulerFieldToSwarm(
    UserCtx    *user,
    Vec         fieldLocal_cellCentered,
    const char *fieldName,
    const char *swarmOutFieldName)
{
    PetscErrorCode ierr;
    DM             fda = user->fda;
    DM             swarm = user->swarm;
    PetscInt       bs;
    DMDALocalInfo  info;
    PetscMPIInt    rank;

    // Intermediate vectors to manage corner data and its ghosts
    Vec cornerGlobal, cornerLocal;

    // Generic pointers to the raw data arrays
    void *cellCenterPtr_read;
    void *cornerPtr_read_with_ghosts;

    // Swarm-related pointers
    PetscInt  *cellIDs = NULL;
    PetscReal *weights = NULL;
    PetscInt  *pids = NULL; // For logging particle IDs in warnings
    void      *swarmOut = NULL;
    PetscInt   nLocal;

    PetscFunctionBegin;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(fda, &info); CHKERRQ(ierr);
    ierr = VecGetBlockSize(fieldLocal_cellCentered, &bs); CHKERRQ(ierr);
    if (bs != 1 && bs != 3) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP, "BlockSize must be 1 or 3.");

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting for field '%s'... \n", fieldName);

    /* STAGE 1: Center-to-Corner Calculation with Communication */

    // (A) Get temporary Vecs from the DMDA pool to manage corner data
    ierr = DMGetGlobalVector(fda, &cornerGlobal); CHKERRQ(ierr);
    ierr = DMGetLocalVector(fda, &cornerLocal); CHKERRQ(ierr);

    ierr = VecSet(cornerGlobal,0.0); CHKERRQ(ierr);
    ierr = VecSet(cornerLocal,0.0); CHKERRQ(ierr);

    // (B) Get a read-only array from the input cell-centered vector
    ierr = DMDAVecGetArrayRead(fda, fieldLocal_cellCentered, &cellCenterPtr_read); CHKERRQ(ierr);

    // (C) Perform the center-to-corner interpolation directly into the global corner vector
    void *cornerPtr_write = NULL;
    ierr = DMDAVecGetArray(fda, cornerLocal, &cornerPtr_write); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d] Starting center-to-corner interpolation for '%s'.\n", rank, fieldName);
        
    // SINGLE, CLEAN CALL SITE: The macro handles the runtime dispatch based on 'bs'.
    ierr = InterpolateFieldFromCenterToCorner(bs, cellCenterPtr_read, cornerPtr_write, user); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d] Finished center-to-corner interpolation for '%s'.\n", rank, fieldName);
    ierr = DMDAVecRestoreArray(fda, cornerLocal, &cornerPtr_write); CHKERRQ(ierr);
    
    ierr = DMDAVecRestoreArrayRead(fda, fieldLocal_cellCentered, &cellCenterPtr_read); CHKERRQ(ierr);

    //////////// DEBUG
    /*
        // Log a synchronized header message from all ranks.
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "DEBUG: Dumping cornerGlobal BEFORE ghost exchange...\n");

    // VecView prints to a PETSc viewer. PETSC_VIEWER_STDOUT_WORLD is a global, synchronized viewer.
    // We only need to call it if logging is active for this function.
    if (is_function_allowed(__func__) && (int)(LOG_DEBUG) <= (int)get_log_level()) {
        ierr = VecView(cornerGlobal, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    }
    
    // Log a synchronized footer message.
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "DEBUG: Finished dumping cornerGlobal.\n");
    */
    ////////// DEBUG

    // ierr = PetscBarrier((PetscObject)cornerGlobal); CHKERRQ(ierr);
    
    // (D) CRITICAL STEP: Communicate the newly computed corner data to fill ghost regions
    LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d] Beginning ghost exchange for corner data...\n", rank);
    ierr = DMLocalToGlobalBegin(fda, cornerLocal, INSERT_VALUES, cornerGlobal); CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(fda, cornerLocal, INSERT_VALUES, cornerGlobal); CHKERRQ(ierr);

    ierr = DMGlobalToLocalBegin(fda, cornerGlobal, INSERT_VALUES, cornerLocal); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(fda, cornerGlobal, INSERT_VALUES, cornerLocal); CHKERRQ(ierr);
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d] Ghost exchange for corner data complete.\n", rank);

    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "DEBUG: About to get array from FDA on all ranks.\n");
       if (is_function_allowed(__func__) && (int)(LOG_DEBUG) <= (int)get_log_level()) {
	 ierr = DMView(user->fda, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
       }
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "DEBUG: Finished DMView.\n");
    
    /////// DEBUG
    /*
        // Log a synchronized header message from all ranks.
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "DEBUG: Dumping cornerLocal AFTER ghost exchange...\n");

    // Here, we want each rank to print its own local vector.
    // PETSC_VIEWER_STDOUT_SELF is the correct tool for this. The LOG_ALLOW_SYNC
    // wrapper will ensure the output from each rank is printed sequentially.
    if (is_function_allowed(__func__) && (int)(LOG_DEBUG) <= (int)get_log_level()) {
        PetscMPIInt rank_d;
        ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank_d); CHKERRQ(ierr);
        
        // Use a synchronized print to create a clean header for each rank's output
        PetscSynchronizedPrintf(PETSC_COMM_WORLD, "\n--- cornerLocal (Rank %d) ---\n", rank_d);
        PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);

        // Print the local vector's contents for this rank
        ierr = VecView(cornerLocal, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);

        // Use a synchronized print to create a clean footer
        PetscSynchronizedPrintf(PETSC_COMM_WORLD, "--- End cornerLocal (Rank %d) ---\n", rank_d);
        PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    }

    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG, "DEBUG: Finished dumping cornerLocal.\n");
    */
    ///// DEBUG
    
    /* STAGE 2: Particle Interpolation using Ghosted Corner Data */

    // (E) Get the local, GHOSTED array of corner data for the final interpolation step
    ierr = DMDAVecGetArrayRead(fda, cornerLocal, &cornerPtr_read_with_ghosts); CHKERRQ(ierr);

    // (F) Retrieve swarm fields for the particle loop
    ierr = DMSwarmGetLocalSize(swarm, &nLocal); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&pids); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, swarmOutFieldName, NULL, NULL, &swarmOut); CHKERRQ(ierr);

    // (G) Loop over each local particle
    for (PetscInt p = 0; p < nLocal; p++) {
        PetscInt iCell_global = cellIDs[3 * p + 0];
        PetscInt jCell_global = cellIDs[3 * p + 1];
        PetscInt kCell_global = cellIDs[3 * p + 2];

        // Convert GLOBAL cell index to LOCAL GHOSTED index for the corner array
        PetscInt iCell_local = iCell_global;//- info.gxs;
        PetscInt jCell_local = jCell_global;//- info.gys;
        PetscInt kCell_local = kCell_global;//- info.gzs;

        // Safety check: Ensure the entire 8-node stencil is within the valid memory
        // region (owned + ghosts) of the local array.

	LOG_ALLOW(LOCAL,LOG_DEBUG," The ghosted max of current rank %d are: %d,%d,%d.\n",rank,info.gxs + info.gxm,info.gys+info.gym,info.gzs+info.gzm);
	
        if (iCell_local < 0 || iCell_local >= info.gxs + info.gxm - 1 ||
            jCell_local < 0 || jCell_local >= info.gys + info.gym - 1 ||
            kCell_local < 0 || kCell_local >= info.gzs + info.gzm - 1)
        {
            LOG_ALLOW(LOCAL, LOG_WARNING,
                      "[Rank %d] Particle PID %lld in global cell (%d,%d,%d) is in an un-interpolatable region (requires ghosts of ghosts or is out of bounds). Zeroing field '%s'.\n",
                      rank, (long long)pids[p], iCell_global, jCell_global, kCell_global, fieldName);

            if (bs == 3) {
                ((PetscReal*)swarmOut)[3*p + 0] = 0.0;
                ((PetscReal*)swarmOut)[3*p + 1] = 0.0;
                ((PetscReal*)swarmOut)[3*p + 2] = 0.0;
            } else {
                 ((PetscReal*)swarmOut)[p] = 0.0;
            }
            continue;
        }

        PetscReal a1 = weights[3 * p + 0];
        PetscReal a2 = weights[3 * p + 1];
        PetscReal a3 = weights[3 * p + 2];

        ierr = InterpolateEulerFieldToSwarmForParticle(
                  fieldName,
                  cornerPtr_read_with_ghosts,
                  iCell_local, jCell_local, kCell_local,
                  a1, a2, a3,
                  swarmOut, p, bs);
        CHKERRQ(ierr);
    }

    /* (H) Restore all retrieved arrays and vectors */
    ierr = DMDAVecRestoreArrayRead(fda, cornerLocal, &cornerPtr_read_with_ghosts); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, swarmOutFieldName, NULL, NULL, &swarmOut); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "weight", NULL, NULL, (void**)&weights); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_pid", NULL, NULL, (void**)&pids); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);

    ierr = DMRestoreLocalVector(fda, &cornerLocal); CHKERRQ(ierr);
    ierr = DMRestoreGlobalVector(fda, &cornerGlobal); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Finished for field '%s'.\n", fieldName);
    PetscFunctionReturn(0);
}

/**
 * @brief Interpolates a cell-centered field (scalar or vector) onto DMSwarm particles,
 *        converting the cell-center data to corner data first, then looping over particles.
 *
 * Steps:
 *   1) Check that the Vec has blockSize=1 or 3 (scalar vs. vector).
 *   2) Map the cell-centered Vec to a local array (fieldGlobal -> localPtr).
 *   3) Allocate a corner array (cornerPtr) via Allocate3DArray(...), sized (zm+1, ym+1, xm+1).
 *   4) Convert from cell-centers to corners via InterpolateFieldFromCenterToCorner(...).
 *   5) Restore the cell-centered array.
 *   6) Retrieve DMSwarm fields: "DMSwarm_CellID", "weight", and swarmOutFieldName.
 *   7) Loop over local particles, clamp i/j/k, skip or zero out if out of range, read (a1,a2,a3).
 *   8) Call InterpolateEulerFieldToSwarmForParticle(...) with cornerPtr to do final Interpolatetion.
 *   9) Restore swarm fields, free the corner array.
 *
 * @param[in]  user              User context with:
 *                                - user->da     (cell-centered DMDA),
 *                                - user->swarm  (DMSwarm).
 * @param[in]  fieldGlobal       Vec with blockSize=1 or 3, storing the cell-centered field.
 * @param[in]  fieldName         Human-readable field name for logging (e.g. "velocity").
 * @param[in]  swarmOutFieldName Name of the DMSwarm field where Interpolatetion results go.
 *
 * @return PetscErrorCode  0 on success, non-zero on error.
 */
/*
PetscErrorCode InterpolateEulerFieldToSwarm(
    UserCtx    *user,
    Vec         fieldGlobal,       //DMDA Vec containing cell‐center data 
    const char *fieldName,         // e.g., "Ucat" 
    const char *swarmOutFieldName) // Name of the output DMSwarm field 
{
  PetscErrorCode ierr;
  DM             fda    = user->fda;      // DM for cell‐center field data 
  //  DM             da     = user->da;        DM for grid information (local indices) 
  DM             swarm  = user->swarm;     DMSwarm for particles 
  PetscInt       bs;                      Block size: 1 (scalar) or 3 (vector) 
  DMDALocalInfo  info;                  // Local grid info 
void          *localPtr   = NULL;       // Pointer to cell‐center data from fda 
void          *cornerPtr  = NULL;       // Will hold the typed corner array 
  void          *swarmOut   = NULL;     // Pointer to the swarm output field 
  PetscInt    *cellIDs    = NULL;     // Particle cell indices from swarm 
  PetscReal     *weights    = NULL;     // Interpolation coefficients from swarm 
  PetscInt       nLocal;
  PetscMPIInt    rank;

  
  
  PetscFunctionBegin;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
  // (A) Check block size and get local domain info 
  ierr = VecGetBlockSize(fieldGlobal, &bs); CHKERRQ(ierr);
  if (bs != 1 && bs != 3) {
    SETERRQ(PETSC_COMM_SELF, PETSC_ERR_SUP,
             "InterpolateEulerFieldToSwarm: blockSize must be 1 or 3, got %d.", (PetscInt)bs);
  }
  ierr = DMDAGetLocalInfo(fda, &info); CHKERRQ(ierr);

  PetscInt xs_node = info.xs;
  PetscInt ys_node = info.ys;
  PetscInt zs_node = info.zs;
  
  LOG_ALLOW(GLOBAL, LOG_INFO,
    "InterpolateEulerFieldToSwarm: Starting with field='%s', blockSize=%d, local domain: (%d x %d x %d)\n",
    fieldName, bs, info.xm, info.ym, info.zm);

  // (B) Map the cell-centered Vec to a local array using the DM attached to fieldGlobal 
  LOG_ALLOW(LOCAL,LOG_DEBUG, "[Rank %d] Dumping DM state (user->fda) BEFORE GetArray:\n", rank);

  if(get_log_level() == LOG_DEBUG && is_function_allowed(__func__)) DMView(user->fda, PETSC_VIEWER_STDOUT_SELF);

  ierr = DMDAVecGetArrayRead(fda, fieldGlobal, &localPtr); CHKERRQ(ierr);

    // Corner domain is larger than cell center domain 
    PetscInt nz = info.zm; 
    PetscInt ny = info.ym;
    PetscInt nx = info.xm;
    LOG_ALLOW(LOCAL, LOG_DEBUG,
      "[Rank %d]  Allocating corner array of size (%d x %d x %d), blockSize = %d\n",
	      rank,nx, ny, nz, bs);
    if (bs == 1) {
      // Declare a typed Pointer for scalar corners 
      PetscReal ***cornerScal = NULL;
      ierr = Allocate3DArray(&cornerScal, nz, ny, nx); CHKERRQ(ierr);
      if (!cornerScal) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "InterpolateEulerFieldToSwarm: 3D Array Allocation for scalar corners failed.\n");
      }
      // Save typed Pointer into cornerPtr so later code can cast appropriately 
      cornerPtr = (void*) cornerScal;
      // (D) Convert cell-center data to corners for scalar field 
      ierr = InterpolateFieldFromCenterToCorner( (PetscReal ***) localPtr, cornerScal, user); CHKERRQ(ierr);
    } else {
      // For vector fields 
      Cmpnts ***cornerVec = NULL;
      ierr = Allocate3DArray(&cornerVec, nz, ny, nx); CHKERRQ(ierr);
      if (!cornerVec) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "InterpolateEulerFieldToSwarm: 3D Array Allocation for vector corners failed.\n");
      }
      cornerPtr = (void*) cornerVec;

      // Comment out the actual call:
      ierr = InterpolateFieldFromCenterToCorner( (Cmpnts ***) localPtr, (Cmpnts ***)cornerPtr, user); CHKERRQ(ierr);
      */
      /*
      // --- DEBUG - BYPASS SECTION ---
      LOG_ALLOW(LOCAL,LOG_DEBUG, "[%d] DEBUG: Bypassing InterpolateFieldFromCenterToCorner call.\n", rank);


      // Instead, populate cornerPtr with fixed values
      // Use local indices 0..nz-1, 0..ny-1, 0..nx-1 because cornerPtr is a local array
      if (bs == 3) {
	Cmpnts ***cornerVec_bypass = (Cmpnts ***)cornerPtr;
	for (PetscInt k_local = 0; k_local < nz; ++k_local) {
          for (PetscInt j_local = 0; j_local < ny; ++j_local) {
	    for (PetscInt i_local = 0; i_local < nx; ++i_local) {
	      cornerVec_bypass[k_local][j_local][i_local].x = 1.0; // Or 0.0, or rank number
	      cornerVec_bypass[k_local][j_local][i_local].y = 2.0;
	      cornerVec_bypass[k_local][j_local][i_local].z = 3.0;
	    }
          }
	}
	LOG_ALLOW(LOCAL,LOG_DEBUG, "[%d] DEBUG: Finished setting fixed values in cornerVec_bypass.\n", rank);
      
      }
      */
      /*
    LOG_ALLOW(LOCAL, LOG_INFO,
      "InterpolateEulerFieldToSwarm: Rank %d Completed center-to-corner Interpolatetion for field='%s'.\n",
	      rank,fieldName);
    }
  // (E) Restore the cell-centered array since we now have corner data 
  ierr = DMDAVecRestoreArrayRead(fda, fieldGlobal, &localPtr); CHKERRQ(ierr);

  // (F) Retrieve swarm fields: cell IDs, weights, and the output field 
  ierr = DMSwarmGetLocalSize(swarm, &nLocal); CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "DMSwarm_CellID",  NULL, NULL, (void**)&cellIDs);  CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, "weight",          NULL, NULL, (void**)&weights);  CHKERRQ(ierr);
  ierr = DMSwarmGetField(swarm, swarmOutFieldName, NULL, NULL, &swarmOut);         CHKERRQ(ierr);
  LOG_ALLOW(GLOBAL, LOG_INFO,
    "InterpolateEulerFieldToSwarm: Found %d local particles to process for field='%s'.\n",
    nLocal, fieldName);

  // (G) Loop over each local particle and perform the final Interpolatetion from corners 
  for (PetscInt p = 0; p < nLocal; p++) {
    PetscInt iCell_global = (PetscInt)cellIDs[3*p + 0];
    PetscInt jCell_global = (PetscInt)cellIDs[3*p + 1];
    PetscInt kCell_global = (PetscInt)cellIDs[3*p + 2];

// --- Convert GLOBAL cell index to LOCAL cell index ---
      // NOTE: This assumes the cell index corresponds directly to the
      //       "lower-left-front" node index for interpolation purposes.
      //       We need the local index relative to the cornerPtr array,
      //       which is indexed from 0 based on owned NODES.
    
      PetscInt iCell_local = iCell_global - xs_node;
      PetscInt jCell_local = jCell_global - ys_node;
      PetscInt kCell_local = kCell_global - zs_node;

      // --- Clamp LOCAL indices and check bounds for the BASE cell index ---
      // This check ensures the base index (i,j,k) for the interpolation
      // is within the bounds of the locally allocated cornerPtr array.
      // The interpolation function will handle the +1 offsets.
    // Boundary clamp: adjust indices to be within [0, mx), [0, my), [0, mz) 
    if (iCell_local >= user->info.mx) iCell_local = user->info.mx - 1;
    if (jCell_local >= user->info.my) jCell_local = user->info.my - 1;
    if (kCell_local >= user->info.mz) kCell_local = user->info.mz - 1;
    
    if (iCell_local < 0 ||
	jCell_local < 0 ||
	kCell_local < 0 ||
        iCell_local >= user->info.xm-1 ||
        jCell_local >= user->info.ym-1 ||
        kCell_local >= user->info.zm-1)
    {

      LOG_ALLOW(LOCAL,LOG_DEBUG,"[Rank %d] Rank/Domain boundary reached(cell wise): (k,j,i) - (%d,%d,%d).\n",rank,kCell_local,jCell_local,iCell_local);
      
      // Out-of-range: set output to zero       
      if (bs == 1) {
        ((PetscReal*)swarmOut)[p] = 0.0;
      } else {
        ((PetscReal*)swarmOut)[3*p + 0] = 0.0;
        ((PetscReal*)swarmOut)[3*p + 1] = 0.0;
        ((PetscReal*)swarmOut)[3*p + 2] = 0.0;
      }
      continue;
    }
    
    
    // Retrieve Interpolatetion coefficients (a1, a2, a3) for this particle 
    PetscReal alpha1 = weights[3*p + 0];
    PetscReal alpha2 = weights[3*p + 1];
    PetscReal alpha3 = weights[3*p + 2];
      */
    /* (Optional) If your final Interpolatetion expects a corner offset, adjust here.
       For example, if the cell center corresponds to the average of corners at (iCell,jCell,kCell)
       and (iCell+1, jCell+1, kCell+1), you might add +1. For now, we pass indices as is.
    */
     /*  
    PetscInt iUse = iCell_local;  // + 1 if required
    PetscInt jUse = jCell_local;  // + 1 if required
    PetscInt kUse = kCell_local;  // + 1 if required
     */
    /* (H) Call the per-particle Interpolatetion function.
       This function will use your _Generic macro (TrilinearInterpolation or PiecewiseLinearInterpolation)
       on the corner data.
    */
     /* 
    ierr = InterpolateEulerFieldToSwarmForParticle(
              fieldName,    // e.g., "Ucat" 
              cornerPtr,    // typed Pointer: (PetscReal***) or (Cmpnts***) 
              iUse, jUse, kUse,
              alpha1, alpha2, alpha3,
              swarmOut,     // Pointer to swarm output array 
              p,            // particle index 
              bs);          // block size: 1 or 3 
    CHKERRQ(ierr);
  }

  // (I) Restore swarm fields 
  ierr = DMSwarmRestoreField(swarm, swarmOutFieldName, NULL, NULL, &swarmOut);        CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "weight",          NULL, NULL, (void**)&weights); CHKERRQ(ierr);
  ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID",  NULL, NULL, (void**)&cellIDs); CHKERRQ(ierr);

  // (J) Deallocate the corner array using the generic deallocation macro 
  if (bs == 1) {
    PetscReal ***cornerScal = (PetscReal ***) cornerPtr;
    ierr = Deallocate3DArray(cornerScal, info.zm, info.ym); CHKERRQ(ierr);
  } else {
    Cmpnts ***cornerVec = (Cmpnts ***) cornerPtr;
    ierr = Deallocate3DArray(cornerVec, info.zm, info.ym); CHKERRQ(ierr);
  }

  LOG_ALLOW(GLOBAL, LOG_INFO,
    "InterpolateEulerFieldToSwarm: Rank %d Completed Interpolatetion of field='%s' for %d local particles.\n",
	    rank,fieldName, nLocal);

  PetscFunctionReturn(0);
}
*/

 /**
 * @brief Interpolates all relevant fields from the DMDA to the DMSwarm.
 *
 * Currently, it Interpolatetes:
 *   - user->Ucat (vector field) into the DMSwarm field "swarmVelocity".
 *
 * To add more fields, duplicate the call to InterpolateOneFieldOverSwarm and provide:
 *   - The global Vec for that field (e.g. user->Tcat for temperature),
 *   - A human-readable field name (for logging),
 *   - A DMSwarm output field name (e.g. "swarmTemperature").
 *
 * @param[in,out] user Pointer to a UserCtx containing:
 *                     - user->da (DM for the grid),
 *                     - user->swarm (DMSwarm for particles),
 *                     - user->Ucat (Vec for the vector field),
 *                     - possibly more fields like user->Tcat, user->Pcat, etc.
 *
 * @return PetscErrorCode Returns 0 on success, non-zero on failure.
 */
PetscErrorCode InterpolateAllFieldsToSwarm(UserCtx *user)
{
  PetscErrorCode ierr;
  PetscMPIInt rank;
  PetscFunctionBegin;

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);

  /* 
     1) Interpolate the 'velocity' field (user->Ucat) into DMSwarm's "swarmVelocity".
        - The function InterpolateOneFieldOverSwarm() loops over all local particles,
          retrieves iCell, jCell, kCell, (a1,a2,a3) from the DMSwarm,
          and calls the appropriate trilinear Interpolatetion routine.

        - "velocity" => just a label used in logging or debugging.
        - "swarmVelocity" => the DMSwarm field name where we store the result
          (assume dof=3 if user->Ucat has blockSize=3).
  */

  LOG_ALLOW(LOCAL, LOG_DEBUG,
		 "InterpolateteAllFieldsToSwarm: Interpolation of ucat to velocity begins on rank %d.\n",rank);
  // Make sure to pass the *GLOBAL* Vector to the function below! 
  ierr = InterpolateEulerFieldToSwarm(user, user->lUcat, 
                                      "Ucat", 
                                      "velocity"); CHKERRQ(ierr);
  /* 
     2) (OPTIONAL) If you have more fields, you can Interpolatete them similarly.

     For example, if user->Tcat is a scalar Vec for "temperature" and the DMSwarm
     has a field "Temperature" (with dof=1), you could do:

         ierr = InterpolateOneFieldOverSwarm(user, user->Tcat,
                                             "temperature",
                                             "swarmTemperature"); CHKERRQ(ierr);

     For pressure:
         ierr = InterpolateOneFieldOverSwarm(user, user->Pcat,
                                             "pressure",
                                             "swarmPressure"); CHKERRQ(ierr);
     
     ...and so forth.
   */

  /* 
     3) Optionally, synchronize or log that all fields are done 
  */
  ierr = MPI_Barrier(PETSC_COMM_WORLD); CHKERRQ(ierr);
  LOG_ALLOW(LOCAL, LOG_INFO,
	    "[rank %d]Completed Interpolateting all fields to the swarm.\n",rank);

  PetscFunctionReturn(0);
}

/////////////////////// Scatter from particles to euler fields

/**
 * @brief Functions for scattering particle data (scalar or vector) onto
 *        Eulerian grid fields by averaging contributions within each cell.
 *
 * This file provides a modular set of functions to perform particle-to-grid
 * projection, specifically calculating cell-averaged quantities from particle properties.
 * It assumes a PETSc environment using DMDA for the grids and DMSwarm for particles.
 *
 * Key Features:
 * - Handles both scalar (DOF=1) and vector (DOF=3) particle fields.
 * - Uses a pre-calculated particle count vector (`ParticleCount`) for normalization.
 * - Implicitly determines the target Eulerian DM (typically `user->da` for scalars,
 *   `user->fda` for vectors) based on standard field names ("P", "Nvert", "Ucat", "Ucont").
 * - Modifies existing Eulerian field vectors (e.g., `user->P`, `user->Ucat`) in place.
 * - Provides a high-level wrapper function (`ScatterAllParticleFieldsToEulerFields`)
 *   to easily scatter a standard set of fields.
 * - Uses only the base `SETERRQ` macro for error reporting to maximize compiler compatibility.
 *
 * Dependencies:
 * - PETSc library (DMDA, DMSwarm, Vec, IS, PetscLog, etc.)
 * - A `UserCtx` struct (defined elsewhere, e.g., "userctx.h") containing pointers
 *   to relevant DMs (`da`, `fda`), Vecs (`ParticleCount`, `P`, `Nvert`, `Ucat`, etc.),
 *   and the `DMSwarm` object (`swarm`).
 * - A custom DMSwarm field named `"DMSwarm_CellID"` (blockSize=3, type=PETSC_INT)
 *   must be registered and populated with the local cell indices for each particle.
 * - Logging infrastructure (`LOG_ALLOW`, etc.) assumed to be defined elsewhere.
 *
 * @defgroup scatter_module Particle-to-Grid Scattering
 * @{
 */

//-----------------------------------------------------------------------------
// Internal Helper Modules (Lower-level building blocks)
//-----------------------------------------------------------------------------
/**
 * @defgroup scatter_module_internal Internal Scattering Helpers
 * @ingroup scatter_module
 * @brief Lower-level functions used by the main scattering routines.
 * @{
 */

/**
 * @brief Determines the target Eulerian DM and expected DOF for scattering a given particle field.
 *
 * Based on hardcoded rules mapping particle field names ("P", "Nvert", "Ucat", "Ucont")
 * to user context DMs (`user->da` or `user->fda`). This function encapsulates the
 * policy of where different fields should be scattered. Modify this function to
 * add rules for custom fields.
 *
 * @param[in] user             Pointer to the UserCtx containing required DMs (`da`, `fda`).
 * @param[in] particleFieldName Name of the particle field.
 * @param[out] targetDM        Pointer to store the determined target DM (`user->da` or `user->fda`).
 * @param[out] expected_dof    Pointer to store the expected DOF (1 or 3) for this field.
 *
 * @return PetscErrorCode Returns 0 on success. Error codes:
 *         - `PETSC_ERR_ARG_NULL` if required inputs are NULL.
 *         - `PETSC_ERR_ARG_WRONG` if `particleFieldName` is not recognized.
 */
PetscErrorCode GetScatterTargetInfo(UserCtx *user, const char *particleFieldName,
                                    DM *targetDM, PetscInt *expected_dof)
{
    char           msg[ERROR_MSG_BUFFER_SIZE]; // Buffer for formatted error messages
    PetscFunctionBeginUser;

    // --- Input Validation ---
    // Check for NULL pointers in essential inputs
    if (!user) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx pointer is NULL.");
    if (!user->da) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx->da is NULL.");
    if (!user->fda) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx->fda is NULL."); // Needed for vector fields
    if (!particleFieldName) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "particleFieldName is NULL.");
    if (!targetDM) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Output targetDM pointer is NULL.");
    if (!expected_dof) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Output expected_dof pointer is NULL.");

    // --- Determine Target DM and DOF based on Field Name ---
    // Compare the input field name with known scalar fields targeting 'da'
    if (strcmp(particleFieldName, "P") == 0 || strcmp(particleFieldName, "Nvert") == 0) {
        *expected_dof = 1;      // Scalar fields have DOF 1
        *targetDM = user->da;   // Target the primary scalar DMDA
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "GetScatterTargetInfo: Field '%s' targets DM 'da' (DOF=1).\n", particleFieldName);
    }
    // Compare with known vector fields targeting 'fda'
    else if (strcmp(particleFieldName, "Ucat") == 0 || strcmp(particleFieldName, "Ucont") == 0) {
        *expected_dof = 3;      // Vector fields have DOF 3
        *targetDM = user->fda;  // Target the vector DMDA (often node-based)
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "GetScatterTargetInfo: Field '%s' targets DM 'fda' (DOF=3).\n", particleFieldName);
    }
    // --- Add rules for other fields here ---
    // else if (strcmp(particleFieldName, "SomeOtherScalar") == 0) { *expected_dof = 1; *targetDM = user->da; }
    // else if (strcmp(particleFieldName, "SomeOtherVector") == 0) { *expected_dof = 3; *targetDM = user->someOtherDM; }
    else {
        // The provided field name doesn't match any known rules
        *targetDM = NULL; // Indicate failure
        *expected_dof = 0;
        // Format the error message manually
        PetscSNPrintf(msg, sizeof(msg), "GetScatterTargetInfo: Field name '%s' is not recognized for automatic DM selection.", particleFieldName);
        // Use SETERRQ with the formatted message and appropriate error code
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, msg); // Use WRONG argument error code
    }

    PetscFunctionReturn(0);
}


/**
 * @brief Accumulates a particle field (scalar or vector) into a target grid sum vector.
 *
 * This function iterates through local particles, identifies their cell using the
 * `"DMSwarm_CellID"` field, and adds the particle's field value (`particleFieldName`)
 * to the corresponding cell location in the `gridSumVec`. It handles both scalar
 * (DOF=1) and vector (DOF=3) fields automatically based on the DOF of `gridSumDM`.
 *
 * IMPORTANT: The caller must ensure `gridSumVec` is zeroed before calling this
 * function if a fresh sum calculation is desired.
 *
 * @param[in] swarm           The DMSwarm containing particles.
 * @param[in] particleFieldName Name of the field on the particles (must match DOF of `gridSumDM`).
 * @param[in] gridSumDM       The DMDA associated with `gridSumVec`. Its DOF determines
 *                            how many components are accumulated.
 * @param[in,out] gridSumVec  The Vec (associated with `gridSumDM`) to accumulate sums into.
 *
 * @return PetscErrorCode 0 on success. Errors if fields don't exist, DMs are incompatible,
 *         or memory access fails.
 */
PetscErrorCode AccumulateParticleField(DM swarm, const char *particleFieldName,
                                       DM gridSumDM, Vec gridSumVec)
{
    PetscErrorCode    ierr;
    PetscInt          dof;                   // DOF determined from gridSumDM
    PetscInt          nlocal, p;             // Local particle count and loop index
    const PetscReal   *particle_arr = NULL;  // Pointer to particle field data array (assuming Real)
    const PetscInt  *cell_id_arr = NULL;   // Pointer to particle cell ID array ("DMSwarm_CellID", Int)
    PetscScalar       *sum_arr_ptr = NULL;   // Pointer to grid sum vector data array (Scalar)
    PetscInt          gxs, gys, gzs;         // Start indices of local ghosted patch (often 0)
    PetscInt          gxm, gym, gzm;         // Dimensions of local ghosted patch (including ghosts)
    PetscMPIInt       rank;                  // MPI rank for logging
    char              msg[ERROR_MSG_BUFFER_SIZE]; // Buffer for formatted error messages

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // --- Get DOF from the target grid DM ---
    ierr = DMDAGetInfo(gridSumDM, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &dof, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
    // Validate that the DOF is supported (currently 1 or 3)
    if (dof != 1 && dof != 3) {
        PetscSNPrintf(msg, sizeof(msg), "AccumulateParticleField: gridSumDM DOF must be 1 or 3, got %d.", dof);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, msg);
    }

    // --- Get Particle Data Arrays ---
    // DMSwarmGetField will return an error if the field doesn't exist, caught by CHKERRQ.
    ierr = DMSwarmGetField(swarm, particleFieldName, NULL, NULL, (void **)&particle_arr); CHKERRQ(ierr);
    ierr = DMSwarmGetField(swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_id_arr); CHKERRQ(ierr);

    // Get number of local particles *after* successfully getting field pointers
    ierr = DMSwarmGetLocalSize(swarm, &nlocal); CHKERRQ(ierr);

    // --- Get Grid Sum Vector Array & Dimensions ---
    ierr = VecGetArray(gridSumVec, &sum_arr_ptr); CHKERRQ(ierr);
    // Get dimensions needed for calculating flat index within the local ghosted array
    ierr = DMDAGetGhostCorners(gridSumDM, &gxs, &gys, &gzs, &gxm, &gym, &gzm); CHKERRQ(ierr);

    // --- Accumulate Locally ---
    LOG_ALLOW(LOCAL, LOG_DEBUG, "AccumulateParticleField (Rank %d): Accumulating '%s' (DOF=%d) from %d particles using CellID field 'DMSwarm_CellID'.\n", rank, particleFieldName, dof, nlocal);
    // Loop over all particles owned by this process
    for (p = 0; p < nlocal; ++p) {
        // Extract local cell indices (relative to start of ghosted patch, [0..gxm-1], etc.)
        // Assumes DMSwarm_CellID stores (i, j, k) contiguously for each particle.
        PetscInt pidx = cell_id_arr[p * 3 + 0]; // Local i-index
        PetscInt pidy = cell_id_arr[p * 3 + 1]; // Local j-index
        PetscInt pidz = cell_id_arr[p * 3 + 2]; // Local k-index

        // Bounds Check: Ensure the particle's cell index is within the valid local ghosted region
        if (pidx >= 0 && pidx < gxm && pidy >= 0 && pidy < gym && pidz >= 0 && pidz < gzm)
        {
             // Calculate the flat 1D index for this cell within the linear ghosted array
             // Uses PETSc's standard C-style row-major ordering (k-slowest, j-middle, i-fastest)
             // Corrected: k (slowest), j, i (fastest)
             PetscInt cell_flat_idx = (pidz * gym + pidy) * gxm + pidx;

             // Calculate the base index for this particle's data in particle_arr
             PetscInt particle_base_idx = p * dof;
             // Calculate the base index for this cell's data in sum_arr_ptr
             PetscInt grid_base_idx = cell_flat_idx * dof;

             // Add particle components to the grid sum vector components
             for (PetscInt c = 0; c < dof; ++c) {
                 sum_arr_ptr[grid_base_idx + c] += particle_arr[particle_base_idx + c];
             }
        } else {
             // Log a warning if a particle's CellID is outside the expected local region.
             // This might indicate particles needing migration or boundary issues.
             LOG_ALLOW(LOCAL, LOG_WARNING, "AccumulateParticleField (Rank %d): Particle %d (field '%s') has out-of-bounds CellID (%d, %d, %d). Ghosted dims: %dx%dx%d. Skipping.\n",
                      rank, p, particleFieldName, pidx, pidy, pidz, gxm, gym, gzm);
        }
    } // End of particle loop

    // --- Restore Access to Arrays ---
    ierr = VecRestoreArray(gridSumVec, &sum_arr_ptr); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, particleFieldName, NULL, NULL, (void **)&particle_arr); CHKERRQ(ierr);
    ierr = DMSwarmRestoreField(swarm, "DMSwarm_CellID", NULL, NULL, (void **)&cell_id_arr); CHKERRQ(ierr);

    // --- Assemble Global Sum Vector ---
    // Crucial for parallel execution: sums contributions for cells shared across process boundaries.
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "AccumulateParticleField: Assembling global sum vector for '%s'.\n", particleFieldName);
    ierr = VecAssemblyBegin(gridSumVec); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(gridSumVec); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}


/**
 * @brief Normalizes a grid vector of sums by a grid vector of counts to produce an average.
 *
 * Calculates `avgVec[i] = sumVec[i] / countVec[i]` for each component of each
 * cell owned by the current process, provided `countVec[i] > 0`. Otherwise, sets `avgVec[i] = 0`.
 * Handles both scalar (DOF=1) and vector (DOF=3) data fields based on the DOF of `dataDM`.
 * Uses DMDA multi-dimensional array accessors (`DMDAVecGetArray...`) for safe and convenient indexing.
 *
 * @param[in] countDM    The DMDA associated with `countVec` (must have DOF=1).
 * @param[in] countVec   The Vec containing particle counts per cell (read-only).
 * @param[in] dataDM     The DMDA associated with `sumVec` and `avgVec` (must have DOF=1 or DOF=3).
 * @param[in] sumVec     The Vec containing the accumulated sums per cell (read-only).
 * @param[in,out] avgVec The Vec where the calculated averages will be stored (overwritten). Must be
 *                       associated with `dataDM`.
 *
 * @return PetscErrorCode 0 on success. Errors on incompatible DMs or memory access failure.
 */
PetscErrorCode NormalizeGridVectorByCount(DM countDM, Vec countVec,
                                          DM dataDM, Vec sumVec, Vec avgVec)
{
    PetscErrorCode ierr;
    PetscInt       data_dof;
    PetscInt       count_dof;
    PetscMPIInt    rank;
    char           msg[ERROR_MSG_BUFFER_SIZE];

    // Pointers for DMDA array accessors - declare specific types
    PetscScalar    ***count_arr_3d = NULL;      // For DOF=1 count vector (3D DMDA)
    PetscScalar    ***sum_arr_scalar = NULL;    // For DOF=1 sum vector (3D DMDA)
    PetscScalar    ***avg_arr_scalar = NULL;    // For DOF=1 avg vector (3D DMDA)
    PetscScalar    ***sum_arr_vector = NULL;   // For DOF=3 sum vector (3D DMDA + DOF)
    PetscScalar    ***avg_arr_vector = NULL;   // For DOF=3 avg vector (3D DMDA + DOF)


    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // --- Validation ---
    ierr = DMDAGetInfo(countDM, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &count_dof, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
    ierr = DMDAGetInfo(dataDM, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &data_dof, NULL, NULL, NULL, NULL, NULL); CHKERRQ(ierr);
    if (count_dof != 1) { PetscSNPrintf(msg, sizeof(msg), "countDM must have DOF=1, got %d.", count_dof); SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, msg); }
    if (data_dof != 1 && data_dof != 3) { PetscSNPrintf(msg, sizeof(msg), "dataDM DOF must be 1 or 3, got %d.", data_dof); SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, msg); }

    // --- Get Array Access using appropriate DMDA accessors ---
    ierr = DMDAVecGetArrayRead(countDM, countVec, &count_arr_3d); CHKERRQ(ierr);

    if (data_dof == 1) {
        ierr = DMDAVecGetArrayRead(dataDM, sumVec, &sum_arr_scalar); CHKERRQ(ierr);
        ierr = DMDAVecGetArray(dataDM, avgVec, &avg_arr_scalar); CHKERRQ(ierr);
    } else { // data_dof == 3
        ierr = DMDAVecGetArrayDOFRead(dataDM, sumVec, &sum_arr_vector); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayDOF(dataDM, avgVec, &avg_arr_vector); CHKERRQ(ierr);
    }

    // Get the corners (global start indices) and dimensions of the *local owned* region
    PetscInt xs, ys, zs, xm, ym, zm;
    ierr = DMDAGetCorners(countDM, &xs, &ys, &zs, &xm, &ym, &zm); CHKERRQ(ierr);

    // --- Normalize Over Owned Cells ---
    LOG_ALLOW(LOCAL, LOG_DEBUG, "NormalizeGridVectorByCount (Rank %d): Normalizing DOF=%d data over owned range [%d:%d, %d:%d, %d:%d].\n",
              rank, data_dof, xs, xs+xm, ys, ys+ym, zs, zs+zm);

    // Loop using GLOBAL indices (i, j, k) over the range owned by this process
    for (PetscInt k = zs; k < zs + zm; ++k) {
        for (PetscInt j = ys; j < ys + ym; ++j) {
            for (PetscInt i = xs; i < xs + xm; ++i) {

                // Access the count using standard 3D indexing
                PetscScalar count = count_arr_3d[k][j][i];

                if (PetscRealPart(count) > 0.5) { // Use tolerance for float comparison
                    if (data_dof == 1) {
                        // Access scalar sum/avg using standard 3D indexing
                        avg_arr_scalar[k][j][i] = sum_arr_scalar[k][j][i] / count;
                    } else { // data_dof == 3
                        // Access vector components using DOF indexing on the last dimension
                        for (PetscInt c = 0; c < data_dof; ++c) {
                            avg_arr_vector[k][j][i * data_dof + c] = sum_arr_vector[k][j][i * data_dof + c] / count;
                        }
                    }
                } else { // count is zero or negative
                    // Set average to zero
                    if (data_dof == 1) {
                        avg_arr_scalar[k][j][i] = 0.0;
                    } else { // data_dof == 3
                        for (PetscInt c = 0; c < data_dof; ++c) {
                            avg_arr_vector[k][j][i * data_dof + c] = 0.0;
                        }
                    }
                } // end if count > 0.5
            } // end i loop
        } // end j loop
    } // end k loop

    // --- Restore Arrays using appropriate functions ---
    ierr = DMDAVecRestoreArrayRead(countDM, countVec, &count_arr_3d); CHKERRQ(ierr);
    if (data_dof == 1) {
        ierr = DMDAVecRestoreArrayRead(dataDM, sumVec, &sum_arr_scalar); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(dataDM, avgVec, &avg_arr_scalar); CHKERRQ(ierr);
    } else { // data_dof == 3
        ierr = DMDAVecRestoreArrayDOFRead(dataDM, sumVec, &sum_arr_vector); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayDOF(dataDM, avgVec, &avg_arr_vector); CHKERRQ(ierr);
    }

    // --- Assemble Final Average Vector ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "NormalizeGridVectorByCount: Assembling final average vector (DOF=%d).\n", data_dof);
    ierr = VecAssemblyBegin(avgVec); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(avgVec); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/** @} */ // End of scatter_module_internal group

//-----------------------------------------------------------------------------
// User-Facing API
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// MODULE 4: Internal Scatter Orchestration Helper - No PetscErrorClear
//-----------------------------------------------------------------------------
/**
 * @brief Internal helper function to orchestrate the scatter operation (accumulate + normalize).
 * @ingroup scatter_module_internal
 *
 * Manages the temporary sum vector and calls the accumulation and normalization
 * functions. Assumes caller determined target DM and DOF. Checks for particle field existence.
 * NOTE: If DMSwarmGetField fails, the subsequent SETERRQ will overwrite the original error.
 *
 * @param[in] user                 Pointer to UserCtx containing swarm, ParticleCount, da.
 * @param[in] particleFieldName    Name of the field in the DMSwarm.
 * @param[in] targetDM             The DMDA where the final average and intermediate sum reside.
 * @param[in] expected_dof         The expected DOF (1 or 3) for the targetDM and field.
 * @param[in,out] eulerFieldAverageVec The pre-created Vec associated with targetDM to store the result.
 *
 * @return PetscErrorCode 0 on success. Errors if particle field doesn't exist or
 *         underlying helpers fail.
 */
static PetscErrorCode ScatterParticleFieldToEulerField_Internal(UserCtx *user,
                                                         const char *particleFieldName,
                                                         DM targetDM,
                                                         PetscInt expected_dof,
                                                         Vec eulerFieldAverageVec)
{
    PetscErrorCode ierr;
    Vec            sumVec = NULL;
    char           msg[ERROR_MSG_BUFFER_SIZE]; // Buffer for formatted error messages

    PetscFunctionBeginUser;

    if (!user || !user->swarm || !user->ParticleCount || !particleFieldName || !targetDM || !eulerFieldAverageVec)
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "NULL input provided to ScatterParticleFieldToEulerField_Internal.");

    // --- Check if Particle Field Exists ---
    // Attempt a GetField call; if it fails, the field doesn't exist.
    // We let CHKERRQ handle the error directly if the field doesn't exist OR
    // we catch it specifically to provide a more tailored message.

    /*
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"Field %s being accessed to check existence \n",particleFieldName);
    ierr = DMSwarmGetField(user->swarm, particleFieldName, NULL, NULL, NULL);
    if (ierr) { // If GetField returns an error
        PetscSNPrintf(msg, sizeof(msg), "Particle field '%s' not found in DMSwarm for scattering.", particleFieldName);
        // Directly set the error, overwriting the one from GetField
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, msg);
    }
    ierr = DMSwarmRestoreField(user->swarm, particleFieldName, NULL, NULL, NULL);
    */
    
    // --- Setup Temporary Sum Vector ---
    ierr = VecDuplicate(eulerFieldAverageVec, &sumVec); CHKERRQ(ierr);
    ierr = VecSet(sumVec, 0.0); CHKERRQ(ierr);
    ierr = PetscSNPrintf(msg, sizeof(msg), "TempSum_%s", particleFieldName); CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject)sumVec, msg); CHKERRQ(ierr);

    // --- Accumulate ---
    // This will call DMSwarmGetField again. If it failed above, it will likely fail here too,
    // unless the error was cleared somehow between the check and here (unlikely).
    // If the check above was skipped (Option 1), this is where the error for non-existent
    // field will be caught by CHKERRQ.
    ierr = AccumulateParticleField(user->swarm, particleFieldName, targetDM, sumVec); CHKERRQ(ierr);

    // Calculate the number of particles per cell.
    ierr = CalculateParticleCountPerCell(user);
    // --- Normalize ---
    ierr = NormalizeGridVectorByCount(user->da, user->ParticleCount, targetDM, sumVec, eulerFieldAverageVec); CHKERRQ(ierr);

    // --- Cleanup ---
    ierr = VecDestroy(&sumVec); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

/**
 * @brief Scatters a particle field (scalar or vector) to the corresponding Eulerian field average.
 * @ingroup scatter_module
 *
 * This is the primary user-facing function for scattering a single field.
 * It determines the target Eulerian DM (`user->da` or `user->fda`) based on the
 * `particleFieldName`, validates that the provided `eulerFieldAverageVec` is compatible
 * with that target DM, and then orchestrates the scatter operation by calling
 * an internal helper function. The final averaged result is stored IN-PLACE in
 * the `eulerFieldAverageVec`.
 *
 * @param[in] user                 Pointer to UserCtx containing `da`, `fda`, `swarm`, `ParticleCount`.
 * @param[in] particleFieldName    Name of the field in the DMSwarm (e.g., "P", "Ucat").
 * @param[in,out] eulerFieldAverageVec Pre-created Vec associated with the correct target DM
 *                                 (implicitly `da` or `fda`). Result stored here.
 *                                 This vector should be zeroed by the caller if desired
 *                                 before this function (e.g., using `VecSet`). The wrapper
 *                                 `ScatterAllParticleFieldsToEulerFields` handles this zeroing.
 *
 * @return PetscErrorCode 0 on success. Errors on NULL input, unrecognized field name,
 *         incompatible target vector, or if the particle field doesn't exist.
 */
PetscErrorCode ScatterParticleFieldToEulerField(UserCtx *user,
                                                const char *particleFieldName,
                                                Vec eulerFieldAverageVec)
{
    PetscErrorCode ierr;
    DM             targetDM = NULL;       // Will point to user->da or user->fda
    PetscInt       expected_dof = 0;      // Will be 1 or 3
    char           msg[ERROR_MSG_BUFFER_SIZE]; // Buffer for formatted error messages

    PetscFunctionBeginUser;

    // --- Essential Input Validation ---
    if (!user) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx pointer is NULL.");
    if (!user->swarm) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx->swarm is NULL.");
    if (!user->ParticleCount) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "UserCtx->ParticleCount is NULL.");
    if (!eulerFieldAverageVec) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Output eulerFieldAverageVec is NULL.");
    // particleFieldName validity checked within GetScatterTargetInfo

    // --- Determine Target DM & DOF using the helper function ---
    ierr = GetScatterTargetInfo(user, particleFieldName, &targetDM, &expected_dof); CHKERRQ(ierr);
    // If GetScatterTargetInfo failed (e.g., unknown name), CHKERRQ would have returned.

    // --- Validate the provided Target Vec's Compatibility ---
    DM vec_dm;
    PetscInt vec_dof;
    // Check that the provided average vector has a DM associated with it
    ierr = VecGetDM(eulerFieldAverageVec, &vec_dm); CHKERRQ(ierr);
    if (!vec_dm) {
         PetscSNPrintf(msg, sizeof(msg), "Provided eulerFieldAverageVec for field '%s' does not have an associated DM.", particleFieldName);
         SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, msg);
    }
    // Get the block size (DOF) of the provided vector
    ierr = VecGetBlockSize(eulerFieldAverageVec, &vec_dof); CHKERRQ(ierr);
    // Compare the vector's associated DM with the one determined by the field name
    if (vec_dm != targetDM) {
        const char *target_dm_name = "targetDM", *vec_dm_name = "vec_dm";
        // Get actual names if possible for a more informative error message
        PetscObjectGetName((PetscObject)targetDM, &target_dm_name);
        PetscObjectGetName((PetscObject)vec_dm, &vec_dm_name);
        PetscSNPrintf(msg, sizeof(msg), "Provided eulerFieldAverageVec associated with DM '%s', but field '%s' requires scatter to DM '%s'.", vec_dm_name, particleFieldName, target_dm_name);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_INCOMP, msg);
    }
    // Compare the vector's DOF with the one expected for the field name
    if (vec_dof != expected_dof) {
         PetscSNPrintf(msg, sizeof(msg), "Field '%s' requires DOF %d, but provided eulerFieldAverageVec has DOF %d.", particleFieldName, expected_dof, vec_dof);
         SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_INCOMP, msg);
    }

    // --- Perform Scatter using Internal Helper ---
    // Log intent before calling the core logic
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "ScatterParticleFieldToEulerField: Scattering field '%s' (DOF=%d).\n", particleFieldName, expected_dof);
    ierr = ScatterParticleFieldToEulerField_Internal(user, // Pass user context
                                                     particleFieldName, // Name of particle field
                                                     targetDM, // Determined target DM (da or fda)
                                                     expected_dof, // Determined DOF (1 or 3)
                                                     eulerFieldAverageVec); // The output vector
    CHKERRQ(ierr); // Handle potential errors from the internal function

    LOG_ALLOW(GLOBAL, LOG_INFO, "ScatterParticleFieldToEulerField: Successfully scattered field '%s'.\n", particleFieldName);
    PetscFunctionReturn(0);
}


/**
 * @brief Scatters a predefined set of particle fields to their corresponding Eulerian fields.
 * @ingroup scatter_module
 *
 * This convenience function calls the unified `ScatterParticleFieldToEulerField`
 * for a standard set of fields (currently just "P", others commented out). It assumes
 * the target Eulerian Vec objects (e.g., `user->P`, `user->Ucat`) exist in the
 * UserCtx structure and are correctly associated with their respective DMs (`user->da`
 * or `user->fda`).
 *
 * **IMPORTANT:** It zeros the target Vecs before scattering. Ensure `user->ParticleCount`
 * has been computed accurately *before* calling this function, reflecting the current
 * particle distribution.
 *
 * @param[in,out] user Pointer to the UserCtx structure containing all required DMs,
 *                     Vecs (`ParticleCount`, target Eulerian fields like `P`, `Ucat`), and `swarm`.
 *
 * @return PetscErrorCode 0 on success. Errors if prerequisites (like ParticleCount)
 *         are missing or if underlying scatter calls fail (e.g., particle field missing).
 */
PetscErrorCode ScatterAllParticleFieldsToEulerFields(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting scattering of specified particle fields to Eulerian grids.\n");

    // --- Pre-computation Check: Ensure Particle Counts are Ready ---
    if (!user->ParticleCount) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "UserCtx->ParticleCount is NULL. Compute counts before calling ScatterAllParticleFieldsToEulerFields.");
    }

    // --- Scatter Particle Field "P" -> Eulerian Field user->P (on da) ---
    // Check if the target Eulerian vector 'user->P' exists.
    if (user->P) {
        
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Scattering particle field 'P' to user->P.\n");
        // Zero the target vector before accumulating the new average for this step/call.

	// Debug Verification ------------------------------------------------
	Vec swarm_P;
	PetscReal Avg_P,Avg_swarm_P;
	     
	ierr = VecMean(user->P,&Avg_P);
	LOG_ALLOW(GLOBAL,LOG_DEBUG," Average of Pressure before scatter: %.4f.\n",Avg_P);
	     
	ierr = DMSwarmCreateGlobalVectorFromField(user->swarm,"P",&swarm_P);
	ierr = VecMean(swarm_P,&Avg_swarm_P);

	LOG_ALLOW(GLOBAL,LOG_DEBUG," Average of Particle Pressure: %.4f.\n",Avg_swarm_P);

	ierr = DMSwarmDestroyGlobalVectorFromField(user->swarm,"P",&swarm_P);
	// Debug----------------------------------------------------------------
	  
	//ierr = VecSet(user->P, 0.0); CHKERRQ(ierr);
        // Call the unified scatter function. It will handle DM determination and validation.
        // It will also error out if the *particle* field "P" doesn't exist in the swarm.
        ierr = ScatterParticleFieldToEulerField(user, "P", user->P); CHKERRQ(ierr);
	ierr = VecMean(user->P,&Avg_P);
	
	LOG_ALLOW(GLOBAL,LOG_DEBUG," Average of Pressure after  scatter: %.4f.\n",Avg_P);
    } else {
        // Only log a warning if the target Eulerian field is missing in the context.
        LOG_ALLOW(GLOBAL, LOG_WARNING, "Skipping scatter for 'P': UserCtx->P is NULL.\n");
    }

    // --- (Commented Out) Scatter Other Fields ---
    // To enable scattering for other fields, uncomment the relevant block
    // AND ensure:
    // 1. The corresponding particle field (e.g., "Nvert", "Ucat") exists in user->swarm.
    // 2. The corresponding target Eulerian Vec (e.g., user->Nvert, user->Ucat) exists in user->ctx.
    // 3. The target Eulerian Vec is associated with the correct DM (da for Nvert, fda for Ucat).

    /*
    // Scatter Particle Field "Nvert" -> Eulerian Field user->Nvert (on da)
    if (user->Nvert) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Scattering particle field 'Nvert' to user->Nvert.\n");
        ierr = VecSet(user->Nvert, 0.0); CHKERRQ(ierr);
        ierr = ScatterParticleFieldToEulerField(user, "Nvert", user->Nvert); CHKERRQ(ierr);
    } else {
         LOG_ALLOW(GLOBAL, LOG_WARNING, "Skipping scatter for 'Nvert': UserCtx->Nvert is NULL.\n");
    }
    */

    /*
    // Scatter Particle Field "Ucat" -> Eulerian Field user->Ucat (on fda)
    if (user->Ucat) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Scattering particle field 'Ucat' to user->Ucat.\n");
        ierr = VecSet(user->Ucat, 0.0); CHKERRQ(ierr);
        ierr = ScatterParticleFieldToEulerField(user, "Ucat", user->Ucat); CHKERRQ(ierr);
    } else {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "Skipping scatter for 'Ucat': UserCtx->Ucat is NULL.\n");
    }
    */

    /*
     // Scatter Particle Field "Ucont" -> Eulerian Field user->Ucont (on fda)
    if (user->Ucont) {
        LOG_ALLOW(GLOBAL, LOG_DEBUG, "Scattering particle field 'Ucont' to user->Ucont.\n");
        ierr = VecSet(user->Ucont, 0.0); CHKERRQ(ierr);
        ierr = ScatterParticleFieldToEulerField(user, "Ucont", user->Ucont); CHKERRQ(ierr);
    } else {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "Skipping scatter for 'Ucont': UserCtx->Ucont is NULL.\n");
    }
    */

    // Add more fields as needed following the pattern above...

    LOG_ALLOW(GLOBAL, LOG_INFO, "Finished scattering specified particle fields.\n");
    PetscFunctionReturn(0);
}

/** @} */ // End of scatter_module group

/**
 * @brief Interpolates a scalar field from corner nodes to all face centers.
 *
 * This routine computes the average of the four corner-node values
 * defining each face of a hexahedral cell:
 *   - X-faces (perpendicular to X): face between (i-1,i) in X-dir
 *   - Y-faces (perpendicular to Y): face between (j-1,j) in Y-dir
 *   - Z-faces (perpendicular to Z): face between (k-1,k) in Z-dir
 *
 * @param[in]  corner_arr  Ghosted node-centered array (global indexing) from user->fda.
 * @param[out] faceX_arr   Local array for X-faces sized [zm][ym][xm+1].
 * @param[out] faceY_arr   Local array for Y-faces sized [zm][ym+1][xm].
 * @param[out] faceZ_arr   Local array for Z-faces sized [zm+1][ym][xm].
 * @param[in]  user        User context containing DMDA 'fda' and GetOwnedCellRange.
 *
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode InterpolateCornerToFaceCenter_Scalar(
    PetscReal ***corner_arr,
    PetscReal ***faceX_arr,
    PetscReal ***faceY_arr,
    PetscReal ***faceZ_arr,
    UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo info;
    ierr = DMDAGetLocalInfo(user->fda, &info); CHKERRQ(ierr);

    // Determine owned-cell ranges based on corner-node ownership
    PetscInt xs, xm, ys, ym, zs, zm;
    ierr = GetOwnedCellRange(&info, 0, &xs, &xm); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&info, 1, &ys, &ym); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&info, 2, &zs, &zm); CHKERRQ(ierr);

    // Global exclusive end indices for cells
    PetscInt xe = xs + xm;
    PetscInt ye = ys + ym;
    PetscInt ze = zs + zm;

    // --- X‐faces: loops k=zs..ze-1, j=ys..ye-1, i=xs..xe (xm+1 faces per row) ---
    for (PetscInt k = zs; k < ze; ++k) {
        PetscInt k_loc = k - zs;
        for (PetscInt j = ys; j < ye; ++j) {
            PetscInt j_loc = j - ys;
            for (PetscInt i = xs; i <= xe; ++i) {
                PetscInt i_loc = i - xs; // 0..xm
                // Average the four corners of the Y-Z face at X = i
                PetscReal sum = corner_arr[k  ][j  ][i]
                              + corner_arr[k+1][j  ][i]
                              + corner_arr[k  ][j+1][i]
                              + corner_arr[k+1][j+1][i];
                faceX_arr[k_loc][j_loc][i_loc] = sum * 0.25;
            }
        }
    }

    // --- Y‐faces: loops k=zs..ze-1, j=ys..ye (ym+1 faces), i=xs..xe-1 ---
    for (PetscInt k = zs; k < ze; ++k) {
        PetscInt k_loc = k - zs;
        for (PetscInt j = ys; j <= ye; ++j) {
            PetscInt j_loc = j - ys; // 0..ym
            for (PetscInt i = xs; i < xe; ++i) {
                PetscInt i_loc = i - xs;
                // Average the four corners of the X-Z face at Y = j
                PetscReal sum = corner_arr[k  ][j][i  ]
                              + corner_arr[k+1][j][i  ]
                              + corner_arr[k  ][j][i+1]
                              + corner_arr[k+1][j][i+1];
                faceY_arr[k_loc][j_loc][i_loc] = sum * 0.25;
            }
        }
    }

    // --- Z‐faces: loops k=zs..ze (zm+1), j=ys..ye-1, i=xs..xe-1 ---
    for (PetscInt k = zs; k <= ze; ++k) {
        PetscInt k_loc = k - zs;
        for (PetscInt j = ys; j < ye; ++j) {
            PetscInt j_loc = j - ys;
            for (PetscInt i = xs; i < xe; ++i) {
                PetscInt i_loc = i - xs;
                // Average the four corners of the X-Y face at Z = k
                PetscReal sum = corner_arr[k][j  ][i  ]
                              + corner_arr[k][j  ][i+1]
                              + corner_arr[k][j+1][i  ]
                              + corner_arr[k][j+1][i+1];
                faceZ_arr[k_loc][j_loc][i_loc] = sum * 0.25;
            }
        }
    }

    return 0;
}

/**
 * @brief Interpolates a vector field from corner nodes to all face centers.
 *
 * Identical to the scalar version, except it averages each component of the
 * Cmpnts struct at the four corner-nodes per face.
 *
 * @param[in]  corner_arr  Ghosted 3-component array (global node indices).
 * @param[out] faceX_arr   Local array of Cmpnts for X-faces sized [zm][ym][xm+1].
 * @param[out] faceY_arr   Local array of Cmpnts for Y-faces sized [zm][ym+1][xm].
 * @param[out] faceZ_arr   Local array of Cmpnts for Z-faces sized [zm+1][ym][xm].
 * @param[in]  user        User context containing DMDA 'fda'.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InterpolateCornerToFaceCenter_Vector(
    Cmpnts ***corner_arr,
    Cmpnts ***faceX_arr,
    Cmpnts ***faceY_arr,
    Cmpnts ***faceZ_arr,
    UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo info;
    PetscMPIInt rank;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
                   "Rank %d starting InterpolateFieldFromCornerToFaceCenter_Vector.\n", rank);
    
    ierr = DMDAGetLocalInfo(user->fda, &info); CHKERRQ(ierr);

    PetscInt xs, xm, ys, ym, zs, zm;
    ierr = GetOwnedCellRange(&info, 0, &xs, &xm); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&info, 1, &ys, &ym); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(&info, 2, &zs, &zm); CHKERRQ(ierr);

    PetscInt xe = xs + xm;
    PetscInt ye = ys + ym;
    PetscInt ze = zs + zm;

    // X-faces
    for (PetscInt k = zs; k < ze; ++k) {
        PetscInt k_loc = k - zs;
        for (PetscInt j = ys; j < ye; ++j) {
            PetscInt j_loc = j - ys;
            for (PetscInt i = xs; i <= xe; ++i) {
                PetscInt i_loc = i - xs;
                Cmpnts sum = {0,0,0};
                sum.x = corner_arr[k  ][j  ][i].x + corner_arr[k+1][j  ][i].x
                      + corner_arr[k  ][j+1][i].x + corner_arr[k+1][j+1][i].x;
                sum.y = corner_arr[k  ][j  ][i].y + corner_arr[k+1][j  ][i].y
                      + corner_arr[k  ][j+1][i].y + corner_arr[k+1][j+1][i].y;
                sum.z = corner_arr[k  ][j  ][i].z + corner_arr[k+1][j  ][i].z
                      + corner_arr[k  ][j+1][i].z + corner_arr[k+1][j+1][i].z;
                faceX_arr[k_loc][j_loc][i_loc].x = sum.x * 0.25;
                faceX_arr[k_loc][j_loc][i_loc].y = sum.y * 0.25;
                faceX_arr[k_loc][j_loc][i_loc].z = sum.z * 0.25;
            }
        }
    }

    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
                   "Rank %d x-face Interpolation complete.\n", rank);
    
    // Y-faces
    for (PetscInt k = zs; k < ze; ++k) {
        PetscInt k_loc = k - zs;
        for (PetscInt j = ys; j <= ye; ++j) {
            PetscInt j_loc = j - ys;
            for (PetscInt i = xs; i < xe; ++i) {
                PetscInt i_loc = i - xs;
                Cmpnts sum = {0,0,0};
                sum.x = corner_arr[k  ][j][i  ].x + corner_arr[k+1][j][i  ].x
                      + corner_arr[k  ][j][i+1].x + corner_arr[k+1][j][i+1].x;
                sum.y = corner_arr[k  ][j][i  ].y + corner_arr[k+1][j][i  ].y
                      + corner_arr[k  ][j][i+1].y + corner_arr[k+1][j][i+1].y;
                sum.z = corner_arr[k  ][j][i  ].z + corner_arr[k+1][j][i  ].z
                      + corner_arr[k  ][j][i+1].z + corner_arr[k+1][j][i+1].z;
                faceY_arr[k_loc][j_loc][i_loc].x = sum.x * 0.25;
                faceY_arr[k_loc][j_loc][i_loc].y = sum.y * 0.25;
                faceY_arr[k_loc][j_loc][i_loc].z = sum.z * 0.25;
            }
        }
    }

    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
                   "Rank %d y-face Interpolation complete.\n", rank);
    
    // Z-faces
    for (PetscInt k = zs; k <= ze; ++k) {
        PetscInt k_loc = k - zs;
        for (PetscInt j = ys; j < ye; ++j) {
            PetscInt j_loc = j - ys;
            for (PetscInt i = xs; i < xe; ++i) {
                PetscInt i_loc = i - xs;
                Cmpnts sum = {0,0,0};
                sum.x = corner_arr[k][j  ][i  ].x + corner_arr[k][j  ][i+1].x
                      + corner_arr[k][j+1][i  ].x + corner_arr[k][j+1][i+1].x;
                sum.y = corner_arr[k][j  ][i  ].y + corner_arr[k][j  ][i+1].y
                      + corner_arr[k][j+1][i  ].y + corner_arr[k][j+1][i+1].y;
                sum.z = corner_arr[k][j  ][i  ].z + corner_arr[k][j  ][i+1].z
                      + corner_arr[k][j+1][i  ].z + corner_arr[k][j+1][i+1].z;
                faceZ_arr[k_loc][j_loc][i_loc].x = sum.x * 0.25;
                faceZ_arr[k_loc][j_loc][i_loc].y = sum.y * 0.25;
                faceZ_arr[k_loc][j_loc][i_loc].z = sum.z * 0.25;
            }
        }
    }

    LOG_ALLOW_SYNC(GLOBAL, LOG_DEBUG,
                   "Rank %d z-face Interpolation complete.\n", rank);
    
    return 0;
}
