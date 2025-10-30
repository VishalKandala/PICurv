#include "Boundaries.h" // The main header for our project
#include <string.h>    // For strcasecmp
#include <ctype.h>     // For isspace

#undef __FUNCT__
#define __FUNCT__ "CanRankServiceInletFace"
/**
 * @brief Determines if the current MPI rank owns any part of the globally defined inlet face,
 *        making it responsible for placing particles on that portion of the surface.
 *
 * The determination is based on the rank's owned nodes (from `DMDALocalInfo`) and
 * the global node counts, in conjunction with the `user->identifiedInletBCFace`.
 * A rank can service an inlet face if it owns the cells adjacent to that global boundary
 * and has a non-zero extent (owns cells) in the tangential dimensions of that face.
 *
 * @param user Pointer to the UserCtx structure, containing `identifiedInletBCFace`.
 * @param info Pointer to the DMDALocalInfo for the current rank's DA (node-based).
 * @param IM_nodes_global Global number of nodes in the I-direction (e.g., user->IM + 1 if user->IM is cell count).
 * @param JM_nodes_global Global number of nodes in the J-direction.
 * @param KM_nodes_global Global number of nodes in the K-direction.
 * @param[out] can_service_inlet_out Pointer to a PetscBool; set to PETSC_TRUE if the rank
 *                                   services (part of) the inlet, PETSC_FALSE otherwise.
 * @return PetscErrorCode 0 on success, non-zero on failure.
 */
PetscErrorCode CanRankServiceInletFace(UserCtx *user, const DMDALocalInfo *info,
                                              PetscInt IM_nodes_global, PetscInt JM_nodes_global, PetscInt KM_nodes_global,
                                              PetscBool *can_service_inlet_out)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank_for_logging; // For detailed debugging logs
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank_for_logging); CHKERRQ(ierr);

    *can_service_inlet_out = PETSC_FALSE; // Default to no service

    if (!user->inletFaceDefined) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d]: Inlet face not defined in user context. Cannot service.\n", rank_for_logging);
        PetscFunctionReturn(0);
    }

    // Get the range of cells owned by this rank in each dimension
    PetscInt owned_start_cell_i, num_owned_cells_on_rank_i;
    PetscInt owned_start_cell_j, num_owned_cells_on_rank_j;
    PetscInt owned_start_cell_k, num_owned_cells_on_rank_k;

    ierr = GetOwnedCellRange(info, 0, &owned_start_cell_i, &num_owned_cells_on_rank_i); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 1, &owned_start_cell_j, &num_owned_cells_on_rank_j); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 2, &owned_start_cell_k, &num_owned_cells_on_rank_k); CHKERRQ(ierr);

    // Determine the global index of the last cell (0-indexed) in each direction.
    // Example: If IM_nodes_global = 11 (nodes 0-10), there are 10 cells (0-9). Last cell index is 9.
    // Formula: global_nodes - 1 (num cells) - 1 (0-indexed) = global_nodes - 2.
    PetscInt last_global_cell_idx_i = (IM_nodes_global > 1) ? (IM_nodes_global - 2) : -1; // -1 if 0 or 1 node (i.e., 0 cells)
    PetscInt last_global_cell_idx_j = (JM_nodes_global > 1) ? (JM_nodes_global - 2) : -1;
    PetscInt last_global_cell_idx_k = (KM_nodes_global > 1) ? (KM_nodes_global - 2) : -1;

    switch (user->identifiedInletBCFace) {
        case BC_FACE_NEG_X: // Inlet on the global I-minimum face (face of cell C_i=0)
            // Rank services if its first owned node is global node 0 (info->xs == 0),
            // and it owns cells in I, J, and K directions.
            if (info->xs == 0 && num_owned_cells_on_rank_i > 0 &&
                num_owned_cells_on_rank_j > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_inlet_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_POS_X: // Inlet on the global I-maximum face (face of cell C_i=last_global_cell_idx_i)
            // Rank services if it owns the last cell in I-direction,
            // and has extent in J and K.
            if (last_global_cell_idx_i >= 0 && /* Check for valid global domain */
                (owned_start_cell_i + num_owned_cells_on_rank_i - 1) == last_global_cell_idx_i && /* Rank's last cell is the global last cell */
                num_owned_cells_on_rank_j > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_inlet_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_NEG_Y:
            if (info->ys == 0 && num_owned_cells_on_rank_j > 0 &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_inlet_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_POS_Y:
            if (last_global_cell_idx_j >= 0 &&
                (owned_start_cell_j + num_owned_cells_on_rank_j - 1) == last_global_cell_idx_j &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_inlet_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_NEG_Z:
            if (info->zs == 0 && num_owned_cells_on_rank_k > 0 &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_j > 0) {
                *can_service_inlet_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_POS_Z:
            if (last_global_cell_idx_k >= 0 &&
                (owned_start_cell_k + num_owned_cells_on_rank_k - 1) == last_global_cell_idx_k &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_j > 0) {
                *can_service_inlet_out = PETSC_TRUE;
            }
            break;
        default:
             LOG_ALLOW(LOCAL, LOG_WARNING, "[Rank %d]: Unknown inlet face %s.\n", rank_for_logging, BCFaceToString((BCFace)user->identifiedInletBCFace));
            break;
    }

      LOG_ALLOW(LOCAL, LOG_TRACE,
      "[Rank %d] Check Service for Inlet %s:\n"
      "    - Local Domain: starts at cell (%d,%d,%d), has (%d,%d,%d) cells.\n"
      "    - Global Domain: has (%d,%d,%d) nodes, so last cell is (%d,%d,%d).\n",
      rank_for_logging,
      BCFaceToString((BCFace)user->identifiedInletBCFace),
      owned_start_cell_i, owned_start_cell_j, owned_start_cell_k,
      num_owned_cells_on_rank_i, num_owned_cells_on_rank_j, num_owned_cells_on_rank_k,
      IM_nodes_global, JM_nodes_global, KM_nodes_global,
      last_global_cell_idx_i, last_global_cell_idx_j, last_global_cell_idx_k);

      LOG_ALLOW(LOCAL, LOG_DEBUG,"[Rank %d] Inlet Face %s Service Check Result: %s\n",
                rank_for_logging,
                BCFaceToString((BCFace)user->identifiedInletBCFace),
                (*can_service_inlet_out) ? "CAN SERVICE" : "CANNOT SERVICE");

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CanRankServiceFace"

/**
 * @brief Determines if the current MPI rank owns any part of a specified global face.
 *
 * This function is a general utility for parallel boundary operations. It checks if the
 * local domain of the current MPI rank is adjacent to a specified global boundary face.
 * A rank "services" a face if it owns the cells adjacent to that global boundary and has
 * a non-zero extent (i.e., owns at least one cell) in the tangential dimensions of that face.
 *
 * @param info              Pointer to the DMDALocalInfo for the current rank's DA.
 * @param IM_nodes_global Global number of nodes in the I-direction (e.g., user->IM + 1 if user->IM is cell count).
 * @param JM_nodes_global Global number of nodes in the J-direction.
 * @param KM_nodes_global Global number of nodes in the K-direction.
 * @param face_id           The specific global face (e.g., BC_FACE_NEG_Z) to check.
 * @param[out] can_service_out Pointer to a PetscBool; set to PETSC_TRUE if the rank
 *                           services the face, PETSC_FALSE otherwise.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode CanRankServiceFace(const DMDALocalInfo *info, PetscInt IM_nodes_global, PetscInt JM_nodes_global, PetscInt KM_nodes_global,
                                  BCFace face_id, PetscBool *can_service_out)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank_for_logging;
    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank_for_logging); CHKERRQ(ierr);

    *can_service_out = PETSC_FALSE; // Default to no service

    // Get the range of cells owned by this rank
    PetscInt owned_start_cell_i, num_owned_cells_on_rank_i;
    PetscInt owned_start_cell_j, num_owned_cells_on_rank_j;
    PetscInt owned_start_cell_k, num_owned_cells_on_rank_k;
    ierr = GetOwnedCellRange(info, 0, &owned_start_cell_i, &num_owned_cells_on_rank_i); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 1, &owned_start_cell_j, &num_owned_cells_on_rank_j); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 2, &owned_start_cell_k, &num_owned_cells_on_rank_k); CHKERRQ(ierr);

    // Determine the global index of the last cell (0-indexed) in each direction.
    PetscInt last_global_cell_idx_i = (IM_nodes_global > 1) ? (IM_nodes_global - 2) : -1;
    PetscInt last_global_cell_idx_j = (JM_nodes_global > 1) ? (JM_nodes_global - 2) : -1;
    PetscInt last_global_cell_idx_k = (KM_nodes_global > 1) ? (KM_nodes_global - 2) : -1;

    switch (face_id) {
        case BC_FACE_NEG_X:
            if (info->xs == 0 && num_owned_cells_on_rank_i > 0 &&
                num_owned_cells_on_rank_j > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_POS_X:
            if (last_global_cell_idx_i >= 0 &&
                (owned_start_cell_i + num_owned_cells_on_rank_i - 1) == last_global_cell_idx_i &&
                num_owned_cells_on_rank_j > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_NEG_Y:
            if (info->ys == 0 && num_owned_cells_on_rank_j > 0 &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_POS_Y:
            if (last_global_cell_idx_j >= 0 &&
                (owned_start_cell_j + num_owned_cells_on_rank_j - 1) == last_global_cell_idx_j &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_k > 0) {
                *can_service_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_NEG_Z:
            if (info->zs == 0 && num_owned_cells_on_rank_k > 0 &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_j > 0) {
                *can_service_out = PETSC_TRUE;
            }
            break;
        case BC_FACE_POS_Z:
            if (last_global_cell_idx_k >= 0 &&
                (owned_start_cell_k + num_owned_cells_on_rank_k - 1) == last_global_cell_idx_k &&
                num_owned_cells_on_rank_i > 0 && num_owned_cells_on_rank_j > 0) {
                *can_service_out = PETSC_TRUE;
            }
            break;
        default:
             LOG_ALLOW(LOCAL, LOG_WARNING, "Rank %d: Unknown face enum %d. \n", rank_for_logging, face_id);
            break;
    }

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d check for face %s: Result=%s. \n",
        rank_for_logging, BCFaceToString((BCFace)face_id), (*can_service_out ? "TRUE" : "FALSE"));

    PROFILE_FUNCTION_END;
        
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GetDeterministicFaceGridLocation"

/**
 * @brief Places particles in a deterministic grid/raster pattern on a specified domain face.
 *
 * This function creates a set of equidistant, parallel lines of particles near the four
 * edges of the face specified by user->identifiedInletBCFace. The number of lines drawn
 * from each edge is hardcoded within this function (default is 2).
 *
 * For example, if grid_layers=2 on face BC_FACE_NEG_X, the function will create particle lines at:
 * - y ~ 0*dy, y ~ 1*dy (parallel to the Z-axis, starting from the J=0 edge)
 * - y ~ y_max, y ~ y_max-dy (parallel to the Z-axis, starting from the J=max edge)
 * - z ~ 0*dz, z ~ 1*dz (parallel to the Y-axis, starting from the K=0 edge)
 * - z ~ z_max, z ~ z_max-dz (parallel to the Y-axis, starting from the K=max edge)
 *
 * The particle's final position is set just inside the target cell face to ensure it is
 * correctly located. The total number of particles (simCtx->np) is distributed as evenly
 * as possible among all generated lines.
 *
 * The function includes extensive validation to stop with an error if the requested grid
 * placement is geometrically impossible (e.g., in a 2D domain or if layers would overlap).
 * It also issues warnings for non-fatal but potentially unintended configurations.
 *
 * @param user Pointer to UserCtx, which must contain a valid identifiedInletBCFace.
 * @param info Pointer to DMDALocalInfo for the current rank's grid layout.
 * @param xs_gnode_rank, ys_gnode_rank, zs_gnode_rank Local starting node indices (incl. ghosts) for the rank's DA.
 * @param IM_cells_global, JM_cells_global, KM_cells_global Global cell counts.
 * @param particle_global_id The unique global ID of the particle being placed (from 0 to np-1).
 * @param[out] ci_metric_lnode_out Local I-node index of the selected cell's origin.
 * @param[out] cj_metric_lnode_out Local J-node index of the selected cell's origin.
 * @param[out] ck_metric_lnode_out Local K-node index of the selected cell's origin.
 * @param[out] xi_metric_logic_out Logical xi-coordinate [0,1] within the cell.
 * @param[out] eta_metric_logic_out Logical eta-coordinate [0,1] within the cell.
 * @param[out] zta_metric_logic_out Logical zta-coordinate [0,1] within the cell.
 * @param[out] placement_successful_out PETSC_TRUE if the point belongs to this rank, PETSC_FALSE otherwise.
 * @return PetscErrorCode
 */
PetscErrorCode GetDeterministicFaceGridLocation(
    UserCtx *user, const DMDALocalInfo *info,
    PetscInt xs_gnode_rank, PetscInt ys_gnode_rank, PetscInt zs_gnode_rank,
    PetscInt IM_cells_global, PetscInt JM_cells_global, PetscInt KM_cells_global,
    PetscInt64 particle_global_id,
    PetscInt *ci_metric_lnode_out, PetscInt *cj_metric_lnode_out, PetscInt *ck_metric_lnode_out,
    PetscReal *xi_metric_logic_out, PetscReal *eta_metric_logic_out, PetscReal *zta_metric_logic_out,
    PetscBool *placement_successful_out)
{
    SimCtx *simCtx = user->simCtx;
    PetscReal global_logic_i, global_logic_j, global_logic_k;
    PetscErrorCode ierr;
    PetscMPIInt rank_for_logging;

    PetscFunctionBeginUser;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank_for_logging); CHKERRQ(ierr);

    *placement_successful_out = PETSC_FALSE; // Default to failure

    // --- Step 1: Configuration and Input Validation ---

    // *** Hardcoded number of grid layers. Change this value to alter the pattern. ***
    const PetscInt grid_layers = 2;

    LOG_ALLOW(LOCAL, LOG_DEBUG,
        "[Rank %d] Placing particle %lld on face %s with grid_layers=%d in global domain (%d,%d,%d) cells.\n",
        rank_for_logging, (long long)particle_global_id, BCFaceToString(user->identifiedInletBCFace), grid_layers,
        IM_cells_global, JM_cells_global, KM_cells_global);

    const char *face_name = BCFaceToString(user->identifiedInletBCFace);

    // Fatal Error Checks: Ensure the requested grid is geometrically possible.
    // The total layers from opposite faces (2 * grid_layers) must be less than the domain size.
    switch (user->identifiedInletBCFace) {
        case BC_FACE_NEG_X: case BC_FACE_POS_X:
            if (JM_cells_global <= 1 || KM_cells_global <= 1) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Cannot place grid on face %s for a 2D/1D domain (J-cells=%d, K-cells=%d).", face_name, JM_cells_global, KM_cells_global);
            if (2 * grid_layers >= JM_cells_global || 2 * grid_layers >= KM_cells_global) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Grid layers (%d) from opposing J/K faces would overlap in this domain (J-cells=%d, K-cells=%d).", grid_layers, JM_cells_global, KM_cells_global);
            break;
        case BC_FACE_NEG_Y: case BC_FACE_POS_Y:
            if (IM_cells_global <= 1 || KM_cells_global <= 1) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Cannot place grid on face %s for a 2D/1D domain (I-cells=%d, K-cells=%d).", face_name, IM_cells_global, KM_cells_global);
            if (2 * grid_layers >= IM_cells_global || 2 * grid_layers >= KM_cells_global) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Grid layers (%d) from opposing I/K faces would overlap in this domain (I-cells=%d, K-cells=%d).", grid_layers, IM_cells_global, KM_cells_global);
            break;
        case BC_FACE_NEG_Z: case BC_FACE_POS_Z:
            if (IM_cells_global <= 1 || JM_cells_global <= 1) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Cannot place grid on face %s for a 2D/1D domain (I-cells=%d, J-cells=%d).", face_name, IM_cells_global, JM_cells_global);
            if (2 * grid_layers >= IM_cells_global || 2 * grid_layers >= JM_cells_global) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Grid layers (%d) from opposing I/J faces would overlap in this domain (I-cells=%d, J-cells=%d).", grid_layers, IM_cells_global, JM_cells_global);
            break;
        default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid identifiedInletBCFace specified: %d", user->identifiedInletBCFace);
    }

    const PetscInt num_lines_total = 4 * grid_layers;
    if (simCtx->np < num_lines_total) {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "Warning: Total particle count (%lld) is less than the number of grid lines requested (%d). Some lines may be empty.\n", (long long)simCtx->np, num_lines_total);
    }
    if (simCtx->np > 0 && simCtx->np % num_lines_total != 0) {
        LOG_ALLOW(GLOBAL, LOG_WARNING, "Warning: Total particle count (%lld) is not evenly divisible by the number of grid lines (%d). Distribution will be uneven.\n", (long long)simCtx->np, num_lines_total);
    }

    // --- Step 2: Map global particle ID to a line and a point on that line ---
    if (simCtx->np == 0) PetscFunctionReturn(0); // Nothing to do

    LOG_ALLOW(LOCAL, LOG_TRACE, "[Rank %d] Distributing %lld particles over %d lines on face %s.\n",
        rank_for_logging, (long long)simCtx->np, num_lines_total, face_name);

    const PetscInt points_per_line = PetscMax(1, simCtx->np / num_lines_total);
    PetscInt line_index = particle_global_id / points_per_line;
    PetscInt point_index_on_line = particle_global_id % points_per_line;
    line_index = PetscMin(line_index, num_lines_total - 1); // Clamp to handle uneven division

    // Decode the line_index into an edge group (0-3) and a layer within that group (0 to grid_layers-1)
    const PetscInt edge_group = line_index / grid_layers;
    const PetscInt layer_index = line_index % grid_layers;

    // --- Step 3: Calculate placement coordinates based on the decoded indices ---
    const PetscReal epsilon = 1.0e-6; // Small offset to keep particles off exact cell boundaries
    const PetscReal layer_spacing_norm_i = (IM_cells_global > 0) ? 1.0 / (PetscReal)IM_cells_global : 0.0;
    const PetscReal layer_spacing_norm_j = (JM_cells_global > 0) ? 1.0 / (PetscReal)JM_cells_global : 0.0;
    const PetscReal layer_spacing_norm_k = (KM_cells_global > 0) ? 1.0 / (PetscReal)KM_cells_global : 0.0;

    PetscReal variable_coord; // The coordinate that varies along a line
    if (points_per_line <= 1) {
        variable_coord = 0.5; // Place single point in the middle
    } else {
        variable_coord = ((PetscReal)point_index_on_line + 0.5)/ (PetscReal)(points_per_line);
    }
    variable_coord = PetscMin(1.0 - epsilon, PetscMax(epsilon, variable_coord)); // Clamp within [eps, 1-eps]

    // Main logic switch to determine the three global logical coordinates
    switch (user->identifiedInletBCFace) {
        case BC_FACE_NEG_X:
            global_logic_i = 0.5 * layer_spacing_norm_i; // Place near the face, in the middle of the first cell
            if (edge_group == 0)      { global_logic_j = (PetscReal)layer_index * layer_spacing_norm_j + epsilon; global_logic_k = variable_coord; }
            else if (edge_group == 1) { global_logic_j = 1.0 - ((PetscReal)layer_index * layer_spacing_norm_j) - epsilon; global_logic_k = variable_coord; }
            else if (edge_group == 2) { global_logic_k = (PetscReal)layer_index * layer_spacing_norm_k + epsilon; global_logic_j = variable_coord; }
            else /* edge_group == 3 */ { global_logic_k = 1.0 - ((PetscReal)layer_index * layer_spacing_norm_k) - epsilon; global_logic_j = variable_coord; }
            break;
        case BC_FACE_POS_X:
            global_logic_i = 1.0 - (0.5 * layer_spacing_norm_i); // Place near the face, in the middle of the last cell
            if (edge_group == 0)      { global_logic_j = (PetscReal)layer_index * layer_spacing_norm_j + epsilon; global_logic_k = variable_coord; }
            else if (edge_group == 1) { global_logic_j = 1.0 - ((PetscReal)layer_index * layer_spacing_norm_j) - epsilon; global_logic_k = variable_coord; }
            else if (edge_group == 2) { global_logic_k = (PetscReal)layer_index * layer_spacing_norm_k + epsilon; global_logic_j = variable_coord; }
            else /* edge_group == 3 */ { global_logic_k = 1.0 - ((PetscReal)layer_index * layer_spacing_norm_k) - epsilon; global_logic_j = variable_coord; }
            break;
        case BC_FACE_NEG_Y:
            global_logic_j = 0.5 * layer_spacing_norm_j;
            if (edge_group == 0)      { global_logic_i = (PetscReal)layer_index * layer_spacing_norm_i + epsilon; global_logic_k = variable_coord; }
            else if (edge_group == 1) { global_logic_i = 1.0 - ((PetscReal)layer_index * layer_spacing_norm_i) - epsilon; global_logic_k = variable_coord; }
            else if (edge_group == 2) { global_logic_k = (PetscReal)layer_index * layer_spacing_norm_k + epsilon; global_logic_i = variable_coord; }
            else /* edge_group == 3 */ { global_logic_k = 1.0 - ((PetscReal)layer_index * layer_spacing_norm_k) - epsilon; global_logic_i = variable_coord; }
            break;
        case BC_FACE_POS_Y:
            global_logic_j = 1.0 - (0.5 * layer_spacing_norm_j);
            if (edge_group == 0)      { global_logic_i = (PetscReal)layer_index * layer_spacing_norm_i + epsilon; global_logic_k = variable_coord; }
            else if (edge_group == 1) { global_logic_i = 1.0 - ((PetscReal)layer_index * layer_spacing_norm_i) - epsilon; global_logic_k = variable_coord; }
            else if (edge_group == 2) { global_logic_k = (PetscReal)layer_index * layer_spacing_norm_k + epsilon; global_logic_i = variable_coord; }
            else /* edge_group == 3 */ { global_logic_k = 1.0 - ((PetscReal)layer_index * layer_spacing_norm_k) - epsilon; global_logic_i = variable_coord; }
            break;
        case BC_FACE_NEG_Z:
            global_logic_k = 0.5 * layer_spacing_norm_k;
            if (edge_group == 0)      { global_logic_i = (PetscReal)layer_index * layer_spacing_norm_i + epsilon; global_logic_j = variable_coord; }
            else if (edge_group == 1) { global_logic_i = 1.0 - ((PetscReal)layer_index * layer_spacing_norm_i) - epsilon; global_logic_j = variable_coord; }
            else if (edge_group == 2) { global_logic_j = (PetscReal)layer_index * layer_spacing_norm_j + epsilon; global_logic_i = variable_coord; }
            else /* edge_group == 3 */ { global_logic_j = 1.0 - ((PetscReal)layer_index * layer_spacing_norm_j) - epsilon; global_logic_i = variable_coord; }
            break;
        case BC_FACE_POS_Z:
            global_logic_k = 1.0 - (0.5 * layer_spacing_norm_k);
            if (edge_group == 0)      { global_logic_i = (PetscReal)layer_index * layer_spacing_norm_i + epsilon; global_logic_j = variable_coord; }
            else if (edge_group == 1) { global_logic_i = 1.0 - ((PetscReal)layer_index * layer_spacing_norm_i) - epsilon; global_logic_j = variable_coord; }
            else if (edge_group == 2) { global_logic_j = (PetscReal)layer_index * layer_spacing_norm_j + epsilon; global_logic_i = variable_coord; }
            else /* edge_group == 3 */ { global_logic_j = 1.0 - ((PetscReal)layer_index * layer_spacing_norm_j) - epsilon; global_logic_i = variable_coord; }
            break;
    }

    LOG_ALLOW(LOCAL, LOG_TRACE,
        "[Rank %d] Particle %lld assigned to line %d (edge group %d, layer %d) with variable_coord=%.4f.\n"
        "    -> Global logical coords: (i,j,k) = (%.6f, %.6f, %.6f)\n",
        rank_for_logging, (long long)particle_global_id, line_index, edge_group, layer_index, variable_coord,
        global_logic_i, global_logic_j, global_logic_k);

    // --- Step 4: Convert global logical coordinate to global cell index and intra-cell logicals ---
    PetscReal global_cell_coord_i = global_logic_i * IM_cells_global;
    PetscInt  I_g = (PetscInt)global_cell_coord_i;
    *xi_metric_logic_out = global_cell_coord_i - I_g;

    PetscReal global_cell_coord_j = global_logic_j * JM_cells_global;
    PetscInt  J_g = (PetscInt)global_cell_coord_j;
    *eta_metric_logic_out = global_cell_coord_j - J_g;

    PetscReal global_cell_coord_k = global_logic_k * KM_cells_global;
    PetscInt  K_g = (PetscInt)global_cell_coord_k;
    *zta_metric_logic_out = global_cell_coord_k - K_g;

    // --- Step 5: Check if this rank owns the target cell and finalize outputs ---
    if ((I_g >= info->xs && I_g < info->xs + info->xm) &&
        (J_g >= info->ys && J_g < info->ys + info->ym) &&
        (K_g >= info->zs && K_g < info->zs + info->zm))
    {
        // Convert global cell index to the local node index for this rank's DA patch
        *ci_metric_lnode_out = (I_g - info->xs) + xs_gnode_rank;
        *cj_metric_lnode_out = (J_g - info->ys) + ys_gnode_rank;
        *ck_metric_lnode_out = (K_g - info->zs) + zs_gnode_rank;
        *placement_successful_out = PETSC_TRUE;
    }

    LOG_ALLOW(LOCAL, LOG_VERBOSE,
        "[Rank %d] Particle %lld placement %s. Local cell origin node: (I,J,K) = (%d,%d,%d), intra-cell logicals: (xi,eta,zta)=(%.6f,%.6f,%.6f)\n",
        rank_for_logging, (long long)particle_global_id,
        (*placement_successful_out ? "SUCCESSFUL" : "NOT ON THIS RANK"),
        *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out,
        *xi_metric_logic_out, *eta_metric_logic_out, *zta_metric_logic_out);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GetRandomFCellAndLogicOnInletFace"

/**
 * @brief Assuming the current rank services the inlet face, this function selects a random
 *        cell (owned by this rank on that face) and random logical coordinates within that cell,
 *        suitable for placing a particle on the inlet surface.
 *
 * It is the caller's responsibility to ensure CanRankServiceInletFace returned true.
 *
 * @param user Pointer to UserCtx.
 * @param info Pointer to DMDALocalInfo for the current rank (node-based).
 * @param xs_gnode, ys_gnode, zs_gnode Local starting node indices (incl. ghosts) for the rank's DA.
 * @param IM_nodes_global, JM_nodes_global, KM_nodes_global Global node counts.
 * @param rand_logic_i_ptr, rand_logic_j_ptr, rand_logic_k_ptr Pointers to RNGs for logical coords.
 * @param[out] ci_metric_lnode_out, cj_metric_lnode_out, ck_metric_lnode_out Local node indices of the selected cell's origin (these are local to the rank's DA including ghosts).
 * @param[out] xi_metric_logic_out, eta_metric_logic_out, zta_metric_logic_out Logical coords [0,1] within the cell.
 * @return PetscErrorCode
 */
PetscErrorCode GetRandomCellAndLogicalCoordsOnInletFace(
    UserCtx *user, const DMDALocalInfo *info,
    PetscInt xs_gnode_rank, PetscInt ys_gnode_rank, PetscInt zs_gnode_rank, // Local starting node index (with ghosts) of the rank's DA patch
    PetscInt IM_nodes_global, PetscInt JM_nodes_global, PetscInt KM_nodes_global,
    PetscRandom *rand_logic_i_ptr, PetscRandom *rand_logic_j_ptr, PetscRandom *rand_logic_k_ptr,
    PetscInt *ci_metric_lnode_out, PetscInt *cj_metric_lnode_out, PetscInt *ck_metric_lnode_out,
    PetscReal *xi_metric_logic_out, PetscReal *eta_metric_logic_out, PetscReal *zta_metric_logic_out)
{
    PetscErrorCode ierr = 0;
    PetscReal r_val_i_sel, r_val_j_sel, r_val_k_sel;
    PetscInt local_cell_idx_on_face_dim1 = 0; // 0-indexed relative to owned cells on face
    PetscInt local_cell_idx_on_face_dim2 = 0;
    PetscMPIInt rank_for_logging; 

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank_for_logging); CHKERRQ(ierr); 

    // Defaults for cell origin node (local index for the rank's DA patch, including ghosts)
    *ci_metric_lnode_out = xs_gnode_rank; *cj_metric_lnode_out = ys_gnode_rank; *ck_metric_lnode_out = zs_gnode_rank;
    // Defaults for logical coordinates
    *xi_metric_logic_out = 0.5; *eta_metric_logic_out = 0.5; *zta_metric_logic_out = 0.5;

    // Get number of cells this rank owns in each dimension (tangential to the face mainly)
    PetscInt owned_start_cell_i, num_owned_cells_on_rank_i;
    PetscInt owned_start_cell_j, num_owned_cells_on_rank_j;
    PetscInt owned_start_cell_k, num_owned_cells_on_rank_k;

    ierr = GetOwnedCellRange(info, 0, &owned_start_cell_i, &num_owned_cells_on_rank_i); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 1, &owned_start_cell_j, &num_owned_cells_on_rank_j); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 2, &owned_start_cell_k, &num_owned_cells_on_rank_k); CHKERRQ(ierr);

    // Index of the last cell (0-indexed) in each global direction
    PetscInt last_global_cell_idx_i = (IM_nodes_global > 1) ? (IM_nodes_global - 2) : -1;
    PetscInt last_global_cell_idx_j = (JM_nodes_global > 1) ? (JM_nodes_global - 2) : -1;
    PetscInt last_global_cell_idx_k = (KM_nodes_global > 1) ? (KM_nodes_global - 2) : -1;
    
    LOG_ALLOW(LOCAL, LOG_TRACE, "Rank %d: Inlet face %s. Owned cells(i,j,k):(%d,%d,%d). GlobNodes(I,J,K):(%d,%d,%d). Rank's DA node starts at (%d,%d,%d).\n",
        rank_for_logging, user->identifiedInletBCFace, num_owned_cells_on_rank_i,num_owned_cells_on_rank_j,num_owned_cells_on_rank_k,
        IM_nodes_global,JM_nodes_global,KM_nodes_global, xs_gnode_rank,ys_gnode_rank,zs_gnode_rank);


    switch (user->identifiedInletBCFace) {
        case BC_FACE_NEG_X: // Particle on -X face of cell C_0 (origin node N_0)
            // Cell origin node is the first owned node in I by this rank (global index info->xs).
            // Its local index within the rank's DA (incl ghosts) is xs_gnode_rank.
            *ci_metric_lnode_out = xs_gnode_rank;
            *xi_metric_logic_out = 1.0e-6;

            // Tangential dimensions are J and K. Select an owned cell randomly on this face.
            // num_owned_cells_on_rank_j/k must be > 0 (checked by CanRankServiceInletFace)
            ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, &r_val_j_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim1 = (PetscInt)(r_val_j_sel * num_owned_cells_on_rank_j); // Index among owned J-cells
            local_cell_idx_on_face_dim1 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim1), num_owned_cells_on_rank_j - 1);
            *cj_metric_lnode_out = ys_gnode_rank + local_cell_idx_on_face_dim1; // Offset from start of rank's J-nodes

            ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, &r_val_k_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim2 = (PetscInt)(r_val_k_sel * num_owned_cells_on_rank_k);
            local_cell_idx_on_face_dim2 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim2), num_owned_cells_on_rank_k - 1);
            *ck_metric_lnode_out = zs_gnode_rank + local_cell_idx_on_face_dim2;

            ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, eta_metric_logic_out); CHKERRQ(ierr);
            ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, zta_metric_logic_out); CHKERRQ(ierr);
            break;

        case BC_FACE_POS_X: // Particle on +X face of cell C_last_I (origin node N_last_I_origin)
            // Origin node of the last I-cell is global_node_idx = last_global_cell_idx_i.
            // Its local index in rank's DA: (last_global_cell_idx_i - info->xs) + xs_gnode_rank
            *ci_metric_lnode_out = xs_gnode_rank + (last_global_cell_idx_i - info->xs);
            *xi_metric_logic_out = 1.0 - 1.0e-6;

            ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, &r_val_j_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim1 = (PetscInt)(r_val_j_sel * num_owned_cells_on_rank_j);
            local_cell_idx_on_face_dim1 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim1), num_owned_cells_on_rank_j - 1);
            *cj_metric_lnode_out = ys_gnode_rank + local_cell_idx_on_face_dim1;

            ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, &r_val_k_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim2 = (PetscInt)(r_val_k_sel * num_owned_cells_on_rank_k);
            local_cell_idx_on_face_dim2 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim2), num_owned_cells_on_rank_k - 1);
            *ck_metric_lnode_out = zs_gnode_rank + local_cell_idx_on_face_dim2;

            ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, eta_metric_logic_out); CHKERRQ(ierr);
            ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, zta_metric_logic_out); CHKERRQ(ierr);
            break;
        // ... (Cases for Y and Z faces, following the same pattern) ...
        case BC_FACE_NEG_Y:
            *cj_metric_lnode_out = ys_gnode_rank;
            *eta_metric_logic_out = 1.0e-6;
            ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val_i_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim1 = (PetscInt)(r_val_i_sel * num_owned_cells_on_rank_i);
            local_cell_idx_on_face_dim1 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim1), num_owned_cells_on_rank_i - 1);
            *ci_metric_lnode_out = xs_gnode_rank + local_cell_idx_on_face_dim1;
            ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, &r_val_k_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim2 = (PetscInt)(r_val_k_sel * num_owned_cells_on_rank_k);
            local_cell_idx_on_face_dim2 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim2), num_owned_cells_on_rank_k - 1);
            *ck_metric_lnode_out = zs_gnode_rank + local_cell_idx_on_face_dim2;
            ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, xi_metric_logic_out); CHKERRQ(ierr);
            ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, zta_metric_logic_out); CHKERRQ(ierr);
            break;
        case BC_FACE_POS_Y:
            *cj_metric_lnode_out = ys_gnode_rank + (last_global_cell_idx_j - info->ys);
            *eta_metric_logic_out = 1.0 - 1.0e-6;
            ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val_i_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim1 = (PetscInt)(r_val_i_sel * num_owned_cells_on_rank_i);
            local_cell_idx_on_face_dim1 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim1), num_owned_cells_on_rank_i - 1);
            *ci_metric_lnode_out = xs_gnode_rank + local_cell_idx_on_face_dim1;
            ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, &r_val_k_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim2 = (PetscInt)(r_val_k_sel * num_owned_cells_on_rank_k);
            local_cell_idx_on_face_dim2 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim2), num_owned_cells_on_rank_k - 1);
            *ck_metric_lnode_out = zs_gnode_rank + local_cell_idx_on_face_dim2;
            ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, xi_metric_logic_out); CHKERRQ(ierr);
            ierr = PetscRandomGetValueReal(*rand_logic_k_ptr, zta_metric_logic_out); CHKERRQ(ierr);
            break;
        case BC_FACE_NEG_Z: // Your example case
            *ck_metric_lnode_out = zs_gnode_rank; // Cell origin is the first owned node in K by this rank
            *zta_metric_logic_out = 1.0e-6;      // Place particle slightly inside this cell from its -Z face
            // Tangential dimensions are I and J
            ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val_i_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim1 = (PetscInt)(r_val_i_sel * num_owned_cells_on_rank_i);
            local_cell_idx_on_face_dim1 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim1), num_owned_cells_on_rank_i - 1);
            *ci_metric_lnode_out = xs_gnode_rank + local_cell_idx_on_face_dim1;

            ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, &r_val_j_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim2 = (PetscInt)(r_val_j_sel * num_owned_cells_on_rank_j);
            local_cell_idx_on_face_dim2 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim2), num_owned_cells_on_rank_j - 1);
            *cj_metric_lnode_out = ys_gnode_rank + local_cell_idx_on_face_dim2;

            ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, xi_metric_logic_out); CHKERRQ(ierr); // Intra-cell logical for I
            ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, eta_metric_logic_out); CHKERRQ(ierr); // Intra-cell logical for J
            break;
        case BC_FACE_POS_Z:
            *ck_metric_lnode_out = zs_gnode_rank + (last_global_cell_idx_k - info->zs);
            *zta_metric_logic_out = 1.0 - 1.0e-6;
            ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, &r_val_i_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim1 = (PetscInt)(r_val_i_sel * num_owned_cells_on_rank_i);
            local_cell_idx_on_face_dim1 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim1), num_owned_cells_on_rank_i - 1);
            *ci_metric_lnode_out = xs_gnode_rank + local_cell_idx_on_face_dim1;
            ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, &r_val_j_sel); CHKERRQ(ierr);
            local_cell_idx_on_face_dim2 = (PetscInt)(r_val_j_sel * num_owned_cells_on_rank_j);
            local_cell_idx_on_face_dim2 = PetscMin(PetscMax(0, local_cell_idx_on_face_dim2), num_owned_cells_on_rank_j - 1);
            *cj_metric_lnode_out = ys_gnode_rank + local_cell_idx_on_face_dim2;
            ierr = PetscRandomGetValueReal(*rand_logic_i_ptr, xi_metric_logic_out); CHKERRQ(ierr);
            ierr = PetscRandomGetValueReal(*rand_logic_j_ptr, eta_metric_logic_out); CHKERRQ(ierr);
            break;
        default:
             SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "GetRandomCellAndLogicOnInletFace: Invalid user->identifiedInletBCFace %d. \n", user->identifiedInletBCFace);
    }

    PetscReal eps = 1.0e-7;
    if (user->identifiedInletBCFace == BC_FACE_NEG_X || user->identifiedInletBCFace == BC_FACE_POS_X) {
        *eta_metric_logic_out = PetscMin(PetscMax(0.0, *eta_metric_logic_out), 1.0 - eps);
        *zta_metric_logic_out = PetscMin(PetscMax(0.0, *zta_metric_logic_out), 1.0 - eps);
    } else if (user->identifiedInletBCFace == BC_FACE_NEG_Y || user->identifiedInletBCFace == BC_FACE_POS_Y) {
        *xi_metric_logic_out  = PetscMin(PetscMax(0.0, *xi_metric_logic_out),  1.0 - eps);
        *zta_metric_logic_out = PetscMin(PetscMax(0.0, *zta_metric_logic_out), 1.0 - eps);
    } else { 
        *xi_metric_logic_out  = PetscMin(PetscMax(0.0, *xi_metric_logic_out),  1.0 - eps);
        *eta_metric_logic_out = PetscMin(PetscMax(0.0, *eta_metric_logic_out), 1.0 - eps);
    }
    
    LOG_ALLOW(LOCAL, LOG_VERBOSE, "Rank %d: Target Cell Node =(%d,%d,%d). (xi,et,zt)=(%.2e,%.2f,%.2f). \n",
        rank_for_logging, *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out,
        *xi_metric_logic_out, *eta_metric_logic_out, *zta_metric_logic_out);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}



#undef __FUNCT__
#define __FUNCT__ "TranslateModernBCsToLegacy"

PetscErrorCode TranslateModernBCsToLegacy(UserCtx *user)
{
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(GLOBAL,LOG_DEBUG," Translating modern BC config to legacy integer codes...\n");

    for (int i = 0; i < 6; i++) {
      BCType modern_type = user->boundary_faces[i].mathematical_type;
      user->bctype[i] = (int)modern_type;

      
      BCFace current_face = (BCFace)i;
      const char* face_str  = BCFaceToString(current_face);
      const char* bc_type_str = BCTypeToString(modern_type);
      LOG_ALLOW(GLOBAL,LOG_TRACE," for face %s(%d), legacy type = %d & modern type = %s .\n",face_str,i,user->bctype[i],bc_type_str);
    }
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"    -> Translation complete.\n");
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "BoundaryCondition_Create"
/**
 * @brief (Private) Creates and configures a specific BoundaryCondition handler object.
 *
 * This function acts as a factory. Based on the requested handler_type, it allocates
 * a BoundaryCondition object and populates it with the correct set of function
 * pointers corresponding to that specific behavior.
 *
 * @param handler_type The specific handler to create (e.g., BC_HANDLER_WALL_NOSLIP).
 * @param[out] new_bc_ptr  A pointer to where the newly created BoundaryCondition
 *                         object's address will be stored.
 * @return PetscErrorCode 0 on success.
 */

PetscErrorCode BoundaryCondition_Create(BCHandlerType handler_type, BoundaryCondition **new_bc_ptr)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    const char* handler_name = BCHandlerTypeToString(handler_type);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Factory called for handler type %s. \n", handler_name);

    ierr = PetscMalloc1(1, new_bc_ptr); CHKERRQ(ierr);
    BoundaryCondition *bc = *new_bc_ptr;

    bc->type        = handler_type;
    bc->priority    = -1;  // Default priority; can be overridden in specific handlers
    bc->data        = NULL;
    bc->Initialize  = NULL;
    bc->PreStep     = NULL;
    bc->Apply       = NULL;
    bc->PostStep    = NULL;
    bc->Destroy     = NULL;

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Allocated generic handler object at address %p.\n", (void*)bc);

    switch (handler_type) {

        case BC_HANDLER_OUTLET_CONSERVATION:
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Dispatching to Create_OutletConservation().\n");
            ierr = Create_OutletConservation(bc); CHKERRQ(ierr);
            break;

        case BC_HANDLER_WALL_NOSLIP:
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Dispatching to Create_WallNoSlip().\n");
            ierr = Create_WallNoSlip(bc); CHKERRQ(ierr);
            break;

        case BC_HANDLER_INLET_CONSTANT_VELOCITY:
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Dispatching to Create_InletConstantVelocity().\n");
	          ierr = Create_InletConstantVelocity(bc); CHKERRQ(ierr);
            break;

        //case BC_HANDLER_INLET_PARABOLIC:
        //    LOG_ALLOW(LOCAL, LOG_DEBUG, "Dispatching to Create_InletParabolicProfile().\n");
	      //    ierr = Create_InletParabolicProfile(bc); CHKERRQ(ierr);
        //    break;
        //Add cases for other handlers here in future phases 
        
        default:
            LOG_ALLOW(GLOBAL, LOG_ERROR, "Handler type (%s) is not recognized or implemented in the factory.\n", handler_name);
            SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE, "Boundary handler type %d (%s) not recognized in factory.\n", handler_type, handler_name);
    }

    if(bc->priority < 0) {
        LOG_ALLOW(GLOBAL, LOG_ERROR, "Handler type %d (%s) did not set a valid priority during creation.\n", handler_type, handler_name);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE, "Boundary handler type %d (%s) did not set a valid priority during creation.\n", handler_type, handler_name);
    }
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Successfully created and configured handler for %s.\n", handler_name);
    PetscFunctionReturn(0);
}

//================================================================================
//
//                       PUBLIC MASTER SETUP FUNCTION
//
//================================================================================
#undef __FUNCT__
#define __FUNCT__ "BoundarySystem_Initialize"
/**
 * @brief Initializes the entire boundary system based on a configuration file.
 */
PetscErrorCode BoundarySystem_Initialize(UserCtx *user, const char *bcs_filename)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting creation and initialization of all boundary handlers.\n");

    // =========================================================================
    // Step 0: Clear any existing boundary handlers (if re-initializing).
    // This ensures no memory leaks if this function is called multiple times.
    // =========================================================================
    for (int i = 0; i < 6; i++) {
        BoundaryFaceConfig *face_cfg = &user->boundary_faces[i];
        if (face_cfg->handler) {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Destroying existing handler on Face %s before re-initialization.\n", BCFaceToString((BCFace)i));
            if (face_cfg->handler->Destroy) {
                ierr = face_cfg->handler->Destroy(face_cfg->handler); CHKERRQ(ierr);
            }
            ierr = PetscFree(face_cfg->handler); CHKERRQ(ierr);
            face_cfg->handler = NULL;
        }
    }
    // =========================================================================

    // Step 0.1: Initiate flux sums to zero
    user->simCtx->FluxInSum = 0.0;
    user->simCtx->FluxOutSum = 0.0;
    user->simCtx->FarFluxInSum = 0.0;
    user->simCtx->FarFluxOutSum = 0.0;
    // =========================================================================

    // Step 1: Parse the configuration file to determine user intent.
    // This function, defined in io.c, populates the configuration enums and parameter
    // lists within the user->boundary_faces array on all MPI ranks.
    ierr = ParseAllBoundaryConditions(user, bcs_filename); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Configuration file '%s' parsed successfully.\n", bcs_filename);

    // Step 2: Create and Initialize the handler object for each of the 6 faces.
    for (int i = 0; i < 6; i++) {
        BoundaryFaceConfig *face_cfg = &user->boundary_faces[i];
        
        const char *face_name = BCFaceToString(face_cfg->face_id);
        const char *type_name = BCTypeToString(face_cfg->mathematical_type);
        const char *handler_name = BCHandlerTypeToString(face_cfg->handler_type);

        LOG_ALLOW(LOCAL, LOG_DEBUG, "Creating handler for Face %s with Type %s and handler '%s'.\n", face_name, type_name,handler_name);

        // Use the private factory to construct the correct handler object based on the parsed type.
        // The factory returns a pointer to the new handler object, which we store in the config struct.
        ierr = BoundaryCondition_Create(face_cfg->handler_type, &face_cfg->handler); CHKERRQ(ierr);

        // Step 3: Call the specific Initialize() method for the newly created handler.
        // This allows the handler to perform its own setup, like reading parameters from the
        // face_cfg->params list and setting the initial field values on its face.
        if (face_cfg->handler && face_cfg->handler->Initialize) {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Calling Initialize() method for handler %s(%s) on Face %s.\n",type_name,handler_name,face_name);
            
            // Prepare the context needed by the Initialize() function.
            BCContext ctx = {
                .user = user,
                .face_id = face_cfg->face_id,
                .global_inflow_sum = &user->simCtx->FluxInSum,  // Global flux sums are not relevant during initialization.
                .global_outflow_sum = &user->simCtx->FluxOutSum,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            ierr = face_cfg->handler->Initialize(face_cfg->handler, &ctx); CHKERRQ(ierr);
        } else {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Handler %s(%s) for Face %s has no Initialize() method, skipping.\n", type_name,handler_name,face_name);
        }
    }
    // =========================================================================
    // NO SYNCHRONIZATION NEEDED HERE
    // =========================================================================
    // Initialize() only reads parameters and allocates memory.
    // It does NOT modify field values (Ucat, Ucont, Ubcs).
    // Field values are set by:
    //   1. Initial conditions (before this function)
    //   2. Apply() during timestepping (after this function)
    // The first call to ApplyBoundaryConditions() will handle synchronization.
    // =========================================================================

    // ====================================================================================
    // --- NEW: Step 4: Synchronize Vectors After Initialization ---
    // This is the CRITICAL fix. The Initialize() calls have modified local vector
    // arrays on some ranks but not others. We must now update the global vector state
    // and then update all local ghost regions to be consistent.
    // ====================================================================================
     
    //LOG_ALLOW(GLOBAL, LOG_DEBUG, "Committing global boundary initializations to local vectors.\n");

    // Commit changes from the global vectors (Ucat, Ucont) to the local vectors (lUcat, lUcont)
    // NOTE: The Apply functions modified Ucat and Ucont via GetArray, which works on the global
    // representation.
    /*    
    ierr = DMGlobalToLocalBegin(user->fda, user->Ucat, INSERT_VALUES, user->lUcat); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->fda, user->Ucat, INSERT_VALUES, user->lUcat); CHKERRQ(ierr);
    
    ierr = DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont); CHKERRQ(ierr);
    ierr = DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont); CHKERRQ(ierr);
    
     // Now, update all local vectors (including ghost cells) from the newly consistent global vectors

    ierr = DMLocalToGlobalBegin(user->fda, user->lUcat, INSERT_VALUES, user->Ucat); CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(user->fda, user->lUcat, INSERT_VALUES, user->Ucat); CHKERRQ(ierr);
    
    ierr = DMLocalToGlobalBegin(user->fda, user->lUcont, INSERT_VALUES, user->Ucont); CHKERRQ(ierr);
    ierr = DMLocalToGlobalEnd(user->fda, user->lUcont, INSERT_VALUES, user->Ucont); CHKERRQ(ierr);
    */
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "All boundary handlers created and initialized successfully.\n");
    PetscFunctionReturn(0);
}


//================================================================================
//
//                      PUBLIC MASTER TIME-STEP FUNCTION
//
//================================================================================

#undef __FUNCT__
#define __FUNCT__ "BoundarySystem_ExecuteStep"
/**
 * @brief Executes all boundary condition handlers in priority order.
 *
 * This function orchestrates the application of boundary conditions across all
 * faces using a priority-based system. Each priority group is executed atomically:
 * handlers at a given priority complete their PreStep, Apply, and PostStep phases,
 * with MPI communication between phases as needed. This ensures proper data flow
 * for boundary conditions that depend on results from other boundaries.
 *
 * Priority execution order (matches legacy):
 *   0 (BC_PRIORITY_INLET):    Inlets - Set inflow, measure flux
 *   1 (BC_PRIORITY_FARFIELD): Farfield - Bidirectional flow, measure flux  
 *   2 (BC_PRIORITY_WALL):     Walls/Symmetry - Set velocity/gradients
 *   3 (BC_PRIORITY_OUTLET):   Outlets - Apply conservation correction
 *
 * NOTE: This function is called INSIDE the ApplyBoundaryConditions iteration loop.
 * It does NOT handle:
 *   - Contra2Cart (done by caller)
 *   - UpdateDummyCells (done by caller)
 *   - DMGlobalToLocal syncs (done by caller)
 *
 * @param user The UserCtx containing boundary configuration and state
 * @return PetscErrorCode 0 on success
 */
PetscErrorCode BoundarySystem_ExecuteStep(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "BoundarySystem_ExecuteStep: Starting.\n");
    
    // =========================================================================
    // PRIORITY 0: INLETS
    // =========================================================================
    {
        PetscReal local_inflow_pre = 0.0;
        PetscReal local_inflow_post = 0.0;
        PetscReal global_inflow_pre = 0.0;
        PetscReal global_inflow_post = 0.0;
        PetscInt  num_handlers = 0;
        
        LOG_ALLOW(LOCAL, LOG_TRACE, "  Priority 0 (INLETS): Begin.\n");
        
        // Phase 1: PreStep - Preparation (e.g., calculate profiles, read files)
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_INLET) continue;
            if (!handler->PreStep) continue;
            
            num_handlers++;
            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = NULL,
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    PreStep: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->PreStep(handler, &ctx, &local_inflow_pre, NULL); CHKERRQ(ierr);
        }
        
        // Optional: Global communication for PreStep (for debugging)
        if (local_inflow_pre != 0.0) {
            ierr = MPI_Allreduce(&local_inflow_pre, &global_inflow_pre, 1, MPIU_REAL,
                                MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
            LOG_ALLOW(GLOBAL, LOG_TRACE, "    PreStep predicted flux: %.6e\n", global_inflow_pre);
        }
        
        // Phase 2: Apply - Set boundary conditions
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_INLET) continue;
            
            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = NULL,
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    Apply: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->Apply(handler, &ctx); CHKERRQ(ierr);
        }
        
        // Phase 3: PostStep - Measure actual flux
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_INLET) continue;
            if (!handler->PostStep) continue;
            
            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = NULL,
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    PostStep: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->PostStep(handler, &ctx, &local_inflow_post, NULL); CHKERRQ(ierr);
        }
        
        // Phase 4: Global communication - Sum flux for other priorities to use
        ierr = MPI_Allreduce(&local_inflow_post, &global_inflow_post, 1, MPIU_REAL,
                            MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
        
        // Store for next priority levels
        user->simCtx->FluxInSum = global_inflow_post;
        
        LOG_ALLOW(GLOBAL, LOG_INFO, 
                  "  Priority 0 (INLETS): %d handler(s), FluxInSum = %.6e\n",
                  num_handlers, global_inflow_post);
    }
    
    // =========================================================================
    // PRIORITY 1: FARFIELD
    // =========================================================================
    {
        PetscReal local_farfield_in_pre = 0.0;
        PetscReal local_farfield_out_pre = 0.0;
        PetscReal local_farfield_in_post = 0.0;
        PetscReal local_farfield_out_post = 0.0;
        PetscReal global_farfield_in_pre = 0.0;
        PetscReal global_farfield_out_pre = 0.0;
        PetscReal global_farfield_in_post = 0.0;
        PetscReal global_farfield_out_post = 0.0;
        PetscInt  num_handlers = 0;
        
        LOG_ALLOW(LOCAL, LOG_TRACE, "  Priority 1 (FARFIELD): Begin.\n");
        
        // Phase 1: PreStep - Analyze flow direction, measure initial flux
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_FARFIELD) continue;
            if (!handler->PreStep) continue;
            
            num_handlers++;
            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = &user->simCtx->FluxInSum,  // Available from Priority 0
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    PreStep: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->PreStep(handler, &ctx, &local_farfield_in_pre, &local_farfield_out_pre);
            CHKERRQ(ierr);
        }
        
        // Phase 2: Global communication (optional, for debugging)
        if (local_farfield_in_pre != 0.0 || local_farfield_out_pre != 0.0) {
            ierr = MPI_Allreduce(&local_farfield_in_pre, &global_farfield_in_pre, 1, MPIU_REAL,
                                MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
            ierr = MPI_Allreduce(&local_farfield_out_pre, &global_farfield_out_pre, 1, MPIU_REAL,
                                MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
            
            LOG_ALLOW(GLOBAL, LOG_DEBUG, 
                      "    Farfield pre-analysis: In=%.6e, Out=%.6e\n",
                      global_farfield_in_pre, global_farfield_out_pre);
        }
        
        // Phase 3: Apply - Set farfield boundary conditions
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_FARFIELD) continue;
            
            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = &user->simCtx->FluxInSum,
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    Apply: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->Apply(handler, &ctx); CHKERRQ(ierr);
        }
        
        // Phase 4: PostStep - Measure actual farfield fluxes
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_FARFIELD) continue;
            if (!handler->PostStep) continue;
            
            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = &user->simCtx->FluxInSum,
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    PostStep: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->PostStep(handler, &ctx, &local_farfield_in_post, &local_farfield_out_post);
            CHKERRQ(ierr);
        }
        
        // Phase 5: Global communication - Store for outlet priority
        if (num_handlers > 0) {
            ierr = MPI_Allreduce(&local_farfield_in_post, &global_farfield_in_post, 1, MPIU_REAL,
                                MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
            ierr = MPI_Allreduce(&local_farfield_out_post, &global_farfield_out_post, 1, MPIU_REAL,
                                MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
            
            // Store for outlet handlers to use
            user->simCtx->FarFluxInSum = global_farfield_in_post;
            user->simCtx->FarFluxOutSum = global_farfield_out_post;
            
            LOG_ALLOW(GLOBAL, LOG_INFO, 
                      "  Priority 1 (FARFIELD): %d handler(s), In=%.6e, Out=%.6e\n",
                      num_handlers, global_farfield_in_post, global_farfield_out_post);
        } else {
            // No farfield handlers - zero out the fluxes
            user->simCtx->FarFluxInSum = 0.0;
            user->simCtx->FarFluxOutSum = 0.0;
        }
    }
    
    // =========================================================================
    // PRIORITY 2: WALLS
    // =========================================================================
    {
        PetscInt num_handlers = 0;
        
        LOG_ALLOW(LOCAL, LOG_TRACE, "  Priority 2 (WALLS): Begin.\n");
        
        // Phase 1: PreStep - Preparation (usually no-op for walls)
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_WALL) continue;
            if (!handler->PreStep) continue;
            
            num_handlers++;
            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = &user->simCtx->FluxInSum,
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    PreStep: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->PreStep(handler, &ctx, NULL, NULL); CHKERRQ(ierr);
        }
        
        // No global communication needed for walls
        
        // Phase 2: Apply - Set boundary conditions
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_WALL) continue;
            
            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = &user->simCtx->FluxInSum,
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    Apply: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->Apply(handler, &ctx); CHKERRQ(ierr);
        }
        
        // Phase 3: PostStep - Post-application processing (usually no-op for walls)
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_WALL) continue;
            if (!handler->PostStep) continue;
            
            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = &user->simCtx->FluxInSum,
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    PostStep: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->PostStep(handler, &ctx, NULL, NULL); CHKERRQ(ierr);
        }
        
        // No global communication needed for walls
        
        LOG_ALLOW(GLOBAL, LOG_INFO, "  Priority 2 (WALLS): %d handler(s) applied.\n",
                  num_handlers);
    }
    
    // =========================================================================
    // PRIORITY 3: OUTLETS
    // =========================================================================
    {
        PetscReal local_outflow_pre = 0.0;
        PetscReal local_outflow_post = 0.0;
        PetscReal global_outflow_pre = 0.0;
        PetscReal global_outflow_post = 0.0;
        PetscInt  num_handlers = 0;
        
        LOG_ALLOW(LOCAL, LOG_TRACE, "  Priority 3 (OUTLETS): Begin.\n");
        
        // Phase 1: PreStep - Measure uncorrected outflow (from ucat)
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_OUTLET) continue;
            if (!handler->PreStep) continue;
            
            num_handlers++;
            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = &user->simCtx->FluxInSum,      // From Priority 0
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    PreStep: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->PreStep(handler, &ctx, NULL, &local_outflow_pre); CHKERRQ(ierr);
        }
        
        // Phase 2: Global communication - Get uncorrected outflow sum
        ierr = MPI_Allreduce(&local_outflow_pre, &global_outflow_pre, 1, MPIU_REAL,
                            MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
        
        // Calculate total inflow (inlet + farfield inflow)
        PetscReal total_inflow = user->simCtx->FluxInSum + user->simCtx->FarFluxInSum;
        
        LOG_ALLOW(GLOBAL, LOG_DEBUG, 
                  "    Uncorrected outflow: %.6e, Total inflow: %.6e (Inlet: %.6e + Farfield: %.6e)\n",
                  global_outflow_pre, total_inflow, user->simCtx->FluxInSum, 
                  user->simCtx->FarFluxInSum);
        
        // Phase 3: Apply - Set corrected boundary conditions
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_OUTLET) continue;
            
            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = &user->simCtx->FluxInSum,      // From Priority 0
                .global_outflow_sum = &global_outflow_pre, // From PreStep above
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum 
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    Apply: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->Apply(handler, &ctx); CHKERRQ(ierr);
        }
        
        // Phase 4: PostStep - Measure corrected outflow (verification)
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_OUTLET) continue;
            if (!handler->PostStep) continue;
            
            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = &user->simCtx->FluxInSum,
                .global_outflow_sum = &global_outflow_pre,
                .global_farfield_inflow_sum = &user->simCtx->FarFluxInSum,
                .global_farfield_outflow_sum = &user->simCtx->FarFluxOutSum
            };
            
            LOG_ALLOW(LOCAL, LOG_TRACE, "    PostStep: Face %d (%s)\n", i, BCFaceToString((BCFace)i));
            ierr = handler->PostStep(handler, &ctx, NULL, &local_outflow_post); CHKERRQ(ierr);
        }
        
        // Phase 5: Global communication - Verify conservation
        ierr = MPI_Allreduce(&local_outflow_post, &global_outflow_post, 1, MPIU_REAL,
                            MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
        
        // Store for legacy compatibility and reporting
        user->simCtx->FluxOutSum = global_outflow_post;
        
        // Conservation check (compare total outflow vs total inflow)
        PetscReal total_outflow = global_outflow_post + user->simCtx->FarFluxOutSum;
        PetscReal flux_error = PetscAbsReal(total_outflow - total_inflow);
        PetscReal relative_error = (total_inflow > 1e-16) ? 
                                   flux_error / total_inflow : flux_error;
        
        LOG_ALLOW(GLOBAL, LOG_INFO, 
                  "  Priority 3 (OUTLETS): %d handler(s), FluxOutSum = %.6e\n",
                  num_handlers, global_outflow_post);
        LOG_ALLOW(GLOBAL, LOG_INFO, 
                  "    Conservation: Total In=%.6e, Total Out=%.6e, Error=%.3e (%.2e)%%)\n",
                  total_inflow, total_outflow, flux_error, relative_error * 100.0);
        
        if (relative_error > 1e-6) {
            LOG_ALLOW(GLOBAL, LOG_WARNING, 
                     "    WARNING: Large mass conservation error (%.2e%%)!\n",
                     relative_error * 100.0);
        }
    }
    
    LOG_ALLOW(LOCAL, LOG_VERBOSE, "Complete.\n");
    
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

//================================================================================
//
//                         PUBLIC MASTER CLEANUP FUNCTION
//
//================================================================================
#undef __FUNCT__
#define __FUNCT__ "BoundarySystem_Destroy"
/**
 * @brief Cleans up and destroys all resources allocated by the boundary system.
 *
 * This function should be called once at the end of the simulation. It iterates
 * through all created handlers and calls their respective Destroy methods to free
 * any privately allocated data (like parameter lists or handler-specific data),
 * and then frees the handler object itself. This prevents memory leaks.
 *
 * @param user The main UserCtx struct containing the boundary system to be destroyed.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode BoundarySystem_Destroy(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting destruction of all boundary handlers. \n");

    for (int i = 0; i < 6; i++) {
        BoundaryFaceConfig *face_cfg = &user->boundary_faces[i];
        const char *face_name = BCFaceToString(face_cfg->face_id);

        // --- Step 1: Free the parameter linked list associated with this face ---
        if (face_cfg->params) {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "  Freeing parameter list for Face %d (%s). \n", i, face_name);
            FreeBC_ParamList(face_cfg->params);
            face_cfg->params = NULL; // Good practice to nullify dangling pointers
        }

        // --- Step 2: Destroy the handler object itself ---
        if (face_cfg->handler) {
            const char *handler_name = BCHandlerTypeToString(face_cfg->handler->type);
            LOG_ALLOW(LOCAL, LOG_DEBUG, "  Destroying handler '%s' on Face %d (%s).\n", handler_name, i, face_name);
            
            // Call the handler's specific cleanup function first, if it exists.
            // This will free any memory stored in the handler's private `data` pointer.
            if (face_cfg->handler->Destroy) {
                ierr = face_cfg->handler->Destroy(face_cfg->handler); CHKERRQ(ierr);
            }

            // Finally, free the generic BoundaryCondition object itself.
            ierr = PetscFree(face_cfg->handler); CHKERRQ(ierr);
            face_cfg->handler = NULL;
        }
    }
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Destruction complete.\n");
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "TransferPeriodicField"
/**
 * @brief (Private) A generic routine to copy data for a single, named field across periodic boundaries.
 *
 * This function encapsulates all logic for a periodic transfer. Given a field name (e.g., "P", "Ucat"),
 * it determines the field's data type (scalar/vector), retrieves the correct DMDA and Vecs from the
 * UserCtx, and then performs the memory copy from the local ghost array to the global array.
 *
 * This must be called AFTER the corresponding local ghost vector has been updated via DMGlobalToLocal.
 *
 * @param user The main UserCtx struct, containing all grid info and field data.
 * @param field_name A string identifier for the field to transfer.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode TransferPeriodicField(UserCtx *user, const char *field_name)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info = user->info;
    PetscInt       xs = info.xs, xe = info.xs + info.xm;
    PetscInt       ys = info.ys, ye = info.ys + info.ym;
    PetscInt       zs = info.zs, ze = info.zs + info.zm;
    PetscInt       mx = info.mx, my = info.my, mz = info.mz;

    // --- Local variables to hold the specific details of the chosen field ---
    DM        dm;
    Vec       global_vec;
    Vec       local_vec;
    PetscInt  dof;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // --- STEP 1: Dispatcher - Set the specific DM, Vecs, and dof based on field_name ---
    if (strcmp(field_name, "Ucat") == 0) {
        dm         = user->fda;
        global_vec = user->Ucat;
        local_vec  = user->lUcat;
        dof        = 3;
    } else if (strcmp(field_name, "P") == 0) {
        dm         = user->da;
        global_vec = user->P;
        local_vec  = user->lP;
        dof        = 1;
    } else if (strcmp(field_name, "Nvert") == 0) {
        dm         = user->da;
        global_vec = user->Nvert;
        local_vec  = user->lNvert;
        dof        = 1;
    } 
    /*
    // Example for future extension:
    else if (strcmp(field_name, "Temperature") == 0) {
        dm         = user->da; // Assuming Temperature is scalar
        global_vec = user->T;
        local_vec  = user->lT;
        dof        = 1;
    }
    */
    else {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_UNKNOWN_TYPE, "Unknown field name '%s' in TransferPeriodicFieldByName.", field_name);
    }

    // --- STEP 2: Execute the copy logic using the dispatched variables ---
    if (dof == 1) { // --- Handle SCALAR fields (PetscReal) ---
        PetscReal ***g_array, ***l_array;
        ierr = DMDAVecGetArray(dm, global_vec, &g_array); CHKERRQ(ierr);
        ierr = DMDAVecGetArray(dm, local_vec, &l_array); CHKERRQ(ierr);

        if (user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC && xs == 0) for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) g_array[k][j][xs] = l_array[k][j][xs-2];
        if (user->boundary_faces[BC_FACE_POS_X].mathematical_type == PERIODIC && xe == mx) for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) g_array[k][j][xe-1] = l_array[k][j][xe+1];
        if (user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC && ys == 0) for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) g_array[k][ys][i] = l_array[k][ys-2][i];
        if (user->boundary_faces[BC_FACE_POS_Y].mathematical_type == PERIODIC && ye == my) for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) g_array[k][ye-1][i] = l_array[k][ye+1][i];
        if (user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC && zs == 0) for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) g_array[zs][j][i] = l_array[zs-2][j][i];
        if (user->boundary_faces[BC_FACE_POS_Z].mathematical_type == PERIODIC && ze == mz) for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) g_array[ze-1][j][i] = l_array[ze+1][j][i];
        
        ierr = DMDAVecRestoreArray(dm, global_vec, &g_array); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(dm, local_vec, &l_array); CHKERRQ(ierr);

    } else { // --- Handle VECTOR fields (Cmpnts) ---
        Cmpnts ***g_array, ***l_array;
        ierr = DMDAVecGetArray(dm, global_vec, &g_array); CHKERRQ(ierr);
        ierr = DMDAVecGetArray(dm, local_vec, &l_array); CHKERRQ(ierr);

        if (user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC && xs == 0) for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) g_array[k][j][xs] = l_array[k][j][xs-2];
        if (user->boundary_faces[BC_FACE_POS_X].mathematical_type == PERIODIC && xe == mx) for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) g_array[k][j][xe-1] = l_array[k][j][xe+1];
        if (user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC && ys == 0) for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) g_array[k][ys][i] = l_array[k][ys-2][i];
        if (user->boundary_faces[BC_FACE_POS_Y].mathematical_type == PERIODIC && ye == my) for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) g_array[k][ye-1][i] = l_array[k][ye+1][i];
        if (user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC && zs == 0) for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) g_array[zs][j][i] = l_array[zs-2][j][i];
        if (user->boundary_faces[BC_FACE_POS_Z].mathematical_type == PERIODIC && ze == mz) for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) g_array[ze-1][j][i] = l_array[ze+1][j][i];
        
        ierr = DMDAVecRestoreArray(dm, global_vec, &g_array); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(dm, local_vec, &l_array); CHKERRQ(ierr);
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ApplyPeriodicBCs"
/**
 * @brief Applies periodic boundary conditions by copying data across domain boundaries for all relevant fields.
 *
 * This function orchestrates the periodic update. It first performs a single, collective
 * ghost-cell exchange for all fields. Then, it calls a generic helper routine to perform
 * the memory copy for each individual field by name.
 *
 * @param user The main UserCtx struct.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ApplyPeriodicBCs(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscBool is_any_periodic = PETSC_FALSE;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    for (int i = 0; i < 6; i++) {
        if (user->boundary_faces[i].mathematical_type == PERIODIC) {
            is_any_periodic = PETSC_TRUE;
            break;
        }
    }

    if (!is_any_periodic) {
        LOG_ALLOW(GLOBAL,LOG_TRACE, "No periodic boundaries defined; skipping ApplyPeriodicBCs.\n");
        PROFILE_FUNCTION_END;
        PetscFunctionReturn(0);
    }

    LOG_ALLOW(GLOBAL, LOG_TRACE, "Applying periodic boundary conditions for all fields.\s");

    // STEP 1: Perform the collective communication for all fields at once.
    ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "P"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "Nvert"); CHKERRQ(ierr);
    /* if (user->solve_temperature) { ierr = UpdateLocalGhosts(user, "Temperature"); CHKERRQ(ierr); } */

    // STEP 2: Call the generic copy routine for each field by name.
    ierr = TransferPeriodicField(user, "Ucat"); CHKERRQ(ierr);
    ierr = TransferPeriodicField(user, "P"); CHKERRQ(ierr);
    ierr = TransferPeriodicField(user, "Nvert"); CHKERRQ(ierr);

    // FUTURE EXTENSION: Adding a new scalar field like Temperature is now trivial.
    /*
    if (user->solve_temperature) {
        ierr = TransferPeriodicField(user, "Temperature"); CHKERRQ(ierr);
    }
    */

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdateDummyCells"
/**
 * @brief Updates the dummy cells (ghost nodes) on the faces of the local domain for NON-PERIODIC boundaries.
 *
 * This function's role is to apply a second-order extrapolation to set the ghost
 * cell values based on the boundary condition value (stored in `ubcs`) and the
 * first interior cell.
 *
 * NOTE: This function deliberately IGNORES periodic boundaries. It is part of a
 * larger workflow where `ApplyPeriodicBCs` handles periodic faces first.
 *
 * CRITICAL DETAIL: This function uses shrunken loop ranges (lxs, lxe, etc.) to
 * intentionally update only the flat part of the faces, avoiding the edges and

 * corners. The edges and corners are then handled separately by `UpdateCornerNodes`.
 * This precisely replicates the logic of the original function.
 *
 * @param user The main UserCtx struct containing all necessary data.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode UpdateDummyCells(UserCtx *user)
{
    PetscErrorCode ierr;
    DM            fda = user->fda;
    DMDALocalInfo info = user->info;
    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    // --- Calculate shrunken loop ranges to avoid edges and corners ---
    PetscInt lxs = xs, lxe = xe;
    PetscInt lys = ys, lye = ye;
    PetscInt lzs = zs, lze = ze;

    if (xs == 0) lxs = xs + 1;
    if (ys == 0) lys = ys + 1;
    if (zs == 0) lzs = zs + 1;

    if (xe == mx) lxe = xe - 1;
    if (ye == my) lye = ye - 1;
    if (ze == mz) lze = ze - 1;

    Cmpnts        ***ucat, ***ubcs;
    PetscFunctionBeginUser;

    ierr = DMDAVecGetArray(fda, user->Bcs.Ubcs, &ubcs); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(fda, user->Ucat, &ucat); CHKERRQ(ierr);

    // -X Face
    if (user->boundary_faces[BC_FACE_NEG_X].mathematical_type != PERIODIC && xs == 0) {
        for (PetscInt k = lzs; k < lze; k++) for (PetscInt j = lys; j < lye; j++) {
            ucat[k][j][xs].x = 2.0 * ubcs[k][j][xs].x - ucat[k][j][xs + 1].x;
            ucat[k][j][xs].y = 2.0 * ubcs[k][j][xs].y - ucat[k][j][xs + 1].y;
            ucat[k][j][xs].z = 2.0 * ubcs[k][j][xs].z - ucat[k][j][xs + 1].z;
        }
    }
    // +X Face
    if (user->boundary_faces[BC_FACE_POS_X].mathematical_type != PERIODIC && xe == mx) {
        for (PetscInt k = lzs; k < lze; k++) for (PetscInt j = lys; j < lye; j++) {
            ucat[k][j][xe-1].x = 2.0 * ubcs[k][j][xe-1].x - ucat[k][j][xe - 2].x;
            ucat[k][j][xe-1].y = 2.0 * ubcs[k][j][xe-1].y - ucat[k][j][xe - 2].y;
            ucat[k][j][xe-1].z = 2.0 * ubcs[k][j][xe-1].z - ucat[k][j][xe - 2].z;
        }
    }

    // -Y Face
    if (user->boundary_faces[BC_FACE_NEG_Y].mathematical_type != PERIODIC && ys == 0) {
        for (PetscInt k = lzs; k < lze; k++) for (PetscInt i = lxs; i < lxe; i++) {
            ucat[k][ys][i].x = 2.0 * ubcs[k][ys][i].x - ucat[k][ys + 1][i].x;
            ucat[k][ys][i].y = 2.0 * ubcs[k][ys][i].y - ucat[k][ys + 1][i].y;
            ucat[k][ys][i].z = 2.0 * ubcs[k][ys][i].z - ucat[k][ys + 1][i].z;
        }
    }
    // +Y Face
    if (user->boundary_faces[BC_FACE_POS_Y].mathematical_type != PERIODIC && ye == my) {
        for (PetscInt k = lzs; k < lze; k++) for (PetscInt i = lxs; i < lxe; i++) {
            ucat[k][ye-1][i].x = 2.0 * ubcs[k][ye-1][i].x - ucat[k][ye-2][i].x;
            ucat[k][ye-1][i].y = 2.0 * ubcs[k][ye-1][i].y - ucat[k][ye-2][i].y;
            ucat[k][ye-1][i].z = 2.0 * ubcs[k][ye-1][i].z - ucat[k][ye-2][i].z;
        }
    }

    // -Z Face
    if (user->boundary_faces[BC_FACE_NEG_Z].mathematical_type != PERIODIC && zs == 0) {
        for (PetscInt j = lys; j < lye; j++) for (PetscInt i = lxs; i < lxe; i++) {
            ucat[zs][j][i].x = 2.0 * ubcs[zs][j][i].x - ucat[zs + 1][j][i].x;
            ucat[zs][j][i].y = 2.0 * ubcs[zs][j][i].y - ucat[zs + 1][j][i].y;
            ucat[zs][j][i].z = 2.0 * ubcs[zs][j][i].z - ucat[zs + 1][j][i].z;
        }
    }
    // +Z Face
    if (user->boundary_faces[BC_FACE_POS_Z].mathematical_type != PERIODIC && ze == mz) {
        for (PetscInt j = lys; j < lye; j++) for (PetscInt i = lxs; i < lxe; i++) {
            ucat[ze-1][j][i].x = 2.0 * ubcs[ze-1][j][i].x - ucat[ze-2][j][i].x;
            ucat[ze-1][j][i].y = 2.0 * ubcs[ze-1][j][i].y - ucat[ze-2][j][i].y;
            ucat[ze-1][j][i].z = 2.0 * ubcs[ze-1][j][i].z - ucat[ze-2][j][i].z;
        }
    }

    ierr = DMDAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(fda, user->Ucat, &ucat); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdateCornerNodes"
/**
 * @brief Updates the corner and edge ghost nodes of the local domain by averaging.
 *
 * This function should be called AFTER the face ghost nodes are finalized by both
 * `ApplyPeriodicBCs` and `UpdateDummyCells`. It resolves the values at shared
 * edges and corners by averaging the values of adjacent, previously-computed
 * ghost nodes.
 *
 * The logic is generic and works correctly regardless of the boundary types on
 * the adjacent faces (e.g., it will correctly average a periodic face neighbor
 * with a wall face neighbor).
 *
 * @param user The main UserCtx struct containing all necessary data.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode UpdateCornerNodes(UserCtx *user)
{
    PetscErrorCode ierr;
    DM            da = user->da, fda = user->fda;
    DMDALocalInfo info = user->info;
    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***ucat;
    PetscReal     ***p;

    PetscFunctionBeginUser;

    ierr = DMDAVecGetArray(fda, user->Ucat, &ucat); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(da, user->P, &p); CHKERRQ(ierr);

    // --- Update Edges and Corners by Averaging ---
    // The order of these blocks ensures that corners (where 3 faces meet) are
    // computed using data from edges (where 2 faces meet), which are computed first.
// Edges connected to the -Z face (k=zs)
  if (zs == 0) {
      if (xs == 0) {
          for (PetscInt j = ys; j < ye; j++) {
              p[zs][j][xs] = 0.5 * (p[zs+1][j][xs] + p[zs][j][xs+1]);
              ucat[zs][j][xs].x = 0.5 * (ucat[zs+1][j][xs].x + ucat[zs][j][xs+1].x);
              ucat[zs][j][xs].y = 0.5 * (ucat[zs+1][j][xs].y + ucat[zs][j][xs+1].y);
              ucat[zs][j][xs].z = 0.5 * (ucat[zs+1][j][xs].z + ucat[zs][j][xs+1].z);
          }
      }
      if (xe == mx) {
          for (PetscInt j = ys; j < ye; j++) {
              p[zs][j][mx-1] = 0.5 * (p[zs+1][j][mx-1] + p[zs][j][mx-2]);
              ucat[zs][j][mx-1].x = 0.5 * (ucat[zs+1][j][mx-1].x + ucat[zs][j][mx-2].x);
              ucat[zs][j][mx-1].y = 0.5 * (ucat[zs+1][j][mx-1].y + ucat[zs][j][mx-2].y);
              ucat[zs][j][mx-1].z = 0.5 * (ucat[zs+1][j][mx-1].z + ucat[zs][j][mx-2].z);
          }
      }
      if (ys == 0) {
          for (PetscInt i = xs; i < xe; i++) {
              p[zs][ys][i] = 0.5 * (p[zs+1][ys][i] + p[zs][ys+1][i]);
              ucat[zs][ys][i].x = 0.5 * (ucat[zs+1][ys][i].x + ucat[zs][ys+1][i].x);
              ucat[zs][ys][i].y = 0.5 * (ucat[zs+1][ys][i].y + ucat[zs][ys+1][i].y);
              ucat[zs][ys][i].z = 0.5 * (ucat[zs+1][ys][i].z + ucat[zs][ys+1][i].z);
          }
      }
      if (ye == my) {
          for (PetscInt i = xs; i < xe; i++) {
              p[zs][my-1][i] = 0.5 * (p[zs+1][my-1][i] + p[zs][my-2][i]);
              ucat[zs][my-1][i].x = 0.5 * (ucat[zs+1][my-1][i].x + ucat[zs][my-2][i].x);
              ucat[zs][my-1][i].y = 0.5 * (ucat[zs+1][my-1][i].y + ucat[zs][my-2][i].y);
              ucat[zs][my-1][i].z = 0.5 * (ucat[zs+1][my-1][i].z + ucat[zs][my-2][i].z);
          }
      }
  }

  // Edges connected to the +Z face (k=ze-1)
  if (ze == mz) {
      if (xs == 0) {
          for (PetscInt j = ys; j < ye; j++) {
              p[mz-1][j][xs] = 0.5 * (p[mz-2][j][xs] + p[mz-1][j][xs+1]);
              ucat[mz-1][j][xs].x = 0.5 * (ucat[mz-2][j][xs].x + ucat[mz-1][j][xs+1].x);
              ucat[mz-1][j][xs].y = 0.5 * (ucat[mz-2][j][xs].y + ucat[mz-1][j][xs+1].y);
              ucat[mz-1][j][xs].z = 0.5 * (ucat[mz-2][j][xs].z + ucat[mz-1][j][xs+1].z);
          }
      }
      if (xe == mx) {
          for (PetscInt j = ys; j < ye; j++) {
              p[mz-1][j][mx-1] = 0.5 * (p[mz-2][j][mx-1] + p[mz-1][j][mx-2]);
              ucat[mz-1][j][mx-1].x = 0.5 * (ucat[mz-2][j][mx-1].x + ucat[mz-1][j][mx-2].x);
              ucat[mz-1][j][mx-1].y = 0.5 * (ucat[mz-2][j][mx-1].y + ucat[mz-1][j][mx-2].y);
              ucat[mz-1][j][mx-1].z = 0.5 * (ucat[mz-2][j][mx-1].z + ucat[mz-1][j][mx-2].z);
          }
      }
      if (ys == 0) {
          for (PetscInt i = xs; i < xe; i++) {
              p[mz-1][ys][i] = 0.5 * (p[mz-2][ys][i] + p[mz-1][ys+1][i]);
              ucat[mz-1][ys][i].x = 0.5 * (ucat[mz-2][ys][i].x + ucat[mz-1][ys+1][i].x);
              ucat[mz-1][ys][i].y = 0.5 * (ucat[mz-2][ys][i].y + ucat[mz-1][ys+1][i].y);
              ucat[mz-1][ys][i].z = 0.5 * (ucat[mz-2][ys][i].z + ucat[mz-1][ys+1][i].z);
          }
      }
      if (ye == my) {
          for (PetscInt i = xs; i < xe; i++) {
              p[mz-1][my-1][i] = 0.5 * (p[mz-2][my-1][i] + p[mz-1][my-2][i]);
              ucat[mz-1][my-1][i].x = 0.5 * (ucat[mz-2][my-1][i].x + ucat[mz-1][my-2][i].x);
              ucat[mz-1][my-1][i].y = 0.5 * (ucat[mz-2][my-1][i].y + ucat[mz-1][my-2][i].y);
              ucat[mz-1][my-1][i].z = 0.5 * (ucat[mz-2][my-1][i].z + ucat[mz-1][my-2][i].z);
          }
      }
  }

  // Remaining edges on the XY plane (that are not on Z faces)
  if (ys == 0) {
      if (xs == 0) {
          for (PetscInt k = zs; k < ze; k++) {
              p[k][ys][xs] = 0.5 * (p[k][ys+1][xs] + p[k][ys][xs+1]);
              ucat[k][ys][xs].x = 0.5 * (ucat[k][ys+1][xs].x + ucat[k][ys][xs+1].x);
              ucat[k][ys][xs].y = 0.5 * (ucat[k][ys+1][xs].y + ucat[k][ys][xs+1].y);
              ucat[k][ys][xs].z = 0.5 * (ucat[k][ys+1][xs].z + ucat[k][ys][xs+1].z);
          }
      }
      if (xe == mx) {
          for (PetscInt k = zs; k < ze; k++) {
              p[k][ys][mx-1] = 0.5 * (p[k][ys+1][mx-1] + p[k][ys][mx-2]);
              ucat[k][ys][mx-1].x = 0.5 * (ucat[k][ys+1][mx-1].x + ucat[k][ys][mx-2].x);
              ucat[k][ys][mx-1].y = 0.5 * (ucat[k][ys+1][mx-1].y + ucat[k][ys][mx-2].y);
              ucat[k][ys][mx-1].z = 0.5 * (ucat[k][ys+1][mx-1].z + ucat[k][ys][mx-2].z);
          }
      }
  }

  if (ye == my) {
      if (xs == 0) {
          for (PetscInt k = zs; k < ze; k++) {
              p[k][my-1][xs] = 0.5 * (p[k][my-2][xs] + p[k][my-1][xs+1]);
              ucat[k][my-1][xs].x = 0.5 * (ucat[k][my-2][xs].x + ucat[k][my-1][xs+1].x);
              ucat[k][my-1][xs].y = 0.5 * (ucat[k][my-2][xs].y + ucat[k][my-1][xs+1].y);
              ucat[k][my-1][xs].z = 0.5 * (ucat[k][my-2][xs].z + ucat[k][my-1][xs+1].z);
          }
      }
      if (xe == mx) {
          for (PetscInt k = zs; k < ze; k++) {
              p[k][my-1][mx-1] = 0.5 * (p[k][my-2][mx-1] + p[k][my-1][mx-2]);
              ucat[k][my-1][mx-1].x = 0.5 * (ucat[k][my-2][mx-1].x + ucat[k][my-1][mx-2].x);
              ucat[k][my-1][mx-1].y = 0.5 * (ucat[k][my-2][mx-1].y + ucat[k][my-1][mx-2].y);
              ucat[k][my-1][mx-1].z = 0.5 * (ucat[k][my-2][mx-1].z + ucat[k][my-1][mx-2].z);
          }
      }
  }

    ierr = DMDAVecRestoreArray(fda, user->Ucat, &ucat); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(da, user->P, &p); CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ApplyWallFunction"
/**
 * @brief Applies wall function modeling to near-wall velocities for all wall-type boundaries.
 *
 * This function implements log-law wall functions to model the near-wall velocity profile
 * without fully resolving the viscous sublayer. It is applicable to ALL wall-type boundaries
 * regardless of their specific boundary condition (no-slip, moving wall, slip, etc.), as
 * determined by the mathematical_type being WALL.
 *
 * MATHEMATICAL BACKGROUND:
 * Wall functions bridge the gap between the wall (y=0) and the first computational cell
 * center by using empirical log-law relationships:
 *   - Viscous sublayer (y+ < 11.81): u+ = y+
 *   - Log-law region (y+ > 11.81): u+ = (1/κ) * ln(E * y+)
 * where u+ = u/u_τ, y+ = y*u_τ/ν, κ = 0.41 (von Karman constant), E = exp(κB)
 *
 * IMPLEMENTATION DETAILS:
 * Unlike standard boundary conditions that set ghost cell values, wall functions:
 *   1. Read velocity from the SECOND interior cell (i±2, j±2, k±2)
 *   2. Compute wall shear stress using log-law
 *   3. Modify velocity at the FIRST interior cell (i±1, j±1, k±1) 
 *   4. Keep ghost cell boundary values (ubcs, ucont) at zero
 *
 * WORKFLOW:
 *   - Called from ApplyBoundaryConditions after standard BC application
 *   - Operates on ucat (Cartesian velocity)
 *   - Updates ustar (friction velocity field) for diagnostics/turbulence models
 *   - Ghost cells remain zero; UpdateDummyCells handles extrapolation afterward
 *
 * GEOMETRIC QUANTITIES:
 *   sb = wall-normal distance from wall to first interior cell center
 *   sc = wall-normal distance from wall to second interior cell center  
 *   These are computed from cell Jacobians (aj) and face area vectors
 *
 * APPLICABILITY:
 *   - Requires simCtx->wallfunction = true
 *   - Only processes faces where mathematical_type == WALL
 *   - Skips solid-embedded cells (nvert >= 0.1)
 *
 * @param user The UserCtx containing all simulation state and geometry
 * @return PetscErrorCode 0 on success
 *
 * @note This function modifies interior cell velocities, NOT ghost cells
 * @note Wall roughness (ks) is currently set to 1e-16 (smooth wall)
 * @see wall_function_loglaw() in wallfunction.c for the actual log-law implementation
 * @see noslip() in wallfunction.c for the initial linear interpolation
 */
PetscErrorCode ApplyWallFunction(UserCtx *user)
{
    PetscErrorCode ierr;
    SimCtx *simCtx = user->simCtx;
    DMDALocalInfo *info = &user->info;
    
    PetscFunctionBeginUser;
    
    // =========================================================================
    // STEP 0: Early exit if wall functions are disabled
    // =========================================================================
    if (!simCtx->wallfunction) {
        PetscFunctionReturn(0);
    }
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Processing wall function boundaries.\n");
    
    // =========================================================================
    // STEP 1: Get read/write access to all necessary field arrays
    // =========================================================================
    Cmpnts ***velocity_cartesian;           // Cartesian velocity (modified)
    Cmpnts ***velocity_contravariant;       // Contravariant velocity (set to zero at walls)
    Cmpnts ***velocity_boundary;            // Boundary condition velocity (kept at zero)
    Cmpnts ***csi, ***eta, ***zet;         // Metric tensor components (face normals)
    PetscReal ***node_vertex_flag;          // Fluid/solid indicator (0=fluid, 1=solid)
    PetscReal ***cell_jacobian;             // Grid Jacobian (1/volume)
    PetscReal ***friction_velocity;         // u_tau (friction velocity field)
    
    ierr = DMDAVecGetArray(user->fda, user->Ucat, &velocity_cartesian); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Ucont, &velocity_contravariant); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Bcs.Ubcs, &velocity_boundary); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi, (const Cmpnts***)&csi); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta, (const Cmpnts***)&eta); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet, (const Cmpnts***)&zet); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da, user->lNvert, (const PetscReal***)&node_vertex_flag); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da, user->lAj, (const PetscReal***)&cell_jacobian); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->lFriction_Velocity, &friction_velocity); CHKERRQ(ierr);
    
    // =========================================================================
    // STEP 2: Define loop bounds (owned portion of the grid for this MPI rank)
    // =========================================================================
    PetscInt grid_start_i = info->xs, grid_end_i = info->xs + info->xm;
    PetscInt grid_start_j = info->ys, grid_end_j = info->ys + info->ym;
    PetscInt grid_start_k = info->zs, grid_end_k = info->zs + info->zm;
    PetscInt grid_size_i = info->mx, grid_size_j = info->my, grid_size_k = info->mz;
    
    // Shrunken loop bounds: exclude domain edges and corners to avoid double-counting
    PetscInt loop_start_i = grid_start_i, loop_end_i = grid_end_i;
    PetscInt loop_start_j = grid_start_j, loop_end_j = grid_end_j;
    PetscInt loop_start_k = grid_start_k, loop_end_k = grid_end_k;
    
    if (grid_start_i == 0) loop_start_i = grid_start_i + 1;
    if (grid_end_i == grid_size_i) loop_end_i = grid_end_i - 1;
    if (grid_start_j == 0) loop_start_j = grid_start_j + 1;
    if (grid_end_j == grid_size_j) loop_end_j = grid_end_j - 1;
    if (grid_start_k == 0) loop_start_k = grid_start_k + 1;
    if (grid_end_k == grid_size_k) loop_end_k = grid_end_k - 1;
    
    // Wall roughness parameter (currently smooth wall)
    const PetscReal wall_roughness_height = 1.e-16;
    
    // =========================================================================
    // STEP 3: Process each of the 6 domain faces
    // =========================================================================
    for (int face_index = 0; face_index < 6; face_index++) {
        BCFace current_face_id = (BCFace)face_index;
        BoundaryFaceConfig *face_config = &user->boundary_faces[current_face_id];
        
        // Only process faces that are mathematical walls (applies to no-slip, moving, slip, etc.)
        if (face_config->mathematical_type != WALL) {
            continue;
        }
        
        // Check if this MPI rank owns part of this face
        PetscBool rank_owns_this_face;
        ierr = CanRankServiceFace(info, user->IM, user->JM, user->KM, 
                                  current_face_id, &rank_owns_this_face); CHKERRQ(ierr);
        
        if (!rank_owns_this_face) {
            continue;
        }
        
        LOG_ALLOW(LOCAL, LOG_TRACE, "Processing Face %d (%s)\n",
                  current_face_id, BCFaceToString(current_face_id));
        
        // =====================================================================
        // Process each face with appropriate indexing
        // =====================================================================
        switch(current_face_id) {
            
            // =================================================================
            // NEGATIVE X FACE (i = 0, normal points in +X direction)
            // =================================================================
            case BC_FACE_NEG_X: {
                if (grid_start_i == 0) {
                    const PetscInt ghost_cell_index = grid_start_i;
                    const PetscInt first_interior_cell = grid_start_i + 1;
                    const PetscInt second_interior_cell = grid_start_i + 2;
                    
                    for (PetscInt k = loop_start_k; k < loop_end_k; k++) {
                        for (PetscInt j = loop_start_j; j < loop_end_j; j++) {
                            
                            // Skip if this is a solid cell (embedded boundary)
                            if (node_vertex_flag[k][j][first_interior_cell] < 0.1) {
                                
                                // Calculate face area from contravariant metric tensor
                                PetscReal face_area = sqrt(
                                    csi[k][j][ghost_cell_index].x * csi[k][j][ghost_cell_index].x + 
                                    csi[k][j][ghost_cell_index].y * csi[k][j][ghost_cell_index].y + 
                                    csi[k][j][ghost_cell_index].z * csi[k][j][ghost_cell_index].z
                                );
                                
                                // Compute wall-normal distances using cell Jacobians
                                // sb = distance from wall to first interior cell center
                                // sc = distance from wall to second interior cell center
                                PetscReal distance_to_first_cell = 0.5 / cell_jacobian[k][j][first_interior_cell] / face_area;
                                PetscReal distance_to_second_cell = 2.0 * distance_to_first_cell + 
                                                                     0.5 / cell_jacobian[k][j][second_interior_cell] / face_area;
                                
                                // Compute unit normal vector pointing INTO the domain
                                PetscReal wall_normal[3];
                                wall_normal[0] = csi[k][j][ghost_cell_index].x / face_area;
                                wall_normal[1] = csi[k][j][ghost_cell_index].y / face_area;
                                wall_normal[2] = csi[k][j][ghost_cell_index].z / face_area;
                                
                                // Define velocities for wall function calculation
                                Cmpnts wall_velocity;      // Ua = velocity at wall (zero for stationary wall)
                                Cmpnts reference_velocity; // Uc = velocity at second interior cell
                                
                                wall_velocity.x = wall_velocity.y = wall_velocity.z = 0.0;
                                reference_velocity = velocity_cartesian[k][j][second_interior_cell];
                                
                                // Step 1: Linear interpolation (provides initial guess)
                                noslip(user, distance_to_second_cell, distance_to_first_cell,
                                      wall_velocity, reference_velocity,
                                      &velocity_cartesian[k][j][first_interior_cell],
                                      wall_normal[0], wall_normal[1], wall_normal[2]);
                                
                                // Step 2: Apply log-law correction (improves near-wall velocity)
                                wall_function_loglaw(user, wall_roughness_height,
                                                    distance_to_second_cell, distance_to_first_cell,
                                                    wall_velocity, reference_velocity,
                                                    &velocity_cartesian[k][j][first_interior_cell],
                                                    &friction_velocity[k][j][first_interior_cell],
                                                    wall_normal[0], wall_normal[1], wall_normal[2]);
                                
                                // Ensure ghost cell BC remains zero (required for proper extrapolation)
                                velocity_boundary[k][j][ghost_cell_index].x = 0.0;
                                velocity_boundary[k][j][ghost_cell_index].y = 0.0;
                                velocity_boundary[k][j][ghost_cell_index].z = 0.0;
                                velocity_contravariant[k][j][ghost_cell_index].x = 0.0;
                            }
                        }
                    }
                }
            } break;
            
            // =================================================================
            // POSITIVE X FACE (i = mx-1, normal points in -X direction)
            // =================================================================
            case BC_FACE_POS_X: {
                if (grid_end_i == grid_size_i) {
                    const PetscInt ghost_cell_index = grid_end_i - 1;
                    const PetscInt first_interior_cell = grid_end_i - 2;
                    const PetscInt second_interior_cell = grid_end_i - 3;
                    
                    for (PetscInt k = loop_start_k; k < loop_end_k; k++) {
                        for (PetscInt j = loop_start_j; j < loop_end_j; j++) {
                            
                            if (node_vertex_flag[k][j][first_interior_cell] < 0.1) {
                                
                                PetscReal face_area = sqrt(
                                    csi[k][j][first_interior_cell].x * csi[k][j][first_interior_cell].x + 
                                    csi[k][j][first_interior_cell].y * csi[k][j][first_interior_cell].y + 
                                    csi[k][j][first_interior_cell].z * csi[k][j][first_interior_cell].z
                                );
                                
                                PetscReal distance_to_first_cell = 0.5 / cell_jacobian[k][j][first_interior_cell] / face_area;
                                PetscReal distance_to_second_cell = 2.0 * distance_to_first_cell + 
                                                                     0.5 / cell_jacobian[k][j][second_interior_cell] / face_area;
                                
                                // Note: Normal flipped for +X face to point INTO domain
                                PetscReal wall_normal[3];
                                wall_normal[0] = -csi[k][j][first_interior_cell].x / face_area;
                                wall_normal[1] = -csi[k][j][first_interior_cell].y / face_area;
                                wall_normal[2] = -csi[k][j][first_interior_cell].z / face_area;
                                
                                Cmpnts wall_velocity, reference_velocity;
                                wall_velocity.x = wall_velocity.y = wall_velocity.z = 0.0;
                                reference_velocity = velocity_cartesian[k][j][second_interior_cell];
                                
                                noslip(user, distance_to_second_cell, distance_to_first_cell,
                                      wall_velocity, reference_velocity,
                                      &velocity_cartesian[k][j][first_interior_cell],
                                      wall_normal[0], wall_normal[1], wall_normal[2]);
                                
                                wall_function_loglaw(user, wall_roughness_height,
                                                    distance_to_second_cell, distance_to_first_cell,
                                                    wall_velocity, reference_velocity,
                                                    &velocity_cartesian[k][j][first_interior_cell],
                                                    &friction_velocity[k][j][first_interior_cell],
                                                    wall_normal[0], wall_normal[1], wall_normal[2]);
                                
                                velocity_boundary[k][j][ghost_cell_index].x = 0.0;
                                velocity_boundary[k][j][ghost_cell_index].y = 0.0;
                                velocity_boundary[k][j][ghost_cell_index].z = 0.0;
                                velocity_contravariant[k][j][first_interior_cell].x = 0.0;
                            }
                        }
                    }
                }
            } break;
            
            // =================================================================
            // NEGATIVE Y FACE (j = 0, normal points in +Y direction)
            // =================================================================
            case BC_FACE_NEG_Y: {
                if (grid_start_j == 0) {
                    const PetscInt ghost_cell_index = grid_start_j;
                    const PetscInt first_interior_cell = grid_start_j + 1;
                    const PetscInt second_interior_cell = grid_start_j + 2;
                    
                    for (PetscInt k = loop_start_k; k < loop_end_k; k++) {
                        for (PetscInt i = loop_start_i; i < loop_end_i; i++) {
                            
                            if (node_vertex_flag[k][first_interior_cell][i] < 0.1) {
                                
                                PetscReal face_area = sqrt(
                                    eta[k][ghost_cell_index][i].x * eta[k][ghost_cell_index][i].x + 
                                    eta[k][ghost_cell_index][i].y * eta[k][ghost_cell_index][i].y + 
                                    eta[k][ghost_cell_index][i].z * eta[k][ghost_cell_index][i].z
                                );
                                
                                PetscReal distance_to_first_cell = 0.5 / cell_jacobian[k][first_interior_cell][i] / face_area;
                                PetscReal distance_to_second_cell = 2.0 * distance_to_first_cell + 
                                                                     0.5 / cell_jacobian[k][second_interior_cell][i] / face_area;
                                
                                PetscReal wall_normal[3];
                                wall_normal[0] = eta[k][ghost_cell_index][i].x / face_area;
                                wall_normal[1] = eta[k][ghost_cell_index][i].y / face_area;
                                wall_normal[2] = eta[k][ghost_cell_index][i].z / face_area;
                                
                                Cmpnts wall_velocity, reference_velocity;
                                wall_velocity.x = wall_velocity.y = wall_velocity.z = 0.0;
                                reference_velocity = velocity_cartesian[k][second_interior_cell][i];
                                
                                noslip(user, distance_to_second_cell, distance_to_first_cell,
                                      wall_velocity, reference_velocity,
                                      &velocity_cartesian[k][first_interior_cell][i],
                                      wall_normal[0], wall_normal[1], wall_normal[2]);
                                
                                wall_function_loglaw(user, wall_roughness_height,
                                                    distance_to_second_cell, distance_to_first_cell,
                                                    wall_velocity, reference_velocity,
                                                    &velocity_cartesian[k][first_interior_cell][i],
                                                    &friction_velocity[k][first_interior_cell][i],
                                                    wall_normal[0], wall_normal[1], wall_normal[2]);
                                
                                velocity_boundary[k][ghost_cell_index][i].x = 0.0;
                                velocity_boundary[k][ghost_cell_index][i].y = 0.0;
                                velocity_boundary[k][ghost_cell_index][i].z = 0.0;
                                velocity_contravariant[k][ghost_cell_index][i].y = 0.0;
                            }
                        }
                    }
                }
            } break;
            
            // =================================================================
            // POSITIVE Y FACE (j = my-1, normal points in -Y direction)
            // =================================================================
            case BC_FACE_POS_Y: {
                if (grid_end_j == grid_size_j) {
                    const PetscInt ghost_cell_index = grid_end_j - 1;
                    const PetscInt first_interior_cell = grid_end_j - 2;
                    const PetscInt second_interior_cell = grid_end_j - 3;
                    
                    for (PetscInt k = loop_start_k; k < loop_end_k; k++) {
                        for (PetscInt i = loop_start_i; i < loop_end_i; i++) {
                            
                            if (node_vertex_flag[k][first_interior_cell][i] < 0.1) {
                                
                                PetscReal face_area = sqrt(
                                    eta[k][first_interior_cell][i].x * eta[k][first_interior_cell][i].x + 
                                    eta[k][first_interior_cell][i].y * eta[k][first_interior_cell][i].y + 
                                    eta[k][first_interior_cell][i].z * eta[k][first_interior_cell][i].z
                                );
                                
                                PetscReal distance_to_first_cell = 0.5 / cell_jacobian[k][first_interior_cell][i] / face_area;
                                PetscReal distance_to_second_cell = 2.0 * distance_to_first_cell + 
                                                                     0.5 / cell_jacobian[k][second_interior_cell][i] / face_area;
                                
                                PetscReal wall_normal[3];
                                wall_normal[0] = -eta[k][first_interior_cell][i].x / face_area;
                                wall_normal[1] = -eta[k][first_interior_cell][i].y / face_area;
                                wall_normal[2] = -eta[k][first_interior_cell][i].z / face_area;
                                
                                Cmpnts wall_velocity, reference_velocity;
                                wall_velocity.x = wall_velocity.y = wall_velocity.z = 0.0;
                                reference_velocity = velocity_cartesian[k][second_interior_cell][i];
                                
                                noslip(user, distance_to_second_cell, distance_to_first_cell,
                                      wall_velocity, reference_velocity,
                                      &velocity_cartesian[k][first_interior_cell][i],
                                      wall_normal[0], wall_normal[1], wall_normal[2]);
                                
                                wall_function_loglaw(user, wall_roughness_height,
                                                    distance_to_second_cell, distance_to_first_cell,
                                                    wall_velocity, reference_velocity,
                                                    &velocity_cartesian[k][first_interior_cell][i],
                                                    &friction_velocity[k][first_interior_cell][i],
                                                    wall_normal[0], wall_normal[1], wall_normal[2]);
                                
                                velocity_boundary[k][ghost_cell_index][i].x = 0.0;
                                velocity_boundary[k][ghost_cell_index][i].y = 0.0;
                                velocity_boundary[k][ghost_cell_index][i].z = 0.0;
                                velocity_contravariant[k][first_interior_cell][i].y = 0.0;
                            }
                        }
                    }
                }
            } break;
            
            // =================================================================
            // NEGATIVE Z FACE (k = 0, normal points in +Z direction)
            // =================================================================
            case BC_FACE_NEG_Z: {
                if (grid_start_k == 0) {
                    const PetscInt ghost_cell_index = grid_start_k;
                    const PetscInt first_interior_cell = grid_start_k + 1;
                    const PetscInt second_interior_cell = grid_start_k + 2;
                    
                    for (PetscInt j = loop_start_j; j < loop_end_j; j++) {
                        for (PetscInt i = loop_start_i; i < loop_end_i; i++) {
                            
                            if (node_vertex_flag[first_interior_cell][j][i] < 0.1) {
                                
                                PetscReal face_area = sqrt(
                                    zet[ghost_cell_index][j][i].x * zet[ghost_cell_index][j][i].x + 
                                    zet[ghost_cell_index][j][i].y * zet[ghost_cell_index][j][i].y + 
                                    zet[ghost_cell_index][j][i].z * zet[ghost_cell_index][j][i].z
                                );
                                
                                PetscReal distance_to_first_cell = 0.5 / cell_jacobian[first_interior_cell][j][i] / face_area;
                                PetscReal distance_to_second_cell = 2.0 * distance_to_first_cell + 
                                                                     0.5 / cell_jacobian[second_interior_cell][j][i] / face_area;
                                
                                PetscReal wall_normal[3];
                                wall_normal[0] = zet[ghost_cell_index][j][i].x / face_area;
                                wall_normal[1] = zet[ghost_cell_index][j][i].y / face_area;
                                wall_normal[2] = zet[ghost_cell_index][j][i].z / face_area;
                                
                                Cmpnts wall_velocity, reference_velocity;
                                wall_velocity.x = wall_velocity.y = wall_velocity.z = 0.0;
                                reference_velocity = velocity_cartesian[second_interior_cell][j][i];
                                
                                noslip(user, distance_to_second_cell, distance_to_first_cell,
                                      wall_velocity, reference_velocity,
                                      &velocity_cartesian[first_interior_cell][j][i],
                                      wall_normal[0], wall_normal[1], wall_normal[2]);
                                
                                wall_function_loglaw(user, wall_roughness_height,
                                                    distance_to_second_cell, distance_to_first_cell,
                                                    wall_velocity, reference_velocity,
                                                    &velocity_cartesian[first_interior_cell][j][i],
                                                    &friction_velocity[first_interior_cell][j][i],
                                                    wall_normal[0], wall_normal[1], wall_normal[2]);
                                
                                velocity_boundary[ghost_cell_index][j][i].x = 0.0;
                                velocity_boundary[ghost_cell_index][j][i].y = 0.0;
                                velocity_boundary[ghost_cell_index][j][i].z = 0.0;
                                velocity_contravariant[ghost_cell_index][j][i].z = 0.0;
                            }
                        }
                    }
                }
            } break;
            
            // =================================================================
            // POSITIVE Z FACE (k = mz-1, normal points in -Z direction)
            // =================================================================
            case BC_FACE_POS_Z: {
                if (grid_end_k == grid_size_k) {
                    const PetscInt ghost_cell_index = grid_end_k - 1;
                    const PetscInt first_interior_cell = grid_end_k - 2;
                    const PetscInt second_interior_cell = grid_end_k - 3;
                    
                    for (PetscInt j = loop_start_j; j < loop_end_j; j++) {
                        for (PetscInt i = loop_start_i; i < loop_end_i; i++) {
                            
                            if (node_vertex_flag[first_interior_cell][j][i] < 0.1) {
                                
                                PetscReal face_area = sqrt(
                                    zet[first_interior_cell][j][i].x * zet[first_interior_cell][j][i].x + 
                                    zet[first_interior_cell][j][i].y * zet[first_interior_cell][j][i].y + 
                                    zet[first_interior_cell][j][i].z * zet[first_interior_cell][j][i].z
                                );
                                
                                PetscReal distance_to_first_cell = 0.5 / cell_jacobian[first_interior_cell][j][i] / face_area;
                                PetscReal distance_to_second_cell = 2.0 * distance_to_first_cell + 
                                                                     0.5 / cell_jacobian[second_interior_cell][j][i] / face_area;
                                
                                PetscReal wall_normal[3];
                                wall_normal[0] = -zet[first_interior_cell][j][i].x / face_area;
                                wall_normal[1] = -zet[first_interior_cell][j][i].y / face_area;
                                wall_normal[2] = -zet[first_interior_cell][j][i].z / face_area;
                                
                                Cmpnts wall_velocity, reference_velocity;
                                wall_velocity.x = wall_velocity.y = wall_velocity.z = 0.0;
                                reference_velocity = velocity_cartesian[second_interior_cell][j][i];
                                
                                noslip(user, distance_to_second_cell, distance_to_first_cell,
                                      wall_velocity, reference_velocity,
                                      &velocity_cartesian[first_interior_cell][j][i],
                                      wall_normal[0], wall_normal[1], wall_normal[2]);
                                
                                wall_function_loglaw(user, wall_roughness_height,
                                                    distance_to_second_cell, distance_to_first_cell,
                                                    wall_velocity, reference_velocity,
                                                    &velocity_cartesian[first_interior_cell][j][i],
                                                    &friction_velocity[first_interior_cell][j][i],
                                                    wall_normal[0], wall_normal[1], wall_normal[2]);
                                
                                velocity_boundary[ghost_cell_index][j][i].x = 0.0;
                                velocity_boundary[ghost_cell_index][j][i].y = 0.0;
                                velocity_boundary[ghost_cell_index][j][i].z = 0.0;
                                velocity_contravariant[first_interior_cell][j][i].z = 0.0;
                            }
                        }
                    }
                }
            } break;
        }
    }
    
    // =========================================================================
    // STEP 4: Restore all arrays and release memory
    // =========================================================================
    ierr = DMDAVecRestoreArray(user->fda, user->Ucat, &velocity_cartesian); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Ucont, &velocity_contravariant); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Bcs.Ubcs, &velocity_boundary); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi, (const Cmpnts***)&csi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta, (const Cmpnts***)&eta); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet, (const Cmpnts***)&zet); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da, user->lNvert, (const PetscReal***)&node_vertex_flag); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da, user->lAj, (const PetscReal***)&cell_jacobian); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, user->lFriction_Velocity, &friction_velocity); CHKERRQ(ierr);
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Complete.\n");
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ApplyBoundaryConditions"
/**
 * @brief Main master function to apply all boundary conditions for a time step.
 *
 * This function orchestrates the entire boundary condition workflow in a specific,
 * dependency-aware order:
 *
 * 1.  **Iterate Non-Periodic BCs:** It then enters an iterative loop to solve for
 *     the non-periodic boundary conditions. This allows complex, coupled conditions
 *     (like mass-conserving outlets that depend on inlet fluxes) to converge.
 *     - `BoundarySystem_ExecuteStep` calculates fluxes and sets boundary values (`ubcs`).
 *     - `UpdateDummyCells` uses these values to update the ghost cells.
 *
 * 2.  **Apply Periodic BCs:** First, it handles all periodic boundaries. This is a
 *     direct, non-iterative data copy that establishes the "wrap-around" state for
 *     all relevant fields. This provides a fixed constraint for the subsequent steps.
 * 
 * 3.  **Update Corner Nodes:** After the loop, `UpdateCornerNodes` is called to
 *     resolve the values at ghost cell edges and corners, using the now-final
 *     values from both the periodic and non-periodic faces.
 *
 * 4.  **Final Ghost Synchronization:** Finally, it performs a global-to-local
 *     update on all key fields to ensure that every processor's ghost-cell data
 *     is fully consistent before the main solver proceeds.
 *
 * @param user The main UserCtx struct containing the complete simulation state.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ApplyBoundaryConditions(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL,LOG_TRACE,"Boundary Condition Application begins.\n");

    // STEP 1: Main iteration loop for applying and converging non-periodic BCs.
    // The number of iterations (e.g., 3) allows information to propagate
    // between coupled boundaries, like an inlet and a conserving outlet.
    for (PetscInt iter = 0; iter < 3; iter++) {
        // (a) Execute the boundary system. This phase calculates fluxes across
        //     the domain and then applies the physical logic for each non-periodic
        //     handler, setting the `ubcs` (boundary value) array.
        ierr = BoundarySystem_ExecuteStep(user); CHKERRQ(ierr);

        LOG_ALLOW(GLOBAL,LOG_VERBOSE,"Boundary Condition Setup Executed.\n");

        // (b) Synchronize the updated ghost cells across all processors to ensure
        //     all ucont values are current before updating the dummy cells.
        ierr = UpdateLocalGhosts(user, "Ucont"); CHKERRQ(ierr);

        // (c) Convert updated Contravariant velocities to Cartesian velocities.
        ierr = Contra2Cart(user); CHKERRQ(ierr);

        // (d) Synchronize the updated Cartesian velocities across all processors
        //     to ensure all ucat values are current before updating the dummy cells.
        ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);

        // (e) If Wall functions are enabled, apply them now to adjust near-wall velocities.
        if(user->simCtx->wallfunction){
          // Apply wall function adjustments to the boundary velocities.
          ierr = ApplyWallFunction(user); CHKERRQ(ierr);

          // Synchronize the updated Cartesian velocities after wall function adjustments.
          ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);

          LOG_ALLOW(GLOBAL,LOG_VERBOSE,"Wall Function Applied at Walls.\n");
        }

        // (f) Update the first layer of ghost cells for non-periodic faces using
        //     the newly computed `ubcs` values.
        ierr = UpdateDummyCells(user); CHKERRQ(ierr);

        LOG_ALLOW(GLOBAL,LOG_VERBOSE,"Dummy Cells/Ghost Cells Updated.\n");

        // (g) Synchronize the updated dummy cells across all processors to ensure
        //     consistency before the next iteration or finalization.
        ierr = UpdateLocalGhosts(user, "P"); CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, "Ucont"); CHKERRQ(ierr);

    }

    // STEP 2: Handle all periodic boundaries. This is a one-time direct copy
    // that sets the absolute constraints for the rest of the solve.
    ierr = ApplyPeriodicBCs(user); CHKERRQ(ierr);


    // STEP 3: Update the corner and edge ghost nodes. This routine calculates
    // values for corners/edges by averaging their neighbors, which have been
    // finalized in the steps above (both periodic and non-periodic).
    ierr = UpdateCornerNodes(user); CHKERRQ(ierr);

    // STEP 4: Final ghost node synchronization. This ensures all changes made
    // to the global vectors are reflected in the local ghost regions of all
    // processors, making the state fully consistent before the next solver stage.
    ierr = UpdateLocalGhosts(user, "P"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "Ucat"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "Ucont"); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

