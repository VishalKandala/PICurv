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


      LOG_ALLOW(LOCAL, LOG_DEBUG,
      "[Rank %d] Check Service for Inlet %s:\n"
      "    - Local Domain: starts at cell (%d,%d,%d), has (%d,%d,%d) cells.\n"
      "    - Global Domain: has (%d,%d,%d) nodes, so last cell is (%d,%d,%d).\n"
      "    - Decision: %s\n",
      rank_for_logging,
      BCFaceToString((BCFace)user->identifiedInletBCFace),
      owned_start_cell_i, owned_start_cell_j, owned_start_cell_k,
      num_owned_cells_on_rank_i, num_owned_cells_on_rank_j, num_owned_cells_on_rank_k,
      IM_nodes_global, JM_nodes_global, KM_nodes_global,
      last_global_cell_idx_i, last_global_cell_idx_j, last_global_cell_idx_k,
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

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d check for face %d: Result=%s. \n",
        rank_for_logging, face_id, (*can_service_out ? "TRUE" : "FALSE"));

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

    LOG_ALLOW(LOCAL, LOG_DEBUG, "[Rank %d] Distributing %lld particles over %d lines on face %s.\n",
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

    LOG_ALLOW(LOCAL, LOG_DEBUG,
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

    LOG_ALLOW(LOCAL, LOG_DEBUG,
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
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Inlet face %s. Owned cells(i,j,k):(%d,%d,%d). GlobNodes(I,J,K):(%d,%d,%d). Rank's DA node starts at (%d,%d,%d).\n",
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
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: Target Cell Node =(%d,%d,%d). (xi,et,zt)=(%.2e,%.2f,%.2f). \n",
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
      LOG_ALLOW(GLOBAL,LOG_DEBUG," for face %s(%d), legacy type = %d & modern type = %s .\n",face_str,i,user->bctype[i],bc_type_str);
    }
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"    -> Translation complete.\n");
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "BoundarySystem_ExecuteStep_Legacy"

/**
 * @brief Acts as a temporary bridge to the legacy boundary condition implementation.
 *
 * This function is the key to our integration strategy. It matches the signature
 * of the modern `BoundarySystem_ExecuteStep` function that SetEulerianFields
 * expects to call.
 *
 * However, instead of containing new handler-based logic, it simply calls the
 * monolithic legacy `FormBCS` function. This allows the modern orchestrator to
 * drive the old solver logic without modification.
 *
 * @param user The UserCtx for a single block.
 * @return PetscErrorCode
 */
PetscErrorCode BoundarySystem_ExecuteStep_Legacy(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    // The sole purpose of this function is to call the old logic for the
    // specific block context that was passed in.
    ierr = FormBCS(user); CHKERRQ(ierr);

    // NOTE: The legacy `main` called Block_Interface_U after the FormBCS loop.
    // We will handle that at a higher level (in our AdvanceSimulation loop)
    // after all blocks have had their BCs applied for the step.
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "InflowFlux"

/**
 * @brief Applies inlet boundary conditions based on the modern BC handling system.
 *
 * This function iterates through all 6 domain faces. For each face identified as an
 * INLET, it applies the velocity profile specified by its assigned handler and
 * parameters (e.g., 'constant_velocity' with vx,vy,vz or 'parabolic' with u_max).
 *
 * It calculates the contravariant flux (Ucont), Cartesian velocity on the face (Ubcs),
 * and the staggered Cartesian velocity (Ucat). It also computes the total incoming
 * flux and area across all MPI ranks.
 *
 * @param user The main UserCtx struct containing the BC configuration and PETSc objects.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode InflowFlux(UserCtx *user) 
{
  PetscErrorCode ierr;
  PetscReal    lFluxIn = 0.0, lAreaIn = 0.0, AreaSumIn;
  Vec          lCoor;
  Cmpnts       ***ucont, ***ubcs, ***ucat, ***coor, ***csi, ***eta, ***zet, ***cent;  
  PetscReal    ***nvert;
  
  DM            da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  // Get context for coordinate transformation if needed by a handler
  SimCtx      *simCtx = user->simCtx;
  PetscReal   CMx_c = simCtx->CMx_c;
  PetscReal   CMy_c = simCtx->CMy_c;
  PetscReal   CMz_c = simCtx->CMz_c;
  
  PetscFunctionBeginUser;

  PROFILE_FUNCTION_BEGIN;

  lxs = xs; lxe = xe; lys = ys; lye = ye; lzs = zs; lze = ze;
  if (xs==0) lxs = xs+1; if (ys==0) lys = ys+1; if (zs==0) lzs = zs+1;
  if (xe==mx) lxe = xe-1; if (ye==my) lye = ye-1; if (ze==mz) lze = ze-1;

  // --- Get PETSc arrays ---
  DMGetCoordinatesLocal(da,&lCoor);
  DMDAVecGetArray(fda, lCoor, &coor); // Use local coordinates
  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecGetArray(fda, user->Ucat,  &ucat);
  DMDAVecGetArray(da,  user->lNvert, &nvert);
  DMDAVecGetArray(fda, user->Cent, &cent);
  DMDAVecGetArray(fda, user->lCsi,  &csi);
  DMDAVecGetArray(fda, user->lEta,  &eta);
  DMDAVecGetArray(fda, user->lZet,  &zet);

  // --- Main Loop over all 6 faces ---
  for (PetscInt fn=0; fn<6; fn++) {
    BoundaryFaceConfig *face_config = &user->boundary_faces[fn];

    if (face_config->mathematical_type != INLET) {
        continue; // Skip non-inlet faces
    }
    
    // This processor only acts if it is on the boundary of the global domain
    PetscBool is_on_boundary = ( (fn==0 && xs==0) || (fn==1 && xe==mx) ||
                                 (fn==2 && ys==0) || (fn==3 && ye==my) ||
                                 (fn==4 && zs==0) || (fn==5 && ze==mz) );
    if (!is_on_boundary) continue;

    LOG_ALLOW(GLOBAL,LOG_DEBUG,"Applying INLET handler for face: %s \n", BCFaceToString(face_config->face_id));

    // --- Loop over the specific face geometry ---
    switch(face_config->face_id) {
      case BC_FACE_NEG_X: // -Xi
      case BC_FACE_POS_X: // +Xi
      {
        PetscReal sign = (face_config->face_id == BC_FACE_NEG_X) ? 1.0 : -1.0;
        PetscInt i = (face_config->face_id == BC_FACE_NEG_X) ? xs : mx - 2;
        for (PetscInt k=lzs; k<lze; k++) {
          for (PetscInt j=lys; j<lye; j++) {
            if ( (sign > 0 && nvert[k][j][i+1] > 0.1) || (sign < 0 && nvert[k][j][i] > 0.1) ) continue;

            PetscReal uin_this_point = 0.0;
            // --- Determine velocity based on the handler for this point ---
            if (face_config->handler_type == BC_HANDLER_INLET_CONSTANT_VELOCITY) {
                PetscBool found;
                ierr = GetBCParamReal(face_config->params, "vx", &uin_this_point, &found); CHKERRQ(ierr);
            } else if(face_config->handler_type== BC_HANDLER_INLET_PARABOLIC){
              PetscBool found;
              PetscReal umax,diameter=1.0;
              ierr = GetBCParamReal(face_config->params,"u_max",&umax,&found); CHKERRQ(ierr);
              ierr = GetBCParamReal(face_config->params,"diameter",&diameter,&found); CHKERRQ(ierr);

              // Radius is in the YZ-plane for an X-face inlet
              PetscReal yc = cent[k][j][i + (sign>0)].y - CMy_c;
              PetscReal zc = cent[k][j][i + (sign>0)].z - CMz_c;
              PetscReal r = sqrt(yc*yc + zc*zc);
              PetscReal r_norm = 2.0 * r / diameter;
              uin_this_point = umax * (1.0 - r_norm * r_norm);
              if (r_norm > 1.0) uin_this_point = 0.0; 
            }
            // Add other X-face handlers like 'else if (handler == ...)' here

            // --- Apply the calculated velocity ---
            PetscReal CellArea = sqrt(csi[k][j][i].z*csi[k][j][i].z + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].x*csi[k][j][i].x);
            lAreaIn += CellArea;
            ucont[k][j][i].x  = sign * uin_this_point * CellArea;
            lFluxIn          += ucont[k][j][i].x;
            ubcs[k][j][i + (sign < 0)].x = sign * uin_this_point * csi[k][j][i].x / CellArea;
            ubcs[k][j][i + (sign < 0)].y = sign * uin_this_point * csi[k][j][i].y / CellArea;
            ubcs[k][j][i + (sign < 0)].z = sign * uin_this_point * csi[k][j][i].z / CellArea;
            ucat[k][j][i + (sign > 0)] = ubcs[k][j][i + (sign < 0)];
          }
        }
      } break;

      case BC_FACE_NEG_Y: // -Eta
      case BC_FACE_POS_Y: // +Eta
      {
        PetscReal sign = (face_config->face_id == BC_FACE_NEG_Y) ? 1.0 : -1.0;
        PetscInt j = (face_config->face_id == BC_FACE_NEG_Y) ? ys : my - 2;
        for (PetscInt k=lzs; k<lze; k++) {
          for (PetscInt i=lxs; i<lxe; i++) {
            if ( (sign > 0 && nvert[k][j+1][i] > 0.1) || (sign < 0 && nvert[k][j][i] > 0.1) ) continue;

            PetscReal uin_this_point = 0.0;
            if (face_config->handler_type == BC_HANDLER_INLET_CONSTANT_VELOCITY) {
                PetscBool found;
                ierr = GetBCParamReal(face_config->params, "vy", &uin_this_point, &found); CHKERRQ(ierr);
            }else if(face_config->handler_type == BC_HANDLER_INLET_PARABOLIC){
              PetscBool found;
              PetscReal umax,diameter=1.0;
              ierr = GetBCParamReal(face_config->params,"umax",&umax,&found); CHKERRQ(ierr);
              ierr = GetBCParamReal(face_config->params,"diameter",&diameter,&found); CHKERRQ(ierr);
              
              // Radius is in the XZ-plane for a Y-face inlet
              PetscReal xc = cent[k][j + (sign>0)][i].x - CMx_c;
              PetscReal zc = cent[k][j + (sign>0)][i].z - CMz_c;
              PetscReal r = sqrt(xc*xc + zc*zc);
              PetscReal r_norm = 2.0 * r / diameter;
              uin_this_point = umax * (1.0 - r_norm * r_norm);
              if (r_norm > 1.0) uin_this_point = 0.0;

            }

            PetscReal CellArea = sqrt(eta[k][j][i].z*eta[k][j][i].z + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].x*eta[k][j][i].x);
            lAreaIn += CellArea;
            ucont[k][j][i].y  = sign * uin_this_point * CellArea;
            lFluxIn          += ucont[k][j][i].y;
            ubcs[k][j + (sign < 0)][i].x = sign * uin_this_point * eta[k][j][i].x / CellArea;
            ubcs[k][j + (sign < 0)][i].y = sign * uin_this_point * eta[k][j][i].y / CellArea;
            ubcs[k][j + (sign < 0)][i].z = sign * uin_this_point * eta[k][j][i].z / CellArea;
            ucat[k][j + (sign > 0)][i] = ubcs[k][j + (sign < 0)][i];
          }
        }
      } break;

      case BC_FACE_NEG_Z: // -Zeta
      case BC_FACE_POS_Z: // +Zeta
      {
        PetscReal sign = (face_config->face_id == BC_FACE_NEG_Z) ? 1.0 : -1.0;
        PetscInt k = (face_config->face_id == BC_FACE_NEG_Z) ? zs : mz - 2;
        for (PetscInt j=lys; j<lye; j++) {
          for (PetscInt i=lxs; i<lxe; i++) {
            if ( (sign > 0 && nvert[k+1][j][i] > 0.1) || (sign < 0 && nvert[k][j][i] > 0.1) ) continue;

            PetscReal uin_this_point = 0.0;
            if (face_config->handler_type == BC_HANDLER_INLET_CONSTANT_VELOCITY) {
                PetscBool found;
                ierr = GetBCParamReal(face_config->params, "vz", &uin_this_point, &found); CHKERRQ(ierr);
            } else if (face_config->handler_type == BC_HANDLER_INLET_PARABOLIC) {
                PetscBool found;
                PetscReal umax, diameter=1.0; // Default diameter for r_norm = 2*r
                ierr = GetBCParamReal(face_config->params, "u_max", &umax, &found); CHKERRQ(ierr);
                ierr = GetBCParamReal(face_config->params, "diameter", &diameter, &found); CHKERRQ(ierr); // Optional
                
                // Radius in the XY-plane for a Z-face inlet.
                PetscReal xc = cent[k + (sign>0)][j][i].x - CMx_c;
                PetscReal yc = cent[k + (sign>0)][j][i].y - CMy_c;
                PetscReal r = sqrt(xc*xc + yc*yc);
                PetscReal r_norm = 2.0 * r / diameter; // Normalized radius
                uin_this_point = umax * (1.0 - r_norm * r_norm);
                if (r_norm > 1.0) uin_this_point = 0.0;
            }

            PetscReal CellArea = sqrt(zet[k][j][i].z*zet[k][j][i].z + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].x*zet[k][j][i].x);
            lAreaIn += CellArea;
            ucont[k][j][i].z  = sign * uin_this_point * CellArea;
            lFluxIn          += ucont[k][j][i].z;
            ubcs[k + (sign < 0)][j][i].x = sign * uin_this_point * zet[k][j][i].x / CellArea;
            ubcs[k + (sign < 0)][j][i].y = sign * uin_this_point * zet[k][j][i].y / CellArea;
            ubcs[k + (sign < 0)][j][i].z = sign * uin_this_point * zet[k][j][i].z / CellArea;
            ucat[k + (sign > 0)][j][i] = ubcs[k + (sign < 0)][j][i];
          }
        }
      } break;
    } // end switch(face_id)
  } // end for(fn)

  // --- Finalize: Sum flux and area from all processes ---
  ierr = MPI_Allreduce(&lFluxIn, &simCtx->FluxInSum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
  ierr = MPI_Allreduce(&lAreaIn, &AreaSumIn, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);
  
  LOG_ALLOW(GLOBAL,LOG_DEBUG,"Inflow: Flux - Area:  %le - %le \n", simCtx->FluxInSum, AreaSumIn);    
  
  // --- Restore PETSc arrays ---
  DMDAVecRestoreArray(fda, lCoor, &coor);
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  DMDAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  DMDAVecRestoreArray(da,  user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->Cent, &cent);
  DMDAVecRestoreArray(fda, user->lCsi,  &csi);
  DMDAVecRestoreArray(fda, user->lEta,  &eta);
  DMDAVecRestoreArray(fda, user->lZet,  &zet);

  // --- Update local vectors for subsequent computations ---
  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  
  PROFILE_FUNCTION_END;

  PetscFunctionReturn(0); 
}

#undef __FUNCT__
#define __FUNCT__ "OutflowFlux"

/**
 * @brief Calculates the total outgoing flux through all OUTLET faces for reporting.
 *
 * NOTE: In a mixed modern/legacy environment, this function is for DIAGNOSTICS ONLY.
 * It reads the contravariant velocities and calculates the total flux passing through
 * faces marked as OUTLET. It does NOT apply any boundary conditions itself, as that
 * is still the responsibility of the legacy FormBCS function.
 *
 * @param user The main UserCtx struct containing BC config and PETSc vectors.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode OutflowFlux(UserCtx *user) {
  
    PetscErrorCode ierr;
    PetscReal      lFluxOut = 0.0;
    Cmpnts         ***ucont;

    DM             fda = user->fda;
    DMDALocalInfo  info = user->info;
    PetscInt       xs = info.xs, xe = info.xs + info.xm;
    PetscInt       ys = info.ys, ye = info.ys + info.ym;
    PetscInt       zs = info.zs, ze = info.zs + info.zm;
    PetscInt       mx = info.mx, my = info.my, mz = info.mz;

    PetscFunctionBeginUser;
    
    PROFILE_FUNCTION_BEGIN;

    ierr = DMDAVecGetArrayRead(fda, user->Ucont, &ucont); CHKERRQ(ierr);
  
    // --- Loop over all 6 faces to find OUTLETS ---
    for (PetscInt fn = 0; fn < 6; fn++) {
        if (user->boundary_faces[fn].mathematical_type != OUTLET) {
            continue; 
        }

        PetscBool is_on_boundary = ( (fn==0 && xs==0) || (fn==1 && xe==mx) ||
                                     (fn==2 && ys==0) || (fn==3 && ye==my) ||
                                     (fn==4 && zs==0) || (fn==5 && ze==mz) );
        if (!is_on_boundary) continue;

        // --- Sum the flux for the appropriate face and component ---
        switch ((BCFace)fn) {
            case BC_FACE_NEG_X: case BC_FACE_POS_X: {
                PetscInt i = (fn == 0) ? xs : mx - 2;
                for (PetscInt k=info.zs; k<info.zs+info.zm; k++) for (PetscInt j=info.ys; j<info.ys+info.ym; j++) {
                    lFluxOut += ucont[k][j][i].x;
                }
            } break;

            case BC_FACE_NEG_Y: case BC_FACE_POS_Y: {
                PetscInt j = (fn == 2) ? ys : my - 2;
                for (PetscInt k=info.zs; k<info.zs+info.zm; k++) for (PetscInt i=info.xs; i<info.xs+info.xm; i++) {
                    lFluxOut += ucont[k][j][i].y;
                }
            } break;

            case BC_FACE_NEG_Z: case BC_FACE_POS_Z: {
                PetscInt k = (fn == 4) ? zs : mz - 2;
                for (PetscInt j=info.ys; j<info.ys+info.ym; j++) for (PetscInt i=info.xs; i<info.xs+info.xm; i++) {
                    lFluxOut += ucont[k][j][i].z;
                }
            } break;
        } // end switch
    } // end for loop

    ierr = DMDAVecRestoreArrayRead(fda, user->Ucont, &ucont); CHKERRQ(ierr);

    // --- Finalize: Sum and store the global total flux ---
    ierr = MPI_Allreduce(&lFluxOut, &user->simCtx->FluxOutSum, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Reported Global FluxOutSum = %.6f\n", user->simCtx->FluxOutSum);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FormBCS"

/* Boundary condition defination (array user->bctype[0-5]):
   0:	interpolation/interface
   -1:  wallfunction
   1:	solid wall (not moving)
   2:	moving solid wall (U=1)
   3:   slip wall/symmetry
   5:	Inlet
   4:	Outlet
   6:   farfield
   7:   periodic
   8:   Characteristic BC
   9:   Analytical Vortex
   10:  Oulet Junction
   11:  Annulus
   12:  Ogrid
   13:  Rheology
   14:  Outlet with Interface
   15:  No Gradient (Similar to Farfield)  
*/

PetscErrorCode FormBCS(UserCtx *user)
{
  DM            da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  PetscInt	i, j, k;

  PetscReal	***nvert,***lnvert; //local working array

  PetscReal	***p,***lp;
  Cmpnts	***ucont, ***ubcs, ***ucat,***lucat, ***csi, ***eta, ***zet;
  Cmpnts	***cent,***centx,***centy,***centz,***coor;
  PetscScalar	FluxIn, FluxOut,ratio;
  PetscScalar   lArea, AreaSum;
 
  PetscScalar   FarFluxIn=0., FarFluxOut=0., FarFluxInSum, FarFluxOutSum;
  PetscScalar   FarAreaIn=0., FarAreaOut=0., FarAreaInSum, FarAreaOutSum;
  PetscScalar   FluxDiff, VelDiffIn, VelDiffOut;
  Cmpnts        V_frame;
 
  PetscReal Un, nx,ny,nz,A;

  SimCtx *simCtx = user->simCtx;  

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  PetscInt	gxs, gxe, gys, gye, gzs, gze;

  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx ) lxe = xe-1;
  if (ye==my ) lye = ye-1;
  if (ze==mz ) lze = ze-1;

  PetscFunctionBeginUser;

  PROFILE_FUNCTION_BEGIN;

  DMDAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);

  DMDAVecGetArray(fda, user->lCsi,  &csi);
  DMDAVecGetArray(fda, user->lEta,  &eta);
  DMDAVecGetArray(fda, user->lZet,  &zet);

  PetscInt ttemp;
  for (ttemp=0; ttemp<3; ttemp++) {
    DMDAVecGetArray(da, user->Nvert, &nvert); 
    DMDAVecGetArray(fda, user->lUcat,  &ucat);
    DMDAVecGetArray(fda, user->Ucont, &ucont);
/* ==================================================================================             */
/*   FAR-FIELD BC */
/* ==================================================================================             */
 

  // reset FAR FLUXES
  FarFluxIn = 0.; FarFluxOut=0.;
  FarAreaIn = 0.; FarAreaOut=0.;

  PetscReal lFlux_abs=0.0,FluxSum_abs=0.0,ratio=0.0;

    V_frame.x=0.;
    V_frame.y=0.;
    V_frame.z=0.;

 

  if (user->bctype[0]==6) {
    if (xs == 0) {
      i= xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i].x = ucat[k][j][i+1].x;
	  ubcs[k][j][i].y = ucat[k][j][i+1].y;
	  ubcs[k][j][i].z = ucat[k][j][i+1].z;
	  ucont[k][j][i].x = ubcs[k][j][i].x * csi[k][j][i].x;
	  FarFluxIn += ucont[k][j][i].x;
	  lFlux_abs += fabs(ucont[k][j][i].x);
	  FarAreaIn += csi[k][j][i].x;
	}
      }
    }
  }
    
  if (user->bctype[1]==6) {
    if (xe==mx) {
      i= xe-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ubcs[k][j][i].x = ucat[k][j][i-1].x;
	  ubcs[k][j][i].y = ucat[k][j][i-1].y;
	  ubcs[k][j][i].z = ucat[k][j][i-1].z;
	  ucont[k][j][i-1].x = ubcs[k][j][i].x * csi[k][j][i-1].x;
	  FarFluxOut += ucont[k][j][i-1].x;
	  lFlux_abs  += fabs(ucont[k][j][i-1].x);
	  FarAreaOut += csi[k][j][i-1].x;
	}
      }
    }
  }

  if (user->bctype[2]==6) {
    if (ys==0) {
      j= ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k][j+1][i].x;
	  ubcs[k][j][i].y = ucat[k][j+1][i].y;
	  ubcs[k][j][i].z = ucat[k][j+1][i].z;
	  ucont[k][j][i].y = ubcs[k][j][i].y * eta[k][j][i].y;
	  FarFluxIn += ucont[k][j][i].y;
	  lFlux_abs += fabs(ucont[k][j][i].y);
	  FarAreaIn += eta[k][j][i].y;
	}
      }
    }
  }
  
  if (user->bctype[3]==6) {
    if (ye==my) {
      j=ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k][j-1][i].x;
	  ubcs[k][j][i].y = ucat[k][j-1][i].y;
	  ubcs[k][j][i].z = ucat[k][j-1][i].z;
	  ucont[k][j-1][i].y = ubcs[k][j][i].y * eta[k][j-1][i].y;
	  FarFluxOut += ucont[k][j-1][i].y;
	  lFlux_abs  += fabs(ucont[k][j-1][i].y);
	  FarAreaOut += eta[k][j-1][i].y;
	}
      }
    }
  }

  if (user->bctype[4]==6 || user->bctype[4]==10) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k+1][j][i].x;
	  ubcs[k][j][i].y = ucat[k+1][j][i].y;
	  ubcs[k][j][i].z = ucat[k+1][j][i].z;
	  ucont[k][j][i].z = ubcs[k][j][i].z * zet[k][j][i].z;
	  FarFluxIn += ucont[k][j][i].z;
	  lFlux_abs += fabs(ucont[k][j][i].z);
	  FarAreaIn += zet[k][j][i].z;
	}
      }
    }
  }

  if (user->bctype[5]==6 || user->bctype[5]==10) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = ucat[k-1][j][i].z;
	  ucont[k-1][j][i].z = ubcs[k][j][i].z * zet[k-1][j][i].z;
	  FarFluxOut += ucont[k-1][j][i].z;
	  lFlux_abs  += fabs(ucont[k-1][j][i].z); 
	  FarAreaOut += zet[k-1][j][i].z;
	}
      }
    }
  }
  
  MPI_Allreduce(&FarFluxIn,&FarFluxInSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&FarFluxOut,&FarFluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&lFlux_abs,&FluxSum_abs,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); 
  MPI_Allreduce(&FarAreaIn,&FarAreaInSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  MPI_Allreduce(&FarAreaOut,&FarAreaOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
 
  if (user->bctype[5]==6 || user->bctype[3]==6 || user->bctype[1]==6) {
  
    ratio=(FarFluxInSum - FarFluxOutSum)/FluxSum_abs;
    if (fabs(FluxSum_abs) <1.e-10) ratio = 0.;
    
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "/FluxSum_abs %le ratio %le \n", FluxSum_abs,ratio);
    
    FluxDiff = 0.5*(FarFluxInSum - FarFluxOutSum) ;
    VelDiffIn  = FluxDiff / FarAreaInSum ;
    
    if (fabs(FarAreaInSum) <1.e-6) VelDiffIn = 0.;

    VelDiffOut  = FluxDiff / FarAreaOutSum ;
  
    if (fabs(FarAreaOutSum) <1.e-6) VelDiffOut = 0.;

    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Far Flux Diff %d %le %le %le %le %le %le %le\n", simCtx->step, FarFluxInSum, FarFluxOutSum, FluxDiff, FarAreaInSum, FarAreaOutSum, VelDiffIn, VelDiffOut);
           
  }
  
  if (user->bctype[5]==10) {
    FluxDiff = simCtx->FluxInSum -( FarFluxOutSum -FarFluxInSum) ;
    VelDiffIn  = 1/3.*FluxDiff / (FarAreaInSum);// +  FarsimCtx->AreaOutSum);
    if (fabs(FarAreaInSum) <1.e-6) VelDiffIn = 0.;

    VelDiffOut  = 2./3.* FluxDiff / (FarAreaOutSum) ;//(FarAreaInSum +  FarsimCtx->AreaOutSum) ;
   
    if (fabs(FarAreaOutSum) <1.e-6) VelDiffOut = 0.;

    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Far Flux Diff %d %le %le %le %le %le %le %le\n", simCtx->step, FarFluxInSum, FarFluxOutSum, FluxDiff, FarAreaInSum, FarAreaOutSum, VelDiffIn, VelDiffOut);
           
  }
  
  
  // scale global mass conservation

  if (user->bctype[5]==6 || user->bctype[5]==10) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ucont[k-1][j][i].z  += ratio*fabs(ucont[k-1][j][i].z);
	  ubcs[k][j][i].z = ucont[k-1][j][i].z/zet[k-1][j][i].z;
	  //  ubcs[k][j][i].z = ucat[k-1][j][i].z + VelDiffOut ;//+ V_frame.z;
	  // ucont[k-1][j][i].z = ubcs[k][j][i].z * zet[k-1][j][i].z;
	}
      }
    }
  }

  if (user->bctype[3]==6) {
    if (ye==my) {
      j=ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {

	  ucont[k][j-1][i].y +=ratio*fabs(ucont[k][j-1][i].y);
	  ubcs[k][j][i].y = ucont[k][j-1][i].y /eta[k][j-1][i].y;
	  //	  ubcs[k][j][i].y = ucat[k][j-1][i].y + VelDiffOut;// + V_frame.y;
	  // ucont[k][j-1][i].y = ubcs[k][j][i].y * eta[k][j-1][i].y;

	}
      }
    }
  }
    
  if (user->bctype[1]==6) {
    if (xe==mx) {
      i= xe-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ucont[k][j][i-1].x +=ratio*fabs(ucont[k][j][i-1].x);
	  ubcs[k][j][i].x = ucont[k][j][i-1].x / csi[k][j][i-1].x ;
	  //  ubcs[k][j][i].x = ucat[k][j][i-1].x + VelDiffOut;// + V_frame.x;
	  // ucont[k][j][i-1].x = ubcs[k][j][i].x * csi[k][j][i-1].x;
	}
      }
    }
  }


  if (user->bctype[0]==6) {
    if (xs == 0) {
      i= xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ucont[k][j][i].x  -=ratio*fabs(ucont[k][j][i].x);
	  ubcs[k][j][i].x = ucont[k][j][i].x / csi[k][j][i].x;
	  // ubcs[k][j][i].x = ucat[k][j][i+1].x - VelDiffIn;// + V_frame.x;
	  // ucont[k][j][i].x = ubcs[k][j][i].x * csi[k][j][i].x;
	}
      }
    }
  }
  

  if (user->bctype[2]==6) {
    if (ys==0) {
      j= ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ucont[k][j][i].y -=ratio*fabs(ucont[k][j][i].y);
	  ubcs[k][j][i].y = ucont[k][j][i].y / eta[k][j][i].y;
	  //	  ubcs[k][j][i].y = ucat[k][j+1][i].y - VelDiffIn;// + V_frame.y;
	  // ucont[k][j][i].y = ubcs[k][j][i].y * eta[k][j][i].y;
	}
      }
    }
  }
  
  
  if (user->bctype[4]==6 || user->bctype[5]==10) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ucont[k][j][i].z -=ratio*fabs(ucont[k][j][i].z);
	  ubcs[k][j][i].z =ucont[k][j][i].z / zet[k][j][i].z;
	  // ubcs[k][j][i].z = ucat[k+1][j][i].z - VelDiffIn;// + V_frame.z;
	  // ucont[k][j][i].z = ubcs[k][j][i].z * zet[k][j][i].z;

	}
      }
    }
  }

//// Amir wall Ogrid
 
if (user->bctype[2]==1 || user->bctype[2]==-1)  {
    if (ys==0) {
    j= ys;
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
        A=sqrt(eta[k][j][i].z*eta[k][j][i].z +
               eta[k][j][i].y*eta[k][j][i].y +
               eta[k][j][i].x*eta[k][j][i].x);
        nx=eta[k][j][i].x/A;
        ny=eta[k][j][i].y/A;
        nz=eta[k][j][i].z/A;
        Un=ucat[k][j+1][i].x*nx+ucat[k][j+1][i].y*ny+ucat[k][j+1][i].z*nz;
        ubcs[k][j][i].x = 0.0;
        ubcs[k][j][i].y = 0.0;
        ubcs[k][j][i].z = 0.0;
        ucont[k][j][i].y = 0.;
      }
    }
    }
 }

/* ==================================================================================             */
/*   SOLID WALL BC (NO-SLIP / NO-PENETRATION) */
/* ==================================================================================             */

// NOTE: This block is added to explicitly handle bctype=1 (solid wall) for all faces.
// It ensures both no-slip (ubcs=0) and no-penetration (ucont_normal=0).
// ubcs is handled by the implicit-zero assumption, but ucont must be set explicitly.

// -X Face (i=0)
if (user->bctype[0]==1 || user->bctype[0]==-1)  {
    if (xs==0) {
      i= xs;
      for (k=lzs; k<lze; k++) {
        for (j=lys; j<lye; j++) {
          ubcs[k][j][i].x = 0.0;
          ubcs[k][j][i].y = 0.0;
          ubcs[k][j][i].z = 0.0;
          ucont[k][j][i].x = 0.0; // Enforce no-penetration
        }
      }
    }
}

// +X Face (i=mx-1)
if (user->bctype[1]==1 || user->bctype[1]==-1)  {
    if (xe==mx) {
      i= xe-1;
      for (k=lzs; k<lze; k++) {
        for (j=lys; j<lye; j++) {
          ubcs[k][j][i].x = 0.0;
          ubcs[k][j][i].y = 0.0;
          ubcs[k][j][i].z = 0.0;
          // The relevant ucont is at the face, index i-1
          ucont[k][j][i-1].x = 0.0; // Enforce no-penetration
        }
      }
    }
}

// -Y Face (j=0)
if (user->bctype[2]==1 || user->bctype[2]==-1)  {
    if (ys==0) {
      j= ys;
      for (k=lzs; k<lze; k++) {
        for (i=lxs; i<lxe; i++) {
          ubcs[k][j][i].x = 0.0;
          ubcs[k][j][i].y = 0.0;
          ubcs[k][j][i].z = 0.0;
          ucont[k][j][i].y = 0.0; // Enforce no-penetration
        }
      }
    }
}

// +Y Face (j=my-1)
if (user->bctype[3]==1 || user->bctype[3]==-1)  {
    if (ye==my) {
      j= ye-1;
      for (k=lzs; k<lze; k++) {
        for (i=lxs; i<lxe; i++) {
          ubcs[k][j][i].x = 0.0;
          ubcs[k][j][i].y = 0.0;
          ubcs[k][j][i].z = 0.0;
          // The relevant ucont is at the face, index j-1
          ucont[k][j-1][i].y = 0.0; // Enforce no-penetration
        }
      }
    }
}

/* Original "Amir wall Ogrid" block can now be removed or commented out
   as its functionality is included above.
if (user->bctype[2]==1 || user->bctype[2]==-1)  { ... }
*/
 
/* ==================================================================================             */
/*   SYMMETRY BC */
/* ==================================================================================             */

  if (user->bctype[0]==3) {
    if (xs==0) {
    i= xs;

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	A=sqrt(csi[k][j][i].z*csi[k][j][i].z +
	       csi[k][j][i].y*csi[k][j][i].y +
	       csi[k][j][i].x*csi[k][j][i].x);
	nx=csi[k][j][i].x/A;
	ny=csi[k][j][i].y/A;
	nz=csi[k][j][i].z/A;
	Un=ucat[k][j][i+1].x*nx+ucat[k][j][i+1].y*ny+ucat[k][j][i+1].z*nz;
	ubcs[k][j][i].x = ucat[k][j][i+1].x-Un*nx;//-V_frame.x;
	ubcs[k][j][i].y = ucat[k][j][i+1].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j][i+1].z-Un*nz;
	ucont[k][j][i].x = 0.;
      }
    }
    }
  }

  if (user->bctype[1]==3) {
    if (xe==mx) {
    i= xe-1;

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	A=sqrt(csi[k][j][i-1].z*csi[k][j][i-1].z +
	       csi[k][j][i-1].y*csi[k][j][i-1].y +
	       csi[k][j][i-1].x*csi[k][j][i-1].x);
	nx=csi[k][j][i-1].x/A;
	ny=csi[k][j][i-1].y/A;
	nz=csi[k][j][i-1].z/A;
	Un=ucat[k][j][i-1].x*nx+ucat[k][j][i-1].y*ny+ucat[k][j][i-1].z*nz;
	ubcs[k][j][i].x = ucat[k][j][i-1].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j][i-1].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j][i-1].z-Un*nz;
	ucont[k][j][i-1].x = 0.;
      }
    }
    }
  }

  if (user->bctype[2]==3) {
    if (ys==0) {
    j= ys;

    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	A=sqrt(eta[k][j][i].z*eta[k][j][i].z +
	       eta[k][j][i].y*eta[k][j][i].y +
	       eta[k][j][i].x*eta[k][j][i].x);
	nx=eta[k][j][i].x/A;
	ny=eta[k][j][i].y/A;
	nz=eta[k][j][i].z/A;
	Un=ucat[k][j+1][i].x*nx+ucat[k][j+1][i].y*ny+ucat[k][j+1][i].z*nz;
	ubcs[k][j][i].x = ucat[k][j+1][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j+1][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j+1][i].z-Un*nz;
	ucont[k][j][i].y = 0.;
      }
    }
    }
  }

  if (user->bctype[3]==3) {
    if (ye==my) {
    j=ye-1;

    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	A=sqrt(eta[k][j-1][i].z*eta[k][j-1][i].z +
	       eta[k][j-1][i].y*eta[k][j-1][i].y +
	       eta[k][j-1][i].x*eta[k][j-1][i].x);
	nx=eta[k][j-1][i].x/A;
	ny=eta[k][j-1][i].y/A;
	nz=eta[k][j-1][i].z/A;
	Un=ucat[k][j-1][i].x*nx+ucat[k][j-1][i].y*ny+ucat[k][j-1][i].z*nz;
	ubcs[k][j][i].x = ucat[k][j-1][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j-1][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j-1][i].z-Un*nz;
	ucont[k][j-1][i].y = 0.;
      }
    }
    }
  }
  

  if (user->bctype[4]==3) {
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	A=sqrt(zet[k][j][i].z*zet[k][j][i].z +
	       zet[k][j][i].y*zet[k][j][i].y +
	       zet[k][j][i].x*zet[k][j][i].x);
	nx=zet[k][j][i].x/A;
	ny=zet[k][j][i].y/A;
	nz=zet[k][j][i].z/A;
	Un=ucat[k+1][j][i].x*nx+ucat[k+1][j][i].y*ny+ucat[k+1][j][i].z*nz;
	ubcs[k][j][i].x = ucat[k+1][j][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k+1][j][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k+1][j][i].z-Un*nz;
	ucont[k][j][i].z = 0.;
	}
      }
    }
  }

  if (user->bctype[5]==3) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	A=sqrt(zet[k-1][j][i].z*zet[k-1][j][i].z +
	       zet[k-1][j][i].y*zet[k-1][j][i].y +
	       zet[k-1][j][i].x*zet[k-1][j][i].x);
	nx=zet[k-1][j][i].x/A;
	ny=zet[k-1][j][i].y/A;
	nz=zet[k-1][j][i].z/A;
	Un=ucat[k-1][j][i].x*nx+ucat[k-1][j][i].y*ny+ucat[k-1][j][i].z*nz;
	ubcs[k][j][i].x = ucat[k-1][j][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k-1][j][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k-1][j][i].z-Un*nz;
	ucont[k-1][j][i].z = 0.;
	}
      }
    }
  }
 
/* ==================================================================================             */
/*     CHARACTERISTIC OUTLET BC :8 */
/* ==================================================================================             */

  if (user->bctype[5]==8) {
    if (ze == mz) {
      k = ze-2;
      FluxOut = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  FluxOut += ucont[k][j][i].z;
	}
      }
    }
    else {
      FluxOut = 0.;
    }
    
    FluxIn = simCtx->FluxInSum + FarFluxInSum;

    MPI_Allreduce(&FluxOut,&simCtx->FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    // PetscGlobalSum(PETSC_COMM_WORLD,&FluxOut, &FluxOutSum);

    //ratio = FluxInSum / FluxOutSum;
    ratio = FluxIn / simCtx->FluxOutSum;
    if (fabs(simCtx->FluxOutSum) < 1.e-6) ratio = 1.;
    //if (fabs(FluxInSum) <1.e-6) ratio = 0.;
    if (fabs(FluxIn) <1.e-6) ratio = 0.;
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Char Ratio %d %le %le %le %le %d %d\n", simCtx->step, ratio, FluxIn, simCtx->FluxOutSum, FarFluxInSum,zs, ze);

    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = ucat[k-1][j][i].y;
	  if (simCtx->step==0 || simCtx->step==1)
	    if (simCtx->inletprofile<0)
	      ubcs[k][j][i].z = -1.;
	    else if (user->bctype[4]==6)
	      ubcs[k][j][i].z = 0.;
	    else
	      ubcs[k][j][i].z = 1.;//ubcs[0][j][i].z;//-1.;//1.;
	  
	  else
	    ucont[k-1][j][i].z = ucont[k-1][j][i].z*ratio;
	  
	  ubcs[k][j][i].z = ucont[k-1][j][i].z / zet[k-1][j][i].z;
	}
      }
    }
  }

  
/* ==================================================================================             */
/*     OUTLETBC :4 */
/* ==================================================================================             */

  
  if (user->bctype[5]==OUTLET || user->bctype[5]==14 || user->bctype[5]==20) {
    lArea=0.;
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"+Zeta Outlet \n");
     LOG_ALLOW(GLOBAL,LOG_DEBUG,"FluxOutSum before FormBCS applied = %.6f \n",simCtx->FluxOutSum);
    if (ze == mz) {
      //    k = ze-3;
      k=ze-1;
      FluxOut = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {

	  if ((nvert[k-1][j][i])<0.1) {
	  FluxOut += (ucat[k-1][j][i].x * (zet[k-1][j][i].x) +
		      ucat[k-1][j][i].y * (zet[k-1][j][i].y) +
		      ucat[k-1][j][i].z * (zet[k-1][j][i].z));

	  lArea += sqrt( (zet[k-1][j][i].x) * (zet[k-1][j][i].x) +
			 (zet[k-1][j][i].y) * (zet[k-1][j][i].y) +
			 (zet[k-1][j][i].z) * (zet[k-1][j][i].z));
	 
	  }
	}
      }
    }
    else {
      FluxOut = 0.;
    }
    
    MPI_Allreduce(&FluxOut,&simCtx->FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea,&AreaSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
   
    user->simCtx->AreaOutSum = AreaSum;

    LOG_ALLOW(GLOBAL,LOG_DEBUG,"AreaOutSum = %.6f | FluxOutSum = %.6f \n",AreaSum,simCtx->FluxOutSum);
    
    if (simCtx->block_number>1 && user->bctype[5]==14) {
      simCtx->FluxOutSum += user->FluxIntfcSum;
      //      AreaSum    += user->AreaIntfcSum;
    }

    FluxIn = simCtx->FluxInSum + FarFluxInSum + user->FluxIntpSum;
    if (user->bctype[5]==20)
      ratio = (FluxIn / simCtx->FluxOutSum);
    else
      ratio = (FluxIn - simCtx->FluxOutSum) / AreaSum;
   
     LOG_ALLOW(GLOBAL,LOG_DEBUG,"Ratio for momentum correction = %.6f \n",ratio);
    
  /*   user->FluxOutSum += ratio*user->simCtx->AreaOutSum; */
    simCtx->FluxOutSum =0.0;
    FluxOut=0.0;
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if ((nvert[k-1][j][i])<0.1) {
	 
	    ubcs[k][j][i].x = ucat[k-1][j][i].x;//+ratio;
	    ubcs[k][j][i].y = ucat[k-1][j][i].y;
	    ubcs[k][j][i].z = ucat[k-1][j][i].z;// + ratio;//*n_z;

	    //  ucont[k-1][j][i].z = ubcs[k][j][i].z * zet[k-1][j][i].z;
	    if (user->bctype[5]==20)
	      ucont[k-1][j][i].z = (ubcs[k][j][i].x * (zet[k-1][j][i].x ) +
				    ubcs[k][j][i].y * (zet[k-1][j][i].y ) +
				    ubcs[k][j][i].z * (zet[k-1][j][i].z ))*ratio;
	    
	    else{
	      ucont[k-1][j][i].z = (ubcs[k][j][i].x * (zet[k-1][j][i].x ) +
				    ubcs[k][j][i].y * (zet[k-1][j][i].y ) +
				    ubcs[k][j][i].z * (zet[k-1][j][i].z ))
		+ ratio * sqrt( (zet[k-1][j][i].x) * (zet[k-1][j][i].x) +
				(zet[k-1][j][i].y) * (zet[k-1][j][i].y) +
				(zet[k-1][j][i].z) * (zet[k-1][j][i].z)); 

	      FluxOut += ucont[k-1][j][i].z;
	    }
	  }//if
	}
      }
    }
    
    MPI_Allreduce(&FluxOut,&simCtx->FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Timestep = %d | FluxInSum = %.6f | FlucOutSum = %.6f | FluxIntfcSum = %.6f | FluxIntpSum = %.6f \n", simCtx->step, simCtx->FluxInSum, simCtx->FluxOutSum, user->FluxIntfcSum,user->FluxIntpSum);

  } else if (user->bctype[5]==2) {
  /* Designed for driven cavity problem (top(k=kmax) wall moving)
   u_x = 1 at k==kmax */
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;// - ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = 1;//sin(2*3.14*simCtx->step*simCtx->dt);//1.;//- ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = 0.;//- ucat[k-1][j][i].z;
	}
      }
    }
  }
  

  /*   OUTLET at k==0 */
  if (user->bctype[4]==OUTLET) {
    lArea=0.;
    if (zs == 0) {
      k = zs;
      //      k= zs + 1;
      FluxOut = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {


	  FluxOut += ucat[k+1][j][i].z * zet[k][j][i].z ;

	  lArea += zet[k][j][i].z;



	}
      }
    }
    else {
      FluxOut = 0.;
    }
    
    FluxIn = simCtx->FluxInSum + FarFluxInSum;

    MPI_Allreduce(&FluxOut,&simCtx->FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea,&AreaSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  

    ratio = (simCtx->FluxInSum - simCtx->FluxOutSum) / AreaSum;
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Ratio b %d  %le %le %le %le %d %d\n", simCtx->step,ratio, simCtx->FluxInSum, simCtx->FluxOutSum, AreaSum,zs, ze);
    
    if (zs==0) {
      k = 0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	 
	  ubcs[k][j][i].x = ucat[k+1][j][i].x;
	  ubcs[k][j][i].y = ucat[k+1][j][i].y;
	  ubcs[k][j][i].z = ucat[k+1][j][i].z;
	  ucont[k][j][i].z = (ubcs[k][j][i].z+ratio) * zet[k][j][i].z;
	}
      }
    }
  }


  
/* ==================================================================================             */
/*     Ogrid :77 */
/* ==================================================================================             */
  /* 
  if (user->bctype[3]=77 && Ogrid)
    {Cmpnts ***coor;
      lArea=0.;
      FluxOut=0.0;
      //    k = ze-3;
      
      Vec Coor; DMGetCoordinatesLocal(da, &Coor); 
      DMDAVecGetArray(fda,Coor,&coor);       
      if (ye==my) {
	j=my-2;
	for (k=lzs; k<lze; k++){
	  for (i=lxs; i<lxe; i++) {
	    
	      FluxOut += ucont[k][j][i].y;
	      lArea += sqrt( (eta[k-1][j][i].x) * (eta[k-1][j][i].x) +
			     (eta[k-1][j][i].y) * (eta[k-1][j][i].y) +
			     (eta[k-1][j][i].z) * (eta[k-1][j][i].z));
	 
	  }
	}
      }
      
      else {
	FluxOut = 0.;
      }
      
      MPI_Allreduce(&FluxOut,&FluxOutSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      MPI_Allreduce(&lArea,&AreaSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
      
      user->FluxOutSum = FluxOutSum;
      simCtx->AreaOutSum = AreaSum;
     export PL="$HOME/CSBL/Codes/fdf_project/pic_les"
 ratio=2*(FluxIn-FluxOutSum)/AreaSum;
       ratio=0.0;   
      if (ye==my){
	j=my-2;
	for (k=lzs; k<lze; k++){      
	  for (i=lxs; i<lxe; i++) {	
	    if ((nvert[k-1][j][i])<0.1 && coor[k][j][i].z >= 800.) {
	      ucont[k][j][i].y *= (1+ratio);
	      
	    }
	  }
	}
      }
            
      if (ye==my){
	j=my-2;
	for (k=lzs; k<lze; k++){      
	  for (i=lxs; i<lxe; i++) {	
	    
	    ubcs[k][j+1][i].z=ucat[k][j][i].z;
	    ubcs[k][j+1][i].y=ucat[k][j][i].y;
	    ubcs[k][j+1][i].x=ucat[k][j][i].x;
	  }
	}
      }
      
      //   LOG_ALLOW(GLOBAL,LOG_DEBUG, "  ratio %le ",ratio);

    }

  */


/* ==================================================================================             */
/*     Channelz */
/* ==================================================================================             */
 // Amir channel flow correction
  if (user->bctype[4]==7 && simCtx->channelz==1) {
 Vec Coor; DMGetCoordinatesLocal(da, &Coor); 
 DMDAVecGetArray(fda,Coor,&coor); 
    Cmpnts  ***uch;
    DMDAVecGetArray(fda, user->Bcs.Uch, &uch);

    lArea=0.0;
   // if (zs==0) {
   // k=0;
    FluxIn=0.0;
    
    double Fluxbcs=0.0, Fluxbcssum,ratiobcs;

   if (zs==0) {
    k=0;
     Fluxbcs=0.0;      
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k][j][i]<0.1){
	    Fluxbcs += ucont[k][j][i].z;

	    lArea +=  sqrt((zet[k][j][i].x) * (zet[k][j][i].x) +
			  (zet[k][j][i].y) * (zet[k][j][i].y) +
			  (zet[k][j][i].z) * (zet[k][j][i].z));
	  }
	}
      }
    }
    

    //int kk=(simCtx->step % (mz-2))+2;
 
    for (k=zs;k<lze;k++){      
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k][j][i]<0.1){
	    FluxIn += ucont[k][j][i].z /((mz)-1);
	  }
	}
      }
    }
   // else {
  //   FluxIn=0.0;
//  }
 
    MPI_Allreduce(&FluxIn,&simCtx->Fluxsum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lArea,&AreaSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&Fluxbcs,&Fluxbcssum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
 
    if (simCtx->step==simCtx->StartStep && simCtx->StartStep > 0 && simCtx->ccc==0) {
      simCtx->ccc=1;
      simCtx->FluxInSum=Fluxbcssum;
//	simCtx->FluxInSum=6.3908; 
      LOG_ALLOW(LOCAL,LOG_DEBUG, "  FluxInSum %le .\n ",simCtx->FluxInSum);
    }
  
	simCtx->FluxInSum=6.35066; 
        ratio=(simCtx->FluxInSum-simCtx->Fluxsum)/AreaSum;
	simCtx->ratio=ratio;
        ratiobcs=(simCtx->FluxInSum-Fluxbcssum)/AreaSum;

    if (zs==0) {
	k=0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k+1][j][i]<0.1){
	    if (simCtx->fish) { 
			 ucont[k][j][i].z+=ratiobcs * /* (1.-exp(-500. * (1.-fabs(coor[k][j][i].y))))  */ sqrt( (zet[k+1][j][i].x) * (zet[k+1][j][i].x) + 
 					   (zet[k+1][j][i].y) * (zet[k+1][j][i].y) + 
 					   (zet[k+1][j][i].z) * (zet[k+1][j][i].z)); 
			}

	    uch[k][j][i].z=ratiobcs * /* (1.-exp(-500. * (1.-fabs(coor[k][j][i].y)))) */   sqrt( (zet[k+1][j][i].x) * (zet[k+1][j][i].x) +
					   (zet[k+1][j][i].y) * (zet[k+1][j][i].y) +
					   (zet[k+1][j][i].z) * (zet[k+1][j][i].z)); 
	  }
	}
      }
    }
  
    if (ze==mz) {
      k=mz-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k-1][j][i]<0.1){
 	    if (simCtx->fish){
		   ucont[k-1][j][i].z+=ratiobcs * /*(1.-exp(-500. * (1.-fabs(coor[k][j][i].y))))  */    sqrt((zet[k-1][j][i].x) * (zet[k-1][j][i].x) + 
 					   (zet[k-1][j][i].y) * (zet[k-1][j][i].y) + 
 					   (zet[k-1][j][i].z) * (zet[k-1][j][i].z));  	
	  }
	    uch[k][j][i].z=ratiobcs *   /*(1.-exp(-500. * (1.-fabs(coor[k][j][i].y)))) */     sqrt( (zet[k+1][j][i].x) * (zet[k+1][j][i].x) +
					   (zet[k+1][j][i].y) * (zet[k+1][j][i].y) +
					   (zet[k+1][j][i].z) * (zet[k+1][j][i].z)); 
	  }
	}
      }
    }
    DMDAVecRestoreArray(fda, user->Bcs.Uch, &uch);
    LOG_ALLOW(GLOBAL,LOG_DEBUG, "Ratio  %le %.15le %.15le  %.15le \n", ratio, ratiobcs, simCtx->FluxInSum,AreaSum);
    DMDAVecRestoreArray(fda,Coor,&coor); 

  ////.................////
   
   
  

    //just for check
/*         if (zs==0) { */
/*       k=10; */
/*       FluxIn=0.0; */
/*       for (j=lys; j<lye; j++) { */
/* 	for (i=lxs; i<lxe; i++) { */
/* 	  if (nvert[k+1][j][i]<0.1){ */
/* 	    FluxIn += ucont[k][j][i].z; */
/* 	    lArea += sqrt((zet[k+1][j][i].x) * (zet[k+1][j][i].x) + */
/* 			  (zet[k+1][j][i].y) * (zet[k+1][j][i].y) + */
/* 			  (zet[k+1][j][i].z) * (zet[k+1][j][i].z)); */
/* 	  } */
/* 	} */
/*       } */
/* 	} */
/*     else { */
/*       FluxIn=0.0; */
/*     } */
/*     // LOG_ALLOW(PETSC_COMM_SELF, "  Fluxsum %le ",FluxIn); */

/*     MPI_Allreduce(&FluxIn,&Fluxsum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */
/*     MPI_Allreduce(&lArea,&AreaSum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); */



  }//channel



  /*  ==================================================================================== */
  /*     Cylinder O-grid */
  /*  ==================================================================================== */
  if (user->bctype[3]==12) {
  /* Designed to test O-grid for flow over a cylinder at jmax velocity is 1 (horizontal) 
   u_x = 1 at k==kmax */
    if (ye==my) {
      j = ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = 1.;
	}
      }
    }
  }
  /*  ==================================================================================== */
  /*     Annulus */
  /*  ==================================================================================== */
  /* designed to test periodic boundary condition for O-grid j=0 rotates */
  DMDAVecGetArray(fda, user->Cent, &cent);
  if (user->bctype[2]==11) {
    if (ys==0){
      j=0;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	
	/*   ubcs[k][j][i].x=0.0; */
	 
/* 	  ubcs[k][j][i].y = -cent[k][j+1][i].z/sqrt(cent[k][j+1][i].z*cent[k][j+1][i].z+cent[k][j+1][i].y*cent[k][j+1][i].y); */
	 
/* 	  ubcs[k][j][i].z =cent[k][j+1][i].y/sqrt(cent[k][j+1][i].z*cent[k][j+1][i].z+cent[k][j+1][i].y*cent[k][j+1][i].y); */
	  ubcs[k][j][i].x = cent[k][j+1][i].y/sqrt(cent[k][j+1][i].x*cent[k][j+1][i].x+cent[k][j+1][i].y*cent[k][j+1][i].y);;
	  ubcs[k][j][i].y =-cent[k][j+1][i].x/sqrt(cent[k][j+1][i].x*cent[k][j+1][i].x+cent[k][j+1][i].y*cent[k][j+1][i].y);
	  ubcs[k][j][i].z =0.0;
	  //  if(k==1)  LOG_ALLOW(PETSC_COMM_SELF, "@ i= %d j=%d k=%d ubcs.y is %le\n",i,j,k,ubcs[k][j][i].y);
	}
      }
    }
  }
  /*  ==================================================================================== */
  /*     Rheology */
  /*  ==================================================================================== */
 
  if(simCtx->rheology && (user->bctype[2]==13 || user->bctype[3]==13 || user->bctype[4]==13 || user->bctype[5]==13)){
      LOG_ALLOW(GLOBAL,LOG_DEBUG, "moving plate velocity for rheology setup is %le \n",simCtx->U_bc);
  }
  
  if (user->bctype[2]==13){
    if (ys==0){
      j=0;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	   ubcs[k][j][i].x = 0.;
	  // ubcs[k][j][i].x = -simCtx->U_bc;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = -simCtx->U_bc;
	  //ubcs[k][j][i].z =0.0;
	}
      }
    }
  }
  if (user->bctype[3]==13){
    if (ye==my){
      j=ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;
	  // ubcs[k][j][i].x = simCtx->U_bc;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = simCtx->U_bc;
	   //ubcs[k][j][i].z =0.0;
	}
      }
    }
  }

 if (user->bctype[4]==13){
    if (zs==0){
      k=0;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x =-simCtx->U_bc;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = 0.;
	}
      }
    }
  }
  if (user->bctype[5]==13){
    if (ze==mz){
      k=ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = simCtx->U_bc;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = 0.;
	}
      }
    }
  }
  DMDAVecRestoreArray(fda, user->Cent, &cent);
  DMDAVecRestoreArray(fda, user->lUcat, &ucat);
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  DMDAVecRestoreArray(da, user->Nvert, &nvert);
  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

 
 
  Contra2Cart(user); // it also does global to local for Ucat
 
/* ==================================================================================             */
/*   WALL FUNCTION */
/* ==================================================================================             */

  if (simCtx->wallfunction && user->bctype[2]==-1) {
  PetscReal ***ustar, ***aj;
  //Mohsen Dec 2015
  Vec Aj  =  user->lAj;
  DMDAVecGetArray(fda, user->Ucat, &ucat);
  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(da, Aj,  &aj);
//  DMDAVecGetArray(da, user->Nvert, &nvert);
	DMDAVecGetArray(da, user->lUstar, &ustar);
 
  // wall function for boundary
  for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
      for (i=lxs; i<lxe; i++) {

	if( nvert[k][j][i]<1.1 &&  user->bctype[2]==-1 && j==1 )
	 {
	  double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
	  double sb, sc;
	  double ni[3], nj[3], nk[3];
	  Cmpnts Uc, Ua;
	  
	  Ua.x = Ua.y = Ua.z = 0;
	  sb = 0.5/aj[k][j][i]/area;
	  
	  
	    sc = 2*sb + 0.5/aj[k][j+1][i]/area;
	    Uc = ucat[k][j+1][i];
	  
	  
	  //Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
	  //if(j==my-2) nj[0]*=-1, nj[1]*=-1, nj[2]*=-1;
	
         double AA=sqrt(eta[k][j][i].z*eta[k][j][i].z +
               eta[k][j][i].y*eta[k][j][i].y +
               eta[k][j][i].x*eta[k][j][i].x);
 	nj[0]=eta[k][j][i].x/AA;
        nj[1]=eta[k][j][i].y/AA;
        nj[2]=eta[k][j][i].z/AA;     
	  noslip (user, sc, sb, Ua, Uc, &ucat[k][j][i], nj[0], nj[1], nj[2]);
	wall_function_loglaw(user, 1.e-16, sc, sb, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], nj[0], nj[1], nj[2]);

	 // nvert[k][j][i]=1.;	/* set nvert to 1 to exclude from rhs */
	// if (k==1) 
	  // LOG_ALLOW(GLOBAL,LOG_DEBUG, " %d   %le   %le  %le   %le   %le   %le   %le   %le   %le\n",i, sb,aj[k][j][i],AA, nj[0], nj[1], nj[2],ucat[k][j][i].x,ucat[k][j][i].y,ucat[k][j][i].z);

	}
      }
  if (ys==0) {
    j= ys;

    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
        ubcs[k][j][i].x = 0.0;
        ubcs[k][j][i].y = 0.0;
        ubcs[k][j][i].z = 0.0;
        ucont[k][j][i].y = 0.;
      }
    }
    }

  DMDAVecRestoreArray(da, Aj,  &aj);
  DMDAVecRestoreArray(da, user->lUstar, &ustar);
  DMDAVecRestoreArray(fda, user->Ucat, &ucat);
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
 // DMDAVecRestoreArray(da, user->Nvert, &nvert);
 // DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
 // DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
//   DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); 
//   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); 

 
  }

/* ==================================================================================             */

  DMDAVecGetArray(fda, user->Ucat, &ucat);
 
/* ==================================================================================             */
 
  // boundary conditions on ghost nodes
  if (xs==0 && user->bctype[0]!=7) {
    i = xs;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i+1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i+1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i+1].z;
      }
    }
  }

  if (xe==mx && user->bctype[0]!=7) {
    i = xe-1;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i-1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i-1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i-1].z;
      }
    }
  }


  if (ys==0 && user->bctype[2]!=7) {
    j = ys;
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j+1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j+1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j+1][i].z;
      }
    }
  }

  if (ye==my && user->bctype[2]!=7) {
    j = ye-1;
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j-1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j-1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j-1][i].z;
      }
    }
  }
 
  if (zs==0 && user->bctype[4]!=7) {
    k = zs;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k+1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k+1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k+1][j][i].z;
      }
    }
  }

  if (ze==mz && user->bctype[4]!=7) {
    k = ze-1;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k-1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k-1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k-1][j][i].z;
      }
    }
  }

  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
 /* ==================================================================================             */
  /*   Periodic BC *///Mohsen
/* ==================================================================================             */
  if (user->bctype[0]==7 || user->bctype[2]==7 || user->bctype[4]==7){
    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    /* /\*   DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert); *\/ */
    /* /\*   DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert); *\/ */
  //Mohsen Dec 2015
    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lNvert, &lnvert);
    DMDAVecGetArray(fda, user->lUcat, &lucat);
    DMDAVecGetArray(da, user->P, &p);
    DMDAVecGetArray(da, user->Nvert, &nvert);
    DMDAVecGetArray(fda, user->Ucat, &ucat);
   
    if (user->bctype[0]==7 || user->bctype[1]==7){
      if (xs==0){
	i=xs;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    if(j>0 && k>0 && j<user->JM && k<user->KM){
	      ucat[k][j][i]=lucat[k][j][i-2];
	      p[k][j][i]=lp[k][j][i-2];
	      nvert[k][j][i]=lnvert[k][j][i-2];
	    }
	  }
	}
      }
    }
    if (user->bctype[2]==7 || user->bctype[3]==7){
      if (ys==0){
	j=ys;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    if(i>0 && k>0 && i<user->IM && k<user->KM){
	      ucat[k][j][i]=lucat[k][j-2][i];
	      p[k][j][i]=lp[k][j-2][i];
	      nvert[k][j][i]=lnvert[k][j-2][i];
	    }
	  }
	}
      }
    }
    if (user->bctype[4]==7 || user->bctype[5]==7){
      if (zs==0){
	k=zs;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    if(i>0 && j>0 && i<user->IM && j<user->JM){
	      ucat[k][j][i]=lucat[k-2][j][i];
	      nvert[k][j][i]=lnvert[k-2][j][i];
	      //amir
		p[k][j][i]=lp[k-2][j][i];
	    }
	  }
	}
      }
    }
    if (user->bctype[0]==7 || user->bctype[1]==7){
      if (xe==mx){
	i=mx-1;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    if(j>0 && k>0 && j<user->JM && k<user->KM){
	      ucat[k][j][i]=lucat[k][j][i+2];
	      p[k][j][i]=lp[k][j][i+2];
	      nvert[k][j][i]=lnvert[k][j][i+2];
	    }
	  }
	}
      }
    }
    if (user->bctype[2]==7 || user->bctype[3]==7){
      if (ye==my){
	j=my-1;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    if(i>0 && k>0 && i<user->IM && k<user->KM){
	      ucat[k][j][i]=lucat[k][j+2][i];
	      p[k][j][i]=lp[k][j+2][i];
	      nvert[k][j][i]=lnvert[k][j+2][i];
	    }
	  }
	}
      }
    }
  
    if (user->bctype[4]==7 || user->bctype[5]==7){
      if (ze==mz){
	k=mz-1;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    if(i>0 && j>0 && i<user->IM && j<user->JM){
	      ucat[k][j][i]=lucat[k+2][j][i];
      	      nvert[k][j][i]=lnvert[k+2][j][i]; 
	      //amir
		p[k][j][i]=lp[k+2][j][i];
	    }
	  }
	}
      }
    }

       


 
    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lNvert, &lnvert);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);
    DMDAVecRestoreArray(da, user->P, &p);
    DMDAVecRestoreArray(da, user->Nvert, &nvert);

  /*  /\*  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); *\/ */
  /* /\*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); *\/ */

  /* /\*   DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P); *\/ */
  /* /\*   DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P); *\/ */

  /* /\*   DMLocalToGlobalBegin(da, user->lNvert, INSERT_VALUES, user->Nvert); *\/ */
  /* /\*   DMLocalToGlobalEnd(da, user->lNvert, INSERT_VALUES, user->Nvert); *\/ */

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
    DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);
}
 // 0 velocity on the corner point

  DMDAVecGetArray(fda, user->Ucat, &ucat);
  DMDAVecGetArray(da, user->P, &p);
  
  if (zs==0) {
    k=0;
    if (xs==0) {
      i=0;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j][i+1].z);
	p[k][j][i]= 0.5*(p[k+1][j][i]+p[k][j][i+1]);
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j][i-1].z);
	p[k][j][i] = 0.5*(p[k+1][j][i]+p[k][j][i-1]);
      }
    }

    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j+1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j+1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j+1][i].z);
	p[k][j][i] = 0.5*(p[k+1][j][i]+p[k][j+1][i]);
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j-1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j-1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j-1][i].z);
	p[k][j][i] = 0.5*(p[k+1][j][i]+p[k][j-1][i]);
      }
    }
  }
 
  if (ze==mz) {
    k=mz-1;
    if (xs==0) {
      i=0;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j][i+1].z);
	p[k][j][i] = 0.5*(p[k-1][j][i]+p[k][j][i+1]);
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j][i-1].z);
	p[k][j][i] = 0.5*(p[k-1][j][i]+p[k][j][i-1]);
      }
    }

    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j+1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j+1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j+1][i].z);
	p[k][j][i] = 0.5*(p[k-1][j][i]+p[k][j+1][i]);
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j-1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j-1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j-1][i].z);
	p[k][j][i] = 0.5*(p[k-1][j][i]+p[k][j-1][i]);
      }
    }
  }
 
  if (ys==0) {
    j=0;
    if (xs==0) {
      i=0;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j+1][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j+1][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j+1][i].z+ucat[k][j][i+1].z);
	p[k][j][i]= 0.5*(p[k][j+1][i]+p[k][j][i+1]);
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j+1][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j+1][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j+1][i].z+ucat[k][j][i-1].z);
	p[k][j][i] = 0.5*(p[k][j+1][i]+p[k][j][i-1]);
      }
    }
  }
 
  if (ye==my) {
    j=my-1;
    if (xs==0) {
      i=0;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j-1][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j-1][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j-1][i].z+ucat[k][j][i+1].z);
	p[k][j][i] = 0.5*(p[k][j-1][i]+p[k][j][i+1]);
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j-1][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j-1][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j-1][i].z+ucat[k][j][i-1].z);
	p[k][j][i] = 0.5*(p[k][j-1][i]+p[k][j][i-1]);
      }
    }
  }

  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  DMDAVecRestoreArray(da, user->P, &p);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);

  //Mohsen Nov 2012
  //Velocity and Presurre at corners for Periodic BC's

  if (user->bctype[0]==7 || user->bctype[2]==7 || user->bctype[4]==7){
  //i-direction

    DMDAVecGetArray(fda, user->lUcat,  &lucat);
    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lNvert, &lnvert);
    DMDAVecGetArray(fda, user->Ucat,  &ucat);
    DMDAVecGetArray(da, user->P, &p);
    DMDAVecGetArray(da, user->Nvert, &nvert);
    
    if (user->bctype[0]==7){
      if (xs==0){
	i=xs;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    ucat[k][j][i]=lucat[k][j][i-2];
	    p[k][j][i]=lp[k][j][i-2];
	    nvert[k][j][i]=lnvert[k][j][i-2];
	  }
	}
      }
    }
    if (user->bctype[1]==7){
      if (xe==mx){
	i=xe-1;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    ucat[k][j][i].x=lucat[k][j][i+2].x;
	    p[k][j][i]=lp[k][j][i+2];
	    nvert[k][j][i]=lnvert[k][j][i+2];
	  }
	}
      }
    }
    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lNvert, &lnvert);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);
    DMDAVecRestoreArray(da, user->P, &p);
    DMDAVecRestoreArray(da, user->Nvert, &nvert);

 /*  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */

/*   DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P); */
/*   DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P); */

/*   DMLocalToGlobalBegin(da, user->lNvert, INSERT_VALUES, user->Nvert); */
/*   DMLocalToGlobalEnd(da, user->lNvert, INSERT_VALUES, user->Nvert); */
  
    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
    DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);
    
    //j-direction
    DMDAVecGetArray(fda, user->lUcat,  &lucat);
    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lNvert, &lnvert);
    DMDAVecGetArray(fda, user->Ucat,  &ucat);
    DMDAVecGetArray(da, user->P, &p);
    DMDAVecGetArray(da, user->Nvert, &nvert);
    
    if (user->bctype[2]==7){
      if (ys==0){
	j=ys;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    ucat[k][j][i]=lucat[k][j-2][i];
	    p[k][j][i]=lp[k][j-2][i];
	    nvert[k][j][i]=lnvert[k][j-2][i];
	  }
	}
      }
    }
    
    if (user->bctype[3]==7){
      if (ye==my){
	j=my-1;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    ucat[k][j][i]=lucat[k][j+2][i];
	    p[k][j][i]=lp[k][j+2][i];
	    nvert[k][j][i]=lnvert[k][j+2][i];
	  }
	}
      }
    }
    
    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lNvert, &lnvert);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);
    DMDAVecRestoreArray(da, user->P, &p);
    DMDAVecRestoreArray(da, user->Nvert, &nvert);
    
/*   DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */

/*   DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P); */
/*   DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P); */

/*   DMLocalToGlobalBegin(da, user->lNvert, INSERT_VALUES, user->Nvert); */
/*   DMLocalToGlobalEnd(da, user->lNvert, INSERT_VALUES, user->Nvert); */

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
    DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);
    
    //k-direction
    DMDAVecGetArray(fda, user->lUcat,  &lucat);
    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lNvert, &lnvert);
    DMDAVecGetArray(fda, user->Ucat,  &ucat);
    DMDAVecGetArray(da, user->P, &p);
    DMDAVecGetArray(da, user->Nvert, &nvert);
    
    if (user->bctype[4]==7){
      if (zs==0){
	k=zs;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    ucat[k][j][i]=lucat[k-2][j][i];
	    nvert[k][j][i]=lnvert[k-2][j][i]; 
	    //amir   
	      p[k][j][i]=lp[k-2][j][i];
	  }
	}
      }
    }
    if (user->bctype[5]==7){
      if (ze==mz){
	k=mz-1;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    ucat[k][j][i]=lucat[k+2][j][i];
	    nvert[k][j][i]=lnvert[k+2][j][i];
	      p[k][j][i]=lp[k+2][j][i];
	  }
	}
      }
    }
    
    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lNvert, &lnvert);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);
    DMDAVecRestoreArray(da, user->P, &p);
    DMDAVecRestoreArray(da, user->Nvert, &nvert);
    
/*   DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */

/*   DMLocalToGlobalBegin(da, user->lP, INSERT_VALUES, user->P); */
/*   DMLocalToGlobalEnd(da, user->lP, INSERT_VALUES, user->P); */

/*   DMLocalToGlobalBegin(da, user->lNvert, INSERT_VALUES, user->Nvert); */
/*   DMLocalToGlobalEnd(da, user->lNvert, INSERT_VALUES, user->Nvert); */

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    
    DMGlobalToLocalBegin(da, user->Nvert, INSERT_VALUES, user->lNvert);
    DMGlobalToLocalEnd(da, user->Nvert, INSERT_VALUES, user->lNvert);
  }
  /* ==================================================================================             */
/*   Analytical Vortex BC */
/* ==================================================================================             */
 
  DMDAVecGetArray(fda, user->Ucat, &ucat);
  DMDAVecGetArray(fda, user->Ucont, &ucont);
  DMDAVecGetArray(fda, user->Cent, &cent); 
  DMDAVecGetArray(fda, user->Centx, &centx);
  DMDAVecGetArray(fda, user->Centy, &centy);
  DMDAVecGetArray(fda, user->Centz, &centz);

  if (user->bctype[0]==9) {
    if (xs==0) {
      i= xs;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ucat[k][j][i].x=-cos(cent[k][j][i+1].x)*sin(cent[k][j][i+1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].y=-sin(cent[k][j][i+1].x)*cos(cent[k][j][i+1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].z =0.0;
	 
	  ucont[k][j][i].x =-(cos(centx[k][j][i].x)*sin(centx[k][j][i].y)*csi[k][j][i].x)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	}
      }
      if (ys==0) {
	j=ys;
	for (k=lzs; k<lze; k++) {
	  ucat[k][j][i].x=cos(cent[k][j+1][i+1].x)*sin(cent[k][j+1][i+1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].y=-sin(cent[k][j+1][i+1].x)*cos(cent[k][j+1][i+1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].z =0.0;
	}
      }
      if (ye==my) {
	j=ye-1;
	for (k=lzs; k<lze; k++) {
	  ucat[k][j][i].x=cos(cent[k][j-1][i+1].x)*sin(cent[k][j-1][i+1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].y=-sin(cent[k][j-1][i+1].x)*cos(cent[k][j-1][i+1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].z =0.0;
	}
      }
    }
  }
  if (user->bctype[1]==9) {
    if (xe==mx) {
      i= xe-1;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  ucat[k][j][i].x=-cos(cent[k][j][i-1].x)*sin(cent[k][j][i-1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].y=-sin(cent[k][j][i-1].x)*cos(cent[k][j][i-1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].z =0.0;
	
	  ucont[k][j][i-1].x =(-cos(centx[k][j][i-1].x)*sin(centx[k][j][i-1].y)*csi[k][j][i-1].x)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	}
      }
      if (ys==0) {
	j=ys;
	for (k=lzs; k<lze; k++) {
	  ucat[k][j][i].x=cos(cent[k][j+1][i-1].x)*sin(cent[k][j+1][i-1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].y=-sin(cent[k][j+1][i-1].x)*cos(cent[k][j+1][i-1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].z =0.0;
	}
      }
      if (ye==my) {
	j=ye-1;
	for (k=lzs; k<lze; k++) {
	  ucat[k][j][i].x=cos(cent[k][j-1][i-1].x)*sin(cent[k][j-1][i-1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].y=-sin(cent[k][j-1][i-1].x)*cos(cent[k][j-1][i-1].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].z =0.0;
	}
      }
    }
  }

  if (user->bctype[2]==9) {
    if (ys==0) {
      j= ys;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ucat[k][j][i].x=cos(cent[k][j+1][i].x)*sin(cent[k][j+1][i].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].y=sin(cent[k][j+1][i].x)*cos(cent[k][j+1][i].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].z =0.0;

	  ucont[k][j][i].y=(sin(centy[k][j][i].x)*cos(centy[k][j][i].y)*eta[k][j][i].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	}
      }
    }
  }
 
 
  if (user->bctype[3]==9) {
    if (ye==my) {
      j= ye-1;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ucat[k][j][i].x=cos(cent[k][j-1][i].x)*sin(cent[k][j-1][i].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].y=sin(cent[k][j-1][i].x)*cos(cent[k][j-1][i].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	  ucat[k][j][i].z =0.0;
	
	  ucont[k][j-1][i].y=(sin(centy[k][j-1][i].x)*cos(centy[k][j-1][i].y)*eta[k][j-1][i].y)*exp(-2.0*simCtx->dt*(simCtx->step+1)/simCtx->ren);
	}
      }
    }
  }
  if (user->bctype[4]==9) {
    if (zs==0) {
      k= zs;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  ucat[k][j][i].x=ucat[k+1][j][i].x;
	  ucat[k][j][i].y=ucat[k+1][j][i].y;
	  ucat[k][j][i].z=ucat[k+1][j][i].z;

	  ucont[k][j][i].z=0.0;
	}
      }
    }
  }
  if (user->bctype[5]==9) {
    if (ze==mz) {
      k= ze-1;
      for (j=ys; j<ye; j++) {
	for (i=xs; i<xe; i++) {
	  ucat[k][j][i].x=ucat[k-1][j][i].x;
	  ucat[k][j][i].y=ucat[k-1][j][i].y;
	  ucat[k][j][i].z=ucat[k-1][j][i].z;

	  ucont[k-1][j][i].z=0.0;
	}
      }
    }
  }

  DMDAVecRestoreArray(fda, user->Cent, &cent);
  DMDAVecRestoreArray(fda, user->Centx, &centx);
  DMDAVecRestoreArray(fda, user->Centy, &centy);
  DMDAVecRestoreArray(fda, user->Centz, &centz);

  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);
  DMDAVecRestoreArray(fda, user->Ucont,  &ucont);
 
  } // ttemp

 
  DMDAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecRestoreArray(fda, user->lCsi,  &csi);
  DMDAVecRestoreArray(fda, user->lEta,  &eta);
  DMDAVecRestoreArray(fda, user->lZet,  &zet);

  DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont); 

  PROFILE_FUNCTION_END;

  PetscFunctionReturn(0);
  }
