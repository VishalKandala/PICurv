#include "Boundaries.h" // The main header for our project
#include <string.h>    // For strcasecmp
#include <ctype.h>     // For isspace

#undef __FUNCT__
#define __FUNCT__ "CanRankServiceInletFace"
/**
 * @brief Internal helper implementation: `CanRankServiceInletFace()`.
 * @details Local to this translation unit.
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
        PROFILE_FUNCTION_END;
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

      LOG_ALLOW(LOCAL, LOG_INFO,"[Rank %d] Inlet Face %s Service Check Result: %s | Owned Cells (I,J,K): (%d,%d,%d) | Starts at Cell (%d,%d,%d)\n",
                rank_for_logging,
                BCFaceToString((BCFace)user->identifiedInletBCFace),
                (*can_service_inlet_out) ? "CAN SERVICE" : "CANNOT SERVICE",
                num_owned_cells_on_rank_i, num_owned_cells_on_rank_j, num_owned_cells_on_rank_k,
                owned_start_cell_i, owned_start_cell_j, owned_start_cell_k);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CanRankServiceFace"

/**
 * @brief Implementation of \ref CanRankServiceFace().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Boundaries.h`.
 * @see CanRankServiceFace()
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
 * @brief Internal helper implementation: `GetDeterministicFaceGridLocation()`.
 * @details Local to this translation unit.
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
    PetscReal global_logic_i = 0.0, global_logic_j = 0.0, global_logic_k = 0.0;
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
    const PetscReal layer_spacing_norm_i = (IM_cells_global > 0) ? 1.0 / (PetscReal)IM_cells_global : 0.0;
    const PetscReal layer_spacing_norm_j = (JM_cells_global > 0) ? 1.0 / (PetscReal)JM_cells_global : 0.0;
    const PetscReal layer_spacing_norm_k = (KM_cells_global > 0) ? 1.0 / (PetscReal)KM_cells_global : 0.0;

    // Grid-aware epsilon: scale with minimum cell size to keep particles away from rank boundaries
    const PetscReal min_layer_spacing = PetscMin(layer_spacing_norm_i, PetscMin(layer_spacing_norm_j, layer_spacing_norm_k));
    const PetscReal epsilon = 0.5 * min_layer_spacing;  // Keep particles 10% of cell width from boundaries

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
        "[Rank %d] Particle %lld placement %s.\n",
        rank_for_logging, (long long)particle_global_id,
        (*placement_successful_out ? "SUCCESSFUL" : "NOT ON THIS RANK"));

    if(*placement_successful_out){    
        LOG_ALLOW(LOCAL,LOG_TRACE,"Local cell origin node: (I,J,K) = (%d,%d,%d), intra-cell logicals: (xi,eta,zta)=(%.6f,%.6f,%.6f)\n",
                    *ci_metric_lnode_out, *cj_metric_lnode_out, *ck_metric_lnode_out,
                    *xi_metric_logic_out, *eta_metric_logic_out, *zta_metric_logic_out);
        }
        
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GetRandomFCellAndLogicOnInletFace"

/**
 * @brief Internal helper implementation: `GetRandomCellAndLogicalCoordsOnInletFace()`.
 * @details Local to this translation unit.
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

    // Get number of cells this rank owns in each dimension (tangential to the face mainly)
    PetscInt owned_start_cell_i, num_owned_cells_on_rank_i;
    PetscInt owned_start_cell_j, num_owned_cells_on_rank_j;
    PetscInt owned_start_cell_k, num_owned_cells_on_rank_k;

    ierr = GetOwnedCellRange(info, 0, &owned_start_cell_i, &num_owned_cells_on_rank_i); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 1, &owned_start_cell_j, &num_owned_cells_on_rank_j); CHKERRQ(ierr);
    ierr = GetOwnedCellRange(info, 2, &owned_start_cell_k, &num_owned_cells_on_rank_k); CHKERRQ(ierr);

    // Defaults for cell origin node (local index for the rank's DA patch, including ghosts)
    *ci_metric_lnode_out = xs_gnode_rank; *cj_metric_lnode_out = ys_gnode_rank; *ck_metric_lnode_out = zs_gnode_rank;
    // Defaults for logical coordinates
    *xi_metric_logic_out = 0.5; *eta_metric_logic_out = 0.5; *zta_metric_logic_out = 0.5;

    // Index of the last cell (0-indexed) in each global direction
    PetscInt last_global_cell_idx_i = (IM_nodes_global > 1) ? (IM_nodes_global - 2) : -1;
    PetscInt last_global_cell_idx_j = (JM_nodes_global > 1) ? (JM_nodes_global - 2) : -1;
    PetscInt last_global_cell_idx_k = (KM_nodes_global > 1) ? (KM_nodes_global - 2) : -1;
    
    LOG_ALLOW(LOCAL, LOG_INFO, "PARTICLE_INIT_DEBUG Rank %d: Inlet face %s.\n"
        "  Owned cells (i,j,k): (%d,%d,%d)\n"
        "  Global nodes (I,J,K): (%d,%d,%d)\n"
        "  info->xs,ys,zs (first owned node GLOBAL): (%d,%d,%d)\n"
        "  info->xm,ym,zm (num owned nodes GLOBAL): (%d,%d,%d)\n"
        "  xs_gnode_rank,ys_gnode_rank,zs_gnode_rank (DMDAGetCorners): (%d,%d,%d)\n"
        "  owned_start_cell (i,j,k) GLOBAL: (%d,%d,%d)\n"
        "  last_global_cell_idx (i,j,k): (%d,%d,%d)\n",
        rank_for_logging, BCFaceToString((BCFace)user->identifiedInletBCFace),
        num_owned_cells_on_rank_i,num_owned_cells_on_rank_j,num_owned_cells_on_rank_k,
        IM_nodes_global,JM_nodes_global,KM_nodes_global,
        info->xs, info->ys, info->zs,
        info->xm, info->ym, info->zm,
        xs_gnode_rank,ys_gnode_rank,zs_gnode_rank,
        owned_start_cell_i, owned_start_cell_j, owned_start_cell_k,
        last_global_cell_idx_i, last_global_cell_idx_j, last_global_cell_idx_k);


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
#define __FUNCT__ "ClassifyMomentumRow"
/**
 * @brief Implementation of \ref ClassifyMomentumRow().
 *
 * Pure index/boundary-type arithmetic: no field reads, no communication.
 * Precedence matters. A conditioned row is reported first because its explicit
 * Dirichlet value is more specific than the homogeneous fallback, and a periodic
 * duplicate is reported before the homogeneous case because the Newton path needs
 * its representative index to build `F = X_dup - X_rep`.
 */
MomentumRowType ClassifyMomentumRow(UserCtx *user, PetscInt i, PetscInt j, PetscInt k,
                                    PetscInt component, PetscInt *ri, PetscInt *rj, PetscInt *rk)
{
    const PetscInt mx = user->info.mx, my = user->info.my, mz = user->info.mz;
    const PetscInt coord[3] = {i, j, k};
    const PetscInt size[3] = {mx, my, mz};
    const BCFace neg_face[3] = {BC_FACE_NEG_X, BC_FACE_NEG_Y, BC_FACE_NEG_Z};
    PetscBool periodic[3], periodic_duplicate = PETSC_FALSE;
    PetscBool residual_zeroed = PETSC_FALSE, conditioned = PETSC_FALSE;

    *ri = i; *rj = j; *rk = k;
    for (PetscInt axis = 0; axis < 3; ++axis) {
        periodic[axis] = (PetscBool)(
            user->boundary_faces[neg_face[axis]].mathematical_type == PERIODIC);
        if (periodic[axis] && coord[axis] == 0) {
            periodic_duplicate = PETSC_TRUE;
            if (axis == 0) *ri = -2;
            else if (axis == 1) *rj = -2;
            else *rk = -2;
        }
        if (periodic[axis] && coord[axis] == size[axis] - 1) {
            periodic_duplicate = PETSC_TRUE;
            if (axis == 0) *ri = mx + 1;
            else if (axis == 1) *rj = my + 1;
            else *rk = mz + 1;
        }

        if (!periodic[axis] && coord[axis] == 0) residual_zeroed = PETSC_TRUE;
        if (coord[axis] == size[axis] - 1) residual_zeroed = PETSC_TRUE;
        if (!periodic[axis] && coord[axis] == size[axis] - 2 && component == axis)
            residual_zeroed = PETSC_TRUE;
    }

    if (!periodic[component] &&
        (coord[component] == 0 || coord[component] == size[component] - 2)) {
        PetscBool tangential_interior = PETSC_TRUE;
        for (PetscInt axis = 0; axis < 3; ++axis) {
            if (axis == component) continue;
            if (coord[axis] < 1 || coord[axis] > size[axis] - 2)
                tangential_interior = PETSC_FALSE;
        }
        conditioned = tangential_interior;
    }

    if (conditioned) return MOM_ROW_FIXED_CONDITIONED;
    if (periodic_duplicate) return MOM_ROW_PERIODIC_DUPLICATE;
    if (residual_zeroed) return MOM_ROW_FIXED_HOMOGENEOUS;
    return MOM_ROW_PHYSICAL;
}

#undef __FUNCT__
#define __FUNCT__ "MomentumRowIsSolidMasked"
/**
 * @brief Implementation of \ref MomentumRowIsSolidMasked().
 * @details Full API contract is documented with the header declaration in
 *          `include/Boundaries.h`.
 * @see MomentumRowIsSolidMasked()
 */
PetscBool MomentumRowIsSolidMasked(const PetscReal ***nvert, PetscInt i, PetscInt j,
                                   PetscInt k, PetscInt component)
{
    /* PICURV_SOLID_THRESHOLD mirrors the 0.1 fluid test ComputeRHS() applies. */
    const PetscReal threshold = 0.1;

    if (nvert == NULL) return PETSC_FALSE;
    /* A solid cell carries no momentum equation in any component. */
    if (nvert[k][j][i] > threshold) return PETSC_TRUE;
    /* A staggered row also loses its equation when the cell it points into is solid. */
    switch (component) {
    case 0: return (PetscBool)(nvert[k][j][i + 1] > threshold);
    case 1: return (PetscBool)(nvert[k][j + 1][i] > threshold);
    case 2: return (PetscBool)(nvert[k + 1][j][i] > threshold);
    default: return PETSC_FALSE;
    }
}

#undef __FUNCT__
#define __FUNCT__ "EnforceRHSBoundaryConditions"
/**
 * @brief Implementation of \ref EnforceRHSBoundaryConditions().
 *
 * The sweep is deliberately expressed over every owned location rather than over
 * the six boundary slabs: restating "which indices can be non-physical" here is
 * exactly the duplication that let the periodic duplicate column go unzeroed.
 * ClassifyMomentumRow() is a handful of integer comparisons and the walk is a
 * single pass with no stencil access, which is negligible next to the several
 * ghosted stencil passes ComputeRHS() has already made over the same range.
 * MomentumNewtonKrylov_ApplyConstraints() walks the same range the same way.
 */
PetscErrorCode EnforceRHSBoundaryConditions(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info = user->info;
    Cmpnts         ***rhs;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // Get a writable pointer to the local data of the global RHS vector.
    ierr = DMDAVecGetArray(user->fda, user->Rhs, &rhs); CHKERRQ(ierr);

    for (PetscInt k = info.zs; k < info.zs + info.zm; k++) {
        for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
            for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
                PetscScalar *row = &rhs[k][j][i].x;   /* .x/.y/.z are contiguous */
                for (PetscInt component = 0; component < 3; component++) {
                    PetscInt ri, rj, rk;
                    if (ClassifyMomentumRow(user, i, j, k, component, &ri, &rj, &rk) != MOM_ROW_PHYSICAL)
                        row[component] = 0.0;
                }
            }
        }
    }

    // --- Release the pointer to the local data ---
    ierr = DMDAVecRestoreArray(user->fda, user->Rhs, &rhs); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL, LOG_TRACE, "Rank %d, Block %d: Finished enforcing RHS boundary conditions.\n",
              user->simCtx->rank, user->_this);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "BoundaryCondition_Create"
/**
 * @brief Internal helper implementation: `BoundaryCondition_Create()`.
 * @details Local to this translation unit.
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
    bc->UpdateUbcs  = NULL;
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

        case BC_HANDLER_PERIODIC_GEOMETRIC:
            LOG_ALLOW(LOCAL,LOG_DEBUG,"Dispatching to Create_PeriodicGeometric().\n");
            ierr = Create_PeriodicGeometric(bc);
            break;
        
        case BC_HANDLER_PERIODIC_DRIVEN_CONSTANT_FLUX:
            LOG_ALLOW(LOCAL,LOG_DEBUG,"Dispatching to Create_PeriodicDrivenConstant().\n");
            ierr = Create_PeriodicDrivenConstant(bc);
            break;
        
        case BC_HANDLER_PERIODIC_DRIVEN_INITIAL_FLUX:
            LOG_ALLOW(LOCAL,LOG_DEBUG,"Dispatching to Create_PeriodicDrivenInitial().\n");
            ierr = Create_PeriodicDrivenInitial(bc);
            break;
                
        case BC_HANDLER_INLET_PARABOLIC:
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Dispatching to Create_InletParabolicProfile().\n");
            ierr = Create_InletParabolicProfile(bc); CHKERRQ(ierr);
            break;

        case BC_HANDLER_INLET_PROFILE_FROM_FILE:
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Dispatching to Create_InletProfileFromFile().\n");
            ierr = Create_InletProfileFromFile(bc); CHKERRQ(ierr);
            break;
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

#undef __FUNCT__
#define __FUNCT__ "BoundarySystem_Validate"
/**
 * @brief Internal helper implementation: `BoundarySystem_Validate()`.
 * @details Local to this translation unit.
 */
PetscErrorCode BoundarySystem_Validate(UserCtx *user)
{
    PetscErrorCode ierr;
    const BCFace neg_faces[3] = {BC_FACE_NEG_X, BC_FACE_NEG_Y, BC_FACE_NEG_Z};
    const BCFace pos_faces[3] = {BC_FACE_POS_X, BC_FACE_POS_Y, BC_FACE_POS_Z};
    const char axis_names[3] = {'X', 'Y', 'Z'};
    DMBoundaryType bx, by, bz;
    PetscBool dm_periodic[3];
    PetscFunctionBeginUser;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Validating parsed boundary condition configuration...\n");
    ierr = DMDAGetInfo(user->da, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                       &bx, &by, &bz, NULL); CHKERRQ(ierr);
    dm_periodic[0] = (PetscBool)(bx == DM_BOUNDARY_PERIODIC);
    dm_periodic[1] = (PetscBool)(by == DM_BOUNDARY_PERIODIC);
    dm_periodic[2] = (PetscBool)(bz == DM_BOUNDARY_PERIODIC);

    // --- Rule Set 1: Geometric periodic faces must be paired and match the DM topology. ---
    for (PetscInt axis = 0; axis < 3; axis++) {
        const PetscBool neg_periodic =
            user->boundary_faces[neg_faces[axis]].mathematical_type == PERIODIC;
        const PetscBool pos_periodic =
            user->boundary_faces[pos_faces[axis]].mathematical_type == PERIODIC;

        PetscCheck(neg_periodic == pos_periodic, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
                   "Configuration Error: Periodic boundaries in the %c direction must be paired; "
                   "%s is %s while %s is %s.",
                   axis_names[axis],
                   BCFaceToString(neg_faces[axis]), neg_periodic ? "PERIODIC" : "not periodic",
                   BCFaceToString(pos_faces[axis]), pos_periodic ? "PERIODIC" : "not periodic");
        PetscCheck(dm_periodic[axis] == neg_periodic, PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,
                   "Configuration Error: The %c-direction DM periodic flag (%d) does not match "
                   "the paired boundary configuration (%s).",
                   axis_names[axis], (int)dm_periodic[axis], neg_periodic ? "PERIODIC" : "not periodic");
    }

    // --- Rule Set 2: Driven Flow Handler Consistency ---
    // This specialized validator will check all rules related to driven flow handlers.
    ierr = Validate_DrivenFlowConfiguration(user); CHKERRQ(ierr);

    // --- Rule Set 3: (Future Extension) Overset Interface Consistency ---
    // ierr = Validate_OversetConfiguration(user); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Boundary configuration is valid.\n");

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
 * @brief Implementation of \ref BoundarySystem_Initialize().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Boundaries.h`.
 * @see BoundarySystem_Initialize()
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

    // Step 1.1: Validate the parsed configuration to ensure there are no Boundary Condition conflicts
    ierr = BoundarySystem_Validate(user); CHKERRQ(ierr);

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

    LOG_ALLOW(GLOBAL, LOG_INFO, "All boundary handlers created and initialized successfully.\n");
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "PropagateBoundaryConfigToCoarserLevels"
/**
 * @brief Internal helper implementation: `PropagateBoundaryConfigToCoarserLevels()`.
 * @details Local to this translation unit.
 */
PetscErrorCode PropagateBoundaryConfigToCoarserLevels(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserMG *usermg = &simCtx->usermg;
    
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Propagating BC configuration from finest to coarser multigrid levels...\n");
    
    // Loop from second-finest down to coarsest
    for (PetscInt level = usermg->mglevels - 2; level >= 0; level--) {
        for (PetscInt bi = 0; bi < simCtx->block_number; bi++) {
            UserCtx *user_coarse = &usermg->mgctx[level].user[bi];
            UserCtx *user_fine   = &usermg->mgctx[level + 1].user[bi];
            
            LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "Rank %d: Copying BC config from level %d to level %d, block %d\n",
                          simCtx->rank, level + 1, level, bi);
            
            // Copy the 6 boundary face configurations
            for (int face_i = 0; face_i < 6; face_i++) {
                user_coarse->boundary_faces[face_i].face_id = user_fine->boundary_faces[face_i].face_id;
                user_coarse->boundary_faces[face_i].mathematical_type = user_fine->boundary_faces[face_i].mathematical_type;
                user_coarse->boundary_faces[face_i].handler_type = user_fine->boundary_faces[face_i].handler_type;
                
                // Copy parameter list (deep copy)
                FreeBC_ParamList(user_coarse->boundary_faces[face_i].params); // Clear any existing
                user_coarse->boundary_faces[face_i].params = NULL;
                
                BC_Param **dst_next = &user_coarse->boundary_faces[face_i].params;
                for (BC_Param *src = user_fine->boundary_faces[face_i].params; src; src = src->next) {
                    BC_Param *new_param;
                    ierr = PetscMalloc1(1, &new_param); CHKERRQ(ierr);
                    ierr = PetscStrallocpy(src->key, &new_param->key); CHKERRQ(ierr);
                    ierr = PetscStrallocpy(src->value, &new_param->value); CHKERRQ(ierr);
                    new_param->next = NULL;
                    *dst_next = new_param;
                    dst_next = &new_param->next;
                }
                
                // IMPORTANT: Do NOT create handler objects for coarser levels
                // Handlers are only needed at finest level for timestepping Apply() calls
                user_coarse->boundary_faces[face_i].handler = NULL;
            }
            
            // Propagate the particle inlet lookup fields to coarse levels as well.
            user_coarse->inletFaceDefined = user_fine->inletFaceDefined;
            user_coarse->identifiedInletBCFace = user_fine->identifiedInletBCFace;
        }
    }
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "BC configuration propagation complete.\n");
    
    PROFILE_FUNCTION_END;
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
 * @brief Implementation of \ref BoundarySystem_ExecuteStep().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Boundaries.h`.
 * @see BoundarySystem_ExecuteStep()
 */
PetscErrorCode BoundarySystem_ExecuteStep(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Starting.\n");
    
    // =========================================================================
    // PRIORITY 0: INLETS
    // =========================================================================
    
        PetscReal local_inflow_pre = 0.0;
        PetscReal local_inflow_post = 0.0;
        PetscReal global_inflow_pre = 0.0;
        PetscReal global_inflow_post = 0.0;
        PetscInt  num_handlers[3] = {0,0,0};
        
        LOG_ALLOW(LOCAL, LOG_TRACE, " (INLETS): Begin.\n");
        
        // Phase 1: PreStep - Preparation (e.g., calculate profiles, read files)
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_INLET) continue;
            if (!handler->PreStep) continue;
            
            num_handlers[0]++;
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
            if(!handler->Apply) continue; // For example Periodic BCs  
            
            num_handlers[1]++;

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
            
            num_handlers[2]++;

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
                  "  (INLETS): %d Prestep(s), %d Application(s), %d Poststep(s), FluxInSum = %.6e\n",
                  num_handlers[0],num_handlers[1],num_handlers[2], global_inflow_post);
    
    // =========================================================================
    // PRIORITY 1: FARFIELD
    // =========================================================================
    
        PetscReal local_farfield_in_pre = 0.0;
        PetscReal local_farfield_out_pre = 0.0;
        PetscReal local_farfield_in_post = 0.0;
        PetscReal local_farfield_out_post = 0.0;
        PetscReal global_farfield_in_pre = 0.0;
        PetscReal global_farfield_out_pre = 0.0;
        PetscReal global_farfield_in_post = 0.0;
        PetscReal global_farfield_out_post = 0.0;
        memset(num_handlers,0,sizeof(num_handlers));
        
        LOG_ALLOW(LOCAL, LOG_TRACE, "  (FARFIELD): Begin.\n");
        
        // Phase 1: PreStep - Analyze flow direction, measure initial flux
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_FARFIELD) continue;
            if (!handler->PreStep) continue;
            
            num_handlers[0]++;
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
            if(!handler->Apply) continue; // For example Periodic BCs            
            
            num_handlers[1]++;

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
            
            num_handlers[2]++;

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
                      "  (FARFIELD): %d Prestep(s), %d Application(s), %d Poststep(s) , InFlux=%.6e, OutFlux=%.6e\n",
                      num_handlers[0],num_handlers[1],num_handlers[2], global_farfield_in_post, global_farfield_out_post);
        } else {
            // No farfield handlers - zero out the fluxes
            user->simCtx->FarFluxInSum = 0.0;
            user->simCtx->FarFluxOutSum = 0.0;
        }
    
    
    // =========================================================================
    // PRIORITY 2: WALLS
    // =========================================================================
    
        memset(num_handlers,0,sizeof(num_handlers));
        
        LOG_ALLOW(LOCAL, LOG_TRACE, "  (WALLS): Begin.\n");
        
        // Phase 1: PreStep - Preparation (usually no-op for walls)
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_WALL) continue;
            if (!handler->PreStep) continue;
            
            num_handlers[0]++;
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
            if(!handler->Apply) continue; // For example Periodic BCs  
            
            num_handlers[1]++;

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
            
            num_handlers[2]++;

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
        
        LOG_ALLOW(GLOBAL, LOG_INFO, "  (WALLS): %d Prestep(s), %d Application(s), %d Poststep(s) applied.\n",
                  num_handlers[0],num_handlers[1],num_handlers[2]);
    
    
    // =========================================================================
    // PRIORITY 3: OUTLETS
    // =========================================================================
    
        PetscReal local_outflow_pre = 0.0;
        PetscReal local_outflow_post = 0.0;
        PetscReal global_outflow_pre = 0.0;
        PetscReal global_outflow_post = 0.0;
        memset(num_handlers,0,sizeof(num_handlers));
        
        LOG_ALLOW(LOCAL, LOG_TRACE, "  (OUTLETS): Begin.\n");
        
        // Phase 1: PreStep - Measure uncorrected outflow (from ucat)
        for (int i = 0; i < 6; i++) {
            BoundaryCondition *handler = user->boundary_faces[i].handler;
            if (!handler || handler->priority != BC_PRIORITY_OUTLET) continue;
            if (!handler->PreStep) continue;
            
            num_handlers[0]++;
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
            if(!handler->Apply) continue; // For example Periodic BCs  
            
            num_handlers[1]++;

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
            
            num_handlers[2]++;

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
        
        // Store for global reporting.
        user->simCtx->FluxOutSum = global_outflow_post;
        
        // Conservation check (compare total outflow vs total inflow)
        PetscReal total_outflow = global_outflow_post + user->simCtx->FarFluxOutSum;
        PetscReal flux_error = PetscAbsReal(total_outflow - total_inflow);
        PetscReal relative_error = (total_inflow > 1e-16) ? 
                                   flux_error / total_inflow : flux_error;
        
        LOG_ALLOW(GLOBAL, LOG_INFO, 
                  "  (OUTLETS): %d Prestep(s), %d Application(s), %d Poststep(s), FluxOutSum = %.6e\n",
                  num_handlers[0],num_handlers[1],num_handlers[2], global_outflow_post);
        LOG_ALLOW(GLOBAL, LOG_INFO, 
                  "    Conservation: Total In=%.6e, Total Out=%.6e, Error=%.3e (%.2e)%%)\n",
                  total_inflow, total_outflow, flux_error, relative_error * 100.0);
        
        if (relative_error > 1e-6) {
            LOG_ALLOW(GLOBAL, LOG_WARNING, 
                     "    WARNING: Large mass conservation error (%.2e%%)!\n",
                     relative_error * 100.0);
        }
    
    
    LOG_ALLOW(LOCAL, LOG_VERBOSE, "Complete.\n");
    
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

// =============================================================================
//
//                   PRIVATE "LIGHT" EXECUTION ENGINE
//
// =============================================================================

#undef __FUNCT__
#define __FUNCT__ "BoundarySystem_RefreshUbcs"
/**
 * @brief Internal helper implementation: `BoundarySystem_RefreshUbcs()`.
 * @details Local to this translation unit.
 */
PetscErrorCode BoundarySystem_RefreshUbcs(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    LOG_ALLOW(GLOBAL, LOG_TRACE, "Refreshing `ubcs` targets for flow-dependent boundaries...\n");

    // Loop through all 6 faces of the domain
    for (int i = 0; i < 6; i++) {
        BoundaryCondition *handler = user->boundary_faces[i].handler;
        
        // THE FILTER:
        // This is the core logic. We only act if a handler exists for the face
        // AND that handler has explicitly implemented the `UpdateUbcs` method.
        if (handler && handler->UpdateUbcs) {
            
            const char *face_name = BCFaceToString((BCFace)i);
            LOG_ALLOW(LOCAL, LOG_TRACE, "  Calling UpdateUbcs() for handler on Face %s.\n", face_name);

            // Prepare the context. For this refresh step, we don't need to pass flux sums.
            BCContext ctx = {
                .user = user,
                .face_id = (BCFace)i,
                .global_inflow_sum = NULL,
                .global_outflow_sum = NULL,
                .global_farfield_inflow_sum = NULL,
                .global_farfield_outflow_sum = NULL
            };
            
            // Call the handler's specific UpdateUbcs function pointer.
            ierr = handler->UpdateUbcs(handler, &ctx); CHKERRQ(ierr);
        }
    }

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
 * @brief Implementation of \ref BoundarySystem_Destroy().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Boundaries.h`.
 * @see BoundarySystem_Destroy()
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
#define __FUNCT__ "TransferPeriodicFieldByDirection"
/**
 * @brief Copies one cell field's wrapped local values onto the owned periodic duplicate plane.
 * @details Handles scalar and three-component storage for one selected logical axis.
 */
static PetscErrorCode TransferPeriodicFieldByDirection(UserCtx *user, FieldId field_id, char direction)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info = user->info;
    PetscInt       xs = info.xs, xe = info.xs + info.xm;
    PetscInt       ys = info.ys, ye = info.ys + info.ym;
    PetscInt       zs = info.zs, ze = info.zs + info.zm;
    PetscInt       mx = info.mx, my = info.my, mz = info.mz;

    FieldView field_view;
    DM        dm;
    Vec       global_vec;
    Vec       local_vec;
    PetscInt  dof;

    PetscFunctionBeginUser;
    PetscCall(FieldGetView(user, field_id, &field_view));
    PetscCheck((field_view.descriptor->capabilities & FIELD_CAPABILITY_PERIODIC_CELL_SYNC) != 0u,
               PETSC_COMM_SELF, PETSC_ERR_SUP,
               "Field '%s' is not registered for periodic cell synchronization.",
               field_view.descriptor->canonical_name);
    dm = field_view.dm;
    global_vec = field_view.global_vec;
    local_vec = field_view.local_vec;
    dof = field_view.descriptor->dof;

    // --- Execute the copy logic based on DoF and Direction ---
    if (dof == 1) { // --- Handle SCALAR fields (PetscReal) ---
        PetscReal ***g_array, ***l_array;
        ierr = DMDAVecGetArray(dm, global_vec, &g_array); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(dm, local_vec, (void*)&l_array); CHKERRQ(ierr); // Use Read for safety

        switch (direction) {
            case 'i':
                if (user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC && xs == 0) for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) g_array[k][j][xs] = l_array[k][j][xs-2];
                if (user->boundary_faces[BC_FACE_POS_X].mathematical_type == PERIODIC && xe == mx) for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) g_array[k][j][xe-1] = l_array[k][j][xe+1];
                break;
            case 'j':
                if (user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC && ys == 0) for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) g_array[k][ys][i] = l_array[k][ys-2][i];
                if (user->boundary_faces[BC_FACE_POS_Y].mathematical_type == PERIODIC && ye == my) for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) g_array[k][ye-1][i] = l_array[k][ye+1][i];
                break;
            case 'k':
                if (user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC && zs == 0) for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) g_array[zs][j][i] = l_array[zs-2][j][i];
                if (user->boundary_faces[BC_FACE_POS_Z].mathematical_type == PERIODIC && ze == mz) for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) g_array[ze-1][j][i] = l_array[ze+1][j][i];
                break;
            default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid direction '%c'", direction);
        }
        ierr = DMDAVecRestoreArray(dm, global_vec, &g_array); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(dm, local_vec, (void*)&l_array); CHKERRQ(ierr);

    } else if (dof == 3) { // --- Handle VECTOR fields (Cmpnts) ---
        Cmpnts ***g_array, ***l_array;
        ierr = DMDAVecGetArray(dm, global_vec, &g_array); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(dm, local_vec, (void*)&l_array); CHKERRQ(ierr);

        switch (direction) {
            case 'i':
                if (user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC && xs == 0) for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) g_array[k][j][xs] = l_array[k][j][xs-2];
                if (user->boundary_faces[BC_FACE_POS_X].mathematical_type == PERIODIC && xe == mx) for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) g_array[k][j][xe-1] = l_array[k][j][xe+1];
                break;
            case 'j':
                if (user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC && ys == 0) for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) g_array[k][ys][i] = l_array[k][ys-2][i];
                if (user->boundary_faces[BC_FACE_POS_Y].mathematical_type == PERIODIC && ye == my) for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) g_array[k][ye-1][i] = l_array[k][ye+1][i];
                break;
            case 'k':
                if (user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC && zs == 0) for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) g_array[zs][j][i] = l_array[zs-2][j][i];
                if (user->boundary_faces[BC_FACE_POS_Z].mathematical_type == PERIODIC && ze == mz) for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) g_array[ze-1][j][i] = l_array[ze+1][j][i];
                break;
            default: SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Invalid direction '%c'", direction);
        }
        ierr = DMDAVecRestoreArray(dm, global_vec, &g_array); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArrayRead(dm, local_vec, (void*)&l_array); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SynchronizePeriodicCellFields"
/**
 * @brief Implementation of \ref SynchronizePeriodicCellFields().
 * @details Full API contract is documented with the header declaration in
 *          `include/Boundaries.h`.
 */
PetscErrorCode SynchronizePeriodicCellFields(UserCtx *user, PetscInt num_fields, const FieldId field_ids[])
{
    PetscErrorCode ierr;
    PetscBool      periodic_i;
    PetscBool      periodic_j;
    PetscBool      periodic_k;

    PetscFunctionBeginUser;

    PetscCheck(num_fields >= 0, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Number of cell fields cannot be negative.");
    if (num_fields == 0) PetscFunctionReturn(0);
    PetscCheck(field_ids != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Cell field-ID array cannot be NULL.");

    periodic_i =
        user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC ||
        user->boundary_faces[BC_FACE_POS_X].mathematical_type == PERIODIC;
    periodic_j =
        user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC ||
        user->boundary_faces[BC_FACE_POS_Y].mathematical_type == PERIODIC;
    periodic_k =
        user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC ||
        user->boundary_faces[BC_FACE_POS_Z].mathematical_type == PERIODIC;

    if (!periodic_i && !periodic_j && !periodic_k) PetscFunctionReturn(0);

    for (PetscInt field = 0; field < num_fields; field++) {
        ierr = UpdateLocalGhosts(user, field_ids[field]); CHKERRQ(ierr);
    }

    if (periodic_i) {
        for (PetscInt field = 0; field < num_fields; field++) {
            ierr = TransferPeriodicFieldByDirection(user, field_ids[field], 'i'); CHKERRQ(ierr);
        }
        for (PetscInt field = 0; field < num_fields; field++) {
            ierr = UpdateLocalGhosts(user, field_ids[field]); CHKERRQ(ierr);
        }
    }

    if (periodic_j) {
        for (PetscInt field = 0; field < num_fields; field++) {
            ierr = TransferPeriodicFieldByDirection(user, field_ids[field], 'j'); CHKERRQ(ierr);
        }
        for (PetscInt field = 0; field < num_fields; field++) {
            ierr = UpdateLocalGhosts(user, field_ids[field]); CHKERRQ(ierr);
        }
    }

    if (periodic_k) {
        for (PetscInt field = 0; field < num_fields; field++) {
            ierr = TransferPeriodicFieldByDirection(user, field_ids[field], 'k'); CHKERRQ(ierr);
        }
        for (PetscInt field = 0; field < num_fields; field++) {
            ierr = UpdateLocalGhosts(user, field_ids[field]); CHKERRQ(ierr);
        }
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GetPersistentFaceField"
/**
 * @brief Resolves one registered persistent single-face-family field.
 */
static PetscErrorCode GetPersistentFaceField(UserCtx *user, FieldId field_id,
                                             char face_direction, DM *dm,
                                             Vec *global_vec, Vec *local_vec,
                                             PetscInt *dof)
{
    FieldView   field_view;
    FieldLayout expected_layout;

    PetscFunctionBeginUser;
    PetscCheck(face_direction == 'i' || face_direction == 'j' || face_direction == 'k',
               PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
               "Invalid face direction '%c'; expected 'i', 'j', or 'k'.", face_direction);

    expected_layout = face_direction == 'i' ? FIELD_LAYOUT_I_FACE :
                      (face_direction == 'j' ? FIELD_LAYOUT_J_FACE : FIELD_LAYOUT_K_FACE);
    PetscCall(FieldGetView(user, field_id, &field_view));
    PetscCheck((field_view.descriptor->capabilities & FIELD_CAPABILITY_PERIODIC_FACE_SYNC) != 0u &&
               field_view.descriptor->layout == expected_layout,
               PETSC_COMM_SELF, PETSC_ERR_SUP,
               "Field '%s' is not registered for %c-face periodic synchronization.",
               field_view.descriptor->canonical_name, face_direction);

    *dm = field_view.dm;
    *global_vec = field_view.global_vec;
    *local_vec = field_view.local_vec;
    *dof = field_view.descriptor->dof;
    PetscFunctionReturn(0);
}

/**
 * @brief Returns whether a registered face field stores physical coordinates.
 */
static PetscErrorCode IsFaceCenterCoordinateField(FieldId field_id, PetscBool *is_coordinate)
{
    const FieldDescriptor *descriptor = NULL;

    PetscFunctionBeginUser;
    PetscCheck(is_coordinate != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Face-coordinate result cannot be NULL.");
    PetscCall(FieldGetDescriptor(field_id, &descriptor));
    *is_coordinate = (PetscBool)((descriptor->capabilities & FIELD_CAPABILITY_PERIODIC_GEOMETRY_SHIFT) != 0u);
    PetscFunctionReturn(0);
}

/**
 * @brief Applies geometric translations to wrapped face-center ghost coordinates.
 */
static PetscErrorCode TranslatePeriodicFaceCenterGhosts(UserCtx *user, Vec local_vec)
{
    DMDALocalInfo info;
    Cmpnts ***array;
    const BCFace negative_faces[3] = {BC_FACE_NEG_X, BC_FACE_NEG_Y, BC_FACE_NEG_Z};
    const BCFace positive_faces[3] = {BC_FACE_POS_X, BC_FACE_POS_Y, BC_FACE_POS_Z};

    PetscFunctionBeginUser;
    PetscCall(DMDAGetLocalInfo(user->fda, &info));
    PetscCall(DMDAVecGetArray(user->fda, local_vec, &array));

    for (PetscInt axis = 0; axis < 3; axis++) {
        const PetscBool active =
            user->boundary_faces[negative_faces[axis]].mathematical_type == PERIODIC ||
            user->boundary_faces[positive_faces[axis]].mathematical_type == PERIODIC;
        Cmpnts translation;

        if (!active) continue;
        PetscCheck(user->periodic_translation_valid[axis], PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
                   "Periodic face-center synchronization requires validated %c-direction geometry.",
                   "XYZ"[axis]);
        translation = user->periodic_translation[axis];

        for (PetscInt k = info.gzs; k < info.gzs + info.gzm; k++) {
            for (PetscInt j = info.gys; j < info.gys + info.gym; j++) {
                for (PetscInt i = info.gxs; i < info.gxs + info.gxm; i++) {
                    PetscReal scale = 0.0;
                    const PetscInt index = axis == 0 ? i : (axis == 1 ? j : k);
                    const PetscInt size = axis == 0 ? info.mx : (axis == 1 ? info.my : info.mz);
                    if (index < 0) scale = -1.0;
                    else if (index >= size) scale = 1.0;
                    if (scale == 0.0) continue;
                    array[k][j][i].x += scale * translation.x;
                    array[k][j][i].y += scale * translation.y;
                    array[k][j][i].z += scale * translation.z;
                }
            }
        }
    }

    PetscCall(DMDAVecRestoreArray(user->fda, local_vec, &array));
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TransferPeriodicFaceFieldByDirection"
/** @brief Transfers one registered face-family field along one periodic axis. */
static PetscErrorCode TransferPeriodicFaceFieldByDirection(UserCtx *user, FieldId field_id,
                                                           char face_direction, char periodic_direction)
{
    PetscErrorCode ierr;
    DMDALocalInfo info = user->info;
    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;
    PetscInt mx = info.mx, my = info.my, mz = info.mz;
    DM dm;
    Vec global_vec, local_vec;
    PetscInt dof;

    PetscFunctionBeginUser;
    PetscCheck(periodic_direction == 'i' || periodic_direction == 'j' || periodic_direction == 'k',
               PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
               "Invalid periodic direction '%c'; expected 'i', 'j', or 'k'.", periodic_direction);
    PetscCall(GetPersistentFaceField(user, field_id, face_direction, &dm, &global_vec, &local_vec, &dof));

    if (dof == 1) {
        PetscReal ***global_array, ***local_array;
        ierr = DMDAVecGetArray(dm, global_vec, &global_array); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(dm, local_vec, &local_array); CHKERRQ(ierr);

        if (periodic_direction == 'i') {
            if (user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC && xs == 0)
                for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) global_array[k][j][0] = local_array[k][j][-2];
            if (user->boundary_faces[BC_FACE_POS_X].mathematical_type == PERIODIC && xe == mx)
                for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) global_array[k][j][mx-1] = local_array[k][j][mx+1];
        } else if (periodic_direction == 'j') {
            if (user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC && ys == 0)
                for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) global_array[k][0][i] = local_array[k][-2][i];
            if (user->boundary_faces[BC_FACE_POS_Y].mathematical_type == PERIODIC && ye == my)
                for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) global_array[k][my-1][i] = local_array[k][my+1][i];
        } else {
            if (user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC && zs == 0)
                for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) global_array[0][j][i] = local_array[-2][j][i];
            if (user->boundary_faces[BC_FACE_POS_Z].mathematical_type == PERIODIC && ze == mz)
                for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) global_array[mz-1][j][i] = local_array[mz+1][j][i];
        }
        ierr = DMDAVecRestoreArrayRead(dm, local_vec, &local_array); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(dm, global_vec, &global_array); CHKERRQ(ierr);
    } else {
        Cmpnts ***global_array, ***local_array;
        ierr = DMDAVecGetArray(dm, global_vec, &global_array); CHKERRQ(ierr);
        ierr = DMDAVecGetArrayRead(dm, local_vec, &local_array); CHKERRQ(ierr);

        if (periodic_direction == 'i') {
            if (user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC && xs == 0)
                for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) global_array[k][j][0] = local_array[k][j][-2];
            if (user->boundary_faces[BC_FACE_POS_X].mathematical_type == PERIODIC && xe == mx)
                for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) global_array[k][j][mx-1] = local_array[k][j][mx+1];
        } else if (periodic_direction == 'j') {
            if (user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC && ys == 0)
                for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) global_array[k][0][i] = local_array[k][-2][i];
            if (user->boundary_faces[BC_FACE_POS_Y].mathematical_type == PERIODIC && ye == my)
                for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) global_array[k][my-1][i] = local_array[k][my+1][i];
        } else {
            if (user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC && zs == 0)
                for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) global_array[0][j][i] = local_array[-2][j][i];
            if (user->boundary_faces[BC_FACE_POS_Z].mathematical_type == PERIODIC && ze == mz)
                for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) global_array[mz-1][j][i] = local_array[mz+1][j][i];
        }
        ierr = DMDAVecRestoreArrayRead(dm, local_vec, &local_array); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(dm, global_vec, &global_array); CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SynchronizePeriodicFaceFields"
// Implements SynchronizePeriodicFaceFields(); the public header owns the
// rendered API contract.
PetscErrorCode SynchronizePeriodicFaceFields(UserCtx *user, char face_direction,
                                             PetscInt num_fields, const FieldId field_ids[])
{
    PetscErrorCode ierr;
    const char periodic_directions[3] = {'i', 'j', 'k'};
    const BCFace negative_faces[3] = {BC_FACE_NEG_X, BC_FACE_NEG_Y, BC_FACE_NEG_Z};
    const BCFace positive_faces[3] = {BC_FACE_POS_X, BC_FACE_POS_Y, BC_FACE_POS_Z};

    PetscFunctionBeginUser;
    PetscCheck(num_fields >= 0, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Number of face fields cannot be negative.");
    if (num_fields == 0) PetscFunctionReturn(0);
    PetscCheck(field_ids != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Face field-ID array cannot be NULL.");

    for (PetscInt field = 0; field < num_fields; field++) {
        DM dm;
        Vec global_vec, local_vec;
        PetscInt dof;
        PetscBool is_coordinate;
        PetscCall(GetPersistentFaceField(user, field_ids[field], face_direction,
                                         &dm, &global_vec, &local_vec, &dof));
        ierr = UpdateLocalGhosts(user, field_ids[field]); CHKERRQ(ierr);
        ierr = IsFaceCenterCoordinateField(field_ids[field], &is_coordinate); CHKERRQ(ierr);
        if (is_coordinate) {
            ierr = TranslatePeriodicFaceCenterGhosts(user, local_vec); CHKERRQ(ierr);
        }
    }

    for (PetscInt direction = 0; direction < 3; direction++) {
        const PetscBool active =
            user->boundary_faces[negative_faces[direction]].mathematical_type == PERIODIC ||
            user->boundary_faces[positive_faces[direction]].mathematical_type == PERIODIC;
        if (!active) continue;

        for (PetscInt field = 0; field < num_fields; field++) {
            ierr = TransferPeriodicFaceFieldByDirection(user, field_ids[field], face_direction,
                                                        periodic_directions[direction]); CHKERRQ(ierr);
        }
        for (PetscInt field = 0; field < num_fields; field++) {
            DM dm;
            Vec global_vec, local_vec;
            PetscInt dof;
            PetscBool is_coordinate;
            ierr = UpdateLocalGhosts(user, field_ids[field]); CHKERRQ(ierr);
            ierr = IsFaceCenterCoordinateField(field_ids[field], &is_coordinate); CHKERRQ(ierr);
            if (is_coordinate) {
                PetscCall(GetPersistentFaceField(user, field_ids[field], face_direction,
                                                 &dm, &global_vec, &local_vec, &dof));
                ierr = TranslatePeriodicFaceCenterGhosts(user, local_vec); CHKERRQ(ierr);
            }
        }
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "GetPersistentStaggeredField"
/**
 * @brief Resolves one registered persistent component-staggered field.
 */
static PetscErrorCode GetPersistentStaggeredField(UserCtx *user, FieldId field_id,
                                                  DM *dm, Vec *global_vec, Vec *local_vec)
{
    FieldView field_view;

    PetscFunctionBeginUser;
    PetscCall(FieldGetView(user, field_id, &field_view));
    PetscCheck((field_view.descriptor->capabilities & FIELD_CAPABILITY_PERIODIC_STAGGERED_SYNC) != 0u &&
               field_view.descriptor->layout == FIELD_LAYOUT_COMPONENT_STAGGERED,
               PETSC_COMM_SELF, PETSC_ERR_SUP,
               "Field '%s' is not registered for component-staggered periodic synchronization.",
               field_view.descriptor->canonical_name);
    *dm = field_view.dm;
    *global_vec = field_view.global_vec;
    *local_vec = field_view.local_vec;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TransferPeriodicStaggeredFieldByDirection"
/**
 * @brief Transfers one component-staggered field along one periodic axis.
 */
static PetscErrorCode TransferPeriodicStaggeredFieldByDirection(UserCtx *user, FieldId field_id,
                                                                char periodic_direction)
{
    PetscErrorCode ierr;
    DMDALocalInfo info = user->info;
    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;
    PetscInt mx = info.mx, my = info.my, mz = info.mz;
    DM dm;
    Vec global_vec, local_vec;
    Cmpnts ***global_array, ***local_array;

    PetscFunctionBeginUser;
    PetscCheck(periodic_direction == 'i' || periodic_direction == 'j' || periodic_direction == 'k',
               PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG,
               "Invalid periodic direction '%c'; expected 'i', 'j', or 'k'.", periodic_direction);
    PetscCall(GetPersistentStaggeredField(user, field_id, &dm, &global_vec, &local_vec));

    ierr = DMDAVecGetArray(dm, global_vec, &global_array); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(dm, local_vec, &local_array); CHKERRQ(ierr);

    if (periodic_direction == 'i') {
        if (user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC && xs == 0)
            for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) global_array[k][j][0] = local_array[k][j][-2];
        if (user->boundary_faces[BC_FACE_POS_X].mathematical_type == PERIODIC && xe == mx)
            for (PetscInt k=zs; k<ze; k++) for (PetscInt j=ys; j<ye; j++) global_array[k][j][mx-1] = local_array[k][j][mx+1];
    } else if (periodic_direction == 'j') {
        if (user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC && ys == 0)
            for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) global_array[k][0][i] = local_array[k][-2][i];
        if (user->boundary_faces[BC_FACE_POS_Y].mathematical_type == PERIODIC && ye == my)
            for (PetscInt k=zs; k<ze; k++) for (PetscInt i=xs; i<xe; i++) global_array[k][my-1][i] = local_array[k][my+1][i];
    } else {
        if (user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC && zs == 0)
            for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) global_array[0][j][i] = local_array[-2][j][i];
        if (user->boundary_faces[BC_FACE_POS_Z].mathematical_type == PERIODIC && ze == mz)
            for (PetscInt j=ys; j<ye; j++) for (PetscInt i=xs; i<xe; i++) global_array[mz-1][j][i] = local_array[mz+1][j][i];
    }

    ierr = DMDAVecRestoreArrayRead(dm, local_vec, &local_array); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(dm, global_vec, &global_array); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SynchronizePeriodicStaggeredFields"
/**
 * @brief Implementation of \ref SynchronizePeriodicStaggeredFields().
 */
PetscErrorCode SynchronizePeriodicStaggeredFields(UserCtx *user, PetscInt num_fields,
                                                  const FieldId field_ids[])
{
    PetscErrorCode ierr;
    const char periodic_directions[3] = {'i', 'j', 'k'};
    const BCFace negative_faces[3] = {BC_FACE_NEG_X, BC_FACE_NEG_Y, BC_FACE_NEG_Z};
    const BCFace positive_faces[3] = {BC_FACE_POS_X, BC_FACE_POS_Y, BC_FACE_POS_Z};

    PetscFunctionBeginUser;
    PetscCheck(num_fields >= 0, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Number of staggered fields cannot be negative.");
    if (num_fields == 0) PetscFunctionReturn(0);
    PetscCheck(field_ids != NULL, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Staggered field-ID array cannot be NULL.");

    for (PetscInt field = 0; field < num_fields; field++) {
        DM dm;
        Vec global_vec, local_vec;
        PetscCall(GetPersistentStaggeredField(user, field_ids[field], &dm, &global_vec, &local_vec));
        ierr = UpdateLocalGhosts(user, field_ids[field]); CHKERRQ(ierr);
    }

    for (PetscInt direction = 0; direction < 3; direction++) {
        const PetscBool active =
            user->boundary_faces[negative_faces[direction]].mathematical_type == PERIODIC ||
            user->boundary_faces[positive_faces[direction]].mathematical_type == PERIODIC;
        if (!active) continue;

        for (PetscInt field = 0; field < num_fields; field++) {
            ierr = TransferPeriodicStaggeredFieldByDirection(user, field_ids[field],
                                                             periodic_directions[direction]); CHKERRQ(ierr);
        }
        for (PetscInt field = 0; field < num_fields; field++) {
            ierr = UpdateLocalGhosts(user, field_ids[field]); CHKERRQ(ierr);
        }
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PreparePeriodicQuickStencilFields"
/**
 * @brief Implementation of \ref PreparePeriodicQuickStencilFields().
 */
PetscErrorCode PreparePeriodicQuickStencilFields(UserCtx *user, Vec local_vector_field,
                                                 Vec local_scalar_field)
{
    DMDALocalInfo info = user->info;
    Cmpnts ***vector_array;
    PetscReal ***scalar_array;
    const PetscInt xs = info.xs, xe = info.xs + info.xm;
    const PetscInt ys = info.ys, ye = info.ys + info.ym;
    const PetscInt zs = info.zs, ze = info.zs + info.zm;
    const PetscInt gxs = info.gxs, gxe = info.gxs + info.gxm;
    const PetscInt gys = info.gys, gye = info.gys + info.gym;
    const PetscInt gzs = info.gzs, gze = info.gzs + info.gzm;

    PetscFunctionBeginUser;
    PetscCheck(local_vector_field && local_scalar_field, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "QUICK stencil repair requires both local vector and scalar fields.");
    PetscCall(DMDAVecGetArray(user->fda, local_vector_field, &vector_array));
    PetscCall(DMDAVecGetArray(user->da, local_scalar_field, &scalar_array));

    if (user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC && xs == 0) {
        for (PetscInt k = gzs; k < gze; k++) for (PetscInt j = gys; j < gye; j++) {
            vector_array[k][j][-1] = vector_array[k][j][-3];
            scalar_array[k][j][-1] = scalar_array[k][j][-3];
        }
    }
    if (user->boundary_faces[BC_FACE_POS_X].mathematical_type == PERIODIC && xe == info.mx) {
        for (PetscInt k = gzs; k < gze; k++) for (PetscInt j = gys; j < gye; j++) {
            vector_array[k][j][info.mx] = vector_array[k][j][info.mx + 2];
            scalar_array[k][j][info.mx] = scalar_array[k][j][info.mx + 2];
        }
    }
    if (user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC && ys == 0) {
        for (PetscInt k = gzs; k < gze; k++) for (PetscInt i = gxs; i < gxe; i++) {
            vector_array[k][-1][i] = vector_array[k][-3][i];
            scalar_array[k][-1][i] = scalar_array[k][-3][i];
        }
    }
    if (user->boundary_faces[BC_FACE_POS_Y].mathematical_type == PERIODIC && ye == info.my) {
        for (PetscInt k = gzs; k < gze; k++) for (PetscInt i = gxs; i < gxe; i++) {
            vector_array[k][info.my][i] = vector_array[k][info.my + 2][i];
            scalar_array[k][info.my][i] = scalar_array[k][info.my + 2][i];
        }
    }
    if (user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC && zs == 0) {
        for (PetscInt j = gys; j < gye; j++) for (PetscInt i = gxs; i < gxe; i++) {
            vector_array[-1][j][i] = vector_array[-3][j][i];
            scalar_array[-1][j][i] = scalar_array[-3][j][i];
        }
    }
    if (user->boundary_faces[BC_FACE_POS_Z].mathematical_type == PERIODIC && ze == info.mz) {
        for (PetscInt j = gys; j < gye; j++) for (PetscInt i = gxs; i < gxe; i++) {
            vector_array[info.mz][j][i] = vector_array[info.mz + 2][j][i];
            scalar_array[info.mz][j][i] = scalar_array[info.mz + 2][j][i];
        }
    }

    PetscCall(DMDAVecRestoreArray(user->da, local_scalar_field, &scalar_array));
    PetscCall(DMDAVecRestoreArray(user->fda, local_vector_field, &vector_array));
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "SynchronizePeriodicLocalStaggeredField"
/**
 * @brief Implementation of \ref SynchronizePeriodicLocalStaggeredField().
 */
PetscErrorCode SynchronizePeriodicLocalStaggeredField(UserCtx *user, Vec local_field)
{
    DMDALocalInfo info = user->info;
    Cmpnts ***array;
    const PetscInt xs = info.xs, xe = info.xs + info.xm;
    const PetscInt ys = info.ys, ye = info.ys + info.ym;
    const PetscInt zs = info.zs, ze = info.zs + info.zm;

    PetscFunctionBeginUser;
    PetscCheck(local_field, PETSC_COMM_SELF, PETSC_ERR_ARG_NULL,
               "Local staggered field cannot be NULL.");
    PetscCall(DMLocalToLocalBegin(user->fda, local_field, INSERT_VALUES, local_field));
    PetscCall(DMLocalToLocalEnd(user->fda, local_field, INSERT_VALUES, local_field));
    PetscCall(DMDAVecGetArray(user->fda, local_field, &array));

    if (user->boundary_faces[BC_FACE_NEG_X].mathematical_type == PERIODIC && xs == 0)
        for (PetscInt k = zs; k < ze; k++) for (PetscInt j = ys; j < ye; j++) array[k][j][0].x = array[k][j][-2].x;
    if (user->boundary_faces[BC_FACE_POS_X].mathematical_type == PERIODIC && xe == info.mx)
        for (PetscInt k = zs; k < ze; k++) for (PetscInt j = ys; j < ye; j++) array[k][j][info.mx - 1].x = array[k][j][info.mx + 1].x;
    if (user->boundary_faces[BC_FACE_NEG_Y].mathematical_type == PERIODIC && ys == 0)
        for (PetscInt k = zs; k < ze; k++) for (PetscInt i = xs; i < xe; i++) array[k][0][i].y = array[k][-2][i].y;
    if (user->boundary_faces[BC_FACE_POS_Y].mathematical_type == PERIODIC && ye == info.my)
        for (PetscInt k = zs; k < ze; k++) for (PetscInt i = xs; i < xe; i++) array[k][info.my - 1][i].y = array[k][info.my + 1][i].y;
    if (user->boundary_faces[BC_FACE_NEG_Z].mathematical_type == PERIODIC && zs == 0)
        for (PetscInt j = ys; j < ye; j++) for (PetscInt i = xs; i < xe; i++) array[0][j][i].z = array[-2][j][i].z;
    if (user->boundary_faces[BC_FACE_POS_Z].mathematical_type == PERIODIC && ze == info.mz)
        for (PetscInt j = ys; j < ye; j++) for (PetscInt i = xs; i < xe; i++) array[info.mz - 1][j][i].z = array[info.mz + 1][j][i].z;

    PetscCall(DMDAVecRestoreArray(user->fda, local_field, &array));
    PetscCall(DMLocalToLocalBegin(user->fda, local_field, INSERT_VALUES, local_field));
    PetscCall(DMLocalToLocalEnd(user->fda, local_field, INSERT_VALUES, local_field));
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ApplyMetricsPeriodicBCs"
/**
 * @brief Internal helper implementation: `ApplyMetricsPeriodicBCs()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ApplyMetricsPeriodicBCs(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    const FieldId cell_fields[] = {FIELD_ID_AJ};
    const FieldId i_face_fields[] = {FIELD_ID_CENTX, FIELD_ID_CSI, FIELD_ID_ICSI,
                                     FIELD_ID_IETA, FIELD_ID_IZET, FIELD_ID_IAJ};
    const FieldId j_face_fields[] = {FIELD_ID_CENTY, FIELD_ID_ETA, FIELD_ID_JCSI,
                                     FIELD_ID_JETA, FIELD_ID_JZET, FIELD_ID_JAJ};
    const FieldId k_face_fields[] = {FIELD_ID_CENTZ, FIELD_ID_ZET, FIELD_ID_KCSI,
                                     FIELD_ID_KETA, FIELD_ID_KZET, FIELD_ID_KAJ};

    ierr = SynchronizePeriodicCellFields(user, 1, cell_fields); CHKERRQ(ierr);
    ierr = SynchronizePeriodicFaceFields(user, 'i', 6, i_face_fields); CHKERRQ(ierr);
    ierr = SynchronizePeriodicFaceFields(user, 'j', 6, j_face_fields); CHKERRQ(ierr);
    ierr = SynchronizePeriodicFaceFields(user, 'k', 6, k_face_fields); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ApplyPeriodicBCs"
/**
 * @brief Internal helper implementation: `ApplyPeriodicBCs()`.
 * @details Local to this translation unit.
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

    LOG_ALLOW(GLOBAL, LOG_TRACE, "Applying periodic boundary conditions for all fields.\n");

    // STEP 1: Synchronize periodic cell-centered fields in deterministic direction order.
    const FieldId cell_fields[] = {FIELD_ID_UCAT, FIELD_ID_P, FIELD_ID_NVERT};
    ierr = SynchronizePeriodicCellFields(user, 3, cell_fields); CHKERRQ(ierr);

    /* A future temperature field must be catalogued before requesting its typed ghost update. */

    // STEP 2: Synchronize persistent staggered endpoints and repair local
    // component-normal ghosts through UpdateLocalGhosts().
    const FieldId staggered_fields[] = {FIELD_ID_UCONT};
    ierr = SynchronizePeriodicStaggeredFields(user, 1, staggered_fields); CHKERRQ(ierr);

    // FUTURE EXTENSION: Add new cell fields through SynchronizePeriodicCellFields().
    /*
    if (user->solve_temperature) {
        const char *temperature_field[] = {"Temperature"};
        ierr = SynchronizePeriodicCellFields(user, 1, temperature_field); CHKERRQ(ierr);
    }
    */

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "UpdateDummyCells"
/**
 * @brief Internal helper implementation: `UpdateDummyCells()`.
 * @details Local to this translation unit.
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
 * @brief Internal helper implementation: `UpdateCornerNodes()`.
 * @details Local to this translation unit.
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

/**
 * @brief Applies the configured wall model at one near-wall cell.
 *
 * The six face cases differ only in which indices they address, so the choice of law
 * is made here once rather than repeated at each of them. Every model writes the
 * corrected near-wall velocity and the friction velocity through the same two
 * pointers, which is what lets one call site serve all of them.
 *
 * Cabot's non-equilibrium form takes a pressure gradient, which is why the caller
 * supplies the read-only pressure array: the gradient is evaluated at the cell rather
 * than assumed zero, since zero would silently reduce the model to its equilibrium form.
 * The other two laws ignore that array.
 */
static PetscErrorCode ApplyWallModelAtCell(UserCtx *user, PetscReal ***pressure,
                                 PetscReal ***wall_eddy_viscosity,
                                 PetscInt i, PetscInt j, PetscInt k,
                                 PetscReal roughness_height,
                                 PetscReal distance_reference, PetscReal distance_boundary,
                                 Cmpnts velocity_wall, Cmpnts velocity_reference,
                                 Cmpnts *velocity_boundary, PetscReal *friction_velocity,
                                 PetscReal normal_x, PetscReal normal_y, PetscReal normal_z)
{
    PetscFunctionBeginUser;
    switch ((WallFunctionModel)user->simCtx->wallfunction) {
    case WALL_FUNCTION_WERNER:
        wall_function(user, distance_reference, distance_boundary,
                      velocity_wall, velocity_reference, velocity_boundary,
                      friction_velocity, normal_x, normal_y, normal_z);
        break;
    case WALL_FUNCTION_CABOT: {
        /* Cabot's departure from an equilibrium profile is driven by the pressure
           gradient in the wall layer, so it is supplied from the resolved pressure at
           this cell rather than assumed zero. Zero would silently reduce the model to
           its equilibrium form, which is not the model the user selected. The pressure
           is the previous projection's, which is the same lag every other explicit use
           of it carries. */
        Cmpnts pressure_gradient = {0.0, 0.0, 0.0};

        if (pressure) {
            PetscCall(ComputeScalarFieldDerivatives(user, i, j, k, pressure, &pressure_gradient));
        }
        wall_function_Cabot(user, roughness_height, distance_reference, distance_boundary,
                            velocity_wall, velocity_reference, velocity_boundary,
                            friction_velocity, normal_x, normal_y, normal_z,
                            pressure_gradient.x, pressure_gradient.y, pressure_gradient.z, 0);
        break;
    }
    case WALL_FUNCTION_LOG_LAW:
    default:
        wall_function_loglaw(user, roughness_height, distance_reference, distance_boundary,
                             velocity_wall, velocity_reference, velocity_boundary,
                             friction_velocity, normal_x, normal_y, normal_z);
        break;
    }

    /* Every face case reaches the model through here, so this is the one place that sees
       both the wall distance and the friction velocity the law produced. y+ is formed
       here for that reason: nothing downstream still holds the distance. */
    {
        WallModelDiagnosticsState *diagnostics = &user->wall_diagnostics;
        const PetscReal            molecular   = 1.0 / user->simCtx->ren;
        const PetscReal            u_tau       = *friction_velocity;
        const PetscReal            y_plus      = u_tau * distance_boundary / molecular;

        /* The viscous operator reaches the wall through a viscosity times a velocity
           gradient, so the modelled stress arrives only if the viscosity is the one that
           reproduces it across this cell: nu_eff = tau_w y / u. Molecular viscosity alone
           delivers a fraction u+/y+ of the stress, which is where most of it was going.
           Formed here because this is the one place holding the stress, the distance and
           the corrected speed together; nothing downstream still has all three. */
        if (wall_eddy_viscosity) {
            const Cmpnts    relative = {velocity_boundary->x - velocity_wall.x,
                                        velocity_boundary->y - velocity_wall.y,
                                        velocity_boundary->z - velocity_wall.z};
            const PetscReal normal_component = relative.x * normal_x + relative.y * normal_y +
                                               relative.z * normal_z;
            const PetscReal tangential[3] = {relative.x - normal_component * normal_x,
                                             relative.y - normal_component * normal_y,
                                             relative.z - normal_component * normal_z};
            const PetscReal speed = PetscSqrtReal(tangential[0] * tangential[0] +
                                                  tangential[1] * tangential[1] +
                                                  tangential[2] * tangential[2]);

            /* A cell the model brought to rest carries no stress to deliver, and the
               quotient below would be meaningless there. */
            if (speed > 1.0e-10) {
                const PetscReal effective = u_tau * u_tau * distance_boundary / speed;

                wall_eddy_viscosity[k][j][i] = PetscMax(effective - molecular, 0.0);
                diagnostics->wall_viscosity_sum += wall_eddy_viscosity[k][j][i];
            } else {
                wall_eddy_viscosity[k][j][i] = 0.0;
            }
        }

        if (diagnostics->cells == 0) {
            diagnostics->friction_velocity_min = u_tau;
            diagnostics->friction_velocity_max = u_tau;
            diagnostics->y_plus_max            = y_plus;
        } else {
            diagnostics->friction_velocity_min = PetscMin(diagnostics->friction_velocity_min, u_tau);
            diagnostics->friction_velocity_max = PetscMax(diagnostics->friction_velocity_max, u_tau);
            diagnostics->y_plus_max            = PetscMax(diagnostics->y_plus_max, y_plus);
        }
        diagnostics->friction_velocity_sum += u_tau;
        diagnostics->friction_velocity_sq  += u_tau * u_tau;
        diagnostics->wall_distance_sum     += distance_boundary;
        diagnostics->y_plus_sum            += y_plus;
        diagnostics->cells                 += 1;
    }
    PetscFunctionReturn(PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "ApplyWallFunction"
/**
 * @brief Internal helper implementation: `ApplyWallFunction()`.
 * @details Local to this translation unit.
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

    /* One pass, one sample set: the state describes this pass and nothing earlier. */
    ierr = PetscMemzero(&user->wall_diagnostics, sizeof(user->wall_diagnostics)); CHKERRQ(ierr);
    /* Cells the model does not reach must not keep a stale wall viscosity from an
       earlier pass, so the field is cleared rather than accumulated into. */
    ierr = VecSet(user->Nu_Wall, 0.0); CHKERRQ(ierr);
    
    // =========================================================================
    // STEP 1: Get read/write access to all necessary field arrays
    // =========================================================================
    Cmpnts ***velocity_cartesian;           // Cartesian velocity (modified)
    Cmpnts ***velocity_contravariant;       // Contravariant velocity (set to zero at walls)
    Cmpnts ***velocity_boundary;            // Boundary condition velocity (kept at zero)
    Cmpnts ***csi, ***eta, ***zet;         // Metric tensor components (face normals)
    PetscReal ***node_vertex_flag;          // Fluid/solid indicator (0=fluid, 1=solid)
    PetscReal ***cell_jacobian;             // Grid Jacobian (1/volume)
    PetscReal ***wall_eddy_viscosity;    // Effective wall eddy viscosity (written)
    PetscReal ***friction_velocity;
    PetscReal ***wall_pressure = NULL;         // u_tau (friction velocity field)
    
    ierr = DMDAVecGetArray(user->fda, user->Ucat, &velocity_cartesian); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Ucont, &velocity_contravariant); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Bcs.Ubcs, &velocity_boundary); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi, (const Cmpnts***)&csi); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta, (const Cmpnts***)&eta); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet, (const Cmpnts***)&zet); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da, user->lNvert, (const PetscReal***)&node_vertex_flag); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da, user->lAj, (const PetscReal***)&cell_jacobian); CHKERRQ(ierr);
    /* The global view, like the Ucat this pass corrects: every cell it writes is one
       this rank owns, and the ghosted image is refreshed by the caller afterwards. */
    ierr = DMDAVecGetArray(user->da, user->Friction_Velocity, &friction_velocity); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->Nu_Wall, &wall_eddy_viscosity); CHKERRQ(ierr);
    /* Read-only pressure for the wall models that need its gradient; Cabot is the only
       one that does, and it reads the previous projection's field. */
    ierr = DMDAVecGetArrayRead(user->da, user->lP, (const PetscReal ***)&wall_pressure); CHKERRQ(ierr);
    
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
    
    // Wall roughness parameter (smooth wall by default, configurable via -wall_roughness).
    const PetscReal wall_roughness_height = user->simCtx->wall_roughness_height;
    
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
                                ierr = ApplyWallModelAtCell(user, wall_pressure, wall_eddy_viscosity,
                                                    first_interior_cell, j, k,
                                                    wall_roughness_height,
                                                    distance_to_second_cell, distance_to_first_cell,
                                                    wall_velocity, reference_velocity,
                                                    &velocity_cartesian[k][j][first_interior_cell],
                                                    &friction_velocity[k][j][first_interior_cell],
                                                    wall_normal[0], wall_normal[1], wall_normal[2]); CHKERRQ(ierr);
                                
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
                                
                                ierr = ApplyWallModelAtCell(user, wall_pressure, wall_eddy_viscosity,
                                                    first_interior_cell, j, k,
                                                    wall_roughness_height,
                                                    distance_to_second_cell, distance_to_first_cell,
                                                    wall_velocity, reference_velocity,
                                                    &velocity_cartesian[k][j][first_interior_cell],
                                                    &friction_velocity[k][j][first_interior_cell],
                                                    wall_normal[0], wall_normal[1], wall_normal[2]); CHKERRQ(ierr);
                                
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
                                
                                ierr = ApplyWallModelAtCell(user, wall_pressure, wall_eddy_viscosity,
                                                    i, first_interior_cell, k,
                                                    wall_roughness_height,
                                                    distance_to_second_cell, distance_to_first_cell,
                                                    wall_velocity, reference_velocity,
                                                    &velocity_cartesian[k][first_interior_cell][i],
                                                    &friction_velocity[k][first_interior_cell][i],
                                                    wall_normal[0], wall_normal[1], wall_normal[2]); CHKERRQ(ierr);
                                
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
                                
                                ierr = ApplyWallModelAtCell(user, wall_pressure, wall_eddy_viscosity,
                                                    i, first_interior_cell, k,
                                                    wall_roughness_height,
                                                    distance_to_second_cell, distance_to_first_cell,
                                                    wall_velocity, reference_velocity,
                                                    &velocity_cartesian[k][first_interior_cell][i],
                                                    &friction_velocity[k][first_interior_cell][i],
                                                    wall_normal[0], wall_normal[1], wall_normal[2]); CHKERRQ(ierr);
                                
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
                                
                                ierr = ApplyWallModelAtCell(user, wall_pressure, wall_eddy_viscosity,
                                                    i, j, first_interior_cell,
                                                    wall_roughness_height,
                                                    distance_to_second_cell, distance_to_first_cell,
                                                    wall_velocity, reference_velocity,
                                                    &velocity_cartesian[first_interior_cell][j][i],
                                                    &friction_velocity[first_interior_cell][j][i],
                                                    wall_normal[0], wall_normal[1], wall_normal[2]); CHKERRQ(ierr);
                                
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
                                
                                ierr = ApplyWallModelAtCell(user, wall_pressure, wall_eddy_viscosity,
                                                    i, j, first_interior_cell,
                                                    wall_roughness_height,
                                                    distance_to_second_cell, distance_to_first_cell,
                                                    wall_velocity, reference_velocity,
                                                    &velocity_cartesian[first_interior_cell][j][i],
                                                    &friction_velocity[first_interior_cell][j][i],
                                                    wall_normal[0], wall_normal[1], wall_normal[2]); CHKERRQ(ierr);
                                
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
    ierr = DMDAVecRestoreArrayRead(user->da, user->lP, (const PetscReal ***)&wall_pressure); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, user->Nu_Wall, &wall_eddy_viscosity); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, user->Friction_Velocity, &friction_velocity); CHKERRQ(ierr);

    /* The correction wrote the global view; refresh the ghosted image here so that every
       caller sees a consistent field rather than each remembering to do it. */
    ierr = UpdateLocalGhosts(user, FIELD_ID_U_TAU); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_NU_WALL); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Complete.\n");
    
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "LogWallModelDiagnostics"
/**
 * @brief Implementation of \ref LogWallModelDiagnostics().
 * @details Full API contract is documented with the header declaration in
 *          `include/Boundaries.h`.
 */
PetscErrorCode LogWallModelDiagnostics(UserCtx *user)
{
    SimCtx   *simCtx = user->simCtx;
    const WallModelDiagnosticsState *state = &user->wall_diagnostics;
    MPI_Comm  comm;

    /* Index 0..3: cell count, u_tau sum, u_tau^2 sum, y+ sum. Index 4..5: wall-distance
       and wall-viscosity sums. One collective rather than six. */
    PetscReal local_sum[6], global_sum[6];
    PetscReal local_max[2], global_max[2];
    PetscReal local_min, global_min;

    PetscFunctionBeginUser;

    if (!simCtx->wallfunction) PetscFunctionReturn(0);

    PetscCall(PetscObjectGetComm((PetscObject)user->da, &comm));

    local_sum[0] = (PetscReal)state->cells;
    local_sum[1] = state->friction_velocity_sum;
    local_sum[2] = state->friction_velocity_sq;
    local_sum[3] = state->y_plus_sum;
    local_sum[4] = state->wall_distance_sum;
    local_sum[5] = state->wall_viscosity_sum;
    local_max[0] = (state->cells > 0) ? state->friction_velocity_max : 0.0;
    local_max[1] = (state->cells > 0) ? state->y_plus_max : 0.0;
    /* A rank that owns no wall face must not win the minimum with a zero it never
       measured, so it contributes the identity instead. */
    local_min    = (state->cells > 0) ? state->friction_velocity_min : PETSC_MAX_REAL;

    PetscCallMPI(MPI_Allreduce(local_sum, global_sum, 6, MPIU_REAL, MPI_SUM, comm));
    PetscCallMPI(MPI_Allreduce(local_max, global_max, 2, MPIU_REAL, MPI_MAX, comm));
    PetscCallMPI(MPI_Allreduce(&local_min, &global_min, 1, MPIU_REAL, MPI_MIN, comm));

    if (simCtx->rank == 0) {
        const PetscReal cells = global_sum[0];
        FILE           *file  = NULL;
        PetscReal       mean, mean_square, variance;

        /* No wall face anywhere is a configuration fact, not a data point: a run with
           wall functions enabled and no WALL boundary would otherwise emit a row of
           zeros that reads like a converged answer. */
        if (cells <= 0.0) {
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "Wall model is enabled but no WALL face was corrected this step; "
                      "no diagnostics written.\n");
            PetscFunctionReturn(0);
        }

        mean        = global_sum[1] / cells;
        mean_square = global_sum[2] / cells;
        variance    = PetscMax(mean_square - mean * mean, 0.0);

        PetscCall(PicurvOpenDiagnosticsCsv(simCtx, "wall_model.csv",
                                           "step,time,wall_cells,u_tau_mean,u_tau_rms,"
                                           "u_tau_min,u_tau_max,y_plus_mean,y_plus_max,"
                                           "wall_distance_mean,nu_wall_over_nu_mean", &file));
        fprintf(file,
                "%d,%.6e,%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",
                (int)simCtx->step, (double)simCtx->ti, (int)cells,
                (double)mean, (double)PetscSqrtReal(variance),
                (double)global_min, (double)global_max[0],
                (double)(global_sum[3] / cells), (double)global_max[1],
                (double)(global_sum[4] / cells),
                (double)(global_sum[5] / cells * simCtx->ren));
        PetscCheck(fclose(file) == 0, PETSC_COMM_SELF, PETSC_ERR_FILE_WRITE,
                   "Unable to close the wall-model diagnostics file.");

        LOG_ALLOW(GLOBAL, LOG_INFO,
                  "  Wall model (%s): u_tau=%.4e, y+ (mean)=%.2f, y+ (max)=%.2f over %d cell(s)\n",
                  WallFunctionModelToString((WallFunctionModel)simCtx->wallfunction),
                  (double)mean, (double)(global_sum[3] / cells),
                  (double)global_max[1], (int)cells);
    }

    /* Whether the first cell sits where the selected law is valid is a property of the
       mesh, not of the configuration, so it cannot be settled before the grid exists.
       It is checked here instead, and a run that stays outside the range is stopped
       rather than left to spend its walltime producing a wall stress the law cannot
       support. */
    {
        /* Consecutive samples tolerated outside the range. A startup transient or a
           brief excursion is not worth ending a run over; a mesh that is simply wrong
           never comes back. */
        const PetscInt  PICURV_WALL_YPLUS_GRACE_SAMPLES = 10;
        const PetscReal y_plus_mean = (global_sum[0] > 0.0) ? global_sum[3] / global_sum[0] : 0.0;
        PetscReal       lower = 0.0, upper = 300.0;
        const char     *range_reason = NULL;

        switch ((WallFunctionModel)simCtx->wallfunction) {
        case WALL_FUNCTION_LOG_LAW:
            /* Below 30 the first cell is under the logarithmic region; above 300 the
               law's own implementation reports no valid branch and the correction
               silently falls back to leaving the velocity alone. */
            lower = 30.0;
            upper = 300.0;
            range_reason = "the logarithmic region the law describes";
            break;
        case WALL_FUNCTION_WERNER:
            /* The two-layer form has a valid branch below y+ = 11.81, so only the upper
               end is a limit. */
            lower = 0.0;
            upper = 300.0;
            range_reason = "the region the power law describes";
            break;
        default:
            /* Cabot integrates across the wall layer, so it has no lower bound worth
               asserting; a first cell far outside the layer leaves its ODE nothing to
               integrate over. */
            lower = 0.0;
            upper = 1000.0;
            range_reason = "the wall layer the model integrates over";
            break;
        }

        if (global_sum[0] > 0.0 && (y_plus_mean < lower || y_plus_mean > upper)) {
            user->wall_yplus_excursions += 1;
            LOG_ALLOW(GLOBAL, LOG_WARNING,
                      "Wall model (%s): first-cell y+ is %.2f, outside %s (%.0f to %.0f). "
                      "The stress it reports is not one the law supports here. Sample %d "
                      "of %d before this run is stopped.\n",
                      WallFunctionModelToString((WallFunctionModel)simCtx->wallfunction),
                      (double)y_plus_mean, range_reason, (double)lower, (double)upper,
                      (int)user->wall_yplus_excursions,
                      (int)PICURV_WALL_YPLUS_GRACE_SAMPLES);

            PetscCheck(user->wall_yplus_excursions < PICURV_WALL_YPLUS_GRACE_SAMPLES,
                       PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
                       "Wall model (%s): first-cell y+ has been outside %s (%.0f to %.0f) "
                       "for %d consecutive samples, most recently %.2f. Refine or coarsen "
                       "the wall-normal spacing so the first cell lands in range, or "
                       "resolve the wall and disable the wall function. Stopping now "
                       "rather than spending the run producing a stress the law cannot "
                       "support.",
                       WallFunctionModelToString((WallFunctionModel)simCtx->wallfunction),
                       range_reason, (double)lower, (double)upper,
                       (int)user->wall_yplus_excursions, (double)y_plus_mean);
        } else {
            /* Consecutive, so a sample back in range clears the count. */
            user->wall_yplus_excursions = 0;
        }
    }

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FinalizePostProjectionCellFields"
/**
 * @brief Implementation of \ref FinalizePostProjectionCellFields().
 * @details Full API contract is documented with the header declaration in
 *          `include/Boundaries.h`.
 */
PetscErrorCode FinalizePostProjectionCellFields(UserCtx *user)
{
    PetscErrorCode ierr;
    const FieldId cell_fields[] = {FIELD_ID_UCAT, FIELD_ID_P};

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Finalizing post-projection cell-centered fields.\n");

    // Ensure flow-dependent Ubcs handlers see the newly reconstructed Ucat.
    ierr = UpdateLocalGhosts(user, FIELD_ID_UCAT); CHKERRQ(ierr);
    ierr = BoundarySystem_RefreshUbcs(user); CHKERRQ(ierr);

    // Establish flat non-periodic faces and periodic endpoints before corners.
    ierr = UpdateDummyCells(user); CHKERRQ(ierr);
    ierr = SynchronizePeriodicCellFields(user, 2, cell_fields); CHKERRQ(ierr);

    // Corner averaging can overwrite periodic endpoints, so restore them after.
    ierr = UpdateCornerNodes(user); CHKERRQ(ierr);
    ierr = SynchronizePeriodicCellFields(user, 2, cell_fields); CHKERRQ(ierr);

    // Synchronize explicitly because the periodic helper is a no-op when every
    // direction is non-periodic.
    ierr = UpdateLocalGhosts(user, FIELD_ID_UCAT); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_P); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Post-projection cell-centered fields finalized.\n");
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ApplyBoundaryConditions"
/**
 * @brief Implementation of \ref ApplyBoundaryConditions().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Boundaries.h`.
 * @see ApplyBoundaryConditions()
 */
PetscErrorCode ApplyBoundaryConditions(UserCtx *user)
{
    PetscErrorCode ierr;
    const FieldId staggered_fields[] = {FIELD_ID_UCONT};
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
        ierr = SynchronizePeriodicStaggeredFields(user, 1, staggered_fields); CHKERRQ(ierr);

        // (c) Convert updated Contravariant velocities to Cartesian velocities.
        ierr = Contra2Cart(user); CHKERRQ(ierr);

        // (d) Synchronize the updated Cartesian velocities across all processors
        //     to ensure all ucat values are current before updating the dummy cells.
        ierr = UpdateLocalGhosts(user, FIELD_ID_UCAT); CHKERRQ(ierr);

        // (e) If Wall functions are enabled, apply them now to adjust near-wall velocities.
        if(user->simCtx->wallfunction){
          // Apply wall function adjustments to the boundary velocities.
          ierr = ApplyWallFunction(user); CHKERRQ(ierr);

          // Synchronize the updated Cartesian velocities after wall function adjustments.
          ierr = UpdateLocalGhosts(user, FIELD_ID_UCAT); CHKERRQ(ierr);

          LOG_ALLOW(GLOBAL,LOG_VERBOSE,"Wall Function Applied at Walls.\n");
        }

        // (f) Update the first layer of ghost cells for non-periodic faces using
        //     the newly computed `ubcs` values.
        ierr = UpdateDummyCells(user); CHKERRQ(ierr);

        LOG_ALLOW(GLOBAL,LOG_VERBOSE,"Dummy Cells/Ghost Cells Updated.\n");

        // (g) Handle all periodic boundaries. This is a parallel direct copy
        // that sets the absolute constraints for the rest of the solve.
        // There is a Ghost update happening inside this function.
        ierr = ApplyPeriodicBCs(user); CHKERRQ(ierr);

        // (h) Update the corner and edge ghost nodes. This routine calculates
        // values for corners/edges by averaging their neighbors, which have been
        // finalized in the steps above (both periodic and non-periodic).
        ierr = UpdateCornerNodes(user); CHKERRQ(ierr);
        
        // (i) Synchronize the updated edge and corner cells across all processors to ensure
        //     consistency before the next iteration or finalization.
        ierr = UpdateLocalGhosts(user, FIELD_ID_P); CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, FIELD_ID_UCAT); CHKERRQ(ierr);
        ierr = SynchronizePeriodicStaggeredFields(user, 1, staggered_fields); CHKERRQ(ierr);

        // (j) Ensure All the corners are synchronized with  a well defined protocol  in case of Periodic boundary conditions
        // To avoid race conditions.
        const FieldId all_fields[] = {FIELD_ID_UCAT, FIELD_ID_P, FIELD_ID_NVERT};
        ierr = SynchronizePeriodicCellFields(user, 3, all_fields); CHKERRQ(ierr);

    }

    // STEP 3: Final ghost node synchronization. This ensures all changes made
    // to the global vectors are reflected in the local ghost regions of all
    // processors, making the state fully consistent before the next solver stage.
    ierr = UpdateLocalGhosts(user, FIELD_ID_P); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_UCAT); CHKERRQ(ierr);
    ierr = SynchronizePeriodicStaggeredFields(user, 1, staggered_fields); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
