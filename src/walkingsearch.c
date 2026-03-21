// walkingsearch.c

#include "walkingsearch.h"
#include "logging.h"
#include <petsc.h>
#include <stdbool.h>
#include <math.h>

// Define maximum traversal steps to prevent infinite loops
#define MAX_TRAVERSAL 1000
#define DISTANCE_THRESHOLD 1e-11
#define REPEAT_COUNT_THRESHOLD 5

#undef __FUNCT__
#define __FUNCT__ "GetCellCharacteristicSize"
/**
 * @brief Implementation of \ref GetCellCharacteristicSize().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/walkingsearch.h`.
 * @see GetCellCharacteristicSize()
 */
PetscErrorCode GetCellCharacteristicSize(const Cell *cell, PetscReal *cellSize)
{
    PetscErrorCode ierr = 0;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    if (!cell || !cellSize) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, 
        "GetCellCharacteristicSize: Null pointer(s).");

    // A simple approach: compute the average of the distances between
    // all pairs of adjacent vertices that share an edge. That gives
    // a typical measure of cell dimension. For a uniform grid, you
    // could do something simpler.

    PetscInt edges[12][2] = {
        {0,1},{1,2},{2,3},{3,0},  // bottom face edges
        {4,5},{5,6},{6,7},{7,4},  // top face edges
        {0,7},{1,6},{2,5},{3,4}   // vertical edges
    };

    PetscReal totalEdgeLen = 0.0;
    for (PetscInt i=0; i<12; i++) {
        PetscInt vA = edges[i][0];
        PetscInt vB = edges[i][1];
        Cmpnts A = cell->vertices[vA];
        Cmpnts B = cell->vertices[vB];
        PetscReal dx = B.x - A.x;
        PetscReal dy = B.y - A.y;
        PetscReal dz = B.z - A.z;
        PetscReal edgeLen = sqrt(dx*dx + dy*dy + dz*dz);
        totalEdgeLen += edgeLen;
    }

    // Average edge length
    *cellSize = totalEdgeLen / 12.0;

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(ierr);
}


#undef __FUNCT__
#define __FUNCT__ "ComputeSignedDistanceToPlane"
/**
 * @brief Internal helper implementation: `ComputeSignedDistanceToPlane()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ComputeSignedDistanceToPlane(const Cmpnts v1, const Cmpnts v2, const Cmpnts v3, const Cmpnts v4,
                                            const Cmpnts cell_centroid,const Cmpnts p_target,
					    PetscReal *d_signed, const PetscReal threshold)
{
    PetscErrorCode ierr = 0;
    PetscMPIInt rank;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    if (d_signed == NULL) {
        ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank); CHKERRQ(ierr);
        LOG_ALLOW(LOCAL, LOG_ERROR, "Output pointer 'd_signed' is NULL on rank %d.\n", rank);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Output pointer 'd_signed' must not be NULL.");
    }

    LOG_ALLOW(LOCAL, LOG_VERBOSE, "Target point: (%.6e, %.6e, %.6e)\n", p_target.x, p_target.y, p_target.z);
    LOG_ALLOW(LOCAL, LOG_VERBOSE, "  Cell Centroid: (%.6e, %.6e, %.6e)\n", cell_centroid.x, cell_centroid.y, cell_centroid.z);
    LOG_ALLOW(LOCAL, LOG_VERBOSE, "  Face Vertices:\n");
    LOG_ALLOW(LOCAL, LOG_VERBOSE, "    v1: (%.6e, %.6e, %.6e)\n", v1.x, v1.y, v1.z);
    LOG_ALLOW(LOCAL, LOG_VERBOSE, "    v2: (%.6e, %.6e, %.6e)\n", v2.x, v2.y, v2.z);
    LOG_ALLOW(LOCAL, LOG_VERBOSE, "    v3: (%.6e, %.6e, %.6e)\n", v3.x, v3.y, v3.z);
    LOG_ALLOW(LOCAL, LOG_VERBOSE, "    v4: (%.6e, %.6e, %.6e)\n", v4.x, v4.y, v4.z);

    // --- Calculate Edge Vectors for Initial Normal Computation ---
    PetscReal edge1_x = v2.x - v1.x;
    PetscReal edge1_y = v2.y - v1.y;
    PetscReal edge1_z = v2.z - v1.z;

    PetscReal edge2_x = v4.x - v1.x;
    PetscReal edge2_y = v4.y - v1.y;
    PetscReal edge2_z = v4.z - v1.z;
    LOG_ALLOW(LOCAL, LOG_VERBOSE, "  Edge1 (v2-v1): (%.6e, %.6e, %.6e)\n", edge1_x, edge1_y, edge1_z);
    LOG_ALLOW(LOCAL, LOG_VERBOSE, "  Edge2 (v4-v1): (%.6e, %.6e, %.6e)\n", edge2_x, edge2_y, edge2_z);

    // --- Compute Initial Normal Vector (Cross Product: edge1 x edge2) ---
    PetscReal normal_x_initial = edge1_y * edge2_z - edge1_z * edge2_y;  
    PetscReal normal_y_initial = edge1_z * edge2_x - edge1_x * edge2_z;
    PetscReal normal_z_initial = edge1_x * edge2_y - edge1_y * edge2_x;
    LOG_ALLOW(LOCAL, LOG_VERBOSE, "  Initial Raw Normal (edge1 x edge2): (%.6e, %.6e, %.6e)\n", normal_x_initial, normal_y_initial, normal_z_initial);

    PetscReal normal_magnitude = sqrt(normal_x_initial * normal_x_initial +
                                   normal_y_initial * normal_y_initial +
                                   normal_z_initial * normal_z_initial);
    LOG_ALLOW(LOCAL, LOG_VERBOSE, "  Initial Normal Magnitude: %.6e\n", normal_magnitude);

    if (normal_magnitude < 1.0e-12) {
        ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank); CHKERRQ(ierr);
        LOG_ALLOW(LOCAL, LOG_ERROR,
                  "Degenerate plane detected on rank %d. Normal magnitude (%.3e) is too small.\n",
                  rank, normal_magnitude);
        LOG_ALLOW(LOCAL, LOG_ERROR, "  Offending vertices for normal: v1(%.3e,%.3e,%.3e), v2(%.3e,%.3e,%.3e), v4(%.3e,%.3e,%.3e)\n",
                  v1.x,v1.y,v1.z, v2.x,v2.y,v2.z, v4.x,v4.y,v4.z);
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_USER, "Degenerate plane detected (normal vector is near zero).");
    }

    // --- Calculate the Centroid of the Four Face Vertices (v1, v2, v3, v4) ---    
    PetscReal face_centroid_x = 0.25 * (v1.x + v2.x + v3.x + v4.x);
    PetscReal face_centroid_y = 0.25 * (v1.y + v2.y + v3.y + v4.y);
    PetscReal face_centroid_z = 0.25 * (v1.z + v2.z + v3.z + v4.z);
    LOG_ALLOW(LOCAL, LOG_VERBOSE, "  Face Centroid: (%.6e, %.6e, %.6e)\n", face_centroid_x, face_centroid_y, face_centroid_z);

    // --- Orient the Normal to Point Outward from the Cell ---
    // Vector from face centroid to cell centroid
    PetscReal vec_fc_to_cc_x = cell_centroid.x - face_centroid_x;
    PetscReal vec_fc_to_cc_y = cell_centroid.y - face_centroid_y;
    PetscReal vec_fc_to_cc_z = cell_centroid.z - face_centroid_z;
    LOG_ALLOW(LOCAL, LOG_VERBOSE, "  Vec (FaceCentroid -> CellCentroid): (%.6e, %.6e, %.6e)\n", vec_fc_to_cc_x, vec_fc_to_cc_y, vec_fc_to_cc_z);

    // Dot product of initial normal with vector from face centroid to cell centroid
    PetscReal dot_prod_orientation = normal_x_initial * vec_fc_to_cc_x +
                                     normal_y_initial * vec_fc_to_cc_y +
                                     normal_z_initial * vec_fc_to_cc_z;
    LOG_ALLOW(LOCAL, LOG_VERBOSE, "  Dot Product for Orientation (N_initial . Vec_FC_to_CC): %.6e\n", dot_prod_orientation);

    PetscReal normal_x = normal_x_initial;
    PetscReal normal_y = normal_y_initial;
    PetscReal normal_z = normal_z_initial;

    // If dot_prod_orientation > 0, initial normal points towards cell interior, so flip it.
    if (dot_prod_orientation > 1.0e-9) { // Use a small epsilon to avoid issues if dot product is extremely close to zero
        normal_x = -normal_x_initial;
        normal_y = -normal_y_initial;
        normal_z = -normal_z_initial;
        LOG_ALLOW(LOCAL, LOG_VERBOSE, "  Initial normal was inward (dot_prod > 0). Flipped normal.\n");
    } else if (dot_prod_orientation == 0.0 && normal_magnitude > 1e-12) {
        // This case is ambiguous or face plane contains cell centroid.
        // This might happen for highly symmetric cells or if face_centroid IS cell_centroid (e.g. 2D cell).
        // For now, we keep the original normal direction based on v1,v2,v4 ordering.
        // A more robust solution for this edge case might be needed if it occurs often.
        ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank); CHKERRQ(ierr);
         LOG_ALLOW(LOCAL, LOG_WARNING, "Rank %d: Dot product for normal orientation is zero. Normal direction might be ambiguous. Keeping initial normal direction from (v2-v1)x(v4-v1).\n", rank);
    }
    LOG_ALLOW(LOCAL, LOG_VERBOSE, "  Oriented Raw Normal: (%.6e, %.6e, %.6e)\n", normal_x, normal_y, normal_z);


    // --- Normalize the (Now Outward-Pointing) Normal Vector ---
    // Note: normal_magnitude was calculated from initial normal.
    // If we flipped, magnitude is the same.
    normal_x /= normal_magnitude;
    normal_y /= normal_magnitude;
    normal_z /= normal_magnitude;
    LOG_ALLOW(LOCAL, LOG_VERBOSE, "  Normalized Outward Normal: (%.6f, %.6f, %.6f)\n", normal_x, normal_y, normal_z);


    // --- Compute Vector from Target Point to Face Centroid ---
    PetscReal vec_p_to_fc_x = face_centroid_x - p_target.x;
    PetscReal vec_p_to_fc_y = face_centroid_y - p_target.y;
    PetscReal vec_p_to_fc_z = face_centroid_z - p_target.z;

    // --- Compute Signed Distance ---
    *d_signed = vec_p_to_fc_x * normal_x +
                vec_p_to_fc_y * normal_y +
                vec_p_to_fc_z * normal_z;

    LOG_ALLOW(LOCAL, LOG_VERBOSE, "  Raw Signed Distance (using outward normal): %.15e\n", *d_signed);

    // --- Apply Threshold ---
    if (fabs(*d_signed) < threshold) {
        LOG_ALLOW(LOCAL, LOG_VERBOSE, "  Distance %.15e is less than threshold %.1e. Setting to 0.0.\n", *d_signed, threshold);
        *d_signed = 0.0;
    }

    LOG_ALLOW(LOCAL, LOG_VERBOSE, "Final Signed Distance: %.15e\n", *d_signed);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "CalculateDistancesToCellFaces"
/**
 * @brief Implementation of \ref CalculateDistancesToCellFaces().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/walkingsearch.h`.
 * @see CalculateDistancesToCellFaces()
 */
PetscErrorCode CalculateDistancesToCellFaces(const Cmpnts p, const Cell *cell, PetscReal *d, const PetscReal threshold)
{
    PetscErrorCode ierr;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // Validate that the 'cell' pointer is not NULL to prevent dereferencing a null pointer.
    if (cell == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "'cell' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "'cell' is NULL.");
    }

    // Validate that the 'd' pointer is not NULL to ensure there is memory allocated for distance storage.
    if (d == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "'d' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "'d' is NULL.");
    }

    // Compute the centroid of the entire cell
    Cmpnts cell_centroid = {0.0, 0.0, 0.0};
    for (int i = 0; i < 8; ++i) {
      cell_centroid.x += cell->vertices[i].x;
      cell_centroid.y += cell->vertices[i].y;
      cell_centroid.z += cell->vertices[i].z;
    }
    cell_centroid.x /= 8.0;
    cell_centroid.y /= 8.0;
    cell_centroid.z /= 8.0;

    LOG_ALLOW(LOCAL,LOG_DEBUG, "Cell Centroid: (%.6e, %.6e, %.6e)\n",
	      cell_centroid.x, cell_centroid.y, cell_centroid.z);    

    
    // Compute the signed distance from point 'p' to the BACK face of the cell.
    // The BACK face is defined by vertices 0, 3, 2, and 1, with its normal vector pointing in the -z direction.
    ierr = ComputeSignedDistanceToPlane(
        cell->vertices[0], // Vertex 0
        cell->vertices[3], // Vertex 3
        cell->vertices[2], // Vertex 2
        cell->vertices[1], // Vertex 1
	cell_centroid,     // Cell centroid 
        p,                  // Target point
        &d[BACK],           // Storage location for the BACK face distance
        threshold           // Threshold for zero distance
    );  CHKERRQ(ierr);

    // Compute the signed distance from point 'p' to the FRONT face of the cell.
    // The FRONT face is defined by vertices 4, 7, 6, and 5, with its normal vector pointing in the +z direction.
    ierr = ComputeSignedDistanceToPlane(
        cell->vertices[4], // Vertex 4
        cell->vertices[7], // Vertex 7
        cell->vertices[6], // Vertex 6
        cell->vertices[5], // Vertex 5
	cell_centroid,     // Cell centroid 
        p,                  // Target point
        &d[FRONT],          // Storage location for the FRONT face distance
        threshold           // Threshold for zero distance
    );  CHKERRQ(ierr);

    // Compute the signed distance from point 'p' to the BOTTOM face of the cell.
    // The BOTTOM face is defined by vertices 0, 1, 6, and 7, with its normal vector pointing in the -y direction.
    ierr = ComputeSignedDistanceToPlane(
        cell->vertices[0], // Vertex 0
        cell->vertices[1], // Vertex 1
        cell->vertices[6], // Vertex 6
        cell->vertices[7], // Vertex 7
	cell_centroid,     // Cell centroid 
        p,                  // Target point
        &d[BOTTOM],         // Storage location for the BOTTOM face distance
        threshold           // Threshold for zero distance
    );  CHKERRQ(ierr);

    // Compute the signed distance from point 'p' to the TOP face of the cell.
    // The TOP face is defined by vertices 3, 4, 5, and 2, with its normal vector pointing in the +y direction.
    ierr = ComputeSignedDistanceToPlane(
        cell->vertices[3], // Vertex 3
        cell->vertices[4], // Vertex 4
        cell->vertices[5], // Vertex 5
        cell->vertices[2], // Vertex 2
	cell_centroid,     // Cell centroid 
	p,                  // Target point
        &d[TOP],            // Storage location for the TOP face distance
        threshold           // Threshold for zero distance
    );  CHKERRQ(ierr);

    // Compute the signed distance from point 'p' to the LEFT face of the cell.
    // The LEFT face is defined by vertices 0, 7, 4, and 3, with its normal vector pointing in the -x direction.
    ierr = ComputeSignedDistanceToPlane(
        cell->vertices[0], // Vertex 0
        cell->vertices[7], // Vertex 7
        cell->vertices[4], // Vertex 4
        cell->vertices[3], // Vertex 3
	cell_centroid,     // Cell centroid 
	p,                  // Target point
        &d[LEFT],           // Storage location for the LEFT face distance
        threshold           // Threshold for zero distance
    );  CHKERRQ(ierr);

    // Compute the signed distance from point 'p' to the RIGHT face of the cell.
    // The RIGHT face is defined by vertices 1, 2, 5, and 6, with its normal vector pointing in the +x direction.
    ierr = ComputeSignedDistanceToPlane(
        cell->vertices[1], // Vertex 1
        cell->vertices[2], // Vertex 2
        cell->vertices[5], // Vertex 5
        cell->vertices[6], // Vertex 6
	cell_centroid,     // Cell centroid 
        p,                  // Target point
        &d[RIGHT],          // Storage location for the RIGHT face distance
        threshold           // Threshold for zero distance
    );  CHKERRQ(ierr);

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Populated d: "
	      "d[LEFT=%d]=%.3e, d[RIGHT=%d]=%.3e, d[BOTTOM=%d]=%.3e, "
	      "d[TOP=%d]=%.3e, d[FRONT=%d]=%.3e, d[BACK=%d]=%.3e\n",
	      LEFT, d[LEFT], RIGHT, d[RIGHT], BOTTOM, d[BOTTOM],
	      TOP, d[TOP], FRONT, d[FRONT], BACK, d[BACK]);
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Raw d: "
	      "[%.3e, %.3e, %.3e, %.3e, %.3e, %.3e]\n",
	      d[0], d[1], d[2], d[3], d[4], d[5]);
    
    PROFILE_FUNCTION_END;      
    PetscFunctionReturn(0); // Indicate successful execution of the function.
}

#undef __FUNCT__
#define __FUNCT__ "DeterminePointPosition"
/**
 * @brief Implementation of \ref DeterminePointPosition().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/walkingsearch.h`.
 * @see DeterminePointPosition()
 */
PetscErrorCode DeterminePointPosition(PetscReal *d, PetscInt *result)
{
    
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    // Validate input pointers
    if (d == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "'d' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input parameter 'd' is NULL.");
    }
    if (result == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "'result' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Output parameter 'result' is NULL.");
    }

    // Initialize flags
    PetscBool isInside = PETSC_TRUE;
    PetscBool isOnBoundary = PETSC_FALSE;
    PetscInt IntersectionCount = 0; // Counts the number of intersections of the point with various planes of the cell.

    // Analyze distances to determine position
    for(int i = 0; i < NUM_FACES; i++) {
        if(d[i] <  0.0) {
            isInside = PETSC_FALSE; // Point is outside in at least one direction
        }
        if(d[i] == 0.0) {
            isOnBoundary = PETSC_TRUE; // Point is exactly at least one face
            IntersectionCount++; 
        }
    }

    // Set the result based on flags
    if(isInside && isOnBoundary) {
      if(IntersectionCount == 1){ 
          *result = 1; // on face
          LOG_ALLOW(LOCAL,LOG_VERBOSE, "Particle is on a face of the cell.\n");
      }
      else if(IntersectionCount == 2){ 
          *result = 2; // on edge
          LOG_ALLOW(LOCAL,LOG_VERBOSE, "Particle is on an edge of the cell.\n");
      }
      else if(IntersectionCount >= 3){ 
          *result = 3; // on corner
          LOG_ALLOW(LOCAL,LOG_VERBOSE, "Particle is on a corner of the cell.\n");
      }
    }
    else if(isInside) {
        *result = 0; // Inside the cell
        LOG_ALLOW(LOCAL,LOG_VERBOSE, "Particle is inside the cell.\n");
    }
    else {
        *result = -1; // Outside the cell
        LOG_ALLOW(LOCAL,LOG_VERBOSE, "Particle is outside the cell.\n");
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0); // Indicate successful execution
}

#undef __FUNCT__
#define __FUNCT__ "GetCellVerticesFromGrid"
/**
 * @brief Implementation of \ref GetCellVerticesFromGrid().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/walkingsearch.h`.
 * @see GetCellVerticesFromGrid()
 */
PetscErrorCode GetCellVerticesFromGrid(Cmpnts ***coor, PetscInt idx, PetscInt idy, PetscInt idz,
                                       Cell *cell)
{

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    // Validate input pointers
    if (coor == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "'coor' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input array 'coor' is NULL.");
    }
    if (cell == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "'cell' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Output parameter 'cell' is NULL.");
    }

    // Assign vertices based on grid indices matching the figure
    // Vertex numbering follows the convention depicted in the figure
    cell->vertices[0] = coor[idz][idy][idx];         // Vertex 0: (i,   j,   k)
    cell->vertices[1] = coor[idz][idy][idx+1];       // Vertex 1: (i+1, j,   k)
    cell->vertices[2] = coor[idz][idy+1][idx+1];     // Vertex 2: (i+1, j+1, k)
    cell->vertices[3] = coor[idz][idy+1][idx];       // Vertex 3: (i,   j+1, k)
    cell->vertices[4] = coor[idz+1][idy+1][idx];     // Vertex 4: (i,   j+1, k+1)
    cell->vertices[5] = coor[idz+1][idy+1][idx+1];   // Vertex 5: (i+1, j+1, k+1)
    cell->vertices[6] = coor[idz+1][idy][idx+1];     // Vertex 6: (i+1, j,   k+1)
    cell->vertices[7] = coor[idz+1][idy][idx];       // Vertex 7: (i,   j,   k+1)

    LOG_ALLOW(LOCAL,LOG_VERBOSE, "Retrieved vertices for cell (%d, %d, %d).\n", idx, idy, idz);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0); // Indicate successful execution
}

#undef __FUNCT__
#define __FUNCT__ "InitializeTraversalParameters"
/**
 * @brief Implementation of \ref InitializeTraversalParameters().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/walkingsearch.h`.
 * @see InitializeTraversalParameters()
 */
PetscErrorCode InitializeTraversalParameters(UserCtx *user, Particle *particle, PetscInt *idx, PetscInt *idy, PetscInt *idz, PetscInt *traversal_steps)
{
    PetscErrorCode ierr;
    DMDALocalInfo info;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // Validate input pointers
    if (user == NULL || particle == NULL || idx == NULL || idy == NULL || idz == NULL || traversal_steps == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "One or more input pointers are NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "One or more input pointers are NULL.");
    }

    // Get grid information (needed for fallback start and potentially validation)
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);

    // --- Check if the particle has a valid previous cell ID ---
    // Assuming particle->cell stores the *global* indices
    if (particle->cell[0] >= 0 && particle->cell[1] >= 0 && particle->cell[2] >= 0) {

      // It sees a valid cell ID exists.
      // Before using it, it explicitly checks if this cell is accessible.
      PetscBool is_handoff_cell_valid;
      ierr = CheckCellWithinLocalGrid(user, particle->cell[0], particle->cell[1], particle->cell[2], &is_handoff_cell_valid); CHKERRQ(ierr);

      if (is_handoff_cell_valid) {
        // FAST PATH: The check passed. The cell is safe. Use it.
        *idx = particle->cell[0];
        *idy = particle->cell[1];
        *idz = particle->cell[2];

	 LOG_ALLOW(LOCAL,LOG_DEBUG, "Particle %lld has previous cell (%d, %d, %d). Starting search there.\n",
		   (long long)particle->PID, *idx, *idy, *idz); // Cast id to long long for printing PetscInt64
	
      } else {
        // SAFE FALLBACK: The check failed! The handoff cell is NOT safe.
        // Discard the unsafe handoff cell and revert to starting at the
        // guaranteed-safe local corner.
        *idx = info.xs;
        *idy = info.ys;
        *idz = info.zs;
      }
    } else {
        // No valid previous cell ID (e.g., -1,-1,-1 or first time).
        // IMPROVED: Use distance-based multi-seed approach to find the closest strategic cell.
        // Sample the 6 face centers + volume center (7 total candidates) and start from the closest.

        Cmpnts ***cent;

        // Get local cell center array (lCent)
        // For cell (i,j,k), the center is at cent[k+1][j+1][i+1]
        ierr = DMDAVecGetArrayRead(user->fda, user->lCent, &cent); CHKERRQ(ierr);

        // Define 7 strategic sampling points: 6 face centers + 1 volume center
        PetscInt candidates[7][3];
        const char* candidate_names[7] = {"X-min face", "X-max face", "Y-min face",
                                          "Y-max face", "Z-min face", "Z-max face", "Volume center"};

        // Calculate middle indices for each dimension
        PetscInt mid_i = info.xs + (info.xm - 1) / 2;
        PetscInt mid_j = info.ys + (info.ym - 1) / 2;
        PetscInt mid_k = info.zs + (info.zm - 1) / 2;
        PetscInt max_i = info.xs + info.xm - 2;
        PetscInt max_j = info.ys + info.ym - 2;
        PetscInt max_k = info.zs + info.zm - 2;

        // Clamp max indices to ensure valid cells
        if (max_i < info.xs) max_i = info.xs;
        if (max_j < info.ys) max_j = info.ys;
        if (max_k < info.zs) max_k = info.zs;

        // X-min face center (i=min, j=mid, k=mid)
        candidates[0][0] = info.xs; candidates[0][1] = mid_j; candidates[0][2] = mid_k;
        // X-max face center (i=max, j=mid, k=mid)
        candidates[1][0] = max_i; candidates[1][1] = mid_j; candidates[1][2] = mid_k;
        // Y-min face center (i=mid, j=min, k=mid)
        candidates[2][0] = mid_i; candidates[2][1] = info.ys; candidates[2][2] = mid_k;
        // Y-max face center (i=mid, j=max, k=mid)
        candidates[3][0] = mid_i; candidates[3][1] = max_j; candidates[3][2] = mid_k;
        // Z-min face center (i=mid, j=mid, k=min)
        candidates[4][0] = mid_i; candidates[4][1] = mid_j; candidates[4][2] = info.zs;
        // Z-max face center (i=mid, j=mid, k=max)
        candidates[5][0] = mid_i; candidates[5][1] = mid_j; candidates[5][2] = max_k;
        // Volume center (i=mid, j=mid, k=mid)
        candidates[6][0] = mid_i; candidates[6][1] = mid_j; candidates[6][2] = mid_k;

        // Find the closest candidate cell
        PetscReal min_distance = PETSC_MAX_REAL;
        PetscInt best_candidate = 6; // Default to volume center

        for (PetscInt i = 0; i < 7; i++) {
            PetscInt ci = candidates[i][0];
            PetscInt cj = candidates[i][1];
            PetscInt ck = candidates[i][2];

            // Get cell center coordinates: cent[k+1][j+1][i+1] for cell (i,j,k)
            Cmpnts cell_center = cent[ck+1][cj+1][ci+1];

            // Calculate Euclidean distance from particle to cell center
            PetscReal dx = particle->loc.x - cell_center.x;
            PetscReal dy = particle->loc.y - cell_center.y;
            PetscReal dz = particle->loc.z - cell_center.z;
            PetscReal distance = sqrt(dx*dx + dy*dy + dz*dz);

            LOG_ALLOW(LOCAL, LOG_DEBUG, "  Candidate %d (%s) at cell (%d,%d,%d): center=(%.3e,%.3e,%.3e), distance=%.3e\n",
                      i, candidate_names[i], ci, cj, ck, cell_center.x, cell_center.y, cell_center.z, distance);

            if (distance < min_distance) {
                min_distance = distance;
                best_candidate = i;
                *idx = ci;
                *idy = cj;
                *idz = ck;
            }
        }

        // Restore array
        ierr = DMDAVecRestoreArrayRead(user->fda, user->lCent, &cent); CHKERRQ(ierr);

        LOG_ALLOW(LOCAL, LOG_INFO, "Particle %lld has no valid previous cell. Multi-seed search selected %s at cell (%d,%d,%d) with distance %.3e.\n",
                  (long long)particle->PID, candidate_names[best_candidate], *idx, *idy, *idz, min_distance);
    }

    // Initialize traversal step counter
    *traversal_steps = 0;

    // Log the chosen starting point
    LOG_ALLOW(LOCAL,LOG_INFO, "Traversal for particle %lld initialized to start at cell (%d, %d, %d).\n",
               (long long)particle->PID, *idx, *idy, *idz);

    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "CheckCellWithinLocalGrid"
/**
 * @brief Implementation of \ref CheckCellWithinLocalGrid().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/walkingsearch.h`.
 * @see CheckCellWithinLocalGrid()
 */
PetscErrorCode CheckCellWithinLocalGrid(UserCtx *user, PetscInt idx, PetscInt idy, PetscInt idz, PetscBool *is_within)
{
    PetscErrorCode ierr;
    DMDALocalInfo info_nodes; // Node information from the DMDA that defines ghost regions (user->fda)

    PetscFunctionBeginUser; // Assuming this is part of your PETSc style
    PROFILE_FUNCTION_BEGIN;

    // Validate inputs
    if (user == NULL || is_within == NULL) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "Input pointer is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input pointer is NULL in CheckCellWithinLocalGrid.");
    }
    if (user->fda == NULL) {
        LOG_ALLOW(LOCAL, LOG_ERROR, "user->fda is NULL.Cannot get ghost info.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "user->fda is NULL. Cannot get ghost info.");
    }

    // Get node info from user->fda (this DMDA has the ghost layer information for nodes)
    ierr = DMDAGetLocalInfo(user->fda, &info_nodes); CHKERRQ(ierr);

    // Determine the range of GLOBAL CELL INDICES that are covered by this rank's ghosted NODAL region.
    // A cell C(i,j,k) has origin node N(i,j,k).
    // The ghosted nodal region starts at global node index info_nodes.gxs and has info_nodes.gxm nodes.

    // Global starting index of the first cell whose origin node is within the ghosted nodal region.
    PetscInt gxs_cell_global_start = info_nodes.gxs;
    PetscInt gys_cell_global_start = info_nodes.gys;
    PetscInt gzs_cell_global_start = info_nodes.gzs;

    // Number of cells that can be formed starting from nodes within the ghosted nodal region.
    // If there are N ghosted nodes (info_nodes.gxm), they can be origins for N-1 cells.
    PetscInt gxm_cell_local_count = (info_nodes.gxm > 0) ? info_nodes.gxm - 1 : 0;
    PetscInt gym_cell_local_count = (info_nodes.gym > 0) ? info_nodes.gym - 1 : 0;
    PetscInt gzm_cell_local_count = (info_nodes.gzm > 0) ? info_nodes.gzm - 1 : 0;

    // Global exclusive end index for cells whose origins are in the ghosted nodal region. 
    PetscInt gxe_cell_global_end_exclusive = gxs_cell_global_start + gxm_cell_local_count;
    PetscInt gye_cell_global_end_exclusive = gys_cell_global_start + gym_cell_local_count;
    PetscInt gze_cell_global_end_exclusive = gzs_cell_global_start + gzm_cell_local_count;

    // Check if the given global cell index (idx, idy, idz) falls within this range.
    // This means the origin node of cell (idx,idy,idz) is within the rank's accessible ghosted node region,
    // and that node can indeed serve as a cell origin (i.e., it's not the very last node in the ghosted region).
    if (idx >= gxs_cell_global_start && idx < gxe_cell_global_end_exclusive &&
        idy >= gys_cell_global_start && idy < gye_cell_global_end_exclusive &&
        idz >= gzs_cell_global_start && idz < gze_cell_global_end_exclusive) {
        *is_within = PETSC_TRUE;
    } else {
        *is_within = PETSC_FALSE;
    }

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Cell (origin node glob idx) (%d, %d, %d) is %s the ghosted local grid (covered cell origins x:[%d..%d), y:[%d..%d), z:[%d..%d)).\n",
              idx, idy, idz, (*is_within) ? "within" : "outside",
              gxs_cell_global_start, gxe_cell_global_end_exclusive,
              gys_cell_global_start, gye_cell_global_end_exclusive,
              gzs_cell_global_start, gze_cell_global_end_exclusive);

    PROFILE_FUNCTION_END;          
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "RetrieveCurrentCell"
/**
 * @brief Implementation of \ref RetrieveCurrentCell().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/walkingsearch.h`.
 * @see RetrieveCurrentCell()
 */
PetscErrorCode RetrieveCurrentCell(UserCtx *user, PetscInt idx, PetscInt idy, PetscInt idz, Cell *cell)
{
    PetscErrorCode ierr;
    Vec Coor;
    Cmpnts ***coor_array;
    PetscMPIInt rank;
    DMDALocalInfo info_nodes;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // Validate input pointers
    if (user == NULL || cell == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "One or more input pointers are NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "One or more input pointers are NULL.");
    }

    // Get local coordinates
    ierr = DMGetCoordinatesLocal(user->da, &Coor); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, Coor, &coor_array); CHKERRQ(ierr);

    // Get the local grid information FOR THE GHOSTED ARRAY from user->fda
    // This info contains the mapping between global and local indices.
    ierr = DMDAGetLocalInfo(user->fda, &info_nodes); CHKERRQ(ierr);

    PetscInt  idx_local = idx; // - info_nodes.gxs;
    PetscInt  idy_local = idy; // - info_nodes.gys;
    PetscInt  idz_local = idz; // - info_nodes.gzs;
    
    LOG_ALLOW(LOCAL,LOG_VERBOSE," [Rank %d] Getting vertex coordinates for cell: %d,%d,%d whose local coordinates are %d,%d,%d.\n",rank,idx,idy,idz,idx_local,idy_local,idz_local);
    
    // Get the current cell's vertices
    ierr = GetCellVerticesFromGrid(coor_array, idx_local, idy_local, idz_local, cell); CHKERRQ(ierr);

    // Restore array
    ierr = DMDAVecRestoreArrayRead(user->fda, Coor, &coor_array); CHKERRQ(ierr);

    // Get MPI rank for debugging
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

    // Debug: Print cell vertices
    LOG_ALLOW(LOCAL,LOG_VERBOSE, "Cell (%d, %d, %d) vertices \n", idx, idy, idz);
    ierr = LOG_CELL_VERTICES(cell, rank); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EvaluateParticlePosition"
/**
 * @brief Implementation of \ref EvaluateParticlePosition().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/walkingsearch.h`.
 * @see EvaluateParticlePosition()
 */
PetscErrorCode EvaluateParticlePosition(const Cell *cell, PetscReal *d, const Cmpnts p, PetscInt *position, const PetscReal threshold)
{
    PetscErrorCode ierr;
    PetscReal cellSize;
    PetscReal cellThreshold;
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // Validate input pointers to ensure they are not NULL, preventing potential segmentation faults.
    if (cell == NULL || position == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "One or more input pointers are NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "One or more input pointers are NULL.");
    }
  
    // Compute a local cell size
    ierr = GetCellCharacteristicSize(cell, &cellSize); CHKERRQ(ierr);

    // scale base threshold by cell size (ex: if threshold = 1e-6 and cell size = 0.01, cellThreshold = 1e-8)
    cellThreshold = threshold*cellSize; 

    // Invoke the function to calculate signed distances from the particle to each face of the cell.
    // The distances are stored in the array pointed to by 'd'.
    ierr = CalculateDistancesToCellFaces(p, cell, d, cellThreshold); CHKERRQ(ierr);
    CHKERRQ(ierr); // Check for errors in distance calculation.


    // Catch degenerate-plane error manually:
    if (ierr == PETSC_ERR_USER) {
        LOG_ALLOW(LOCAL, LOG_WARNING,
                  "Skipping cell due to degenerate face.\n");
        // We can set *position = -1 here
        *position = -1; // treat as outside
        return 0; // not a fatal error, just skip
    } else {
        CHKERRQ(ierr);
    }


    // Debugging output: Print the computed distances to each face for verification purposes.
    LOG_ALLOW(LOCAL,LOG_DEBUG, "Face Distances:\n");
    if(get_log_level() == LOG_DEBUG && is_function_allowed(__func__)) LOG_FACE_DISTANCES(d); 
    CHKERRQ(ierr); // Check for errors in printing distances.

    // Determine the particle's position relative to the cell based on the computed distances.
    // The function sets the value pointed to by 'position' accordingly:
    // 0 for inside, 1 (2 or 3) for on the boundary, and -1 for outside.
    ierr = DeterminePointPosition(d,position); 
    CHKERRQ(ierr); // Check for errors in position determination.

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0); // Indicate successful execution of the function.
}

#undef __FUNCT__
#define __FUNCT__ "UpdateCellIndicesBasedOnDistances"
/**
 * @brief Internal helper implementation: `UpdateCellIndicesBasedOnDistances()`.
 * @details Local to this translation unit.
 */
PetscErrorCode UpdateCellIndicesBasedOnDistances( PetscReal d[NUM_FACES], PetscInt *idx, PetscInt *idy, PetscInt *idz)
{

   PetscFunctionBeginUser;
   PROFILE_FUNCTION_BEGIN;
  /*
    PetscInt cxm,cxs;  // maximum & minimum cell ID in x
    PetscInt cym,cys;  // maximum & minimum cell ID in y
    PetscInt czm,czs;  // maximum & minimum cell ID in z

    cxs = info->xs; cxm = cxs + info->xm - 2;
    cys = info->ys; cym = cys + info->ym - 2;
    czs = info->zs; czm = czs + info->zm - 2; 
  */

  //    LOG_ALLOW(LOCAL, LOG_DEBUG, "Received d: "
  //    "d[LEFT=%d]=%.3e, d[RIGHT=%d]=%.3e, d[BOTTOM=%d]=%.3e, "
  //	      "d[TOP=%d]=%.3e, d[FRONT=%d]=%.3e, d[BACK=%d]=%.3e\n",
		//	      LEFT, d[LEFT], RIGHT, d[RIGHT], BOTTOM, d[BOTTOM],
		//  TOP, d[TOP], FRONT, d[FRONT], BACK, d[BACK]);
    //  LOG_ALLOW(LOCAL, LOG_DEBUG, "Raw d: "
  //	      "[%.3e, %.3e, %.3e, %.3e, %.3e, %.3e]\n",
  //	      d[0], d[1], d[2], d[3], d[4], d[5]);
    
    // Validate input pointers
    if (d == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "'d' is NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Input array 'd' is NULL.");
    }
    if (idx == NULL || idy == NULL || idz == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "One or more index pointers are NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "One or more index pointers are NULL.");
    }

    // Debug: Print current face distances
    LOG_ALLOW(LOCAL,LOG_DEBUG, "Current Face Distances:\n");
    if(get_log_level() == LOG_DEBUG && is_function_allowed(__func__)) LOG_FACE_DISTANCES(d);

    // Update k-direction based on FRONT and BACK distances
    if (d[FRONT] < 0.0) {
      LOG_ALLOW(LOCAL,LOG_VERBOSE, "Condition met: d[FRONT] < 0.0, incrementing idz.\n");
        (*idz) += 1;
    }
    else if(d[BACK] < 0.0){
      LOG_ALLOW(LOCAL,LOG_VERBOSE, "Condition met: d[BACK] < 0.0, decrementing idz.\n");
        (*idz) -= 1;
    }

    // Update i-direction based on LEFT and RIGHT distances
    if (d[LEFT] < 0.0) {
      LOG_ALLOW(LOCAL,LOG_VERBOSE, "Condition met: d[LEFT] < 0.0, decrementing idx.\n");
        (*idx) -= 1;
    }
    else if (d[RIGHT] < 0.0) {
      LOG_ALLOW(LOCAL,LOG_VERBOSE, "Condition met: d[RIGHT] < 0.0, incrementing idx.\n");
        (*idx) += 1;
    }

    // Update j-direction based on BOTTOM and TOP distances
    if (d[BOTTOM] < 0.0) {
      LOG_ALLOW(LOCAL,LOG_VERBOSE, "Condition met: d[BOTTOM] < 0.0, decrementing idy.\n");
        (*idy) -= 1;
    }
    else if (d[TOP] < 0.0) {
      LOG_ALLOW(LOCAL,LOG_VERBOSE, "Condition met: d[TOP] < 0.0, incrementing idy.\n");
        (*idy) += 1;
    }

    /*
    // The 'cell' corners you can reference go from [xs .. xs+xm-1], but
    // to form a valid cell in x, you need (idx+1) in range, so max is (xs+xm-2).
    *idx = PetscMax(cxs,               PetscMin(*idx, cxm));
    *idy = PetscMax(cys,               PetscMin(*idy, cym));
    *idz = PetscMax(czs,               PetscMin(*idz, czm));
    */

    LOG_ALLOW(LOCAL,LOG_DEBUG, "Updated Indices  - idx, idy, idz: %d, %d, %d\n", *idx, *idy, *idz);

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0); // Indicate successful execution
}


#undef __FUNCT__
#define __FUNCT__ "FinalizeTraversal"
/**
 * @brief Implementation of \ref FinalizeTraversal().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/walkingsearch.h`.
 * @see FinalizeTraversal()
 */
PetscErrorCode FinalizeTraversal(UserCtx *user, Particle *particle, PetscInt traversal_steps, PetscBool cell_found, PetscInt idx, PetscInt idy, PetscInt idz)
{
    PetscFunctionBeginUser;
    // Validate input pointers
    if (user == NULL || particle == NULL) {
      LOG_ALLOW(LOCAL,LOG_ERROR, "One or more input pointers are NULL.\n");
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "One or more input pointers are NULL.");
    }

    if (cell_found) {
      LOG_ALLOW(LOCAL,LOG_INFO, "Particle located in cell (%d, %d, %d) after %d traversal steps.\n",
            idx, idy, idz, traversal_steps);
    }
    else {
      LOG_ALLOW(LOCAL,LOG_WARNING, "Particle could not be located within the grid after %d traversal steps.\n", (PetscInt)traversal_steps);
        particle->cell[0] = -1;
        particle->cell[1] = -1;
        particle->cell[2] = -1;
    }

    LOG_ALLOW(LOCAL, LOG_INFO, "Completed final traversal sync across all ranks.\n");


    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "FindOwnerOfCell"
/**
 * @brief Internal helper implementation: `FindOwnerOfCell()`.
 * @details Local to this translation unit.
 */
PetscErrorCode FindOwnerOfCell(UserCtx *user, PetscInt i, PetscInt j, PetscInt k, PetscMPIInt *owner_rank)
{
    PetscErrorCode ierr;
    PetscMPIInt size;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(ierr);

    // --- 1. Input Validation ---
    if (!user || !user->RankCellInfoMap) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "UserCtx or RankCellInfoMap is not initialized in FindOwnerOfCell.");
    }
    if (!owner_rank) {
        SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_NULL, "Output pointer owner_rank is NULL in FindOwnerOfCell.");
    }

    // --- 2. Linear Search through the Decomposition Map ---
    // Initialize to a "not found" state.
    *owner_rank = -1;

    // Loop through the map, which contains the ownership info for every rank 'r'.
    for (PetscMPIInt r = 0; r < size; ++r) {
        const RankCellInfo *info = &user->RankCellInfoMap[r];

        // A rank owns a cell if the cell's index is within its start (inclusive)
        // and end (exclusive) range for all three dimensions.
        if ((i >= info->xs_cell && i < info->xs_cell + info->xm_cell) &&
            (j >= info->ys_cell && j < info->ys_cell + info->ym_cell) &&
            (k >= info->zs_cell && k < info->zs_cell + info->zm_cell))
        {
            *owner_rank = r; // We found the owner.
            break;           // The search is over, exit the loop.
        }
    }

    // --- 3. Logging for Diagnostics ---
    // This is extremely useful for debugging particle migration issues.
    if (*owner_rank == -1) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "No owner found for global cell (%d, %d, %d). It is likely outside the domain.\n", i, j, k);
    } else {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "Owner of cell (%d, %d, %d) is Rank %d.\n", i, j, k, *owner_rank);
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "LocateParticleOrFindMigrationTarget"

/**
 * @brief Implementation of \ref LocateParticleOrFindMigrationTarget().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/walkingsearch.h`.
 * @see LocateParticleOrFindMigrationTarget()
 */
PetscErrorCode LocateParticleOrFindMigrationTarget(UserCtx *user,
                                                        Particle *particle,
                                                        ParticleLocationStatus *status_out)
{
    PetscErrorCode ierr;
    PetscMPIInt    rank;
    
    // --- Search State Variables ---
    PetscInt  idx, idy, idz;           // Current search cell global indices
    PetscInt  traversal_steps;         // Counter to prevent infinite loops
    PetscBool search_concluded = PETSC_FALSE; // Flag to terminate the main while loop

    // --- Oscillation/Stuck Loop Detection Variables ---
    PetscInt  repeatedIndexCount = 0;
    PetscInt  prevIdx = PETSC_MIN_INT, prevIdy = PETSC_MIN_INT, prevIdz = PETSC_MIN_INT;
    PetscInt  last_position_result = -999;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);
    if (user && user->simCtx) {
        user->simCtx->searchMetrics.searchAttempts++;
    }

    // --- 1. Initialize the Search ---
    ierr = InitializeTraversalParameters(user, particle, &idx, &idy, &idz, &traversal_steps); CHKERRQ(ierr);

    LOG_ALLOW(LOCAL,LOG_VERBOSE, " The Threshold for considering a particle to be at a face is %.16f.\n",DISTANCE_THRESHOLD);
    
    LOG_ALLOW(LOCAL,LOG_TRACE," [PID %lld]Traversal Initiated at : i = %d, j = %d, k = %d.\n",(long long)particle->PID,idx,idy,idz); 
    
    // --- 2. Main Walking Search Loop ---
    while (!search_concluded && traversal_steps < MAX_TRAVERSAL) {
        traversal_steps++;

        // --- 2a. GLOBAL Domain Boundary Check ---
        // IMPROVED: Instead of failing immediately when hitting a boundary, clamp the indices
        // to the boundary and continue searching. This allows the search to explore other
        // directions that might lead to the particle.
        PetscBool hit_boundary = PETSC_FALSE;
        if (idx < 0) {
            idx = 0;
            hit_boundary = PETSC_TRUE;
            user->simCtx->searchMetrics.boundaryClampCount++;
            LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %lld]: Hit global i-min boundary, clamped to idx=0.\n", (long long)particle->PID);
        }
        if (idx >= (user->IM - 1)) {
            idx = user->IM - 2;
            hit_boundary = PETSC_TRUE;
            user->simCtx->searchMetrics.boundaryClampCount++;
            LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %lld]: Hit global i-max boundary, clamped to idx=%d.\n", (long long)particle->PID, idx);
        }
        if (idy < 0) {
            idy = 0;
            hit_boundary = PETSC_TRUE;
            user->simCtx->searchMetrics.boundaryClampCount++;
            LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %lld]: Hit global j-min boundary, clamped to idy=0.\n", (long long)particle->PID);
        }
        if (idy >= (user->JM - 1)) {
            idy = user->JM - 2;
            hit_boundary = PETSC_TRUE;
            user->simCtx->searchMetrics.boundaryClampCount++;
            LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %lld]: Hit global j-max boundary, clamped to idy=%d.\n", (long long)particle->PID, idy);
        }
        if (idz < 0) {
            idz = 0;
            hit_boundary = PETSC_TRUE;
            user->simCtx->searchMetrics.boundaryClampCount++;
            LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %lld]: Hit global k-min boundary, clamped to idz=0.\n", (long long)particle->PID);
        }
        if (idz >= (user->KM - 1)) {
            idz = user->KM - 2;
            hit_boundary = PETSC_TRUE;
            user->simCtx->searchMetrics.boundaryClampCount++;
            LOG_ALLOW(LOCAL, LOG_DEBUG, "[PID %lld]: Hit global k-max boundary, clamped to idz=%d.\n", (long long)particle->PID, idz);
        }

        if (hit_boundary) {
            LOG_ALLOW(LOCAL, LOG_INFO, "[PID %lld]: Hit GLOBAL domain boundary. Clamped to boundary cell (%d,%d,%d) and continuing search.\n",
                      (long long)particle->PID, idx, idy, idz);
        }

        // --- 2b. LOCAL GHOST REGION CHECK (PREVENTS SEGV) ---
        // Before trying to access cell data, check if we have it.
        PetscBool is_cell_local;
        ierr = CheckCellWithinLocalGrid(user, idx, idy, idz, &is_cell_local); CHKERRQ(ierr);
        if (!is_cell_local) {
            // We have walked outside the local rank's ghost region.
            // This definitively means the particle belongs to another rank.
            // Conclude the search here; the current (idx,idy,idz) is the handoff cell.
            LOG_ALLOW(LOCAL, LOG_INFO, "[PID %lld]: Walked outside local ghost region to cell (%d,%d,%d). Concluding search for handoff.\n",
                      (long long)particle->PID, idx, idy, idz);
            search_concluded = PETSC_TRUE;
            continue; // Skip the rest of the loop; proceed to finalization.
        }
	
        // --- 2c. Stuck Loop Detection & Enhanced Tie-Breaker ---
        if (idx == prevIdx && idy == prevIdy && idz == prevIdz) {
            repeatedIndexCount++;
            if (repeatedIndexCount > REPEAT_COUNT_THRESHOLD) {
                // Only apply tie-breaker if we are stuck for the right reason (on a boundary)
                if (last_position_result >= 1) {
                    user->simCtx->searchMetrics.tieBreakCount++;
                    LOG_ALLOW(LOCAL, LOG_WARNING, "[PID %lld]: Stuck on boundary of cell (%d,%d,%d) for %d steps. Applying enhanced tie-breaker.\n",
                              (long long)particle->PID, idx, idy, idz, repeatedIndexCount);
                    
                    // Re-evaluate at the stuck cell to get definitive weights
                    Cell final_cell;
                    PetscReal final_d[NUM_FACES];
                    PetscInt final_position; // Dummy variable
                    
                    ierr = RetrieveCurrentCell(user, idx, idy, idz, &final_cell); CHKERRQ(ierr);
                    ierr = EvaluateParticlePosition(&final_cell, final_d, particle->loc, &final_position, DISTANCE_THRESHOLD);
                    
                    if (ierr == 0) { // If evaluation succeeded
                        ierr = UpdateParticleWeights(final_d, particle); CHKERRQ(ierr);
                        search_concluded = PETSC_TRUE; // Conclude search, accepting this cell.
                    } else { // Evaluation failed (e.g., degenerate cell)
                        LOG_ALLOW(LOCAL, LOG_ERROR, "[PID %lld]: Tie-breaker failed during final evaluation at cell (%d,%d,%d). Search fails.\n",
                                  (long long)particle->PID, idx, idy, idz);
                        idx = -1; // Invalidate result
                        search_concluded = PETSC_TRUE;
                    }
                } else { // Stuck for the wrong reason (not on a boundary)
                    LOG_ALLOW(LOCAL, LOG_ERROR, "[PID %lld]: Search is stuck at cell (%d,%d,%d) but not on a boundary. This indicates a logic error. Failing search.\n",
                              (long long)particle->PID, idx, idy, idz);
                    idx = -1; // Invalidate result
                    search_concluded = PETSC_TRUE;
                }
                if(search_concluded) continue;
            }
        } else {
            repeatedIndexCount = 0;
        }
        prevIdx = idx; prevIdy = idy; prevIdz = idz;

        // --- 2d. Geometric Evaluation ---
        Cell      current_cell;
        PetscReal distances[NUM_FACES];
        PetscInt  position_in_cell;

        ierr = RetrieveCurrentCell(user, idx, idy, idz, &current_cell); CHKERRQ(ierr);
        ierr = EvaluateParticlePosition(&current_cell, distances, particle->loc, &position_in_cell, DISTANCE_THRESHOLD); CHKERRQ(ierr);
        last_position_result = position_in_cell;

        // --- 2e. Decision Making ---
        if (position_in_cell >= 0) { // Particle is INSIDE or ON THE BOUNDARY
            search_concluded = PETSC_TRUE;
            ierr = UpdateParticleWeights(distances, particle); CHKERRQ(ierr);
        } else { // Particle is OUTSIDE
            ierr = UpdateCellIndicesBasedOnDistances(distances, &idx, &idy, &idz); CHKERRQ(ierr);
        }
    }

    user->simCtx->searchMetrics.traversalStepsSum += traversal_steps;
    user->simCtx->searchMetrics.maxTraversalSteps = PetscMax(user->simCtx->searchMetrics.maxTraversalSteps,
                                                             traversal_steps);

    // --- 3. Finalize and Determine Actionable Status ---
    if (idx == -1 || (!search_concluded && traversal_steps >= MAX_TRAVERSAL)) {
        if (idx != -1) {
            LOG_ALLOW(LOCAL, LOG_ERROR, "[PID %lld]: Search FAILED, exceeded MAX_TRAVERSAL limit of %d.\n",
                      (long long)particle->PID, MAX_TRAVERSAL);
        }
        *status_out = LOST;
        particle->cell[0] = -1; particle->cell[1] = -1; particle->cell[2] = -1;
    } else {
        // Search succeeded in finding a candidate cell, now determine its owner.
        PetscMPIInt owner_rank;
        ierr = FindOwnerOfCell(user, idx, idy, idz, &owner_rank); CHKERRQ(ierr);

	LOG_ALLOW(LOCAL,LOG_TRACE," [PID %ld] Owner rank : %d.\n",particle->PID,owner_rank);
	
        // Always update the particle's cell index. It's a good guess for the receiving rank.
        particle->cell[0] = idx; particle->cell[1] = idy; particle->cell[2] = idz;

        if (owner_rank == rank) {
            *status_out = ACTIVE_AND_LOCATED;
        } else if (owner_rank != -1) {
            // Particle belongs to another rank. Return the direct, actionable status.
            *status_out = MIGRATING_OUT;
            particle->destination_rank = owner_rank;
        } else { // Found a valid index, but no owner in the map.
            *status_out = LOST;
            particle->cell[0] = -1; particle->cell[1] = -1; particle->cell[2] = -1;
        }
    }

    LOG_ALLOW(LOCAL,LOG_DEBUG,"[Rank %d][PID %ld] Search complete.\n",rank,particle->PID);

    // --- 4. Report the Final Outcome ---
    ierr = ReportSearchOutcome(particle, *status_out, traversal_steps); CHKERRQ(ierr);
    
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ReportSearchOutcome"
/**
 * @brief Implementation of \ref ReportSearchOutcome().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/walkingsearch.h`.
 * @see ReportSearchOutcome()
 */
PetscErrorCode ReportSearchOutcome(const Particle *particle,
                                          ParticleLocationStatus status,
                                          PetscInt traversal_steps)
{
    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    switch (status) {
        case ACTIVE_AND_LOCATED:
            LOG_ALLOW(LOCAL, LOG_INFO, "Search SUCCESS [PID %lld]: Located in global cell (%d, %d, %d) after %d steps.\n",
                      (long long)particle->PID, particle->cell[0], particle->cell[1], particle->cell[2], traversal_steps);
            break;
        case MIGRATING_OUT:
            LOG_ALLOW(LOCAL, LOG_INFO, "Search SUCCESS [PID %lld]: Identified for migration to Rank %d. Handoff cell is (%d, %d, %d) after %d steps.\n",
                      (long long)particle->PID, particle->destination_rank, particle->cell[0], particle->cell[1], particle->cell[2], traversal_steps);
            break;
        case LOST:
            LOG_ALLOW(LOCAL, LOG_WARNING, "Search FAILED [PID %lld]: Particle is LOST after %d steps.\n",
                      (long long)particle->PID, traversal_steps);
            break;
        default:
            LOG_ALLOW(LOCAL, LOG_WARNING, "Search ended with unexpected status %d for PID %lld.\n", status, (long long)particle->PID);
            break;
    }
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
