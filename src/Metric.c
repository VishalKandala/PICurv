/* Metric.c ------------------------------------------------------------------
 * Utility routines for curvilinear-grid metric operations used by the
 * particle-swarm module.
 *
 *  –  Logical (xi,eta,zta) → physical (x,y,z) mapping via trilinear blend.
 *  –  Jacobian matrix and determinant for contravariant velocity conversion.
 *
 * The only data this file needs from the application is the DMDA that stores
 * vertex coordinates in the usual PETSc coordinate DM (user->da) and the
 * coordinate array type `Cmpnts` (three-component struct {x,y,z}).
 * ---------------------------------------------------------------------------*/

#include <petsc.h>
#include "Metric.h"          /* forward declarations + Cmpnts + UserCtx */


/* ------------------------------------------------------------------------- */
/**
 * @brief Implementation of \ref MetricGetCellVertices().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Metric.h`.
 * @see MetricGetCellVertices()
 */
PetscErrorCode MetricGetCellVertices(UserCtx *user,
                                     const Cmpnts ***X,   /* coord array */
                                     PetscInt i,PetscInt j,PetscInt k,
                                     Cmpnts V[8])
{
  (void)user;
  PetscFunctionBeginUser;
  PROFILE_HELPER_BEGIN;
  for (PetscInt c = 0; c < 8; ++c) {
    PetscInt ii = i + ((c & 1) ? 1 : 0);
    PetscInt jj = j + ((c & 2) ? 1 : 0);
    PetscInt kk = k + ((c & 4) ? 1 : 0);
    LOG_LOOP_ALLOW(GLOBAL, LOG_VERBOSE,i+j+k,10," ii: %d,jj:%d,kk:%d - Retrieved.\n",ii,jj,kk);
    V[c] = X[kk][jj][ii];
  }
  PROFILE_HELPER_END;
  PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref MetricFacePlane().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Metric.h`.
 * @see MetricFacePlane()
 */
PetscBool MetricFacePlane(const Cmpnts v1, const Cmpnts v2, const Cmpnts v3, const Cmpnts v4,
                          const Cmpnts cell_centroid, FacePlane *plane, PetscBool *degenerate)
{
    PROFILE_HELPER_BEGIN;
    const PetscReal edge1_x = v2.x - v1.x, edge1_y = v2.y - v1.y, edge1_z = v2.z - v1.z;
    const PetscReal edge2_x = v4.x - v1.x, edge2_y = v4.y - v1.y, edge2_z = v4.z - v1.z;
    PetscReal normal_x = edge1_y * edge2_z - edge1_z * edge2_y;
    PetscReal normal_y = edge1_z * edge2_x - edge1_x * edge2_z;
    PetscReal normal_z = edge1_x * edge2_y - edge1_y * edge2_x;
    const PetscReal normal_magnitude = sqrt(normal_x * normal_x + normal_y * normal_y + normal_z * normal_z);
    PetscBool ambiguous = PETSC_FALSE;

    plane->centroid.x = 0.25 * (v1.x + v2.x + v3.x + v4.x);
    plane->centroid.y = 0.25 * (v1.y + v2.y + v3.y + v4.y);
    plane->centroid.z = 0.25 * (v1.z + v2.z + v3.z + v4.z);
    *degenerate = (PetscBool)(normal_magnitude < 1.0e-12);
    if (!*degenerate) {
        // Orient outward: flip a normal that points from the face toward the cell centroid.
        const PetscReal dot_prod_orientation = normal_x * (cell_centroid.x - plane->centroid.x) +
                                               normal_y * (cell_centroid.y - plane->centroid.y) +
                                               normal_z * (cell_centroid.z - plane->centroid.z);
        if (dot_prod_orientation > 1.0e-9) {
            normal_x = -normal_x; normal_y = -normal_y; normal_z = -normal_z;
        } else if (dot_prod_orientation == 0.0) {
            ambiguous = PETSC_TRUE; // Keep the (v2-v1)x(v4-v1) direction
        }
        normal_x /= normal_magnitude;
        normal_y /= normal_magnitude;
        normal_z /= normal_magnitude;
    }
    plane->normal.x = normal_x; plane->normal.y = normal_y; plane->normal.z = normal_z;
    PROFILE_HELPER_END;
    return ambiguous;
}

/**
 * @brief Implementation of \ref MetricBuildCellFacePlanes().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Metric.h`.
 * @see MetricBuildCellFacePlanes()
 */
void MetricBuildCellFacePlanes(const Cell *cell, CellFacePlanes *planes)
{
    PROFILE_HELPER_BEGIN;
    // Face vertex orders: each face's v1, v2, v4 span its normal, as the walking search has always used.
    static const int face_vertices[NUM_FACES][4] = {
        [LEFT] = {0, 7, 4, 3}, [RIGHT] = {1, 2, 5, 6}, [BOTTOM] = {0, 1, 6, 7},
        [TOP] = {3, 4, 5, 2}, [FRONT] = {4, 7, 6, 5}, [BACK] = {0, 3, 2, 1}};
    Cmpnts   cell_centroid = {0.0, 0.0, 0.0};
    PetscInt mask = 0;

    for (int i = 0; i < 8; ++i) {
      cell_centroid.x += cell->vertices[i].x;
      cell_centroid.y += cell->vertices[i].y;
      cell_centroid.z += cell->vertices[i].z;
    }
    cell_centroid.x /= 8.0;
    cell_centroid.y /= 8.0;
    cell_centroid.z /= 8.0;
    for (int f = 0; f < NUM_FACES; ++f) {
        const int *v = face_vertices[f];
        PetscBool  degenerate;
        if (MetricFacePlane(cell->vertices[v[0]], cell->vertices[v[1]], cell->vertices[v[2]], cell->vertices[v[3]],
                            cell_centroid, &planes->face[f], &degenerate)) {
            LOG_ALLOW(LOCAL, LOG_WARNING, "Face %d: dot product for normal orientation is zero; keeping the (v2-v1)x(v4-v1) direction.\n", f);
        }
        if (degenerate) mask |= (PetscInt)1 << f;
    }
    planes->degenerate = mask;
    planes->size = 0.0;
    PROFILE_HELPER_END;
}

/** @brief The catalogued fields MetricComputeFacePlanes() owns. */
static const FieldId kFacePlaneFields[] = {
    FIELD_ID_PLANE_NORMAL_I, FIELD_ID_PLANE_NORMAL_J, FIELD_ID_PLANE_NORMAL_K,
    FIELD_ID_PLANE_CENTROID_I, FIELD_ID_PLANE_CENTROID_J, FIELD_ID_PLANE_CENTROID_K, FIELD_ID_SEARCH_CELL_SIZE,
    FIELD_ID_DUAL_PLANE_NORMAL_I, FIELD_ID_DUAL_PLANE_NORMAL_J, FIELD_ID_DUAL_PLANE_NORMAL_K,
    FIELD_ID_DUAL_PLANE_CENTROID_I, FIELD_ID_DUAL_PLANE_CENTROID_J, FIELD_ID_DUAL_PLANE_CENTROID_K};

/** @brief Offset from a face's lower-corner origin to where its family stores it (see Metric.h). */
static const PetscInt kPrimalFaceShift[3][3] = {{0, 1, 1}, {1, 0, 1}, {1, 1, 0}};

/**
 * @brief Implementation of \ref MetricComputeFacePlanes().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Metric.h`.
 * @see MetricComputeFacePlanes()
 */
PetscErrorCode MetricComputeFacePlanes(UserCtx *user)
{
    const DMDALocalInfo *info = &user->info;
    const PetscInt       gs[3] = {info->gxs, info->gys, info->gzs};
    const PetscInt       ge[3] = {info->gxs + info->gxm, info->gys + info->gym, info->gzs + info->gzm};
    Vec                  Coor;
    const PetscScalar   *src_flat[2] = {NULL, NULL};
    Vec                 *globals[2][6] = {
        {&user->PlaneNormalI, &user->PlaneNormalJ, &user->PlaneNormalK,
         &user->PlaneCentroidI, &user->PlaneCentroidJ, &user->PlaneCentroidK},
        {&user->DualPlaneNormalI, &user->DualPlaneNormalJ, &user->DualPlaneNormalK,
         &user->DualPlaneCentroidI, &user->DualPlaneCentroidJ, &user->DualPlaneCentroidK}};
    Vec                 *locals[2][6] = {
        {&user->lPlaneNormalI, &user->lPlaneNormalJ, &user->lPlaneNormalK,
         &user->lPlaneCentroidI, &user->lPlaneCentroidJ, &user->lPlaneCentroidK},
        {&user->lDualPlaneNormalI, &user->lDualPlaneNormalJ, &user->lDualPlaneNormalK,
         &user->lDualPlaneCentroidI, &user->lDualPlaneCentroidJ, &user->lDualPlaneCentroidK}};
    Cmpnts              *out[2][6];
    PetscReal           *size = NULL;

    PetscFunctionBeginUser;
    if (user->SearchCellSize) PetscFunctionReturn(0);
    PROFILE_FUNCTION_BEGIN;
    for (int which = 0; which < 2; ++which) {
        for (int f = 0; f < 6; ++f) {
            PetscCall(DMCreateGlobalVector(user->fda, globals[which][f]));
            PetscCall(DMCreateLocalVector(user->fda, locals[which][f]));
            PetscCall(VecSet(*locals[which][f], 0.0)); // A zero normal marks a face no plane was built for
            PetscCall(VecGetArray(*locals[which][f], (PetscScalar **)&out[which][f]));
        }
    }
    PetscCall(DMCreateGlobalVector(user->da, &user->SearchCellSize));
    PetscCall(DMCreateLocalVector(user->da, &user->lSearchCellSize));
    PetscCall(VecSet(user->lSearchCellSize, 0.0));
    PetscCall(VecGetArray(user->lSearchCellSize, &size));

    PetscCall(DMGetCoordinatesLocal(user->da, &Coor));
    PetscCall(VecGetArrayRead(Coor, &src_flat[0]));
    PetscCall(VecGetArrayRead(user->lCent, &src_flat[1]));
    const PetscBool periodic[3] = {(PetscBool)(user->simCtx->i_periodic != 0), (PetscBool)(user->simCtx->j_periodic != 0),
                                   (PetscBool)(user->simCtx->k_periodic != 0)};
    const PetscInt  mx[3] = {info->mx, info->my, info->mz};
    for (int which = 0; which < 2; ++which) {
        const Cmpnts *X = (const Cmpnts *)src_flat[which];
        // Points a cell may use, as the walk and the interpolation use them: physical nodes 0..mx-2
        // (the last DA point, and its wrapped ghost copies, are a dummy); cell centres 1..mx-2 across
        // a wall and 0..mx-1 on a periodic axis, where the boundary centres are the seam images.
        // No plane is oriented by a cell that touches any other point.
        PetscInt lo[3], hi[3];
        for (int a = 0; a < 3; ++a) {
            lo[a] = PetscMax(gs[a], (which && !periodic[a]) ? 1 : 0);
            hi[a] = PetscMin(ge[a] - 1, (which && periodic[a]) ? mx[a] - 1 : mx[a] - 2);
        }
        #define PICURV_AT(i, j, k) X[MetricCellFacePlanesIndex(info, (i), (j), (k))]
        for (PetscInt k = gs[2]; k < ge[2]; ++k) {
            for (PetscInt j = gs[1]; j < ge[1]; ++j) {
                for (PetscInt i = gs[0]; i < ge[0]; ++i) {
                    const PetscInt o[3] = {i, j, k};
                    for (PetscInt a = 0; a < 3; ++a) {
                        const PetscInt b = (a + 1) % 3, c = (a + 2) % 3;
                        // The face's own points, and the cell that orients it: the one above (origin o) when all
                        // its points are usable, else the one below. Without either the face keeps no plane.
                        if (o[b] < lo[b] || o[b] + 1 > hi[b] || o[c] < lo[c] || o[c] + 1 > hi[c]) continue;
                        PetscInt  cell[3] = {i, j, k};
                        PetscBool above = (PetscBool)(o[a] >= lo[a] && o[a] + 1 <= hi[a]);
                        if (!above) {
                            if (o[a] - 1 < lo[a] || o[a] > hi[a]) continue;
                            cell[a] -= 1;
                        }
                        Cell hex; // Vertex numbering of GetCellVerticesFromGrid()
                        hex.vertices[0] = PICURV_AT(cell[0],     cell[1],     cell[2]);
                        hex.vertices[1] = PICURV_AT(cell[0] + 1, cell[1],     cell[2]);
                        hex.vertices[2] = PICURV_AT(cell[0] + 1, cell[1] + 1, cell[2]);
                        hex.vertices[3] = PICURV_AT(cell[0],     cell[1] + 1, cell[2]);
                        hex.vertices[4] = PICURV_AT(cell[0],     cell[1] + 1, cell[2] + 1);
                        hex.vertices[5] = PICURV_AT(cell[0] + 1, cell[1] + 1, cell[2] + 1);
                        hex.vertices[6] = PICURV_AT(cell[0] + 1, cell[1],     cell[2] + 1);
                        hex.vertices[7] = PICURV_AT(cell[0],     cell[1],     cell[2] + 1);
                        Cmpnts centroid = {0.0, 0.0, 0.0};
                        for (int v = 0; v < 8; ++v) {
                            centroid.x += hex.vertices[v].x;
                            centroid.y += hex.vertices[v].y;
                            centroid.z += hex.vertices[v].z;
                        }
                        centroid.x /= 8.0; centroid.y /= 8.0; centroid.z /= 8.0;
                        // The face at o in the vertex order of the minus face (LEFT, BOTTOM, BACK) of cell o.
                        PetscInt e1[3] = {0, 0, 0}, e2[3] = {0, 0, 0};
                        if (a == 0) { e1[2] = 1; e2[1] = 1; }      // LEFT: o, o+k, o+j+k, o+j
                        else if (a == 1) { e1[0] = 1; e2[2] = 1; } // BOTTOM: o, o+i, o+i+k, o+k
                        else { e1[1] = 1; e2[0] = 1; }             // BACK: o, o+j, o+i+j, o+i
                        FacePlane plane;
                        PetscBool degenerate;
                        if (MetricFacePlane(PICURV_AT(i, j, k), PICURV_AT(i + e1[0], j + e1[1], k + e1[2]),
                                            PICURV_AT(i + e1[0] + e2[0], j + e1[1] + e2[1], k + e1[2] + e2[2]),
                                            PICURV_AT(i + e2[0], j + e2[1], k + e2[2]), centroid, &plane, &degenerate)) {
                            LOG_ALLOW(LOCAL, LOG_WARNING, "Face family %d at (%d,%d,%d): ambiguous orientation; keeping (v2-v1)x(v4-v1).\n",
                                      (int)a, (int)i, (int)j, (int)k);
                        }
                        if (degenerate) continue;
                        // Stored pointing toward +a: the minus face of the cell above, the plus face of the one below.
                        const PetscReal sign = above ? -1.0 : 1.0;
                        const PetscInt *sh = which ? (const PetscInt[3]){0, 0, 0} : kPrimalFaceShift[a];
                        const PetscInt  at = MetricCellFacePlanesIndex(info, i + sh[0], j + sh[1], k + sh[2]);
                        out[which][a][at].x = sign * plane.normal.x;
                        out[which][a][at].y = sign * plane.normal.y;
                        out[which][a][at].z = sign * plane.normal.z;
                        out[which][3 + a][at] = plane.centroid;
                    }
                    if (which == 0 && i >= lo[0] && i + 1 <= hi[0] && j >= lo[1] && j + 1 <= hi[1] && k >= lo[2] && k + 1 <= hi[2]) {
                        Cell hex;
                        hex.vertices[0] = PICURV_AT(i,     j,     k);
                        hex.vertices[1] = PICURV_AT(i + 1, j,     k);
                        hex.vertices[2] = PICURV_AT(i + 1, j + 1, k);
                        hex.vertices[3] = PICURV_AT(i,     j + 1, k);
                        hex.vertices[4] = PICURV_AT(i,     j + 1, k + 1);
                        hex.vertices[5] = PICURV_AT(i + 1, j + 1, k + 1);
                        hex.vertices[6] = PICURV_AT(i + 1, j,     k + 1);
                        hex.vertices[7] = PICURV_AT(i,     j,     k + 1);
                        PetscCall(GetCellCharacteristicSize(&hex, &size[MetricCellFacePlanesIndex(info, i + 1, j + 1, k + 1)]));
                    }
                }
            }
        }
        #undef PICURV_AT
    }
    PetscCall(VecRestoreArrayRead(user->lCent, &src_flat[1]));
    PetscCall(VecRestoreArrayRead(Coor, &src_flat[0]));

    // Owned values become the global fields; the catalogued ghost update then gives every rank the
    // owner's planes, which on every face a particle can reach are the ones it computed itself
    // (the periodic face repair writes only the ghosts beyond the physical range, which none reads).
    PetscCall(VecRestoreArray(user->lSearchCellSize, &size));
    for (int which = 0; which < 2; ++which) {
        for (int f = 0; f < 6; ++f) PetscCall(VecRestoreArray(*locals[which][f], (PetscScalar **)&out[which][f]));
    }
    for (size_t f = 0; f < sizeof(kFacePlaneFields) / sizeof(kFacePlaneFields[0]); ++f) {
        FieldView view;
        const PetscScalar *local = NULL;
        PetscScalar       *global = NULL;
        const PetscInt     dof = (kFacePlaneFields[f] == FIELD_ID_SEARCH_CELL_SIZE) ? 1 : 3;
        PetscCall(FieldGetView(user, kFacePlaneFields[f], &view));
        // Copy the owned block (DMLocalToGlobal cannot insert on a periodic axis held by one rank).
        PetscCall(VecGetArrayRead(view.local_vec, &local));
        PetscCall(VecGetArray(view.global_vec, &global));
        for (PetscInt k = info->zs, g = 0; k < info->zs + info->zm; ++k) {
            for (PetscInt j = info->ys; j < info->ys + info->ym; ++j) {
                for (PetscInt i = info->xs; i < info->xs + info->xm; ++i, ++g) {
                    PetscCall(PetscArraycpy(&global[dof * g], &local[dof * MetricCellFacePlanesIndex(info, i, j, k)], dof));
                }
            }
        }
        PetscCall(VecRestoreArray(view.global_vec, &global));
        PetscCall(VecRestoreArrayRead(view.local_vec, &local));
        PetscCall(UpdateLocalGhosts(user, kFacePlaneFields[f]));
    }
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref MetricDestroyFacePlanes().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Metric.h`.
 * @see MetricDestroyFacePlanes()
 */
PetscErrorCode MetricDestroyFacePlanes(UserCtx *user)
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    for (size_t f = 0; f < sizeof(kFacePlaneFields) / sizeof(kFacePlaneFields[0]); ++f) {
        const FieldDescriptor *descriptor = NULL;
        PetscCall(FieldGetDescriptor(kFacePlaneFields[f], &descriptor));
        PetscCall(VecDestroy((Vec *)((char *)user + descriptor->global_vec_offset)));
        PetscCall(VecDestroy((Vec *)((char *)user + descriptor->local_vec_offset)));
    }
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref MetricGetFacePlaneArrays().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Metric.h`.
 * @see MetricGetFacePlaneArrays()
 */
PetscErrorCode MetricGetFacePlaneArrays(UserCtx *user, PetscBool dual, FacePlaneArrays *arrays)
{
    Vec normals[3], centroids[3];

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCheck(user->lSearchCellSize, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
               "Particle face planes are not built: MetricComputeFacePlanes() runs in CalculateAllGridMetrics() "
               "for a particle run, collectively, once the grid is final.");
    if (dual) {
        normals[0] = user->lDualPlaneNormalI; normals[1] = user->lDualPlaneNormalJ; normals[2] = user->lDualPlaneNormalK;
        centroids[0] = user->lDualPlaneCentroidI; centroids[1] = user->lDualPlaneCentroidJ; centroids[2] = user->lDualPlaneCentroidK;
    } else {
        normals[0] = user->lPlaneNormalI; normals[1] = user->lPlaneNormalJ; normals[2] = user->lPlaneNormalK;
        centroids[0] = user->lPlaneCentroidI; centroids[1] = user->lPlaneCentroidJ; centroids[2] = user->lPlaneCentroidK;
    }
    arrays->dual = dual;
    for (int a = 0; a < 3; ++a) {
        PetscCall(VecGetArrayRead(normals[a], (const PetscScalar **)&arrays->normal[a]));
        PetscCall(VecGetArrayRead(centroids[a], (const PetscScalar **)&arrays->centroid[a]));
    }
    arrays->size = NULL;
    if (!dual) PetscCall(VecGetArrayRead(user->lSearchCellSize, &arrays->size));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref MetricRestoreFacePlaneArrays().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Metric.h`.
 * @see MetricRestoreFacePlaneArrays()
 */
PetscErrorCode MetricRestoreFacePlaneArrays(UserCtx *user, FacePlaneArrays *arrays)
{
    Vec normals[3], centroids[3];

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    if (arrays->dual) {
        normals[0] = user->lDualPlaneNormalI; normals[1] = user->lDualPlaneNormalJ; normals[2] = user->lDualPlaneNormalK;
        centroids[0] = user->lDualPlaneCentroidI; centroids[1] = user->lDualPlaneCentroidJ; centroids[2] = user->lDualPlaneCentroidK;
    } else {
        normals[0] = user->lPlaneNormalI; normals[1] = user->lPlaneNormalJ; normals[2] = user->lPlaneNormalK;
        centroids[0] = user->lPlaneCentroidI; centroids[1] = user->lPlaneCentroidJ; centroids[2] = user->lPlaneCentroidK;
    }
    for (int a = 0; a < 3; ++a) {
        PetscCall(VecRestoreArrayRead(normals[a], (const PetscScalar **)&arrays->normal[a]));
        PetscCall(VecRestoreArrayRead(centroids[a], (const PetscScalar **)&arrays->centroid[a]));
    }
    if (!arrays->dual) PetscCall(VecRestoreArrayRead(user->lSearchCellSize, &arrays->size));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref MetricGatherCellFacePlanes().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Metric.h`.
 * @see MetricGatherCellFacePlanes()
 */
void MetricGatherCellFacePlanes(const FacePlaneArrays *arrays, const DMDALocalInfo *info,
                                PetscInt i, PetscInt j, PetscInt k, CellFacePlanes *planes)
{
    PROFILE_HELPER_BEGIN;
    static const PetscInt dual_shift[3] = {0, 0, 0};
    static const Face     minus_face[3] = {LEFT, BOTTOM, BACK}, plus_face[3] = {RIGHT, TOP, FRONT};

    planes->degenerate = 0;
    for (PetscInt a = 0; a < 3; ++a) {
        const PetscInt *sh = arrays->dual ? dual_shift : kPrimalFaceShift[a];
        for (int side = 0; side < 2; ++side) {
            const PetscInt  at = MetricCellFacePlanesIndex(info, i + sh[0] + (a == 0 ? side : 0),
                                                           j + sh[1] + (a == 1 ? side : 0), k + sh[2] + (a == 2 ? side : 0));
            const Face      f  = side ? plus_face[a] : minus_face[a];
            const Cmpnts    n  = arrays->normal[a][at];
            FacePlane      *p  = &planes->face[f];
            // Stored toward +a: outward for the plus face, reversed (exactly) for the minus face.
            if (side) p->normal = n;
            else { p->normal.x = -n.x; p->normal.y = -n.y; p->normal.z = -n.z; }
            p->centroid = arrays->centroid[a][at];
            if (n.x == 0.0 && n.y == 0.0 && n.z == 0.0) planes->degenerate |= (PetscInt)1 << f;
        }
    }
    planes->size = arrays->dual ? 0.0 : arrays->size[MetricCellFacePlanesIndex(info, i + 1, j + 1, k + 1)];
    PROFILE_HELPER_END;
}


/* ------------------------------------------------------------------------- */
/**
 * @brief Blend eight corner values at the supplied trilinear reference coordinates.
 */
static inline void TrilinearBlend(const Cmpnts V[8],
                                  PetscReal xi,PetscReal eta,PetscReal zta,
                                  Cmpnts *Xp)
{
  PetscReal x=0,y=0,z=0;
  for (PetscInt c=0;c<8;++c) {
    PetscReal N = ((c&1)?xi : 1.0-xi ) *
                  ((c&2)?eta: 1.0-eta) *
                  ((c&4)?zta: 1.0-zta);
    x += N * V[c].x;
    y += N * V[c].y;
    z += N * V[c].z;
  }
  Xp->x = x; Xp->y = y; Xp->z = z;
}


/* ------------------------------------------------------------------------- */
/**
 * @brief Implementation of \ref MetricLogicalToPhysical().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Metric.h`.
 * @see MetricLogicalToPhysical()
 */
PetscErrorCode MetricLogicalToPhysical(UserCtx  *user,
                                       const Cmpnts ***X,
                                       PetscInt i,PetscInt j,PetscInt k,
                                       PetscReal xi,PetscReal eta,PetscReal zta,
                                       Cmpnts *Xp)
{
  PetscErrorCode ierr;
  Cmpnts V[8];
  PetscFunctionBeginUser;
  
  PROFILE_HELPER_BEGIN;

  ierr = MetricGetCellVertices(user,X,i,j,k,V); CHKERRQ(ierr);
  TrilinearBlend(V,xi,eta,zta,Xp);

  PROFILE_FUNCTION_END;

  PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------------- */
/**
 * @brief Implementation of \ref MetricJacobian().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Metric.h`.
 * @see MetricJacobian()
 */
PetscErrorCode MetricJacobian(UserCtx *user,
                              const Cmpnts ***X,
                              PetscInt i,PetscInt j,PetscInt k,
                              PetscReal xi,PetscReal eta,PetscReal zta,
                              PetscReal J[3][3], PetscReal *detJ)
{
  PetscErrorCode ierr;
  Cmpnts V[8];
  PetscFunctionBeginUser;

  PROFILE_HELPER_BEGIN;

  ierr = MetricGetCellVertices(user,X,i,j,k,V); CHKERRQ(ierr);

  /* derivatives of the [0,1] trilinear shape functions TrilinearBlend uses */
  PetscReal dN_dXi[8], dN_dEta[8], dN_dZta[8];
  for (PetscInt c=0;c<8;++c) {
    PetscReal sx = (c & 1) ?  1.0 : -1.0;
    PetscReal sy = (c & 2) ?  1.0 : -1.0;
    PetscReal sz = (c & 4) ?  1.0 : -1.0;
    dN_dXi [c] = sx * ( (c&2?eta:1-eta) ) * ( (c&4?zta:1-zta) );
    dN_dEta[c] = sy * ( (c&1?xi :1-xi ) ) * ( (c&4?zta:1-zta) );
    dN_dZta[c] = sz * ( (c&1?xi :1-xi ) ) * ( (c&2?eta:1-eta) );
  }

  /* assemble Jacobian */
  PetscReal x_xi=0,y_xi=0,z_xi=0,
            x_eta=0,y_eta=0,z_eta=0,
            x_zta=0,y_zta=0,z_zta=0;
  for (PetscInt c=0;c<8;++c) {
    x_xi  += dN_dXi [c]*V[c].x;  y_xi  += dN_dXi [c]*V[c].y;  z_xi  += dN_dXi [c]*V[c].z;
    x_eta += dN_dEta[c]*V[c].x;  y_eta += dN_dEta[c]*V[c].y;  z_eta += dN_dEta[c]*V[c].z;
    x_zta += dN_dZta[c]*V[c].x;  y_zta += dN_dZta[c]*V[c].y;  z_zta += dN_dZta[c]*V[c].z;
  }

  J[0][0]=x_xi;  J[0][1]=x_eta;  J[0][2]=x_zta;
  J[1][0]=y_xi;  J[1][1]=y_eta;  J[1][2]=y_zta;
  J[2][0]=z_xi;  J[2][1]=z_eta;  J[2][2]=z_zta;

  if (detJ) {
    *detJ = x_xi*(y_eta*z_zta - y_zta*z_eta)
          - x_eta*(y_xi*z_zta - y_zta*z_xi)
          + x_zta*(y_xi*z_eta - y_eta*z_xi);
  }

  PROFILE_FUNCTION_END;

  PetscFunctionReturn(0);
}


/* ------------------------------------------------------------------------- */
/**
 * @brief Implementation of \ref MetricVelocityContravariant().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Metric.h`.
 * @see MetricVelocityContravariant()
 */
PetscErrorCode MetricVelocityContravariant(const PetscReal J[3][3], PetscReal detJ,
                                           const PetscReal u[3], PetscReal uc[3])
{
  PetscFunctionBeginUser;

  PROFILE_HELPER_BEGIN;

  /* contravariant basis vectors (row of adjugate(J)) divided by detJ */
  PetscReal gxi[3]  = {  J[1][1]*J[2][2]-J[1][2]*J[2][1],
                        -J[0][1]*J[2][2]+J[0][2]*J[2][1],
                         J[0][1]*J[1][2]-J[0][2]*J[1][1] };
  PetscReal geta[3] = { -J[1][0]*J[2][2]+J[1][2]*J[2][0],
                         J[0][0]*J[2][2]-J[0][2]*J[2][0],
                        -J[0][0]*J[1][2]+J[0][2]*J[1][0] };
  PetscReal gzta[3] = {  J[1][0]*J[2][1]-J[1][1]*J[2][0],
                        -J[0][0]*J[2][1]+J[0][1]*J[2][0],
                         J[0][0]*J[1][1]-J[0][1]*J[1][0] };

  PetscReal invDet = 1.0 / detJ;
  for (int d=0; d<3; ++d) { gxi[d]  *= invDet; geta[d] *= invDet; gzta[d] *= invDet; }

  uc[0] = gxi [0]*u[0] + gxi [1]*u[1] + gxi [2]*u[2];
  uc[1] = geta[0]*u[0] + geta[1]*u[1] + geta[2]*u[2];
  uc[2] = gzta[0]*u[0] + gzta[1]*u[1] + gzta[2]*u[2];

  PROFILE_FUNCTION_END;

  PetscFunctionReturn(0);
}

/////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Internal helper implementation: `InvertCovariantMetricTensor()`.
 * @details Local to this translation unit.
 */
PetscErrorCode InvertCovariantMetricTensor(double covariantTensor[3][3], double contravariantTensor[3][3])
{
	PetscFunctionBeginUser;
	PROFILE_HELPER_BEGIN;

	const double a11=covariantTensor[0][0], a12=covariantTensor[0][1], a13=covariantTensor[0][2];
	const double a21=covariantTensor[1][0], a22=covariantTensor[1][1], a23=covariantTensor[1][2];
	const double a31=covariantTensor[2][0], a32=covariantTensor[2][1], a33=covariantTensor[2][2];

	double det = a11*(a33*a22-a32*a23) - a21*(a33*a12-a32*a13) + a31*(a23*a12-a22*a13);

	/* Degeneracy has to be judged against the rows' own magnitude, not against a fixed
	   number. The rows here are face-area vectors, so the determinant carries their
	   product - for an orthogonal cell it is exactly the square of the cell volume, and
	   it shrinks as the mesh is refined while the cell stays perfectly well formed. An
	   absolute floor therefore rejects fine grids: a wall-resolved cell of 5e-9 volume
	   has a determinant near 3e-17. The ratio below is 1 for an orthogonal cell and
	   approaches 0 only as the rows become coplanar, which is the actual failure. */
	const double row_scale = sqrt(a11*a11 + a12*a12 + a13*a13) *
	                         sqrt(a21*a21 + a22*a22 + a23*a23) *
	                         sqrt(a31*a31 + a32*a32 + a33*a33);

	if (row_scale <= 0.0 || fabs(det) <= 1.0e-10 * row_scale) {
		SETERRQ(PETSC_COMM_SELF, PETSC_ERR_MAT_LU_ZRPVT,
		        "Matrix is singular: |det| = %g is degenerate against the rows' own scale "
		        "%g (ratio %g); the rows are coplanar or one of them vanishes.",
		        (double)fabs(det), (double)row_scale,
		        (double)(row_scale > 0.0 ? fabs(det)/row_scale : 0.0));
	}

	contravariantTensor[0][0] = (a33*a22-a32*a23)/det;
    contravariantTensor[0][1] = -(a33*a12-a32*a13)/det;
    contravariantTensor[0][2] = (a23*a12-a22*a13)/det;
	contravariantTensor[1][0] = -(a33*a21-a31*a23)/det;
    contravariantTensor[1][1] = (a33*a11-a31*a13)/det;
    contravariantTensor[1][2] = -(a23*a11-a21*a13)/det;
	contravariantTensor[2][0] = (a32*a21-a31*a22)/det;
    contravariantTensor[2][1] = -(a32*a11-a31*a12)/det;
    contravariantTensor[2][2] = (a22*a11-a21*a12)/det;

	PROFILE_HELPER_END;
	PetscFunctionReturn(0);
}


/**
 * @brief Internal helper implementation: `CalculateFaceNormalAndArea()`.
 * @details Local to this translation unit.
 */
PetscErrorCode CalculateFaceNormalAndArea(Cmpnts csi, Cmpnts eta, Cmpnts zet, double ni[3], double nj[3], double nk[3], double *Ai, double *Aj, double *Ak)
{
	PetscFunctionBeginUser;
	PROFILE_HELPER_BEGIN;
	double g[3][3];
	double G[3][3];
	
	g[0][0]=csi.x, g[0][1]=csi.y, g[0][2]=csi.z;
	g[1][0]=eta.x, g[1][1]=eta.y, g[1][2]=eta.z;
	g[2][0]=zet.x, g[2][1]=zet.y, g[2][2]=zet.z;
	
	PetscCall(InvertCovariantMetricTensor(g, G));
	double xcsi=G[0][0], ycsi=G[1][0], zcsi=G[2][0];
	double xeta=G[0][1], yeta=G[1][1], zeta=G[2][1];
	double xzet=G[0][2], yzet=G[1][2], zzet=G[2][2];
	      
	double nx_i = xcsi, ny_i = ycsi, nz_i = zcsi;
	double nx_j = xeta, ny_j = yeta, nz_j = zeta;
	double nx_k = xzet, ny_k = yzet, nz_k = zzet;
	      
	double sum_i=sqrt(nx_i*nx_i+ny_i*ny_i+nz_i*nz_i);
	double sum_j=sqrt(nx_j*nx_j+ny_j*ny_j+nz_j*nz_j);
	double sum_k=sqrt(nx_k*nx_k+ny_k*ny_k+nz_k*nz_k);

	*Ai = sqrt( g[0][0]*g[0][0] + g[0][1]*g[0][1] + g[0][2]*g[0][2] );	// area
	*Aj = sqrt( g[1][0]*g[1][0] + g[1][1]*g[1][1] + g[1][2]*g[1][2] );
	*Ak =sqrt( g[2][0]*g[2][0] + g[2][1]*g[2][1] + g[2][2]*g[2][2] );
		
	nx_i /= sum_i, ny_i /= sum_i, nz_i /= sum_i;
	nx_j /= sum_j, ny_j /= sum_j, nz_j /= sum_j;
	nx_k /= sum_k, ny_k /= sum_k, nz_k /= sum_k;

	ni[0] = nx_i, ni[1] = ny_i, ni[2] = nz_i;
	nj[0] = nx_j, nj[1] = ny_j, nj[2] = nz_j;
	nk[0] = nx_k, nk[1] = ny_k, nk[2] = nz_k;

	PROFILE_FUNCTION_END;
	PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `ComputeCellCharacteristicLengthScale()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ComputeCellCharacteristicLengthScale(PetscReal ajc, Cmpnts csi, Cmpnts eta, Cmpnts zet, double *dx, double *dy, double *dz)
{
		PetscFunctionBeginUser;
		PROFILE_HELPER_BEGIN;
        double ni[3], nj[3], nk[3];
        double Li, Lj, Lk;
        double Ai, Aj, Ak;
        double vol = 1./ajc;

        PetscCall(CalculateFaceNormalAndArea(csi, eta, zet, ni, nj, nk, &Ai, &Aj, &Ak));
        Li = vol / Ai;
        Lj = vol / Aj;
        Lk = vol / Ak;

        // Length scale vector = di * ni_vector + dj * nj_vector + dk * nk_vector
        *dx = fabs( Li * ni[0] + Lj * nj[0] + Lk * nk[0] );
        *dy = fabs( Li * ni[1] + Lj * nj[1] + Lk * nk[1] );
        *dz = fabs( Li * ni[2] + Lj * nj[2] + Lk * nk[2] );

		PROFILE_FUNCTION_END;
		PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref ComputeCellDirectionalExtents().
 * @details Full API contract is documented with the header declaration in
 *          `include/Metric.h`.
 * @see ComputeCellDirectionalExtents()
 */
PetscErrorCode ComputeCellDirectionalExtents(PetscReal ajc, Cmpnts csi, Cmpnts eta, Cmpnts zet,
                                             double *l_xi, double *l_eta, double *l_zeta)
{
    const double area_xi   = sqrt(csi.x*csi.x + csi.y*csi.y + csi.z*csi.z);
    const double area_eta  = sqrt(eta.x*eta.x + eta.y*eta.y + eta.z*eta.z);
    const double area_zeta = sqrt(zet.x*zet.x + zet.y*zet.y + zet.z*zet.z);

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCheck(ajc > 0.0, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Cell Jacobian must be positive to define cell extents; received %g.", (double)ajc);
    PetscCheck(area_xi > 0.0 && area_eta > 0.0 && area_zeta > 0.0, PETSC_COMM_SELF,
               PETSC_ERR_ARG_OUTOFRANGE,
               "A cell face has zero area (%g, %g, %g), so the cell has no extent across it.",
               area_xi, area_eta, area_zeta);

    /* Volume over the area of a face pair is the distance between those faces: the
       cell's extent in that grid direction. It uses nothing but the cell's own metrics,
       so it is the same wherever and however the cell sits in space. */
    *l_xi   = (1.0/ajc)/area_xi;
    *l_eta  = (1.0/ajc)/area_eta;
    *l_zeta = (1.0/ajc)/area_zeta;
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref ComputeCellEdgeVectors().
 * @details Full API contract is documented with the header declaration in
 *          `include/Metric.h`.
 * @see ComputeCellEdgeVectors()
 */
PetscErrorCode ComputeCellEdgeVectors(PetscReal ajc, Cmpnts csi, Cmpnts eta, Cmpnts zet,
                                      Cmpnts edges[3])
{
    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCheck(ajc > 0.0, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
               "Cell Jacobian must be positive to define cell edges; received %g.", (double)ajc);
    /* The face-area vectors are J grad(xi), J grad(eta), J grad(zeta), with J the cell
       volume. The covariant basis dx/dxi equals J grad(eta) x grad(zeta), so it is
       (eta x zet) / J - the cell's edge along xi, direction and length together. */
    edges[0].x = ajc * (eta.y*zet.z - eta.z*zet.y);
    edges[0].y = ajc * (eta.z*zet.x - eta.x*zet.z);
    edges[0].z = ajc * (eta.x*zet.y - eta.y*zet.x);
    edges[1].x = ajc * (zet.y*csi.z - zet.z*csi.y);
    edges[1].y = ajc * (zet.z*csi.x - zet.x*csi.z);
    edges[1].z = ajc * (zet.x*csi.y - zet.y*csi.x);
    edges[2].x = ajc * (csi.y*eta.z - csi.z*eta.y);
    edges[2].y = ajc * (csi.z*eta.x - csi.x*eta.z);
    edges[2].z = ajc * (csi.x*eta.y - csi.y*eta.x);
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}


/* -------------------------------------------------------------------------- */
/**
 * @brief Internal helper implementation: `CheckAndFixGridOrientation()`.
 * @details Local to this translation unit.
 */
PetscErrorCode CheckAndFixGridOrientation(UserCtx *user)
{
    PetscErrorCode ierr;
    PetscReal      aj_min, aj_max;
    PetscMPIInt    rank;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    /* ---------------- step 1: global extrema of Aj ---------------- */
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(ierr);

    ierr = VecMin(user->Aj, NULL, &aj_min); CHKERRQ(ierr);  /* already global */
    ierr = VecMax(user->Aj, NULL, &aj_max); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO,
        "[orientation] Global Aj range: [%.3e , %.3e]\n",
        (double)aj_min, (double)aj_max);

    /* ---------------- step 2: detect malformed mesh ---------------- */
    if (aj_min < 0.0 && aj_max > 0.0)
        SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER,
            "Mixed Jacobian signs detected – grid is topologically inconsistent.");

    /* A uniformly left-handed grid is refused. The former "repair" negated Csi, Eta,
       Zet and Aj on this level only, left the face and cell-centre metrics as they
       were, and then reset the orientation flag to +1; a mirrored duct solved that way
       reached forty times its bulk velocity. Renumbering one logical axis is exact. */
    PetscInt orientation = +1;
    if (aj_max < 0.0) {
        SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER,
                "The grid is left-handed: every cell Jacobian is negative, as a mirrored or "
                "oddly permuted grid is. The solver needs right-handed logical axes. Renumber "
                "one logical axis (grid.gen: append reverse:axis=i|j|k to the transforms); "
                "that axis's two faces trade names.");
    }

    /* ---------------- step 4: store result in UserCtx -------------- */
    user->GridOrientation = orientation;

    if (!rank)
        LOG_ALLOW(LOCAL, LOG_INFO, "[orientation] Grid confirmed right-handed.\n");

    PROFILE_FUNCTION_END;        

    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `ApplyPeriodicCorrectionsToCellCentersAndSpacing()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ApplyPeriodicCorrectionsToCellCentersAndSpacing(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info = user->info;
    PetscInt       xs = info.xs, xe = info.xs + info.xm;
    PetscInt       ys = info.ys, ye = info.ys + info.ym;
    PetscInt       zs = info.zs, ze = info.zs + info.zm;
    PetscInt       mx = info.mx, my = info.my, mz = info.mz;
    const BCFace   neg_faces[3] = {BC_FACE_NEG_X, BC_FACE_NEG_Y, BC_FACE_NEG_Z};
    Cmpnts       ***cent, ***lcent;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;

    // The dummy slot before the first cell of a periodic axis holds the last cell's centre
    // shifted back by the validated seam translation T, and the slot after the last cell holds
    // the first cell's centre shifted forward by T; exact for any translated period. Axes are
    // corrected in turn, refreshing ghosts between them, so edge and corner dummies collect the
    // translation of every periodic axis they lie across.
    for (PetscInt axis = 0; axis < 3; axis++) {
        if (user->boundary_faces[neg_faces[axis]].mathematical_type != PERIODIC) continue;
        PetscCheck(user->periodic_translation_valid[axis], PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE,
                   "Periodic cell-centre correction on axis %d needs the validated seam translation.", (int)axis);
        const Cmpnts T = user->periodic_translation[axis];
        ierr = UpdateLocalGhosts(user, FIELD_ID_CENT); CHKERRQ(ierr);
        ierr = DMDAVecGetArray(user->fda, user->Cent, &cent); CHKERRQ(ierr);
        ierr = DMDAVecGetArray(user->fda, user->lCent, &lcent); CHKERRQ(ierr);
        for (PetscInt k = zs; k < ze; k++) {
            for (PetscInt j = ys; j < ye; j++) {
                for (PetscInt i = xs; i < xe; i++) {
                    const PetscInt idx[3] = {i, j, k}, n[3] = {mx, my, mz};
                    PetscInt src_i = i, src_j = j, src_k = k;
                    PetscReal sign;
                    if (idx[axis] == 0) sign = -1.0;              // before the first cell: last cell - T
                    else if (idx[axis] == n[axis] - 1) sign = 1.0; // after the last cell: first cell + T
                    else continue;
                    // Wrapped local ghost of the source cell: index -2 is the last physical cell,
                    // n + 1 the first.
                    if (axis == 0) src_i = (sign < 0.0) ? -2 : mx + 1;
                    if (axis == 1) src_j = (sign < 0.0) ? -2 : my + 1;
                    if (axis == 2) src_k = (sign < 0.0) ? -2 : mz + 1;
                    cent[k][j][i].x = lcent[src_k][src_j][src_i].x + sign * T.x;
                    cent[k][j][i].y = lcent[src_k][src_j][src_i].y + sign * T.y;
                    cent[k][j][i].z = lcent[src_k][src_j][src_i].z + sign * T.z;
                }
            }
        }
        ierr = DMDAVecRestoreArray(user->fda, user->lCent, &lcent); CHKERRQ(ierr);
        ierr = DMDAVecRestoreArray(user->fda, user->Cent, &cent); CHKERRQ(ierr);
    }

    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `ApplyPeriodicCorrectionsToIFaceCenter()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ApplyPeriodicCorrectionsToIFaceCenter(UserCtx *user)
{
    const FieldId fields[] = {FIELD_ID_CENTX};

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(SynchronizePeriodicFaceFields(user, 'i', 1, fields));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}


/**
 * @brief Internal helper implementation: `ApplyPeriodicCorrectionsToJFaceCenter()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ApplyPeriodicCorrectionsToJFaceCenter(UserCtx *user)
{
    const FieldId fields[] = {FIELD_ID_CENTY};

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(SynchronizePeriodicFaceFields(user, 'j', 1, fields));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}


/**
 * @brief Internal helper implementation: `ApplyPeriodicCorrectionsToKFaceCenter()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ApplyPeriodicCorrectionsToKFaceCenter(UserCtx *user)
{
    const FieldId fields[] = {FIELD_ID_CENTZ};

    PetscFunctionBeginUser;
    PROFILE_HELPER_BEGIN;
    PetscCall(SynchronizePeriodicFaceFields(user, 'k', 1, fields));
    PROFILE_HELPER_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `ComputeFaceMetrics()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ComputeFaceMetrics(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    Cmpnts       ***csi_arr, ***eta_arr, ***zet_arr;
    Cmpnts       ***nodal_coords_arr;
    Vec            localCoords_from_dm;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting calculation and update for Csi, Eta, Zet.\n");

    ierr = DMDAGetLocalInfo(user->fda, &info); CHKERRQ(ierr);

    // --- 1. Get Nodal Physical Coordinates (Local Ghosted Array directly) ---
    ierr = DMGetCoordinatesLocal(user->da, &localCoords_from_dm); CHKERRQ(ierr);
    if (!localCoords_from_dm) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE, "DMGetCoordinatesLocal failed to return a coordinate vector. \n");
    ierr = DMDAVecGetArrayRead(user->fda, localCoords_from_dm, &nodal_coords_arr); CHKERRQ(ierr);

    // --- 2. Get arrays for output global Vecs from UserCtx ---
    ierr = DMDAVecGetArray(user->fda, user->Csi, &csi_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Eta, &eta_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Zet, &zet_arr); CHKERRQ(ierr);

    // Define owned node ranges (global indices)
    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;

    // Global domain dimensions (total number of nodes)
    PetscInt mx = info.mx;
    PetscInt my = info.my;
    PetscInt mz = info.mz;

    // --- 3. Calculate Csi, Eta, Zet for INTERIOR Stencils ---
    // Start loops from 1 if at global boundary 0 to ensure k_node-1 etc. are valid.
    PetscInt k_loop_start = (zs == 0) ? zs + 1 : zs;
    PetscInt j_loop_start = (ys == 0) ? ys + 1 : ys;
    PetscInt i_loop_start = (xs == 0) ? xs + 1 : xs;

   // These represent the surface area of the curvilinear cell face and the normal rotated such that the direction of increasing coordinate is maintained.
   // The metric vectors (Csi, Eta, Zet) are defined to point in the direction of their corresponding increasing computational coordinate.

    // Calculate Csi
    for (PetscInt k_node = k_loop_start; k_node < ze; ++k_node) {
        for (PetscInt j_node = j_loop_start; j_node < ye; ++j_node) {
            for (PetscInt i_node = xs; i_node < xe; ++i_node) {
	      
                PetscReal dx_deta = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].x + nodal_coords_arr[k_node-1][j_node][i_node].x - nodal_coords_arr[k_node][j_node-1][i_node].x - nodal_coords_arr[k_node-1][j_node-1][i_node].x);
                PetscReal dy_deta = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].y + nodal_coords_arr[k_node-1][j_node][i_node].y - nodal_coords_arr[k_node][j_node-1][i_node].y - nodal_coords_arr[k_node-1][j_node-1][i_node].y);
                PetscReal dz_deta = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].z + nodal_coords_arr[k_node-1][j_node][i_node].z - nodal_coords_arr[k_node][j_node-1][i_node].z - nodal_coords_arr[k_node-1][j_node-1][i_node].z);
                PetscReal dx_dzeta = 0.5 * (nodal_coords_arr[k_node][j_node-1][i_node].x + nodal_coords_arr[k_node][j_node][i_node].x - nodal_coords_arr[k_node-1][j_node-1][i_node].x - nodal_coords_arr[k_node-1][j_node][i_node].x);
                PetscReal dy_dzeta = 0.5 * (nodal_coords_arr[k_node][j_node-1][i_node].y + nodal_coords_arr[k_node][j_node][i_node].y - nodal_coords_arr[k_node-1][j_node-1][i_node].y - nodal_coords_arr[k_node-1][j_node][i_node].y);
                PetscReal dz_dzeta = 0.5 * (nodal_coords_arr[k_node][j_node-1][i_node].z + nodal_coords_arr[k_node][j_node][i_node].z - nodal_coords_arr[k_node-1][j_node-1][i_node].z - nodal_coords_arr[k_node-1][j_node][i_node].z);

		csi_arr[k_node][j_node][i_node].x = dy_deta * dz_dzeta - dz_deta * dy_dzeta;
                csi_arr[k_node][j_node][i_node].y = dz_deta * dx_dzeta - dx_deta * dz_dzeta;
                csi_arr[k_node][j_node][i_node].z = dx_deta * dy_dzeta - dy_deta * dx_dzeta;
            }
        }
    }

    // Calculate Eta
    for (PetscInt k_node = k_loop_start; k_node < ze; ++k_node) {
        for (PetscInt j_node = ys; j_node < ye; ++j_node) {
            for (PetscInt i_node = i_loop_start; i_node < xe; ++i_node) {
	      
                PetscReal dx_dxi = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].x + nodal_coords_arr[k_node-1][j_node][i_node].x - nodal_coords_arr[k_node][j_node][i_node-1].x - nodal_coords_arr[k_node-1][j_node][i_node-1].x);
                PetscReal dy_dxi = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].y + nodal_coords_arr[k_node-1][j_node][i_node].y - nodal_coords_arr[k_node][j_node][i_node-1].y - nodal_coords_arr[k_node-1][j_node][i_node-1].y);
                PetscReal dz_dxi = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].z + nodal_coords_arr[k_node-1][j_node][i_node].z - nodal_coords_arr[k_node][j_node][i_node-1].z - nodal_coords_arr[k_node-1][j_node][i_node-1].z);
                PetscReal dx_dzeta = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].x + nodal_coords_arr[k_node][j_node][i_node-1].x - nodal_coords_arr[k_node-1][j_node][i_node].x - nodal_coords_arr[k_node-1][j_node][i_node-1].x);
                PetscReal dy_dzeta = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].y + nodal_coords_arr[k_node][j_node][i_node-1].y - nodal_coords_arr[k_node-1][j_node][i_node].y - nodal_coords_arr[k_node-1][j_node][i_node-1].y);
                PetscReal dz_dzeta = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].z + nodal_coords_arr[k_node][j_node][i_node-1].z - nodal_coords_arr[k_node-1][j_node][i_node].z - nodal_coords_arr[k_node-1][j_node][i_node-1].z);

		eta_arr[k_node][j_node][i_node].x = dy_dzeta * dz_dxi  - dz_dzeta * dy_dxi;
                eta_arr[k_node][j_node][i_node].y = dz_dzeta * dx_dxi  - dx_dzeta * dz_dxi;
                eta_arr[k_node][j_node][i_node].z = dx_dzeta * dy_dxi  - dy_dzeta * dx_dxi;
            }
        }
    }

    // Calculate Zet
    for (PetscInt k_node = zs; k_node < ze; ++k_node) {
        for (PetscInt j_node = j_loop_start; j_node < ye; ++j_node) {
            for (PetscInt i_node = i_loop_start; i_node < xe; ++i_node) {
	      
                PetscReal dx_dxi = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].x + nodal_coords_arr[k_node][j_node-1][i_node].x - nodal_coords_arr[k_node][j_node][i_node-1].x - nodal_coords_arr[k_node][j_node-1][i_node-1].x);
                PetscReal dy_dxi = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].y + nodal_coords_arr[k_node][j_node-1][i_node].y - nodal_coords_arr[k_node][j_node][i_node-1].y - nodal_coords_arr[k_node][j_node-1][i_node-1].y);
                PetscReal dz_dxi = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].z + nodal_coords_arr[k_node][j_node-1][i_node].z - nodal_coords_arr[k_node][j_node][i_node-1].z - nodal_coords_arr[k_node][j_node-1][i_node-1].z);
                PetscReal dx_deta = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].x + nodal_coords_arr[k_node][j_node][i_node-1].x - nodal_coords_arr[k_node][j_node-1][i_node].x - nodal_coords_arr[k_node][j_node-1][i_node-1].x);
                PetscReal dy_deta = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].y + nodal_coords_arr[k_node][j_node][i_node-1].y - nodal_coords_arr[k_node][j_node-1][i_node].y - nodal_coords_arr[k_node][j_node-1][i_node-1].y);
                PetscReal dz_deta = 0.5 * (nodal_coords_arr[k_node][j_node][i_node].z + nodal_coords_arr[k_node][j_node][i_node-1].z - nodal_coords_arr[k_node][j_node-1][i_node].z - nodal_coords_arr[k_node][j_node-1][i_node-1].z);

		zet_arr[k_node][j_node][i_node].x = dy_dxi * dz_deta - dz_dxi * dy_deta;
                zet_arr[k_node][j_node][i_node].y = dz_dxi * dx_deta - dx_dxi * dz_deta;
                zet_arr[k_node][j_node][i_node].z = dx_dxi * dy_deta - dy_dxi * dx_deta;
            }
        }
    }

    // --- 4. Boundary Extrapolation ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Extrapolating boundary values for Csi, Eta, Zet.\n");
    PetscInt i_bnd, j_bnd, k_bnd;

    if (xs == 0) { // If this rank owns the global i=0 boundary
        i_bnd = 0;
        for (k_bnd = zs; k_bnd < ze; ++k_bnd) {
            for (j_bnd = ys; j_bnd < ye; ++j_bnd) {
                if (i_bnd + 1 < mx) {
                    eta_arr[k_bnd][j_bnd][i_bnd] = eta_arr[k_bnd][j_bnd][i_bnd+1];
                    zet_arr[k_bnd][j_bnd][i_bnd] = zet_arr[k_bnd][j_bnd][i_bnd+1];
                }
            }
        }
    }
    if (xe == mx) { // If this rank owns the global i=mx-1 boundary
        i_bnd = mx - 1;
        for (k_bnd = zs; k_bnd < ze; ++k_bnd) {
            for (j_bnd = ys; j_bnd < ye; ++j_bnd) {
                if (i_bnd - 1 >= 0) {
                    eta_arr[k_bnd][j_bnd][i_bnd] = eta_arr[k_bnd][j_bnd][i_bnd-1];
                    zet_arr[k_bnd][j_bnd][i_bnd] = zet_arr[k_bnd][j_bnd][i_bnd-1];
                }
            }
        }
    }
    if (ys == 0) {
        j_bnd = 0;
        for (k_bnd = zs; k_bnd < ze; ++k_bnd) {
            for (i_bnd = xs; i_bnd < xe; ++i_bnd) {
                if (j_bnd + 1 < my) {
                    csi_arr[k_bnd][j_bnd][i_bnd] = csi_arr[k_bnd][j_bnd+1][i_bnd];
                    zet_arr[k_bnd][j_bnd][i_bnd] = zet_arr[k_bnd][j_bnd+1][i_bnd];
                }
            }
        }
    }
    if (ye == my) {
        j_bnd = my - 1;
        for (k_bnd = zs; k_bnd < ze; ++k_bnd) {
            for (i_bnd = xs; i_bnd < xe; ++i_bnd) {
                if (j_bnd - 1 >= 0) {
                    csi_arr[k_bnd][j_bnd][i_bnd] = csi_arr[k_bnd][j_bnd-1][i_bnd];
                    zet_arr[k_bnd][j_bnd][i_bnd] = zet_arr[k_bnd][j_bnd-1][i_bnd];
                }
            }
        }
    }
    if (zs == 0) {
        k_bnd = 0;
        for (j_bnd = ys; j_bnd < ye; ++j_bnd) {
            for (i_bnd = xs; i_bnd < xe; ++i_bnd) {
                if (k_bnd + 1 < mz) {
                    csi_arr[k_bnd][j_bnd][i_bnd] = csi_arr[k_bnd+1][j_bnd][i_bnd];
                    eta_arr[k_bnd][j_bnd][i_bnd] = eta_arr[k_bnd+1][j_bnd][i_bnd];
                }
            }
        }
    }
    if (ze == mz) {
        k_bnd = mz - 1;
        for (j_bnd = ys; j_bnd < ye; ++j_bnd) {
            for (i_bnd = xs; i_bnd < xe; ++i_bnd) {
                if (k_bnd - 1 >= 0) {
                    csi_arr[k_bnd][j_bnd][i_bnd] = csi_arr[k_bnd-1][j_bnd][i_bnd];
                    eta_arr[k_bnd][j_bnd][i_bnd] = eta_arr[k_bnd-1][j_bnd][i_bnd];
                }
            }
        }
    }

    if (info.xs==0 && info.ys==0 && info.zs==0) {
      PetscReal dot = zet_arr[0][0][0].z;  /* dot with global +z */
      LOG_ALLOW(GLOBAL,LOG_DEBUG,"Zet(k=0)·ez = %.3f   (should be >0 for right-handed grid)\n", dot);
    }
    
    // --- 5. Restore all arrays ---
    ierr = DMDAVecRestoreArrayRead(user->fda, localCoords_from_dm, &nodal_coords_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Csi, &csi_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Eta, &eta_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Zet, &zet_arr); CHKERRQ(ierr);

    // --- 6. Assemble Global Vectors ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Assembling global Csi, Eta, Zet.\n");
    ierr = VecAssemblyBegin(user->Csi); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->Csi); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->Eta); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->Eta); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->Zet); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->Zet); CHKERRQ(ierr);

    // --- 7. Update Local Ghosted Versions ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Updating local lCsi, lEta, lZet.\n");
    ierr = UpdateLocalGhosts(user, FIELD_ID_CSI); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_ETA); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_ZET); CHKERRQ(ierr);
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Completed calculation, extrapolation, and update for Csi, Eta, Zet.\n");

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}


/**
 * @brief Implementation of \ref ComputeCellCenteredJacobianInverse().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/Metric.h`.
 * @see ComputeCellCenteredJacobianInverse()
 */
PetscErrorCode ComputeCellCenteredJacobianInverse(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    PetscScalar  ***aj_arr;
    Cmpnts       ***nodal_coords_arr;
    Vec            localCoords_from_dm;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting calculation, extrapolation, and update for Aj.\n");

    // --- 1. Get Nodal Coordinates and Output Array ---
    ierr = DMGetCoordinatesLocal(user->da, &localCoords_from_dm); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, localCoords_from_dm, &nodal_coords_arr); CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->Aj, &aj_arr); CHKERRQ(ierr);

    // Define owned node ranges (global indices)
    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;

    // Global domain dimensions (total number of nodes)
    PetscInt mx = info.mx;
    PetscInt my = info.my;
    PetscInt mz = info.mz;

    // --- 2. Calculate Aj for INTERIOR Stencils ---
    
    PetscInt k_start_node = (zs == 0) ? zs + 1 : zs;
    PetscInt j_start_node = (ys == 0) ? ys + 1 : ys;
    PetscInt i_start_node = (xs == 0) ? xs + 1 : xs;

    PetscInt k_end_node = (ze == mz) ? ze - 1 : ze;
    PetscInt j_end_node = (ye == my) ? ye - 1 : ye;
    PetscInt i_end_node = (xe == mx) ? xe - 1 : xe;

    for (PetscInt k_node = k_start_node; k_node < k_end_node; ++k_node) {
        for (PetscInt j_node = j_start_node; j_node < j_end_node; ++j_node) {
            for (PetscInt i_node = i_start_node; i_node < i_end_node; ++i_node) {
	      
                PetscReal dx_dxi = 0.25 * ( (nodal_coords_arr[k_node][j_node][i_node].x + nodal_coords_arr[k_node][j_node-1][i_node].x + nodal_coords_arr[k_node-1][j_node][i_node].x + nodal_coords_arr[k_node-1][j_node-1][i_node].x) - (nodal_coords_arr[k_node][j_node][i_node-1].x + nodal_coords_arr[k_node][j_node-1][i_node-1].x + nodal_coords_arr[k_node-1][j_node][i_node-1].x + nodal_coords_arr[k_node-1][j_node-1][i_node-1].x) );

		PetscReal dy_dxi = 0.25 * ( (nodal_coords_arr[k_node][j_node][i_node].y + nodal_coords_arr[k_node][j_node-1][i_node].y + nodal_coords_arr[k_node-1][j_node][i_node].y + nodal_coords_arr[k_node-1][j_node-1][i_node].y) - (nodal_coords_arr[k_node][j_node][i_node-1].y + nodal_coords_arr[k_node][j_node-1][i_node-1].y + nodal_coords_arr[k_node-1][j_node][i_node-1].y + nodal_coords_arr[k_node-1][j_node-1][i_node-1].y) );

		PetscReal dz_dxi = 0.25 * ( (nodal_coords_arr[k_node][j_node][i_node].z + nodal_coords_arr[k_node][j_node-1][i_node].z + nodal_coords_arr[k_node-1][j_node][i_node].z + nodal_coords_arr[k_node-1][j_node-1][i_node].z) - (nodal_coords_arr[k_node][j_node][i_node-1].z + nodal_coords_arr[k_node][j_node-1][i_node-1].z + nodal_coords_arr[k_node-1][j_node][i_node-1].z + nodal_coords_arr[k_node-1][j_node-1][i_node-1].z) );

		PetscReal dx_deta = 0.25 * ( (nodal_coords_arr[k_node][j_node][i_node].x + nodal_coords_arr[k_node][j_node][i_node-1].x + nodal_coords_arr[k_node-1][j_node][i_node].x + nodal_coords_arr[k_node-1][j_node][i_node-1].x) - (nodal_coords_arr[k_node][j_node-1][i_node].x + nodal_coords_arr[k_node][j_node-1][i_node-1].x + nodal_coords_arr[k_node-1][j_node-1][i_node].x + nodal_coords_arr[k_node-1][j_node-1][i_node-1].x) );

		PetscReal dy_deta = 0.25 * ( (nodal_coords_arr[k_node][j_node][i_node].y + nodal_coords_arr[k_node][j_node][i_node-1].y + nodal_coords_arr[k_node-1][j_node][i_node].y + nodal_coords_arr[k_node-1][j_node][i_node-1].y) - (nodal_coords_arr[k_node][j_node-1][i_node].y + nodal_coords_arr[k_node][j_node-1][i_node-1].y + nodal_coords_arr[k_node-1][j_node-1][i_node].y + nodal_coords_arr[k_node-1][j_node-1][i_node-1].y) );

		PetscReal dz_deta = 0.25 * ( (nodal_coords_arr[k_node][j_node][i_node].z + nodal_coords_arr[k_node][j_node][i_node-1].z + nodal_coords_arr[k_node-1][j_node][i_node].z + nodal_coords_arr[k_node-1][j_node][i_node-1].z) - (nodal_coords_arr[k_node][j_node-1][i_node].z + nodal_coords_arr[k_node][j_node-1][i_node-1].z + nodal_coords_arr[k_node-1][j_node-1][i_node].z + nodal_coords_arr[k_node-1][j_node-1][i_node-1].z) );

		PetscReal dx_dzeta = 0.25 * ( (nodal_coords_arr[k_node][j_node][i_node].x + nodal_coords_arr[k_node][j_node-1][i_node].x + nodal_coords_arr[k_node][j_node][i_node-1].x + nodal_coords_arr[k_node][j_node-1][i_node-1].x) - (nodal_coords_arr[k_node-1][j_node][i_node].x + nodal_coords_arr[k_node-1][j_node-1][i_node].x + nodal_coords_arr[k_node-1][j_node][i_node-1].x + nodal_coords_arr[k_node-1][j_node-1][i_node-1].x) );

		PetscReal dy_dzeta = 0.25 * ( (nodal_coords_arr[k_node][j_node][i_node].y + nodal_coords_arr[k_node][j_node-1][i_node].y + nodal_coords_arr[k_node][j_node][i_node-1].y + nodal_coords_arr[k_node][j_node-1][i_node-1].y) - (nodal_coords_arr[k_node-1][j_node][i_node].y + nodal_coords_arr[k_node-1][j_node-1][i_node].y + nodal_coords_arr[k_node-1][j_node][i_node-1].y + nodal_coords_arr[k_node-1][j_node-1][i_node-1].y) );

		PetscReal dz_dzeta = 0.25 * ( (nodal_coords_arr[k_node][j_node][i_node].z + nodal_coords_arr[k_node][j_node-1][i_node].z + nodal_coords_arr[k_node][j_node][i_node-1].z + nodal_coords_arr[k_node][j_node-1][i_node-1].z) - (nodal_coords_arr[k_node-1][j_node][i_node].z + nodal_coords_arr[k_node-1][j_node-1][i_node].z + nodal_coords_arr[k_node-1][j_node][i_node-1].z + nodal_coords_arr[k_node-1][j_node-1][i_node-1].z) );

		PetscReal jacobian_det = dx_dxi * (dy_deta * dz_dzeta - dz_deta * dy_dzeta) - dy_dxi * (dx_deta * dz_dzeta - dz_deta * dx_dzeta) + dz_dxi * (dx_deta * dy_dzeta - dy_deta * dx_dzeta);
                if (PetscAbsReal(jacobian_det) < 1.0e-18) { SETERRQ(PETSC_COMM_SELF, PETSC_ERR_FLOP_COUNT, "Jacobian is near zero..."); }
                aj_arr[k_node][j_node][i_node] = 1.0 / jacobian_det;
            }
        }
    }

    // --- 4. Boundary Extrapolation for Aj ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Extrapolating boundary values for Aj. \n");
    PetscInt i_bnd, j_bnd, k_bnd;

    if (xs == 0) {
        i_bnd = 0;
        for (k_bnd = zs; k_bnd < ze; ++k_bnd) {
            for (j_bnd = ys; j_bnd < ye; ++j_bnd) {
                if (i_bnd + 1 < mx) aj_arr[k_bnd][j_bnd][i_bnd] = aj_arr[k_bnd][j_bnd][i_bnd+1];
            }
        }
    }
    if (xe == mx) {
        i_bnd = mx - 1;
        for (k_bnd = zs; k_bnd < ze; ++k_bnd) {
            for (j_bnd = ys; j_bnd < ye; ++j_bnd) {
                if (i_bnd - 1 >= 0) aj_arr[k_bnd][j_bnd][i_bnd] = aj_arr[k_bnd][j_bnd][i_bnd-1];
            }
        }
    }
    // (Similar extrapolation blocks for Y and Z boundaries for aj_arr)
    if (ys == 0) {
        j_bnd = 0;
        for (k_bnd = zs; k_bnd < ze; ++k_bnd) {
            for (i_bnd = xs; i_bnd < xe; ++i_bnd) {
                 if (j_bnd + 1 < my) aj_arr[k_bnd][j_bnd][i_bnd] = aj_arr[k_bnd][j_bnd+1][i_bnd];
            }
        }
    }
    if (ye == my) {
        j_bnd = my - 1;
        for (k_bnd = zs; k_bnd < ze; ++k_bnd) {
            for (i_bnd = xs; i_bnd < xe; ++i_bnd) {
                 if (j_bnd - 1 >= 0) aj_arr[k_bnd][j_bnd][i_bnd] = aj_arr[k_bnd][j_bnd-1][i_bnd];
            }
        }
    }
    if (zs == 0) {
        k_bnd = 0;
        for (j_bnd = ys; j_bnd < ye; ++j_bnd) {
            for (i_bnd = xs; i_bnd < xe; ++i_bnd) {
                if (k_bnd + 1 < mz) aj_arr[k_bnd][j_bnd][i_bnd] = aj_arr[k_bnd+1][j_bnd][i_bnd];
            }
        }
    }
    if (ze == mz) {
        k_bnd = mz - 1;
        for (j_bnd = ys; j_bnd < ye; ++j_bnd) {
            for (i_bnd = xs; i_bnd < xe; ++i_bnd) {
                 if (k_bnd - 1 >= 0) aj_arr[k_bnd][j_bnd][i_bnd] = aj_arr[k_bnd-1][j_bnd][i_bnd];
            }
        }
    }

    // --- 5. Restore arrays ---
    ierr = DMDAVecRestoreArrayRead(user->fda, localCoords_from_dm, &nodal_coords_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, user->Aj, &aj_arr); CHKERRQ(ierr);

    // --- 6. Assemble Global Vector ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Assembling global Aj.\n");
    ierr = VecAssemblyBegin(user->Aj); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(user->Aj); CHKERRQ(ierr);

    // --- 7. Update Local Ghosted Version ---
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Updating local lAj.\n");
    ierr = UpdateLocalGhosts(user, FIELD_ID_AJ); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Completed calculation, extrapolation, and update for Aj.\n");
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


/**
 * @brief Internal helper implementation: `ComputeCellCentersAndSpacing()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ComputeCellCentersAndSpacing(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    Vec            lCoords;
    const Cmpnts ***coor;
    Cmpnts       ***cent, ***gs;
    PetscReal    xcp, ycp, zcp, xcm, ycm, zcm;
    PetscInt     xs,ys,zs,xe,ye,ze,mx,my,mz;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Computing cell centers and spacing for level %d block %d...\n", user->simCtx->rank, user->thislevel, user->_this);

    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(user->da, &lCoords); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, lCoords, &coor); CHKERRQ(ierr);

    ierr = DMDAVecGetArray(user->fda, user->Cent, &cent); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->GridSpace, &gs); CHKERRQ(ierr);

    xs = info.xs; xe = info.xs + info.xm;
    ys = info.ys; ye = info.ys + info.ym;
    zs = info.zs; ze = info.zs + info.zm;
    mx = info.mx; my =  info.my; mz = info.mz;

    PetscInt k_start_node = (zs == 0) ? zs + 1 : zs;
    PetscInt j_start_node = (ys == 0) ? ys + 1 : ys;
    PetscInt i_start_node = (xs == 0) ? xs + 1 : xs;

    PetscInt k_end_node = (ze == mz) ? ze - 1 : ze;
    PetscInt j_end_node = (ye == my) ? ye - 1 : ye;
    PetscInt i_end_node = (xe == mx) ? xe - 1 : xe;    
    
    // Loop over the interior OWNED cells (stencil requires i-1, j-1, k-1)
    for (PetscInt k=k_start_node; k<k_end_node; k++) {
        for (PetscInt j=j_start_node; j<j_end_node; j++) {
            for (PetscInt i=i_start_node; i<i_end_node; i++) {
                // Calculate cell center as the average of its 8 corner nodes
                cent[k][j][i].x = 0.125 * (coor[k][j][i].x + coor[k][j-1][i].x + coor[k-1][j][i].x + coor[k-1][j-1][i].x + coor[k][j][i-1].x + coor[k][j-1][i-1].x + coor[k-1][j][i-1].x + coor[k-1][j-1][i-1].x);
                cent[k][j][i].y = 0.125 * (coor[k][j][i].y + coor[k][j-1][i].y + coor[k-1][j][i].y + coor[k-1][j-1][i].y + coor[k][j][i-1].y + coor[k][j-1][i-1].y + coor[k-1][j][i-1].y + coor[k-1][j-1][i-1].y);
                cent[k][j][i].z = 0.125 * (coor[k][j][i].z + coor[k][j-1][i].z + coor[k-1][j][i].z + coor[k-1][j-1][i].z + coor[k][j][i-1].z + coor[k][j-1][i-1].z + coor[k-1][j][i-1].z + coor[k-1][j-1][i-1].z);

                // Calculate Grid Spacing in i-direction (distance between i-face centers)
                xcp = 0.25 * (coor[k][j][i].x + coor[k][j-1][i].x + coor[k-1][j-1][i].x + coor[k-1][j][i].x);
                ycp = 0.25 * (coor[k][j][i].y + coor[k][j-1][i].y + coor[k-1][j-1][i].y + coor[k-1][j][i].y);
                zcp = 0.25 * (coor[k][j][i].z + coor[k][j-1][i].z + coor[k-1][j-1][i].z + coor[k-1][j][i].z);
                xcm = 0.25 * (coor[k][j][i-1].x + coor[k][j-1][i-1].x + coor[k-1][j-1][i-1].x + coor[k-1][j][i-1].x);
                ycm = 0.25 * (coor[k][j][i-1].y + coor[k][j-1][i-1].y + coor[k-1][j-1][i-1].y + coor[k-1][j][i-1].y);
                zcm = 0.25 * (coor[k][j][i-1].z + coor[k][j-1][i-1].z + coor[k-1][j-1][i-1].z + coor[k-1][j][i-1].z);
                gs[k][j][i].x = PetscSqrtReal(PetscSqr(xcp-xcm) + PetscSqr(ycp-ycm) + PetscSqr(zcp-zcm));
       
                // Calculate Grid Spacing in j-direction (distance between j-face centers)
                xcp = 0.25 * (coor[k][j][i].x + coor[k][j][i-1].x + coor[k-1][j][i].x + coor[k-1][j][i-1].x);
                ycp = 0.25 * (coor[k][j][i].y + coor[k][j][i-1].y + coor[k-1][j][i].y + coor[k-1][j][i-1].y);
                zcp = 0.25 * (coor[k][j][i].z + coor[k][j][i-1].z + coor[k-1][j][i].z + coor[k-1][j][i-1].z);
                xcm = 0.25 * (coor[k][j-1][i].x + coor[k][j-1][i-1].x + coor[k-1][j-1][i].x + coor[k-1][j-1][i-1].x);
                ycm = 0.25 * (coor[k][j-1][i].y + coor[k][j-1][i-1].y + coor[k-1][j-1][i].y + coor[k-1][j-1][i-1].y);
                zcm = 0.25 * (coor[k][j-1][i].z + coor[k][j-1][i-1].z + coor[k-1][j-1][i].z + coor[k-1][j-1][i-1].z);
                gs[k][j][i].y = PetscSqrtReal(PetscSqr(xcp-xcm) + PetscSqr(ycp-ycm) + PetscSqr(zcp-zcm));

                // Calculate Grid Spacing in k-direction (distance between k-face centers)
                xcp = 0.25 * (coor[k][j][i].x + coor[k][j][i-1].x + coor[k][j-1][i].x + coor[k][j-1][i-1].x);
                ycp = 0.25 * (coor[k][j][i].y + coor[k][j][i-1].y + coor[k][j-1][i].y + coor[k][j-1][i-1].y);
                zcp = 0.25 * (coor[k][j][i].z + coor[k][j][i-1].z + coor[k][j-1][i].z + coor[k][j-1][i-1].z);
                xcm = 0.25 * (coor[k-1][j][i].x + coor[k-1][j][i-1].x + coor[k-1][j-1][i].x + coor[k-1][j-1][i-1].x);
                ycm = 0.25 * (coor[k-1][j][i].y + coor[k-1][j][i-1].y + coor[k-1][j-1][i].y + coor[k-1][j-1][i-1].y);
                zcm = 0.25 * (coor[k-1][j][i].z + coor[k-1][j-1][i-1].z + coor[k-1][j-1][i].z + coor[k-1][j-1][i-1].z);
                gs[k][j][i].z = PetscSqrtReal(PetscSqr(xcp-xcm) + PetscSqr(ycp-ycm) + PetscSqr(zcp-zcm));
            }
        }
    }
    
    ierr = DMDAVecRestoreArrayRead(user->fda, lCoords, &coor); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Cent, &cent); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->GridSpace, &gs); CHKERRQ(ierr);

    // Assemble and update ghost regions for the new data
    ierr = VecAssemblyBegin(user->Cent); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->Cent); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->GridSpace); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->GridSpace); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_CENT); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_GRID_SPACE); CHKERRQ(ierr);

    ierr = ApplyPeriodicCorrectionsToCellCentersAndSpacing(user); CHKERRQ(ierr);

    // Final assembly and ghost update after corrections
    ierr = VecAssemblyBegin(user->Cent); CHKERRQ(ierr); 
    ierr = VecAssemblyEnd(user->Cent); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_CENT); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `ComputeIFaceMetrics()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ComputeIFaceMetrics(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    Vec            lCoords;
    const Cmpnts ***coor;
    Cmpnts       ***centx; //***gs;
    const Cmpnts ***centx_const;
    Cmpnts       ***icsi, ***ieta, ***izet;
    PetscScalar  ***iaj;
    PetscReal      dxdc, dydc, dzdc, dxde, dyde, dzde, dxdz, dydz, dzdz;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Computing i-face metrics for level %d block %d...\n", user->simCtx->rank, user->thislevel, user->_this);

    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    PetscInt xs = info.xs, xe = info.xs + info.xm, mx = info.mx;
    PetscInt ys = info.ys, ye = info.ys + info.ym, my = info.my;
    PetscInt zs = info.zs, ze = info.zs + info.zm, mz = info.mz;
    PetscInt lxe = xe;
    PetscInt lys = ys; PetscInt lye = ye;
    PetscInt lzs = zs; PetscInt lze = ze;

    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;

    if (xe==mx) lxe=xe-1;
    if (ye==my) lye=ye-1;
    if (ze==mz) lze=ze-1;

    // --- Part 1: Calculate the location of i-face centers (Centx) ---
    ierr = DMGetCoordinatesLocal(user->da, &lCoords); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, lCoords, &coor); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Centx, &centx); CHKERRQ(ierr);
    //  ierr = DMDAVecGetArray(user->fda, user->lGridSpace,&gs); CHKERRQ(ierr);

    //LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d:   Calculating i-face centers (Centx) with i[%d,%d], j[%d,%d], k[%d,%d] ...\n", user->simCtx->rank,gxs,gxe,gys,gye,gzs,gze);

    // Populate only owned physical face centers. Periodic endpoint and ghost
    // coordinates are established by the canonical face-field synchronizer.
    for (PetscInt k = PetscMax(zs, 1); k < PetscMin(ze, mz - 1); k++) {
        for (PetscInt j = PetscMax(ys, 1); j < PetscMin(ye, my - 1); j++) {
            for (PetscInt i = xs; i < PetscMin(xe, mx - 1); i++) {
                //----- DEBUG ------
                //LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d:     Calculating i-face center at (k=%d, j=%d, i=%d)\n", user->simCtx->rank, k, j, i);
                //LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d:       Using corner nodes: (%f,%f,%f), (%f,%f,%f), (%f,%f,%f), (%f,%f,%f)\n", user->simCtx->rank,
                //          coor[k][j][i].x, coor[k][j][i].y, coor[k][j][i].z,
                //          coor[k-1][j][i].x, coor[k-1][j][i].y, coor[k-1][j][i].z,
                //          coor[k][j-1][i].x, coor[k][j-1][i].y, coor[k][j-1][i].z,
                //          coor[k-1][j-1][i].x, coor[k-1][j-1][i].y, coor[k-1][j-1][i].z);

                centx[k][j][i].x = 0.25 * (coor[k][j][i].x + coor[k-1][j][i].x + coor[k][j-1][i].x + coor[k-1][j-1][i].x);
                centx[k][j][i].y = 0.25 * (coor[k][j][i].y + coor[k-1][j][i].y + coor[k][j-1][i].y + coor[k-1][j-1][i].y);
                centx[k][j][i].z = 0.25 * (coor[k][j][i].z + coor[k-1][j][i].z + coor[k][j-1][i].z + coor[k-1][j-1][i].z);
            
                //LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d:       Calculated i-face center: (%f,%f,%f)\n", user->simCtx->rank, centx[k][j][i].x, centx[k][j][i].y, centx[k][j][i].z);
            }
        }
    }

    //LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d: i-face center coordinates calculated. \n", user->simCtx->rank);
    /*
    if(xs==0){
      for(PetscInt k=gzs+1;k < gze; k++){
	for(PetscInt j=gys+1;j < gye; j++){
	  PetscInt i=0; 
	  centx[k][j][i-1].x=centx[k][j][i].x-gs[k][j][i-2].x;
	  centx[k][j][i-1].y=centx[k][j][i].y;
	  centx[k][j][i-1].z=centx[k][j][i].z;
	}
      }
    }
    if (xe==mx){
      for(PetscInt k=gzs+1; k<gze; k++) {
	for (PetscInt j=gys+1; j<gye;j++) {
	  PetscInt i=mx-1; 
	  centx[k][j][i].x=centx[k][j][i-1].x+gs[k][j][i+2].x;
	  centx[k][j][i].y=centx[k][j][i-1].y;
	  centx[k][j][i].z=centx[k][j][i-1].z;
	}
      }
    }
    */
    
    ierr = DMDAVecRestoreArrayRead(user->fda, lCoords, &coor); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Centx, &centx); CHKERRQ(ierr);

    // ierr = DMDAVecRestoreArray(user->fda, user->lGridSpace,&gs); CHKERRQ(ierr);

    {
        const FieldId face_centers[] = {FIELD_ID_CENTX};
        ierr = SynchronizePeriodicFaceFields(user, 'i', 1, face_centers); CHKERRQ(ierr);
    }

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d:   i-face centers (Centx) calculated and ghosts updated.\n", user->simCtx->rank);

    // --- Part 2: Calculate metrics using face-centered coordinates ---
    ierr = DMDAVecGetArrayRead(user->fda, user->lCentx, &centx_const); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->ICsi, &icsi); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->IEta, &ieta); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->IZet, &izet); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->IAj, &iaj); CHKERRQ(ierr);

    // Loop over the OWNED region where we will store the final metrics
    for (PetscInt k=lzs; k<lze; k++) {
        for (PetscInt j=lys; j<lye; j++) {
            for (PetscInt i=xs; i<lxe; i++) {

                // --- Stencil Logic for d/dcsi (derivative in i-direction) ---
                if(i == 0 && user->boundary_faces[BC_FACE_NEG_X].mathematical_type != PERIODIC){
                        // Forward difference at the domain's min-i boundary
                        dxdc = centx_const[k][j][i+1].x - centx_const[k][j][i].x;
                        dydc = centx_const[k][j][i+1].y - centx_const[k][j][i].y;
                        dzdc = centx_const[k][j][i+1].z - centx_const[k][j][i].z;
                } else if (i == mx - 2 && user->boundary_faces[BC_FACE_POS_X].mathematical_type != PERIODIC) { 
                    // Backward difference at the domain's max-i boundary
                    dxdc = centx_const[k][j][i].x - centx_const[k][j][i-1].x;
                    dydc = centx_const[k][j][i].y - centx_const[k][j][i-1].y;
                    dzdc = centx_const[k][j][i].z - centx_const[k][j][i-1].z;
                } else { // Central difference in the interior  (or if PERIODIC BCs)
                    dxdc = 0.5 * (centx_const[k][j][i+1].x - centx_const[k][j][i-1].x);
                    dydc = 0.5 * (centx_const[k][j][i+1].y - centx_const[k][j][i-1].y);
                    dzdc = 0.5 * (centx_const[k][j][i+1].z - centx_const[k][j][i-1].z);
                }

                // --- Stencil Logic for d/deta (derivative in j-direction) ---
                if (j == 1 && user->boundary_faces[BC_FACE_NEG_Y].mathematical_type != PERIODIC) { 
                    // Forward difference
                    dxde = centx_const[k][j+1][i].x - centx_const[k][j][i].x;
                    dyde = centx_const[k][j+1][i].y - centx_const[k][j][i].y;
                    dzde = centx_const[k][j+1][i].z - centx_const[k][j][i].z;
                } else if (j == my - 2 && user->boundary_faces[BC_FACE_POS_Y].mathematical_type != PERIODIC) { 
                    // Backward difference
                    dxde = centx_const[k][j][i].x - centx_const[k][j-1][i].x;
                    dyde = centx_const[k][j][i].y - centx_const[k][j-1][i].y;
                    dzde = centx_const[k][j][i].z - centx_const[k][j-1][i].z;
                } else { // Central difference (interior or PERIODIC)
                    dxde = 0.5 * (centx_const[k][j+1][i].x - centx_const[k][j-1][i].x);
                    dyde = 0.5 * (centx_const[k][j+1][i].y - centx_const[k][j-1][i].y);
                    dzde = 0.5 * (centx_const[k][j+1][i].z - centx_const[k][j-1][i].z);
                }

                // --- Stencil Logic for d/dzeta (derivative in k-direction) ---
                if (k == 1 && user->boundary_faces[BC_FACE_NEG_Z].mathematical_type != PERIODIC) { 
                    // Forward difference
                    dxdz = centx_const[k+1][j][i].x - centx_const[k][j][i].x;
                    dydz = centx_const[k+1][j][i].y - centx_const[k][j][i].y;
                    dzdz = centx_const[k+1][j][i].z - centx_const[k][j][i].z;
                } else if (k == mz - 2 && user->boundary_faces[BC_FACE_POS_Z].mathematical_type != PERIODIC) { 
                    // Backward difference
                    dxdz = centx_const[k][j][i].x - centx_const[k-1][j][i].x;
                    dydz = centx_const[k][j][i].y - centx_const[k-1][j][i].y;
                    dzdz = centx_const[k][j][i].z - centx_const[k-1][j][i].z;
                } else { // Central difference (Interior + PERIODIC)
                    dxdz = 0.5 * (centx_const[k+1][j][i].x - centx_const[k-1][j][i].x);
                    dydz = 0.5 * (centx_const[k+1][j][i].y - centx_const[k-1][j][i].y);
                    dzdz = 0.5 * (centx_const[k+1][j][i].z - centx_const[k-1][j][i].z);
                }

                // --- Metric calculations (identical to legacy FormMetrics) ---
                icsi[k][j][i].x = dyde * dzdz - dzde * dydz;
                icsi[k][j][i].y = -dxde * dzdz + dzde * dxdz;
                icsi[k][j][i].z = dxde * dydz - dyde * dxdz;

                ieta[k][j][i].x = dydz * dzdc - dzdz * dydc;
                ieta[k][j][i].y = -dxdz * dzdc + dzdz * dxdc;
                ieta[k][j][i].z = dxdz * dydc - dydz * dxdc;

                izet[k][j][i].x = dydc * dzde - dzdc * dyde;
                izet[k][j][i].y = -dxdc * dzde + dzdc * dxde;
                izet[k][j][i].z = dxdc * dyde - dydc * dxde;

                iaj[k][j][i] = dxdc * icsi[k][j][i].x + dydc * icsi[k][j][i].y + dzdc * icsi[k][j][i].z;
                if (PetscAbsScalar(iaj[k][j][i]) > 1e-12) {
                    iaj[k][j][i] = 1.0 / iaj[k][j][i];
                }
            }
        }
    }

    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCentx, &centx_const); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->ICsi, &icsi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->IEta, &ieta); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->IZet, &izet); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, user->IAj, &iaj); CHKERRQ(ierr);
    
    // --- Part 3: Assemble global vectors and update local ghosts ---
    ierr = VecAssemblyBegin(user->ICsi); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->ICsi); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->IEta); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->IEta); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->IZet); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->IZet); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->IAj); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->IAj); CHKERRQ(ierr);

    ierr = UpdateLocalGhosts(user, FIELD_ID_ICSI); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_IETA); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_IZET); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_IAJ); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `ComputeJFaceMetrics()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ComputeJFaceMetrics(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    Vec            lCoords;
    const Cmpnts ***coor;
    Cmpnts       ***centy; //***gs;
    const Cmpnts ***centy_const;
    Cmpnts       ***jcsi, ***jeta, ***jzet;
    PetscScalar  ***jaj;
    PetscReal      dxdc, dydc, dzdc, dxde, dyde, dzde, dxdz, dydz, dzdz;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Computing j-face metrics for level %d block %d...\n", user->simCtx->rank, user->thislevel, user->_this);

    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    PetscInt xs = info.xs, xe = info.xs + info.xm, mx = info.mx;
    PetscInt ys = info.ys, ye = info.ys + info.ym, my = info.my;
    PetscInt zs = info.zs, ze = info.zs + info.zm, mz = info.mz;
    PetscInt lxs = xs; PetscInt lxe = xe;
    PetscInt lye = ye;
    PetscInt lzs = zs; PetscInt lze = ze;

    if (xs==0) lxs = xs+1;
    if (zs==0) lzs = zs+1;

    if (xe==mx) lxe=xe-1;
    if (ye==my) lye=ye-1;
    if (ze==mz) lze=ze-1;

    // --- Part 1: Calculate the location of i-face centers (Centx) ---
    ierr = DMGetCoordinatesLocal(user->da, &lCoords); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, lCoords, &coor); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Centy, &centy); CHKERRQ(ierr);
    //  ierr = DMDAVecGetArray(user->fda, user->lGridSpace,&gs); CHKERRQ(ierr);

    for (PetscInt k = PetscMax(zs, 1); k < PetscMin(ze, mz - 1); k++) {
        for (PetscInt j = ys; j < PetscMin(ye, my - 1); j++) {
            for (PetscInt i = PetscMax(xs, 1); i < PetscMin(xe, mx - 1); i++) {
                centy[k][j][i].x = 0.25 * (coor[k][j][i].x + coor[k-1][j][i].x + coor[k][j][i-1].x + coor[k-1][j][i-1].x);
                centy[k][j][i].y = 0.25 * (coor[k][j][i].y + coor[k-1][j][i].y + coor[k][j][i-1].y + coor[k-1][j][i-1].y);
                centy[k][j][i].z = 0.25 * (coor[k][j][i].z + coor[k-1][j][i].z + coor[k][j][i-1].z + coor[k-1][j][i-1].z);
            }
        }
    }

    /*
    if(ys==0){
      for(PetscInt k=gzs+1;k < gze; k++){
	for(PetscInt i=gxs+1;j < gxe; i++){
	  PetscInt j=0; 
	  centy[k][j-1][i].x=centy[k][j][i].x;
	  centy[k][j-1][i].y=centy[k][j][i].y-gs[k][j-2][i].y;
	  centy[k][j-1][i].z=centy[k][j][i].z;
	}
      }
    }
    if (ye==my){
      for(PetscInt k=gzs+1; k<gze; k++) {
	for (PetscInt i=gxs+1; j<gxe;i++) {
	  PetscInt j=my-1; 
	  centy[k][j][i].x=centy[k][j-1][i].x
	  centy[k][j][i].y=centy[k][j-1][i].y+gs[k][j+2][i].y;
	  centy[k][j][i].z=centy[k][j-1][i].z;
	}
      }
    }
    */
    
    ierr = DMDAVecRestoreArrayRead(user->fda, lCoords, &coor); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Centy, &centy); CHKERRQ(ierr);
    // ierr = DMDAVecRestoreArray(user->fda, user->lGridSpace,&gs); CHKERRQ(ierr);

    {
        const FieldId face_centers[] = {FIELD_ID_CENTY};
        ierr = SynchronizePeriodicFaceFields(user, 'j', 1, face_centers); CHKERRQ(ierr);
    }

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d:   j-face centers (Centx) calculated and ghosts updated.\n", user->simCtx->rank);

    // --- Part 2: Calculate metrics using face-centered coordinates ---
    ierr = DMDAVecGetArrayRead(user->fda, user->lCenty, &centy_const); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->JCsi, &jcsi); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->JEta, &jeta); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->JZet, &jzet); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->JAj, &jaj); CHKERRQ(ierr);

    // Loop over the OWNED region where we will store the final metrics
    for (PetscInt k=lzs; k<lze; k++) {
        for (PetscInt j=ys; j<lye; j++) {
            for (PetscInt i=lxs; i<lxe; i++) {

                // --- Stencil Logic for d/dcsi (derivative in i-direction) ---
                if (i == 1 && user->boundary_faces[BC_FACE_NEG_X].mathematical_type != PERIODIC) { 
                    // Forward difference at the domain's min-i boundary
                    dxdc = centy_const[k][j][i+1].x - centy_const[k][j][i].x;
                    dydc = centy_const[k][j][i+1].y - centy_const[k][j][i].y;
                    dzdc = centy_const[k][j][i+1].z - centy_const[k][j][i].z;
                } else if (i == mx - 2 && user->boundary_faces[BC_FACE_POS_X].mathematical_type != PERIODIC) { 
                    // Backward difference at the domain's max-i boundary
                    dxdc = centy_const[k][j][i].x - centy_const[k][j][i-1].x;
                    dydc = centy_const[k][j][i].y - centy_const[k][j][i-1].y;
                    dzdc = centy_const[k][j][i].z - centy_const[k][j][i-1].z;
                } else { // Central difference in the interior or PERIODIC
                    dxdc = 0.5 * (centy_const[k][j][i+1].x - centy_const[k][j][i-1].x);
                    dydc = 0.5 * (centy_const[k][j][i+1].y - centy_const[k][j][i-1].y);
                    dzdc = 0.5 * (centy_const[k][j][i+1].z - centy_const[k][j][i-1].z);
                }

                // --- Stencil Logic for d/deta (derivative in j-direction) ---
                if (j == 0 && user->boundary_faces[BC_FACE_NEG_Y].mathematical_type != PERIODIC) { 
                    // Forward difference
                    dxde = centy_const[k][j+1][i].x - centy_const[k][j][i].x;
                    dyde = centy_const[k][j+1][i].y - centy_const[k][j][i].y;
                    dzde = centy_const[k][j+1][i].z - centy_const[k][j][i].z;
                } else if (j == my - 2 && user->boundary_faces[BC_FACE_POS_Y].mathematical_type != PERIODIC) {
                    // Backward difference
                    dxde = centy_const[k][j][i].x - centy_const[k][j-1][i].x;
                    dyde = centy_const[k][j][i].y - centy_const[k][j-1][i].y;
                    dzde = centy_const[k][j][i].z - centy_const[k][j-1][i].z;
                } else { // Central difference (interior or PERIODIC)
                    dxde = 0.5 * (centy_const[k][j+1][i].x - centy_const[k][j-1][i].x);
                    dyde = 0.5 * (centy_const[k][j+1][i].y - centy_const[k][j-1][i].y);
                    dzde = 0.5 * (centy_const[k][j+1][i].z - centy_const[k][j-1][i].z);
                }

                // --- Stencil Logic for d/dzeta (derivative in k-direction) ---
                if (k == 1 && user->boundary_faces[BC_FACE_NEG_Z].mathematical_type != PERIODIC) { 
                    // Forward difference
                    dxdz = centy_const[k+1][j][i].x - centy_const[k][j][i].x;
                    dydz = centy_const[k+1][j][i].y - centy_const[k][j][i].y;
                    dzdz = centy_const[k+1][j][i].z - centy_const[k][j][i].z;
                } else if (k == mz - 2 && user->boundary_faces[BC_FACE_POS_Z].mathematical_type != PERIODIC) { 
                    // Backward difference
                    dxdz = centy_const[k][j][i].x - centy_const[k-1][j][i].x;
                    dydz = centy_const[k][j][i].y - centy_const[k-1][j][i].y;
                    dzdz = centy_const[k][j][i].z - centy_const[k-1][j][i].z;
                } else { // Central difference (Interior or PERIODIC)
                    dxdz = 0.5 * (centy_const[k+1][j][i].x - centy_const[k-1][j][i].x);
                    dydz = 0.5 * (centy_const[k+1][j][i].y - centy_const[k-1][j][i].y);
                    dzdz = 0.5 * (centy_const[k+1][j][i].z - centy_const[k-1][j][i].z);
                }

                // --- Metric calculations (identical to legacy FormMetrics) ---
                jcsi[k][j][i].x = dyde * dzdz - dzde * dydz;
                jcsi[k][j][i].y = -dxde * dzdz + dzde * dxdz;
                jcsi[k][j][i].z = dxde * dydz - dyde * dxdz;

                jeta[k][j][i].x = dydz * dzdc - dzdz * dydc;
                jeta[k][j][i].y = -dxdz * dzdc + dzdz * dxdc;
                jeta[k][j][i].z = dxdz * dydc - dydz * dxdc;

                jzet[k][j][i].x = dydc * dzde - dzdc * dyde;
                jzet[k][j][i].y = -dxdc * dzde + dzdc * dxde;
                jzet[k][j][i].z = dxdc * dyde - dydc * dxde;

                jaj[k][j][i] = dxdc * jcsi[k][j][i].x + dydc * jcsi[k][j][i].y + dzdc * jcsi[k][j][i].z;
                if (PetscAbsScalar(jaj[k][j][i]) > 1e-12) {
                    jaj[k][j][i] = 1.0 / jaj[k][j][i];
                }
            }
        }
    }

    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCenty, &centy_const); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->JCsi, &jcsi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->JEta, &jeta); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->JZet, &jzet); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, user->JAj, &jaj); CHKERRQ(ierr);
    
    // --- Part 3: Assemble global vectors and update local ghosts ---
    ierr = VecAssemblyBegin(user->JCsi); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->JCsi); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->JEta); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->JEta); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->JZet); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->JZet); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->JAj); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->JAj); CHKERRQ(ierr);

    ierr = UpdateLocalGhosts(user, FIELD_ID_JCSI); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_JETA); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_JZET); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_JAJ); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `ComputeKFaceMetrics()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ComputeKFaceMetrics(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    Vec            lCoords;
    const Cmpnts ***coor;
    Cmpnts       ***centz; //***gs;
    const Cmpnts ***centz_const;
    Cmpnts       ***kcsi, ***keta, ***kzet;
    PetscScalar  ***kaj;
    PetscReal      dxdc, dydc, dzdc, dxde, dyde, dzde, dxdz, dydz, dzdz;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Computing k-face metrics for level %d block %d...\n", user->simCtx->rank, user->thislevel, user->_this);

    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    PetscInt xs = info.xs, xe = info.xs + info.xm, mx = info.mx;
    PetscInt ys = info.ys, ye = info.ys + info.ym, my = info.my;
    PetscInt zs = info.zs, ze = info.zs + info.zm, mz = info.mz;
    PetscInt lxs = xs; PetscInt lxe = xe;
    PetscInt lys = ys; PetscInt lye = ye;
    PetscInt lze = ze;

    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;

    if (xe==mx) lxe=xe-1;
    if (ye==my) lye=ye-1;
    if (ze==mz) lze=ze-1;

    // --- Part 1: Calculate the location of i-face centers (Centx) ---
    ierr = DMGetCoordinatesLocal(user->da, &lCoords); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, lCoords, &coor); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Centz, &centz); CHKERRQ(ierr);
    //  ierr = DMDAVecGetArray(user->fda, user->lGridSpace,&gs); CHKERRQ(ierr);

    for (PetscInt k = zs; k < PetscMin(ze, mz - 1); k++) {
        for (PetscInt j = PetscMax(ys, 1); j < PetscMin(ye, my - 1); j++) {
            for (PetscInt i = PetscMax(xs, 1); i < PetscMin(xe, mx - 1); i++) {
                centz[k][j][i].x = 0.25 * (coor[k][j][i].x + coor[k][j-1][i].x + coor[k][j][i-1].x + coor[k][j-1][i-1].x);
                centz[k][j][i].y = 0.25 * (coor[k][j][i].y + coor[k][j-1][i].y + coor[k][j][i-1].y + coor[k][j-1][i-1].y);
                centz[k][j][i].z = 0.25 * (coor[k][j][i].z + coor[k][j-1][i].z + coor[k][j][i-1].z + coor[k][j-1][i-1].z);
            }
        }
    }

    /*
    if(zs==0){
      for(PetscInt j=gys+1;j < gye; j++){
	for(PetscInt i=gxs+1;j < gxe; i++){
	  PetscInt k=0; 
	  centz[k-1][j][i].x=centz[k][j][i].x;
	  centz[k-1][j][i].y=centz[k][j][i].y;
	  centz[k-1][j][i].z=centz[k][j][i].z-gs[k-2][j][i].z;
	}
      }
    }
    if (ze==mz){
      for(PetscInt j=gys+1; j<gye; j++) {
	for (PetscInt i=gxs+1; j<gxe;i++) {
	  PetscInt k=mz-1; 
	  centy[k][j][i].x=centy[k-1][j][i].x
	  centy[k][j][i].y=centy[k-1][j][i].y;
	  centz[k][j][i].z=centz[k-1][j][i].z+gs[k+2][j][1].z;
	}
      }
    }
    */
    
    ierr = DMDAVecRestoreArrayRead(user->fda, lCoords, &coor); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->Centz, &centz); CHKERRQ(ierr);
    // ierr = DMDAVecRestoreArray(user->fda, user->lGridSpace,&gs); CHKERRQ(ierr);

    {
        const FieldId face_centers[] = {FIELD_ID_CENTZ};
        ierr = SynchronizePeriodicFaceFields(user, 'k', 1, face_centers); CHKERRQ(ierr);
    }

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d:   k-face centers (Centx) calculated and ghosts updated.\n", user->simCtx->rank);

    // --- Part 2: Calculate metrics using face-centered coordinates ---
    ierr = DMDAVecGetArrayRead(user->fda, user->lCentz, &centz_const); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->KCsi, &kcsi); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->KEta, &keta); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->KZet, &kzet); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->KAj, &kaj); CHKERRQ(ierr);

    // Loop over the OWNED region where we will store the final metrics
    for (PetscInt k=zs; k<lze; k++) {
        for (PetscInt j=lys; j<lye; j++) {
            for (PetscInt i=lxs; i<lxe; i++) {

                // --- Stencil Logic for d/dcsi (derivative in i-direction) ---
                if (i == 1 && user->boundary_faces[BC_FACE_NEG_X].mathematical_type != PERIODIC) { 
                    // Forward difference at the domain's min-i boundary
                    dxdc = centz_const[k][j][i+1].x - centz_const[k][j][i].x;
                    dydc = centz_const[k][j][i+1].y - centz_const[k][j][i].y;
                    dzdc = centz_const[k][j][i+1].z - centz_const[k][j][i].z;
                } else if (i == mx - 2 && user->boundary_faces[BC_FACE_POS_X].mathematical_type != PERIODIC) { 
                    // Backward difference at the domain's max-i boundary
                    dxdc = centz_const[k][j][i].x - centz_const[k][j][i-1].x;
                    dydc = centz_const[k][j][i].y - centz_const[k][j][i-1].y;
                    dzdc = centz_const[k][j][i].z - centz_const[k][j][i-1].z;
                } else { // Central difference in the interior (or PERIODIC)
                    dxdc = 0.5 * (centz_const[k][j][i+1].x - centz_const[k][j][i-1].x);
                    dydc = 0.5 * (centz_const[k][j][i+1].y - centz_const[k][j][i-1].y);
                    dzdc = 0.5 * (centz_const[k][j][i+1].z - centz_const[k][j][i-1].z);
                }

                // --- Stencil Logic for d/deta (derivative in j-direction) ---
                if (j == 1 && user->boundary_faces[BC_FACE_NEG_Y].mathematical_type != PERIODIC) { 
                    // Forward difference
                    dxde = centz_const[k][j+1][i].x - centz_const[k][j][i].x;
                    dyde = centz_const[k][j+1][i].y - centz_const[k][j][i].y;
                    dzde = centz_const[k][j+1][i].z - centz_const[k][j][i].z;
                } else if (j == my - 2 && user->boundary_faces[BC_FACE_POS_Y].mathematical_type != PERIODIC) { 
                    // Backward difference
                    dxde = centz_const[k][j][i].x - centz_const[k][j-1][i].x;
                    dyde = centz_const[k][j][i].y - centz_const[k][j-1][i].y;
                    dzde = centz_const[k][j][i].z - centz_const[k][j-1][i].z;
                } else { // Central difference (interior or PERIODIC)
                    dxde = 0.5 * (centz_const[k][j+1][i].x - centz_const[k][j-1][i].x);
                    dyde = 0.5 * (centz_const[k][j+1][i].y - centz_const[k][j-1][i].y);
                    dzde = 0.5 * (centz_const[k][j+1][i].z - centz_const[k][j-1][i].z);
                }

                // --- Stencil Logic for d/dzeta (derivative in k-direction) ---
                if (k == 0 && user->boundary_faces[BC_FACE_NEG_Z].mathematical_type != PERIODIC) { 
                    // Forward difference
                    dxdz = centz_const[k+1][j][i].x - centz_const[k][j][i].x;
                    dydz = centz_const[k+1][j][i].y - centz_const[k][j][i].y;
                    dzdz = centz_const[k+1][j][i].z - centz_const[k][j][i].z;
                } else if (k == mz - 2 && user->boundary_faces[BC_FACE_POS_Z].mathematical_type != PERIODIC) { 
                    // Backward difference
                    dxdz = centz_const[k][j][i].x - centz_const[k-1][j][i].x;
                    dydz = centz_const[k][j][i].y - centz_const[k-1][j][i].y;
                    dzdz = centz_const[k][j][i].z - centz_const[k-1][j][i].z;
                } else { // Central difference (Interior or PERIODIC)
                    dxdz = 0.5 * (centz_const[k+1][j][i].x - centz_const[k-1][j][i].x);
                    dydz = 0.5 * (centz_const[k+1][j][i].y - centz_const[k-1][j][i].y);
                    dzdz = 0.5 * (centz_const[k+1][j][i].z - centz_const[k-1][j][i].z);
                }

                // --- Metric calculations (identical to legacy FormMetrics) ---
                kcsi[k][j][i].x = dyde * dzdz - dzde * dydz;
                kcsi[k][j][i].y = -dxde * dzdz + dzde * dxdz;
                kcsi[k][j][i].z = dxde * dydz - dyde * dxdz;

                keta[k][j][i].x = dydz * dzdc - dzdz * dydc;
                keta[k][j][i].y = -dxdz * dzdc + dzdz * dxdc;
                keta[k][j][i].z = dxdz * dydc - dydz * dxdc;

                kzet[k][j][i].x = dydc * dzde - dzdc * dyde;
                kzet[k][j][i].y = -dxdc * dzde + dzdc * dxde;
                kzet[k][j][i].z = dxdc * dyde - dydc * dxde;

                kaj[k][j][i] = dxdc * kcsi[k][j][i].x + dydc * kcsi[k][j][i].y + dzdc * kcsi[k][j][i].z;
                if (PetscAbsScalar(kaj[k][j][i]) > 1e-12) {
                    kaj[k][j][i] = 1.0 / kaj[k][j][i];
                }
            }
        }
    }

    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCentz, &centz_const); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->KCsi, &kcsi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->KEta, &keta); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->KZet, &kzet); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, user->KAj, &kaj); CHKERRQ(ierr);
    
    // --- Part 3: Assemble global vectors and update local ghosts ---
    ierr = VecAssemblyBegin(user->KCsi); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->KCsi); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->KEta); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->KEta); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->KZet); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->KZet); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->KAj); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->KAj); CHKERRQ(ierr);

    ierr = UpdateLocalGhosts(user, FIELD_ID_KCSI); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_KETA); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_KZET); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, FIELD_ID_KAJ); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

/**
 * @brief Convert logical cell indices into the flattened global cell identifier.
 */
static PetscInt Gidx(PetscInt i, PetscInt j, PetscInt k, UserCtx *user)
{
  PROFILE_HELPER_BEGIN;
  PetscInt nidx;
  DMDALocalInfo	info = user->info;

  PetscInt	mx = info.mx, my = info.my;
  
  AO ao;
  DMDAGetAO(user->da, &ao);
  nidx=i+j*mx+k*mx*my;
  
  AOApplicationToPetsc(ao,1,&nidx);
  
  PROFILE_HELPER_END;
  return (nidx);
}

/**
 * @brief Internal helper implementation: `ComputeMetricsDivergence()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ComputeMetricsDivergence(UserCtx *user)
{
  DM            da = user->da, fda = user->fda;
  DMDALocalInfo info = user->info;
  PetscInt      xs = info.xs, xe = info.xs + info.xm;
  PetscInt      ys = info.ys, ye = info.ys + info.ym;
  PetscInt      zs = info.zs, ze = info.zs + info.zm;
  PetscInt      mx = info.mx, my = info.my, mz = info.mz;
  PetscInt      lxs, lys, lzs, lxe, lye, lze;
  PetscInt      i, j, k;
  Vec           Div;
  PetscReal     ***div, ***aj;
  Cmpnts        ***csi, ***eta, ***zet;
  PetscReal     maxdiv;

  PetscFunctionBeginUser;

  PROFILE_FUNCTION_BEGIN;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs == 0) lxs = xs + 1;
  if (ys == 0) lys = ys + 1;
  if (zs == 0) lzs = zs + 1;

  if (xe == mx) lxe = xe - 1;
  if (ye == my) lye = ye - 1;
  if (ze == mz) lze = ze - 1;

  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
  DMDAVecGetArray(da, user->lAj, &aj);

  VecDuplicate(user->P, &Div);
  VecSet(Div, 0.);
  DMDAVecGetArray(da, Div, &div);

  for (k = lzs; k < lze; k++) {
    for (j = lys; j < lye; j++) {
      for (i = lxs; i < lxe; i++) {
        PetscReal divergence = (csi[k][j][i].x - csi[k][j][i-1].x +
                                eta[k][j][i].x - eta[k][j-1][i].x +
                                zet[k][j][i].x - zet[k-1][j][i].x +
                                csi[k][j][i].y - csi[k][j][i-1].y +
                                eta[k][j][i].y - eta[k][j-1][i].y +
                                zet[k][j][i].y - zet[k-1][j][i].y +
                                csi[k][j][i].z - csi[k][j][i-1].z +
                                eta[k][j][i].z - eta[k][j-1][i].z +
                                zet[k][j][i].z - zet[k-1][j][i].z) * aj[k][j][i];
        div[k][j][i] = fabs(divergence);
      }
    }
  }

  DMDAVecRestoreArray(da, Div, &div);

  PetscInt MaxFlatIndex = -1;
  VecMax(Div, &MaxFlatIndex, &maxdiv);
  LOG_ALLOW(GLOBAL,LOG_INFO,"The Maximum Metric Divergence is %e at flat index %" PetscInt_FMT ".\n",maxdiv,MaxFlatIndex);

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (Gidx(i,j,k,user) == MaxFlatIndex) {
	  LOG_ALLOW(GLOBAL,LOG_INFO,"The Maximum Metric Divergence(%e) is at location [%d][%d][%d]. \n", maxdiv,(int)k,(int)j,(int)i);
	}
      }
    }
  }
 
  
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  DMDAVecRestoreArray(da, user->lAj, &aj);
  VecDestroy(&Div);


  PROFILE_FUNCTION_END;

  PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `ComputeMetricNorms()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ComputeMetricNorms(UserCtx *user)
{
  
  DMDALocalInfo info = user->info;
  PetscInt      xs = info.xs, xe = info.xs + info.xm;
  PetscInt      ys = info.ys, ye = info.ys + info.ym;
  PetscInt      zs = info.zs, ze = info.zs + info.zm;
  PetscInt      i, j, k;
  
  PetscFunctionBeginUser;

  PROFILE_FUNCTION_BEGIN;

  PetscReal CsiMax, EtaMax, ZetMax;
  PetscReal ICsiMax, IEtaMax, IZetMax;
  PetscReal JCsiMax, JEtaMax, JZetMax;
  PetscReal KCsiMax, KEtaMax, KZetMax;
  PetscReal AjMax, IAjMax, JAjMax, KAjMax;

  PetscInt CsiMaxArg, EtaMaxArg, ZetMaxArg;
  PetscInt ICsiMaxArg, IEtaMaxArg, IZetMaxArg;
  PetscInt JCsiMaxArg, JEtaMaxArg, JZetMaxArg;
  PetscInt KCsiMaxArg, KEtaMaxArg, KZetMaxArg;
  PetscInt AjMaxArg, IAjMaxArg, JAjMaxArg, KAjMaxArg;  

  // Max Values
  VecMax(user->lCsi,&CsiMaxArg,&CsiMax);
  VecMax(user->lEta,&EtaMaxArg,&EtaMax);
  VecMax(user->lZet,&ZetMaxArg,&ZetMax);

  VecMax(user->lICsi,&ICsiMaxArg,&ICsiMax);
  VecMax(user->lIEta,&IEtaMaxArg,&IEtaMax);
  VecMax(user->lIZet,&IZetMaxArg,&IZetMax);

  VecMax(user->lJCsi,&JCsiMaxArg,&JCsiMax);
  VecMax(user->lJEta,&JEtaMaxArg,&JEtaMax);
  VecMax(user->lJZet,&JZetMaxArg,&JZetMax);

  VecMax(user->lKCsi,&KCsiMaxArg,&KCsiMax);
  VecMax(user->lKEta,&KEtaMaxArg,&KEtaMax);
  VecMax(user->lKZet,&KZetMaxArg,&KZetMax);

  VecMax(user->lAj,&AjMaxArg,&AjMax);
  VecMax(user->lIAj,&IAjMaxArg,&IAjMax);
  VecMax(user->lJAj,&JAjMaxArg,&JAjMax);
  VecMax(user->lKAj,&KAjMaxArg,&KAjMax);

  VecMax(user->lAj,&AjMaxArg,&AjMax);
  VecMax(user->lIAj,&IAjMaxArg,&IAjMax);
  VecMax(user->lJAj,&JAjMaxArg,&JAjMax);
  VecMax(user->lKAj,&KAjMaxArg,&KAjMax);

  LOG_ALLOW(GLOBAL,LOG_INFO," Metric Norms for MG level %d .\n",user->thislevel);
  
  LOG_ALLOW(GLOBAL,LOG_INFO,"The Max Metric Values are: CsiMax = %le, EtaMax = %le, ZetMax = %le.\n",CsiMax,EtaMax,ZetMax);
  LOG_ALLOW(GLOBAL,LOG_INFO,"The Max Metric Values are: ICsiMax = %le, IEtaMax = %le, IZetMax = %le.\n",ICsiMax,IEtaMax,IZetMax);
  LOG_ALLOW(GLOBAL,LOG_INFO,"The Max Metric Values are: JCsiMax = %le, JEtaMax = %le, JZetMax = %le.\n",JCsiMax,JEtaMax,JZetMax);
  LOG_ALLOW(GLOBAL,LOG_INFO,"The Max Metric Values are: KCsiMax = %le, KEtaMax = %le, KZetMax = %le.\n",KCsiMax,KEtaMax,KZetMax);
  LOG_ALLOW(GLOBAL,LOG_INFO,"The Max Volumes(Inverse) are: Aj = %le, IAj = %le, JAj = %le, KAj = %le.\n",AjMax,IAjMax,JAjMax,KAjMax);

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (Gidx(i,j,k,user) == CsiMaxArg) {
	  LOG_ALLOW(GLOBAL,LOG_INFO,"Max Csi = %le is at [%d][%d][%d] \n", CsiMax,k,j,i);
	}
	if (Gidx(i,j,k,user) == EtaMaxArg) {
	  LOG_ALLOW(GLOBAL,LOG_INFO,"Max Eta = %le is at [%d][%d][%d] \n", EtaMax,k,j,i);
	}
	if (Gidx(i,j,k,user) == ZetMaxArg) {
	  LOG_ALLOW(GLOBAL,LOG_INFO,"Max Zet = %le is at [%d][%d][%d] \n", ZetMax,k,j,i);
	}
	if (Gidx(i,j,k,user) == ICsiMaxArg) {
	  LOG_ALLOW(GLOBAL,LOG_INFO,"Max ICsi = %le is at [%d][%d][%d] \n", ICsiMax,k,j,i);
	}
	if (Gidx(i,j,k,user) == IEtaMaxArg) {
	  LOG_ALLOW(GLOBAL,LOG_INFO,"Max IEta = %le is at [%d][%d][%d] \n", IEtaMax,k,j,i);
	}
	if (Gidx(i,j,k,user) == IZetMaxArg) {
	  LOG_ALLOW(GLOBAL,LOG_INFO,"Max IZet = %le is at [%d][%d][%d] \n", IZetMax,k,j,i);
	}
	if (Gidx(i,j,k,user) == JCsiMaxArg) {
	  LOG_ALLOW(GLOBAL,LOG_INFO,"Max JCsi = %le is at [%d][%d][%d] \n", JCsiMax,k,j,i);
	}
	if (Gidx(i,j,k,user) == JEtaMaxArg) {
	  LOG_ALLOW(GLOBAL,LOG_INFO,"Max JEta = %le is at [%d][%d][%d] \n", JEtaMax,k,j,i);
	}
	if (Gidx(i,j,k,user) == JZetMaxArg) {
	  LOG_ALLOW(GLOBAL,LOG_INFO,"Max JZet = %le is at [%d][%d][%d] \n", JZetMax,k,j,i);
	}
	if (Gidx(i,j,k,user) == KCsiMaxArg) {
	  LOG_ALLOW(GLOBAL,LOG_INFO,"Max KCsi = %le is at [%d][%d][%d] \n", KCsiMax,k,j,i);
	}
	if (Gidx(i,j,k,user) == KEtaMaxArg) {
	  LOG_ALLOW(GLOBAL,LOG_INFO,"Max KEta = %le is at [%d][%d][%d] \n", KEtaMax,k,j,i);
	}
	if (Gidx(i,j,k,user) == KZetMaxArg) {
	  LOG_ALLOW(GLOBAL,LOG_INFO,"Max KZet = %le is at [%d][%d][%d] \n", KZetMax,k,j,i);
	}
	if (Gidx(i,j,k,user) == AjMaxArg) {
	  LOG_ALLOW(GLOBAL,LOG_INFO,"Max Aj = %le is at [%d][%d][%d] \n", AjMax,k,j,i);
	}
	if (Gidx(i,j,k,user) == IAjMaxArg) {
	  LOG_ALLOW(GLOBAL,LOG_INFO,"Max IAj = %le is at [%d][%d][%d] \n", IAjMax,k,j,i);
	}
	if (Gidx(i,j,k,user) == JAjMaxArg) {
	  LOG_ALLOW(GLOBAL,LOG_INFO,"Max JAj = %le is at [%d][%d][%d] \n", JAjMax,k,j,i);
	}
	if (Gidx(i,j,k,user) == KAjMaxArg) {
	  LOG_ALLOW(GLOBAL,LOG_INFO,"Max KAj = %le is at [%d][%d][%d] \n", KAjMax,k,j,i);
	}
      }
    }
  }

  /*
  VecView(user->lCsi,PETSC_VIEWER_STDOUT_WORLD);
  VecView(user->lEta,PETSC_VIEWER_STDOUT_WORLD);
  VecView(user->lZet,PETSC_VIEWER_STDOUT_WORLD);
  */
  
  PROFILE_FUNCTION_END;

  PetscFunctionReturn(0);
}

/**
 * @brief Internal helper implementation: `CalculateAllGridMetrics()`.
 * @details Local to this translation unit.
 */
PetscErrorCode CalculateAllGridMetrics(SimCtx *simCtx)
{
    PetscErrorCode ierr;
    UserMG         *usermg = &simCtx->usermg;
    MGCtx          *mgctx = usermg->mgctx;
    PetscInt       nblk = simCtx->block_number;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(GLOBAL, LOG_INFO, "Calculating grid metrics for all levels and blocks...\n");

    // Loop through all levels and all blocks
    for (PetscInt level = usermg->mglevels -1 ; level >=0; level--) {
        for (PetscInt bi = 0; bi < nblk; bi++) {
            UserCtx *user = &mgctx[level].user[bi];
            LOG_ALLOW_SYNC(LOCAL, LOG_DEBUG, "Rank %d: Calculating metrics for level %d, block %d\n", simCtx->rank, level, bi);

            // Call the modern, modular helper functions for each UserCtx.
            // These functions are self-contained and operate on the data within the provided context.
            ierr = ComputeFaceMetrics(user); CHKERRQ(ierr);
            ierr = ComputeCellCenteredJacobianInverse(user); CHKERRQ(ierr);
	        ierr = CheckAndFixGridOrientation(user); CHKERRQ(ierr);
	        ierr = ComputeCellCentersAndSpacing(user); CHKERRQ(ierr);
            ierr = ComputeIFaceMetrics(user); CHKERRQ(ierr);
            ierr = ComputeJFaceMetrics(user); CHKERRQ(ierr);
            ierr = ComputeKFaceMetrics(user); CHKERRQ(ierr); 

            // Apply Periodic Boundary Condition Adjustments if necessary
            ierr = ApplyMetricsPeriodicBCs(user); CHKERRQ(ierr);
            // Diagnostics
	    ierr = ComputeMetricNorms(user);
            if (level == usermg->mglevels - 1) {
                ierr = ComputeMetricsDivergence(user); CHKERRQ(ierr);
            }
        }
    }

    // Face planes for particle search and interpolation, from the final coordinates and lCent.
    if (simCtx->np > 0) {
        for (PetscInt bi = 0; bi < nblk; bi++) {
            UserCtx *user = &mgctx[usermg->mglevels - 1].user[bi];
            ierr = MetricComputeFacePlanes(user); CHKERRQ(ierr);
        }
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Grid metrics calculation complete.\n");

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------------- */
/* End of Metric.c */
