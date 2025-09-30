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


#undef __FUNCT__
#define __FUNCT__ "MetricGetCellVertices"
/* ------------------------------------------------------------------------- */
/**
 * @brief Extract the eight vertex coordinates of the hexahedral cell (i,j,k).
 *
 * Vertices are returned in the standard trilinear ordering: bit 0 → x-corner,
 * bit 1 → y-corner, bit 2 → z-corner.  (000 = origin of the cell, 111 = far
 * corner.)
 */
PetscErrorCode MetricGetCellVertices(UserCtx *user,
                                     const Cmpnts ***X,   /* coord array */
                                     PetscInt i,PetscInt j,PetscInt k,
                                     Cmpnts V[8])
{
  PetscFunctionBeginUser;
  for (PetscInt c = 0; c < 8; ++c) {
    PetscInt ii = i + ((c & 1) ? 1 : 0);
    PetscInt jj = j + ((c & 2) ? 1 : 0);
    PetscInt kk = k + ((c & 4) ? 1 : 0);
    LOG_LOOP_ALLOW(GLOBAL, LOG_DEBUG,i+j+k,10," ii: %d,jj:%d,kk:%d - Retrieved.\n",ii,jj,kk);
    V[c] = X[kk][jj][ii];
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "TrilinearBlend"

/* ------------------------------------------------------------------------- */
/**
 * @brief Map logical coordinates to physical space using trilinear basis.
 *
 * @param[in]   V   Array of the eight vertex coordinates (MetricGetCellVertices).
 * @param[in]   xi,eta,zta  Logical coordinates in [0,1].
 * @param[out]  Xp  Physical coordinate.
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


#undef _FUNCT__
#define __FUNCT__ "MetricLogicalToPhysical"
/* ------------------------------------------------------------------------- */
/**
 * @brief Public wrapper: map (cell index, ξ,η,ζ) to (x,y,z).
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
  
  PROFILE_FUNCTION_BEGIN;

  ierr = MetricGetCellVertices(user,X,i,j,k,V); CHKERRQ(ierr);
  TrilinearBlend(V,xi,eta,zta,Xp);

  PROFILE_FUNCTION_END;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MetricJacobian"
/* ------------------------------------------------------------------------- */
/**
 * @brief Compute Jacobian matrix and its determinant at (xi,eta,zta).
 *
 *        J = [ x_ξ  x_η  x_ζ ]
 *            [ y_ξ  y_η  y_ζ ]
 *            [ z_ξ  z_η  z_ζ ]
 *
 * This is handy for converting physical velocities (u,v,w) into contravariant
 * components and for volume weighting.
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

  PROFILE_FUNCTION_BEGIN;

  ierr = MetricGetCellVertices(user,X,i,j,k,V); CHKERRQ(ierr);

  /* derivatives of trilinear shape functions */
  PetscReal dN_dXi[8], dN_dEta[8], dN_dZta[8];
  for (PetscInt c=0;c<8;++c) {
    PetscReal sx = (c & 1) ?  1.0 : -1.0;
    PetscReal sy = (c & 2) ?  1.0 : -1.0;
    PetscReal sz = (c & 4) ?  1.0 : -1.0;
    dN_dXi [c] = 0.125 * sx * ( (c&2?eta:1-eta) ) * ( (c&4?zta:1-zta) );
    dN_dEta[c] = 0.125 * sy * ( (c&1?xi :1-xi ) ) * ( (c&4?zta:1-zta) );
    dN_dZta[c] = 0.125 * sz * ( (c&1?xi :1-xi ) ) * ( (c&2?eta:1-eta) );
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


#undef __FUNCT__
#define __FUNCT__ "MetricVelocityContravariant"
/* ------------------------------------------------------------------------- */
/**
 * @brief Convert physical velocity (u,v,w) to contravariant components (u^xi, u^eta, u^zta).
 */
PetscErrorCode MetricVelocityContravariant(const PetscReal J[3][3], PetscReal detJ,
                                           const PetscReal u[3], PetscReal uc[3])
{
  PetscFunctionBeginUser;

  PROFILE_FUNCTION_BEGIN;

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

#undef __FUNCT__
#define __FUNCT__ "CheckAndFixGridOrientation"
/* -------------------------------------------------------------------------- */
/**
 * @brief Ensure a **right-handed** metric basis (`Csi`, `Eta`, `Zet`) and a
 *        **positive Jacobian** (`Aj`) over the whole domain.
 *
 * The metric-generation kernels are completely algebraic, so they will happily
 * deliver a *left-handed* basis if the mesh file enumerates nodes in the
 * opposite ζ-direction.  
 * This routine makes the orientation explicit and—if needed—repairs it
 * **once per run**:
 *
 * | Step | Action |
 * |------|--------|
 * | 1 | Compute global `Aj_min`, `Aj_max`.                          |
 * | 2 | **Mixed signs** (`Aj_min < 0 && Aj_max > 0`) &rarr; abort: the mesh is topologically inconsistent. |
 * | 3 | **All negative** (`Aj_max < 0`) &rarr; flip <br>`Csi`, `Eta`, `Zet`, `Aj` & update local ghosts. |
 * | 4 | Store `user->orientation = ±1` so BC / IC routines can apply sign-aware logic if they care about inlet direction. |
 *
 * @param[in,out] user  Fully initialised #UserCtx that already contains  
 *                      `Csi`, `Eta`, `Zet`, `Aj`, their **local** ghosts, and
 *                      valid distributed DMs.
 *
 * @return `0` on success or a PETSc error code on failure.
 *
 * @note  Call **immediately after** `ComputeCellCenteredJacobianInverse()` and
 *        before any routine that differentiates or applies BCs.
 *
 * @author vishal kandala
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

    /* Default: grid is right-handed unless proven otherwise */
    PetscInt orientation = +1;

    /* ---------------- step 3: repair left-handed mesh -------------- */
    if (aj_max < 0.0) {                  /* entire domain has Aj < 0   */
        orientation = -1;

        if (!rank)
            LOG_ALLOW(LOCAL, LOG_INFO,
                "[orientation] Detected left-handed grid – flipping metric vectors\n");

        /* Flip sign of *all* metric vectors and Aj                     */
        ierr = VecScale(user->Csi, -1.0); CHKERRQ(ierr);
        ierr = VecScale(user->Eta, -1.0); CHKERRQ(ierr);
        ierr = VecScale(user->Zet, -1.0); CHKERRQ(ierr);
        ierr = VecScale(user->Aj , -1.0); CHKERRQ(ierr);

        /* Local ghost regions now stale – refresh                      */
        ierr = UpdateLocalGhosts(user, "Csi"); CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, "Eta"); CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, "Zet"); CHKERRQ(ierr);
        ierr = UpdateLocalGhosts(user, "Aj");  CHKERRQ(ierr);

        /* Sanity print: Aj must be > 0 now                             */
        ierr = VecMin(user->Aj, NULL, &aj_min); CHKERRQ(ierr);
        ierr = VecMax(user->Aj, NULL, &aj_max); CHKERRQ(ierr);

        if (aj_min <= 0.0)
            SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER,
                "Failed to flip grid orientation – Aj still non-positive.");
	else if (aj_min && aj_max > 0.0)
	  orientation = +1;
    }

    /* ---------------- step 4: store result in UserCtx -------------- */
    user->GridOrientation = orientation;

    if (!rank)
        LOG_ALLOW(LOCAL, LOG_INFO,
            "[orientation] Grid confirmed %s-handed after flip (orientation=%+d)\n",
            (orientation>0) ? "right" : "left", orientation);

    PROFILE_FUNCTION_END;        

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeFaceMetrics"
/**
 * @brief Computes the primary face metric components (Csi, Eta, Zet), including
 *        boundary extrapolation, and stores them in the corresponding global Vec
 *        members of the UserCtx structure (user->Csi, user->Eta, user->Zet).
 *
 *        This is a self-contained routine that performs the following steps:
 *        1. Obtains local ghosted nodal coordinates using DMGetCoordinatesLocal.
 *        2. Calculates metrics for INTERIOR faces where finite difference stencils are valid.
 *        3. EXTRAPOLATES metrics for faces on the physical domain boundaries by copying
 *           from the nearest computed interior face.
 *        4. Assembles the global `user->Csi`, `user->Eta`, `user->Zet` Vecs.
 *        5. Updates the local ghosted `user->lCsi`, `user->lEta`, `user->lZet` Vecs.
 *
 * @param[in,out] user             Pointer to the UserCtx structure.
 *
 * @return PetscErrorCode 0 on success.
 *
 * @note
 *  - This function is a complete "compute and make ready" unit for Csi, Eta, and Zet.
 *  - It's recommended to call `VecZeroEntries` on user->Csi, Eta, Zet before this
 *    if they might contain old data.
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
    ierr = UpdateLocalGhosts(user, "Csi"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "Eta"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "Zet"); CHKERRQ(ierr);
    
    LOG_ALLOW(GLOBAL, LOG_INFO, "Completed calculation, extrapolation, and update for Csi, Eta, Zet.\n");

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ComputeCellCenteredJacobianInverse"
/**
 * @brief Calculates the cell-centered inverse Jacobian determinant (1/J), including
 *        boundary extrapolation, stores it in `user->Aj`, assembles `user->Aj`, and
 *        updates `user->lAj`.
 *
 * @param[in,out] user             Pointer to the UserCtx structure.
 *
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeCellCenteredJacobianInverse(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    PetscScalar  ***aj_arr;
    Cmpnts       ***nodal_coords_arr;
    Vec            localCoords_from_dm;

    PetscFunctionBeginUser;
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
    ierr = UpdateLocalGhosts(user, "Aj"); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL, LOG_INFO, "Completed calculation, extrapolation, and update for Aj.\n");
    PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ComputeCellCentersAndSpacing"
/**
 * @brief Computes the physical location of cell centers and the spacing between them.
 *
 * This function calculates two key geometric properties from the nodal coordinates:
 * 1.  `Cent`: A vector field storing the (x,y,z) coordinates of the center of each grid cell.
 * 2.  `GridSpace`: A vector field storing the physical distance between adjacent
 *     cell centers in the i, j, and k computational directions.
 *
 * It is a direct adaptation of the corresponding logic from the legacy `FormMetrics`.
 *
 * @param user The UserCtx for a specific grid level. The function populates `user->Cent` and `user->GridSpace`.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode ComputeCellCentersAndSpacing(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    Vec            lCoords;
    const Cmpnts ***coor;
    Cmpnts       ***cent, ***gs;
    PetscReal      xcp, ycp, zcp, xcm, ycm, zcm;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Computing cell centers and spacing for level %d block %d...\n", user->simCtx->rank, user->thislevel, user->_this);

    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocal(user->da, &lCoords); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, lCoords, &coor); CHKERRQ(ierr);

    ierr = DMDAVecGetArray(user->fda, user->Cent, &cent); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->GridSpace, &gs); CHKERRQ(ierr);
    
    // Loop over the interior OWNED cells (stencil requires i-1, j-1, k-1)
    for (PetscInt k=info.zs+1; k<info.zs+info.zm; k++) {
        for (PetscInt j=info.ys+1; j<info.ys+info.ym; j++) {
            for (PetscInt i=info.xs+1; i<info.xs+info.xm; i++) {
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
    ierr = UpdateLocalGhosts(user, "Cent"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "GridSpace"); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeIFaceMetrics"
/**
 * @brief Computes metrics centered on constant-i faces (i-faces).
 *
 * This function calculates the metric terms (`ICsi`, `IEta`, `IZet`) and the
 * inverse Jacobian (`IAj`) located at the geometric center of each constant-i
 * face. This is a critical step for staggered-grid finite difference schemes.
 *
 * The process is a direct and faithful refactoring of the corresponding logic
 * from the legacy `FormMetrics` function:
 * 1.  It first calculates the physical (x,y,z) coordinates of the center of
 *     each i-face and stores them in the `user->Centx` vector.
 * 2.  It then uses a boundary-aware, second-order finite difference stencil on
 *     the `Centx` field to compute the derivatives (e.g., d(x)/d(csi)).
 *     - Central differences are used in the grid interior.
 *     - One-sided differences are used at the physical domain boundaries.
 * 3.  Finally, these derivatives are used to compute the final metric terms and
 *     the inverse Jacobian, which are stored in their respective `Vec` objects.
 *
 * @param user The UserCtx for a specific grid level. This function populates
 *             the `user->ICsi`, `user->IEta`, `user->IZet`, and `user->IAj` vectors.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
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
    PetscInt gxs = info.gxs, gxe = info.gxs + info.gxm;
    PetscInt gys = info.gys, gye = info.gys + info.gym;
    PetscInt gzs = info.gzs, gze = info.gzs + info.gzm;
    
    PetscInt lxs = xs; PetscInt lxe = xe;
    PetscInt lys = ys; PetscInt lye = ye;
    PetscInt lzs = zs; PetscInt lze = ze;

    if (xs==0) lxs = xs+1;
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

    // Loop over the ghosted region to calculate all local face centers
    for (PetscInt k = gzs + 1; k < gze; k++) {
        for (PetscInt j = gys + 1; j < gye; j++) {
            for (PetscInt i = gxs; i < gxe; i++) {
                centx[k][j][i].x = 0.25 * (coor[k][j][i].x + coor[k-1][j][i].x + coor[k][j-1][i].x + coor[k-1][j-1][i].x);
                centx[k][j][i].y = 0.25 * (coor[k][j][i].y + coor[k-1][j][i].y + coor[k][j-1][i].y + coor[k-1][j-1][i].y);
                centx[k][j][i].z = 0.25 * (coor[k][j][i].z + coor[k-1][j][i].z + coor[k][j-1][i].z + coor[k-1][j-1][i].z);
            }
        }
    }

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

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d:   i-face centers (Centx) calculated and ghosts updated.\n", user->simCtx->rank);

    // --- Part 2: Calculate metrics using face-centered coordinates ---
    ierr = DMDAVecGetArrayRead(user->fda, user->Centx, &centx_const); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->ICsi, &icsi); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->IEta, &ieta); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->IZet, &izet); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->IAj, &iaj); CHKERRQ(ierr);

    // Loop over the OWNED region where we will store the final metrics
    for (PetscInt k=lzs; k<lze; k++) {
        for (PetscInt j=lys; j<lye; j++) {
            for (PetscInt i=xs; i<lxe; i++) {

                // --- Stencil Logic for d/dcsi (derivative in i-direction) ---
                if (i == 0) { // Forward difference at the domain's min-i boundary
                    dxdc = centx_const[k][j][i+1].x - centx_const[k][j][i].x;
                    dydc = centx_const[k][j][i+1].y - centx_const[k][j][i].y;
                    dzdc = centx_const[k][j][i+1].z - centx_const[k][j][i].z;
                } else if (i == mx - 2) { // Backward difference at the domain's max-i boundary
                    dxdc = centx_const[k][j][i].x - centx_const[k][j][i-1].x;
                    dydc = centx_const[k][j][i].y - centx_const[k][j][i-1].y;
                    dzdc = centx_const[k][j][i].z - centx_const[k][j][i-1].z;
                } else { // Central difference in the interior
                    dxdc = 0.5 * (centx_const[k][j][i+1].x - centx_const[k][j][i-1].x);
                    dydc = 0.5 * (centx_const[k][j][i+1].y - centx_const[k][j][i-1].y);
                    dzdc = 0.5 * (centx_const[k][j][i+1].z - centx_const[k][j][i-1].z);
                }

                // --- Stencil Logic for d/deta (derivative in j-direction) ---
                if (j == 1) { // Forward difference
                    dxde = centx_const[k][j+1][i].x - centx_const[k][j][i].x;
                    dyde = centx_const[k][j+1][i].y - centx_const[k][j][i].y;
                    dzde = centx_const[k][j+1][i].z - centx_const[k][j][i].z;
                } else if (j == my - 2) { // Backward difference
                    dxde = centx_const[k][j][i].x - centx_const[k][j-1][i].x;
                    dyde = centx_const[k][j][i].y - centx_const[k][j-1][i].y;
                    dzde = centx_const[k][j][i].z - centx_const[k][j-1][i].z;
                } else { // Central difference
                    dxde = 0.5 * (centx_const[k][j+1][i].x - centx_const[k][j-1][i].x);
                    dyde = 0.5 * (centx_const[k][j+1][i].y - centx_const[k][j-1][i].y);
                    dzde = 0.5 * (centx_const[k][j+1][i].z - centx_const[k][j-1][i].z);
                }

                // --- Stencil Logic for d/dzeta (derivative in k-direction) ---
                if (k == 1) { // Forward difference
                    dxdz = centx_const[k+1][j][i].x - centx_const[k][j][i].x;
                    dydz = centx_const[k+1][j][i].y - centx_const[k][j][i].y;
                    dzdz = centx_const[k+1][j][i].z - centx_const[k][j][i].z;
                } else if (k == mz - 2) { // Backward difference
                    dxdz = centx_const[k][j][i].x - centx_const[k-1][j][i].x;
                    dydz = centx_const[k][j][i].y - centx_const[k-1][j][i].y;
                    dzdz = centx_const[k][j][i].z - centx_const[k-1][j][i].z;
                } else { // Central difference
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

    ierr = DMDAVecRestoreArrayRead(user->fda, user->Centx, &centx_const); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->ICsi, &icsi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->IEta, &ieta); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->IZet, &izet); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, user->IAj, &iaj); CHKERRQ(ierr);
    
    // --- Part 3: Assemble global vectors and update local ghosts ---
    ierr = VecAssemblyBegin(user->ICsi); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->ICsi); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->IEta); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->IEta); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->IZet); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->IZet); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->IAj); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->IAj); CHKERRQ(ierr);

    ierr = UpdateLocalGhosts(user, "ICsi"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "IEta"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "IZet"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "IAj"); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeJFaceMetrics"
/**
 * @brief Computes metrics centered on constant-j faces (j-faces).
 *
 * This function calculates the metric terms (`JCsi`, `JEta`, `JZet`) and the
 * inverse Jacobian (`JAj`) located at the geometric center of each constant-j
 * face. This is a critical step for staggered-grid finite difference schemes.
 *
 * The process is a direct and faithful refactoring of the corresponding logic
 * from the legacy `FormMetrics` function:
 * 1.  It first calculates the physical (x,y,z) coordinates of the center of
 *     each i-face and stores them in the `user->Centy` vector.
 * 2.  It then uses a boundary-aware, second-order finite difference stencil on
 *     the `Centy` field to compute the derivatives (e.g., d(x)/d(csi)).
 *     - Central differences are used in the grid interior.
 *     - One-sided differences are used at the physical domain boundaries.
 * 3.  Finally, these derivatives are used to compute the final metric terms and
 *     the inverse Jacobian, which are stored in their respective `Vec` objects.
 *
 * @param user The UserCtx for a specific grid level. This function populates
 *             the `user->JCsi`, `user->JEta`, `user->JZet`, and `user->JAj` vectors.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
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
    PetscInt gxs = info.gxs, gxe = info.gxs + info.gxm;
    PetscInt gys = info.gys, gye = info.gys + info.gym;
    PetscInt gzs = info.gzs, gze = info.gzs + info.gzm;
    
    PetscInt lxs = xs; PetscInt lxe = xe;
    PetscInt lys = ys; PetscInt lye = ye;
    PetscInt lzs = zs; PetscInt lze = ze;

    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;

    if (xe==mx) lxe=xe-1;
    if (ye==my) lye=ye-1;
    if (ze==mz) lze=ze-1;

    // --- Part 1: Calculate the location of i-face centers (Centx) ---
    ierr = DMGetCoordinatesLocal(user->da, &lCoords); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, lCoords, &coor); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Centy, &centy); CHKERRQ(ierr);
    //  ierr = DMDAVecGetArray(user->fda, user->lGridSpace,&gs); CHKERRQ(ierr);

    // Loop over the ghosted region to calculate all local face centers
    for (PetscInt k = gzs + 1; k < gze; k++) {
        for (PetscInt j = gys; j < gye; j++) {
            for (PetscInt i = gxs + 1; i < gxe; i++) {
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

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d:   j-face centers (Centx) calculated and ghosts updated.\n", user->simCtx->rank);

    // --- Part 2: Calculate metrics using face-centered coordinates ---
    ierr = DMDAVecGetArrayRead(user->fda, user->Centy, &centy_const); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->JCsi, &jcsi); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->JEta, &jeta); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->JZet, &jzet); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->JAj, &jaj); CHKERRQ(ierr);

    // Loop over the OWNED region where we will store the final metrics
    for (PetscInt k=lzs; k<lze; k++) {
        for (PetscInt j=ys; j<lye; j++) {
            for (PetscInt i=lxs; i<lxe; i++) {

                // --- Stencil Logic for d/dcsi (derivative in i-direction) ---
                if (i == 1) { // Forward difference at the domain's min-i boundary
                    dxdc = centy_const[k][j][i+1].x - centy_const[k][j][i].x;
                    dydc = centy_const[k][j][i+1].y - centy_const[k][j][i].y;
                    dzdc = centy_const[k][j][i+1].z - centy_const[k][j][i].z;
                } else if (i == mx - 2) { // Backward difference at the domain's max-i boundary
                    dxdc = centy_const[k][j][i].x - centy_const[k][j][i-1].x;
                    dydc = centy_const[k][j][i].y - centy_const[k][j][i-1].y;
                    dzdc = centy_const[k][j][i].z - centy_const[k][j][i-1].z;
                } else { // Central difference in the interior
                    dxdc = 0.5 * (centy_const[k][j][i+1].x - centy_const[k][j][i-1].x);
                    dydc = 0.5 * (centy_const[k][j][i+1].y - centy_const[k][j][i-1].y);
                    dzdc = 0.5 * (centy_const[k][j][i+1].z - centy_const[k][j][i-1].z);
                }

                // --- Stencil Logic for d/deta (derivative in j-direction) ---
                if (j == 0) { // Forward difference
                    dxde = centy_const[k][j+1][i].x - centy_const[k][j][i].x;
                    dyde = centy_const[k][j+1][i].y - centy_const[k][j][i].y;
                    dzde = centy_const[k][j+1][i].z - centy_const[k][j][i].z;
                } else if (j == my - 2) { // Backward difference
                    dxde = centy_const[k][j][i].x - centy_const[k][j-1][i].x;
                    dyde = centy_const[k][j][i].y - centy_const[k][j-1][i].y;
                    dzde = centy_const[k][j][i].z - centy_const[k][j-1][i].z;
                } else { // Central difference
                    dxde = 0.5 * (centy_const[k][j+1][i].x - centy_const[k][j-1][i].x);
                    dyde = 0.5 * (centy_const[k][j+1][i].y - centy_const[k][j-1][i].y);
                    dzde = 0.5 * (centy_const[k][j+1][i].z - centy_const[k][j-1][i].z);
                }

                // --- Stencil Logic for d/dzeta (derivative in k-direction) ---
                if (k == 1) { // Forward difference
                    dxdz = centy_const[k+1][j][i].x - centy_const[k][j][i].x;
                    dydz = centy_const[k+1][j][i].y - centy_const[k][j][i].y;
                    dzdz = centy_const[k+1][j][i].z - centy_const[k][j][i].z;
                } else if (k == mz - 2) { // Backward difference
                    dxdz = centy_const[k][j][i].x - centy_const[k-1][j][i].x;
                    dydz = centy_const[k][j][i].y - centy_const[k-1][j][i].y;
                    dzdz = centy_const[k][j][i].z - centy_const[k-1][j][i].z;
                } else { // Central difference
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

    ierr = DMDAVecRestoreArrayRead(user->fda, user->Centy, &centy_const); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->JCsi, &jcsi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->JEta, &jeta); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->JZet, &jzet); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, user->JAj, &jaj); CHKERRQ(ierr);
    
    // --- Part 3: Assemble global vectors and update local ghosts ---
    ierr = VecAssemblyBegin(user->JCsi); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->JCsi); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->JEta); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->JEta); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->JZet); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->JZet); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->JAj); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->JAj); CHKERRQ(ierr);

    ierr = UpdateLocalGhosts(user, "JCsi"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "JEta"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "JZet"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "JAj"); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeJFaceMetrics"
/**
 * @brief Computes metrics centered on constant-k faces (k-faces).
 *
 * This function calculates the metric terms (`KCsi`, `KEta`, `KZet`) and the
 * inverse Jacobian (`KAj`) located at the geometric center of each constant-k
 * face. This is a critical step for staggered-grid finite difference schemes.
 *
 * The process is a direct and faithful refactoring of the corresponding logic
 * from the legacy `FormMetrics` function:
 * 1.  It first calculates the physical (x,y,z) coordinates of the center of
 *     each i-face and stores them in the `user->Centz` vector.
 * 2.  It then uses a boundary-aware, second-order finite difference stencil on
 *     the `Centz` field to compute the derivatives (e.g., d(x)/d(csi)).
 *     - Central differences are used in the grid interior.
 *     - One-sided differences are used at the physical domain boundaries.
 * 3.  Finally, these derivatives are used to compute the final metric terms and
 *     the inverse Jacobian, which are stored in their respective `Vec` objects.
 *
 * @param user The UserCtx for a specific grid level. This function populates
 *             the `user->KCsi`, `user->KEta`, `user->KZet`, and `user->KAj` vectors.
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
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
    PetscInt gxs = info.gxs, gxe = info.gxs + info.gxm;
    PetscInt gys = info.gys, gye = info.gys + info.gym;
    PetscInt gzs = info.gzs, gze = info.gzs + info.gzm;
    
    PetscInt lxs = xs; PetscInt lxe = xe;
    PetscInt lys = ys; PetscInt lye = ye;
    PetscInt lzs = zs; PetscInt lze = ze;

    if (xs==0) lxs = xs+1;
    if (ys==0) lys = ys+1;
    if (zs==0) lzs = zs+1;

    if (xe==mx) lxe=xe-1;
    if (ye==my) lye=ye-1;
    if (ze==mz) lze=ze-1;

    // --- Part 1: Calculate the location of i-face centers (Centx) ---
    ierr = DMGetCoordinatesLocal(user->da, &lCoords); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, lCoords, &coor); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->Centz, &centz); CHKERRQ(ierr);
    //  ierr = DMDAVecGetArray(user->fda, user->lGridSpace,&gs); CHKERRQ(ierr);

    // Loop over the ghosted region to calculate all local face centers
    for (PetscInt k = gzs; k < gze; k++) {
        for (PetscInt j = gys; j < gye; j++) {
            for (PetscInt i = gxs + 1; i < gxe; i++) {
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

    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d:   k-face centers (Centx) calculated and ghosts updated.\n", user->simCtx->rank);

    // --- Part 2: Calculate metrics using face-centered coordinates ---
    ierr = DMDAVecGetArrayRead(user->fda, user->Centz, &centz_const); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->KCsi, &kcsi); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->KEta, &keta); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->fda, user->KZet, &kzet); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, user->KAj, &kaj); CHKERRQ(ierr);

    // Loop over the OWNED region where we will store the final metrics
    for (PetscInt k=zs; k<lze; k++) {
        for (PetscInt j=lys; j<lye; j++) {
            for (PetscInt i=lxs; i<lxe; i++) {

                // --- Stencil Logic for d/dcsi (derivative in i-direction) ---
                if (i == 1) { // Forward difference at the domain's min-i boundary
                    dxdc = centz_const[k][j][i+1].x - centz_const[k][j][i].x;
                    dydc = centz_const[k][j][i+1].y - centz_const[k][j][i].y;
                    dzdc = centz_const[k][j][i+1].z - centz_const[k][j][i].z;
                } else if (i == mx - 2) { // Backward difference at the domain's max-i boundary
                    dxdc = centz_const[k][j][i].x - centz_const[k][j][i-1].x;
                    dydc = centz_const[k][j][i].y - centz_const[k][j][i-1].y;
                    dzdc = centz_const[k][j][i].z - centz_const[k][j][i-1].z;
                } else { // Central difference in the interior
                    dxdc = 0.5 * (centz_const[k][j][i+1].x - centz_const[k][j][i-1].x);
                    dydc = 0.5 * (centz_const[k][j][i+1].y - centz_const[k][j][i-1].y);
                    dzdc = 0.5 * (centz_const[k][j][i+1].z - centz_const[k][j][i-1].z);
                }

                // --- Stencil Logic for d/deta (derivative in j-direction) ---
                if (j == 1) { // Forward difference
                    dxde = centz_const[k][j+1][i].x - centz_const[k][j][i].x;
                    dyde = centz_const[k][j+1][i].y - centz_const[k][j][i].y;
                    dzde = centz_const[k][j+1][i].z - centz_const[k][j][i].z;
                } else if (j == my - 2) { // Backward difference
                    dxde = centz_const[k][j][i].x - centz_const[k][j-1][i].x;
                    dyde = centz_const[k][j][i].y - centz_const[k][j-1][i].y;
                    dzde = centz_const[k][j][i].z - centz_const[k][j-1][i].z;
                } else { // Central difference
                    dxde = 0.5 * (centz_const[k][j+1][i].x - centz_const[k][j-1][i].x);
                    dyde = 0.5 * (centz_const[k][j+1][i].y - centz_const[k][j-1][i].y);
                    dzde = 0.5 * (centz_const[k][j+1][i].z - centz_const[k][j-1][i].z);
                }

                // --- Stencil Logic for d/dzeta (derivative in k-direction) ---
                if (k == 0) { // Forward difference
                    dxdz = centz_const[k+1][j][i].x - centz_const[k][j][i].x;
                    dydz = centz_const[k+1][j][i].y - centz_const[k][j][i].y;
                    dzdz = centz_const[k+1][j][i].z - centz_const[k][j][i].z;
                } else if (k == mz - 2) { // Backward difference
                    dxdz = centz_const[k][j][i].x - centz_const[k-1][j][i].x;
                    dydz = centz_const[k][j][i].y - centz_const[k-1][j][i].y;
                    dzdz = centz_const[k][j][i].z - centz_const[k-1][j][i].z;
                } else { // Central difference
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

    ierr = DMDAVecRestoreArrayRead(user->fda, user->Centz, &centz_const); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->KCsi, &kcsi); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->KEta, &keta); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->fda, user->KZet, &kzet); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(user->da, user->KAj, &kaj); CHKERRQ(ierr);
    
    // --- Part 3: Assemble global vectors and update local ghosts ---
    ierr = VecAssemblyBegin(user->KCsi); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->KCsi); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->KEta); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->KEta); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->KZet); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->KZet); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(user->KAj); CHKERRQ(ierr); ierr = VecAssemblyEnd(user->KAj); CHKERRQ(ierr);

    ierr = UpdateLocalGhosts(user, "KCsi"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "KEta"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "KZet"); CHKERRQ(ierr);
    ierr = UpdateLocalGhosts(user, "KAj"); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MetricsDivergence"
/**
 * @brief Performs a diagnostic check on the divergence of the face area metric vectors.
 *
 * For a closed cell, the sum of the face area vectors should be zero (Gauss's
 * divergence theorem). This function computes a measure of this divergence and
 * reports the maximum value over the domain. A small value indicates a
 * well-formed grid. This is a direct adaptation of the legacy function.
 *
 * @param user The UserCtx for a specific grid level (typically the finest).
 * @return PetscErrorCode 0 on success, or a PETSc error code on failure.
 */
PetscErrorCode MetricsDivergence(UserCtx *user)
{
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    Vec            Div;
    PetscReal    ***div_arr;
    const Cmpnts ***csi_arr, ***eta_arr, ***zet_arr;
    const PetscScalar ***aj_arr;
    PetscReal      maxdiv;
    PetscInt       max_idx;

    PetscFunctionBeginUser;

    PROFILE_FUNCTION_BEGIN;

    LOG_ALLOW(LOCAL, LOG_INFO, "Rank %d: Performing metric divergence check for level %d block %d...\n", user->simCtx->rank, user->thislevel, user->_this);

    ierr = DMDAGetLocalInfo(user->da, &info); CHKERRQ(ierr);

    ierr = DMDAVecGetArrayRead(user->fda, user->lCsi, &csi_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lEta, &eta_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->fda, user->lZet, &zet_arr); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(user->da, user->lAj, &aj_arr); CHKERRQ(ierr);
    
    ierr = VecDuplicate(user->P, &Div); CHKERRQ(ierr);
    ierr = VecSet(Div, 0.0); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(user->da, Div, &div_arr); CHKERRQ(ierr);
    
    // Loop over the interior OWNED cells
    for (PetscInt k=info.xs+1; k<info.xs+info.xm; k++) {
        for (PetscInt j=info.ys+1; j<info.ys+info.ym; j++) {
            for (PetscInt i=info.xs+1; i<info.xs+info.xm; i++) {
                PetscReal div_x = (csi_arr[k][j][i].x - csi_arr[k][j][i-1].x) + (eta_arr[k][j][i].x - eta_arr[k][j-1][i].x) + (zet_arr[k][j][i].x - zet_arr[k-1][j][i].x);
                PetscReal div_y = (csi_arr[k][j][i].y - csi_arr[k][j][i-1].y) + (eta_arr[k][j][i].y - eta_arr[k][j-1][i].y) + (zet_arr[k][j][i].y - zet_arr[k-1][j][i].y);
                PetscReal div_z = (csi_arr[k][j][i].z - csi_arr[k][j][i-1].z) + (eta_arr[k][j][i].z - eta_arr[k][j-1][i].z) + (zet_arr[k][j][i].z - zet_arr[k-1][j][i].z);
                div_arr[k][j][i] = PetscAbsScalar((div_x + div_y + div_z) * aj_arr[k][j][i]);
            }
        }
    }

    ierr = DMDAVecRestoreArray(user->da, Div, &div_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lCsi, &csi_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lEta, &eta_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->fda, user->lZet, &zet_arr); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArrayRead(user->da, user->lAj, &aj_arr); CHKERRQ(ierr);
    
    ierr = VecMax(Div, &max_idx, &maxdiv); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL, LOG_INFO, "Maximum metric divergence: %e\n", maxdiv);

    ierr = VecDestroy(&Div); CHKERRQ(ierr);

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}

static PetscInt Gidx(PetscInt i, PetscInt j, PetscInt k, UserCtx *user)

{
  PetscInt nidx;
  DMDALocalInfo	info = user->info;

  PetscInt	mx = info.mx, my = info.my;
  
  AO ao;
  DMDAGetAO(user->da, &ao);
  nidx=i+j*mx+k*mx*my;
  
  AOApplicationToPetsc(ao,1,&nidx);
  
  return (nidx);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeMetricsDivergence"
/**
 * @brief Computes the divergence of the grid metrics and identifies the maximum value.
 *
 * This function serves as a diagnostic tool to assess the quality of the grid
 * metrics. It calculates the divergence of the face metrics (Csi, Eta, Zet)
 * and scales it by the inverse of the cell Jacobian. The maximum divergence
 * value is then located, and its grid coordinates are printed to the console,
 * helping to pinpoint areas of potential grid quality issues.
 *
 * @param user The UserCtx, containing all necessary grid data.
 * @return PetscErrorCode
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

  PetscReal MaxFlatIndex;
  VecMax(Div, &MaxFlatIndex, &maxdiv);
  LOG_ALLOW(GLOBAL,LOG_INFO,"The Maximum Metric Divergence is %e at flat index %d.\n",maxdiv,MaxFlatIndex);

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (Gidx(i,j,k,user) == MaxFlatIndex) {
	  LOG_ALLOW(GLOBAL,LOG_INFO,"The Maximum Metric Divergence(%e) is at location [%d][%d][%d]. \n", maxdiv,k,j,i);
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

#undef __FUNCT__
#define __FUNCT__ "ComputeMetricNorms"
/**
 * @brief Computes the max-min values of the grid metrics.
 *
 * This function serves as a diagnostic tool to assess the quality of the grid
 * metrics. It calculates the bounds of the face metrics (Csi, Eta, Zet).
 *
 * @param user The UserCtx, containing all necessary grid data.
 * @return PetscErrorCode
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

#undef __FUNCT__
#define __FUNCT__ "ComputeAllGridMetrics"
/**
 * @brief Orchestrates the calculation of all grid metrics.
 *
 * This function iterates through every UserCtx in the multigrid and multi-block
 * hierarchy. For each context, it calls a series of modern, modular helper
 * functions to compute the face metrics (Csi, Eta, Zet), the cell-centered
 * inverse Jacobian (Aj), and to validate the grid's orientation.
 *
 * @param simCtx The master SimCtx, containing the configured UserCtx hierarchy.
 * @return PetscErrorCode
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

            // Diagnostics
	    ierr = ComputeMetricNorms(user);
            if (level == usermg->mglevels - 1) {
                ierr = ComputeMetricsDivergence(user); CHKERRQ(ierr);
            }
        }
    }

    LOG_ALLOW(GLOBAL, LOG_INFO, "Grid metrics calculation complete.\n");

    PROFILE_FUNCTION_END;

    PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------------- */
/* End of Metric.c */

