#include "rhs.h"
#include "poisson.h"   /* the face gradient shared with the projection */
#include "verification_sources.h"
#include "les.h"   /* SymTensor kernels for the gradient (Clark) closure term */

/** @brief Owned-cell bounds of the residual: interior cells, excluding the dummy layers. */
static inline void RHSInteriorBounds(const DMDALocalInfo *info, PetscInt lo[3], PetscInt hi[3])
{
    const PetscInt m[3] = {info->mx, info->my, info->mz};
    const PetscInt start[3] = {info->xs, info->ys, info->zs};
    const PetscInt end[3] = {info->xs + info->xm, info->ys + info->ym, info->zs + info->zm};

    for (PetscInt a = 0; a < 3; a++) {
        lo[a] = (start[a] == 0) ? 1 : start[a];
        hi[a] = (end[a] == m[a]) ? m[a] - 1 : end[a];
    }
}

/** Cells whose `nvert` exceeds this are blanked inner-body cells, which convect as fluid. */
#define RHS_INNER_BLANK 7.0

/** @brief Whether a cell takes part in convection: fluid, or a blanked inner-body cell. */
static inline PetscBool RHSConvectsAsFluid(PetscReal nvert)
{
    return (PetscBool)(nvert < 0.1 || nvert > RHS_INNER_BLANK);
}

/** @brief Cartesian velocity at cell @p c displaced by @p d along axis @p n. */
static inline Cmpnts RHSVelocityAt(const Cmpnts ***ucat, const PetscInt c[3], PetscInt n, PetscInt d)
{
    PetscInt cc[3] = {c[0], c[1], c[2]};

    cc[n] += d;
    return ucat[cc[2]][cc[1]][cc[0]];
}

/**
 * @brief QUICK flux through the face between cell @p c and `c + e_n`, upwinded by the sign of
 *        @p ucon. @p m1 and @p p2 are the offsets used for the far upwind cells, -1 and +2 in
 *        the interior; a one-sided face substitutes 0 or +1 for the side it cannot reach.
 */
static inline Cmpnts RHSQuickFlux(const Cmpnts ***ucat, const PetscInt c[3], PetscInt n,
                                  PetscReal ucon, PetscInt m1, PetscInt p2)
{
    const PetscReal coef = 0.125;
    const PetscReal up = ucon + fabs(ucon), um = ucon - fabs(ucon);
    const Cmpnts    u0 = RHSVelocityAt(ucat, c, n, 0), u1 = RHSVelocityAt(ucat, c, n, 1);
    const Cmpnts    um1 = RHSVelocityAt(ucat, c, n, m1), up2 = RHSVelocityAt(ucat, c, n, p2);
    Cmpnts          f;

    f.x = um * (coef * (-up2.x -2. * u1.x +3. * u0.x) +u1.x) +
          up * (coef * (-um1.x -2. * u0.x +3. * u1.x) +u0.x);
    f.y = um * (coef * (-up2.y -2. * u1.y +3. * u0.y) +u1.y) +
          up * (coef * (-um1.y -2. * u0.y +3. * u1.y) +u0.y);
    f.z = um * (coef * (-up2.z -2. * u1.z +3. * u0.z) +u1.z) +
          up * (coef * (-um1.z -2. * u0.z +3. * u1.z) +u0.z);
    return f;
}

/** @brief Central flux through the face between cell @p c and `c + e_n`. */
static inline Cmpnts RHSCentralFlux(const Cmpnts ***ucat, const PetscInt c[3], PetscInt n, PetscReal ucon)
{
    const Cmpnts u0 = RHSVelocityAt(ucat, c, n, 0), u1 = RHSVelocityAt(ucat, c, n, 1);
    Cmpnts       f;

    f.x = ucon * ( u0.x + u1.x );
    f.y = ucon * ( u0.y + u1.y );
    f.z = ucon * ( u0.z + u1.z );
    return f;
}

/**
 * @brief Convective flux of Cartesian momentum through the face between cell @p c and
 *        `c + e_n`, carried by half the contravariant flux on that face.
 *
 * Interior faces use the central flux when @p central (central differencing or LES) and QUICK
 * otherwise; with central differencing the two end faces of an axis are central as well. A face
 * at the low end of an axis or behind a solid cell is QUICK one-sided toward +n, and one at the
 * high end or ahead of a solid cell one-sided toward -n, except across a periodic seam, which
 * keeps the full stencil. A face meeting none of these carries no flux.
 */
static inline Cmpnts RHSConvectiveFaceFlux(const Cmpnts ***ucat, const PetscReal ***nvert,
                                           const PetscInt c[3], PetscInt n, PetscReal ucon,
                                           PetscInt m, PetscBool periodic, PetscBool central)
{
    const PetscInt  s = c[n];
    const PetscReal nv_0 = nvert[c[2]][c[1]][c[0]];
    const PetscReal nv_p = nvert[c[2] + (n == 2)][c[1] + (n == 1)][c[0] + (n == 0)];
    const Cmpnts    zero = {0.0, 0.0, 0.0};
    /* The cell below exists in the local array everywhere but at s = 0 on a bounded axis. */
    const PetscBool has_m = (PetscBool)(s > 0 || periodic);
    const PetscReal nv_m = has_m ? nvert[c[2] - (n == 2)][c[1] - (n == 1)][c[0] - (n == 0)] : 0.0;
    const PetscBool full_seam = (PetscBool)(periodic && nv_m < 0.1 && nv_p < 0.1);

    if (s > 0 && s < m - 2 && RHSConvectsAsFluid(nv_p) && RHSConvectsAsFluid(nv_m)) {
        return central ? RHSCentralFlux(ucat, c, n, ucon) : RHSQuickFlux(ucat, c, n, ucon, -1, 2);
    }
    if (central && (s == 0 || s == m - 2) && RHSConvectsAsFluid(nv_p) && RHSConvectsAsFluid(nv_0)) {
        return RHSCentralFlux(ucat, c, n, ucon);
    }
    if (s < m - 2 && (s == 0 || nv_m > 0.1)) {
        return RHSQuickFlux(ucat, c, n, ucon, (full_seam && s == 0) ? -1 : 0, 2);
    }
    if (s > 0 && (s == m - 2 || nv_p > 0.1)) {
        return RHSQuickFlux(ucat, c, n, ucon, -1, (full_seam && s == m - 2) ? 2 : 1);
    }
    return zero;
}

/**
 * @brief Implementation of \ref ComputeConvection().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/rhs.h`.
 * @see ComputeConvection()
 */
PetscErrorCode ComputeConvection(UserCtx *user, Vec Ucont, Vec Ucat, Vec Conv)
{
    const SimCtx        *simCtx = user->simCtx;
    const DMDALocalInfo *info = &user->info;
    const PetscInt       m[3] = {info->mx, info->my, info->mz};
    const PetscBool      central = (PetscBool)(simCtx->les || simCtx->central);
    PetscBool            periodic[3];
    PetscInt             lo[3], hi[3];
    Vec                  Fp[3];
    Cmpnts            ***fp[3], ***conv;
    const Cmpnts      ***ucont, ***ucat;
    const PetscReal   ***nvert;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    RHSInteriorBounds(info, lo, hi);
    PoissonOperator_PeriodicAxes(user, periodic);
    PetscCall(PreparePeriodicQuickStencilFields(user, Ucat, user->lNvert));
    PetscCall(VecSet(Conv, 0.0));
    /* Face fluxes per axis, from fda's pool and zeroed as fresh vectors would be: faces that
       carry no flux must read as zero. */
    for (PetscInt n = 0; n < 3; n++) {
        PetscCall(DMGetLocalVector(user->fda, &Fp[n]));
        PetscCall(VecZeroEntries(Fp[n]));
        PetscCall(DMDAVecGetArray(user->fda, Fp[n], &fp[n]));
    }
    PetscCall(DMDAVecGetArrayRead(user->fda, Ucont, (void *)&ucont));
    PetscCall(DMDAVecGetArrayRead(user->fda, Ucat, (void *)&ucat));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecGetArray(user->fda, Conv, &conv));

    /* Fluxes on every face of the owned cells: along n, from the face below the first cell. */
    for (PetscInt n = 0; n < 3; n++) {
        PetscInt flo[3] = {lo[0], lo[1], lo[2]};

        flo[n] -= 1;
        for (PetscInt k = flo[2]; k < hi[2]; k++) {
            for (PetscInt j = flo[1]; j < hi[1]; j++) {
                for (PetscInt i = flo[0]; i < hi[0]; i++) {
                    const PetscInt  c[3] = {i, j, k};
                    const PetscReal ucon = (&ucont[k][j][i].x)[n] * 0.5;

                    fp[n][k][j][i] = RHSConvectiveFaceFlux(ucat, (const PetscReal ***)nvert, c, n, ucon,
                                                           m[n], periodic[n], central);
                }
            }
        }
    }

    /* Convective term in Cartesian components: the net flux out of each cell. */
    for (PetscInt k = lo[2]; k < hi[2]; k++) {
        for (PetscInt j = lo[1]; j < hi[1]; j++) {
            for (PetscInt i = lo[0]; i < hi[0]; i++) {
                conv[k][j][i].x =
                    fp[0][k][j][i].x - fp[0][k][j][i-1].x +
                    fp[1][k][j][i].x - fp[1][k][j-1][i].x +
                    fp[2][k][j][i].x - fp[2][k-1][j][i].x;
                conv[k][j][i].y =
                    fp[0][k][j][i].y - fp[0][k][j][i-1].y +
                    fp[1][k][j][i].y - fp[1][k][j-1][i].y +
                    fp[2][k][j][i].y - fp[2][k-1][j][i].y;
                conv[k][j][i].z =
                    fp[0][k][j][i].z - fp[0][k][j][i-1].z +
                    fp[1][k][j][i].z - fp[1][k][j-1][i].z +
                    fp[2][k][j][i].z - fp[2][k-1][j][i].z;
            }
        }
    }

    PetscCall(DMDAVecRestoreArray(user->fda, Conv, &conv));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, Ucat, (void *)&ucat));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, Ucont, (void *)&ucont));
    for (PetscInt n = 0; n < 3; n++) {
        PetscCall(DMDAVecRestoreArray(user->fda, Fp[n], &fp[n]));
        PetscCall(DMRestoreLocalVector(user->fda, &Fp[n]));
    }
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Convective term calculated .\n");
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}


/**
 * @brief Implementation of \ref Viscous().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/rhs.h`.
 * @see Viscous()
 */

PetscErrorCode Viscous(UserCtx *user, Vec Ucont, Vec Ucat, Vec Visc)
{
  
  Vec		Csi = user->lCsi, Eta = user->lEta, Zet = user->lZet;

  Cmpnts	***ucont, ***ucat;

  Cmpnts	***csi, ***eta, ***zet;
  Cmpnts	***icsi, ***ieta, ***izet;
  Cmpnts	***jcsi, ***jeta, ***jzet;
  Cmpnts	***kcsi, ***keta, ***kzet;

  PetscReal	***nvert;

  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  PetscInt	i, j, k;
  Vec		Fp1, Fp2, Fp3;
  Cmpnts	***fp1, ***fp2, ***fp3;
  Cmpnts	***visc;
  PetscReal	***aj, ***iaj, ***jaj, ***kaj;

  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal     ajc;

  PetscReal	dudc, dude, dudz, dvdc, dvde, dvdz, dwdc, dwde, dwdz;
  PetscReal	csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;
  PetscReal	g11, g21, g31;
  PetscReal	r11, r21, r31, r12, r22, r32, r13, r23, r33;

  PetscScalar	solid,innerblank;

  // --- CONTEXT ACQUISITION BLOCK ---
  // Get the master simulation context from the UserCtx.
  SimCtx *simCtx = user->simCtx;

  // Create local variables to mirror the legacy globals for minimal code changes.
  const LESModelType les = simCtx->les;
  const	PetscReal ren = simCtx->ren;
  const PetscInt clark = simCtx->les_gradient_model;
  solid = 0.5;
  innerblank = 7.;

  PROFILE_FUNCTION_BEGIN;
  
  DMDAVecGetArray(fda, Ucont, &ucont);
  DMDAVecGetArray(fda, Ucat,  &ucat);
  DMDAVecGetArray(fda, Visc,  &visc);

  DMDAVecGetArray(fda, Csi, &csi);
  DMDAVecGetArray(fda, Eta, &eta);
  DMDAVecGetArray(fda, Zet, &zet);

  DMDAVecGetArray(fda, user->lICsi, &icsi);
  DMDAVecGetArray(fda, user->lIEta, &ieta);
  DMDAVecGetArray(fda, user->lIZet, &izet);

  DMDAVecGetArray(fda, user->lJCsi, &jcsi);
  DMDAVecGetArray(fda, user->lJEta, &jeta);
  DMDAVecGetArray(fda, user->lJZet, &jzet);

  DMDAVecGetArray(fda, user->lKCsi, &kcsi);
  DMDAVecGetArray(fda, user->lKEta, &keta);
  DMDAVecGetArray(fda, user->lKZet, &kzet);

  DMDAVecGetArray(da, user->lNvert, &nvert);

  // Pooled work vectors, zeroed as the fresh duplicates they replace were (see ComputeRHS()).
  DMGetLocalVector(fda, &Fp1); VecZeroEntries(Fp1);
  DMGetLocalVector(fda, &Fp2); VecZeroEntries(Fp2);
  DMGetLocalVector(fda, &Fp3); VecZeroEntries(Fp3);

  DMDAVecGetArray(fda, Fp1, &fp1);
  DMDAVecGetArray(fda, Fp2, &fp2);
  DMDAVecGetArray(fda, Fp3, &fp3);

  DMDAVecGetArray(da, user->lAj, &aj);

  info = user->info;
  
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;

  /* First we calculate the flux on cell surfaces. Stored on the upper integer
     node. For example, along i direction, the flux are stored at node 0:mx-2*/
  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;


  if (xe==mx) lxe=xe-1;
  if (ye==my) lye=ye-1;
  if (ze==mz) lze=ze-1;
 
  VecSet(Visc,0.0);
  
  PetscReal ***lnu_t;
 
  PetscReal ***lnu_wall = NULL;

  if(les) {
    DMDAVecGetArray(da, user->lNu_t, &lnu_t);
  }
  /* A wall model replaces the subgrid viscosity at its own face rather than adding to
     it, so the flux there carries the stress the model computed instead of the
     molecular fraction of it. */
  if (user->simCtx->wallfunction) {
    DMDAVecGetArray(da, user->lNu_Wall, &lnu_wall);
  }

  /* The visc flux on each surface center is stored at previous integer node */

  DMDAVecGetArray(da, user->lIAj, &iaj);
 /*  for (k=zs; k<ze; k++) { */
/*     for (j=ys; j<ye; j++) { */
/*       for (i=xs; i<xe; i++) { */
/* 	if (i==1 && (j==0 ||j==1 || j==2) && (k==21 || k==22|| k==20)) */
/* 	  PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d u is %.15le v is %.15le  w is %.15le \n",i,j,k,ucat[k][j][i].x,ucat[k][j][i].y,ucat[k][j][i].z ); */
/*       } */
/*     } */
/*   } */
  // i direction
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs-1; i<lxe; i++) {
	
	dudc = ucat[k][j][i+1].x - ucat[k][j][i].x;
	dvdc = ucat[k][j][i+1].y - ucat[k][j][i].y;
	dwdc = ucat[k][j][i+1].z - ucat[k][j][i].z;
	
	if ((nvert[k][j+1][i  ]> solid && nvert[k][j+1][i  ]<innerblank)  ||
	    (nvert[k][j+1][i+1]> solid && nvert[k][j+1][i+1]<innerblank)) {
	  dude = (ucat[k][j  ][i+1].x + ucat[k][j  ][i].x -
		  ucat[k][j-1][i+1].x - ucat[k][j-1][i].x) * 0.5;
	  dvde = (ucat[k][j  ][i+1].y + ucat[k][j  ][i].y -
		  ucat[k][j-1][i+1].y - ucat[k][j-1][i].y) * 0.5;
	  dwde = (ucat[k][j  ][i+1].z + ucat[k][j  ][i].z -
		  ucat[k][j-1][i+1].z - ucat[k][j-1][i].z) * 0.5;
	}
	else if  ((nvert[k][j-1][i  ]> solid && nvert[k][j-1][i  ]<innerblank)  ||
		  (nvert[k][j-1][i+1]> solid && nvert[k][j-1][i+1]<innerblank)) {
	  dude = (ucat[k][j+1][i+1].x + ucat[k][j+1][i].x -
		  ucat[k][j  ][i+1].x - ucat[k][j  ][i].x) * 0.5;
	  dvde = (ucat[k][j+1][i+1].y + ucat[k][j+1][i].y -
		  ucat[k][j  ][i+1].y - ucat[k][j  ][i].y) * 0.5;
	  dwde = (ucat[k][j+1][i+1].z + ucat[k][j+1][i].z -
		  ucat[k][j  ][i+1].z - ucat[k][j  ][i].z) * 0.5;
	}
	else {
	  dude = (ucat[k][j+1][i+1].x + ucat[k][j+1][i].x -
		  ucat[k][j-1][i+1].x - ucat[k][j-1][i].x) * 0.25;
	  dvde = (ucat[k][j+1][i+1].y + ucat[k][j+1][i].y -
		  ucat[k][j-1][i+1].y - ucat[k][j-1][i].y) * 0.25;
	  dwde = (ucat[k][j+1][i+1].z + ucat[k][j+1][i].z -
		  ucat[k][j-1][i+1].z - ucat[k][j-1][i].z) * 0.25;
	}
	
	if ((nvert[k+1][j][i  ]> solid && nvert[k+1][j][i  ]<innerblank)||
	    (nvert[k+1][j][i+1]> solid && nvert[k+1][j][i+1]<innerblank)) {
	  dudz = (ucat[k  ][j][i+1].x + ucat[k  ][j][i].x -
		  ucat[k-1][j][i+1].x - ucat[k-1][j][i].x) * 0.5;
	  dvdz = (ucat[k  ][j][i+1].y + ucat[k  ][j][i].y -
		  ucat[k-1][j][i+1].y - ucat[k-1][j][i].y) * 0.5;
	  dwdz = (ucat[k  ][j][i+1].z + ucat[k  ][j][i].z -
		  ucat[k-1][j][i+1].z - ucat[k-1][j][i].z) * 0.5;
	}
	else if ((nvert[k-1][j][i  ]> solid && nvert[k-1][j][i  ]<innerblank) ||
		 (nvert[k-1][j][i+1]> solid && nvert[k-1][j][i+1]<innerblank)) {

	  dudz = (ucat[k+1][j][i+1].x + ucat[k+1][j][i].x -
		  ucat[k  ][j][i+1].x - ucat[k  ][j][i].x) * 0.5;
	  dvdz = (ucat[k+1][j][i+1].y + ucat[k+1][j][i].y -
		  ucat[k  ][j][i+1].y - ucat[k  ][j][i].y) * 0.5;
	  dwdz = (ucat[k+1][j][i+1].z + ucat[k+1][j][i].z -
		  ucat[k  ][j][i+1].z - ucat[k  ][j][i].z) * 0.5;
	}
	else {
	  dudz = (ucat[k+1][j][i+1].x + ucat[k+1][j][i].x -
		  ucat[k-1][j][i+1].x - ucat[k-1][j][i].x) * 0.25;
	  dvdz = (ucat[k+1][j][i+1].y + ucat[k+1][j][i].y -
		  ucat[k-1][j][i+1].y - ucat[k-1][j][i].y) * 0.25;
	  dwdz = (ucat[k+1][j][i+1].z + ucat[k+1][j][i].z -
		  ucat[k-1][j][i+1].z - ucat[k-1][j][i].z) * 0.25;
	}

 	csi0 = icsi[k][j][i].x;
	csi1 = icsi[k][j][i].y;
	csi2 = icsi[k][j][i].z;

	eta0 = ieta[k][j][i].x;
	eta1 = ieta[k][j][i].y;
	eta2 = ieta[k][j][i].z;
	
	zet0 = izet[k][j][i].x;
	zet1 = izet[k][j][i].y;
	zet2 = izet[k][j][i].z;

	g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
	g21 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
	g31 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;

	r11 = dudc * csi0 + dude * eta0 + dudz * zet0;
	r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;
	r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;

	r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
	r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
	r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;

	r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
	r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
	r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;

	ajc = iaj[k][j][i]; 

	double nu = 1./ren, nu_t=0;
	
	if( les ) {
	  //nu_t = pow( 0.5 * ( sqrt(lnu_t[k][j][i]) + sqrt(lnu_t[k][j][i+1]) ), 2.0) * Sabs;
	  nu_t = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j][i+1]);
	  /* Zero is right for a wall-resolved run, where the subgrid stress vanishes at
	     the wall. With a wall model the same face is where the modelled stress has to
	     enter, so it takes the model's effective viscosity instead. */
	  if ( BoundaryFaceHasWallShear(&user->boundary_faces[BC_FACE_NEG_X]) && i==0 ) nu_t = lnu_wall ? lnu_wall[k][j][1] : 0.0;
	  if ( BoundaryFaceHasWallShear(&user->boundary_faces[BC_FACE_POS_X]) && i==mx-2 ) nu_t = lnu_wall ? lnu_wall[k][j][mx-2] : 0.0;    
	  fp1[k][j][i].x = (g11 * dudc + g21 * dude + g31 * dudz + r11 * csi0 + r21 * csi1 + r31 * csi2) * ajc * (nu_t); 
	  fp1[k][j][i].y = (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * csi0 + r22 * csi1 + r32 * csi2) * ajc * (nu_t);
	  fp1[k][j][i].z = (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * csi0 + r23 * csi1 + r33 * csi2) * ajc * (nu_t);
	}
	else {
	  fp1[k][j][i].x = 0;
	  fp1[k][j][i].y = 0;
	  fp1[k][j][i].z = 0;
	}

	fp1[k][j][i].x += (g11 * dudc + g21 * dude + g31 * dudz+ r11 * csi0 + r21 * csi1 + r31 * csi2 ) * ajc * (nu);
	fp1[k][j][i].y += (g11 * dvdc + g21 * dvde + g31 * dvdz+ r12 * csi0 + r22 * csi1 + r32 * csi2 ) * ajc * (nu);
	fp1[k][j][i].z += (g11 * dwdc + g21 * dwde + g31 * dwdz+ r13 * csi0 + r23 * csi1 + r33 * csi2 ) * ajc * (nu);


	if(clark) {
	  /* dudc, dude, dudz are differences between neighbouring cells along each grid
	     direction, so each already equals that direction's spacing times the physical
	     derivative. Their outer products are therefore Delta_k^2 du_i/dx_k du_j/dx_k
	     as they stand. Multiplying by a squared cell extent as well - which this term
	     did, carried over from the legacy - made the stress scale as Delta^4 and left
	     it smaller than intended by Delta^2, 1e-4 to 1e-6 on a production mesh; the
	     extents were also Cartesian, so they did not belong to the directions they
	     were paired with once the grid turned. */
	  
	  /* The gradient-model tensor is the sum over the three computational
	     directions of the per-cell velocity difference's outer product with itself;
	     the cell's own extent in each direction is already inside that difference.
	     It is symmetric by construction, so it is carried as one. */
	  const Cmpnts grad_csi = { dudc, dvdc, dwdc };
	  const Cmpnts grad_eta = { dude, dvde, dwde };
	  const Cmpnts grad_zet = { dudz, dvdz, dwdz };
	  const SymTensor gradient_tensor =
	      SymTensorCombine(1.0, SymTensorSelfOuter(grad_csi), 1.0,
	                       SymTensorCombine(1.0, SymTensorSelfOuter(grad_eta),
	                                        1.0, SymTensorSelfOuter(grad_zet)));
	  
	  {
	    const Cmpnts face_normal = { csi0, csi1, csi2 };
	    const Cmpnts gradient_flux = SymTensorTimesVector(gradient_tensor, face_normal);

	    fp1[k][j][i].x -= gradient_flux.x / 12.;
	    fp1[k][j][i].y -= gradient_flux.y / 12.;
	    fp1[k][j][i].z -= gradient_flux.z / 12.;
	  }
	}

      }
    }
  }
  DMDAVecRestoreArray(da, user->lIAj, &iaj);  
 

  // j direction
  DMDAVecGetArray(da, user->lJAj, &jaj);
  for (k=lzs; k<lze; k++) {
    for (j=lys-1; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
 
	if ((nvert[k][j  ][i+1]> solid && nvert[k][j  ][i+1]<innerblank)||
	    (nvert[k][j+1][i+1]> solid && nvert[k][j+1][i+1]<innerblank)) {
	  dudc = (ucat[k][j+1][i  ].x + ucat[k][j][i  ].x -
		  ucat[k][j+1][i-1].x - ucat[k][j][i-1].x) * 0.5;
	  dvdc = (ucat[k][j+1][i  ].y + ucat[k][j][i  ].y -
		  ucat[k][j+1][i-1].y - ucat[k][j][i-1].y) * 0.5;
	  dwdc = (ucat[k][j+1][i  ].z + ucat[k][j][i  ].z -
		  ucat[k][j+1][i-1].z - ucat[k][j][i-1].z) * 0.5;
	}
	else if ((nvert[k][j  ][i-1]> solid && nvert[k][j  ][i-1]<innerblank) ||
		 (nvert[k][j+1][i-1]> solid && nvert[k][j+1][i-1]<innerblank)) {
	  dudc = (ucat[k][j+1][i+1].x + ucat[k][j][i+1].x -
		  ucat[k][j+1][i  ].x - ucat[k][j][i  ].x) * 0.5;
	  dvdc = (ucat[k][j+1][i+1].y + ucat[k][j][i+1].y -
		  ucat[k][j+1][i  ].y - ucat[k][j][i  ].y) * 0.5;
	  dwdc = (ucat[k][j+1][i+1].z + ucat[k][j][i+1].z -
		  ucat[k][j+1][i  ].z - ucat[k][j][i  ].z) * 0.5;
	}
	else {
	  dudc = (ucat[k][j+1][i+1].x + ucat[k][j][i+1].x -
		  ucat[k][j+1][i-1].x - ucat[k][j][i-1].x) * 0.25;
	  dvdc = (ucat[k][j+1][i+1].y + ucat[k][j][i+1].y -
		  ucat[k][j+1][i-1].y - ucat[k][j][i-1].y) * 0.25;
	  dwdc = (ucat[k][j+1][i+1].z + ucat[k][j][i+1].z -
		  ucat[k][j+1][i-1].z - ucat[k][j][i-1].z) * 0.25;
	}

	dude = ucat[k][j+1][i].x - ucat[k][j][i].x;
	dvde = ucat[k][j+1][i].y - ucat[k][j][i].y;
	dwde = ucat[k][j+1][i].z - ucat[k][j][i].z;

	if ((nvert[k+1][j  ][i]> solid && nvert[k+1][j  ][i]<innerblank)||
	    (nvert[k+1][j+1][i]> solid && nvert[k+1][j+1][i]<innerblank)) {
	  dudz = (ucat[k  ][j+1][i].x + ucat[k  ][j][i].x -
		  ucat[k-1][j+1][i].x - ucat[k-1][j][i].x) * 0.5;
	  dvdz = (ucat[k  ][j+1][i].y + ucat[k  ][j][i].y -
		  ucat[k-1][j+1][i].y - ucat[k-1][j][i].y) * 0.5;
	  dwdz = (ucat[k  ][j+1][i].z + ucat[k  ][j][i].z -
		  ucat[k-1][j+1][i].z - ucat[k-1][j][i].z) * 0.5;
	}
	else if ((nvert[k-1][j  ][i]> solid && nvert[k-1][j  ][i]<innerblank)||
		 (nvert[k-1][j+1][i]> solid && nvert[k-1][j+1][i]<innerblank)) {
	  dudz = (ucat[k+1][j+1][i].x + ucat[k+1][j][i].x -
		  ucat[k  ][j+1][i].x - ucat[k  ][j][i].x) * 0.5;
	  dvdz = (ucat[k+1][j+1][i].y + ucat[k+1][j][i].y -
		  ucat[k  ][j+1][i].y - ucat[k  ][j][i].y) * 0.5;
	  dwdz = (ucat[k+1][j+1][i].z + ucat[k+1][j][i].z -
		  ucat[k  ][j+1][i].z - ucat[k  ][j][i].z) * 0.5;
	}
	else {
	  dudz = (ucat[k+1][j+1][i].x + ucat[k+1][j][i].x -
		  ucat[k-1][j+1][i].x - ucat[k-1][j][i].x) * 0.25;
	  dvdz = (ucat[k+1][j+1][i].y + ucat[k+1][j][i].y -
		  ucat[k-1][j+1][i].y - ucat[k-1][j][i].y) * 0.25;
	  dwdz = (ucat[k+1][j+1][i].z + ucat[k+1][j][i].z -
		  ucat[k-1][j+1][i].z - ucat[k-1][j][i].z) * 0.25;
	}

 	csi0 = jcsi[k][j][i].x;
	csi1 = jcsi[k][j][i].y;
	csi2 = jcsi[k][j][i].z;

	eta0 = jeta[k][j][i].x;
	eta1 = jeta[k][j][i].y;
	eta2 = jeta[k][j][i].z;
	
	zet0 = jzet[k][j][i].x;
	zet1 = jzet[k][j][i].y;
	zet2 = jzet[k][j][i].z;


	g11 = csi0 * eta0 + csi1 * eta1 + csi2 * eta2;
	g21 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
	g31 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;

	r11 = dudc * csi0 + dude * eta0 + dudz * zet0;
	r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;
	r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;

	r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
	r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
	r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;

	r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
	r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
	r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;

	//	if (i==1 && j==0 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d dvdc is %.15le dvde is %.15le dvdz is %.15le \n",i,j,k,dvdc,dvde,dvdz);
	//	if (i==1 && j==0 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d dwdc is %.15le dwde is %.15le dwdz is %.15le \n",i,j,k,dwdc,dwde,dwdz);
	//	if (i==1 && j==0 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d jcsi is %.15le jeta is %.15le jzet is %.15le \n",i,j,k,jcsi[k][j][i].z,jeta[k][j][i].z,jzet[k][j][i].z);
	//	if (i==1 && j==0 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i=%d j=%d k=%d r13 is %.15le r23 is %.15le r33 is %.15le \n",i,j,k,r13,r23,r33);



	ajc = jaj[k][j][i];

	double nu = 1./ren, nu_t = 0;
		
	if( les ) {
	  //nu_t = pow( 0.5 * ( sqrt(lnu_t[k][j][i]) + sqrt(lnu_t[k][j+1][i]) ), 2.0) * Sabs;
	  nu_t = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j+1][i]);
	  /* Zero is right for a wall-resolved run, where the subgrid stress vanishes at
	     the wall. With a wall model the same face is where the modelled stress has to
	     enter, so it takes the model's effective viscosity instead. */
	  if ( BoundaryFaceHasWallShear(&user->boundary_faces[BC_FACE_NEG_Y]) && j==0 ) nu_t = lnu_wall ? lnu_wall[k][1][i] : 0.0;
	  if ( BoundaryFaceHasWallShear(&user->boundary_faces[BC_FACE_POS_Y]) && j==my-2 ) nu_t = lnu_wall ? lnu_wall[k][my-2][i] : 0.0;
		
	  fp2[k][j][i].x = (g11 * dudc + g21 * dude + g31 * dudz + r11 * eta0 + r21 * eta1 + r31 * eta2) * ajc * (nu_t);
	  fp2[k][j][i].y = (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * eta0 + r22 * eta1 + r32 * eta2) * ajc * (nu_t);
	  fp2[k][j][i].z = (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * eta0 + r23 * eta1 + r33 * eta2) * ajc * (nu_t);
	}
	else {
	  fp2[k][j][i].x = 0;
	  fp2[k][j][i].y = 0;
	  fp2[k][j][i].z = 0;
	}
		
	fp2[k][j][i].x += (g11 * dudc + g21 * dude + g31 * dudz+ r11 * eta0 + r21 * eta1 + r31 * eta2 ) * ajc * (nu);
	fp2[k][j][i].y += (g11 * dvdc + g21 * dvde + g31 * dvdz+ r12 * eta0 + r22 * eta1 + r32 * eta2 ) * ajc * (nu);
	fp2[k][j][i].z += (g11 * dwdc + g21 * dwde + g31 * dwdz+ r13 * eta0 + r23 * eta1 + r33 * eta2 ) * ajc * (nu);
		
	if(clark) {
	  /* dudc, dude, dudz are differences between neighbouring cells along each grid
	     direction, so each already equals that direction's spacing times the physical
	     derivative. Their outer products are therefore Delta_k^2 du_i/dx_k du_j/dx_k
	     as they stand. Multiplying by a squared cell extent as well - which this term
	     did, carried over from the legacy - made the stress scale as Delta^4 and left
	     it smaller than intended by Delta^2, 1e-4 to 1e-6 on a production mesh; the
	     extents were also Cartesian, so they did not belong to the directions they
	     were paired with once the grid turned. */
			
	  /* The gradient-model tensor is the sum over the three computational
	     directions of the per-cell velocity difference's outer product with itself;
	     the cell's own extent in each direction is already inside that difference.
	     It is symmetric by construction, so it is carried as one. */
	  const Cmpnts grad_csi = { dudc, dvdc, dwdc };
	  const Cmpnts grad_eta = { dude, dvde, dwde };
	  const Cmpnts grad_zet = { dudz, dvdz, dwdz };
	  const SymTensor gradient_tensor =
	      SymTensorCombine(1.0, SymTensorSelfOuter(grad_csi), 1.0,
	                       SymTensorCombine(1.0, SymTensorSelfOuter(grad_eta),
	                                        1.0, SymTensorSelfOuter(grad_zet)));
			
	  {
	    const Cmpnts face_normal = { eta0, eta1, eta2 };
	    const Cmpnts gradient_flux = SymTensorTimesVector(gradient_tensor, face_normal);

	    fp2[k][j][i].x -= gradient_flux.x / 12.;
	    fp2[k][j][i].y -= gradient_flux.y / 12.;
	    fp2[k][j][i].z -= gradient_flux.z / 12.;
	  }
	}
      }
    }
  }

  DMDAVecRestoreArray(da, user->lJAj, &jaj);
  // k direction

  DMDAVecGetArray(da, user->lKAj, &kaj);
  for (k=lzs-1; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if ((nvert[k  ][j][i+1]> solid && nvert[k  ][j][i+1]<innerblank)||
	    (nvert[k+1][j][i+1]> solid && nvert[k+1][j][i+1]<innerblank)) {
	  dudc = (ucat[k+1][j][i  ].x + ucat[k][j][i  ].x -
		  ucat[k+1][j][i-1].x - ucat[k][j][i-1].x) * 0.5;
	  dvdc = (ucat[k+1][j][i  ].y + ucat[k][j][i  ].y -
		  ucat[k+1][j][i-1].y - ucat[k][j][i-1].y) * 0.5;
	  dwdc = (ucat[k+1][j][i  ].z + ucat[k][j][i  ].z -
		  ucat[k+1][j][i-1].z - ucat[k][j][i-1].z) * 0.5;
	}
	else if ((nvert[k  ][j][i-1]> solid && nvert[k  ][j][i-1]<innerblank) ||
		 (nvert[k+1][j][i-1]> solid && nvert[k+1][j][i-1]<innerblank)) {
	  dudc = (ucat[k+1][j][i+1].x + ucat[k][j][i+1].x -
		  ucat[k+1][j][i  ].x - ucat[k][j][i  ].x) * 0.5;
	  dvdc = (ucat[k+1][j][i+1].y + ucat[k][j][i+1].y -
		  ucat[k+1][j][i  ].y - ucat[k][j][i  ].y) * 0.5;
	  dwdc = (ucat[k+1][j][i+1].z + ucat[k][j][i+1].z -
		  ucat[k+1][j][i  ].z - ucat[k][j][i  ].z) * 0.5;
	}
	else {
	  dudc = (ucat[k+1][j][i+1].x + ucat[k][j][i+1].x -
		  ucat[k+1][j][i-1].x - ucat[k][j][i-1].x) * 0.25;
	  dvdc = (ucat[k+1][j][i+1].y + ucat[k][j][i+1].y -
		  ucat[k+1][j][i-1].y - ucat[k][j][i-1].y) * 0.25;
	  dwdc = (ucat[k+1][j][i+1].z + ucat[k][j][i+1].z -
		  ucat[k+1][j][i-1].z - ucat[k][j][i-1].z) * 0.25;
	}
	  
	if ((nvert[k  ][j+1][i]> solid && nvert[k  ][j+1][i]<innerblank)||
	    (nvert[k+1][j+1][i]> solid && nvert[k+1][j+1][i]<innerblank)) {
	  dude = (ucat[k+1][j  ][i].x + ucat[k][j  ][i].x -
		  ucat[k+1][j-1][i].x - ucat[k][j-1][i].x) * 0.5;
	  dvde = (ucat[k+1][j  ][i].y + ucat[k][j  ][i].y -
		  ucat[k+1][j-1][i].y - ucat[k][j-1][i].y) * 0.5;
	  dwde = (ucat[k+1][j  ][i].z + ucat[k][j  ][i].z -
		  ucat[k+1][j-1][i].z - ucat[k][j-1][i].z) * 0.5;
	}
	else if ((nvert[k  ][j-1][i]> solid && nvert[k  ][j-1][i]<innerblank) ||
		 (nvert[k+1][j-1][i]> solid && nvert[k+1][j-1][i]<innerblank)){
	  dude = (ucat[k+1][j+1][i].x + ucat[k][j+1][i].x -
		  ucat[k+1][j  ][i].x - ucat[k][j  ][i].x) * 0.5;
	  dvde = (ucat[k+1][j+1][i].y + ucat[k][j+1][i].y -
		  ucat[k+1][j  ][i].y - ucat[k][j  ][i].y) * 0.5;
	  dwde = (ucat[k+1][j+1][i].z + ucat[k][j+1][i].z -
		  ucat[k+1][j  ][i].z - ucat[k][j  ][i].z) * 0.5;
	}
	else {
	  dude = (ucat[k+1][j+1][i].x + ucat[k][j+1][i].x -
		  ucat[k+1][j-1][i].x - ucat[k][j-1][i].x) * 0.25;
	  dvde = (ucat[k+1][j+1][i].y + ucat[k][j+1][i].y -
		  ucat[k+1][j-1][i].y - ucat[k][j-1][i].y) * 0.25;
	  dwde = (ucat[k+1][j+1][i].z + ucat[k][j+1][i].z -
		  ucat[k+1][j-1][i].z - ucat[k][j-1][i].z) * 0.25;
	}

	dudz = ucat[k+1][j][i].x - ucat[k][j][i].x;
	dvdz = ucat[k+1][j][i].y - ucat[k][j][i].y;
	dwdz = ucat[k+1][j][i].z - ucat[k][j][i].z;


 	csi0 = kcsi[k][j][i].x;
	csi1 = kcsi[k][j][i].y;
	csi2 = kcsi[k][j][i].z;

	eta0 = keta[k][j][i].x;
	eta1 = keta[k][j][i].y;
	eta2 = keta[k][j][i].z;
	
	zet0 = kzet[k][j][i].x;
	zet1 = kzet[k][j][i].y;
	zet2 = kzet[k][j][i].z;


	g11 = csi0 * zet0 + csi1 * zet1 + csi2 * zet2;
	g21 = eta0 * zet0 + eta1 * zet1 + eta2 * zet2;
	g31 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

	r11 = dudc * csi0 + dude * eta0 + dudz * zet0;
	r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;
	r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;

	r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
	r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
	r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;

	r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
	r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
	r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;

	ajc = kaj[k][j][i];

	double nu = 1./ren, nu_t =0;
		
	if( les ) {
	  //nu_t = pow( 0.5 * ( sqrt(lnu_t[k][j][i]) + sqrt(lnu_t[k+1][j][i]) ), 2.0) * Sabs;
	  nu_t = 0.5 * (lnu_t[k][j][i] + lnu_t[k+1][j][i]);
	  /* Zero is right for a wall-resolved run, where the subgrid stress vanishes at
	     the wall. With a wall model the same face is where the modelled stress has to
	     enter, so it takes the model's effective viscosity instead. */
	  if ( BoundaryFaceHasWallShear(&user->boundary_faces[BC_FACE_NEG_Z]) && k==0 ) nu_t = lnu_wall ? lnu_wall[1][j][i] : 0.0;
	  if ( BoundaryFaceHasWallShear(&user->boundary_faces[BC_FACE_POS_Z]) && k==mz-2 ) nu_t = lnu_wall ? lnu_wall[mz-2][j][i] : 0.0;
		
	  fp3[k][j][i].x = (g11 * dudc + g21 * dude + g31 * dudz + r11 * zet0 + r21 * zet1 + r31 * zet2) * ajc * (nu_t);
	  fp3[k][j][i].y = (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * zet0 + r22 * zet1 + r32 * zet2) * ajc * (nu_t);
	  fp3[k][j][i].z = (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * zet0 + r23 * zet1 + r33 * zet2) * ajc * (nu_t);
	}
	else {
	  fp3[k][j][i].x = 0;
	  fp3[k][j][i].y = 0;
	  fp3[k][j][i].z = 0;
	}
	fp3[k][j][i].x += (g11 * dudc + g21 * dude + g31 * dudz + r11 * zet0 + r21 * zet1 + r31 * zet2) * ajc * (nu);//
	fp3[k][j][i].y += (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * zet0 + r22 * zet1 + r32 * zet2) * ajc * (nu);//
	fp3[k][j][i].z += (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * zet0 + r23 * zet1 + r33 * zet2) * ajc * (nu);//

	if(clark) {
	  /* dudc, dude, dudz are differences between neighbouring cells along each grid
	     direction, so each already equals that direction's spacing times the physical
	     derivative. Their outer products are therefore Delta_k^2 du_i/dx_k du_j/dx_k
	     as they stand. Multiplying by a squared cell extent as well - which this term
	     did, carried over from the legacy - made the stress scale as Delta^4 and left
	     it smaller than intended by Delta^2, 1e-4 to 1e-6 on a production mesh; the
	     extents were also Cartesian, so they did not belong to the directions they
	     were paired with once the grid turned. */
			
	  /* The gradient-model tensor is the sum over the three computational
	     directions of the per-cell velocity difference's outer product with itself;
	     the cell's own extent in each direction is already inside that difference.
	     It is symmetric by construction, so it is carried as one. */
	  const Cmpnts grad_csi = { dudc, dvdc, dwdc };
	  const Cmpnts grad_eta = { dude, dvde, dwde };
	  const Cmpnts grad_zet = { dudz, dvdz, dwdz };
	  const SymTensor gradient_tensor =
	      SymTensorCombine(1.0, SymTensorSelfOuter(grad_csi), 1.0,
	                       SymTensorCombine(1.0, SymTensorSelfOuter(grad_eta),
	                                        1.0, SymTensorSelfOuter(grad_zet)));
			
	  {
	    const Cmpnts face_normal = { zet0, zet1, zet2 };
	    const Cmpnts gradient_flux = SymTensorTimesVector(gradient_tensor, face_normal);

	    fp3[k][j][i].x -= gradient_flux.x / 12.;
	    fp3[k][j][i].y -= gradient_flux.y / 12.;
	    fp3[k][j][i].z -= gradient_flux.z / 12.;
	  }
	}
      }
    }
  }

  DMDAVecRestoreArray(da, user->lKAj, &kaj);
  
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	visc[k][j][i].x =
	  (fp1[k][j][i].x - fp1[k][j][i-1].x +
	   fp2[k][j][i].x - fp2[k][j-1][i].x +
	   fp3[k][j][i].x - fp3[k-1][j][i].x);
	
	visc[k][j][i].y =
	  (fp1[k][j][i].y - fp1[k][j][i-1].y +
	   fp2[k][j][i].y - fp2[k][j-1][i].y +
	   fp3[k][j][i].y - fp3[k-1][j][i].y);

	visc[k][j][i].z =
	  (fp1[k][j][i].z - fp1[k][j][i-1].z +
	   fp2[k][j][i].z - fp2[k][j-1][i].z +
	   fp3[k][j][i].z - fp3[k-1][j][i].z);
			
      }
    }
  }
/*   for (k=zs; k<ze; k++) { */
/*     for (j=ys; j<ye; j++) { */
/*       for (i=xs; i<xe; i++) { */
/* 	if (i==1 && j==1 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d  fp1.z is %.15le \n",i,j,k,fp1[k][j][i].z); */
/* 	if (i==0 && j==1 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d  fp1.z is %.15le \n",i,j,k,fp1[k][j][i].z); */
/* 	if (i==1 && j==1 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d  fp2.z is %.15le \n",i,j,k,fp2[k][j][i].z); */
/* 	if (i==1 && j==0 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d  fp2.z is %.15le \n",i,j,k,fp2[k][j][i].z); */
/* 	if (i==1 && j==1 && k==21) PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d  fp3.z is %.15le \n",i,j,k,fp3[k][j][i].z); */
/* 	if (i==1 && j==1 && k==20) PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d  fp3.z is %.15le \n",i,j,k,fp3[k][j][i].z); */
	 
/*       } */
/*     } */
/*   }  */
  DMDAVecRestoreArray(fda, Ucont, &ucont);
  DMDAVecRestoreArray(fda, Ucat,  &ucat);
  DMDAVecRestoreArray(fda, Visc,  &visc);

  DMDAVecRestoreArray(fda, Csi, &csi);
  DMDAVecRestoreArray(fda, Eta, &eta);
  DMDAVecRestoreArray(fda, Zet, &zet);

  DMDAVecRestoreArray(fda, Fp1, &fp1);
  DMDAVecRestoreArray(fda, Fp2, &fp2);
  DMDAVecRestoreArray(fda, Fp3, &fp3);

  DMDAVecRestoreArray(da, user->lAj, &aj);

  DMDAVecRestoreArray(fda, user->lICsi, &icsi);
  DMDAVecRestoreArray(fda, user->lIEta, &ieta);
  DMDAVecRestoreArray(fda, user->lIZet, &izet);

  DMDAVecRestoreArray(fda, user->lJCsi, &jcsi);
  DMDAVecRestoreArray(fda, user->lJEta, &jeta);
  DMDAVecRestoreArray(fda, user->lJZet, &jzet);

  DMDAVecRestoreArray(fda, user->lKCsi, &kcsi);
  DMDAVecRestoreArray(fda, user->lKEta, &keta);
  DMDAVecRestoreArray(fda, user->lKZet, &kzet);

  DMDAVecRestoreArray(da, user->lNvert, &nvert);

  if(les) {
    DMDAVecRestoreArray(da, user->lNu_t, &lnu_t);
  }
  /* Opened independently of the turbulence model, so closed independently of it. */
  if (lnu_wall) DMDAVecRestoreArray(da, user->lNu_Wall, &lnu_wall);
 

  DMRestoreLocalVector(fda, &Fp1);
  DMRestoreLocalVector(fda, &Fp2);
  DMRestoreLocalVector(fda, &Fp3);


  LOG_ALLOW(GLOBAL,LOG_DEBUG,"Viscous terms calculated .\n");

  PROFILE_FUNCTION_END;
  
  return(0);
}

/**
 * @brief Internal helper implementation: `ComputeBodyForces()`.
 * @details Local to this translation unit.
 */
PetscErrorCode ComputeBodyForces(UserCtx *user, Vec Rct)
{
	PetscErrorCode ierr;
    PetscFunctionBeginUser;
	PROFILE_HELPER_BEGIN;

    // --- 1. Apply momentum source for driven channel/pipe flows ---
    // This function will internally check if a driven flow BC is active.
    ierr = ComputeDrivenChannelFlowSource(user, Rct); CHKERRQ(ierr);

    // --- 2. (Future Extension) Apply gravitational force ---
    // if (user->simCtx->gravityEnabled) {
    //     ierr = ApplyGravitationalForce(user, Rhs); CHKERRQ(ierr);
    // }
    //
    // Adding a body force here: see the contract in include/BodyForces.h.
    // In particular, this function runs once per RESIDUAL EVALUATION, not once
    // per timestep, so any force carrying state across calls (a filter, ramp,
    // moving average, or integral term) must gate its update on simCtx->step.

    PROFILE_HELPER_END;
    PetscFunctionReturn(0); 
}

/**
 * @brief Contravariant cell-centre components of the spatial terms, `Rct = aj F_n . (Visc - Conv)`,
 *        with each base vector `F_n` averaged from its two faces to the cell centre.
 */
static PetscErrorCode RHSContravariantSpatialTerms(UserCtx *user, Vec Conv, Vec Visc, Vec Rct)
{
    PetscInt         lo[3], hi[3];
    const Cmpnts  ***csi, ***eta, ***zet, ***conv, ***visc;
    const PetscReal ***aj;
    Cmpnts        ***rct;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    RHSInteriorBounds(&user->info, lo, hi);
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lCsi, (void *)&csi));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lEta, (void *)&eta));
    PetscCall(DMDAVecGetArrayRead(user->fda, user->lZet, (void *)&zet));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lAj, (void *)&aj));
    PetscCall(DMDAVecGetArrayRead(user->fda, Conv, (void *)&conv));
    PetscCall(DMDAVecGetArrayRead(user->fda, Visc, (void *)&visc));
    PetscCall(DMDAVecGetArray(user->fda, Rct, &rct));
    for (PetscInt k = lo[2]; k < hi[2]; k++) {
        for (PetscInt j = lo[1]; j < hi[1]; j++) {
            for (PetscInt i = lo[0]; i < hi[0]; i++) {
                const Cmpnts rc = {visc[k][j][i].x - conv[k][j][i].x,
                                   visc[k][j][i].y - conv[k][j][i].y,
                                   visc[k][j][i].z - conv[k][j][i].z};
                rct[k][j][i].x = aj[k][j][i] *
                    (0.5 * (csi[k][j][i].x + csi[k][j][i-1].x) * rc.x +
                     0.5 * (csi[k][j][i].y + csi[k][j][i-1].y) * rc.y +
                     0.5 * (csi[k][j][i].z + csi[k][j][i-1].z) * rc.z);
                rct[k][j][i].y = aj[k][j][i] *
                    (0.5 * (eta[k][j][i].x + eta[k][j-1][i].x) * rc.x +
                     0.5 * (eta[k][j][i].y + eta[k][j-1][i].y) * rc.y +
                     0.5 * (eta[k][j][i].z + eta[k][j-1][i].z) * rc.z);
                rct[k][j][i].z = aj[k][j][i] *
                    (0.5 * (zet[k][j][i].x + zet[k-1][j][i].x) * rc.x +
                     0.5 * (zet[k][j][i].y + zet[k-1][j][i].y) * rc.y +
                     0.5 * (zet[k][j][i].z + zet[k-1][j][i].z) * rc.z);
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, Rct, &rct));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, Visc, (void *)&visc));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, Conv, (void *)&conv));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lAj, (void *)&aj));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lZet, (void *)&zet));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lEta, (void *)&eta));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, user->lCsi, (void *)&csi));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Face residual from the ghosted cell-centre terms: the average of @p Rct over the two
 *        cells beside each face, less the pressure-gradient flux through it.
 *
 * The flux is the projection's face gradient (PoissonOperator_FaceGradientFlux()) applied to
 * `P`. Every interior face is set, and on a periodic axis also the seam face at index 0.
 */
static PetscErrorCode RHSFaceTermsWithPressureGradient(UserCtx *user, Vec Rct, Vec Rhs)
{
    const DMDALocalInfo *info = &user->info;
    const PetscInt       m[3] = {info->mx, info->my, info->mz};
    const PetscInt       start[3] = {info->xs, info->ys, info->zs};
    PetscBool            periodic[3];
    PetscInt             interior_lo[3], interior_hi[3];
    PoissonFaceMetrics   metrics;
    const PetscReal   ***nvert, ***p;
    const Cmpnts      ***rct;
    Cmpnts            ***rhs;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    RHSInteriorBounds(info, interior_lo, interior_hi);
    PoissonOperator_PeriodicAxes(user, periodic);
    PetscCall(PoissonOperator_GetFaceMetrics(user, &metrics));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecGetArrayRead(user->da, user->lP, (void *)&p));
    PetscCall(DMDAVecGetArrayRead(user->fda, Rct, (void *)&rct));
    PetscCall(DMDAVecGetArray(user->fda, Rhs, &rhs));
    for (PetscInt n = 0; n < 3; n++) {
        PetscInt lo[3] = {interior_lo[0], interior_lo[1], interior_lo[2]};
        PetscInt across[3] = {0, 0, 0};

        if (periodic[n] && start[n] == 0) lo[n] = 0;
        across[n] = 1;
        for (PetscInt k = lo[2]; k < interior_hi[2]; k++) {
            for (PetscInt j = lo[1]; j < interior_hi[1]; j++) {
                for (PetscInt i = lo[0]; i < interior_hi[0]; i++) {
                    const PetscInt   c[3] = {i, j, k};
                    const PetscReal *own = &rct[k][j][i].x;
                    const PetscReal *next = &rct[k + across[2]][j + across[1]][i + across[0]].x;

                    (&rhs[k][j][i].x)[n] = 0.5 * (own[n] + next[n]) -
                        PoissonOperator_FaceGradientFlux(&metrics, nvert, p, c, n, m, periodic);
                }
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, Rhs, &rhs));
    PetscCall(DMDAVecRestoreArrayRead(user->fda, Rct, (void *)&rct));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lP, (void *)&p));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(PoissonOperator_RestoreFaceMetrics(user, &metrics));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Zeroes the residual rows that carry no unknown: the component a 2D run disables, every
 *        component of a solid cell, and each face component whose far cell is solid.
 */
static PetscErrorCode RHSZeroInactiveRows(UserCtx *user, Vec Rhs)
{
    const PetscInt     two_d = user->simCtx->TwoD;
    PetscInt           lo[3], hi[3];
    const PetscReal ***nvert;
    Cmpnts          ***rhs;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    RHSInteriorBounds(&user->info, lo, hi);
    PetscCall(DMDAVecGetArrayRead(user->da, user->lNvert, (void *)&nvert));
    PetscCall(DMDAVecGetArray(user->fda, Rhs, &rhs));
    for (PetscInt k = lo[2]; k < hi[2]; k++) {
        for (PetscInt j = lo[1]; j < hi[1]; j++) {
            for (PetscInt i = lo[0]; i < hi[0]; i++) {
                if (two_d == 1)      rhs[k][j][i].x = 0.0;
                else if (two_d == 2) rhs[k][j][i].y = 0.0;
                else if (two_d == 3) rhs[k][j][i].z = 0.0;
                if (nvert[k][j][i] > 0.1) {
                    rhs[k][j][i].x = 0.0; rhs[k][j][i].y = 0.0; rhs[k][j][i].z = 0.0;
                }
                if (nvert[k][j][i+1] > 0.1) rhs[k][j][i].x = 0.0;
                if (nvert[k][j+1][i] > 0.1) rhs[k][j][i].y = 0.0;
                if (nvert[k+1][j][i] > 0.1) rhs[k][j][i].z = 0.0;
            }
        }
    }
    PetscCall(DMDAVecRestoreArray(user->fda, Rhs, &rhs));
    PetscCall(DMDAVecRestoreArrayRead(user->da, user->lNvert, (void *)&nvert));
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}

/**
 * @brief Implementation of \ref ComputeRHS().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/rhs.h`.
 * @see ComputeRHS()
 */
PetscErrorCode ComputeRHS(UserCtx *user, Vec Rhs)
{
    const FieldId cell_fields[] = {FIELD_ID_UCAT};
    Vec           Conv, Visc, Rct;

    PetscFunctionBeginUser;
    PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d, Block %" PetscInt_FMT ": Computing RHS (FormFunction1)...\n",
              user->simCtx->rank, user->_this);
    /* Work vectors from fda's pool, zeroed as fresh vectors would be: the residual is evaluated
       many times per step, and allocating them each time costs more than the zeroing. */
    PetscCall(DMGetLocalVector(user->fda, &Conv)); PetscCall(VecZeroEntries(Conv));
    PetscCall(DMGetLocalVector(user->fda, &Visc)); PetscCall(VecZeroEntries(Visc));
    PetscCall(DMGetLocalVector(user->fda, &Rct));  PetscCall(VecZeroEntries(Rct));

    /* Cartesian velocity, with its periodic copies and ghosts, from the contravariant flux. */
    PetscCall(Contra2Cart(user));
    PetscCall(SynchronizePeriodicCellFields(user, 1, cell_fields));
    PetscCall(UpdateLocalGhosts(user, FIELD_ID_UCAT));

    /* Spatial terms in Cartesian components at cell centres, then contravariant. */
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Calculating convective terms...\n");
    PetscCall(ComputeConvection(user, user->lUcont, user->lUcat, Conv));
    if (!user->simCtx->invicid) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "  Calculating viscous terms...\n");
        PetscCall(Viscous(user, user->lUcont, user->lUcat, Visc));
    }
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Converting Cartesian RHS to Contravariant RHS...\n");
    PetscCall(RHSContravariantSpatialTerms(user, Conv, Visc, Rct));
    PetscCall(ComputeBodyForces(user, Rct));
    PetscCall(SynchronizePeriodicLocalStaggeredField(user, Rct));

    /* Face residual: averaged spatial terms less the pressure gradient, then inactive rows. */
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Adding pressure gradient term to RHS...\n");
    PetscCall(RHSFaceTermsWithPressureGradient(user, Rct, Rhs));
    PetscCall(RHSZeroInactiveRows(user, Rhs));
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Final cleanup and edge-cases complete .\n");

    PetscCall(DMRestoreLocalVector(user->fda, &Rct));
    PetscCall(DMRestoreLocalVector(user->fda, &Visc));
    PetscCall(DMRestoreLocalVector(user->fda, &Conv));
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d, Block %" PetscInt_FMT ": RHS computation complete.\n",
              user->simCtx->rank, user->_this);
    PROFILE_FUNCTION_END;
    PetscFunctionReturn(0);
}
