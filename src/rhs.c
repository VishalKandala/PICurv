#include "rhs.h"

#include "solvers.h" // Or a new "rhs.h"
#include "logging.h"

#include "rhs.h"
#include "logging.h"

// Forward declarations for the legacy helper functions if they are not in a header.
// It is best practice to move these to rhs.h if they are in rhs.c.

#undef __FUNCT__
#define __FUNCT__ "CalculateCovariantMetrics"
PetscErrorCode CalculateCovariantMetrics(double g[3][3], double G[3][3])
{
	PetscFunctionBeginUser;
	PROFILE_FUNCTION_BEGIN;
	/*
		| csi.x  csi.y csi.z |-1		| x.csi  x.eta x.zet | 
		| eta.x eta.y eta.z |	 =	| y.csi   y.eta  y.zet |
		| zet.x zet.y zet.z |		| z.csi  z.eta z.zet |
	
	*/
	const double a11=g[0][0], a12=g[0][1], a13=g[0][2];
	const double a21=g[1][0], a22=g[1][1], a23=g[1][2];
	const double a31=g[2][0], a32=g[2][1], a33=g[2][2];

	double det= a11*(a33*a22-a32*a23)-a21*(a33*a12-a32*a13)+a31*(a23*a12-a22*a13);
	
	G[0][0] = (a33*a22-a32*a23)/det,	G[0][1] = - (a33*a12-a32*a13)/det, 	G[0][2] = (a23*a12-a22*a13)/det;
	G[1][0] = -(a33*a21-a31*a23)/det, G[1][1] = (a33*a11-a31*a13)/det,	G[1][2] = - (a23*a11-a21*a13)/det;
	G[2][0] = (a32*a21-a31*a22)/det,	G[2][1] = - (a32*a11-a31*a12)/det,	G[2][2] = (a22*a11-a21*a12)/det;
	
	PROFILE_FUNCTION_END;
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "CalculateNormalAndArea"
PetscErrorCode CalculateNormalAndArea(Cmpnts csi, Cmpnts eta, Cmpnts zet, double ni[3], double nj[3], double nk[3], double *Ai, double *Aj, double *Ak)
{
	PetscFunctionBeginUser;
	PROFILE_FUNCTION_BEGIN;
	double g[3][3];
	double G[3][3];
	
	g[0][0]=csi.x, g[0][1]=csi.y, g[0][2]=csi.z;
	g[1][0]=eta.x, g[1][1]=eta.y, g[1][2]=eta.z;
	g[2][0]=zet.x, g[2][1]=zet.y, g[2][2]=zet.z;
	
	CalculateCovariantMetrics(g, G);
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

#undef __FUNCT__
#define __FUNCT__ "Calculatedxdydz"
PetscErrorCode Calculatedxdydz(PetscReal ajc, Cmpnts csi, Cmpnts eta, Cmpnts zet, double *dx, double *dy, double *dz)
{
		PetscFunctionBeginUser;
		PROFILE_FUNCTION_BEGIN;
        double ni[3], nj[3], nk[3];
        double Li, Lj, Lk;
        double Ai, Aj, Ak;
        double vol = 1./ajc;

        CalculateNormalAndArea(csi, eta, zet, ni, nj, nk, &Ai, &Aj, &Ak);
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

#undef __FUNCT__
#define __FUNCT__ "Convection"
PetscErrorCode Convection(UserCtx *user, Vec Ucont, Vec Ucat, Vec Conv)
{

  // --- CONTEXT ACQUISITION BLOCK ---
  // Get the master simulation context from the UserCtx.
  SimCtx *simCtx = user->simCtx;

  // Create local variables to mirror the legacy globals for minimal code changes.
  const PetscInt les = simCtx->les;
  const PetscInt central = simCtx->central; // Get this from SimCtx now
  // --- END CONTEXT ACQUISITION BLOCK ---
 
  Cmpnts	***ucont, ***ucat;
  DM		da = user->da, fda = user->fda;
  DMDALocalInfo	info;
  PetscInt	xs, xe, ys, ye, zs, ze; // Local grid information
  PetscInt	mx, my, mz; // Dimensions in three directions
  PetscInt	i, j, k;
  Vec		Fp1, Fp2, Fp3;
  Cmpnts	***fp1, ***fp2, ***fp3;
  Cmpnts	***conv;
 
  PetscReal	ucon, up, um;
  PetscReal	coef = 0.125, innerblank=7.;

  PetscInt	lxs, lxe, lys, lye, lzs, lze, gxs, gxe, gys, gye, gzs,gze;

  PetscReal	***nvert,***aj;

  PROFILE_FUNCTION_BEGIN;

  DMDAGetLocalInfo(da, &info);
  mx = info.mx; my = info.my; mz = info.mz;
  xs = info.xs; xe = xs + info.xm;
  ys = info.ys; ye = ys + info.ym;
  zs = info.zs; ze = zs + info.zm;
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  DMDAVecGetArray(fda, Ucont, &ucont);
  DMDAVecGetArray(fda, Ucat,  &ucat);
  DMDAVecGetArray(fda, Conv,  &conv);
  DMDAVecGetArray(da, user->lAj, &aj);

  VecDuplicate(Ucont, &Fp1);
  VecDuplicate(Ucont, &Fp2);
  VecDuplicate(Ucont, &Fp3);

  DMDAVecGetArray(fda, Fp1, &fp1);
  DMDAVecGetArray(fda, Fp2, &fp2);
  DMDAVecGetArray(fda, Fp3, &fp3);

  DMDAVecGetArray(da, user->lNvert, &nvert);
  

  /* We have two different sets of node: 1. grid node, the physical points
     where grid lines intercross; 2. storage node, where we store variables.
     All node without explicitly specified as "grid node" refers to
     storage node.

     The integer node is defined at cell center while half node refers to
     the actual grid node. (The reason to choose this arrangement is we need
     ghost node, which is half node away from boundaries, to specify boundary
     conditions. By using this storage arrangement, the actual storage need
     is (IM+1) * (JM + 1) * (KM+1) where IM, JM, & KM refer to the number of
     grid nodes along i, j, k directions.)

     DA, the data structure used to define the storage of 3D arrays, is defined
     as mx * my * mz. mx = IM+1, my = JM+1, mz = KM+1.

     Staggered grid arrangement is used in this solver.
     Pressure is stored at interger node (hence the cell center) and volume
     fluxes defined on the center of each surface of a given control volume
     is stored on the cloest upper integer node. */

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
  
  //Mohsen Sep 2012//
/* First update the computational ghost points velocity for periodic boundary conditions
 just for this subroutine because of Quick scheme for velocity deravatives */
 /*  for (k=zs; k<ze; k++) { */
/*     for (j=ys; j<ye; j++) { */
/*       for (i=xs; i<xe; i++) { */
/* 	if (i==1 && (j==1) && (k==0 || k==1)) */
/* 	  PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d u is %.15le v is %.15le  w is %.15le \n",i,j,k,ucat[k][j][i].x,ucat[k][j][i].y,ucat[k][j][i].z ); */
/*       } */
/*     } */
/*   } */
  if (user->bctype[0]==7 || user->bctype[1]==7){
    if (xs==0){
      i=xs-1;
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  ucat[k][j][i].x=ucat[k][j][i-2].x;
	  ucat[k][j][i].y=ucat[k][j][i-2].y;
	  ucat[k][j][i].z=ucat[k][j][i-2].z;
	  nvert[k][j][i]=nvert[k][j][i-2];
	}
      }
    }
    if (xe==mx){
      i=mx;
      for (k=gzs; k<gze; k++) {
	for (j=gys; j<gye; j++) {
	  ucat[k][j][i].x=ucat[k][j][i+2].x;
	  ucat[k][j][i].y=ucat[k][j][i+2].y;
	  ucat[k][j][i].z=ucat[k][j][i+2].z;
	  nvert[k][j][i]=nvert[k][j][i+2];
	}
      }
    }
  }  
  if (user->bctype[2]==7 || user->bctype[3]==7){
    if (ys==0){
      j=ys-1;
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i].x=ucat[k][j-2][i].x;
	  ucat[k][j][i].y=ucat[k][j-2][i].y;
	  ucat[k][j][i].z=ucat[k][j-2][i].z;
	  nvert[k][j][i]=nvert[k][j-2][i];
	}
      }
    }
    if (ye==my){
      j=my;
      for (k=gzs; k<gze; k++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i].x=ucat[k][j+2][i].x;
	  ucat[k][j][i].y=ucat[k][j+2][i].y;
	  ucat[k][j][i].z=ucat[k][j+2][i].z;
	  nvert[k][j][i]=nvert[k][j+2][i];
	}
      }
    }
  }
  if (user->bctype[4]==7 || user->bctype[5]==7){
    if (zs==0){
      k=zs-1;
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i].x=ucat[k-2][j][i].x;
	  ucat[k][j][i].y=ucat[k-2][j][i].y;
	  ucat[k][j][i].z=ucat[k-2][j][i].z;
	  nvert[k][j][i]=nvert[k-2][j][i];
	}
      }
    }
    if (ze==mz){
      k=mz;
      for (j=gys; j<gye; j++) {
	for (i=gxs; i<gxe; i++) {
	  ucat[k][j][i].x=ucat[k+2][j][i].x;
	  ucat[k][j][i].y=ucat[k+2][j][i].y;
	  ucat[k][j][i].z=ucat[k+2][j][i].z;
	  nvert[k][j][i]=nvert[k+2][j][i];
	}
      }
    }
  }

  VecSet(Conv, 0.0);
 
 /* Calculating the convective terms on cell centers.
    First calcualte the contribution from i direction
    The flux is evaluated by QUICK scheme */

  for (k=lzs; k<lze; k++){
    for (j=lys; j<lye; j++){
      for (i=lxs-1; i<lxe; i++){
       
       
	ucon = ucont[k][j][i].x * 0.5;
       
	up = ucon + fabs(ucon);
	um = ucon - fabs(ucon);
	
	if (i>0  && i<mx-2 &&
		 (nvert[k][j][i+1] < 0.1 || nvert[k][j][i+1]>innerblank) &&
		 (nvert[k][j][i-1] < 0.1 || nvert[k][j][i-1]>innerblank)) { // interial nodes
	  if ((les || central)) {
	    fp1[k][j][i].x =  ucon * ( ucat[k][j][i].x + ucat[k][j][i+1].x );
	    fp1[k][j][i].y =  ucon * ( ucat[k][j][i].y + ucat[k][j][i+1].y );
	    fp1[k][j][i].z =  ucon * ( ucat[k][j][i].z + ucat[k][j][i+1].z );

	  } else {
	  fp1[k][j][i].x =  
	    um * (coef * (-ucat[k][j][i+2].x -2.* ucat[k][j][i+1].x +3.* ucat[k][j][i  ].x) +ucat[k][j][i+1].x) +
	    up * (coef * (-ucat[k][j][i-1].x -2.* ucat[k][j][i  ].x +3.* ucat[k][j][i+1].x) +ucat[k][j][i  ].x);
	  fp1[k][j][i].y = 
	    um * (coef * (-ucat[k][j][i+2].y -2.* ucat[k][j][i+1].y +3.* ucat[k][j][i  ].y) +ucat[k][j][i+1].y) +
	    up * (coef * (-ucat[k][j][i-1].y -2.* ucat[k][j][i  ].y +3.* ucat[k][j][i+1].y) +ucat[k][j][i  ].y);
	  fp1[k][j][i].z = 
	    um * (coef * (-ucat[k][j][i+2].z -2.* ucat[k][j][i+1].z +3.* ucat[k][j][i  ].z) +ucat[k][j][i+1].z) +
	    up * (coef * (-ucat[k][j][i-1].z -2.* ucat[k][j][i  ].z +3.* ucat[k][j][i+1].z) +ucat[k][j][i  ].z);
	  }
	}
	else  if ((les || central) && (i==0 || i==mx-2) &&
		  (nvert[k][j][i+1] < 0.1 || nvert[k][j][i+1]>innerblank) &&
		  (nvert[k][j][i  ] < 0.1 || nvert[k][j][i  ]>innerblank)) 
	  {
	    fp1[k][j][i].x =  ucon * ( ucat[k][j][i].x + ucat[k][j][i+1].x );
	    fp1[k][j][i].y =  ucon * ( ucat[k][j][i].y + ucat[k][j][i+1].y );
	    fp1[k][j][i].z =  ucon * ( ucat[k][j][i].z + ucat[k][j][i+1].z );
	  }
	else if (i==0 ||(nvert[k][j][i-1] > 0.1) ) {
	  if (user->bctype[0]==7 && i==0 && (nvert[k][j][i-1]<0.1 && nvert[k][j][i+1]<0.1)){//Mohsen Feb 12 
	    fp1[k][j][i].x =  
	      um * (coef * (-ucat[k][j][i+2].x -2.* ucat[k][j][i+1].x +3.* ucat[k][j][i  ].x) +ucat[k][j][i+1].x) +
	      up * (coef * (-ucat[k][j][i-1].x -2.* ucat[k][j][i  ].x +3.* ucat[k][j][i+1].x) +ucat[k][j][i  ].x);
	    fp1[k][j][i].y = 
	      um * (coef * (-ucat[k][j][i+2].y -2.* ucat[k][j][i+1].y +3.* ucat[k][j][i  ].y) +ucat[k][j][i+1].y) +
	      up * (coef * (-ucat[k][j][i-1].y -2.* ucat[k][j][i  ].y +3.* ucat[k][j][i+1].y) +ucat[k][j][i  ].y);
	    fp1[k][j][i].z = 
	      um * (coef * (-ucat[k][j][i+2].z -2.* ucat[k][j][i+1].z +3.* ucat[k][j][i  ].z) +ucat[k][j][i+1].z) +
	      up * (coef * (-ucat[k][j][i-1].z -2.* ucat[k][j][i  ].z +3.* ucat[k][j][i+1].z) +ucat[k][j][i  ].z);
	  }else{
	    fp1[k][j][i].x = 
	      um * (coef * (-ucat[k][j][i+2].x -2.* ucat[k][j][i+1].x +3.* ucat[k][j][i  ].x) +ucat[k][j][i+1].x) +
	      up * (coef * (-ucat[k][j][i  ].x -2.* ucat[k][j][i  ].x +3.* ucat[k][j][i+1].x) +ucat[k][j][i  ].x);
	    fp1[k][j][i].y = 
	      um * (coef * (-ucat[k][j][i+2].y -2.* ucat[k][j][i+1].y +3.* ucat[k][j][i  ].y) +ucat[k][j][i+1].y) +
	      up * (coef * (-ucat[k][j][i  ].y -2.* ucat[k][j][i  ].y +3.* ucat[k][j][i+1].y) +ucat[k][j][i  ].y);
	    fp1[k][j][i].z = 
	      um * (coef * (-ucat[k][j][i+2].z -2.* ucat[k][j][i+1].z +3.* ucat[k][j][i  ].z) +ucat[k][j][i+1].z) +
	      up * (coef * (-ucat[k][j][i  ].z -2.* ucat[k][j][i  ].z +3.* ucat[k][j][i+1].z) +ucat[k][j][i  ].z);	  
	  }
	}
	else if (i==mx-2 ||(nvert[k][j][i+1]) > 0.1) {
	  if (user->bctype[0]==7 && i==mx-2 &&(nvert[k][j][i-1]<0.1 && nvert[k][j][i+1]<0.1)){//Mohsen Feb 12 
	    fp1[k][j][i].x =  
	    um * (coef * (-ucat[k][j][i+2].x -2.* ucat[k][j][i+1].x +3.* ucat[k][j][i  ].x) +ucat[k][j][i+1].x) +
	    up * (coef * (-ucat[k][j][i-1].x -2.* ucat[k][j][i  ].x +3.* ucat[k][j][i+1].x) +ucat[k][j][i  ].x);
	  fp1[k][j][i].y = 
	    um * (coef * (-ucat[k][j][i+2].y -2.* ucat[k][j][i+1].y +3.* ucat[k][j][i  ].y) +ucat[k][j][i+1].y) +
	    up * (coef * (-ucat[k][j][i-1].y -2.* ucat[k][j][i  ].y +3.* ucat[k][j][i+1].y) +ucat[k][j][i  ].y);
	  fp1[k][j][i].z = 
	    um * (coef * (-ucat[k][j][i+2].z -2.* ucat[k][j][i+1].z +3.* ucat[k][j][i  ].z) +ucat[k][j][i+1].z) +
	    up * (coef * (-ucat[k][j][i-1].z -2.* ucat[k][j][i  ].z +3.* ucat[k][j][i+1].z) +ucat[k][j][i  ].z);
	  }else{
	  fp1[k][j][i].x = 
	    um * (coef * (-ucat[k][j][i+1].x -2. * ucat[k][j][i+1].x +3. * ucat[k][j][i  ].x) +ucat[k][j][i+1].x) +
	    up * (coef * (-ucat[k][j][i-1].x -2. * ucat[k][j][i  ].x +3. * ucat[k][j][i+1].x) +ucat[k][j][i  ].x);
	  fp1[k][j][i].y = 
	    um * (coef * (-ucat[k][j][i+1].y -2. * ucat[k][j][i+1].y +3. * ucat[k][j][i  ].y) +ucat[k][j][i+1].y) +
	    up * (coef * (-ucat[k][j][i-1].y -2. * ucat[k][j][i  ].y +3. * ucat[k][j][i+1].y) +ucat[k][j][i  ].y);
	  fp1[k][j][i].z = 
	    um * (coef * (-ucat[k][j][i+1].z -2. * ucat[k][j][i+1].z +3. * ucat[k][j][i  ].z) +ucat[k][j][i+1].z) +
	    up * (coef * (-ucat[k][j][i-1].z -2. * ucat[k][j][i  ].z +3. * ucat[k][j][i+1].z) +ucat[k][j][i  ].z);
	  }
	}
      }
    }
  }

  /* j direction */
  for (k=lzs; k<lze; k++) {
    for(j=lys-1; j<lye; j++) {
      for(i=lxs; i<lxe; i++) {
	ucon = ucont[k][j][i].y * 0.5;

	up = ucon + fabs(ucon);
	um = ucon - fabs(ucon);

	if (j>0 && j<my-2 &&
	    (nvert[k][j+1][i] < 0.1 || nvert[k][j+1][i] > innerblank) &&
	    (nvert[k][j-1][i] < 0.1 || nvert[k][j-1][i] > innerblank))  {
	  if ((les || central)) {
	    fp2[k][j][i].x =  ucon * ( ucat[k][j][i].x + ucat[k][j+1][i].x );
	    fp2[k][j][i].y =  ucon * ( ucat[k][j][i].y + ucat[k][j+1][i].y );
	    fp2[k][j][i].z =  ucon * ( ucat[k][j][i].z + ucat[k][j+1][i].z );

	  } else {
	  fp2[k][j][i].x = 
	    um * (coef * (-ucat[k][j+2][i].x -2. * ucat[k][j+1][i].x +3. * ucat[k][j  ][i].x) +ucat[k][j+1][i].x) +		     
	    up * (coef * (-ucat[k][j-1][i].x -2. * ucat[k][j  ][i].x +3. * ucat[k][j+1][i].x) +ucat[k][j  ][i].x);		     
	  fp2[k][j][i].y = 				     
	    um * (coef * (-ucat[k][j+2][i].y -2. * ucat[k][j+1][i].y +3. * ucat[k][j  ][i].y) +ucat[k][j+1][i].y) +		     
	    up * (coef * (-ucat[k][j-1][i].y -2. * ucat[k][j  ][i].y +3. * ucat[k][j+1][i].y) +ucat[k][j  ][i].y);		     
	  fp2[k][j][i].z = 				     
	    um * (coef * (-ucat[k][j+2][i].z -2. * ucat[k][j+1][i].z +3. * ucat[k][j  ][i].z) +ucat[k][j+1][i].z) +		     
	    up * (coef * (-ucat[k][j-1][i].z -2. * ucat[k][j  ][i].z +3. * ucat[k][j+1][i].z) +ucat[k][j  ][i].z);
	  }
	}
	else  if ((les || central) && (j==0 || i==my-2) &&
		  (nvert[k][j+1][i] < 0.1 || nvert[k][j+1][i]>innerblank) &&
		  (nvert[k][j  ][i] < 0.1 || nvert[k][j  ][i]>innerblank)) 
	  {
	    fp2[k][j][i].x =  ucon * ( ucat[k][j][i].x + ucat[k][j+1][i].x );
	    fp2[k][j][i].y =  ucon * ( ucat[k][j][i].y + ucat[k][j+1][i].y );
	    fp2[k][j][i].z =  ucon * ( ucat[k][j][i].z + ucat[k][j+1][i].z );
	  }
	else if (j==0 || (nvert[k][j-1][i]) > 0.1) {
	  if (user->bctype[2]==7 && j==0 && (nvert[k][j-1][i]<0.1 && nvert[k][j+1][i]<0.1 )){//Mohsen Feb 12 //
	    fp2[k][j][i].x = 
	    um * (coef * (-ucat[k][j+2][i].x -2. * ucat[k][j+1][i].x +3. * ucat[k][j  ][i].x) +ucat[k][j+1][i].x) +		     
	    up * (coef * (-ucat[k][j-1][i].x -2. * ucat[k][j  ][i].x +3. * ucat[k][j+1][i].x) +ucat[k][j  ][i].x);		     
	  fp2[k][j][i].y = 				     
	    um * (coef * (-ucat[k][j+2][i].y -2. * ucat[k][j+1][i].y +3. * ucat[k][j  ][i].y) +ucat[k][j+1][i].y) +		     
	    up * (coef * (-ucat[k][j-1][i].y -2. * ucat[k][j  ][i].y +3. * ucat[k][j+1][i].y) +ucat[k][j  ][i].y);		     
	  fp2[k][j][i].z = 				     
	    um * (coef * (-ucat[k][j+2][i].z -2. * ucat[k][j+1][i].z +3. * ucat[k][j  ][i].z) +ucat[k][j+1][i].z) +		     
	    up * (coef * (-ucat[k][j-1][i].z -2. * ucat[k][j  ][i].z +3. * ucat[k][j+1][i].z) +ucat[k][j  ][i].z);
	  }else{
	  fp2[k][j][i].x = 
	    um * (coef * (-ucat[k][j+2][i].x -2. * ucat[k][j+1][i].x +3. * ucat[k][j  ][i].x) +ucat[k][j+1][i].x) +		     
	    up * (coef * (-ucat[k][j  ][i].x -2. * ucat[k][j  ][i].x +3. * ucat[k][j+1][i].x) +ucat[k][j  ][i].x);		     
	  fp2[k][j][i].y = 				     
	    um * (coef * (-ucat[k][j+2][i].y -2. * ucat[k][j+1][i].y +3. * ucat[k][j  ][i].y) +ucat[k][j+1][i].y) +		     
	    up * (coef * (-ucat[k][j  ][i].y -2. * ucat[k][j  ][i].y +3. * ucat[k][j+1][i].y) +ucat[k][j  ][i].y);		     
	  fp2[k][j][i].z = 				     
	    um * (coef * (-ucat[k][j+2][i].z -2. * ucat[k][j+1][i].z +3. * ucat[k][j  ][i].z) +ucat[k][j+1][i].z) +		     
	    up * (coef * (-ucat[k][j  ][i].z -2. * ucat[k][j  ][i].z +3. * ucat[k][j+1][i].z) +ucat[k][j  ][i].z);
	}
	}
	else if (j==my-2 ||(nvert[k][j+1][i]) > 0.1) {
	  if (user->bctype[2]==7 && j==my-2 && (nvert[k][j-1][i]<0.1 && nvert[k][j+1][i]<0.1 )){//Mohsen Feb 12// 
	    fp2[k][j][i].x = 
	    um * (coef * (-ucat[k][j+2][i].x -2. * ucat[k][j+1][i].x +3. * ucat[k][j  ][i].x) +ucat[k][j+1][i].x) +		     
	    up * (coef * (-ucat[k][j-1][i].x -2. * ucat[k][j  ][i].x +3. * ucat[k][j+1][i].x) +ucat[k][j  ][i].x);		     
	  fp2[k][j][i].y = 				     
	    um * (coef * (-ucat[k][j+2][i].y -2. * ucat[k][j+1][i].y +3. * ucat[k][j  ][i].y) +ucat[k][j+1][i].y) +		     
	    up * (coef * (-ucat[k][j-1][i].y -2. * ucat[k][j  ][i].y +3. * ucat[k][j+1][i].y) +ucat[k][j  ][i].y);		     
	  fp2[k][j][i].z = 				     
	    um * (coef * (-ucat[k][j+2][i].z -2. * ucat[k][j+1][i].z +3. * ucat[k][j  ][i].z) +ucat[k][j+1][i].z) +		     
	    up * (coef * (-ucat[k][j-1][i].z -2. * ucat[k][j  ][i].z +3. * ucat[k][j+1][i].z) +ucat[k][j  ][i].z);
	  }else{
	  fp2[k][j][i].x = 
	    um * (coef * (-ucat[k][j+1][i].x -2. * ucat[k][j+1][i].x +3. * ucat[k][j  ][i].x) +ucat[k][j+1][i].x) +		     
	    up * (coef * (-ucat[k][j-1][i].x -2. * ucat[k][j  ][i].x +3. * ucat[k][j+1][i].x) +ucat[k][j  ][i].x);		     
	  fp2[k][j][i].y = 				     
	    um * (coef * (-ucat[k][j+1][i].y -2. * ucat[k][j+1][i].y +3. * ucat[k][j  ][i].y) +ucat[k][j+1][i].y) +		     
	    up * (coef * (-ucat[k][j-1][i].y -2. * ucat[k][j  ][i].y +3. * ucat[k][j+1][i].y) +ucat[k][j][i  ].y);		     
	  fp2[k][j][i].z = 				     
	    um * (coef * (-ucat[k][j+1][i].z -2. * ucat[k][j+1][i].z +3. * ucat[k][j  ][i].z) +ucat[k][j+1][i].z) +		     
	    up * (coef * (-ucat[k][j-1][i].z -2. * ucat[k][j  ][i].z +3. * ucat[k][j+1][i].z) +ucat[k][j][i  ].z);
	  }
	}
      }
    }
  }
 

  /* k direction */
  for (k=lzs-1; k<lze; k++) {
    for(j=lys; j<lye; j++) {
      for(i=lxs; i<lxe; i++) {
	ucon = ucont[k][j][i].z * 0.5;
       
	up = ucon + fabs(ucon);
	um = ucon - fabs(ucon);

	if (k>0 && k<mz-2 &&
	    (nvert[k+1][j][i] < 0.1 || nvert[k+1][j][i] > innerblank) &&
	    (nvert[k-1][j][i] < 0.1 || nvert[k-1][j][i] > innerblank)) {
	  if ((les || central)) {
	    fp3[k][j][i].x =  ucon * ( ucat[k][j][i].x + ucat[k+1][j][i].x );
	    fp3[k][j][i].y =  ucon * ( ucat[k][j][i].y + ucat[k+1][j][i].y );
	    fp3[k][j][i].z =  ucon * ( ucat[k][j][i].z + ucat[k+1][j][i].z );

	  } else {
	  fp3[k][j][i].x = 
	    um * (coef * (-ucat[k+2][j][i].x -2. * ucat[k+1][j][i].x +3. * ucat[k  ][j][i].x) +ucat[k+1][j][i].x)   +		     
	    up * (coef * (-ucat[k-1][j][i].x -2. * ucat[k  ][j][i].x +3. * ucat[k+1][j][i].x) +ucat[k  ][j][i].x);  		     
	  fp3[k][j][i].y = 	       			     
	    um * (coef * (-ucat[k+2][j][i].y -2. * ucat[k+1][j][i].y +3. * ucat[k  ][j][i].y) +ucat[k+1][j][i].y)   +		     
	    up * (coef * (-ucat[k-1][j][i].y -2. * ucat[k  ][j][i].y +3. * ucat[k+1][j][i].y) +ucat[k  ][j][i].y);  		     
	  fp3[k][j][i].z = 	       			     
	    um * (coef * (-ucat[k+2][j][i].z -2. * ucat[k+1][j][i].z +3. * ucat[k  ][j][i].z) +ucat[k+1][j][i].z)   +		     
	    up * (coef * (-ucat[k-1][j][i].z -2. * ucat[k  ][j][i].z +3. * ucat[k+1][j][i].z) +ucat[k  ][j][i].z);
	  }
	}
	else if ((les || central) && (k==0 || k==mz-2) &&
		 (nvert[k+1][j][i] < 0.1 || nvert[k+1][j][i]>innerblank) &&
		 (nvert[k  ][j][i] < 0.1 || nvert[k  ][j][i]>innerblank)) 
	  {
	    fp3[k][j][i].x =  ucon * ( ucat[k][j][i].x + ucat[k+1][j][i].x );
	    fp3[k][j][i].y =  ucon * ( ucat[k][j][i].y + ucat[k+1][j][i].y );
	    fp3[k][j][i].z =  ucon * ( ucat[k][j][i].z + ucat[k+1][j][i].z );
	  }
	else if (k<mz-2 && (k==0 ||(nvert[k-1][j][i]) > 0.1)) {
	  if(user->bctype[4]==7 && k==0 && (nvert[k-1][j][i]<0.1 && nvert[k+1][j][i]<0.1)){//Mohsen Feb 12//
	    fp3[k][j][i].x = 
	    um * (coef * (-ucat[k+2][j][i].x -2. * ucat[k+1][j][i].x +3. * ucat[k  ][j][i].x) +ucat[k+1][j][i].x)   +		     
	    up * (coef * (-ucat[k-1][j][i].x -2. * ucat[k  ][j][i].x +3. * ucat[k+1][j][i].x) +ucat[k  ][j][i].x);  		     
	  fp3[k][j][i].y = 	       			     
	    um * (coef * (-ucat[k+2][j][i].y -2. * ucat[k+1][j][i].y +3. * ucat[k  ][j][i].y) +ucat[k+1][j][i].y)   +		     
	    up * (coef * (-ucat[k-1][j][i].y -2. * ucat[k  ][j][i].y +3. * ucat[k+1][j][i].y) +ucat[k  ][j][i].y);  		     
	  fp3[k][j][i].z = 	       			     
	    um * (coef * (-ucat[k+2][j][i].z -2. * ucat[k+1][j][i].z +3. * ucat[k  ][j][i].z) +ucat[k+1][j][i].z)   +		     
	    up * (coef * (-ucat[k-1][j][i].z -2. * ucat[k  ][j][i].z +3. * ucat[k+1][j][i].z) +ucat[k  ][j][i].z);
	  }else{
	  fp3[k][j][i].x = 
	    um * (coef * (-ucat[k+2][j][i].x -2. * ucat[k+1][j][i].x +3. * ucat[k  ][j][i].x) +ucat[k+1][j][i].x)   +		     
	    up * (coef * (-ucat[k  ][j][i].x -2. * ucat[k  ][j][i].x +3. * ucat[k+1][j][i].x) +ucat[k][j][i  ].x);  		     
	  fp3[k][j][i].y = 	       			     
	    um * (coef * (-ucat[k+2][j][i].y -2. * ucat[k+1][j][i].y +3. * ucat[k  ][j][i].y) +ucat[k+1][j][i].y)   +		     
	    up * (coef * (-ucat[k  ][j][i].y -2. * ucat[k  ][j][i].y +3. * ucat[k+1][j][i].y) +ucat[k][j][i  ].y);  		     
	  fp3[k][j][i].z = 	       			     
	    um * (coef * (-ucat[k+2][j][i].z -2. * ucat[k+1][j][i].z +3. * ucat[k  ][j][i].z) +ucat[k+1][j][i].z)   +		     
	    up * (coef * (-ucat[k  ][j][i].z -2. * ucat[k  ][j][i].z +3. * ucat[k+1][j][i].z) +ucat[k][j][i  ].z);
	  }
	}
	else if (k>0 && (k==mz-2 ||(nvert[k+1][j][i]) > 0.1)) {
	  if (user->bctype[4]==7 && k==mz-2 && (nvert[k-1][j][i]<0.1 && nvert[k+1][j][i]<0.1)){//Mohsen Feb 12//
	    fp3[k][j][i].x = 
	    um * (coef * (-ucat[k+2][j][i].x -2. * ucat[k+1][j][i].x +3. * ucat[k  ][j][i].x) +ucat[k+1][j][i].x)   +		     
	    up * (coef * (-ucat[k-1][j][i].x -2. * ucat[k  ][j][i].x +3. * ucat[k+1][j][i].x) +ucat[k  ][j][i].x);  		     
	  fp3[k][j][i].y = 	       			     
	    um * (coef * (-ucat[k+2][j][i].y -2. * ucat[k+1][j][i].y +3. * ucat[k  ][j][i].y) +ucat[k+1][j][i].y)   +		     
	    up * (coef * (-ucat[k-1][j][i].y -2. * ucat[k  ][j][i].y +3. * ucat[k+1][j][i].y) +ucat[k  ][j][i].y);  		     
	  fp3[k][j][i].z = 	       			     
	    um * (coef * (-ucat[k+2][j][i].z -2. * ucat[k+1][j][i].z +3. * ucat[k  ][j][i].z) +ucat[k+1][j][i].z)   +		     
	    up * (coef * (-ucat[k-1][j][i].z -2. * ucat[k  ][j][i].z +3. * ucat[k+1][j][i].z) +ucat[k  ][j][i].z);
	  }else{
	  fp3[k][j][i].x = 
	    um * (coef * (-ucat[k+1][j][i].x -2. * ucat[k+1][j][i].x +3. * ucat[k  ][j][i].x) +ucat[k+1][j][i].x)   +		     
	    up * (coef * (-ucat[k-1][j][i].x -2. * ucat[k  ][j][i].x +3. * ucat[k+1][j][i].x) +ucat[k][j][i  ].x);  		     
	  fp3[k][j][i].y = 	       			     
	    um * (coef * (-ucat[k+1][j][i].y -2. * ucat[k+1][j][i].y +3. * ucat[k  ][j][i].y) +ucat[k+1][j][i].y)   +		     
	    up * (coef * (-ucat[k-1][j][i].y -2. * ucat[k  ][j][i].y +3. * ucat[k+1][j][i].y) +ucat[k][j][i  ].y);  		     
	  fp3[k][j][i].z = 	       			     
	    um * (coef * (-ucat[k+1][j][i].z -2. * ucat[k+1][j][i].z +3. * ucat[k  ][j][i].z) +ucat[k+1][j][i].z)   +		     
	    up * (coef * (-ucat[k-1][j][i].z -2. * ucat[k  ][j][i].z +3. * ucat[k+1][j][i].z) +ucat[k][j][i  ].z);
	  }
	}
      }
    }
  }
 
  /* Calculate the convective terms under cartesian coordinates */

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	conv[k][j][i].x = 
	  fp1[k][j][i].x - fp1[k][j][i-1].x +
	  fp2[k][j][i].x - fp2[k][j-1][i].x +
	  fp3[k][j][i].x - fp3[k-1][j][i].x;

	conv[k][j][i].y = 
	  fp1[k][j][i].y - fp1[k][j][i-1].y +
	  fp2[k][j][i].y - fp2[k][j-1][i].y +
	  fp3[k][j][i].y - fp3[k-1][j][i].y;

	conv[k][j][i].z =
	  fp1[k][j][i].z - fp1[k][j][i-1].z +
	  fp2[k][j][i].z - fp2[k][j-1][i].z +
	  fp3[k][j][i].z - fp3[k-1][j][i].z;
      }
    }
  } 
  /* for (k=zs; k<ze; k++) { */
/*     for (j=ys; j<ye; j++) { */
/*       for (i=xs; i<xe; i++) { */
/* 	if (i==1 && (j==1) && (k==1 || k==21 || k==22|| k==200)) */
/* 	  PetscPrintf(PETSC_COMM_SELF, "@ i= %d j=%d k=%d conv.y is %.15le  conv.z is %.15le \n",i,j,k,conv[k][j][i].y,conv[k][j][i].z); */
/*       } */
/*     } */
/*   }  */

  DMDAVecRestoreArray(fda, Ucont, &ucont);
  DMDAVecRestoreArray(fda, Ucat,  &ucat);
  DMDAVecRestoreArray(fda, Conv,  &conv);
  DMDAVecRestoreArray(da, user->lAj, &aj);

  DMDAVecRestoreArray(fda, Fp1, &fp1);
  DMDAVecRestoreArray(fda, Fp2, &fp2);
  DMDAVecRestoreArray(fda, Fp3, &fp3);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);

  VecDestroy(&Fp1);
  VecDestroy(&Fp2);
  VecDestroy(&Fp3);


  LOG_ALLOW(GLOBAL,LOG_DEBUG,"Convective term calculated .\n");
  
  PROFILE_FUNCTION_END;
  return (0);
}


#undef __FUNCT__
#define __FUNCT__ "Viscous"
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
  const PetscInt les = simCtx->les;
  const PetscInt rans = simCtx->rans;
  const PetscInt ti = simCtx->step; // Assuming simCtx->step is the new integer time counter
  const	PetscReal ren = simCtx->ren;
  const PetscInt clark = simCtx->clark;
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

  VecDuplicate(Ucont, &Fp1);
  VecDuplicate(Ucont, &Fp2);
  VecDuplicate(Ucont, &Fp3);

  DMDAVecGetArray(fda, Fp1, &fp1);
  DMDAVecGetArray(fda, Fp2, &fp2);
  DMDAVecGetArray(fda, Fp3, &fp3);

  DMDAVecGetArray(da, user->lAj, &aj);

  DMDAGetLocalInfo(da, &info);
  
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
 
  if(les) {
    DMDAVecGetArray(da, user->lNu_t, &lnu_t);
  } else if (rans) {
   
    DMDAVecGetArray(da, user->lNu_t, &lnu_t);
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
	
	if( les || (rans && ti>0) ) {
	  //nu_t = pow( 0.5 * ( sqrt(lnu_t[k][j][i]) + sqrt(lnu_t[k][j][i+1]) ), 2.0) * Sabs;
	  nu_t = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j][i+1]);
	  if ( (user->bctype[0]==1 && i==0) || (user->bctype[1]==1 && i==mx-2) ) nu_t=0;    
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
	  double dc, de, dz;
	  Calculatedxdydz (ajc, csi[k][j][i], eta[k][j][i], zet[k][j][i], &dc, &de, &dz);
	  double dc2=dc*dc, de2=de*de, dz2=dz*dz;
	  
	  double t11 = ( dudc * dudc * dc2 + dude * dude * de2 + dudz * dudz * dz2 );
	  double t12 = ( dudc * dvdc * dc2 + dude * dvde * de2 + dudz * dvdz * dz2 );
	  double t13 = ( dudc * dwdc * dc2 + dude * dwde * de2 + dudz * dwdz * dz2 );
	  double t21 = t12;
	  double t22 = ( dvdc * dvdc * dc2 + dvde * dvde * de2 + dvdz * dvdz * dz2 );
	  double t23 = ( dvdc * dwdc * dc2 + dvde * dwde * de2 + dvdz * dwdz * dz2 );
	  double t31 = t13;
	  double t32 = t23;
	  double t33 = ( dwdc * dwdc * dc2 + dwde * dwde * de2 + dwdz * dwdz * dz2 );
	  
	  fp1[k][j][i].x -= ( t11 * csi0 + t12 * csi1 + t13 * csi2 ) / 12.;
	  fp1[k][j][i].y -= ( t21 * csi0 + t22 * csi1 + t23 * csi2 ) / 12.;
	  fp1[k][j][i].z -= ( t31 * csi0 + t32 * csi1 + t33 * csi2 ) / 12.;
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
		
	if( les || (rans && ti>0) ) {
	  //nu_t = pow( 0.5 * ( sqrt(lnu_t[k][j][i]) + sqrt(lnu_t[k][j+1][i]) ), 2.0) * Sabs;
	  nu_t = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j+1][i]);
	  if ( (user->bctype[2]==1 && j==0) || (user->bctype[3]==1 && j==my-2) ) nu_t=0;
		
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
	  double dc, de, dz;
	  Calculatedxdydz(ajc, csi[k][j][i], eta[k][j][i], zet[k][j][i], &dc, &de, &dz);
	  double dc2=dc*dc, de2=de*de, dz2=dz*dz;
			
	  double t11 = ( dudc * dudc * dc2 + dude * dude * de2 + dudz * dudz * dz2 );
	  double t12 = ( dudc * dvdc * dc2 + dude * dvde * de2 + dudz * dvdz * dz2 );
	  double t13 = ( dudc * dwdc * dc2 + dude * dwde * de2 + dudz * dwdz * dz2 );
	  double t21 = t12;
	  double t22 = ( dvdc * dvdc * dc2 + dvde * dvde * de2 + dvdz * dvdz * dz2 );
	  double t23 = ( dvdc * dwdc * dc2 + dvde * dwde * de2 + dvdz * dwdz * dz2 );
	  double t31 = t13;
	  double t32 = t23;
	  double t33 = ( dwdc * dwdc * dc2 + dwde * dwde * de2 + dwdz * dwdz * dz2 );
			
	  fp2[k][j][i].x -= ( t11 * eta0 + t12 * eta1 + t13 * eta2 ) / 12.;
	  fp2[k][j][i].y -= ( t21 * eta0 + t22 * eta1 + t23 * eta2 ) / 12.;
	  fp2[k][j][i].z -= ( t31 * eta0 + t32 * eta1 + t33 * eta2 ) / 12.;
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
		
	if( les || (rans && ti>0) ) {
	  //nu_t = pow( 0.5 * ( sqrt(lnu_t[k][j][i]) + sqrt(lnu_t[k+1][j][i]) ), 2.0) * Sabs;
	  nu_t = 0.5 * (lnu_t[k][j][i] + lnu_t[k+1][j][i]);
	  if ( (user->bctype[4]==1 && k==0) || (user->bctype[5]==1 && k==mz-2) ) nu_t=0;
		
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
	  double dc, de, dz;
	  Calculatedxdydz(ajc, csi[k][j][i], eta[k][j][i], zet[k][j][i], &dc, &de, &dz);
	  double dc2=dc*dc, de2=de*de, dz2=dz*dz;
			
	  double t11 = ( dudc * dudc * dc2 + dude * dude * de2 + dudz * dudz * dz2 );
	  double t12 = ( dudc * dvdc * dc2 + dude * dvde * de2 + dudz * dvdz * dz2 );
	  double t13 = ( dudc * dwdc * dc2 + dude * dwde * de2 + dudz * dwdz * dz2 );
	  double t21 = t12;
	  double t22 = ( dvdc * dvdc * dc2 + dvde * dvde * de2 + dvdz * dvdz * dz2 );
	  double t23 = ( dvdc * dwdc * dc2 + dvde * dwde * de2 + dvdz * dwdz * dz2 );
	  double t31 = t13;
	  double t32 = t23;
	  double t33 = ( dwdc * dwdc * dc2 + dwde * dwde * de2 + dwdz * dwdz * dz2 );
			
	  fp3[k][j][i].x -= ( t11 * zet0 + t12 * zet1 + t13 * zet2 ) / 12.;
	  fp3[k][j][i].y -= ( t21 * zet0 + t22 * zet1 + t23 * zet2 ) / 12.;
	  fp3[k][j][i].z -= ( t31 * zet0 + t32 * zet1 + t33 * zet2 ) / 12.;
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
  } else if (rans) {
  
    DMDAVecRestoreArray(da, user->lNu_t, &lnu_t);
  }
 

  VecDestroy(&Fp1);
  VecDestroy(&Fp2);
  VecDestroy(&Fp3);


  LOG_ALLOW(GLOBAL,LOG_DEBUG,"Viscous terms calculated .\n");

  PROFILE_FUNCTION_END;
  
  return(0);
}

#undef __FUNCT__
#define __FUNCT__ "ComputeRHS"
/**
 * @brief Computes the Right-Hand-Side (RHS) of the momentum equations.
 *
 * This function orchestrates the calculation of the spatial discretization of the
 * convective and diffusive terms. It calls specialized helper functions
 * (`Convection`, `Viscous`) to compute these terms and then combines them with
 * the pressure gradient to form the final RHS vector.
 *
 * This is a minimally-edited version of the legacy kernel. It operates on a
 * single UserCtx (one grid block) and retrieves all global parameters
 * (Re, rans, les, etc.) from the master SimCtx.
 *
 * @param user The UserCtx for the specific block being computed.
 * @param Rhs  The PETSc Vec where the calculated RHS will be stored.
 * @return PetscErrorCode 0 on success.
 */
PetscErrorCode ComputeRHS(UserCtx *user, Vec Rhs)
{
    PetscErrorCode ierr;
    SimCtx *simCtx = user->simCtx;
    DM da = user->da, fda = user->fda;
    DMDALocalInfo info = user->info;
    PetscInt i,j,k;
    PetscReal dpdc, dpde,dpdz;
    // --- Local Grid Indices and Parameters ---
    PetscInt xs = info.xs, xe = xs + info.xm, mx = info.mx;
    PetscInt ys = info.ys, ye = ys + info.ym, my = info.my;
    PetscInt zs = info.zs, ze = zs + info.zm, mz = info.mz;
    PetscInt lxs = (xs==0) ? xs+1 : xs;
    PetscInt lys = (ys==0) ? ys+1 : ys;
    PetscInt lzs = (zs==0) ? zs+1 : zs;
    PetscInt lxe = (xe==mx) ? xe-1 : xe;
    PetscInt lye = (ye==my) ? ye-1 : ye;
    PetscInt lze = (ze==mz) ? ze-1 : ze;
    PetscInt mz_end = (user->bctype[5] == 8) ? mz - 2 : mz - 3;

    // --- Array Pointers ---
    Cmpnts ***csi, ***eta, ***zet, ***icsi, ***ieta, ***izet, ***jcsi, ***jeta, ***jzet, ***kcsi, ***keta, ***kzet;
    PetscReal ***p, ***iaj, ***jaj, ***kaj, ***aj, ***nvert, ***nvert_o;
    Cmpnts ***rhs, ***rc, ***rct,***ucont;
    
    // --- Temporary Vectors ---
    Vec Conv, Visc, Rc, Rct;

    PetscFunctionBeginUser;
	PROFILE_FUNCTION_BEGIN;
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d, Block %d: Computing RHS (FormFunction1)...\n",
              simCtx->rank, user->_this);

    // --- Get all necessary array pointers ---
    ierr = DMDAVecGetArrayRead(fda, user->lCsi, &csi); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(fda, user->lEta, &eta); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(fda, user->lZet, &zet); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da,  user->lAj,  &aj);  CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(fda, user->lICsi, &icsi); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(fda, user->lIEta, &ieta); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(fda, user->lIZet, &izet); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(fda, user->lJCsi, &jcsi); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(fda, user->lJEta, &jeta); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(fda, user->lJZet, &jzet); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(fda, user->lKCsi, &kcsi); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(fda, user->lKEta, &keta); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(fda, user->lKZet, &kzet); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da, user->lIAj, &iaj); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da, user->lJAj, &jaj); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da, user->lKAj, &kaj); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da, user->lP, &p); CHKERRQ(ierr);
    ierr = DMDAVecGetArrayRead(da, user->lNvert, &nvert); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(fda, Rhs, &rhs); CHKERRQ(ierr);

    // --- Create temporary work vectors ---
    ierr = VecDuplicate(user->lUcont, &Rc); CHKERRQ(ierr);
    ierr = VecDuplicate(Rc, &Rct); CHKERRQ(ierr);
    ierr = VecDuplicate(Rct, &Conv); CHKERRQ(ierr);
    ierr = VecDuplicate(Rct, &Visc); CHKERRQ(ierr);

    // ========================================================================
    //   CORE LOGIC (UNCHANGED FROM LEGACY CODE)
    // ========================================================================

    // 1. Obtain Cartesian velocity from Contravariant velocity
    ierr = Contra2Cart(user); CHKERRQ(ierr);

    // 2. Compute Convective term
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Calculating convective terms...\n");
    if (simCtx->moveframe || simCtx->rotateframe) {
      // ierr = Convection_MV(user, user->lUcont, user->lUcat, Conv); CHKERRQ(ierr);
    } else {
        ierr = Convection(user, user->lUcont, user->lUcat, Conv); CHKERRQ(ierr);
    }

    // 3. Compute Viscous term
    if (simCtx->invicid) {
        ierr = VecSet(Visc, 0.0); CHKERRQ(ierr);
    } else {
        LOG_ALLOW(LOCAL, LOG_DEBUG, "  Calculating viscous terms...\n");
        ierr = Viscous(user, user->lUcont, user->lUcat, Visc); CHKERRQ(ierr);
    }

    // 4. Combine terms to get Cartesian RHS: Rc = Visc - Conv
    ierr = VecWAXPY(Rc, -1.0, Conv, Visc); CHKERRQ(ierr);

    // 5. Convert Cartesian RHS (Rc) to Contravariant RHS (Rct)
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Converting Cartesian RHS to Contravariant RHS...\n");
    ierr = DMDAVecGetArray(fda, Rct, &rct); CHKERRQ(ierr);
    ierr = DMDAVecGetArray(fda, Rc, &rc); CHKERRQ(ierr);
    
    for (k = lzs; k < lze; k++) {
        for (j = lys; j < lye; j++) {
            for (i = lxs; i < lxe; i++) {
                rct[k][j][i].x = aj[k][j][i] *
                    (0.5 * (csi[k][j][i].x + csi[k][j][i-1].x) * rc[k][j][i].x +
                     0.5 * (csi[k][j][i].y + csi[k][j][i-1].y) * rc[k][j][i].y +
                     0.5 * (csi[k][j][i].z + csi[k][j][i-1].z) * rc[k][j][i].z);
                rct[k][j][i].y = aj[k][j][i] *
                    (0.5 * (eta[k][j][i].x + eta[k][j-1][i].x) * rc[k][j][i].x +
                     0.5 * (eta[k][j][i].y + eta[k][j-1][i].y) * rc[k][j][i].y +
                     0.5 * (eta[k][j][i].z + eta[k][j-1][i].z) * rc[k][j][i].z);
                rct[k][j][i].z = aj[k][j][i] *
                    (0.5 * (zet[k][j][i].x + zet[k-1][j][i].x) * rc[k][j][i].x +
                     0.5 * (zet[k][j][i].y + zet[k-1][j][i].y) * rc[k][j][i].y +
                     0.5 * (zet[k][j][i].z + zet[k-1][j][i].z) * rc[k][j][i].z);
            }
        }
    }
    ierr = DMDAVecRestoreArray(fda, Rct, &rct); CHKERRQ(ierr);
    ierr = DMDAVecRestoreArray(fda, Rc, &rc); CHKERRQ(ierr);

    PetscBarrier(NULL);
 
    DMLocalToLocalBegin(fda, Rct, INSERT_VALUES, Rct);
    DMLocalToLocalEnd(fda, Rct, INSERT_VALUES, Rct);

    DMDAVecGetArray(fda, Rct, &rct);

    if (user->bctype[0]==7 || user->bctype[1]==7){
      if (xs==0){
	i=xs;
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    rct[k][j][i].x=rct[k][j][i-2].x;
	  }
	}
      }

      if (xe==mx){
	i=mx-1;
	for (k=lzs; k<lze; k++) {
	  for (j=lys; j<lye; j++) {
	    rct[k][j][i].x=rct[k][j][i+2].x;
	  }
	}
      }
    }

    if (user->bctype[2]==7 || user->bctype[3]==7){
      if (ys==0){
	j=ys;
	for (k=lzs; k<lze; k++) {
	  for (i=lxs; i<lxe; i++) {
	    rct[k][j][i].y=rct[k][j-2][i].y;
	  }
	}
      }

      if (ye==my){
	j=my-1;
	for (k=lzs; k<lze; k++) {
	  for (i=lxs; i<lxe; i++) {
	    rct[k][j][i].y=rct[k][j+2][i].y;
	  }
	}
      }
    }
    if (user->bctype[4]==7 || user->bctype[5]==7){
      if (zs==0){
	k=zs;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    rct[k][j][i].z=rct[k-2][j][i].z;
	  }
	}
      }
      if (ze==mz){
	k=mz-1;
	for (j=lys; j<lye; j++) {
	  for (i=lxs; i<lxe; i++) {
	    rct[k][j][i].z=rct[k+2][j][i].z;
	  }
	}
      }
    }
    
    // 6. Add Pressure Gradient Term and Finalize RHS
    // This involves calculating pressure derivatives (dpdc, dpde, dpdz) and using
    // them to adjust the contravariant RHS. The full stencil logic is preserved.
    LOG_ALLOW(LOCAL, LOG_DEBUG, "  Adding pressure gradient term to RHS...\n");
    
    ierr = DMDAVecRestoreArray(fda, Rct, &rct); CHKERRQ(ierr);
    
    ierr = DMLocalToLocalBegin(fda, Rct, INSERT_VALUES, Rct); CHKERRQ(ierr);
    ierr = DMLocalToLocalEnd(fda, Rct, INSERT_VALUES, Rct); CHKERRQ(ierr);
    
    ierr = DMDAVecGetArray(fda, Rct, &rct); CHKERRQ(ierr);

    for (k = lzs; k < lze; k++) {
        for (j = lys; j < lye; j++) {
            for (i = lxs; i < lxe; i++) {
                PetscReal dpdc, dpde, dpdz;
	dpdc = p[k][j][i+1] - p[k][j][i];

	if ((j==my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i]+nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
	    dpde = (p[k][j  ][i] + p[k][j  ][i+1] -
		    p[k][j-1][i] - p[k][j-1][i+1]) * 0.5;
	  }
	}
	else if ((j==my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i]+nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1) {
	    dpde = (p[k][j  ][i] + p[k][j  ][i+1] -
		    p[k][j-1][i] - p[k][j-1][i+1]) * 0.5;
	  }
	}
	else if ((j == 1 && user->bctype[2]!=7) || nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		    p[k][j  ][i] - p[k][j  ][i+1]) * 0.5;
	  }
	}
	else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		    p[k][j  ][i] - p[k][j  ][i+1]) * 0.5;
	  }
	}
	else {
	  dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		  p[k][j-1][i] - p[k][j-1][i+1]) * 0.25;
	}

	if ((k == mz-2 && user->bctype[4]!=7) || nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
	    dpdz = (p[k  ][j][i] + p[k  ][j][i+1] -
		    p[k-1][j][i] - p[k-1][j][i+1]) * 0.5;
	  }
	}
	else  if ((k == mz-2 || k==1) &&  user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1) {
	    dpdz = (p[k  ][j][i] + p[k  ][j][i+1] -
		    p[k-1][j][i] - p[k-1][j][i+1]) * 0.5;
	  }
	}
	else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		    p[k  ][j][i] - p[k  ][j][i+1]) * 0.5;
	  }
	}
	else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		    p[k  ][j][i] - p[k  ][j][i+1]) * 0.5;
	  }
	}
	else {
	  dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		  p[k-1][j][i] - p[k-1][j][i+1]) * 0.25;
	}
	    
	rhs[k][j][i].x =0.5 * (rct[k][j][i].x + rct[k][j][i+1].x);
	
	
	rhs[k][j][i].x -=
	  (dpdc * (icsi[k][j][i].x * icsi[k][j][i].x +
		   icsi[k][j][i].y * icsi[k][j][i].y +
		   icsi[k][j][i].z * icsi[k][j][i].z)+
	   dpde * (ieta[k][j][i].x * icsi[k][j][i].x +
		   ieta[k][j][i].y * icsi[k][j][i].y +
		   ieta[k][j][i].z * icsi[k][j][i].z)+
	   dpdz * (izet[k][j][i].x * icsi[k][j][i].x +
		   izet[k][j][i].y * icsi[k][j][i].y +
		   izet[k][j][i].z * icsi[k][j][i].z)) * iaj[k][j][i];

	if ((i == mx-2 && user->bctype[0]!=7) || nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
	    dpdc = (p[k][j][i  ] + p[k][j+1][i  ] -
		    p[k][j][i-1] - p[k][j+1][i-1]) * 0.5;
	  }
	}
	else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1) {
	    dpdc = (p[k][j][i  ] + p[k][j+1][i  ] -
		    p[k][j][i-1] - p[k][j+1][i-1]) * 0.5;
	  }
	}
	else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		    p[k][j][i  ] - p[k][j+1][i  ]) * 0.5;
	  }
	}
	else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		    p[k][j][i  ] - p[k][j+1][i  ]) * 0.5;
	  }
	}
	else {
	  dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		  p[k][j][i-1] - p[k][j+1][i-1]) * 0.25;
	}
	
	dpde = p[k][j+1][i] - p[k][j][i];
	
	if ((k == mz-2 && user->bctype[4]!=7)|| nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
	    dpdz = (p[k  ][j][i] + p[k  ][j+1][i] -
		    p[k-1][j][i] - p[k-1][j+1][i]) * 0.5;
	  }
	}
	else if ((k == mz-2 || k==1) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1) {
	    dpdz = (p[k  ][j][i] + p[k  ][j+1][i] -
		    p[k-1][j][i] - p[k-1][j+1][i]) * 0.5;
	  }
	}
	else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		    p[k  ][j][i] - p[k  ][j+1][i]) * 0.5;
	  }
	}
	else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		    p[k  ][j][i] - p[k  ][j+1][i]) * 0.5;
	  }
	}
	else {
	  dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		  p[k-1][j][i] - p[k-1][j+1][i]) * 0.25;
	}

	rhs[k][j][i].y =0.5 * (rct[k][j][i].y + rct[k][j+1][i].y);
	  
	  
	rhs[k][j][i].y -=
	  (dpdc * (jcsi[k][j][i].x * jeta[k][j][i].x +
		   jcsi[k][j][i].y * jeta[k][j][i].y +
		   jcsi[k][j][i].z * jeta[k][j][i].z) +
	   dpde * (jeta[k][j][i].x * jeta[k][j][i].x +
		   jeta[k][j][i].y * jeta[k][j][i].y +
		   jeta[k][j][i].z * jeta[k][j][i].z) +
	   dpdz * (jzet[k][j][i].x * jeta[k][j][i].x +
		   jzet[k][j][i].y * jeta[k][j][i].y +
		   jzet[k][j][i].z * jeta[k][j][i].z)) * jaj[k][j][i];

	if ((i == mx-2 && user->bctype[0]!=7)|| nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
	    dpdc = (p[k][j][i  ] + p[k+1][j][i  ] -
		    p[k][j][i-1] - p[k+1][j][i-1]) * 0.5;
	  }
	}
	else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1) {
	    dpdc = (p[k][j][i  ] + p[k+1][j][i  ] -
		    p[k][j][i-1] - p[k+1][j][i-1]) * 0.5;
	  }
	}
	else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		    p[k][j][i  ] - p[k+1][j][i  ]) * 0.5;
	  }
	}
	else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		    p[k][j][i  ] - p[k+1][j][i  ]) * 0.5;
	  }
	}
	else {
	  dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		  p[k][j][i-1] - p[k+1][j][i-1]) * 0.25;
	}

	if ((j == my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
	    dpde = (p[k][j  ][i] + p[k+1][j  ][i] -
		    p[k][j-1][i] - p[k+1][j-1][i]) * 0.5;
	  }
	}
	else  if ((j == my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1) {
	    dpde = (p[k][j  ][i] + p[k+1][j  ][i] -
		    p[k][j-1][i] - p[k+1][j-1][i]) * 0.5;
	  }
	}
	else if ((j == 1 && user->bctype[2]!=7)|| nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		    p[k][j  ][i] - p[k+1][j  ][i]) * 0.5;
	  }
	}
	else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		    p[k][j  ][i] - p[k+1][j  ][i]) * 0.5;
	  }
	}
	else {
	  dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		  p[k][j-1][i] - p[k+1][j-1][i]) * 0.25;
	}

	dpdz = (p[k+1][j][i] - p[k][j][i]);
	
	rhs[k][j][i].z =0.5 * (rct[k][j][i].z + rct[k+1][j][i].z);
	
	rhs[k][j][i].z -=
	  (dpdc * (kcsi[k][j][i].x * kzet[k][j][i].x +
		   kcsi[k][j][i].y * kzet[k][j][i].y +
		   kcsi[k][j][i].z * kzet[k][j][i].z) +
	   dpde * (keta[k][j][i].x * kzet[k][j][i].x +
		   keta[k][j][i].y * kzet[k][j][i].y +
		   keta[k][j][i].z * kzet[k][j][i].z) +
	   dpdz * (kzet[k][j][i].x * kzet[k][j][i].x +
		   kzet[k][j][i].y * kzet[k][j][i].y +
		   kzet[k][j][i].z * kzet[k][j][i].z)) * kaj[k][j][i];
	
      }
    }
  }
  
  
  //Mohsen March 2012//
  
  // rhs.x at boundaries for periodic bc at i direction//
  if (user->bctype[0]==7 && xs==0){
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {	  
	i=xs;

	dpdc = p[k][j][i+1] - p[k][j][i];
	
	if ((j==my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i]+nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
	    dpde = (p[k][j  ][i] + p[k][j  ][i+1] -
		    p[k][j-1][i] - p[k][j-1][i+1]) * 0.5;
	  }
	}
	else if ((j==my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i]+nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1) {
	    dpde = (p[k][j  ][i] + p[k][j  ][i+1] -
		    p[k][j-1][i] - p[k][j-1][i+1]) * 0.5;
	  }
	}
	else if ((j == 1 && user->bctype[2]!=7) || nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		    p[k][j  ][i] - p[k][j  ][i+1]) * 0.5;
	  }
	}
	else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		    p[k][j  ][i] - p[k][j  ][i+1]) * 0.5;
	  }
	}
	else {
	  dpde = (p[k][j+1][i] + p[k][j+1][i+1] -
		  p[k][j-1][i] - p[k][j-1][i+1]) * 0.25;
	}

	if ((k == mz-2 && user->bctype[4]!=7) || nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
	    dpdz = (p[k  ][j][i] + p[k  ][j][i+1] -
		    p[k-1][j][i] - p[k-1][j][i+1]) * 0.5;
	  }
	}
	else  if ((k == mz-2 || k==1) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1) {
	    dpdz = (p[k  ][j][i] + p[k  ][j][i+1] -
		    p[k-1][j][i] - p[k-1][j][i+1]) * 0.5;
	  }
	}
	else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		    p[k  ][j][i] - p[k  ][j][i+1]) * 0.5;
	  }
	}
	else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		    p[k  ][j][i] - p[k  ][j][i+1]) * 0.5;
	  }
	}
	else {
	  dpdz = (p[k+1][j][i] + p[k+1][j][i+1] -
		  p[k-1][j][i] - p[k-1][j][i+1]) * 0.25;
	}

	rhs[k][j][i].x =0.5 * (rct[k][j][i].x + rct[k][j][i+1].x);
	rhs[k][j][i].x -=
	  (dpdc * (icsi[k][j][i].x * icsi[k][j][i].x +
		   icsi[k][j][i].y * icsi[k][j][i].y +
		   icsi[k][j][i].z * icsi[k][j][i].z)+
	   dpde * (ieta[k][j][i].x * icsi[k][j][i].x +
		   ieta[k][j][i].y * icsi[k][j][i].y +
		   ieta[k][j][i].z * icsi[k][j][i].z)+
	   dpdz * (izet[k][j][i].x * icsi[k][j][i].x +
		   izet[k][j][i].y * icsi[k][j][i].y +
		   izet[k][j][i].z * icsi[k][j][i].z)) * iaj[k][j][i];
      }
    }
  }
  
// rhs.y at boundaries for periodic bc at j direction//
  if (user->bctype[2]==7 && ys==0){
    for (k=lzs; k<lze; k++) {
      for (i=lxs; i<lxe; i++) {	  

	j=ys;

	if ((i == mx-2 && user->bctype[0]!=7) || nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
	    dpdc = (p[k][j][i  ] + p[k][j+1][i  ] -
		    p[k][j][i-1] - p[k][j+1][i-1]) * 0.5;
	  }
	}
	else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1) {
	    dpdc = (p[k][j][i  ] + p[k][j+1][i  ] -
		    p[k][j][i-1] - p[k][j+1][i-1]) * 0.5;
	  }
	}
	else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		    p[k][j][i  ] - p[k][j+1][i  ]) * 0.5;
	  }
	}
	else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		    p[k][j][i  ] - p[k][j+1][i  ]) * 0.5;
	  }
	}
	else {
	  dpdc = (p[k][j][i+1] + p[k][j+1][i+1] -
		  p[k][j][i-1] - p[k][j+1][i-1]) * 0.25;
	}

	dpde = p[k][j+1][i] - p[k][j][i];

	if ((k == mz-2 && user->bctype[4]!=7)|| nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
	    dpdz = (p[k  ][j][i] + p[k  ][j+1][i] -
		    p[k-1][j][i] - p[k-1][j+1][i]) * 0.5;
	  }
	}
	else if ((k == mz-2 || k==1 ) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1) {
	    dpdz = (p[k  ][j][i] + p[k  ][j+1][i] -
		    p[k-1][j][i] - p[k-1][j+1][i]) * 0.5;
	  }
	}
	else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		    p[k  ][j][i] - p[k  ][j+1][i]) * 0.5;
	  }
	}
	else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
	  if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		    p[k  ][j][i] - p[k  ][j+1][i]) * 0.5;
	  }
	}
	else {
	  dpdz = (p[k+1][j][i] + p[k+1][j+1][i] -
		  p[k-1][j][i] - p[k-1][j+1][i]) * 0.25;
	}

	rhs[k][j][i].y =0.5 * (rct[k][j][i].y + rct[k][j+1][i].y);
	 
       	rhs[k][j][i].y -=
	  (dpdc * (jcsi[k][j][i].x * jeta[k][j][i].x +
		   jcsi[k][j][i].y * jeta[k][j][i].y +
		   jcsi[k][j][i].z * jeta[k][j][i].z)+
	   dpde * (jeta[k][j][i].x * jeta[k][j][i].x +
		   jeta[k][j][i].y * jeta[k][j][i].y +
		   jeta[k][j][i].z * jeta[k][j][i].z)+
	   dpdz * (jzet[k][j][i].x * jeta[k][j][i].x +
		   jzet[k][j][i].y * jeta[k][j][i].y +
		   jzet[k][j][i].z * jeta[k][j][i].z)) * jaj[k][j][i];

      }
    }
  }
  
  // rhs.z at boundaries for periodic bc at k direction//
  if (user->bctype[4]==7&& zs==0){
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

	k=zs; 

	if ((i == mx-2 && user->bctype[0]!=7)|| nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
	    dpdc = (p[k][j][i  ] + p[k+1][j][i  ] -
		    p[k][j][i-1] - p[k+1][j][i-1]) * 0.5;
	  }
	}
	else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
	  if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1) {
	    dpdc = (p[k][j][i  ] + p[k+1][j][i  ] -
		    p[k][j][i-1] - p[k+1][j][i-1]) * 0.5;
	  }
	}
	else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		    p[k][j][i  ] - p[k+1][j][i  ]) * 0.5;
	  }
	}
	else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
	  if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
	    dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		    p[k][j][i  ] - p[k+1][j][i  ]) * 0.5;
	  }
	}
	else {
	  dpdc = (p[k][j][i+1] + p[k+1][j][i+1] -
		  p[k][j][i-1] - p[k+1][j][i-1]) * 0.25;
	}

	if ((j == my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
	    dpde = (p[k][j  ][i] + p[k+1][j  ][i] -
		    p[k][j-1][i] - p[k+1][j-1][i]) * 0.5;
	  }
	}
	else  if ((j == my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1) {
	    dpde = (p[k][j  ][i] + p[k+1][j  ][i] -
		    p[k][j-1][i] - p[k+1][j-1][i]) * 0.5;
	  }
	}
	else if ((j == 1 && user->bctype[2]!=7)|| nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		    p[k][j  ][i] - p[k+1][j  ][i]) * 0.5;
	  }
	}
	else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
	  if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
	    dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		    p[k][j  ][i] - p[k+1][j  ][i]) * 0.5;
	  }
	}
	else {
	  dpde = (p[k][j+1][i] + p[k+1][j+1][i] -
		  p[k][j-1][i] - p[k+1][j-1][i]) * 0.25;
	}

	dpdz = (p[k+1][j][i] - p[k][j][i]);
	
	rhs[k][j][i].z =0.5 * (rct[k][j][i].z + rct[k+1][j][i].z);
	
	rhs[k][j][i].z -=
	  (dpdc * (kcsi[k][j][i].x * kzet[k][j][i].x +
		   kcsi[k][j][i].y * kzet[k][j][i].y +
		   kcsi[k][j][i].z * kzet[k][j][i].z)+
	   dpde * (keta[k][j][i].x * kzet[k][j][i].x +
		   keta[k][j][i].y * kzet[k][j][i].y +
		   keta[k][j][i].z * kzet[k][j][i].z)+
	   dpdz * (kzet[k][j][i].x * kzet[k][j][i].x +
		   kzet[k][j][i].y * kzet[k][j][i].y +
		   kzet[k][j][i].z * kzet[k][j][i].z)) * kaj[k][j][i];
	
      }
    }
  }
    
    ierr = DMDAVecRestoreArray(fda, Rct, &rct); CHKERRQ(ierr);

    LOG_ALLOW(GLOBAL,LOG_DEBUG,"Pressure Gradient added to RHS .\n");
    PetscInt TwoD = simCtx->TwoD;

    LOG_ALLOW(GLOBAL,LOG_DEBUG,"Final cleanup and edge-cases initiated .\n");
    
    // 7. Final clean-up for immersed boundaries and 2D cases
      for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (TwoD==1)
	  rhs[k][j][i].x =0.;
	else if (TwoD==2)
	  rhs[k][j][i].y =0.;
	else if (TwoD==3)
	  rhs[k][j][i].z =0.;
	
	if (nvert[k][j][i]>0.1) {
	  rhs[k][j][i].x = 0;
	  rhs[k][j][i].y = 0;
	  rhs[k][j][i].z = 0;
	}
	if (nvert[k][j][i+1]>0.1) {
	  rhs[k][j][i].x=0;
	}
	if (nvert[k][j+1][i]>0.1) {
	  rhs[k][j][i].y=0;
	}
	if (nvert[k+1][j][i]>0.1) {
	  rhs[k][j][i].z=0;
	}
      }
    }
  }
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"Final cleanup and edge-cases complete .\n");

      // ========================================================================

    // --- Restore all PETSc array pointers ---
    // DMDAVecRestoreArray(fda, user->lUcont, &ucont);   
    DMDAVecRestoreArray(fda, Rhs, &rhs);
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"Rhs restored successfully! .\n");
    
    DMDAVecRestoreArray(fda, user->lCsi, &csi);
    DMDAVecRestoreArray(fda, user->lEta, &eta);
    DMDAVecRestoreArray(fda, user->lZet, &zet);
    DMDAVecRestoreArray(da,  user->lAj,  &aj);
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"Face metrics restored successfully! .\n");
    
    DMDAVecRestoreArray(fda, user->lICsi, &icsi);
    DMDAVecRestoreArray(fda, user->lIEta, &ieta);
    DMDAVecRestoreArray(fda, user->lIZet, &izet);
    DMDAVecRestoreArray(da, user->lIAj, &iaj);
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"I Face metrics restored successfully! .\n");
    
    DMDAVecRestoreArray(fda, user->lJCsi, &jcsi);
    DMDAVecRestoreArray(fda, user->lJEta, &jeta);
    DMDAVecRestoreArray(fda, user->lJZet, &jzet);
    DMDAVecRestoreArray(da, user->lJAj, &jaj);  
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"J Face metrics restored successfully! .\n");
  
    DMDAVecRestoreArray(fda, user->lKCsi, &kcsi);
    DMDAVecRestoreArray(fda, user->lKEta, &keta);
    DMDAVecRestoreArray(fda, user->lKZet, &kzet);
    DMDAVecRestoreArray(da, user->lKAj, &kaj);
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"K Face metrics restored successfully! .\n");
  
    DMDAVecRestoreArray(da, user->lP, &p);
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"Pressure restored successfully! .\n");

    DMDAVecRestoreArray(da, user->lNvert, &nvert);
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"Nvert restored successfully! .\n");

    LOG_ALLOW(GLOBAL,LOG_DEBUG,"Cell Centered scalars restored successfully! .\n");
    
    // --- Destroy temporary work vectors ---
    ierr = VecDestroy(&Conv); CHKERRQ(ierr);
    ierr = VecDestroy(&Visc); CHKERRQ(ierr);
    ierr = VecDestroy(&Rc); CHKERRQ(ierr);
    ierr = VecDestroy(&Rct); CHKERRQ(ierr);
    LOG_ALLOW(GLOBAL,LOG_DEBUG,"Temporary work vectors destroyed successfully! .\n");
    
    LOG_ALLOW(LOCAL, LOG_DEBUG, "Rank %d, Block %d: RHS computation complete.\n",
              simCtx->rank, user->_this);


	PROFILE_FUNCTION_END;		  
    PetscFunctionReturn(0);
}
