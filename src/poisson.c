// In src/poisson.c (or wherever PoissonSolver_MG is)
#include "poisson.h" // The new header for this file
#include "logging.h"

static void Calculate_Covariant_metrics(double g[3][3], double G[3][3])
{
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
};

static void Calculate_normal(Cmpnts csi, Cmpnts eta, Cmpnts zet, double ni[3], double nj[3], double nk[3])
{
	double g[3][3];
	double G[3][3];
	
	g[0][0]=csi.x, g[0][1]=csi.y, g[0][2]=csi.z;
	g[1][0]=eta.x, g[1][1]=eta.y, g[1][2]=eta.z;
	g[2][0]=zet.x, g[2][1]=zet.y, g[2][2]=zet.z;
	
	Calculate_Covariant_metrics(g, G);
	double xcsi=G[0][0], ycsi=G[1][0], zcsi=G[2][0];
	double xeta=G[0][1], yeta=G[1][1], zeta=G[2][1];
	double xzet=G[0][2], yzet=G[1][2], zzet=G[2][2];
	      
	double nx_i = xcsi, ny_i = ycsi, nz_i = zcsi;
	double nx_j = xeta, ny_j = yeta, nz_j = zeta;
	double nx_k = xzet, ny_k = yzet, nz_k = zzet;
	      
	double sum_i=sqrt(nx_i*nx_i+ny_i*ny_i+nz_i*nz_i);
	double sum_j=sqrt(nx_j*nx_j+ny_j*ny_j+nz_j*nz_j);
	double sum_k=sqrt(nx_k*nx_k+ny_k*ny_k+nz_k*nz_k);

//	*Ai = sqrt( g[0][0]*g[0][0] + g[0][1]*g[0][1] + g[0][2]*g[0][2] );	// area
//	*Aj = sqrt( g[1][0]*g[1][0] + g[1][1]*g[1][1] + g[1][2]*g[1][2] );
//	*Ak =sqrt( g[2][0]*g[2][0] + g[2][1]*g[2][1] + g[2][2]*g[2][2] );
		
	nx_i /= sum_i, ny_i /= sum_i, nz_i /= sum_i;
	nx_j /= sum_j, ny_j /= sum_j, nz_j /= sum_j;
	nx_k /= sum_k, ny_k /= sum_k, nz_k /= sum_k;

	ni[0] = nx_i, ni[1] = ny_i, ni[2] = nz_i;
	nj[0] = nx_j, nj[1] = ny_j, nj[2] = nz_j;
	nk[0] = nx_k, nk[1] = ny_k, nk[2] = nz_k;
}

static PetscErrorCode GhostNodeVelocity(UserCtx *user)
{

  // --- CONTEXT ACQUISITION BLOCK ---
  // Get the master simulation context from the UserCtx.
  SimCtx *simCtx = user->simCtx;

  // Create local variables to mirror the legacy globals for minimal code changes.
  PetscInt wallfunction = simCtx->wallfunction;
  const PetscInt les = simCtx->les;
  const PetscInt rans = simCtx->rans;
  const PetscInt ti = simCtx->step;
  const PetscReal U_bc = simCtx->U_bc;
  // --- END CONTEXT ACQUISITION BLOCK ---
  
  DM            da = user->da, fda = user->fda;
  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt	i, j, k;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  Cmpnts        ***ubcs, ***ucat,***lucat, ***csi, ***eta, ***zet,***cent;


  PetscReal Un, nx,ny,nz,A;

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

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  DMDAVecGetArray(fda, user->lCsi,  &csi);
  DMDAVecGetArray(fda, user->lEta,  &eta);
  DMDAVecGetArray(fda, user->lZet,  &zet);
 
  
/* ==================================================================================             */
/*   WALL FUNCTION */
/* ==================================================================================             */
  if (wallfunction) {
    PetscReal ***ustar, ***aj, ***nvert;
  DMDAVecGetArray(fda, user->Ucat, &ucat);
  DMDAVecGetArray(da, user->lUstar, &ustar);
  DMDAVecGetArray(da, user->lAj,  &aj);
  DMDAVecGetArray(da, user->lNvert,  &nvert);

  // wall function for boundary
  //Dec 2015
  for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
      for (i=lxs; i<lxe; i++) {
	if( nvert[k][j][i]<1.1 && ( (user->bctype[0]==-1 && i==1) || (user->bctype[1]==-1 &&  i==mx-2) ) ) {
	  double area = sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
	  double sb, sc; 
	  double ni[3], nj[3], nk[3];
	  Cmpnts Uc, Ua;
	  
	  Ua.x = Ua.y = Ua.z = 0;
	  sb = 0.5/aj[k][j][i]/area;
	  
	  if(i==1) {
	    sc = 2*sb + 0.5/aj[k][j][i+1]/area;
	    Uc = ucat[k][j][i+1];
	  }
	  else {
	    sc = 2*sb + 0.5/aj[k][j][i-1]/area;
	    Uc = ucat[k][j][i-1];
	  }
	  
	  Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
	  if(i==mx-2) ni[0]*=-1, ni[1]*=-1, ni[2]*=-1;
	  
	  //if(i==1) printf("%f %f, %f %f %f, %e %e %e, %e\n", sc, sb, Ua.x, Ua.y, Ua.z, Uc.x, Uc.y, Uc.z, ucat[k][j][i+2].z);
	  wall_function (user, sc, sb, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], ni[0], ni[1], ni[2]);
	  nvert[k][j][i]=1;	/* set nvert to 1 to exclude from rhs */
	  //	  ucat[k][j][i] = lucat[k][j][i];	  
	}
	
	if( nvert[k][j][i]<1.1 && ( (user->bctype[2]==-1 && j==1) || (user->bctype[3]==-1 &&  j==my-2) ) ) {
	  double area = sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
	  double sb, sc; 
	  double ni[3], nj[3], nk[3];
	  Cmpnts Uc, Ua;
	  
	  Ua.x = Ua.y = Ua.z = 0;
	  sb = 0.5/aj[k][j][i]/area;
	  
	  if(j==1) {
	    sc = 2*sb + 0.5/aj[k][j+1][i]/area;
	    Uc = ucat[k][j+1][i];
	  }
	  else {
	    sc = 2*sb + 0.5/aj[k][j-1][i]/area;
	    Uc = ucat[k][j-1][i];
	  }
	  
	  Calculate_normal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
	  //if(j==my-2) nj[0]*=-1, nj[1]*=-1, nj[2]*=-1;
	  
	  wall_function (user, sc, sb, Ua, Uc, &ucat[k][j][i], &ustar[k][j][i], nj[0], nj[1], nj[2]);
	  nvert[k][j][i]=1;	/* set nvert to 1 to exclude from rhs */
	  //	  ucat[k][j][i] = lucat[k][j][i];
	}
	
      }
  DMDAVecRestoreArray(da, user->lAj,  &aj);
  DMDAVecRestoreArray(da, user->lNvert,  &nvert);
  DMDAVecRestoreArray(da, user->lUstar, &ustar);   
  DMDAVecRestoreArray(fda, user->Ucat, &ucat);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
 /*  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */

  } // if wallfunction

/* ==================================================================================             */
/*   Driven Cavity  */
/* ==================================================================================             */

  DMDAVecGetArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecGetArray(fda, user->Ucat,  &ucat);

  /* Designed for driven cavity problem (top(k=kmax) wall moving)
   u_x = 1 at k==kmax */
  if (user->bctype[5]==2) {
    if (ze==mz) {
      k = ze-1;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;// - ucat[k-1][j][i].x;
	  ubcs[k][j][i].y = 1.;//- ucat[k-1][j][i].y;
	  ubcs[k][j][i].z = 0.;//- ucat[k-1][j][i].z;
	}
      }
    }
  }
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
/* designed to test periodic boundary condition for c grid j=0 rotates */
 
  if (user->bctype[2]==11) {
    DMDAVecGetArray(fda, user->Cent, &cent); 
   if (ys==0){
      j=0;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x =cent[k][j+1][i].y/sqrt(cent[k][j+1][i].x*cent[k][j+1][i].x+cent[k][j+1][i].y*cent[k][j+1][i].y);;
	  // ubcs[k][j][i].x=0.0;
	  ubcs[k][j][i].y =-cent[k][j+1][i].x/sqrt(cent[k][j+1][i].x*cent[k][j+1][i].x+cent[k][j+1][i].y*cent[k][j+1][i].y);
	  // ubcs[k][j][i].y = -cent[k][j+1][i].z/sqrt(cent[k][j+1][i].z*cent[k][j+1][i].z+cent[k][j+1][i].y*cent[k][j+1][i].y);
	  ubcs[k][j][i].z = 0.0;
	  // ubcs[k][j][i].z =cent[k][j+1][i].y/sqrt(cent[k][j+1][i].z*cent[k][j+1][i].z+cent[k][j+1][i].y*cent[k][j+1][i].y);
	}
      }
    }
   DMDAVecRestoreArray(fda, user->Cent,  &cent);
  }

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
	Un=ucat[k][j][i].x*nx+ucat[k][j][i].y*ny+ucat[k][j][i].z*nz;
	ubcs[k][j][i].x = ucat[k][j][i].x-Un*nx;
	ubcs[k][j][i].y = ucat[k][j][i].y-Un*ny;
	ubcs[k][j][i].z = ucat[k][j][i].z-Un*nz;
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
	}
      }
    }
  }
/*  ==================================================================================== */
  /*     Rheology */
  /*  ==================================================================================== */
 
  // PetscPrintf(PETSC_COMM_WORLD, "moving plate velocity for rheology setup is %le \n",U_bc);

  if (user->bctype[2]==13){
    if (ys==0){
      j=0;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  ubcs[k][j][i].x = 0.;
	  //ubcs[k][j][i].x = -U_bc;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z =-U_bc;
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
	  //ubcs[k][j][i].x = U_bc;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = U_bc;
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
	  ubcs[k][j][i].x =-U_bc;
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
	  ubcs[k][j][i].x = U_bc;
	  ubcs[k][j][i].y = 0.;
	  ubcs[k][j][i].z = 0.;
	}
      }
    }
  }
/* ==================================================================================             */
  // boundary conditions on ghost nodes
/* ==================================================================================             */
//Mohsen Aug 2012
  if (xs==0 && user->bctype[0]!=7) {
    i = xs;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i+1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i+1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i+1].z;
      }
    }
  }

  if (xe==mx && user->bctype[0]!=7) {
    i = xe-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j][i-1].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j][i-1].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j][i-1].z;
      }
    }
  }

  if (ys==0 && user->bctype[2]!=7) {
    j = ys;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j+1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j+1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j+1][i].z;
      }
    }
  }

  if (ye==my && user->bctype[2]!=7) {
    j = ye-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k][j-1][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k][j-1][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k][j-1][i].z;
      }
    }
  }

  if (zs==0 && user->bctype[4]!=7) {
    k = zs;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k+1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k+1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k+1][j][i].z;
      }
    }
  }

  if (ze==mz && user->bctype[4]!=7) {
    k = ze-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 2 * ubcs[k][j][i].x - ucat[k-1][j][i].x;
	ucat[k][j][i].y = 2 * ubcs[k][j][i].y - ucat[k-1][j][i].y;
	ucat[k][j][i].z = 2 * ubcs[k][j][i].z - ucat[k-1][j][i].z;
      }
    }
  }
  DMDAVecRestoreArray(fda, user->Bcs.Ubcs, &ubcs);
  DMDAVecRestoreArray(fda, user->Ucat,&ucat);
/* ==================================================================================             */
/*   Periodic BC Mohsen Aug 2012
/* ==================================================================================             */

  if (user->bctype[0]==7 || user->bctype[2]==7 ||  user->bctype[4]==7 ){

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
 
    DMDAVecGetArray(fda, user->lUcat, &lucat);
    DMDAVecGetArray(fda, user->Ucat, &ucat);
   
    if (user->bctype[0]==7 || user->bctype[1]==7){
      if (xs==0){
	i=xs;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    if(k>0 && k<user->KM && j>0 && j<user->JM){
	      ucat[k][j][i]=lucat[k][j][i-2];
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
	    if(k>0 && k<user->KM && i>0 && i<user->IM){
	      ucat[k][j][i]=lucat[k][j-2][i];
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
	    if(j>0 && j<user->JM && i>0 && i<user->IM){
	      ucat[k][j][i]=lucat[k-2][j][i];
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
	    if(k>0 && k<user->KM && j>0 && j<user->JM){
	      ucat[k][j][i]=lucat[k][j][i+2];
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
	    if(k>0 && k<user->KM && i>0 && i<user->IM){
	    ucat[k][j][i]=lucat[k][j+2][i];
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
	    if(j>0 && j<user->JM && i>0 && i<user->IM){
	      ucat[k][j][i].x=lucat[k+2][j][i].x;
	    }
	  }
	}
      }
    }
    
    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
 /*  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
  }

  // velocity on the corner point
  DMDAVecGetArray(fda, user->Ucat,  &ucat);

  if (zs==0 ) {
    k=0;
    if (xs==0) {
      i=0;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j][i+1].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j][i+1].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j][i+1].z);
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j][i-1].z);
      }
    }

    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j+1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j+1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j+1][i].z);
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k+1][j][i].x+ucat[k][j-1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k+1][j][i].y+ucat[k][j-1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k+1][j][i].z+ucat[k][j-1][i].z);
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
      }
    }
    if (xe == mx) {
      i=mx-1;
      for (j=ys; j<ye; j++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j][i-1].z);
      }
    }
    if (ys==0) {
      j=0;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j+1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j+1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j+1][i].z);
      }
    }

    if (ye==my) {
      j=my-1;
      for (i=xs; i<xe; i++) {
	ucat[k][j][i].x = 0.5*(ucat[k-1][j][i].x+ucat[k][j-1][i].x);
	ucat[k][j][i].y = 0.5*(ucat[k-1][j][i].y+ucat[k][j-1][i].y);
	ucat[k][j][i].z = 0.5*(ucat[k-1][j][i].z+ucat[k][j-1][i].z);
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
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j+1][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j+1][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j+1][i].z+ucat[k][j][i-1].z);
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
      }
    }

    if (xe==mx) {
      i=mx-1;
      for (k=zs; k<ze; k++) {
	ucat[k][j][i].x = 0.5*(ucat[k][j-1][i].x+ucat[k][j][i-1].x);
	ucat[k][j][i].y = 0.5*(ucat[k][j-1][i].y+ucat[k][j][i-1].y);
	ucat[k][j][i].z = 0.5*(ucat[k][j-1][i].z+ucat[k][j][i-1].z);
      }
    }
  }

  DMDAVecRestoreArray(fda, user->Ucat,  &ucat);

  DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  if (user->bctype[0]==7 || user->bctype[2]==7 ||  user->bctype[4]==7 ){
    //i-direction
    DMDAVecGetArray(fda, user->lUcat, &lucat);
    DMDAVecGetArray(fda, user->Ucat, &ucat);
 
    if (user->bctype[0]==7){
      if (xs==0){
	i=xs;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    ucat[k][j][i]=lucat[k][j][i-2];
	  }
	}
      }
    }
    if (user->bctype[1]==7){
      if (xe==mx){
	i=xe-1;
	for (k=zs; k<ze; k++) {
	  for (j=ys; j<ye; j++) {
	    ucat[k][j][i]=lucat[k][j][i+2];
	  }
	}
      }
    }
    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);

 /*  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  //j-direction
  
    DMDAVecGetArray(fda, user->lUcat, &lucat);
    DMDAVecGetArray(fda, user->Ucat, &ucat);
 
    if (user->bctype[2]==7){
      if (ys==0){
	j=ys;
	for (k=zs; k<ze; k++) {
	  for (i=xs; i<xe; i++) {
	    ucat[k][j][i]=lucat[k][j-2][i];
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
	  }
	}
      }
    }

    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);
  
 /*  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);

  //k-direction
  
    DMDAVecGetArray(fda, user->lUcat, &lucat);
    DMDAVecGetArray(fda, user->Ucat, &ucat);
 
    if (user->bctype[4]==7){
      if (zs==0){
	k=zs;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	  ucat[k][j][i]=lucat[k-2][j][i];
	  }
	}
      }
    }
    if (user->bctype[5]==7){
      if (ze==mz){
	k=mz-1;
	for (j=ys; j<ye; j++) {
	  for (i=xs; i<xe; i++) {
	    ucat[k][j][i].x=lucat[k+2][j][i].x;
	  }
	}
      }
    }
   
    DMDAVecRestoreArray(fda, user->lUcat, &lucat);
    DMDAVecRestoreArray(fda, user->Ucat, &ucat);

 /*  DMLocalToGlobalBegin(fda, user->lUcat, INSERT_VALUES, user->Ucat); */
/*   DMLocalToGlobalEnd(fda, user->lUcat, INSERT_VALUES, user->Ucat); */

    DMGlobalToLocalBegin(fda, user->Ucat, INSERT_VALUES, user->lUcat);
    DMGlobalToLocalEnd(fda, user->Ucat, INSERT_VALUES, user->lUcat);
  }
 
  DMDAVecRestoreArray(fda, user->lCsi,  &csi);
  DMDAVecRestoreArray(fda, user->lEta,  &eta);
  DMDAVecRestoreArray(fda, user->lZet,  &zet);

  return(0);
  
}

#define GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user) \
  if ((user->isc)) { \
    ic = i; \
    ia = 0; \
  } \
  else { \
    ic = (i+1) / 2; \
    ia = (i - 2 * (ic)) == 0 ? 1 : -1; \
    if (i==1 || i==mx-2) ia = 0; \
  }\
  if ((user->jsc)) { \
    jc = j; \
    ja = 0; \
  } \
  else { \
    jc = (j+1) / 2; \
    ja = (j - 2 * (jc)) == 0 ? 1 : -1; \
    if (j==1 || j==my-2) ja = 0; \
  } \
  if ((user->ksc)) { \
    kc = k; \
    ka = 0; \
  } \
  else { \
    kc = (k+1) / 2; \
    ka = (k - 2 * (kc)) == 0 ? 1 : -1; \
    if (k==1 || k==mz-2) ka = 0; \
  } \
  if (ka==-1 && nvert_c[kc-1][jc][ic] > 0.1) ka = 0; \
  else if (ka==1 && nvert_c[kc+1][jc][ic] > 0.1) ka = 0; \
  if (ja==-1 && nvert_c[kc][jc-1][ic] > 0.1) ja = 0; \
  else if (ja==1 && nvert_c[kc][jc+1][ic] > 0.1) ja = 0; \
  if (ia==-1 && nvert_c[kc][jc][ic-1] > 0.1) ia = 0; \
  else if (ia==1 && nvert_c[kc][jc][ic+1] > 0.1) ia = 0;

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

static PetscErrorCode GridRestriction(PetscInt i, PetscInt j, PetscInt k,
			       PetscInt *ih, PetscInt *jh, PetscInt *kh,
			       UserCtx *user)
{
  if ((user->isc)) {
    *ih = i;
  }
  else {
    *ih = 2 * i;
  }

  if ((user->jsc)) {
    *jh = j;
  }
  else {
    *jh = 2 * j;
  }

  if ((user->ksc)) {
    *kh = k;
  }
  else {
    *kh = 2 * k;
  }

  return 0;
}


//Cell Center
#define CP  0 
// Face Centers
#define EP  1
#define WP  2
#define NP  3
#define SP  4
#define TP  5
#define BP  6

// Edge Centers
#define NE  7
#define SE  8
#define NW  9
#define SW  10
#define TN  11
#define BN  12
#define TS  13
#define BS  14
#define TE  15
#define BE  16
#define TW  17
#define BW  18

/**
 * @brief Performs the projection step to enforce an incompressible velocity field.
 *
 * @details
 * This function executes the final "correction" stage of a fractional-step (projection)
 * method for solving the incompressible Navier-Stokes equations. After solving the
 * Pressure Poisson Equation to get a pressure correction `Phi`, this function uses
 * `Phi` to correct the intermediate velocity field, making it divergence-free.
 *
 * The main steps are:
 * 1.  **Calculate Pressure Gradient**: At each velocity point (on the cell faces), it
 *     computes the gradient of the pressure correction field (`∇Φ`). This is done using
 *     finite differences on a generalized curvilinear grid.
 *
 * 2.  **Correct Velocity Field**: It updates the contravariant velocity components `Ucont`
 *     by subtracting the scaled pressure gradient. The update rule is:
 *     `U_new = U_intermediate - Δt * G * ∇Φ`, where `G` is the metric tensor `g_ij`
 *     that transforms the gradient into contravariant components.
 *
 * 3.  **Boundary-Aware Gradients**: The gradient calculation is highly sensitive to
 *     boundaries. The code uses complex, dynamic stencils that switch from centered
 *     differences to one-sided differences if a standard stencil would use a point
 *     inside an immersed solid (`nvert > 0.1`). This preserves accuracy at the
 *     fluid-solid interface.
 *
 * 4.  **Periodic Boundary Updates**: Includes dedicated loops to correctly calculate
 *     the pressure gradient at the domain edges for periodic boundary conditions.
 *
 * 5.  **Finalization**: After correcting `Ucont`, it calls helper functions to convert
 *     the velocity back to Cartesian coordinates (`Contra2Cart`) and update all
 *     ghost cell values.
 *
 * @param[in, out] user Pointer to the UserCtx struct, containing simulation context,
 *                      grid metrics, and the velocity/pressure fields.
 *
 * @return PetscErrorCode 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode Projection(UserCtx *user)
{
  PetscFunctionBeginUser;
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Entering Projection step to correct velocity field.\n");

  //================================================================================
  // Section 1: Initialization and Data Acquisition
  //================================================================================

  // --- Get simulation and grid context ---
  SimCtx *simCtx = user->simCtx;
  DM da = user->da, fda = user->fda;
  DMDALocalInfo info = user->info;

  // --- Grid dimensions ---
  PetscInt mx = info.mx, my = info.my, mz = info.mz;
  PetscInt xs = info.xs, xe = info.xs + info.xm;
  PetscInt ys = info.ys, ye = info.ys + info.ym;
  PetscInt zs = info.zs, ze = info.zs + info.zm;

  // --- Loop bounds (excluding outer ghost layers) ---
  PetscInt lxs = (xs == 0) ? xs + 1 : xs;
  PetscInt lxe = (xe == mx) ? xe - 1 : xe;
  PetscInt lys = (ys == 0) ? ys + 1 : ys;
  PetscInt lye = (ye == my) ? ye - 1 : ye;
  PetscInt lzs = (zs == 0) ? zs + 1 : zs;
  PetscInt lze = (ze == mz) ? ze - 1 : ze;

  // --- Get direct pointer access to grid metric and field data ---
  Cmpnts ***icsi, ***ieta, ***izet, ***jcsi, ***jeta, ***jzet, ***kcsi, ***keta, ***kzet;
  PetscReal ***iaj, ***jaj, ***kaj, ***p, ***nvert;
  Cmpnts ***ucont;
  DMDAVecGetArray(fda, user->lICsi, &icsi); DMDAVecGetArray(fda, user->lIEta, &ieta); DMDAVecGetArray(fda, user->lIZet, &izet);
  DMDAVecGetArray(fda, user->lJCsi, &jcsi); DMDAVecGetArray(fda, user->lJEta, &jeta); DMDAVecGetArray(fda, user->lJZet, &jzet);
  DMDAVecGetArray(fda, user->lKCsi, &kcsi); DMDAVecGetArray(fda, user->lKEta, &keta); DMDAVecGetArray(fda, user->lKZet, &kzet);
  DMDAVecGetArray(da, user->lIAj, &iaj); DMDAVecGetArray(da, user->lJAj, &jaj); DMDAVecGetArray(da, user->lKAj, &kaj);
  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(da, user->lPhi, &p); // Note: using lPhi, which is the pressure correction
  //DMDAVecGetArray(da,user->lP,&p);
  DMDAVecGetArray(fda, user->Ucont, &ucont);

  // --- Constants for clarity ---
  const PetscReal IBM_FLUID_THRESHOLD = 0.1;
  const PetscReal scale = simCtx->dt * simCtx->st / COEF_TIME_ACCURACY;

  LOG_ALLOW(GLOBAL,LOG_DEBUG," Starting velocity correction: Scale = %le .\n",scale);

  //================================================================================
  // Section 2: Correct Velocity Components
  //================================================================================
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Calculating pressure gradients and correcting velocity components.\n");
  
  // --- Main loop over interior domain points ---
  for (PetscInt k = lzs; k < lze; k++) {
    for (PetscInt j = lys; j < lye; j++) {
      for (PetscInt i = lxs; i < lxe; i++) {

        // --- Correct U_contravariant (x-component of velocity) ---
        PetscInt i_end = (user->bctype[0] == 7) ? mx - 1 : mx - 2;
        if (i < i_end) {
	  
          if (!(nvert[k][j][i] > IBM_FLUID_THRESHOLD || nvert[k][j][i + 1] > IBM_FLUID_THRESHOLD)) {
            // Compute pressure derivatives (dp/d_csi, dp/d_eta, dp/d_zet) at the i-face

	    PetscReal dpdc = p[k][j][i + 1] - p[k][j][i];
            PetscReal dpde = 0.0, dpdz = 0.0;

            // Boundary-aware stencil for dp/d_eta
	      if ((j==my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i]+nvert[k][j+1][i+1] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
		  dpde = (p[k][j][i] + p[k][j][i+1] -
			  p[k][j-1][i] - p[k][j-1][i+1]) * 0.5;
		}
	      }

	      else if ((j==my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i]+nvert[k][j+1][i+1] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1) { dpde = (p[k][j][i] + p[k][j][i+1] - p[k][j-1][i] - p[k][j-1][i+1]) * 0.5; }
	      }

	      else if ((j == 1 && user->bctype[2]!=7) || nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) { dpde = (p[k][j+1][i] + p[k][j+1][i+1] - p[k][j][i] - p[k][j][i+1]) * 0.5; }
	      }

	      else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) { dpde = (p[k][j+1][i] + p[k][j+1][i+1] - p[k][j][i] - p[k][j][i+1]) * 0.5; }
	      }

	      else { dpde = (p[k][j+1][i] + p[k][j+1][i+1] - p[k][j-1][i] - p[k][j-1][i+1]) * 0.25; }

            // Boundary-aware stencil for dp/d_zet
	      if ((k == mz-2 && user->bctype[4]!=7) || nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) { dpdz = (p[k][j][i] + p[k][j][i+1] - p[k-1][j][i] - p[k-1][j][i+1]) * 0.5; }
	      }

	      else  if ((k == mz-2 || k==1) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1) { dpdz = (p[k][j][i] + p[k][j][i+1] - p[k-1][j][i] - p[k-1][j][i+1]) * 0.5; }
	      }

	      else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) { dpdz = (p[k+1][j][i] + p[k+1][j][i+1] - p[k][j][i] - p[k][j][i+1]) * 0.5; }
	      }

	      else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) { dpdz = (p[k+1][j][i] + p[k+1][j][i+1] - p[k][j][i] - p[k][j][i+1]) * 0.5; }
	      }

	      else { dpdz = (p[k+1][j][i] + p[k+1][j][i+1] - p[k-1][j][i] - p[k-1][j][i+1]) * 0.25; }

            // Apply the correction: U_new = U_old - dt * (g11*dpdc + g12*dpde + g13*dpdz)

            
	      
            PetscReal grad_p_x = (dpdc * (icsi[k][j][i].x * icsi[k][j][i].x + icsi[k][j][i].y * icsi[k][j][i].y
					  + icsi[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i] +
                                dpde * (ieta[k][j][i].x * icsi[k][j][i].x + ieta[k][j][i].y * icsi[k][j][i].y
					+ ieta[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i] +
                                dpdz * (izet[k][j][i].x * icsi[k][j][i].x + izet[k][j][i].y * icsi[k][j][i].y
					+ izet[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i]);

	    PetscReal correction = grad_p_x*scale;
	    //LOG_LOOP_ALLOW_EXACT(GLOBAL,LOG_DEBUG,k,5," Flux correction in Csi Direction: %le.\n",correction);
	    ucont[k][j][i].x -= correction;

          }
        }

        // --- Correct V_contravariant (y-component of velocity) ---
        PetscInt j_end = (user->bctype[2] == 7) ? my - 1 : my - 2;
        if (j < j_end) {
          if (!(nvert[k][j][i] > IBM_FLUID_THRESHOLD || nvert[k][j + 1][i] > IBM_FLUID_THRESHOLD)) {
            PetscReal dpdc = 0.0, dpde = 0.0, dpdz = 0.0;
            dpde = p[k][j + 1][i] - p[k][j][i];

            // Boundary-aware stencil for dp/d_csi
	      if ((i == mx-2 && user->bctype[0]!=7) || nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) { dpdc = (p[k][j][i] + p[k][j+1][i] - p[k][j][i-1] - p[k][j+1][i-1]) * 0.5; }
	      } else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1) { dpdc = (p[k][j][i] + p[k][j+1][i] - p[k][j][i-1] - p[k][j+1][i-1]) * 0.5; }
	      } else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) { dpdc = (p[k][j][i+1] + p[k][j+1][i+1] - p[k][j][i] - p[k][j+1][i]) * 0.5; }
	      } else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) { dpdc = (p[k][j][i+1] + p[k][j+1][i+1] - p[k][j][i] - p[k][j+1][i]) * 0.5; }
	      } else { dpdc = (p[k][j][i+1] + p[k][j+1][i+1] - p[k][j][i-1] - p[k][j+1][i-1]) * 0.25; }

            // Boundary-aware stencil for dp/d_zet
	      if ((k == mz-2 && user->bctype[4]!=7)|| nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) { dpdz = (p[k][j][i] + p[k][j+1][i] - p[k-1][j][i] - p[k-1][j+1][i]) * 0.5; }
	      } else if ((k == mz-2 || k==1 ) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1) { dpdz = (p[k][j][i] + p[k][j+1][i] - p[k-1][j][i] - p[k-1][j+1][i]) * 0.5; }
	      } else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) { dpdz = (p[k+1][j][i] + p[k+1][j+1][i] - p[k][j][i] - p[k][j+1][i]) * 0.5; }
	      } else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) { dpdz = (p[k+1][j][i] + p[k+1][j+1][i] - p[k][j][i] - p[k][j+1][i]) * 0.5; }
	      } else { dpdz = (p[k+1][j][i] + p[k+1][j+1][i] - p[k-1][j][i] - p[k-1][j+1][i]) * 0.25; }

            PetscReal grad_p_y = (dpdc * (jcsi[k][j][i].x * jeta[k][j][i].x + jcsi[k][j][i].y * jeta[k][j][i].y + jcsi[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i] +
                                dpde * (jeta[k][j][i].x * jeta[k][j][i].x + jeta[k][j][i].y * jeta[k][j][i].y + jeta[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i] +
                                dpdz * (jzet[k][j][i].x * jeta[k][j][i].x + jzet[k][j][i].y * jeta[k][j][i].y + jzet[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i]);

	    PetscReal correction = grad_p_y*scale;
	    //LOG_LOOP_ALLOW_EXACT(GLOBAL,LOG_DEBUG,k,5," Flux correction in Eta Direction: %le.\n",correction);
            ucont[k][j][i].y -= correction;
          }
        }
        
        // --- Correct W_contravariant (z-component of velocity) ---
        PetscInt k_end = (user->bctype[4] == 7) ? mz - 1 : mz - 2;
        if (k < k_end) {
          if (!(nvert[k][j][i] > IBM_FLUID_THRESHOLD || nvert[k + 1][j][i] > IBM_FLUID_THRESHOLD)) {
            PetscReal dpdc = 0.0, dpde = 0.0, dpdz = 0.0;
            dpdz = p[k + 1][j][i] - p[k][j][i];
            
            // Boundary-aware stencil for dp/d_csi
	      if ((i == mx-2 && user->bctype[0]!=7)|| nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) { dpdc = (p[k][j][i] + p[k+1][j][i] - p[k][j][i-1] - p[k+1][j][i-1]) * 0.5; }
	      } else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1) { dpdc = (p[k][j][i] + p[k+1][j][i] - p[k][j][i-1] - p[k+1][j][i-1]) * 0.5; }
	      } else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) { dpdc = (p[k][j][i+1] + p[k+1][j][i+1] - p[k][j][i] - p[k+1][j][i]) * 0.5; }
	      } else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) { dpdc = (p[k][j][i+1] + p[k+1][j][i+1] - p[k][j][i] - p[k+1][j][i]) * 0.5; }
	      } else { dpdc = (p[k][j][i+1] + p[k+1][j][i+1] - p[k][j][i-1] - p[k+1][j][i-1]) * 0.25; }

            // Boundary-aware stencil for dp/d_eta
	      if ((j == my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) { dpde = (p[k][j][i] + p[k+1][j][i] - p[k][j-1][i] - p[k+1][j-1][i]) * 0.5; }
	      } else  if ((j == my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1) { dpde = (p[k][j][i] + p[k+1][j][i] - p[k][j-1][i] - p[k+1][j-1][i]) * 0.5; }
	      } else if ((j == 1 && user->bctype[2]!=7)|| nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) { dpde = (p[k][j+1][i] + p[k+1][j+1][i] - p[k][j][i] - p[k+1][j][i]) * 0.5; }
	      } else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) { dpde = (p[k][j+1][i] + p[k+1][j+1][i] - p[k][j][i] - p[k+1][j][i]) * 0.5; }
	      } else { dpde = (p[k][j+1][i] + p[k+1][j+1][i] - p[k][j-1][i] - p[k+1][j-1][i]) * 0.25; }

            PetscReal grad_p_z = (dpdc * (kcsi[k][j][i].x * kzet[k][j][i].x + kcsi[k][j][i].y * kzet[k][j][i].y + kcsi[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i] +
                                dpde * (keta[k][j][i].x * kzet[k][j][i].x + keta[k][j][i].y * kzet[k][j][i].y + keta[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i] +
                                dpdz * (kzet[k][j][i].x * kzet[k][j][i].x + kzet[k][j][i].y * kzet[k][j][i].y + kzet[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i]);

	    // ========================= DEBUG PRINT  =========================
	    LOG_LOOP_ALLOW_EXACT(GLOBAL, LOG_DEBUG, k, 5,
				 "[k=%d, j=%d, i=%d] ---- Neighbor Pressures ----\n"
				 "  Central Z-Neighbors: p[k+1][j][i] = %g | p[k][j][i] = %g\n"
				 "  Eta-Stencil (Y-dir): p[k][j-1][i] = %g, p[k+1][j-1][i] = %g | p[k][j+1][i] = %g, p[k+1][j+1][i] = %g\n"
				 "  Csi-Stencil (X-dir): p[k][j][i-1] = %g, p[k+1][j][i-1] = %g | p[k][j][i+1] = %g, p[k+1][j][i+1] = %g\n",
				 k, j, i,
				 p[k + 1][j][i], p[k][j][i],
				 p[k][j - 1][i], p[k + 1][j - 1][i], p[k][j + 1][i], p[k + 1][j + 1][i],
				 p[k][j][i - 1], p[k + 1][j][i - 1], p[k][j][i + 1], p[k + 1][j][i + 1]);
	    // ======================= END DEBUG PRINT =======================

	    LOG_LOOP_ALLOW_EXACT(GLOBAL,LOG_DEBUG,k,5," dpdc: %le | dpde: %le | dpdz: %le.\n",dpdc,dpde,dpdz);	    
	    PetscReal correction = grad_p_z*scale;
	    //LOG_LOOP_ALLOW_EXACT(GLOBAL,LOG_DEBUG,k,5," Flux correction in Zet Direction: %le.\n",correction);
            ucont[k][j][i].z -= correction;
          }
        }
      }
    }
  }

  // --- Explicit correction for periodic boundaries (if necessary) ---
  // The main loop handles the interior, but this handles the first physical layer at periodic boundaries.
  // Note: This logic is largely duplicated from the main loop and could be merged, but is preserved for fidelity.
  if (user->bctype[0] == 7 && xs == 0) {
    for (PetscInt k=lzs; k<lze; k++) {
      for (PetscInt j=lys; j<lye; j++) {
	PetscInt i=xs;
	
	PetscReal dpdc = p[k][j][i+1] - p[k][j][i];
	
	PetscReal dpde = 0.;
	PetscReal dpdz = 0.;
		
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
	
	
	
	if (!(nvert[k][j][i] + nvert[k][j][i+1])) {
	  ucont[k][j][i].x -=
	    (dpdc * (icsi[k][j][i].x * icsi[k][j][i].x +
		     icsi[k][j][i].y * icsi[k][j][i].y +
		     icsi[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i] +
	     dpde * (ieta[k][j][i].x * icsi[k][j][i].x +
		     ieta[k][j][i].y * icsi[k][j][i].y +
		     ieta[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i] +
	     dpdz * (izet[k][j][i].x * icsi[k][j][i].x +
		     izet[k][j][i].y * icsi[k][j][i].y +
		     izet[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i])
	    * scale;
	  
	}
      }
    } 
  }
  if (user->bctype[2] == 7 && ys == 0) {

    for (PetscInt k=lzs; k<lze; k++) {
      for (PetscInt i=lxs; i<lxe; i++) {
	PetscInt j=ys;
	
	PetscReal dpdc = 0.;
	PetscReal dpdz = 0.;
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
	
	PetscReal dpde = p[k][j+1][i] - p[k][j][i];
	
	if ((k == mz-2 && user->bctype[4]!=7)|| nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
	  if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
	    dpdz = (p[k  ][j][i] + p[k  ][j+1][i] -
		    p[k-1][j][i] - p[k-1][j+1][i]) * 0.5;
	  }
	}
	else if ((k == mz-2 || k==1 ) && user->bctype[0]==7 && nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
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
		
	if (!(nvert[k][j][i] + nvert[k][j+1][i])) {
	  ucont[k][j][i].y -=
	    (dpdc * (jcsi[k][j][i].x * jeta[k][j][i].x +
		     jcsi[k][j][i].y * jeta[k][j][i].y +
		     jcsi[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i] +
	     dpde * (jeta[k][j][i].x * jeta[k][j][i].x +
		     jeta[k][j][i].y * jeta[k][j][i].y +
		     jeta[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i] +
	     dpdz * (jzet[k][j][i].x * jeta[k][j][i].x +
		     jzet[k][j][i].y * jeta[k][j][i].y +
		     jzet[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i])
	    * scale;
	}
      }
    }
  }

  if (user->bctype[4] == 7 && zs == 0) {
    for (PetscInt j=lys; j<lye; j++) {
      for (PetscInt i=lxs; i<lxe; i++) {
	
	PetscInt k=zs;
	PetscReal dpdc = 0.;
	PetscReal dpde = 0.;

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
	
	PetscReal dpdz = p[k+1][j][i] - p[k][j][i];
	
	if (!(nvert[k][j][i] + nvert[k+1][j][i])) {
	  
	  ucont[k][j][i].z -=
	   (dpdc * (kcsi[k][j][i].x * kzet[k][j][i].x +
	  		       kcsi[k][j][i].y * kzet[k][j][i].y +
	  		       kcsi[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i] +
	  	       dpde * (keta[k][j][i].x * kzet[k][j][i].x +
	  		       keta[k][j][i].y * kzet[k][j][i].y +
	  		       keta[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i] +
	  	       dpdz * (kzet[k][j][i].x * kzet[k][j][i].x +
	  		       kzet[k][j][i].y * kzet[k][j][i].y +
	  		       kzet[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i])
	                    * scale;
	  
	}
      }
    }
  }
  
  // --- Optional: Enforce specific mass flow rate for channel flow simulations ---
  if (simCtx->channelz) {
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Applying channel flow mass flux correction.\n");
    // DMDAVecGetArray(fda, user->lCsi, &csi); DMDAVecGetArray(fda, user->lEta, &eta); DMDAVecGetArray(fda, user->lZet, &zet);
    //DMDAVecGetArray(da, user->lAj, &aj);
    // This entire block was commented out in the original code. It calculates the current
    // total flux in the z-direction, compares it to a target flux (user->simCtx->FluxInSum),
    // and computes a uniform velocity correction `ratio` to enforce the target.
    // The implementation for applying the `ratio` was also commented out.
    // Preserving this as-is.
  } 

  //================================================================================
  // Section 3: Finalization and Cleanup
  //================================================================================

  // --- Restore access to all PETSc vector arrays ---
  DMDAVecRestoreArray(fda, user->Ucont, &ucont);
  // DMDAVecRestoreArray(fda, user->lCsi, &csi); DMDAVecRestoreArray(fda, user->lEta, &eta); DMDAVecRestoreArray(fda, user->lZet, &zet);
  //DMDAVecRestoreArray(da, user->lAj, &aj);
  DMDAVecRestoreArray(fda, user->lICsi, &icsi); DMDAVecRestoreArray(fda, user->lIEta, &ieta); DMDAVecRestoreArray(fda, user->lIZet, &izet);
  DMDAVecRestoreArray(fda, user->lJCsi, &jcsi); DMDAVecRestoreArray(fda, user->lJEta, &jeta); DMDAVecRestoreArray(fda, user->lJZet, &jzet);
  DMDAVecRestoreArray(fda, user->lKCsi, &kcsi); DMDAVecRestoreArray(fda, user->lKEta, &keta); DMDAVecRestoreArray(fda, user->lKZet, &kzet);
  DMDAVecRestoreArray(da, user->lIAj, &iaj); DMDAVecRestoreArray(da, user->lJAj, &jaj); DMDAVecRestoreArray(da, user->lKAj, &kaj);
  DMDAVecRestoreArray(da, user->lP, &p);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
 
  // --- Update ghost cells for the newly corrected velocity field ---
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Updating ghost cells for corrected velocity.\n");
  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

  // --- Convert velocity to Cartesian and update ghost nodes ---
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Converting velocity to Cartesian and finalizing ghost nodes.\n");
  Contra2Cart(user);
  GhostNodeVelocity(user);

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Exiting Projection step.\n");
  PetscFunctionReturn(0);
}

/**
 * @brief Updates the pressure field and handles periodic boundary conditions.
 *
 * @details
 * This function is a core part of the projection method in a CFD solver. It is
 * called after the Pressure Poisson Equation has been solved to obtain a pressure
 * correction field, `Phi`.
 *
 * The function performs two main operations:
 *
 * 1.  **Core Pressure Update**: It updates the main pressure field `P` by adding the
 *     pressure correction `Phi` at every grid point in the fluid domain. The
 *     operation is `P_new = P_old + Phi`.
 *
 * 2.  **Periodic Boundary Synchronization**: If any of the domain boundaries are
 *     periodic (`bctype == 7`), this function manually updates the pressure values
 *     in the ghost cells. It copies values from the physical cells on one side of
 *     the domain to the corresponding ghost cells on the opposite side. This ensures
 *     that subsequent calculations involving pressure gradients are correct across
 *     periodic boundaries. The refactored version consolidates the original code's
 *     redundant updates into a single, efficient pass.
 *
 * @param[in, out] user Pointer to the UserCtx struct, containing simulation context and
 *                      the pressure vectors `P` and `Phi`.
 *
 * @return PetscErrorCode 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode UpdatePressure(UserCtx *user)
{
  PetscFunctionBeginUser;
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Entering UpdatePressure.\n");

  //================================================================================
  // Section 1: Initialization and Data Acquisition
  //================================================================================
  DM da = user->da;
  DMDALocalInfo info = user->info;

  // Local grid extents for the main update loop
  PetscInt xs = info.xs, xe = info.xs + info.xm;
  PetscInt ys = info.ys, ye = info.ys + info.ym;
  PetscInt zs = info.zs, ze = info.zs + info.zm;
  PetscInt mx = info.mx, my = info.my, mz = info.mz;

  PetscInt i,j,k;

  // --- Get direct pointer access to PETSc vector data for performance ---
  PetscReal ***p, ***phi;
  DMDAVecGetArray(da, user->P, &p);
  DMDAVecGetArray(da, user->Phi, &phi);

  //================================================================================
  // Section 2: Core Pressure Update
  //================================================================================
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Performing core pressure update (P_new = P_old + Phi).\n");
  for (PetscInt k = zs; k < ze; k++) {
    for (PetscInt j = ys; j < ye; j++) {
      for (PetscInt i = xs; i < xe; i++) {
        // This is the fundamental pressure update in a projection method.
        p[k][j][i] += phi[k][j][i];
      }
    }
  }

  // Restore arrays now that the core computation is done.
  DMDAVecRestoreArray(da, user->Phi, &phi);
  DMDAVecRestoreArray(da, user->P, &p);

  
  //================================================================================
  // Section 3: Handle Periodic Boundary Condition Synchronization
  //================================================================================
  // This block is executed only if at least one boundary is periodic.
  // The original code contained many redundant Get/Restore and update calls.
  // This refactored version performs the same effective logic but in a single,
  // efficient pass, which is numerically equivalent and much cleaner.
  if (user->bctype[0] == 7 || user->bctype[1] == 7 ||
      user->bctype[2] == 7 || user->bctype[3] == 7 ||
      user->bctype[4] == 7 || user->bctype[5] == 7)
  {
    LOG_ALLOW(GLOBAL, LOG_DEBUG, "Synchronizing ghost cells for periodic boundaries.\n");

    // First, update the local vectors (including ghost regions) with the new global data.
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalBegin(da, user->Phi, INSERT_VALUES, user->lPhi);
    DMGlobalToLocalEnd(da, user->Phi, INSERT_VALUES, user->lPhi);

    // Get pointers to all necessary local and global arrays ONCE.
    PetscReal ***lp, ***lphi;
    DMDAVecGetArray(da, user->lP, &lp);
    DMDAVecGetArray(da, user->lPhi, &lphi);
    DMDAVecGetArray(da, user->P, &p);
    DMDAVecGetArray(da, user->Phi, &phi);

    // --- X-Direction Periodic Update ---
    if (user->bctype[0] == 7 || user->bctype[1] == 7) {
      // Update left boundary physical cells from right boundary ghost cells
      if (xs == 0) {
        PetscInt i = 0;
        for (k = zs; k < ze; k++) {
          for (j = ys; j < ye; j++) {
            // Note: Accessing lp[...][i-2] reads from the ghost cell region.
            p[k][j][i] = lp[k][j][i - 2];
            phi[k][j][i] = lphi[k][j][i - 2];
          }
        }
      }
      // Update right boundary physical cells from left boundary ghost cells
      if (xe == mx) {
        PetscInt i = mx - 1;
        for (k = zs; k < ze; k++) {
          for (j = ys; j < ye; j++) {
            p[k][j][i] = lp[k][j][i + 2];
            phi[k][j][i] = lphi[k][j][i + 2];
          }
        }
      }
    }

    // --- Y-Direction Periodic Update ---
    if (user->bctype[2] == 7 || user->bctype[3] == 7) {
      // Update bottom boundary
      if (ys == 0) {
        PetscInt j = 0;
        for (k = zs; k < ze; k++) {
          for (i = xs; i < xe; i++) {
            p[k][j][i] = lp[k][j - 2][i];
            phi[k][j][i] = lphi[k][j - 2][i];
          }
        }
      }
      // Update top boundary
      if (ye == my) {
        PetscInt j = my - 1;
        for (k = zs; k < ze; k++) {
          for (i = xs; i < xe; i++) {
            p[k][j][i] = lp[k][j + 2][i];
            phi[k][j][i] = lphi[k][j + 2][i];
          }
        }
      }
    }

    // --- Z-Direction Periodic Update ---
    if (user->bctype[4] == 7 || user->bctype[5] == 7) {
      // Update front boundary
      if (zs == 0) {
        PetscInt k = 0;
        for (j = ys; j < ye; j++) {
          for (i = xs; i < xe; i++) {
            p[k][j][i] = lp[k - 2][j][i];
            phi[k][j][i] = lphi[k - 2][j][i];
          }
        }
      }
      // Update back boundary
      if (ze == mz) {
        PetscInt k = mz - 1;
        for (j = ys; j < ye; j++) {
          for (i = xs; i < xe; i++) {
            p[k][j][i] = lp[k + 2][j][i];
            phi[k][j][i] = lphi[k + 2][j][i];
          }
        }
      }
    }
    
    // Restore all arrays ONCE.
    DMDAVecRestoreArray(da, user->lP, &lp);
    DMDAVecRestoreArray(da, user->lPhi, &lphi);
    DMDAVecRestoreArray(da, user->P, &p);
    DMDAVecRestoreArray(da, user->Phi, &phi);

    // After manually updating the physical boundary cells, we must call
    // DMGlobalToLocal again to ensure all processes have the updated ghost
    // values for the *next* function that needs them.
    DMGlobalToLocalBegin(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalEnd(da, user->P, INSERT_VALUES, user->lP);
    DMGlobalToLocalBegin(da, user->Phi, INSERT_VALUES, user->lPhi);
    DMGlobalToLocalEnd(da, user->Phi, INSERT_VALUES, user->lPhi);
  }
  
  //================================================================================
  // Section 4: Final Cleanup (pointers already restored)
  //================================================================================

  UpdateLocalGhosts(user,"P");
  UpdateLocalGhosts(user,"Phi");
  
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Exiting UpdatePressure.\n");
  PetscFunctionReturn(0);
}

static PetscErrorCode mymatmultadd(Mat mat, Vec v1, Vec v2)
{

  Vec vt;
  VecDuplicate(v2, &vt);
  MatMult(mat, v1, vt);
  VecAYPX(v2, 1., vt);
  VecDestroy(&vt);
  return(0);
}

PetscErrorCode PoissonNullSpaceFunction(MatNullSpace nullsp,Vec X, void *ctx)
{
  UserCtx *user = (UserCtx*)ctx;

  DM da = user->da;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal	***x, ***nvert;
  PetscInt	i, j, k;

/*   /\* First remove a constant from the Vec field X *\/ */


  /* Then apply boundary conditions */
  DMDAVecGetArray(da, X, &x);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscReal lsum, sum;
  PetscReal  lnum, num;

  if (user->multinullspace) PetscPrintf(PETSC_COMM_WORLD, "MultiNullSpace!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  if (!user->multinullspace) {
    lsum = 0;
    lnum = 0;
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k][j][i] < 0.1) {
	    lsum += x[k][j][i];
	    lnum ++;
	  }
	}
      }
    }

    MPI_Allreduce(&lsum,&sum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lnum,&num,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
   /*  PetscGlobalSum(&lsum, &sum, PETSC_COMM_WORLD); */
/*     PetscGlobalSum(&lnum, &num, PETSC_COMM_WORLD); */
    sum = sum / (-1.0 * num);

    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  if (nvert[k][j][i] < 0.1) {
	    x[k][j][i] +=sum;
	  }
	}
      }
    }
  }
  else {
    lsum = 0;
    lnum = 0;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	for (k=lzs; k<lze; k++) {
	  if (k<user->KSKE[2*(j*mx+i)] && nvert[k][j][i]<0.1) {
	    lsum += x[k][j][i];
	    lnum ++;
	  }
	}
      }
    }
    MPI_Allreduce(&lsum,&sum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lnum,&num,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
   /*  PetscGlobalSum(&lsum, &sum, PETSC_COMM_WORLD); */
/*     PetscGlobalSum(&lnum, &num, PETSC_COMM_WORLD); */
    sum /= -num;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	for (k=lzs; k<lze; k++) {
	  if (k<user->KSKE[2*(j*mx+i)] && nvert[k][j][i]<0.1) {
	    x[k][j][i] += sum;
	  }
	}
      }
    }

    lsum = 0;
    lnum = 0;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	for (k=lzs; k<lze; k++) {
	  if (k>=user->KSKE[2*(j*mx+i)] && nvert[k][j][i]<0.1) {
	    lsum += x[k][j][i];
	    lnum ++;
	  }
	}
      }
    }
    MPI_Allreduce(&lsum,&sum,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(&lnum,&num,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
   /*  PetscGlobalSum(&lsum, &sum, PETSC_COMM_WORLD); */
/*     PetscGlobalSum(&lnum, &num, PETSC_COMM_WORLD); */
    sum /= -num;
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	for (k=lzs; k<lze; k++) {
	  if (k>=user->KSKE[2*(j*mx+i)] && nvert[k][j][i]<0.1) {
	    x[k][j][i] += sum;
	  }
	}
      }
    }
			   
  } //if multinullspace
  if (zs == 0) {
    k = 0;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (ze == mz) {
    k = mz-1;
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (ys == 0) {
    j = 0;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (ye == my) {
    j = my-1;
    for (k=zs; k<ze; k++) {
      for (i=xs; i<xe; i++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (xs == 0) {
    i = 0;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	x[k][j][i] = 0.;
      }
    }
  }

  if (xe == mx) {
    i = mx-1;
    for (k=zs; k<ze; k++) {
      for (j=ys; j<ye; j++) {
	x[k][j][i] = 0.;
      }
    }
  }

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (nvert[k][j][i] > 0.1)
	  x[k][j][i] = 0.;
      }
    }
  }
  DMDAVecRestoreArray(da, X, &x);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);

  return 0;
}

PetscErrorCode MyInterpolation(Mat A, Vec X, Vec F)
{
  UserCtx *user;

  MatShellGetContext(A, (void**)&user);

  
  
  DM	da = user->da;

  DM	da_c = *user->da_c;

  DMDALocalInfo	info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  
  PetscReal ***f, ***x, ***nvert, ***nvert_c;
  PetscInt i, j, k, ic, jc, kc, ia, ja, ka;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;


  DMDAVecGetArray(da,   F, &f);


  Vec lX;
  DMCreateLocalVector(da_c, &lX);
 
  DMGlobalToLocalBegin(da_c, X, INSERT_VALUES, lX);
  DMGlobalToLocalEnd(da_c, X, INSERT_VALUES, lX);  
  DMDAVecGetArray(da_c, lX, &x);

  DMDAVecGetArray(da, user->lNvert, &nvert);
  DMDAVecGetArray(da_c, *(user->lNvert_c), &nvert_c);
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {

	GridInterpolation(i, j, k, ic, jc, kc, ia, ja, ka, user);

	  f[k][j][i] = (x[kc   ][jc   ][ic   ] * 9 +
			x[kc   ][jc+ja][ic   ] * 3 +
			x[kc   ][jc   ][ic+ia] * 3 +
			x[kc   ][jc+ja][ic+ia]) * 3./64. +
	    (x[kc+ka][jc   ][ic   ] * 9 +
	     x[kc+ka][jc+ja][ic   ] * 3 +
	     x[kc+ka][jc   ][ic+ia] * 3 +
	     x[kc+ka][jc+ja][ic+ia]) /64.;
      }
    }
  }

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {

	if (i==0) {
	  f[k][j][i] = 0.;//-f[k][j][i+1];
	}
	else if (i==mx-1) {
	  f[k][j][i] = 0.;//-f[k][j][i-1];
	}
	else if (j==0) {
	  f[k][j][i] = 0.;//-f[k][j+1][i];
	}
	else if (j==my-1) {
	  f[k][j][i] = 0.;//-f[k][j-1][i];
	}
	else if (k==0) {
	  f[k][j][i] = 0.;//-f[k+1][j][i];
	}
	else if (k==mz-1) {
	  f[k][j][i] = 0.;//-f[k-1][j][i];
	}
	if (nvert[k][j][i] > 0.1) f[k][j][i] = 0.;

      }
    }
  }

  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(da_c, *(user->lNvert_c), &nvert_c);

  DMDAVecRestoreArray(da_c, lX, &x);
 
  VecDestroy(&lX);
  DMDAVecRestoreArray(da,   F,  &f);



  return 0;
 
}

static PetscErrorCode RestrictResidual_SolidAware(Mat A, Vec X, Vec F)
{
  UserCtx *user;
  MatShellGetContext(A, (void**)&user);
  
  DM	da = user->da;
  DM	da_f = *user->da_f;

  DMDALocalInfo	info;
  DMDAGetLocalInfo(da, &info);
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  
  PetscReal ***f, ***x, ***nvert;
  PetscInt i, j, k, ih, jh, kh, ia, ja, ka;

  DMDAVecGetArray(da,   F, &f);

  Vec lX;
  DMCreateLocalVector(da_f, &lX);
  DMGlobalToLocalBegin(da_f, X, INSERT_VALUES, lX);
  DMGlobalToLocalEnd(da_f, X, INSERT_VALUES, lX);  
  DMDAVecGetArray(da_f, lX, &x);

  DMDAVecGetArray(da, user->lNvert, &nvert);

  PetscReal ***nvert_f;
  DMDAVecGetArray(da_f, user->user_f->lNvert, &nvert_f);

  if ((user->isc)) ia = 0;
  else ia = 1;

  if ((user->jsc)) ja = 0;
  else ja = 1;

  if ((user->ksc)) ka = 0;
  else ka = 1;

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
        // --- CORRECTED LOGIC ---
        // First, check if the current point is a boundary point.
        // If it is, it does not contribute to the coarse grid residual.
        if (i==0 || i==mx-1 || j==0 || j==my-1 || k==0 || k==mz-1 || nvert[k][j][i] > 0.1) {
            f[k][j][i] = 0.0;
        } 
        // Only if it's a true interior fluid point, perform the restriction.
        else {
            GridRestriction(i, j, k, &ih, &jh, &kh, user);
            f[k][j][i] = 0.125 *
              (x[kh   ][jh   ][ih   ] * PetscMax(0., 1 - nvert_f[kh   ][jh   ][ih   ]) +
               x[kh   ][jh   ][ih-ia] * PetscMax(0., 1 - nvert_f[kh   ][jh   ][ih-ia]) +
               x[kh   ][jh-ja][ih   ] * PetscMax(0., 1 - nvert_f[kh   ][jh-ja][ih   ]) +
               x[kh-ka][jh   ][ih   ] * PetscMax(0., 1 - nvert_f[kh-ka][jh   ][ih   ]) +
               x[kh   ][jh-ja][ih-ia] * PetscMax(0., 1 - nvert_f[kh   ][jh-ja][ih-ia]) +
               x[kh-ka][jh-ja][ih   ] * PetscMax(0., 1 - nvert_f[kh-ka][jh-ja][ih   ]) +
               x[kh-ka][jh   ][ih-ia] * PetscMax(0., 1 - nvert_f[kh-ka][jh   ][ih-ia]) +
               x[kh-ka][jh-ja][ih-ia] * PetscMax(0., 1 - nvert_f[kh-ka][jh-ja][ih-ia]));
        }
      }
    }
  }

  DMDAVecRestoreArray(da_f, user->user_f->lNvert, &nvert_f);
  DMDAVecRestoreArray(da_f, lX, &x);
  VecDestroy(&lX);
  DMDAVecRestoreArray(da,   F,  &f);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);

  return 0;
}

PetscErrorCode MyRestriction(Mat A, Vec X, Vec F)
{
  UserCtx *user;

  MatShellGetContext(A, (void**)&user);

  
  DM	da = user->da;

  DM	da_f = *user->da_f;

  DMDALocalInfo	info;
  DMDAGetLocalInfo(da, &info);
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
  //  PetscInt	lxs, lxe, lys, lye, lzs, lze;
  
  PetscReal ***f, ***x, ***nvert;
  PetscInt i, j, k, ih, jh, kh, ia, ja, ka;

  DMDAVecGetArray(da,   F, &f);

  Vec lX;
 
  DMCreateLocalVector(da_f, &lX);
  DMGlobalToLocalBegin(da_f, X, INSERT_VALUES, lX);
  DMGlobalToLocalEnd(da_f, X, INSERT_VALUES, lX);  
  DMDAVecGetArray(da_f, lX, &x);

  DMDAVecGetArray(da, user->lNvert, &nvert);

  PetscReal ***nvert_f;
  DMDAVecGetArray(da_f, user->user_f->lNvert, &nvert_f);

  if ((user->isc)) ia = 0;
  else ia = 1;

  if ((user->jsc)) ja = 0;
  else ja = 1;

  if ((user->ksc)) ka = 0;
  else ka = 1;

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	if (k==0) {
	  f[k][j][i] = 0.;
	}
	else if (k==mz-1) {
	  f[k][j][i] = 0.;
	}
	else if (j==0) {
	  f[k][j][i] = 0.;
	}
	else if (j==my-1) {
	  f[k][j][i] = 0.;
	}
	else if (i==0) {
	  f[k][j][i] = 0.;
	}
	else if (i==mx-1) {
	  f[k][j][i] = 0.;
	}
	else {
	  GridRestriction(i, j, k, &ih, &jh, &kh, user);
	  f[k][j][i] = 0.125 *
	    (x[kh   ][jh   ][ih   ] * PetscMax(0., 1 - nvert_f[kh   ][jh   ][ih   ]) +
	     x[kh   ][jh   ][ih-ia] * PetscMax(0., 1 - nvert_f[kh   ][jh   ][ih-ia]) +
	     x[kh   ][jh-ja][ih   ] * PetscMax(0., 1 - nvert_f[kh   ][jh-ja][ih   ]) +
	     x[kh-ka][jh   ][ih   ] * PetscMax(0., 1 - nvert_f[kh-ka][jh   ][ih   ]) +
	     x[kh   ][jh-ja][ih-ia] * PetscMax(0., 1 - nvert_f[kh   ][jh-ja][ih-ia]) +
	     x[kh-ka][jh-ja][ih   ] * PetscMax(0., 1 - nvert_f[kh-ka][jh-ja][ih   ]) +
	     x[kh-ka][jh   ][ih-ia] * PetscMax(0., 1 - nvert_f[kh-ka][jh   ][ih-ia]) +
	     x[kh-ka][jh-ja][ih-ia] * PetscMax(0., 1 - nvert_f[kh-ka][jh-ja][ih-ia]));



	  if (nvert[k][j][i] > 0.1) f[k][j][i] = 0.;
	}
      }
    }
  }


  DMDAVecRestoreArray(da_f, user->user_f->lNvert, &nvert_f);

  DMDAVecRestoreArray(da_f, lX, &x);
  VecDestroy(&lX);
 
  DMDAVecRestoreArray(da,   F,  &f);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);


  return 0;
}

/**
 * @brief Assembles the Left-Hand Side (LHS) matrix for the Pressure Poisson Equation.
 *
 * @details
 * This function constructs the sparse matrix `A` (the LHS) for the linear system `Ax = B`,
 * which is the Pressure Poisson Equation (PPE). The matrix `A` represents a discrete
 * version of the negative Laplacian operator (-∇²), tailored for a general curvilinear,
 * staggered grid.
 *
 * The assembly process is highly complex and follows these main steps:
 *
 * 1.  **Matrix Initialization**: On the first call, it allocates a sparse PETSc matrix `A`
 *     pre-allocating space for a 19-point stencil per row. On subsequent calls, it
 *     simply zeroes out the existing matrix.
 *
 * 2.  **Metric Tensor Calculation**: It first computes the 9 components of the metric
 *     tensor (`g11`, `g12`, ..., `g33`) at the cell faces. These `g_ij` coefficients
 *     are essential for defining the Laplacian operator in generalized curvilinear
 *     coordinates and account for grid stretching and non-orthogonality.
 *
 * 3.  **Stencil Assembly Loop**: The function iterates through every grid point `(i,j,k)`.
 *     For each point, it determines the matrix row entries based on its status:
 *     - **Boundary/Solid Point**: If the point is on a domain boundary or inside an
 *       immersed solid (`nvert > 0.1`), it sets an identity row (`A(row,row) = 1`),
 *       effectively removing it from the linear solve.
 *     - **Fluid Point**: For a fluid point, it computes the 19 coefficients of the
 *       finite volume stencil. This involves summing contributions from each of the
 *       six faces of the control volume around the point.
 *
 * 4.  **Boundary-Aware Stencils**: The stencil calculation is critically dependent on the
 *     state of neighboring cells. The code contains intricate logic to check if neighbors
 *     are fluid or solid. If a standard centered-difference stencil would cross into a
 *     solid, the scheme is automatically adapted to a one-sided difference to maintain
 *     accuracy at the fluid-solid interface.
 *
 * 5.  **Matrix Finalization**: After all rows corresponding to the local processor's
 *     grid points have been set, `MatAssemblyBegin/End` is called to finalize the
 *     global matrix, making it ready for use by a linear solver.
 *
 * @param[in, out] user Pointer to the UserCtx struct, containing all simulation context,
 *                      grid data, and the matrix `A`.
 *
 * @return PetscErrorCode 0 on success, or a non-zero error code on failure.
 */
PetscErrorCode PoissonLHSNew(UserCtx *user)
{
  PetscFunctionBeginUser;
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Entering PoissonLHSNew to assemble Laplacian matrix.\n");
  PetscErrorCode ierr;
  //================================================================================
  // Section 1: Initialization and Data Acquisition
  //================================================================================

  
  // --- Get simulation and grid context ---
  SimCtx *simCtx = user->simCtx;
  DM da = user->da, fda = user->fda;
  DMDALocalInfo info = user->info;
  PetscInt IM = user->IM, JM = user->JM, KM = user->KM;
  PetscInt i,j,k;

  // --- Grid dimensions ---
  PetscInt mx = info.mx, my = info.my, mz = info.mz;
  PetscInt xs = info.xs, xe = info.xs + info.xm;
  PetscInt ys = info.ys, ye = info.ys + info.ym;
  PetscInt zs = info.zs, ze = info.zs + info.zm;
  PetscInt gxs = info.gxs, gxe = gxs + info.gxm;
  PetscInt gys = info.gys, gye = gys + info.gym;
  PetscInt gzs = info.gzs, gze = gzs + info.gzm;

  // --- Define constants for clarity ---
  const PetscReal IBM_FLUID_THRESHOLD = 0.1;

  // --- Allocate the LHS matrix A on the first call ---
  if (!user->assignedA) {
    LOG_ALLOW(GLOBAL, LOG_INFO, "First call: Creating LHS matrix 'A' with 19-point stencil preallocation.\n");
    PetscInt N = mx * my * mz; // Total size
    PetscInt M;                // Local size
    VecGetLocalSize(user->Phi, &M);
    // Create a sparse AIJ matrix, preallocating for 19 non-zeros per row (d=diagonal, o=off-diagonal)
    MatCreateAIJ(PETSC_COMM_WORLD, M, M, N, N, 19, PETSC_NULL, 19, PETSC_NULL, &(user->A));
    user->assignedA = PETSC_TRUE;
  }

  // Zero out matrix entries from the previous solve
  MatZeroEntries(user->A);

  // --- Get direct pointer access to grid metric data ---
  Cmpnts ***csi, ***eta, ***zet, ***icsi, ***ieta, ***izet, ***jcsi, ***jeta, ***jzet, ***kcsi, ***keta, ***kzet;
  PetscReal ***aj, ***iaj, ***jaj, ***kaj, ***nvert;
  DMDAVecGetArray(fda, user->lCsi, &csi); DMDAVecGetArray(fda, user->lEta, &eta); DMDAVecGetArray(fda, user->lZet, &zet);
  DMDAVecGetArray(fda, user->lICsi, &icsi); DMDAVecGetArray(fda, user->lIEta, &ieta); DMDAVecGetArray(fda, user->lIZet, &izet);
  DMDAVecGetArray(fda, user->lJCsi, &jcsi); DMDAVecGetArray(fda, user->lJEta, &jeta); DMDAVecGetArray(fda, user->lJZet, &jzet);
  DMDAVecGetArray(fda, user->lKCsi, &kcsi); DMDAVecGetArray(fda, user->lKEta, &keta); DMDAVecGetArray(fda, user->lKZet, &kzet);
  DMDAVecGetArray(da, user->lAj, &aj); DMDAVecGetArray(da, user->lIAj, &iaj); DMDAVecGetArray(da, user->lJAj, &jaj); DMDAVecGetArray(da, user->lKAj, &kaj);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  // --- Create temporary vectors for the metric tensor components G_ij ---
  Vec G11, G12, G13, G21, G22, G23, G31, G32, G33;
  PetscReal ***g11, ***g12, ***g13, ***g21, ***g22, ***g23, ***g31, ***g32, ***g33;
  VecDuplicate(user->lAj, &G11); VecDuplicate(user->lAj, &G12); VecDuplicate(user->lAj, &G13);
  VecDuplicate(user->lAj, &G21); VecDuplicate(user->lAj, &G22); VecDuplicate(user->lAj, &G23);
  VecDuplicate(user->lAj, &G31); VecDuplicate(user->lAj, &G32); VecDuplicate(user->lAj, &G33);
  DMDAVecGetArray(da, G11, &g11); DMDAVecGetArray(da, G12, &g12); DMDAVecGetArray(da, G13, &g13);
  DMDAVecGetArray(da, G21, &g21); DMDAVecGetArray(da, G22, &g22); DMDAVecGetArray(da, G23, &g23);
  DMDAVecGetArray(da, G31, &g31); DMDAVecGetArray(da, G32, &g32); DMDAVecGetArray(da, G33, &g33);

  //================================================================================
  // Section 2: Pre-compute Metric Tensor Coefficients (g_ij)
  //================================================================================
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Pre-computing metric tensor components (g_ij).\n");
  for (k = gzs; k < gze; k++) {
    for (j = gys; j < gye; j++) {
      for (i = gxs; i < gxe; i++) {
        // These coefficients represent the dot products of the grid's contravariant base vectors,
        // scaled by face area. They are the core of the Laplacian operator on a curvilinear grid.
        if(i>-1 && j>-1 && k>-1 && i<IM+1 && j<JM+1 && k<KM+1){
            g11[k][j][i] = (icsi[k][j][i].x * icsi[k][j][i].x + icsi[k][j][i].y * icsi[k][j][i].y + icsi[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i];
            g12[k][j][i] = (ieta[k][j][i].x * icsi[k][j][i].x + ieta[k][j][i].y * icsi[k][j][i].y + ieta[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i];
            g13[k][j][i] = (izet[k][j][i].x * icsi[k][j][i].x + izet[k][j][i].y * icsi[k][j][i].y + izet[k][j][i].z * icsi[k][j][i].z) * iaj[k][j][i];
            g21[k][j][i] = (jcsi[k][j][i].x * jeta[k][j][i].x + jcsi[k][j][i].y * jeta[k][j][i].y + jcsi[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i];
            g22[k][j][i] = (jeta[k][j][i].x * jeta[k][j][i].x + jeta[k][j][i].y * jeta[k][j][i].y + jeta[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i];
            g23[k][j][i] = (jzet[k][j][i].x * jeta[k][j][i].x + jzet[k][j][i].y * jeta[k][j][i].y + jzet[k][j][i].z * jeta[k][j][i].z) * jaj[k][j][i];
            g31[k][j][i] = (kcsi[k][j][i].x * kzet[k][j][i].x + kcsi[k][j][i].y * kzet[k][j][i].y + kcsi[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i];
            g32[k][j][i] = (keta[k][j][i].x * kzet[k][j][i].x + keta[k][j][i].y * kzet[k][j][i].y + keta[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i];
            g33[k][j][i] = (kzet[k][j][i].x * kzet[k][j][i].x + kzet[k][j][i].y * kzet[k][j][i].y + kzet[k][j][i].z * kzet[k][j][i].z) * kaj[k][j][i];
        }
      }
    }
  }

  //================================================================================
  // Section 3: Assemble the LHS Matrix A
  //================================================================================
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Assembling the LHS matrix A using a 19-point stencil.\n");

  // --- Define domain boundaries for stencil logic, accounting for periodic BCs ---
  PetscInt x_str, x_end, y_str, y_end, z_str, z_end;
  if (user->bctype[0] == 7) { x_end = mx - 1; x_str = 0; }
  else { x_end = mx - 2; x_str = 1; }
  if (user->bctype[2] == 7) { y_end = my - 1; y_str = 0; }
  else { y_end = my - 2; y_str = 1; }
  if (user->bctype[4] == 7) { z_end = mz - 1; z_str = 0; }
  else { z_end = mz - 2; z_str = 1; }

  // --- Main assembly loop over all local grid points ---
  for (k = zs; k < ze; k++) {
    for (j = ys; j < ye; j++) {
      for (i = xs; i < xe; i++) {
        PetscScalar vol[19]; // Holds the 19 stencil coefficient values for the current row
        PetscInt idx[19];    // Holds the 19 global column indices for the current row
        PetscInt row = Gidx(i, j, k, user); // Global index for the current row

        // --- Handle Domain Boundary and Immersed Solid Points ---
        // For these points, we don't solve the Poisson equation. We set an identity
        // row (A_ii = 1) to effectively fix the pressure value (usually to 0).
        if (i == 0 || i == mx - 1 || j == 0 || j == my - 1 || k == 0 || k == mz - 1 || nvert[k][j][i] > IBM_FLUID_THRESHOLD) {
          vol[CP] = 1.0;
          idx[CP] = row;
          MatSetValues(user->A, 1, &row, 1, &idx[CP], &vol[CP], INSERT_VALUES);
        }
        // --- Handle Fluid Points ---
        else {
          for (PetscInt m = 0; m < 19; m++) {
            vol[m] = 0.0;
          }

          /************************************************************************
           * EAST FACE CONTRIBUTION (between i and i+1)
           ************************************************************************/
          if (nvert[k][j][i + 1] < IBM_FLUID_THRESHOLD && i != x_end) { // East neighbor is fluid
            // Primary derivative term: d/d_csi (g11 * dP/d_csi)
            vol[CP] -= g11[k][j][i];
            vol[EP] += g11[k][j][i];

            // Cross-derivative term: d/d_csi (g12 * dP/d_eta).
            // This requires an average of dP/d_eta. If a neighbor is solid, the stencil
            // dynamically switches to a one-sided difference to avoid using solid points.
	      if ((j == my-2 && user->bctype[2]!=7) || nvert[k][j+1][i] + nvert[k][j+1][i+1] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
		  vol[CP] += g12[k][j][i] * 0.5; vol[EP] += g12[k][j][i] * 0.5;
		  vol[SP] -= g12[k][j][i] * 0.5; vol[SE] -= g12[k][j][i] * 0.5;
		}
	      }
	      else if ((j == my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i] + nvert[k][j+1][i+1] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k][j-1][i+1] < 0.1) {
		  vol[CP] += g12[k][j][i] * 0.5; vol[EP] += g12[k][j][i] * 0.5;
		  vol[SP] -= g12[k][j][i] * 0.5; vol[SE] -= g12[k][j][i] * 0.5;
		}
	      }
	      else if ((j == 1 && user->bctype[2]!=7) || nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
		  vol[NP] += g12[k][j][i] * 0.5; vol[NE] += g12[k][j][i] * 0.5;
		  vol[CP] -= g12[k][j][i] * 0.5; vol[EP] -= g12[k][j][i] * 0.5;
		}
	      }
	      else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k][j-1][i+1] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k][j+1][i+1] < 0.1) {
		  vol[NP] += g12[k][j][i] * 0.5; vol[NE] += g12[k][j][i] * 0.5;
		  vol[CP] -= g12[k][j][i] * 0.5; vol[EP] -= g12[k][j][i] * 0.5;
		}
	      }
	      else { // Centered difference
		vol[NP] += g12[k][j][i] * 0.25; vol[NE] += g12[k][j][i] * 0.25;
		vol[SP] -= g12[k][j][i] * 0.25; vol[SE] -= g12[k][j][i] * 0.25;
	      }
            
            // Cross-derivative term: d/d_csi (g13 * dP/d_zet)
	      if ((k == mz-2 && user->bctype[4]!=7) || nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
		  vol[CP] += g13[k][j][i] * 0.5; vol[EP] += g13[k][j][i] * 0.5;
		  vol[BP] -= g13[k][j][i] * 0.5; vol[BE] -= g13[k][j][i] * 0.5;
		}
	      }
	      else if ((k == mz-2 || k==1) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j][i+1] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j][i+1] < 0.1) {
		  vol[CP] += g13[k][j][i] * 0.5; vol[EP] += g13[k][j][i] * 0.5;
		  vol[BP] -= g13[k][j][i] * 0.5; vol[BE] -= g13[k][j][i] * 0.5;
		}
	      }
	      else if ((k == 1 && user->bctype[4]!=7) || nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
		  vol[TP] += g13[k][j][i] * 0.5; vol[TE] += g13[k][j][i] * 0.5;
		  vol[CP] -= g13[k][j][i] * 0.5; vol[EP] -= g13[k][j][i] * 0.5;
		}
	      }
	      else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j][i+1] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j][i+1] < 0.1) {
		  vol[TP] += g13[k][j][i] * 0.5; vol[TE] += g13[k][j][i] * 0.5;
		  vol[CP] -= g13[k][j][i] * 0.5; vol[EP] -= g13[k][j][i] * 0.5;
		}
	      }
	      else { // Centered difference
		vol[TP] += g13[k][j][i] * 0.25; vol[TE] += g13[k][j][i] * 0.25;
		vol[BP] -= g13[k][j][i] * 0.25; vol[BE] -= g13[k][j][i] * 0.25;
	      }
          }

          /************************************************************************
           * WEST FACE CONTRIBUTION (between i-1 and i)
           ************************************************************************/
          if (nvert[k][j][i-1] < IBM_FLUID_THRESHOLD && i != x_str) {
	      vol[CP] -= g11[k][j][i-1];
	      vol[WP] += g11[k][j][i-1];

	      if ((j == my-2 && user->bctype[2]!=7) || nvert[k][j+1][i] + nvert[k][j+1][i-1] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k][j-1][i-1] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
		  vol[CP] -= g12[k][j][i-1] * 0.5; vol[WP] -= g12[k][j][i-1] * 0.5;
		  vol[SP] += g12[k][j][i-1] * 0.5; vol[SW] += g12[k][j][i-1] * 0.5;
		}
	      }
	      else if ((j == my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i] + nvert[k][j+1][i-1] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k][j-1][i-1] < 0.1) {
		  vol[CP] -= g12[k][j][i-1] * 0.5; vol[WP] -= g12[k][j][i-1] * 0.5;
		  vol[SP] += g12[k][j][i-1] * 0.5; vol[SW] += g12[k][j][i-1] * 0.5;
		}
	      }
	      else if ((j == 1 && user->bctype[2]!=7)|| nvert[k][j-1][i] + nvert[k][j-1][i-1] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k][j+1][i-1] < 0.1) {
		  vol[NP] -= g12[k][j][i-1] * 0.5; vol[NW] -= g12[k][j][i-1] * 0.5;
		  vol[CP] += g12[k][j][i-1] * 0.5; vol[WP] += g12[k][j][i-1] * 0.5;
		}
	      }
	      else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k][j-1][i-1] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k][j+1][i-1] < 0.1) {
		  vol[NP] -= g12[k][j][i-1] * 0.5; vol[NW] -= g12[k][j][i-1] * 0.5;
		  vol[CP] += g12[k][j][i-1] * 0.5; vol[WP] += g12[k][j][i-1] * 0.5;
		}
	      }
	      else {
		vol[NP] -= g12[k][j][i-1] * 0.25; vol[NW] -= g12[k][j][i-1] * 0.25;
		vol[SP] += g12[k][j][i-1] * 0.25; vol[SW] += g12[k][j][i-1] * 0.25;
	      }

	      if ((k == mz-2 && user->bctype[4]!=7) || nvert[k+1][j][i] + nvert[k+1][j][i-1] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j][i-1] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
		  vol[CP] -= g13[k][j][i-1] * 0.5; vol[WP] -= g13[k][j][i-1] * 0.5;
		  vol[BP] += g13[k][j][i-1] * 0.5; vol[BW] += g13[k][j][i-1] * 0.5;
		}
	      }
	      else if ((k == mz-2 || k==1) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j][i-1] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j][i-1] < 0.1) {
		  vol[CP] -= g13[k][j][i-1] * 0.5; vol[WP] -= g13[k][j][i-1] * 0.5;
		  vol[BP] += g13[k][j][i-1] * 0.5; vol[BW] += g13[k][j][i-1] * 0.5;
		}
	      }
	      else if ((k == 1 && user->bctype[4]!=7) || nvert[k-1][j][i] + nvert[k-1][j][i-1] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j][i-1] < 0.1) {
		  vol[TP] -= g13[k][j][i-1] * 0.5; vol[TW] -= g13[k][j][i-1] * 0.5;
		  vol[CP] += g13[k][j][i-1] * 0.5; vol[WP] += g13[k][j][i-1] * 0.5;
		}
	      }
	      else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j][i-1] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j][i-1] < 0.1) {
		  vol[TP] -= g13[k][j][i-1] * 0.5; vol[TW] -= g13[k][j][i-1] * 0.5;
		  vol[CP] += g13[k][j][i-1] * 0.5; vol[WP] += g13[k][j][i-1] * 0.5;
		}
	      }
	      else {
		vol[TP] -= g13[k][j][i-1] * 0.25; vol[TW] -= g13[k][j][i-1] * 0.25;
		vol[BP] += g13[k][j][i-1] * 0.25; vol[BW] += g13[k][j][i-1] * 0.25;
	      }
          }

          /************************************************************************
           * NORTH FACE CONTRIBUTION (between j and j+1)
           ************************************************************************/
          if (nvert[k][j+1][i] < IBM_FLUID_THRESHOLD && j != y_end) {
	      if ((i == mx-2 && user->bctype[0]!=7)|| nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
		  vol[CP] += g21[k][j][i] * 0.5; vol[NP] += g21[k][j][i] * 0.5;
		  vol[WP] -= g21[k][j][i] * 0.5; vol[NW] -= g21[k][j][i] * 0.5;
		}
	      }
	      else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k][j+1][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k][j+1][i-1] < 0.1) {
		  vol[CP] += g21[k][j][i] * 0.5; vol[NP] += g21[k][j][i] * 0.5;
		  vol[WP] -= g21[k][j][i] * 0.5; vol[NW] -= g21[k][j][i] * 0.5;
		}
	      }
	      else if ((i == 1 && user->bctype[0]!=7) || nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
		  vol[EP] += g21[k][j][i] * 0.5; vol[NE] += g21[k][j][i] * 0.5;
		  vol[CP] -= g21[k][j][i] * 0.5; vol[NP] -= g21[k][j][i] * 0.5;
		}
	      }
	      else if ((i == 1 || i==mx-2) && user->bctype[0]==7 &&  nvert[k][j][i-1] + nvert[k][j+1][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k][j+1][i+1] < 0.1) {
		  vol[EP] += g21[k][j][i] * 0.5; vol[NE] += g21[k][j][i] * 0.5;
		  vol[CP] -= g21[k][j][i] * 0.5; vol[NP] -= g21[k][j][i] * 0.5;
		}
	      }
	      else {
		vol[EP] += g21[k][j][i] * 0.25; vol[NE] += g21[k][j][i] * 0.25;
		vol[WP] -= g21[k][j][i] * 0.25; vol[NW] -= g21[k][j][i] * 0.25;
	      }

	      vol[CP] -= g22[k][j][i];
	      vol[NP] += g22[k][j][i];

	      if ((k == mz-2 && user->bctype[4]!=7)|| nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
		  vol[CP] += g23[k][j][i] * 0.5; vol[NP] += g23[k][j][i] * 0.5;
		  vol[BP] -= g23[k][j][i] * 0.5; vol[BN] -= g23[k][j][i] * 0.5;
		}
	      }
	      else if ((k == mz-2 || k==1 ) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j+1][i] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j+1][i] < 0.1) {
		  vol[CP] += g23[k][j][i] * 0.5; vol[NP] += g23[k][j][i] * 0.5;
		  vol[BP] -= g23[k][j][i] * 0.5; vol[BN] -= g23[k][j][i] * 0.5;
		}
	      }
	      else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
		  vol[TP] += g23[k][j][i] * 0.5; vol[TN] += g23[k][j][i] * 0.5;
		  vol[CP] -= g23[k][j][i] * 0.5; vol[NP] -= g23[k][j][i] * 0.5;
		}
	      }
	      else if ((k == 1 || k==mz-2 ) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j+1][i] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j+1][i] < 0.1) {
		  vol[TP] += g23[k][j][i] * 0.5; vol[TN] += g23[k][j][i] * 0.5;
		  vol[CP] -= g23[k][j][i] * 0.5; vol[NP] -= g23[k][j][i] * 0.5;
		}
	      }
	      else {
		vol[TP] += g23[k][j][i] * 0.25; vol[TN] += g23[k][j][i] * 0.25;
		vol[BP] -= g23[k][j][i] * 0.25; vol[BN] -= g23[k][j][i] * 0.25;
	      }
          }

          /************************************************************************
           * SOUTH FACE CONTRIBUTION (between j-1 and j)
           ************************************************************************/
          if (nvert[k][j-1][i] < IBM_FLUID_THRESHOLD && j != y_str) {
	      if ((i == mx-2 && user->bctype[0]!=7) || nvert[k][j][i+1] + nvert[k][j-1][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k][j-1][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
		  vol[CP] -= g21[k][j-1][i] * 0.5; vol[SP] -= g21[k][j-1][i] * 0.5;
		  vol[WP] += g21[k][j-1][i] * 0.5; vol[SW] += g21[k][j-1][i] * 0.5;
		}
	      }
	      else  if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k][j-1][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k][j-1][i-1] < 0.1) {
		  vol[CP] -= g21[k][j-1][i] * 0.5; vol[SP] -= g21[k][j-1][i] * 0.5;
		  vol[WP] += g21[k][j-1][i] * 0.5; vol[SW] += g21[k][j-1][i] * 0.5;
		}
	      }
	      else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k][j-1][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k][j-1][i+1] < 0.1) {
		  vol[EP] -= g21[k][j-1][i] * 0.5; vol[SE] -= g21[k][j-1][i] * 0.5;
		  vol[CP] += g21[k][j-1][i] * 0.5; vol[SP] += g21[k][j-1][i] * 0.5;
		}
	      }
	      else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k][j-1][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k][j-1][i+1] < 0.1) {
		  vol[EP] -= g21[k][j-1][i] * 0.5; vol[SE] -= g21[k][j-1][i] * 0.5;
		  vol[CP] += g21[k][j-1][i] * 0.5; vol[SP] += g21[k][j-1][i] * 0.5;
		}
	      }
	      else {
		vol[EP] -= g21[k][j-1][i] * 0.25; vol[SE] -= g21[k][j-1][i] * 0.25;
		vol[WP] += g21[k][j-1][i] * 0.25; vol[SW] += g21[k][j-1][i] * 0.25;
	      }

	      vol[CP] -= g22[k][j-1][i];
	      vol[SP] += g22[k][j-1][i];

	      if ((k == mz-2 && user->bctype[4]!=7)|| nvert[k+1][j][i] + nvert[k+1][j-1][i] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j-1][i] < 0.1 && (k!=1 || (k==1 && user->bctype[4]==7))) {
		  vol[CP] -= g23[k][j-1][i] * 0.5; vol[SP] -= g23[k][j-1][i] * 0.5;
		  vol[BP] += g23[k][j-1][i] * 0.5; vol[BS] += g23[k][j-1][i] * 0.5;
		}
	      }
	      else if ((k == mz-2 || k==1) && user->bctype[4]==7 && nvert[k+1][j][i] + nvert[k+1][j-1][i] > 0.1) {
		if (nvert[k-1][j][i] + nvert[k-1][j-1][i] < 0.1 ) {
		  vol[CP] -= g23[k][j-1][i] * 0.5; vol[SP] -= g23[k][j-1][i] * 0.5;
		  vol[BP] += g23[k][j-1][i] * 0.5; vol[BS] += g23[k][j-1][i] * 0.5;
		}
	      }
	      else if ((k == 1 && user->bctype[4]!=7)|| nvert[k-1][j][i] + nvert[k-1][j-1][i] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j-1][i] < 0.1) {
		  vol[TP] -= g23[k][j-1][i] * 0.5; vol[TS] -= g23[k][j-1][i] * 0.5;
		  vol[CP] += g23[k][j-1][i] * 0.5; vol[SP] += g23[k][j-1][i] * 0.5;
		}
	      }
	      else if ((k == 1 || k==mz-2) && user->bctype[4]==7 && nvert[k-1][j][i] + nvert[k-1][j-1][i] > 0.1) {
		if (nvert[k+1][j][i] + nvert[k+1][j-1][i] < 0.1) {
		  vol[TP] -= g23[k][j-1][i] * 0.5; vol[TS] -= g23[k][j-1][i] * 0.5;
		  vol[CP] += g23[k][j-1][i] * 0.5; vol[SP] += g23[k][j-1][i] * 0.5;
		}
	      }
	      else {
		vol[TP] -= g23[k][j-1][i] * 0.25; vol[TS] -= g23[k][j-1][i] * 0.25;
		vol[BP] += g23[k][j-1][i] * 0.25; vol[BS] += g23[k][j-1][i] * 0.25;
	      }
          }

          /************************************************************************
           * TOP FACE CONTRIBUTION (between k and k+1)
           ************************************************************************/
          if (nvert[k+1][j][i] < IBM_FLUID_THRESHOLD && k != z_end) {
	      if ((i == mx-2 && user->bctype[0]!=7)|| nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
		  vol[CP] += g31[k][j][i] * 0.5; vol[TP] += g31[k][j][i] * 0.5;
		  vol[WP] -= g31[k][j][i] * 0.5; vol[TW] -= g31[k][j][i] * 0.5;
		}
	      }
	      else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k+1][j][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k+1][j][i-1] < 0.1) {
		  vol[CP] += g31[k][j][i] * 0.5; vol[TP] += g31[k][j][i] * 0.5;
		  vol[WP] -= g31[k][j][i] * 0.5; vol[TW] -= g31[k][j][i] * 0.5;
		}
	      }
	      else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
		  vol[EP] += g31[k][j][i] * 0.5; vol[TE] += g31[k][j][i] * 0.5;
		  vol[CP] -= g31[k][j][i] * 0.5; vol[TP] -= g31[k][j][i] * 0.5;
		}
	      }
	      else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k+1][j][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k+1][j][i+1] < 0.1) {
		  vol[EP] += g31[k][j][i] * 0.5; vol[TE] += g31[k][j][i] * 0.5;
		  vol[CP] -= g31[k][j][i] * 0.5; vol[TP] -= g31[k][j][i] * 0.5;
		}
	      }
	      else {
		vol[EP] += g31[k][j][i] * 0.25; vol[TE] += g31[k][j][i] * 0.25;
		vol[WP] -= g31[k][j][i] * 0.25; vol[TW] -= g31[k][j][i] * 0.25;
	      }

	      if ((j == my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
		  vol[CP] += g32[k][j][i] * 0.5; vol[TP] += g32[k][j][i] * 0.5;
		  vol[SP] -= g32[k][j][i] * 0.5; vol[TS] -= g32[k][j][i] * 0.5;
		}
	      }
	      else if ((j == my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i] + nvert[k+1][j+1][i] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k+1][j-1][i] < 0.1) {
		  vol[CP] += g32[k][j][i] * 0.5; vol[TP] += g32[k][j][i] * 0.5;
		  vol[SP] -= g32[k][j][i] * 0.5; vol[TS] -= g32[k][j][i] * 0.5;
		}
	      }
	      else if ((j == 1 && user->bctype[2]!=7)|| nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
		  vol[NP] += g32[k][j][i] * 0.5; vol[TN] += g32[k][j][i] * 0.5;
		  vol[CP] -= g32[k][j][i] * 0.5; vol[TP] -= g32[k][j][i] * 0.5;
		}
	      }
	      else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k+1][j-1][i] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k+1][j+1][i] < 0.1) {
		  vol[NP] += g32[k][j][i] * 0.5; vol[TN] += g32[k][j][i] * 0.5;
		  vol[CP] -= g32[k][j][i] * 0.5; vol[TP] -= g32[k][j][i] * 0.5;
		}
	      }
	      else {
		vol[NP] += g32[k][j][i] * 0.25; vol[TN] += g32[k][j][i] * 0.25;
		vol[SP] -= g32[k][j][i] * 0.25; vol[TS] -= g32[k][j][i] * 0.25;
	      }

	      vol[CP] -= g33[k][j][i];
	      vol[TP] += g33[k][j][i];
          }

          /************************************************************************
           * BOTTOM FACE CONTRIBUTION (between k-1 and k)
           ************************************************************************/
          if (nvert[k-1][j][i] < IBM_FLUID_THRESHOLD && k != z_str) {
	      if ((i == mx-2 && user->bctype[0]!=7)|| nvert[k][j][i+1] + nvert[k-1][j][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k-1][j][i-1] < 0.1 && (i!=1 || (i==1 && user->bctype[0]==7))) {
		  vol[CP] -= g31[k-1][j][i] * 0.5; vol[BP] -= g31[k-1][j][i] * 0.5;
		  vol[WP] += g31[k-1][j][i] * 0.5; vol[BW] += g31[k-1][j][i] * 0.5;
		}
	      }
	      else if ((i == mx-2 || i==1) && user->bctype[0]==7 && nvert[k][j][i+1] + nvert[k-1][j][i+1] > 0.1) {
		if (nvert[k][j][i-1] + nvert[k-1][j][i-1] < 0.1) {
		  vol[CP] -= g31[k-1][j][i] * 0.5; vol[BP] -= g31[k-1][j][i] * 0.5;
		  vol[WP] += g31[k-1][j][i] * 0.5; vol[BW] += g31[k-1][j][i] * 0.5;
		}
	      }
	      else if ((i == 1 && user->bctype[0]!=7)|| nvert[k][j][i-1] + nvert[k-1][j][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k-1][j][i+1] < 0.1) {
		  vol[EP] -= g31[k-1][j][i] * 0.5; vol[BE] -= g31[k-1][j][i] * 0.5;
		  vol[CP] += g31[k-1][j][i] * 0.5; vol[BP] += g31[k-1][j][i] * 0.5;
		}
	      }
	      else if ((i == 1 || i==mx-2) && user->bctype[0]==7 && nvert[k][j][i-1] + nvert[k-1][j][i-1] > 0.1) {
		if (nvert[k][j][i+1] + nvert[k-1][j][i+1] < 0.1) {
		  vol[EP] -= g31[k-1][j][i] * 0.5; vol[BE] -= g31[k-1][j][i] * 0.5;
		  vol[CP] += g31[k-1][j][i] * 0.5; vol[BP] += g31[k-1][j][i] * 0.5;
		}
	      }
	      else {
		vol[EP] -= g31[k-1][j][i] * 0.25; vol[BE] -= g31[k-1][j][i] * 0.25;
		vol[WP] += g31[k-1][j][i] * 0.25; vol[BW] += g31[k-1][j][i] * 0.25;
	      }

	      if ((j == my-2 && user->bctype[2]!=7)|| nvert[k][j+1][i] + nvert[k-1][j+1][i] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k-1][j-1][i] < 0.1 && (j!=1 || (j==1 && user->bctype[2]==7))) {
		  vol[CP] -= g32[k-1][j][i] * 0.5; vol[BP] -= g32[k-1][j][i] * 0.5;
		  vol[SP] += g32[k-1][j][i] * 0.5; vol[BS] += g32[k-1][j][i] * 0.5;
		}
	      }
	      else if ((j == my-2 || j==1) && user->bctype[2]==7 && nvert[k][j+1][i] + nvert[k-1][j+1][i] > 0.1) {
		if (nvert[k][j-1][i] + nvert[k-1][j-1][i] < 0.1) {
		  vol[CP] -= g32[k-1][j][i] * 0.5; vol[BP] -= g32[k-1][j][i] * 0.5;
		  vol[SP] += g32[k-1][j][i] * 0.5; vol[BS] += g32[k-1][j][i] * 0.5;
		}
	      }
	      else if ((j == 1 && user->bctype[2]!=7)|| nvert[k][j-1][i] + nvert[k-1][j-1][i] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k-1][j+1][i] < 0.1) {
		  vol[NP] -= g32[k-1][j][i] * 0.5; vol[BN] -= g32[k-1][j][i] * 0.5;
		  vol[CP] += g32[k-1][j][i] * 0.5; vol[BP] += g32[k-1][j][i] * 0.5;
		}
	      }
	      else if ((j == 1 || j==my-2) && user->bctype[2]==7 && nvert[k][j-1][i] + nvert[k-1][j-1][i] > 0.1) {
		if (nvert[k][j+1][i] + nvert[k-1][j+1][i] < 0.1) {
		  vol[NP] -= g32[k-1][j][i] * 0.5; vol[BN] -= g32[k-1][j][i] * 0.5;
		  vol[CP] += g32[k-1][j][i] * 0.5; vol[BP] += g32[k-1][j][i] * 0.5;
		}
	      }
	      else {
		vol[NP] -= g32[k-1][j][i] * 0.25; vol[BN] -= g32[k-1][j][i] * 0.25;
		vol[SP] += g32[k-1][j][i] * 0.25; vol[BS] += g32[k-1][j][i] * 0.25;
	      }

	      vol[CP] -= g33[k-1][j][i];
	      vol[BP] += g33[k-1][j][i];
          }

          // --- Final scaling and insertion into the matrix ---

          // Scale all stencil coefficients by the negative cell volume (-aj).
          for (PetscInt m = 0; m < 19; m++) {
            vol[m] *= -aj[k][j][i];
          }

          // Set the global column indices for the 19 stencil points, handling periodic BCs.
          idx[CP] = Gidx(i, j, k, user);
          if (user->bctype[0]==7 && i==mx-2) idx[EP] = Gidx(1, j, k, user); else idx[EP] = Gidx(i+1, j, k, user);
          if (user->bctype[0]==7 && i==1) idx[WP] = Gidx(mx-2, j, k, user); else idx[WP] = Gidx(i-1, j, k, user);
          if (user->bctype[2]==7 && j==my-2) idx[NP] = Gidx(i, 1, k, user); else idx[NP] = Gidx(i, j+1, k, user);
          if (user->bctype[2]==7 && j==1) idx[SP] = Gidx(i, my-2, k, user); else idx[SP] = Gidx(i, j-1, k, user);
          if (user->bctype[4]==7 && k==mz-2) idx[TP] = Gidx(i, j, 1, user); else idx[TP] = Gidx(i, j, k+1, user);
          if (user->bctype[4]==7 && k==1) idx[BP] = Gidx(i, j, mz-2, user); else idx[BP] = Gidx(i, j, k-1, user);
          if (user->bctype[0]==7 && user->bctype[2]==7 && i==mx-2 && j==my-2) idx[NE] = Gidx(1, 1, k, user); else if (user->bctype[0]==7 && i==mx-2) idx[NE] = Gidx(1, j+1, k, user); else if (user->bctype[2]==7 && j==my-2) idx[NE] = Gidx(i+1, 1, k, user); else idx[NE] = Gidx(i+1, j+1, k, user);
          if (user->bctype[0]==7 && user->bctype[2]==7 && i==mx-2 && j==1) idx[SE] = Gidx(1, my-2, k, user); else if (user->bctype[0]==7 && i==mx-2) idx[SE] = Gidx(1, j-1, k, user); else if (user->bctype[2]==7 && j==1) idx[SE] = Gidx(i+1, my-2, k, user); else idx[SE] = Gidx(i+1, j-1, k, user);
          if (user->bctype[0]==7 && user->bctype[2]==7 && i==1 && j==my-2) idx[NW] = Gidx(mx-2, 1, k, user); else if (user->bctype[0]==7 && i==1) idx[NW] = Gidx(mx-2, j+1, k, user); else if (user->bctype[2]==7 && j==my-2) idx[NW] = Gidx(i-1, 1, k, user); else idx[NW] = Gidx(i-1, j+1, k, user);
          if (user->bctype[0]==7 && user->bctype[2]==7 && i==1 && j==1) idx[SW] = Gidx(mx-2, my-2, k, user); else if (user->bctype[0]==7 && i==1) idx[SW] = Gidx(mx-2, j-1, k, user); else if (user->bctype[2]==7 && j==1) idx[SW] = Gidx(i-1, my-2, k, user); else idx[SW] = Gidx(i-1, j-1, k, user);
          if (user->bctype[2]==7 && user->bctype[4]==7 && j==my-2 && k==mz-2) idx[TN] = Gidx(i, 1, 1, user); else if (user->bctype[2]==7 && j==my-2) idx[TN] = Gidx(i, 1, k+1, user); else if (user->bctype[4]==7 && k==mz-2) idx[TN] = Gidx(i, j+1, 1, user); else idx[TN] = Gidx(i, j+1, k+1, user);
          if (user->bctype[2]==7 && user->bctype[4]==7 && j==my-2 && k==1) idx[BN] = Gidx(i, 1, mz-2, user); else if(user->bctype[2]==7 && j==my-2) idx[BN] = Gidx(i, 1, k-1, user); else if (user->bctype[4]==7 && k==1) idx[BN] = Gidx(i, j+1, mz-2, user); else idx[BN] = Gidx(i, j+1, k-1, user);
          if (user->bctype[2]==7 && user->bctype[4]==7 && j==1 && k==mz-2) idx[TS] = Gidx(i, my-2, 1, user); else if (user->bctype[2]==7 && j==1) idx[TS] = Gidx(i, my-2, k+1, user); else if (user->bctype[4]==7 && k==mz-2) idx[TS] = Gidx(i, j-1, 1, user); else idx[TS] = Gidx(i, j-1, k+1, user);
          if (user->bctype[2]==7 && user->bctype[4]==7 && j==1 && k==1) idx[BS] = Gidx(i, my-2, mz-2, user); else if (user->bctype[2]==7 && j==1) idx[BS] = Gidx(i, my-2, k-1, user); else if (user->bctype[4]==7 && k==1) idx[BS] = Gidx(i, j-1, mz-2, user); else idx[BS] = Gidx(i, j-1, k-1, user);
          if (user->bctype[0]==7 && user->bctype[4]==7 && i==mx-2 && k==mz-2) idx[TE] = Gidx(1, j, 1, user); else if(user->bctype[0]==7 && i==mx-2) idx[TE] = Gidx(1, j, k+1, user); else if(user->bctype[4]==7 && k==mz-2) idx[TE] = Gidx(i+1, j, 1, user); else idx[TE] = Gidx(i+1, j, k+1, user);
          if (user->bctype[0]==7 && user->bctype[4]==7 && i==mx-2 && k==1) idx[BE] = Gidx(1, j, mz-2, user); else if(user->bctype[0]==7 && i==mx-2) idx[BE] = Gidx(1, j, k-1, user); else if(user->bctype[4]==7 && k==1) idx[BE] = Gidx(i+1, j, mz-2, user); else idx[BE] = Gidx(i+1, j, k-1, user);
          if (user->bctype[0]==7 && user->bctype[4]==7 && i==1 && k==mz-2) idx[TW] = Gidx(mx-2, j, 1, user); else if(user->bctype[0]==7 && i==1) idx[TW] = Gidx(mx-2, j, k+1, user); else if (user->bctype[4]==7 && k==mz-2) idx[TW] = Gidx(i-1, j, 1, user); else idx[TW] = Gidx(i-1, j, k+1, user);
          if (user->bctype[0]==7 && user->bctype[4]==7 && i==1 && k==1) idx[BW] = Gidx(mx-2, j, mz-2, user); else if (user->bctype[0]==7 && i==1) idx[BW] = Gidx(mx-2, j, k-1, user); else if (user->bctype[4]==7 && k==1) idx[BW] = Gidx(i-1, j, mz-2, user); else idx[BW] = Gidx(i-1, j, k-1, user);

          // Insert the computed row into the matrix A.
          MatSetValues(user->A, 1, &row, 19, idx, vol, INSERT_VALUES);
        }
      }
    }
  }

  //================================================================================
  // Section 4: Finalize Matrix and Cleanup
  //================================================================================
  
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Finalizing matrix assembly.\n");
  MatAssemblyBegin(user->A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(user->A, MAT_FINAL_ASSEMBLY);

  PetscReal max_A;

  ierr = MatNorm(user->A,NORM_INFINITY,&max_A);CHKERRQ(ierr);

  LOG_ALLOW(GLOBAL,LOG_DEBUG," Max value in A matrix for level %d =  %le.\n",user->thislevel,max_A);

  // if (get_log_level() >= LOG_DEBUG) {
  //  ierr = MatView(user->A,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  // }
  
  // --- Restore access to all PETSc vectors and destroy temporary ones ---
  DMDAVecRestoreArray(da, G11, &g11); DMDAVecRestoreArray(da, G12, &g12); DMDAVecRestoreArray(da, G13, &g13);
  DMDAVecRestoreArray(da, G21, &g21); DMDAVecRestoreArray(da, G22, &g22); DMDAVecRestoreArray(da, G23, &g23);
  DMDAVecRestoreArray(da, G31, &g31); DMDAVecRestoreArray(da, G32, &g32); DMDAVecRestoreArray(da, G33, &g33);
  
  VecDestroy(&G11); VecDestroy(&G12); VecDestroy(&G13);
  VecDestroy(&G21); VecDestroy(&G22); VecDestroy(&G23);
  VecDestroy(&G31); VecDestroy(&G32); VecDestroy(&G33);

  DMDAVecRestoreArray(fda, user->lCsi, &csi); DMDAVecRestoreArray(fda, user->lEta, &eta); DMDAVecRestoreArray(fda, user->lZet, &zet);
  DMDAVecRestoreArray(fda, user->lICsi, &icsi); DMDAVecRestoreArray(fda, user->lIEta, &ieta); DMDAVecRestoreArray(fda, user->lIZet, &izet);
  DMDAVecRestoreArray(fda, user->lJCsi, &jcsi); DMDAVecRestoreArray(fda, user->lJEta, &jeta); DMDAVecRestoreArray(fda, user->lJZet, &jzet);
  DMDAVecRestoreArray(fda, user->lKCsi, &kcsi); DMDAVecRestoreArray(fda, user->lKEta, &keta); DMDAVecRestoreArray(fda, user->lKZet, &kzet);
  DMDAVecRestoreArray(da, user->lAj, &aj); DMDAVecRestoreArray(da, user->lIAj, &iaj); DMDAVecRestoreArray(da, user->lJAj, &jaj); DMDAVecRestoreArray(da, user->lKAj, &kaj);
  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Exiting PoissonLHSNew.\n");
  PetscFunctionReturn(0);
}

PetscErrorCode PoissonRHS(UserCtx *user, Vec B)
{
  DMDALocalInfo info = user->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;
 
  PetscInt      i, j, k;
  PetscReal	***nvert, ***aj, ***rb, dt = user->simCtx->dt;
  struct Components{
    PetscReal x;
    PetscReal y;
    PetscReal z;
  } *** ucont;

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Entering PoissonRHS to compute pressure equation RHS.\n");

  DMDAVecGetArray(user->da, B, &rb);
  DMDAVecGetArray(user->fda, user->lUcont, &ucont);
  DMDAVecGetArray(user->da, user->lNvert, &nvert);
  DMDAVecGetArray(user->da, user->lAj, &aj);


   LOG_ALLOW(GLOBAL, LOG_DEBUG, "Computing RHS values for each cell.\n");
  
  for (k=zs; k<ze; k++) { 
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {

	if (i==0 || i==mx-1 || j==0 || j==my-1 ||  k==0 || k==mz-1) {
	  rb[k][j][i] = 0.;
	}
	else if (nvert[k][j][i] > 0.1) {
	  rb[k][j][i] = 0;
	}
	else {
	  rb[k][j][i] = -(ucont[k][j][i].x - ucont[k][j][i-1].x +
			  ucont[k][j][i].y - ucont[k][j-1][i].y +
			  ucont[k][j][i].z - ucont[k-1][j][i].z) / dt
 	    * aj[k][j][i] / user->simCtx->st * COEF_TIME_ACCURACY;
	 
	}
      }
    }
  }


  // --- Check the solvability condition for the Poisson equation ---
  // The global sum of the RHS (proportional to the total divergence) must be zero.
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Verifying solvability condition (sum of RHS terms).\n");
  PetscReal lsum=0., sum=0.;

  for (k=zs; k<ze; k++) {
    for (j=ys; j<ye; j++) {
      for (i=xs; i<xe; i++) {
	
	lsum += rb[k][j][i] / aj[k][j][i]* dt/COEF_TIME_ACCURACY;

      }
    }
  }
  
  MPI_Allreduce(&lsum,&sum,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);

  LOG_ALLOW(GLOBAL, LOG_INFO, "Global Sum of RHS (Divergence Check): %le\n", sum);

  user->simCtx->summationRHS = sum;
	
  DMDAVecRestoreArray(user->fda, user->lUcont, &ucont);
  DMDAVecRestoreArray(user->da, user->lNvert, &nvert);
  DMDAVecRestoreArray(user->da, user->lAj, &aj);
  DMDAVecRestoreArray(user->da, B, &rb);
 
  return 0;
}

PetscErrorCode VolumeFlux_rev(UserCtx *user, PetscReal *ibm_Flux, 
			      PetscReal *ibm_Area, PetscInt flg)
{

  // --- CONTEXT ACQUISITION BLOCK ---
  // Get the master simulation context from the UserCtx.
  SimCtx *simCtx = user->simCtx;

  // Create a local variable to mirror the legacy global for minimal code changes.
  const PetscInt NumberOfBodies = simCtx->NumberOfBodies;
  // --- END CONTEXT ACQUISITION BLOCK ---
  
  DM	da = user->da, fda = user->fda;

  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt i, j, k;
  PetscInt	lxs, lys, lzs, lxe, lye, lze;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  PetscReal ***nvert, ibmval=1.5;
  Cmpnts ***ucor, ***csi, ***eta, ***zet;
  DMDAVecGetArray(fda, user->Ucont, &ucor);
  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  PetscReal libm_Flux, libm_area;
  libm_Flux = 0;
  libm_area = 0;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > ibmval-0.4 && nvert[k][j][i+1] < ibmval && i < mx-2) {
	    libm_Flux += ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] > ibmval-0.4 && nvert[k][j+1][i] < ibmval && j < my-2) {
	    libm_Flux += ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] > ibmval-0.4 && nvert[k+1][j][i] < ibmval && k < mz-2) {
	    libm_Flux += ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

	if (nvert[k][j][i] > ibmval-0.4 && nvert[k][j][i] < ibmval) {
	  if (nvert[k][j][i+1] < 0.1 && i < mx-2) {
	    libm_Flux -= ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < my-2) {
	    libm_Flux -= ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < mz-2) {
	    libm_Flux -= ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

      }
    }
  }

  MPI_Allreduce(&libm_Flux, ibm_Flux,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&libm_area, ibm_Area,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);

 /*  PetscGlobalSum(&libm_Flux, ibm_Flux, PETSC_COMM_WORLD); */
/*   PetscGlobalSum(&libm_area, ibm_Area, PETSC_COMM_WORLD); */
  PetscPrintf(PETSC_COMM_WORLD, "IBMFlux REV %le %le\n", *ibm_Flux, *ibm_Area);

  PetscReal correction;

  if (*ibm_Area > 1.e-15) {
    if (flg)
      correction = (*ibm_Flux + user->FluxIntpSum) / *ibm_Area;
    else
      correction = *ibm_Flux / *ibm_Area;
  }
  else {
    correction = 0;
  }

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > ibmval-0.4 && nvert[k][j][i+1] < ibmval && i < mx-2) {
	    ucor[k][j][i].x -= sqrt(csi[k][j][i].x * csi[k][j][i].x +
				    csi[k][j][i].y * csi[k][j][i].y +
				    csi[k][j][i].z * csi[k][j][i].z) *
				    correction;
			  
	  }
	  if (nvert[k][j+1][i] > ibmval-0.4 && nvert[k][j+1][i] < ibmval && j < my-2) {
	    ucor[k][j][i].y -= sqrt(eta[k][j][i].x * eta[k][j][i].x + 
				    eta[k][j][i].y * eta[k][j][i].y +
				    eta[k][j][i].z * eta[k][j][i].z) *
				    correction;
	  }
	  if (nvert[k+1][j][i] > ibmval-0.4 && nvert[k+1][j][i] < ibmval && k < mz-2) {
	    ucor[k][j][i].z -= sqrt(zet[k][j][i].x * zet[k][j][i].x +
				    zet[k][j][i].y * zet[k][j][i].y +
				    zet[k][j][i].z * zet[k][j][i].z) *
				    correction;
	  }
	}

	if (nvert[k][j][i] > ibmval-0.4 && nvert[k][j][i] < ibmval) {
	  if (nvert[k][j][i+1] < 0.1 && i < mx-2) {
	    ucor[k][j][i].x += sqrt(csi[k][j][i].x * csi[k][j][i].x +
				    csi[k][j][i].y * csi[k][j][i].y +
				    csi[k][j][i].z * csi[k][j][i].z) *
				    correction;
			  
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < my-2) {
	    ucor[k][j][i].y += sqrt(eta[k][j][i].x * eta[k][j][i].x +
				    eta[k][j][i].y * eta[k][j][i].y +
				    eta[k][j][i].z * eta[k][j][i].z) *
				    correction;
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < mz-2) {
	    ucor[k][j][i].z += sqrt(zet[k][j][i].x * zet[k][j][i].x +
				    zet[k][j][i].y * zet[k][j][i].y +
				    zet[k][j][i].z * zet[k][j][i].z) *
				    correction;
	  }
	}

      }
    }
  }
  


  libm_Flux = 0;
  libm_area = 0;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > ibmval-0.4 && nvert[k][j][i+1] < ibmval && i < mx-2) {
	    libm_Flux += ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] > ibmval-0.4 && nvert[k][j+1][i] < ibmval && j < my-2) {
	    libm_Flux += ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] > ibmval-0.4 && nvert[k+1][j][i] < ibmval && k < mz-2) {
	    libm_Flux += ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

	if (nvert[k][j][i] > ibmval-0.4 && nvert[k][j][i] < ibmval) {
	  if (nvert[k][j][i+1] < 0.1 && i < mx-2) {
	    libm_Flux -= ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < my-2) {
	    libm_Flux -= ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < mz-2) {
	    libm_Flux -= ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

      }
    }
  }

  MPI_Allreduce(&libm_Flux, ibm_Flux,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&libm_area, ibm_Area,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);

 /*  PetscGlobalSum(&libm_Flux, ibm_Flux, PETSC_COMM_WORLD); */
/*   PetscGlobalSum(&libm_area, ibm_Area, PETSC_COMM_WORLD); */
  PetscPrintf(PETSC_COMM_WORLD, "IBMFlux22 REV %le %le\n", *ibm_Flux, *ibm_Area);

  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  DMDAVecRestoreArray(fda, user->Ucont, &ucor);

  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  return 0;
}


PetscErrorCode VolumeFlux(UserCtx *user, PetscReal *ibm_Flux, PetscReal *ibm_Area, PetscInt flg)
{
  // --- CONTEXT ACQUISITION BLOCK ---
  // Get the master simulation context from the UserCtx.
  SimCtx *simCtx = user->simCtx;

  // Create local variables to mirror the legacy globals for minimal code changes.
  const PetscInt NumberOfBodies = simCtx->NumberOfBodies;
  const PetscInt channelz = simCtx->channelz;
  const PetscInt fish = simCtx->fish;
  // --- END CONTEXT ACQUISITION BLOCK ---

  DM	da = user->da, fda = user->fda;

  DMDALocalInfo	info = user->info;

  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt i, j, k,ibi;
  PetscInt	lxs, lys, lzs, lxe, lye, lze;

  PetscInt	gxs, gxe, gys, gye, gzs, gze;
  
  gxs = info.gxs; gxe = gxs + info.gxm;
  gys = info.gys; gye = gys + info.gym;
  gzs = info.gzs; gze = gzs + info.gzm;

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;
  
  PetscReal epsilon=1.e-8;
  PetscReal ***nvert, ibmval=1.9999;
  
  struct Components {
    PetscReal x;
    PetscReal y;
    PetscReal z;
  }***ucor, ***lucor,***csi, ***eta, ***zet;
 

  PetscInt xend=mx-2 ,yend=my-2,zend=mz-2;

  if (user->bctype[0]==7) xend=mx-1;
  if (user->bctype[2]==7) yend=my-1;
  if (user->bctype[4]==7) zend=mz-1;

  DMDAVecGetArray(fda, user->Ucont, &ucor);
  DMDAVecGetArray(fda, user->lCsi, &csi);
  DMDAVecGetArray(fda, user->lEta, &eta);
  DMDAVecGetArray(fda, user->lZet, &zet);
  DMDAVecGetArray(da, user->lNvert, &nvert);

  PetscReal libm_Flux, libm_area, libm_Flux_abs=0., ibm_Flux_abs;
  libm_Flux = 0;
  libm_area = 0;

  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Entering VolumeFlux to enforce no-penetration condition.\n");

  //Mohsen March 2017
  PetscReal *lIB_Flux, *lIB_area,*IB_Flux,*IB_Area;
  if (NumberOfBodies > 1) { 
  
    lIB_Flux=(PetscReal *)calloc(NumberOfBodies,sizeof(PetscReal));
    lIB_area=(PetscReal *)calloc(NumberOfBodies,sizeof(PetscReal));
    IB_Flux=(PetscReal *)calloc(NumberOfBodies,sizeof(PetscReal));
    IB_Area=(PetscReal *)calloc(NumberOfBodies,sizeof(PetscReal));
  

    for (ibi=0; ibi<NumberOfBodies; ibi++) {
      lIB_Flux[ibi]=0.0;
      lIB_area[ibi]=0.0;
      IB_Flux[ibi]=0.0;
      IB_Area[ibi]=0.0;
    }
  }


  //================================================================================
  // PASS 1: Calculate Uncorrected Flux and Area
  // This pass measures the total fluid "leakage" across the immersed boundary.
  //================================================================================
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Pass 1: Measuring uncorrected flux and area.\n");
  
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > 0.1 && nvert[k][j][i+1] < ibmval && i < xend) {
	    
	    if (fabs(ucor[k][j][i].x)>epsilon) {
	      libm_Flux += ucor[k][j][i].x;
	      if (flg==3)
		libm_Flux_abs += fabs(ucor[k][j][i].x)/sqrt(csi[k][j][i].x * csi[k][j][i].x +
							    csi[k][j][i].y * csi[k][j][i].y +
							    csi[k][j][i].z * csi[k][j][i].z);
	      else
		libm_Flux_abs += fabs(ucor[k][j][i].x);
	      
	      libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
				csi[k][j][i].y * csi[k][j][i].y +
				csi[k][j][i].z * csi[k][j][i].z);
	      
	      if (NumberOfBodies > 1) {
		
		ibi=(int)((nvert[k][j][i+1]-1.0)*1001);
		lIB_Flux[ibi] += ucor[k][j][i].x;
		lIB_area[ibi] += sqrt(csi[k][j][i].x * csi[k][j][i].x +
				      csi[k][j][i].y * csi[k][j][i].y +
				      csi[k][j][i].z * csi[k][j][i].z);
	      }
	    } else
	      ucor[k][j][i].x=0.;
	    
	  }
	  if (nvert[k][j+1][i] > 0.1 && nvert[k][j+1][i] < ibmval && j < yend) {
	    
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	      libm_Flux += ucor[k][j][i].y;
	      if (flg==3)
		libm_Flux_abs += fabs(ucor[k][j][i].y)/sqrt(eta[k][j][i].x * eta[k][j][i].x +
							    eta[k][j][i].y * eta[k][j][i].y +
							    eta[k][j][i].z * eta[k][j][i].z);
	      else
		libm_Flux_abs += fabs(ucor[k][j][i].y);
	      libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
				eta[k][j][i].y * eta[k][j][i].y +
				eta[k][j][i].z * eta[k][j][i].z);
	      if (NumberOfBodies > 1) {
		
		ibi=(int)((nvert[k][j+1][i]-1.0)*1001);	

		lIB_Flux[ibi] += ucor[k][j][i].y;
		lIB_area[ibi] += sqrt(eta[k][j][i].x * eta[k][j][i].x +
				      eta[k][j][i].y * eta[k][j][i].y +
				      eta[k][j][i].z * eta[k][j][i].z);
	      }
	    } else
	      ucor[k][j][i].y=0.;
	  }
	  if (nvert[k+1][j][i] > 0.1 && nvert[k+1][j][i] < ibmval && k < zend) {
	    
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	      libm_Flux += ucor[k][j][i].z;
	      if (flg==3)
		libm_Flux_abs += fabs(ucor[k][j][i].z)/sqrt(zet[k][j][i].x * zet[k][j][i].x +
							    zet[k][j][i].y * zet[k][j][i].y +
							    zet[k][j][i].z * zet[k][j][i].z);
	      else
		libm_Flux_abs += fabs(ucor[k][j][i].z);
	      libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
				zet[k][j][i].y * zet[k][j][i].y +
				zet[k][j][i].z * zet[k][j][i].z);
	      
	      if (NumberOfBodies > 1) {
		
		ibi=(int)((nvert[k+1][j][i]-1.0)*1001);
		lIB_Flux[ibi] += ucor[k][j][i].z;
		lIB_area[ibi] += sqrt(zet[k][j][i].x * zet[k][j][i].x +
				      zet[k][j][i].y * zet[k][j][i].y +
				      zet[k][j][i].z * zet[k][j][i].z);
	      }
	    }else
	      ucor[k][j][i].z=0.;
	  }
	}

	if (nvert[k][j][i] > 0.1 && nvert[k][j][i] < ibmval) {
	  
	  if (nvert[k][j][i+1] < 0.1 && i < xend) {
	    if (fabs(ucor[k][j][i].x)>epsilon) {
	      libm_Flux -= ucor[k][j][i].x;
	      if (flg==3)
		libm_Flux_abs += fabs(ucor[k][j][i].x)/sqrt(csi[k][j][i].x * csi[k][j][i].x +
							    csi[k][j][i].y * csi[k][j][i].y +
							    csi[k][j][i].z * csi[k][j][i].z);
	      else
		libm_Flux_abs += fabs(ucor[k][j][i].x);
	      libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
				csi[k][j][i].y * csi[k][j][i].y +
				csi[k][j][i].z * csi[k][j][i].z);
	      if (NumberOfBodies > 1) {
		
		ibi=(int)((nvert[k][j][i]-1.0)*1001);
		lIB_Flux[ibi] -= ucor[k][j][i].x;
		lIB_area[ibi] += sqrt(csi[k][j][i].x * csi[k][j][i].x +
				      csi[k][j][i].y * csi[k][j][i].y +
				      csi[k][j][i].z * csi[k][j][i].z);
	      }
			  
	    }else
	      ucor[k][j][i].x=0.;
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < yend) {
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	      libm_Flux -= ucor[k][j][i].y;
	      if (flg==3)
		libm_Flux_abs += fabs(ucor[k][j][i].y)/ sqrt(eta[k][j][i].x * eta[k][j][i].x +
							     eta[k][j][i].y * eta[k][j][i].y +
							     eta[k][j][i].z * eta[k][j][i].z);
	      else
		libm_Flux_abs += fabs(ucor[k][j][i].y);
	      libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
				eta[k][j][i].y * eta[k][j][i].y +
				eta[k][j][i].z * eta[k][j][i].z);
	      if (NumberOfBodies > 1) {

		ibi=(int)((nvert[k][j][i]-1.0)*1001);
		lIB_Flux[ibi] -= ucor[k][j][i].y;
		lIB_area[ibi] += sqrt(eta[k][j][i].x * eta[k][j][i].x +
				      eta[k][j][i].y * eta[k][j][i].y +
				      eta[k][j][i].z * eta[k][j][i].z);
	      }
	    }else
	      ucor[k][j][i].y=0.;
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < zend) {
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	      libm_Flux -= ucor[k][j][i].z;
	      if (flg==3)
		libm_Flux_abs += fabs(ucor[k][j][i].z)/sqrt(zet[k][j][i].x * zet[k][j][i].x +
							    zet[k][j][i].y * zet[k][j][i].y +
							    zet[k][j][i].z * zet[k][j][i].z);
	      else
		libm_Flux_abs += fabs(ucor[k][j][i].z);
	      libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
				zet[k][j][i].y * zet[k][j][i].y +
				zet[k][j][i].z * zet[k][j][i].z);
	      if (NumberOfBodies > 1) {
		
		ibi=(int)((nvert[k][j][i]-1.0)*1001);
		lIB_Flux[ibi] -= ucor[k][j][i].z;
		lIB_area[ibi] += sqrt(zet[k][j][i].x * zet[k][j][i].x +
				      zet[k][j][i].y * zet[k][j][i].y +
				      zet[k][j][i].z * zet[k][j][i].z);
	      }
	    }else
	      ucor[k][j][i].z=0.;
	  }
	}
	
      }
    }
  }
  
  MPI_Allreduce(&libm_Flux, ibm_Flux,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&libm_Flux_abs, &ibm_Flux_abs,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&libm_area, ibm_Area,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);

  if (NumberOfBodies > 1) { 
    MPI_Allreduce(lIB_Flux,IB_Flux,NumberOfBodies,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    MPI_Allreduce(lIB_area,IB_Area,NumberOfBodies,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
  }

  PetscReal correction;

  PetscReal *Correction;
  if (NumberOfBodies > 1) {
      Correction=(PetscReal *)calloc(NumberOfBodies,sizeof(PetscReal));
      for (ibi=0; ibi<NumberOfBodies; ibi++) Correction[ibi]=0.0;
  }

  if (*ibm_Area > 1.e-15) {
    if (flg>1) 
      correction = (*ibm_Flux + user->FluxIntpSum)/ ibm_Flux_abs;
    else if (flg)
      correction = (*ibm_Flux + user->FluxIntpSum) / *ibm_Area;
    else
      correction = *ibm_Flux / *ibm_Area;
    if (NumberOfBodies > 1) 
      for (ibi=0; ibi<NumberOfBodies; ibi++) if (IB_Area[ibi]>1.e-15) Correction[ibi] = IB_Flux[ibi] / IB_Area[ibi];
  }
  else {
    correction = 0;
  }
  // --- Log the uncorrected results and calculated correction ---
  LOG_ALLOW(GLOBAL, LOG_INFO, "IBM Uncorrected Flux: %g, Area: %g, Correction: %g\n", *ibm_Flux, *ibm_Area, correction);
  if  (NumberOfBodies>1){
    for (ibi=0; ibi<NumberOfBodies; ibi++)   LOG_ALLOW(GLOBAL, LOG_INFO, "  [Body %d] Uncorrected Flux: %g, Area: %g, Correction: %g\n", ibi, IB_Flux[ibi], IB_Area[ibi], Correction[ibi]);
  }

  //================================================================================
  // PASS 2: Apply Correction to Velocity Field
  // This pass modifies the velocity at the boundary to enforce zero net flux.
  //================================================================================
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Pass 2: Applying velocity corrections at the boundary.\n");
  
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > 0.1 && nvert[k][j][i+1] <ibmval && i < xend) {
	    if (fabs(ucor[k][j][i].x)>epsilon){
	      if (flg==3) 
		ucor[k][j][i].x -=correction*fabs(ucor[k][j][i].x)/
		  sqrt(csi[k][j][i].x * csi[k][j][i].x +
		       csi[k][j][i].y * csi[k][j][i].y +
		       csi[k][j][i].z * csi[k][j][i].z);
	      else if (flg==2) 
		ucor[k][j][i].x -=correction*fabs(ucor[k][j][i].x);
	      else if (NumberOfBodies > 1) {
		ibi=(int)((nvert[k][j][i+1]-1.0)*1001);
		ucor[k][j][i].x -= sqrt(csi[k][j][i].x * csi[k][j][i].x +
					csi[k][j][i].y * csi[k][j][i].y +
					csi[k][j][i].z * csi[k][j][i].z) *
		  Correction[ibi];
	      }
	      else
		ucor[k][j][i].x -= sqrt(csi[k][j][i].x * csi[k][j][i].x +
					csi[k][j][i].y * csi[k][j][i].y +
					csi[k][j][i].z * csi[k][j][i].z) *
		  correction;
	    }
	  }
	  if (nvert[k][j+1][i] > 0.1 && nvert[k][j+1][i] < ibmval && j < yend) {
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	      if (flg==3) 
		ucor[k][j][i].y -=correction*fabs(ucor[k][j][i].y)/
		  sqrt(eta[k][j][i].x * eta[k][j][i].x + 
		       eta[k][j][i].y * eta[k][j][i].y +
		       eta[k][j][i].z * eta[k][j][i].z);
	      else if (flg==2) 
		ucor[k][j][i].y -=correction*fabs(ucor[k][j][i].y);
	      else if (NumberOfBodies > 1) {
		ibi=(int)((nvert[k][j+1][i]-1.0)*1001);
		ucor[k][j][i].y -= sqrt(eta[k][j][i].x * eta[k][j][i].x +
					eta[k][j][i].y * eta[k][j][i].y +
					eta[k][j][i].z * eta[k][j][i].z) *
		  Correction[ibi];
	      }
	      else
		ucor[k][j][i].y -= sqrt(eta[k][j][i].x * eta[k][j][i].x + 
					eta[k][j][i].y * eta[k][j][i].y +
					eta[k][j][i].z * eta[k][j][i].z) *
		  correction;
	    }
	  }
	  if (nvert[k+1][j][i] > 0.1 && nvert[k+1][j][i] < ibmval && k < zend) {
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	      if (flg==3) 
		ucor[k][j][i].z -= correction*fabs(ucor[k][j][i].z)/
		  sqrt(zet[k][j][i].x * zet[k][j][i].x +
		       zet[k][j][i].y * zet[k][j][i].y +
		       zet[k][j][i].z * zet[k][j][i].z);
	      else if (flg==2) 
		ucor[k][j][i].z -= correction*fabs(ucor[k][j][i].z);
	      else if (NumberOfBodies > 1) {
		ibi=(int)((nvert[k+1][j][i]-1.0)*1001);
		ucor[k][j][i].z -= sqrt(zet[k][j][i].x * zet[k][j][i].x +
					zet[k][j][i].y * zet[k][j][i].y +
					zet[k][j][i].z * zet[k][j][i].z) *
		  Correction[ibi];
	      }
	      else
		ucor[k][j][i].z -= sqrt(zet[k][j][i].x * zet[k][j][i].x +
					zet[k][j][i].y * zet[k][j][i].y +
					zet[k][j][i].z * zet[k][j][i].z) *
		  correction;
	    }
	  }
	}

	if (nvert[k][j][i] > 0.1 && nvert[k][j][i] < ibmval) {
	  if (nvert[k][j][i+1] < 0.1 && i < xend) {
	    if (fabs(ucor[k][j][i].x)>epsilon) {
	      if (flg==3) 
		ucor[k][j][i].x += correction*fabs(ucor[k][j][i].x)/
		  sqrt(csi[k][j][i].x * csi[k][j][i].x +
		       csi[k][j][i].y * csi[k][j][i].y +
		       csi[k][j][i].z * csi[k][j][i].z);
	      else if (flg==2) 
		ucor[k][j][i].x += correction*fabs(ucor[k][j][i].x);
	      else if (NumberOfBodies > 1) {
		ibi=(int)((nvert[k][j][i]-1.0)*1001);
		ucor[k][j][i].x += sqrt(csi[k][j][i].x * csi[k][j][i].x +
					csi[k][j][i].y * csi[k][j][i].y +
					csi[k][j][i].z * csi[k][j][i].z) *
		  Correction[ibi];
	      }
	      else
		ucor[k][j][i].x += sqrt(csi[k][j][i].x * csi[k][j][i].x +
					csi[k][j][i].y * csi[k][j][i].y +
					csi[k][j][i].z * csi[k][j][i].z) *
		  correction;
	    }
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < yend) {
	    if (fabs(ucor[k][j][i].y)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].y +=correction*fabs(ucor[k][j][i].y)/
		sqrt(eta[k][j][i].x * eta[k][j][i].x +
		     eta[k][j][i].y * eta[k][j][i].y +
		     eta[k][j][i].z * eta[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].y +=correction*fabs(ucor[k][j][i].y);
	    else if (NumberOfBodies > 1) {
	      ibi=(int)((nvert[k][j][i]-1.0)*1001);
	      ucor[k][j][i].y += sqrt(eta[k][j][i].x * eta[k][j][i].x +
				      eta[k][j][i].y * eta[k][j][i].y +
				      eta[k][j][i].z * eta[k][j][i].z) *
		Correction[ibi];
	    }
	    else
	    ucor[k][j][i].y += sqrt(eta[k][j][i].x * eta[k][j][i].x +
				    eta[k][j][i].y * eta[k][j][i].y +
				    eta[k][j][i].z * eta[k][j][i].z) *
				    correction;
	    }
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < zend) {
	    if (fabs(ucor[k][j][i].z)>epsilon) {
	    if (flg==3) 
	      ucor[k][j][i].z += correction*fabs(ucor[k][j][i].z)/
		sqrt(zet[k][j][i].x * zet[k][j][i].x +
		     zet[k][j][i].y * zet[k][j][i].y +
		     zet[k][j][i].z * zet[k][j][i].z);
	    else if (flg==2) 
	      ucor[k][j][i].z += correction*fabs(ucor[k][j][i].z);
	    else if (NumberOfBodies > 1) {
	      ibi=(int)((nvert[k][j][i]-1.0)*1001);
	      ucor[k][j][i].z += sqrt(zet[k][j][i].x * zet[k][j][i].x +
				      zet[k][j][i].y * zet[k][j][i].y +
				      zet[k][j][i].z * zet[k][j][i].z) *
		Correction[ibi];
	    }
	    else
	      ucor[k][j][i].z += sqrt(zet[k][j][i].x * zet[k][j][i].x +
				      zet[k][j][i].y * zet[k][j][i].y +
				      zet[k][j][i].z * zet[k][j][i].z) *
		correction;
	    }
	  }
	}
	
      }
    }
  }
  
  //================================================================================
  // PASS 3: Verification
  // This optional pass recalculates the flux to confirm the correction was successful.
  //================================================================================
  LOG_ALLOW(GLOBAL, LOG_DEBUG, "Pass 3: Verifying corrected flux.\n");
  
  libm_Flux = 0;
  libm_area = 0;
  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] > 0.1 && nvert[k][j][i+1] < ibmval && i < xend) {
	    libm_Flux += ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] > 0.1 && nvert[k][j+1][i] < ibmval && j < yend) {
	    libm_Flux += ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] > 0.1 && nvert[k+1][j][i] < ibmval && k < zend) {
	    libm_Flux += ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

	if (nvert[k][j][i] > 0.1 && nvert[k][j][i] < ibmval) {
	  if (nvert[k][j][i+1] < 0.1 && i < xend) {
	    libm_Flux -= ucor[k][j][i].x;
	    libm_area += sqrt(csi[k][j][i].x * csi[k][j][i].x +
			  csi[k][j][i].y * csi[k][j][i].y +
			  csi[k][j][i].z * csi[k][j][i].z);
			  
	  }
	  if (nvert[k][j+1][i] < 0.1 && j < yend) {
	    libm_Flux -= ucor[k][j][i].y;
	    libm_area += sqrt(eta[k][j][i].x * eta[k][j][i].x +
			  eta[k][j][i].y * eta[k][j][i].y +
			  eta[k][j][i].z * eta[k][j][i].z);
	  }
	  if (nvert[k+1][j][i] < 0.1 && k < zend) {
	    libm_Flux -= ucor[k][j][i].z;
	    libm_area += sqrt(zet[k][j][i].x * zet[k][j][i].x +
			  zet[k][j][i].y * zet[k][j][i].y +
			  zet[k][j][i].z * zet[k][j][i].z);
	  }
	}

      }
    }
  }

  MPI_Allreduce(&libm_Flux, ibm_Flux,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);
  MPI_Allreduce(&libm_area, ibm_Area,1,MPI_DOUBLE,MPI_SUM, PETSC_COMM_WORLD);

 /*  PetscGlobalSum(&libm_Flux, ibm_Flux, PETSC_COMM_WORLD); */
/*   PetscGlobalSum(&libm_area, ibm_Area, PETSC_COMM_WORLD); */
    LOG_ALLOW(GLOBAL, LOG_INFO, "IBM Corrected (Verified) Flux: %g, Area: %g\n", *ibm_Flux, *ibm_Area);


  if (user->bctype[0]==7 || user->bctype[1]==7){
    if (xe==mx){
      i=mx-2;
      for (k=lzs; k<lze; k++) {
	for (j=lys; j<lye; j++) {
	  // if(j>0 && k>0 && j<user->JM && k<user->KM){
	    if ((nvert[k][j][i]>ibmval && nvert[k][j][i+1]<0.1) || (nvert[k][j][i]<0.1 && nvert[k][j][i+1]>ibmval)) ucor[k][j][i].x=0.0;
	    
	    // }
	}
      }
    }
  }

  if (user->bctype[2]==7 || user->bctype[3]==7){
    if (ye==my){
      j=my-2;
      for (k=lzs; k<lze; k++) {
	for (i=lxs; i<lxe; i++) {
	  // if(i>0 && k>0 && i<user->IM && k<user->KM){
	    if ((nvert[k][j][i]>ibmval && nvert[k][j+1][i]<0.1) || (nvert[k][j][i]<0.1 && nvert[k][j+1][i]>ibmval)) ucor[k][j][i].y=0.0;
	    // }
	}
      }
    }
  }

  if (user->bctype[4]==7 || user->bctype[5]==7){
    if (ze==mz){
      k=mz-2;
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  // if(i>0 && j>0 && i<user->IM && j<user->JM){
	    if ((nvert[k][j][i]>ibmval && nvert[k+1][j][i]<0.1) || (nvert[k][j][i]<0.1 && nvert[k+1][j][i]>ibmval)) ucor[k][j][i].z=0.0;
	    // }
	}
      }
    }
  }


  DMDAVecRestoreArray(da, user->lNvert, &nvert);
  DMDAVecRestoreArray(fda, user->lCsi, &csi);
  DMDAVecRestoreArray(fda, user->lEta, &eta);
  DMDAVecRestoreArray(fda, user->lZet, &zet);
  DMDAVecRestoreArray(fda, user->Ucont, &ucor);

  DMGlobalToLocalBegin(fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(fda, user->Ucont, INSERT_VALUES, user->lUcont);

  /* periodci boundary condition update corrected flux */
  //Mohsen Dec 2015
  DMDAVecGetArray(fda, user->lUcont, &lucor);
  DMDAVecGetArray(fda, user->Ucont, &ucor);
  
  if (user->bctype[0]==7 || user->bctype[1]==7){
    if (xs==0){
      i=xs;
      for (k=zs; k<ze; k++) {
	for (j=ys; j<ye; j++) {
	  if(j>0 && k>0 && j<user->JM && k<user->KM){
	    ucor[k][j][i].x=lucor[k][j][i-2].x;  
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
	    ucor[k][j][i].y=lucor[k][j-2][i].y;
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
	    ucor[k][j][i].z=lucor[k-2][j][i].z;
	  }
	}
      }
    }
  }
  DMDAVecRestoreArray(fda, user->lUcont, &lucor);
  DMDAVecRestoreArray(fda, user->Ucont, &ucor);

 /*  DMLocalToGlobalBegin(user->fda, user->lUcont, INSERT_VALUES, user->Ucont); */
/*   DMLocalToGlobalEnd(user->fda, user->lUcont, INSERT_VALUES, user->Ucont); */
  DMGlobalToLocalBegin(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);
  DMGlobalToLocalEnd(user->fda, user->Ucont, INSERT_VALUES, user->lUcont);

  if (NumberOfBodies > 1) {
    free(lIB_Flux);
    free(lIB_area);
    free(IB_Flux);
    free(IB_Area);
    free(Correction);
  }

 LOG_ALLOW(GLOBAL, LOG_DEBUG, "Exiting VolumeFlux.\n");
  
  return 0;
}

PetscErrorCode FullyBlocked(UserCtx *user)
{
  DM da = user->da;
  Vec nNvert;
  DMDALocalInfo info = user->info;

  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt i, j, k;

  PetscInt *KSKE = user->KSKE;
  PetscReal ***nvert;
  PetscBool *Blocked;

  DMDACreateNaturalVector(da, &nNvert);
  DMDAGlobalToNaturalBegin(da, user->Nvert, INSERT_VALUES, nNvert);
  DMDAGlobalToNaturalEnd(da, user->Nvert, INSERT_VALUES, nNvert);

  VecScatter ctx;
  Vec Zvert;
  VecScatterCreateToZero(nNvert, &ctx, &Zvert);

  VecScatterBegin(ctx, nNvert, Zvert, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(ctx, nNvert, Zvert, INSERT_VALUES, SCATTER_FORWARD);

  VecScatterDestroy(&ctx);
  VecDestroy(&nNvert);

  PetscInt rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

  if (!rank) {

    VecGetArray3d(Zvert, mz, my, mx, 0, 0, 0, &nvert);
    PetscMalloc(mx*my*sizeof(PetscBool), &Blocked);
    for (j=1; j<my-1; j++) {
      for (i=1; i<mx-1; i++) {
	Blocked[j*mx+i] = PETSC_FALSE;
	for (k=0; k<mz; k++) {
	  if (nvert[k][j][i] > 0.1) {
	    if (!Blocked[j*mx+i]) {
	      KSKE[2*(j*mx+i)] = k;
	      Blocked[j*mx+i] = PETSC_TRUE;
	    }
	    else {
	      KSKE[2*(j*mx+i)] = PetscMin(KSKE[2*(j*mx+i)], k);
	    }
	  }
	}
      }
    }


    user->multinullspace = PETSC_TRUE;
    for (j=1; j<my-1; j++) {
      for (i=1; i<mx-1; i++) {
	if (!Blocked[j*mx+i]) {
	  user->multinullspace = PETSC_FALSE;
	  break;
	}
      }
    }
    PetscFree(Blocked);
    VecRestoreArray3d(Zvert, mz, my, mx, 0, 0, 0, &nvert);
    MPI_Bcast(&user->multinullspace, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    if (user->multinullspace) {
      MPI_Bcast(user->KSKE, 2*mx*my, MPI_INT, 0, PETSC_COMM_WORLD);

    }
  }
  else {
    MPI_Bcast(&user->multinullspace, 1, MPI_INT, 0, PETSC_COMM_WORLD);
    if (user->multinullspace) {
      MPI_Bcast(user->KSKE, 2*mx*my, MPI_INT, 0, PETSC_COMM_WORLD);
    }
  }



  VecDestroy(&Zvert);
  return 0;
}

PetscErrorCode MyNvertRestriction(UserCtx *user_h, UserCtx *user_c)
{
  //  DA		da = user_c->da, fda = user_c->fda;



  DMDALocalInfo	info = user_c->info;
  PetscInt	xs = info.xs, xe = info.xs + info.xm;
  PetscInt  	ys = info.ys, ye = info.ys + info.ym;
  PetscInt	zs = info.zs, ze = info.zs + info.zm;
  PetscInt	mx = info.mx, my = info.my, mz = info.mz;

  PetscInt i,j,k;
  PetscInt ih, jh, kh, ia, ja, ka;
  PetscInt	lxs, lxe, lys, lye, lzs, lze;

  PetscReal ***nvert, ***nvert_h;

  DMDAVecGetArray(user_h->da, user_h->lNvert, &nvert_h);
  DMDAVecGetArray(user_c->da, user_c->Nvert, &nvert);

  lxs = xs; lxe = xe;
  lys = ys; lye = ye;
  lzs = zs; lze = ze;

  if (xs==0) lxs = xs+1;
  if (ys==0) lys = ys+1;
  if (zs==0) lzs = zs+1;

  if (xe==mx) lxe = xe-1;
  if (ye==my) lye = ye-1;
  if (ze==mz) lze = ze-1;

  if ((user_c->isc)) ia = 0;
  else ia = 1;

  if ((user_c->jsc)) ja = 0;
  else ja = 1;

  if ((user_c->ksc)) ka = 0;
  else ka = 1;

  VecSet(user_c->Nvert, 0.);
  if (user_c->thislevel > 0) {
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  GridRestriction(i, j, k, &ih, &jh, &kh, user_c);
	  if (nvert_h[kh   ][jh   ][ih   ] *
	      nvert_h[kh   ][jh   ][ih-ia] *
	      nvert_h[kh   ][jh-ja][ih   ] *
	      nvert_h[kh-ka][jh   ][ih   ] *
	      nvert_h[kh   ][jh-ja][ih-ia] *
	      nvert_h[kh-ka][jh   ][ih-ia] *
	      nvert_h[kh-ka][jh-ja][ih   ] *
	      nvert_h[kh-ka][jh-ja][ih-ia] > 0.1) {
	    nvert[k][j][i] = PetscMax(1., nvert[k][j][i]);
	  }
	}
      }
    }
  }
  else {
    for (k=lzs; k<lze; k++) {
      for (j=lys; j<lye; j++) {
	for (i=lxs; i<lxe; i++) {
	  GridRestriction(i, j, k, &ih, &jh, &kh, user_c);
	  if (nvert_h[kh   ][jh   ][ih   ] *
	      nvert_h[kh   ][jh   ][ih-ia] *
	      nvert_h[kh   ][jh-ja][ih   ] *
	      nvert_h[kh-ka][jh   ][ih   ] *
	      nvert_h[kh   ][jh-ja][ih-ia] *
	      nvert_h[kh-ka][jh   ][ih-ia] *
	      nvert_h[kh-ka][jh-ja][ih   ] *
	      nvert_h[kh-ka][jh-ja][ih-ia] > 0.1) {
	    nvert[k][j][i] = PetscMax(1., nvert[k][j][i]);
	  }
	}
      }
    }
  }
  DMDAVecRestoreArray(user_h->da, user_h->lNvert, &nvert_h);
  DMDAVecRestoreArray(user_c->da, user_c->Nvert, &nvert);

  DMGlobalToLocalBegin(user_c->da, user_c->Nvert, INSERT_VALUES, user_c->lNvert);
  DMGlobalToLocalEnd(user_c->da, user_c->Nvert, INSERT_VALUES, user_c->lNvert);
  //Mohsen Dec 2015
  DMDAVecGetArray(user_c->da, user_c->lNvert, &nvert);
  DMDAVecGetArray(user_c->da, user_c->Nvert, &nvert_h);

  for (k=lzs; k<lze; k++) {
    for (j=lys; j<lye; j++) {
      for (i=lxs; i<lxe; i++) {
	if (nvert_h[k][j][i] < 0.1) {
	  if (nvert[k][j][i+1] + nvert[k][j][i-1] > 1.1 &&
	      nvert[k][j+1][i] + nvert[k][j-1][i] > 1.1 &&
	      nvert[k+1][j][i] + nvert[k-1][j][i] > 1.1) {
	    nvert_h[k][j][i] = 1.;
	  }
	}
      }
    }
  }

  DMDAVecRestoreArray(user_c->da, user_c->lNvert, &nvert);
  DMDAVecRestoreArray(user_c->da, user_c->Nvert, &nvert_h);
  DMGlobalToLocalBegin(user_c->da, user_c->Nvert, INSERT_VALUES, user_c->lNvert);
  DMGlobalToLocalEnd(user_c->da, user_c->Nvert, INSERT_VALUES, user_c->lNvert);
 /*  DMLocalToGlobalBegin(user_c->da, user_c->lNvert, INSERT_VALUES, user_c->Nvert); */
/*   DMLocalToGlobalEnd(user_c->da, user_c->lNvert, INSERT_VALUES, user_c->Nvert); */
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "PoissonSolver_MG"
PetscErrorCode PoissonSolver_MG(UserMG *usermg)
{
    // --- CONTEXT ACQUISITION BLOCK ---
    // Get the master simulation context from the first block's UserCtx on the finest level.
    // This provides access to all former global variables.
    SimCtx *simCtx = usermg->mgctx[0].user[0].simCtx;

    // Create local variables to mirror the legacy globals for minimal code changes.
    const PetscInt block_number = simCtx->block_number;
    const PetscInt immersed = simCtx->immersed;
    const PetscInt MHV = simCtx->MHV;
    const PetscInt LV = simCtx->LV;
    PetscMPIInt rank = simCtx->rank;
    // --- END CONTEXT ACQUISITION BLOCK ---

    PetscErrorCode ierr;
    PetscInt l, bi;
    MGCtx *mgctx = usermg->mgctx;
    KSP mgksp, subksp;
    PC mgpc, subpc;
    UserCtx *user;

    PetscFunctionBeginUser; // Moved to after variable declarations
    LOG_ALLOW(GLOBAL, LOG_INFO, "Starting Multigrid Poisson Solve...\n");

    for (bi = 0; bi < block_number; bi++) {
        
        // ====================================================================
        //   SECTION: Immersed Boundary Specific Setup (Conditional)
        // ====================================================================
        if (immersed) {
            LOG_ALLOW(LOCAL, LOG_DEBUG, "Block %d: Performing IBM pre-solve setup (Nvert restriction, etc.).\n", bi);
            for (l = usermg->mglevels - 1; l > 0; l--) {
                mgctx[l].user[bi].multinullspace = PETSC_FALSE;
                MyNvertRestriction(&mgctx[l].user[bi], &mgctx[l-1].user[bi]);
            }
            // Coarsest level check for disconnected domains due to IBM
            l = 0;
            user = mgctx[l].user;
            ierr = PetscMalloc1(user[bi].info.mx * user[bi].info.my * 2, &user[bi].KSKE); CHKERRQ(ierr);
            FullyBlocked(&user[bi]);
        }
        

        l = usermg->mglevels - 1;
        user = mgctx[l].user;
        
        // --- 1. Compute RHS of the Poisson Equation ---
        LOG_ALLOW(LOCAL, LOG_DEBUG, "Block %d: Computing Poisson RHS...\n", bi);
        ierr = VecDuplicate(user[bi].P, &user[bi].B); CHKERRQ(ierr);
        
        PetscReal ibm_Flux, ibm_Area;
        PetscInt flg = immersed - 1;

        // Calculate volume flux source terms (often from IBM)
        VolumeFlux(&user[bi], &ibm_Flux, &ibm_Area, flg);
        if (MHV || LV) {
            flg = ((MHV > 1 || LV) && bi == 0) ? 1 : 0;
            VolumeFlux_rev(&user[bi], &ibm_Flux, &ibm_Area, flg);
        }
        // Calculate the main divergence term
        PoissonRHS(&user[bi], user[bi].B);
        
        // --- 2. Assemble LHS Matrix (Laplacian) on all MG levels ---
        LOG_ALLOW(LOCAL, LOG_DEBUG, "Block %d: Assembling Poisson LHS on all levels...\n", bi);
        for (l = usermg->mglevels - 1; l >= 0; l--) {
            user = mgctx[l].user;
	    LOG_ALLOW(GLOBAL,LOG_DEBUG," Calculating LHS for Level %d.\n",l);
            PoissonLHSNew(&user[bi]);
        }

        // --- 3. Setup PETSc KSP and PCMG (Multigrid Preconditioner) ---
        LOG_ALLOW(LOCAL, LOG_DEBUG, "Block %d: Configuring KSP and PCMG...\n", bi);

	ierr = KSPCreate(PETSC_COMM_WORLD, &mgksp); CHKERRQ(ierr);
        ierr = KSPAppendOptionsPrefix(mgksp, "ps_"); CHKERRQ(ierr);

	// =======================================================================
        DualMonitorCtx *monctx;
        char           filen[128];
        PetscBool      has_monitor_option;

        // 1. Allocate the context and set it up.
        ierr = PetscNew(&monctx); CHKERRQ(ierr);

	monctx->step = simCtx->step;
	monctx->block_id = bi;
	monctx->file_handle = NULL;

	// Only rank 0 handles the file.
        if (!rank) {
	  sprintf(filen, "logs/Poisson_Solver_Convergence_History_Block_%d.log", bi);
	  // On the very first step of the entire simulation, TRUNCATE the file.
	  if (simCtx->step == simCtx->StartStep) {
	    monctx->file_handle = fopen(filen, "w");
	  } else { // For all subsequent steps, APPEND to the file.
	    monctx->file_handle = fopen(filen, "a");
	  }
 
            if (monctx->file_handle) {
                PetscFPrintf(PETSC_COMM_SELF, monctx->file_handle, "--- Convergence for Timestep %d, Block %d ---\n", (int)simCtx->step, bi);
            } else {
                SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_FILE_OPEN, "Could not open KSP monitor log file: %s", filen);
            }
        }

        ierr = PetscOptionsHasName(NULL, NULL, "-ps_ksp_pic_monitor_true_residual", &has_monitor_option); CHKERRQ(ierr);
        monctx->log_to_console = has_monitor_option;

        ierr = KSPMonitorSet(mgksp, DualKSPMonitor, monctx, DualMonitorDestroy); CHKERRQ(ierr);
        // =======================================================================
	
        ierr = KSPGetPC(mgksp, &mgpc); CHKERRQ(ierr);
        ierr = PCSetType(mgpc, PCMG); CHKERRQ(ierr);

	ierr = PCMGSetLevels(mgpc, usermg->mglevels, PETSC_NULL); CHKERRQ(ierr);
        ierr = PCMGSetCycleType(mgpc, PC_MG_CYCLE_V); CHKERRQ(ierr);
        ierr = PCMGSetType(mgpc, PC_MG_MULTIPLICATIVE); CHKERRQ(ierr);

        // --- 4. Define Restriction and Interpolation Operators for MG ---
        for (l = usermg->mglevels - 1; l > 0; l--) {

         // Get stable pointers directly from the main mgctx array.
         // These pointers point to memory that will persist.
	  UserCtx *fine_user_ctx   = &mgctx[l].user[bi];
	  UserCtx *coarse_user_ctx = &mgctx[l-1].user[bi];

	  // --- Configure the context pointers ---
	  // The coarse UserCtx needs to know about the fine grid for restriction.
	  coarse_user_ctx->da_f     = &(fine_user_ctx->da);
	  coarse_user_ctx->user_f   = fine_user_ctx;

	  // The fine UserCtx needs to know about the coarse grid for interpolation.
	  fine_user_ctx->da_c       = &(coarse_user_ctx->da);
	  fine_user_ctx->user_c     = coarse_user_ctx;
	  fine_user_ctx->lNvert_c   = &(coarse_user_ctx->lNvert);

	  // --- Get matrix dimensions ---
	  PetscInt m_c = (coarse_user_ctx->info.xm * coarse_user_ctx->info.ym * coarse_user_ctx->info.zm);
	  PetscInt m_f = (fine_user_ctx->info.xm * fine_user_ctx->info.ym * fine_user_ctx->info.zm);
	  PetscInt M_c = (coarse_user_ctx->info.mx * coarse_user_ctx->info.my * coarse_user_ctx->info.mz);
	  PetscInt M_f = (fine_user_ctx->info.mx * fine_user_ctx->info.my * fine_user_ctx->info.mz);

	  LOG_ALLOW(GLOBAL,LOG_DEBUG,"level = %d; m_c = %d; m_f = %d; M_c = %d; M_f = %d.\n",l,m_c,m_f,M_c,M_f);
	  // --- Create the MatShell objects ---
	  // Pass the STABLE pointer coarse_user_ctx as the context for restriction.
	  ierr = MatCreateShell(PETSC_COMM_WORLD, m_c, m_f, M_c, M_f, (void*)coarse_user_ctx, &fine_user_ctx->MR); CHKERRQ(ierr);
    
	  // Pass the STABLE pointer fine_user_ctx as the context for interpolation.
	  ierr = MatCreateShell(PETSC_COMM_WORLD, m_f, m_c, M_f, M_c, (void*)fine_user_ctx, &fine_user_ctx->MP); CHKERRQ(ierr);
    
	  // --- Set the operations for the MatShells ---
	  ierr = MatShellSetOperation(fine_user_ctx->MR, MATOP_MULT, (void(*)(void))RestrictResidual_SolidAware); CHKERRQ(ierr);
	  ierr = MatShellSetOperation(fine_user_ctx->MP, MATOP_MULT, (void(*)(void))MyInterpolation); CHKERRQ(ierr);
    
	  // --- Register the operators with PCMG ---
	  ierr = PCMGSetRestriction(mgpc, l, fine_user_ctx->MR); CHKERRQ(ierr);
	  ierr = PCMGSetInterpolation(mgpc, l, fine_user_ctx->MP); CHKERRQ(ierr);
	  
        }
        
        // --- 5. Configure Solvers on Each MG Level ---
        for (l = usermg->mglevels - 1; l >= 0; l--) {
            user = mgctx[l].user;
            if (l > 0) { // Smoother for fine levels
                ierr = PCMGGetSmoother(mgpc, l, &subksp); CHKERRQ(ierr);
            } else { // Direct or iterative solver for the coarsest level
                ierr = PCMGGetCoarseSolve(mgpc, &subksp); CHKERRQ(ierr);
                ierr = KSPSetTolerances(subksp, 1.e-8, PETSC_DEFAULT, PETSC_DEFAULT, 40); CHKERRQ(ierr);
            } 
            
            ierr = KSPSetOperators(subksp, user[bi].A, user[bi].A); CHKERRQ(ierr);
	    ierr = KSPSetFromOptions(subksp); CHKERRQ(ierr); 
            ierr = KSPGetPC(subksp, &subpc); CHKERRQ(ierr);
	    ierr = PCSetType(subpc, PCBJACOBI); CHKERRQ(ierr);
	    
	    KSP *subsubksp;
	    PC subsubpc;
	    PetscInt nlocal;

	    // This logic is required for both the smoother and the coarse solve
	    // since both use PCBJACOBI.
	    ierr = KSPSetUp(subksp); CHKERRQ(ierr); // Set up KSP to allow access to sub-KSPs
	    ierr = PCBJacobiGetSubKSP(subpc, &nlocal, NULL, &subsubksp); CHKERRQ(ierr);

	    for (PetscInt abi = 0; abi < nlocal; abi++) {
	      ierr = KSPGetPC(subsubksp[abi], &subsubpc); CHKERRQ(ierr);
	      // Add the critical shift amount
	      ierr = PCFactorSetShiftAmount(subsubpc, 1.e-10); CHKERRQ(ierr);
	    }
	    
	    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &user[bi].nullsp); CHKERRQ(ierr);
            ierr = MatNullSpaceSetFunction(user[bi].nullsp, PoissonNullSpaceFunction, &user[bi]); CHKERRQ(ierr);
            ierr = MatSetNullSpace(user[bi].A, user[bi].nullsp); CHKERRQ(ierr);
            
            ierr = PCMGSetResidual(mgpc, l, PCMGResidualDefault, user[bi].A); CHKERRQ(ierr);
            ierr = KSPSetUp(subksp); CHKERRQ(ierr);

            if (l < usermg->mglevels - 1) {
                ierr = MatCreateVecs(user[bi].A, &user[bi].R, PETSC_NULL); CHKERRQ(ierr);
                ierr = PCMGSetRhs(mgpc, l, user[bi].R); CHKERRQ(ierr);
            }
        }

        // --- 6. Set Final KSP Operators and Solve ---
        l = usermg->mglevels - 1;
        user = mgctx[l].user;
        
        LOG_ALLOW(LOCAL, LOG_DEBUG, "Block %d: Setting KSP operators and solving...\n", bi);
        ierr = KSPSetOperators(mgksp, user[bi].A, user[bi].A); CHKERRQ(ierr);
        ierr = MatSetNullSpace(user[bi].A, user[bi].nullsp); CHKERRQ(ierr);
        ierr = KSPSetFromOptions(mgksp); CHKERRQ(ierr);
        ierr = KSPSetUp(mgksp); CHKERRQ(ierr);
        ierr = KSPSolve(mgksp, user[bi].B, user[bi].Phi); CHKERRQ(ierr);

        // --- 7. Cleanup for this block ---
        for (l = usermg->mglevels - 1; l >= 0; l--) {
            user = mgctx[l].user;
            MatNullSpaceDestroy(&user[bi].nullsp);
            MatDestroy(&user[bi].A);
            user[bi].assignedA = PETSC_FALSE;
            if (l > 0) {
                MatDestroy(&user[bi].MR);
                MatDestroy(&user[bi].MP);
            } else if (l==0 && immersed) {
                PetscFree(user[bi].KSKE);
            }
            if (l < usermg->mglevels - 1) {
                VecDestroy(&user[bi].R);
            }
        }
        
        KSPDestroy(&mgksp);
        VecDestroy(&mgctx[usermg->mglevels-1].user[bi].B);
        
    } // End of loop over blocks

    LOG_ALLOW(GLOBAL, LOG_INFO, "Multigrid Poisson Solve complete.\n");
    PetscFunctionReturn(0);
}
