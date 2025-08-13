#include "wallfunction.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* double integrate_F(double nu, double utau, double ya, double yb); */
/* double integrate_Fy(double nu, double utau, double ya, double yb); */
/* double sign(double a); */

#define KAPPA 0.41
#define LOGLAW_B 5.5

void noslip (UserCtx *user, double sc, double sb, Cmpnts Ua, Cmpnts Uc, 
	     Cmpnts *Ub, double nx, double ny, double nz)
{
	double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
	double un = u_c * nx + v_c * ny + w_c * nz;
	double ut = u_c - un * nx, vt = v_c - un * ny, wt = w_c - un * nz;
	double ut_mag = sqrt( ut*ut + vt*vt + wt*wt );
	
	(*Ub).x = sb/sc * u_c;
	(*Ub).y = sb/sc * v_c;
	(*Ub).z = sb/sc * w_c;
	
	(*Ub).x += Ua.x;
	(*Ub).y += Ua.y;
	(*Ub).z += Ua.z;
}

void freeslip (UserCtx *user, double sc, double sb, Cmpnts Ua, Cmpnts Uc, 
	       Cmpnts *Ub, double nx, double ny, double nz)
{
  //  in the Normal direction interpolation 
  double Uan = Ua.x * nx + Ua.y * ny + Ua.z * nz;
  double Ucn = Uc.x * nx + Uc.y * ny + Uc.z * nz;
  double Ubn = Uan + (Ucn - Uan) * sb/sc;
  //   in the Tangential direction Ub_t = Uc_t
  double ut = Uc.x - Ucn * nx, vt = Uc.y - Ucn * ny, wt = Uc.z - Ucn * nz;
  
	(*Ub).x = ut;
	(*Ub).y = vt;
	(*Ub).z = wt;
	
	(*Ub).x += Ubn*nx;
	(*Ub).y += Ubn*ny;
	(*Ub).z += Ubn*nz;
}

/* ali added the Wall function for rough and smooth walls using Log-low   */ 
double E_coeff (double utau, double ks, double nu)
{
  double kplus=utau*ks/nu, dB;
  if(kplus<=2.25) dB = 0.0;
  else if(kplus>2.25 && kplus<90.) dB = (LOGLAW_B-8.5+1./KAPPA*log(kplus))*sin(0.4258*(log(kplus)-0.811));
  else if(kplus>=90.) dB = LOGLAW_B-8.5+1./KAPPA*log(kplus);
  return exp(KAPPA*(LOGLAW_B-dB));
}

double u_hydset_roughness(double nu, double y, double utau, double ks)
{
  double y0plus=11.53, y1plus=300,f;
  double yplus = utau*y/nu;
 
  if(yplus<=y0plus){f= utau * yplus;}
  else if (yplus<=y1plus) {f= utau/KAPPA*log(E_coeff(utau,ks,nu)*yplus);}
  else {f=-1.;}
  return f;
}

double f_hydset(double nu, double u, double y, double utau0, double ks) 
{
double y0plus=11.53, f;
double yplus=utau0*y/nu;
 if (yplus<=y0plus) {f= utau0*yplus-u;}
 else {f= utau0*(1./KAPPA*log(E_coeff(utau0,ks,nu)*yplus))-u;}
return f;
}

double df_hydset (double nu, double u, double y, double utau0, double ks)
{
double eps=1.e-7;
return (f_hydset(nu, u, y, utau0 + eps, ks) - f_hydset(nu, u, y, utau0 - eps, ks))/(2.*eps);
}


double find_utau_hydset(double nu,double u, double y, double utau_guess, double ks)
{
 double utau,utau0 = utau_guess;  
 int ii;
 for(ii=0;ii<30;ii++){
 utau=utau0-f_hydset(nu,u,y,utau0,ks)/df_hydset(nu,u,y,utau0,ks);
 if (fabs(utau-utau0)<1.e-7)break;
 utau0=utau;
  }
 if (ii==30) LOG_ALLOW(LOCAL,LOG_DEBUG, "iter utau diff %d %le %le  %le\n",ii,utau,utau0,fabs(utau-utau0));
 //printf("iter utau diff %d %le %le  %le\n",ii,utau,utau0,fabs(utau-utau0));
return utau;
};

//printf("iter  yplus  utau diff %d %le %le  %le\n",i,yplus,utau, fabs(u-uc));


double nu_t(double yplus)// in fact, nu_t/nu
{
	return KAPPA * yplus * pow ( 1. - exp( - yplus / 19. ) , 2.0 );
};
double integrate_1(double nu,double y,double utau, int m)
{
	int ny=30;	//number of integartion points
	int i;
	double nut,f=0.0,yplus;
	double dy=y/(ny-1);
	double F[ny];
	for (i=0;i<ny;i++){
		yplus=(i*dy)*utau/nu;
		nut=nu_t(yplus);
		if (m==0){
		 F[i]=1.0/(1.0+nut)/nu;
		}	else
		{
			F[i]=dy*i/(1.0+nut)/nu;
		}
	}
	for (i=1;i<ny-1;i++){
		f+=(F[i-1]+2*F[i]+F[i+1])*0.25*dy;
	}
	f+=.1667*dy*(2*F[0]+F[1]+2*F[ny-1]+F[ny-2]);
	return f;
}
double taw(double nu, double utau, double y, double u, double dpdt)
{
	double yplus = y*utau/nu;
	double sum=0.0;
	double f1 = integrate_1(nu,y,utau,0);
	double f2 = integrate_1(nu,y,utau,1);
	return 1/f1*(u-dpdt*f2); 
}
double u_Cabot(double nu, double y, double utau, double dpdt, double taw)
{
	double f1 = integrate_1(nu,y,utau,0);
	double f2 = integrate_1(nu,y,utau,1);
	return  taw*f1+dpdt*f2;
};

double u_Werner(double nu, double y, double utau)
{
	double yplus = utau * y / nu; 
	double A=8.3, B=1./7.;
	
	if(yplus<=11.81) return yplus*utau;
	else return A * pow(yplus, B) * utau;
};

double f_Werner(double nu, double u, double y, double utau)
{
	double ym=11.81*nu/utau;
	double A=8.3, B=1./7.;
	
	if( fabs(u) <= nu/(2*ym) * pow(A, 2./(1.-B) ) ) return utau*utau - u/y*nu;
	else return utau*utau -  u/fabs(u) * pow( 0.5*(1-B) * pow(A, (1+B)/(1-B)) * pow(nu/y, 1+B) + (1+B)/A * pow(nu/y, B) * fabs(u), 2/(1+B) );
}

double df_Werner(double nu, double u, double y, double utau)
{
	double eps=1.e-7;
	return ( f_Werner(nu, u, y, utau+eps) - f_Werner(nu, u, y, utau-eps) ) / ( 2*eps ) ;
}

double f_Cabot(double nu, double u, double y, double utau, double dpdt, double dpdtn)
{
	return utau - sqrt(fabs(taw( nu, utau, y, u, dpdt )));
}

double df_Cabot(double nu, double u, double y, double utau, double dpdt, double dpdtn)
{
	double eps=1.e-7;
	return ( f_Cabot(nu, u, y, utau+eps, dpdt,dpdtn) - f_Cabot(nu, u, y, utau-eps, dpdt, dpdtn)) / ( 2*eps ) ;
}

/*
	integrate F dy 	= integrate [ 1/ (nut + nu) ] dy
				= integrate [ 1/ (nu) / {1 + k*yp( 1- exp(-yp/A) )^2 }  ] dy
				= integrate [ 1/ utau / {1 + k*yp( 1- exp(-yp/A) )^2 }  ] d(yp)
				= 1/ utau * integrate [  1 / {1 + k*yp( 1- exp(-yp/A)^2 ) }  ] d(yp)

	integrate Fy dy	= y/ utau * integrate [  1 / {1 + k*yp( 1- exp(-yp/A)^2 ) }  ] d(yp)
				= nu/utau^2 * integrate [  yp / {1 + k*yp( 1- exp(-yp/A)^2 ) }  ] d(yp)
*/


void find_utau_Cabot(double nu, double u, double y, double guess, double dpdt, double dpdtn, double *utau, double *taw1, double *taw2)
{
	double x, x0=guess;
	
	int i;
	
	for(i=0; i<30; i++) {
		x =fabs (x0 - f_Cabot(nu, u, y, x0, dpdt, dpdtn)/df_Cabot(nu, u, y, x0, dpdt, dpdtn));
		if( fabs(x0 - x) < 1.e-10) break;
		x0 = x;
	}
	
	if( fabs(x0 - x) > 1.e-5 && i>=29 ) printf("\n!!!!!!!! Iteration Failed !!!!!!!!\n");
		*utau=x;
		*taw1=taw( nu, x, y, u, dpdt );
	    *taw2=taw( nu, x, y, 0.0, dpdtn );

};

double find_utau_Werner(double nu, double u, double y, double guess)
{
	double x, x0=guess;
	
	int i;
	
	for(i=0; i<20; i++) {
		x = x0 - f_Werner(nu, u, y, x0)/df_Werner(nu, u, y, x0);
		if( fabs(x0 - x) < 1.e-7 ) break;
		x0 = x;
	}
	
	if( fabs(x0 - x) > 1.e-5 && i>=19 ) printf("\n!!!!!!!! Iteration Failed !!!!!!!!\n");
		
	return x;
};

double sign(double a)
{
	if(a>0) return 1;
	else if(a<0) return -1;
	else return 0;
}


double u_loglaw(double y, double utau, double roughness)
{
	return utau *  1./KAPPA * log( (roughness+y) / roughness ) ;
}

double find_utau_loglaw(double u, double y, double roughness)
{
	return KAPPA * u / log( (y+roughness) / roughness);
};

void wall_function (UserCtx *user, double sc, double sb, 
		    Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, 
		    double nx, double ny, double nz)
{
	SimCtx *simCtx = user->simCtx;
        double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
	double un = u_c * nx + v_c * ny + w_c * nz;
	double ut = u_c - un * nx, vt = v_c - un * ny, wt = w_c - un * nz;
	double ut_mag = sqrt( ut*ut + vt*vt + wt*wt );
	//*ustar = find_utau_Cabot(1./simCtx->ren, ut_mag, sc, 0.01, 0);
	double utau=0.05;
	double ut_mag_modeled = u_Werner(1./simCtx->ren, sb, utau);
	//double ut_mag_modeled = u_Cabot(1./simCtx->ren, sb, *ustar, 0);
	
	if(ut_mag>1.e-10) {
		ut *= ut_mag_modeled/ut_mag;
		vt *= ut_mag_modeled/ut_mag;
		wt *= ut_mag_modeled/ut_mag;
	}
	else ut=vt=wt=0;
					
	// u = ut + (u.n)n
	(*Ub).x = ut + sb/sc * un * nx;
	(*Ub).y = vt + sb/sc * un * ny;
	(*Ub).z = wt + sb/sc * un * nz;
	
	(*Ub).x += Ua.x;
	(*Ub).y += Ua.y;
	(*Ub).z += Ua.z;
}

void wall_function_loglaw (UserCtx *user, double ks, double sc, double sb, 
			   Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, 
			   double nx, double ny, double nz)
{
	SimCtx *simCtx = user->simCtx;
        double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
	//	double ua_n= Ua.x*nx + Ua.y* ny + Ua.z* nz;
	double un = u_c * nx + v_c * ny + w_c * nz;
	double ut = u_c - un * nx, vt = v_c - un * ny, wt = w_c - un * nz;
	double ut_mag = sqrt( ut*ut + vt*vt + wt*wt );
	//*ustar = find_utau_Cabot_roughness(1./simCtx->ren, ut_mag, sc, 0.01, 0, ks);
        *ustar = find_utau_hydset(1./simCtx->ren, ut_mag, sc, 0.001, ks);
	//double ut_mag_modeled = u_Cabot_roughness(1./simCtx->ren, sb, *ustar, 0, ks);
	double ut_mag_modeled = u_hydset_roughness(1./simCtx->ren, sb, *ustar, ks);
	
	if (ut_mag_modeled<0.) {
	  //	  freeslip(user, sc,sb,Ua, Uc,Ub,nx,ny,nz);
	  ut_mag_modeled=ut_mag;
	  //*ustar =0.;
	  //  LOG_ALLOW(LOCAL,LOG_DEBUG, "Y+>300!!!!!!!!!!!!!!! %le ",*ustar);
	}
	
	//PetscReal yp=sc* (*ustar)*simCtx->ren;
	//if (yp>300) LOG_ALLOW(LOCAL,LOG_DEBUG, "  Y+>300!!!!!!!!!!!!!!!  %le ",yp);
	
	if(ut_mag>1.e-10) {
		ut *= ut_mag_modeled/ut_mag;
		vt *= ut_mag_modeled/ut_mag;
		wt *= ut_mag_modeled/ut_mag;
	}
	else ut=vt=wt=0;
					
	// u = ut + (u.n)n
	(*Ub).x = ut + (sb/sc * un ) * nx;
	(*Ub).y = vt + (sb/sc * un ) * ny;
	(*Ub).z = wt + (sb/sc * un ) * nz;
	
	(*Ub).x += Ua.x;
	(*Ub).y += Ua.y;
	(*Ub).z += Ua.z;
	
}

void wall_function_Cabot (UserCtx *user, double ks, double sc, double sb, 
			   Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, 
			   double nx, double ny, double nz, double dpdx, double dpdy, double dpdz, int count)
	    {
	        SimCtx *simCtx = user->simCtx;
		double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
	   //	double ua_n= Ua.x*nx + Ua.y* ny + Ua.z* nz;
	        double un = u_c * nx + v_c * ny + w_c * nz;
	        double ut = u_c - un * nx, vt = v_c - un * ny, wt = w_c - un * nz;
			double ut_mag = sqrt( ut*ut + vt*vt + wt*wt );	   	   
			double dpdt = dpdx*ut/ut_mag + dpdy*vt/ut_mag + dpdz*wt/ut_mag;
			double dpdtn = sqrt(dpdx*dpdx+dpdy*dpdy+dpdz*dpdz-pow(dpdx*ut/ut_mag,2.0) - pow(dpdx*nx,2.0) - 
			pow(dpdy*vt/ut_mag,2.0) - pow(dpdy*ny,2.0) - pow(dpdz*nz,2.0)-pow(dpdz*wt/ut_mag,2.0));
			double utau,taw1,taw2;
			dpdtn=0.0; //dpdt=0.0;
		       if (count ==0 ||count > 4){
			find_utau_Cabot(1./simCtx->ren, ut_mag, sc, 0.01, dpdt, dpdtn, &utau , &taw1, &taw2);
			*ustar=utau;}
			//*ustar=0.045;
			double ut_mag_modeled = u_Cabot(1./simCtx->ren, sb, *ustar, dpdt, taw1);

				   	if(ut_mag>1.e-10) {
		ut *= ut_mag_modeled/ut_mag;
		vt *= ut_mag_modeled/ut_mag;
		wt *= ut_mag_modeled/ut_mag;
	     }
	      else ut=vt=wt=0;
					
	// u = ut + (u.n)n
	(*Ub).x = ut + (sb/sc * un ) * nx;
	(*Ub).y = vt + (sb/sc * un ) * ny;
	(*Ub).z = wt + (sb/sc * un ) * nz;
	
	(*Ub).x += Ua.x;
	(*Ub).y += Ua.y;
	(*Ub).z += Ua.z;
				   
	   }

