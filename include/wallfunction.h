#ifndef WALLFUNCTION_H
#define WALLFUNCTION_H

#include <petscpf.h>
#include <petscdmswarm.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <petsctime.h>
#include <petscsys.h>
#include <petscdmcomposite.h>
#include <petscsystypes.h>

// Include additional headers
#include "variables.h"         // Shared type definitions
#include "ParticleSwarm.h"  // Particle swarm functions
#include "walkingsearch.h"  // Particle location functions
#include "grid.h"           // Grid functions
#include "logging.h"        // Logging macros
#include "io.h"             // Data Input and Output functions
#include "interpolation.h"  // Interpolation routines
#include "AnalyticalSolution.h" // Analytical Solution for testing
#include "ParticleMotion.h" // Functions related to motion of particles
#include "Boundaries.h"     //  Functions related to Boundary condition

void noslip (UserCtx *user, double sc, double sb, Cmpnts Ua, Cmpnts Uc, 
	     Cmpnts *Ub, double nx, double ny, double nz);

void freeslip (UserCtx *user, double sc, double sb, Cmpnts Ua, Cmpnts Uc, 
	       Cmpnts *Ub, double nx, double ny, double nz);

double E_coeff (double utau, double ks, double nu);

double u_hydset_roughness(double nu, double y, double utau, double ks);

double f_hydset(double nu, double u, double y, double utau0, double ks);

double df_hydset (double nu, double u, double y, double utau0, double ks);

double find_utau_hydset(double nu,double u, double y, double utau_guess, double ks);

double nu_t(double yplus);

double integrate_1(double nu,double y,double utau, int m);

double taw(double nu, double utau, double y, double u, double dpdt);

double u_Cabot(double nu, double y, double utau, double dpdt, double taw);

double u_Werner(double nu, double y, double utau);

double f_Werner(double nu, double u, double y, double utau);

double df_Werner(double nu, double u, double y, double utau);

double f_Cabot(double nu, double u, double y, double utau, double dpdt, double dpdtn);

double df_Cabot(double nu, double u, double y, double utau, double dpdt, double dpdtn);

void find_utau_Cabot(double nu, double u, double y, double guess, double dpdt, double dpdtn, double *utau, double *taw1, double *taw2);

double find_utau_Werner(double nu, double u, double y, double guess);

double sign(double a);

double u_loglaw(double y, double utau, double roughness);

double find_utau_loglaw(double u, double y, double roughness);

void wall_function (UserCtx *user, double sc, double sb, 
		    Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, 
		    double nx, double ny, double nz);

void wall_function_loglaw (UserCtx *user, double ks, double sc, double sb, 
			   Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, 
			   double nx, double ny, double nz);

void wall_function_Cabot (UserCtx *user, double ks, double sc, double sb, 
			   Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, 
			  double nx, double ny, double nz, double dpdx, double dpdy, double dpdz, int count);


#endif // WALLFUNCTION_H
