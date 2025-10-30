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

// ===========================================================================
// Function Declarations
// ===========================================================================

/**
 * @brief Applies no-slip wall boundary condition with linear interpolation
 *
 * This function enforces a no-slip boundary condition (zero velocity at the wall)
 * by linearly interpolating between the wall velocity (typically zero) and the
 * velocity at a reference point in the flow.
 *
 * MATHEMATICAL FORMULATION:
 *   For a point at distance sb from the wall, with a reference velocity Uc
 *   at distance sc from the wall:
 *     U_boundary = U_wall + (U_reference - U_wall) * (sb / sc)
 *
 * PHYSICAL INTERPRETATION:
 *   This provides a first-order approximation assuming linear velocity variation
 *   in the near-wall region, which is valid in the viscous sublayer.
 *
 * @param[in]  user              Simulation context (unused but required for interface)
 * @param[in]  distance_reference Wall-normal distance to reference point (sc)
 * @param[in]  distance_boundary  Wall-normal distance to boundary point (sb)
 * @param[in]  velocity_wall      Velocity at the wall (Ua), typically zero
 * @param[in]  velocity_reference Velocity at reference point (Uc)
 * @param[out] velocity_boundary  Computed velocity at boundary point (Ub)
 * @param[in]  normal_x           X-component of wall normal vector
 * @param[in]  normal_y           Y-component of wall normal vector
 * @param[in]  normal_z           Z-component of wall normal vector
 *
 * @note For moving walls, velocity_wall would be non-zero
 * @note The normal vector should point INTO the fluid domain
 */
void noslip(UserCtx *user, 
            double distance_reference, double distance_boundary,
            Cmpnts velocity_wall, Cmpnts velocity_reference,
            Cmpnts *velocity_boundary,
            double normal_x, double normal_y, double normal_z);

/**
 * @brief Applies free-slip wall boundary condition
 *
 * Free-slip conditions allow tangential flow but enforce zero normal velocity.
 * This is appropriate for inviscid walls or symmetry planes where there is
 * no shear stress but flow cannot penetrate the boundary.
 *
 * DECOMPOSITION:
 *   Velocity is decomposed into normal and tangential components:
 *     U = U_n * n + U_t
 *   where U_n = U · n (normal component)
 *         U_t = U - U_n * n (tangential component)
 *
 * BOUNDARY CONDITIONS:
 *   - Normal component: Interpolated from interior (∂U_n/∂n ≠ 0)
 *   - Tangential component: Extrapolated from interior (∂U_t/∂n = 0)
 *
 * @param[in]  user              Simulation context
 * @param[in]  distance_reference Wall-normal distance to reference point
 * @param[in]  distance_boundary  Wall-normal distance to boundary point
 * @param[in]  velocity_wall      Velocity at the wall
 * @param[in]  velocity_reference Velocity at reference point
 * @param[out] velocity_boundary  Computed velocity at boundary point
 * @param[in]  normal_x           X-component of wall normal vector
 * @param[in]  normal_y           Y-component of wall normal vector
 * @param[in]  normal_z           Z-component of wall normal vector
 */
void freeslip(UserCtx *user,
              double distance_reference, double distance_boundary,
              Cmpnts velocity_wall, Cmpnts velocity_reference,
              Cmpnts *velocity_boundary,
              double normal_x, double normal_y, double normal_z);
			  
/**
 * @brief Computes roughness-modified log-law coefficient E
 *
 * The coefficient E accounts for wall roughness effects on the log-law:
 *   u+ = (1/κ) ln(E * y+)
 *
 * ROUGHNESS REGIMES:
 *   1. Hydraulically smooth (ks+ < 2.25): E = exp(κ*B), no roughness effect
 *   2. Transitional (2.25 < ks+ < 90): Smooth interpolation between regimes
 *   3. Fully rough (ks+ > 90): E modified by roughness height
 * where ks+ = u_τ * ks / ν (roughness Reynolds number)
 *
 * @param[in] friction_velocity Friction velocity u_τ
 * @param[in] roughness_height  Equivalent sand grain roughness height ks
 * @param[in] kinematic_viscosity Kinematic viscosity ν
 * @return    Roughness-modified log-law coefficient E
 *
 * @note This follows the Nikuradse sand grain roughness correlation
 */
double E_coeff(double friction_velocity, double roughness_height, double kinematic_viscosity);

/**
 * @brief Computes velocity from log-law for rough walls
 *
 * Calculates the tangential velocity at a given wall distance using
 * the roughness-modified log-law of the wall.
 *
 * APPLICABLE REGIMES:
 *   - Viscous sublayer (y+ < 11.53): u+ = y+
 *   - Log-layer (11.53 < y+ < 300): u+ = (1/κ) ln(E * y+)
 *   - Outer layer (y+ > 300): Returns -1 (invalid, wake region)
 *
 * @param[in] kinematic_viscosity Kinematic viscosity ν
 * @param[in] wall_distance       Distance from wall y
 * @param[in] friction_velocity   Friction velocity u_τ
 * @param[in] roughness_height    Equivalent sand grain roughness ks
 * @return    Tangential velocity u_t, or -1 if outside valid range
 *
 * @note The upper limit y+ < 300 ensures we stay within the log-layer
 * @note For y+ > 300, wake effects become important and log-law is invalid
 */
double u_hydset_roughness(double kinematic_viscosity, double wall_distance,
                          double friction_velocity, double roughness_height);

/**
 * @brief Residual function for friction velocity equation (log-law with roughness)
 *
 * This function computes the residual for Newton-Raphson iteration:
 *   f(u_τ) = u_predicted(u_τ) - u_known
 * where u_predicted comes from the log-law or linear law depending on y+.
 *
 * @param[in] kinematic_viscosity Kinematic viscosity
 * @param[in] known_velocity      Known velocity at distance y
 * @param[in] wall_distance       Distance from wall
 * @param[in] friction_velocity_guess Current guess for u_τ
 * @param[in] roughness_height    Wall roughness height
 * @return    Residual value (should be zero at solution)
 */
double f_hydset(double kinematic_viscosity, double known_velocity, 
                double wall_distance, double friction_velocity_guess, 
                double roughness_height);
				
/**
 * @brief Numerical derivative of residual function
 *
 * Computes df/du_τ using finite differences for Newton-Raphson iteration.
 *
 * @param[in] kinematic_viscosity Kinematic viscosity
 * @param[in] known_velocity      Known velocity
 * @param[in] wall_distance       Distance from wall
 * @param[in] friction_velocity_guess Current guess for u_τ
 * @param[in] roughness_height    Wall roughness
 * @return    Derivative df/du_τ
 */
double df_hydset(double kinematic_viscosity, double known_velocity,
                 double wall_distance, double friction_velocity_guess,
                 double roughness_height);
				 
/**
 * @brief Solves for friction velocity using Newton-Raphson iteration
 *
 * Given a known velocity at a known distance from the wall, this function
 * iteratively solves for the friction velocity u_τ that satisfies the
 * roughness-modified log-law or linear law.
 *
 * ALGORITHM:
 *   Newton-Raphson: u_τ^(n+1) = u_τ^n - f(u_τ^n) / f'(u_τ^n)
 *   Convergence criterion: |u_τ^(n+1) - u_τ^n| < 1e-7
 *
 * @param[in] kinematic_viscosity Kinematic viscosity
 * @param[in] known_velocity      Velocity at reference point
 * @param[in] wall_distance       Distance from wall to reference point
 * @param[in] initial_guess       Initial guess for u_τ
 * @param[in] roughness_height    Wall roughness height
 * @return    Converged friction velocity u_τ
 *
 * @note Maximum 30 iterations; warns if convergence not achieved
 * @note Typical convergence in 3-5 iterations for good initial guess
 */
double find_utau_hydset(double kinematic_viscosity, double known_velocity,
                        double wall_distance, double initial_guess,
                        double roughness_height);

						
/**
 * @brief Computes turbulent eddy viscosity ratio (ν_t / ν)
 *
 * Uses the mixing length model with van Driest damping:
 *   ν_t / ν = κ * y+ * [1 - exp(-y+ / A+)]²
 * where A+ ≈ 19 is the damping coefficient.
 *
 * PHYSICAL INTERPRETATION:
 *   - Near wall (y+ → 0): ν_t → 0 (viscous effects dominate)
 *   - Far from wall (y+ >> A+): ν_t ~ κ * y+ (mixing length theory)
 *   - Damping function smoothly transitions between these limits
 *
 * @param[in] yplus Normalized wall distance y+ = y * u_τ / ν
 * @return    Eddy viscosity ratio ν_t / ν
 *
 * @note This is the standard mixing length model with van Driest damping
 * @note Valid in the inner layer (y+ < 50-100)
 */
double nu_t(double yplus);

/**
 * @brief Integrates eddy viscosity profile from wall to distance y
 *
 * Computes integrals needed for non-equilibrium wall functions:
 *   If mode=0: ∫[0 to y] dy / (ν + ν_t)
 *   If mode=1: ∫[0 to y] y dy / (ν + ν_t)
 *
 * These integrals appear in the solution of the momentum equation
 * with pressure gradients in the near-wall region.
 *
 * NUMERICAL METHOD:
 *   - Trapezoidal rule with 30 integration points
 *   - Higher-order correction at endpoints for accuracy
 *
 * @param[in] kinematic_viscosity Kinematic viscosity ν
 * @param[in] wall_distance       Distance from wall y
 * @param[in] friction_velocity   Friction velocity u_τ
 * @param[in] integration_mode    0 for F integral, 1 for Fy integral
 * @return    Value of the integral
 *
 * @note Used in Cabot wall function for pressure gradient effects
 */
double integrate_1(double kinematic_viscosity, double wall_distance,
                   double friction_velocity, int integration_mode);

/**
 * @brief Computes wall shear stress with pressure gradient effects
 *
 * Solves the integrated momentum equation in the near-wall region:
 *   τ_w = (u - dp/dx * F_y) / F_1
 * where F_1 and F_y are integrals of the effective viscosity profile.
 *
 * This accounts for non-equilibrium effects due to streamwise pressure gradients.
 *
 * @param[in] kinematic_viscosity Kinematic viscosity
 * @param[in] friction_velocity   Friction velocity u_τ
 * @param[in] wall_distance       Distance from wall
 * @param[in] velocity            Velocity at wall_distance
 * @param[in] pressure_gradient_tangent Tangential pressure gradient dp/ds
 * @return    Wall shear stress τ_w
 *
 * @note This is more accurate than equilibrium wall functions in separated flows
 */
double taw(double kinematic_viscosity, double friction_velocity,
           double wall_distance, double velocity, double pressure_gradient_tangent);
		   
/**
 * @brief Computes velocity using Cabot wall function
 *
 * Reconstructs velocity from wall shear stress and pressure gradient:
 *   u = τ_w * F1 + (dp/dx) * Fy
 *
 * @param[in] kinematic_viscosity Kinematic viscosity
 * @param[in] wall_distance       Distance from wall
 * @param[in] friction_velocity   Friction velocity
 * @param[in] pressure_gradient_tangent Tangential pressure gradient
 * @param[in] wall_shear_stress   Wall shear stress
 * @return    Velocity at wall_distance
 */
double u_Cabot(double kinematic_viscosity, double wall_distance,
               double friction_velocity, double pressure_gradient_tangent,
               double wall_shear_stress);
			   
/**
 * @brief Residual function for Cabot wall function
 *
 * Computes: f(u_τ) = u_τ - sqrt(|τ_w(u_τ)|)
 * This enforces consistency between u_τ and τ_w.
 *
 * @param[in] kinematic_viscosity Kinematic viscosity
 * @param[in] velocity            Velocity
 * @param[in] wall_distance       Distance from wall
 * @param[in] friction_velocity_guess Guess for u_τ
 * @param[in] pressure_gradient_tangent Tangential pressure gradient
 * @param[in] pressure_gradient_normal Normal pressure gradient (currently unused)
 * @return    Residual value
 */
double f_Cabot(double kinematic_viscosity, double velocity, double wall_distance,
               double friction_velocity_guess, double pressure_gradient_tangent,
               double pressure_gradient_normal);
			   
/**
 * @brief Numerical derivative for Cabot wall function
 */
double df_Cabot(double kinematic_viscosity, double velocity, double wall_distance,
                double friction_velocity_guess, double pressure_gradient_tangent,
                double pressure_gradient_normal);
				
/**
 * @brief Solves for friction velocity using Cabot wall function
 *
 * This non-equilibrium wall function accounts for pressure gradient effects,
 * making it more accurate in separated or strongly accelerating/decelerating flows.
 *
 * @param[in]  kinematic_viscosity Kinematic viscosity
 * @param[in]  velocity            Velocity at reference point
 * @param[in]  wall_distance       Distance from wall
 * @param[in]  initial_guess       Initial guess for u_τ
 * @param[in]  pressure_gradient_tangent Tangential pressure gradient
 * @param[in]  pressure_gradient_normal Normal pressure gradient
 * @param[out] friction_velocity   Converged friction velocity
 * @param[out] wall_shear_velocity Wall shear stress for velocity
 * @param[out] wall_shear_normal   Wall shear stress for normal pressure gradient
 */
void find_utau_Cabot(double kinematic_viscosity, double velocity, double wall_distance,
                     double initial_guess, double pressure_gradient_tangent,
                     double pressure_gradient_normal, double *friction_velocity,
                     double *wall_shear_velocity, double *wall_shear_normal);
					 
/**
 * @brief Computes velocity using Werner-Wengle wall function
 *
 * Algebraic wall function that provides explicit relation:
 *   u+ = y+                    for y+ < 11.81 (viscous sublayer)
 *   u+ = A * (y+)^B            for y+ ≥ 11.81 (power law)
 * where A = 8.3, B = 1/7 are empirical constants.
 *
 * ADVANTAGES:
 *   - Explicit (no iteration required)
 *   - Computationally cheap
 *   - Reasonable accuracy for equilibrium boundary layers
 *
 * LIMITATIONS:
 *   - Less accurate than log-law for high Reynolds numbers
 *   - No pressure gradient effects
 *
 * @param[in] kinematic_viscosity Kinematic viscosity
 * @param[in] wall_distance       Distance from wall
 * @param[in] friction_velocity   Friction velocity
 * @return    Tangential velocity
 */
double u_Werner(double kinematic_viscosity, double wall_distance,
                double friction_velocity);					 


/**
 * @brief Residual function for Werner-Wengle iteration
 *
 * Computes residual: f(u_τ) = u_τ² - g(u, y, ν)
 * where g is derived from the velocity profile inversion.
 *
 * @param[in] kinematic_viscosity Kinematic viscosity
 * @param[in] velocity            Known velocity
 * @param[in] wall_distance       Distance from wall
 * @param[in] friction_velocity   Guess for friction velocity
 * @return    Residual value
 */
double f_Werner(double kinematic_viscosity, double velocity,
                double wall_distance, double friction_velocity);
				
/**
 * @brief Numerical derivative for Werner-Wengle iteration
 */
double df_Werner(double kinematic_viscosity, double velocity,
                 double wall_distance, double friction_velocity);
				 
/**
 * @brief Solves for friction velocity using Werner-Wengle wall function
 *
 * @param[in] kinematic_viscosity Kinematic viscosity
 * @param[in] velocity            Velocity at reference point
 * @param[in] wall_distance       Distance from wall
 * @param[in] initial_guess       Initial guess for u_τ
 * @return    Converged friction velocity
 */
double find_utau_Werner(double kinematic_viscosity, double velocity,
                        double wall_distance, double initial_guess);
						
/**
 * @brief Computes velocity using simple log-law (smooth wall with roughness offset)
 *
 * Simple logarithmic profile:
 *   u = (u_τ / κ) * ln((y + y0) / y0)
 * where y0 is a roughness length scale.
 *
 * @param[in] wall_distance     Distance from wall
 * @param[in] friction_velocity Friction velocity
 * @param[in] roughness_length  Roughness length scale y0
 * @return    Velocity
 *
 * @note This is a simplified model; use u_hydset_roughness for better accuracy
 */
double u_loglaw(double wall_distance, double friction_velocity, double roughness_length);

/**
 * @brief Solves for friction velocity using simple log-law (explicit formula)
 *
 * Explicit inversion: u_τ = κ * u / ln((y + y0) / y0)
 *
 * @param[in] velocity         Known velocity
 * @param[in] wall_distance    Distance from wall
 * @param[in] roughness_length Roughness length scale
 * @return    Friction velocity
 */
double find_utau_loglaw(double velocity, double wall_distance, double roughness_length);

/**
 * @brief Applies standard wall function with Werner-Wengle model
 *
 * This is a high-level interface that:
 *   1. Decomposes velocity into normal/tangential components
 *   2. Applies wall function to tangential velocity
 *   3. Interpolates normal velocity
 *   4. Reconstructs total velocity
 *
 * @param[in]  user              Simulation context
 * @param[in]  distance_reference Distance to reference point
 * @param[in]  distance_boundary  Distance to boundary point
 * @param[in]  velocity_wall      Wall velocity
 * @param[in]  velocity_reference Reference velocity
 * @param[out] velocity_boundary  Output boundary velocity
 * @param[out] friction_velocity  Output friction velocity
 * @param[in]  normal_x           X-component of wall normal
 * @param[in]  normal_y           Y-component of wall normal
 * @param[in]  normal_z           Z-component of wall normal
 */
void wall_function(UserCtx *user, double distance_reference, double distance_boundary,
                   Cmpnts velocity_wall, Cmpnts velocity_reference,
                   Cmpnts *velocity_boundary, PetscReal *friction_velocity,
                   double normal_x, double normal_y, double normal_z);

/**
 * @brief Applies log-law wall function with roughness correction
 *
 * This is the recommended wall function interface for most applications.
 * It uses the roughness-corrected log-law which is accurate for both
 * smooth and rough walls.
 *
 * @param[in]  user              Simulation context
 * @param[in]  roughness_height  Wall roughness height ks
 * @param[in]  distance_reference Distance to reference point
 * @param[in]  distance_boundary  Distance to boundary point
 * @param[in]  velocity_wall      Wall velocity (typically zero)
 * @param[in]  velocity_reference Reference velocity
 * @param[out] velocity_boundary  Output boundary velocity
 * @param[out] friction_velocity  Output friction velocity u_τ
 * @param[in]  normal_x           X-component of wall normal
 * @param[in]  normal_y           Y-component of wall normal
 * @param[in]  normal_z           Z-component of wall normal
 *
 * @note Sets velocity to reference value if y+ > 300 (outside log-layer)
 */
void wall_function_loglaw(UserCtx *user, double roughness_height,
                          double distance_reference, double distance_boundary,
                          Cmpnts velocity_wall, Cmpnts velocity_reference,
                          Cmpnts *velocity_boundary, PetscReal *friction_velocity,
                          double normal_x, double normal_y, double normal_z);
						  
/**
 * @brief Applies Cabot non-equilibrium wall function with pressure gradients
 *
 * This wall function accounts for pressure gradient effects and is most
 * accurate in separated or strongly accelerating/decelerating flows.
 *
 * @param[in]  user               Simulation context
 * @param[in]  roughness_height   Wall roughness (currently unused in this version)
 * @param[in]  distance_reference Distance to reference point
 * @param[in]  distance_boundary  Distance to boundary point
 * @param[in]  velocity_wall      Wall velocity
 * @param[in]  velocity_reference Reference velocity
 * @param[out] velocity_boundary  Output boundary velocity
 * @param[out] friction_velocity  Output friction velocity
 * @param[in]  normal_x           X-component of wall normal
 * @param[in]  normal_y           Y-component of wall normal
 * @param[in]  normal_z           Z-component of wall normal
 * @param[in]  pressure_gradient_x X-component of pressure gradient
 * @param[in]  pressure_gradient_y Y-component of pressure gradient
 * @param[in]  pressure_gradient_z Z-component of pressure gradient
 * @param[in]  iteration_count    Current iteration count
 *
 * @note Only updates friction velocity every 5 iterations for stability
 * @note Normal pressure gradient currently set to zero
 */
void wall_function_Cabot(UserCtx *user, double roughness_height,
                         double distance_reference, double distance_boundary,
                         Cmpnts velocity_wall, Cmpnts velocity_reference,
                         Cmpnts *velocity_boundary, PetscReal *friction_velocity,
                         double normal_x, double normal_y, double normal_z,
                         double pressure_gradient_x, double pressure_gradient_y,
                         double pressure_gradient_z, int iteration_count);						  

#endif // WALLFUNCTION_H
