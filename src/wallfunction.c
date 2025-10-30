/**
 * @file wallfunction.c
 * @brief Wall function implementations for near-wall turbulence modeling
 *
 * This file contains various wall function models that bridge the gap between
 * the wall and the first computational grid point. Wall functions allow simulations
 * to avoid resolving the viscous sublayer, reducing computational cost while
 * maintaining reasonable accuracy for turbulent wall-bounded flows.
 *
 * PHYSICAL BACKGROUND:
 * The near-wall region in turbulent flows is characterized by three layers:
 *   1. Viscous sublayer (y+ < 5): Dominated by viscous effects, u+ = y+
 *   2. Buffer layer (5 < y+ < 30): Transition region
 *   3. Log-law region (y+ > 30): Inertial layer, u+ = (1/κ) ln(y+) + B
 * where:
 *   u+ = u / u_τ (normalized velocity)
 *   y+ = y * u_τ / ν (normalized wall distance)
 *   u_τ = sqrt(τ_w / ρ) (friction velocity)
 *   κ ≈ 0.41 (von Karman constant)
 *   B ≈ 5.5 (log-law intercept for smooth walls)
 *
 * IMPLEMENTED MODELS:
 *   - No-slip: Linear interpolation for stationary walls
 *   - Free-slip: Normal interpolation, tangential extrapolation
 *   - Werner-Wengle: Algebraic wall function
 *   - Log-law: Standard equilibrium wall function with roughness
 *   - Cabot: Non-equilibrium wall function with pressure gradient effects
 *
 * @author Original implementation adapted from legacy code
 * @date Enhanced with comprehensive documentation
 */

#include "wallfunction.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// ============================================================================
//                        PHYSICAL CONSTANTS
// ============================================================================

/** @brief von Karman constant (universal turbulence constant) */
#define KAPPA 0.41

/** @brief Log-law intercept constant B for smooth walls */
#define LOGLAW_B 5.5

/** @brief Viscous sublayer thickness y+ threshold */
#define VISCOUS_SUBLAYER_YPLUS 11.81

/** @brief Smooth-to-rough transition y+ threshold */
#define ROUGHNESS_TRANSITION_YPLUS 2.25

/** @brief Fully rough regime y+ threshold */
#define FULLY_ROUGH_YPLUS 90.0

/** @brief Eddy viscosity damping coefficient (van Driest damping) */
#define DAMPING_COEFFICIENT 19.0

// ============================================================================
//                   BASIC BOUNDARY CONDITION FUNCTIONS
// ============================================================================

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
            double normal_x, double normal_y, double normal_z)
{
    // Compute velocity difference between reference point and wall
    double delta_u = velocity_reference.x - velocity_wall.x;
    double delta_v = velocity_reference.y - velocity_wall.y;
    double delta_w = velocity_reference.z - velocity_wall.z;
    
    // Linear interpolation factor
    double interpolation_factor = distance_boundary / distance_reference;
    
    // Apply linear interpolation
    (*velocity_boundary).x = interpolation_factor * delta_u;
    (*velocity_boundary).y = interpolation_factor * delta_v;
    (*velocity_boundary).z = interpolation_factor * delta_w;
    
    // Add wall velocity (shift to absolute frame)
    (*velocity_boundary).x += velocity_wall.x;
    (*velocity_boundary).y += velocity_wall.y;
    (*velocity_boundary).z += velocity_wall.z;
}

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
              double normal_x, double normal_y, double normal_z)
{
    // Extract normal components of velocities
    double wall_normal_velocity = velocity_wall.x * normal_x + 
                                   velocity_wall.y * normal_y + 
                                   velocity_wall.z * normal_z;
    
    double reference_normal_velocity = velocity_reference.x * normal_x + 
                                        velocity_reference.y * normal_y + 
                                        velocity_reference.z * normal_z;
    
    // Interpolate normal velocity component
    double boundary_normal_velocity = wall_normal_velocity + 
        (reference_normal_velocity - wall_normal_velocity) * (distance_boundary / distance_reference);
    
    // Extract tangential velocity components (extrapolated from reference point)
    double tangential_u = velocity_reference.x - reference_normal_velocity * normal_x;
    double tangential_v = velocity_reference.y - reference_normal_velocity * normal_y;
    double tangential_w = velocity_reference.z - reference_normal_velocity * normal_z;
    
    // Reconstruct total velocity: U = U_t + U_n * n
    (*velocity_boundary).x = tangential_u + boundary_normal_velocity * normal_x;
    (*velocity_boundary).y = tangential_v + boundary_normal_velocity * normal_y;
    (*velocity_boundary).z = tangential_w + boundary_normal_velocity * normal_z;
}

// ============================================================================
//              ROUGHNESS-CORRECTED LOG-LAW FUNCTIONS
// ============================================================================

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
double E_coeff(double friction_velocity, double roughness_height, double kinematic_viscosity)
{
    // Compute roughness Reynolds number
    double roughness_reynolds = friction_velocity * roughness_height / kinematic_viscosity;
    
    double roughness_correction;
    
    if (roughness_reynolds <= ROUGHNESS_TRANSITION_YPLUS) {
        // Hydraulically smooth regime: no correction needed
        roughness_correction = 0.0;
    }
    else if (roughness_reynolds < FULLY_ROUGH_YPLUS) {
        // Transitional regime: smooth interpolation using sine function
        // This provides a gradual transition between smooth and rough wall behavior
        double transition_factor = log(roughness_reynolds) - log(ROUGHNESS_TRANSITION_YPLUS);
        double max_correction = LOGLAW_B - 8.5 + (1.0 / KAPPA) * log(roughness_reynolds);
        
        roughness_correction = max_correction * 
            sin(0.4258 * (log(roughness_reynolds) - 0.811));
    }
    else {
        // Fully rough regime: maximum correction
        roughness_correction = LOGLAW_B - 8.5 + (1.0 / KAPPA) * log(roughness_reynolds);
    }
    
    // Return E = exp[κ(B - ΔB)]
    return exp(KAPPA * (LOGLAW_B - roughness_correction));
}

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
                          double friction_velocity, double roughness_height)
{
    // Normalized wall distance
    double yplus = friction_velocity * wall_distance / kinematic_viscosity;
    
    // Threshold values for regime transitions
    const double viscous_limit = VISCOUS_SUBLAYER_YPLUS;
    const double loglayer_limit = 300.0;
    
    double tangential_velocity;
    
    if (yplus <= viscous_limit) {
        // Viscous sublayer: linear velocity profile
        tangential_velocity = friction_velocity * yplus;
    }
    else if (yplus <= loglayer_limit) {
        // Log-layer: logarithmic velocity profile with roughness correction
        double E = E_coeff(friction_velocity, roughness_height, kinematic_viscosity);
        tangential_velocity = (friction_velocity / KAPPA) * log(E * yplus);
    }
    else {
        // Outside valid range (wake region or beyond)
        tangential_velocity = -1.0;
    }
    
    return tangential_velocity;
}

// ============================================================================
//         NEWTON-RAPHSON SOLVER FOR FRICTION VELOCITY (LOG-LAW)
// ============================================================================

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
                double roughness_height)
{
    double yplus = friction_velocity_guess * wall_distance / kinematic_viscosity;
    double residual;
    
    if (yplus <= VISCOUS_SUBLAYER_YPLUS) {
        // Viscous sublayer: u = u_τ * y+
        residual = friction_velocity_guess * yplus - known_velocity;
    }
    else {
        // Log-layer: u = (u_τ / κ) * ln(E * y+)
        double E = E_coeff(friction_velocity_guess, roughness_height, kinematic_viscosity);
        residual = friction_velocity_guess * (1.0 / KAPPA * log(E * yplus)) - known_velocity;
    }
    
    return residual;
}

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
                 double roughness_height)
{
    const double perturbation = 1.e-7;
    
    double f_plus = f_hydset(kinematic_viscosity, known_velocity, wall_distance,
                            friction_velocity_guess + perturbation, roughness_height);
    
    double f_minus = f_hydset(kinematic_viscosity, known_velocity, wall_distance,
                             friction_velocity_guess - perturbation, roughness_height);
    
    return (f_plus - f_minus) / (2.0 * perturbation);
}

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
                        double roughness_height)
{
    double friction_velocity = initial_guess;
    double friction_velocity_old = initial_guess;
    
    const int max_iterations = 30;
    const double convergence_tolerance = 1.e-7;
    
    int iteration;
    for (iteration = 0; iteration < max_iterations; iteration++) {
        // Newton-Raphson update
        double residual = f_hydset(kinematic_viscosity, known_velocity, wall_distance,
                                  friction_velocity_old, roughness_height);
        double derivative = df_hydset(kinematic_viscosity, known_velocity, wall_distance,
                                     friction_velocity_old, roughness_height);
        
        friction_velocity = friction_velocity_old - residual / derivative;
        
        // Check convergence
        if (fabs(friction_velocity - friction_velocity_old) < convergence_tolerance) {
            break;
        }
        
        friction_velocity_old = friction_velocity;
    }
    
    // Warn if convergence not achieved
    if (iteration == max_iterations) {
        LOG_ALLOW(LOCAL, LOG_DEBUG, 
                 "WARNING: u_tau iteration did not converge after %d iterations. "
                 "Final difference: %le\n", 
                 iteration, fabs(friction_velocity - friction_velocity_old));
    }
    
    return friction_velocity;
}

// ============================================================================
//               TURBULENT EDDY VISCOSITY MODELS
// ============================================================================

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
double nu_t(double yplus)
{
    // van Driest damping function: [1 - exp(-y+ / A+)]²
    double damping_function = 1.0 - exp(-yplus / DAMPING_COEFFICIENT);
    
    // Mixing length model: ν_t / ν = κ * y+ * D²
    return KAPPA * yplus * damping_function * damping_function;
}

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
                   double friction_velocity, int integration_mode)
{
    const int num_points = 30;
    double step_size = wall_distance / (num_points - 1);
    double integral_values[num_points];
    
    // Compute integrand at each point
    for (int i = 0; i < num_points; i++) {
        double current_distance = i * step_size;
        double yplus = current_distance * friction_velocity / kinematic_viscosity;
        double eddy_viscosity_ratio = nu_t(yplus);
        
        if (integration_mode == 0) {
            // Mode 0: integrand = 1 / [(1 + ν_t/ν) * ν]
            integral_values[i] = 1.0 / ((1.0 + eddy_viscosity_ratio) * kinematic_viscosity);
        }
        else {
            // Mode 1: integrand = y / [(1 + ν_t/ν) * ν]
            integral_values[i] = current_distance / 
                                ((1.0 + eddy_viscosity_ratio) * kinematic_viscosity);
        }
    }
    
    // Trapezoidal rule integration with corrected endpoints
    double integral_sum = 0.0;
    
    // Interior points: weight = (f[i-1] + 2*f[i] + f[i+1]) * dy/4
    for (int i = 1; i < num_points - 1; i++) {
        integral_sum += (integral_values[i-1] + 2.0 * integral_values[i] + 
                        integral_values[i+1]) * 0.25 * step_size;
    }
    
    // Endpoint correction for higher accuracy
    integral_sum += 0.1667 * step_size * 
                   (2.0 * integral_values[0] + integral_values[1] +
                    2.0 * integral_values[num_points-1] + integral_values[num_points-2]);
    
    return integral_sum;
}

// ============================================================================
//            CABOT NON-EQUILIBRIUM WALL FUNCTION
// ============================================================================

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
           double wall_distance, double velocity, double pressure_gradient_tangent)
{
    // Compute integrated viscosity profiles
    double F1 = integrate_1(kinematic_viscosity, wall_distance, friction_velocity, 0);
    double Fy = integrate_1(kinematic_viscosity, wall_distance, friction_velocity, 1);
    
    // Apply momentum balance: τ_w * F1 = u - (dp/dx) * Fy
    return (1.0 / F1) * (velocity - pressure_gradient_tangent * Fy);
}

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
               double wall_shear_stress)
{
    double F1 = integrate_1(kinematic_viscosity, wall_distance, friction_velocity, 0);
    double Fy = integrate_1(kinematic_viscosity, wall_distance, friction_velocity, 1);
    
    return wall_shear_stress * F1 + pressure_gradient_tangent * Fy;
}

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
               double pressure_gradient_normal)
{
    double wall_shear = taw(kinematic_viscosity, friction_velocity_guess,
                           wall_distance, velocity, pressure_gradient_tangent);
    
    return friction_velocity_guess - sqrt(fabs(wall_shear));
}

/**
 * @brief Numerical derivative for Cabot wall function
 */
double df_Cabot(double kinematic_viscosity, double velocity, double wall_distance,
                double friction_velocity_guess, double pressure_gradient_tangent,
                double pressure_gradient_normal)
{
    const double perturbation = 1.e-7;
    
    double f_plus = f_Cabot(kinematic_viscosity, velocity, wall_distance,
                           friction_velocity_guess + perturbation,
                           pressure_gradient_tangent, pressure_gradient_normal);
    
    double f_minus = f_Cabot(kinematic_viscosity, velocity, wall_distance,
                            friction_velocity_guess - perturbation,
                            pressure_gradient_tangent, pressure_gradient_normal);
    
    return (f_plus - f_minus) / (2.0 * perturbation);
}

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
                     double *wall_shear_velocity, double *wall_shear_normal)
{
    double current_guess = initial_guess;
    double new_guess;
    
    const int max_iterations = 30;
    const double convergence_tolerance = 1.e-10;
    
    int iteration;
    for (iteration = 0; iteration < max_iterations; iteration++) {
        double residual = f_Cabot(kinematic_viscosity, velocity, wall_distance,
                                 current_guess, pressure_gradient_tangent,
                                 pressure_gradient_normal);
        
        double derivative = df_Cabot(kinematic_viscosity, velocity, wall_distance,
                                    current_guess, pressure_gradient_tangent,
                                    pressure_gradient_normal);
        
        new_guess = fabs(current_guess - residual / derivative);
        
        if (fabs(current_guess - new_guess) < convergence_tolerance) {
            break;
        }
        
        current_guess = new_guess;
    }
    
    if (fabs(current_guess - new_guess) > 1.e-5 && iteration >= 29) {
        printf("\n!!!!!!!! Cabot Iteration Failed !!!!!!!!\n");
    }
    
    *friction_velocity = new_guess;
    *wall_shear_velocity = taw(kinematic_viscosity, new_guess, wall_distance,
                              velocity, pressure_gradient_tangent);
    *wall_shear_normal = taw(kinematic_viscosity, new_guess, wall_distance,
                            0.0, pressure_gradient_normal);
}

// ============================================================================
//            WERNER-WENGLE ALGEBRAIC WALL FUNCTION
// ============================================================================

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
                double friction_velocity)
{
    double yplus = friction_velocity * wall_distance / kinematic_viscosity;
    
    // Werner-Wengle constants
    const double power_law_coefficient = 8.3;
    const double power_law_exponent = 1.0 / 7.0;
    
    double velocity;
    
    if (yplus <= VISCOUS_SUBLAYER_YPLUS) {
        // Viscous sublayer: u+ = y+
        velocity = yplus * friction_velocity;
    }
    else {
        // Power-law region: u+ = A * (y+)^B
        velocity = power_law_coefficient * pow(yplus, power_law_exponent) * friction_velocity;
    }
    
    return velocity;
}

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
                double wall_distance, double friction_velocity)
{
    const double power_law_coefficient = 8.3;
    const double power_law_exponent = 1.0 / 7.0;
    
    // Transition point between viscous and power-law regions
    double transition_distance = VISCOUS_SUBLAYER_YPLUS * kinematic_viscosity / friction_velocity;
    
    // Transition velocity
    double transition_velocity = kinematic_viscosity / (2.0 * transition_distance) *
        pow(power_law_coefficient, 2.0 / (1.0 - power_law_exponent));
    
    double residual;
    
    if (fabs(velocity) <= transition_velocity) {
        // Viscous sublayer regime
        residual = friction_velocity * friction_velocity - 
                  velocity / wall_distance * kinematic_viscosity;
    }
    else {
        // Power-law regime (more complex inversion formula)
        double term1 = 0.5 * (1.0 - power_law_exponent) *
                      pow(power_law_coefficient, (1.0 + power_law_exponent) / 
                          (1.0 - power_law_exponent)) *
                      pow(kinematic_viscosity / wall_distance, 1.0 + power_law_exponent);
        
        double term2 = (1.0 + power_law_exponent) / power_law_coefficient *
                      pow(kinematic_viscosity / wall_distance, power_law_exponent) *
                      fabs(velocity);
        
        residual = friction_velocity * friction_velocity -
                  velocity / fabs(velocity) * pow(term1 + term2, 2.0 / (1.0 + power_law_exponent));
    }
    
    return residual;
}

/**
 * @brief Numerical derivative for Werner-Wengle iteration
 */
double df_Werner(double kinematic_viscosity, double velocity,
                 double wall_distance, double friction_velocity)
{
    const double perturbation = 1.e-7;
    
    double f_plus = f_Werner(kinematic_viscosity, velocity, wall_distance,
                            friction_velocity + perturbation);
    
    double f_minus = f_Werner(kinematic_viscosity, velocity, wall_distance,
                             friction_velocity - perturbation);
    
    return (f_plus - f_minus) / (2.0 * perturbation);
}

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
                        double wall_distance, double initial_guess)
{
    double current_guess = initial_guess;
    double new_guess;
    
    const int max_iterations = 20;
    const double convergence_tolerance = 1.e-7;
    
    int iteration;
    for (iteration = 0; iteration < max_iterations; iteration++) {
        double residual = f_Werner(kinematic_viscosity, velocity, wall_distance, current_guess);
        double derivative = df_Werner(kinematic_viscosity, velocity, wall_distance, current_guess);
        
        new_guess = current_guess - residual / derivative;
        
        if (fabs(current_guess - new_guess) < convergence_tolerance) {
            break;
        }
        
        current_guess = new_guess;
    }
    
    if (fabs(current_guess - new_guess) > 1.e-5 && iteration >= 19) {
        printf("\n!!!!!!!! Werner-Wengle Iteration Failed !!!!!!!!\n");
    }
    
    return new_guess;
}

// ============================================================================
//              SIMPLE LOG-LAW WALL FUNCTION (NO ROUGHNESS)
// ============================================================================

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
double u_loglaw(double wall_distance, double friction_velocity, double roughness_length)
{
    return friction_velocity * (1.0 / KAPPA) * log((roughness_length + wall_distance) / roughness_length);
}

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
double find_utau_loglaw(double velocity, double wall_distance, double roughness_length)
{
    return KAPPA * velocity / log((wall_distance + roughness_length) / roughness_length);
}

// ============================================================================
//                  UTILITY FUNCTIONS
// ============================================================================

/**
 * @brief Returns the sign of a number
 *
 * @param[in] value Input value
 * @return    +1 if positive, -1 if negative, 0 if zero
 */
static double sign(double value)
{
    if (value > 0) return 1.0;
    else if (value < 0) return -1.0;
    else return 0.0;
}

// ============================================================================
//          HIGH-LEVEL WALL FUNCTION INTERFACE FUNCTIONS
// ============================================================================

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
                   double normal_x, double normal_y, double normal_z)
{
    SimCtx *simCtx = user->simCtx;
    double kinematic_viscosity = 1.0 / simCtx->ren;
    
    // Decompose velocity into components
    double delta_u = velocity_reference.x - velocity_wall.x;
    double delta_v = velocity_reference.y - velocity_wall.y;
    double delta_w = velocity_reference.z - velocity_wall.z;
    
    double normal_velocity = delta_u * normal_x + delta_v * normal_y + delta_w * normal_z;
    
    double tangential_u = delta_u - normal_velocity * normal_x;
    double tangential_v = delta_v - normal_velocity * normal_y;
    double tangential_w = delta_w - normal_velocity * normal_z;
    
    double tangential_magnitude = sqrt(tangential_u * tangential_u +
                                      tangential_v * tangential_v +
                                      tangential_w * tangential_w);
    
    // Apply Werner-Wengle wall function
    double utau_guess = 0.05;
    double tangential_modeled = u_Werner(kinematic_viscosity, distance_boundary, utau_guess);
    
    // Scale tangential components
    if (tangential_magnitude > 1.e-10) {
        tangential_u *= tangential_modeled / tangential_magnitude;
        tangential_v *= tangential_modeled / tangential_magnitude;
        tangential_w *= tangential_modeled / tangential_magnitude;
    }
    else {
        tangential_u = tangential_v = tangential_w = 0.0;
    }
    
    // Reconstruct velocity: U = U_t + U_n * n
    (*velocity_boundary).x = tangential_u + (distance_boundary / distance_reference) * normal_velocity * normal_x;
    (*velocity_boundary).y = tangential_v + (distance_boundary / distance_reference) * normal_velocity * normal_y;
    (*velocity_boundary).z = tangential_w + (distance_boundary / distance_reference) * normal_velocity * normal_z;
    
    // Add wall velocity
    (*velocity_boundary).x += velocity_wall.x;
    (*velocity_boundary).y += velocity_wall.y;
    (*velocity_boundary).z += velocity_wall.z;
}

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
                          double normal_x, double normal_y, double normal_z)
{
    SimCtx *simCtx = user->simCtx;
    double kinematic_viscosity = 1.0 / simCtx->ren;
    
    // Decompose velocity
    double delta_u = velocity_reference.x - velocity_wall.x;
    double delta_v = velocity_reference.y - velocity_wall.y;
    double delta_w = velocity_reference.z - velocity_wall.z;
    
    double normal_velocity = delta_u * normal_x + delta_v * normal_y + delta_w * normal_z;
    
    double tangential_u = delta_u - normal_velocity * normal_x;
    double tangential_v = delta_v - normal_velocity * normal_y;
    double tangential_w = delta_w - normal_velocity * normal_z;
    
    double tangential_magnitude = sqrt(tangential_u * tangential_u +
                                      tangential_v * tangential_v +
                                      tangential_w * tangential_w);
    
    // Solve for friction velocity
    *friction_velocity = find_utau_hydset(kinematic_viscosity, tangential_magnitude,
                                         distance_reference, 0.001, roughness_height);
    
    // Compute wall function velocity
    double tangential_modeled = u_hydset_roughness(kinematic_viscosity, distance_boundary,
                                                   *friction_velocity, roughness_height);
    
    // Check if outside valid range (y+ > 300)
    if (tangential_modeled < 0.0) {
        tangential_modeled = tangential_magnitude;
    }
    
    // Scale tangential components
    if (tangential_magnitude > 1.e-10) {
        tangential_u *= tangential_modeled / tangential_magnitude;
        tangential_v *= tangential_modeled / tangential_magnitude;
        tangential_w *= tangential_modeled / tangential_magnitude;
    }
    else {
        tangential_u = tangential_v = tangential_w = 0.0;
    }
    
    // Reconstruct total velocity
    (*velocity_boundary).x = tangential_u + (distance_boundary / distance_reference) * normal_velocity * normal_x;
    (*velocity_boundary).y = tangential_v + (distance_boundary / distance_reference) * normal_velocity * normal_y;
    (*velocity_boundary).z = tangential_w + (distance_boundary / distance_reference) * normal_velocity * normal_z;
    
    // Add wall velocity
    (*velocity_boundary).x += velocity_wall.x;
    (*velocity_boundary).y += velocity_wall.y;
    (*velocity_boundary).z += velocity_wall.z;
}

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
                         double pressure_gradient_z, int iteration_count)
{
    SimCtx *simCtx = user->simCtx;
    double kinematic_viscosity = 1.0 / simCtx->ren;
    
    // Decompose velocity
    double delta_u = velocity_reference.x - velocity_wall.x;
    double delta_v = velocity_reference.y - velocity_wall.y;
    double delta_w = velocity_reference.z - velocity_wall.z;
    
    double normal_velocity = delta_u * normal_x + delta_v * normal_y + delta_w * normal_z;
    
    double tangential_u = delta_u - normal_velocity * normal_x;
    double tangential_v = delta_v - normal_velocity * normal_y;
    double tangential_w = delta_w - normal_velocity * normal_z;
    
    double tangential_magnitude = sqrt(tangential_u * tangential_u +
                                      tangential_v * tangential_v +
                                      tangential_w * tangential_w);
    
    // Compute tangential pressure gradient
    double pressure_gradient_tangent = 0.0;
    if (tangential_magnitude > 1.e-10) {
        pressure_gradient_tangent = (pressure_gradient_x * tangential_u +
                                    pressure_gradient_y * tangential_v +
                                    pressure_gradient_z * tangential_w) / tangential_magnitude;
    }
    
    // Normal pressure gradient (currently set to zero)
    double pressure_gradient_normal = 0.0;
    
    // Solve for friction velocity (only on first iteration or periodically)
    if (iteration_count == 0 || iteration_count > 4) {
        double utau, wall_shear1, wall_shear2;
        find_utau_Cabot(kinematic_viscosity, tangential_magnitude, distance_reference,
                       0.01, pressure_gradient_tangent, pressure_gradient_normal,
                       &utau, &wall_shear1, &wall_shear2);
        *friction_velocity = utau;
    }
    
    // Compute wall function velocity with pressure gradient correction
    double tangential_modeled = u_Cabot(kinematic_viscosity, distance_boundary,
                                       *friction_velocity, pressure_gradient_tangent,
                                       taw(kinematic_viscosity, *friction_velocity,
                                           distance_reference, tangential_magnitude,
                                           pressure_gradient_tangent));
    
    // Scale tangential components
    if (tangential_magnitude > 1.e-10) {
        tangential_u *= tangential_modeled / tangential_magnitude;
        tangential_v *= tangential_modeled / tangential_magnitude;
        tangential_w *= tangential_modeled / tangential_magnitude;
    }
    else {
        tangential_u = tangential_v = tangential_w = 0.0;
    }
    
    // Reconstruct total velocity
    (*velocity_boundary).x = tangential_u + (distance_boundary / distance_reference) * normal_velocity * normal_x;
    (*velocity_boundary).y = tangential_v + (distance_boundary / distance_reference) * normal_velocity * normal_y;
    (*velocity_boundary).z = tangential_w + (distance_boundary / distance_reference) * normal_velocity * normal_z;
    
    // Add wall velocity
    (*velocity_boundary).x += velocity_wall.x;
    (*velocity_boundary).y += velocity_wall.y;
    (*velocity_boundary).z += velocity_wall.z;
}