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
 * @brief Internal helper implementation: `noslip()`.
 * @details Local to this translation unit.
 */
void noslip(UserCtx *user, 
            double distance_reference, double distance_boundary,
            Cmpnts velocity_wall, Cmpnts velocity_reference,
            Cmpnts *velocity_boundary,
            double normal_x, double normal_y, double normal_z)
{
    (void)user;
    (void)normal_x;
    (void)normal_y;
    (void)normal_z;
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
 * @brief Internal helper implementation: `freeslip()`.
 * @details Local to this translation unit.
 */
void freeslip(UserCtx *user,
              double distance_reference, double distance_boundary,
              Cmpnts velocity_wall, Cmpnts velocity_reference,
              Cmpnts *velocity_boundary,
              double normal_x, double normal_y, double normal_z)
{
    (void)user;
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
 * @brief Internal helper implementation: `E_coeff()`.
 * @details Local to this translation unit.
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
 * @brief Internal helper implementation: `u_hydset_roughness()`.
 * @details Local to this translation unit.
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
 * @brief Internal helper implementation: `f_hydset()`.
 * @details Local to this translation unit.
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
 * @brief Implementation of \ref df_hydset().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/wallfunction.h`.
 * @see df_hydset()
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
 * @brief Implementation of \ref find_utau_hydset().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/wallfunction.h`.
 * @see find_utau_hydset()
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
 * @brief Internal helper implementation: `nu_t()`.
 * @details Local to this translation unit.
 */
double nu_t(double yplus)
{
    // van Driest damping function: [1 - exp(-y+ / A+)]²
    double damping_function = 1.0 - exp(-yplus / DAMPING_COEFFICIENT);
    
    // Mixing length model: ν_t / ν = κ * y+ * D²
    return KAPPA * yplus * damping_function * damping_function;
}

/**
 * @brief Implementation of \ref integrate_1().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/wallfunction.h`.
 * @see integrate_1()
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
 * @brief Implementation of \ref taw().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/wallfunction.h`.
 * @see taw()
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
 * @brief Implementation of \ref u_Cabot().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/wallfunction.h`.
 * @see u_Cabot()
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
 * @brief Internal helper implementation: `f_Cabot()`.
 * @details Local to this translation unit.
 */
double f_Cabot(double kinematic_viscosity, double velocity, double wall_distance,
               double friction_velocity_guess, double pressure_gradient_tangent,
               double pressure_gradient_normal)
{
    (void)pressure_gradient_normal;
    double wall_shear = taw(kinematic_viscosity, friction_velocity_guess,
                           wall_distance, velocity, pressure_gradient_tangent);
    
    return friction_velocity_guess - sqrt(fabs(wall_shear));
}

/**
 * @brief Implementation of \ref df_Cabot().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/wallfunction.h`.
 * @see df_Cabot()
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
 * @brief Implementation of \ref find_utau_Cabot().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/wallfunction.h`.
 * @see find_utau_Cabot()
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
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "Cabot wall-function iteration reached its limit without the requested tolerance.\n");
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
 * @brief Internal helper implementation: `u_Werner()`.
 * @details Local to this translation unit.
 */
double u_Werner(double kinematic_viscosity, double wall_distance,
                double friction_velocity)
{
    /* The pointwise profile. Not the production path: the correction uses the cell
       integral of this profile, `u_Werner_explicit()`, which is what the Werner-Wengle
       model actually specifies. Retained for reference and unit coverage. */
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

#undef __FUNCT__
#define __FUNCT__ "taw_Werner"
/**
 * @brief Implementation of \ref taw_Werner().
 * @details Full API contract is documented with the header declaration in
 *          `include/wallfunction.h`.
 */
double taw_Werner(double kinematic_viscosity, double velocity, double wall_distance)
{
    const double A = 8.3;
    const double B = 1.0 / 7.0;
    const double speed = fabs(velocity);
    const double nu_over_y = kinematic_viscosity / wall_distance;

    /* The switch is written in terms of the geometry alone. Expressing it through the
       friction velocity, as the iterative residual does, makes the branch depend on the
       answer and forces a Newton solve the closed form does not need. */
    const double transition_speed = 0.5 * nu_over_y * pow(A, 2.0 / (1.0 - B));

    if (speed <= transition_speed) {
        /* Integrating u+ = y+ across the cell leaves the factor of two that separates
           this from the pointwise relation tau_w = nu u / y. */
        return 2.0 * kinematic_viscosity * speed / wall_distance;
    }

    {
        const double term1 = 0.5 * (1.0 - B) * pow(A, (1.0 + B) / (1.0 - B)) *
                             pow(nu_over_y, 1.0 + B);
        const double term2 = (1.0 + B) / A * pow(nu_over_y, B) * speed;

        return pow(term1 + term2, 2.0 / (1.0 + B));
    }
}

#undef __FUNCT__
#define __FUNCT__ "u_Werner_explicit"
/**
 * @brief Implementation of \ref u_Werner_explicit().
 * @details Full API contract is documented with the header declaration in
 *          `include/wallfunction.h`.
 */
double u_Werner_explicit(double kinematic_viscosity, double wall_distance,
                         double wall_shear_stress)
{
    const double A = 8.3;
    const double B = 1.0 / 7.0;
    const double stress = fabs(wall_shear_stress);
    const double nu_over_y = kinematic_viscosity / wall_distance;

    /* The same threshold as taw_Werner(), carried over to the stress side. The two
       branches meet here exactly, so selecting on the stress cannot land on the wrong
       side of the velocity switch. */
    const double transition_stress = nu_over_y * nu_over_y * pow(A, 2.0 / (1.0 - B));

    if (stress <= transition_stress) {
        return stress * wall_distance / (2.0 * kinematic_viscosity);
    }

    {
        const double offset = 0.5 * (1.0 - B) * pow(A, (1.0 + B) / (1.0 - B)) *
                              pow(nu_over_y, 1.0 + B);
        /* Non-negative by construction above the threshold - it is zero only in the
           limit of the switch - so the max guards roundoff, not a modelling case. */
        const double bracket = PetscMax(pow(stress, 0.5 * (1.0 + B)) - offset, 0.0);

        return A / (1.0 + B) * pow(1.0 / nu_over_y, B) * bracket;
    }
}

/**
 * @brief Internal helper implementation: `f_Werner()`.
 * @details Local to this translation unit.
 *
 *          Not wired into the solver. Retained with `df_Werner()` and
 *          `find_utau_Werner()` for reference and for their unit coverage; the
 *          production path is the closed-form `taw_Werner()` / `u_Werner_explicit()`
 *          pair. Two things keep this residual out of it. Its branch threshold is
 *          written through the friction velocity it is solving for, so the branch can
 *          flip mid-iteration and a Newton solve becomes necessary where the model's
 *          own inversion needs none. And its branches invert different relations - the
 *          sublayer branch is pointwise (tau_w = nu u / y), the power branch is the cell
 *          integral - so it reconstructs `u_Werner()` only below y+ = 11.81, which is
 *          the range its unit test happens to sample.
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
 * @brief Implementation of \ref df_Werner().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/wallfunction.h`.
 * @see df_Werner()
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
 * @brief Implementation of \ref find_utau_Werner().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/wallfunction.h`.
 * @see find_utau_Werner()
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
        LOG_ALLOW(GLOBAL, LOG_WARNING,
                  "Werner-Wengle wall-function iteration reached its limit without the requested tolerance.\n");
    }
    
    return new_guess;
}

// ============================================================================
//              SIMPLE LOG-LAW WALL FUNCTION (NO ROUGHNESS)
// ============================================================================

/**
 * @brief Implementation of \ref u_loglaw().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/wallfunction.h`.
 * @see u_loglaw()
 */
double u_loglaw(double wall_distance, double friction_velocity, double roughness_length)
{
    return friction_velocity * (1.0 / KAPPA) * log((roughness_length + wall_distance) / roughness_length);
}

/**
 * @brief Internal helper implementation: `find_utau_loglaw()`.
 * @details Local to this translation unit.
 */
double find_utau_loglaw(double velocity, double wall_distance, double roughness_length)
{
    return KAPPA * velocity / log((wall_distance + roughness_length) / roughness_length);
}

// ============================================================================
//                  UTILITY FUNCTIONS
// ============================================================================

// ============================================================================
//          HIGH-LEVEL WALL FUNCTION INTERFACE FUNCTIONS
// ============================================================================

/**
 * @brief Implementation of \ref wall_function().
 * @details Full API contract (arguments, ownership, side effects) is documented with
 *          the header declaration in `include/wallfunction.h`.
 * @see wall_function()
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
    
    /* Werner-Wengle inverts in closed form, so the wall stress comes straight from the
       reference cell's speed with no inner iteration, and the boundary cell's velocity
       comes back through the exact inverse of that same relation. Using two relations
       here - the cell integral one way and the pointwise profile the other - is what
       previously drove the corrected cell faster than the reference cell it sits inside.
       Same two distances, in the same roles, as the log-law and Cabot paths. */
    double wall_shear_stress = taw_Werner(kinematic_viscosity, tangential_magnitude,
                                          distance_reference);
    double tangential_modeled = u_Werner_explicit(kinematic_viscosity, distance_boundary,
                                                  wall_shear_stress);
    if (friction_velocity) {
        *friction_velocity = (PetscReal)sqrt(wall_shear_stress);
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
 * @brief Internal helper implementation: `wall_function_loglaw()`.
 * @details Local to this translation unit.
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
 * @brief Internal helper implementation: `wall_function_Cabot()`.
 * @details Local to this translation unit.
 */
void wall_function_Cabot(UserCtx *user, double roughness_height,
                         double distance_reference, double distance_boundary,
                         Cmpnts velocity_wall, Cmpnts velocity_reference,
                         Cmpnts *velocity_boundary, PetscReal *friction_velocity,
                         double normal_x, double normal_y, double normal_z,
                         double pressure_gradient_x, double pressure_gradient_y,
                         double pressure_gradient_z, int iteration_count)
{
    (void)roughness_height;
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
