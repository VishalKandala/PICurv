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
 * IMPLEMENTED MODELS (each a profile and the friction velocity that inverts it, which
 * ApplyWallFunction() turns into the wall shear stress):
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

/** @brief Werner-Wengle power-law coefficient A and exponent B in u+ = A (y+)^B. */
#define WERNER_A 8.3
#define WERNER_B (1.0 / 7.0)

/** @brief Werner-Wengle sublayer edge: where u+ = y+ meets u+ = A (y+)^B, at A^(1/(1-B)).
 *  The commonly quoted 11.81 is this value rounded, and rounding it leaves the two branches
 *  apart by about 2e-5 of u+, so the profile and its inverse jump there. */
#define WERNER_SUBLAYER_YPLUS (pow(WERNER_A, 1.0 / (1.0 - WERNER_B)))

/** @brief Smooth-to-rough transition y+ threshold */
#define ROUGHNESS_TRANSITION_YPLUS 2.25

/** @brief Fully rough regime y+ threshold */
#define FULLY_ROUGH_YPLUS 90.0

/** @brief Eddy viscosity damping coefficient (van Driest damping) */
#define DAMPING_COEFFICIENT 19.0

// ============================================================================
//                   BASIC BOUNDARY CONDITION FUNCTIONS
// ============================================================================

// ============================================================================
//              ROUGHNESS-CORRECTED LOG-LAW FUNCTIONS
// ============================================================================

/**
 * @brief Implementation of \ref E_coeff().
 * @see E_coeff()
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
 * @brief The log law's sublayer edge for a given E: the y+ where u+ = y+ meets
 *        u+ = ln(E y+)/kappa.
 *
 * Switching at a fixed 11.81 instead left the two branches apart by about 2.4% of u+ for a
 * smooth wall (E = exp(kappa B)), so the profile jumped and the friction velocity it implies
 * had two roots in a band of U y/nu. Switching where they meet keeps the profile continuous
 * and increasing, so each reference speed has one friction velocity. Newton from 11.81
 * converges in a few steps for every E the roughness function produces.
 */
static inline double LogLawSublayerEdge(double E)
{
    double y = 11.81;

    for (int iteration = 0; iteration < 50; iteration++) {
        const double g = y - log(E * y) / KAPPA;
        const double dg = 1.0 - 1.0 / (KAPPA * y);
        const double step = g / dg;

        y -= step;
        if (fabs(step) < 1.e-12 * y) break;
    }
    return y;
}

/**
 * @brief Implementation of \ref u_hydset_roughness().
 * @see u_hydset_roughness()
 */
double u_hydset_roughness(double kinematic_viscosity, double wall_distance,
                          double friction_velocity, double roughness_height)
{
    // Normalized wall distance
    double yplus = friction_velocity * wall_distance / kinematic_viscosity;
    const double E = E_coeff(friction_velocity, roughness_height, kinematic_viscosity);
    
    double tangential_velocity;
    
    /* The forward form of exactly the relation find_utau_hydset() inverts, at every y+:
       the logarithmic branch continues above y+ = 300 as it does in the inverse. Whether
       the first cell lies where the law describes the flow is the range check in
       LogWallModelDiagnostics(), not a property of the profile. */
    if (yplus <= LogLawSublayerEdge(E)) {
        // Viscous sublayer: linear velocity profile
        tangential_velocity = friction_velocity * yplus;
    }
    else {
        // Log-layer: logarithmic velocity profile with roughness correction
        tangential_velocity = (friction_velocity / KAPPA) * log(E * yplus);
    }
    
    return tangential_velocity;
}

// ============================================================================
//         NEWTON-RAPHSON SOLVER FOR FRICTION VELOCITY (LOG-LAW)
// ============================================================================

/**
 * @brief Implementation of \ref f_hydset().
 * @see f_hydset()
 */
double f_hydset(double kinematic_viscosity, double known_velocity, 
                double wall_distance, double friction_velocity_guess, 
                double roughness_height)
{
    double yplus = friction_velocity_guess * wall_distance / kinematic_viscosity;
    double E = E_coeff(friction_velocity_guess, roughness_height, kinematic_viscosity);
    double residual;
    
    if (yplus <= LogLawSublayerEdge(E)) {
        // Viscous sublayer: u = u_τ * y+
        residual = friction_velocity_guess * yplus - known_velocity;
    }
    else {
        // Log-layer: u = (u_τ / κ) * ln(E * y+)
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
        if (derivative == 0.0 || !isfinite(derivative)) break;
        
        friction_velocity = friction_velocity_old - residual / derivative;
        /* The friction velocity is positive. The linear branch's residual u_tau^2 y/nu - U
           has a mirror root at -u_tau, which an undamped step from a small guess can reach
           (a rough wall made it common); halving toward the last positive iterate keeps the
           iteration on the physical root. */
        if (!(friction_velocity > 0.0)) friction_velocity = 0.5 * friction_velocity_old;
        
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
 * @brief Implementation of \ref nu_t().
 * @see nu_t()
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
    /* Composite Simpson in s, with y = wall_distance * s^2. The integrand falls off like
       1/y+ across the log layer, so uniform points in y left the near-wall rise under-sampled:
       30 of them were about 33% off at y+ = 1000. Clustering them toward the wall as s^2 does
       makes the mapped integrand smooth at every y+ the model is allowed to reach. */
    const int    intervals = 128;                 /* even */
    const double ds = 1.0 / intervals;
    double       integral_sum = 0.0;

    for (int i = 0; i <= intervals; i++) {
        const double s_value = i * ds;
        const double current_distance = wall_distance * s_value * s_value;
        const double jacobian = 2.0 * wall_distance * s_value;      /* dy/ds */
        const double yplus = current_distance * friction_velocity / kinematic_viscosity;
        const double denominator = (1.0 + nu_t(fabs(yplus))) * kinematic_viscosity;
        const double integrand = ((integration_mode == 0) ? 1.0 : current_distance) / denominator * jacobian;
        const double weight = (i == 0 || i == intervals) ? 1.0 : ((i % 2) ? 4.0 : 2.0);

        integral_sum += weight * integrand;
    }
    return integral_sum * ds / 3.0;
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
 * @brief Implementation of \ref f_Cabot().
 * @see f_Cabot()
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
        /* Per wall cell and step, so reported only the first time on each rank. */
        static PetscBool reported = PETSC_FALSE;
        if (!reported) {
            reported = PETSC_TRUE;
            LOG_ALLOW(LOCAL, LOG_WARNING,
                      "Cabot wall-function iteration reached its limit without the requested tolerance "
                      "(reported once per rank; later occurrences are not).\n");
        }
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
 * @brief Implementation of \ref u_Werner().
 * @see u_Werner()
 */
double u_Werner(double kinematic_viscosity, double wall_distance,
                double friction_velocity)
{
    /* The pointwise profile. ApplyWallFunction() takes the friction velocity from its exact
       inverse, find_utau_Werner(), at the reference cell. */
    double yplus = friction_velocity * wall_distance / kinematic_viscosity;
    
    // Werner-Wengle constants
    const double power_law_coefficient = WERNER_A;
    const double power_law_exponent = WERNER_B;
    
    double velocity;
    
    if (yplus <= WERNER_SUBLAYER_YPLUS) {
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
    PROFILE_HELPER_BEGIN;
    const double A = WERNER_A;
    const double B = WERNER_B;
    const double speed = fabs(velocity);
    const double nu_over_y = kinematic_viscosity / wall_distance;

    /* The switch is written in terms of the geometry alone. Expressing it through the
       friction velocity, as the iterative residual does, makes the branch depend on the
       answer and forces a Newton solve the closed form does not need. */
    const double transition_speed = 0.5 * nu_over_y * pow(A, 2.0 / (1.0 - B));

    if (speed <= transition_speed) {
        /* Integrating u+ = y+ across the cell leaves the factor of two that separates
           this from the pointwise relation tau_w = nu u / y. */
        PROFILE_HELPER_END;
        return 2.0 * kinematic_viscosity * speed / wall_distance;
    }

    {
        const double term1 = 0.5 * (1.0 - B) * pow(A, (1.0 + B) / (1.0 - B)) *
                             pow(nu_over_y, 1.0 + B);
        const double term2 = (1.0 + B) / A * pow(nu_over_y, B) * speed;

        PROFILE_HELPER_END;
        return pow(term1 + term2, 2.0 / (1.0 + B));
    }
    PROFILE_HELPER_END;
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
    PROFILE_HELPER_BEGIN;
    const double A = WERNER_A;
    const double B = WERNER_B;
    const double stress = fabs(wall_shear_stress);
    const double nu_over_y = kinematic_viscosity / wall_distance;

    /* The same threshold as taw_Werner(), carried over to the stress side. The two
       branches meet here exactly, so selecting on the stress cannot land on the wrong
       side of the velocity switch. */
    const double transition_stress = nu_over_y * nu_over_y * pow(A, 2.0 / (1.0 - B));

    if (stress <= transition_stress) {
        PROFILE_HELPER_END;
        return stress * wall_distance / (2.0 * kinematic_viscosity);
    }

    {
        const double offset = 0.5 * (1.0 - B) * pow(A, (1.0 + B) / (1.0 - B)) *
                              pow(nu_over_y, 1.0 + B);
        /* Non-negative by construction above the threshold - it is zero only in the
           limit of the switch - so the max guards roundoff, not a modelling case. */
        const double bracket = PetscMax(pow(stress, 0.5 * (1.0 + B)) - offset, 0.0);

        PROFILE_HELPER_END;
        return A / (1.0 + B) * pow(1.0 / nu_over_y, B) * bracket;
    }
    PROFILE_HELPER_END;
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
    const double A = WERNER_A;
    const double B = WERNER_B;
    const double speed = fabs(velocity);
    const double y_over_nu = wall_distance / kinematic_viscosity;
    (void)initial_guess;

    /* The exact inverse of `u_Werner()`, branch for branch. Along that profile
       U y / nu = (y+)^2 in the sublayer, so the branch it was evaluated on is decided by
       the data alone: the sublayer holds while U y / nu <= WERNER_SUBLAYER_YPLUS^2, the
       same switch `u_Werner()` applies to y+. Both branches then invert explicitly. */
    if (speed * y_over_nu <= WERNER_SUBLAYER_YPLUS * WERNER_SUBLAYER_YPLUS) {
        return sqrt(speed / y_over_nu);                           /* U = u_tau^2 y / nu */
    }
    return pow(speed / (A * pow(y_over_nu, B)), 1.0 / (1.0 + B)); /* U = A u_tau^(1+B) (y/nu)^B */
}

// ============================================================================
//                  UTILITY FUNCTIONS
// ============================================================================

// ============================================================================
//          HIGH-LEVEL WALL FUNCTION INTERFACE FUNCTIONS
// ============================================================================

