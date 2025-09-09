#pragma once

#include <cmath>
#include <vector>
#include <array>

namespace lbm {
namespace utils {

/**
 * @brief Mathematical utility functions for LBM simulations
 */
class MathUtils {
public:
    /**
     * @brief Calculate equilibrium distribution function for D2Q9 lattice
     * @param direction Velocity direction index (0-8)
     * @param density Local fluid density
     * @param velocity_x Velocity component in x-direction
     * @param velocity_y Velocity component in y-direction
     * @return Equilibrium distribution value
     */
    static double calculateEquilibrium(int direction, double density, 
                                     double velocity_x, double velocity_y);

    /**
     * @brief Calculate macroscopic density from distribution functions
     * @param distributions Array of distribution functions [Q]
     * @return Macroscopic density
     */
    static double calculateDensity(const double* distributions);

    /**
     * @brief Calculate macroscopic velocity from distribution functions
     * @param distributions Array of distribution functions [Q]
     * @param density Macroscopic density
     * @return Pair of (velocity_x, velocity_y)
     */
    static std::pair<double, double> calculateVelocity(const double* distributions, double density);

    /**
     * @brief Calculate gradient using isotropic finite differences
     * @param field 2D field data [Nx][Ny]
     * @param Nx Grid size in x-direction
     * @param Ny Grid size in y-direction
     * @param x Grid position x
     * @param y Grid position y
     * @return Pair of (gradient_x, gradient_y)
     */
    static std::pair<double, double> calculateGradient(double** field, int Nx, int Ny, int x, int y);

    /**
     * @brief Calculate Laplacian using isotropic finite differences
     * @param field 2D field data [Nx][Ny]
     * @param Nx Grid size in x-direction
     * @param Ny Grid size in y-direction
     * @param x Grid position x
     * @param y Grid position y
     * @return Laplacian value
     */
    static double calculateLaplacian(double** field, int Nx, int Ny, int x, int y);

    /**
     * @brief Clamp value to specified range
     * @param value Value to clamp
     * @param min_val Minimum allowed value
     * @param max_val Maximum allowed value
     * @return Clamped value
     */
    template<typename T>
    static T clamp(T value, T min_val, T max_val) {
        return std::max(min_val, std::min(max_val, value));
    }

    /**
     * @brief Check if value is finite (not NaN or infinity)
     * @param value Value to check
     * @return true if finite, false otherwise
     */
    static bool isFinite(double value) {
        return std::isfinite(value);
    }

    /**
     * @brief Safe division with fallback value
     * @param numerator Numerator
     * @param denominator Denominator
     * @param fallback_value Value to return if denominator is too small
     * @param epsilon Minimum allowed denominator magnitude
     * @return Division result or fallback value
     */
    static double safeDivision(double numerator, double denominator, 
                              double fallback_value = 0.0, double epsilon = 1e-12);

    /**
     * @brief Calculate analytical Poiseuille flow velocity profile
     * @param y Grid position in y-direction
     * @param channel_height Channel height
     * @param body_force Applied body force
     * @param viscosity Kinematic viscosity
     * @return Analytical velocity
     */
    static double analyticalPoiseuilleVelocity(int y, double channel_height, 
                                              double body_force, double viscosity);

    /**
     * @brief Calculate L2 norm of the difference between two fields
     * @param field1 First field [Nx*Ny]
     * @param field2 Second field [Nx*Ny]
     * @param size Total number of grid points
     * @return L2 norm of difference
     */
    static double calculateL2Norm(const double* field1, const double* field2, int size);

    /**
     * @brief Calculate maximum absolute difference between two fields
     * @param field1 First field [Nx*Ny]
     * @param field2 Second field [Nx*Ny]
     * @param size Total number of grid points
     * @return Maximum absolute difference
     */
    static double calculateMaxDifference(const double* field1, const double* field2, int size);

private:
    // LBM constants
    static constexpr int Q = 9;
    static constexpr double cs2 = 1.0 / 3.0;
    static constexpr int ex[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    static constexpr int ey[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
    static constexpr double w[Q] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 
                                   1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
};

} // namespace utils
} // namespace lbm