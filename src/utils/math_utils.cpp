#include "../../include/utils/math_utils.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lbm {
namespace utils {

double MathUtils::calculateEquilibrium(int direction, double density, 
                                      double velocity_x, double velocity_y) {
    if (direction < 0 || direction >= Q) {
        throw std::invalid_argument("Invalid direction index");
    }
    
    if (!isFinite(density) || !isFinite(velocity_x) || !isFinite(velocity_y)) {
        throw std::invalid_argument("Non-finite input values");
    }
    
    double eu = ex[direction] * velocity_x + ey[direction] * velocity_y;
    double u2 = velocity_x * velocity_x + velocity_y * velocity_y;
    
    return w[direction] * density * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u2);
}

double MathUtils::calculateDensity(const double* distributions) {
    if (!distributions) {
        throw std::invalid_argument("Null pointer for distributions");
    }
    
    double density = 0.0;
    for (int i = 0; i < Q; i++) {
        if (!isFinite(distributions[i])) {
            throw std::runtime_error("Non-finite distribution function detected");
        }
        density += distributions[i];
    }
    
    return density;
}

std::pair<double, double> MathUtils::calculateVelocity(const double* distributions, double density) {
    if (!distributions) {
        throw std::invalid_argument("Null pointer for distributions");
    }
    
    if (!isFinite(density) || std::abs(density) < 1e-12) {
        return {0.0, 0.0};  // Return zero velocity for near-zero density
    }
    
    double velocity_x = 0.0, velocity_y = 0.0;
    
    for (int i = 0; i < Q; i++) {
        if (!isFinite(distributions[i])) {
            throw std::runtime_error("Non-finite distribution function detected");
        }
        velocity_x += distributions[i] * ex[i];
        velocity_y += distributions[i] * ey[i];
    }
    
    return {velocity_x / density, velocity_y / density};
}

std::pair<double, double> MathUtils::calculateGradient(double** field, int Nx, int Ny, int x, int y) {
    if (!field) {
        throw std::invalid_argument("Null pointer for field");
    }
    
    if (x < 0 || x >= Nx || y < 0 || y >= Ny) {
        throw std::out_of_range("Grid coordinates out of bounds");
    }
    
    double grad_x = 0.0, grad_y = 0.0;
    
    // Use isotropic finite differences with D2Q9 stencil
    for (int i = 1; i < Q; i++) {  // Skip i=0 (rest particle)
        int x_nb = (x + ex[i] + Nx) % Nx;  // Periodic in x
        int y_nb = y + ey[i];
        
        if (y_nb >= 0 && y_nb < Ny) {  // Check y bounds
            if (!isFinite(field[x_nb][y_nb])) {
                throw std::runtime_error("Non-finite field value detected");
            }
            grad_x += w[i] * ex[i] * field[x_nb][y_nb] / cs2;
            grad_y += w[i] * ey[i] * field[x_nb][y_nb] / cs2;
        }
    }
    
    return {grad_x, grad_y};
}

double MathUtils::calculateLaplacian(double** field, int Nx, int Ny, int x, int y) {
    if (!field) {
        throw std::invalid_argument("Null pointer for field");
    }
    
    if (x < 0 || x >= Nx || y < 0 || y >= Ny) {
        throw std::out_of_range("Grid coordinates out of bounds");
    }
    
    double laplacian = 0.0;
    
    // Use isotropic finite differences with D2Q9 stencil
    for (int i = 1; i < Q; i++) {  // Skip i=0 (rest particle)
        int x_nb = (x + ex[i] + Nx) % Nx;  // Periodic in x
        int y_nb = y + ey[i];
        
        if (y_nb >= 0 && y_nb < Ny) {  // Check y bounds
            if (!isFinite(field[x_nb][y_nb]) || !isFinite(field[x][y])) {
                throw std::runtime_error("Non-finite field value detected");
            }
            laplacian += 2.0 * w[i] * (field[x_nb][y_nb] - field[x][y]) / cs2;
        }
    }
    
    return laplacian;
}

double MathUtils::safeDivision(double numerator, double denominator, 
                              double fallback_value, double epsilon) {
    if (!isFinite(numerator) || !isFinite(denominator)) {
        return fallback_value;
    }
    
    if (std::abs(denominator) < epsilon) {
        return fallback_value;
    }
    
    double result = numerator / denominator;
    return isFinite(result) ? result : fallback_value;
}

double MathUtils::analyticalPoiseuilleVelocity(int y, double channel_height, 
                                              double body_force, double viscosity) {
    if (channel_height <= 0 || viscosity <= 0) {
        throw std::invalid_argument("Channel height and viscosity must be positive");
    }
    
    double H = channel_height;
    double y_center = y - H / 2.0;
    
    return body_force * H * H / (8.0 * viscosity) * (1.0 - 4.0 * y_center * y_center / (H * H));
}

double MathUtils::calculateL2Norm(const double* field1, const double* field2, int size) {
    if (!field1 || !field2) {
        throw std::invalid_argument("Null pointer for field data");
    }
    
    if (size <= 0) {
        throw std::invalid_argument("Size must be positive");
    }
    
    double sum_squared_diff = 0.0;
    
    for (int i = 0; i < size; i++) {
        if (!isFinite(field1[i]) || !isFinite(field2[i])) {
            throw std::runtime_error("Non-finite field values detected");
        }
        double diff = field1[i] - field2[i];
        sum_squared_diff += diff * diff;
    }
    
    return std::sqrt(sum_squared_diff / size);
}

double MathUtils::calculateMaxDifference(const double* field1, const double* field2, int size) {
    if (!field1 || !field2) {
        throw std::invalid_argument("Null pointer for field data");
    }
    
    if (size <= 0) {
        throw std::invalid_argument("Size must be positive");
    }
    
    double max_diff = 0.0;
    
    for (int i = 0; i < size; i++) {
        if (!isFinite(field1[i]) || !isFinite(field2[i])) {
            throw std::runtime_error("Non-finite field values detected");
        }
        double diff = std::abs(field1[i] - field2[i]);
        max_diff = std::max(max_diff, diff);
    }
    
    return max_diff;
}

} // namespace utils
} // namespace lbm