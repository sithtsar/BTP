#ifndef INTERFACE_DETECTOR_H
#define INTERFACE_DETECTOR_H

#include "core/lattice.h"
#include "core/fluid_state.h"
#include <array>
#include <cmath>
#include <vector>
#include <algorithm>

namespace elbm {

/**
 * Interface Detection and Analysis for Two-Phase Flows
 *
 * Provides utilities for:
 * - Detecting interface regions between immiscible fluids
 * - Computing interface properties (position, width, normal vectors)
 * - Separating bulk phases from interface for diagnostics
 */
class InterfaceDetector {
public:
    /**
     * Node classification based on density
     */
    enum class Phase {
        LIQUID,    // High density phase
        GAS,       // Low density phase
        INTERFACE, // Transitional region
        INVALID    // Solid or invalid
    };

    /**
     * Interface statistics
     */
    struct InterfaceStats {
        double mean_position_x;  // Mean x-coordinate of interface nodes
        double mean_position_y;  // Mean y-coordinate of interface nodes
        double radius;           // Effective radius (for droplets)
        double width;            // Interface width (10%-90% density transition)
        int num_interface_nodes; // Count of interface nodes
        double spurious_vel_interface; // Maximum velocity at interface
    };

    /**
     * Constructor
     * @param rho_liquid - Density of liquid phase
     * @param rho_gas - Density of gas phase
     * @param interface_threshold - Fraction of density range to define interface (default 0.1 = 10%)
     */
    InterfaceDetector(double rho_liquid, double rho_gas, double interface_threshold = 0.1)
        : rho_liquid_(rho_liquid)
        , rho_gas_(rho_gas)
        , interface_threshold_(interface_threshold) {

        // Compute density thresholds for phase classification
        double rho_range = std::abs(rho_liquid - rho_gas);
        threshold_ = interface_threshold * rho_range;

        // Interface region: gas + threshold < rho < liquid - threshold
        rho_interface_lower_ = std::min(rho_liquid, rho_gas) + threshold_;
        rho_interface_upper_ = std::max(rho_liquid, rho_gas) - threshold_;
    }

    /**
     * Classify a node based on its density
     */
    Phase classify(double rho) const {
        if (rho < rho_interface_lower_) {
            return Phase::GAS;
        } else if (rho > rho_interface_upper_) {
            return Phase::LIQUID;
        } else {
            return Phase::INTERFACE;
        }
    }

    /**
     * Check if a node is at the interface
     */
    bool isInterface(double rho) const {
        return (rho >= rho_interface_lower_) && (rho <= rho_interface_upper_);
    }

    /**
     * Compute density gradient using finite differences
     * @param grid - Lattice grid
     * @param x, y - Coordinates
     * @return Gradient vector [∂ρ/∂x, ∂ρ/∂y]
     */
    template<int D, int Q>
    std::array<double, D> computeDensityGradient(const LatticeGrid<D, Q>& grid, int x, int y) const {
        static_assert(D == 2, "InterfaceDetector currently supports 2D only");

        int nx = grid.nx();
        int ny = grid.ny();

        // Central differences with periodic boundaries
        int xp = (x + 1) % nx;
        int xm = (x - 1 + nx) % nx;
        int yp = (y + 1) % ny;
        int ym = (y - 1 + ny) % ny;

        double grad_x = 0.5 * (grid(xp, y).fluid.rho - grid(xm, y).fluid.rho);
        double grad_y = 0.5 * (grid(x, yp).fluid.rho - grid(x, ym).fluid.rho);

        return {grad_x, grad_y};
    }

    /**
     * Compute interface normal vector (points from gas to liquid)
     * Normal is -∇ρ / |∇ρ|
     */
    template<int D>
    std::array<double, D> computeInterfaceNormal(const std::array<double, D>& gradient) const {
        double grad_mag = 0.0;
        for (int d = 0; d < D; ++d) {
            grad_mag += gradient[d] * gradient[d];
        }
        grad_mag = std::sqrt(grad_mag);

        std::array<double, D> normal;
        if (grad_mag > 1e-12) {
            for (int d = 0; d < D; ++d) {
                normal[d] = -gradient[d] / grad_mag;
            }
        } else {
            normal.fill(0.0);
        }

        return normal;
    }

    /**
     * Analyze interface properties for a droplet configuration
     * @param grid - Lattice grid
     * @param center_x, center_y - Expected droplet center (for radius calculation)
     * @return Interface statistics
     */
    template<int D, int Q>
    InterfaceStats analyzeInterface(const LatticeGrid<D, Q>& grid,
                                   int center_x, int center_y) const {
        InterfaceStats stats;
        stats.mean_position_x = 0.0;
        stats.mean_position_y = 0.0;
        stats.num_interface_nodes = 0;
        stats.spurious_vel_interface = 0.0;

        std::vector<double> radii;
        radii.reserve(100);

        int nx = grid.nx();
        int ny = grid.ny();

        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                const auto& node = grid(x, y);

                if (node.is_solid) continue;

                if (isInterface(node.fluid.rho)) {
                    stats.mean_position_x += x;
                    stats.mean_position_y += y;
                    stats.num_interface_nodes++;

                    // Track spurious velocity at interface
                    double vel_mag = 0.0;
                    for (int d = 0; d < D; ++d) {
                        vel_mag += node.fluid.u[d] * node.fluid.u[d];
                    }
                    vel_mag = std::sqrt(vel_mag);
                    stats.spurious_vel_interface = std::max(stats.spurious_vel_interface, vel_mag);

                    // Compute distance from center for radius estimation
                    double dx = x - center_x;
                    double dy = y - center_y;
                    double r = std::sqrt(dx*dx + dy*dy);
                    radii.push_back(r);
                }
            }
        }

        if (stats.num_interface_nodes > 0) {
            stats.mean_position_x /= stats.num_interface_nodes;
            stats.mean_position_y /= stats.num_interface_nodes;

            // Compute mean radius
            double mean_radius = 0.0;
            for (double r : radii) {
                mean_radius += r;
            }
            stats.radius = mean_radius / radii.size();

            // Compute interface width (std dev of radii gives approximate width)
            double variance = 0.0;
            for (double r : radii) {
                double diff = r - stats.radius;
                variance += diff * diff;
            }
            stats.width = std::sqrt(variance / radii.size());
        } else {
            stats.radius = 0.0;
            stats.width = 0.0;
        }

        return stats;
    }

    /**
     * Compute average pressure in each phase
     * @return [P_liquid, P_gas, P_interface]
     */
    template<int D, int Q>
    std::array<double, 3> computePhaseAveragedPressure(const LatticeGrid<D, Q>& grid) const {
        double P_liquid = 0.0;
        double P_gas = 0.0;
        double P_interface = 0.0;

        int count_liquid = 0;
        int count_gas = 0;
        int count_interface = 0;

        for (const auto& node : grid) {
            if (node.is_solid) continue;

            Phase phase = classify(node.fluid.rho);
            double pressure = node.fluid.p;

            switch (phase) {
                case Phase::LIQUID:
                    P_liquid += pressure;
                    count_liquid++;
                    break;
                case Phase::GAS:
                    P_gas += pressure;
                    count_gas++;
                    break;
                case Phase::INTERFACE:
                    P_interface += pressure;
                    count_interface++;
                    break;
                default:
                    break;
            }
        }

        if (count_liquid > 0) P_liquid /= count_liquid;
        if (count_gas > 0) P_gas /= count_gas;
        if (count_interface > 0) P_interface /= count_interface;

        return {P_liquid, P_gas, P_interface};
    }

    /**
     * Compute Laplace pressure and surface tension
     * Uses Laplace law: ΔP = σ/R
     * @param P_liquid - Pressure inside droplet
     * @param P_gas - Pressure outside droplet
     * @param radius - Droplet radius
     * @return [laplace_pressure, surface_tension]
     */
    std::array<double, 2> computeLaplaceParameters(double P_liquid, double P_gas, double radius) const {
        double laplace_pressure = P_liquid - P_gas;
        double surface_tension = (radius > 1e-10) ? laplace_pressure * radius : 0.0;
        return {laplace_pressure, surface_tension};
    }

    // Getters
    double getRhoLiquid() const { return rho_liquid_; }
    double getRhoGas() const { return rho_gas_; }
    double getThreshold() const { return threshold_; }

private:
    double rho_liquid_;            // Liquid phase density
    double rho_gas_;               // Gas phase density
    double interface_threshold_;   // Fraction of density range
    double threshold_;             // Absolute density threshold
    double rho_interface_lower_;   // Lower bound of interface region
    double rho_interface_upper_;   // Upper bound of interface region
};

} // namespace elbm

#endif // INTERFACE_DETECTOR_H
