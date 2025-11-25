#ifndef VDW_FORCE_H
#define VDW_FORCE_H

#include "core/lattice.h"
#include "core/fluid_state.h"
#include <array>
#include <cmath>

namespace elbm {

/**
 * Van der Waals-Korteweg Force Calculator
 *
 * Implements the free-energy approach for two-phase immiscible flows
 * Force: F = -κ·ρ·∇∇²ρ + ∇P_non-ideal
 *
 * Based on:
 * - Swift et al. (1995) PRL 75:830
 * - Mazloomi et al. (2015) PRL 114:174502 (ELBM extension)
 */
template<typename Lattice, int NX, int NY>
class VanDerWaalsForce {
public:
    VanDerWaalsForce(double kappa, double a, double b)
        : kappa_(kappa), a_(a), b_(b) {}

    /**
     * Compute van der Waals force at a node
     * F = -κ·ρ·∇(∇²ρ)
     */
    std::array<double, Lattice::D> computeForce(
        const LatticeGrid<Lattice::D, Lattice::Q>& grid,
        int x, int y) const {

        std::array<double, Lattice::D> force;
        force.fill(0.0);

        double rho = grid(x, y).fluid.rho;

        // Compute Laplacian of density using finite differences
        double laplacian_rho = computeLaplacian(grid, x, y);

        // Compute gradient of Laplacian
        auto grad_laplacian = computeGradient(grid, x, y,
            [this](const auto& g, int i, int j) { return computeLaplacian(g, i, j); });

        // F = -κ·ρ·∇(∇²ρ)
        force[0] = -kappa_ * rho * grad_laplacian[0];
        force[1] = -kappa_ * rho * grad_laplacian[1];

        return force;
    }

    /**
     * Compute pressure from van der Waals equation of state
     * P = ρ·cs² - a·ρ² + b·ρ³
     */
    double computePressure(double rho) const {
        constexpr double cs2 = 1.0 / 3.0;
        return rho * cs2 - a_ * rho * rho + b_ * rho * rho * rho;
    }

    /**
     * Get parameters for coexistence curve
     * For phase separation: ρ_liquid and ρ_gas
     */
    std::pair<double, double> getCoexistenceDensities() const {
        // Simplified: for demo purposes
        // Proper calculation requires solving van der Waals EOS
        double rho_gas = 0.2;
        double rho_liquid = 2.0;
        return {rho_gas, rho_liquid};
    }

private:
    double kappa_;  // Capillary coefficient (interface width parameter)
    double a_;      // Attractive interaction parameter
    double b_;      // Repulsive interaction parameter

    /**
     * Compute Laplacian using 9-point stencil (for D2Q9)
     * ∇²ρ = (ρ_E + ρ_W + ρ_N + ρ_S - 4ρ_C) / Δx²
     *       + 0.5*(ρ_NE + ρ_NW + ρ_SE + ρ_SW - 4ρ_C) / Δx²
     */
    double computeLaplacian(
        const LatticeGrid<Lattice::D, Lattice::Q>& grid,
        int x, int y) const {

        double rho_C = grid(x, y).fluid.rho;

        // Periodic boundary conditions
        int xp = (x + 1) % NX;
        int xm = (x - 1 + NX) % NX;
        int yp = (y + 1) % NY;
        int ym = (y - 1 + NY) % NY;

        // Cardinal directions (weight 1)
        double rho_E = grid(xp, y).fluid.rho;
        double rho_W = grid(xm, y).fluid.rho;
        double rho_N = grid(x, yp).fluid.rho;
        double rho_S = grid(x, ym).fluid.rho;

        // Diagonal directions (weight 0.5)
        double rho_NE = grid(xp, yp).fluid.rho;
        double rho_NW = grid(xm, yp).fluid.rho;
        double rho_SE = grid(xp, ym).fluid.rho;
        double rho_SW = grid(xm, ym).fluid.rho;

        // 9-point Laplacian
        double laplacian = (rho_E + rho_W + rho_N + rho_S - 4.0 * rho_C)
                         + 0.5 * (rho_NE + rho_NW + rho_SE + rho_SW - 4.0 * rho_C);

        return laplacian; // Δx² = 1 in lattice units
    }

    /**
     * Compute gradient using central differences
     * Generic template to compute gradient of any field
     */
    template<typename FieldFunc>
    std::array<double, Lattice::D> computeGradient(
        const LatticeGrid<Lattice::D, Lattice::Q>& grid,
        int x, int y,
        FieldFunc field) const {

        std::array<double, Lattice::D> grad;

        // Periodic boundary conditions
        int xp = (x + 1) % NX;
        int xm = (x - 1 + NX) % NX;
        int yp = (y + 1) % NY;
        int ym = (y - 1 + NY) % NY;

        // Central difference (2nd order accurate)
        grad[0] = (field(grid, xp, y) - field(grid, xm, y)) / 2.0;
        grad[1] = (field(grid, x, yp) - field(grid, x, ym)) / 2.0;

        return grad;
    }
};

} // namespace elbm

#endif // VDW_FORCE_H
