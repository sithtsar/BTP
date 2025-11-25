#ifndef COLOR_GRADIENT_H
#define COLOR_GRADIENT_H

#include "core/lattice.h"
#include "core/fluid_state.h"
#include <array>
#include <cmath>
#include <algorithm>

namespace elbm {

/**
 * Color-Gradient Two-Phase Model
 *
 * Implements the Rothman-Keller color-gradient method for immiscible two-phase flows.
 * Based on ETH Zurich research: "Thermodynamic Consistency, Interface Stability,
 * and Active Matter Dynamics" (2024).
 *
 * Key Features:
 * - Tracks red and blue distribution functions separately
 * - Three-step algorithm: Collision → Perturbation → Recoloring
 * - Independently tunable surface tension
 * - Minimal spurious currents
 * - Compatible with ELBM stabilization
 *
 * References:
 * - Rothman & Keller (1988) - Original color-gradient model
 * - Latva-Kokko & Rothman (2005) - Diffusion properties
 * - ETH Research Paper - Entropic stabilization at interfaces
 */

template<typename Lattice>
class ColorGradientModel {
public:
    using DistArray = std::array<double, Lattice::Q>;
    using VecArray = std::array<double, Lattice::D>;

    /**
     * Color-Gradient node state
     * Stores red and blue distributions separately
     */
    struct ColorNode {
        DistArray f_red;      // Red fluid distribution
        DistArray f_blue;     // Blue fluid distribution
        double rho_red;       // Red density
        double rho_blue;      // Blue density
        double rho_total;     // Total density
        VecArray u;           // Velocity (shared)
        double color_field;   // ρ^N = (ρ^R - ρ^B)/(ρ^R + ρ^B)
        VecArray color_gradient;  // ∇ρ^N
        VecArray interface_normal; // n = -∇ρ^N / |∇ρ^N|

        ColorNode() {
            f_red.fill(0.0);
            f_blue.fill(0.0);
            rho_red = 0.0;
            rho_blue = 0.0;
            rho_total = 0.0;
            u.fill(0.0);
            color_field = 0.0;
            color_gradient.fill(0.0);
            interface_normal.fill(0.0);
        }
    };

    /**
     * Constructor
     * @param surface_tension - Surface tension parameter σ
     * @param beta_recolor - Recoloring strength parameter (default: 0.9, increased from 0.7)
     * CORRECTED: Increased default from 0.7 to 0.9 for stronger phase separation
     */
    ColorGradientModel(double surface_tension, double beta_recolor = 0.99)
        : sigma_(surface_tension)
        , beta_recolor_(beta_recolor) {

        // Compute perturbation amplitude A from desired surface tension
        // Increased amplitude for effective surface tension: A = σ × 10
        perturbation_amplitude_ = surface_tension * 10.0;  // = 0.01 for σ=0.001
    }

    /**
     * Initialize droplet configuration
     * @param grid - Grid of color nodes
     * @param center_x, center_y - Droplet center
     * @param radius - Droplet radius
     * @param rho_inside - Density inside droplet (red phase)
     * @param rho_outside - Density outside droplet (blue phase)
     */
    template<int D, int Q>
    void initializeDroplet(std::vector<std::vector<ColorNode>>& grid,
                          int nx, int ny,
                          int center_x, int center_y, double radius,
                          double rho_inside, double rho_outside) {

        static_assert(D == 2 && Q == 9, "Color-Gradient currently supports D2Q9 only");

        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                auto& node = grid[y][x];

                // Distance from center
                double dx = x - center_x;
                double dy = y - center_y;
                double r = std::sqrt(dx*dx + dy*dy);

                // Smooth interface using tanh profile
                // Further optimized: Width = 0.5 for very sharp interface
                double interface_width = 0.5;
                double rho_red = rho_inside * 0.5 * (1.0 - std::tanh(2.0 * (r - radius) / interface_width));
                double rho_blue = rho_outside * 0.5 * (1.0 + std::tanh(2.0 * (r - radius) / interface_width));

                node.rho_red = rho_red;
                node.rho_blue = rho_blue;
                node.rho_total = rho_red + rho_blue;
                node.u.fill(0.0);

                // Initialize color field
                if (node.rho_total > 1e-12) {
                    node.color_field = (rho_red - rho_blue) / node.rho_total;
                } else {
                    node.color_field = 0.0;
                }

                // Initialize distributions at equilibrium
                // Distribute total mass into red and blue based on color field
                for (int i = 0; i < Q; ++i) {
                    double f_eq_total = computeEquilibrium<D, Q>(i, node.rho_total, node.u);

                    // Simple initial distribution: proportional to density
                    node.f_red[i] = (rho_red / node.rho_total) * f_eq_total;
                    node.f_blue[i] = (rho_blue / node.rho_total) * f_eq_total;
                }
            }
        }
    }

    /**
     * Compute equilibrium distribution (2nd order)
     * Uses standard BGK polynomial equilibrium for initial state
     */
    template<int D, int Q>
    double computeEquilibrium(int i, double rho, const VecArray& u) const {
        constexpr double cs2 = 1.0 / 3.0;  // c_s² for D2Q9
        constexpr double cs4 = cs2 * cs2;

        double ci_dot_u = 0.0;
        for (int d = 0; d < D; ++d) {
            ci_dot_u += Lattice::c[i][d] * u[d];
        }

        double u_sq = 0.0;
        for (int d = 0; d < D; ++d) {
            u_sq += u[d] * u[d];
        }

        double f_eq = Lattice::weight(i) * rho * (
            1.0 + ci_dot_u / cs2 +
            ci_dot_u * ci_dot_u / (2.0 * cs4) -
            u_sq / (2.0 * cs2)
        );

        return f_eq;
    }

    /**
     * Compute color field and gradient for all nodes
     */
    template<int D, int Q>
    void computeColorFieldAndGradient(std::vector<std::vector<ColorNode>>& grid, int nx, int ny) {
        // Compute color field for all nodes
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                auto& node = grid[y][x];

                if (node.rho_total > 1e-12) {
                    node.color_field = (node.rho_red - node.rho_blue) / node.rho_total;
                } else {
                    node.color_field = 0.0;
                }
            }
        }

        // Compute gradients using central differences with periodic BC
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                auto& node = grid[y][x];

                int xp = (x + 1) % nx;
                int xm = (x - 1 + nx) % nx;
                int yp = (y + 1) % ny;
                int ym = (y - 1 + ny) % ny;

                // Central difference gradient
                node.color_gradient[0] = 0.5 * (grid[y][xp].color_field - grid[y][xm].color_field);
                node.color_gradient[1] = 0.5 * (grid[yp][x].color_field - grid[ym][x].color_field);

                // Compute interface normal: n = -∇ρ^N / |∇ρ^N|
                double grad_mag = std::sqrt(
                    node.color_gradient[0] * node.color_gradient[0] +
                    node.color_gradient[1] * node.color_gradient[1]
                );

                if (grad_mag > 1e-12) {
                    node.interface_normal[0] = -node.color_gradient[0] / grad_mag;
                    node.interface_normal[1] = -node.color_gradient[1] / grad_mag;
                } else {
                    node.interface_normal[0] = 0.0;
                    node.interface_normal[1] = 0.0;
                }
            }
        }
    }

    /**
     * Step 2: Compute perturbation for surface tension
     *
     * From ETH paper: Perturbation creates anisotropic stress at interface
     * Formula: Δf_i = A·w_i·|∇ρ^N|·[(c_i·n)² - c_s²]
     *
     * This creates tension perpendicular to interface normal
     * Returns the perturbation array for thermodynamic consistency in ELBM
     */
    template<int D, int Q>
    DistArray computePerturbation(ColorNode& node) const {
        constexpr double cs2 = 1.0 / 3.0;
        DistArray delta_surf{};
        delta_surf.fill(0.0);

        double grad_mag = std::sqrt(
            node.color_gradient[0] * node.color_gradient[0] +
            node.color_gradient[1] * node.color_gradient[1]
        );

        // CRITICAL FIX: Remove threshold check - perturbation formula naturally scales with grad_mag
        // Threshold created death spiral: weak gradient → no perturbation → weaker gradient
        // Formula includes grad_mag factor, so automatically weak in bulk regions
        // Only skip if gradient is exactly zero (division by zero in normal vector)

        for (int i = 0; i < Q; ++i) {
            // Compute c_i · n
            double ci_dot_n = 0.0;
            for (int d = 0; d < D; ++d) {
                ci_dot_n += Lattice::c[i][d] * node.interface_normal[d];
            }

            // Perturbation: A·w_i·|∇ρ^N|·[(c_i·n)² - c_s²]
            // Positive sign for surface tension (Gunn-Rothman formula)
            double perturbation = perturbation_amplitude_ * Lattice::weight(i) *
                                  grad_mag * (ci_dot_n * ci_dot_n - cs2);

            delta_surf[i] = perturbation;
        }

        return delta_surf;
    }

    /**
     * Step 3: Recoloring - Maximize phase separation
     *
     * From ETH paper (Latva-Kokko & Rothman 2005):
     * Redistributes f into f^R and f^B to maximize color separation
     *
     * Formula:
     * f_i^R = (ρ^R/ρ)·f_i + β·ρ^R·ρ^B·cos(φ_i)·|f_i^eq|/ρ²
     * f_i^B = (ρ^B/ρ)·f_i - β·ρ^R·ρ^B·cos(φ_i)·|f_i^eq|/ρ²
     *
     * where cos(φ_i) = (c_i · n) / |c_i|
     */
    template<int D, int Q>
    void recolor(ColorNode& node, const DistArray& f_total) {
        if (node.rho_total < 1e-12) return;

        double rho_red = node.rho_red;
        double rho_blue = node.rho_blue;
        double rho = node.rho_total;

        // Compute |∇ρ^N| to detect interface
        double grad_mag = std::sqrt(
            node.color_gradient[0] * node.color_gradient[0] +
            node.color_gradient[1] * node.color_gradient[1]
        );

        // CRITICAL FIX: ALWAYS apply recoloring to prevent interface diffusion death spiral
        // Previous threshold check (grad_mag < 1e-9) created vicious cycle:
        //   diffusion → weak gradient → no recoloring → more diffusion → death
        // Now: continuous recoloring maintains interface sharpness even with weak gradients

        // For very weak gradients (true bulk regions), normal vector becomes undefined
        // In this case, recoloring reduces to simple proportional split (cos_phi ≈ 0)
        // So no explicit threshold check needed - algorithm handles it naturally
        for (int i = 0; i < Q; ++i) {
            // Compute cos(φ_i) = (c_i · n) / |c_i|
            double ci_dot_n = 0.0;
            double ci_mag = 0.0;
            for (int d = 0; d < D; ++d) {
                ci_dot_n += Lattice::c[i][d] * node.interface_normal[d];
                ci_mag += Lattice::c[i][d] * Lattice::c[i][d];
            }
            ci_mag = std::sqrt(ci_mag);

            double cos_phi = (ci_mag > 1e-12) ? (ci_dot_n / ci_mag) : 0.0;

            // Equilibrium distribution magnitude
            double f_eq = computeEquilibrium<D, Q>(i, rho, node.u);

            // Recoloring formula
            double recolor_term = beta_recolor_ * rho_red * rho_blue * cos_phi * std::abs(f_eq) / (rho * rho);

            node.f_red[i] = (rho_red / rho) * f_total[i] + recolor_term;
            node.f_blue[i] = (rho_blue / rho) * f_total[i] - recolor_term;

            // Ensure positivity
            node.f_red[i] = std::max(node.f_red[i], 0.0);
            node.f_blue[i] = std::max(node.f_blue[i], 0.0);
        }
    }

    /**
     * Compute macroscopic quantities from colored distributions
     */
    template<int D, int Q>
    void computeMacroscopic(ColorNode& node) {
        node.rho_red = 0.0;
        node.rho_blue = 0.0;
        node.u.fill(0.0);

        // Sum densities
        for (int i = 0; i < Q; ++i) {
            node.rho_red += node.f_red[i];
            node.rho_blue += node.f_blue[i];
        }

        node.rho_total = node.rho_red + node.rho_blue;

        // Compute velocity (momentum / total density)
        if (node.rho_total > 1e-12) {
            for (int d = 0; d < D; ++d) {
                double momentum = 0.0;
                for (int i = 0; i < Q; ++i) {
                    momentum += (node.f_red[i] + node.f_blue[i]) * Lattice::c[i][d];
                }
                node.u[d] = momentum / node.rho_total;
            }
        }

        // Update color field
        if (node.rho_total > 1e-12) {
            node.color_field = (node.rho_red - node.rho_blue) / node.rho_total;
        } else {
            node.color_field = 0.0;
        }
    }

    // Getters
    double getSurfaceTension() const { return sigma_; }
    double getBetaRecolor() const { return beta_recolor_; }
    double getPerturbationAmplitude() const { return perturbation_amplitude_; }

private:
    double sigma_;                   // Surface tension
    double beta_recolor_;            // Recoloring strength (0-1)
    double perturbation_amplitude_;  // Perturbation amplitude A
};

} // namespace elbm

#endif // COLOR_GRADIENT_H
