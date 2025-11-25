#ifndef ELBM_SOLVER_H
#define ELBM_SOLVER_H

#include "core/lattice.h"
#include "core/fluid_state.h"
#include "solvers/equilibrium.h"
#include <cmath>
#include <algorithm>
#include <limits>

namespace elbm {

template<typename Lattice>
class ELBMSolver {
public:
    ELBMSolver(double viscosity, double dt = 1.0, double dx = 1.0,
               int max_iter = 100, double tol = 1e-10)
        : dt_(dt), dx_(dx), max_iter_(max_iter), tol_(tol) {
        // Viscosity will be controlled by alpha and beta
        // ν = cs² (1/(αβ) - 1/2) δt
        target_viscosity_ = viscosity;

        // Calculate beta from target viscosity, assuming alpha converges to 2 in smooth flows.
        // This is the key to matching BGK results in simple shear/channel flows.
        beta_ = 0.5 / (target_viscosity_ / Lattice::cs2 + 0.5);
    }



    // Compute H-function (discrete entropy)
    // H = Σ f_i · ln(f_i / w_i)
    double computeH(const std::array<double, Lattice::Q>& f) const {
        constexpr double eps = 1e-12; // Conservative threshold
        double H = 0.0;

        for (int i = 0; i < Lattice::Q; ++i) {
            if (f[i] > eps) {
                const double log_term = std::log(f[i] / Lattice::weight(i));

                // Safety check for numerical issues
                if (std::isfinite(log_term)) {
                    H += f[i] * log_term;
                }
            }
        }

        // Check for NaN/Inf in result
        if (!std::isfinite(H)) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        return H;
    }

    // Compute bounds for alpha parameter
    // From thesis eq 3.46: Ensure f + α·Δ ≥ 0 for all i
    // When Δ_i < 0: α ≤ f_i/|Δ_i| (upper bound)
    // When Δ_i > 0: α ≥ -f_i/Δ_i (lower bound, usually satisfied by α≥0)
    void computeAlphaBounds(const std::array<double, Lattice::Q>& f,
                           const std::array<double, Lattice::Q>& f_eq,
                           double& alpha_min, double& alpha_max) const {
        alpha_min = 0.0;
        alpha_max = 2.0;

        constexpr double eps = 1e-14;

        for (int i = 0; i < Lattice::Q; ++i) {
            const double delta_i = f_eq[i] - f[i];

            if (std::abs(delta_i) > eps) {
                // For f + α·Δ ≥ 0: α ≥ -f/Δ (if Δ>0) or α ≤ -f/Δ (if Δ<0)
                const double bound = -f[i] / delta_i;

                if (delta_i < 0.0) {
                    // Negative delta: α ≤ f_i/|Δ_i| = -f_i/Δ_i (UPPER bound)
                    alpha_max = std::min(alpha_max, bound);
                } else {
                    // Positive delta: α ≥ -f_i/Δ_i (LOWER bound)
                    alpha_min = std::max(alpha_min, bound);
                }
            }
        }

        // Ensure valid bounds
        alpha_min = std::max(0.0, alpha_min);
        alpha_max = std::min(2.0, alpha_max);

        // Ensure min < max
        if (alpha_min >= alpha_max) {
            alpha_min = 0.0;
            alpha_max = 2.0;
        }
    }

    // Solve for alpha parameter using Newton-Raphson
    // Find α such that H(f + α·Δ + Δ_surf) = H(f), where Δ = f_eq - f
    // Δ_surf is the surface tension perturbation (for Color-Gradient thermodynamic consistency)
    double solveAlpha(const std::array<double, Lattice::Q>& f,
                      const std::array<double, Lattice::Q>& f_eq,
                      const std::array<double, Lattice::Q>& delta_surf = {}) const {
        // Check if perturbation is provided (non-empty)
        bool has_perturbation = !delta_surf.empty();

        // Compute initial H
        const double H0 = computeH(f);

        // Compute delta
        std::array<double, Lattice::Q> delta;
        for (int i = 0; i < Lattice::Q; ++i) {
            delta[i] = f_eq[i] - f[i];
        }

        // Check if already at equilibrium
        double delta_norm = 0.0;
        for (int i = 0; i < Lattice::Q; ++i) {
            delta_norm += delta[i] * delta[i];
        }
        if (delta_norm < 1e-16) {
            return 2.0; // Already at equilibrium
        }

        // Compute bounds
        double alpha_min, alpha_max;
        computeAlphaBounds(f, f_eq, alpha_min, alpha_max);

        // Check boundary cases
        std::array<double, Lattice::Q> f_test;

        // Test alpha_max
        bool alpha_max_valid = true;
        for (int i = 0; i < Lattice::Q; ++i) {
            f_test[i] = f[i] + alpha_max * delta[i];
            if (has_perturbation) {
                f_test[i] += delta_surf[i];
            }
            if (f_test[i] < 0.0) {
                alpha_max_valid = false;
                break;
            }
        }

        if (alpha_max_valid) {
            const double H_max = computeH(f_test);
            if (std::abs(H_max - H0) < tol_) {
                return alpha_max;
            }
        }

        // Newton-Raphson iteration
        double alpha = 1.0; // Initial guess

        for (int iter = 0; iter < max_iter_; ++iter) {
            // Clamp alpha to bounds
            alpha = std::clamp(alpha, alpha_min + 1e-10, alpha_max - 1e-10);

            // Compute f_alpha = f + alpha * delta + delta_surf
            for (int i = 0; i < Lattice::Q; ++i) {
                f_test[i] = f[i] + alpha * delta[i];
                if (has_perturbation) {
                    f_test[i] += delta_surf[i];
                }
            }

            // Compute H(alpha) and its derivative
            double H_alpha = computeH(f_test);
            double dH_dalpha = 0.0;

            for (int i = 0; i < Lattice::Q; ++i) {
                if (f_test[i] > 1e-16) {
                    dH_dalpha += delta[i] * (std::log(f_test[i] / Lattice::weight(i)) + 1.0);
                }
            }

            // Check convergence
            const double residual = H_alpha - H0;
            if (std::abs(residual) < tol_) {
                return alpha;
            }

            // Newton update
            if (std::abs(dH_dalpha) > 1e-16) {
                alpha -= residual / dH_dalpha;
            } else {
                // Bisection fallback
                if (residual > 0.0) {
                    alpha_max = alpha;
                    alpha = 0.5 * (alpha_min + alpha_max);
                } else {
                    alpha_min = alpha;
                    alpha = 0.5 * (alpha_min + alpha_max);
                }
            }
        }

        // If not converged, return safe value
        return std::clamp(alpha, alpha_min, alpha_max);
    }

    // Two-step entropic collision
    // delta_surf: surface tension perturbation (for Color-Gradient thermodynamic consistency)
    void collide(Node<Lattice::D, Lattice::Q>& node,
                 const std::array<double, Lattice::D>& force = {},
                 const std::array<double, Lattice::Q>& delta_surf = {}) {
        auto& df = node.df;
        auto& fluid = node.fluid;

        // Step 1: Compute macroscopic variables
        computeMacroscopic<Lattice>(df.f, fluid);

        // Step 2: Compute equilibrium distribution
        // TEMPORARY: Use BGK equilibrium for Color-Gradient stability testing
        // TODO: Switch back to EntropicEquilibrium after stabilization
        BGKEquilibrium<Lattice>::compute(fluid, df.f_eq);

        // Step 3: Solve for alpha (iso-entropy step)
        // Include surface tension perturbation for thermodynamic consistency
        const double alpha = solveAlpha(df.f, df.f_eq, delta_surf);
        node.alpha = alpha;

        // Step 4: Compute auxiliary distribution f*
        for (int i = 0; i < Lattice::Q; ++i) {
            df.f_star[i] = df.f[i] + alpha * (df.f_eq[i] - df.f[i]);
        }

        // Safety check: Ensure f* is valid and positive
        bool f_star_valid = true;
        constexpr double eps = 1e-12;
        for (int i = 0; i < Lattice::Q; ++i) {
            if (!std::isfinite(df.f_star[i]) || df.f_star[i] < -eps) {
                f_star_valid = false;
                break;
            }
            // Clamp small negative values to zero
            if (df.f_star[i] < 0.0) {
                df.f_star[i] = 0.0;
            }
        }

        // If f* is invalid, fall back to BGK collision
        if (!f_star_valid) {
            for (int i = 0; i < Lattice::Q; ++i) {
                df.f_new[i] = df.f[i] - 0.5 * (df.f[i] - df.f_eq[i]);
            }
            node.alpha = 1.0;
            node.beta = 0.5;
            return;
        }

        // Step 5: Dissipation step with beta
        // f_new = (1-β)f + β·f*
        for (int i = 0; i < Lattice::Q; ++i) {
            df.f_new[i] = (1.0 - beta_) * df.f[i] + beta_ * df.f_star[i];

            // Final safety check
            if (!std::isfinite(df.f_new[i])) {
                df.f_new[i] = df.f[i]; // Keep old value
            }
        }

        // Step 6: Add body force term if present (Guo forcing scheme)
        bool has_force = false;
        if (!force.empty()) {
            for (int d = 0; d < Lattice::D; ++d) {
                if (std::abs(force[d]) > 1e-14) {
                    has_force = true;
                    break;
                }
            }
        }

        if (has_force) {
            // Use average alpha for forcing (alpha from current step)
            const double omega_eff = 1.0 / (node.alpha * beta_);

            for (int i = 0; i < Lattice::Q; ++i) {
                // Compute ci · u and ci · F
                double ci_dot_u = 0.0;
                double ci_dot_F = 0.0;

                if constexpr (Lattice::D == 2) {
                    ci_dot_u = Lattice::cx(i) * fluid.u[0] + Lattice::cy(i) * fluid.u[1];
                    ci_dot_F = Lattice::cx(i) * force[0] + Lattice::cy(i) * force[1];
                } else if constexpr (Lattice::D == 3) {
                    ci_dot_u = Lattice::cx(i) * fluid.u[0] + Lattice::cy(i) * fluid.u[1] +
                              Lattice::cz(i) * fluid.u[2];
                    ci_dot_F = Lattice::cx(i) * force[0] + Lattice::cy(i) * force[1] +
                              Lattice::cz(i) * force[2];
                }

                // u · F
                double u_dot_F = 0.0;
                for (int d = 0; d < Lattice::D; ++d) {
                    u_dot_F += fluid.u[d] * force[d];
                }

                // Guo forcing term
                double force_term = Lattice::weight(i) * (1.0 - 0.5 * omega_eff) *
                                   (ci_dot_F / Lattice::cs2 +
                                    (ci_dot_u * ci_dot_F) / (Lattice::cs2 * Lattice::cs2) -
                                    u_dot_F / Lattice::cs2);

                df.f_new[i] += force_term;
            }
        }

        // Update beta for next iteration
        node.beta = beta_;
    }

    // Streaming step with selectable periodic boundaries
    void stream(LatticeGrid<Lattice::D, Lattice::Q>& grid,
                bool periodic_x = false, bool periodic_y = false, bool periodic_z = false) {
        const int nx = grid.nx();
        const int ny = grid.ny();
        const int nz = (Lattice::D == 3) ? grid.nz() : 1;

        // Create temporary storage for streaming
        LatticeGrid<Lattice::D, Lattice::Q> temp(nx, ny, nz);

        // Copy post-collision distributions to temporary grid
        #pragma omp parallel for collapse(3)
        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    temp(x, y, z).df.f = grid(x, y, z).df.f_new;
                }
            }
        }

        // Stream to neighbor nodes
        #pragma omp parallel for collapse(3)
        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    // Skip streaming from solid nodes
                    if (grid(x, y, z).is_solid) continue;

                    for (int i = 0; i < Lattice::Q; ++i) {
                        int xn, yn, zn;

                        if constexpr (Lattice::D == 2) {
                            xn = x + Lattice::cx(i);
                            yn = y + Lattice::cy(i);
                            zn = 0;
                            if (periodic_x) xn = (xn + nx) % nx;
                            if (periodic_y) yn = (yn + ny) % ny;
                        } else {
                            xn = x + Lattice::cx(i);
                            yn = y + Lattice::cy(i);
                            zn = z + Lattice::cz(i);
                            if (periodic_x) xn = (xn + nx) % nx;
                            if (periodic_y) yn = (yn + ny) % ny;
                            if (periodic_z) zn = (zn + nz) % nz;
                        }

                        // Skip if out of bounds (for non-periodic directions)
                        if (xn < 0 || xn >= nx || yn < 0 || yn >= ny) continue;
                        if constexpr (Lattice::D == 3) {
                            if (zn < 0 || zn >= nz) continue;
                        }

                        // Skip streaming to solid nodes
                        if (grid(xn, yn, zn).is_solid) continue;

                        // Stream to neighbor
                        grid(xn, yn, zn).df.f[i] = temp(x, y, z).df.f[i];
                    }
                }
            }
        }
    }

    // Get effective viscosity
    double viscosity(double alpha) const {
        // ν = cs² (1/(αβ) - 1/2) δt
        return Lattice::cs2 * (1.0 / (alpha * beta_) - 0.5) * dt_;
    }

    double beta() const { return beta_; }
    void setBeta(double beta) { beta_ = beta; }

private:
    double dt_;                               // Time step
    double dx_;                               // Lattice spacing
    double beta_;                             // Dissipation parameter
    double target_viscosity_;                 // Target viscosity
    int max_iter_;                            // Max iterations for alpha solver
    double tol_;                              // Tolerance for alpha solver
};

} // namespace elbm

#endif // ELBM_SOLVER_H
