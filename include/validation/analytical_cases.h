#ifndef ANALYTICAL_CASES_H
#define ANALYTICAL_CASES_H

#include "core/lattice.h"
#include "core/fluid_state.h"
#include <cmath>
#include <complex>
#include <array>

namespace elbm {
namespace validation {

/**
 * Couette Flow: Flow between two parallel plates with relative motion
 * Analytical solution: u(y) = U * y / H
 * where U is the top plate velocity and H is the channel height
 */
template<typename Lattice>
class CouetteFlow {
public:
    CouetteFlow(int nx, int ny, double U_wall, double viscosity)
        : nx_(nx), ny_(ny), U_wall_(U_wall), viscosity_(viscosity) {}

    // Initialize velocity field
    void initialize(LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                auto& node = grid(x, y);

                // Linear velocity profile
                const double u_x = U_wall_ * static_cast<double>(y) / (ny_ - 1);

                node.fluid.rho = 1.0;
                node.fluid.u[0] = u_x;
                node.fluid.u[1] = 0.0;

                // Initialize distributions to equilibrium
                BGKEquilibrium<Lattice>::compute(node.fluid, node.df.f);
            }
        }
    }

    // Compute analytical solution
    std::vector<double> analytical_solution() const {
        std::vector<double> u_analytical(ny_);
        for (int y = 0; y < ny_; ++y) {
            u_analytical[y] = U_wall_ * static_cast<double>(y) / (ny_ - 1);
        }
        return u_analytical;
    }

    // Compute L2 error
    double compute_error(const LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        auto u_exact = analytical_solution();
        double error = 0.0;

        const int x_mid = nx_ / 2;

        for (int y = 0; y < ny_; ++y) {
            const double u_sim = grid(x_mid, y).fluid.u[0];
            const double diff = u_sim - u_exact[y];
            error += diff * diff;
        }

        return std::sqrt(error / ny_);
    }

private:
    int nx_, ny_;
    double U_wall_;
    double viscosity_;
};

/**
 * Poiseuille Flow: Pressure-driven flow between two parallel plates
 * Analytical solution: u(y) = (dp/dx) * y * (H - y) / (2 * μ)
 * where dp/dx is the pressure gradient, H is channel height, μ is dynamic viscosity
 */
template<typename Lattice>
class PoiseuilleFlow {
public:
    PoiseuilleFlow(int nx, int ny, double dp_dx, double viscosity)
        : nx_(nx), ny_(ny), dp_dx_(dp_dx), viscosity_(viscosity) {}

    // Initialize from rest (let body force develop the flow)
    // This is the standard approach for body-force driven Poiseuille validation
    void initialize(LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                auto& node = grid(x, y);

                // Start from rest - body force will develop the parabolic profile
                node.fluid.rho = 1.0;
                node.fluid.u[0] = 0.0;  // Zero initial velocity
                node.fluid.u[1] = 0.0;

                // Initialize distributions to equilibrium at rest
                BGKEquilibrium<Lattice>::compute(node.fluid, node.df.f);
            }
        }
    }

    // Compute analytical solution
    std::vector<double> analytical_solution() const {
        std::vector<double> u_analytical(ny_);
        const double H = ny_ - 1;
        const double rho = 1.0;

        for (int y = 0; y < ny_; ++y) {
            const double y_pos = static_cast<double>(y);
            u_analytical[y] = -dp_dx_ * y_pos * (H - y_pos) / (2.0 * rho * viscosity_);
        }
        return u_analytical;
    }

    // Compute L2 error
    double compute_error(const LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        auto u_exact = analytical_solution();
        double error = 0.0;

        const int x_mid = nx_ / 2;

        for (int y = 0; y < ny_; ++y) {
            const double u_sim = grid(x_mid, y).fluid.u[0];
            const double diff = u_sim - u_exact[y];
            error += diff * diff;
        }

        return std::sqrt(error / ny_);
    }

private:
    int nx_, ny_;
    double dp_dx_;
    double viscosity_;
};

/**
 * Taylor-Green Vortex: 2D decaying vortex for testing viscosity extraction
 * Analytical solution:
 *   u_x(x,y,t) = -U0 * cos(kx) * sin(ky) * exp(-2k²νt)
 *   u_y(x,y,t) =  U0 * sin(kx) * cos(ky) * exp(-2k²νt)
 *   p(x,y,t) = -ρ0 * U0² / 4 * (cos(2kx) + cos(2ky)) * exp(-4k²νt)
 */
template<typename Lattice>
class TaylorGreenVortex {
public:
    TaylorGreenVortex(int nx, int ny, double U0, double viscosity, double k = 1.0)
        : nx_(nx), ny_(ny), U0_(U0), viscosity_(viscosity), k_(k) {}

    // Initialize velocity field at t=0
    void initialize(LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        const double Lx = 2.0 * M_PI / k_;
        const double Ly = 2.0 * M_PI / k_;

        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                auto& node = grid(x, y);

                // Physical coordinates
                const double x_pos = Lx * x / nx_;
                const double y_pos = Ly * y / ny_;

                // Velocity field at t=0
                const double u_x = -U0_ * std::cos(k_ * x_pos) * std::sin(k_ * y_pos);
                const double u_y =  U0_ * std::sin(k_ * x_pos) * std::cos(k_ * y_pos);

                node.fluid.rho = 1.0;
                node.fluid.u[0] = u_x;
                node.fluid.u[1] = u_y;

                // Initialize distributions to equilibrium
                BGKEquilibrium<Lattice>::compute(node.fluid, node.df.f);
            }
        }
    }

    // Compute analytical solution at time t
    void analytical_solution(double t, std::vector<std::array<double, 2>>& u_exact) const {
        u_exact.resize(nx_ * ny_);

        const double Lx = 2.0 * M_PI / k_;
        const double Ly = 2.0 * M_PI / k_;
        const double decay = std::exp(-2.0 * k_ * k_ * viscosity_ * t);

        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const double x_pos = Lx * x / nx_;
                const double y_pos = Ly * y / ny_;

                const int idx = x + nx_ * y;
                u_exact[idx][0] = -U0_ * std::cos(k_ * x_pos) * std::sin(k_ * y_pos) * decay;
                u_exact[idx][1] =  U0_ * std::sin(k_ * x_pos) * std::cos(k_ * y_pos) * decay;
            }
        }
    }

    // Compute kinetic energy (should decay as exp(-4k²νt))
    double compute_kinetic_energy(const LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        double energy = 0.0;

        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const auto& u = grid(x, y).fluid.u;
                energy += u[0] * u[0] + u[1] * u[1];
            }
        }

        return 0.5 * energy / (nx_ * ny_);
    }

    // Extract viscosity from energy decay
    double extract_viscosity(double E0, double E1, double dt) const {
        // E(t) = E0 * exp(-4k²νt)
        // ν = -ln(E1/E0) / (4k²dt)
        if (E1 > 0 && E0 > 0 && E1 < E0) {
            return -std::log(E1 / E0) / (4.0 * k_ * k_ * dt);
        }
        return -1.0; // Invalid
    }

    // Compute L2 error
    double compute_error(const LatticeGrid<Lattice::D, Lattice::Q>& grid, double t) const {
        std::vector<std::array<double, 2>> u_exact;
        analytical_solution(t, u_exact);

        double error = 0.0;

        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const int idx = x + nx_ * y;
                const auto& u_sim = grid(x, y).fluid.u;

                const double diff_x = u_sim[0] - u_exact[idx][0];
                const double diff_y = u_sim[1] - u_exact[idx][1];

                error += diff_x * diff_x + diff_y * diff_y;
            }
        }

        return std::sqrt(error / (nx_ * ny_));
    }

private:
    int nx_, ny_;
    double U0_;
    double viscosity_;
    double k_;
};

/**
 * Womersley Flow: Oscillating pressure-driven flow between parallel plates
 * Analytical solution for oscillatory flow (pulsatile flow)
 * u(y,t) = Re[A * (1 - cosh(λy)/cosh(λH/2)) * exp(iωt)]
 * where λ = sqrt(iω/ν), ω is angular frequency
 */
template<typename Lattice>
class WomersleyFlow {
public:
    WomersleyFlow(int nx, int ny, double amplitude, double frequency, double viscosity)
        : nx_(nx), ny_(ny), amplitude_(amplitude), frequency_(frequency),
          viscosity_(viscosity), omega_(2.0 * M_PI * frequency) {}

    // Initialize from rest
    void initialize(LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                auto& node = grid(x, y);
                node.fluid.rho = 1.0;
                node.fluid.u[0] = 0.0;
                node.fluid.u[1] = 0.0;
                BGKEquilibrium<Lattice>::compute(node.fluid, node.df.f);
            }
        }
    }

    // Compute analytical solution at time t
    std::vector<double> analytical_solution(double t) const {
        std::vector<double> u_analytical(ny_);
        const double H = ny_ - 1;

        // Womersley number: α = H * sqrt(ω/ν)
        const double alpha = H * std::sqrt(omega_ / viscosity_);

        // Complex wave number
        const std::complex<double> i(0.0, 1.0);
        const std::complex<double> lambda = std::sqrt(i * omega_ / viscosity_);

        for (int y = 0; y < ny_; ++y) {
            const double y_pos = static_cast<double>(y) - H / 2.0;

            // Complex velocity
            const std::complex<double> u_complex = amplitude_ *
                (1.0 - std::cosh(lambda * y_pos) / std::cosh(lambda * H / 2.0)) *
                std::exp(i * omega_ * t);

            u_analytical[y] = u_complex.real();
        }
        return u_analytical;
    }

    double compute_error(const LatticeGrid<Lattice::D, Lattice::Q>& grid, double t) const {
        auto u_exact = analytical_solution(t);
        double error = 0.0;
        const int x_mid = nx_ / 2;

        for (int y = 0; y < ny_; ++y) {
            const double u_sim = grid(x_mid, y).fluid.u[0];
            const double diff = u_sim - u_exact[y];
            error += diff * diff;
        }
        return std::sqrt(error / ny_);
    }

private:
    int nx_, ny_;
    double amplitude_;
    double frequency_;
    double viscosity_;
    double omega_;
};

/**
 * Hagen-Poiseuille Flow: Pressure-driven flow in a circular pipe
 * Analytical solution: u(r) = (dp/dx) * (R² - r²) / (4μ)
 * where R is pipe radius, r is radial distance from center, μ is dynamic viscosity
 */
template<typename Lattice>
class HagenPoiseuilleFlow {
public:
    HagenPoiseuilleFlow(int nx, int ny, double dp_dx, double viscosity)
        : nx_(nx), ny_(ny), dp_dx_(dp_dx), viscosity_(viscosity) {
        // Compute effective radius (inscribed circle in domain)
        const double cx = nx / 2.0;
        const double cy = ny / 2.0;
        radius_ = std::min(cx, cy) - 1.0;
    }

    // Initialize from rest (or with analytical profile if use_analytical=true)
    void initialize(LatticeGrid<Lattice::D, Lattice::Q>& grid, bool use_analytical = false) const {
        const double cx = nx_ / 2.0;
        const double cy = ny_ / 2.0;

        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                auto& node = grid(x, y);

                // Check if inside pipe
                const double dx = x - cx;
                const double dy = y - cy;
                const double r = std::sqrt(dx * dx + dy * dy);

                if (r <= radius_) {
                    // Fluid node
                    node.fluid.rho = 1.0;
                    if (use_analytical) {
                        node.fluid.u[0] = analytical_velocity(r);
                    } else {
                        node.fluid.u[0] = 0.0;
                    }
                    node.fluid.u[1] = 0.0;
                    node.is_solid = false;

                    // Initialize distributions with equilibrium
                    BGKEquilibrium<Lattice>::compute(node.fluid, node.df.f);
                } else {
                    // Solid node (wall)
                    node.is_solid = true;
                    node.fluid.rho = 1.0;
                    node.fluid.u[0] = 0.0;
                    node.fluid.u[1] = 0.0;
                    // Initialize solid node distributions to zero (will be handled by bounce-back)
                    for (int i = 0; i < Lattice::Q; ++i) {
                        node.df.f[i] = 0.0;
                    }
                }
            }
        }
    }

    // Compute analytical solution
    double analytical_velocity(double r) const {
        const double rho = 1.0;
        if (r <= radius_) {
            return -dp_dx_ * (radius_ * radius_ - r * r) / (4.0 * rho * viscosity_);
        }
        return 0.0;
    }

    // Compute L2 error (radial profile)
    double compute_error(const LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        const double cx = nx_ / 2.0;
        const double cy = ny_ / 2.0;

        double error = 0.0;
        int count = 0;

        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const double dx = x - cx;
                const double dy = y - cy;
                const double r = std::sqrt(dx * dx + dy * dy);

                if (r <= radius_ && !grid(x, y).is_solid) {
                    const double u_sim = grid(x, y).fluid.u[0];
                    const double u_exact = analytical_velocity(r);
                    const double diff = u_sim - u_exact;
                    error += diff * diff;
                    count++;
                }
            }
        }

        return count > 0 ? std::sqrt(error / count) : 0.0;
    }

private:
    int nx_, ny_;
    double dp_dx_;
    double viscosity_;
    double radius_;
};

/**
 * Stokes Flow: Low Reynolds number flow (creeping flow)
 * Simple shear flow validation: du/dy = const
 * Analytical solution: u(y) = γ * y, where γ is shear rate
 */
template<typename Lattice>
class StokesShearFlow {
public:
    StokesShearFlow(int nx, int ny, double shear_rate, double viscosity)
        : nx_(nx), ny_(ny), shear_rate_(shear_rate), viscosity_(viscosity) {}

    void initialize(LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                auto& node = grid(x, y);

                // Linear velocity profile
                const double u_x = shear_rate_ * static_cast<double>(y);

                node.fluid.rho = 1.0;
                node.fluid.u[0] = u_x;
                node.fluid.u[1] = 0.0;

                BGKEquilibrium<Lattice>::compute(node.fluid, node.df.f);
            }
        }
    }

    std::vector<double> analytical_solution() const {
        std::vector<double> u_analytical(ny_);
        for (int y = 0; y < ny_; ++y) {
            u_analytical[y] = shear_rate_ * static_cast<double>(y);
        }
        return u_analytical;
    }

    // Verify constant shear stress: τ = μ * du/dy = μ * γ
    double compute_shear_stress(const LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        const int x_mid = nx_ / 2;
        const double rho = 1.0;

        // Compute du/dy numerically
        double du_dy = 0.0;
        int count = 0;

        for (int y = 1; y < ny_ - 1; ++y) {
            const double u_above = grid(x_mid, y + 1).fluid.u[0];
            const double u_below = grid(x_mid, y - 1).fluid.u[0];
            du_dy += (u_above - u_below) / 2.0;
            count++;
        }

        du_dy /= count;
        return rho * viscosity_ * du_dy;
    }

    double theoretical_shear_stress() const {
        return viscosity_; // μ * γ with ρ = 1
    }

    double compute_error(const LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        auto u_exact = analytical_solution();
        double error = 0.0;
        const int x_mid = nx_ / 2;

        for (int y = 0; y < ny_; ++y) {
            const double u_sim = grid(x_mid, y).fluid.u[0];
            const double diff = u_sim - u_exact[y];
            error += diff * diff;
        }
        return std::sqrt(error / ny_);
    }

private:
    int nx_, ny_;
    double shear_rate_;
    double viscosity_;
};

/**
 * Kolmogorov Flow: 2D sinusoidal forcing pattern
 * Used to study transition to turbulence
 * Analytical solution (steady state): u(y) = (F/k²ν) * sin(ky)
 * where F is forcing amplitude, k is wave number, ν is viscosity
 */
template<typename Lattice>
class KolmogorovFlow {
public:
    KolmogorovFlow(int nx, int ny, double force_amplitude, double viscosity, double k = 1.0)
        : nx_(nx), ny_(ny), force_amplitude_(force_amplitude),
          viscosity_(viscosity), k_(k) {}

    void initialize(LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                auto& node = grid(x, y);
                node.fluid.rho = 1.0;
                node.fluid.u[0] = 0.0;
                node.fluid.u[1] = 0.0;
                BGKEquilibrium<Lattice>::compute(node.fluid, node.df.f);
            }
        }
    }

    std::vector<double> analytical_solution() const {
        std::vector<double> u_analytical(ny_);
        const double L = 2.0 * M_PI / k_;

        for (int y = 0; y < ny_; ++y) {
            const double y_pos = L * y / ny_;
            u_analytical[y] = (force_amplitude_ / (k_ * k_ * viscosity_)) *
                              std::sin(k_ * y_pos);
        }
        return u_analytical;
    }

    // Get forcing term for y position
    double get_force(int y) const {
        const double L = 2.0 * M_PI / k_;
        const double y_pos = L * y / ny_;
        return force_amplitude_ * std::sin(k_ * y_pos);
    }

    double compute_error(const LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        auto u_exact = analytical_solution();
        double error = 0.0;
        const int x_mid = nx_ / 2;

        for (int y = 0; y < ny_; ++y) {
            const double u_sim = grid(x_mid, y).fluid.u[0];
            const double diff = u_sim - u_exact[y];
            error += diff * diff;
        }
        return std::sqrt(error / ny_);
    }

private:
    int nx_, ny_;
    double force_amplitude_;
    double viscosity_;
    double k_;
};

} // namespace validation
} // namespace elbm

#endif // ANALYTICAL_CASES_H
