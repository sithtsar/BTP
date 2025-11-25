#ifndef BGK_SOLVER_H
#define BGK_SOLVER_H

#include "core/lattice.h"
#include "core/fluid_state.h"
#include "solvers/equilibrium.h"
#include <cmath>

namespace elbm {

template<typename Lattice>
class BGKSolver {
public:
    BGKSolver(double viscosity, double dt = 1.0, double dx = 1.0)
        : dt_(dt), dx_(dx) {
        // Compute relaxation time from viscosity
        // ν = cs² (τ - Δt/2)
        tau_ = viscosity / Lattice::cs2 + 0.5 * dt_;
        omega_ = dt_ / tau_;

        // Check stability condition: Δt/τ > 1/2
        if (omega_ <= 0.5) {
            throw std::runtime_error("Stability condition violated: Δt/τ must be > 0.5");
        }


    }

    // Collision step with optional external force
    void collide(Node<Lattice::D, Lattice::Q>& node,
                 const std::array<double, Lattice::D>& force = {}) {
        auto& df = node.df;
        auto& fluid = node.fluid;

        // Compute macroscopic variables
        computeMacroscopic<Lattice>(df.f, fluid);

        // Compute equilibrium distribution
        BGKEquilibrium<Lattice>::compute(fluid, df.f_eq);

        // BGK collision
        for (int i = 0; i < Lattice::Q; ++i) {
            df.f_new[i] = df.f[i] - omega_ * (df.f[i] - df.f_eq[i]);
        }

        // Add body force term if present (Guo forcing scheme)
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
            for (int i = 0; i < Lattice::Q; ++i) {
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

                double u_dot_F = 0.0;
                for (int d = 0; d < Lattice::D; ++d) {
                    u_dot_F += fluid.u[d] * force[d];
                }

                double force_term = Lattice::weight(i) * (1.0 - 0.5 * omega_) *
                                   (ci_dot_F / Lattice::cs2 +
                                    (ci_dot_u * ci_dot_F) / (Lattice::cs2 * Lattice::cs2) -
                                    u_dot_F / Lattice::cs2);

                df.f_new[i] += force_term;
            }
        }
    }

    // Compute H-function for entropy monitoring (needed by EntropyMonitor)
    double computeH(const std::array<double, Lattice::Q>& f) const {
        constexpr double eps = 1e-12;
        double H = 0.0;

        for (int i = 0; i < Lattice::Q; ++i) {
            if (f[i] > eps) {
                const double log_term = std::log(f[i] / Lattice::weight(i));
                if (std::isfinite(log_term)) {
                    H += f[i] * log_term;
                }
            }
        }

        return std::isfinite(H) ? H : 0.0;
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

    // Get viscosity
    double viscosity() const {
        return Lattice::cs2 * (tau_ - 0.5 * dt_);
    }

    // Get relaxation parameters
    double tau() const { return tau_; }
    double omega() const { return omega_; }

private:
    double dt_;                               // Time step
    double dx_;                               // Lattice spacing
    double tau_;                              // Relaxation time
    double omega_;                            // Collision frequency (1/tau)
};

} // namespace elbm

#endif // BGK_SOLVER_H
