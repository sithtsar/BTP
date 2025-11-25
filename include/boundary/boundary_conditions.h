#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H

#include "core/lattice.h"
#include "core/fluid_state.h"
#include <vector>

namespace elbm {

// Boundary condition types
enum class BCType {
    PERIODIC,
    BOUNCE_BACK,
    PRESSURE,
    VELOCITY,
    OPEN
};

// Bounce-back boundary condition (no-slip wall)
template<typename Lattice>
class BounceBackBC {
public:
    static void apply(Node<Lattice::D, Lattice::Q>& node) {
        // Swap populations in opposite directions
        auto& f = node.df.f;
        auto f_temp = f;

        for (int i = 0; i < Lattice::Q; ++i) {
            const int opp = Lattice::opposite[i];
            f[i] = f_temp[opp];
        }
    }
};

// Pressure boundary condition (Zou-He scheme)
template<typename Lattice>
class PressureBC {
public:
    static void applyLeft(Node<Lattice::D, Lattice::Q>& node, double rho_target) {
        // For D2Q9 at left boundary (x = 0)
        // Unknown distributions: f1, f5, f8 (pointing right)
        // Known distributions: all others

        auto& f = node.df.f;
        const auto& u = node.fluid.u;

        if constexpr (Lattice::D == 2) {
            // Compute velocity from known distributions
            const double rho = rho_target;
            const double uy = (f[2] + f[5] + f[6] - f[4] - f[7] - f[8]) / rho;
            const double ux = 1.0 - (f[0] + f[2] + f[4] + 2.0 * (f[3] + f[6] + f[7])) / rho;

            // Set unknown populations
            f[1] = f[3] + (2.0 / 3.0) * rho * ux;
            f[5] = f[7] - 0.5 * (f[2] - f[4]) + 0.5 * rho * uy + (1.0 / 6.0) * rho * ux;
            f[8] = f[6] + 0.5 * (f[2] - f[4]) - 0.5 * rho * uy + (1.0 / 6.0) * rho * ux;
        }
    }

    static void applyRight(Node<Lattice::D, Lattice::Q>& node, double rho_target) {
        // For D2Q9 at right boundary (x = nx-1)
        // Unknown distributions: f3, f6, f7 (pointing left)

        auto& f = node.df.f;

        if constexpr (Lattice::D == 2) {
            const double rho = rho_target;
            const double uy = (f[2] + f[5] + f[6] - f[4] - f[7] - f[8]) / rho;
            const double ux = -1.0 + (f[0] + f[2] + f[4] + 2.0 * (f[1] + f[5] + f[8])) / rho;

            // Set unknown populations
            f[3] = f[1] - (2.0 / 3.0) * rho * ux;
            f[6] = f[8] - 0.5 * (f[2] - f[4]) + 0.5 * rho * uy - (1.0 / 6.0) * rho * ux;
            f[7] = f[5] + 0.5 * (f[2] - f[4]) - 0.5 * rho * uy - (1.0 / 6.0) * rho * ux;
        }
    }
};

// Velocity boundary condition (Zou-He scheme)
template<typename Lattice>
class VelocityBC {
public:
    static void applyLeft(Node<Lattice::D, Lattice::Q>& node,
                         const std::array<double, Lattice::D>& u_target) {
        // For D2Q9 at left boundary (x = 0)
        // Known velocity, compute density and unknown populations
        auto& f = node.df.f;

        if constexpr (Lattice::D == 2) {
            const double ux = u_target[0];
            const double uy = u_target[1];

            // Compute density from known populations
            const double rho = (f[0] + f[2] + f[4] + 2.0 * (f[3] + f[6] + f[7])) / (1.0 - ux);

            // Set unknown populations (f1, f5, f8 pointing right)
            f[1] = f[3] + (2.0 / 3.0) * rho * ux;
            f[5] = f[7] - 0.5 * (f[2] - f[4]) + 0.5 * rho * uy + (1.0 / 6.0) * rho * ux;
            f[8] = f[6] + 0.5 * (f[2] - f[4]) - 0.5 * rho * uy + (1.0 / 6.0) * rho * ux;
        }
    }

    static void applyTop(Node<Lattice::D, Lattice::Q>& node,
                        const std::array<double, Lattice::D>& u_target) {
        // For D2Q9 at top boundary (y = ny-1)
        // Unknown distributions: f4, f7, f8 (pointing down)

        auto& f = node.df.f;

        if constexpr (Lattice::D == 2) {
            const double ux = u_target[0];
            const double uy = u_target[1];

            // Compute density
            const double rho = (f[0] + f[1] + f[3] + 2.0 * (f[2] + f[5] + f[6])) /
                              (1.0 + uy);

            // Set unknown populations
            f[4] = f[2] - (2.0 / 3.0) * rho * uy;
            f[7] = f[5] + 0.5 * (f[1] - f[3]) - 0.5 * rho * ux - (1.0 / 6.0) * rho * uy;
            f[8] = f[6] - 0.5 * (f[1] - f[3]) + 0.5 * rho * ux - (1.0 / 6.0) * rho * uy;
        }
    }

    static void applyBottom(Node<Lattice::D, Lattice::Q>& node,
                           const std::array<double, Lattice::D>& u_target) {
        // For D2Q9 at bottom boundary (y = 0)
        // Unknown distributions: f2, f5, f6 (pointing up)

        auto& f = node.df.f;

        if constexpr (Lattice::D == 2) {
            const double ux = u_target[0];
            const double uy = u_target[1];

            // Compute density
            const double rho = (f[0] + f[1] + f[3] + 2.0 * (f[4] + f[7] + f[8])) /
                              (1.0 - uy);

            // Set unknown populations
            f[2] = f[4] + (2.0 / 3.0) * rho * uy;
            f[5] = f[7] - 0.5 * (f[1] - f[3]) + 0.5 * rho * ux + (1.0 / 6.0) * rho * uy;
            f[6] = f[8] + 0.5 * (f[1] - f[3]) - 0.5 * rho * ux + (1.0 / 6.0) * rho * uy;
        }
    }
};

// Boundary condition manager
template<typename Lattice>
class BoundaryManager {
public:
    void applyBoundaries(LatticeGrid<Lattice::D, Lattice::Q>& grid) {
        const int nx = grid.nx();
        const int ny = grid.ny();

        // Apply left boundary (inlet: velocity or pressure) - exclude corners
        if (left_bc_type_ == BCType::PRESSURE) {
            for (int y = 1; y < ny - 1; ++y) {
                PressureBC<Lattice>::applyLeft(grid(0, y), left_rho_);
            }
        } else if (left_bc_type_ == BCType::VELOCITY) {
            for (int y = 1; y < ny - 1; ++y) {
                VelocityBC<Lattice>::applyLeft(grid(0, y), left_velocity_);
            }
        }

        // Apply right boundary (outlet: pressure or open) - exclude corners
        if (right_bc_type_ == BCType::PRESSURE) {
            for (int y = 1; y < ny - 1; ++y) {
                PressureBC<Lattice>::applyRight(grid(nx - 1, y), right_rho_);
            }
        }

        // Apply top boundary (wall or velocity) - exclude corners
        if (top_bc_type_ == BCType::BOUNCE_BACK) {
            for (int x = 1; x < nx - 1; ++x) {
                BounceBackBC<Lattice>::apply(grid(x, ny - 1));
            }
        } else if (top_bc_type_ == BCType::VELOCITY) {
            for (int x = 1; x < nx - 1; ++x) {
                VelocityBC<Lattice>::applyTop(grid(x, ny - 1), top_velocity_);
            }
        }

        // Apply bottom boundary (wall or velocity) - exclude corners
        if (bottom_bc_type_ == BCType::BOUNCE_BACK) {
            for (int x = 1; x < nx - 1; ++x) {
                BounceBackBC<Lattice>::apply(grid(x, 0));
            }
        } else if (bottom_bc_type_ == BCType::VELOCITY) {
            for (int x = 1; x < nx - 1; ++x) {
                VelocityBC<Lattice>::applyBottom(grid(x, 0), bottom_velocity_);
            }
        }

        // Handle corner nodes explicitly
        // Bottom-left corner (0, 0)
        if (bottom_bc_type_ == BCType::BOUNCE_BACK || left_bc_type_ == BCType::VELOCITY) {
            BounceBackBC<Lattice>::apply(grid(0, 0));
        }
        // Bottom-right corner (nx-1, 0)
        if (bottom_bc_type_ == BCType::BOUNCE_BACK || right_bc_type_ == BCType::PRESSURE) {
            BounceBackBC<Lattice>::apply(grid(nx - 1, 0));
        }
        // Top-left corner (0, ny-1)
        if (top_bc_type_ == BCType::BOUNCE_BACK || left_bc_type_ == BCType::VELOCITY) {
            BounceBackBC<Lattice>::apply(grid(0, ny - 1));
        }
        // Top-right corner (nx-1, ny-1)
        if (top_bc_type_ == BCType::BOUNCE_BACK || right_bc_type_ == BCType::PRESSURE) {
            BounceBackBC<Lattice>::apply(grid(nx - 1, ny - 1));
        }
    }

    // Setters for boundary conditions
    void setLeftBC(BCType type, double rho = 1.0) {
        left_bc_type_ = type;
        left_rho_ = rho;
    }

    void setLeftBC(BCType type, const std::array<double, Lattice::D>& u) {
        left_bc_type_ = type;
        left_velocity_ = u;
    }

    void setRightBC(BCType type, double rho = 1.0) {
        right_bc_type_ = type;
        right_rho_ = rho;
    }

    void setTopBC(BCType type, const std::array<double, Lattice::D>& u = {}) {
        top_bc_type_ = type;
        top_velocity_ = u;
    }

    void setBottomBC(BCType type, const std::array<double, Lattice::D>& u = {}) {
        bottom_bc_type_ = type;
        bottom_velocity_ = u;
    }

private:
    BCType left_bc_type_ = BCType::PERIODIC;
    BCType right_bc_type_ = BCType::PERIODIC;
    BCType top_bc_type_ = BCType::PERIODIC;
    BCType bottom_bc_type_ = BCType::PERIODIC;

    double left_rho_ = 1.0;
    double right_rho_ = 1.0;
    std::array<double, Lattice::D> left_velocity_ = {};
    std::array<double, Lattice::D> top_velocity_ = {};
    std::array<double, Lattice::D> bottom_velocity_ = {};
};

} // namespace elbm

#endif // BOUNDARY_CONDITIONS_H
