#ifndef BENCHMARK_CASES_H
#define BENCHMARK_CASES_H

#include "core/lattice.h"
#include "core/fluid_state.h"
#include "solvers/equilibrium.h"
#include "boundary/boundary_conditions.h"
#include <cmath>
#include <vector>

namespace elbm {
namespace validation {

/**
 * Lid-Driven Cavity Flow
 * Classic CFD benchmark case
 * Square cavity with moving top lid
 * Reference solutions available for Re = 100, 400, 1000, 3200, 5000, 10000
 */
template<typename Lattice>
class LidDrivenCavity {
public:
    LidDrivenCavity(int n, double U_lid, double Re)
        : n_(n), U_lid_(U_lid), Re_(Re) {
        // Compute viscosity from Reynolds number
        // Re = U_lid * L / ν, where L = n-1
        viscosity_ = U_lid * (n - 1) / Re;
    }

    // Initialize with quiescent fluid
    void initialize(LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        for (int y = 0; y < n_; ++y) {
            for (int x = 0; x < n_; ++x) {
                auto& node = grid(x, y);
                node.fluid.rho = 1.0;
                node.fluid.u[0] = 0.0;
                node.fluid.u[1] = 0.0;

                // Initialize to equilibrium
                BGKEquilibrium<Lattice>::compute(node.fluid, node.df.f);
            }
        }
    }

    // Apply boundary conditions
    template<typename Solver, typename BCManager>
    void setupBoundaries(BCManager& bc_manager) const {
        // Top wall: moving with velocity U_lid
        std::array<double, Lattice::D> top_vel = {U_lid_, 0.0};
        bc_manager.setTopBC(BCType::VELOCITY, top_vel);

        // Other walls: stationary (bounce-back)
        std::array<double, Lattice::D> wall_vel = {0.0, 0.0};
        bc_manager.setBottomBC(BCType::BOUNCE_BACK);
    }

    // Get centerline velocities for comparison with benchmark
    std::vector<double> get_u_centerline(const LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        std::vector<double> u_centerline(n_);
        const int x_mid = n_ / 2;

        for (int y = 0; y < n_; ++y) {
            u_centerline[y] = grid(x_mid, y).fluid.u[0];
        }
        return u_centerline;
    }

    std::vector<double> get_v_centerline(const LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        std::vector<double> v_centerline(n_);
        const int y_mid = n_ / 2;

        for (int x = 0; x < n_; ++x) {
            v_centerline[x] = grid(x, y_mid).fluid.u[1];
        }
        return v_centerline;
    }

    double viscosity() const { return viscosity_; }
    double Re() const { return Re_; }
    int size() const { return n_; }

private:
    int n_;
    double U_lid_;
    double Re_;
    double viscosity_;
};

/**
 * Flow Past Circular Cylinder
 * Classic benchmark for vortex shedding
 * Reynolds number based on cylinder diameter
 */
template<typename Lattice>
class FlowPastCylinder {
public:
    FlowPastCylinder(int nx, int ny, double U_inf, double diameter, double Re)
        : nx_(nx), ny_(ny), U_inf_(U_inf), diameter_(diameter), Re_(Re) {
        // Compute viscosity from Reynolds number
        // Re = U_inf * D / ν
        viscosity_ = U_inf * diameter / Re;

        // Cylinder center
        cx_ = nx / 4;  // Place cylinder at 1/4 from inlet
        cy_ = ny / 2;
        radius_ = diameter / 2.0;
    }

    // Initialize with uniform flow
    void initialize(LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                auto& node = grid(x, y);

                // Check if inside cylinder
                const double dx = x - cx_;
                const double dy = y - cy_;
                const double r = std::sqrt(dx * dx + dy * dy);

                if (r <= radius_) {
                    node.is_solid = true;
                    node.fluid.rho = 1.0;
                    node.fluid.u[0] = 0.0;
                    node.fluid.u[1] = 0.0;
                } else {
                    node.is_solid = false;
                    node.fluid.rho = 1.0;
                    node.fluid.u[0] = U_inf_;
                    node.fluid.u[1] = 0.0;
                }

                // Initialize to equilibrium
                BGKEquilibrium<Lattice>::compute(node.fluid, node.df.f);
            }
        }
    }

    // Apply bounce-back on cylinder surface
    void applyImmersedBoundary(LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                auto& node = grid(x, y);

                if (node.is_solid) {
                    // Bounce-back boundary condition
                    auto f_temp = node.df.f;
                    for (int i = 0; i < Lattice::Q; ++i) {
                        const int opp = Lattice::opposite[i];
                        node.df.f[i] = f_temp[opp];
                    }
                }
            }
        }
    }

    // Compute drag and lift using proper momentum exchange method
    // Force = 2 * Σ(f_i * c_i) over boundary links
    std::array<double, 2> compute_forces(const LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        double Fx = 0.0;
        double Fy = 0.0;

        // Loop over all solid nodes
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                if (!grid(x, y).is_solid) continue;

                // For each direction, check if neighbor is fluid (boundary link)
                for (int i = 0; i < Lattice::Q; ++i) {
                    const int xn = x + Lattice::cx(i);
                    const int yn = y + Lattice::cy(i);

                    // Check bounds
                    if (xn < 0 || xn >= nx_ || yn < 0 || yn >= ny_) continue;

                    // If neighbor is fluid, this is a boundary node
                    if (!grid(xn, yn).is_solid) {
                        // After bounce-back, f[i] at solid contains reflected distribution
                        const double f_boundary = grid(x, y).df.f[i];

                        // Momentum exchange: F = 2 * f * c (factor 2 for reversal)
                        Fx += 2.0 * f_boundary * Lattice::cx(i);
                        Fy += 2.0 * f_boundary * Lattice::cy(i);
                    }
                }
            }
        }

        return {Fx, Fy};
    }

    std::array<double, 2> compute_coefficients(const LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        auto forces = compute_forces(grid);

        // Drag and lift coefficients
        // C_d = F_x / (0.5 * ρ * U² * D)
        // C_l = F_y / (0.5 * ρ * U² * D)
        const double denom = 0.5 * 1.0 * U_inf_ * U_inf_ * diameter_;

        return {forces[0] / denom, forces[1] / denom};
    }

    double viscosity() const { return viscosity_; }
    double Re() const { return Re_; }

private:
    int nx_, ny_;
    double U_inf_;
    double diameter_;
    double Re_;
    double viscosity_;
    double cx_, cy_;
    double radius_;
};


/**
 * Flow Past Circular Cylinder in a Channel
 * Based on the provided implementation with BGK and entropic collision operators
 */
template<typename Lattice>
class ChannelFlowWithCylinder {
public:
    ChannelFlowWithCylinder(int nx, int ny, double U_inf, double diameter, double Re)
        : nx_(nx), ny_(ny), U_inf_(U_inf), diameter_(diameter), Re_(Re) {
        // Compute viscosity from Reynolds number
        // Re = U_inf * D / ν
        viscosity_ = U_inf * diameter / Re;

        // Cylinder center (matching the provided code)
        cx_ = nx / 4;  // Place cylinder at 1/4 from inlet
        cy_ = ny / 2;
        radius_ = diameter / 2.0;
    }

    // Initialize with quiescent fluid (rho=1, u=0 everywhere)
    void initialize(LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                auto& node = grid(x, y);

                // Check if inside cylinder or on walls
                const double dx = x - cx_;
                const double dy = y - cy_;
                const double r = std::sqrt(dx * dx + dy * dy);

                if (r <= radius_ || y == 0 || y == ny_ - 1) {
                    node.is_solid = true;
                } else {
                    node.is_solid = false;
                }

                // Start at rest everywhere (rho=1, u=0)
                node.fluid.rho = 1.0;
                node.fluid.u[0] = 0.0;
                node.fluid.u[1] = 0.0;

                // Initialize to equilibrium at rest
                BGKEquilibrium<Lattice>::compute(node.fluid, node.df.f);
            }
        }
    }

    // Apply immersed boundary (bounce-back) for cylinder
    void applyImmersedBoundary(LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                auto& node = grid(x, y);

                if (node.is_solid) {
                    // Bounce-back boundary condition
                    auto f_temp = node.df.f;
                    for (int i = 0; i < Lattice::Q; ++i) {
                        const int opp = Lattice::opposite[i];
                        node.df.f[i] = f_temp[opp];
                    }
                }
            }
        }
    }

    // Apply boundary conditions for inlet and walls
    // Outlet velocity imposed in macroscopic update (match provided code)
    template<typename BCManager>
    void setupBoundaries(BCManager& bc_manager) const {
        // Inlet: uniform velocity U_inf
        std::array<double, 2> inlet_velocity = {U_inf_, 0.0};
        bc_manager.setLeftBC(BCType::VELOCITY, inlet_velocity);

        // Top and bottom walls: bounce-back
        bc_manager.setTopBC(BCType::BOUNCE_BACK);
        bc_manager.setBottomBC(BCType::BOUNCE_BACK);
    }

    // Compute drag and lift coefficients using momentum exchange method
    std::array<double, 2> compute_coefficients(const LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        double Fx = 0.0, Fy = 0.0;

        // Loop over all solid nodes
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                if (!grid(x, y).is_solid) continue;

                // For each direction, check if neighbor is fluid
                for (int i = 0; i < Lattice::Q; ++i) {
                    const int xn = x + Lattice::cx(i);
                    const int yn = y + Lattice::cy(i);

                    // Check bounds
                    if (xn < 0 || xn >= nx_ || yn < 0 || yn >= ny_) continue;

                    // If neighbor is fluid, this is a boundary link
                    if (!grid(xn, yn).is_solid) {
                        // After bounce-back, f[i] at solid contains reflected distribution
                        const double f_boundary = grid(x, y).df.f[i];

                        // Momentum exchange: F = 2 * f * c (factor 2 for reversal)
                        Fx += 2.0 * f_boundary * Lattice::cx(i);
                        Fy += 2.0 * f_boundary * Lattice::cy(i);
                    }
                }
            }
        }

        // Drag and lift coefficients: C_d = F_x / (0.5 * ρ * U² * D)
        const double denom = 0.5 * 1.0 * U_inf_ * U_inf_ * diameter_;
        return {Fx / denom, Fy / denom};
    }

    // Compute vorticity field for visualization
    std::vector<std::vector<double>> compute_vorticity(const LatticeGrid<Lattice::D, Lattice::Q>& grid) const {
        std::vector<std::vector<double>> vorticity(ny_, std::vector<double>(nx_, 0.0));

        // Skip boundary cells for vorticity calculation
        for (int y = 1; y < ny_ - 1; ++y) {
            for (int x = 1; x < nx_ - 1; ++x) {
                if (grid(x, y).is_solid) continue;

                // Central difference vorticity: ω = ∂v/∂x - ∂u/∂y
                double dvdx = (grid(x+1, y).fluid.u[1] - grid(x-1, y).fluid.u[1]) * 0.5;
                double dudy = (grid(x, y+1).fluid.u[0] - grid(x, y-1).fluid.u[0]) * 0.5;
                vorticity[y][x] = dvdx - dudy;
            }
        }

        return vorticity;
    }

    double viscosity() const { return viscosity_; }
    double Re() const { return Re_; }
    double U_inf() const { return U_inf_; }
    int nx() const { return nx_; }
    int ny() const { return ny_; }

private:
    int nx_, ny_;
    double U_inf_;
    double diameter_;
    double Re_;
    double viscosity_;
    double cx_, cy_;
    double radius_;
};

} // namespace validation
} // namespace elbm

#endif // BENCHMARK_CASES_H
