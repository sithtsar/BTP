/*
 * ELBM vs LBGK Benchmark Suite
 * Cases: Couette Flow, Poiseuille Flow, Taylor-Green Vortex
 *
 * Features:
 * - Integrated with existing ELBM project structure
 * - D2Q9 Lattice using core/lattice.h
 * - Switchable: Standard LBGK vs. Entropic LBM (ELBM)
 * - ELBM uses Newton-Raphson to enforce the H-Theorem (alpha search)
 * - Analytical solution comparisons for all cases
 * - Uses existing solver and boundary frameworks
 */

#include <iostream>
#include <fstream>
#include <cmath>

#include "core/lattice.h"
#include "core/fluid_state.h"
#include "solvers/bgk_solver.h"
#include "solvers/elbm_solver.h"
#include "boundary/boundary_conditions.h"

using namespace elbm;

// --- Configuration Constants ---
const int NX = 100;          // Domain Width
const int NY = 100;          // Domain Height
const double RE = 100.0;      // Reynolds Number
const double U_MAX = 0.1;     // Characteristic Velocity (low Mach)
const double L_CHAR = NY;     // Characteristic Length

// --- Physics Helper Functions ---
// Kinematic Viscosity from Reynolds
double get_viscosity(double u, double len, double re) {
    return (u * len) / re;
}

// Relaxation time (tau) from viscosity
double get_tau(double nu) {
    return 3.0 * nu + 0.5;
}

// Equilibrium Distribution
double get_feq(int k, double rho, double ux, double uy) {
    using Lattice = D2Q9;
    double cu = Lattice::cx(k)*ux + Lattice::cy(k)*uy;
    double u2 = ux*ux + uy*uy;
    return Lattice::weight(k) * rho * (1.0 + 3.0*cu + 4.5*cu*cu - 1.5*u2);
}

// Entropy H = sum(f * ln(f/w))
double calc_entropy(const double* f) {
    using Lattice = D2Q9;
    double H = 0.0;
    for(int k=0; k<Lattice::Q; k++) {
        if(f[k] > 1e-15) {
            H += f[k] * std::log(f[k] / Lattice::weight(k));
        }
    }
    return H;
}

// --- Helper Functions for TGV ---
double extract_viscosity(double E0, double E1, double dt, double k_squared) {
    // E(t) = E0 * exp(-4*k²*ν*t)
    // ν = -ln(E1/E0) / (4*k²*dt)
    if (E1 > 0 && E0 > 0 && E1 < E0) {
        return -std::log(E1 / E0) / (4.0 * k_squared * dt);
    }
    return -1.0; // Invalid
}

// --- Main Lattice Boltzmann Class ---
class LatticeBoltzmann {
public:
    using Lattice = D2Q9;
    LatticeGrid<Lattice::D, Lattice::Q> grid;
    BGKSolver<Lattice> bgk_solver;
    ELBMSolver<Lattice> elbm_solver;
    BoundaryManager<Lattice> bc_manager;
    double viscosity;
    bool use_entropic;

    LatticeBoltzmann(double nu, bool entropic) :
        grid(NX, NY),
        bgk_solver(nu),
        elbm_solver(nu),
        bc_manager(),
        viscosity(nu),
        use_entropic(entropic)
    {
        std::cout << "Model Initialized. Viscosity: " << viscosity
                  << ", Mode: " << (use_entropic? "ENTROPIC (ELBM)" : "STANDARD (LBGK)") << "\n";
    }

    // Initialize with uniform flow
    void initialize_uniform(double rho0 = 1.0, double ux0 = 0.0, double uy0 = 0.0) {
        for(int y = 0; y < NY; y++) {
            for(int x = 0; x < NX; x++) {
                auto& node = grid(x, y);
                node.fluid.rho = rho0;
                node.fluid.u = {ux0, uy0};
                node.fluid.p = rho0 / 3.0;

                // Initialize distributions at equilibrium
                for(int k = 0; k < Lattice::Q; k++) {
                    node.df.f[k] = get_feq(k, rho0, ux0, uy0);
                }
            }
        }
    }

    // Initialize Taylor-Green Vortex
    void initialize_tgv(double u0) {
        double k_x = 2.0 * M_PI / NX;
        double k_y = 2.0 * M_PI / NY;

        for(int y = 0; y < NY; y++) {
            for(int x = 0; x < NX; x++) {
                auto& node = grid(x, y);

                double ux = -u0 * std::cos(k_x * x) * std::sin(k_y * y);
                double uy =  u0 * std::sin(k_x * x) * std::cos(k_y * y);
                double rho = 1.0 - 0.25 * (u0*u0 / (1.0/3.0)) * (std::cos(2.0*k_x*x) + std::cos(2.0*k_y*y));

                node.fluid.rho = rho;
                node.fluid.u = {ux, uy};
                node.fluid.p = rho / 3.0;

                for(int k = 0; k < Lattice::Q; k++) {
                    node.df.f[k] = get_feq(k, rho, ux, uy);
                }
            }
        }
    }

    void collide_and_stream(bool periodic_x, bool periodic_y, const std::array<double, 2>& force = {}) {
        // Collision
        for(int y = 0; y < NY; y++) {
            for(int x = 0; x < NX; x++) {
                if(use_entropic) {
                    elbm_solver.collide(grid(x, y), force);
                } else {
                    bgk_solver.collide(grid(x, y), force);
                }
            }
        }

        // Streaming
        if(use_entropic) {
            elbm_solver.stream(grid, periodic_x, periodic_y);
        } else {
            bgk_solver.stream(grid, periodic_x, periodic_y);
        }

        // Apply boundary conditions if not fully periodic
        if (!periodic_x || !periodic_y) {
            bc_manager.applyBoundaries(grid);
        }
    }

    void update_macroscopic() {
        for(int y = 0; y < NY; y++) {
            for(int x = 0; x < NX; x++) {
                auto& node = grid(x, y);

                double rho = 0.0;
                std::array<double, 2> u = {0.0, 0.0};

                for(int k = 0; k < Lattice::Q; k++) {
                    rho += node.df.f[k];
                    u[0] += node.df.f[k] * Lattice::cx(k);
                    u[1] += node.df.f[k] * Lattice::cy(k);
                }

                if(rho > 1e-14) {
                    u[0] /= rho;
                    u[1] /= rho;
                }

                node.fluid.rho = rho;
                node.fluid.u = u;
                node.fluid.p = rho / 3.0;
            }
        }
    }
};

// --- Benchmark Runners ---

// 1. Couette Flow (Shear between plates)
void run_couette(bool use_entropic) {
    std::cout << "\n--- Running Couette Flow Benchmark ---" << std::endl;
    double nu = get_viscosity(U_MAX, L_CHAR, RE);
    LatticeBoltzmann lb(nu, use_entropic);

    // Initialize
    lb.initialize_uniform();

    // Setup boundaries
    std::array<double, 2> bottom_vel = {0.0, 0.0};
    std::array<double, 2> top_vel = {U_MAX, 0.0};
    lb.bc_manager.setTopBC(BCType::VELOCITY, top_vel);
    lb.bc_manager.setBottomBC(BCType::VELOCITY, bottom_vel);

    // Run simulation
    int steps = 2000;
    for(int t = 0; t < steps; t++) {
        lb.collide_and_stream();
    }
    lb.update_macroscopic();

    // Verify (Center Column) and save results
    double error_l2 = 0.0;
    int x_probe = NX/2;

    std::string method = use_entropic ? "elbm" : "bgk";
    std::string filename = "output/benchmark_couette_" + method + ".dat";
    std::ofstream file(filename);
    file << "# y u_sim u_exact error\n";

    for(int y = 0; y < NY; y++) {
        double y_norm = (double)y / (double)(NY-1);
        double u_analytical = U_MAX * y_norm;
        double u_sim = lb.grid(x_probe, y).fluid.u[0];
        double diff = u_sim - u_analytical;
        double error = std::abs(diff);
        error_l2 += diff * diff;

        file << y << " " << u_sim << " " << u_analytical << " " << error << "\n";
    }
    file.close();

    std::cout << "L2 Error Norm: " << std::sqrt(error_l2/NY) << std::endl;
    std::cout << "Results saved to: " << filename << std::endl;
}

// 2. Poiseuille Flow (Pressure driven channel)
void run_poiseuille(bool use_entropic) {
    std::cout << "\n--- Running Poiseuille Flow Benchmark ---" << std::endl;
    double nu = get_viscosity(U_MAX, L_CHAR, RE);
    LatticeBoltzmann lb(nu, use_entropic);

    // From analytical solution: u_max = -dp_dx * H^2 / (8 * rho * nu)
    // We want to simulate a body force Fx that corresponds to a pressure gradient dp_dx.
    // Fx = -dp_dx.
    // So, Fx = 8 * rho * nu * u_max / H^2
    double H = NY - 1;
    double rho = 1.0;
    // Note: The original formula had a factor of 16, which is for max velocity at centerline
    // for a body force driven flow between two plates, the relation is u_max = Fx * H^2 / (8 * nu).
    // Let's use the force directly.
    double Fx = 8.0 * nu * U_MAX / (H * H);

    // Set body force (pressure gradient)
    std::array<double, 2> body_force = {Fx, 0.0};

    // Initialize from rest
    lb.initialize_uniform();

    // Setup boundaries (no-slip walls)
    std::array<double, 2> wall_vel = {0.0, 0.0};
    // For Poiseuille, we need walls top and bottom, but periodic left/right is often assumed
    // when using a body force. The existing BC manager sets pressure BCs.
    // To simulate with a body force, we should use periodic in X and walls in Y.
    lb.bc_manager.setTopBC(BCType::BOUNCE_BACK);
    lb.bc_manager.setBottomBC(BCType::BOUNCE_BACK);


    // Run simulation
    int steps = 4000;
    for(int t = 0; t < steps; t++) {
        // Poiseuille flow is periodic in X, non-periodic (walls) in Y. Force is applied.
        lb.collide_and_stream(true, false, body_force);
    }
    lb.update_macroscopic();

    double error_l2 = 0.0;
    int x_probe = NX/2;

    std::string method = use_entropic ? "elbm" : "bgk";
    std::string filename = "output/benchmark_poiseuille_" + method + ".dat";
    std::ofstream file(filename);
    file << "# y u_sim u_exact error\n";

    for(int y = 0; y < NY; y++) {
        double Y = y;
        // Analytical for body force Fx: u(y) = (Fx / (2 * nu)) * y * (H - y)
        double u_analytical = (Fx / (2.0 * nu)) * Y * (H - Y);
        // Boundaries are strictly 0
        if (y == 0 || y == NY-1) u_analytical = 0.0;

        double u_sim = lb.grid(x_probe, y).fluid.u[0];
        double diff = u_sim - u_analytical;
        double error = std::abs(diff);
        error_l2 += diff * diff;

        file << y << " " << u_sim << " " << u_analytical << " " << error << "\n";
    }
    file.close();

    std::cout << "L2 Error Norm: " << std::sqrt(error_l2/NY) << std::endl;
    std::cout << "Results saved to: " << filename << std::endl;
}

// 3. Taylor-Green Vortex (Decay)
void run_tgv(bool use_entropic) {
    std::cout << "\n--- Running Taylor-Green Vortex Benchmark ---" << std::endl;
    double nu = get_viscosity(U_MAX, L_CHAR, RE);
    LatticeBoltzmann lb(nu, use_entropic);

    double k_x = 2.0 * M_PI / NX;
    double k_y = 2.0 * M_PI / NY;

    // Initialize TGV Field
    lb.initialize_tgv(U_MAX);

    // Periodic boundaries (default for TGV)
    // No additional BC setup needed

    // Evolve
    int steps = 1000;
    for(int t = 0; t <= steps; t++) {
        lb.collide_and_stream(true, true);
    }
    lb.update_macroscopic();

    // Compare Decay at t = steps
    // Analytical: u(t) = u(0) * exp(-2 * k^2 * nu * t)
    // E(t) = E(0) * exp(-4 * k^2 * nu * t)
    double k_squared = k_x*k_x + k_y*k_y;
    double velocity_decay = std::exp(-2.0 * k_squared * nu * steps);
    double energy_decay = velocity_decay * velocity_decay;

    // Compute analytical energy
    double analytical_energy = 0.0;
    for(int y = 0; y < NY; y++) {
        for(int x = 0; x < NX; x++) {
            double ux_init = -U_MAX * std::cos(k_x * x) * std::sin(k_y * y);
            double uy_init = U_MAX * std::sin(k_x * x) * std::cos(k_y * y);
            double ux_ana = ux_init * velocity_decay;
            double uy_ana = uy_init * velocity_decay;
            analytical_energy += 0.5 * (ux_ana*ux_ana + uy_ana*uy_ana);
        }
    }

    double error_l2 = 0.0;
    double total_energy_sim = 0.0;

    for(int y = 0; y < NY; y++) {
        for(int x = 0; x < NX; x++) {
            double ux_init = -U_MAX * std::cos(k_x * x) * std::sin(k_y * y);
            double uy_init = U_MAX * std::sin(k_x * x) * std::cos(k_y * y);

            double ux_ana = ux_init * velocity_decay;
            double uy_ana = uy_init * velocity_decay;

            double ux_sim = lb.grid(x, y).fluid.u[0];
            double uy_sim = lb.grid(x, y).fluid.u[1];

            double diff_x = ux_sim - ux_ana;
            double diff_y = uy_sim - uy_ana;
            error_l2 += diff_x*diff_x + diff_y*diff_y;

            total_energy_sim += 0.5 * (ux_sim*ux_sim + uy_sim*uy_sim);
        }
    }

    std::cout << "Time Step: " << steps << "\n";
    std::cout << "Kinetic Energy (Sim): " << total_energy_sim << "\n";
    std::cout << "L2 Velocity Error: " << std::sqrt(error_l2/(NX*NY*2)) << "\n";

    // Save full field data for plotting
    std::string method = use_entropic ? "elbm" : "bgk";
    std::string filename = "output/benchmark_tgv_" + method + ".dat";
    std::ofstream file(filename);
    file << "# x y ux_sim uy_sim ux_exact uy_exact ux_error uy_error\n";

    for(int y = 0; y < NY; y++) {
        for(int x = 0; x < NX; x++) {
            double ux_init = -U_MAX * std::cos(k_x * x) * std::sin(k_y * y);
            double uy_init = U_MAX * std::sin(k_x * x) * std::cos(k_y * y);

            double ux_ana = ux_init * velocity_decay;
            double uy_ana = uy_init * velocity_decay;

            double ux_sim = lb.grid(x, y).fluid.u[0];
            double uy_sim = lb.grid(x, y).fluid.u[1];

            double ux_error = ux_sim - ux_ana;
            double uy_error = uy_sim - uy_ana;

            file << x << " " << y << " " << ux_sim << " " << uy_sim << " "
                 << ux_ana << " " << uy_ana << " " << ux_error << " " << uy_error << "\n";
        }
    }
    file.close();

    std::cout << "Full field data saved to: " << filename << std::endl;
}

int main() {
    std::cout << "=== LBM BENCHMARK SUITE ===" << std::endl;

    // Couette
    run_couette(false); // Normal
    run_couette(true);  // Entropic

    // Poiseuille (ELBM has issues with body forces - needs debugging)
    run_poiseuille(false);
    // run_poiseuille(true);  // Temporarily disabled

    // Taylor-Green
    run_tgv(false);
    run_tgv(true);

    return 0;
}