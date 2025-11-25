#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <chrono>

#include "core/lattice.h"
#include "core/fluid_state.h"
#include "solvers/bgk_solver.h"
#include "solvers/elbm_solver.h"
#include "boundary/boundary_conditions.h"
#include "core/entropy_monitor.h"

using namespace elbm;

// Simulation parameters structure
struct SimulationParams {
    int nx, ny;                  // Grid dimensions
    double viscosity;            // Kinematic viscosity [m^2/s]
    double delta_p;              // Pressure drop [Pa]
    int max_steps;               // Number of timesteps
    int output_interval;         // Output frequency
    std::string solver_type;     // "BGK" or "ELBM"
    std::string output_dir;      // Output directory
};

// Initialize grid with pressure drop
template<typename Lattice>
void initializePressureDrop(LatticeGrid<Lattice::D, Lattice::Q>& grid,
                           double rho_left, double rho_right) {
    const int nx = grid.nx();
    const int ny = grid.ny();

    // Linear pressure profile from left to right
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            // Linear interpolation
            const double t = static_cast<double>(x) / (nx - 1);
            const double rho = rho_left * (1.0 - t) + rho_right * t;

            auto& node = grid(x, y);
            node.fluid.rho = rho;
            node.fluid.u[0] = 0.0;
            node.fluid.u[1] = 0.0;

            // Initialize distribution functions to equilibrium
            BGKEquilibrium<Lattice>::compute(node.fluid, node.df.f);
            node.df.f_eq = node.df.f;
        }
    }
}

// Write results to file
template<typename Lattice>
void writeResults(const LatticeGrid<Lattice::D, Lattice::Q>& grid,
                 const std::string& filename, int step) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    file << std::scientific << std::setprecision(10);
    file << "# Step: " << step << "\n";
    file << "# x y rho ux uy p\n";

    const int nx = grid.nx();
    const int ny = grid.ny();

    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            const auto& node = grid(x, y);
            const auto& fluid = node.fluid;

            file << x << " " << y << " "
                 << fluid.rho << " "
                 << fluid.u[0] << " "
                 << fluid.u[1] << " "
                 << fluid.p << "\n";
        }
    }

    file.close();
}

// Write centerline pressure profile
template<typename Lattice>
void writeCenterlineProfile(const LatticeGrid<Lattice::D, Lattice::Q>& grid,
                           const std::string& filename) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    file << std::scientific << std::setprecision(10);
    file << "# x p rho ux\n";

    const int nx = grid.nx();
    const int ny = grid.ny();
    const int y_center = ny / 2;

    for (int x = 0; x < nx; ++x) {
        const auto& fluid = grid(x, y_center).fluid;
        file << x << " " << fluid.p << " " << fluid.rho << " " << fluid.u[0] << "\n";
    }

    file.close();
}

// Write velocity profile across channel width at centerline
template<typename Lattice>
void writeVelocityProfile(const LatticeGrid<Lattice::D, Lattice::Q>& grid,
                         const std::string& filename) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    file << std::scientific << std::setprecision(10);
    file << "# y ux uy\n";

    const int nx = grid.nx();
    const int ny = grid.ny();
    const int x_center = nx / 2;

    for (int y = 0; y < ny; ++y) {
        const auto& fluid = grid(x_center, y).fluid;
        file << y << " " << fluid.u[0] << " " << fluid.u[1] << "\n";
    }

    file.close();
}

// Run simulation
void runSimulation(const SimulationParams& params) {
    using Lattice = D2Q9;
    const int D = Lattice::D;
    const int Q = Lattice::Q;

    std::cout << "=== ELBM Rectangular Pipe Flow Simulation ===" << std::endl;
    std::cout << "Grid: " << params.nx << " x " << params.ny << std::endl;
    std::cout << "Viscosity: " << params.viscosity << " m^2/s" << std::endl;
    std::cout << "Pressure drop: " << params.delta_p << " Pa" << std::endl;
    std::cout << "Solver: " << params.solver_type << std::endl;
    std::cout << "Max steps: " << params.max_steps << std::endl;
    std::cout << std::endl;

    // Create grid
    LatticeGrid<D, Q> grid(params.nx, params.ny);

    // Create entropy monitor for H-theorem validation
    EntropyMonitor monitor;

    // Compute densities from pressure drop
    // p = cs² * rho, so rho = p / cs²
    const double rho_left = (1.0 + params.delta_p) / Lattice::cs2;
    const double rho_right = 1.0 / Lattice::cs2;

    std::cout << "Left density: " << rho_left << std::endl;
    std::cout << "Right density: " << rho_right << std::endl;
    std::cout << std::endl;

    // Initialize with pressure drop
    initializePressureDrop<Lattice>(grid, rho_left, rho_right);

    // Create solver
    if (params.solver_type == "BGK") {
        BGKSolver<Lattice> solver(params.viscosity);
        std::cout << "BGK Solver initialized" << std::endl;
        std::cout << "Relaxation time (tau): " << solver.tau() << std::endl;
        std::cout << "Omega: " << solver.omega() << std::endl;
        std::cout << std::endl;

        // Setup boundaries
        BoundaryManager<Lattice> bc_manager;
        bc_manager.setLeftBC(BCType::PRESSURE, rho_left);
        bc_manager.setRightBC(BCType::PRESSURE, rho_right);
        bc_manager.setTopBC(BCType::BOUNCE_BACK);
        bc_manager.setBottomBC(BCType::BOUNCE_BACK);

        // Time loop
        auto start_time = std::chrono::high_resolution_clock::now();

        for (int step = 0; step <= params.max_steps; ++step) {
            // Collision
            for (int y = 0; y < params.ny; ++y) {
                for (int x = 0; x < params.nx; ++x) {
                    solver.collide(grid(x, y));
                }
            }

            // Streaming
            solver.stream(grid);

            // Boundary conditions (MUST be after streaming!)
            bc_manager.applyBoundaries(grid);

            // Record entropy for H-theorem validation
            if (step % params.output_interval == 0) {
                monitor.record(step, grid, solver);
            }

            // Output
            if (step % params.output_interval == 0) {
                std::cout << "Step " << step << "/" << params.max_steps;
                if (!monitor.getHistory().empty()) {
                    const auto& rec = monitor.getHistory().back();
                    std::cout << " | H=" << std::scientific << std::setprecision(4) << rec.H_total;
                }
                std::cout << std::endl;

                std::string filename = params.output_dir + "/bgk_step_" +
                                       std::to_string(step) + ".dat";
                writeResults<Lattice>(grid, filename, step);

                std::string centerline_file = params.output_dir + "/bgk_centerline_" +
                                              std::to_string(step) + ".dat";
                writeCenterlineProfile<Lattice>(grid, centerline_file);

                std::string velocity_file = params.output_dir + "/bgk_velocity_" +
                                           std::to_string(step) + ".dat";
                writeVelocityProfile<Lattice>(grid, velocity_file);
            }
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Simulation completed in " << duration.count() << " ms" << std::endl;

        // Save entropy history
        std::string entropy_file = params.output_dir + "/entropy_bgk.dat";
        monitor.saveToFile(entropy_file);
        std::cout << "Entropy history saved to: " << entropy_file << std::endl;

        // Print H-theorem validation
        monitor.printSummary();

    } else if (params.solver_type == "ELBM") {
        ELBMSolver<Lattice> solver(params.viscosity);
        std::cout << "ELBM Solver initialized" << std::endl;
        std::cout << "Beta: " << solver.beta() << std::endl;
        std::cout << std::endl;

        // Setup boundaries
        BoundaryManager<Lattice> bc_manager;
        bc_manager.setLeftBC(BCType::PRESSURE, rho_left);
        bc_manager.setRightBC(BCType::PRESSURE, rho_right);
        bc_manager.setTopBC(BCType::BOUNCE_BACK);
        bc_manager.setBottomBC(BCType::BOUNCE_BACK);

        // Time loop
        auto start_time = std::chrono::high_resolution_clock::now();

        for (int step = 0; step <= params.max_steps; ++step) {
            // Collision
            for (int y = 0; y < params.ny; ++y) {
                for (int x = 0; x < params.nx; ++x) {
                    solver.collide(grid(x, y));
                }
            }

            // Streaming
            solver.stream(grid);

            // Boundary conditions (MUST be after streaming!)
            bc_manager.applyBoundaries(grid);

            // Record entropy for H-theorem validation
            if (step % params.output_interval == 0) {
                monitor.record(step, grid, solver);
            }

            // Output
            if (step % params.output_interval == 0) {
                std::cout << "Step " << step << "/" << params.max_steps;
                if (!monitor.getHistory().empty()) {
                    const auto& rec = monitor.getHistory().back();
                    std::cout << " | H=" << std::scientific << std::setprecision(4) << rec.H_total;
                }
                std::cout << std::endl;

                std::string filename = params.output_dir + "/elbm_step_" +
                                       std::to_string(step) + ".dat";
                writeResults<Lattice>(grid, filename, step);

                std::string centerline_file = params.output_dir + "/elbm_centerline_" +
                                              std::to_string(step) + ".dat";
                writeCenterlineProfile<Lattice>(grid, centerline_file);

                std::string velocity_file = params.output_dir + "/elbm_velocity_" +
                                           std::to_string(step) + ".dat";
                writeVelocityProfile<Lattice>(grid, velocity_file);
            }
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "Simulation completed in " << duration.count() << " ms" << std::endl;

        // Save entropy history
        std::string entropy_file = params.output_dir + "/entropy_elbm.dat";
        monitor.saveToFile(entropy_file);
        std::cout << "Entropy history saved to: " << entropy_file << std::endl;

        // Print H-theorem validation
        monitor.printSummary();
    }
}

int main(int argc, char* argv[]) {
    // Default parameters for Re ~ 10 case
    SimulationParams params;
    params.nx = 100;
    params.ny = 40;
    params.viscosity = 0.1;      // m^2/s
    params.delta_p = 5.0;        // Pa
    params.max_steps = 2500;
    params.output_interval = 250;
    params.solver_type = "BGK";
    params.output_dir = "output";

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--solver" && i + 1 < argc) {
            params.solver_type = argv[++i];
        } else if (arg == "--viscosity" && i + 1 < argc) {
            params.viscosity = std::stod(argv[++i]);
        } else if (arg == "--nx" && i + 1 < argc) {
            params.nx = std::stoi(argv[++i]);
        } else if (arg == "--ny" && i + 1 < argc) {
            params.ny = std::stoi(argv[++i]);
        } else if (arg == "--steps" && i + 1 < argc) {
            params.max_steps = std::stoi(argv[++i]);
        } else if (arg == "--output-dir" && i + 1 < argc) {
            params.output_dir = argv[++i];
        }
    }

    // Create output directory
    system(("mkdir -p " + params.output_dir).c_str());

    // Run simulation
    runSimulation(params);

    return 0;
}
