/**
 * Two-Phase Immiscible Flow: H-Theorem Validation
 *
 * Demonstrates H-theorem validity (dH/dt ≤ 0) for two-phase flows
 * Compares BGK and ELBM collision operators on a two-phase density configuration
 */

#include "core/lattice.h"
#include "core/fluid_state.h"
#include "core/entropy_monitor.h"
#include "core/interface_detector.h"
#include "solvers/bgk_solver.h"
#include "solvers/elbm_solver.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>

using namespace elbm;

// Lattice parameters
constexpr int NX = 128;
constexpr int NY = 128;
constexpr int TIMESTEPS = 5000;
constexpr int OUTPUT_INTERVAL = 250;
constexpr int DIAGNOSTIC_INTERVAL = 50; // Record diagnostics more frequently

// Physical parameters
constexpr double RHO_LIQUID = 1.2;
constexpr double RHO_GAS = 0.8;
constexpr double DROPLET_RADIUS = 30.0;
constexpr double VISCOSITY = 0.1;

using Lattice2D = D2Q9;

/**
 * Initialize droplet with smooth interface (no forces, just initial condition)
 */
template<int D, int Q>
void initializeDroplet(LatticeGrid<D, Q>& grid, double R) {
    int cx = NX / 2;
    int cy = NY / 2;

    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            auto& node = grid(x, y);

            // Distance from center
            double dx = x - cx;
            double dy = y - cy;
            double r = std::sqrt(dx*dx + dy*dy);

            // Smooth interface using tanh profile
            double interface_width = 5.0;
            double rho = RHO_GAS + 0.5 * (RHO_LIQUID - RHO_GAS) *
                         (1.0 - std::tanh(2.0 * (r - R) / interface_width));

            node.fluid.rho = rho;
            node.fluid.u.fill(0.0);

            // Initialize distributions at equilibrium
            // Use entropic equilibrium for thermodynamic consistency (esp. for ELBM)
            EntropicEquilibrium<Lattice2D>::compute(node.fluid, node.df.f);
        }
    }
}

/**
 * Save grid state
 */
template<int D, int Q>
void saveState(const LatticeGrid<D, Q>& grid, const std::string& filename) {
    std::ofstream file(filename);
    file << "# x y rho ux uy p\n";
    file << std::scientific << std::setprecision(8);

    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            const auto& node = grid(x, y);
            file << x << " " << y << " "
                 << node.fluid.rho << " "
                 << node.fluid.u[0] << " "
                 << node.fluid.u[1] << " "
                 << node.fluid.p << "\n";
        }
    }
    file.close();
}

/**
 * Compute macroscopic variables
 */
template<int D, int Q>
void computeMacro(LatticeGrid<D, Q>& grid) {
    for (auto& node : grid) {
        if (node.is_solid) continue;

        double rho = 0.0;
        std::array<double, D> u;
        u.fill(0.0);

        for (int i = 0; i < Q; ++i) {
            rho += node.df.f[i];
            u[0] += node.df.f[i] * Lattice2D::cx(i);
            u[1] += node.df.f[i] * Lattice2D::cy(i);
        }

        if (rho > 1e-14) {
            for (int d = 0; d < D; ++d) {
                u[d] /= rho;
            }
        }

        node.fluid.rho = rho;
        node.fluid.u = u;
        node.fluid.p = rho / 3.0;
    }
}

/**
 * Streaming step
 */
template<int D, int Q>
void stream(LatticeGrid<D, Q>& grid) {
    LatticeGrid<D, Q> temp(NX, NY);

    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            for (int i = 0; i < Q; ++i) {
                int xn = (x + Lattice2D::cx(i) + NX) % NX;
                int yn = (y + Lattice2D::cy(i) + NY) % NY;
                temp(xn, yn).df.f[i] = grid(x, y).df.f_new[i];
            }
        }
    }

    for (int y = 0; y < NY; ++y) {
        for (int x = 0; x < NX; ++x) {
            grid(x, y).df.f = temp(x, y).df.f;
        }
    }
}

/**
 * Run BGK simulation
 */
void runBGK() {
    std::cout << "\n=== Running BGK Simulation ===\n";

    LatticeGrid<2, 9> grid(NX, NY);
    BGKSolver<Lattice2D> solver(VISCOSITY);
    EntropyMonitor monitor;
    InterfaceDetector detector(RHO_LIQUID, RHO_GAS, 0.1);

    initializeDroplet(grid, DROPLET_RADIUS);
    std::cout << "Two-phase configuration initialized\n";
    std::cout << "Recording diagnostics every " << DIAGNOSTIC_INTERVAL << " steps\n";

    int center_x = NX / 2;
    int center_y = NY / 2;

    for (int t = 0; t <= TIMESTEPS; ++t) {
        // Simple collision (no forces)
        for (int y = 0; y < NY; ++y) {
            for (int x = 0; x < NX; ++x) {
                solver.collide(grid(x, y));
            }
        }

        stream(grid);
        computeMacro(grid);

        // Record diagnostics more frequently for time-series analysis
        if (t % DIAGNOSTIC_INTERVAL == 0) {
            monitor.recordWithDiagnostics(t, grid, solver, detector, center_x, center_y);
        }

        // Save spatial snapshots less frequently
        if (t % OUTPUT_INTERVAL == 0) {
            saveState(grid, "output/bgk_t" + std::to_string(t) + ".dat");

            if (!monitor.getHistory().empty()) {
                const auto& rec = monitor.getHistory().back();
                const auto& m = rec.macro;

                std::cout << "\n  t=" << std::setw(5) << t << " | Convergence Progress:\n";
                std::cout << "    Entropy:   H=" << std::fixed << std::setprecision(4) << rec.H_total
                          << " (range: " << rec.H_min << " - " << rec.H_max << ")\n";
                std::cout << "    Pressure:  ΔP=" << std::scientific << std::setprecision(3) << m.laplace_pressure
                          << " | σ=" << m.surface_tension
                          << " | P_in=" << m.pressure_liquid
                          << " | P_out=" << m.pressure_gas << "\n";
                std::cout << "    Interface: R=" << std::fixed << std::setprecision(2) << m.interface_radius
                          << " | width=" << m.interface_width
                          << " | nodes=" << m.num_interface_nodes << "\n";
                std::cout << "    Spurious:  global=" << std::scientific << std::setprecision(3) << rec.spurious_vel
                          << " | interface=" << m.spurious_vel_interface << "\n";
                std::cout << "    Conservation: ΔM=" << std::abs(m.total_mass - (NX*NY*RHO_LIQUID + NX*NY*RHO_GAS)/2.0)
                          << " | p_mag=" << m.total_momentum_mag << "\n";
            }
        }
    }

    // Save both formats
    monitor.saveToFile("output/entropy_bgk.dat");
    monitor.saveDiagnosticsToCSV("output/bgk_diagnostics.csv");
    monitor.printSummary();

    std::cout << "\nDiagnostics saved to: output/bgk_diagnostics.csv\n";
}

/**
 * Run ELBM simulation
 */
void runELBM() {
    std::cout << "\n=== Running ELBM Simulation ===\n";

    LatticeGrid<2, 9> grid(NX, NY);
    ELBMSolver<Lattice2D> solver(VISCOSITY);
    EntropyMonitor monitor;
    InterfaceDetector detector(RHO_LIQUID, RHO_GAS, 0.1);

    initializeDroplet(grid, DROPLET_RADIUS);
    std::cout << "Two-phase configuration initialized\n";
    std::cout << "Recording diagnostics every " << DIAGNOSTIC_INTERVAL << " steps\n";

    int center_x = NX / 2;
    int center_y = NY / 2;

    // Track alpha statistics
    double alpha_sum = 0.0;
    double alpha_min = 2.0;
    double alpha_max = 0.0;
    int alpha_count = 0;

    for (int t = 0; t <= TIMESTEPS; ++t) {
        // ELBM collision (no forces)
        for (int y = 0; y < NY; ++y) {
            for (int x = 0; x < NX; ++x) {
                solver.collide(grid(x, y));

                // Track alpha statistics
                double alpha = grid(x, y).alpha;
                alpha_sum += alpha;
                alpha_min = std::min(alpha_min, alpha);
                alpha_max = std::max(alpha_max, alpha);
                alpha_count++;
            }
        }

        stream(grid);
        computeMacro(grid);

        // Record diagnostics more frequently for time-series analysis
        if (t % DIAGNOSTIC_INTERVAL == 0) {
            monitor.recordWithDiagnostics(t, grid, solver, detector, center_x, center_y);
        }

        // Save spatial snapshots less frequently
        if (t % OUTPUT_INTERVAL == 0) {
            saveState(grid, "output/elbm_t" + std::to_string(t) + ".dat");

            if (!monitor.getHistory().empty()) {
                const auto& rec = monitor.getHistory().back();
                const auto& m = rec.macro;

                double alpha_mean = alpha_sum / alpha_count;

                std::cout << "\n  t=" << std::setw(5) << t << " | Convergence Progress:\n";
                std::cout << "    Entropy:   H=" << std::fixed << std::setprecision(4) << rec.H_total
                          << " (range: " << rec.H_min << " - " << rec.H_max << ")\n";
                std::cout << "    Pressure:  ΔP=" << std::scientific << std::setprecision(3) << m.laplace_pressure
                          << " | σ=" << m.surface_tension
                          << " | P_in=" << m.pressure_liquid
                          << " | P_out=" << m.pressure_gas << "\n";
                std::cout << "    Interface: R=" << std::fixed << std::setprecision(2) << m.interface_radius
                          << " | width=" << m.interface_width
                          << " | nodes=" << m.num_interface_nodes << "\n";
                std::cout << "    Spurious:  global=" << std::scientific << std::setprecision(3) << rec.spurious_vel
                          << " | interface=" << m.spurious_vel_interface << "\n";
                std::cout << "    Alpha:     mean=" << std::fixed << std::setprecision(4) << alpha_mean
                          << " | range=[" << alpha_min << ", " << alpha_max << "]"
                          << " | center=" << grid(NX/2, NY/2).alpha << "\n";
                std::cout << "    Conservation: ΔM=" << std::scientific << std::setprecision(3)
                          << std::abs(m.total_mass - (NX*NY*RHO_LIQUID + NX*NY*RHO_GAS)/2.0)
                          << " | p_mag=" << m.total_momentum_mag << "\n";

                // Reset alpha statistics for next interval
                alpha_sum = 0.0;
                alpha_min = 2.0;
                alpha_max = 0.0;
                alpha_count = 0;
            }
        }
    }

    // Save both formats
    monitor.saveToFile("output/entropy_elbm.dat");
    monitor.saveDiagnosticsToCSV("output/elbm_diagnostics.csv");
    monitor.printSummary();

    std::cout << "\nDiagnostics saved to: output/elbm_diagnostics.csv\n";
}

int main() {
    std::cout << "========================================================\n";
    std::cout << "Two-Phase Immiscible Flow: H-Theorem Validation\n";
    std::cout << "========================================================\n";
    std::cout << "Problem: Two-phase density field\n";
    std::cout << "Domain: " << NX << "×" << NY << "\n";
    std::cout << "Phases:\n";
    std::cout << "  - Higher density: ρ = " << RHO_LIQUID << "\n";
    std::cout << "  - Lower density:  ρ = " << RHO_GAS << "\n";
    std::cout << "Density ratio: " << RHO_LIQUID/RHO_GAS << ":1\n";
    std::cout << "Timesteps: " << TIMESTEPS << "\n";
    std::cout << "\nObjective: Validate dH/dt ≤ 0 for both BGK and ELBM\n";
    std::cout << "========================================================\n";

    runBGK();
    runELBM();

    std::cout << "\n=== Results Summary ===\n";
    std::cout << "Both methods should satisfy H-theorem (dH/dt ≤ 0)\n";
    std::cout << "ELBM should show:\n";
    std::cout << "  - Similar or lower spurious currents\n";
    std::cout << "  - Stable entropy decrease\n";
    std::cout << "  - α parameter close to 2.0 (minimal correction needed)\n";
    std::cout << "\nVisualize: python plotting/plot_twophase_results.py\n";

    return 0;
}
