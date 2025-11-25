#ifndef ENTROPY_MONITOR_H
#define ENTROPY_MONITOR_H

#include "core/fluid_state.h"
#include "core/interface_detector.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

namespace elbm {

/**
 * Entropy Monitor for H-Theorem Validation
 *
 * Tracks the discrete entropy H(t) = Σ f_i · ln(f_i / w_i)
 * Validates H-theorem: dH/dt ≤ 0
 * Also tracks comprehensive macroscopic diagnostics for two-phase flows
 */
class EntropyMonitor {
public:
    /**
     * Macroscopic diagnostics for two-phase flows
     */
    struct MacroscopicDiagnostics {
        // Pressures
        double pressure_liquid;
        double pressure_gas;
        double pressure_interface;
        double laplace_pressure;    // ΔP = P_liquid - P_gas
        double surface_tension;     // σ = ΔP * R

        // Interface properties
        double interface_radius;
        double interface_width;
        double interface_position_x;
        double interface_position_y;
        int num_interface_nodes;
        double spurious_vel_interface;

        // Conservation laws
        double total_mass;
        double total_momentum_x;
        double total_momentum_y;
        double total_momentum_mag;
        double kinetic_energy;

        // Initialize with zeros
        MacroscopicDiagnostics()
            : pressure_liquid(0), pressure_gas(0), pressure_interface(0)
            , laplace_pressure(0), surface_tension(0)
            , interface_radius(0), interface_width(0)
            , interface_position_x(0), interface_position_y(0)
            , num_interface_nodes(0), spurious_vel_interface(0)
            , total_mass(0), total_momentum_x(0), total_momentum_y(0)
            , total_momentum_mag(0), kinetic_energy(0) {}
    };

    struct Record {
        int timestep;
        double H_total;      // Total entropy (sum over all nodes)
        double H_min;        // Minimum H across grid
        double H_max;        // Maximum H across grid
        double H_mean;       // Mean H
        double spurious_vel; // Maximum spurious velocity

        // Per-phase entropy tracking (NEW: addresses H-convergence issue)
        double H_liquid;     // Total H in liquid phase
        double H_gas;        // Total H in gas phase
        double H_interface;  // Total H at interface
        double H_liquid_mean;   // Mean H per node in liquid
        double H_gas_mean;      // Mean H per node in gas
        double H_interface_mean; // Mean H per node at interface
        int nodes_liquid;    // Number of liquid nodes
        int nodes_gas;       // Number of gas nodes
        int nodes_interface; // Number of interface nodes

        MacroscopicDiagnostics macro; // Macroscopic diagnostics
    };

    EntropyMonitor() = default;

    /**
     * Record entropy for current timestep
     */
    template<int D, int Q, typename Solver>
    void record(int timestep,
                const LatticeGrid<D, Q>& grid,
                const Solver& solver) {

        Record rec;
        rec.timestep = timestep;
        rec.H_total = 0.0;
        rec.H_min = std::numeric_limits<double>::max();
        rec.H_max = std::numeric_limits<double>::lowest();
        rec.spurious_vel = 0.0;

        int count = 0;
        for (int y = 0; y < grid.ny(); ++y) {
            for (int x = 0; x < grid.nx(); ++x) {
                const auto& node = grid(x, y);

                if (node.is_solid) continue;

                // Compute H for this node
                double H = solver.computeH(node.df.f);

                rec.H_total += H;
                rec.H_min = std::min(rec.H_min, H);
                rec.H_max = std::max(rec.H_max, H);

                // Track maximum velocity (spurious currents)
                double vel_mag = 0.0;
                for (int d = 0; d < D; ++d) {
                    vel_mag += node.fluid.u[d] * node.fluid.u[d];
                }
                vel_mag = std::sqrt(vel_mag);
                rec.spurious_vel = std::max(rec.spurious_vel, vel_mag);

                count++;
            }
        }

        rec.H_mean = rec.H_total / count;
        history_.push_back(rec);
    }

    /**
     * Record entropy + macroscopic diagnostics for two-phase flows
     * @param timestep - Current timestep
     * @param grid - Lattice grid
     * @param solver - LBM solver (BGK or ELBM)
     * @param detector - Interface detector
     * @param center_x, center_y - Droplet center for interface analysis
     */
    template<int D, int Q, typename Solver>
    void recordWithDiagnostics(int timestep,
                              const LatticeGrid<D, Q>& grid,
                              const Solver& solver,
                              const InterfaceDetector& detector,
                              int center_x, int center_y) {

        Record rec;
        rec.timestep = timestep;
        rec.H_total = 0.0;
        rec.H_min = std::numeric_limits<double>::max();
        rec.H_max = std::numeric_limits<double>::lowest();
        rec.spurious_vel = 0.0;

        // Per-phase entropy initialization
        rec.H_liquid = 0.0;
        rec.H_gas = 0.0;
        rec.H_interface = 0.0;
        rec.nodes_liquid = 0;
        rec.nodes_gas = 0;
        rec.nodes_interface = 0;

        // Basic entropy and velocity tracking with phase separation
        int count = 0;
        for (const auto& node : grid) {
            if (node.is_solid) continue;

            // Compute H for this node
            double H = solver.computeH(node.df.f);

            rec.H_total += H;
            rec.H_min = std::min(rec.H_min, H);
            rec.H_max = std::max(rec.H_max, H);

            // Classify node and accumulate per-phase entropy
            InterfaceDetector::Phase phase = detector.classify(node.fluid.rho);
            switch (phase) {
                case InterfaceDetector::Phase::LIQUID:
                    rec.H_liquid += H;
                    rec.nodes_liquid++;
                    break;
                case InterfaceDetector::Phase::GAS:
                    rec.H_gas += H;
                    rec.nodes_gas++;
                    break;
                case InterfaceDetector::Phase::INTERFACE:
                    rec.H_interface += H;
                    rec.nodes_interface++;
                    break;
                default:
                    break;
            }

            // Track maximum velocity (spurious currents)
            double vel_mag = 0.0;
            for (int d = 0; d < D; ++d) {
                vel_mag += node.fluid.u[d] * node.fluid.u[d];
            }
            vel_mag = std::sqrt(vel_mag);
            rec.spurious_vel = std::max(rec.spurious_vel, vel_mag);

            count++;
        }

        rec.H_mean = rec.H_total / count;

        // Compute mean H per node for each phase
        rec.H_liquid_mean = (rec.nodes_liquid > 0) ? rec.H_liquid / rec.nodes_liquid : 0.0;
        rec.H_gas_mean = (rec.nodes_gas > 0) ? rec.H_gas / rec.nodes_gas : 0.0;
        rec.H_interface_mean = (rec.nodes_interface > 0) ? rec.H_interface / rec.nodes_interface : 0.0;

        // Compute macroscopic diagnostics
        computeMacroscopicDiagnostics(grid, detector, center_x, center_y, rec.macro);

        history_.push_back(rec);
    }

    /**
     * Compute macroscopic diagnostics for two-phase flows
     */
    template<int D, int Q>
    void computeMacroscopicDiagnostics(const LatticeGrid<D, Q>& grid,
                                      const InterfaceDetector& detector,
                                      int center_x, int center_y,
                                      MacroscopicDiagnostics& diag) {

        // Compute phase-averaged pressures
        auto pressures = detector.computePhaseAveragedPressure(grid);
        diag.pressure_liquid = pressures[0];
        diag.pressure_gas = pressures[1];
        diag.pressure_interface = pressures[2];

        // Analyze interface properties
        auto interface_stats = detector.analyzeInterface(grid, center_x, center_y);
        diag.interface_radius = interface_stats.radius;
        diag.interface_width = interface_stats.width;
        diag.interface_position_x = interface_stats.mean_position_x;
        diag.interface_position_y = interface_stats.mean_position_y;
        diag.num_interface_nodes = interface_stats.num_interface_nodes;
        diag.spurious_vel_interface = interface_stats.spurious_vel_interface;

        // Compute Laplace pressure and surface tension
        auto laplace = detector.computeLaplaceParameters(
            diag.pressure_liquid, diag.pressure_gas, diag.interface_radius);
        diag.laplace_pressure = laplace[0];
        diag.surface_tension = laplace[1];

        // Compute conservation laws
        diag.total_mass = 0.0;
        diag.total_momentum_x = 0.0;
        diag.total_momentum_y = 0.0;
        diag.kinetic_energy = 0.0;

        for (const auto& node : grid) {
            if (node.is_solid) continue;

            diag.total_mass += node.fluid.rho;
            diag.total_momentum_x += node.fluid.rho * node.fluid.u[0];
            diag.total_momentum_y += node.fluid.rho * node.fluid.u[1];

            double vel_sq = 0.0;
            for (int d = 0; d < D; ++d) {
                vel_sq += node.fluid.u[d] * node.fluid.u[d];
            }
            diag.kinetic_energy += 0.5 * node.fluid.rho * vel_sq;
        }

        diag.total_momentum_mag = std::sqrt(
            diag.total_momentum_x * diag.total_momentum_x +
            diag.total_momentum_y * diag.total_momentum_y);
    }

    /**
     * Check if H-theorem is satisfied
     * Returns true if dH/dt ≤ 0 for all recorded timesteps
     */
    bool verifyHTheorem() const {
        if (history_.size() < 2) return true;

        for (size_t i = 1; i < history_.size(); ++i) {
            double dH = history_[i].H_total - history_[i-1].H_total;
            if (dH > 1e-10) { // Small tolerance for numerical errors
                return false;
            }
        }
        return true;
    }

    /**
     * Get maximum violation of H-theorem
     * Returns max(H(t+1) - H(t)) over all timesteps
     */
    double getMaxViolation() const {
        if (history_.size() < 2) return 0.0;

        double max_violation = 0.0;
        for (size_t i = 1; i < history_.size(); ++i) {
            double dH = history_[i].H_total - history_[i-1].H_total;
            max_violation = std::max(max_violation, dH);
        }
        return max_violation;
    }

    /**
     * Save history to file (entropy only - legacy format)
     */
    void saveToFile(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        // Header
        file << "# H-function history for H-theorem validation\n";
        file << "# timestep H_total H_min H_max H_mean spurious_vel\n";
        file << std::scientific << std::setprecision(12);

        // Data
        for (const auto& rec : history_) {
            file << rec.timestep << " "
                 << rec.H_total << " "
                 << rec.H_min << " "
                 << rec.H_max << " "
                 << rec.H_mean << " "
                 << rec.spurious_vel << "\n";
        }

        file.close();
    }

    /**
     * Save comprehensive diagnostics to CSV file for time-series analysis
     */
    void saveDiagnosticsToCSV(const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }

        // CSV Header
        file << "timestep,H_total,H_min,H_max,H_mean,spurious_vel,"
             << "H_liquid,H_gas,H_interface,H_liquid_mean,H_gas_mean,H_interface_mean,nodes_liquid,nodes_gas,nodes_interface,"
             << "P_liquid,P_gas,P_interface,laplace_pressure,surface_tension,"
             << "interface_radius,interface_width,interface_pos_x,interface_pos_y,num_interface_nodes,spurious_vel_interface,"
             << "total_mass,momentum_x,momentum_y,momentum_mag,kinetic_energy\n";

        file << std::scientific << std::setprecision(12);

        // Data rows
        for (const auto& rec : history_) {
            const auto& m = rec.macro; // alias for macroscopic diagnostics

            file << rec.timestep << ","
                 << rec.H_total << ","
                 << rec.H_min << ","
                 << rec.H_max << ","
                 << rec.H_mean << ","
                 << rec.spurious_vel << ","
                 << rec.H_liquid << ","
                 << rec.H_gas << ","
                 << rec.H_interface << ","
                 << rec.H_liquid_mean << ","
                 << rec.H_gas_mean << ","
                 << rec.H_interface_mean << ","
                 << rec.nodes_liquid << ","
                 << rec.nodes_gas << ","
                 << rec.nodes_interface << ","
                 << m.pressure_liquid << ","
                 << m.pressure_gas << ","
                 << m.pressure_interface << ","
                 << m.laplace_pressure << ","
                 << m.surface_tension << ","
                 << m.interface_radius << ","
                 << m.interface_width << ","
                 << m.interface_position_x << ","
                 << m.interface_position_y << ","
                 << m.num_interface_nodes << ","
                 << m.spurious_vel_interface << ","
                 << m.total_mass << ","
                 << m.total_momentum_x << ","
                 << m.total_momentum_y << ","
                 << m.total_momentum_mag << ","
                 << m.kinetic_energy << "\n";
        }

        file.close();
    }

    /**
     * Print summary statistics
     */
    void printSummary() const {
        if (history_.empty()) {
            std::cout << "No entropy data recorded.\n";
            return;
        }

        const auto& first = history_.front();
        const auto& last = history_.back();

        std::cout << "\n=== Entropy Monitor Summary ===\n";
        std::cout << "Timesteps: " << first.timestep << " -> " << last.timestep << "\n";
        std::cout << "\nGlobal Entropy:\n";
        std::cout << "  Initial H: " << std::scientific << std::setprecision(6) << first.H_total << "\n";
        std::cout << "  Final H:   " << last.H_total << "\n";
        std::cout << "  ΔH:        " << (last.H_total - first.H_total) << "\n";
        std::cout << "  H-theorem: " << (verifyHTheorem() ? "SATISFIED ✓" : "VIOLATED ✗") << "\n";
        std::cout << "  Max violation: " << getMaxViolation() << "\n";

        // Per-phase entropy (if tracked)
        if (last.nodes_liquid > 0 || last.nodes_gas > 0) {
            std::cout << "\nPer-Phase Entropy (Mean H per node):\n";
            std::cout << "  Liquid phase:    H_mean = " << std::fixed << std::setprecision(4)
                      << last.H_liquid_mean << " (" << last.nodes_liquid << " nodes)\n";
            std::cout << "  Gas phase:       H_mean = " << last.H_gas_mean
                      << " (" << last.nodes_gas << " nodes)\n";
            std::cout << "  Interface:       H_mean = " << last.H_interface_mean
                      << " (" << last.nodes_interface << " nodes)\n";

            std::cout << "\n  IMPORTANT: Different phases have DIFFERENT entropy densities!\n";
            std::cout << "  This is thermodynamically correct - H_liquid ≠ H_gas at equilibrium.\n";
            std::cout << "  The H-theorem (dH/dt ≤ 0) applies to TOTAL entropy, not per-phase values.\n";
        }

        std::cout << "\nFinal spurious velocity: " << std::scientific << last.spurious_vel << "\n";
        std::cout << "================================\n\n";
    }

    /**
     * Get history for external analysis
     */
    const std::vector<Record>& getHistory() const {
        return history_;
    }

    /**
     * Clear history
     */
    void clear() {
        history_.clear();
    }

private:
    std::vector<Record> history_;
};

} // namespace elbm

#endif // ENTROPY_MONITOR_H
