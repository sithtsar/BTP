#include "../../include/analysis/h_theorem.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

namespace lbm {
namespace analysis {

HTheoremAnalyzer::HTheoremAnalyzer(int Nx, int Ny, const AnalysisConfig& config)
    : Nx_(Nx), Ny_(Ny), config_(config), previous_h_value_(0.0), 
      has_previous_data_(false), last_analysis_timestep_(-1) {
    
    // Allocate working arrays
    h_field_flat_ = std::make_unique<double[]>(Nx * Ny);
    feq_temp_ = std::make_unique<double[]>(Q);
    
    // Initialize evolution data
    evolution_data_.timesteps.reserve(1000);
    evolution_data_.h_values.reserve(1000);
    evolution_data_.h_changes.reserve(1000);
    evolution_data_.entropy_production_rates.reserve(1000);
    evolution_data_.monotonicity_violations.reserve(1000);
    
    if (config_.enable_statistical_analysis) {
        evolution_data_.h_variances.reserve(1000);
        evolution_data_.kinetic_energies.reserve(1000);
        evolution_data_.total_masses.reserve(1000);
    }
}

HTheoremAnalyzer::HTheoremAnalyzer(int Nx, int Ny)
    : HTheoremAnalyzer(Nx, Ny, AnalysisConfig{}) {
}

HTheoremAnalyzer::HTheoremMetrics HTheoremAnalyzer::calculateHFunctionStandardBGK(
    double*** f, double** rho, double** ux, double** uy, int timestep) {
    
    HTheoremMetrics metrics;
    metrics.timestep = timestep;
    
    double total_h = 0.0;
    double total_mass = 0.0;
    double total_kinetic_energy = 0.0;
    double h_min = std::numeric_limits<double>::max();
    double h_max = std::numeric_limits<double>::lowest();
    
    std::vector<double> h_values;
    if (config_.enable_statistical_analysis) {
        h_values.reserve(Nx_ * Ny_);
    }
    
    // Calculate H-function at each grid point
    for (int x = 0; x < Nx_; ++x) {
        for (int y = 0; y < Ny_; ++y) {
            // Calculate equilibrium distributions
            for (int q = 0; q < Q; ++q) {
                feq_temp_[q] = calculateEquilibrium(q, rho[x][y], ux[x][y], uy[x][y]);
            }
            
            // Calculate local H-function using standard BGK form
            double local_h = 0.0;
            for (int q = 0; q < Q; ++q) {
                if (f[q][x][y] > 1e-15 && feq_temp_[q] > 1e-15) {
                    local_h += f[q][x][y] * std::log(f[q][x][y] / feq_temp_[q]);
                }
            }
            
            // Check for numerical anomalies
            if (std::isnan(local_h) || std::isinf(local_h)) {
                metrics.has_nan_values = true;
                local_h = 0.0; // Set to zero to avoid contaminating the sum
            }
            
            total_h += local_h;
            total_mass += rho[x][y];
            total_kinetic_energy += 0.5 * rho[x][y] * (ux[x][y] * ux[x][y] + uy[x][y] * uy[x][y]);
            
            h_min = std::min(h_min, local_h);
            h_max = std::max(h_max, local_h);
            
            if (config_.enable_statistical_analysis) {
                h_values.push_back(local_h);
            }
            
            // Store in flat array for spatial analysis
            h_field_flat_[x * Ny_ + y] = local_h;
        }
    }
    
    // Set basic metrics
    metrics.h_function_value = total_h;
    metrics.total_mass = total_mass;
    metrics.kinetic_energy = total_kinetic_energy;
    metrics.h_function_min = h_min;
    metrics.h_function_max = h_max;
    
    // Calculate H-function change and monotonicity
    if (has_previous_data_) {
        metrics.h_function_change = total_h - previous_h_value_;
        metrics.monotonic_decrease = (metrics.h_function_change <= config_.monotonicity_tolerance);
        
        if (!metrics.monotonic_decrease && config_.log_violations) {
            logMonotonicityViolation(timestep, total_h, previous_h_value_);
        }
        
        // Calculate entropy production rate (assuming unit time step)
        metrics.entropy_production = calculateEntropyProductionRate(total_h, previous_h_value_, 1.0);
    } else {
        metrics.h_function_change = 0.0;
        metrics.monotonic_decrease = true;
        metrics.entropy_production = 0.0;
    }
    
    // Calculate spatial statistics
    if (config_.enable_statistical_analysis && !h_values.empty()) {
        double mean = std::accumulate(h_values.begin(), h_values.end(), 0.0) / h_values.size();
        double variance = 0.0;
        for (double h : h_values) {
            variance += (h - mean) * (h - mean);
        }
        variance /= h_values.size();
        metrics.h_function_variance = variance;
        
        // Calculate maximum local change
        if (has_previous_data_) {
            double max_change = 0.0;
            for (int i = 0; i < Nx_ * Ny_; ++i) {
                // This is a simplified version - in practice, we'd need to store previous field
                max_change = std::max(max_change, std::abs(h_field_flat_[i]));
            }
            metrics.max_local_change = max_change;
        }
    }
    
    // Update previous data
    previous_h_value_ = total_h;
    has_previous_data_ = true;
    last_analysis_timestep_ = timestep;
    
    return metrics;
}

HTheoremAnalyzer::HTheoremMetrics HTheoremAnalyzer::calculateHFunctionEntropicBGK(
    double*** f, double** rho, double** ux, double** uy, int timestep) {
    
    HTheoremMetrics metrics;
    metrics.timestep = timestep;
    
    double total_h = 0.0;
    double total_mass = 0.0;
    double total_kinetic_energy = 0.0;
    double h_min = std::numeric_limits<double>::max();
    double h_max = std::numeric_limits<double>::lowest();
    
    std::vector<double> h_values;
    if (config_.enable_statistical_analysis) {
        h_values.reserve(Nx_ * Ny_);
    }
    
    // Calculate H-function at each grid point using entropic form
    for (int x = 0; x < Nx_; ++x) {
        for (int y = 0; y < Ny_; ++y) {
            // For entropic BGK, use logarithmic form of H-function
            double local_h = 0.0;
            for (int q = 0; q < Q; ++q) {
                if (f[q][x][y] > 1e-15) {
                    local_h += f[q][x][y] * std::log(f[q][x][y] / w[q]);
                }
            }
            
            // Check for numerical anomalies
            if (std::isnan(local_h) || std::isinf(local_h)) {
                metrics.has_nan_values = true;
                local_h = 0.0;
            }
            
            total_h += local_h;
            total_mass += rho[x][y];
            total_kinetic_energy += 0.5 * rho[x][y] * (ux[x][y] * ux[x][y] + uy[x][y] * uy[x][y]);
            
            h_min = std::min(h_min, local_h);
            h_max = std::max(h_max, local_h);
            
            if (config_.enable_statistical_analysis) {
                h_values.push_back(local_h);
            }
            
            h_field_flat_[x * Ny_ + y] = local_h;
        }
    }
    
    // Set basic metrics
    metrics.h_function_value = total_h;
    metrics.total_mass = total_mass;
    metrics.kinetic_energy = total_kinetic_energy;
    metrics.h_function_min = h_min;
    metrics.h_function_max = h_max;
    
    // Calculate H-function change and monotonicity
    if (has_previous_data_) {
        metrics.h_function_change = total_h - previous_h_value_;
        // For entropic BGK, H should decrease monotonically
        metrics.monotonic_decrease = (metrics.h_function_change <= config_.monotonicity_tolerance);
        
        if (!metrics.monotonic_decrease && config_.log_violations) {
            logMonotonicityViolation(timestep, total_h, previous_h_value_);
        }
        
        metrics.entropy_production = calculateEntropyProductionRate(total_h, previous_h_value_, 1.0);
    } else {
        metrics.h_function_change = 0.0;
        metrics.monotonic_decrease = true;
        metrics.entropy_production = 0.0;
    }
    
    // Calculate spatial statistics
    if (config_.enable_statistical_analysis && !h_values.empty()) {
        double mean = std::accumulate(h_values.begin(), h_values.end(), 0.0) / h_values.size();
        double variance = 0.0;
        for (double h : h_values) {
            variance += (h - mean) * (h - mean);
        }
        variance /= h_values.size();
        metrics.h_function_variance = variance;
    }
    
    // Update previous data
    previous_h_value_ = total_h;
    has_previous_data_ = true;
    last_analysis_timestep_ = timestep;
    
    return metrics;
}

void HTheoremAnalyzer::recordHEvolution(const HTheoremMetrics& metrics) {
    evolution_data_.timesteps.push_back(metrics.timestep);
    evolution_data_.h_values.push_back(metrics.h_function_value);
    evolution_data_.h_changes.push_back(metrics.h_function_change);
    evolution_data_.entropy_production_rates.push_back(metrics.entropy_production);
    evolution_data_.monotonicity_violations.push_back(!metrics.monotonic_decrease);
    
    if (config_.enable_statistical_analysis) {
        evolution_data_.h_variances.push_back(metrics.h_function_variance);
        evolution_data_.kinetic_energies.push_back(metrics.kinetic_energy);
        evolution_data_.total_masses.push_back(metrics.total_mass);
    }
    
    // Store metrics for history
    metrics_history_.push_back(metrics);
    
    // Update evolution statistics
    updateEvolutionStatistics();
}

bool HTheoremAnalyzer::verifyMonotonicDecrease(double tolerance) const {
    if (evolution_data_.h_values.size() < 2) {
        return true; // Not enough data to verify
    }
    
    for (size_t i = 1; i < evolution_data_.h_values.size(); ++i) {
        double change = evolution_data_.h_values[i] - evolution_data_.h_values[i-1];
        if (change > tolerance) {
            return false;
        }
    }
    
    return true;
}

HTheoremAnalyzer::HTheoremMetrics HTheoremAnalyzer::getLatestMetrics() const {
    if (metrics_history_.empty()) {
        return HTheoremMetrics{}; // Return default-constructed metrics
    }
    return metrics_history_.back();
}

double HTheoremAnalyzer::calculateEntropyProductionRate(double current_h, double previous_h, double dt) const {
    return -(current_h - previous_h) / dt; // Negative because H should decrease
}

HTheoremAnalyzer::SpatialStatistics HTheoremAnalyzer::analyzeSpatialDistribution(double** h_field) const {
    SpatialStatistics stats;
    
    std::vector<double> values;
    values.reserve(Nx_ * Ny_);
    
    for (int x = 0; x < Nx_; ++x) {
        for (int y = 0; y < Ny_; ++y) {
            values.push_back(h_field[x][y]);
        }
    }
    
    if (values.empty()) {
        return stats;
    }
    
    // Calculate mean
    stats.mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
    
    // Calculate variance
    double variance_sum = 0.0;
    for (double val : values) {
        variance_sum += (val - stats.mean) * (val - stats.mean);
    }
    stats.variance = variance_sum / values.size();
    stats.standard_deviation = std::sqrt(stats.variance);
    
    // Find min and max
    auto minmax = std::minmax_element(values.begin(), values.end());
    stats.min_value = *minmax.first;
    stats.max_value = *minmax.second;
    
    return stats;
}

bool HTheoremAnalyzer::writeEvolutionData(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }
    
    // Write header
    file << "# H-theorem evolution data\n";
    file << "# Columns: timestep, h_value, h_change, entropy_production, monotonic_violation";
    if (config_.enable_statistical_analysis) {
        file << ", h_variance, kinetic_energy, total_mass";
    }
    file << "\n";
    
    // Write data
    for (size_t i = 0; i < evolution_data_.timesteps.size(); ++i) {
        file << evolution_data_.timesteps[i] << ","
             << std::scientific << std::setprecision(12)
             << evolution_data_.h_values[i] << ","
             << evolution_data_.h_changes[i] << ","
             << evolution_data_.entropy_production_rates[i] << ","
             << (evolution_data_.monotonicity_violations[i] ? 1 : 0);
        
        if (config_.enable_statistical_analysis && i < evolution_data_.h_variances.size()) {
            file << "," << evolution_data_.h_variances[i]
                 << "," << evolution_data_.kinetic_energies[i]
                 << "," << evolution_data_.total_masses[i];
        }
        file << "\n";
    }
    
    return true;
}

std::string HTheoremAnalyzer::generateDiagnosticReport() const {
    std::ostringstream report;
    
    report << "=== H-theorem Analysis Diagnostic Report ===\n\n";
    
    if (evolution_data_.timesteps.empty()) {
        report << "No data available for analysis.\n";
        return report.str();
    }
    
    // Basic statistics
    report << "Data Summary:\n";
    report << "  Total timesteps analyzed: " << evolution_data_.timesteps.size() << "\n";
    report << "  Timestep range: " << evolution_data_.timesteps.front() 
           << " to " << evolution_data_.timesteps.back() << "\n";
    
    // H-function evolution
    report << "\nH-function Evolution:\n";
    report << "  Initial H-value: " << std::scientific << std::setprecision(6) 
           << evolution_data_.h_values.front() << "\n";
    report << "  Final H-value: " << evolution_data_.h_values.back() << "\n";
    report << "  Total change: " << (evolution_data_.h_values.back() - evolution_data_.h_values.front()) << "\n";
    
    // Monotonicity analysis
    int violation_count = std::count(evolution_data_.monotonicity_violations.begin(),
                                   evolution_data_.monotonicity_violations.end(), true);
    report << "\nMonotonicity Analysis:\n";
    report << "  Monotonicity violations: " << violation_count << "\n";
    report << "  Violation percentage: " << std::fixed << std::setprecision(2)
           << (100.0 * violation_count / evolution_data_.timesteps.size()) << "%\n";
    report << "  Overall monotonic: " << (violation_count == 0 ? "YES" : "NO") << "\n";
    
    // Entropy production
    if (!evolution_data_.entropy_production_rates.empty()) {
        double avg_entropy_production = std::accumulate(
            evolution_data_.entropy_production_rates.begin(),
            evolution_data_.entropy_production_rates.end(), 0.0) / 
            evolution_data_.entropy_production_rates.size();
        
        report << "\nEntropy Production:\n";
        report << "  Average entropy production rate: " << std::scientific 
               << avg_entropy_production << "\n";
    }
    
    // Statistical analysis
    if (config_.enable_statistical_analysis && !evolution_data_.h_variances.empty()) {
        report << "\nStatistical Analysis:\n";
        report << "  Average H-function variance: " << std::scientific
               << std::accumulate(evolution_data_.h_variances.begin(),
                                evolution_data_.h_variances.end(), 0.0) / 
                                evolution_data_.h_variances.size() << "\n";
    }
    
    return report.str();
}

void HTheoremAnalyzer::reset() {
    evolution_data_ = HEvolutionData{};
    metrics_history_.clear();
    previous_h_value_ = 0.0;
    has_previous_data_ = false;
    last_analysis_timestep_ = -1;
}

double HTheoremAnalyzer::calculateEquilibrium(int direction, double rho, double ux, double uy) const {
    double uxeq = ux * ex[direction];
    double uyeq = uy * ey[direction];
    double u2 = ux * ux + uy * uy;
    double usq = uxeq + uyeq;
    
    return w[direction] * rho * (1.0 + usq / cs2 + 0.5 * usq * usq / (cs2 * cs2) - 0.5 * u2 / cs2);
}

double HTheoremAnalyzer::calculateLocalHFunction(const double* f, const double* feq) const {
    double h = 0.0;
    for (int q = 0; q < Q; ++q) {
        if (f[q] > 1e-15 && feq[q] > 1e-15) {
            h += f[q] * std::log(f[q] / feq[q]);
        }
    }
    return h;
}

double HTheoremAnalyzer::calculateLocalHFunctionLog(const double* f) const {
    double h = 0.0;
    for (int q = 0; q < Q; ++q) {
        if (f[q] > 1e-15) {
            h += f[q] * std::log(f[q] / w[q]);
        }
    }
    return h;
}

bool HTheoremAnalyzer::checkNumericalAnomalies(double** h_field) const {
    for (int x = 0; x < Nx_; ++x) {
        for (int y = 0; y < Ny_; ++y) {
            if (std::isnan(h_field[x][y]) || std::isinf(h_field[x][y])) {
                return true;
            }
        }
    }
    return false;
}

void HTheoremAnalyzer::updateEvolutionStatistics() {
    if (evolution_data_.h_values.size() < 2) {
        return;
    }
    
    // Update violation count
    evolution_data_.monotonicity_violation_count = std::count(
        evolution_data_.monotonicity_violations.begin(),
        evolution_data_.monotonicity_violations.end(), true);
    
    // Update overall monotonicity
    evolution_data_.overall_monotonic = (evolution_data_.monotonicity_violation_count == 0);
    
    // Calculate average entropy production
    if (!evolution_data_.entropy_production_rates.empty()) {
        evolution_data_.average_entropy_production = std::accumulate(
            evolution_data_.entropy_production_rates.begin(),
            evolution_data_.entropy_production_rates.end(), 0.0) / 
            evolution_data_.entropy_production_rates.size();
    }
    
    // Estimate H-function decay rate (simple linear fit)
    if (evolution_data_.h_values.size() >= 10) {
        double initial_h = evolution_data_.h_values.front();
        double final_h = evolution_data_.h_values.back();
        double time_span = evolution_data_.timesteps.back() - evolution_data_.timesteps.front();
        
        if (time_span > 0 && initial_h != 0) {
            evolution_data_.h_function_decay_rate = (final_h - initial_h) / (initial_h * time_span);
        }
    }
}

void HTheoremAnalyzer::logMonotonicityViolation(int timestep, double h_current, double h_previous) const {
    if (config_.log_violations) {
        std::cout << "WARNING: H-theorem monotonicity violation at timestep " << timestep
                  << ": H increased from " << std::scientific << h_previous
                  << " to " << h_current << " (change: " << (h_current - h_previous) << ")\n";
    }
}

} // namespace analysis
} // namespace lbm