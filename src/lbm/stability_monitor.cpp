#include "lbm/stability_monitor.hpp"
#include <cmath>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <iostream>

namespace lbm {

StabilityMonitor::StabilityMonitor(const StabilityConfig& config)
    : config_(config), issue_count_(0), initial_mass_(0.0), initial_mass_set_(false) {
    stability_history_.reserve(1000); // Reserve space for efficiency
    error_log_.reserve(100);
}

StabilityMetrics StabilityMonitor::checkStability(const LatticeData& data, int timestep) {
    StabilityMetrics metrics;
    
    // Check for NaN values
    if (config_.check_nan) {
        metrics.nan_count = checkForNaN(data);
        metrics.no_nan_values = (metrics.nan_count == 0);
        if (!metrics.no_nan_values) {
            metrics.error_message += "NaN values detected (" + std::to_string(metrics.nan_count) + "); ";
        }
    }
    
    // Check for infinity values
    if (config_.check_inf) {
        metrics.inf_count = checkForInfinity(data);
        metrics.no_inf_values = (metrics.inf_count == 0);
        if (!metrics.no_inf_values) {
            metrics.error_message += "Infinity values detected (" + std::to_string(metrics.inf_count) + "); ";
        }
    }
    
    // Check density bounds
    if (config_.check_density_bounds) {
        checkDensityBounds(data, metrics);
    }
    
    // Check velocity bounds
    if (config_.check_velocity_bounds) {
        checkVelocityBounds(data, metrics);
    }
    
    // Check mass conservation
    if (config_.check_mass_conservation) {
        checkMassConservation(data, metrics);
    }
    
    // Log issues if any
    if (!metrics.isStable()) {
        issue_count_++;
        std::string severity = "WARNING";
        if (!metrics.no_nan_values || !metrics.no_inf_values) {
            severity = "CRITICAL";
        }
        logStabilityIssue(metrics.error_message, timestep, severity);
    }
    
    // Store in history if real-time monitoring is enabled
    if (config_.real_time_monitoring) {
        stability_history_.push_back(metrics);
        
        // Limit history size to prevent memory issues
        if (stability_history_.size() > 10000) {
            stability_history_.erase(stability_history_.begin(), 
                                   stability_history_.begin() + 1000);
        }
    }
    
    return metrics;
}

bool StabilityMonitor::shouldCheck(int timestep) const {
    if (timestep < 0) return true; // Always check if timestep not provided
    return (timestep % config_.check_interval == 0);
}

void StabilityMonitor::logStabilityIssue(const std::string& message, int timestep, 
                                        const std::string& severity) {
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    
    std::stringstream log_entry;
    log_entry << "[" << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S") << "] ";
    log_entry << "[" << severity << "] ";
    if (timestep >= 0) {
        log_entry << "[t=" << timestep << "] ";
    }
    log_entry << message;
    
    error_log_.push_back(log_entry.str());
    
    // Limit log size
    if (error_log_.size() > 1000) {
        error_log_.erase(error_log_.begin(), error_log_.begin() + 100);
    }
    
    // Print to console if detailed logging is enabled
    if (config_.detailed_logging) {
        std::cerr << log_entry.str() << std::endl;
    }
}

bool StabilityMonitor::hasStabilityIssues() const {
    return issue_count_ > 0;
}

int StabilityMonitor::getIssueCount() const {
    return issue_count_;
}

void StabilityMonitor::reset() {
    stability_history_.clear();
    error_log_.clear();
    issue_count_ = 0;
    initial_mass_ = 0.0;
    initial_mass_set_ = false;
}

void StabilityMonitor::updateConfig(const StabilityConfig& config) {
    config_ = config;
}

const StabilityConfig& StabilityMonitor::getConfig() const {
    return config_;
}

const std::vector<StabilityMetrics>& StabilityMonitor::getStabilityHistory() const {
    return stability_history_;
}

std::string StabilityMonitor::generateDiagnosticReport() const {
    std::stringstream report;
    
    report << "=== Stability Monitor Diagnostic Report ===\n";
    report << "Total issues detected: " << issue_count_ << "\n";
    report << "History entries: " << stability_history_.size() << "\n";
    report << "Log entries: " << error_log_.size() << "\n\n";
    
    // Configuration summary
    report << "Configuration:\n";
    report << "  Max velocity: " << config_.max_velocity << "\n";
    report << "  Min density: " << config_.min_density << "\n";
    report << "  Max density: " << config_.max_density << "\n";
    report << "  Mass conservation tolerance: " << config_.mass_conservation_tolerance << "\n";
    report << "  Check interval: " << config_.check_interval << "\n\n";
    
    // Recent stability metrics
    if (!stability_history_.empty()) {
        const auto& latest = stability_history_.back();
        report << "Latest Stability Metrics:\n";
        report << "  Stable: " << (latest.isStable() ? "YES" : "NO") << "\n";
        report << "  Max velocity: " << latest.max_velocity << "\n";
        report << "  Min density: " << latest.min_density << "\n";
        report << "  Max density: " << latest.max_density << "\n";
        report << "  Total mass: " << latest.total_mass << "\n";
        report << "  Mass change: " << latest.mass_change << "\n";
        report << "  NaN count: " << latest.nan_count << "\n";
        report << "  Inf count: " << latest.inf_count << "\n";
        report << "  Negative density count: " << latest.negative_density_count << "\n";
        report << "  Excessive velocity count: " << latest.excessive_velocity_count << "\n\n";
    }
    
    // Recent error log entries
    if (!error_log_.empty()) {
        report << "Recent Error Log (last 10 entries):\n";
        int start = std::max(0, static_cast<int>(error_log_.size()) - 10);
        for (int i = start; i < static_cast<int>(error_log_.size()); ++i) {
            report << "  " << error_log_[i] << "\n";
        }
    }
    
    return report.str();
}

int StabilityMonitor::checkForNaN(const LatticeData& data) const {
    int nan_count = 0;
    
    // Check density field
    for (int i = 0; i < data.nx; ++i) {
        for (int j = 0; j < data.ny; ++j) {
            if (isNaN(data.density[i][j])) nan_count++;
            if (isNaN(data.velocity_x[i][j])) nan_count++;
            if (isNaN(data.velocity_y[i][j])) nan_count++;
            
            // Check distribution functions
            for (int k = 0; k < 9; ++k) {
                if (isNaN(data.distribution[i][j][k])) nan_count++;
            }
        }
    }
    
    return nan_count;
}

int StabilityMonitor::checkForInfinity(const LatticeData& data) const {
    int inf_count = 0;
    
    // Check density field
    for (int i = 0; i < data.nx; ++i) {
        for (int j = 0; j < data.ny; ++j) {
            if (isInfinite(data.density[i][j])) inf_count++;
            if (isInfinite(data.velocity_x[i][j])) inf_count++;
            if (isInfinite(data.velocity_y[i][j])) inf_count++;
            
            // Check distribution functions
            for (int k = 0; k < 9; ++k) {
                if (isInfinite(data.distribution[i][j][k])) inf_count++;
            }
        }
    }
    
    return inf_count;
}

void StabilityMonitor::checkDensityBounds(const LatticeData& data, StabilityMetrics& metrics) const {
    metrics.min_density = std::numeric_limits<double>::max();
    metrics.max_density = std::numeric_limits<double>::lowest();
    metrics.negative_density_count = 0;
    
    for (int i = 0; i < data.nx; ++i) {
        for (int j = 0; j < data.ny; ++j) {
            double density = data.density[i][j];
            
            metrics.min_density = std::min(metrics.min_density, density);
            metrics.max_density = std::max(metrics.max_density, density);
            
            if (density < 0.0) {
                metrics.negative_density_count++;
            }
            
            if (density < config_.min_density || density > config_.max_density) {
                metrics.density_positive = false;
            }
        }
    }
    
    if (metrics.negative_density_count > 0) {
        metrics.density_positive = false;
        metrics.error_message += "Negative densities detected (" + 
                                std::to_string(metrics.negative_density_count) + "); ";
    }
    
    if (metrics.min_density < config_.min_density) {
        metrics.error_message += "Density below minimum (" + 
                                std::to_string(metrics.min_density) + " < " + 
                                std::to_string(config_.min_density) + "); ";
    }
    
    if (metrics.max_density > config_.max_density) {
        metrics.error_message += "Density above maximum (" + 
                                std::to_string(metrics.max_density) + " > " + 
                                std::to_string(config_.max_density) + "); ";
    }
}

void StabilityMonitor::checkVelocityBounds(const LatticeData& data, StabilityMetrics& metrics) const {
    metrics.max_velocity = 0.0;
    metrics.excessive_velocity_count = 0;
    
    for (int i = 0; i < data.nx; ++i) {
        for (int j = 0; j < data.ny; ++j) {
            double vx = data.velocity_x[i][j];
            double vy = data.velocity_y[i][j];
            double velocity_magnitude = std::sqrt(vx * vx + vy * vy);
            
            metrics.max_velocity = std::max(metrics.max_velocity, velocity_magnitude);
            
            if (velocity_magnitude > config_.max_velocity) {
                metrics.excessive_velocity_count++;
                metrics.velocity_bounded = false;
            }
        }
    }
    
    if (!metrics.velocity_bounded) {
        metrics.error_message += "Excessive velocities detected (" + 
                                std::to_string(metrics.excessive_velocity_count) + 
                                " points, max=" + std::to_string(metrics.max_velocity) + "); ";
    }
}

void StabilityMonitor::checkMassConservation(const LatticeData& data, StabilityMetrics& metrics) {
    metrics.total_mass = calculateTotalMass(data);
    
    if (!initial_mass_set_) {
        initial_mass_ = metrics.total_mass;
        initial_mass_set_ = true;
        metrics.mass_change = 0.0;
    } else {
        metrics.mass_change = std::abs(metrics.total_mass - initial_mass_);
        
        if (metrics.mass_change > config_.mass_conservation_tolerance) {
            metrics.mass_conserved = false;
            metrics.error_message += "Mass not conserved (change=" + 
                                    std::to_string(metrics.mass_change) + "); ";
        }
    }
}

double StabilityMonitor::calculateTotalMass(const LatticeData& data) const {
    double total_mass = 0.0;
    
    for (int i = 0; i < data.nx; ++i) {
        for (int j = 0; j < data.ny; ++j) {
            total_mass += data.density[i][j];
        }
    }
    
    return total_mass;
}

bool StabilityMonitor::isNaN(double value) const {
    return std::isnan(value);
}

bool StabilityMonitor::isInfinite(double value) const {
    return std::isinf(value);
}

} // namespace lbm