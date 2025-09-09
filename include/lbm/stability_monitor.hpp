#pragma once

#include <string>
#include <vector>
#include <chrono>
#include <memory>

namespace lbm {

/**
 * @brief Structure to hold lattice data for stability monitoring
 */
struct LatticeData {
    std::vector<std::vector<double>> density;     // 2D density field
    std::vector<std::vector<double>> velocity_x;  // 2D x-velocity field
    std::vector<std::vector<double>> velocity_y;  // 2D y-velocity field
    std::vector<std::vector<std::vector<double>>> distribution; // 3D distribution functions [x][y][direction]
    
    int nx, ny;  // Grid dimensions
    
    LatticeData(int nx, int ny) : nx(nx), ny(ny) {
        density.resize(nx, std::vector<double>(ny, 0.0));
        velocity_x.resize(nx, std::vector<double>(ny, 0.0));
        velocity_y.resize(nx, std::vector<double>(ny, 0.0));
        distribution.resize(nx, std::vector<std::vector<double>>(ny, std::vector<double>(9, 0.0)));
    }
};

/**
 * @brief Structure containing comprehensive stability metrics
 */
struct StabilityMetrics {
    bool density_positive = true;
    bool velocity_bounded = true;
    bool no_nan_values = true;
    bool no_inf_values = true;
    bool mass_conserved = true;
    
    double max_velocity = 0.0;
    double min_density = 0.0;
    double max_density = 0.0;
    double total_mass = 0.0;
    double mass_change = 0.0;
    
    int nan_count = 0;
    int inf_count = 0;
    int negative_density_count = 0;
    int excessive_velocity_count = 0;
    
    std::string error_message;
    
    bool isStable() const {
        return density_positive && velocity_bounded && no_nan_values && 
               no_inf_values && mass_conserved;
    }
};

/**
 * @brief Configuration parameters for stability monitoring
 */
struct StabilityConfig {
    double max_velocity = 0.1;           // Maximum allowed velocity magnitude
    double min_density = 1e-6;           // Minimum allowed density
    double max_density = 1e6;            // Maximum allowed density
    double mass_conservation_tolerance = 1e-10;  // Tolerance for mass conservation
    
    bool check_nan = true;               // Enable NaN detection
    bool check_inf = true;               // Enable infinity detection
    bool check_mass_conservation = true; // Enable mass conservation check
    bool check_velocity_bounds = true;   // Enable velocity bounds check
    bool check_density_bounds = true;    // Enable density bounds check
    
    int check_interval = 100;            // Check every N timesteps
    bool real_time_monitoring = true;    // Enable real-time monitoring
    bool detailed_logging = false;       // Enable detailed error logging
};

/**
 * @brief Comprehensive numerical stability monitoring system for LBM simulations
 * 
 * This class provides real-time monitoring of numerical stability indicators
 * including density positivity, velocity bounds, NaN/infinity detection,
 * and mass conservation verification.
 */
class StabilityMonitor {
public:
    /**
     * @brief Constructor with configuration
     * @param config Stability monitoring configuration
     */
    explicit StabilityMonitor(const StabilityConfig& config = StabilityConfig{});
    
    /**
     * @brief Destructor
     */
    ~StabilityMonitor() = default;
    
    /**
     * @brief Check stability of lattice data
     * @param data Lattice data to analyze
     * @param timestep Current simulation timestep
     * @return Comprehensive stability metrics
     */
    StabilityMetrics checkStability(const LatticeData& data, int timestep = -1);
    
    /**
     * @brief Check if monitoring should be performed at this timestep
     * @param timestep Current timestep
     * @return True if monitoring should be performed
     */
    bool shouldCheck(int timestep) const;
    
    /**
     * @brief Log a stability issue with detailed information
     * @param message Error message
     * @param timestep Current timestep
     * @param severity Severity level (INFO, WARNING, ERROR, CRITICAL)
     */
    void logStabilityIssue(const std::string& message, int timestep = -1, 
                          const std::string& severity = "ERROR");
    
    /**
     * @brief Check if any stability issues have been detected
     * @return True if stability issues exist
     */
    bool hasStabilityIssues() const;
    
    /**
     * @brief Get the number of stability issues detected
     * @return Number of issues
     */
    int getIssueCount() const;
    
    /**
     * @brief Reset stability monitoring state
     */
    void reset();
    
    /**
     * @brief Update configuration
     * @param config New configuration
     */
    void updateConfig(const StabilityConfig& config);
    
    /**
     * @brief Get current configuration
     * @return Current configuration
     */
    const StabilityConfig& getConfig() const;
    
    /**
     * @brief Get stability history
     * @return Vector of historical stability metrics
     */
    const std::vector<StabilityMetrics>& getStabilityHistory() const;
    
    /**
     * @brief Generate comprehensive diagnostic report
     * @return Detailed diagnostic string
     */
    std::string generateDiagnosticReport() const;

private:
    StabilityConfig config_;
    std::vector<StabilityMetrics> stability_history_;
    std::vector<std::string> error_log_;
    int issue_count_;
    double initial_mass_;
    bool initial_mass_set_;
    
    /**
     * @brief Check for NaN values in lattice data
     * @param data Lattice data to check
     * @return Number of NaN values found
     */
    int checkForNaN(const LatticeData& data) const;
    
    /**
     * @brief Check for infinity values in lattice data
     * @param data Lattice data to check
     * @return Number of infinity values found
     */
    int checkForInfinity(const LatticeData& data) const;
    
    /**
     * @brief Check density bounds and positivity
     * @param data Lattice data to check
     * @param metrics Metrics structure to update
     */
    void checkDensityBounds(const LatticeData& data, StabilityMetrics& metrics) const;
    
    /**
     * @brief Check velocity bounds
     * @param data Lattice data to check
     * @param metrics Metrics structure to update
     */
    void checkVelocityBounds(const LatticeData& data, StabilityMetrics& metrics) const;
    
    /**
     * @brief Check mass conservation
     * @param data Lattice data to check
     * @param metrics Metrics structure to update
     */
    void checkMassConservation(const LatticeData& data, StabilityMetrics& metrics);
    
    /**
     * @brief Calculate total mass in the system
     * @param data Lattice data
     * @return Total mass
     */
    double calculateTotalMass(const LatticeData& data) const;
    
    /**
     * @brief Check if a double value is NaN
     * @param value Value to check
     * @return True if NaN
     */
    bool isNaN(double value) const;
    
    /**
     * @brief Check if a double value is infinite
     * @param value Value to check
     * @return True if infinite
     */
    bool isInfinite(double value) const;
};

} // namespace lbm