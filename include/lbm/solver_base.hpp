#pragma once

#include <memory>
#include <string>
#include <vector>
#include "stability_monitor.hpp"
#include "../utils/validation.hpp"

namespace lbm {

/**
 * @brief Base class for all LBM solvers providing common interface and functionality
 * 
 * This abstract base class defines the common interface for all LBM solvers
 * and provides shared functionality like stability monitoring and validation.
 */
class SolverBase {
public:
    /**
     * @brief Virtual destructor for proper cleanup of derived classes
     */
    virtual ~SolverBase() = default;

    /**
     * @brief Initialize the solver with grid and physical parameters
     * @throws std::runtime_error if initialization fails
     */
    virtual void initialize() = 0;

    /**
     * @brief Perform one time step of the simulation
     * @throws std::runtime_error if step fails due to numerical instability
     */
    virtual void step() = 0;

    /**
     * @brief Check if the simulation is numerically stable
     * @return true if stable, false otherwise
     */
    virtual bool isStable() const = 0;

    /**
     * @brief Write simulation output data
     * @param timestep Current simulation timestep
     * @throws std::runtime_error if output writing fails
     */
    virtual void writeOutput(int timestep) = 0;

    /**
     * @brief Get the current simulation time step
     * @return Current timestep
     */
    int getCurrentTimestep() const { return current_timestep_; }

    /**
     * @brief Get the grid dimensions
     * @return Pair of (Nx, Ny) grid dimensions
     */
    std::pair<int, int> getGridDimensions() const { return {Nx_, Ny_}; }

    /**
     * @brief Get stability monitoring results
     * @return Reference to stability monitor
     */
    const StabilityMonitor& getStabilityMonitor() const { return stability_monitor_; }

    /**
     * @brief Enable or disable stability monitoring
     * @param enabled Whether to enable stability monitoring
     */
    void setStabilityMonitoring(bool enabled) { stability_monitoring_enabled_ = enabled; }

    /**
     * @brief Set the output directory for simulation data
     * @param output_dir Path to output directory
     */
    void setOutputDirectory(const std::string& output_dir) { output_directory_ = output_dir; }

    // Public methods for testing
    /**
     * @brief Check stability (public version for testing)
     * @return true if stable, false if unstable
     */
    bool testCheckStability() { return checkAndHandleStability(); }

    /**
     * @brief Validate data (public version for testing)
     * @param data_description Description of data being validated
     * @return true if validation passes
     */
    bool testValidateData(const std::string& data_description) { return validateData(data_description); }

protected:
    /**
     * @brief Protected constructor for derived classes
     * @param Nx Grid size in x-direction
     * @param Ny Grid size in y-direction
     */
    SolverBase(int Nx, int Ny);

    /**
     * @brief Update the current timestep counter
     * @param timestep New timestep value
     */
    void updateTimestep(int timestep) { current_timestep_ = timestep; }

    /**
     * @brief Get the current timestep (protected non-const version)
     * @return Current timestep value
     */
    int getCurrentTimestep() { return current_timestep_; }

    /**
     * @brief Check stability and handle instabilities
     * @return true if stable, false if unstable
     */
    bool checkAndHandleStability();

    /**
     * @brief Validate simulation data
     * @param data_description Description of data being validated
     * @return true if validation passes
     */
    bool validateData(const std::string& data_description);

    // Grid dimensions
    int Nx_, Ny_;
    
    // Current simulation state
    int current_timestep_;
    
    // Monitoring and validation
    StabilityMonitor stability_monitor_;
    DataValidator validator_;
    bool stability_monitoring_enabled_;
    
    // Output configuration
    std::string output_directory_;
    
    // LBM constants
    static constexpr int Q = 9;  // Number of velocity directions
    static constexpr double cs2 = 1.0 / 3.0;  // Speed of sound squared
    
    // Velocity directions (D2Q9)
    static constexpr int ex[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    static constexpr int ey[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
    static constexpr double w[Q] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 
                                   1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
};

} // namespace lbm