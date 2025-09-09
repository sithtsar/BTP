#include "../../include/lbm/solver_base.hpp"
#include <stdexcept>
#include <iostream>

namespace lbm {

SolverBase::SolverBase(int Nx, int Ny) 
    : Nx_(Nx), Ny_(Ny), current_timestep_(0), 
      stability_monitoring_enabled_(true), output_directory_("data") {
    
    if (Nx <= 0 || Ny <= 0) {
        throw std::invalid_argument("Grid dimensions must be positive");
    }
    
    // Initialize stability monitor with default settings
    StabilityConfig stability_config;
    stability_config.max_velocity = 0.1;
    stability_config.min_density = 1e-6;
    stability_config.max_density = 1e6;
    stability_config.check_interval = 100;
    stability_monitor_.updateConfig(stability_config);
    
    // Initialize data validator with default rules
    ValidationRules rules;
    rules.max_velocity = 0.1;
    rules.min_density = 1e-6;
    rules.max_density = 1e6;
    rules.check_nan = true;
    rules.check_inf = true;
    rules.check_mass_conservation = true;
    rules.mass_conservation_tolerance = 1e-10;
    validator_.updateRules(rules);
}

bool SolverBase::checkAndHandleStability() {
    if (!stability_monitoring_enabled_) {
        return true;
    }
    
    // This method should be called by derived classes with their specific data
    // The actual stability check implementation depends on the specific solver type
    
    if (stability_monitor_.hasStabilityIssues()) {
        std::cerr << "Stability issues detected at timestep " << current_timestep_ << std::endl;
        
        // Get the latest stability metrics for detailed reporting
        auto history = stability_monitor_.getStabilityHistory();
        if (!history.empty()) {
            auto metrics = history.back();
            
            if (!metrics.density_positive) {
                std::cerr << "  - Negative density detected (min: " << metrics.min_density << ")" << std::endl;
            }
            if (!metrics.velocity_bounded) {
                std::cerr << "  - Velocity out of bounds (max: " << metrics.max_velocity << ")" << std::endl;
            }
            if (!metrics.no_nan_values) {
                std::cerr << "  - NaN values detected" << std::endl;
            }
            if (!metrics.no_inf_values) {
                std::cerr << "  - Infinite values detected" << std::endl;
            }
        }
        
        return false;
    }
    
    return true;
}

bool SolverBase::validateData(const std::string& data_description) {
    // This is a placeholder - derived classes should implement specific validation
    // by calling validator_.validate() with their actual data
    
    std::cout << "Validating " << data_description << " at timestep " << current_timestep_ << std::endl;
    return true;
}

} // namespace lbm