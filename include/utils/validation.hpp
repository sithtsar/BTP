#pragma once

#include "lbm/stability_monitor.hpp"
#include <string>
#include <vector>
#include <functional>
#include <memory>

namespace lbm {

/**
 * @brief Validation severity levels
 */
enum class ValidationLevel {
    BASIC,          // Basic checks (NaN, infinity, bounds)
    COMPREHENSIVE,  // Includes physics consistency checks
    STRICT          // Includes advanced mathematical validations
};

/**
 * @brief Validation result status
 */
enum class ValidationStatus {
    PASSED,         // All validations passed
    WARNING,        // Minor issues detected
    FAILED,         // Critical issues detected
    ERROR           // Validation process failed
};

/**
 * @brief Individual validation rule result
 */
struct ValidationRuleResult {
    std::string rule_name;
    ValidationStatus status;
    std::string message;
    double severity_score;  // 0.0 (no issue) to 1.0 (critical)
    
    ValidationRuleResult(const std::string& name, ValidationStatus stat, 
                        const std::string& msg, double score = 0.0)
        : rule_name(name), status(stat), message(msg), severity_score(score) {}
};

/**
 * @brief Comprehensive validation result
 */
struct ValidationResult {
    ValidationStatus overall_status;
    std::vector<ValidationRuleResult> rule_results;
    double overall_severity_score;
    std::string summary_message;
    int warnings_count;
    int failures_count;
    int errors_count;
    
    ValidationResult() 
        : overall_status(ValidationStatus::PASSED), overall_severity_score(0.0),
          warnings_count(0), failures_count(0), errors_count(0) {}
    
    bool isValid() const {
        return overall_status == ValidationStatus::PASSED || 
               overall_status == ValidationStatus::WARNING;
    }
    
    bool hasCriticalIssues() const {
        return overall_status == ValidationStatus::FAILED || 
               overall_status == ValidationStatus::ERROR;
    }
};

/**
 * @brief Configuration for validation rules
 */
struct ValidationRules {
    // Basic validation parameters
    double max_velocity = 0.1;
    double min_density = 1e-6;
    double max_density = 1e6;
    double mass_conservation_tolerance = 1e-10;
    
    // Physics validation parameters
    double max_mach_number = 0.1;              // Maximum Mach number for LBM validity
    double min_reynolds_number = 1.0;          // Minimum Reynolds number
    double max_reynolds_number = 1000.0;       // Maximum Reynolds number
    double viscosity_tolerance = 1e-8;         // Tolerance for viscosity calculations
    
    // Advanced validation parameters
    double equilibrium_tolerance = 1e-6;       // Tolerance for equilibrium distribution
    double symmetry_tolerance = 1e-8;          // Tolerance for symmetry checks
    double momentum_conservation_tolerance = 1e-10; // Momentum conservation tolerance
    
    // Validation control flags
    bool check_nan = true;
    bool check_inf = true;
    bool check_bounds = true;
    bool check_mass_conservation = true;
    bool check_momentum_conservation = true;
    bool check_physics_consistency = true;
    bool check_equilibrium_validity = true;
    bool check_symmetry_properties = true;
    bool check_boundary_conditions = true;
    
    // Performance and memory checks
    bool check_memory_usage = false;
    bool check_computational_efficiency = false;
    size_t max_memory_usage_mb = 1024;         // Maximum memory usage in MB
    double max_computation_time_ratio = 2.0;   // Max computation time vs theoretical
};

/**
 * @brief Custom validation rule function type
 */
using ValidationRuleFunction = std::function<ValidationRuleResult(const LatticeData&, const ValidationRules&)>;

/**
 * @brief Comprehensive data validation system for LBM simulations
 * 
 * This class provides a flexible validation pipeline with configurable rules
 * and different validation levels. It can perform basic stability checks,
 * physics consistency validation, and advanced mathematical verifications.
 */
class DataValidator {
public:
    /**
     * @brief Constructor with validation rules
     * @param rules Validation configuration
     */
    explicit DataValidator(const ValidationRules& rules = ValidationRules{});
    
    /**
     * @brief Destructor
     */
    ~DataValidator() = default;
    
    /**
     * @brief Validate lattice data with specified validation level
     * @param data Lattice data to validate
     * @param level Validation level (BASIC, COMPREHENSIVE, STRICT)
     * @param timestep Current timestep (optional, for context)
     * @return Comprehensive validation result
     */
    ValidationResult validate(const LatticeData& data, 
                            ValidationLevel level = ValidationLevel::COMPREHENSIVE,
                            int timestep = -1);
    
    /**
     * @brief Add custom validation rule
     * @param rule_name Name of the custom rule
     * @param rule_function Function implementing the validation logic
     * @param level Minimum validation level for this rule
     */
    void addCustomRule(const std::string& rule_name, 
                      ValidationRuleFunction rule_function,
                      ValidationLevel level = ValidationLevel::COMPREHENSIVE);
    
    /**
     * @brief Remove custom validation rule
     * @param rule_name Name of the rule to remove
     */
    void removeCustomRule(const std::string& rule_name);
    
    /**
     * @brief Update validation rules
     * @param rules New validation configuration
     */
    void updateRules(const ValidationRules& rules);
    
    /**
     * @brief Get current validation rules
     * @return Current validation configuration
     */
    const ValidationRules& getRules() const;
    
    /**
     * @brief Generate detailed diagnostic report
     * @param result Validation result to analyze
     * @return Formatted diagnostic report
     */
    std::string generateDiagnosticReport(const ValidationResult& result) const;
    
    /**
     * @brief Generate summary report for multiple validation results
     * @param results Vector of validation results
     * @return Summary report
     */
    std::string generateSummaryReport(const std::vector<ValidationResult>& results) const;
    
    /**
     * @brief Check if validation should be performed at this timestep
     * @param timestep Current timestep
     * @param check_interval Interval for validation checks
     * @return True if validation should be performed
     */
    static bool shouldValidate(int timestep, int check_interval = 100);
    
    /**
     * @brief Create validation rules from stability configuration
     * @param stability_config Stability monitor configuration
     * @return Corresponding validation rules
     */
    static ValidationRules fromStabilityConfig(const StabilityConfig& stability_config);

private:
    ValidationRules rules_;
    std::vector<std::pair<std::string, ValidationRuleFunction>> custom_rules_;
    std::vector<ValidationLevel> custom_rule_levels_;
    
    // Built-in validation rule implementations
    ValidationRuleResult validateBasicStability(const LatticeData& data) const;
    ValidationRuleResult validateNaNValues(const LatticeData& data) const;
    ValidationRuleResult validateInfinityValues(const LatticeData& data) const;
    ValidationRuleResult validateDensityBounds(const LatticeData& data) const;
    ValidationRuleResult validateVelocityBounds(const LatticeData& data) const;
    ValidationRuleResult validateMassConservation(const LatticeData& data) const;
    ValidationRuleResult validateMomentumConservation(const LatticeData& data) const;
    ValidationRuleResult validatePhysicsConsistency(const LatticeData& data) const;
    ValidationRuleResult validateMachNumber(const LatticeData& data) const;
    ValidationRuleResult validateReynoldsNumber(const LatticeData& data) const;
    ValidationRuleResult validateEquilibriumDistribution(const LatticeData& data) const;
    ValidationRuleResult validateSymmetryProperties(const LatticeData& data) const;
    ValidationRuleResult validateBoundaryConditions(const LatticeData& data) const;
    
    // Utility methods
    double calculateMachNumber(const LatticeData& data) const;
    double calculateReynoldsNumber(const LatticeData& data) const;
    double calculateViscosity(const LatticeData& data) const;
    bool checkEquilibriumDistribution(const LatticeData& data, int i, int j) const;
    double calculateEquilibrium(int direction, double density, double vx, double vy) const;
    
    // Result processing methods
    void processValidationResults(ValidationResult& result) const;
    ValidationStatus determineOverallStatus(const std::vector<ValidationRuleResult>& results) const;
    double calculateOverallSeverity(const std::vector<ValidationRuleResult>& results) const;
    std::string generateSummaryMessage(const ValidationResult& result) const;
    
    // Constants for LBM calculations
    static constexpr double LATTICE_WEIGHTS[9] = {
        4.0/9.0,  1.0/9.0,  1.0/9.0,  1.0/9.0,  1.0/9.0,
        1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0
    };
    
    static constexpr int LATTICE_VELOCITIES[9][2] = {
        {0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1},
        {1, 1}, {-1, 1}, {-1, -1}, {1, -1}
    };
    
    static constexpr double SPEED_OF_SOUND_SQUARED = 1.0/3.0;
};

} // namespace lbm