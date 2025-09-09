#include "utils/validation.hpp"
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <numeric>

namespace lbm {

DataValidator::DataValidator(const ValidationRules& rules) : rules_(rules) {
    custom_rules_.reserve(10);
    custom_rule_levels_.reserve(10);
}

ValidationResult DataValidator::validate(const LatticeData& data, ValidationLevel level, int timestep) {
    ValidationResult result;
    
    try {
        // Basic validation (always performed)
        if (rules_.check_nan) {
            result.rule_results.push_back(validateNaNValues(data));
        }
        
        if (rules_.check_inf) {
            result.rule_results.push_back(validateInfinityValues(data));
        }
        
        if (rules_.check_bounds) {
            result.rule_results.push_back(validateDensityBounds(data));
            result.rule_results.push_back(validateVelocityBounds(data));
        }
        
        if (rules_.check_mass_conservation) {
            result.rule_results.push_back(validateMassConservation(data));
        }
        
        // Comprehensive validation
        if (level >= ValidationLevel::COMPREHENSIVE) {
            if (rules_.check_momentum_conservation) {
                result.rule_results.push_back(validateMomentumConservation(data));
            }
            
            if (rules_.check_physics_consistency) {
                result.rule_results.push_back(validatePhysicsConsistency(data));
                result.rule_results.push_back(validateMachNumber(data));
                result.rule_results.push_back(validateReynoldsNumber(data));
            }
            
            if (rules_.check_boundary_conditions) {
                result.rule_results.push_back(validateBoundaryConditions(data));
            }
        }
        
        // Strict validation
        if (level >= ValidationLevel::STRICT) {
            if (rules_.check_equilibrium_validity) {
                result.rule_results.push_back(validateEquilibriumDistribution(data));
            }
            
            if (rules_.check_symmetry_properties) {
                result.rule_results.push_back(validateSymmetryProperties(data));
            }
        }
        
        // Apply custom rules based on validation level
        for (size_t i = 0; i < custom_rules_.size(); ++i) {
            if (level >= custom_rule_levels_[i]) {
                try {
                    auto custom_result = custom_rules_[i].second(data, rules_);
                    result.rule_results.push_back(custom_result);
                } catch (const std::exception& e) {
                    result.rule_results.emplace_back(
                        custom_rules_[i].first,
                        ValidationStatus::ERROR,
                        "Custom rule execution failed: " + std::string(e.what()),
                        1.0
                    );
                }
            }
        }
        
        // Process results to determine overall status
        processValidationResults(result);
        
    } catch (const std::exception& e) {
        result.overall_status = ValidationStatus::ERROR;
        result.summary_message = "Validation process failed: " + std::string(e.what());
        result.errors_count = 1;
    }
    
    return result;
}

void DataValidator::addCustomRule(const std::string& rule_name, 
                                 ValidationRuleFunction rule_function,
                                 ValidationLevel level) {
    // Check if rule already exists
    auto it = std::find_if(custom_rules_.begin(), custom_rules_.end(),
                          [&rule_name](const auto& pair) { return pair.first == rule_name; });
    
    if (it != custom_rules_.end()) {
        // Update existing rule
        size_t index = std::distance(custom_rules_.begin(), it);
        custom_rules_[index].second = rule_function;
        custom_rule_levels_[index] = level;
    } else {
        // Add new rule
        custom_rules_.emplace_back(rule_name, rule_function);
        custom_rule_levels_.push_back(level);
    }
}

void DataValidator::removeCustomRule(const std::string& rule_name) {
    auto it = std::find_if(custom_rules_.begin(), custom_rules_.end(),
                          [&rule_name](const auto& pair) { return pair.first == rule_name; });
    
    if (it != custom_rules_.end()) {
        size_t index = std::distance(custom_rules_.begin(), it);
        custom_rules_.erase(it);
        custom_rule_levels_.erase(custom_rule_levels_.begin() + index);
    }
}

void DataValidator::updateRules(const ValidationRules& rules) {
    rules_ = rules;
}

const ValidationRules& DataValidator::getRules() const {
    return rules_;
}

std::string DataValidator::generateDiagnosticReport(const ValidationResult& result) const {
    std::stringstream report;
    
    report << "=== Data Validation Diagnostic Report ===\n";
    report << "Overall Status: ";
    
    switch (result.overall_status) {
        case ValidationStatus::PASSED:
            report << "PASSED";
            break;
        case ValidationStatus::WARNING:
            report << "WARNING";
            break;
        case ValidationStatus::FAILED:
            report << "FAILED";
            break;
        case ValidationStatus::ERROR:
            report << "ERROR";
            break;
    }
    
    report << "\nOverall Severity Score: " << std::fixed << std::setprecision(4) 
           << result.overall_severity_score << "\n";
    report << "Summary: " << result.summary_message << "\n\n";
    
    report << "Issue Counts:\n";
    report << "  Warnings: " << result.warnings_count << "\n";
    report << "  Failures: " << result.failures_count << "\n";
    report << "  Errors: " << result.errors_count << "\n\n";
    
    report << "Individual Rule Results:\n";
    for (const auto& rule_result : result.rule_results) {
        report << "  [" << rule_result.rule_name << "] ";
        
        switch (rule_result.status) {
            case ValidationStatus::PASSED:
                report << "PASSED";
                break;
            case ValidationStatus::WARNING:
                report << "WARNING";
                break;
            case ValidationStatus::FAILED:
                report << "FAILED";
                break;
            case ValidationStatus::ERROR:
                report << "ERROR";
                break;
        }
        
        report << " (severity: " << std::fixed << std::setprecision(3) 
               << rule_result.severity_score << ")\n";
        
        if (!rule_result.message.empty()) {
            report << "    " << rule_result.message << "\n";
        }
    }
    
    return report.str();
}

std::string DataValidator::generateSummaryReport(const std::vector<ValidationResult>& results) const {
    if (results.empty()) {
        return "No validation results to summarize.\n";
    }
    
    std::stringstream report;
    
    report << "=== Validation Summary Report ===\n";
    report << "Total validation runs: " << results.size() << "\n\n";
    
    // Count status occurrences
    int passed = 0, warnings = 0, failures = 0, errors = 0;
    double avg_severity = 0.0;
    
    for (const auto& result : results) {
        switch (result.overall_status) {
            case ValidationStatus::PASSED: passed++; break;
            case ValidationStatus::WARNING: warnings++; break;
            case ValidationStatus::FAILED: failures++; break;
            case ValidationStatus::ERROR: errors++; break;
        }
        avg_severity += result.overall_severity_score;
    }
    
    avg_severity /= results.size();
    
    report << "Status Distribution:\n";
    report << "  Passed: " << passed << " (" << std::fixed << std::setprecision(1) 
           << (100.0 * passed / results.size()) << "%)\n";
    report << "  Warnings: " << warnings << " (" << std::fixed << std::setprecision(1) 
           << (100.0 * warnings / results.size()) << "%)\n";
    report << "  Failures: " << failures << " (" << std::fixed << std::setprecision(1) 
           << (100.0 * failures / results.size()) << "%)\n";
    report << "  Errors: " << errors << " (" << std::fixed << std::setprecision(1) 
           << (100.0 * errors / results.size()) << "%)\n\n";
    
    report << "Average Severity Score: " << std::fixed << std::setprecision(4) 
           << avg_severity << "\n";
    
    return report.str();
}

bool DataValidator::shouldValidate(int timestep, int check_interval) {
    if (timestep < 0) return true;
    return (timestep % check_interval == 0);
}

ValidationRules DataValidator::fromStabilityConfig(const StabilityConfig& stability_config) {
    ValidationRules rules;
    
    rules.max_velocity = stability_config.max_velocity;
    rules.min_density = stability_config.min_density;
    rules.max_density = stability_config.max_density;
    rules.mass_conservation_tolerance = stability_config.mass_conservation_tolerance;
    
    rules.check_nan = stability_config.check_nan;
    rules.check_inf = stability_config.check_inf;
    rules.check_bounds = true;
    rules.check_mass_conservation = stability_config.check_mass_conservation;
    
    return rules;
}

// Private method implementations

ValidationRuleResult DataValidator::validateNaNValues(const LatticeData& data) const {
    int nan_count = 0;
    
    for (int i = 0; i < data.nx; ++i) {
        for (int j = 0; j < data.ny; ++j) {
            if (std::isnan(data.density[i][j])) nan_count++;
            if (std::isnan(data.velocity_x[i][j])) nan_count++;
            if (std::isnan(data.velocity_y[i][j])) nan_count++;
            
            for (int k = 0; k < 9; ++k) {
                if (std::isnan(data.distribution[i][j][k])) nan_count++;
            }
        }
    }
    
    if (nan_count == 0) {
        return ValidationRuleResult("NaN Detection", ValidationStatus::PASSED, 
                                   "No NaN values detected", 0.0);
    } else {
        double severity = std::min(1.0, nan_count / (data.nx * data.ny * 11.0));
        return ValidationRuleResult("NaN Detection", ValidationStatus::FAILED,
                                   "Detected " + std::to_string(nan_count) + " NaN values", 
                                   severity);
    }
}

ValidationRuleResult DataValidator::validateInfinityValues(const LatticeData& data) const {
    int inf_count = 0;
    
    for (int i = 0; i < data.nx; ++i) {
        for (int j = 0; j < data.ny; ++j) {
            if (std::isinf(data.density[i][j])) inf_count++;
            if (std::isinf(data.velocity_x[i][j])) inf_count++;
            if (std::isinf(data.velocity_y[i][j])) inf_count++;
            
            for (int k = 0; k < 9; ++k) {
                if (std::isinf(data.distribution[i][j][k])) inf_count++;
            }
        }
    }
    
    if (inf_count == 0) {
        return ValidationRuleResult("Infinity Detection", ValidationStatus::PASSED,
                                   "No infinity values detected", 0.0);
    } else {
        double severity = std::min(1.0, inf_count / (data.nx * data.ny * 11.0));
        return ValidationRuleResult("Infinity Detection", ValidationStatus::FAILED,
                                   "Detected " + std::to_string(inf_count) + " infinity values",
                                   severity);
    }
}

ValidationRuleResult DataValidator::validateDensityBounds(const LatticeData& data) const {
    double min_density = std::numeric_limits<double>::max();
    double max_density = std::numeric_limits<double>::lowest();
    int violations = 0;
    
    for (int i = 0; i < data.nx; ++i) {
        for (int j = 0; j < data.ny; ++j) {
            double density = data.density[i][j];
            min_density = std::min(min_density, density);
            max_density = std::max(max_density, density);
            
            if (density < rules_.min_density || density > rules_.max_density) {
                violations++;
            }
        }
    }
    
    if (violations == 0) {
        return ValidationRuleResult("Density Bounds", ValidationStatus::PASSED,
                                   "All densities within bounds [" + 
                                   std::to_string(min_density) + ", " + 
                                   std::to_string(max_density) + "]", 0.0);
    } else {
        double severity = std::min(1.0, violations / (data.nx * data.ny * 1.0));
        std::string message = std::to_string(violations) + " density violations. " +
                             "Range: [" + std::to_string(min_density) + ", " + 
                             std::to_string(max_density) + "], " +
                             "Expected: [" + std::to_string(rules_.min_density) + ", " + 
                             std::to_string(rules_.max_density) + "]";
        
        ValidationStatus status = (severity > 0.1) ? ValidationStatus::FAILED : ValidationStatus::WARNING;
        return ValidationRuleResult("Density Bounds", status, message, severity);
    }
}

ValidationRuleResult DataValidator::validateVelocityBounds(const LatticeData& data) const {
    double max_velocity = 0.0;
    int violations = 0;
    
    for (int i = 0; i < data.nx; ++i) {
        for (int j = 0; j < data.ny; ++j) {
            double vx = data.velocity_x[i][j];
            double vy = data.velocity_y[i][j];
            double velocity_magnitude = std::sqrt(vx * vx + vy * vy);
            
            max_velocity = std::max(max_velocity, velocity_magnitude);
            
            if (velocity_magnitude > rules_.max_velocity) {
                violations++;
            }
        }
    }
    
    if (violations == 0) {
        return ValidationRuleResult("Velocity Bounds", ValidationStatus::PASSED,
                                   "All velocities within bounds (max: " + 
                                   std::to_string(max_velocity) + ")", 0.0);
    } else {
        double severity = std::min(1.0, violations / (data.nx * data.ny * 1.0));
        std::string message = std::to_string(violations) + " velocity violations. " +
                             "Max velocity: " + std::to_string(max_velocity) + ", " +
                             "Limit: " + std::to_string(rules_.max_velocity);
        
        ValidationStatus status = (severity > 0.1) ? ValidationStatus::FAILED : ValidationStatus::WARNING;
        return ValidationRuleResult("Velocity Bounds", status, message, severity);
    }
}

ValidationRuleResult DataValidator::validateMassConservation(const LatticeData& data) const {
    static double initial_mass = -1.0;
    static bool first_call = true;
    
    double total_mass = 0.0;
    for (int i = 0; i < data.nx; ++i) {
        for (int j = 0; j < data.ny; ++j) {
            total_mass += data.density[i][j];
        }
    }
    
    if (first_call) {
        initial_mass = total_mass;
        first_call = false;
        return ValidationRuleResult("Mass Conservation", ValidationStatus::PASSED,
                                   "Initial mass established: " + std::to_string(total_mass), 0.0);
    }
    
    double mass_change = std::abs(total_mass - initial_mass);
    double relative_change = mass_change / initial_mass;
    
    if (mass_change <= rules_.mass_conservation_tolerance) {
        return ValidationRuleResult("Mass Conservation", ValidationStatus::PASSED,
                                   "Mass conserved (change: " + std::to_string(mass_change) + ")", 0.0);
    } else {
        double severity = std::min(1.0, relative_change * 1000.0); // Scale for visibility
        std::string message = "Mass not conserved. Change: " + std::to_string(mass_change) +
                             " (relative: " + std::to_string(relative_change) + ")";
        
        ValidationStatus status = (relative_change > 0.01) ? ValidationStatus::FAILED : ValidationStatus::WARNING;
        return ValidationRuleResult("Mass Conservation", status, message, severity);
    }
}

ValidationRuleResult DataValidator::validateMomentumConservation(const LatticeData& data) const {
    double total_momentum_x = 0.0;
    double total_momentum_y = 0.0;
    
    for (int i = 0; i < data.nx; ++i) {
        for (int j = 0; j < data.ny; ++j) {
            total_momentum_x += data.density[i][j] * data.velocity_x[i][j];
            total_momentum_y += data.density[i][j] * data.velocity_y[i][j];
        }
    }
    
    double momentum_magnitude = std::sqrt(total_momentum_x * total_momentum_x + 
                                        total_momentum_y * total_momentum_y);
    
    // For Poiseuille flow, momentum should be approximately conserved
    // This is a simplified check - in practice, boundary conditions affect momentum
    if (momentum_magnitude < rules_.momentum_conservation_tolerance * data.nx * data.ny) {
        return ValidationRuleResult("Momentum Conservation", ValidationStatus::PASSED,
                                   "Momentum approximately conserved", 0.0);
    } else {
        return ValidationRuleResult("Momentum Conservation", ValidationStatus::WARNING,
                                   "Momentum magnitude: " + std::to_string(momentum_magnitude), 0.3);
    }
}

ValidationRuleResult DataValidator::validatePhysicsConsistency(const LatticeData& data) const {
    // Check basic physics consistency
    int inconsistencies = 0;
    
    for (int i = 0; i < data.nx; ++i) {
        for (int j = 0; j < data.ny; ++j) {
            // Check if density and distribution functions are consistent
            double density_from_dist = 0.0;
            for (int k = 0; k < 9; ++k) {
                density_from_dist += data.distribution[i][j][k];
            }
            
            double density_error = std::abs(density_from_dist - data.density[i][j]);
            if (density_error > 1e-10) {
                inconsistencies++;
            }
        }
    }
    
    if (inconsistencies == 0) {
        return ValidationRuleResult("Physics Consistency", ValidationStatus::PASSED,
                                   "Density and distribution functions consistent", 0.0);
    } else {
        double severity = std::min(1.0, inconsistencies / (data.nx * data.ny * 1.0));
        return ValidationRuleResult("Physics Consistency", ValidationStatus::WARNING,
                                   std::to_string(inconsistencies) + " consistency violations", severity);
    }
}

ValidationRuleResult DataValidator::validateMachNumber(const LatticeData& data) const {
    double max_mach = calculateMachNumber(data);
    
    if (max_mach <= rules_.max_mach_number) {
        return ValidationRuleResult("Mach Number", ValidationStatus::PASSED,
                                   "Max Mach number: " + std::to_string(max_mach), 0.0);
    } else {
        double severity = std::min(1.0, max_mach / rules_.max_mach_number - 1.0);
        return ValidationRuleResult("Mach Number", ValidationStatus::WARNING,
                                   "Max Mach number (" + std::to_string(max_mach) + 
                                   ") exceeds limit (" + std::to_string(rules_.max_mach_number) + ")",
                                   severity);
    }
}

ValidationRuleResult DataValidator::validateReynoldsNumber(const LatticeData& data) const {
    double reynolds = calculateReynoldsNumber(data);
    
    if (reynolds >= rules_.min_reynolds_number && reynolds <= rules_.max_reynolds_number) {
        return ValidationRuleResult("Reynolds Number", ValidationStatus::PASSED,
                                   "Reynolds number: " + std::to_string(reynolds), 0.0);
    } else {
        double severity = 0.3; // Moderate severity for Reynolds number issues
        std::string message = "Reynolds number (" + std::to_string(reynolds) + 
                             ") outside expected range [" + std::to_string(rules_.min_reynolds_number) +
                             ", " + std::to_string(rules_.max_reynolds_number) + "]";
        return ValidationRuleResult("Reynolds Number", ValidationStatus::WARNING, message, severity);
    }
}

ValidationRuleResult DataValidator::validateEquilibriumDistribution(const LatticeData& data) const {
    int violations = 0;
    
    for (int i = 0; i < data.nx; ++i) {
        for (int j = 0; j < data.ny; ++j) {
            if (!checkEquilibriumDistribution(data, i, j)) {
                violations++;
            }
        }
    }
    
    if (violations == 0) {
        return ValidationRuleResult("Equilibrium Distribution", ValidationStatus::PASSED,
                                   "All distributions within equilibrium tolerance", 0.0);
    } else {
        double severity = std::min(1.0, violations / (data.nx * data.ny * 1.0));
        return ValidationRuleResult("Equilibrium Distribution", ValidationStatus::WARNING,
                                   std::to_string(violations) + " equilibrium violations", severity);
    }
}

ValidationRuleResult DataValidator::validateSymmetryProperties(const LatticeData& data) const {
    // Check for basic symmetry properties in the flow field
    // This is a simplified check - real symmetry validation would be more complex
    
    int center_x = data.nx / 2;
    int violations = 0;
    
    for (int i = 1; i < center_x; ++i) {
        for (int j = 0; j < data.ny; ++j) {
            int mirror_i = data.nx - 1 - i;
            
            double density_diff = std::abs(data.density[i][j] - data.density[mirror_i][j]);
            double vx_sum = std::abs(data.velocity_x[i][j] + data.velocity_x[mirror_i][j]); // Should be opposite
            double vy_diff = std::abs(data.velocity_y[i][j] - data.velocity_y[mirror_i][j]); // Should be same
            
            if (density_diff > rules_.symmetry_tolerance ||
                vx_sum > rules_.symmetry_tolerance ||
                vy_diff > rules_.symmetry_tolerance) {
                violations++;
            }
        }
    }
    
    if (violations == 0) {
        return ValidationRuleResult("Symmetry Properties", ValidationStatus::PASSED,
                                   "Flow field exhibits expected symmetry", 0.0);
    } else {
        double severity = std::min(1.0, violations / (center_x * data.ny * 1.0));
        return ValidationRuleResult("Symmetry Properties", ValidationStatus::WARNING,
                                   std::to_string(violations) + " symmetry violations", severity);
    }
}

ValidationRuleResult DataValidator::validateBoundaryConditions(const LatticeData& data) const {
    int violations = 0;
    
    // Check boundary conditions (simplified - assumes no-slip walls at top/bottom)
    for (int i = 0; i < data.nx; ++i) {
        // Bottom boundary (j = 0)
        if (std::abs(data.velocity_x[i][0]) > 1e-10 || std::abs(data.velocity_y[i][0]) > 1e-10) {
            violations++;
        }
        
        // Top boundary (j = ny-1)
        if (std::abs(data.velocity_x[i][data.ny-1]) > 1e-10 || std::abs(data.velocity_y[i][data.ny-1]) > 1e-10) {
            violations++;
        }
    }
    
    if (violations == 0) {
        return ValidationRuleResult("Boundary Conditions", ValidationStatus::PASSED,
                                   "Boundary conditions properly enforced", 0.0);
    } else {
        double severity = std::min(1.0, violations / (2.0 * data.nx));
        return ValidationRuleResult("Boundary Conditions", ValidationStatus::WARNING,
                                   std::to_string(violations) + " boundary condition violations", severity);
    }
}

// Utility method implementations

double DataValidator::calculateMachNumber(const LatticeData& data) const {
    double max_velocity = 0.0;
    
    for (int i = 0; i < data.nx; ++i) {
        for (int j = 0; j < data.ny; ++j) {
            double vx = data.velocity_x[i][j];
            double vy = data.velocity_y[i][j];
            double velocity_magnitude = std::sqrt(vx * vx + vy * vy);
            max_velocity = std::max(max_velocity, velocity_magnitude);
        }
    }
    
    double speed_of_sound = std::sqrt(SPEED_OF_SOUND_SQUARED);
    return max_velocity / speed_of_sound;
}

double DataValidator::calculateReynoldsNumber(const LatticeData& data) const {
    // Simplified Reynolds number calculation for channel flow
    double max_velocity = 0.0;
    
    for (int i = 0; i < data.nx; ++i) {
        for (int j = 0; j < data.ny; ++j) {
            double vx = data.velocity_x[i][j];
            double vy = data.velocity_y[i][j];
            double velocity_magnitude = std::sqrt(vx * vx + vy * vy);
            max_velocity = std::max(max_velocity, velocity_magnitude);
        }
    }
    
    double characteristic_length = data.ny; // Channel height
    double viscosity = calculateViscosity(data);
    
    return max_velocity * characteristic_length / viscosity;
}

double DataValidator::calculateViscosity(const LatticeData& data) const {
    // Simplified viscosity calculation - in practice this would depend on tau
    // For now, use a typical value for LBM
    return 1.0/6.0; // Corresponds to tau = 1.0
}

bool DataValidator::checkEquilibriumDistribution(const LatticeData& data, int i, int j) const {
    double density = data.density[i][j];
    double vx = data.velocity_x[i][j];
    double vy = data.velocity_y[i][j];
    
    for (int k = 0; k < 9; ++k) {
        double expected_eq = calculateEquilibrium(k, density, vx, vy);
        double actual = data.distribution[i][j][k];
        
        if (std::abs(actual - expected_eq) > rules_.equilibrium_tolerance) {
            return false;
        }
    }
    
    return true;
}

double DataValidator::calculateEquilibrium(int direction, double density, double vx, double vy) const {
    double cx = LATTICE_VELOCITIES[direction][0];
    double cy = LATTICE_VELOCITIES[direction][1];
    double weight = LATTICE_WEIGHTS[direction];
    
    double cu = cx * vx + cy * vy;
    double u_sqr = vx * vx + vy * vy;
    
    return weight * density * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u_sqr);
}

void DataValidator::processValidationResults(ValidationResult& result) const {
    result.overall_status = determineOverallStatus(result.rule_results);
    result.overall_severity_score = calculateOverallSeverity(result.rule_results);
    
    // Count different types of issues
    for (const auto& rule_result : result.rule_results) {
        switch (rule_result.status) {
            case ValidationStatus::WARNING:
                result.warnings_count++;
                break;
            case ValidationStatus::FAILED:
                result.failures_count++;
                break;
            case ValidationStatus::ERROR:
                result.errors_count++;
                break;
            default:
                break;
        }
    }
    
    result.summary_message = generateSummaryMessage(result);
}

ValidationStatus DataValidator::determineOverallStatus(const std::vector<ValidationRuleResult>& results) const {
    bool has_error = false;
    bool has_failure = false;
    bool has_warning = false;
    
    for (const auto& result : results) {
        switch (result.status) {
            case ValidationStatus::ERROR:
                has_error = true;
                break;
            case ValidationStatus::FAILED:
                has_failure = true;
                break;
            case ValidationStatus::WARNING:
                has_warning = true;
                break;
            default:
                break;
        }
    }
    
    if (has_error) return ValidationStatus::ERROR;
    if (has_failure) return ValidationStatus::FAILED;
    if (has_warning) return ValidationStatus::WARNING;
    return ValidationStatus::PASSED;
}

double DataValidator::calculateOverallSeverity(const std::vector<ValidationRuleResult>& results) const {
    if (results.empty()) return 0.0;
    
    double max_severity = 0.0;
    double avg_severity = 0.0;
    
    for (const auto& result : results) {
        max_severity = std::max(max_severity, result.severity_score);
        avg_severity += result.severity_score;
    }
    
    avg_severity /= results.size();
    
    // Return weighted combination of max and average severity
    return 0.7 * max_severity + 0.3 * avg_severity;
}

std::string DataValidator::generateSummaryMessage(const ValidationResult& result) const {
    std::stringstream message;
    
    switch (result.overall_status) {
        case ValidationStatus::PASSED:
            message << "All validation checks passed successfully";
            break;
        case ValidationStatus::WARNING:
            message << "Validation completed with " << result.warnings_count << " warning(s)";
            break;
        case ValidationStatus::FAILED:
            message << "Validation failed with " << result.failures_count << " critical issue(s)";
            if (result.warnings_count > 0) {
                message << " and " << result.warnings_count << " warning(s)";
            }
            break;
        case ValidationStatus::ERROR:
            message << "Validation process encountered " << result.errors_count << " error(s)";
            break;
    }
    
    return message.str();
}

} // namespace lbm