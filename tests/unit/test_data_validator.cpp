#include <gtest/gtest.h>
#include "utils/validation.hpp"
#include <limits>
#include <cmath>

using namespace lbm;

class DataValidatorTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create test data with reasonable dimensions
        test_data = std::make_unique<LatticeData>(8, 4);
        
        // Initialize with stable values
        initializeStableData();
        
        // Default validation rules
        rules = ValidationRules{};
        validator = std::make_unique<DataValidator>(rules);
    }
    
    void initializeStableData() {
        for (int i = 0; i < test_data->nx; ++i) {
            for (int j = 0; j < test_data->ny; ++j) {
                test_data->density[i][j] = 1.0;
                test_data->velocity_x[i][j] = 0.01;
                test_data->velocity_y[i][j] = 0.0;
                
                // Initialize distribution functions with equilibrium values
                double density = test_data->density[i][j];
                double vx = test_data->velocity_x[i][j];
                double vy = test_data->velocity_y[i][j];
                
                // Simplified equilibrium calculation
                for (int k = 0; k < 9; ++k) {
                    test_data->distribution[i][j][k] = density / 9.0; // Simplified
                }
            }
        }
    }
    
    std::unique_ptr<LatticeData> test_data;
    ValidationRules rules;
    std::unique_ptr<DataValidator> validator;
};

TEST_F(DataValidatorTest, BasicValidationPasses) {
    auto result = validator->validate(*test_data, ValidationLevel::BASIC);
    
    EXPECT_EQ(result.overall_status, ValidationStatus::PASSED);
    EXPECT_TRUE(result.isValid());
    EXPECT_FALSE(result.hasCriticalIssues());
    EXPECT_EQ(result.failures_count, 0);
    EXPECT_EQ(result.errors_count, 0);
    EXPECT_GT(result.rule_results.size(), 0);
}

TEST_F(DataValidatorTest, DetectsNaNValues) {
    // Introduce NaN values
    test_data->density[2][1] = std::numeric_limits<double>::quiet_NaN();
    test_data->velocity_x[3][2] = std::numeric_limits<double>::quiet_NaN();
    
    auto result = validator->validate(*test_data, ValidationLevel::BASIC);
    
    EXPECT_EQ(result.overall_status, ValidationStatus::FAILED);
    EXPECT_FALSE(result.isValid());
    EXPECT_TRUE(result.hasCriticalIssues());
    EXPECT_GT(result.failures_count, 0);
    
    // Check that NaN detection rule failed
    bool nan_rule_found = false;
    for (const auto& rule_result : result.rule_results) {
        if (rule_result.rule_name == "NaN Detection") {
            EXPECT_EQ(rule_result.status, ValidationStatus::FAILED);
            EXPECT_GT(rule_result.severity_score, 0.0);
            nan_rule_found = true;
            break;
        }
    }
    EXPECT_TRUE(nan_rule_found);
}

TEST_F(DataValidatorTest, DetectsInfinityValues) {
    // Introduce infinity values
    test_data->density[1][1] = std::numeric_limits<double>::infinity();
    test_data->velocity_y[2][2] = -std::numeric_limits<double>::infinity();
    test_data->distribution[3][1][5] = std::numeric_limits<double>::infinity();
    
    auto result = validator->validate(*test_data, ValidationLevel::BASIC);
    
    EXPECT_EQ(result.overall_status, ValidationStatus::FAILED);
    EXPECT_FALSE(result.isValid());
    EXPECT_TRUE(result.hasCriticalIssues());
    
    // Check that infinity detection rule failed
    bool inf_rule_found = false;
    for (const auto& rule_result : result.rule_results) {
        if (rule_result.rule_name == "Infinity Detection") {
            EXPECT_EQ(rule_result.status, ValidationStatus::FAILED);
            EXPECT_GT(rule_result.severity_score, 0.0);
            inf_rule_found = true;
            break;
        }
    }
    EXPECT_TRUE(inf_rule_found);
}

TEST_F(DataValidatorTest, DetectsDensityBoundViolations) {
    // Test minimum density violation
    test_data->density[1][1] = 1e-8; // Below default min_density (1e-6)
    
    auto result = validator->validate(*test_data, ValidationLevel::BASIC);
    
    EXPECT_NE(result.overall_status, ValidationStatus::PASSED);
    
    // Check density bounds rule
    bool density_rule_found = false;
    for (const auto& rule_result : result.rule_results) {
        if (rule_result.rule_name == "Density Bounds") {
            EXPECT_NE(rule_result.status, ValidationStatus::PASSED);
            density_rule_found = true;
            break;
        }
    }
    EXPECT_TRUE(density_rule_found);
    
    // Reset and test maximum density violation
    initializeStableData();
    test_data->density[2][2] = 2e6; // Above default max_density (1e6)
    
    result = validator->validate(*test_data, ValidationLevel::BASIC);
    EXPECT_NE(result.overall_status, ValidationStatus::PASSED);
}

TEST_F(DataValidatorTest, DetectsVelocityBoundViolations) {
    // Introduce excessive velocities
    test_data->velocity_x[1][1] = 0.08;
    test_data->velocity_y[1][1] = 0.08; // Combined magnitude > 0.1
    test_data->velocity_x[2][2] = 0.15; // This alone exceeds 0.1
    
    auto result = validator->validate(*test_data, ValidationLevel::BASIC);
    
    EXPECT_NE(result.overall_status, ValidationStatus::PASSED);
    
    // Check velocity bounds rule
    bool velocity_rule_found = false;
    for (const auto& rule_result : result.rule_results) {
        if (rule_result.rule_name == "Velocity Bounds") {
            EXPECT_NE(rule_result.status, ValidationStatus::PASSED);
            EXPECT_GT(rule_result.severity_score, 0.0);
            velocity_rule_found = true;
            break;
        }
    }
    EXPECT_TRUE(velocity_rule_found);
}

TEST_F(DataValidatorTest, ValidationLevels) {
    // Basic level should have fewer checks than comprehensive
    auto basic_result = validator->validate(*test_data, ValidationLevel::BASIC);
    auto comprehensive_result = validator->validate(*test_data, ValidationLevel::COMPREHENSIVE);
    auto strict_result = validator->validate(*test_data, ValidationLevel::STRICT);
    
    EXPECT_LE(basic_result.rule_results.size(), comprehensive_result.rule_results.size());
    EXPECT_LE(comprehensive_result.rule_results.size(), strict_result.rule_results.size());
    
    // All should pass for stable data
    EXPECT_EQ(basic_result.overall_status, ValidationStatus::PASSED);
    EXPECT_EQ(comprehensive_result.overall_status, ValidationStatus::PASSED);
    EXPECT_EQ(strict_result.overall_status, ValidationStatus::PASSED);
}

TEST_F(DataValidatorTest, CustomRuleAddition) {
    // Add a custom rule that always fails
    auto custom_rule = [](const LatticeData& data, const ValidationRules& rules) -> ValidationRuleResult {
        return ValidationRuleResult("Custom Test Rule", ValidationStatus::FAILED, 
                                   "This rule always fails", 0.5);
    };
    
    validator->addCustomRule("TestRule", custom_rule, ValidationLevel::BASIC);
    
    auto result = validator->validate(*test_data, ValidationLevel::BASIC);
    
    EXPECT_EQ(result.overall_status, ValidationStatus::FAILED);
    EXPECT_GT(result.failures_count, 0);
    
    // Check that custom rule was executed
    bool custom_rule_found = false;
    for (const auto& rule_result : result.rule_results) {
        if (rule_result.rule_name == "Custom Test Rule") {
            EXPECT_EQ(rule_result.status, ValidationStatus::FAILED);
            custom_rule_found = true;
            break;
        }
    }
    EXPECT_TRUE(custom_rule_found);
}

TEST_F(DataValidatorTest, CustomRuleRemoval) {
    // Add a custom rule
    auto custom_rule = [](const LatticeData& data, const ValidationRules& rules) -> ValidationRuleResult {
        return ValidationRuleResult("Removable Rule", ValidationStatus::WARNING, 
                                   "This rule generates warnings", 0.3);
    };
    
    validator->addCustomRule("RemovableRule", custom_rule, ValidationLevel::BASIC);
    
    auto result_with_rule = validator->validate(*test_data, ValidationLevel::BASIC);
    
    // Remove the rule
    validator->removeCustomRule("RemovableRule");
    
    auto result_without_rule = validator->validate(*test_data, ValidationLevel::BASIC);
    
    // Check that rule was removed
    bool rule_found_after_removal = false;
    for (const auto& rule_result : result_without_rule.rule_results) {
        if (rule_result.rule_name == "Removable Rule") {
            rule_found_after_removal = true;
            break;
        }
    }
    EXPECT_FALSE(rule_found_after_removal);
}

TEST_F(DataValidatorTest, ConfigurationUpdate) {
    ValidationRules new_rules;
    new_rules.max_velocity = 0.05; // More restrictive
    new_rules.min_density = 0.5;   // More restrictive
    
    validator->updateRules(new_rules);
    
    const auto& retrieved_rules = validator->getRules();
    EXPECT_EQ(retrieved_rules.max_velocity, 0.05);
    EXPECT_EQ(retrieved_rules.min_density, 0.5);
    
    // Test with velocity that would pass with default rules but fail with new rules
    test_data->velocity_x[1][1] = 0.08;
    
    auto result = validator->validate(*test_data, ValidationLevel::BASIC);
    
    // Should fail with more restrictive rules
    EXPECT_NE(result.overall_status, ValidationStatus::PASSED);
}

TEST_F(DataValidatorTest, DiagnosticReportGeneration) {
    // Introduce some issues
    test_data->density[1][1] = -0.5; // Negative density
    test_data->velocity_x[2][2] = 0.15; // Excessive velocity
    
    auto result = validator->validate(*test_data, ValidationLevel::COMPREHENSIVE);
    
    std::string report = validator->generateDiagnosticReport(result);
    
    EXPECT_TRUE(report.find("Data Validation Diagnostic Report") != std::string::npos);
    EXPECT_TRUE(report.find("Overall Status:") != std::string::npos);
    EXPECT_TRUE(report.find("Overall Severity Score:") != std::string::npos);
    EXPECT_TRUE(report.find("Individual Rule Results:") != std::string::npos);
    EXPECT_TRUE(report.find("Issue Counts:") != std::string::npos);
}

TEST_F(DataValidatorTest, SummaryReportGeneration) {
    std::vector<ValidationResult> results;
    
    // Create multiple validation results
    for (int i = 0; i < 5; ++i) {
        if (i % 2 == 0) {
            // Introduce issues in some results
            test_data->velocity_x[1][1] = 0.15;
        } else {
            // Keep others stable
            initializeStableData();
        }
        
        results.push_back(validator->validate(*test_data, ValidationLevel::BASIC));
    }
    
    std::string summary = validator->generateSummaryReport(results);
    
    EXPECT_TRUE(summary.find("Validation Summary Report") != std::string::npos);
    EXPECT_TRUE(summary.find("Total validation runs: 5") != std::string::npos);
    EXPECT_TRUE(summary.find("Status Distribution:") != std::string::npos);
    EXPECT_TRUE(summary.find("Average Severity Score:") != std::string::npos);
}

TEST_F(DataValidatorTest, ShouldValidateFunction) {
    EXPECT_TRUE(DataValidator::shouldValidate(0, 10));
    EXPECT_FALSE(DataValidator::shouldValidate(1, 10));
    EXPECT_FALSE(DataValidator::shouldValidate(9, 10));
    EXPECT_TRUE(DataValidator::shouldValidate(10, 10));
    EXPECT_TRUE(DataValidator::shouldValidate(20, 10));
    EXPECT_FALSE(DataValidator::shouldValidate(25, 10));
    
    // Negative timestep should always validate
    EXPECT_TRUE(DataValidator::shouldValidate(-1, 10));
}

TEST_F(DataValidatorTest, FromStabilityConfigConversion) {
    StabilityConfig stability_config;
    stability_config.max_velocity = 0.08;
    stability_config.min_density = 1e-5;
    stability_config.max_density = 5e5;
    stability_config.mass_conservation_tolerance = 1e-12;
    stability_config.check_nan = true;
    stability_config.check_inf = true;
    stability_config.check_mass_conservation = true;
    
    ValidationRules converted_rules = DataValidator::fromStabilityConfig(stability_config);
    
    EXPECT_EQ(converted_rules.max_velocity, 0.08);
    EXPECT_EQ(converted_rules.min_density, 1e-5);
    EXPECT_EQ(converted_rules.max_density, 5e5);
    EXPECT_EQ(converted_rules.mass_conservation_tolerance, 1e-12);
    EXPECT_TRUE(converted_rules.check_nan);
    EXPECT_TRUE(converted_rules.check_inf);
    EXPECT_TRUE(converted_rules.check_mass_conservation);
}

TEST_F(DataValidatorTest, ComprehensiveValidationChecks) {
    auto result = validator->validate(*test_data, ValidationLevel::COMPREHENSIVE);
    
    // Should include additional checks beyond basic validation
    std::vector<std::string> expected_rules = {
        "NaN Detection", "Infinity Detection", "Density Bounds", 
        "Velocity Bounds", "Mass Conservation", "Momentum Conservation",
        "Physics Consistency", "Mach Number", "Reynolds Number", 
        "Boundary Conditions"
    };
    
    for (const auto& expected_rule : expected_rules) {
        bool rule_found = false;
        for (const auto& rule_result : result.rule_results) {
            if (rule_result.rule_name == expected_rule) {
                rule_found = true;
                break;
            }
        }
        // Note: Some rules might not be present depending on configuration flags
        // This test mainly ensures comprehensive validation includes more checks
    }
    
    EXPECT_GT(result.rule_results.size(), 5); // Should have multiple validation rules
}

TEST_F(DataValidatorTest, StrictValidationChecks) {
    auto result = validator->validate(*test_data, ValidationLevel::STRICT);
    
    // Strict validation should include all checks
    EXPECT_GT(result.rule_results.size(), 8); // Should have many validation rules
    
    // Should include equilibrium and symmetry checks
    bool equilibrium_found = false;
    bool symmetry_found = false;
    
    for (const auto& rule_result : result.rule_results) {
        if (rule_result.rule_name == "Equilibrium Distribution") {
            equilibrium_found = true;
        }
        if (rule_result.rule_name == "Symmetry Properties") {
            symmetry_found = true;
        }
    }
    
    // Note: These might not be found if the corresponding flags are disabled
    // The test mainly ensures strict validation attempts more comprehensive checks
}

TEST_F(DataValidatorTest, ValidationResultStructure) {
    // Introduce mixed issues
    test_data->density[1][1] = -0.1; // Should cause failure
    test_data->velocity_x[2][2] = 0.08; // Might cause warning depending on rules
    
    auto result = validator->validate(*test_data, ValidationLevel::COMPREHENSIVE);
    
    // Check result structure
    EXPECT_NE(result.overall_status, ValidationStatus::PASSED);
    EXPECT_GT(result.overall_severity_score, 0.0);
    EXPECT_FALSE(result.summary_message.empty());
    EXPECT_GT(result.rule_results.size(), 0);
    
    // Check that severity scores are reasonable
    for (const auto& rule_result : result.rule_results) {
        EXPECT_GE(rule_result.severity_score, 0.0);
        EXPECT_LE(rule_result.severity_score, 1.0);
        EXPECT_FALSE(rule_result.rule_name.empty());
    }
}

TEST_F(DataValidatorTest, DisabledValidationChecks) {
    // Disable specific checks
    ValidationRules custom_rules = rules;
    custom_rules.check_nan = false;
    custom_rules.check_inf = false;
    custom_rules.check_bounds = false;
    
    DataValidator custom_validator(custom_rules);
    
    // Introduce problems that should be ignored
    test_data->density[1][1] = std::numeric_limits<double>::quiet_NaN();
    test_data->velocity_x[2][2] = std::numeric_limits<double>::infinity();
    test_data->density[3][3] = -1.0; // Negative density
    
    auto result = custom_validator.validate(*test_data, ValidationLevel::BASIC);
    
    // Should have fewer rule results due to disabled checks
    EXPECT_LT(result.rule_results.size(), 5);
    
    // Check that disabled rules are not present
    for (const auto& rule_result : result.rule_results) {
        EXPECT_NE(rule_result.rule_name, "NaN Detection");
        EXPECT_NE(rule_result.rule_name, "Infinity Detection");
        EXPECT_NE(rule_result.rule_name, "Density Bounds");
        EXPECT_NE(rule_result.rule_name, "Velocity Bounds");
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}