#include <gtest/gtest.h>
#include "lbm/stability_monitor.hpp"
#include <limits>
#include <cmath>

using namespace lbm;

class StabilityMonitorTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create test data with reasonable dimensions
        test_data = std::make_unique<LatticeData>(10, 5);
        
        // Initialize with stable values
        for (int i = 0; i < test_data->nx; ++i) {
            for (int j = 0; j < test_data->ny; ++j) {
                test_data->density[i][j] = 1.0;
                test_data->velocity_x[i][j] = 0.01;
                test_data->velocity_y[i][j] = 0.01;
                
                // Initialize distribution functions
                for (int k = 0; k < 9; ++k) {
                    test_data->distribution[i][j][k] = 0.1;
                }
            }
        }
        
        // Default configuration
        config = StabilityConfig{};
        monitor = std::make_unique<StabilityMonitor>(config);
    }
    
    std::unique_ptr<LatticeData> test_data;
    StabilityConfig config;
    std::unique_ptr<StabilityMonitor> monitor;
};

TEST_F(StabilityMonitorTest, StableDataReturnsStableMetrics) {
    auto metrics = monitor->checkStability(*test_data);
    
    EXPECT_TRUE(metrics.isStable());
    EXPECT_TRUE(metrics.density_positive);
    EXPECT_TRUE(metrics.velocity_bounded);
    EXPECT_TRUE(metrics.no_nan_values);
    EXPECT_TRUE(metrics.no_inf_values);
    EXPECT_TRUE(metrics.mass_conserved);
    EXPECT_EQ(metrics.nan_count, 0);
    EXPECT_EQ(metrics.inf_count, 0);
    EXPECT_EQ(metrics.negative_density_count, 0);
    EXPECT_EQ(metrics.excessive_velocity_count, 0);
}

TEST_F(StabilityMonitorTest, DetectsNaNValues) {
    // Introduce NaN values
    test_data->density[5][2] = std::numeric_limits<double>::quiet_NaN();
    test_data->velocity_x[3][1] = std::numeric_limits<double>::quiet_NaN();
    
    auto metrics = monitor->checkStability(*test_data);
    
    EXPECT_FALSE(metrics.isStable());
    EXPECT_FALSE(metrics.no_nan_values);
    EXPECT_EQ(metrics.nan_count, 2);
    EXPECT_TRUE(metrics.error_message.find("NaN values detected") != std::string::npos);
}

TEST_F(StabilityMonitorTest, DetectsInfinityValues) {
    // Introduce infinity values
    test_data->density[2][3] = std::numeric_limits<double>::infinity();
    test_data->velocity_y[4][1] = -std::numeric_limits<double>::infinity();
    test_data->distribution[1][2][5] = std::numeric_limits<double>::infinity();
    
    auto metrics = monitor->checkStability(*test_data);
    
    EXPECT_FALSE(metrics.isStable());
    EXPECT_FALSE(metrics.no_inf_values);
    EXPECT_EQ(metrics.inf_count, 3);
    EXPECT_TRUE(metrics.error_message.find("Infinity values detected") != std::string::npos);
}

TEST_F(StabilityMonitorTest, DetectsNegativeDensity) {
    // Introduce negative densities
    test_data->density[1][1] = -0.5;
    test_data->density[2][2] = -0.1;
    test_data->density[3][3] = -1e-8;
    
    auto metrics = monitor->checkStability(*test_data);
    
    EXPECT_FALSE(metrics.isStable());
    EXPECT_FALSE(metrics.density_positive);
    EXPECT_EQ(metrics.negative_density_count, 3);
    EXPECT_TRUE(metrics.error_message.find("Negative densities detected") != std::string::npos);
}

TEST_F(StabilityMonitorTest, DetectsDensityOutOfBounds) {
    // Test minimum density violation
    test_data->density[1][1] = 1e-8; // Below default min_density (1e-6)
    
    auto metrics = monitor->checkStability(*test_data);
    
    EXPECT_FALSE(metrics.isStable());
    EXPECT_FALSE(metrics.density_positive);
    EXPECT_TRUE(metrics.error_message.find("Density below minimum") != std::string::npos);
    
    // Reset and test maximum density violation
    test_data->density[1][1] = 1.0;
    test_data->density[2][2] = 2e6; // Above default max_density (1e6)
    
    metrics = monitor->checkStability(*test_data);
    
    EXPECT_FALSE(metrics.isStable());
    EXPECT_FALSE(metrics.density_positive);
    EXPECT_TRUE(metrics.error_message.find("Density above maximum") != std::string::npos);
}

TEST_F(StabilityMonitorTest, DetectsExcessiveVelocity) {
    // Introduce excessive velocities
    test_data->velocity_x[1][1] = 0.08; // Combined magnitude will exceed 0.1
    test_data->velocity_y[1][1] = 0.08;
    test_data->velocity_x[2][2] = 0.15; // This alone exceeds 0.1
    test_data->velocity_y[2][2] = 0.0;
    
    auto metrics = monitor->checkStability(*test_data);
    
    EXPECT_FALSE(metrics.isStable());
    EXPECT_FALSE(metrics.velocity_bounded);
    EXPECT_EQ(metrics.excessive_velocity_count, 2);
    EXPECT_GT(metrics.max_velocity, 0.1);
    EXPECT_TRUE(metrics.error_message.find("Excessive velocities detected") != std::string::npos);
}

TEST_F(StabilityMonitorTest, DetectsMassConservationViolation) {
    // First check establishes initial mass
    auto metrics1 = monitor->checkStability(*test_data);
    EXPECT_TRUE(metrics1.mass_conserved);
    EXPECT_EQ(metrics1.mass_change, 0.0);
    
    // Modify total mass significantly
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 3; ++j) {
            test_data->density[i][j] = 2.0; // Double the density in part of domain
        }
    }
    
    auto metrics2 = monitor->checkStability(*test_data);
    
    EXPECT_FALSE(metrics2.isStable());
    EXPECT_FALSE(metrics2.mass_conserved);
    EXPECT_GT(metrics2.mass_change, config.mass_conservation_tolerance);
    EXPECT_TRUE(metrics2.error_message.find("Mass not conserved") != std::string::npos);
}

TEST_F(StabilityMonitorTest, ConfigurableParameters) {
    // Create custom configuration
    StabilityConfig custom_config;
    custom_config.max_velocity = 0.05;
    custom_config.min_density = 0.5;
    custom_config.max_density = 2.0;
    custom_config.mass_conservation_tolerance = 1e-12;
    
    StabilityMonitor custom_monitor(custom_config);
    
    // Test with velocity that would be acceptable with default config
    test_data->velocity_x[1][1] = 0.08;
    test_data->velocity_y[1][1] = 0.0;
    
    auto metrics = custom_monitor.checkStability(*test_data);
    
    EXPECT_FALSE(metrics.isStable());
    EXPECT_FALSE(metrics.velocity_bounded);
    EXPECT_EQ(metrics.excessive_velocity_count, 1);
}

TEST_F(StabilityMonitorTest, CheckIntervalFunctionality) {
    config.check_interval = 5;
    StabilityMonitor interval_monitor(config);
    
    EXPECT_TRUE(interval_monitor.shouldCheck(0));
    EXPECT_FALSE(interval_monitor.shouldCheck(1));
    EXPECT_FALSE(interval_monitor.shouldCheck(4));
    EXPECT_TRUE(interval_monitor.shouldCheck(5));
    EXPECT_TRUE(interval_monitor.shouldCheck(10));
    EXPECT_FALSE(interval_monitor.shouldCheck(13));
}

TEST_F(StabilityMonitorTest, IssueCountingAndLogging) {
    EXPECT_FALSE(monitor->hasStabilityIssues());
    EXPECT_EQ(monitor->getIssueCount(), 0);
    
    // Introduce instability
    test_data->density[1][1] = -0.5;
    monitor->checkStability(*test_data);
    
    EXPECT_TRUE(monitor->hasStabilityIssues());
    EXPECT_EQ(monitor->getIssueCount(), 1);
    
    // Introduce another instability
    test_data->velocity_x[2][2] = 0.2;
    monitor->checkStability(*test_data);
    
    EXPECT_EQ(monitor->getIssueCount(), 2);
}

TEST_F(StabilityMonitorTest, ResetFunctionality) {
    // Introduce instability and check
    test_data->density[1][1] = -0.5;
    monitor->checkStability(*test_data);
    
    EXPECT_TRUE(monitor->hasStabilityIssues());
    EXPECT_EQ(monitor->getIssueCount(), 1);
    EXPECT_FALSE(monitor->getStabilityHistory().empty());
    
    // Reset and verify
    monitor->reset();
    
    EXPECT_FALSE(monitor->hasStabilityIssues());
    EXPECT_EQ(monitor->getIssueCount(), 0);
    EXPECT_TRUE(monitor->getStabilityHistory().empty());
}

TEST_F(StabilityMonitorTest, ConfigurationUpdate) {
    StabilityConfig new_config;
    new_config.max_velocity = 0.2;
    new_config.check_interval = 50;
    
    monitor->updateConfig(new_config);
    
    const auto& retrieved_config = monitor->getConfig();
    EXPECT_EQ(retrieved_config.max_velocity, 0.2);
    EXPECT_EQ(retrieved_config.check_interval, 50);
}

TEST_F(StabilityMonitorTest, DiagnosticReportGeneration) {
    // Introduce some instabilities
    test_data->density[1][1] = -0.5;
    test_data->velocity_x[2][2] = 0.15;
    
    monitor->checkStability(*test_data, 100);
    
    std::string report = monitor->generateDiagnosticReport();
    
    EXPECT_TRUE(report.find("Stability Monitor Diagnostic Report") != std::string::npos);
    EXPECT_TRUE(report.find("Total issues detected: 1") != std::string::npos);
    EXPECT_TRUE(report.find("Configuration:") != std::string::npos);
    EXPECT_TRUE(report.find("Latest Stability Metrics:") != std::string::npos);
    EXPECT_TRUE(report.find("Stable: NO") != std::string::npos);
}

TEST_F(StabilityMonitorTest, StabilityHistoryTracking) {
    config.real_time_monitoring = true;
    StabilityMonitor history_monitor(config);
    
    // Perform multiple checks
    for (int i = 0; i < 5; ++i) {
        history_monitor.checkStability(*test_data, i);
    }
    
    const auto& history = history_monitor.getStabilityHistory();
    EXPECT_EQ(history.size(), 5);
    
    // All should be stable
    for (const auto& metrics : history) {
        EXPECT_TRUE(metrics.isStable());
    }
}

TEST_F(StabilityMonitorTest, DisabledChecks) {
    // Disable specific checks
    config.check_nan = false;
    config.check_inf = false;
    config.check_velocity_bounds = false;
    
    StabilityMonitor selective_monitor(config);
    
    // Introduce problems that should be ignored
    test_data->density[1][1] = std::numeric_limits<double>::quiet_NaN();
    test_data->velocity_x[2][2] = std::numeric_limits<double>::infinity();
    test_data->velocity_y[3][3] = 0.2; // Excessive velocity
    
    auto metrics = selective_monitor.checkStability(*test_data);
    
    // Should still detect density issues (not disabled)
    EXPECT_FALSE(metrics.isStable()); // Due to NaN in density affecting density_positive
    EXPECT_EQ(metrics.nan_count, 0); // NaN check disabled
    EXPECT_EQ(metrics.inf_count, 0); // Inf check disabled
    EXPECT_EQ(metrics.excessive_velocity_count, 0); // Velocity check disabled
}

// Test edge cases and boundary conditions
TEST_F(StabilityMonitorTest, EdgeCaseValues) {
    // Test with very small positive values
    test_data->density[1][1] = 1e-15;
    test_data->velocity_x[2][2] = 1e-10;
    
    auto metrics = monitor->checkStability(*test_data);
    
    // Should detect density below minimum
    EXPECT_FALSE(metrics.isStable());
    EXPECT_FALSE(metrics.density_positive);
    
    // Test with values exactly at boundaries
    test_data->density[1][1] = config.min_density; // Exactly at minimum
    test_data->velocity_x[2][2] = config.max_velocity; // Exactly at maximum
    test_data->velocity_y[2][2] = 0.0;
    
    metrics = monitor->checkStability(*test_data);
    
    // Should be stable (boundary values are acceptable)
    EXPECT_TRUE(metrics.isStable());
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}