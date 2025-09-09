#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "../../include/lbm/single_phase.hpp"
#include "../../include/analysis/h_theorem.hpp"

using namespace lbm;

class HTheoremVerificationTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Set up configuration for H-theorem testing
        config.Nx = 20;
        config.Ny = 21;
        config.tau = 1.0;
        config.rho0 = 1.0;
        config.gravity = 1e-6;
        config.max_timesteps = 3000;
        config.output_interval = 500;
        config.use_entropic_bgk = true;  // Start with entropic for H-theorem
        config.max_velocity_limit = 0.1;
        config.min_density_limit = 1e-6;
        config.stability_check_interval = 100;
        config.output_prefix = "h_theorem_test";
        config.write_analytical_comparison = true;
        config.enable_h_theorem_analysis = true;
        config.h_analysis_interval = 50;
        
        test_output_dir = "h_theorem_verification_test";
        
        // Clean up any existing test files
        if (std::filesystem::exists(test_output_dir)) {
            std::filesystem::remove_all(test_output_dir);
        }
        std::filesystem::create_directories(test_output_dir);
    }
    
    void TearDown() override {
        // Clean up test files
        if (std::filesystem::exists(test_output_dir)) {
            std::filesystem::remove_all(test_output_dir);
        }
    }
    
    struct HEvolutionData {
        std::vector<int> timesteps;
        std::vector<double> h_values;
        std::vector<double> h_changes;
        std::vector<double> entropy_production;
        std::vector<bool> monotonic_violations;
    };
    
    HEvolutionData readHEvolutionData(const std::string& filename) {
        HEvolutionData data;
        std::ifstream file(filename);
        std::string line;
        
        // Skip header
        std::getline(file, line);
        
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string token;
            std::vector<std::string> tokens;
            
            while (std::getline(iss, token, ',')) {
                tokens.push_back(token);
            }
            
            if (tokens.size() >= 5) {
                data.timesteps.push_back(std::stoi(tokens[0]));
                data.h_values.push_back(std::stod(tokens[1]));
                data.h_changes.push_back(std::stod(tokens[2]));
                data.entropy_production.push_back(std::stod(tokens[3]));
                data.monotonic_violations.push_back(std::stoi(tokens[4]) != 0);
            }
        }
        
        return data;
    }
    
    int countMonotonicityViolations(const HEvolutionData& data) {
        return std::count(data.monotonic_violations.begin(), data.monotonic_violations.end(), true);
    }
    
    double calculateHDecayRate(const HEvolutionData& data) {
        if (data.h_values.size() < 10) return 0.0;
        
        // Fit exponential decay: H(t) = H0 * exp(-r*t)
        // Using log-linear regression: log(H) = log(H0) - r*t
        
        std::vector<double> log_h;
        std::vector<double> timesteps_double;
        
        for (size_t i = 0; i < data.h_values.size(); i++) {
            if (data.h_values[i] > 0) {
                log_h.push_back(std::log(data.h_values[i]));
                timesteps_double.push_back(static_cast<double>(data.timesteps[i]));
            }
        }
        
        if (log_h.size() < 5) return 0.0;
        
        // Simple linear regression
        double n = log_h.size();
        double sum_t = 0.0, sum_log_h = 0.0, sum_t_log_h = 0.0, sum_t_sqr = 0.0;
        
        for (size_t i = 0; i < log_h.size(); i++) {
            sum_t += timesteps_double[i];
            sum_log_h += log_h[i];
            sum_t_log_h += timesteps_double[i] * log_h[i];
            sum_t_sqr += timesteps_double[i] * timesteps_double[i];
        }
        
        double slope = (n * sum_t_log_h - sum_t * sum_log_h) / (n * sum_t_sqr - sum_t * sum_t);
        return -slope;  // Decay rate is negative of slope
    }
    
    double calculateAverageEntropyProduction(const HEvolutionData& data) {
        if (data.entropy_production.empty()) return 0.0;
        
        double sum = 0.0;
        int count = 0;
        for (double ep : data.entropy_production) {
            if (std::isfinite(ep)) {
                sum += ep;
                count++;
            }
        }
        
        return (count > 0) ? sum / count : 0.0;
    }
    
    SinglePhaseSolver::SimulationConfig config;
    std::string test_output_dir;
};

TEST_F(HTheoremVerificationTest, EntropicBGKMonotonicity) {
    // Test that entropic BGK satisfies H-theorem monotonicity
    config.use_entropic_bgk = true;
    SinglePhaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    
    EXPECT_TRUE(solver.runSimulation());
    
    // Check that H-evolution file was created
    std::string h_evolution_file = test_output_dir + "/" + config.output_prefix + "_h_evolution.csv";
    EXPECT_TRUE(std::filesystem::exists(h_evolution_file)) 
        << "H-evolution file not created";
    
    // Read H-evolution data
    auto h_data = readHEvolutionData(h_evolution_file);
    EXPECT_GT(h_data.timesteps.size(), 0) << "No H-theorem data found";
    
    // Check monotonicity: H should not increase (allow small numerical violations)
    int violations = countMonotonicityViolations(h_data);
    double violation_rate = static_cast<double>(violations) / h_data.timesteps.size();
    
    EXPECT_LT(violation_rate, 0.01) 
        << "Too many monotonicity violations for entropic BGK: " 
        << violations << " out of " << h_data.timesteps.size();
    
    // Check that H-function decreases overall
    if (h_data.h_values.size() >= 2) {
        double initial_h = h_data.h_values[0];
        double final_h = h_data.h_values.back();
        EXPECT_LT(final_h, initial_h) 
            << "H-function should decrease overall";
    }
}

TEST_F(HTheoremVerificationTest, StandardBGKMonotonicityComparison) {
    // Compare H-theorem behavior between standard and entropic BGK
    HEvolutionData standard_data, entropic_data;
    
    // Run standard BGK
    config.use_entropic_bgk = false;
    SinglePhaseSolver standard_solver(config);
    standard_solver.setOutputDirectory(test_output_dir + "/standard");
    std::filesystem::create_directories(test_output_dir + "/standard");
    
    EXPECT_TRUE(standard_solver.runSimulation());
    
    std::string standard_h_file = test_output_dir + "/standard/" + config.output_prefix + "_h_evolution.csv";
    EXPECT_TRUE(std::filesystem::exists(standard_h_file));
    standard_data = readHEvolutionData(standard_h_file);
    
    // Run entropic BGK
    config.use_entropic_bgk = true;
    SinglePhaseSolver entropic_solver(config);
    entropic_solver.setOutputDirectory(test_output_dir + "/entropic");
    std::filesystem::create_directories(test_output_dir + "/entropic");
    
    EXPECT_TRUE(entropic_solver.runSimulation());
    
    std::string entropic_h_file = test_output_dir + "/entropic/" + config.output_prefix + "_h_evolution.csv";
    EXPECT_TRUE(std::filesystem::exists(entropic_h_file));
    entropic_data = readHEvolutionData(entropic_h_file);
    
    // Compare violation rates
    int standard_violations = countMonotonicityViolations(standard_data);
    int entropic_violations = countMonotonicityViolations(entropic_data);
    
    double standard_rate = static_cast<double>(standard_violations) / standard_data.timesteps.size();
    double entropic_rate = static_cast<double>(entropic_violations) / entropic_data.timesteps.size();
    
    // Entropic BGK should have fewer or equal violations
    EXPECT_LE(entropic_rate, standard_rate + 0.005) 
        << "Entropic BGK should not have more violations than standard BGK";
    
    // Both should show overall H-function decrease
    if (standard_data.h_values.size() >= 2) {
        EXPECT_LT(standard_data.h_values.back(), standard_data.h_values[0])
            << "Standard BGK H-function should decrease overall";
    }
    
    if (entropic_data.h_values.size() >= 2) {
        EXPECT_LT(entropic_data.h_values.back(), entropic_data.h_values[0])
            << "Entropic BGK H-function should decrease overall";
    }
}

TEST_F(HTheoremVerificationTest, EntropyProductionPositivity) {
    // Test that entropy production is non-negative (H-theorem requirement)
    config.use_entropic_bgk = true;
    SinglePhaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    
    EXPECT_TRUE(solver.runSimulation());
    
    std::string h_evolution_file = test_output_dir + "/" + config.output_prefix + "_h_evolution.csv";
    auto h_data = readHEvolutionData(h_evolution_file);
    
    // Check that entropy production is mostly non-negative
    int negative_production_count = 0;
    for (double ep : h_data.entropy_production) {
        if (std::isfinite(ep) && ep < -1e-12) {  // Allow small numerical errors
            negative_production_count++;
        }
    }
    
    double negative_rate = static_cast<double>(negative_production_count) / h_data.entropy_production.size();
    EXPECT_LT(negative_rate, 0.05) 
        << "Too many negative entropy production values: " << negative_production_count;
    
    // Average entropy production should be positive
    double avg_entropy_production = calculateAverageEntropyProduction(h_data);
    EXPECT_GT(avg_entropy_production, 0.0) 
        << "Average entropy production should be positive";
}

TEST_F(HTheoremVerificationTest, HFunctionDecayRate) {
    // Test exponential decay characteristics of H-function
    config.use_entropic_bgk = true;
    config.max_timesteps = 5000;  // Longer simulation for decay analysis
    
    SinglePhaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    
    EXPECT_TRUE(solver.runSimulation());
    
    std::string h_evolution_file = test_output_dir + "/" + config.output_prefix + "_h_evolution.csv";
    auto h_data = readHEvolutionData(h_evolution_file);
    
    // Calculate decay rate
    double decay_rate = calculateHDecayRate(h_data);
    EXPECT_GT(decay_rate, 0.0) << "H-function should exhibit exponential decay";
    
    // Decay rate should be reasonable (not too fast or too slow)
    EXPECT_LT(decay_rate, 1.0) << "H-function decay rate too high";
    EXPECT_GT(decay_rate, 1e-6) << "H-function decay rate too low";
    
    // Check that later times show smaller H-function changes
    if (h_data.h_changes.size() > 100) {
        double early_avg_change = 0.0, late_avg_change = 0.0;
        int early_count = 0, late_count = 0;
        
        for (size_t i = 10; i < h_data.h_changes.size() / 3; i++) {
            early_avg_change += std::abs(h_data.h_changes[i]);
            early_count++;
        }
        
        for (size_t i = 2 * h_data.h_changes.size() / 3; i < h_data.h_changes.size(); i++) {
            late_avg_change += std::abs(h_data.h_changes[i]);
            late_count++;
        }
        
        early_avg_change /= early_count;
        late_avg_change /= late_count;
        
        EXPECT_LT(late_avg_change, early_avg_change) 
            << "H-function changes should decrease over time";
    }
}

TEST_F(HTheoremVerificationTest, ParameterDependentBehavior) {
    // Test H-theorem behavior under different parameters
    std::vector<double> tau_values = {0.6, 1.0, 1.5, 2.0};
    std::vector<double> decay_rates;
    std::vector<double> violation_rates;
    
    for (double tau : tau_values) {
        config.tau = tau;
        config.use_entropic_bgk = true;
        config.max_timesteps = 2000;
        
        std::string subdir = test_output_dir + "/tau_" + std::to_string(int(tau * 10));
        std::filesystem::create_directories(subdir);
        
        SinglePhaseSolver solver(config);
        solver.setOutputDirectory(subdir);
        
        EXPECT_TRUE(solver.runSimulation()) << "Simulation failed for tau=" << tau;
        
        std::string h_file = subdir + "/" + config.output_prefix + "_h_evolution.csv";
        auto h_data = readHEvolutionData(h_file);
        
        double decay_rate = calculateHDecayRate(h_data);
        int violations = countMonotonicityViolations(h_data);
        double violation_rate = static_cast<double>(violations) / h_data.timesteps.size();
        
        decay_rates.push_back(decay_rate);
        violation_rates.push_back(violation_rate);
        
        // Each parameter set should satisfy basic H-theorem requirements
        EXPECT_GT(decay_rate, 0.0) << "No decay detected for tau=" << tau;
        EXPECT_LT(violation_rate, 0.02) << "Too many violations for tau=" << tau;
    }
    
    // All decay rates should be positive and reasonable
    for (size_t i = 0; i < decay_rates.size(); i++) {
        EXPECT_GT(decay_rates[i], 1e-7) << "Decay rate too small for tau=" << tau_values[i];
        EXPECT_LT(decay_rates[i], 1.0) << "Decay rate too large for tau=" << tau_values[i];
    }
}

TEST_F(HTheoremVerificationTest, NumericalPrecisionEffects) {
    // Test H-theorem under different numerical precision conditions
    config.use_entropic_bgk = true;
    
    // Test with different tolerances/intervals
    std::vector<int> analysis_intervals = {10, 50, 100};
    
    for (int interval : analysis_intervals) {
        config.h_analysis_interval = interval;
        config.max_timesteps = 2000;
        
        std::string subdir = test_output_dir + "/interval_" + std::to_string(interval);
        std::filesystem::create_directories(subdir);
        
        SinglePhaseSolver solver(config);
        solver.setOutputDirectory(subdir);
        
        EXPECT_TRUE(solver.runSimulation()) 
            << "Simulation failed for interval=" << interval;
        
        std::string h_file = subdir + "/" + config.output_prefix + "_h_evolution.csv";
        auto h_data = readHEvolutionData(h_file);
        
        // Should still satisfy H-theorem regardless of analysis frequency
        int violations = countMonotonicityViolations(h_data);
        double violation_rate = static_cast<double>(violations) / h_data.timesteps.size();
        
        EXPECT_LT(violation_rate, 0.02) 
            << "H-theorem violations for analysis interval=" << interval;
        
        double decay_rate = calculateHDecayRate(h_data);
        EXPECT_GT(decay_rate, 0.0) 
            << "No decay detected for interval=" << interval;
    }
}

TEST_F(HTheoremVerificationTest, EquilibriumApproach) {
    // Test that H-function approaches minimum as system reaches equilibrium
    config.use_entropic_bgk = true;
    config.max_timesteps = 8000;  // Long simulation to approach equilibrium
    config.h_analysis_interval = 100;
    
    SinglePhaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    
    EXPECT_TRUE(solver.runSimulation());
    
    std::string h_evolution_file = test_output_dir + "/" + config.output_prefix + "_h_evolution.csv";
    auto h_data = readHEvolutionData(h_evolution_file);
    
    // Check that H-function levels off at later times
    if (h_data.h_values.size() > 50) {
        // Compare last quarter with third quarter
        size_t quarter_size = h_data.h_values.size() / 4;
        size_t third_quarter_start = 2 * quarter_size;
        size_t fourth_quarter_start = 3 * quarter_size;
        
        double third_quarter_avg = 0.0, fourth_quarter_avg = 0.0;
        
        for (size_t i = third_quarter_start; i < fourth_quarter_start; i++) {
            third_quarter_avg += h_data.h_values[i];
        }
        third_quarter_avg /= quarter_size;
        
        for (size_t i = fourth_quarter_start; i < h_data.h_values.size(); i++) {
            fourth_quarter_avg += h_data.h_values[i];
        }
        fourth_quarter_avg /= (h_data.h_values.size() - fourth_quarter_start);
        
        // Should approach a steady value
        double relative_change = std::abs(fourth_quarter_avg - third_quarter_avg) / third_quarter_avg;
        EXPECT_LT(relative_change, 0.1) 
            << "H-function not approaching equilibrium, relative change: " << relative_change;
    }
    
    // Check that entropy production approaches zero
    if (h_data.entropy_production.size() > 20) {
        double final_avg_production = 0.0;
        int count = 0;
        
        for (size_t i = h_data.entropy_production.size() - 20; i < h_data.entropy_production.size(); i++) {
            if (std::isfinite(h_data.entropy_production[i])) {
                final_avg_production += h_data.entropy_production[i];
                count++;
            }
        }
        
        if (count > 0) {
            final_avg_production /= count;
            
            // Should approach zero as system reaches equilibrium (but might not be exactly zero due to boundaries)
            EXPECT_LT(final_avg_production, calculateAverageEntropyProduction(h_data) * 0.5)
                << "Entropy production not decreasing toward equilibrium";
        }
    }
}

TEST_F(HTheoremVerificationTest, ConsistencyWithAnalyticalSolution) {
    // Test H-theorem consistency with analytical solution approach
    config.use_entropic_bgk = true;
    config.max_timesteps = 4000;
    
    SinglePhaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    
    EXPECT_TRUE(solver.runSimulation());
    
    // Check velocity convergence
    std::string final_velocity_file = test_output_dir + "/" + config.output_prefix + "_velocity_t4000.csv";
    EXPECT_TRUE(std::filesystem::exists(final_velocity_file));
    
    // Read and analyze velocity convergence
    std::ifstream file(final_velocity_file);
    std::string line;
    std::getline(file, line);  // Skip header
    
    double max_relative_error = 0.0;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        std::vector<std::string> tokens;
        
        while (std::getline(iss, token, ',')) {
            tokens.push_back(token);
        }
        
        if (tokens.size() >= 3) {
            double numerical = std::stod(tokens[1]);
            double analytical = std::stod(tokens[2]);
            
            if (std::abs(analytical) > 1e-10) {
                double relative_error = std::abs(numerical - analytical) / std::abs(analytical);
                max_relative_error = std::max(max_relative_error, relative_error);
            }
        }
    }
    
    // Good convergence to analytical solution should correlate with proper H-theorem behavior
    EXPECT_LT(max_relative_error, 0.05) << "Poor convergence to analytical solution";
    
    // Check H-theorem data
    std::string h_evolution_file = test_output_dir + "/" + config.output_prefix + "_h_evolution.csv";
    auto h_data = readHEvolutionData(h_evolution_file);
    
    int violations = countMonotonicityViolations(h_data);
    double violation_rate = static_cast<double>(violations) / h_data.timesteps.size();
    
    // Good analytical convergence should correspond to good H-theorem behavior
    EXPECT_LT(violation_rate, 0.01) 
        << "H-theorem violations despite good analytical convergence";
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}