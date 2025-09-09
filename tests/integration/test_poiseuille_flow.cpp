#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <cmath>
#include "../../include/lbm/single_phase.hpp"

using namespace lbm;

class PoiseuilleFlowIntegrationTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Set up realistic Poiseuille flow parameters
        config.Nx = 20;
        config.Ny = 21;
        config.tau = 1.0;
        config.rho0 = 1.0;
        config.gravity = 1e-6;
        config.max_timesteps = 5000;
        config.output_interval = 1000;
        config.use_entropic_bgk = false;
        config.max_velocity_limit = 0.1;
        config.min_density_limit = 1e-6;
        config.stability_check_interval = 100;
        config.output_prefix = "integration_test";
        config.write_analytical_comparison = true;
        
        test_output_dir = "integration_test_output";
        
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
    
    struct ErrorMetrics {
        double l1_error;
        double l2_error;
        double max_error;
        double relative_error;
    };
    
    ErrorMetrics calculateErrorMetrics(const std::vector<double>& numerical,
                                     const std::vector<double>& analytical) {
        ErrorMetrics metrics = {0.0, 0.0, 0.0, 0.0};
        
        if (numerical.size() != analytical.size()) {
            return metrics;
        }
        
        double sum_abs_error = 0.0;
        double sum_sqr_error = 0.0;
        double max_analytical = 0.0;
        
        for (size_t i = 0; i < numerical.size(); i++) {
            double error = std::abs(numerical[i] - analytical[i]);
            sum_abs_error += error;
            sum_sqr_error += error * error;
            metrics.max_error = std::max(metrics.max_error, error);
            max_analytical = std::max(max_analytical, std::abs(analytical[i]));
        }
        
        metrics.l1_error = sum_abs_error / numerical.size();
        metrics.l2_error = std::sqrt(sum_sqr_error / numerical.size());
        metrics.relative_error = (max_analytical > 0) ? metrics.l2_error / max_analytical : 0.0;
        
        return metrics;
    }
    
    std::vector<std::pair<double, double>> readVelocityProfile(const std::string& filename) {
        std::vector<std::pair<double, double>> profile;
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
            
            if (tokens.size() >= 3) {
                // double y = std::stod(tokens[0]);  // Not used in this function
                double ux_numerical = std::stod(tokens[1]);
                double ux_analytical = std::stod(tokens[2]);
                profile.push_back({ux_numerical, ux_analytical});
            }
        }
        
        return profile;
    }
    
    SinglePhaseSolver::SimulationConfig config;
    std::string test_output_dir;
};

TEST_F(PoiseuilleFlowIntegrationTest, StandardBGKConvergence) {
    // Test complete simulation workflow with standard BGK
    config.use_entropic_bgk = false;
    SinglePhaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    
    // Run complete simulation
    EXPECT_TRUE(solver.runSimulation());
    
    // Check that all expected output files were created
    std::vector<int> expected_timesteps = {0, 1000, 2000, 3000, 4000, 5000};
    for (int t : expected_timesteps) {
        std::string filename = test_output_dir + "/integration_test_velocity_t" + std::to_string(t) + ".csv";
        EXPECT_TRUE(std::filesystem::exists(filename)) 
            << "Missing output file: " << filename;
    }
    
    // Analyze final results
    std::string final_output = test_output_dir + "/integration_test_velocity_t5000.csv";
    auto velocity_profile = readVelocityProfile(final_output);
    EXPECT_GT(velocity_profile.size(), 0) << "Could not read velocity profile";
    
    // Extract numerical and analytical velocities
    std::vector<double> numerical, analytical;
    for (const auto& point : velocity_profile) {
        numerical.push_back(point.first);
        analytical.push_back(point.second);
    }
    
    // Calculate error metrics
    auto metrics = calculateErrorMetrics(numerical, analytical);
    
    // Standard BGK should converge with reasonable accuracy
    EXPECT_LT(metrics.relative_error, 0.05) 
        << "Standard BGK relative error too high: " << metrics.relative_error;
    EXPECT_LT(metrics.l2_error, 1e-4) 
        << "Standard BGK L2 error too high: " << metrics.l2_error;
    
    // Check velocity profile shape
    double max_numerical = *std::max_element(numerical.begin(), numerical.end());
    double max_analytical = *std::max_element(analytical.begin(), analytical.end());
    EXPECT_NEAR(max_numerical, max_analytical, max_analytical * 0.05)
        << "Maximum velocity mismatch";
    
    // Check boundary conditions (walls should have zero velocity)
    EXPECT_NEAR(numerical.front(), 0.0, 1e-6) << "Wall velocity not zero";
    EXPECT_NEAR(numerical.back(), 0.0, 1e-6) << "Wall velocity not zero";
}

TEST_F(PoiseuilleFlowIntegrationTest, EntropicBGKConvergence) {
    // Test complete simulation workflow with entropic BGK
    config.use_entropic_bgk = true;
    SinglePhaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    
    // Run complete simulation
    EXPECT_TRUE(solver.runSimulation());
    
    // Check that all expected output files were created including H-function files
    std::vector<int> expected_timesteps = {0, 1000, 2000, 3000, 4000, 5000};
    for (int t : expected_timesteps) {
        std::string velocity_file = test_output_dir + "/integration_test_velocity_t" + std::to_string(t) + ".csv";
        std::string h_file = test_output_dir + "/integration_test_h_y_t" + std::to_string(t) + ".csv";
        
        EXPECT_TRUE(std::filesystem::exists(velocity_file)) 
            << "Missing velocity file: " << velocity_file;
        EXPECT_TRUE(std::filesystem::exists(h_file)) 
            << "Missing H-function file: " << h_file;
    }
    
    // Analyze final results
    std::string final_output = test_output_dir + "/integration_test_velocity_t5000.csv";
    auto velocity_profile = readVelocityProfile(final_output);
    EXPECT_GT(velocity_profile.size(), 0) << "Could not read velocity profile";
    
    // Extract numerical and analytical velocities
    std::vector<double> numerical, analytical;
    for (const auto& point : velocity_profile) {
        numerical.push_back(point.first);
        analytical.push_back(point.second);
    }
    
    // Calculate error metrics
    auto metrics = calculateErrorMetrics(numerical, analytical);
    
    // Entropic BGK should also converge with good accuracy
    EXPECT_LT(metrics.relative_error, 0.05) 
        << "Entropic BGK relative error too high: " << metrics.relative_error;
    EXPECT_LT(metrics.l2_error, 1e-4) 
        << "Entropic BGK L2 error too high: " << metrics.l2_error;
}

TEST_F(PoiseuilleFlowIntegrationTest, ConvergenceComparison) {
    // Compare convergence of standard vs entropic BGK
    ErrorMetrics standard_metrics, entropic_metrics;
    
    // Run standard BGK
    config.use_entropic_bgk = false;
    SinglePhaseSolver standard_solver(config);
    standard_solver.setOutputDirectory(test_output_dir + "/standard");
    std::filesystem::create_directories(test_output_dir + "/standard");
    
    EXPECT_TRUE(standard_solver.runSimulation());
    
    // Analyze standard BGK results
    std::string standard_output = test_output_dir + "/standard/integration_test_velocity_t5000.csv";
    auto standard_profile = readVelocityProfile(standard_output);
    
    std::vector<double> standard_numerical, standard_analytical;
    for (const auto& point : standard_profile) {
        standard_numerical.push_back(point.first);
        standard_analytical.push_back(point.second);
    }
    standard_metrics = calculateErrorMetrics(standard_numerical, standard_analytical);
    
    // Run entropic BGK
    config.use_entropic_bgk = true;
    SinglePhaseSolver entropic_solver(config);
    entropic_solver.setOutputDirectory(test_output_dir + "/entropic");
    std::filesystem::create_directories(test_output_dir + "/entropic");
    
    EXPECT_TRUE(entropic_solver.runSimulation());
    
    // Analyze entropic BGK results
    std::string entropic_output = test_output_dir + "/entropic/integration_test_velocity_t5000.csv";
    auto entropic_profile = readVelocityProfile(entropic_output);
    
    std::vector<double> entropic_numerical, entropic_analytical;
    for (const auto& point : entropic_profile) {
        entropic_numerical.push_back(point.first);
        entropic_analytical.push_back(point.second);
    }
    entropic_metrics = calculateErrorMetrics(entropic_numerical, entropic_analytical);
    
    // Both methods should achieve similar accuracy
    EXPECT_LT(std::abs(standard_metrics.relative_error - entropic_metrics.relative_error), 0.02)
        << "Large difference in convergence accuracy between BGK methods";
    
    // Both should meet minimum accuracy requirements
    EXPECT_LT(standard_metrics.relative_error, 0.1) << "Standard BGK accuracy insufficient";
    EXPECT_LT(entropic_metrics.relative_error, 0.1) << "Entropic BGK accuracy insufficient";
}

TEST_F(PoiseuilleFlowIntegrationTest, GridResolutionScaling) {
    // Test how accuracy scales with grid resolution
    std::vector<int> grid_sizes = {11, 21, 41};
    std::vector<ErrorMetrics> metrics_by_resolution;
    
    for (int Ny : grid_sizes) {
        config.Ny = Ny;
        config.use_entropic_bgk = false;
        config.max_timesteps = 3000;  // Shorter for multiple runs
        
        std::string subdir = test_output_dir + "/grid_" + std::to_string(Ny);
        std::filesystem::create_directories(subdir);
        
        SinglePhaseSolver solver(config);
        solver.setOutputDirectory(subdir);
        
        EXPECT_TRUE(solver.runSimulation()) << "Simulation failed for Ny=" << Ny;
        
        // Analyze results
        std::string output_file = subdir + "/integration_test_velocity_t3000.csv";
        auto profile = readVelocityProfile(output_file);
        
        std::vector<double> numerical, analytical;
        for (const auto& point : profile) {
            numerical.push_back(point.first);
            analytical.push_back(point.second);
        }
        
        auto metrics = calculateErrorMetrics(numerical, analytical);
        metrics_by_resolution.push_back(metrics);
    }
    
    // Error should generally decrease with increasing resolution
    EXPECT_GT(metrics_by_resolution[0].l2_error, metrics_by_resolution[2].l2_error * 0.5)
        << "No clear improvement with increased resolution";
    
    // All resolutions should converge
    for (size_t i = 0; i < metrics_by_resolution.size(); i++) {
        EXPECT_LT(metrics_by_resolution[i].relative_error, 0.15)
            << "Poor convergence for grid size " << grid_sizes[i];
    }
}

TEST_F(PoiseuilleFlowIntegrationTest, ViscosityAccuracy) {
    // Test accuracy for different viscosity values (tau values)
    std::vector<double> tau_values = {0.6, 1.0, 1.5, 2.0};
    
    for (double tau : tau_values) {
        config.tau = tau;
        config.use_entropic_bgk = false;
        config.max_timesteps = 4000;
        
        std::string subdir = test_output_dir + "/tau_" + std::to_string(int(tau*10));
        std::filesystem::create_directories(subdir);
        
        SinglePhaseSolver solver(config);
        solver.setOutputDirectory(subdir);
        
        EXPECT_TRUE(solver.runSimulation()) << "Simulation failed for tau=" << tau;
        
        // Analyze results
        std::string output_file = subdir + "/integration_test_velocity_t4000.csv";
        auto profile = readVelocityProfile(output_file);
        
        std::vector<double> numerical, analytical;
        for (const auto& point : profile) {
            numerical.push_back(point.first);
            analytical.push_back(point.second);
        }
        
        auto metrics = calculateErrorMetrics(numerical, analytical);
        
        // Should maintain good accuracy across viscosity range
        EXPECT_LT(metrics.relative_error, 0.08) 
            << "Poor accuracy for tau=" << tau << ", error=" << metrics.relative_error;
        
        // Check that maximum velocity scales correctly with viscosity
        double max_numerical = *std::max_element(numerical.begin(), numerical.end());
        double max_analytical = *std::max_element(analytical.begin(), analytical.end());
        
        EXPECT_NEAR(max_numerical, max_analytical, max_analytical * 0.05)
            << "Viscosity scaling incorrect for tau=" << tau;
    }
}

TEST_F(PoiseuilleFlowIntegrationTest, StabilityUnderParameterVariation) {
    // Test stability with various parameter combinations
    struct TestParams {
        double tau;
        double gravity;
        bool entropic;
        std::string name;
    };
    
    std::vector<TestParams> test_cases = {
        {0.55, 1e-6, false, "low_viscosity_standard"},
        {2.5, 1e-6, false, "high_viscosity_standard"},
        {0.55, 1e-5, false, "high_force_standard"},
        {1.0, 1e-7, false, "low_force_standard"},
        {0.55, 1e-6, true, "low_viscosity_entropic"},
        {2.5, 1e-6, true, "high_viscosity_entropic"},
        {0.55, 1e-5, true, "high_force_entropic"}
    };
    
    for (const auto& params : test_cases) {
        config.tau = params.tau;
        config.gravity = params.gravity;
        config.use_entropic_bgk = params.entropic;
        config.max_timesteps = 2000;
        
        std::string subdir = test_output_dir + "/" + params.name;
        std::filesystem::create_directories(subdir);
        
        SinglePhaseSolver solver(config);
        solver.setOutputDirectory(subdir);
        
        EXPECT_TRUE(solver.runSimulation()) 
            << "Simulation failed for " << params.name;
        
        // Check final stability
        EXPECT_TRUE(solver.isStable()) 
            << "Solver unstable at end for " << params.name;
        
        // Check that reasonable results were produced
        std::string output_file = subdir + "/integration_test_velocity_t2000.csv";
        EXPECT_TRUE(std::filesystem::exists(output_file))
            << "No output file for " << params.name;
        
        auto profile = readVelocityProfile(output_file);
        EXPECT_GT(profile.size(), 0) << "Empty profile for " << params.name;
        
        // Check that velocities are reasonable
        std::vector<double> numerical;
        for (const auto& point : profile) {
            numerical.push_back(point.first);
        }
        
        double max_vel = *std::max_element(numerical.begin(), numerical.end());
        EXPECT_GT(max_vel, 0.0) << "No flow detected for " << params.name;
        EXPECT_LT(max_vel, 0.3) << "Unrealistic velocities for " << params.name;
    }
}

TEST_F(PoiseuilleFlowIntegrationTest, LongTermStability) {
    // Test long-term stability and steady-state behavior
    config.max_timesteps = 10000;
    config.output_interval = 2000;
    config.use_entropic_bgk = false;
    
    SinglePhaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    
    EXPECT_TRUE(solver.runSimulation());
    
    // Check that solver remained stable throughout
    EXPECT_TRUE(solver.isStable());
    
    // Analyze convergence to steady state
    std::vector<int> timesteps = {4000, 6000, 8000, 10000};
    std::vector<std::vector<double>> profiles;
    
    for (int t : timesteps) {
        std::string filename = test_output_dir + "/integration_test_velocity_t" + std::to_string(t) + ".csv";
        auto profile = readVelocityProfile(filename);
        
        std::vector<double> numerical;
        for (const auto& point : profile) {
            numerical.push_back(point.first);
        }
        profiles.push_back(numerical);
    }
    
    // Check that profiles are converging (changes between late times should be small)
    for (size_t i = 1; i < profiles.size(); i++) {
        auto metrics = calculateErrorMetrics(profiles[i], profiles[i-1]);
        EXPECT_LT(metrics.relative_error, 0.01) 
            << "Large changes in velocity profile between timesteps " 
            << timesteps[i-1] << " and " << timesteps[i];
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}