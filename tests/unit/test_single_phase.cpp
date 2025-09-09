#include <gtest/gtest.h>
#include <filesystem>
#include "../../include/lbm/single_phase.hpp"

using namespace lbm;

class SinglePhaseTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Set up basic configuration
        config.Nx = 10;
        config.Ny = 11;
        config.tau = 1.0;
        config.rho0 = 1.0;
        config.gravity = 1e-6;
        config.max_timesteps = 100;
        config.output_interval = 50;
        config.use_entropic_bgk = false;
        config.max_velocity_limit = 0.1;
        config.min_density_limit = 1e-6;
        config.stability_check_interval = 10;
        config.output_prefix = "test";
        config.write_analytical_comparison = true;
        
        test_output_dir = "test_single_phase_output";
        
        // Clean up any existing test files
        if (std::filesystem::exists(test_output_dir)) {
            std::filesystem::remove_all(test_output_dir);
        }
    }
    
    void TearDown() override {
        // Clean up test files
        if (std::filesystem::exists(test_output_dir)) {
            std::filesystem::remove_all(test_output_dir);
        }
    }
    
    SinglePhaseSolver::SimulationConfig config;
    std::string test_output_dir;
};

TEST_F(SinglePhaseTest, Construction) {
    // Test valid construction
    EXPECT_NO_THROW(SinglePhaseSolver solver(config));
    
    SinglePhaseSolver solver(config);
    EXPECT_EQ(solver.getGridDimensions().first, config.Nx);
    EXPECT_EQ(solver.getGridDimensions().second, config.Ny);
    EXPECT_EQ(solver.getConfig().tau, config.tau);
    
    // Test invalid configurations
    auto invalid_config = config;
    invalid_config.tau = 0.3;  // Too small for stability
    EXPECT_THROW(SinglePhaseSolver invalid_solver(invalid_config), std::invalid_argument);
    
    invalid_config = config;
    invalid_config.rho0 = -1.0;  // Negative density
    EXPECT_THROW(SinglePhaseSolver invalid_solver(invalid_config), std::invalid_argument);
    
    invalid_config = config;
    invalid_config.max_timesteps = 0;  // Invalid timesteps
    EXPECT_THROW(SinglePhaseSolver invalid_solver(invalid_config), std::invalid_argument);
}

TEST_F(SinglePhaseTest, Initialization) {
    SinglePhaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    
    // Test initialization
    EXPECT_NO_THROW(solver.initialize());
    
    // Check that macroscopic variables are accessible after initialization
    double density = solver.getDensity(5, 5);
    EXPECT_GT(density, 0.0);
    EXPECT_NEAR(density, config.rho0, 1e-10);
    
    auto velocity = solver.getVelocity(5, 5);
    EXPECT_TRUE(std::isfinite(velocity.first));
    EXPECT_TRUE(std::isfinite(velocity.second));
    
    // Test that solver is initially stable
    EXPECT_TRUE(solver.isStable());
}

TEST_F(SinglePhaseTest, Timestepping) {
    SinglePhaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    solver.initialize();
    
    // Test initial timestep
    EXPECT_EQ(static_cast<const SinglePhaseSolver&>(solver).getCurrentTimestep(), 0);
    
    // Test stepping
    EXPECT_NO_THROW(solver.step());
    EXPECT_EQ(static_cast<const SinglePhaseSolver&>(solver).getCurrentTimestep(), 1);
    
    EXPECT_NO_THROW(solver.step());
    EXPECT_EQ(static_cast<const SinglePhaseSolver&>(solver).getCurrentTimestep(), 2);
    
    // Test that solver remains stable after a few steps
    for (int i = 0; i < 10; i++) {
        solver.step();
    }
    EXPECT_TRUE(solver.isStable());
}

TEST_F(SinglePhaseTest, MacroscopicVariables) {
    SinglePhaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    solver.initialize();
    
    // Test density access
    for (int x = 0; x < config.Nx; x++) {
        for (int y = 0; y < config.Ny; y++) {
            double density = solver.getDensity(x, y);
            EXPECT_GT(density, 0.0);
            EXPECT_TRUE(std::isfinite(density));
        }
    }
    
    // Test velocity access
    for (int x = 0; x < config.Nx; x++) {
        for (int y = 0; y < config.Ny; y++) {
            auto velocity = solver.getVelocity(x, y);
            EXPECT_TRUE(std::isfinite(velocity.first));
            EXPECT_TRUE(std::isfinite(velocity.second));
        }
    }
    
    // Test out-of-bounds access
    EXPECT_EQ(solver.getDensity(-1, 5), 0.0);
    EXPECT_EQ(solver.getDensity(config.Nx, 5), 0.0);
    EXPECT_EQ(solver.getDensity(5, -1), 0.0);
    EXPECT_EQ(solver.getDensity(5, config.Ny), 0.0);
    
    auto invalid_velocity = solver.getVelocity(-1, 5);
    EXPECT_EQ(invalid_velocity.first, 0.0);
    EXPECT_EQ(invalid_velocity.second, 0.0);
}

TEST_F(SinglePhaseTest, AnalyticalComparison) {
    SinglePhaseSolver solver(config);
    
    // Test analytical velocity calculation
    for (int y = 0; y < config.Ny; y++) {
        double analytical = solver.getAnalyticalVelocity(y);
        EXPECT_TRUE(std::isfinite(analytical));
        EXPECT_GE(analytical, 0.0);  // Should be non-negative for Poiseuille flow
    }
    
    // Test that analytical velocity is maximum at channel center
    int center_y = config.Ny / 2;
    double center_velocity = solver.getAnalyticalVelocity(center_y);
    double wall_velocity1 = solver.getAnalyticalVelocity(0);
    double wall_velocity2 = solver.getAnalyticalVelocity(config.Ny - 1);
    
    EXPECT_GT(center_velocity, wall_velocity1);
    EXPECT_GT(center_velocity, wall_velocity2);
    EXPECT_NEAR(wall_velocity1, 0.0, 1e-10);
    EXPECT_NEAR(wall_velocity2, 0.0, 1e-10);
}

TEST_F(SinglePhaseTest, HFunctionCalculation) {
    // Test with entropic BGK
    config.use_entropic_bgk = true;
    SinglePhaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    solver.initialize();
    
    // Test H-function access
    for (int x = 0; x < config.Nx; x++) {
        for (int y = 0; y < config.Ny; y++) {
            double H = solver.getHFunction(x, y);
            EXPECT_TRUE(std::isfinite(H));
        }
    }
    
    // Test out-of-bounds access
    EXPECT_EQ(solver.getHFunction(-1, 5), 0.0);
    EXPECT_EQ(solver.getHFunction(config.Nx, 5), 0.0);
}

TEST_F(SinglePhaseTest, OutputGeneration) {
    SinglePhaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    solver.initialize();
    
    // Test output writing
    EXPECT_NO_THROW(solver.writeOutput(0));
    
    // Check that output files are created
    std::string velocity_file = test_output_dir + "/test_velocity_t0.csv";
    EXPECT_TRUE(std::filesystem::exists(velocity_file));
    
    // Test with entropic BGK (should create H-function file)
    config.use_entropic_bgk = true;
    SinglePhaseSolver entropic_solver(config);
    entropic_solver.setOutputDirectory(test_output_dir);
    entropic_solver.initialize();
    entropic_solver.writeOutput(0);
    
    std::string h_file = test_output_dir + "/test_h_y_t0.csv";
    EXPECT_TRUE(std::filesystem::exists(h_file));
}

TEST_F(SinglePhaseTest, StabilityMonitoring) {
    SinglePhaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    solver.initialize();
    
    // Test initial stability
    EXPECT_TRUE(solver.isStable());
    
    // Run a few timesteps and check stability
    for (int i = 0; i < 20; i++) {
        solver.step();
        EXPECT_TRUE(solver.isStable());  // Should remain stable for normal parameters
    }
    
    // Test stability monitor access
    const auto& monitor = solver.getStabilityMonitor();
    EXPECT_FALSE(monitor.hasStabilityIssues());
}

TEST_F(SinglePhaseTest, CollisionOperators) {
    // Test standard BGK
    config.use_entropic_bgk = false;
    SinglePhaseSolver standard_solver(config);
    standard_solver.setOutputDirectory(test_output_dir);
    standard_solver.initialize();
    
    for (int i = 0; i < 10; i++) {
        EXPECT_NO_THROW(standard_solver.step());
    }
    EXPECT_TRUE(standard_solver.isStable());
    
    // Test entropic BGK
    config.use_entropic_bgk = true;
    SinglePhaseSolver entropic_solver(config);
    entropic_solver.setOutputDirectory(test_output_dir);
    entropic_solver.initialize();
    
    for (int i = 0; i < 10; i++) {
        EXPECT_NO_THROW(entropic_solver.step());
    }
    EXPECT_TRUE(entropic_solver.isStable());
}

TEST_F(SinglePhaseTest, FullSimulation) {
    // Test short simulation
    config.max_timesteps = 50;
    config.output_interval = 25;
    
    SinglePhaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    
    // Test successful simulation
    EXPECT_TRUE(solver.runSimulation());
    
    // Check that output files were created
    EXPECT_TRUE(std::filesystem::exists(test_output_dir + "/test_velocity_t0.csv"));
    EXPECT_TRUE(std::filesystem::exists(test_output_dir + "/test_velocity_t25.csv"));
    EXPECT_TRUE(std::filesystem::exists(test_output_dir + "/test_velocity_t50.csv"));
}

TEST_F(SinglePhaseTest, ErrorHandling) {
    SinglePhaseSolver solver(config);
    
    // Test stepping without initialization
    EXPECT_THROW(solver.step(), std::runtime_error);
    
    // Test invalid output directory
    solver.setOutputDirectory("/invalid/path/that/cannot/be/created");
    EXPECT_THROW(solver.initialize(), std::runtime_error);
}

TEST_F(SinglePhaseTest, MemoryManagement) {
    // Test that multiple solvers can be created and destroyed without issues
    for (int i = 0; i < 5; i++) {
        auto solver = std::make_unique<SinglePhaseSolver>(config);
        solver->setOutputDirectory(test_output_dir);
        solver->initialize();
        
        for (int j = 0; j < 10; j++) {
            solver->step();
        }
        
        EXPECT_TRUE(solver->isStable());
        // Solver will be automatically destroyed at end of scope
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}