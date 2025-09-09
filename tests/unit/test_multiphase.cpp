#include <gtest/gtest.h>
#include <filesystem>
#include "../../include/lbm/multiphase.hpp"

using namespace lbm;

class MultiphaseTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Set up basic configuration
        config.Nx = 16;
        config.Ny = 8;
        config.phi_L = 0.0;
        config.phi_H = 1.0;
        config.phi_0 = 0.5;
        config.rho_L = 1.0;
        config.rho_H = 10.0;  // Smaller density ratio for testing
        config.mu_L = 0.01;
        config.mu_H = 0.1;
        config.xi = 2.0;
        config.sigma = 0.01;
        config.tau_phi = 0.7;
        config.gravity = 1e-6;
        config.max_timesteps = 100;
        config.output_interval = 50;
        config.max_velocity_limit = 0.1;
        config.min_density_limit = 1e-6;
        config.max_density_ratio = 1e3;
        config.min_tau_limit = 0.51;
        config.max_tau_limit = 5.0;
        config.stability_check_interval = 10;
        config.interface_thickness_tolerance = 0.5;
        config.phase_conservation_tolerance = 1e-3;
        config.interface_position_tolerance = 1.0;
        config.output_prefix = "test_multiphase";
        config.write_full_fields = true;
        config.write_averaged_profiles = true;
        config.write_analytical_comparison = true;
        
        test_output_dir = "test_multiphase_output";
        
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
    
    MultiphaseSolver::SimulationConfig config;
    std::string test_output_dir;
};

TEST_F(MultiphaseTest, Construction) {
    // Test valid construction
    EXPECT_NO_THROW(MultiphaseSolver solver(config));
    
    MultiphaseSolver solver(config);
    EXPECT_EQ(solver.getGridDimensions().first, config.Nx);
    EXPECT_EQ(solver.getGridDimensions().second, config.Ny);
    EXPECT_EQ(solver.getConfig().rho_L, config.rho_L);
    EXPECT_EQ(solver.getConfig().rho_H, config.rho_H);
    
    // Test invalid configurations
    auto invalid_config = config;
    invalid_config.rho_L = -1.0;  // Negative density
    EXPECT_THROW(MultiphaseSolver invalid_solver(invalid_config), std::invalid_argument);
    
    invalid_config = config;
    invalid_config.mu_H = 0.0;  // Zero viscosity
    EXPECT_THROW(MultiphaseSolver invalid_solver(invalid_config), std::invalid_argument);
    
    invalid_config = config;
    invalid_config.tau_phi = 0.3;  // Too small for stability
    EXPECT_THROW(MultiphaseSolver invalid_solver(invalid_config), std::invalid_argument);
    
    invalid_config = config;
    invalid_config.xi = 0.0;  // Zero interface thickness
    EXPECT_THROW(MultiphaseSolver invalid_solver(invalid_config), std::invalid_argument);
    
    invalid_config = config;
    invalid_config.max_timesteps = 0;  // Invalid timesteps
    EXPECT_THROW(MultiphaseSolver invalid_solver(invalid_config), std::invalid_argument);
}

TEST_F(MultiphaseTest, Initialization) {
    MultiphaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    
    // Test initialization
    EXPECT_NO_THROW(solver.initialize());
    
    // Check that phase field is initialized properly
    for (int x = 0; x < config.Nx; x++) {
        for (int y = 0; y < config.Ny; y++) {
            double phi = solver.getPhaseField(x, y);
            EXPECT_GE(phi, config.phi_L);
            EXPECT_LE(phi, config.phi_H);
            EXPECT_TRUE(std::isfinite(phi));
        }
    }
    
    // Check that density field is initialized
    for (int x = 0; x < config.Nx; x++) {
        for (int y = 0; y < config.Ny; y++) {
            double density = solver.getDensity(x, y);
            EXPECT_GT(density, 0.0);
            EXPECT_TRUE(std::isfinite(density));
        }
    }
    
    // Test that solver is initially stable
    EXPECT_TRUE(solver.isStable());
}

TEST_F(MultiphaseTest, PhaseFieldInitialization) {
    MultiphaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    solver.initialize();
    
    // Check that interface is roughly in the middle
    int center_y = config.Ny / 2;
    
    // Check that lower half tends toward light phase
    double avg_phi_lower = 0.0;
    int count_lower = 0;
    for (int x = 0; x < config.Nx; x++) {
        for (int y = 0; y < center_y; y++) {
            avg_phi_lower += solver.getPhaseField(x, y);
            count_lower++;
        }
    }
    avg_phi_lower /= count_lower;
    EXPECT_LT(avg_phi_lower, config.phi_0);  // Should be closer to light phase
    
    // Check that upper half tends toward heavy phase
    double avg_phi_upper = 0.0;
    int count_upper = 0;
    for (int x = 0; x < config.Nx; x++) {
        for (int y = center_y; y < config.Ny; y++) {
            avg_phi_upper += solver.getPhaseField(x, y);
            count_upper++;
        }
    }
    avg_phi_upper /= count_upper;
    EXPECT_GT(avg_phi_upper, config.phi_0);  // Should be closer to heavy phase
}

TEST_F(MultiphaseTest, Timestepping) {
    MultiphaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    solver.initialize();
    
    // Test initial timestep
    EXPECT_EQ(static_cast<const MultiphaseSolver&>(solver).getCurrentTimestep(), 0);
    
    // Test stepping
    EXPECT_NO_THROW(solver.step());
    EXPECT_EQ(static_cast<const MultiphaseSolver&>(solver).getCurrentTimestep(), 1);
    
    EXPECT_NO_THROW(solver.step());
    EXPECT_EQ(static_cast<const MultiphaseSolver&>(solver).getCurrentTimestep(), 2);
    
    // Test that solver remains stable after a few steps
    for (int i = 0; i < 10; i++) {
        solver.step();
    }
    EXPECT_TRUE(solver.isStable());
}

TEST_F(MultiphaseTest, MacroscopicVariables) {
    MultiphaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    solver.initialize();
    
    // Test field access
    for (int x = 0; x < config.Nx; x++) {
        for (int y = 0; y < config.Ny; y++) {
            // Phase field
            double phi = solver.getPhaseField(x, y);
            EXPECT_GE(phi, config.phi_L);
            EXPECT_LE(phi, config.phi_H);
            EXPECT_TRUE(std::isfinite(phi));
            
            // Density
            double density = solver.getDensity(x, y);
            EXPECT_GT(density, 0.0);
            EXPECT_TRUE(std::isfinite(density));
            
            // Velocity
            auto velocity = solver.getVelocity(x, y);
            EXPECT_TRUE(std::isfinite(velocity.first));
            EXPECT_TRUE(std::isfinite(velocity.second));
            
            // Pressure
            double pressure = solver.getPressure(x, y);
            EXPECT_TRUE(std::isfinite(pressure));
        }
    }
    
    // Test out-of-bounds access
    EXPECT_EQ(solver.getPhaseField(-1, 5), 0.0);
    EXPECT_EQ(solver.getDensity(config.Nx, 5), 0.0);
    
    auto invalid_velocity = solver.getVelocity(5, -1);
    EXPECT_EQ(invalid_velocity.first, 0.0);
    EXPECT_EQ(invalid_velocity.second, 0.0);
    
    EXPECT_EQ(solver.getPressure(5, config.Ny), 0.0);
}

TEST_F(MultiphaseTest, InterfaceTracking) {
    MultiphaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    solver.initialize();
    
    // Test interface position detection
    for (int x = 0; x < config.Nx; x++) {
        double interface_pos = solver.getInterfacePosition(x);
        EXPECT_GE(interface_pos, 0.0);
        EXPECT_LT(interface_pos, config.Ny);
    }
    
    // Test out-of-bounds interface position
    EXPECT_EQ(solver.getInterfacePosition(-1), -1.0);
    EXPECT_EQ(solver.getInterfacePosition(config.Nx), -1.0);
}

TEST_F(MultiphaseTest, PhaseConservation) {
    MultiphaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    solver.initialize();
    
    // Get initial phase masses
    auto initial_masses = solver.calculatePhaseMasses();
    EXPECT_GT(initial_masses.first, 0.0);   // Light phase mass
    EXPECT_GT(initial_masses.second, 0.0);  // Heavy phase mass
    
    // Run a few timesteps
    for (int i = 0; i < 20; i++) {
        solver.step();
    }
    
    // Check that masses are approximately conserved
    auto current_masses = solver.calculatePhaseMasses();
    double light_change = std::abs(current_masses.first - initial_masses.first) / initial_masses.first;
    double heavy_change = std::abs(current_masses.second - initial_masses.second) / initial_masses.second;
    
    // Allow some tolerance for numerical errors
    EXPECT_LT(light_change, 0.1);  // 10% tolerance
    EXPECT_LT(heavy_change, 0.1);  // 10% tolerance
}

TEST_F(MultiphaseTest, AnalyticalComparison) {
    MultiphaseSolver solver(config);
    
    // Test analytical velocity calculation
    for (int y = 0; y < config.Ny; y++) {
        double analytical = solver.getAnalyticalVelocity(y);
        EXPECT_TRUE(std::isfinite(analytical));
        EXPECT_GE(analytical, 0.0);  // Should be non-negative for Poiseuille flow
    }
    
    // Test that analytical velocity profile makes sense
    // (maximum somewhere in the middle, zero at walls)
    double wall_velocity1 = solver.getAnalyticalVelocity(0);
    double wall_velocity2 = solver.getAnalyticalVelocity(config.Ny - 1);
    double mid_velocity = solver.getAnalyticalVelocity(config.Ny / 2);
    
    EXPECT_NEAR(wall_velocity1, 0.0, 1e-10);
    EXPECT_NEAR(wall_velocity2, 0.0, 1e-10);
    EXPECT_GT(mid_velocity, 0.0);
}

TEST_F(MultiphaseTest, OutputGeneration) {
    MultiphaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    solver.initialize();
    
    // Test output writing
    EXPECT_NO_THROW(solver.writeOutput(0));
    
    // Check that output files are created
    std::string full_file = test_output_dir + "/test_multiphase_t0.csv";
    std::string avg_file = test_output_dir + "/test_multiphase_avg_t0.csv";
    
    EXPECT_TRUE(std::filesystem::exists(full_file));
    EXPECT_TRUE(std::filesystem::exists(avg_file));
}

TEST_F(MultiphaseTest, StabilityMonitoring) {
    MultiphaseSolver solver(config);
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

TEST_F(MultiphaseTest, FullSimulation) {
    // Test short simulation
    config.max_timesteps = 50;
    config.output_interval = 25;
    
    MultiphaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    
    // Test successful simulation
    EXPECT_TRUE(solver.runSimulation());
    
    // Check that output files were created
    EXPECT_TRUE(std::filesystem::exists(test_output_dir + "/test_multiphase_t0.csv"));
    EXPECT_TRUE(std::filesystem::exists(test_output_dir + "/test_multiphase_t25.csv"));
    EXPECT_TRUE(std::filesystem::exists(test_output_dir + "/test_multiphase_t50.csv"));
    
    EXPECT_TRUE(std::filesystem::exists(test_output_dir + "/test_multiphase_avg_t0.csv"));
    EXPECT_TRUE(std::filesystem::exists(test_output_dir + "/test_multiphase_avg_t25.csv"));
    EXPECT_TRUE(std::filesystem::exists(test_output_dir + "/test_multiphase_avg_t50.csv"));
}

TEST_F(MultiphaseTest, ErrorHandling) {
    MultiphaseSolver solver(config);
    
    // Test stepping without initialization
    EXPECT_THROW(solver.step(), std::runtime_error);
    
    // Test invalid output directory
    solver.setOutputDirectory("/invalid/path/that/cannot/be/created");
    EXPECT_THROW(solver.initialize(), std::runtime_error);
}

TEST_F(MultiphaseTest, MemoryManagement) {
    // Test that multiple solvers can be created and destroyed without issues
    for (int i = 0; i < 3; i++) {  // Fewer iterations due to larger memory footprint
        auto solver = std::make_unique<MultiphaseSolver>(config);
        solver->setOutputDirectory(test_output_dir);
        solver->initialize();
        
        for (int j = 0; j < 5; j++) {
            solver->step();
        }
        
        EXPECT_TRUE(solver->isStable());
        // Solver will be automatically destroyed at end of scope
    }
}

TEST_F(MultiphaseTest, DensityRatioHandling) {
    // Test with extreme density ratio
    config.rho_H = 1000.0 * config.rho_L;
    config.max_density_ratio = 1e4;  // Allow higher ratio
    
    MultiphaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    solver.initialize();
    
    // Check that density interpolation works correctly
    for (int x = 0; x < config.Nx; x++) {
        for (int y = 0; y < config.Ny; y++) {
            double density = solver.getDensity(x, y);
            EXPECT_GE(density, config.rho_L);
            EXPECT_LE(density, config.rho_H);
        }
    }
    
    // Run a few steps to test stability with high density ratio
    for (int i = 0; i < 10; i++) {
        solver.step();
    }
    
    // Should still be stable (though this is challenging numerically)
    EXPECT_TRUE(solver.isStable());
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}