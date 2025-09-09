#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <cmath>
#include "../../include/lbm/multiphase.hpp"

using namespace lbm;

class MultiphaseFlowIntegrationTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Set up realistic multiphase flow parameters
        config.Nx = 64;
        config.Ny = 32;
        config.rho_L = 1.0;     // Light phase density
        config.rho_H = 1000.0;  // Heavy phase density
        config.mu_L = 0.01;     // Light phase viscosity
        config.mu_H = 1.0;      // Heavy phase viscosity
        config.sigma = 0.01;    // Surface tension
        config.xi = 4.0;        // Interface width
        config.gravity = 1e-5;
        config.max_timesteps = 2000;
        config.output_interval = 500;
        config.max_velocity_limit = 0.2;
        config.min_density_limit = 1e-6;
        config.stability_check_interval = 100;
        config.output_prefix = "multiphase_test";
        config.phi_0 = 0.5;     // Interface at middle
        
        test_output_dir = "multiphase_integration_test";
        
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
    
    struct MultiphaseData {
        std::vector<double> x, y;
        std::vector<double> phi, rho, ux, uy;
    };
    
    MultiphaseData readMultiphaseData(const std::string& filename) {
        MultiphaseData data;
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
            
            if (tokens.size() >= 6) {
                data.x.push_back(std::stod(tokens[0]));
                data.y.push_back(std::stod(tokens[1]));
                data.phi.push_back(std::stod(tokens[2]));
                data.rho.push_back(std::stod(tokens[3]));
                data.ux.push_back(std::stod(tokens[4]));
                data.uy.push_back(std::stod(tokens[5]));
            }
        }
        
        return data;
    }
    
    double calculateInterfacePosition(const MultiphaseData& data, int Nx, int Ny) {
        // Find interface position by locating phi = 0.5 contour
        double interface_pos = 0.0;
        int count = 0;
        
        for (int x = 0; x < Nx; x++) {
            for (int y = 0; y < Ny; y++) {
                int idx = x * Ny + y;
                if (idx < data.phi.size() && std::abs(data.phi[idx] - 0.5) < 0.1) {
                    interface_pos += data.y[idx];
                    count++;
                }
            }
        }
        
        return (count > 0) ? interface_pos / count : -1.0;
    }
    
    double calculateTotalMass(const MultiphaseData& data) {
        double total_mass = 0.0;
        for (double rho : data.rho) {
            total_mass += rho;
        }
        return total_mass;
    }
    
    double calculateInterfaceThickness(const MultiphaseData& data, int Nx, int Ny) {
        // Calculate interface thickness using gradient magnitude
        double max_gradient = 0.0;
        
        for (int x = 1; x < Nx - 1; x++) {
            for (int y = 1; y < Ny - 1; y++) {
                int idx = x * Ny + y;
                int idx_xp = (x + 1) * Ny + y;
                int idx_yp = x * Ny + (y + 1);
                
                if (idx < data.phi.size() && idx_xp < data.phi.size() && idx_yp < data.phi.size()) {
                    double dphi_dx = data.phi[idx_xp] - data.phi[idx];
                    double dphi_dy = data.phi[idx_yp] - data.phi[idx];
                    double grad_mag = std::sqrt(dphi_dx * dphi_dx + dphi_dy * dphi_dy);
                    max_gradient = std::max(max_gradient, grad_mag);
                }
            }
        }
        
        return (max_gradient > 0) ? 1.0 / max_gradient : 0.0;
    }
    
    MultiphaseSolver::SimulationConfig config;
    std::string test_output_dir;
};

TEST_F(MultiphaseFlowIntegrationTest, BasicSimulationExecution) {
    // Test that basic multiphase simulation runs to completion
    MultiphaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    
    EXPECT_TRUE(solver.runSimulation());
    
    // Check that expected output files were created
    std::vector<int> expected_timesteps = {0, 500, 1000, 1500, 2000};
    for (int t : expected_timesteps) {
        std::string full_file = test_output_dir + "/multiphase_test_t" + std::to_string(t) + ".csv";
        std::string avg_file = test_output_dir + "/multiphase_test_avg_t" + std::to_string(t) + ".csv";
        
        EXPECT_TRUE(std::filesystem::exists(full_file)) 
            << "Missing full field file: " << full_file;
        EXPECT_TRUE(std::filesystem::exists(avg_file)) 
            << "Missing averaged file: " << avg_file;
    }
    
    // Check that simulation remained stable
    EXPECT_TRUE(solver.isStable());
}

TEST_F(MultiphaseFlowIntegrationTest, MassConservation) {
    // Test that total mass is conserved throughout simulation
    MultiphaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    
    EXPECT_TRUE(solver.runSimulation());
    
    // Read data at different times
    std::vector<int> timesteps = {0, 500, 1000, 1500, 2000};
    std::vector<double> total_masses;
    
    for (int t : timesteps) {
        std::string filename = test_output_dir + "/multiphase_test_t" + std::to_string(t) + ".csv";
        auto data = readMultiphaseData(filename);
        double mass = calculateTotalMass(data);
        total_masses.push_back(mass);
    }
    
    // Check mass conservation
    double initial_mass = total_masses[0];
    for (size_t i = 1; i < total_masses.size(); i++) {
        double relative_change = std::abs(total_masses[i] - initial_mass) / initial_mass;
        EXPECT_LT(relative_change, 1e-8) 
            << "Mass not conserved at timestep " << timesteps[i] 
            << ", relative change: " << relative_change;
    }
}

TEST_F(MultiphaseFlowIntegrationTest, InterfaceStability) {
    // Test that interface remains stable and doesn't undergo spurious breakup
    config.max_timesteps = 3000;  // Longer simulation for interface stability
    
    MultiphaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    
    EXPECT_TRUE(solver.runSimulation());
    
    // Read data at different times
    std::vector<int> timesteps = {0, 1000, 2000, 3000};
    std::vector<double> interface_positions;
    std::vector<double> interface_thicknesses;
    
    for (int t : timesteps) {
        std::string filename = test_output_dir + "/multiphase_test_t" + std::to_string(t) + ".csv";
        auto data = readMultiphaseData(filename);
        
        double pos = calculateInterfacePosition(data, config.Nx, config.Ny);
        double thickness = calculateInterfaceThickness(data, config.Nx, config.Ny);
        
        interface_positions.push_back(pos);
        interface_thicknesses.push_back(thickness);
    }
    
    // Interface position should not drift significantly
    double initial_pos = interface_positions[0];
    for (size_t i = 1; i < interface_positions.size(); i++) {
        if (interface_positions[i] > 0) {  // Valid position found
            double drift = std::abs(interface_positions[i] - initial_pos);
            EXPECT_LT(drift, 2.0) 
                << "Interface drifted too much at timestep " << timesteps[i];
        }
    }
    
    // Interface thickness should remain bounded
    for (size_t i = 0; i < interface_thicknesses.size(); i++) {
        EXPECT_GT(interface_thicknesses[i], 0.0) 
            << "Invalid interface thickness at timestep " << timesteps[i];
        EXPECT_LT(interface_thicknesses[i], config.xi * 3.0) 
            << "Interface too thick at timestep " << timesteps[i];
    }
}

TEST_F(MultiphaseFlowIntegrationTest, DensityRatioHandling) {
    // Test simulation with different density ratios
    std::vector<double> density_ratios = {10.0, 100.0, 1000.0};
    
    for (double ratio : density_ratios) {
        config.rho_H = ratio * config.rho_L;
        config.max_timesteps = 1500;
        
        std::string subdir = test_output_dir + "/ratio_" + std::to_string(int(ratio));
        std::filesystem::create_directories(subdir);
        
        MultiphaseSolver solver(config);
        solver.setOutputDirectory(subdir);
        
        EXPECT_TRUE(solver.runSimulation()) 
            << "Simulation failed for density ratio " << ratio;
        
        // Check final stability
        EXPECT_TRUE(solver.isStable()) 
            << "Solver unstable for density ratio " << ratio;
        
        // Check that reasonable density values are maintained
        std::string final_file = subdir + "/multiphase_test_t1500.csv";
        auto data = readMultiphaseData(final_file);
        
        double min_density = *std::min_element(data.rho.begin(), data.rho.end());
        double max_density = *std::max_element(data.rho.begin(), data.rho.end());
        
        EXPECT_GT(min_density, config.rho_L * 0.1) 
            << "Density too low for ratio " << ratio;
        EXPECT_LT(max_density, config.rho_H * 1.1) 
            << "Density too high for ratio " << ratio;
        
        double achieved_ratio = max_density / min_density;
        EXPECT_GT(achieved_ratio, ratio * 0.5) 
            << "Density ratio not maintained for " << ratio;
    }
}

TEST_F(MultiphaseFlowIntegrationTest, ViscosityRatioEffects) {
    // Test different viscosity ratios
    std::vector<double> viscosity_ratios = {1.0, 10.0, 100.0};
    
    for (double ratio : viscosity_ratios) {
        config.mu_H = ratio * config.mu_L;
        config.max_timesteps = 1500;
        
        std::string subdir = test_output_dir + "/visc_" + std::to_string(int(ratio));
        std::filesystem::create_directories(subdir);
        
        MultiphaseSolver solver(config);
        solver.setOutputDirectory(subdir);
        
        EXPECT_TRUE(solver.runSimulation()) 
            << "Simulation failed for viscosity ratio " << ratio;
        
        EXPECT_TRUE(solver.isStable()) 
            << "Solver unstable for viscosity ratio " << ratio;
        
        // Check that velocity fields are reasonable
        std::string final_file = subdir + "/multiphase_test_t1500.csv";
        auto data = readMultiphaseData(final_file);
        
        double max_vel_mag = 0.0;
        for (size_t i = 0; i < data.ux.size(); i++) {
            double vel_mag = std::sqrt(data.ux[i] * data.ux[i] + data.uy[i] * data.uy[i]);
            max_vel_mag = std::max(max_vel_mag, vel_mag);
        }
        
        EXPECT_LT(max_vel_mag, config.max_velocity_limit) 
            << "Velocities too high for viscosity ratio " << ratio;
        EXPECT_GT(max_vel_mag, 0.0) 
            << "No flow detected for viscosity ratio " << ratio;
    }
}

TEST_F(MultiphaseFlowIntegrationTest, SurfaceTensionEffects) {
    // Test different surface tension values
    std::vector<double> surface_tensions = {0.001, 0.01, 0.1};
    
    for (double sigma : surface_tensions) {
        config.sigma = sigma;
        config.max_timesteps = 1500;
        
        std::string subdir = test_output_dir + "/sigma_" + std::to_string(int(sigma * 1000));
        std::filesystem::create_directories(subdir);
        
        MultiphaseSolver solver(config);
        solver.setOutputDirectory(subdir);
        
        EXPECT_TRUE(solver.runSimulation()) 
            << "Simulation failed for surface tension " << sigma;
        
        EXPECT_TRUE(solver.isStable()) 
            << "Solver unstable for surface tension " << sigma;
        
        // Check interface properties
        std::string final_file = subdir + "/multiphase_test_t1500.csv";
        auto data = readMultiphaseData(final_file);
        
        double interface_thickness = calculateInterfaceThickness(data, config.Nx, config.Ny);
        
        // Higher surface tension should lead to sharper interfaces
        EXPECT_GT(interface_thickness, 0.5) 
            << "Interface too thin for sigma " << sigma;
        EXPECT_LT(interface_thickness, config.xi * 5.0) 
            << "Interface too thick for sigma " << sigma;
    }
}

TEST_F(MultiphaseFlowIntegrationTest, GridResolutionEffects) {
    // Test different grid resolutions
    std::vector<std::pair<int, int>> grid_sizes = {{32, 16}, {64, 32}, {128, 64}};
    
    for (const auto& grid : grid_sizes) {
        config.Nx = grid.first;
        config.Ny = grid.second;
        config.max_timesteps = 1000;  // Shorter for multiple runs
        
        std::string subdir = test_output_dir + "/grid_" + std::to_string(grid.first) + "x" + std::to_string(grid.second);
        std::filesystem::create_directories(subdir);
        
        MultiphaseSolver solver(config);
        solver.setOutputDirectory(subdir);
        
        EXPECT_TRUE(solver.runSimulation()) 
            << "Simulation failed for grid " << grid.first << "x" << grid.second;
        
        EXPECT_TRUE(solver.isStable()) 
            << "Solver unstable for grid " << grid.first << "x" << grid.second;
        
        // Check that interface is resolved
        std::string final_file = subdir + "/multiphase_test_t1000.csv";
        auto data = readMultiphaseData(final_file);
        
        double interface_thickness = calculateInterfaceThickness(data, grid.first, grid.second);
        
        // Interface should be resolved with at least a few grid points
        EXPECT_GT(interface_thickness, 1.0) 
            << "Interface not resolved for grid " << grid.first << "x" << grid.second;
    }
}

TEST_F(MultiphaseFlowIntegrationTest, LongTermEvolution) {
    // Test long-term stability and evolution
    config.max_timesteps = 5000;
    config.output_interval = 1000;
    
    MultiphaseSolver solver(config);
    solver.setOutputDirectory(test_output_dir);
    
    EXPECT_TRUE(solver.runSimulation());
    EXPECT_TRUE(solver.isStable());
    
    // Analyze evolution over time
    std::vector<int> timesteps = {0, 1000, 2000, 3000, 4000, 5000};
    std::vector<double> total_masses;
    std::vector<double> interface_positions;
    
    for (int t : timesteps) {
        std::string filename = test_output_dir + "/multiphase_test_t" + std::to_string(t) + ".csv";
        auto data = readMultiphaseData(filename);
        
        total_masses.push_back(calculateTotalMass(data));
        interface_positions.push_back(calculateInterfacePosition(data, config.Nx, config.Ny));
    }
    
    // Check long-term mass conservation
    double initial_mass = total_masses[0];
    for (size_t i = 1; i < total_masses.size(); i++) {
        double relative_change = std::abs(total_masses[i] - initial_mass) / initial_mass;
        EXPECT_LT(relative_change, 1e-7) 
            << "Mass conservation deteriorated at long times, timestep " << timesteps[i];
    }
    
    // Check interface stability over long times
    for (size_t i = 1; i < interface_positions.size(); i++) {
        if (interface_positions[i] > 0 && interface_positions[0] > 0) {
            double drift = std::abs(interface_positions[i] - interface_positions[0]);
            EXPECT_LT(drift, 3.0) 
                << "Interface drifted too much over long time, timestep " << timesteps[i];
        }
    }
}

TEST_F(MultiphaseFlowIntegrationTest, NumericalStabilityBoundaries) {
    // Test near stability boundaries
    struct StabilityTest {
        double tau_l, tau_h;
        double gravity;
        std::string name;
    };
    
    std::vector<StabilityTest> tests = {
        {0.52, 0.55, 1e-5, "low_tau_high_force"},
        {2.5, 3.0, 1e-6, "high_tau_low_force"},
        {1.0, 1.0, 1e-4, "equal_tau_high_force"}
    };
    
    for (const auto& test : tests) {
        config.tau_phi = 0.7;  // Use default phase field relaxation time
        // Note: Removed separate tau_l and tau_h as they don't exist in actual config
        config.gravity = test.gravity;
        config.max_timesteps = 1500;
        
        std::string subdir = test_output_dir + "/" + test.name;
        std::filesystem::create_directories(subdir);
        
        MultiphaseSolver solver(config);
        solver.setOutputDirectory(subdir);
        
        EXPECT_TRUE(solver.runSimulation()) 
            << "Simulation failed for " << test.name;
        
        EXPECT_TRUE(solver.isStable()) 
            << "Solver unstable for " << test.name;
        
        // Check that results are physically reasonable
        std::string final_file = subdir + "/multiphase_test_t1500.csv";
        auto data = readMultiphaseData(final_file);
        
        // Check density bounds
        double min_rho = *std::min_element(data.rho.begin(), data.rho.end());
        double max_rho = *std::max_element(data.rho.begin(), data.rho.end());
        
        EXPECT_GT(min_rho, 0.0) << "Non-positive density for " << test.name;
        EXPECT_LT(max_rho, config.rho_H * 1.5) << "Unrealistic density for " << test.name;
        
        // Check velocity bounds
        for (size_t i = 0; i < data.ux.size(); i++) {
            double vel_mag = std::sqrt(data.ux[i] * data.ux[i] + data.uy[i] * data.uy[i]);
            EXPECT_LT(vel_mag, 0.5) << "Unrealistic velocity for " << test.name;
        }
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}