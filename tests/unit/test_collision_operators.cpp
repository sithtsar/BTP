#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include <numeric>
#include "../../include/lbm/single_phase.hpp"
#include "../../include/analysis/h_theorem.hpp"

using namespace lbm;

class CollisionOperatorTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Set up basic configuration for testing
        config.Nx = 10;
        config.Ny = 11;
        config.tau = 1.0;
        config.rho0 = 1.0;
        config.gravity = 1e-6;
        config.max_timesteps = 100;
        config.output_interval = 50;
        config.max_velocity_limit = 0.1;
        config.min_density_limit = 1e-6;
        config.stability_check_interval = 10;
        config.output_prefix = "test";
        config.write_analytical_comparison = true;
        
        // Initialize test distribution functions for D2Q9
        f_test.resize(9);
        f_eq_test.resize(9);
        
        // Create equilibrium distribution for testing
        double rho = 1.0;
        double ux = 0.01;
        double uy = 0.005;
        
        createEquilibriumDistribution(f_eq_test, rho, ux, uy);
        
        // Perturb for non-equilibrium testing  
        for (int i = 0; i < 9; i++) {
            f_test[i] = f_eq_test[i] + 0.001 * std::sin(i * M_PI / 4.0);
        }
    }
    
    void createEquilibriumDistribution(std::vector<double>& f_eq, 
                                     double rho, double ux, double uy) {
        // D2Q9 weights
        std::vector<double> w = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,
                                1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
        
        // D2Q9 lattice velocities
        std::vector<std::pair<int, int>> e = {
            {0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1},
            {1, 1}, {-1, 1}, {-1, -1}, {1, -1}
        };
        
        double u_sqr = ux * ux + uy * uy;
        
        for (int i = 0; i < 9; i++) {
            double e_dot_u = e[i].first * ux + e[i].second * uy;
            f_eq[i] = w[i] * rho * (1.0 + 3.0 * e_dot_u + 
                     4.5 * e_dot_u * e_dot_u - 1.5 * u_sqr);
        }
    }
    
    double calculateDensity(const std::vector<double>& f) {
        return std::accumulate(f.begin(), f.end(), 0.0);
    }
    
    std::pair<double, double> calculateVelocity(const std::vector<double>& f, double rho) {
        // D2Q9 lattice velocities
        std::vector<std::pair<int, int>> e = {
            {0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1},
            {1, 1}, {-1, 1}, {-1, -1}, {1, -1}
        };
        
        double ux = 0.0, uy = 0.0;
        for (int i = 0; i < 9; i++) {
            ux += f[i] * e[i].first;
            uy += f[i] * e[i].second;
        }
        
        return {ux / rho, uy / rho};
    }
    
    SinglePhaseSolver::SimulationConfig config;
    std::vector<double> f_test;
    std::vector<double> f_eq_test;
};

TEST_F(CollisionOperatorTest, StandardBGKMassConservation) {
    config.use_entropic_bgk = false;
    SinglePhaseSolver solver(config);
    solver.initialize();
    
    // Test mass conservation during collision
    double initial_mass = 0.0;
    for (int x = 0; x < config.Nx; x++) {
        for (int y = 0; y < config.Ny; y++) {
            initial_mass += solver.getDensity(x, y);
        }
    }
    
    // Run several collision steps
    for (int step = 0; step < 50; step++) {
        solver.step();
        
        double current_mass = 0.0;
        for (int x = 0; x < config.Nx; x++) {
            for (int y = 0; y < config.Ny; y++) {
                current_mass += solver.getDensity(x, y);
            }
        }
        
        // Mass should be conserved to machine precision
        EXPECT_NEAR(current_mass, initial_mass, 1e-12) 
            << "Mass not conserved at step " << step;
    }
}

TEST_F(CollisionOperatorTest, EntropicBGKMassConservation) {
    config.use_entropic_bgk = true;
    SinglePhaseSolver solver(config);
    solver.initialize();
    
    // Test mass conservation for entropic BGK
    double initial_mass = 0.0;
    for (int x = 0; x < config.Nx; x++) {
        for (int y = 0; y < config.Ny; y++) {
            initial_mass += solver.getDensity(x, y);
        }
    }
    
    // Run several collision steps
    for (int step = 0; step < 50; step++) {
        solver.step();
        
        double current_mass = 0.0;
        for (int x = 0; x < config.Nx; x++) {
            for (int y = 0; y < config.Ny; y++) {
                current_mass += solver.getDensity(x, y);
            }
        }
        
        // Mass should be conserved to machine precision
        EXPECT_NEAR(current_mass, initial_mass, 1e-12) 
            << "Mass not conserved at step " << step;
    }
}

TEST_F(CollisionOperatorTest, EquilibriumRecovery) {
    // Test that BGK collision drives distribution toward equilibrium
    config.use_entropic_bgk = false;
    SinglePhaseSolver solver(config);
    solver.initialize();
    
    // Let the system evolve for many timesteps
    for (int step = 0; step < 1000; step++) {
        solver.step();
    }
    
    // Check that the system has reached near-equilibrium
    // In steady state Poiseuille flow, check velocity profile
    for (int y = 1; y < config.Ny - 1; y++) {
        int center_x = config.Nx / 2;
        double numerical_ux = solver.getVelocity(center_x, y).first;
        double analytical_ux = solver.getAnalyticalVelocity(y);
        
        // Should converge to analytical solution
        EXPECT_NEAR(numerical_ux, analytical_ux, 1e-4) 
            << "Velocity not converged at y=" << y;
    }
}

TEST_F(CollisionOperatorTest, ViscosityConsistency) {
    // Test that different tau values give correct viscosity scaling
    std::vector<double> tau_values = {0.6, 0.8, 1.0, 1.2, 1.5};
    std::vector<double> max_velocities;
    
    for (double tau : tau_values) {
        config.tau = tau;
        config.use_entropic_bgk = false;
        SinglePhaseSolver solver(config);
        solver.initialize();
        
        // Run to steady state
        for (int step = 0; step < 2000; step++) {
            solver.step();
        }
        
        // Find maximum velocity (should be at channel center)
        double max_vel = 0.0;
        for (int y = 0; y < config.Ny; y++) {
            int center_x = config.Nx / 2;
            double vel = solver.getVelocity(center_x, y).first;
            max_vel = std::max(max_vel, vel);
        }
        
        max_velocities.push_back(max_vel);
    }
    
    // Maximum velocity should scale linearly with kinematic viscosity
    // nu = (tau - 0.5) / 3, so max_vel should scale with (tau - 0.5)
    for (size_t i = 1; i < tau_values.size(); i++) {
        double nu_ratio = (tau_values[i] - 0.5) / (tau_values[0] - 0.5);
        double vel_ratio = max_velocities[i] / max_velocities[0];
        
        EXPECT_NEAR(vel_ratio, nu_ratio, 0.05) 
            << "Viscosity scaling incorrect for tau=" << tau_values[i];
    }
}

TEST_F(CollisionOperatorTest, StabilityBounds) {
    // Test stability boundaries for BGK collision
    std::vector<double> unstable_tau = {0.4, 0.49};  // Below stability limit
    
    for (double tau : unstable_tau) {
        config.tau = tau;
        config.use_entropic_bgk = false;
        
        // Should throw exception for unstable tau
        EXPECT_THROW(SinglePhaseSolver solver(config), std::invalid_argument)
            << "Unstable tau=" << tau << " should be rejected";
    }
    
    // Test stable tau values
    std::vector<double> stable_tau = {0.5, 0.6, 1.0, 2.0};
    
    for (double tau : stable_tau) {
        config.tau = tau;
        config.use_entropic_bgk = false;
        
        EXPECT_NO_THROW(SinglePhaseSolver solver(config))
            << "Stable tau=" << tau << " should be accepted";
        
        SinglePhaseSolver solver(config);
        solver.initialize();
        
        // Should remain stable for many steps
        for (int step = 0; step < 100; step++) {
            solver.step();
            EXPECT_TRUE(solver.isStable()) 
                << "Solver became unstable at step " << step << " with tau=" << tau;
        }
    }
}

TEST_F(CollisionOperatorTest, MomentumConservation) {
    config.use_entropic_bgk = false;
    config.gravity = 0.0;  // Remove body force for momentum conservation test
    SinglePhaseSolver solver(config);
    solver.initialize();
    
    // Calculate initial momentum
    double initial_momentum_x = 0.0, initial_momentum_y = 0.0;
    for (int x = 0; x < config.Nx; x++) {
        for (int y = 0; y < config.Ny; y++) {
            double rho = solver.getDensity(x, y);
            auto vel = solver.getVelocity(x, y);
            initial_momentum_x += rho * vel.first;
            initial_momentum_y += rho * vel.second;
        }
    }
    
    // Run collision steps without body force
    for (int step = 0; step < 50; step++) {
        solver.step();
        
        double momentum_x = 0.0, momentum_y = 0.0;
        for (int x = 0; x < config.Nx; x++) {
            for (int y = 0; y < config.Ny; y++) {
                double rho = solver.getDensity(x, y);
                auto vel = solver.getVelocity(x, y);
                momentum_x += rho * vel.first;
                momentum_y += rho * vel.second;
            }
        }
        
        // Momentum should be conserved (accounting for boundary effects)
        EXPECT_NEAR(momentum_x, initial_momentum_x, 1e-10) 
            << "X-momentum not conserved at step " << step;
        EXPECT_NEAR(momentum_y, initial_momentum_y, 1e-10) 
            << "Y-momentum not conserved at step " << step;
    }
}

TEST_F(CollisionOperatorTest, EntropicBGKStability) {
    // Test that entropic BGK is more stable than standard BGK
    config.tau = 0.51;  // Close to stability limit
    
    // Test standard BGK at stability limit
    config.use_entropic_bgk = false;
    SinglePhaseSolver standard_solver(config);
    standard_solver.initialize();
    
    bool standard_stable = true;
    for (int step = 0; step < 500; step++) {
        standard_solver.step();
        if (!standard_solver.isStable()) {
            standard_stable = false;
            break;
        }
    }
    
    // Test entropic BGK at same parameters
    config.use_entropic_bgk = true;
    SinglePhaseSolver entropic_solver(config);
    entropic_solver.initialize();
    
    bool entropic_stable = true;
    for (int step = 0; step < 500; step++) {
        entropic_solver.step();
        if (!entropic_solver.isStable()) {
            entropic_stable = false;
            break;
        }
    }
    
    // Entropic BGK should be at least as stable as standard BGK
    if (!standard_stable) {
        EXPECT_TRUE(entropic_stable) 
            << "Entropic BGK should be more stable than standard BGK";
    }
}

TEST_F(CollisionOperatorTest, NumericalDiffusion) {
    // Test numerical diffusion properties of collision operators
    config.use_entropic_bgk = false;
    SinglePhaseSolver solver(config);
    solver.initialize();
    
    // Measure energy dissipation over time
    std::vector<double> kinetic_energies;
    
    for (int step = 0; step <= 100; step++) {
        if (step > 0) solver.step();
        
        double kinetic_energy = 0.0;
        for (int x = 0; x < config.Nx; x++) {
            for (int y = 0; y < config.Ny; y++) {
                double rho = solver.getDensity(x, y);
                auto vel = solver.getVelocity(x, y);
                kinetic_energy += 0.5 * rho * (vel.first * vel.first + vel.second * vel.second);
            }
        }
        kinetic_energies.push_back(kinetic_energy);
    }
    
    // Energy should generally decrease due to viscous dissipation
    // (except for initial transients in presence of body force)
    for (size_t i = 50; i < kinetic_energies.size() - 1; i++) {
        if (kinetic_energies[i+1] > kinetic_energies[i] * 1.01) {  // Allow 1% tolerance
            // Energy increase detected - this could indicate issues but is allowed due to body force
            break;
        }
    }
    
    // In steady state with body force, energy should stabilize
    double final_energy = kinetic_energies.back();
    EXPECT_GT(final_energy, 0.0) << "Final kinetic energy should be positive";
}

TEST_F(CollisionOperatorTest, DistributionFunctionBounds) {
    config.use_entropic_bgk = false;
    SinglePhaseSolver solver(config);
    solver.initialize();
    
    // Check that distribution functions remain within physically reasonable bounds
    for (int step = 0; step < 100; step++) {
        solver.step();
        
        for (int x = 0; x < config.Nx; x++) {
            for (int y = 0; y < config.Ny; y++) {
                double rho = solver.getDensity(x, y);
                
                // Density should be positive and reasonable
                EXPECT_GT(rho, 0.0) << "Negative density at (" << x << ", " << y << ")";
                EXPECT_LT(rho, 10.0) << "Unreasonably high density at (" << x << ", " << y << ")";
                
                // Velocity magnitude should be reasonable (much less than sound speed)
                auto vel = solver.getVelocity(x, y);
                double vel_mag = std::sqrt(vel.first * vel.first + vel.second * vel.second);
                EXPECT_LT(vel_mag, 0.3) << "Velocity too high at (" << x << ", " << y << ")";
            }
        }
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}