#include <gtest/gtest.h>
#include <cmath>
#include <memory>
#include <fstream>
#include "../../include/analysis/h_theorem.hpp"

using namespace lbm::analysis;

class HTheoremAnalyzerTest : public ::testing::Test {
protected:
    void SetUp() override {
        Nx = 10;
        Ny = 10;
        
        // Allocate test data arrays
        allocateTestArrays();
        
        // Initialize with simple test data
        initializeTestData();
        
        // Create analyzer with default config
        HTheoremAnalyzer::AnalysisConfig config;
        config.enable_statistical_analysis = true;
        config.enable_monotonicity_check = true;
        config.monotonicity_tolerance = 1e-12;
        
        analyzer = std::make_unique<HTheoremAnalyzer>(Nx, Ny, config);
    }
    
    void TearDown() override {
        deallocateTestArrays();
    }
    
    void allocateTestArrays() {
        // Allocate distribution functions
        f = new double**[Q];
        for (int q = 0; q < Q; ++q) {
            f[q] = new double*[Nx];
            for (int x = 0; x < Nx; ++x) {
                f[q][x] = new double[Ny];
            }
        }
        
        // Allocate macroscopic variables
        rho = new double*[Nx];
        ux = new double*[Nx];
        uy = new double*[Nx];
        for (int x = 0; x < Nx; ++x) {
            rho[x] = new double[Ny];
            ux[x] = new double[Ny];
            uy[x] = new double[Ny];
        }
    }
    
    void deallocateTestArrays() {
        // Deallocate distribution functions
        for (int q = 0; q < Q; ++q) {
            for (int x = 0; x < Nx; ++x) {
                delete[] f[q][x];
            }
            delete[] f[q];
        }
        delete[] f;
        
        // Deallocate macroscopic variables
        for (int x = 0; x < Nx; ++x) {
            delete[] rho[x];
            delete[] ux[x];
            delete[] uy[x];
        }
        delete[] rho;
        delete[] ux;
        delete[] uy;
    }
    
    void initializeTestData() {
        // Initialize with equilibrium distributions
        for (int x = 0; x < Nx; ++x) {
            for (int y = 0; y < Ny; ++y) {
                rho[x][y] = 1.0;
                ux[x][y] = 0.01 * x / Nx;  // Small velocity gradient
                uy[x][y] = 0.0;
                
                // Calculate equilibrium distributions
                for (int q = 0; q < Q; ++q) {
                    f[q][x][y] = calculateEquilibrium(q, rho[x][y], ux[x][y], uy[x][y]);
                }
            }
        }
    }
    
    double calculateEquilibrium(int direction, double rho_val, double ux_val, double uy_val) {
        static constexpr int ex[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
        static constexpr int ey[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
        static constexpr double w[9] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 
                                       1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
        static constexpr double cs2 = 1.0 / 3.0;
        
        double uxeq = ux_val * ex[direction];
        double uyeq = uy_val * ey[direction];
        double u2 = ux_val * ux_val + uy_val * uy_val;
        double usq = uxeq + uyeq;
        
        return w[direction] * rho_val * (1.0 + usq / cs2 + 0.5 * usq * usq / (cs2 * cs2) - 0.5 * u2 / cs2);
    }
    
    void perturbDistributions(double perturbation_strength) {
        // Add small perturbations to test non-equilibrium behavior
        for (int x = 0; x < Nx; ++x) {
            for (int y = 0; y < Ny; ++y) {
                for (int q = 0; q < Q; ++q) {
                    f[q][x][y] += perturbation_strength * (0.5 - static_cast<double>(rand()) / RAND_MAX);
                    f[q][x][y] = std::max(f[q][x][y], 1e-15); // Ensure positivity
                }
            }
        }
    }
    
    static constexpr int Q = 9;
    int Nx, Ny;
    double*** f;
    double** rho;
    double** ux;
    double** uy;
    std::unique_ptr<HTheoremAnalyzer> analyzer;
};

TEST_F(HTheoremAnalyzerTest, ConstructorInitialization) {
    EXPECT_EQ(analyzer->getConfig().enable_statistical_analysis, true);
    EXPECT_EQ(analyzer->getConfig().enable_monotonicity_check, true);
    EXPECT_DOUBLE_EQ(analyzer->getConfig().monotonicity_tolerance, 1e-12);
    
    // Check that evolution data is initially empty
    const auto& evolution_data = analyzer->getEvolutionData();
    EXPECT_TRUE(evolution_data.timesteps.empty());
    EXPECT_TRUE(evolution_data.h_values.empty());
}

TEST_F(HTheoremAnalyzerTest, EquilibriumHFunctionStandardBGK) {
    // For equilibrium distributions, H-function should be zero (or very close)
    auto metrics = analyzer->calculateHFunctionStandardBGK(f, rho, ux, uy, 0);
    
    EXPECT_EQ(metrics.timestep, 0);
    EXPECT_NEAR(metrics.h_function_value, 0.0, 1e-10);
    EXPECT_FALSE(metrics.has_nan_values);
    EXPECT_FALSE(metrics.has_inf_values);
    EXPECT_GT(metrics.total_mass, 0.0);
    EXPECT_GE(metrics.kinetic_energy, 0.0);
}

TEST_F(HTheoremAnalyzerTest, EquilibriumHFunctionEntropicBGK) {
    // For equilibrium distributions, entropic H-function should be well-defined
    auto metrics = analyzer->calculateHFunctionEntropicBGK(f, rho, ux, uy, 0);
    
    EXPECT_EQ(metrics.timestep, 0);
    EXPECT_FALSE(metrics.has_nan_values);
    EXPECT_FALSE(metrics.has_inf_values);
    EXPECT_GT(metrics.total_mass, 0.0);
    EXPECT_GE(metrics.kinetic_energy, 0.0);
    
    // For equilibrium, H-function should be finite and well-defined
    EXPECT_FALSE(std::isnan(metrics.h_function_value));
    EXPECT_FALSE(std::isinf(metrics.h_function_value));
}

TEST_F(HTheoremAnalyzerTest, NonEquilibriumHFunction) {
    // Perturb distributions away from equilibrium
    perturbDistributions(0.01);
    
    auto metrics = analyzer->calculateHFunctionStandardBGK(f, rho, ux, uy, 0);
    
    // H-function should be positive for non-equilibrium distributions
    EXPECT_GT(metrics.h_function_value, 0.0);
    EXPECT_FALSE(metrics.has_nan_values);
    EXPECT_FALSE(metrics.has_inf_values);
}

TEST_F(HTheoremAnalyzerTest, HFunctionEvolutionTracking) {
    // Calculate H-function at multiple timesteps
    auto metrics1 = analyzer->calculateHFunctionStandardBGK(f, rho, ux, uy, 0);
    analyzer->recordHEvolution(metrics1);
    
    // Perturb slightly to simulate evolution
    perturbDistributions(0.005);
    auto metrics2 = analyzer->calculateHFunctionStandardBGK(f, rho, ux, uy, 1);
    analyzer->recordHEvolution(metrics2);
    
    // Check evolution data
    const auto& evolution_data = analyzer->getEvolutionData();
    EXPECT_EQ(evolution_data.timesteps.size(), 2);
    EXPECT_EQ(evolution_data.h_values.size(), 2);
    EXPECT_EQ(evolution_data.h_changes.size(), 2);
    
    EXPECT_EQ(evolution_data.timesteps[0], 0);
    EXPECT_EQ(evolution_data.timesteps[1], 1);
    
    // Second timestep should have calculated change
    EXPECT_NE(metrics2.h_function_change, 0.0);
}

TEST_F(HTheoremAnalyzerTest, MonotonicityVerification) {
    // Create a sequence where H decreases
    std::vector<double> h_values = {10.0, 8.0, 6.0, 4.0, 2.0};
    
    // Manually populate evolution data
    for (size_t i = 0; i < h_values.size(); ++i) {
        HTheoremAnalyzer::HTheoremMetrics metrics;
        metrics.timestep = static_cast<int>(i);
        metrics.h_function_value = h_values[i];
        if (i > 0) {
            metrics.h_function_change = h_values[i] - h_values[i-1];
            metrics.monotonic_decrease = (metrics.h_function_change <= 0.0);
        }
        analyzer->recordHEvolution(metrics);
    }
    
    EXPECT_TRUE(analyzer->verifyMonotonicDecrease(1e-12));
}

TEST_F(HTheoremAnalyzerTest, MonotonicityViolationDetection) {
    // Create a sequence with a violation
    std::vector<double> h_values = {10.0, 8.0, 9.0, 7.0, 5.0}; // Increase at index 2
    
    for (size_t i = 0; i < h_values.size(); ++i) {
        HTheoremAnalyzer::HTheoremMetrics metrics;
        metrics.timestep = static_cast<int>(i);
        metrics.h_function_value = h_values[i];
        if (i > 0) {
            metrics.h_function_change = h_values[i] - h_values[i-1];
            metrics.monotonic_decrease = (metrics.h_function_change <= 0.0);
        }
        analyzer->recordHEvolution(metrics);
    }
    
    EXPECT_FALSE(analyzer->verifyMonotonicDecrease(1e-12));
    
    const auto& evolution_data = analyzer->getEvolutionData();
    EXPECT_GT(evolution_data.monotonicity_violation_count, 0);
    EXPECT_FALSE(evolution_data.overall_monotonic);
}

TEST_F(HTheoremAnalyzerTest, EntropyProductionCalculation) {
    double current_h = 5.0;
    double previous_h = 6.0;
    double dt = 1.0;
    
    double entropy_production = analyzer->calculateEntropyProductionRate(current_h, previous_h, dt);
    
    // Entropy production should be positive when H decreases
    EXPECT_GT(entropy_production, 0.0);
    EXPECT_DOUBLE_EQ(entropy_production, 1.0); // (6.0 - 5.0) / 1.0
}

TEST_F(HTheoremAnalyzerTest, SpatialStatisticsAnalysis) {
    // Create a test H-function field
    double** h_field = new double*[Nx];
    for (int x = 0; x < Nx; ++x) {
        h_field[x] = new double[Ny];
        for (int y = 0; y < Ny; ++y) {
            h_field[x][y] = 1.0 + 0.1 * x + 0.05 * y; // Linear variation
        }
    }
    
    auto stats = analyzer->analyzeSpatialDistribution(h_field);
    
    EXPECT_GT(stats.mean, 0.0);
    EXPECT_GE(stats.variance, 0.0);
    EXPECT_GE(stats.standard_deviation, 0.0);
    EXPECT_LE(stats.min_value, stats.max_value);
    EXPECT_LE(stats.min_value, stats.mean);
    EXPECT_GE(stats.max_value, stats.mean);
    
    // Clean up
    for (int x = 0; x < Nx; ++x) {
        delete[] h_field[x];
    }
    delete[] h_field;
}

TEST_F(HTheoremAnalyzerTest, NumericalAnomalyDetection) {
    // Introduce NaN values
    f[0][0][0] = std::numeric_limits<double>::quiet_NaN();
    
    auto metrics = analyzer->calculateHFunctionStandardBGK(f, rho, ux, uy, 0);
    
    EXPECT_TRUE(metrics.has_nan_values);
}

TEST_F(HTheoremAnalyzerTest, ConfigurationManagement) {
    HTheoremAnalyzer::AnalysisConfig new_config;
    new_config.enable_statistical_analysis = false;
    new_config.monotonicity_tolerance = 1e-6;
    new_config.output_prefix = "test_output";
    
    analyzer->setConfig(new_config);
    
    const auto& config = analyzer->getConfig();
    EXPECT_FALSE(config.enable_statistical_analysis);
    EXPECT_DOUBLE_EQ(config.monotonicity_tolerance, 1e-6);
    EXPECT_EQ(config.output_prefix, "test_output");
}

TEST_F(HTheoremAnalyzerTest, ResetFunctionality) {
    // Add some data
    auto metrics = analyzer->calculateHFunctionStandardBGK(f, rho, ux, uy, 0);
    analyzer->recordHEvolution(metrics);
    
    EXPECT_FALSE(analyzer->getEvolutionData().timesteps.empty());
    
    // Reset
    analyzer->reset();
    
    EXPECT_TRUE(analyzer->getEvolutionData().timesteps.empty());
    EXPECT_TRUE(analyzer->getEvolutionData().h_values.empty());
}

TEST_F(HTheoremAnalyzerTest, DiagnosticReportGeneration) {
    // Add some test data
    for (int t = 0; t < 5; ++t) {
        HTheoremAnalyzer::HTheoremMetrics metrics;
        metrics.timestep = t;
        metrics.h_function_value = 10.0 - t * 2.0; // Decreasing H
        metrics.h_function_change = (t > 0) ? -2.0 : 0.0;
        metrics.monotonic_decrease = true;
        metrics.entropy_production = 2.0;
        analyzer->recordHEvolution(metrics);
    }
    
    std::string report = analyzer->generateDiagnosticReport();
    
    EXPECT_FALSE(report.empty());
    EXPECT_NE(report.find("H-theorem Analysis Diagnostic Report"), std::string::npos);
    EXPECT_NE(report.find("Total timesteps analyzed: 5"), std::string::npos);
    EXPECT_NE(report.find("Monotonicity violations: 0"), std::string::npos);
}

TEST_F(HTheoremAnalyzerTest, MassConservationTracking) {
    auto metrics = analyzer->calculateHFunctionStandardBGK(f, rho, ux, uy, 0);
    
    // Total mass should equal sum of densities
    double expected_mass = 0.0;
    for (int x = 0; x < Nx; ++x) {
        for (int y = 0; y < Ny; ++y) {
            expected_mass += rho[x][y];
        }
    }
    
    EXPECT_NEAR(metrics.total_mass, expected_mass, 1e-12);
}

TEST_F(HTheoremAnalyzerTest, KineticEnergyCalculation) {
    auto metrics = analyzer->calculateHFunctionStandardBGK(f, rho, ux, uy, 0);
    
    // Calculate expected kinetic energy
    double expected_ke = 0.0;
    for (int x = 0; x < Nx; ++x) {
        for (int y = 0; y < Ny; ++y) {
            expected_ke += 0.5 * rho[x][y] * (ux[x][y] * ux[x][y] + uy[x][y] * uy[x][y]);
        }
    }
    
    EXPECT_NEAR(metrics.kinetic_energy, expected_ke, 1e-12);
}

// Test fixture for file I/O operations
class HTheoremAnalyzerFileTest : public HTheoremAnalyzerTest {
protected:
    void SetUp() override {
        HTheoremAnalyzerTest::SetUp();
        test_filename = "test_h_evolution.csv";
    }
    
    void TearDown() override {
        // Clean up test file
        std::remove(test_filename.c_str());
        HTheoremAnalyzerTest::TearDown();
    }
    
    std::string test_filename;
};

TEST_F(HTheoremAnalyzerFileTest, WriteEvolutionData) {
    // Add some test data
    for (int t = 0; t < 3; ++t) {
        HTheoremAnalyzer::HTheoremMetrics metrics;
        metrics.timestep = t;
        metrics.h_function_value = 10.0 - t;
        metrics.h_function_change = (t > 0) ? -1.0 : 0.0;
        metrics.entropy_production = 1.0;
        metrics.h_function_variance = 0.1 * t;
        analyzer->recordHEvolution(metrics);
    }
    
    bool success = analyzer->writeEvolutionData(test_filename);
    EXPECT_TRUE(success);
    
    // Verify file exists and has content
    std::ifstream file(test_filename);
    EXPECT_TRUE(file.is_open());
    
    std::string line;
    int line_count = 0;
    while (std::getline(file, line)) {
        line_count++;
    }
    
    // Should have header lines plus data lines
    EXPECT_GE(line_count, 5); // 2 header lines + 3 data lines
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}