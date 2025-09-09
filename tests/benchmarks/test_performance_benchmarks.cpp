#include <gtest/gtest.h>
#include <chrono>
#include <vector>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include "../../include/lbm/single_phase.hpp"
#include "../../include/lbm/multiphase.hpp"

using namespace lbm;
using namespace std::chrono;

class PerformanceBenchmarkTest : public ::testing::Test {
protected:
    void SetUp() override {
        test_output_dir = "benchmark_results";
        
        // Clean up and create output directory
        if (std::filesystem::exists(test_output_dir)) {
            std::filesystem::remove_all(test_output_dir);
        }
        std::filesystem::create_directories(test_output_dir);
        
        // Initialize basic configurations
        setupSinglePhaseConfig();
        setupMultiphaseConfig();
    }
    
    void TearDown() override {
        // Keep benchmark results for analysis
        // std::filesystem::remove_all(test_output_dir);
    }
    
    void setupSinglePhaseConfig() {
        single_config.Nx = 32;
        single_config.Ny = 33;
        single_config.tau = 1.0;
        single_config.rho0 = 1.0;
        single_config.gravity = 1e-6;
        single_config.max_timesteps = 1000;
        single_config.output_interval = 1000;  // Minimal output for benchmarking
        single_config.use_entropic_bgk = false;
        single_config.max_velocity_limit = 0.1;
        single_config.min_density_limit = 1e-6;
        single_config.stability_check_interval = 100;
        single_config.output_prefix = "benchmark";
        single_config.write_analytical_comparison = false;
        single_config.enable_h_theorem_analysis = false;
    }
    
    void setupMultiphaseConfig() {
        multi_config.Nx = 64;
        multi_config.Ny = 32;
        multi_config.rho_L = 1.0;
        multi_config.rho_H = 1000.0;
        multi_config.mu_L = 0.01;
        multi_config.mu_H = 1.0;
        multi_config.sigma = 0.01;
        multi_config.xi = 4.0;
        multi_config.gravity = 1e-5;
        multi_config.max_timesteps = 500;  // Shorter for heavy multiphase benchmarks
        multi_config.output_interval = 500;
        multi_config.max_velocity_limit = 0.2;
        multi_config.min_density_limit = 1e-6;
        multi_config.stability_check_interval = 100;
        multi_config.output_prefix = "multiphase_benchmark";
        multi_config.phi_0 = 0.5;
    }
    
    struct BenchmarkResult {
        std::string test_name;
        double runtime_seconds;
        double timesteps_per_second;
        double memory_usage_mb;
        int grid_points;
        int timesteps;
        double efficiency_score;  // timesteps * grid_points / runtime
    };
    
    double measureRuntime(std::function<bool()> simulation_func) {
        auto start = high_resolution_clock::now();
        bool success = simulation_func();
        auto end = high_resolution_clock::now();
        
        EXPECT_TRUE(success) << "Simulation failed during benchmark";
        
        auto duration = duration_cast<microseconds>(end - start);
        return duration.count() / 1e6;  // Convert to seconds
    }
    
    void saveBenchmarkResults(const std::vector<BenchmarkResult>& results) {
        std::ofstream file(test_output_dir + "/benchmark_results.csv");
        file << "Test Name,Runtime (s),Timesteps/s,Memory (MB),Grid Points,Timesteps,Efficiency Score\n";
        
        for (const auto& result : results) {
            file << result.test_name << ","
                 << std::fixed << std::setprecision(6) << result.runtime_seconds << ","
                 << std::fixed << std::setprecision(2) << result.timesteps_per_second << ","
                 << std::fixed << std::setprecision(2) << result.memory_usage_mb << ","
                 << result.grid_points << ","
                 << result.timesteps << ","
                 << std::scientific << std::setprecision(2) << result.efficiency_score << "\n";
        }
    }
    
    double estimateMemoryUsage(int Nx, int Ny, bool is_multiphase = false) {
        // Rough memory estimation
        int grid_points = Nx * Ny;
        double mb_per_grid_point = is_multiphase ? 0.5 : 0.2;  // Rough estimates
        return grid_points * mb_per_grid_point;
    }
    
    SinglePhaseSolver::SimulationConfig single_config;
    MultiphaseSolver::SimulationConfig multi_config;
    std::string test_output_dir;
    std::vector<BenchmarkResult> all_results;
};

TEST_F(PerformanceBenchmarkTest, SinglePhaseStandardBGKBaseline) {
    // Baseline performance test for standard BGK
    single_config.use_entropic_bgk = false;
    single_config.max_timesteps = 2000;
    
    auto simulation_func = [&]() {
        SinglePhaseSolver solver(single_config);
        solver.setOutputDirectory(test_output_dir + "/sp_standard");
        return solver.runSimulation();
    };
    
    double runtime = measureRuntime(simulation_func);
    
    BenchmarkResult result;
    result.test_name = "SinglePhase_StandardBGK_Baseline";
    result.runtime_seconds = runtime;
    result.timesteps_per_second = single_config.max_timesteps / runtime;
    result.memory_usage_mb = estimateMemoryUsage(single_config.Nx, single_config.Ny);
    result.grid_points = single_config.Nx * single_config.Ny;
    result.timesteps = single_config.max_timesteps;
    result.efficiency_score = (result.timesteps * result.grid_points) / runtime;
    
    all_results.push_back(result);
    
    // Performance expectations
    EXPECT_LT(runtime, 30.0) << "Standard BGK baseline too slow";
    EXPECT_GT(result.timesteps_per_second, 50.0) << "Insufficient timestep throughput";
    
    std::cout << "Standard BGK Baseline: " << result.timesteps_per_second 
              << " timesteps/sec, " << runtime << "s total\n";
}

TEST_F(PerformanceBenchmarkTest, SinglePhaseEntropicBGKComparison) {
    // Compare entropic BGK performance to standard BGK
    single_config.use_entropic_bgk = true;
    single_config.max_timesteps = 2000;
    
    auto simulation_func = [&]() {
        SinglePhaseSolver solver(single_config);
        solver.setOutputDirectory(test_output_dir + "/sp_entropic");
        return solver.runSimulation();
    };
    
    double runtime = measureRuntime(simulation_func);
    
    BenchmarkResult result;
    result.test_name = "SinglePhase_EntropicBGK";
    result.runtime_seconds = runtime;
    result.timesteps_per_second = single_config.max_timesteps / runtime;
    result.memory_usage_mb = estimateMemoryUsage(single_config.Nx, single_config.Ny);
    result.grid_points = single_config.Nx * single_config.Ny;
    result.timesteps = single_config.max_timesteps;
    result.efficiency_score = (result.timesteps * result.grid_points) / runtime;
    
    all_results.push_back(result);
    
    // Entropic BGK should be slower but not excessively so
    EXPECT_LT(runtime, 60.0) << "Entropic BGK too slow";
    EXPECT_GT(result.timesteps_per_second, 30.0) << "Entropic BGK throughput too low";
    
    // Compare with standard BGK if available
    auto standard_it = std::find_if(all_results.begin(), all_results.end(),
        [](const BenchmarkResult& r) { return r.test_name.find("StandardBGK") != std::string::npos; });
    
    if (standard_it != all_results.end()) {
        double performance_ratio = result.timesteps_per_second / standard_it->timesteps_per_second;
        EXPECT_GT(performance_ratio, 0.3) << "Entropic BGK more than 3x slower than standard";
        EXPECT_LT(performance_ratio, 2.0) << "Entropic BGK unexpectedly faster than standard";
    }
    
    std::cout << "Entropic BGK: " << result.timesteps_per_second 
              << " timesteps/sec, " << runtime << "s total\n";
}

TEST_F(PerformanceBenchmarkTest, GridSizeScaling) {
    // Test performance scaling with grid size
    std::vector<std::pair<int, int>> grid_sizes = {
        {16, 17}, {32, 33}, {64, 65}, {128, 129}
    };
    
    for (const auto& grid : grid_sizes) {
        single_config.Nx = grid.first;
        single_config.Ny = grid.second;
        single_config.use_entropic_bgk = false;
        single_config.max_timesteps = 1000;  // Fixed timesteps for scaling test
        
        std::string test_name = "GridScaling_" + std::to_string(grid.first) + "x" + std::to_string(grid.second);
        
        auto simulation_func = [&]() {
            SinglePhaseSolver solver(single_config);
            solver.setOutputDirectory(test_output_dir + "/" + test_name);
            return solver.runSimulation();
        };
        
        double runtime = measureRuntime(simulation_func);
        
        BenchmarkResult result;
        result.test_name = test_name;
        result.runtime_seconds = runtime;
        result.timesteps_per_second = single_config.max_timesteps / runtime;
        result.memory_usage_mb = estimateMemoryUsage(grid.first, grid.second);
        result.grid_points = grid.first * grid.second;
        result.timesteps = single_config.max_timesteps;
        result.efficiency_score = (result.timesteps * result.grid_points) / runtime;
        
        all_results.push_back(result);
        
        std::cout << "Grid " << grid.first << "x" << grid.second 
                  << ": " << result.timesteps_per_second << " timesteps/sec\n";
    }
    
    // Check scaling behavior
    if (all_results.size() >= 2) {
        // Find grid scaling results
        std::vector<BenchmarkResult*> scaling_results;
        for (auto& result : all_results) {
            if (result.test_name.find("GridScaling") != std::string::npos) {
                scaling_results.push_back(&result);
            }
        }
        
        if (scaling_results.size() >= 2) {
            // Sort by grid points
            std::sort(scaling_results.begin(), scaling_results.end(),
                [](const BenchmarkResult* a, const BenchmarkResult* b) {
                    return a->grid_points < b->grid_points;
                });
            
            // Check that performance scales reasonably with grid size
            for (size_t i = 1; i < scaling_results.size(); i++) {
                double grid_ratio = static_cast<double>(scaling_results[i]->grid_points) / scaling_results[0]->grid_points;
                double time_ratio = scaling_results[i]->runtime_seconds / scaling_results[0]->runtime_seconds;
                
                // Runtime should scale roughly linearly with grid points (O(N))
                EXPECT_LT(time_ratio, grid_ratio * 2.0) 
                    << "Performance scaling worse than O(N^2) from " 
                    << scaling_results[0]->test_name << " to " << scaling_results[i]->test_name;
                EXPECT_GT(time_ratio, grid_ratio * 0.5) 
                    << "Unexpectedly good scaling from " 
                    << scaling_results[0]->test_name << " to " << scaling_results[i]->test_name;
            }
        }
    }
}

TEST_F(PerformanceBenchmarkTest, MultiphaseBaseline) {
    // Baseline multiphase performance
    auto simulation_func = [&]() {
        MultiphaseSolver solver(multi_config);
        solver.setOutputDirectory(test_output_dir + "/mp_baseline");
        return solver.runSimulation();
    };
    
    double runtime = measureRuntime(simulation_func);
    
    BenchmarkResult result;
    result.test_name = "Multiphase_Baseline";
    result.runtime_seconds = runtime;
    result.timesteps_per_second = multi_config.max_timesteps / runtime;
    result.memory_usage_mb = estimateMemoryUsage(multi_config.Nx, multi_config.Ny, true);
    result.grid_points = multi_config.Nx * multi_config.Ny;
    result.timesteps = multi_config.max_timesteps;
    result.efficiency_score = (result.timesteps * result.grid_points) / runtime;
    
    all_results.push_back(result);
    
    // Multiphase should be slower than single-phase but still reasonable
    EXPECT_LT(runtime, 120.0) << "Multiphase baseline too slow";
    EXPECT_GT(result.timesteps_per_second, 5.0) << "Multiphase throughput too low";
    
    std::cout << "Multiphase Baseline: " << result.timesteps_per_second 
              << " timesteps/sec, " << runtime << "s total\n";
}

TEST_F(PerformanceBenchmarkTest, DensityRatioPerformanceImpact) {
    // Test performance impact of different density ratios
    std::vector<double> density_ratios = {10.0, 100.0, 1000.0};
    
    for (double ratio : density_ratios) {
        multi_config.rho_H = ratio * multi_config.rho_L;
        multi_config.max_timesteps = 300;  // Shorter for multiple runs
        
        std::string test_name = "Multiphase_DensityRatio_" + std::to_string(int(ratio));
        
        auto simulation_func = [&]() {
            MultiphaseSolver solver(multi_config);
            solver.setOutputDirectory(test_output_dir + "/" + test_name);
            return solver.runSimulation();
        };
        
        double runtime = measureRuntime(simulation_func);
        
        BenchmarkResult result;
        result.test_name = test_name;
        result.runtime_seconds = runtime;
        result.timesteps_per_second = multi_config.max_timesteps / runtime;
        result.memory_usage_mb = estimateMemoryUsage(multi_config.Nx, multi_config.Ny, true);
        result.grid_points = multi_config.Nx * multi_config.Ny;
        result.timesteps = multi_config.max_timesteps;
        result.efficiency_score = (result.timesteps * result.grid_points) / runtime;
        
        all_results.push_back(result);
        
        // Performance shouldn't degrade drastically with density ratio
        EXPECT_LT(runtime, 60.0) << "Performance too poor for density ratio " << ratio;
        EXPECT_GT(result.timesteps_per_second, 3.0) 
            << "Throughput too low for density ratio " << ratio;
        
        std::cout << "Density ratio " << ratio << ": " << result.timesteps_per_second 
                  << " timesteps/sec\n";
    }
}

TEST_F(PerformanceBenchmarkTest, MemoryUsageEstimation) {
    // Test memory usage scaling
    std::vector<std::pair<int, int>> grids = {{32, 33}, {64, 65}, {128, 129}};
    
    for (const auto& grid : grids) {
        single_config.Nx = grid.first;
        single_config.Ny = grid.second;
        single_config.max_timesteps = 100;  // Short run for memory test
        single_config.use_entropic_bgk = false;
        
        auto simulation_func = [&]() {
            SinglePhaseSolver solver(single_config);
            solver.setOutputDirectory(test_output_dir + "/memory_test");
            return solver.runSimulation();
        };
        
        double runtime = measureRuntime(simulation_func);
        
        // Check that simulation completes without memory issues
        EXPECT_GT(runtime, 0.0) << "Simulation failed for grid " << grid.first << "x" << grid.second;
        
        int grid_points = grid.first * grid.second;
        double estimated_memory = estimateMemoryUsage(grid.first, grid.second);
        
        // Memory usage should scale roughly linearly with grid points
        EXPECT_LT(estimated_memory, grid_points * 1.0) << "Memory estimate too high";
        EXPECT_GT(estimated_memory, grid_points * 0.01) << "Memory estimate too low";
        
        std::cout << "Grid " << grid.first << "x" << grid.second 
                  << ": ~" << estimated_memory << " MB estimated\n";
    }
}

TEST_F(PerformanceBenchmarkTest, StabilityMonitoringOverhead) {
    // Test performance impact of stability monitoring
    single_config.max_timesteps = 1500;
    
    // Test with frequent stability checks
    single_config.stability_check_interval = 10;
    single_config.use_entropic_bgk = false;
    
    auto frequent_check_func = [&]() {
        SinglePhaseSolver solver(single_config);
        solver.setOutputDirectory(test_output_dir + "/frequent_checks");
        return solver.runSimulation();
    };
    
    double runtime_frequent = measureRuntime(frequent_check_func);
    
    // Test with infrequent stability checks
    single_config.stability_check_interval = 500;
    
    auto infrequent_check_func = [&]() {
        SinglePhaseSolver solver(single_config);
        solver.setOutputDirectory(test_output_dir + "/infrequent_checks");
        return solver.runSimulation();
    };
    
    double runtime_infrequent = measureRuntime(infrequent_check_func);
    
    // Record both results
    BenchmarkResult frequent_result;
    frequent_result.test_name = "StabilityChecks_Frequent";
    frequent_result.runtime_seconds = runtime_frequent;
    frequent_result.timesteps_per_second = single_config.max_timesteps / runtime_frequent;
    frequent_result.memory_usage_mb = estimateMemoryUsage(single_config.Nx, single_config.Ny);
    frequent_result.grid_points = single_config.Nx * single_config.Ny;
    frequent_result.timesteps = single_config.max_timesteps;
    frequent_result.efficiency_score = (frequent_result.timesteps * frequent_result.grid_points) / runtime_frequent;
    
    BenchmarkResult infrequent_result;
    infrequent_result.test_name = "StabilityChecks_Infrequent";
    infrequent_result.runtime_seconds = runtime_infrequent;
    infrequent_result.timesteps_per_second = single_config.max_timesteps / runtime_infrequent;
    infrequent_result.memory_usage_mb = estimateMemoryUsage(single_config.Nx, single_config.Ny);
    infrequent_result.grid_points = single_config.Nx * single_config.Ny;
    infrequent_result.timesteps = single_config.max_timesteps;
    infrequent_result.efficiency_score = (infrequent_result.timesteps * infrequent_result.grid_points) / runtime_infrequent;
    
    all_results.push_back(frequent_result);
    all_results.push_back(infrequent_result);
    
    // Stability monitoring should have minimal overhead
    double overhead_ratio = runtime_frequent / runtime_infrequent;
    EXPECT_LT(overhead_ratio, 1.2) << "Stability monitoring overhead too high";
    
    std::cout << "Stability monitoring overhead: " << (overhead_ratio - 1.0) * 100 << "%\n";
}

TEST_F(PerformanceBenchmarkTest, HTheoremAnalysisOverhead) {
    // Test performance impact of H-theorem analysis
    single_config.max_timesteps = 1500;
    single_config.use_entropic_bgk = true;
    
    // Test without H-theorem analysis
    single_config.enable_h_theorem_analysis = false;
    
    auto no_h_analysis_func = [&]() {
        SinglePhaseSolver solver(single_config);
        solver.setOutputDirectory(test_output_dir + "/no_h_analysis");
        return solver.runSimulation();
    };
    
    double runtime_no_h = measureRuntime(no_h_analysis_func);
    
    // Test with H-theorem analysis
    single_config.enable_h_theorem_analysis = true;
    single_config.h_analysis_interval = 50;
    
    auto with_h_analysis_func = [&]() {
        SinglePhaseSolver solver(single_config);
        solver.setOutputDirectory(test_output_dir + "/with_h_analysis");
        return solver.runSimulation();
    };
    
    double runtime_with_h = measureRuntime(with_h_analysis_func);
    
    // Record results
    BenchmarkResult no_h_result;
    no_h_result.test_name = "HTheorem_Disabled";
    no_h_result.runtime_seconds = runtime_no_h;
    no_h_result.timesteps_per_second = single_config.max_timesteps / runtime_no_h;
    no_h_result.memory_usage_mb = estimateMemoryUsage(single_config.Nx, single_config.Ny);
    no_h_result.grid_points = single_config.Nx * single_config.Ny;
    no_h_result.timesteps = single_config.max_timesteps;
    no_h_result.efficiency_score = (no_h_result.timesteps * no_h_result.grid_points) / runtime_no_h;
    
    BenchmarkResult with_h_result;
    with_h_result.test_name = "HTheorem_Enabled";
    with_h_result.runtime_seconds = runtime_with_h;
    with_h_result.timesteps_per_second = single_config.max_timesteps / runtime_with_h;
    with_h_result.memory_usage_mb = estimateMemoryUsage(single_config.Nx, single_config.Ny);
    with_h_result.grid_points = single_config.Nx * single_config.Ny;
    with_h_result.timesteps = single_config.max_timesteps;
    with_h_result.efficiency_score = (with_h_result.timesteps * with_h_result.grid_points) / runtime_with_h;
    
    all_results.push_back(no_h_result);
    all_results.push_back(with_h_result);
    
    // H-theorem analysis should have reasonable overhead
    double overhead_ratio = runtime_with_h / runtime_no_h;
    EXPECT_LT(overhead_ratio, 1.5) << "H-theorem analysis overhead too high";
    
    std::cout << "H-theorem analysis overhead: " << (overhead_ratio - 1.0) * 100 << "%\n";
}

TEST_F(PerformanceBenchmarkTest, ComprehensivePerformanceReport) {
    // This test runs last and generates a comprehensive report
    
    // Ensure we have results to report
    EXPECT_GT(all_results.size(), 0) << "No benchmark results to report";
    
    // Save all results
    saveBenchmarkResults(all_results);
    
    // Generate summary statistics
    if (!all_results.empty()) {
        auto best_single_phase = std::max_element(all_results.begin(), all_results.end(),
            [](const BenchmarkResult& a, const BenchmarkResult& b) {
                return (a.test_name.find("SinglePhase") != std::string::npos && 
                       b.test_name.find("SinglePhase") != std::string::npos) ? 
                       a.timesteps_per_second < b.timesteps_per_second : false;
            });
        
        auto best_multiphase = std::max_element(all_results.begin(), all_results.end(),
            [](const BenchmarkResult& a, const BenchmarkResult& b) {
                return (a.test_name.find("Multiphase") != std::string::npos && 
                       b.test_name.find("Multiphase") != std::string::npos) ? 
                       a.timesteps_per_second < b.timesteps_per_second : false;
            });
        
        // Performance regression checks
        double min_acceptable_single_phase_perf = 30.0;  // timesteps/sec
        double min_acceptable_multiphase_perf = 3.0;     // timesteps/sec
        
        bool found_good_single_phase = false;
        bool found_good_multiphase = false;
        
        for (const auto& result : all_results) {
            if (result.test_name.find("SinglePhase") != std::string::npos && 
                result.timesteps_per_second >= min_acceptable_single_phase_perf) {
                found_good_single_phase = true;
            }
            if (result.test_name.find("Multiphase") != std::string::npos && 
                result.timesteps_per_second >= min_acceptable_multiphase_perf) {
                found_good_multiphase = true;
            }
        }
        
        EXPECT_TRUE(found_good_single_phase) 
            << "No single-phase configuration achieved minimum performance threshold";
        EXPECT_TRUE(found_good_multiphase) 
            << "No multiphase configuration achieved minimum performance threshold";
        
        std::cout << "\n=== PERFORMANCE BENCHMARK SUMMARY ===\n";
        std::cout << "Total tests run: " << all_results.size() << "\n";
        
        if (best_single_phase != all_results.end() && 
            best_single_phase->test_name.find("SinglePhase") != std::string::npos) {
            std::cout << "Best single-phase: " << best_single_phase->test_name 
                      << " (" << best_single_phase->timesteps_per_second << " timesteps/sec)\n";
        }
        
        if (best_multiphase != all_results.end() && 
            best_multiphase->test_name.find("Multiphase") != std::string::npos) {
            std::cout << "Best multiphase: " << best_multiphase->test_name 
                      << " (" << best_multiphase->timesteps_per_second << " timesteps/sec)\n";
        }
        
        std::cout << "Results saved to: " << test_output_dir << "/benchmark_results.csv\n";
        std::cout << "======================================\n\n";
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}