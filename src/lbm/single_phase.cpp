#include "../../include/lbm/single_phase.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace lbm {

SinglePhaseSolver::SinglePhaseSolver(const SimulationConfig& config)
    : SolverBase(config.Nx, config.Ny), config_(config), f_(nullptr), f_temp_(nullptr),
      rho_(nullptr), ux_(nullptr), uy_(nullptr), is_stable_(true),
      last_stability_check_(0), memory_allocated_(false) {
    
    // Validate configuration
    if (config_.tau <= 0.5) {
        throw std::invalid_argument("Relaxation time tau must be > 0.5 for stability");
    }
    
    if (config_.rho0 <= 0) {
        throw std::invalid_argument("Reference density must be positive");
    }
    
    if (config_.max_timesteps <= 0) {
        throw std::invalid_argument("Maximum timesteps must be positive");
    }
    
    // Calculate derived parameters
    omega_ = 1.0 / config_.tau;
    cs2_ = 1.0 / 3.0;
    nu_ = cs2_ * (config_.tau - 0.5);
    
    // Configure stability monitor
    StabilityConfig stability_config;
    stability_config.max_velocity = config_.max_velocity_limit;
    stability_config.min_density = config_.min_density_limit;
    stability_config.max_density = 10.0 * config_.rho0;  // Allow 10x reference density
    stability_config.check_interval = config_.stability_check_interval;
    stability_monitor_.updateConfig(stability_config);
    
    // Initialize H-theorem analyzer if enabled
    if (config_.enable_h_theorem_analysis) {
        analysis::HTheoremAnalyzer::AnalysisConfig h_config;
        h_config.enable_spatial_analysis = true;
        h_config.enable_monotonicity_check = true;
        h_config.enable_entropy_production = true;
        h_config.enable_statistical_analysis = true;
        h_config.analysis_interval = config_.h_analysis_interval;
        h_config.output_prefix = config_.output_prefix + "_h_theorem";
        
        h_analyzer_ = std::make_unique<analysis::HTheoremAnalyzer>(Nx_, Ny_, h_config);
    }
    
    std::cout << "SinglePhaseSolver initialized:" << std::endl;
    std::cout << "  Grid: " << Nx_ << "x" << Ny_ << std::endl;
    std::cout << "  tau: " << config_.tau << ", omega: " << omega_ << std::endl;
    std::cout << "  Kinematic viscosity: " << nu_ << std::endl;
    std::cout << "  Collision operator: " << (config_.use_entropic_bgk ? "Entropic BGK" : "Standard BGK") << std::endl;
    std::cout << "  H-theorem analysis: " << (config_.enable_h_theorem_analysis ? "Enabled" : "Disabled") << std::endl;
}

SinglePhaseSolver::~SinglePhaseSolver() {
    deallocateMemory();
}

void SinglePhaseSolver::allocateMemory() {
    if (memory_allocated_) {
        return;  // Already allocated
    }
    
    try {
        // Allocate distribution functions
        f_ = new double**[Q];
        f_temp_ = new double**[Q];
        for (int i = 0; i < Q; i++) {
            f_[i] = new double*[Nx_];
            f_temp_[i] = new double*[Nx_];
            for (int x = 0; x < Nx_; x++) {
                f_[i][x] = new double[Ny_];
                f_temp_[i][x] = new double[Ny_];
            }
        }
        
        // Allocate macroscopic variables
        rho_ = new double*[Nx_];
        ux_ = new double*[Nx_];
        uy_ = new double*[Nx_];
        for (int x = 0; x < Nx_; x++) {
            rho_[x] = new double[Ny_];
            ux_[x] = new double[Ny_];
            uy_[x] = new double[Ny_];
        }
        
        memory_allocated_ = true;
        
    } catch (const std::bad_alloc& e) {
        deallocateMemory();  // Clean up partial allocation
        throw std::runtime_error("Failed to allocate memory for solver arrays");
    }
}

void SinglePhaseSolver::deallocateMemory() {
    if (!memory_allocated_) {
        return;  // Nothing to deallocate
    }
    
    // Deallocate distribution functions
    if (f_) {
        for (int i = 0; i < Q; i++) {
            if (f_[i]) {
                for (int x = 0; x < Nx_; x++) {
                    delete[] f_[i][x];
                }
                delete[] f_[i];
            }
        }
        delete[] f_;
        f_ = nullptr;
    }
    
    if (f_temp_) {
        for (int i = 0; i < Q; i++) {
            if (f_temp_[i]) {
                for (int x = 0; x < Nx_; x++) {
                    delete[] f_temp_[i][x];
                }
                delete[] f_temp_[i];
            }
        }
        delete[] f_temp_;
        f_temp_ = nullptr;
    }
    
    // Deallocate macroscopic variables
    if (rho_) {
        for (int x = 0; x < Nx_; x++) {
            delete[] rho_[x];
        }
        delete[] rho_;
        rho_ = nullptr;
    }
    
    if (ux_) {
        for (int x = 0; x < Nx_; x++) {
            delete[] ux_[x];
        }
        delete[] ux_;
        ux_ = nullptr;
    }
    
    if (uy_) {
        for (int x = 0; x < Nx_; x++) {
            delete[] uy_[x];
        }
        delete[] uy_;
        uy_ = nullptr;
    }
    
    memory_allocated_ = false;
}

void SinglePhaseSolver::initialize() {
    allocateMemory();
    initializeDistributions();
    updateMacroscopicVariables();
    
    // Validate output directory
    if (!utils::IOUtils::validateOutputDirectory(output_directory_)) {
        throw std::runtime_error("Cannot create or access output directory: " + output_directory_);
    }
    
    std::cout << "Solver initialized successfully" << std::endl;
}

void SinglePhaseSolver::initializeDistributions() {
    for (int x = 0; x < Nx_; x++) {
        for (int y = 0; y < Ny_; y++) {
            for (int i = 0; i < Q; i++) {
                f_[i][x][y] = w[i] * config_.rho0;
                f_temp_[i][x][y] = 0.0;
            }
        }
    }
}

void SinglePhaseSolver::step() {
    if (!memory_allocated_) {
        throw std::runtime_error("Solver not initialized - call initialize() first");
    }
    
    // Perform LBM steps
    collision();
    streaming();
    applyBoundaryConditions();
    updateMacroscopicVariables();
    
    // Update timestep
    updateTimestep(getCurrentTimestep() + 1);
    
    // Perform H-theorem analysis if enabled
    if (config_.enable_h_theorem_analysis && h_analyzer_ && 
        getCurrentTimestep() % config_.h_analysis_interval == 0) {
        performHTheoremAnalysis(getCurrentTimestep());
    }
    
    // Check stability periodically
    if (getCurrentTimestep() % config_.stability_check_interval == 0) {
        is_stable_ = checkNumericalStability();
        last_stability_check_ = getCurrentTimestep();
        
        if (!is_stable_) {
            std::cerr << "Numerical instability detected at timestep " << getCurrentTimestep() << std::endl;
        }
    }
}

bool SinglePhaseSolver::isStable() const {
    return is_stable_;
}

void SinglePhaseSolver::collision() {
    for (int x = 0; x < Nx_; x++) {
        for (int y = 0; y < Ny_; y++) {
            if (config_.use_entropic_bgk) {
                entropicBGKCollision(x, y);
            } else {
                standardBGKCollision(x, y);
            }
        }
    }
}

void SinglePhaseSolver::standardBGKCollision(int x, int y) {
    // Calculate macroscopic variables
    double rho_local = 0.0;
    double ux_local = 0.0, uy_local = 0.0;
    
    for (int i = 0; i < Q; i++) {
        rho_local += f_[i][x][y];
        ux_local += f_[i][x][y] * ex[i];
        uy_local += f_[i][x][y] * ey[i];
    }
    
    // Prevent division by zero
    if (rho_local < config_.min_density_limit) {
        rho_local = config_.rho0;
        ux_local = 0.0;
        uy_local = 0.0;
    } else {
        ux_local /= rho_local;
        uy_local /= rho_local;
    }
    
    // Add body force to velocity
    ux_local += 0.5 * config_.gravity / rho_local;
    
    // Calculate equilibrium distributions
    double feq[Q];
    for (int i = 0; i < Q; i++) {
        feq[i] = utils::MathUtils::calculateEquilibrium(i, rho_local, ux_local, uy_local);
    }
    
    // Calculate force terms
    for (int i = 0; i < Q; i++) {
        double force_term = calculateForceTerm(i, ux_local, uy_local);
        
        // BGK collision with force
        f_temp_[i][x][y] = f_[i][x][y] - omega_ * (f_[i][x][y] - feq[i]) + force_term;
        
        // Ensure non-negative distributions
        if (f_temp_[i][x][y] < 0.0) {
            f_temp_[i][x][y] = 1e-12;  // Small positive value
        }
    }
}

void SinglePhaseSolver::entropicBGKCollision(int x, int y) {
    // Calculate macroscopic variables
    double rho_local = 0.0;
    double ux_local = 0.0, uy_local = 0.0;
    
    for (int i = 0; i < Q; i++) {
        rho_local += f_[i][x][y];
        ux_local += f_[i][x][y] * ex[i];
        uy_local += f_[i][x][y] * ey[i];
    }
    
    // Prevent division by zero
    if (rho_local < config_.min_density_limit) {
        rho_local = config_.rho0;
        ux_local = 0.0;
        uy_local = 0.0;
    } else {
        ux_local /= rho_local;
        uy_local /= rho_local;
    }
    
    // Add body force to velocity
    ux_local += 0.5 * config_.gravity / rho_local;
    
    // Calculate equilibrium distributions
    double feq[Q];
    for (int i = 0; i < Q; i++) {
        feq[i] = utils::MathUtils::calculateEquilibrium(i, rho_local, ux_local, uy_local);
    }
    
    // Find optimal alpha parameter for entropy conservation
    double alpha = findOptimalAlpha(x, y, feq);
    
    // Apply entropic collision with force terms
    for (int i = 0; i < Q; i++) {
        double force_term = calculateForceTerm(i, ux_local, uy_local);
        
        // Entropic BGK collision
        f_temp_[i][x][y] = f_[i][x][y] + alpha * omega_ * (feq[i] - f_[i][x][y]) + force_term;
        
        // Ensure non-negative distributions
        if (f_temp_[i][x][y] < 0.0) {
            f_temp_[i][x][y] = 1e-12;  // Small positive value
        }
    }
}

double SinglePhaseSolver::calculateForceTerm(int direction, double ux, double uy) const {
    double eu = ex[direction] * ux + ey[direction] * uy;
    double F_dot_e = config_.gravity * ex[direction];
    
    return w[direction] * (1.0 - 0.5 * omega_) * 
           (3.0 * F_dot_e + 9.0 * F_dot_e * eu - 3.0 * config_.gravity * ux);
}

double SinglePhaseSolver::findOptimalAlpha(int x, int y, const double* feq) const {
    // Calculate current entropy H = sum_i f_i * ln(f_i/w_i)
    double H_current = calculateHFunction(x, y);
    
    // Use bisection method to find alpha such that entropy is conserved
    double alpha_min = 0.0, alpha_max = 2.0;
    double alpha = 1.0;
    const int max_iter = 20;
    const double tolerance = 1e-8;
    
    for (int iter = 0; iter < max_iter; iter++) {
        alpha = 0.5 * (alpha_min + alpha_max);
        
        // Calculate H at intermediate state
        double H_test = 0.0;
        for (int i = 0; i < Q; i++) {
            double f_test = f_[i][x][y] + alpha * (feq[i] - f_[i][x][y]);
            if (f_test > 1e-12) {
                H_test += f_test * std::log(f_test / w[i]);
            }
        }
        
        double diff = H_test - H_current;
        if (std::abs(diff) < tolerance) {
            break;
        }
        
        if (diff > 0) {
            alpha_max = alpha;
        } else {
            alpha_min = alpha;
        }
    }
    
    return alpha;
}

void SinglePhaseSolver::streaming() {
    // Create temporary array for streaming
    static double*** f_new = nullptr;
    static bool first_call = true;
    
    if (first_call) {
        f_new = new double**[Q];
        for (int i = 0; i < Q; i++) {
            f_new[i] = new double*[Nx_];
            for (int x = 0; x < Nx_; x++) {
                f_new[i][x] = new double[Ny_];
            }
        }
        first_call = false;
    }
    
    // Initialize f_new to zero
    for (int i = 0; i < Q; i++) {
        for (int x = 0; x < Nx_; x++) {
            for (int y = 0; y < Ny_; y++) {
                f_new[i][x][y] = 0.0;
            }
        }
    }
    
    // Stream from all points
    for (int x = 0; x < Nx_; x++) {
        for (int y = 0; y < Ny_; y++) {
            for (int i = 0; i < Q; i++) {
                int x_next = (x + ex[i] + Nx_) % Nx_;  // Periodic in x
                int y_next = y + ey[i];
                
                if (y_next >= 0 && y_next < Ny_) {
                    f_new[i][x_next][y_next] = f_temp_[i][x][y];
                }
            }
        }
    }
    
    // Copy back to f array
    for (int i = 0; i < Q; i++) {
        for (int x = 0; x < Nx_; x++) {
            for (int y = 0; y < Ny_; y++) {
                f_[i][x][y] = f_new[i][x][y];
            }
        }
    }
}

void SinglePhaseSolver::applyBoundaryConditions() {
    // Bottom wall (y=0) - bounce back
    for (int x = 0; x < Nx_; x++) {
        f_[2][x][0] = f_[4][x][0];  // north <- south
        f_[5][x][0] = f_[7][x][0];  // NE <- SW
        f_[6][x][0] = f_[8][x][0];  // NW <- SE
    }
    
    // Top wall (y=Ny-1) - bounce back
    for (int x = 0; x < Nx_; x++) {
        f_[4][x][Ny_-1] = f_[2][x][Ny_-1];  // south <- north
        f_[7][x][Ny_-1] = f_[5][x][Ny_-1];  // SW <- NE
        f_[8][x][Ny_-1] = f_[6][x][Ny_-1];  // SE <- NW
    }
}

void SinglePhaseSolver::updateMacroscopicVariables() {
    for (int x = 0; x < Nx_; x++) {
        for (int y = 0; y < Ny_; y++) {
            // Calculate density
            rho_[x][y] = 0.0;
            for (int i = 0; i < Q; i++) {
                rho_[x][y] += f_[i][x][y];
            }
            
            // Calculate velocity
            ux_[x][y] = 0.0;
            uy_[x][y] = 0.0;
            
            if (rho_[x][y] > config_.min_density_limit) {
                for (int i = 0; i < Q; i++) {
                    ux_[x][y] += f_[i][x][y] * ex[i];
                    uy_[x][y] += f_[i][x][y] * ey[i];
                }
                ux_[x][y] /= rho_[x][y];
                uy_[x][y] /= rho_[x][y];
            }
        }
    }
}

bool SinglePhaseSolver::checkNumericalStability() {
    // Create lattice data structure for stability checking
    LatticeData lattice_data(Nx_, Ny_);
    
    // Copy data to lattice structure
    for (int x = 0; x < Nx_; ++x) {
        for (int y = 0; y < Ny_; ++y) {
            lattice_data.density[x][y] = rho_[x][y];
            lattice_data.velocity_x[x][y] = ux_[x][y];
            lattice_data.velocity_y[x][y] = uy_[x][y];
            
            for (int i = 0; i < Q; ++i) {
                lattice_data.distribution[x][y][i] = f_[i][x][y];
            }
        }
    }
    
    // Check stability using monitor
    StabilityMetrics metrics = stability_monitor_.checkStability(lattice_data, getCurrentTimestep());
    
    return metrics.isStable();
}

double SinglePhaseSolver::calculateHFunction(int x, int y) const {
    double H = 0.0;
    for (int i = 0; i < Q; i++) {
        if (f_[i][x][y] > 1e-12) {
            H += f_[i][x][y] * std::log(f_[i][x][y] / w[i]);
        }
    }
    return H;
}

double SinglePhaseSolver::getDensity(int x, int y) const {
    if (x < 0 || x >= Nx_ || y < 0 || y >= Ny_ || !memory_allocated_) {
        return 0.0;
    }
    return rho_[x][y];
}

std::pair<double, double> SinglePhaseSolver::getVelocity(int x, int y) const {
    if (x < 0 || x >= Nx_ || y < 0 || y >= Ny_ || !memory_allocated_) {
        return {0.0, 0.0};
    }
    return {ux_[x][y], uy_[x][y]};
}

double SinglePhaseSolver::getAnalyticalVelocity(int y) const {
    return utils::MathUtils::analyticalPoiseuilleVelocity(y, Ny_ - 1, config_.gravity, nu_);
}

double SinglePhaseSolver::getHFunction(int x, int y) const {
    if (x < 0 || x >= Nx_ || y < 0 || y >= Ny_ || !memory_allocated_) {
        return 0.0;
    }
    return calculateHFunction(x, y);
}

void SinglePhaseSolver::writeOutput(int timestep) {
    writeVelocityProfile(timestep);
    if (config_.use_entropic_bgk) {
        writeHFunctionData(timestep);
    }
    if (config_.enable_h_theorem_analysis && h_analyzer_) {
        writeHTheoremResults(timestep);
    }
}

void SinglePhaseSolver::writeVelocityProfile(int timestep) {
    std::string filename = output_directory_ + "/" + 
                          utils::IOUtils::generateTimestepFilename(config_.output_prefix + "_velocity", timestep);
    
    std::vector<int> y_positions;
    std::vector<double> velocities;
    std::vector<double> analytical_velocities;
    
    for (int y = 0; y < Ny_; y++) {
        // Calculate average velocity across x-direction
        double avg_ux = 0.0;
        for (int x = 0; x < Nx_; x++) {
            avg_ux += ux_[x][y];
        }
        avg_ux /= Nx_;
        
        y_positions.push_back(y);
        velocities.push_back(avg_ux);
        
        if (config_.write_analytical_comparison) {
            analytical_velocities.push_back(getAnalyticalVelocity(y));
        }
    }
    
    utils::IOUtils::writeVelocityProfile(filename, y_positions, velocities, analytical_velocities);
}

void SinglePhaseSolver::writeHFunctionData(int timestep) {
    std::string filename = output_directory_ + "/" + 
                          utils::IOUtils::generateTimestepFilename(config_.output_prefix + "_h_y", timestep);
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Warning: Could not open H-function output file: " << filename << std::endl;
        return;
    }
    
    file << "y,H_y\n";
    
    for (int y = 0; y < Ny_; y++) {
        double H_y = 0.0;
        for (int x = 0; x < Nx_; x++) {
            H_y += calculateHFunction(x, y);
        }
        file << y << "," << utils::IOUtils::formatDouble(H_y) << "\n";
    }
    
    file.close();
}

bool SinglePhaseSolver::runSimulation() {
    try {
        initialize();
        
        std::cout << "Starting simulation with " << config_.max_timesteps << " timesteps..." << std::endl;
        std::cout << "Expected max velocity: " << getAnalyticalVelocity(Ny_/2) << std::endl;
        
        for (int t = 0; t <= config_.max_timesteps; t++) {
            step();
            
            if (!isStable()) {
                std::cerr << "Simulation terminated due to numerical instability at timestep " << t << std::endl;
                return false;
            }
            
            if (t % config_.output_interval == 0) {
                auto history = stability_monitor_.getStabilityHistory();
                double max_vel = history.empty() ? 0.0 : history.back().max_velocity;
                utils::IOUtils::printProgress(t, config_.max_timesteps, 
                                            "Max velocity: " + std::to_string(max_vel));
                writeOutput(t);
            }
        }
        
        std::cout << "Simulation completed successfully!" << std::endl;
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Simulation failed with error: " << e.what() << std::endl;
        return false;
    }
}

void SinglePhaseSolver::performHTheoremAnalysis(int timestep) {
    if (!h_analyzer_ || !memory_allocated_) {
        return;
    }
    
    analysis::HTheoremAnalyzer::HTheoremMetrics metrics;
    
    if (config_.use_entropic_bgk) {
        metrics = h_analyzer_->calculateHFunctionEntropicBGK(f_, rho_, ux_, uy_, timestep);
    } else {
        metrics = h_analyzer_->calculateHFunctionStandardBGK(f_, rho_, ux_, uy_, timestep);
    }
    
    // Record metrics for evolution tracking
    h_analyzer_->recordHEvolution(metrics);
    
    // Check for monotonicity violations and log warnings
    if (!metrics.monotonic_decrease && timestep > 0) {
        std::cout << "WARNING: H-theorem monotonicity violation detected at timestep " 
                  << timestep << std::endl;
    }
    
    // Check for numerical anomalies
    if (metrics.has_nan_values || metrics.has_inf_values) {
        std::cerr << "ERROR: Numerical anomalies detected in H-function at timestep " 
                  << timestep << std::endl;
        is_stable_ = false;
    }
}

void SinglePhaseSolver::writeHTheoremResults(int timestep) {
    if (!h_analyzer_) {
        return;
    }
    
    // Write evolution data periodically
    if (timestep % (config_.output_interval * 10) == 0 || timestep == config_.max_timesteps) {
        std::string evolution_filename = output_directory_ + "/" + 
                                       config_.output_prefix + "_h_evolution.csv";
        h_analyzer_->writeEvolutionData(evolution_filename);
        
        // Write diagnostic report
        std::string report_filename = output_directory_ + "/" + 
                                    config_.output_prefix + "_h_report.txt";
        std::ofstream report_file(report_filename);
        if (report_file.is_open()) {
            report_file << h_analyzer_->generateDiagnosticReport();
            report_file.close();
        }
    }
}

analysis::HTheoremAnalyzer::HTheoremMetrics SinglePhaseSolver::getLatestHMetrics() const {
    if (!h_analyzer_) {
        return analysis::HTheoremAnalyzer::HTheoremMetrics{};
    }
    return h_analyzer_->getLatestMetrics();
}

} // namespace lbm