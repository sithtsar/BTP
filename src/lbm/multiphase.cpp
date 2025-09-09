#include "../../include/lbm/multiphase.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <limits>

namespace lbm {

MultiphaseSolver::MultiphaseSolver(const SimulationConfig& config)
    : SolverBase(config.Nx, config.Ny), config_(config), h_(nullptr), h_temp_(nullptr),
      g_(nullptr), g_temp_(nullptr), phi_(nullptr), rho_(nullptr), mu_(nullptr),
      tau_(nullptr), ux_(nullptr), uy_(nullptr), p_(nullptr), grad_phi_x_(nullptr),
      grad_phi_y_(nullptr), lapl_phi_(nullptr), is_stable_(true), interface_stable_(true),
      phase_conserved_(true), last_stability_check_(0), memory_allocated_(false) {
    
    // Validate configuration
    if (config_.rho_L <= 0 || config_.rho_H <= 0) {
        throw std::invalid_argument("Fluid densities must be positive");
    }
    
    if (config_.mu_L <= 0 || config_.mu_H <= 0) {
        throw std::invalid_argument("Fluid viscosities must be positive");
    }
    
    if (config_.tau_phi <= 0.5) {
        throw std::invalid_argument("Phase field relaxation time must be > 0.5 for stability");
    }
    
    if (config_.xi <= 0 || config_.sigma <= 0) {
        throw std::invalid_argument("Interface thickness and surface tension must be positive");
    }
    
    if (config_.max_timesteps <= 0) {
        throw std::invalid_argument("Maximum timesteps must be positive");
    }
    
    // Calculate derived parameters
    config_.beta = 12.0 * config_.sigma / config_.xi;
    config_.kappa = 3.0 * config_.sigma * config_.xi / 2.0;
    cs2_ = 1.0 / 3.0;
    M_ = config_.tau_phi * cs2_;
    
    // Configure stability monitor
    StabilityConfig stability_config;
    stability_config.max_velocity = config_.max_velocity_limit;
    stability_config.min_density = config_.min_density_limit;
    stability_config.max_density = config_.max_density_ratio * config_.rho_L;
    stability_config.check_interval = config_.stability_check_interval;
    stability_monitor_.updateConfig(stability_config);
    
    std::cout << "MultiphaseSolver initialized:" << std::endl;
    std::cout << "  Grid: " << Nx_ << "x" << Ny_ << std::endl;
    std::cout << "  Density ratio: " << config_.rho_H / config_.rho_L << std::endl;
    std::cout << "  Viscosity ratio: " << config_.mu_H / config_.mu_L << std::endl;
    std::cout << "  Interface thickness: " << config_.xi << std::endl;
    std::cout << "  Surface tension: " << config_.sigma << std::endl;
}

MultiphaseSolver::~MultiphaseSolver() {
    deallocateMemory();
}

void MultiphaseSolver::allocateMemory() {
    if (memory_allocated_) {
        return;  // Already allocated
    }
    
    try {
        // Allocate distribution functions
        h_ = new double**[Q];
        h_temp_ = new double**[Q];
        g_ = new double**[Q];
        g_temp_ = new double**[Q];
        
        for (int i = 0; i < Q; i++) {
            h_[i] = new double*[Nx_];
            h_temp_[i] = new double*[Nx_];
            g_[i] = new double*[Nx_];
            g_temp_[i] = new double*[Nx_];
            
            for (int x = 0; x < Nx_; x++) {
                h_[i][x] = new double[Ny_];
                h_temp_[i][x] = new double[Ny_];
                g_[i][x] = new double[Ny_];
                g_temp_[i][x] = new double[Ny_];
            }
        }
        
        // Allocate macroscopic variables
        phi_ = new double*[Nx_];
        rho_ = new double*[Nx_];
        mu_ = new double*[Nx_];
        tau_ = new double*[Nx_];
        ux_ = new double*[Nx_];
        uy_ = new double*[Nx_];
        p_ = new double*[Nx_];
        
        // Allocate gradient fields
        grad_phi_x_ = new double*[Nx_];
        grad_phi_y_ = new double*[Nx_];
        lapl_phi_ = new double*[Nx_];
        
        for (int x = 0; x < Nx_; x++) {
            phi_[x] = new double[Ny_];
            rho_[x] = new double[Ny_];
            mu_[x] = new double[Ny_];
            tau_[x] = new double[Ny_];
            ux_[x] = new double[Ny_];
            uy_[x] = new double[Ny_];
            p_[x] = new double[Ny_];
            grad_phi_x_[x] = new double[Ny_];
            grad_phi_y_[x] = new double[Ny_];
            lapl_phi_[x] = new double[Ny_];
        }
        
        memory_allocated_ = true;
        
    } catch (const std::bad_alloc& e) {
        deallocateMemory();  // Clean up partial allocation
        throw std::runtime_error("Failed to allocate memory for multiphase solver arrays");
    }
}

void MultiphaseSolver::deallocateMemory() {
    if (!memory_allocated_) {
        return;  // Nothing to deallocate
    }
    
    // Deallocate distribution functions
    if (h_) {
        for (int i = 0; i < Q; i++) {
            if (h_[i]) {
                for (int x = 0; x < Nx_; x++) {
                    delete[] h_[i][x];
                }
                delete[] h_[i];
            }
        }
        delete[] h_;
        h_ = nullptr;
    }
    
    if (h_temp_) {
        for (int i = 0; i < Q; i++) {
            if (h_temp_[i]) {
                for (int x = 0; x < Nx_; x++) {
                    delete[] h_temp_[i][x];
                }
                delete[] h_temp_[i];
            }
        }
        delete[] h_temp_;
        h_temp_ = nullptr;
    }
    
    if (g_) {
        for (int i = 0; i < Q; i++) {
            if (g_[i]) {
                for (int x = 0; x < Nx_; x++) {
                    delete[] g_[i][x];
                }
                delete[] g_[i];
            }
        }
        delete[] g_;
        g_ = nullptr;
    }
    
    if (g_temp_) {
        for (int i = 0; i < Q; i++) {
            if (g_temp_[i]) {
                for (int x = 0; x < Nx_; x++) {
                    delete[] g_temp_[i][x];
                }
                delete[] g_temp_[i];
            }
        }
        delete[] g_temp_;
        g_temp_ = nullptr;
    }
    
    // Deallocate field arrays
    auto deallocate_field = [this](double**& field) {
        if (field) {
            for (int x = 0; x < Nx_; x++) {
                delete[] field[x];
            }
            delete[] field;
            field = nullptr;
        }
    };
    
    deallocate_field(phi_);
    deallocate_field(rho_);
    deallocate_field(mu_);
    deallocate_field(tau_);
    deallocate_field(ux_);
    deallocate_field(uy_);
    deallocate_field(p_);
    deallocate_field(grad_phi_x_);
    deallocate_field(grad_phi_y_);
    deallocate_field(lapl_phi_);
    
    memory_allocated_ = false;
}

void MultiphaseSolver::initialize() {
    allocateMemory();
    initializePhaseField();
    initializeHydrodynamics();
    updateFluidProperties();
    calculateGradients();
    
    // Calculate initial phase masses for conservation checking
    initial_phase_masses_ = calculatePhaseMasses();
    
    // Validate output directory
    if (!utils::IOUtils::validateOutputDirectory(output_directory_)) {
        throw std::runtime_error("Cannot create or access output directory: " + output_directory_);
    }
    
    std::cout << "Multiphase solver initialized successfully" << std::endl;
    std::cout << "Initial light phase mass: " << initial_phase_masses_.first << std::endl;
    std::cout << "Initial heavy phase mass: " << initial_phase_masses_.second << std::endl;
}

void MultiphaseSolver::initializePhaseField() {
    for (int x = 0; x < Nx_; x++) {
        for (int y = 0; y < Ny_; y++) {
            // Create horizontal layers: light fluid below, heavy fluid above
            double y_interface = Ny_ / 2.0;
            phi_[x][y] = config_.phi_0 + (config_.phi_H - config_.phi_L) / 2.0 * 
                        std::tanh(2.0 * (y - y_interface) / config_.xi);
            
            // Clamp to valid range
            phi_[x][y] = utils::MathUtils::clamp(phi_[x][y], config_.phi_L, config_.phi_H);
            
            // Initialize phase field distributions
            for (int i = 0; i < Q; i++) {
                double Gamma_i = w[i] * (1.0 + 3.0 * (ex[i] * 0.0 + ey[i] * 0.0) / cs2_);
                h_[i][x][y] = phi_[x][y] * Gamma_i;
                h_temp_[i][x][y] = 0.0;
            }
        }
    }
}

void MultiphaseSolver::initializeHydrodynamics() {
    for (int x = 0; x < Nx_; x++) {
        for (int y = 0; y < Ny_; y++) {
            ux_[x][y] = 0.0;
            uy_[x][y] = 0.0;
            p_[x][y] = 0.0;
            
            // Initialize velocity distributions
            for (int i = 0; i < Q; i++) {
                double Gamma_i = w[i] * (1.0 + 3.0 * (ex[i] * ux_[x][y] + ey[i] * uy_[x][y]) / cs2_ +
                                        4.5 * std::pow(ex[i] * ux_[x][y] + ey[i] * uy_[x][y], 2) / (cs2_ * cs2_) -
                                        1.5 * (ux_[x][y] * ux_[x][y] + uy_[x][y] * uy_[x][y]) / cs2_);
                g_[i][x][y] = p_[x][y] * w[i] + (Gamma_i - w[i]);
                g_temp_[i][x][y] = 0.0;
            }
        }
    }
}

void MultiphaseSolver::updateFluidProperties() {
    for (int x = 0; x < Nx_; x++) {
        for (int y = 0; y < Ny_; y++) {
            // Clamp phi to valid range to prevent instability
            phi_[x][y] = utils::MathUtils::clamp(phi_[x][y], config_.phi_L, config_.phi_H);
            
            // Linear interpolation for density
            double phi_normalized = (phi_[x][y] - config_.phi_L) / (config_.phi_H - config_.phi_L);
            rho_[x][y] = config_.rho_L + phi_normalized * (config_.rho_H - config_.rho_L);
            
            // Ensure minimum density to prevent division by zero
            rho_[x][y] = std::max(rho_[x][y], config_.min_density_limit);
            
            // Linear interpolation for dynamic viscosity
            mu_[x][y] = config_.mu_L + phi_normalized * (config_.mu_H - config_.mu_L);
            
            // Ensure minimum viscosity
            mu_[x][y] = std::max(mu_[x][y], 0.1 * config_.mu_L);
            
            // Calculate relaxation time with safety bounds
            tau_[x][y] = mu_[x][y] / (rho_[x][y] * cs2_) + 0.5;
            tau_[x][y] = utils::MathUtils::clamp(tau_[x][y], config_.min_tau_limit, config_.max_tau_limit);
        }
    }
}

void MultiphaseSolver::calculateGradients() {
    for (int x = 0; x < Nx_; x++) {
        for (int y = 0; y < Ny_; y++) {
            try {
                auto gradient = utils::MathUtils::calculateGradient(phi_, Nx_, Ny_, x, y);
                grad_phi_x_[x][y] = gradient.first;
                grad_phi_y_[x][y] = gradient.second;
                
                lapl_phi_[x][y] = utils::MathUtils::calculateLaplacian(phi_, Nx_, Ny_, x, y);
                
            } catch (const std::exception& e) {
                // Handle boundary cases gracefully
                grad_phi_x_[x][y] = 0.0;
                grad_phi_y_[x][y] = 0.0;
                lapl_phi_[x][y] = 0.0;
            }
        }
    }
}

double MultiphaseSolver::calculateChemicalPotential(int x, int y) const {
    double phi_val = phi_[x][y];
    double mu_phi = 4.0 * config_.beta * (phi_val - config_.phi_L) * 
                    (phi_val - config_.phi_H) * (phi_val - config_.phi_0) - 
                    config_.kappa * lapl_phi_[x][y];
    return mu_phi;
}

void MultiphaseSolver::step() {
    if (!memory_allocated_) {
        throw std::runtime_error("Solver not initialized - call initialize() first");
    }
    
    // Calculate gradients
    calculateGradients();
    
    // Update fluid properties
    updateFluidProperties();
    
    // Collision steps
    phaseFieldCollision();
    hydrodynamicCollision();
    
    // Streaming
    streaming();
    
    // Boundary conditions
    applyBoundaryConditions();
    
    // Update macroscopic variables
    updateMacroscopicVariables();
    
    // Update timestep
    updateTimestep(getCurrentTimestep() + 1);
    
    // Check stability periodically
    if (getCurrentTimestep() % config_.stability_check_interval == 0) {
        is_stable_ = checkNumericalStability();
        interface_stable_ = checkInterfaceStability();
        phase_conserved_ = checkPhaseConservation();
        last_stability_check_ = getCurrentTimestep();
        
        if (!is_stable_ || !interface_stable_ || !phase_conserved_) {
            std::cerr << "Stability issues detected at timestep " << getCurrentTimestep() << std::endl;
            if (!is_stable_) std::cerr << "  - Numerical instability" << std::endl;
            if (!interface_stable_) std::cerr << "  - Interface instability" << std::endl;
            if (!phase_conserved_) std::cerr << "  - Phase conservation violation" << std::endl;
        }
    }
}

bool MultiphaseSolver::isStable() const {
    return is_stable_ && interface_stable_ && phase_conserved_;
}

void MultiphaseSolver::phaseFieldCollision() {
    for (int x = 0; x < Nx_; x++) {
        for (int y = 0; y < Ny_; y++) {
            // Calculate forcing term
            double grad_phi_mag = std::sqrt(grad_phi_x_[x][y] * grad_phi_x_[x][y] + 
                                           grad_phi_y_[x][y] * grad_phi_y_[x][y]);
            
            double force_factor = 0.0;
            if (grad_phi_mag > 1e-12) {
                force_factor = (1.0 - 4.0 * std::pow(phi_[x][y] - config_.phi_0, 2)) / config_.xi;
            }
            
            for (int i = 0; i < Q; i++) {
                // Equilibrium distribution
                double Gamma_i = w[i] * (1.0 + 3.0 * (ex[i] * ux_[x][y] + ey[i] * uy_[x][y]) / cs2_ +
                                        4.5 * std::pow(ex[i] * ux_[x][y] + ey[i] * uy_[x][y], 2) / (cs2_ * cs2_) -
                                        1.5 * (ux_[x][y] * ux_[x][y] + uy_[x][y] * uy_[x][y]) / cs2_);
                double h_eq = phi_[x][y] * Gamma_i;
                
                // Forcing term
                double F_phi = 0.0;
                if (grad_phi_mag > 1e-12) {
                    F_phi = force_factor * w[i] * (ex[i] * grad_phi_x_[x][y] + ey[i] * grad_phi_y_[x][y]) / grad_phi_mag;
                }
                
                double h_eq_bar = h_eq - 0.5 * F_phi;
                
                // Collision
                h_temp_[i][x][y] = h_[i][x][y] - (h_[i][x][y] - h_eq_bar) / (config_.tau_phi + 0.5) + F_phi;
                
                // Ensure stability
                if (!utils::MathUtils::isFinite(h_temp_[i][x][y])) {
                    h_temp_[i][x][y] = h_eq;
                }
            }
        }
    }
}

void MultiphaseSolver::hydrodynamicCollision() {
    for (int x = 0; x < Nx_; x++) {
        for (int y = 0; y < Ny_; y++) {
            // Surface tension force
            double mu_phi = calculateChemicalPotential(x, y);
            double F_s_x = mu_phi * grad_phi_x_[x][y];
            double F_s_y = mu_phi * grad_phi_y_[x][y];
            
            // Body force
            double F_b_x = rho_[x][y] * config_.gravity;
            double F_b_y = 0.0;
            
            // Total force
            double F_total_x = F_s_x + F_b_x;
            double F_total_y = F_s_y + F_b_y;
            
            for (int i = 0; i < Q; i++) {
                // Equilibrium distribution
                double Gamma_i = w[i] * (1.0 + 3.0 * (ex[i] * ux_[x][y] + ey[i] * uy_[x][y]) / cs2_ +
                                        4.5 * std::pow(ex[i] * ux_[x][y] + ey[i] * uy_[x][y], 2) / (cs2_ * cs2_) -
                                        1.5 * (ux_[x][y] * ux_[x][y] + uy_[x][y] * uy_[x][y]) / cs2_);
                double g_eq = p_[x][y] * w[i] + (Gamma_i - w[i]);
                
                // Forcing term with safety check
                double F_i = 0.0;
                if (rho_[x][y] > config_.min_density_limit) {
                    F_i = w[i] * (ex[i] * F_total_x + ey[i] * F_total_y) / (rho_[x][y] * cs2_);
                }
                
                double g_eq_bar = g_eq - 0.5 * F_i;
                
                // BGK collision with safety check
                if (tau_[x][y] > 0.5) {
                    g_temp_[i][x][y] = g_[i][x][y] - (g_[i][x][y] - g_eq_bar) / tau_[x][y] + F_i;
                } else {
                    g_temp_[i][x][y] = g_[i][x][y]; // No collision if tau is too small
                }
                
                // Ensure stability
                if (!utils::MathUtils::isFinite(g_temp_[i][x][y])) {
                    g_temp_[i][x][y] = g_eq;
                }
            }
        }
    }
}

void MultiphaseSolver::streaming() {
    // Phase field streaming
    static double*** h_new = nullptr;
    static bool first_call_h = true;
    
    if (first_call_h) {
        h_new = new double**[Q];
        for (int i = 0; i < Q; i++) {
            h_new[i] = new double*[Nx_];
            for (int x = 0; x < Nx_; x++) {
                h_new[i][x] = new double[Ny_];
            }
        }
        first_call_h = false;
    }
    
    // Initialize h_new to zero
    for (int i = 0; i < Q; i++) {
        for (int x = 0; x < Nx_; x++) {
            for (int y = 0; y < Ny_; y++) {
                h_new[i][x][y] = 0.0;
            }
        }
    }
    
    // Stream phase field
    for (int x = 0; x < Nx_; x++) {
        for (int y = 0; y < Ny_; y++) {
            for (int i = 0; i < Q; i++) {
                int x_next = (x + ex[i] + Nx_) % Nx_;
                int y_next = y + ey[i];
                
                if (y_next >= 0 && y_next < Ny_) {
                    h_new[i][x_next][y_next] = h_temp_[i][x][y];
                }
            }
        }
    }
    
    // Copy back
    for (int i = 0; i < Q; i++) {
        for (int x = 0; x < Nx_; x++) {
            for (int y = 0; y < Ny_; y++) {
                h_[i][x][y] = h_new[i][x][y];
            }
        }
    }
    
    // Hydrodynamic streaming
    static double*** g_new = nullptr;
    static bool first_call_g = true;
    
    if (first_call_g) {
        g_new = new double**[Q];
        for (int i = 0; i < Q; i++) {
            g_new[i] = new double*[Nx_];
            for (int x = 0; x < Nx_; x++) {
                g_new[i][x] = new double[Ny_];
            }
        }
        first_call_g = false;
    }
    
    // Initialize g_new to zero
    for (int i = 0; i < Q; i++) {
        for (int x = 0; x < Nx_; x++) {
            for (int y = 0; y < Ny_; y++) {
                g_new[i][x][y] = 0.0;
            }
        }
    }
    
    // Stream velocity field
    for (int x = 0; x < Nx_; x++) {
        for (int y = 0; y < Ny_; y++) {
            for (int i = 0; i < Q; i++) {
                int x_next = (x + ex[i] + Nx_) % Nx_;
                int y_next = y + ey[i];
                
                if (y_next >= 0 && y_next < Ny_) {
                    g_new[i][x_next][y_next] = g_temp_[i][x][y];
                }
            }
        }
    }
    
    // Copy back
    for (int i = 0; i < Q; i++) {
        for (int x = 0; x < Nx_; x++) {
            for (int y = 0; y < Ny_; y++) {
                g_[i][x][y] = g_new[i][x][y];
            }
        }
    }
}

void MultiphaseSolver::applyBoundaryConditions() {
    // No-slip walls at top and bottom
    for (int x = 0; x < Nx_; x++) {
        // Bottom wall (y=0)
        for (int i = 0; i < Q; i++) {
            if (ey[i] > 0) {  // Bounce back upward moving distributions
                int i_opp = -1;
                for (int j = 0; j < Q; j++) {
                    if (ex[j] == -ex[i] && ey[j] == -ey[i]) {
                        i_opp = j;
                        break;
                    }
                }
                if (i_opp >= 0) {
                    h_[i_opp][x][0] = h_[i][x][0];
                    g_[i_opp][x][0] = g_[i][x][0];
                }
            }
        }
        
        // Top wall (y=Ny-1)
        for (int i = 0; i < Q; i++) {
            if (ey[i] < 0) {  // Bounce back downward moving distributions
                int i_opp = -1;
                for (int j = 0; j < Q; j++) {
                    if (ex[j] == -ex[i] && ey[j] == -ey[i]) {
                        i_opp = j;
                        break;
                    }
                }
                if (i_opp >= 0) {
                    h_[i_opp][x][Ny_-1] = h_[i][x][Ny_-1];
                    g_[i_opp][x][Ny_-1] = g_[i][x][Ny_-1];
                }
            }
        }
    }
}

void MultiphaseSolver::updateMacroscopicVariables() {
    for (int x = 0; x < Nx_; x++) {
        for (int y = 0; y < Ny_; y++) {
            // Update phase field
            phi_[x][y] = 0.0;
            for (int i = 0; i < Q; i++) {
                phi_[x][y] += h_[i][x][y];
            }
            
            // Clamp phase field to valid range
            phi_[x][y] = utils::MathUtils::clamp(phi_[x][y], config_.phi_L, config_.phi_H);
            
            // Update pressure
            p_[x][y] = 0.0;
            for (int i = 0; i < Q; i++) {
                p_[x][y] += g_[i][x][y];
            }
            
            // Update velocity
            double vel_x = 0.0, vel_y = 0.0;
            for (int i = 0; i < Q; i++) {
                vel_x += g_[i][x][y] * ex[i];
                vel_y += g_[i][x][y] * ey[i];
            }
            
            // Add force correction to velocity with safety check
            double F_total_x = rho_[x][y] * config_.gravity;  // Simplified total force
            double F_total_y = 0.0;
            
            if (rho_[x][y] > config_.min_density_limit) {
                ux_[x][y] = vel_x + 0.5 * F_total_x / rho_[x][y];
                uy_[x][y] = vel_y + 0.5 * F_total_y / rho_[x][y];
            } else {
                ux_[x][y] = 0.0;
                uy_[x][y] = 0.0;
            }
            
            // Clamp velocities to reasonable bounds to prevent instability
            ux_[x][y] = utils::MathUtils::clamp(ux_[x][y], -config_.max_velocity_limit, config_.max_velocity_limit);
            uy_[x][y] = utils::MathUtils::clamp(uy_[x][y], -config_.max_velocity_limit, config_.max_velocity_limit);
        }
    }
}

bool MultiphaseSolver::checkNumericalStability() {
    // Create lattice data structure for stability checking
    LatticeData lattice_data(Nx_, Ny_);
    
    // Copy data to lattice structure  
    for (int x = 0; x < Nx_; ++x) {
        for (int y = 0; y < Ny_; ++y) {
            lattice_data.density[x][y] = rho_[x][y];
            lattice_data.velocity_x[x][y] = ux_[x][y];
            lattice_data.velocity_y[x][y] = uy_[x][y];
            
            for (int i = 0; i < Q; ++i) {
                lattice_data.distribution[x][y][i] = h_[i][x][y];
            }
        }
    }
    
    // Check stability using monitor
    StabilityMetrics metrics = stability_monitor_.checkStability(lattice_data, getCurrentTimestep());
    
    return metrics.isStable();
}

bool MultiphaseSolver::checkInterfaceStability() {
    // Check interface thickness and position stability
    double interface_thickness_sum = 0.0;
    int interface_count = 0;
    
    for (int x = 0; x < Nx_; x++) {
        double interface_y = getInterfacePosition(x);
        if (interface_y >= 0 && interface_y < Ny_) {
            // Calculate local interface thickness
            double max_gradient = 0.0;
            for (int y = std::max(0, static_cast<int>(interface_y) - 5); 
                 y < std::min(Ny_, static_cast<int>(interface_y) + 5); y++) {
                double grad_mag = std::sqrt(grad_phi_x_[x][y] * grad_phi_x_[x][y] + 
                                           grad_phi_y_[x][y] * grad_phi_y_[x][y]);
                max_gradient = std::max(max_gradient, grad_mag);
            }
            
            if (max_gradient > 1e-12) {
                double thickness = (config_.phi_H - config_.phi_L) / max_gradient;
                interface_thickness_sum += thickness;
                interface_count++;
            }
        }
    }
    
    if (interface_count > 0) {
        double avg_thickness = interface_thickness_sum / interface_count;
        double thickness_deviation = std::abs(avg_thickness - config_.xi) / config_.xi;
        
        if (thickness_deviation > config_.interface_thickness_tolerance) {
            return false;
        }
    }
    
    return true;
}

bool MultiphaseSolver::checkPhaseConservation() {
    auto current_masses = calculatePhaseMasses();
    
    double light_mass_change = std::abs(current_masses.first - initial_phase_masses_.first) / 
                              initial_phase_masses_.first;
    double heavy_mass_change = std::abs(current_masses.second - initial_phase_masses_.second) / 
                              initial_phase_masses_.second;
    
    return (light_mass_change < config_.phase_conservation_tolerance) && 
           (heavy_mass_change < config_.phase_conservation_tolerance);
}

double MultiphaseSolver::getPhaseField(int x, int y) const {
    if (x < 0 || x >= Nx_ || y < 0 || y >= Ny_ || !memory_allocated_) {
        return 0.0;
    }
    return phi_[x][y];
}

double MultiphaseSolver::getDensity(int x, int y) const {
    if (x < 0 || x >= Nx_ || y < 0 || y >= Ny_ || !memory_allocated_) {
        return 0.0;
    }
    return rho_[x][y];
}

std::pair<double, double> MultiphaseSolver::getVelocity(int x, int y) const {
    if (x < 0 || x >= Nx_ || y < 0 || y >= Ny_ || !memory_allocated_) {
        return {0.0, 0.0};
    }
    return {ux_[x][y], uy_[x][y]};
}

double MultiphaseSolver::getPressure(int x, int y) const {
    if (x < 0 || x >= Nx_ || y < 0 || y >= Ny_ || !memory_allocated_) {
        return 0.0;
    }
    return p_[x][y];
}

double MultiphaseSolver::getAnalyticalVelocity(int y) const {
    // Analytical solution for two-layer Poiseuille flow
    double H = Ny_ - 1;
    double y_interface = H / 2.0;
    
    if (y <= y_interface) {
        // Light fluid (lower layer)
        double C1 = config_.gravity * config_.rho_L * (H - y_interface) / (2.0 * config_.mu_L);
        double C2 = config_.gravity * config_.rho_H * (H - y_interface) * y_interface / (2.0 * config_.mu_H) - 
                    config_.gravity * config_.rho_L * y_interface * y_interface / (2.0 * config_.mu_L);
        return C1 * y + C2;
    } else {
        // Heavy fluid (upper layer)
        double C3 = config_.gravity * config_.rho_H / (2.0 * config_.mu_H);
        double C4 = -config_.gravity * config_.rho_H * H * H / (2.0 * config_.mu_H) + 
                    config_.gravity * config_.rho_L * (H - y_interface) * y_interface / (2.0 * config_.mu_L) + 
                    config_.gravity * config_.rho_H * (H - y_interface) * y_interface / (2.0 * config_.mu_H);
        return C3 * (H * H - y * y) + C4;
    }
}

double MultiphaseSolver::getInterfacePosition(int x) const {
    if (x < 0 || x >= Nx_ || !memory_allocated_) {
        return -1.0;
    }
    
    // Find y-position where phi is closest to phi_0
    double min_diff = std::numeric_limits<double>::max();
    int interface_y = -1;
    
    for (int y = 0; y < Ny_; y++) {
        double diff = std::abs(phi_[x][y] - config_.phi_0);
        if (diff < min_diff) {
            min_diff = diff;
            interface_y = y;
        }
    }
    
    return static_cast<double>(interface_y);
}

std::pair<double, double> MultiphaseSolver::calculatePhaseMasses() const {
    if (!memory_allocated_) {
        return {0.0, 0.0};
    }
    
    double light_mass = 0.0;
    double heavy_mass = 0.0;
    
    for (int x = 0; x < Nx_; x++) {
        for (int y = 0; y < Ny_; y++) {
            double phi_normalized = (phi_[x][y] - config_.phi_L) / (config_.phi_H - config_.phi_L);
            light_mass += (1.0 - phi_normalized);
            heavy_mass += phi_normalized;
        }
    }
    
    return {light_mass, heavy_mass};
}

void MultiphaseSolver::writeOutput(int timestep) {
    if (config_.write_full_fields) {
        writeFullFieldData(timestep);
    }
    if (config_.write_averaged_profiles) {
        writeAveragedProfiles(timestep);
    }
}

void MultiphaseSolver::writeFullFieldData(int timestep) {
    std::string filename = output_directory_ + "/" + 
                          utils::IOUtils::generateTimestepFilename(config_.output_prefix, timestep);
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Warning: Could not open output file: " << filename << std::endl;
        return;
    }
    
    file << "x,y,phi,rho,ux,uy,p";
    if (config_.write_analytical_comparison) {
        file << ",analytical_ux";
    }
    file << "\n";
    
    for (int y = 0; y < Ny_; y++) {
        for (int x = 0; x < Nx_; x++) {
            file << x << "," << y << "," 
                 << utils::IOUtils::formatDouble(phi_[x][y]) << ","
                 << utils::IOUtils::formatDouble(rho_[x][y]) << ","
                 << utils::IOUtils::formatDouble(ux_[x][y]) << ","
                 << utils::IOUtils::formatDouble(uy_[x][y]) << ","
                 << utils::IOUtils::formatDouble(p_[x][y]);
            
            if (config_.write_analytical_comparison) {
                file << "," << utils::IOUtils::formatDouble(getAnalyticalVelocity(y));
            }
            file << "\n";
        }
    }
    
    file.close();
}

void MultiphaseSolver::writeAveragedProfiles(int timestep) {
    std::string filename = output_directory_ + "/" + 
                          utils::IOUtils::generateTimestepFilename(config_.output_prefix + "_avg", timestep);
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Warning: Could not open averaged output file: " << filename << std::endl;
        return;
    }
    
    file << "y,avg_phi,avg_rho,avg_ux";
    if (config_.write_analytical_comparison) {
        file << ",analytical_ux";
    }
    file << "\n";
    
    for (int y = 0; y < Ny_; y++) {
        double avg_phi = 0.0, avg_rho = 0.0, avg_ux = 0.0;
        
        for (int x = 0; x < Nx_; x++) {
            avg_phi += phi_[x][y];
            avg_rho += rho_[x][y];
            avg_ux += ux_[x][y];
        }
        
        avg_phi /= Nx_;
        avg_rho /= Nx_;
        avg_ux /= Nx_;
        
        // Validate data before output
        if (!utils::MathUtils::isFinite(avg_phi)) avg_phi = config_.phi_0;
        if (!utils::MathUtils::isFinite(avg_rho)) avg_rho = config_.rho_L;
        if (!utils::MathUtils::isFinite(avg_ux)) avg_ux = 0.0;
        
        file << y << "," 
             << utils::IOUtils::formatDouble(avg_phi) << ","
             << utils::IOUtils::formatDouble(avg_rho) << ","
             << utils::IOUtils::formatDouble(avg_ux);
        
        if (config_.write_analytical_comparison) {
            file << "," << utils::IOUtils::formatDouble(getAnalyticalVelocity(y));
        }
        file << "\n";
    }
    
    file.close();
}

bool MultiphaseSolver::runSimulation() {
    try {
        initialize();
        
        std::cout << "Starting multiphase simulation with " << config_.max_timesteps << " timesteps..." << std::endl;
        std::cout << "Interface thickness: " << config_.xi << std::endl;
        std::cout << "Surface tension: " << config_.sigma << std::endl;
        
        for (int t = 0; t <= config_.max_timesteps; t++) {
            step();
            
            if (!isStable()) {
                std::cerr << "Simulation terminated due to instability at timestep " << t << std::endl;
                return false;
            }
            
            if (t % config_.output_interval == 0) {
                auto current_masses = calculatePhaseMasses();
                auto history = stability_monitor_.getStabilityHistory();
                double max_vel = history.empty() ? 0.0 : history.back().max_velocity;
                utils::IOUtils::printProgress(t, config_.max_timesteps, 
                                            "Max velocity: " + std::to_string(max_vel) +
                                            ", Light mass: " + std::to_string(current_masses.first));
                writeOutput(t);
            }
        }
        
        std::cout << "Multiphase simulation completed successfully!" << std::endl;
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Multiphase simulation failed with error: " << e.what() << std::endl;
        return false;
    }
}

} // namespace lbm