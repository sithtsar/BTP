#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <vector>
#include <algorithm>

// Grid parameters
const int Nx = 256;
const int Ny = 64;
const int Q = 9;

// Physical parameters
const double cs2 = 1.0 / 3.0;
const double cs = std::sqrt(cs2);

// Phase field parameters
const double phi_L = 0.0;    // Light fluid phase field value
const double phi_H = 1.0;    // Heavy fluid phase field value  
const double phi_0 = 0.5;    // Interface location

// Fluid properties
const double rho_L = 1.0;    // Light fluid density
const double rho_H = 1000.0; // Heavy fluid density (high density ratio)
const double mu_L = 0.01;    // Light fluid viscosity
const double mu_H = 1.0;     // Heavy fluid viscosity

// Surface tension parameters
const double xi = 4.0;       // Interface thickness
const double sigma = 0.01;   // Surface tension coefficient
const double beta = 12.0 * sigma / xi;
const double kappa = 3.0 * sigma * xi / 2.0;

// Simulation parameters
const double tau_phi = 0.7;  // Phase field relaxation time
const double M = tau_phi * cs2; // Mobility
const double gravity = 1e-6;       // Gravitational acceleration
const int maxT = 20000;
const int output_interval = 1000;

// LBM parameters
const int ex[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int ey[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
const double w[Q] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 
                     1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

// Distribution functions
double h[Q][Nx][Ny];      // Phase field distribution
double h_temp[Q][Nx][Ny];
double g[Q][Nx][Ny];      // Velocity-based distribution  
double g_temp[Q][Nx][Ny];

// Macroscopic variables
double phi[Nx][Ny];       // Phase field
double rho[Nx][Ny];       // Density
double mu[Nx][Ny];        // Dynamic viscosity
double tau[Nx][Ny];       // Relaxation time
double ux[Nx][Ny];        // x-velocity
double uy[Nx][Ny];        // y-velocity
double p[Nx][Ny];         // Pressure

// Gradient arrays
double grad_phi_x[Nx][Ny];
double grad_phi_y[Nx][Ny];
double lapl_phi[Nx][Ny];

// Initialize phase field with a layered configuration
void initialize_phase_field() {
    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            // Create horizontal layers: light fluid below, heavy fluid above
            double y_interface = Ny / 2.0;
            phi[x][y] = phi_0 + (phi_H - phi_L) / 2.0 * 
                       std::tanh(2.0 * (y - y_interface) / xi);
            
            // Initialize phase field distributions
            for (int i = 0; i < Q; i++) {
                double Gamma_i = w[i] * (1.0 + 3.0 * (ex[i] * 0.0 + ey[i] * 0.0) / cs2);
                h[i][x][y] = phi[x][y] * Gamma_i;
            }
        }
    }
}

// Initialize hydrodynamic fields
void initialize_hydrodynamics() {
    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            ux[x][y] = 0.0;
            uy[x][y] = 0.0;
            p[x][y] = 0.0;
            
            // Initialize velocity distributions
            for (int i = 0; i < Q; i++) {
                double Gamma_i = w[i] * (1.0 + 3.0 * (ex[i] * ux[x][y] + ey[i] * uy[x][y]) / cs2 +
                                        4.5 * std::pow(ex[i] * ux[x][y] + ey[i] * uy[x][y], 2) / (cs2 * cs2) -
                                        1.5 * (ux[x][y] * ux[x][y] + uy[x][y] * uy[x][y]) / cs2);
                g[i][x][y] = p[x][y] * w[i] + (Gamma_i - w[i]);
            }
        }
    }
}

// Update density and viscosity from phase field
void update_properties() {
    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            // Clamp phi to valid range to prevent instability
            phi[x][y] = std::max(phi_L, std::min(phi_H, phi[x][y]));
            
            // Linear interpolation for density
            rho[x][y] = rho_L + (phi[x][y] - phi_L) * (rho_H - rho_L) / (phi_H - phi_L);
            
            // Ensure minimum density to prevent division by zero
            rho[x][y] = std::max(rho[x][y], 0.1 * rho_L);
            
            // Linear interpolation for dynamic viscosity
            mu[x][y] = mu_L + (phi[x][y] - phi_L) * (mu_H - mu_L) / (phi_H - phi_L);
            
            // Ensure minimum viscosity
            mu[x][y] = std::max(mu[x][y], 0.1 * mu_L);
            
            // Calculate relaxation time with safety bounds
            tau[x][y] = mu[x][y] / (rho[x][y] * cs2) + 0.5;
            tau[x][y] = std::max(0.51, std::min(5.0, tau[x][y])); // Stability bounds
        }
    }
}

// Calculate gradients using isotropic finite differences
void calculate_gradients() {
    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            grad_phi_x[x][y] = 0.0;
            grad_phi_y[x][y] = 0.0;
            lapl_phi[x][y] = 0.0;
            
            for (int i = 1; i < Q; i++) {  // Skip i=0 (rest particle)
                int x_nb = (x + ex[i] + Nx) % Nx;
                int y_nb = y + ey[i];
                
                if (y_nb >= 0 && y_nb < Ny) {
                    // Gradient calculation
                    grad_phi_x[x][y] += w[i] * ex[i] * phi[x_nb][y_nb] / cs2;
                    grad_phi_y[x][y] += w[i] * ey[i] * phi[x_nb][y_nb] / cs2;
                    
                    // Laplacian calculation
                    lapl_phi[x][y] += 2.0 * w[i] * (phi[x_nb][y_nb] - phi[x][y]) / cs2;
                }
            }
        }
    }
}

// Calculate chemical potential
double chemical_potential(int x, int y) {
    double mu_phi = 4.0 * beta * (phi[x][y] - phi_L) * (phi[x][y] - phi_H) * (phi[x][y] - phi_0) 
                    - kappa * lapl_phi[x][y];
    return mu_phi;
}

// Phase field collision step
void phase_field_collision() {
    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            // Calculate forcing term
            double grad_phi_mag = std::sqrt(grad_phi_x[x][y] * grad_phi_x[x][y] + 
                                           grad_phi_y[x][y] * grad_phi_y[x][y]);
            
            double force_factor = 0.0;
            if (grad_phi_mag > 1e-12) {
                force_factor = (1.0 - 4.0 * std::pow(phi[x][y] - phi_0, 2)) / xi;
            }
            
            for (int i = 0; i < Q; i++) {
                // Equilibrium distribution
                double Gamma_i = w[i] * (1.0 + 3.0 * (ex[i] * ux[x][y] + ey[i] * uy[x][y]) / cs2 +
                                        4.5 * std::pow(ex[i] * ux[x][y] + ey[i] * uy[x][y], 2) / (cs2 * cs2) -
                                        1.5 * (ux[x][y] * ux[x][y] + uy[x][y] * uy[x][y]) / cs2);
                double h_eq = phi[x][y] * Gamma_i;
                
                // Forcing term
                double F_phi = 0.0;
                if (grad_phi_mag > 1e-12) {
                    F_phi = force_factor * w[i] * (ex[i] * grad_phi_x[x][y] + ey[i] * grad_phi_y[x][y]) / grad_phi_mag;
                }
                
                double h_eq_bar = h_eq - 0.5 * F_phi;
                
                // Collision
                h_temp[i][x][y] = h[i][x][y] - (h[i][x][y] - h_eq_bar) / (tau_phi + 0.5) + F_phi;
            }
        }
    }
}

// Hydrodynamic collision step
void hydrodynamic_collision() {
    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            // Surface tension force
            double mu_phi = chemical_potential(x, y);
            double F_s_x = mu_phi * grad_phi_x[x][y];
            double F_s_y = mu_phi * grad_phi_y[x][y];
            
            // Body force
            double F_b_x = rho[x][y] * gravity;
            double F_b_y = 0.0;
            
            // Pressure force (simplified - can be enhanced)
            double F_p_x = 0.0; 
            double F_p_y = 0.0;
            
            // Total force
            double F_total_x = F_s_x + F_b_x + F_p_x;
            double F_total_y = F_s_y + F_b_y + F_p_y;
            
            for (int i = 0; i < Q; i++) {
                // Equilibrium distribution  
                double Gamma_i = w[i] * (1.0 + 3.0 * (ex[i] * ux[x][y] + ey[i] * uy[x][y]) / cs2 +
                                        4.5 * std::pow(ex[i] * ux[x][y] + ey[i] * uy[x][y], 2) / (cs2 * cs2) -
                                        1.5 * (ux[x][y] * ux[x][y] + uy[x][y] * uy[x][y]) / cs2);
                double g_eq = p[x][y] * w[i] + (Gamma_i - w[i]);
                
                // Forcing term with safety check
                double F_i = 0.0;
                if (rho[x][y] > 1e-10) {
                    F_i = w[i] * (ex[i] * F_total_x + ey[i] * F_total_y) / (rho[x][y] * cs2);
                }
                
                double g_eq_bar = g_eq - 0.5 * F_i;
                
                // BGK collision with safety check
                if (tau[x][y] > 0.5) {
                    g_temp[i][x][y] = g[i][x][y] - (g[i][x][y] - g_eq_bar) / tau[x][y] + F_i;
                } else {
                    g_temp[i][x][y] = g[i][x][y]; // No collision if tau is too small
                }
            }
        }
    }
}

// Streaming step
void streaming() {
    // Phase field streaming
    static double h_new[Q][Nx][Ny];
    for (int i = 0; i < Q; i++) {
        for (int x = 0; x < Nx; x++) {
            for (int y = 0; y < Ny; y++) {
                h_new[i][x][y] = 0.0;
            }
        }
    }
    
    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            for (int i = 0; i < Q; i++) {
                int x_next = (x + ex[i] + Nx) % Nx;
                int y_next = y + ey[i];
                
                if (y_next >= 0 && y_next < Ny) {
                    h_new[i][x_next][y_next] = h_temp[i][x][y];
                }
            }
        }
    }
    
    // Copy back
    for (int i = 0; i < Q; i++) {
        for (int x = 0; x < Nx; x++) {
            for (int y = 0; y < Ny; y++) {
                h[i][x][y] = h_new[i][x][y];
            }
        }
    }
    
    // Hydrodynamic streaming
    static double g_new[Q][Nx][Ny];
    for (int i = 0; i < Q; i++) {
        for (int x = 0; x < Nx; x++) {
            for (int y = 0; y < Ny; y++) {
                g_new[i][x][y] = 0.0;
            }
        }
    }
    
    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            for (int i = 0; i < Q; i++) {
                int x_next = (x + ex[i] + Nx) % Nx;
                int y_next = y + ey[i];
                
                if (y_next >= 0 && y_next < Ny) {
                    g_new[i][x_next][y_next] = g_temp[i][x][y];
                }
            }
        }
    }
    
    // Copy back
    for (int i = 0; i < Q; i++) {
        for (int x = 0; x < Nx; x++) {
            for (int y = 0; y < Ny; y++) {
                g[i][x][y] = g_new[i][x][y];
            }
        }
    }
}

// Apply boundary conditions
void boundary_conditions() {
    // No-slip walls at top and bottom
    for (int x = 0; x < Nx; x++) {
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
                    h[i_opp][x][0] = h[i][x][0];
                    g[i_opp][x][0] = g[i][x][0];
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
                    h[i_opp][x][Ny-1] = h[i][x][Ny-1];
                    g[i_opp][x][Ny-1] = g[i][x][Ny-1];
                }
            }
        }
    }
}

// Update macroscopic variables
void update_macroscopic() {
    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            // Update phase field
            phi[x][y] = 0.0;
            for (int i = 0; i < Q; i++) {
                phi[x][y] += h[i][x][y];
            }
            
            // Update pressure
            p[x][y] = 0.0;
            for (int i = 0; i < Q; i++) {
                p[x][y] += g[i][x][y];
            }
            
            // Update velocity
            double vel_x = 0.0, vel_y = 0.0;
            for (int i = 0; i < Q; i++) {
                vel_x += g[i][x][y] * ex[i];
                vel_y += g[i][x][y] * ey[i];
            }
            
            // Add force correction to velocity with safety check
            double F_total_x = rho[x][y] * gravity;  // Simplified total force
            double F_total_y = 0.0;
            
            if (rho[x][y] > 1e-10) {
                ux[x][y] = vel_x + 0.5 * F_total_x / rho[x][y];
                uy[x][y] = vel_y + 0.5 * F_total_y / rho[x][y];
            } else {
                ux[x][y] = 0.0;
                uy[x][y] = 0.0;
            }
            
            // Clamp velocities to reasonable bounds to prevent instability
            ux[x][y] = std::max(-1.0, std::min(1.0, ux[x][y]));
            uy[x][y] = std::max(-1.0, std::min(1.0, uy[x][y]));
        }
    }
}

// Analytical solution for two-layer Poiseuille flow
double analytical_velocity_twolayer(int y) {
    double H = Ny - 1;
    double y_interface = H / 2.0;
    
    if (y <= y_interface) {
        // Light fluid (lower layer)
        double C1 = gravity * rho_L * (H - y_interface) / (2.0 * mu_L);
        double C2 = gravity * rho_H * (H - y_interface) * y_interface / (2.0 * mu_H) - 
                    gravity * rho_L * y_interface * y_interface / (2.0 * mu_L);
        return C1 * y + C2;
    } else {
        // Heavy fluid (upper layer)  
        double C3 = gravity * rho_H / (2.0 * mu_H);
        double C4 = -gravity * rho_H * H * H / (2.0 * mu_H) + 
                    gravity * rho_L * (H - y_interface) * y_interface / (2.0 * mu_L) + 
                    gravity * rho_H * (H - y_interface) * y_interface / (2.0 * mu_H);
        return C3 * (H * H - y * y) + C4;
    }
}

// Output data
void output_data(int t) {
    std::string filename = "data/multiphase_t" + std::to_string(t) + ".csv";
    std::ofstream file(filename);
    
    file << "x,y,phi,rho,ux,uy,p,analytical_ux\n";
    
    for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
            double analytical = analytical_velocity_twolayer(y);
            file << x << "," << y << "," << phi[x][y] << "," << rho[x][y] << ","
                 << ux[x][y] << "," << uy[x][y] << "," << p[x][y] << "," 
                 << analytical << "\n";
        }
    }
    
    file.close();
    
    // Output average profiles
    std::string avg_filename = "data/multiphase_avg_t" + std::to_string(t) + ".csv";
    std::ofstream avg_file(avg_filename);
    
    avg_file << "y,avg_phi,avg_rho,avg_ux,analytical_ux\n";
    
    for (int y = 0; y < Ny; y++) {
        double avg_phi_y = 0.0, avg_rho_y = 0.0, avg_ux_y = 0.0;
        
        for (int x = 0; x < Nx; x++) {
            avg_phi_y += phi[x][y];
            avg_rho_y += rho[x][y];
            avg_ux_y += ux[x][y];
        }
        
        avg_phi_y /= Nx;
        avg_rho_y /= Nx;
        avg_ux_y /= Nx;
        
        // Validate data before output
        if (!std::isfinite(avg_phi_y)) avg_phi_y = 0.0;
        if (!std::isfinite(avg_rho_y)) avg_rho_y = rho_L;
        if (!std::isfinite(avg_ux_y)) avg_ux_y = 0.0;
        
        double analytical = analytical_velocity_twolayer(y);
        
        avg_file << y << "," << avg_phi_y << "," << avg_rho_y << ","
                 << avg_ux_y << "," << analytical << "\n";
    }
    
    avg_file.close();
}

int main() {
    std::cout << "Multiphase LBM Simulation Starting..." << std::endl;
    std::cout << "Grid: " << Nx << "x" << Ny << std::endl;
    std::cout << "Density ratio: " << rho_H/rho_L << std::endl;
    std::cout << "Viscosity ratio: " << mu_H/mu_L << std::endl;
    std::cout << "Interface thickness: " << xi << std::endl;
    
    // Initialize
    initialize_phase_field();
    initialize_hydrodynamics();
    update_properties();
    
    // Main time loop
    for (int t = 0; t <= maxT; t++) {
        // Calculate gradients
        calculate_gradients();
        
        // Update fluid properties
        update_properties();
        
        // Collision steps
        phase_field_collision();
        hydrodynamic_collision();
        
        // Streaming
        streaming();
        
        // Boundary conditions
        boundary_conditions();
        
        // Update macroscopic variables
        update_macroscopic();
        
        // Output
        if (t % output_interval == 0) {
            std::cout << "Time step: " << t << std::endl;
            output_data(t);
        }
    }
    
    std::cout << "Simulation completed!" << std::endl;
    return 0;
}