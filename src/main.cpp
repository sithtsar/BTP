#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>

const int Nx = 20;
const int Ny = 21;
const int Q = 9;
const double tau = 1.0;
const double omega = 1.0 / tau;
const double g = 1e-6;
const int maxT = 10000;
const int output_interval = 1000;
const double rho0 = 1.0;
const double cs2 = 1.0 / 3.0;
const double nu = cs2 * (tau - 0.5);

const int ex[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int ey[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
const double w[Q] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 
                     1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

double f[Q][Nx][Ny];
double f_temp[Q][Nx][Ny];

void initialize() {
    for (int i = 0; i < Q; i++) {
        for (int x = 0; x < Nx; x++) {
            for (int y = 0; y < Ny; y++) {
                f[i][x][y] = w[i] * rho0;
            }
        }
    }
}

double analytical_velocity(int y) {
    double H = Ny - 1;
    return g * H * H / (8.0 * nu) * (1.0 - 4.0 * (y - H/2.0) * (y - H/2.0) / (H * H));
}

void collision(bool entropic) {
    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            double rho = 0.0;
            for (int i = 0; i < Q; i++) {
                rho += f[i][x][y];
            }
            
            // Prevent division by zero
            if (rho < 1e-10) {
                rho = rho0;
            }
            
            double ux = 0.0, uy = 0.0;
            for (int i = 0; i < Q; i++) {
                ux += f[i][x][y] * ex[i];
                uy += f[i][x][y] * ey[i];
            }
            ux /= rho;
            uy /= rho;
            
            // Add body force to velocity
            ux += 0.5 * g / rho;
            uy += 0.0;
            
            double u2 = ux * ux + uy * uy;
            
            double feq[Q];
            for (int i = 0; i < Q; i++) {
                double eu = ex[i] * ux + ey[i] * uy;
                feq[i] = w[i] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u2);
            }
            
            // Improved force term implementation
            double force_term[Q];
            for (int i = 0; i < Q; i++) {
                double eu = ex[i] * ux + ey[i] * uy;
                double F_dot_e = g * ex[i];
                force_term[i] = w[i] * (1.0 - 0.5 * omega) * 
                               (3.0 * F_dot_e + 9.0 * F_dot_e * eu - 3.0 * g * ux);
            }
            
            if (!entropic) {
                for (int i = 0; i < Q; i++) {
                    f_temp[i][x][y] = f[i][x][y] - omega * (f[i][x][y] - feq[i]) + force_term[i];
                }
            } else {
                // Proper Entropic BGK (ELBGK) implementation
                double alpha = 1.0;
                
                // Calculate current entropy H = sum_i f_i * ln(f_i/w_i)
                double H_current = 0.0;
                for (int i = 0; i < Q; i++) {
                    if (f[i][x][y] > 1e-12) {
                        H_current += f[i][x][y] * std::log(f[i][x][y] / w[i]);
                    }
                }
                
                // Find alpha such that entropy is conserved during over-relaxation
                // We need to solve: H(f + alpha*(feq - f)) = H(f)
                // Use bisection method to find alpha
                double alpha_min = 0.0, alpha_max = 2.0;
                int max_iter = 20;
                double tolerance = 1e-8;
                
                for (int iter = 0; iter < max_iter; iter++) {
                    alpha = 0.5 * (alpha_min + alpha_max);
                    
                    // Calculate H at intermediate state
                    double H_test = 0.0;
                    for (int i = 0; i < Q; i++) {
                        double f_test = f[i][x][y] + alpha * (feq[i] - f[i][x][y]);
                        if (f_test > 1e-12) {
                            H_test += f_test * std::log(f_test / w[i]);
                        }
                    }
                    
                    double diff = H_test - H_current;
                    if (std::abs(diff) < tolerance) break;
                    
                    if (diff > 0) {
                        alpha_max = alpha;
                    } else {
                        alpha_min = alpha;
                    }
                }
                
                // Apply entropic collision with found alpha
                for (int i = 0; i < Q; i++) {
                    f_temp[i][x][y] = f[i][x][y] + alpha * omega * (feq[i] - f[i][x][y]) + force_term[i];
                }
            }
        }
    }
}

void streaming() {
    static double f_new[Q][Nx][Ny];
    
    // Initialize f_new to zero
    for (int i = 0; i < Q; i++) {
        for (int x = 0; x < Nx; x++) {
            for (int y = 0; y < Ny; y++) {
                f_new[i][x][y] = 0.0;
            }
        }
    }
    
    // Stream from all interior and boundary points
    for (int x = 0; x < Nx; x++) {
        for (int y = 0; y < Ny; y++) {
            for (int i = 0; i < Q; i++) {
                int x_next = (x + ex[i] + Nx) % Nx;
                int y_next = y + ey[i];
                
                if (y_next >= 0 && y_next < Ny) {
                    f_new[i][x_next][y_next] = f_temp[i][x][y];
                }
            }
        }
    }
    
    // Copy back to f array
    for (int i = 0; i < Q; i++) {
        for (int x = 0; x < Nx; x++) {
            for (int y = 0; y < Ny; y++) {
                f[i][x][y] = f_new[i][x][y];
            }
        }
    }
}

void boundary_conditions() {
    // Bottom wall (y=0) - bounce back
    for (int x = 0; x < Nx; x++) {
        f[2][x][0] = f[4][x][0];  // north <- south
        f[5][x][0] = f[7][x][0];  // NE <- SW
        f[6][x][0] = f[8][x][0];  // NW <- SE
    }
    
    // Top wall (y=Ny-1) - bounce back  
    for (int x = 0; x < Nx; x++) {
        f[4][x][Ny-1] = f[2][x][Ny-1];  // south <- north
        f[7][x][Ny-1] = f[5][x][Ny-1];  // SW <- NE
        f[8][x][Ny-1] = f[6][x][Ny-1];  // SE <- NW
    }
}

void output_data(int t, const std::string& mode) {
    std::ofstream vel_file("data/" + mode + "_velocity_t" + std::to_string(t) + ".csv");
    std::ofstream h_file("data/" + mode + "_h_y_t" + std::to_string(t) + ".csv");
    
    vel_file << "y,avg_ux,analytic_ux\n";
    h_file << "y,H_y\n";
    
    for (int y = 0; y < Ny; y++) {
        double avg_ux = 0.0;
        double H_y = 0.0;
        
        for (int x = 0; x < Nx; x++) {
            double rho = 0.0;
            double ux = 0.0;
            
            for (int i = 0; i < Q; i++) {
                rho += f[i][x][y];
                ux += f[i][x][y] * ex[i];
            }
            
            // Prevent division by zero
            if (rho > 1e-10) {
                ux /= rho;
            } else {
                ux = 0.0;
            }
            avg_ux += ux;
            
            // Calculate H-function: H = sum_i f_i * ln(f_i/w_i)
            // Entropy S = -H should be maximized (H should be minimized)
            for (int i = 0; i < Q; i++) {
                if (f[i][x][y] > 1e-12) {
                    H_y += f[i][x][y] * std::log(f[i][x][y] / w[i]);
                }
            }
        }
        avg_ux /= Nx;
        
        double analytic = analytical_velocity(y);
        vel_file << y << "," << avg_ux << "," << analytic << "\n";
        h_file << y << "," << H_y << "\n";
    }
    
    vel_file.close();
    h_file.close();
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: ./lbm_sim [standard|entropic]" << std::endl;
        return 1;
    }
    
    std::string mode(argv[1]);
    bool entropic = (mode == "entropic");
    
    std::cout << "Running LBM simulation in " << mode << " mode..." << std::endl;
    std::cout << "Grid: " << Nx << "x" << Ny << ", tau=" << tau << ", g=" << g << std::endl;
    std::cout << "Analytical nu=" << nu << std::endl;
    std::cout << "Expected max velocity: " << g * (Ny-1) * (Ny-1) / (8.0 * nu) << std::endl;
    
    initialize();
    
    for (int t = 0; t <= maxT; t++) {
        collision(entropic);
        streaming();
        boundary_conditions();
        
        if (t % output_interval == 0) {
            std::cout << "Time step: " << t << std::endl;
            output_data(t, mode);
        }
    }
    
    std::cout << "Simulation completed. Data saved to data/ directory." << std::endl;
    return 0;
}