/*
 * Entropic Lattice Boltzmann Model (ELBM) with Color-Gradient for Multiphase Flow
 * Case: Static Droplet in 2D
 *
 * Features:
 * - D2Q9 Lattice
 * - Entropic Stabilizer (Newton-Raphson root finding for alpha)
 * - Latva-Kokko Color Gradient (Recoloring + Perturbation)
 * - Exact H-Theorem validity check
 *
 * References:
 *  Entropic Stabilizers
 *  Color Gradient Method
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <iomanip>
#include <algorithm>

// --- Simulation Constants ---
const int NX = 64;
const int NY = 64;
const int MAX_ITER = 2000;
const int OUTPUT_FREQ = 100;

// Physical Parameters
const double REYNOLDS = 100.0;
const double VISCOSITY = 0.05; // Base viscosity
const double TAU_BASE = 3.0 * VISCOSITY + 0.5;
const double SURFACE_TENSION = 0.01;
const double BETA_SEG = 0.7;   // Segregation parameter
const double DENSITY_RED = 1.0;
const double DENSITY_BLUE = 1.0;

// Lattice D2Q9 Constants
const int Q = 9;
const int cx[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int cy[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
const double w[Q] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,
                     1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

// --- Data Structures ---
struct Node {
    double f_red[Q];
    double f_blue[Q];
    double f_new_red[Q];
    double f_new_blue[Q];
    double rho;
    double rho_red;
    double rho_blue;
    double ux;
    double uy;
    double phi; // Phase field: (rho_r - rho_b) / (rho_r + rho_b)
};

std::vector<Node> grid(NX * NY);

// --- Helper Functions ---

inline int idx(int x, int y) {
    // Periodic boundaries
    if (x < 0) x += NX;
    if (x >= NX) x -= NX;
    if (y < 0) y += NY;
    if (y >= NY) y -= NY;
    return y * NX + x;
}

// Compute Macroscopic variables (Density, Velocity, Color Field)
void compute_macroscopic(int i) {
    double rho_r = 0.0, rho_b = 0.0;
    double mom_x = 0.0, mom_y = 0.0;

    for (int k = 0; k < Q; k++) {
        rho_r += grid[i].f_red[k];
        rho_b += grid[i].f_blue[k];
        double f_tot = grid[i].f_red[k] + grid[i].f_blue[k];
        mom_x += f_tot * cx[k];
        mom_y += f_tot * cy[k];
    }

    grid[i].rho_red = rho_r;
    grid[i].rho_blue = rho_b;
    grid[i].rho = rho_r + rho_b;

    // Check for division by zero (though rho should be ~1.0)
    if (grid[i].rho > 1e-9) {
        grid[i].ux = mom_x / grid[i].rho;
        grid[i].uy = mom_y / grid[i].rho;
        grid[i].phi = (rho_r - rho_b) / grid[i].rho;
    } else {
        grid[i].ux = 0.0;
        grid[i].uy = 0.0;
        grid[i].phi = 0.0;
    }
}

// Equilibrium Distribution (Standard Polynomial)
double get_feq(int k, double rho, double ux, double uy) {
    double cu = cx[k] * ux + cy[k] * uy;
    double u2 = ux * ux + uy * uy;
    return w[k] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
}

// Entropy Function H = sum [f * ln(f/w)]
double calculate_entropy(const double* f, const double* weights) {
    double H = 0.0;
    for (int k = 0; k < Q; k++) {
        if (f[k] > 1e-12) {
            H += f[k] * std::log(f[k] / weights[k]);
        }
    }
    return H;
}

// Entropic Stabilizer: Find alpha such that H(f + alpha(feq - f)) = H(f)
// Returns the stabilization parameter alpha.
double find_entropic_alpha(int i, double feq[Q], double f_total[Q]) {
    double delta[Q];
    for(int k=0; k<Q; k++) delta[k] = feq[k] - f_total[k];

    // Check if we are already at equilibrium (delta ~ 0)
    double max_delta = 0.0;
    for(int k=0; k<Q; k++) max_delta = std::max(max_delta, std::abs(delta[k]));
    if(max_delta < 1e-7) return 2.0;

    // Target Entropy (Current state)
    double H_in = calculate_entropy(f_total, w);

    // Initial Guess (Standard LBM = 2.0)
    double alpha = 2.0;

    // Newton-Raphson Solver
    // Function g(alpha) = H(f + alpha*delta) - H_in = 0
    for(int iter=0; iter<10; iter++) {
        double H_curr = 0.0;
        double dH_dalpha = 0.0;
        bool possible = true;

        for(int k=0; k<Q; k++) {
            double f_cand = f_total[k] + alpha * delta[k];
            if(f_cand <= 1e-12) { possible = false; break; } // Positivity check
            double term = std::log(f_cand / w[k]);
            H_curr += f_cand * term;
            dH_dalpha += delta[k] * (term + 1.0); // derivative of x*ln(x) is ln(x)+1
        }

        if (!possible) {
            alpha *= 0.5; // Backtrack if negativity occurs
            continue;
        }

        double g_val = H_curr - H_in;

        if(std::abs(g_val) < 1e-6) break; // Converged

        // Update alpha
        if(std::abs(dH_dalpha) < 1e-9) break;
        double change = g_val / dH_dalpha;

        // Clamp change to prevent wild oscillations
        change = std::max(-0.5, std::min(0.5, change));

        alpha = alpha - change;
    }

    // Safety: Entropic LBM requires alpha >= 1 (usually around 2)
    // If alpha < 1, it means the mirror state is "close", system is far from eq.
    // We strictly limit to valid physics range.
    return std::max(1.0, std::min(3.0, alpha));
}

// --- Main Simulation Loop ---

void initialize() {
    double center_x = NX / 2.0;
    double center_y = NY / 2.0;
    double radius = NY / 4.0;

    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int i = idx(x, y);
            double dist = std::sqrt(std::pow(x - center_x, 2) + std::pow(y - center_y, 2));

            // Initialization: Red Droplet, Blue Background
            // Tanh profile for smooth interface
            double phase = 0.5 * (1.0 - std::tanh((dist - radius) / 1.5));

            grid[i].rho_red = DENSITY_RED * phase;
            grid[i].rho_blue = DENSITY_BLUE * (1.0 - phase);
            grid[i].rho = grid[i].rho_red + grid[i].rho_blue;
            grid[i].ux = 0.0;
            grid[i].uy = 0.0;

            for (int k = 0; k < Q; k++) {
                double feq = get_feq(k, grid[i].rho, 0, 0);
                // Distribute eq according to density fraction
                grid[i].f_red[k] = feq * (grid[i].rho_red / grid[i].rho);
                grid[i].f_blue[k] = feq * (grid[i].rho_blue / grid[i].rho);
            }

            // Compute macroscopic quantities including phi
            compute_macroscopic(i);
        }
    }
}

void collision_step() {
    // 1. Calculate Gradients of Phase Field (Color Gradient)
    std::vector<double> GradX(NX*NY, 0.0);
    std::vector<double> GradY(NX*NY, 0.0);

    for(int y=0; y<NY; y++) {
        for(int x=0; x<NX; x++) {
            int i = idx(x,y);
            // Isotropic Gradient (using neighbors)
            for(int k=1; k<Q; k++) {
                int nb = idx(x + cx[k], y + cy[k]);
                double phi_nb = grid[nb].phi;
                GradX[i] += w[k] * cx[k] * phi_nb * 3.0;
                GradY[i] += w[k] * cy[k] * phi_nb * 3.0;
            }
        }
    }

    // 2. Local Collision Operations
    for(int i=0; i < NX*NY; i++) {
        compute_macroscopic(i);

        // Sum total f
        double f_total[Q];
        for(int k=0; k<Q; k++) f_total[k] = grid[i].f_red[k] + grid[i].f_blue[k];

        // --- Step A: Entropic Relaxation ---
        // Calculate Equilibrium
        double feq[Q];
        for(int k=0; k<Q; k++) feq[k] = get_feq(k, grid[i].rho, grid[i].ux, grid[i].uy);

        // Find Entropic Stabilizer Alpha
        double alpha = find_entropic_alpha(i, feq, f_total);

        // Standard ELBM two-step collision:
        // Step 1: Iso-entropic relaxation to f*
        double f_star[Q];
        for(int k=0; k<Q; k++) {
            f_star[k] = f_total[k] + alpha * (feq[k] - f_total[k]);
        }

        // Step 2: Dissipative relaxation with beta
        double beta = 1.0 / (2.0 * TAU_BASE);
        double f_post[Q];
        for(int k=0; k<Q; k++) {
            f_post[k] = (1.0 - beta) * f_total[k] + beta * f_star[k];
        }

        // --- Step B: Perturbation (Surface Tension) ---
        double grad_mag = std::sqrt(GradX[i]*GradX[i] + GradY[i]*GradY[i]);
        if (grad_mag > 1e-6) {
            double nx = GradX[i] / grad_mag;
            double ny = GradY[i] / grad_mag;

            // Parameter A for tension sigma: A = 9/2 * sigma / tau
            double A = (9.0 * SURFACE_TENSION) / (2.0 * TAU_BASE);

            for(int k=0; k<Q; k++) {
                double ci_n = cx[k]*nx + cy[k]*ny;
                double term = A * grad_mag * w[k] * (ci_n * ci_n - 1.0/3.0); // Simplified perturbation
                f_post[k] += term;
            }
        }

        // --- Step C: Recoloring (Segregation) ---
        double rho_tot = grid[i].rho;
        double ratio_r = grid[i].rho_red / rho_tot;
        double ratio_b = grid[i].rho_blue / rho_tot;

        double nx = (grad_mag > 1e-9)? GradX[i]/grad_mag : 0.0;
        double ny = (grad_mag > 1e-9)? GradY[i]/grad_mag : 0.0;

        for(int k=0; k<Q; k++) {
            double ci_dot_n = cx[k]*nx + cy[k]*ny;
            double cos_theta = ci_dot_n / std::sqrt(cx[k]*cx[k] + cy[k]*cy[k] + 1e-12);

            // Standard Latva-Kokko recoloring
            double f_eq = get_feq(k, rho_tot, grid[i].ux, grid[i].uy);
            double term_seg = BETA_SEG * ratio_r * ratio_b * cos_theta * std::abs(f_eq) / rho_tot;

            // Recoloring
            grid[i].f_new_red[k] = ratio_r * f_post[k] + term_seg;
            grid[i].f_new_blue[k] = ratio_b * f_post[k] - term_seg;
        }
    }
}

void streaming_step() {
    for (int y = 0; y < NY; y++) {
        for (int x = 0; x < NX; x++) {
            int i = idx(x, y);
            for (int k = 0; k < Q; k++) {
                // Pull from neighbor (Target based streaming)
                int prev_x = x - cx[k];
                int prev_y = y - cy[k];
                int i_prev = idx(prev_x, prev_y);

                grid[i].f_red[k] = grid[i_prev].f_new_red[k];
                grid[i].f_blue[k] = grid[i_prev].f_new_blue[k];
            }
        }
    }
}

// Write simple legacy VTK
void write_vtk(int step) {
    std::string filename = "droplet_" + std::to_string(step) + ".vtk";
    std::ofstream out(filename);

    out << "# vtk DataFile Version 3.0\n";
    out << "ELBM Droplet\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << NX << " " << NY << " 1\n";
    out << "ORIGIN 0 0 0\n";
    out << "SPACING 1 1 1\n";
    out << "POINT_DATA " << NX * NY << "\n";

    out << "SCALARS density double 1\n";
    out << "LOOKUP_TABLE default\n";
    for(const auto& n : grid) out << n.rho << "\n";

    out << "SCALARS phase double 1\n";
    out << "LOOKUP_TABLE default\n";
    for(const auto& n : grid) out << n.phi << "\n";

    out << "VECTORS velocity double\n";
    for(const auto& n : grid) out << n.ux << " " << n.uy << " 0.0\n";

    out.close();
}

double get_global_entropy() {
    double H = 0.0;
    for(int i=0; i<NX*NY; i++) {
        double f_tot[Q];
        for(int k=0; k<Q; k++) f_tot[k] = grid[i].f_red[k] + grid[i].f_blue[k];
        H += calculate_entropy(f_tot, w);
    }
    return H;
}

int main() {
    std::cout << "Initializing ELBM Static Droplet Simulation..." << std::endl;
    initialize();

    for (int t = 0; t <= MAX_ITER; t++) {
        if (t % OUTPUT_FREQ == 0) {
            double H = get_global_entropy();
            std::cout << "Step: " << t << " | Global H: " << std::fixed << std::setprecision(6) << H << std::endl;
            write_vtk(t);
        }

        collision_step();
        streaming_step();
    }

    std::cout << "Simulation Complete." << std::endl;
    return 0;
}