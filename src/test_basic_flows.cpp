/*
 * UNIFIED LATTICE BOLTZMANN SOLVER: BGK & ENTROPIC
 * Verifying Couette, Poiseuille, and Taylor-Green Vortex Flows
 *
 * Features:
 * - D2Q9 Lattice
 * - BGK and Entropic (ELBM) Collision Operators
 * - Newton-Raphson Alpha Search for ELBM
 * - Guo Forcing Scheme
 * - Zou-He and Bounce-Back Boundary Conditions
 * - L2 Error Norm Verification
 *
 * Author: CFD Research Specialist
 * Date: October 27, 2023
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sstream>

// --- PHYSICAL CONSTANTS ---
const int Q = 9;
const int nx_default = 64;
const int ny_default = 64;

// D2Q9 Weights and Vectors
const double w[Q] = {
    4.0/9.0,
    1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,
    1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0
};

const int cx[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
const int cy[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
const int opp[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6}; // Opposite direction index

// --- CONFIGURATION STRUCTS ---
enum FlowType { COUETTE, POISEUILLE, TAYLOR_GREEN };
enum CollisionModel { BGK, ENTROPIC };

struct SimulationConfig {
    int nx;
    int ny;
    double tau;           // Relaxation time
    double u0;            // Characteristic velocity (Lid or TGV max)
    double force_x;       // Body force for Poiseuille
    int max_steps;
    int check_interval;
    FlowType flow;
    CollisionModel collision;
    std::string output_prefix;
};

// --- HELPER FUNCTIONS ---

// 1. Equilibrium Distribution (Low Mach Expansion)
double get_feq(int i, double rho, double ux, double uy) {
    double cu = cx[i]*ux + cy[i]*uy;
    double u2 = ux*ux + uy*uy;
    double cs2 = 1.0/3.0;

    // 2nd order expansion
    return w[i] * rho * (1.0 + cu/cs2 + (cu*cu)/(2.0*cs2*cs2) - u2/(2.0*cs2));
}

// 2. Entropy Calculation: H = Sum f * ln(f/w)
double calculate_H(const std::vector<double>& f_local) {
    double H = 0.0;
    for(int i=0; i<Q; ++i) {
        if(f_local[i] > 1e-16) {
            H += f_local[i] * std::log(f_local[i] / w[i]);
        }
    }
    return H;
}

// 3. Newton-Raphson Solver for Alpha (ELBM)
// Target: H(f + alpha(feq - f)) - H(f) = 0
double find_alpha_entropy(const std::vector<double>& f, const std::vector<double>& feq, double H_target) {
    double alpha = 2.0; // Initial guess (BGK limit)
    const int max_iter = 20;
    const double tol = 1e-6;

    for(int iter=0; iter<max_iter; ++iter) {
        double H_new = 0.0;
        double dH_dalpha = 0.0;
        bool negative_pop = false;

        for(int i=0; i<Q; ++i) {
            double delta = feq[i] - f[i];
            double f_next = f[i] + alpha * delta;

            if(f_next <= 0.0) {
                negative_pop = true;
                break;
            }

            double log_term = std::log(f_next / w[i]);
            H_new += f_next * log_term;
            dH_dalpha += delta * (log_term + 1.0);
        }

        if (negative_pop) {
            alpha *= 0.5; // Backtrack if we hit negative populations
            continue;
        }

        double residual = H_new - H_target;
        if(std::abs(residual) < tol) return alpha;

        // Newton Step
        if(std::abs(dH_dalpha) < 1e-10) break; // Avoid div by zero
        alpha = alpha - residual / dH_dalpha;

        // Safety Clamping
        if(alpha < 1.0) alpha = 1.0;
        if(alpha > 10.0) alpha = 10.0;
    }
    return alpha;
}

// --- SOLVER CLASS ---
class LBMSolver {
private:
    SimulationConfig cfg;
    std::vector<double> f;       // Current populations
    std::vector<double> f_new;   // Next populations (buffer)

    // Macroscopic Fields for Analysis
    std::vector<double> rho_field;
    std::vector<double> ux_field;
    std::vector<double> uy_field;

    // Indexing: Row-Major (y, x, i)
    int idx(int x, int y, int i) const { return (y * cfg.nx + x) * Q + i; }
    int scalar_idx(int x, int y) const { return y * cfg.nx + x; }

public:
    LBMSolver(SimulationConfig c) : cfg(c) {
        int size = cfg.nx * cfg.ny * Q;
        f.resize(size);
        f_new.resize(size);
        rho_field.resize(cfg.nx * cfg.ny);
        ux_field.resize(cfg.nx * cfg.ny);
        uy_field.resize(cfg.nx * cfg.ny);

        initialize();
    }

    void initialize() {
        for(int y=0; y<cfg.ny; ++y) {
            for(int x=0; x<cfg.nx; ++x) {
                double rho = 1.0;
                double ux = 0.0;
                double uy = 0.0;

                // TGV Initialization (Exact Solution at t=0)
                if(cfg.flow == TAYLOR_GREEN) {
                    double kx = 2.0 * M_PI * x / (double)cfg.nx;
                    double ky = 2.0 * M_PI * y / (double)cfg.ny;
                    ux = -cfg.u0 * cos(kx) * sin(ky);
                    uy =  cfg.u0 * sin(kx) * cos(ky);
                    // Pressure P = -rho u0^2 / 4 (cos2x + cos2y)
                    double P = -0.25 * 1.0 * cfg.u0 * cfg.u0 * (cos(2.0*kx) + cos(2.0*ky));
                    rho = 1.0 + 3.0 * P;
                }

                // Couette Initialization
                if(cfg.flow == COUETTE) {
                    double H = cfg.ny - 1;
                    if(y == H) ux = cfg.u0;
                }

                for(int i=0; i<Q; ++i) {
                    f[idx(x,y,i)] = get_feq(i, rho, ux, uy);
                }

                // Initialize fields
                rho_field[scalar_idx(x,y)] = rho;
                ux_field[scalar_idx(x,y)] = ux;
                uy_field[scalar_idx(x,y)] = uy;
            }
        }
    }

    void collide_and_stream() {
        // Main LBM Loop

        for(int y=0; y<cfg.ny; ++y) {
            for(int x=0; x<cfg.nx; ++x) {

                // 1. Compute Moments (rho, u)
                double rho = 0.0;
                double ux = 0.0;
                double uy = 0.0;
                std::vector<double> f_curr(Q);

                for(int i=0; i<Q; ++i) {
                    f_curr[i] = f[idx(x,y,i)];
                    rho += f_curr[i];
                    ux += f_curr[i] * cx[i];
                    uy += f_curr[i] * cy[i];
                }

                // 2. Force Correction (Guo Scheme Step 1: Velocity shift)
                double Fx = 0.0;
                if(cfg.flow == POISEUILLE) Fx = cfg.force_x;

                double ux_phys = ux / rho + (0.5 * Fx) / rho;
                double uy_phys = uy / rho;

                // Save for analysis
                rho_field[scalar_idx(x,y)] = rho;
                ux_field[scalar_idx(x,y)] = ux_phys;
                uy_field[scalar_idx(x,y)] = uy_phys;

                // 3. Equilibrium and Collision
                std::vector<double> feq(Q);
                for(int i=0; i<Q; ++i) feq[i] = get_feq(i, rho, ux_phys, uy_phys);

                double beta = 1.0 / (2.0 * cfg.tau);
                double alpha = 2.0; // Standard BGK value

                if(cfg.collision == ENTROPIC) {
                    double H_curr = calculate_H(f_curr);
                    alpha = find_alpha_entropy(f_curr, feq, H_curr);
                }

                // 4. Guo Forcing Source Term S_i
                std::vector<double> S(Q, 0.0);
                if(cfg.flow == POISEUILLE) {
                    double cs2 = 1.0/3.0;
                    for(int i=0; i<Q; ++i) {
                        double ci_minus_u_x = cx[i] - ux_phys;
                        double ci_dot_u = cx[i]*ux_phys + cy[i]*uy_phys;

                        double term1 = (ci_minus_u_x * Fx) / cs2;
                        double term2 = (ci_dot_u * (cx[i] * Fx)) / (cs2*cs2);

                        S[i] = w[i] * (term1 + term2);
                    }
                }

                // 5. Post-Collision and Streaming
                double force_factor = (1.0 - 0.5/cfg.tau);

                for(int i=0; i<Q; ++i) {
                    double omega = alpha * beta;

                    if(cfg.collision == BGK) omega = 1.0 / cfg.tau;

                    double f_val = f_curr[i] + omega * (feq[i] - f_curr[i]) + force_factor * S[i];

                    // STREAMING STEP
                    int next_x = (x + cx[i] + cfg.nx) % cfg.nx;
                    int next_y = (y + cy[i] + cfg.ny) % cfg.ny;

                    f_new[idx(next_x, next_y, i)] = f_val;
                }
            }
        }
    }

    void apply_boundary_conditions() {
        if(cfg.flow == TAYLOR_GREEN) return; // Fully periodic

        // --- TOP WALL (Y = NY-1) ---
        double u_wall_x = (cfg.flow == COUETTE)? cfg.u0 : 0.0;
        int y_top = cfg.ny - 1;

        for(int x=0; x<cfg.nx; ++x) {
            double f0 = f_new[idx(x, y_top, 0)];
            double f1 = f_new[idx(x, y_top, 1)];
            double f2 = f_new[idx(x, y_top, 2)];
            double f3 = f_new[idx(x, y_top, 3)];
            double f5 = f_new[idx(x, y_top, 5)];
            double f6 = f_new[idx(x, y_top, 6)];

            double rho_wall = f0 + f1 + f3 + 2.0*(f2 + f5 + f6);

            // Zou-He Updates
            f_new[idx(x, y_top, 4)] = f2;
            f_new[idx(x, y_top, 7)] = f5 + 0.5*(f1 - f3) - 0.5*rho_wall*u_wall_x;
            f_new[idx(x, y_top, 8)] = f6 - 0.5*(f1 - f3) + 0.5*rho_wall*u_wall_x;
        }

        // --- BOTTOM WALL (Y = 0) ---
        int y_bot = 0;

        for(int x=0; x<cfg.nx; ++x) {
            double f0 = f_new[idx(x, y_bot, 0)];
            double f1 = f_new[idx(x, y_bot, 1)];
            double f3 = f_new[idx(x, y_bot, 3)];
            double f4 = f_new[idx(x, y_bot, 4)];
            double f7 = f_new[idx(x, y_bot, 7)];
            double f8 = f_new[idx(x, y_bot, 8)];

            double rho_bot = f0 + f1 + f3 + 2.0*(f4 + f7 + f8);

            f_new[idx(x, y_bot, 2)] = f4;
            f_new[idx(x, y_bot, 5)] = f7 - 0.5*(f1 - f3);
            f_new[idx(x, y_bot, 6)] = f8 + 0.5*(f1 - f3);
        }
    }

    void step() {
        collide_and_stream();
        apply_boundary_conditions();
        f = f_new;
    }

    // --- ANALYTICAL VERIFICATION ---
    double compute_error(int timestep) {
        double error_L2 = 0.0;
        double norm_L2 = 0.0;
        double nu = (cfg.tau - 0.5)/3.0;

        for(int y=0; y<cfg.ny; ++y) {
            for(int x=0; x<cfg.nx; ++x) {
                double u_exact = 0.0;
                double u_num = ux_field[scalar_idx(x,y)];

                // Analytical Solutions
                if(cfg.flow == COUETTE) {
                    u_exact = cfg.u0 * ((double)y / (cfg.ny-1));
                }
                else if(cfg.flow == POISEUILLE) {
                    double H = cfg.ny - 1.0;
                    double yy = (double)y;
                    u_exact = (cfg.force_x / (2.0 * 1.0 * nu)) * yy * (H - yy);
                }
                else if(cfg.flow == TAYLOR_GREEN) {
                    double k = 2.0*M_PI/cfg.nx;
                    double decay = std::exp(-2.0 * nu * k * k * timestep);
                    u_exact = -cfg.u0 * cos(k*x) * sin(k*y) * decay;
                }

                error_L2 += (u_num - u_exact)*(u_num - u_exact);
                norm_L2 += u_exact * u_exact;
            }
        }

        if(norm_L2 < 1e-9) norm_L2 = 1.0;
        return std::sqrt(error_L2 / norm_L2);
    }

    void save_results(const std::string& filename) {
        std::ofstream file(filename);
        file << std::scientific << std::setprecision(8);
        file << "# x y rho ux uy\n";

        for(int y=0; y<cfg.ny; ++y) {
            for(int x=0; x<cfg.nx; ++x) {
                file << x << " " << y << " "
                     << rho_field[scalar_idx(x,y)] << " "
                     << ux_field[scalar_idx(x,y)] << " "
                     << uy_field[scalar_idx(x,y)] << "\n";
            }
        }
        file.close();
    }

    void save_profile(const std::string& filename, int timestep) {
        std::ofstream file(filename);
        file << std::scientific << std::setprecision(8);
        file << "# y ux u_exact error\n";

        double nu = (cfg.tau - 0.5)/3.0;
        int x_mid = cfg.nx / 2;

        for(int y=0; y<cfg.ny; ++y) {
            double u_num = ux_field[scalar_idx(x_mid, y)];
            double u_exact = 0.0;

            if(cfg.flow == COUETTE) {
                u_exact = cfg.u0 * ((double)y / (cfg.ny-1));
            }
            else if(cfg.flow == POISEUILLE) {
                double H = cfg.ny - 1.0;
                double yy = (double)y;
                u_exact = (cfg.force_x / (2.0 * 1.0 * nu)) * yy * (H - yy);
            }

            double error = std::abs(u_num - u_exact);
            file << y << " " << u_num << " " << u_exact << " " << error << "\n";
        }
        file.close();
    }

    // Compute global entropy for H-theorem validation
    struct EntropyStats {
        double H_total;
        double H_mean;
        double max_velocity;
    };

    EntropyStats compute_global_entropy() const {
        EntropyStats stats;
        stats.H_total = 0.0;
        stats.max_velocity = 0.0;
        int node_count = 0;

        for(int y=0; y<cfg.ny; ++y) {
            for(int x=0; x<cfg.nx; ++x) {
                // Get distributions at this node
                std::vector<double> f_local(Q);
                for(int q=0; q<Q; ++q) {
                    f_local[q] = f[idx(x, y, q)];
                }

                // Compute H for this node
                double H_node = calculate_H(f_local);
                stats.H_total += H_node;
                node_count++;

                // Track max velocity
                double vel_mag = std::sqrt(
                    ux_field[scalar_idx(x, y)] * ux_field[scalar_idx(x, y)] +
                    uy_field[scalar_idx(x, y)] * uy_field[scalar_idx(x, y)]
                );
                stats.max_velocity = std::max(stats.max_velocity, vel_mag);
            }
        }

        stats.H_mean = stats.H_total / node_count;
        return stats;
    }
};

void run_simulation(SimulationConfig config) {
    std::string flow_name;
    switch(config.flow) {
        case COUETTE: flow_name = "Couette"; break;
        case POISEUILLE: flow_name = "Poiseuille"; break;
        case TAYLOR_GREEN: flow_name = "Taylor-Green"; break;
    }

    std::string collision_name = (config.collision == BGK) ? "BGK" : "ELBM";

    std::cout << "\n=== " << flow_name << " Flow | " << collision_name << " ===" << std::endl;
    std::cout << "Grid: " << config.nx << "x" << config.ny << std::endl;
    std::cout << "tau: " << config.tau << ", nu: " << (config.tau - 0.5)/3.0 << std::endl;

    LBMSolver solver(config);

    // Open entropy tracking file
    std::stringstream entropy_ss;
    entropy_ss << config.output_prefix << "_entropy_" << collision_name << ".dat";
    std::ofstream entropy_file(entropy_ss.str());
    entropy_file << "# timestep H_total H_mean max_velocity\n";
    entropy_file << std::scientific << std::setprecision(12);

    // Timesteps for saving snapshots (for temporal evolution plots)
    std::vector<int> snapshot_times = {0, 5000, 10000, 15000, 20000};

    for(int i=0; i<=config.max_steps; ++i) {
        solver.step();

        // Save snapshots at specific timesteps for temporal evolution plots
        if(std::find(snapshot_times.begin(), snapshot_times.end(), i) != snapshot_times.end()) {
            std::stringstream snapshot_ss;
            snapshot_ss << config.output_prefix << "_" << collision_name << "_t"
                       << std::setfill('0') << std::setw(5) << i << ".dat";
            solver.save_results(snapshot_ss.str());
        }

        if(i % config.check_interval == 0) {
            double error = solver.compute_error(i);

            // Compute global H-function using public method
            auto entropy_stats = solver.compute_global_entropy();

            entropy_file << i << " "
                        << entropy_stats.H_total << " "
                        << entropy_stats.H_mean << " "
                        << entropy_stats.max_velocity << "\n";
            entropy_file.flush();

            std::cout << "Step " << std::setw(6) << i
                     << " | L2 error: " << std::scientific << std::setprecision(6)
                     << error
                     << " | H: " << entropy_stats.H_total << std::endl;
        }
    }

    entropy_file.close();

    // Save final results
    std::stringstream ss;
    ss << config.output_prefix << "_" << collision_name << ".dat";
    solver.save_results(ss.str());

    std::stringstream ss2;
    ss2 << config.output_prefix << "_profile_" << collision_name << ".dat";
    solver.save_profile(ss2.str(), config.max_steps);

    std::cout << "Results saved to: " << ss.str() << std::endl;
    std::cout << "=== Simulation Complete ===" << std::endl;
}

int main(int argc, char* argv[]) {
    std::string test_case = "all";
    if(argc > 1) {
        test_case = argv[1];
    }

    std::cout << "=== LBM Analytical Validation: BGK & ELBM ===" << std::endl;
    std::cout << "Test case: " << test_case << std::endl;

    if(test_case == "couette" || test_case == "all") {
        // Couette Flow - Both BGK and ELBM
        SimulationConfig config;
        config.nx = 64;
        config.ny = 64;
        config.tau = 0.8;
        config.u0 = 0.1;
        config.force_x = 0.0;
        config.max_steps = 20000;
        config.check_interval = 2000;
        config.flow = COUETTE;
        config.output_prefix = "output/analytical_validation/couette";

        config.collision = BGK;
        run_simulation(config);

        config.collision = ENTROPIC;
        run_simulation(config);
    }

    if(test_case == "poiseuille" || test_case == "all") {
        // Poiseuille Flow - Both BGK and ELBM
        SimulationConfig config;
        config.nx = 100;
        config.ny = 50;
        config.tau = 0.8;
        config.u0 = 0.0;
        config.force_x = 1e-5;
        config.max_steps = 20000;
        config.check_interval = 2000;
        config.flow = POISEUILLE;
        config.output_prefix = "output/analytical_validation/poiseuille";

        config.collision = BGK;
        run_simulation(config);

        config.collision = ENTROPIC;
        run_simulation(config);
    }

    if(test_case == "taylor" || test_case == "all") {
        // Taylor-Green Vortex - Both BGK and ELBM
        SimulationConfig config;
        config.nx = 64;
        config.ny = 64;
        config.tau = 0.8;
        config.u0 = 0.05;
        config.force_x = 0.0;
        config.max_steps = 10000;
        config.check_interval = 1000;
        config.flow = TAYLOR_GREEN;
        config.output_prefix = "output/analytical_validation/taylor";

        config.collision = BGK;
        run_simulation(config);

        config.collision = ENTROPIC;
        run_simulation(config);
    }

    std::cout << "\n=== All Tests Complete ===" << std::endl;
    return 0;
}
