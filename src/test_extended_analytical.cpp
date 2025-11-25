/**
 * Extended Analytical Test Cases - BGK vs ELBM vs Analytical
 * Tests: Hagen-Poiseuille (circular pipe) and Stokes shear flow
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>

#include "core/lattice.h"
#include "core/fluid_state.h"
#include "solvers/bgk_solver.h"
#include "solvers/elbm_solver.h"
#include "boundary/boundary_conditions.h"
#include "validation/analytical_cases.h"

using namespace elbm;
using namespace elbm::validation;

// Run Hagen-Poiseuille flow with a given solver
template<typename Solver>
void run_hagen_poiseuille(const std::string& solver_name, Solver& solver,
                          const HagenPoiseuilleFlow<D2Q9>& pipe,
                          const std::array<double, 2>& force,
                          int nx, int ny, int steps) {
    using Lattice = D2Q9;

    LatticeGrid<Lattice::D, Lattice::Q> grid(nx, ny);
    pipe.initialize(grid);

    std::cout << "  Running " << solver_name << " solver..." << std::endl;

    // Get circle center and radius for bounce-back
    const double cx = nx / 2.0;
    const double cy = ny / 2.0;
    const double radius = std::min(cx, cy) - 1.0;

    for (int step = 0; step < steps; ++step) {
        // Collision
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                if (!grid(x, y).is_solid) {
                    solver.collide(grid(x, y), force);
                }
            }
        }

        // Streaming
        solver.stream(grid);

        // Apply bounce-back at circular wall (AFTER streaming)
        // For solid wall nodes, we need to reflect incoming distributions
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                if (grid(x, y).is_solid) {
                    // Bounce-back: reverse directions
                    // For D2Q9: opposites are [0,2,1,4,3,6,5,8,7]
                    auto& f = grid(x, y).df.f;
                    std::array<double, 9> f_temp = f;

                    f[0] = f_temp[0];  // Rest particle stays
                    f[1] = f_temp[3];  // East <- West
                    f[2] = f_temp[4];  // North <- South
                    f[3] = f_temp[1];  // West <- East
                    f[4] = f_temp[2];  // South <- North
                    f[5] = f_temp[7];  // NE <- SW
                    f[6] = f_temp[8];  // NW <- SE
                    f[7] = f_temp[5];  // SW <- NE
                    f[8] = f_temp[6];  // SE <- NW
                }
            }
        }

        if (step % 10000 == 0 && step > 0) {
            double error = pipe.compute_error(grid);
            std::cout << "    Step " << step << ", L2 error: " << error << std::endl;
        }
    }

    double final_error = pipe.compute_error(grid);
    std::cout << "  Final L2 error: " << final_error << std::endl;

    // Save radial profile
    std::string filename = "output/hagen_poiseuille_" + solver_name + ".dat";
    std::ofstream file(filename);
    file << "# r u_sim u_exact error\n";

    // Sample along x-axis through center (reuse cx, cy, radius from above)
    for (int x = 0; x < nx; ++x) {
        const double r = std::abs(x - cx);
        if (r <= radius && !grid(x, ny/2).is_solid) {
            const double u_sim = grid(x, ny/2).fluid.u[0];
            const double u_exact = pipe.analytical_velocity(r);
            const double error = std::abs(u_sim - u_exact);
            file << r << " " << u_sim << " " << u_exact << " " << error << "\n";
        }
    }
    file.close();
    std::cout << "  Results saved to: " << filename << "\n";
}

// Run Stokes shear flow with a given solver
template<typename Solver>
void run_stokes_shear(const std::string& solver_name, Solver& solver,
                     const StokesShearFlow<D2Q9>& stokes,
                     int nx, int ny, int steps) {
// ... same as before ...
    for (int step = 0; step < steps; ++step) {
        // Collision
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                solver.collide(grid(x, y));
            }
        }
// ... same as before ...
}

void test_hagen_poiseuille_flow() {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "Testing Hagen-Poiseuille Flow - BGK vs ELBM" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    using Lattice = D2Q9;
    const int nx = 64;
    const int ny = 64;
    // Use smaller pressure gradient and larger viscosity for stability
    const double dp_dx = -0.0001;  // Reduced from -0.001
    const double viscosity = 0.15;  // Increased from 0.1
    const int steps = 30000;  // Reduced for faster testing

    std::cout << "Grid: " << nx << " x " << ny << std::endl;
    std::cout << "Pressure gradient: " << dp_dx << std::endl;
    std::cout << "Viscosity: " << viscosity << std::endl;
    std::cout << "Steps: " << steps << std::endl;
    std::cout << std::endl;

    HagenPoiseuilleFlow<Lattice> pipe(nx, ny, dp_dx, viscosity);
    const std::array<double, 2> force = {-dp_dx, 0.0};

    // Test with BGK
    {
        std::cout << "BGK Solver:\n";
        BGKSolver<Lattice> solver(viscosity);
        run_hagen_poiseuille("bgk", solver, pipe, force, nx, ny, steps);
    }

    // Test with ELBM
    {
        std::cout << "\nELBM Solver:\n";
        ELBMSolver<Lattice> solver(viscosity);
        run_hagen_poiseuille("elbm", solver, pipe, force, nx, ny, steps);
    }
}

void test_stokes_shear_flow() {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "Testing Stokes Shear Flow - BGK vs ELBM" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    using Lattice = D2Q9;
    const int nx = 100;
    const int ny = 40;
    const double shear_rate = 0.01;
    const double viscosity = 0.1;
    const int steps = 20000;

    std::cout << "Grid: " << nx << " x " << ny << std::endl;
    std::cout << "Shear rate: " << shear_rate << std::endl;
    std::cout << "Viscosity: " << viscosity << std::endl;
    std::cout << "Steps: " << steps << std::endl;
    std::cout << std::endl;

    StokesShearFlow<Lattice> stokes(nx, ny, shear_rate, viscosity);

    // Test with BGK
    {
        std::cout << "BGK Solver:\n";
        BGKSolver<Lattice> solver(viscosity);
        run_stokes_shear("bgk", solver, stokes, nx, ny, steps);
    }

    // Test with ELBM
    {
        std::cout << "\nELBM Solver:\n";
        ELBMSolver<Lattice> solver(viscosity);
        run_stokes_shear("elbm", solver, stokes, nx, ny, steps);
    }
}

int main() {
    std::cout << "\nExtended Analytical Validation Tests" << std::endl;
    std::cout << "BGK vs ELBM vs Analytical Solutions" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    try {
        test_hagen_poiseuille_flow();
        test_stokes_shear_flow();

        std::cout << "\n" << std::string(60, '=') << std::endl;
        std::cout << "All Extended Tests Completed Successfully" << std::endl;
        std::cout << std::string(60, '=') << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
