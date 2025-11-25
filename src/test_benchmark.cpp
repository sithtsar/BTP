#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <string>
#include <chrono>

#include "core/lattice.h"
#include "core/fluid_state.h"
#include "core/entropy_monitor.h"
#include "solvers/bgk_solver.h"
#include "solvers/elbm_solver.h"
#include "solvers/equilibrium.h"
#include "boundary/boundary_conditions.h"
#include "validation/benchmark_cases.h"

using namespace elbm;
using namespace elbm::validation;

// ... (keep existing test_lid_driven_cavity and test_flow_past_cylinder functions) ...
void test_lid_driven_cavity(int Re) {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "Testing Lid-Driven Cavity Flow - Re = " << Re << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    using Lattice = D2Q9;
    const int n = 128;  // Grid size
    const double U_lid = 0.1;
    const int steps = 50000;
    const int output_interval = 5000;

    // Create cavity case
    LidDrivenCavity<Lattice> cavity(n, U_lid, Re);

    // Create grid
    LatticeGrid<Lattice::D, Lattice::Q> grid(n, n);

    // Initialize
    cavity.initialize(grid);

    // Create solver
    BGKSolver<Lattice> solver(cavity.viscosity());

    // Setup boundaries
    BoundaryManager<Lattice> bc_manager;
    cavity.setupBoundaries<BGKSolver<Lattice>>(bc_manager);

    std::cout << "Grid: " << n << " x " << n << std::endl;
    std::cout << "Lid velocity: " << U_lid << std::endl;
    std::cout << "Reynolds number: " << Re << std::endl;
    std::cout << "Viscosity: " << cavity.viscosity() << std::endl;
    std::cout << "Steps: " << steps << std::endl;
    std::cout << std::endl;

    // Output file for centerline velocities
    std::string output_file = "output/cavity_re" + std::to_string(Re) + "_profiles.dat";
    std::ofstream file(output_file);
    file << "# step y u_x_centerline v_y_centerline\n";

    // Run simulation
    for (int step = 0; step <= steps; ++step) {
        // Collision
        for (int y = 0; y < n; ++y) {
            for (int x = 0; x < n; ++x) {
                solver.collide(grid(x, y));
            }
        }

        // Boundary conditions
        bc_manager.applyBoundaries(grid);

        // Streaming
        solver.stream(grid);

        if (step % output_interval == 0) {
            std::cout << "Step " << step << "/" << steps << std::endl;

            // Get centerline velocities
            auto u_center = cavity.get_u_centerline(grid);
            auto v_center = cavity.get_v_centerline(grid);

            // Save to file
            for (int i = 0; i < n; ++i) {
                file << step << " " << i << " " << u_center[i] << " " << v_center[i] << "\n";
            }
            file << "\n";  // Blank line between timesteps
        }
    }

    file.close();
    std::cout << "\nResults saved to: " << output_file << std::endl;

    // Save final field
    std::string field_file = "output/cavity_re" + std::to_string(Re) + "_final.dat";
    std::ofstream field(field_file);
    field << "# x y rho ux uy p\n";

    for (int y = 0; y < n; ++y) {
        for (int x = 0; x < n; ++x) {
            const auto& fluid = grid(x, y).fluid;
            field << x << " " << y << " "
                 << fluid.rho << " "
                 << fluid.u[0] << " "
                 << fluid.u[1] << " "
                 << fluid.p << "\n";
        }
    }
    field.close();
    std::cout << "Final field saved to: " << field_file << std::endl;
}

void test_flow_past_cylinder(int Re) {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "Testing Flow Past Cylinder - Re = " << Re << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    using Lattice = D2Q9;
    const int nx = 400;
    const int ny = 120;
    const double U_inf = 0.1;
    const double diameter = 20.0;
    const int steps = 10000;
    const int output_interval = 100;

    // Create cylinder case
    FlowPastCylinder<Lattice> cylinder(nx, ny, U_inf, diameter, Re);

    // Create grid
    LatticeGrid<Lattice::D, Lattice::Q> grid(nx, ny);

    // Initialize
    cylinder.initialize(grid);

    // Create solver
    BGKSolver<Lattice> solver(cylinder.viscosity());

    std::cout << "Grid: " << nx << " x " << ny << std::endl;
    std::cout << "Freestream velocity: " << U_inf << std::endl;
    std::cout << "Cylinder diameter: " << diameter << std::endl;
    std::cout << "Reynolds number: " << Re << std::endl;
    std::cout << "Viscosity: " << cylinder.viscosity() << std::endl;
    std::cout << "Steps: " << steps << std::endl;
    std::cout << std::endl;

    // Output file for forces
    std::string force_file = "output/cylinder_re" + std::to_string(Re) + "_forces.dat";
    std::ofstream file(force_file);
    file << "# step Cd Cl\n";

    // Setup boundary conditions (inlet/outlet/walls)
    BoundaryManager<Lattice> bc_manager;
    std::array<double, 2> inlet_velocity = {U_inf, 0.0};
    bc_manager.setLeftBC(BCType::VELOCITY, inlet_velocity);
    bc_manager.setRightBC(BCType::PRESSURE, 1.0);
    bc_manager.setTopBC(BCType::BOUNCE_BACK);
    bc_manager.setBottomBC(BCType::BOUNCE_BACK);

    // Run simulation
    for (int step = 0; step <= steps; ++step) {
        // Step 1: Collision (skip solid nodes)
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                if (!grid(x, y).is_solid) {
                    solver.collide(grid(x, y));
                }
            }
        }

        // Step 2: Streaming (all nodes)
        solver.stream(grid);

        // Step 3: Apply cylinder immersed boundary (AFTER streaming)
        cylinder.applyImmersedBoundary(grid);

        // Step 4: Apply inlet/outlet and wall boundaries
        bc_manager.applyBoundaries(grid);

        if (step % output_interval == 0) {
            auto coeffs = cylinder.compute_coefficients(grid);
            file << step << " " << coeffs[0] << " " << coeffs[1] << "\n";

            if (step % 1000 == 0) {
                std::cout << "Step " << step << ", Cd = " << coeffs[0]
                         << ", Cl = " << coeffs[1] << std::endl;
            }
        }
    }

    file.close();
    std::cout << "\nForce coefficients saved to: " << force_file << std::endl;

    // Save final velocity field
    std::string field_file = "output/cylinder_re" + std::to_string(Re) + "_final.dat";
    std::ofstream field(field_file);
    field << "# x y rho ux uy p is_solid\n";

    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            const auto& node = grid(x, y);
            field << x << " " << y << " "
                 << node.fluid.rho << " "
                 << node.fluid.u[0] << " "
                 << node.fluid.u[1] << " "
                 << node.fluid.p << " "
                 << node.is_solid << "\n";
        }
    }
    field.close();
    std::cout << "Final field saved to: " << field_file << std::endl;
}

template<typename SolverType>
void run_channel_flow_simulation(int Re, const std::string& solver_name) {
    using Lattice = D2Q9;
    const int nx = 200;  // Match provided code dimensions
    const int ny = 50;
    const double diameter = 10.0;  // radius = ny/10 = 5, diameter = 10
    const int steps = 20000;  // Match provided code
    const int output_interval = 1000;

    // Fixed inlet velocity, vary viscosity to change Re
    // Re = U * D / ν  =>  ν = U * D / Re
    const double U_inf = 0.1;  // Fixed inlet velocity
    const double viscosity = U_inf * diameter / Re;

    // Check Mach number constraint (Ma = U/cs < 0.3 for incompressibility)
    const double cs = 1.0 / std::sqrt(3.0);  // Speed of sound in lattice units
    const double Ma = U_inf / cs;

    if (Ma > 0.3) {
        std::cerr << "WARNING: Mach number " << Ma << " exceeds 0.3 for Re=" << Re
                  << ". Results may be affected by compressibility." << std::endl;
    }

    ChannelFlowWithCylinder<Lattice> channel_flow(nx, ny, U_inf, diameter, Re);
    LatticeGrid<Lattice::D, Lattice::Q> grid(nx, ny);
    channel_flow.initialize(grid);

    SolverType solver(viscosity);  // Use viscosity directly

    std::cout << "Grid: " << nx << " x " << ny << std::endl;
    std::cout << "Freestream velocity: " << U_inf << " (varied)" << std::endl;
    std::cout << "Cylinder diameter: " << diameter << std::endl;
    std::cout << "Reynolds number: " << Re << std::endl;
    std::cout << "Viscosity: " << viscosity << " (fixed)" << std::endl;
    std::cout << "Mach number: " << Ma << std::endl;
    std::cout << "Steps: " << steps << std::endl;

    // Setup boundary conditions
    BoundaryManager<Lattice> bc_manager;
    channel_flow.setupBoundaries(bc_manager);

    for (int step = 0; step <= steps; ++step) {
        // Step 1: Collision (skip solid nodes)
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                if (!grid(x, y).is_solid) {
                    solver.collide(grid(x, y));
                }
            }
        }

        // Step 2: Streaming
        solver.stream(grid);

        // Step 3: Apply cylinder immersed boundary (AFTER streaming)
        channel_flow.applyImmersedBoundary(grid);

        // Step 4: Apply inlet/outlet and wall boundaries (AFTER streaming)
        bc_manager.applyBoundaries(grid);

        // Step 5: Update macroscopic quantities after boundary conditions
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                if (!grid(x, y).is_solid) {
                    computeMacroscopic<Lattice>(grid(x, y).df.f, grid(x, y).fluid);
                    // Impose outlet velocity (match provided code)
                    if (x == nx - 1) {
                        grid(x, y).fluid.u[0] = U_inf;
                        grid(x, y).fluid.u[1] = 0.0;
                    }
                }
            }
        }

        if (step % output_interval == 0) {
            std::cout << "Step " << step << "/" << steps << std::endl;
        }
    }



    // Final macroscopic update for all nodes (including boundaries) before output
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            computeMacroscopic<Lattice>(grid(x, y).df.f, grid(x, y).fluid);
        }
    }

    // Compute vorticity
    auto vorticity = channel_flow.compute_vorticity(grid);

    // Output CSV format matching provided code
    std::string csv_file = "LBM_output_Re" + std::to_string(Re) + "_" + solver_name + ".csv";
    std::ofstream out(csv_file);
    out << "x,y,u,v,speed,vorticity\n";

    // Skip boundary cells for output (match provided code)
    for (int x = 1; x < nx - 1; ++x) {
        for (int y = 1; y < ny - 1; ++y) {
            if (grid(x, y).is_solid) continue;  // Skip solid cells

            double u = grid(x, y).fluid.u[0];
            double v = grid(x, y).fluid.u[1];
            double speed = std::sqrt(u * u + v * v);
            double vort = vorticity[y][x];

            out << x << "," << y << ","
                << u << "," << v << ","
                << speed << "," << vort << "\n";
        }
    }
    out.close();
    std::cout << "Re=" << Re << " done, output written to file: " << csv_file << std::endl;
}


void test_channel_flow_with_cylinder(int Re, const std::string& solver_type) {
    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << "Testing Channel Flow with Cylinder - Re = " << Re << " (" << solver_type << ")" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    using Lattice = D2Q9;
    if (solver_type == "BGK") {
        run_channel_flow_simulation<BGKSolver<Lattice>>(Re, "BGK");
    } else {
        run_channel_flow_simulation<ELBMSolver<Lattice>>(Re, "ELBM");
    }
}

int main(int argc, char* argv[]) {
    system("mkdir -p output");

    std::cout << "==================================================" << std::endl;
    std::cout << "ELBM Benchmark Test Suite" << std::endl;
    std::cout << "==================================================" << std::endl;

    // Check if specific test configuration is requested
    if (argc > 1) {
        std::string test_type = argv[1];

        if (test_type == "cylinder" && argc >= 3) {
            // Usage: ./test_benchmark cylinder <Re> [0=BGK|1=ELBM]
            int Re = std::stoi(argv[2]);
            int solver_choice = (argc > 3) ? std::stoi(argv[3]) : -1;

            auto start = std::chrono::high_resolution_clock::now();

            if (solver_choice == 0 || solver_choice == -1) {
                test_channel_flow_with_cylinder(Re, "BGK");
            }
            if (solver_choice == 1 || solver_choice == -1) {
                test_channel_flow_with_cylinder(Re, "ELBM");
            }

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
            std::cerr << "Total time: " << duration.count() << "s" << std::endl;

            return 0;
        }
    }

    // Default: Run all tests
    // Channel flow with cylinder - BGK and ELBM comparison at multiple Re
    // Match provided code Re values
    std::vector<int> reynolds_numbers = {10, 100, 1000};

    for (int Re : reynolds_numbers) {
        test_channel_flow_with_cylinder(Re, "BGK");
        test_channel_flow_with_cylinder(Re, "ELBM");
    }

    std::cout << "\n==================================================" << std::endl;
    std::cout << "All benchmark tests complete!" << std::endl;
    std::cout << "==================================================" << std::endl;

    return 0;
}

