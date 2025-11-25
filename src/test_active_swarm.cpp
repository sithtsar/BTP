#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <random>

#include "core/lattice.h"
#include "core/fluid_state.h"
#include "solvers/elbm_solver.h"
#include "boundary/boundary_conditions.h"
#include "active_matter/active_particles.h"

using namespace elbm;

// --- Simulation Parameters ---
const int LX = 400;
const int LY = 100;
const int NUM_PARTICLES = 500;
const double SWIM_SPEED = 0.1;  // Reasonable swim speed
const double TUMBLE_RATE = 0.1;
const double COUPLING_STRENGTH = 0.01;  // Reduced for stability
const double VISCOSITY = 1.0/6.0;
const int MAX_STEPS = 2000;
const int OUTPUT_INTERVAL = 100;

// Utility to write fluid and particle data
template<typename Lattice>
void write_simulation_data(const LatticeGrid<Lattice::D, Lattice::Q>& grid,
                          const ActiveSwarm& swarm,
                          const std::string& filename) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error opening " << filename << std::endl;
        return;
    }
    file << std::scientific << std::setprecision(6);

    // Write fluid data
    file << "# Fluid data: x y u_x u_y rho\n";
    for (int y = 0; y < LY; ++y) {
        for (int x = 0; x < LX; ++x) {
            const auto& node = grid(x, y);
            file << x << " " << y << " "
                 << node.fluid.u[0] << " " << node.fluid.u[1] << " "
                 << node.fluid.rho << "\n";
        }
        file << "\n";
    }

    // Write particle data
    file << "# Particle data: x y vx vy orientation\n";
    const auto& particles = swarm.get_particles();
    for (const auto& particle : particles) {
        file << particle.position[0] << " " << particle.position[1] << " "
             << particle.velocity[0] << " " << particle.velocity[1] << " "
             << particle.orientation << "\n";
    }
}

int main() {
    std::cout << "--- Active Particle Swarm Simulation ---" << std::endl;
    std::cout << "Grid: " << LX << "x" << LY << ", Particles: " << NUM_PARTICLES << std::endl;
    std::cout << "Swim speed: " << SWIM_SPEED << ", Tumble rate: " << TUMBLE_RATE << std::endl;

    // Random number generator
    std::mt19937 gen(42);

    // 1. Create simulation components
    using Lattice = D2Q9;
    LatticeGrid<Lattice::D, Lattice::Q> grid(LX, LY);
    ELBMSolver<Lattice> solver(VISCOSITY);
    ActiveSwarm swarm(LX, LY, NUM_PARTICLES, SWIM_SPEED, TUMBLE_RATE);

    // 2. Initialize fluid at rest
    for (auto& node : grid) {
        node.fluid.rho = 1.0;
        node.fluid.u[0] = 0.0;
        node.fluid.u[1] = 0.0;
        BGKEquilibrium<Lattice>::compute(node.fluid, node.df.f);
    }

    // 3. Initialize particles in a cluster
    swarm.initialize_cluster(LX/2.0, LY/2.0, 20.0, gen);

    // --- Main Simulation Loop ---
    for (int step = 0; step <= MAX_STEPS; ++step) {
        if (step % OUTPUT_INTERVAL == 0) {
            std::cout << "Step: " << step << "/" << MAX_STEPS << std::endl;
            std::string filename = "./output/active_swarm_step_" + std::to_string(step) + ".dat";
            write_simulation_data<Lattice>(grid, swarm, filename);
        }

        // 1. Compute forces from particles on fluid
        std::vector<double> force_x, force_y;
        swarm.compute_fluid_forces(force_x, force_y, LX, LY, COUPLING_STRENGTH);

        // 2. LBM Collision (with particle forces)
        #pragma omp parallel for collapse(2)
        for (int y = 0; y < LY; ++y) {
            for (int x = 0; x < LX; ++x) {
                const size_t idx = y * LX + x;
                std::array<double, 2> local_force = {force_x[idx], force_y[idx]};
                solver.collide(grid(x,y), local_force);
            }
        }

        // 3. LBM Streaming (Periodic in X)
        solver.stream(grid, /*periodic_x=*/true, /*periodic_y=*/false);

        // 4. Apply Bounce-Back Boundary Conditions for Top/Bottom walls
        for (int x = 0; x < LX; ++x) {
            // Top wall (y = LY - 1)
            BounceBackBC<Lattice>::apply(grid(x, LY - 1));
            // Bottom wall (y = 0)
            BounceBackBC<Lattice>::apply(grid(x, 0));
        }

        // 5. Update particle positions and orientations
        const double dt = 0.01;  // Standard time step
        swarm.update_particles(dt, grid, gen, SWIM_SPEED);

        // 6. Apply boundary conditions to particles
        swarm.apply_periodic_boundaries(LX, LY);
    }

    std::cout << "Simulation finished." << std::endl;

    return 0;
}