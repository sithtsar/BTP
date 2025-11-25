#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#include "core/lattice.h"
#include "core/fluid_state.h"
#include "solvers/elbm_solver.h"
#include "boundary/boundary_conditions.h"
#include "active_matter/nematic_tensor.h"

using namespace elbm;

// --- Benchmark Parameters (from user prompt section 7.1) ---
const int LX = 400;
const int LY = 100;
const double ZETA = 0.05;      // Activity (above threshold)
const double K_ELASTIC = 0.04;
const double GAMMA_ROT = 0.3;   // Rotational Diffusion (not used yet)
const double VISCOSITY = 1.0/6.0; // Corresponds to tau = 1.0
const int MAX_STEPS = 1000;
const int OUTPUT_INTERVAL = 100;

// Utility to write fluid and nematic fields
template<typename Lattice>
void write_fields(const LatticeGrid<Lattice::D, Lattice::Q>& grid, 
                  const ActiveNematic<Lattice>& nematic, 
                  const std::string& filename) {
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error opening " << filename << std::endl;
        return;
    }
    file << std::scientific << std::setprecision(6);
    file << "# x y u_x u_y rho Q_xx Q_xy\n";

    const auto& Q_xx = nematic.get_Q_xx();
    const auto& Q_xy = nematic.get_Q_xy();

    for (int y = 0; y < LY; ++y) {
        for (int x = 0; x < LX; ++x) {
            const auto& node = grid(x, y);
            const size_t idx = y * LX + x;
            file << x << " " << y << " "
                 << node.fluid.u[0] << " " << node.fluid.u[1] << " "
                 << node.fluid.rho << " "
                 << Q_xx[idx] << " " << Q_xy[idx] << "\n";
        }
        file << "\n";
    }
}

int main() {
    using Lattice = D2Q9;

    std::cout << "--- Active Nematic Channel Flow Simulation ---" << std::endl;
    std::cout << "Grid: " << LX << "x" << LY << ", Activity: " << ZETA << std::endl;

    // 1. Create simulation components
    LatticeGrid<Lattice::D, Lattice::Q> grid(LX, LY);
    ELBMSolver<Lattice> solver(VISCOSITY);
    ActiveNematic<Lattice> nematic(LX, LY, 1, ZETA, K_ELASTIC, GAMMA_ROT);

    // 2. Initialization
    // Initialize fluid at rest with uniform density
    for (auto& node : grid) {
        node.fluid.rho = 1.0;
        BGKEquilibrium<Lattice>::compute(node.fluid, node.df.f);
    }
    
    // Initialize nematic field (homeotropic + noise)
    nematic.initialize(1.0, 1e-4);

    // Set top and bottom boundaries as solid walls for bounce-back
    // Note: The BounceBackBC class swaps populations. We need to apply it
    // after streaming to the boundary nodes.
    // The current BoundaryManager is too complex. A simpler approach is to
    // have a list of boundary nodes and apply the BC directly.
    std::vector<Node<Lattice::D, Lattice::Q>*> top_wall_nodes;
    std::vector<Node<Lattice::D, Lattice::Q>*> bottom_wall_nodes;
    for(int x = 0; x < LX; ++x) {
        // Here we mark the fluid nodes adjacent to the wall, not the wall itself.
        // The streaming will send populations to these nodes, which then bounce back.
        // A simpler way for this problem is just to apply bounce back on the f_new
        // populations of the wall-adjacent nodes before they are streamed.
        // However, the standard is to stream and then apply BC.
        // Let's use a simpler direct application loop for now.
    }


    // --- Main Simulation Loop ---
    for (int step = 0; step <= MAX_STEPS; ++step) {
        if (step % OUTPUT_INTERVAL == 0) {
            std::cout << "Step: " << step << "/" << MAX_STEPS << std::endl;
            std::string filename = "../output/active_nematic_step_" + std::to_string(step) + ".dat";
            write_fields(grid, nematic, filename);
        }

        // 1. Compute Active Force from nematic field
        nematic.compute_active_force();
        const auto& force_x = nematic.get_force_x();
        const auto& force_y = nematic.get_force_y();

        // 2. LBM Collision (with active force)
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
        // This must be done after streaming
        for (int x = 0; x < LX; ++x) {
            // Top wall (y = LY - 1)
            BounceBackBC<Lattice>::apply(grid(x, LY - 1));
            // Bottom wall (y = 0)
            BounceBackBC<Lattice>::apply(grid(x, 0));
        }

        // 5. Update Q-Tensor (Beris-Edwards equation)
        const double dt = 0.01;  // Much smaller time step for Q-tensor
        nematic.update_Q_tensor(grid, dt);

        // 6. Apply Q-tensor boundary conditions (homeotropic anchoring)
        // At walls (y=0 and y=LY-1): director along y, so Q_xx = -0.5, Q_xy = 0
        auto& Q_xx = nematic.Q_xx();
        auto& Q_xy = nematic.Q_xy();
        for (int x = 0; x < LX; ++x) {
            // Bottom wall (y=0)
            const size_t idx_bottom = 0 * LX + x;
            Q_xx[idx_bottom] = -0.5;
            Q_xy[idx_bottom] = 0.0;

            // Top wall (y=LY-1)
            const size_t idx_top = (LY-1) * LX + x;
            Q_xx[idx_top] = -0.5;
            Q_xy[idx_top] = 0.0;
        }

    }

    std::cout << "Simulation finished." << std::endl;

    return 0;
}
