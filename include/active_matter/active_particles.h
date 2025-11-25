#ifndef ACTIVE_PARTICLES_H
#define ACTIVE_PARTICLES_H

#include <vector>
#include <array>
#include <random>
#include <cmath>
#include <algorithm>
#include "core/lattice.h"
#include "core/fluid_state.h"

namespace elbm {

struct ActiveParticle {
    std::array<double, 2> position;    // x, y coordinates
    std::array<double, 2> velocity;    // vx, vy components
    double orientation;                // angle in radians
    double tumble_rate;                // rate of direction changes
    double swim_speed;                 // swimming speed
    double radius;                     // particle radius for collision/interaction

    ActiveParticle() : position{0.0, 0.0}, velocity{0.0, 0.0}, orientation{0.0},
                        tumble_rate{0.0}, swim_speed{0.0}, radius{1.0} {}

    ActiveParticle(double x, double y, double speed, double tumble, double rad = 1.0)
        : position{x, y}, velocity{0.0, 0.0}, orientation{0.0},
          tumble_rate{tumble}, swim_speed{speed}, radius{rad} {}

    void update_velocity() {
        velocity[0] = swim_speed * std::cos(orientation);
        velocity[1] = swim_speed * std::sin(orientation);
    }
};

class ActiveSwarm {
public:
    ActiveSwarm(int nx, int ny, int num_particles, double swim_speed = 0.1,
                double tumble_rate = 0.1, double particle_radius = 1.0);

    // Initialize particles randomly in domain
    void initialize_random(std::mt19937& gen);

    // Initialize particles in a specific pattern (e.g., center cluster)
    void initialize_cluster(double center_x, double center_y, double cluster_radius, std::mt19937& gen);

    // Update particle positions and orientations (run-and-tumble dynamics)
    void update_particles(double dt, const LatticeGrid<2, 9>& grid, std::mt19937& gen, double swim_speed);

    // Apply periodic boundary conditions
    void apply_periodic_boundaries(int nx, int ny);

    // Compute forces on fluid from particles (momentum coupling)
    void compute_fluid_forces(std::vector<double>& force_x, std::vector<double>& force_y,
                             int nx, int ny, double coupling_strength = 1.0);

    // Get particle data for output/visualization
    const std::vector<ActiveParticle>& get_particles() const { return particles_; }
    std::vector<ActiveParticle>& get_particles() { return particles_; }

    // Add chemotaxis towards a chemical gradient
    void add_chemotaxis(const std::vector<double>& chemical_field, int nx, int ny,
                       double chemotactic_strength = 0.1);

private:
    int num_particles_;
    double swim_speed_;
    double tumble_rate_;
    double particle_radius_;
    std::vector<ActiveParticle> particles_;

    // Helper function to compute tumbling probability
    bool should_tumble(double dt, std::mt19937& gen);

    // Helper function to generate new random orientation
    double random_orientation(std::mt19937& gen);

    // Interpolate fluid velocity at particle position
    std::array<double, 2> interpolate_velocity(double x, double y,
                                              const LatticeGrid<2, 9>& grid) const;
};

} // namespace elbm

#endif // ACTIVE_PARTICLES_H