#include "active_matter/active_particles.h"
#include <iostream>

namespace elbm {

ActiveSwarm::ActiveSwarm(int nx, int ny, int num_particles, double swim_speed,
                        double tumble_rate, double particle_radius)
    : num_particles_(num_particles), swim_speed_(swim_speed),
      tumble_rate_(tumble_rate), particle_radius_(particle_radius),
      particles_(num_particles) {}

void ActiveSwarm::initialize_random(std::mt19937& gen) {
    std::uniform_real_distribution<> x_dist(0.0, 400.0); // Assuming domain size
    std::uniform_real_distribution<> y_dist(0.0, 100.0);
    std::uniform_real_distribution<> angle_dist(0.0, 2.0 * M_PI);

    for (auto& particle : particles_) {
        particle.position[0] = x_dist(gen);
        particle.position[1] = y_dist(gen);
        particle.orientation = angle_dist(gen);
        particle.update_velocity();
    }
}

void ActiveSwarm::initialize_cluster(double center_x, double center_y,
                                   double cluster_radius, std::mt19937& gen) {
    std::uniform_real_distribution<> radius_dist(0.0, cluster_radius);
    std::uniform_real_distribution<> angle_dist(0.0, 2.0 * M_PI);
    std::uniform_real_distribution<> orient_dist(0.0, 2.0 * M_PI);

    for (auto& particle : particles_) {
        double radius = radius_dist(gen);
        double angle = angle_dist(gen);

        particle.position[0] = center_x + radius * std::cos(angle);
        particle.position[1] = center_y + radius * std::sin(angle);
        particle.orientation = orient_dist(gen);
        particle.update_velocity();
    }
}

void ActiveSwarm::update_particles(double dt, const LatticeGrid<2, 9>& grid, std::mt19937& gen, double swim_speed) {
    for (auto& particle : particles_) {
        // Run-and-tumble dynamics
        if (should_tumble(dt, gen)) {
            particle.orientation = random_orientation(gen);
        }

        // Add fluid velocity to particle motion (advection)
        // TODO: Fix interpolate_velocity - for now use zero fluid velocity
        auto fluid_vel = std::array<double, 2>{0.0, 0.0};

        // Update particle velocity (swimming + fluid advection)
        particle.velocity[0] = swim_speed * std::cos(particle.orientation) + fluid_vel[0];
        particle.velocity[1] = swim_speed * std::sin(particle.orientation) + fluid_vel[1];

        // Update position
        particle.position[0] += particle.velocity[0] * dt;
        particle.position[1] += particle.velocity[1] * dt;
    }
}

void ActiveSwarm::apply_periodic_boundaries(int nx, int ny) {
    for (auto& particle : particles_) {
        // Periodic in x
        if (particle.position[0] < 0) particle.position[0] += nx;
        if (particle.position[0] >= nx) particle.position[0] -= nx;

        // Bounce off top/bottom walls
        if (particle.position[1] < particle.radius) {
            particle.position[1] = particle.radius;
            particle.velocity[1] = -particle.velocity[1]; // Reflect
        }
        if (particle.position[1] >= ny - particle.radius) {
            particle.position[1] = ny - particle.radius;
            particle.velocity[1] = -particle.velocity[1]; // Reflect
        }
    }
}

void ActiveSwarm::compute_fluid_forces(std::vector<double>& force_x, std::vector<double>& force_y,
                                     int nx, int ny, double coupling_strength) {
    // Initialize forces to zero
    force_x.assign(nx * ny, 0.0);
    force_y.assign(nx * ny, 0.0);

    // Distribute particle momentum to nearby fluid nodes
    for (const auto& particle : particles_) {
        int x = static_cast<int>(particle.position[0]);
        int y = static_cast<int>(particle.position[1]);

        // Simple bilinear interpolation to neighboring nodes
        double fx = x;
        double fy = y;
        double wx = 1.0 - (particle.position[0] - fx);
        double wy = 1.0 - (particle.position[1] - fy);

        // Particle force on fluid (momentum transfer)
        double particle_force_x = -coupling_strength * particle.velocity[0];
        double particle_force_y = -coupling_strength * particle.velocity[1];

        // Distribute to 4 nearest nodes
        for (int dy = 0; dy <= 1; ++dy) {
            for (int dx = 0; dx <= 1; ++dx) {
                int nx_pos = (x + dx) % nx;
                int ny_pos = std::clamp(y + dy, 0, ny - 1);

                double weight = (dx == 0 ? wx : 1.0 - wx) * (dy == 0 ? wy : 1.0 - wy);
                size_t idx = nx_pos + ny_pos * nx;

                force_x[idx] += weight * particle_force_x;
                force_y[idx] += weight * particle_force_y;
            }
        }
    }
}

void ActiveSwarm::add_chemotaxis(const std::vector<double>& chemical_field, int nx, int ny,
                               double chemotactic_strength) {
    for (auto& particle : particles_) {
        int x = static_cast<int>(particle.position[0]);
        int y = static_cast<int>(particle.position[1]);

        // Compute chemical gradient at particle position
        // Simple central difference
        int xp = (x + 1) % nx;
        int xm = (x - 1 + nx) % nx;
        int yp = std::min(y + 1, ny - 1);
        int ym = std::max(y - 1, 0);

        double c_center = chemical_field[x + y * nx];
        double c_right = chemical_field[xp + y * nx];
        double c_left = chemical_field[xm + y * nx];
        double c_up = chemical_field[x + yp * nx];
        double c_down = chemical_field[x + ym * nx];

        double dc_dx = 0.5 * (c_right - c_left);
        double dc_dy = 0.5 * (c_up - c_down);

        // Chemotactic torque (bias orientation towards gradient)
        double torque = chemotactic_strength * (dc_dx * std::cos(particle.orientation) +
                                               dc_dy * std::sin(particle.orientation));

        particle.orientation += torque * 0.01; // Small time step for orientation
    }
}

bool ActiveSwarm::should_tumble(double dt, std::mt19937& gen) {
    std::uniform_real_distribution<> dist(0.0, 1.0);
    double prob = tumble_rate_ * dt;
    return dist(gen) < prob;
}

double ActiveSwarm::random_orientation(std::mt19937& gen) {
    std::uniform_real_distribution<> dist(0.0, 2.0 * M_PI);
    return dist(gen);
}

std::array<double, 2> ActiveSwarm::interpolate_velocity(double x, double y,
                                                      const LatticeGrid<2, 9>& grid) const {
    int ix = static_cast<int>(x);
    int iy = static_cast<int>(y);

    // Clamp to grid boundaries
    ix = std::max(0, std::min(ix, grid.nx() - 1));
    iy = std::max(0, std::min(iy, grid.ny() - 1));

    // Simple nearest neighbor interpolation for now
    const auto& node = grid(ix, iy);
    return {node.fluid.u[0], node.fluid.u[1]};
}

} // namespace elbm