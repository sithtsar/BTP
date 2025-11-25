#ifndef FLUID_STATE_H
#define FLUID_STATE_H

#include <array>
#include <vector>

namespace elbm {

// Macroscopic fluid variables
template<int D>
struct FluidState {
    double rho;                      // Density
    std::array<double, D> u;         // Velocity
    double p;                        // Pressure

    FluidState() : rho(1.0), p(0.0) {
        u.fill(0.0);
    }
};

// Distribution functions for a single node
template<int Q>
struct DistributionFunction {
    std::array<double, Q> f;         // Current distribution
    std::array<double, Q> f_eq;      // Equilibrium distribution
    std::array<double, Q> f_star;    // Auxiliary distribution (for ELBM)
    std::array<double, Q> f_new;     // Post-collision distribution

    DistributionFunction() {
        f.fill(0.0);
        f_eq.fill(0.0);
        f_star.fill(0.0);
        f_new.fill(0.0);
    }
};

// Grid node containing all relevant information
template<int D, int Q>
struct Node {
    DistributionFunction<Q> df;      // Distribution functions
    FluidState<D> fluid;             // Macroscopic variables
    bool is_solid;                   // Boundary flag
    double alpha;                    // Entropic parameter
    double beta;                     // Dissipation parameter

    Node() : is_solid(false), alpha(2.0), beta(1.0) {}
};

// Complete lattice grid
template<int D, int Q>
class LatticeGrid {
public:
    LatticeGrid(int nx, int ny, int nz = 1)
        : nx_(nx), ny_(ny), nz_(nz) {
        size_t total = nx * ny * nz;
        nodes_.resize(total);
    }

    // Access node by indices
    Node<D, Q>& operator()(int x, int y, int z = 0) {
        return nodes_[index(x, y, z)];
    }

    const Node<D, Q>& operator()(int x, int y, int z = 0) const {
        return nodes_[index(x, y, z)];
    }

    // Dimensions
    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int nz() const { return nz_; }
    size_t size() const { return nodes_.size(); }

    // Iterator access
    auto begin() { return nodes_.begin(); }
    auto end() { return nodes_.end(); }
    auto begin() const { return nodes_.begin(); }
    auto end() const { return nodes_.end(); }

private:
    int nx_, ny_, nz_;
    std::vector<Node<D, Q>> nodes_;

    size_t index(int x, int y, int z) const {
        return x + nx_ * (y + ny_ * z);
    }
};

} // namespace elbm

#endif // FLUID_STATE_H
