#ifndef EQUILIBRIUM_H
#define EQUILIBRIUM_H

#include "core/lattice.h"
#include "core/fluid_state.h"
#include <cmath>

namespace elbm {

// Standard BGK equilibrium distribution (2nd order Hermite expansion)
template<typename Lattice>
class BGKEquilibrium {
public:
    static void compute(const FluidState<Lattice::D>& fluid,
                       std::array<double, Lattice::Q>& f_eq) {
        const double rho = fluid.rho;
        const auto& u = fluid.u;

        // Pre-compute u^2
        double u2 = 0.0;
        for (int d = 0; d < Lattice::D; ++d) {
            u2 += u[d] * u[d];
        }

        // Compute equilibrium for each direction
        for (int i = 0; i < Lattice::Q; ++i) {
            // Compute ci · u
            double ci_dot_u = 0.0;
            if constexpr (Lattice::D == 2) {
                ci_dot_u = Lattice::cx(i) * u[0] + Lattice::cy(i) * u[1];
            } else if constexpr (Lattice::D == 3) {
                ci_dot_u = Lattice::cx(i) * u[0] + Lattice::cy(i) * u[1] +
                          Lattice::cz(i) * u[2];
            }

            const double ci_dot_u2 = ci_dot_u * ci_dot_u;

            // Standard BGK equilibrium
            f_eq[i] = Lattice::weight(i) * rho * (
                1.0 + ci_dot_u / Lattice::cs2 +
                ci_dot_u2 / (2.0 * Lattice::cs2 * Lattice::cs2) -
                u2 / (2.0 * Lattice::cs2)
            );
        }
    }
};

// Entropic equilibrium distribution
template<typename Lattice>
class EntropicEquilibrium {
public:
    static void compute(const FluidState<Lattice::D>& fluid,
                       std::array<double, Lattice::Q>& f_eq) {
        const double rho = fluid.rho;
        const auto& u = fluid.u;

        // For each direction, compute the entropic equilibrium
        // f_i^eq = w_i * rho * ∏(2 - √(1+u_j²)) * ((2u_j + √(1+3u_j²))/(1-u_j))^(c_ij/c)

        // Pre-compute terms for each dimension
        std::array<double, Lattice::D> term1, term2;
        for (int d = 0; d < Lattice::D; ++d) {
            const double uj = u[d];
            const double uj2 = uj * uj;

            term1[d] = 2.0 - std::sqrt(1.0 + uj2);

            // Numerator and denominator for power term
            const double numer = 2.0 * uj + std::sqrt(1.0 + 3.0 * uj2);
            const double denom = 1.0 - uj;
            term2[d] = numer / denom;
        }

        // Compute equilibrium for each direction
        for (int i = 0; i < Lattice::Q; ++i) {
            double prod = rho;

            // Product over dimensions
            for (int d = 0; d < Lattice::D; ++d) {
                prod *= term1[d];
            }

            // Power terms based on velocity components
            if constexpr (Lattice::D == 2) {
                const int cx = Lattice::cx(i);
                const int cy = Lattice::cy(i);
                if (cx != 0) prod *= std::pow(term2[0], cx);
                if (cy != 0) prod *= std::pow(term2[1], cy);
            } else if constexpr (Lattice::D == 3) {
                const int cx = Lattice::cx(i);
                const int cy = Lattice::cy(i);
                const int cz = Lattice::cz(i);
                if (cx != 0) prod *= std::pow(term2[0], cx);
                if (cy != 0) prod *= std::pow(term2[1], cy);
                if (cz != 0) prod *= std::pow(term2[2], cz);
            }

            f_eq[i] = Lattice::weight(i) * prod;
        }
    }
};

// Compute macroscopic variables from distribution functions
template<typename Lattice>
void computeMacroscopic(const std::array<double, Lattice::Q>& f,
                       FluidState<Lattice::D>& fluid) {
    // Compute density
    fluid.rho = 0.0;
    for (int i = 0; i < Lattice::Q; ++i) {
        fluid.rho += f[i];
    }

    // Compute velocity
    fluid.u.fill(0.0);
    for (int i = 0; i < Lattice::Q; ++i) {
        if constexpr (Lattice::D == 2) {
            fluid.u[0] += f[i] * Lattice::cx(i);
            fluid.u[1] += f[i] * Lattice::cy(i);
        } else if constexpr (Lattice::D == 3) {
            fluid.u[0] += f[i] * Lattice::cx(i);
            fluid.u[1] += f[i] * Lattice::cy(i);
            fluid.u[2] += f[i] * Lattice::cz(i);
        }
    }

    // Normalize velocity by density
    if (fluid.rho > 1e-12) {
        for (int d = 0; d < Lattice::D; ++d) {
            fluid.u[d] /= fluid.rho;
        }
    }

    // Compute pressure (isothermal equation of state)
    fluid.p = Lattice::cs2 * fluid.rho;
}

} // namespace elbm

#endif // EQUILIBRIUM_H
