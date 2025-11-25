#ifndef NEMATIC_TENSOR_H
#define NEMATIC_TENSOR_H

#include "core/lattice.h"
#include "core/fluid_state.h"
#include <vector>
#include <cmath>
#include <random>

namespace elbm {

template<class Lattice>
class ActiveNematic {
public:
    ActiveNematic(int nx, int ny, int nz, double zeta, double k_elastic, double gamma_rot)
        : nx_(nx), ny_(ny), nz_(nz), zeta_(zeta), k_elastic_(k_elastic), gamma_rot_(gamma_rot),
          Q_xx_(nx * ny * nz), Q_xy_(nx * ny * nz),
          force_x_(nx * ny * nz), force_y_(nx * ny * nz)
    {
        // For 3D, we'd need more components
        static_assert(Lattice::D == 2, "ActiveNematic currently only supports 2D");
    }

    void initialize(double initial_s, double noise_strength) {
        std::mt19937 gen(1337); // Fixed seed for reproducibility
        std::uniform_real_distribution<> dis(-noise_strength, noise_strength);

        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const size_t idx = index(x, y, 0);
                // Initialize with homeotropic alignment (director along y) + noise
                // n = (nx, ny) = (small_noise, 1)
                // Q_xx = S * (nx*nx - 1/2) ~= S * (-1/2)
                // Q_xy = S * (nx*ny) ~= S * small_noise
                double director_x = dis(gen);
                double director_y = 1.0;
                double norm = std::sqrt(director_x*director_x + director_y*director_y);
                director_x /= norm;
                director_y /= norm;

                Q_xx_[idx] = initial_s * (director_x * director_x - 0.5);
                Q_xy_[idx] = initial_s * director_x * director_y;
            }
        }
    }

    // Corresponds to F_active = -zeta * div(Q) from the document
    void compute_active_force() {
        #pragma omp parallel for collapse(2)
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const size_t idx = index(x, y, 0);

                // Periodic boundaries in X, solid walls in Y
                const size_t ip = index((x + 1) % nx_, y, 0);
                const size_t im = index((x - 1 + nx_) % nx_, y, 0);
                const size_t jp = index(x, std::min(y + 1, ny_ - 1), 0);
                const size_t jm = index(x, std::max(y - 1, 0), 0);
                
                // Central Difference for divergence
                // div(Q)_x = d(Qxx)/dx + d(Qxy)/dy
                // div(Q)_y = d(Qxy)/dx + d(Qyy)/dy
                // In 2D, Qyy = -Qxx because Q is traceless
                const double dQxx_dx = (Q_xx_[ip] - Q_xx_[im]) * 0.5;
                const double dQxy_dy = (Q_xy_[jp] - Q_xy_[jm]) * 0.5;
                const double dQxy_dx = (Q_xy_[ip] - Q_xy_[im]) * 0.5;
                const double dQyy_dy = -(Q_xx_[jp] - Q_xx_[jm]) * 0.5;

                force_x_[idx] = -zeta_ * (dQxx_dx + dQxy_dy);
                force_y_[idx] = -zeta_ * (dQxy_dx + dQyy_dy);
            }
        }
    }

    // Compute molecular field H = -δF/δQ from Landau-de Gennes free energy
    void compute_molecular_field(std::vector<double>& H_xx, std::vector<double>& H_xy) const {
        H_xx.resize(nx_ * ny_ * nz_);
        H_xy.resize(nx_ * ny_ * nz_);

        #pragma omp parallel for collapse(2)
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const size_t idx = index(x, y, 0);

                // Compute Q invariants
                const double Q_xx = Q_xx_[idx];
                const double Q_xy = Q_xy_[idx];
                const double Q_yy = -Q_xx;  // Traceless condition
                const double Tr_Q2 = Q_xx*Q_xx + 2.0*Q_xy*Q_xy + Q_yy*Q_yy;
                const double Tr_Q3 = Q_xx*(Q_xx*Q_xx + 2.0*Q_xy*Q_xy) + Q_yy*(Q_yy*Q_yy + 2.0*Q_xy*Q_xy);

                // Bulk free energy derivatives (simplified Landau expansion)
                // dF_bulk/dQ_xx = A0*(1-γ/3)*2*Q_xx - A0*γ/3 * 3*Q_xx² + A0*γ/4 * 4*Q_xx*Tr_Q2
                // Using A0=1, γ=0.5 for simplicity (can be made parameters later)
                const double A0 = 1.0;
                const double gamma = 0.5;
                const double dF_dQxx_bulk = A0 * (1.0 - gamma/3.0) * 2.0 * Q_xx
                                          - A0 * gamma/3.0 * 3.0 * Q_xx*Q_xx
                                          + A0 * gamma/4.0 * 4.0 * Q_xx * Tr_Q2;
                const double dF_dQxy_bulk = A0 * (1.0 - gamma/3.0) * 4.0 * Q_xy
                                          - A0 * gamma/3.0 * 6.0 * Q_xx * Q_xy
                                          + A0 * gamma/4.0 * 8.0 * Q_xy * Tr_Q2;

                // Elastic energy derivatives (∇²Q terms)
                // Need second derivatives of Q
                const size_t ip = index((x + 1) % nx_, y, 0);
                const size_t im = index((x - 1 + nx_) % nx_, y, 0);
                const size_t jp = index(x, std::min(y + 1, ny_ - 1), 0);
                const size_t jm = index(x, std::max(y - 1, 0), 0);

                // Laplacian of Q_xx and Q_xy
                const double lap_Qxx = (Q_xx_[ip] - 2.0*Q_xx_[idx] + Q_xx_[im]) +
                                      (Q_xx_[jp] - 2.0*Q_xx_[idx] + Q_xx_[jm]);
                const double lap_Qxy = (Q_xy_[ip] - 2.0*Q_xy_[idx] + Q_xy_[im]) +
                                      (Q_xy_[jp] - 2.0*Q_xy_[idx] + Q_xy_[jm]);

                // Elastic contribution: K * ∇²Q
                const double dF_dQxx_elastic = k_elastic_ * lap_Qxx;
                const double dF_dQxy_elastic = k_elastic_ * lap_Qxy;

                // Total molecular field H = -δF/δQ
                H_xx[idx] = -(dF_dQxx_bulk + dF_dQxx_elastic);
                H_xy[idx] = -(dF_dQxy_bulk + dF_dQxy_elastic);
            }
        }
    }

    // Compute velocity gradients for co-rotation term
    void compute_velocity_gradients(const LatticeGrid<Lattice::D, Lattice::Q>& grid,
                                   std::vector<double>& dux_dx, std::vector<double>& dux_dy,
                                   std::vector<double>& duy_dx, std::vector<double>& duy_dy) const {
        dux_dx.resize(nx_ * ny_ * nz_);
        dux_dy.resize(nx_ * ny_ * nz_);
        duy_dx.resize(nx_ * ny_ * nz_);
        duy_dy.resize(nx_ * ny_ * nz_);

        #pragma omp parallel for collapse(2)
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const size_t idx = index(x, y, 0);

                // Central differences for velocity gradients
                const size_t ip = index((x + 1) % nx_, y, 0);
                const size_t im = index((x - 1 + nx_) % nx_, y, 0);
                const size_t jp = index(x, std::min(y + 1, ny_ - 1), 0);
                const size_t jm = index(x, std::max(y - 1, 0), 0);

                // Convert indices back to coordinates for grid access
                const int x_ip = ip % nx_;
                const int y_ip = ip / nx_;
                const int x_im = im % nx_;
                const int y_im = im / nx_;
                const int x_jp = jp % nx_;
                const int y_jp = jp / nx_;
                const int x_jm = jm % nx_;
                const int y_jm = jm / nx_;

                const double ux_ip = grid(x_ip, y_ip).fluid.u[0];
                const double ux_im = grid(x_im, y_im).fluid.u[0];
                const double ux_jp = grid(x_jp, y_jp).fluid.u[0];
                const double ux_jm = grid(x_jm, y_jm).fluid.u[0];

                const double uy_ip = grid(x_ip, y_ip).fluid.u[1];
                const double uy_im = grid(x_im, y_im).fluid.u[1];
                const double uy_jp = grid(x_jp, y_jp).fluid.u[1];
                const double uy_jm = grid(x_jm, y_jm).fluid.u[1];

                dux_dx[idx] = 0.5 * (ux_ip - ux_im);
                dux_dy[idx] = 0.5 * (ux_jp - ux_jm);
                duy_dx[idx] = 0.5 * (uy_ip - uy_im);
                duy_dy[idx] = 0.5 * (uy_jp - uy_jm);
            }
        }
    }

    // Beris-Edwards Q-tensor evolution equation
    void update_Q_tensor(const LatticeGrid<Lattice::D, Lattice::Q>& grid, double dt, double xi = 0.7) {
        // Temporary arrays for new Q values
        std::vector<double> Q_xx_new = Q_xx_;
        std::vector<double> Q_xy_new = Q_xy_;

        // Compute molecular field H
        std::vector<double> H_xx, H_xy;
        compute_molecular_field(H_xx, H_xy);

        // Compute velocity gradients
        std::vector<double> dux_dx, dux_dy, duy_dx, duy_dy;
        compute_velocity_gradients(grid, dux_dx, dux_dy, duy_dx, duy_dy);

        #pragma omp parallel for collapse(2)
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const size_t idx = index(x, y, 0);

                // Current Q values
                const double Q_xx = Q_xx_[idx];
                const double Q_xy = Q_xy_[idx];
                const double Q_yy = -Q_xx;  // Traceless

                // Velocity and its gradients at this point
                const double ux = grid(x, y).fluid.u[0];
                const double uy = grid(x, y).fluid.u[1];
                const double duxdx = dux_dx[idx];
                const double duxdy = dux_dy[idx];
                const double duydx = duy_dx[idx];
                const double duydy = duy_dy[idx];

                // Strain rate tensor E_ij = 0.5*(∂u_i/∂x_j + ∂u_j/∂x_i)
                const double E_xx = duxdx;
                const double E_xy = 0.5 * (duxdy + duydx);
                const double E_yy = duydy;

                // Vorticity tensor Ω_ij = 0.5*(∂u_i/∂x_j - ∂u_j/∂x_i)
                const double Omega_xy = 0.5 * (duxdy - duydx);

                // Advection term: u · ∇Q
                const size_t ip = index((x + 1) % nx_, y, 0);
                const size_t im = index((x - 1 + nx_) % nx_, y, 0);
                const size_t jp = index(x, std::min(y + 1, ny_ - 1), 0);
                const size_t jm = index(x, std::max(y - 1, 0), 0);

                const double dQxx_dx = 0.5 * (Q_xx_[ip] - Q_xx_[im]);
                const double dQxx_dy = 0.5 * (Q_xx_[jp] - Q_xx_[jm]);
                const double dQxy_dx = 0.5 * (Q_xy_[ip] - Q_xy_[im]);
                const double dQxy_dy = 0.5 * (Q_xy_[jp] - Q_xy_[jm]);

                const double adv_Qxx = ux * dQxx_dx + uy * dQxx_dy;
                const double adv_Qxy = ux * dQxy_dx + uy * dQxy_dy;

                // Co-rotation term: S(Q, u) where S = (ξE + Ω)(Q + I/2) + (Q + I/2)(ξE - Ω) - 2ξ(Q + I/2)Tr(Q E)
                // For 2D, I/2 = diag(0.5, 0.5), but since Q is traceless, Q + I/d = Q + diag(0.5, 0.5)
                const double Q_plus_xx = Q_xx + 0.5;
                const double Q_plus_xy = Q_xy;
                const double Q_plus_yx = Q_xy;  // Symmetric
                const double Q_plus_yy = Q_yy + 0.5;

                // ξE + Ω and ξE - Ω matrices
                const double A_xx = xi * E_xx;
                const double A_xy = xi * E_xy + Omega_xy;
                const double A_yx = xi * E_xy - Omega_xy;
                const double A_yy = xi * E_yy;

                // (ξE + Ω)(Q + I/2)
                const double temp1_xx = A_xx * Q_plus_xx + A_xy * Q_plus_yx;
                const double temp1_xy = A_xx * Q_plus_xy + A_xy * Q_plus_yy;

                // (Q + I/2)(ξE - Ω)
                const double temp2_xx = Q_plus_xx * A_xx + Q_plus_yx * A_yx;
                const double temp2_xy = Q_plus_xx * A_xy + Q_plus_yx * A_yy;

                // S = temp1 + temp2 - 2ξ(Q + I/2)Tr(Q E)
                const double Tr_Q_E = Q_xx * E_xx + 2.0 * Q_xy * E_xy + Q_yy * E_yy;
                const double correction = 2.0 * xi * Tr_Q_E;

                const double S_xx = temp1_xx + temp2_xx - correction * Q_plus_xx;
                const double S_xy = temp1_xy + temp2_xy - correction * Q_plus_xy;

                // Co-rotation contribution
                const double corot_Qxx = S_xx;
                const double corot_Qxy = S_xy;

                // Relaxation term: Γ H
                const double relax_Qxx = gamma_rot_ * H_xx[idx];
                const double relax_Qxy = gamma_rot_ * H_xy[idx];

                // Total derivative: dQ/dt = -Advection + Co-rotation + Relaxation
                const double dQxx_dt = -adv_Qxx + corot_Qxx + relax_Qxx;
                const double dQxy_dt = -adv_Qxy + corot_Qxy + relax_Qxy;

                // Euler integration
                Q_xx_new[idx] = Q_xx + dt * dQxx_dt;
                Q_xy_new[idx] = Q_xy + dt * dQxy_dt;
            }
        }

        // Update Q tensors
        Q_xx_ = Q_xx_new;
        Q_xy_ = Q_xy_new;
    }

    const std::vector<double>& get_force_x() const { return force_x_; }
    const std::vector<double>& get_force_y() const { return force_y_; }
    const std::vector<double>& get_Q_xx() const { return Q_xx_; }
    const std::vector<double>& get_Q_xy() const { return Q_xy_; }

    // Non-const access for boundary conditions
    std::vector<double>& Q_xx() { return Q_xx_; }
    std::vector<double>& Q_xy() { return Q_xy_; }

private:
    int nx_, ny_, nz_;
    double zeta_;        // Activity
    double k_elastic_;   // Elastic constant
    double gamma_rot_;   // Rotational diffusivity

    std::vector<double> Q_xx_, Q_xy_;
    std::vector<double> force_x_, force_y_;

    size_t index(int x, int y, int z) const {
        return x + nx_ * (y + ny_ * z);
    }
};

} // namespace elbm

#endif // NEMATIC_TENSOR_H
