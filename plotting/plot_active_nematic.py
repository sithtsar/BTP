#!/usr/bin/env python3
"""
Quick visualization of active nematic simulation results
"""
import numpy as np
import matplotlib.pyplot as plt

def load_data(filename):
    """Load simulation data from .dat file"""
    data = np.loadtxt(filename, skiprows=1)
    x = data[:, 0].astype(int)
    y = data[:, 1].astype(int)
    ux = data[:, 2]
    uy = data[:, 3]
    rho = data[:, 4]
    Qxx = data[:, 5]
    Qxy = data[:, 6]

    # Reshape to grid
    nx = len(np.unique(x))
    ny = len(np.unique(y))

    ux_grid = ux.reshape((ny, nx))
    uy_grid = uy.reshape((ny, nx))
    rho_grid = rho.reshape((ny, nx))
    Qxx_grid = Qxx.reshape((ny, nx))
    Qxy_grid = Qxy.reshape((ny, nx))

    return ux_grid, uy_grid, rho_grid, Qxx_grid, Qxy_grid

def plot_results():
    """Plot the simulation results"""
    # Load final state
    ux, uy, rho, Qxx, Qxy = load_data('output/active_nematic_step_1000.dat')

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Velocity magnitude
    vel_mag = np.sqrt(ux**2 + uy**2)
    im1 = axes[0,0].imshow(vel_mag, origin='lower', cmap='viridis')
    axes[0,0].set_title('Velocity Magnitude')
    plt.colorbar(im1, ax=axes[0,0])

    # Horizontal velocity
    im2 = axes[0,1].imshow(ux, origin='lower', cmap='RdBu_r')
    axes[0,1].set_title('Horizontal Velocity u_x')
    plt.colorbar(im2, ax=axes[0,1])

    # Vertical velocity
    im3 = axes[0,2].imshow(uy, origin='lower', cmap='RdBu_r')
    axes[0,2].set_title('Vertical Velocity u_y')
    plt.colorbar(im3, ax=axes[0,2])

    # Density
    im4 = axes[1,0].imshow(rho, origin='lower', cmap='plasma')
    axes[1,0].set_title('Density ρ')
    plt.colorbar(im4, ax=axes[1,0])

    # Q_xx (order parameter)
    im5 = axes[1,1].imshow(Qxx, origin='lower', cmap='RdYlBu_r')
    axes[1,1].set_title('Q_xx (Order Parameter)')
    plt.colorbar(im5, ax=axes[1,1])

    # Q_xy
    im6 = axes[1,2].imshow(Qxy, origin='lower', cmap='RdYlBu_r')
    axes[1,2].set_title('Q_xy')
    plt.colorbar(im6, ax=axes[1,2])

    plt.tight_layout()
    plt.savefig('active_nematic_results.png', dpi=150, bbox_inches='tight')
    plt.show()

    # Plot velocity profiles
    fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    # Horizontal velocity profile at center
    center_y = uy.shape[0] // 2
    ax1.plot(ux[center_y, :], label='u_x')
    ax1.set_xlabel('x position')
    ax1.set_ylabel('Horizontal velocity')
    ax1.set_title('Horizontal Velocity Profile (y = center)')
    ax1.grid(True)

    # Vertical velocity profile
    ax2.plot(uy[:, ux.shape[1]//2], label='u_y')
    ax2.set_xlabel('y position')
    ax2.set_ylabel('Vertical velocity')
    ax2.set_title('Vertical Velocity Profile (x = center)')
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig('active_nematic_profiles.png', dpi=150, bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    plot_results()