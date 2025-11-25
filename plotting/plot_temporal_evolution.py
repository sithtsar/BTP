#!/usr/bin/env python3
"""
Plot temporal evolution of physical quantities for analytical validation cases.

Shows:
1. Couette Flow: Velocity profiles u(y) at multiple timesteps
2. Poiseuille Flow: Velocity profiles u(y) at multiple timesteps
3. Taylor-Green Vortex: Vorticity magnitude field at multiple timesteps
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm

# Snapshot timesteps (must match test_basic_flows.cpp)
TIMESTEPS = [0, 5000, 10000, 15000, 20000]

# Color schemes for temporal evolution
COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
LINESTYLES = ['-', '--', '-.', ':', '-']

def load_snapshot(filepath):
    """Load a snapshot file (x y rho ux uy format)"""
    if not Path(filepath).exists():
        print(f"Warning: {filepath} not found")
        return None
    data = np.loadtxt(filepath, comments='#')
    return {
        'x': data[:, 0],
        'y': data[:, 1],
        'rho': data[:, 2],
        'ux': data[:, 3],
        'uy': data[:, 4]
    }

def get_grid_shape(data):
    """Determine grid dimensions from flat data"""
    x_unique = np.unique(data['x'])
    y_unique = np.unique(data['y'])
    return len(x_unique), len(y_unique)

def reshape_to_grid(data, nx, ny):
    """Reshape flat arrays to 2D grid"""
    return {
        'x': data['x'].reshape(ny, nx),
        'y': data['y'].reshape(ny, nx),
        'rho': data['rho'].reshape(ny, nx),
        'ux': data['ux'].reshape(ny, nx),
        'uy': data['uy'].reshape(ny, nx)
    }

def plot_couette_temporal_evolution(collision='BGK', output_dir='figures/analytical_validation'):
    """Plot Couette flow velocity profiles at multiple timesteps"""

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Load data at each timestep
    for i, t in enumerate(TIMESTEPS):
        filepath = f"output/analytical_validation/couette_{collision}_t{t:05d}.dat"
        data = load_snapshot(filepath)

        if data is None:
            continue

        # Get grid dimensions and extract centerline profile
        nx, ny = get_grid_shape(data)
        grid = reshape_to_grid(data, nx, ny)

        # Extract centerline (middle x position)
        x_mid = nx // 2
        y_profile = grid['y'][:, x_mid]
        ux_profile = grid['ux'][:, x_mid]

        # Analytical solution
        u_exact = 0.1 * (y_profile / (ny - 1))

        # Plot numerical profile
        ax1.plot(ux_profile, y_profile,
                color=COLORS[i], linestyle=LINESTYLES[i], linewidth=2,
                label=f't = {t}', alpha=0.8)

        # Plot error
        error = np.abs(ux_profile - u_exact)
        ax2.plot(error, y_profile,
                color=COLORS[i], linestyle=LINESTYLES[i], linewidth=2,
                label=f't = {t}', alpha=0.8)

    # Plot analytical solution on first panel
    y_analytical = np.linspace(0, ny-1, 100)
    u_analytical = 0.1 * (y_analytical / (ny - 1))
    ax1.plot(u_analytical, y_analytical, 'k-', linewidth=3,
            label='Analytical', alpha=0.5)

    # Format left panel (velocity profiles)
    ax1.set_xlabel('u(y)', fontsize=12)
    ax1.set_ylabel('y', fontsize=12)
    ax1.set_title(f'Couette Flow: {collision}\nVelocity Profile Evolution', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10, framealpha=0.9)
    ax1.grid(True, alpha=0.3)

    # Format right panel (error)
    ax2.set_xlabel('|u - u_exact|', fontsize=12)
    ax2.set_ylabel('y', fontsize=12)
    ax2.set_title('Absolute Error vs Analytical', fontsize=13, fontweight='bold')
    ax2.legend(fontsize=10, framealpha=0.9)
    ax2.grid(True, alpha=0.3)
    ax2.set_xscale('log')

    plt.tight_layout()

    output_path = Path(output_dir) / f'couette_temporal_{collision}.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

def plot_poiseuille_temporal_evolution(collision='BGK', output_dir='figures/analytical_validation'):
    """Plot Poiseuille flow velocity profiles at multiple timesteps"""

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Physical parameters (must match test_basic_flows.cpp)
    tau = 0.8
    nu = (tau - 0.5) / 3.0
    force_x = 1e-5

    # Load data at each timestep
    for i, t in enumerate(TIMESTEPS):
        filepath = f"output/analytical_validation/poiseuille_{collision}_t{t:05d}.dat"
        data = load_snapshot(filepath)

        if data is None:
            continue

        # Get grid dimensions and extract centerline profile
        nx, ny = get_grid_shape(data)
        grid = reshape_to_grid(data, nx, ny)

        # Extract centerline (middle x position)
        x_mid = nx // 2
        y_profile = grid['y'][:, x_mid]
        ux_profile = grid['ux'][:, x_mid]

        # Analytical solution (parabolic profile)
        H = ny - 1.0
        rho = 1.0  # Lattice density
        u_exact = (force_x / (2.0 * rho * nu)) * y_profile * (H - y_profile)

        # Plot numerical profile
        ax1.plot(ux_profile, y_profile,
                color=COLORS[i], linestyle=LINESTYLES[i], linewidth=2,
                label=f't = {t}', alpha=0.8)

        # Plot error
        error = np.abs(ux_profile - u_exact)
        ax2.plot(error, y_profile,
                color=COLORS[i], linestyle=LINESTYLES[i], linewidth=2,
                label=f't = {t}', alpha=0.8)

    # Plot analytical solution on first panel
    y_analytical = np.linspace(0, ny-1, 100)
    u_analytical = (force_x / (2.0 * rho * nu)) * y_analytical * (H - y_analytical)
    ax1.plot(u_analytical, y_analytical, 'k-', linewidth=3,
            label='Analytical', alpha=0.5)

    # Format left panel (velocity profiles)
    ax1.set_xlabel('u(y)', fontsize=12)
    ax1.set_ylabel('y', fontsize=12)
    ax1.set_title(f'Poiseuille Flow: {collision}\nVelocity Profile Evolution', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10, framealpha=0.9)
    ax1.grid(True, alpha=0.3)

    # Format right panel (error)
    ax2.set_xlabel('|u - u_exact|', fontsize=12)
    ax2.set_ylabel('y', fontsize=12)
    ax2.set_title('Absolute Error vs Analytical', fontsize=13, fontweight='bold')
    ax2.legend(fontsize=10, framealpha=0.9)
    ax2.grid(True, alpha=0.3)
    ax2.set_xscale('log')

    plt.tight_layout()

    output_path = Path(output_dir) / f'poiseuille_temporal_{collision}.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

def compute_vorticity(ux, uy, dx=1.0, dy=1.0):
    """Compute vorticity (curl of velocity field) using central differences"""
    # ω_z = ∂v/∂x - ∂u/∂y
    dudy = np.gradient(uy, dy, axis=0)
    dvdx = np.gradient(ux, dx, axis=1)
    return dvdx - dudy

def plot_taylor_green_temporal_evolution(collision='BGK', output_dir='figures/analytical_validation'):
    """Plot Taylor-Green vortex vorticity field at multiple timesteps"""

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Fewer timesteps for TGV (shorter simulation)
    tgv_timesteps = [0, 2500, 5000, 7500, 10000]

    fig = plt.figure(figsize=(16, 6))
    gs = gridspec.GridSpec(1, 5, figure=fig, wspace=0.3)

    for idx, t in enumerate(tgv_timesteps):
        filepath = f"output/analytical_validation/taylor_{collision}_t{t:05d}.dat"
        data = load_snapshot(filepath)

        if data is None:
            continue

        # Get grid dimensions and reshape
        nx, ny = get_grid_shape(data)
        grid = reshape_to_grid(data, nx, ny)

        # Compute vorticity
        vorticity = compute_vorticity(grid['ux'], grid['uy'])

        # Create subplot
        ax = fig.add_subplot(gs[0, idx])

        # Plot vorticity field with diverging colormap
        vmax = np.max(np.abs(vorticity))
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

        im = ax.contourf(grid['x'], grid['y'], vorticity,
                        levels=20, cmap='RdBu_r', norm=norm)

        ax.set_aspect('equal')
        ax.set_title(f't = {t}', fontsize=11, fontweight='bold')
        ax.set_xlabel('x', fontsize=10)
        if idx == 0:
            ax.set_ylabel('y', fontsize=10)
        else:
            ax.set_yticklabels([])

        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Vorticity', fontsize=9)

    fig.suptitle(f'Taylor-Green Vortex: {collision}\nVorticity Field Evolution',
                fontsize=14, fontweight='bold', y=1.02)

    output_path = Path(output_dir) / f'taylor_temporal_{collision}.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()

def plot_all_temporal_evolution():
    """Generate all temporal evolution plots"""

    print("\n=== Generating Temporal Evolution Plots ===\n")

    # Couette Flow
    print("Couette Flow:")
    plot_couette_temporal_evolution('BGK')
    plot_couette_temporal_evolution('ELBM')

    # Poiseuille Flow
    print("\nPoiseuille Flow:")
    plot_poiseuille_temporal_evolution('BGK')
    plot_poiseuille_temporal_evolution('ELBM')

    # Taylor-Green Vortex
    print("\nTaylor-Green Vortex:")
    plot_taylor_green_temporal_evolution('BGK')
    plot_taylor_green_temporal_evolution('ELBM')

    print("\n=== All temporal evolution plots generated ===")

if __name__ == "__main__":
    plot_all_temporal_evolution()
