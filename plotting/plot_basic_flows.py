#!/usr/bin/env python3
"""
Plot results from unified LBM test cases
Compares BGK and ELBM solvers for:
- Couette Flow
- Poiseuille Flow
- Taylor-Green Vortex
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
plt.rcParams['figure.figsize'] = (15, 10)
plt.rcParams['font.size'] = 10

def load_profile_data(filename):
    """Load profile data (y, ux, u_exact, error)"""
    if not Path(filename).exists():
        print(f"Warning: {filename} not found")
        return None
    data = np.loadtxt(filename, skiprows=1)
    return {
        'y': data[:, 0],
        'u_sim': data[:, 1],
        'u_exact': data[:, 2],
        'error': data[:, 3]
    }

def load_field_data(filename):
    """Load full field data (x, y, rho, ux, uy)"""
    if not Path(filename).exists():
        print(f"Warning: {filename} not found")
        return None
    data = np.loadtxt(filename, skiprows=1)
    return {
        'x': data[:, 0],
        'y': data[:, 1],
        'rho': data[:, 2],
        'ux': data[:, 3],
        'uy': data[:, 4]
    }

def plot_couette():
    """Plot Couette flow profiles"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # Load data
    bgk = load_profile_data('output/couette_profile_BGK.dat')
    elbm = load_profile_data('output/couette_profile_ELBM.dat')

    if bgk is None or elbm is None:
        print("Couette data not found")
        return

    # Plot 1: Velocity profiles
    ax = axes[0]
    ax.plot(bgk['u_exact'], bgk['y'], 'k--', linewidth=2.5, label='Analytical', alpha=0.9, zorder=1)
    ax.plot(bgk['u_sim'], bgk['y'], 'bo-', linewidth=1.5, markersize=4, markevery=5,
            label='BGK', alpha=0.7, zorder=3)
    ax.plot(elbm['u_sim'], elbm['y'], 'r^--', linewidth=1.5, markersize=4, markevery=5,
            label='ELBM', alpha=0.7, zorder=2)
    ax.set_xlabel('Velocity u(y)')
    ax.set_ylabel('Height y')
    ax.set_title('Couette Flow: Velocity Profile')
    ax.legend(loc='best', framealpha=0.9)
    ax.grid(True, alpha=0.3)

    # Plot 2: Error comparison
    ax = axes[1]
    ax.semilogy(bgk['y'], bgk['error'], 'b-', linewidth=1.5, label='BGK Error', alpha=0.8)
    ax.semilogy(elbm['y'], elbm['error'], 'r-', linewidth=1.5, label='ELBM Error', alpha=0.8)
    ax.set_xlabel('Height y')
    ax.set_ylabel('Absolute Error')
    ax.set_title('Point-wise Error')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 3: Relative deviation
    ax = axes[2]
    bgk_rel = (bgk['u_sim'] - bgk['u_exact']) / (bgk['u_exact'] + 1e-10) * 100
    elbm_rel = (elbm['u_sim'] - elbm['u_exact']) / (elbm['u_exact'] + 1e-10) * 100
    ax.plot(bgk_rel, bgk['y'], 'b-', linewidth=2, label='BGK', alpha=0.8)
    ax.plot(elbm_rel, elbm['y'], 'r--', linewidth=2, label='ELBM', alpha=0.8)
    ax.axvline(0, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel('Relative Error (%)')
    ax.set_ylabel('Height y')
    ax.set_title('Relative Deviation from Analytical')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('figures/couette_comparison.png', dpi=300, bbox_inches='tight')
    print("Saved: figures/couette_comparison.png")

def plot_poiseuille():
    """Plot Poiseuille flow profiles"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # Load data
    bgk = load_profile_data('output/poiseuille_profile_BGK.dat')
    elbm = load_profile_data('output/poiseuille_profile_ELBM.dat')

    if bgk is None or elbm is None:
        print("Poiseuille data not found")
        return

    # Plot 1: Velocity profiles
    ax = axes[0]
    ax.plot(bgk['u_exact'], bgk['y'], 'k--', linewidth=2.5, label='Analytical', alpha=0.9, zorder=1)
    ax.plot(bgk['u_sim'], bgk['y'], 'bs-', linewidth=1.5, markersize=4, markevery=4,
            label='BGK', alpha=0.7, zorder=3)
    ax.plot(elbm['u_sim'], elbm['y'], 'r^--', linewidth=1.5, markersize=4, markevery=4,
            label='ELBM', alpha=0.7, zorder=2)
    ax.set_xlabel('Velocity u(y)')
    ax.set_ylabel('Height y')
    ax.set_title('Poiseuille Flow: Velocity Profile')
    ax.legend(loc='best', framealpha=0.9)
    ax.grid(True, alpha=0.3)

    # Plot 2: Error comparison
    ax = axes[1]
    ax.semilogy(bgk['y'], bgk['error'], 'b-', linewidth=2, label='BGK Error', alpha=0.8)
    ax.semilogy(elbm['y'], elbm['error'], 'r--', linewidth=2, label='ELBM Error', alpha=0.8)
    ax.set_xlabel('Height y')
    ax.set_ylabel('Absolute Error')
    ax.set_title('Point-wise Error')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 3: Relative deviation
    ax = axes[2]
    # Avoid division by zero at walls
    mask_bgk = np.abs(bgk['u_exact']) > 1e-10
    mask_elbm = np.abs(elbm['u_exact']) > 1e-10

    bgk_rel = np.zeros_like(bgk['u_exact'])
    bgk_rel[mask_bgk] = (bgk['u_sim'][mask_bgk] - bgk['u_exact'][mask_bgk]) / bgk['u_exact'][mask_bgk] * 100

    elbm_rel = np.zeros_like(elbm['u_exact'])
    elbm_rel[mask_elbm] = (elbm['u_sim'][mask_elbm] - elbm['u_exact'][mask_elbm]) / elbm['u_exact'][mask_elbm] * 100

    ax.plot(bgk_rel, bgk['y'], 'b-', linewidth=2, label='BGK', alpha=0.8)
    ax.plot(elbm_rel, elbm['y'], 'r--', linewidth=2, label='ELBM', alpha=0.8)
    ax.axvline(0, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel('Relative Error (%)')
    ax.set_ylabel('Height y')
    ax.set_title('Relative Deviation from Analytical')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('figures/poiseuille_comparison.png', dpi=300, bbox_inches='tight')
    print("Saved: figures/poiseuille_comparison.png")

def plot_taylor_green():
    """Plot Taylor-Green vortex fields"""
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))

    # Load data
    bgk = load_field_data('output/taylor_BGK.dat')
    elbm = load_field_data('output/taylor_ELBM.dat')

    if bgk is None or elbm is None:
        print("Taylor-Green data not found")
        return

    # Reshape data to 2D grid
    nx = int(np.max(bgk['x'])) + 1
    ny = int(np.max(bgk['y'])) + 1

    def reshape_field(data):
        return {
            'ux': data['ux'].reshape(ny, nx),
            'uy': data['uy'].reshape(ny, nx),
            'rho': data['rho'].reshape(ny, nx),
            'speed': np.sqrt(data['ux']**2 + data['uy']**2).reshape(ny, nx)
        }

    bgk_2d = reshape_field(bgk)
    elbm_2d = reshape_field(elbm)

    # BGK plots
    # Velocity magnitude
    im = axes[0, 0].contourf(bgk_2d['speed'], levels=20, cmap='viridis')
    axes[0, 0].set_title('BGK: Velocity Magnitude')
    axes[0, 0].set_xlabel('x')
    axes[0, 0].set_ylabel('y')
    plt.colorbar(im, ax=axes[0, 0])

    # Streamlines
    X, Y = np.meshgrid(range(nx), range(ny))
    axes[0, 1].streamplot(X, Y, bgk_2d['ux'], bgk_2d['uy'],
                          color=bgk_2d['speed'], cmap='plasma', density=1.5)
    axes[0, 1].set_title('BGK: Streamlines')
    axes[0, 1].set_xlabel('x')
    axes[0, 1].set_ylabel('y')

    # Vorticity
    uy_dx = np.gradient(bgk_2d['uy'], axis=1)
    ux_dy = np.gradient(bgk_2d['ux'], axis=0)
    vorticity_bgk = uy_dx - ux_dy

    im = axes[0, 2].contourf(vorticity_bgk, levels=20, cmap='RdBu_r')
    axes[0, 2].set_title('BGK: Vorticity')
    axes[0, 2].set_xlabel('x')
    axes[0, 2].set_ylabel('y')
    plt.colorbar(im, ax=axes[0, 2])

    # ELBM plots
    # Velocity magnitude
    im = axes[1, 0].contourf(elbm_2d['speed'], levels=20, cmap='viridis')
    axes[1, 0].set_title('ELBM: Velocity Magnitude')
    axes[1, 0].set_xlabel('x')
    axes[1, 0].set_ylabel('y')
    plt.colorbar(im, ax=axes[1, 0])

    # Streamlines
    axes[1, 1].streamplot(X, Y, elbm_2d['ux'], elbm_2d['uy'],
                          color=elbm_2d['speed'], cmap='plasma', density=1.5)
    axes[1, 1].set_title('ELBM: Streamlines')
    axes[1, 1].set_xlabel('x')
    axes[1, 1].set_ylabel('y')

    # Vorticity
    uy_dx = np.gradient(elbm_2d['uy'], axis=1)
    ux_dy = np.gradient(elbm_2d['ux'], axis=0)
    vorticity_elbm = uy_dx - ux_dy

    im = axes[1, 2].contourf(vorticity_elbm, levels=20, cmap='RdBu_r')
    axes[1, 2].set_title('ELBM: Vorticity')
    axes[1, 2].set_xlabel('x')
    axes[1, 2].set_ylabel('y')
    plt.colorbar(im, ax=axes[1, 2])

    plt.tight_layout()
    plt.savefig('figures/taylor_green_fields.png', dpi=300, bbox_inches='tight')
    print("Saved: figures/taylor_green_fields.png")

    # Plot difference
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    diff_ux = elbm_2d['ux'] - bgk_2d['ux']
    diff_uy = elbm_2d['uy'] - bgk_2d['uy']
    diff_speed = elbm_2d['speed'] - bgk_2d['speed']

    im = axes[0].contourf(diff_ux, levels=20, cmap='RdBu_r')
    axes[0].set_title('Difference: ELBM - BGK (ux)')
    plt.colorbar(im, ax=axes[0])

    im = axes[1].contourf(diff_uy, levels=20, cmap='RdBu_r')
    axes[1].set_title('Difference: ELBM - BGK (uy)')
    plt.colorbar(im, ax=axes[1])

    im = axes[2].contourf(diff_speed, levels=20, cmap='RdBu_r')
    axes[2].set_title('Difference: ELBM - BGK (speed)')
    plt.colorbar(im, ax=axes[2])

    plt.tight_layout()
    plt.savefig('figures/taylor_green_difference.png', dpi=300, bbox_inches='tight')
    print("Saved: figures/taylor_green_difference.png")

def plot_summary():
    """Create summary comparison plot"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Couette profile
    ax = axes[0, 0]
    bgk = load_profile_data('output/couette_profile_BGK.dat')
    elbm = load_profile_data('output/couette_profile_ELBM.dat')
    if bgk and elbm:
        ax.plot(bgk['u_exact'], bgk['y'], 'k--', linewidth=2.5, label='Analytical', alpha=0.9, zorder=1)
        ax.plot(bgk['u_sim'], bgk['y'], 'bo-', linewidth=1.5, markersize=3, markevery=5,
                label='BGK', alpha=0.7, zorder=3)
        ax.plot(elbm['u_sim'], elbm['y'], 'r^--', linewidth=1.5, markersize=3, markevery=5,
                label='ELBM', alpha=0.7, zorder=2)
        ax.set_xlabel('Velocity u(y)')
        ax.set_ylabel('Height y')
        ax.set_title('Couette Flow')
        ax.legend(framealpha=0.9)
        ax.grid(True, alpha=0.3)

    # Poiseuille profile
    ax = axes[0, 1]
    bgk = load_profile_data('output/poiseuille_profile_BGK.dat')
    elbm = load_profile_data('output/poiseuille_profile_ELBM.dat')
    if bgk and elbm:
        ax.plot(bgk['u_exact'], bgk['y'], 'k--', linewidth=2.5, label='Analytical', alpha=0.9, zorder=1)
        ax.plot(bgk['u_sim'], bgk['y'], 'bs-', linewidth=1.5, markersize=3, markevery=4,
                label='BGK', alpha=0.7, zorder=3)
        ax.plot(elbm['u_sim'], elbm['y'], 'r^--', linewidth=1.5, markersize=3, markevery=4,
                label='ELBM', alpha=0.7, zorder=2)
        ax.set_xlabel('Velocity u(y)')
        ax.set_ylabel('Height y')
        ax.set_title('Poiseuille Flow')
        ax.legend(framealpha=0.9)
        ax.grid(True, alpha=0.3)

    # Taylor-Green velocity field (BGK)
    ax = axes[1, 0]
    bgk = load_field_data('output/taylor_BGK.dat')
    if bgk:
        nx = int(np.max(bgk['x'])) + 1
        ny = int(np.max(bgk['y'])) + 1
        speed = np.sqrt(bgk['ux']**2 + bgk['uy']**2).reshape(ny, nx)
        im = ax.contourf(speed, levels=20, cmap='viridis')
        ax.set_title('Taylor-Green: BGK')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plt.colorbar(im, ax=ax, label='Speed')

    # Taylor-Green velocity field (ELBM)
    ax = axes[1, 1]
    elbm = load_field_data('output/taylor_ELBM.dat')
    if elbm:
        nx = int(np.max(elbm['x'])) + 1
        ny = int(np.max(elbm['y'])) + 1
        speed = np.sqrt(elbm['ux']**2 + elbm['uy']**2).reshape(ny, nx)
        im = ax.contourf(speed, levels=20, cmap='viridis')
        ax.set_title('Taylor-Green: ELBM')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plt.colorbar(im, ax=ax, label='Speed')

    plt.tight_layout()
    plt.savefig('figures/analytical_validation_summary.png', dpi=300, bbox_inches='tight')
    print("Saved: figures/analytical_validation_summary.png")

def main():
    # Create figures directory
    Path('figures').mkdir(exist_ok=True)

    print("=" * 60)
    print("Plotting Analytical Validation Results")
    print("=" * 60)

    # Generate plots
    print("\n1. Plotting Couette flow...")
    plot_couette()

    print("\n2. Plotting Poiseuille flow...")
    plot_poiseuille()

    print("\n3. Plotting Taylor-Green vortex...")
    plot_taylor_green()

    print("\n4. Creating summary plot...")
    plot_summary()

    print("\n" + "=" * 60)
    print("All plots generated successfully!")
    print("=" * 60)

    # Display plots
    if '--show' in sys.argv:
        plt.show()

if __name__ == '__main__':
    main()
