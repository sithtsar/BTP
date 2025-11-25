#!/usr/bin/env python3
"""
Quick verification script to show a few frames of the active swarm animation
"""
import numpy as np
import matplotlib.pyplot as plt
import os

def load_swarm_data(filename):
    """Load simulation data from .dat file"""
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()

        # Find where particle data starts
        fluid_end = -1
        for i, line in enumerate(lines):
            if line.startswith('# Particle data'):
                fluid_end = i
                break

        if fluid_end == -1:
            raise ValueError("Could not find particle data section")

        # Load fluid data
        fluid_data = np.loadtxt(lines[1:fluid_end], skiprows=0)
        x_f = fluid_data[:, 0].astype(int)
        y_f = fluid_data[:, 1].astype(int)
        ux = fluid_data[:, 2]
        uy = fluid_data[:, 3]
        rho = fluid_data[:, 4]

        # Get grid dimensions
        nx = len(np.unique(x_f))
        ny = len(np.unique(y_f))

        # Reshape fluid data
        ux_grid = ux.reshape((ny, nx))
        uy_grid = uy.reshape((ny, nx))
        rho_grid = rho.reshape((ny, nx))

        # Load particle data
        particle_data = np.loadtxt(lines[fluid_end+1:])
        if len(particle_data.shape) == 1:  # Single particle
            particle_data = particle_data.reshape(1, -1)

        px = particle_data[:, 0]
        py = particle_data[:, 1]
        pvx = particle_data[:, 2]
        pvy = particle_data[:, 3]
        porient = particle_data[:, 4]

        return ux_grid, uy_grid, rho_grid, px, py, pvx, pvy, porient

    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None

def show_key_frames():
    """Show key frames from the animation to verify it works"""
    steps_to_show = [0, 500, 1000, 1500, 2000]

    fig, axes = plt.subplots(1, len(steps_to_show), figsize=(20, 4))

    im = None
    for i, step in enumerate(steps_to_show):
        filename = f'output/active_swarm_step_{step}.dat'

        if os.path.exists(filename):
            data = load_swarm_data(filename)
            if data is not None:
                ux, uy, rho, px, py, pvx, pvy, porient = data

                # Plot horizontal velocity
                im = axes[i].imshow(ux, origin='lower',
                                   extent=(0, ux.shape[1], 0, ux.shape[0]),
                                   cmap='RdBu_r', vmin=-0.03, vmax=0.03, aspect='auto')

                # Plot particles
                axes[i].scatter(px, py, c='red', s=20, alpha=0.8, edgecolors='black', linewidth=0.5)

                # Add orientation arrows for some particles
                for j in range(0, len(px), 10):  # Every 10th particle
                    dx = 4 * np.cos(porient[j])
                    dy = 4 * np.sin(porient[j])
                    axes[i].arrow(px[j], py[j], dx, dy, head_width=1.5, head_length=2,
                                 fc='yellow', ec='yellow', alpha=0.9, linewidth=1)

                axes[i].set_title(f'Step {step}', fontsize=12)
                axes[i].set_xlabel('X')
                axes[i].set_ylabel('Y')

    # Add colorbar if we have at least one valid frame
    if im is not None:
        cbar = plt.colorbar(im, ax=axes, shrink=0.8, aspect=30)
        cbar.set_label('Horizontal Velocity u_x')

    plt.suptitle('Active Swarm: Evolution of Fluid Perturbations Over Time', fontsize=16)
    plt.tight_layout()
    plt.savefig('animation_key_frames.png', dpi=150, bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    show_key_frames()