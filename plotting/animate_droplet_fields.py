#!/usr/bin/env python3
"""
Spatial Field Animation for Two-Phase Droplet Simulations

Creates animated visualizations of:
- Density field with interface contour
- Velocity field (quiver plot with magnitude colormap)
- Combined view (density + velocity overlay)

Usage:
    python plotting/animate_droplet_fields.py output/elbm_t*.dat --output figures/elbm_animation.mp4
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
from pathlib import Path
import sys
import glob
import re
import argparse

def load_snapshot(filename):
    """Load a single snapshot .dat file"""
    data = np.loadtxt(filename, skiprows=1)
    x = data[:, 0].astype(int)
    y = data[:, 1].astype(int)
    rho = data[:, 2]
    ux = data[:, 3]
    uy = data[:, 4]
    p = data[:, 5]

    nx = x.max() + 1
    ny = y.max() + 1

    # Reshape into 2D grids
    rho_grid = rho.reshape((ny, nx))
    ux_grid = ux.reshape((ny, nx))
    uy_grid = uy.reshape((ny, nx))
    p_grid = p.reshape((ny, nx))

    return {
        'rho': rho_grid,
        'ux': ux_grid,
        'uy': uy_grid,
        'p': p_grid,
        'nx': nx,
        'ny': ny
    }

def extract_timestep(filename):
    """Extract timestep number from filename like elbm_t250.dat"""
    match = re.search(r'_t(\d+)\.dat', filename)
    return int(match.group(1)) if match else 0

def load_all_snapshots(pattern):
    """Load all snapshots matching the pattern, sorted by timestep"""
    files = sorted(glob.glob(pattern), key=extract_timestep)
    print(f"Found {len(files)} snapshot files")

    snapshots = []
    timesteps = []
    for f in files:
        snapshots.append(load_snapshot(f))
        timesteps.append(extract_timestep(f))

    return snapshots, timesteps

def create_density_animation(snapshots, timesteps, output_file):
    """Create animation of density field with interface contour"""
    print(f"Creating density animation: {output_file}")

    fig, ax = plt.subplots(figsize=(10, 10))

    # Find global min/max for consistent colorbar
    rho_min = min(s['rho'].min() for s in snapshots)
    rho_max = max(s['rho'].max() for s in snapshots)

    def update(frame):
        ax.clear()
        data = snapshots[frame]
        t = timesteps[frame]

        # Plot density field
        im = ax.imshow(data['rho'], origin='lower', cmap='RdBu_r',
                      vmin=rho_min, vmax=rho_max, interpolation='bilinear')

        # Add interface contour (at mid-density)
        rho_mid = 0.5 * (rho_min + rho_max)
        ax.contour(data['rho'], levels=[rho_mid], colors='black',
                  linewidths=2, origin='lower')

        ax.set_title(f'Density Field (t = {t})', fontsize=16, weight='bold')
        ax.set_xlabel('X', fontsize=14)
        ax.set_ylabel('Y', fontsize=14)

        # Colorbar
        if frame == 0:
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label('Density ρ', fontsize=14)

        return [im]

    anim = animation.FuncAnimation(fig, update, frames=len(snapshots),
                                  interval=200, blit=False)

    # Save animation
    Writer = animation.writers['pillow'] if output_file.endswith('.gif') else animation.writers['ffmpeg']
    writer = Writer(fps=5, bitrate=1800)
    anim.save(output_file, writer=writer)
    plt.close()
    print(f"  → Saved to {output_file}")

def create_velocity_animation(snapshots, timesteps, output_file):
    """Create animation of velocity field"""
    print(f"Creating velocity animation: {output_file}")

    fig, ax = plt.subplots(figsize=(10, 10))

    # Find global velocity magnitude range
    vel_max = max(np.sqrt(s['ux']**2 + s['uy']**2).max() for s in snapshots)

    def update(frame):
        ax.clear()
        data = snapshots[frame]
        t = timesteps[frame]

        vel_mag = np.sqrt(data['ux']**2 + data['uy']**2)

        # Plot velocity magnitude as background
        im = ax.imshow(vel_mag, origin='lower', cmap='viridis',
                      vmin=0, vmax=vel_max, interpolation='bilinear')

        # Overlay velocity vectors (subsample for visibility)
        nx, ny = data['nx'], data['ny']
        skip = max(nx // 20, 1)
        X, Y = np.meshgrid(np.arange(0, nx, skip), np.arange(0, ny, skip))
        U = data['ux'][::skip, ::skip]
        V = data['uy'][::skip, ::skip]

        ax.quiver(X, Y, U, V, color='white', alpha=0.7, scale=0.5, width=0.003)

        ax.set_title(f'Velocity Field (t = {t})', fontsize=16, weight='bold')
        ax.set_xlabel('X', fontsize=14)
        ax.set_ylabel('Y', fontsize=14)

        # Colorbar
        if frame == 0:
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label('|u|', fontsize=14)

        return [im]

    anim = animation.FuncAnimation(fig, update, frames=len(snapshots),
                                  interval=200, blit=False)

    # Save animation
    Writer = animation.writers['pillow'] if output_file.endswith('.gif') else animation.writers['ffmpeg']
    writer = Writer(fps=5, bitrate=1800)
    anim.save(output_file, writer=writer)
    plt.close()
    print(f"  → Saved to {output_file}")

def create_combined_animation(snapshots, timesteps, output_file):
    """Create animation with density + velocity overlay"""
    print(f"Creating combined animation: {output_file}")

    fig, ax = plt.subplots(figsize=(12, 10))

    # Find global ranges
    rho_min = min(s['rho'].min() for s in snapshots)
    rho_max = max(s['rho'].max() for s in snapshots)
    vel_max = max(np.sqrt(s['ux']**2 + s['uy']**2).max() for s in snapshots)

    def update(frame):
        ax.clear()
        data = snapshots[frame]
        t = timesteps[frame]

        # Plot density as background
        im = ax.imshow(data['rho'], origin='lower', cmap='RdBu_r',
                      vmin=rho_min, vmax=rho_max, interpolation='bilinear',
                      alpha=0.7)

        # Add interface contour
        rho_mid = 0.5 * (rho_min + rho_max)
        ax.contour(data['rho'], levels=[rho_mid], colors='black',
                  linewidths=3, origin='lower')

        # Overlay velocity vectors
        nx, ny = data['nx'], data['ny']
        skip = max(nx // 15, 1)
        X, Y = np.meshgrid(np.arange(0, nx, skip), np.arange(0, ny, skip))
        U = data['ux'][::skip, ::skip]
        V = data['uy'][::skip, ::skip]

        vel_mag_sub = np.sqrt(U**2 + V**2)
        ax.quiver(X, Y, U, V, vel_mag_sub, cmap='plasma',
                 scale=1.0, width=0.004, alpha=0.8)

        ax.set_title(f'Density + Velocity Fields (t = {t})',
                    fontsize=16, weight='bold')
        ax.set_xlabel('X', fontsize=14)
        ax.set_ylabel('Y', fontsize=14)

        # Add text annotation with max velocity
        max_vel = vel_mag_sub.max() if vel_mag_sub.size > 0 else 0
        ax.text(0.02, 0.98, f'Max |u| = {max_vel:.2e}',
               transform=ax.transAxes, fontsize=12,
               verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        return [im]

    anim = animation.FuncAnimation(fig, update, frames=len(snapshots),
                                  interval=200, blit=False)

    # Save animation
    Writer = animation.writers['pillow'] if output_file.endswith('.gif') else animation.writers['ffmpeg']
    writer = Writer(fps=5, bitrate=1800)
    anim.save(output_file, writer=writer)
    plt.close()
    print(f"  → Saved to {output_file}")

def create_comparison_figure(snapshots, timesteps, output_file, indices=[0, -1]):
    """Create static comparison figure (initial vs final state)"""
    print(f"Creating comparison figure: {output_file}")

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Initial vs Final State Comparison', fontsize=18, weight='bold')

    for row, idx in enumerate(indices):
        data = snapshots[idx]
        t = timesteps[idx]

        # Find ranges
        rho_min, rho_max = data['rho'].min(), data['rho'].max()
        vel_mag = np.sqrt(data['ux']**2 + data['uy']**2)
        vel_max = vel_mag.max()

        # Density
        im1 = axes[row, 0].imshow(data['rho'], origin='lower', cmap='RdBu_r',
                                  vmin=rho_min, vmax=rho_max)
        rho_mid = 0.5 * (rho_min + rho_max)
        axes[row, 0].contour(data['rho'], levels=[rho_mid], colors='black',
                            linewidths=2, origin='lower')
        axes[row, 0].set_title(f'Density (t = {t})', fontsize=14)
        plt.colorbar(im1, ax=axes[row, 0])

        # Velocity magnitude
        im2 = axes[row, 1].imshow(vel_mag, origin='lower', cmap='viridis',
                                  vmin=0, vmax=vel_max)
        axes[row, 1].set_title(f'Velocity Magnitude (t = {t})', fontsize=14)
        plt.colorbar(im2, ax=axes[row, 1])

        # Combined
        axes[row, 2].imshow(data['rho'], origin='lower', cmap='RdBu_r',
                           vmin=rho_min, vmax=rho_max, alpha=0.6)
        axes[row, 2].contour(data['rho'], levels=[rho_mid], colors='black',
                            linewidths=2, origin='lower')
        nx, ny = data['nx'], data['ny']
        skip = max(nx // 12, 1)
        X, Y = np.meshgrid(np.arange(0, nx, skip), np.arange(0, ny, skip))
        U = data['ux'][::skip, ::skip]
        V = data['uy'][::skip, ::skip]
        axes[row, 2].quiver(X, Y, U, V, color='white', alpha=0.8, scale=0.5)
        axes[row, 2].set_title(f'Combined (t = {t})', fontsize=14)

        for ax in axes[row, :]:
            ax.set_xlabel('X')
            ax.set_ylabel('Y')

    plt.tight_layout()
    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  → Saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Animate two-phase droplet spatial fields')
    parser.add_argument('pattern', help='File pattern for snapshots (e.g., output/elbm_t*.dat)')
    parser.add_argument('--output', '-o', default='figures', help='Output directory or file')
    parser.add_argument('--type', choices=['density', 'velocity', 'combined', 'all'], default='all',
                       help='Animation type to create')
    parser.add_argument('--format', choices=['mp4', 'gif'], default='mp4',
                       help='Output format')
    args = parser.parse_args()

    # Load snapshots
    snapshots, timesteps = load_all_snapshots(args.pattern)

    if len(snapshots) == 0:
        print(f"ERROR: No snapshots found matching pattern: {args.pattern}")
        sys.exit(1)

    print(f"Loaded {len(snapshots)} snapshots (timesteps {timesteps[0]} - {timesteps[-1]})")

    # Determine output directory/files
    if args.output.endswith(('.' + args.format)):
        output_file = args.output
        output_dir = str(Path(output_file).parent)
    else:
        output_dir = args.output
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # Extract solver name from pattern
        solver_name = 'droplet'
        if 'bgk' in args.pattern.lower():
            solver_name = 'bgk'
        elif 'elbm' in args.pattern.lower():
            solver_name = 'elbm'

        output_file = str(Path(output_dir) / f'{solver_name}_animation.{args.format}')

    # Create animations
    if args.type in ['density', 'all']:
        density_file = output_file.replace('_animation', '_density')
        create_density_animation(snapshots, timesteps, density_file)

    if args.type in ['velocity', 'all']:
        velocity_file = output_file.replace('_animation', '_velocity')
        create_velocity_animation(snapshots, timesteps, velocity_file)

    if args.type in ['combined', 'all']:
        combined_file = output_file.replace('_animation', '_combined')
        create_combined_animation(snapshots, timesteps, combined_file)

    # Always create static comparison figure
    comparison_file = str(Path(output_dir) / f'{solver_name}_comparison.png').replace('_animation', '')
    create_comparison_figure(snapshots, timesteps, comparison_file)

    print("\n✓ Animation complete!")

if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("Usage: python animate_droplet_fields.py output/elbm_t*.dat --output figures/elbm_animation.mp4")
        print("   or: python animate_droplet_fields.py output/bgk_t*.dat -o figures/")
        print("\nOptions:")
        print("  --type density|velocity|combined|all  (default: all)")
        print("  --format mp4|gif  (default: mp4)")
        sys.exit(0)

    main()
