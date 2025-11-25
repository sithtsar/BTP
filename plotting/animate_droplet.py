#!/usr/bin/env python3
"""
Create an animated GIF of the Color-Gradient ELBM droplet evolution

Usage:
    python plotting/animate_droplet.py [vtk_pattern] [output_file]

Requirements:
    - numpy
    - matplotlib
    - imageio
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
import glob
import re
import sys
from pathlib import Path

# Animation functionality removed - using snapshots instead

def read_vtk_file(filename):
    """Read VTK file and extract data"""
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Parse dimensions
    dim_line = [l for l in lines if l.startswith('DIMENSIONS')][0]
    nx, ny, nz = map(int, dim_line.split()[1:4])

    # Find data sections
    density_start = None
    phase_start = None

    for i, line in enumerate(lines):
        if 'SCALARS density' in line:
            density_start = i + 2
        elif 'SCALARS phase' in line:
            phase_start = i + 2

    # Read phase data
    phase = []
    if phase_start:
        for i in range(phase_start, phase_start + nx*ny):
            phase.append(float(lines[i].strip()))
    phase = np.array(phase).reshape(ny, nx)

    return {
        'nx': nx, 'ny': ny,
        'phase': phase
    }

def extract_timestep(filename):
    """Extract timestep from filename"""
    match = re.search(r'droplet_(\d+)\.vtk', filename)
    return int(match.group(1)) if match else 0

def create_snapshots(vtk_files, output_dir):
    """Create individual snapshots of droplet evolution"""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Sort files by timestep
    vtk_files.sort(key=extract_timestep)

    # Select representative timesteps
    all_timesteps = [extract_timestep(f) for f in vtk_files]
    selected_timesteps = [0, 200, 500, 1000, 1500, 2000]
    available_files = {extract_timestep(f): f for f in vtk_files}

    for t in selected_timesteps:
        if t in available_files:
            filename = available_files[t]
            data = read_vtk_file(filename)

            # Create figure
            fig, ax = plt.subplots(figsize=(8, 8))
            im = ax.imshow(data['phase'], cmap='RdBu_r', vmin=-1, vmax=1,
                          extent=(0, data['nx'], 0, data['ny']), origin='lower')
            ax.set_title(f'Color-Gradient ELBM Droplet\nTime Step: {t}', fontsize=14, weight='bold')
            ax.set_xlabel('x')
            ax.set_ylabel('y')

            # Add colorbar
            cbar = plt.colorbar(im, ax=ax, shrink=0.8)
            cbar.set_label('Phase Field φ')

            # Add expected droplet circle
            center_x, center_y = data['nx']/2, data['ny']/2
            radius = data['ny']/4.0
            circle = Circle((center_x, center_y), radius, fill=False,
                          edgecolor='black', linestyle='--', linewidth=2, alpha=0.7)
            ax.add_patch(circle)
            ax.text(center_x, center_y + radius + 2, f'Expected R = {radius:.1f}',
                   ha='center', va='bottom', fontsize=10, weight='bold')

            # Save figure
            output_file = output_path / f'droplet_snapshot_t{t:04d}.png'
            plt.savefig(output_file, dpi=150, bbox_inches='tight')
            plt.close()

            print(f"Saved snapshot for timestep {t}")

    print(f"✓ Snapshots saved to: {output_path}")

def main():
    vtk_pattern = "droplet_*.vtk"
    output_file = "figures/02_colorgradient/droplet_evolution.gif"

    if len(sys.argv) > 1:
        vtk_pattern = sys.argv[1]
    if len(sys.argv) > 2:
        output_file = sys.argv[2]

    # Find VTK files
    vtk_files = glob.glob(vtk_pattern)
    if not vtk_files:
        print(f"No VTK files found matching pattern: {vtk_pattern}")
        sys.exit(1)

    print(f"Found {len(vtk_files)} VTK files for animation")

    # Create output directory if needed
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Create snapshots
    create_snapshots(vtk_files, Path(output_file).parent)

if __name__ == "__main__":
    main()