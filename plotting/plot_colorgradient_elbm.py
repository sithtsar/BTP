#!/usr/bin/env python3
"""
Color-Gradient ELBM Static Droplet Analysis

Creates comprehensive plots for the Color-Gradient ELBM simulation:
1. Droplet shape evolution (phase field snapshots)
2. Entropy evolution over time
3. Radial profiles of density and phase field
4. Velocity field evolution
5. Interface properties analysis

Usage:
    python plotting/plot_colorgradient_elbm.py [vtk_pattern] [output_dir]

Arguments:
    vtk_pattern: Glob pattern for VTK files (default: "droplet_*.vtk")
    output_dir: Output directory for plots (default: "figures/02_colorgradient")
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
import glob
import re
import sys
from matplotlib.patches import Circle

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
    velocity_start = None

    for i, line in enumerate(lines):
        if 'SCALARS density' in line:
            density_start = i + 2  # Skip LOOKUP_TABLE line
        elif 'SCALARS phase' in line:
            phase_start = i + 2
        elif 'VECTORS velocity' in line:
            velocity_start = i + 1

    # Read density data
    density = []
    if density_start:
        for i in range(density_start, density_start + nx*ny):
            density.append(float(lines[i].strip()))
    density = np.array(density).reshape(ny, nx)

    # Read phase data
    phase = []
    if phase_start:
        for i in range(phase_start, phase_start + nx*ny):
            phase.append(float(lines[i].strip()))
    phase = np.array(phase).reshape(ny, nx)

    # Read velocity data
    velocity = []
    if velocity_start:
        for i in range(velocity_start, velocity_start + nx*ny):
            vx, vy, vz = map(float, lines[i].strip().split())
            velocity.append([vx, vy])
    velocity = np.array(velocity).reshape(ny, nx, 2)

    return {
        'nx': nx, 'ny': ny,
        'density': density,
        'phase': phase,
        'velocity': velocity
    }

def extract_timestep(filename):
    """Extract timestep from filename"""
    match = re.search(r'droplet_(\d+)\.vtk', filename)
    return int(match.group(1)) if match else 0

def calculate_entropy(density, phase):
    """Calculate Boltzmann H-function: H = sum f ln f"""
    # For multiphase systems, compute entropy from total distribution
    # H = sum f ln f (Boltzmann entropy, negative for physical systems)

    # Normalize density to probability distribution
    f_total = density / np.sum(density)

    # Compute Boltzmann entropy: H = sum f ln f
    # Add small epsilon to avoid log(0)
    eps = 1e-12
    H = np.sum(f_total * np.log(np.maximum(f_total, eps)))

    return H

def calculate_radial_profile(data, center_x, center_y, max_radius):
    """Calculate radial profile from center"""
    ny, nx = data.shape
    radii = []
    values = []

    for y in range(ny):
        for x in range(nx):
            r = np.sqrt((x - center_x)**2 + (y - center_y)**2)
            if r <= max_radius:
                radii.append(r)
                values.append(data[y, x])

    # Bin by radius
    r_bins = np.linspace(0, max_radius, 50)
    profile = []

    for i in range(len(r_bins)-1):
        mask = (np.array(radii) >= r_bins[i]) & (np.array(radii) < r_bins[i+1])
        if np.any(mask):
            profile.append(np.mean(np.array(values)[mask]))
        else:
            profile.append(0.0)

    return r_bins[:-1], profile

def plot_droplet_snapshots(vtk_files, output_dir):
    """Plot droplet shape evolution snapshots"""
    fig, axes = plt.subplots(2, 4, figsize=(18, 9))
    axes = axes.flatten()

    # Select representative timesteps
    timesteps = [0, 200, 500, 800, 1000, 1200, 1500, 2000]
    available_files = {extract_timestep(f): f for f in vtk_files}

    im = None
    for i, t in enumerate(timesteps):
        if t in available_files and i < len(axes):
            data = read_vtk_file(available_files[t])

            # Plot phase field
            im = axes[i].imshow(data['phase'], cmap='RdBu_r', vmin=-1, vmax=1,
                              extent=(0, data['nx'], 0, data['ny']), origin='lower')
            axes[i].set_title(f't = {t}', fontsize=12, weight='bold')
            axes[i].set_xlabel('x', fontsize=10)
            axes[i].set_ylabel('y', fontsize=10)

            # Add circle showing expected droplet radius
            center_x, center_y = data['nx']/2, data['ny']/2
            radius = data['ny']/4.0
            circle = Circle((center_x, center_y), radius, fill=False,
                          edgecolor='black', linestyle='--', linewidth=1, alpha=0.7)
            axes[i].add_patch(circle)

            # Only show x/y labels on edge subplots
            if i < 4:  # Top row
                axes[i].set_xlabel('')
            if i % 4 != 0:  # Not left column
                axes[i].set_ylabel('')

    # Colorbar
    if im is not None:
        cbar_ax = fig.add_axes((0.92, 0.15, 0.02, 0.7))
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_label('Phase Field φ', fontsize=12)

    fig.suptitle('Color-Gradient ELBM: Droplet Shape Evolution', fontsize=16, weight='bold', y=0.98)
    plt.subplots_adjust(left=0.05, right=0.9, bottom=0.08, top=0.92, wspace=0.1, hspace=0.2)

    output_path = Path(output_dir) / 'droplet_evolution.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved droplet evolution to: {output_path}")
    plt.close()

def plot_entropy_evolution(vtk_files, output_dir):
    """Plot entropy evolution over time"""
    fig, ax = plt.subplots(figsize=(12, 6))

    timesteps = []
    entropies = []

    for filename in sorted(vtk_files, key=extract_timestep):
        t = extract_timestep(filename)
        data = read_vtk_file(filename)

        # Boltzmann H-function: H = sum f ln f (negative for physical systems)
        entropy = calculate_entropy(data['density'], data['phase'])

        timesteps.append(t)
        entropies.append(entropy)

    # Plot entropy evolution
    ax.plot(timesteps, entropies, 'b-', linewidth=2, marker='o', markersize=4, alpha=0.8)
    ax.set_xlabel('Time Step', fontsize=12)
    ax.set_ylabel('Boltzmann Entropy H (nats)', fontsize=12)
    ax.set_title('H-Theorem Validation: Entropy Evolution', fontsize=14, weight='bold')
    ax.grid(True, alpha=0.3)

    # Add H-theorem annotation
    ax.text(0.02, 0.98, 'H-Theorem: dH/dt ≤ 0\n(Entropy decreases or stays constant)',
            transform=ax.transAxes, fontsize=11, weight='bold',
            verticalalignment='top', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

    # Add annotations for initial and final values
    initial_entropy = entropies[0]
    final_entropy = entropies[-1]
    delta_H = final_entropy - initial_entropy

    ax.text(0.02, 0.05, f'H_initial = {initial_entropy:.4f}\nH_final = {final_entropy:.4f}\nΔH = {delta_H:.6f}',
            transform=ax.transAxes, fontsize=11, weight='bold',
            verticalalignment='bottom', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

    # Ensure y-axis shows negative values properly
    y_min = min(entropies) * 1.05  # 5% margin
    y_max = max(entropies) * 0.95  # 5% margin
    ax.set_ylim(y_min, y_max)

    plt.tight_layout()

    output_path = Path(output_dir) / 'entropy_evolution.png'
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    print(f"Saved entropy evolution to: {output_path}")
    plt.close()

def plot_radial_profiles(vtk_files, output_dir):
    """Plot radial profiles of density and phase field"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # Select a few timesteps
    timesteps = [0, 500, 1000, 2000]
    available_files = {extract_timestep(f): f for f in vtk_files}
    colors = ['blue', 'green', 'orange', 'red']

    center_x, center_y = 32, 32  # For 64x64 grid
    max_radius = 25

    for i, t in enumerate(timesteps):
        if t in available_files:
            data = read_vtk_file(available_files[t])
            color = colors[i]

            # Density profile
            r_dens, dens_profile = calculate_radial_profile(data['density'], center_x, center_y, max_radius)
            ax1.plot(r_dens, dens_profile, color=color, linewidth=2, label=f't = {t}', alpha=0.8)

            # Phase profile
            r_phi, phi_profile = calculate_radial_profile(data['phase'], center_x, center_y, max_radius)
            ax2.plot(r_phi, phi_profile, color=color, linewidth=2, label=f't = {t}', alpha=0.8)

    # Expected droplet radius
    expected_radius = 16  # ny/4 = 64/4 = 16
    for ax in [ax1, ax2]:
        ax.axvline(x=expected_radius, color='black', linestyle='--', alpha=0.7,
                  label=f'Expected R = {expected_radius}')
        ax.grid(True, alpha=0.3)

    ax1.set_xlabel('Radius (lattice units)', fontsize=12)
    ax1.set_ylabel('Density ρ', fontsize=12)
    ax1.set_title('A) Radial Density Profile', fontsize=14, weight='bold')
    ax1.legend(loc='upper right')
    ax1.set_xlim(0, max_radius)
    ax1.set_ylim(0.95, 1.05)

    ax2.set_xlabel('Radius (lattice units)', fontsize=12)
    ax2.set_ylabel('Phase Field φ', fontsize=12)
    ax2.set_title('B) Radial Phase Field Profile', fontsize=14, weight='bold')
    ax2.legend(loc='lower right')
    ax2.set_xlim(0, max_radius)
    ax2.set_ylim(-1.1, 1.1)

    fig.suptitle('Color-Gradient ELBM: Radial Profiles', fontsize=16, weight='bold', y=0.98)
    plt.subplots_adjust(left=0.08, right=0.95, bottom=0.12, top=0.88, wspace=0.25)

    output_path = Path(output_dir) / 'radial_profiles.png'
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    print(f"Saved radial profiles to: {output_path}")
    plt.close()

def plot_velocity_evolution(vtk_files, output_dir):
    """Plot velocity field evolution"""
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes = axes.flatten()

    # Select timesteps
    timesteps = [0, 500, 1000, 1500, 2000]
    available_files = {extract_timestep(f): f for f in vtk_files}

    im = None
    for i, t in enumerate(timesteps[:5]):  # Only 5 subplots
        if t in available_files:
            data = read_vtk_file(available_files[t])

            # Plot velocity magnitude
            vel_mag = np.sqrt(data['velocity'][:,:,0]**2 + data['velocity'][:,:,1]**2)
            im = axes[i].imshow(vel_mag, cmap='viridis', extent=(0, data['nx'], 0, data['ny']),
                              origin='lower', vmin=0, vmax=0.01)  # Small velocities expected
            axes[i].set_title(f't = {t}', fontsize=12, weight='bold')
            axes[i].set_xlabel('x', fontsize=10)
            axes[i].set_ylabel('y', fontsize=10)

            # Only show x/y labels on edge subplots
            if i < 3:  # Top row
                axes[i].set_xlabel('')
            if i % 3 != 0:  # Not left column
                axes[i].set_ylabel('')

    # Colorbar
    if im is not None:
        cbar_ax = fig.add_axes((0.92, 0.15, 0.02, 0.6))
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_label('Velocity Magnitude', fontsize=12)

    # Leave last subplot empty or add summary
    axes[-1].text(0.5, 0.5, 'Static Droplet:\nMinimal Velocities\nExpected',
                  transform=axes[-1].transAxes, fontsize=12, ha='center', va='center',
                  bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.7))
    axes[-1].set_xlim(0, 1)
    axes[-1].set_ylim(0, 1)
    axes[-1].set_xticks([])
    axes[-1].set_yticks([])

    fig.suptitle('Color-Gradient ELBM: Velocity Field Evolution', fontsize=16, weight='bold', y=0.98)
    plt.subplots_adjust(left=0.07, right=0.9, bottom=0.08, top=0.92, wspace=0.15, hspace=0.25)

    output_path = Path(output_dir) / 'velocity_evolution.png'
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    print(f"Saved velocity evolution to: {output_path}")
    plt.close()

def plot_interface_analysis(vtk_files, output_dir):
    """Analyze interface properties"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    timesteps = []
    interface_areas = []
    max_velocities = []

    for filename in sorted(vtk_files, key=extract_timestep):
        t = extract_timestep(filename)
        data = read_vtk_file(filename)

        # Calculate interface area (cells with significant phase gradient)
        grad_phi_x = np.gradient(data['phase'], axis=1)
        grad_phi_y = np.gradient(data['phase'], axis=0)
        grad_mag = np.sqrt(grad_phi_x**2 + grad_phi_y**2)
        interface_area = np.sum(grad_mag > 0.1)

        # Maximum velocity
        vel_mag = np.sqrt(data['velocity'][:,:,0]**2 + data['velocity'][:,:,1]**2)
        max_vel = np.max(vel_mag)

        timesteps.append(t)
        interface_areas.append(interface_area)
        max_velocities.append(max_vel)

    # Interface area evolution
    ax1.plot(timesteps, interface_areas, 'b-', linewidth=2, marker='o', markersize=4, alpha=0.8)
    ax1.set_xlabel('Time Step', fontsize=12)
    ax1.set_ylabel('Interface Area (cells)', fontsize=12)
    ax1.set_title('A) Interface Area Evolution', fontsize=14, weight='bold')
    ax1.grid(True, alpha=0.3)

    # Expected interface area for a circle
    expected_radius = 16
    expected_perimeter = 2 * np.pi * expected_radius
    ax1.axhline(y=expected_perimeter, color='red', linestyle='--', alpha=0.7,
               label=f'Expected: {expected_perimeter:.0f}')
    ax1.legend(loc='upper right')
    ax1.set_ylim(bottom=0)

    # Maximum velocity evolution
    ax2.semilogy(timesteps, max_velocities, 'r-', linewidth=2, marker='s', markersize=4, alpha=0.8)
    ax2.set_xlabel('Time Step', fontsize=12)
    ax2.set_ylabel('Maximum Velocity (log scale)', fontsize=12)
    ax2.set_title('B) Maximum Velocity Evolution', fontsize=14, weight='bold')
    ax2.grid(True, alpha=0.3)

    # Add annotation about spurious currents
    final_max_vel = max_velocities[-1]
    ax2.text(0.05, 0.95, f'Final max velocity: {final_max_vel:.2e}',
            transform=ax2.transAxes, fontsize=11, weight='bold',
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

    fig.suptitle('Color-Gradient ELBM: Interface Analysis', fontsize=16, weight='bold', y=0.98)
    plt.subplots_adjust(left=0.08, right=0.95, bottom=0.12, top=0.88, wspace=0.25)

    output_path = Path(output_dir) / 'interface_analysis.png'
    plt.savefig(output_path, dpi=200, bbox_inches='tight')
    print(f"Saved interface analysis to: {output_path}")
    plt.close()

def create_summary_statistics(vtk_files):
    """Print summary statistics"""
    print("\n" + "="*60)
    print("COLOR-GRADIENT ELBM STATIC DROPLET ANALYSIS")
    print("="*60)

    # Load initial and final data
    initial_file = [f for f in vtk_files if extract_timestep(f) == 0][0]
    final_file = sorted(vtk_files, key=extract_timestep)[-1]

    initial_data = read_vtk_file(initial_file)
    final_data = read_vtk_file(final_file)

    print(f"\nDomain: {initial_data['nx']} × {initial_data['ny']}")
    print(f"Total timesteps: {len(vtk_files)}")
    print(f"Simulation time: {extract_timestep(final_file)} steps")

    # Interface analysis
    grad_phi_x = np.gradient(final_data['phase'], axis=1)
    grad_phi_y = np.gradient(final_data['phase'], axis=0)
    grad_mag = np.sqrt(grad_phi_x**2 + grad_phi_y**2)
    interface_cells = np.sum(grad_mag > 0.1)

    expected_radius = initial_data['ny'] / 4.0
    expected_perimeter = 2 * np.pi * expected_radius

    print("\nInterface Analysis:")
    print(f"  Expected radius: {expected_radius:.1f}")
    print(f"  Expected perimeter: {expected_perimeter:.1f}")
    print(f"  Final interface cells: {interface_cells}")
    print(f"  Interface quality: {interface_cells/expected_perimeter:.2f}")

    # Velocity analysis
    vel_mag = np.sqrt(final_data['velocity'][:,:,0]**2 + final_data['velocity'][:,:,1]**2)
    max_vel = np.max(vel_mag)
    mean_vel = np.mean(vel_mag)

    print("\nVelocity Analysis:")
    print(f"  Maximum velocity: {max_vel:.2e}")
    print(f"  Mean velocity: {mean_vel:.2e}")

    # Phase field analysis
    phi_center = final_data['phase'][32, 32]  # Center of 64x64 grid
    phi_range = np.max(final_data['phase']) - np.min(final_data['phase'])

    print("\nPhase Field Analysis:")
    print(f"  Center phase: {phi_center:.4f}")
    print(f"  Phase range: {phi_range:.4f}")

    print("\n✓ Analysis complete!")
    print("="*60 + "\n")

def main():
    # Parse arguments
    vtk_pattern = "droplet_*.vtk"
    output_dir = "figures/02_colorgradient"

    if len(sys.argv) > 1:
        vtk_pattern = sys.argv[1]
    if len(sys.argv) > 2:
        output_dir = sys.argv[2]

    # Find VTK files
    vtk_files = glob.glob(vtk_pattern)
    if not vtk_files:
        print(f"No VTK files found matching pattern: {vtk_pattern}")
        sys.exit(1)

    vtk_files.sort(key=extract_timestep)
    print(f"Found {len(vtk_files)} VTK files")

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Generate plots
    print("\nGenerating plots...")

    create_summary_statistics(vtk_files)
    plot_droplet_snapshots(vtk_files, output_dir)
    plot_entropy_evolution(vtk_files, output_dir)
    plot_radial_profiles(vtk_files, output_dir)
    plot_velocity_evolution(vtk_files, output_dir)
    plot_interface_analysis(vtk_files, output_dir)

    print("\n✓ All plots generated!")
    print(f"  → Output directory: {output_dir}")
    print(f"  → Generated {len(list(output_path.glob('*.png')))} plot files")

if __name__ == "__main__":
    main()