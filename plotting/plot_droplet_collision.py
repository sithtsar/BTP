"""
Droplet Collision Analysis for Color-Gradient ELBM

Analyzes the droplet collision simulation results, checking:
1. Momentum conservation
2. Collision dynamics and coalescence
3. Interface evolution
4. Velocity field evolution
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import glob
import re
from pathlib import Path

def load_vtk_data(filename):
    """Load VTK data from collision simulation."""
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Parse dimensions
    dims_line = [l for l in lines if l.startswith('DIMENSIONS')][0]
    nx, ny, nz = map(int, dims_line.split()[1:4])

    # Find data sections
    data_start = 0
    for i, line in enumerate(lines):
        if line.strip() == 'POINT_DATA {}'.format(nx*ny):
            data_start = i + 1
            break

    # Extract data
    data_lines = lines[data_start:]

    # Parse density
    density_start = data_lines.index('SCALARS density float 1\n') + 2
    density = []
    for i in range(density_start, density_start + nx*ny):
        density.append(float(data_lines[i].strip()))
    density = np.array(density).reshape(ny, nx)

    # Parse phi
    phi_start = data_lines.index('SCALARS phi float 1\n') + 2
    phi = []
    for i in range(phi_start, phi_start + nx*ny):
        phi.append(float(data_lines[i].strip()))
    phi = np.array(phi).reshape(ny, nx)

    # Parse velocity
    vel_start = data_lines.index('VECTORS velocity float\n') + 1
    ux = []
    uy = []
    for i in range(vel_start, vel_start + nx*ny):
        parts = data_lines[i].strip().split()
        ux.append(float(parts[0]))
        uy.append(float(parts[1]))
    ux = np.array(ux).reshape(ny, nx)
    uy = np.array(uy).reshape(ny, nx)

    return density, phi, ux, uy

def compute_momentum(density, ux, uy):
    """Compute total momentum."""
    px = np.sum(density * ux)
    py = np.sum(density * uy)
    return px, py

def analyze_collision():
    """Main analysis function."""
    print("="*60)
    print("DROPLET COLLISION ANALYSIS")
    print("="*60)

    # Find all VTK files
    def extract_number(filename):
        match = re.search(r'(\d+)', filename)
        return int(match.group(1)) if match else 0

    vtk_files = sorted(glob.glob("droplet_collision_*.vtk"), key=extract_number)

    if not vtk_files:
        print("No VTK files found!")
        return

    print(f"Found {len(vtk_files)} VTK files")

    # Load initial state
    density0, phi0, ux0, uy0 = load_vtk_data(vtk_files[0])
    px0, py0 = compute_momentum(density0, ux0, uy0)

    print("Initial state:")
    print(".6f")
    print(".6f")

    # Analyze all timesteps
    timesteps = []
    px_history = []
    py_history = []
    px_error = []
    py_error = []

    snapshots = []

    for i, vtk_file in enumerate(vtk_files):
        match = re.search(r'(\d+)', vtk_file)
        timestep = int(match.group(1)) if match else 0
        timesteps.append(timestep)

        density, phi, ux, uy = load_vtk_data(vtk_file)
        px, py = compute_momentum(density, ux, uy)

        px_history.append(px)
        py_history.append(py)

        # Compute relative errors
        if abs(px0) > 1e-10:
            px_err = abs(px - px0) / abs(px0)
        else:
            px_err = abs(px)
        px_error.append(px_err)

        if abs(py0) > 1e-10:
            py_err = abs(py - py0) / abs(py0)
        else:
            py_err = abs(py)
        py_error.append(py_err)

        # Save snapshots for animation
        if i % 5 == 0:  # Every 5th frame for animation
            snapshots.append((timestep, density, phi, ux, uy))

        if i % 10 == 0:
            print("6d")

    # Create output directory
    output_dir = Path("figures/03_droplet_collision")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Plot momentum conservation
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))

    # Momentum components
    ax1.plot(timesteps, px_history, 'b-', linewidth=2, label='p_x')
    ax1.plot(timesteps, py_history, 'r-', linewidth=2, label='p_y')
    ax1.axhline(px0, color='b', linestyle='--', alpha=0.7, label='.4f')
    ax1.axhline(py0, color='r', linestyle='--', alpha=0.7, label='.4f')
    ax1.set_xlabel('Timestep')
    ax1.set_ylabel('Momentum')
    ax1.set_title('Momentum Components', fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Momentum errors
    ax2.semilogy(timesteps, px_error, 'b-', linewidth=2, label='p_x error')
    ax2.semilogy(timesteps, py_error, 'r-', linewidth=2, label='p_y error')
    ax2.axhline(1e-6, color='k', linestyle='--', linewidth=1, label='Target: < 1e-6')
    ax2.set_xlabel('Timestep')
    ax2.set_ylabel('Relative Momentum Error')
    ax2.set_title('Momentum Conservation Error', fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # Initial density field
    im1 = ax3.imshow(density0.T, origin='lower', cmap='viridis', extent=[0, density0.shape[1], 0, density0.shape[0]])
    ax3.set_title('Initial Density Field', fontweight='bold')
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    plt.colorbar(im1, ax=ax3, label='ρ')

    # Final density field
    density_final, _, _, _ = load_vtk_data(vtk_files[-1])
    im2 = ax4.imshow(density_final.T, origin='lower', cmap='viridis', extent=[0, density_final.shape[1], 0, density_final.shape[0]])
    ax4.set_title('Final Density Field', fontweight='bold')
    ax4.set_xlabel('x')
    ax4.set_ylabel('y')
    plt.colorbar(im2, ax=ax4, label='ρ')

    plt.tight_layout()
    plt.savefig(output_dir / 'droplet_collision_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Create collision sequence animation
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))

    def animate(frame_idx):
        timestep, density, phi, ux, uy = snapshots[frame_idx]

        # Clear axes
        for ax in [ax1, ax2, ax3, ax4]:
            ax.clear()

        # Density
        im1 = ax1.imshow(density.T, origin='lower', cmap='viridis',
                        extent=[0, density.shape[1], 0, density.shape[0]], vmin=0.8, vmax=1.2)
        ax1.set_title('.1f', fontweight='bold')
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')

        # Phase field
        im2 = ax2.imshow(phi.T, origin='lower', cmap='RdBu_r',
                        extent=[0, phi.shape[1], 0, phi.shape[0]], vmin=-1, vmax=1)
        ax2.set_title('Phase Field φ', fontweight='bold')
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')

        # Velocity magnitude
        vel_mag = np.sqrt(ux**2 + uy**2)
        im3 = ax3.imshow(vel_mag.T, origin='lower', cmap='plasma',
                        extent=[0, vel_mag.shape[1], 0, vel_mag.shape[0]])
        ax3.set_title('Velocity Magnitude', fontweight='bold')
        ax3.set_xlabel('x')
        ax3.set_ylabel('y')

        # Velocity vectors (subsampled)
        stride = 4
        x_vec = np.arange(0, ux.shape[1], stride)
        y_vec = np.arange(0, ux.shape[0], stride)
        X_vec, Y_vec = np.meshgrid(x_vec, y_vec)
        U_vec = ux[::stride, ::stride]
        V_vec = uy[::stride, ::stride]

        ax4.quiver(X_vec, Y_vec, U_vec, V_vec, scale=10, alpha=0.7)
        ax4.set_title('Velocity Field', fontweight='bold')
        ax4.set_xlabel('x')
        ax4.set_ylabel('y')
        ax4.set_xlim(0, ux.shape[1])
        ax4.set_ylim(0, ux.shape[0])

        return [ax1, ax2, ax3, ax4]

    # Create animation
    anim = FuncAnimation(fig, animate, frames=len(snapshots), interval=200, blit=False)
    anim.save(output_dir / 'droplet_collision_animation.gif', writer='pillow', fps=5, dpi=100)
    plt.close()

    # Print summary
    final_px_err = px_error[-1]
    final_py_err = py_error[-1]
    max_px_err = max(px_error)
    max_py_err = max(py_error)

    print("\n" + "="*60)
    print("COLLISION ANALYSIS SUMMARY")
    print("="*60)

    print("\nMomentum Conservation:")
    print(f"  Final p_x error: {final_px_err:.6e}")
    print(f"  Final p_y error: {final_py_err:.6e}")
    print(f"  Max p_x error: {max_px_err:.6e}")
    print(f"  Max p_y error: {max_py_err:.6e}")

    print("\nSimulation Parameters:")
    print(f"  Domain: {density0.shape[1]} × {density0.shape[0]}")
    print(f"  Total timesteps: {len(vtk_files)}")
    print(f"  Initial momentum: ({px0:.6f}, {py0:.6f})")

    print("\nValidation Results:")
    if final_px_err < 1e-4 and final_py_err < 1e-4:
        print("  ✓ EXCELLENT: Momentum conserved to < 1e-4")
    elif final_px_err < 1e-2 and final_py_err < 1e-2:
        print("  ✓ GOOD: Momentum conserved to < 1e-2")
    elif final_px_err < 1e-1 and final_py_err < 1e-1:
        print("  ⚠ ACCEPTABLE: Momentum conserved to < 1e-1")
    else:
        print("  ✗ POOR: Momentum conservation failed")

    print(f"\nPlots saved to: {output_dir}/")
    print("="*60)

if __name__ == "__main__":
    analyze_collision()