#!/usr/bin/env python3
"""
Create animated GIF of active nematic spontaneous flow evolution
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import FancyArrow
import glob
import os

def load_data(filename):
    """Load simulation data from .dat file"""
    try:
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

        return ux_grid, uy_grid, rho_grid, Qxx_grid, Qxy_grid, nx, ny
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None

def create_frame(step, fig, ax, data_cache):
    """Create a single frame of the animation"""
    # Clear the entire figure to prevent accumulation
    fig.clear()
    ax = fig.add_subplot(111)

    # Load data for this step
    filename = f"output/active_nematic_step_{step}.dat"
    if filename not in data_cache:
        data = load_data(filename)
        if data is None:
            return
        data_cache[filename] = data

    ux, uy, rho, Qxx, Qxy, nx, ny = data_cache[filename]

    # Calculate velocity magnitude
    vel_mag = np.sqrt(ux**2 + uy**2)

    # Plot velocity magnitude as background
    im = ax.imshow(vel_mag, origin='lower', extent=[0, nx-1, 0, ny-1],
                   cmap='viridis', alpha=0.8, vmin=0, vmax=0.003)

    # Add velocity vectors (downsampled for clarity)
    step_size = 20  # Sample every 20th point
    y_coords, x_coords = np.mgrid[step_size//2:ny:step_size, step_size//2:nx:step_size]

    # Scale vectors for better visibility
    scale_factor = 4000  # Moderate scaling
    ax.quiver(x_coords, y_coords,
              ux[::step_size, ::step_size] * scale_factor,
              uy[::step_size, ::step_size] * scale_factor,
              scale=1, color='white', alpha=0.95, width=0.005,
              headwidth=3, headlength=4, headaxislength=3)

    # Add Q-tensor visualization (director field)
    # Calculate director angle from Q-tensor
    # For 2D nematics: n = (cosθ, sinθ) where θ = 0.5 * atan2(2Qxy, Qxx)
    theta = 0.5 * np.arctan2(2*Qxy, Qxx - Qxy)  # Director angle

    # Plot director field (subsampled but more visible)
    director_step = 12  # Smaller step for more directors
    y_dir, x_dir = np.mgrid[director_step//2:ny:director_step, director_step//2:nx:director_step]

    # Show directors everywhere (not just high order regions)
    Q_magnitude = np.sqrt(Qxx**2 + Qxy**2)

    director_length = 8  # Good visibility
    for i in range(len(y_dir)):
        for j in range(len(x_dir)):
            yi, xj = int(y_dir[i, j]), int(x_dir[i, j])  # Fixed indexing
            angle = theta[yi, xj]
            # Scale director length by order parameter magnitude
            length_scale = max(0.4, min(1.0, Q_magnitude[yi, xj] * 2))
            dx = director_length * length_scale * np.cos(angle)
            dy = director_length * length_scale * np.sin(angle)
            ax.arrow(xj, yi, dx, dy, head_width=2, head_length=3,
                    fc='red', ec='red', alpha=0.9, linewidth=1.5)

    # Styling
    ax.set_xlabel('x position', fontsize=12)
    ax.set_ylabel('y position', fontsize=12)
    ax.set_title(f'Active Nematic Spontaneous Flow\nStep: {step}', fontsize=14, pad=20)

    # Add legend text
    ax.text(0.02, 0.98, 'White arrows: Fluid velocity\nRed bars: Rod directors',
            transform=ax.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Add colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('Velocity Magnitude', fontsize=10)

    # Add boundary lines
    ax.axhline(y=0, color='black', linewidth=2, alpha=0.8)
    ax.axhline(y=ny-1, color='black', linewidth=2, alpha=0.8)

    ax.set_xlim(0, nx-1)
    ax.set_ylim(0, ny-1)

    return ax

def create_animation():
    """Create the animated GIF"""
    # Get all available time steps
    files = glob.glob("output/active_nematic_step_*.dat")
    steps = []
    for f in files:
        step_str = f.replace("output/active_nematic_step_", "").replace(".dat", "")
        steps.append(int(step_str))
    steps.sort()

    print(f"Creating animation with {len(steps)} frames: {steps}")

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8), dpi=100)
    fig.patch.set_facecolor('white')

    # Initialize ax for the first frame
    ax = create_frame(steps[0], fig, ax, {})

    # Cache for loaded data
    data_cache = {}

    # Create animation
    def animate_frame(step_idx):
        step = steps[step_idx]
        new_ax = create_frame(step, fig, ax, data_cache)
        if new_ax is not None:
            return new_ax.get_children()  # Return all artists in the axes
        return []

    anim = animation.FuncAnimation(fig, animate_frame,
                                  frames=len(steps),
                                  interval=800,  # 800ms between frames
                                  blit=False, repeat=True)

    # Save as GIF
    print("Saving animation as GIF...")
    anim.save('active_nematic_evolution.gif',
              writer='pillow', fps=1, dpi=100)

    print("Animation saved as 'active_nematic_evolution.gif'")

    # Also save as MP4 for higher quality
    try:
        print("Saving MP4 version...")
        anim.save('active_nematic_evolution.mp4',
                  writer='ffmpeg', fps=1, dpi=150,
                  extra_args=['-vcodec', 'libx264'])
        print("MP4 saved as 'active_nematic_evolution.mp4'")
    except Exception as e:
        print(f"MP4 creation failed (ffmpeg not available): {e}")

if __name__ == '__main__':
    # Check if output directory exists
    if not os.path.exists('output'):
        print("Error: output directory not found!")
        exit(1)

    # Check if we have data files
    files = glob.glob("output/active_nematic_step_*.dat")
    if len(files) == 0:
        print("Error: No simulation data files found in output/")
        exit(1)

    create_animation()