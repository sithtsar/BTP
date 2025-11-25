#!/usr/bin/env python3
"""
Visualization of active particle swarm simulation results
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
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

def plot_swarm_snapshot(step):
    """Plot a single snapshot of the swarm simulation"""
    filename = f'output/active_swarm_step_{step}.dat'
    if not os.path.exists(filename):
        print(f"File {filename} not found")
        return

    data = load_swarm_data(filename)
    if data is None:
        return

    ux, uy, rho, px, py, pvx, pvy, porient = data

    fig, axes = plt.subplots(2, 2, figsize=(15, 10))

    # Velocity magnitude
    vel_mag = np.sqrt(ux**2 + uy**2)
    im1 = axes[0,0].imshow(vel_mag, origin='lower', extent=[0, ux.shape[1], 0, ux.shape[0]],
                           cmap='viridis', aspect='auto')
    axes[0,0].set_title(f'Fluid Velocity Magnitude (Step {step})')
    plt.colorbar(im1, ax=axes[0,0])

    # Plot particles as arrows showing orientation
    for i in range(len(px)):
        # Particle position
        axes[0,0].plot(px[i], py[i], 'ro', markersize=2, alpha=0.7)
        # Orientation arrow
        dx = 3 * np.cos(porient[i])
        dy = 3 * np.sin(porient[i])
        axes[0,0].arrow(px[i], py[i], dx, dy, head_width=1, head_length=1,
                        fc='red', ec='red', alpha=0.5, linewidth=0.5)

    # Horizontal velocity
    im2 = axes[0,1].imshow(ux, origin='lower', extent=[0, ux.shape[1], 0, ux.shape[0]],
                           cmap='RdBu_r', aspect='auto')
    axes[0,1].set_title('Horizontal Velocity u_x')
    plt.colorbar(im2, ax=axes[0,1])

    # Vertical velocity
    im3 = axes[1,0].imshow(uy, origin='lower', extent=[0, ux.shape[1], 0, ux.shape[0]],
                           cmap='RdBu_r', aspect='auto')
    axes[1,0].set_title('Vertical Velocity u_y')
    plt.colorbar(im3, ax=axes[1,0])

    # Density
    im4 = axes[1,1].imshow(rho, origin='lower', extent=[0, ux.shape[1], 0, ux.shape[0]],
                           cmap='plasma', aspect='auto')
    axes[1,1].set_title('Fluid Density ρ')
    plt.colorbar(im4, ax=axes[1,1])

    # Plot particles on all subplots
    for ax in axes.flat:
        ax.plot(px, py, 'ro', markersize=1, alpha=0.6)

    plt.tight_layout()
    plt.savefig(f'active_swarm_step_{step}.png', dpi=150, bbox_inches='tight')
    plt.show()

def create_animation():
    """Create animation of the swarm simulation showing fluid perturbations"""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    def animate(frame):
        step = frame * 50  # Every 50 steps for smoother animation
        filename = f'output/active_swarm_step_{step}.dat'

        if not os.path.exists(filename):
            return []

        data = load_swarm_data(filename)
        if data is None:
            return []

        ux, uy, rho, px, py, pvx, pvy, porient = data

        # Clear previous plots
        for ax in axes.flat:
            ax.clear()

        # Left plot: Horizontal velocity field showing side-to-side perturbations
        im1 = axes[0].imshow(ux, origin='lower',
                            extent=[0, ux.shape[1], 0, ux.shape[0]],
                            cmap='RdBu_r', vmin=-0.02, vmax=0.02, aspect='auto', animated=True)
        axes[0].set_title(f'Horizontal Fluid Velocity (Step {step})', fontsize=14)
        axes[0].set_xlabel('X position')
        axes[0].set_ylabel('Y position')

        # Add colorbar
        cbar1 = plt.colorbar(im1, ax=axes[0], shrink=0.8)
        cbar1.set_label('Velocity u_x')

        # Plot particles with orientation arrows
        axes[0].scatter(px, py, c='red', s=30, alpha=0.8, edgecolors='black', linewidth=0.5)

        # Add orientation arrows for particles (showing swimming direction)
        for i in range(len(px)):
            if i % 10 == 0:  # Show arrows for every 10th particle to avoid clutter
                dx = 4 * np.cos(porient[i])
                dy = 4 * np.sin(porient[i])
                axes[0].arrow(px[i], py[i], dx, dy, head_width=1.5, head_length=2,
                              fc='yellow', ec='yellow', alpha=0.9, linewidth=1)

        # Right plot: Velocity magnitude with streamlines
        vel_mag = np.sqrt(ux**2 + uy**2)
        im2 = axes[1].imshow(vel_mag, origin='lower',
                            extent=[0, ux.shape[1], 0, ux.shape[0]],
                            cmap='viridis', vmin=0, vmax=0.05, aspect='auto', animated=True)
        axes[1].set_title(f'Fluid Velocity Magnitude (Step {step})', fontsize=14)
        axes[1].set_xlabel('X position')
        axes[1].set_ylabel('Y position')

        # Add colorbar
        cbar2 = plt.colorbar(im2, ax=axes[1], shrink=0.8)
        cbar2.set_label('Velocity magnitude')

        # Plot particles
        axes[1].scatter(px, py, c='red', s=20, alpha=0.7, edgecolors='white', linewidth=0.3)

        # Add some streamlines to show flow direction (sample every 10 points)
        y_sample = np.arange(10, ux.shape[0]-10, 10)
        x_sample = np.arange(10, ux.shape[1]-10, 10)
        axes[1].streamplot(x_sample, y_sample, ux[y_sample[:, None], x_sample],
                          uy[y_sample[:, None], x_sample],
                          density=1.5, color='white', linewidth=0.5, arrowsize=0.8)

        plt.tight_layout()
        return [im1, im2]

    # Create animation with more frames for smoother motion
    frames = 41  # 0, 50, 100, ..., 2000
    anim = animation.FuncAnimation(fig, animate, frames=frames, interval=200, blit=True)

    # Save as GIF
    anim.save('active_swarm_animation.gif', writer='pillow', fps=5, dpi=100)
    print("Animation saved as 'active_swarm_animation.gif'")

    # Also save as MP4 if ffmpeg is available
    try:
        anim.save('active_swarm_animation.mp4', writer='ffmpeg', fps=10, dpi=150)
        print("Animation also saved as 'active_swarm_animation.mp4'")
    except:
        print("MP4 save failed - ffmpeg not available")

    plt.show()

def plot_particle_trajectories():
    """Plot particle trajectories over time"""
    steps = [0, 500, 1000, 1500, 2000]
    colors = ['blue', 'green', 'red', 'orange', 'purple']

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    for i, step in enumerate(steps):
        filename = f'output/active_swarm_step_{step}.dat'
        if os.path.exists(filename):
            data = load_swarm_data(filename)
            if data is not None:
                ux, uy, rho, px, py, pvx, pvy, porient = data
                ax.scatter(px, py, c=colors[i], s=10, alpha=0.7,
                          label=f'Step {step}', edgecolors='black', linewidth=0.5)

    ax.set_xlabel('X position')
    ax.set_ylabel('Y position')
    ax.set_title('Particle Trajectories Over Time')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.savefig('particle_trajectories.png', dpi=150, bbox_inches='tight')
    plt.show()

def create_simple_animation():
    """Create a simple animation focusing on horizontal velocity perturbations"""
    fig, (ax, cax) = plt.subplots(1, 2, figsize=(14, 6),
                                  gridspec_kw={'width_ratios': [1, 0.05]})

    def animate(frame):
        step = frame * 50  # Every 50 steps
        filename = f'output/active_swarm_step_{step}.dat'

        if not os.path.exists(filename):
            return []

        data = load_swarm_data(filename)
        if data is None:
            return []

        ux, uy, rho, px, py, pvx, pvy, porient = data

        ax.clear()

        # Show horizontal velocity field
        im = ax.imshow(ux, origin='lower',
                      extent=(0, ux.shape[1], 0, ux.shape[0]),
                      cmap='RdBu_r', vmin=-0.03, vmax=0.03, aspect='auto', animated=True)

        # Plot particles with swimming direction
        ax.scatter(px, py, c='red', s=40, alpha=0.9, edgecolors='black', linewidth=1)

        # Add orientation arrows (sample every 5th particle)
        for i in range(0, len(px), 5):
            dx = 6 * np.cos(porient[i])
            dy = 6 * np.sin(porient[i])
            ax.arrow(px[i], py[i], dx, dy, head_width=2, head_length=3,
                    fc='yellow', ec='yellow', alpha=1.0, linewidth=1.5)

        ax.set_title(f'Active Swarm: Fluid Horizontal Velocity & Particle Motion (Step {step})', fontsize=16)
        ax.set_xlabel('X position', fontsize=12)
        ax.set_ylabel('Y position', fontsize=12)

        # Add text showing time
        ax.text(0.02, 0.98, f't = {step}', transform=ax.transAxes,
               fontsize=14, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        # Update colorbar (create it only once)
        if frame == 0:
            cbar = plt.colorbar(im, cax=cax)
            cbar.set_label('Horizontal Velocity u_x', fontsize=12)

        return [im]

    # Create animation
    frames = 41  # 0 to 2000 in steps of 50
    anim = animation.FuncAnimation(fig, animate, frames=frames, interval=150, blit=True)

    # Save animation
    anim.save('active_swarm_simple.gif', writer='pillow', fps=6, dpi=120)
    print("Simple animation saved as 'active_swarm_simple.gif'")

    try:
        anim.save('active_swarm_simple.mp4', writer='ffmpeg', fps=10, dpi=150)
        print("Simple animation also saved as 'active_swarm_simple.mp4'")
    except:
        print("MP4 save failed - ffmpeg not available")

    plt.show()

def create_quick_animation():
    """Create a quick animation showing particle motion"""
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    def animate(frame):
        step = frame * 200  # Every 200 steps for faster animation
        filename = f'output/active_swarm_step_{step}.dat'

        if not os.path.exists(filename):
            return []

        data = load_swarm_data(filename)
        if data is None:
            return []

        ux, uy, rho, px, py, pvx, pvy, porient = data

        ax.clear()

        # Plot a subset of the velocity field (to show perturbations)
        # Show only a central region
        x_start, x_end = 150, 250
        y_start, y_end = 30, 70

        ux_subset = ux[y_start:y_end, x_start:x_end]
        extent = (x_start, x_end, y_start, y_end)

        im = ax.imshow(ux_subset, origin='lower', extent=extent,
                      cmap='RdBu_r', vmin=-0.01, vmax=0.01, aspect='auto', animated=True)

        # Plot particles in this region
        mask = (px >= x_start) & (px < x_end) & (py >= y_start) & (py < y_end)
        ax.scatter(px[mask], py[mask], c='red', s=30, alpha=0.8, edgecolors='black', linewidth=1)

        # Add orientation arrows for some particles
        count = 0
        for i in range(len(px)):
            if mask[i] and count < 10:  # Show arrows for up to 10 particles in view
                dx = 2 * np.cos(porient[i])
                dy = 2 * np.sin(porient[i])
                ax.arrow(px[i], py[i], dx, dy, head_width=1, head_length=1.5,
                        fc='yellow', ec='yellow', alpha=0.9, linewidth=1)
                count += 1

        ax.set_title(f'Active Particles & Fluid Perturbations (Step {step})', fontsize=14)
        ax.set_xlabel('X position')
        ax.set_ylabel('Y position')

        # Add colorbar
        if frame == 0:
            cbar = plt.colorbar(im, ax=ax, shrink=0.8)
            cbar.set_label('Horizontal Velocity u_x')

        return [im]

    # Create quick animation
    frames = 11  # 0, 200, 400, ..., 2000
    anim = animation.FuncAnimation(fig, animate, frames=frames, interval=300, blit=True)

    anim.save('active_swarm_quick.gif', writer='pillow', fps=3, dpi=100)
    print("Quick animation saved as 'active_swarm_quick.gif'")

    plt.show()

def show_motion_comparison():
    """Show comparison of particle positions at different times to demonstrate motion"""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    steps = [0, 1000, 2000]
    titles = ['Initial (Step 0)', 'Middle (Step 1000)', 'Final (Step 2000)']

    im = None
    for i, (step, title) in enumerate(zip(steps, titles)):
        filename = f'output/active_swarm_step_{step}.dat'

        if os.path.exists(filename):
            data = load_swarm_data(filename)
            if data is not None:
                ux, uy, rho, px, py, pvx, pvy, porient = data

                # Plot velocity field
                im = axes[i].imshow(ux, origin='lower',
                                   extent=(0, ux.shape[1], 0, ux.shape[0]),
                                   cmap='RdBu_r', vmin=-0.01, vmax=0.01, aspect='auto')

                # Plot particles
                axes[i].scatter(px, py, c='red', s=20, alpha=0.7, edgecolors='black', linewidth=0.5)

                # Add some orientation arrows
                for j in range(0, len(px), 20):  # Every 20th particle
                    dx = 3 * np.cos(porient[j])
                    dy = 3 * np.sin(porient[j])
                    axes[i].arrow(px[j], py[j], dx, dy, head_width=1, head_length=2,
                                 fc='yellow', ec='yellow', alpha=0.8, linewidth=1)

                axes[i].set_title(f'{title}\nParticles Moving!', fontsize=12)
                axes[i].set_xlabel('X position')
                axes[i].set_ylabel('Y position')

    # Add colorbar if we have data
    if im is not None:
        cbar = plt.colorbar(im, ax=axes, shrink=0.8, aspect=30)
        cbar.set_label('Horizontal Velocity u_x')

    plt.suptitle('Active Particle Swarm: Demonstrating Motion Over Time', fontsize=16)
    plt.tight_layout()
    plt.savefig('particle_motion_demo.png', dpi=150, bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    # Plot final state
    plot_swarm_snapshot(2000)

    # Plot trajectories
    plot_particle_trajectories()

    # Show motion comparison
    print("Creating motion comparison plot...")
    show_motion_comparison()