#!/usr/bin/env python3
"""
Visualization for Extended Analytical Test Cases
Plots results from Womersley, Hagen-Poiseuille, Stokes, and Kolmogorov flows
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
plt.rcParams['figure.figsize'] = (12, 10)
plt.rcParams['font.size'] = 10

def plot_womersley(filename='output/womersley_profile.dat', save_dir='figures/extended_analytical'):
    """Plot Womersley flow profile"""
    if not Path(filename).exists():
        print(f"Warning: {filename} not found")
        return

    data = np.loadtxt(filename)
    y = data[:, 0]
    u_x = data[:, 1]

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    ax.plot(u_x, y, 'b-', linewidth=2, label='LBM simulation')
    ax.set_xlabel('Velocity u(y)')
    ax.set_ylabel('y position')
    ax.set_title('Womersley Flow: Oscillatory Pressure-Driven Flow')
    ax.grid(True, alpha=0.3)
    ax.legend()

    Path(save_dir).mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(f'{save_dir}/womersley_profile.png', dpi=150, bbox_inches='tight')
    print(f"Saved: {save_dir}/womersley_profile.png")
    plt.close()

def plot_hagen_poiseuille(filename='output/hagen_poiseuille_radial.dat', save_dir='figures/extended_analytical'):
    """Plot Hagen-Poiseuille (circular pipe) radial profile"""
    if not Path(filename).exists():
        print(f"Warning: {filename} not found")
        return

    data = np.loadtxt(filename)
    r = data[:, 0]
    u_analytical = data[:, 1]
    u_simulated = data[:, 2]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Profile comparison
    ax1.plot(r, u_analytical, 'r--', linewidth=2, label='Analytical', alpha=0.7)
    ax1.plot(r, u_simulated, 'b-', linewidth=2, label='LBM simulation')
    ax1.set_xlabel('Radial distance r')
    ax1.set_ylabel('Velocity u(r)')
    ax1.set_title('Hagen-Poiseuille Flow: Velocity Profile')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Error plot
    error = np.abs(u_simulated - u_analytical)
    ax2.plot(r, error, 'g-', linewidth=2)
    ax2.set_xlabel('Radial distance r')
    ax2.set_ylabel('Absolute Error |u_sim - u_exact|')
    ax2.set_title('Point-wise Error')
    ax2.grid(True, alpha=0.3)

    l2_error = np.sqrt(np.mean(error**2))
    ax2.text(0.05, 0.95, f'L2 error: {l2_error:.6f}',
             transform=ax2.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    Path(save_dir).mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(f'{save_dir}/hagen_poiseuille_radial.png', dpi=150, bbox_inches='tight')
    print(f"Saved: {save_dir}/hagen_poiseuille_radial.png")
    plt.close()

def plot_stokes_shear(filename='output/stokes_shear_profile.dat', save_dir='figures/extended_analytical'):
    """Plot Stokes shear flow profile"""
    if not Path(filename).exists():
        print(f"Warning: {filename} not found")
        return

    data = np.loadtxt(filename)
    y = data[:, 0]
    u_x = data[:, 1]

    # Analytical solution (linear)
    shear_rate = (u_x[-1] - u_x[0]) / (y[-1] - y[0])
    u_analytical = shear_rate * y

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Profile comparison
    ax1.plot(y, u_analytical, 'r--', linewidth=2, label='Analytical (linear)', alpha=0.7)
    ax1.plot(y, u_x, 'b-', linewidth=2, label='LBM simulation')
    ax1.set_xlabel('y position')
    ax1.set_ylabel('Velocity u(y)')
    ax1.set_title('Stokes Shear Flow: Linear Velocity Profile')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Shear rate verification
    du_dy = np.gradient(u_x, y)
    ax2.plot(y[1:-1], du_dy[1:-1], 'g-', linewidth=2, label='Numerical du/dy')
    ax2.axhline(y=shear_rate, color='r', linestyle='--', linewidth=2, label=f'Theory: {shear_rate:.6f}')
    ax2.set_xlabel('y position')
    ax2.set_ylabel('Shear rate du/dy')
    ax2.set_title('Shear Rate Verification (Should be Constant)')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    Path(save_dir).mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(f'{save_dir}/stokes_shear_profile.png', dpi=150, bbox_inches='tight')
    print(f"Saved: {save_dir}/stokes_shear_profile.png")
    plt.close()

def plot_kolmogorov(filename='output/kolmogorov_profile.dat', save_dir='figures/extended_analytical'):
    """Plot Kolmogorov flow profile"""
    if not Path(filename).exists():
        print(f"Warning: {filename} not found")
        return

    data = np.loadtxt(filename)
    y = data[:, 0]
    u_x = data[:, 1]

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    ax.plot(y, u_x, 'b-', linewidth=2, label='LBM simulation')
    ax.set_xlabel('y position')
    ax.set_ylabel('Velocity u(y)')
    ax.set_title('Kolmogorov Flow: Sinusoidal Forcing Pattern')
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Try to fit sinusoid
    try:
        from scipy.optimize import curve_fit

        def sinusoid(x, A, k, phi, offset):
            return A * np.sin(k * x + phi) + offset

        popt, _ = curve_fit(sinusoid, y, u_x, p0=[np.max(u_x), 2*np.pi/len(y), 0, np.mean(u_x)])
        u_fit = sinusoid(y, *popt)

        ax.plot(y, u_fit, 'r--', linewidth=2, alpha=0.7, label=f'Fit: A={popt[0]:.6f}, k={popt[1]:.6f}')
        ax.legend()
    except:
        pass

    Path(save_dir).mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(f'{save_dir}/kolmogorov_profile.png', dpi=150, bbox_inches='tight')
    print(f"Saved: {save_dir}/kolmogorov_profile.png")
    plt.close()

def create_summary_figure(save_dir='figures/extended_analytical'):
    """Create a 2x2 summary figure with all cases"""
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

    # Womersley
    if Path('output/womersley_profile.dat').exists():
        ax1 = fig.add_subplot(gs[0, 0])
        data = np.loadtxt('output/womersley_profile.dat')
        y, u_x = data[:, 0], data[:, 1]
        ax1.plot(u_x, y, 'b-', linewidth=2)
        ax1.set_xlabel('Velocity')
        ax1.set_ylabel('y position')
        ax1.set_title('(a) Womersley Flow\n(Oscillatory)')
        ax1.grid(True, alpha=0.3)

    # Hagen-Poiseuille
    if Path('output/hagen_poiseuille_radial.dat').exists():
        ax2 = fig.add_subplot(gs[0, 1])
        data = np.loadtxt('output/hagen_poiseuille_radial.dat')
        r, u_analytical, u_simulated = data[:, 0], data[:, 1], data[:, 2]
        ax2.plot(r, u_analytical, 'r--', linewidth=2, label='Analytical', alpha=0.7)
        ax2.plot(r, u_simulated, 'b-', linewidth=2, label='Simulation')
        ax2.set_xlabel('Radial distance r')
        ax2.set_ylabel('Velocity')
        ax2.set_title('(b) Hagen-Poiseuille\n(Circular Pipe)')
        ax2.grid(True, alpha=0.3)
        ax2.legend()

    # Stokes Shear
    if Path('output/stokes_shear_profile.dat').exists():
        ax3 = fig.add_subplot(gs[1, 0])
        data = np.loadtxt('output/stokes_shear_profile.dat')
        y, u_x = data[:, 0], data[:, 1]
        shear_rate = (u_x[-1] - u_x[0]) / (y[-1] - y[0])
        u_analytical = shear_rate * y
        ax3.plot(y, u_analytical, 'r--', linewidth=2, label='Analytical', alpha=0.7)
        ax3.plot(y, u_x, 'b-', linewidth=2, label='Simulation')
        ax3.set_xlabel('y position')
        ax3.set_ylabel('Velocity')
        ax3.set_title('(c) Stokes Shear Flow\n(Creeping Flow)')
        ax3.grid(True, alpha=0.3)
        ax3.legend()

    # Kolmogorov
    if Path('output/kolmogorov_profile.dat').exists():
        ax4 = fig.add_subplot(gs[1, 1])
        data = np.loadtxt('output/kolmogorov_profile.dat')
        y, u_x = data[:, 0], data[:, 1]
        ax4.plot(y, u_x, 'b-', linewidth=2)
        ax4.set_xlabel('y position')
        ax4.set_ylabel('Velocity')
        ax4.set_title('(d) Kolmogorov Flow\n(Turbulence Transition)')
        ax4.grid(True, alpha=0.3)

    fig.suptitle('Extended Analytical Validation Cases', fontsize=16, fontweight='bold')

    Path(save_dir).mkdir(parents=True, exist_ok=True)
    plt.savefig(f'{save_dir}/extended_analytical_summary.png', dpi=150, bbox_inches='tight')
    print(f"Saved: {save_dir}/extended_analytical_summary.png")
    plt.close()

def main():
    print("=" * 60)
    print("Extended Analytical Cases Visualization")
    print("=" * 60)

    save_dir = 'figures/extended_analytical'

    print("\nGenerating individual plots...")
    plot_womersley(save_dir=save_dir)
    plot_hagen_poiseuille(save_dir=save_dir)
    plot_stokes_shear(save_dir=save_dir)
    plot_kolmogorov(save_dir=save_dir)

    print("\nGenerating summary figure...")
    create_summary_figure(save_dir=save_dir)

    print("\n" + "=" * 60)
    print("Visualization complete!")
    print(f"All figures saved to: {save_dir}/")
    print("=" * 60)

if __name__ == '__main__':
    main()
