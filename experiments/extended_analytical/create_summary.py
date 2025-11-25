#!/usr/bin/env python3
"""
Create summary figure combining all extended analytical cases
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def main():
    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.3)

    # Load data
    womersley_data = np.loadtxt('output/womersley_analytical.dat')
    hagen_data = np.loadtxt('output/hagen_poiseuille_analytical.dat')
    stokes_data = np.loadtxt('output/stokes_shear_analytical.dat')
    kolm_data = np.loadtxt('output/kolmogorov_analytical.dat')

    # 1. Womersley Flow - show 4 time snapshots
    ax1 = fig.add_subplot(gs[0, 0])
    y_wom = womersley_data[:, 0]
    for i in [1, 3, 5, 7]:  # Select 4 time points
        u_wom = womersley_data[:, i]
        ax1.plot(u_wom, y_wom, linewidth=2, label=f't={i-1}', alpha=0.7)
    ax1.set_xlabel('Velocity')
    ax1.set_ylabel('y position')
    ax1.set_title('(a) Womersley Flow\n(Oscillatory, α=31.7)', fontsize=11, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8)
    ax1.axvline(x=0, color='k', linestyle='--', alpha=0.2)

    # 2. Hagen-Poiseuille Flow
    ax2 = fig.add_subplot(gs[0, 1])
    r_hagen = hagen_data[:, 0]
    u_hagen = hagen_data[:, 1]
    ax2.plot(r_hagen, u_hagen, 'b-', linewidth=2.5)
    ax2.axvline(x=30, color='r', linestyle='--', linewidth=2, alpha=0.5, label='Wall (r=30)')
    ax2.set_xlabel('Radial distance r')
    ax2.set_ylabel('Velocity u(r)')
    ax2.set_title('(b) Hagen-Poiseuille Flow\n(Circular Pipe, Parabolic)', fontsize=11, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_ylim(bottom=-0.2)

    # 3. Stokes Shear Flow
    ax3 = fig.add_subplot(gs[0, 2])
    y_stokes = stokes_data[:, 0]
    u_stokes = stokes_data[:, 1]
    du_dy = stokes_data[:, 2]
    ax3.plot(y_stokes, u_stokes, 'b-', linewidth=2.5, label='Velocity')
    ax3_twin = ax3.twinx()
    ax3_twin.plot(y_stokes[5:-5], du_dy[5:-5], 'g--', linewidth=2, alpha=0.7, label='du/dy')
    ax3_twin.axhline(y=0.01, color='r', linestyle=':', linewidth=2, label='γ=0.01')
    ax3.set_xlabel('y position')
    ax3.set_ylabel('Velocity u(y)', color='b')
    ax3_twin.set_ylabel('Shear rate du/dy', color='g')
    ax3.set_title('(c) Stokes Shear Flow\n(Creeping Flow, Constant Shear)', fontsize=11, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.tick_params(axis='y', labelcolor='b')
    ax3_twin.tick_params(axis='y', labelcolor='g')
    lines1, labels1 = ax3.get_legend_handles_labels()
    lines2, labels2 = ax3_twin.get_legend_handles_labels()
    ax3.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=8)

    # 4. Kolmogorov Flow
    ax4 = fig.add_subplot(gs[1, 0])
    y_kolm = kolm_data[:, 0]
    u_kolm = kolm_data[:, 1]
    ax4.plot(y_kolm, u_kolm, 'b-', linewidth=2.5)
    ax4.set_xlabel('y position')
    ax4.set_ylabel('Velocity u(y)')
    ax4.set_title('(d) Kolmogorov Flow\n(Sinusoidal Forcing, k=0.049)', fontsize=11, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.axhline(y=0, color='k', linestyle='--', alpha=0.2)

    # 5. Hagen-Poiseuille 2D cross-section
    ax5 = fig.add_subplot(gs[1, 1])
    R = 30.0
    theta = np.linspace(0, 2*np.pi, 100)
    x_circle = R * np.cos(theta)
    y_circle = R * np.sin(theta)

    x_grid = np.linspace(-35, 35, 100)
    y_grid = np.linspace(-35, 35, 100)
    X, Y = np.meshgrid(x_grid, y_grid)
    R_grid = np.sqrt(X**2 + Y**2)

    # Recreate velocity field
    U = -0.001 * (R**2 - R_grid**2) / (4.0 * 0.1)
    U[R_grid > R] = 0

    contour = ax5.contourf(X, Y, U, levels=20, cmap='viridis')
    ax5.plot(x_circle, y_circle, 'r-', linewidth=2.5, label='Pipe wall')
    ax5.set_xlabel('x position')
    ax5.set_ylabel('y position')
    ax5.set_title('(e) Hagen-Poiseuille\n2D Cross-Section', fontsize=11, fontweight='bold')
    ax5.set_aspect('equal')
    cbar = plt.colorbar(contour, ax=ax5, label='Velocity')
    cbar.ax.tick_params(labelsize=8)
    ax5.legend(fontsize=9)

    # 6. Kolmogorov 2D field
    ax6 = fig.add_subplot(gs[1, 2])
    nx = 128
    X_kolm, Y_kolm = np.meshgrid(np.arange(nx), y_kolm)
    U_kolm_2d = np.tile(u_kolm[:, np.newaxis], (1, nx))

    contour2 = ax6.contourf(X_kolm, Y_kolm, U_kolm_2d, levels=20, cmap='RdBu_r')
    ax6.set_xlabel('x position')
    ax6.set_ylabel('y position')
    ax6.set_title('(f) Kolmogorov Flow\n2D Velocity Field', fontsize=11, fontweight='bold')
    ax6.set_aspect('equal')
    cbar2 = plt.colorbar(contour2, ax=ax6, label='Velocity')
    cbar2.ax.tick_params(labelsize=8)

    fig.suptitle('Extended Analytical Validation Cases - Complete Suite',
                 fontsize=18, fontweight='bold', y=0.98)

    # Add text box with summary
    textstr = 'Analytical Solutions Implemented:\n'
    textstr += '• Womersley: Pulsatile flow (cardiovascular)\n'
    textstr += '• Hagen-Poiseuille: Circular pipe (3D geometry)\n'
    textstr += '• Stokes: Creeping flow (low Re)\n'
    textstr += '• Kolmogorov: Turbulence transition'

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    fig.text(0.99, 0.02, textstr, fontsize=9, verticalalignment='bottom',
             horizontalalignment='right', bbox=props, family='monospace')

    plt.tight_layout(rect=[0, 0.03, 1, 0.96])

    # Save
    output_dir = Path('figures/extended_analytical')
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / 'extended_analytical_complete_summary.png', dpi=200, bbox_inches='tight')
    print(f"\n{'='*70}")
    print(f"Complete summary figure saved:")
    print(f"  {output_dir / 'extended_analytical_complete_summary.png'}")
    print(f"{'='*70}\n")

    plt.close()

if __name__ == '__main__':
    main()
