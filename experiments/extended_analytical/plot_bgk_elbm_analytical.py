#!/usr/bin/env python3
"""
Plot BGK vs ELBM vs Analytical Solutions
3-way comparison for extended analytical test cases
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def plot_hagen_poiseuille():
    """Plot Hagen-Poiseuille flow comparison"""
    try:
        bgk_data = np.loadtxt('output/hagen_poiseuille_bgk.dat')
        elbm_data = np.loadtxt('output/hagen_poiseuille_elbm.dat')

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Radial profile
        ax = axes[0]
        ax.plot(bgk_data[:, 0], bgk_data[:, 2], 'k--', linewidth=2, label='Analytical', alpha=0.7)
        ax.plot(bgk_data[:, 0], bgk_data[:, 1], 'b-', linewidth=1.5, label='BGK', alpha=0.8)
        ax.plot(elbm_data[:, 0], elbm_data[:, 1], 'r-', linewidth=1.5, label='ELBM', alpha=0.8)
        ax.set_xlabel('Radial distance r', fontsize=12)
        ax.set_ylabel('Velocity u(r)', fontsize=12)
        ax.set_title('Hagen-Poiseuille Flow: Radial Profile', fontsize=13, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        # Error comparison
        ax = axes[1]
        ax.semilogy(bgk_data[:, 0], np.abs(bgk_data[:, 3]), 'b-', linewidth=1.5, label='BGK error')
        ax.semilogy(elbm_data[:, 0], np.abs(elbm_data[:, 3]), 'r-', linewidth=1.5, label='ELBM error')
        ax.set_xlabel('Radial distance r', fontsize=12)
        ax.set_ylabel('Absolute Error', fontsize=12)
        ax.set_title('Error Comparison', fontsize=13, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        output_dir = Path('figures/extended_analytical')
        output_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_dir / 'hagen_poiseuille_comparison.png', dpi=200, bbox_inches='tight')
        print(f"Saved: {output_dir / 'hagen_poiseuille_comparison.png'}")
        plt.close()

    except Exception as e:
        print(f"Error plotting Hagen-Poiseuille: {e}")

def plot_stokes_shear():
    """Plot Stokes shear flow comparison"""
    try:
        bgk_data = np.loadtxt('output/stokes_shear_bgk.dat')
        elbm_data = np.loadtxt('output/stokes_shear_elbm.dat')

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Velocity profile
        ax = axes[0]
        ax.plot(bgk_data[:, 2], bgk_data[:, 0], 'k--', linewidth=2, label='Analytical', alpha=0.7)
        ax.plot(bgk_data[:, 1], bgk_data[:, 0], 'b-', linewidth=1.5, label='BGK', alpha=0.8)
        ax.plot(elbm_data[:, 1], elbm_data[:, 0], 'r-', linewidth=1.5, label='ELBM', alpha=0.8)
        ax.set_xlabel('Velocity u(y)', fontsize=12)
        ax.set_ylabel('y position', fontsize=12)
        ax.set_title('Stokes Shear Flow: Velocity Profile', fontsize=13, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        # Error comparison
        ax = axes[1]
        ax.semilogy(bgk_data[:, 0], np.abs(bgk_data[:, 3]), 'b-', linewidth=1.5, label='BGK error')
        ax.semilogy(elbm_data[:, 0], np.abs(elbm_data[:, 3]), 'r-', linewidth=1.5, label='ELBM error')
        ax.set_xlabel('y position', fontsize=12)
        ax.set_ylabel('Absolute Error', fontsize=12)
        ax.set_title('Error Comparison', fontsize=13, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        output_dir = Path('figures/extended_analytical')
        output_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_dir / 'stokes_shear_comparison.png', dpi=200, bbox_inches='tight')
        print(f"Saved: {output_dir / 'stokes_shear_comparison.png'}")
        plt.close()

    except Exception as e:
        print(f"Error plotting Stokes shear: {e}")

def create_combined_summary():
    """Create combined summary figure"""
    try:
        bgk_hp = np.loadtxt('output/hagen_poiseuille_bgk.dat')
        elbm_hp = np.loadtxt('output/hagen_poiseuille_elbm.dat')
        bgk_ss = np.loadtxt('output/stokes_shear_bgk.dat')
        elbm_ss = np.loadtxt('output/stokes_shear_elbm.dat')

        fig = plt.figure(figsize=(16, 10))
        gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

        # Hagen-Poiseuille profile
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(bgk_hp[:, 0], bgk_hp[:, 2], 'k--', linewidth=2.5, label='Analytical', alpha=0.7)
        ax1.plot(bgk_hp[:, 0], bgk_hp[:, 1], 'b-', linewidth=1.5, label='BGK', alpha=0.8)
        ax1.plot(elbm_hp[:, 0], elbm_hp[:, 1], 'r-', linewidth=1.5, label='ELBM', alpha=0.8)
        ax1.set_xlabel('Radial distance r', fontsize=11)
        ax1.set_ylabel('Velocity u(r)', fontsize=11)
        ax1.set_title('(a) Hagen-Poiseuille: Radial Profile', fontsize=12, fontweight='bold')
        ax1.legend(fontsize=9, loc='upper right')
        ax1.grid(True, alpha=0.3)

        # Hagen-Poiseuille error
        ax2 = fig.add_subplot(gs[0, 1])
        ax2.semilogy(bgk_hp[:, 0], np.abs(bgk_hp[:, 3]), 'b-', linewidth=1.5, label='BGK', alpha=0.8)
        ax2.semilogy(elbm_hp[:, 0], np.abs(elbm_hp[:, 3]), 'r-', linewidth=1.5, label='ELBM', alpha=0.8)
        ax2.set_xlabel('Radial distance r', fontsize=11)
        ax2.set_ylabel('Absolute Error', fontsize=11)
        ax2.set_title('(b) Hagen-Poiseuille: Error', fontsize=12, fontweight='bold')
        ax2.legend(fontsize=9)
        ax2.grid(True, alpha=0.3)

        # Stokes profile
        ax3 = fig.add_subplot(gs[1, 0])
        ax3.plot(bgk_ss[:, 2], bgk_ss[:, 0], 'k--', linewidth=2.5, label='Analytical', alpha=0.7)
        ax3.plot(bgk_ss[:, 1], bgk_ss[:, 0], 'b-', linewidth=1.5, label='BGK', alpha=0.8)
        ax3.plot(elbm_ss[:, 1], elbm_ss[:, 0], 'r-', linewidth=1.5, label='ELBM', alpha=0.8)
        ax3.set_xlabel('Velocity u(y)', fontsize=11)
        ax3.set_ylabel('y position', fontsize=11)
        ax3.set_title('(c) Stokes Shear: Velocity Profile', fontsize=12, fontweight='bold')
        ax3.legend(fontsize=9, loc='upper left')
        ax3.grid(True, alpha=0.3)

        # Stokes error
        ax4 = fig.add_subplot(gs[1, 1])
        ax4.semilogy(bgk_ss[:, 0], np.abs(bgk_ss[:, 3]), 'b-', linewidth=1.5, label='BGK', alpha=0.8)
        ax4.semilogy(elbm_ss[:, 0], np.abs(elbm_ss[:, 3]), 'r-', linewidth=1.5, label='ELBM', alpha=0.8)
        ax4.set_xlabel('y position', fontsize=11)
        ax4.set_ylabel('Absolute Error', fontsize=11)
        ax4.set_title('(d) Stokes Shear: Error', fontsize=12, fontweight='bold')
        ax4.legend(fontsize=9)
        ax4.grid(True, alpha=0.3)

        fig.suptitle('Extended Analytical Validation: BGK vs ELBM vs Analytical',
                     fontsize=16, fontweight='bold', y=0.995)

        # Add summary text
        textstr = 'Comparison Summary:\\n'
        textstr += f'• Hagen-Poiseuille (Circular Pipe)\\n'
        textstr += f'• Stokes Shear (Linear Shear Flow)\\n'
        textstr += f'• Both BGK and ELBM compared to analytical'

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        fig.text(0.99, 0.02, textstr, fontsize=9, verticalalignment='bottom',
                 horizontalalignment='right', bbox=props, family='monospace')

        plt.tight_layout(rect=[0, 0.03, 1, 0.98])

        output_dir = Path('figures/extended_analytical')
        output_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_dir / 'bgk_elbm_analytical_summary.png', dpi=200, bbox_inches='tight')
        print(f"\\nSaved: {output_dir / 'bgk_elbm_analytical_summary.png'}")
        plt.close()

    except Exception as e:
        print(f"Error creating combined summary: {e}")

def main():
    print("\\n" + "="*70)
    print("Plotting BGK vs ELBM vs Analytical Comparisons")
    print("="*70 + "\\n")

    plot_hagen_poiseuille()
    plot_stokes_shear()
    create_combined_summary()

    print("\\n" + "="*70)
    print("All comparison plots generated successfully")
    print("="*70 + "\\n")

if __name__ == '__main__':
    main()
