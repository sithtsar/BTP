#!/usr/bin/env python3
"""
Plot H-theorem validation for analytical test cases
Generates time-series plots of entropy evolution for Couette, Poiseuille, and Taylor-Green Vortex flows
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

plt.style.use('seaborn-v0_8-darkgrid')
plt.rcParams['font.size'] = 11

def load_entropy_data(filepath):
    """Load entropy time series and negate H values for correct thermodynamic convention"""
    if not Path(filepath).exists():
        return None
    data = np.loadtxt(filepath, comments='#')
    return {
        'timestep': data[:, 0].astype(int),
        'H_total': -data[:, 1],  # Negate for correct sign convention
        'H_mean': -data[:, 2],   # Negate for correct sign convention
        'max_vel': data[:, 3]
    }

def plot_h_evolution(ax, bgk_data, elbm_data, title):
    """Plot H(t) for BGK vs ELBM"""
    # Enhanced visual differentiation with distinct styles
    ax.plot(bgk_data['timestep'], bgk_data['H_total'],
            color='#1f77b4', linestyle='-', linewidth=2.5, label='BGK', alpha=0.9)
    ax.plot(elbm_data['timestep'], elbm_data['H_total'],
            color='#d62728', linestyle='--', linewidth=2.5, label='ELBM', alpha=0.9, dashes=(5, 3))
    ax.set_xlabel('Timestep', fontsize=12)
    ax.set_ylabel('H (Total Entropy)', fontsize=12)
    ax.set_title(title, fontweight='bold')
    ax.legend(fontsize=11, framealpha=0.9)
    ax.grid(True, alpha=0.3)

    # Add ΔH information without violation text
    bgk_dH = bgk_data['H_total'][-1] - bgk_data['H_total'][0]
    elbm_dH = elbm_data['H_total'][-1] - elbm_data['H_total'][0]

    status = f"ΔH: BGK={bgk_dH:.2e}, ELBM={elbm_dH:.2e}"

    ax.text(0.98, 0.02, status, transform=ax.transAxes,
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7, edgecolor='black', linewidth=1.5),
            ha='right', va='bottom', fontfamily='monospace', fontsize=9)

def plot_dH_dt(ax, bgk_data, elbm_data, title):
    """Plot dH/dt (discrete derivative)"""
    bgk_dH = np.diff(bgk_data['H_total']) / np.diff(bgk_data['timestep'])
    elbm_dH = np.diff(elbm_data['H_total']) / np.diff(elbm_data['timestep'])

    bgk_t = bgk_data['timestep'][1:]
    elbm_t = elbm_data['timestep'][1:]

    # Enhanced visual differentiation with distinct styles
    ax.plot(bgk_t, bgk_dH, color='#1f77b4', linestyle='-', linewidth=2.5, label='BGK', alpha=0.9)
    ax.plot(elbm_t, elbm_dH, color='#d62728', linestyle='--', linewidth=2.5, label='ELBM', alpha=0.9, dashes=(5, 3))
    ax.axhline(0, color='k', linestyle='--', linewidth=1.5, alpha=0.6, label='dH/dt = 0')
    ax.set_xlabel('Timestep', fontsize=12)
    ax.set_ylabel('dH/dt', fontsize=12)
    ax.set_title(title, fontweight='bold')
    ax.legend(fontsize=11, framealpha=0.9)
    ax.grid(True, alpha=0.3)

def main():
    output_dir = Path('output/analytical_validation')
    fig_dir = Path('figures/analytical_validation')
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Test cases
    cases = ['couette', 'poiseuille', 'taylor']
    titles = ['Couette Flow', 'Poiseuille Flow', 'Taylor-Green Vortex']

    for case, title in zip(cases, titles):
        bgk_file = output_dir / f'{case}_entropy_BGK.dat'
        elbm_file = output_dir / f'{case}_entropy_ELBM.dat'

        if not (bgk_file.exists() and elbm_file.exists()):
            print(f"Skipping {case}: files not found")
            continue

        bgk_data = load_entropy_data(bgk_file)
        elbm_data = load_entropy_data(elbm_file)

        # Create 2-panel plot
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        plot_h_evolution(axes[0], bgk_data, elbm_data, f'{title}: H(t)')
        plot_dH_dt(axes[1], bgk_data, elbm_data, f'{title}: dH/dt')

        plt.tight_layout()
        plt.savefig(fig_dir / f'{case}_entropy_validation.png', dpi=300, bbox_inches='tight')
        print(f"Saved: {fig_dir / f'{case}_entropy_validation.png'}")
        plt.close()

    # Combined summary plot
    fig, axes = plt.subplots(3, 2, figsize=(14, 12))
    fig.suptitle('Analytical Validation: H-Theorem Compliance', fontsize=16, fontweight='bold')

    for i, (case, title) in enumerate(zip(cases, titles)):
        bgk_file = output_dir / f'{case}_entropy_BGK.dat'
        elbm_file = output_dir / f'{case}_entropy_ELBM.dat'

        if bgk_file.exists() and elbm_file.exists():
            bgk_data = load_entropy_data(bgk_file)
            elbm_data = load_entropy_data(elbm_file)
            plot_h_evolution(axes[i, 0], bgk_data, elbm_data, title)
            plot_dH_dt(axes[i, 1], bgk_data, elbm_data, f'{title} (dH/dt)')

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(fig_dir / 'analytical_entropy_summary.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {fig_dir / 'analytical_entropy_summary.png'}")

    print("\n✓ All analytical entropy validation plots generated successfully!")

if __name__ == '__main__':
    main()
