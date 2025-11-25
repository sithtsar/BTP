#!/usr/bin/env python3
"""
Plot H-theorem validation for cylinder flow benchmarks
Compares BGK vs ELBM at Re = 10, 50, 100
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

plt.style.use('seaborn-v0_8-darkgrid')
plt.rcParams['font.size'] = 11

def load_entropy_history(filepath):
    """Load entropy history from EntropyMonitor output"""
    if not Path(filepath).exists():
        return None
    data = np.loadtxt(filepath, comments='#')
    return {
        'timestep': data[:, 0].astype(int),
        'H_total': data[:, 1],
        'H_min': data[:, 2],
        'H_max': data[:, 3],
        'H_mean': data[:, 4],
        'spurious_vel': data[:, 5]
    }

def plot_combined_reynolds(re_numbers):
    """Combined plot showing all Reynolds numbers"""
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('Cylinder Flow: H-Theorem Validation Across Reynolds Numbers',
                 fontsize=16, fontweight='bold')

    colors_bgk = ['#1f77b4', '#ff7f0e', '#2ca02c']
    colors_elbm = ['#d62728', '#9467bd', '#8c564b']

    # Plot 1: H(t) for all Re
    ax = axes[0, 0]
    for i, re in enumerate(re_numbers):
        bgk_file = Path(f'output/channel_cylinder/re{re}_BGK_entropy.dat')
        elbm_file = Path(f'output/channel_cylinder/re{re}_ELBM_entropy.dat')

        if bgk_file.exists():
            data = load_entropy_history(bgk_file)
            ax.plot(data['timestep'], data['H_total'],
                   color=colors_bgk[i], linestyle='-', linewidth=2,
                   label=f'BGK Re={re}', alpha=0.7)

        if elbm_file.exists():
            data = load_entropy_history(elbm_file)
            ax.plot(data['timestep'], data['H_total'],
                   color=colors_elbm[i], linestyle='--', linewidth=2,
                   label=f'ELBM Re={re}', alpha=0.7)

    ax.set_xlabel('Timestep')
    ax.set_ylabel('H (Total Entropy)')
    ax.set_title('Entropy Evolution: All Re')
    ax.legend(ncol=2, fontsize=9)
    ax.grid(True, alpha=0.3)

    # Plot 2: Spurious currents
    ax = axes[0, 1]
    for i, re in enumerate(re_numbers):
        bgk_file = Path(f'output/channel_cylinder/re{re}_BGK_entropy.dat')
        elbm_file = Path(f'output/channel_cylinder/re{re}_ELBM_entropy.dat')

        if bgk_file.exists():
            data = load_entropy_history(bgk_file)
            ax.semilogy(data['timestep'], data['spurious_vel'],
                       color=colors_bgk[i], linestyle='-', linewidth=2,
                       label=f'BGK Re={re}', alpha=0.7)

        if elbm_file.exists():
            data = load_entropy_history(elbm_file)
            ax.semilogy(data['timestep'], data['spurious_vel'],
                       color=colors_elbm[i], linestyle='--', linewidth=2,
                       label=f'ELBM Re={re}', alpha=0.7)

    ax.set_xlabel('Timestep')
    ax.set_ylabel('Max Velocity (Spurious)')
    ax.set_title('Spurious Currents')
    ax.legend(ncol=2, fontsize=9)
    ax.grid(True, alpha=0.3)

    # Plot 3: dH/dt for lowest Re
    ax = axes[1, 0]
    re_low = re_numbers[0]
    bgk_file = Path(f'output/channel_cylinder/re{re_low}_BGK_entropy.dat')
    elbm_file = Path(f'output/channel_cylinder/re{re_low}_ELBM_entropy.dat')

    if bgk_file.exists():
        data = load_entropy_history(bgk_file)
        dH_bgk = np.diff(data['H_total']) / np.diff(data['timestep'])
        ax.plot(data['timestep'][1:], dH_bgk, 'b-', linewidth=2, label='BGK', alpha=0.7)

    if elbm_file.exists():
        data = load_entropy_history(elbm_file)
        dH_elbm = np.diff(data['H_total']) / np.diff(data['timestep'])
        ax.plot(data['timestep'][1:], dH_elbm, 'r-', linewidth=2, label='ELBM', alpha=0.7)

    ax.axhline(0, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel('Timestep')
    ax.set_ylabel('dH/dt')
    ax.set_title(f'Entropy Rate: Re={re_low}')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 4: dH/dt for highest Re
    ax = axes[1, 1]
    re_high = re_numbers[-1]
    bgk_file = Path(f'output/channel_cylinder/re{re_high}_BGK_entropy.dat')
    elbm_file = Path(f'output/channel_cylinder/re{re_high}_ELBM_entropy.dat')

    if bgk_file.exists():
        data = load_entropy_history(bgk_file)
        dH_bgk = np.diff(data['H_total']) / np.diff(data['timestep'])
        ax.plot(data['timestep'][1:], dH_bgk, 'b-', linewidth=2, label='BGK', alpha=0.7)

    if elbm_file.exists():
        data = load_entropy_history(elbm_file)
        dH_elbm = np.diff(data['H_total']) / np.diff(data['timestep'])
        ax.plot(data['timestep'][1:], dH_elbm, 'r-', linewidth=2, label='ELBM', alpha=0.7)

    ax.axhline(0, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel('Timestep')
    ax.set_ylabel('dH/dt')
    ax.set_title(f'Entropy Rate: Re={re_high}')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    fig_dir = Path('figures/channel_cylinder')
    fig_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(fig_dir / 'entropy_validation_combined.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {fig_dir / 'entropy_validation_combined.png'}")

def plot_individual_reynolds(re):
    """Detailed plot for single Reynolds number"""
    bgk_file = Path(f'output/channel_cylinder/re{re}_BGK_entropy.dat')
    elbm_file = Path(f'output/channel_cylinder/re{re}_ELBM_entropy.dat')

    if not (bgk_file.exists() and elbm_file.exists()):
        print(f"Skipping Re={re}: files not found")
        return

    bgk_data = load_entropy_history(bgk_file)
    elbm_data = load_entropy_history(elbm_file)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'Cylinder Flow Entropy Analysis: Re={re}', fontsize=16, fontweight='bold')

    # H(t)
    ax = axes[0, 0]
    ax.plot(bgk_data['timestep'], bgk_data['H_total'], 'b-', linewidth=2, label='BGK')
    ax.plot(elbm_data['timestep'], elbm_data['H_total'], 'r-', linewidth=2, label='ELBM')
    ax.set_xlabel('Timestep')
    ax.set_ylabel('H (Total Entropy)')
    ax.set_title('Entropy Evolution')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # dH/dt
    ax = axes[0, 1]
    dH_bgk = np.diff(bgk_data['H_total']) / np.diff(bgk_data['timestep'])
    dH_elbm = np.diff(elbm_data['H_total']) / np.diff(elbm_data['timestep'])
    ax.plot(bgk_data['timestep'][1:], dH_bgk, 'b-', linewidth=2, label='BGK')
    ax.plot(elbm_data['timestep'][1:], dH_elbm, 'r-', linewidth=2, label='ELBM')
    ax.axhline(0, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel('Timestep')
    ax.set_ylabel('dH/dt')
    ax.set_title('Entropy Rate (should be ≤ 0)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Spurious currents
    ax = axes[1, 0]
    ax.semilogy(bgk_data['timestep'], bgk_data['spurious_vel'], 'b-', linewidth=2, label='BGK')
    ax.semilogy(elbm_data['timestep'], elbm_data['spurious_vel'], 'r-', linewidth=2, label='ELBM')
    ax.set_xlabel('Timestep')
    ax.set_ylabel('Max Velocity')
    ax.set_title('Spurious Currents')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # H range
    ax = axes[1, 1]
    ax.plot(bgk_data['timestep'], bgk_data['H_max'], 'b-', linewidth=1.5, label='BGK max', alpha=0.7)
    ax.plot(bgk_data['timestep'], bgk_data['H_mean'], 'b--', linewidth=2, label='BGK mean')
    ax.plot(bgk_data['timestep'], bgk_data['H_min'], 'b-', linewidth=1.5, label='BGK min', alpha=0.7)
    ax.plot(elbm_data['timestep'], elbm_data['H_max'], 'r-', linewidth=1.5, label='ELBM max', alpha=0.7)
    ax.plot(elbm_data['timestep'], elbm_data['H_mean'], 'r--', linewidth=2, label='ELBM mean')
    ax.plot(elbm_data['timestep'], elbm_data['H_min'], 'r-', linewidth=1.5, label='ELBM min', alpha=0.7)
    ax.set_xlabel('Timestep')
    ax.set_ylabel('H per node')
    ax.set_title('Entropy Distribution (min/mean/max)')
    ax.legend(ncol=2, fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    fig_dir = Path('figures/channel_cylinder')
    fig_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(fig_dir / f're_{re}_entropy_detailed.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {fig_dir / f're_{re}_entropy_detailed.png'}")

def main():
    re_numbers = [10, 50, 100]

    # Combined plot
    plot_combined_reynolds(re_numbers)

    # Individual plots
    for re in re_numbers:
        plot_individual_reynolds(re)

    print("\n✓ All cylinder entropy validation plots generated successfully!")

if __name__ == '__main__':
    main()
