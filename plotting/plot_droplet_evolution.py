#!/usr/bin/env python3
"""
Comprehensive Time-Series Plotting for Two-Phase Droplet Simulations

Generates multi-panel figures showing:
1. Per-phase entropy evolution (H_liquid, H_gas, H_interface)
2. Pressure evolution (P_liquid, P_gas, ΔP)
3. Surface tension and Laplace pressure
4. Spurious currents (global and interface-localized)
5. Interface properties (radius, width)
6. Conservation laws (mass, momentum, energy)

Usage:
    python plotting/plot_droplet_evolution.py output/bgk_diagnostics.csv output/elbm_diagnostics.csv
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
import sys

def load_diagnostics(filename):
    """Load diagnostics CSV file"""
    return pd.read_csv(filename)

def plot_entropy_evolution(ax, bgk_data, elbm_data):
    """Panel A: Per-phase entropy evolution"""
    ax.plot(bgk_data['timestep'], bgk_data['H_liquid_mean'],
            'r-', linewidth=2, label='BGK Liquid', alpha=0.7)
    ax.plot(bgk_data['timestep'], bgk_data['H_gas_mean'],
            'b-', linewidth=2, label='BGK Gas', alpha=0.7)
    ax.plot(bgk_data['timestep'], bgk_data['H_interface_mean'],
            'g-', linewidth=2, label='BGK Interface', alpha=0.7)

    ax.plot(elbm_data['timestep'], elbm_data['H_liquid_mean'],
            'r--', linewidth=2, label='ELBM Liquid', alpha=0.7)
    ax.plot(elbm_data['timestep'], elbm_data['H_gas_mean'],
            'b--', linewidth=2, label='ELBM Gas', alpha=0.7)
    ax.plot(elbm_data['timestep'], elbm_data['H_interface_mean'],
            'g--', linewidth=2, label='ELBM Interface', alpha=0.7)

    ax.set_xlabel('Time Step', fontsize=12)
    ax.set_ylabel('Mean H (entropy density)', fontsize=12)
    ax.set_title('A) Per-Phase Entropy Evolution', fontsize=14, weight='bold')
    ax.legend(fontsize=9, ncol=2)
    ax.grid(True, alpha=0.3)

    # Add annotation explaining different H-values
    ax.text(0.95, 0.05, 'Note: H_liquid ≠ H_gas is physically correct!',
            transform=ax.transAxes, fontsize=9,
            horizontalalignment='right', verticalalignment='bottom',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

def plot_pressure_evolution(ax, bgk_data, elbm_data):
    """Panel B: Pressure evolution"""
    ax.plot(bgk_data['timestep'], bgk_data['P_liquid'],
            'r-', linewidth=2, label='BGK P_in', alpha=0.7)
    ax.plot(bgk_data['timestep'], bgk_data['P_gas'],
            'b-', linewidth=2, label='BGK P_out', alpha=0.7)
    ax.plot(bgk_data['timestep'], bgk_data['laplace_pressure'],
            'k-', linewidth=2, label='BGK ΔP', alpha=0.7)

    ax.plot(elbm_data['timestep'], elbm_data['P_liquid'],
            'r--', linewidth=2, label='ELBM P_in', alpha=0.7)
    ax.plot(elbm_data['timestep'], elbm_data['P_gas'],
            'b--', linewidth=2, label='ELBM P_out', alpha=0.7)
    ax.plot(elbm_data['timestep'], elbm_data['laplace_pressure'],
            'k--', linewidth=2, label='ELBM ΔP', alpha=0.7)

    ax.set_xlabel('Time Step', fontsize=12)
    ax.set_ylabel('Pressure (lattice units)', fontsize=12)
    ax.set_title('B) Pressure Evolution', fontsize=14, weight='bold')
    ax.legend(fontsize=9, ncol=2)
    ax.grid(True, alpha=0.3)

def plot_surface_tension(ax, bgk_data, elbm_data):
    """Panel C: Surface tension evolution"""
    ax.plot(bgk_data['timestep'], bgk_data['surface_tension'],
            'purple', linewidth=2, label='BGK σ', alpha=0.7)
    ax.plot(elbm_data['timestep'], elbm_data['surface_tension'],
            'purple', linewidth=2, linestyle='--', label='ELBM σ', alpha=0.7)

    # Also plot Laplace pressure (should converge to σ/R)
    ax2 = ax.twinx()
    ax2.plot(bgk_data['timestep'], bgk_data['laplace_pressure'],
             'orange', linewidth=1.5, alpha=0.5, label='BGK ΔP')
    ax2.plot(elbm_data['timestep'], elbm_data['laplace_pressure'],
             'orange', linewidth=1.5, linestyle='--', alpha=0.5, label='ELBM ΔP')
    ax2.set_ylabel('Laplace Pressure ΔP', fontsize=11, color='orange')
    ax2.tick_params(axis='y', labelcolor='orange')

    ax.set_xlabel('Time Step', fontsize=12)
    ax.set_ylabel('Surface Tension σ', fontsize=12, color='purple')
    ax.set_title('C) Surface Tension & Laplace Pressure', fontsize=14, weight='bold')
    ax.tick_params(axis='y', labelcolor='purple')
    ax.legend(loc='upper left', fontsize=9)
    ax2.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)

def plot_spurious_currents(ax, bgk_data, elbm_data):
    """Panel D: Spurious currents (log scale)"""
    ax.semilogy(bgk_data['timestep'], bgk_data['spurious_vel'],
                'r-', linewidth=2, label='BGK Global', alpha=0.7)
    ax.semilogy(bgk_data['timestep'], bgk_data['spurious_vel_interface'],
                'r:', linewidth=2, label='BGK Interface', alpha=0.7)

    ax.semilogy(elbm_data['timestep'], elbm_data['spurious_vel'],
                'b--', linewidth=2, label='ELBM Global', alpha=0.7)
    ax.semilogy(elbm_data['timestep'], elbm_data['spurious_vel_interface'],
                'b-.', linewidth=2, label='ELBM Interface', alpha=0.7)

    ax.set_xlabel('Time Step', fontsize=12)
    ax.set_ylabel('Spurious Velocity (log scale)', fontsize=12)
    ax.set_title('D) Spurious Currents Evolution', fontsize=14, weight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which='both')

    # Add annotation showing ELBM improvement
    if len(elbm_data) > 0:
        bgk_final = bgk_data['spurious_vel'].iloc[-1]
        elbm_final = elbm_data['spurious_vel'].iloc[-1]
        improvement = bgk_final / elbm_final
        ax.text(0.05, 0.05, f'ELBM improvement: {improvement:.1f}×',
                transform=ax.transAxes, fontsize=10, weight='bold',
                horizontalalignment='left', verticalalignment='bottom',
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.7))

def plot_interface_properties(ax, bgk_data, elbm_data):
    """Panel E: Interface properties"""
    ax.plot(bgk_data['timestep'], bgk_data['interface_radius'],
            'purple', linewidth=2, label='BGK Radius', alpha=0.7)
    ax.plot(elbm_data['timestep'], elbm_data['interface_radius'],
            'purple', linewidth=2, linestyle='--', label='ELBM Radius', alpha=0.7)

    # Plot interface width on secondary axis
    ax2 = ax.twinx()
    ax2.plot(bgk_data['timestep'], bgk_data['interface_width'],
             'teal', linewidth=2, alpha=0.5, label='BGK Width')
    ax2.plot(elbm_data['timestep'], elbm_data['interface_width'],
             'teal', linewidth=2, linestyle='--', alpha=0.5, label='ELBM Width')
    ax2.set_ylabel('Interface Width', fontsize=11, color='teal')
    ax2.tick_params(axis='y', labelcolor='teal')

    ax.set_xlabel('Time Step', fontsize=12)
    ax.set_ylabel('Interface Radius', fontsize=12, color='purple')
    ax.set_title('E) Interface Radius & Width', fontsize=14, weight='bold')
    ax.tick_params(axis='y', labelcolor='purple')
    ax.legend(loc='upper left', fontsize=9)
    ax2.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)

def plot_conservation_laws(ax, bgk_data, elbm_data):
    """Panel F: Conservation laws"""
    # Compute conservation errors (relative to initial values)
    bgk_mass_err = np.abs(bgk_data['total_mass'] - bgk_data['total_mass'].iloc[0]) / bgk_data['total_mass'].iloc[0]
    elbm_mass_err = np.abs(elbm_data['total_mass'] - elbm_data['total_mass'].iloc[0]) / elbm_data['total_mass'].iloc[0]

    ax.semilogy(bgk_data['timestep'], bgk_mass_err,
                'r-', linewidth=2, label='BGK Mass Error', alpha=0.7)
    ax.semilogy(elbm_data['timestep'], elbm_mass_err,
                'r--', linewidth=2, label='ELBM Mass Error', alpha=0.7)

    ax.semilogy(bgk_data['timestep'], bgk_data['momentum_mag'],
                'b-', linewidth=2, label='BGK |Momentum|', alpha=0.7)
    ax.semilogy(elbm_data['timestep'], elbm_data['momentum_mag'],
                'b--', linewidth=2, label='ELBM |Momentum|', alpha=0.7)

    # Plot kinetic energy on secondary axis
    ax2 = ax.twinx()
    ax2.plot(bgk_data['timestep'], bgk_data['kinetic_energy'],
             'green', linewidth=2, alpha=0.5, label='BGK KE')
    ax2.plot(elbm_data['timestep'], elbm_data['kinetic_energy'],
             'green', linewidth=2, linestyle='--', alpha=0.5, label='ELBM KE')
    ax2.set_ylabel('Kinetic Energy', fontsize=11, color='green')
    ax2.tick_params(axis='y', labelcolor='green')
    ax2.set_yscale('log')

    ax.set_xlabel('Time Step', fontsize=12)
    ax.set_ylabel('Mass Error & Momentum (log scale)', fontsize=12)
    ax.set_title('F) Conservation Laws', fontsize=14, weight='bold')
    ax.legend(loc='upper left', fontsize=9)
    ax2.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3, which='both')

def create_comprehensive_figure(bgk_data, elbm_data, output_dir):
    """Create 6-panel comprehensive figure"""
    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.3, wspace=0.3)

    # Create 6 subplots
    ax1 = fig.add_subplot(gs[0, 0])  # Entropy
    ax2 = fig.add_subplot(gs[0, 1])  # Pressure
    ax3 = fig.add_subplot(gs[1, 0])  # Surface tension
    ax4 = fig.add_subplot(gs[1, 1])  # Spurious currents
    ax5 = fig.add_subplot(gs[2, 0])  # Interface properties
    ax6 = fig.add_subplot(gs[2, 1])  # Conservation laws

    # Plot each panel
    plot_entropy_evolution(ax1, bgk_data, elbm_data)
    plot_pressure_evolution(ax2, bgk_data, elbm_data)
    plot_surface_tension(ax3, bgk_data, elbm_data)
    plot_spurious_currents(ax4, bgk_data, elbm_data)
    plot_interface_properties(ax5, bgk_data, elbm_data)
    plot_conservation_laws(ax6, bgk_data, elbm_data)

    # Overall title
    fig.suptitle('Two-Phase Droplet Evolution: BGK vs ELBM Comparison',
                 fontsize=16, weight='bold', y=0.995)

    # Save figure
    output_path = Path(output_dir) / 'droplet_time_evolution.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved comprehensive figure to: {output_path}")

    # Also save as PDF for publication quality
    output_path_pdf = Path(output_dir) / 'droplet_time_evolution.pdf'
    plt.savefig(output_path_pdf, bbox_inches='tight')
    print(f"Saved PDF version to: {output_path_pdf}")

    plt.close()

def create_individual_plots(bgk_data, elbm_data, output_dir):
    """Create individual plots for each panel (easier to view/share)"""
    output_dir = Path(output_dir)

    # Panel A: Entropy
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_entropy_evolution(ax, bgk_data, elbm_data)
    plt.savefig(output_dir / 'entropy_evolution.png', dpi=200, bbox_inches='tight')
    plt.close()

    # Panel B: Pressure
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_pressure_evolution(ax, bgk_data, elbm_data)
    plt.savefig(output_dir / 'pressure_evolution.png', dpi=200, bbox_inches='tight')
    plt.close()

    # Panel C: Surface tension
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_surface_tension(ax, bgk_data, elbm_data)
    plt.savefig(output_dir / 'surface_tension_evolution.png', dpi=200, bbox_inches='tight')
    plt.close()

    # Panel D: Spurious currents
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_spurious_currents(ax, bgk_data, elbm_data)
    plt.savefig(output_dir / 'spurious_currents_evolution.png', dpi=200, bbox_inches='tight')
    plt.close()

    # Panel E: Interface
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_interface_properties(ax, bgk_data, elbm_data)
    plt.savefig(output_dir / 'interface_evolution.png', dpi=200, bbox_inches='tight')
    plt.close()

    # Panel F: Conservation
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_conservation_laws(ax, bgk_data, elbm_data)
    plt.savefig(output_dir / 'conservation_evolution.png', dpi=200, bbox_inches='tight')
    plt.close()

    print(f"Saved 6 individual plots to: {output_dir}")

def print_final_statistics(bgk_data, elbm_data):
    """Print summary statistics for comparison"""
    print("\n" + "="*60)
    print("FINAL STATE COMPARISON (BGK vs ELBM)")
    print("="*60)

    print("\nEntropy (Mean H per node):")
    print(f"  BGK  - Liquid: {bgk_data['H_liquid_mean'].iloc[-1]:.4f}, Gas: {bgk_data['H_gas_mean'].iloc[-1]:.4f}")
    print(f"  ELBM - Liquid: {elbm_data['H_liquid_mean'].iloc[-1]:.4f}, Gas: {elbm_data['H_gas_mean'].iloc[-1]:.4f}")

    print("\nPressure:")
    print(f"  BGK  - ΔP: {bgk_data['laplace_pressure'].iloc[-1]:.6e}")
    print(f"  ELBM - ΔP: {elbm_data['laplace_pressure'].iloc[-1]:.6e}")

    print("\nSurface Tension:")
    print(f"  BGK  - σ: {bgk_data['surface_tension'].iloc[-1]:.6e}")
    print(f"  ELBM - σ: {elbm_data['surface_tension'].iloc[-1]:.6e}")

    print("\nSpurious Currents:")
    print(f"  BGK  - Global: {bgk_data['spurious_vel'].iloc[-1]:.6e}, Interface: {bgk_data['spurious_vel_interface'].iloc[-1]:.6e}")
    print(f"  ELBM - Global: {elbm_data['spurious_vel'].iloc[-1]:.6e}, Interface: {elbm_data['spurious_vel_interface'].iloc[-1]:.6e}")
    improvement = bgk_data['spurious_vel'].iloc[-1] / elbm_data['spurious_vel'].iloc[-1]
    print(f"  ELBM Improvement Factor: {improvement:.1f}×")

    print("\nInterface:")
    print(f"  BGK  - Radius: {bgk_data['interface_radius'].iloc[-1]:.2f}, Width: {bgk_data['interface_width'].iloc[-1]:.2f}")
    print(f"  ELBM - Radius: {elbm_data['interface_radius'].iloc[-1]:.2f}, Width: {elbm_data['interface_width'].iloc[-1]:.2f}")

    print("\nConservation (Final values):")
    bgk_mass_err = abs(bgk_data['total_mass'].iloc[-1] - bgk_data['total_mass'].iloc[0]) / bgk_data['total_mass'].iloc[0]
    elbm_mass_err = abs(elbm_data['total_mass'].iloc[-1] - elbm_data['total_mass'].iloc[0]) / elbm_data['total_mass'].iloc[0]
    print(f"  BGK  - Mass error: {bgk_mass_err:.6e}, Momentum: {bgk_data['momentum_mag'].iloc[-1]:.6e}")
    print(f"  ELBM - Mass error: {elbm_mass_err:.6e}, Momentum: {elbm_data['momentum_mag'].iloc[-1]:.6e}")

    print("="*60 + "\n")

def main():
    if len(sys.argv) < 3:
        print("Usage: python plot_droplet_evolution.py <bgk_diagnostics.csv> <elbm_diagnostics.csv> [output_dir]")
        print("Example: python plotting/plot_droplet_evolution.py output/bgk_diagnostics.csv output/elbm_diagnostics.csv figures/")
        sys.exit(1)

    bgk_file = sys.argv[1]
    elbm_file = sys.argv[2]
    output_dir = sys.argv[3] if len(sys.argv) > 3 else 'figures'

    # Create output directory if needed
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Load data
    print(f"Loading BGK data from: {bgk_file}")
    bgk_data = load_diagnostics(bgk_file)
    print(f"  → Loaded {len(bgk_data)} time steps")

    print(f"Loading ELBM data from: {elbm_file}")
    elbm_data = load_diagnostics(elbm_file)
    print(f"  → Loaded {len(elbm_data)} time steps")

    # Print statistics
    print_final_statistics(bgk_data, elbm_data)

    # Create plots
    print("\nGenerating plots...")
    create_comprehensive_figure(bgk_data, elbm_data, output_dir)
    create_individual_plots(bgk_data, elbm_data, output_dir)

    print("\n✓ Plotting complete!")
    print(f"  → Comprehensive figure: {output_dir}/droplet_time_evolution.png")
    print(f"  → Individual plots: {output_dir}/*.png")

if __name__ == "__main__":
    main()
