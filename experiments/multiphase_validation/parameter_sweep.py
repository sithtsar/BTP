"""
Parameter Sweep for Spurious Current Optimization

Systematically tests different parameter combinations to minimize
spurious currents in static droplet simulations.

Parameters varied:
1. Interaction strength (g)
2. Relaxation parameter (omega)
3. Density ratio (rho_liquid/rho_gas)
4. Time steps (equilibration time)

Objective: Find parameter combination with spurious currents < 1e-3
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import json
from datetime import datetime

sys.path.append(str(Path(__file__).parent))

from static_droplet_mrt import StaticDroplet


def run_parameter_sweep(
    g_values: list = None,
    omega_values: list = None,
    density_ratios: list = None,
    time_steps: int = 20000,
    N: int = 64,
    droplet_radius: float = 15.0
):
    """
    Run parameter sweep to optimize spurious currents.

    Args:
        g_values: List of interaction strengths to test
        omega_values: List of relaxation parameters to test
        density_ratios: List of density ratios to test
        time_steps: Number of equilibration steps
        N: Domain size
        droplet_radius: Droplet radius

    Returns:
        Dictionary with results
    """
    # Default parameter ranges
    if g_values is None:
        g_values = [-5.0, -4.7, -4.5, -4.2, -4.0]

    if omega_values is None:
        omega_values = [0.8, 0.9, 1.0, 1.1, 1.2]

    if density_ratios is None:
        density_ratios = [5.0, 10.0, 14.0, 20.0]

    print(f"\n{'='*70}")
    print("PARAMETER SWEEP FOR SPURIOUS CURRENT OPTIMIZATION")
    print(f"{'='*70}")
    print(f"Testing:")
    print(f"  g values: {g_values}")
    print(f"  omega values: {omega_values}")
    print(f"  density ratios: {density_ratios}")
    print(f"  time steps: {time_steps}")
    print(f"  domain size: {N}×{N}")
    print(f"{'='*70}\n")

    results = []
    total_runs = len(g_values) * len(omega_values) * len(density_ratios)
    run_count = 0

    # Sweep over all parameter combinations
    for g in g_values:
        for omega in omega_values:
            for ratio in density_ratios:
                run_count += 1

                # Calculate densities from ratio
                rho_liquid = ratio * 0.15
                rho_gas = 0.15

                print(f"\n{'='*70}")
                print(f"Run {run_count}/{total_runs}")
                print(f"Parameters: g={g:.2f}, ω={omega:.2f}, ratio={ratio:.1f}")
                print(f"{'='*70}")

                try:
                    # Create simulation
                    sim = StaticDroplet(
                        N=N,
                        droplet_radius=droplet_radius,
                        rho_liquid=rho_liquid,
                        rho_gas=rho_gas,
                        omega=omega,
                        g_interaction=g,
                        rho0=1.0,
                        use_mrt=True,
                        target='cpu'
                    )

                    # Initialize
                    sim.initialize()

                    # Run simulation
                    sim.run(time_steps=time_steps, print_interval=5000)

                    # Compute metrics
                    spurious = sim.compute_spurious_currents()
                    pressure_data = sim.compute_laplace_pressure()

                    rho = sim.get_density_field()

                    result = {
                        'g': g,
                        'omega': omega,
                        'density_ratio': ratio,
                        'rho_liquid': rho_liquid,
                        'rho_gas': rho_gas,
                        'spurious_currents': spurious,
                        'delta_P': pressure_data['delta_P'],
                        'P_inside': pressure_data['P_inside'],
                        'P_outside': pressure_data['P_outside'],
                        'rho_min': rho.min(),
                        'rho_max': rho.max(),
                        'success': True
                    }

                    print(f"✓ Success: spurious = {spurious:.6e}, ΔP = {pressure_data['delta_P']:.6f}")

                except Exception as e:
                    print(f"✗ Failed: {str(e)}")
                    result = {
                        'g': g,
                        'omega': omega,
                        'density_ratio': ratio,
                        'rho_liquid': rho_liquid,
                        'rho_gas': rho_gas,
                        'spurious_currents': np.inf,
                        'delta_P': np.nan,
                        'P_inside': np.nan,
                        'P_outside': np.nan,
                        'rho_min': np.nan,
                        'rho_max': np.nan,
                        'success': False,
                        'error': str(e)
                    }

                results.append(result)

    return results


def analyze_and_plot_results(results: list, save_dir: str = "figures/multiphase"):
    """
    Analyze parameter sweep results and create visualization.

    Args:
        results: List of result dictionaries
        save_dir: Directory to save figures
    """
    Path(save_dir).mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*70}")
    print("ANALYZING PARAMETER SWEEP RESULTS")
    print(f"{'='*70}")

    # Convert to arrays for analysis
    successful = [r for r in results if r['success']]

    if len(successful) == 0:
        print("✗ No successful runs!")
        return

    print(f"Successful runs: {len(successful)}/{len(results)}")

    # Find best configuration
    best = min(successful, key=lambda x: x['spurious_currents'])

    print(f"\n{'='*70}")
    print("BEST CONFIGURATION")
    print(f"{'='*70}")
    print(f"g = {best['g']:.2f}")
    print(f"ω = {best['omega']:.2f}")
    print(f"Density ratio = {best['density_ratio']:.1f}")
    print(f"Spurious currents = {best['spurious_currents']:.6e}")
    print(f"ΔP = {best['delta_P']:.6f}")
    print(f"Status: {'✓ EXCELLENT' if best['spurious_currents'] < 1e-3 else '✓ GOOD' if best['spurious_currents'] < 0.01 else '⚠ ACCEPTABLE'}")
    print(f"{'='*70}\n")

    # Create visualization
    fig = plt.figure(figsize=(16, 12))

    # Extract data
    g_vals = np.array([r['g'] for r in successful])
    omega_vals = np.array([r['omega'] for r in successful])
    ratio_vals = np.array([r['density_ratio'] for r in successful])
    spurious_vals = np.array([r['spurious_currents'] for r in successful])
    delta_p_vals = np.array([r['delta_P'] for r in successful])

    # 1. Spurious currents vs g (for each omega)
    ax1 = plt.subplot(2, 3, 1)
    unique_omegas = np.unique(omega_vals)
    for omega in unique_omegas:
        mask = omega_vals == omega
        if np.any(mask):
            g_subset = g_vals[mask]
            spurious_subset = spurious_vals[mask]
            # Sort for plotting
            sort_idx = np.argsort(g_subset)
            ax1.semilogy(g_subset[sort_idx], spurious_subset[sort_idx],
                        'o-', label=f'ω={omega:.2f}', linewidth=2, markersize=6)

    ax1.axhline(1e-3, color='k', linestyle='--', linewidth=1, label='Target: 1e-3')
    ax1.set_xlabel('Interaction Strength g', fontsize=11)
    ax1.set_ylabel('Spurious Currents', fontsize=11)
    ax1.set_title('Effect of g on Spurious Currents', fontsize=12, fontweight='bold')
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    # 2. Spurious currents vs omega (for each g)
    ax2 = plt.subplot(2, 3, 2)
    unique_gs = np.unique(g_vals)
    for g in unique_gs:
        mask = g_vals == g
        if np.any(mask):
            omega_subset = omega_vals[mask]
            spurious_subset = spurious_vals[mask]
            sort_idx = np.argsort(omega_subset)
            ax2.semilogy(omega_subset[sort_idx], spurious_subset[sort_idx],
                        'o-', label=f'g={g:.2f}', linewidth=2, markersize=6)

    ax2.axhline(1e-3, color='k', linestyle='--', linewidth=1, label='Target: 1e-3')
    ax2.set_xlabel('Relaxation Parameter ω', fontsize=11)
    ax2.set_ylabel('Spurious Currents', fontsize=11)
    ax2.set_title('Effect of ω on Spurious Currents', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)

    # 3. Spurious currents vs density ratio
    ax3 = plt.subplot(2, 3, 3)
    unique_ratios = np.unique(ratio_vals)
    for ratio in unique_ratios:
        mask = ratio_vals == ratio
        if np.any(mask):
            g_subset = g_vals[mask]
            spurious_subset = spurious_vals[mask]
            sort_idx = np.argsort(g_subset)
            ax3.semilogy(g_subset[sort_idx], spurious_subset[sort_idx],
                        'o-', label=f'ratio={ratio:.1f}', linewidth=2, markersize=6)

    ax3.axhline(1e-3, color='k', linestyle='--', linewidth=1, label='Target: 1e-3')
    ax3.set_xlabel('Interaction Strength g', fontsize=11)
    ax3.set_ylabel('Spurious Currents', fontsize=11)
    ax3.set_title('Effect of Density Ratio', fontsize=12, fontweight='bold')
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3)

    # 4. Laplace pressure vs g
    ax4 = plt.subplot(2, 3, 4)
    for omega in unique_omegas:
        mask = omega_vals == omega
        if np.any(mask):
            g_subset = g_vals[mask]
            dp_subset = delta_p_vals[mask]
            sort_idx = np.argsort(g_subset)
            ax4.plot(g_subset[sort_idx], dp_subset[sort_idx],
                    'o-', label=f'ω={omega:.2f}', linewidth=2, markersize=6)

    ax4.set_xlabel('Interaction Strength g', fontsize=11)
    ax4.set_ylabel('Laplace Pressure ΔP', fontsize=11)
    ax4.set_title('Effect of g on Surface Tension', fontsize=12, fontweight='bold')
    ax4.legend(fontsize=9)
    ax4.grid(True, alpha=0.3)

    # 5. 2D heatmap: spurious currents vs (g, omega)
    ax5 = plt.subplot(2, 3, 5)

    # Create grid for heatmap
    g_unique = np.unique(g_vals)
    omega_unique = np.unique(omega_vals)

    # Average over density ratios
    spurious_grid = np.zeros((len(g_unique), len(omega_unique)))
    for i, g in enumerate(g_unique):
        for j, omega in enumerate(omega_unique):
            mask = (g_vals == g) & (omega_vals == omega)
            if np.any(mask):
                spurious_grid[i, j] = np.mean(spurious_vals[mask])
            else:
                spurious_grid[i, j] = np.nan

    im = ax5.imshow(np.log10(spurious_grid), cmap='RdYlGn_r', aspect='auto',
                    extent=[omega_unique.min(), omega_unique.max(),
                           g_unique.min(), g_unique.max()],
                    origin='lower')
    ax5.set_xlabel('Relaxation Parameter ω', fontsize=11)
    ax5.set_ylabel('Interaction Strength g', fontsize=11)
    ax5.set_title('log₁₀(Spurious Currents)', fontsize=12, fontweight='bold')
    plt.colorbar(im, ax=ax5, label='log₁₀(spurious)')

    # 6. Best configurations table
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('off')

    # Top 5 configurations
    top5 = sorted(successful, key=lambda x: x['spurious_currents'])[:5]

    table_text = "TOP 5 CONFIGURATIONS\n"
    table_text += "="*45 + "\n"
    for idx, config in enumerate(top5, 1):
        table_text += f"{idx}. g={config['g']:.2f}, ω={config['omega']:.2f}, "
        table_text += f"ratio={config['density_ratio']:.1f}\n"
        table_text += f"   Spurious: {config['spurious_currents']:.6e}\n"
        table_text += f"   ΔP: {config['delta_P']:.6f}\n"
        table_text += "\n"

    ax6.text(0.05, 0.95, table_text, transform=ax6.transAxes,
            fontsize=9, verticalalignment='top',
            fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(f"{save_dir}/parameter_sweep_results.png",
               dpi=300, bbox_inches='tight')
    plt.close()

    print(f"✓ Parameter sweep plot saved to: {save_dir}/parameter_sweep_results.png")

    # Save results to JSON
    output_file = Path("output/multiphase/parameter_sweep_results.json")
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, 'w') as f:
        json.dump({
            'timestamp': datetime.now().isoformat(),
            'best_config': best,
            'all_results': results
        }, f, indent=2, default=str)

    print(f"✓ Results saved to: {output_file}")


if __name__ == "__main__":
    # Run parameter sweep
    print("Starting parameter sweep (this will take a while)...")

    results = run_parameter_sweep(
        g_values=[-5.0, -4.7, -4.5, -4.2],
        omega_values=[0.8, 1.0, 1.2],
        density_ratios=[5.0, 10.0, 14.0],
        time_steps=20000,
        N=64,
        droplet_radius=15.0
    )

    # Analyze results
    analyze_and_plot_results(results)

    print("\n" + "="*70)
    print("Parameter sweep completed!")
    print("="*70)
