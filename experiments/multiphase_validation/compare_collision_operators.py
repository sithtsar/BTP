"""
Compare Collision Operators: MRT vs SRT

Direct comparison of Multiple Relaxation Time (MRT) and Single Relaxation Time (SRT/BGK)
collision operators for Shan-Chen multi-phase simulations.

Metrics compared:
1. Spurious currents (lower is better)
2. Numerical stability
3. Computational performance
4. Interface quality
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import time

sys.path.append(str(Path(__file__).parent))

from static_droplet_mrt import StaticDroplet


def run_comparison(
    N: int = 64,
    droplet_radius: float = 15.0,
    time_steps: int = 15000,
    omega: float = 1.0,
    g_interaction: float = -4.7
):
    """
    Run MRT vs SRT comparison on static droplet.

    Args:
        N: Domain size
        droplet_radius: Droplet radius
        time_steps: Number of simulation steps
        omega: Relaxation parameter
        g_interaction: Shan-Chen interaction strength

    Returns:
        Dictionary with comparison results
    """
    print(f"\n{'='*70}")
    print("COLLISION OPERATOR COMPARISON: MRT vs SRT")
    print(f"{'='*70}\n")

    results = {}

    # Test 1: MRT
    print(f"{'='*70}")
    print("TEST 1: MRT (Multiple Relaxation Time)")
    print(f"{'='*70}\n")

    try:
        start_time = time.time()

        sim_mrt = StaticDroplet(
            N=N,
            droplet_radius=droplet_radius,
            rho_liquid=2.1,
            rho_gas=0.15,
            omega=omega,
            g_interaction=g_interaction,
            use_mrt=True,
            target='cpu'
        )

        sim_mrt.initialize()
        sim_mrt.run(time_steps=time_steps, print_interval=3000)

        mrt_time = time.time() - start_time

        # Compute metrics
        spurious_mrt = sim_mrt.compute_spurious_currents()
        pressure_mrt = sim_mrt.compute_laplace_pressure()
        rho_mrt = sim_mrt.get_density_field()

        results['MRT'] = {
            'spurious_currents': spurious_mrt,
            'delta_P': pressure_mrt['delta_P'],
            'P_inside': pressure_mrt['P_inside'],
            'P_outside': pressure_mrt['P_outside'],
            'wall_time': mrt_time,
            'rho_min': rho_mrt.min(),
            'rho_max': rho_mrt.max(),
            'success': True,
            'simulation': sim_mrt
        }

        print(f"\n✓ MRT Results:")
        print(f"  Spurious currents: {spurious_mrt:.6e}")
        print(f"  Laplace pressure:  {pressure_mrt['delta_P']:.6f}")
        print(f"  Wall time:         {mrt_time:.2f} seconds")

    except Exception as e:
        print(f"\n✗ MRT failed: {str(e)}")
        results['MRT'] = {'success': False, 'error': str(e)}

    # Test 2: SRT
    print(f"\n{'='*70}")
    print("TEST 2: SRT (Single Relaxation Time / BGK)")
    print(f"{'='*70}\n")

    try:
        start_time = time.time()

        sim_srt = StaticDroplet(
            N=N,
            droplet_radius=droplet_radius,
            rho_liquid=2.1,
            rho_gas=0.15,
            omega=omega,
            g_interaction=g_interaction,
            use_mrt=False,  # Use SRT
            target='cpu'
        )

        sim_srt.initialize()
        sim_srt.run(time_steps=time_steps, print_interval=3000)

        srt_time = time.time() - start_time

        # Compute metrics
        spurious_srt = sim_srt.compute_spurious_currents()
        pressure_srt = sim_srt.compute_laplace_pressure()
        rho_srt = sim_srt.get_density_field()

        results['SRT'] = {
            'spurious_currents': spurious_srt,
            'delta_P': pressure_srt['delta_P'],
            'P_inside': pressure_srt['P_inside'],
            'P_outside': pressure_srt['P_outside'],
            'wall_time': srt_time,
            'rho_min': rho_srt.min(),
            'rho_max': rho_srt.max(),
            'success': True,
            'simulation': sim_srt
        }

        print(f"\n✓ SRT Results:")
        print(f"  Spurious currents: {spurious_srt:.6e}")
        print(f"  Laplace pressure:  {pressure_srt['delta_P']:.6f}")
        print(f"  Wall time:         {srt_time:.2f} seconds")

    except Exception as e:
        print(f"\n✗ SRT failed: {str(e)}")
        results['SRT'] = {'success': False, 'error': str(e)}

    return results


def plot_comparison(results: dict, save_dir: str = "figures/multiphase"):
    """
    Create comparison visualization.

    Args:
        results: Dictionary from run_comparison()
        save_dir: Directory to save figures
    """
    Path(save_dir).mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*70}")
    print("GENERATING COMPARISON PLOTS")
    print(f"{'='*70}")

    if not results['MRT']['success'] or not results['SRT']['success']:
        print("✗ Cannot create comparison (one or both methods failed)")
        return

    # Create comprehensive comparison figure
    fig = plt.figure(figsize=(16, 10))

    # 1. Density field comparison
    ax1 = plt.subplot(2, 3, 1)
    rho_mrt = results['MRT']['simulation'].get_density_field()
    im1 = ax1.imshow(rho_mrt.T, origin='lower', cmap='viridis')
    ax1.set_title('MRT: Density Field', fontsize=12, fontweight='bold')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    plt.colorbar(im1, ax=ax1, label='ρ')

    ax2 = plt.subplot(2, 3, 2)
    rho_srt = results['SRT']['simulation'].get_density_field()
    im2 = ax2.imshow(rho_srt.T, origin='lower', cmap='viridis')
    ax2.set_title('SRT: Density Field', fontsize=12, fontweight='bold')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    plt.colorbar(im2, ax=ax2, label='ρ')

    # 2. Difference map
    ax3 = plt.subplot(2, 3, 3)
    diff = np.abs(rho_mrt - rho_srt)
    im3 = ax3.imshow(diff.T, origin='lower', cmap='hot')
    ax3.set_title(f'|Difference| (max={diff.max():.4f})', fontsize=12, fontweight='bold')
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    plt.colorbar(im3, ax=ax3, label='|ρ_MRT - ρ_SRT|')

    # 3. Radial profiles
    ax4 = plt.subplot(2, 3, 4)
    r_mrt, rho_avg_mrt, rho_std_mrt = results['MRT']['simulation'].compute_radial_density_profile()
    r_srt, rho_avg_srt, rho_std_srt = results['SRT']['simulation'].compute_radial_density_profile()

    ax4.plot(r_mrt, rho_avg_mrt, 'b-', linewidth=2, label='MRT')
    ax4.plot(r_srt, rho_avg_srt, 'r--', linewidth=2, label='SRT')
    ax4.fill_between(r_mrt, rho_avg_mrt - rho_std_mrt, rho_avg_mrt + rho_std_mrt,
                     alpha=0.2, color='blue')
    ax4.fill_between(r_srt, rho_avg_srt - rho_std_srt, rho_avg_srt + rho_std_srt,
                     alpha=0.2, color='red')

    ax4.set_xlabel('Radial Distance r', fontsize=11)
    ax4.set_ylabel('Density ρ', fontsize=11)
    ax4.set_title('Radial Density Profiles', fontsize=12, fontweight='bold')
    ax4.legend(fontsize=10)
    ax4.grid(True, alpha=0.3)

    # 4. Metrics comparison (bar chart)
    ax5 = plt.subplot(2, 3, 5)

    metrics = ['Spurious\nCurrents', 'Laplace\nPressure ΔP', 'Wall Time\n(seconds)']
    mrt_vals = [
        results['MRT']['spurious_currents'],
        results['MRT']['delta_P'],
        results['MRT']['wall_time']
    ]
    srt_vals = [
        results['SRT']['spurious_currents'],
        results['SRT']['delta_P'],
        results['SRT']['wall_time']
    ]

    x = np.arange(len(metrics))
    width = 0.35

    # Normalize for visualization (different scales)
    mrt_norm = [mrt_vals[0] / max(mrt_vals[0], srt_vals[0]),
                mrt_vals[1] / max(mrt_vals[1], srt_vals[1]),
                mrt_vals[2] / max(mrt_vals[2], srt_vals[2])]
    srt_norm = [srt_vals[0] / max(mrt_vals[0], srt_vals[0]),
                srt_vals[1] / max(mrt_vals[1], srt_vals[1]),
                srt_vals[2] / max(mrt_vals[2], srt_vals[2])]

    rects1 = ax5.bar(x - width/2, mrt_norm, width, label='MRT', color='blue', alpha=0.7)
    rects2 = ax5.bar(x + width/2, srt_norm, width, label='SRT', color='red', alpha=0.7)

    ax5.set_ylabel('Normalized Value', fontsize=11)
    ax5.set_title('Performance Metrics (Normalized)', fontsize=12, fontweight='bold')
    ax5.set_xticks(x)
    ax5.set_xticklabels(metrics, fontsize=9)
    ax5.legend(fontsize=10)
    ax5.grid(True, alpha=0.3, axis='y')

    # 5. Summary table
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('off')

    table_text = "COMPARISON SUMMARY\n"
    table_text += "="*50 + "\n\n"

    table_text += "Spurious Currents:\n"
    table_text += f"  MRT: {results['MRT']['spurious_currents']:.6e}\n"
    table_text += f"  SRT: {results['SRT']['spurious_currents']:.6e}\n"
    if results['MRT']['spurious_currents'] < results['SRT']['spurious_currents']:
        improvement = (results['SRT']['spurious_currents'] - results['MRT']['spurious_currents']) / results['SRT']['spurious_currents'] * 100
        table_text += f"  ✓ MRT is {improvement:.1f}% better\n\n"
    else:
        improvement = (results['MRT']['spurious_currents'] - results['SRT']['spurious_currents']) / results['MRT']['spurious_currents'] * 100
        table_text += f"  ✓ SRT is {improvement:.1f}% better\n\n"

    table_text += "Laplace Pressure ΔP:\n"
    table_text += f"  MRT: {results['MRT']['delta_P']:.6f}\n"
    table_text += f"  SRT: {results['SRT']['delta_P']:.6f}\n"
    table_text += f"  Difference: {abs(results['MRT']['delta_P'] - results['SRT']['delta_P']):.6f}\n\n"

    table_text += "Computational Performance:\n"
    table_text += f"  MRT: {results['MRT']['wall_time']:.2f} seconds\n"
    table_text += f"  SRT: {results['SRT']['wall_time']:.2f} seconds\n"
    speedup = results['MRT']['wall_time'] / results['SRT']['wall_time']
    if speedup > 1:
        table_text += f"  SRT is {speedup:.2f}× faster\n\n"
    else:
        table_text += f"  MRT is {1/speedup:.2f}× faster\n\n"

    table_text += "Recommendation:\n"
    if results['MRT']['spurious_currents'] < 0.01 and results['MRT']['spurious_currents'] < results['SRT']['spurious_currents']:
        table_text += "  ✓ Use MRT for better accuracy\n"
    elif results['SRT']['wall_time'] < results['MRT']['wall_time'] * 0.7 and results['SRT']['spurious_currents'] < 0.05:
        table_text += "  ✓ Use SRT for faster computation\n"
    else:
        table_text += "  ⚠ Trade-off: MRT more stable,\n    SRT potentially faster\n"

    ax6.text(0.05, 0.95, table_text, transform=ax6.transAxes,
            fontsize=9, verticalalignment='top',
            fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(f"{save_dir}/collision_operator_comparison.png",
               dpi=300, bbox_inches='tight')
    plt.close()

    print(f"✓ Comparison plot saved to: {save_dir}/collision_operator_comparison.png")


if __name__ == "__main__":
    # Run comparison
    results = run_comparison(
        N=64,
        droplet_radius=15.0,
        time_steps=15000,
        omega=1.0,
        g_interaction=-4.7
    )

    # Plot comparison
    plot_comparison(results)

    print("\n" + "="*70)
    print("Collision operator comparison completed!")
    print("="*70)
