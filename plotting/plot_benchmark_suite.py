#!/usr/bin/env python3
"""
Benchmark Suite Results Plotting
LBGK vs ELBM comparison for Couette, Poiseuille, and TGV
Publication-quality figures for fluid dynamics validation
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def load_benchmark_data(filename):
    """Load benchmark result data"""
    if not Path(filename).exists():
        return None
    return np.loadtxt(filename, comments='#')

def plot_couette_comparison():
    """Plot Couette flow BGK vs ELBM comparison"""
    # Load data
    bgk_data = load_benchmark_data('output/benchmark_couette_bgk.dat')
    elbm_data = load_benchmark_data('output/benchmark_couette_elbm.dat')

    if bgk_data is None:
        print("Warning: Couette BGK data not found")
        return None

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    y = bgk_data[:, 0]
    u_bgk = bgk_data[:, 1]
    u_exact = bgk_data[:, 2]
    error_bgk = bgk_data[:, 3]

    # Velocity profiles comparison
    ax1.plot(u_exact, y, 'k-', label='Analytical', linewidth=3, alpha=0.8)
    ax1.plot(u_bgk, y, 'b-o', label='BGK', markersize=4, linewidth=2, alpha=0.8)

    error_elbm = None
    if elbm_data is not None:
        u_elbm = elbm_data[:, 1]
        error_elbm = elbm_data[:, 3]
        ax1.plot(u_elbm, y, 'r-s', label='ELBM', markersize=4, linewidth=2, alpha=0.8)

    ax1.set_xlabel('Velocity u_x', fontsize=14, fontweight='bold')
    ax1.set_ylabel('y position', fontsize=14, fontweight='bold')
    ax1.set_title('Couette Flow: Velocity Profiles', fontsize=16, fontweight='bold')
    ax1.legend(fontsize=12, framealpha=0.9)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(-0.01, 0.11)

    # BGK error
    ax2.plot(error_bgk, y, 'b-', linewidth=2.5, label='BGK Error')
    ax2.fill_betweenx(y, 0, error_bgk, alpha=0.3, color='blue')
    ax2.set_xlabel('Absolute Error', fontsize=14, fontweight='bold')
    ax2.set_ylabel('y position', fontsize=14, fontweight='bold')
    ax2.set_title('BGK Error Distribution', fontsize=16, fontweight='bold')
    ax2.grid(True, alpha=0.3)

    # ELBM error (if available)
    if elbm_data is not None and error_elbm is not None:
        ax3.plot(error_elbm, y, 'r-', linewidth=2.5, label='ELBM Error')
        ax3.fill_betweenx(y, 0, error_elbm, alpha=0.3, color='red')
        ax3.set_xlabel('Absolute Error', fontsize=14, fontweight='bold')
        ax3.set_ylabel('y position', fontsize=14, fontweight='bold')
        ax3.set_title('ELBM Error Distribution', fontsize=16, fontweight='bold')
        ax3.grid(True, alpha=0.3)

        # Error comparison
        ax4.plot(error_bgk, y, 'b-', linewidth=2, label='BGK', alpha=0.8)
        ax4.plot(error_elbm, y, 'r-', linewidth=2, label='ELBM', alpha=0.8)
        ax4.set_xlabel('Absolute Error', fontsize=14, fontweight='bold')
        ax4.set_ylabel('y position', fontsize=14, fontweight='bold')
        ax4.set_title('Error Comparison', fontsize=16, fontweight='bold')
        ax4.legend(fontsize=12)
        ax4.grid(True, alpha=0.3)
    else:
        ax3.text(0.5, 0.5, 'ELBM data not available', transform=ax3.transAxes,
                ha='center', va='center', fontsize=14)
        ax4.text(0.5, 0.5, 'ELBM data not available', transform=ax4.transAxes,
                ha='center', va='center', fontsize=14)

    plt.tight_layout()
    return fig

def plot_poiseuille_comparison():
    """Plot Poiseuille flow BGK vs ELBM comparison"""
    # Load data
    bgk_data = load_benchmark_data('output/benchmark_poiseuille_bgk.dat')
    elbm_data = load_benchmark_data('output/benchmark_poiseuille_elbm.dat')

    if bgk_data is None:
        print("Warning: Poiseuille BGK data not found")
        return None

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    y = bgk_data[:, 0]
    u_bgk = bgk_data[:, 1]
    u_exact = bgk_data[:, 2]
    error_bgk = bgk_data[:, 3]

    # Velocity profiles comparison
    ax1.plot(u_exact, y, 'k-', label='Analytical', linewidth=3, alpha=0.8)
    ax1.plot(u_bgk, y, 'b-o', label='BGK', markersize=4, linewidth=2, alpha=0.8)

    error_elbm = None
    if elbm_data is not None:
        u_elbm = elbm_data[:, 1]
        error_elbm = elbm_data[:, 3]
        ax1.plot(u_elbm, y, 'r-s', label='ELBM', markersize=4, linewidth=2, alpha=0.8)

    ax1.set_xlabel('Velocity u_x', fontsize=14, fontweight='bold')
    ax1.set_ylabel('y position', fontsize=14, fontweight='bold')
    ax1.set_title('Poiseuille Flow: Velocity Profiles', fontsize=16, fontweight='bold')
    ax1.legend(fontsize=12, framealpha=0.9)
    ax1.grid(True, alpha=0.3)

    # BGK error
    ax2.plot(error_bgk, y, 'b-', linewidth=2.5, label='BGK Error')
    ax2.fill_betweenx(y, 0, error_bgk, alpha=0.3, color='blue')
    ax2.set_xlabel('Absolute Error', fontsize=14, fontweight='bold')
    ax2.set_ylabel('y position', fontsize=14, fontweight='bold')
    ax2.set_title('BGK Error Distribution', fontsize=16, fontweight='bold')
    ax2.grid(True, alpha=0.3)

    # ELBM error (if available)
    if elbm_data is not None and error_elbm is not None:
        ax3.plot(error_elbm, y, 'r-', linewidth=2.5, label='ELBM Error')
        ax3.fill_betweenx(y, 0, error_elbm, alpha=0.3, color='red')
        ax3.set_xlabel('Absolute Error', fontsize=14, fontweight='bold')
        ax3.set_ylabel('y position', fontsize=14, fontweight='bold')
        ax3.set_title('ELBM Error Distribution', fontsize=16, fontweight='bold')
        ax3.grid(True, alpha=0.3)

        # Error comparison
        ax4.plot(error_bgk, y, 'b-', linewidth=2, label='BGK', alpha=0.8)
        ax4.plot(error_elbm, y, 'r-', linewidth=2, label='ELBM', alpha=0.8)
        ax4.set_xlabel('Absolute Error', fontsize=14, fontweight='bold')
        ax4.set_ylabel('y position', fontsize=14, fontweight='bold')
        ax4.set_title('Error Comparison', fontsize=16, fontweight='bold')
        ax4.legend(fontsize=12)
        ax4.grid(True, alpha=0.3)
    else:
        ax3.text(0.5, 0.5, 'ELBM data not available', transform=ax3.transAxes,
                ha='center', va='center', fontsize=14)
        ax4.text(0.5, 0.5, 'ELBM data not available', transform=ax4.transAxes,
                ha='center', va='center', fontsize=14)

    plt.tight_layout()
    return fig

def plot_tgv_comparison():
    """Plot TGV velocity field comparison"""
    # Load data
    bgk_data = load_benchmark_data('output/benchmark_tgv_bgk.dat')
    elbm_data = load_benchmark_data('output/benchmark_tgv_elbm.dat')

    if bgk_data is None:
        print("Warning: TGV BGK data not found")
        return None

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    # Reshape data to grid
    nx, ny = 100, 100  # From benchmark constants
    x = bgk_data[:, 0].reshape((ny, nx))
    y = bgk_data[:, 1].reshape((ny, nx))

    # BGK velocity magnitude
    ux_bgk = bgk_data[:, 2].reshape((ny, nx))
    uy_bgk = bgk_data[:, 3].reshape((ny, nx))
    vel_mag_bgk = np.sqrt(ux_bgk**2 + uy_bgk**2)

    im1 = ax1.contourf(x, y, vel_mag_bgk, levels=20, cmap='viridis')
    ax1.set_xlabel('x', fontsize=14, fontweight='bold')
    ax1.set_ylabel('y', fontsize=14, fontweight='bold')
    ax1.set_title('TGV: BGK Velocity Magnitude', fontsize=16, fontweight='bold')
    plt.colorbar(im1, ax=ax1, label='|u|')

    # ELBM velocity magnitude (if available)
    if elbm_data is not None:
        ux_elbm = elbm_data[:, 2].reshape((ny, nx))
        uy_elbm = elbm_data[:, 3].reshape((ny, nx))
        vel_mag_elbm = np.sqrt(ux_elbm**2 + uy_elbm**2)

        im2 = ax2.contourf(x, y, vel_mag_elbm, levels=20, cmap='plasma')
        ax2.set_xlabel('x', fontsize=14, fontweight='bold')
        ax2.set_ylabel('y', fontsize=14, fontweight='bold')
        ax2.set_title('TGV: ELBM Velocity Magnitude', fontsize=16, fontweight='bold')
        plt.colorbar(im2, ax=ax2, label='|u|')

        # Error field
        ux_exact = elbm_data[:, 4].reshape((ny, nx))
        uy_exact = elbm_data[:, 5].reshape((ny, nx))
        vel_mag_exact = np.sqrt(ux_exact**2 + uy_exact**2)
        error_field = np.abs(vel_mag_elbm - vel_mag_exact)

        im3 = ax3.contourf(x, y, error_field, levels=20, cmap='RdYlBu_r')
        ax3.set_xlabel('x', fontsize=14, fontweight='bold')
        ax3.set_ylabel('y', fontsize=14, fontweight='bold')
        ax3.set_title('TGV: ELBM Velocity Error', fontsize=16, fontweight='bold')
        plt.colorbar(im3, ax=ax3, label='Error')

        # Velocity vectors (subsampled)
        stride = 8
        x_vec = x[::stride, ::stride]
        y_vec = y[::stride, ::stride]
        ux_vec = ux_elbm[::stride, ::stride]
        uy_vec = uy_elbm[::stride, ::stride]

        ax4.quiver(x_vec, y_vec, ux_vec, uy_vec, scale=5, alpha=0.7)
        ax4.set_xlabel('x', fontsize=14, fontweight='bold')
        ax4.set_ylabel('y', fontsize=14, fontweight='bold')
        ax4.set_title('TGV: ELBM Velocity Field', fontsize=16, fontweight='bold')
    else:
        for ax in [ax2, ax3, ax4]:
            ax.text(0.5, 0.5, 'ELBM data not available', transform=ax.transAxes,
                   ha='center', va='center', fontsize=14)

    plt.tight_layout()
    return fig

def create_results_summary():
    """Create summary table of benchmark results"""
    results = {}

    # Load and compute L2 errors
    test_cases = ['couette', 'poiseuille', 'tgv']

    for case in test_cases:
        for method in ['bgk', 'elbm']:
            filename = f'output/benchmark_{case}_{method}.dat'
            data = load_benchmark_data(filename)
            if data is not None:
                if case == 'tgv':
                    # For TGV, compute L2 error from error columns
                    ux_error = data[:, 6]  # ux_error column
                    uy_error = data[:, 7]  # uy_error column
                    l2_error = np.sqrt(np.mean(ux_error**2 + uy_error**2))
                else:
                    # For Couette/Poiseuille, use error column
                    l2_error = np.sqrt(np.mean(data[:, 3]**2))

                results[f'{case}_{method}'] = l2_error
            else:
                results[f'{case}_{method}'] = None

    # Create summary plot
    fig, ax = plt.subplots(figsize=(12, 8))

    cases = ['Couette Flow', 'Poiseuille Flow', 'Taylor-Green Vortex']
    methods = ['BGK', 'ELBM']
    colors = ['blue', 'red']

    x = np.arange(len(cases))
    width = 0.35

    for i, method in enumerate(['bgk', 'elbm']):
        values = []
        for case in ['couette', 'poiseuille', 'tgv']:
            key = f'{case}_{method}'
            val = results.get(key)
            values.append(val if val is not None else 0)

        bars = ax.bar(x + i*width, values, width, label=methods[i],
                     color=colors[i], alpha=0.7)

        # Add value labels on bars
        for bar, val in zip(bars, values):
            if val > 0:
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{val:.2e}', ha='center', va='bottom', fontsize=10)

    ax.set_xlabel('Test Case', fontsize=14, fontweight='bold')
    ax.set_ylabel('L2 Error Norm', fontsize=14, fontweight='bold')
    ax.set_title('Benchmark Suite: L2 Error Comparison', fontsize=16, fontweight='bold')
    ax.set_xticks(x + width/2)
    ax.set_xticklabels(cases, fontsize=12)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')

    # Add table below plot
    table_data = []
    for case in cases:
        row = [case]
        for method in ['bgk', 'elbm']:
            key = f'{case.lower().split()[0]}_{method}'
            val = results.get(key)
            row.append(f'{val:.2e}' if val is not None else 'N/A')
        table_data.append(row)

    # Create table as separate figure to avoid bbox issues
    fig_table, ax_table = plt.subplots(figsize=(10, 4))
    ax_table.axis('off')
    table = ax_table.table(cellText=table_data,
                          colLabels=['Test Case', 'BGK L2 Error', 'ELBM L2 Error'],
                          loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.2, 1.5)
    ax_table.set_title('Benchmark Results Summary', fontsize=16, fontweight='bold', pad=20)

    return fig, fig_table

def main():
    """Generate all benchmark comparison plots"""
    print("Generating benchmark suite plots...")

    # Create output directory
    output_dir = Path("figures/benchmark_suite")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate individual plots
    print("Creating Couette flow comparison...")
    couette_fig = plot_couette_comparison()
    if couette_fig:
        couette_fig.savefig(output_dir / "couette_comparison.png", dpi=300, bbox_inches='tight')
        plt.close(couette_fig)

    print("Creating Poiseuille flow comparison...")
    poiseuille_fig = plot_poiseuille_comparison()
    if poiseuille_fig:
        poiseuille_fig.savefig(output_dir / "poiseuille_comparison.png", dpi=300, bbox_inches='tight')
        plt.close(poiseuille_fig)

    print("Creating TGV comparison...")
    tgv_fig = plot_tgv_comparison()
    if tgv_fig:
        tgv_fig.savefig(output_dir / "tgv_comparison.png", dpi=300, bbox_inches='tight')
        plt.close(tgv_fig)

    print("Creating results summary...")
    summary_fig, table_fig = create_results_summary()
    summary_fig.savefig(output_dir / "benchmark_summary.png", dpi=300, bbox_inches='tight')
    table_fig.savefig(output_dir / "benchmark_table.png", dpi=300, bbox_inches='tight')
    plt.close(summary_fig)
    plt.close(table_fig)

    print(f"\nBenchmark plots saved to: {output_dir}/")
    print("Generated files:")
    print("  - couette_comparison.png")
    print("  - poiseuille_comparison.png")
    print("  - tgv_comparison.png")
    print("  - benchmark_summary.png")
    print("  - benchmark_table.png")

if __name__ == "__main__":
    main()