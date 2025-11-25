"""
Plot lbmpy Validation Results
==============================
Generates comparison plots for Couette, Poiseuille, and Taylor-Green vortex
comparing C++ and lbmpy implementations.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Add parent directory to path for imports
sys.path.append(str(Path(__file__).parent.parent.parent))


def plot_couette_comparison():
    """Plot Couette flow C++ vs lbmpy comparison."""
    cpp_file = Path("../../output/couette_results.dat")
    lbmpy_file = Path("../../output/couette_lbmpy.dat")

    if not cpp_file.exists() or not lbmpy_file.exists():
        print(f"Warning: Missing data files for Couette flow")
        return

    cpp_data = np.loadtxt(cpp_file)
    lbmpy_data = np.loadtxt(lbmpy_file)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Plot 1: Velocity profiles
    ax = axes[0]
    ax.plot(cpp_data[:, 1], cpp_data[:, 0], 'o-', label='C++ simulation',
            markersize=4, alpha=0.7, linewidth=1.5)
    ax.plot(lbmpy_data[:, 1], lbmpy_data[:, 0], 's-', label='lbmpy simulation',
            markersize=4, alpha=0.7, linewidth=1.5)
    ax.plot(cpp_data[:, 2], cpp_data[:, 0], 'k--', label='Analytical', linewidth=2)
    ax.set_xlabel('Velocity u', fontsize=12)
    ax.set_ylabel('y position', fontsize=12)
    ax.set_title('Couette Flow: Velocity Profiles', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Plot 2: Error comparison
    ax = axes[1]
    cpp_error = np.abs(cpp_data[:, 1] - cpp_data[:, 2])
    lbmpy_error = np.abs(lbmpy_data[:, 1] - lbmpy_data[:, 2])
    ax.semilogy(cpp_error, cpp_data[:, 0], 'o-', label='C++ error',
                markersize=4, alpha=0.7, linewidth=1.5)
    ax.semilogy(lbmpy_error, lbmpy_data[:, 0], 's-', label='lbmpy error',
                markersize=4, alpha=0.7, linewidth=1.5)
    ax.set_xlabel('Absolute Error', fontsize=12)
    ax.set_ylabel('y position', fontsize=12)
    ax.set_title('Error Comparison', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')

    # Plot 3: Direct C++ vs lbmpy comparison
    ax = axes[2]
    ax.plot(cpp_data[:, 1], lbmpy_data[:, 1], 'o', markersize=5, alpha=0.7)
    # Perfect agreement line
    min_val = min(cpp_data[:, 1].min(), lbmpy_data[:, 1].min())
    max_val = max(cpp_data[:, 1].max(), lbmpy_data[:, 1].max())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', linewidth=2, label='Perfect agreement')
    ax.set_xlabel('C++ velocity', fontsize=12)
    ax.set_ylabel('lbmpy velocity', fontsize=12)
    ax.set_title('C++ vs lbmpy Direct Comparison', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.axis('equal')

    # Compute and display L2 errors
    cpp_l2 = np.sqrt(np.mean(cpp_error**2))
    lbmpy_l2 = np.sqrt(np.mean(lbmpy_error**2))
    diff_l2 = np.sqrt(np.mean((cpp_data[:, 1] - lbmpy_data[:, 1])**2))

    fig.suptitle(f'Couette Flow | C++ L2: {cpp_l2:.2e} | lbmpy L2: {lbmpy_l2:.2e} | Difference: {diff_l2:.2e}',
                 fontsize=13, fontweight='bold', y=1.02)

    plt.tight_layout()
    output_file = '../../figures/analytical_validation/couette_lbmpy_comparison.png'
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()


def plot_poiseuille_comparison():
    """Plot Poiseuille flow C++ vs lbmpy comparison."""
    cpp_file = Path("../../output/poiseuille_results.dat")
    lbmpy_file = Path("../../output/poiseuille_lbmpy.dat")

    if not cpp_file.exists() or not lbmpy_file.exists():
        print(f"Warning: Missing data files for Poiseuille flow")
        return

    cpp_data = np.loadtxt(cpp_file)
    lbmpy_data = np.loadtxt(lbmpy_file)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Plot 1: Velocity profiles
    ax = axes[0]
    ax.plot(cpp_data[:, 1], cpp_data[:, 0], 'o-', label='C++ simulation',
            markersize=4, alpha=0.7, linewidth=1.5)
    ax.plot(lbmpy_data[:, 1], lbmpy_data[:, 0], 's-', label='lbmpy simulation',
            markersize=4, alpha=0.7, linewidth=1.5)
    ax.plot(cpp_data[:, 2], cpp_data[:, 0], 'k--', label='Analytical', linewidth=2)
    ax.set_xlabel('Velocity u', fontsize=12)
    ax.set_ylabel('y position', fontsize=12)
    ax.set_title('Poiseuille Flow: Velocity Profiles', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Plot 2: Error comparison
    ax = axes[1]
    cpp_error = np.abs(cpp_data[:, 1] - cpp_data[:, 2])
    lbmpy_error = np.abs(lbmpy_data[:, 1] - lbmpy_data[:, 2])
    ax.semilogy(cpp_error, cpp_data[:, 0], 'o-', label='C++ error',
                markersize=4, alpha=0.7, linewidth=1.5)
    ax.semilogy(lbmpy_error, lbmpy_data[:, 0], 's-', label='lbmpy error',
                markersize=4, alpha=0.7, linewidth=1.5)
    ax.set_xlabel('Absolute Error', fontsize=12)
    ax.set_ylabel('y position', fontsize=12)
    ax.set_title('Error Comparison', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')

    # Plot 3: Direct C++ vs lbmpy comparison
    ax = axes[2]
    ax.plot(cpp_data[:, 1], lbmpy_data[:, 1], 'o', markersize=5, alpha=0.7)
    # Perfect agreement line
    min_val = min(cpp_data[:, 1].min(), lbmpy_data[:, 1].min())
    max_val = max(cpp_data[:, 1].max(), lbmpy_data[:, 1].max())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', linewidth=2, label='Perfect agreement')
    ax.set_xlabel('C++ velocity', fontsize=12)
    ax.set_ylabel('lbmpy velocity', fontsize=12)
    ax.set_title('C++ vs lbmpy Direct Comparison', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.axis('equal')

    # Compute and display L2 errors
    cpp_l2 = np.sqrt(np.mean(cpp_error**2))
    lbmpy_l2 = np.sqrt(np.mean(lbmpy_error**2))
    diff_l2 = np.sqrt(np.mean((cpp_data[:, 1] - lbmpy_data[:, 1])**2))

    fig.suptitle(f'Poiseuille Flow | C++ L2: {cpp_l2:.2e} | lbmpy L2: {lbmpy_l2:.2e} | Difference: {diff_l2:.2e}',
                 fontsize=13, fontweight='bold', y=1.02)

    plt.tight_layout()
    output_file = '../../figures/analytical_validation/poiseuille_lbmpy_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()


def plot_taylor_green_comparison():
    """Plot Taylor-Green vortex C++ vs lbmpy comparison."""
    cpp_file = Path("../../output/tg_energy.dat")
    lbmpy_file = Path("../../output/tg_lbmpy.dat")

    if not cpp_file.exists() or not lbmpy_file.exists():
        print(f"Warning: Missing data files for Taylor-Green vortex")
        return

    cpp_data = np.loadtxt(cpp_file)
    lbmpy_data = np.loadtxt(lbmpy_file)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Plot 1: Energy decay
    ax = axes[0]
    ax.semilogy(cpp_data[:, 1], cpp_data[:, 2], 'o-', label='C++ simulation',
                markersize=5, alpha=0.7, linewidth=1.5)
    ax.semilogy(lbmpy_data[:, 1], lbmpy_data[:, 2], 's-', label='lbmpy simulation',
                markersize=5, alpha=0.7, linewidth=1.5)

    # Analytical decay
    E0_cpp = cpp_data[0, 2]
    t = cpp_data[:, 1]
    nu = 0.01
    k = 2.0 * np.pi / 64
    E_analytical = E0_cpp * np.exp(-4.0 * nu * k**2 * t)
    ax.semilogy(t, E_analytical, 'k--', label='Analytical', linewidth=2)

    ax.set_xlabel('Time steps', fontsize=12)
    ax.set_ylabel('Kinetic Energy', fontsize=12)
    ax.set_title('Taylor-Green Vortex: Energy Decay', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, which='both')

    # Plot 2: Energy difference
    ax = axes[1]
    # Interpolate lbmpy data to C++ time points
    lbmpy_energy_interp = np.interp(cpp_data[:, 1], lbmpy_data[:, 1], lbmpy_data[:, 2])
    energy_diff = np.abs(cpp_data[:, 2] - lbmpy_energy_interp)
    ax.semilogy(cpp_data[:, 1], energy_diff, 'o-', markersize=5, alpha=0.7, linewidth=1.5)
    ax.set_xlabel('Time steps', fontsize=12)
    ax.set_ylabel('|E_cpp - E_lbmpy|', fontsize=12)
    ax.set_title('Energy Difference', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, which='both')

    # Plot 3: Viscosity extraction over time
    ax = axes[2]

    # Extract viscosity at each time point for both
    def extract_viscosity(data):
        E0 = data[0, 2]
        times = data[1:, 1]
        energies = data[1:, 2]
        nu_extracted = -np.log(energies / E0) / (4.0 * k**2 * times)
        return times, nu_extracted

    t_cpp, nu_cpp = extract_viscosity(cpp_data)
    t_lbmpy, nu_lbmpy = extract_viscosity(lbmpy_data)

    ax.plot(t_cpp, nu_cpp, 'o-', label='C++ extracted', markersize=4, alpha=0.7, linewidth=1.5)
    ax.plot(t_lbmpy, nu_lbmpy, 's-', label='lbmpy extracted', markersize=4, alpha=0.7, linewidth=1.5)
    ax.axhline(y=nu, color='k', linestyle='--', linewidth=2, label='Expected')
    ax.set_xlabel('Time steps', fontsize=12)
    ax.set_ylabel('Extracted viscosity', fontsize=12)
    ax.set_title('Viscosity Extraction', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0.008, 0.012])

    # Compute final viscosity errors
    nu_cpp_final = nu_cpp[-1]
    nu_lbmpy_final = nu_lbmpy[-1]
    error_cpp = abs(nu_cpp_final - nu) / nu * 100
    error_lbmpy = abs(nu_lbmpy_final - nu) / nu * 100

    fig.suptitle(f'Taylor-Green Vortex | C++ ν error: {error_cpp:.2f}% | lbmpy ν error: {error_lbmpy:.2f}%',
                 fontsize=13, fontweight='bold', y=1.02)

    plt.tight_layout()
    output_file = '../../figures/analytical_validation/taylor_green_lbmpy_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()


def main():
    """Generate all lbmpy comparison plots."""
    print("="*70)
    print("Generating lbmpy Validation Comparison Plots")
    print("="*70)
    print()

    # Ensure output directory exists
    Path("../../figures/analytical_validation").mkdir(parents=True, exist_ok=True)

    # Generate plots
    plot_couette_comparison()
    plot_poiseuille_comparison()
    plot_taylor_green_comparison()

    print()
    print("="*70)
    print("All lbmpy comparison plots generated!")
    print("Location: figures/analytical_validation/")
    print("="*70)


if __name__ == "__main__":
    main()
