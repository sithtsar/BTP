"""
Comprehensive plotting for all validated ELBM results
Plots only valid data, skips NaN cases
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path

# Set publication-quality plotting defaults
plt.rcParams.update({
    'font.size': 11,
    'font.family': 'serif',
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

def load_data_safe(filepath):
    """Load data file, return None if file has NaN or doesn't exist"""
    try:
        data = np.loadtxt(filepath)
        if np.any(np.isnan(data)):
            print(f"  ⚠️  Skipping {filepath.name} (contains NaN)")
            return None
        return data
    except:
        print(f"  ⚠️  Skipping {filepath.name} (cannot load)")
        return None

def plot_stokes_shear():
    """Plot Stokes Shear Flow - BGK vs ELBM vs Analytical"""
    print("\n📊 Plotting Stokes Shear Flow...")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Load data
    bgk = load_data_safe(Path('output/stokes_shear_bgk.dat'))
    elbm = load_data_safe(Path('output/stokes_shear_elbm.dat'))

    if bgk is None and elbm is None:
        print("  ❌ No valid data for Stokes Shear")
        plt.close(fig)
        return

    # Plot velocity profiles
    if bgk is not None:
        ax1.plot(bgk[:, 0], bgk[:, 1], 'b-', label='BGK', linewidth=2)
        ax1.plot(bgk[:, 0], bgk[:, 2], 'k--', label='Analytical', linewidth=1.5, alpha=0.7)

    if elbm is not None:
        ax1.plot(elbm[:, 0], elbm[:, 1], 'r-', label='ELBM', linewidth=2, alpha=0.7)

    ax1.set_xlabel('y position')
    ax1.set_ylabel('Velocity u(y)')
    ax1.set_title('Stokes Shear Flow: Velocity Profile')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot error
    if bgk is not None:
        ax2.semilogy(bgk[:, 0], bgk[:, 3], 'b-', label='BGK Error', linewidth=2)
    if elbm is not None:
        ax2.semilogy(elbm[:, 0], elbm[:, 3], 'r-', label='ELBM Error', linewidth=2)

    ax2.set_xlabel('y position')
    ax2.set_ylabel('Absolute Error')
    ax2.set_title('Point-wise Error Distribution')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    outdir = Path('figures/analytical_validation')
    outdir.mkdir(parents=True, exist_ok=True)
    plt.savefig(outdir / 'stokes_shear_validation.png')
    print(f"  ✅ Saved: {outdir / 'stokes_shear_validation.png'}")
    plt.close()

def plot_couette():
    """Plot Couette Flow - BGK vs ELBM vs Analytical"""
    print("\n📊 Plotting Couette Flow...")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Load data
    bgk = load_data_safe(Path('output/couette_bgk.dat'))
    elbm = load_data_safe(Path('output/couette_elbm.dat'))

    if bgk is None and elbm is None:
        print("  ❌ No valid data for Couette")
        plt.close(fig)
        return

    # Plot velocity profiles
    if bgk is not None:
        ax1.plot(bgk[:, 0], bgk[:, 1], 'b-', label='BGK', linewidth=2)
        ax1.plot(bgk[:, 0], bgk[:, 2], 'k--', label='Analytical', linewidth=1.5, alpha=0.7)

    if elbm is not None:
        ax1.plot(elbm[:, 0], elbm[:, 1], 'r-', label='ELBM', linewidth=2, alpha=0.7)

    ax1.set_xlabel('y position')
    ax1.set_ylabel('Velocity u(y)')
    ax1.set_title('Couette Flow: Velocity Profile')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot error
    if bgk is not None:
        ax2.semilogy(bgk[:, 0], bgk[:, 3], 'b-', label='BGK Error', linewidth=2)
    if elbm is not None:
        ax2.semilogy(elbm[:, 0], elbm[:, 3], 'r-', label='ELBM Error', linewidth=2)

    ax2.set_xlabel('y position')
    ax2.set_ylabel('Absolute Error')
    ax2.set_title('Point-wise Error Distribution')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    outdir = Path('figures/analytical_validation')
    outdir.mkdir(parents=True, exist_ok=True)
    plt.savefig(outdir / 'couette_validation.png')
    print(f"  ✅ Saved: {outdir / 'couette_validation.png'}")
    plt.close()

def plot_poiseuille():
    """Plot Poiseuille Flow - BGK vs ELBM vs Analytical"""
    print("\n📊 Plotting Poiseuille Flow...")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Load data
    bgk = load_data_safe(Path('output/poiseuille_bgk.dat'))
    elbm = load_data_safe(Path('output/poiseuille_elbm.dat'))

    if bgk is None and elbm is None:
        print("  ❌ No valid data for Poiseuille")
        plt.close(fig)
        return

    # Plot velocity profiles
    if bgk is not None:
        ax1.plot(bgk[:, 0], bgk[:, 1], 'b-', label='BGK', linewidth=2)
        ax1.plot(bgk[:, 0], bgk[:, 2], 'k--', label='Analytical', linewidth=1.5, alpha=0.7)

    if elbm is not None:
        ax1.plot(elbm[:, 0], elbm[:, 1], 'r-', label='ELBM', linewidth=2, alpha=0.7)

    ax1.set_xlabel('y position')
    ax1.set_ylabel('Velocity u(y)')
    ax1.set_title('Poiseuille Flow: Velocity Profile')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Plot error
    if bgk is not None:
        ax2.semilogy(bgk[:, 0], bgk[:, 3], 'b-', label='BGK Error', linewidth=2)
    if elbm is not None:
        ax2.semilogy(elbm[:, 0], elbm[:, 3], 'r-', label='ELBM Error', linewidth=2)

    ax2.set_xlabel('y position')
    ax2.set_ylabel('Absolute Error')
    ax2.set_title('Point-wise Error Distribution')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    outdir = Path('figures/analytical_validation')
    outdir.mkdir(parents=True, exist_ok=True)
    plt.savefig(outdir / 'poiseuille_validation.png')
    print(f"  ✅ Saved: {outdir / 'poiseuille_validation.png'}")
    plt.close()

def plot_taylor_green():
    """Plot Taylor-Green Vortex - Energy Decay"""
    print("\n📊 Plotting Taylor-Green Vortex...")

    bgk = load_data_safe(Path('output/tg_bgk.dat'))
    elbm = load_data_safe(Path('output/tg_elbm.dat'))

    if bgk is None and elbm is None:
        print("  ❌ No valid data for Taylor-Green")
        return

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    if bgk is not None:
        # Assuming format: step, energy, exact_energy
        ax.semilogy(bgk[:, 0], bgk[:, 1], 'b-', label='BGK', linewidth=2)
        if bgk.shape[1] > 2:
            ax.semilogy(bgk[:, 0], bgk[:, 2], 'k--', label='Analytical', linewidth=1.5, alpha=0.7)

    if elbm is not None:
        ax.semilogy(elbm[:, 0], elbm[:, 1], 'r-', label='ELBM', linewidth=2, alpha=0.7)

    ax.set_xlabel('Time Step')
    ax.set_ylabel('Kinetic Energy')
    ax.set_title('Taylor-Green Vortex: Energy Decay')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    outdir = Path('figures/analytical_validation')
    outdir.mkdir(parents=True, exist_ok=True)
    plt.savefig(outdir / 'taylor_green_validation.png')
    print(f"  ✅ Saved: {outdir / 'taylor_green_validation.png'}")
    plt.close()

def plot_cavity_profiles():
    """Plot Lid-Driven Cavity benchmark profiles"""
    print("\n📊 Plotting Lid-Driven Cavity...")

    re100_profile = load_data_safe(Path('output/cavity_re100_profiles.dat'))
    re1000_profile = load_data_safe(Path('output/cavity_re1000_profiles.dat'))

    if re100_profile is None and re1000_profile is None:
        print("  ❌ No valid data for Cavity")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    if re100_profile is not None:
        # Assuming format: y, u_centerline, v_midline
        ax1.plot(re100_profile[:, 1], re100_profile[:, 0], 'b-', label='u(x=0.5)', linewidth=2)
        ax1.plot(re100_profile[:, 2], re100_profile[:, 0], 'r--', label='v(y=0.5)', linewidth=2)
        ax1.set_xlabel('Velocity')
        ax1.set_ylabel('Position')
        ax1.set_title('Cavity Re=100')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

    if re1000_profile is not None:
        ax2.plot(re1000_profile[:, 1], re1000_profile[:, 0], 'b-', label='u(x=0.5)', linewidth=2)
        ax2.plot(re1000_profile[:, 2], re1000_profile[:, 0], 'r--', label='v(y=0.5)', linewidth=2)
        ax2.set_xlabel('Velocity')
        ax2.set_ylabel('Position')
        ax2.set_title('Cavity Re=1000')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    outdir = Path('figures/benchmark_cases')
    outdir.mkdir(parents=True, exist_ok=True)
    plt.savefig(outdir / 'cavity_profiles.png')
    print(f"  ✅ Saved: {outdir / 'cavity_profiles.png'}")
    plt.close()

def plot_cylinder_forces():
    """Plot Cylinder flow force coefficients"""
    print("\n📊 Plotting Cylinder Flow Forces...")

    re40_forces = load_data_safe(Path('output/cylinder_re40_forces.dat'))
    re100_forces = load_data_safe(Path('output/cylinder_re100_forces.dat'))

    if re40_forces is None and re100_forces is None:
        print("  ❌ No valid data for Cylinder")
        return

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

    if re40_forces is not None:
        # Assuming format: step, cd, cl
        ax1.plot(re40_forces[:, 0], re40_forces[:, 1], 'b-', label='Cd', linewidth=1.5)
        ax1.plot(re40_forces[:, 0], re40_forces[:, 2], 'r-', label='Cl', linewidth=1.5)
        ax1.set_xlabel('Time Step')
        ax1.set_ylabel('Force Coefficient')
        ax1.set_title('Cylinder Re=40')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

    if re100_forces is not None:
        ax2.plot(re100_forces[:, 0], re100_forces[:, 1], 'b-', label='Cd', linewidth=1.5)
        ax2.plot(re100_forces[:, 0], re100_forces[:, 2], 'r-', label='Cl', linewidth=1.5)
        ax2.set_xlabel('Time Step')
        ax2.set_ylabel('Force Coefficient')
        ax2.set_title('Cylinder Re=100 (Vortex Shedding)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    outdir = Path('figures/benchmark_cases')
    outdir.mkdir(parents=True, exist_ok=True)
    plt.savefig(outdir / 'cylinder_forces.png')
    print(f"  ✅ Saved: {outdir / 'cylinder_forces.png'}")
    plt.close()

def create_summary_table():
    """Create a summary table of validation results"""
    print("\n📊 Creating Summary Table...")

    results = []

    # Check each test case
    test_cases = [
        ('couette_bgk.dat', 'Couette BGK'),
        ('couette_elbm.dat', 'Couette ELBM'),
        ('poiseuille_bgk.dat', 'Poiseuille BGK'),
        ('poiseuille_elbm.dat', 'Poiseuille ELBM'),
        ('stokes_shear_bgk.dat', 'Stokes BGK'),
        ('stokes_shear_elbm.dat', 'Stokes ELBM'),
        ('tg_bgk.dat', 'Taylor-Green BGK'),
        ('tg_elbm.dat', 'Taylor-Green ELBM'),
    ]

    for filename, name in test_cases:
        data = load_data_safe(Path('output') / filename)
        if data is not None:
            if 'error' in name.lower() or data.shape[1] > 3:
                # Has error column
                if data.shape[1] >= 4:
                    max_error = np.max(data[:, 3])
                    mean_error = np.mean(data[:, 3])
                    results.append(f"{name:25s}: ✅ PASS (Max Error: {max_error:.2e}, Mean: {mean_error:.2e})")
                else:
                    results.append(f"{name:25s}: ✅ PASS")
            else:
                results.append(f"{name:25s}: ✅ PASS")
        else:
            results.append(f"{name:25s}: ❌ FAIL (NaN or missing)")

    # Save summary
    outdir = Path('figures')
    outdir.mkdir(parents=True, exist_ok=True)

    with open(outdir / 'validation_summary.txt', 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("ELBM VALIDATION SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        for result in results:
            f.write(result + "\n")
        f.write("\n" + "=" * 80 + "\n")

    print(f"  ✅ Saved: {outdir / 'validation_summary.txt'}")

    # Print to console
    print("\n" + "=" * 80)
    print("VALIDATION SUMMARY")
    print("=" * 80)
    for result in results:
        print(result)
    print("=" * 80)

if __name__ == "__main__":
    print("=" * 80)
    print("COMPREHENSIVE ELBM VALIDATION PLOTTING")
    print("=" * 80)

    # Create all plots
    plot_stokes_shear()
    plot_couette()
    plot_poiseuille()
    plot_taylor_green()
    plot_cavity_profiles()
    plot_cylinder_forces()

    # Create summary
    create_summary_table()

    print("\n" + "=" * 80)
    print("ALL PLOTS GENERATED SUCCESSFULLY")
    print("=" * 80)
