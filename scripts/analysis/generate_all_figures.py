#!/usr/bin/env python3
"""
Generate all figures from Yu thesis comparing BGK and ELBM

Original Figures (3.2-3.19):
Figures 3.2-3.7:   Case 1 (Re~10)
Figures 3.8-3.13:  Case 2 (Re~100)
Figures 3.14-3.19: Case 3 (Re~1000)

Each case has:
- 2 centerline profiles (BGK, ELBM)
- 4 centerplane contours at different timesteps (LBGK/ELBM at 1250, 2500)

Additional Figures:
- Velocity profile evolution (BGK and ELBM for each case)
- H-function entropy evolution and H-theorem validation
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

sys.path.append('plotting')
import plot_utils

def load_velocity_profile(filepath):
    """Load velocity profile data (y ux uy format)"""
    data = np.loadtxt(filepath, comments='#')
    return {
        'y': data[:, 0],
        'ux': data[:, 1],
        'uy': data[:, 2]
    }

def load_entropy_data(filepath):
    """Load entropy history data"""
    data = np.loadtxt(filepath, comments='#')
    return {
        'timestep': data[:, 0],
        'H_total': data[:, 1],
        'H_min': data[:, 2],
        'H_max': data[:, 3],
        'H_mean': data[:, 4],
        'spurious_vel': data[:, 5]
    }

def generate_case_figures(case_name, case_dir, re_value, viscosity, figure_start):
    """Generate all 6 figures for a case"""

    print(f"\n{'='*60}")
    print(f"Generating figures for {case_name} (Re~{re_value})")
    print(f"{'='*60}")

    output_dir = Path("figures") / case_name
    output_dir.mkdir(parents=True, exist_ok=True)

    bgk_dir = Path("output") / f"{case_dir}_bgk"
    elbm_dir = Path("output") / f"{case_dir}_elbm"

    timesteps = [0, 1250, 2500]
    figure_num = figure_start

    # Figure X.Y and X.Y+1: Centerline pressure profiles
    print(f"\nFigure 3.{figure_num}: Centerline Pressure Profile (LBGK)")
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = ['blue', 'green', 'red']
    labels = ['t=0', 't=1250', 't=2500']

    for i, step in enumerate(timesteps):
        bgk_file = bgk_dir / f"bgk_centerline_{step}.dat"
        if bgk_file.exists():
            data = plot_utils.load_centerline_data(bgk_file)
            ax.plot(data['x'], data['p'], color=colors[i],
                   linewidth=2, label=labels[i], marker='o', markersize=3)

    ax.set_xlabel('x (lattice units)', fontsize=14)
    ax.set_ylabel('Pressure (lattice units)', fontsize=14)
    ax.set_title(f'Centerline Pressure Profile (LBGK) - {case_name}\nν={viscosity} m²/s, Re~{re_value}',
                fontsize=14, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / f"figure_3_{figure_num}.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: figure_3_{figure_num}.png")

    figure_num += 1

    print(f"\nFigure 3.{figure_num}: Centerline Pressure Profile (ELBM)")
    fig, ax = plt.subplots(figsize=(10, 6))

    for i, step in enumerate(timesteps):
        elbm_file = elbm_dir / f"elbm_centerline_{step}.dat"
        if elbm_file.exists():
            data = plot_utils.load_centerline_data(elbm_file)
            ax.plot(data['x'], data['p'], color=colors[i],
                   linewidth=2, label=labels[i], marker='s', markersize=3)

    ax.set_xlabel('x (lattice units)', fontsize=14)
    ax.set_ylabel('Pressure (lattice units)', fontsize=14)
    ax.set_title(f'Centerline Pressure Profile (ELBM) - {case_name}\nν={viscosity} m²/s, Re~{re_value}',
                fontsize=14, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / f"figure_3_{figure_num}.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: figure_3_{figure_num}.png")

    figure_num += 1

    # Figures X.Y+2 through X.Y+5: Centerplane contours at different timesteps
    for step in [1250, 2500]:
        for solver, solver_dir in [("LBGK", bgk_dir), ("ELBM", elbm_dir)]:
            print(f"\nFigure 3.{figure_num}: Centerplane Pressure Profile after {step} steps ({solver})")

            data_file = solver_dir / f"{'bgk' if solver == 'LBGK' else 'elbm'}_step_{step}.dat"

            if data_file.exists():
                data = plot_utils.load_full_field_data(data_file)

                # Determine grid dimensions
                nx = len(np.unique(data['x']))
                ny = len(np.unique(data['y']))

                # Reshape to grid
                p_grid = data['p'].reshape((ny, nx))

                fig, ax = plt.subplots(figsize=(12, 4))
                im = ax.contourf(p_grid, levels=20, cmap='viridis')
                ax.set_xlabel('x (lattice units)', fontsize=14)
                ax.set_ylabel('y (lattice units)', fontsize=14)
                ax.set_title(f'Centerplane Pressure Profile after {step} timesteps ({solver})\n{case_name}, ν={viscosity} m²/s, Re~{re_value}',
                           fontsize=14, fontweight='bold')
                cbar = plt.colorbar(im, ax=ax)
                cbar.set_label('Pressure (lattice units)', fontsize=12)
                plt.tight_layout()
                plt.savefig(output_dir / f"figure_3_{figure_num}.png", dpi=300, bbox_inches='tight')
                plt.close()
                print(f"  Saved: figure_3_{figure_num}.png")
            else:
                print(f"  Warning: {data_file} not found, skipping...")

            figure_num += 1

    return figure_num

def plot_velocity_profiles(case_name, case_dir, re_value, viscosity):
    """Generate velocity profile evolution figures for BGK and ELBM"""

    output_dir = Path("figures") / case_name
    output_dir.mkdir(parents=True, exist_ok=True)

    bgk_dir = Path("output") / f"{case_dir}_bgk"
    elbm_dir = Path("output") / f"{case_dir}_elbm"

    timesteps = [0, 1250, 2500]
    colors = ['blue', 'green', 'red']
    labels = ['t=0', 't=1250', 't=2500']

    # BGK velocity profiles
    fig, ax = plt.subplots(figsize=(10, 6))
    for i, step in enumerate(timesteps):
        velocity_file = bgk_dir / f"bgk_velocity_{step}.dat"
        if velocity_file.exists():
            data = load_velocity_profile(velocity_file)
            ax.plot(data['ux'], data['y'], color=colors[i], linewidth=2,
                   label=labels[i], marker='o', markersize=3)

    ax.set_xlabel('Velocity u(y)', fontsize=14)
    ax.set_ylabel('y (lattice units)', fontsize=14)
    ax.set_title(f'Velocity Profile Evolution - BGK\n{case_name}, ν={viscosity} m²/s, Re~{re_value}',
                fontsize=14, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / f"velocity_profiles_{case_name.lower().replace(' ', '_').replace('~', '')}_bgk.png",
               dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: velocity_profiles_{case_name.lower().replace(' ', '_').replace('~', '')}_bgk.png")

    # ELBM velocity profiles
    fig, ax = plt.subplots(figsize=(10, 6))
    for i, step in enumerate(timesteps):
        velocity_file = elbm_dir / f"elbm_velocity_{step}.dat"
        if velocity_file.exists():
            data = load_velocity_profile(velocity_file)
            ax.plot(data['ux'], data['y'], color=colors[i], linewidth=2,
                   label=labels[i], marker='s', markersize=3)

    ax.set_xlabel('Velocity u(y)', fontsize=14)
    ax.set_ylabel('y (lattice units)', fontsize=14)
    ax.set_title(f'Velocity Profile Evolution - ELBM\n{case_name}, ν={viscosity} m²/s, Re~{re_value}',
                fontsize=14, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / f"velocity_profiles_{case_name.lower().replace(' ', '_').replace('~', '')}_elbm.png",
               dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: velocity_profiles_{case_name.lower().replace(' ', '_').replace('~', '')}_elbm.png")

def plot_entropy_evolution(case_name, case_dir, re_value, viscosity):
    """Generate entropy evolution comparison figure"""

    output_dir = Path("figures") / case_name
    output_dir.mkdir(parents=True, exist_ok=True)

    bgk_entropy_file = Path("output") / f"{case_dir}_bgk" / "entropy_bgk.dat"
    elbm_entropy_file = Path("output") / f"{case_dir}_elbm" / "entropy_elbm.dat"

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # BGK entropy evolution
    if bgk_entropy_file.exists():
        bgk_data = load_entropy_data(bgk_entropy_file)
        ax1.plot(bgk_data['timestep'], bgk_data['H_total'], 'b-', linewidth=2, label='H_total')
        ax1.fill_between(bgk_data['timestep'], bgk_data['H_min'], bgk_data['H_max'],
                        alpha=0.3, color='blue', label='H range')
        ax1.set_xlabel('Timestep', fontsize=12)
        ax1.set_ylabel('Entropy H', fontsize=12)
        ax1.set_title(f'H-function Evolution - BGK\n{case_name}, ν={viscosity} m²/s, Re~{re_value}',
                     fontsize=13, fontweight='bold')
        ax1.legend(fontsize=10)
        ax1.grid(True, alpha=0.3)

    # ELBM entropy evolution
    if elbm_entropy_file.exists():
        elbm_data = load_entropy_data(elbm_entropy_file)
        ax2.plot(elbm_data['timestep'], elbm_data['H_total'], 'r-', linewidth=2, label='H_total')
        ax2.fill_between(elbm_data['timestep'], elbm_data['H_min'], elbm_data['H_max'],
                        alpha=0.3, color='red', label='H range')
        ax2.set_xlabel('Timestep', fontsize=12)
        ax2.set_ylabel('Entropy H', fontsize=12)
        ax2.set_title(f'H-function Evolution - ELBM\n{case_name}, ν={viscosity} m²/s, Re~{re_value}',
                     fontsize=13, fontweight='bold')
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3)

    plt.suptitle(f'H-theorem Validation: Entropy Evolution\n{case_name} (Re~{re_value})',
                fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(output_dir / f"entropy_evolution_{case_name.lower().replace(' ', '_').replace('~', '')}.png",
               dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: entropy_evolution_{case_name.lower().replace(' ', '_').replace('~', '')}.png")

def generate_comparison_summary():
    """Generate a summary comparison figure across all cases"""
    print(f"\n{'='*60}")
    print("Generating Summary Comparison Figure")
    print(f"{'='*60}")

    fig, axes = plt.subplots(3, 2, figsize=(16, 12))

    cases = [
        ("case1_re10", 10, 0.1, "Case 1: Re~10"),
        ("case2_re100", 100, 0.01, "Case 2: Re~100"),
        ("case3_re1000", 1000, 0.001, "Case 3: Re~1000")
    ]

    for i, (case_dir, re, viscosity, title) in enumerate(cases):
        # Final timestep (2500)
        bgk_file = Path("output") / f"{case_dir}_bgk" / "bgk_centerline_2500.dat"
        elbm_file = Path("output") / f"{case_dir}_elbm" / "elbm_centerline_2500.dat"

        if bgk_file.exists() and elbm_file.exists():
            bgk_data = plot_utils.load_centerline_data(bgk_file)
            elbm_data = plot_utils.load_centerline_data(elbm_file)

            # BGK
            ax = axes[i, 0]
            ax.plot(bgk_data['x'], bgk_data['p'], 'b-', linewidth=2)
            ax.set_xlabel('x (lattice units)', fontsize=12)
            ax.set_ylabel('Pressure', fontsize=12)
            ax.set_title(f'{title} - LBGK (t=2500)', fontweight='bold', fontsize=12)
            ax.grid(True, alpha=0.3)

            # ELBM
            ax = axes[i, 1]
            ax.plot(elbm_data['x'], elbm_data['p'], 'r-', linewidth=2)
            ax.set_xlabel('x (lattice units)', fontsize=12)
            ax.set_ylabel('Pressure', fontsize=12)
            ax.set_title(f'{title} - ELBM (t=2500)', fontweight='bold', fontsize=12)
            ax.grid(True, alpha=0.3)

    plt.suptitle('BGK vs ELBM Stability Comparison\nFinal Timestep (t=2500)',
                fontsize=16, fontweight='bold')
    plt.tight_layout()

    output_dir = Path("figures")
    output_dir.mkdir(exist_ok=True)
    plt.savefig(output_dir / "summary_comparison.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("  Saved: summary_comparison.png")

def main():
    print("="*60)
    print("Generating All 19 Figures from Yu Thesis")
    print("="*60)

    # Create figures directory
    Path("figures").mkdir(exist_ok=True)

    # Case 1: Re~10 (Figures 3.2-3.7)
    figure_num = generate_case_figures(
        case_name="Case 1 Re~10",
        case_dir="case1_re10",
        re_value=10,
        viscosity=0.1,
        figure_start=2
    )

    # Velocity profiles and entropy for Case 1
    plot_velocity_profiles("Case 1 Re~10", "case1_re10", 10, 0.1)
    plot_entropy_evolution("Case 1 Re~10", "case1_re10", 10, 0.1)

    # Case 2: Re~100 (Figures 3.8-3.13)
    figure_num = generate_case_figures(
        case_name="Case 2 Re~100",
        case_dir="case2_re100",
        re_value=100,
        viscosity=0.01,
        figure_start=figure_num
    )

    # Velocity profiles and entropy for Case 2
    plot_velocity_profiles("Case 2 Re~100", "case2_re100", 100, 0.01)
    plot_entropy_evolution("Case 2 Re~100", "case2_re100", 100, 0.01)

    # Case 3: Re~1000 (Figures 3.14-3.19)
    figure_num = generate_case_figures(
        case_name="Case 3 Re~1000",
        case_dir="case3_re1000",
        re_value=1000,
        viscosity=0.001,
        figure_start=figure_num
    )

    # Velocity profiles and entropy for Case 3
    plot_velocity_profiles("Case 3 Re~1000", "case3_re1000", 1000, 0.001)
    plot_entropy_evolution("Case 3 Re~1000", "case3_re1000", 1000, 0.001)

    # Generate summary comparison
    generate_comparison_summary()

    print("\n" + "="*60)
    print("All figures generated successfully!")
    print(f"Original figures: {figure_num - 2 + 1}")
    print(f"Additional figures: 9 (3 velocity profile pairs + 3 entropy plots)")
    print(f"Total figures: {figure_num - 2 + 1 + 9}")
    print("Output directory: figures/")
    print("="*60)

if __name__ == "__main__":
    main()
