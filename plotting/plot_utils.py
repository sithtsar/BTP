"""
Plotting utilities for ELBM simulation results
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def load_centerline_data(filepath):
    """Load centerline pressure profile data"""
    data = np.loadtxt(filepath, comments='#')
    return {
        'x': data[:, 0],
        'p': data[:, 1],
        'rho': data[:, 2],
        'ux': data[:, 3]
    }

def load_full_field_data(filepath):
    """Load full field data"""
    data = np.loadtxt(filepath, comments='#')
    return {
        'x': data[:, 0],
        'y': data[:, 1],
        'rho': data[:, 2],
        'ux': data[:, 3],
        'uy': data[:, 4],
        'p': data[:, 5]
    }

def plot_centerline_comparison(bgk_data, elbm_data, title, save_path=None):
    """Plot centerline pressure profile comparison"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # BGK
    ax1.plot(bgk_data['x'], bgk_data['p'], 'b-', linewidth=2)
    ax1.set_xlabel('x (lattice units)', fontsize=12)
    ax1.set_ylabel('Pressure (lattice units)', fontsize=12)
    ax1.set_title(f'{title} - BGK', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)

    # ELBM
    ax2.plot(elbm_data['x'], elbm_data['p'], 'r-', linewidth=2)
    ax2.set_xlabel('x (lattice units)', fontsize=12)
    ax2.set_ylabel('Pressure (lattice units)', fontsize=12)
    ax2.set_title(f'{title} - ELBM', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    return fig

def plot_centerplane_contour(data, nx, ny, title, save_path=None):
    """Plot centerplane pressure contour"""
    # Reshape data to grid
    p_grid = data['p'].reshape((ny, nx))

    fig, ax = plt.subplots(figsize=(12, 4))

    im = ax.contourf(p_grid, levels=20, cmap='viridis')
    ax.set_xlabel('x (lattice units)', fontsize=12)
    ax.set_ylabel('y (lattice units)', fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Pressure (lattice units)', fontsize=12)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    return fig

def plot_time_evolution(data_dict, title, ylabel='Pressure', save_path=None):
    """Plot time evolution of centerline profiles"""
    fig, ax = plt.subplots(figsize=(10, 6))

    colors = plt.cm.viridis(np.linspace(0, 1, len(data_dict)))

    for i, (label, data) in enumerate(data_dict.items()):
        ax.plot(data['x'], data['p'], color=colors[i],
               linewidth=2, label=label, alpha=0.8)

    ax.set_xlabel('x (lattice units)', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    return fig

def plot_bgk_vs_elbm_comparison(bgk_files, elbm_files, steps, case_name, save_dir=None):
    """Create comprehensive comparison plots for a test case"""

    figures = []

    # Plot each timestep
    for step in steps:
        bgk_file = bgk_files[step]
        elbm_file = elbm_files[step]

        bgk_data = load_centerline_data(bgk_file)
        elbm_data = load_centerline_data(elbm_file)

        title = f'{case_name} - Step {step}'
        save_path = None
        if save_dir:
            save_path = Path(save_dir) / f'{case_name.lower().replace(" ", "_")}_step_{step}.png'

        fig = plot_centerline_comparison(bgk_data, elbm_data, title, save_path)
        figures.append(fig)

    return figures

def plot_stability_comparison(case_data, save_path=None):
    """Plot stability comparison across different Reynolds numbers"""
    fig, axes = plt.subplots(3, 2, figsize=(14, 12))

    cases = ['Re~10', 'Re~100', 'Re~1000']

    for i, (case_name, data) in enumerate(case_data.items()):
        bgk_data = data['bgk']
        elbm_data = data['elbm']

        # BGK
        ax = axes[i, 0]
        ax.plot(bgk_data['x'], bgk_data['p'], 'b-', linewidth=2)
        ax.set_xlabel('x (lattice units)')
        ax.set_ylabel('Pressure')
        ax.set_title(f'{cases[i]} - BGK', fontweight='bold')
        ax.grid(True, alpha=0.3)

        # ELBM
        ax = axes[i, 1]
        ax.plot(elbm_data['x'], elbm_data['p'], 'r-', linewidth=2)
        ax.set_xlabel('x (lattice units)')
        ax.set_ylabel('Pressure')
        ax.set_title(f'{cases[i]} - ELBM', fontweight='bold')
        ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    return fig
