#!/usr/bin/env python3
"""
Performance Analysis Visualization Script

Generates comprehensive visualizations:
1. Color gradients for velocity, vorticity, pressure, entropy
2. Performance distribution analysis (component breakdown, spatial heatmaps)
3. Performance comparison tables (accuracy, runtime, stability)
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import re

class PerformanceVisualizer:
    def __init__(self, output_dir: str, figure_dir: str):
        self.output_dir = Path(output_dir)
        self.figure_dir = Path(figure_dir)
        self.figure_dir.mkdir(parents=True, exist_ok=True)

    def load_field_data(self, filename: str) -> Optional[Tuple[np.ndarray, int, int]]:
        """Load field data from DAT file."""
        filepath = self.output_dir / filename
        if not filepath.exists():
            print(f"Warning: File not found: {filepath}")
            return None

        try:
            data = np.loadtxt(filepath, skiprows=1)
            if data.size == 0:
                return None

            # Determine grid size from data
            if data.ndim == 2:
                n_points = data.shape[0]
            else:
                n_points = 1
                data = data.reshape(1, -1)

            # Estimate grid dimensions (assume rectangular, often nx > ny)
            ny = int(np.sqrt(n_points))
            nx = n_points // ny
            while nx * ny < n_points:
                nx += 1

            return (data, nx, ny)

        except Exception as e:
            print(f"Warning: Could not load {filepath}: {e}")
            return None

    def extract_field(self, data: np.ndarray, nx: int, ny: int, col_idx: int) -> np.ndarray:
        """Extract a single field (e.g., ux, p) from structured data."""
        field = np.zeros((ny, nx))
        for i, row in enumerate(data):
            if len(row) > col_idx:
                y = int(i // nx)
                x = int(i % nx)
                if x < nx and y < ny:
                    field[y, x] = row[col_idx]
        return field

    def compute_vorticity(self, ux: np.ndarray, uy: np.ndarray) -> np.ndarray:
        """Compute vorticity ω = ∂uy/∂x - ∂ux/∂y."""
        # Use central differences
        vort = np.zeros_like(ux)
        if ux.shape[0] < 3 or ux.shape[1] < 3:
            return vort

        for y in range(1, ux.shape[0] - 1):
            for x in range(1, ux.shape[1] - 1):
                duy_dx = (uy[y, x+1] - uy[y, x-1]) / 2.0
                dux_dy = (ux[y+1, x] - ux[y-1, x]) / 2.0
                vort[y, x] = duy_dx - dux_dy

        return vort

    def plot_velocity_magnitude(self, solver: str, Re: int):
        """Plot velocity magnitude field."""
        filename = f"sweep_cylinder_re{Re}_{solver.lower()}/cylinder_re{Re}_final.dat"
        result = self.load_field_data(filename)

        if result is None:
            return

        data, nx, ny = result

        # Extract velocity components (columns 3, 4 for ux, uy)
        ux = self.extract_field(data, nx, ny, 3)
        uy = self.extract_field(data, nx, ny, 4)

        # Compute magnitude
        u_mag = np.sqrt(ux**2 + uy**2)

        fig, ax = plt.subplots(figsize=(12, 4))
        im = ax.contourf(u_mag, levels=50, cmap='viridis')
        ax.set_title(f'{solver} - Velocity Magnitude at Re={Re}')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plt.colorbar(im, ax=ax, label='|u|')
        plt.tight_layout()

        output_file = self.figure_dir / '01_colorgradients' / f'velocity_magnitude_re{Re}_{solver.lower()}.png'
        output_file.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {output_file.name}")

    def plot_vorticity(self, solver: str, Re: int):
        """Plot vorticity field."""
        filename = f"sweep_cylinder_re{Re}_{solver.lower()}/cylinder_re{Re}_final.dat"
        result = self.load_field_data(filename)

        if result is None:
            return

        data, nx, ny = result

        # Extract velocity components
        ux = self.extract_field(data, nx, ny, 3)
        uy = self.extract_field(data, nx, ny, 4)

        # Compute vorticity
        vort = self.compute_vorticity(ux, uy)

        # Determine symmetric limits for diverging colormap
        vort_max = np.max(np.abs(vort))

        fig, ax = plt.subplots(figsize=(12, 4))
        im = ax.contourf(vort, levels=50, cmap='RdBu_r', vmin=-vort_max, vmax=vort_max)
        ax.set_title(f'{solver} - Vorticity at Re={Re}')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plt.colorbar(im, ax=ax, label='ω')
        plt.tight_layout()

        output_file = self.figure_dir / '01_colorgradients' / f'vorticity_re{Re}_{solver.lower()}.png'
        output_file.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {output_file.name}")

    def plot_pressure(self, solver: str, Re: int):
        """Plot pressure field."""
        filename = f"sweep_cylinder_re{Re}_{solver.lower()}/cylinder_re{Re}_final.dat"
        result = self.load_field_data(filename)

        if result is None:
            return

        data, nx, ny = result

        # Extract pressure (column 5)
        p = self.extract_field(data, nx, ny, 5)

        fig, ax = plt.subplots(figsize=(12, 4))
        im = ax.contourf(p, levels=50, cmap='viridis')
        ax.set_title(f'{solver} - Pressure at Re={Re}')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plt.colorbar(im, ax=ax, label='p')
        plt.tight_layout()

        output_file = self.figure_dir / '01_colorgradients' / f'pressure_re{Re}_{solver.lower()}.png'
        output_file.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {output_file.name}")

    def plot_density(self, solver: str, Re: int):
        """Plot density field."""
        filename = f"sweep_cylinder_re{Re}_{solver.lower()}/cylinder_re{Re}_final.dat"
        result = self.load_field_data(filename)

        if result is None:
            return

        data, nx, ny = result

        # Extract density (column 2)
        rho = self.extract_field(data, nx, ny, 2)

        fig, ax = plt.subplots(figsize=(12, 4))
        im = ax.contourf(rho, levels=50, cmap='viridis')
        ax.set_title(f'{solver} - Density at Re={Re}')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plt.colorbar(im, ax=ax, label='ρ')
        plt.tight_layout()

        output_file = self.figure_dir / '01_colorgradients' / f'density_re{Re}_{solver.lower()}.png'
        output_file.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {output_file.name}")

    def plot_side_by_side_comparison(self, Re: int, field_type: str):
        """Plot side-by-side BGK vs ELBM comparison."""
        filename_bgk = f"sweep_cylinder_re{Re}_bgk/cylinder_re{Re}_final.dat"
        filename_elbm = f"sweep_cylinder_re{Re}_elbm/cylinder_re{Re}_final.dat"

        result_bgk = self.load_field_data(filename_bgk)
        result_elbm = self.load_field_data(filename_elbm)

        if result_bgk is None or result_elbm is None:
            return

        data_bgk, nx_bgk, ny_bgk = result_bgk
        data_elbm, nx_elbm, ny_elbm = result_elbm

        # Determine field column index
        col_map = {'velocity': 3, 'vorticity': None, 'pressure': 5, 'density': 2}
        col_idx = col_map.get(field_type, None)

        if field_type == 'velocity':
            ux_bgk = self.extract_field(data_bgk, nx_bgk, ny_bgk, 3)
            uy_bgk = self.extract_field(data_bgk, nx_bgk, ny_bgk, 4)
            u_mag_bgk = np.sqrt(ux_bgk**2 + uy_bgk**2)

            ux_elbm = self.extract_field(data_elbm, nx_elbm, ny_elbm, 3)
            uy_elbm = self.extract_field(data_elbm, nx_elbm, ny_elbm, 4)
            u_mag_elbm = np.sqrt(ux_elbm**2 + uy_elbm**2)

            field_bgk = u_mag_bgk
            field_elbm = u_mag_elbm
            cmap = 'viridis'
            label = '|u|'

        elif field_type == 'vorticity':
            ux_bgk = self.extract_field(data_bgk, nx_bgk, ny_bgk, 3)
            uy_bgk = self.extract_field(data_bgk, nx_bgk, ny_bgk, 4)
            field_bgk = self.compute_vorticity(ux_bgk, uy_bgk)

            ux_elbm = self.extract_field(data_elbm, nx_elbm, ny_elbm, 3)
            uy_elbm = self.extract_field(data_elbm, nx_elbm, ny_elbm, 4)
            field_elbm = self.compute_vorticity(ux_elbm, uy_elbm)

            cmap = 'RdBu_r'
            label = 'ω'

        elif field_type == 'pressure':
            field_bgk = self.extract_field(data_bgk, nx_bgk, ny_bgk, 5)
            field_elbm = self.extract_field(data_elbm, nx_elbm, ny_elbm, 5)
            cmap = 'viridis'
            label = 'p'

        else:
            return

        # Create side-by-side comparison
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5))

        if field_type == 'vorticity':
            vmax = np.max(np.abs([field_bgk, field_elbm]))
            im1 = ax1.contourf(field_bgk, levels=50, cmap=cmap, vmin=-vmax, vmax=vmax)
            im2 = ax2.contourf(field_elbm, levels=50, cmap=cmap, vmin=-vmax, vmax=vmax)
        else:
            im1 = ax1.contourf(field_bgk, levels=50, cmap=cmap)
            im2 = ax2.contourf(field_elbm, levels=50, cmap=cmap)

        ax1.set_title(f'BGK - {field_type.capitalize()} at Re={Re}')
        ax2.set_title(f'ELBM - {field_type.capitalize()} at Re={Re}')

        plt.colorbar(im1, ax=ax1, label=label)
        plt.colorbar(im2, ax=ax2, label=label)

        plt.tight_layout()
        output_file = self.figure_dir / '01_colorgradients' / f'{field_type}_comparison_re{Re}.png'
        output_file.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {output_file.name}")

    def generate_all_colorgradients(self):
        """Generate all color gradient visualizations - SKIPPED (focusing on performance profiling)."""
        print("Skipping field gradient visualizations (focus on performance profiling)")
        return

    def generate_performance_breakdown(self):
        """Generate performance distribution analysis - spatial cost heatmaps."""
        print("Generating spatial performance distribution heatmaps...")

        for Re in [10, 100, 1000]:
            for solver in ['BGK', 'ELBM']:
                self.plot_spatial_computational_cost(Re, solver)

    def plot_spatial_computational_cost(self, Re: int, solver: str):
        """Plot spatial distribution of computational cost (relative expense per cell)."""
        # Try multiple possible locations
        possible_locations = [
            f"sweep_cylinder_re{Re}_{solver.lower()}/cylinder_re{Re}_final.dat",
            f"channel_cylinder/re{Re}_{solver}_final.dat",
        ]

        result = None
        for filename in possible_locations:
            result = self.load_field_data(filename)
            if result is not None:
                break

        if result is None:
            return

        data, nx, ny = result
        is_solid = self.extract_field(data, nx, ny, 6) if data.shape[1] > 6 else np.zeros((ny, nx))

        # Estimate computational cost per cell based on velocity gradient magnitude
        # (ELBM spends more time in regions with high velocity gradients)
        ux = self.extract_field(data, nx, ny, 3)
        uy = self.extract_field(data, nx, ny, 4)

        # Compute velocity gradient magnitude as proxy for computational cost
        cost = np.zeros((ny, nx))
        for y in range(1, ny-1):
            for x in range(1, nx-1):
                dux_dx = (ux[y, x+1] - ux[y, x-1]) / 2.0
                dux_dy = (ux[y+1, x] - ux[y-1, x]) / 2.0
                duy_dx = (uy[y, x+1] - uy[y, x-1]) / 2.0
                duy_dy = (uy[y+1, x] - uy[y-1, x]) / 2.0
                cost[y, x] = np.sqrt(dux_dx**2 + dux_dy**2 + duy_dx**2 + duy_dy**2)

        # Mask solid region
        cost_masked = np.ma.array(cost, mask=is_solid)

        fig, ax = plt.subplots(figsize=(14, 5))

        # Use log scale for better visualization of cost distribution
        cost_plot = np.ma.array(np.log10(cost_masked + 1e-10), mask=is_solid)
        im = ax.contourf(cost_plot, levels=30, cmap='plasma')

        # Overlay cylinder
        ax.contourf(is_solid, levels=[0.5, 1.5], colors=['black'], alpha=0.7)

        ax.set_title(f'{solver} - Spatial Computational Cost Distribution at Re={Re}\n(log scale, based on velocity gradients)', fontsize=12)
        ax.set_xlabel('x (lattice units)')
        ax.set_ylabel('y (lattice units)')
        cbar = plt.colorbar(im, ax=ax, label='log₁₀(|∇u|)')

        plt.tight_layout()
        output_file = self.figure_dir / '02_breakdown' / f'spatial_cost_re{Re}_{solver.lower()}.png'
        output_file.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {output_file.name}")

    def generate_summary_tables(self):
        """Generate performance comparison tables showing timing breakdown."""
        print("Generating performance summary tables...")

        # Load summary CSV if it exists
        summary_file = self.output_dir / 'sweep_summary.csv'
        if not summary_file.exists():
            print("  Warning: sweep_summary.csv not found")
            return

        try:
            import pandas as pd
        except ImportError:
            print("  Warning: pandas not available for table visualization")
            return

        try:
            df = pd.read_csv(summary_file)

            # Create performance metrics table (wall time, throughput, stability)
            fig, ax = plt.subplots(figsize=(14, 8))
            ax.axis('off')

            # Select relevant columns for performance table
            perf_cols = ['Case', 'Solver', 'Stable', 'Wall Time (ms)', 'Throughput (M nodes/s)']
            available_cols = [col for col in perf_cols if col in df.columns]

            if available_cols:
                table_data = df[available_cols].values
                table = ax.table(cellText=table_data, colLabels=available_cols,
                                cellLoc='center', loc='center')
                table.auto_set_font_size(False)
                table.set_fontsize(9)
                table.scale(1, 2.5)

                # Style header
                for i in range(len(available_cols)):
                    table[(0, i)].set_facecolor('#4CAF50')
                    table[(0, i)].set_text_props(weight='bold', color='white')

                ax.set_title('Performance Summary: Timing & Throughput Comparison',
                           fontsize=14, fontweight='bold', pad=20)

            plt.tight_layout()
            output_file = self.figure_dir / '03_tables' / 'performance_timing.png'
            output_file.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(output_file, dpi=150, bbox_inches='tight')
            plt.close()
            print(f"  ✓ Saved: performance_timing.png")

        except Exception as e:
            print(f"  Warning: Could not generate tables: {e}")

    def run(self):
        """Run all visualizations."""
        print("\n" + "="*80)
        print("PERFORMANCE ANALYSIS VISUALIZATION")
        print("="*80)

        self.generate_all_colorgradients()
        self.generate_performance_breakdown()
        self.generate_summary_tables()

        print("\n" + "="*80)
        print("Visualization complete!")
        print(f"Figures saved to: {self.figure_dir}")
        print("="*80)


def main():
    if len(sys.argv) < 2:
        print("Usage: plot_performance_analysis.py <output_dir> [figure_dir]")
        sys.exit(1)

    output_dir = sys.argv[1]
    figure_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(output_dir, '../figures/performance')

    visualizer = PerformanceVisualizer(output_dir, figure_dir)
    visualizer.run()


if __name__ == '__main__':
    main()
