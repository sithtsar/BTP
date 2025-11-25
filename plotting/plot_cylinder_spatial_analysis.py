#!/usr/bin/env python3
"""
Cylinder Flow Spatial Analysis Visualization

Generates spatial distribution visualizations:
1. Vorticity heatmaps around cylinder
2. Alpha-solver iteration count distribution
3. Performance concentration analysis
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path
from typing import Optional, Tuple

class CylinderSpatialAnalyzer:
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

            # Determine grid size
            n_points = data.shape[0] if data.ndim == 2 else 1
            ny = int(np.sqrt(n_points))
            nx = n_points // ny
            while nx * ny < n_points:
                nx += 1

            return (data, nx, ny)

        except Exception as e:
            print(f"Warning: Could not load {filepath}: {e}")
            return None

    def extract_field(self, data: np.ndarray, nx: int, ny: int, col_idx: int) -> np.ndarray:
        """Extract a single field column."""
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
        vort = np.zeros_like(ux)
        if ux.shape[0] < 3 or ux.shape[1] < 3:
            return vort

        for y in range(1, ux.shape[0] - 1):
            for x in range(1, ux.shape[1] - 1):
                duy_dx = (uy[y, x+1] - uy[y, x-1]) / 2.0
                dux_dy = (ux[y+1, x] - ux[y-1, x]) / 2.0
                vort[y, x] = duy_dx - dux_dy

        return vort

    def get_cylinder_mask(self, data: np.ndarray, nx: int, ny: int) -> np.ndarray:
        """Extract cylinder geometry from is_solid column."""
        is_solid = np.zeros((ny, nx))
        if data.shape[1] <= 6:
            return is_solid

        for i, row in enumerate(data):
            y = int(i // nx)
            x = int(i % nx)
            if x < nx and y < ny and len(row) > 6:
                is_solid[y, x] = row[6]

        return is_solid

    def plot_vorticity_heatmap(self, Re: int, solver: str):
        """Plot vorticity heatmap with cylinder overlay."""
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

        # Extract velocity and cylinder mask
        ux = self.extract_field(data, nx, ny, 3)
        uy = self.extract_field(data, nx, ny, 4)
        is_solid = self.get_cylinder_mask(data, nx, ny)

        # Compute vorticity
        vort = self.compute_vorticity(ux, uy)

        # Mask solid region
        vort_masked = np.ma.array(vort, mask=is_solid)

        fig, ax = plt.subplots(figsize=(14, 5))

        # Plot vorticity with symmetric colormap
        vort_max = np.max(np.abs(vort_masked))
        im = ax.contourf(vort_masked, levels=50, cmap='RdBu_r', vmin=-vort_max, vmax=vort_max)

        # Overlay cylinder as black region
        solid_contour = ax.contourf(is_solid, levels=[0.5, 1.5], colors=['black'], alpha=0.7)

        ax.set_title(f'{solver} - Vorticity Distribution at Re={Re} (Cylinder in Black)', fontsize=12)
        ax.set_xlabel('x (lattice units)')
        ax.set_ylabel('y (lattice units)')
        plt.colorbar(im, ax=ax, label='Vorticity ω')

        plt.tight_layout()
        output_file = self.figure_dir / f'vorticity_spatial_re{Re}_{solver.lower()}.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {output_file.name}")

    def plot_velocity_streamlines(self, Re: int, solver: str):
        """Plot velocity field with streamlines."""
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

        # Extract velocity
        ux = self.extract_field(data, nx, ny, 3)
        uy = self.extract_field(data, nx, ny, 4)
        is_solid = self.get_cylinder_mask(data, nx, ny)

        fig, ax = plt.subplots(figsize=(14, 5))

        # Plot velocity magnitude as background
        u_mag = np.sqrt(ux**2 + uy**2)
        u_mag_masked = np.ma.array(u_mag, mask=is_solid)

        im = ax.contourf(u_mag_masked, levels=30, cmap='viridis')

        # Create coordinate grids for streamlines
        y_coords, x_coords = np.meshgrid(np.arange(ny), np.arange(nx), indexing='ij')

        # Plot streamlines (subsample for clarity)
        skip = max(ny // 10, 2)
        ax.quiver(x_coords[::skip, ::skip], y_coords[::skip, ::skip],
                 ux[::skip, ::skip], uy[::skip, ::skip],
                 alpha=0.5, scale=50)

        # Overlay cylinder
        ax.contourf(is_solid, levels=[0.5, 1.5], colors=['black'], alpha=0.7)

        ax.set_title(f'{solver} - Velocity Field at Re={Re}', fontsize=12)
        ax.set_xlabel('x (lattice units)')
        ax.set_ylabel('y (lattice units)')
        plt.colorbar(im, ax=ax, label='|u|')

        plt.tight_layout()
        output_file = self.figure_dir / f'velocity_field_re{Re}_{solver.lower()}.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {output_file.name}")

    def plot_comparison_heatmap(self, Re: int, field_type: str = 'vorticity'):
        """Plot side-by-side comparison of spatial distributions."""
        # Try multiple possible locations
        bgk_locations = [
            f"sweep_cylinder_re{Re}_bgk/cylinder_re{Re}_final.dat",
            f"channel_cylinder/re{Re}_BGK_final.dat",
        ]
        elbm_locations = [
            f"sweep_cylinder_re{Re}_elbm/cylinder_re{Re}_final.dat",
            f"channel_cylinder/re{Re}_ELBM_final.dat",
        ]

        result_bgk = None
        for filename in bgk_locations:
            result_bgk = self.load_field_data(filename)
            if result_bgk is not None:
                break

        result_elbm = None
        for filename in elbm_locations:
            result_elbm = self.load_field_data(filename)
            if result_elbm is not None:
                break

        if result_bgk is None or result_elbm is None:
            return

        data_bgk, nx_bgk, ny_bgk = result_bgk
        data_elbm, nx_elbm, ny_elbm = result_elbm

        if field_type == 'vorticity':
            ux_bgk = self.extract_field(data_bgk, nx_bgk, ny_bgk, 3)
            uy_bgk = self.extract_field(data_bgk, nx_bgk, ny_bgk, 4)
            vort_bgk = self.compute_vorticity(ux_bgk, uy_bgk)

            ux_elbm = self.extract_field(data_elbm, nx_elbm, ny_elbm, 3)
            uy_elbm = self.extract_field(data_elbm, nx_elbm, ny_elbm, 4)
            vort_elbm = self.compute_vorticity(ux_elbm, uy_elbm)

            is_solid_bgk = self.get_cylinder_mask(data_bgk, nx_bgk, ny_bgk)
            is_solid_elbm = self.get_cylinder_mask(data_elbm, nx_elbm, ny_elbm)

            field_bgk = np.ma.array(vort_bgk, mask=is_solid_bgk)
            field_elbm = np.ma.array(vort_elbm, mask=is_solid_elbm)

            title_suffix = 'Vorticity'
            cmap = 'RdBu_r'
            label = 'ω'

        else:
            return

        # Create comparison
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5))

        vmax = np.max(np.abs([field_bgk, field_elbm]))
        im1 = ax1.contourf(field_bgk, levels=50, cmap=cmap, vmin=-vmax, vmax=vmax)
        im2 = ax2.contourf(field_elbm, levels=50, cmap=cmap, vmin=-vmax, vmax=vmax)

        # Overlay cylinders
        ax1.contourf(is_solid_bgk, levels=[0.5, 1.5], colors=['black'], alpha=0.7)
        ax2.contourf(is_solid_elbm, levels=[0.5, 1.5], colors=['black'], alpha=0.7)

        ax1.set_title(f'BGK - {title_suffix} at Re={Re}')
        ax2.set_title(f'ELBM - {title_suffix} at Re={Re}')

        plt.colorbar(im1, ax=ax1, label=label)
        plt.colorbar(im2, ax=ax2, label=label)

        plt.tight_layout()
        output_file = self.figure_dir / f'{field_type}_comparison_spatial_re{Re}.png'
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  ✓ Saved: {output_file.name}")

    def run(self):
        """Generate all spatial analysis visualizations."""
        print("\n" + "="*80)
        print("CYLINDER SPATIAL ANALYSIS")
        print("="*80)

        print("\nGenerating vorticity heatmaps...")
        for Re in [10, 100, 1000]:
            for solver in ['BGK', 'ELBM']:
                self.plot_vorticity_heatmap(Re, solver)

        print("\nGenerating velocity field visualizations...")
        for Re in [10, 100, 1000]:
            for solver in ['BGK', 'ELBM']:
                self.plot_velocity_streamlines(Re, solver)

        print("\nGenerating spatial comparison plots...")
        for Re in [10, 100, 1000]:
            self.plot_comparison_heatmap(Re, 'vorticity')

        print("\n" + "="*80)
        print("Spatial analysis complete!")
        print(f"Figures saved to: {self.figure_dir}")
        print("="*80)


def main():
    if len(sys.argv) < 2:
        print("Usage: plot_cylinder_spatial_analysis.py <output_dir> [figure_dir]")
        sys.exit(1)

    output_dir = sys.argv[1]
    figure_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(output_dir, '../figures/performance')

    analyzer = CylinderSpatialAnalyzer(output_dir, figure_dir)
    analyzer.run()


if __name__ == '__main__':
    main()
