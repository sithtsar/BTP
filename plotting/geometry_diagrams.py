#!/usr/bin/env python3
"""
Geometry Diagram Generator for ELBM Test Cases

Creates publication-quality domain schematics showing:
- Domain dimensions
- Boundary condition types
- Key parameters (Re, ν, U, etc.)
- Timestep information
- Coordinate system
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle
import numpy as np
from pathlib import Path

class GeometryDiagram:
    """Base class for creating domain geometry diagrams"""

    def __init__(self, figsize=(10, 6)):
        self.fig, self.ax = plt.subplots(1, 1, figsize=figsize)
        self.ax.set_aspect('equal')
        self.ax.axis('off')

    def add_boundary(self, x0, y0, x1, y1, bc_type, label=None, color='blue'):
        """Add boundary condition visualization"""
        bc_styles = {
            'velocity': {'linestyle': '-', 'linewidth': 4, 'color': 'red'},
            'pressure': {'linestyle': '-', 'linewidth': 4, 'color': 'blue'},
            'periodic': {'linestyle': '--', 'linewidth': 3, 'color': 'green'},
            'bounce_back': {'linestyle': '-', 'linewidth': 5, 'color': 'black'}
        }

        style = bc_styles.get(bc_type, {'linestyle': '-', 'linewidth': 2, 'color': color})
        self.ax.plot([x0, x1], [y0, y1], **style, zorder=10)

        if label:
            mid_x, mid_y = (x0 + x1) / 2, (y0 + y1) / 2

            if abs(y1 - y0) < 0.1:  # Horizontal boundary
                offset_y = 0.05 if y0 < 0.5 else -0.05
                va = 'bottom' if y0 < 0.5 else 'top'
                self.ax.text(mid_x, mid_y + offset_y, label,
                           ha='center', va=va, fontsize=11, fontweight='bold',
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
            else:  # Vertical boundary
                offset_x = -0.08 if x0 < 0.5 else 0.08
                ha = 'right' if x0 < 0.5 else 'left'
                self.ax.text(mid_x + offset_x, mid_y, label,
                           ha=ha, va='center', fontsize=11, fontweight='bold',
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    def add_domain_box(self, nx, ny):
        """Add domain rectangle with dimensions"""
        rect = Rectangle((0, 0), 1, 1, fill=False, edgecolor='gray',
                        linewidth=2, linestyle='--', zorder=1)
        self.ax.add_patch(rect)

        # Dimension annotations
        self.ax.annotate('', xy=(1.05, 0), xytext=(1.05, 1),
                        arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
        self.ax.text(1.12, 0.5, f'ny = {ny}', ha='left', va='center',
                    fontsize=10, rotation=90)

        self.ax.annotate('', xy=(0, -0.05), xytext=(1, -0.05),
                        arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
        self.ax.text(0.5, -0.12, f'nx = {nx}', ha='center', va='top', fontsize=10)

    def add_coordinate_system(self, x=0.05, y=0.05, length=0.08):
        """Add coordinate system arrows"""
        self.ax.arrow(x, y, length, 0, head_width=0.02, head_length=0.02,
                     fc='black', ec='black', linewidth=2, zorder=20)
        self.ax.text(x + length + 0.02, y, 'x', fontsize=12, fontweight='bold')

        self.ax.arrow(x, y, 0, length, head_width=0.02, head_length=0.02,
                     fc='black', ec='black', linewidth=2, zorder=20)
        self.ax.text(x, y + length + 0.02, 'y', fontsize=12, fontweight='bold')

    def add_parameter_box(self, params, x=0.02, y=0.97):
        """Add parameter information box"""
        param_text = '\n'.join([f'{k}: {v}' for k, v in params.items()])
        self.ax.text(x, y, param_text, transform=self.ax.transAxes,
                    ha='left', va='top', fontsize=10, fontfamily='monospace',
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow',
                             alpha=0.9, edgecolor='black', linewidth=1.5))

    def add_title(self, title):
        """Add centered title"""
        self.ax.text(0.5, 1.15, title, ha='center', va='bottom',
                    fontsize=14, fontweight='bold', transform=self.ax.transAxes)

    def save(self, filepath):
        """Save diagram to file"""
        self.ax.set_xlim(-0.2, 1.2)
        self.ax.set_ylim(-0.2, 1.3)
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved geometry diagram: {filepath}")


def create_couette_diagram():
    """Couette flow: Moving top wall, stationary bottom"""
    diag = GeometryDiagram()

    diag.add_domain_box(nx=64, ny=64)
    diag.add_boundary(0, 1, 1, 1, 'velocity', label='Moving Wall (u=U₀)')
    diag.add_boundary(0, 0, 1, 0, 'bounce_back', label='Stationary Wall (u=0)')
    diag.add_boundary(0, 0, 0, 1, 'periodic', label='Periodic')
    diag.add_boundary(1, 0, 1, 1, 'periodic', label='Periodic')

    # Velocity profile schematic
    y_profile = np.linspace(0, 1, 50)
    u_profile = 0.15 * y_profile / 1.0  # Linear
    u_profile += 0.5
    diag.ax.plot(u_profile, y_profile, 'r--', linewidth=2, alpha=0.7, label='u(y) profile')

    diag.add_coordinate_system()

    params = {
        'Grid': '64 × 64',
        'U₀': '0.1',
        'ν': '0.1',
        'τ': '0.8',
        'Steps': '20,000',
        'Re': '~6.4'
    }
    diag.add_parameter_box(params)
    diag.add_title('Couette Flow: Analytical Validation')

    fig_dir = Path('figures/geometry_diagrams')
    fig_dir.mkdir(parents=True, exist_ok=True)
    diag.save(fig_dir / 'couette_geometry.png')


def create_poiseuille_diagram():
    """Poiseuille flow: Pressure-driven channel"""
    diag = GeometryDiagram()

    diag.add_domain_box(nx=100, ny=50)
    diag.add_boundary(0, 1, 1, 1, 'bounce_back', label='No-slip Wall')
    diag.add_boundary(0, 0, 1, 0, 'bounce_back', label='No-slip Wall')
    diag.add_boundary(0, 0, 0, 1, 'periodic', label='Periodic')
    diag.add_boundary(1, 0, 1, 1, 'periodic', label='Periodic')

    # Body force arrow
    arrow = FancyArrowPatch((0.2, 0.5), (0.8, 0.5),
                           arrowstyle='->', mutation_scale=30,
                           linewidth=3, color='orange', zorder=15)
    diag.ax.add_patch(arrow)
    diag.ax.text(0.5, 0.55, 'Body Force (Fₓ)', ha='center', fontsize=11,
                fontweight='bold', color='orange')

    # Parabolic velocity profile
    y_profile = np.linspace(0, 1, 50)
    u_profile = 0.2 * 4 * y_profile * (1 - y_profile)  # Parabolic
    u_profile += 0.5
    diag.ax.plot(u_profile, y_profile, 'r--', linewidth=2, alpha=0.7)

    diag.add_coordinate_system()

    params = {
        'Grid': '100 × 50',
        'Fₓ': '1×10⁻⁵',
        'ν': '0.1',
        'τ': '0.8',
        'Steps': '20,000',
        'Max u': '~0.006'
    }
    diag.add_parameter_box(params)
    diag.add_title('Poiseuille Flow: Pressure-Driven Channel')

    fig_dir = Path('figures/geometry_diagrams')
    fig_dir.mkdir(parents=True, exist_ok=True)
    diag.save(fig_dir / 'poiseuille_geometry.png')


def create_taylor_green_diagram():
    """Taylor-Green vortex: Decaying periodic vortices"""
    diag = GeometryDiagram(figsize=(10, 10))

    diag.add_domain_box(nx=64, ny=64)
    diag.add_boundary(0, 1, 1, 1, 'periodic', label='Periodic')
    diag.add_boundary(0, 0, 1, 0, 'periodic', label='Periodic')
    diag.add_boundary(0, 0, 0, 1, 'periodic', label='Periodic')
    diag.add_boundary(1, 0, 1, 1, 'periodic', label='Periodic')

    # Vortex visualization
    x = np.linspace(0, 1, 20)
    y = np.linspace(0, 1, 20)
    X, Y = np.meshgrid(x, y)

    U = -0.05 * np.cos(2 * np.pi * X) * np.sin(2 * np.pi * Y)
    V = 0.05 * np.sin(2 * np.pi * X) * np.cos(2 * np.pi * Y)

    diag.ax.streamplot(X, Y, U, V, density=1.5, color='blue',
                      linewidth=1, arrowsize=1.2, zorder=5)

    # Mark vortex centers
    centers = [(0.25, 0.25), (0.75, 0.25), (0.25, 0.75), (0.75, 0.75)]
    for (cx, cy) in centers:
        diag.ax.plot(cx, cy, 'ro', markersize=8, zorder=10)

    diag.add_coordinate_system()

    params = {
        'Grid': '64 × 64',
        'U₀': '0.05',
        'ν': '0.1',
        'τ': '0.8',
        'Steps': '10,000',
        'k': '2π/L'
    }
    diag.add_parameter_box(params, y=0.05)

    # Add decay equation
    diag.ax.text(0.98, 0.97, 'E(t) = E₀ exp(-2νk²t)',
                transform=diag.ax.transAxes,
                ha='right', va='top', fontsize=11, fontfamily='serif',
                bbox=dict(boxstyle='round', facecolor='lightcyan', alpha=0.9))

    diag.add_title('Taylor-Green Vortex: Viscous Decay Validation')

    fig_dir = Path('figures/geometry_diagrams')
    fig_dir.mkdir(parents=True, exist_ok=True)
    diag.save(fig_dir / 'taylor_green_geometry.png')


def create_cylinder_diagram(Re):
    """Channel flow with cylinder"""
    diag = GeometryDiagram()

    diag.add_domain_box(nx=400, ny=120)
    diag.add_boundary(0, 0, 0, 1, 'velocity', label=f'Inlet (U∞)')
    diag.add_boundary(1, 0, 1, 1, 'pressure', label='Outlet (P=P₀)')
    diag.add_boundary(0, 1, 1, 1, 'bounce_back', label='Wall')
    diag.add_boundary(0, 0, 1, 0, 'bounce_back', label='Wall')

    # Cylinder
    cylinder_x = 0.33
    cylinder_y = 0.5
    cylinder_r = 20.0 / 400.0

    circle = Circle((cylinder_x, cylinder_y), cylinder_r,
                   facecolor='gray', edgecolor='black', linewidth=2, zorder=10)
    diag.ax.add_patch(circle)
    diag.ax.text(cylinder_x, cylinder_y, 'Cylinder\nD=20',
                ha='center', va='center', fontsize=9, fontweight='bold')

    # Flow direction arrow
    arrow = FancyArrowPatch((0.05, 0.5), (0.25, 0.5),
                           arrowstyle='->', mutation_scale=25,
                           linewidth=2.5, color='red', zorder=8)
    diag.ax.add_patch(arrow)

    diag.add_coordinate_system()

    # Compute velocity from Re
    U_inf = (Re * 0.02) / 20.0

    params = {
        'Grid': '400 × 120',
        'Re': str(Re),
        'U∞': f'{U_inf:.4f}',
        'D': '20',
        'ν': '0.02 (fixed)',
        'Steps': '10,000',
        'Scheme': 'Immersed Boundary'
    }
    diag.add_parameter_box(params)
    diag.add_title(f'Cylinder Flow Benchmark: Re = {Re}')

    fig_dir = Path('figures/geometry_diagrams')
    fig_dir.mkdir(parents=True, exist_ok=True)
    diag.save(fig_dir / f'cylinder_re{Re}_geometry.png')


def create_twophase_diagram():
    """Two-phase droplet relaxation"""
    diag = GeometryDiagram(figsize=(10, 10))

    diag.add_domain_box(nx=128, ny=128)
    diag.add_boundary(0, 1, 1, 1, 'periodic', label='Periodic')
    diag.add_boundary(0, 0, 1, 0, 'periodic', label='Periodic')
    diag.add_boundary(0, 0, 0, 1, 'periodic', label='Periodic')
    diag.add_boundary(1, 0, 1, 1, 'periodic', label='Periodic')

    # Droplet (liquid phase)
    droplet = Circle((0.5, 0.5), 30.0/128.0,
                    facecolor='blue', edgecolor='darkblue',
                    linewidth=2, alpha=0.6, zorder=10)
    diag.ax.add_patch(droplet)
    diag.ax.text(0.5, 0.5, 'Liquid\nρ=1.2', ha='center', va='center',
                fontsize=11, fontweight='bold', color='white')

    # Gas phase label
    diag.ax.text(0.15, 0.15, 'Gas\nρ=0.8', ha='center', va='center',
                fontsize=11, fontweight='bold',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7))

    # Interface annotation
    diag.ax.annotate('', xy=(0.5 + 30.0/128.0, 0.5), xytext=(0.75, 0.5),
                    arrowprops=dict(arrowstyle='->', lw=1.5, color='red'))
    diag.ax.text(0.78, 0.5, 'Interface\nR₀=30', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

    diag.add_coordinate_system()

    params = {
        'Grid': '128 × 128',
        'ρ_liquid': '1.2',
        'ρ_gas': '0.8',
        'Ratio': '1.5:1',
        'R₀': '30',
        'ν': '0.1',
        'Steps': '5,000',
        'Validation': 'H-theorem'
    }
    diag.add_parameter_box(params)
    diag.add_title('Two-Phase Droplet: H-Theorem Validation')

    fig_dir = Path('figures/geometry_diagrams')
    fig_dir.mkdir(parents=True, exist_ok=True)
    diag.save(fig_dir / 'twophase_geometry.png')


def main():
    """Generate all geometry diagrams"""
    print("Generating geometry diagrams for all test cases...")

    create_couette_diagram()
    create_poiseuille_diagram()
    create_taylor_green_diagram()

    for re in [10, 50, 100]:
        create_cylinder_diagram(re)

    create_twophase_diagram()

    print("\nAll geometry diagrams generated successfully!")
    print("Location: figures/geometry_diagrams/")

if __name__ == '__main__':
    main()
