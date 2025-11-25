#!/usr/bin/env python3
"""
Hagen-Poiseuille Flow: Pressure-driven flow in circular pipe
Analytical solution validation
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def hagen_poiseuille_analytical(r, R, dp_dx, rho, viscosity):
    """
    Hagen-Poiseuille analytical solution
    u(r) = -(dp/dx) * (R² - r²) / (4μρ)
    """
    mu = rho * viscosity
    u = -dp_dx * (R**2 - r**2) / (4.0 * mu)
    u[r > R] = 0  # Zero outside pipe
    return u

def main():
    # Parameters
    R = 30.0  # Pipe radius
    dp_dx = -0.001  # Pressure gradient
    rho = 1.0
    viscosity = 0.1

    print(f"Hagen-Poiseuille Flow Analysis")
    print(f"=" * 60)
    print(f"Pipe radius: {R}")
    print(f"Pressure gradient: {dp_dx}")
    print(f"Density: {rho}")
    print(f"Viscosity: {viscosity}")
    print(f"=" * 60)

    # Radial coordinate
    r = np.linspace(0, R + 5, 200)

    # Analytical solution
    u_analytical = hagen_poiseuille_analytical(r, R, dp_dx, rho, viscosity)

    # Maximum velocity (at r=0)
    u_max = hagen_poiseuille_analytical(np.array([0]), R, dp_dx, rho, viscosity)[0]
    print(f"\nMaximum velocity (centerline): {u_max:.6f}")

    # Volume flow rate
    Q_analytical = np.pi * R**4 * abs(dp_dx) / (8.0 * rho * viscosity)
    print(f"Volume flow rate: {Q_analytical:.6f}")

    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Radial profile
    ax1.plot(r, u_analytical, 'b-', linewidth=2, label='Analytical')
    ax1.axvline(x=R, color='r', linestyle='--', linewidth=2, alpha=0.5, label=f'Wall (r={R})')
    ax1.set_xlabel('Radial distance r')
    ax1.set_ylabel('Velocity u(r)')
    ax1.set_title('Velocity Profile (Parabolic)')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # 2D cross-section visualization
    theta = np.linspace(0, 2*np.pi, 100)
    x_circle = R * np.cos(theta)
    y_circle = R * np.sin(theta)

    # Create velocity field
    x_grid = np.linspace(-R-5, R+5, 100)
    y_grid = np.linspace(-R-5, R+5, 100)
    X, Y = np.meshgrid(x_grid, y_grid)
    R_grid = np.sqrt(X**2 + Y**2)
    U = hagen_poiseuille_analytical(R_grid.flatten(), R, dp_dx, rho, viscosity).reshape(X.shape)

    contour = ax2.contourf(X, Y, U, levels=20, cmap='viridis')
    ax2.plot(x_circle, y_circle, 'r-', linewidth=2, label='Pipe wall')
    ax2.set_xlabel('x position')
    ax2.set_ylabel('y position')
    ax2.set_title('2D Cross-section (Velocity Contours)')
    ax2.set_aspect('equal')
    plt.colorbar(contour, ax=ax2, label='Velocity')
    ax2.legend()

    plt.tight_layout()

    # Save
    output_dir = Path('figures/extended_analytical')
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / 'hagen_poiseuille_analytical.png', dpi=150, bbox_inches='tight')
    print(f"\nFigure saved: {output_dir / 'hagen_poiseuille_analytical.png'}")

    # Save profile data
    output_data = Path('output')
    output_data.mkdir(parents=True, exist_ok=True)

    with open(output_data / 'hagen_poiseuille_analytical.dat', 'w') as f:
        f.write(f"# Hagen-Poiseuille flow analytical solution\n")
        f.write(f"# R = {R}, dp_dx = {dp_dx}, viscosity = {viscosity}\n")
        f.write(f"# r u_analytical\n")
        for ri, ui in zip(r, u_analytical):
            f.write(f"{ri:.6f} {ui:.10e}\n")

    print(f"Data saved: {output_data / 'hagen_poiseuille_analytical.dat'}")
    plt.close()

if __name__ == '__main__':
    main()
