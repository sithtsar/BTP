#!/usr/bin/env python3
"""
Stokes Shear Flow and Kolmogorov Flow
Analytical solutions for both cases
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def stokes_shear_analytical(y, shear_rate):
    """
    Stokes shear flow: u(y) = γ * y
    """
    return shear_rate * y

def kolmogorov_analytical(y, L, force_amplitude, k, viscosity):
    """
    Kolmogorov flow: u(y) = (F/k²ν) * sin(ky)
    """
    y_phys = L * y / len(y)
    return (force_amplitude / (k**2 * viscosity)) * np.sin(k * y_phys)

def main():
    print("="*60)
    print("Extended Analytical Solutions: Stokes & Kolmogorov Flows")
    print("="*60)

    # Create figure with 2x2 layout
    fig = plt.figure(figsize=(16, 12))

    # ===== STOKES SHEAR FLOW =====
    print("\n1. Stokes Shear Flow")
    print("-" * 40)

    H = 40.0
    shear_rate = 0.01
    viscosity_stokes = 0.1
    rho = 1.0

    y_stokes = np.linspace(0, H, 100)
    u_stokes = stokes_shear_analytical(y_stokes, shear_rate)

    # Theoretical shear stress
    tau_theory = rho * viscosity_stokes * shear_rate

    print(f"Channel height: {H}")
    print(f"Shear rate γ: {shear_rate}")
    print(f"Viscosity: {viscosity_stokes}")
    print(f"Theoretical shear stress τ = μγ: {tau_theory:.6f}")

    # Plot 1: Velocity profile
    ax1 = plt.subplot(2, 2, 1)
    ax1.plot(y_stokes, u_stokes, 'b-', linewidth=2, label='Analytical')
    ax1.set_xlabel('y position')
    ax1.set_ylabel('Velocity u(y)')
    ax1.set_title('(a) Stokes Shear Flow:\nLinear Velocity Profile')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Plot 2: Shear rate (should be constant)
    ax2 = plt.subplot(2, 2, 2)
    du_dy = np.gradient(u_stokes, y_stokes)
    ax2.plot(y_stokes, du_dy, 'g-', linewidth=2, label='Numerical du/dy')
    ax2.axhline(y=shear_rate, color='r', linestyle='--', linewidth=2, label=f'Theory: {shear_rate:.6f}')
    ax2.set_xlabel('y position')
    ax2.set_ylabel('Shear rate du/dy')
    ax2.set_title('(b) Shear Rate Verification\n(Constant Throughout)')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # ===== KOLMOGOROV FLOW =====
    print("\n2. Kolmogorov Flow")
    print("-" * 40)

    ny_kolm = 128
    force_amplitude = 0.0001
    viscosity_kolm = 0.01
    k = 2.0 * np.pi / ny_kolm  # One wavelength
    L = 2.0 * np.pi / k

    y_kolm = np.arange(ny_kolm)
    u_kolm = kolmogorov_analytical(y_kolm, L, force_amplitude, k, viscosity_kolm)

    print(f"Domain size: {ny_kolm}")
    print(f"Wave number k: {k:.6f}")
    print(f"Force amplitude: {force_amplitude}")
    print(f"Viscosity: {viscosity_kolm}")
    print(f"Max velocity: {np.max(np.abs(u_kolm)):.6e}")

    # Plot 3: Velocity profile
    ax3 = plt.subplot(2, 2, 3)
    ax3.plot(y_kolm, u_kolm, 'b-', linewidth=2, label='Analytical')
    ax3.set_xlabel('y position')
    ax3.set_ylabel('Velocity u(y)')
    ax3.set_title('(c) Kolmogorov Flow:\nSinusoidal Forcing Pattern')
    ax3.grid(True, alpha=0.3)
    ax3.legend()

    # Plot 4: 2D visualization of forcing pattern
    ax4 = plt.subplot(2, 2, 4)
    nx_kolm = 128
    X, Y = np.meshgrid(np.arange(nx_kolm), y_kolm)
    U_field = np.tile(u_kolm[:, np.newaxis], (1, nx_kolm))

    contour = ax4.contourf(X, Y, U_field, levels=20, cmap='RdBu_r')
    ax4.set_xlabel('x position')
    ax4.set_ylabel('y position')
    ax4.set_title('(d) Kolmogorov Flow:\n2D Velocity Field')
    ax4.set_aspect('equal')
    plt.colorbar(contour, ax=ax4, label='Velocity')

    fig.suptitle('Extended Analytical Cases: Stokes & Kolmogorov Flows',
                 fontsize=16, fontweight='bold')
    plt.tight_layout()

    # Save figure
    output_dir = Path('figures/extended_analytical')
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / 'stokes_kolmogorov_analytical.png', dpi=150, bbox_inches='tight')
    print(f"\nFigure saved: {output_dir / 'stokes_kolmogorov_analytical.png'}")

    # Save data
    output_data = Path('output')
    output_data.mkdir(parents=True, exist_ok=True)

    # Stokes data
    with open(output_data / 'stokes_shear_analytical.dat', 'w') as f:
        f.write(f"# Stokes shear flow analytical solution\n")
        f.write(f"# shear_rate = {shear_rate}, tau_theory = {tau_theory:.6f}\n")
        f.write(f"# y u du_dy\n")
        for yi, ui, di in zip(y_stokes, u_stokes, du_dy):
            f.write(f"{yi:.6f} {ui:.10e} {di:.10e}\n")

    # Kolmogorov data
    with open(output_data / 'kolmogorov_analytical.dat', 'w') as f:
        f.write(f"# Kolmogorov flow analytical solution\n")
        f.write(f"# k = {k:.6f}, force_amplitude = {force_amplitude}\n")
        f.write(f"# y u\n")
        for yi, ui in zip(y_kolm, u_kolm):
            f.write(f"{yi:.6f} {ui:.10e}\n")

    print(f"Data saved: {output_data / 'stokes_shear_analytical.dat'}")
    print(f"Data saved: {output_data / 'kolmogorov_analytical.dat'}")

    plt.close()

    print("\n" + "="*60)
    print("All analytical solutions generated successfully!")
    print("="*60)

if __name__ == '__main__':
    main()
