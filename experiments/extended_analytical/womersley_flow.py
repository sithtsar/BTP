#!/usr/bin/env python3
"""
Womersley Flow: Oscillatory pressure-driven flow between parallel plates
Analytical solution validation
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def womersley_analytical(y, H, t, amplitude, frequency, viscosity):
    """
    Womersley analytical solution
    u(y,t) = Re[A * (1 - cosh(λy)/cosh(λH/2)) * exp(iωt)]
    """
    omega = 2 * np.pi * frequency

    # Complex wave number
    lambda_val = np.sqrt(1j * omega / viscosity)

    # Center coordinates
    y_centered = y - H / 2.0

    # Complex velocity
    u_complex = amplitude * (1.0 - np.cosh(lambda_val * y_centered) / np.cosh(lambda_val * H / 2.0)) * np.exp(1j * omega * t)

    return np.real(u_complex)

def main():
    # Parameters
    H = 40.0  # Channel height
    amplitude = 0.01
    frequency = 0.01  # Hz
    viscosity = 0.1

    # Womersley number
    omega = 2 * np.pi * frequency
    alpha = H * np.sqrt(omega / viscosity)

    print(f"Womersley Flow Analysis")
    print(f"=" * 60)
    print(f"Channel height: {H}")
    print(f"Oscillation amplitude: {amplitude}")
    print(f"Frequency: {frequency} Hz")
    print(f"Viscosity: {viscosity}")
    print(f"Womersley number α: {alpha:.4f}")
    print(f"=" * 60)

    # Grid
    y = np.linspace(0, H-1, 100)

    # Time points over one period
    period = 1.0 / frequency
    n_snapshots = 8
    times = np.linspace(0, period, n_snapshots, endpoint=False)

    # Create figure
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    for idx, t in enumerate(times):
        u = womersley_analytical(y, H, t, amplitude, frequency, viscosity)

        ax = axes[idx]
        ax.plot(u, y, 'b-', linewidth=2)
        ax.set_xlabel('Velocity')
        ax.set_ylabel('y position')
        ax.set_title(f't = {t/period:.3f}T')
        ax.grid(True, alpha=0.3)
        ax.axvline(x=0, color='k', linestyle='--', alpha=0.3)

    fig.suptitle(f'Womersley Flow (α = {alpha:.2f}): Velocity Profiles Over One Period',
                 fontsize=14, fontweight='bold')

    plt.tight_layout()

    # Save
    output_dir = Path('figures/extended_analytical')
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / 'womersley_flow_analytical.png', dpi=150, bbox_inches='tight')
    print(f"\nFigure saved: {output_dir / 'womersley_flow_analytical.png'}")

    # Also save profile data
    output_data = Path('output')
    output_data.mkdir(parents=True, exist_ok=True)

    with open(output_data / 'womersley_analytical.dat', 'w') as f:
        f.write(f"# Womersley flow analytical solution\n")
        f.write(f"# alpha = {alpha:.4f}\n")
        f.write(f"# y ")
        for i, t in enumerate(times):
            f.write(f"u_t{i} ")
        f.write("\n")

        for yi in y:
            f.write(f"{yi:.6f} ")
            for t in times:
                u = womersley_analytical(np.array([yi]), H, t, amplitude, frequency, viscosity)[0]
                f.write(f"{u:.10e} ")
            f.write("\n")

    print(f"Data saved: {output_data / 'womersley_analytical.dat'}")
    plt.close()

if __name__ == '__main__':
    main()
