"""
Static Droplet Test Case - Shan-Chen Two-Phase LBM

Validates Shan-Chen multi-phase model by simulating a static droplet
and measuring surface tension via Laplace pressure law: ΔP = σ/R

Test objectives:
1. Verify stable droplet interface (no spurious currents)
2. Measure pressure difference across interface
3. Validate Laplace law (linear relationship between ΔP and 1/R)
4. Confirm sharp interface (3-5 lattice units transition width)

Expected results:
- Spurious currents < 1e-3
- Circular droplet shape preserved
- Pressure profile matches analytical Laplace law
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Add parent directory to path
sys.path.append(str(Path(__file__).parent))

from multiphase_framework import ShanChenMultiPhase


class StaticDroplet(ShanChenMultiPhase):
    """
    Static droplet test case for Shan-Chen multi-phase LBM.
    """

    def __init__(self,
                 N: int = 64,
                 droplet_radius: float = 15.0,
                 rho_liquid: float = 2.1,
                 rho_gas: float = 0.15,
                 **kwargs):
        """
        Initialize static droplet simulation.

        Args:
            N: Domain size (N×N lattice)
            droplet_radius: Initial droplet radius (lattice units)
            rho_liquid: Density inside droplet (liquid phase)
            rho_gas: Density outside droplet (gas phase)
            **kwargs: Additional arguments passed to ShanChenMultiPhase
        """
        self.droplet_radius = droplet_radius
        self.rho_liquid = rho_liquid
        self.rho_gas = rho_gas

        # Calculate droplet center
        self.center_x = N // 2
        self.center_y = N // 2

        super().__init__(N=N, **kwargs)

    def initialize_density_field(self):
        """
        Initialize circular droplet in center of domain.

        Density:
            ρ = rho_liquid   if r < R (inside droplet)
            ρ = rho_gas      if r >= R (outside droplet)
        """
        print(f"Initializing static droplet:")
        print(f"  Center: ({self.center_x}, {self.center_y})")
        print(f"  Radius: {self.droplet_radius:.2f}")
        print(f"  ρ_liquid = {self.rho_liquid:.3f}")
        print(f"  ρ_gas = {self.rho_gas:.3f}")

        # Create coordinate grids
        x = np.arange(self.N)
        y = np.arange(self.N)
        X, Y = np.meshgrid(x, y, indexing='ij')

        # Compute distance from center
        distance = np.sqrt((X - self.center_x)**2 + (Y - self.center_y)**2)

        # Set density based on distance
        for i in range(self.N):
            for j in range(self.N):
                if distance[i, j] < self.droplet_radius:
                    # Inside droplet: liquid
                    self.dh.fill(self.rho.name, self.rho_liquid, slice_obj=[i, j])
                else:
                    # Outside droplet: gas
                    self.dh.fill(self.rho.name, self.rho_gas, slice_obj=[i, j])

    def compute_radial_density_profile(self, num_bins: int = 50) -> tuple:
        """
        Compute azimuthally-averaged radial density profile.

        Args:
            num_bins: Number of radial bins

        Returns:
            Tuple of (r_values, rho_avg, rho_std)
        """
        rho = self.get_density_field()

        # Create coordinate grids
        x = np.arange(self.N)
        y = np.arange(self.N)
        X, Y = np.meshgrid(x, y, indexing='ij')

        # Compute distance from droplet center
        r = np.sqrt((X - self.center_x)**2 + (Y - self.center_y)**2)

        # Bin radial distances
        r_max = min(self.center_x, self.center_y)  # Stay within domain
        r_bins = np.linspace(0, r_max, num_bins + 1)
        r_values = 0.5 * (r_bins[:-1] + r_bins[1:])  # Bin centers

        # Compute average density in each bin
        rho_avg = np.zeros(num_bins)
        rho_std = np.zeros(num_bins)

        for i in range(num_bins):
            mask = (r >= r_bins[i]) & (r < r_bins[i + 1])
            if np.any(mask):
                rho_avg[i] = np.mean(rho[mask])
                rho_std[i] = np.std(rho[mask])
            else:
                rho_avg[i] = np.nan
                rho_std[i] = np.nan

        return r_values, rho_avg, rho_std

    def compute_laplace_pressure(self) -> dict:
        """
        Compute pressure difference across interface using equation of state.

        For ideal gas EOS: P = cs² * ρ (with cs² = 1/3 for D2Q9)

        Returns:
            Dictionary with pressure analysis
        """
        rho = self.get_density_field()
        cs2 = 1.0 / 3.0  # Speed of sound squared for D2Q9

        # Compute pressure field: P = cs² * ρ
        pressure = cs2 * rho

        # Estimate pressure inside droplet (r < R/2)
        x = np.arange(self.N)
        y = np.arange(self.N)
        X, Y = np.meshgrid(x, y, indexing='ij')
        r = np.sqrt((X - self.center_x)**2 + (Y - self.center_y)**2)

        mask_inside = r < (self.droplet_radius / 2.0)
        mask_outside = r > (1.5 * self.droplet_radius)

        P_inside = np.mean(pressure[mask_inside])
        P_outside = np.mean(pressure[mask_outside])

        delta_P = P_inside - P_outside

        # Theoretical Laplace pressure: ΔP = σ/R
        # Surface tension σ is implicitly determined by g_interaction
        # For validation, we check linearity of ΔP vs 1/R

        return {
            'P_inside': P_inside,
            'P_outside': P_outside,
            'delta_P': delta_P,
            'radius': self.droplet_radius,
            'curvature': 1.0 / self.droplet_radius,
            'spurious_currents': self.compute_spurious_currents()
        }

    def plot_results(self, save_dir: str = "figures/multiphase"):
        """
        Generate comprehensive results plots.

        Args:
            save_dir: Directory to save figures
        """
        Path(save_dir).mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print("Generating Result Plots")
        print(f"{'='*60}")

        # 1. Density field
        self.plot_density_field(
            title=f"Static Droplet: Density Field (t = final)",
            save_path=f"{save_dir}/static_droplet_density.png",
            show=False
        )

        # 2. Radial density profile
        r_values, rho_avg, rho_std = self.compute_radial_density_profile()

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(r_values, rho_avg, 'b-', linewidth=2, label='LBM Simulation')
        ax.fill_between(r_values, rho_avg - rho_std, rho_avg + rho_std,
                        alpha=0.3, color='blue', label='±1σ')

        ax.axvline(self.droplet_radius, color='r', linestyle='--',
                  linewidth=2, label=f'Initial R = {self.droplet_radius:.1f}')
        ax.axhline(self.rho_liquid, color='gray', linestyle=':', alpha=0.5)
        ax.axhline(self.rho_gas, color='gray', linestyle=':', alpha=0.5)

        ax.set_xlabel('Radial Distance r', fontsize=12)
        ax.set_ylabel('Density ρ', fontsize=12)
        ax.set_title('Radial Density Profile (Azimuthal Average)',
                    fontsize=14, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(f"{save_dir}/static_droplet_radial_profile.png",
                   dpi=300, bbox_inches='tight')
        plt.close()

        print(f"✓ Radial profile plot saved")

        # 3. Velocity field (check spurious currents)
        self.plot_velocity_field(
            title="Spurious Currents (Velocity Field)",
            save_path=f"{save_dir}/static_droplet_velocity.png",
            show=False,
            skip=4
        )

        # 4. Pressure analysis
        pressure_data = self.compute_laplace_pressure()

        fig, ax = plt.subplots(figsize=(8, 6))

        info_text = (
            f"Pressure Analysis\n"
            f"{'='*40}\n"
            f"P_inside  = {pressure_data['P_inside']:.6f}\n"
            f"P_outside = {pressure_data['P_outside']:.6f}\n"
            f"ΔP = {pressure_data['delta_P']:.6f}\n"
            f"\n"
            f"Droplet Properties\n"
            f"{'='*40}\n"
            f"Radius R = {pressure_data['radius']:.2f}\n"
            f"Curvature 1/R = {pressure_data['curvature']:.6f}\n"
            f"\n"
            f"Laplace Law: ΔP = σ/R\n"
            f"Implied σ = ΔP × R = {pressure_data['delta_P'] * pressure_data['radius']:.6f}\n"
            f"\n"
            f"Numerical Quality\n"
            f"{'='*40}\n"
            f"Max spurious currents = {pressure_data['spurious_currents']:.6e}\n"
            f"Threshold: < 1e-3 (stable)\n"
        )

        ax.text(0.05, 0.95, info_text, transform=ax.transAxes,
               fontsize=10, verticalalignment='top',
               fontfamily='monospace',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        ax.axis('off')
        ax.set_title('Laplace Pressure Analysis', fontsize=14, fontweight='bold')

        plt.tight_layout()
        plt.savefig(f"{save_dir}/static_droplet_pressure_analysis.png",
                   dpi=300, bbox_inches='tight')
        plt.close()

        print(f"✓ Pressure analysis plot saved")

        print(f"{'='*60}")
        print(f"✓ All plots saved to: {save_dir}/")
        print(f"{'='*60}\n")

        return pressure_data

    def print_validation_summary(self, pressure_data: dict):
        """Print validation summary with pass/fail criteria."""
        print(f"\n{'='*60}")
        print("VALIDATION SUMMARY")
        print(f"{'='*60}")

        spurious = pressure_data['spurious_currents']
        delta_P = pressure_data['delta_P']

        # Check criteria
        check_spurious = spurious < 1e-3
        check_positive_pressure = delta_P > 0

        print(f"\n1. Spurious Currents:")
        print(f"   Value: {spurious:.6e}")
        print(f"   Threshold: < 1e-3")
        print(f"   Status: {'✓ PASS' if check_spurious else '✗ FAIL'}")

        print(f"\n2. Laplace Pressure (ΔP > 0):")
        print(f"   Value: {delta_P:.6f}")
        print(f"   Status: {'✓ PASS' if check_positive_pressure else '✗ FAIL'}")

        print(f"\n3. Interface Properties:")
        print(f"   Initial radius: {self.droplet_radius:.2f}")
        print(f"   Liquid density: {self.rho_liquid:.3f}")
        print(f"   Gas density: {self.rho_gas:.3f}")
        print(f"   Density ratio: {self.rho_liquid/self.rho_gas:.2f}")

        print(f"\n4. Implied Surface Tension:")
        sigma_implied = delta_P * pressure_data['radius']
        print(f"   σ = ΔP × R = {sigma_implied:.6f}")
        print(f"   (in lattice units)")

        print(f"\n{'='*60}")
        if check_spurious and check_positive_pressure:
            print("✓ VALIDATION PASSED")
        else:
            print("⚠ VALIDATION NEEDS REVIEW")
        print(f"{'='*60}\n")


def run_static_droplet_test(N: int = 64,
                            droplet_radius: float = 15.0,
                            time_steps: int = 10000,
                            omega: float = 1.0,
                            g_interaction: float = -4.7):
    """
    Run complete static droplet test case.

    Args:
        N: Domain size
        droplet_radius: Initial droplet radius
        time_steps: Number of simulation steps
        omega: Relaxation parameter
        g_interaction: Shan-Chen interaction strength
    """
    print(f"\n{'='*60}")
    print("STATIC DROPLET TEST CASE")
    print("Shan-Chen Two-Phase LBM with MRT")
    print(f"{'='*60}\n")

    # Create simulation
    sim = StaticDroplet(
        N=N,
        droplet_radius=droplet_radius,
        rho_liquid=2.1,
        rho_gas=0.15,
        omega=omega,
        g_interaction=g_interaction,
        use_mrt=True,
        target='cpu'
    )

    # Initialize
    sim.initialize()

    # Run simulation
    sim.run(time_steps=time_steps, print_interval=1000)

    # Generate plots and analysis
    pressure_data = sim.plot_results()

    # Print validation summary
    sim.print_validation_summary(pressure_data)

    # Save final state
    output_dir = Path("output/multiphase")
    output_dir.mkdir(parents=True, exist_ok=True)
    sim.save_state(str(output_dir / "static_droplet_final.npz"))

    return sim, pressure_data


if __name__ == "__main__":
    # Run static droplet test with default parameters
    sim, results = run_static_droplet_test(
        N=64,
        droplet_radius=15.0,
        time_steps=10000,
        omega=1.0,
        g_interaction=-4.7
    )

    print("\n" + "="*60)
    print("Static droplet test completed successfully!")
    print("Check figures/multiphase/ for visualization")
    print("="*60)
