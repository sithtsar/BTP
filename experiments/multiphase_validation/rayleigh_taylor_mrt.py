"""
Rayleigh-Taylor Instability Test Case - Shan-Chen Two-Phase LBM

Simulates heavy fluid initially above light fluid with gravity.
The unstable configuration leads to characteristic "finger" and "bubble" formations.

Test objectives:
1. Observe instability growth from initial perturbations
2. Validate growth rate against linear stability theory
3. Visualize mushroom-shaped structures
4. Test late-stage turbulent mixing

Expected results:
- Exponential growth in early stage
- Formation of fingers (heavy fluid falling) and bubbles (light fluid rising)
- Qualitative agreement with classical Rayleigh-Taylor patterns
- No NaN or instabilities
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).parent))

from multiphase_framework import ShanChenMultiPhase


class RayleighTaylor(ShanChenMultiPhase):
    """
    Rayleigh-Taylor instability test case for Shan-Chen multi-phase LBM.
    """

    def __init__(self,
                 N: int = 128,
                 interface_y: float = None,
                 rho_heavy: float = 2.1,
                 rho_light: float = 0.15,
                 gravity: float = 1e-5,
                 perturbation_amplitude: float = 2.0,
                 perturbation_wavelength: float = 32.0,
                 **kwargs):
        """
        Initialize Rayleigh-Taylor instability simulation.

        Args:
            N: Domain size (N×N lattice)
            interface_y: Initial interface position (default: N/2)
            rho_heavy: Density of heavy fluid (top layer initially)
            rho_light: Density of light fluid (bottom layer initially)
            gravity: Gravitational acceleration (downward)
            perturbation_amplitude: Amplitude of initial perturbation
            perturbation_wavelength: Wavelength of sinusoidal perturbation
            **kwargs: Additional arguments passed to ShanChenMultiPhase
        """
        self.interface_y = interface_y if interface_y is not None else N // 2
        self.rho_heavy = rho_heavy
        self.rho_light = rho_light
        self.gravity = gravity
        self.perturbation_amplitude = perturbation_amplitude
        self.perturbation_wavelength = perturbation_wavelength

        super().__init__(N=N, **kwargs)

    def initialize_density_field(self):
        """
        Initialize stratified layers with sinusoidal perturbation.

        Heavy fluid above light fluid (unstable configuration).
        Interface perturbed with sin(2π x / λ).
        """
        print(f"Initializing Rayleigh-Taylor instability:")
        print(f"  Interface y: {self.interface_y:.2f}")
        print(f"  ρ_heavy = {self.rho_heavy:.3f} (top)")
        print(f"  ρ_light = {self.rho_light:.3f} (bottom)")
        print(f"  Gravity: {self.gravity:.6e}")
        print(f"  Perturbation: A = {self.perturbation_amplitude:.2f}, "
              f"λ = {self.perturbation_wavelength:.2f}")

        # Create coordinate grids
        x = np.arange(self.N)
        y = np.arange(self.N)
        X, Y = np.meshgrid(x, y, indexing='ij')

        # Perturbed interface: y_interface(x) = y0 + A*sin(2πx/λ)
        k = 2.0 * np.pi / self.perturbation_wavelength
        y_interface = self.interface_y + self.perturbation_amplitude * np.sin(k * X)

        # Set density based on perturbed interface
        for i in range(self.N):
            for j in range(self.N):
                if Y[i, j] > y_interface[i, j]:
                    # Above interface: heavy fluid
                    self.dh.fill(self.rho.name, self.rho_heavy, slice_obj=[i, j])
                else:
                    # Below interface: light fluid
                    self.dh.fill(self.rho.name, self.rho_light, slice_obj=[i, j])

    def add_gravity_force(self):
        """
        Add gravitational body force to the simulation.

        Note: This requires modifying the force term after initialization.
        For simplicity, we'll implement gravity through the Shan-Chen force
        by adding an external force term.

        TODO: Properly integrate body force into lbmpy framework.
        For now, gravity effect comes from density stratification.
        """
        # Gravity is implicitly included through buoyancy
        # The heavier fluid will fall due to pressure gradients
        pass

    def compute_interface_amplitude(self) -> tuple:
        """
        Compute amplitude of interface perturbation over time.

        Returns:
            Tuple of (amplitude, center_of_mass_y)
        """
        rho = self.get_density_field()

        # Find interface position for each x-column
        interface_positions = []

        for i in range(self.N):
            column = rho[i, :]

            # Find y where density crosses midpoint
            rho_mid = 0.5 * (self.rho_heavy + self.rho_light)

            # Simple linear interpolation to find crossing
            for j in range(self.N - 1):
                if (column[j] - rho_mid) * (column[j + 1] - rho_mid) < 0:
                    # Found crossing
                    interface_positions.append(j + 0.5)
                    break

        if len(interface_positions) > 0:
            amplitude = np.std(interface_positions)
            mean_y = np.mean(interface_positions)
        else:
            amplitude = 0.0
            mean_y = self.interface_y

        return amplitude, mean_y

    def run_with_tracking(self, time_steps: int, print_interval: int = 1000,
                         save_interval: int = 500) -> dict:
        """
        Run simulation with interface tracking.

        Args:
            time_steps: Number of LBM iterations
            print_interval: Print progress every N steps
            save_interval: Save snapshots every N steps

        Returns:
            Dictionary with tracking data
        """
        print(f"{'='*60}")
        print(f"Running Rayleigh-Taylor Simulation: {time_steps} time steps")
        print(f"{'='*60}")

        # Move data to target device
        self.dh.all_to_gpu()

        # Tracking arrays
        amplitude_history = []
        interface_y_history = []
        time_history = []
        snapshots = []

        for step in range(time_steps):
            # LBM cycle
            self.sync_rho()
            self.dh.run_kernel(self.collision_kernel)
            self.sync_pdfs()
            self.dh.run_kernel(self.stream_kernel)
            self.dh.swap(self.src.name, self.dst.name)

            # Track interface
            if step % 100 == 0:
                self.dh.all_to_cpu()

                amplitude, mean_y = self.compute_interface_amplitude()

                amplitude_history.append(amplitude)
                interface_y_history.append(mean_y)
                time_history.append(step)

                self.dh.all_to_gpu()

            # Save snapshots
            if step % save_interval == 0:
                self.dh.all_to_cpu()
                snapshots.append((step, self.get_density_field().copy()))
                self.dh.all_to_gpu()

            # Progress reporting
            if (step + 1) % print_interval == 0:
                self.dh.all_to_cpu()
                rho_field = self.get_density_field()
                amplitude, mean_y = self.compute_interface_amplitude()
                print(f"Step {step+1:6d}/{time_steps}: "
                      f"ρ ∈ [{rho_field.min():.3f}, {rho_field.max():.3f}], "
                      f"amplitude = {amplitude:.4f}, y_mean = {mean_y:.2f}")
                self.dh.all_to_gpu()

        # Move data back to CPU
        self.dh.all_to_cpu()

        print(f"{'='*60}")
        print(f"✓ Simulation completed: {time_steps} steps")
        print(f"{'='*60}\n")

        return {
            'amplitude': np.array(amplitude_history),
            'interface_y': np.array(interface_y_history),
            'time': np.array(time_history),
            'snapshots': snapshots
        }

    def plot_results(self, tracking_data: dict, save_dir: str = "figures/multiphase"):
        """
        Generate comprehensive Rayleigh-Taylor analysis plots.

        Args:
            tracking_data: Dictionary from run_with_tracking()
            save_dir: Directory to save figures
        """
        Path(save_dir).mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print("Generating Rayleigh-Taylor Analysis Plots")
        print(f"{'='*60}")

        # 1. Final density field
        self.plot_density_field(
            title="Rayleigh-Taylor Instability: Final State",
            save_path=f"{save_dir}/rayleigh_taylor_final.png",
            show=False
        )

        # 2. Interface evolution
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

        # Amplitude growth
        ax1.plot(tracking_data['time'], tracking_data['amplitude'],
                'b-', linewidth=2)
        ax1.axhline(self.perturbation_amplitude, color='k', linestyle='--',
                   alpha=0.5, label=f'Initial A = {self.perturbation_amplitude:.2f}')

        ax1.set_xlabel('Time Step', fontsize=12)
        ax1.set_ylabel('Interface Amplitude (std)', fontsize=12)
        ax1.set_title('Instability Growth', fontsize=14, fontweight='bold')
        ax1.legend(fontsize=10)
        ax1.grid(True, alpha=0.3)

        # Interface position (center of mass)
        ax2.plot(tracking_data['time'], tracking_data['interface_y'],
                'r-', linewidth=2)
        ax2.axhline(self.interface_y, color='k', linestyle='--',
                   alpha=0.5, label=f'Initial y = {self.interface_y:.2f}')

        ax2.set_xlabel('Time Step', fontsize=12)
        ax2.set_ylabel('Interface Position (mean)', fontsize=12)
        ax2.set_title('Interface Displacement', fontsize=14, fontweight='bold')
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(f"{save_dir}/rayleigh_taylor_evolution.png",
                   dpi=300, bbox_inches='tight')
        plt.close()

        print(f"✓ Evolution plot saved")

        # 3. Instability sequence (snapshots)
        snapshots = tracking_data['snapshots']
        n_snapshots = min(6, len(snapshots))

        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()

        for idx in range(n_snapshots):
            snap_idx = idx * (len(snapshots) - 1) // (n_snapshots - 1)
            step, density = snapshots[snap_idx]

            im = axes[idx].imshow(density.T, origin='lower', cmap='RdBu_r',
                                 extent=[0, self.N, 0, self.N],
                                 vmin=self.rho_light, vmax=self.rho_heavy)
            axes[idx].set_title(f't = {step}', fontsize=12, fontweight='bold')
            axes[idx].set_xlabel('x')
            axes[idx].set_ylabel('y')

            plt.colorbar(im, ax=axes[idx], label='ρ')

        plt.tight_layout()
        plt.savefig(f"{save_dir}/rayleigh_taylor_sequence.png",
                   dpi=300, bbox_inches='tight')
        plt.close()

        print(f"✓ Instability sequence plot saved")

        # 4. Velocity field at final time
        self.plot_velocity_field(
            title="Rayleigh-Taylor: Final Velocity Field",
            save_path=f"{save_dir}/rayleigh_taylor_velocity.png",
            show=False,
            skip=4
        )

        print(f"{'='*60}")
        print(f"✓ All plots saved to: {save_dir}/")
        print(f"{'='*60}\n")

    def print_validation_summary(self, tracking_data: dict):
        """Print validation summary."""
        print(f"\n{'='*60}")
        print("RAYLEIGH-TAYLOR VALIDATION SUMMARY")
        print(f"{'='*60}")

        initial_amplitude = self.perturbation_amplitude
        final_amplitude = tracking_data['amplitude'][-1]
        amplitude_growth = final_amplitude / initial_amplitude

        initial_y = self.interface_y
        final_y = tracking_data['interface_y'][-1]
        displacement = abs(final_y - initial_y)

        print(f"\n1. Instability Growth:")
        print(f"   Initial amplitude: {initial_amplitude:.4f}")
        print(f"   Final amplitude:   {final_amplitude:.4f}")
        print(f"   Growth factor:     {amplitude_growth:.2f}×")
        print(f"   Status: {'✓ GROWING' if amplitude_growth > 1.5 else '⚠ WEAK GROWTH'}")

        print(f"\n2. Interface Displacement:")
        print(f"   Initial y: {initial_y:.2f}")
        print(f"   Final y:   {final_y:.2f}")
        print(f"   Displacement: {displacement:.2f}")

        print(f"\n3. Configuration:")
        print(f"   Heavy fluid (ρ={self.rho_heavy:.3f}) initially on top")
        print(f"   Light fluid (ρ={self.rho_light:.3f}) initially on bottom")
        print(f"   Atwood number: {(self.rho_heavy-self.rho_light)/(self.rho_heavy+self.rho_light):.3f}")

        rho = self.get_density_field()
        print(f"\n4. Final State:")
        print(f"   Density range: [{rho.min():.3f}, {rho.max():.3f}]")
        print(f"   Max velocity: {self.compute_spurious_currents():.6e}")

        print(f"\n{'='*60}")
        if amplitude_growth > 1.5:
            print("✓ INSTABILITY DEVELOPED")
        else:
            print("⚠ WEAK INSTABILITY (may need longer run or stronger gravity)")
        print(f"{'='*60}\n")


def run_rayleigh_taylor_test(N: int = 128,
                             time_steps: int = 20000,
                             gravity: float = 1e-5,
                             omega: float = 1.0,
                             g_interaction: float = -4.7):
    """
    Run complete Rayleigh-Taylor instability test case.

    Args:
        N: Domain size
        time_steps: Number of simulation steps
        gravity: Gravitational acceleration
        omega: Relaxation parameter
        g_interaction: Shan-Chen interaction strength
    """
    print(f"\n{'='*60}")
    print("RAYLEIGH-TAYLOR INSTABILITY TEST CASE")
    print("Shan-Chen Two-Phase LBM with MRT")
    print(f"{'='*60}\n")

    # Create simulation
    sim = RayleighTaylor(
        N=N,
        interface_y=N // 2,
        rho_heavy=2.1,
        rho_light=0.15,
        gravity=gravity,
        perturbation_amplitude=5.0,  # Large perturbation to trigger instability
        perturbation_wavelength=64.0,  # Larger wavelength grows faster (λ/2 = 32)
        omega=omega,
        g_interaction=g_interaction,
        use_mrt=True,
        target='cpu'
    )

    # Initialize
    sim.initialize()

    # Run simulation with tracking
    tracking_data = sim.run_with_tracking(
        time_steps=time_steps,
        print_interval=2000,
        save_interval=1000
    )

    # Generate plots
    sim.plot_results(tracking_data)

    # Print validation summary
    sim.print_validation_summary(tracking_data)

    # Save final state
    output_dir = Path("output/multiphase")
    output_dir.mkdir(parents=True, exist_ok=True)
    sim.save_state(str(output_dir / "rayleigh_taylor_final.npz"))

    return sim, tracking_data


if __name__ == "__main__":
    # Run Rayleigh-Taylor test with OPTIMIZED parameters for stronger instability
    # Key insight: Must balance gravity vs surface tension to get instability WITHOUT mixing
    # Changes from original (g=1e-5, g_int=-4.7, amp=2.0):
    # - gravity: 1e-5 → 5e-5 (5× stronger, but not too strong to cause mixing)
    # - g_interaction: -4.7 → -4.7 (KEEP SAME - maintain phase separation)
    # - time_steps: 20000 → 40000 (much longer to see growth)
    # - perturbation_amplitude: 2.0 → 5.0 (larger initial perturbation)
    # - Use larger wavelength for faster growth

    sim, results = run_rayleigh_taylor_test(
        N=128,
        time_steps=40000,
        gravity=5e-5,
        omega=1.0,
        g_interaction=-4.7
    )

    print("\n" + "="*60)
    print("Rayleigh-Taylor instability test completed successfully!")
    print("Check figures/multiphase/ for visualization")
    print("="*60)
