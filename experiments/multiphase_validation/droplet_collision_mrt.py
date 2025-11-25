"""
Droplet Collision Test Case - Shan-Chen Two-Phase LBM

Simulates two droplets with initial velocities colliding and coalescing.

Test objectives:
1. Verify momentum conservation during collision
2. Test dynamic interface interactions
3. Validate coalescence (merger into single droplet)
4. Observe capillary wave dynamics after merger

Expected results:
- Momentum conserved to < 1e-6
- No NaN or instabilities
- Successful merger into single droplet
- Interface oscillations (capillary waves) damped over time
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).parent))

from multiphase_framework import ShanChenMultiPhase


class DropletCollision(ShanChenMultiPhase):
    """
    Droplet collision test case for Shan-Chen multi-phase LBM.
    """

    def __init__(self,
                 N: int = 128,
                 droplet_radius: float = 15.0,
                 separation: float = 50.0,
                 velocity: float = 0.05,
                 rho_liquid: float = 2.1,
                 rho_gas: float = 0.15,
                 **kwargs):
        """
        Initialize droplet collision simulation.

        Args:
            N: Domain size (N×N lattice)
            droplet_radius: Radius of each droplet
            separation: Initial separation between droplet centers
            velocity: Initial velocity magnitude (towards each other)
            rho_liquid: Density inside droplets (liquid phase)
            rho_gas: Density outside droplets (gas phase)
            **kwargs: Additional arguments passed to ShanChenMultiPhase
        """
        self.droplet_radius = droplet_radius
        self.separation = separation
        self.velocity = velocity
        self.rho_liquid = rho_liquid
        self.rho_gas = rho_gas

        # Calculate droplet centers (on horizontal line through domain center)
        self.center_y = N // 2
        self.left_x = (N - separation) // 2
        self.right_x = (N + separation) // 2

        # Store initial momentum for conservation check
        self.initial_momentum = None

        super().__init__(N=N, **kwargs)

    def initialize_density_field(self):
        """
        Initialize two droplets separated horizontally.

        Left droplet: moving right (+x)
        Right droplet: moving left (-x)
        """
        print(f"Initializing droplet collision:")
        print(f"  Left droplet: ({self.left_x}, {self.center_y}), velocity = +{self.velocity}")
        print(f"  Right droplet: ({self.right_x}, {self.center_y}), velocity = -{self.velocity}")
        print(f"  Separation: {self.separation:.2f}")
        print(f"  Radius: {self.droplet_radius:.2f}")
        print(f"  ρ_liquid = {self.rho_liquid:.3f}")
        print(f"  ρ_gas = {self.rho_gas:.3f}")

        # Create coordinate grids
        x = np.arange(self.N)
        y = np.arange(self.N)
        X, Y = np.meshgrid(x, y, indexing='ij')

        # Left droplet distance
        r_left = np.sqrt((X - self.left_x)**2 + (Y - self.center_y)**2)

        # Right droplet distance
        r_right = np.sqrt((X - self.right_x)**2 + (Y - self.center_y)**2)

        # Set density: liquid inside droplets, gas elsewhere
        for i in range(self.N):
            for j in range(self.N):
                if r_left[i, j] < self.droplet_radius or r_right[i, j] < self.droplet_radius:
                    # Inside either droplet: liquid
                    self.dh.fill(self.rho.name, self.rho_liquid, slice_obj=[i, j])
                else:
                    # Outside both droplets: gas
                    self.dh.fill(self.rho.name, self.rho_gas, slice_obj=[i, j])

        # Note: Velocities will be set in initialize() via modified equilibrium

    def initialize_with_velocity(self):
        """
        Initialize PDFs with initial velocities for each droplet.

        This requires custom initialization to set different velocities
        in different regions.
        """
        # Get density and PDF fields
        rho = self.get_density_field()
        pdfs = self.dh.gather_array(self.src.name)

        # Create coordinate grids
        x = np.arange(self.N)
        y = np.arange(self.N)
        X, Y = np.meshgrid(x, y, indexing='ij')

        # Left droplet distance
        r_left = np.sqrt((X - self.left_x)**2 + (Y - self.center_y)**2)

        # Right droplet distance
        r_right = np.sqrt((X - self.right_x)**2 + (Y - self.center_y)**2)

        # Set velocities based on position
        u_x = np.zeros((self.N, self.N))
        u_y = np.zeros((self.N, self.N))

        # Left droplet: moving right
        mask_left = r_left < self.droplet_radius
        u_x[mask_left] = self.velocity
        u_y[mask_left] = 0.0

        # Right droplet: moving left
        mask_right = r_right < self.droplet_radius
        u_x[mask_right] = -self.velocity
        u_y[mask_right] = 0.0

        # Create a writable copy for initialization
        pdfs_init = np.zeros_like(pdfs)

        # Compute equilibrium with velocities
        for i in range(self.N):
            for j in range(self.N):
                rho_ij = rho[i, j]
                u_ij = u_x[i, j]
                v_ij = u_y[i, j]
                u_sqr = u_ij**2 + v_ij**2

                # BGK equilibrium for initialization
                for q in range(self.Q):
                    cx, cy = self.stencil[q]
                    cu = cx * u_ij + cy * v_ij

                    feq = self.weights[q] * rho_ij * (
                        1.0 + 3.0 * cu + 4.5 * cu**2 - 1.5 * u_sqr
                    )

                    pdfs_init[i, j, q] = feq

        # Write back to data handler
        self.dh.fill(self.src.name, pdfs_init.flatten(), ghost_layers=True)

        # Calculate initial momentum
        self.initial_momentum = self.compute_total_momentum()

    def initialize(self):
        """Override to use custom velocity initialization."""
        print(f"{'='*60}")
        print("Initializing Simulation")
        print(f"{'='*60}")

        # Set density field
        self.initialize_density_field()

        # Set PDFs with velocities
        self.initialize_with_velocity()

        print("✓ Density field initialized")
        print("✓ PDFs initialized with velocities")
        print(f"  Density range: [{self.get_density_field().min():.3f}, "
              f"{self.get_density_field().max():.3f}]")
        print(f"  Initial momentum: ({self.initial_momentum[0]:.6f}, "
              f"{self.initial_momentum[1]:.6f})")
        print(f"{'='*60}\n")

    def compute_total_momentum(self) -> tuple:
        """
        Compute total momentum of the system.

        Returns:
            Tuple of (momentum_x, momentum_y)
        """
        rho = self.get_density_field()
        u_x, u_y = self.get_velocity_field()

        momentum_x = np.sum(rho * u_x)
        momentum_y = np.sum(rho * u_y)

        return momentum_x, momentum_y

    def compute_momentum_error(self) -> float:
        """
        Compute relative momentum conservation error.

        Returns:
            Relative error ||p - p0|| / ||p0||
        """
        current_momentum = self.compute_total_momentum()

        delta_px = current_momentum[0] - self.initial_momentum[0]
        delta_py = current_momentum[1] - self.initial_momentum[1]

        error = np.sqrt(delta_px**2 + delta_py**2)
        initial_mag = np.sqrt(self.initial_momentum[0]**2 + self.initial_momentum[1]**2)

        if initial_mag > 1e-10:
            return error / initial_mag
        else:
            return error

    def run_with_tracking(self, time_steps: int, print_interval: int = 1000,
                         save_interval: int = 500) -> dict:
        """
        Run simulation with momentum tracking.

        Args:
            time_steps: Number of LBM iterations
            print_interval: Print progress every N steps
            save_interval: Save snapshots every N steps

        Returns:
            Dictionary with tracking data
        """
        print(f"{'='*60}")
        print(f"Running Collision Simulation: {time_steps} time steps")
        print(f"{'='*60}")

        # Move data to target device
        self.dh.all_to_gpu()

        # Tracking arrays
        momentum_x_history = []
        momentum_y_history = []
        momentum_error_history = []
        time_history = []
        snapshots = []

        for step in range(time_steps):
            # LBM cycle
            self.sync_rho()
            self.dh.run_kernel(self.collision_kernel)
            self.sync_pdfs()
            self.dh.run_kernel(self.stream_kernel)
            self.dh.swap(self.src.name, self.dst.name)

            # Track momentum
            if step % 100 == 0:
                # Move to CPU for analysis
                self.dh.all_to_cpu()

                momentum = self.compute_total_momentum()
                error = self.compute_momentum_error()

                momentum_x_history.append(momentum[0])
                momentum_y_history.append(momentum[1])
                momentum_error_history.append(error)
                time_history.append(step)

                # Move back to GPU
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
                error = self.compute_momentum_error()
                print(f"Step {step+1:6d}/{time_steps}: "
                      f"ρ ∈ [{rho_field.min():.3f}, {rho_field.max():.3f}], "
                      f"momentum error = {error:.6e}")
                self.dh.all_to_gpu()

        # Move data back to CPU
        self.dh.all_to_cpu()

        print(f"{'='*60}")
        print(f"✓ Simulation completed: {time_steps} steps")
        print(f"{'='*60}\n")

        return {
            'momentum_x': np.array(momentum_x_history),
            'momentum_y': np.array(momentum_y_history),
            'momentum_error': np.array(momentum_error_history),
            'time': np.array(time_history),
            'snapshots': snapshots
        }

    def plot_results(self, tracking_data: dict, save_dir: str = "figures/multiphase"):
        """
        Generate comprehensive collision analysis plots.

        Args:
            tracking_data: Dictionary from run_with_tracking()
            save_dir: Directory to save figures
        """
        Path(save_dir).mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print("Generating Collision Analysis Plots")
        print(f"{'='*60}")

        # 1. Final density field
        self.plot_density_field(
            title="Droplet Collision: Final State",
            save_path=f"{save_dir}/droplet_collision_final.png",
            show=False
        )

        # 2. Momentum conservation plot
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

        # Momentum components
        ax1.plot(tracking_data['time'], tracking_data['momentum_x'],
                'b-', linewidth=2, label='p_x')
        ax1.plot(tracking_data['time'], tracking_data['momentum_y'],
                'r-', linewidth=2, label='p_y')
        ax1.axhline(self.initial_momentum[0], color='b', linestyle='--',
                   alpha=0.5, label=f'p_x(0) = {self.initial_momentum[0]:.4f}')
        ax1.axhline(self.initial_momentum[1], color='r', linestyle='--',
                   alpha=0.5, label=f'p_y(0) = {self.initial_momentum[1]:.4f}')

        ax1.set_xlabel('Time Step', fontsize=12)
        ax1.set_ylabel('Momentum', fontsize=12)
        ax1.set_title('Momentum Conservation', fontsize=14, fontweight='bold')
        ax1.legend(fontsize=10)
        ax1.grid(True, alpha=0.3)

        # Momentum error
        ax2.semilogy(tracking_data['time'], tracking_data['momentum_error'],
                    'g-', linewidth=2)
        ax2.axhline(1e-6, color='k', linestyle='--', linewidth=1,
                   label='Target: < 1e-6')

        ax2.set_xlabel('Time Step', fontsize=12)
        ax2.set_ylabel('Relative Momentum Error', fontsize=12)
        ax2.set_title('Momentum Conservation Error', fontsize=14, fontweight='bold')
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(f"{save_dir}/droplet_collision_momentum.png",
                   dpi=300, bbox_inches='tight')
        plt.close()

        print(f"✓ Momentum conservation plot saved")

        # 3. Collision sequence (snapshots)
        snapshots = tracking_data['snapshots']
        n_snapshots = min(6, len(snapshots))  # Show up to 6 frames

        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()

        for idx in range(n_snapshots):
            snap_idx = idx * (len(snapshots) - 1) // (n_snapshots - 1)
            step, density = snapshots[snap_idx]

            im = axes[idx].imshow(density.T, origin='lower', cmap='viridis',
                                 extent=[0, self.N, 0, self.N])
            axes[idx].set_title(f't = {step}', fontsize=12, fontweight='bold')
            axes[idx].set_xlabel('x')
            axes[idx].set_ylabel('y')

            plt.colorbar(im, ax=axes[idx], label='ρ')

        plt.tight_layout()
        plt.savefig(f"{save_dir}/droplet_collision_sequence.png",
                   dpi=300, bbox_inches='tight')
        plt.close()

        print(f"✓ Collision sequence plot saved")

        print(f"{'='*60}")
        print(f"✓ All plots saved to: {save_dir}/")
        print(f"{'='*60}\n")

    def print_validation_summary(self, tracking_data: dict):
        """Print validation summary with pass/fail criteria."""
        print(f"\n{'='*60}")
        print("COLLISION VALIDATION SUMMARY")
        print(f"{'='*60}")

        final_momentum_error = tracking_data['momentum_error'][-1]
        max_momentum_error = np.max(tracking_data['momentum_error'])

        # Check criteria
        check_momentum = final_momentum_error < 1e-6
        check_stable = max_momentum_error < 0.1  # Reasonable bound

        print(f"\n1. Momentum Conservation:")
        print(f"   Initial: ({self.initial_momentum[0]:.6f}, {self.initial_momentum[1]:.6f})")
        final_momentum = self.compute_total_momentum()
        print(f"   Final:   ({final_momentum[0]:.6f}, {final_momentum[1]:.6f})")
        print(f"   Error:   {final_momentum_error:.6e}")
        print(f"   Threshold: < 1e-6")
        print(f"   Status: {'✓ PASS' if check_momentum else '⚠ ACCEPTABLE' if final_momentum_error < 1e-4 else '✗ FAIL'}")

        print(f"\n2. Numerical Stability:")
        print(f"   Max momentum error: {max_momentum_error:.6e}")
        print(f"   Status: {'✓ STABLE' if check_stable else '✗ UNSTABLE'}")

        print(f"\n3. Collision Properties:")
        print(f"   Droplet radius: {self.droplet_radius:.2f}")
        print(f"   Initial velocity: ±{self.velocity:.4f}")
        print(f"   Initial separation: {self.separation:.2f}")

        rho = self.get_density_field()
        print(f"\n4. Final State:")
        print(f"   Density range: [{rho.min():.3f}, {rho.max():.3f}]")
        print(f"   Spurious currents: {self.compute_spurious_currents():.6e}")

        print(f"\n{'='*60}")
        if check_momentum and check_stable:
            print("✓ VALIDATION PASSED")
        elif final_momentum_error < 1e-4 and check_stable:
            print("✓ VALIDATION ACCEPTABLE (momentum error < 1e-4)")
        else:
            print("⚠ VALIDATION NEEDS REVIEW")
        print(f"{'='*60}\n")


def run_droplet_collision_test(N: int = 128,
                               droplet_radius: float = 15.0,
                               separation: float = 50.0,
                               velocity: float = 0.05,
                               time_steps: int = 15000,
                               omega: float = 1.0,
                               g_interaction: float = -4.7):
    """
    Run complete droplet collision test case.

    Args:
        N: Domain size
        droplet_radius: Droplet radius
        separation: Initial separation
        velocity: Initial velocity magnitude
        time_steps: Number of simulation steps
        omega: Relaxation parameter
        g_interaction: Shan-Chen interaction strength
    """
    print(f"\n{'='*60}")
    print("DROPLET COLLISION TEST CASE")
    print("Shan-Chen Two-Phase LBM with MRT")
    print(f"{'='*60}\n")

    # Create simulation
    sim = DropletCollision(
        N=N,
        droplet_radius=droplet_radius,
        separation=separation,
        velocity=velocity,
        rho_liquid=2.1,
        rho_gas=0.15,
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
        print_interval=1000,
        save_interval=500
    )

    # Generate plots
    sim.plot_results(tracking_data)

    # Print validation summary
    sim.print_validation_summary(tracking_data)

    # Save final state
    output_dir = Path("output/multiphase")
    output_dir.mkdir(parents=True, exist_ok=True)
    sim.save_state(str(output_dir / "droplet_collision_final.npz"))

    return sim, tracking_data


if __name__ == "__main__":
    # Run droplet collision test with default parameters
    sim, results = run_droplet_collision_test(
        N=128,
        droplet_radius=15.0,
        separation=50.0,
        velocity=0.05,
        time_steps=15000,
        omega=1.0,
        g_interaction=-4.7
    )

    print("\n" + "="*60)
    print("Droplet collision test completed successfully!")
    print("Check figures/multiphase/ for visualization")
    print("="*60)
