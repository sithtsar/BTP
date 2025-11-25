"""
Shan-Chen Multi-Phase LBM Framework using lbmpy

This module provides a base class for implementing Shan-Chen two-phase
(liquid-gas) lattice Boltzmann simulations using lbmpy with MRT collision operator.

Based on:
- Shan & Chen (1993) - Lattice Boltzmann model for simulating flows with multiple phases
- lbmpy documentation: https://deepwiki.com/lssfau/lbmpy
"""

import numpy as np
import sympy as sp
from abc import ABC, abstractmethod
from typing import Tuple, Dict, Optional
import matplotlib.pyplot as plt

try:
    import pystencils as ps
    from pystencils import create_data_handling, create_kernel, CreateKernelConfig
    import lbmpy
    from lbmpy import LBMConfig, Method, ForceModel
    from lbmpy.enums import Stencil
    from lbmpy.stencils import LBStencil
    from lbmpy.creationfunctions import create_lb_update_rule, create_lb_method
    from lbmpy.maxwellian_equilibrium import get_weights
    from lbmpy.updatekernels import create_stream_pull_with_output_kernel
    from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
except ImportError as e:
    print(f"Error importing required packages: {e}")
    print("Please ensure lbmpy, pystencils, and sympy are installed:")
    print("  uv pip install lbmpy pystencils sympy")
    raise


class ShanChenMultiPhase(ABC):
    """
    Base class for Shan-Chen two-phase simulations using lbmpy.

    This implements the pseudo-potential approach where phase separation
    occurs through non-local density-dependent interactions.
    """

    def __init__(self,
                 N: int = 64,
                 omega: float = 1.0,
                 g_interaction: float = -4.7,
                 rho0: float = 1.0,
                 use_mrt: bool = True,
                 target: str = 'cpu'):
        """
        Initialize Shan-Chen multi-phase simulation.

        Args:
            N: Domain size (N×N lattice)
            omega: Relaxation parameter (for SRT) or base for MRT
            g_interaction: Interaction strength (negative for phase separation)
            rho0: Reference density for psi function
            use_mrt: Use MRT collision operator (True) or SRT (False)
            target: Computational target ('cpu' or 'gpu')
        """
        self.N = N
        self.omega = omega
        self.g_interaction = g_interaction
        self.rho0 = rho0
        self.use_mrt = use_mrt

        # Setup stencil (D2Q9 for 2D simulations)
        self.stencil = LBStencil(Stencil.D2Q9)
        self.D = self.stencil.D  # Dimensionality
        self.Q = self.stencil.Q  # Number of velocities

        # Get lattice weights
        self.weights = get_weights(self.stencil)

        # Setup data handling with periodic boundaries
        self.target = ps.Target.GPU if target == 'gpu' else ps.Target.CPU
        self.dh = create_data_handling(
            domain_size=(N,) * self.D,
            periodicity=True,
            default_target=self.target
        )

        # Add arrays for PDFs (src/dst for double buffering)
        self.src = self.dh.add_array('src', values_per_cell=self.Q)
        self.dst = self.dh.add_array_like('dst', 'src')

        # Add array for density field
        self.rho = self.dh.add_array('rho')

        # Velocity field (for output/analysis)
        self.u = self.dh.add_array('u', values_per_cell=self.D)

        # Create symbolic expressions for Shan-Chen force
        self._setup_force_term()

        # Create LBM collision and streaming kernels
        self._create_kernels()

        # Initialize synchronization functions
        self.sync_pdfs = self.dh.synchronization_function([self.src.name])
        self.sync_rho = self.dh.synchronization_function([self.rho.name])

        print(f"{'='*60}")
        print(f"Shan-Chen Multi-Phase Simulation Initialized")
        print(f"{'='*60}")
        print(f"Domain size:         {N}×{N}")
        print(f"Stencil:             D{self.D}Q{self.Q}")
        print(f"Collision operator:  {'MRT' if use_mrt else 'SRT'}")
        print(f"Relaxation (omega):  {omega:.4f}")
        print(f"Interaction (g):     {g_interaction:.4f}")
        print(f"Reference density:   {rho0:.4f}")
        print(f"Target:              {target.upper()}")
        print(f"{'='*60}\n")

    def _setup_force_term(self):
        """
        Setup Shan-Chen pseudo-potential force term.

        Force formula:
            F = -g * psi(rho) * sum_i[ w_i * psi(rho_i) * c_i ]

        where psi(rho) = rho0 * (1 - exp(-rho/rho0))
        """
        # Create symbolic density field access
        rho_field = self.rho.center_vector[0]

        # Define psi function symbolically
        def psi_symbolic(dens):
            return self.rho0 * (1.0 - sp.exp(-dens / self.rho0))

        # Zero vector for force accumulation
        zero_vec = sp.Matrix([0] * self.D)

        # Compute force: sum over neighbors
        force_sum = zero_vec
        for direction, weight in zip(self.stencil, self.weights):
            # Access density at neighbor
            rho_neighbor = self.rho[direction]

            # Add contribution: w_i * psi(rho_i) * c_i
            force_sum += weight * psi_symbolic(rho_neighbor) * sp.Matrix(direction)

        # Final force: -g * psi(rho_center) * sum
        self.force = -self.g_interaction * psi_symbolic(rho_field) * force_sum

        print(f"✓ Shan-Chen force term configured")
        print(f"  Pseudo-potential: psi(ρ) = {self.rho0}*(1 - exp(-ρ/{self.rho0}))")
        print(f"  Force model: Guo (for Galilean invariance)")

    def _create_kernels(self):
        """Create LBM collision and streaming kernels using lbmpy."""

        # Configure LBM method
        if self.use_mrt:
            method = Method.MRT
            print(f"\n✓ Using MRT collision operator")
        else:
            method = Method.SRT
            print(f"\n✓ Using SRT (BGK) collision operator")

        # LBM configuration with Shan-Chen force
        lbm_config = LBMConfig(
            stencil=self.stencil,
            method=method,
            relaxation_rate=self.omega,
            compressible=True,  # Required for phase separation
            force_model=ForceModel.GUO,  # Guo forcing scheme
            force=self.force,
            kernel_type='collide_only'
        )

        # Create collision update rule
        collision_rule = create_lb_update_rule(
            lbm_config=lbm_config,
            optimization={'symbolic_field': self.src}
        )

        # Store the LB method (needed for streaming)
        self.lb_method = collision_rule.method

        # Create streaming kernel with density output
        stream_rule = create_stream_pull_with_output_kernel(
            self.lb_method,
            self.src,
            self.dst,
            {'density': self.rho}
        )

        # Kernel compilation configuration
        config = CreateKernelConfig(
            target=self.target,
            cpu_openmp=False
        )

        # Compile kernels
        self.collision_kernel = create_kernel(collision_rule, config=config).compile()
        self.stream_kernel = create_kernel(stream_rule, config=config).compile()

        print(f"✓ Collision kernel compiled")
        print(f"✓ Streaming kernel compiled")

        # Create initialization kernel (without force for initial equilibrium)
        method_no_force = create_lb_method(
            lbm_config=LBMConfig(
                stencil=self.stencil,
                relaxation_rate=self.omega,
                compressible=True
            )
        )

        init_assignments = macroscopic_values_setter(
            method_no_force,
            velocity=tuple([0] * self.D),
            pdfs=self.src.center_vector,
            density=self.rho.center
        )

        self.init_kernel = create_kernel(
            init_assignments,
            ghost_layers=0,
            config=config
        ).compile()

        print(f"✓ Initialization kernel compiled\n")

    @abstractmethod
    def initialize_density_field(self):
        """
        Initialize the density field (must be implemented by subclasses).

        This should set self.rho to desired initial configuration
        (e.g., droplet, stratified layers, random fluctuations).
        """
        pass

    def initialize(self):
        """
        Initialize simulation: set density field and PDFs to equilibrium.
        """
        print(f"{'='*60}")
        print("Initializing Simulation")
        print(f"{'='*60}")

        # Set density field (implemented by subclass)
        self.initialize_density_field()

        # Initialize PDFs to equilibrium based on density
        self.dh.run_kernel(self.init_kernel)

        print("✓ Density field initialized")
        print("✓ PDFs initialized to equilibrium")
        print(f"  Density range: [{self.get_density_field().min():.3f}, "
              f"{self.get_density_field().max():.3f}]")
        print(f"{'='*60}\n")

    def run(self, time_steps: int, print_interval: int = 1000):
        """
        Run simulation for specified number of time steps.

        Args:
            time_steps: Number of LBM iterations
            print_interval: Print progress every N steps
        """
        print(f"{'='*60}")
        print(f"Running Simulation: {time_steps} time steps")
        print(f"{'='*60}")

        # Move data to target device (GPU if configured)
        self.dh.all_to_gpu()

        for step in range(time_steps):
            # LBM cycle
            self.sync_rho()                      # Sync density ghost layers
            self.dh.run_kernel(self.collision_kernel)  # Collision with force
            self.sync_pdfs()                     # Sync PDF ghost layers
            self.dh.run_kernel(self.stream_kernel)     # Streaming + compute density
            self.dh.swap(self.src.name, self.dst.name) # Swap buffers

            # Progress reporting
            if (step + 1) % print_interval == 0:
                rho_field = self.get_density_field()
                print(f"Step {step+1:6d}/{time_steps}: "
                      f"ρ ∈ [{rho_field.min():.3f}, {rho_field.max():.3f}]")

        # Move data back to CPU
        self.dh.all_to_cpu()

        print(f"{'='*60}")
        print(f"✓ Simulation completed: {time_steps} steps")
        print(f"{'='*60}\n")

    def get_density_field(self) -> np.ndarray:
        """
        Get current density field as numpy array.

        Returns:
            Density field (N×N)
        """
        return self.dh.gather_array(self.rho.name)

    def get_velocity_field(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute velocity field from current PDFs.

        Returns:
            Tuple of (u_x, u_y) velocity components
        """
        # Get PDFs and density
        pdfs = self.dh.gather_array(self.src.name)
        rho = self.get_density_field()

        # Compute velocity: u = (1/rho) * sum_i(c_i * f_i)
        u_x = np.zeros_like(rho)
        u_y = np.zeros_like(rho)

        for i in range(self.Q):
            cx, cy = self.stencil[i]
            u_x += cx * pdfs[..., i]
            u_y += cy * pdfs[..., i]

        u_x /= rho
        u_y /= rho

        return u_x, u_y

    def compute_spurious_currents(self) -> float:
        """
        Compute maximum spurious currents (velocities near interface).

        Lower values indicate better numerical stability.
        Typical threshold: < 1e-3 for stable droplet.

        Returns:
            Maximum velocity magnitude
        """
        u_x, u_y = self.get_velocity_field()
        velocity_magnitude = np.sqrt(u_x**2 + u_y**2)
        return np.max(velocity_magnitude)

    def plot_density_field(self,
                          title: str = "Density Field",
                          save_path: Optional[str] = None,
                          show: bool = True):
        """
        Plot density field as heatmap.

        Args:
            title: Plot title
            save_path: Optional path to save figure
            show: Whether to display figure
        """
        rho = self.get_density_field()

        fig, ax = plt.subplots(figsize=(8, 7))

        im = ax.imshow(
            rho.T,  # Transpose for correct orientation
            origin='lower',
            cmap='viridis',
            extent=[0, self.N, 0, self.N]
        )

        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Density ρ', fontsize=12)

        ax.set_xlabel('x', fontsize=12)
        ax.set_ylabel('y', fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold')

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ Figure saved to: {save_path}")

        if show:
            plt.show()
        else:
            plt.close()

    def plot_velocity_field(self,
                           title: str = "Velocity Field",
                           save_path: Optional[str] = None,
                           show: bool = True,
                           skip: int = 4):
        """
        Plot velocity field as quiver plot over density.

        Args:
            title: Plot title
            save_path: Optional path to save figure
            show: Whether to display figure
            skip: Plot every Nth velocity vector (for clarity)
        """
        rho = self.get_density_field()
        u_x, u_y = self.get_velocity_field()

        fig, ax = plt.subplots(figsize=(8, 7))

        # Background: density field
        im = ax.imshow(
            rho.T,
            origin='lower',
            cmap='gray',
            alpha=0.5,
            extent=[0, self.N, 0, self.N]
        )

        # Foreground: velocity vectors
        x = np.arange(0, self.N, skip)
        y = np.arange(0, self.N, skip)
        X, Y = np.meshgrid(x, y)

        U = u_x[::skip, ::skip].T
        V = u_y[::skip, ::skip].T

        ax.quiver(X, Y, U, V, color='red', alpha=0.8)

        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Density ρ', fontsize=12)

        ax.set_xlabel('x', fontsize=12)
        ax.set_ylabel('y', fontsize=12)
        ax.set_title(title, fontsize=14, fontweight='bold')

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ Figure saved to: {save_path}")

        if show:
            plt.show()
        else:
            plt.close()

    def save_state(self, filename: str):
        """
        Save current simulation state to file.

        Args:
            filename: Output filename (.npz format)
        """
        rho = self.get_density_field()
        u_x, u_y = self.get_velocity_field()

        np.savez(
            filename,
            density=rho,
            velocity_x=u_x,
            velocity_y=u_y,
            parameters={
                'N': self.N,
                'omega': self.omega,
                'g_interaction': self.g_interaction,
                'rho0': self.rho0,
                'use_mrt': self.use_mrt
            }
        )

        print(f"✓ State saved to: {filename}")


# Example usage (to be removed in production)
if __name__ == "__main__":
    print("This is a base class module. Use specific test cases:")
    print("  - static_droplet_mrt.py")
    print("  - droplet_collision_mrt.py")
    print("  - rayleigh_taylor_mrt.py")
    print("  - spinodal_decomposition_mrt.py")
