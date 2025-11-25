"""
Entropic LBM (ELBM) Multi-Phase Framework using Shan-Chen Model

CRITICAL NOTE: Due to lbmpy limitations, we cannot use ForceModel.GUO with entropic methods.
Instead, we use ForceModel.LUO which is compatible with entropic collision operators.

Key differences from MRT/SRT implementation:
1. Force model: GUO → LUO (compatible with entropic)
2. Collision operator: TRT_KBC_N4 with entropic=True
3. Entropy maximization: Ensures unconditional stability
4. Two relaxation rates: omega_s (viscosity) and omega_h (entropy-optimized)

Reference:
- lbmpy documentation on entropic methods
- Shan & Chen (1993) for pseudo-potential model
- Ansumali & Karlin (2000) for entropic LBM
"""

import numpy as np
import sympy as sp
from abc import ABC, abstractmethod
import matplotlib.pyplot as plt
from pathlib import Path

from lbmpy.creationfunctions import create_lb_collision_rule, create_lb_method
from lbmpy.enums import Method, Stencil, ForceModel
from lbmpy.stencils import LBStencil
from lbmpy.maxwellian_equilibrium import get_weights
from lbmpy.advanced_streaming import is_inplace

from pystencils import create_kernel, Field
from pystencils.datahandling import create_data_handling, SerialDataHandling


class ShanChenEntropicMultiPhase(ABC):
    """
    Base class for Shan-Chen multi-phase simulations using Entropic LBM.

    Key features:
    - Entropic collision operator (unconditionally stable)
    - Shan-Chen pseudo-potential forces
    - ForceModel.LUO (compatible with entropic methods)
    - Automatic entropy maximization via alpha parameter
    """

    def __init__(self,
                 N: int = 64,
                 omega: float = 1.0,
                 g_interaction: float = -4.7,
                 rho0: float = 1.0,
                 entropic_newton_iters: int = 3,
                 target: str = 'cpu'):
        """
        Initialize Shan-Chen multi-phase ELBM simulation.

        Args:
            N: Domain size (N×N lattice)
            omega: Shear relaxation rate (determines viscosity)
            g_interaction: Shan-Chen interaction strength (negative = attraction)
            rho0: Reference density for pseudo-potential
            entropic_newton_iters: Newton iterations for entropy optimization
            target: 'cpu' or 'gpu'
        """
        self.N = N
        self.omega_s = omega  # Shear viscosity rate
        self.g_interaction = g_interaction
        self.rho0 = rho0
        self.entropic_iters = entropic_newton_iters
        self.target = target

        # Setup lattice
        self.stencil = LBStencil(Stencil.D2Q9)
        self.q = len(self.stencil)
        self.dim = 2

        # Get weights from stencil using get_weights()
        self.weights = get_weights(self.stencil)

        # Fields for pystencils
        self.dh = None
        self.pdf_field = None
        self.velocity_field = None
        self.density_field = None

        # Compiled kernels
        self.collision_kernel = None
        self.stream_kernel = None

        # Create symbolic force term (Shan-Chen)
        self._setup_force_term()

        # Create LBM kernels with ELBM
        self._create_entropic_kernels()

    def _setup_force_term(self):
        """
        Setup Shan-Chen pseudo-potential force term.

        Force: F = -g * psi(rho_center) * sum_i[ w_i * psi(rho_i) * c_i ]
        Pseudo-potential: psi(rho) = rho0 * (1 - exp(-rho / rho0))
        """
        print("✓ Setting up Shan-Chen force term for ELBM")
        print(f"  Pseudo-potential: psi(ρ) = {self.rho0}*(1 - exp(-ρ/{self.rho0}))")
        print(f"  Force model: LUO (compatible with entropic)")

        # Symbolic density field
        rho_field = sp.Symbol('rho')

        # Pseudo-potential function
        def psi_symbolic(dens):
            return self.rho0 * (1.0 - sp.exp(-dens / self.rho0))

        # Access to neighbor densities
        from pystencils.field import Field
        self.rho = Field.create_generic('rho', spatial_dimensions=2, index_dimensions=0)

        # Compute interaction force
        zero_vec = sp.Matrix([0, 0])
        force_sum = zero_vec

        for direction, weight in zip(self.stencil, self.weights):
            rho_neighbor = self.rho[direction]
            c_i = sp.Matrix(direction)
            force_sum += weight * psi_symbolic(rho_neighbor) * c_i

        # Final force: F = -g * psi(rho_center) * sum
        self.force = -self.g_interaction * psi_symbolic(rho_field) * force_sum

    def _create_entropic_kernels(self):
        """
        Create ELBM collision and streaming kernels.

        Uses TRT_KBC_N4 method with entropic=True and ForceModel.LUO.
        The free relaxation rate (omega_h) is automatically determined
        by entropy maximization.
        """
        print("\n✓ Creating Entropic LBM kernels")
        print(f"  Method: TRT_KBC_N4 (entropic)")
        print(f"  Omega_s: {self.omega_s} (shear viscosity)")
        print(f"  Omega_h: FREE (entropy-optimized)")
        print(f"  Entropy Newton iterations: {self.entropic_iters}")

        # Create LB method with entropic collision
        lb_method = create_lb_method(
            stencil=self.stencil,
            method=Method.TRT_KBC_N4,
            relaxation_rates=[self.omega_s, sp.Symbol("omega_free")],
            compressible=True,
            zero_centered=False,  # Required for entropic methods
            force_model=ForceModel.LUO,  # Compatible with entropic
            force=tuple(self.force),
        )

        # PDF field (population field)
        self.src = Field.create_generic('src', spatial_dimensions=2,
                                       index_dimensions=1,
                                       index_shape=(self.q,))
        self.dst = Field.create_generic('dst', spatial_dimensions=2,
                                       index_dimensions=1,
                                       index_shape=(self.q,))

        # Create collision rule with ELBM
        collision_rule = create_lb_collision_rule(
            lb_method=lb_method,
            src_field=self.src,
            dst_field=self.dst,
            entropic=True,  # Enable entropy maximization
            entropic_newton_iterations=self.entropic_iters,
        )

        print("✓ Entropic collision rule created")

        # Compile collision kernel
        self.collision_kernel = create_kernel(
            collision_rule,
            target=self.target,
            cpu_openmp=False
        ).compile()

        print("✓ Collision kernel compiled")

        # Create streaming kernel (simple advection)
        from lbmpy.updatekernels import create_stream_pull_with_output_kernel

        self.stream_kernel = create_stream_pull_with_output_kernel(
            lb_method.stencil,
            src_field=self.src,
            dst_field=self.dst,
            output={
                'density': self.rho,
                'velocity': Field.create_generic('u', spatial_dimensions=2,
                                                index_dimensions=1,
                                                index_shape=(2,))
            }
        )

        self.stream_kernel = create_kernel(
            self.stream_kernel,
            target=self.target,
            cpu_openmp=False
        ).compile()

        print("✓ Streaming kernel compiled")

    def initialize(self):
        """Initialize simulation data structures and fields."""
        print(f"\n{'='*60}")
        print("Initializing Entropic Multi-Phase Simulation")
        print(f"{'='*60}")

        # Create data handling
        self.dh = create_data_handling(
            domain_size=(self.N, self.N),
            periodicity=(True, True),
            default_ghost_layers=1
        )

        # Add fields
        self.dh.add_array('src', values_per_cell=self.q)
        self.dh.add_array('dst', values_per_cell=self.q)
        self.dh.add_array('rho', values_per_cell=1)
        self.dh.add_array('u', values_per_cell=2)

        print(f"Domain size:         {self.N}×{self.N}")
        print(f"Stencil:             D2Q9")
        print(f"Collision operator:  ELBM (TRT_KBC_N4)")
        print(f"Omega_s:             {self.omega_s:.4f}")
        print(f"Interaction (g):     {self.g_interaction:.4f}")
        print(f"Reference density:   {self.rho0:.4f}")
        print(f"Target:              {self.target.upper()}")
        print(f"{'='*60}\n")

        # Initialize density field (subclass implements this)
        self.initialize_density_field()

        # Initialize PDFs to equilibrium
        self._initialize_pdfs()

        print("✓ Initialization complete\n")

    @abstractmethod
    def initialize_density_field(self):
        """Initialize density field (implemented by subclasses)."""
        pass

    def _initialize_pdfs(self):
        """Initialize PDFs to equilibrium based on density field."""
        print("Initializing PDFs to equilibrium...")

        # Get density array
        rho = self.dh.cpu_arrays['rho']
        pdfs = self.dh.cpu_arrays['src']

        # Initialize with equilibrium (zero velocity)
        for i in range(self.N):
            for j in range(self.N):
                density = rho[i, j]

                # Equilibrium with zero velocity
                for q in range(self.q):
                    w = float(self.weights[q])
                    pdfs[i, j, q] = w * density

        # Copy to dst
        self.dh.cpu_arrays['dst'][:] = pdfs[:]

        # Get density range for validation
        density_min = rho.min()
        density_max = rho.max()

        print(f"✓ PDFs initialized to equilibrium")
        print(f"  Density range: [{density_min:.3f}, {density_max:.3f}]")

    def get_density_field(self) -> np.ndarray:
        """Get current density field."""
        return self.dh.cpu_arrays['rho'][:, :, 0]

    def get_velocity_field(self) -> tuple:
        """Get current velocity field (u_x, u_y)."""
        u = self.dh.cpu_arrays['u']
        return u[:, :, 0], u[:, :, 1]

    def sync_pdfs(self):
        """Sync PDF field for force calculation."""
        self.dh.to_cpu('src')

    def sync_rho(self):
        """Sync density field."""
        self.dh.to_cpu('rho')

    def run(self, time_steps: int, print_interval: int = 1000):
        """
        Run ELBM simulation.

        Args:
            time_steps: Number of time steps
            print_interval: Print progress every N steps
        """
        print(f"{'='*60}")
        print(f"Running ELBM Simulation: {time_steps} time steps")
        print(f"{'='*60}")

        # Move data to target device
        if self.target == 'gpu':
            self.dh.all_to_gpu()

        for step in range(time_steps):
            # ELBM cycle (collision includes entropy maximization)
            self.sync_rho()
            self.dh.run_kernel(self.collision_kernel)

            self.sync_pdfs()
            self.dh.run_kernel(self.stream_kernel)

            # Swap src/dst
            self.dh.swap(self.src.name, self.dst.name)

            # Print progress
            if (step + 1) % print_interval == 0:
                self.dh.to_cpu('rho')
                rho = self.get_density_field()
                print(f"Step {step+1:6d}/{time_steps}: "
                      f"ρ ∈ [{rho.min():.3f}, {rho.max():.3f}]")

        # Move data back to CPU
        self.dh.all_to_cpu()

        print(f"{'='*60}")
        print(f"✓ ELBM simulation completed: {time_steps} steps")
        print(f"{'='*60}\n")

    def compute_spurious_currents(self) -> float:
        """Compute maximum spurious velocity magnitude."""
        u_x, u_y = self.get_velocity_field()
        velocity_magnitude = np.sqrt(u_x**2 + u_y**2)
        return velocity_magnitude.max()

    def plot_density_field(self, title: str = "Density Field",
                          save_path: str = None, show: bool = True):
        """Plot density field."""
        rho = self.get_density_field()

        plt.figure(figsize=(10, 8))
        plt.imshow(rho.T, origin='lower', cmap='viridis',
                  extent=[0, self.N, 0, self.N])
        plt.colorbar(label='Density ρ')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(title, fontsize=14, fontweight='bold')
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ Figure saved to: {save_path}")

        if show:
            plt.show()
        else:
            plt.close()

    def save_state(self, filename: str):
        """Save simulation state."""
        rho = self.get_density_field()
        u_x, u_y = self.get_velocity_field()

        np.savez(filename,
                density=rho,
                velocity_x=u_x,
                velocity_y=u_y,
                omega=self.omega_s,
                g_interaction=self.g_interaction,
                N=self.N)

        print(f"✓ State saved to: {filename}")
