"""
lbmpy Framework for ELBM/BGK/Analytical Comparison

This module provides base classes for implementing lattice Boltzmann simulations
using lbmpy with both BGK and ELBM collision operators, and comparing results
with analytical solutions.
"""

import numpy as np
from abc import ABC, abstractmethod
from typing import Tuple, Callable, Dict, Optional
import matplotlib.pyplot as plt


class LBMPySimulation(ABC):
    """Base class for lbmpy simulations with ELBM/BGK comparison."""

    def __init__(self, nx: int, ny: int, nu: float, u_max: float = 0.1):
        """
        Initialize simulation parameters.

        Args:
            nx: Number of lattice nodes in x-direction
            ny: Number of lattice nodes in y-direction
            nu: Kinematic viscosity
            u_max: Maximum velocity (for stability)
        """
        self.nx = nx
        self.ny = ny
        self.nu = nu
        self.u_max = u_max

        # Lattice parameters (D2Q9)
        self.cs2 = 1.0 / 3.0  # Speed of sound squared
        self.q = 9  # Number of discrete velocities

        # Calculate relaxation parameter
        self.tau = 3.0 * nu + 0.5
        self.omega = 1.0 / self.tau

        # Lattice weights for D2Q9
        self.w = np.array([
            4/9,                           # 0: center
            1/9, 1/9, 1/9, 1/9,           # 1-4: cardinal directions
            1/36, 1/36, 1/36, 1/36        # 5-8: diagonal directions
        ])

        # Lattice velocities for D2Q9
        self.c = np.array([
            [0, 0],    # 0
            [1, 0],    # 1: east
            [0, 1],    # 2: north
            [-1, 0],   # 3: west
            [0, -1],   # 4: south
            [1, 1],    # 5: northeast
            [-1, 1],   # 6: northwest
            [-1, -1],  # 7: southwest
            [1, -1]    # 8: southeast
        ])

        # Initialize distribution functions
        self.f = np.zeros((self.q, ny, nx))
        self.f_eq = np.zeros((self.q, ny, nx))

        # Macroscopic fields
        self.rho = np.ones((ny, nx))
        self.u = np.zeros((ny, nx))
        self.v = np.zeros((ny, nx))

    def equilibrium_bgk(self, rho: np.ndarray, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        """
        Calculate BGK equilibrium distribution (2nd order Hermite expansion).

        Args:
            rho: Density field
            u: x-velocity field
            v: y-velocity field

        Returns:
            Equilibrium distribution f_eq[q, ny, nx]
        """
        f_eq = np.zeros((self.q, self.ny, self.nx))
        u_sqr = u**2 + v**2

        for i in range(self.q):
            cu = self.c[i, 0] * u + self.c[i, 1] * v
            f_eq[i] = self.w[i] * rho * (
                1.0
                + 3.0 * cu
                + 4.5 * cu**2
                - 1.5 * u_sqr
            )

        return f_eq

    def equilibrium_entropic(self, rho: np.ndarray, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        """
        Calculate entropic equilibrium distribution.

        Based on: f_i^eq = w_i * rho * ∏(2-√(1+u_j²)) * ((2u_j+√(1+3u_j²))/(1-u_j))^(c_ij/c)

        Args:
            rho: Density field
            u: x-velocity field
            v: y-velocity field

        Returns:
            Equilibrium distribution f_eq[q, ny, nx]
        """
        f_eq = np.zeros((self.q, self.ny, self.nx))

        # Prefactor: ∏(2-√(1+u_j²))
        eps = 1e-10  # Small epsilon for numerical stability
        prefactor = (2.0 - np.sqrt(1.0 + u**2 + eps)) * (2.0 - np.sqrt(1.0 + v**2 + eps))

        for i in range(self.q):
            # Power term for each component
            power_u = self.c[i, 0]
            power_v = self.c[i, 1]

            # Avoid division by zero
            denom_u = np.maximum(1.0 - u, eps)
            denom_v = np.maximum(1.0 - v, eps)

            term_u = ((2.0 * u + np.sqrt(1.0 + 3.0 * u**2)) / denom_u) ** power_u
            term_v = ((2.0 * v + np.sqrt(1.0 + 3.0 * v**2)) / denom_v) ** power_v

            f_eq[i] = self.w[i] * rho * prefactor * term_u * term_v

        return f_eq

    def compute_macroscopic(self):
        """Compute macroscopic density and velocity from distribution functions."""
        self.rho = np.sum(self.f, axis=0)
        self.u = np.sum(self.f * self.c[:, 0, np.newaxis, np.newaxis], axis=0) / self.rho
        self.v = np.sum(self.f * self.c[:, 1, np.newaxis, np.newaxis], axis=0) / self.rho

    def stream(self):
        """Streaming step: propagate distributions along lattice velocities."""
        for i in range(self.q):
            self.f[i] = np.roll(self.f[i], self.c[i, 0], axis=1)  # x-direction
            self.f[i] = np.roll(self.f[i], self.c[i, 1], axis=0)  # y-direction

    def collide_bgk(self):
        """BGK collision step."""
        self.f_eq = self.equilibrium_bgk(self.rho, self.u, self.v)
        self.f = self.f - self.omega * (self.f - self.f_eq)

    def collide_elbm(self):
        """
        Entropic LBM collision step (simplified two-step process).

        Step 1: Iso-entropy relaxation (find α such that H(f*) = H(f))
        Step 2: Dissipation relaxation (β-relaxation for viscosity)
        """
        # Use entropic equilibrium
        self.f_eq = self.equilibrium_entropic(self.rho, self.u, self.v)

        # Compute delta
        delta = self.f_eq - self.f

        # Find alpha bounds to ensure positivity
        alpha = self.find_alpha_isoentropic(self.f, delta)

        # Step 1: Iso-entropy relaxation
        f_star = self.f + alpha * delta

        # Step 2: Dissipation (use beta = omega for simplicity)
        beta = self.omega
        self.f = (1.0 - beta) * self.f + beta * f_star

    def find_alpha_isoentropic(self, f: np.ndarray, delta: np.ndarray,
                                max_iter: int = 10, tol: float = 1e-6) -> float:
        """
        Find alpha parameter that satisfies H(f*) = H(f) using Newton's method.

        Simplified version: uses average alpha across all lattice nodes.

        Args:
            f: Current distribution
            delta: f_eq - f
            max_iter: Maximum Newton iterations
            tol: Convergence tolerance

        Returns:
            alpha: Relaxation parameter
        """
        # Compute alpha bounds for positivity
        eps = 1e-10
        alpha_min = 0.0
        alpha_max = 2.0

        # Find bounds where f + alpha*delta > 0
        for i in range(self.q):
            mask_neg = delta[i] < 0
            if np.any(mask_neg):
                alpha_max = min(alpha_max, np.min(f[i][mask_neg] / np.abs(delta[i][mask_neg])))

            mask_pos = delta[i] > 0
            if np.any(mask_pos):
                alpha_min = max(alpha_min, np.min(-f[i][mask_pos] / delta[i][mask_pos]))

        # Ensure valid bounds
        alpha_max = min(alpha_max, 2.0)
        alpha_min = max(alpha_min, 0.0)

        if alpha_max <= alpha_min:
            return alpha_min

        # Initial guess
        alpha = 0.5 * (alpha_min + alpha_max)

        # Newton iteration to find alpha such that H(f*) = H(f)
        for _ in range(max_iter):
            f_star = f + alpha * delta

            # Ensure positivity
            f_star = np.maximum(f_star, eps)

            # Compute H(f*) - H(f)
            H_f = self.compute_entropy(f)
            H_fstar = self.compute_entropy(f_star)

            dH = H_fstar - H_f

            if abs(dH) < tol:
                break

            # Newton step: approximate derivative
            dalpha = 0.001
            f_test = f + (alpha + dalpha) * delta
            f_test = np.maximum(f_test, eps)
            H_test = self.compute_entropy(f_test)
            dHdalpha = (H_test - H_fstar) / dalpha

            if abs(dHdalpha) > eps:
                alpha = alpha - dH / dHdalpha
                alpha = np.clip(alpha, alpha_min, alpha_max)

        return float(alpha)

    def compute_entropy(self, f: np.ndarray) -> float:
        """
        Compute discrete entropy: H = Σ f_i * ln(f_i / w_i)

        Args:
            f: Distribution function

        Returns:
            Total entropy
        """
        eps = 1e-10
        H = 0.0
        for i in range(self.q):
            f_safe = np.maximum(f[i], eps)
            H += np.sum(f_safe * np.log(f_safe / self.w[i]))
        return H

    @abstractmethod
    def apply_boundary_conditions(self):
        """Apply boundary conditions (must be implemented by subclasses)."""
        pass

    @abstractmethod
    def analytical_solution(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute analytical solution (must be implemented by subclasses).

        Returns:
            Tuple of (y_coords, u_analytical)
        """
        pass

    def initialize(self):
        """Initialize distribution functions to equilibrium."""
        self.f = self.equilibrium_bgk(self.rho, self.u, self.v)

    def run_bgk(self, max_steps: int, check_interval: int = 100,
                convergence_threshold: float = 1e-3) -> Dict:
        """
        Run BGK simulation until convergence.

        Args:
            max_steps: Maximum number of time steps
            check_interval: Steps between convergence checks
            convergence_threshold: L2 error threshold for convergence

        Returns:
            Dictionary with results
        """
        print(f"\n{'='*60}")
        print("Running BGK Simulation")
        print(f"{'='*60}")

        self.initialize()
        y_analytical, u_analytical = self.analytical_solution()
        errors = []

        for step in range(max_steps):
            self.collide_bgk()
            self.stream()
            self.apply_boundary_conditions()
            self.compute_macroscopic()

            if step % check_interval == 0:
                u_sim = self.extract_profile()
                l2_error = np.sqrt(np.mean((u_sim - u_analytical)**2))
                errors.append((step, l2_error))

                print(f"Step {step:6d}: L2 error = {l2_error:.6f}")

                if l2_error < convergence_threshold:
                    print(f"✓ Converged at step {step}")
                    break

        u_final = self.extract_profile()

        return {
            'method': 'BGK',
            'y': y_analytical,
            'u': u_final,
            'errors': errors,
            'final_step': step,
            'omega': self.omega,
            'tau': self.tau
        }

    def run_elbm(self, max_steps: int, check_interval: int = 100,
                 convergence_threshold: float = 1e-3) -> Dict:
        """
        Run ELBM simulation until convergence.

        Args:
            max_steps: Maximum number of time steps
            check_interval: Steps between convergence checks
            convergence_threshold: L2 error threshold for convergence

        Returns:
            Dictionary with results
        """
        print(f"\n{'='*60}")
        print("Running ELBM Simulation")
        print(f"{'='*60}")

        self.initialize()
        y_analytical, u_analytical = self.analytical_solution()
        errors = []

        for step in range(max_steps):
            self.collide_elbm()
            self.stream()
            self.apply_boundary_conditions()
            self.compute_macroscopic()

            if step % check_interval == 0:
                u_sim = self.extract_profile()
                l2_error = np.sqrt(np.mean((u_sim - u_analytical)**2))
                errors.append((step, l2_error))

                print(f"Step {step:6d}: L2 error = {l2_error:.6f}")

                if l2_error < convergence_threshold:
                    print(f"✓ Converged at step {step}")
                    break

        u_final = self.extract_profile()

        return {
            'method': 'ELBM',
            'y': y_analytical,
            'u': u_final,
            'errors': errors,
            'final_step': step,
            'omega': self.omega,
            'tau': self.tau
        }

    @abstractmethod
    def extract_profile(self) -> np.ndarray:
        """
        Extract velocity profile for comparison (must be implemented by subclasses).

        Returns:
            1D velocity profile
        """
        pass

    def plot_comparison(self, bgk_results: Dict, elbm_results: Dict,
                       save_path: Optional[str] = None):
        """
        Generate 3-way comparison plot: BGK, ELBM, and Analytical.

        Args:
            bgk_results: Results dictionary from run_bgk()
            elbm_results: Results dictionary from run_elbm()
            save_path: Optional path to save figure
        """
        y_analytical, u_analytical = self.analytical_solution()

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        # Plot 1: Velocity profiles
        ax1.plot(y_analytical, u_analytical, 'k--', linewidth=2, label='Analytical', alpha=0.7)
        ax1.plot(bgk_results['y'], bgk_results['u'], 'b-', linewidth=1.5, label='BGK')
        ax1.plot(elbm_results['y'], elbm_results['u'], 'r-', linewidth=1.5, label='ELBM')

        ax1.set_xlabel('y', fontsize=12)
        ax1.set_ylabel('u', fontsize=12)
        ax1.set_title(f'{self.__class__.__name__}: Velocity Profiles', fontsize=14, fontweight='bold')
        ax1.legend(fontsize=10)
        ax1.grid(True, alpha=0.3)

        # Plot 2: Convergence history
        if bgk_results['errors']:
            bgk_steps, bgk_errors = zip(*bgk_results['errors'])
            ax2.semilogy(bgk_steps, bgk_errors, 'b-', linewidth=1.5, label='BGK')

        if elbm_results['errors']:
            elbm_steps, elbm_errors = zip(*elbm_results['errors'])
            ax2.semilogy(elbm_steps, elbm_errors, 'r-', linewidth=1.5, label='ELBM')

        ax2.axhline(y=1e-3, color='k', linestyle='--', linewidth=1, label='Threshold (1e-3)')
        ax2.set_xlabel('Time Step', fontsize=12)
        ax2.set_ylabel('L2 Error', fontsize=12)
        ax2.set_title('Convergence History', fontsize=14, fontweight='bold')
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"\n✓ Figure saved to: {save_path}")

        plt.show()

    def print_summary(self, bgk_results: Dict, elbm_results: Dict):
        """Print summary statistics comparing BGK and ELBM."""
        y_analytical, u_analytical = self.analytical_solution()

        bgk_error = np.sqrt(np.mean((bgk_results['u'] - u_analytical)**2))
        elbm_error = np.sqrt(np.mean((elbm_results['u'] - u_analytical)**2))

        print(f"\n{'='*60}")
        print("SUMMARY")
        print(f"{'='*60}")
        print(f"{'Method':<15} {'Final L2 Error':<20} {'Steps':<10}")
        print(f"{'-'*60}")
        print(f"{'BGK':<15} {bgk_error:<20.6e} {bgk_results['final_step']:<10}")
        print(f"{'ELBM':<15} {elbm_error:<20.6e} {elbm_results['final_step']:<10}")
        print(f"{'Analytical':<15} {'0.0':<20} {'N/A':<10}")
        print(f"{'='*60}")

        if bgk_error < 1e-3 and elbm_error < 1e-3:
            print("✓ Both methods converged successfully!")
        elif bgk_error < 1e-3:
            print("✓ BGK converged, ELBM needs more steps")
        elif elbm_error < 1e-3:
            print("✓ ELBM converged, BGK needs more steps")
        else:
            print("⚠ Neither method fully converged")
