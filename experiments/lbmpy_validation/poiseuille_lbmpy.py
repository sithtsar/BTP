"""
Poiseuille Flow Validation using lbmpy
======================================
Implements body-force driven channel flow with parabolic velocity profile.
Compares BGK, ELBM, and analytical solution: u(y) = -(dp/dx)*y*(H-y)/(2*rho*nu)
"""

import numpy as np
import pystencils as ps
import sympy as sp
from lbmpy.advanced_streaming import LBMPeriodicityHandling
from lbmpy.analytical_solutions import poiseuille_flow
from lbmpy.boundaries import NoSlip
from lbmpy.boundaries.boundaryhandling import LatticeBoltzmannBoundaryHandling
from lbmpy.creationfunctions import create_lb_update_rule, LBMConfig, LBMOptimisation
from lbmpy.enums import Method, ForceModel, Stencil
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
from lbmpy.stencils import LBStencil
from lbmpy.relaxationrates import relaxation_rate_from_lattice_viscosity
import matplotlib.pyplot as plt
from pathlib import Path


def poiseuille_flow_validation(method_type='BGK', target=ps.Target.CPU, stencil_name=Stencil.D2Q9):
    """
    Run Poiseuille flow validation with lbmpy.

    Parameters:
    -----------
    method_type : str
        'BGK' or 'ELBM'
    target : ps.Target
        Computation target (CPU or GPU)
    stencil_name : Stencil
        Lattice stencil to use

    Returns:
    --------
    dict : Results containing velocity profiles and errors
    """
    # Physical parameters (matching C++ implementation)
    nx = 100
    ny = 40
    rho_0 = 1.0
    viscosity = 0.1  # kinematic viscosity
    dp_dx = -1e-5  # pressure gradient (negative: decreases in +x)
    time_steps = 20000

    # Convert pressure gradient to force density
    # F = -dp/dx (force per unit volume)
    ext_force_density = -dp_dx

    # Lattice setup
    stencil = LBStencil(stencil_name)
    omega = relaxation_rate_from_lattice_viscosity(viscosity)

    # Domain size and periodicity (periodic in x, non-periodic in y)
    L = (nx, ny)
    periodicity = [True, False]

    # Data handling
    dh = ps.create_data_handling(L, periodicity=periodicity, default_target=target)

    # Add arrays
    src = dh.add_array('src', values_per_cell=stencil.Q)
    dst = dh.add_array_like('dst', 'src')
    u = dh.add_array('u', values_per_cell=stencil.D)

    # Initialize
    dh.fill(src.name, 0.0, ghost_layers=True)
    dh.fill(dst.name, 0.0, ghost_layers=True)
    dh.fill(u.name, 0.0, ghost_layers=True)

    # LBM configuration with body force
    if method_type == 'ELBM':
        omega_free = sp.Symbol('omega_free')
        lbm_config = LBMConfig(
            stencil=stencil,
            relaxation_rates=[omega, omega_free],
            method=Method.TRT,
            entropic=True,
            entropic_newton_iterations=3,
            zero_centered=False,
            compressible=True,
            force_model=ForceModel.LUO,  # LUO compatible with entropic
            force=tuple([ext_force_density, 0.0]),
            output={'velocity': u}
        )
    else:  # BGK
        lbm_config = LBMConfig(
            stencil=stencil,
            relaxation_rate=omega,
            method=Method.SRT,
            compressible=True,
            force_model=ForceModel.GUO,
            force=tuple([ext_force_density, 0.0]),
            output={'velocity': u}
        )

    lbm_opt = LBMOptimisation(symbolic_field=src, symbolic_temporary_field=dst)
    update = create_lb_update_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)
    method = update.method

    # Compile kernel
    config = ps.CreateKernelConfig(cpu_openmp=False, target=target)
    stream_collide = ps.create_kernel(update, config=config).compile()

    # Initialize PDFs to equilibrium at rest
    init = macroscopic_values_setter(
        method,
        velocity=(0.0, 0.0),
        pdfs=src.center_vector,
        density=rho_0
    )
    init_kernel = ps.create_kernel(init).compile()
    dh.run_kernel(init_kernel)

    # Boundary handling
    lbbh = LatticeBoltzmannBoundaryHandling(method, dh, src.name, target=target)

    # No-slip walls on top and bottom
    noslip = NoSlip()
    wall_thickness = 1
    lbbh.set_boundary(noslip, ps.make_slice[:, :wall_thickness])
    lbbh.set_boundary(noslip, ps.make_slice[:, -wall_thickness:])

    # Periodic boundary handling for x-direction
    sync_pdfs = LBMPeriodicityHandling(stencil, dh, src.name)

    print("="*60)
    print(f"Poiseuille Flow Validation - {method_type} (lbmpy)")
    print("="*60)
    print(f"Grid: {nx} x {ny}")
    print(f"Pressure gradient (dp/dx): {dp_dx}")
    print(f"Body force density: {ext_force_density}")
    print(f"Viscosity: {viscosity}")
    print(f"Relaxation rate (omega): {omega:.6f}")
    print(f"Method: {method_type}")
    print(f"Time steps: {time_steps}")
    print()

    # Time loop with proper convergence checking
    last_max_vel = -1
    for step in range(time_steps):
        # Apply periodic BC (before everything)
        sync_pdfs()

        # Apply wall BC
        lbbh()

        # Collision and streaming
        dh.run_kernel(stream_collide)

        # Swap arrays
        dh.swap(src.name, dst.name)

        # Progress output and convergence check
        if step % 100 == 0 and step > 0:
            dh.to_cpu(u.name)
            velocity = dh.gather_array(u.name)

            # Average over x-direction (should be uniform for Poiseuille)
            u_profile = np.mean(velocity[:, :, 0], axis=0)

            max_vel = np.max(u_profile)
            u_min, u_max = np.min(u_profile), np.max(u_profile)

            if step % 1000 == 0:
                print(f"Step {step:5d}: u_range = [{u_min:.6e}, {u_max:.6e}]")

            # Check convergence (velocity change < 5e-6)
            if last_max_vel > 0 and abs(max_vel / last_max_vel - 1) < 5e-6:
                print(f"Converged at step {step}")
                break
            last_max_vel = max_vel

    # Extract final velocity profile
    dh.to_cpu(u.name)
    velocity = dh.gather_array(u.name)

    # Average over x-direction
    u_profile = np.mean(velocity[:, :, 0], axis=0)

    # Correct for F/2 term from Guo forcing scheme
    u_profile -= ext_force_density / (2.0 * rho_0)

    # Remove wall regions for comparison
    u_profile_inner = u_profile[wall_thickness:-wall_thickness]
    y_inner = np.arange(len(u_profile_inner))

    # Analytical solution using lbmpy's poiseuille_flow function
    # poiseuille_flow(y_pos, channel_width, force_density, dynamic_viscosity)
    H = ny - 2 * wall_thickness  # actual channel width
    mid = (H - 1) / 2.0
    y_centered = y_inner - mid
    dynamic_viscosity = rho_0 * viscosity  # mu = rho * nu
    u_analytical = poiseuille_flow(y_centered, H, ext_force_density, dynamic_viscosity)

    # Compute L2 error
    error = np.sqrt(np.mean((u_profile_inner - u_analytical)**2))

    print()
    print(f"Final L2 error: {error:.6e}")
    print()

    # Prepare full profiles for output (including walls)
    y_full = np.arange(ny)
    u_analytical_full = np.zeros(ny)
    u_analytical_full[wall_thickness:-wall_thickness] = u_analytical

    # Save results
    results = {
        'method': method_type,
        'y': y_full,
        'u_sim': u_profile,
        'u_analytical': u_analytical_full,
        'error': error,
        'params': {
            'nx': nx,
            'ny': ny,
            'dp_dx': dp_dx,
            'viscosity': viscosity,
            'time_steps': time_steps,
            'omega': omega,
            'force_density': ext_force_density
        }
    }

    return results


def plot_three_way_comparison(bgk_results, elbm_results):
    """Plot BGK, ELBM, and analytical solutions together."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    y = bgk_results['y']
    u_analytical = bgk_results['u_analytical']
    u_bgk = bgk_results['u_sim']
    u_elbm = elbm_results['u_sim']

    # Plot 1: Velocity profiles
    ax1.plot(u_analytical, y, 'k--', linewidth=2, label='Analytical', alpha=0.8)
    ax1.plot(u_bgk, y, 'b-', linewidth=1.5, label='BGK', alpha=0.8)
    ax1.plot(u_elbm, y, 'r-', linewidth=1.5, label='ELBM', alpha=0.8)
    ax1.set_xlabel('Velocity u', fontsize=12)
    ax1.set_ylabel('y position', fontsize=12)
    ax1.set_title('Poiseuille Flow: Velocity Profiles', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    # Plot 2: Errors
    error_bgk = np.abs(u_bgk - u_analytical)
    error_elbm = np.abs(u_elbm - u_analytical)

    # Filter out zero errors at walls
    mask = (error_bgk > 1e-12) & (error_elbm > 1e-12)
    y_filtered = y[mask]
    error_bgk_filtered = error_bgk[mask]
    error_elbm_filtered = error_elbm[mask]

    ax2.semilogy(error_bgk_filtered, y_filtered, 'b-', linewidth=1.5, label='BGK error', alpha=0.8)
    ax2.semilogy(error_elbm_filtered, y_filtered, 'r-', linewidth=1.5, label='ELBM error', alpha=0.8)
    ax2.set_xlabel('Absolute Error', fontsize=12)
    ax2.set_ylabel('y position', fontsize=12)
    ax2.set_title('Error Comparison', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3, which='both')

    # Add L2 errors in title
    l2_bgk = bgk_results['error']
    l2_elbm = elbm_results['error']
    fig.suptitle(f'Poiseuille Flow | BGK L2: {l2_bgk:.2e} | ELBM L2: {l2_elbm:.2e}',
                 fontsize=13, fontweight='bold', y=1.02)

    plt.tight_layout()

    # Save figure - use Path(__file__) to find repo root
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent.parent
    output_dir = repo_root / 'figures' / 'lbmpy_validation'
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / 'poiseuille_3way_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved: {output_file}")
    plt.close()


if __name__ == "__main__":
    import sys
    import os
    from pathlib import Path

    # Use Path(__file__) to find repo root and create directories
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent.parent
    output_dir = repo_root / 'output' / 'lbmpy_validation'
    figures_dir = repo_root / 'figures' / 'lbmpy_validation'
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "="*70)
    print("POISEUILLE FLOW: BGK vs ELBM vs ANALYTICAL")
    print("="*70 + "\n")

    # Run BGK validation
    try:
        bgk_results = poiseuille_flow_validation(method_type='BGK')
        print(f"✓ BGK validation complete (L2 error: {bgk_results['error']:.6e})")
    except Exception as e:
        print(f"✗ BGK validation FAILED with exception: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    # Run ELBM validation
    try:
        elbm_results = poiseuille_flow_validation(method_type='ELBM')
        print(f"✓ ELBM validation complete (L2 error: {elbm_results['error']:.6e})")
    except Exception as e:
        print(f"✗ ELBM validation FAILED with exception: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    # Plot comparison
    plot_three_way_comparison(bgk_results, elbm_results)

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"BGK L2 error:  {bgk_results['error']:.6e}")
    print(f"ELBM L2 error: {elbm_results['error']:.6e}")

    if bgk_results['error'] < 0.01 and elbm_results['error'] < 0.01:
        print("\n✓ All methods converged successfully!")
        sys.exit(0)
    else:
        print("\n⚠ Some methods did not fully converge")
        sys.exit(1)
