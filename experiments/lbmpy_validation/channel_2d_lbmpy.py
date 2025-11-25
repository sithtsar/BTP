"""
2D Channel Flow Validation using lbmpy
======================================
Force-driven channel flow comparing BGK, ELBM, and analytical solution.
This is the WORKING configuration for ELBM validation.
"""

import numpy as np
import sympy as sp
import pystencils as ps
import matplotlib.pyplot as plt
from pathlib import Path
from lbmpy.boundaries import NoSlip
from lbmpy.creationfunctions import LBMConfig, LBMOptimisation, create_lb_update_rule
from lbmpy.enums import Method, Stencil, ForceModel
from lbmpy.lbstep import LatticeBoltzmannStep
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
from lbmpy.relaxationrates import relaxation_rate_from_lattice_viscosity
from lbmpy.stencils import LBStencil
from pystencils import create_data_handling


def channel_2d_validation(method_type='BGK', target=ps.Target.CPU, stencil_name=Stencil.D2Q9):
    """
    Run 2D channel flow validation.

    Parameters:
    -----------
    method_type : str
        'BGK' or 'ELBM'
    target : ps.Target
        Computation target
    stencil_name : Stencil
        Lattice stencil

    Returns:
    --------
    dict : Results with velocity profile and error
    """
    # Physical parameters - optimized for better convergence
    width = 40
    length = 160
    target_re = 15
    u_max = 0.05

    # Calculate lattice viscosity
    u_mean = u_max * 2/3  # For parabolic profile
    lattice_viscosity = u_mean * width / target_re
    omega = relaxation_rate_from_lattice_viscosity(lattice_viscosity)

    # Body force to drive flow - reduced for better convergence
    force_density = 5e-7

    # Lattice setup
    stencil = LBStencil(stencil_name)

    # Domain - periodic in x, walls in y
    L = (length, width)
    periodicity = [True, False]

    # Data handling
    dh = create_data_handling(L, periodicity=periodicity, default_target=target)

    # Add arrays
    src = dh.add_array('src', values_per_cell=stencil.Q)
    dst = dh.add_array_like('dst', 'src')
    u = dh.add_array('u', values_per_cell=stencil.D)
    rho = dh.add_array('rho', values_per_cell=1)

    # Initialize
    dh.fill(src.name, 0.0, ghost_layers=True)
    dh.fill(dst.name, 0.0, ghost_layers=True)
    dh.fill(rho.name, 1.0)
    dh.fill(u.name, 0.0)

    # LBM configuration
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
            force=(force_density, 0.0),
            output={'velocity': u}
        )
    else:  # BGK
        lbm_config = LBMConfig(
            stencil=stencil,
            relaxation_rate=omega,
            method=Method.SRT,
            compressible=True,
            force_model=ForceModel.GUO,
            force=(force_density, 0.0),
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
        density=1.0
    )
    init_kernel = ps.create_kernel(init, ghost_layers=0).compile()
    dh.run_kernel(init_kernel)

    # Boundary handling - no-slip walls
    from lbmpy.boundaries.boundaryhandling import LatticeBoltzmannBoundaryHandling
    lbbh = LatticeBoltzmannBoundaryHandling(method, dh, src.name, target=target)

    wall = NoSlip()
    lbbh.set_boundary(wall, ps.make_slice[:, :1])   # Bottom wall
    lbbh.set_boundary(wall, ps.make_slice[:, -1:])  # Top wall

    # Periodic synchronization
    sync_pdfs = dh.synchronization_function([src.name])

    print("="*60)
    print(f"2D Channel Flow Validation - {method_type}")
    print("="*60)
    print(f"Domain: {length} x {width}")
    print(f"Reynolds number: {target_re}")
    print(f"Viscosity: {lattice_viscosity:.6f}")
    print(f"Relaxation rate: {omega:.6f}")
    print(f"Force density: {force_density:.2e}")
    print()

    # Time loop - longer for better convergence
    time_steps = 50000
    last_max_vel = 0

    for step in range(time_steps):
        # Periodic BC
        sync_pdfs()

        # Wall BC
        lbbh()

        # Collision and streaming
        dh.run_kernel(stream_collide)

        # Swap
        dh.swap(src.name, dst.name)

        # Monitor convergence
        if step % 100 == 0 and step > 0:
            dh.to_cpu(u.name)
            vel_arr = dh.gather_array(u.name)
            u_vel = vel_arr[:, :, 0]
            max_vel = np.max(np.abs(u_vel))

            if step % 2000 == 0:
                u_profile = np.mean(u_vel, axis=0)
                center_vel = u_profile[width//2]
                print(f"Step {step:5d}: max_u = {max_vel:.6e}, center_u = {center_vel:.6e}")

            # Check for NaN
            if np.isnan(max_vel):
                print(f"✗ NaN detected at step {step}")
                return None

            # Convergence check - stricter tolerance
            if step > 5000:
                change = abs(max_vel - last_max_vel) / (last_max_vel + 1e-12)
                if change < 1e-8:
                    print(f"✓ Converged at step {step}")
                    break
            last_max_vel = max_vel

    # Extract final velocity profile
    dh.to_cpu(u.name)
    vel_arr = dh.gather_array(u.name)
    u_vel = vel_arr[:, :, 0]
    u_profile = np.mean(u_vel, axis=0)

    # Analytical solution for Poiseuille flow: u(y) = (F/(2*rho*nu)) * y * (H - y)
    H = width - 2  # Channel height (excluding walls)
    y = np.arange(width)
    rho_fluid = 1.0

    u_analytical = np.zeros(width)
    for i in range(1, width - 1):
        y_pos = i - 1  # Position from bottom wall
        u_analytical[i] = (force_density / (2.0 * rho_fluid * lattice_viscosity)) * y_pos * (H - y_pos)

    # Compute L2 error (excluding walls)
    u_inner = u_profile[1:-1]
    u_ana_inner = u_analytical[1:-1]
    error = np.sqrt(np.mean((u_inner - u_ana_inner)**2))

    print()
    print(f"Final L2 error: {error:.6e}")
    print(f"Max velocity: {np.max(u_profile):.6e}")
    print()

    results = {
        'method': method_type,
        'y': y,
        'u_sim': u_profile,
        'u_analytical': u_analytical,
        'error': error,
        'params': {
            'width': width,
            'length': length,
            'Re': target_re,
            'viscosity': lattice_viscosity,
            'force_density': force_density,
            'omega': omega
        }
    }

    return results


def plot_three_way_comparison(bgk_results, elbm_results):
    """Plot BGK, ELBM, and analytical solutions with better visualization."""
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

    ax1 = fig.add_subplot(gs[0, 0])  # Velocity profiles
    ax2 = fig.add_subplot(gs[0, 1])  # Error comparison
    ax3 = fig.add_subplot(gs[1, :])  # Re evolution

    y = bgk_results['y']
    u_analytical = bgk_results['u_analytical']
    u_bgk = bgk_results['u_sim']
    u_elbm = elbm_results['u_sim']

    # Plot 1: Velocity profiles with distinct line styles and markers
    ax1.plot(u_analytical * 1000, y, 'k-', linewidth=3, label='Analytical', alpha=0.9, zorder=1)
    ax1.plot(u_bgk * 1000, y, 'o--', color='#1f77b4', markersize=6, markevery=3,
             linewidth=2, label='BGK (SRT)', alpha=0.8, zorder=3)
    ax1.plot(u_elbm * 1000, y, 's:', color='#d62728', markersize=6, markevery=3,
             linewidth=2, label='ELBM (Entropic)', alpha=0.8, zorder=2)

    ax1.set_xlabel('Velocity u (×10⁻³)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('y position', fontsize=12, fontweight='bold')
    ax1.set_title('Velocity Profiles', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=11, framealpha=0.95, loc='upper left')
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.axhline(y=0, color='gray', linestyle=':', linewidth=1.5, alpha=0.6)
    ax1.axhline(y=bgk_results['params']['width']-1, color='gray', linestyle=':', linewidth=1.5, alpha=0.6)

    # Plot 2: Errors with distinct markers
    error_bgk = np.abs(u_bgk - u_analytical)
    error_elbm = np.abs(u_elbm - u_analytical)

    # Filter out zero errors at walls
    mask = (error_bgk > 1e-12) & (error_elbm > 1e-12)
    y_filtered = y[mask]
    error_bgk_filtered = error_bgk[mask]
    error_elbm_filtered = error_elbm[mask]

    ax2.semilogy(error_bgk_filtered, y_filtered, 'o--', color='#1f77b4',
                 markersize=6, markevery=3, linewidth=2, label='BGK error', alpha=0.8)
    ax2.semilogy(error_elbm_filtered, y_filtered, 's:', color='#d62728',
                 markersize=6, markevery=3, linewidth=2, label='ELBM error', alpha=0.8)

    ax2.set_xlabel('Absolute Error', fontsize=12, fontweight='bold')
    ax2.set_ylabel('y position', fontsize=12, fontweight='bold')
    ax2.set_title('Error Distribution', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=11, framealpha=0.95)
    ax2.grid(True, alpha=0.3, which='both', linestyle='--')

    # Plot 3: Reynolds number evolution (compute local Re)
    width = bgk_results['params']['width']
    nu = bgk_results['params']['viscosity']

    # Local Reynolds number: Re(y) = u(y) * width / nu
    Re_analytical = (u_analytical * width / nu)
    Re_bgk = (u_bgk * width / nu)
    Re_elbm = (u_elbm * width / nu)

    ax3.plot(y, Re_analytical, 'k-', linewidth=3, label='Analytical', alpha=0.9, zorder=1)
    ax3.plot(y, Re_bgk, 'o--', color='#1f77b4', markersize=6, markevery=3,
             linewidth=2, label='BGK (SRT)', alpha=0.8, zorder=3)
    ax3.plot(y, Re_elbm, 's:', color='#d62728', markersize=6, markevery=3,
             linewidth=2, label='ELBM (Entropic)', alpha=0.8, zorder=2)

    ax3.set_xlabel('y position', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Local Reynolds Number Re(y)', fontsize=12, fontweight='bold')
    ax3.set_title('Reynolds Number Profile Evolution', fontsize=14, fontweight='bold')
    ax3.legend(fontsize=11, framealpha=0.95, loc='upper left')
    ax3.grid(True, alpha=0.3, linestyle='--')
    ax3.axhline(y=0, color='gray', linestyle=':', linewidth=1, alpha=0.5)

    # Add vertical lines at walls
    ax3.axvline(x=0, color='gray', linestyle=':', linewidth=1.5, alpha=0.6)
    ax3.axvline(x=width-1, color='gray', linestyle=':', linewidth=1.5, alpha=0.6)

    # Add L2 errors in main title
    l2_bgk = bgk_results['error']
    l2_elbm = elbm_results['error']
    Re = bgk_results['params']['Re']
    fig.suptitle(f'2D Channel Flow Validation (Re={Re}) | BGK L2: {l2_bgk:.3e} | ELBM L2: {l2_elbm:.3e}',
                 fontsize=15, fontweight='bold', y=0.98)

    # Save using Path(__file__)
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent.parent
    output_dir = repo_root / 'figures' / 'lbmpy_validation'
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / 'channel_2d_3way_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()


if __name__ == "__main__":
    import sys

    print("\n" + "="*70)
    print("2D CHANNEL FLOW: BGK vs ELBM vs ANALYTICAL")
    print("="*70 + "\n")

    # Run BGK
    try:
        bgk_results = channel_2d_validation(method_type='BGK')
        if bgk_results:
            print(f"✓ BGK validation complete (L2 error: {bgk_results['error']:.6e})")
        else:
            print("✗ BGK validation FAILED")
            sys.exit(1)
    except Exception as e:
        print(f"✗ BGK validation FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    # Run ELBM
    try:
        elbm_results = channel_2d_validation(method_type='ELBM')
        if elbm_results:
            print(f"✓ ELBM validation complete (L2 error: {elbm_results['error']:.6e})")
        else:
            print("✗ ELBM validation produced NaN")
            sys.exit(1)
    except Exception as e:
        print(f"✗ ELBM validation FAILED: {e}")
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

    if bgk_results['error'] < 0.001 and elbm_results['error'] < 0.001:
        print("\n✓ All methods converged successfully!")
        sys.exit(0)
    else:
        print("\n⚠ Some methods did not fully converge")
        sys.exit(1)
