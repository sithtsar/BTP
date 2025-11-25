"""
2D Channel Flow: Reynolds Number Evolution
===========================================
Shows how the parabolic velocity profile changes with increasing Re.
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


def channel_2d_at_re(target_re, method_type='BGK', target=ps.Target.CPU, stencil_name=Stencil.D2Q9):
    """
    Run 2D channel flow at specific Reynolds number.

    Parameters:
    -----------
    target_re : float
        Target Reynolds number
    method_type : str
        'BGK' or 'ELBM'
    target : ps.Target
        Computation target
    stencil_name : Stencil
        Lattice stencil

    Returns:
    --------
    dict : Results with velocity profile and parameters
    """
    # Physical parameters - reduced resolution for higher Re
    width = 32
    length = 128

    # For force-driven flow: Re = u_max * H / nu
    # We'll adjust force to achieve target Re while keeping nu moderate
    lattice_viscosity = 0.1  # Fixed viscosity
    omega = relaxation_rate_from_lattice_viscosity(lattice_viscosity)

    # Calculate force needed for target Re
    # u_max = target_re * nu / H, and u_max = F * H^2 / (8 * rho * nu)
    # So: F = 8 * rho * nu * target_re * nu / (H * H^2) = 8 * rho * nu^2 * target_re / H^3
    H = width - 2  # Channel height (excluding walls)
    rho_fluid = 1.0
    force_density = 8.0 * rho_fluid * lattice_viscosity**2 * target_re / (H**3)

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
            force_model=ForceModel.LUO,
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

    print(f"Running Re={target_re} ({method_type})...", end=' ', flush=True)

    # Time loop - reduced for faster runs at higher Re
    time_steps = 30000
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
        if step % 500 == 0 and step > 0:
            dh.to_cpu(u.name)
            vel_arr = dh.gather_array(u.name)
            u_vel = vel_arr[:, :, 0]
            max_vel = np.max(np.abs(u_vel))

            # Check for NaN
            if np.isnan(max_vel):
                print(f"NaN at step {step}")
                return None

            # Convergence check
            if step > 5000:
                change = abs(max_vel - last_max_vel) / (last_max_vel + 1e-12)
                if change < 1e-8:
                    print(f"converged at step {step}")
                    break
            last_max_vel = max_vel

    # Extract final velocity profile
    dh.to_cpu(u.name)
    vel_arr = dh.gather_array(u.name)
    u_vel = vel_arr[:, :, 0]
    u_profile = np.mean(u_vel, axis=0)

    # Analytical solution
    y = np.arange(width)
    u_analytical = np.zeros(width)
    for i in range(1, width - 1):
        y_pos = i - 1
        u_analytical[i] = (force_density / (2.0 * rho_fluid * lattice_viscosity)) * y_pos * (H - y_pos)

    # Compute L2 error
    u_inner = u_profile[1:-1]
    u_ana_inner = u_analytical[1:-1]
    error = np.sqrt(np.mean((u_inner - u_ana_inner)**2))

    print(f"L2 error: {error:.6e}")

    results = {
        'method': method_type,
        'y': y,
        'u_sim': u_profile,
        'u_analytical': u_analytical,
        'error': error,
        'Re': target_re,
        'params': {
            'width': width,
            'length': length,
            'viscosity': lattice_viscosity,
            'force_density': force_density,
            'omega': omega
        }
    }

    return results


def plot_re_evolution(all_results):
    """Plot how velocity profile evolves with Reynolds number."""
    fig = plt.figure(figsize=(18, 6))
    gs = fig.add_gridspec(1, 3, hspace=0.25, wspace=0.3)

    ax1 = fig.add_subplot(gs[0])  # BGK profiles at different Re
    ax2 = fig.add_subplot(gs[1])  # ELBM profiles at different Re
    ax3 = fig.add_subplot(gs[2])  # Error comparison

    # Color map for different Re values
    colors = plt.cm.viridis(np.linspace(0, 1, len(all_results['BGK'])))

    re_values = sorted(all_results['BGK'].keys())

    # Plot 1: BGK velocity profiles at different Re
    for i, re in enumerate(re_values):
        results = all_results['BGK'][re]
        y = results['y']
        u_sim = results['u_sim']
        u_analytical = results['u_analytical']

        # Plot analytical (dashed) and simulation (solid with markers)
        ax1.plot(u_analytical * 1000, y, '--', color=colors[i], linewidth=1.5, alpha=0.5)
        ax1.plot(u_sim * 1000, y, 'o-', color=colors[i], markersize=4, markevery=3,
                linewidth=2, label=f'Re={re}', alpha=0.8)

    ax1.set_xlabel('Velocity u (×10⁻³)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('y position', fontsize=12, fontweight='bold')
    ax1.set_title('BGK: Velocity Profile Evolution with Re', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10, framealpha=0.95, loc='upper left')
    ax1.grid(True, alpha=0.3, linestyle='--')

    # Plot 2: ELBM velocity profiles at different Re
    for i, re in enumerate(re_values):
        results = all_results['ELBM'][re]
        y = results['y']
        u_sim = results['u_sim']
        u_analytical = results['u_analytical']

        # Plot analytical (dashed) and simulation (solid with markers)
        ax2.plot(u_analytical * 1000, y, '--', color=colors[i], linewidth=1.5, alpha=0.5)
        ax2.plot(u_sim * 1000, y, 's-', color=colors[i], markersize=4, markevery=3,
                linewidth=2, label=f'Re={re}', alpha=0.8)

    ax2.set_xlabel('Velocity u (×10⁻³)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('y position', fontsize=12, fontweight='bold')
    ax2.set_title('ELBM: Velocity Profile Evolution with Re', fontsize=13, fontweight='bold')
    ax2.legend(fontsize=10, framealpha=0.95, loc='upper left')
    ax2.grid(True, alpha=0.3, linestyle='--')

    # Plot 3: L2 error vs Re
    bgk_errors = [all_results['BGK'][re]['error'] for re in re_values]
    elbm_errors = [all_results['ELBM'][re]['error'] for re in re_values]

    ax3.semilogy(re_values, bgk_errors, 'o--', color='#1f77b4', markersize=8,
                linewidth=2, label='BGK (SRT)', alpha=0.8)
    ax3.semilogy(re_values, elbm_errors, 's:', color='#d62728', markersize=8,
                linewidth=2, label='ELBM (Entropic)', alpha=0.8)

    ax3.set_xlabel('Reynolds Number', fontsize=12, fontweight='bold')
    ax3.set_ylabel('L2 Error', fontsize=12, fontweight='bold')
    ax3.set_title('Convergence Error vs Reynolds Number', fontsize=13, fontweight='bold')
    ax3.legend(fontsize=11, framealpha=0.95)
    ax3.grid(True, alpha=0.3, which='both', linestyle='--')

    # Main title
    fig.suptitle('2D Channel Flow: Parabolic Profile Evolution with Reynolds Number',
                 fontsize=15, fontweight='bold', y=0.98)

    # Save using Path(__file__)
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent.parent
    output_dir = repo_root / 'figures' / 'lbmpy_validation'
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / 'channel_2d_re_evolution.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved: {output_file}")
    plt.close()


if __name__ == "__main__":
    import sys

    print("\n" + "="*70)
    print("2D CHANNEL FLOW: REYNOLDS NUMBER EVOLUTION")
    print("="*70 + "\n")

    # Reynolds numbers to test - higher range to test stability
    re_values = [20, 30, 40, 50, 60]

    all_results = {
        'BGK': {},
        'ELBM': {}
    }

    # Run BGK for all Re values
    print("BGK Simulations:")
    print("-" * 70)
    for re in re_values:
        try:
            results = channel_2d_at_re(target_re=re, method_type='BGK')
            if results:
                all_results['BGK'][re] = results
            else:
                print(f"✗ BGK at Re={re} FAILED")
                sys.exit(1)
        except Exception as e:
            print(f"✗ BGK at Re={re} FAILED: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)

    print()

    # Run ELBM for all Re values
    print("ELBM Simulations:")
    print("-" * 70)
    for re in re_values:
        try:
            results = channel_2d_at_re(target_re=re, method_type='ELBM')
            if results:
                all_results['ELBM'][re] = results
            else:
                print(f"✗ ELBM at Re={re} FAILED")
                sys.exit(1)
        except Exception as e:
            print(f"✗ ELBM at Re={re} FAILED: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)

    # Plot comparison
    plot_re_evolution(all_results)

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    for re in re_values:
        bgk_err = all_results['BGK'][re]['error']
        elbm_err = all_results['ELBM'][re]['error']
        print(f"Re={re:2d}: BGK L2={bgk_err:.3e}, ELBM L2={elbm_err:.3e}")

    print("\n✓ All simulations completed successfully!")
    sys.exit(0)
