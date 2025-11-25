"""
Couette Flow Validation using lbmpy - STABLE VERSION
=====================================================
Lower wall velocity and higher viscosity for ELBM stability.
"""

import numpy as np
import pystencils as ps
from lbmpy.stencils import LBStencil
from lbmpy.enums import Stencil, Method
from lbmpy.boundaries import UBB
from lbmpy.boundaries.boundaryhandling import LatticeBoltzmannBoundaryHandling
from lbmpy.creationfunctions import LBMConfig, LBMOptimisation, create_lb_update_rule
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
from lbmpy.relaxationrates import relaxation_rate_from_lattice_viscosity
import sympy as sp
import matplotlib.pyplot as plt
from pathlib import Path


def couette_flow_validation_stable(method_type='BGK', target=ps.Target.CPU, stencil_name=Stencil.D2Q9):
    """
    Run Couette flow validation with STABLE parameters for ELBM.
    """
    # REDUCED parameters for stability
    nx = 100
    ny = 40
    U_wall = 0.01  # Reduced from 0.1
    viscosity = 0.15  # Increased from 0.1
    rho_0 = 1.0
    time_steps = 30000  # More time for convergence

    # Lattice setup
    stencil = LBStencil(stencil_name)
    omega = relaxation_rate_from_lattice_viscosity(viscosity)

    # Domain size and periodicity
    L = (nx, ny)
    periodicity = [True, False]

    # Data handling
    dh = ps.create_data_handling(L, periodicity=periodicity, default_target=target)

    src = dh.add_array('src', values_per_cell=stencil.Q)
    dst = dh.add_array_like('dst', 'src')
    rho = dh.add_array('rho', values_per_cell=1)
    u = dh.add_array('u', values_per_cell=stencil.D)

    # Initialize
    dh.fill(src.name, 0.0, ghost_layers=True)
    dh.fill(dst.name, 0.0, ghost_layers=True)
    dh.fill(rho.name, rho_0)
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
            output={'velocity': u}
        )
    else:  # BGK
        lbm_config = LBMConfig(
            stencil=stencil,
            relaxation_rate=omega,
            method=Method.SRT,
            compressible=True,
            output={'velocity': u}
        )

    lbm_opt = LBMOptimisation(symbolic_field=src, symbolic_temporary_field=dst)
    update = create_lb_update_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)
    method = update.method

    # Compile kernels
    config = ps.CreateKernelConfig(cpu_openmp=False, target=target)
    stream_collide_kernel = ps.create_kernel(update, config=config).compile()

    # Initialize PDFs to equilibrium
    init = macroscopic_values_setter(
        method,
        velocity=(0.0, 0.0),
        pdfs=src.center_vector,
        density=rho_0
    )
    init_kernel = ps.create_kernel(init, ghost_layers=0).compile()
    dh.run_kernel(init_kernel)

    # Boundary handling
    lbbh = LatticeBoltzmannBoundaryHandling(method, dh, src.name, target=target)

    wall_thickness = 1
    vel_top = sp.Matrix([U_wall, 0.0])
    vel_bottom = sp.Matrix([0.0, 0.0])

    lbbh.set_boundary(UBB(velocity=vel_top), ps.make_slice[:, -wall_thickness:])
    lbbh.set_boundary(UBB(velocity=vel_bottom), ps.make_slice[:, :wall_thickness])

    sync_pdfs = dh.synchronization_function([src.name])

    print("="*60)
    print(f"Couette Flow Validation - {method_type} (STABLE)")
    print("="*60)
    print(f"Grid: {nx} x {ny}")
    print(f"Top wall velocity: {U_wall} (REDUCED)")
    print(f"Viscosity: {viscosity} (INCREASED)")
    print(f"Relaxation rate (omega): {omega:.6f}")
    print(f"Time steps: {time_steps}")
    print()

    # Time loop
    last_max_vel = -1
    for step in range(time_steps):
        sync_pdfs()
        lbbh()
        dh.run_kernel(stream_collide_kernel)
        dh.swap(src.name, dst.name)

        if step % 100 == 0 and step > 0:
            dh.to_cpu(u.name)
            velocity = dh.gather_array(u.name)
            u_mid = velocity[nx//2, :, 0]
            u_min, u_max = np.min(u_mid), np.max(u_mid)
            max_vel = u_max

            if step % 1000 == 0:
                print(f"Step {step:5d}: u_range = [{u_min:.6f}, {u_max:.6f}]")

            # Check for NaN
            if np.isnan(max_vel):
                print(f"✗ NaN detected at step {step} - simulation unstable!")
                return None

            # Convergence check
            if last_max_vel > 0 and abs(max_vel / last_max_vel - 1) < 5e-6:
                print(f"Converged at step {step}")
                break
            last_max_vel = max_vel

    # Extract velocity profile
    dh.to_cpu(u.name)
    velocity = dh.gather_array(u.name)
    u_profile = np.mean(velocity[:, :, 0], axis=0)

    # Analytical solution
    H = ny - 1
    y = np.arange(ny)
    u_analytical = U_wall * y / H

    # L2 error
    error = np.sqrt(np.mean((u_profile - u_analytical)**2))

    print()
    print(f"Final L2 error: {error:.6e}")
    print()

    results = {
        'method': method_type,
        'y': y,
        'u_sim': u_profile,
        'u_analytical': u_analytical,
        'error': error,
        'params': {
            'nx': nx,
            'ny': ny,
            'U_wall': U_wall,
            'viscosity': viscosity,
            'time_steps': time_steps,
            'omega': omega
        }
    }

    return results


if __name__ == "__main__":
    import sys

    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent.parent

    print("\n" + "="*70)
    print("COUETTE FLOW (STABLE): BGK vs ELBM vs ANALYTICAL")
    print("="*70 + "\n")

    # Run BGK
    try:
        bgk_results = couette_flow_validation_stable(method_type='BGK')
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
        elbm_results = couette_flow_validation_stable(method_type='ELBM')
        if elbm_results:
            print(f"✓ ELBM validation complete (L2 error: {elbm_results['error']:.6e})")
        else:
            print("✗ ELBM validation produced NaN - still unstable")
            sys.exit(1)
    except Exception as e:
        print(f"✗ ELBM validation FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    y = bgk_results['y']
    u_analytical = bgk_results['u_analytical']
    u_bgk = bgk_results['u_sim']
    u_elbm = elbm_results['u_sim']

    ax1.plot(u_analytical, y, 'k--', linewidth=2, label='Analytical', alpha=0.8)
    ax1.plot(u_bgk, y, 'b-', linewidth=1.5, label='BGK', alpha=0.8)
    ax1.plot(u_elbm, y, 'r-', linewidth=1.5, label='ELBM', alpha=0.8)
    ax1.set_xlabel('Velocity u', fontsize=12)
    ax1.set_ylabel('y position', fontsize=12)
    ax1.set_title('Couette Flow (Stable): Velocity Profiles', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)

    error_bgk = np.abs(u_bgk - u_analytical)
    error_elbm = np.abs(u_elbm - u_analytical)

    ax2.semilogy(error_bgk, y, 'b-', linewidth=1.5, label='BGK error', alpha=0.8)
    ax2.semilogy(error_elbm, y, 'r-', linewidth=1.5, label='ELBM error', alpha=0.8)
    ax2.set_xlabel('Absolute Error', fontsize=12)
    ax2.set_ylabel('y position', fontsize=12)
    ax2.set_title('Error Comparison', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3, which='both')

    l2_bgk = bgk_results['error']
    l2_elbm = elbm_results['error']
    fig.suptitle(f'Couette Flow (Stable) | BGK L2: {l2_bgk:.2e} | ELBM L2: {l2_elbm:.2e}',
                 fontsize=13, fontweight='bold', y=1.02)

    plt.tight_layout()

    output_dir = repo_root / 'figures' / 'lbmpy_validation'
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / 'couette_stable_3way_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved: {output_file}")

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"BGK L2 error:  {bgk_results['error']:.6e}")
    print(f"ELBM L2 error: {elbm_results['error']:.6e}")
    print("="*70)
