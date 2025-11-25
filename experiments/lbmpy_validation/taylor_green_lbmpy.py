"""
Taylor-Green Vortex Validation using lbmpy
==========================================
Implements decaying 2D vortex for viscosity extraction validation.
Compares BGK, ELBM, and analytical energy decay: E(t) = E0 * exp(-4*nu*k^2*t)
"""

import numpy as np
import pystencils as ps
import sympy as sp
from lbmpy.stencils import LBStencil
from lbmpy.enums import Stencil, Method
from lbmpy.creationfunctions import LBMConfig, LBMOptimisation, create_lb_update_rule, create_lb_method
from lbmpy.macroscopic_value_kernels import compile_macroscopic_values_setter
from lbmpy.relaxationrates import relaxation_rate_from_lattice_viscosity
import matplotlib.pyplot as plt
from pathlib import Path


def taylor_green_vortex_validation(method_type='BGK', target=ps.Target.CPU, stencil_name=Stencil.D2Q9):
    """
    Run Taylor-Green vortex validation with lbmpy.

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
    dict : Results containing energy decay and viscosity extraction
    """
    # Physical parameters (matching C++ implementation)
    nx = 64
    ny = 64
    U0 = 0.1
    viscosity = 0.01
    k = 2.0 * np.pi / nx
    rho_0 = 1.0
    time_steps = 1000
    output_interval = 100

    # Lattice setup
    stencil = LBStencil(stencil_name)
    omega = relaxation_rate_from_lattice_viscosity(viscosity)

    # Domain size (fully periodic)
    L = (nx, ny)
    periodicity = [True, True]

    # Data handling
    dh = ps.create_data_handling(L, periodicity=periodicity, default_target=target)

    # Add arrays
    src = dh.add_array('src', values_per_cell=stencil.Q)
    dst = dh.add_array_like('dst', 'src')
    u = dh.add_array('u', values_per_cell=stencil.D)
    rho = dh.add_array('rho', values_per_cell=1)

    # Initialize
    dh.fill(src.name, 0.0, ghost_layers=True)
    dh.fill(dst.name, 0.0, ghost_layers=True)

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
            output={'velocity': u, 'density': rho}
        )
    else:  # BGK
        lbm_config = LBMConfig(
            stencil=stencil,
            relaxation_rate=omega,
            method=Method.SRT,
            compressible=True,
            output={'velocity': u, 'density': rho}
        )

    lbm_opt = LBMOptimisation(symbolic_field=src, symbolic_temporary_field=dst)
    update = create_lb_update_rule(lbm_config=lbm_config, lbm_optimisation=lbm_opt)
    method = update.method

    # Compile kernel
    config = ps.CreateKernelConfig(cpu_openmp=False, target=target)
    stream_collide = ps.create_kernel(update, config=config).compile()

    print("="*60)
    print(f"Taylor-Green Vortex Validation - {method_type} (lbmpy)")
    print("="*60)
    print(f"Grid: {nx} x {ny}")
    print(f"Initial velocity amplitude: {U0}")
    print(f"Viscosity: {viscosity}")
    print(f"Wave number k: {k}")
    print(f"Relaxation rate (omega): {omega:.6f}")
    print(f"Method: {method_type}")
    print(f"Time steps: {time_steps}")
    print()

    # Initialize Taylor-Green vortex velocity field
    print("Initializing Taylor-Green vortex...")

    # Get array shapes (including ghost layers)
    rho_arr = dh.cpu_arrays[rho.name]
    u_arr = dh.cpu_arrays[u.name]

    # Determine ghost layer size (typically 1)
    ghost_layers = 1
    inner_slice = (slice(ghost_layers, -ghost_layers), slice(ghost_layers, -ghost_layers))

    # Create spatially varying density and velocity fields
    density_field = np.ones((nx, ny))
    velocity_field = np.zeros((nx, ny, stencil.D))

    for i in range(nx):
        for j in range(ny):
            x = float(i)
            y = float(j)
            velocity_field[i, j, 0] = -U0 * np.cos(k * x) * np.sin(k * y)  # u_x
            velocity_field[i, j, 1] = U0 * np.sin(k * x) * np.cos(k * y)   # u_y

    # Fill the data handling arrays (only inner region, not ghost layers)
    if len(rho_arr.shape) == 3:
        rho_arr[inner_slice[0], inner_slice[1], 0] = density_field
    else:
        rho_arr[inner_slice] = density_field

    # Velocity array
    u_arr[inner_slice[0], inner_slice[1], :] = velocity_field

    # Use compile_macroscopic_values_setter for proper initialization
    quantities_to_set = {
        'density': dh.cpu_arrays[rho.name],
        'velocity': dh.cpu_arrays[u.name]
    }

    setter_kernel_compiled = compile_macroscopic_values_setter(
        method,
        pdf_arr=dh.cpu_arrays[src.name],
        quantities_to_set=quantities_to_set,
        ghost_layers=0
    )

    # Run the setter to initialize PDFs
    setter_kernel_compiled(dh.cpu_arrays[src.name])

    print("Initialization complete.")
    print()

    # Compute initial kinetic energy
    def compute_kinetic_energy():
        dh.to_cpu(u.name)
        dh.to_cpu(rho.name)
        vel = dh.gather_array(u.name)
        density = dh.gather_array(rho.name)

        u_sq = vel[:, :, 0]**2 + vel[:, :, 1]**2
        # Handle both 2D and 3D density arrays
        if len(density.shape) == 3:
            rho_field = density[:, :, 0]
        else:
            rho_field = density
        KE = 0.5 * np.sum(rho_field * u_sq)
        return KE / (nx * ny)  # Normalized

    E0 = compute_kinetic_energy()
    print(f"Initial kinetic energy: {E0:.6e}")
    print()

    # Storage for energy history
    energy_history = []
    time_history = []
    energy_history.append(E0)
    time_history.append(0.0)

    # Periodic synchronization
    sync_pdfs = dh.synchronization_function([src.name])

    # Time loop
    for step in range(1, time_steps + 1):
        # Periodic BC
        sync_pdfs()

        # Collision and streaming (no wall boundaries - fully periodic)
        dh.run_kernel(stream_collide)

        # Swap arrays
        dh.swap(src.name, dst.name)

        # Output
        if step % output_interval == 0:
            E = compute_kinetic_energy()
            t = float(step)

            energy_history.append(E)
            time_history.append(t)

            # Extract viscosity from energy decay
            # For 2D Taylor-Green: E(t) = E0 * exp(-4 * nu * k^2 * t)
            # (factor of 4 instead of 2 for 2D vortex decay)
            # nu = -ln(E/E0) / (4 * k^2 * t)
            if E > 0 and E0 > 0:
                nu_extracted = -np.log(E / E0) / (4.0 * k**2 * t)
                error_percent = abs(nu_extracted - viscosity) / viscosity * 100.0

                print(f"Step {step:5d}: KE = {E:.6e}, "
                      f"nu_extracted = {nu_extracted:.6f}, "
                      f"error = {error_percent:.2f}%")

    print()

    # Final viscosity extraction
    E_final = energy_history[-1]
    t_final = time_history[-1]
    if E_final > 0 and E0 > 0:
        nu_extracted_final = -np.log(E_final / E0) / (4.0 * k**2 * t_final)
        error_percent_final = abs(nu_extracted_final - viscosity) / viscosity * 100.0
        print(f"Final viscosity extraction:")
        print(f"  Expected: {viscosity:.6f}")
        print(f"  Extracted: {nu_extracted_final:.6f}")
        print(f"  Error: {error_percent_final:.2f}%")
    else:
        nu_extracted_final = 0.0
        error_percent_final = 100.0
        print("WARNING: Energy decay is problematic")

    print()

    # Save results
    results = {
        'method': method_type,
        'time': np.array(time_history),
        'energy': np.array(energy_history),
        'nu_extracted': nu_extracted_final,
        'nu_expected': viscosity,
        'error_percent': error_percent_final,
        'E0': E0,
        'k': k,
        'params': {
            'nx': nx,
            'ny': ny,
            'U0': U0,
            'viscosity': viscosity,
            'k': k,
            'time_steps': time_steps,
            'omega': omega
        }
    }

    return results


def plot_three_way_comparison(bgk_results, elbm_results):
    """Plot BGK, ELBM, and analytical energy decay."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    t_bgk = bgk_results['time']
    E_bgk = bgk_results['energy']
    t_elbm = elbm_results['time']
    E_elbm = elbm_results['energy']

    # Analytical solution
    E0 = bgk_results['E0']
    nu = bgk_results['nu_expected']
    k = bgk_results['k']
    t_analytical = np.linspace(0, max(t_bgk[-1], t_elbm[-1]), 100)
    E_analytical = E0 * np.exp(-4.0 * nu * k**2 * t_analytical)

    # Plot 1: Energy decay
    ax1.semilogy(t_analytical, E_analytical, 'k--', linewidth=2, label='Analytical', alpha=0.8)
    ax1.semilogy(t_bgk, E_bgk, 'bo-', linewidth=1.5, markersize=6, label='BGK', alpha=0.8)
    ax1.semilogy(t_elbm, E_elbm, 'rs-', linewidth=1.5, markersize=6, label='ELBM', alpha=0.8)
    ax1.set_xlabel('Time steps', fontsize=12)
    ax1.set_ylabel('Kinetic Energy', fontsize=12)
    ax1.set_title('Taylor-Green Vortex: Energy Decay', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3, which='both')

    # Plot 2: Viscosity extraction over time
    def extract_viscosity_history(time, energy, E0, k):
        nu_hist = []
        for i in range(1, len(time)):
            t = time[i]
            E = energy[i]
            if E > 0 and E0 > 0 and t > 0:
                nu_ext = -np.log(E / E0) / (4.0 * k**2 * t)
                nu_hist.append(nu_ext)
            else:
                nu_hist.append(0)
        return np.array(nu_hist)

    nu_bgk_hist = extract_viscosity_history(t_bgk, E_bgk, E0, k)
    nu_elbm_hist = extract_viscosity_history(t_elbm, E_elbm, E0, k)

    ax2.plot(t_bgk[1:], nu_bgk_hist, 'b-', linewidth=1.5, label='BGK extracted', alpha=0.8)
    ax2.plot(t_elbm[1:], nu_elbm_hist, 'r-', linewidth=1.5, label='ELBM extracted', alpha=0.8)
    ax2.axhline(y=nu, color='k', linestyle='--', linewidth=2, label='Expected')
    ax2.set_xlabel('Time steps', fontsize=12)
    ax2.set_ylabel('Extracted viscosity', fontsize=12)
    ax2.set_title('Viscosity Extraction', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim([0.008, 0.012])

    # Add errors in title
    error_bgk = bgk_results['error_percent']
    error_elbm = elbm_results['error_percent']
    fig.suptitle(f'Taylor-Green Vortex | BGK error: {error_bgk:.2f}% | ELBM error: {error_elbm:.2f}%',
                 fontsize=13, fontweight='bold', y=1.02)

    plt.tight_layout()

    # Save figure - use Path(__file__) to find repo root
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent.parent
    output_dir = repo_root / 'figures' / 'lbmpy_validation'
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / 'taylor_green_3way_comparison.png'
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
    print("TAYLOR-GREEN VORTEX: BGK vs ELBM vs ANALYTICAL")
    print("="*70 + "\n")

    # Run BGK validation
    try:
        bgk_results = taylor_green_vortex_validation(method_type='BGK')
        print(f"✓ BGK validation complete (error: {bgk_results['error_percent']:.2f}%)")
    except Exception as e:
        print(f"✗ BGK validation FAILED with exception: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    # Run ELBM validation
    try:
        elbm_results = taylor_green_vortex_validation(method_type='ELBM')
        print(f"✓ ELBM validation complete (error: {elbm_results['error_percent']:.2f}%)")
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
    print(f"BGK viscosity error:  {bgk_results['error_percent']:.2f}%")
    print(f"ELBM viscosity error: {elbm_results['error_percent']:.2f}%")

    if bgk_results['error_percent'] < 1.0 and elbm_results['error_percent'] < 1.0:
        print("\n✓ All methods converged successfully!")
        sys.exit(0)
    else:
        print("\n⚠ Some methods did not fully converge")
        sys.exit(1)
