"""
2D Entropic Channel Flow Test
==============================
Simple 2D channel flow at low Re to test entropic stability.
"""

import numpy as np
import sympy as sp
import pystencils as ps
from lbmpy.boundaries import NoSlip, UBB
from lbmpy.creationfunctions import LBMConfig
from lbmpy.enums import Method, Stencil, ForceModel
from lbmpy.lbstep import LatticeBoltzmannStep
from lbmpy.relaxationrates import relaxation_rate_from_lattice_viscosity
from lbmpy.stencils import LBStencil
from pystencils import create_data_handling

print("="*60)
print("2D ENTROPIC CHANNEL FLOW TEST")
print("="*60)

# Physical parameters - very low Re for stability
width = 40
length = 160
target_re = 15
u_max = 0.05

print(f"Domain: {length} x {width}")
print(f"Target Re: {target_re}")
print(f"Max velocity: {u_max}")

# Calculate lattice viscosity
u_mean = u_max * 2/3  # For parabolic profile
lattice_viscosity = u_mean * width / target_re
omega = relaxation_rate_from_lattice_viscosity(lattice_viscosity)

print(f"Lattice viscosity: {lattice_viscosity:.6f}")
print(f"Relaxation rate: {omega:.6f}")

# Body force to drive flow (very small)
force_density = 1e-6

print(f"Force density: {force_density:.2e}")
print()

# Configure entropic LBM
try:
    print("Creating LBM configuration with entropic=True...")
    omega_free = sp.Symbol('omega_free')

    lbm_config = LBMConfig(
        stencil=LBStencil(Stencil.D2Q9),
        compressible=True,
        method=Method.TRT,
        relaxation_rates=[omega, omega_free],
        entropic=True,
        entropic_newton_iterations=3,
        zero_centered=False,
        force_model=ForceModel.LUO,  # LUO compatible with entropic
        force=(force_density, 0.0)
    )
    print("✓ LBM config created successfully")
except Exception as e:
    print(f"✗ Failed to create LBM config: {e}")
    import sys
    sys.exit(1)

# Create data handling - 2D domain, periodic in x
try:
    print("Creating data handling...")
    dh = create_data_handling(
        domain_size=(length, width),
        periodicity=(True, False),  # Periodic in x, walls in y
        parallel=False
    )
    print(f"✓ Data handling created: {dh.dim}D domain")

    print("Creating Lattice Boltzmann step...")
    sc = LatticeBoltzmannStep(
        data_handling=dh,
        lbm_config=lbm_config,
        name="entropic_channel"
    )
    print("✓ LB step created successfully")
except Exception as e:
    print(f"✗ Failed to create simulation: {e}")
    import traceback
    traceback.print_exc()
    import sys
    sys.exit(1)

# Set up boundaries - no-slip walls at top and bottom
try:
    print("Setting up boundaries...")
    bh = sc.boundary_handling
    wall = NoSlip()

    # Bottom wall (y=0)
    bh.set_boundary(wall, ps.make_slice[:, :1])

    # Top wall (y=width-1)
    bh.set_boundary(wall, ps.make_slice[:, -1:])

    print("✓ Boundaries configured (no-slip top and bottom, periodic x)")
except Exception as e:
    print(f"✗ Failed to set boundaries: {e}")
    import traceback
    traceback.print_exc()
    import sys
    sys.exit(1)

# Run simulation
print()
print("Starting simulation...")
print("-"*60)

try:
    last_max_vel = 0
    for i in range(200):  # Run 200 batches of 100 steps = 20,000 steps
        sc.run(100)

        # Extract velocity array properly
        vel_arr = sc.velocity[:, :]
        u_vel = vel_arr[:, :, 0]  # x-component
        max_vel = np.max(np.abs(u_vel))

        if (i + 1) % 10 == 0:
            # Average profile over x (should be uniform due to periodicity)
            u_profile = np.mean(u_vel, axis=0)
            center_vel = u_profile[width//2]
            print(f"Batch {i+1:3d}: steps = {sc.time_steps_run:6d}, "
                  f"max_u = {max_vel:.6e}, center_u = {center_vel:.6e}")

            # Check convergence (after initial development)
            if i > 20:
                change = abs(max_vel - last_max_vel) / (last_max_vel + 1e-12)
                if change < 1e-6:
                    print(f"✓ Converged at step {sc.time_steps_run} (change: {change:.2e})")
                    break
            last_max_vel = max_vel

        # Check for NaN
        if np.isnan(max_vel):
            print("✗ NaN detected - simulation unstable!")
            import sys
            sys.exit(1)

    print()
    print("="*60)
    print("✓ ENTROPIC 2D SIMULATION COMPLETED SUCCESSFULLY")
    print("="*60)
    print(f"Total time steps: {sc.time_steps_run}")

    # Extract final velocity profile
    vel_arr = sc.velocity[:, :]
    u_vel = vel_arr[:, :, 0]
    u_profile = np.mean(u_vel, axis=0)

    print(f"Final max velocity: {np.max(u_profile):.6e}")
    print(f"Final center velocity: {u_profile[width//2]:.6e}")

    # Check if profile is parabolic-like
    if u_profile[0] < 1e-8 and u_profile[-1] < 1e-8 and u_profile[width//2] > 1e-8:
        print("✓ Velocity profile looks parabolic (zero at walls, max at center)")

    import sys
    sys.exit(0)

except Exception as e:
    print(f"✗ Simulation failed: {e}")
    import traceback
    traceback.print_exc()
    import sys
    sys.exit(1)
