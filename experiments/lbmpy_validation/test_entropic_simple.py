"""
Simple Entropic LBM Test
========================
Test basic entropic method functionality with a simple pipe flow setup.
"""

import numpy as np
from lbmpy.boundaries import FixedDensity, NoSlip
from lbmpy.creationfunctions import LBMConfig
from lbmpy.enums import Method, Stencil
from lbmpy.lbstep import LatticeBoltzmannStep
from lbmpy.relaxationrates import relaxation_rate_from_lattice_viscosity
from lbmpy.stencils import LBStencil
from pystencils import create_data_handling
from pystencils.slicing import slice_from_direction

print("="*60)
print("SIMPLE ENTROPIC LBM TEST")
print("="*60)

# Physical parameters
diameter = 20
length = 16 * diameter
u_max = 0.07
re = 150

print(f"Domain: {length} x {diameter} x {diameter}")
print(f"Target Re: {re}")
print(f"Max velocity: {u_max}")

# Calculate lattice viscosity and relaxation rate
u_mean = u_max / 2
lattice_viscosity = u_mean * diameter / re
omega = relaxation_rate_from_lattice_viscosity(lattice_viscosity)

print(f"Lattice viscosity: {lattice_viscosity:.6f}")
print(f"Relaxation rate: {omega:.6f}")
print()

# Configure entropic LBM
try:
    print("Creating LBM configuration with entropic=True...")
    # Use TRT with two relaxation rates for entropic method
    import sympy as sp
    omega_free = sp.Symbol('omega_free')  # Free parameter determined by entropy condition

    lbm_config = LBMConfig(
        stencil=LBStencil(Stencil.D3Q19),
        compressible=True,
        method=Method.TRT,  # Two-Relaxation-Time required for entropic
        relaxation_rates=[omega, omega_free],  # Two rates: shear and bulk
        entropic=True,  # Enable entropic method
        entropic_newton_iterations=3,  # Newton iterations for entropy optimization
        zero_centered=False  # Required for entropic methods
    )
    print("✓ LBM config created successfully")
except Exception as e:
    print(f"✗ Failed to create LBM config: {e}")
    import sys
    sys.exit(1)

# Create data handling and scenario
try:
    print("Creating data handling...")
    dh = create_data_handling(domain_size=(length, diameter, diameter), parallel=False)
    print(f"✓ Data handling created: {dh.dim}D domain")

    print("Creating Lattice Boltzmann step...")
    sc = LatticeBoltzmannStep(
        data_handling=dh,
        lbm_config=lbm_config,
        name="entropic_pipe"
    )
    print("✓ LB step created successfully")
except Exception as e:
    print(f"✗ Failed to create simulation: {e}")
    import traceback
    traceback.print_exc()
    import sys
    sys.exit(1)

# Set up boundaries
def pipe_geometry(x, y, z):
    """Define circular pipe geometry"""
    radius = diameter / 2
    y_mid = diameter / 2
    z_mid = diameter / 2
    return (y - y_mid)**2 + (z - z_mid)**2 > radius**2

try:
    print("Setting up boundaries...")
    bh = sc.boundary_handling
    wall = NoSlip()
    bh.set_boundary(wall, mask_callback=pipe_geometry)

    # Add outflow boundary
    outflow = FixedDensity(1.0)
    bh.set_boundary(outflow, slice_from_direction('E', 3))
    print("✓ Boundaries configured")
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
    for i in range(10):  # Run 10 batches of 100 steps
        sc.run(100)
        # Extract velocity array properly
        vel_arr = sc.velocity[:, :, :]
        max_vel = np.max(np.abs(vel_arr))
        print(f"Batch {i+1:2d}: Time steps = {sc.time_steps_run:5d}, Max velocity = {max_vel:.6f}")

        # Check for NaN
        if np.isnan(max_vel):
            print("✗ NaN detected - simulation unstable!")
            break
    else:
        print()
        print("="*60)
        print("✓ ENTROPIC SIMULATION COMPLETED SUCCESSFULLY")
        print("="*60)
        print(f"Total time steps: {sc.time_steps_run}")
        vel_arr = sc.velocity[:, :, :]
        print(f"Final max velocity: {np.max(np.abs(vel_arr)):.6f}")
        import sys
        sys.exit(0)

except Exception as e:
    print(f"✗ Simulation failed: {e}")
    import traceback
    traceback.print_exc()
    import sys
    sys.exit(1)
