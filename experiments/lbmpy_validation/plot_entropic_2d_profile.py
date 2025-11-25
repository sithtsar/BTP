"""
Visualize 2D Entropic Channel Flow Profile
==========================================
"""

import numpy as np
import sympy as sp
import pystencils as ps
import matplotlib.pyplot as plt
from pathlib import Path
from lbmpy.boundaries import NoSlip
from lbmpy.creationfunctions import LBMConfig
from lbmpy.enums import Method, Stencil, ForceModel
from lbmpy.lbstep import LatticeBoltzmannStep
from lbmpy.relaxationrates import relaxation_rate_from_lattice_viscosity
from lbmpy.stencils import LBStencil
from pystencils import create_data_handling

# Same parameters as test
width = 40
length = 160
target_re = 15
u_max = 0.05
u_mean = u_max * 2/3
lattice_viscosity = u_mean * width / target_re
omega = relaxation_rate_from_lattice_viscosity(lattice_viscosity)
force_density = 1e-6

# Setup simulation
omega_free = sp.Symbol('omega_free')
lbm_config = LBMConfig(
    stencil=LBStencil(Stencil.D2Q9),
    compressible=True,
    method=Method.TRT,
    relaxation_rates=[omega, omega_free],
    entropic=True,
    entropic_newton_iterations=3,
    zero_centered=False,
    force_model=ForceModel.LUO,
    force=(force_density, 0.0)
)

dh = create_data_handling(
    domain_size=(length, width),
    periodicity=(True, False),
    parallel=False
)

sc = LatticeBoltzmannStep(
    data_handling=dh,
    lbm_config=lbm_config,
    name="entropic_channel"
)

# Boundaries
bh = sc.boundary_handling
wall = NoSlip()
bh.set_boundary(wall, ps.make_slice[:, :1])
bh.set_boundary(wall, ps.make_slice[:, -1:])

# Run to steady state
print("Running simulation to steady state...")
for i in range(200):
    sc.run(100)
    if (i + 1) % 50 == 0:
        print(f"  Step {sc.time_steps_run}")

# Extract velocity profile
vel_arr = sc.velocity[:, :]
u_vel = vel_arr[:, :, 0]
u_profile = np.mean(u_vel, axis=0)

# Analytical solution for Poiseuille flow
# u(y) = (F/(2*rho*nu)) * y * (H - y)
H = width - 2  # Channel height (excluding walls)
y = np.arange(width)
rho = 1.0

# Adjust for wall positions
y_inner = np.arange(1, width - 1)
y_centered = y_inner - width/2

# Analytical profile
u_analytical = np.zeros(width)
for i in range(1, width - 1):
    y_pos = i - 1  # Position from bottom wall
    u_analytical[i] = (force_density / (2.0 * rho * lattice_viscosity)) * y_pos * (H - y_pos)

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Left: Velocity profile
ax1.plot(u_profile * 1000, y, 'r-', linewidth=2, label='ELBM', alpha=0.8)
ax1.plot(u_analytical * 1000, y, 'k--', linewidth=2, label='Analytical', alpha=0.6)
ax1.set_xlabel('Velocity u (×10⁻³)', fontsize=12)
ax1.set_ylabel('y position', fontsize=12)
ax1.set_title('2D Entropic Channel Flow Profile', fontsize=14, fontweight='bold')
ax1.legend(fontsize=11)
ax1.grid(True, alpha=0.3)
ax1.axhline(y=0, color='gray', linestyle=':', linewidth=1, alpha=0.5)
ax1.axhline(y=width-1, color='gray', linestyle=':', linewidth=1, alpha=0.5)

# Right: Error
error = np.abs(u_profile - u_analytical)
error_inner = error[1:-1]
y_inner = y[1:-1]

ax2.semilogy(error_inner, y_inner, 'b-', linewidth=2, alpha=0.8)
ax2.set_xlabel('Absolute Error', fontsize=12)
ax2.set_ylabel('y position', fontsize=12)
ax2.set_title('Error Distribution', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3, which='both')

# Compute L2 error
l2_error = np.sqrt(np.mean(error_inner**2))
fig.suptitle(f'2D Entropic Channel Flow | Re={target_re} | L2 Error: {l2_error:.2e}',
             fontsize=13, fontweight='bold', y=0.98)

plt.tight_layout()

# Save using Path(__file__) to find repo root
script_dir = Path(__file__).resolve().parent
repo_root = script_dir.parent.parent
output_dir = repo_root / 'figures' / 'lbmpy_validation'
output_dir.mkdir(parents=True, exist_ok=True)
output_file = output_dir / 'entropic_2d_channel_profile.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\n✓ Saved: {output_file}")

# Print summary
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
print(f"Domain: {length} x {width}")
print(f"Reynolds number: {target_re}")
print(f"Max velocity: {np.max(u_profile):.6e}")
print(f"L2 error: {l2_error:.6e}")
print(f"Wall velocities: u(0)={u_profile[0]:.2e}, u(H)={u_profile[-1]:.2e}")
print("="*60)
