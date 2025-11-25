"""
Verify the analytical solution for plane Poiseuille flow.
"""
import numpy as np
import matplotlib.pyplot as plt

# Parameters
H = 30  # Channel height (interior, excluding walls)
rho = 1.0
nu = 0.1  # kinematic viscosity
F = 1e-4  # body force

# Analytical solution: u(y) = F/(2*rho*nu) * y * (H - y)
y = np.linspace(0, H, 100)
u_analytical = (F / (2.0 * rho * nu)) * y * (H - y)

# Maximum velocity should be at y = H/2
y_max = H / 2
u_max_analytical = (F / (2.0 * rho * nu)) * y_max * (H - y_max)
u_max_formula = F * H**2 / (8 * rho * nu)  # Alternative formula

print("Analytical Solution Verification")
print("=" * 50)
print(f"Channel height: H = {H}")
print(f"Body force: F = {F}")
print(f"Density: ρ = {rho}")
print(f"Kinematic viscosity: ν = {nu}")
print(f"Dynamic viscosity: μ = ρ*ν = {rho * nu}")
print()
print(f"Maximum velocity (from profile): {u_max_analytical:.6e}")
print(f"Maximum velocity (from formula):  {u_max_formula:.6e}")
print(f"Location of u_max: y = {y_max} (should be H/2 = {H/2})")
print()

# Check that maximum is at center
u_at_center = (F / (2.0 * rho * nu)) * (H/2) * (H - H/2)
print(f"Velocity at center: {u_at_center:.6e}")
print(f"Matches u_max: {np.isclose(u_at_center, u_max_analytical)}")
print()

# Reynolds number
Re = u_max_analytical * H / nu
print(f"Reynolds number: Re = u_max * H / ν = {Re:.2f}")
print()

# Plot
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(u_analytical * 1000, y, 'k-', linewidth=2, label='Analytical Solution')
ax.axhline(y=H/2, color='r', linestyle='--', alpha=0.5, label=f'Centerline (y={H/2})')
ax.axvline(x=u_max_analytical * 1000, color='r', linestyle=':', alpha=0.5, label=f'u_max={u_max_analytical*1000:.3f}×10⁻³')
ax.set_xlabel('Velocity u (×10⁻³)', fontsize=12)
ax.set_ylabel('y position', fontsize=12)
ax.set_title(f'Plane Poiseuille Flow: u(y) = F/(2ρν) · y(H-y)\nRe={Re:.1f}', fontsize=13)
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('/tmp/verify_analytical_poiseuille.png', dpi=150)
print(f"✓ Saved: /tmp/verify_analytical_poiseuille.png")

# Verify the formula for different Re values
print("\nRe scaling verification:")
print("-" * 50)
re_values = [10, 20, 30, 40, 50]
for target_re in re_values:
    # For given Re, what force is needed?
    # Re = u_max * H / nu
    # u_max = F * H^2 / (8 * rho * nu)
    # So: Re = [F * H^2 / (8 * rho * nu)] * H / nu = F * H^3 / (8 * rho * nu^2)
    # Therefore: F = 8 * rho * nu^2 * Re / H^3
    F_needed = 8 * rho * nu**2 * target_re / H**3
    u_max_expected = F_needed * H**2 / (8 * rho * nu)
    re_check = u_max_expected * H / nu
    print(f"Re={target_re:2d}: F={F_needed:.6e}, u_max={u_max_expected:.6e}, Re_check={re_check:.2f}")
