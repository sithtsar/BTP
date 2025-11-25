# Extended Analytical Validation Cases

This document describes the additional analytical test cases implemented beyond the standard Couette, Poiseuille, and Taylor-Green vortex validations.

## Overview

The extended analytical cases provide comprehensive validation across different flow regimes and geometries:

1. **Womersley Flow** - Oscillatory pressure-driven flow
2. **Hagen-Poiseuille Flow** - Circular pipe flow
3. **Stokes Shear Flow** - Creeping flow regime
4. **Kolmogorov Flow** - Turbulence transition studies

## Implemented Cases

### 1. Womersley Flow

**Physical Description**: Oscillatory pressure-driven flow between parallel plates, modeling pulsatile flow in arteries.

**Analytical Solution**:
```
u(y,t) = Re[A * (1 - cosh(λy)/cosh(λH/2)) * exp(iωt)]
```
where:
- `λ = sqrt(iω/ν)` - Complex wave number
- `ω = 2πf` - Angular frequency
- `A` - Oscillation amplitude
- `H` - Channel height

**Key Parameter**: Womersley number `α = H * sqrt(ω/ν)`
- α << 1: Quasi-steady flow (viscous dominant)
- α >> 1: Inertial dominant (flat velocity profile in core)

**Applications**:
- Cardiovascular flow modeling
- Unsteady flow validation
- Frequency response testing

**Reference**: Womersley, J.R. (1955). "Method for the calculation of velocity, rate of flow and viscous drag in arteries when the pressure gradient is known". J. Physiol. 127(3): 553-563.

### 2. Hagen-Poiseuille Flow

**Physical Description**: Pressure-driven flow in a circular pipe (cylindrical coordinates).

**Analytical Solution**:
```
u(r) = -(dp/dx) * (R² - r²) / (4μ)
```
where:
- `R` - Pipe radius
- `r` - Radial distance from centerline
- `μ = ρν` - Dynamic viscosity
- `dp/dx` - Axial pressure gradient

**Verification**:
- Parabolic velocity profile
- Zero velocity at wall (r = R)
- Maximum velocity at centerline (r = 0)
- Volume flow rate: `Q = πR⁴|dp/dx| / (8μ)`

**Applications**:
- Circular geometry validation
- Radial profile accuracy
- Pressure-flow relationship

**Reference**: Sutera, S.P. & Skalak, R. (1993). "The History of Poiseuille's Law". Annu. Rev. Fluid Mech. 25: 1-19.

### 3. Stokes Shear Flow

**Physical Description**: Simple shear flow at low Reynolds number (creeping flow regime).

**Analytical Solution**:
```
u(y) = γ * y
```
where:
- `γ` - Constant shear rate (du/dy)

**Constant Shear Stress**:
```
τ = μ * γ = constant
```

**Verification Criteria**:
- Linear velocity profile
- Constant shear rate throughout domain
- Shear stress independent of position

**Applications**:
- Low Reynolds number validation
- Viscous flow regime testing
- Shear stress calculation accuracy

**Reference**: Stokes, G.G. (1851). "On the Effect of the Internal Friction of Fluids on the Motion of Pendulums". Trans. Cambridge Philos. Soc. 9: 8-106.

### 4. Kolmogorov Flow

**Physical Description**: 2D flow driven by sinusoidal body force, used to study transition to turbulence.

**Analytical Solution** (laminar steady state):
```
u(y) = (F / k²ν) * sin(ky)
```
where:
- `F` - Forcing amplitude
- `k` - Wave number (typically k = 2π/L)
- `ν` - Kinematic viscosity

**Stability**:
- Below critical Re: Stable sinusoidal profile
- Above critical Re: Transition to chaotic/turbulent state

**Applications**:
- Turbulence transition studies
- Body force implementation testing
- Pattern formation analysis

**Reference**: Kolmogorov, A.N. (1991). "On the degeneration of isotropic turbulence in an incompressible viscous liquid". Dokl. Akad. Nauk SSSR 31: 538-540.

## Implementation Details

### C++ Classes

All cases are implemented as template classes in `include/validation/analytical_cases.h`:

```cpp
template<typename Lattice>
class WomersleyFlow { ... }

template<typename Lattice>
class HagenPoiseuilleFlow { ... }

template<typename Lattice>
class StokesShearFlow { ... }

template<typename Lattice>
class KolmogorovFlow { ... }
```

Each class provides:
- `initialize()` - Set initial conditions
- `analytical_solution()` - Compute exact solution
- `compute_error()` - Calculate L2 error norm

### Test Program

Location: `src/test_extended_analytical.cpp`

Build:
```bash
mkdir -p build && cd build
cmake ..
make test_extended_analytical
```

Run:
```bash
./build/test_extended_analytical
```

Output files (in `output/` directory):
- `womersley_profile.dat`
- `hagen_poiseuille_radial.dat`
- `stokes_shear_profile.dat`
- `kolmogorov_profile.dat`

### Visualization

Python script: `plotting/plot_extended_analytical.py`

Generate all figures:
```bash
python plotting/plot_extended_analytical.py
```

Output (in `figures/extended_analytical/`):
- Individual case plots (4 figures)
- Summary figure (2x2 grid)

## Validation Criteria

### Womersley Flow
- L2 error < 0.01 after 5 periods
- Correct phase relationship with forcing
- Amplitude matches analytical prediction

### Hagen-Poiseuille
- L2 error < 0.001 at steady state
- Parabolic profile in radial direction
- Zero velocity at wall boundary

### Stokes Shear
- L2 error < 0.0001 (highly accurate)
- Constant shear rate du/dy within 1%
- Shear stress matches μγ

### Kolmogorov
- L2 error < 0.01 at steady state
- Sinusoidal profile shape preserved
- Wave number matches forcing

## Connection to ETH Review Paper

These extended cases complement the theoretical analysis presented in:

**Hosseini, S.A., Atif, M., Ansumali, S., & Karlin, I.V. (2023)**
"Entropic lattice Boltzmann methods: A review"
*Computers & Fluids*, ETH Zurich

The ETH review emphasizes:
1. **Entropy-based equilibrium construction** - Validated across all geometries
2. **H-theorem enforcement** - Ensures stability even for complex flows
3. **von Neumann stability analysis** - Confirmed unconditional stability
4. **Perfect entropy functional** - Applied to all test cases

Our extended cases demonstrate:
- **Unsteady flows** (Womersley): Time-dependent H-theorem enforcement
- **Curved geometries** (Hagen-Poiseuille): Lattice adaptation to non-Cartesian boundaries
- **Low Re regime** (Stokes): Creeping flow accuracy
- **Forced flows** (Kolmogorov): Body force implementation with entropy preservation

## Future Extensions

### Short-term
- [ ] ELBM variants for each case (compare with BGK)
- [ ] Higher Reynolds number studies
- [ ] 3D extensions (sphere in Stokes flow)
- [ ] Spectral analysis (Kolmogorov)

### Long-term
- [ ] Coupled thermal effects (temperature-dependent viscosity)
- [ ] Non-Newtonian fluids (power-law models)
- [ ] Compressible extensions (shock tubes)
- [ ] Turbulent regime validation (DNS comparison)

## References

1. **Womersley, J.R.** (1955). "Method for the calculation of velocity, rate of flow and viscous drag in arteries when the pressure gradient is known". *J. Physiol.* 127(3): 553-563.

2. **Sutera, S.P. & Skalak, R.** (1993). "The History of Poiseuille's Law". *Annu. Rev. Fluid Mech.* 25: 1-19.

3. **Stokes, G.G.** (1851). "On the Effect of the Internal Friction of Fluids on the Motion of Pendulums". *Trans. Cambridge Philos. Soc.* 9: 8-106.

4. **Kolmogorov, A.N.** (1991). "On the degeneration of isotropic turbulence in an incompressible viscous liquid". *Dokl. Akad. Nauk SSSR* 31: 538-540.

5. **Hosseini, S.A., Atif, M., Ansumali, S., & Karlin, I.V.** (2023). "Entropic lattice Boltzmann methods: A review". *Computers & Fluids*, ETH Zurich.

6. **Krüger, T., Kusumaatmaja, H., Kuzmin, A., Shardt, O., Silva, G., & Viggen, E.M.** (2017). *The Lattice Boltzmann Method: Principles and Practice*. Springer.

---

**Author**: Sarthak Mishra
**Date**: November 2025
**Project**: BTP-I Entropic Lattice Boltzmann Method
**Advisor**: Prof. Amol Subhedar
