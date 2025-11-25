# Multi-Phase Fluid Simulation using lbmpy

## Overview

This directory contains implementations of **Shan-Chen two-phase** (liquid-gas) lattice Boltzmann simulations using lbmpy with MRT collision operators. The Shan-Chen model uses pseudo-potential interactions to achieve phase separation and simulate surface tension effects.

## Theory: Shan-Chen Pseudo-Potential Model

### Basic Concept

The Shan-Chen model introduces phase separation through a **non-local, density-dependent interaction force**:

```
F = -g * psi(ρ) * Σ_i [ w_i * psi(ρ_i) * c_i ]
```

where:
- `g < 0`: Interaction strength (negative for attraction → phase separation)
- `psi(ρ) = ρ₀ * (1 - exp(-ρ/ρ₀))`: Pseudo-potential function
- `w_i`: Lattice weights from D2Q9 stencil
- `ρ_i`: Density at neighboring node i
- `c_i`: Lattice velocity direction i

### Phase Separation Mechanism

1. **High density regions** (liquid): Large `psi(ρ)` → strong attractive force → densification
2. **Low density regions** (gas): Small `psi(ρ)` → weak force → rarefaction
3. **Interface**: Sharp transition where forces balance
4. **Surface tension**: Emerges naturally from force imbalance at curved interfaces

### Laplace Pressure Law

For a static droplet with radius R:

```
ΔP = P_inside - P_outside = σ/R
```

where σ is the surface tension (implicitly determined by `g`).

This relationship validates:
- ΔP should be positive (pressure higher inside droplet)
- ΔP ∝ 1/R (linear with curvature)
- Different `g` values → different surface tensions

### Collision Operators

**MRT (Multiple Relaxation Time)**:
- Relaxes different moments at different rates
- Better stability than BGK/SRT
- Reduces spurious currents near interfaces
- Default choice for multi-phase simulations

**SRT (Single Relaxation Time / BGK)**:
- Simpler, single relaxation parameter ω
- Faster computation
- Higher spurious currents
- Suitable for testing/prototyping

**Entropic LBM** (Future):
- Unconditionally stable via H-theorem
- Requires custom implementation for multi-phase
- Research-grade extension

## Directory Structure

```
multiphase_validation/
├── README.md                          # This file
├── multiphase_framework.py            # Base class for Shan-Chen LBM
├── static_droplet_mrt.py              # Phase 1: Static droplet validation
├── droplet_collision_mrt.py           # Phase 2: Dynamic collision
├── rayleigh_taylor_mrt.py             # Phase 3: Instability test
├── spinodal_decomposition_mrt.py      # Phase 4: Phase separation
├── plot_multiphase_results.py         # Visualization utilities
└── compare_collision_operators.py     # SRT vs MRT comparison
```

## Implementation Files

### 1. `multiphase_framework.py`

**Base class**: `ShanChenMultiPhase`

Provides core functionality:
- Shan-Chen force calculation (symbolic with sympy)
- MRT/SRT collision kernel creation (via lbmpy)
- Streaming kernel with density output
- Periodic boundary conditions
- Data handling (CPU/GPU support via pystencils)
- Visualization methods (density field, velocity field)

**Key Methods**:
- `initialize_density_field()`: **Abstract** - implement in subclasses
- `initialize()`: Set PDFs to equilibrium based on density
- `run(time_steps)`: Execute LBM simulation loop
- `get_density_field()`: Extract current density as numpy array
- `get_velocity_field()`: Compute velocities from PDFs
- `compute_spurious_currents()`: Measure numerical artifacts
- `plot_density_field()`: Visualize density distribution
- `plot_velocity_field()`: Visualize velocity vectors
- `save_state()`: Save simulation state to .npz file

### 2. `static_droplet_mrt.py`

**Purpose**: Validate Shan-Chen model via static droplet test

**Test Objectives**:
1. ✓ Stable droplet interface (minimize spurious currents)
2. ✓ Measure Laplace pressure: ΔP = σ/R
3. ✓ Verify phase separation (liquid core, gas outside)
4. ✓ Confirm sharp interface (3-5 lattice units transition)

**Parameters** (default):
```python
N = 64                 # Domain size (64×64)
droplet_radius = 15.0  # Initial radius
rho_liquid = 2.1       # Liquid phase density
rho_gas = 0.15         # Gas phase density
omega = 1.0            # Relaxation parameter
g_interaction = -4.7   # Shan-Chen interaction strength
time_steps = 10000     # Equilibration time
```

**Expected Results**:
- Spurious currents: < 1e-3 (for stable droplet)
  - **Note**: Achieving this requires careful parameter tuning
  - Higher values (0.01-0.1) are common but indicate numerical artifacts
- Positive Laplace pressure: ΔP > 0
- Circular droplet shape preserved
- Density ratio preserved: ρ_liquid / ρ_gas ≈ 14

**Output**:
- `figures/multiphase/static_droplet_density.png`: Density field heatmap
- `figures/multiphase/static_droplet_radial_profile.png`: Azimuthal average
- `figures/multiphase/static_droplet_velocity.png`: Spurious current visualization
- `figures/multiphase/static_droplet_pressure_analysis.png`: Laplace law validation
- `output/multiphase/static_droplet_final.npz`: Final state

**Usage**:
```bash
uv run python experiments/multiphase_validation/static_droplet_mrt.py
```

### 3. `droplet_collision_mrt.py`

**Purpose**: Test dynamic interactions between droplets

**Scenario**: Two droplets with initial velocities collide and coalesce

**Validation Criteria**:
- Momentum conservation: < 1e-6 deviation
- Successful coalescence (single merged droplet)
- No NaN or instabilities
- Interface dynamics (capillary waves after merger)

**Status**: ⚠️ Pending implementation

### 4. `rayleigh_taylor_mrt.py`

**Purpose**: Validate gravity-driven instability

**Scenario**: Heavy fluid initially above light fluid → instability growth

**Validation Criteria**:
- Instability develops (mushroom structures)
- Growth rate matches linear stability theory
- Late-stage turbulent mixing observed

**Status**: ⚠️ Pending implementation

### 5. `spinodal_decomposition_mrt.py`

**Purpose**: Test spontaneous phase separation

**Scenario**: Uniform density with random fluctuations → domain formation

**Validation Criteria**:
- Domains emerge from noise
- Power-law coarsening dynamics
- Equilibrium: two separated phases

**Status**: ⚠️ Pending implementation

## Running the Simulations

### Prerequisites

Ensure packages are installed (already in `pyproject.toml`):
```bash
uv pip install lbmpy pystencils sympy numpy matplotlib
```

### Quick Start

**Run static droplet test**:
```bash
uv run python experiments/multiphase_validation/static_droplet_mrt.py
```

**Check results**:
```bash
ls figures/multiphase/
# static_droplet_density.png
# static_droplet_radial_profile.png
# static_droplet_velocity.png
# static_droplet_pressure_analysis.png
```

### Parameter Tuning Guide

#### Reducing Spurious Currents

Spurious currents are parasitic velocities near interfaces. To minimize:

1. **Increase equilibration time**:
   ```python
   time_steps = 20000  # Double the default
   ```

2. **Adjust interaction strength** (trade-off with interface width):
   ```python
   g_interaction = -4.5  # Less negative → weaker force → lower currents
   ```

3. **Use smaller density ratio**:
   ```python
   rho_liquid = 1.5
   rho_gas = 0.3
   # Ratio: 5 instead of 14
   ```

4. **Refine MRT relaxation rates** (advanced):
   ```python
   # Custom relaxation rates for moment space
   # Requires deeper lbmpy knowledge
   ```

#### Surface Tension Control

Surface tension σ is implicitly set by `g_interaction`:

- **Higher |g|** → Stronger attraction → Higher σ → Sharper interface
- **Lower |g|** → Weaker attraction → Lower σ → Diffuse interface

Typical range: `-5.0 < g < -4.0`

To measure σ from simulations:
```python
sigma_implied = delta_P * radius
```

Run multiple simulations with different radii and fit: ΔP vs 1/R

#### Domain Size and Resolution

- **Minimum**: `N = 4 * droplet_radius` (avoid periodic boundary effects)
- **Recommended**: `N = 64` for radius ~15
- **High-resolution**: `N = 128` or `N = 256` (slower but more accurate)

Interface width typically ~3-5 lattice units, so ensure:
```
droplet_radius >> 5 (for well-resolved droplet)
N - droplet_radius >> 5 (to avoid domain boundary)
```

## Common Issues & Solutions

### Issue 1: High Spurious Currents (> 0.01)

**Symptoms**:
- Validation fails with "NEEDS REVIEW"
- Velocity field shows circulation patterns near interface

**Solutions**:
1. Run longer (20k-50k steps instead of 10k)
2. Reduce density ratio (e.g., ρ_liquid/ρ_gas = 5 instead of 14)
3. Use `g = -4.5` instead of `-4.7`
4. Check MRT is enabled: `use_mrt=True`

### Issue 2: Droplet Evaporates or Explodes

**Symptoms**:
- Density field becomes uniform
- NaN values appear

**Solutions**:
1. Check `g < 0` (negative for attraction)
2. Reduce `|g|` if too strong (try `-3.0` to `-4.5` range)
3. Ensure `omega` is reasonable (`0.5 < omega < 1.5`)
4. Verify initial density ratio not extreme (< 20)

### Issue 3: Non-Circular Droplet

**Symptoms**:
- Droplet elongates or develops lobes

**Solutions**:
1. Increase domain size (`N` too small → periodic BC artifacts)
2. Run longer (interface still relaxing)
3. Check initialization (ensure perfect circle in `initialize_density_field()`)

### Issue 4: Import Errors

**Error**: `ImportError: cannot import name 'Stencil' from 'pystencils.stencil'`

**Solution**:
```python
# WRONG:
from pystencils.stencil import Stencil

# CORRECT:
from lbmpy.enums import Stencil
```

Use imports from `lbmpy.enums` and `lbmpy.stencils`, not `pystencils.stencil`.

## Validation Criteria Summary

### Static Droplet

| Metric | Threshold | Status (Current) |
|--------|-----------|-----------------|
| Spurious currents | < 1e-3 | ✗ 0.105 (needs tuning) |
| Laplace pressure ΔP | > 0 | ✓ 0.495 |
| Interface sharpness | 3-5 lattice units | ✓ (visual) |
| Circular shape | Preserved | ✓ (visual) |

**Interpretation**:
- ✓ **Phase separation works**: Droplet forms and maintains interface
- ✓ **Surface tension present**: Positive Laplace pressure observed
- ⚠️ **Spurious currents high**: Requires parameter tuning for quantitative accuracy
- ✓ **Qualitative success**: Model demonstrates correct physics

### Next Steps for Validation

1. **Parameter sweep**: Test g ∈ [-5.0, -4.0], omega ∈ [0.8, 1.2]
2. **Longer runs**: 50k steps to ensure equilibration
3. **Multi-radius test**: Run R = [10, 15, 20, 25] and fit ΔP vs 1/R
4. **Compare with lbmpy examples**: Check `tests/shan_chen/test_shan_chen_two_phase.py`

## Performance Notes

### CPU vs GPU

- **CPU** (default): Easier debugging, slower for large domains
- **GPU**: 10-100× faster, requires CUDA/ROCm, change `target='gpu'`

### Typical Runtimes (CPU, M1 Mac)

| Configuration | Time Steps | Wall Time |
|--------------|------------|-----------|
| N=64, 10k steps | 10,000 | ~30 seconds |
| N=64, 50k steps | 50,000 | ~2.5 minutes |
| N=128, 10k steps | 10,000 | ~3 minutes |
| N=256, 10k steps | 10,000 | ~25 minutes |

### Optimization Tips

1. Use GPU target for large domains (N > 128)
2. Reduce `print_interval` to avoid I/O overhead
3. Save state only at end (not every N steps)
4. Use `show=False` in plotting to avoid display overhead

## Future Extensions

### Phase 2-4: Additional Test Cases

Once static droplet is fully tuned, implement:

1. **Droplet collision**: Test dynamic interactions, momentum conservation
2. **Rayleigh-Taylor**: Validate instability growth rates
3. **Spinodal decomposition**: Test phase separation kinetics

### Entropic Multi-Phase

Extend to entropic collision operator:
- Define H-function for two-phase system
- Implement alpha-solver with Shan-Chen forces
- Test unconditional stability at high density ratios

### C++ Port

For production performance:
- Port framework to C++ (follow `include/solvers/` structure)
- Create `multiphase_bgk_solver.h` and `multiphase_elbm_solver.h`
- Extend `fluid_state.h` for multi-component systems
- Add `test_multiphase.cpp` executable

## References

### Papers

1. **Shan & Chen (1993)**: "Lattice Boltzmann model for simulating flows with multiple phases and components"
   - Original pseudo-potential model
   - PRE 47(3), 1815-1819

2. **Shan & Chen (1994)**: "Simulation of nonideal gases and liquid-gas phase transitions by the lattice Boltzmann equation"
   - Extended model with equation of state
   - PRE 49(4), 2941-2948

3. **Yuan & Schaefer (2006)**: "Equations of state in a lattice Boltzmann model"
   - Comprehensive review of pseudo-potential models
   - Physics of Fluids 18, 042101

### lbmpy Documentation

- **DeepWiki**: https://deepwiki.com/lssfau/lbmpy
- **Section 5.1**: Shan-Chen Two-Phase Model
- **Section 5.2**: Shan-Chen Two-Component Model
- **Tutorial Notebooks**:
  - `doc/notebooks/08_tutorial_shanchen_twophase.ipynb`
  - `doc/notebooks/09_tutorial_shanchen_twocomponent.ipynb`

### Test Cases

- `tests/shan_chen/test_shan_chen_two_phase.py`: Reference implementation
- `tests/shan_chen/test_shan_chen_two_component.py`: Binary mixture example

## Contact & Support

For issues or questions:
- Check lbmpy documentation: https://deepwiki.com/lssfau/lbmpy
- Review existing validation scripts in `experiments/lbmpy_validation/`
- Consult project README.md in repository root
- Refer to CLAUDE.md for project context

---

**Last Updated**: 2025-11-15
**Status**: Phase 1 (Static Droplet) Implemented ✓ | Phases 2-4 Pending ⚠️
