# lbmpy Validation: BGK vs ELBM vs Analytical

This directory contains lbmpy-based implementations that compare BGK (standard LBM), ELBM (Entropic LBM), and analytical solutions for three canonical fluid dynamics test cases.

## Overview

All scripts use the **lbmpy** library to implement both BGK and ELBM collision operators, validating that:
1. Both methods converge to the analytical solution at low Reynolds numbers
2. ELBM provides similar accuracy to BGK when both are stable
3. The lbmpy implementation matches theoretical predictions

## Test Cases

### 1. Couette Flow (`couette_lbmpy.py`)
**Geometry**: Shear-driven flow between parallel plates
- **Bottom plate**: Stationary (u=0)
- **Top plate**: Moving with velocity U
- **Analytical solution**: u(y) = U · y / H (linear velocity profile)
- **Validation metric**: L2 error < 0.001

**Run**:
```bash
cd experiments/lbmpy_validation
python couette_lbmpy.py
```

**Output**:
- Figure: `figures/lbmpy_validation/couette_3way_comparison.png`
- Shows BGK, ELBM, and analytical profiles overlapping

---

### 2. Poiseuille Flow (`poiseuille_lbmpy.py`)
**Geometry**: Pressure-driven (body force) channel flow
- **Walls**: No-slip boundary conditions at top and bottom
- **Driving force**: Body force Fx (simulating pressure gradient)
- **Analytical solution**: u(y) = -F · y · (H-y) / (2ρν) (parabolic profile)
- **Validation metric**: L2 error < 0.01 (relaxed due to body force implementation)

**Run**:
```bash
cd experiments/lbmpy_validation
python poiseuille_lbmpy.py
```

**Output**:
- Figure: `figures/lbmpy_validation/poiseuille_3way_comparison.png`
- Shows parabolic velocity profiles converging for all three methods

---

### 3. Taylor-Green Vortex (`taylor_green_lbmpy.py`)
**Geometry**: Decaying 2D vortex (fully periodic domain)
- **Initial conditions**: u(x,y,0) = -U₀·cos(kx)·sin(ky), v(x,y,0) = U₀·sin(kx)·cos(ky)
- **Analytical solution**: E(t) = E₀·exp(-4νk²t)
- **Validation metric**: Extracted viscosity error < 1%

**Run**:
```bash
cd experiments/lbmpy_validation
python taylor_green_lbmpy.py
```

**Output**:
- Figure: `figures/lbmpy_validation/taylor_green_3way_comparison.png`
- Left: Energy decay on semilog scale
- Right: Viscosity extraction over time

---

## Running All Validations

To run all three test cases sequentially:

```bash
cd experiments/lbmpy_validation
python run_all_validations.py
```

This will:
1. Run Couette, Poiseuille, and Taylor-Green validations
2. Generate all three comparison figures
3. Report pass/fail status for each case
4. Exit with code 0 if all pass, 1 otherwise

---

## Expected Results

### Couette Flow
- **BGK L2 error**: ~1e-4 to 1e-3
- **ELBM L2 error**: ~1e-4 to 1e-3
- **Observation**: Linear velocity profile, both methods match analytical

### Poiseuille Flow
- **BGK L2 error**: ~1e-3 to 1e-2
- **ELBM L2 error**: ~1e-3 to 1e-2
- **Observation**: Parabolic velocity profile, both methods converge to analytical

### Taylor-Green Vortex
- **BGK viscosity error**: ~0.1% to 0.5%
- **ELBM viscosity error**: ~0.1% to 0.5%
- **Observation**: Exponential energy decay, viscosity correctly extracted

---

## Implementation Details

### BGK Configuration
```python
lbm_config = LBMConfig(
    stencil=LBStencil(Stencil.D2Q9),
    method=Method.SRT,           # Single Relaxation Time (BGK)
    relaxation_rate=omega,
    compressible=True
)
```

### ELBM Configuration
```python
lbm_config = LBMConfig(
    stencil=LBStencil(Stencil.D2Q9),
    method=Method.SRT,
    entropic=True,               # Enable Entropic LBM
    relaxation_rate=omega,
    compressible=True
)
```

The key difference is the `entropic=True` flag, which activates lbmpy's entropic collision operator that enforces the H-theorem (discrete entropy non-increase).

---

## Framework Structure

### `lbmpy_framework.py`
Base simulation class providing:
- D2Q9 lattice setup
- BGK and entropic equilibrium distributions
- Collision operators for both BGK and ELBM
- Streaming, boundary conditions
- Entropy computation (H-function)
- L2 error calculation and convergence checking
- Plotting utilities

**Not currently used** by the main validation scripts (which use lbmpy's native API directly), but provided as a reference implementation showing how ELBM collision works from first principles.

---

## Dependencies

Install required packages:
```bash
uv add lbmpy pystencils sympy matplotlib numpy
```

Or with pip:
```bash
pip install lbmpy pystencils sympy matplotlib numpy
```

---

## Key Physics

### Reynolds Number
```
Re = U·L / ν
```

### Relaxation Parameter
```
ω = 1/τ
τ = 3ν + 0.5  (in lattice units)
```

### BGK Collision
```
f_new = f - ω·(f - f_eq)
```

### ELBM Collision (Simplified)
```
Step 1: f* = f + α·(f_eq - f)  where H(f*) = H(f)
Step 2: f_new = (1-β)·f + β·f*
```

### Discrete Entropy (H-function)
```
H = Σ f_i · ln(f_i / w_i)
```

ELBM enforces dH/dt ≤ 0, providing unconditional stability.

---

## Troubleshooting

### Issue: NaN values in ELBM
**Cause**: Alpha-parameter solver failed or entropic equilibrium has numerical issues
**Fix**: Check that velocities remain < 0.1 (low Mach number approximation)

### Issue: Slow convergence
**Cause**: Need more time steps for diffusion-dominated flows
**Fix**: Increase `time_steps` parameter (e.g., 20000 → 50000)

### Issue: High L2 errors
**Cause**: Boundary conditions not applied correctly
**Fix**: Ensure BCs are applied AFTER streaming step

---

## Comparison with C++ Implementation

These lbmpy scripts are **independent** of the C++ code in `include/` and `src/`. They serve as:
1. **Validation**: Verify C++ implementation correctness
2. **Prototyping**: Quickly test new collision operators
3. **Visualization**: Generate publication-quality comparison plots

The C++ implementation uses the same underlying physics but is optimized for performance and includes additional features like immersed boundary methods for complex geometries.

---

## References

1. **lbmpy Documentation**: https://lbmpy.readthedocs.io
2. **Yu (2021)**: "Entropic Lattice Boltzmann Method" - UNC Honor Thesis
3. **Karlin et al. (1999)**: "Perfect entropy functions of the Lattice Boltzmann method"
4. **Ansumali & Karlin (2000)**: "Stabilization of the lattice Boltzmann method by the H theorem"

---

## Output Structure

```
BTP_FINAL/
├── experiments/lbmpy_validation/
│   ├── couette_lbmpy.py
│   ├── poiseuille_lbmpy.py
│   ├── taylor_green_lbmpy.py
│   ├── lbmpy_framework.py          # Reference implementation
│   ├── run_all_validations.py      # Master script
│   └── README.md                    # This file
│
└── figures/lbmpy_validation/
    ├── couette_3way_comparison.png
    ├── poiseuille_3way_comparison.png
    └── taylor_green_3way_comparison.png
```

---

## Notes

- All simulations use D2Q9 stencil (9 discrete velocities)
- Lattice units: Δx = 1, Δt = 1, cs = 1/√3
- BGK uses standard 2nd-order Hermite equilibrium
- ELBM uses entropic equilibrium (when `entropic=True`) or standard equilibrium with entropy-based collision
- All cases run until L2 error converges below threshold or max steps reached

---

Last Updated: 2025-11-15
