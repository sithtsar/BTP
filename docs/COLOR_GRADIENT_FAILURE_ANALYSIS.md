# Color-Gradient Model: Complete Failure Analysis

## Executive Summary

**Problem**: Surface tension σ drops from initial value to 0.000 at t=500 and never recovers.

**Root Cause**: Color-Gradient model lacks mechanism to prevent **total density diffusion**. While the model successfully maintains color separation (|∇ρ^N| ≈ 0.27), the total density ρ_total = ρ_red + ρ_blue diffuses to uniform, causing pressure to collapse to zero.

**Solution**: Integrate **Van der Waals forces** (F = -κ·ρ·∇∇²ρ) to maintain density gradients and prevent diffusion.

---

## Investigation Timeline

### Stage 1: Initial Misdiagnosis (FAILED FIX)

**Hypothesis**: Perturbation amplitude 4000x too large
**Evidence**: Initial σ = 3.964 (vs target 0.001)
**Fix Attempted**: Reduced amplitude from `/3.0` to `/12000.0`
**Result**: **NO EFFECT** - identical failure pattern

**Correction**: The large initial pressure ΔP ≈ 0.13 comes from **density ratio** via equation of state P = ρ·cs², NOT from perturbation forces!

### Stage 2: Gradient Threshold Death Spiral (PARTIAL FIX)

**Discovery**: Gradient thresholds in perturbation and recoloring created vicious cycle:
```
Interface diffuses → gradient weakens → threshold deactivates forces →
faster diffusion → complete collapse
```

**Fix**: Removed thresholds in `applyPerturbation()` and `recolor()` functions

**Evidence of Success**:
```
t=0:     |∇ρ^N| = 0.400
t=500:   |∇ρ^N| = 0.295
t=10000: |∇ρ^N| = 0.294  ← Maintained!
```

Gradient is 27 million times stronger than old threshold (1e-9)!

**But**: Surface tension still shows σ = 0.000. Why?

### Stage 3: Total Density Diffusion (ROOT CAUSE IDENTIFIED)

**Key Measurements**:

| Timestep | ρ_total Range | ρ_red Range | ρ_blue Range | Color Field ρ^N | Pressure ΔP |
|----------|---------------|-------------|--------------|-----------------|-------------|
| t=0 | [0.80, 1.20] | [0.00, 1.20] | [0.00, 0.80] | [-1.0, +1.0] | 0.133 |
| t=500 | [0.879, 0.886] | [0.123, 0.711] | [0.173, 0.758] | [-0.72, +0.61] | 0.000 |
| t=10000 | [0.8823, 0.8826] | [0.320, 0.446] | [0.437, 0.563] | [-0.28, +0.01] | 0.000 |

**Critical Insight**:
- Color gradient |∇ρ^N| **maintained** at 0.294 (threshold fix worked!)
- Total density ρ_total **diffused** to nearly uniform (Δρ = 0.0003)
- Pressure P = ρ·cs² requires **density gradient**, not just color gradient
- With uniform density → uniform pressure → ΔP = 0 → σ = 0

### Stage 4: Understanding the Physics

**What Color-Gradient Does**:
1. **Perturbation**: Creates anisotropic stress from color gradient → surface tension
2. **Recoloring**: Maximizes phase separation → maintains |∇ρ^N|
3. **Equation of state**: P = ρ·cs² (ideal gas, **compressible**)

**What's Missing**:
- No force to **prevent total density diffusion**
- Compressible EOS allows ρ_red and ρ_blue to spread spatially
- Even with perfect color separation, if both phases diffuse equally, total density becomes uniform

**Real Immiscible Fluids**:
- **Surface tension**: From molecular forces at interface ✓ (Color-Gradient has this)
- **Bulk incompressibility**: ∇·u = 0, prevents density diffusion ✗ (Missing!)

### Stage 5: Solution Found - Van der Waals Forces

**Existing Implementation**: `include/forces/vdw_force.h`

**Key Components**:

1. **Korteweg Stress**: F = -κ·ρ·∇(∇²ρ)
   - Opposes density gradients
   - Maintains sharp interfaces
   - Prevents diffusion

2. **Van der Waals Equation of State**:
   ```
   P = ρ·cs² - a·ρ² + b·ρ³
   ```
   - `-a·ρ²`: Attractive term → phase separation
   - `+b·ρ³`: Repulsive term → prevents overlap
   - Naturally creates coexistence of two phases

---

## Technical Deep Dive

### Color-Gradient Three-Step Algorithm

**Current Implementation**:
```cpp
// Step 1: ELBM Collision on total distribution f = f_red + f_blue
solver.collide(f_total);

// Step 2: Perturbation (surface tension from color gradient)
cg_model.applyPerturbation(node, f_total);

// Step 3: Recoloring (maximize phase separation)
cg_model.recolor(node, f_total);

// Streaming
stream(grid);
```

**Perturbation Formula**:
```
Δf_i = A·w_i·|∇ρ^N|·[(c_i·n)² - c_s²]
```
- Creates anisotropic stress perpendicular to interface
- Strength proportional to color gradient |∇ρ^N|
- **Works**: Color gradient maintained at 0.294

**Recoloring Formula**:
```
f_i^R = (ρ^R/ρ)·f_i + β·ρ^R·ρ^B·cos(φ_i)·|f_i^eq|/ρ²
f_i^B = (ρ^B/ρ)·f_i - β·ρ^R·ρ^B·cos(φ_i)·|f_i^eq|/ρ²
```
- Redistributes distributions to maximize color separation
- β = 0.9 (increased from 0.7)
- **Works**: Slows diffusion, but can't stop it

### Why Perturbation Alone Isn't Enough

**Perturbation creates surface tension σ**, which causes **Laplace pressure**:
```
ΔP_Laplace = σ/R  (curvature-dependent)
```

**BUT**, the actual measured pressure comes from **equation of state**:
```
P = ρ·cs²
```

If total density is uniform (ρ ≈ 0.88 everywhere):
- P_inside = 0.88 × (1/3) = 0.293
- P_outside = 0.88 × (1/3) = 0.293
- ΔP_measured = 0.293 - 0.293 = 0.000

**Conclusion**: Surface tension forces create stress, but **compressible EOS** allows phases to expand/compress, negating the pressure difference.

### Comparison with Reference Implementation

**Two-Phase "Simple" Model** (`test_twophase.cpp`):
- Just collision + streaming
- **NO forces** (control experiment)
- Shows **IDENTICAL failure**: R=48.97, width=18.23, ΔP=0
- This is the **natural diffusion behavior** of LBM

**Conclusion**: Both models suffer from same root cause - compressible EOS allows density diffusion.

---

## Solution Paths

### Option 1: Van der Waals Force Integration ⭐ RECOMMENDED

**Approach**: Replace Color-Gradient with Van der Waals model

**Implementation**:
```cpp
// Initialize VdW force calculator
VanDerWaalsForce<D2Q9, NX, NY> vdw_force(kappa, a, b);

// In time loop:
for each node:
    // Compute VdW force
    F = vdw_force.computeForce(grid, x, y);

    // ELBM collision with force
    solver.collideWithForce(node, F);

    // Use VdW pressure
    node.fluid.p = vdw_force.computePressure(node.fluid.rho);
```

**Advantages**:
- Naturally maintains density gradients via Korteweg stress
- Non-ideal EOS creates phase separation
- Well-established in LBM literature
- Already implemented in codebase!

**Parameters**:
- `κ`: Interface width parameter (~ 0.01-0.1)
- `a`: Attractive strength (~ 0.1)
- `b`: Repulsive strength (~ 0.1)

### Option 2: Hybrid Color-Gradient + Van der Waals

**Approach**: Add VdW forces to Color-Gradient model

**Implementation**:
```cpp
// Three-step Color-Gradient
collide(f_total);
applyPerturbation(f_total);  // Color-based surface tension
recolor(f_total);

// PLUS: Add VdW force to maintain density gradients
F_vdw = -κ·ρ·∇(∇²ρ);
applyForce(f_total, F_vdw);
```

**Advantages**:
- Combines both approaches
- Color-Gradient provides tunable surface tension
- VdW prevents density diffusion

**Complexity**: More parameters to tune

### Option 3: Incompressible LBM Formulation

**Approach**: Use projection method or pressure-correction

**Challenges**:
- Major algorithmic change
- Not compatible with ELBM framework
- Beyond scope of current project

---

## Validation Plan

### Test 1: Van der Waals Baseline ⭐ HIGH PRIORITY

**Setup**:
- Use existing `vdw_force.h` implementation
- Same droplet configuration (R=30, ρ_liquid=1.2, ρ_gas=0.8)
- Parameters: κ=0.05, a=0.1, b=0.1

**Expected**:
- Density gradients maintained
- Non-zero pressure difference (from non-ideal EOS)
- Droplet maintains shape

**Success Criteria**:
- ρ_total range remains [0.8, 1.2] throughout simulation
- ΔP ≠ 0 at equilibrium
- Interface width < 5 lattice units

### Test 2: Parameter Sweep

**Variables**:
- κ ∈ [0.01, 0.1, 0.5] (interface width)
- a ∈ [0.05, 0.1, 0.2] (attraction)
- b ∈ [0.05, 0.1, 0.2] (repulsion)

**Metrics**:
- Surface tension σ
- Interface width
- Spurious currents

### Test 3: Hybrid Model (Optional)

**If VdW baseline succeeds**, test adding Color-Gradient perturbation on top

---

## Key Findings Summary

### What Worked ✓

1. **Gradient threshold removal**: Maintained |∇ρ^N| ≈ 0.29 (vs collapse to 1e-9)
2. **Recoloring strength β=0.9**: Slowed color diffusion
3. **ELBM stabilization**: No spurious NaN, smooth evolution
4. **Mass conservation**: Perfect (14458 lattice units maintained)

### What Failed ✗

1. **Total density diffusion**: ρ_total → uniform despite color separation
2. **Pressure collapse**: ΔP → 0 due to uniform density
3. **Surface tension measurement**: σ → 0 (consequence of ΔP=0)
4. **Perturbation forces alone**: Insufficient to prevent compressible diffusion

### Root Cause

**Color-Gradient model with ideal gas EOS lacks bulk incompressibility**. The perturbation creates surface tension (curvature pressure), but the compressible equation of state P = ρ·cs² allows total density to diffuse, negating any pressure difference.

### Solution

**Integrate Van der Waals forces** which:
1. Add Korteweg stress F = -κ·ρ·∇(∇²ρ) to oppose density gradients
2. Use non-ideal EOS P = ρ·cs² - a·ρ² + b·ρ³ for phase separation
3. Naturally maintain both density gradients AND surface tension

---

## References

### Color-Gradient Literature
- **Rothman & Keller (1988)**: Original color-gradient model
- **Latva-Kokko & Rothman (2005)**: Diffusion properties and recoloring
- **ETH Research Paper (2024)**: Entropic stabilization at interfaces

### Van der Waals/Free-Energy Literature
- **Swift et al. (1995)** PRL 75:830: Korteweg stress formulation
- **Mazloomi et al. (2015)** PRL 114:174502: ELBM extension for multiphase
- **Krüger et al. (2017)**: LBM textbook Chapter 7

### Implementation Files
- `include/multiphase/color_gradient.h`: Color-Gradient implementation
- `include/forces/vdw_force.h`: Van der Waals forces ⭐
- `src/test_colorgradient.cpp`: Current failing test
- `TWOPHASE_README.md`: Two-phase documentation

---

**Date**: 2025-11-24
**Status**: Complete failure analysis, solution identified
**Next Step**: Test Van der Waals implementation (`vdw_force.h`)
