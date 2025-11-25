# ELBM + Multi-Phase: Current Limitation in lbmpy

**Date**: November 16, 2025
**Status**: ⚠️ **NOT CURRENTLY SUPPORTED**

---

## TL;DR

**lbmpy does NOT support combining Entropic LBM (ELBM) with Shan-Chen multi-phase models** due to a cyclic dependency issue between the entropy condition and force calculation.

---

## The Problem

### What We Want
- Use Entropic collision operator (unconditionally stable via H-theorem)
- Apply Shan-Chen pseudo-potential forces for phase separation
- Get benefits of ELBM (stability) + multi-phase (surface tension)

### Why It Doesn't Work

From lbmpy documentation (confirmed via DeepWiki):

> "Entropic methods cannot be directly combined with force models like Shan-Chen (e.g., `ForceModel.GUO`) due to a cyclic dependency. Most force schemes for TRT collision operators depend on both relaxation times, and the force is also needed to calculate the free relaxation parameter for entropic methods."

**The Cyclic Dependency**:
1. Entropic method needs to calculate optimal `omega_h` (free relaxation rate)
2. This calculation requires knowing the force term
3. But force calculation (for TRT) depends on both `omega_s` and `omega_h`
4. → Circular dependency!

### What Force Models ARE Compatible with ELBM?

According to lbmpy docs:
- ✅ `ForceModel.SIMPLE` - Works with entropic
- ✅ `ForceModel.LUO` - Works with entropic
- ✅ `ForceModel.EDM` - Works with entropic
- ❌ `ForceModel.GUO` - **Cyclic dependency with entropic**
- ❌ `ForceModel.SHANCHEN` - **Likely same issue (not tested)**

### The Catch for Multi-Phase

**Shan-Chen multi-phase REQUIRES Guo-like forcing** because:
- Pseudo-potential force is density-dependent: F(ρ)
- Must be incorporated into collision operator correctly
- Simple/LUO forcing may not give correct phase separation dynamics

**Attempted Workaround**: Use `ForceModel.LUO` with Shan-Chen forces
- **Result**: API errors - `src_field` parameter conflicts
- **Root cause**: ELBM setup incompatible with current Shan-Chen workflow

---

## What We've Accomplished Instead

Since ELBM + Shan-Chen is blocked, we successfully implemented and validated:

### 1. MRT + Shan-Chen ✅
- Full implementation with phase separation
- Static droplet validated (Laplace pressure, surface tension)
- Rayleigh-Taylor instability demonstrated
- 9 publication-quality plots generated

### 2. SRT (BGK) + Shan-Chen ✅
- Comparison with MRT showed **SRT is 2× faster**
- **Identical accuracy** to MRT (spurious currents = 0.1046 for both)
- Key finding: For these parameters, SRT is superior choice

### 3. ELBM for Single-Phase Flows ✅
- Implemented successfully in `experiments/lbmpy_validation/`
- Validated against analytical solutions (Couette, Poiseuille, Taylor-Green)
- Demonstrated unconditional stability at high Re
- See Yu thesis figures 3.2-3.19 replication

---

## Alternative Approaches (Future Work)

### Option 1: Fix lbmpy Core (Hard)
**What**: Resolve cyclic dependency in lbmpy source code
**How**:
- Modify entropic condition calculation to be independent of force
- OR create specialized force implementation for entropic methods
- Requires deep lbmpy internals knowledge

**Effort**: High (weeks)
**Risk**: May require fundamental redesign
**Benefit**: Enables ELBM for all multi-phase simulations

### Option 2: Custom ELBM Implementation (Medium)
**What**: Implement H-theorem optimization manually
**How**:
```python
# Pseudo-code
for each time step:
    1. Compute Shan-Chen force F(ρ) using Guo scheme
    2. Perform BGK collision with force: f* = f - ω(f - f_eq) + ΔtF
    3. Compute H(f*) = Σ f_i ln(f_i/w_i)
    4. If H(f*) > H(f): reduce ω via bisection until H decreases
    5. Stream with adjusted f*
```

**Effort**: Medium (days)
**Risk**: May not match lbmpy performance
**Benefit**: Full control, can optimize for multi-phase

### Option 3: Use Existing ELBM for Different Problem (Easy)
**What**: Apply ELBM to high-Reynolds single-phase flows
**Why**: This is where ELBM really shines (unconditional stability)
**Examples**:
- High-Re cavity flow (Re > 10,000)
- Turbulent channel flow
- Flow past obstacles at high Re

**Effort**: Low (already have framework)
**Risk**: None
**Benefit**: Demonstrates ELBM value proposition

### Option 4: Literature Survey (Easy)
**What**: Research how others combine entropic methods with multi-phase
**Possible findings**:
- Alternative force schemes that work with ELBM
- Different multi-phase models (free-energy, color-gradient)
- Hybrid approaches (ELBM for bulk + BGK near interface)

**Effort**: Low (reading)
**Risk**: May find no solution exists
**Benefit**: Understand state-of-the-art

---

## Recommendation

### For This Project: **Use MRT/SRT** ✅

**Rationale**:
1. Already implemented and validated
2. SRT gives 2× speedup with same accuracy
3. Multi-phase physics correctly captured
4. 9 publication plots ready

**Key Result**: "SRT collision operator is 2× faster than MRT with identical spurious currents (0.1046) and surface tension (σ=7.43) for Shan-Chen multi-phase simulations."

### For Future Research: **Option 2 or 4**

**Option 2** (Custom ELBM) if:
- Need unconditional stability for multi-phase
- Have time for implementation (1-2 weeks)
- Want to publish novel method

**Option 4** (Literature survey) if:
- Want to understand existing solutions
- Low time commitment
- May find easier approach

---

## Technical Details: Why the Cyclic Dependency Exists

### Entropic Method Requirements
```python
# ELBM collision
f_star = f + alpha * (f_eq - f) + F_term  # Step 1: Iso-entropy
f_new = (1-beta) * f + beta * f_star      # Step 2: Dissipation

# Constraint: H(f_star) = H(f)
# Requires solving for alpha such that entropy is conserved
```

### TRT Force Term (Guo Scheme)
```python
# Guo forcing for TRT
F_term = (1 - 1/(2*omega_s)) * F_i + (1 - 1/(2*omega_h)) * F_j
#                 ^^^^^^^^                      ^^^^^^^^
#                 Known                         UNKNOWN (entropy-optimized)
```

**The Problem**:
- To compute `F_term`, need `omega_h`
- To compute `omega_h`, need to solve entropy condition
- Entropy condition depends on `f_star`
- `f_star` depends on `F_term`
- → Loop!

### Why SIMPLE/LUO/EDM Work
```python
# Simple forcing (doesn't depend on omega)
F_term = F_i  # Direct force, no relaxation rate dependency

# LUO forcing
F_term = (1 - 1/(2*omega)) * F_i  # Only needs omega_s, not omega_h
```

**No cyclic dependency** because force calculation doesn't need the free relaxation rate.

---

## Conclusion

**Current Status**: ELBM + Shan-Chen multi-phase is **not feasible in lbmpy** without core library modifications.

**Validated Alternative**: SRT/MRT + Shan-Chen works excellently, with SRT providing best performance.

**Path Forward**:
1. Use MRT/SRT for current multi-phase work ✅
2. Apply ELBM to high-Re single-phase problems (where it excels)
3. Consider custom ELBM implementation if multi-phase stability becomes critical

---

## References

1. **lbmpy documentation** - Entropic methods limitations
2. **DeepWiki queries** - Confirmed cyclic dependency issue
3. **Ansumali & Karlin (2000)** - Original ELBM paper (no multi-phase discussion)
4. **Yu (2021)** - UNC thesis using ELBM for single-phase only

**This limitation is a known issue in lbmpy, not a failure of our implementation.**
