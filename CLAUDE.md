# ELBM Project Context for Claude

## Project Overview

This project replicates Keming Yu's 2021 honor thesis "Entropic Lattice Boltzmann Method" from UNC Chapel Hill. The goal is to implement and validate the entropic LBM approach, which provides unconditional stability compared to standard BGK methods.

## Key Mathematical Foundations

### H-Function (Discrete Entropy) - Equation 3.32-3.33
```
H = Σ f_i · ln(f_i / w_i)
```
where:
- `f_i`: Distribution function at velocity direction i
- `w_i`: Lattice weight for direction i

### H-Theorem
```
dH/dt ≤ 0
```
Entropy must not increase - this is the stability guarantee.

### Two-Step Entropic Collision

**Step 1: Iso-Entropy (α-relaxation)** - Equation 3.40-3.42
```
f* = f + α*(f_eq - f)
```
where α is solved to satisfy: `H(f*) = H(f)`

**Step 2: Dissipation (β-relaxation)** - Equation 3.43
```
f_new = (1-β)f + β·f*
```

### Alpha-Parameter Bounds - Equation 3.46
```
α ∈ [α_min, α_max]
```
where bounds ensure f + α·Δ ≥ 0 for all i:
- When Δ_i < 0: α ≤ f_i/|Δ_i| (upper bound)
- When Δ_i > 0: α ≥ -f_i/Δ_i (lower bound, usually 0)

### Viscosity Relation - Equation 3.47
```
ν = cs²(1/(αβ) - 1/2)Δt
```
where cs² = 1/3 for D2Q9 lattice.

### Equilibrium Distributions

**BGK Equilibrium** (2nd order Hermite expansion):
```
f_i^eq = w_i·ρ[1 + (c_i·u)/cs² + (c_i·u)²/(2cs⁴) - u²/(2cs²)]
```

**Entropic Equilibrium** - Equation 3.39:
```
f_i^eq = w_i·ρ·∏(2-√(1+u_j²))·((2u_j+√(1+3u_j²))/(1-u_j))^(c_ij/c)
```

## Critical Implementation Notes

### Bugs Found & Fixed

1. **Wrong Equilibrium in ELBM** (CATASTROPHIC)
   - Was using BGK equilibrium instead of entropic
   - Caused complete ELBM failure with NaN cascade
   - **Fix**: Use EntropicEquilibrium in ELBM solver
   - **Location**: `include/solvers/elbm_solver.h` line 168

2. **Alpha Bounds Reversed** (CRITICAL)
   - Min/max logic was inverted
   - Alpha-solver searched wrong range
   - **Fix**: Corrected bound computation logic
   - **Location**: `include/solvers/elbm_solver.h` lines 49-87

3. **Boundary Condition Timing** (CATASTROPHIC)
   - BCs applied BEFORE streaming instead of AFTER
   - Streaming overwrites BC values → no boundary enforcement
   - **Fix**: Move bc_manager.applyBoundaries() after solver.stream()
   - **Location**: `src/test_analytical.cpp` lines 58-61 and 140-143

4. **Poiseuille Formula Error**
   - Used kinematic viscosity ν instead of dynamic viscosity μ = ρ·ν
   - **Fix**: Multiply by density in formula
   - **Location**: `include/validation/analytical_cases.h` lines 96 and 116

### Coordinate System

- **y = 0**: Bottom boundary
- **y = ny-1**: Top boundary
- **x = 0**: Left boundary
- **x = nx-1**: Right boundary
- **x increases**: Left to right
- **y increases**: Bottom to top

### Lattice Units

All simulations use lattice units where:
- Δx = 1 (lattice spacing)
- Δt = 1 (time step)
- cs = 1/√3 (speed of sound)
- cs² = 1/3

### Boundary Conditions

**Zou-He Scheme**: Used for velocity and pressure BCs
- Extrapolates unknown distributions from known ones
- Enforces mass/momentum conservation

**Bounce-Back**: Used for no-slip walls
- Reflects distributions: f_i(wall) = f_opp(wall)
- Simple but effective for stationary walls

**Operation Order (CRITICAL)**:
```
1. Collision
2. Streaming
3. Boundary conditions ← MUST BE LAST!
```

## Test Cases & Expected Results

### Rectangular Pipe Flow (Primary Validation)

**Purpose**: Reproduce Yu thesis Figures 3.2-3.19

**Cases**:
- Case 1: Re~10 (ν=0.1) → Both BGK/ELBM stable
- Case 2: Re~100 (ν=0.01) → BGK diffusive, ELBM stable
- Case 3: Re~1000 (ν=0.001) → BGK fails, ELBM stable

**Expected**: 19 figures showing progressive divergence between methods

### Analytical Validation

**Couette Flow**:
- Analytical: u(y) = U·y/H (linear)
- Expected L2 error: < 0.001 after fixes
- Convergence: ~20,000 steps

**Poiseuille Flow**:
- Analytical: u(y) = -(dp/dx)·y·(H-y)/(2ρν) (parabolic)
- Expected L2 error: < 0.001 after fixes
- Requires: Body force to sustain flow

**Taylor-Green Vortex**:
- Analytical: E(t) = E₀·exp(-4k²νt)
- Expected: ν extraction error < 1%
- Validates: Chapman-Enskog viscosity

### Benchmark Cases

**Lid-Driven Cavity**:
- Re=100: Primary vortex, steady state
- Re=1000: Secondary vortices in corners
- Reference: Ghia et al. (1982)

**Flow Past Cylinder** (In Progress):
- Re=40: Steady flow, Cd ≈ 1.5-1.6
- Re=100: Vortex shedding, Cd ≈ 1.3-1.4
- Issues: Immersed boundary needs refinement

## Repository Structure

```
BTP_FINAL/
├── include/                  # Header-only C++ implementation
│   ├── core/                # Lattice structures, fluid state
│   ├── solvers/             # BGK, ELBM, equilibrium
│   ├── boundary/            # Boundary conditions
│   └── validation/          # Test case implementations
├── src/                     # C++ executable programs
│   ├── main.cpp            # Rectangular pipe (primary validation)
│   ├── test_analytical.cpp # Analytical validation (BGK + ELBM)
│   └── test_benchmark.cpp  # Benchmark cases
├── scripts/                 # Automation scripts
│   ├── build.sh            # Build system wrapper
│   ├── run_all_cases.sh    # Run all Re cases
│   ├── verify_results.sh   # Verify simulation outputs
│   └── generate_all_figures.py # Generate all 19 thesis figures
├── plotting/                # Python visualization utilities
│   ├── plot_utils.py
│   ├── plot_analytical_3way.py
│   └── plot_analytical_detailed.py
├── experiments/             # Validation experiments
│   ├── lbmpy_validation/   # LBMpy reference implementations
│   │   ├── couette_lbmpy.py
│   │   ├── poiseuille_lbmpy.py
│   │   ├── taylor_green_lbmpy.py
│   │   ├── compare_with_cpp.py
│   │   └── plot_lbmpy_results.py
│   └── notebooks/          # Marimo comparison notebooks
│       └── cpp_vs_lbmpy_comparison.py
├── notebooks/               # Interactive Marimo notebooks
│   ├── 01_pressure_profiles.py
│   ├── 02_analytical_validation.py
│   └── 03_benchmark_cases.py
├── docs/                    # Documentation
│   ├── papers/             # Original research papers
│   │   ├── ELBM-merged.pdf # Yu thesis (PRIMARY REFERENCE)
│   │   ├── PhysRevLett.114.174502.pdf
│   │   └── duncan-orkwis-2025...pdf
│   ├── ocr_output/         # Structured thesis extraction
│   └── presentation.tex    # LaTeX presentation
├── tools/                   # Development tools
│   ├── ocr_extraction.py   # Mistral OCR thesis extraction
│   └── README.md           # Tool documentation
├── figures/                 # Generated plots (~19 figures)
│   ├── Case 1 Re~10/       # Figures 3.2-3.7
│   ├── Case 2 Re~100/      # Figures 3.8-3.13
│   ├── Case 3 Re~1000/     # Figures 3.14-3.19
│   ├── thesis_replication/
│   ├── analytical_validation/
│   └── benchmark_cases/
├── output/                  # Simulation results (.dat files, gitignored)
├── logs/                    # Simulation logs
├── build/                   # Build artifacts (gitignored)
├── CMakeLists.txt          # CMake build configuration
├── pyproject.toml          # Python dependencies (UV)
├── uv.lock                 # UV lockfile
├── CLAUDE.md               # This file
├── README.md               # User guide
└── RESULTS.md              # Validation results
```

## Common Issues & Solutions

### Issue: ELBM Producing NaN
**Causes**:
1. Using BGK equilibrium in ELBM → Use entropic equilibrium
2. Wrong alpha bounds → Check bound computation logic
3. H-function log(0) → Add epsilon threshold

### Issue: Analytical Tests Show Flat Profiles
**Causes**:
1. **BC timing wrong** (most common) → Apply BCs AFTER streaming
2. Not enough convergence time → Use τ ≈ H²/ν steps
3. Missing body force for Poiseuille → Add setBodyForce()

### Issue: Cylinder Flow NaN
**Causes**:
1. No inlet/outlet BCs → Add velocity inlet, pressure outlet
2. Streaming corrupts solid nodes → Apply bounce-back AFTER streaming
3. Wrong force calculation → Use momentum exchange method

## Performance Expectations

| Solver | Speed | Max Stable Re | Use Case |
|--------|-------|---------------|----------|
| BGK | 1.0x (baseline) | ~100-200 | Low Re, quick iterations |
| ELBM | 0.09x (~11x slower) | Unlimited | High Re, stability critical |

## Key Papers & References

1. **Yu (2021)** - PRIMARY: Honor thesis on ELBM
2. **Ansumali & Karlin (2000)** - H-theorem stabilization
3. **Karlin et al. (1999)** - Perfect entropy functions
4. **Krüger et al. (2017)** - LBM textbook
5. **Ghia et al. (1982)** - Cavity benchmark reference

## Build & Run Commands

```bash
# Build everything
./scripts/build.sh

# Run primary validation (thesis replication)
./scripts/run_all_cases.sh           # Creates all 19 figures

# Run analytical tests
./build/test_analytical              # Couette, Poiseuille, Taylor-Green

# Run benchmarks
./build/test_benchmark               # Cavity, cylinder

# Generate all figures
python scripts/generate_all_figures.py

# Verify results
./scripts/verify_results.sh

# Visualize
source .venv/bin/activate
cd notebooks && marimo edit 01_pressure_profiles.py
```

## Current Status

✅ **Completed**:
- All 19 thesis figures reproduced
- ELBM stability demonstrated at Re~1000
- Analytical tests passing (with fixes)
- Lid-driven cavity validated
- Interactive notebooks created
- Comprehensive documentation

⚠️ **In Progress**:
- Cylinder flow immersed boundary refinement

## Notes for Future Development

- **Optimization**: α-solver could use better initial guess (previous step's α)
- **Parallelization**: Collision step is embarrassingly parallel (OpenMP ready)
- **3D**: D3Q19 implementation exists but untested
- **Multiphase**: Framework ready for van der Waals extension
- **GPU**: Structure suitable for CUDA/HIP port

## Quick Reference: Key Equations

```
BGK Collision:     f_new = f - ω(f - f_eq)
ELBM Collision:    f_new = (1-β)f + β·f*, where H(f*) = H(f)
Viscosity:         ν = cs²(τ - Δt/2)       [BGK]
                   ν = cs²(1/(αβ) - 1/2)Δt [ELBM]
Stability:         ω > 1/2                 [BGK]
                   Always stable           [ELBM]
```

---

Last Updated: 2025-11-14
Project Status: Core Implementation Complete
