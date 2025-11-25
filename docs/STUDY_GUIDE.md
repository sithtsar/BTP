# Comprehensive Study Guide: Entropic Lattice Boltzmann Method

**Author**: Sarthak Mishra  
**Project**: BTP-I Entropic Lattice Boltzmann Method  
**Advisor**: Prof. Amol Subhedar  
**Date**: November 2025

---

## Table of Contents

1. [Mathematical Foundations](#1-mathematical-foundations)
2. [Lattice Structures](#2-lattice-structures)
3. [Equilibrium Distributions](#3-equilibrium-distributions)
4. [Collision Operators](#4-collision-operators)
5. [Boundary Conditions](#5-boundary-conditions)
6. [Flow Geometries and Test Cases](#6-flow-geometries-and-test-cases)
7. [Extended Analytical Cases](#7-extended-analytical-cases)
8. [Multiphase Flow](#8-multiphase-flow)
9. [Implementation Details](#9-implementation-details)
10. [Key Derivations](#10-key-derivations)

---

## 1. Mathematical Foundations

### 1.1 Lattice Boltzmann Equation

The fundamental LBM equation discretizes the Boltzmann equation in velocity, space, and time:

$$
f_i(\mathbf{x} + \mathbf{c}_i \Delta t, t + \Delta t) = f_i(\mathbf{x}, t) + \Omega_i(\mathbf{x}, t)
$$

where:
- $f_i(\mathbf{x}, t)$: Distribution function at position $\mathbf{x}$, time $t$, velocity direction $i$
- $\mathbf{c}_i$: Discrete velocity vector for direction $i$
- $\Omega_i$: Collision operator
- $\Delta t$: Time step (typically $\Delta t = 1$ in lattice units)

### 1.2 Macroscopic Variables

From distribution functions, we compute macroscopic quantities:

**Density**:
$$
\rho(\mathbf{x}, t) = \sum_{i=0}^{Q-1} f_i(\mathbf{x}, t)
$$

**Velocity**:
$$
\mathbf{u}(\mathbf{x}, t) = \frac{1}{\rho(\mathbf{x}, t)} \sum_{i=0}^{Q-1} f_i(\mathbf{x}, t) \mathbf{c}_i
$$

**Pressure** (isothermal equation of state):
$$
p(\mathbf{x}, t) = c_s^2 \rho(\mathbf{x}, t)
$$

where $c_s$ is the speed of sound in the lattice.

### 1.3 H-Function (Discrete Entropy)

The H-function is the discrete analog of entropy:

$$
H = \sum_{i=0}^{Q-1} f_i \ln\left(\frac{f_i}{w_i}\right)
$$

where:
- $f_i$: Distribution function in direction $i$
- $w_i$: Lattice weight for direction $i$

**Key Properties**:
- $H \geq 0$ (non-negative)
- Minimum when $f_i = w_i \rho$ (equilibrium)
- Convex function (guarantees unique minimum)

### 1.4 H-Theorem

The H-theorem states that entropy must not increase:

$$
\frac{dH}{dt} \leq 0
$$

This is the **fundamental stability guarantee** for ELBM. If entropy never increases, the system cannot become unstable.

### 1.5 Viscosity Relations

**BGK Viscosity**:
$$
\nu = c_s^2 \left(\tau - \frac{\Delta t}{2}\right) = c_s^2 \left(\frac{1}{\omega} - \frac{1}{2}\right) \Delta t
$$

where:
- $\tau$: Relaxation time
- $\omega = 1/\tau$: Relaxation parameter

**ELBM Viscosity**:
$$
\nu = c_s^2 \left(\frac{1}{\alpha \beta} - \frac{1}{2}\right) \Delta t
$$

where:
- $\alpha$: Iso-entropy parameter (solved to preserve entropy)
- $\beta$: Dissipation parameter (controls viscosity)

### 1.6 Stability Conditions

**BGK Stability**:
- Requires: $\omega < 2$ or $\tau > \Delta t/2$
- Fails at high Reynolds numbers (Re > 100-200)
- Numerical diffusion at low viscosity

**ELBM Stability**:
- **Unconditionally stable** (no upper limit on Re)
- H-theorem enforcement prevents instability
- Works even at vanishing viscosity

---

## 2. Lattice Structures

### 2.1 D2Q9 Lattice (2D, 9 Velocities)

**Discrete Velocities**:
$$
\mathbf{c}_i = \begin{cases}
(0, 0) & i=0 \text{ (rest)} \\
(\pm 1, 0), (0, \pm 1) & i=1,2,3,4 \text{ (cardinal)} \\
(\pm 1, \pm 1) & i=5,6,7,8 \text{ (diagonal)}
\end{cases}
$$

**Lattice Weights**:
$$
w_i = \begin{cases}
4/9 & i=0 \\
1/9 & i=1,2,3,4 \\
1/36 & i=5,6,7,8
\end{cases}
$$

**Speed of Sound**:
$$
c_s = \frac{1}{\sqrt{3}}, \quad c_s^2 = \frac{1}{3}
$$

**Opposite Directions** (for bounce-back):
$$
\text{opposite} = [0, 3, 4, 1, 2, 7, 8, 5, 6]
$$

### 2.2 D3Q19 Lattice (3D, 19 Velocities)

**Discrete Velocities**:
- Rest: $(0,0,0)$
- Face directions: $(\pm 1,0,0), (0,\pm 1,0), (0,0,\pm 1)$
- Edge directions: $(\pm 1,\pm 1,0), (\pm 1,0,\pm 1), (0,\pm 1,\pm 1)$

**Lattice Weights**:
$$
w_i = \begin{cases}
1/3 & i=0 \text{ (rest)} \\
1/18 & i=1-6 \text{ (face)} \\
1/36 & i=7-18 \text{ (edge)}
\end{cases}
$$

**Speed of Sound**: Same as D2Q9: $c_s^2 = 1/3$

### 2.3 Lattice Units

All simulations use **lattice units**:
- $\Delta x = 1$ (lattice spacing)
- $\Delta t = 1$ (time step)
- $c_s = 1/\sqrt{3}$ (speed of sound)
- $c_s^2 = 1/3$

**Conversion to Physical Units**:
- Physical length: $L_{phys} = L_{lattice} \times \Delta x_{phys}$
- Physical time: $t_{phys} = t_{lattice} \times \Delta t_{phys}$
- Physical velocity: $u_{phys} = u_{lattice} \times \frac{\Delta x_{phys}}{\Delta t_{phys}}$

---

## 3. Equilibrium Distributions

### 3.1 BGK Equilibrium (Polynomial)

The BGK equilibrium is a **2nd order Hermite expansion**:

$$
f_i^{eq} = w_i \rho \left[1 + \frac{\mathbf{c}_i \cdot \mathbf{u}}{c_s^2} + \frac{(\mathbf{c}_i \cdot \mathbf{u})^2}{2c_s^4} - \frac{\mathbf{u}^2}{2c_s^2}\right]
$$

**Expanded Form** (for D2Q9):
$$
f_i^{eq} = w_i \rho \left[1 + 3(\mathbf{c}_i \cdot \mathbf{u}) + \frac{9}{2}(\mathbf{c}_i \cdot \mathbf{u})^2 - \frac{3}{2}\mathbf{u}^2\right]
$$

**Properties**:
- Moment-matching approach (matches first two moments)
- **Does NOT minimize entropy functional**
- Simple to compute
- Root cause of BGK instabilities

### 3.2 Entropic Equilibrium

The entropic equilibrium **minimizes the H-function** subject to mass and momentum constraints:

$$
f_i^{eq} = w_i \rho \prod_{j=1}^{D} \left(2 - \sqrt{1 + u_j^2}\right) \left(\frac{2u_j + \sqrt{1 + 3u_j^2}}{1 - u_j}\right)^{c_{ij}/c}
$$

where:
- $D$: Number of spatial dimensions
- $u_j$: Velocity component in dimension $j$
- $c_{ij}$: Component of $\mathbf{c}_i$ in dimension $j$
- $c$: Lattice speed (typically $c = 1$)

**For D2Q9**:
$$
f_i^{eq} = w_i \rho (2 - \sqrt{1 + u_x^2})(2 - \sqrt{1 + u_y^2}) \left(\frac{2u_x + \sqrt{1 + 3u_x^2}}{1 - u_x}\right)^{c_{ix}} \left(\frac{2u_y + \sqrt{1 + 3u_y^2}}{1 - u_y}\right)^{c_{iy}}
$$

**Properties**:
- **Minimizes entropy functional** $H$
- Guarantees H-theorem: $dH/dt \leq 0$
- More complex than BGK but necessary for stability
- Self-adjusting sound speed prevents lattice violations

**Key Difference from BGK**:
- BGK: Moment-matching (matches 2nd moment)
- Entropic: Entropy-minimization (does NOT match 2nd moment exactly)
- Pressure tensor deviations negligible for $Ma < 0.5$

---

## 4. Collision Operators

### 4.1 BGK Collision Operator

The BGK (Bhatnagar-Gross-Krook) collision is the simplest:

$$
f_i^{new} = f_i - \omega (f_i - f_i^{eq})
$$

where $\omega = 1/\tau$ is the relaxation parameter.

**Stability Constraint**:
$$
\omega < 2 \quad \text{or} \quad \tau > \frac{\Delta t}{2}
$$

**Limitations**:
- Fails at high Reynolds numbers (Re > 100-200)
- Numerical diffusion at low viscosity
- No entropy guarantee

### 4.2 ELBM Two-Step Collision

The Entropic LBM uses a **two-step collision process**:

#### Step 1: Iso-Entropy Relaxation (α-relaxation)

$$
f^* = f + \alpha (f^{eq} - f)
$$

where $\alpha$ is solved such that:
$$
H(f^*) = H(f)
$$

This step **preserves entropy** (no dissipation).

#### Step 2: Dissipation Relaxation (β-relaxation)

$$
f^{new} = (1 - \beta) f + \beta f^*
$$

This step introduces **viscous dissipation** controlled by $\beta$.

**Combined**:
$$
f^{new} = (1 - \beta) f + \beta [f + \alpha (f^{eq} - f)] = f + \alpha \beta (f^{eq} - f)
$$

### 4.3 Alpha-Parameter Solver

The $\alpha$ parameter is found using **Newton-Raphson iteration**:

**Objective**: Solve $H(f + \alpha \Delta) = H(f)$ where $\Delta = f^{eq} - f$

**Bounds** (to ensure $f + \alpha \Delta \geq 0$):
- When $\Delta_i < 0$: $\alpha \leq \frac{f_i}{|\Delta_i|}$ (upper bound)
- When $\Delta_i > 0$: $\alpha \geq -\frac{f_i}{\Delta_i}$ (lower bound, usually 0)

**Bounds Computation**:
$$
\alpha_{min} = \max_i \left\{0, -\frac{f_i}{\Delta_i} \text{ for } \Delta_i > 0\right\}
$$

$$
\alpha_{max} = \min_i \left\{2, \frac{f_i}{|\Delta_i|} \text{ for } \Delta_i < 0\right\}
$$

**Newton-Raphson Update**:
$$
\alpha_{n+1} = \alpha_n - \frac{H(f + \alpha_n \Delta) - H(f)}{dH/d\alpha}
$$

where:
$$
\frac{dH}{d\alpha} = \sum_i \Delta_i \left[\ln\left(\frac{f_i + \alpha \Delta_i}{w_i}\right) + 1\right]
$$

**Convergence**: Typically 10-15 iterations, tolerance $10^{-10}$

### 4.4 Viscosity Control in ELBM

The effective relaxation parameter in ELBM is:
$$
\omega_{eff} = \frac{1}{\alpha \beta}
$$

Therefore:
$$
\nu = c_s^2 \left(\frac{1}{\alpha \beta} - \frac{1}{2}\right) \Delta t
$$

**Key Insight**: By independently controlling $\alpha$ and $\beta$, we can achieve arbitrarily low viscosity while maintaining stability.

---

## 5. Boundary Conditions

### 5.1 Operation Order (CRITICAL!)

The **correct order** is:
```
1. Collision
2. Streaming
3. Boundary Conditions ← MUST BE LAST!
```

**Why?** Streaming moves distributions to neighbor nodes. If BCs are applied before streaming, the streaming step will overwrite the boundary values.

### 5.2 Bounce-Back (No-Slip Walls)

**Principle**: Reflect distributions in opposite directions.

$$
f_i(\text{wall}) = f_{\bar{i}}(\text{wall})
$$

where $\bar{i}$ is the opposite direction of $i$.

**Implementation**:
```cpp
for (int i = 0; i < Q; ++i) {
    f[i] = f_old[opposite[i]];
}
```

**Properties**:
- Simple and effective
- Enforces zero velocity at wall
- Works for stationary walls
- 2nd order accurate

### 5.3 Zou-He Scheme (Pressure/Velocity Boundaries)

The Zou-He scheme extrapolates **unknown distributions** from known ones while enforcing mass and momentum conservation.

#### Pressure Boundary (Left, x=0)

**Unknown distributions**: $f_1, f_5, f_8$ (pointing right)

**Known**: Density $\rho$ (from pressure: $p = c_s^2 \rho$)

**Compute velocity from known distributions**:
$$
u_y = \frac{f_2 + f_5 + f_6 - f_4 - f_7 - f_8}{\rho}
$$

$$
u_x = 1 - \frac{f_0 + f_2 + f_4 + 2(f_3 + f_6 + f_7)}{\rho}
$$

**Set unknown populations**:
$$
f_1 = f_3 + \frac{2}{3} \rho u_x
$$

$$
f_5 = f_7 - \frac{1}{2}(f_2 - f_4) + \frac{1}{2} \rho u_y + \frac{1}{6} \rho u_x
$$

$$
f_8 = f_6 + \frac{1}{2}(f_2 - f_4) - \frac{1}{2} \rho u_y + \frac{1}{6} \rho u_x
$$

#### Velocity Boundary (Left, x=0)

**Known**: Velocity $\mathbf{u} = (u_x, u_y)$

**Compute density from known distributions**:
$$
\rho = \frac{f_0 + f_2 + f_4 + 2(f_3 + f_6 + f_7)}{1 - u_x}
$$

**Set unknown populations**: Same as pressure boundary.

#### Top/Bottom Boundaries

Similar logic applies, but with different unknown distributions:
- **Top** (y = ny-1): Unknown are $f_4, f_7, f_8$ (pointing down)
- **Bottom** (y = 0): Unknown are $f_2, f_5, f_6$ (pointing up)

### 5.4 Periodic Boundaries

**Principle**: Wrap around domain boundaries.

**Implementation** (during streaming):
```cpp
x_new = (x + c_x[i] + nx) % nx;
y_new = (y + c_y[i] + ny) % ny;
```

**Use Cases**:
- Fully periodic domains (Taylor-Green vortex)
- Inlet/outlet for channel flows (with pressure BCs)

---

## 6. Flow Geometries and Test Cases

### 6.1 Rectangular Channel Flow (Primary Validation)

**Geometry**: 100 × 40 lattice units

**Boundary Conditions**:
- Left (x=0): Pressure inlet ($\rho_{in} = 1.0 + \Delta p / c_s^2$)
- Right (x=99): Pressure outlet ($\rho_{out} = 1.0$)
- Top/Bottom (y=0, y=39): No-slip walls (bounce-back)

**Pressure Drop**: $\Delta p = 5.0$ Pa

**Three Reynolds Number Cases**:

| Case | Re | $\nu$ (lattice) | Expected Behavior |
|------|----|-----------------|-------------------|
| 1 | ~10 | 0.1 | Both BGK/ELBM stable |
| 2 | ~100 | 0.01 | BGK diffusive, ELBM stable |
| 3 | ~1000 | 0.001 | **BGK fails, ELBM stable** |

**Key Results**:
- Case 1: Both methods stable, comparable
- Case 2: BGK shows excessive diffusion, ELBM maintains pressure gradient
- Case 3: BGK produces NaN, ELBM remains stable

### 6.2 Couette Flow

**Geometry**: Flow between two parallel plates

**Boundary Conditions**:
- Bottom (y=0): Stationary wall ($u = 0$)
- Top (y=H): Moving wall ($u = U$)
- Left/Right: Periodic

**Analytical Solution**:
$$
u(y) = U \frac{y}{H}
$$

**Linear velocity profile** from bottom to top.

**Validation**: L2 error < 0.001

### 6.3 Poiseuille Flow

**Geometry**: Pressure-driven flow between parallel plates

**Boundary Conditions**:
- Top/Bottom: No-slip walls
- Left/Right: Periodic (or pressure BCs)
- **Body force**: $F_x = -\frac{dp}{dx}$ to sustain flow

**Analytical Solution**:
$$
u(y) = -\frac{dp}{dx} \frac{y(H-y)}{2\mu}
$$

where $\mu = \rho \nu$ is the **dynamic viscosity**.

**Parabolic velocity profile** with maximum at centerline.

**Key Point**: Must use **dynamic viscosity** $\mu = \rho \nu$, not kinematic $\nu$!

**Validation**: L2 error < 0.01

### 6.4 Taylor-Green Vortex

**Geometry**: Fully periodic domain

**Initial Conditions**:
$$
u_x(x,y,0) = -U_0 \cos(kx) \sin(ky)
$$

$$
u_y(x,y,0) = U_0 \sin(kx) \cos(ky)
$$

**Analytical Solution** (decaying vortex):
$$
u_x(x,y,t) = -U_0 \cos(kx) \sin(ky) e^{-2k^2 \nu t}
$$

$$
u_y(x,y,t) = U_0 \sin(kx) \cos(ky) e^{-2k^2 \nu t}
$$

**Kinetic Energy Decay**:
$$
E(t) = E_0 e^{-4k^2 \nu t}
$$

**Validation**: Extract viscosity from energy decay, error < 1%

### 6.5 Lid-Driven Cavity

**Geometry**: Square cavity

**Boundary Conditions**:
- Top wall: Moving with velocity $U$ (velocity BC)
- Other walls: No-slip (bounce-back)

**Key Features**:
- Primary vortex in center
- Secondary corner vortices at high Re
- Benchmark: Ghia et al. (1982)

**Reynolds Numbers**:
- Re = 100: Primary vortex
- Re = 1000: Secondary corner vortices

### 6.6 Flow Past Cylinder

**Geometry**: Cylinder in channel

**Boundary Conditions**:
- Inlet: Velocity profile
- Outlet: Pressure outlet
- Cylinder: Bounce-back (immersed boundary)
- Top/Bottom: No-slip walls

**Key Metrics**:
- Drag coefficient: $C_d = \frac{2F_d}{\rho U^2 D}$
- Lift coefficient: $C_l = \frac{2F_l}{\rho U^2 D}$
- Strouhal number (for unsteady flow)

**Reynolds Numbers**:
- Re = 40: Steady flow, $C_d \approx 1.5-1.6$
- Re = 100: Vortex shedding, $C_d \approx 1.3-1.4$

---

## 7. Extended Analytical Cases

### 7.1 Womersley Flow

**Physical Description**: Oscillatory pressure-driven flow (pulsatile flow in arteries)

**Analytical Solution**:
$$
u(y,t) = \text{Re}\left[A \left(1 - \frac{\cosh(\lambda y)}{\cosh(\lambda H/2)}\right) e^{i\omega t}\right]
$$

where:
- $\lambda = \sqrt{i\omega/\nu}$: Complex wave number
- $\omega = 2\pi f$: Angular frequency
- $A$: Oscillation amplitude
- $H$: Channel height

**Womersley Number**:
$$
\alpha = H \sqrt{\frac{\omega}{\nu}}
$$

- $\alpha \ll 1$: Quasi-steady flow (viscous dominant)
- $\alpha \gg 1$: Inertial dominant (flat velocity profile in core)

**Applications**: Cardiovascular flow modeling, unsteady flow validation

### 7.2 Hagen-Poiseuille Flow

**Physical Description**: Pressure-driven flow in circular pipe

**Analytical Solution** (cylindrical coordinates):
$$
u(r) = -\frac{dp}{dx} \frac{R^2 - r^2}{4\mu}
$$

where:
- $R$: Pipe radius
- $r$: Radial distance from centerline
- $\mu = \rho \nu$: Dynamic viscosity

**Volume Flow Rate**:
$$
Q = \frac{\pi R^4 |dp/dx|}{8\mu}
$$

**Key Features**:
- Parabolic velocity profile in radial direction
- Zero velocity at wall ($r = R$)
- Maximum velocity at centerline ($r = 0$)

### 7.3 Stokes Shear Flow

**Physical Description**: Simple shear flow at low Reynolds number

**Analytical Solution**:
$$
u(y) = \gamma y
$$

where $\gamma$ is the constant shear rate ($du/dy = \gamma$).

**Constant Shear Stress**:
$$
\tau = \mu \gamma = \text{constant}
$$

**Validation Criteria**:
- Linear velocity profile
- Constant shear rate throughout domain
- Shear stress independent of position

### 7.4 Kolmogorov Flow

**Physical Description**: 2D flow driven by sinusoidal body force

**Analytical Solution** (laminar steady state):
$$
u(y) = \frac{F}{k^2 \nu} \sin(ky)
$$

where:
- $F$: Forcing amplitude
- $k$: Wave number (typically $k = 2\pi/L$)
- $\nu$: Kinematic viscosity

**Stability**:
- Below critical Re: Stable sinusoidal profile
- Above critical Re: Transition to chaotic/turbulent state

**Applications**: Turbulence transition studies, body force implementation testing

---

## 8. Multiphase Flow

### 8.1 Shan-Chen Pseudo-Potential Model

**Basic Concept**: Phase separation through non-local, density-dependent interaction force.

**Interaction Force**:
$$
\mathbf{F}(\mathbf{x}) = -G \psi(\rho(\mathbf{x})) \sum_{i} w_i \psi(\rho(\mathbf{x} + \mathbf{c}_i)) \mathbf{c}_i
$$

where:
- $G < 0$: Interaction strength (negative for attraction → phase separation)
- $\psi(\rho)$: Pseudo-potential function
- $w_i$: Lattice weights

**Pseudo-Potential**:
$$
\psi(\rho) = \rho_0 \left(1 - e^{-\rho/\rho_0}\right)
$$

### 8.2 Phase Separation Mechanism

1. **High density regions** (liquid): Large $\psi(\rho)$ → strong attractive force → densification
2. **Low density regions** (gas): Small $\psi(\rho)$ → weak force → rarefaction
3. **Interface**: Sharp transition where forces balance
4. **Surface tension**: Emerges naturally from force imbalance at curved interfaces

### 8.3 Laplace Pressure Law

For a static droplet with radius $R$:

$$
\Delta P = P_{inside} - P_{outside} = \frac{\sigma}{R}
$$

where $\sigma$ is the surface tension (implicitly determined by $G$).

**Validation**:
- $\Delta P > 0$ (pressure higher inside droplet)
- $\Delta P \propto 1/R$ (linear with curvature)
- Different $G$ values → different surface tensions

### 8.4 Collision Operators for Multiphase

**MRT (Multiple Relaxation Time)**:
- Relaxes different moments at different rates
- Better stability than BGK/SRT
- Reduces spurious currents near interfaces
- Default choice for multi-phase simulations

**SRT (Single Relaxation Time / BGK)**:
- Simpler, single relaxation parameter $\omega$
- Faster computation
- Higher spurious currents
- Suitable for testing/prototyping

**Key Finding**: SRT is 2× faster than MRT with identical accuracy for static/quasi-static flows.

### 8.5 Rayleigh-Taylor Instability

**Configuration**: Heavy fluid initially above light fluid

**Bond Number**:
$$
\text{Bo} = \frac{\Delta \rho \cdot g \cdot L^2}{\sigma}
$$

- Bo < 1: Surface-tension dominated (stable stratified layers)
- Bo ~ 1-10: Transition regime
- Bo > 10: Gravity dominated (turbulent mixing)

**Results**: Gravitational settling with stable stratified layers for Bo < 1.

---

## 9. Implementation Details

### 9.1 Coordinate System

- **y = 0**: Bottom boundary
- **y = ny-1**: Top boundary
- **x = 0**: Left boundary
- **x = nx-1**: Right boundary
- **x increases**: Left to right
- **y increases**: Bottom to top

### 9.2 Critical Implementation Bugs (Fixed)

#### Bug #1: Wrong Equilibrium in ELBM
- **Issue**: Used BGK equilibrium instead of entropic
- **Impact**: Complete ELBM failure (NaN cascade)
- **Fix**: Use `EntropicEquilibrium` in ELBM solver

#### Bug #2: Alpha Bounds Reversed
- **Issue**: Min/max logic inverted
- **Impact**: Alpha-solver searched incorrect range
- **Fix**: Corrected bound computation logic

#### Bug #3: Boundary Condition Timing
- **Issue**: BCs applied BEFORE streaming
- **Impact**: Streaming overwrites boundary values
- **Fix**: Move `bc_manager.applyBoundaries()` AFTER `solver.stream()`

#### Bug #4: Poiseuille Formula Error
- **Issue**: Used kinematic viscosity $\nu$ instead of dynamic $\mu = \rho \nu$
- **Impact**: Incorrect parabolic profile magnitude
- **Fix**: Multiply by density in formula

### 9.3 Alpha-Solver Implementation

**Algorithm**:
1. Compute initial H: $H_0 = H(f)$
2. Compute delta: $\Delta = f^{eq} - f$
3. Compute bounds: $[\alpha_{min}, \alpha_{max}]$
4. Check boundary cases (test $\alpha_{max}$)
5. Newton-Raphson iteration:
   - Compute $H(\alpha)$ and $dH/d\alpha$
   - Update: $\alpha_{n+1} = \alpha_n - \frac{H(\alpha_n) - H_0}{dH/d\alpha}$
   - Clamp to bounds
6. Return converged $\alpha$

**Convergence**: Typically 10-15 iterations, tolerance $10^{-10}$

### 9.4 Performance

| Solver | Time (2500 steps) | Speedup | Max Stable Re |
|--------|-------------------|---------|---------------|
| BGK | 0.4 seconds | 1.0× | ~100-200 |
| ELBM | 5.0 seconds | 0.09× | **Unlimited** |

**Trade-off**: ELBM is ~11× slower due to $\alpha$-solver, but enables arbitrarily high Re.

---

## 10. Key Derivations

### 10.1 Chapman-Enskog Expansion

The Chapman-Enskog expansion shows that LBM recovers the Navier-Stokes equations in the limit of small Knudsen number.

**Key Result**: The viscosity is related to the relaxation time:
$$
\nu = c_s^2 \left(\tau - \frac{\Delta t}{2}\right)
$$

**Derivation** (brief):
1. Expand distribution function: $f_i = f_i^{(0)} + \epsilon f_i^{(1)} + \ldots$
2. $f_i^{(0)}$ is the equilibrium distribution
3. $f_i^{(1)}$ gives the viscous stress tensor
4. Match with Navier-Stokes: $\tau_{ij} = \mu \left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right)$
5. Result: $\mu = \rho c_s^2 (\tau - \Delta t/2)$

### 10.2 Entropic Equilibrium Derivation

The entropic equilibrium minimizes the H-function subject to constraints.

**Lagrangian**:
$$
\mathcal{L} = H - \lambda_0 \left(\sum_i f_i - \rho\right) - \sum_j \lambda_j \left(\sum_i f_i c_{ij} - \rho u_j\right)
$$

**Variation**:
$$
\frac{\partial \mathcal{L}}{\partial f_i} = \ln\left(\frac{f_i}{w_i}\right) + 1 - \lambda_0 - \sum_j \lambda_j c_{ij} = 0
$$

**Solution**:
$$
f_i^{eq} = w_i \exp\left(\lambda_0 - 1 + \sum_j \lambda_j c_{ij}\right)
$$

**For D2Q9**, this leads to the product form:
$$
f_i^{eq} = w_i \rho \prod_j \left(2 - \sqrt{1 + u_j^2}\right) \left(\frac{2u_j + \sqrt{1 + 3u_j^2}}{1 - u_j}\right)^{c_{ij}/c}
$$

### 10.3 H-Theorem for ELBM

**Theorem**: If $f^{eq}$ minimizes $H$, then $dH/dt \leq 0$.

**Proof Sketch**:
1. Compute $dH/dt$ from collision:
   $$
   \frac{dH}{dt} = \sum_i \frac{\partial H}{\partial f_i} \frac{df_i}{dt}
   $$
2. For iso-entropy step: $H(f^*) = H(f)$ → no change
3. For dissipation step: $f^{new} = (1-\beta)f + \beta f^*$
4. Since $f^*$ is closer to equilibrium than $f$, and equilibrium minimizes $H$:
   $$
   H(f^{new}) \leq H(f)
   $$
5. Therefore: $dH/dt \leq 0$ ✓

### 10.4 Viscosity Relation Derivation

**Starting Point**: Effective relaxation parameter
$$
\omega_{eff} = \frac{1}{\alpha \beta}
$$

**From Chapman-Enskog**:
$$
\nu = c_s^2 \left(\tau - \frac{\Delta t}{2}\right) = c_s^2 \left(\frac{1}{\omega} - \frac{1}{2}\right) \Delta t
$$

**Substitute**:
$$
\nu = c_s^2 \left(\frac{1}{\alpha \beta} - \frac{1}{2}\right) \Delta t
$$

**Key Insight**: By controlling $\alpha$ and $\beta$ independently, we can achieve any viscosity while maintaining stability.

### 10.5 Zou-He Boundary Condition Derivation

**Principle**: Enforce mass and momentum conservation at boundary.

**Mass Conservation**:
$$
\rho = \sum_i f_i
$$

**Momentum Conservation**:
$$
\rho u_x = \sum_i f_i c_{ix}
$$

$$
\rho u_y = \sum_i f_i c_{iy}
$$

**For Left Boundary (x=0)**:
- Unknown: $f_1, f_5, f_8$ (pointing right)
- Known: All others

**From mass conservation**:
$$
\rho = f_0 + f_1 + f_2 + f_3 + f_4 + f_5 + f_6 + f_7 + f_8
$$

**From x-momentum**:
$$
\rho u_x = f_1 - f_3 + f_5 - f_6 - f_7 + f_8
$$

**From y-momentum**:
$$
\rho u_y = f_2 - f_4 + f_5 + f_6 - f_7 - f_8
$$

**Solve for unknowns**: This gives the Zou-He formulas.

---

## Quick Reference: Key Equations

### Fundamental Equations

**Lattice Boltzmann Equation**:
$$
f_i(\mathbf{x} + \mathbf{c}_i \Delta t, t + \Delta t) = f_i(\mathbf{x}, t) + \Omega_i(\mathbf{x}, t)
$$

**H-Function**:
$$
H = \sum_i f_i \ln\left(\frac{f_i}{w_i}\right)
$$

**H-Theorem**:
$$
\frac{dH}{dt} \leq 0
$$

### Collision Operators

**BGK**:
$$
f_i^{new} = f_i - \omega (f_i - f_i^{eq})
$$

**ELBM**:
$$
f^* = f + \alpha (f^{eq} - f), \quad H(f^*) = H(f)
$$

$$
f^{new} = (1 - \beta) f + \beta f^*
$$

### Viscosity

**BGK**:
$$
\nu = c_s^2 \left(\frac{1}{\omega} - \frac{1}{2}\right) \Delta t
$$

**ELBM**:
$$
\nu = c_s^2 \left(\frac{1}{\alpha \beta} - \frac{1}{2}\right) \Delta t
$$

### Stability

**BGK**: $\omega < 2$ or $\tau > \Delta t/2$  
**ELBM**: **Unconditionally stable**

---

## Study Tips for Presentation

1. **Start with the Big Picture**: Why ELBM? (Stability at high Re)
2. **Mathematical Foundation**: H-theorem is the key concept
3. **Two-Step Collision**: Explain why we need both steps
4. **Alpha-Solver**: Mention Newton-Raphson but don't go into details
5. **Boundary Conditions**: Emphasize the correct order (collision → streaming → BC)
6. **Results**: Focus on the three Re cases showing progressive divergence
7. **Extended Analysis**: Mention the four cases briefly
8. **Multiphase**: Highlight the Shan-Chen model and key findings

**Common Questions to Prepare For**:
- Why is ELBM more stable than BGK? (H-theorem enforcement)
- What is the computational cost? (~11× slower but enables high Re)
- How does the alpha-solver work? (Newton-Raphson to preserve entropy)
- What are the limitations? (Computational cost, complexity)
- Future work? (ML integration, 3D, optimization)

---

**End of Study Guide**

*Good luck with your presentation!*

