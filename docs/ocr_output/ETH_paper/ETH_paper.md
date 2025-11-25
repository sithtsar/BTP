---

<!-- Page 0 -->

# Entropic lattice Boltzmann methods: A review 

S.A. Hosseini ${ }^{\mathrm{a}}$, M. Atif ${ }^{\mathrm{b}}$, S. Ansumali ${ }^{\mathrm{c}}$, I.V. Karlin ${ }^{\mathrm{a}, *}$<br>${ }^{a}$ Department of Mechanical and Process Engineering, ETH Zurich, Zurich, 8092, Switzerland<br>${ }^{\mathrm{b}}$ Computational Science Initiative, Brookhaven National Laboratory, Upton, NY, 11973, USA<br>${ }^{c}$ Engineering Mechanics Unit, Jawaharlal Nehru Centre for Advanced Scientific Research, Jakkur, Bangalore, 560064, India

## A R T I C L E I N F O

MSC:
0000
1111
Keywords:
Lattice Boltzmann method
Entropy
Compressible flows
Multiphase flows
Incompressible flows

## A B S T R A C T

In the late 90's and early 2000's the concept of a discrete H theorem and Lyapunov functionals as a way to ensure stability of lattice Boltzmann solvers was a shift of paradigm in the construction of discrete kinetic solvers and opened the door for new discussions and perspectives on the matter. The entropic construction proposed to reorganize the relaxation collision operator by changing both the equilibrium attractor and relaxation process by introducing a discrete entropy functional and enforcing an H-theorem. The concept has proven to be effective in stabilizing lattice Boltzmann solvers in a variety of different area of applications ranging from isothermal weakly compressible, to fully compressible and multi-phase flows. Here we review basic building blocks of the entropic lattice Boltzmann method and discuss its extension to multiphase and compressible flows.

## 1. Introduction

The lattice Boltzmann method, is a numerical method intended for hydrodynamic simulations which came about in the late 80's [1], initially to cure the statistical noise problem plaguing its ancestor the lattice gas automata [2] by replacing the boolean occupation numbers with ensemble-averaged populations [3]. This first work was then followed by a large number of publications proposing other possible reductions of degrees of freedom in the lattice gas automata with a focus on hydrodynamics, culminating with the reduction of the lattice Boltzmann collision matrix into a diagonal one [4] leading to the now widely used lattice Bhatnagar-Gross-Krook (LBGK) model [5]. The resulting discrete system has turned, over the past decades, into a very popular alternative to classical discrete solvers for the Navier-Stokes equations especially in the incompressible flow limit. This rapid growth in popularity can be explained by a number of observations: The final product is a very simple system of coupled algebraic equations and as noted by S. Succi, "it is difficult to imagine something simpler than this equation to recover the Navier-Stokes equations" [3]. Being, at least initially, for the most part used in the incompressible limit it offers a degree of locality of operators that is incomparable to alternatives such as the Navier-Stokes-Poisson equations and contrary to other fully hyperbolic approximations like Chorin's artificial compressibility method [6] it allows for accurate transient flow simulations.

Early on after the establishment of the LBGK as an efficient tool for hydrodynamic simulations, it was noted that in its original form the discrete solver is quite sensitive to both maximum non-dimensional
velocity and viscosity as low non-dimensional values of viscosity or larger values of velocity lead to instabilities preventing large Reynolds number simulations at reasonable cost. This was the starting point for two decades of research on different strategies to enhance stability of the LBGK; Most of the approaches developed in that context [7-10], although formulated from different perspectives, relied on control/tuning of the relaxation rates of different moments of the distribution function. As such they can be categorized as sub-classes of the multiple relaxation time approach [7].

A change of paradigm was proposed in the form of the entropic LBGK; The entropic LBGK in its modern form relies on two main ingredients, namely a construction of the equilibrium state guaranteeing minimization of an entropy functional and dynamic control over the relaxation rate to enforce H-theorem through the streaming-collision process.

In [11] the authors pointed out that the classical moment-matchingbased approaches to the construction of the discrete equilibrium state fail to address entropy and therefor lack an H-theorem which can be the root of observed instabilities. To that end, they proposed a variational approach to the construction of the discrete equilibrium function minimizing a concave function of the discrete distribution functions under constraints on conserved moments. This construction led to the first entropic quadratic discrete equilibrium which was shown to be the minimizer of the Tsallis entropy with $\varphi=3 / 2$ [12]. Further research showed that in this variational approach to the construction

[^0]
[^0]:    * Corresponding author.

    E-mail address: ikarlin@ethz.ch (I.V. Karlin).
    https://doi.org/10.1016/j.compfluid.2023.105884
    Received 22 March 2023; Accepted 24 March 2023
    Available online 30 March 2023
    0045-7930/© 2023 The Authors. Published by Elsevier Ltd. This is an open access article under the CC BY license (http://creativecommons.org/licenses/by/4.0/).

---

<!-- Page 1 -->

of the discrete equilibrium two main challenges had to be considered: (a) the subspace of hydrodynamic parameters range accessible to the discrete equilibrium and (b) Errors in higher-order moments appearing at the Euler-Navier-Stokes level not explicitly covered by the initially set-out constraints [13]. For instance the Tsallis-type entropy-based construction led to the equilibria being only well-defined for a bounded domain of hydrodynamic parameters. In [13], the authors provided a possible answer to improve upon the previously-cited two points by noting that contrary to the continuous Boltzmann case where the entropy function is unique, and similar to ergodic Markov chains that support an infinite number of convex entropy functions here the form of the entropy function can be considered a degree of freedom of the construction. This ultimately led to the now-widely used form of the discrete entropy function which was presented initially in [14] as the perfect entropy functional. This choice of entropy functional was shown to recover the correct equilibrium pressure tensor in the limit of vanishing velocity and with deviations that were fourth-order in Mach number. While this discrete equilibrium was introduced without deriving explicit expressions, a closed form solution for the 1-D case was presented in [15] and extended to multiple dimensions and thermal flows in [16,17]. Another ingredient of the entropic lattice Boltzmann method guaranteeing non-linear stability is dynamic control over the maximum relaxation path-length set by imposing equality of entropy between post-relaxation and pre-collision populations in the limit of zero viscosity. This step boils down to finding non-trivial roots of the entropy estimate equation which can be solved either through numerical iterative solvers, or via approximations such as [18] or more advanced strategies such as the essentially entropic method, discussed in detail in this review. In the present contribution we will present a non-exhaustive review of the family of entropic lattice Boltzmann solvers and discuss challenges and most recent progress in that area. The review will start with Section 2, with a detailed discussion centered around entropic lattice Boltzmann for isothermal flows covering different aspects such as construction of the discrete equilibrium with its hydrodynamic/numerical properties and dynamic entropy-constrained relaxation and its variants such as the multiple relaxation time entropic collision operator. This will be followed by Section 3 covering the topic of entropic lattice Boltzmann models for non-ideal fluids with a brief overview of the theory of lattice Boltzmann for non-ideal fluids and different entropic formulations. The last section, i.e. Section 4, will discuss entropic models for compressible flows with a focus on higher order lattices covering more advanced optimization strategies such as the concept of shifted lattices.

## 2. Entropic lattice Boltzmann method for isothermal flows

### 2.1. Brief refresher on the lattice Boltzmann method

The lattice Boltzmann method can be considered to be a discrete solver for the Boltzmann equation with the BGK approximation [5] for the collision operator:
$\partial_{t} f(\boldsymbol{v}, \boldsymbol{r}, t)+\boldsymbol{v} \cdot \nabla f(\boldsymbol{v}, \boldsymbol{r}, t)=\frac{1}{t}\left[f^{\mathrm{eq}}\left(\rho, \boldsymbol{u}, k_{B} T / m\right)-f(\boldsymbol{v}, \boldsymbol{r}, t)\right]$,
where $f$ is the probability distribution function, $f^{\text {eq }}$ the equilibrium distribution function, $\boldsymbol{v}$ the particle speed, $\boldsymbol{u}$ the fluid velocity, $\rho$ the density and $T$ the temperature. $t$ is the relaxation time, $r$ and $t$ the coordinates in space and time, $k_{B}$ the Boltzmann constant and $m$ the molecular mass. The equilibrium distribution function, also known as the Maxwell-Boltzmann distribution function, is found as the state minimizing entropy under constraints on conserved moments of the distribution function, i.e. density, momentum and energy:
$f^{\text {eq }}\left(\rho, \boldsymbol{u}, k_{B} T / m\right)=\frac{\rho}{\left(2 \pi k_{B} T / m\right)^{d / 2}} \exp \left[-\frac{m\left(\boldsymbol{v}-\boldsymbol{u}\right)^{2}}{2 k_{B} T}\right]$.
Here $D$ denotes the physical dimension of the considered system. To arrive to a set of simple algebraic time-evolution equations two main ingredients are applied to Eq. (1): (a) Discretization of the space of particle velocities and (b) Discretization in physical space and time.

### 2.1.1. Discretization of particles speed space and discrete equilibrium state construction

In the discretization of the particles speed space, construction of the discrete equilibrium is the first and most important step. In the classical approach to the construction of the discrete equilibrium state, while a variety of realizations can be found in the literature, recovery of a set of moments of the equilibrium function needed to correctly recover the target macroscopic system of equations itself determined from a multi-scale perturbation analysis is the main paradigm. The number of constraints, i.e. number of moments needed, determines the minimum number of discrete velocities:
$\int_{V} \phi(\boldsymbol{v}) f^{\mathrm{eq}} d \boldsymbol{v}=\sum_{i=1}^{Q} \phi\left(\boldsymbol{c}_{i}\right) f_{i}^{\mathrm{eq}}$,
where $\phi$ is a vector of size $Q$, i.e. number of moments constraints, with elements corresponding to linearly independent polynomials of $\boldsymbol{v}$ and $\boldsymbol{c}_{i}$ the set of discrete particle velocities. In addition to constraints on the moments, properties such as symmetry and isotropy have to be considered when choosing the discrete velocities. For Lagrangian algorithms like the lattice Boltzmann method, to avoid interpolation operations for reconstruction of distribution functions on the grid points the discrete velocities lattice must be space-filling.

One systematic approach to the construction of a discrete system satisfying all of the above-listed constraints is the combined use of Hermite expansions and Gauss-Hermite quadrature rules. The latter quadrature allows to estimate integrals of the following form [19]:
$\int_{V} P(v) \exp \left(-v^{2}\right) d v \approx \sum_{i=1}^{Q} w_{i} f\left(c_{i}\right)$,
where $P(v)$ is a polynomial function of $v$. Designating the order of the polynomial by $M$, the above approximation becomes exact if $2 Q-1 \geq$ $M$ and the discrete velocities $c_{i}$ are set to the roots of the Hermite polynomial of order $Q$. Ideally, given the form of the Maxwell-Boltzmann equilibrium distribution function and the Gauss-Hermite quadrature one could directly obtain an integral of the form of Eq. (4) for moments of the equilibrium via a simple change of variable $\eta=\frac{1 / 2 \pi}{1 / 2 Q}$. However this would mean that the discrete velocities $c_{i}$ would now become function of local velocity and temperature, e.g. for a third-order quadrature $c_{i} \in\left\{u-\sqrt{3 k_{B} T / m}, u, u+\sqrt{3 k_{B} T / m}\right\}$. This is not convenient as in the context of a Lagrangian discretization in space, discrete distribution functions would not always propagate on-lattice. To overcome that issue the equilibrium distribution function is expanded around a reference state, i.e. static reference frame at a given reference temperature; One famous and widely used results of this construction approach is the second-order polynomial form of the equilibrium distribution function:
$f_{i}^{\mathrm{eq}}=w_{i} \rho\left(1+\frac{\boldsymbol{c}_{i} \cdot \boldsymbol{u}}{\zeta^{2}}+\frac{1}{\zeta^{2}}\left(\boldsymbol{c}_{i} \cdot \boldsymbol{u}\right)^{2}-\frac{1}{\zeta^{2}} \boldsymbol{u}^{2}\right)+\mathcal{O}\left(\boldsymbol{u}^{3}\right)$.
Here $\zeta$ is the sound speed at the lattice reference temperature and $w_{i}$ are weights stemming from the Gauss-Hermite quadrature.

### 2.1.2. Discretization in physical space and time

After discretization in phase-space one ends up with a set of $Q$ coupled $D+1$-dimensional hyperbolic partial differential equations (PDE):
$\partial_{t} f_{i}(\boldsymbol{r}, t)+\boldsymbol{c}_{i} \cdot \nabla f_{i}(\boldsymbol{r}, t)=\frac{1}{\tau}\left[f_{i}^{\mathrm{eq}}\left(\rho, \boldsymbol{u}, \zeta^{2}\right)-f_{i}(\boldsymbol{r}, t)\right]$.
This system of PDEs can be discretized using a variety of approaches such as the classical finite volumes or finite differences approximations. In the context of the lattice Boltzmann method each PDE is integrated along its characteristics leading to a set of coupled implicit algebraic equations:
$\int_{t}^{t+\delta t}\left(\partial_{t^{\prime}} f_{i}+\boldsymbol{c}_{i} \cdot \nabla f_{i}\right) d t^{\prime}=$

---

<!-- Page 2 -->

$$
\int_{t}^{t+\delta t} \frac{1}{\tau}\left[f_{i}^{\mathrm{eq}}\left(\rho, \boldsymbol{u}, \varsigma^{2}\right)-f_{i}\left(\boldsymbol{r}, t^{\prime}\right)\right] d t^{\prime}
$$

The advection speeds in each PDE being constant, the integration along characteristics yields an exact estimate of the advection term:
$\int_{t}^{t+\delta t}\left(\partial_{t^{\prime}} f_{i}+\boldsymbol{c}_{i} \cdot \nabla f_{i}\right) d t^{\prime}=f_{i}\left(\boldsymbol{r}+\boldsymbol{c}_{i} \delta t, t+\delta t\right)-f_{i}(\boldsymbol{r}, t)$,
while a lower order approximation is used to evaluate the collision term. To minimize dissipation a second-order scheme based on the trapezoidal rule is employed leading to:
$\int_{t}^{t+\delta t} \frac{1}{\tau}\left[f_{i}^{\mathrm{eq}}\left(\boldsymbol{r}, t^{\prime}\right)-f_{i}\left(\boldsymbol{r}, t^{\prime}\right)\right] d t^{\prime}=\frac{\delta t}{2 \tau}\left[f_{i}^{\mathrm{eq}}(\boldsymbol{r}, t)-f_{i}(\boldsymbol{r}, t)\right]$
$+\frac{\delta t}{2 \tau}\left[f_{i}^{\mathrm{eq}}\left(\boldsymbol{r}+\boldsymbol{c}_{i} \delta t, t+\delta t\right)-f_{i}\left(\boldsymbol{r}+\boldsymbol{c}_{i} \delta t, t+\delta t\right)\right]+\mathcal{O}\left(\delta t^{3}\right)$
which after a change of variables,
$\bar{f}_{i}=f_{i}-\frac{\delta t}{2 \tau}\left(f_{i}^{\mathrm{eq}}-f_{i}\right)$,
in combination with Eqs. (7) and (8) leads to the classical collisionstreaming equations:
$\bar{f}_{i}\left(\boldsymbol{r}+\boldsymbol{c}_{i} \delta t, t+\delta t\right)-\bar{f}_{i}(\boldsymbol{r}, t)=\frac{\delta t}{\bar{\tau}}\left[f_{i}^{\mathrm{eq}}(\boldsymbol{r}, t)-\bar{f}_{i}(\boldsymbol{r}, t)\right]$,
with
$\bar{\tau}=\tau+\delta t / 2$.
Note that for the sake of simplicity the overbar for the redefined distribution function will be omitted in the remainder of the article and the redefined relaxation coefficient noted as:
$\beta=\frac{\delta t}{2 \tau+\delta t}$.

### 2.1.3. Hydrodynamic limit: Macroscopic balance equations

At the macroscopic level, the leading order dynamics of this system of equations can be easily derived by subsequent application of a Taylor expansion around $(\boldsymbol{r}, t)$ and a multi-scale perturbation analysis based on a smallness parameter $\varepsilon \propto \delta t / T$ with $T$ the flow characteristic time. In the case of first-neighbor lattices, i.e. D1Q3, D2Q9 and D3Q27, the following macroscopic equations are recovered at order $\varepsilon+\varepsilon^{2}$ :

$$
\begin{aligned}
\partial_{t} \rho+\nabla \cdot \rho \boldsymbol{u} & =0 \\
\partial_{t}(\rho \boldsymbol{u})+\boldsymbol{\nabla} \cdot \rho \boldsymbol{u} \otimes \boldsymbol{u}+\boldsymbol{\nabla} \cdot \boldsymbol{T} & =0
\end{aligned}
$$

with:
$\boldsymbol{T}=P \boldsymbol{I}-\mu\left(\nabla \boldsymbol{u}+\nabla \boldsymbol{u}^{\dagger}-\frac{2}{D} \nabla \cdot \boldsymbol{u} \boldsymbol{I}\right)-\eta \boldsymbol{\nabla} \cdot \boldsymbol{u} \boldsymbol{I}$,
where $P$ is the pressure, $\mu$ the dynamic viscosity and $\eta$ the bulk viscosity. The dynamic viscosity is tied to the relaxation time via:
$\tau=\frac{\mu}{P}$,
while the recovered bulk viscosity is $\eta=\frac{2 \mu}{D}$ and the pressure $P=\rho \varsigma^{2}$ with the sound speed tied to the grid and time-step sizes, $\delta r$ and $\delta t$ as:
$\varsigma^{2}=\frac{\delta t^{2}}{3 \delta r^{2}}$,
meaning the reference temperature for a given choice of time-step and grid size can be found as:
$\frac{T}{m}=\frac{\delta t^{2}}{3 \delta r^{2} k_{B}}$.
With this system of equations, under diffusive scaling, i.e. $\delta t \propto \delta r^{2}$ and equivalent to $U^{\prime} / \varsigma \propto 1 / \delta t$ where $U^{\prime}$ is the characteristic convective velocity, in the limit of $\delta t \rightarrow 0$ one converges to the incompressible Navier-Stokes equations. Under convective scaling on the other hand, i.e. $\delta t \propto \delta r$ and equivalent to $U^{\prime} / \varsigma \propto \varepsilon s t$, one recovers the isothermal compressible Navier-Stokes equations instead.

### 2.2. Construction of entropic discrete equilibria

The first change of paradigm in the context of the entropic lattice Boltzmann method concerns the construction of the equilibrium state. In the context of the entropic lattice Boltzmann method the discrete equilibrium state is found as the minimizer of a Lyapunov functional $H$ :
$H=\sum_{i=1}^{Q} h_{i}\left(f_{i}\right)$
where $h_{i}$ are convex functions, under mass and momentum conservation constraints:

$$
\begin{aligned}
& \sum_{i=1}^{Q} f_{i}^{\mathrm{eq}}=\rho \\
& \sum_{i=1}^{Q} c_{i} f_{i}^{\mathrm{eq}}=\rho \boldsymbol{u}
\end{aligned}
$$

which is different from the polynomial construction which also imposes a constraint on the second-order moment of the equilibrium distribution function. Introducing the corresponding Lagrange multipliers, $\lambda_{0}, \lambda_{a}$ :
$\delta \sum_{i=1}^{Q}\left(h_{i}\left(f_{i}^{\mathrm{eq}}\right)-\lambda_{0} f_{i}^{\mathrm{eq}}-\lambda_{a} c_{i a} f_{i}^{\mathrm{eq}}\right)=0$,
yields
$h_{i}^{\prime}\left(f_{i}^{\mathrm{eq}}\right)=\lambda_{0}+\lambda_{a} c_{i a}$.
Defining the inverse of $h_{i}^{\prime}\left(f_{i}^{\mathrm{eq}}\right)$ as $\mu_{i}=\left[h_{i}^{\prime}\left(f_{i}^{\mathrm{eq}}\right)\right]^{-1}$, which must exist due to the convexity of $h_{i}$ the formal solution of the minimization problem reads:
$f_{i}^{\mathrm{eq}}=\mu_{i}\left(\lambda_{0}+\lambda_{a} c_{i a}\right)$.
Taking the $H$-function to be $[14,16]$ :
$H=\sum_{i=1}^{Q} f_{i} \ln \left(f_{i} / w_{i}\right)$,
and considering stencils with $3^{D}$ discrete velocities the corresponding equilibria can be expressed as:
$f_{i}^{\mathrm{eq}}=w_{i} \exp \left(\lambda_{0}\right) \prod_{a=1}^{D} \exp \left(c_{i a} \lambda_{a}\right)$.
Introducing the following changes of variables, $X=\exp \left(-\lambda_{0}\right)$ and $Z_{a}=\exp \left(\lambda_{a}\right)$ the equilibrium can be re-written as:
$f_{i}^{\mathrm{eq}}=w_{i} X^{-1} \prod_{a=1}^{D} Z_{a}^{c_{i a}}$.
Here $X$ and $Z_{a}$, i.e. the corresponding Lagrange multipliers, can be obtained by writing the constraints on the moments:

$$
\begin{aligned}
& \rho X=\sum_{i=1}^{Q} w_{i} \prod_{a=1}^{D} Z_{a}^{c_{i a}} \\
& \rho u_{p} X=\sum_{i=1}^{Q} w_{i} c_{i p} \prod_{a=1}^{D} Z_{a}^{c_{i a}}
\end{aligned}
$$

Solving this system of equations results in:

$$
\begin{aligned}
& Z_{a}=\frac{2 u_{a}+\sqrt{\left(u_{a} / \varsigma\right)^{2}+1}}{1-u_{a}} \\
& X^{-1}=\rho \prod_{a=1}^{D}\left(2-\sqrt{\left(u_{a} / \varsigma\right)^{2}+1}\right)
\end{aligned}
$$

---

<!-- Page 3 -->

![img-0.jpeg](img-0.jpeg)

**Fig. 1.** H-function of Eq. (33) shown with black dashed line as a function of *f*<sup>*</sup><sub>1</sub>. The red symbol represents the equilibrium state of Eq. (31) while the blue symbols show the equilibrium of Eq. (5).

leading to the following final expression for the discrete entropic equilibrium:

$$f_i^{(0)} = w_i \rho \prod_{a=1}^{D} \left( 2 - \sqrt{(u_a/\zeta)^2 + 1} \right) \left( \frac{2u_a + \sqrt{(u_a/\zeta)^2 + 1}}{1 - u_a} \right)^{1/a}.\tag{31}$$

where in 1-D *w*<sup>*</sup><sub>i</sub> ∈ {1/6, 2/3, 1/6}. In D-dimensional form, *w*<sup>*</sup><sub>i</sub> are the same as those obtained from the Gauss–Hermite quadrature for discrete velocities placed at the third-order Hermite polynomial roots.

To better illustrate the concept of entropy and equilibrium as its minimizer let us consider the simple case of a D1Q3 stencil under fixed density ρ = 1 and velocity *u*<sup>*</sup><sub>x</sub> = 0.2. Having fixed the values of the collisional invariants, only one free parameter is left, i.e. one of the populations:

$$f_0 = \rho - (f_1 + f_{-1}),\tag{32a}$$

$$f_{-1} = f_1 - \rho u_x. \tag{32b}$$

The discrete entropy of Eq. (25) can therefore be re-written as:

$$\begin{split}
H &= \frac{1}{2} \ln(3\rho) - 3f_1 + 0.5(3\rho u_x)(\rho - 2f_1 + \rho u_x) \\
&+ (f_1 - \rho u_x) \ln(6f_1 - 6\rho u_x) + f_1 \ln(6f_1).
\end{split} \tag{33}$$

The discrete entropy function for that case is shown in Fig. 1. This figure confirms the convexity of the discrete entropy function and that the entropic equilibrium does indeed minimize it. The second-order polynomial on the other hand, as observed in that plot does not guarantee minimization of the entropy functional.

### 2.3. Properties of discrete entropic equilibria

#### 2.3.1. Equilibrium pressure and sound speed

As observed in the previous section, to enforce the minimization of entropy, the strict condition on the second-order moment of the entropic equilibrium was lifted. Although the introduction of weights into the equilibrium distribution function allows to enforce constraints on the second-order moments for zero velocity one expects deviation for non-zero velocities. To that end let us first consider the diagonal components of the second-order moments tensor in the co-moving reference frame. We consider the co-moving reference frame here to focus on the recovered thermodynamic pressure. For the exact Maxwell–Boltzmann distribution at a reference pressure leading to a sound speed ζ it is readily computed as:

$$\overline{\Pi_{xx}^{200}} + \overline{\Pi_{yy}^{200}} = 2\rho \zeta^2. \tag{34}$$

For polynomial equilibria like the one shown in Eq. (5) the correct form of the thermodynamic pressure is recovered as the condition is explicitly taken into account in the construction. For the entropic equilibrium of Eq. (31) one recovers the following form:

$$\overline{\Pi_{xx}^{200}} + \overline{\Pi_{yy}^{200}} = 2\rho \zeta^2 + P_{xx}^* + P_{yy}^*,\tag{35}$$

with:

$$P_{xx}^* = -\rho \zeta^2 \left[ (u_x/\zeta)^2 - 2\sqrt{(u_a/\zeta)^2 + 1} + 2 \right]. \tag{36}$$

The deviations are illustrated in Fig. 2 via the following normalized variable:

$$\delta = \frac{P_{xx}^* + P_{yy}^*}{2\rho \zeta^2}. \tag{37}$$

Accordingly on observes that the entropic equilibrium pressure tensor admits deviation from the Maxwell–Boltzmann pressure tensor that are not isotropic. However, deviations are negligible for non-dimensional velocities as large as 0.5 which is well above the domain of validity of the weakly compressible flow assumption making them immaterial.

Going further in the analysis, a characteristics analysis of the Euler-level balance equations for conserved moments points to the existence of two non-symmetrical sound propagation speeds [20]:

$$\begin{aligned}
c_x^{\text{entropic,}*}+ &= \frac{u_x + \zeta \sqrt{2 \sqrt{(u_x/\zeta)^2 + 1 - 1}}}{\sqrt{(u_x/\zeta)^2 + 1}} - u_x,\\c_0^2 &u_x - \zeta \sqrt{2 \sqrt{(u_x/\zeta)^2 + 1 - 1}}} \\
c_x^{\text{entropic,}*}- &= \frac{u_x - \zeta \sqrt{2 \sqrt{(u_x/\zeta)^2 + 1 - 1}}}{\sqrt{(u_x/\zeta)^2 + 1}} - u_x.
\end{aligned} \tag{38b}$$

The propagation speed of these two eigen-modes are shown in Fig. 3 as a function of speed. Keeping in mind that the fastest eigen-modes in the system cannot propagate faster than the lattice links and only considering two *physical* eigen-modes, i.e. *c*<sup>*</sup><sub>x</sub> and *c*<sup>*</sup><sub>y</sub>, one arrives to the following condition on linear stability:

$$\max\{|u_x + c_x^*|, |u_x + c_y^^*|) \le \frac{\delta r}{\delta t},\tag{39}$$

which is satisfied unconditionally by the entropic equilibrium as the sound speed self-adjusts as a function of local velocity. For the polynomial equilibria, sound speed being constant, one recovers the following maximum convection velocity:

$$|u_x^{\text{max}}| = \frac{\delta r}{\delta t} - \zeta = 0.4226 \frac{\delta r}{\delta t},\tag{40}$$

which as will be shown in the next section through stability analyses, is indeed the maximum reachable velocity.

The slightly modified equilibrium pressure also affects dissipation at the Navier–Stokes level. The altered dissipation behavior is discussed in the next section.

#### 2.3.2. Viscous dissipation at Navier–Stokes level

To illustrate the impact of the entropic equilibrium on dissipation we will focus on the simple case of a 1-D system with three discrete velocities. Performing the classical Chapman–Enskog analysis, at order ε<sup>2</sup> the following momentum balance equation is recovered:

$$\delta_r^{(2)} \rho u_x + \delta_x (1 - \beta) \Pi_{xx}^{(1)} = 0,\tag{41}$$

where after some algebra the non-equilibrium stress tensor *Π*<sub>xx</sub><sup>(1)</sup> for any diagonal equilibrium pressure of the form:

$$P_{xx} = \rho \zeta^2 + P_{xx}^*,\tag{42}$$

results in:

$$\Pi_{xx}^{(1)} = -\frac{1}{2\beta} \left[ 2A\rho \zeta^2 \delta_x u_x + B\delta_x \rho \right],\tag{43}$$

where [20]

$$A = 1 - \frac{3u_x^2}{2\zeta^2} - \frac{3}{2\rho \zeta^2} u_x \delta_{u_x} P_{xx}^* - \frac{(\delta_{u_x} P_{xx}^*)^2}{2\rho^2 \zeta^2} - \frac{\delta_x P_{xx}^*}{2\zeta^2},\tag{44}$$

![img-0.jpeg](images/page3-img-0.jpeg.jpeg)



---

<!-- Page 4 -->

![img-1.jpeg](img-1.jpeg)

**Fig. 2.** Illustration of deviations in trace of pressure tensor for the entropic equilibrium. Left: Normalized deviation in 2-D. Right: Normalized error in 1-D.

![img-2.jpeg](img-2.jpeg)

**Fig. 3.** (Left) Sound speed for entropic and polynomial equilibria as a function of velocity *u<sub>x</sub>*. (Right) Comparison of the speed of fastest propagating eigen-modes: (blue dashed line) polynomial and (red line) entropic equilibria. *Source:* Plots reproduced from [20].

and

$$B = -3u_x \partial_x P_{xx}^* - \frac{\zeta^2}{\rho} \partial_{u_x} P_{xx}^* - \frac{\partial_{u_x} P_{xx}^* \partial_x P_{xx}^*}{\rho} - u_x^3. \tag{45}$$

For the case of the polynomial equilibrium one gets:

$$A = 1 - \frac{3u_x^2}{2\zeta^2}, \tag{46}$$

and:

$$B = -u_x^3, \tag{47}$$

which, neglecting the second term, indicates that the corresponding partial differential equation is only dissipative for:

$$\left| u_x \right| \le \frac{2\zeta}{3}. \tag{48}$$

In the case of the entropic equilibrium [20]:

$$A = 1 - \frac{3u_x^2}{2\zeta^2} + \frac{(u_x/\zeta)^2 + 3 (u_x/\zeta)^4 - 2 \sqrt{(u_x/\zeta)^2 + 1 + 2}}{2(u_x/\zeta)^2 + 2} \tag{49}$$

and:

$$B = 0. \tag{50}$$

Clearly for the entropic equilibrium, in the limit of *u<sub>x</sub> → δr/δt*, *A* goes to zero. The effective viscosities of Eqs. (46) and (49) are compared and shown in Fig. 4. Going even further and comparing the contribution from the second term, i.e. *B*, it is observed that the entropic equilibrium does not admit such a deviation.

### 2.3.3. Linear stability: Entropic vs polynomial equilibria

To systematically quantify the linear stability and spectral properties of the entropic equilibrium we consider the von Neumann approach. It aims at studying the time evolution of a perturbation *f<sub>i</sub><sup>a</sup>* that is injected into the linearized discrete time evolution equations. Briefly put, the perturbation is expanded as a combination of standing waves, whose propagation speed and attenuation rate is obtained. A positive attenuation rate will result in a growth of the signal at the corresponding wave-length and linear instability of the solver for the set of parameters considered, here *β*, Ma = ||*u*||/*ζ*. On the other hand, a scheme is linearly stable if the attenuation rate remains negative for all wave-numbers.

Furthermore, the spectral behavior and accuracy can be readily analyzed by comparing the spectral dispersions and dissipations to the theoretical modes obtained from the linearized Navier–Stokes equations. The Navier–Stokes theoretical modes for an isothermal flow can be expressed as [21]:

$$\omega^{\text{linear}} = u \cdot k - iv k^2, \tag{51a}$$

$$\omega^{\text{acoustic}} = (u \pm \zeta) \cdot k - i \left( \frac{D - 1}{D} v + \frac{\eta}{2\rho} \right) k^2, \tag{51b}$$

where *D* is the physical dimension of the system and *k* the wavenumber vector. As a consequence, the von Neumann stability analysis

![img-1.jpeg](images/page4-img-1.jpeg.jpeg)



![img-2.jpeg](images/page4-img-2.jpeg.jpeg)



---

<!-- Page 5 -->

![img-3.jpeg](img-3.jpeg)

Fig. 4. Comparison of the effective bulk viscosity: (blue) polynomial and (red) entropic equilibria.

*Source:* Plot reproduced from [\[20\]](#page-11-13).

can be used to evaluate the spectral behavior and linear stability domain of a solver for a given set of parameters. As such it can be perceived as a tool to objectively evaluate the stabilization properties of different collision models, on the basis of necessary conditions. The latter comes from the fact that the analysis relies on a linearization step and as such gets the sufficient condition for stability only under the linear regime assumption (small amplitude perturbations). It has been widely used in the past to evaluate the stability properties of the lattice Boltzmann method. Details on this analysis approach can be found in [\[7,21\]](#page-11-14),[28], among other sources.

The equilibrium state is one of the most important components of a kinetic scheme and controls, for the most part, the leading-order dynamics of the system, but also the behavior of higher-order terms. The effects of the equilibrium in the limit of vanishing wave-number were studied in the previous sections. Here using the von Neumann formalism we look at the effect of the equilibrium on the linear stability domain. The eigenvalue problem of the von Neumann equations is solved for different values of non-dimensional viscosities, over the entire wave-number space, i.e. \( k_x \) and \( k_y \) with a resolution of 100 points in each direction. The highest Mach number resulting in negative dissipation rates over all wave-numbers is retained as the linear stability limit. These limits, as obtained for different polynomial equilibria, are shown in Fig. 5. Looking at those results a number of points are worth noting: For all of these equilibria, regardless of the value of the non-dimensional viscosity, the maximum stable Mach number never goes beyond \( Ma = \sqrt{3} - 1 \approx 0.732 \). This confirms the observation in Eq. (40):

$$
Ma_{p}^{max} = \frac{|a|_{p}^{max}|}{\zeta} = \sqrt{3} - 1. \tag{52}
$$

Furthermore, while the addition of third-order components appears not to have a large effect on the stability domain, moving to the product form extends it.

Apart from extending the linear stability domain, the addition of the fourth-order component results in more isotropic behavior especially for small values of the non-dimensional viscosity. The directional stability domains obtained with different orders of the EDF are shown in Fig. 6. For the entropic equilibrium on the other hand, the scheme is found to come with unconditional linear stability for all values of the Mach number supported by the stencil, i.e. \( Ma = \sqrt{3} \), even for vanishing viscosities. The stability domain is shown in Fig. 7. This in turn confirms the effectiveness of the discrete equilibrium construction approach in guaranteeing linear stability (by enforcing a discrete H-theorem).

### 2.4. Entropy-constrained relaxation towards the equilibrium attractor

Apart from the form of the discrete equilibrium state, the entropic lattice Boltzmann method differs from the LBGK in the collision operator. The idea is to implement the discrete-time analog of Boltzmann's H-theorem. To illustrate that, let us re-write the discrete time-evolution equation as:

$$
f_i(\mathbf{r} + \mathbf{c}_i \delta t, t + \delta t) - f_i(\mathbf{r}, t) = \alpha \beta [f_i^{\text{eq}}(\mathbf{r}, t) - f_i(\mathbf{r}, t)],\tag{53}
$$

where we have introduced a parameters \( \alpha \), which we shall refer to as the over-relaxation parameter and is operationally available as the positive root of the following entropy condition:

$$
H(f + \alpha (f^{\text{eq}} - f)) = H(f). \tag{54}
$$

This condition on the over-relaxation parameter guarantees that even in the limit of zero viscosity, i.e. \( \beta \rightarrow 1 \), the entropy does not increase. In practice, this is achieved by solving Eq. (54) equivalent to finding the intersection of the line starting from \( f \) moving in the direction of \( f^{\text{eq}} \) with the iso-entropy surface \( H = H(f) \) which gives the convexity of the entropy functional is guaranteed to exist. This process is illustrated in Fig. 8. Here note a few implications of the entropy condition of Eq. (54):

- **Over-relaxation**: Thanks to convexity of the entropy functional, solution to Eq. (54) always leads to over-relaxation, i.e. \( \alpha > 1 \);
- **Mirror symmetry**: Let \( f \) be a population vector, and \( f(\alpha) = f + \alpha (f^{\text{eq}} - f) \) its entropic mirror state, with the same value of the entropy, \( H(f(\alpha)) = H(f) \). If the entropy estimate is applied to \( f(\alpha) \) instead of \( f \), then the initial state is recovered in the form:

$$
f = f(\alpha) + \alpha^t (f^{\text{eq}} - f(\alpha)) \tag{55}
$$

with another over-relaxation \( \alpha^t > 1 \) which satisfies a duality relation:

$$
\alpha^t \alpha = \alpha^t + \alpha. \tag{56}
$$

- **Resolved simulation limit**: Whenever the simulation is resolved (populations stay close to the local equilibrium), the maximal over-relaxation parameter \( \alpha \) becomes fixed automatically to the value \( \alpha = 2 \) (and so is also the mirror value, \( \alpha^t = 2 \)). Then the entropic relaxation converges to the LBGK over-relaxation. This is a direct implication of the built-in discrete-time \( H \)-theorem. A resolved simulation, on the kinetic level, is characterized by the fact that all populations are asymptotically close to the local equilibrium. Then, the entropy function can be well-represented by a second-order approximation; For:

$$
\delta f = f - f^{\text{eq}}, |\delta f/f^{\text{eq}}| \ll 1, \tag{57}
$$

one gets:

$$
H(f) \approx H^{\text{eq}} + \frac{1}{2} \sum_{i=1}^O \delta f_i^2/f_i^{\text{eq}}. \tag{58}
$$

The levels of the entropy are then asymptotically close to the levels of the above quadratic form. Under such condition the entropy estimate of Eq. (54) results in \( \alpha = 2 \).

Whenever the grid is coarsened, levels of the entropy functions beyond the quadratic approximations start to be explored by populations at some sites. The entropy estimate then delivers \( \alpha < 2 \) for some of them, and \( \alpha > 2 \) for some others. In this situation, stipulating the LBGK updates with fixed \( \alpha = 2 \) on the sites requiring \( \alpha < 2 \), amounts to violating the \( H \)-theorem. The domain of validity of the quadratic approximation and the \( H \)-theorem does not hold. Of course, no such violation would occur at sites with \( \alpha > 2 \); Nevertheless, it is clear that lack of knowledge

![img-3.jpeg](images/page5-img-3.jpeg.jpeg)



---

<!-- Page 6 -->

![img-4.jpeg](img-4.jpeg)

**Fig. 5.** Linear stability domains of SRT collision operator with (from left to right) second order, third order and product form equilibria.
*Source:* Reproduced from [28].

![img-5.jpeg](img-5.jpeg)

**Fig. 6.** Illustration of anisotropy of linear stability domains for EDFs of orders (from left to right) two, three and four, and for seven different non-dimensional kinematic viscosities, i.e. (—)5 × 10<sup>−3</sup>, (—)1 × 10<sup>−3</sup>, (—)5 × 10<sup>−3</sup>, (—)0.01, (—)0.05, (—)0.1, (—)0.5.
*Source:* Reproduced from [28].

![img-6.jpeg](img-6.jpeg)

**Fig. 7.** Illustration of linear stability domain for the entropic EDF for seven different non-dimensional kinematic viscosities, i.e. (—)5 × 10<sup>−3</sup>, (—)1 × 10<sup>−3</sup>, (—)5 × 10<sup>−3</sup>, (—)0.01, (—)0.05, (—)0.1, (—)0.5.
*Source:* Reproduced from [26,28].

of this repartition between "forgiving" (α > 2) and "unforgiving" sites (α < 2) is permanently exposed to instabilities. Entropic LBM exploits the self-adaptive mechanism of stabilization by choosing automatically the over-relaxation α at each node, which guarantees the H-theorem at all sites and all discrete time-steps. When the grid is coarsened, over-relaxation α becomes "smeared" in an interval, [α<sub>min</sub>, α<sub>max</sub>], with 1 < α<sub>min</sub> < 2, and α<sub>max</sub> > 2. The self-adapted over-relaxation set up by Eq. (54), results in two oppositely directed effects: If α < 2, the relaxation parameter αβ is less than the corresponding standard LBGK parameter 2β (set by the viscosity), and hence the entropic relaxation will tend to smoothen any flow perturbation. On the other hand, if α > 2, the flow perturbation is enhanced.

### 2.5. Subgrid-scale model interpretation of the entropic stabilizer

The entropic stabilizer, as introduced in the previous section modifies the effective viscosity locally so as to ensure decay of the discrete H-function. Considering that the modified viscosity is introduced at the level of the Navier–Stokes viscosity, and not under the form of some hyper-viscosity, it is clear that in regions where the stabilizer is active a form of subgrid-scale model is applied. This observation was first made numerically in [29] where the authors showed that the entropic model allows for stable under-resolved simulations correctly capturing dynamics of large-scale eddy structures, opening the door for an interpretation of the entropic model as a class of large eddy simulation (LES) tool. To explicitly evaluate the effect of the entropic stabilizer on effective viscosity a closed-form expression for the maximum pathlength α is needed. Later on in [18], to minimize the computational effort of solving the implicit system of Eq. (54) the authors derived a closed-form approximation for α assuming distributions close to the equilibrium distribution. A Taylor expansion of the discrete H-function, i.e. Eq. (25), near equilibrium led to the following expression:

$$
\alpha = 2 - \frac{4a_2}{a_1} + \frac{16a_2^2}{a_1^2} - \frac{8a_3}{a_1} + \frac{80a_2a_3}{a_1^2} - \frac{80a_3^3}{a_1^3} - \frac{16a_4}{a_1}, \tag{59}
$$

with:

$$
a_n = \frac{(-1)^{n-1}}{n(n+1)} \sum_{i=1}^{N} \frac{(f_i - f_i^{\text{eq}})^{n+1}}{f_i^n}. \tag{60}
$$

Writing the distribution function as \( f_i = f_i^{\text{eq}} + ε f_i^{(1)} \) and using the Chapman–Enskog expansion to evaluate the first-order non-equilibrium part of the distribution function:

$$
f_i^{(1)} = -\frac{\rho u_i}{2\beta} H_2 : \left(\nabla u + \nabla u^1\right). \tag{61}
$$

![img-4.jpeg](images/page6-img-4.jpeg.jpeg)



![img-5.jpeg](images/page6-img-5.jpeg.jpeg)



![img-6.jpeg](images/page6-img-6.jpeg.jpeg)



---

<!-- Page 7 -->

![img-7.jpeg](img-7.jpeg)

**Fig. 8.** Illustration of entropic relaxation process.

where *H*<sup>2</sup> is the second-order Hermite polynomials tensor, the leading-order dependence on the strain tensor is obvious. In [30] the authors carried out this development and found the following leading-order expression for the maximum path-length:

$$
\alpha = 2 + \frac{2\tilde{\mu} t}{3\tilde{\beta}} \frac{S_{ij} S_{ik} S_{jk}}{S_{lm} S_{lm}}, \tag{62}
$$

where *S* = (*P**I* − *T*)/2*μ*, which as observed in [31,32] is similar to a Smagorinsky subgrid-scale model for the turbulent viscosity ν<sub>T</sub> [33], i.e.

$$
\nu_T = \left(C_t \Delta\right)^2 \sqrt{S_{ij} S_{ji}}, \tag{63}
$$

where *C*<sub>*t*</sub> is the Smagorinsky constant. At the difference of the Smagorinsky, or any of the other well-known models, the turbulent viscosity recovered by the entropic model is not positive definite and as such can admit, in certain situations, energy transfer from under-resolved small scales to larger resolved scales [32,34]. The effect of this specific form of the turbulent viscosity was checked in detail in [34] where the authors compared simulations using the Smagorinsky model, the entropic model with the entropic stabilizer and lattice BGK model with an effective viscosity set to the form recovered in Eq. (62). It was observed that while the last two solvers led to similar results, the entropic subgrid closure led to better agreement with direct numerical simulation data as compared to the positive-definite Smagorinsky closure.

### 2.6. An efficient approach to compute the over-relaxation path-length: Essentially entropic formulation

Typically, the *H* theorem in entropic LBM is numerically enforced by iteratively solving the nonlinear equation (54) with a combination of bisection and Newton–Raphson methods. In Refs. [32,35] it was demonstrated that an analytic solution can be obtained by assuming a natural criterion of negative entropy change, thereby reducing the problem to solving the inequality

$$
\Delta H \equiv H(f + \alpha(f^{\text{eq}} - f)) - H(f) \leq 0. \tag{64}
$$

As any inequality accepts multiple solutions, we search for the solution corresponding to the maximal path length α such that Δ*H* → 0. Similar to the case with entropic LBM, the analytic solution should also reduce to LBGK (α = 2) close to equilibrium.

The inequality Eq. (64) is solved by first recasting it (using Eq. (25) and the entropic equilibrium from Eq. (27)) to

$$
\begin{split}
\Delta H & \equiv \sum_{a=1}^{D} f_i(1 + \alpha \tilde{z}_i) \log(1 + \alpha \tilde{z}_i) \\
& - \alpha \beta \sum_{a=1}^{D} f_i(1 + z_i) \log(1 + z_i) \leq 0,
\end{split} \tag{65}
$$

where *z*<sub>*i*</sub> = *f*<sup>eq</sup><sub>*i*</sub>/*f*<sub>*i*</sub> − 1 is the dimensionless departure from the equilibrium and *z̃*<sub>*i*</sub> = αβ*z*<sub>*i*</sub>. This is followed by creating a framework for constructing Padé approximants on logarithms via quadrature on appropriate convex function. For example, using only the convexity arguments one can construct both upper and lower bounds on any convex function via

$$
F\left(\frac{a + b}{2}\right) \leq \frac{1}{b - a} \int_0^b F(x) \, dx \leq \frac{1}{2} [F(a) + F(b)],
$$

where *a*, *b* are positive real numbers with *b* > *a*. The upper bound is the integral approximated via the trapezoid rule and the lower bound is the integral approximated as the value of the function at the midpoint rule. Integrating the convex function *F*(*z*) = 1/(1 + *z*) from 0 to *z* one obtains the bounds on log(1 + *x*) as

$$
\frac{2x}{2 + x} \leq \log(1 + x) \leq \frac{y(2 + y)}{2(1 + y)}. \tag{68}
$$

The above inequality provides loose bounds on log(1 + *x*), however, one can employ higher Gauss–Legendre and Newton–Cotes quadratures to construct tighter upper and lower bounds on log(1 + *x*).

Exploiting the above bounds on logarithm and neglecting some negative-definite terms that contribute to entropy production, one can derive a quadratic in α. The computationally inexpensive positive root of the quadratic is readily available as 2*c*/(b + √b<sup>2</sup> − 4*ac*) with

$$
a = \sum_{a=1}^{D} f_{i, -} \frac{z_i^2}{2}, \quad b = \sum_{a=1}^{D} f_i \frac{z_i^2}{2}, \quad c = \sum_{a=1}^{D} f_i \frac{2z_i^2}{2 + z_i}. \tag{69}
$$

where *f*<sub>*i*,−</sub>*, *z*<sub>*i*,−</sub> are the *f*<sub>*i*</sub>, *z*<sub>*i*</sub> for which *z*<sub>*i*</sub> < 0. In practice, the bounds on the logarithm from the trapezoid and midpoint rule lead to a diffusive numerical scheme, thus, one needs to employ bounds from higher quadrature rules.

### 2.7. Entropy-controlled modulation of relaxation rates of specific eigenmodes: The multiple-relaxation time entropic model

As a more recent iteration of entropic LBM proposed first in [36], the authors proposed a new approach to construct an optimal mirror state as the maximizer of the discrete entropy function under the constraint of over-relaxation of the hydrodynamic stresses. To that end it is noted that all discrete populations can be re-written as a combination of their moments:

$$
f_i = k_i + s_i + h_i, \tag{70}
$$

where conserved moments under collision are grouped into *k*<sub>*i*</sub>, moments of order two controlling the stress tensor in *s*<sub>*i*</sub> and all higher order moments not appearing in the previous terms in *h*<sub>*i*</sub>. Using that representation a modified mirror state can be defined as:

$$
f_i^{\text{mirr}} = k_i + \left[ 2s_i^{\text{eq}} - s_i \right] + \left[ (1 - \gamma)h_i + \gamma h_i^{\text{eq}} \right], \tag{71}
$$

where γ is a parameter to be determined later. It is clear that for any value of γ this construction of the mirror state recovers the Navier–Stokes equation with the correct viscosity.

![img-7.jpeg](images/page7-img-7.jpeg.jpeg)



---

<!-- Page 8 -->

![img-8.jpeg](img-8.jpeg)

Fig. 9. Illustration of applications of the entropic multiple relaxation time model: (left) Flow over the SD7003 at Re = 6 × 10<sup>4</sup> and (right) Cold flow in valve–piston assembly. Source: Images reproduced from [39,40].

The parameter γ, contrary to classical multiple-relaxation-time approaches is determined by maximizing the entropy of the post-collision state ƒ* which translates into solving the following equation:

$$\sum_{i=1}^{N} \Delta h_i \ln \left(1 + \frac{(1 - \beta\gamma)\Delta h_i - (2\beta - 1)\Delta s_i}{f_i^{\text{eq}}} \right) = 0. \tag{72}$$

An explicit form of this constraint can be derived using a Taylor expansion of $\Delta h_i/f_i^{\text{eq}}$ and $\Delta s_i/f_i^{\text{eq}}$, i.e. $\ln(1 + x) = x - x^2/2 + \dots$. Keeping only first order contributions:

$$\ln \left(1 + \frac{(1 - \beta\gamma)\Delta h_i - (2\beta - 1)\Delta s_i}{f_i^{\text{eq}}} \right) \approx \left(1 - \beta\gamma\right)\frac{\Delta h_i}{f_i^{\text{eq}}} - (2\beta - 1)\frac{\Delta s_i}{f_i^{\text{eq}}}, \tag{73}$$

which leads to the following condition on γ:

$$\gamma = \frac{1}{\beta} - \left(2 - \frac{1}{\beta}\right)\frac{\langle\Delta h|\Delta s\rangle}{\langle\Delta h|\Delta h\rangle}, \tag{74}$$

where $\langle X|Y\rangle$ is the entropic inner product defined as [37]:

$$\langle X|Y\rangle = \sum_{i=1}^{N} \frac{X_iY_i}{f_i^{\text{eq}}}.\tag{75}$$

An extension of this approach to arbitrary number of ghost moments stabilizer γ<sub>m</sub> was proposed in [37]. In that case the mirror state is defined as:

$$f_i^{\text{mirr}} = k_i + \left[2s_i^{\text{eq}} - s_i\right] + \sum_{m=1}^{K} \left[(1 - \gamma_m)h_{mi} + \gamma_mh_{mi}^{\text{eq}}\right], \tag{76}$$

and the stabilizer γ<sub>m</sub> is computed as [37]:

$$\gamma_m = \frac{1}{\beta} - \left(2 - \frac{1}{\beta}\right)\sum_{n=1}^{K} C_{mn}^{-1}\langle\Delta s|\Delta h_n\rangle,\tag{77}$$

where the correlation matrix C is defined as $C_{mn} = \langle\Delta h_m|\Delta h_n\rangle$. A similar approach extending the original two-relaxation time entropic formulation to independent relaxation rates for all modes was also proposed and discussed in [38].

The multiple relaxation time entropic model has been successfully used to model a wide range of highly turbulent flows with complex geometries, as illustrated in Fig. 9.

### 3. Extension to non-ideal fluids

In the context of the present contribution, by non-ideal fluids, we refer to the diffuse interface van der Waals fluid leading to the Korteweg stress tensor which is obtained by constrained minimization problem of the van der Waals fluid free energy subject to global conservation of mass leading to [41],

$$T_K = \nabla \otimes \frac{\partial\mathcal{L}}{\partial (\nabla \rho)} - \mathcal{L} \mathbf{I},\tag{78}$$

where $\mathbf{I}$ is unit tensor and $\mathcal{L}$ is the Lagrange function,

$$\mathcal{L} = \mathcal{A} + \frac{1}{2} \kappa |\nabla \rho|^2 - \lambda \rho,\tag{79}$$

$\mathcal{A}$ the bulk free energy, and $\lambda$, corresponding to the chemical potential, is,

$$\lambda = \frac{\partial \mathcal{A}}{\partial \rho} - \kappa \nabla^2 \rho.\tag{80}$$

Here $\kappa$ is the capillary coefficient. The Korteweg's stress tensor is then obtained as [42]:

$$T_K = \left(P - \kappa \rho \nabla^2 \rho - \frac{1}{2} \kappa |\nabla \rho|^2\right) \mathbf{I} + \kappa \nabla \rho \otimes \nabla \rho,\tag{81}$$

where

$$P = \rho \frac{\partial \mathcal{A}}{\partial \rho} - \mathcal{A},\tag{82}$$

is the thermodynamic pressure, or equation of state, different from the ideal gas equation of state. Including the non-ideal stress tensor into the usual balance equations the target hydrodynamic system is:

$$\partial_t \rho + \nabla \cdot \rho \mathbf{u} = 0,\tag{83}$$

$$\partial_t \rho \mathbf{u} + \nabla \cdot \rho \mathbf{u} \otimes \mathbf{u} + \nabla \cdot \mathbf{T} = 0,\tag{84}$$

where $\mathbf{u}$ is the fluid velocity and the stress tensor $\mathbf{T}$ is

$$T = T_K + T_{\text{NS}}.\tag{85}$$

The Navier–Stokes viscous stress tensor reads,

$$T_{\text{NS}} = -2\mu S - \eta (\nabla \cdot \mathbf{u}) \mathbf{I},\tag{86}$$

where $S$ is the trace-free rate-of-strain tensor,

$$S = \frac{1}{2} |\nabla \mathbf{u} + \nabla \mathbf{u} - \frac{2}{D} (\nabla \cdot \mathbf{u}) \mathbf{I} \tag{87}$$

and $\mu$ and $\eta$ are the dynamic and the bulk viscosity, respectively.

Eqs. (83) together with Eqs. (81), (85), (86) and (87) describe the dynamics of the van der Waals fluid targeted by non-ideal lattice Boltzmann models in the hydrodynamic limit. This system of equation can also be recovered from the kinetic theory:

$$\partial_t f + \mathbf{v} \cdot \nabla f = -\frac{1}{\tau} \left(f - f^{\text{eq}}\right) - \frac{1}{\rho} \frac{\partial f^{\text{eq}}}{\partial \mathbf{u}} \cdot \mathbf{F}_{\text{nloc}},\tag{88}$$

with:

$$\mathbf{F}_{\text{nloc}} = \iiint \nabla V \left( |\mathbf{r} - \mathbf{r}_1 | \right) f_2(\mathbf{r}, \mathbf{v}, \mathbf{r}_1, \mathbf{v}_1, t) d\mathbf{v}_1 d\mathbf{r}_1 d\mathbf{v},\tag{89}$$

where $f_2$ is the two-particle distribution function and $V$ is the interaction potential. Using the Enskog–Vlasov approximation and taking into account the form-invariance of the hydrodynamic equations with respect to pressure partition between the BGK collision term and external body force a family of models taylored to lattice Boltzmann models is recovered:

$$\partial_t f + \mathbf{v} \cdot \nabla f = -\frac{1}{\tau} \left(f - f^{\text{eq}}\right) - \frac{1}{\rho} \frac{\partial f^{\text{eq}}}{\partial \mathbf{u}} \cdot \left[\nabla \left(P - P_0\right) - \kappa \rho \nabla \nabla^2 \rho\right].\tag{90}$$

Here $P_0$ is the pressure carried by the BGK collision operator, a tunable parameter usually set to $P_0 = \rho \varsigma^2$. For detailed discussion on that topic we refer readers to [43]. In the next subsections we will discuss two main approaches to introduce the new stress tensor into LBM and different entropic realizations proposed in the literature.

### 3.1. ELBM for non-ideal fluids

#### 3.1.1. Entropic free pressure approach

For the LB solver to recover the target system of macroscopic equations, the second-order moment of the discrete equilibrium function must recover the Korteweg stress tensor. One way to recover this tensor is to introduce it into the equilibrium pressure [44]. To that end using

![img-8.jpeg](images/page8-img-8.jpeg.jpeg)



---

<!-- Page 9 -->

the entropic construction, for a D1Q3 lattice, introduction of this new pressure boils down to:
$\sum_{i=1}^{Q} c_{i}^{2} \mu_{i}\left(\lambda_{0}+\lambda_{x} c_{i}\right)-\rho u_{x}^{2}-P-\kappa \rho \frac{d^{2} \rho}{d x^{2}}+\frac{\kappa}{2}\left(\frac{d \rho}{d x}\right)^{2}=0$.
Contrary to the ideal gas case, here due to highly non-linear nature of the system it will be extremely difficult to find a closed-form solution. To simplify the process, the pressure is separated into two contributions as:
$P+\kappa \rho \frac{d^{2} \rho}{d x^{2}}-\frac{\kappa}{2}\left(\frac{d \rho}{d x}\right)^{2}=\rho \varsigma^{2}+K$,
and it is assumed that the non-ideal contribution, $K$ is a constant. This assumption is however mathematically untrue but serves as an approximation and allows us to proceed with the solution to find an entropy function. After some algebra one gets:
$\lambda_{0}=4 \exp (X-1)-K$
$\lambda_{x}=\exp (X-1)+\frac{1}{2} K$
which after inversion lead to the following entropy functions:
$h_{0}=(X+K / 4) \ln \left(\frac{X+K}{4}\right)$,
$h_{x}=(X-K / 2) \ln (X-K / 2)$.
Finally expanding the Lagrange multipliers in powers of $u_{x}$ up to order two one gets:
$f_{0}^{\mathrm{eq}}=w_{0} \rho\left(1-\frac{u_{x}^{2}}{2 \varsigma^{2}}\right)-K$
$f_{ \pm}^{\mathrm{eq}}=w_{ \pm} \rho\left(1+\frac{u_{x}^{2}}{\varsigma^{2}}+\frac{c_{ \pm} u_{x}}{\varsigma^{2}}\right)+\frac{K}{2}$.
Compared to the pressure-based construction of Swift et al. [45] it was observed that the present entropic construction allowed for lower kinematic viscosities [44]. It was also observed that the use of larger lattices, mainly to remove the Galilean-variant error in the diagonal third-order moments of the equilibrium allowed to reach density ratios of the order of 50 . Note that the now first-order in Mach Galileanvariant error can also be removed via an appropriate correction term. Please see for further discussions $[43,46]$.

### 3.1.2. Entropic free energy forcing scheme

The forcing approach for non-ideal fluids consists in using the classical discrete isothermal BGK collision operator with corresponding equilibrium distribution function with an ideal equation of state and accounting for deviations from that reference pressure, $\rho \varsigma^{2}$, via an external body force defined as [47]:
$\boldsymbol{F}=\boldsymbol{\nabla} \cdot\left(P-\rho \varsigma^{2}\right) \boldsymbol{I}-\kappa \rho \boldsymbol{\nabla} \boldsymbol{\nabla}^{2} \rho$.
changing the evolution equation into:
$f_{i}\left(\boldsymbol{r}+c_{i} \delta t, t+\delta t\right)-f_{i}(\boldsymbol{r}, t)=\alpha \beta\left[f_{i}^{\mathrm{eq}}(\boldsymbol{r}, t)-f_{i}(\boldsymbol{r}, t)\right]+\mathscr{G}_{i}(\boldsymbol{F})$,
where the last term on the RHS is a general term representing the implementation of the body force in the lattice Boltzmann solver. Here for the sake of readability we will assume the exact difference approximation, i.e.
$\mathscr{G}_{i}(\boldsymbol{F})=f_{i}^{\mathrm{eq}}\left(\rho, \boldsymbol{u}+\frac{\boldsymbol{F} \delta t}{\rho}\right)-f_{i}^{\mathrm{eq}}(\rho, \boldsymbol{u})$.
A comprehensive review of different body force realizations can be found in [43]. A simple multiscale analysis shows that under the condition that:
$\sum_{i=1}^{Q} \mathscr{G}_{i}(\boldsymbol{F})=0$,
$\sum_{i=1}^{Q} c_{i} \mathscr{G}_{i}(\boldsymbol{F})=\boldsymbol{F}$,
one recovers the Euler equation with the Korteweg stress tensor.
Note that in the presence of external body force the condition on the maximum path-length $\alpha$ of Eq. (54) changes into [47,48]:
$H\left(f^{\prime}+\alpha\left(f^{\mathrm{eq}}\left(\rho, \boldsymbol{u}+\frac{\boldsymbol{F} \delta t}{\rho}\right)-f^{\prime}\right)\right)=H\left(f^{\prime}\right)$
with:
$f_{i}^{\prime}=f_{i}+\left[f_{i}^{\mathrm{eq}}\left(\rho, \boldsymbol{u}+\frac{\boldsymbol{F} \delta t}{\rho}\right)-f_{i}^{\mathrm{eq}}(\rho, \boldsymbol{u})\right]$.
The first realization of the multiple relaxation time entropic model for non-ideal fluids was proposed in [49]. The authors used the exact difference method to introduce the external body force and redefined the entropic dot product for the evaluation of the stabilizer $\gamma$ so as to take into account the external body force:
$\langle X \mid Y\rangle=\sum_{i=1}^{Q} \frac{X_{i} Y_{i}}{f_{i}^{\mathrm{eq}}\left(\rho, \boldsymbol{u}+\boldsymbol{F} \delta t / \rho\right)}$.
A similar route was undertaken later on in [50] using a similar approach, i.e. exact difference forcing and raw moments-based multiple relaxation time entropic collision operator. The multiple relaxation time entropic model was also applied to the non-ideal fluid model of [46] in [51]. At the difference of previous realizations the central Hermite polynomials were used as the projection space for the relaxation process. This specific choice of moments space was based on the orthogonality of the central Hermite polynomials in the conoving reference frame with respect to the weighted dot product which is equivalent to the non-shifted entropic dot product. As such it was expected that the stabilizer $\gamma$ would force ghost moments to equilibrium everywhere except for interfacial regions where the body force is pronounced. This was readily confirmed by simulations of a drop impact on a thin liquid film, illustrated in Fig. 10. A comparative study of the stability of the classical SRT and entropic multiple relaxation time collision operators in [51] showed that the latter allowed for pronounced extension of the domain compared to the former. Results based on a 2-D case involving impact of a droplet onto a liquid film are shown in Fig. 11. Furthermore, it was observed that the use of the proposed realization of the entropic multiple relaxation time model allowed for faster convergence of spurious currents to their steady values as shown in Fig. 12 and, in general, reduced spurious currents at curved interfaces especially at lower viscosities, illustrated in Fig. 13, opening the door for higher Reynolds number simulations. Based on these observations, the entropic multiple relaxation time model has been used in a number of publications to model high Reynolds number multiphase flows, see for instance $[52,53]$.

## 4. Entropic lattice Boltzmann for compressible flows

In this section we consider lattice Boltzmann models for compressible flows of the following form [54]:

$$
\begin{aligned}
& f_{i}\left(\boldsymbol{r}+c_{i} \delta t, t+\delta t\right)-f_{i}(\boldsymbol{r}, t)=2 \beta_{1}\left[f_{i}^{\mathrm{eq}}(\boldsymbol{r}, t)-f_{i}(\boldsymbol{r}, t)\right] \\
& +2\left(\beta_{1}-\beta_{2}\right)\left[f_{i}^{\prime}(\boldsymbol{r}, t)-f_{i}^{\mathrm{eq}}(\boldsymbol{r}, t)\right] \\
& g_{i}\left(\boldsymbol{r}+c_{i} \delta t, t+\delta t\right)-f_{i}(\boldsymbol{r}, t)=2 \beta_{1}\left[g_{i}^{\mathrm{eq}}(\boldsymbol{r}, t)-g_{i}(\boldsymbol{r}, t)\right] \\
& +2\left(\beta_{1}-\beta_{2}\right)\left[g_{i}^{\prime}(\boldsymbol{r}, t)-g_{i}^{\mathrm{eq}}(\boldsymbol{r}, t)\right]
\end{aligned}
$$

where the second population, $g_{i}$ is introduced to transport energy contributions from non-translational degrees of freedom and:
$g_{i}^{\mathrm{eq}}=\left(2 C_{c}-D\right) T f_{i}^{\mathrm{eq}}$.
where $C_{c}$ is the specific heat capacity at constant volume. The collisional invariants are computed as:

$$
\rho=\sum_{i=1}^{Q} f_{i}
$$

---

<!-- Page 10 -->

![img-9.jpeg](img-9.jpeg)

**Fig. 10.** Snapshots at (from left to right) t=0.1, 0.3, 0.5 and 0.7 ms of simulations of 3-D drop impact on liquid film. The right-most column shows snapshots at t=0.7 ms from experiments. The rows (from top to bottom) correspond to We=82, 167 and 328. The planes behind iso-surfaces representing the liquid show the distribution of density and *t/ρr* on the central plane.

*Source:* Figure is reproduced from [51].

![img-10.jpeg](img-10.jpeg)

**Fig. 11.** Stability domain as obtained for the 2-D drop impact on liquid film simulation as a function of non-dimensional impact velocity *U<sub>0</sub>* and viscosity. The drop radius is set to *R<sub>0</sub>* = 50*λv*. Single relaxation time collision operator stability domain are shown with red circular symbols while those of the entropic model are shown with blue square markers.

*Source:* Figure is reproduced from [51].

![img-11.jpeg](img-11.jpeg)

**Fig. 12.** (right) Time-evolution of maximum spurious currents along with the density and velocity fields at the converged state as obtained from simulations with SRT and E-MRT models at ν = 0.03 and *T<sub>c</sub>* = 0.59.

*Source:* Figure is reproduced from [51].

![img-12.jpeg](img-12.jpeg)

**Fig. 13.** Maximum spurious currents for different viscosities at *T<sub>c</sub>* = 0.59 as obtained from 2-D drop simulations with (red circular markers) single relaxation time and (square blue markers) entropic collision operators.

*Source:* Figure is reproduced from [51].

$$
\rho \mathbf{u} = \sum_{i=1}^{Q} \mathbf{c}_i f_i \tag{105b}
$$

$$
2\rho E^{\text{tot}} = 2C_v \rho T + \rho \mathbf{u}^2 = \sum_{i=1}^{Q} \mathbf{c}_i^2 f_i + g_i \tag{105c}
$$

where *E<sup>tot</sup>* is the total energy for a polyatomic gas. The second terms on the left hand-side of Eqs. (103a) and (103b) are introduced to lift the unity Prandtl limitation of the BGK collision operator. The quasi-equilibrium, *f<sub>i</sub><sup>+</sup>* is computed, via a Grad's approximation, as:

$$
f_i^+ = f_i^{\text{eq}} + w_i \frac{(\widehat{\Pi}_{\alpha \beta \gamma} - \Pi_{\alpha \beta \gamma}^{\text{eq}})(\mathbf{c}_{i\alpha} \mathbf{c}_{i\beta} \mathbf{c}_{i\gamma} - 3\mathbf{c}_{i\gamma} T \delta_{\alpha \beta})}{6T^3}, \text{ if Pr} \leq 1
$$

$$
f_i^+ = f_i^{\text{eq}} + w_i \frac{(\widehat{\Pi}_{\alpha \beta}^+ - \Pi_{\alpha \beta}^{\text{eq}})(\mathbf{c}_{i\alpha} \mathbf{c}_{i\beta} - T \delta_{\alpha \beta})}{2T^2}, \text{ if Pr} > 1
$$

where:

$$
\Pi_{\alpha \beta}^{+} = \sum_{i=1}^{Q} \left[ (\mathbf{c}_{i\alpha} - u_{\alpha}) (\mathbf{c}_{i\beta} - u_{\beta}) - \frac{2}{D} \partial_{\gamma} u_{\gamma} \delta_{\alpha \beta} \right] f_i,
$$

while the quasi-equilibrium *g<sub>i</sub><sup>+</sup>* is:

$$
g_i^+ = g_i^{\text{eq}} + w_i \frac{(q_{\alpha} - q_{\alpha}^{\text{eq}}) c_{i\alpha}}{T}, \text{ if Pr} \leq 1
$$

$$
g_i^+ = g_i^{\text{eq}} + w_i (K - \rho (2C_v - D) T), \text{ if Pr} > 1
$$

with:

$$
q_{\alpha} = \sum_{i=1}^{Q} (\mathbf{c}_{i\alpha} - u_{\alpha}) g_i,
$$

and

$$
K = \sum_{i=1}^{Q} g_i.
$$

Upon proper construction of the discrete equilibrium *f<sub>i</sub><sup>+</sup>* , the system of equations presented here recovers, in the hydrodynamic limit, the Navier–Stokes–Fourier system of equations [54,55]:

$$
\partial_{\nu} \rho + \nabla \cdot \rho \mathbf{u} = 0, \tag{111a}
$$

$$
\partial_{\nu} \rho \mathbf{u} + \nabla \cdot \rho \mathbf{u} \otimes \mathbf{u} + \nabla \cdot \mathbf{T}_{\text{NS}} = 0, \tag{111b}
$$

$$
\partial_{\nu} \rho C_v T + \nabla \cdot \rho C_{\rho} T - \mathbf{u} \cdot \nabla + \nabla \mathbf{u} : \mathbf{T}_{\text{NS}} - \nabla \cdot \lambda \nabla T = 0. \tag{111c}
$$

We will discuss the construction of the discrete equilibrium for higher order lattices intended for compressible flows in the next sections.

![img-9.jpeg](images/page10-img-9.jpeg.jpeg)



![img-10.jpeg](images/page10-img-10.jpeg.jpeg)



![img-11.jpeg](images/page10-img-11.jpeg.jpeg)



![img-12.jpeg](images/page10-img-12.jpeg.jpeg)



---

<!-- Page 11 -->

### 4.1. Entropic equilibrium construction

The entropic construction of the discrete equilibrium state introduced for isothermal models, see Eq. (22), can be reformulated in a more general form as a minimization problem subject to *M* constraints:

$$
\delta H + \delta t \sum_{i=1}^{M} \lambda_m \Pi_m = 0. \tag{112}
$$

The formal solution of this constrained minimization leads to a function of the following form:

$$
f_i^{\text{eq}} = \rho w_i \exp \left[ \sum_{m=1}^{M} \lambda_m \left( \sum_{j=1}^{Q} \frac{\partial \Pi_m}{\partial f_j} \right) \right]. \tag{113}
$$

Note that other forms of the minimizer without the weights *w<sub>i</sub>* have also been proposed and used in the literature [56], most notably for entropic Grad moments methods [56,57]. For instance, a model imposing only constraints on collisional invariants, i.e.

$$
\sum_{i=1}^{Q} f_i^{\text{eq}} = \rho, \tag{114a}
$$

$$
\sum_{i=1}^{Q} e_i^j f_i^{\text{eq}} = \rho u, \tag{114b}
$$

$$
\sum_{i=1}^{Q} e_i^2 f_i^{\text{eq}} = \rho (u^2 + D r T), \tag{114c}
$$

would lead to the following discrete equilibrium [58]:

$$
f_i^{\text{eq}} = \rho w_i \exp \left[ \lambda_0 + \lambda_0 c_{i0} + \lambda_2 e_i^2 \right]. \tag{115}
$$

It is interesting to note that while, for the most part, entropic equilibria construction has been done by enforcing constraints on collisional invariants, one may reduce higher-order moments error by adding corresponding constraint in Eq. (112). This is sometimes referred to as *guiding* the equilibrium and corresponding discrete equilibria are referred to as *guided equilibria* [59,60]. In the context of the lattice Boltzmann method, this extension of constraints was discussed for the first time in [13] through the concept of auxiliary and target equilibria. There, auxiliary equilibria were constructed by enforcing constraints on collisional invariants and target equilibria, a combination of auxiliary equilibria and additional degrees of freedom, by enforcing constraints on higher order moments.

Once the form of the equilibrium distribution function has been determined, its construction consists of finding the expression of the different Lagrangian multiplicators. This is done by introducing back the discrete equilibrium into the set of constraints which would lead to a system of *M* equations with *M* unknowns, i.e. the Lagrange multipliers, to be determined. While an analytical expression was derived for the isothermal case with *D + 1* constraints, for larger systems no such solutions exist. In the absence of a closed form solution one can use numerical methods such as Newton iterations to find the Lagrange multipliers at every grid-point and every time-step [54].

### 4.2. Higher-order lattice construction

As shown in previous sections, one systematic approach to choose an optimal set of discrete velocities is to rely on the Gauss–Hermite quadrature and roots of Hermite polynomials. However, apart from the third-order quadrature leading to the *D∂Q3<sup>δ</sup>* lattices, all other higher order quadratures result in off-lattice propagation of some of the discrete distribution functions. In [61], starting from a set of discrete velocities the authors proposed an approach to find a reference temperature and corresponding weights. This is achieved through the *closure relation* and *matching* conditions. For a set of discrete velocities

![img-13.jpeg](img-13.jpeg)

**Fig. 14.** Illustration of the *D2Q49* lattice.

*V* with *Q* vectors *c<sub>i</sub>*, the *Q*<sup>th</sup> power of *c<sub>i</sub>* can be written as a linear combination of lower order odd-powers from *Q − 2* to 1, i.e.

$$
c_i^Q = a_{Q-2} e_i^{Q-2} + a_{Q-4} e_i^{Q-4} + \dots + a_1 c_i. \tag{116}
$$

For instance, in the case of the *D1Q3* lattice one has *c<sub>i</sub><sup>3</sup>* = *c<sub>i</sub>*. This essentially means that the moment of order *Q* is not an independent moment and cannot be set at one's will. The only possibility is to set the linear in *u* term of the *Q*<sup>th</sup> order to its Maxwell–Boltzmann counter-part and in doing so determine the reference temperature, which is referred to as the *matching* condition. Consider for instance the *D1Q3* lattice again. The third order moment is going to be *u<sub>x</sub>* while the Maxwell–Boltzmann distribution leads to *u<sub>x</sub><sup>3</sup>* + 3*T<sub>0</sub>u<sub>x</sub>*. To match the linear term one must have 3*T<sub>0</sub> = 1. Note that not any choice of lattice admits a reference temperature. For example the velocity set *V* = {−2,−1,0,+1,+2} will lead to a closure relation of the form *c<sub>i</sub><sup>3</sup>* = 5*c<sub>i</sub><sup>3</sup>* − 4*c<sub>i</sub>* and a matching condition, 15*T<sub>0</sub><sup>2</sup>* − 15*T<sub>0</sub> + 4 = 0, which does not admit any solutions. This explains why the shortest admissible five-velocity lattice is *V* = {−3,−1,0,+1,+3} with *T<sub>0</sub>* = 1 ± √2/5. Once the reference temperature is determined, the weights are readily found by matching the moments of the discrete equilibrium at *ρ* = 1 and *u<sub>x</sub>* = 0 to their Maxwell–Boltzmann counter-parts. Considering the condition of positivity of the weights one also finds the range of temperature that can be covered by the chosen system of discrete velocities. The closure relations and reference temperatures of a number of 1-D lattice are summarized in Table 1.

One successful example of such lattices is the *D2Q49* shown in Fig. 14. The closure relation for the 1-D set is

$$
c_i^2 = 14c_i^5 - 49c_i^3 + 36c_i. \tag{117}
$$

The 1-D weights read:

$$
w_0 = \frac{36 - 49T + 42T^2 - 15T^3}{36}, \tag{118a}
$$

$$
w_{\pm 1} = \frac{T(12 - 13T + 5T^2)}{16}, \tag{118b}
$$

$$
w_{\pm 2} = \frac{T(-3 + 10T - 5T^2)}{40}, \tag{118c}
$$

$$
w_{\pm 3} = \frac{T(4 - 15T + 15T^2)}{720}, \tag{118d}
$$

which lead to *T<sub>min</sub>* = 1 − √2/5 and *T<sub>max</sub>* = 1 + √2/5. Note that the range of accessible temperatures can be further extended by changing the ratio of the largest and shortest discrete velocities, here ±3 and ±1. In [54] the author also proposed pruning strategies to reduce the number of discrete velocities in 2- and 3-D, leading to the *D3Q39* lattice, which reduces the discrete velocities by one order of magnitude compared to the tensor product of the *D1Q7*, i.e. *D3Q343*.

![img-13.jpeg](images/page11-img-13.jpeg.jpeg)



---

<!-- Page 12 -->

Table 1 One-dimensional Maxwell lattices with odd number of integer-valued velocities, Q = 3, 5, 7, 9, 11. Second column: Lattice vectors; Third column: Closure relation, defining the reference temperature T_{0} through the matching condition (fourth column).

|  Q | ν | Closure | T_{0}  |
| --- | --- | --- | --- |
|  3 | {0, ±1} | c_{1}^{3} = c_{1} | 1/3  |
|  5 | {0, ±1, ±3} | c_{1}^{4} = 10c_{1}^{4} - 9c_{1} | 1 ± √2/5  |
|  7 | {0, ±1, ±2, ±3} | c_{1}^{5} = 14c_{1}^{5} - 49c_{1}^{5} + 36c_{1} | 0.697953  |
|  9 | {0, ±1, ±2, ±3, ±5} | c_{1}^{6} = 39c_{1}^{6} - 399c_{1}^{6} + 1261c_{1}^{6} - 900c_{1} | 0.756081, 2.175382  |
|  11 | {0, ±1, ±2, ±3, ±4, ±5} | c_{1}^{11} = 55c_{1}^{6} - 1023c_{1}^{6} + 7645c_{1}^{6} - 21076c_{1}^{6} + 14400c_{1} | 1.062794  |

![img-14.jpeg](img-14.jpeg)

Fig. 15. Illustration of the D2Q49 lattice with a shift of U_{x} = δr/δt.

### 4.3. Strategies to extend operation range: shifted lattices

As observed for both isothermal and compressible models, errors in higher-order moments scale with the deviations of local temperature and velocity from the lattice reference temperature and velocity. For all symmetric lattices considered up to that point the lattice reference velocity is U = 0. In [62] the authors proposed to challenge the idea of a reference frame at rest by introducing a non-zero shift U. It was noted that the discrete entropy functional as defined in Eq. (25) is uniquely defined by the weights w_{i}. The weights of a lattice with Q discrete velocities, as shown in the previous section, are determined by matching the first Q moments of the Maxwell–Boltzmann equilibrium distribution function at temperature T and u_{x} = 0:

$$\sum_{i=1}^{Q} \phi(c_i) w_i(0, T) = \int \phi(v)f^{\text{MB}}(0, T) dv. \tag{119}$$

It was shown through the Galilean-invariance of the moments of the Maxwell–Boltzmann distribution function and the binomial theorem that the weights are also Galilean-invariant and therefore untouched by the change of reference frame. The immediate consequences of that observation are: (a) construction of 3-D lattices via tensorial product of the 1-D lattice remains as before, (b) assuming U = kδr/δt with k ∈ Z the propagation remains on-lattice and (c) the discrete entropy functional is Galilean invariant and therefore equilibrium populations are form invariant under the shift of reference frame. This point along with the effect of the shift on operation range has also been discussed for standard isothermal lattices [26]. The process of changing the reference frame and resulting discrete lattice is illustrated in Fig. 15 through the D2Q49 lattice. The use of the shifted lattice along with the entropic equilibrium has been successfully used to model a wide variety of high Mach number flows as illustrated in Fig. 16.

### 5. Conclusion

While the majority of efforts in improving the LBGK in the early 2000s up to the 2010s was focused multiple relaxation time formulations, the concept of a discrete H theorem and Lyapunov functionals to restore stability in the LBGK provided a different perspective. The concept of entropic LBM's affected the original BGK at two different levels: (a) construction of the discrete equilibrium state and (b) relaxation process. At the discrete equilibrium level, the attractor is constructed by considering conditional minimization of a convex Lyapunov functional leading to a considerably wider domain of linear stability. At the level of the relaxation operator, the introduction of the dynamic parameter a prevents the discrete distribution functions from going beyond their mirror state and in doing so guarantees non-linear stability. A testament to the advantages of the entropic construction is that over the past two decades it has been successfully extended to different areas of application such as multiphase and compressible flows, also showing the robustness of the concept.

![img-15.jpeg](img-15.jpeg)

Fig. 16. Drag coefficient c_{d} as a function of the free stream Mach number for the Busemann biplane simulations. Inset: snapshots of the pressure distribution around the biplane for three different Mach numbers: Ma = 1.5, top; Ma = 1.7, bottom left; Ma = 2.0, bottom right. Source: Figure is reproduced from [58].

### Declaration of competing interest

The authors declare that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

### Data availability

No data was used for the research described in the article.

### Acknowledgments

S. A. H. and I. V. K. gratefully acknowledge support by European Research Council (ERC) Advanced Grant 834763-PonD. Computational resources at the Swiss National Super Computing Center CSCS were provided under the grants s1066 and s1222.

### References

- [1] McNamara GR, Zanetti G. Use of the boltzmann equation to simulate lattice-gas automata. Phys Rev Lett 1988;61(20):2332, Publisher: APS.
- [2] Frisch U, Hasslacher B, Pomeau Y. Lattice-gas automata for the navier–stokes equation. Phys Rev Lett 1986;56(14):1505–8. http://dx.doi.org/10.1103/PhysRevLett.56.1505, URL https://link.aps.org/doi/10.1103/PhysRevLett.56.1505.

![img-14.jpeg](images/page12-img-14.jpeg.jpeg)



![img-15.jpeg](images/page12-img-15.jpeg.jpeg)



---

<!-- Page 13 -->

[3] Succi S. The lattice Boltzmann equation for fluid dynamics and beyond. Clarendon; 2002, URL http://physicstoday.sci tation.org/doi/10.1063/1.1537916.
[4] Chen S, Chen H, Martnez D, Matthaeus W. Lattice boltzmann model for simulation of magnetohydrodynamics. Phys Rev Lett 1991;67(27):3776-9. http://dx.doi.org/10.1103/PhysRevLett.67.3776, URL https://link.aps.org/doi/ 10.1103/PhysRevLett.67.3776.
[5] Bhatnagar PL, Gross EP, Krook M. A model for collision processes in gases. I. small amplitude processes in charged and neutral one-component systems. Phys Rev 1954;94(3):511-25. http://dx.doi.org/10.1103/PhysRev.94.511, URL https://link.aps.org/doi/10.1103/PhysRev.94.511.
[6] Chorin AJ. A numerical method for solving incompressible viscous flow problems. J Comput Phys 1997;135(2):118-25. http://dx.doi.org/10.1006/jcph.1997.5716, URL https://linkinghub.elsevier.com/retrieve/pii/S0021999197957168.
[7] Lallemand P, Luo L-S. Theory of the lattice Boltzmann method: Dispersion, dissipation, isotropy, galllean invariance, and stability. Phys Rev E 2000;61(6):6546-62. http://dx.doi.org/10.1103/PhysRevE.61.6546, URL https: //link.aps.org/doi/10.1103/PhysRevE.61.6546.
[8] Geier M, Greiner A, Korvink JG. Cascaded digital lattice Boltzmann automata for high revmolds number flow. Phys Rev E 2006;73(6):066705, Publisher: APS.
[9] Geier M, Schönherr M, Pasquali A, Krafczyk M. The cumulant lattice Boltzmann equation in three dimensions: Theory and validation. Comput Math Appl 2015;70(4):507-47. http://dx.doi.org/10.1016/j.camwa.2015.05.001, URL https: //linkinghub.elsevier.com/retrieve/pii/S0898122115002126.
[10] Latt J, Chopard B. Lattice Boltzmann method with regularized pre-collision distribution functions. Math Comput Simulation 2006;72(2-6):165-8, Publisher: Elsevier.
[11] Karlin I, Gorban A, Dukek G, Nonnenmacher T. Dynamic correction to moment approximations. Phys Rev E 1998;57(2):1668-72. http://dx.doi.org/10.1103/ PhysRevE.57.1668, URL https://link.aps.org/doi/10.1103/PhysRevE.57.1668.
[12] Tsallis C. Possible generalization of Boltzmann-Gibbs statistics. J Stat Phys 1988;52(1-2):479-87. http://dx.doi.org/10.1007/BF01016429, URL http://link. springer.com/10.1007/BF01016429.
[13] Karlin IV, Succi S. Equilibria for discrete kinetic equations. Phys Rev E 1998;58(4):R4053-6. http://dx.doi.org/10.1103/PhysRevE.58.R4053, URL https: //link.aps.org/doi/10.1103/PhysRevE.58.R4053.
[14] Karlin IV, Ferrante A, Öttinger HC. Perfect entropy functions of the lattice Boltzmann method. Europhys Lett 1999;47(2):182-8. http://dx.doi.org/10.1209/ epl/11999-00370-1, URL https://eopscience.iop.org/article/10.1209/epi/1199900370-1.
[15] Ansumali S, Karlin IV. Entropy function approach to the lattice Boltzmann method. J Stat Phys 2002;107(1/2):291-308. http://dx.doi.org/10.1023/A: 1014575024265, URL http://link.springer.com/10.1023/A:1014575024265.
[16] Ansumali S, Karlin IV, Öttinger HC. Minimal entropic kinetic models for hydrodynamics. Europhys Lett 2003;63(6):798, Publisher: IOP Publishing.
[17] Ansumali S, Karlin IV. Consistent lattice Boltzmann method. Phys Rev Lett 2005;95(26):260605. http://dx.doi.org/10.1103/PhysRevLett.95.260605, URL https://link.aps.org/doi/10.1103/PhysRevLett.95.260605.
[18] Chikatamarla SS, Ansumali S, Karlin IV. Entropic lattice Boltzmann models for hydrodynamics in three dimensions. Phys Rev Lett 2006;97(1):010201, Publisher: APS.
[19] Methods of numerical integration. Elsevier; 1984, http://dx.doi.org/10. 1016/C2013-0-10566-1, URL https://linkinghub.elsevier.com/retrieve/pii/ C20130105661.
[20] Hosseini SA, Karlin IV. Entropic equilibrium for the lattice Boltzmann method: hydrodynamics and numerical properties. 2023, http://dx.doi.org/10.48550/ ARXIV.2303.08163, URL https://arxiv.org/abs/2303.08163, Publisher: arXiv Version Number: 1.
[21] Chávez-Modena M, Ferrer E, Rubio G. Improving the stability of multiplerelaxation lattice Boltzmann methods with central moments. Comput \& Fluids 2018;172:397-409. http://dx.doi.org/10.1016/j.compfluid.2018.03.084, URL https://linkinghub.elsevier.com/retrieve/pii/S0045793018301889.
[22] Sterling JD, Chen S. Stability analysis of lattice Boltzmann methods. J Comput Phys 1996;123(1):196-206. http://dx.doi.org/10.1006/jcph.1996.0016, URL https://linkinghub.elsevier.com/retrieve/pii/S0021999196900169.
[23] Worthing RA, Mozer J, Seeley G. Stability of lattice Boltzmann methods in hydrodynamic regimes. Phys Rev E 1997;56(2):2243-53. http://dx.doi.org/10.1103/ PhysRevE.56.2243, URL https://link.aps.org/doi/10.1103/PhysRevE.56.2243.
[24] Hosseini SA, Darabiha N, Thévenin D, Eshghinejadfard A. Stability limits of the single relaxation-time advection-diffusion lattice Boltzmann scheme. Internat J Modern Phys C 2017;28(12):1750141. http://dx.doi.org/10.1142/ S0129183117501418, URL https://www.worldscientific.com/doi/abs/10.1142/ S0129183117501418.
[25] Hosseini SA, Coreixas C, Darabiha N, Thévenin D. Stability of the lattice kinetic scheme and choice of the free relaxation parameter. Phys Rev E 2019;99(6):063305. http://dx.doi.org/10.1103/PhysRevE.99.063305, URL https: //link.aps.org/doi/10.1103/PhysRevE.99.063305.
[26] Hosseini SA, Coreixas C, Darabiha N, Thévenin D. Extensive analysis of the lattice Boltzmann method on shifted stencils. Phys Rev E 2019;100(6):063301. http://dx.doi.org/10.1103/PhysRevE.100.063301, URL https://link.aps.org/doi/ 10.1103/PhysRevE.100.063301.
[27] Wissocq G, Sagaut P, Boussuge J-F. An extended spectral analysis of the lattice Boltzmann method: modal interactions and stability issues. J Comput Phys 2019;380:311-33. http://dx.doi.org/10.1016/j.jcp.2018.12.015, URL https: //linkinghub.elsevier.com/retrieve/pii/S0021999118308118.
[28] Hosseini SA. Development of a lattice Boltzmann-based numerical method for the simulation od reacting flows [Ph.D. thesis], Universite Paris-Saclay/Otto-von-Guericke-Universität Magdeburg; 2020, URL https://opendata.uni-halle.de/ /handle/1981185920/33674.
[29] Karlin IV, Ansumali S, DE Angelis E, Öttinger HC, Succi S. Entropic lattice Boltzmann method for large scale turbulence simulation. 2003, http: //dx.doi.org/10.48550/ARXIV.CO9D-MAT-9306003, URL https://arxiv.org/abs/ cond-mat/0306003, Publisher: arXiv Version Number: 1.
[30] Malaspinas O, Deville M, Chopard B. Towards a physical interpretation of the entropic lattice Boltzmann method. Phys Rev E 2008;78(6):066705. http://dx.doi.org/10.1103/PhysRevE.78.066705, URL https://link.aps.org/doi/ 10.1103/PhysRevE.78.066705.
[31] Malaspinas O. Increasing stability and accuracy of the lattice Boltzmann scheme: recurricity and regularization. 2015, http://dx.doi.org/10.48550/ARXIV.1505. 06900, URL https://arxiv.org/abs/1505.06900, Publisher: arXiv.
[32] Atif M, Kolluru PK, Ansumali S. Essentially entropic lattice Boltzmann model: Theory and simulations. Phys Rev E 2022;106(5):055307. http://dx.doi.org/10. 1103/PhysRevE.106.055307, URL https://link.aps.org/doi/10.1103/PhysRevE. 106.055307.
[33] Smagorinsky J, Manabe S, Holloway JL. Numerical results from a nine-level general circulation model of the atmosphere. Mon Weather Rev 93(12):727-68.
[34] Buzzicotti M, Tauzin G. Inertial range statistics of the entropic lattice Boltzmann method in three-dimensional turbulence. Phys Rev E 2021;104(1):015302. http://dx.doi.org/10.1103/PhysRevE.104.015302, URL https://link.aps.org/doi/ 10.1103/PhysRevE.104.015302.
[35] Atif M, Kolluru PK, Thantanapally C, Ansumali S. Essentially entropic lattice Boltzmann model. Phys Rev Lett 2017;119(24):240602. http://dx.doi. org/10.1103/PhysRevLett.119.240602, URL https://link.aps.org/doi/10.1103/ PhysRevLett.119.240602.
[36] Karlin IV, Bösch F, Chikatamarla S. Gibbs' principle for the lattice-kinetic theory of fluid dynamics. Phys Rev E 2014;90(3):031302, Publisher: APS.
[37] Karlin I, Bösch F, Chikatamarla S, Succi S. Entropy-assisted computing of lowdissipative systems. Entropy 2015;17(12):8099-110. http://dx.doi.org/10.3390/ e17127867, URL http://www.mdpi.com/1099-4300/17/12/7867.
[38] Wang L. Enhanced multi-relaxation-time lattice Boltzmann model by entropic stabilizers. Phys Rev E 2020;102(2):023307. http://dx.doi.org/10.1103/PhysRevE. 102.023307, URL https://link.aps.org/doi/10.1103/PhysRevE.102.023307.
[39] Dorschner R, Chikatamarla SS, Karlin IV. Transitional flows with the entropic lattice Boltzmann method. J Fluid Mech 2017;824:388-412. http://dx. doi.org/10.1017/jfm.2017.356, URL https://www.cambridge.org/core/product/ identifier/S0022112017003561/type/journal_article.
[40] Dorschner R, Bösch F, Chikatamarla SS, Bouboucha K, Karlin IV. Entropic multi-relaxation time lattice Boltzmann model for complex flows. J Fluid Mech 2016;801:623-51. http://dx.doi.org/10.1017/jfm.2016.448, URL https://www.cambridge.org/core/product/identifier/S0022112016004481/ type/journal_article.
[41] Anderson DM, McFadden GB, Wheeler AA. Diffuse-interface methods in fluid mechanics. Annu Rev Fluid Mech 1998;30(1):139-65, Publisher: Annual Reviews 4139 El Camino Way, PO Box 10139, Palo Alto, CA 94303-0139, USA.
[42] Korteweg DJ. Sur la forme que prennent les équations du mouvement des fluides si l'on tient compte des forces capillaires causées par des variations de densité considérables mais continues et sur la théorie de la capillarité dans l'hypothèse d'une variation continue de la densité. Arch Névél Sci Exactes Nat 1901;6(265).
[43] Hosseini SA, Karlin IV. Lattice Boltzmann for non-ideal fluids: fundamentals and practice. 2023, http://dx.doi.org/10.48550/ARXIV.2301.02011, URL https: //arxiv.org/abs/2301.02011, Publisher: arXiv Version Number: 1.
[44] Moquélien AM. Entropic lattice Boltzmann method for two-phase flows [Ph.D. thesis], ETH Zurich; 2016, http://dx.doi.org/10.3929/ETHZ-A-010809769, URL http://hdl.handle.net/20.500.11850/126173, Artwork Size: 189 p. Medium: application/pdf/4spec: 189 p.
[45] Swift MR, Osborn WA, Yeomans JM. Lattice Boltzmann simulation of nonideal fluids. Phys Rev Lett 1995;75(5):830-3. http://dx.doi.org/10.1103/PhysRevLett. 75.830, URL https://link.aps.org/doi/10.1103/PhysRevLett.75.830.
[46] Hosseini S, Dorschner R, Karlin I. Towards a consistent lattice Boltzmann model for two-phase fluids. J Fluid Mech 2022;953:A4. http://dx.doi.org/10. 1017/jfm.2022.867, URL https://www.cambridge.org/core/product/identifier/ S0022112022008679/type/journal_article.
[47] Mazloomi A, Chikatamarla S, Karlin I. Entropic lattice Boltzmann method for multiphase flows. Phys Rev Lett 2015;114(7):174502.

---

<!-- Page 14 -->

[48] Mazloomi Moqaddam A, Chikatamarla SS, Karlin IV. Simulation of binary droplet collisions with the entropic lattice Boltzmann method. Phys Fluids 2016;28(2):022106. http://dx.doi.org/10.1063/1.4942017, URL http://aip. scitation.org/doi/10.1063/1.4942017.
[49] Bösch F, Dorschner B, Karlin I. Entropic multi-relaxation free-energy lattice Boltzmann model for two-phase flows. Europhys Lett 2018;122(1):14002, Publisher: IOP Publishing.
[50] Qin F, Mazloomi Moqaddam A, Kang Q, Derome D, Carmeliet J. Entropic multiple-relaxation-time multirange pseudopotential lattice Boltzmann model for two-phase flow. Phys Fluids 2018;30(3):032104. http://dx.doi.org/10.1063/1. 5016965, URL http://aip.scitation.org/doi/10.1063/1.5016965.
[51] Hosseini SA, Dorschner B, Karlin IV. Entropic multi-relaxation-time lattice Boltzmann model for large density ratio two-phase flows. Commun Comput Phys 2022. URL https://arxiv.org/abs/2201, Publisher: arXiv Version Number: 2.
[52] Luo KH, Fei L, Wang G. A unified lattice Boltzmann model and application to multiphase flows. Phil Trans R Soc A 2021;379(2208):20200397. http://dx. doi.org/10.1098/rsta.2020.0397, URL https://royalsocietypublishing.org/doi/10. 1098/rsta.2020.0397.
[53] Wang G, Fei L, Lei T, Wang Q, Luo KH. Droplet impact on a heated porous plate above the leidenfrost temperature: A lattice Boltzmann study. Phys Fluids 2022;34(9):093319. http://dx.doi.org/10.1063/5.0118079, URL https://aip. scitation.org/doi/10.1063/5.0118079.
[54] Frapolli N. Entropic lattice Boltzmann models for thermal and compressible flows [Ph.D. thesis], ETH Zurich; 2017, http://dx.doi.org/10.3929/ETHZ-A010890892, URL http://hdl.handle.net/20.500.11850/130664, Artwork Size: 273 p. Medium: application/pdf Pages: 273 p.
[55] Frapolli N, Chikatamarla SS, Karlin IV. Entropic lattice Boltzmann model for compressible flows. Phys Rev E 2015;92(6):061301. http://dx.doi.org/10. 1103/PhysRevE.92.061301, URL https://link.aps.org/doi/10.1103/PhysRevE.92. 061301 .
[56] Öttinger HC, Struchtrup H, Torrilhon M. Formulation of moment equations for rarefied gases within two frameworks of non-equilibrium thermodynamics: RET and GENERIC. Phil Trans R Soc A 2020;378(2170):20190174. http://dx. doi.org/10.1098/rsta.2019.0174, URL https://royalsocietypublishing.org/doi/10. 1098/rsta.2019.0174.
[57] Levermore CD. Moment closure hierarchies for kinetic theories. J Stat Phys 1996;83(5-6):1021-65. http://dx.doi.org/10.1007/BF02179552, URL http ://link.springer.com/10.1007/BF02179552.
[58] Frapolli N, Chikatamarla S, Karlin I. Theory, analysis, and applications of the entropic lattice Boltzmann model for compressible flows. Entropy 2020;22(3):370. http://dx.doi.org/10.3390/e22030370, URL https://www.mdpi.com/1099-4300/ 22/3/370.
[59] Prasianakis NI, Karlin IV. Lattice Boltzmann method for thermal flow simulation on standard lattices. Phys Rev E 2007;76(1):016702. http://dx.doi.org/10. 1103/PhysRevE.76.016702, URL https://link.aps.org/doi/10.1103/PhysRevE.76. 016702 .
[60] Latt J, Coreixas C, Beny J, Parmigiani A. Efficient supersonic flow simulations using lattice Boltzmann methods based on numerical equilibria. Phil Trans R Soc A 2020;378(2175):20190559. http://dx.doi.org/10.1098/rsta.2019.0559, URL https://royalsocietypublishing.org/doi/10.1098/rsta.2019.0559.
[61] Karlin I, Asinari P. Factorization symmetry in the lattice Boltzmann method. Physica A 2010;389(8):1530-48. http://dx.doi.org/10.1016/j.physa.2009.12. 032, URL https://linkinghub.elsevier.com/retrieve/pii/S0378437109010462.
[62] Frapolli N, Chikatamarla S, Karlin I. Lattice kinetic theory in a comoving galilean reference frame. Phys Rev Lett 2016;117(1):010604. http:// dx.doi.org/10.1103/PhysRevLett.117.010604, URL https://link.aps.org/doi/10. 1103/PhysRevLett.117.010604.