---

<!-- Page 0 -->

# Entropic Lattice Boltzmann Method 

Keming Yu

Approved by Sorin Mitran (Thesis Advisor), Boyce Eugene Griffith (Defense Committee Member), Shahar Kovalsky (Defense Committee Member)

A thesis presented for the degree of Bachelor of Science in Mathematics
![img-0.jpeg](img-0.jpeg)

Department of Mathematics
University of North Carolina in Chapel Hill
Apr. 2021

![img-0.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page0-img-0.jpeg.jpeg)



---

<!-- Page 1 -->

# Entropic Lattice Boltzmann Method 

Keming Yu


#### Abstract

Lattice Boltzmann methods are kinetic descriptions of fluid flow that are efficiently implemented through a stream and collide approach. The collision operation is typically an approximation of the microscopic physics, with the BGK linear approximation a widely used choice. It has a wide range of application in computational fluid dynamics.

The honor thesis starts with an introduction to the kinetic theory of gas. First, we introduce the Liouville Equation that describes the evolution of particle distribution function of a system in phase space. With an analysis of the Liouville Equation, BBGKY hierarchy is introduced for the reduced distribution evolution. Using the Bogoliubov Hypothesis, we are able to close the hierarchy equation and derive the Boltzmann's Equation for kinetic theory of gas. An analysis of Boltzmann's Equation is presented, with the center of analysis being the Chapman-Enskog analysis that reproduces the Navier-Stokes Equation in its expansion.

Then, the Lattice Boltzmann Method is introduced by describing collision and streaming steps. The formulation of Lattice Boltzmann Method is centered around discretizing the velocity space by using Gauss-Hermite Quadrature that gives a truncated version of the Maxwell-Boltzmann distribution. Under this truncation, Lattice Boltzmann Method becomes essentially a forward Euler scheme with BGK operator under a simplification of collision physics that is described by a limited direction of possible microscopic velocities.

The BGK relaxation approximation however suffers from known deficiencies leading to instability at high Reynolds numbers that are the result of not satisfying a microscopic H-theorem. Entropic lattice methods seek an alternative relaxation procedure that is guaranteed to be unconditionally stable. They however typically are computationally expensive requiring solution of a nonlinear equation at each lattice site and each time step. This honor thesis further introduces the possibility of an alternative relaxation approach based upon geodesic transport on a statistical manifold within the overall framework of information geometry.

---

<!-- Page 2 -->

# Chapter 1 

## Kinetic Theory of Gas

For the discussion of kinetic theory of gas, we follow a standard textbook: [10] Richard L. Liboff. Kinetic theory: classical, quantum, and relativistic descriptions.

### 1.1 The Livouille Equation

### 1.1.1 Classical Mechanics Review

To fully determine a state of a given system, a set of independent generalized coordiantes is employed, where the degree of freedom is the minimum number to determine this given system. For a system of $N$ degree of freedom, the set of generalized coordiates is $\left\{q_{k} \mid k=1 \ldots . N\right\}$ or can be denoted as a vector. Correspondingly, the generalzied velocities are calculated with a differentiation in time.

$$
\begin{aligned}
& \boldsymbol{q}=\left(q_{1}, q_{2}, \ldots q_{N}\right) \\
& \dot{\boldsymbol{q}}=\left(q_{1}^{\prime}, q_{2}, \ldots q_{N}^{\prime}\right)
\end{aligned}
$$

The system with a conservative force field and kinetic energy is characterized by its Lagrangian.

$$
L(\boldsymbol{q}, \dot{\boldsymbol{q}})=T(\boldsymbol{q}, \dot{\boldsymbol{q}})-V(\boldsymbol{q})
$$

The Hamiltonian Principle dictates the dynamics of motion by minimizing or maximizing the action of the system.

$$
\delta \int_{1}^{2} L(\boldsymbol{q}, \dot{\boldsymbol{q}}, t) \mathrm{dt}=0
$$

From Hamiltonian Principle, a differential equation can be obtained for the equation of motion. This is the Lagrange's Equation

$$
\frac{d}{d t}\left(\frac{\partial L}{\partial \dot{q}_{l}}\right)-\frac{\partial L}{\partial q_{l}}=0 \quad(l=1 \ldots . N)
$$

After a Legendre transform of $(\boldsymbol{q}, \dot{\boldsymbol{q}}) \rightarrow(\boldsymbol{q}, \boldsymbol{p})$, where $p_{l}=\frac{\partial L}{\partial q_{l}}$ is the canonical momentum, the transform takes Lagrangian to a Hamiltonian.

$$
H=\sum_{l} \frac{\partial L}{\partial \dot{q}_{l}} \dot{q}_{l}-L
$$

Then the coordinates, velocities and momentum could have the following expression by differentiating the Hamiltonian

$$
\dot{q}_{l}=\frac{\partial H}{\partial p_{l}}, p_{l}=-\frac{\partial H}{\partial q_{l}}, \frac{\partial H}{\partial t}=-\frac{\partial L}{\partial t}
$$

This is to be further observed that if generalized coordinates are constant in time, it follows that Lagrangian is not dependent on time explicitly. For such system, the Hamitonian is constant in time and can be identified as energy.

$$
\frac{d}{d t}\left(\sum_{l} \frac{\partial L}{\partial \dot{q}_{l}} \dot{q}_{l}-L\right)=\frac{d}{d t} H=0 \quad H=E=\text { constant }
$$

---

<!-- Page 3 -->

System constant of motion and of degree of freedom $N$ lives in a $2 N$ Cartesian space with $q_{l}$ and $p_{l}$ as variables. This is the $\Gamma$-space. Dynamical variables live in $2 N+1$ Cartesian space, the $\bar{\Gamma}$-space. Let $u(\boldsymbol{q}, \boldsymbol{p}, t)$ denotes a dynamical variable that lives in $\bar{\Gamma}$-space. Using Hamilton's Equations 1.7, we can obtain the equation of motion for $u(\boldsymbol{q}, \boldsymbol{p}, t)$.

$$
\frac{d u}{d t}=\sum_{l}\left(\frac{\partial u}{\partial q_{l}} \frac{\partial H}{\partial p_{l}}-\frac{\partial u}{\partial p_{l}} \frac{\partial H}{\partial q_{l}}\right)+\frac{\partial u}{\partial t}=[u, H]+\frac{\partial u}{\partial t}
$$

Here, $[u, H]$ is a Poisson Bracket. If dynamical variable $u(\boldsymbol{q}, \boldsymbol{p}, t)$ is only a function of Hamiltonian $H$, the equation above goes to zero and $u(\boldsymbol{q}, \boldsymbol{p}, t)$ is constant of motion.

# 1.1.2 Canonical Tranformation 

A transformation $(q, p) \rightarrow\left(q^{\prime}, p^{\prime}\right)$ is canonical if before and after the transformation, the Hamiltonian Principle is satisifed. Therefore, for a canoncial transformation, the following relationship beteween the Lagrangian before and after the transformation is true.

$$
L\left(q^{\prime}, p^{\prime}\right)=L(q, p)+\frac{d G_{1}(q, p, t)}{d t}
$$

A function $G_{1}\left(q, q^{\prime}, t\right)$ must exist if the transformation is canonical. We can also express 1.10 in terms of the Hamiltonian of the system to observe the following quality of function $G_{1}$.

$$
p=\frac{\partial G_{1}}{\partial q}, p^{\prime}=-\frac{\partial G_{1}}{\partial q^{\prime}}, H-H^{\prime}=\frac{\partial G}{\partial t}
$$

An example of relation 1.11 is the exchange transformation, given by generating function below and the following antisymmetric transformation.

$$
\begin{array}{r}
G_{1}\left(q, q^{\prime}\right)=-q q^{\prime} \\
p^{\prime}=q, q^{\prime}=-p
\end{array}
$$

Another equivalent definition for canonical transformation is that new variables' Hamiltonian also satisfy the Hamiltonian equation 1.7. If a Legendre transformation is perfromed on generating function $G_{1}$ to new generating function $G_{2}=\sum p^{\prime} q^{\prime}+G_{1}$ have new relations.

$$
p=\frac{\partial G_{2}\left(q, p^{\prime}\right)}{\partial q}, q^{\prime}=\frac{\partial G_{2}\left(q, p^{\prime}\right)}{\partial p^{\prime}}
$$

Canonical invariants are quantities that remain invariant under canonical transformations. An example of such invariant is the above mentioned Poisson's Bracket. Canonical transformtion also satisfies the group property where the product of two canonical transformations is also canonical. If a variable commutes with Hamiltonian in terms of Poisson Bracket, then H must must be independent of that variable. It follows for such variable having no explicit dependence on time for its Hamiltonian, it must be constant in time.

### 1.1.3 Liouville Theorem

The Liouville Theorem states that the Jacobian of canonical transformation must be unity.

---

<!-- Page 4 -->

The Geometric significance of Liouville Theorem is volume elements in the $\Gamma$-space remain the same before and after the canonical transformation. Now, consider the action integral from time $t$ to time $t+T$ :

$$
S(t, T)=\int_{t}^{t+T} L(q, p) d t
$$

Converting to Hamiltonian and differentiating this integral, we obtain

$$
\frac{d S}{d t}=\sum p^{\prime} q^{\prime}-\sum p \dot{q}+\left(H-H^{\prime}\right)
$$

where we denote $q^{\prime}=q(t+T)$ and similar for momentum. If Hamiltonian is constant in time $H=H^{\prime}$, we know action would be a generating function because it is a differential. The conclusion of the above derivation is

1. Action $S$ is a generating function for physical motion in time.
2. Differential motion in time is a canonical transformation.
3. Following the group property of canonical transformation, extended motion in time must also be a canonical transformation.

# 1.1.4 Liouville's Equation 

The state of system is a point in the $\Gamma$-space, and as time evolves, it moves along certain trajectory in the $\Gamma$-space. Now consider large replica of the system. This abstract collection of system is called ensemble. If there are $\mathcal{N}$ systems in the ensemble, there will be $\mathcal{N}$ points in the $\Gamma$-space. Since a system with $N$ degree of freedom is specified with $2 N$ variables, the trajectories of system in the ensemble cannot cross each other in the $\Gamma$-space. We introduce a function $D(q, p, t)$ as

$$
D(q, p, t)=\frac{d \mathcal{N}}{d \Omega}
$$

where $d \Omega$ is a volume element about a point in the $\Gamma$-space. Therefore, function $D$ represents the density of systems. Consider an ensemble in a continously closed surface in the $\Gamma$-space. Since the trajectory never cross each other, the points in the surface will remain in the surface. Therefore, $d \mathcal{N}$ must be invariant under a canonical transformation. With Liouville's Theorem, we also know that $d \Omega$ is invariant under canonical transformation. Therefore, their ratio $\frac{d \mathcal{N}}{d \Omega}$ is also invariant under canonical transformation. Then, on a system trajectory, functon $D$ is constant. This is the Liouville's Equation.

$$
\frac{d D}{d t}=0
$$

Or, for any dynamical variable like $D$, we refer to equation 1.9

$$
\frac{\partial D}{\partial t}+[D, H]=0
$$

From the above, we know that a variable constant of motion must solve the Liouville's Equation. Therefore, knowing the solution to Liouville's Equation is equivalent of knowing all the orbits of systems in the ensemble. Another derivation of Liouville's Equation can come from fluid dynamics interpretation of ensemble in $\Gamma$-space. An equivalent way of writing equation 1.20:

$$
\frac{\partial D}{\partial t}+\boldsymbol{u} \cdot \nabla D=0
$$

where $\boldsymbol{u}$ is the velocity of system points. The above "continuity equation" can be solved in terms of method of characteristics, which further confirmed solution of Liouville's Equation gives the knowlege of all orbits of systems in the ensemble. The Liouville Equation can be solved in terms of Taylor expansion of $D$ around initial condition $D_{0}(q, p)$.

$$
D(q, p, \Delta t)=\left\{1+\Delta t[H, D]+\frac{(\Delta t)^{2}}{2}[H,[H, D]]+\cdots\right\} D_{0}(p, q)
$$

---

<!-- Page 5 -->

This solution can be seen as how $D$ evolves in time for a point in the hyper-plane $(q . p)$ in $\tilde{\Gamma}$-space. Liouville's Equation, as mentioned above, can be solved by method of characteristics. For orbits:

$$
q=q_{0}+\tilde{q}(t), p=p_{0}+\tilde{p}(t)
$$

Function $\tilde{q}(t)$ and $\tilde{p}(t)$ vanishes at $t=0$ because the initial conditions are known. Function $D$ is constant along system trajectories. Thus, the solution of Liouville's Equation is

$$
D(q, p, t)=D_{0}[q-\tilde{q}(t), p-\tilde{p}(t)]
$$

We may also write Liouville's Equation in the form of Schrodinger's Equation:

$$
i \frac{\partial D}{\partial t}=i[H, D]=\hat{\Lambda} D
$$

Where the Liouville Operator $\hat{\Lambda}$ is defined as

$$
\hat{\Lambda}=i[H .]=i \sum\left(\frac{\partial H}{\partial q} \frac{\partial}{\partial p}-\frac{\partial H}{\partial p} \frac{\partial}{\partial q}\right)
$$

It can be proved that $\hat{\Lambda}$ is a Hermitian operator. Therefore, it has real eigenvalues and has orthogonal eigenfunctions. Two things are of importance here. First, this property is general for any conservative particle interactions because we did not specifiy the form of interaction here. Second, the eigenvalue of $\hat{\Lambda}$ is real implies $[H . \quad]$ is purely imaginery. Therfore, Liouville Equation must have oscillatory solutions.

# 1.1.5 Eigenfunction Expansions and Resolvent 

With Liouville's Equation in the form of equation 1.25, we multiply through an integration factor to rewrite the equation to be

$$
\frac{\partial}{\partial t}\left(\exp \left(i \int_{0}^{t} d t^{\prime} \hat{\Lambda}\right) D\right)=0
$$

An integration gives solution

$$
D(q, p, t)=\exp \left(-i \int_{0}^{t} d t^{\prime} \hat{\Lambda}\right) D(q, p, 0)
$$

The above solution can be Taylor expanded to achieve the same form as equation 1.22 with a short time interval $\Delta t$. For a longer time interval, we consider the eigenproperties of Liouville operator $\hat{\Lambda}$.

$$
\hat{\Lambda} \psi_{n}=\omega_{n} \psi_{n}
$$

Assuming that these eigenfunctions $\psi_{n}$ span the Hilbert space, we have

$$
D(q, p, 0)=\sum_{n}<\psi_{n} \mid D(q, p, 0)>\psi_{n}
$$

Combining with equation 1.28 , the solution is

$$
D(q, p, t)=\sum_{n}<\psi_{n} \mid D(q, p, 0)>\exp \left(-i \int_{0}^{t} d t^{\prime} \omega_{n}\right) \psi_{n}
$$

We can derive the special case for an ideal gas with $N$ particles in a box of length $L$. Ideal gas has no interaction between particles. Therefore, the Hamiltonian of the particle is purley kinetic. For this case, the Liouville operator $\hat{\Lambda}_{0}$ is

$$
\hat{\Lambda}_{0}=-i \sum_{s} \frac{\partial}{\partial \boldsymbol{p}_{s}}\left(\frac{\boldsymbol{p}_{s}^{2}}{2 m}\right) \cdot \frac{\partial}{\partial \boldsymbol{x}_{s}}=-i \sum_{s} \boldsymbol{v}_{s} \cdot \frac{\partial}{\partial \boldsymbol{x}_{s}}
$$

---

<!-- Page 6 -->

Here, we used velocity $\boldsymbol{v}_{s}=\boldsymbol{p}_{s} / m$. The eigenvalues and eigenvectors can be obtained:

$$
\begin{gathered}
\omega_{(\boldsymbol{k})}=\sum_{s=1}^{N} \boldsymbol{v}_{s} \cdot \boldsymbol{k}_{s} \\
\psi_{(\boldsymbol{k})}=\frac{1}{L^{3 / 2}} \exp \left(i \sum_{s} \boldsymbol{k}_{s} \cdot \boldsymbol{x}_{s}\right)
\end{gathered}
$$

We eigenexpand the initial condition $D\left(\boldsymbol{x}^{N}, \boldsymbol{p}^{N}, 0\right)$ where $\boldsymbol{x}^{N}, \boldsymbol{p}^{N}$ are $3 N$ dimensional vectors

$$
\begin{gathered}
D\left(\boldsymbol{x}^{N}, \boldsymbol{p}^{N}, 0\right)=\sum_{(\boldsymbol{k})} D_{(\boldsymbol{k})}\left(\boldsymbol{p}^{N}\right) \psi_{(\boldsymbol{k})}\left(\boldsymbol{x}^{N}\right) \\
D_{(\boldsymbol{k})}\left(\boldsymbol{p}^{N}\right)=\frac{1}{L^{3 / 2}} \int d \boldsymbol{x}^{N} \exp \left(-i \sum_{s} \boldsymbol{k}_{s} \cdot \boldsymbol{x}_{s}\right) D\left(\boldsymbol{x}^{N}, \boldsymbol{p}^{N}, 0\right)
\end{gathered}
$$

The full solution for $D\left(\boldsymbol{x}^{N}, \boldsymbol{p}^{N}, t\right)$ as ideal gas is

$$
D\left(\boldsymbol{x}^{N}, \boldsymbol{p}^{N}, t\right)=\frac{1}{L^{3 / 2}} \sum_{\langle\boldsymbol{k}\rangle} D_{(\boldsymbol{k})}\left(\boldsymbol{p}^{N}\right) e^{i \sum_{s} \boldsymbol{k}_{s} \cdot\left(\boldsymbol{x}_{s}-\boldsymbol{v}_{s} t\right)}
$$

For the case of a free particle, we have the following relation

$$
-i t \hat{\Lambda}_{0}=-\sum_{s} t \boldsymbol{v}_{s} \cdot \frac{\partial}{\partial \boldsymbol{x}_{s}}
$$

Consder particle 1 in the ensemble, we have the following relation

$$
\exp \left(-t \boldsymbol{v}_{1} \cdot \frac{\partial}{\partial \boldsymbol{x}_{1}}\right) D\left(\boldsymbol{x}_{1}, \boldsymbol{v}_{1}, 0\right)=D\left(\boldsymbol{x}_{1}-t \boldsymbol{v}_{1}, \boldsymbol{v}_{1}, 0\right)
$$

Therefore operator $\exp \left(-i t \hat{\Lambda}_{0}\right)$ propagtes particle backward in time. It is also called free particle operator. At last, Liouville's Equation can be solved in terms of Laplace Transform.

$$
\tilde{D}(s)=\int_{0}^{\infty} d t e^{-s t} D(t)
$$

In the $s$ domain, the Liouville Eqation is

$$
-i \tilde{D}(0)=(\hat{\Lambda}-i s) \tilde{D}(s)
$$

We define the Resolvent operator as $\hat{R}=(\hat{\Lambda}-i s)^{-1}$. The solution of the Liouville Equation can be obatined by taking the inverse transform

$$
D(t)=\frac{1}{2 \pi i} \int_{\gamma-i \infty}^{\gamma+i \infty} d s e^{s t} \tilde{D}(s)
$$

where the line $s=\gamma$ in the complex $s$ plane lies to the right of all singularties of operator $\hat{R}$.

# 1.1.6 Distribution Functions 

Consider a system with constant Hamiltonian with energy $E$. Quantum Mechanics reveals there is an uncertainty around $E$, and this small spread of $E$ generates an energy shell in the $\Gamma$-space. System of particles lies in this energy shell of volume $\Omega$. Following definition 1.18

$$
\mathcal{N}=\int_{\Omega} D d p d q=\text { number of occupie states in } \Omega
$$

In a volume $\Delta \Omega$, the probability that an unspecified state is occupied in $\Delta \Omega$ is

$$
\frac{\Delta \Omega}{\Omega}=\frac{\int_{\Delta \Omega} D d p d q}{\int_{\Omega} D d p d q}=\int_{\Delta \Omega} f_{N} d p d q
$$

with quantity $f_{N}$ defined as above. If we take $\Delta \Omega$ to be infinitesimal, we can see that

$$
f_{N}(q, p, t)=C D(q, p, t)
$$

---

<!-- Page 7 -->

where $C$ is a constant. Therefore, it is also obvious that function $f_{N}(q, p, t)$ also solves Liouville's Equation. The function $f_{N}(q, p, t)$ has the physical significance of joint probability of particle $n$ on the point $\left(\boldsymbol{x}_{n}, \boldsymbol{p}_{n}\right)$ for $1 \leqslant n \leqslant N$ particles in the system. $f_{N}$ is called the $N$-body joint probability density for the $N$ body system. Comparing with $D$ as an ensemble of system points, $f_{N}$ describes a single system. As a probability distribution, $f_{N}$ must be normalized

$$
\int_{\Omega} f_{N} d p d q=1
$$

It follows if $G$ is any dynamical variable, the expectation or average of $G$ is

$$
\langle G\rangle=\int_{\Omega} f_{N} G d p d q
$$

We introduce the notation

$$
d s=d \boldsymbol{x}_{s} d \boldsymbol{p}_{s}
$$

For a system of $N$ particles. Consider subsystem with $s<N$. The probability of finding a subsystem in the phasevolume $d 1 d 2 \ldots . . d s$ from probability distribution $f_{N}(1,2, \ldots \ldots, N)$ is

$$
f_{s}(1,2, \ldots \ldots, s) d 1 d 2 \ldots . . d s
$$

To construct $f_{s}$ from $f_{N}$, the information on state of particles $s+1, \ldots \ldots, N$ must be integrated out from $f_{N}(1,2, \ldots \ldots, N)$, which gives

$$
f_{s}(1,2, \ldots \ldots, s)=\int f_{N}(1,2, \ldots \ldots, N) d(s+1) \ldots \ldots d N
$$

This is the reduced distribution. Apart from joint-probability distributions, we also encounter conditional probability distributions. For a particle 1 is phase volume $d 1$ about the phase point 1 , given particles $2 \ldots . . N$ are in volume $d 2 \ldots . . d N$, the probability distribution is given by the Bayes' Formula. This could be also generalized to multiple particles.

$$
h_{n}(1,2 \ldots n \mid n+1 \ldots N)=\frac{f_{N}(1, \ldots, N)}{f_{N-n}(n+1, \ldots, N)}
$$

The $s$ - tuple distribution is defined as following. The below product represents the probable number of $s$ - tuple of particles such that one of the particles is in volume $d 1$ about point 1 , another in volume $d 2$ about point 2 and so on.

$$
F_{s}(1, \ldots, s) d 1 \ldots d s
$$

The relation between $F_{s}$ and $f_{s}$ is

$$
F_{s}=s!\frac{N!}{s!(N-s)!} f_{s}=\frac{N!}{(N-s)!} f_{s}
$$

The average kinetic and potential energies are best written in the form of $F_{s}$. Since $F_{s}$ does not care about state of individual particles, it must be symmetric. Therefore, it follows $f_{s}$ must also be properly symmetrized.

# 1.1.7 Markov Process 

Introduce the phase vector $\boldsymbol{z}=(1,2, \ldots, N)$. In terms of phase vectors, the conditional distribution functions may be written in the form of products

$$
\prod\left(z, t \mid z_{0}, t_{0}\right) d z
$$

which represents the probability in finding the system in the state $d z$ granted it was in state $z_{0}$ at time $z_{0}$. There are 3 properties that has to be satisfied for this two-time distribution.

1. System evolves and has intermediate states between $\left(t-t_{0}\right)$. The probability is normalized to 1 for all points.

$$
\int_{\forall z} \prod\left(z, t \mid z_{0}, t_{0}\right) d z=1
$$

---

<!-- Page 8 -->

2. System does not change for a time interval of zero.

$$
\prod\left(\boldsymbol{z}, t_{0} \mid \boldsymbol{z}_{0}, t_{0}\right)=\delta\left(\boldsymbol{z}-\boldsymbol{z}_{0}\right)
$$

3. $\prod\left(\boldsymbol{z}, t \mid \boldsymbol{z}_{0}, t_{0}\right)$ obeys Chapman-Kolmogorov Equation.

$$
\prod\left(\boldsymbol{z}, t_{0} \mid \boldsymbol{z}_{0}, t_{0}\right)=\int \prod\left(\boldsymbol{z}, t_{0} \mid \boldsymbol{z}^{\prime}, t^{\prime}\right) \prod\left(\boldsymbol{z}^{\prime}, t^{\prime} \mid \boldsymbol{z}_{0}, t_{0}\right)
$$

The Chapman-Kolmogorov contains an apprxomation such that the system for a later state only depends on the state immediately before it. Such process is called Markov Process. Further noted that this particular path from $\boldsymbol{z}_{0}$ to $\boldsymbol{z}$ is uncertain. When the product $\prod\left(\boldsymbol{z}, t_{0} \mid \boldsymbol{z}_{0}, t_{0}\right)$ equals the joint probability distribution $f_{N}$, the path is unique when solved by Liouville's Equation. Therefore, a differential equation equivalent to Liouville Equation for $\prod\left(\boldsymbol{z}, t_{0} \mid \boldsymbol{z}_{0}, t_{0}\right)$ derived from CK equation could be applied to random or quantum processes more effectively. For a system whose behavior is homogeneous in time and for a countable set, we have the equivalent expression for CK equation where $\{z\} \rightarrow\{l\}, l=0,1,2,3, \ldots$

$$
\prod\left(l \mid l_{0} ; t+\Delta t\right)=\sum_{\forall j} \prod(l \mid j ; \Delta t) \prod\left(j \mid l_{0} ; t\right)
$$

We can derive a canonical form for master equation

$$
\frac{\partial}{\partial t} \prod\left(l \mid l_{0} ; t\right)=\sum_{\forall j}\left[w_{l j} \prod\left(j \mid l_{0} ; t\right)-w_{j l} \prod\left(l \mid l_{0} ; t\right)\right]
$$

Where $w_{l j}$ denotes the probability rate at which transition from state $j$ to state $l$ occur. The above equation could be applied to discussion of random walks where transitional probability rate is

$$
w_{j l}=\frac{1}{2}\left(\delta_{j, l+1}+\delta_{j, l-1}\right)
$$

Here, the particle moved through displacement $l$ in total $n$ step. The probability of being at position $l$ in total $n$ step is

$$
P(l, n)=\frac{1}{2}[P(l+1, n+1), P(l-1, n-1)]-P(l, n-1)
$$

The moments could be further computed

$$
\begin{gathered}
M_{1}(n)=\sum_{l} l P(n, l)=0 \\
M_{2}(n)=\sum_{l} l^{2} P(n, l)=n
\end{gathered}
$$

# 1.1.8 Central Limit Theorem 

A random variable is a mapping $\xi$ from a sample space to the real line. We can associate a probability $P(x)$ with this mapping $\xi$. A fundamental property of it would be normalization $\sum P(\xi)=1$. It has also expectation $\langle\xi\rangle$ and variance $D(\xi)=\left\langle\xi^{2}-\langle\xi\rangle^{2}\right\rangle$. The characteristic function $\phi(a)$ for $P(x=\xi)$ where $x$ are the elements of the real line is

$$
\phi(a)=\sum_{\forall x} P(x) e^{i a x}
$$

It could be observed that $x$ and $a$ form a Fourier pair. The the first derivative of $\ln \phi(a)$ with respect to $a$ gives $i<\xi>$, whereas the second derivative gives $-D(\xi)$. We have the Taylor expansion form for $\ln \phi(a)$.

$$
\ln \phi(a)=0+i<\xi>a-\frac{1}{2} D(\xi) a^{2}+\cdots
$$

For statistically independent $n$ random variables $\xi_{r}$ with $r=1,2, \ldots n$, we have

$$
\phi(a)=\prod_{r=1}^{n} \phi_{r}\left(\xi_{r}\right), \xi=\sum_{r=1}^{n} \xi_{r}
$$

---

<!-- Page 9 -->

The central limit theorem address the siuation when $n \rightarrow \infty$ and the probability distribution are all the same for $n$ variables. We can find the expectation and variance for $n$ random variables by differentiating $\ln \phi(a)=\ln \prod_{r=1}^{n} \phi_{r}\left(\xi_{r}\right)$.

$$
\begin{gathered}
<\xi>(\xi)=-i \frac{d}{d a} \ln \prod_{r=1}^{n} \phi_{r}\left(\xi_{r}\right)=\sum_{r=1}^{n}<\xi_{r}>\left(\xi_{r}\right)=n<\xi_{r}> \\
D(\xi)=-\frac{d^{2}}{d a^{2}} \ln \prod_{r=1}^{n} \phi_{r}\left(\xi_{r}\right)=\sum_{r=1}^{n} D\left(\xi_{r}\right)=n D\left(\xi_{r}\right)
\end{gathered}
$$

Now we can take the inverse Fourier Transform for $\phi(a)$ to find the probability distribution for large $n$, where $l$ below is the displacement for random walk. Notice below, only when $n \rightarrow \infty$, the integration limit for Fourier Transform could be taken from $-\infty$ to $\infty$.

$$
\begin{gathered}
P(l, n)=\frac{1}{2 \pi} \int_{-\infty}^{\infty} d a e^{-i a l} \phi(a) \\
=\frac{1}{2 \pi} \int_{-\infty}^{\infty} d a e^{-i a l} \exp \left(i a<\xi>(\xi)-\frac{1}{2} D(\xi) a^{2}+\cdots\right)
\end{gathered}
$$

We have result

$$
P(l, n)=\frac{1}{(2 \pi D)^{1 / 2}} \exp \left\{\frac{(l-<\xi>)^{2}}{2 D}\right\}
$$

For a large $n$, introduce a continous variable $x=l \delta$ for number of intervals $\delta$. The probility that the net displacement occur in the interval $(x, x+\Delta x)$ at the $n$ th step is $P(x, n) \Delta x=P(l, n) \Delta x / \delta$. Therefore, for a continous limit, we have

$$
P(x, t)=\frac{1}{(2 \pi \sigma)^{1 / 2}} e^{-x^{2} / 2 \sigma^{2}}
$$

where $\sigma^{2}=\delta^{2} t$. The above Gaussian Distribution maintains normalization. The variance $D(x)$ is occur in the expression as $D(x)=\sigma^{2}$. Hence, $\sigma$ is standard deviation of the probability distribution. Notice that $\sigma^{2} \propto t$, which coincides with the result for random walk 1.58.

# 1.2 Analysis of Liouville's Equation 

### 1.2.1 BBGKY Hierachy

We can rewrite the Liouville Equation in terms of another operator $\hat{L}_{N}$

$$
\frac{\partial f_{N}}{\partial t}-\left[\hat{L}_{N}, f_{N}\right]=0
$$

where this new Liouville Operator can be casted in terms of Hamiltonian $H$

$$
\begin{gathered}
H=\sum_{i=1}^{N} \frac{p_{i}^{2}}{2 m}+\sum_{i<j}^{N} \sum \Phi\left(\left|\boldsymbol{x}_{i}-\boldsymbol{x}_{j}\right|\right) \\
\hat{L}_{N}=-\sum_{l=1}^{N} \frac{\boldsymbol{p}_{l}}{m} \cdot \frac{\partial}{\partial \boldsymbol{x}_{l}}+\sum_{i<j}^{N} \sum \sum_{l=1}^{N} \frac{\partial}{\partial \boldsymbol{x}_{l}} \Phi\left(\left|\boldsymbol{x}_{i}-\boldsymbol{x}_{j}\right|\right) \cdot \frac{\partial}{\partial \boldsymbol{p}_{l}}
\end{gathered}
$$

For a two particle interaction potential $\Phi\left(\left|\boldsymbol{x}_{i}-\boldsymbol{x}_{j}\right|\right)=\Phi_{i j}$. The first term implies a kinetic energy operator $\hat{K}_{l}$. We can simplify the 2nd term by noticing the terms in the summation is only nonzero if $l=i$ or $l=j$. Therefore, we can define operator $\hat{O}_{i j}$ for these terms

$$
\hat{O}_{i j}=\frac{\partial}{\partial \boldsymbol{x}_{i}} \Phi\left(\left|\boldsymbol{x}_{i}-\boldsymbol{x}_{j}\right|\right) \cdot \frac{\partial}{\partial \boldsymbol{p}_{i}}+\frac{\partial}{\partial \boldsymbol{x}_{j}} \Phi\left(\left|\boldsymbol{x}_{i}-\boldsymbol{x}_{j}\right|\right) \cdot \frac{\partial}{\partial \boldsymbol{p}_{j}}
$$

---

<!-- Page 10 -->

For a two particle potential, its force $\boldsymbol{G}_{i j}$ is antisymmetric under particle exchange. Therefore, we have the equality

$$
\boldsymbol{G}_{i j}=-\frac{\partial}{\partial \boldsymbol{x}_{i}} \Phi_{i j}=\frac{\partial}{\partial \boldsymbol{x}_{j}} \Phi_{i j}
$$

We can further simplify $\hat{O}_{i j}$ to be

$$
\hat{O}_{i j}=\frac{\partial}{\partial \boldsymbol{x}_{i}} \Phi_{i j} \cdot\left(\frac{\partial}{\partial \boldsymbol{p}_{i}}-\frac{\partial}{\partial \boldsymbol{p}_{j}}\right)=-\boldsymbol{G}_{i j} \cdot\left(\frac{\partial}{\partial \boldsymbol{p}_{i}}-\frac{\partial}{\partial \boldsymbol{p}_{j}}\right)
$$

Now we have a compact form for $\hat{L}_{N}$

$$
\hat{L}_{N}=-\sum_{l=1}^{N} \hat{K}_{l}+\sum_{i<j}^{N} \sum \hat{O}_{i j}
$$

we are interested in an operator for reduced distribution $f_{s}$ where $s<N$. Hence, the operator $\hat{L}_{N}$ can be partitioned

$$
\hat{L}_{N}=\hat{L}_{s}+\hat{L}_{N, s+1}
$$

where the reduced s-particle Liouville Operator is

$$
\hat{L}_{s}=-\sum_{l=1}^{s} \hat{K}_{l}+\sum_{i<j}^{s} \sum \hat{O}_{i j}
$$

and the remainder Liouville Operator is

$$
\hat{L}_{N, s+1}=-\sum_{l=1}^{s} \hat{K}_{l}+\sum_{i=1}^{s} \sum_{j=s+1}^{N} \hat{O}_{i j}+\sum_{i=s+1}^{N} \sum_{j=s+1(i<j)}^{N} \hat{O}_{i j}
$$

With 1.48 and 1.73, we can rewrite the Liouville Equation into

$$
\left(\frac{\partial}{\partial t}-\hat{L}_{s}\right) f_{s}=\int \mathrm{d}(s+1) \ldots \mathrm{dN}\left(\hat{L}_{N, s+1}\right) f_{N}
$$

On the right hand side, the intergal is

$$
\int \mathrm{d}(s+1) \ldots \mathrm{dN}\left(-\sum_{l=1}^{s} \hat{K}_{l}+\sum_{i=1}^{s} \sum_{j=s+1}^{N} \hat{O}_{i j}+\sum_{i=s+1}^{N} \sum_{j=s+1(i<j)}^{N} \hat{O}_{i j}\right) f_{N}
$$

The first and the third terms are surface terms and will vanish. Equation 1.75 becomes

$$
\left(\frac{\partial}{\partial t}-\hat{L}_{s}\right) f_{s}=-\sum_{i=1}^{s} \frac{\partial}{\partial \boldsymbol{p}_{i}} \cdot \int \mathrm{~d}(s+1) \ldots \mathrm{dN} \sum_{j=s+1}^{N} \boldsymbol{G}_{i j} f_{N}
$$

The above equation can be further simplified for a distribution $f_{N}$ that is symmetric under the exchange of particles. In such fluid, the particles would be identical. This makes the terms in the above 2nd summation along $j$ to be the same for all $j$. Then, we set this index to be $s+1$ and the rest integral with repect to $s+2 \ldots . N$ becomes 1 . Now, we have

$$
\left(\frac{\partial}{\partial t}-\hat{L}_{s}\right) f_{s}=-(N-s) \sum_{i=1}^{s} \frac{\partial}{\partial \boldsymbol{p}_{i}} \cdot \int \mathrm{~d}(s+1) \boldsymbol{G}_{i, s+1} f_{s+1}
$$

The above $N$ coupled equations are called BBGKY equations. They are also called hierachy, abbreviated by $\mathrm{BY}_{s}$ for the $s$ th equation in the hierachy. There are 3 basic properties for these equations

1. The $N$ th equation in the hierachy is the Liouville Equation.
2. The LHS of these equation does not equal to zero means they are not constant along trajectories in the $\tilde{\Gamma}$-space. They do not vanishes due to the interaction represented by the RHS of these equation.

---

<!-- Page 11 -->

3. The first equation $\mathrm{BY}_{1}$ is the generic form of all kinetic equations

$$
\frac{\partial f_{1}}{\partial t}+\frac{\boldsymbol{p}_{1}}{m} \cdot \frac{\partial}{\partial \boldsymbol{x}_{1}} f_{1}=-\hat{A}_{1} f_{2}(1,2)=\hat{J} f_{1}(1)
$$

If the second equality is true by a transformation to close $\mathrm{BY}_{1}, \mathrm{BY}_{1}$ will have the form of kinetic equations, where $\hat{J}$ is often called the collision integral.
The first equation and the second equation $\mathrm{BY}_{1}$ and $\mathrm{BY}_{2}$ can be used to retrieve the conservation of energy.

# 1.2.2 Correlation Expansions: Vlasov Limit 

We first non-dimensionalize $\mathrm{BY}_{s}$. For a fluid with $N$ particles in a volume $V$, the number density $n_{0}$ is

$$
n_{0}=\frac{N}{V}
$$

For this fluid, we can assign a mean thermal speed $C$ and a corresponding temperature $T$

$$
m C^{2}=k_{B} T
$$

where $k_{B}$ is Boltzmann's constant. The strength of potential $\Phi_{0}$ and a characteristic scale length $r_{0}$ are defined by

$$
G_{i j}=\frac{\Phi_{0}}{r_{0}} \hat{G}_{i j}
$$

Now, the force $\hat{G}_{i j}$ is non-dimensional. We renormalize $f_{s}$ to $F_{s}$ such that the integral of $F_{1}$ gives the number density.

$$
F_{s}=V^{s} f_{s}
$$

The full non-dimensionalization scheme is

$$
\begin{array}{ll}
x=r_{0} \tilde{x} & p=m C \tilde{p} \\
t=\frac{r_{0}}{C} \tilde{t} & F_{s}=(m C)^{-3 s} \tilde{F}_{s}
\end{array}
$$

The hierachy can be rewritten in terms of these variables:

$$
\left(\frac{\partial}{\partial t}+\sum_{l=1}^{s} \hat{K}_{l}-\alpha \sum_{i<j}^{s} \sum_{\hat{O}_{i j}}\right) F_{s}=\frac{\alpha}{\gamma} \hat{I}_{s} F_{s+1}
$$

where operator $\hat{I}_{s}$, coefficient $\alpha$ and $\gamma$ are defined as

$$
\begin{gathered}
\hat{I}_{s}=-N \sum_{i=1}^{s} \frac{\partial}{\partial \hat{\boldsymbol{p}}_{i}} \cdot \int \mathrm{~d}(s+1) \boldsymbol{G}_{i, s+1} \\
\alpha=\frac{\Phi_{0}}{k_{B} T} \\
\gamma^{-1}=n_{0} r_{0}^{3}
\end{gathered}
$$

Consider the coefficient $\frac{\alpha}{\gamma}$, if $r_{0}$ characterize the length of the entire volume, the number of particles $\mathcal{N}_{0}$ in a sphere would be $\mathcal{N}_{0}=n_{0} V=n_{0} r_{0}^{3}$. This factor becomes the ratio of potential energy to kinetic energy in the volume.

$$
\frac{\alpha}{\gamma}=\frac{\mathcal{N}_{0}^{2} \Phi_{0}}{\mathcal{N}_{0} k_{B} T} \sim \frac{\left\langle E_{\Phi}\right\rangle}{\left\langle E_{K}\right\rangle}
$$

This terms therefore characterize the relative strength of interaction of particles. If $\alpha / \gamma \geq 1$, the fluid is strongly coupled. If $\alpha / \gamma \ll 1$, the fluid is weakly coupled. Returning to the discussion of correlations, for two statistically dependent particles (1,2), we can write their correlation function $C_{2}(1,2)$ as

$$
f_{2}(1,2)=f_{1}(1) f_{1}(2)+C_{2}(1,2)
$$

---

<!-- Page 12 -->

If the two particles are statistically independent, $C_{2}=0$ and we say they are uncorrelated. We will use $C_{n}(1, \ldots n)$ to describe the correlation function between $n$ particles. First, consider the limit for long range interaction where $\alpha=\frac{\Phi_{0}}{k_{B} T} \ll 1$ and $\gamma^{-1}=n_{0} c_{0}^{3} \gg 1$. In this case, the correlation between particles is also small. Setting $\alpha \xrightarrow{\gamma} \alpha, \gamma^{-1} \rightarrow \varepsilon^{-1} \gamma^{-1}$ and

$$
\begin{gathered}
f_{2}=f_{1} f_{1}+\varepsilon C_{2} \\
f_{3}=f_{1} f_{1} f_{1}+\varepsilon \sum_{P} f_{1} C_{2}+\varepsilon^{2} C_{3}
\end{gathered}
$$

Returning to using $F_{s}$, we obtain new forms of $\mathrm{BY}_{1}$ and $\mathrm{BY}_{2}$ from equation 1.79

$$
\begin{gathered}
\left(\frac{\partial}{\partial t}+\hat{\kappa}_{1}\right) F_{1}=\frac{\alpha}{\gamma} \hat{I}_{1}\left[F_{1}(1) F_{1}(2)+\varepsilon C_{2}(1,2)\right] \\
\left(\frac{\partial}{\partial t}+\hat{\kappa}_{2}-\varepsilon \alpha \hat{O}_{12}\right)\left[F_{1}(1) F_{1}(2)+\varepsilon C_{2}(1,2)\right]=\frac{\alpha}{\gamma} \hat{I}_{2}\left[F_{1}(1) F_{1}(2) F_{1}(3)+\varepsilon F_{1}(1) C_{2}(2,3)\right. \\
\left.+\varepsilon F_{1}(2) C_{2}(1,3)+\varepsilon F_{1}(3) C_{2}(1,2)+\varepsilon^{2} C_{3}(1,2,3)\right)
\end{gathered}
$$

Here, we defined $\hat{\kappa}_{s}=\sum_{l=1}^{s} \hat{K}_{l}$.Similarly, if all the terms with higher order of $\varepsilon$ can be discarded, we have a general form for $\mathrm{BY}_{s}$

$$
\left(\frac{\partial}{\partial t}+\hat{\kappa}_{s}\right) \prod_{i=1}^{s} F_{1}(i)=\frac{\alpha}{\gamma} \hat{I}_{1} \prod_{i=1}^{s+1} F_{1}(i)
$$

There are $N$ equations for one single unknown $F_{1}$. Therefore, these $N$ equations are redundant. To establish the redundancy, we must have

$$
\sum \prod_{i} F_{1}\left[\left(\frac{\partial}{\partial t}+\hat{K}_{l}\right) F_{1}(l)-\frac{\alpha}{\gamma} \hat{I}_{s} F_{1}(l) F_{1}(s+1)\right]=0
$$

The only equation we have left is the first equation in equation 1.83, which now is written in dimensional form:

$$
\left(\frac{\partial}{\partial t}+\boldsymbol{v} \cdot \frac{\partial}{\partial x}\right) F(\boldsymbol{x}, \boldsymbol{v}, t)=-\frac{n_{0}}{m} \frac{\partial}{\partial \boldsymbol{v}} \cdot \int \mathrm{~d} \boldsymbol{x}^{\prime} \mathrm{d} \boldsymbol{v}^{\prime} \boldsymbol{G}\left(\boldsymbol{x}, \boldsymbol{x}^{\prime}\right) F(\boldsymbol{x}, \boldsymbol{v}, t) F\left(\boldsymbol{x}^{\prime}, \boldsymbol{v}^{\prime}, t\right)
$$

This is the Vlasov Equation. Defining the number density

$$
n\left(\boldsymbol{x}^{\prime}, t\right)=n_{0} \int F\left(\boldsymbol{x}^{\prime}, \boldsymbol{v}^{\prime}, t\right) \mathrm{d} \boldsymbol{v}^{\prime}
$$

and the mean field strength

$$
\boldsymbol{G}(\boldsymbol{x}, t)=\int n\left(\boldsymbol{x}^{\prime}, t\right) \boldsymbol{G}\left(\boldsymbol{x}, \boldsymbol{x}^{\prime}\right) \mathrm{d} \boldsymbol{x}^{\prime}
$$

The equation can be rewritten as

$$
\left(\frac{\partial}{\partial t}+\boldsymbol{v} \cdot \frac{\partial}{\partial x}+\frac{\boldsymbol{G}}{m} \cdot \frac{\partial}{\partial \boldsymbol{v}}\right) F(\boldsymbol{x}, \boldsymbol{v}, t)=0
$$

The equation describe particle distribution $F(\boldsymbol{x}, \boldsymbol{v}, t)$ developing under the instantaneous average of all two particle forces in the fluid. This equation has similar form to the one particle Liouville Equation. The solution of above equation can be obtained by the method of characteristics, which is given by

$$
F=F\left[\frac{p^{2}}{2 m}+\Phi(\boldsymbol{x})\right]
$$

where $\boldsymbol{G}(\boldsymbol{x})=-\frac{\partial}{\partial \boldsymbol{x}} \Phi(\boldsymbol{x})$ and $\Phi(\boldsymbol{x})$ is a known function.

# 1.2.3 Prigogine Analysis 

Now we introduce a perturbation solution to the Liouville's Equation. Define the perturbation as

$$
\hat{L}_{N}=\hat{L}_{0}+\delta \hat{L}
$$

---

<!-- Page 13 -->

where $\hat{L}_{0}$ is a free particle operator and the perturbation to it is

$$
\delta \hat{L}=\sum_{l<j}^{N} \sum \frac{\partial}{\partial \boldsymbol{x}_{i}} \Phi_{i j} \cdot\left[\frac{\partial}{\partial \boldsymbol{p}_{i}}-\frac{\partial}{\partial \boldsymbol{p}_{j}}\right]
$$

The eigenfunctions of free particles operator were discussed before:

$$
\begin{gathered}
\psi_{(\boldsymbol{k})}=L^{-3 N / 2} \exp \left(i \sum \boldsymbol{k}_{l} \cdot \boldsymbol{x}_{l}\right) \equiv|(\boldsymbol{k})\rangle \\
\omega_{(\boldsymbol{k})}=\sum \boldsymbol{k}_{l} \cdot \boldsymbol{v}_{l}
\end{gathered}
$$

where $\boldsymbol{v}_{l}=\boldsymbol{p}_{l} / m$. In this case, the solution with respect to operator $\hat{L}_{N}$ can be written in the form of eigenfunction expansion

$$
f_{N}(1, \ldots, N)=\sum_{\langle\boldsymbol{k}\rangle} a_{(\boldsymbol{k})}\left(\boldsymbol{p}^{N}, t\right) \psi_{(\boldsymbol{k})}\left(\boldsymbol{x}^{N}\right) e^{-i \omega_{(\boldsymbol{k})} t}
$$

The terms in the series can be regrouped such that coefficients in the new series have special coefficients.

$$
\begin{gathered}
f_{N}=\frac{1}{V^{N}}\left[a_{0}\left(\boldsymbol{p}^{N}, t\right)+\frac{1}{V} \sum_{j=1}^{N} \sum_{\boldsymbol{k}(j)} a_{1}\left(\boldsymbol{k}_{j}, \boldsymbol{p}^{N}, t\right) e^{i \boldsymbol{k}_{j} \cdot \boldsymbol{x}_{i}} e^{-i \omega_{j} t}\right] \\
+\frac{1}{V^{N}}\left[\frac{1}{V^{2}} \sum \sum_{j<l} \sum_{\boldsymbol{k}_{j} \boldsymbol{k}_{l}}\left[\left(1-\delta_{\boldsymbol{k}_{j}+\boldsymbol{k}_{l}}\right)+\delta_{\boldsymbol{k}_{j}+\boldsymbol{k}_{l}}\right] a_{2}\left(\boldsymbol{k}_{j}, \boldsymbol{k}_{l}, \boldsymbol{p}^{N}, t\right) e^{i\left[\boldsymbol{k}_{j} \cdot \boldsymbol{x}_{j}+\boldsymbol{k}_{l} \cdot \boldsymbol{x}_{l}\right]} e^{-i \omega_{j} t}\right]
\end{gathered}
$$

where $\hat{V}=V /(2 \pi)^{3}$. The integral of coefficients $a_{0}, a_{1}, a_{2}$ has the special significance as the following

$$
\begin{gathered}
1=\int a_{0}\left(\boldsymbol{p}^{N}, t\right) \mathrm{d} \boldsymbol{p}^{N} \\
n_{1}(\boldsymbol{x})=\frac{N}{V}\left[1+\int a_{1}\left(\boldsymbol{k}, \boldsymbol{p}^{N}, t\right) e^{i \boldsymbol{k} \cdot \boldsymbol{x}} \mathrm{~d} \boldsymbol{k} \mathrm{~d} \boldsymbol{p}^{N}\right] \\
n_{2}\left(\boldsymbol{x}, \boldsymbol{x}^{\prime}\right)=\frac{N(N-1)}{2 V^{2}}\left[1+\int a_{2}\left(\boldsymbol{k},-\boldsymbol{k}, \boldsymbol{p}^{N}, t\right) e^{i \boldsymbol{k} \cdot\left(\boldsymbol{x}-\boldsymbol{x}^{\prime}\right)} \mathrm{d} \boldsymbol{k} \mathrm{~d} \boldsymbol{p}^{N}\right]
\end{gathered}
$$

Therefore, $a_{0}$ has the meaning of $N$ particle momentum distribution. For a homogeneous system, the intergals in the 2 nd and 3 rd equation goes to zero. In this case, $n_{1}$ has the significance of number density and $n_{2}$ has the significance of pair distribution. With the reorganized form of $f_{N}$ in equation 1.93, we can obtain equations of motion for these coefficients $a_{(\boldsymbol{k})}$. First we have to Fourier Expand the potential in the perturbed operator

$$
\begin{gathered}
\delta \hat{L}=\sum_{l<n}^{N} \sum \frac{\partial}{\partial \boldsymbol{x}_{l}} \Phi_{l n} \cdot\left[\frac{\partial}{\partial \boldsymbol{p}_{l}}-\frac{\partial}{\partial \boldsymbol{p}_{n}}\right] \\
=\sum_{l<n}^{N} \sum \sum_{\boldsymbol{K}} \phi_{\boldsymbol{K}} i \boldsymbol{K} e^{i \boldsymbol{K} \cdot\left(\boldsymbol{x}_{l}-\boldsymbol{x}_{n}\right)} \cdot\left[\frac{\partial}{\partial \boldsymbol{p}_{l}}-\frac{\partial}{\partial \boldsymbol{p}_{n}}\right] \\
\equiv \sum_{l<n}^{N} \sum \sum_{\boldsymbol{K}} \phi_{\boldsymbol{K}} i \boldsymbol{K} e^{i \boldsymbol{K} \cdot\left(\boldsymbol{x}_{l}-\boldsymbol{x}_{n}\right)} \cdot \hat{\boldsymbol{\theta}}_{l n}
\end{gathered}
$$

Inserting equation 1.93 into equation 1.66

$$
\begin{gathered}
\frac{\partial}{\partial t} a_{(\boldsymbol{k})}=\sum_{\langle\boldsymbol{k}^{\prime}\rangle} e^{i \omega_{(\boldsymbol{k})} t}\langle\langle\boldsymbol{k})| \delta \hat{L}|\langle\boldsymbol{k}^{\prime}\rangle\rangle\left(e^{i \omega_{\left(\boldsymbol{k}^{\prime}\right)}^{t} t} a_{\left(\boldsymbol{k}^{\prime}\right)}\right) \\
=\frac{1}{V} \sum_{\left\langle\boldsymbol{k}^{\prime}\right\rangle} \sum_{l<n}^{N} \sum \sum_{\boldsymbol{K}} e^{i \omega_{(\boldsymbol{k})} t}\left|\langle\boldsymbol{k}\rangle\right| e^{i \boldsymbol{K} \cdot\left(\boldsymbol{x}_{l}-\boldsymbol{x}_{n}\right)}\left|\left(\boldsymbol{k}^{\prime}\right\rangle\right\rangle\rangle \phi_{\boldsymbol{K}} i \boldsymbol{K} \cdot \hat{\boldsymbol{\theta}}_{l n}\left(e^{i \omega_{\left(\boldsymbol{k}^{\prime}\right)}^{t} t} a_{\left(\boldsymbol{k}^{\prime}\right)}\right)
\end{gathered}
$$

---

<!-- Page 14 -->

the matrix element can be evaluated

$$
\left\langle(\boldsymbol{k})\right| e^{i \boldsymbol{K} \cdot\left(\boldsymbol{x}_{l}-\boldsymbol{x}_{n}\right)}\left|\left(\boldsymbol{k}^{\prime}\right)\right\rangle=\delta\left(\boldsymbol{K}-\boldsymbol{k}_{l}-\boldsymbol{k}_{l}^{\prime}\right) \delta\left(\boldsymbol{k}_{\mathrm{n}}^{\prime}-\boldsymbol{k}_{n}-\boldsymbol{K}\right) \prod_{j \neq l, n}^{N} \delta\left(\boldsymbol{k}_{j}^{\prime}-\boldsymbol{k}_{j}\right)
$$

# 1.2.4 Bogoliubov Hypothesis 

Bogoliubov Hypothesis address the behavior of an enclosed gas in three different time intervals.

1. Interval $\tau_{1}$ : The amount of time for two molecules in each other's interaction domain.
2. Interval $\tau_{2}$ : Mean time interval between collisions.
3. Interval $\tau_{3}$ : The amount of time for a molecule to traverse the entire closed domain.

For a mole of gas in a macroscopic container, the following order of magnitude is true:

$$
\tau_{1} \ll \tau_{2} \ll \tau_{3}
$$

Three corresponding displacement can be labeled base on these three intervals:

1. Displacement $\lambda_{1}$ : range of interaction
2. Displacement $\lambda_{2}$ : mean free path
3. Displacement $\lambda_{3}$ : edge length of the entire domain

In the three time intervals, the behavior of molecules can be described in three stages.

1. The initial stage $0<t \leq \tau_{1}$ : There exists no collisions between molecules. No equilibrating force exists. Newtonian mechanics require the information at least $3 N$ independent variables to track the motion of individual particles.
2. The kinetic stage $\tau_{1}<t \leq \tau_{2}$ : Collisions exist and there is a tendency towards equilibrium. Hypothesis: The particles are interrelated and are best seen by particles' probability distributions. All $s$ - particle distribution is a functional of $f_{1}$

$$
f(s)=f_{s}\left(1 \ldots s, f_{1}\right)
$$

with time dependence full contained in $f_{1}$
3. The hydrodynamic stage $\tau_{2}<t$ : Equilibrium reached. Hypothesis: all particle distributions depend only on macroscopic variables for fluids. (density, macroscopic velocity, temperatures)
The initial stage essentially formulates an $n$-body problem, whereas the hydrodynamic stage could be described in terms of fluid dynamical equations and thermodynamics. We are particularly interested in the behavior of kinetic stage. Consider the hierachy in the thermodynamical limit

$$
\lim _{N \rightarrow \infty, V \rightarrow \infty}\left(\frac{N-s}{N}\right)=\lim \left(\frac{N}{V}\right)=\frac{1}{v}
$$

where $v$ is the specific volume. In this limit, the hierachy becomes

$$
\left(\frac{\partial}{\partial t}-\hat{L}_{s}\right) F_{s}-\frac{1}{v} \sum_{i=1}^{s} \int \mathrm{~d}(s+1) \hat{O}_{i, s+1} F_{s+1}=0
$$

Here, we reuse the renormalized distribution $F_{s}$ from 1.79. Since the hypothesis in 1.97 must also hold true for $F_{s}$, we can write

$$
\frac{\partial F_{s}}{\partial t}=\left\langle\frac{\delta F_{s}}{\delta F_{1}}, \frac{\partial F_{1}}{\partial t}\right\rangle
$$

The above bracket denotes inner product with respect to all $F_{1}$ functions. For statistically independent molecules $F_{s}=\prod_{i=1}^{s} F_{1}(s)$, hence

$$
\frac{\partial F_{s}}{\partial t}=\sum_{k} \frac{\partial F_{1}(k)}{\partial t} \prod_{l \neq k} F_{1}(l)
$$

---

<!-- Page 15 -->

To discuss the behavior of rare gas, the following expansion is introduced

$$
F_{s}\left(1, \ldots s ; F_{1}\right)=F_{s}^{(0)}+\frac{1}{v} F_{s}^{(1)}+\frac{1}{v^{2}} F_{s}^{(2)}+\ldots
$$

Substituting this for $F_{2}$ in $\mathrm{BY}_{1}$ (equation 1.78)

$$
\begin{gathered}
\frac{\partial F_{1}}{\partial t}=-\frac{\boldsymbol{p}_{1}}{m} \cdot \frac{\partial F_{1}}{\partial \boldsymbol{x}_{1}}+\frac{1}{v} \hat{I}_{12}\left(F_{2}^{(0)}+\frac{1}{v} F_{2}^{(1)}+\cdots\right) \\
\equiv A^{(0)}+\frac{1}{v} A^{(1)}+\cdots \\
\hat{I}_{12}=\int \mathrm{d} 2 \hat{O}_{12}
\end{gathered}
$$

Using the inner product notations in equation 1.99 , to obtain an equation for entire hierachy

$$
\begin{gathered}
\frac{\partial F_{s}}{\partial t}=\left\langle\frac{\delta F_{s}^{(0)}}{\delta F_{1}}, A^{(0)}\right\rangle+\frac{1}{v}\left\{\left\langle\frac{\delta F_{s}^{(1)}}{\delta F_{1}}, A^{(0)}\right\rangle,\left\langle\frac{\delta F_{s}^{(0)}}{\delta F_{1}}, A^{(1)}\right\rangle\right\}+\cdots \\
\equiv \hat{D}^{(0)} F_{s}^{(0)}+\frac{1}{v}\left[\hat{D}^{(0)} F_{s}^{(1)}+\hat{D}^{(1)} F_{s}^{(0)}\right]+\cdots
\end{gathered}
$$

The above equation must be equivalent to inserting expansion 1.101 into the hierachy 1.98.

$$
\frac{\partial F_{s}}{\partial t}=\hat{L}_{s}\left(F_{s}^{(0)}+\frac{1}{v} F_{s}^{(1)}+\cdots\right)+\frac{1}{v} \sum_{i=1}^{s} \hat{I}_{i, s+1}\left(F_{s+1}^{(0)}+\frac{1}{v} F_{s+1}^{(1)}\right)
$$

We obtain a sequence of equations:

$$
\begin{gathered}
\hat{D}^{(0)} F_{s}^{(0)}=\hat{L}_{s} F_{s}^{(0)} \\
\ldots \\
\sum_{k=0}^{n} \hat{D}^{(k)} F_{s}^{(n-k)}=\hat{L}_{s} F_{s}^{(n)}+\sum_{i=1}^{s} \hat{I}_{i, s+1} F_{s+1}^{(n-1)}
\end{gathered}
$$

To obtain a kinetic equation from equation 1.102 to $O(1 / v)$, we have to close the equation with a expression for $F_{2}^{(0)}\left(F_{1}\right)$. Therefore, from the above, the below equation has to be solved with appropriate boundary conditions.

$$
\hat{D}^{(0)} F_{2}^{(0)}=\hat{L}_{s} F_{2}^{(0)}
$$

Here, we take this boundary condition to be when carried sufficiently far back in time, the particles become uncorrelated again. We found before that $\exp \left(t \hat{L}_{s}\right)$ propagates particles phase variables backward in time. We can label it to be $\hat{\Delta}_{-t}^{(s)}$ and employ it on the boundary condition.

$$
\lim _{t \rightarrow \infty} \hat{\Delta}_{-t}^{(s)} F_{s}\left(1 \ldots s ; F_{1}\right)=\lim _{t \rightarrow \infty} \hat{\Delta}_{-t}^{(s)} \prod_{k=1}^{s} F_{1}(k)
$$

Since this boundary condition holds for $F_{1}$, it must hold for $\lim _{t \rightarrow \infty} \hat{\Delta}_{-t}^{(1)} F_{1}(k)$ where $\hat{\Delta}_{-t}^{(1)}$ is the free 1 particle propagator :

$$
\begin{gathered}
\lim _{t \rightarrow \infty} \hat{\Delta}_{-t}^{(s)} F_{s}\left(1 \ldots s ; \hat{\Delta}_{-t}^{(1)} F_{1}\right)=\lim _{t \rightarrow \infty} \hat{\Delta}_{-t}^{(s)} \prod_{k=1}^{s} \hat{\Delta}_{-t}^{(1)} F_{1}(k) \\
=\lim _{t \rightarrow \infty} \hat{\Delta}_{-t}^{(s)} \prod_{k=1}^{s} F_{1}\left(\boldsymbol{x}_{k}+\frac{\boldsymbol{p}_{k}}{m} t, \boldsymbol{p}_{k}\right) \\
=\lim _{t \rightarrow \infty} \prod_{k=1}^{s} F_{1}\left(\hat{\Delta}_{-t}^{(s)} \boldsymbol{x}_{k}+\frac{\boldsymbol{p}_{k}}{m} t, \hat{\Delta}_{-t}^{(s)} \boldsymbol{p}_{k}\right) \\
==\lim _{t \rightarrow \infty} \prod_{k=1}^{s} F_{1}\left(\boldsymbol{x}_{k}^{(s)}, \boldsymbol{p}_{k}^{(s)}\right)
\end{gathered}
$$

---

<!-- Page 16 -->

where $\boldsymbol{x}_{k}^{(s)}, \boldsymbol{p}_{k}^{(s)}$ are the phase values for the $k$ th particles traveled back in time to $t=-\infty$. Returning to equation 1.104

$$
\hat{D}^{(0)} F_{s}^{(0)}=\left\langle\frac{\delta F_{s}^{(0)}}{\delta \hat{F}_{1}}, A^{(0)}\left(F_{1}\right)\right\rangle=\left\langle\frac{\delta F_{s}^{(0)}}{\delta \hat{F}_{1}}, \hat{L}_{1} F_{1}\right\rangle
$$

Again, since the above equation holds for $F_{1}$, it must hold for $\hat{\Delta}_{-1}^{(1)} F_{1}(k)$, we have

$$
\hat{D}^{(0)} F_{s}^{(0)}=\left\langle\frac{\delta F_{s}^{(0)}}{\delta \hat{\Delta}_{-1}^{(1)} F_{1}}, \hat{L}_{1} \hat{\Delta}_{-1}^{(1)} F_{1}\right\rangle
$$

After some computation and inserting the correct boundary condition, this give solution

$$
F_{s}^{(0)}\left(1 \ldots s ; F_{1}\right)=\prod_{k=1}^{s} F_{1}\left(\boldsymbol{x}_{k}^{(s)}, \boldsymbol{p}_{k}^{(s)}\right)
$$

Therefore, the lowest order distribution is given by a product of one particle distribution with phase variables propagated backward in time and coupled through the subspace $s$ - particle Hamiltonian. The final step in obtaining a closed kinetic equation is to obtain a consistent expression for $\hat{L}_{2}$. Starting with the definition introduced in equation 1.73

$$
\hat{L}_{2}=-\hat{\kappa}_{2}+\hat{O}_{12}
$$

where

$$
\hat{\kappa}_{2}=\frac{\boldsymbol{p}_{1}}{m} \cdot \frac{\partial}{\partial \boldsymbol{x}_{1}}+\frac{\boldsymbol{p}_{2}}{m} \cdot \frac{\partial}{\partial \boldsymbol{x}_{2}}
$$

Apply this operator on $F_{2}^{(0)}$ :

$$
\hat{O}_{12} F_{2}^{(0)}=\hat{\kappa}_{2} F_{2}^{(0)}+\hat{L}_{2} F_{2}^{(0)}
$$

Also recalling the definition of $A^{(0)}$ from 1.102 , we have

$$
A^{(0)}=-\frac{\boldsymbol{p}}{m} \cdot \frac{\partial}{\partial \boldsymbol{x}} F_{1}(\boldsymbol{x}, \boldsymbol{p})
$$

Therefore, we have
$\hat{L}_{2} F_{2}^{(0)}=-\frac{\boldsymbol{p}_{1}^{(2)}}{m} F_{1}\left(\boldsymbol{x}_{2}^{(2)}, \boldsymbol{p}_{2}^{(2)}\right) \cdot \frac{\partial}{\partial \boldsymbol{x}_{1}^{(2)}} F_{1}\left(\boldsymbol{x}_{1}^{(2)}, \boldsymbol{p}_{1}^{(2)}\right)-\frac{\boldsymbol{p}_{2}^{(2)}}{m} F_{1}\left(\boldsymbol{x}_{1}^{(2)}, \boldsymbol{p}_{1}^{(2)}\right) \cdot \frac{\partial}{\partial \boldsymbol{x}_{2}^{(2)}} F_{1}\left(\boldsymbol{x}_{2}^{(2)}, \boldsymbol{p}_{2}^{(2)}\right)$
we perform a coordinate transformation to find an alternative expression for $\hat{\kappa}_{2}$ :

$$
\begin{gathered}
\hat{\kappa}_{2}=\boldsymbol{g} \cdot \frac{\partial}{\partial \boldsymbol{r}}+\frac{\boldsymbol{p}_{1}}{m} \cdot \frac{\partial}{\partial \boldsymbol{x}_{1}} \equiv \hat{K}_{B}+\hat{K}_{1} \\
\boldsymbol{r}=\boldsymbol{x}_{2}-\boldsymbol{x}_{1} \\
\boldsymbol{g}=\frac{\boldsymbol{p}_{2}-\boldsymbol{p}_{1}}{m}
\end{gathered}
$$

Insert this back into $\mathrm{BY}_{1}$ for equation 1.102

$$
\frac{\partial F_{1}}{\partial t}+\frac{\boldsymbol{p}_{1}}{m} \cdot \frac{\partial F_{1}}{\partial \boldsymbol{x}_{1}}=\frac{1}{v} \int \mathrm{~d} 2 \hat{K}_{B} F_{2}^{(0)}+\frac{1}{v} \int \mathrm{~d} 2\left(\hat{L}_{2}-\hat{K}_{1}\right) F_{2}^{(0)}
$$

If $F_{2}^{(0)}$ is homogeneous (independent of $\boldsymbol{x}_{1}$ and $\boldsymbol{x}_{2}$ ) in the interaction domain, the second intergal evaluates to zero. The first integral becomes a collision integral when we insert the expression for $F_{2}^{(0)}$ :

$$
\int \mathrm{d} 2 \hat{K}_{B} F_{2}^{(0)}=\int \mathrm{d} \boldsymbol{x}_{2} \mathrm{~d} \boldsymbol{p}_{2} \boldsymbol{g} \cdot \frac{\partial}{\partial \boldsymbol{r}}\left\{F_{1}\left(\boldsymbol{p}_{1}^{2}\right) F_{1}\left(\boldsymbol{p}_{1}^{2}\right)\right\}
$$

The intergration range would correspond to the result before collision and after the collision. Label the result of collision as $\dot{\boldsymbol{p}}_{1}^{\prime}$ and $\dot{\boldsymbol{p}}_{2}^{\prime}$, we can evalue the integral for its spatial part within a cylinder or radius $b$ :

$$
\int \mathrm{d} \boldsymbol{x}_{2} \boldsymbol{g} \cdot \frac{\partial}{\partial \boldsymbol{r}}\left\{F_{1}\left(\boldsymbol{p}_{1}^{2}\right) F_{1}\left(\boldsymbol{p}_{1}^{2}\right)\right\}=2 \pi \int b \mathrm{~d} b g F_{1}\left(\boldsymbol{p}_{1}^{\prime}\right) F_{1}\left(\boldsymbol{p}_{1}^{\prime}\right)-F_{1}\left(\boldsymbol{p}_{1}\right) F_{1}\left(\boldsymbol{p}_{1}\right)
$$

---

<!-- Page 17 -->

we obtain the final equation

$$
\frac{\partial F_{1}}{\partial t}+\frac{\boldsymbol{p}_{1}}{m} \frac{\partial F_{1}}{\partial \boldsymbol{x}_{1}}=\frac{2 \pi}{v} \int \mathrm{~d} \boldsymbol{p}_{2} \int b \mathrm{~d} b g F_{1}\left(\boldsymbol{p}_{1}^{\prime}\right) F_{1}\left(\boldsymbol{p}_{1}^{\prime}\right)-F_{1}\left(\boldsymbol{p}_{1}\right) F_{1}\left(\boldsymbol{p}_{1}\right)
$$

This is the Boltzmann's Equation.

# 1.3 Boltzmann Equation and Fluid Dynamics 

### 1.3.1 Scattering

We know the Hamiltonian of two particle interaction in a central potential can be reduced to a onebody Hamiltonian with reduced mass with center-of-mass absorbed in the potential as a constant

$$
H=\frac{p^{2}}{2 \mu}+V(r)
$$

where $\mu=\frac{m_{1} m_{2}}{m_{1}+m_{2}}$ is the reduced mass, $\boldsymbol{p}=\mu \dot{\boldsymbol{r}}$ is the reduced momentum and $\boldsymbol{r}$ is the position of center of mass. In a central potential, the angular momentum is conserved and therefore particles stay in a plane of constant orientation. The degree of freedom is reduced from 3 to 2 . The Hamiltonian can be written in terms of polar coordiante variables $\phi$ and $r$

$$
H=\frac{p_{r}^{2}}{2 \mu}+\frac{L^{2}}{2 \mu r}+V(r)
$$

where $L$ is a constant total angular momentum. Solving the set of Hamiltonian Equations gives $\phi(r)$

$$
\phi(r)=\int \frac{L \mathrm{dr} / r^{2}}{\sqrt{2 \mu(E-V(r))-(L / r)^{2}}}+\text { constant }
$$

With the convention $V(\infty)=0$, we are addressing unbounded/scattering trajectories here. Now, consider the scattering event with two particles sufferring in a collision. We can quantify the range of interaction by $r_{0}$, the range of interaction, where $V(r)=0$ when $r>r_{0}$. Before and After the collision the interaction is minimal, there is no particle interactions, $V(r)=0$ and $r>r_{0}$. Assuming the interaction occurs in a interval around $t=0$. Let the relative velocity of two particles at $t= \pm \infty$ be $\boldsymbol{g}$ and $\boldsymbol{g}^{\prime}$. With conservation of energy: $\frac{1}{2} \mu g^{2}=\frac{1}{2} \mu g^{\prime 2}$ and conservation of momentum, we have

$$
\boldsymbol{g}^{\prime}=\boldsymbol{g}-2 \boldsymbol{\alpha}(\boldsymbol{\alpha} \cdot \boldsymbol{g})=(\widehat{I}-2 \widehat{\boldsymbol{\alpha} \boldsymbol{\alpha}}) \boldsymbol{g}=\widehat{S} \boldsymbol{g}
$$

where $\boldsymbol{\alpha}$ is the unit apsidal vector. The matrix $S$ is called the scattering matrix. The scattering matrix must be symmetric to be consistent with the time reversibility of physical laws $\left(\boldsymbol{g}=\widehat{S}^{-1} \boldsymbol{g}^{\prime}\right)$. It is clear that the magnitude of relative velocity of two particles is conserved in the collision. Also writing $L=\mu g s=\mu g^{\prime} s^{\prime}$ for total angular momentum and $s$ be an impact parameter, these are all constants and the integral for $\phi(r)$ can be consequently performed: Let angle $\boldsymbol{\psi}$ defined as

$$
\boldsymbol{\psi}=\phi\left(r_{\min }\right)-\phi(\infty)=\int_{r_{\min }}^{\infty} \frac{L \mathrm{dr} / r^{2}}{\sqrt{2 \mu(E-V(r))-(L / r)^{2}}}
$$

with the scattering angle $\theta$

$$
2 \boldsymbol{\psi}+\theta=\pi
$$

The total energy and total angular momentum is related by

$$
E=\frac{L^{2}}{2 \mu s}
$$

Let the form of potential be

$$
V(r)=K r^{-N}=K u^{N}
$$

the result of simplifying integration

$$
\boldsymbol{\psi}(s)=\int_{0}^{0} \frac{s \mathrm{du}}{\sqrt{1-s^{2} u^{2}-\left(K u^{N} / E\right)}}
$$

---

<!-- Page 18 -->

With the denominator vanishing at the apse:

$$
1-s^{2} \ddot{u}^{2}-\frac{V(\ddot{u})}{E}=0
$$

Defining the nondimensional inverse radius $\beta=s u$ and nondimensional impact parameter $b=$ $s(E / K)^{1 / N}$.

$$
\begin{gathered}
\psi(b)=\int_{0}^{\beta} \frac{\mathrm{d} \beta}{\sqrt{1-\beta^{2}-(\beta / b)^{N}}} \\
1-\bar{\beta}^{2}-\left(\frac{\bar{\beta}}{b}\right)^{N}=0
\end{gathered}
$$

Differential cross section can be defined by imagine a uniform beam of particles of Energy $E$ and intensity $I$ passing through a scatterer. The number of particles scattered into the element of solid angle $\mathrm{d} \Omega$ about $\Omega$ is proportional to the incident intensity $I$ and the element of solid angle $\mathrm{d} \Omega$. With a proportionality factor $\sigma$, the number of particles deflected into element $\mathrm{d} \Omega$ per second is

$$
I \sigma(\Omega) \mathrm{d} \Omega=I \mathrm{~d} \phi s \mathrm{~d} s
$$

This factor $\sigma$ is the differential cross section. is the The azimuthal angle $\phi$ which locates a section of the incident beam is the same angle in the spherical coordinate system around a scatterer, which gives $\mathrm{d} \Omega=\mathrm{d} \cos (\mathrm{d} \phi)$. Now inserting this into equation 1.125 , we have the classical equation for scattering cross section

$$
\sigma(E, \theta)=s(E, \theta) \frac{\operatorname{ds}(E, \theta)}{\sin \theta \mathrm{d} \theta}
$$

where the functional form of $s(E, \theta)$ can be obtained from the equation 1.120. On the other hand, the total cross section $\sigma_{T}$ relates to the total number of particles that are scattered out of the incident beam

$$
\sigma_{T}=\int_{4 \pi} \sigma \mathrm{~d} \Omega=\pi r_{0}^{2}
$$

It represents the obstructional area that scatterer presents to the incident beam. Another important factor occuring in the Boltzmann equation is $g \sigma \mathrm{~d} \cos (\theta)$ :

$$
g \sigma \mathrm{~d} \cos (\theta)=g s \mathrm{~d} s=\left(\frac{2 K}{\mu}\right)^{2 / N} g^{(N-4) / N} b \mathrm{~d} b
$$

The case when $N=4$ is the potential for Maxwell molecules. These particles have factor $g \sigma \mathrm{~d} \cos (\theta)$ independent of $g$. A fact that is closely related to the linearized Boltzmann Equation.

# 1.3.2 The Boltzmann Equation 

We now present an another way to derive the Boltzmann's Equation that is closer to Boltzmann's work and more physically revealing. Consider the limit of no interaction, particles are mutually independent and $F(\boldsymbol{x}, \boldsymbol{v}, t)$ statifies the one-particle Liouville Equation. This equation states that the net number of particles that enters the phase element $\delta \boldsymbol{x} \delta \boldsymbol{v}$ following a particle's trajectory in a time interval $\delta t$ is zero:

$$
\delta R=\delta \boldsymbol{x} \delta \boldsymbol{v} \delta t\left(\frac{\partial F}{\partial t}+\boldsymbol{v} \cdot \frac{\partial F}{\partial \boldsymbol{x}}+\frac{\boldsymbol{K}}{m} \cdot \frac{\partial F}{\partial \boldsymbol{v}}\right)=0
$$

where $\mathbf{K}$ is the externally supported force field. Now consider the interaction of particles with an interaction range of $r_{0}$. When $r>r_{0}$, interaction vanishes. When particles enter the interaction region $r \lesssim r_{0}$, they experience a collision. Let the mean distance between collisions, the mean free path, be $l$. The first important criterion for the below derivation is that $l \gg r_{0}$. In the present discussion, this phenomenon can be characterized by the net rate at which collision increase or decrease the number of particles in the phase volume. We may label $\delta R=\delta R_{+}-\delta R_{-}$, where $\delta R_{+}$ is the number of particles entering the phase volume and $\delta R_{-}$is the number of particles leaving the phase volume.

---

<!-- Page 19 -->

First, we consider $\delta R_{-}$. Velocity of all particles could be divided into two groups, the small velocity group that has velocity of interval $\delta \boldsymbol{v}$ about $\boldsymbol{v}$ and all particles with other velocities $\boldsymbol{v}_{1}$. The number of particles removed is simply the number of collisions that happen between particles of velocity about $\boldsymbol{v}$ and all other particles (velocity $\boldsymbol{v}_{1}$ ). Hence, in counting $\delta R_{-}$, all pairs between $\boldsymbol{v}$ and $\boldsymbol{v}_{1}$ must be considered. There are two properties of such pairs:

1. One particle is in phase element $\delta \boldsymbol{v} \delta \boldsymbol{x}$ about $(\boldsymbol{v}, \boldsymbol{x})$ and another particle is in $\delta \boldsymbol{v}_{1} \delta \boldsymbol{x}_{1}$ about $\left(\boldsymbol{v}_{1}, \boldsymbol{x}_{1}\right)$
2. They undergo a collision

By definition, the total number of particles that have such properties is (with phase variable $\boldsymbol{z}=(\boldsymbol{x}, \boldsymbol{v})$

$$
\delta R_{-}=\int_{1} F_{2}\left(\boldsymbol{z}, \boldsymbol{z}_{1}\right) \delta \boldsymbol{v}_{1} \delta \boldsymbol{x}_{1} \delta \boldsymbol{v} \delta \boldsymbol{x}
$$

To calculate $\delta R_{+}$, we notice the above intergal calculates how many particles were sent into $\delta \boldsymbol{v}$ about $\boldsymbol{v}$. Therefore, how many particles were sent out from $\delta \boldsymbol{v}$ about $\boldsymbol{v}$ could be given by an inverse collision

$$
\delta R_{+}=\int_{1} F_{2}\left(\boldsymbol{z}^{\prime}, \boldsymbol{z}_{1}^{\prime}\right) \delta \boldsymbol{v}_{1}^{\prime} \delta \boldsymbol{x}_{1}^{\prime} \delta \boldsymbol{v}^{\prime} \delta \boldsymbol{x}^{\prime}=\int_{1} F_{2}\left(\boldsymbol{z}^{\prime}, \boldsymbol{z}_{1}^{\prime}\right) \delta \boldsymbol{v}_{1} \delta \boldsymbol{x}_{1} \delta \boldsymbol{v} \delta \boldsymbol{x}
$$

where the primed variables denotes the variables after the collision. The second equality is true to preserve the time reversibility of collision. At last, we realize $\delta \boldsymbol{x}_{1}=g \delta t s \mathrm{ds} \mathrm{d} \phi=g \sigma \mathrm{~d} \Omega$. Inserting all of the above to equation 1.128 , we have

$$
\left(\frac{\partial F}{\partial t}+\boldsymbol{v} \cdot \frac{\partial F}{\partial \boldsymbol{x}}+\frac{\boldsymbol{K}}{m} \cdot \frac{\partial F}{\partial \boldsymbol{v}}\right)=\int_{1}\left[F_{2}\left(\boldsymbol{z}, \boldsymbol{z}_{1}\right)-F_{2}\left(\boldsymbol{z}^{\prime}, \boldsymbol{z}_{1}^{\prime}\right)\right] \partial \boldsymbol{v}_{1} g \sigma \mathrm{~d} \Omega
$$

If we can further assume $F_{2}$ is homogeneous and molecular chaos (reduced distributions are not correlated), we have $F_{2}\left(\boldsymbol{z}, \boldsymbol{z}_{1}\right)=F_{2}\left(\boldsymbol{v}, \boldsymbol{v}_{1}\right)=\frac{N-1}{N} F(\boldsymbol{v}) F\left(\boldsymbol{v}_{1}\right)$ and similar for $F_{2}\left(\boldsymbol{z}^{\prime}, \boldsymbol{z}_{1}^{\prime}\right)=\frac{N-1}{N} F\left(\boldsymbol{v}^{\prime}\right) F\left(\boldsymbol{v}_{1}^{\prime}\right)$. Now for a large ensemble of particles we neglect 1 compare to $N$. We have the Boltzmann's Equation again

$$
\left(\frac{\partial F}{\partial t}+\boldsymbol{v} \cdot \frac{\partial F}{\partial \boldsymbol{x}}+\frac{\boldsymbol{K}}{m} \cdot \frac{\partial F}{\partial \boldsymbol{v}}\right)=\iint\left[F^{\prime} F_{1}^{\prime}-F_{1} F\right] \mathrm{d} \boldsymbol{v}_{1} g \sigma \mathrm{~d} \Omega
$$

where we introduced expression $F_{1}=F\left(\boldsymbol{v}_{1}\right), F^{\prime}=F\left(\boldsymbol{v}^{\prime}\right), F_{1}^{\prime}=\left(\boldsymbol{v}_{1}^{\prime}\right)$. Note that the conservation derived from kinematics of two particle collision is still true:

$$
\begin{aligned}
\boldsymbol{v}_{1}^{\prime} & =\boldsymbol{v}_{1}+\boldsymbol{\alpha}(\boldsymbol{\alpha} \cdot \boldsymbol{g}) \\
\boldsymbol{v}^{\prime} & =\boldsymbol{v}-\boldsymbol{\alpha}(\boldsymbol{\alpha} \cdot \boldsymbol{g})
\end{aligned}
$$

A summary of all assumptions used to derive the Boltzmann Equation from the above reasoning

1. Range of interaction/mean free path $\ll 1$
2. Particles trajectories are rectilinear before and after collision.
3. $F$ is homogeneous in the entire course of collision
4. Molecular chaos

# 1.3.3 Fluid Dynamics and Boltzmann H-theorem 

We first wish to identifies the quantities (mass, momentum, energy, etc) that are perserved in the collision could still preserve itself under the effect of collision integral. The collision integral $\hat{J}(F)$ occurs in the left side of Boltzmann's equation.

$$
\hat{J}(F)=\iint\left[F^{\prime} F_{1}^{\prime}-F_{1} F\right]_{1} \mathrm{~d} \boldsymbol{v}_{1} g \sigma \mathrm{~d} \Omega
$$

Let the linear operator $\hat{I}[\phi(\boldsymbol{v})]$ be defined as

$$
\hat{I}[\phi(\boldsymbol{v})]=\iint \hat{J}(F) \phi(\boldsymbol{v}) \mathrm{d} \boldsymbol{v}=\iiint\left[F^{\prime} F_{1}^{\prime}-F_{1} F\right] \mathrm{d} \boldsymbol{v}_{1} g \sigma \mathrm{~d} \Omega \mathrm{~d} \boldsymbol{v} \phi(\boldsymbol{v})
$$

---

<!-- Page 20 -->

Changing variables $\left(\boldsymbol{v}, \boldsymbol{v}_{1}\right) \rightarrow\left(\boldsymbol{v}_{1}, \boldsymbol{v}\right)$,

$$
\hat{I}[\phi]=\hat{I}\left[\phi_{1}\right]
$$

Next, we consider

$$
\hat{I}\left[\phi^{\prime}\right]=\iiint\left[F^{\prime} F_{1}^{\prime}-F_{1} F\right] \mathrm{d} \boldsymbol{v}_{1} g \sigma \mathrm{~d} \Omega \mathrm{~d} \boldsymbol{v} \phi\left(\boldsymbol{v}^{\prime}\right)
$$

Another change of variables $\left(\boldsymbol{v}, \boldsymbol{v}_{1}\right) \rightarrow\left(\boldsymbol{v}^{\prime}, \boldsymbol{v}_{1}^{\prime}\right)$ would carry a Jacobian. However, Liouville Theorem states $\mathrm{d} \boldsymbol{v}_{1} \mathrm{~d} \boldsymbol{v}=\mathrm{d} \boldsymbol{v}_{1}^{\prime} \mathrm{d} \boldsymbol{v}^{\prime}$. Now we have

$$
\hat{I}\left[\phi^{\prime}\right]=\iiint\left[F_{1} F-F^{\prime} F_{1}^{\prime}\right] \mathrm{d} \boldsymbol{v}_{1}^{\prime} g^{\prime} \sigma^{\prime} \mathrm{d} \Omega^{\prime} \mathrm{d} \boldsymbol{v}^{\prime} \phi\left(\boldsymbol{v}^{\prime}\right)=-\hat{I}[\phi]
$$

Finally, a change of varible $\left(\boldsymbol{v}_{1}^{\prime}, \boldsymbol{v}^{\prime}\right) \rightarrow\left(\boldsymbol{v}^{\prime}, \boldsymbol{v}_{1}^{\prime}\right)$ gives

$$
\hat{I}\left[\phi^{\prime}\right]=\hat{I}\left[\phi_{1}^{\prime}\right]
$$

In summary, the following relation is true

$$
\hat{I}[\phi]=\frac{1}{4} \hat{I}\left[\phi+\phi_{1}-\phi^{\prime}-\phi_{1}^{\prime}\right]
$$

Now, if $\hat{I}[\phi]=0$, the quantity $\phi$ must preserve in a collision. This property is called collisional invariant. Three fundamental invariants are

$$
\hat{I}[1]=\hat{I}[\boldsymbol{v}]=\hat{I}\left[v^{2}\right]=0
$$

To obtain a series of equation that characterize macroscopic variables, we would like to integrate the Boltzmann Equation with respect to different moments. The right hand side after integration gives

$$
\int \mathrm{d} \boldsymbol{v}\left(\frac{\partial F}{\partial t}+\boldsymbol{v} \cdot \frac{\partial F}{\partial \boldsymbol{x}}\right)\left(\begin{array}{c}
1 \\
m \boldsymbol{v} \\
\frac{m v^{2}}{2}
\end{array}\right)=\left(\begin{array}{c}
\frac{\partial u}{\partial t}+\nabla \cdot n \boldsymbol{u} \\
\frac{\partial p \boldsymbol{u}}{\partial t}+\nabla \cdot \hat{P} \\
\frac{\partial e_{K}}{\partial t}+\nabla \cdot \boldsymbol{q}
\end{array}\right)
$$

and the lefthand side goes to zero from the property of collisional invariants. The above three equations are contiuity equation, conservation of momentum and conservation of energy. The macroscopic varibles are $n$ number density, $\boldsymbol{u}$ macroscopic mean velocity, $\hat{P}$ pressure tensor, $e_{K}$ kinetic energy density and $\boldsymbol{q}$ heat flow vector.

Moving to the subject of irreversibility: an irreversible law is characterized by an equation that has a solution peculiar to a particular direction of time. The Boltzmann's equation is such an irreversible equation. Consider the transformation $(\boldsymbol{x}, \boldsymbol{v}, t) \rightarrow(-\boldsymbol{x},-\boldsymbol{v},-t)$ : The left hand side of the equation switch sign, but the collision integral on the right hand side maintains its sign because the integration measure of $\hat{J}(f)$ is positive and $f \geq 0$. To formalize this irreversibility, we define Boltzmann's entropy

$$
\mathcal{H}=\int f \ln f \mathrm{dxdp}
$$

We multiply by $(1+\ln f)$ and integrate with respect to entire phase space volume for both side of Boltzmann's Equation

$$
\frac{\partial}{\partial t} \int \mathrm{~d} \boldsymbol{x} \mathrm{~d} \boldsymbol{v} f \ln f+\int \mathrm{d} \boldsymbol{x} \mathrm{~d} \boldsymbol{v} \frac{\partial}{\partial \boldsymbol{x}} \cdot \boldsymbol{v} \ln f+\frac{\boldsymbol{K}}{m} \cdot \frac{\partial}{\partial \boldsymbol{v}} \int \mathrm{~d} \boldsymbol{x} \mathrm{~d} \boldsymbol{v} f \ln f=\int \mathrm{d} \boldsymbol{x} \mathrm{~d} \boldsymbol{v} \hat{I}(1+\ln f)
$$

Dropping all the surface terms and passing to large volumes:

$$
\frac{1}{m} \frac{d}{d t} \mathcal{H}=\int \mathrm{d} \boldsymbol{x} \hat{I}(1+\ln f)
$$

With property of collisional invariant, we can write

$$
\begin{aligned}
4 \hat{I}(1+\ln f) & =\hat{I}(1+\ln f)+\hat{I}\left(1+\ln f_{1}\right)-\hat{I}\left(1+\ln f^{\prime}\right)-\hat{I}\left(1+\ln f_{1}^{\prime}\right) \\
& =\hat{I}\left(\ln \frac{f_{1} f}{f_{1}^{\prime} f^{\prime}}\right)=-\hat{I}\left(\ln \frac{f_{1}^{\prime} f^{\prime}}{f_{1} f}\right)
\end{aligned}
$$

---

<!-- Page 21 -->

Therefore we have the following evaluation for $4 \hat{I}(1+\ln f)$

$$
4 \hat{I}(1+\ln f)=-\int \mathrm{d} \boldsymbol{v} \int \mathrm{~d} \boldsymbol{v}_{1}\left[f^{\prime} f_{1}^{\prime}-f_{1} f\right] g \sigma \mathrm{~d} \Omega\left(\ln \frac{f_{1}^{\prime} f^{\prime}}{f_{1} f}\right)
$$

We notice the integrand of the above equation must be non-negative. Therefore $\hat{I}(1+\ln f) \leqslant 0$. Inserting this into equation 1.143, we have

$$
\frac{d}{d t} \mathcal{H} \leqslant 0
$$

This is Boltzmann's H-theorem. It states for arbitrary initial condition, $\mathcal{H}$ decreases until

$$
f^{\prime} f_{1}^{\prime}=f_{1} f
$$

and then it will remain constant. Using the method of Lagrange Multiplier, we can find the solution of this minimization problem:

$$
F_{0}(\boldsymbol{v})=n f_{0}(\boldsymbol{v})=\frac{n}{(2 \pi R T)^{3 / 2}} \exp \left(-\frac{(\boldsymbol{v}-\boldsymbol{u})^{2}}{2 R T}\right)
$$

where $R=k_{B} / m$ is the gas constant. This is the Maxwell-Boltzmann Distribution. Although the above distribution is a solution to the statstically balanced equation, there is another more general and physically relevant solution to the problem. It is obtained by observed that the above distribution always holds for particle distribution, temeperature and macroscopic velocity at a particular location and time $(\boldsymbol{x}, t)$.

$$
F_{0}(\boldsymbol{x}, \boldsymbol{v}, t)=\frac{n(\boldsymbol{x}, t)}{(2 \pi R T(\boldsymbol{x}, t))^{3 / 2}} \exp \left(-\frac{(\boldsymbol{v}-\boldsymbol{u}(\boldsymbol{x}, t))^{2}}{2 R T(\boldsymbol{x}, t)}\right)
$$

To distinguish the solutions, the latter is called local Maxwellian and the former is called absolute Maxwellian. Both of them could be termed as the equlibrium distribution. However, the thermal equilibrium for a force-free fluid implies all the macroscopic variables are constants. In this state, the system is better to be described with absolute Maxwellian. Prior to the thermal equilibrium stage, the system is in a state of local Maxwellian.

# 1.3.4 The Chapman-Enskog Expansion 

The collision integral can be written as

$$
\hat{J}(F \mid F)=-F(\boldsymbol{v}) \iint \mathrm{d} \boldsymbol{v}_{1} g \sigma \mathrm{~d} \Omega F\left(\boldsymbol{v}_{1}\right)+\iint \mathrm{d} \boldsymbol{v}_{1} g \sigma \mathrm{~d} \Omega F^{\prime} F_{1}^{\prime}
$$

we can identify that the quantity

$$
\nu(\boldsymbol{v})=\iint \mathrm{d} \boldsymbol{v}_{1} g \sigma \mathrm{~d} \Omega F\left(\boldsymbol{v}_{1}\right)
$$

represents a collsional frequency $\nu(\boldsymbol{v})$. Chapman Enskog expansion is relevant to the domain of large collisional frequency/small mean free path (with $C \sim l v$ being the thermal speed). The first step of Chapman Enskog Expansion is to write the Boltzmann's Equation in the following form

$$
\left(\frac{\partial}{\partial t}+\hat{\mathcal{D}}\right) F=\frac{1}{\varepsilon} \hat{I} F
$$

where $\varepsilon$ is a small dimensionless parameter. It will be eventually set to 1 . Here the operator $\hat{\mathcal{D}}$ is

$$
\hat{\mathcal{D}}=\boldsymbol{v} \cdot \frac{\partial F}{\partial \boldsymbol{x}}+\frac{\boldsymbol{K}}{m} \cdot \frac{\partial F}{\partial \boldsymbol{v}}
$$

The second step is to introduce the expansion

$$
F=F^{(0)}+\varepsilon F^{(1)}+\varepsilon^{2} F^{(2)}+\cdots
$$

---

<!-- Page 22 -->

where $F$ satisfy the moment relations

$$
\begin{aligned}
\int F \mathrm{~d} \boldsymbol{v} & =n \\
\int F \boldsymbol{v} \mathrm{~d} \boldsymbol{v} & =n \boldsymbol{u} \\
\int F \frac{m(\boldsymbol{v}-\boldsymbol{u})^{2}}{2} \mathrm{~d} \boldsymbol{v}=\int F \frac{m c^{2}}{2} \mathrm{~d} \boldsymbol{v} & =\frac{3}{2} k_{B} n T
\end{aligned}
$$

The third step is to stipulate that macroscopic variables $(n, \boldsymbol{u}, T)$ are all in the order of $O(1)$. For higher order terms $i>0$ in the series $F^{(i)}$, they only contribute to higher moments such as $\boldsymbol{Q}$ and $\tilde{P}$.

$$
\begin{aligned}
& \int F^{(0)}\left(\begin{array}{c}
1 \\
\boldsymbol{v} \\
c^{2}
\end{array}\right) \mathrm{d} \boldsymbol{v}=\left(\begin{array}{c}
n \\
n \boldsymbol{u} \\
\frac{3}{2} k_{B} n T / m
\end{array}\right) \\
& \int F^{(i)}\left(\begin{array}{c}
1 \\
\boldsymbol{v} \\
c^{2}
\end{array}\right) \mathrm{d} \boldsymbol{v}=\left(\begin{array}{l}
0 \\
0 \\
0
\end{array}\right) \\
& \boldsymbol{Q}=\sum_{l} \varepsilon^{l} \boldsymbol{Q}^{(l)}=\frac{1}{2} \sum_{l} \varepsilon^{l} \int F^{(l)} \boldsymbol{c} m c^{2} \mathrm{~d} \boldsymbol{c} \\
& \tilde{P}=\sum_{l} \varepsilon^{l} \tilde{P}^{(l)}=\sum_{l} \varepsilon^{l} \int F^{(l)} m \widetilde{\boldsymbol{c}} \mathrm{~d} \boldsymbol{c}
\end{aligned}
$$

Now we can expand $\hat{\mathcal{D}}$ and $\hat{J}$ in the powers of $\varepsilon$ :

$$
\hat{\mathcal{D}} F=\hat{\mathcal{D}} F^{(0)}+\varepsilon \hat{\mathcal{D}} F^{(1)}+\cdots
$$

Define the ordered operator

$$
\hat{J}^{(s)}\left(F^{(0)}, F^{(1)}, \ldots, F^{(s)}\right)=\sum_{n}^{n+l=s} \sum_{l} \hat{J}\left(F^{(l)} \mid F^{(n)}\right)
$$

we have the expansion for collision integral $\hat{J}$

$$
\hat{J}(F \mid F)=\hat{J}^{(0)}\left(F^{(0)}\right)+\hat{J}^{(1)}\left(F^{(0)}, F^{(1)}\right)+\hat{J}^{(2)}\left(F^{(0)}, F^{(1)}, F^{(2)}\right)+\cdots
$$

The fourth step is to expand the time derivative. Chapman-Enskog expansion stipulate that the time dependence of $F$ is solely dependent on the hydrodynamic variables $n, \boldsymbol{u}, T$. The time derivative is expanded as

$$
\frac{\partial}{\partial t}=\frac{\partial_{0}}{\partial t}+\varepsilon \frac{\partial_{1}}{\partial t}+\varepsilon^{2} \frac{\partial_{2}}{\partial t}+\cdots
$$

the physical meaning of different order of time derivative is that lowest order term vary most rapidly and the higher order terms are more slowly varying. The significance of different order of time derivative can be found by integrating different orders of $\varepsilon$ equations to obtain corresponding order of conservation equations. $(r>0)$

$$
\begin{aligned}
\frac{\partial_{0} n}{\partial t} & =-\nabla \cdot n \boldsymbol{u} \\
\frac{\partial_{r} n}{\partial t} & =0 \\
\frac{\partial_{0} \boldsymbol{u}}{\partial t} & =-(\boldsymbol{u} \cdot \nabla) \boldsymbol{u}-\frac{1}{\rho} \nabla \cdot \tilde{P}^{(0)}+\tilde{\boldsymbol{K}} \\
\frac{\partial_{\rho} \boldsymbol{u}}{\partial t} & =-\frac{1}{\rho} \nabla \cdot \tilde{P}^{(r)} \\
\frac{\partial_{0} T}{\partial t} & =-(\boldsymbol{u} \cdot \nabla) T-\frac{2}{3} \frac{1}{n k_{B}}\left[\tilde{P}^{(0)}: \widetilde{\nabla \boldsymbol{u}}+\nabla \cdot \boldsymbol{Q}^{(0)}\right] \\
\frac{\partial_{r} T}{\partial t} & =-\frac{2}{3} \frac{1}{n k_{B}}\left[\tilde{P}^{(r)}: \widetilde{\nabla \boldsymbol{u}}+\nabla \cdot \boldsymbol{Q}^{(r)}\right]
\end{aligned}
$$

---

<!-- Page 23 -->

Now, substituting all above expansions and matching terms with different orders of $\varepsilon$, we have

$$
\begin{aligned}
0 & =\hat{J}^{(0)}\left(F^{(0)}\right) \\
\left(\frac{\partial_{0}}{\partial t}+\hat{\mathcal{D}}\right) F^{(0)} & =\hat{J}^{(1)}\left(F^{(0)}, F^{(1)}\right) \\
\left(\frac{\partial_{0}}{\partial t}+\hat{\mathcal{D}}\right) F^{(0)}+\frac{\partial_{1}}{\partial t} F^{(0)} & =\hat{J}^{(2)}\left(F^{(0)}, F^{(1)}, F^{(2)}\right)
\end{aligned}
$$

The first order solution to the above equation is simply solving $\hat{J}^{(0)}\left(F^{(0)} \mid F^{(0)}\right)=0$. For a state that is no longer evolved by the collision integral, the distribution must be in equilibrium. Therefore, for the first order solution, the result is calculated before to be local Maxwellian (equation 1.147). We can calculate the moments of local Maxwellian to find the macroscopic quantities: $(\boldsymbol{c}=\boldsymbol{v}-\boldsymbol{u})$

$$
\begin{aligned}
n(\boldsymbol{x}, t) & =\int F^{0}(\boldsymbol{x}, \boldsymbol{v}, t) \mathrm{d} \boldsymbol{v} \\
n(\boldsymbol{x}, t) \boldsymbol{u}(\boldsymbol{x}, t) & =\int F^{0}(\boldsymbol{x}, \boldsymbol{v}, t) \boldsymbol{v} \mathrm{d} \boldsymbol{v} \\
\frac{3}{2} n(\boldsymbol{x}, t) k_{B} T(\boldsymbol{x}, t) & =\int F^{0}(\boldsymbol{x}, \boldsymbol{c}, t) \frac{m c^{2}}{2} \mathrm{~d} \boldsymbol{c} \\
\boldsymbol{Q} & =\int F^{0}(\boldsymbol{x}, \boldsymbol{c}, t) \frac{m c^{2}}{2} \boldsymbol{c} \mathrm{~d} \boldsymbol{c}=0 \\
\tilde{P} & =\int F^{0}(\boldsymbol{x}, \boldsymbol{c}, t) m \overleftrightarrow{\boldsymbol{c}} \mathrm{d} \boldsymbol{c}=\tilde{I} n(\boldsymbol{x}, t) k_{B} T(\boldsymbol{x}, t)
\end{aligned}
$$

Therefore, to the lowest order, the heat flow is none and the pressure tensor is purely diagonal. Substituting these into the conservation laws, they give us the Euler's equations.

$$
\begin{aligned}
\frac{\partial n}{\partial t}+\nabla \cdot n \boldsymbol{u} & =0 \\
\rho\left(\frac{\partial}{\partial t}+\boldsymbol{u} \cdot \nabla\right) \boldsymbol{u}+\nabla p & =\rho \boldsymbol{K} \\
\left(\frac{\partial}{\partial t}+\boldsymbol{u} \cdot \nabla\right) \frac{p}{n^{5 / 3}} & =0
\end{aligned}
$$

For the second order solution, we introduce the function $F^{(1)}=F^{(0)} \Phi$. Using equation 1.155, we have

$$
\hat{\square} \Phi=\frac{1}{F^{(0)}} \hat{J}^{(1)}\left(F^{(0)}, F^{(0)} \Phi\right)=\frac{1}{F^{(0)}}\left[\hat{J}\left(F^{(0)} \mid F^{(0)} \Phi\right)+\hat{J}\left(F^{(0)} \Phi \mid F^{(0)}\right)\right]
$$

The second order equation can be written as

$$
\frac{1}{F^{(0)}}\left(\frac{\partial_{0}}{\partial t}+\hat{\mathcal{D}}\right) F^{(0)}=\hat{\square} \Phi
$$

Terms on the left gives

$$
\begin{aligned}
\frac{1}{F^{(0)}} \frac{\partial_{0} F^{(0)}}{\partial t} & =\left[\frac{1}{n} \frac{\partial_{0} n}{\partial t}+2 \boldsymbol{\xi} \cdot \frac{\partial_{0} \boldsymbol{u}}{\partial t}+\left(\xi^{2}-\frac{3}{2}\right) \frac{1}{T} \frac{\partial_{0} T}{\partial t}\right] \\
\frac{1}{F^{(0)}} \boldsymbol{v} \cdot \frac{\partial_{0} F^{(0)}}{\partial \boldsymbol{x}} & =\boldsymbol{v} \cdot\left[\frac{1}{n} \frac{\partial_{0} n}{\partial \boldsymbol{x}}+2 \boldsymbol{\xi} \cdot \frac{\partial_{0} \boldsymbol{u}}{\partial \boldsymbol{x}}+\left(\xi^{2}-\frac{3}{2}\right) \frac{1}{T} \frac{\partial_{0} T}{\partial \boldsymbol{x}}\right]
\end{aligned}
$$

where

$$
\xi^{2}=\frac{c^{2}}{2 R T}
$$

Replacing these time derivatives to expressions on the RHS of equations 1.158

$$
\sqrt{2 R T}\left(\xi^{2}-\frac{5}{2}\right) \boldsymbol{\xi} \cdot \nabla \ln T+2\left(\overline{\boldsymbol{\xi} \boldsymbol{\xi}}-\frac{1}{3} \xi^{2} \tilde{I}\right) \cdot \widehat{\nabla \boldsymbol{u}}=\hat{\square} \Phi
$$

---

<!-- Page 24 -->

this is a linear inhomogeneous integral equation for distribution $\Phi$. If this is solved, the distribution $F$ is known to the second order with $F=F^{(0)}(1+\Phi)$. From the structure of operator $\hat{\square}$, the homogeneous solution should be any linear combination with three summational invariants:

$$
\Phi_{h}=\alpha+\boldsymbol{\beta} \cdot m \boldsymbol{c}+\frac{1}{2} \gamma m c^{2}
$$

On the other hand, noticing the form of the inhomogeneous part of the equation, the particular solution of the equation may take form of

$$
\Phi_{i}=\boldsymbol{A}(\xi) \cdot(2 R T)^{1 / 2} \nabla \ln T+2 \overline{\widehat{B(\xi)}}: \widehat{\nabla \boldsymbol{u}}
$$

Inserting this back into the equation 1.164 gives two separate integral equation for vector field $\boldsymbol{A}(\xi)$ and tensor field $\overline{\widehat{B(\xi)}}$.

$$
\begin{aligned}
& \hat{\square} \boldsymbol{A}(\xi)=\left(\xi^{2}-\frac{5}{2}\right) \boldsymbol{\xi} \\
& \widehat{\square} \overline{\widehat{B(\xi)}}=\left(\overline{\boldsymbol{\xi} \boldsymbol{\xi}}-\frac{1}{3} \xi^{2} \widehat{I}\right)
\end{aligned}
$$

Note that the only variables in $\boldsymbol{A}(\xi)$ and $\overline{\widehat{B(\xi)}}$ are $\boldsymbol{\xi}, n$ and $\boldsymbol{u}$. The only vector/tensor field that can form these elements is $\boldsymbol{\xi}$ and $\left(\overline{\boldsymbol{\xi} \boldsymbol{\xi}}-\frac{1}{3} \xi^{2} \widehat{I}\right)$ itself. Hence, we have the form of the vector/tensor field as

$$
\begin{aligned}
\boldsymbol{A}(\xi) & =\mathcal{A}\left(\xi^{2}\right) \boldsymbol{\xi} \\
\overline{\widehat{B(\xi)}} & =\mathcal{B}\left(\xi^{2}\right)\left(\overline{\boldsymbol{\xi} \boldsymbol{\xi}}-\frac{1}{3} \xi^{2} \widehat{I}\right)
\end{aligned}
$$

Now, we can form the total solution with scalar function $\mathcal{A}\left(\xi^{2}\right)$ and $\mathcal{B}\left(\xi^{2}\right)$. There are several constriants (the higher order equation does not contribute to the macroscopic quantities $n, \boldsymbol{u}, T$ in equation 1.153) that could help determine the coefficients in the homogeneous solution.

$$
\begin{array}{r}
\int F^{0}\left(\alpha+\frac{1}{2} \gamma m c^{2}\right) \mathrm{d} \boldsymbol{c}=0 \\
\int F^{0}\left(\mathcal{A}\left(\xi^{2}\right) \nabla \ln T+\boldsymbol{\beta} m\right) m c^{2} \mathrm{~d} \boldsymbol{c}=0 \\
\int F^{0}\left(\alpha+\frac{1}{2} \gamma m c^{2}\right) \frac{m c^{2}}{2} \mathrm{~d} \boldsymbol{c}=0
\end{array}
$$

The first and the third constriants implies that $\alpha=\gamma=0$. The second constriant implies that $\boldsymbol{\beta}$ is in the same direction as $\nabla \ln T$ and can be absorbed into the latter. Hence, the entire homogeneous part is dropped. We have the full solution of Boltzmann Equation to the second order.

$$
F=F^{(0)}\left[1+(2 R T)^{1 / 2} \mathcal{A}\left(\xi^{2}\right) \boldsymbol{\xi} \cdot \nabla \ln T+2 \mathcal{B}\left(\xi^{2}\right)\left(\overline{\boldsymbol{\xi} \boldsymbol{\xi}}-\frac{1}{3} \xi^{2} \widehat{I}\right): \widehat{\nabla \boldsymbol{u}}\right]
$$

Each iteration of Chapman Enskog Expansion gives a more detailed set of hydrodynamical equations, better suited to higher order fluctuations in the fluid. The first iteration as noted above gives the Euler's Equations. The second iteration replace the momentum equation by the new Navier-Stokes Equation. Using equation 1.153, we have an expression for stress tensor $\tilde{P}$ in the conservation equations 1.141

$$
\tilde{P}=2 R T \int F \overline{\tilde{\boldsymbol{\xi}}} \mathrm{d} \boldsymbol{c}=\tilde{I} p+\tilde{P}^{(1)}=\tilde{I} p+m \int F^{(0)} \Phi \overline{\boldsymbol{c}} \mathrm{d} \boldsymbol{c}
$$

The first term in the integral contribute to the diagonal pressure tensor. The second term in equation integrates to zero. The remaning third term contribute to the integral and were explicitly written and labelled as $\tilde{P}^{(1)}$. After some tensor integrals, we have expression

$$
\tilde{P}^{(1)}=\frac{4 k_{B} T}{5}\left(\tilde{\Lambda}-\frac{1}{3} \tilde{I} \operatorname{Tr} \widehat{\nabla \boldsymbol{u}}\right) \int F^{0} \overline{\widehat{B(\xi)}}: \widehat{\square} \overline{\widehat{B(\xi)}} \mathrm{d} \boldsymbol{c}
$$

---

<!-- Page 25 -->

where $\tilde{\Lambda}$ is the symmetric strain tensor in fluid dynamics.

$$
\Lambda_{i k}=\frac{1}{2}\left(\frac{\partial u_{i}}{\partial x_{k}}+\frac{\partial u_{k}}{\partial x_{i}}\right)
$$

We recognize the coefficient of shear viscosity from the expression

$$
\eta=-\frac{2}{5} k_{B} T \int F^{0} \overline{\bar{B}(\xi)}: \bar{\square} \overline{\bar{B}(\xi)} \mathrm{d} \boldsymbol{c}
$$

With this result we can write the momentum equation

$$
\rho\left(\frac{\partial \boldsymbol{u}}{\partial t}+(\boldsymbol{u} \cdot \nabla) \boldsymbol{u}\right)=-\nabla p+\eta \nabla^{2} \boldsymbol{u}-\frac{\eta}{3} \nabla(\nabla \cdot \boldsymbol{u})
$$

With incompressiblility condition $(\nabla \cdot \boldsymbol{u})=0$, the last term goes to zero. We demonstrated that the second order solution of Chapman Enskog Expansion can be used to reconstruct Navier-Stokes equation.

---

<!-- Page 26 -->

# Chapter 2 

## Lattice Boltzmann Method

Since the inventio of lattice boltzmann method, a range of literature has been dedicated to this field. Textbooks are also written to explain the numerical scheme and its related applications. Here we primarily follow two textbooks: [12] S. Succi. The lattice Boltzmann equation for fluid dynamics and beyond. and [8] Timm Krüger, Halim Kusumaatmaja, Alexandr Kuzmin, Orest Shardt, Goncalo Silva, and Erlend Magnus Viggen. The Lattice Boltzmann Method: Principles and Practice.

### 2.1 Lattice Gas Celluar Automata Method

The predecessor of LBM is developed by Uriel Frish, Brosl Hasslacher and Yves Pomeau and termed as Lattice Gas Celluar Automata Method (LGCA)[5]. The discovery here is that a simple cellular automaton(FHP automaton) obeying only conservation laws at a microscopic level was able to reproduce the complexity of a real fluid flow. We will start by consider such automaton as a regular lattice with hexagonal symmetry with each site surrounded by six neighbors identified by six connecting vectors $\vec{c}_{i}=c_{i a}, i=1, \ldots, 6$ where the index $a$ scanning all spatial dimensions.

Each lattice site hosts up to six particles with following prescriptions

1. All particle have the same mass $m=1$.
2. Particles can move only along one of the six directions defined by $c_{i a}$.
3. In one time-cycle the particles hop to the nearest neighbor pointed by the vector $c_{i a}$. Both larger jumps and smaller jumps are forbidden. Hence, all of the particles would have the same energy.
4. Exclusion principle: no particles sitting on the same site can move along the same direction $c_{i a}$.
In a real gas, classical molecules can move along any direction and any speed without limitation of exclusion principle, whereas here they are confined in a six-barred cage, allowed to move along only six monochromatic beams and subject to exclusion principle. However, as we will see, this crude image can simulate realistic hydrodynamics.

With the above prescriptions, the state of each site could be specified by a yes or no saying whether there is a particle on the site. This dichotomic situation could be coded with a single binary digit per site and direction so that the entire state of the lattice gas could be specified by $6 N$ bits, whereas $N$ is the number of lattice sites. Introducing an occupation number $n_{i}$

$$
\begin{array}{ll}
n_{i}(\vec{x}, t)=0 & \text { particle absence at site } \vec{x} \text { and time } t \\
n_{i}(\vec{x}, t)=1 & \text { particle presence at site } \vec{x} \text { and time } t
\end{array}
$$

The collection of occpuation numbers $n_{i}(\vec{x}, t)$ over the entire lattice defines a $6 N$-dimensional timedependent Boolean field whose evolution takes place in a Boolean phase-space consisiting of $2^{6 N}$ discrete states. The microdynamics of the Boolean field cannot be expected to reproduce the true molecular dynmaics to any reasonable degree of microscopic accuracy. However, as has been known since Gibbs, many different microscopic systems can give rise to the same macroscopic dynamics, an it could be hoped that the macroscopic dynamics of the lattice Boolean field can reproduce the real-life hydordynamics even if its microdynamics does not.

The evolution rules of automatons have two mechanisms based on hydrodynamics: Free-steaming and Collisions. Free-streaming consists of particle transfer from site to site according to the set of discrete speeds $c_{i a}$. A particle at site $\vec{x}$ at time $t$ with speed $c_{i a}$ will move to site $\vec{x}+c_{i a}$ at time $t+1$.

$$
n_{i}\left(\vec{x}+c_{i a}, t+1\right)=n_{i}(\vec{x}, t)
$$

---

<!-- Page 27 -->

This defines the discrete free-streaming operator as

$$
\Delta_{i} n_{i}=n_{i}\left(\vec{x}+c_{i a}, t+1\right)-n_{i}(\vec{x}, t)
$$

This is a direct translation of the Boltzmann free-streaming opeator $D_{t}=\partial_{t}+v_{a} \partial_{a}$ to a discrete lattice in which space-time are discretized according to the synchronous "light-cone" rule:

$$
\Delta x_{i a}=c_{i a} \Delta t
$$

The relation between discrete and continuous streaming operator is

$$
\Delta_{i}=e^{D_{i}}-1
$$

where $D_{i}=\partial_{t}+c_{i a} \partial_{a}$ is the generator of space-time translations along the $i-$ th direction. The magnitude $c$ plays the role of "light speed" in the discrete world where no signal can propagte faster than $c$ in lattice.

Once on the same site, particles interact and reshuffle their momenta as long as it is in a direction allowed by the lattice. This mimics the real-life collisions of real gas. Any such collisions must 1 , conserve the particle number, and 2 , conserve the total momentum. Once a sufficiently large group of particles and collisions are considered, these properties will be necessary to achieve hydrodynamic behavior. However, because Navier Stokes equation has rotational invariance, the hexagonal structure becomes necessary when any other simpler lattice shape(square, for example) does not achieve this rotational invariance. This requirement is relevant to the number of spacetime resolutions. In general, for a $D$ dimensional lattice with $b$ discrete speeds, the following tensorial identity must be ensured

$$
\left[\sum_{i=1}^{b} c_{i a} c_{i b}\left(c_{i c} c_{i d}-\frac{c^{2}}{D} \delta_{c d}\right)\right] u_{c} u_{d}=u_{a} u_{b}
$$

for any choice of dyadic $u_{a} u_{b}$ to provide enough symmetry. Intuitively, the point is that as opposed to a scalar field, a fourth-order tensor field are demanding when they sense more details of the space-time fabric. It is harder to make them "believe" that a square and a circle are the same thing. Surprisingly and pleasingly, a hexagon does serve this purpose.

Back to the discussion of coliision. Symbollically, its effect on the occupation number is a change from $n_{i}$ to $n_{i}^{\prime}$ on the same site

$$
n_{i}^{\prime}-n_{i}=C_{i}(\bar{n})
$$

where $\bar{n}=\left[n_{1}, n_{2}, \ldots n_{b}\right\rangle$ denotes the set of occupation numbers at a given lattice site. To formalize the expression of $C_{i}$, we can label the phase-space by a bit of string $\bar{s}=\left[s_{1}, s_{2}, \ldots s_{b}\right\rangle$ spanning the set of all possible $2^{b}$ states at a given lattice site. For example, numbering discrete speeds 1-6 counterclockwise starting from the rightward propagation, $c_{1 s}=1, c_{1 y}=0$, the pre and post collisional states read $\bar{s}=|100100\rangle$ and $\vec{s}=|010010\rangle$ respectively.

It is natural to define a transitional matrix $A\left(\bar{s}, \bar{s}^{\prime}\right)$ flagging all permissible collisions from source state $\bar{s}$ to destination state $\bar{s}^{\prime}$ as follows

$$
\begin{array}{ll}
A_{i}\left(\bar{s}, \bar{s}^{\prime}\right)=0 & \text { collisions forbidden } \\
A_{i}\left(\bar{s}, \bar{s}^{\prime}\right)=1 & \text { collisions allowed }
\end{array}
$$

Collisions are allowed if they obey conservation laws. This transitional matrix obeys the semidetailed balance condition:

$$
\sum_{\bar{s}} A\left(\bar{s}, \bar{s}^{\prime}\right)=1
$$

It means that every destination state necessarily comes from a source state within the phase-space of the automaton. This condition does not imply a one-to-one source-destination relationship. Define the probability to have $\bar{s}$ as input state with occpuation

$$
P(\bar{s}, \bar{n})=\prod_{i=1}^{b} n_{i}^{s_{i}} \hat{n}_{i}^{\bar{s}_{i}}
$$

---

<!-- Page 28 -->

where the accent $\tilde{n}$ denotes complement to one $\left(\tilde{n}_{i}=1-n_{i}\right)$ as befits to fermionic degrees of freedom. Example: the probability of occupying state $\bar{s}=[100100\rangle$ as an input string is given by $P[100100\rangle=n_{1} \tilde{n_{2}} \tilde{n_{3}} n_{4} \tilde{n_{5}} \tilde{n_{6}}$. This probability is always zero except when a particle with speed $\vec{c}_{1}$ AND a particle with speed $\vec{c}_{4}$ are sitting simulatenously on the node. Here, the "particle absence" is equivalent with "hole presence", echoing the particle-hole symmetry with fermionic matter. With these preparations, the collision operator could be recast in the traditional gain minus loss form

$$
C_{i}=\sum_{\bar{s}, \bar{s}^{\prime}}\left(s_{i}^{\prime}-s_{i}\right) P(\bar{s}, \bar{n}) A\left(\bar{s}, \bar{s}^{\prime}\right)
$$

Here we can check that the value of $C_{i}$ is either -1 (annihilation), 0 (no action) and +1 (generation). With this definition, the collision operator can obey the property of preserving the Boolean nature of the occupation numbers. This can be verified by the Boolean-breaking occurrences, $n_{i}=0$, $C_{i}=-1$ or $n_{i}=1, C_{i}=+1$ will never occur. To sum up, the LGCA update rule reads as follows

$$
\Delta_{i} n_{i}=C_{i}
$$

This represents the microdynamic equation for the Boolean lattice gas, the analogue of Newton equations for real molecules. This equation consitutes the starting point of a lattice BBGKY hierachy, ending up with the Navier-Stokes Equation. At each level, one formulates a lattice counterpart of the various approximations pertaining to the four levels of the hierachy. The fundations of these approximtion procedure are the conservation laws of classical mechanics

1. Mass and Momentum Conservation

$$
\begin{gathered}
\sum_{i} C_{i}=0 \\
\sum_{i} c_{i a} C_{i}=0
\end{gathered}
$$

2. Angular Momentum Conservation (rotational invariance)

$$
\sum_{i} c_{i a} c_{i b} c_{i c} c_{i d}=\frac{b c^{4}}{D(D+2)}\left(\delta_{a c} \delta_{b d}+\delta_{a d} \delta_{b c}-\frac{2}{D} \delta_{a b} \delta_{c d}\right)
$$

The standard bottom-up procedure will take us from many-body particle dynamics to continumm fluid-like equations. This is done by three formal steps familar from classical statistical mechanics.

1. From Newton-Hamilton to Liouville

With ergodicity assumption, the deterministic Newton-Hamilton equations of motion are replaced by probabilistic Liouville Equation for $N$-body distribution function $f_{N}\left(\overrightarrow{x_{1}}, \overrightarrow{v_{1}}, \ldots\right.$, $\overrightarrow{x_{N}}, \overrightarrow{v_{N}} ; t)$. Symbollically:

$$
\ddot{z}(t)=\left[\overrightarrow{x_{1}}, \overrightarrow{v_{1}}, \ldots, \overrightarrow{x_{N}}, \overrightarrow{v_{N}}\right] \rightarrow f_{N}\left(\overrightarrow{x_{1}}, \overrightarrow{v_{1}}, \ldots, \overrightarrow{x_{N}}, \overrightarrow{v_{N}} ; t\right)
$$

2. From Liouville to Boltzmann

With diluteness assumption, high-order distribution functions are expressed in terms of the low-order ones by integration over many-body phase-space coordinates. The procedure ends up at the lowest level (Boltzmann) with single-body distribution function. Hierachy is closed by a truncation in the diluteness parameter $(s / d)^{3}$. ( $s$ is the effective molecular diamater and $d$ is the mean intermolecular distance. Symbollically:

$$
f_{N}\left(\overrightarrow{x_{1}}, \overrightarrow{v_{1}}, \ldots, \overrightarrow{x_{N}},overrightarrow{v_{N}} ; t\right) \rightarrow f\left(\overrightarrow{x_{1}}, \overrightarrow{v_{1}}, t\right)
$$

3. From Boltzmann to Navier-Stokes

With integration over momentum space degrees of freedom, a set of partial differential equations for the space-time depedent moments of the Boltzmann distribution function is obtained. For $m=1, \ldots M$, number of moments, the set of moment equations takes the form of an $M$-dimensional generalized continuity equation

$$
\partial_{t} \rho_{m}+\partial_{a} J_{a m}=C_{m}
$$

---

<!-- Page 29 -->

where

$$
\begin{gathered}
\rho_{m}=\int f \phi_{m}(\vec{v}) \mathrm{d} \vec{v} \\
J_{a m}=\int f v_{a} \phi_{m}(\vec{v}) \mathrm{d} \vec{v} \\
C_{m}=\int C[f] \phi_{m}(\vec{v}) \mathrm{d} \vec{v}
\end{gathered}
$$

with $\rho_{m}$, the generalized densities, $J_{a m}$, the generalized currents and $C_{m}$, the effect of intermolecular interactions. Here, $\left\{\phi_{m}\right\}$ is a complete set of basis functions in the momentum space, typically Hermite polynomials. Equation 2.18 does not close the hierachy. Instead, it is closed by splitting up the one particle distribution into a local equilibrium and nonequilibrium components:

$$
f=f^{e}+f^{\mathrm{ne}}
$$

with the assumption that $f^{\mathrm{ne}} \sim O(\mathrm{Kn}) f^{e}$, Kn being the Knudsen number. To leading order $O(\mathrm{Kn})$, this delivers the Euler Equations for inviscid flows. On the next order $O\left(\mathrm{Kn}^{2}\right)$, the result is desired Navier-Stokes Equation.
The same steps are involved in the process of deriving the lattice Navier-Stokes equations from lattice BBGKY, which is to be handled with great care since contiuum symmetries are always at risk of being broken by the lattice discreteness. The most dangerous effects bear upon the discrete lattice equilibria. Since LGCA obeys exclusion principle, lattice quilibria are expressed by a Fermi - Dirac distribution

$$
f_{i}^{\mathrm{eq}}=\frac{\rho / b}{1+e^{\Phi_{i}}}
$$

where $\Phi_{i}$ is a linear combination of the collisional invariants, mass and momentums. Restricting to isothermal ideal fluids only, where energy and mass are the same, we have

$$
\Phi_{i}=A+B c_{i a} u_{a}
$$

Using the method of Lagrangian Multipliers, we would like to obtain a closed, analytical expression of the parameters $A$ and $B$ by minmizing the Boltzmann Entropy.

$$
H=-\int \tilde{f} \ln \tilde{f} \mathrm{~d} \vec{v}
$$

where $\tilde{f}=f /(1-f)$ for fermions. However, the form of Fermi-Dirac distribution leads to two nonlinear equations when applying the conservation relations in the discrete case, whereas in a continuous momentum space, this would not be a problem. In the discrete momentum space, the exponential function plays no special role, we would like to replace it by other functional form like polynomials. It turns out to be natural to expand the expoential in powers of the flow field (by Mach number $\mathrm{Ma}=u / c_{s}$ ) Conservation constaints can be then imposed order-by-order in the perturbation series and since Navier Stokes Equation exhibits a quadratic nonlinearity, we should expect a second-order expansion should occur. The specific LGCA equilibria expression is

$$
f_{i}^{\mathrm{eq}}=\frac{\rho}{b}\left[1+u_{i}+G q_{i}\right]
$$

where

$$
\begin{gathered}
u_{i}=\frac{c_{i a} u_{a}}{c_{s}^{2}} \\
q_{i}=\frac{Q_{i a b} u_{a} u_{b}}{2 c_{s}^{4}}
\end{gathered}
$$

with

$$
Q_{i a b}=c_{i a} c_{i b}-c_{s}^{2} \delta_{a b}
$$

is the projetor along the $i$-th direction. Here,

$$
c_{s}=\frac{c}{\sqrt{D}}
$$

---

<!-- Page 30 -->

is the lattice sound speed. Notice the Galilean invariance is broken by this manner with

$$
G(\rho)=\frac{1-2 d}{1-d} \neq 1
$$

where $d=\rho / b$ is reduced fluid density. Now, equation 2.17 becomes a condition to match the tensor equality

$$
P_{a b}^{e}=\int f^{e} v_{a} v_{b} \mathrm{~d} \vec{v}=\rho\left(u_{a} u_{b}+c_{s}^{2} \delta_{a b}\right)
$$

This completes our formulation of LGCA. Its implementation has two major advantages:

1. Exact Computation (No round-off errors): The entire update scheme implies that the corresponding algorithm could be implemented in pure Boolean logic without the round-off error of floating point representation of real numbers. It also eases out the problem of numerical drifts in long-time simulations in fluid dynamics and statistical mechanics.
2. Virtually unlimited parallelism: LGCA are fairly well positioned in terms of low communicativity between processors since the collision step, the most time-consuming operation, is completely local and requires no communication among the processors.
On the other hand, LGCA has a number of flaws that ultimately ground it almost to a halt in early 90 s. We here briefly discusses these diseases and how to cure them: This formal derivation of "Navier-Stokes look-alike" equations from LGCA microdynamics highlight two basic diseases:
3. Lack of Galiean Invariance
4. Anomalous velocity depedence on the fluid pressure

These track to the fact that only a finite number of speeds are allowed. As a result, the equilibria can only be deinfed with a perturbative expansion in the flow field This fails to reproduce the Navier-Stokes inertial and pressure tensors. It delivers instead the following
where

$$
\sum_{i} f_{i}^{e} c_{i a} c_{i b}=\rho\left[g u_{a} u_{b}+c_{s}^{2}\left(1-g M^{2}\right) \delta_{a b}\right]
$$

is the Galiean breaking factor.
Additional drawback of LGCA that comes are

1. Statistical Noises
2. Low Reynolds' Number
3. Expoential Complexity
4. Spurious Invariants

# 2.2 Lattice Boltzmann Method Overview 

Lattice Boltzmann Method(LBM) was first developed to cure the diease of statistical noises in LGCA[11]. However, it proved itself to be a more effective method than LGCA in terms of reproducing the correct hydrodynamics and still possesses the advantages of LGCA, especially on its compatibility with parallel computation. Normally, conventional solvers of Navier Stokes equation include the boundary information exchange of both diffusion equation and advection equation. Solving these two equation separately is also not an easy task. This requirement of separation greatly adds to the burden of computation. Later in the chapter, we will see that LBM is essentially equivalent to a forward Euler's method that although being unstable in some regime, it is highly efficient in terms of computation and parallelization.

The basic quantity of LBM is the discrete velocity function $f_{i}(\boldsymbol{x}, t)$. It represents the density of particles in the phase space with velocity $\boldsymbol{c}_{i}=\left\{c_{i x}, c_{i y}, c_{i z}\right\}$ at position $\boldsymbol{x}$ and time $t$. It satisfies the basic conservation relation

$$
\rho(\boldsymbol{x}, t)=\sum_{i} f_{i}(\boldsymbol{x}, t) \quad \rho \boldsymbol{u}(\boldsymbol{x}, t)=\sum_{i} \boldsymbol{c}_{i} f_{i}(\boldsymbol{x}, t)
$$

---

<!-- Page 31 -->

Unlike the usual continuous population, $f_{i}(\boldsymbol{x}, t)$ is discrete in space, time and velocity. For space and time, the discretization is carried out with time step $\Delta t$ and lattice spacing $\Delta x$. In traditional LBM literatures, the most common choice for them is a lattice unit: $\Delta t=\Delta x=1$. As long as the non-dimensional numbers in the system solved numerically and the real simulated system is the same, the behavior of two system will be the same, according to the law of similarity.

Discrete velocity $\boldsymbol{c}_{i}$ shows how populations are discretized in velocity space. They essentially restrict the motion of population in only the direction of $\boldsymbol{c}_{i}$. Together, they come with a weight and forms a discrete velocity set $\left\{\boldsymbol{c}_{i}, w_{i}\right\}$. These sets are denoted in literature as DdQq, where $d$ is the number of physical dimensions and $q$ is the number of discrete velocities.

By discretizing the Boltzmann's equation in velocity, space, and time, we obtain the lattice Boltzmann equation:

$$
f\left(\boldsymbol{x}+\boldsymbol{c}_{i} \Delta t, t+\Delta t\right)=f_{i}(\boldsymbol{x}, t)+\Omega_{i}(\boldsymbol{x}, t)
$$

This shows the particle distribution moves from a point $(\boldsymbol{x}, t)$ to another adjacent point $\left(\boldsymbol{x}+\boldsymbol{c}_{i} \Delta t\right.$, $t+\Delta t)$ but at the same time is affect by a particle collision $\Omega_{i}(\boldsymbol{x}, t)$. For many applications, we resort to the most simplest collision integral $\Omega_{i}(\boldsymbol{x}, t)$, the Bhatnagar-Gross-Krook(BGK) operator.

$$
\Omega_{i}(\boldsymbol{x}, t)=-\frac{f_{i}-f_{i}^{\mathrm{eq}}}{\tau} \Delta t
$$

BGK operator essentially relaxes the population $f_{i}$ towards an equilibrium $f_{i}^{\text {eq }}$ at a rate determined by $\tau$. This equilibrium is given by

$$
f_{i}^{\mathrm{eq}}=w_{i} \rho\left(1+\frac{\boldsymbol{u} \cdot \boldsymbol{c}_{i}}{c_{e}^{2}}+\frac{\left(\boldsymbol{u} \cdot \boldsymbol{c}_{i}\right)^{2}}{2 c_{e}^{4}}-\frac{\boldsymbol{u} \cdot \boldsymbol{u}}{2 c_{e}^{2}}\right)
$$

where the weights $w_{i}$ are contained in the discrete velocity set. With Chapman-Enskog analysis, it can be proved that with a correct chocie of $\tau$, equation 2.35 reproduce Navier-Stokes equation. This choice is related to the kinematic viscosity $\nu$

$$
\nu=c_{e}^{2}\left(\tau-\frac{\Delta t}{2}\right)
$$

In summary, there are two parts in LBM that are performed in succession

1. Collision(Relaxation):

$$
f_{i}^{\star}(\boldsymbol{x}, t)=f_{i}(\boldsymbol{x}, t)-\frac{f_{i}(\boldsymbol{x}, t)-f_{i}^{\mathrm{eq}}(\boldsymbol{x}, t)}{\tau} \Delta t
$$

where $f_{i}^{\star}(\boldsymbol{x}, t)$ is a post-collision state.
2. Streaming:

$$
f\left(\boldsymbol{x}+\boldsymbol{c}_{i} \Delta t, t+\Delta t\right)=f_{i}^{\star}(\boldsymbol{x}, t)
$$

it brings the post-collision population to the adjacent cell.

# 2.3 Discretization of Velocity Space 

To derive the above numerical scheme, the most important step is discretization of velocity space. A continuous Boltzmann equation is difficult to solve because it is an advection equation of 7 variables in 3D. A correct discretization that still reproduces the hydrodynamics can reduce this deegree of freedom.

In order to do this, we notice the similarity between Maxwell-Boltzmann Distribution and the generating function of hermite polynomials. It can be shown that Maxwell-Boltzmann Distribution can be reduced to a truncated sum of Hermite Polynomials while retaining the correct macroscopic moments. Then, the macroscopic moments can be evaluated by the method of Gaussian-Hermite Quadrature which uses a discrete sum over the polynomial integrand at specific abscissae. We will show that these abscissae corresponds to the discrete velocities and their corresponding weight when used for evaluating macroscopic moments corresponds to the weights in the discrete velocity set.

---

<!-- Page 32 -->

To do this, we first perform a non-dimensionalization of force-free Boltzmann Equation

$$
\frac{\partial f}{\partial t}+\xi_{\alpha} \frac{\partial f}{\partial \alpha}=\Omega(f)
$$

and non-dimensionalized Maxwell-Boltzmann Distribution is

$$
f^{\mathrm{eq}}(\rho, \boldsymbol{u}, \theta, \boldsymbol{\xi})=\frac{\rho}{(2 \pi \theta)^{d / 2}} e^{-(\boldsymbol{\xi}-\boldsymbol{u})^{2} / 2 \theta}
$$

where all quantities $\rho, \boldsymbol{u}, \theta, \boldsymbol{\xi}$ (density, bulk velocity, temperature, micoscopic velocity) are nondimensional. Compare this equilibrium distribution with the weight(generating) function of Hermite Polynomials

$$
\omega(x)=\frac{1}{\sqrt{2 \pi}} e^{-x^{2} / 2}
$$

We can express Maxwell-Boltzmann Distribution in terms of (multi-dimensional) weight function

$$
f^{\mathrm{eq}}(\rho, \boldsymbol{u}, \theta, \boldsymbol{\xi})=\frac{\rho}{(2 \pi \theta)^{d / 2}} e^{-(\boldsymbol{\xi}-\boldsymbol{u})^{2} / 2 \theta}=\frac{\rho}{\theta^{d / 2}} \omega\left(\frac{\boldsymbol{\xi}-\boldsymbol{u}}{\sqrt{\theta}}\right)
$$

Now we can apply the Hermite expansion to equilibirum distribution

$$
f^{\mathrm{eq}}(\rho, \boldsymbol{u}, \theta, \boldsymbol{\xi})=\omega(\boldsymbol{\xi}) \sum_{n=0}^{\infty} \frac{1}{n!} \boldsymbol{a}^{(n), \mathrm{eq}}(\rho, \boldsymbol{u}, \theta) \cdot \boldsymbol{H}^{(n)}(\boldsymbol{\xi})
$$

where $\boldsymbol{H}^{(n)}(\boldsymbol{\xi})$ is a multi-dimensional Hermite Polynomial

$$
\boldsymbol{H}^{(n)}(\boldsymbol{\xi})=(-1)^{n} \frac{1}{\omega(\boldsymbol{\xi})} \nabla^{(n)} \omega(\boldsymbol{\xi})
$$

and the coefficients $\boldsymbol{a}^{(n), \mathrm{eq}}(\rho, \boldsymbol{u}, \theta)$ is given by

$$
\boldsymbol{a}^{(n), \mathrm{eq}}(\rho, \boldsymbol{u}, \theta)=\int f^{\mathrm{eq}}(\rho, \boldsymbol{u}, \theta, \boldsymbol{\xi}) \boldsymbol{H}^{(n)}(\boldsymbol{\xi}) \mathrm{d}^{d} \xi=\rho \int \omega(\boldsymbol{\eta}) \boldsymbol{H}^{(n)}(\sqrt{\theta} \boldsymbol{\eta}+\boldsymbol{u}) \mathrm{d}^{d} \eta
$$

Now, computing these integrals, we find that the first few coefficient corresponds to the conserved moments

$$
\begin{aligned}
a^{(0), \mathrm{eq}} & =\rho \\
a^{(1), \mathrm{eq}} & =\rho u_{\alpha} \\
a^{(2), \mathrm{eq}} & =\rho\left(u_{\alpha} u_{\beta}+(\theta-1) \delta_{\alpha \beta}\right)
\end{aligned}
$$

Hence, we do not need the full form of Maxwell-Boltzmann Equation to preserve the conserved moments. A truncated version of it up to the third coefficients is enough to reproduces density, momentum and energy.

$$
\begin{aligned}
f^{\mathrm{eq}}(\rho, \boldsymbol{u}, \theta, \boldsymbol{\xi}) & =\omega(\boldsymbol{\xi}) \rho\left(1+\xi_{\alpha} u_{\alpha}+\left(u_{\alpha} u_{\beta}+(\theta-1) \delta_{\alpha \beta}\right)\left(\xi_{\alpha} \xi_{\beta}-\delta_{\alpha \beta}\right)\right) \\
& =\omega(\boldsymbol{\xi}) \rho Q(\boldsymbol{u}, \theta, \boldsymbol{\xi})
\end{aligned}
$$

Using this new definition of equilibrium distribution, we can compute the conserved moments again. However, since it now takes a polynomial form, the intergal can be evaluated by using GaussianHermite Quadrature.

$$
\begin{aligned}
\boldsymbol{a}^{(n), \mathrm{eq}}(\rho, \boldsymbol{u}, \theta) & =\int f^{\mathrm{eq}}(\rho, \boldsymbol{u}, \theta, \boldsymbol{\xi}) \boldsymbol{H}^{(n)}(\boldsymbol{\xi}) \mathrm{d}^{d} \xi \\
& =\int \omega(\boldsymbol{\xi}) \rho Q(\boldsymbol{u}, \theta, \boldsymbol{\xi}) \boldsymbol{H}^{(n)}(\boldsymbol{\xi}) \mathrm{d}^{d} \xi \\
& =\rho \sum_{i=1}^{n} w_{i} Q\left(\boldsymbol{\xi}_{i}\right) \boldsymbol{H}^{(n)}\left(\boldsymbol{\xi}_{i}\right)
\end{aligned}
$$

---

<!-- Page 33 -->

By the rules of Gaussian Quadrature, for a polynomial of order $N$, it needs at least $n \geqslant(N+$ 1) $/ 2$ terms to calculate the moments exactly. Therefore, the choice of $n$ affects the number of abscissae $\xi_{i}$, which is at least 3 since $Q$ is of order 2 . These abscissae are then given by the roots of correpsonding order Hermite Polynomials $H^{(n)}(\xi)$. Define $f_{i}^{\mathrm{eq}}\left(\rho, \boldsymbol{u}, \theta, \boldsymbol{\xi}_{i}\right)$ as the distribution evaluated at abscissae $\boldsymbol{\xi}_{i}$, now we have a discrete set of equilibirum distributions that also satisfies the conservation laws

$$
f_{i}^{\mathrm{eq}}\left(\rho, \boldsymbol{u}, \theta, \boldsymbol{\xi}_{i}\right)=w_{i} \rho\left(1+\xi_{i \alpha} u_{\alpha}+\left(u_{\alpha} u_{\beta}+(\theta-1) \delta_{\alpha \beta}\right)\left(\xi_{i \alpha} \xi_{i \beta}-\delta_{\alpha \beta}\right)\right)
$$

Many of these abscissae carry a factor of $\sqrt{3}$. We can introduce a new set of particle velocities to eliminate this. These are the set of discrete velocities.

$$
\boldsymbol{c}_{i}=\frac{\boldsymbol{\xi}_{i}}{\sqrt{3}}
$$

Using the new definition, we obtain the form of equilibrium distribution mentioned in the previous section:

$$
f_{i}^{\mathrm{eq}}=w_{i} \rho\left(1+\frac{c_{i \alpha} u_{\alpha}}{c_{s}^{2}}+\frac{u_{\alpha} u_{\beta}\left(c_{i \alpha} c_{i \beta}-c_{s}^{2} \delta_{\alpha \beta}\right)}{2 c_{s}^{4}}\right)
$$

where $c_{s}$ is the speed of sound. In a similar manner, we can also discretize the particle distribution

$$
f_{i}(\boldsymbol{x}, t)=\frac{w_{i}}{\omega\left(\boldsymbol{c}_{i}\right)} f\left(\boldsymbol{x}, \boldsymbol{c}_{i}, t\right)
$$

Substituting this back to the boltzmann equation, we have the discrete version of Boltzmann's Equation

$$
\partial_{t} f_{i}+c_{i \alpha} \partial_{\alpha} f_{i}=\Omega\left(f_{i}\right)
$$

The choice of discrete velocity set often depends on the application. However, it must be rotationally invariant to be physical. This limits the choices of discrete velocities by adding new conditions. If the conditions are satisfied, the extra degree of freedom can be used for other purposes. The most common chocie of them are $D 1 Q 3, D 2 Q 9, D 3 Q 19$.

# 2.4 Discretization in Space and Time 

For the discrete version of Boltzmann Equation, we can solve it using the Method of Characteristics. The solution is

$$
f_{i}\left(\boldsymbol{x}+\boldsymbol{c}_{i} \Delta t, t+\Delta t\right)-f_{i}(\boldsymbol{x}, t)=\int_{0}^{\Delta t} \Omega_{i}\left(\boldsymbol{x}+\boldsymbol{c}_{i} \zeta, t+\zeta\right) \mathrm{d} \zeta
$$

To integrate the right hand side of the equation, LBM chooses a explicit Forward Euler scheme for the space-time integration. This greatly simplifies the LBM scheme but also puts stability criterion on the numerical method.

$$
f_{i}\left(\boldsymbol{x}+\boldsymbol{c}_{i} \Delta t, t+\Delta t\right)-f_{i}(\boldsymbol{x}, t)=\Delta t \Omega_{i}(\boldsymbol{x}, t)
$$

this is the Lattice Boltzmann Equation. Substituting the BGK operator in the equation, we obatin the LBGK equation

$$
f_{i}\left(\boldsymbol{x}+\boldsymbol{c}_{i} \Delta t, t+\Delta t\right)-f_{i}(\boldsymbol{x}, t)=-\frac{\Delta t}{\tau}\left(f_{i}(\boldsymbol{x}, t)-f_{i}^{\mathrm{eq}}(\boldsymbol{x}, t)\right)
$$

Depending on the value of $\Delta t / \tau$, the relaxation can be characterized by one of the following ways

1. Under-relaxation: $\Delta t / \tau>1$, the population $f_{i}$ decays exponetially towards the equilibrium distribution.
2. Full-relaxation: $\Delta t / \tau=1, f_{i}$ decays directly to $f_{i}^{\text {eq }}$.
3. Over-relaxation: $1 / 2<\Delta t / \tau<1$, the population $f_{i}$ oscillates around the equilibrium distribution but also decays exponentially in terms of amplitude.

---

<!-- Page 34 -->

We can also see from here that $1 / 2<\Delta t / \tau$ is a necessary condition for numerical stability. For cases $1 / 2>\Delta t / \tau, f_{i}$ will oscillates with an increasing amplitude. Now, we have successfully discretized Boltzmann Equation in velocity, space and time.

# 2.5 Chapman Enskog Analysis of LBM 

We now need to prove that LBE does simulate the behavior of fluids. This is done with ChapmanEnskog expansion, which is very similar to the procedure described in Chapter 1. We start with a perturbative expansion of $f_{i}$ with Knudsen's number Kn as the expansion parameter $\varepsilon$

$$
f_{i}=f_{i}^{\mathrm{eq}}+\varepsilon f_{i}^{(1)}+\varepsilon^{2} f_{i}^{(2)}+\cdots
$$

Similarly, we require that the higher order expansion does not contribute to mass and momentum individually at each order. The property of collision invariant also required that BGK collision does not alter the mass and momentum as well. These are the constriants of the problem. Let $f_{i}^{\text {neq }}=f_{i}-f_{i}^{\text {eq }}$

$$
\begin{aligned}
& \sum_{i} f_{i}^{\text {neq }}=0 \quad \sum_{i} \boldsymbol{c}_{i} f_{i}^{\text {neq }}=\mathbf{0} \\
& \sum_{i} f_{i}^{(n)}=0 \quad \sum_{i} \boldsymbol{c}_{i} f_{i}^{(n)}=\mathbf{0} \quad \text { for all } n \geqslant 1
\end{aligned}
$$

Taylor expanding the LBE, we obtain

$$
\Delta t\left(\partial_{t}+c_{i \alpha} \partial_{\alpha}\right) f_{i}+\frac{\Delta t^{2}}{2}\left(\partial_{t}+c_{i \alpha} \partial_{\alpha}\right)^{2} f_{i}+O\left(\Delta t^{3}\right)=-\frac{\Delta t}{\tau} f_{i}^{\text {neq }}
$$

We will discard the higher order terms on the basis that $\Delta t^{n}\left(\partial_{t}+c_{i \alpha} \partial_{\alpha}\right)^{n} f_{i}$ scale with $O\left(\mathrm{Kn}^{n}\right)$. Since finally in the chapman enskog expansion we will prove that the first order expansion is sufficient to reproduce Navier-Stokes equation. This omission of higher order terms in time can be justifed. The second order derivative term can also be subtracted form the equation by applying $\Delta t\left(\partial_{t}+c_{i \alpha} \partial_{\alpha}\right)$ on to the equation itself.

$$
\Delta t\left(\partial_{t}+c_{i \alpha} \partial_{\alpha}\right) f_{i}=-\frac{\Delta t}{\tau} f_{i}^{\text {neq }}+\Delta t\left(\partial_{t}+c_{i \alpha} \partial_{\alpha}\right) \frac{\Delta t}{2 \tau} f_{i}^{\text {neq }}
$$

Now we perform a multiscale expansion on the time derivative of the equation in interest.

$$
\Delta t \partial_{t} f_{i}=\Delta t\left(\varepsilon \partial_{t}^{(1)} f_{i}+\varepsilon^{2} \partial_{t}^{(2)} f_{i}+\cdots\right) \quad \Delta t c_{i \alpha} \partial_{\alpha} f_{i}=\Delta t\left(\varepsilon c_{i \alpha} \partial_{\alpha}^{(1)}\right) f_{i}
$$

Different order of Kn would now separate the equation

$$
\begin{aligned}
O(\varepsilon) & :\left(\partial_{t}^{(1)}+c_{i \alpha} \partial_{\alpha}^{(1)}\right) f_{i}^{\mathrm{eq}}=-\frac{1}{\tau} f_{i}^{(1)} \\
O\left(\varepsilon^{2}\right) & : \partial_{t}^{(2)} f_{i}^{\mathrm{eq}}+\left(\partial_{t}^{(1)}+c_{i \alpha} \partial_{\alpha}^{(1)}\right)\left(1-\frac{\Delta t}{2 \tau}\right) f_{i}^{(1)}=-\frac{1}{\tau} f_{i}^{(2)}
\end{aligned}
$$

Taking zeroth to second moments to the $O(\varepsilon)$ equation, we have the $O(\varepsilon)$ moments equation

$$
\begin{aligned}
\partial_{t}^{(1)} \rho+\partial_{\gamma}^{(1)}\left(\rho u_{\gamma}\right) & =0 \\
\partial_{t}^{(1)} \rho u_{\alpha}+\partial_{\beta}^{(1)} \Pi_{\alpha \beta}^{\mathrm{eq}} & =0 \\
\partial_{t}^{(1)} \Pi_{\alpha \beta}^{\mathrm{eq}}+\partial_{\gamma}^{(1)} \Pi_{\alpha \beta \gamma}^{\mathrm{eq}} & =-\frac{1}{\tau} \Pi_{\alpha \beta}^{(1)}
\end{aligned}
$$

The moments defined are

$$
\begin{aligned}
\Pi_{\alpha \beta}^{\mathrm{eq}} & =\sum_{i} c_{i \alpha} c_{i \beta} f_{i}^{\mathrm{eq}}=\rho u_{\alpha} u_{\beta}+\rho c_{s}^{2} \delta_{\alpha \beta} \\
\Pi_{\alpha \beta \gamma}^{\mathrm{eq}} & =\sum_{i} c_{i \alpha} c_{i \beta} c_{i \gamma} f_{i}^{\mathrm{eq}}=\rho c_{s}^{2}\left(u_{\alpha} \delta_{\beta \gamma}+u_{\beta} \delta_{\alpha \gamma}+u_{\gamma} \delta_{\alpha \beta}\right) \\
\Pi_{\alpha \beta}^{(1)} & =\sum_{i} c_{i \alpha} c_{i \beta} f_{i}^{(1)}
\end{aligned}
$$

---

<!-- Page 35 -->

Here, only keeping terms in the first order of perturbative expansion and reversing it allows us to get Euler's equation sets: the first one is continuity equation and the second one is Euler's momentum equation. Going to the next order of perturbative expansion, we can find the first and the second moment of $O\left(\varepsilon^{2}\right)$ equations

$$
\begin{aligned}
\partial_{t}^{(2)} \rho & =0 \\
\partial_{t}^{(2)} \rho u_{\alpha}+\partial_{\beta}^{(1)}\left(1-\frac{\Delta t}{2 \tau}\right) \Pi_{\alpha \beta}^{(1)} & =0
\end{aligned}
$$

These set of equation can be seen as correction to continuity equation and Euler's momentum equation shown in the first order moment equations. As expected, the second order does not contribute to the continuity equation, but the Euler's momentum equation gets a correction with known tensor $\Pi_{\alpha \beta}^{(1)}$. Summarizing all of the above, we have

$$
\begin{aligned}
\left(\varepsilon \partial_{t}^{(1)}+\varepsilon^{2} \partial_{t}^{(2)}\right) \rho+\varepsilon \partial_{\gamma}^{(1)}\left(\rho u_{\gamma}\right) & =0 \\
\left(\varepsilon \partial_{t}^{(1)}+\varepsilon^{2} \partial_{t}^{(2)}\right) \rho u_{\alpha}+\varepsilon \partial_{\beta}^{(1)} \Pi_{\alpha \beta}^{\mathrm{eq}} & =-\varepsilon^{2} \partial_{\beta}^{(1)}\left(1-\frac{\Delta t}{2 \tau}\right) \Pi_{\alpha \beta}^{(1)}
\end{aligned}
$$

Now an expression of $\Pi_{\alpha \beta}^{(1)}$ is needed to close the expansion. This is done by expanding isothermal equation of state and equilibrium distribution up to the order of $O\left(u^{2}\right)$

$$
\Pi_{\alpha \beta}^{(1)}=-\rho c_{s}^{2} \tau\left(\partial_{\beta}^{(1)} u_{\alpha}+\partial_{\alpha}^{(1)} u_{\beta}\right)+\tau \partial_{\gamma}^{(1)} \rho u_{\alpha} u_{\beta} u_{\gamma}
$$

In this expansion, the $O\left(u^{3}\right)$ order terms are neglected reasonably if we can assume that $u^{2} \ll c_{s}^{2}$. This is only true for a weakly compressible case $\mathrm{Ma}^{2} \ll 1$. Hence, LBGK is only valid for weakly compressible fluids. With this expression, the expansion is closed and we have found Navier-Stokes Equation by reversing the perturbative expansion

$$
\begin{aligned}
\partial_{t} \rho+\partial_{\gamma}\left(\rho u_{\gamma}\right) & =0 \\
\partial_{t}\left(\rho u_{\alpha}\right)+\partial_{\beta}\left(\rho u_{\alpha} u_{\beta}\right) & =-\partial_{\alpha} p+\partial_{\beta}\left[\eta\left(\partial_{\beta} u_{\alpha}+\partial_{\alpha} u_{\beta}\right)\right]
\end{aligned}
$$

where we have defined

$$
p=\rho c_{s}^{2} \quad \eta=\rho c_{s}^{2}\left(\tau-\frac{\Delta t}{2}\right)
$$

From the Chapman-Enskog analysis, we know that LBGK is only valid for weakly compressible fluids and its relaxation is controlled by viscosity and therefore the Reynolds' number of the problem. Since relaxation essentially determines the numerical stability of the scheme, we can say that LBGK's numerical stability is limited by Reynolds' number in the problem.

# 2.6 Stability Analysis of LBGK 

There are many literatures dedicated for the numerical stability analysis of LBM and LBGK in particular(see, for example, [9] ). Here, we state the two sufficient conditions for stability of LBM

1. $\Delta t / \tau>1 / 2$ as shown in section 2.4
2. Non-negativity of all equilibrium distribution $f_{i}^{\text {eq }}$

In particular, we notice that satisfying condition 1 does not automatically imply the validity of condition 2 and vice versa. In the over-relaxation region, even both positive $f_{i}^{\text {eq }}$ and $f_{i}$ cannot ensure the positivity of the post-collision state $f_{i}^{s}$. However, for high Reynolds number or low viscosity simulation, over-relaxation region is necessary. Therefore, usually LBGK fails at high Reynolds number. To see the reason of numerical stability of LBGK and how to modify the collision, we move to the next chapter on the Entropic Scheme of LBM.

---

<!-- Page 36 -->

# Chapter 3 

## Entropic methods

The first portion of the chapter which focuses on the discussion of relationship between numerical instability and H-theorem follows a textbook[12] and a published lecture note[13]. The discussion of multiple relaxation time follows a Lattice Boltzmann Method textbook[8]. The rest of the discussion on Entropic Lattice Boltzman Method is a summary of a variety of journal articles.

### 3.1 BGK and Entropy

Boltzmann's H-theorem states that quantity $H$ must decrease for an isolated system until it reaches its minimum. The quantity $H$ is a representation of entropy in thermodyanmics, whereas the second law of thermodyanmics comes as a consequence of H-theorem. $H$ is defined continuously as

The H-theorem then states

$$
H=\int f \ln f \mathrm{~d}^{d} v
$$

$$
\frac{d H}{d t} \leqslant 0
$$

For an ideal isolated gas, when $H$ reaches it minimum, the equilibirum distribution of the system is Maxwellian. In the lattice Boltzmann context where discrete velocity can only take a number of certain direction, the integration over velocity space becomes a sum over all possible microscopic velocities.

$$
H=\sum_{i=1}^{N} f_{i} \ln \left(f_{i}\right)
$$

It is well known that numerical stability is a issue for LBGK method on high Reynold's number. The reason behind instabilities is in LBGK's drastic simplification of the real physical process. It cannot ensure the validity of $H$ theorem during relaxation. The LBGK scheme is (time step $\Delta t$ is unity)

$$
f_{i}^{n+1}=\left(1-\frac{1}{\tau}\right) f_{i}^{n}+\frac{1}{\tau} f_{i}^{\mathrm{eq}, n}
$$

with the upscript representing the time step. The expression of equilibirum distribution is

$$
f_{i}^{\mathrm{eq}}=w_{i} \rho\left(1+\frac{c_{i \alpha} u_{\alpha}}{c_{s}^{2}}+\frac{u_{\alpha} u_{\beta}\left(c_{i \alpha} c_{i \beta}-c_{s}^{2} \delta_{\alpha \beta}\right)}{2 c_{s}^{4}}\right)
$$

H-theorem requires

$$
H_{i}^{n+1}-H^{n} \leqslant 0
$$

It could be shown that the relaxation scheme specified in 3.4 cannot ensure the validity of 3.6.

Consider the $H$-theorem for the BGK equation in the continuum. The simplest candidate for a BGK $H$-function is the "Euclidean distance":

$$
h[f]=f^{2}
$$

---

<!-- Page 37 -->

the BGK equation on convective equation has the form of

$$
\partial_{t} f+v_{a} \partial_{a} f=\omega(f-g)
$$

Taking dot product with $-f$ over the velocity space on both side:

$$
\partial_{t} s+\partial_{a} S_{a}=\omega[(f, f)-(f, g)]
$$

where the local entropy density $s$ and the corresponding entropy flux $S_{a}$ are

$$
\begin{aligned}
s(x, t) & =-\frac{1}{2}(f, f)=\int_{v} \frac{1}{2} f^{2} \mathrm{dv} \\
S_{a}(x, t) & =-\int_{v} v_{a} f^{2} \mathrm{dv}
\end{aligned}
$$

Here, $g$ is the attractor of the collision operator. It would be a target equilibrium of the relaxtation. A local $H$-theorem implies that the local entropy production

$$
\sigma=\omega[(f, f)-(f, g)]
$$

is positive for any functions of $f$ and $g$. The first term is clearly positive whereas the second term has no single sign a priori and prevents general conclusion on the existence of $H$-theorem. One exception is the case where target equilibrium coincides with the global uniform equilibrium $g=\rho / b$. In this case, $g$ reduces to the first term of the Hermite expansion $g=\rho h_{0}(v), h_{0}(v)=1$, and since Hermite basis are orthogonal, we have

$$
(f, g)=(g, g)
$$

The entropy production will be

$$
\sigma=\omega[(f-g, f-g)]
$$

This is clearly positive-definite. We can extend the idea to a more general attractor $g=g_{0} h_{0}(v)+$ $g_{1} h_{1}(v)$ and $f=f_{0} h_{0}(v)+f_{1} h_{1}(v)+\ldots$ (rest of the hermite expansion). We have the scalar product

$$
(g, g)=g_{0}^{2}+g_{0}^{2} \quad(f, g)=f_{0} g_{0}+f_{1} g_{1}
$$

Condition in 3.13 is satisfied if and only if $f_{0}=g_{0}$ and $f_{1}=g_{1}$. This is the usual requirement that actual distribution and target equilibirum have to share the same density $\rho_{0}$ and current density $\rho u_{a}$. Moving to the LBGK world, we have similar expression for 3.12 on global entropy. This is more complicated since now we have $f_{i}^{n}(\boldsymbol{x}, t)$ and $f_{i}^{n+1}=\left(\boldsymbol{x}+\boldsymbol{c}_{i}, t+1\right)$ which contribute separately to the local entropy production.

$$
\Delta S=-\left(H_{2}(t+1)-H_{2}(t)\right)=2 \omega \sum_{x} \sum_{i} 2\left(f_{i}^{n}-g_{i}^{n}\right)\left(f_{i}^{n}+f_{i}^{n+1}\right)
$$

where periodic boundary condition is assumed. Using LBGK update, this becomes (dropping all subscript)

$$
\Delta S=-\omega \sum_{x} \sum_{i}\left[\omega\left(f_{i}-g_{i}\right)^{2}-2 f_{i} g_{i}\right]
$$

Invoking a condition very similar to $3.13:\left(f_{i}, g_{i}\right)=\left(g_{i}, g_{i}\right)$

$$
\Delta S=\omega(\omega-2) \sum_{x} \sum_{i}\left(f_{i}-g_{i}\right)^{2}
$$

From here, we know the $H$-theorem is ensured if $0<\omega<2$. This is the same condition derived from the linear analysis of LBGK equation where it produces a stability condition similar to the forward Euler scheme

$$
|1-\omega|<1
$$

Notice this conclusion from 3.18 does not apply in the over relaxation regime where $1<\omega<2$. From the LBGK update, it is clear that a positive pair $f_{i}^{n}, g_{i}^{n}$ does not lead to a positive $f_{i}^{n+1}$. This is precisely the region that is most relevant to turbulent flows: for $\omega=1, \nu=c_{s}^{2} / 2$ which is well above 0.1 and needs an unaffordable resolution to reach for such high Reynolds' number.

---

<!-- Page 38 -->

Nonlinear stability is more complicated: it holds in a restricted strip $0<\omega<\omega_{N L}[f]<2$, where the threshold $\omega_{N L}[f]$ is an unknown functional of the solution $f$. It is almost impossible to get an exact form of $\omega_{N L}$ but the actual value can be computed by requiring $H$-function or, more generally, Lyapunov functional $\Lambda[f]$ should not decrease by collisions. Actual value of the relaxation parameter can be obtained by solving the (nonlinear) algebraic problem

$$
\Lambda[f]=\Lambda\left[f^{*}\right]
$$

where $f^{*}$ is the post-collision state such that entropy is conserved. The update $f^{\prime}$ then can be obtained by a linear interpolation between $f$ and $f^{*}$.

$$
f^{\prime}=\beta f+(1-\beta) f^{*}
$$

this new relaxation scheme leads to the entropic LBE methods.

# 3.2 MRT 

One way of improving numerical stability is to construct collision operators that have more degress of freedom than BGK. These are multiple-relaxation-time method(MRT) and two-relaxation-time method(TRT). They stem from the observation that not all population relaxes towards equilibrium in the same rate. In BGK, there is only one collision time $\tau$ to control the relaxation of all population. It is then possible to propose different collision time for different population. Let $\omega=1 / \tau$ be the BGK collision rate and $\omega_{i}$ be the collision rate of each population $f_{i}$. They still needs to satisfy conservation of density of and momentums.

$$
\begin{aligned}
\sum \Omega_{i}=-\sum \omega_{i}\left(f_{i}-f_{\mathrm{eq}}\right) & =0 \\
\sum \Omega_{i} \boldsymbol{c}_{i}=-\sum \omega_{i}\left(f_{i}-f_{\mathrm{eq}}\right) \boldsymbol{c}_{i} & =\mathbf{0}
\end{aligned}
$$

For a model with distinct value for each microscopic population (called MRT-L model), the situation might be complicated. Another idea is to use the velocity moments.

$$
\begin{aligned}
\sum \Omega_{i}=-\sum \omega\left(f_{i}-f_{\mathrm{eq}}\right) & =0 \\
\sum \Omega_{i} \boldsymbol{c}_{i}=-\sum \omega^{\prime}\left(f_{i}-f_{\mathrm{eq}}\right) \boldsymbol{c}_{i} & =\mathbf{0}
\end{aligned}
$$

In principle, each moment (density, momentum, and etc) can be relaxed with their individual rate ( $\omega, \omega^{\prime}$ for density and momentum here) such that their relaxation rate can be individually controlled. This is the basic idea behind MRT. Hence, to operate relaxation of each moment, MRT collisions have to be performed in moment space. The detail steps would be 1)mapping the populations $f_{i}$ to moment space, 2) perform the collision, 3) map the relaxed moment back to the population space. To start this, we notice moments can be represented as certain summations over populations with Hermite polynomials

$$
\boldsymbol{a}^{(n)}=\sum_{i} \boldsymbol{H}_{i}^{(n)} f_{i}
$$

A similar definition could be devised for moment $m_{k}$ in a DdQq velocity set.

$$
m_{k}=\sum_{i=0}^{q-1} M_{k i} f_{i}
$$

where $M$ is a $q \times q$ matrix. Generally, $q$ moments $m_{k}$ can be generated from $q$ populations $f_{i}$ through linear transformation from population space to the moment space. By carefully choosing the transformation $M$, the obtained moments $m_{k}$ can be made directly corresponding to the hydrodynamic moments. From this ,we can write MRT collision operator in detail: from the BGK collision we have

$$
\boldsymbol{f}^{n+1}-\boldsymbol{f}^{n}=-\omega\left(\boldsymbol{f}^{n}-\boldsymbol{f}^{o q}\right) \Delta t
$$

---

<!-- Page 39 -->

The expression is not changed when multiply the right side by $\boldsymbol{I}=\boldsymbol{M}^{-1} \boldsymbol{M}$, assuming it can be inverted:

$$
\begin{aligned}
\boldsymbol{f}^{n+1}-\boldsymbol{f}^{n} & =-\boldsymbol{M}^{-1} \boldsymbol{M} \omega\left(\boldsymbol{f}^{n}-\boldsymbol{f}^{\mathrm{eq}}\right) \Delta t \\
& =-\boldsymbol{M}^{-1} \omega\left(\boldsymbol{M} \boldsymbol{f}^{n}-\boldsymbol{M} \boldsymbol{f}^{\mathrm{eq}}\right) \Delta t \\
& =-\boldsymbol{M}^{-1} \omega \boldsymbol{I}\left(\boldsymbol{m}^{n}-\boldsymbol{m}^{\mathrm{eq}}\right) \Delta t \\
& =-\boldsymbol{M}^{-1} \boldsymbol{S}\left(\boldsymbol{m}^{n}-\boldsymbol{m}^{\mathrm{eq}}\right) \Delta t
\end{aligned}
$$

Where we introduced $\boldsymbol{S}=\omega \boldsymbol{I}$. The above equation has three parts: $\boldsymbol{S}\left(\boldsymbol{m}^{n}-\boldsymbol{m}^{\mathrm{eq}}\right)$ is the collision in momentum space, $\boldsymbol{M}^{-1}$ transform from moment space to population space and the left hand side is the usual streaming step. Now, we can introduce different relaxation rate for each moment by redefining the diagonal matrix $\boldsymbol{S}$ :

$$
S_{k k}=\omega_{k}
$$

where $k$ goes from 0 to $q-1$. There are two ways to construct the transformation matrix $\boldsymbol{M}$. The first is to use Hermite Polynomials since we have shown in the discretization of velocity space that all moments can be expressed in terms of Hermite Polynomials. Another more widely used method is Gram-Schmidt procedure. It starts with vectors of known moments. The next step is to take a combination of velocity vectors $c_{i \alpha}$ of appropriate order and find coefficients to construct vector that is orthogonal to previous one. Repeating this, we can find matrix $M_{k i}$ with row $\boldsymbol{M}_{k}$. To summarize, the entire MRT algorithm is

1. Compute conserved moments

$$
\rho=\sum_{i} f_{i}, \rho \boldsymbol{u}=\sum_{i} f_{i} \boldsymbol{c}_{i}
$$

2. Transform to moment space

$$
\boldsymbol{m}=\boldsymbol{M} \boldsymbol{f}
$$

3. Compute equilibrium moments (the method to compute coefficients will be discussed later)

$$
m_{k}^{\mathrm{eq}}=\rho \sum_{l, m, n} a_{k, l m n} u_{x}^{l} u_{y}^{m} u_{z}^{n}
$$

4. Collide

$$
m_{k}^{\prime}=m_{k}-\omega_{k}\left(m_{k}-m_{k}^{\mathrm{eq}}\right) \Delta t
$$

5. Transform to population space

$$
f_{i}^{\prime}=\sum M_{i k} m_{k}
$$

6. Stream

$$
\boldsymbol{f}^{n+1}=\boldsymbol{f}^{\prime}
$$

# 3.3 Basic Entropic Scheme 

To start the entropic scheme (first described in [2]), we distinguish between three stages of relaxation:

1. Start of the relaxation, when population starts moving away from non-equilibrium: $f_{0}$
2. Intermediate stage
3. Final regression towards equilibrium $f^{0}$

Now, let the collision integral computed at initial stages be $Q_{0}=Q\left(f_{0}\right)$ and let

$$
f(\Gamma, a)=f_{0}(\Gamma)+a Q_{0}(\Gamma)
$$

---

<!-- Page 40 -->

be all the intermediate steps toward equilibrium. The approximation is only valid if the scalar variable $a \geqslant 0$ is not too large to exceed a certain limit $a^{*}$. This limit is important for entropic schemes. Now, let $S(f)$ be the entropy at state $f(\Gamma)$, and let $S(\Gamma, a)$ be the entropy at state $f(\Gamma, a)$. From $H$-theorem, the state is only accessible if the entropy production grows in $\left[0, a^{*}\right][3]$ :

1. Entropy $S(a)$ increase from 0 to $a^{*}$
2. Entropy $S(a)$ decrease as $a$ exceeds $a^{*}$

Determining the values for $a^{*}$ is the primary objective of the method. To do this, entropic method introduces two novel features:

1. Definition of the local equilibrium as the minimizer of the $H$-function.
2. Computation of the relaxation parameter to satisfy the monoticity of $H$-function.

# 3.3.1 $\boldsymbol{H}$-function 

We first must construct the $H$-function for entropic version of LBGK, which is based on the following observation:
"If discrete velocities are constructed from the zeros of Hermite polynomials, the method of discrete velocity is essentially the same as Grad's moment method based on the expansion of distribution function around a fixed Maxwellian distribution."[13][1]

To link discrete model to entropic Grad's moment method[7], first we note the continuous Boltzmann's $H$-function is a convex function given by

$$
H=\int F(\boldsymbol{x}, \boldsymbol{c}, t) \ln F(\boldsymbol{x}, \boldsymbol{c}, t) \mathrm{d} \boldsymbol{c}
$$

For 1D isothermal flows, a discrete form of the $H$-function could be directly computed from above[6]. For higher dimensional flows, this is accomplished by using Gauss-Hermite quadrature formulas. (we know that Hermite coefficients are directly related to the macroscopic variables) Therefore, we can obtain

$$
H_{\left(\boldsymbol{w}_{i}, \boldsymbol{c}_{i}\right)}=\sum_{i=1}^{b} f_{i} \ln \left(\frac{f_{i}}{w_{i}}\right)
$$

where weights $w_{i}$ is associated with the $i$-th particles with discrete velocity $\boldsymbol{c}_{i}$ and the total number of links is $b$. Furthermore, in $D$-dimensions, the one particle distribution function $F(\boldsymbol{x}, \boldsymbol{c}, t)$ is related to the $i$-th particle distribution function by

$$
f_{i}(\boldsymbol{x}, t)=w_{i}\left(2 \pi T_{0}\right)^{D / 2} \exp \left[\frac{\boldsymbol{c}_{i}^{2}}{2 T_{0}}\right] F(\boldsymbol{x}, \boldsymbol{c}, t)
$$

where $T_{0}$ is the reference temperature.

### 3.3.2 Equilibrium distribution

The equilibrium distribution is obtained from solving the minimization problem[1]:

$$
f_{i}^{\mathrm{eq}}=w_{i} \rho \prod_{j=1}^{D}\left(2-\sqrt{1+u_{j}^{2}}\right)\left(\frac{2 u_{j}+\sqrt{1+3 u_{j}^{2}}}{1-u_{j}}\right)^{\boldsymbol{c}_{i j} / \boldsymbol{c}}
$$

where $j$ is the index of spatial directions.

### 3.3.3 The relaxation procedure

The monoticity constraint on $H$-function is imposed with a two-step constraint:

1. In the first step, population is changed based on the bare collision $\boldsymbol{\Delta}=\boldsymbol{f}^{\mathrm{eq}}-\boldsymbol{f}$ such that the $H$-function value does not change.

---

<!-- Page 41 -->

2. In the section step, dissipation is introduced and value of $H$ is decreased.

A collision integral is defined on the kinetic equation:

$$
\boldsymbol{\Delta}(\boldsymbol{f})=\boldsymbol{f}^{\mathrm{eq}}-\boldsymbol{f}
$$

it is admissible if it conserves mass and momentum, and entropy production inequality is satisfied:

$$
\langle\nabla H \mid \Delta\rangle \leqslant 0
$$

where $\nabla H=\left\{\partial H / \partial x_{i}\right\}_{i=1, \ldots, D}$ is the spatial gradient of the enrtopy function. Equality hold if $\boldsymbol{f}^{\mathrm{eq}}=\boldsymbol{f}$. We can use this to obtain a constant value for entropy production at each time step. Introduce an auxiliary population $\boldsymbol{f}^{*}$ :

$$
\boldsymbol{f}^{*}=\boldsymbol{f}_{\alpha^{*}}=\boldsymbol{f}+\alpha^{*} \Delta
$$

where the scalar parameter $\alpha^{*}$ is the upper limit for the updating rule. It is the solution of the equation:

$$
H(\boldsymbol{f})=H(\boldsymbol{f}+\alpha \Delta)
$$

We can characterize solutions of the equation to the following:

1. $\alpha_{1}=0$
2. $\alpha_{2}$ can assume the following values:
i. the degenerate case $\alpha_{1}=\alpha_{2}=0$, occurs only if $\boldsymbol{f}=\boldsymbol{f}^{\mathrm{eq}}$.
ii. the boundary case, where

$$
\alpha_{2}=\alpha^{*}=\min _{i=0, \ldots, b-1 ; \Delta_{i}<0}\left\{\frac{\boldsymbol{f}_{i}}{\left|\Delta_{i}\right|}\right\}
$$

iii. the standard case, wherre the $H$-theorem(3.42) equation is solved with a root finder. Now, from the auxilary population, the collision is set as

$$
\boldsymbol{f}(\beta)=(1-\beta) \boldsymbol{f}+\beta \boldsymbol{f}^{*}
$$

where $\beta \in[0,1]$ is a fixed parameter that controls the viscosity coefficient in ths Navier-Stokes equation. It is given by [1]

$$
\nu=c_{s}^{2}\left(\frac{1}{\alpha \beta}-\frac{1}{2}\right) \delta t
$$

These two steps describe ELBM and give it full control of viscosity,

# 3.3.4 Remarks on the progression of $\boldsymbol{H}$-function 

The population $\boldsymbol{f}_{\alpha}$ can be expressd in a simple way around the local minimum of the $H$-function:

$$
\boldsymbol{f}_{\alpha}=\boldsymbol{f}+\alpha \Delta=\boldsymbol{f}+\alpha\left(\boldsymbol{f}^{\mathrm{eq}}-\boldsymbol{f}\right)=(1-\alpha) \boldsymbol{f}+\alpha \boldsymbol{f}^{\mathrm{eq}}
$$

Geometrical intepretation of the monotonicity of the $H$-function around its local minium gives the following remarks:

1. The old population $\boldsymbol{f}$ limits from below the new value.
2. The equilbrium $\boldsymbol{f}^{\text {eq }}$ fixes a central value for $\alpha$.

The lower and higher bounds, given for the populations. are associated with the value of parameter $\alpha$, i.e. $\alpha=0$ and $\alpha=1$ respectively, which implies $\alpha \in[0, a]$ with $a>1$. We can divide all the cases of possible value of $\alpha$ :

1. $\alpha=0$, static solution, we have $\boldsymbol{f}_{\alpha}=\boldsymbol{f}$.
2. $\alpha \in(0,1)$, the value of $H$-function is decreasing.

---

<!-- Page 42 -->

3. $\alpha=1$, a population that is the same as $\boldsymbol{f}^{\text {eq }}$. Local minimum of $H$-function is achieved.
4. $\alpha \in(1,2)$, the value of $H$-function is increasing.
5. $\alpha=2$, corresponding to a solution of equation 3.42 that gives the following value for population $\boldsymbol{f}_{\alpha}$ :

$$
\boldsymbol{f}_{\alpha}=\boldsymbol{f}^{\mathrm{eq}}-\boldsymbol{f}^{\mathrm{neq}}>0, \text { if } \boldsymbol{f}^{\mathrm{eq}}>\boldsymbol{f}^{\mathrm{neq}}
$$

For each value of $\alpha \in[0, a]$, we also have a different progression of $H$-function:

1. $H$ is locally quadratic: value $\alpha=1$ corresponds to the abscissa of the vertex of parabola
2. $H$ is not locally quadratic: $\alpha<2$ or $\alpha>2$, which matches to under and over relaxation of the scheme.

# 3.4 Basic Information Geometry 

Before introducing a new entropic scheme, the distinguishability of probability distributions has to be addressed first, which is in the realm of information geometry. Here we follow a conference paper[4]. More specifically, the notion of distances between two different probability distributions has to be introduced. For a parametric family of probability distributions, it is a set of distributions $p_{\theta}(x)$ labeled by $\theta=\left(\theta^{1} \ldots \theta^{n}\right)$. This forms a statistical manifold, namely, a space where each point, labeled by $\theta$, represents a probability distribution. Unlike generic manifolds, statistical manifolds possess a natural notion of distance - the information metric. This is inevitable resulting from the definition. In other words, Geometry is instrinsic to the structure of statistical manifolds.

The distance $\mathrm{d} l$ between two neighboring points $\theta$ and $\theta+\mathrm{d} \theta$ is given by Pythagora's theorem which, written in the form of a metric tensor $g_{a b}$, is

$$
\mathrm{d} l^{2}=g_{a b} \mathrm{~d} \theta^{a} \mathrm{~d} \theta^{b}
$$

N.Cecov states that the metric $g_{a b}$ on the manifold of probability distribution is essentially unique: up to an overall scale factor there is only one metric that takes into account the fact that these are not distances between simple structureless dots but distances between probability distributions. One result from this has to be emphasized: having a notion of distance means there is a notion of volume, and this implies there is an unique and objective notion of distribution that is uniform over the space of parameters - equal volumes are assigned equal probabilities.

An $n$ - dimensional manifold $\mathcal{M}$ is a smooth, possibly curved, space that is locally like $\mathbb{R}^{n}$. Therefore, it is possible to set up a coordinate frame (a mapping $\mathcal{M} \rightarrow \mathbb{R}^{n}$ ) such that each point $\theta$ is identified with a coordinate. For information manifolds, each point would be a probability distribution $p_{\theta}(x)$. A very convenient notion is $p_{\theta}(x)=p(x \mid \theta)$ : for example, we have for multinomial distribution, multivariate Gaussian and the canonical distribution:

$$
p\left(\left\{n_{i}\right\} \mid \theta\right)=\frac{N!}{n_{1}!n_{2}!\ldots n_{m}!}\left(\theta^{1}\right)^{n_{1}}\left(\theta^{2}\right)^{n_{2}} \ldots\left(\theta^{m}\right)^{n_{m}}
$$

where $N=\sum_{i=1}^{m} n_{i}$ and $1=\sum_{i=1}^{m} \theta_{i}$. This forms a statistical manifold of dimension $(m-1)$ called a Simplex, $S_{m-1}$. The parameters $\theta$ are a convenient choice of coordinates.

$$
p(x \mid \mu, \sigma)=\frac{1}{\left(2 \pi \sigma^{2}\right)^{n / 2}} \exp \left[-\frac{1}{2 \sigma^{2}} \sum_{a=1}^{n}\left(x^{a}-\mu^{a}\right)^{2}\right]
$$

where means are $\mu^{a}, a=1, \ldots n$ and variance is $\sigma^{2}$. It forms a statistical manifold of dimension $(n+1)$ with coordinates $\theta=\left(\mu^{1} \ldots \mu^{n}, \sigma^{2}\right)$.

$$
p(i \mid F)=\frac{1}{Z} e^{-\lambda_{k} f_{i}^{k}}
$$

---

<!-- Page 43 -->

This is derived from maximizing the Shanon entropy $S[p]$ subject to the constraints on the expected values of $n$ functions $f_{i}^{k}=f^{k}\left(x_{i}\right)$, labeled by superscript $k=1,2 \ldots n$.

$$
\left\langle f^{k}\right\rangle=\sum_{i} p_{i} f_{i}^{k}=F_{i}^{k}
$$

They form an $n$-dimensional statistical manifold. The coordinates can be either $F=\left(F^{1} \ldots F^{n}\right)$ or Lagrange multipliers $\lambda=\left(\lambda^{1} \ldots \lambda^{n}\right)$.

The basic idea of differential geometry is that curved spaces are locally flat: curvature effects can be ignored in a locally flat region. Hence, within the close vicinity of any point $x$ we can always transform from the original coordinate $x^{a}$ to new coordinates $\hat{x}^{\alpha}=\hat{x}^{\alpha}\left(x^{1} \ldots x^{n}\right)$ that we declare to be local Cartesian (here denoted by $\hat{x}^{\alpha}$ ). An infinitesimal displacement is given by

$$
\mathrm{d} \hat{x}^{\alpha}=X_{a}^{\alpha} \mathrm{d} x^{a} \quad \text { where } \quad X_{a}^{\alpha}=\frac{\partial \hat{x}^{\alpha}}{\partial x^{a}}
$$

and the corresponding infinitesimal distance is

$$
\mathrm{d} l^{2}=\delta_{\alpha \beta} \mathrm{d} \hat{x}^{\alpha} \mathrm{d} \hat{x}^{\beta}
$$

Changing back to the original frame

$$
\mathrm{d} l^{2}=\delta_{\alpha \beta} X_{a}^{\alpha} X_{a}^{\beta} \mathrm{d} x^{a} \mathrm{~d} x^{b}
$$

defining the quantity $g_{a b}=\delta_{\alpha \beta} X_{a}^{\alpha} X_{a}^{\beta}$, we can retrieve the formula 3.47. Under a coordinate transformation:

$$
g_{a b}=X_{a}^{a^{\prime}} X_{b}^{b^{\prime}} g_{a^{\prime} b^{\prime}} \quad \text { where } \quad X_{a}^{a^{\prime}}=\frac{\partial x^{a^{\prime}}}{\partial x^{a}}
$$

so that infinitesimal distance $\mathrm{d} l$ is invariant under the transformation of coordinates. To find the length of a finite length curve $x(\lambda)$, we integrate:

$$
l=\int_{\lambda_{1}}^{\lambda_{2}} \mathrm{~d} l=\int_{\lambda_{1}}^{\lambda_{2}}\left(g_{a b} \frac{d x^{a}}{d \lambda} \frac{d x^{b}}{d \lambda}\right)^{1 / 2} \mathrm{~d} \lambda
$$

With a definition of distance, we can also measure other quanities like angle, areas and volumes. To find an expression for $n$-dimensional volume, we transform to the Cartesian frame

$$
\mathrm{d} V_{n}=\mathrm{d} \hat{x}_{1} \mathrm{~d} \hat{x}_{2} \ldots \mathrm{~d} \hat{x}_{n}
$$

and then transform back to the original coordinates:

$$
\mathrm{d} V_{n}=\left|\frac{\mathrm{d} \hat{x}}{\mathrm{~d} x}\right| \mathrm{d} x_{1} \mathrm{~d} x_{2} \ldots \mathrm{~d} x_{n}=\left|\operatorname{det} X_{a}^{\alpha}\right| \mathrm{d}^{n} x
$$

to calculate the Jacobian of transformation $\operatorname{det}\left(X_{a}^{\alpha}\right)$, we use the definition of metric $g_{a b}$ again. Taking the determinant, we have

$$
\operatorname{det}\left(g_{a b}\right)=\operatorname{det}\left(\delta_{\alpha \beta} X_{a}^{\alpha} X_{a}^{\beta}\right)=\left[\operatorname{det}\left(X_{a}^{\alpha}\right)\right]^{2}
$$

let $g=\operatorname{det}\left(g_{a b}\right)$. We have an expression for the volume element with respect to the metric $g_{a b}(x)$ :

$$
\mathrm{d} V_{n}=g^{1 / 2}(x) \mathrm{d}^{n} x
$$

If we have equal probabilities for equal volumes, the following is true:

$$
p(x) \mathrm{d}^{n} x \propto g^{1 / 2}(x) \mathrm{d}^{n} x
$$

Our goal is to derive $g_{a b}$ corresponding to $p(x \mid \theta)$ in a way that illuminates the meaning of information metric, its interpretation and how it is used. We give two such derivations.

---

<!-- Page 44 -->

# 3.4.1 Derivations from Distinguishability 

Consider the difference between two probability distributions

$$
\Delta=\frac{p(x \mid \theta+\mathrm{d} \theta)-p(x \mid \theta)}{p(x \mid \theta)}=\frac{\partial \log p(x \mid \theta)}{\partial \theta^{a}} \mathrm{~d} \theta^{a}
$$

One of the candidate can be the expectation of these relative differences $\langle\Delta\rangle$. However, they would vanish under averaging. The variance does not vanish from averaging:

$$
\mathrm{d} l^{2}=\left\langle\Delta^{2}\right\rangle=\int \mathrm{d} x p(x \mid \theta) \frac{\partial \log p(x \mid \theta)}{\partial \theta^{a}} \frac{\partial \log p(x \mid \theta)}{\partial \theta^{b}} \mathrm{~d} \theta^{a} \mathrm{~d} \theta^{b}
$$

this is the metric that we seek. A small distance means that the difference between probability distribution is small and therefore indistinguishable. It suggests introducing the metric $g_{a b}$

$$
g_{a b}=\int \mathrm{d} x p(x \mid \theta) \frac{\partial \log p(x \mid \theta)}{\partial \theta^{a}} \frac{\partial \log p(x \mid \theta)}{\partial \theta^{b}}
$$

called the Fisher information matrix. To introduce the notion of distance, normally we say two points are difficult to distinguish because they are close. We can invert this by saying two point $\theta$ and $\theta+\mathrm{d} \theta$ are close whenever they are hard to distinguish. Furthermore, the variance $\mathrm{d} l^{2}=\left\langle\Delta^{2}\right\rangle$ is always positive and vanishes only if $\mathrm{d} \theta$ is zero. Thus, it is natural to introduce $g_{a b}$ as the metric tensor of a Riemannian Space. This is the information metric.

### 3.4.2 Derivation from entropy

Here, we use relative entropy $S[p, q]$ as a tool for updating probabilities from a prior $q$ to posteior $p$ when new information in the forms of constriant becomes available. Therefore, $S[p, q]$ can be used to rank those distributions $p$ relative to $q$ so that the preferred posterior is that maximizes $S[p, q]$ under constraints. The functional form of $S[p, q]$ is derived from the very conservative design that recognize the value of information: what has been learned in the past is valuable and should not be discarded unless rendered obsolete by new information. This is the Principle of Minimal Updating: beliefs should be revised only to the extent required by the new evidence. From this, those distributions $p$ that has a higher entropy is closer to $q$ because they reflect a less drastic revisions to our beliefs.

The term closer means there is a connection between entropy and distance. However, since entropy is not reflective $S[p, q] \neq S[q, p]$, it cannot be used as distance directly. Since length is a local concept whereas entropy is a non-local concept, if there is any relation between distance and entropy, it should be a relation between two infinitesimally distribution $q$ and $p=q+\mathrm{d} q$.

Now consider the entropy of one distribution $p(x \mid \theta^{\prime})$ relative to another $p(x \mid \theta)$

$$
S\left(\theta, \theta^{\prime}\right)=-\int \mathrm{d} x p(x \mid \theta^{\prime}) \log \frac{p(x \mid \theta^{\prime})}{p(x \mid \theta)}
$$

we want to Taylor expand this entropy to see how it changes when $\theta^{\prime}=\theta+\mathrm{d} \theta$ is in the close vicinity of $\theta$. Using Gibb's inequality, the first non-vanishing term in the Taylor expansion is in second order:

$$
S(\theta+\mathrm{d} \theta, \theta)=\left.\frac{1}{2} \frac{\partial^{2} S\left(\theta, \theta^{\prime}\right)}{\partial \theta^{\prime a} \partial \theta^{\prime b}}\right|_{\theta^{\prime}=\theta} \mathrm{d} \theta^{a} \mathrm{~d} \theta^{b}+\cdots \leqslant 0
$$

this suggests distance can be defined as

$$
S(\theta+\mathrm{d} \theta, \theta)=-\frac{1}{2} \mathrm{~d} l^{2}
$$

which means a calculation of the second order derivative gives the information metric

$$
\left.\frac{\partial^{2} S\left(\theta, \theta^{\prime}\right)}{\partial \theta^{\prime a} \partial \theta^{\prime b}}\right|_{\theta^{\prime}=\theta}=\int \mathrm{d} x p(x \mid \theta) \frac{\partial \log p(x \mid \theta)}{\partial \theta^{a}} \frac{\partial \log p(x \mid \theta)}{\partial \theta^{b}}=g_{a b}
$$

---

<!-- Page 45 -->

# 3.4.3 Futher Works 

The design of information geometry and its application on LBM is in the frontier of the research field. To see how information geometry hinted us on a new relaxation scheme, we look at the diagram below
![img-1.jpeg](img-1.jpeg)

Figure 3.1.
In the usual entropy scheme that is discussed in detail before, we have done the following on the information manifold: choosing a direction that points to an equal entropy state $f^{*}$, this direction is characterized by parameter $\alpha$. Since diffusivity/viscosity prevents us to remain in an equal entropy state, we have to choose to advance for a particular distance in this direction. This is characterized by parameter $\beta$. Hence, we have the total relaxation/collision step as

$$
f_{i}^{n+1}=f_{i}+\alpha \beta\left(f_{i}^{n}-f_{i}^{\text {eq }}\right)
$$

From the least action principle and the maximum entropy principle (second law of thermodynamics), we know that in actual physics, population always move in a direction that maximizes the entropy (or equivalently, minimizes the $H$ function). Therefore, on information manifold, now a new relaxation procedure would tell us to move the population in the direction that produces the maximal entropy change and finally reach the equilibrium distribution. This direction is always orthogonal to the equal $H$-function or entropy curve above.

Since this particular normal curve in desire that is perpendicular to the contour line of equal entropy and connects population $f$ and $f^{\text {eq }}$ cannot be a straight line on the information manifold, finding a particular direction would be hard. Therefore, it is now possible to perform a coordinate transformation such that in the new geometry, the equal entropy line would be straight. That will also striaghten the normal curve, making the direction easy to find in the transformed coordinate system. After finding this particular direction, a reverse transformation of the coordinate system allows us to perform relaxation on the original infomation manifold. This direction is the information geodesic. Hence, the new entropic update would be an information geodesic update. The next step of the post-honor thesis research would be realizing this new entropic scheme.

### 3.5 Numerical Results

To prove the validity of the entropic scheme, we perform numerical simulation of a flow in a rectangular pipe. In this geometry, a initial pressure drop will allow fluid to flow from left side of the pipe to the right side of the pipe. The pressure profile will gradually flatten if viscosity is

![img-1.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page45-img-1.jpeg.jpeg)



---

<!-- Page 46 -->

prominent in the fluid. On the other hand, if viscosity is low or the Reynold's number is high, the fluid would not loses its drop on its pressure profile if it is stable physically. Fluid dynamics states that physically the fluid is stable for Reynolds number about 2300 . We will see that the traditional LBGK is not able to simulate this stable fluid flow for Reynolds number about $10^{2}$ or greater. On the other hand, ELBM is unconditionally stable.

1. For initial condition that allows physical simulation of both LBGK and ELBM: $\nu=0.1 \mathrm{~m}^{2} / \mathrm{s}$ $\Delta p=5.0 \mathrm{pa} \operatorname{Re} \sim 10$

- Centerline Pressure Profile (LBGK)
![img-2.jpeg](img-2.jpeg)

Figure 3.2.

- Centerline Pressure Profile (ELBM) (Since we did not implement a correct set of boundary condition exchange for ELBM, the boundary behavior in ELBM is not physical. In the interior of the pipe, however, ELBM shows equivalency to the LBGK.)
![img-3.jpeg](img-3.jpeg)

Figure 3.3.

- Centerplane Pressure Profile after 1250 time steps (LBGK)
![img-4.jpeg](img-4.jpeg)

Figure 3.4.

- Centerplane Pressure Profile after 1250 time steps (ELBM) (Here the incorrect boundary condition for ELBM made perturbation on pressure more prominent than LBGK)

![img-2.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page46-img-2.jpeg.jpeg)



![img-3.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page46-img-3.jpeg.jpeg)



![img-4.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page46-img-4.jpeg.jpeg)



---

<!-- Page 47 -->

![img-5.jpeg](img-5.jpeg)

Figure 3.5.

- Centerplane Pressure Profile after 2500 time steps (LBGK)
![img-6.jpeg](img-6.jpeg)

Figure 3.6.

- Centerplane Pressure Profile after 2500 time steps (ELBM)
![img-7.jpeg](img-7.jpeg)

Figure 3.7.
2. Start of breakdown of LBGK: $\nu=0.01 \mathrm{~m}^{2} / \mathrm{s} \Delta p=5.0 \mathrm{pa} \operatorname{Re} \sim 100$

- Centerline Pressure Profile (LBGK) (Here we see even the viscosity is small, the pressure profile is flattened because of the internal diffusivity of the numerical simulation. This is the feature of numerical instability in LBGK. It prevents a correct reporduction of the fluid dynamical behavior of the fluid.)
![img-8.jpeg](img-8.jpeg)

Figure 3.8.

![img-5.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page47-img-5.jpeg.jpeg)



![img-6.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page47-img-6.jpeg.jpeg)



![img-7.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page47-img-7.jpeg.jpeg)



![img-8.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page47-img-8.jpeg.jpeg)



---

<!-- Page 48 -->

- Centerline Pressure Profile (ELBM) (ELBM is, on the other hand, irresponsive to the change to a smaller viscosity. Diffusion is minimal in the simulation)
![img-9.jpeg](img-9.jpeg)

Figure 3.9.

- Centerplane Pressure Profile after 1250 time steps (LBGK) (LBGK simulation is not physical now)
![img-10.jpeg](img-10.jpeg)

Figure 3.10.

- Centerplane Pressure Profile after 1250 time steps (ELBM) (ELBM simulations retains a reasonable fluid behavior.)
![img-11.jpeg](img-11.jpeg)

Figure 3.11.

- Centerplane Pressure Profile after 2500 time steps (LBGK)
![img-12.jpeg](img-12.jpeg)

Figure 3.12.

![img-9.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page48-img-9.jpeg.jpeg)



![img-10.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page48-img-10.jpeg.jpeg)



![img-11.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page48-img-11.jpeg.jpeg)



![img-12.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page48-img-12.jpeg.jpeg)



---

<!-- Page 49 -->

- Centerplane Pressure Profile after 2500 time steps (ELBM)
![img-13.jpeg](img-13.jpeg)

Figure 3.13.
3. Total Breakdonw of LBGK, Stable Simulation for ELBM: $\nu=0.001 \mathrm{~m}^{2} / \mathrm{s} \Delta p=5.0 \mathrm{pa} \operatorname{Re} \sim$ 1000

- Centerline Pressure Profile (LBGK)
![img-14.jpeg](img-14.jpeg)

Figure 3.14.

- Centerline Pressure Profile (ELBM)
![img-15.jpeg](img-15.jpeg)

Figure 3.15.

- Centerplane Pressure Profile after 1250 time steps (LBGK)
![img-16.jpeg](img-16.jpeg)

Figure 3.16.

![img-13.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page49-img-13.jpeg.jpeg)



![img-14.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page49-img-14.jpeg.jpeg)



![img-15.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page49-img-15.jpeg.jpeg)



![img-16.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page49-img-16.jpeg.jpeg)



---

<!-- Page 50 -->

- Centerplane Pressure Profile after 1250 time steps (ELBM)
![img-17.jpeg](img-17.jpeg)

Figure 3.17.

- Centerplane Pressure Profile after 2500 time steps (LBGK)
![img-18.jpeg](img-18.jpeg)

Figure 3.18.

- Centerplane Pressure Profile after 2500 time steps (ELBM)
![img-19.jpeg](img-19.jpeg)

Figure 3.19.

![img-17.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page50-img-17.jpeg.jpeg)



![img-18.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page50-img-18.jpeg.jpeg)



![img-19.jpeg](/Users/sarthakmishra/Documents/Repos/BTP_FINAL/scripts/ocr_output/images/page50-img-19.jpeg.jpeg)



---

<!-- Page 51 -->

.

---

<!-- Page 52 -->

# Bibliography 

[1] S Ansumali, I. V Karlin, and H. C Öttinger. Minimal entropic kinetic models for hydrodynamics. Europhysics Letters (EPL), 63(6):798-804, sep 2003.
[2] Santosh Ansumali and Iliya V. Karlin. Stabilization of the lattice Boltzmann method by the $\boldsymbol{h}$ theorem: a numerical test. Physical Review E, 62(6):7999-8003, dec 2000.
[3] Santosh Ansumali and Iliya V. Karlin. Single relaxation time model for entropic lattice Boltzmann methods. Physical Review E, 65(5):56312, may 2002.
[4] Ariel Caticha. The basics of information geometry. Pages 15-26. Clos Lucé, Amboise, France, 2015.
[5] Uriel Frisch. Lattice gas automata for the Navier-Stokes equations. a new approach to hydrodynamics and turbulence. Physica Scripta, 40(3):423-423, sep 1989.
[6] I. V Karlin, A Ferrante, and H. C Öttinger. Perfect entropy functions of the Lattice Boltzmann method. Europhysics Letters (EPL), 47(2):182-188, jul 1999.
[7] A.M. Kogan. Derivation of Grad's type equations and study of their relaxation properties by the method of maximization of entropy. Journal of Applied Mathematics and Mechanics, 29(1):130-142, jan 1965.
[8] Timm Krüger, Halim Kusumaatmaja, Alexandr Kuzmin, Orest Shardt, Goncalo Silva, and Erlend Magnus Viggen. The Lattice Boltzmann Method: Principles and Practice. Graduate Texts in Physics. Springer International Publishing, Cham, 2017.
[9] Pierre Lallemand and Li-Shi Luo. Theory of the lattice Boltzmann method: Dispersion, dissipation, isotropy, Galilean invariance, and stability. Physical Review E, 61(6):6546-6562, jun 2000.
[10] Richard L. Liboff. Kinetic theory: classical, quantum, and relativistic descriptions. Graduate texts in contemporary physics. Springer, New York, 3rd ed edition, 2003.
[11] Y. H Qian, D D'Humières, and P Lallemand. Lattice BGK Models for Navier-Stokes Equation. Europhysics Letters (EPL), 17(6):479-484, feb 1992.
[12] S. Succi. The lattice Boltzmann equation for fluid dynamics and beyond. Numerical mathematics and scientific computation. Clarendon Press : Oxford University Press, Oxford : New York, 2001.
[13] Francesca Tosi and Sauro Succi. An Introduction to Entropic Lattice Boltzmann Scheme. SIMAI e-Lecture Notes, 1(0), may 2008. Number: 0.