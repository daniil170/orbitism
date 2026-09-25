# Numerical Methods: Integration Schemes

This document describes the numerical integration methodologies implemented in OrbitSim:
1. **Explicit Euler** (Milestone M0)
2. **Classical 4th-Order Runge-Kutta (RK4)** (Milestone M5)

---

## 1. Overview and Rationale

To solve the system of ordinary differential equations (ODEs):

$$\frac{d\mathbf{s}}{dt} = \mathbf{f}(t, \mathbf{s})$$

we implement numerical integration methods that step the system forward in time. Explicit Euler serves as our baseline first-order integrator, while classical RK4 provides a high-order benchmark.

---

## 2. Derivation from Taylor Series (Explicit Euler)

Expanding the state vector $\mathbf{s}(t + \Delta t)$ around $t$ using Taylor series:

$$\mathbf{s}(t + \Delta t) = \mathbf{s}(t) + \frac{d\mathbf{s}}{dt}\Delta t + \mathcal{O}(\Delta t^2)$$

Truncating after the first derivative term gives the discrete update rule:

$$\mathbf{s}_{n+1} = \mathbf{s}_n + \mathbf{f}(t_n, \mathbf{s}_n) \cdot \Delta t$$

---

## 3. Explicit Update Equations (Explicit Euler)

For our 2D orbital state $\mathbf{s} = [x, y, v_x, v_y]^T$, the component-wise updates are:

$$x_{n+1} = x_n + v_{x, n} \cdot \Delta t$$

$$y_{n+1} = y_n + v_{y, n} \cdot \Delta t$$

$$v_{x, n+1} = v_{x, n} + a_{x, n} \cdot \Delta t$$

$$v_{y, n+1} = v_{y, n} + a_{y, n} \cdot \Delta t$$

where accelerations are evaluated at the beginning of the step:

$$a_{x, n} = -\mu \frac{x_n}{r_n^3}, \quad a_{y, n} = -\mu \frac{y_n}{r_n^3}$$

---

## 4. Error Characteristics (Explicit Euler)

* **Local Truncation Error (LTE)**: $\mathcal{O}(\Delta t^2)$ per step.
* **Global Truncation Error (GTE)**: $\mathcal{O}(\Delta t)$ over a fixed time interval $[0, T]$.
* **Order of Accuracy**: First order ($p = 1$).

### Physical Implications for Orbital Mechanics

The Explicit Euler method is **non-symplectic** and **non-conservative**. For a central gravitational field, it systematically introduces positive energy error at each step, causing the simulated satellite to artificially gain energy and spiral outwards over time.

---

## 5. Architectural Decoupling

Both integrator interfaces (`src/integrators/euler.py` and `src/integrators/rk4.py`) are strictly decoupled from orbital physics:
* Integrators accept any callable `derivative_fn(t, state) -> dstate/dt`.
* They contain no domain-specific knowledge of gravity, satellite masses, or orbital radii.
* Any autonomous or non-autonomous ODE system conforming to this signature can be stepped generically.

---

## 6. Classical 4th-Order Runge-Kutta Method (RK4)

### Formulation and Butcher Tableau

The classical 4th-order Runge-Kutta method (often referred to as Kutta's method or simply RK4) evaluates the derivative vector $\mathbf{f}(t, \mathbf{s})$ four times per step to cancel Taylor series error terms up to 4th order.

The Butcher tableau for classical RK4 is:

$$\begin{array}{c|cccc}
0 & & & & \\
1/2 & 1/2 & & & \\
1/2 & 0 & 1/2 & & \\
1 & 0 & 0 & 1 & \\
\hline
& 1/6 & 1/3 & 1/3 & 1/6
\end{array}$$

### Stage Equations

For state vector $\mathbf{s}_n$ at time $t_n$ and step size $\Delta t$:

1. **Stage 1 (Initial Slope)**:
   $$\mathbf{k}_1 = \mathbf{f}(t_n, \mathbf{s}_n)$$

2. **Stage 2 (Midpoint Slope from $\mathbf{k}_1$)**:
   $$\mathbf{k}_2 = \mathbf{f}\left(t_n + \frac{\Delta t}{2}, \mathbf{s}_n + \frac{\Delta t}{2} \mathbf{k}_1\right)$$

3. **Stage 3 (Midpoint Slope from $\mathbf{k}_2$)**:
   $$\mathbf{k}_3 = \mathbf{f}\left(t_n + \frac{\Delta t}{2}, \mathbf{s}_n + \frac{\Delta t}{2} \mathbf{k}_2\right)$$

4. **Stage 4 (Endpoint Slope from $\mathbf{k}_3$)**:
   $$\mathbf{k}_4 = \mathbf{f}\left(t_n + \Delta t, \mathbf{s}_n + \Delta t\,\mathbf{k}_3\right)$$

5. **State Update (Simpson's Quadrature Weighted Average)**:
   $$\mathbf{s}_{n+1} = \mathbf{s}_n + \frac{\Delta t}{6}\left(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4\right)$$

### Error Characteristics and Convergence Order

* **Local Truncation Error (LTE)**: $\mathcal{O}(\Delta t^5)$ per step.
* **Global Truncation Error (GTE)**: $\mathcal{O}(\Delta t^4)$ over a fixed duration $[0, T]$.
* **Order of Accuracy**: Fourth order ($p = 4$).

### Physical Implications for Orbital Mechanics

* **Higher-Order Truncation Suppression**: By canceling error terms through 4th order, RK4 suppresses error accumulation by multiple orders of magnitude relative to Explicit Euler ($> 10^4\times$ to $10^9\times$ depending on $\Delta t$).
* **Non-Symplectic Nature**: Like Explicit Euler, classical RK4 is non-symplectic and does not strictly preserve the Hamiltonian or phase-space symplectic 2-form ($d\mathbf{p} \wedge d\mathbf{q} \neq \text{const}$). However, unlike Euler's large positive energy injection, RK4 introduces a tiny negative energy dissipation ($\delta\varepsilon < 0$), causing an exceedingly slow inward spiral.
* **Absence of Phase-Winding**: Over multi-orbit spans ($100T$), the energy error remains so small ($|\delta\varepsilon| \le 2.5 \times 10^{-5}$) that the semi-major axis drift $|\Delta a| < 0.2\text{ km}$, keeping the cumulative phase lag $|\Delta\theta| \ll 2\pi$. As a result, RK4 never enters the geometric phase-wound regime that caused non-monotonic Cartesian error in Explicit Euler.

---

## 7. Velocity Verlet Integrator (Milestone M6)

### Formulation and Algorithm

The Velocity Verlet method is a canonical, time-reversible, second-order **symplectic integrator** designed specifically for separable Hamiltonian systems $\mathcal{H}(\mathbf{q}, \mathbf{p}) = \mathcal{T}(\mathbf{p}) + \mathcal{V}(\mathbf{q})$.

The discrete step equations advancing state from $t_n$ to $t_{n+1} = t_n + \Delta t$ are:

1. **Position Update**:
   $$\mathbf{r}_{n+1} = \mathbf{r}_n + \mathbf{v}_n \Delta t + \frac{1}{2} \mathbf{a}(\mathbf{r}_n) \Delta t^2$$

2. **Acceleration Evaluation**:
   $$\mathbf{a}_{n+1} = \mathbf{a}(\mathbf{r}_{n+1}) = -\mu \frac{\mathbf{r}_{n+1}}{r_{n+1}^3}$$

3. **Velocity Update**:
   $$\mathbf{v}_{n+1} = \mathbf{v}_n + \frac{1}{2}\left[\mathbf{a}(\mathbf{r}_n) + \mathbf{a}_{n+1}\right] \Delta t$$

### Error Characteristics and Symplectic Invariance

* **Local Truncation Error (LTE)**: $\mathcal{O}(\Delta t^3)$ for position, $\mathcal{O}(\Delta t^3)$ for velocity.
* **Global Truncation Error (GTE)**: $\mathcal{O}(\Delta t^2)$ over fixed duration $[0, T]$.
* **Order of Accuracy**: Second order ($p = 2$).
* **Symplecticity**: The Jacobian map $J = \partial(\mathbf{r}_{n+1}, \mathbf{v}_{n+1}) / \partial(\mathbf{r}_n, \mathbf{v}_n)$ satisfies $\det(J) \equiv 1$ and preserves the symplectic 2-form $\omega = d\mathbf{r} \wedge d\mathbf{v}$.
* **Shadow Hamiltonian**: By backward error analysis, Velocity Verlet exactly solves a modified Hamiltonian $\widetilde{\mathcal{H}} = \mathcal{H} + \mathcal{O}(\Delta t^2)$. Consequently, orbital energy $\varepsilon$ does not secularly drift, but oscillates within a bounded envelope of width $\mathcal{O}(\Delta t^2)$ indefinitely.
* **Exact Angular Momentum Conservation**: For any central force $\mathbf{a}(\mathbf{r}) \parallel \mathbf{r}$, the torque vanishes identically:
  $$\mathbf{r}_{n+1} \times \mathbf{v}_{n+1} = \mathbf{r}_n \times \mathbf{v}_n$$
  Specific angular momentum $h_z$ is an exact discrete invariant preserved down to floating-point roundoff.

---

## 8. Yoshida 4th-Order Symplectic Integrator (Milestone M7)

### Formulation and Composition Principle

The Yoshida 4th-order integrator (Yoshida, 1990) is constructed via a symmetric composition of a symmetric, time-reversible, second-order symplectic base integrator $S_2(h)$ (such as Velocity Verlet / leapfrog):

$$S_4(h) = S_2(w_1 h) \circ S_2(w_0 h) \circ S_2(w_1 h)$$

Because $S_2(h)$ is symmetric ($S_2(-h)^{-1} = S_2(h)$), its asymptotic expansion in terms of Lie differential operators contains only odd powers of the timestep $h$:

$$S_2(h) = \exp\left(h D + h^3 D_3 + h^5 D_5 + \dots\right)$$

Applying the symmetric composition of three substeps with weights $(w_1, w_0, w_1)$ yields the effective Lie generator:

$$S_4(h) = \exp\left[(2w_1 + w_0) h D + (2w_1^3 + w_0^3) h^3 D_3 + \mathcal{O}(h^5)\right]$$

### Yoshida Coefficients Origin

To achieve 4th-order consistency, the composition weights must satisfy two algebraic conditions:

1. **Consistency Condition** (order 1 recovery):
   $$2 w_1 + w_0 = 1$$

2. **Order Condition** (elimination of the 3rd-order error term):
   $$2 w_1^3 + w_0^3 = 0 \implies w_0 = -2^{1/3} w_1$$

Substituting $w_0$ into the consistency condition yields:

$$w_1 = \frac{1}{2 - 2^{1/3}} \approx 1.3512071919596576$$

$$w_0 = -\frac{2^{1/3}}{2 - 2^{1/3}} \approx -1.7024143839193153$$

### The Negative Substep ($w_0 h$)

The intermediate substep weight $w_0 < 0$ is strictly negative. This is not an error:
* By the Suzuki (1991) theorem, no symmetric composition of order $p \ge 3$ with real coefficients can have all positive weights.
* The negative substep represents an exact backwards-in-time integration step.
* Because the underlying Velocity Verlet scheme is unconditionally time-reversible ($S_2(-h) \circ S_2(h) = I$), stepping backwards by $w_0 h$ is mathematically rigorous, numerically stable, and exact.

### Computational Cost & Force Evaluations

In our decoupled architecture (`src/integrators/yoshida4.py`), the integrator is implemented as a pure composition of three stateless `velocity_verlet_step` calls:
* **Substep 1**: advances by $w_1 h$ (2 force evaluations: $a(r_0)$ and $a(r_1)$).
* **Substep 2**: advances by $w_0 h$ (2 force evaluations: $a(r_1)$ and $a(r_2)$).
* **Substep 3**: advances by $w_1 h$ (2 force evaluations: $a(r_2)$ and $a(r_3)$).
* **Total Force Evaluations**: exactly **6 force evaluations per outer timestep**.

#### Optimization Analysis:
* At intermediate substep boundaries, $a(r_1)$ and $a(r_2)$ are evaluated twice (at the end of one substep and the beginning of the next).
* In autonomous systems $\mathbf{a} = \mathbf{a}(\mathbf{r})$, these evaluations could theoretically be cached to reduce the cost to 4 evaluations per outer step (or 3 across outer steps using FSAL).
* In M7, to preserve stateless architectural purity, avoid stateful side-effects, and strictly verify the mathematical composition $S_2 \circ S_2 \circ S_2$, caching is intentionally avoided (no premature optimization).

### Invariants and Shadow Hamiltonian

* **Local Truncation Error (LTE)**: $\mathcal{O}(h^5)$ per step.
* **Global Truncation Error (GTE)**: $\mathcal{O}(h^4)$ over fixed duration $[0, T]$.
* **Symplecticity & Phase-Space Volume**: As a composition of canonical symplectic mappings, the discrete update $S_4(h)$ strictly preserves the canonical symplectic 2-form $\omega = d\mathbf{p} \wedge d\mathbf{q}$ and phase-space volume ($\det J \equiv 1$). This is an exact mathematical property of the map.
* **Shadow Hamiltonian & Energy**: By backward error analysis, Yoshida 4 exactly tracks the trajectory of a perturbed shadow Hamiltonian $\widetilde{\mathcal{H}} = \mathcal{H} + \mathcal{O}(h^4)$. Consequently, energy error does not experience secular growth, but oscillates within a bounded envelope. On an unperturbed circular orbit, because $r_0$ is a stationary minimum of the effective potential $V_{\text{eff}}(r)$, the energy deviation scales as $(\Delta r)^2 \sim \mathcal{O}(h^8)$ for large steps before reaching the machine precision roundoff floor ($\sim 10^{-14}$).
* **Central Force Angular Momentum**: Specific angular momentum $h_z$ is theoretically an exact invariant for central forces in each symmetric Verlet substep; in numerical computation it is preserved down to double-precision floating-point roundoff ($|\delta h| < 10^{-13}$).
