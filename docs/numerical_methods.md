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

