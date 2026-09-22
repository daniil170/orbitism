# Numerical Methods: Explicit Euler Integration

This document describes the numerical integration methodology implemented in OrbitSim (Milestone M0).

---

## 1. Overview and Rationale

To solve the system of ordinary differential equations (ODEs):

$$\frac{d\mathbf{s}}{dt} = \mathbf{f}(t, \mathbf{s})$$

we implement the **Explicit Euler method** (forward Euler). It serves as the baseline first-order integrator for our numerical benchmarking pipeline.

---

## 2. Derivation from Taylor Series

Expanding the state vector $\mathbf{s}(t + \Delta t)$ around $t$ using Taylor series:

$$\mathbf{s}(t + \Delta t) = \mathbf{s}(t) + \frac{d\mathbf{s}}{dt}\Delta t + \mathcal{O}(\Delta t^2)$$

Truncating after the first derivative term gives the discrete update rule:

$$\mathbf{s}_{n+1} = \mathbf{s}_n + \mathbf{f}(t_n, \mathbf{s}_n) \cdot \Delta t$$

---

## 3. Explicit Update Equations

For our 2D orbital state $\mathbf{s} = [x, y, v_x, v_y]^T$, the component-wise updates are:

$$x_{n+1} = x_n + v_{x, n} \cdot \Delta t$$

$$y_{n+1} = y_n + v_{y, n} \cdot \Delta t$$

$$v_{x, n+1} = v_{x, n} + a_{x, n} \cdot \Delta t$$

$$v_{y, n+1} = v_{y, n} + a_{y, n} \cdot \Delta t$$

where accelerations are evaluated at the beginning of the step:

$$a_{x, n} = -\mu \frac{x_n}{r_n^3}, \quad a_{y, n} = -\mu \frac{y_n}{r_n^3}$$

---

## 4. Error Characteristics

* **Local Truncation Error (LTE)**: $\mathcal{O}(\Delta t^2)$ per step.
* **Global Truncation Error (GTE)**: $\mathcal{O}(\Delta t)$ over a fixed time interval $[0, T]$.
* **Order of Accuracy**: First order ($p = 1$).

### Physical Implications for Orbital Mechanics

The Explicit Euler method is **non-symplectic** and **non-conservative**. For a central gravitational field, it systematically introduces positive energy error at each step, causing the simulated satellite to artificially gain energy and spiral outwards over time.

This known limitation makes Explicit Euler an ideal benchmark to contrast against higher-order methods (e.g. RK4) and symplectic integrators (e.g. Verlet/Leapfrog) in future project phases.

---

## 5. Architectural Decoupling

The integrator interface in `src/integrators/euler.py` is strictly decoupled from orbital physics:
* The integrator accepts any callable `derivative_fn(t, state) -> dstate/dt`.
* It does not contain domain-specific knowledge of gravity, satellite masses, or orbital radii.
* Any autonomous or non-autonomous ODE system conforming to this signature can be stepped with this integrator.

