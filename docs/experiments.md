# OrbitSim Numerical Experiments and Validation

This document outlines the scientific methodology for running numerical experiments and validating orbital integration schemes within OrbitSim.

---

## 1. Scientific Metrics and Invariants

In an unperturbed two-body Keplerian orbit, several physical quantities are analytical constants of motion. Tracking the deviation of these invariants provides a quantitative measure of numerical solver performance:

### Specific Orbital Energy ($\varepsilon$)

$$\varepsilon = \frac{v^2}{2} - \frac{\mu}{r} = \text{const}$$

* **Energy Drift**: $\Delta \varepsilon(t) = \varepsilon(t) - \varepsilon(0)$.
* **Relative Energy Error**: $\frac{|\varepsilon(t) - \varepsilon(0)|}{|\varepsilon(0)|}$.

### Specific Angular Momentum ($h$)

$$\mathbf{h} = \mathbf{r} \times \mathbf{v} \implies h_z = x v_y - y v_x = \text{const}$$

* **Angular Momentum Drift**: $\Delta h_z(t) = h_z(t) - h_z(0)$.

### Orbital Geometry

* **Radius Drift**: For a circular orbit, $r(t) - r(0)$ should ideally be $0$.
* **Eccentricity Vector**: $\mathbf{e} = \frac{\mathbf{v} \times \mathbf{h}}{\mu} - \frac{\mathbf{r}}{r}$.

---

## 2. Methodology for Integrator Benchmarking

Every experiment will follow a standardized protocol:

1. **Parameter Sweep**:
   - Vary timestep $\Delta t \in \{60.0\text{ s}, 30.0\text{ s}, 15.0\text{ s}, 7.5\text{ s}, 3.75\text{ s}, 1.875\text{ s}\}$.
2. **Duration**:
   - Short term (1 orbit $\approx 5553.62\text{ s}$).
   - Medium term (10 orbits).
   - Long term (100+ orbits).
3. **Outputs**:
   - Tabular trajectory and summary metrics saved in `results/`.
   - Convergence verification against analytical solutions.

---

## 3. M1 — Explicit Euler Convergence Experiment

### Objective
Empirically investigate the convergence rate and error scaling of the Explicit Euler integrator as a function of timestep $\Delta t$ on an analytical 2D circular Low Earth Orbit (LEO).

### Mathematical Reference Solution
For an unperturbed circular orbit with radius $r_0 = R_E + 400\text{ km} = 6778.137\text{ km}$ and Earth gravitational parameter $\mu = 398600.435\text{ km}^3/\text{s}^2$:

* Mean motion:
  $$n = \sqrt{\frac{\mu}{r_0^3}}$$
* Orbital period:
  $$T = \frac{2\pi}{n} \approx 5553.624319\text{ s} \quad (\approx 92.56\text{ min})$$
* Exact state vector at time $t \ge 0$:
  $$x_{\text{exact}}(t) = r_0 \cos(nt)$$
  $$y_{\text{exact}}(t) = r_0 \sin(nt)$$
  $$v_{x, \text{exact}}(t) = -r_0 n \sin(nt)$$
  $$v_{y, \text{exact}}(t) = r_0 n \cos(nt)$$

### Experiment Parameters and Interval
* **Initial State**: $x_0 = r_0$, $y_0 = 0$, $v_{x0} = 0$, $v_{y0} = \sqrt{\mu / r_0}$.
* **Time Interval**: Exactly one orbital period $[0, T]$ ($t_{\text{final}} = T$).
* **Timesteps**:
  $$\Delta t \in \{60.0\text{ s}, 30.0\text{ s}, 15.0\text{ s}, 7.5\text{ s}, 3.75\text{ s}, 1.875\text{ s}\}$$
* **Step Shortening**: Because $T$ is generally not an integer multiple of $\Delta t$, the final step is shortened:
  $$\Delta t_{\text{step}} = \min(\Delta t, T - t)$$
  guaranteeing that the simulation terminates precisely at $t = T$.

### Error Definitions
Position and velocity are evaluated separately due to differing physical dimensions:

* **Position Error**:
  $$e_r(t) = \sqrt{(x_{\text{num}}(t) - x_{\text{exact}}(t))^2 + (y_{\text{num}}(t) - y_{\text{exact}}(t))^2} \quad [\text{km}]$$
* **Velocity Error**:
  $$e_v(t) = \sqrt{(v_{x, \text{num}}(t) - v_{x, \text{exact}}(t))^2 + (v_{y, \text{num}}(t) - v_{y, \text{exact}}(t))^2} \quad [\text{km/s}]$$

For each run, we record:
* Final errors: $e_r(T)$, $e_v(T)$.
* Maximum errors over the period: $\max_{t \in [0, T]} e_r(t)$, $\max_{t \in [0, T]} e_v(t)$.

### Convergence Order Formulation
Between two consecutive timesteps $\Delta t_1 > \Delta t_2$, empirical convergence order $p$ is calculated as:

$$p = \frac{\ln(E(\Delta t_1) / E(\Delta t_2))}{\ln(\Delta t_1 / \Delta t_2)}$$

evaluated separately for position error and velocity error.

### Physical Invariant Diagnostics
* **Specific Orbital Energy**:
  $$\varepsilon(t) = \frac{v(t)^2}{2} - \frac{\mu}{r(t)}, \quad \varepsilon_0 = -\frac{\mu}{2 r_0}$$
  $$\text{Relative Energy Error} = \frac{\varepsilon(t) - \varepsilon_0}{|\varepsilon_0|}$$
* **Specific Angular Momentum**:
  $$h_z(t) = x(t) v_y(t) - y(t) v_x(t), \quad h_0 = r_0 \sqrt{\frac{\mu}{r_0}}$$
  $$\text{Relative Angular Momentum Error} = \frac{h_z(t) - h_0}{|h_0|}$$

### Scientific Distinction: Hypothesis vs Direct Measurement
* **Direct Measurements**:
  - Exact Euclidean position errors $e_r(T)$ and $\max e_r(t)$ [km].
  - Exact Euclidean velocity errors $e_v(T)$ and $\max e_v(t)$ [km/s].
  - Fractional drift in specific energy and angular momentum.
  - Empirical slope $p$ between discrete runs.
* **Scientific Hypotheses / Theoretical Expectations**:
  - Explicit Euler is theoretically first order ($p = 1$); empirical $p$ is expected to asymptotically approach $1$ as $\Delta t \to 0$.
  - Explicit Euler is non-symplectic and known to systematically inflate orbital energy in central force fields over time.
