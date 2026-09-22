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
  - Explicit Euler is non-symplectic and not energy-preserving for central force fields.

---

## 4. M2 — Time Evolution of Explicit Euler Error

### Scientific Question
How does the numerical error of the Explicit Euler integrator evolve over time throughout one complete orbital period ($t \in [0, T]$) for different numerical timesteps $\Delta t$?

Unlike M1 (which investigated the final error $E_{\text{final}} = f(\Delta t)$ at $t = T$), M2 investigates the discrete error trajectory across sampled time steps:
$$E = f(t, \Delta t), \quad t \in [0, T]$$

### Setup & Initial Conditions
* **Physical Model**: 2D classical two-body problem (Earth fixed at origin).
* **Gravitational Parameter**: $\mu = 398600.435\text{ km}^3/\text{s}^2$.
* **Reference Orbit**: Circular Low Earth Orbit (LEO) at altitude $h = 400.0\text{ km}$ above equatorial radius $R_E = 6378.137\text{ km}$.
* **Radius**: $r_0 = R_E + h = 6778.137\text{ km}$.
* **Mean Motion**: $n = \sqrt{\mu / r_0^3} \approx 0.001131358\text{ rad/s}$.
* **Analytical Period**: $T = \frac{2\pi}{n} \approx 5553.624318623783\text{ s} \quad (\approx 92.56\text{ min})$.
* **Initial State**: $x_0 = r_0$, $y_0 = 0\text{ km}$, $v_{x0} = 0\text{ km/s}$, $v_{y0} = \sqrt{\mu / r_0} \approx 7.668558109995\text{ km/s}$.
* **Timesteps**:
  $$\Delta t \in \{60.0\text{ s}, 30.0\text{ s}, 15.0\text{ s}, 7.5\text{ s}, 3.75\text{ s}, 1.875\text{ s}\}$$
* **Step Shortening**: $\Delta t_{\text{step}} = \min(\Delta t, T - t)$ ensures simulation terminates exactly at $t = T$ without overshoot.

### Error and Diagnostic Formulations
At each discrete timestamp $t_k \in [0, T]$, the numerical state $[x_{\text{num}}, y_{\text{num}}, v_{x,\text{num}}, v_{y,\text{num}}]$ is evaluated against the independent exact analytical state:

* **Analytical Reference State**:
  $$x_{\text{exact}}(t) = r_0 \cos(nt), \quad y_{\text{exact}}(t) = r_0 \sin(nt)$$
  $$v_{x, \text{exact}}(t) = -r_0 n \sin(nt), \quad v_{y, \text{exact}}(t) = r_0 n \cos(nt)$$

* **Position Error**:
  $$e_r(t) = \sqrt{(x_{\text{num}}(t) - x_{\text{exact}}(t))^2 + (y_{\text{num}}(t) - y_{\text{exact}}(t))^2} \quad [\text{km}]$$

* **Velocity Error**:
  $$e_v(t) = \sqrt{(v_{x,\text{num}}(t) - v_{x,\text{exact}}(t))^2 + (v_{y,\text{num}}(t) - v_{y,\text{exact}}(t))^2} \quad [\text{km/s}]$$

* **Specific Orbital Energy Diagnostic**:
  $$\varepsilon(t) = \frac{v(t)^2}{2} - \frac{\mu}{r(t)}, \quad \varepsilon_0 = -\frac{\mu}{2 r_0}$$
  $$\delta\varepsilon(t) = \frac{\varepsilon(t) - \varepsilon_0}{|\varepsilon_0|}$$

* **Specific Angular Momentum Diagnostic**:
  $$h_z(t) = x(t) v_y(t) - y(t) v_x(t), \quad h_0 = r_0 \sqrt{\frac{\mu}{r_0}}$$
  $$\delta h(t) = \frac{h_z(t) - h_0}{|h_0|}$$

### Scientific Classification of Results

#### 1. Tested Facts
* Evaluated across 6 discrete timestep values ($\Delta t \in \{60.0, 30.0, 15.0, 7.5, 3.75, 1.875\}$ s) over one orbital period $T = 5553.624318623783$ s.
* Dataset contains exactly 5840 rows without any NaN/Inf values or duplicate timestamps.
* Simulation terminates precisely at $t = T$ without step overshoot ($t \le T$).
* For the investigated circular orbit, one orbital period, and the tested timestep values, the recorded position and velocity errors were numerically non-decreasing at every saved discrete time point.
* Maximum position and velocity errors equal final errors at $t = T$ for all tested timesteps.
* Independent arithmetic verification matches CSV recorded values to machine floating-point precision.

#### 2. Experimental Observations
* Reducing $\Delta t$ systematically reduces final and maximum position, velocity, energy, and angular momentum errors.
* For this configuration and timeframe, relative error in specific orbital energy $\delta\varepsilon(t)$ starts at zero, remains positive for $t > 0$, is monotonically non-decreasing on the sampled discrete time grid, and reaches its maximum at $t = T$.
* For this configuration and timeframe, relative error in specific angular momentum $\delta h(t)$ starts at zero, remains positive for $t > 0$, is monotonically non-decreasing on the sampled discrete time grid, and reaches its maximum at $t = T$.

#### 3. Theoretical Expectations
* Explicit Euler is a first-order integration method; global truncation error is expected to scale as $\mathcal{O}(\Delta t)$ under standard convergence conditions.
* Explicit Euler is non-symplectic and not energy-preserving for Hamiltonian orbital systems; energy conservation is not theoretically guaranteed.

#### 4. Open Questions
* Error evolution and monotonicity behavior over long-term integration horizons ($t \gg T$, multiple orbits).
* Error trajectory behavior for eccentric orbits ($e > 0$) where velocity and distance vary along the trajectory.
* Boundedness or secular growth of energy error over multi-orbit time scales.

### Main Conclusion of M2
For the investigated circular 2D orbit, over one analytical orbital period and for $\Delta t \in \{60.0, 30.0, 15.0, 7.5, 3.75, 1.875\}$ s, the recorded position and velocity errors were numerically non-decreasing at all saved time points and reached their maximum at the final time $T$. Reducing $\Delta t$ systematically reduced the final error. The observed positive growth of energy and angular-momentum relative errors is an experimental result for this configuration and time interval, not a universal property claimed for Explicit Euler.
