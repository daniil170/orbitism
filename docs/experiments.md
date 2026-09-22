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

---

## 5. M3 — Long-Term Stability of Explicit Euler

### Goal
Evaluate the long-term error trajectory and diagnostic invariant drift of the Explicit Euler integration scheme over multi-orbit time scales ($t \le 100 T$).

### Research Question
How do numerical timestep $\Delta t$ and integration duration ($t \in [0, 100 T]$) influence position error $e_r(t)$, velocity error $e_v(t)$, relative energy error $\delta\varepsilon(t)$, and relative angular momentum error $\delta h(t)$ of the Explicit Euler method in a 2D two-body central force field?

### Hypothesis
1. Explicit Euler is a non-symplectic numerical integration method and does not conserve physical energy invariants. Over long integration horizons ($t \gg T$), relative specific orbital energy error $\delta\varepsilon(t)$ and position error $e_r(t)$ are hypothesized to exhibit non-decreasing growth.
2. Decreasing the numerical timestep $\Delta t$ systematically reduces the magnitude of position, velocity, energy, and angular momentum errors accumulated across the $100 T$ integration horizon.

### Setup
* **Physical Model**: 2D classical two-body Keplerian central gravity field (Earth fixed at origin).
* **Equations of Motion**:
  $$\frac{d\mathbf{r}}{dt} = \mathbf{v}, \quad \frac{d\mathbf{v}}{dt} = -\frac{\mu}{r^3} \mathbf{r}$$

### Initial Conditions
* **Gravitational Parameter**: $\mu = 398600.435\text{ km}^3/\text{s}^2$.
* **Reference Orbital Altitude**: $h = 400.0\text{ km}$ above equatorial surface ($R_E = 6378.137\text{ km}$).
* **Initial Orbital Radius**: $r_0 = R_E + h = 6778.137\text{ km}$.
* **Initial State Vector**:
  $$x_0 = r_0 = 6778.137\text{ km}, \quad y_0 = 0.0\text{ km}$$
  $$v_{x0} = 0.0\text{ km/s}, \quad v_{y0} = \sqrt{\frac{\mu}{r_0}} \approx 7.668558109995\text{ km/s}$$
* **Analytical Period**: $T = 5553.624318623783\text{ s}$ ($\approx 92.56\text{ min}$).

### Numerical Method
* **Scheme**: First-order Explicit Euler method:
  $$\mathbf{r}_{k+1} = \mathbf{r}_k + \mathbf{v}_k \Delta t, \quad \mathbf{v}_{k+1} = \mathbf{v}_k + \mathbf{a}(\mathbf{r}_k) \Delta t$$
* **Boundary Shortening**: Adaptive last step $\Delta t_{\text{step}} = \min(\Delta t, t_{\text{final}} - t)$ to guarantee exact termination at target final time horizon $t_{\text{final}} = N \cdot T$.

### Timesteps
Investigated discrete timesteps (identical to M1 and M2):
$$\Delta t \in \{60.0\text{ s}, 30.0\text{ s}, 15.0\text{ s}, 7.5\text{ s}, 3.75\text{ s}, 1.875\text{ s}\}$$

### Duration
* **Maximum Horizon**: $t_{\text{final}} = 100 T \approx 555362.431862\text{ s}$ ($\approx 6.43\text{ days} \approx 154.27\text{ hours}$).
* **Checkpoint Evaluation Horizons**: $1 T, 5 T, 10 T, 20 T, 50 T, 100 T$.

### Metrics
Strictly standard OrbitSim diagnostic suite (identical to M1 and M2):
* **Euclidean Position Error**:
  $$e_r(t) = \sqrt{(x_{\text{num}}(t) - x_{\text{exact}}(t))^2 + (y_{\text{num}}(t) - y_{\text{exact}}(t))^2} \quad [\text{km}]$$
* **Euclidean Velocity Error**:
  $$e_v(t) = \sqrt{(v_{x,\text{num}}(t) - v_{x,\text{exact}}(t))^2 + (v_{y,\text{num}}(t) - v_{y,\text{exact}}(t))^2} \quad [\text{km/s}]$$
* **Relative Specific Orbital Energy Error**:
  $$\delta\varepsilon(t) = \frac{\varepsilon(t) - \varepsilon_0}{|\varepsilon_0|}, \quad \varepsilon(t) = \frac{v(t)^2}{2} - \frac{\mu}{r(t)}, \quad \varepsilon_0 = -\frac{\mu}{2 r_0}$$
* **Relative Specific Angular Momentum Error**:
  $$\delta h(t) = \frac{h_z(t) - h_0}{|h_0|}, \quad h_z(t) = x(t) v_y(t) - y(t) v_x(t), \quad h_0 = r_0 v_{y0}$$

### Reference
Exact analytical circular Keplerian orbit at discrete timestamp $t$:
$$x_{\text{exact}}(t) = r_0 \cos(nt), \quad y_{\text{exact}}(t) = r_0 \sin(nt)$$
$$v_{x,\text{exact}}(t) = -r_0 n \sin(nt), \quad v_{y,\text{exact}}(t) = r_0 n \cos(nt)$$
where mean motion $n = \sqrt{\frac{\mu}{r_0^3}} \approx 0.001131358\text{ rad/s}$.

### Implementation Details
* **Module Location**: `experiments/exp03_euler_long_term.py`
* **Infrastructure Reuse**:
  - `src/simulation/simulator.py`: Reuses `Simulator` class for stepping and adaptive boundary termination (`t_final = 100 * T`).
  - `src/simulation/metrics.py`: Reuses `compute_error_history` for position error, velocity error, relative orbital energy error, and angular momentum error calculation.
  - `src/physics/gravity.py`: Reuses `MU_EARTH`, `R_EARTH`, and `create_circular_orbit_state`.
  - `src/physics/analytical.py`: Reuses `circular_orbit_period`.

### Data Products
1. **Full Trajectory Dataset**: `results/exp03_euler_long_term.csv`
   - **Schema**: `dt_s`, `step`, `t_s`, `x_km`, `y_km`, `vx_km_s`, `vy_km_s`, `position_error_km`, `velocity_error_km_s`, `energy_error`, `angular_momentum_error`
2. **Checkpoint Diagnostic Summary**: `results/exp03_euler_long_term_summary.csv`
   - **Schema**: `dt_s`, `checkpoint_T`, `time_s`, `position_error_km`, `velocity_error_km_s`, `energy_error`, `angular_momentum_error`
   - **Horizons**: $1T, 5T, 10T, 20T, 50T, 100T$

### Testing Description
* **Test Module**: `tests/test_long_term_stability.py`
* **Verified Facts**:
  1. `test_initial_error_at_t0_is_zero`: Initial state error at $t=0$ equals zero within float precision.
  2. `test_long_term_history_starts_at_t0`: Trajectory history starts at $t=0$ for all timesteps.
  3. `test_long_term_history_ends_at_100T_without_overshoot`: Integration terminates exactly at $t=100T$ with no step exceeding $100T + 10^{-12}\text{ s}$.
  4. `test_long_term_history_strictly_monotonic`: Timestamps are strictly monotonically increasing ($t_{k+1} > t_k$).
  5. `test_long_term_no_nan_or_inf_in_states_and_metrics`: States, position/velocity errors, energy errors, and angular momentum errors contain no NaN or Inf values.
  6. `test_long_term_all_checkpoints_exist`: Summary extraction extracts exact entries for all 6 requested checkpoints ($1T, 5T, 10T, 20T, 50T, 100T$) per timestep.
  7. `test_independent_verification_of_diagnostics`: Independent arithmetic evaluation of analytical reference states and invariant errors matches production module output to machine precision.
* **Regression Tests**:
  1. `test_regression_against_m2_trajectory_steps`: Un-shortened trajectory states during the first orbital period match M2 trajectory data to floating point precision.

### Limitations
1. Restricted to standard unperturbed 2D Keplerian motion (no $J_2$ Earth oblateness, no 3D inclination, no atmospheric drag, no lunar/solar gravity).
2. Unperturbed circular reference solution is valid only for ideal two-body motion.
3. Memory and disk storage considerations when exporting un-sampled full trajectories for fine timesteps over 100 orbits.

### Open Questions
1. How does the accumulated position error $e_r(100T)$ scale with timestep $\Delta t$?
2. Does relative energy error $\delta\varepsilon(t)$ grow linearly or non-linearly with time over $100 T$?

---

## 6. M4 — Phase Error and Radial Drift Decomposition

### Goal
Decompose the long-term Euclidean position error of the Explicit Euler integrator into orthogonal radial drift and along-track phase components over multi-orbit time scales ($t \in [0, 100 T]$). Quantitatively determine whether Euclidean error non-monotonicity is driven by radial drift, along-track phase error, or a scale-dependent regime transition.

### Mathematical Formulation
The unperturbed circular analytical state at time $t$ provides an orthonormal rotating frame basis $\{\hat{\mathbf{e}}_r(t), \hat{\mathbf{e}}_\theta(t)\}$:

$$\hat{\mathbf{e}}_r(t) = \begin{bmatrix} \cos(nt) \\ \sin(nt) \end{bmatrix}, \quad \hat{\mathbf{e}}_\theta(t) = \begin{bmatrix} -\sin(nt) \\ \cos(nt) \end{bmatrix}$$

Projecting the displacement vector $\Delta\mathbf{r}(t) = \mathbf{r}_{\text{num}}(t) - \mathbf{r}_{\text{exact}}(t)$ onto this basis yields:

* **Radial Error Projection**:
  $$e_R(t) = \Delta\mathbf{r}(t) \cdot \hat{\mathbf{e}}_r(t) = (x_{\text{num}} - x_{\text{exact}})\cos(nt) + (y_{\text{num}} - y_{\text{exact}})\sin(nt) \quad [\text{km}]$$
* **Along-Track Error Projection**:
  $$e_T(t) = \Delta\mathbf{r}(t) \cdot \hat{\mathbf{e}}_\theta(t) = -(x_{\text{num}} - x_{\text{exact}})\sin(nt) + (y_{\text{num}} - y_{\text{exact}})\cos(nt) \quad [\text{km}]$$
* **Exact Pythagorean Identity**:
  $$e_r^2(t) \equiv e_R^2(t) + e_T^2(t) \implies e_r(t) = \sqrt{e_R^2(t) + e_T^2(t)}$$
* **Variance Partition Fractions**:
  $$\rho_R(t) = \frac{e_R^2(t)}{e_r^2(t)}, \quad \rho_T(t) = \frac{e_T^2(t)}{e_r^2(t)}, \quad \rho_R(t) + \rho_T(t) \equiv 1.0 \quad (t > 0)$$
* **Dominance Ratio**:
  $$\chi(t) = \frac{|e_T(t)|}{|e_R(t)|}$$
* **Polar Coordinate Invariants**:
  $$\Delta r(t) = r_{\text{num}}(t) - r_0, \quad \Delta\theta(t) = \theta_{\text{unwrapped}}(t) - nt$$

### Implementation Details
* **Metrics Extension**: `src/simulation/metrics.py` extended with `radial_position_error`, `along_track_position_error`, `radial_drift`, `unwrap_phase`, `angular_phase_error`, `radial_error_fraction`, `along_track_error_fraction`, `dominance_ratio`, and `compute_decomposed_error_history`.
* **Experiment Runner**: `experiments/exp04_phase_radial_decomposition.py` reuses existing `Simulator` and Explicit Euler logic without modification.
* **Test Suite**: `tests/test_phase_radial_decomposition.py` with 13 comprehensive tests covering facts, identities, unwrap correctness, and M3 regression.

### Checkpoint Summary Results

| $\Delta t$ (s) | Ckpt | $e_r$ (km) | $e_R$ (km) | $e_T$ (km) | $\rho_R$ | $\rho_T$ | $\chi$ | Mode | $\Delta r$ (km) | $\Delta\theta$ (rad) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 60.000 | 1T | 16076.51 | -11814.33 | -10903.03 | 0.5401 | 0.4599 | 0.92 | radial | 5231.83 | -2.0035 |
| 60.000 | 5T | 8437.60 | 8416.87 | 591.09 | 0.9951 | 0.0049 | 0.07 | radial | 8428.36 | -18.8107 |
| 60.000 | 10T | 12177.48 | 7156.27 | 9852.86 | 0.3453 | 0.6547 | 1.38 | phase | 10287.81 | -43.3668 |
| 60.000 | 20T | 36268.21 | -19934.46 | -30298.52 | 0.3021 | 0.6979 | 1.52 | phase | 26253.50 | -96.2282 |
| 60.000 | 50T | 51019.15 | -30672.95 | -40769.16 | 0.3614 | 0.6386 | 1.33 | phase | 40477.41 | -265.9947 |
| 60.000 | 100T | 61782.18 | 4585.69 | -61611.77 | 0.0055 | 0.9945 | 13.44 | phase | 55872.85 | -560.5919 |
| 30.000 | 1T | 10232.37 | -4484.14 | -9197.49 | 0.1920 | 0.8080 | 2.05 | phase | 2701.12 | -1.3264 |
| 30.000 | 5T | 18366.42 | -13125.26 | -12847.30 | 0.5107 | 0.4893 | 0.98 | radial | 7551.52 | -14.5960 |
| 30.000 | 10T | 14552.91 | -2852.10 | 14270.70 | 0.0384 | 0.9616 | 5.00 | phase | 8022.76 | -36.3968 |
| 30.000 | 20T | 26796.61 | -26629.40 | 2988.84 | 0.9876 | 0.0124 | 0.11 | radial | 13296.87 | -84.9724 |
| 30.000 | 50T | 25822.48 | 13171.76 | -22210.47 | 0.2602 | 0.7398 | 1.69 | phase | 23076.57 | -245.8832 |
| 30.000 | 100T | 32930.88 | 21011.73 | -25356.46 | 0.4071 | 0.5929 | 1.21 | phase | 30841.37 | -528.5272 |
| 15.000 | 1T | 5932.93 | -1081.58 | -5833.51 | 0.0332 | 0.9668 | 5.39 | phase | 1375.43 | -0.7973 |
| 15.000 | 5T | 13541.95 | -8678.55 | 10395.53 | 0.4107 | 0.5893 | 1.20 | phase | 3789.68 | -10.8148 |
| 15.000 | 10T | 18362.81 | -16059.52 | 8904.20 | 0.7649 | 0.2351 | 0.55 | radial | 6083.77 | -29.0390 |
| 15.000 | 20T | 22978.33 | -22917.80 | 1666.79 | 0.9947 | 0.0053 | 0.07 | radial | 9447.36 | -72.3595 |
| 15.000 | 50T | 22491.36 | -1510.74 | -22440.57 | 0.0045 | 0.9955 | 14.85 | phase | 16272.34 | -221.2517 |
| 15.000 | 100T | 20598.18 | 20346.75 | -3208.50 | 0.9757 | 0.0243 | 0.16 | radial | 20535.85 | -490.2062 |
| 7.500 | 1T | 3210.43 | -24.97 | -3210.33 | 0.0001 | 0.9999 | 128.57 | phase | 699.26 | -0.4438 |
| 7.500 | 5T | 8765.13 | -2463.89 | -8411.70 | 0.0790 | 0.9210 | 3.41 | phase | 2675.40 | -7.3801 |
| 7.500 | 10T | 17478.70 | -16964.81 | -4207.15 | 0.9421 | 0.0579 | 0.25 | radial | 4243.14 | -21.5995 |
| 7.500 | 20T | 15171.83 | -6947.06 | -13487.88 | 0.2097 | 0.7903 | 1.94 | phase | 6710.80 | -58.1320 |
| 7.500 | 50T | 24234.59 | -24132.39 | -2223.27 | 0.9916 | 0.0084 | 0.09 | radial | 10717.95 | -191.5097 |
| 7.500 | 100T | 26967.08 | -20408.11 | -17627.62 | 0.5727 | 0.4273 | 0.86 | radial | 15504.34 | -442.0520 |
| 3.750 | 1T | 1672.96 | 157.22 | -1665.56 | 0.0088 | 0.9912 | 10.59 | phase | 354.41 | -0.2357 |
| 3.750 | 5T | 11314.48 | -7889.17 | 8110.38 | 0.4862 | 0.5138 | 1.03 | phase | 1407.99 | -4.5762 |
| 3.750 | 10T | 14282.78 | -11940.42 | -7837.35 | 0.6989 | 0.3011 | 0.66 | radial | 2606.60 | -14.7196 |
| 3.750 | 20T | 7338.75 | 1350.58 | 7213.41 | 0.0339 | 0.9661 | 5.34 | phase | 4089.67 | -43.2565 |
| 3.750 | 50T | 7321.84 | 7316.42 | -281.46 | 0.9985 | 0.0015 | 0.04 | radial | 7319.23 | -157.0996 |
| 3.750 | 100T | 10859.66 | 10174.96 | 3795.04 | 0.8779 | 0.1221 | 0.37 | radial | 10594.54 | -383.0541 |
| 1.875 | 1T | 853.37 | 127.41 | -843.80 | 0.0223 | 0.9777 | 6.62 | phase | 178.78 | -0.1216 |
| 1.875 | 5T | 13900.91 | -13387.09 | -3744.47 | 0.9274 | 0.0726 | 0.28 | radial | 817.87 | -2.6261 |
| 1.875 | 10T | 14896.84 | -14734.08 | -2196.09 | 0.9783 | 0.0217 | 0.15 | radial | 1475.33 | -9.1555 |
| 1.875 | 20T | 13434.78 | -10343.35 | 8573.71 | 0.5927 | 0.4073 | 0.83 | radial | 2507.29 | -29.4511 |
| 1.875 | 50T | 7209.34 | 2738.15 | -6669.12 | 0.1443 | 0.8557 | 2.44 | phase | 4842.40 | -119.9918 |
| 1.875 | 100T | 7319.22 | 7315.09 | 245.94 | 0.9989 | 0.0011 | 0.03 | radial | 7317.23 | -314.1418 |

### Scientific Classification of Results

#### 1. Tested Facts
* **Exact Orthogonality**: The rotating frame projections satisfy $e_r^2(t) = e_R^2(t) + e_T^2(t)$ across every discrete step of all simulations within floating-point tolerance ($< 10^{-10}\text{ km}^2$).
* **Variance Partition**: The variance fractions satisfy $\rho_R(t) + \rho_T(t) \equiv 1.0$ for all $t > 0$ within $10^{-14}$.
* **Monotonic Polar Invariants**: Radial drift $\Delta r(t)$ is strictly monotonically increasing ($\Delta r > 0$), and unwrapped angular phase error $\Delta\theta(t)$ is strictly monotonically decreasing ($\Delta\theta < 0$, persistent phase lag) across the entire 100-orbit integration for all tested timesteps.
* **Non-Monotonic Cartesian Projections**: Unlike the polar invariants, $e_R(t)$ and $e_T(t)$ oscillate and repeatedly change signs between positive and negative values.

#### 2. Experimental Observations
* **Short-Term Phase Dominance ($t \le 1T$)**: For fine timesteps ($\Delta t \le 15\text{ s}$), along-track phase error accounts for $> 96\%$ of total position variance ($\rho_T > 0.96$, $\chi \gg 1$) at $t = 1T$.
* **Long-Term Alternation**: Over multi-orbit horizons ($t \ge 5T$), neither component dominates monotonically. Both $\rho_R$ and $\rho_T$ periodically alternate between values close to $1.0$ and values close to $0.0$.
* **Phase Slip Frequency**: Finer timesteps experience slower secular phase slip. For example, over $100 T$:
  * $\Delta t = 60.0\text{ s} \implies \Delta\theta(100T) = -560.59\text{ rad} \approx -89.2\text{ revolutions}$.
  * $\Delta t = 1.875\text{ s} \implies \Delta\theta(100T) = -314.14\text{ rad} \approx -50.0\text{ revolutions}$.

#### 3. Scientific Interpretation and Resolution of Hypotheses
* **Supported Hypothesis — $H_C$ (Scale-Dependent Regime Transition and Phase-Wound Coupling)**:
  The experimental evidence refutes both pure phase dominance ($H_A$) and pure radial dominance ($H_B$). The error trajectory exhibits a scale-dependent transition driven by Keplerian orbital mechanics:
  1. Artificial energy injection increases orbital energy $\varepsilon(t)$, enlarging the numerical semi-major axis $a(t) > r_0$.
  2. By Kepler's third law, mean motion decreases: $n(a) = \sqrt{\mu / a^3} < n_0$. The satellite falls behind the reference body at rate $\Delta n(t) < 0$.
  3. Continuous integration produces a cumulative phase lag $\Delta\theta(t) < 0$.
  4. The Euclidean distance between the numerical and reference positions is governed by the law of cosines:
     $$e_r^2(t) = (\Delta r(t))^2 + 4 r_0 (r_0 + \Delta r(t))\sin^2\left(\frac{\Delta\theta(t)}{2}\right)$$
  5. As $\Delta\theta(t)$ winds through multiples of $\pi$, the position error vector rotates in the $\{\hat{\mathbf{e}}_r, \hat{\mathbf{e}}_\theta\}$ plane, causing $e_r(t)$ to oscillate between conjunction minima ($e_r \approx \Delta r$) and opposition maxima ($e_r \approx 2r_0 + \Delta r$).

* **Resolution of the M3 Convergence Paradox (Why smaller $\Delta t$ had larger $e_r$ at $t = 5T$)**:
  * At $t = 5T$, for $\Delta t = 60\text{ s}$, the phase error is $\Delta\theta = -18.8107\text{ rad} \approx -2.9938 \times (2\pi)$. The numerical satellite has completed almost exactly 3 extra phase laps, placing it in near-conjunction ($\Delta\theta \pmod{2\pi} \approx 2.2^\circ$). Along-track displacement vanishes ($e_T = 591.09\text{ km}$), leaving only the radial separation ($e_r \approx e_R \approx 8437.60\text{ km}$).
  * For $\Delta t = 1.875\text{ s}$, the phase error is $\Delta\theta = -2.6261\text{ rad} \approx -150.46^\circ$. The satellite is near orbital opposition (opposite side of Earth), where chord distance across the diameter maximizes the Euclidean error ($e_r \approx 13900.91\text{ km}$).
  * Thus, $\Delta t = 1.875\text{ s}$ has substantially *smaller* true physical drift ($\Delta r = 817.87\text{ km}$ vs. $8428.36\text{ km}$) and *smaller* phase slip ($|\Delta\theta| = 2.63\text{ rad}$ vs. $18.81\text{ rad}$), but orbital geometry across the diameter creates an apparent inversion when evaluating scalar Euclidean distance.

### Limitations
1. **Unperturbed 2D Keplerian Approximation**: Excludes $J_2$, out-of-plane cross-track error, and atmospheric drag.
2. **Integrator Specificity**: Results reflect the non-symplectic, non-conservative properties of the Explicit Euler integrator. Symplectic and higher-order Runge-Kutta integrators require separate investigation.

### Open Questions
1. How do higher-order non-symplectic integrators (e.g., Runge-Kutta 4th order) scale in radial drift versus phase error over $100 T$?
2. Does a symplectic integrator (e.g., Verlet / Leapfrog) eliminate the secular quadratic phase divergence by preserving an exact shadow Hamiltonian?
