# Milestone M8 — Systematic Integrator Comparison Methodology

## 1. Overview and Research Question

Milestone M8 establishes a rigorous, unified experimental methodology for systematically comparing the four numerical integrators developed across Milestones M1 through M7 of OrbitSim:
1. **Explicit Euler** (1st-order, non-symplectic, single-step)
2. **Velocity Verlet** (2nd-order, symplectic, time-reversible)
3. **Classical Runge-Kutta 4 (RK4)** (4th-order, non-symplectic, multi-stage)
4. **Yoshida 4** (4th-order, symplectic, symmetric composition)

### Core Research Question

> **"How do Explicit Euler, Velocity Verlet, Classical RK4, and Yoshida 4 compare under both matched-timestep ($\Delta t$) and matched-computational-budget (force evaluations) paradigms, and how does the choice of benchmark metric (instantaneous Cartesian truncation error vs. geometric invariant preservation vs. along-track phase drift) alter the scientific assessment of solver suitability for orbital mechanics?"**

### Guiding Principles

* **No Pre-ordained Winner**: The benchmark does not aggregate metrics into a single scalar score, nor does it declare a universally "best" method. Different orbital applications require different trade-offs (e.g. short-term landing precision vs. multi-decade ephemeris stability).
* **Exact Analytical Ground Truth**: The benchmark evaluates numerical errors against the exact closed-form Keplerian circular solution, not against another numerical integrator.
* **Separation of Discretization Error, Phase Drift, and Geometric Invariants**: Instantaneous Euclidean position error combines both shape distortion and along-track phase lag. Invariants reflect preservation of underlying phase-space geometry.

---

## 2. Hypotheses and Research Questions

The systematic comparison investigates five specific scientific questions:

1. **Short-Term Convergence Orders**:
   * *Hypothesis*: Over a single orbital period ($1T$), empirical convergence orders $p$ computed via dyadic timestep halving will match theoretical expectations: $p \approx 1$ (Euler), $p \approx 2$ (Velocity Verlet), $p \approx 4$ (RK4), and $p \approx 4$ (Yoshida 4).

2. **Cost-Normalized Efficiency (Pareto Trade-offs)**:
   * *Hypothesis*: When normalized to equal numbers of gravitational force evaluations, high-order methods (RK4 and Yoshida 4) will achieve lower absolute position errors than low-order methods (Euler and Velocity Verlet) at fine-to-moderate budgets, but the computational overhead per step will alter the cost-accuracy crossover points.

3. **Symplectic Invariant Preservation vs. High-Order Dissipation**:
   * *Hypothesis*: For the autonomous Hamiltonian system under consideration and fixed timestep, symplectic methods (Velocity Verlet and Yoshida 4) under appropriate backward error analysis conditions admit interpretation via a nearby modified Hamiltonian and typically exhibit bounded energy error over long horizons, in contrast to the typical secular energy drift of non-symplectic methods (Euler and RK4). This does not imply exact conservation of original energy and is not an unconditional guarantee for arbitrary problems, step sizes, or horizons.

4. **Phase Error vs. Radial Drift Decomposition**:
   * *Hypothesis*: Across stable integrators, total Euclidean position error over multi-orbit horizons ($t \gg T$) is dominated by along-track phase error ($e_T(t) \propto t \Delta t^p$), while radial drift ($e_R(t)$) remains bounded for symplectic methods but grows for non-symplectic methods.

5. **Angular Momentum Conservation**:
   * *Hypothesis*: Continuous physical conservation of angular momentum ($\frac{d\mathbf{h}}{dt} = \mathbf{r} \times \mathbf{F} = 0$) arises from central force symmetry. However, in discrete time, numerical integrators do not necessarily conserve the discrete cross product $\mathbf{r} \times \mathbf{v}$ (for instance, Explicit Euler generates an $\mathcal{O}(h^2 (\mathbf{v}_n \times \mathbf{a}_n))$ step error). M8 empirically measures angular momentum error $\delta h(t)$ across all methods rather than assuming its conservation.

---

## 3. Mathematical Reference Solution

### Orbital Parameters and Constants

To ensure consistency across M1–M8, all simulations use standard OrbitSim physical constants:
* Central body: Earth
* Gravitational parameter: $\mu = 398600.435\text{ km}^3/\text{s}^2$
* Mean equatorial radius: $R_E = 6378.137\text{ km}$
* Circular orbit altitude: $h_0 = 400.0\text{ km}$
* Orbital radius: $r_0 = R_E + h_0 = 6778.137\text{ km}$

### Exact Analytical State

For a circular orbit in the $xy$-plane with initial condition at the positive $x$-axis:
$$\mathbf{r}_0 = \begin{bmatrix} r_0 \\ 0 \end{bmatrix}, \quad \mathbf{v}_0 = \begin{bmatrix} 0 \\ v_0 \end{bmatrix}, \quad v_0 = \sqrt{\frac{\mu}{r_0}} \approx 7.672624\text{ km/s}$$

The orbital kinematics are known in closed analytical form for all $t \ge 0$:
* **Mean motion**:
  $$n = \sqrt{\frac{\mu}{r_0^3}} \approx 1.131370 \times 10^{-3}\text{ rad/s}$$
* **Orbital period**:
  $$T = \frac{2\pi}{n} \approx 5553.624319\text{ s} \quad (\approx 92.5604\text{ min})$$
* **Exact position**:
  $$x_{\text{exact}}(t) = r_0 \cos(nt), \quad y_{\text{exact}}(t) = r_0 \sin(nt)$$
* **Exact velocity**:
  $$v_{x,\text{exact}}(t) = -r_0 n \sin(nt), \quad v_{y,\text{exact}}(t) = r_0 n \cos(nt)$$
* **Exact specific orbital energy**:
  $$\varepsilon_0 = -\frac{\mu}{2 r_0} \approx -29.4035\text{ km}^2/\text{s}^2$$
* **Exact specific angular momentum**:
  $$h_0 = r_0 v_0 = \sqrt{\mu r_0} \approx 52006.126\text{ km}^2/\text{s}$$

### Rationale for Analytical Ground Truth

Using an analytical reference solution eliminates the need for a numerical "surrogate truth" (e.g. an 8th-order integrator with tiny $\Delta t$). The analytical solution:
1. Has **zero discretization truncation error** for all $t$.
2. Incurs only IEEE 754 double-precision roundoff in trigonometric evaluations ($\sim 10^{-16}$).
3. Introduces no numerical phase lag or numerical energy damping into the error baseline.

---

## 4. Simulation Horizons

The benchmark evaluates performance across two distinct physical timescales:

| Horizon | Duration ($t_f$) | Orbits | Primary Scientific Purpose |
| :--- | :--- | :--- | :--- |
| **Short-Term** | $1T \approx 5553.62\text{ s}$ | $1.0$ | One orbital period, used primarily to characterize discretization error and empirical convergence order. |
| **Long-Term** | $100T \approx 555255\text{ s}$ | $100.0$ | Invariant drift vs. boundedness, symplectic stability, secular phase accumulation, structural failure (divergence). |

### Checkpoint Schedule

For long-term runs, metrics are logged at discrete checkpoints:
$$\mathcal{T}_{\text{ckpt}} = \{1.0T, 5.0T, 10.0T, 20.0T, 50.0T, 100.0T\}$$
This schedule captures the transition from initial phase error accumulation to multi-orbit secular drift.

---

## 5. Timestep Grid ($\Delta t$)

To allow uniform comparison across all four methods, a common dyadic halving grid is defined:

$$\Delta t \in \{60.0\text{ s}, 30.0\text{ s}, 15.0\text{ s}, 7.5\text{ s}, 3.75\text{ s}, 1.875\text{ s}\}$$

### Rationale

* **Dyadic Ratios**: Successive step sizes satisfy $\Delta t_k / \Delta t_{k+1} = 2.0$, allowing direct convergence order computation:
  $$p = \log_2\left(\frac{e(\Delta t_k)}{e(\Delta t_{k+1})}\right)$$
* **Numerical Regime**:
  * For circular LEO ($T \approx 5554\text{ s}$), $\Delta t = 60\text{ s}$ corresponds to $\approx 92.5$ steps per orbit.
  * $\Delta t = 60\text{ s}$ is included because it lies within the previously tested regime for the current two-body model.
* **Dynamic Range**: Spans a factor of 32 in timestep size, providing 5 distinct convergence intervals for regression analysis.

---

## 6. Error and Fidelity Metrics

Let $\mathbf{s}_k = [x_k, y_k, v_{x,k}, v_{y,k}]^T$ denote the numerical state at time $t_k$, and $\mathbf{s}_{\text{exact}}(t_k) = [x_{\text{exact}}(t_k), y_{\text{exact}}(t_k), v_{x,\text{exact}}(t_k), v_{y,\text{exact}}(t_k)]^T$ the exact analytical state.

### 6.1 Instantaneous Kinematic Errors

* **Euclidean Position Error**:
  $$e_r(t_k) = \|\mathbf{r}_k - \mathbf{r}_{\text{exact}}(t_k)\|_2 = \sqrt{(x_k - x_{\text{exact}}(t_k))^2 + (y_k - y_{\text{exact}}(t_k))^2} \quad [\text{km}]$$
* **Euclidean Velocity Error**:
  $$e_v(t_k) = \|\mathbf{v}_k - \mathbf{v}_{\text{exact}}(t_k)\|_2 = \sqrt{(v_{x,k} - v_{x,\text{exact}}(t_k))^2 + (v_{y,k} - v_{y,\text{exact}}(t_k))^2} \quad [\text{km/s}]$$

### 6.2 Trajectory Aggregates vs. Checkpoints

To avoid conflating worst-case errors with final-state values:
* **Trajectory Maximum Error**:
  $$e_{r,\max} = \max_{0 \le t_k \le t_f} e_r(t_k), \quad e_{v,\max} = \max_{0 \le t_k \le t_f} e_v(t_k)$$
* **Final Checkpoint Error**:
  $$e_{r,\text{final}} = e_r(t_f), \quad e_{v,\text{final}} = e_v(t_f)$$

### 6.3 Rotating Frame Decomposition (Radial vs. Along-Track)

In orbital mechanics, errors naturally split into orbital shape changes (radial) and timing errors (along-track):
* **Radial unit vector**: $\mathbf{e}_r(t) = [\cos(nt), \sin(nt)]^T$
* **Along-track unit vector**: $\mathbf{e}_\theta(t) = [-\sin(nt), \cos(nt)]^T$
* **Radial Position Error Projection**:
  $$e_R(t_k) = (\mathbf{r}_k - \mathbf{r}_{\text{exact}}(t_k)) \cdot \mathbf{e}_r(t_k) \quad [\text{km}]$$
* **Along-Track Position Error Projection**:
  $$e_T(t_k) = (\mathbf{r}_k - \mathbf{r}_{\text{exact}}(t_k)) \cdot \mathbf{e}_\theta(t_k) \quad [\text{km}]$$
* **Radial Drift**:
  $$\Delta r(t_k) = r_k - r_0 = \sqrt{x_k^2 + y_k^2} - r_0 \quad [\text{km}]$$
* **Unwrapped Phase Error**:
  $$\Delta \theta(t_k) = \theta_{\text{unwrapped}}(t_k) - n t_k \quad [\text{rad}]$$
* **Dominance Ratio**:
  $$\chi(t_k) = \frac{|e_T(t_k)|}{|e_R(t_k)|}$$

---

## 7. Physical Invariants

### 7.1 Mathematical Definitions

* **Specific Orbital Energy**:
  $$\varepsilon_k = \frac{1}{2} v_k^2 - \frac{\mu}{r_k}, \quad \delta \varepsilon(t_k) = \frac{\varepsilon_k - \varepsilon_0}{|\varepsilon_0|}$$
* **Specific Angular Momentum**:
  $$h_{z,k} = x_k v_{y,k} - y_k v_{x,k}, \quad \delta h(t_k) = \frac{h_{z,k} - h_0}{|h_0|}$$

### 7.2 Strict Classification of Invariant Behavior

Scientific analysis must strictly distinguish between four distinct physical and numerical aspects:

1. **Continuous Physical Conservation**:
   * For the continuous two-body Keplerian problem, angular momentum is strictly conserved ($\frac{d\mathbf{h}}{dt} = \mathbf{r} \times \mathbf{F}(\mathbf{r}) = \mathbf{0}$) because gravitational acceleration is purely central. Energy is strictly conserved ($\frac{d\varepsilon}{dt} = 0$) because the Hamiltonian is autonomous and time-independent.
2. **Discrete Numerical Integrator Behavior**:
   * *Energy Invariance*: For the autonomous Hamiltonian system under consideration and fixed timestep, symplectic methods (Velocity Verlet, Yoshida 4) under appropriate backward error analysis conditions admit interpretation via a nearby modified Hamiltonian $\tilde{H} = H + \Delta t^p H_p + \dots$ and typically exhibit bounded energy error over long horizons, in contrast to the typical secular energy drift of non-symplectic methods (Euler, RK4). This does not imply exact conservation of original continuous energy and is not an unconditional guarantee for arbitrary problems, step sizes, or simulation horizons.
   * *Angular Momentum Invariance*: Continuous central force symmetry does not automatically imply exact discrete angular momentum preservation for an arbitrary discrete algorithm. For example, a single step of Explicit Euler updates state as $\mathbf{r}_{n+1} = \mathbf{r}_n + h \mathbf{v}_n$ and $\mathbf{v}_{n+1} = \mathbf{v}_n + h \mathbf{a}_n$. The resulting discrete cross product satisfies $\mathbf{r}_{n+1} \times \mathbf{v}_{n+1} = \mathbf{r}_n \times \mathbf{v}_n + h^2 (\mathbf{v}_n \times \mathbf{a}_n)$, yielding a systematic $\mathcal{O}(h^2)$ angular momentum error per step. M8 must empirically measure and track $\delta h(t)$ across each method without presupposing exact conservation.
3. **Floating-Point Roundoff Limits**:
   * Double-precision IEEE 754 arithmetic introduces a finite precision floor ($\epsilon_{\text{mach}} \approx 2.22 \times 10^{-16}$). Even when a numerical map algebraically preserves an invariant (e.g. Velocity Verlet on central forces), the accumulated deviation will be bounded by roundoff growth rather than truncation. In non-conservative solvers, truncation error completely dwarfs roundoff.
4. **Observed Experimental Characterization**:
   * The benchmark explicitly characterizes the empirical time evolution of invariants as:
     * **Secular Drift**: $|\delta \varepsilon(t)| \propto t^\alpha$ or $|\delta h(t)| \propto t^\beta$.
     * **Bounded Oscillation**: Invariant error remains strictly bounded within an amplitude envelope with zero secular trend.
     * **Roundoff-Limited**: Invariant error stays near the float64 arithmetic precision boundary.

---

## 8. Cost Normalization Model

### 8.1 Force Evaluation Proxy

Wall-clock execution time varies based on CPU architecture, background operating system tasks, and compiler optimizations. Therefore, **the number of gravitational acceleration (force) evaluations** is adopted as the primary, hardware-independent metric of computational cost.

This is an implementation-level force-evaluation cost proxy for the current codebase, not a universal model of wall-clock computational cost.

In the current modular implementation of OrbitSim:
* **Explicit Euler**: 1 evaluation per external step ($\mathbf{a}_n = \mathbf{a}(\mathbf{r}_n)$).
* **Velocity Verlet**: 2 evaluations per external step (un-cached evaluation at $\mathbf{r}_n$ and $\mathbf{r}_{n+1}$).
* **Classical RK4**: 4 evaluations per external step (stages $k_1, k_2, k_3, k_4$).
* **Yoshida 4**: 6 evaluations per external step (3 Velocity Verlet substeps $\times$ 2 evaluations each).

$$\mathbf{C}_{\text{step}} = [1, 2, 4, 6] \quad \text{for } [\text{Euler, Verlet, RK4, Yoshida 4}]$$

*Note: In specialized production codes, Velocity Verlet can be cached to 1 eval/step and Yoshida 4 to 3 evals/step across steps. However, in accordance with the project guidelines, cost accounting reflects the actual, current un-cached codebase.*

### 8.2 Matched Computational Budgets

To compare methods fairly at equal computational cost over $100T$ ($t_f \approx 555255\text{ s}$), the total number of force evaluations is held constant across solvers:
$$N_{\text{evals}} = N_{\text{steps}} \times \mathbf{C}_{\text{step}} \approx \frac{t_f}{\Delta t} \times \mathbf{C}_{\text{step}} \implies \Delta t_{\text{method}} = \Delta t_{\text{base}} \times \mathbf{C}_{\text{step}}$$

Four matched evaluation budgets are defined:

| Budget Name | Target Force Evals | Euler $\Delta t$ (s) | Verlet $\Delta t$ (s) | RK4 $\Delta t$ (s) | Yoshida 4 $\Delta t$ (s) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Budget 1 (Low Cost)** | $\approx 55,525$ | $10.0$ | $20.0$ | $40.0$ | $60.0$ |
| **Budget 2 (Medium Cost)** | $\approx 111,050$ | $5.0$ | $10.0$ | $20.0$ | $30.0$ |
| **Budget 3 (High Cost)** | $\approx 222,100$ | $2.5$ | $5.0$ | $10.0$ | $15.0$ |
| **Budget 4 (Ultra-High Cost)** | $\approx 444,200$ | $1.25$ | $2.5$ | $5.0$ | $7.5$ |

---

## 9. Comparison Rules and Trade-Off Dimensions

To ensure scientific integrity:

1. **Dual-Axis Comparison**:
   * Every metric must be presented under **Fixed Timestep** (mathematical behavior) and **Fixed Evaluation Budget** (computational efficiency).
2. **Multi-Faceted Assessment**:
   * Evaluation must consider all four primary dimensions:
     1. Discretization accuracy ($e_r, e_v$).
     2. Invariant preservation ($\delta \varepsilon, \delta h$).
     3. Geometric stability ($e_R$ vs. $e_T$).
     4. Computational complexity (force evaluations per unit simulated time).
3. **No Aggregate Ranking**:
   * No solver will be labeled "best". The benchmark must characterize the Pareto frontier:
     * *Euler*: Minimal implementation complexity, lowest cost per step, but mathematically unstable for orbital propagation.
     * *RK4*: Highest short-term point accuracy at large timesteps, but non-conservative and dissipative over long horizons.
     * *Velocity Verlet*: Low cost per step (2 evals), strictly symplectic with bounded energy, but moderate 2nd-order truncation error.
     * *Yoshida 4*: 4th-order discretization accuracy combined with strict symplecticity, but higher per-step cost (6 evals).

---

## 10. Planned CSV Output Schemas

Four standardized CSV schemas are defined for Milestone M8:

### 10.1 Convergence Study (`exp08_convergence.csv`)
Logs 1-orbit short-term convergence data:
`method,dt,num_steps,force_evals,e_r,e_v,order_p_r,ratio_r,max_de,max_dh,final_dr,final_dtheta,max_e_r,max_e_v,order_p_v,final_de,final_dh,max_dr`

### 10.2 Long-Term Fixed Timestep (`exp08_long_term_fixed_dt.csv`)
Logs 100-orbit performance at fixed $\Delta t \in \{60.0, 30.0, 15.0\}\text{ s}$:
`method,dt,evals_per_step,total_steps,actual_force_evals,max_position_error,final_position_error,max_velocity_error,final_velocity_error,max_radial_deviation,final_radial_deviation,final_phase_error,max_relative_energy_error,final_relative_energy_error,max_relative_angular_momentum_error`

### 10.3 Cost-Normalized Comparison (`exp08_cost_normalized.csv`)
Logs 100-orbit performance across matched evaluation budgets:
`budget_name,target_evals,method,dt,evals_per_step,total_steps,actual_force_evals,max_position_error,final_position_error,max_radial_deviation,final_radial_deviation,final_phase_error,max_relative_energy_error,final_relative_energy_error,max_relative_angular_momentum_error`

### 10.4 Invariant Quality Classification (`exp08_invariant_quality.csv`)
Logs qualitative and quantitative invariant metrics:
`method,dt,max_de_trajectory,final_de_checkpoint,de_character,max_dh_trajectory,final_dh_checkpoint,dh_character`

---

## 11. Limitations and Reproducibility

### Current Model Limitations

1. **Unperturbed Two-Body Dynamics**: Benchmarked exclusively on planar, circular Keplerian orbits ($\ddot{\mathbf{r}} = -\mu \mathbf{r} / r^3$). Perturbations ($J_2$, atmospheric drag, third bodies) and eccentric orbits are reserved for future milestones.
2. **Fixed Timesteps**: All solvers operate with constant $\Delta t$; adaptive step-size control is not benchmarked.
3. **Hardware Independence**: Timing measurements are excluded from quantitative claims in favor of deterministic force evaluation counts.

### Reproducibility Protocol

All results are reproducible by executing the standardized suite with fixed initial conditions:
* Initial state: $\mathbf{r}_0 = [6778.137, 0.0]^T\text{ km}$, $\mathbf{v}_0 = [0.0, \sqrt{\mu/r_0}]^T\text{ km/s}$
* Floating-point environment: Standard IEEE 754 double precision (64-bit float).

---

## 12. Executed Milestone M8.1: Analytical Convergence Study Results

### 12.1 Execution Parameters
* **Execution Timestamp**: 2026-09-26 (calibrated on baseline commit `6c2b577`)
* **Reference Solution**: Exact closed-form Keplerian circular solution ($r_0 = 6778.137\text{ km}$, $n \approx 1.131370 \times 10^{-3}\text{ rad/s}$, $T \approx 5553.624319\text{ s}$)
* **Number of Solvers**: Exactly 4 (`Explicit Euler`, `Velocity Verlet`, `Classical RK4`, `Yoshida 4`)
* **Timestep Grid**: Exactly 6 timesteps per solver: $\Delta t \in \{60.0, 30.0, 15.0, 7.5, 3.75, 1.875\}\text{ s}$
* **Total Executed Simulations**: 24 runs (4 solvers $\times$ 6 timesteps)
* **Integration Horizon**: Exactly 1 orbital period ($t_f = T \approx 5553.624319\text{ s}$), identical across all runs
* **Generated Output Artifacts**:
  * `results/exp08_convergence.csv` (24 data rows + header)
  * `results/exp08_convergence_position.png` (log-log position error vs $\Delta t$)
  * `results/exp08_convergence_velocity.png` (log-log velocity error vs $\Delta t$)
  * `results/exp08_convergence_invariants.png` (energy and angular momentum errors vs $\Delta t$)

### 12.2 Measured Convergence Data Table

| Method | $\Delta t$ [s] | Steps | Evals | $e_r(T)$ [km] | $e_v(T)$ [km/s] | Order $p_r$ | Order $p_v$ | $\max |\delta\varepsilon|$ | $\max |\delta h|$ | Final $\Delta\theta$ [rad] |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Euler** | 60.0 | 93 | 93 | $1.5978 \times 10^{4}$ | $1.1313 \times 10^{1}$ | — | — | $3.6382 \times 10^{-1}$ | $2.3854 \times 10^{-1}$ | $1.9855 \times 10^{0}$ |
| **Euler** | 30.0 | 186 | 186 | $1.0243 \times 10^{4}$ | $8.7365 \times 10^{0}$ | 0.64 | 0.37 | $2.4212 \times 10^{-1}$ | $1.4653 \times 10^{-1}$ | $1.3281 \times 10^{0}$ |
| **Euler** | 15.0 | 371 | 371 | $5.9398 \times 10^{3}$ | $5.7017 \times 10^{0}$ | 0.79 | 0.62 | $1.5179 \times 10^{-1}$ | $8.5569 \times 10^{-2}$ | $7.9830 \times 10^{-1}$ |
| **Euler** | 7.5 | 741 | 741 | $3.2143 \times 10^{3}$ | $3.3052 \times 10^{0}$ | 0.89 | 0.79 | $8.8250 \times 10^{-2}$ | $4.7260 \times 10^{-2}$ | $4.4433 \times 10^{-1}$ |
| **Euler** | 3.75 | 1481 | 1481 | $1.6729 \times 10^{3}$ | $1.7867 \times 10^{0}$ | 0.94 | 0.89 | $4.8234 \times 10^{-2}$ | $2.5025 \times 10^{-2}$ | $2.3568 \times 10^{-1}$ |
| **Euler** | 1.875 | 2962 | 2962 | $8.5333 \times 10^{2}$ | $9.2990 \times 10^{-1}$ | 0.97 | 0.94 | $2.5315 \times 10^{-2}$ | $1.2903 \times 10^{-2}$ | $1.2158 \times 10^{-1}$ |
| **Velocity Verlet** | 60.0 | 93 | 186 | $6.5437 \times 10^{1}$ | $7.3781 \times 10^{-2}$ | — | — | $5.2714 \times 10^{-6}$ | $4.1994 \times 10^{-16}$ | $9.6542 \times 10^{-3}$ |
| **Velocity Verlet** | 30.0 | 186 | 372 | $1.6352 \times 10^{1}$ | $1.8490 \times 10^{-2}$ | 2.00 | 2.00 | $3.3118 \times 10^{-7}$ | $1.1198 \times 10^{-15}$ | $2.4124 \times 10^{-3}$ |
| **Velocity Verlet** | 15.0 | 371 | 742 | $4.0893 \times 10^{0}$ | $4.6242 \times 10^{-3}$ | 2.00 | 2.00 | $2.0727 \times 10^{-8}$ | $1.2598 \times 10^{-15}$ | $6.0330 \times 10^{-4}$ |
| **Velocity Verlet** | 7.5 | 741 | 1482 | $1.0223 \times 10^{0}$ | $1.1562 \times 10^{-3}$ | 2.00 | 2.00 | $1.2958 \times 10^{-9}$ | $1.3998 \times 10^{-15}$ | $1.5083 \times 10^{-4}$ |
| **Velocity Verlet** | 3.75 | 1481 | 2962 | $2.5553 \times 10^{-1}$ | $2.8909 \times 10^{-4}$ | 2.00 | 2.00 | $8.0990 \times 10^{-11}$ | $4.7593 \times 10^{-15}$ | $3.7699 \times 10^{-5}$ |
| **Velocity Verlet** | 1.875 | 2962 | 5924 | $6.3883 \times 10^{-2}$ | $7.2273 \times 10^{-5}$ | 2.00 | 2.00 | $5.0647 \times 10^{-12}$ | $2.0997 \times 10^{-15}$ | $9.4249 \times 10^{-6}$ |
| **RK4** | 60.0 | 93 | 372 | $2.8783 \times 10^{-2}$ | $3.2435 \times 10^{-5}$ | — | — | $2.5036 \times 10^{-7}$ | $1.2518 \times 10^{-7}$ | $4.2391 \times 10^{-6}$ |
| **RK4** | 30.0 | 186 | 744 | $1.5464 \times 10^{-3}$ | $1.7482 \times 10^{-6}$ | 4.22 | 4.21 | $7.8581 \times 10^{-9}$ | $3.9291 \times 10^{-9}$ | $2.2801 \times 10^{-7}$ |
| **RK4** | 15.0 | 371 | 1484 | $8.8815 \times 10^{-5}$ | $1.0043 \times 10^{-7}$ | 4.12 | 4.12 | $2.4552 \times 10^{-10}$ | $1.2276 \times 10^{-10}$ | $1.3101 \times 10^{-8}$ |
| **RK4** | 7.5 | 741 | 2964 | $5.3059 \times 10^{-6}$ | $6.0003 \times 10^{-9}$ | 4.07 | 4.06 | $7.6783 \times 10^{-12}$ | $3.8391 \times 10^{-12}$ | $7.8276 \times 10^{-10}$ |
| **RK4** | 3.75 | 1481 | 5924 | $3.2393 \times 10^{-7}$ | $3.6649 \times 10^{-10}$ | 4.03 | 4.03 | $2.4455 \times 10^{-13}$ | $1.2234 \times 10^{-13}$ | $4.7809 \times 10^{-11}$ |
| **RK4** | 1.875 | 2962 | 11848 | $1.9976 \times 10^{-8}$ | $2.2575 \times 10^{-11}$ | 4.02 | 4.02 | $7.6121 \times 10^{-15}$ | $3.7795 \times 10^{-15}$ | $2.9372 \times 10^{-12}$ |
| **Yoshida 4** | 60.0 | 93 | 558 | $8.3884 \times 10^{-1}$ | $9.4614 \times 10^{-4}$ | — | — | $4.7836 \times 10^{-10}$ | $4.1994 \times 10^{-16}$ | $1.2376 \times 10^{-4}$ |
| **Yoshida 4** | 30.0 | 186 | 1116 | $5.2354 \times 10^{-2}$ | $5.9210 \times 10^{-5}$ | 4.00 | 4.00 | $1.8699 \times 10^{-12}$ | $1.3998 \times 10^{-15}$ | $7.7240 \times 10^{-6}$ |
| **Yoshida 4** | 15.0 | 371 | 2226 | $3.2715 \times 10^{-3}$ | $3.6999 \times 10^{-6}$ | 4.00 | 4.00 | $1.2324 \times 10^{-14}$ | $3.6395 \times 10^{-15}$ | $4.8265 \times 10^{-7}$ |
| **Yoshida 4** | 7.5 | 741 | 4446 | $2.0446 \times 10^{-4}$ | $2.3123 \times 10^{-7}$ | 4.00 | 4.00 | $5.5580 \times 10^{-15}$ | $2.7996 \times 10^{-15}$ | $3.0164 \times 10^{-8}$ |
| **Yoshida 4** | 3.75 | 1481 | 8886 | $1.2778 \times 10^{-5}$ | $1.4456 \times 10^{-8}$ | 4.00 | 4.00 | $2.3803 \times 10^{-14}$ | $1.1898 \times 10^{-14}$ | $1.8852 \times 10^{-9}$ |
| **Yoshida 4** | 1.875 | 2962 | 17772 | $7.9828 \times 10^{-7}$ | $9.0321 \times 10^{-10}$ | 4.00 | 4.00 | $1.6070 \times 10^{-14}$ | $8.1188 \times 10^{-15}$ | $1.1783 \times 10^{-10}$ |

### 12.3 Comparison: Theoretical Expectations vs. Observed Results

1. **Explicit Euler**:
   * *Theoretical Expectation*: $\mathcal{O}(\Delta t)$ asymptotic convergence ($p = 1$).
   * *Observed Result*: Empirical order $p_r$ transitions from $0.64$ at coarse steps ($\Delta t = 60 \to 30\text{ s}$) to $0.97$ at fine steps ($\Delta t = 3.75 \to 1.875\text{ s}$). Large secular drifts observed in energy ($\max |\delta\varepsilon| \approx 0.36$ at $60\text{ s}$) and angular momentum ($\max |\delta h| \approx 0.24$ at $60\text{ s}$).
2. **Velocity Verlet**:
   * *Theoretical Expectation*: $\mathcal{O}(\Delta t^2)$ convergence ($p = 2$). Symplectic preservation of modified Hamiltonian. Continuous angular momentum symmetry preserved in discrete updates.
   * *Observed Result*: Empirical order is strictly $p_r = 2.00$ and $p_v = 2.00$ across all 5 halving intervals. Relative angular momentum deviation $\max |\delta h|$ remains bounded between $4.2 \times 10^{-16}$ and $4.8 \times 10^{-15}$ (at the IEEE 754 float64 machine epsilon limit). Energy error scales as $\Delta t^2$ ($5.27 \times 10^{-6}$ at $60\text{ s}$ down to $5.06 \times 10^{-12}$ at $1.875\text{ s}$).
3. **Classical RK4**:
   * *Theoretical Expectation*: $\mathcal{O}(\Delta t^4)$ convergence ($p = 4$). Non-symplectic multi-stage formulation.
   * *Observed Result*: Empirical order $p_r$ is $4.22$ at $60 \to 30\text{ s}$ and steadily approaches $4.02$ at $3.75 \to 1.875\text{ s}$. Achieves the smallest absolute Cartesian position error at equal $\Delta t$ ($1.99 \times 10^{-8}\text{ km}$ at $1.875\text{ s}$). Angular momentum error $\max |\delta h|$ is not machine-limited; it exhibits non-symplectic truncation error scaling as $\mathcal{O}(\Delta t^4)$ ($1.25 \times 10^{-7}$ at $60\text{ s}$ down to $3.78 \times 10^{-15}$ at $1.875\text{ s}$).
4. **Yoshida 4**:
   * *Theoretical Expectation*: $\mathcal{O}(\Delta t^4)$ convergence ($p = 4$). Symmetric symplectic composition.
   * *Observed Result*: Empirical order is strictly $p_r = 4.00$ and $p_v = 4.00$ across all halving intervals. Angular momentum error $\max |\delta h|$ is strictly machine-precision bounded ($\sim 10^{-16} - 10^{-14}$) across all timesteps. Relative energy error is bounded and reaches the float64 roundoff floor ($\sim 10^{-14}$) for $\Delta t \le 15\text{ s}$.

### 12.4 M8.1 Experiment Scope and Limitations
* **Single Orbital Horizon**: The M8.1 benchmark covers strictly one orbital period ($1T$). Long-term secular error growth, stability horizons, and structural divergence over $100T$ are part of Milestone M8.2 and were not evaluated in this experiment.
* **Matched-Timestep vs. Matched-Cost**: The results above reflect identical $\Delta t$ grids. Because the four solvers incur different numbers of force evaluations per step ($1, 2, 4, 6$), equal-$\Delta t$ comparisons do not reflect equal computational cost. Cost-normalized evaluations are deferred to M8.2.
