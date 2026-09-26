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

### Rationale and Endpoint Handling

* **Dyadic Ratios**: Successive step sizes satisfy $\Delta t_k / \Delta t_{k+1} = 2.0$, allowing direct convergence order computation:
  $$p = \log_2\left(\frac{e(\Delta t_k)}{e(\Delta t_{k+1})}\right)$$
* **Nominal Timestep Definition**: In convergence analysis, $\Delta t$ denotes the **nominal integration step**.
* **Common Physical Endpoint**: Each run uses the nominal timestep $\Delta t$ for all full steps and a shortened final step when necessary so that the trajectory terminates at the common physical time $T$. Because $T \approx 5553.624319\text{ s}$ is non-integer in seconds, the final step $\Delta t_{\text{last}} = T - N_{\text{full}} \Delta t < \Delta t$ hits $t_f = T$ exactly, preventing artificial time-offset phase errors.
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
Logs 1-orbit short-term convergence data across all 4 solvers and 6 nominal timesteps (24 rows total):
`method,dt,num_steps,force_evals,e_r,e_v,order_p_r,ratio_r,max_de,max_dh,final_dr,final_dtheta,max_e_r,max_e_v,order_p_v,final_de,final_dh,max_dr`

* **`method`**: Numerical integrator identifier (`Euler`, `Velocity Verlet`, `RK4`, `Yoshida 4`).
* **`dt`**: Nominal integration timestep $\Delta t$ in seconds. All full steps use $\Delta t$, with a shortened step at the end to terminate at $T$.
* **`num_steps`**: Total discrete integration steps taken to reach physical time $T = 5553.624319\text{ s}$.
* **`force_evals`**: Total gravitational force evaluations ($N_{\text{steps}} \times \mathbf{C}_{\text{step}}$).
* **`e_r`**: **Final Euclidean position error** at physical time $T$: $\|\mathbf{r}_{\text{num}}(T) - \mathbf{r}_{\text{exact}}(T)\|_2$ [km].
* **`e_v`**: **Final Euclidean velocity error** at physical time $T$: $\|\mathbf{v}_{\text{num}}(T) - \mathbf{v}_{\text{exact}}(T)\|_2$ [km/s].
* **`order_p_r`**: Empirical convergence order calculated from final position error $e_r$ between adjacent nominal timesteps: $\log_2(e_r(\Delta t) / e_r(\Delta t / 2))$.
* **`ratio_r`**: Error reduction ratio between adjacent nominal timesteps: $e_r(\Delta t) / e_r(\Delta t / 2)$.
* **`max_de`**: Maximum absolute relative specific orbital energy error over the solver's discrete trajectory: $\max_{t_k} |(\varepsilon(t_k) - \varepsilon_0)/\varepsilon_0|$.
* **`max_dh`**: Maximum absolute relative specific angular momentum error over the solver's discrete trajectory: $\max_{t_k} |(h_z(t_k) - h_0)/h_0|$.
* **`final_dr`**: Final absolute radial deviation at physical time $T$: $|r_{\text{num}}(T) - r_0|$ [km].
* **`final_dtheta`**: Final along-track unwrapped phase error at physical time $T$: $|\Delta\theta(T)| = |\theta_{\text{unwrapped}}(T) - n T|$ [rad].
* **`max_e_r`**: **Maximum Euclidean position error** evaluated over the solver's own discrete trajectory: $\max_{t_k} \|\mathbf{r}_{\text{num}}(t_k) - \mathbf{r}_{\text{exact}}(t_k)\|_2$ [km].
* **`max_e_v`**: **Maximum Euclidean velocity error** evaluated over the solver's own discrete trajectory: $\max_{t_k} \|\mathbf{v}_{\text{num}}(t_k) - \mathbf{v}_{\text{exact}}(t_k)\|_2$ [km/s].
* **`order_p_v`**: Empirical convergence order calculated from final velocity error $e_v$ between adjacent nominal timesteps: $\log_2(e_v(\Delta t) / e_v(\Delta t / 2))$.
* **`final_de`**: Final signed relative specific orbital energy error at physical time $T$: $(\varepsilon(T) - \varepsilon_0)/|\varepsilon_0|$.
* **`final_dh`**: Final signed relative specific angular momentum error at physical time $T$: $(h_z(T) - h_0)/|h_0|$.
* **`max_dr`**: Maximum absolute radial deviation evaluated over the solver's own discrete trajectory: $\max_{t_k} |r_{\text{num}}(t_k) - r_0|$ [km].


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
* **Step Composition**: Each run executes $N_{\text{full}} = \lfloor T / \Delta t \rfloor$ nominal steps of size $\Delta t$, followed by a single shortened step $\Delta t_{\text{last}} = T - N_{\text{full}} \Delta t$ terminating precisely at $t_f = T$ (e.g. 92 steps of $60.0\text{ s}$ + 1 step of $33.624\text{ s}$).
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

---

## 13. Executed Milestone M8.2a: Long-Term Fixed-dt Comparison Results (100 Orbits)

### 13.1 Execution Parameters
* **Execution Timestamp**: 2026-09-26
* **Simulation Horizon**: Exactly 100 orbital periods:
  $$t_f = 100 \times T = 555362.4318623783\text{ s} \quad (\approx 6.4278\text{ days})$$
* **Solvers Evaluated**: 4 methods (`Explicit Euler`, `Velocity Verlet`, `Classical RK4`, `Yoshida 4`).
* **Nominal Timestep Grid**: $\Delta t \in \{60.0, 30.0, 15.0\}\text{ s}$.
* **Total Executed Simulations**: 12 runs (4 solvers $\times$ 3 timesteps).
* **Step Composition**: Each run integrates $N_{\text{full}} = \lfloor 100 T / \Delta t \rfloor$ nominal full steps, followed by a single shortened final step:
  $$\Delta t_{\text{last}} = 100 T - 555360.0\text{ s} \approx 2.431862378\text{ s}$$
  All 12 runs terminate at $t_f = 100 T$ with $|t_{\text{final}} - 100 T| < 10^{-11}\text{ s}$.
  * $\Delta t = 60.0\text{ s}$: 9256 full steps + 1 step of $2.431862\text{ s}$ = 9257 steps.
  * $\Delta t = 30.0\text{ s}$: 18512 full steps + 1 step of $2.431862\text{ s}$ = 18513 steps.
  * $\Delta t = 15.0\text{ s}$: 37024 full steps + 1 step of $2.431862\text{ s}$ = 37025 steps.
* **Discrete Checkpoints**: Sampled at $\{1.0T, 5.0T, 10.0T, 20.0T, 50.0T, 100.0T\}$. At each checkpoint, analytical reference kinematics are evaluated at the actual sampled discrete timestamp $t_{\text{actual}}$.
* **Generated Output Artifacts**:
  * `results/exp08_long_term_fixed_dt.csv` (12 trajectory-wide summary rows + header)
  * `results/exp08_long_term_checkpoints.csv` (72 checkpoint rows + header)
  * `results/exp08_long_term_position_error.png`
  * `results/exp08_long_term_energy_deviation.png`
  * `results/exp08_long_term_angular_momentum.png`
  * `results/exp08_long_term_phase_error.png`

### 13.2 Summary Results Table (100 Orbits)

| Method | $\Delta t$ [s] | Steps | Evals | Max $e_r$ [km] | Final $e_r$ [km] | Final $\Delta\theta$ [rad] | $\max |\delta\varepsilon|$ | Final $\delta\varepsilon$ | $\max |\delta h|$ | Final $\delta h$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Euler** | 60.0 | 9257 | 9257 | $7.1477 \times 10^{4}$ | $6.1782 \times 10^{4}$ | $-560.59$ | $0.8566$ | $+0.8566$ | $1.4543$ | $+1.4543$ |
| **Euler** | 30.0 | 18513 | 18513 | $4.7068 \times 10^{4}$ | $3.2931 \times 10^{4}$ | $-528.53$ | $0.8050$ | $+0.8050$ | $1.2342$ | $+1.2342$ |
| **Euler** | 15.0 | 37025 | 37025 | $3.5824 \times 10^{4}$ | $2.0598 \times 10^{4}$ | $-490.21$ | $0.7518$ | $+0.7518$ | $1.0024$ | $+1.0024$ |
| **Velocity Verlet** | 60.0 | 9257 | 18514 | $6.2835 \times 10^{3}$ | $6.2806 \times 10^{3}$ | $-0.9634$ | $5.2719 \times 10^{-6}$ | $+7.99 \times 10^{-8}$ | $5.18 \times 10^{-15}$ | $+3.50 \times 10^{-15}$ |
| **Velocity Verlet** | 30.0 | 18513 | 37026 | $1.6314 \times 10^{3}$ | $1.6307 \times 10^{3}$ | $-0.2412$ | $3.3120 \times 10^{-7}$ | $+3.31 \times 10^{-10}$ | $8.26 \times 10^{-15}$ | $-2.52 \times 10^{-15}$ |
| **Velocity Verlet** | 15.0 | 37025 | 74050 | $4.0892 \times 10^{2}$ | $4.0874 \times 10^{2}$ | $-0.0603$ | $2.0727 \times 10^{-8}$ | $+1.64 \times 10^{-12}$ | $8.68 \times 10^{-15}$ | $+2.10 \times 10^{-15}$ |
| **RK4** | 60.0 | 9257 | 37028 | $8.2499 \times 10^{1}$ | $8.2499 \times 10^{1}$ | $+0.0122$ | $2.5183 \times 10^{-5}$ | $-2.5183 \times 10^{-5}$ | $1.2591 \times 10^{-5}$ | $-1.2591 \times 10^{-5}$ |
| **RK4** | 30.0 | 18513 | 74052 | $2.6410 \times 10^{0}$ | $2.6410 \times 10^{0}$ | $+3.90 \times 10^{-4}$ | $7.8632 \times 10^{-7}$ | $-7.8632 \times 10^{-7}$ | $3.9316 \times 10^{-7}$ | $-3.9316 \times 10^{-7}$ |
| **RK4** | 15.0 | 37025 | 148100 | $8.6568 \times 10^{-2}$ | $8.6568 \times 10^{-2}$ | $+1.28 \times 10^{-5}$ | $2.4568 \times 10^{-8}$ | $-2.4568 \times 10^{-8}$ | $1.2284 \times 10^{-8}$ | $-1.2284 \times 10^{-8}$ |
| **Yoshida 4** | 60.0 | 9257 | 55542 | $8.3837 \times 10^{1}$ | $8.3832 \times 10^{1}$ | $+0.0124$ | $4.7840 \times 10^{-10}$ | $+7.01 \times 10^{-15}$ | $1.526 \times 10^{-14}$ | $-1.34 \times 10^{-14}$ |
| **Yoshida 4** | 30.0 | 18513 | 111078 | $5.2353 \times 10^{0}$ | $5.2351 \times 10^{0}$ | $+7.72 \times 10^{-4}$ | $1.9291 \times 10^{-12}$ | $+6.02 \times 10^{-14}$ | $3.094 \times 10^{-14}$ | $+3.00 \times 10^{-14}$ |
| **Yoshida 4** | 15.0 | 37025 | 222150 | $3.2714 \times 10^{-1}$ | $3.2712 \times 10^{-1}$ | $+4.83 \times 10^{-5}$ | $6.7300 \times 10^{-14}$ | $+5.78 \times 10^{-14}$ | $3.024 \times 10^{-14}$ | $+2.90 \times 10^{-14}$ |

*Note: All trajectory-wide maxima (`max_*`) are evaluated strictly over the discrete sampled states along each solver's trajectory, not over continuous time.*

### 13.3 Scientific Observations and Analysis

1. **Explicit Euler — Strong Unbounded Growth**:
   * Euler exhibits strong unbounded growth of trajectory error, energy deviation, and radial error over the tested horizon ($\delta\varepsilon \approx +0.86$, $\delta h \approx +1.45$, and $\Delta r > 20000\text{ km}$ at $\Delta t = 60\text{ s}$).
   * Massive along-track phase errors ($|\Delta\theta| > 490\text{ rad}$) and trajectory breakdown occur across all tested timesteps.
   * This confirms the known numerical instability of Explicit Euler for orbital propagation without assuming an unverified specific growth law.

2. **Velocity Verlet — Symplectic Boundedness vs. Phase Drift**:
   * *Energy*: No secular trend was detected over the tested 100-orbit horizon; the relative energy deviation remained bounded and oscillatory within an $\mathcal{O}(\Delta t^2)$ envelope ($\max |\delta\varepsilon| \approx 5.27 \times 10^{-6}$ at $60\text{ s}$, $2.07 \times 10^{-8}$ at $15\text{ s}$).
   * *Angular Momentum*: In the conducted experiments, relative angular-momentum deviation for Velocity Verlet remained at the float64 numerical-precision level ($\max |\delta h| < 8.7 \times 10^{-15}$). This is consistent with central force symmetry and the structure of the discrete update.
   * *Radial Deviation*: Orbital radial deviation remains comparatively small over 100 orbits ($\Delta r(100T) \approx 0.23\text{ km}$ at $60\text{ s}$, $5.4 \times 10^{-5}\text{ km}$ at $15\text{ s}$).
   * *Phase Error vs. Radial Error*: The dominant contribution to Cartesian position error is along-track phase drift, while radial deviation remains comparatively small. Phase error accumulates over time, reaching $|\Delta\theta| \approx 0.963\text{ rad}$ at $60\text{ s}$, leading to a Cartesian position error of $e_r \approx 6280\text{ km}$ (where chord distance $2 r_0 \sin(|\Delta\theta| / 2) \approx 6280.5\text{ km}$).

3. **Classical RK4 — Dissipative Secular Drift and Radial Bias**:
   * *Energy and Invariants*: RK4 exhibits a measurable secular drift in energy and angular momentum over the tested 100-orbit horizon ($\delta\varepsilon(100T) \approx -2.52 \times 10^{-5}$, $\delta h(100T) \approx -1.26 \times 10^{-5}$ at $\Delta t = 60\text{ s}$), accompanied by a small monotonic inward radial bias ($\Delta r(100T) \approx -0.171\text{ km} = -171\text{ m}$ at $\Delta t = 60\text{ s}$).
   * *No Center Fall / Extrapolation*: The observed inward radial deviation is bounded to 171 m after 100 orbits at $\Delta t = 60\text{ s}$ and does not constitute a rapid collapse or spiral into the center on the tested horizon.
   * *Position Error*: Due to 4th-order local truncation error, Cartesian position error remains relatively small over 100 orbits ($e_r \approx 82.5\text{ km}$ at $60\text{ s}$, $0.0866\text{ km}$ at $15\text{ s}$).

4. **Yoshida 4 — Fourth-Order Symplectic Fidelity**:
   * *Energy*: Energy error remained bounded and oscillatory within the expected 4th-order modified Hamiltonian envelope ($\max |\delta\varepsilon| \approx 4.78 \times 10^{-10}$ at $60\text{ s}$); no secular trend was detected over the tested 100-orbit horizon. At $\Delta t = 15\text{ s}$, energy error reached the float64 arithmetic precision limit ($\sim 6.7 \times 10^{-14}$). (Note: $\mathcal{O}(\Delta t^4)$ describes the expected modified-Hamiltonian scaling structure; the special empirical $\mathcal{O}(\Delta t^8)$ scaling identified in M7 under specific conditions is not extrapolated to M8.2a).
   * *Angular Momentum*: In the conducted experiments, relative angular-momentum deviation for Yoshida 4 remained at the float64 numerical-precision level ($\max |\delta h| < 3.1 \times 10^{-14}$), consistent with central force symmetry.
   * *Radial Deviation and Units*: The maximum radial deviation across the entire 100-orbit trajectory remained within $148\text{ m}$ at $\Delta t = 60\text{ s}$ and $0.58\text{ m}$ at $\Delta t = 15\text{ s}$. At the final endpoint $t = 100T$, the radial deviation is $\Delta r(100T) \approx 7.0 \times 10^{-6}\text{ km} \approx 7.0\text{ mm}$ at $\Delta t = 60\text{ s}$, and $\sim 1.3 \times 10^{-9}\text{ km} \approx 1.3\ \mu\text{m}$ at $\Delta t = 15\text{ s}$. No unverified sub-millimeter claim is made for $\Delta t = 60\text{ s}$.
   * *Phase Error vs. Radial Error*: The dominant contribution to Cartesian position error is along-track phase drift, while radial deviation remains comparatively small ($e_r \approx 83.8\text{ km}$ at $60\text{ s}$, $0.327\text{ km}$ at $15\text{ s}$, consistent with $e_r \approx r_0 |\Delta\theta|$).

### 13.4 M8.2a Scope and Limitations
* **Fixed-dt Only**: All runs use fixed nominal $\Delta t$. Because each step of Yoshida 4 requires 6 force evaluations vs. 4 for RK4, 2 for Verlet, and 1 for Euler, equal-$\Delta t$ comparisons do not normalize for computational work.
* **M8.2b Transition**: Cost-normalized evaluation across matched force evaluation budgets is documented below in Section 14.

---

## 14. Executed Milestone M8.2b: Cost-Normalized Comparison Results (100 Orbits)

### 14.1 Objective and Methodology
* **Core Research Question**: How do Explicit Euler, Velocity Verlet, Classical RK4, and Yoshida 4 perform over a 100-orbit horizon ($t_f = 100 T \approx 555362.432\text{ s}$) when provided with an identical computational budget, measured by the total number of gravitational acceleration evaluations?
* **Computational Cost Proxy**: In the modular, un-cached implementations in this codebase, force evaluations per outer timestep are:
  * Explicit Euler: 1 eval/step
  * Velocity Verlet: 2 evals/step
  * Classical RK4: 4 evals/step
  * Yoshida 4: 6 evals/step (3 substeps $\times$ 2 evals/substep)
  * *Important Note*: This is an implementation-level proxy reflecting force evaluation count, NOT a hardware-dependent wall-clock benchmark.
* **Budget Selection and Step Accounting**:
  To guarantee exact mathematical integer step counts without remainder or step truncation across all 4 solvers, budgets $B$ are selected to be common multiples of $\text{lcm}(1, 2, 4, 6) = 12$:
  $$B \in \{27768, 55536, 111072\}$$
  For each budget $B$ and method:
  $$N_{\text{steps}} = \frac{B}{\text{evals\_per\_step}}, \quad \Delta t = \frac{t_f}{N_{\text{steps}}} = \frac{100 T}{N_{\text{steps}}}$$
  * *Exact Budget Equality*: $\text{actual\_force\_evals} = N_{\text{steps}} \times \text{evals\_per\_step} \equiv B$ exactly.
  * *Uniform Step Size*: All steps within a run are identical; no shortened final step is needed.
  * *Final Time Precision*: Terminal physical timestamp satisfies $|t_{\text{final}} - 100 T| < 10^{-11}\text{ s}$ across all runs.
* **Generated Output Artifacts**:
  * `results/exp08_cost_normalized.csv` (12 data rows + header)
  * `results/exp08_cost_normalized_position_error.png`
  * `results/exp08_cost_normalized_velocity_error.png`
  * `results/exp08_cost_normalized_energy.png`
  * `results/exp08_cost_normalized_angular_momentum.png`
  * `results/exp08_cost_normalized_phase.png`

### 14.2 Cost-Normalized Results Table (100 Orbits)

| Budget $B$ | Method | Steps | $\Delta t$ [s] | Final $e_r$ [km] | Max $e_r$ [km] | Final $\Delta r$ [km] | Final $\Delta\theta$ [rad] | $\max |\delta\varepsilon|$ | Final $\delta\varepsilon$ | $\max |\delta h|$ | Final $\delta h$ |
| :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **27,768** | Euler | 27768 | 20.000 | $3.152 \times 10^{4}$ | $3.922 \times 10^{4}$ | $2.076 \times 10^{4}$ | $-506.84$ | $7.756 \times 10^{-1}$ | $+7.756 \times 10^{-1}$ | $1.099 \times 10^{0}$ | $+1.099 \times 10^{0}$ |
| **27,768** | Verlet | 13884 | 40.000 | $2.883 \times 10^{3}$ | $2.884 \times 10^{3}$ | $0.020$ | $-0.4286$ | $1.045 \times 10^{-6}$ | $+3.034 \times 10^{-9}$ | $5.319 \times 10^{-15}$ | $+9.799 \times 10^{-16}$ |
| **27,768** | RK4 | 6942 | 80.000 | $3.458 \times 10^{2}$ | $3.458 \times 10^{2}$ | $0.720$ | $+0.0510$ | $1.062 \times 10^{-4}$ | $-1.062 \times 10^{-4}$ | $5.312 \times 10^{-5}$ | $-5.312 \times 10^{-5}$ |
| **27,768** | Yoshida 4 | 4628 | 120.001 | $1.344 \times 10^{3}$ | $1.344 \times 10^{3}$ | $0.029$ | $+0.1987$ | $1.225 \times 10^{-7}$ | $+1.485 \times 10^{-9}$ | $5.879 \times 10^{-15}$ | $-3.779 \times 10^{-15}$ |
| **55,536** | Euler | 55536 | 10.000 | $2.431 \times 10^{4}$ | $3.131 \times 10^{4}$ | $1.608 \times 10^{4}$ | $-463.31$ | $7.164 \times 10^{-1}$ | $+7.164 \times 10^{-1}$ | $8.759 \times 10^{-1}$ | $+8.759 \times 10^{-1}$ |
| **55,536** | Verlet | 27768 | 20.000 | $7.264 \times 10^{2}$ | $7.267 \times 10^{2}$ | $3.12 \times 10^{-4}$ | $-0.1072$ | $6.549 \times 10^{-8}$ | $+1.180 \times 10^{-11}$ | $7.979 \times 10^{-15}$ | $+1.680 \times 10^{-15}$ |
| **55,536** | RK4 | 13884 | 40.000 | $1.099 \times 10^{1}$ | $1.099 \times 10^{1}$ | $0.022$ | $+0.0016$ | $3.314 \times 10^{-6}$ | $-3.314 \times 10^{-6}$ | $1.657 \times 10^{-6}$ | $-1.657 \times 10^{-6}$ |
| **55,536** | Yoshida 4 | 9256 | 60.000 | $8.383 \times 10^{1}$ | $8.384 \times 10^{1}$ | $7.31 \times 10^{-6}$ | $+0.0124$ | $4.784 \times 10^{-10}$ | $+2.247 \times 10^{-14}$ | $7.279 \times 10^{-15}$ | $-5.599 \times 10^{-16}$ |
| **111,072** | Euler | 111072 | 5.000 | $1.369 \times 10^{4}$ | $2.617 \times 10^{4}$ | $1.247 \times 10^{4}$ | $-408.91$ | $6.449 \times 10^{-1}$ | $+6.449 \times 10^{-1}$ | $6.778 \times 10^{-1}$ | $+6.778 \times 10^{-1}$ |
| **111,072** | Verlet | 55536 | 10.000 | $1.817 \times 10^{2}$ | $1.818 \times 10^{2}$ | $4.87 \times 10^{-6}$ | $-0.0268$ | $4.095 \times 10^{-9}$ | $+2.127 \times 10^{-14}$ | $2.170 \times 10^{-14}$ | $-1.246 \times 10^{-14}$ |
| **111,072** | RK4 | 27768 | 20.000 | $3.563 \times 10^{-1}$ | $3.563 \times 10^{-1}$ | $7.02 \times 10^{-4}$ | $+5.26 \times 10^{-5}$ | $1.035 \times 10^{-7}$ | $-1.035 \times 10^{-7}$ | $5.177 \times 10^{-8}$ | $-5.177 \times 10^{-8}$ |
| **111,072** | Yoshida 4 | 18512 | 30.000 | $5.235 \times 10^{0}$ | $5.235 \times 10^{0}$ | $1.79 \times 10^{-9}$ | $+7.72 \times 10^{-4}$ | $1.901 \times 10^{-12}$ | $+9.908 \times 10^{-15}$ | $1.694 \times 10^{-14}$ | $+4.759 \times 10^{-15}$ |

*Note: All trajectory-wide maxima (`max_*`) are evaluated strictly over discrete sampled states along each solver's trajectory.*

### 14.3 Multi-Dimensional Trade-Off Analysis (No Ranking/Scores)

Cost normalization reveals clear physical and computational trade-offs depending on the target metric:

1. **Cartesian Trajectory Discrepancy ($e_r$):**
   * Within the scope of the three tested computational budgets, four solvers, and the 100-orbit horizon, Classical RK4 exhibited the smallest final Cartesian position error $e_r(100T)$ among the evaluated methods ($345.8\text{ km}$ vs $1344\text{ km}$ and $2883\text{ km}$ at 27k; $10.99\text{ km}$ vs $83.83\text{ km}$ at 55k; $0.356\text{ km}$ vs $5.24\text{ km}$ at 111k). This is an empirical observation specific to this setup and horizon, not a universal property of RK4 across all dynamical problems.
   * At matched force-evaluation budgets, Yoshida 4 requires 6 force evaluations per outer step compared to 4 for Classical RK4. Consequently, its timestep in this experiment is $1.5\times$ larger than that of RK4 ($\Delta t_{\text{Yoshida}} = 1.5 \Delta t_{\text{RK4}}$). Because along-track phase drift scales with step size ($\Delta\theta \propto t \cdot \Delta t^4$), this larger step size directly contributes to a larger Cartesian position error for Yoshida 4 compared to RK4 on the evaluated 100-orbit horizon, illustrating that symplecticity does not inherently imply smaller Cartesian trajectory error at equal function-evaluation budgets.

2. **Invariant Preservation (Energy and Angular Momentum):**
   * *Angular Momentum*: In the conducted experiments, relative angular-momentum deviation for Velocity Verlet and Yoshida 4 remained at the IEEE 754 float64 numerical-precision level ($\max |\delta h| \le 2.2 \times 10^{-14}$) across all tested budgets. This reflects the central-force symmetry of the problem and discrete update structure, without claiming universal conservation in arbitrary settings.
   * *Energy*: In the conducted 100-orbit fixed-step experiment on the 400 km circular two-body orbit, relative energy deviation for Velocity Verlet and Yoshida 4 remained bounded and oscillatory ($\max |\delta\varepsilon| \sim 1.9 \times 10^{-12}$ for Yoshida 4 at 111k evals); no secular trend was detected over the tested horizon.
   * *Classical RK4*: In the conducted experiments over the 100-orbit horizon, Classical RK4 exhibited a systematic decrease in both energy and angular momentum (final $\delta\varepsilon \approx -1.06 \times 10^{-4}$ at 27k, $-3.31 \times 10^{-6}$ at 55k, and $-1.04 \times 10^{-7}$ at 111k; final $\delta h \approx -5.31 \times 10^{-5}$ at 27k, $-5.18 \times 10^{-8}$ at 111k). This is an empirical observation of the scheme's truncation characteristics under this fixed-step Keplerian configuration.

3. **Radial Shape vs. Along-Track Phase Drift:**
   * For Velocity Verlet and Yoshida 4 across all tested budgets, the dominant contribution to Cartesian position error is associated with along-track phase drift, while radial deviation remains substantially smaller (e.g. at 55k evals, Yoshida 4 has final $\Delta r \approx 7.3\text{ mm}$ vs $e_r \approx 83.8\text{ km}$; Verlet has final $\Delta r \approx 0.31\text{ m}$ vs $e_r \approx 726\text{ km}$).
   * For Classical RK4, a small negative (inward) radial bias was observed on the 100-orbit horizon ($\Delta r(100T) \approx -720\text{ m}$ at 27k, $-22\text{ m}$ at 55k, $-0.70\text{ m}$ at 111k), without implying rapid collapse or spiral into the center.

4. **Explicit Euler:**
   * Euler demonstrated strong unbounded growth of position error, energy deviation, and angular-momentum deviation over the tested 100-orbit horizon across all budgets ($\Delta r > 12000\text{ km}$, $\delta\varepsilon > 0.64$, $\delta h > 0.67$). While Euler is a first-order non-symplectic scheme, receiving $6\times$ more integration steps than Yoshida 4 does not prevent large deviations under the tested conditions.

### 14.4 Scope and Limitations
* **Computational Cost Proxy**: Normalization is based strictly on force evaluations in the current modular, un-cached implementation. It is an implementation-level proxy and NOT a wall-clock benchmark, NOT a universal cost of the methods, and NOT a hardware-independent performance guarantee. No conclusions regarding CPU cache, SIMD, or wall-clock runtime are made.
* **Strict Problem Scope**: Results and conclusions are strictly limited to:
  * 2D planar motion;
  * Unperturbed Newtonian two-body Earth-satellite model;
  * Circular orbit at 400 km altitude ($r_0 = 6778.137\text{ km}$);
  * Fixed-timestep integration;
  * Standard IEEE 754 float64 arithmetic precision;
  * Integration horizon of exactly 100 orbital periods ($100T \approx 555362.432\text{ s}$);
  * Current implementation-level solver routines.
* **No Generalization**: Conclusions cannot be extrapolated to:
  * Eccentric or highly eccentric orbits;
  * Perturbed dynamics ($J_2$, atmospheric drag, third-body luni-solar perturbations);
  * 3D orbital dynamics;
  * Adaptive timestepping or embedded Runge-Kutta / symplectic algorithms;
  * Alternative or optimized implementations of the same methods;
  * Physical hardware-dependent execution speed.
* **Milestone M8.3 Status**: Milestone M8.3 is deferred and was **NOT** executed in this phase.
