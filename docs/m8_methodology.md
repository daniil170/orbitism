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
   * *Hypothesis*: Over long durations ($100T$), non-symplectic integrators (Euler and RK4) will exhibit secular drift in specific orbital energy ($\delta \varepsilon(t) \propto t$), whereas symplectic integrators (Velocity Verlet and Yoshida 4) will constrain energy error within bounded oscillatory bands ($\max |\delta \varepsilon(t)| \le C \Delta t^p$) without secular growth.

4. **Phase Error vs. Radial Drift Decomposition**:
   * *Hypothesis*: Across all stable integrators, total Euclidean position error over multi-orbit horizons ($t \gg T$) is dominated by along-track phase error ($e_T(t) \propto t \Delta t^p$), while radial drift ($e_R(t)$) remains bounded for symplectic methods but grows for non-symplectic methods.

5. **Angular Momentum Conservation**:
   * *Hypothesis*: Because Earth two-body gravity is an exact central force ($\mathbf{r} \times \mathbf{F} = 0$), angular momentum error $\delta h(t)$ is governed primarily by Cartesian roundoff accumulation rather than discretization truncation in all symmetric/consistent formulations.

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
| **Short-Term** | $1T \approx 5553.62\text{ s}$ | $1.0$ | Discretization error convergence rates ($p$), local truncation scaling, clean asymptotic slope measurement. |
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
* **Numerical Stability**:
  * For circular LEO ($T \approx 5554\text{ s}$), $\Delta t = 60\text{ s}$ corresponds to $\approx 92.5$ steps per orbit.
  * All four integrators remain stable (non-exploding within 1 orbit) at $\Delta t = 60\text{ s}$.
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

Scientific analysis must strictly distinguish between three distinct physical/numerical regimes:

1. **Mathematical Conservation Property**:
   * *Symplectic Integrators (Verlet, Yoshida 4)*: Preserve the symplectic 2-form $d\mathbf{q} \wedge d\mathbf{p}$. They exactly conserve a perturbed Hamiltonian $\tilde{H} = H + \Delta t^p H_p + \dots$. Energy error is mathematically bounded for exponentially long times.
   * *Non-Symplectic Integrators (Euler, RK4)*: Do not preserve phase space volume or symplectic 2-form. Truncation error drives artificial dissipation or artificial energy pumping.
2. **Floating-Point Roundoff Floor**:
   * For double-precision IEEE 754 arithmetic, machine epsilon is $\epsilon_{\text{mach}} \approx 2.22 \times 10^{-16}$.
   * Angular momentum conservation in a central force field is preserved to near machine precision ($\sim 10^{-15} - 10^{-14}$) subject only to Cartesian roundoff.
   * For Yoshida 4 at small $\Delta t$, energy error reaches a roundoff noise floor of $\sim 10^{-14} - 10^{-13}$.
3. **Observed Experimental Character**:
   * Classify energy error time evolution into:
     * **Secular Drift**: $|\delta \varepsilon(t)| \propto t^\alpha$ ($\alpha \ge 1$, typical of Euler and RK4).
     * **Bounded Oscillation**: $|\delta \varepsilon(t)| \le M(\Delta t)$, zero mean secular slope (typical of Verlet and Yoshida 4).

---

## 8. Cost Normalization Model

### 8.1 Force Evaluation Proxy

Wall-clock execution time varies based on CPU architecture, background operating system tasks, and compiler optimizations. Therefore, **the number of gravitational acceleration (force) evaluations** is adopted as the primary, hardware-independent metric of computational cost.

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
`method,dt,num_steps,force_evals,e_r,e_v,order_p_r,ratio_r,max_de,max_dh,final_dr,final_dtheta`

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
