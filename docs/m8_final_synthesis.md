# Milestone M8: Final Scientific Synthesis

## Systematic Comparison of Numerical Integrators in Two-Body Orbital Mechanics

---

### Executive Summary

Milestone M8 completes the systematic numerical and geometric comparison of four foundational integration algorithms implemented in OrbitSim:
1. **Explicit Euler** (1st-order, single-step, non-symplectic);
2. **Velocity Verlet** (2nd-order, two-stage, symplectic, time-reversible);
3. **Classical Runge-Kutta 4 (RK4)** (4th-order, four-stage, non-symplectic);
4. **Yoshida 4** (4th-order, three-stage Velocity Verlet composition, symplectic, symmetric).

This benchmark was conducted across four distinct experimental investigations:
* **M8.1**: Analytical Reference + Discretization Convergence ($1T$ horizon, $\Delta t \in [1.875, 60.0]\text{ s}$);
* **M8.2a**: Long-Term Fixed-Timestep Comparison ($100T$ horizon, $\Delta t \in [15.0, 60.0]\text{ s}$);
* **M8.2b**: Cost-Normalized Comparison ($100T$ horizon, matched force-evaluation budgets);
* **M8.3**: Orbital Period Error and Phase Drift Analysis ($100T$ horizon, prograde positive x-axis event detection via cubic Hermite root finding).

All simulations were benchmarked against the exact closed-form analytical Keplerian reference solution on an unperturbed geocentric circular orbit at an altitude of $400\text{ km}$.

---

## 1. Scientific Objective

The overarching research objective of Milestone M8 was formulated as follows:

> *"How do Explicit Euler, Velocity Verlet, Classical RK4, and Yoshida 4 compare under both matched-timestep ($\Delta t$) and matched-computational-budget (force evaluations) paradigms, and how does the choice of benchmark metric (instantaneous Cartesian truncation error vs. geometric invariant preservation vs. along-track phase drift) alter the scientific assessment of solver suitability for orbital mechanics?"*

### Core Guiding Principles

1. **Absence of a Pre-Ordained "Winner"**: M8 is not a competitive benchmark designed to declare a single "best" algorithm. In orbital mechanics, different applications present fundamentally irreconcilable requirements (e.g., short-term rendezvous navigation requiring minimal instantaneous Cartesian error vs. multi-year ephemeris propagation requiring strict energy and phase boundedness).
2. **Exact Analytical Baseline**: Numerical trajectories are evaluated strictly against the exact closed-form Keplerian circular solution, completely eliminating reliance on surrogate numerical approximations.
3. **Multi-Faceted Assessment**: Performance is evaluated across four distinct physical axes:
   * Local and asymptotic discretization truncation accuracy ($e_r, e_v$);
   * First-integral preservation (orbital energy $\varepsilon$ and specific angular momentum $h$);
   * Geometric orbital stability (radial deviation $\Delta r$ vs. along-track transverse error $e_T$);
   * Along-track phase drift ($\Delta\theta(t)$) and orbital period error ($\delta T$).
4. **Dual-Axis Evaluation**: Every solver is characterized both at **equal step size** $\Delta t$ (mathematical map properties) and at **equal computational effort** (cost-normalized force evaluation budgets).

---

## 2. Experimental Scope and Parameterization

All conclusions established in Milestone M8 apply strictly within the following experimental domain:

| Parameter | Value | Scientific Rationale |
| :--- | :--- | :--- |
| **Physical Model** | Planar 2D Unperturbed Two-Body Gravity ($\ddot{\mathbf{r}} = -\mu \mathbf{r} / r^3$) | Isolates algorithmic discretization effects from environmental perturbations. |
| **Central Body** | Earth ($\mu = 398600.435\text{ km}^3/\text{s}^2$, $R_E = 6378.137\text{ km}$) | Canonical astrodynamics standard parameters. |
| **Orbit Geometry** | Circular ($e = 0.0$), Altitude $h = 400.0\text{ km}$ ($r_0 = 6778.137\text{ km}$) | Provides exact analytical closed-form reference for all $t \ge 0$. |
| **Initial State** | $\mathbf{r}_0 = [6778.137, 0.0]^T\text{ km}$, $\mathbf{v}_0 = [0.0, \sqrt{\mu/r_0}]^T \approx [0.0, 7.668631]^T\text{ km/s}$ | Prograde circular motion starting on positive x-axis. |
| **Analytical Period** | $T_{\text{exact}} = 2\pi \sqrt{r_0^3 / \mu} \approx 5553.624318623783\text{ s}$ ($\approx 92.5604\text{ min}$) | Primary fundamental timescale for orbital revolutions. |
| **Mean Motion** | $n = \sqrt{\mu / r_0^3} \approx 1.131370197 \times 10^{-3}\text{ rad/s}$ | Analytical angular velocity. |
| **Simulation Horizons** | Short-term: $1T$ ($\approx 5553.62\text{ s}$); Long-term: $100T$ ($\approx 555362.432\text{ s} \approx 6.43\text{ days}$) | Short-term isolates local order; long-term reveals secular drift. |
| **Floating-Point Standard** | IEEE 754 Double Precision (`float64`, 53-bit mantissa) | Machine epsilon $\epsilon_{\text{mach}} \approx 2.22 \times 10^{-16}$. |
| **Force Evaluations / Step** | Euler: 1, Velocity Verlet: 2, Classical RK4: 4, Yoshida 4: 6 | Modular un-cached implementation-level cost proxy. |

---

## 3. Overview of Constituent Experiments

```
                                  Milestone M8
                      Systematic Integrator Comparison
                                     │
         ┌───────────────────────────┼───────────────────────────┐
         ▼                           ▼                           ▼
       M8.1                        M8.2                        M8.3
Convergence & Reference     Long-Term & Cost-Normalized   Period Error & Phase Drift
(1 Orbit, dt ∈ [1.875, 60]s) (100 Orbits, Matched Budgets)(100 Orbits, Hermite Events)
         │                           │                           │
         │                  ┌────────┴────────┐                  │
         │                  ▼                 ▼                  │
         │                M8.2a             M8.2b                │
         │             Fixed-dt        Cost-Normalized           │
         │           (100 Orbits)      (Matched Evals)           │
         │                  │                 │                  │
         └──────────────────┴────────┬────────┴──────────────────┘
                                     ▼
                            Final M8 Synthesis
```

### 3.1 M8.1: Analytical Reference + Discretization Convergence
* **Focus**: Convergence order $p$ and truncation error over one complete orbital period ($1T$).
* **Timestep Grid**: $\Delta t \in \{60.0, 30.0, 15.0, 7.5, 3.75, 1.875\}\text{ s}$.
* **Key Metrics**: Terminal position error $e_r(T)$, velocity error $e_v(T)$, empirical order $p = \log_2(E(2h)/E(h))$, maximum energy deviation $\max |\delta\varepsilon|$, maximum angular momentum deviation $\max |\delta h|$.
* **Execution**: 24 runs (4 solvers $\times$ 6 timesteps), strictly hitting $t_f = T_{\text{exact}}$ via shortened final step where necessary.

### 3.2 M8.2a: Long-Term Fixed-Timestep Comparison
* **Focus**: Secular drift vs. bounded oscillations over a 100-orbit horizon ($100T \approx 6.43\text{ days}$).
* **Timestep Grid**: $\Delta t \in \{60.0, 30.0, 15.0\}\text{ s}$.
* **Checkpoints**: $\mathcal{T} = \{1.0, 5.0, 10.0, 20.0, 50.0, 100.0\}T$.
* **Key Metrics**: Cartesian error $e_r(t)$, radial deviation $\Delta r(t) = r(t) - r_0$, transverse error $e_T(t)$, phase error $\Delta\theta(t)$, energy error $\delta\varepsilon(t)$, angular momentum error $\delta h(t)$.
* **Execution**: 12 runs (4 solvers $\times$ 3 timesteps), 72 logged checkpoint records.

### 3.3 M8.2b: Cost-Normalized Comparison
* **Focus**: Performance comparison under equal computational effort, defined by fixed budgets of gravitational acceleration (force) evaluations.
* **Evaluation Budgets**: $B \in \{27768, 55536, 111072\}$ force evaluations over $100T$.
* **Matched Timesteps**: Derived from $B = N_{\text{steps}} \times \mathbf{C}_{\text{step}}$, where $\mathbf{C}_{\text{step}} = [1, 2, 4, 6]$ for [Euler, Verlet, RK4, Yoshida 4]:
  * Budget $27768$: Euler $\Delta t = 20\text{ s}$, Verlet $40\text{ s}$, RK4 $80\text{ s}$, Yoshida 4 $120\text{ s}$.
  * Budget $55536$: Euler $\Delta t = 10\text{ s}$, Verlet $20\text{ s}$, RK4 $40\text{ s}$, Yoshida 4 $60\text{ s}$.
  * Budget $111072$: Euler $\Delta t = 5\text{ s}$, Verlet $10\text{ s}$, RK4 $20\text{ s}$, Yoshida 4 $30\text{ s}$.
* **Execution**: 12 runs (4 solvers $\times$ 3 budgets), integer step counts, identical terminal physical duration.

### 3.4 M8.3: Orbital Period Error and Phase Drift Analysis
* **Focus**: Sub-step event detection of orbital revolutions, empirical scaling of period error $\delta T$, phase drift rate $|\dot{\theta}_{\text{drift}}|$, and decomposition of Cartesian error into chord error vs. radial drift.
* **Event Definition**: Crossing of the positive x-axis in the prograde direction ($y_i \le 0 < y_{i+1}$, $x > 0$, $v_y > 0$).
* **Interpolation**: Cubic Hermite polynomial root finding with bracket-preserving bisection fallback ($|P(s)| < 10^{-12}$).
* **Timestep Grid**: $\Delta t \in \{60.0, 30.0, 15.0, 7.5\}\text{ s}$ over $\ge 100.5T$ with strictly constant step size (no shortened step).
* **Checkpoints**: Revolutions $m \in \{1, 5, 10, 20, 50, 100\}$.
* **Execution**: 16 runs (4 solvers $\times$ 4 timesteps), 96 logged checkpoint records.

---

## 4. Key Numerical Observations

### 4.1 Discretization Convergence Orders

Empirical convergence orders $p = \log_2(E(2h)/E(h))$ obtained across the experimental series:

| Solver | Formal Order | Observed Order (M8.1, $1T$) | Observed Period Order $p_T$ (M8.3, $100T$) | Observed Phase Order $p_\theta$ (M8.3, $100T$) | Characteristic Behavior |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Explicit Euler** | 1 | $0.64 \to 0.97$ | N/A (diverged) | $0.08 \to 0.15$ | Approaches 1st-order on short arcs; sub-linear in multi-orbit phase due to exponential radial expansion. |
| **Velocity Verlet** | 2 | $2.00$ | $2.00$ | $2.00$ | Strictly quadratic scaling across all timestep intervals and horizons. |
| **Classical RK4** | 4 | $4.02 \to 4.22$ | $4.87 \to 4.97$ | $4.87 \to 4.97$ | 4th-order on short arcs; exhibits super-convergence on circular orbit due to leading symmetric truncation cancellation. |
| **Yoshida 4** | 4 | $4.00$ | $4.00$ | $4.00$ | Strictly quartic scaling across all timestep intervals, metrics, and horizons. |

### 4.2 Physical Invariant Preservation

#### Specific Angular Momentum ($\delta h = (h - h_0)/|h_0|$)
* **Velocity Verlet & Yoshida 4**: Specific angular momentum is preserved to machine precision ($|\delta h| \le 2.2 \times 10^{-14}$) across all timesteps, budgets, and horizons. This exact discrete conservation is a direct mathematical consequence of central forces ($\mathbf{r} \times \mathbf{a} = \mathbf{0}$) combined with the discrete symplectic kick-drift update structure.
* **Classical RK4**: Exhibits systematic, irreversible secular dissipation of angular momentum. At $100T$ with $\Delta t = 60\text{ s}$, $\delta h \approx -1.26 \times 10^{-5}$; at $\Delta t = 15\text{ s}$, $\delta h \approx -1.23 \times 10^{-8}$.
* **Explicit Euler**: Exhibits rapid secular growth of angular momentum ($\delta h > +0.67$ at $100T$ across all tested timesteps).

#### Specific Orbital Energy ($\delta\varepsilon = (\varepsilon - \varepsilon_0)/|\varepsilon_0|$)
* **Velocity Verlet & Yoshida 4**: Energy error remains strictly bounded within oscillatory envelopes with zero detectable secular trend over the entire 100-orbit horizon. In accordance with backward error analysis for autonomous Hamiltonian systems, the numerical solution traces the level sets of a nearby modified Hamiltonian $\tilde{H} = H + \mathcal{O}(\Delta t^p)$.
  * Velocity Verlet: Maximum deviation $|\delta\varepsilon| \sim 5.3 \times 10^{-6}$ at $\Delta t = 60\text{ s}$; $\sim 8.1 \times 10^{-11}$ at $\Delta t = 3.75\text{ s}$.
  * Yoshida 4: Maximum deviation $|\delta\varepsilon| \sim 4.8 \times 10^{-10}$ at $\Delta t = 60\text{ s}$; $\sim 1.9 \times 10^{-12}$ at $\Delta t = 30\text{ s}$.
* **Classical RK4**: Exhibits secular, monotonic energy loss (numerical dissipation). Over 100 orbits at $\Delta t = 60\text{ s}$, energy error decreases steadily to $\delta\varepsilon \approx -2.52 \times 10^{-5}$ (Budget 27k: $\delta\varepsilon \approx -1.06 \times 10^{-4}$ at $\Delta t = 80\text{ s}$).
* **Explicit Euler**: Exhibits massive secular energy gain ($\delta\varepsilon > +0.64$ at $100T$ across all tested timesteps).

### 4.3 Orbital Geometry: Radial Deviation vs. Along-Track Phase Drift

A central finding of Milestone M8 is the quantitative decomposition of Cartesian position error $e_r(t) = \|\mathbf{r}_{\text{num}}(t) - \mathbf{r}_{\text{exact}}(t)\|_2$:

$$e_r(t)^2 \approx e_R(t)^2 + e_T(t)^2$$

Where $e_R$ is the radial error and $e_T$ is the along-track transverse error. Under near-circular geometry ($r \approx r_0$), the along-track error is related to phase error $\Delta\theta = \theta_{\text{num}} - n t$ by the geometric chord formula:

$$e_{\text{chord}}(\Delta\theta) = 2 r_0 \left|\sin\left(\frac{\Delta\theta}{2}\right)\right|$$

Empirical results at the 100-orbit checkpoint ($t \approx 555362.432\text{ s}$, $\Delta t = 60.0\text{ s}$):

| Solver | Cartesian Error $e_r$ (km) | Geometric Chord Error $e_{\text{chord}}$ (km) | Ratio $e_{\text{chord}} / e_r$ | Radial Error $|\Delta r|$ (km) | Transverse Error $|e_T|$ (km) | Radial Dominance Ratio $|\Delta r| / e_r$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Euler** | $61765.052$ | $8658.998$ | **0.140** | $55873.786$ | $26268.490$ | **0.905** |
| **Velocity Verlet** | $6280.591$ | $6280.483$ | **0.99998** | $0.234$ | $5566.000$ | **$3.7 \times 10^{-5}$** |
| **Classical RK4** | $82.498$ | $82.499$ | **1.00001** | $0.171$ | $82.496$ | **$2.1 \times 10^{-3}$** |
| **Yoshida 4** | $83.832$ | $83.832$ | **1.00000** | $0.00001$ ($1\text{ cm}$) | $83.831$ | **$1.2 \times 10^{-7}$** |

#### Scientific Deductions from Geometry Decomposition:
1. **Geometric Fidelity of Symplectic and Conservative Schemes**: For Velocity Verlet, Yoshida 4, and Classical RK4, the numerical orbit remains almost perfectly circular with minimal radial drift ($|\Delta r| \le 0.23\text{ km}$). Total Cartesian position error at $100T$ is virtually $100\%$ attributable to along-track phase drift ($e_{\text{chord}} / e_r \approx 1.0$).
2. **Breakdown for Explicit Euler**: For Explicit Euler, radial drift accounts for over $90\%$ of the error magnitude ($|\Delta r| > 55000\text{ km}$). The circular chord model completely breaks down ($e_{\text{chord}} / e_r = 0.14$).

### 4.4 Orbital Period Error and Phase Drift Direction

At each positive x-axis prograde crossing event $t^*_m$, the true unwrapped angle is $\theta_{\text{num}}(t^*_m) \equiv 2\pi m$. The relative cumulative period error $\delta T_{\text{cum}} = (t^*_m/m - T_{\text{exact}})/T_{\text{exact}}$ satisfies the mathematical identity:

$$\delta T_{\text{cum}}(m) \equiv -\frac{\Delta\theta(t^*_m)}{2\pi m}$$

Results at 100 orbits ($m = 100$, $\Delta t = 60.0\text{ s}$):
* **Velocity Verlet**: $\delta T_{\text{cum}} \approx +1.532 \times 10^{-3}$ ($\Delta T \approx +8.51\text{ s/orbit}$). Because the numerical period is slightly longer than analytical, the numerical satellite **lags** behind the reference position ($\Delta\theta(100T) \approx -0.963\text{ rad} \approx -55.2^\circ$).
* **Classical RK4**: $\delta T_{\text{cum}} \approx -1.937 \times 10^{-5}$ ($\Delta T \approx -0.108\text{ s/orbit}$). The numerical period is slightly shorter than analytical; the satellite **leads** ($\Delta\theta(100T) \approx +0.0122\text{ rad} \approx +0.70^\circ$).
* **Yoshida 4**: $\delta T_{\text{cum}} \approx -1.968 \times 10^{-5}$ ($\Delta T \approx -0.109\text{ s/orbit}$). The numerical period is slightly shorter than analytical; the satellite **leads** ($\Delta\theta(100T) \approx +0.0124\text{ rad} \approx +0.71^\circ$).
* **Explicit Euler**: Due to outward spiraling, the orbital period expands continuously ($T \propto r^{3/2}$). In $t = 100 T_{\text{exact}}$, Euler completes only 10 revolutions ($\Delta t = 60\text{ s}$) to 29 revolutions ($\Delta t = 7.5\text{ s}$); the 100th crossing is never reached.

---

## 5. Important Conceptual Distinctions

Milestone M8 establishes several critical distinctions that are essential for the rigorous interpretation of numerical algorithms in astrodynamics:

### 5.1 Cartesian Truncation Error $\ne$ Invariant Preservation

A numerical scheme may achieve low Cartesian position error over finite time while systematically violating Hamiltonian invariant structure, or vice versa:
* **Case 1 (Low Cartesian Error with Systematic Dissipation)**: At Budget $111072$ ($100T$), Classical RK4 achieves a Cartesian position error of only $e_r \approx 0.356\text{ km}$, yet exhibits continuous secular dissipation of energy ($\delta\varepsilon \approx -1.04 \times 10^{-7}$) and angular momentum ($\delta h \approx -5.18 \times 10^{-8}$).
* **Case 2 (Strict Invariant Boundedness with Large Cartesian Error)**: At $\Delta t = 60\text{ s}$ ($100T$), Velocity Verlet maintains energy error strictly bounded within $|\delta\varepsilon| \le 5.2 \times 10^{-6}$ and angular momentum to machine precision ($|\delta h| < 10^{-14}$), yet its Cartesian position error grows to $e_r \approx 6280.6\text{ km}$ due to along-track phase drift.

### 5.2 Invariant Preservation $\ne$ Minimal Phase Drift

Symplectic integrators are guaranteed by backward error analysis to preserve a modified Hamiltonian $\tilde{H} = H + \Delta t^p H_p + \dots$ over exponentially long times (for analytic vector fields). This guarantees:
1. Bounded energy oscillations;
2. Conservation of orbital shape (radial boundedness, no artificial spiral).

**Crucially, symplecticity does NOT imply zero along-track phase drift.** The discrete map modifies the fundamental frequency of the system:
$$\tilde{\omega} = \omega_0 \left(1 + C_p \Delta t^p + \mathcal{O}(\Delta t^{p+2})\right)$$
This frequency perturbation generates a secular phase error that accumulates linearly in time:
$$\Delta\theta(t) \approx (\tilde{\omega} - \omega_0) t \propto \Delta t^p \cdot t$$
Both Velocity Verlet ($p=2$) and Yoshida 4 ($p=4$) exhibit strictly linear secular phase accumulation in time, despite exact geometric stability.

### 5.3 Dissipative High-Order vs. Symplectic High-Order Comparison

Comparing Classical RK4 (4th-order, non-symplectic) against Yoshida 4 (4th-order, symplectic):
* **At Equal Timestep ($\Delta t = 60\text{ s}$)**: Over 100 orbits, RK4 yields Cartesian error $e_r \approx 82.5\text{ km}$ and phase error $\Delta\theta \approx +0.0122\text{ rad}$, which are slightly smaller than Yoshida 4 ($e_r \approx 83.8\text{ km}$, $\Delta\theta \approx +0.0124\text{ rad}$), due to smaller leading truncation error coefficients in the standard RK4 Butcher tableau.
* **At Equal Computational Budget (Budget 55k, $\approx 55536$ evaluations)**:
  * RK4 operates at $\Delta t = 40.0\text{ s}$ ($4\text{ evals/step}$) and achieves $e_r \approx 10.99\text{ km}$.
  * Yoshida 4 operates at $\Delta t = 60.0\text{ s}$ ($6\text{ evals/step}$) and achieves $e_r \approx 83.83\text{ km}$.
  * However, Yoshida 4 bounds energy deviation to $|\delta\varepsilon| \le 4.8 \times 10^{-10}$ and preserves $h$ to $7.3 \times 10^{-15}$, whereas RK4 suffers secular dissipation ($\delta\varepsilon \approx -3.3 \times 10^{-6}$, $\delta h \approx -1.7 \times 10^{-6}$).
* **Long-Term Extrapolation**: On horizons spanning thousands of orbits, RK4's secular dissipation causes cumulative spiral-in and ephemeris degradation, whereas Yoshida 4 remains bounded on the modified Hamiltonian torus. On shorter horizons with high force-evaluation cost constraints, RK4 provides smaller point-to-point Cartesian error.

---

## 6. Interpretation of Cost-Normalized vs. Fixed-$\Delta t$ Paradigms

Milestone M8 evaluated solvers across two complementary axes that answer fundamentally different scientific questions:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              EVALUATION AXIS                                │
├──────────────────────────────────────┬──────────────────────────────────────┤
│        Fixed Timestep (Δt)           │       Cost-Normalized Budget (B)     │
├──────────────────────────────────────┼──────────────────────────────────────┤
│ Scientific Question:                 │ Scientific Question:                 │
│ "What are the mathematical error     │ "What is the achievable physical     │
│ scaling properties of the discrete   │ fidelity for a given computational   │
│ map per unit integration step?"      │ investment in force evaluations?"    │
├──────────────────────────────────────┼──────────────────────────────────────┤
│ Primary Application:                 │ Primary Application:                 │
│ Verifying formal convergence orders, │ Practical algorithm selection under  │
│ local truncation error, and          │ compute budgets (e.g. real-time,     │
│ asymptotic behavior.                 │ embedded, multi-body simulations).   │
└──────────────────────────────────────┴──────────────────────────────────────┘
```

### Force Evaluations as an Implementation Proxy

* In the current codebase, computational cost is quantified by the number of gravitational acceleration evaluations:
  $$\mathbf{C}_{\text{step}} = [1, 2, 4, 6] \quad \text{for } [\text{Euler, Verlet, RK4, Yoshida 4}]$$
* **Strict Qualification**: This metric is an **implementation-level proxy** for this decoupled, modular codebase. It is:
  * **NOT equivalent to wall-clock execution time**;
  * **NOT a universal model of computational complexity**;
  * **NOT a hardware-independent performance guarantee**.
* In specialized production implementations, Velocity Verlet can be cached to 1 evaluation/step and Yoshida 4 to 3 evaluations/step by storing accelerations across step boundaries. In accordance with OrbitSim guidelines, cost accounting reflects the current decoupled modular architecture.

---

## 7. Comprehensive Summary of Evaluated Solvers

The table below provides a strictly descriptive synthesis of the four solvers benchmarked in Milestone M8. No columns are included for "winner", "rank", or "recommendation".

| Property / Metric | Explicit Euler | Velocity Verlet | Classical RK4 | Yoshida 4 |
| :--- | :--- | :--- | :--- | :--- |
| **Formal Order ($p$)** | 1 | 2 | 4 | 4 |
| **Symplectic Structure** | Non-symplectic | Symplectic | Non-symplectic | Symplectic |
| **Time Reversibility** | No | Yes | No | Yes |
| **Evaluations / Step** | 1 | 2 | 4 | 6 |
| **Observed $1T$ Order** | $p \approx 0.94 - 0.97$ | $p = 2.00$ | $p \approx 4.02 - 4.22$ | $p = 4.00$ |
| **Observed Period Order** | Sub-linear ($p \le 0.15$) | $p = 2.00$ | $p \approx 4.87 - 4.97$ | $p = 4.00$ |
| **Energy Behavior ($100T$)** | Secular growth ($\delta\varepsilon > +0.64$) | Bounded oscillations ($|\delta\varepsilon| \le 5.3 \times 10^{-6}$) | Secular dissipation ($\delta\varepsilon \sim -10^{-4} \dots -10^{-7}$) | Bounded oscillations ($|\delta\varepsilon| \le 1.9 \times 10^{-12}$) |
| **Angular Momentum ($100T$)** | Secular growth ($\delta h > +0.67$) | Machine precision ($|\delta h| \le 2.2 \times 10^{-14}$) | Secular dissipation ($\delta h \sim -10^{-5} \dots -10^{-8}$) | Machine precision ($|\delta h| \le 2.2 \times 10^{-14}$) |
| **Phase Error Direction** | Extreme Lag ($\Delta\theta \ll 0$) | Lag ($\Delta\theta < 0$, $\delta T > 0$) | Lead ($\Delta\theta > 0$, $\delta T < 0$) | Lead ($\Delta\theta > 0$, $\delta T < 0$) |
| **Radial Stability ($100T$)** | Severe spiral outward ($\Delta r > 15000\text{ km}$) | Strictly bounded ($|\Delta r| \le 0.23\text{ km}$) | Small inward bias ($|\Delta r| \le 0.17\text{ km}$) | Strictly bounded ($|\Delta r| \le 1\text{ cm}$) |
| **Cartesian Error ($100T$, $\Delta t = 60\text{ s}$)** | $61765.1\text{ km}$ | $6280.6\text{ km}$ | $82.5\text{ km}$ | $83.8\text{ km}$ |
| **Cartesian Error ($100T$, Budget 55k)** | $24314.1\text{ km}$ ($\Delta t = 10\text{ s}$) | $726.4\text{ km}$ ($\Delta t = 20\text{ s}$) | $10.99\text{ km}$ ($\Delta t = 40\text{ s}$) | $83.83\text{ km}$ ($\Delta t = 60\text{ s}$) |
| **Observed Operational Characteristics** | Unsuitable for orbital mechanics; rapid secular departure from physical state. | Robust structural stability; machine-precision angular momentum; moderate phase drift. | High point-to-point Cartesian precision over short/medium horizons; irreversible secular dissipation. | Simultaneous 4th-order Cartesian accuracy, strict energy boundedness, and machine-precision angular momentum. |

---

## 8. Substantial Scientific Limitations

The conclusions established in Milestone M8 are bound by the following structural limitations:

1. **Two-Dimensional Planar Dynamics**: Benchmarked exclusively in 2D Cartesian coordinates ($x, y$). Out-of-plane orbital dynamics, orbital inclination, nodal precession, and 3D perturbations were not modeled.
2. **Unperturbed Central-Force Field**: Dynamics were governed strictly by idealized point-mass Newtonian gravity ($\ddot{\mathbf{r}} = -\mu \mathbf{r} / r^3$). No geopotential harmonics ($J_2, J_3, J_4$), atmospheric drag, solar radiation pressure, or third-body gravitational perturbations (Moon, Sun) were included.
3. **Single Orbit Baseline (Circular, 400 km Altitude)**: Evaluated exclusively on a circular orbit ($e = 0.0$) at $r_0 = 6778.137\text{ km}$. Highly eccentric orbits ($e \to 1$), hyperbolic escape trajectories, and parabolic orbits were not investigated. On eccentric orbits, varying velocity along the orbit significantly alters local truncation error distributions.
4. **Strictly Constant Timestepping**: Solvers operated with fixed $\Delta t$. Adaptive step-size control (e.g. embedded Runge-Kutta-Fehlberg or adaptive symplectic schemes) was not evaluated.
5. **Finite Time Horizon (100 Orbital Periods)**: Long-term comparisons spanned $100T \approx 6.43\text{ days}$. Behavior over multi-year, decadal, or secular resonance timescales was not simulated.
6. **Finite Timestep Exploration**: Step sizes were restricted to $\Delta t \in [1.875, 120.0]\text{ s}$. Asymptotic behavior as $\Delta t \to 0$ or stability limits near the linear stability boundary ($\Delta t \to \Delta t_{\text{crit}}$) were not fully mapped.
7. **Arithmetic Precision Limits**: Computations utilized standard IEEE 754 double precision (`float64`). Truncation errors below $\sim 10^{-14}$ encounter the floating-point roundoff floor.
8. **Implementation-Specific Cost Accounting**: Force evaluation counts reflect the current un-cached modular software architecture and do not capture hardware-specific cache locality, memory bandwidth, or instruction-level parallelism.

---

## 9. What Milestone M8 Does NOT Establish

To maintain rigorous scientific boundaries, it is explicitly declared that Milestone M8 does **NOT** establish:

1. **A Universal "Best" Integrator**: M8 does not declare any method universally superior. Solver selection is strictly application-dependent.
2. **A Global Solver Ranking or Hierarchy**: No scalar score or tier list is valid. A method that excels in short-term landing guidance (e.g. RK4) may be inappropriate for multi-century planetary ephemerides, where symplectic schemes are required.
3. **Universal Stability Guarantees**: Demonstrating bounded energy over 100 orbits does not constitute a mathematical proof of stability for all initial conditions, arbitrary step sizes, or indefinite time horizons.
4. **Wall-Clock Performance Superiority**: Force evaluation counts do not directly translate into wall-clock speed across different hardware architectures, compilers, or programming languages.
5. **That Symplectic Methods Always Outperform Runge-Kutta Schemes**: On finite horizons, RK4 demonstrated smaller Cartesian error and phase drift than Velocity Verlet across all tested step sizes, and smaller Cartesian error than Yoshida 4 under equal force evaluation budgets.
6. **That Runge-Kutta Methods Always Outperform Symplectic Methods**: On long horizons where invariant preservation is paramount, RK4's secular energy and angular momentum dissipation causes irreversible structural degradation that symplectic methods avoid.
7. **That Observed Convergence Orders Generalize to Arbitrary Systems**: The empirical convergence rates observed here reflect the specific smoothness and central symmetry of unperturbed circular Keplerian motion and cannot be assumed to hold unconditionally in stiff, chaotic, or non-smooth dynamical systems.

---

### Conclusion and Transition to Milestone M9

Milestone M8 formally establishes the empirical foundations of fixed-step numerical integration in OrbitSim. With discretization error scaling, invariant preservation characteristics, and phase drift mechanisms rigorously quantified and documented, the project is mathematically prepared to transition to subsequent developmental milestones (including elliptical orbits, perturbation modeling, and adaptive integration techniques).
