# OrbitSim

A scientific computing and orbital simulation engine designed for satellite trajectory modeling and rigorous investigation of numerical integration methods, precision, stability, and computational performance.

---

## Goal

OrbitSim investigates the numerical mechanics of orbital dynamics:
- Error propagation and long-term numerical stability across integration schemes.
- Conservation of physical invariants (total energy, angular momentum) in discrete time.
- Accuracy vs. performance trade-offs between low-order, high-order, and symplectic solvers.
- Scientific correctness and reproducibility over decorative graphics.

---

## Current Model

In the current development phase, OrbitSim models **2D planar satellite motion** around Earth:
- Classical two-body gravitational model.
- Earth positioned at the origin $(0, 0)$ as an unmoving central body.
- Satellite treated as a point mass.
- Unperturbed central gravity field (no atmospheric drag, no lunar/solar perturbations, no $J_2$).
- Circular Low Earth Orbit (LEO) at an altitude of $400\text{ km}$ ($r_0 = 6778.137\text{ km}$).

---

## Physics

The equations of motion are governed by Newton's law of universal gravitation:

$$\ddot{\mathbf{r}} = -\frac{\mu}{r^3}\mathbf{r}$$

* Gravitational parameter of Earth: $\mu = 398600.435\text{ km}^3/\text{s}^2$
* Equatorial radius of Earth: $R_E = 6378.137\text{ km}$
* State vector: $\mathbf{s} = [x, y, v_x, v_y]^T$ (km, km/s)
* Derivative: $\frac{d\mathbf{s}}{dt} = [v_x, v_y, a_x, a_y]^T$

---

## Numerical Methods

* **Explicit Euler (Forward Euler)**: First-order numerical integrator ($p = 1$) serving as baseline solver.
* **Classical Runge-Kutta (RK4)**: Fourth-order integrator ($p = 4$) providing high-order benchmark accuracy.
* **Velocity Verlet**: Second-order symplectic integrator ($p = 2$) preserving phase-space geometry and invariants.
* **Yoshida 4th-Order Symplectic**: Fourth-order symplectic integrator ($p = 4$) via symmetric composition of Velocity Verlet steps, combining high-order accuracy with exact phase-space symplecticity.

All integrators are strictly decoupled from physical equations and operate generically on state derivative callables.

---

## Project Structure

```text
orbitsim/
│
├── src/
│   ├── __init__.py
│   ├── physics/
│   │   ├── __init__.py
│   │   └── gravity.py           # Gravitational acceleration, derivatives, initial conditions
│   │
│   ├── integrators/
│   │   ├── __init__.py
│   │   └── euler.py             # Explicit Euler numerical integrator
│   │
│   ├── simulation/
│   │   ├── __init__.py
│   │   └── simulator.py         # Trajectory runner and simulation orchestrator
│   │
│   └── models/
│       ├── __init__.py
│       └── state.py             # State2D and StateDerivative2D representations
│
├── tests/
│   ├── __init__.py
│   ├── test_gravity.py          # Singularity, force directions, inverse-square invariants
│   ├── test_euler.py            # Decoupled step verification, analytical ODE steps
│   ├── test_initial_conditions.py # Circular orbital velocity and perpendicularity
│   └── test_simulator.py        # Simulation execution, trajectory length, validations
│
├── experiments/
│   └── README.md                # Reproducible numerical experiments guide
│
├── results/
│   └── .gitkeep                 # Output artifacts directory
│
├── docs/
│   ├── mathematical_model.md    # Formal mathematical specifications (LaTeX)
│   ├── numerical_methods.md     # Integrator formulation and truncation errors
│   └── experiments.md           # Error metrics and validation protocols
│
├── README.md                    # Project documentation
├── .gitignore                   # Git exclusion rules
├── requirements.txt             # Testing dependencies (pytest)
└── pyproject.toml               # Package configuration and pytest settings
```

---

## Development

### 1. Create Virtual Environment

Using Python 3 (3.9+ recommended):

```bash
python3 -m venv .venv
source .venv/bin/activate
```

### 2. Install Dependencies

OrbitSim core math relies exclusively on Python standard library (`math`). Install development/testing tools:

```bash
pip install -r requirements.txt
```

### 3. Run Tests

Execute the automated test suite via `pytest`:

```bash
pytest
```

Or run with verbose output:

```bash
pytest -v
```

---

## Current Status

**M7 — Fourth-Order Symplectic Integration**
- **M0**: 2D two-body equations, physical constants, and singularity handling established.
- **M1**: Explicit Euler convergence study ($p = 1.0$) across timestep refinements.
- **M2**: Time evolution of numerical error and physical invariants over 1 orbital period.
- **M3**: Long-term stability analysis of Explicit Euler over 100 orbital periods ($100T$).
- **M4**: Error decomposition into radial drift and along-track phase error; proof of phase lag dominance.
- **M5**: Classical RK4 4th-order integrator implementation, convergence ($p = 4.0$), and secular dissipation benchmark.
- **M6**: Symplectic Velocity Verlet integrator, machine-precision angular momentum preservation, bounded energy oscillations, and long-term stability validation.
- **M7**: Fourth-order symplectic integrator (Yoshida 1990) via symmetric composition; confirmed $p = 4.0$ convergence, bounded $\mathcal{O}(\Delta t^4)$ energy oscillations, machine-precision angular momentum, and cost-normalized 100T benchmark across four solvers.

---

## Roadmap

- Two-body model
- Euler integration
- Numerical error analysis
- RK4
- Solver comparison
- Long-term stability analysis
- Validation
- 3D orbital mechanics
- Perturbations
- Performance optimization
- Visualization
- Web interface

