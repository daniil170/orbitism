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

In the current development phase (Milestone M0), OrbitSim models **2D planar satellite motion** around Earth:
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

* **Explicit Euler (Forward Euler)**: First-order numerical integrator ($p = 1$) serving as the baseline solver.
  $$\mathbf{s}_{n+1} = \mathbf{s}_n + \mathbf{f}(t_n, \mathbf{s}_n) \cdot \Delta t$$
  The integrator is decoupled from the orbital physics domain and accepts arbitrary state derivative functions.

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

**M0 — Mathematical foundation / initial Euler prototype**
- 2D two-body equations and physical constants established.
- Singularity handling and invariant tests implemented.
- Decoupled Explicit Euler integrator operational.
- Initial simulation pipeline functional.
- Zero external runtime dependencies.

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

