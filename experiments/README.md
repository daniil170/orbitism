# OrbitSim Experiments

This directory contains standalone, reproducible numerical experiments for analyzing orbital simulation dynamics, numerical integrator precision, and numerical errors.

## Experiment Principles

1. **Reproducibility**: Every experiment must have deterministic initial conditions and fixed parameters (timestep, duration, integrator).
2. **Decoupled Scripts**: Experiments import the core engine from `src/` and output summary metrics into `results/`.
3. **Hypothesis-Driven**: Each experiment investigates a specific scientific question (e.g., convergence rate, energy drift, phase deviation).

## Implemented Experiments (M1)

- **`exp01_euler_convergence.py`**: Explicit Euler convergence order study for a circular Low Earth Orbit (LEO, altitude 400 km) over one complete analytical period ($T \approx 5553.62$ s). Evaluates position and velocity errors across $\Delta t \in \{60, 30, 15, 7.5, 3.75, 1.875\}$ s, empirical convergence orders, and energy/angular momentum drift. Output is written to `results/exp01_euler_convergence.csv`.

## Planned Experiments (M2+)

- `exp02_euler_energy_drift`: Long-term orbital energy and eccentricity growth in Euler integration over 10 to 100 orbits.
- `exp03_euler_vs_analytical`: Analytical phase vs radius error decomposition.
