# OrbitSim Experiments

This directory contains standalone, reproducible numerical experiments for analyzing orbital simulation dynamics, numerical integrator precision, and numerical errors.

## Experiment Principles

1. **Reproducibility**: Every experiment must have deterministic initial conditions and fixed parameters (timestep, duration, integrator).
2. **Decoupled Scripts**: Experiments import the core engine from `src/` and output metrics/logs into `results/`.
3. **Hypothesis-Driven**: Each experiment investigates a specific question (e.g., energy conservation drift as a function of timestep $\Delta t$, phase errors over multiple orbits).

## Planned Experiments (M1+)

- `exp01_euler_step_convergence`: Truncation error convergence rate analysis for Explicit Euler.
- `exp02_euler_energy_drift`: Long-term orbital energy and eccentricity growth in Euler integration.
- `exp03_euler_vs_analytical`: Trajectory deviation from the analytical Keplerian orbit.

