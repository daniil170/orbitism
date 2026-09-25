"""Experiment M7: Fourth-Order Symplectic Integration.

Investigates whether the 4th-order symplectic Yoshida integrator achieves simultaneous
high discretization accuracy and long-term geometric fidelity for conservative two-body motion.

Evaluates:
1. 1-orbit convergence rates across dt in [1.875, 60.0] seconds.
2. 100-orbit long-term behavior comparing:
   - Explicit Euler (order 1, non-symplectic)
   - Classical RK4 (order 4, non-symplectic)
   - Velocity Verlet (order 2, symplectic)
   - Yoshida 4 (order 4, symplectic)
3. Cost-normalized comparison (matched force evaluations vs fixed timestep).
"""

import csv
import math
from pathlib import Path
import sys
import time as pytime
from typing import Any, Dict, List, Optional, Tuple, Type, Union

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.integrators.euler import EulerIntegrator
from src.integrators.rk4 import RK4Integrator
from src.integrators.velocity_verlet import VelocityVerletIntegrator
from src.integrators.yoshida4 import Yoshida4Integrator
from src.models.state import State2D
from src.physics.gravity import (
    MU_EARTH,
    R_EARTH,
    create_circular_orbit_state,
)
from src.physics.analytical import (
    circular_orbit_period,
    circular_orbit_exact_state,
)
from src.simulation.simulator import Simulator
from src.simulation.metrics import (
    DecomposedErrorRecord,
    compute_decomposed_error_history,
    convergence_order,
    dominance_ratio,
    position_error,
    velocity_error,
)

DEFAULT_DT_VALUES: List[float] = [60.0, 30.0, 15.0, 7.5, 3.75, 1.875]
DEFAULT_ALTITUDE_KM: float = 400.0
DEFAULT_NUM_ORBITS: float = 100.0
DEFAULT_CHECKPOINTS_T: List[float] = [1.0, 5.0, 10.0, 20.0, 50.0, 100.0]

FORCE_EVALS_PER_STEP: Dict[str, int] = {
    "Euler": 1,
    "Velocity Verlet": 2,
    "RK4": 4,
    "Yoshida 4": 6,
}

SUMMARY_FIELDNAMES: List[str] = [
    "method",
    "dt",
    "orbit_count",
    "time",
    "e_r",
    "e_v",
    "e_R",
    "e_T",
    "Δr",
    "Δθ",
    "ρ_R",
    "ρ_T",
    "χ",
    "δε",
    "δh",
    "force_evals",
]

LONG_TERM_COMPARISON_FIELDNAMES: List[str] = [
    "method",
    "dt",
    "evals_per_step",
    "total_force_evals",
    "max_rel_energy_deviation",
    "final_rel_energy_deviation",
    "max_rel_angular_momentum_deviation",
    "max_radial_deviation",
    "final_radial_deviation",
    "final_phase_error",
    "final_state_error",
    "final_velocity_error",
]

COST_NORMALIZED_FIELDNAMES: List[str] = [
    "budget_name",
    "target_evals",
    "method",
    "dt",
    "evals_per_step",
    "total_steps",
    "actual_force_evals",
    "max_rel_energy_deviation",
    "final_rel_energy_deviation",
    "max_rel_angular_momentum_deviation",
    "max_radial_deviation",
    "final_radial_deviation",
    "final_phase_error",
    "final_state_error",
]

CONVERGENCE_FIELDNAMES: List[str] = [
    "method",
    "dt",
    "num_steps",
    "force_evals",
    "e_r",
    "e_v",
    "order_p_r",
    "ratio_r",
    "max_de",
    "max_dh",
    "final_dr",
    "final_dtheta",
]


def run_single_simulation(
    integrator_cls: Type[Any],
    dt: float,
    num_orbits: float = DEFAULT_NUM_ORBITS,
    altitude_km: float = DEFAULT_ALTITUDE_KM,
    mu: float = MU_EARTH,
    r_earth: float = R_EARTH,
) -> Tuple[List[DecomposedErrorRecord], float, int]:
    """Execute a single orbital simulation and return error history, orbital period, and total steps."""
    r0 = r_earth + altitude_km
    t_orbit = circular_orbit_period(r0, mu=mu)
    t_final = num_orbits * t_orbit
    initial_state = create_circular_orbit_state(altitude_km=altitude_km, mu=mu, r_earth=r_earth)

    sim = Simulator(initial_state=initial_state, mu=mu, integrator=integrator_cls)
    trajectory = sim.run(dt=dt, t_final=t_final)
    times = sim.times

    history = compute_decomposed_error_history(times=times, trajectory=trajectory, r0=r0, mu=mu)
    total_steps = len(trajectory) - 1
    return history, t_orbit, total_steps


def run_convergence_study(
    altitude_km: float = DEFAULT_ALTITUDE_KM,
    dt_values: Optional[List[float]] = None,
    mu: float = MU_EARTH,
    r_earth: float = R_EARTH,
) -> List[Dict[str, Any]]:
    """Run 1-orbit convergence study for Yoshida 4 across decreasing timesteps."""
    if dt_values is None:
        dt_values = DEFAULT_DT_VALUES

    r0 = r_earth + altitude_km
    t_orbit = circular_orbit_period(r0, mu=mu)
    exact_state = circular_orbit_exact_state(t_orbit, r0=r0, mu=mu)
    initial_state = create_circular_orbit_state(altitude_km=altitude_km, mu=mu, r_earth=r_earth)

    records = []
    prev_dt: Optional[float] = None
    prev_er: Optional[float] = None

    for dt in dt_values:
        sim = Simulator(initial_state=initial_state, mu=mu, integrator=Yoshida4Integrator)
        trajectory = sim.run(dt=dt, t_final=t_orbit)
        history = compute_decomposed_error_history(times=sim.times, trajectory=trajectory, r0=r0, mu=mu)

        final_sim_state = trajectory[-1]
        er = position_error(final_sim_state, exact_state)
        ev = velocity_error(final_sim_state, exact_state)

        num_steps = len(trajectory) - 1
        force_evals = num_steps * FORCE_EVALS_PER_STEP["Yoshida 4"]

        if prev_dt is not None and prev_er is not None and er > 0.0:
            order_p = convergence_order(prev_er, er, prev_dt, dt)
            ratio = prev_er / er
        else:
            order_p = 0.0
            ratio = 1.0

        max_de = max(abs(h.energy_error) for h in history)
        max_dh = max(abs(h.angular_momentum_error) for h in history)
        final_dr = abs(history[-1].delta_r)
        final_dtheta = abs(history[-1].delta_theta)

        records.append({
            "method": "Yoshida 4",
            "dt": dt,
            "num_steps": num_steps,
            "force_evals": force_evals,
            "e_r": er,
            "e_v": ev,
            "order_p_r": order_p,
            "ratio_r": ratio,
            "max_de": max_de,
            "max_dh": max_dh,
            "final_dr": final_dr,
            "final_dtheta": final_dtheta,
        })

        prev_dt = dt
        prev_er = er

    return records


def run_long_term_grid_experiment(
    dt_values: Optional[List[float]] = None,
    num_orbits: float = DEFAULT_NUM_ORBITS,
    altitude_km: float = DEFAULT_ALTITUDE_KM,
    mu: float = MU_EARTH,
    r_earth: float = R_EARTH,
) -> Tuple[Dict[str, Dict[float, List[DecomposedErrorRecord]]], float]:
    """Run 100-orbit simulations for Euler, RK4, Velocity Verlet, and Yoshida 4 across standard dt values."""
    if dt_values is None:
        dt_values = DEFAULT_DT_VALUES

    methods: List[Tuple[str, Type[Any]]] = [
        ("Euler", EulerIntegrator),
        ("RK4", RK4Integrator),
        ("Velocity Verlet", VelocityVerletIntegrator),
        ("Yoshida 4", Yoshida4Integrator),
    ]

    all_histories: Dict[str, Dict[float, List[DecomposedErrorRecord]]] = {}
    t_orbit = 0.0

    for method_name, integrator_cls in methods:
        start_t = pytime.perf_counter()
        all_histories[method_name] = {}
        for dt in dt_values:
            history, t_orbit, _ = run_single_simulation(
                integrator_cls=integrator_cls,
                dt=dt,
                num_orbits=num_orbits,
                altitude_km=altitude_km,
                mu=mu,
                r_earth=r_earth,
            )
            all_histories[method_name][dt] = history
        elapsed = pytime.perf_counter() - start_t
        print(f"[{method_name}] Completed 100T across {len(dt_values)} timesteps in {elapsed:.2f} s")

    return all_histories, t_orbit


def extract_checkpoint_summaries(
    all_histories: Dict[str, Dict[float, List[DecomposedErrorRecord]]],
    orbital_period: float,
    checkpoints_T: Optional[List[float]] = None,
) -> List[Dict[str, Any]]:
    """Extract metrics at discrete checkpoints (1T, 5T, 10T, 20T, 50T, 100T)."""
    if checkpoints_T is None:
        checkpoints_T = DEFAULT_CHECKPOINTS_T

    summary_records: List[Dict[str, Any]] = []

    for method_name, histories_by_dt in all_histories.items():
        evals_per_step = FORCE_EVALS_PER_STEP[method_name]
        for dt, history in histories_by_dt.items():
            for checkpoint_T in checkpoints_T:
                target_t = checkpoint_T * orbital_period
                # Find closest record index
                closest_idx = min(range(len(history)), key=lambda i: abs(history[i].t - target_t))
                record = history[closest_idx]
                chi = dominance_ratio(record.e_R, record.e_T)
                force_evals = closest_idx * evals_per_step

                summary_records.append({
                    "method": method_name,
                    "dt": dt,
                    "orbit_count": checkpoint_T,
                    "time": record.t,
                    "e_r": record.position_error,
                    "e_v": record.velocity_error,
                    "e_R": record.e_R,
                    "e_T": record.e_T,
                    "Δr": record.delta_r,
                    "Δθ": record.delta_theta,
                    "ρ_R": record.fraction_R,
                    "ρ_T": record.fraction_T,
                    "χ": chi,
                    "δε": record.energy_error,
                    "δh": record.angular_momentum_error,
                    "force_evals": force_evals,
                })

    return summary_records


def compute_long_term_comparison_records(
    all_histories: Dict[str, Dict[float, List[DecomposedErrorRecord]]],
) -> List[Dict[str, Any]]:
    """Compute overall 100T trajectory-wide max/final statistics across all methods and dt values."""
    records: List[Dict[str, Any]] = []

    for method_name, histories_by_dt in all_histories.items():
        evals_per_step = FORCE_EVALS_PER_STEP[method_name]
        for dt, history in histories_by_dt.items():
            total_steps = len(history) - 1
            total_force_evals = total_steps * evals_per_step

            max_de = max(abs(h.energy_error) for h in history)
            final_de = history[-1].energy_error
            max_dh = max(abs(h.angular_momentum_error) for h in history)
            max_dr = max(abs(h.delta_r) for h in history)
            final_dr = abs(history[-1].delta_r)
            final_dtheta = history[-1].delta_theta
            final_er = history[-1].position_error
            final_ev = history[-1].velocity_error

            records.append({
                "method": method_name,
                "dt": dt,
                "evals_per_step": evals_per_step,
                "total_force_evals": total_force_evals,
                "max_rel_energy_deviation": max_de,
                "final_rel_energy_deviation": final_de,
                "max_rel_angular_momentum_deviation": max_dh,
                "max_radial_deviation": max_dr,
                "final_radial_deviation": final_dr,
                "final_phase_error": final_dtheta,
                "final_state_error": final_er,
                "final_velocity_error": final_ev,
            })

    return records


def run_cost_normalized_comparison(
    altitude_km: float = DEFAULT_ALTITUDE_KM,
    num_orbits: float = DEFAULT_NUM_ORBITS,
    mu: float = MU_EARTH,
    r_earth: float = R_EARTH,
) -> List[Dict[str, Any]]:
    """Run equal-cost comparison by matching the number of force evaluations across methods.

    Evaluations per step:
        Euler: 1
        Velocity Verlet: 2
        RK4: 4
        Yoshida 4: 6

    Matched budgets over 100 orbits (T_final ~ 555255 s):
        Budget 1 (~55,525 evals):  Y4 (60s), RK4 (40s),  VV (20s),  Euler (10s)
        Budget 2 (~111,050 evals): Y4 (30s), RK4 (20s),  VV (10s),  Euler (5s)
        Budget 3 (~222,100 evals): Y4 (15s), RK4 (10s),  VV (5s),   Euler (2.5s)
        Budget 4 (~444,200 evals): Y4 (7.5s), RK4 (5s),  VV (2.5s), Euler (1.25s)
    """
    budgets = [
        {
            "name": "Budget_1_55k_evals",
            "target_evals": 55525,
            "methods": [
                ("Yoshida 4", Yoshida4Integrator, 60.0),
                ("RK4", RK4Integrator, 40.0),
                ("Velocity Verlet", VelocityVerletIntegrator, 20.0),
                ("Euler", EulerIntegrator, 10.0),
            ],
        },
        {
            "name": "Budget_2_111k_evals",
            "target_evals": 111050,
            "methods": [
                ("Yoshida 4", Yoshida4Integrator, 30.0),
                ("RK4", RK4Integrator, 20.0),
                ("Velocity Verlet", VelocityVerletIntegrator, 10.0),
                ("Euler", EulerIntegrator, 5.0),
            ],
        },
        {
            "name": "Budget_3_222k_evals",
            "target_evals": 222100,
            "methods": [
                ("Yoshida 4", Yoshida4Integrator, 15.0),
                ("RK4", RK4Integrator, 10.0),
                ("Velocity Verlet", VelocityVerletIntegrator, 5.0),
                ("Euler", EulerIntegrator, 2.5),
            ],
        },
        {
            "name": "Budget_4_444k_evals",
            "target_evals": 444200,
            "methods": [
                ("Yoshida 4", Yoshida4Integrator, 7.5),
                ("RK4", RK4Integrator, 5.0),
                ("Velocity Verlet", VelocityVerletIntegrator, 2.5),
                ("Euler", EulerIntegrator, 1.25),
            ],
        },
    ]

    records: List[Dict[str, Any]] = []

    for b in budgets:
        b_name = b["name"]
        target_evals = b["target_evals"]
        for method_name, integrator_cls, dt in b["methods"]:
            history, _, total_steps = run_single_simulation(
                integrator_cls=integrator_cls,
                dt=dt,
                num_orbits=num_orbits,
                altitude_km=altitude_km,
                mu=mu,
                r_earth=r_earth,
            )
            evals_per_step = FORCE_EVALS_PER_STEP[method_name]
            actual_evals = total_steps * evals_per_step

            max_de = max(abs(h.energy_error) for h in history)
            final_de = history[-1].energy_error
            max_dh = max(abs(h.angular_momentum_error) for h in history)
            max_dr = max(abs(h.delta_r) for h in history)
            final_dr = abs(history[-1].delta_r)
            final_dtheta = history[-1].delta_theta
            final_er = history[-1].position_error

            records.append({
                "budget_name": b_name,
                "target_evals": target_evals,
                "method": method_name,
                "dt": dt,
                "evals_per_step": evals_per_step,
                "total_steps": total_steps,
                "actual_force_evals": actual_evals,
                "max_rel_energy_deviation": max_de,
                "final_rel_energy_deviation": final_de,
                "max_rel_angular_momentum_deviation": max_dh,
                "max_radial_deviation": max_dr,
                "final_radial_deviation": final_dr,
                "final_phase_error": final_dtheta,
                "final_state_error": final_er,
            })

    return records


def save_csv(records: List[Dict[str, Any]], fieldnames: List[str], filepath: Path) -> None:
    """Save records to CSV with UTF-8 encoding."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for r in records:
            writer.writerow(r)


def save_trajectories_csv(
    all_histories: Dict[str, Dict[float, List[DecomposedErrorRecord]]],
    orbital_period: float,
    filepath: Path,
) -> None:
    """Save full step trajectory dataset for all methods and timesteps to CSV."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_FIELDNAMES, lineterminator="\n")
        writer.writeheader()

        for method_name, histories_by_dt in all_histories.items():
            evals_per_step = FORCE_EVALS_PER_STEP[method_name]
            for dt, history in histories_by_dt.items():
                for idx, record in enumerate(history):
                    chi = dominance_ratio(record.e_R, record.e_T)
                    writer.writerow({
                        "method": method_name,
                        "dt": dt,
                        "orbit_count": record.t / orbital_period,
                        "time": record.t,
                        "e_r": record.position_error,
                        "e_v": record.velocity_error,
                        "e_R": record.e_R,
                        "e_T": record.e_T,
                        "Δr": record.delta_r,
                        "Δθ": record.delta_theta,
                        "ρ_R": record.fraction_R,
                        "ρ_T": record.fraction_T,
                        "χ": chi,
                        "δε": record.energy_error,
                        "δh": record.angular_momentum_error,
                        "force_evals": idx * evals_per_step,
                    })


def main() -> None:
    """Run full Milestone M7 experimental pipeline."""
    print("=" * 70)
    print("OrbitSim — Milestone M7: Fourth-Order Symplectic Integration")
    print("=" * 70)

    results_dir = PROJECT_ROOT / "results"
    conv_path = results_dir / "exp07_convergence.csv"
    summary_path = results_dir / "exp07_fourth_order_symplectic_summary.csv"
    long_term_path = results_dir / "exp07_long_term_comparison.csv"
    cost_norm_path = results_dir / "exp07_cost_normalized_comparison.csv"
    trajectories_path = results_dir / "exp07_fourth_order_symplectic_integrators.csv"

    print("\nPhase 1: Running Yoshida 4 Convergence Study (1 Orbit)...")
    conv_records = run_convergence_study()
    save_csv(conv_records, CONVERGENCE_FIELDNAMES, conv_path)
    print(f"Convergence results saved to {conv_path}")
    for r in conv_records:
        print(f"  dt={r['dt']:5.3f}s: e_r={r['e_r']:11.4e} km, ratio={r['ratio_r']:5.2f}, order_p={r['order_p_r']:4.2f}")

    print("\nPhase 2: Running 100-Orbit Grid (Euler, RK4, Velocity Verlet, Yoshida 4)...")
    all_histories, t_orbit = run_long_term_grid_experiment()

    print("\nPhase 3: Extracting Checkpoint Summaries (1T, 5T, 10T, 20T, 50T, 100T)...")
    summary_records = extract_checkpoint_summaries(all_histories, orbital_period=t_orbit)
    save_csv(summary_records, SUMMARY_FIELDNAMES, summary_path)
    print(f"Checkpoint summary saved to {summary_path} ({len(summary_records)} rows)")

    print("\nPhase 4: Computing Full Trajectory Long-Term Statistics...")
    long_term_records = compute_long_term_comparison_records(all_histories)
    save_csv(long_term_records, LONG_TERM_COMPARISON_FIELDNAMES, long_term_path)
    print(f"Long-term comparison saved to {long_term_path} ({len(long_term_records)} rows)")

    print("\nPhase 5: Running Cost-Normalized Comparison (Matched Force Evaluation Budgets)...")
    cost_records = run_cost_normalized_comparison()
    save_csv(cost_records, COST_NORMALIZED_FIELDNAMES, cost_norm_path)
    print(f"Cost-normalized comparison saved to {cost_norm_path} ({len(cost_records)} rows)")

    print(f"\nPhase 6: Saving Full Trajectory Dataset to {trajectories_path}...")
    save_trajectories_csv(all_histories, orbital_period=t_orbit, filepath=trajectories_path)
    print("Full trajectories dataset saved.")

    print("\n" + "=" * 70)
    print("M7 Experiment Execution Completed Successfully.")
    print("=" * 70)


if __name__ == "__main__":
    main()
