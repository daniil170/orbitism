"""Experiment M8: Systematic Accuracy, Cost, and Invariant Comparison.

Milestone M8.1: Analytical Reference + Convergence Experiment across 4 solvers.
Milestone M8.2: Long-Term & Cost-Normalized Comparison (Architecture Skeleton - Deferred).

Research Question (M8.1):
    "How does the numerical error of Explicit Euler, Classical RK4, Velocity Verlet,
    and Yoshida 4 scale with timestep dt when integrating a known circular two-body orbit
    against the exact closed-form analytical Keplerian reference solution?"

Solvers:
    1. Explicit Euler (order 1, non-symplectic)
    2. Velocity Verlet (order 2, symplectic)
    3. Classical RK4 (order 4, non-symplectic)
    4. Yoshida 4 (order 4, symplectic)
"""

import csv
from dataclasses import dataclass
import math
import os
from pathlib import Path
import sys
from typing import Any, Callable, Dict, List, Optional, Tuple, Type

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.integrators.euler import EulerIntegrator
from src.integrators.rk4 import RK4Integrator
from src.integrators.velocity_verlet import VelocityVerletIntegrator
from src.integrators.yoshida4 import Yoshida4Integrator
from src.models.state import State2D
from src.physics.gravity import MU_EARTH, R_EARTH, create_circular_orbit_state
from src.physics.analytical import (
    circular_orbit_exact_state,
    circular_orbit_period,
    specific_orbital_energy,
    specific_angular_momentum,
)
from src.simulation.simulator import Simulator
from src.simulation.metrics import (
    compute_decomposed_error_history,
    convergence_order,
    position_error,
    velocity_error,
)

# Standard experimental grid
DEFAULT_ALTITUDE_KM: float = 400.0
DEFAULT_NUM_ORBITS: float = 100.0
DEFAULT_DT_VALUES: List[float] = [60.0, 30.0, 15.0, 7.5, 3.75, 1.875]

# Force evaluations per outer timestep in CURRENT decoupled implementations:
# Note: This is an implementation-level force-evaluation cost proxy for the
# current codebase, not a universal model of wall-clock computational cost.
# These values reflect the current modular architecture (un-cached)
# and do not represent theoretical lower bounds for all possible implementations.
CURRENT_IMPLEMENTATION_EVALS_PER_STEP: Dict[str, int] = {
    "Euler": 1,
    "Velocity Verlet": 2,
    "RK4": 4,
    "Yoshida 4": 6,
}

# Solvers compared in M8
SOLVER_REGISTRY: Dict[str, Type[Any]] = {
    "Euler": EulerIntegrator,
    "Velocity Verlet": VelocityVerletIntegrator,
    "RK4": RK4Integrator,
    "Yoshida 4": Yoshida4Integrator,
}

# Standard output schemas for M8
PLANNED_CONVERGENCE_COLUMNS: List[str] = [
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
    "max_e_r",
    "max_e_v",
    "order_p_v",
    "final_de",
    "final_dh",
    "max_dr",
]

DEFAULT_LONG_TERM_DT_VALUES: List[float] = [60.0, 30.0, 15.0]
DEFAULT_CHECKPOINTS_T: List[float] = [1.0, 5.0, 10.0, 20.0, 50.0, 100.0]

PLANNED_LONG_TERM_FIXED_DT_COLUMNS: List[str] = [
    "method",
    "dt",
    "evals_per_step",
    "total_steps",
    "actual_force_evals",
    "max_position_error",
    "final_position_error",
    "max_velocity_error",
    "final_velocity_error",
    "max_radial_deviation",
    "final_radial_deviation",
    "final_phase_error",
    "max_relative_energy_error",
    "final_relative_energy_error",
    "max_relative_angular_momentum_error",
    "final_relative_angular_momentum_error",
]

PLANNED_LONG_TERM_CHECKPOINTS_COLUMNS: List[str] = [
    "method",
    "dt",
    "orbit_count",
    "target_time",
    "actual_time",
    "position_error",
    "velocity_error",
    "radial_deviation",
    "phase_error",
    "relative_energy_error",
    "relative_angular_momentum_error",
    "force_evals",
]


PLANNED_ACCURACY_VS_COST_COLUMNS: List[str] = [
    "budget_name",
    "target_evals",
    "method",
    "dt",
    "evals_per_step",
    "total_steps",
    "actual_force_evals",
    "max_position_error",
    "final_position_error",
    "max_radial_deviation",
    "final_radial_deviation",
    "final_phase_error",
    "max_relative_energy_error",
    "final_relative_energy_error",
    "max_relative_angular_momentum_error",
]

PLANNED_INVARIANT_QUALITY_COLUMNS: List[str] = [
    "method",
    "dt",
    "max_de_trajectory",
    "final_de_checkpoint",
    "de_character",  # 'secular_drift' vs 'bounded_oscillation'
    "max_dh_trajectory",
    "final_dh_checkpoint",
    "dh_character",  # 'roundoff_limited' vs 'secular_drift'
]


def compute_force_evals_per_step(method: str) -> int:
    """Return number of force evaluations per external integration step."""
    if method not in CURRENT_IMPLEMENTATION_EVALS_PER_STEP:
        raise KeyError(f"Unknown solver method: {method}")
    return CURRENT_IMPLEMENTATION_EVALS_PER_STEP[method]


def compute_total_force_evals(method: str, num_steps: int) -> int:
    """Compute total force evaluations given solver method and number of integration steps."""
    return num_steps * compute_force_evals_per_step(method)


def compute_dt_for_budget(method: str, t_final: float, target_evals: int) -> float:
    """Compute required timestep dt to achieve target force evaluation budget over t_final."""
    cost_per_step = compute_force_evals_per_step(method)
    target_steps = target_evals / cost_per_step
    return t_final / target_steps


@dataclass(frozen=True)
class CostBudgetDefinition:
    """Specification of an equal-cost evaluation budget across solvers."""

    budget_name: str
    target_evals: int
    dt_assignments: Dict[str, float]


# Planned matched budgets for 100 orbits (T_final ~ 555255 s)
PLANNED_MATCHED_BUDGETS: List[CostBudgetDefinition] = [
    CostBudgetDefinition(
        budget_name="Budget_55k",
        target_evals=55525,
        dt_assignments={"Euler": 10.0, "Velocity Verlet": 20.0, "RK4": 40.0, "Yoshida 4": 60.0},
    ),
    CostBudgetDefinition(
        budget_name="Budget_111k",
        target_evals=111050,
        dt_assignments={"Euler": 5.0, "Velocity Verlet": 10.0, "RK4": 20.0, "Yoshida 4": 30.0},
    ),
    CostBudgetDefinition(
        budget_name="Budget_222k",
        target_evals=222100,
        dt_assignments={"Euler": 2.5, "Velocity Verlet": 5.0, "RK4": 10.0, "Yoshida 4": 15.0},
    ),
    CostBudgetDefinition(
        budget_name="Budget_444k",
        target_evals=444200,
        dt_assignments={"Euler": 1.25, "Velocity Verlet": 2.5, "RK4": 5.0, "Yoshida 4": 7.5},
    ),
]


# ==============================================================================
# M8.1: ANALYTICAL REFERENCE + CONVERGENCE STUDY
# ==============================================================================

def run_convergence_study(
    dt_values: Optional[List[float]] = None,
    altitude_km: float = DEFAULT_ALTITUDE_KM,
    mu: float = MU_EARTH,
    r_earth: float = R_EARTH,
) -> List[Dict[str, Any]]:
    """Execute M8.1 analytical convergence study across 4 solvers over 1 orbital period.

    Compares Explicit Euler, Velocity Verlet, RK4, and Yoshida 4 against the exact
    closed-form analytical Keplerian circular orbit solution.

    Args:
        dt_values: List of timestep sizes in seconds (defaults to DEFAULT_DT_VALUES).
        altitude_km: Orbit altitude in km (defaults to 400 km).
        mu: Gravitational parameter in km^3/s^2.
        r_earth: Central body radius in km.

    Returns:
        List of 24 metric records (4 solvers x 6 timesteps).
    """
    if dt_values is None:
        dt_values = list(DEFAULT_DT_VALUES)

    r0 = r_earth + altitude_km
    t_final = circular_orbit_period(r0, mu=mu)
    initial_state = create_circular_orbit_state(
        altitude_km=altitude_km, mu=mu, r_earth=r_earth
    )

    records: List[Dict[str, Any]] = []

    for method_name, integrator_cls in SOLVER_REGISTRY.items():
        evals_per_step = CURRENT_IMPLEMENTATION_EVALS_PER_STEP[method_name]
        prev_dt: Optional[float] = None
        prev_er: Optional[float] = None
        prev_ev: Optional[float] = None

        for dt in dt_values:
            sim = Simulator(initial_state=initial_state, mu=mu, integrator=integrator_cls)
            trajectory = sim.run(dt=dt, t_final=t_final)
            times = sim.times

            num_steps = len(trajectory) - 1
            force_evals = num_steps * evals_per_step

            history = compute_decomposed_error_history(
                times=times, trajectory=trajectory, r0=r0, mu=mu
            )

            # Verification of final physical time
            final_record = history[-1]
            assert abs(final_record.t - t_final) < 1e-11, (
                f"Physical final time mismatch: {final_record.t} vs {t_final}"
            )

            # Errors against exact analytical solution
            final_er = final_record.position_error
            final_ev = final_record.velocity_error
            max_er = max(h.position_error for h in history)
            max_ev = max(h.velocity_error for h in history)

            # Empirical convergence orders
            if prev_dt is not None and prev_er is not None and final_er > 0.0:
                order_p_r = convergence_order(prev_er, final_er, prev_dt, dt)
                ratio_r = prev_er / final_er
            else:
                order_p_r = 0.0
                ratio_r = 1.0

            if prev_dt is not None and prev_ev is not None and final_ev > 0.0:
                order_p_v = convergence_order(prev_ev, final_ev, prev_dt, dt)
            else:
                order_p_v = 0.0

            # Invariant errors
            max_de = max(abs(h.energy_error) for h in history)
            final_de = final_record.energy_error
            max_dh = max(abs(h.angular_momentum_error) for h in history)
            final_dh = final_record.angular_momentum_error

            # Radial and along-track phase metrics
            max_dr = max(abs(h.delta_r) for h in history)
            final_dr = abs(final_record.delta_r)
            final_dtheta = abs(final_record.delta_theta)

            record = {
                "method": method_name,
                "dt": dt,
                "num_steps": num_steps,
                "force_evals": force_evals,
                "e_r": final_er,
                "e_v": final_ev,
                "order_p_r": order_p_r,
                "ratio_r": ratio_r,
                "max_de": max_de,
                "max_dh": max_dh,
                "final_dr": final_dr,
                "final_dtheta": final_dtheta,
                "max_e_r": max_er,
                "max_e_v": max_ev,
                "order_p_v": order_p_v,
                "final_de": final_de,
                "final_dh": final_dh,
                "max_dr": max_dr,
            }
            records.append(record)

            prev_dt = dt
            prev_er = final_er
            prev_ev = final_ev

    return records


def save_convergence_csv(records: List[Dict[str, Any]], filepath: Path) -> None:
    """Save convergence benchmark records to CSV."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=PLANNED_CONVERGENCE_COLUMNS, lineterminator="\n"
        )
        writer.writeheader()
        for r in records:
            writer.writerow(r)


def plot_convergence_results(
    records: List[Dict[str, Any]], output_dir: Path
) -> List[Path]:
    """Generate convergence analysis graphs for M8.1.

    Creates:
    1. exp08_convergence_position.png (log-log position error vs dt)
    2. exp08_convergence_velocity.png (log-log velocity error vs dt)
    3. exp08_convergence_invariants.png (energy and angular momentum errors vs dt)

    Args:
        records: List of convergence study records.
        output_dir: Directory where plots will be saved.

    Returns:
        List of paths to generated image files.
    """
    os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_ROOT / ".cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_dir.mkdir(parents=True, exist_ok=True)
    generated_plots: List[Path] = []

    by_method: Dict[str, List[Dict[str, Any]]] = {}
    for r in records:
        by_method.setdefault(r["method"], []).append(r)

    colors = {
        "Euler": "#d95f02",
        "Velocity Verlet": "#7570b3",
        "RK4": "#1b9e77",
        "Yoshida 4": "#e7298a",
    }
    markers = {
        "Euler": "s",
        "Velocity Verlet": "^",
        "RK4": "o",
        "Yoshida 4": "D",
    }

    # 1. Position Error vs dt
    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
    for method_name, rows in by_method.items():
        dts = [r["dt"] for r in rows]
        ers = [r["e_r"] for r in rows]
        ax.loglog(
            dts,
            ers,
            label=f"{method_name}",
            color=colors.get(method_name, "black"),
            marker=markers.get(method_name, "o"),
            linewidth=1.8,
            markersize=6,
        )

    # Reference slope lines
    dt_ref = [1.875, 60.0]
    ax.loglog(dt_ref, [50.0 * (dt / 1.875) ** 1 for dt in dt_ref], "k:", alpha=0.5, label=r"$\mathcal{O}(\Delta t)$ slope")
    ax.loglog(dt_ref, [0.05 * (dt / 1.875) ** 2 for dt in dt_ref], "k--", alpha=0.5, label=r"$\mathcal{O}(\Delta t^2)$ slope")
    ax.loglog(dt_ref, [5e-7 * (dt / 1.875) ** 4 for dt in dt_ref], "k-.", alpha=0.5, label=r"$\mathcal{O}(\Delta t^4)$ slope")

    ax.set_xlabel(r"Timestep $\Delta t$ [s]", fontsize=12)
    ax.set_ylabel(r"Final Position Error $e_r(T)$ [km]", fontsize=12)
    ax.set_title("M8.1: Discretization Error vs. Timestep (1 Orbit, Circular Reference)", fontsize=13)
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(loc="upper left", framealpha=0.9, fontsize=10)
    fig.tight_layout()

    pos_path = output_dir / "exp08_convergence_position.png"
    fig.savefig(pos_path)
    plt.close(fig)
    generated_plots.append(pos_path)

    # 2. Velocity Error vs dt
    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
    for method_name, rows in by_method.items():
        dts = [r["dt"] for r in rows]
        evs = [r["e_v"] for r in rows]
        ax.loglog(
            dts,
            evs,
            label=f"{method_name}",
            color=colors.get(method_name, "black"),
            marker=markers.get(method_name, "o"),
            linewidth=1.8,
            markersize=6,
        )

    ax.set_xlabel(r"Timestep $\Delta t$ [s]", fontsize=12)
    ax.set_ylabel(r"Final Velocity Error $e_v(T)$ [km/s]", fontsize=12)
    ax.set_title("M8.1: Velocity Error vs. Timestep (1 Orbit, Circular Reference)", fontsize=13)
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(loc="upper left", framealpha=0.9, fontsize=10)
    fig.tight_layout()

    vel_path = output_dir / "exp08_convergence_velocity.png"
    fig.savefig(vel_path)
    plt.close(fig)
    generated_plots.append(vel_path)

    # 3. Invariant Errors vs dt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5), dpi=150)
    for method_name, rows in by_method.items():
        dts = [r["dt"] for r in rows]
        des = [r["max_de"] for r in rows]
        dhs = [max(r["max_dh"], 1e-16) for r in rows]

        ax1.loglog(
            dts,
            des,
            label=method_name,
            color=colors.get(method_name, "black"),
            marker=markers.get(method_name, "o"),
            linewidth=1.8,
        )
        ax2.loglog(
            dts,
            dhs,
            label=method_name,
            color=colors.get(method_name, "black"),
            marker=markers.get(method_name, "o"),
            linewidth=1.8,
        )

    ax1.set_xlabel(r"Timestep $\Delta t$ [s]", fontsize=11)
    ax1.set_ylabel(r"$\max |\delta \varepsilon|$", fontsize=11)
    ax1.set_title("Relative Energy Error", fontsize=12)
    ax1.grid(True, which="both", ls="--", alpha=0.4)
    ax1.legend(loc="upper left", fontsize=9)

    ax2.set_xlabel(r"Timestep $\Delta t$ [s]", fontsize=11)
    ax2.set_ylabel(r"$\max |\delta h|$", fontsize=11)
    ax2.set_title("Relative Angular Momentum Error", fontsize=12)
    ax2.grid(True, which="both", ls="--", alpha=0.4)
    ax2.legend(loc="upper left", fontsize=9)

    fig.suptitle("M8.1: Invariant Conservation vs. Timestep (1 Orbit)", fontsize=13)
    fig.tight_layout()

    inv_path = output_dir / "exp08_convergence_invariants.png"
    fig.savefig(inv_path)
    plt.close(fig)
    generated_plots.append(inv_path)

    return generated_plots


def print_convergence_table(records: List[Dict[str, Any]], t_final: float) -> None:
    """Print ASCII formatted table of convergence results."""
    print("=" * 125)
    print("OrbitSim M8.1: Analytical Reference + Convergence Experiment")
    print(f"Physical Duration: T = {t_final:.6f} s (~{t_final / 60.0:.2f} min), Altitude = 400 km")
    print("=" * 125)
    header = (
        f"{'Method':<16} | {'dt (s)':>7} | {'Steps':>5} | {'Evals':>6} | "
        f"{'e_r (km)':>12} | {'e_v (km/s)':>12} | {'Order p_r':>9} | {'Order p_v':>9} | "
        f"{'Max |δε|':>11} | {'Max |δh|':>11} | {'Final Δθ (rad)':>14}"
    )
    print(header)
    print("-" * 125)
    for r in records:
        p_r_str = f"{r['order_p_r']:.2f}" if r["order_p_r"] > 0 else "N/A"
        p_v_str = f"{r['order_p_v']:.2f}" if r["order_p_v"] > 0 else "N/A"
        print(
            f"{r['method']:<16} | {r['dt']:>7.3f} | {r['num_steps']:>5d} | {r['force_evals']:>6d} | "
            f"{r['e_r']:>12.4e} | {r['e_v']:>12.4e} | {p_r_str:>9} | {p_v_str:>9} | "
            f"{r['max_de']:>11.4e} | {r['max_dh']:>11.4e} | {r['final_dtheta']:>14.4e}"
        )
    print("=" * 125)


# ==============================================================================
# M8.2a: LONG-TERM FIXED-DT COMPARISON (100 ORBITS)
# ==============================================================================

def run_long_term_fixed_dt_study(
    dt_values: Optional[List[float]] = None,
    num_orbits: float = DEFAULT_NUM_ORBITS,
    checkpoints_T: Optional[List[float]] = None,
    altitude_km: float = DEFAULT_ALTITUDE_KM,
    mu: float = MU_EARTH,
    r_earth: float = R_EARTH,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Execute M8.2a long-term fixed-dt comparison across 4 solvers over 100 orbits.

    Evaluates Explicit Euler, Velocity Verlet, RK4, and Yoshida 4 across nominal
    timesteps dt in [60.0, 30.0, 15.0] s.

    Args:
        dt_values: Nominal timesteps in seconds (defaults to [60.0, 30.0, 15.0]).
        num_orbits: Number of orbital periods (defaults to 100.0).
        checkpoints_T: Checkpoint horizons in orbital periods [1.0, 5.0, 10.0, 20.0, 50.0, 100.0].
        altitude_km: Orbit altitude in km (defaults to 400.0 km).
        mu: Gravitational parameter in km^3/s^2.
        r_earth: Central body radius in km.

    Returns:
        Tuple of (summary_records, checkpoint_records).
    """
    if dt_values is None:
        dt_values = list(DEFAULT_LONG_TERM_DT_VALUES)
    if checkpoints_T is None:
        checkpoints_T = list(DEFAULT_CHECKPOINTS_T)

    r0 = r_earth + altitude_km
    t_orbit = circular_orbit_period(r0, mu=mu)
    t_final = num_orbits * t_orbit
    initial_state = create_circular_orbit_state(
        altitude_km=altitude_km, mu=mu, r_earth=r_earth
    )

    summary_records: List[Dict[str, Any]] = []
    checkpoint_records: List[Dict[str, Any]] = []

    for method_name, integrator_cls in SOLVER_REGISTRY.items():
        evals_per_step = CURRENT_IMPLEMENTATION_EVALS_PER_STEP[method_name]
        for dt in dt_values:
            sim = Simulator(initial_state=initial_state, mu=mu, integrator=integrator_cls)
            trajectory = sim.run(dt=dt, t_final=t_final)
            times = sim.times

            total_steps = len(trajectory) - 1
            actual_force_evals = total_steps * evals_per_step

            history = compute_decomposed_error_history(
                times=times, trajectory=trajectory, r0=r0, mu=mu
            )

            # Verification of final physical time
            final_record = history[-1]
            assert abs(final_record.t - t_final) < 1e-11, (
                f"Physical final time mismatch: {final_record.t} vs {t_final}"
            )

            # Extract metrics at discrete checkpoint times
            for checkpoint_T in checkpoints_T:
                target_t = checkpoint_T * t_orbit
                closest_idx = min(range(len(history)), key=lambda i: abs(history[i].t - target_t))
                rec = history[closest_idx]
                checkpoint_evals = closest_idx * evals_per_step

                checkpoint_records.append({
                    "method": method_name,
                    "dt": dt,
                    "orbit_count": checkpoint_T,
                    "target_time": target_t,
                    "actual_time": rec.t,
                    "position_error": rec.position_error,
                    "velocity_error": rec.velocity_error,
                    "radial_deviation": rec.delta_r,
                    "phase_error": rec.delta_theta,
                    "relative_energy_error": rec.energy_error,
                    "relative_angular_momentum_error": rec.angular_momentum_error,
                    "force_evals": checkpoint_evals,
                })

            # Trajectory-wide maxima over discrete sampled states
            max_pos_err = max(h.position_error for h in history)
            max_vel_err = max(h.velocity_error for h in history)
            max_rad_dev = max(abs(h.delta_r) for h in history)
            max_rel_energy_err = max(abs(h.energy_error) for h in history)
            max_rel_ang_mom_err = max(abs(h.angular_momentum_error) for h in history)

            summary_records.append({
                "method": method_name,
                "dt": dt,
                "evals_per_step": evals_per_step,
                "total_steps": total_steps,
                "actual_force_evals": actual_force_evals,
                "max_position_error": max_pos_err,
                "final_position_error": final_record.position_error,
                "max_velocity_error": max_vel_err,
                "final_velocity_error": final_record.velocity_error,
                "max_radial_deviation": max_rad_dev,
                "final_radial_deviation": abs(final_record.delta_r),
                "final_phase_error": final_record.delta_theta,
                "max_relative_energy_error": max_rel_energy_err,
                "final_relative_energy_error": final_record.energy_error,
                "max_relative_angular_momentum_error": max_rel_ang_mom_err,
                "final_relative_angular_momentum_error": final_record.angular_momentum_error,
            })

    return summary_records, checkpoint_records


def save_long_term_csv(records: List[Dict[str, Any]], filepath: Path) -> None:
    """Save 100-orbit fixed-dt summary records to CSV."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=PLANNED_LONG_TERM_FIXED_DT_COLUMNS, lineterminator="\n"
        )
        writer.writeheader()
        for r in records:
            writer.writerow(r)


def save_checkpoints_csv(records: List[Dict[str, Any]], filepath: Path) -> None:
    """Save 100-orbit checkpoint records to CSV."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=PLANNED_LONG_TERM_CHECKPOINTS_COLUMNS, lineterminator="\n"
        )
        writer.writeheader()
        for r in records:
            writer.writerow(r)


def plot_long_term_results(
    checkpoint_records: List[Dict[str, Any]], output_dir: Path
) -> List[Path]:
    """Generate 100-orbit fixed-dt analysis plots for M8.2a.

    Creates:
    1. exp08_long_term_position_error.png (position error vs checkpoint orbit)
    2. exp08_long_term_energy_deviation.png (relative energy error vs checkpoint orbit)
    3. exp08_long_term_angular_momentum.png (relative angular momentum error vs checkpoint orbit)
    4. exp08_long_term_phase_error.png (phase error vs checkpoint orbit)

    Args:
        checkpoint_records: List of checkpoint records.
        output_dir: Destination directory.

    Returns:
        List of generated image paths.
    """
    os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_ROOT / ".cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_dir.mkdir(parents=True, exist_ok=True)
    generated_plots: List[Path] = []

    dt_values = [60.0, 30.0, 15.0]
    methods = ["Euler", "Velocity Verlet", "RK4", "Yoshida 4"]
    colors = {
        "Euler": "#d95f02",
        "Velocity Verlet": "#7570b3",
        "RK4": "#1b9e77",
        "Yoshida 4": "#e7298a",
    }
    markers = {
        "Euler": "s",
        "Velocity Verlet": "^",
        "RK4": "o",
        "Yoshida 4": "D",
    }

    by_method_dt: Dict[Tuple[str, float], List[Dict[str, Any]]] = {}
    for r in checkpoint_records:
        by_method_dt.setdefault((r["method"], r["dt"]), []).append(r)

    metrics_to_plot = [
        ("position_error", "exp08_long_term_position_error.png", "Position Error [km]", "M8.2a: Position Error vs. Orbit Count (100 Orbits)", True),
        ("relative_energy_error", "exp08_long_term_energy_deviation.png", "Relative Energy Error |δε|", "M8.2a: Relative Energy Deviation vs. Orbit Count", True),
        ("relative_angular_momentum_error", "exp08_long_term_angular_momentum.png", "Relative Angular Momentum Error |δh|", "M8.2a: Relative Angular Momentum Deviation vs. Orbit Count", True),
        ("phase_error", "exp08_long_term_phase_error.png", "Along-Track Phase Error |Δθ| [rad]", "M8.2a: Along-Track Phase Error vs. Orbit Count", True),
    ]

    for field_name, filename, ylabel, title, use_log in metrics_to_plot:
        fig, axes = plt.subplots(1, 3, figsize=(16, 5), dpi=150, sharey=True)
        for ax_idx, dt in enumerate(dt_values):
            ax = axes[ax_idx]
            for method in methods:
                series = by_method_dt.get((method, dt), [])
                if not series:
                    continue
                orbits = [s["orbit_count"] for s in series]
                raw_vals = [abs(s[field_name]) for s in series]
                vals = [max(v, 1e-16) for v in raw_vals]

                if use_log:
                    ax.semilogy(
                        orbits,
                        vals,
                        label=method,
                        color=colors[method],
                        marker=markers[method],
                        linewidth=1.8,
                        markersize=6,
                    )
                else:
                    ax.plot(
                        orbits,
                        vals,
                        label=method,
                        color=colors[method],
                        marker=markers[method],
                        linewidth=1.8,
                        markersize=6,
                    )

            ax.set_title(f"Nominal $\\Delta t = {dt}$ s", fontsize=12)
            ax.set_xlabel("Orbit Count", fontsize=11)
            ax.grid(True, which="both", ls="--", alpha=0.4)
            if ax_idx == 0:
                ax.set_ylabel(ylabel, fontsize=11)
                ax.legend(loc="best", fontsize=9, framealpha=0.9)

        fig.suptitle(title, fontsize=14)
        fig.tight_layout()

        out_path = output_dir / filename
        fig.savefig(out_path)
        plt.close(fig)
        generated_plots.append(out_path)

    return generated_plots


def print_long_term_table(records: List[Dict[str, Any]], t_final: float) -> None:
    """Print ASCII formatted table of 100-orbit fixed-dt summary results."""
    print("=" * 140)
    print("OrbitSim M8.2a: Long-Term Fixed-dt Comparison (100 Orbits)")
    print(f"Physical Duration: T_final = {t_final:.3f} s (~{t_final / 86400.0:.2f} days), Altitude = 400 km")
    print("=" * 140)
    header = (
        f"{'Method':<16} | {'dt (s)':>7} | {'Steps':>6} | {'Evals':>7} | "
        f"{'Max Pos Err (km)':>17} | {'Final Pos Err (km)':>18} | "
        f"{'Final Δθ (rad)':>15} | {'Max |δε|':>11} | {'Final δε':>11} | {'Max |δh|':>11}"
    )
    print(header)
    print("-" * 140)
    for r in records:
        print(
            f"{r['method']:<16} | {r['dt']:>7.3f} | {r['total_steps']:>6d} | {r['actual_force_evals']:>7d} | "
            f"{r['max_position_error']:>17.4e} | {r['final_position_error']:>18.4e} | "
            f"{r['final_phase_error']:>15.4e} | {r['max_relative_energy_error']:>11.4e} | "
            f"{r['final_relative_energy_error']:>11.4e} | {r['max_relative_angular_momentum_error']:>11.4e}"
        )
    print("=" * 140)


def run_planned_accuracy_vs_timestep() -> List[Dict[str, Any]]:
    """Execute Study 1 (M8.1): Error metrics vs. timestep across 4 solvers.

    Returns:
        List of convergence records for M8.1.
    """
    return run_convergence_study()


def run_planned_long_term_fixed_dt() -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Execute Study 2a (M8.2a): Long-term fixed-dt comparison across 4 solvers over 100 orbits.

    Returns:
        Tuple of (summary_records, checkpoint_records).
    """
    return run_long_term_fixed_dt_study()


def run_planned_accuracy_vs_cost() -> None:
    """Planned study 2b (M8.2b): Error metrics vs. force evaluations across 4 solvers.

    TODO (M8.2b): Execute matched budget runs and generate Pareto-optimal efficiency curves.
    """
    raise NotImplementedError(
        "M8.2b execution pending. Cost-normalized comparison is deferred."
    )


def run_planned_invariant_quality_comparison() -> None:
    """Planned study 3: Detailed classification of energy and angular momentum behavior.

    TODO: Classify bounded oscillations vs secular dissipation across solvers.
    """
    raise NotImplementedError(
        "Execution pending. Invariant quality classification is deferred."
    )


def main() -> None:
    """Execute Milestone M8.1 and M8.2a experimental pipelines."""
    print("=" * 80)
    print("OrbitSim — Milestone M8: Systematic Solver Benchmark")
    print("=" * 80)

    results_dir = PROJECT_ROOT / "results"

    # Phase 1: M8.1 Convergence Study
    print("\n>>> Phase 1: Running Milestone M8.1 (Analytical Reference + Convergence)...")
    r0 = R_EARTH + DEFAULT_ALTITUDE_KM
    t_orbit = circular_orbit_period(r0, mu=MU_EARTH)
    conv_records = run_convergence_study()
    print_convergence_table(conv_records, t_final=t_orbit)

    conv_csv = results_dir / "exp08_convergence.csv"
    save_convergence_csv(conv_records, conv_csv)
    print(f"Convergence results saved to: {conv_csv.resolve()}")

    print("Generating M8.1 convergence plots...")
    conv_plots = plot_convergence_results(conv_records, results_dir)
    for p in conv_plots:
        print(f"  Plot: {p.resolve()}")

    # Phase 2: M8.2a Long-Term Fixed-dt Study
    print("\n>>> Phase 2: Running Milestone M8.2a (Long-Term Fixed-dt, 100 Orbits)...")
    t_100 = 100.0 * t_orbit
    summary_records, checkpoint_records = run_long_term_fixed_dt_study()
    print_long_term_table(summary_records, t_final=t_100)

    long_term_csv = results_dir / "exp08_long_term_fixed_dt.csv"
    save_long_term_csv(summary_records, long_term_csv)
    print(f"Long-term summary saved to: {long_term_csv.resolve()} ({len(summary_records)} rows)")

    checkpoints_csv = results_dir / "exp08_long_term_checkpoints.csv"
    save_checkpoints_csv(checkpoint_records, checkpoints_csv)
    print(f"Checkpoints saved to: {checkpoints_csv.resolve()} ({len(checkpoint_records)} rows)")

    print("Generating M8.2a long-term plots...")
    lt_plots = plot_long_term_results(checkpoint_records, results_dir)
    for p in lt_plots:
        print(f"  Plot: {p.resolve()}")

    print("\n" + "=" * 80)
    print("Milestones M8.1 and M8.2a Execution Complete.")
    print("Milestone M8.2b (Cost-Normalized Comparison) is deferred and NOT executed.")
    print("=" * 80)


if __name__ == "__main__":
    main()
