"""Experiment M3: Long-Term Stability of Explicit Euler.

Investigates how numerical error and physical diagnostic invariants evolve over
multi-orbit time scales (0 <= t <= 100T) across multiple timesteps dt:
    error = f(t, dt) for t in [0, 100T]
"""

import csv
import math
from pathlib import Path
import sys
from typing import Dict, List, Optional, Tuple

# Ensure project root is in sys.path for direct script execution
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.models.state import State2D
from src.physics.gravity import (
    MU_EARTH,
    R_EARTH,
    create_circular_orbit_state,
)
from src.physics.analytical import circular_orbit_period
from src.simulation.simulator import Simulator
from src.simulation.metrics import (
    ErrorRecord,
    compute_error_history,
)

DEFAULT_DT_VALUES: List[float] = [60.0, 30.0, 15.0, 7.5, 3.75, 1.875]
DEFAULT_ALTITUDE_KM: float = 400.0
DEFAULT_NUM_ORBITS: float = 100.0
DEFAULT_CHECKPOINTS_T: List[float] = [1.0, 5.0, 10.0, 20.0, 50.0, 100.0]


def run_euler_long_term_experiment(
    dt_values: Optional[List[float]] = None,
    num_orbits: float = DEFAULT_NUM_ORBITS,
    altitude_km: float = DEFAULT_ALTITUDE_KM,
    mu: float = MU_EARTH,
    r_earth: float = R_EARTH,
) -> Tuple[Dict[float, List[ErrorRecord]], float]:
    """Run Explicit Euler long-term stability study over multi-orbit time scales (100T).

    For each timestep dt in dt_values, integrates the circular orbit from
    t = 0 to t = 100 * T (where T is analytical period) using step shortening at the end
    dt_step = min(dt, t_final - t) to prevent overshoot. Computes the complete error
    and invariant history at every discrete step.

    Args:
        dt_values: Timesteps in seconds to evaluate.
        num_orbits: Integration duration in units of orbital period T (default: 100.0).
        altitude_km: Orbital altitude in km.
        mu: Gravitational parameter in km^3/s^2.
        r_earth: Central body radius in km.

    Returns:
        Tuple of (experiment_histories, orbital_period_T), where experiment_histories
        maps timestep dt -> List[ErrorRecord].
    """
    if dt_values is None:
        dt_values = DEFAULT_DT_VALUES

    r0 = r_earth + altitude_km
    t_orbit = circular_orbit_period(r0, mu=mu)
    t_final = num_orbits * t_orbit

    initial_state = create_circular_orbit_state(
        altitude_km=altitude_km, mu=mu, r_earth=r_earth
    )

    experiment_histories: Dict[float, List[ErrorRecord]] = {}

    for dt in dt_values:
        sim = Simulator(initial_state=initial_state, mu=mu)
        trajectory = sim.run(dt=dt, t_final=t_final)
        times = sim.times

        history = compute_error_history(
            times=times,
            trajectory=trajectory,
            r0=r0,
            mu=mu,
        )
        experiment_histories[dt] = history

    return experiment_histories, t_orbit


def save_long_term_trajectory_to_csv(
    experiment_histories: Dict[float, List[ErrorRecord]],
    filepath: Path,
) -> None:
    """Save full time evolution histories for all timesteps to a single CSV file.

    Args:
        experiment_histories: Mapping of dt -> List[ErrorRecord].
        filepath: Destination CSV Path.
    """
    filepath.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "dt_s",
        "step",
        "t_s",
        "x_km",
        "y_km",
        "vx_km_s",
        "vy_km_s",
        "position_error_km",
        "velocity_error_km_s",
        "energy_error",
        "angular_momentum_error",
    ]

    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=fieldnames, lineterminator="\n"
        )
        writer.writeheader()

        for dt, history in experiment_histories.items():
            for step_idx, record in enumerate(history):
                writer.writerow(
                    {
                        "dt_s": dt,
                        "step": step_idx,
                        "t_s": record.t,
                        "x_km": record.state.x,
                        "y_km": record.state.y,
                        "vx_km_s": record.state.vx,
                        "vy_km_s": record.state.vy,
                        "position_error_km": record.position_error,
                        "velocity_error_km_s": record.velocity_error,
                        "energy_error": record.energy_error,
                        "angular_momentum_error": record.angular_momentum_error,
                    }
                )


def extract_checkpoint_summary(
    experiment_histories: Dict[float, List[ErrorRecord]],
    orbital_period: float,
    checkpoints_T: Optional[List[float]] = None,
) -> List[Dict[str, float]]:
    """Extract diagnostic metrics at discrete orbital period checkpoints (1T, 5T, 10T, 20T, 50T, 100T).

    Strategy:
        For each requested checkpoint N * T, selects the record in history with minimum
        absolute timestamp difference |t - N * T|.

    Args:
        experiment_histories: Mapping of dt -> List[ErrorRecord].
        orbital_period: Analytical orbital period T in seconds.
        checkpoints_T: List of checkpoint multiples of T (default: [1.0, 5.0, 10.0, 20.0, 50.0, 100.0]).

    Returns:
        List of dictionaries containing checkpoint metrics.
    """
    if checkpoints_T is None:
        checkpoints_T = DEFAULT_CHECKPOINTS_T

    summary_records: List[Dict[str, float]] = []

    for dt, history in experiment_histories.items():
        for checkpoint_T in checkpoints_T:
            target_t = checkpoint_T * orbital_period
            closest_record = min(history, key=lambda r: abs(r.t - target_t))
            summary_records.append(
                {
                    "dt_s": dt,
                    "checkpoint_T": checkpoint_T,
                    "time_s": closest_record.t,
                    "position_error_km": closest_record.position_error,
                    "velocity_error_km_s": closest_record.velocity_error,
                    "energy_error": closest_record.energy_error,
                    "angular_momentum_error": closest_record.angular_momentum_error,
                }
            )

    return summary_records


def save_checkpoint_summary_to_csv(
    summary_records: List[Dict[str, float]],
    filepath: Path,
) -> None:
    """Save diagnostic checkpoint summary table to a CSV file.

    Args:
        summary_records: List of dictionary records created by extract_checkpoint_summary.
        filepath: Destination CSV Path.
    """
    filepath.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "dt_s",
        "checkpoint_T",
        "time_s",
        "position_error_km",
        "velocity_error_km_s",
        "energy_error",
        "angular_momentum_error",
    ]

    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=fieldnames, lineterminator="\n"
        )
        writer.writeheader()
        for rec in summary_records:
            writer.writerow(rec)


def print_long_term_summary(
    summary_records: List[Dict[str, float]],
    orbital_period: float,
) -> None:
    """Print formatted summary tables showing error progression across orbital checkpoints.

    Args:
        summary_records: List of diagnostic checkpoint dictionary records.
        orbital_period: Analytical orbital period T in seconds.
    """
    print("=" * 115)
    print("OrbitSim M3: Long-Term Stability of Explicit Euler (100 Orbits)")
    print(f"Orbital Period T = {orbital_period:.6f} s (~{orbital_period / 60.0:.2f} min)")
    print("=" * 115)

    current_dt = None
    for rec in summary_records:
        dt = rec["dt_s"]
        if dt != current_dt:
            current_dt = dt
            print(f"\n--- Timestep dt = {dt:.3f} s ---")
            headers = [
                "Checkpoint",
                "t (s)",
                "Pos Err (km)",
                "Vel Err (km/s)",
                "Energy Err",
                "AngMom Err",
            ]
            print(
                f"{headers[0]:>10} | {headers[1]:>12} | {headers[2]:>16} | "
                f"{headers[3]:>16} | {headers[4]:>14} | {headers[5]:>14}"
            )
            print("-" * 92)

        checkpoint_str = f"{rec['checkpoint_T']:.0f}T"
        print(
            f"{checkpoint_str:>10} | {rec['time_s']:>12.3f} | "
            f"{rec['position_error_km']:>16.6e} | "
            f"{rec['velocity_error_km_s']:>16.6e} | "
            f"{rec['energy_error']:>14.6e} | "
            f"{rec['angular_momentum_error']:>14.6e}"
        )

    print("\n" + "=" * 115)


def main() -> None:
    """Execute M3 long-term stability experiment, display summary, and export CSVs."""
    print("Running M3 Explicit Euler long-term integration (100 orbits)...")
    experiment_histories, t_orbit = run_euler_long_term_experiment()

    summary_records = extract_checkpoint_summary(experiment_histories, orbital_period=t_orbit)
    print_long_term_summary(summary_records, orbital_period=t_orbit)

    output_trajectory_csv = Path("results/exp03_euler_long_term.csv")
    save_long_term_trajectory_to_csv(experiment_histories, output_trajectory_csv)
    print(f"\nFull trajectory history saved to: {output_trajectory_csv.resolve()}")

    output_summary_csv = Path("results/exp03_euler_long_term_summary.csv")
    save_checkpoint_summary_to_csv(summary_records, output_summary_csv)
    print(f"Checkpoint summary table saved to: {output_summary_csv.resolve()}")


if __name__ == "__main__":
    main()
