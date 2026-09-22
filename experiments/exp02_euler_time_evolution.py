"""Experiment M2: Time Evolution of Explicit Euler Error.

Investigates how numerical error and physical diagnostic invariants evolve over
time throughout one complete orbital period across multiple timesteps dt:
    error = f(t, dt) for t in [0, T]
"""

import csv
import math
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional

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


def run_euler_time_evolution_experiment(
    dt_values: Optional[List[float]] = None,
    altitude_km: float = DEFAULT_ALTITUDE_KM,
    mu: float = MU_EARTH,
    r_earth: float = R_EARTH,
) -> Dict[float, List[ErrorRecord]]:
    """Run Explicit Euler time evolution study over one full orbital period.

    For each timestep dt in dt_values, integrates the circular orbit from
    t = 0 to t = T (where T is analytical period) using step shortening at the end
    dt_step = min(dt, T - t) to prevent overshoot. Computes the complete error
    and invariant history at every discrete step.

    Args:
        dt_values: Timesteps in seconds to evaluate.
        altitude_km: Orbital altitude in km.
        mu: Gravitational parameter in km^3/s^2.
        r_earth: Central body radius in km.

    Returns:
        Dictionary mapping timestep dt -> List[ErrorRecord].
    """
    if dt_values is None:
        dt_values = DEFAULT_DT_VALUES

    r0 = r_earth + altitude_km
    t_final = circular_orbit_period(r0, mu=mu)
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

    return experiment_histories


def save_time_evolution_to_csv(
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


def print_time_evolution_summary(
    experiment_histories: Dict[float, List[ErrorRecord]],
    t_final: float,
) -> None:
    """Print formatted summary tables showing error progression along the orbit.

    Args:
        experiment_histories: Mapping of dt -> List[ErrorRecord].
        t_final: Analytical orbital period in seconds.
    """
    print("=" * 115)
    print("OrbitSim M2: Time Evolution of Explicit Euler Error")
    print(f"Orbital Period T = {t_final:.6f} s (~{t_final / 60.0:.2f} min)")
    print("=" * 115)

    checkpoints_fraction = [0.0, 0.25, 0.5, 0.75, 1.0]

    for dt, history in experiment_histories.items():
        print(f"\n--- Timestep dt = {dt:.3f} s (Total Steps: {len(history) - 1}) ---")
        headers = [
            "Fraction T",
            "t (s)",
            "Pos Err (km)",
            "Vel Err (km/s)",
            "Energy Err",
            "AngMom Err",
        ]
        print(
            f"{headers[0]:>10} | {headers[1]:>10} | {headers[2]:>16} | "
            f"{headers[3]:>16} | {headers[4]:>14} | {headers[5]:>14}"
        )
        print("-" * 90)

        # Find closest record to each checkpoint fraction
        for frac in checkpoints_fraction:
            target_t = frac * t_final
            closest_record = min(history, key=lambda r: abs(r.t - target_t))
            print(
                f"{frac:>10.2f} | {closest_record.t:>10.3f} | "
                f"{closest_record.position_error:>16.6e} | "
                f"{closest_record.velocity_error:>16.6e} | "
                f"{closest_record.energy_error:>14.6e} | "
                f"{closest_record.angular_momentum_error:>14.6e}"
            )

    print("\n" + "=" * 115)


def main() -> None:
    """Execute M2 time evolution experiment, display summary, and export CSV."""
    r0 = R_EARTH + DEFAULT_ALTITUDE_KM
    t_final = circular_orbit_period(r0, mu=MU_EARTH)

    experiment_histories = run_euler_time_evolution_experiment()
    print_time_evolution_summary(experiment_histories, t_final=t_final)

    output_csv = Path("results/exp02_euler_time_evolution.csv")
    save_time_evolution_to_csv(experiment_histories, output_csv)
    print(f"\nFull time evolution history saved to: {output_csv.resolve()}")


if __name__ == "__main__":
    main()
