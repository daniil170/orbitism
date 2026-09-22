"""Experiment M5: Long-Term Stability of RK4.

Investigates how numerical error, physical diagnostic invariants, and rotating-frame
decomposed error components evolve over multi-orbit time scales (0 <= t <= 100T)
using the classical 4th-order Runge-Kutta (RK4) integrator across multiple timesteps dt:
    error = f(t, dt) for t in [0, 100T]
"""

import csv
import math
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional, Tuple

# Ensure project root is in sys.path for direct script execution
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.integrators.rk4 import RK4Integrator
from src.models.state import State2D
from src.physics.gravity import (
    MU_EARTH,
    R_EARTH,
    create_circular_orbit_state,
)
from src.physics.analytical import circular_orbit_period
from src.simulation.simulator import Simulator
from src.simulation.metrics import (
    DecomposedErrorRecord,
    compute_decomposed_error_history,
    dominance_ratio,
)

DEFAULT_DT_VALUES: List[float] = [60.0, 30.0, 15.0, 7.5, 3.75, 1.875]
DEFAULT_ALTITUDE_KM: float = 400.0
DEFAULT_NUM_ORBITS: float = 100.0
DEFAULT_CHECKPOINTS_T: List[float] = [1.0, 5.0, 10.0, 20.0, 50.0, 100.0]


def run_rk4_long_term_experiment(
    dt_values: Optional[List[float]] = None,
    num_orbits: float = DEFAULT_NUM_ORBITS,
    altitude_km: float = DEFAULT_ALTITUDE_KM,
    mu: float = MU_EARTH,
    r_earth: float = R_EARTH,
) -> Tuple[Dict[float, List[DecomposedErrorRecord]], float]:
    """Run RK4 long-term stability and error decomposition experiment over 100T.

    For each timestep dt in dt_values, integrates the circular orbit from
    t = 0 to t = 100 * T (where T is analytical period) using terminal step shortening
    dt_step = min(dt, t_final - t) to prevent overshoot. Computes the complete
    decomposed error and invariant history at every discrete step.

    Args:
        dt_values: Timesteps in seconds to evaluate.
        num_orbits: Integration duration in units of orbital period T (default: 100.0).
        altitude_km: Orbital altitude in km.
        mu: Gravitational parameter in km^3/s^2.
        r_earth: Central body radius in km.

    Returns:
        Tuple of (experiment_histories, orbital_period_T), where experiment_histories
        maps timestep dt -> List[DecomposedErrorRecord].
    """
    if dt_values is None:
        dt_values = DEFAULT_DT_VALUES

    r0 = r_earth + altitude_km
    t_orbit = circular_orbit_period(r0, mu=mu)
    t_final = num_orbits * t_orbit

    initial_state = create_circular_orbit_state(
        altitude_km=altitude_km, mu=mu, r_earth=r_earth
    )

    experiment_histories: Dict[float, List[DecomposedErrorRecord]] = {}

    for dt in dt_values:
        sim = Simulator(initial_state=initial_state, mu=mu, integrator=RK4Integrator)
        trajectory = sim.run(dt=dt, t_final=t_final)
        times = sim.times

        history = compute_decomposed_error_history(
            times=times,
            trajectory=trajectory,
            r0=r0,
            mu=mu,
        )
        experiment_histories[dt] = history

    return experiment_histories, t_orbit


def save_rk4_trajectory_to_csv(
    experiment_histories: Dict[float, List[DecomposedErrorRecord]],
    filepath: Path,
) -> None:
    """Save full time evolution decomposed histories for all timesteps to CSV.

    Columns 1-11 match the core metric schema of exp03/exp04.
    Columns 12-19 contain the decomposed orthogonal metrics and fractions.

    Args:
        experiment_histories: Mapping of dt -> List[DecomposedErrorRecord].
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
        "r_num_km",
        "delta_r_km",
        "theta_unwrapped_rad",
        "delta_theta_rad",
        "e_R_km",
        "e_T_km",
        "fraction_R",
        "fraction_T",
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
                        "r_num_km": record.r_num,
                        "delta_r_km": record.delta_r,
                        "theta_unwrapped_rad": record.theta_unwrapped,
                        "delta_theta_rad": record.delta_theta,
                        "e_R_km": record.e_R,
                        "e_T_km": record.e_T,
                        "fraction_R": record.fraction_R,
                        "fraction_T": record.fraction_T,
                    }
                )


def extract_rk4_checkpoint_summary(
    experiment_histories: Dict[float, List[DecomposedErrorRecord]],
    orbital_period: float,
    checkpoints_T: Optional[List[float]] = None,
) -> List[Dict[str, Any]]:
    """Extract decomposed metrics at discrete checkpoints (1T, 5T, 10T, 20T, 50T, 100T).

    Args:
        experiment_histories: Mapping of dt -> List[DecomposedErrorRecord].
        orbital_period: Analytical orbital period T in seconds.
        checkpoints_T: List of checkpoint multiples of T.

    Returns:
        List of dictionaries containing checkpoint metrics and dominance classifications.
    """
    if checkpoints_T is None:
        checkpoints_T = DEFAULT_CHECKPOINTS_T

    summary_records: List[Dict[str, Any]] = []

    for dt, history in experiment_histories.items():
        for checkpoint_T in checkpoints_T:
            target_t = checkpoint_T * orbital_period
            closest_record = min(history, key=lambda r: abs(r.t - target_t))

            chi = dominance_ratio(closest_record.e_R, closest_record.e_T)
            if closest_record.fraction_T > 0.5:
                dominant_mode = "phase"
            elif closest_record.fraction_R > 0.5:
                dominant_mode = "radial"
            else:
                dominant_mode = "equipartition"

            summary_records.append(
                {
                    "dt_s": dt,
                    "checkpoint_T": checkpoint_T,
                    "time_s": closest_record.t,
                    "position_error_km": closest_record.position_error,
                    "velocity_error_km_s": closest_record.velocity_error,
                    "energy_error": closest_record.energy_error,
                    "angular_momentum_error": closest_record.angular_momentum_error,
                    "delta_r_km": closest_record.delta_r,
                    "delta_theta_rad": closest_record.delta_theta,
                    "e_R_km": closest_record.e_R,
                    "e_T_km": closest_record.e_T,
                    "fraction_R": closest_record.fraction_R,
                    "fraction_T": closest_record.fraction_T,
                    "dominance_ratio": chi,
                    "dominant_mode": dominant_mode,
                }
            )

    return summary_records


def save_rk4_summary_to_csv(
    summary_records: List[Dict[str, Any]],
    filepath: Path,
) -> None:
    """Save decomposed checkpoint summary to CSV file.

    Args:
        summary_records: List of summary dictionaries.
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
        "delta_r_km",
        "delta_theta_rad",
        "e_R_km",
        "e_T_km",
        "fraction_R",
        "fraction_T",
        "dominance_ratio",
        "dominant_mode",
    ]

    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=fieldnames, lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(summary_records)


def print_rk4_long_term_summary(
    summary_records: List[Dict[str, Any]],
    orbital_period: float,
) -> None:
    """Print formatted summary tables showing error progression across orbital checkpoints.

    Args:
        summary_records: List of summary records created by extract_rk4_checkpoint_summary.
        orbital_period: Analytical orbital period T in seconds.
    """
    print("=" * 135)
    print("OrbitSim M5: Long-Term Stability of RK4 (100 Orbits)")
    print(f"Orbital Period T = {orbital_period:.6f} s (~{orbital_period / 60.0:.2f} min)")
    print("=" * 135)

    current_dt = None
    for rec in summary_records:
        dt = rec["dt_s"]
        if dt != current_dt:
            current_dt = dt
            print(f"\n--- Timestep dt = {dt:.3f} s ---")
            headers = [
                "Ckpt",
                "t (s)",
                "e_r (km)",
                "e_v (km/s)",
                "e_R (km)",
                "e_T (km)",
                "rho_R",
                "rho_T",
                "chi",
                "Mode",
                "Delta r (km)",
                "Delta theta (rad)",
            ]
            print(
                f"{headers[0]:>5} | {headers[1]:>10} | {headers[2]:>12} | {headers[3]:>12} | "
                f"{headers[4]:>12} | {headers[5]:>12} | {headers[6]:>6} | {headers[7]:>6} | "
                f"{headers[8]:>6} | {headers[9]:>6} | {headers[10]:>12} | {headers[11]:>17}"
            )
            print("-" * 135)

        ckpt_str = f"{rec['checkpoint_T']:.0f}T"
        print(
            f"{ckpt_str:>5} | {rec['time_s']:>10.1f} | "
            f"{rec['position_error_km']:>12.4e} | {rec['velocity_error_km_s']:>12.4e} | "
            f"{rec['e_R_km']:>12.4e} | {rec['e_T_km']:>12.4e} | "
            f"{rec['fraction_R']:>6.4f} | {rec['fraction_T']:>6.4f} | "
            f"{rec['dominance_ratio']:>6.2f} | {rec['dominant_mode']:>6} | "
            f"{rec['delta_r_km']:>12.4e} | {rec['delta_theta_rad']:>17.8e}"
        )

    print("\n" + "=" * 135)


def main() -> None:
    """Execute M5 RK4 long-term stability experiment, display summary, and export CSVs."""
    print("Running M5 RK4 long-term integration (100 orbits)...")
    experiment_histories, t_orbit = run_rk4_long_term_experiment()

    summary_records = extract_rk4_checkpoint_summary(
        experiment_histories, orbital_period=t_orbit
    )
    print_rk4_long_term_summary(summary_records, orbital_period=t_orbit)

    output_trajectory_csv = Path("results/exp05_rk4_long_term.csv")
    save_rk4_trajectory_to_csv(experiment_histories, output_trajectory_csv)
    print(f"\nFull trajectory history saved to: {output_trajectory_csv.resolve()}")

    output_summary_csv = Path("results/exp05_rk4_long_term_summary.csv")
    save_rk4_summary_to_csv(summary_records, output_summary_csv)
    print(f"Checkpoint summary table saved to: {output_summary_csv.resolve()}")


if __name__ == "__main__":
    main()
