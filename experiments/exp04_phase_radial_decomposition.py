"""Experiment M4: Phase Error and Radial Drift Decomposition.

Investigates whether long-term Euclidean position error of the Explicit Euler
integrator over multi-orbit time scales (0 <= t <= 100T) is dominated by:
1. Radial drift e_R(t) / Delta r(t)
2. Along-track phase error e_T(t) / Delta theta(t)
3. Or a transition between regimes.
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


def run_phase_radial_decomposition_experiment(
    dt_values: Optional[List[float]] = None,
    num_orbits: float = DEFAULT_NUM_ORBITS,
    altitude_km: float = DEFAULT_ALTITUDE_KM,
    mu: float = MU_EARTH,
    r_earth: float = R_EARTH,
) -> Tuple[Dict[float, List[DecomposedErrorRecord]], float]:
    """Run Explicit Euler phase error and radial drift decomposition experiment over 100T.

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
        sim = Simulator(initial_state=initial_state, mu=mu)
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


def save_decomposed_trajectory_to_csv(
    experiment_histories: Dict[float, List[DecomposedErrorRecord]],
    filepath: Path,
) -> None:
    """Save full time evolution decomposed histories for all timesteps to CSV.

    Columns 1-11 are backwards compatible with exp03_euler_long_term.csv.
    Columns 12-19 contain the decomposed metrics.

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


def extract_decomposed_checkpoint_summary(
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


def save_decomposed_summary_to_csv(
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


def main() -> None:
    """Execute M4 Phase Error and Radial Drift Decomposition experiment."""
    results_dir = PROJECT_ROOT / "results"
    traj_filepath = results_dir / "exp04_phase_radial_decomposition.csv"
    summary_filepath = results_dir / "exp04_phase_radial_summary.csv"

    print("=" * 70)
    print("OrbitSim — Experiment M4: Phase Error and Radial Drift Decomposition")
    print("=" * 70)
    print(f"Investigated timesteps dt (s): {DEFAULT_DT_VALUES}")
    print(f"Integration duration: {DEFAULT_NUM_ORBITS} orbits (100T)")
    print(f"Checkpoints: {DEFAULT_CHECKPOINTS_T} T")
    print("-" * 70)

    print("Running simulations and computing error decompositions...")
    histories, t_orbit = run_phase_radial_decomposition_experiment()

    print(f"Saving full decomposed trajectory to: {traj_filepath}")
    save_decomposed_trajectory_to_csv(histories, traj_filepath)

    print(f"Extracting checkpoint summaries and saving to: {summary_filepath}")
    summaries = extract_decomposed_checkpoint_summary(histories, t_orbit)
    save_decomposed_summary_to_csv(summaries, summary_filepath)

    print("\nCheckpoint Summary Results:")
    print("-" * 115)
    print(
        f"{'dt (s)':>8} | {'ckpt':>5} | {'e_r (km)':>12} | {'e_R (km)':>12} | {'e_T (km)':>12} | "
        f"{'rho_R':>7} | {'rho_T':>7} | {'chi':>7} | {'mode':>8}"
    )
    print("-" * 115)
    for s in summaries:
        print(
            f"{s['dt_s']:8.3f} | {s['checkpoint_T']:4.0f}T | {s['position_error_km']:12.2f} | "
            f"{s['e_R_km']:12.2f} | {s['e_T_km']:12.2f} | {s['fraction_R']:7.4f} | "
            f"{s['fraction_T']:7.4f} | {s['dominance_ratio']:7.2f} | {s['dominant_mode']:>8}"
        )
    print("-" * 115)
    print("Experiment M4 completed successfully.")


if __name__ == "__main__":
    main()
