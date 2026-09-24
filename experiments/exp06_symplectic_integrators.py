"""Experiment M6: Symplectic Integrators and Long-Term Energy Preservation.

Investigates whether the symplectic Velocity Verlet integrator preserves long-term orbital
structure, invariants (energy, angular momentum), radial stability, and along-track phase
better than non-symplectic integrators (Explicit Euler and classical RK4) over 100 orbital
periods across multiple timesteps dt in [1.875, 60.0] seconds.
"""

import csv
import math
from pathlib import Path
import sys
import time as pytime
from typing import Any, Dict, List, Optional, Tuple, Type, Union

# Ensure project root is in sys.path for direct script execution
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.integrators.euler import EulerIntegrator
from src.integrators.rk4 import RK4Integrator
from src.integrators.velocity_verlet import VelocityVerletIntegrator
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

REQUIRED_CSV_FIELDNAMES: List[str] = [
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
]


def run_integrator_long_term(
    integrator_cls: Union[
        Type[EulerIntegrator],
        Type[RK4Integrator],
        Type[VelocityVerletIntegrator],
    ],
    dt_values: Optional[List[float]] = None,
    num_orbits: float = DEFAULT_NUM_ORBITS,
    altitude_km: float = DEFAULT_ALTITUDE_KM,
    mu: float = MU_EARTH,
    r_earth: float = R_EARTH,
) -> Tuple[Dict[float, List[DecomposedErrorRecord]], float]:
    """Run long-term stability and decomposed error evaluation for a given integrator.

    Args:
        integrator_cls: Integrator class (EulerIntegrator, RK4Integrator, or VelocityVerletIntegrator).
        dt_values: List of timesteps in seconds to evaluate.
        num_orbits: Number of orbital periods T to integrate (default: 100.0).
        altitude_km: Orbital altitude in km (default: 400.0).
        mu: Gravitational parameter in km^3/s^2.
        r_earth: Equatorial Earth radius in km.

    Returns:
        Tuple of (histories_by_dt, orbital_period_T), where histories_by_dt maps dt -> List[DecomposedErrorRecord].
    """
    if dt_values is None:
        dt_values = DEFAULT_DT_VALUES

    r0 = r_earth + altitude_km
    t_orbit = circular_orbit_period(r0, mu=mu)
    t_final = num_orbits * t_orbit

    initial_state = create_circular_orbit_state(
        altitude_km=altitude_km, mu=mu, r_earth=r_earth
    )

    histories: Dict[float, List[DecomposedErrorRecord]] = {}

    for dt in dt_values:
        sim = Simulator(initial_state=initial_state, mu=mu, integrator=integrator_cls)
        trajectory = sim.run(dt=dt, t_final=t_final)
        times = sim.times

        history = compute_decomposed_error_history(
            times=times,
            trajectory=trajectory,
            r0=r0,
            mu=mu,
        )
        histories[dt] = history

    return histories, t_orbit


def run_m6_full_experiment(
    dt_values: Optional[List[float]] = None,
    num_orbits: float = DEFAULT_NUM_ORBITS,
    altitude_km: float = DEFAULT_ALTITUDE_KM,
    mu: float = MU_EARTH,
    r_earth: float = R_EARTH,
) -> Tuple[Dict[str, Dict[float, List[DecomposedErrorRecord]]], float]:
    """Run identical long-term simulations across Euler, RK4, and Velocity Verlet.

    Returns:
        Tuple of (all_histories, orbital_period_T), where all_histories maps
        method_name -> {dt -> List[DecomposedErrorRecord]}.
    """
    if dt_values is None:
        dt_values = DEFAULT_DT_VALUES

    method_configs = [
        ("Euler", EulerIntegrator),
        ("RK4", RK4Integrator),
        ("Velocity Verlet", VelocityVerletIntegrator),
    ]

    all_histories: Dict[str, Dict[float, List[DecomposedErrorRecord]]] = {}
    t_orbit = 0.0

    for method_name, integrator_cls in method_configs:
        start_t = pytime.perf_counter()
        histories, t_orbit = run_integrator_long_term(
            integrator_cls=integrator_cls,
            dt_values=dt_values,
            num_orbits=num_orbits,
            altitude_km=altitude_km,
            mu=mu,
            r_earth=r_earth,
        )
        elapsed = pytime.perf_counter() - start_t
        print(f"[{method_name}] Completed 100T across {len(dt_values)} timesteps in {elapsed:.2f} s")
        all_histories[method_name] = histories

    return all_histories, t_orbit


def extract_m6_checkpoint_summary(
    all_histories: Dict[str, Dict[float, List[DecomposedErrorRecord]]],
    orbital_period: float,
    checkpoints_T: Optional[List[float]] = None,
) -> List[Dict[str, Any]]:
    """Extract metrics at discrete checkpoints (1T, 5T, 10T, 20T, 50T, 100T) for all methods.

    Returns:
        List of dictionaries with exact required column schema.
    """
    if checkpoints_T is None:
        checkpoints_T = DEFAULT_CHECKPOINTS_T

    summary_records: List[Dict[str, Any]] = []

    for method_name, histories_by_dt in all_histories.items():
        for dt, history in histories_by_dt.items():
            for checkpoint_T in checkpoints_T:
                target_t = checkpoint_T * orbital_period
                closest_record = min(history, key=lambda r: abs(r.t - target_t))
                chi = dominance_ratio(closest_record.e_R, closest_record.e_T)

                summary_records.append(
                    {
                        "method": method_name,
                        "dt": dt,
                        "orbit_count": checkpoint_T,
                        "time": closest_record.t,
                        "e_r": closest_record.position_error,
                        "e_v": closest_record.velocity_error,
                        "e_R": closest_record.e_R,
                        "e_T": closest_record.e_T,
                        "Δr": closest_record.delta_r,
                        "Δθ": closest_record.delta_theta,
                        "ρ_R": closest_record.fraction_R,
                        "ρ_T": closest_record.fraction_T,
                        "χ": chi,
                        "δε": closest_record.energy_error,
                        "δh": closest_record.angular_momentum_error,
                    }
                )

    return summary_records


def save_m6_summary_to_csv(
    summary_records: List[Dict[str, Any]],
    filepath: Path,
) -> None:
    """Save checkpoint summary records to CSV file using required column names."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=REQUIRED_CSV_FIELDNAMES, lineterminator="\n"
        )
        writer.writeheader()
        for record in summary_records:
            writer.writerow(record)


def save_m6_trajectories_to_csv(
    all_histories: Dict[str, Dict[float, List[DecomposedErrorRecord]]],
    orbital_period: float,
    filepath: Path,
) -> None:
    """Save full step trajectory dataset for all methods and timesteps to CSV.

    Uses exact required column names and UTF-8 encoding.
    """
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=REQUIRED_CSV_FIELDNAMES, lineterminator="\n"
        )
        writer.writeheader()

        for method_name, histories_by_dt in all_histories.items():
            for dt, history in histories_by_dt.items():
                for record in history:
                    chi = dominance_ratio(record.e_R, record.e_T)
                    writer.writerow(
                        {
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
                        }
                    )


def main() -> None:
    """Execute M6 full experiment pipeline and generate datasets."""
    print("=" * 70)
    print("OrbitSim — Milestone M6: Symplectic Integrator Pipeline")
    print("=" * 70)

    results_dir = PROJECT_ROOT / "results"
    summary_path = results_dir / "exp06_symplectic_summary.csv"
    trajectories_path = results_dir / "exp06_symplectic_integrators.csv"

    print("\nPhase 4: Running 100T simulations for Euler, RK4, and Velocity Verlet...")
    all_histories, t_orbit = run_m6_full_experiment()

    print(f"\nPhase 5: Extracting checkpoint summaries at 1T, 5T, 10T, 20T, 50T, 100T...")
    summary_records = extract_m6_checkpoint_summary(all_histories, orbital_period=t_orbit)

    print(f"Writing summary dataset to {summary_path} ({len(summary_records)} rows)...")
    save_m6_summary_to_csv(summary_records, summary_path)

    print(f"Writing full trajectory dataset to {trajectories_path}...")
    save_m6_trajectories_to_csv(all_histories, orbital_period=t_orbit, filepath=trajectories_path)

    print("\nDataset generation completed successfully.")
    print("=" * 70)


if __name__ == "__main__":
    main()
