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
from src.physics.analytical import (
    circular_orbit_period,
    circular_orbit_exact_state,
    specific_orbital_energy,
    specific_angular_momentum,
)
from src.simulation.simulator import Simulator
from src.simulation.metrics import (
    position_error,
    velocity_error,
    convergence_order,
)

DEFAULT_DT_VALUES: List[float] = [60.0, 30.0, 15.0, 7.5, 3.75, 1.875]
DEFAULT_ALTITUDE_KM: float = 400.0


def run_euler_convergence_experiment(
    dt_values: Optional[List[float]] = None,
    altitude_km: float = DEFAULT_ALTITUDE_KM,
    mu: float = MU_EARTH,
    r_earth: float = R_EARTH,
) -> List[Dict[str, Any]]:
    """Run Explicit Euler convergence study over one full orbital period.

    Args:
        dt_values: Timesteps in seconds to evaluate.
        altitude_km: Orbital altitude in km.
        mu: Gravitational parameter in km^3/s^2.
        r_earth: Central body radius in km.

    Returns:
        List of result dictionaries containing summary metrics for each dt.
    """
    if dt_values is None:
        dt_values = DEFAULT_DT_VALUES

    r0 = r_earth + altitude_km
    t_final = circular_orbit_period(r0, mu=mu)
    initial_state = create_circular_orbit_state(
        altitude_km=altitude_km, mu=mu, r_earth=r_earth
    )

    ref_energy = specific_orbital_energy(initial_state, mu=mu)
    ref_angular_momentum = specific_angular_momentum(initial_state)

    results: List[Dict[str, Any]] = []

    for dt in dt_values:
        sim = Simulator(initial_state=initial_state, mu=mu)
        trajectory = sim.run(dt=dt, t_final=t_final)
        times = sim.times

        steps = len(trajectory) - 1

        max_pos_err = 0.0
        max_vel_err = 0.0
        max_energy_err = 0.0
        max_ang_mom_err = 0.0

        for t_k, s_k in zip(times, trajectory):
            exact_k = circular_orbit_exact_state(t_k, r0=r0, mu=mu)
            e_r = position_error(s_k, exact_k)
            e_v = velocity_error(s_k, exact_k)

            energy_k = specific_orbital_energy(s_k, mu=mu)
            rel_energy_err = (energy_k - ref_energy) / abs(ref_energy)

            h_k = specific_angular_momentum(s_k)
            rel_ang_err = (h_k - ref_angular_momentum) / abs(ref_angular_momentum)

            if e_r > max_pos_err:
                max_pos_err = e_r
            if e_v > max_vel_err:
                max_vel_err = e_v
            if abs(rel_energy_err) > max_energy_err:
                max_energy_err = abs(rel_energy_err)
            if abs(rel_ang_err) > max_ang_mom_err:
                max_ang_mom_err = abs(rel_ang_err)

        final_exact = circular_orbit_exact_state(t_final, r0=r0, mu=mu)
        final_state = trajectory[-1]
        final_pos_err = position_error(final_state, final_exact)
        final_vel_err = velocity_error(final_state, final_exact)

        final_energy = specific_orbital_energy(final_state, mu=mu)
        final_energy_err = (final_energy - ref_energy) / abs(ref_energy)

        final_h = specific_angular_momentum(final_state)
        final_ang_err = (final_h - ref_angular_momentum) / abs(ref_angular_momentum)

        row: Dict[str, Any] = {
            "dt_s": dt,
            "steps": steps,
            "final_position_error_km": final_pos_err,
            "max_position_error_km": max_pos_err,
            "final_velocity_error_km_s": final_vel_err,
            "max_velocity_error_km_s": max_vel_err,
            "position_order": None,
            "velocity_order": None,
            "max_energy_error": max_energy_err,
            "final_energy_error": final_energy_err,
            "max_angular_momentum_error": max_ang_mom_err,
            "final_angular_momentum_error": final_ang_err,
        }
        results.append(row)

    # Compute convergence orders between adjacent timesteps
    for i in range(1, len(results)):
        prev = results[i - 1]
        curr = results[i]
        curr["position_order"] = convergence_order(
            error_1=prev["max_position_error_km"],
            error_2=curr["max_position_error_km"],
            dt_1=prev["dt_s"],
            dt_2=curr["dt_s"],
        )
        curr["velocity_order"] = convergence_order(
            error_1=prev["max_velocity_error_km_s"],
            error_2=curr["max_velocity_error_km_s"],
            dt_1=prev["dt_s"],
            dt_2=curr["dt_s"],
        )

    return results


def save_results_to_csv(
    results: List[Dict[str, Any]], filepath: Path
) -> None:
    """Save experiment summary metrics to CSV file.

    Args:
        results: List of metric rows.
        filepath: Destination Path.
    """
    filepath.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "dt_s",
        "steps",
        "final_position_error_km",
        "max_position_error_km",
        "final_velocity_error_km_s",
        "max_velocity_error_km_s",
        "position_order",
        "velocity_order",
        "max_energy_error",
        "final_energy_error",
        "max_angular_momentum_error",
        "final_angular_momentum_error",
    ]

    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=fieldnames, lineterminator="\n"
        )
        writer.writeheader()

        for row in results:
            formatted_row = dict(row)
            if formatted_row["position_order"] is None:
                formatted_row["position_order"] = "N/A"
            if formatted_row["velocity_order"] is None:
                formatted_row["velocity_order"] = "N/A"
            writer.writerow(formatted_row)


def print_results_table(results: List[Dict[str, Any]], t_final: float) -> None:
    """Print formatted ASCII summary table of experiment results.

    Args:
        results: List of metric rows.
        t_final: Total simulation duration in seconds.
    """
    print("=" * 115)
    print("OrbitSim M1: Explicit Euler Convergence Experiment")
    print(f"Orbital Period T = {t_final:.6f} s (~{t_final / 60.0:.2f} min)")
    print("=" * 115)

    headers = [
        "dt (s)",
        "Steps",
        "Final Pos Err (km)",
        "Max Pos Err (km)",
        "Max Vel Err (km/s)",
        "Pos Ord",
        "Vel Ord",
        "Max Energy Err",
        "Max AngMom Err",
    ]
    header_line = (
        f"{headers[0]:>8} | {headers[1]:>6} | {headers[2]:>18} | "
        f"{headers[3]:>18} | {headers[4]:>18} | {headers[5]:>8} | "
        f"{headers[6]:>8} | {headers[7]:>14} | {headers[8]:>14}"
    )
    print(header_line)
    print("-" * 115)

    for r in results:
        pos_ord_str = (
            f"{r['position_order']:.4f}"
            if r["position_order"] is not None
            else "N/A"
        )
        vel_ord_str = (
            f"{r['velocity_order']:.4f}"
            if r["velocity_order"] is not None
            else "N/A"
        )

        print(
            f"{r['dt_s']:>8.3f} | {r['steps']:>6d} | "
            f"{r['final_position_error_km']:>18.6e} | "
            f"{r['max_position_error_km']:>18.6e} | "
            f"{r['max_velocity_error_km_s']:>18.6e} | "
            f"{pos_ord_str:>8} | {vel_ord_str:>8} | "
            f"{r['max_energy_error']:>14.6e} | {r['max_angular_momentum_error']:>14.6e}"
        )
    print("=" * 115)


def main() -> None:
    """Execute M1 convergence experiment, print table, and persist CSV results."""
    r0 = R_EARTH + DEFAULT_ALTITUDE_KM
    t_final = circular_orbit_period(r0, mu=MU_EARTH)

    results = run_euler_convergence_experiment()
    print_results_table(results, t_final=t_final)

    output_csv = Path("results/exp01_euler_convergence.csv")
    save_results_to_csv(results, output_csv)
    print(f"\nResults saved successfully to: {output_csv.resolve()}")


if __name__ == "__main__":
    main()
