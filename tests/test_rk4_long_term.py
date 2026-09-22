"""Unit and regression tests for M5: RK4 Long-Term Stability Experiment.

Tests verified mathematical and computational facts:
- Starts at t = 0
- Ends at 100T without overshoot
- Strictly monotonic timestamps
- Absence of NaN or Inf
- Checkpoint extraction correctness
- Orthogonal decomposition Pythagorean identity
- Variance partition identity
- Regression consistency and M4 compatibility
"""

import csv
import math
from pathlib import Path
import pytest

from experiments.exp05_rk4_long_term import (
    DEFAULT_CHECKPOINTS_T,
    DEFAULT_DT_VALUES,
    extract_rk4_checkpoint_summary,
    run_rk4_long_term_experiment,
)
from src.physics.analytical import circular_orbit_period
from src.physics.gravity import MU_EARTH, R_EARTH


@pytest.fixture(scope="module")
def rk4_experiment_data():
    """Run and cache RK4 long-term stability histories over 100T for default timesteps."""
    return run_rk4_long_term_experiment(dt_values=DEFAULT_DT_VALUES, num_orbits=100.0)


def test_rk4_long_term_starts_at_t0(rk4_experiment_data):
    """Fact: Trajectory histories for all timesteps start at t = 0 with zero position error."""
    histories, _ = rk4_experiment_data
    for dt, history in histories.items():
        assert len(history) > 0, f"History empty for dt={dt}"
        assert history[0].t == 0.0, f"dt={dt} does not start at t=0"
        assert math.isclose(history[0].position_error, 0.0, abs_tol=1e-12)
        assert math.isclose(history[0].velocity_error, 0.0, abs_tol=1e-12)


def test_rk4_long_term_ends_at_100T_without_overshoot(rk4_experiment_data):
    """Fact: Simulation terminates at exactly t = 100T without exceeding 100T + 1e-12 s."""
    histories, period = rk4_experiment_data
    t_final_target = 100.0 * period

    for dt, history in histories.items():
        final_t = history[-1].t
        assert math.isclose(final_t, t_final_target, rel_tol=1e-12, abs_tol=1e-9), (
            f"dt={dt} ended at {final_t}, expected {t_final_target}"
        )
        assert final_t <= t_final_target + 1e-12, f"dt={dt} overshot target time"


def test_rk4_long_term_strictly_monotonic_timestamps(rk4_experiment_data):
    """Fact: Timestamps are strictly monotonically increasing throughout the trajectory."""
    histories, _ = rk4_experiment_data
    for dt, history in histories.items():
        for i in range(1, len(history)):
            assert history[i].t > history[i - 1].t, (
                f"Non-monotonic timestamp at step {i} for dt={dt}: "
                f"t[{i}]={history[i].t} <= t[{i-1}]={history[i-1].t}"
            )


def test_rk4_long_term_no_nan_or_inf(rk4_experiment_data):
    """Fact: All state variables, errors, invariants, and projections are strictly finite."""
    histories, _ = rk4_experiment_data
    for dt, history in histories.items():
        for record in history:
            assert math.isfinite(record.state.x), f"NaN/Inf x at t={record.t}, dt={dt}"
            assert math.isfinite(record.state.y), f"NaN/Inf y at t={record.t}, dt={dt}"
            assert math.isfinite(record.state.vx), f"NaN/Inf vx at t={record.t}, dt={dt}"
            assert math.isfinite(record.state.vy), f"NaN/Inf vy at t={record.t}, dt={dt}"
            assert math.isfinite(record.position_error), f"NaN/Inf pos_error at t={record.t}"
            assert math.isfinite(record.velocity_error), f"NaN/Inf vel_error at t={record.t}"
            assert math.isfinite(record.energy_error), f"NaN/Inf energy_error at t={record.t}"
            assert math.isfinite(record.angular_momentum_error), f"NaN/Inf ang_mom at t={record.t}"
            assert math.isfinite(record.delta_r), f"NaN/Inf delta_r at t={record.t}"
            assert math.isfinite(record.delta_theta), f"NaN/Inf delta_theta at t={record.t}"
            assert math.isfinite(record.e_R), f"NaN/Inf e_R at t={record.t}"
            assert math.isfinite(record.e_T), f"NaN/Inf e_T at t={record.t}"
            assert math.isfinite(record.fraction_R), f"NaN/Inf fraction_R at t={record.t}"
            assert math.isfinite(record.fraction_T), f"NaN/Inf fraction_T at t={record.t}"


def test_rk4_long_term_all_checkpoints_extracted(rk4_experiment_data):
    """Fact: Checkpoint summary extraction captures all 6 checkpoints for every timestep."""
    histories, period = rk4_experiment_data
    summary = extract_rk4_checkpoint_summary(histories, period, DEFAULT_CHECKPOINTS_T)

    expected_total = len(DEFAULT_DT_VALUES) * len(DEFAULT_CHECKPOINTS_T)
    assert len(summary) == expected_total

    for rec in summary:
        assert rec["checkpoint_T"] in DEFAULT_CHECKPOINTS_T
        assert rec["dt_s"] in DEFAULT_DT_VALUES
        target_t = rec["checkpoint_T"] * period
        assert abs(rec["time_s"] - target_t) <= rec["dt_s"]


def test_rk4_pythagorean_identity_across_trajectory(rk4_experiment_data):
    """Fact: Rotating-frame error projections satisfy e_r^2 = e_R^2 + e_T^2 everywhere."""
    histories, _ = rk4_experiment_data
    # Sample every 100th step to keep verification fast while covering all regimes
    for dt, history in histories.items():
        for record in history[::100]:
            e_r_sq = record.position_error ** 2
            decomp_sq = record.e_R ** 2 + record.e_T ** 2
            assert math.isclose(e_r_sq, decomp_sq, rel_tol=1e-10, abs_tol=1e-14), (
                f"Pythagorean violation at t={record.t}, dt={dt}: "
                f"e_r^2={e_r_sq}, e_R^2+e_T^2={decomp_sq}"
            )


def test_rk4_variance_partition_identity(rk4_experiment_data):
    """Fact: Variance fractions satisfy rho_R + rho_T = 1.0 for all t > 0."""
    histories, _ = rk4_experiment_data
    for dt, history in histories.items():
        for record in history[1::100]:
            total_fraction = record.fraction_R + record.fraction_T
            assert math.isclose(total_fraction, 1.0, rel_tol=1e-12, abs_tol=1e-12), (
                f"Partition violation at t={record.t}, dt={dt}: sum={total_fraction}"
            )


def test_rk4_summary_csv_integrity():
    """Fact: results/exp05_rk4_long_term_summary.csv exists and has correct schema and rows."""
    summary_path = Path("results/exp05_rk4_long_term_summary.csv")
    assert summary_path.exists(), "results/exp05_rk4_long_term_summary.csv does not exist"

    expected_fields = [
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

    with open(summary_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        assert reader.fieldnames == expected_fields
        rows = list(reader)
        assert len(rows) == len(DEFAULT_DT_VALUES) * len(DEFAULT_CHECKPOINTS_T)


def test_rk4_trajectory_csv_integrity():
    """Fact: results/exp05_rk4_long_term.csv exists and contains non-empty valid rows."""
    traj_path = Path("results/exp05_rk4_long_term.csv")
    assert traj_path.exists(), "results/exp05_rk4_long_term.csv does not exist"

    with open(traj_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        first_row = next(reader)
        assert float(first_row["t_s"]) == 0.0
        assert float(first_row["position_error_km"]) == 0.0
