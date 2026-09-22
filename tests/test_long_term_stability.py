"""Unit and regression tests for M3: Long-Term Stability of Explicit Euler.

Verifies verified empirical facts, trajectory time boundaries over 100 orbits,
non-NaN/Inf numerical state integrity, checkpoint existence, independent formula
evaluations, and regression consistency against M2 single-orbit results.
"""

import csv
import math
from pathlib import Path
import pytest

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
from src.simulation.metrics import (
    position_error,
    velocity_error,
)
from experiments.exp03_euler_long_term import (
    DEFAULT_DT_VALUES,
    DEFAULT_CHECKPOINTS_T,
    run_euler_long_term_experiment,
    extract_checkpoint_summary,
)


@pytest.fixture(scope="module")
def m3_experiment_data():
    """Compute and cache M3 error histories and period for default timesteps over 100T."""
    return run_euler_long_term_experiment(dt_values=DEFAULT_DT_VALUES, num_orbits=100.0)


# ==============================================================================
# 1. TESTED FACTS
# ==============================================================================

def test_initial_error_at_t0_is_zero():
    """Fact 1: At t = 0, initial state error relative to exact circular orbit must be zero."""
    r0 = R_EARTH + 400.0
    initial_state = create_circular_orbit_state(altitude_km=400.0, mu=MU_EARTH, r_earth=R_EARTH)
    exact_state_t0 = circular_orbit_exact_state(0.0, r0=r0, mu=MU_EARTH)

    pos_err = position_error(initial_state, exact_state_t0)
    vel_err = velocity_error(initial_state, exact_state_t0)

    assert math.isclose(pos_err, 0.0, abs_tol=1e-12)
    assert math.isclose(vel_err, 0.0, abs_tol=1e-12)


def test_long_term_history_starts_at_t0(m3_experiment_data):
    """Fact 2: First entry in trajectory history for all runs must be at t = 0 with zero error."""
    histories, _ = m3_experiment_data
    for dt, history in histories.items():
        first_rec = history[0]
        assert first_rec.t == 0.0
        assert math.isclose(first_rec.position_error, 0.0, abs_tol=1e-12)
        assert math.isclose(first_rec.velocity_error, 0.0, abs_tol=1e-12)
        assert math.isclose(first_rec.energy_error, 0.0, abs_tol=1e-12)
        assert math.isclose(first_rec.angular_momentum_error, 0.0, abs_tol=1e-12)


def test_long_term_history_ends_at_100T_without_overshoot(m3_experiment_data):
    """Fact 3: Final entry must terminate precisely at t = 100 * T with no overshoot."""
    histories, t_orbit = m3_experiment_data
    t_100T = 100.0 * t_orbit

    for dt, history in histories.items():
        last_rec = history[-1]
        assert math.isclose(last_rec.t, t_100T, abs_tol=1e-11)

        # No step exceeds t_100T
        for rec in history:
            assert rec.t <= t_100T + 1e-12


def test_long_term_history_strictly_monotonic(m3_experiment_data):
    """Fact 4: Timestamps across all integration steps must be strictly increasing: t_{k+1} > t_k."""
    histories, _ = m3_experiment_data
    for dt, history in histories.items():
        assert len(history) > 1
        for i in range(1, len(history)):
            assert history[i].t > history[i - 1].t


def test_long_term_no_nan_or_inf_in_states_and_metrics(m3_experiment_data):
    """Fact 5: States and diagnostic metrics must contain no NaN or Inf floating point values."""
    histories, _ = m3_experiment_data
    for dt, history in histories.items():
        for rec in history:
            assert not math.isnan(rec.t) and not math.isinf(rec.t)
            assert not math.isnan(rec.state.x) and not math.isinf(rec.state.x)
            assert not math.isnan(rec.state.y) and not math.isinf(rec.state.y)
            assert not math.isnan(rec.state.vx) and not math.isinf(rec.state.vx)
            assert not math.isnan(rec.state.vy) and not math.isinf(rec.state.vy)
            assert not math.isnan(rec.position_error) and not math.isinf(rec.position_error)
            assert not math.isnan(rec.velocity_error) and not math.isinf(rec.velocity_error)
            assert not math.isnan(rec.energy_error) and not math.isinf(rec.energy_error)
            assert not math.isnan(rec.angular_momentum_error) and not math.isinf(rec.angular_momentum_error)


def test_long_term_all_checkpoints_exist(m3_experiment_data):
    """Fact 6: Diagnostic checkpoint summary must extract entries for all requested horizons (1T, 5T, 10T, 20T, 50T, 100T)."""
    histories, t_orbit = m3_experiment_data
    summary_records = extract_checkpoint_summary(histories, orbital_period=t_orbit)

    # Expected count: len(DEFAULT_DT_VALUES) * len(DEFAULT_CHECKPOINTS_T)
    expected_count = len(DEFAULT_DT_VALUES) * len(DEFAULT_CHECKPOINTS_T)
    assert len(summary_records) == expected_count

    for checkpoint in DEFAULT_CHECKPOINTS_T:
        for dt in DEFAULT_DT_VALUES:
            matching = [
                rec for rec in summary_records
                if rec["dt_s"] == dt and rec["checkpoint_T"] == checkpoint
            ]
            assert len(matching) == 1
            rec = matching[0]
            # Time must be close to checkpoint target time
            assert abs(rec["time_s"] - checkpoint * t_orbit) <= dt


def test_independent_verification_of_metrics(m3_experiment_data):
    """Fact 7: Independently recalculate state errors and invariants using raw mathematical formulas.

    Verifies sampled checkpoint states without using the production metrics module,
    confirming exact mathematical correctness of recorded diagnostic data.
    """
    histories, t_orbit = m3_experiment_data
    r0 = R_EARTH + 400.0
    mu = MU_EARTH
    n_mean = math.sqrt(mu / (r0**3))
    epsilon0 = -mu / (2.0 * r0)
    h0 = r0 * math.sqrt(mu / r0)

    # Check discrete records at index 0, 100, 1000, and final step for dt = 60.0 and dt = 1.875
    sample_indices = [0, 100, 1000, -1]

    for dt in [60.0, 1.875]:
        history = histories[dt]
        for idx in sample_indices:
            rec = history[idx]
            t = rec.t
            s = rec.state

            # Independent analytical reference state
            nt = n_mean * t
            x_exact = r0 * math.cos(nt)
            y_exact = r0 * math.sin(nt)
            vx_exact = -r0 * n_mean * math.sin(nt)
            vy_exact = r0 * n_mean * math.cos(nt)

            # Independent position & velocity Euclidean error
            expected_pos_err = math.sqrt((s.x - x_exact) ** 2 + (s.y - y_exact) ** 2)
            expected_vel_err = math.sqrt((s.vx - vx_exact) ** 2 + (s.vy - vy_exact) ** 2)

            # Independent specific orbital energy & relative error
            r = math.sqrt(s.x**2 + s.y**2)
            v_sq = s.vx**2 + s.vy**2
            energy = 0.5 * v_sq - mu / r
            expected_energy_err = (energy - epsilon0) / abs(epsilon0)

            # Independent specific angular momentum & relative error
            hz = s.x * s.vy - s.y * s.vx
            expected_ang_err = (hz - h0) / abs(h0)

            assert math.isclose(rec.position_error, expected_pos_err, abs_tol=1e-12)
            assert math.isclose(rec.velocity_error, expected_vel_err, abs_tol=1e-12)
            assert math.isclose(rec.energy_error, expected_energy_err, abs_tol=1e-12)
            assert math.isclose(rec.angular_momentum_error, expected_ang_err, abs_tol=1e-12)


# ==============================================================================
# 2. REGRESSION TESTS
# ==============================================================================

def test_regression_against_m2_trajectory_steps(m3_experiment_data):
    """Regression Test: Un-shortened trajectory steps in M3 must match M2 trajectory states to float precision."""
    m2_csv_path = Path(__file__).resolve().parent.parent / "results" / "exp02_euler_time_evolution.csv"
    assert m2_csv_path.exists(), "M2 baseline CSV must exist for regression verification."

    # Index M2 rows by (dt, step)
    m2_data = {}
    with open(m2_csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            dt = float(row["dt_s"])
            step = int(row["step"])
            m2_data[(dt, step)] = row

    histories, t_orbit = m3_experiment_data

    for dt, history in histories.items():
        # Check first 10 steps (which are strictly within the first period and un-shortened)
        for step_idx in range(10):
            key = (dt, step_idx)
            assert key in m2_data, f"Step {step_idx} for dt={dt} must exist in M2 data."
            m2_row = m2_data[key]
            m3_rec = history[step_idx]

            assert math.isclose(m3_rec.t, float(m2_row["t_s"]), abs_tol=1e-12)
            assert math.isclose(m3_rec.state.x, float(m2_row["x_km"]), abs_tol=1e-12)
            assert math.isclose(m3_rec.state.y, float(m2_row["y_km"]), abs_tol=1e-12)
            assert math.isclose(m3_rec.state.vx, float(m2_row["vx_km_s"]), abs_tol=1e-12)
            assert math.isclose(m3_rec.state.vy, float(m2_row["vy_km_s"]), abs_tol=1e-12)
            assert math.isclose(m3_rec.position_error, float(m2_row["position_error_km"]), abs_tol=1e-12)
            assert math.isclose(m3_rec.velocity_error, float(m2_row["velocity_error_km_s"]), abs_tol=1e-12)
            assert math.isclose(m3_rec.energy_error, float(m2_row["energy_error"]), abs_tol=1e-12)
            assert math.isclose(m3_rec.angular_momentum_error, float(m2_row["angular_momentum_error"]), abs_tol=1e-12)

