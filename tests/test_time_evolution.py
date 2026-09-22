"""Unit and regression tests for M2: Time Evolution of Explicit Euler Error.

Verifies mathematical consistency of time histories, boundary conditions,
monotonicity, regression against M1 results, and independent diagnostic evaluations.
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
from src.simulation.simulator import Simulator
from src.simulation.metrics import (
    ErrorRecord,
    compute_error_history,
    position_error,
    velocity_error,
)
from experiments.exp02_euler_time_evolution import (
    DEFAULT_DT_VALUES,
    run_euler_time_evolution_experiment,
)


@pytest.fixture(scope="module")
def m2_histories():
    """Compute and cache M2 error histories for all default timesteps."""
    return run_euler_time_evolution_experiment(dt_values=DEFAULT_DT_VALUES)


def test_initial_error_at_t0_is_zero():
    """Test 1: At t = 0, position and velocity errors must be zero within float precision."""
    r0 = R_EARTH + 400.0
    initial_state = create_circular_orbit_state(altitude_km=400.0, mu=MU_EARTH, r_earth=R_EARTH)
    exact_state_t0 = circular_orbit_exact_state(0.0, r0=r0, mu=MU_EARTH)

    pos_err = position_error(initial_state, exact_state_t0)
    vel_err = velocity_error(initial_state, exact_state_t0)

    assert math.isclose(pos_err, 0.0, abs_tol=1e-12)
    assert math.isclose(vel_err, 0.0, abs_tol=1e-12)


def test_time_history_starts_at_t0(m2_histories):
    """Test 2: First entry in error history must be at t = 0 with zero errors."""
    for dt, history in m2_histories.items():
        first_record = history[0]
        assert first_record.t == 0.0
        assert math.isclose(first_record.position_error, 0.0, abs_tol=1e-12)
        assert math.isclose(first_record.velocity_error, 0.0, abs_tol=1e-12)
        assert math.isclose(first_record.energy_error, 0.0, abs_tol=1e-12)
        assert math.isclose(first_record.angular_momentum_error, 0.0, abs_tol=1e-12)


def test_time_history_ends_at_period_without_overshoot(m2_histories):
    """Test 3: Final entry must be at t = T within tolerance, with no t > T."""
    r0 = R_EARTH + 400.0
    t_period = circular_orbit_period(r0, mu=MU_EARTH)

    for dt, history in m2_histories.items():
        # Check termination time
        last_record = history[-1]
        assert math.isclose(last_record.t, t_period, abs_tol=1e-11)

        # Strict no-overshoot check across all steps
        for rec in history:
            assert rec.t <= t_period + 1e-12


def test_time_history_strictly_monotonic(m2_histories):
    """Test 4: Timestamps must be strictly increasing: t_0 < t_1 < ... < t_N."""
    for dt, history in m2_histories.items():
        assert len(history) > 1
        for i in range(1, len(history)):
            assert history[i].t > history[i - 1].t


def test_regression_against_m1_final_errors(m2_histories):
    """Test 5: Final position and velocity errors at t = T must match M1 baseline values."""
    m1_csv_path = Path(__file__).resolve().parent.parent / "results" / "exp01_euler_convergence.csv"
    assert m1_csv_path.exists(), "M1 baseline results CSV must exist for regression verification."

    m1_data = {}
    with open(m1_csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            dt = float(row["dt_s"])
            m1_data[dt] = {
                "final_position_error_km": float(row["final_position_error_km"]),
                "final_velocity_error_km_s": float(row["final_velocity_error_km_s"]),
                "steps": int(row["steps"]),
            }

    for dt, history in m2_histories.items():
        assert dt in m1_data, f"dt={dt} must be present in M1 results."
        last_record = history[-1]
        expected_pos_err = m1_data[dt]["final_position_error_km"]
        expected_vel_err = m1_data[dt]["final_velocity_error_km_s"]
        expected_steps = m1_data[dt]["steps"]

        # Number of integration steps is len(history) - 1
        assert len(history) - 1 == expected_steps
        assert math.isclose(last_record.position_error, expected_pos_err, rel_tol=1e-12)
        assert math.isclose(last_record.velocity_error, expected_vel_err, rel_tol=1e-12)


def test_independent_verification_of_diagnostics(m2_histories):
    """Test 6: Independently recalculate errors and diagnostics using raw formulas.

    Verifies intermediate checkpoints without calling the production metric functions,
    ensuring independent validation of position error, velocity error, energy error,
    and angular momentum error.
    """
    r0 = R_EARTH + 400.0
    mu = MU_EARTH
    n_mean = math.sqrt(mu / (r0**3))
    epsilon0 = -mu / (2.0 * r0)
    h0 = r0 * math.sqrt(mu / r0)

    # Sample intermediate steps from dt = 60.0 and dt = 30.0 runs
    sample_indices = [0, 5, 20, 50, -1]

    for dt in [60.0, 30.0]:
        history = m2_histories[dt]
        for idx in sample_indices:
            rec = history[idx]
            t = rec.t
            s = rec.state

            # Independent analytical state computation
            nt = n_mean * t
            x_exact = r0 * math.cos(nt)
            y_exact = r0 * math.sin(nt)
            vx_exact = -r0 * n_mean * math.sin(nt)
            vy_exact = r0 * n_mean * math.cos(nt)

            # Independent position & velocity error
            expected_pos_err = math.sqrt((s.x - x_exact) ** 2 + (s.y - y_exact) ** 2)
            expected_vel_err = math.sqrt(
                (s.vx - vx_exact) ** 2 + (s.vy - vy_exact) ** 2
            )

            # Independent specific energy & relative error
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
