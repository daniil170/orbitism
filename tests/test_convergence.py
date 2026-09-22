"""Unit tests for analytical circular orbit solutions and convergence metrics."""

import math
import pytest

from src.models.state import State2D
from src.physics.gravity import (
    MU_EARTH,
    R_EARTH,
    create_circular_orbit_state,
)
from src.physics.analytical import (
    circular_orbit_mean_motion,
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


def test_analytical_circular_solution_at_t0():
    """Analytical circular state at t=0 must match initial circular orbit state."""
    r0 = R_EARTH + 400.0
    initial_state = create_circular_orbit_state(altitude_km=400.0)
    exact_t0 = circular_orbit_exact_state(t=0.0, r0=r0, mu=MU_EARTH)

    assert math.isclose(exact_t0.x, initial_state.x, abs_tol=1e-12)
    assert math.isclose(exact_t0.y, initial_state.y, abs_tol=1e-12)
    assert math.isclose(exact_t0.vx, initial_state.vx, abs_tol=1e-12)
    assert math.isclose(exact_t0.vy, initial_state.vy, abs_tol=1e-12)


def test_analytical_circular_solution_at_period_returns_to_initial():
    """Analytical circular solution at t=T must return to initial conditions within float precision."""
    r0 = R_EARTH + 400.0
    period = circular_orbit_period(r0, mu=MU_EARTH)
    initial_state = create_circular_orbit_state(altitude_km=400.0)
    exact_t_final = circular_orbit_exact_state(t=period, r0=r0, mu=MU_EARTH)

    assert math.isclose(exact_t_final.x, initial_state.x, rel_tol=1e-12)
    assert math.isclose(exact_t_final.y, initial_state.y, abs_tol=1e-10)
    assert math.isclose(exact_t_final.vx, initial_state.vx, abs_tol=1e-10)
    assert math.isclose(exact_t_final.vy, initial_state.vy, rel_tol=1e-12)


def test_position_and_velocity_error_against_self():
    """Error of an exact state evaluated against itself must be zero."""
    state = State2D(x=6778.137, y=0.0, vx=0.0, vy=7.67)
    e_r = position_error(state, state)
    e_v = velocity_error(state, state)

    assert math.isclose(e_r, 0.0, abs_tol=1e-15)
    assert math.isclose(e_v, 0.0, abs_tol=1e-15)


def test_convergence_order_synthetic_case():
    """Convergence order helper must return exact theoretical slope for synthetic rates.

    For p=1: halving dt halves error (order = 1.0).
    For p=2: halving dt quarters error (order = 2.0).
    """
    # 1st order synthetic check
    p1 = convergence_order(error_1=2.0, error_2=1.0, dt_1=10.0, dt_2=5.0)
    assert math.isclose(p1, 1.0, rel_tol=1e-14)

    # 2nd order synthetic check
    p2 = convergence_order(error_1=4.0, error_2=1.0, dt_1=10.0, dt_2=5.0)
    assert math.isclose(p2, 2.0, rel_tol=1e-14)


def test_simulation_ends_exactly_at_t_final():
    """Simulator running with t_final must finish precisely at t_final."""
    initial_state = create_circular_orbit_state(altitude_km=400.0)
    sim = Simulator(initial_state=initial_state)

    t_final = 100.25  # Non-integer duration
    dt = 10.0
    sim.run(dt=dt, t_final=t_final)

    assert math.isclose(sim.times[-1], t_final, abs_tol=1e-12)
    assert sim.times[0] == 0.0


def test_simulation_shortens_last_timestep_when_not_multiple():
    """When duration is not an integer multiple of dt, the final step must be shortened."""
    initial_state = create_circular_orbit_state(altitude_km=400.0)
    sim = Simulator(initial_state=initial_state)

    t_final = 25.0
    dt = 10.0
    # Expected times: [0.0, 10.0, 20.0, 25.0]
    sim.run(dt=dt, t_final=t_final)

    assert len(sim.times) == 4
    assert math.isclose(sim.times[1] - sim.times[0], 10.0, abs_tol=1e-12)
    assert math.isclose(sim.times[2] - sim.times[1], 10.0, abs_tol=1e-12)
    assert math.isclose(sim.times[3] - sim.times[2], 5.0, abs_tol=1e-12)
    assert math.isclose(sim.times[-1], 25.0, abs_tol=1e-12)
