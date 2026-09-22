"""Unit tests for initial conditions and circular orbit parameters."""

import math
from src.physics.gravity import (
    MU_EARTH,
    R_EARTH,
    circular_orbit_velocity,
    create_circular_orbit_state,
    gravitational_acceleration,
)


def test_circular_orbit_velocity_formula():
    """Initial circular orbit velocity must match v = sqrt(mu / r)."""
    r0 = R_EARTH + 400.0  # 6778.137 km
    expected_v = math.sqrt(MU_EARTH / r0)
    calculated_v = circular_orbit_velocity(r0, mu=MU_EARTH)

    assert math.isclose(calculated_v, expected_v, rel_tol=1e-15)
    # Physical sanity check: LEO orbital speed is approximately ~7.67 km/s
    assert 7.6 < calculated_v < 7.8


def test_create_circular_orbit_state_400km():
    """Verify initial state vector for 400 km circular orbit."""
    state = create_circular_orbit_state(altitude_km=400.0)
    expected_r = R_EARTH + 400.0

    # Position
    assert math.isclose(state.x, expected_r, rel_tol=1e-15)
    assert math.isclose(state.y, 0.0, abs_tol=1e-15)
    assert math.isclose(state.r, expected_r, rel_tol=1e-15)

    # Velocity
    assert math.isclose(state.vx, 0.0, abs_tol=1e-15)
    expected_vy = math.sqrt(MU_EARTH / expected_r)
    assert math.isclose(state.vy, expected_vy, rel_tol=1e-15)
    assert math.isclose(state.v, expected_vy, rel_tol=1e-15)


def test_perpendicularity_of_velocity_and_radius():
    """In a circular orbit at initial point, velocity must be strictly perpendicular to radius.

    Invariant: r0 . v0 = x0 * vx0 + y0 * vy0 == 0.
    """
    state = create_circular_orbit_state(altitude_km=400.0)
    dot_product = state.x * state.vx + state.y * state.vy
    assert math.isclose(dot_product, 0.0, abs_tol=1e-12)


def test_centripetal_acceleration_balance():
    """For a circular orbit, centripetal acceleration v^2 / r must exactly equal gravitational acceleration mu / r^2."""
    state = create_circular_orbit_state(altitude_km=400.0)
    ax, ay = gravitational_acceleration(state.x, state.y)
    a_grav_mag = math.hypot(ax, ay)

    a_centripetal = (state.v ** 2) / state.r
    assert math.isclose(a_grav_mag, a_centripetal, rel_tol=1e-14)

