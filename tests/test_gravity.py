"""Unit tests for gravitational physics model."""

import math
import pytest

from src.models.state import State2D
from src.physics.gravity import (
    MU_EARTH,
    gravitational_acceleration,
    two_body_derivatives,
)


def test_singularity_at_coordinate_origin_raises_error():
    """In the coordinate origin (r=0), the model must not silently produce NaN/Inf."""
    with pytest.raises(ValueError, match="singularity|radius cannot be zero"):
        gravitational_acceleration(0.0, 0.0)


def test_acceleration_direction_along_x_axis():
    """For position [x, 0] with x > 0, acceleration must point in negative x direction."""
    x = 7000.0
    ax, ay = gravitational_acceleration(x, 0.0)
    assert ax < 0.0, f"Expected ax < 0 for positive x, got {ax}"
    assert math.isclose(ay, 0.0, abs_tol=1e-12), f"Expected ay == 0, got {ay}"

    # Verify symmetry for negative x
    ax_neg, ay_neg = gravitational_acceleration(-x, 0.0)
    assert ax_neg > 0.0, f"Expected ax > 0 for negative x, got {ax_neg}"
    assert math.isclose(ay_neg, 0.0, abs_tol=1e-12)
    assert math.isclose(ax, -ax_neg, rel_tol=1e-14)


def test_acceleration_direction_along_y_axis():
    """For position [0, y] with y > 0, acceleration must point in negative y direction."""
    y = 7000.0
    ax, ay = gravitational_acceleration(0.0, y)
    assert ay < 0.0, f"Expected ay < 0 for positive y, got {ay}"
    assert math.isclose(ax, 0.0, abs_tol=1e-12), f"Expected ax == 0, got {ax}"

    # Verify symmetry for negative y
    ax_neg, ay_neg = gravitational_acceleration(0.0, -y)
    assert ay_neg > 0.0, f"Expected ay > 0 for negative y, got {ay_neg}"
    assert math.isclose(ax_neg, 0.0, abs_tol=1e-12)
    assert math.isclose(ay, -ay_neg, rel_tol=1e-14)


def test_central_force_collinearity_invariant():
    """Gravitational acceleration must be central (antiparallel to radius vector).

    Invariant: The 2D cross product r x a = x * ay - y * ax must equal 0.
    """
    test_positions = [
        (3000.0, 4000.0),
        (-5000.0, 2000.0),
        (-4000.0, -4000.0),
        (1000.0, -6000.0),
    ]
    for x, y in test_positions:
        ax, ay = gravitational_acceleration(x, y)
        cross_product = x * ay - y * ax
        assert math.isclose(cross_product, 0.0, abs_tol=1e-12), (
            f"Cross product r x a should be 0, got {cross_product} for ({x}, {y})"
        )
        # Dot product r . a must be strictly negative (attractive force)
        dot_product = x * ax + y * ay
        assert dot_product < 0.0, f"Dot product r . a should be negative, got {dot_product}"


def test_inverse_square_law_invariant():
    """Doubling the distance must reduce acceleration magnitude by exactly a factor of 4."""
    r1 = 7000.0
    r2 = 2.0 * r1

    ax1, _ = gravitational_acceleration(r1, 0.0)
    ax2, _ = gravitational_acceleration(r2, 0.0)

    ratio = abs(ax1) / abs(ax2)
    assert math.isclose(ratio, 4.0, rel_tol=1e-12), f"Expected ratio 4.0, got {ratio}"


def test_two_body_derivatives_consistency():
    """Verify that two_body_derivatives correctly packs [vx, vy, ax, ay]."""
    state = State2D(x=6778.137, y=0.0, vx=0.0, vy=7.67)
    deriv = two_body_derivatives(0.0, state, mu=MU_EARTH)

    assert deriv.vx == state.vx
    assert deriv.vy == state.vy
    expected_ax, expected_ay = gravitational_acceleration(state.x, state.y, mu=MU_EARTH)
    assert math.isclose(deriv.ax, expected_ax, rel_tol=1e-14)
    assert math.isclose(deriv.ay, expected_ay, rel_tol=1e-14)

