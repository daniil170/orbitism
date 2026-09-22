"""Gravitational physics model for classical two-body problem."""

import math
from typing import Tuple, Union

from src.models.state import State2D, StateDerivative2D

# Physical constants for Earth
MU_EARTH: float = 398600.435  # Standard gravitational parameter, km^3/s^2
R_EARTH: float = 6378.137    # Equatorial radius of Earth, km


def gravitational_acceleration(
    x: float, y: float, mu: float = MU_EARTH
) -> Tuple[float, float]:
    """Compute 2D gravitational acceleration towards central body at origin.

    Formula:
        r = sqrt(x^2 + y^2)
        ax = -mu * x / r^3
        ay = -mu * y / r^3

    Args:
        x: Position x-coordinate in km.
        y: Position y-coordinate in km.
        mu: Gravitational parameter in km^3/s^2 (defaults to Earth's mu).

    Returns:
        Tuple (ax, ay) of acceleration components in km/s^2.

    Raises:
        ValueError: If radius is zero (singularity at coordinate origin).
    """
    r = math.hypot(x, y)
    if r == 0.0:
        raise ValueError(
            "Gravitational singularity: radius cannot be zero at coordinate origin."
        )

    r3 = r * r * r
    ax = -mu * x / r3
    ay = -mu * y / r3
    return (ax, ay)


def two_body_derivatives(
    t: float, state: Union[State2D, Tuple[float, float, float, float]], mu: float = MU_EARTH
) -> StateDerivative2D:
    """Compute state derivative dstate/dt = [vx, vy, ax, ay] for two-body problem.

    Args:
        t: Time parameter (seconds). Not explicitly used in autonomous gravity field.
        state: Current kinematic state (x, y, vx, vy) in km and km/s.
        mu: Gravitational parameter in km^3/s^2.

    Returns:
        StateDerivative2D instance with (vx, vy, ax, ay).
    """
    if isinstance(state, State2D):
        x, y, vx, vy = state.x, state.y, state.vx, state.vy
    else:
        x, y, vx, vy = state[0], state[1], state[2], state[3]

    ax, ay = gravitational_acceleration(x, y, mu=mu)
    return StateDerivative2D(vx=vx, vy=vy, ax=ax, ay=ay)


def circular_orbit_velocity(r: float, mu: float = MU_EARTH) -> float:
    """Compute the circular orbital velocity at distance r from central body.

    Formula:
        v = sqrt(mu / r)

    Args:
        r: Distance from center of mass in km.
        mu: Gravitational parameter in km^3/s^2.

    Returns:
        Circular orbital speed in km/s.

    Raises:
        ValueError: If radius r <= 0.
    """
    if r <= 0.0:
        raise ValueError(f"Radius must be strictly positive, got {r} km.")
    return math.sqrt(mu / r)


def create_circular_orbit_state(
    altitude_km: float = 400.0,
    mu: float = MU_EARTH,
    r_earth: float = R_EARTH,
) -> State2D:
    """Create initial state for a circular orbit at a specified altitude.

    Coordinate conventions:
        Satellite is placed on the +X axis: x0 = r0, y0 = 0.
        Velocity is perpendicular along +Y axis: vx0 = 0, vy0 = sqrt(mu / r0).

    Args:
        altitude_km: Altitude above Earth equatorial surface in km (default: 400 km).
        mu: Gravitational parameter in km^3/s^2.
        r_earth: Equatorial radius of Earth in km.

    Returns:
        State2D initial state vector [x0, y0, vx0, vy0].
    """
    r0 = r_earth + altitude_km
    vy0 = circular_orbit_velocity(r0, mu=mu)
    return State2D(x=r0, y=0.0, vx=0.0, vy=vy0)

