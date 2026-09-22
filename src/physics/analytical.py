"""Analytical solutions and Keplerian invariants for 2D orbital motion."""

import math
from src.models.state import State2D
from src.physics.gravity import MU_EARTH


def circular_orbit_mean_motion(r0: float, mu: float = MU_EARTH) -> float:
    """Compute mean motion n = sqrt(mu / r0^3) for a circular orbit.

    Args:
        r0: Orbital radius in km.
        mu: Gravitational parameter in km^3/s^2.

    Returns:
        Mean motion in rad/s.
    """
    if r0 <= 0.0:
        raise ValueError(f"Orbital radius must be positive, got {r0} km.")
    return math.sqrt(mu / (r0**3))


def circular_orbit_period(r0: float, mu: float = MU_EARTH) -> float:
    """Compute orbital period T = 2 * pi / n for a circular orbit.

    Args:
        r0: Orbital radius in km.
        mu: Gravitational parameter in km^3/s^2.

    Returns:
        Orbital period in seconds.
    """
    n = circular_orbit_mean_motion(r0, mu=mu)
    return 2.0 * math.pi / n


def circular_orbit_exact_state(
    t: float, r0: float, mu: float = MU_EARTH
) -> State2D:
    """Evaluate exact analytical state on a 2D circular orbit at time t.

    Formulas:
        n = sqrt(mu / r0^3)
        x(t)  =  r0 * cos(n * t)
        y(t)  =  r0 * sin(n * t)
        vx(t) = -r0 * n * sin(n * t)
        vy(t) =  r0 * n * cos(n * t)

    Args:
        t: Time in seconds.
        r0: Orbital radius in km.
        mu: Gravitational parameter in km^3/s^2.

    Returns:
        Exact State2D at time t.
    """
    n = circular_orbit_mean_motion(r0, mu=mu)
    nt = n * t
    x = r0 * math.cos(nt)
    y = r0 * math.sin(nt)
    vx = -r0 * n * math.sin(nt)
    vy = r0 * n * math.cos(nt)
    return State2D(x=x, y=y, vx=vx, vy=vy)


def specific_orbital_energy(state: State2D, mu: float = MU_EARTH) -> float:
    """Compute specific orbital energy epsilon = v^2 / 2 - mu / r.

    Args:
        state: Kinematic state.
        mu: Gravitational parameter in km^3/s^2.

    Returns:
        Specific orbital energy in km^2/s^2.
    """
    if state.r <= 0.0:
        raise ValueError("Orbital radius must be positive.")
    return 0.5 * (state.v**2) - (mu / state.r)


def specific_angular_momentum(state: State2D) -> float:
    """Compute 2D specific angular momentum h_z = x * vy - y * vx.

    Args:
        state: Kinematic state.

    Returns:
        Z-component of specific angular momentum in km^2/s.
    """
    return state.x * state.vy - state.y * state.vx
