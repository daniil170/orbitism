"""Error and convergence metrics for numerical simulation benchmarking."""

import math
from src.models.state import State2D


def position_error(state_num: State2D, state_exact: State2D) -> float:
    """Compute Euclidean distance error between numerical and exact position.

    Formula:
        e_r = sqrt((x_num - x_exact)^2 + (y_num - y_exact)^2)

    Args:
        state_num: Numerical State2D.
        state_exact: Reference analytical State2D.

    Returns:
        Position error in kilometers (km).
    """
    return math.hypot(state_num.x - state_exact.x, state_num.y - state_exact.y)


def velocity_error(state_num: State2D, state_exact: State2D) -> float:
    """Compute Euclidean magnitude error between numerical and exact velocity.

    Formula:
        e_v = sqrt((vx_num - vx_exact)^2 + (vy_num - vy_exact)^2)

    Args:
        state_num: Numerical State2D.
        state_exact: Reference analytical State2D.

    Returns:
        Velocity error in kilometers per second (km/s).
    """
    return math.hypot(
        state_num.vx - state_exact.vx, state_num.vy - state_exact.vy
    )


def convergence_order(
    error_1: float, error_2: float, dt_1: float, dt_2: float
) -> float:
    """Compute empirical convergence order p between adjacent timesteps.

    Formula:
        p = log(error_1 / error_2) / log(dt_1 / dt_2)

    Args:
        error_1: Error with timestep dt_1 (dt_1 > dt_2).
        error_2: Error with timestep dt_2.
        dt_1: Coarser timestep in seconds.
        dt_2: Finer timestep in seconds.

    Returns:
        Empirical convergence order p.

    Raises:
        ValueError: If errors or timesteps are non-positive, or dt_1 == dt_2.
    """
    if error_1 <= 0.0 or error_2 <= 0.0:
        raise ValueError(
            f"Errors must be strictly positive, got error_1={error_1}, error_2={error_2}."
        )
    if dt_1 <= 0.0 or dt_2 <= 0.0:
        raise ValueError(
            f"Timesteps must be strictly positive, got dt_1={dt_1}, dt_2={dt_2}."
        )
    if dt_1 == dt_2:
        raise ValueError("Timesteps dt_1 and dt_2 must be different.")

    return math.log(error_1 / error_2) / math.log(dt_1 / dt_2)
