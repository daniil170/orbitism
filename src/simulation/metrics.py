"""Error and convergence metrics for numerical simulation benchmarking."""

from dataclasses import dataclass
import math
from typing import Iterator, List, Sequence, Tuple

from src.models.state import State2D
from src.physics.gravity import MU_EARTH
from src.physics.analytical import (
    circular_orbit_exact_state,
    specific_orbital_energy,
    specific_angular_momentum,
)


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


@dataclass(frozen=True)
class ErrorRecord:
    """Record of kinematic state, errors, and physical invariants at time t.

    Attributes:
        t: Simulation timestamp in seconds.
        state: Kinematic state State2D at timestamp t.
        position_error: Euclidean position error e_r in km.
        velocity_error: Euclidean velocity error e_v in km/s.
        energy_error: Relative specific orbital energy error (epsilon - epsilon0) / |epsilon0|.
        angular_momentum_error: Relative specific angular momentum error (h - h0) / |h0|.
    """

    t: float
    state: State2D
    position_error: float
    velocity_error: float
    energy_error: float
    angular_momentum_error: float

    def __iter__(self) -> Iterator:
        """Allow tuple unpacking (t, state, pos_err, vel_err, energy_err, ang_mom_err)."""
        return iter(
            (
                self.t,
                self.state,
                self.position_error,
                self.velocity_error,
                self.energy_error,
                self.angular_momentum_error,
            )
        )

    def to_tuple(self) -> Tuple[float, State2D, float, float, float, float]:
        """Convert record to a tuple."""
        return (
            self.t,
            self.state,
            self.position_error,
            self.velocity_error,
            self.energy_error,
            self.angular_momentum_error,
        )


def relative_energy_error(
    state: State2D, ref_energy: float, mu: float = MU_EARTH
) -> float:
    """Compute relative error in specific orbital energy: (epsilon - epsilon0) / |epsilon0|.

    Args:
        state: Current kinematic state.
        ref_energy: Reference specific orbital energy epsilon0 in km^2/s^2.
        mu: Gravitational parameter in km^3/s^2.

    Returns:
        Relative energy error (dimensionless).
    """
    if abs(ref_energy) == 0.0:
        raise ValueError("Reference energy magnitude cannot be zero.")
    energy = specific_orbital_energy(state, mu=mu)
    return (energy - ref_energy) / abs(ref_energy)


def relative_angular_momentum_error(
    state: State2D, ref_angular_momentum: float
) -> float:
    """Compute relative error in specific angular momentum: (h_z - h0) / |h0|.

    Args:
        state: Current kinematic state.
        ref_angular_momentum: Reference specific angular momentum h0 in km^2/s.

    Returns:
        Relative angular momentum error (dimensionless).
    """
    if abs(ref_angular_momentum) == 0.0:
        raise ValueError("Reference angular momentum magnitude cannot be zero.")
    h = specific_angular_momentum(state)
    return (h - ref_angular_momentum) / abs(ref_angular_momentum)


def compute_error_history(
    times: Sequence[float],
    trajectory: Sequence[State2D],
    r0: float,
    mu: float = MU_EARTH,
) -> List[ErrorRecord]:
    """Compute full time evolution of errors and diagnostic invariants along trajectory.

    For each timestamp t_k and state s_k, evaluates the exact reference circular
    orbit state at time t_k and calculates position error, velocity error,
    energy error, and angular momentum error.

    Args:
        times: Sequence of timestamps [t_0, t_1, ..., t_N] in seconds.
        trajectory: Sequence of numerical states [s_0, s_1, ..., s_N].
        r0: Reference orbital radius in km.
        mu: Gravitational parameter in km^3/s^2.

    Returns:
        List of ErrorRecord instances from t_0 to t_N.

    Raises:
        ValueError: If lengths of times and trajectory do not match or are empty.
    """
    if len(times) != len(trajectory):
        raise ValueError(
            f"Length mismatch: {len(times)} timestamps vs {len(trajectory)} states."
        )
    if not times:
        raise ValueError("Trajectory history cannot be empty.")

    ref_energy = -mu / (2.0 * r0)
    ref_angular_momentum = r0 * math.sqrt(mu / r0)

    history: List[ErrorRecord] = []
    for t_k, state_k in zip(times, trajectory):
        exact_k = circular_orbit_exact_state(t_k, r0=r0, mu=mu)
        e_r = position_error(state_k, exact_k)
        e_v = velocity_error(state_k, exact_k)
        e_energy = relative_energy_error(state_k, ref_energy=ref_energy, mu=mu)
        e_ang = relative_angular_momentum_error(
            state_k, ref_angular_momentum=ref_angular_momentum
        )

        history.append(
            ErrorRecord(
                t=t_k,
                state=state_k,
                position_error=e_r,
                velocity_error=e_v,
                energy_error=e_energy,
                angular_momentum_error=e_ang,
            )
        )

    return history
