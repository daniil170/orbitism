"""Error and convergence metrics for numerical simulation benchmarking."""

from dataclasses import dataclass
import math
from typing import Iterator, List, Sequence, Tuple

from src.models.state import State2D
from src.physics.gravity import MU_EARTH
from src.physics.analytical import (
    circular_orbit_exact_state,
    circular_orbit_mean_motion,
    circular_orbit_period,
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


def radial_position_error(
    state_num: State2D,
    state_exact: State2D,
    t: float,
    r0: float,
    mu: float = MU_EARTH,
) -> float:
    """Compute radial component of position error in rotating frame: e_R.

    Formula:
        e_R(t) = Delta r(t) . e_r(t) = (x_num - x_exact) * cos(nt) + (y_num - y_exact) * sin(nt)

    Args:
        state_num: Numerical State2D.
        state_exact: Reference analytical State2D.
        t: Simulation time in seconds.
        r0: Reference orbital radius in km.
        mu: Gravitational parameter in km^3/s^2.

    Returns:
        Radial position error in km.
    """
    n = circular_orbit_mean_motion(r0, mu=mu)
    nt = n * t
    dx = state_num.x - state_exact.x
    dy = state_num.y - state_exact.y
    return dx * math.cos(nt) + dy * math.sin(nt)


def along_track_position_error(
    state_num: State2D,
    state_exact: State2D,
    t: float,
    r0: float,
    mu: float = MU_EARTH,
) -> float:
    """Compute along-track (transverse) component of position error in rotating frame: e_T.

    Formula:
        e_T(t) = Delta r(t) . e_theta(t) = -(x_num - x_exact) * sin(nt) + (y_num - y_exact) * cos(nt)

    Args:
        state_num: Numerical State2D.
        state_exact: Reference analytical State2D.
        t: Simulation time in seconds.
        r0: Reference orbital radius in km.
        mu: Gravitational parameter in km^3/s^2.

    Returns:
        Along-track position error in km.
    """
    n = circular_orbit_mean_motion(r0, mu=mu)
    nt = n * t
    dx = state_num.x - state_exact.x
    dy = state_num.y - state_exact.y
    return -dx * math.sin(nt) + dy * math.cos(nt)


def radial_drift(state_num: State2D, r0: float) -> float:
    """Compute radial drift Delta r = r_num - r0.

    Args:
        state_num: Numerical State2D.
        r0: Reference circular orbital radius in km.

    Returns:
        Radial drift in km.
    """
    return state_num.r - r0


def unwrap_phase(angles: Sequence[float]) -> List[float]:
    """Unwrap 1D sequence of radian angles across [-pi, pi] branch cuts.

    Adjusts cumulative phase so that absolute jumps between consecutive points
    satisfy |delta_theta| <= pi.

    Args:
        angles: Sequence of raw radian angles (e.g. from atan2).

    Returns:
        List of continuous unwrapped radian angles.
    """
    if not angles:
        return []
    unwrapped = [angles[0]]
    for i in range(1, len(angles)):
        d_theta = angles[i] - angles[i - 1]
        d_theta = (d_theta + math.pi) % (2.0 * math.pi) - math.pi
        unwrapped.append(unwrapped[-1] + d_theta)
    return unwrapped


def angular_phase_error(
    theta_unwrapped: float,
    t: float,
    r0: float,
    mu: float = MU_EARTH,
) -> float:
    """Compute unwrapped angular phase error Delta theta = theta_unwrapped - n * t.

    Args:
        theta_unwrapped: Continuous unwrapped phase in radians.
        t: Simulation time in seconds.
        r0: Reference circular orbital radius in km.
        mu: Gravitational parameter in km^3/s^2.

    Returns:
        Angular phase error in radians.
    """
    n = circular_orbit_mean_motion(r0, mu=mu)
    return theta_unwrapped - n * t


def radial_error_fraction(e_R: float, e_r: float) -> float:
    """Compute fraction of position error variance along radial direction: rho_R.

    Formula:
        rho_R = e_R^2 / e_r^2 (0.0 if e_r <= 1e-15)

    Args:
        e_R: Radial position error in km.
        e_r: Total Euclidean position error in km.

    Returns:
        Dimensionless variance fraction in [0, 1].
    """
    if e_r <= 1e-15:
        return 0.0
    return (e_R**2) / (e_r**2)


def along_track_error_fraction(e_T: float, e_r: float) -> float:
    """Compute fraction of position error variance along along-track direction: rho_T.

    Formula:
        rho_T = e_T^2 / e_r^2 (0.0 if e_r <= 1e-15)

    Args:
        e_T: Along-track position error in km.
        e_r: Total Euclidean position error in km.

    Returns:
        Dimensionless variance fraction in [0, 1].
    """
    if e_r <= 1e-15:
        return 0.0
    return (e_T**2) / (e_r**2)


def dominance_ratio(e_R: float, e_T: float) -> float:
    """Compute ratio of along-track to radial error magnitude: chi = |e_T| / |e_R|.

    Args:
        e_R: Radial position error in km.
        e_T: Along-track position error in km.

    Returns:
        Dimensionless ratio chi. Returns inf if |e_R| == 0 and |e_T| > 0; 1.0 if both 0.
    """
    abs_r = abs(e_R)
    abs_t = abs(e_T)
    if abs_r <= 1e-15:
        return float("inf") if abs_t > 1e-15 else 1.0
    return abs_t / abs_r


@dataclass(frozen=True)
class DecomposedErrorRecord:
    """Record of kinematic state, decomposition metrics, and physical invariants at time t.

    Attributes:
        t: Simulation timestamp in seconds.
        state: Kinematic state State2D at timestamp t.
        position_error: Euclidean position error e_r in km.
        velocity_error: Euclidean velocity error e_v in km/s.
        energy_error: Relative specific orbital energy error.
        angular_momentum_error: Relative specific angular momentum error.
        r_num: Instantaneous numerical orbital radius in km.
        delta_r: Radial drift r_num - r0 in km.
        theta_unwrapped: Continuous unwrapped phase in radians.
        delta_theta: Unwrapped phase error theta_unwrapped - n * t in radians.
        e_R: Radial position error projection in km.
        e_T: Along-track position error projection in km.
        fraction_R: Radial error variance fraction rho_R in [0, 1].
        fraction_T: Along-track error variance fraction rho_T in [0, 1].
    """

    t: float
    state: State2D
    position_error: float
    velocity_error: float
    energy_error: float
    angular_momentum_error: float
    r_num: float
    delta_r: float
    theta_unwrapped: float
    delta_theta: float
    e_R: float
    e_T: float
    fraction_R: float
    fraction_T: float


def compute_decomposed_error_history(
    times: Sequence[float],
    trajectory: Sequence[State2D],
    r0: float,
    mu: float = MU_EARTH,
) -> List[DecomposedErrorRecord]:
    """Compute full time evolution of decomposed errors and invariants along trajectory.

    Args:
        times: Sequence of timestamps [t_0, t_1, ..., t_N] in seconds.
        trajectory: Sequence of numerical states [s_0, s_1, ..., s_N].
        r0: Reference orbital radius in km.
        mu: Gravitational parameter in km^3/s^2.

    Returns:
        List of DecomposedErrorRecord instances from t_0 to t_N.

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

    # First pass: compute raw atan2 angles
    raw_angles = [math.atan2(s.y, s.x) for s in trajectory]
    unwrapped_angles = unwrap_phase(raw_angles)

    history: List[DecomposedErrorRecord] = []
    for t_k, state_k, theta_u in zip(times, trajectory, unwrapped_angles):
        exact_k = circular_orbit_exact_state(t_k, r0=r0, mu=mu)
        e_r = position_error(state_k, exact_k)
        e_v = velocity_error(state_k, exact_k)
        e_energy = relative_energy_error(state_k, ref_energy=ref_energy, mu=mu)
        e_ang = relative_angular_momentum_error(
            state_k, ref_angular_momentum=ref_angular_momentum
        )

        r_num_k = state_k.r
        dr_k = radial_drift(state_k, r0)
        d_theta_k = angular_phase_error(theta_u, t_k, r0, mu=mu)
        e_R_k = radial_position_error(state_k, exact_k, t_k, r0, mu=mu)
        e_T_k = along_track_position_error(state_k, exact_k, t_k, r0, mu=mu)
        frac_R = radial_error_fraction(e_R_k, e_r)
        frac_T = along_track_error_fraction(e_T_k, e_r)

        history.append(
            DecomposedErrorRecord(
                t=t_k,
                state=state_k,
                position_error=e_r,
                velocity_error=e_v,
                energy_error=e_energy,
                angular_momentum_error=e_ang,
                r_num=r_num_k,
                delta_r=dr_k,
                theta_unwrapped=theta_u,
                delta_theta=d_theta_k,
                e_R=e_R_k,
                e_T=e_T_k,
                fraction_R=frac_R,
                fraction_T=frac_T,
            )
        )

    return history


def cubic_hermite_root(
    y0: float,
    y1: float,
    dy0: float,
    dy1: float,
    dt: float,
    tol: float = 1e-12,
    max_iter: int = 50,
) -> float:
    """Find root s* in [0, 1] of cubic Hermite polynomial P(s) = 0.

    Interpolates y(t) across a step [t_0, t_0 + dt] with normalized parameter
    s = (t - t_0) / dt in [0, 1], where y(0) = y0, y(1) = y1, y'(0) = dt * dy0,
    y'(1) = dt * dy1.

    Uses Newton-Raphson iterations with bracket-preserving bisection fallback.

    Args:
        y0: Value at s=0.
        y1: Value at s=1.
        dy0: Derivative dy/dt at s=0.
        dy1: Derivative dy/dt at s=1.
        dt: Interval step size in seconds (must be positive).
        tol: Convergence tolerance for |P(s)| and interval width.
        max_iter: Maximum number of iterations.

    Returns:
        Root s* in [0, 1].

    Raises:
        ValueError: If dt <= 0 or endpoints do not bracket zero.
        RuntimeError: If root finding fails to converge within max_iter.
    """
    if dt <= 0:
        raise ValueError(f"Timestep dt must be positive, got {dt}.")

    if abs(y0) < 1e-15:
        return 0.0
    if abs(y1) < 1e-15:
        return 1.0

    if y0 * y1 > 0:
        raise ValueError(
            f"Endpoints do not bracket zero: y0={y0}, y1={y1}."
        )

    def p(s: float) -> float:
        s2 = s * s
        s3 = s2 * s
        return (1.0 - 3.0 * s2 + 2.0 * s3) * y0 + s2 * (3.0 - 2.0 * s) * y1 + dt * (
            s * (1.0 - s) ** 2 * dy0 + s2 * (s - 1.0) * dy1
        )

    def dp(s: float) -> float:
        s2 = s * s
        return 6.0 * (s2 - s) * y0 + 6.0 * (s - s2) * y1 + dt * (
            (3.0 * s2 - 4.0 * s + 1.0) * dy0 + (3.0 * s2 - 2.0 * s) * dy1
        )

    a, b = 0.0, 1.0
    p_a = y0

    # Initial guess via linear interpolation
    s = -y0 / (y1 - y0)
    s = max(0.0, min(1.0, s))

    for _ in range(max_iter):
        val = p(s)
        if abs(val) < tol:
            return s

        if p_a * val <= 0.0:
            b = s
        else:
            a = s
            p_a = val

        if (b - a) < tol:
            return 0.5 * (a + b)

        d = dp(s)
        if abs(d) > 1e-15:
            s_new = s - val / d
            if a < s_new < b:
                s = s_new
                continue

        s = 0.5 * (a + b)

    final_val = p(s)
    if abs(final_val) <= 10.0 * tol or (b - a) <= 10.0 * tol:
        return s

    raise RuntimeError(
        f"Cubic Hermite root solver failed to converge: s={s}, P(s)={final_val}, "
        f"bracket=[{a}, {b}], tol={tol}, max_iter={max_iter}."
    )


@dataclass(frozen=True)
class OrbitCrossingEvent:
    """Record of an orbital crossing of the positive x-axis in the prograde direction.

    Event condition: y_i <= 0 < y_{i+1} with x > 0 and v_y > 0.

    Attributes:
        orbit_index: 1-based completed orbit count m = 1, 2, ...
        t_cross: Precise crossing timestamp in seconds.
        state_cross: State2D at crossing timestamp.
        period_instantaneous: Duration of this revolution: t_cross(m) - t_cross(m-1) [s].
        period_cumulative: Cumulative mean period: t_cross(m) / m [s].
        delta_T_cumulative: Relative period error: (period_cumulative - T_exact) / T_exact.
        step_index: Discrete step index preceding the crossing.
    """

    orbit_index: int
    t_cross: float
    state_cross: State2D
    period_instantaneous: float
    period_cumulative: float
    delta_T_cumulative: float
    step_index: int


def detect_positive_x_crossings(
    times: Sequence[float],
    trajectory: Sequence[State2D],
    r0: float,
    mu: float = MU_EARTH,
) -> List[OrbitCrossingEvent]:
    """Detect crossings of the positive x-axis in the prograde direction.

    Identifies events where y crosses 0 with positive velocity (y_i <= 0 < y_{i+1},
    x > 0, v_y > 0). Uses cubic Hermite interpolation to find sub-step crossing
    timestamps and states.

    Initial state at t=0 on positive x-axis (t*_0 = 0) is known analytically;
    the first detected crossing corresponds to m=1 (one complete revolution).

    Args:
        times: Sequence of discrete timestamps in seconds.
        trajectory: Sequence of numerical states State2D.
        r0: Reference circular orbital radius in km.
        mu: Gravitational parameter in km^3/s^2.

    Returns:
        List of OrbitCrossingEvent instances for m = 1, 2, ...

    Raises:
        ValueError: If input sequences have unequal length or fewer than 2 points.
    """
    if len(times) != len(trajectory):
        raise ValueError(
            f"Length mismatch: {len(times)} timestamps vs {len(trajectory)} states."
        )
    if len(times) < 2:
        raise ValueError("Trajectory must have at least 2 points for event detection.")

    t_exact = circular_orbit_period(r0, mu=mu)
    events: List[OrbitCrossingEvent] = []
    t_prev_cross = 0.0
    orbit_count = 0

    for i in range(len(trajectory) - 1):
        s0 = trajectory[i]
        s1 = trajectory[i + 1]
        t0 = times[i]
        t1 = times[i + 1]
        dt = t1 - t0

        if i == 0 and abs(s0.y) < 1e-14:
            continue

        if s0.y <= 0.0 < s1.y and s0.x > 0.0 and s1.x > 0.0 and s0.vy > 0.0 and s1.vy > 0.0:
            s_root = cubic_hermite_root(
                y0=s0.y,
                y1=s1.y,
                dy0=s0.vy,
                dy1=s1.vy,
                dt=dt,
            )
            t_cross = t0 + s_root * dt
            orbit_count += 1

            s2 = s_root * s_root
            s3 = s2 * s_root
            h00 = 1.0 - 3.0 * s2 + 2.0 * s3
            h10 = s_root * (1.0 - s_root) ** 2
            h01 = 3.0 * s2 - 2.0 * s3
            h11 = s2 * (s_root - 1.0)

            x_cross = h00 * s0.x + h01 * s1.x + dt * (h10 * s0.vx + h11 * s1.vx)
            y_cross = 0.0

            dh00 = 6.0 * (s2 - s_root)
            dh10 = 3.0 * s2 - 4.0 * s_root + 1.0
            dh01 = 6.0 * (s_root - s2)
            dh11 = 3.0 * s2 - 2.0 * s_root

            vx_cross = (dh00 * s0.x + dh01 * s1.x) / dt + dh10 * s0.vx + dh11 * s1.vx
            vy_cross = (dh00 * s0.y + dh01 * s1.y) / dt + dh10 * s0.vy + dh11 * s1.vy

            state_cross = State2D(x=x_cross, y=y_cross, vx=vx_cross, vy=vy_cross)

            t_inst = t_cross - t_prev_cross
            t_cum = t_cross / orbit_count
            delta_t_cum = (t_cum - t_exact) / t_exact

            events.append(
                OrbitCrossingEvent(
                    orbit_index=orbit_count,
                    t_cross=t_cross,
                    state_cross=state_cross,
                    period_instantaneous=t_inst,
                    period_cumulative=t_cum,
                    delta_T_cumulative=delta_t_cum,
                    step_index=i,
                )
            )
            t_prev_cross = t_cross

    return events


def phase_derived_chord_error(delta_theta: float, r0: float) -> float:
    """Compute geometric chord error between two points on circular orbit: 2*r0*|sin(Delta_theta / 2)|.

    This serves as a reference geometric chord model under zero radial error.

    Args:
        delta_theta: Angular phase difference in radians.
        r0: Reference orbital radius in km.

    Returns:
        Chord error in km.
    """
    return 2.0 * r0 * abs(math.sin(delta_theta / 2.0))


def phase_derived_along_track_displacement(delta_theta: float, r0: float) -> float:
    """Compute first-order along-track arc displacement diagnostic: r0 * Delta_theta.

    Args:
        delta_theta: Angular phase error in radians.
        r0: Reference orbital radius in km.

    Returns:
        Along-track displacement in km.
    """
    return r0 * delta_theta
