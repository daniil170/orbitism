"""Velocity Verlet numerical integrator.

Decoupled from physical domain: operates generically on state vectors and
arbitrary derivative vector functions dstate/dt = f(t, state).
"""

import math
from typing import Callable, Sequence, Union

from src.models.state import State2D, StateDerivative2D

DerivativeCallable = Callable[[float, State2D], StateDerivative2D]
GenericDerivativeCallable = Callable[[float, Sequence[float]], Sequence[float]]


def velocity_verlet_step(
    t: float,
    state: Union[State2D, Sequence[float]],
    derivative_fn: Union[DerivativeCallable, GenericDerivativeCallable],
    dt: float,
) -> Union[State2D, list[float]]:
    """Perform a single step of the Velocity Verlet method.

    Formulas:
        r_{n+1} = r_n + v_n * dt + 0.5 * a(r_n) * dt^2
        a_{n+1} = a(r_{n+1})
        v_{n+1} = v_n + 0.5 * (a_n + a_{n+1}) * dt

    For 2D kinematic state:
        x_{n+1}  = x_n  + vx_n * dt + 0.5 * ax_n * dt^2
        y_{n+1}  = y_n  + vy_n * dt + 0.5 * ay_n * dt^2
        ax_{n+1}, ay_{n+1} = a(x_{n+1}, y_{n+1})
        vx_{n+1} = vx_n + 0.5 * (ax_n + ax_{n+1}) * dt
        vy_{n+1} = vy_n + 0.5 * (ay_n + ay_{n+1}) * dt

    Args:
        t: Current time in seconds.
        state: Current state (State2D instance or sequence of numbers).
        derivative_fn: Function returning dstate/dt given (t, state).
        dt: Timestep (delta t) in seconds.

    Returns:
        New state after step dt, in same type structure as input.
    """
    if isinstance(state, State2D):
        # Initial derivative evaluation at (t, state)
        d1 = derivative_fn(t, state)
        ax_n = d1.ax if isinstance(d1, StateDerivative2D) else d1[2]
        ay_n = d1.ay if isinstance(d1, StateDerivative2D) else d1[3]

        # Position update
        dt_sq_half = 0.5 * (dt**2)
        x_next = state.x + state.vx * dt + ax_n * dt_sq_half
        y_next = state.y + state.vy * dt + ay_n * dt_sq_half

        # Intermediate velocity for derivative call at t + dt
        vx_half = state.vx + 0.5 * ax_n * dt
        vy_half = state.vy + 0.5 * ay_n * dt
        state_for_accel = State2D(x=x_next, y=y_next, vx=vx_half, vy=vy_half)

        # Acceleration update at new position
        d2 = derivative_fn(t + dt, state_for_accel)
        ax_next = d2.ax if isinstance(d2, StateDerivative2D) else d2[2]
        ay_next = d2.ay if isinstance(d2, StateDerivative2D) else d2[3]

        # Velocity update
        dt_half = 0.5 * dt
        vx_next = state.vx + (ax_n + ax_next) * dt_half
        vy_next = state.vy + (ay_n + ay_next) * dt_half

        return State2D(x=x_next, y=y_next, vx=vx_next, vy=vy_next)

    # Generic sequence fallback: assumes first half is positions q, second half is velocities v
    s_vals = [float(v) for v in state]
    dim = len(s_vals)
    k = dim // 2

    q_n = s_vals[:k]
    v_n = s_vals[k:]

    d1 = derivative_fn(t, s_vals)
    d1_vals = (
        d1.to_list() if isinstance(d1, StateDerivative2D) else [float(x) for x in d1]
    )
    a_n = d1_vals[k:]

    dt_sq_half = 0.5 * (dt**2)
    q_next = [q_n[i] + v_n[i] * dt + a_n[i] * dt_sq_half for i in range(k)]
    v_half = [v_n[i] + 0.5 * a_n[i] * dt for i in range(k)]

    intermediate_seq = q_next + v_half
    d2 = derivative_fn(t + dt, intermediate_seq)
    d2_vals = (
        d2.to_list() if isinstance(d2, StateDerivative2D) else [float(x) for x in d2]
    )
    a_next = d2_vals[k:]

    dt_half = 0.5 * dt
    v_next = [v_n[i] + (a_n[i] + a_next[i]) * dt_half for i in range(k)]

    return q_next + v_next


class VelocityVerletIntegrator:
    """Velocity Verlet numerical integrator class with standard interface."""

    def __init__(self, derivative_fn: Union[DerivativeCallable, GenericDerivativeCallable]):
        self.derivative_fn = derivative_fn

    def step(
        self,
        t: float,
        state: Union[State2D, Sequence[float]],
        dt: float,
    ) -> Union[State2D, list[float]]:
        """Advance state by timestep dt."""
        return velocity_verlet_step(
            t=t, state=state, derivative_fn=self.derivative_fn, dt=dt
        )
