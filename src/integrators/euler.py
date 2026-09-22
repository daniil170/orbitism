"""Explicit Euler numerical integrator.

Decoupled from physical domain: operates generically on state vectors and
arbitrary derivative vector functions dstate/dt = f(t, state).
"""

from typing import Callable, Sequence, Union

from src.models.state import State2D, StateDerivative2D

DerivativeCallable = Callable[[float, State2D], StateDerivative2D]
GenericDerivativeCallable = Callable[[float, Sequence[float]], Sequence[float]]


def euler_step(
    t: float,
    state: Union[State2D, Sequence[float]],
    derivative_fn: Union[DerivativeCallable, GenericDerivativeCallable],
    dt: float,
) -> Union[State2D, list[float]]:
    """Perform a single step of the Explicit Euler method.

    Formula:
        state_{n+1} = state_n + f(t_n, state_n) * dt

    For 2D kinematic state:
        x_{n+1}  = x_n  + vx_n * dt
        y_{n+1}  = y_n  + vy_n * dt
        vx_{n+1} = vx_n + ax_n * dt
        vy_{n+1} = vy_n + ay_n * dt

    Args:
        t: Current time.
        state: Current state (State2D instance or sequence of numbers).
        derivative_fn: Function returning dstate/dt given (t, state).
        dt: Timestep (delta t).

    Returns:
        New state after step dt, in same type structure as input.
    """
    deriv = derivative_fn(t, state)

    if isinstance(state, State2D):
        if isinstance(deriv, StateDerivative2D):
            vx, vy, ax, ay = deriv.vx, deriv.vy, deriv.ax, deriv.ay
        else:
            vx, vy, ax, ay = deriv[0], deriv[1], deriv[2], deriv[3]

        return State2D(
            x=state.x + vx * dt,
            y=state.y + vy * dt,
            vx=state.vx + ax * dt,
            vy=state.vy + ay * dt,
        )

    # Generic sequence fallback
    s_vals = list(state)
    d_vals = list(deriv) if not isinstance(deriv, StateDerivative2D) else deriv.to_list()
    return [s + d * dt for s, d in zip(s_vals, d_vals)]


class EulerIntegrator:
    """Explicit Euler integrator class with standard interface."""

    def __init__(self, derivative_fn: Union[DerivativeCallable, GenericDerivativeCallable]):
        self.derivative_fn = derivative_fn

    def step(
        self,
        t: float,
        state: Union[State2D, Sequence[float]],
        dt: float,
    ) -> Union[State2D, list[float]]:
        """Advance state by timestep dt."""
        return euler_step(t=t, state=state, derivative_fn=self.derivative_fn, dt=dt)

