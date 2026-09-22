"""Classical 4th-order Runge-Kutta (RK4) numerical integrator.

Decoupled from physical domain: operates generically on state vectors and
arbitrary derivative vector functions dstate/dt = f(t, state).
"""

from typing import Callable, Sequence, Union

from src.models.state import State2D, StateDerivative2D

DerivativeCallable = Callable[[float, State2D], StateDerivative2D]
GenericDerivativeCallable = Callable[[float, Sequence[float]], Sequence[float]]


def rk4_step(
    t: float,
    state: Union[State2D, Sequence[float]],
    derivative_fn: Union[DerivativeCallable, GenericDerivativeCallable],
    dt: float,
) -> Union[State2D, list[float]]:
    """Perform a single step of the classical 4th-order Runge-Kutta method (RK4).

    Formulas:
        k1 = f(t, s)
        k2 = f(t + dt / 2, s + dt / 2 * k1)
        k3 = f(t + dt / 2, s + dt / 2 * k2)
        k4 = f(t + dt, s + dt * k3)
        s_{next} = s + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)

    Args:
        t: Current time in seconds.
        state: Current state (State2D instance or sequence of numbers).
        derivative_fn: Function returning dstate/dt given (t, state).
        dt: Timestep (delta t) in seconds.

    Returns:
        New state after step dt, in same type structure as input.
    """
    if isinstance(state, State2D):
        # Stage 1: evaluate derivative at initial state
        d1 = derivative_fn(t, state)
        k1_vx = d1.vx if isinstance(d1, StateDerivative2D) else d1[0]
        k1_vy = d1.vy if isinstance(d1, StateDerivative2D) else d1[1]
        k1_ax = d1.ax if isinstance(d1, StateDerivative2D) else d1[2]
        k1_ay = d1.ay if isinstance(d1, StateDerivative2D) else d1[3]

        # Stage 2: evaluate derivative at midpoint using k1
        t_mid = t + 0.5 * dt
        s2 = State2D(
            x=state.x + 0.5 * dt * k1_vx,
            y=state.y + 0.5 * dt * k1_vy,
            vx=state.vx + 0.5 * dt * k1_ax,
            vy=state.vy + 0.5 * dt * k1_ay,
        )
        d2 = derivative_fn(t_mid, s2)
        k2_vx = d2.vx if isinstance(d2, StateDerivative2D) else d2[0]
        k2_vy = d2.vy if isinstance(d2, StateDerivative2D) else d2[1]
        k2_ax = d2.ax if isinstance(d2, StateDerivative2D) else d2[2]
        k2_ay = d2.ay if isinstance(d2, StateDerivative2D) else d2[3]

        # Stage 3: evaluate derivative at midpoint using k2
        s3 = State2D(
            x=state.x + 0.5 * dt * k2_vx,
            y=state.y + 0.5 * dt * k2_vy,
            vx=state.vx + 0.5 * dt * k2_ax,
            vy=state.vy + 0.5 * dt * k2_ay,
        )
        d3 = derivative_fn(t_mid, s3)
        k3_vx = d3.vx if isinstance(d3, StateDerivative2D) else d3[0]
        k3_vy = d3.vy if isinstance(d3, StateDerivative2D) else d3[1]
        k3_ax = d3.ax if isinstance(d3, StateDerivative2D) else d3[2]
        k3_ay = d3.ay if isinstance(d3, StateDerivative2D) else d3[3]

        # Stage 4: evaluate derivative at endpoint using k3
        t_end = t + dt
        s4 = State2D(
            x=state.x + dt * k3_vx,
            y=state.y + dt * k3_vy,
            vx=state.vx + dt * k3_ax,
            vy=state.vy + dt * k3_ay,
        )
        d4 = derivative_fn(t_end, s4)
        k4_vx = d4.vx if isinstance(d4, StateDerivative2D) else d4[0]
        k4_vy = d4.vy if isinstance(d4, StateDerivative2D) else d4[1]
        k4_ax = d4.ax if isinstance(d4, StateDerivative2D) else d4[2]
        k4_ay = d4.ay if isinstance(d4, StateDerivative2D) else d4[3]

        # Combine stages using Simpson's rule weighting
        dt_sixth = dt / 6.0
        return State2D(
            x=state.x + dt_sixth * (k1_vx + 2.0 * k2_vx + 2.0 * k3_vx + k4_vx),
            y=state.y + dt_sixth * (k1_vy + 2.0 * k2_vy + 2.0 * k3_vy + k4_vy),
            vx=state.vx + dt_sixth * (k1_ax + 2.0 * k2_ax + 2.0 * k3_ax + k4_ax),
            vy=state.vy + dt_sixth * (k1_ay + 2.0 * k2_ay + 2.0 * k3_ay + k4_ay),
        )

    # Generic sequence fallback
    s_vals = [float(v) for v in state]

    def _get_deriv_list(d: Union[StateDerivative2D, Sequence[float]]) -> list[float]:
        if isinstance(d, StateDerivative2D):
            return d.to_list()
        return [float(v) for v in d]

    d1 = derivative_fn(t, s_vals)
    k1 = _get_deriv_list(d1)

    t_mid = t + 0.5 * dt
    s2 = [s + 0.5 * dt * k for s, k in zip(s_vals, k1)]
    d2 = derivative_fn(t_mid, s2)
    k2 = _get_deriv_list(d2)

    s3 = [s + 0.5 * dt * k for s, k in zip(s_vals, k2)]
    d3 = derivative_fn(t_mid, s3)
    k3 = _get_deriv_list(d3)

    t_end = t + dt
    s4 = [s + dt * k for s, k in zip(s_vals, k3)]
    d4 = derivative_fn(t_end, s4)
    k4 = _get_deriv_list(d4)

    dt_sixth = dt / 6.0
    return [
        s + dt_sixth * (k1_i + 2.0 * k2_i + 2.0 * k3_i + k4_i)
        for s, k1_i, k2_i, k3_i, k4_i in zip(s_vals, k1, k2, k3, k4)
    ]


class RK4Integrator:
    """Classical 4th-order Runge-Kutta (RK4) integrator class with standard interface."""

    def __init__(self, derivative_fn: Union[DerivativeCallable, GenericDerivativeCallable]):
        self.derivative_fn = derivative_fn

    def step(
        self,
        t: float,
        state: Union[State2D, Sequence[float]],
        dt: float,
    ) -> Union[State2D, list[float]]:
        """Advance state by timestep dt using classical RK4."""
        return rk4_step(t=t, state=state, derivative_fn=self.derivative_fn, dt=dt)
