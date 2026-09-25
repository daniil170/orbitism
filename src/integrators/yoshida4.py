"""Yoshida 4th-order symplectic numerical integrator.

Constructed via symmetric composition of three 2nd-order symmetric symplectic steps
(Velocity Verlet / leapfrog):
    S4(h) = S2(w1 * h) ∘ S2(w0 * h) ∘ S2(w1 * h)

Coefficients origin (Yoshida, 1990):
    w1 = 1 / (2 - 2^(1/3))
    w0 = -2^(1/3) / (2 - 2^(1/3))

Properties:
    1) Consistency: 2 * w1 + w0 = 1
    2) Order condition: 2 * w1^3 + w0^3 = 0, eliminating the 3rd-order Lie error term
       and yielding global order 4 (error O(h^4)).
    3) Symplecticity: each substep preserves the canonical symplectic form dq ∧ dp,
       so their composition strictly preserves the symplectic 2-form.
    4) Time-reversibility: S4(-h) ∘ S4(h) = I.

Computational cost:
    Composed of 3 independent Velocity Verlet substeps.
    Each substep performs 2 evaluations of the derivative function (acceleration):
        Substep 1: at t and t + w1 * h (2 evals)
        Substep 2: at t + w1 * h and t + (w1 + w0) * h (2 evals)
        Substep 3: at t + (w1 + w0) * h and t + (2*w1 + w0) * h = t + h (2 evals)
    Total: exactly 6 force evaluations per outer timestep without caching.
"""

import math
from typing import Callable, Sequence, Union

from src.integrators.velocity_verlet import velocity_verlet_step
from src.models.state import State2D, StateDerivative2D

DerivativeCallable = Callable[[float, State2D], StateDerivative2D]
GenericDerivativeCallable = Callable[[float, Sequence[float]], Sequence[float]]

# Yoshida composition coefficients
_CBRT_2 = 2.0 ** (1.0 / 3.0)
W1: float = 1.0 / (2.0 - _CBRT_2)
W0: float = -_CBRT_2 / (2.0 - _CBRT_2)


def yoshida4_step(
    t: float,
    state: Union[State2D, Sequence[float]],
    derivative_fn: Union[DerivativeCallable, GenericDerivativeCallable],
    dt: float,
) -> Union[State2D, list[float]]:
    """Perform a single step of the Yoshida 4th-order symplectic integrator.

    Composition:
        h1 = w1 * dt
        h0 = w0 * dt
        s1 = S2(t, s, h1)
        s2 = S2(t + h1, s1, h0)
        s3 = S2(t + h1 + h0, s2, h1)

    Args:
        t: Current time in seconds.
        state: Current state (State2D instance or sequence of numbers).
        derivative_fn: Function returning dstate/dt given (t, state).
        dt: Timestep (delta t) in seconds.

    Returns:
        New state after step dt, in same type structure as input.
    """
    h1 = W1 * dt
    h0 = W0 * dt

    # Substep 1: S2(w1 * dt) starting at t
    s1 = velocity_verlet_step(t=t, state=state, derivative_fn=derivative_fn, dt=h1)

    # Substep 2: S2(w0 * dt) starting at t + h1 (negative substep)
    t2 = t + h1
    s2 = velocity_verlet_step(t=t2, state=s1, derivative_fn=derivative_fn, dt=h0)

    # Substep 3: S2(w1 * dt) starting at t + h1 + h0
    t3 = t2 + h0
    s3 = velocity_verlet_step(t=t3, state=s2, derivative_fn=derivative_fn, dt=h1)

    return s3


class Yoshida4Integrator:
    """Yoshida 4th-order symplectic numerical integrator class with standard interface."""

    def __init__(self, derivative_fn: Union[DerivativeCallable, GenericDerivativeCallable]):
        self.derivative_fn = derivative_fn

    def step(
        self,
        t: float,
        state: Union[State2D, Sequence[float]],
        dt: float,
    ) -> Union[State2D, list[float]]:
        """Advance state by timestep dt using Yoshida 4th-order method."""
        return yoshida4_step(
            t=t, state=state, derivative_fn=self.derivative_fn, dt=dt
        )
