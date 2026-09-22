"""Unit tests for Explicit Euler integrator."""

import math
from src.integrators.euler import euler_step, EulerIntegrator
from src.models.state import State2D, StateDerivative2D


def test_euler_step_linear_constant_velocity():
    """Explicit Euler single step for constant velocity motion (ax=0, ay=0).

    x(t + dt) = x(t) + vx * dt
    """
    state = State2D(x=10.0, y=20.0, vx=3.0, vy=-4.0)
    dt = 2.0

    # Constant velocity derivative function: acceleration is zero
    def const_vel_deriv(t, s):
        return StateDerivative2D(vx=s.vx, vy=s.vy, ax=0.0, ay=0.0)

    next_state = euler_step(0.0, state, const_vel_deriv, dt)

    assert math.isclose(next_state.x, 10.0 + 3.0 * 2.0, rel_tol=1e-12)
    assert math.isclose(next_state.y, 20.0 + (-4.0) * 2.0, rel_tol=1e-12)
    assert math.isclose(next_state.vx, 3.0, rel_tol=1e-12)
    assert math.isclose(next_state.vy, -4.0, rel_tol=1e-12)


def test_euler_step_with_acceleration():
    """Explicit Euler single step with known non-zero acceleration.

    vx_{n+1} = vx_n + ax * dt
    x_{n+1}  = x_n  + vx_n * dt
    """
    state = State2D(x=0.0, y=0.0, vx=10.0, vy=0.0)
    dt = 0.5
    constant_ax = -2.0
    constant_ay = 1.0

    def constant_accel_deriv(t, s):
        return StateDerivative2D(vx=s.vx, vy=s.vy, ax=constant_ax, ay=constant_ay)

    next_state = euler_step(0.0, state, constant_accel_deriv, dt)

    # In explicit Euler: position uses velocity at time n, not n+1
    assert math.isclose(next_state.x, 0.0 + 10.0 * 0.5, rel_tol=1e-12)
    assert math.isclose(next_state.y, 0.0 + 0.0 * 0.5, rel_tol=1e-12)
    # Velocity update
    assert math.isclose(next_state.vx, 10.0 + constant_ax * 0.5, rel_tol=1e-12)
    assert math.isclose(next_state.vy, 0.0 + constant_ay * 0.5, rel_tol=1e-12)


def test_euler_integrator_decoupled_from_domain():
    """EulerIntegrator should work with any generic 1D/scalar or list ODE.

    For example, radioactive/exponential decay: dy/dt = -lambda * y.
    Analytical Euler step: y_{n+1} = y_n * (1 - lambda * dt).
    """
    decay_rate = 0.5
    dt = 0.1
    y0 = [100.0]

    def decay_deriv(t, state_seq):
        return [-decay_rate * state_seq[0]]

    integrator = EulerIntegrator(decay_deriv)
    y1 = integrator.step(0.0, y0, dt)

    expected = 100.0 * (1.0 - decay_rate * dt)
    assert math.isclose(y1[0], expected, rel_tol=1e-12)

