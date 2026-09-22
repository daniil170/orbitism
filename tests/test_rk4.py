"""Unit tests for classical 4th-order Runge-Kutta (RK4) integrator."""

import math
import pytest

from src.integrators.euler import EulerIntegrator, euler_step
from src.integrators.rk4 import RK4Integrator, rk4_step
from src.models.state import State2D, StateDerivative2D
from src.physics.analytical import (
    circular_orbit_exact_state,
    circular_orbit_period,
    specific_orbital_energy,
)
from src.physics.gravity import (
    MU_EARTH,
    R_EARTH,
    create_circular_orbit_state,
    two_body_derivatives,
)
from src.simulation.metrics import position_error, relative_energy_error, velocity_error
from src.simulation.simulator import Simulator


def test_rk4_step_linear_constant_velocity():
    """RK4 must reproduce constant-velocity motion exactly (ax=0, ay=0).

    For the ODE:
        dx/dt = vx, dy/dt = vy, dvx/dt = 0, dvy/dt = 0
    The exact solution is linear:
        x(t + dt) = x(t) + vx * dt
        y(t + dt) = y(t) + vy * dt
    Since all derivatives of order >= 2 vanish identically, RK4 must match
    the analytical solution down to floating-point precision.
    """
    state = State2D(x=100.0, y=-50.0, vx=7.5, vy=-3.2)
    dt = 12.0

    def const_vel_deriv(t, s):
        return StateDerivative2D(vx=s.vx, vy=s.vy, ax=0.0, ay=0.0)

    next_state = rk4_step(0.0, state, const_vel_deriv, dt)

    assert math.isclose(next_state.x, 100.0 + 7.5 * dt, rel_tol=1e-14)
    assert math.isclose(next_state.y, -50.0 + (-3.2) * dt, rel_tol=1e-14)
    assert math.isclose(next_state.vx, 7.5, rel_tol=1e-14)
    assert math.isclose(next_state.vy, -3.2, rel_tol=1e-14)


def test_rk4_step_more_accurate_than_euler_single_step():
    """One RK4 step must be substantially more accurate than one Euler step.

    Tested on scalar exponential decay:
        dy/dt = -lambda * y, y(0) = y0
        y_exact(dt) = y0 * exp(-lambda * dt)

    Local truncation error:
        Euler: O(dt^2)
        RK4:   O(dt^5)
    """
    decay_rate = 0.5
    y0 = [100.0]
    dt = 0.2

    def decay_deriv(t, s):
        return [-decay_rate * s[0]]

    # Analytical solution
    exact_y = y0[0] * math.exp(-decay_rate * dt)

    # Euler step
    euler_y = euler_step(0.0, y0, decay_deriv, dt)[0]
    euler_error = abs(euler_y - exact_y)

    # RK4 step
    rk4_y = rk4_step(0.0, y0, decay_deriv, dt)[0]
    rk4_error = abs(rk4_y - exact_y)

    assert rk4_error < euler_error
    # RK4 LTE is O(dt^5) vs Euler O(dt^2); error ratio must be dramatic
    assert rk4_error < 0.01 * euler_error


def test_rk4_returns_finite_values():
    """RK4Integrator must produce strictly finite values (no NaN or Inf)."""
    initial_state = create_circular_orbit_state(altitude_km=400.0)
    sim = Simulator(
        initial_state=initial_state,
        dt=60.0,
        num_steps=100,
        integrator=RK4Integrator,
    )
    trajectory = sim.run()

    for state in trajectory:
        assert math.isfinite(state.x), f"Non-finite x: {state.x}"
        assert math.isfinite(state.y), f"Non-finite y: {state.y}"
        assert math.isfinite(state.vx), f"Non-finite vx: {state.vx}"
        assert math.isfinite(state.vy), f"Non-finite vy: {state.vy}"


def test_rk4_preserves_orbit_substantially_better_than_euler_short_horizon():
    """RK4 must preserve orbit substantially better than Euler over one orbit.

    Compares position error, velocity error, and energy conservation
    between Euler and RK4 for dt = 60.0 s over one orbital period T.
    """
    altitude_km = 400.0
    r0 = R_EARTH + altitude_km
    t_period = circular_orbit_period(r0, mu=MU_EARTH)
    initial_state = create_circular_orbit_state(altitude_km=altitude_km)
    dt = 60.0

    # Run Euler
    sim_euler = Simulator(
        initial_state=initial_state,
        dt=dt,
        t_final=t_period,
        integrator=EulerIntegrator,
    )
    traj_euler = sim_euler.run()

    # Run RK4
    sim_rk4 = Simulator(
        initial_state=initial_state,
        dt=dt,
        t_final=t_period,
        integrator=RK4Integrator,
    )
    traj_rk4 = sim_rk4.run()

    # Exact state at final time
    exact_state = circular_orbit_exact_state(t_period, r0=r0, mu=MU_EARTH)

    err_pos_euler = position_error(traj_euler[-1], exact_state)
    err_pos_rk4 = position_error(traj_rk4[-1], exact_state)

    err_vel_euler = velocity_error(traj_euler[-1], exact_state)
    err_vel_rk4 = velocity_error(traj_rk4[-1], exact_state)

    ref_energy = -MU_EARTH / (2.0 * r0)
    err_energy_euler = relative_energy_error(traj_euler[-1], ref_energy=ref_energy, mu=MU_EARTH)
    err_energy_rk4 = relative_energy_error(traj_rk4[-1], ref_energy=ref_energy, mu=MU_EARTH)

    # RK4 error must be substantially smaller without arbitrary fragile thresholds
    assert err_pos_rk4 < err_pos_euler
    assert err_vel_rk4 < err_vel_euler
    assert abs(err_energy_rk4) < abs(err_energy_euler)

    # Order of magnitude check: RK4 should outperform Euler by at least 100x at dt=60s
    assert err_pos_rk4 < 0.01 * err_pos_euler
    assert abs(err_energy_rk4) < 0.01 * abs(err_energy_euler)


def test_rk4_generic_sequence_ode():
    """RK4Integrator should work with generic list/sequence ODEs."""
    omega = 2.0
    dt = 0.05
    y0 = [1.0, 0.0]  # [x, vx] for simple harmonic oscillator: d2x/dt2 = -omega^2 x

    def sho_deriv(t, s):
        return [s[1], -(omega**2) * s[0]]

    integrator = RK4Integrator(sho_deriv)
    y1 = integrator.step(0.0, y0, dt)

    exact_x = math.cos(omega * dt)
    exact_vx = -omega * math.sin(omega * dt)

    assert math.isclose(y1[0], exact_x, rel_tol=1e-5)
    assert math.isclose(y1[1], exact_vx, rel_tol=1e-5)


def test_simulator_supports_rk4_instance_and_class():
    """Simulator should accept RK4Integrator class or instance."""
    initial_state = create_circular_orbit_state(altitude_km=400.0)

    # Pass class
    sim_cls = Simulator(initial_state=initial_state, dt=1.0, num_steps=5, integrator=RK4Integrator)
    traj_cls = sim_cls.run()
    assert isinstance(sim_cls.integrator, RK4Integrator)

    # Pass instance
    custom_fn = lambda t, s: two_body_derivatives(t, s, mu=MU_EARTH)
    inst = RK4Integrator(custom_fn)
    sim_inst = Simulator(initial_state=initial_state, dt=1.0, num_steps=5, integrator=inst)
    traj_inst = sim_inst.run()
    assert sim_inst.integrator is inst

    # Trajectories must match
    for s1, s2 in zip(traj_cls, traj_inst):
        assert math.isclose(s1.x, s2.x, rel_tol=1e-12)
        assert math.isclose(s1.y, s2.y, rel_tol=1e-12)


def test_simulator_default_remains_euler():
    """Simulator without explicit integrator argument must default to EulerIntegrator."""
    sim = Simulator()
    assert isinstance(sim.integrator, EulerIntegrator)
