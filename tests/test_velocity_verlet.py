"""Mathematical and numerical tests for Velocity Verlet integrator."""

import math
import pytest

try:
    import numpy as np
except ImportError:
    class _NPCompat:
        """Lightweight compatibility wrapper when numpy is not installed."""
        @staticmethod
        def isfinite(val: float) -> bool:
            return math.isfinite(val)

    np = _NPCompat()

from src.integrators.euler import EulerIntegrator
from src.integrators.velocity_verlet import (
    VelocityVerletIntegrator,
    velocity_verlet_step,
)
from src.models.state import State2D, StateDerivative2D
from src.physics.analytical import (
    circular_orbit_exact_state,
    circular_orbit_period,
)
from src.physics.gravity import (
    MU_EARTH,
    R_EARTH,
    create_circular_orbit_state,
    two_body_derivatives,
)
from src.simulation.metrics import convergence_order, position_error
from src.simulation.simulator import Simulator


def test_t1_single_step_correctness_constant_acceleration():
    """T1: Verify single-step correctness against known analytical solution.

    For constant acceleration a = (gx, gy):
        dx/dt = vx, dvx/dt = gx
        dy/dt = vy, dvy/dt = gy
    The exact solution is quadratic in time:
        x(t + dt) = x(t) + vx * dt + 0.5 * gx * dt^2
        vx(t + dt) = vx(t) + gx * dt
    Velocity Verlet integrates constant acceleration motion with zero truncation error.
    """
    gx, gy = -9.81, 3.5
    state = State2D(x=10.0, y=20.0, vx=2.5, vy=-1.5)
    dt = 0.5

    def const_accel_deriv(t, s):
        return StateDerivative2D(vx=s.vx, vy=s.vy, ax=gx, ay=gy)

    next_state = velocity_verlet_step(0.0, state, const_accel_deriv, dt)

    expected_x = 10.0 + 2.5 * dt + 0.5 * gx * (dt**2)
    expected_y = 20.0 + (-1.5) * dt + 0.5 * gy * (dt**2)
    expected_vx = 2.5 + gx * dt
    expected_vy = -1.5 + gy * dt

    assert math.isclose(next_state.x, expected_x, rel_tol=1e-14, abs_tol=1e-14)
    assert math.isclose(next_state.y, expected_y, rel_tol=1e-14, abs_tol=1e-14)
    assert math.isclose(next_state.vx, expected_vx, rel_tol=1e-14, abs_tol=1e-14)
    assert math.isclose(next_state.vy, expected_vy, rel_tol=1e-14, abs_tol=1e-14)


def test_t1_single_step_harmonic_oscillator():
    """T1: Verify single-step correctness against harmonic oscillator formula.

    For d2x/dt2 = -omega^2 * x with x0 = 1.0, v0 = 0.0:
        a0 = -omega^2 * x0
        x1 = x0 + v0 * dt + 0.5 * a0 * dt^2 = 1.0 - 0.5 * omega^2 * dt^2
        a1 = -omega^2 * x1
        v1 = v0 + 0.5 * (a0 + a1) * dt
    """
    omega = 1.5
    dt = 0.1
    x0, v0 = 1.0, 0.0
    state = [x0, v0]

    def sho_deriv(t, s):
        return [s[1], -(omega**2) * s[0]]

    next_state = velocity_verlet_step(0.0, state, sho_deriv, dt)

    expected_x1 = x0 + v0 * dt - 0.5 * (omega**2) * x0 * (dt**2)
    expected_a1 = -(omega**2) * expected_x1
    expected_v1 = v0 + 0.5 * (-(omega**2) * x0 + expected_a1) * dt

    assert math.isclose(next_state[0], expected_x1, rel_tol=1e-14)
    assert math.isclose(next_state[1], expected_v1, rel_tol=1e-14)


def test_t2_finite_value_validation():
    """T2: Require np.isfinite(...) for all state variables over an orbit."""
    initial_state = create_circular_orbit_state(altitude_km=400.0)
    sim = Simulator(
        initial_state=initial_state,
        dt=60.0,
        num_steps=100,
        integrator=VelocityVerletIntegrator,
    )
    trajectory = sim.run()

    for s in trajectory:
        assert np.isfinite(s.x), f"Non-finite x: {s.x}"
        assert np.isfinite(s.y), f"Non-finite y: {s.y}"
        assert np.isfinite(s.vx), f"Non-finite vx: {s.vx}"
        assert np.isfinite(s.vy), f"Non-finite vy: {s.vy}"


def test_t3_expected_convergence_order():
    """T3: Empirical convergence order p must satisfy 1.8 <= p <= 2.2."""
    altitude_km = 400.0
    r0 = R_EARTH + altitude_km
    t_period = circular_orbit_period(r0, mu=MU_EARTH)
    initial_state = create_circular_orbit_state(altitude_km=altitude_km)

    dt_1 = 30.0
    dt_2 = 15.0

    sim_1 = Simulator(
        initial_state=initial_state,
        dt=dt_1,
        t_final=t_period,
        integrator=VelocityVerletIntegrator,
    )
    traj_1 = sim_1.run()

    sim_2 = Simulator(
        initial_state=initial_state,
        dt=dt_2,
        t_final=t_period,
        integrator=VelocityVerletIntegrator,
    )
    traj_2 = sim_2.run()

    exact_state = circular_orbit_exact_state(t_period, r0=r0, mu=MU_EARTH)
    err_1 = position_error(traj_1[-1], exact_state)
    err_2 = position_error(traj_2[-1], exact_state)

    p = convergence_order(err_1, err_2, dt_1, dt_2)

    assert 1.8 <= p <= 2.2, f"Expected 1.8 <= p <= 2.2, got p={p:.4f}"


def test_t4_time_reversibility():
    """T4: Integrate N forward steps then N backward steps; state_final ≈ state_initial."""
    initial_state = create_circular_orbit_state(altitude_km=400.0)
    deriv_fn = lambda t, s: two_body_derivatives(t, s, mu=MU_EARTH)
    integrator = VelocityVerletIntegrator(deriv_fn)

    dt = 30.0
    n_steps = 100

    current_state = initial_state
    current_time = 0.0

    # Step forward N steps
    for _ in range(n_steps):
        current_state = integrator.step(current_time, current_state, dt)
        current_time += dt

    # Step backward N steps with negative timestep -dt
    for _ in range(n_steps):
        current_state = integrator.step(current_time, current_state, -dt)
        current_time -= dt

    # Check that position and velocity return to initial state within numerical roundoff tolerance
    assert math.isclose(current_state.x, initial_state.x, rel_tol=1e-10, abs_tol=1e-8)
    assert math.isclose(current_state.y, initial_state.y, rel_tol=1e-10, abs_tol=1e-8)
    assert math.isclose(current_state.vx, initial_state.vx, rel_tol=1e-10, abs_tol=1e-8)
    assert math.isclose(current_state.vy, initial_state.vy, rel_tol=1e-10, abs_tol=1e-8)


def test_t5_long_run_robustness_100T():
    """T5: Run 100T for largest timestep (60s); require no NaN, no Inf, positive radius."""
    altitude_km = 400.0
    r0 = R_EARTH + altitude_km
    t_period = circular_orbit_period(r0, mu=MU_EARTH)
    initial_state = create_circular_orbit_state(altitude_km=altitude_km)

    dt = 60.0
    t_final = 100.0 * t_period

    sim = Simulator(
        initial_state=initial_state,
        dt=dt,
        t_final=t_final,
        integrator=VelocityVerletIntegrator,
    )
    trajectory = sim.run()

    assert len(trajectory) > 0

    for s in trajectory:
        assert math.isfinite(s.x), f"Non-finite x: {s.x}"
        assert math.isfinite(s.y), f"Non-finite y: {s.y}"
        assert math.isfinite(s.vx), f"Non-finite vx: {s.vx}"
        assert math.isfinite(s.vy), f"Non-finite vy: {s.vy}"

        r = math.hypot(s.x, s.y)
        assert r > 0.0, f"Non-positive radius: {r}"
        assert r > R_EARTH * 0.5, f"Unphysical orbit collapse: r={r}"


def test_simulator_supports_velocity_verlet_instance_and_class():
    """Simulator should accept VelocityVerletIntegrator class or instance."""
    initial_state = create_circular_orbit_state(altitude_km=400.0)

    # Pass class
    sim_cls = Simulator(
        initial_state=initial_state,
        dt=1.0,
        num_steps=5,
        integrator=VelocityVerletIntegrator,
    )
    traj_cls = sim_cls.run()
    assert isinstance(sim_cls.integrator, VelocityVerletIntegrator)

    # Pass instance
    custom_fn = lambda t, s: two_body_derivatives(t, s, mu=MU_EARTH)
    inst = VelocityVerletIntegrator(custom_fn)
    sim_inst = Simulator(
        initial_state=initial_state, dt=1.0, num_steps=5, integrator=inst
    )
    traj_inst = sim_inst.run()
    assert sim_inst.integrator is inst

    for s1, s2 in zip(traj_cls, traj_inst):
        assert math.isclose(s1.x, s2.x, rel_tol=1e-12)
        assert math.isclose(s1.y, s2.y, rel_tol=1e-12)
