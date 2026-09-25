"""Mathematical and numerical tests for Yoshida 4th-order symplectic integrator."""

import math
from typing import Sequence, Union
import pytest

try:
    import numpy as np
except ImportError:
    class _NPCompat:
        @staticmethod
        def isfinite(val: float) -> bool:
            return math.isfinite(val)

    np = _NPCompat()

from src.integrators.velocity_verlet import VelocityVerletIntegrator
from src.integrators.yoshida4 import (
    W0,
    W1,
    Yoshida4Integrator,
    yoshida4_step,
)
from src.models.state import State2D, StateDerivative2D
from src.physics.analytical import (
    circular_orbit_exact_state,
    circular_orbit_period,
    specific_orbital_energy,
    specific_angular_momentum,
)
from src.physics.gravity import (
    MU_EARTH,
    R_EARTH,
    create_circular_orbit_state,
    two_body_derivatives,
)
from src.simulation.metrics import (
    convergence_order,
    position_error,
    compute_decomposed_error_history,
)
from src.simulation.simulator import Simulator


def test_yoshida4_coefficients_mathematical_properties():
    """Verify algebraic properties of Yoshida composition coefficients."""
    # Consistency condition: sum of step weights must equal 1
    sum_w = 2.0 * W1 + W0
    assert math.isclose(sum_w, 1.0, rel_tol=1e-15, abs_tol=1e-15)

    # Order condition: elimination of the 3rd-order Lie error term
    cubic_sum = 2.0 * (W1**3) + (W0**3)
    assert abs(cubic_sum) < 1e-15

    # Negative substep confirmation
    assert W1 > 0.0, "w1 should be positive"
    assert W0 < 0.0, "w0 must be negative (time-reversal substep)"


def test_yoshida4_force_evaluations_count():
    """Verify that Yoshida 4 executes exactly 6 force evaluations per outer timestep."""
    eval_count = 0

    def counting_derivative(t: float, s: State2D) -> StateDerivative2D:
        nonlocal eval_count
        eval_count += 1
        return two_body_derivatives(t, s, mu=MU_EARTH)

    state = create_circular_orbit_state(altitude_km=400.0)
    _ = yoshida4_step(0.0, state, counting_derivative, dt=30.0)

    assert eval_count == 6, f"Expected exactly 6 force evaluations, got {eval_count}"


def test_yoshida4_single_step_finite_and_types():
    """Verify single step execution, type preservation, and finiteness."""
    state = create_circular_orbit_state(altitude_km=400.0)
    integrator = Yoshida4Integrator(
        lambda t, s: two_body_derivatives(t, s, mu=MU_EARTH)
    )

    next_state = integrator.step(0.0, state, dt=60.0)

    assert isinstance(next_state, State2D)
    assert math.isfinite(next_state.x)
    assert math.isfinite(next_state.y)
    assert math.isfinite(next_state.vx)
    assert math.isfinite(next_state.vy)


def test_yoshida4_sequence_support():
    """Verify generic sequence fallback matches State2D execution."""
    state = create_circular_orbit_state(altitude_km=400.0)
    seq_state = [state.x, state.y, state.vx, state.vy]

    def deriv_state2d(t: float, s: State2D) -> StateDerivative2D:
        return two_body_derivatives(t, s, mu=MU_EARTH)

    def deriv_seq(t: float, s: Sequence[float]) -> list[float]:
        d = two_body_derivatives(
            t, State2D(x=s[0], y=s[1], vx=s[2], vy=s[3]), mu=MU_EARTH
        )
        return d.to_list()

    next_s = yoshida4_step(0.0, state, deriv_state2d, dt=30.0)
    next_seq = yoshida4_step(0.0, seq_state, deriv_seq, dt=30.0)

    assert isinstance(next_seq, list)
    assert math.isclose(next_s.x, next_seq[0], rel_tol=1e-14)
    assert math.isclose(next_s.y, next_seq[1], rel_tol=1e-14)
    assert math.isclose(next_s.vx, next_seq[2], rel_tol=1e-14)
    assert math.isclose(next_s.vy, next_seq[3], rel_tol=1e-14)


def test_yoshida4_convergence_order():
    """Verify 4th-order convergence on Keplerian orbit (error ratio ~ 16, order ~ 4)."""
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
        integrator=Yoshida4Integrator,
    )
    traj_1 = sim_1.run()

    sim_2 = Simulator(
        initial_state=initial_state,
        dt=dt_2,
        t_final=t_period,
        integrator=Yoshida4Integrator,
    )
    traj_2 = sim_2.run()

    exact_state = circular_orbit_exact_state(t_period, r0=r0, mu=MU_EARTH)
    err_1 = position_error(traj_1[-1], exact_state)
    err_2 = position_error(traj_2[-1], exact_state)

    ratio = err_1 / err_2
    p = convergence_order(err_1, err_2, dt_1, dt_2)

    # 2^4 = 16.0; require order p between 3.8 and 4.2
    assert 15.0 <= ratio <= 17.0, f"Expected ratio ~ 16.0, got {ratio:.3f}"
    assert 3.8 <= p <= 4.2, f"Expected order ~ 4.0, got {p:.3f}"


def test_yoshida4_time_reversibility():
    """Verify time-reversibility: N forward steps then N backward steps return to initial state."""
    initial_state = create_circular_orbit_state(altitude_km=400.0)
    deriv_fn = lambda t, s: two_body_derivatives(t, s, mu=MU_EARTH)
    integrator = Yoshida4Integrator(deriv_fn)

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

    assert math.isclose(current_state.x, initial_state.x, rel_tol=1e-10, abs_tol=1e-8)
    assert math.isclose(current_state.y, initial_state.y, rel_tol=1e-10, abs_tol=1e-8)
    assert math.isclose(current_state.vx, initial_state.vx, rel_tol=1e-10, abs_tol=1e-8)
    assert math.isclose(current_state.vy, initial_state.vy, rel_tol=1e-10, abs_tol=1e-8)


def test_yoshida4_physical_invariants_conservation():
    """Verify conservation of angular momentum to machine precision and bounded energy."""
    altitude_km = 400.0
    r0 = R_EARTH + altitude_km
    t_period = circular_orbit_period(r0, mu=MU_EARTH)
    initial_state = create_circular_orbit_state(altitude_km=altitude_km)

    # 10 orbits at dt = 30s
    sim = Simulator(
        initial_state=initial_state,
        dt=30.0,
        t_final=10.0 * t_period,
        integrator=Yoshida4Integrator,
    )
    trajectory = sim.run()

    h0 = specific_angular_momentum(initial_state)
    e0 = specific_orbital_energy(initial_state, mu=MU_EARTH)

    max_delta_h_rel = 0.0
    max_delta_e_rel = 0.0

    for s in trajectory:
        h = specific_angular_momentum(s)
        e = specific_orbital_energy(s, mu=MU_EARTH)

        dh_rel = abs((h - h0) / h0)
        de_rel = abs((e - e0) / e0)

        max_delta_h_rel = max(max_delta_h_rel, dh_rel)
        max_delta_e_rel = max(max_delta_e_rel, de_rel)

    # Symplectic integrators in central force preserve angular momentum to machine precision
    assert max_delta_h_rel < 1e-13, f"Angular momentum drift too large: {max_delta_h_rel}"
    # Energy must oscillate within tight bound without secular drift
    assert max_delta_e_rel < 1e-9, f"Energy drift too large: {max_delta_e_rel}"


def test_yoshida4_simulator_integration():
    """Verify Simulator accepts Yoshida4Integrator as class and instance."""
    initial_state = create_circular_orbit_state(altitude_km=400.0)

    # Pass class
    sim_cls = Simulator(
        initial_state=initial_state,
        dt=10.0,
        num_steps=5,
        integrator=Yoshida4Integrator,
    )
    traj_cls = sim_cls.run()
    assert isinstance(sim_cls.integrator, Yoshida4Integrator)

    # Pass instance
    custom_fn = lambda t, s: two_body_derivatives(t, s, mu=MU_EARTH)
    inst = Yoshida4Integrator(custom_fn)
    sim_inst = Simulator(
        initial_state=initial_state, dt=10.0, num_steps=5, integrator=inst
    )
    traj_inst = sim_inst.run()
    assert sim_inst.integrator is inst

    for s1, s2 in zip(traj_cls, traj_inst):
        assert math.isclose(s1.x, s2.x, rel_tol=1e-14)
        assert math.isclose(s1.y, s2.y, rel_tol=1e-14)
