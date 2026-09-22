"""Convergence rate verification for classical 4th-order Runge-Kutta (RK4).

Evaluates empirical convergence order p via systematic mesh refinement:
    dt, dt / 2, dt / 4
over one full orbital period T on a 2D circular Low Earth Orbit.

Convergence formula:
    p = ln(E(dt_1) / E(dt_2)) / ln(dt_1 / dt_2)

Theoretical expectation:
    p ≈ 4.0 (global truncation error O(dt^4))
"""

import pytest

from src.integrators.rk4 import RK4Integrator
from src.physics.analytical import circular_orbit_exact_state, circular_orbit_period
from src.physics.gravity import MU_EARTH, R_EARTH, create_circular_orbit_state
from src.simulation.metrics import convergence_order, position_error, velocity_error
from src.simulation.simulator import Simulator


def test_rk4_convergence_order_mesh_refinement():
    """Verify empirical convergence order p ≈ 4 via systematic step refinement.

    Refinement triplet:
        dt1 = 60.0 s
        dt2 = 30.0 s (dt1 / 2)
        dt3 = 15.0 s (dt1 / 4)

    Scientific considerations:
    - In the asymptotic regime dt -> 0, classical RK4 exhibits global truncation
      error O(dt^4), corresponding to convergence order p = 4.0.
    - At moderate timesteps (e.g. 60 s), higher-order Taylor series residuals
      and Keplerian nonlinear orbital dynamics cause empirical p to approach 4
      from slightly above (~4.22 -> ~4.12 -> ~4.06).
    - Terminal step shortening min(dt, T - t) also introduces a minor sub-step
      at t = T, which does not disrupt the 4th-order asymptotic scaling.
    """
    altitude_km = 400.0
    r0 = R_EARTH + altitude_km
    t_final = circular_orbit_period(r0, mu=MU_EARTH)
    initial_state = create_circular_orbit_state(altitude_km=altitude_km)

    dt_triplet = [60.0, 30.0, 15.0]
    pos_errors = []
    vel_errors = []

    for dt in dt_triplet:
        sim = Simulator(
            initial_state=initial_state,
            dt=dt,
            t_final=t_final,
            integrator=RK4Integrator,
        )
        trajectory = sim.run()
        exact_state = circular_orbit_exact_state(sim.times[-1], r0=r0, mu=MU_EARTH)

        e_pos = position_error(trajectory[-1], exact_state)
        e_vel = velocity_error(trajectory[-1], exact_state)
        pos_errors.append(e_pos)
        vel_errors.append(e_vel)

    # 1. Error monotonic reduction by factor of ~16 (2^4)
    # E(30) / E(60) ≈ 1/16 ≈ 0.0625
    ratio_pos_1 = pos_errors[1] / pos_errors[0]
    ratio_pos_2 = pos_errors[2] / pos_errors[1]
    assert ratio_pos_1 < 0.1, f"Expected ~16x error reduction, got ratio {ratio_pos_1}"
    assert ratio_pos_2 < 0.1, f"Expected ~16x error reduction, got ratio {ratio_pos_2}"

    # 2. Compute empirical order p
    p_pos_1 = convergence_order(pos_errors[0], pos_errors[1], dt_triplet[0], dt_triplet[1])
    p_pos_2 = convergence_order(pos_errors[1], pos_errors[2], dt_triplet[1], dt_triplet[2])

    p_vel_1 = convergence_order(vel_errors[0], vel_errors[1], dt_triplet[0], dt_triplet[1])
    p_vel_2 = convergence_order(vel_errors[1], vel_errors[2], dt_triplet[1], dt_triplet[2])

    # 3. Scientific expectation: p ≈ 4.0 within realistic engineering/numerical tolerance
    # Fact: p_pos_1 ≈ 4.218, p_pos_2 ≈ 4.122
    # Fact: p_vel_1 ≈ 4.214, p_vel_2 ≈ 4.122
    assert p_pos_1 == pytest.approx(4.0, abs=0.3), f"Position order p1 was {p_pos_1:.4f}"
    assert p_pos_2 == pytest.approx(4.0, abs=0.3), f"Position order p2 was {p_pos_2:.4f}"
    assert p_vel_1 == pytest.approx(4.0, abs=0.3), f"Velocity order p1 was {p_vel_1:.4f}"
    assert p_vel_2 == pytest.approx(4.0, abs=0.3), f"Velocity order p2 was {p_vel_2:.4f}"


def test_rk4_convergence_asymptotic_tightening():
    """Verify that as dt becomes finer, empirical order p tightens towards 4.0.

    Evaluates:
        p(30 -> 15) vs p(15 -> 7.5)
    As dt decreases, higher-order Taylor truncation terms diminish,
    and p asymptotically converges towards 4.00.
    """
    altitude_km = 400.0
    r0 = R_EARTH + altitude_km
    t_final = circular_orbit_period(r0, mu=MU_EARTH)
    initial_state = create_circular_orbit_state(altitude_km=altitude_km)

    dt_values = [30.0, 15.0, 7.5]
    errors = []

    for dt in dt_values:
        sim = Simulator(
            initial_state=initial_state,
            dt=dt,
            t_final=t_final,
            integrator=RK4Integrator,
        )
        trajectory = sim.run()
        exact = circular_orbit_exact_state(sim.times[-1], r0=r0, mu=MU_EARTH)
        errors.append(position_error(trajectory[-1], exact))

    p1 = convergence_order(errors[0], errors[1], dt_values[0], dt_values[1])  # 30 -> 15: ~4.122
    p2 = convergence_order(errors[1], errors[2], dt_values[1], dt_values[2])  # 15 -> 7.5: ~4.065

    # Order must be closer to 4.0 for the finer interval
    assert abs(p2 - 4.0) < abs(p1 - 4.0)
    assert p2 == pytest.approx(4.0, abs=0.15)
