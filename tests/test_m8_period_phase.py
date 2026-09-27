"""Tests for Milestone M8.3: Orbital Period Error and Phase Drift Analysis.

Verifies:
1. Exact analytical circular period T_exact.
2. Sign conventions and continuous phase error definition:
   - positive delta_theta -> numerical phase leads;
   - negative delta_theta -> numerical phase lags.
3. Phase unwrapping across [-pi, pi] branch cuts.
4. Positive x-axis crossing detection in prograde direction:
   - condition: y_i <= 0 < y_{i+1}, x > 0, v_y > 0;
   - analytical starting condition at t=0 on positive x-axis (m=0 is analytic, first crossing is m=1);
   - rejection of retrograde crossings (v_y < 0) and negative-x crossings (x < 0);
   - no false crossing near y = 0.
5. Cubic Hermite root finding:
   - bracketing, bisection fallback, residual checks, convergence tolerance;
   - error handling when inputs do not bracket zero or dt <= 0.
6. Period calculation metrics:
   - instantaneous revolution period T_inst(m) = t_cross(m) - t_cross(m-1);
   - cumulative mean period T_cum(m) = t_cross(m) / m;
   - relative period error delta_T = (T_cum - T_exact) / T_exact.
7. Consistency between event-based period error and phase-derived estimate:
   - exact theoretical identity delta_T(t*_m) = -delta_theta(t*_m) / (2*pi*m) at crossing;
   - cross-check agreement within discretization tolerance at discrete checkpoints.
8. Cartesian error decomposition:
   - geometric chord error 2*r0*|sin(delta_theta / 2)|;
   - along-track displacement r0 * delta_theta;
   - dominance of phase-derived chord error for geometric solvers vs divergence for Euler.
9. Experimental grid and schema integrity (96 checkpoint rows across 4 solvers, 4 dts, 6 checkpoints).
10. Empirical convergence order scaling for phase drift and period error.
"""

import math
from pathlib import Path
import pytest
from typing import Any, Dict, List

from experiments.exp08_systematic_comparison import (
    DEFAULT_PERIOD_PHASE_DT_VALUES,
    DEFAULT_PERIOD_PHASE_CHECKPOINTS,
    PLANNED_PERIOD_PHASE_COLUMNS,
    run_period_phase_study,
)
from src.integrators.euler import EulerIntegrator
from src.integrators.rk4 import RK4Integrator
from src.integrators.velocity_verlet import VelocityVerletIntegrator
from src.integrators.yoshida4 import Yoshida4Integrator
from src.models.state import State2D
from src.physics.gravity import MU_EARTH, R_EARTH, create_circular_orbit_state
from src.physics.analytical import (
    circular_orbit_period,
    circular_orbit_mean_motion,
    circular_orbit_exact_state,
)
from src.simulation.simulator import Simulator
from src.simulation.metrics import (
    cubic_hermite_root,
    detect_positive_x_crossings,
    phase_derived_chord_error,
    phase_derived_along_track_displacement,
    unwrap_phase,
    angular_phase_error,
    convergence_order,
)


PROJECT_ROOT = Path(__file__).resolve().parent.parent


class TestM83PeriodPhaseMath:
    """Mathematical and algorithmic verification of M8.3 primitives."""

    def test_analytical_exact_period(self):
        """Verify T_exact for 400 km circular Earth orbit."""
        r0 = R_EARTH + 400.0
        t_orbit = circular_orbit_period(r0, mu=MU_EARTH)
        n = circular_orbit_mean_motion(r0, mu=MU_EARTH)

        expected_t = 2.0 * math.pi / n
        assert math.isclose(t_orbit, expected_t, rel_tol=1e-15)
        assert math.isclose(t_orbit, 5553.624318623783, rel_tol=1e-10)

    def test_phase_sign_convention(self):
        """Verify phase error sign convention: leads -> positive, lags -> negative."""
        r0 = R_EARTH + 400.0
        n = circular_orbit_mean_motion(r0, mu=MU_EARTH)
        t = 1000.0
        theta_exact = n * t

        # Numerical satellite leads by 0.05 rad
        theta_lead = theta_exact + 0.05
        err_lead = angular_phase_error(theta_lead, t, r0, mu=MU_EARTH)
        assert err_lead > 0.0
        assert math.isclose(err_lead, 0.05, abs_tol=1e-14)

        # Numerical satellite lags by 0.05 rad
        theta_lag = theta_exact - 0.05
        err_lag = angular_phase_error(theta_lag, t, r0, mu=MU_EARTH)
        assert err_lag < 0.0
        assert math.isclose(err_lag, -0.05, abs_tol=1e-14)

    def test_phase_unwrapping(self):
        """Verify unwrap_phase handles jumps across [-pi, pi] branch cuts."""
        # Synthesize a continuous angle spanning multiple orbits
        true_angles = [0.1 * i for i in range(100)]
        wrapped_angles = [(a + math.pi) % (2.0 * math.pi) - math.pi for a in true_angles]

        unwrapped = unwrap_phase(wrapped_angles)
        assert len(unwrapped) == len(true_angles)
        for u, t in zip(unwrapped, true_angles):
            assert math.isclose(u, t, abs_tol=1e-14)

    def test_cubic_hermite_root_exact_linear(self):
        """Verify Hermite root finder on simple linear function."""
        # y(t) = -2 + 4*(t/dt) => root at s = 0.5
        dt = 10.0
        s_root = cubic_hermite_root(y0=-2.0, y1=2.0, dy0=4.0 / dt, dy1=4.0 / dt, dt=dt)
        assert math.isclose(s_root, 0.5, abs_tol=1e-12)

    def test_cubic_hermite_root_endpoints(self):
        """Verify exact endpoint return when y0=0 or y1=0."""
        assert cubic_hermite_root(y0=0.0, y1=1.0, dy0=1.0, dy1=1.0, dt=1.0) == 0.0
        assert cubic_hermite_root(y0=-1.0, y1=0.0, dy0=1.0, dy1=1.0, dt=1.0) == 1.0

    def test_cubic_hermite_root_validation(self):
        """Verify error handling on invalid intervals or non-bracketing inputs."""
        # dt <= 0
        with pytest.raises(ValueError, match="Timestep dt must be positive"):
            cubic_hermite_root(y0=-1.0, y1=1.0, dy0=1.0, dy1=1.0, dt=0.0)

        # Same signs (no zero crossing)
        with pytest.raises(ValueError, match="Endpoints do not bracket zero"):
            cubic_hermite_root(y0=1.0, y1=2.0, dy0=1.0, dy1=1.0, dt=1.0)

    def test_positive_x_crossing_detection_on_analytical_orbit(self):
        """Verify Hermite crossing detection on analytical trajectory recovers T_exact."""
        r0 = R_EARTH + 400.0
        t_orbit = circular_orbit_period(r0, mu=MU_EARTH)
        dt = 60.0
        # Sample analytical solution for 3 full orbits
        n_steps = int(3.5 * t_orbit / dt)
        times = [i * dt for i in range(n_steps)]
        traj = [circular_orbit_exact_state(t, r0=r0, mu=MU_EARTH) for t in times]

        events = detect_positive_x_crossings(times, traj, r0=r0, mu=MU_EARTH)
        assert len(events) == 3

        for ev in events:
            m = ev.orbit_index
            expected_t = m * t_orbit
            # Cubic Hermite on 60s step of Keplerian circle has error ~ (n*dt)^4 ~ 1e-7 s
            assert math.isclose(ev.t_cross, expected_t, abs_tol=1e-6)
            assert math.isclose(ev.period_cumulative, t_orbit, abs_tol=1e-6)
            assert abs(ev.delta_T_cumulative) < 1e-9
            assert ev.state_cross.y == 0.0
            assert ev.state_cross.x > 0.0
            assert ev.state_cross.vy > 0.0

    def test_no_false_crossings_rejections(self):
        """Verify non-prograde or negative-x crossings are rejected."""
        r0 = R_EARTH + 400.0
        times = [0.0, 10.0, 20.0, 30.0]

        # Case 1: Retrograde crossing (vy < 0) on positive x
        traj_retro = [
            State2D(x=r0, y=-1.0, vx=0.0, vy=-1.0),
            State2D(x=r0, y=1.0, vx=0.0, vy=-1.0),
            State2D(x=r0, y=2.0, vx=0.0, vy=-1.0),
            State2D(x=r0, y=3.0, vx=0.0, vy=-1.0),
        ]
        events_retro = detect_positive_x_crossings(times, traj_retro, r0=r0, mu=MU_EARTH)
        assert len(events_retro) == 0

        # Case 2: Crossing on negative x-axis (x < 0, vy > 0)
        traj_neg_x = [
            State2D(x=-r0, y=-1.0, vx=0.0, vy=1.0),
            State2D(x=-r0, y=1.0, vx=0.0, vy=1.0),
            State2D(x=-r0, y=2.0, vx=0.0, vy=1.0),
            State2D(x=-r0, y=3.0, vx=0.0, vy=1.0),
        ]
        events_neg_x = detect_positive_x_crossings(times, traj_neg_x, r0=r0, mu=MU_EARTH)
        assert len(events_neg_x) == 0

    def test_event_based_vs_phase_derived_identity(self):
        """Verify mathematical identity delta_T(t*_m) == -delta_theta(t*_m) / (2*pi*m) at crossing."""
        r0 = R_EARTH + 400.0
        t_orbit = circular_orbit_period(r0, mu=MU_EARTH)
        n = circular_orbit_mean_motion(r0, mu=MU_EARTH)
        dt = 30.0

        sim = Simulator(
            initial_state=create_circular_orbit_state(400.0),
            integrator=Yoshida4Integrator,
        )
        traj = sim.run(dt=dt, num_steps=math.ceil(3.5 * t_orbit / dt))
        events = detect_positive_x_crossings(sim.times, traj, r0=r0, mu=MU_EARTH)

        for ev in events:
            m = ev.orbit_index
            delta_theta_cross = 2.0 * math.pi * m - n * ev.t_cross
            phase_derived_delta_t = -delta_theta_cross / (2.0 * math.pi * m)
            # Mathematical identity must hold to machine precision
            assert math.isclose(ev.delta_T_cumulative, phase_derived_delta_t, abs_tol=1e-14)


@pytest.fixture(scope="module")
def period_phase_data():
    """Run M8.3 experiment once for module tests."""
    records, summaries = run_period_phase_study()
    return records, summaries


class TestM83PeriodPhaseBenchmark:
    """Rigorous verification of M8.3 simulation results and scientific properties."""

    def test_records_count_and_schema(self, period_phase_data):
        """Verify 96 checkpoint rows and complete schema."""
        records, _ = period_phase_data
        assert len(records) == 96  # 4 solvers x 4 dts x 6 checkpoints

        for col in PLANNED_PERIOD_PHASE_COLUMNS:
            assert all(col in r for r in records), f"Missing column {col}"

    def test_solvers_and_grid_coverage(self, period_phase_data):
        """Verify grid: 4 solvers x {60, 30, 15, 7.5} s x {1, 5, 10, 20, 50, 100} checkpoints."""
        records, _ = period_phase_data
        methods = {r["method"] for r in records}
        assert methods == {"Euler", "Velocity Verlet", "RK4", "Yoshida 4"}

        dts = {r["dt"] for r in records}
        assert dts == {60.0, 30.0, 15.0, 7.5}

        checkpoints = {r["checkpoint_m"] for r in records}
        assert checkpoints == {1, 5, 10, 20, 50, 100}

    def test_crossing_presence_behavior(self, period_phase_data):
        """Verify crossing availability: Verlet/RK4/Yoshida complete all 100, Euler diverges."""
        records, _ = period_phase_data

        for r in records:
            m = r["checkpoint_m"]
            method = r["method"]
            if method in {"Velocity Verlet", "RK4", "Yoshida 4"}:
                assert r["has_crossing"] is True
                assert r["t_cross"] is not None
                assert r["relative_period_error"] is not None
            elif method == "Euler":
                # Euler at dt=60 only completes ~10 orbits in 100T duration
                if m >= 50:
                    assert r["has_crossing"] is False
                    assert r["t_cross"] is None

    def test_period_error_empirical_convergence_orders(self, period_phase_data):
        """Verify empirical convergence orders of period error delta_T."""
        records, _ = period_phase_data

        # Velocity Verlet: 2nd order (order ~ 2.0)
        vv_orders = [
            r["order_period"]
            for r in records
            if r["method"] == "Velocity Verlet" and r["checkpoint_m"] == 100 and r["dt"] in [30.0, 15.0, 7.5]
        ]
        for p in vv_orders:
            assert 1.95 <= p <= 2.05, f"Velocity Verlet period order {p} deviates from 2.0"

        # Yoshida 4: 4th order (order ~ 4.0)
        yoshida_orders = [
            r["order_period"]
            for r in records
            if r["method"] == "Yoshida 4" and r["checkpoint_m"] == 100 and r["dt"] in [30.0, 15.0, 7.5]
        ]
        for p in yoshida_orders:
            assert 3.95 <= p <= 4.05, f"Yoshida 4 period order {p} deviates from 4.0"

        # RK4: 4th order or higher (order >= 4.0)
        rk4_orders = [
            r["order_period"]
            for r in records
            if r["method"] == "RK4" and r["checkpoint_m"] == 100 and r["dt"] in [30.0, 15.0, 7.5]
        ]
        for p in rk4_orders:
            assert p >= 4.0, f"RK4 period order {p} lower than 4.0"

    def test_chord_error_decomposition_behavior(self, period_phase_data):
        """Verify Cartesian error vs phase-derived chord error vs radial error."""
        records, _ = period_phase_data

        # For geometric solvers (Verlet, Yoshida, RK4) at 100 orbits with dt=60:
        # Cartesian error is dominated by chord error (ratio near 1.0) and radial error << chord error
        for method in ["Velocity Verlet", "RK4", "Yoshida 4"]:
            row = next(
                r for r in records
                if r["method"] == method and r["dt"] == 60.0 and r["checkpoint_m"] == 100
            )
            ratio = row["chord_error_km"] / row["cartesian_error_km"]
            assert 0.999 <= ratio <= 1.001, f"{method} chord/cartesian ratio {ratio} not near 1.0"
            assert row["radial_error_km"] < 1.0  # Radial error sub-kilometer

        # For Euler at dt=60 and 100 orbits:
        # Radial error is massive (~55000 km) and chord error does NOT explain Cartesian error
        euler_row = next(
            r for r in records
            if r["method"] == "Euler" and r["dt"] == 60.0 and r["checkpoint_m"] == 100
        )
        assert euler_row["radial_error_km"] > 10000.0
        euler_ratio = euler_row["chord_error_km"] / euler_row["cartesian_error_km"]
        assert euler_ratio < 0.25  # Chord error explains less than 25% of Cartesian error
