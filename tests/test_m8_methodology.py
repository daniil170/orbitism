"""Tests for Milestone M8: Systematic Comparison methodology, reference solutions, and cost accounting.

Validates:
1. Exact analytical Keplerian circular reference solution and properties.
2. Invariant conservation of the analytical reference state.
3. Kinematic and decomposition error metric calculations.
4. Cost accounting models and matched evaluation budget calculations.
5. Schema consistency across planned M8 experimental artifacts.
"""

import math
import pytest

from src.models.state import State2D
from src.physics.gravity import MU_EARTH, R_EARTH
from src.physics.analytical import (
    circular_orbit_mean_motion,
    circular_orbit_period,
    circular_orbit_exact_state,
    specific_orbital_energy,
    specific_angular_momentum,
)
from src.simulation.metrics import (
    position_error,
    velocity_error,
    relative_energy_error,
    relative_angular_momentum_error,
    radial_drift,
    angular_phase_error,
    radial_position_error,
    along_track_position_error,
    radial_error_fraction,
    along_track_error_fraction,
    dominance_ratio,
)
from experiments.exp08_systematic_comparison import (
    CURRENT_IMPLEMENTATION_EVALS_PER_STEP,
    compute_force_evals_per_step,
    compute_total_force_evals,
    compute_dt_for_budget,
    PLANNED_MATCHED_BUDGETS,
    PLANNED_CONVERGENCE_COLUMNS,
    PLANNED_LONG_TERM_FIXED_DT_COLUMNS,
    PLANNED_ACCURACY_VS_COST_COLUMNS,
    PLANNED_INVARIANT_QUALITY_COLUMNS,
)


# ==============================================================================
# 1. ANALYTICAL REFERENCE SOLUTION TESTS
# ==============================================================================

class TestAnalyticalReferenceSolution:
    """Verifies that the exact circular Keplerian orbit satisfies all physical properties."""

    @pytest.fixture
    def orbit_params(self):
        r0 = R_EARTH + 400.0  # 6778.137 km
        mu = MU_EARTH         # 398600.435 km^3/s^2
        n = math.sqrt(mu / (r0**3))
        t_period = 2.0 * math.pi / n
        v0 = math.sqrt(mu / r0)
        return {"r0": r0, "mu": mu, "n": n, "t_period": t_period, "v0": v0}

    def test_mean_motion_and_period(self, orbit_params):
        r0 = orbit_params["r0"]
        mu = orbit_params["mu"]
        n_calc = circular_orbit_mean_motion(r0, mu=mu)
        t_calc = circular_orbit_period(r0, mu=mu)

        assert math.isclose(n_calc, orbit_params["n"], rel_tol=1e-14)
        assert math.isclose(t_calc, orbit_params["t_period"], rel_tol=1e-14)
        assert 5553.0 < t_calc < 5554.0  # ~92.56 minutes

    def test_exact_state_quadrants(self, orbit_params):
        r0 = orbit_params["r0"]
        mu = orbit_params["mu"]
        t_period = orbit_params["t_period"]
        v0 = orbit_params["v0"]

        # t = 0: (r0, 0), (0, v0)
        s0 = circular_orbit_exact_state(0.0, r0=r0, mu=mu)
        assert math.isclose(s0.x, r0, abs_tol=1e-9)
        assert math.isclose(s0.y, 0.0, abs_tol=1e-9)
        assert math.isclose(s0.vx, 0.0, abs_tol=1e-9)
        assert math.isclose(s0.vy, v0, abs_tol=1e-9)

        # t = T/4: (0, r0), (-v0, 0)
        s_quarter = circular_orbit_exact_state(0.25 * t_period, r0=r0, mu=mu)
        assert math.isclose(s_quarter.x, 0.0, abs_tol=1e-8)
        assert math.isclose(s_quarter.y, r0, abs_tol=1e-8)
        assert math.isclose(s_quarter.vx, -v0, abs_tol=1e-8)
        assert math.isclose(s_quarter.vy, 0.0, abs_tol=1e-8)

        # t = T/2: (-r0, 0), (0, -v0)
        s_half = circular_orbit_exact_state(0.5 * t_period, r0=r0, mu=mu)
        assert math.isclose(s_half.x, -r0, abs_tol=1e-8)
        assert math.isclose(s_half.y, 0.0, abs_tol=1e-8)
        assert math.isclose(s_half.vx, 0.0, abs_tol=1e-8)
        assert math.isclose(s_half.vy, -v0, abs_tol=1e-8)

        # t = T: completes full orbit, returns to s0
        s_full = circular_orbit_exact_state(t_period, r0=r0, mu=mu)
        assert math.isclose(s_full.x, s0.x, abs_tol=1e-8)
        assert math.isclose(s_full.y, s0.y, abs_tol=1e-8)
        assert math.isclose(s_full.vx, s0.vx, abs_tol=1e-8)
        assert math.isclose(s_full.vy, s0.vy, abs_tol=1e-8)

    def test_exact_state_invariants_constancy(self, orbit_params):
        r0 = orbit_params["r0"]
        mu = orbit_params["mu"]
        t_period = orbit_params["t_period"]

        ref_energy = -mu / (2.0 * r0)
        ref_h = r0 * orbit_params["v0"]

        test_times = [0.0, 0.1 * t_period, 0.33 * t_period, 1.0 * t_period, 10.5 * t_period]
        for t in test_times:
            state = circular_orbit_exact_state(t, r0=r0, mu=mu)
            assert math.isclose(state.r, r0, rel_tol=1e-14)
            assert math.isclose(state.v, orbit_params["v0"], rel_tol=1e-14)

            e = specific_orbital_energy(state, mu=mu)
            h = specific_angular_momentum(state)
            assert math.isclose(e, ref_energy, rel_tol=1e-14)
            assert math.isclose(h, ref_h, rel_tol=1e-14)

            # Invariant errors against reference should be zero to machine precision
            de = relative_energy_error(state, ref_energy=ref_energy, mu=mu)
            dh = relative_angular_momentum_error(state, ref_angular_momentum=ref_h)
            assert abs(de) < 1e-14
            assert abs(dh) < 1e-14


# ==============================================================================
# 2. METRIC CALCULATION TESTS
# ==============================================================================

class TestMetricCalculations:
    """Verifies kinematic and decomposed error metrics."""

    def test_position_and_velocity_error(self):
        s_exact = State2D(x=10.0, y=0.0, vx=0.0, vy=5.0)
        s_num = State2D(x=13.0, y=4.0, vx=0.3, vy=5.4)

        # dr = sqrt(3^2 + 4^2) = 5.0
        # dv = sqrt(0.3^2 + 0.4^2) = 0.5
        assert math.isclose(position_error(s_num, s_exact), 5.0, rel_tol=1e-14)
        assert math.isclose(velocity_error(s_num, s_exact), 0.5, rel_tol=1e-14)

    def test_radial_and_along_track_projections(self):
        r0 = 7000.0
        mu = MU_EARTH
        t = 0.0  # at t = 0, e_r = [1, 0], e_theta = [0, 1]

        s_exact = State2D(x=r0, y=0.0, vx=0.0, vy=math.sqrt(mu / r0))
        # Add +2 km radial shift and +5 km along-track shift
        s_num = State2D(x=r0 + 2.0, y=5.0, vx=0.0, vy=math.sqrt(mu / r0))

        e_R = radial_position_error(s_num, s_exact, t=t, r0=r0, mu=mu)
        e_T = along_track_position_error(s_num, s_exact, t=t, r0=r0, mu=mu)
        e_tot = position_error(s_num, s_exact)

        assert math.isclose(e_R, 2.0, rel_tol=1e-12)
        assert math.isclose(e_T, 5.0, rel_tol=1e-12)
        assert math.isclose(e_tot, math.sqrt(2.0**2 + 5.0**2), rel_tol=1e-12)

        # Variance fractions
        rho_R = radial_error_fraction(e_R, e_tot)
        rho_T = along_track_error_fraction(e_T, e_tot)
        assert math.isclose(rho_R + rho_T, 1.0, rel_tol=1e-12)
        assert math.isclose(dominance_ratio(e_R, e_T), 2.5, rel_tol=1e-12)

    def test_phase_and_radial_drift(self):
        r0 = 7000.0
        mu = MU_EARTH
        n = circular_orbit_mean_motion(r0, mu=mu)

        s_num = State2D(x=7010.0, y=0.0, vx=0.0, vy=7.5)
        assert math.isclose(radial_drift(s_num, r0=r0), 10.0, rel_tol=1e-12)

        # Phase error test
        t = 100.0
        expected_theta = n * t + 0.05  # 0.05 rad phase advance
        d_theta = angular_phase_error(theta_unwrapped=expected_theta, t=t, r0=r0, mu=mu)
        assert math.isclose(d_theta, 0.05, rel_tol=1e-12)


# ==============================================================================
# 3. COST ACCOUNTING AND BUDGET TESTS
# ==============================================================================

class TestCostAccounting:
    """Verifies cost models, evaluation multipliers, and matched budget assignments."""

    def test_evals_per_step_model(self):
        expected_evals = {
            "Euler": 1,
            "Velocity Verlet": 2,
            "RK4": 4,
            "Yoshida 4": 6,
        }
        for method, evals in expected_evals.items():
            assert compute_force_evals_per_step(method) == evals
            assert CURRENT_IMPLEMENTATION_EVALS_PER_STEP[method] == evals

        with pytest.raises(KeyError):
            compute_force_evals_per_step("NonExistentMethod")

    def test_compute_total_force_evals(self):
        assert compute_total_force_evals("Euler", 1000) == 1000
        assert compute_total_force_evals("Velocity Verlet", 1000) == 2000
        assert compute_total_force_evals("RK4", 1000) == 4000
        assert compute_total_force_evals("Yoshida 4", 1000) == 6000

    def test_compute_dt_for_budget(self):
        t_final = 60000.0
        target_evals = 6000

        # Euler: 6000 evals / 1 = 6000 steps => dt = 60000 / 6000 = 10 s
        dt_euler = compute_dt_for_budget("Euler", t_final=t_final, target_evals=target_evals)
        assert math.isclose(dt_euler, 10.0, rel_tol=1e-12)

        # Velocity Verlet: 6000 evals / 2 = 3000 steps => dt = 60000 / 3000 = 20 s
        dt_vv = compute_dt_for_budget("Velocity Verlet", t_final=t_final, target_evals=target_evals)
        assert math.isclose(dt_vv, 20.0, rel_tol=1e-12)

        # RK4: 6000 evals / 4 = 1500 steps => dt = 60000 / 1500 = 40 s
        dt_rk4 = compute_dt_for_budget("RK4", t_final=t_final, target_evals=target_evals)
        assert math.isclose(dt_rk4, 40.0, rel_tol=1e-12)

        # Yoshida 4: 6000 evals / 6 = 1000 steps => dt = 60000 / 1000 = 60 s
        dt_y4 = compute_dt_for_budget("Yoshida 4", t_final=t_final, target_evals=target_evals)
        assert math.isclose(dt_y4, 60.0, rel_tol=1e-12)

    def test_planned_matched_budgets_integrity(self):
        required_solvers = {"Euler", "Velocity Verlet", "RK4", "Yoshida 4"}
        for budget in PLANNED_MATCHED_BUDGETS:
            assert budget.budget_name.startswith("Budget_")
            assert budget.target_evals > 0
            assert set(budget.dt_assignments.keys()) == required_solvers

            # Verify that timestep ratios match cost multipliers:
            # dt_VV / dt_Euler = 2 / 1 = 2
            # dt_RK4 / dt_Euler = 4 / 1 = 4
            # dt_Y4 / dt_Euler = 6 / 1 = 6
            dt_e = budget.dt_assignments["Euler"]
            dt_vv = budget.dt_assignments["Velocity Verlet"]
            dt_rk4 = budget.dt_assignments["RK4"]
            dt_y4 = budget.dt_assignments["Yoshida 4"]

            assert math.isclose(dt_vv / dt_e, 2.0, rel_tol=1e-9)
            assert math.isclose(dt_rk4 / dt_e, 4.0, rel_tol=1e-9)
            assert math.isclose(dt_y4 / dt_e, 6.0, rel_tol=1e-9)


# ==============================================================================
# 4. CSV SCHEMA SPECIFICATION TESTS
# ==============================================================================

class TestCSVSchemas:
    """Verifies planned CSV column definitions for M8 benchmark outputs."""

    def test_convergence_schema(self):
        assert "method" in PLANNED_CONVERGENCE_COLUMNS
        assert "dt" in PLANNED_CONVERGENCE_COLUMNS
        assert "order_p_r" in PLANNED_CONVERGENCE_COLUMNS
        assert "force_evals" in PLANNED_CONVERGENCE_COLUMNS
        assert len(PLANNED_CONVERGENCE_COLUMNS) >= 10

    def test_long_term_fixed_dt_schema(self):
        assert "method" in PLANNED_LONG_TERM_FIXED_DT_COLUMNS
        assert "dt" in PLANNED_LONG_TERM_FIXED_DT_COLUMNS
        assert "actual_force_evals" in PLANNED_LONG_TERM_FIXED_DT_COLUMNS
        assert "max_relative_energy_error" in PLANNED_LONG_TERM_FIXED_DT_COLUMNS
        assert "final_phase_error" in PLANNED_LONG_TERM_FIXED_DT_COLUMNS

    def test_cost_normalized_schema(self):
        assert "budget_name" in PLANNED_ACCURACY_VS_COST_COLUMNS
        assert "target_evals" in PLANNED_ACCURACY_VS_COST_COLUMNS
        assert "actual_force_evals" in PLANNED_ACCURACY_VS_COST_COLUMNS
        assert "max_position_error" in PLANNED_ACCURACY_VS_COST_COLUMNS

    def test_invariant_quality_schema(self):
        assert "method" in PLANNED_INVARIANT_QUALITY_COLUMNS
        assert "de_character" in PLANNED_INVARIANT_QUALITY_COLUMNS
        assert "dh_character" in PLANNED_INVARIANT_QUALITY_COLUMNS
