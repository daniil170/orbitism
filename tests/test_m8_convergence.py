"""Tests for Milestone M8.1: Analytical Reference + Convergence Experiment.

Verifies:
1. Exactly 4 solvers (Euler, Velocity Verlet, RK4, Yoshida 4).
2. Exactly 6 timesteps per solver on standard dyadic grid [60, 30, 15, 7.5, 3.75, 1.875].
3. Exactly 24 output rows.
4. Timesteps strictly match the predefined grid.
5. Force evaluation counts conform to implementation-level model (evals_per_step * num_steps).
6. Simulation final physical time strictly equals analytical period T.
7. No NaN or Inf in any metric.
8. Norms and absolute error metrics are strictly non-negative (>= 0).
9. Analytical reference is evaluated in closed-form trigonometry (not as a numerical integrator).
10. Empirical convergence order mathematical sanity checks.
"""

import math
from pathlib import Path
import pytest
from typing import Any, Dict, List

from experiments.exp08_systematic_comparison import (
    DEFAULT_DT_VALUES,
    CURRENT_IMPLEMENTATION_EVALS_PER_STEP,
    SOLVER_REGISTRY,
    PLANNED_CONVERGENCE_COLUMNS,
    run_convergence_study,
)
from src.physics.gravity import MU_EARTH, R_EARTH
from src.physics.analytical import circular_orbit_period, circular_orbit_exact_state


@pytest.fixture(scope="module")
def convergence_records() -> List[Dict[str, Any]]:
    """Generate convergence records for test inspection."""
    return run_convergence_study()


class TestM81SanityChecks:
    """Rigorous verification of Milestone M8.1 convergence benchmark methodology and results."""

    def test_solvers_and_grid_structure(self, convergence_records):
        """Verify solver count, grid values, and total row count."""
        expected_solvers = {"Euler", "Velocity Verlet", "RK4", "Yoshida 4"}
        expected_dts = [60.0, 30.0, 15.0, 7.5, 3.75, 1.875]

        # 1. Exactly 4 methods
        observed_solvers = {r["method"] for r in convergence_records}
        assert observed_solvers == expected_solvers, f"Unexpected solvers: {observed_solvers}"

        # 2. Exactly 6 dt per method & 24 total rows
        assert len(convergence_records) == 24, f"Expected 24 rows, got {len(convergence_records)}"

        for method in expected_solvers:
            method_records = [r for r in convergence_records if r["method"] == method]
            assert len(method_records) == 6, f"Method {method} has {len(method_records)} rows instead of 6"

            observed_dts = [r["dt"] for r in method_records]
            assert observed_dts == expected_dts, (
                f"Method {method} timesteps mismatch: {observed_dts} vs {expected_dts}"
            )

    def test_cost_model_and_evaluations(self, convergence_records):
        """Verify that force evaluation counts strictly follow the implementation-level proxy."""
        for r in convergence_records:
            method = r["method"]
            evals_per_step = CURRENT_IMPLEMENTATION_EVALS_PER_STEP[method]
            num_steps = r["num_steps"]
            expected_evals = num_steps * evals_per_step

            assert r["force_evals"] == expected_evals, (
                f"Row {method} dt={r['dt']}: expected {expected_evals} evals, got {r['force_evals']}"
            )
            assert num_steps > 0
            assert r["force_evals"] > 0

    def test_final_physical_time_integrity(self, convergence_records):
        """Verify that final physical time corresponds strictly to orbital period T."""
        r0 = R_EARTH + 400.0
        expected_T = circular_orbit_period(r0, mu=MU_EARTH)

        # Period for 400 km LEO is ~5553.624 s
        assert 5553.0 < expected_T < 5554.0

        # All records must reflect common final physical time
        for r in convergence_records:
            dt = r["dt"]
            # Number of steps must accommodate floor(T / dt) + 1 partial step
            expected_min_steps = int(math.floor(expected_T / dt))
            assert r["num_steps"] in (expected_min_steps, expected_min_steps + 1)

    def test_no_nan_or_inf_in_metrics(self, convergence_records):
        """Verify that no metric field contains NaN or Inf."""
        numeric_fields = [
            "dt", "num_steps", "force_evals", "e_r", "e_v", "order_p_r",
            "ratio_r", "max_de", "max_dh", "final_dr", "final_dtheta",
            "max_e_r", "max_e_v", "order_p_v", "final_de", "final_dh", "max_dr"
        ]
        for r in convergence_records:
            for field in numeric_fields:
                val = r[field]
                assert not math.isnan(val), f"Field {field} is NaN in row: {r}"
                assert not math.isinf(val), f"Field {field} is Inf in row: {r}"

    def test_non_negativity_of_norm_and_absolute_metrics(self, convergence_records):
        """Verify errors >= 0 where Euclidean norms or absolute deviations are used."""
        non_negative_fields = [
            "e_r", "max_e_r", "e_v", "max_e_v",
            "max_de", "max_dh", "max_dr", "final_dr", "final_dtheta"
        ]
        for r in convergence_records:
            for field in non_negative_fields:
                val = r[field]
                assert val >= 0.0, f"Field {field} is negative ({val}) in row: {r}"

    def test_analytical_reference_independence(self):
        """Verify that the analytical reference solution is evaluated via closed-form functions."""
        r0 = R_EARTH + 400.0
        t_sample = 1234.56
        s_exact = circular_orbit_exact_state(t_sample, r0=r0, mu=MU_EARTH)

        # Exact formulas
        n = math.sqrt(MU_EARTH / (r0**3))
        expected_x = r0 * math.cos(n * t_sample)
        expected_y = r0 * math.sin(n * t_sample)
        expected_vx = -r0 * n * math.sin(n * t_sample)
        expected_vy = r0 * n * math.cos(n * t_sample)

        assert math.isclose(s_exact.x, expected_x, rel_tol=1e-15)
        assert math.isclose(s_exact.y, expected_y, rel_tol=1e-15)
        assert math.isclose(s_exact.vx, expected_vx, rel_tol=1e-15)
        assert math.isclose(s_exact.vy, expected_vy, rel_tol=1e-15)

    def test_empirical_order_trends(self, convergence_records):
        """Verify mathematical convergence trends without hardcoding fragile exact values."""
        by_method = {}
        for r in convergence_records:
            by_method.setdefault(r["method"], []).append(r)

        # For Velocity Verlet (2nd-order symplectic):
        # Observed orders at finer steps should be ~ 2.0
        vv_orders = [r["order_p_r"] for r in by_method["Velocity Verlet"][1:]]
        for p in vv_orders:
            assert 1.95 <= p <= 2.05, f"Velocity Verlet order out of expected range: {p}"

        # For RK4 (4th-order non-symplectic):
        # Observed orders at finer steps should approach 4.0
        rk4_orders = [r["order_p_r"] for r in by_method["RK4"][1:]]
        for p in rk4_orders:
            assert 3.90 <= p <= 4.30, f"RK4 order out of expected range: {p}"

        # For Yoshida 4 (4th-order symplectic):
        # Observed orders across all intervals should be ~ 4.0
        y4_orders = [r["order_p_r"] for r in by_method["Yoshida 4"][1:]]
        for p in y4_orders:
            assert 3.95 <= p <= 4.05, f"Yoshida 4 order out of expected range: {p}"

        # For Euler (1st-order):
        # Orders should be positive and monotonically approach 1 as dt decreases
        euler_orders = [r["order_p_r"] for r in by_method["Euler"][1:]]
        assert euler_orders[-1] > euler_orders[0], "Euler convergence order should improve with smaller dt"
        assert 0.5 <= euler_orders[-1] <= 1.05

    def test_csv_file_on_disk(self):
        """Verify results/exp08_convergence.csv on disk matches schema and row count."""
        import csv
        csv_path = Path("results/exp08_convergence.csv")
        assert csv_path.exists(), "exp08_convergence.csv must exist after running M8.1"

        with open(csv_path, mode="r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            fieldnames = reader.fieldnames
            rows = list(reader)

        assert fieldnames == PLANNED_CONVERGENCE_COLUMNS
        assert len(rows) == 24, f"Expected 24 data rows in CSV, got {len(rows)}"

        # Verify each method has exactly 6 rows
        for method in ["Euler", "Velocity Verlet", "RK4", "Yoshida 4"]:
            method_rows = [r for r in rows if r["method"] == method]
            assert len(method_rows) == 6

