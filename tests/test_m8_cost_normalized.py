"""Tests for Milestone M8.2b: Cost-Normalized Solver Comparison (100 Orbits).

Verifies:
1. Solvers and budget grid completeness (4 methods x 3 budgets = 12 runs).
2. Strict budget equality: force_evals_actual == requested_budget.
3. Step accounting: num_steps == budget // evals_per_step and dt == final_time / num_steps.
4. Exact final time reaching 100*T within < 1e-11 s.
5. Absence of NaN/Inf and valid non-negative error norms.
6. CSV file schema and row count on disk.
7. Physical qualitative trade-offs without ranking or "winner" assertions.
"""

import csv
import math
from pathlib import Path
import pytest

from src.physics.gravity import MU_EARTH, R_EARTH
from src.physics.analytical import circular_orbit_period
from experiments.exp08_systematic_comparison import (
    CURRENT_IMPLEMENTATION_EVALS_PER_STEP,
    DEFAULT_ALTITUDE_KM,
    DEFAULT_COST_NORMALIZED_BUDGETS,
    DEFAULT_NUM_ORBITS,
    PLANNED_COST_NORMALIZED_COLUMNS,
    SOLVER_REGISTRY,
    run_cost_normalized_study,
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def cost_normalized_data():
    """Run or load M8.2b cost-normalized study once for module tests."""
    csv_path = PROJECT_ROOT / "results" / "exp08_cost_normalized.csv"
    if csv_path.exists():
        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            records = []
            for row in reader:
                parsed = {
                    "method": row["method"],
                    "force_eval_budget": int(row["force_eval_budget"]),
                    "force_evals_actual": int(row["force_evals_actual"]),
                    "num_steps": int(row["num_steps"]),
                    "dt": float(row["dt"]),
                    "final_time": float(row["final_time"]),
                    "final_position_error_km": float(row["final_position_error_km"]),
                    "final_velocity_error_km_s": float(row["final_velocity_error_km_s"]),
                    "max_position_error_km": float(row["max_position_error_km"]),
                    "max_velocity_error_km_s": float(row["max_velocity_error_km_s"]),
                    "final_radial_error_km": float(row["final_radial_error_km"]),
                    "max_radial_error_km": float(row["max_radial_error_km"]),
                    "final_phase_error_rad": float(row["final_phase_error_rad"]),
                    "max_abs_phase_error_rad": float(row["max_abs_phase_error_rad"]),
                    "final_relative_energy_error": float(row["final_relative_energy_error"]),
                    "max_abs_relative_energy_error": float(row["max_abs_relative_energy_error"]),
                    "final_relative_angular_momentum_error": float(row["final_relative_angular_momentum_error"]),
                    "max_abs_relative_angular_momentum_error": float(row["max_abs_relative_angular_momentum_error"]),
                }
                records.append(parsed)
            return records
    return run_cost_normalized_study()


class TestM8CostNormalized:
    """Test suite for Milestone M8.2b cost-normalized comparison."""

    def test_solvers_and_budgets_presence(self, cost_normalized_data):
        """Verify all 4 solvers and 3 budgets are evaluated (12 total rows)."""
        records = cost_normalized_data
        assert len(records) == 12

        observed_methods = {r["method"] for r in records}
        expected_methods = set(SOLVER_REGISTRY.keys())
        assert observed_methods == expected_methods

        observed_budgets = {r["force_eval_budget"] for r in records}
        expected_budgets = set(DEFAULT_COST_NORMALIZED_BUDGETS)
        assert observed_budgets == expected_budgets

        # Each method must have exactly 3 runs
        for method in expected_methods:
            method_rows = [r for r in records if r["method"] == method]
            assert len(method_rows) == 3
            assert {r["force_eval_budget"] for r in method_rows} == expected_budgets

    def test_strict_budget_equality_and_steps(self, cost_normalized_data):
        """Verify force_evals_actual == requested_budget and correct step formula."""
        r0 = R_EARTH + DEFAULT_ALTITUDE_KM
        t_orbit = circular_orbit_period(r0, mu=MU_EARTH)
        t_final = DEFAULT_NUM_ORBITS * t_orbit

        for r in cost_normalized_data:
            method = r["method"]
            budget = r["force_eval_budget"]
            evals_per_step = CURRENT_IMPLEMENTATION_EVALS_PER_STEP[method]

            # Budget must be exactly divisible by evals_per_step
            assert budget % evals_per_step == 0

            # Step count must equal budget // evals_per_step
            expected_steps = budget // evals_per_step
            assert r["num_steps"] == expected_steps

            # Actual force evaluations must match requested budget exactly
            assert r["force_evals_actual"] == budget
            assert r["force_evals_actual"] == r["num_steps"] * evals_per_step

            # Timestep must equal final_time / num_steps
            expected_dt = t_final / expected_steps
            assert math.isclose(r["dt"], expected_dt, rel_tol=1e-12)

    def test_final_time_precision(self, cost_normalized_data):
        """Verify all runs terminate exactly at 100*T with error < 1e-11 s."""
        r0 = R_EARTH + DEFAULT_ALTITUDE_KM
        t_orbit = circular_orbit_period(r0, mu=MU_EARTH)
        t_final = DEFAULT_NUM_ORBITS * t_orbit

        for r in cost_normalized_data:
            assert abs(r["final_time"] - t_final) < 1e-11

    def test_numeric_sanity_and_validity(self, cost_normalized_data):
        """Verify no NaN/Inf, non-negativity of norms, and max >= final."""
        for r in cost_normalized_data:
            for k, v in r.items():
                if isinstance(v, (int, float)):
                    assert not math.isnan(v), f"NaN in {k}"
                    assert not math.isinf(v), f"Inf in {k}"

            assert r["max_position_error_km"] >= 0.0
            assert r["final_position_error_km"] >= 0.0
            assert r["max_velocity_error_km_s"] >= 0.0
            assert r["final_velocity_error_km_s"] >= 0.0
            assert r["max_radial_error_km"] >= 0.0
            assert r["final_radial_error_km"] >= 0.0
            assert r["max_abs_relative_energy_error"] >= 0.0
            assert r["max_abs_relative_angular_momentum_error"] >= 0.0
            assert r["max_abs_phase_error_rad"] >= 0.0

            # Maxima must be greater than or equal to final values
            assert r["max_position_error_km"] >= r["final_position_error_km"] - 1e-12
            assert r["max_velocity_error_km_s"] >= r["final_velocity_error_km_s"] - 1e-12
            assert r["max_radial_error_km"] >= r["final_radial_error_km"] - 1e-12
            assert r["max_abs_relative_energy_error"] >= abs(r["final_relative_energy_error"]) - 1e-15
            assert (
                r["max_abs_relative_angular_momentum_error"]
                >= abs(r["final_relative_angular_momentum_error"]) - 1e-15
            )
            assert r["max_abs_phase_error_rad"] >= abs(r["final_phase_error_rad"]) - 1e-12

    def test_csv_file_on_disk(self):
        """Verify that exp08_cost_normalized.csv exists and has expected schema."""
        csv_path = PROJECT_ROOT / "results" / "exp08_cost_normalized.csv"
        if not csv_path.exists():
            pytest.skip("exp08_cost_normalized.csv not generated yet")

        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            assert reader.fieldnames == PLANNED_COST_NORMALIZED_COLUMNS
            rows = list(reader)
            assert len(rows) == 12

    def test_qualitative_tradeoffs_without_ranking(self, cost_normalized_data):
        """Verify characteristic solver behaviors across matched budgets."""
        by_method = {}
        for r in cost_normalized_data:
            by_method.setdefault(r["method"], []).append(r)

        # 1. Euler diverges across all budgets
        for r in by_method["Euler"]:
            assert r["max_position_error_km"] > 5000.0
            assert r["max_radial_error_km"] > 5000.0
            assert r["max_abs_relative_energy_error"] > 0.05

        # 2. Velocity Verlet preserves angular momentum to machine precision
        for r in by_method["Velocity Verlet"]:
            assert r["max_abs_relative_angular_momentum_error"] < 1e-13
            assert r["max_abs_relative_energy_error"] < 1e-4

        # 3. RK4 exhibits secular dissipation in energy and angular momentum
        for r in by_method["RK4"]:
            assert r["final_relative_energy_error"] < 0.0
            assert r["final_relative_angular_momentum_error"] < 0.0

        # 4. Yoshida 4 preserves angular momentum to machine precision and bounds energy
        for r in by_method["Yoshida 4"]:
            assert r["max_abs_relative_angular_momentum_error"] < 1e-13
            assert r["max_abs_relative_energy_error"] < 1e-6
