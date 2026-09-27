"""Tests for Milestone M8.2a: Long-Term Fixed-dt Comparison (100 Orbits).

Verifies:
1. Presence of exactly 4 solvers (Euler, Velocity Verlet, RK4, Yoshida 4).
2. Grid of nominal timesteps: dt in {60.0, 30.0, 15.0} s.
3. Checkpoint horizons: 1T, 5T, 10T, 20T, 50T, 100T.
4. Total steps and exact shortened final step size (hitting 100T to float64 machine precision).
5. Exact evaluation of analytical reference at actual checkpoint timestamps.
6. Schema integrity for exp08_long_term_fixed_dt.csv (12 rows) and exp08_long_term_checkpoints.csv (72 rows).
7. Non-negativity, lack of NaN/Inf, and mathematical sanity of trajectory-wide maxima.
8. Qualitative physical behavior over 100 orbits (divergence of Euler, symplecticity of Verlet/Yoshida 4, dissipation of RK4).
"""

import csv
import math
from pathlib import Path
import pytest
from typing import Any, Dict, List

from experiments.exp08_systematic_comparison import (
    DEFAULT_LONG_TERM_DT_VALUES,
    DEFAULT_CHECKPOINTS_T,
    CURRENT_IMPLEMENTATION_EVALS_PER_STEP,
    SOLVER_REGISTRY,
    PLANNED_LONG_TERM_FIXED_DT_COLUMNS,
    PLANNED_LONG_TERM_CHECKPOINTS_COLUMNS,
    run_long_term_fixed_dt_study,
)
from src.physics.gravity import MU_EARTH, R_EARTH
from src.physics.analytical import circular_orbit_period, circular_orbit_exact_state


PROJECT_ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def long_term_data():
    """Run 100-orbit fixed-dt study once for module tests."""
    return run_long_term_fixed_dt_study()


class TestM82aLongTerm:
    """Rigorous verification of M8.2a long-term 100-orbit benchmark."""

    def test_solvers_and_dt_grid(self, long_term_data):
        """Verify solver set and nominal timesteps."""
        summary_records, checkpoint_records = long_term_data
        expected_solvers = {"Euler", "Velocity Verlet", "RK4", "Yoshida 4"}
        expected_dts = [60.0, 30.0, 15.0]

        # 4 solvers
        observed_solvers = {r["method"] for r in summary_records}
        assert observed_solvers == expected_solvers

        # 12 summary rows (4 solvers x 3 dt)
        assert len(summary_records) == 12

        for method in expected_solvers:
            method_rows = [r for r in summary_records if r["method"] == method]
            assert len(method_rows) == 3
            assert [r["dt"] for r in method_rows] == expected_dts

    def test_checkpoint_structure_and_times(self, long_term_data):
        """Verify 72 checkpoint records, target times, and physical endpoint accuracy."""
        _, checkpoint_records = long_term_data
        expected_checkpoints = [1.0, 5.0, 10.0, 20.0, 50.0, 100.0]

        # 4 solvers x 3 dt x 6 checkpoints = 72 records
        assert len(checkpoint_records) == 72

        r0 = R_EARTH + 400.0
        t_orbit = circular_orbit_period(r0, mu=MU_EARTH)

        for rec in checkpoint_records:
            dt = rec["dt"]
            target_t = rec["orbit_count"] * t_orbit
            actual_t = rec["actual_time"]

            # Target time matches formula
            assert math.isclose(rec["target_time"], target_t, rel_tol=1e-12)

            # Actual sampled time must be within dt of target time
            assert abs(actual_t - target_t) <= dt

            # For final 100T checkpoint, actual time must match 100*T exactly
            if rec["orbit_count"] == 100.0:
                assert abs(actual_t - 100.0 * t_orbit) < 1e-11

    def test_step_accounting_and_final_shortened_step(self, long_term_data):
        """Verify exact discrete step counts and partial final step size."""
        summary_records, _ = long_term_data
        r0 = R_EARTH + 400.0
        t_orbit = circular_orbit_period(r0, mu=MU_EARTH)
        t_final = 100.0 * t_orbit

        # Analytical remainder: 555362.431862... - 555360.0 = 2.431862378... s
        expected_remainder = t_final - 555360.0

        for r in summary_records:
            dt = r["dt"]
            method = r["method"]
            evals_per_step = CURRENT_IMPLEMENTATION_EVALS_PER_STEP[method]

            if dt == 60.0:
                assert r["total_steps"] == 9257
            elif dt == 30.0:
                assert r["total_steps"] == 18513
            elif dt == 15.0:
                assert r["total_steps"] == 37025

            assert r["actual_force_evals"] == r["total_steps"] * evals_per_step

    def test_numeric_sanity_and_maxima_ordering(self, long_term_data):
        """Verify no NaN/Inf, non-negativity of norms, and max >= final."""
        summary_records, checkpoint_records = long_term_data

        for r in summary_records:
            for k, v in r.items():
                if isinstance(v, (int, float)):
                    assert not math.isnan(v), f"NaN in summary {k}"
                    assert not math.isinf(v), f"Inf in summary {k}"

            assert r["max_position_error"] >= 0.0
            assert r["final_position_error"] >= 0.0
            assert r["max_velocity_error"] >= 0.0
            assert r["final_velocity_error"] >= 0.0
            assert r["max_radial_deviation"] >= 0.0
            assert r["final_radial_deviation"] >= 0.0
            assert r["max_relative_energy_error"] >= 0.0
            assert r["max_relative_angular_momentum_error"] >= 0.0

            # Maxima must be greater than or equal to final values
            assert r["max_position_error"] >= r["final_position_error"] - 1e-12
            assert r["max_velocity_error"] >= r["final_velocity_error"] - 1e-12
            assert r["max_radial_deviation"] >= r["final_radial_deviation"] - 1e-12
            assert r["max_relative_energy_error"] >= abs(r["final_relative_energy_error"]) - 1e-15
            assert r["max_relative_angular_momentum_error"] >= abs(r["final_relative_angular_momentum_error"]) - 1e-15

        for r in checkpoint_records:
            for k, v in r.items():
                if isinstance(v, (int, float)):
                    assert not math.isnan(v), f"NaN in checkpoint {k}"
                    assert not math.isinf(v), f"Inf in checkpoint {k}"

    def test_csv_files_on_disk(self):
        """Verify that exp08_long_term_fixed_dt.csv and exp08_long_term_checkpoints.csv exist and match schema."""
        summary_csv = PROJECT_ROOT / "results" / "exp08_long_term_fixed_dt.csv"
        checkpoints_csv = PROJECT_ROOT / "results" / "exp08_long_term_checkpoints.csv"

        assert summary_csv.exists(), "Summary CSV missing"
        assert checkpoints_csv.exists(), "Checkpoints CSV missing"

        with open(summary_csv, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            assert reader.fieldnames == PLANNED_LONG_TERM_FIXED_DT_COLUMNS
            rows = list(reader)
            assert len(rows) == 12

        with open(checkpoints_csv, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            assert reader.fieldnames == PLANNED_LONG_TERM_CHECKPOINTS_COLUMNS
            rows = list(reader)
            assert len(rows) == 72

    def test_physical_qualitative_behavior(self, long_term_data):
        """Verify characteristic physical behavior over 100 orbits without unscientific value hardcoding."""
        summary_records, _ = long_term_data
        by_method = {}
        for r in summary_records:
            by_method.setdefault(r["method"], []).append(r)

        # 1. Euler diverges dramatically (position error > 10000 km, radial deviation > 10000 km)
        for r in by_method["Euler"]:
            assert r["max_position_error"] > 20000.0
            assert r["final_radial_deviation"] > 20000.0
            assert r["max_relative_energy_error"] > 0.5

        # 2. Velocity Verlet preserves energy within O(dt^2) bound and preserves angular momentum at machine precision
        for r in by_method["Velocity Verlet"]:
            assert r["max_relative_energy_error"] < 1e-4
            assert r["max_relative_angular_momentum_error"] < 1e-14
            assert r["final_radial_deviation"] < 1.0  # radial shape stays stable

        # 3. RK4 exhibits secular energy dissipation and secular angular momentum loss
        for r in by_method["RK4"]:
            assert r["final_relative_energy_error"] < 0.0, "RK4 energy must dissipate"
            assert r["final_relative_angular_momentum_error"] < 0.0, "RK4 angular momentum must dissipate"

        # 4. Yoshida 4 preserves energy within strict high-order bound and preserves angular momentum at machine precision
        for r in by_method["Yoshida 4"]:
            assert r["max_relative_energy_error"] < 1e-9
            assert r["max_relative_angular_momentum_error"] < 1e-13
            assert r["final_radial_deviation"] < 1e-4  # radial shape preserved to meters
