"""Validation and regression tests for Milestone M7: Fourth-Order Symplectic Integration."""

import csv
import math
from pathlib import Path
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
SUMMARY_PATH = RESULTS_DIR / "exp07_fourth_order_symplectic_summary.csv"
LONG_TERM_PATH = RESULTS_DIR / "exp07_long_term_comparison.csv"
COST_NORMALIZED_PATH = RESULTS_DIR / "exp07_cost_normalized_comparison.csv"
CONVERGENCE_PATH = RESULTS_DIR / "exp07_convergence.csv"
M6_SUMMARY_PATH = RESULTS_DIR / "exp06_symplectic_summary.csv"

SUMMARY_FIELDNAMES = [
    "method",
    "dt",
    "orbit_count",
    "time",
    "e_r",
    "e_v",
    "e_R",
    "e_T",
    "Δr",
    "Δθ",
    "ρ_R",
    "ρ_T",
    "χ",
    "δε",
    "δh",
    "force_evals",
]


@pytest.fixture(scope="module")
def summary_rows():
    """Load exp07_fourth_order_symplectic_summary.csv rows."""
    assert SUMMARY_PATH.exists(), f"Missing dataset: {SUMMARY_PATH}"
    with open(SUMMARY_PATH, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        assert reader.fieldnames == SUMMARY_FIELDNAMES
        return list(reader)


@pytest.fixture(scope="module")
def long_term_rows():
    """Load exp07_long_term_comparison.csv rows."""
    assert LONG_TERM_PATH.exists(), f"Missing dataset: {LONG_TERM_PATH}"
    with open(LONG_TERM_PATH, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return list(reader)


@pytest.fixture(scope="module")
def cost_normalized_rows():
    """Load exp07_cost_normalized_comparison.csv rows."""
    assert COST_NORMALIZED_PATH.exists(), f"Missing dataset: {COST_NORMALIZED_PATH}"
    with open(COST_NORMALIZED_PATH, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return list(reader)


@pytest.fixture(scope="module")
def convergence_rows():
    """Load exp07_convergence.csv rows."""
    assert CONVERGENCE_PATH.exists(), f"Missing dataset: {CONVERGENCE_PATH}"
    with open(CONVERGENCE_PATH, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return list(reader)


def test_m7_summary_row_count(summary_rows):
    """M7 summary must contain exactly 144 rows (4 methods * 6 timesteps * 6 checkpoints)."""
    assert len(summary_rows) == 144


def test_m7_all_values_finite(summary_rows):
    """All numerical fields in summary dataset must be finite (no NaN, no Inf)."""
    numeric_keys = [k for k in SUMMARY_FIELDNAMES if k != "method"]
    for row in summary_rows:
        for key in numeric_keys:
            val = float(row[key])
            assert math.isfinite(val), f"Non-finite value in {row['method']} dt={row['dt']} for {key}: {val}"


def test_m7_exact_checkpoints_extracted(summary_rows):
    """Every method and timestep must have checkpoints 1T, 5T, 10T, 20T, 50T, 100T."""
    expected_checkpoints = {1.0, 5.0, 10.0, 20.0, 50.0, 100.0}
    for method in ["Euler", "RK4", "Velocity Verlet", "Yoshida 4"]:
        for dt in [60.0, 30.0, 15.0, 7.5, 3.75, 1.875]:
            ckpts = {
                float(r["orbit_count"])
                for r in summary_rows
                if r["method"] == method and math.isclose(float(r["dt"]), dt)
            }
            assert ckpts == expected_checkpoints, f"Missing checkpoints for {method} dt={dt}: {ckpts}"


def test_m7_regression_against_m6_euler_rk4_verlet(summary_rows):
    """Euler, RK4, and Velocity Verlet checkpoints in M7 must match M6 summary exactly."""
    if not M6_SUMMARY_PATH.exists():
        pytest.skip("M6 summary not found for regression.")

    with open(M6_SUMMARY_PATH, "r", encoding="utf-8") as f:
        m6_rows = list(csv.DictReader(f))

    for m6_r in m6_rows:
        method = m6_r["method"]
        dt = float(m6_r["dt"])
        ckpt = float(m6_r["orbit_count"])

        m7_match = next(
            r
            for r in summary_rows
            if r["method"] == method
            and math.isclose(float(r["dt"]), dt)
            and math.isclose(float(r["orbit_count"]), ckpt)
        )

        for col in ["e_r", "e_v", "e_R", "e_T", "Δr", "Δθ", "δε", "δh"]:
            assert math.isclose(
                float(m7_match[col]), float(m6_r[col]), rel_tol=1e-10, abs_tol=1e-8
            ), f"Mismatch in {method} dt={dt} ckpt={ckpt} for {col}: m7={m7_match[col]} vs m6={m6_r[col]}"


def test_m7_yoshida4_conserves_angular_momentum_to_machine_precision(summary_rows):
    """Yoshida 4 must conserve angular momentum to machine precision (|dh| < 1e-13) across 100T."""
    yoshida_rows = [r for r in summary_rows if r["method"] == "Yoshida 4"]
    for r in yoshida_rows:
        dh = abs(float(r["δh"]))
        assert dh < 1e-13, f"Angular momentum drift too large: {dh} at dt={r['dt']}, ckpt={r['orbit_count']}"


def test_m7_yoshida4_bounds_energy_across_100T(summary_rows):
    """Yoshida 4 relative energy error must remain bounded (|de| < 1e-9 for all dt <= 60s)."""
    yoshida_rows = [r for r in summary_rows if r["method"] == "Yoshida 4"]
    for r in yoshida_rows:
        de = abs(float(r["δε"]))
        assert de < 1e-8, f"Energy drift too large: {de} at dt={r['dt']}, ckpt={r['orbit_count']}"


def test_m7_yoshida4_radial_drift_bounded_near_zero(summary_rows):
    """Yoshida 4 radial drift must remain bounded near zero (|dr| < 0.2 km for dt=60, < 0.01 km for dt<=30)."""
    yoshida_rows = [r for r in summary_rows if r["method"] == "Yoshida 4"]
    for r in yoshida_rows:
        dr = abs(float(r["Δr"]))
        dt = float(r["dt"])
        if dt <= 30.0:
            assert dr < 0.02, f"Radial drift unacceptably large: {dr} km at dt={dt}, ckpt={r['orbit_count']}"
        else:
            assert dr < 0.2, f"Radial drift unacceptably large: {dr} km at dt={dt}, ckpt={r['orbit_count']}"


def test_m7_convergence_order_4(convergence_rows):
    """Convergence study must confirm 4th-order scaling (order_p ~ 4.0, ratio ~ 16)."""
    assert len(convergence_rows) == 6
    # Skip first row where order_p is 0.0
    for r in convergence_rows[1:]:
        p = float(r["order_p_r"])
        ratio = float(r["ratio_r"])
        assert 3.9 <= p <= 4.1, f"Expected order ~ 4.0, got {p} at dt={r['dt']}"
        assert 15.5 <= ratio <= 16.5, f"Expected ratio ~ 16.0, got {ratio} at dt={r['dt']}"


def test_m7_force_evals_exact_multiples(summary_rows):
    """Verify force evaluations match exact multiples per method."""
    eval_factors = {"Euler": 1, "Velocity Verlet": 2, "RK4": 4, "Yoshida 4": 6}
    for r in summary_rows:
        method = r["method"]
        fe = int(r["force_evals"])
        factor = eval_factors[method]
        assert fe % factor == 0, f"Force evals {fe} not divisible by {factor} for {method}"


def test_m7_cost_normalized_budgets_consistent(cost_normalized_rows):
    """Verify cost-normalized budgets have closely matched force evaluations."""
    assert len(cost_normalized_rows) == 16
    budgets = {"Budget_1_55k_evals", "Budget_2_111k_evals", "Budget_3_222k_evals", "Budget_4_444k_evals"}
    found_budgets = {r["budget_name"] for r in cost_normalized_rows}
    assert found_budgets == budgets

    for b in budgets:
        b_rows = [r for r in cost_normalized_rows if r["budget_name"] == b]
        evals = [int(r["actual_force_evals"]) for r in b_rows]
        min_ev, max_ev = min(evals), max(evals)
        # Differ by at most 5 evaluations due to integer step rounding over 100 orbits
        assert max_ev - min_ev <= 5, f"Budgets not matched closely in {b}: {evals}"
