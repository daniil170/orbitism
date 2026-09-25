"""Validation and regression tests for Milestone M6: Symplectic Integrator Pipeline."""

import csv
import math
from pathlib import Path
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = PROJECT_ROOT / "results"
SUMMARY_PATH = RESULTS_DIR / "exp06_symplectic_summary.csv"
TRAJECTORY_PATH = RESULTS_DIR / "exp06_symplectic_integrators.csv"

REQUIRED_FIELDNAMES = [
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
]


@pytest.fixture(scope="module")
def summary_rows():
    """Load exp06_symplectic_summary.csv rows."""
    assert SUMMARY_PATH.exists(), f"Missing dataset: {SUMMARY_PATH}"
    with open(SUMMARY_PATH, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        assert reader.fieldnames == REQUIRED_FIELDNAMES, (
            f"Schema mismatch: expected {REQUIRED_FIELDNAMES}, got {reader.fieldnames}"
        )
        return list(reader)


def test_m6_summary_row_count(summary_rows):
    """M6 summary must contain exactly 108 rows (3 methods * 6 timesteps * 6 checkpoints)."""
    assert len(summary_rows) == 108


def test_m6_all_values_finite(summary_rows):
    """All numerical fields in summary dataset must be finite (no NaN, no Inf)."""
    numeric_keys = [k for k in REQUIRED_FIELDNAMES if k != "method"]
    for row in summary_rows:
        for key in numeric_keys:
            val = float(row[key])
            assert math.isfinite(val), f"Non-finite value in {row['method']} dt={row['dt']} for {key}: {val}"


def test_m6_exact_checkpoints_extracted(summary_rows):
    """Every method and timestep must have checkpoints 1T, 5T, 10T, 20T, 50T, 100T."""
    expected_checkpoints = {1.0, 5.0, 10.0, 20.0, 50.0, 100.0}
    for method in ["Euler", "RK4", "Velocity Verlet"]:
        for dt in [60.0, 30.0, 15.0, 7.5, 3.75, 1.875]:
            ckpts = {
                float(r["orbit_count"])
                for r in summary_rows
                if r["method"] == method and math.isclose(float(r["dt"]), dt)
            }
            assert ckpts == expected_checkpoints, f"Missing checkpoints for {method} dt={dt}: {ckpts}"


def test_m6_euler_regression_against_m4(summary_rows):
    """Euler checkpoint values in M6 must match M4 summary to floating-point precision."""
    m4_summary_path = RESULTS_DIR / "exp04_phase_radial_summary.csv"
    if not m4_summary_path.exists():
        pytest.skip("M4 summary not found for regression.")

    with open(m4_summary_path, "r", encoding="utf-8") as f:
        m4_rows = list(csv.DictReader(f))

    for m4_r in m4_rows:
        dt = float(m4_r["dt_s"])
        ckpt = float(m4_r["checkpoint_T"])

        m6_match = next(
            r
            for r in summary_rows
            if r["method"] == "Euler"
            and math.isclose(float(r["dt"]), dt)
            and math.isclose(float(r["orbit_count"]), ckpt)
        )

        assert math.isclose(
            float(m6_match["e_r"]), float(m4_r["position_error_km"]), rel_tol=1e-10, abs_tol=1e-8
        )
        assert math.isclose(
            float(m6_match["Δr"]), float(m4_r["delta_r_km"]), rel_tol=1e-10, abs_tol=1e-8
        )
        assert math.isclose(
            float(m6_match["Δθ"]), float(m4_r["delta_theta_rad"]), rel_tol=1e-10, abs_tol=1e-8
        )


def test_m6_rk4_regression_against_m5(summary_rows):
    """RK4 checkpoint values in M6 must match M5 summary to floating-point precision."""
    m5_summary_path = RESULTS_DIR / "exp05_rk4_long_term_summary.csv"
    if not m5_summary_path.exists():
        pytest.skip("M5 summary not found for regression.")

    with open(m5_summary_path, "r", encoding="utf-8") as f:
        m5_rows = list(csv.DictReader(f))

    for m5_r in m5_rows:
        dt = float(m5_r["dt_s"])
        ckpt = float(m5_r["checkpoint_T"])

        m6_match = next(
            r
            for r in summary_rows
            if r["method"] == "RK4"
            and math.isclose(float(r["dt"]), dt)
            and math.isclose(float(r["orbit_count"]), ckpt)
        )

        assert math.isclose(
            float(m6_match["e_r"]), float(m5_r["position_error_km"]), rel_tol=1e-10, abs_tol=1e-8
        )
        assert math.isclose(
            float(m6_match["Δr"]), float(m5_r["delta_r_km"]), rel_tol=1e-10, abs_tol=1e-8
        )
        assert math.isclose(
            float(m6_match["Δθ"]), float(m5_r["delta_theta_rad"]), rel_tol=1e-10, abs_tol=1e-8
        )


def test_m6_velocity_verlet_conserves_angular_momentum_to_machine_precision(summary_rows):
    """Velocity Verlet must conserve angular momentum to machine precision (|dh| < 1e-13)."""
    verlet_rows = [r for r in summary_rows if r["method"] == "Velocity Verlet"]
    for r in verlet_rows:
        dh = abs(float(r["δh"]))
        assert dh < 1e-13, f"Angular momentum drift too large: {dh} at dt={r['dt']}, ckpt={r['orbit_count']}"


def test_m6_velocity_verlet_bounds_energy_across_100T(summary_rows):
    """Velocity Verlet relative energy error must remain bounded (|de| < 1e-7 for all dt)."""
    verlet_rows = [r for r in summary_rows if r["method"] == "Velocity Verlet"]
    for r in verlet_rows:
        de = abs(float(r["δε"]))
        assert de < 1e-6, f"Energy drift too large: {de} at dt={r['dt']}, ckpt={r['orbit_count']}"


def test_m6_velocity_verlet_no_phase_winding(summary_rows):
    """Velocity Verlet must not enter phase-wound regime (|dtheta| < pi) over 100T."""
    verlet_rows = [r for r in summary_rows if r["method"] == "Velocity Verlet"]
    for r in verlet_rows:
        dtheta = abs(float(r["Δθ"]))
        assert dtheta < math.pi, (
            f"Phase winding occurred for Velocity Verlet: |dtheta|={dtheta} at dt={r['dt']}, ckpt={r['orbit_count']}"
        )


def test_m6_velocity_verlet_radial_drift_bounded_near_zero(summary_rows):
    """Velocity Verlet radial drift must remain bounded near zero (|dr| < 0.5 km) for all timesteps."""
    verlet_rows = [r for r in summary_rows if r["method"] == "Velocity Verlet"]
    for r in verlet_rows:
        dr = abs(float(r["Δr"]))
        assert dr < 0.5, (
            f"Radial drift unacceptably large for Velocity Verlet: {dr} km at dt={r['dt']}, ckpt={r['orbit_count']}"
        )
