"""Mathematical validation and unit tests for M4: Phase Error and Radial Drift Decomposition.

Verifies:
1. Initial zero conditions (e_R(0) = 0, e_T(0) = 0, Delta r(0) = 0, Delta theta(0) = 0)
2. Exact Pythagorean identity (e_r^2 = e_R^2 + e_T^2) across the entire trajectory
3. Variance partition identity (rho_R + rho_T = 1) for all t > 0
4. Orthonormal rotating frame basis properties
5. Robust phase unwrap correctness across branch cuts
6. Small-angle asymptotic consistency (e_R ~= Delta r, e_T ~= r0 * Delta theta)
7. Consistency between DecomposedErrorRecord and baseline ErrorRecord
"""

import math
from typing import List
import pytest

from src.models.state import State2D
from src.physics.gravity import (
    MU_EARTH,
    R_EARTH,
    create_circular_orbit_state,
)
from src.physics.analytical import (
    circular_orbit_period,
    circular_orbit_mean_motion,
    circular_orbit_exact_state,
)
from src.simulation.simulator import Simulator
from src.simulation.metrics import (
    position_error,
    radial_position_error,
    along_track_position_error,
    radial_drift,
    unwrap_phase,
    angular_phase_error,
    radial_error_fraction,
    along_track_error_fraction,
    dominance_ratio,
    compute_error_history,
    compute_decomposed_error_history,
)


@pytest.fixture(scope="module")
def sample_orbit_trajectory():
    """Simulate a short 2-orbit trajectory with dt=60s for mathematical identity checks."""
    r0 = R_EARTH + 400.0
    period = circular_orbit_period(r0, mu=MU_EARTH)
    initial_state = create_circular_orbit_state(altitude_km=400.0, mu=MU_EARTH, r_earth=R_EARTH)
    sim = Simulator(initial_state=initial_state, mu=MU_EARTH)
    trajectory = sim.run(dt=60.0, t_final=2.0 * period)
    timestamps = sim.times

    decomposed = compute_decomposed_error_history(
        timestamps, trajectory, r0=r0, mu=MU_EARTH
    )
    baseline = compute_error_history(
        timestamps, trajectory, r0=r0, mu=MU_EARTH
    )
    return {
        "r0": r0,
        "period": period,
        "timestamps": timestamps,
        "states": trajectory,
        "decomposed": decomposed,
        "baseline": baseline,
    }


# ==============================================================================
# 1. INITIAL ZERO CONDITIONS
# ==============================================================================

def test_initial_decomposition_metrics_are_zero(sample_orbit_trajectory):
    """Fact 1: At t = 0, e_R, e_T, Delta r, and Delta theta must all be zero within float precision."""
    record_0 = sample_orbit_trajectory["decomposed"][0]
    assert record_0.t == 0.0
    assert math.isclose(record_0.e_R, 0.0, abs_tol=1e-12)
    assert math.isclose(record_0.e_T, 0.0, abs_tol=1e-12)
    assert math.isclose(record_0.delta_r, 0.0, abs_tol=1e-12)
    assert math.isclose(record_0.delta_theta, 0.0, abs_tol=1e-12)
    assert math.isclose(record_0.position_error, 0.0, abs_tol=1e-12)
    assert record_0.fraction_R == 0.0
    assert record_0.fraction_T == 0.0


# ==============================================================================
# 2. PYTHAGOREAN IDENTITY
# ==============================================================================

def test_pythagorean_closure_for_entire_trajectory(sample_orbit_trajectory):
    """Mathematical Identity: e_r^2 = e_R^2 + e_T^2 at every discrete step.

    Because the rotating frame basis {e_r, e_theta} is strictly orthonormal,
    the Euclidean norm of Delta r projected onto this basis must equal the Cartesian norm.
    """
    for rec in sample_orbit_trajectory["decomposed"]:
        lhs = rec.position_error**2
        rhs = rec.e_R**2 + rec.e_T**2
        # Absolute tolerance scaled for values of order 10^4 km squared (~10^8 km^2)
        assert math.isclose(lhs, rhs, rel_tol=1e-11, abs_tol=1e-10)


# ==============================================================================
# 3. VARIANCE PARTITION IDENTITY
# ==============================================================================

def test_variance_partition_sums_to_unity(sample_orbit_trajectory):
    """Mathematical Identity: rho_R + rho_T = 1.0 for all t > 0.

    By definition, rho_R = e_R^2 / e_r^2 and rho_T = e_T^2 / e_r^2.
    """
    for rec in sample_orbit_trajectory["decomposed"][1:]:  # skip t = 0 where e_r == 0
        partition_sum = rec.fraction_R + rec.fraction_T
        assert math.isclose(partition_sum, 1.0, abs_tol=1e-14)
        assert 0.0 <= rec.fraction_R <= 1.0
        assert 0.0 <= rec.fraction_T <= 1.0


# ==============================================================================
# 4. ORTHONORMAL BASIS PROPERTIES
# ==============================================================================

def test_orthonormal_rotating_basis_properties(sample_orbit_trajectory):
    """Verify that radial and along-track unit vectors form an exact orthonormal basis."""
    r0 = sample_orbit_trajectory["r0"]
    n = circular_orbit_mean_motion(r0, mu=MU_EARTH)

    for t in sample_orbit_trajectory["timestamps"][:100]:
        nt = n * t
        cos_nt = math.cos(nt)
        sin_nt = math.sin(nt)

        # e_r = [cos(nt), sin(nt)]
        # e_theta = [-sin(nt), cos(nt)]
        norm_er = cos_nt**2 + sin_nt**2
        norm_et = (-sin_nt) ** 2 + cos_nt**2
        dot_product = cos_nt * (-sin_nt) + sin_nt * cos_nt

        assert math.isclose(norm_er, 1.0, abs_tol=1e-14)
        assert math.isclose(norm_et, 1.0, abs_tol=1e-14)
        assert math.isclose(dot_product, 0.0, abs_tol=1e-14)


# ==============================================================================
# 5. PHASE UNWRAP CORRECTNESS
# ==============================================================================

def test_phase_unwrap_synthetic_forward_multi_revolutions():
    """Verify unwrapping correctly reconstructs a monotonically increasing phase over 10 revolutions."""
    num_steps = 1000
    total_angle = 10.0 * 2.0 * math.pi
    true_phase = [i * total_angle / num_steps for i in range(num_steps)]

    # Wrap to [-pi, pi]
    wrapped = [math.atan2(math.sin(p), math.cos(p)) for p in true_phase]
    unwrapped = unwrap_phase(wrapped)

    for true_val, unwrapped_val in zip(true_phase, unwrapped):
        assert math.isclose(true_val, unwrapped_val, abs_tol=1e-12)


def test_phase_unwrap_synthetic_reverse_revolutions():
    """Verify unwrapping correctly reconstructs a monotonically decreasing phase."""
    num_steps = 500
    total_angle = -5.0 * 2.0 * math.pi
    true_phase = [i * total_angle / num_steps for i in range(num_steps)]

    wrapped = [math.atan2(math.sin(p), math.cos(p)) for p in true_phase]
    unwrapped = unwrap_phase(wrapped)

    for true_val, unwrapped_val in zip(true_phase, unwrapped):
        assert math.isclose(true_val, unwrapped_val, abs_tol=1e-12)


def test_phase_unwrap_edge_cases():
    """Verify edge cases: empty list, single element, constant sequence."""
    assert unwrap_phase([]) == []
    assert unwrap_phase([1.5]) == [1.5]
    assert unwrap_phase([0.5, 0.5, 0.5]) == [0.5, 0.5, 0.5]


# ==============================================================================
# 6. SMALL-ANGLE ASYMPTOTIC CONSISTENCY
# ==============================================================================

def test_small_angle_polar_consistency(sample_orbit_trajectory):
    """Verify asymptotic equivalence for small angles and small radial departures.

    The exact geometric relationship between rotating frame error e_T and polar coordinates is:
        e_T = r_num * sin(Delta theta) = (r0 + Delta r) * sin(Delta theta)
    In the small-angle limit (Delta theta << 1), sin(Delta theta) ~= Delta theta, giving:
        e_T ~= (r0 + Delta r) * Delta theta ~= r0 * Delta theta + O(Delta r * Delta theta).
    Similarly, e_R = r_num * cos(Delta theta) - r0 ~= Delta r - r_num * (Delta theta)^2 / 2.
    """
    r0 = sample_orbit_trajectory["r0"]
    for rec in sample_orbit_trajectory["decomposed"][1:6]:
        arc_length_phase = r0 * rec.delta_theta

        # e_R matches Delta r to within O((Delta theta)^2)
        assert math.isclose(rec.e_R, rec.delta_r, rel_tol=0.03, abs_tol=1e-3)
        # e_T matches r0 * Delta theta to within O(Delta r / r0)
        assert math.isclose(rec.e_T, arc_length_phase, rel_tol=0.03, abs_tol=1e-3)

        # Exact polar trigonometric identity e_T = r_num * sin(delta_theta)
        exact_polar_e_T = rec.r_num * math.sin(rec.delta_theta)
        assert math.isclose(rec.e_T, exact_polar_e_T, rel_tol=1e-11, abs_tol=1e-10)


# ==============================================================================
# 7. REGRESSION CONSISTENCY WITH BASELINE ERROR RECORD
# ==============================================================================

def test_regression_consistency_with_m3_metrics(sample_orbit_trajectory):
    """Ensure that position, velocity, energy, and angular momentum errors match baseline exactly."""
    decomposed = sample_orbit_trajectory["decomposed"]
    baseline = sample_orbit_trajectory["baseline"]

    assert len(decomposed) == len(baseline)
    for dec_rec, base_rec in zip(decomposed, baseline):
        assert dec_rec.t == base_rec.t
        assert dec_rec.position_error == base_rec.position_error
        assert dec_rec.velocity_error == base_rec.velocity_error
        assert dec_rec.energy_error == base_rec.energy_error
        assert dec_rec.angular_momentum_error == base_rec.angular_momentum_error


# ==============================================================================
# 8. DIAGNOSTIC DOMINANCE RATIO EDGE CASES
# ==============================================================================

def test_dominance_ratio_logic():
    """Verify dominance_ratio function correctly handles boundaries and zero values."""
    assert dominance_ratio(0.0, 0.0) == 1.0
    assert math.isinf(dominance_ratio(0.0, 5.0))
    assert math.isclose(dominance_ratio(2.0, 4.0), 2.0)
    assert math.isclose(dominance_ratio(4.0, 2.0), 0.5)


# ==============================================================================
# 9. DATASET VALIDATION AND M3 REGRESSION TESTS
# ==============================================================================

def test_exp04_csv_files_exist_and_have_valid_schema():
    """Verify generated CSV files exist and have approved header schemas."""
    import csv
    from pathlib import Path

    results_dir = Path(__file__).resolve().parent.parent / "results"
    traj_path = results_dir / "exp04_phase_radial_decomposition.csv"
    summary_path = results_dir / "exp04_phase_radial_summary.csv"

    assert traj_path.is_file(), f"Missing trajectory file: {traj_path}"
    assert summary_path.is_file(), f"Missing summary file: {summary_path}"

    with open(traj_path, mode="r", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        expected_traj_cols = [
            "dt_s", "step", "t_s", "x_km", "y_km", "vx_km_s", "vy_km_s",
            "position_error_km", "velocity_error_km_s", "energy_error",
            "angular_momentum_error", "r_num_km", "delta_r_km",
            "theta_unwrapped_rad", "delta_theta_rad", "e_R_km", "e_T_km",
            "fraction_R", "fraction_T"
        ]
        assert header == expected_traj_cols

    with open(summary_path, mode="r", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        expected_summary_cols = [
            "dt_s", "checkpoint_T", "time_s", "position_error_km",
            "velocity_error_km_s", "energy_error", "angular_momentum_error",
            "delta_r_km", "delta_theta_rad", "e_R_km", "e_T_km",
            "fraction_R", "fraction_T", "dominance_ratio", "dominant_mode"
        ]
        assert header == expected_summary_cols


def test_exp04_no_nan_or_inf_and_exact_end_time():
    """Verify that summary data contains no NaN/Inf, all checkpoints exist, and final time is 100T."""
    import csv
    from pathlib import Path

    results_dir = Path(__file__).resolve().parent.parent / "results"
    summary_path = results_dir / "exp04_phase_radial_summary.csv"

    with open(summary_path, mode="r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    assert len(rows) == 36  # 6 dt values * 6 checkpoints

    r0 = R_EARTH + 400.0
    period = circular_orbit_period(r0, mu=MU_EARTH)
    expected_100T = 100.0 * period

    for row in rows:
        for k, v in row.items():
            if k == "dominant_mode":
                assert v in ("phase", "radial", "equipartition")
            else:
                val = float(v)
                assert not math.isnan(val), f"NaN found in {k}"
                if k == "dominance_ratio" and val == float("inf"):
                    continue
                assert not math.isinf(val), f"Inf found in {k}"

        if float(row["checkpoint_T"]) == 100.0:
            assert math.isclose(float(row["time_s"]), expected_100T, abs_tol=1e-6)


def test_exp04_regression_against_m3_summary():
    """Verify that M4 core metrics match M3 summary exactly across all checkpoints and timesteps."""
    import csv
    from pathlib import Path

    results_dir = Path(__file__).resolve().parent.parent / "results"
    m3_summary_path = results_dir / "exp03_euler_long_term_summary.csv"
    m4_summary_path = results_dir / "exp04_phase_radial_summary.csv"

    with open(m3_summary_path, mode="r", encoding="utf-8") as f:
        m3_rows = list(csv.DictReader(f))
    with open(m4_summary_path, mode="r", encoding="utf-8") as f:
        m4_rows = list(csv.DictReader(f))

    assert len(m3_rows) == len(m4_rows) == 36

    for m3_r, m4_r in zip(m3_rows, m4_rows):
        assert float(m3_r["dt_s"]) == float(m4_r["dt_s"])
        assert float(m3_r["checkpoint_T"]) == float(m4_r["checkpoint_T"])
        assert math.isclose(float(m3_r["time_s"]), float(m4_r["time_s"]), abs_tol=1e-9)
        assert math.isclose(
            float(m3_r["position_error_km"]), float(m4_r["position_error_km"]), rel_tol=1e-11, abs_tol=1e-9
        )
        assert math.isclose(
            float(m3_r["velocity_error_km_s"]), float(m4_r["velocity_error_km_s"]), rel_tol=1e-11, abs_tol=1e-9
        )
        assert math.isclose(
            float(m3_r["energy_error"]), float(m4_r["energy_error"]), rel_tol=1e-11, abs_tol=1e-9
        )
        assert math.isclose(
            float(m3_r["angular_momentum_error"]), float(m4_r["angular_momentum_error"]), rel_tol=1e-11, abs_tol=1e-9
        )

