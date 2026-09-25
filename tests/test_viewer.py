"""Unit tests for the 2D orbital trajectory viewer."""

import os
from pathlib import Path
import tempfile

import matplotlib.pyplot as plt
import pytest

from src.models.state import State2D
from src.physics.gravity import R_EARTH
from src.visualization.viewer import (
    OrbitViewer,
    simulate_methods_for_viewing,
    plot_static_orbit,
    create_orbit_animation,
    export_html_viewer,
)


def test_simulate_methods_for_viewing_produces_valid_trajectories():
    results = simulate_methods_for_viewing(
        methods=["Euler", "RK4", "Velocity Verlet"],
        num_orbits=0.25,
        dt=60.0,
    )
    assert set(results.keys()) == {"Euler", "RK4", "Velocity Verlet"}
    for name, tdata in results.items():
        assert len(tdata.trajectory) > 5
        assert len(tdata.times) == len(tdata.trajectory)
        assert tdata.times[0] == 0.0
        # Check starting state matches 400km circular orbit
        first = tdata.trajectory[0]
        assert abs(first.r - (R_EARTH + 400.0)) < 1e-3


def test_plot_static_orbit_renders_earth_and_equal_aspect():
    results = simulate_methods_for_viewing(
        methods=["Velocity Verlet"],
        num_orbits=0.1,
        dt=60.0,
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        png_path = Path(tmpdir) / "test_orbit.png"
        fig, ax = plot_static_orbit(results, save_path=png_path)
        assert png_path.exists()
        assert png_path.stat().st_size > 0
        # Check aspect ratio (matplotlib represents 'equal' as 1.0)
        assert ax.get_aspect() in ("equal", 1.0)
        plt.close(fig)


def test_export_html_viewer_generates_valid_html_canvas():
    results = simulate_methods_for_viewing(
        methods=["Euler", "Velocity Verlet"],
        num_orbits=0.1,
        dt=60.0,
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        html_path = Path(tmpdir) / "test_viewer.html"
        out_file = export_html_viewer(results, output_path=html_path)
        assert out_file.exists()
        content = out_file.read_text(encoding="utf-8")
        assert "canvas id=\"orbitCanvas\"" in content
        assert "R_EARTH" in content
        assert "Velocity Verlet" in content
        assert "Euler" in content
        assert "Telemetry HUD" in content


def test_orbit_viewer_class_orchestration():
    viewer = OrbitViewer(methods=["RK4"], num_orbits=0.1, dt=60.0)
    data = viewer.simulate()
    assert "RK4" in data
    with tempfile.TemporaryDirectory() as tmpdir:
        html_path = Path(tmpdir) / "orbit.html"
        viewer.export_html(html_path)
        assert html_path.exists()
