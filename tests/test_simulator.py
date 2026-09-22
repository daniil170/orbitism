"""Unit tests for the Simulator runner."""

import pytest

from src.models.state import State2D
from src.physics.gravity import create_circular_orbit_state
from src.simulation.simulator import Simulator


def test_simulator_trajectory_length_and_initial_condition():
    """Simulator running N steps must return a trajectory with N + 1 states."""
    initial_state = create_circular_orbit_state(altitude_km=400.0)
    dt = 1.0  # seconds
    num_steps = 100

    sim = Simulator()
    trajectory = sim.run(initial_state=initial_state, dt=dt, num_steps=num_steps)

    assert len(trajectory) == num_steps + 1
    assert len(sim.times) == num_steps + 1

    # First state must match initial conditions exactly
    assert trajectory[0] == initial_state
    assert sim.times[0] == 0.0
    assert sim.times[-1] == pytest.approx(100.0)

    # Trajectory attribute on simulator must match returned trajectory
    assert sim.trajectory == trajectory


def test_simulator_constructor_arguments():
    """Simulator initialized via constructor arguments should run without additional args."""
    initial_state = create_circular_orbit_state(altitude_km=400.0)
    sim = Simulator(initial_state=initial_state, dt=0.5, num_steps=10)
    trajectory = sim.run()

    assert len(trajectory) == 11
    assert trajectory[0] == initial_state


def test_simulator_invalid_inputs_raise_errors():
    """Simulator must validate inputs and reject negative steps or non-positive dt."""
    sim = Simulator()
    valid_state = create_circular_orbit_state(altitude_km=400.0)

    with pytest.raises(ValueError, match="Initial state"):
        sim.run(initial_state=None, dt=1.0, num_steps=10)

    with pytest.raises(ValueError, match="Timestep dt must be positive"):
        sim.run(initial_state=valid_state, dt=-1.0, num_steps=10)

    with pytest.raises(ValueError, match="Number of steps"):
        sim.run(initial_state=valid_state, dt=1.0, num_steps=-5)

