"""Simulation runner for orbital trajectories."""

from typing import Callable, List, Optional, Union

from src.integrators.euler import EulerIntegrator
from src.models.state import State2D, StateDerivative2D
from src.physics.gravity import two_body_derivatives, MU_EARTH


class Simulator:
    """Simulates orbital motion using numerical integration.

    Steps through time using a specified integrator (defaults to Explicit Euler)
    and equations of motion (defaults to Earth two-body gravity).
    """

    def __init__(
        self,
        initial_state: Optional[State2D] = None,
        dt: Optional[float] = None,
        num_steps: Optional[int] = None,
        t_final: Optional[float] = None,
        mu: float = MU_EARTH,
        derivative_fn: Optional[Callable[[float, State2D], StateDerivative2D]] = None,
        integrator: Optional[Union[EulerIntegrator, Any]] = None,
    ):
        """Initialize simulator instance.

        Args:
            initial_state: Optional initial kinematic state.
            dt: Optional timestep in seconds.
            num_steps: Optional number of integration steps.
            t_final: Optional total simulation duration in seconds.
            mu: Gravitational parameter in km^3/s^2.
            derivative_fn: Optional custom derivative function dstate/dt = f(t, state).
            integrator: Optional integrator instance or integrator class (defaults to EulerIntegrator).
        """
        self.initial_state = initial_state
        self.dt = dt
        self.num_steps = num_steps
        self.t_final = t_final
        self.mu = mu

        self.derivative_fn = derivative_fn or (
            lambda t, state: two_body_derivatives(t, state, mu=self.mu)
        )
        if integrator is None:
            self.integrator = EulerIntegrator(self.derivative_fn)
        elif isinstance(integrator, type):
            self.integrator = integrator(self.derivative_fn)
        else:
            self.integrator = integrator

        self.trajectory: List[State2D] = []
        self.times: List[float] = []

    def run(
        self,
        initial_state: Optional[State2D] = None,
        dt: Optional[float] = None,
        num_steps: Optional[int] = None,
        t_final: Optional[float] = None,
    ) -> List[State2D]:
        """Execute simulation and return trajectory.

        Args:
            initial_state: Starting State2D (overrides constructor value if provided).
            dt: Timestep in seconds (overrides constructor value if provided).
            num_steps: Number of steps (overrides constructor value if provided).
            t_final: Final simulation time in seconds (overrides constructor value if provided).

        Returns:
            List of State2D instances [state_0, state_1, ..., state_N].

        Raises:
            ValueError: If initial_state, dt, or duration parameters are missing or invalid.
        """
        state = initial_state if initial_state is not None else self.initial_state
        step_dt = dt if dt is not None else self.dt
        steps = num_steps if num_steps is not None else self.num_steps
        final_time = t_final if t_final is not None else self.t_final

        if state is None:
            raise ValueError("Initial state must be specified.")
        if step_dt is None or step_dt <= 0:
            raise ValueError(f"Timestep dt must be positive, got {step_dt}.")
        if steps is None and final_time is None:
            raise ValueError("Either number of steps or t_final must be specified.")
        if steps is not None and steps < 0:
            raise ValueError(f"Number of steps must be non-negative, got {steps}.")
        if final_time is not None and final_time < 0:
            raise ValueError(f"Final time t_final must be non-negative, got {final_time}.")

        current_state = state
        current_time = 0.0

        self.trajectory = [current_state]
        self.times = [current_time]

        if final_time is not None:
            while current_time < final_time - 1e-12:
                dt_step = min(step_dt, final_time - current_time)
                current_state = self.integrator.step(
                    t=current_time,
                    state=current_state,
                    dt=dt_step,
                )
                current_time += dt_step
                if abs(current_time - final_time) < 1e-12:
                    current_time = final_time
                self.trajectory.append(current_state)
                self.times.append(current_time)
        else:
            for _ in range(steps):
                current_state = self.integrator.step(
                    t=current_time,
                    state=current_state,
                    dt=step_dt,
                )
                current_time += step_dt
                self.trajectory.append(current_state)
                self.times.append(current_time)

        return self.trajectory
