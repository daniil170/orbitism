"""State representations for 2D orbital dynamics."""

from dataclasses import dataclass
import math
from typing import Sequence, Tuple


@dataclass(frozen=True)
class State2D:
    """Represents the 2D kinematic state of a body.

    Coordinate system: Geocentric Cartesian 2D.
    Units:
        x, y: kilometers (km)
        vx, vy: kilometers per second (km/s)
    """

    x: float
    y: float
    vx: float
    vy: float

    @property
    def r(self) -> float:
        """Distance from coordinate origin in km."""
        return math.hypot(self.x, self.y)

    @property
    def v(self) -> float:
        """Magnitude of velocity vector in km/s."""
        return math.hypot(self.vx, self.vy)

    def to_tuple(self) -> Tuple[float, float, float, float]:
        """Convert state to a tuple (x, y, vx, vy)."""
        return (self.x, self.y, self.vx, self.vy)

    def to_list(self) -> list[float]:
        """Convert state to a list [x, y, vx, vy]."""
        return [self.x, self.y, self.vx, self.vy]

    @classmethod
    def from_sequence(cls, seq: Sequence[float]) -> "State2D":
        """Create a State2D instance from a sequence of 4 floats [x, y, vx, vy]."""
        if len(seq) != 4:
            raise ValueError(f"Expected 4 elements [x, y, vx, vy], got {len(seq)}")
        return cls(x=float(seq[0]), y=float(seq[1]), vx=float(seq[2]), vy=float(seq[3]))


@dataclass(frozen=True)
class StateDerivative2D:
    """Represents the time derivative of a 2D state vector.

    Units:
        vx, vy: km/s (rate of position change)
        ax, ay: km/s^2 (rate of velocity change / gravitational acceleration)
    """

    vx: float
    vy: float
    ax: float
    ay: float

    def to_tuple(self) -> Tuple[float, float, float, float]:
        """Convert derivative to a tuple (vx, vy, ax, ay)."""
        return (self.vx, self.vy, self.ax, self.ay)

    def to_list(self) -> list[float]:
        """Convert derivative to a list [vx, vy, ax, ay]."""
        return [self.vx, self.vy, self.ax, self.ay]

    @classmethod
    def from_sequence(cls, seq: Sequence[float]) -> "StateDerivative2D":
        """Create a StateDerivative2D instance from a sequence of 4 floats [vx, vy, ax, ay]."""
        if len(seq) != 4:
            raise ValueError(f"Expected 4 elements [vx, vy, ax, ay], got {len(seq)}")
        return cls(vx=float(seq[0]), vy=float(seq[1]), ax=float(seq[2]), ay=float(seq[3]))

