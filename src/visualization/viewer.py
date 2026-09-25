"""2D orbital trajectory visualizer and animator for OrbitSim.

Displays 2D orbital motion around Earth:
- Earth centered at (0, 0) with physical radius R_E = 6378.137 km.
- Trajectory path curves.
- Instantaneous satellite position marker.
- Direction of motion velocity vector arrow.
- Coordinate axes with 1:1 equal aspect ratio.
- Real-time simulation metrics HUD (time, x, y, radius, altitude).
- Multi-method visual comparison (Euler, RK4, Velocity Verlet).
"""

from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.patches as patches

from src.models.state import State2D
from src.physics.gravity import R_EARTH, MU_EARTH, create_circular_orbit_state
from src.physics.analytical import circular_orbit_period
from src.integrators.euler import EulerIntegrator
from src.integrators.rk4 import RK4Integrator
from src.integrators.velocity_verlet import VelocityVerletIntegrator
from src.integrators.yoshida4 import Yoshida4Integrator
from src.simulation.simulator import Simulator

METHOD_COLORS: Dict[str, str] = {
    "Euler": "#e74c3c",           # Red
    "RK4": "#2ecc71",             # Green
    "Velocity Verlet": "#3498db", # Blue
    "Yoshida 4": "#9b59b6",       # Purple
}

METHOD_CLASSES = {
    "Euler": EulerIntegrator,
    "RK4": RK4Integrator,
    "Velocity Verlet": VelocityVerletIntegrator,
    "Yoshida 4": Yoshida4Integrator,
}


@dataclass
class TrajectoryData:
    """Container for simulated trajectory and timestamps."""

    method_name: str
    times: List[float]
    trajectory: List[State2D]
    color: str = "#3498db"


def simulate_methods_for_viewing(
    methods: Optional[Sequence[str]] = None,
    num_orbits: float = 2.0,
    dt: float = 30.0,
    altitude_km: float = 400.0,
    mu: float = MU_EARTH,
    r_earth: float = R_EARTH,
) -> Dict[str, TrajectoryData]:
    """Run standard OrbitSim simulations for visual comparison without mutating core.

    Args:
        methods: List of method names to simulate ('Euler', 'RK4', 'Velocity Verlet').
        num_orbits: Number of orbital periods to run.
        dt: Timestep in seconds.
        altitude_km: Initial circular orbit altitude in km.
        mu: Gravitational parameter.
        r_earth: Central body radius.

    Returns:
        Dictionary mapping method name to TrajectoryData.
    """
    if methods is None:
        methods = ["Euler", "RK4", "Velocity Verlet"]

    r0 = r_earth + altitude_km
    single_period = circular_orbit_period(r0, mu=mu)
    t_final = num_orbits * single_period
    init_state = create_circular_orbit_state(
        altitude_km=altitude_km, mu=mu, r_earth=r_earth
    )

    results: Dict[str, TrajectoryData] = {}
    for m in methods:
        if m not in METHOD_CLASSES:
            raise ValueError(f"Unknown integrator method '{m}'. Available: {list(METHOD_CLASSES.keys())}")
        cls = METHOD_CLASSES[m]
        sim = Simulator(initial_state=init_state, dt=dt, t_final=t_final, mu=mu, integrator=cls)
        traj = sim.run()
        results[m] = TrajectoryData(
            method_name=m,
            times=list(sim.times),
            trajectory=traj,
            color=METHOD_COLORS.get(m, "#9b59b6"),
        )

    return results


def plot_static_orbit(
    data: Union[TrajectoryData, Dict[str, TrajectoryData]],
    show_arrows: bool = True,
    r_earth: float = R_EARTH,
    title: Optional[str] = None,
    save_path: Optional[Union[str, Path]] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """Create static 2D planar plot of orbital trajectories.

    Args:
        data: Single TrajectoryData or dict of TrajectoryData.
        show_arrows: Whether to draw velocity direction vector arrows.
        r_earth: Earth radius in km.
        title: Optional plot title.
        save_path: Optional file path to save image.

    Returns:
        (fig, ax) Matplotlib figure and axes.
    """
    if isinstance(data, TrajectoryData):
        trajectories = {data.method_name: data}
    else:
        trajectories = data

    fig, ax = plt.subplots(figsize=(8, 8), dpi=100)

    # 1. Earth at origin
    earth_circle = patches.Circle(
        (0.0, 0.0),
        r_earth,
        facecolor="#2980b9",
        edgecolor="#1a5276",
        alpha=0.85,
        zorder=2,
        label=f"Earth (R_E = {r_earth:.0f} km)",
    )
    ax.add_patch(earth_circle)
    ax.plot(0.0, 0.0, "w+", markersize=10, markeredgewidth=1.5, zorder=3)

    # Track bounding box for equal aspect view
    max_coord = r_earth * 1.5

    # 2. Plot each trajectory
    for name, tdata in trajectories.items():
        xs = [s.x for s in tdata.trajectory]
        ys = [s.y for s in tdata.trajectory]
        max_coord = max(max_coord, max(abs(x) for x in xs), max(abs(y) for y in ys))

        # Trail path
        ax.plot(xs, ys, color=tdata.color, linewidth=1.5, label=f"{name} trajectory", zorder=4)

        # Initial point
        ax.scatter([xs[0]], [ys[0]], color=tdata.color, s=40, marker="o", edgecolors="black", zorder=5)

        # Final / current satellite point
        final_state = tdata.trajectory[-1]
        ax.scatter(
            [final_state.x],
            [final_state.y],
            color=tdata.color,
            s=80,
            marker="*",
            edgecolors="black",
            zorder=6,
            label=f"{name} current (r={final_state.r:.0f} km)",
        )

        # Velocity arrow
        if show_arrows and final_state.v > 0:
            scale_len = r_earth * 0.25
            vx_unit = final_state.vx / final_state.v
            vy_unit = final_state.vy / final_state.v
            ax.annotate(
                "",
                xy=(final_state.x + vx_unit * scale_len, final_state.y + vy_unit * scale_len),
                xytext=(final_state.x, final_state.y),
                arrowprops=dict(arrowstyle="->", color=tdata.color, lw=2.0),
                zorder=7,
            )

    # 3. Coordinate axes and styling
    ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5, zorder=1)
    ax.axvline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.5, zorder=1)
    ax.set_aspect("equal", adjustable="box")

    pad = max_coord * 1.15
    ax.set_xlim(-pad, pad)
    ax.set_ylim(-pad, pad)

    ax.set_xlabel("X (km)", fontsize=11)
    ax.set_ylabel("Y (km)", fontsize=11)
    ax.set_title(title or "OrbitSim — 2D Orbital Trajectory", fontsize=12, pad=12)
    ax.grid(True, linestyle=":", alpha=0.6, zorder=1)
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)

    fig.tight_layout()

    if save_path:
        fig.savefig(save_path)

    return fig, ax


def create_orbit_animation(
    data: Union[TrajectoryData, Dict[str, TrajectoryData]],
    r_earth: float = R_EARTH,
    frames: int = 150,
    interval_ms: int = 40,
    title: Optional[str] = None,
) -> Tuple[plt.Figure, FuncAnimation]:
    """Create animated 2D orbit playback using FuncAnimation.

    Args:
        data: TrajectoryData instance or dict of TrajectoryData.
        r_earth: Central body radius in km.
        frames: Number of animation frames to render.
        interval_ms: Milliseconds per frame.
        title: Animation title.

    Returns:
        (fig, anim) Tuple of figure and animation object.
    """
    if isinstance(data, TrajectoryData):
        trajectories = {data.method_name: data}
    else:
        trajectories = data

    fig, ax = plt.subplots(figsize=(8, 8), dpi=100)

    # Earth circle
    earth_circle = patches.Circle(
        (0.0, 0.0),
        r_earth,
        facecolor="#2980b9",
        edgecolor="#1a5276",
        alpha=0.85,
        zorder=2,
        label=f"Earth (R_E = {r_earth:.0f} km)",
    )
    ax.add_patch(earth_circle)
    ax.plot(0.0, 0.0, "w+", markersize=10, markeredgewidth=1.5, zorder=3)

    # Determine maximum extent
    max_coord = r_earth * 1.5
    total_steps = 0
    for tdata in trajectories.values():
        total_steps = max(total_steps, len(tdata.trajectory))
        for s in tdata.trajectory:
            max_coord = max(max_coord, abs(s.x), abs(s.y))

    pad = max_coord * 1.15
    ax.set_xlim(-pad, pad)
    ax.set_ylim(-pad, pad)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("X (km)", fontsize=11)
    ax.set_ylabel("Y (km)", fontsize=11)
    ax.set_title(title or "OrbitSim — 2D Trajectory Animation", fontsize=12, pad=12)
    ax.grid(True, linestyle=":", alpha=0.6, zorder=1)

    # Graphical elements per method
    lines = {}
    dots = {}
    arrows = {}
    for name, tdata in trajectories.items():
        (line,) = ax.plot([], [], color=tdata.color, linewidth=1.5, label=f"{name}", zorder=4)
        (dot,) = ax.plot([], [], color=tdata.color, marker="o", markersize=8, markeredgecolor="black", zorder=6)
        lines[name] = line
        dots[name] = dot

    hud_text = ax.text(
        0.03,
        0.97,
        "",
        transform=ax.transAxes,
        verticalalignment="top",
        fontsize=9,
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.85, edgecolor="#ccc"),
        zorder=10,
    )

    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)

    step_indices = [int(i * (total_steps - 1) / max(1, frames - 1)) for i in range(frames)]

    arrow_patch = [None]  # container to store temporary arrow annotation

    def init():
        for line in lines.values():
            line.set_data([], [])
        for dot in dots.values():
            dot.set_data([], [])
        hud_text.set_text("")
        return list(lines.values()) + list(dots.values()) + [hud_text]

    def update(frame_idx):
        idx = step_indices[frame_idx]
        hud_lines = []

        if arrow_patch[0] is not None:
            arrow_patch[0].remove()
            arrow_patch[0] = None

        first_data = next(iter(trajectories.values()))
        curr_t = first_data.times[min(idx, len(first_data.times) - 1)]
        hud_lines.append(f"Simulation Time: {curr_t:8.1f} s")

        scale_len = r_earth * 0.25

        for name, tdata in trajectories.items():
            curr_idx = min(idx, len(tdata.trajectory) - 1)
            sub_traj = tdata.trajectory[: curr_idx + 1]
            xs = [s.x for s in sub_traj]
            ys = [s.y for s in sub_traj]
            lines[name].set_data(xs, ys)

            curr_state = tdata.trajectory[curr_idx]
            dots[name].set_data([curr_state.x], [curr_state.y])

            r = curr_state.r
            alt = r - r_earth
            hud_lines.append(
                f"{name:15s}: x={curr_state.x:8.1f} km, y={curr_state.y:8.1f} km, r={r:8.1f} km (alt: {alt:6.1f} km)"
            )

        # Draw velocity arrow for primary (or first) method
        primary_name = "Velocity Verlet" if "Velocity Verlet" in trajectories else next(iter(trajectories.keys()))
        primary_idx = min(idx, len(trajectories[primary_name].trajectory) - 1)
        p_state = trajectories[primary_name].trajectory[primary_idx]
        if p_state.v > 0:
            vx_u = p_state.vx / p_state.v
            vy_u = p_state.vy / p_state.v
            arrow_patch[0] = ax.annotate(
                "",
                xy=(p_state.x + vx_u * scale_len, p_state.y + vy_u * scale_len),
                xytext=(p_state.x, p_state.y),
                arrowprops=dict(arrowstyle="->", color=trajectories[primary_name].color, lw=2.0),
                zorder=7,
            )

        hud_text.set_text("\n".join(hud_lines))
        return list(lines.values()) + list(dots.values()) + [hud_text]

    anim = FuncAnimation(fig, update, frames=frames, init_func=init, interval=interval_ms, blit=False)
    fig.tight_layout()
    return fig, anim


def export_html_viewer(
    data: Union[TrajectoryData, Dict[str, TrajectoryData]],
    output_path: Union[str, Path] = "results/orbit_viewer.html",
    r_earth: float = R_EARTH,
) -> Path:
    """Generate self-contained interactive 2D HTML/Canvas viewer.

    Requires zero browser dependencies: pure HTML5 Canvas + JS.
    Includes play/pause, speed controls, scrubber, equal scale, and live HUD.

    Args:
        data: TrajectoryData or dict of TrajectoryData.
        output_path: Output HTML file path.
        r_earth: Central body radius in km.

    Returns:
        Path to generated HTML file.
    """
    if isinstance(data, TrajectoryData):
        trajectories = {data.method_name: data}
    else:
        trajectories = data

    export_dict = {}
    for name, tdata in trajectories.items():
        # Subsample to at most 1200 points for smooth browser rendering
        step = max(1, len(tdata.trajectory) // 1000)
        sub_times = tdata.times[::step]
        sub_states = tdata.trajectory[::step]
        export_dict[name] = {
            "color": tdata.color,
            "times": [round(t, 2) for t in sub_times],
            "x": [round(s.x, 2) for s in sub_states],
            "y": [round(s.y, 2) for s in sub_states],
            "vx": [round(s.vx, 4) for s in sub_states],
            "vy": [round(s.vy, 4) for s in sub_states],
            "r": [round(s.r, 2) for s in sub_states],
        }

    json_payload = json.dumps(export_dict)

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>OrbitSim — 2D Interactive Orbit Viewer</title>
  <style>
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, monospace;
      margin: 0; padding: 20px;
      background: #11141a; color: #ecf0f1;
      display: flex; flex-direction: column; align-items: center;
    }}
    h1 {{ margin: 0 0 8px 0; font-size: 20px; font-weight: 600; color: #5dade2; }}
    .subtitle {{ font-size: 13px; color: #95a5a6; margin-bottom: 16px; }}
    .main-container {{
      display: flex; gap: 20px; background: #1a1e27; border-radius: 8px;
      padding: 16px; box-shadow: 0 4px 16px rgba(0,0,0,0.5);
    }}
    canvas {{
      background: #0b0d13; border: 1px solid #2c3e50; border-radius: 4px;
      cursor: crosshair;
    }}
    .panel {{
      width: 320px; display: flex; flex-direction: column; gap: 14px;
      font-size: 13px;
    }}
    .hud {{
      background: #0e1219; border: 1px solid #2a3342; border-radius: 4px;
      padding: 12px; font-family: monospace; font-size: 12px; line-height: 1.5;
    }}
    .hud-title {{ font-weight: bold; color: #3498db; margin-bottom: 6px; font-family: sans-serif; }}
    .method-row {{
      display: flex; align-items: center; gap: 8px; margin: 4px 0;
    }}
    .legend-box {{
      width: 12px; height: 12px; border-radius: 2px; display: inline-block;
    }}
    .controls {{
      display: flex; flex-direction: column; gap: 8px;
    }}
    .btn-row {{ display: flex; gap: 8px; }}
    button {{
      flex: 1; padding: 8px; border: none; border-radius: 4px;
      background: #34495e; color: #ecf0f1; font-weight: bold; cursor: pointer;
      transition: background 0.2s;
    }}
    button:hover {{ background: #4a627a; }}
    button.active {{ background: #2980b9; }}
    label {{ display: flex; justify-content: space-between; font-size: 12px; color: #bdc3c7; }}
    input[type=range] {{ width: 100%; }}
  </style>
</head>
<body>
  <h1>OrbitSim — 2D Orbital Dynamics Viewer</h1>
  <div class="subtitle">Interactive 2D Planar Verification (Earth centered at origin, 1:1 equal scale)</div>

  <div class="main-container">
    <canvas id="orbitCanvas" width="600" height="600"></canvas>

    <div class="panel">
      <div class="controls">
        <div class="btn-row">
          <button id="playBtn" class="active">Play</button>
          <button id="pauseBtn">Pause</button>
          <button id="resetBtn">Reset</button>
        </div>
        <label>Animation Progress <span id="frameLabel">0%</span></label>
        <input type="range" id="scrubber" min="0" max="1000" value="0">

        <label>Playback Speed <span id="speedLabel">1.0x</span></label>
        <input type="range" id="speedSlider" min="0.2" max="5.0" step="0.2" value="1.0">
      </div>

      <div class="hud">
        <div class="hud-title">Methods & Active Trails</div>
        <div id="methodToggles"></div>
      </div>

      <div class="hud" id="metricsHud">
        <div class="hud-title">Telemetry HUD</div>
        <div>Time: <span id="timeVal">0.0</span> s</div>
        <div id="methodTelemetry"></div>
      </div>

      <div style="font-size: 11px; color: #7f8c8d; line-height: 1.4;">
        * Scale: 1:1 isometric projection.<br>
        * Blue circle: Earth (R_E = {r_earth:.1f} km).<br>
        * Arrows indicate instantaneous velocity direction.
      </div>
    </div>
  </div>

  <script>
    const data = {json_payload};
    const R_EARTH = {r_earth};
    const canvas = document.getElementById("orbitCanvas");
    const ctx = canvas.getContext("2d");

    let isPlaying = true;
    let playbackSpeed = 1.0;
    let progress = 0.0; // 0 to 1
    const methodVisibility = {{}};

    // Populate method toggles
    const togglesDiv = document.getElementById("methodToggles");
    for (const m in data) {{
      methodVisibility[m] = true;
      const row = document.createElement("div");
      row.className = "method-row";
      row.innerHTML = `<span class="legend-box" style="background:${{data[m].color}}"></span>
        <input type="checkbox" id="chk_${{m}}" checked>
        <label for="chk_${{m}}" style="cursor:pointer; color:#ecf0f1;">${{m}}</label>`;
      togglesDiv.appendChild(row);
      document.getElementById(`chk_${{m}}`).addEventListener("change", (e) => {{
        methodVisibility[m] = e.target.checked;
      }});
    }}

    // Compute max bounding extent
    let maxR = R_EARTH * 1.5;
    for (const m in data) {{
      for (const r of data[m].r) {{
        if (r > maxR) maxR = r;
      }}
    }}
    const worldRadius = maxR * 1.15;

    function toCanvas(x, y) {{
      const cx = canvas.width / 2;
      const cy = canvas.height / 2;
      const scale = (canvas.width / 2) / worldRadius;
      return [cx + x * scale, cy - y * scale];
    }}

    function drawGrid() {{
      ctx.strokeStyle = "#1e2633";
      ctx.lineWidth = 1;
      const step = 2000;
      for (let w = -100000; w <= 100000; w += step) {{
        const [x1, y1] = toCanvas(w, -worldRadius);
        const [x2, y2] = toCanvas(w, worldRadius);
        ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x2, y2); ctx.stroke();
        const [hx1, hy1] = toCanvas(-worldRadius, w);
        const [hx2, hy2] = toCanvas(worldRadius, w);
        ctx.beginPath(); ctx.moveTo(hx1, hy1); ctx.lineTo(hx2, hy2); ctx.stroke();
      }}
      // Axes
      ctx.strokeStyle = "#475569";
      ctx.lineWidth = 1.5;
      const [ax0, ay0] = toCanvas(-worldRadius, 0);
      const [ax1, ay1] = toCanvas(worldRadius, 0);
      ctx.beginPath(); ctx.moveTo(ax0, ay0); ctx.lineTo(ax1, ay1); ctx.stroke();
      const [bx0, by0] = toCanvas(0, -worldRadius);
      const [bx1, by1] = toCanvas(0, worldRadius);
      ctx.beginPath(); ctx.moveTo(bx0, by0); ctx.lineTo(bx1, by1); ctx.stroke();
    }}

    function drawEarth() {{
      const [cx, cy] = toCanvas(0, 0);
      const scale = (canvas.width / 2) / worldRadius;
      ctx.beginPath();
      ctx.arc(cx, cy, R_EARTH * scale, 0, 2 * Math.PI);
      ctx.fillStyle = "#1b4f72";
      ctx.fill();
      ctx.strokeStyle = "#5dade2";
      ctx.lineWidth = 1.5;
      ctx.stroke();

      // Earth center marker
      ctx.strokeStyle = "#ffffff";
      ctx.beginPath();
      ctx.moveTo(cx - 5, cy); ctx.lineTo(cx + 5, cy);
      ctx.moveTo(cx, cy - 5); ctx.lineTo(cx, cy + 5);
      ctx.stroke();
    }}

    function render() {{
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      drawGrid();
      drawEarth();

      let telemetryHtml = "";
      let currentTime = 0.0;

      for (const m in data) {{
        if (!methodVisibility[m]) continue;
        const d = data[m];
        const n = d.x.length;
        const currIdx = Math.min(n - 1, Math.floor(progress * (n - 1)));
        currentTime = d.times[currIdx];

        // Draw trajectory trail
        ctx.beginPath();
        for (let i = 0; i <= currIdx; i++) {{
          const [cx, cy] = toCanvas(d.x[i], d.y[i]);
          if (i === 0) ctx.moveTo(cx, cy);
          else ctx.lineTo(cx, cy);
        }}
        ctx.strokeStyle = d.color;
        ctx.lineWidth = 2.0;
        ctx.stroke();

        // Draw initial point marker
        const [ix, iy] = toCanvas(d.x[0], d.y[0]);
        ctx.fillStyle = d.color;
        ctx.beginPath(); ctx.arc(ix, iy, 3, 0, 2*Math.PI); ctx.fill();

        // Draw current satellite marker
        const currX = d.x[currIdx];
        const currY = d.y[currIdx];
        const [sx, sy] = toCanvas(currX, currY);
        ctx.fillStyle = d.color;
        ctx.beginPath(); ctx.arc(sx, sy, 6, 0, 2*Math.PI); ctx.fill();
        ctx.strokeStyle = "#ffffff"; ctx.lineWidth = 1.5; ctx.stroke();

        // Velocity vector arrow
        const vx = d.vx[currIdx];
        const vy = d.vy[currIdx];
        const vMag = Math.hypot(vx, vy);
        if (vMag > 0) {{
          const arrowLen = R_EARTH * 0.25;
          const [tipX, tipY] = toCanvas(currX + (vx / vMag) * arrowLen, currY + (vy / vMag) * arrowLen);
          ctx.beginPath(); ctx.moveTo(sx, sy); ctx.lineTo(tipX, tipY);
          ctx.strokeStyle = d.color; ctx.lineWidth = 2.5; ctx.stroke();
        }}

        telemetryHtml += `<div style="color:${{d.color}}; margin-top:4px;"><strong>${{m}}</strong>: ` +
          `r=${{d.r[currIdx].toFixed(1)}} km (alt: ${{(d.r[currIdx] - R_EARTH).toFixed(1)}} km)<br>` +
          `&nbsp;&nbsp;x=${{currX.toFixed(1)}}, y=${{currY.toFixed(1)}} km</div>`;
      }}

      document.getElementById("timeVal").innerText = currentTime.toFixed(1);
      document.getElementById("methodTelemetry").innerHTML = telemetryHtml;
      document.getElementById("scrubber").value = Math.floor(progress * 1000);
      document.getElementById("frameLabel").innerText = Math.floor(progress * 100) + "%";
    }}

    function animate() {{
      if (isPlaying) {{
        progress += 0.001 * playbackSpeed;
        if (progress > 1.0) progress = 0.0;
      }}
      render();
      requestAnimationFrame(animate);
    }}

    // UI listeners
    document.getElementById("playBtn").addEventListener("click", () => {{
      isPlaying = true;
      document.getElementById("playBtn").classList.add("active");
      document.getElementById("pauseBtn").classList.remove("active");
    }});
    document.getElementById("pauseBtn").addEventListener("click", () => {{
      isPlaying = false;
      document.getElementById("pauseBtn").classList.add("active");
      document.getElementById("playBtn").classList.remove("active");
    }});
    document.getElementById("resetBtn").addEventListener("click", () => {{
      progress = 0.0;
      render();
    }});
    document.getElementById("scrubber").addEventListener("input", (e) => {{
      progress = e.target.value / 1000.0;
      render();
    }});
    document.getElementById("speedSlider").addEventListener("input", (e) => {{
      playbackSpeed = parseFloat(e.target.value);
      document.getElementById("speedLabel").innerText = playbackSpeed.toFixed(1) + "x";
    }});

    render();
    requestAnimationFrame(animate);
  </script>
</body>
</html>
"""
    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text(html_content, encoding="utf-8")
    return out_file


class OrbitViewer:
    """Convenient orchestrator for 2D OrbitSim visual simulations."""

    def __init__(
        self,
        methods: Optional[Sequence[str]] = None,
        num_orbits: float = 2.0,
        dt: float = 30.0,
        altitude_km: float = 400.0,
    ):
        self.methods = methods or ["Euler", "RK4", "Velocity Verlet"]
        self.num_orbits = num_orbits
        self.dt = dt
        self.altitude_km = altitude_km
        self.data: Dict[str, TrajectoryData] = {}

    def simulate(self) -> Dict[str, TrajectoryData]:
        """Run trajectory calculations."""
        self.data = simulate_methods_for_viewing(
            methods=self.methods,
            num_orbits=self.num_orbits,
            dt=self.dt,
            altitude_km=self.altitude_km,
        )
        return self.data

    def show_static(self, save_path: Optional[Union[str, Path]] = None) -> Tuple[plt.Figure, plt.Axes]:
        """Render static 2D comparison plot."""
        if not self.data:
            self.simulate()
        return plot_static_orbit(self.data, save_path=save_path)

    def show_animation(
        self, frames: int = 150, interval_ms: int = 40
    ) -> Tuple[plt.Figure, FuncAnimation]:
        """Render animated orbit."""
        if not self.data:
            self.simulate()
        return create_orbit_animation(self.data, frames=frames, interval_ms=interval_ms)

    def export_html(self, output_path: Union[str, Path] = "results/orbit_viewer.html") -> Path:
        """Export interactive HTML viewer."""
        if not self.data:
            self.simulate()
        return export_html_viewer(self.data, output_path=output_path)
