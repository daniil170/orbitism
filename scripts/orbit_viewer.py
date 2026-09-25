#!/usr/bin/env python3
"""CLI utility to launch the OrbitSim 2D orbital trajectory viewer.

Usage examples:
    # Generate and open interactive HTML/Canvas 2D viewer in browser (default):
    python scripts/orbit_viewer.py

    # Generate static 2D comparison PNG plot:
    python scripts/orbit_viewer.py --mode static --output results/orbit_comparison.png

    # Simulate 5 orbits with dt = 60s to clearly see Euler's outward spiral vs Verlet's stability:
    python scripts/orbit_viewer.py --orbits 5 --dt 60 --open
"""

import argparse
from pathlib import Path
import subprocess
import sys

# Ensure project root is in sys.path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.visualization.viewer import OrbitViewer, simulate_methods_for_viewing, plot_static_orbit, export_html_viewer


def parse_args():
    parser = argparse.ArgumentParser(description="OrbitSim 2D Orbital Viewer")
    parser.add_argument(
        "--mode",
        choices=["html", "static", "animate"],
        default="html",
        help="Visualization mode: 'html' (interactive browser viewer), 'static' (PNG image), 'animate' (matplotlib animation)",
    )
    parser.add_argument(
        "--orbits",
        type=float,
        default=2.0,
        help="Number of orbital periods to simulate (default: 2.0)",
    )
    parser.add_argument(
        "--dt",
        type=float,
        default=30.0,
        help="Timestep in seconds (default: 30.0)",
    )
    parser.add_argument(
        "--methods",
        nargs="+",
        default=["Euler", "RK4", "Velocity Verlet"],
        help="Methods to include (default: Euler RK4 'Velocity Verlet')",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Custom output file path",
    )
    parser.add_argument(
        "--open",
        action="store_true",
        help="Automatically open output file with system viewer (macOS 'open')",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    print(f"=== OrbitSim 2D Orbit Viewer ===")
    print(f"Methods:    {args.methods}")
    print(f"Orbits:     {args.orbits}")
    print(f"Timestep:   {args.dt} s")
    print(f"Mode:       {args.mode}")

    viewer = OrbitViewer(
        methods=args.methods,
        num_orbits=args.orbits,
        dt=args.dt,
    )
    print("Simulating trajectories with core numerical integrators...")
    viewer.simulate()
    print("Simulation complete.")

    if args.mode == "html":
        out_path = Path(args.output) if args.output else Path("results/orbit_viewer.html")
        viewer.export_html(output_path=out_path)
        print(f"Interactive HTML viewer generated: {out_path.resolve()}")
        if args.open:
            subprocess.run(["open", str(out_path.resolve())], check=False)
            print("Opened in default browser.")

    elif args.mode == "static":
        out_path = Path(args.output) if args.output else Path("results/orbit_comparison.png")
        viewer.show_static(save_path=out_path)
        print(f"Static 2D trajectory plot saved: {out_path.resolve()}")
        if args.open:
            subprocess.run(["open", str(out_path.resolve())], check=False)

    elif args.mode == "animate":
        import matplotlib.pyplot as plt
        fig, anim = viewer.show_animation()
        if args.output:
            out_path = Path(args.output)
            if out_path.suffix == ".html":
                with open(out_path, "w", encoding="utf-8") as f:
                    f.write(anim.to_jshtml())
                print(f"Animation HTML saved: {out_path.resolve()}")
            else:
                anim.save(str(out_path))
                print(f"Animation saved: {out_path.resolve()}")
            if args.open:
                subprocess.run(["open", str(out_path.resolve())], check=False)
        else:
            print("Displaying animation window (close window to exit)...")
            plt.show()


if __name__ == "__main__":
    main()
