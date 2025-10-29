"""Wrapper to run umbrella sampling pipeline and export visualization assets."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

from Chronosfold.umbrella_suite.config import UmbrellaPipelineConfig
from Chronosfold.umbrella_suite.pipeline import UmbrellaSamplingPipeline
from Chronosfold.umbrella_suite.visualization import plot_umbrella_diagnostics

logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Umbrella sampling visualization wrapper")
    parser.add_argument(
        "--structure",
        type=Path,
        default=Path("Chronosfold/WNK/5DRB.pdb"),
        help="Input PDB file",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("umbrella_results/wnk_videos"),
        help="Target directory for outputs",
    )
    parser.add_argument(
        "--animation",
        action="store_true",
        help="Generate MP4 animation of CV histograms",
    )
    parser.add_argument(
        "--frames",
        type=int,
        default=120,
        help="Frames for animation",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=24,
        help="Frames per second for animation",
    )
    return parser.parse_args()


def build_config(structure: Path, out_dir: Path) -> UmbrellaPipelineConfig:
    centers = [8.0, 9.0, 10.0, 11.0]
    return UmbrellaPipelineConfig(
        structure_file=structure,
        protein_selection=["A:245:CA", "A:292:CA"],
        window_centers=centers,
        force_constant=12.0,
        simulation_time_ps=25.0,
        output_dir=out_dir,
        temperature_kelvin=300.0,
        batch_size=2,
    )


def generate_animation(out_dir: Path, windows, fps: int, frames: int) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    centers = [w.center for w in windows]
    histograms = [w.histogram_counts for w in windows]
    edges = windows[0].histogram_edges
    bar_container = ax.bar(edges[:-1], histograms[0], width=np.diff(edges), align="edge", alpha=0.6)
    ax.set_ylim(0, max(hist.max() for hist in histograms) * 1.2)
    ax.set_xlabel("CV (Å)")
    ax.set_ylabel("Probability")

    def update(frame: int):
        idx = frame % len(histograms)
        for bar, height in zip(bar_container, histograms[idx]):
            bar.set_height(height)
        ax.set_title(f"Umbrella histogram window={centers[idx]:.1f} Å")
        return bar_container

    writer = animation.FFMpegWriter(fps=fps)
    mp4_path = out_dir / "umbrella_histograms.mp4"
    with writer.saving(fig, str(mp4_path), dpi=200):
        for frame in range(frames):
            update(frame)
            writer.grab_frame()
    plt.close(fig)
    LOGGER.info("Animation saved to %s", mp4_path)


def main() -> None:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    config = build_config(args.structure, args.out)
    pipeline = UmbrellaSamplingPipeline(config)
    result = pipeline.run()
    windows = result["windows"]
    pmf_df = result["pmf"]

    diagnostics_path = args.out / "umbrella_diagnostics.png"
    plot_umbrella_diagnostics(windows, pmf_df, title_suffix="WNK", save_path=diagnostics_path)
    LOGGER.info("Diagnostics saved to %s", diagnostics_path)

    if args.animation and windows:
        generate_animation(args.out, windows, fps=args.fps, frames=args.frames)


if __name__ == "__main__":
    main()
