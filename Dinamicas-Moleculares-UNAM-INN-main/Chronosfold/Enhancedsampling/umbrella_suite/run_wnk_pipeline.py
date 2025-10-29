"""Command line runner for the WNK umbrella sampling pilot."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import List

from .analysis import compute_pmf
from .config import UmbrellaPipelineConfig
from .pipeline import UmbrellaSamplingPipeline
from .visualization import plot_umbrella_diagnostics

logging.basicConfig(level=logging.INFO)
LOGGER = logging.getLogger(__name__)


DEFAULT_SELECTION = [
    "A:245:CA",  # beta3 reference residue
    "A:292:CA",  # helix alphaC representative
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run WNK umbrella sampling pilot")
    parser.add_argument(
        "--structure",
        type=Path,
        default=Path("Chronosfold/WNK/5DRB.pdb"),
        help="Path to the WNK PDB structure",
    )
    parser.add_argument(
        "--windows",
        type=int,
        default=6,
        help="Number of umbrella windows between 8 and 14 Å",
    )
    parser.add_argument(
        "--time-ps",
        type=float,
        default=50.0,
        help="Simulation time per window in picoseconds",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("umbrella_results/wnk_pilot"),
        help="Output directory for results",
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=Path("umbrella_results/wnk_pilot/wnk_pilot_diagnostics.png"),
        help="Path to save the diagnostics plot",
    )
    parser.add_argument(
        "--synthetic",
        dest="synthetic",
        action="store_true",
        help="Force synthetic dataset generation (default)",
    )
    parser.add_argument(
        "--no-synthetic",
        dest="synthetic",
        action="store_false",
        help="Attempt real OpenMM execution if dependencies are available",
    )
    parser.set_defaults(synthetic=True)
    return parser.parse_args()


def build_window_centers(count: int) -> List[float]:
    count = max(2, min(count, 12))
    return [8.0 + i * (6.0 / max(count - 1, 1)) for i in range(count)]


def main() -> None:
    args = parse_args()
    centers = build_window_centers(args.windows)

    config = UmbrellaPipelineConfig(
        structure_file=args.structure,
        protein_selection=DEFAULT_SELECTION,
        window_centers=centers,
        force_constant=12.0,
        simulation_time_ps=args.time_ps,
        output_dir=args.out,
        temperature_kelvin=300.0,
        batch_size=3,
        synthetic_if_missing=args.synthetic,
    )

    pipeline = UmbrellaSamplingPipeline(config)
    result = pipeline.run(force_synthetic=args.synthetic)
    windows = result["windows"]
    pmf_df = result["pmf"]

    LOGGER.info("Generated %d windows", len(windows))
    LOGGER.info("PMF range: %.2f kcal/mol", pmf_df["pmf"].max() if not pmf_df.empty else 0.0)

    plot_umbrella_diagnostics(windows, pmf_df, title_suffix="WNK pilot", save_path=args.plot)
    LOGGER.info("Diagnostics plot saved to %s", args.plot)
    LOGGER.info("Sample metadata written to %s", args.out / "umbrella_metadata.json")


if __name__ == "__main__":
    main()
