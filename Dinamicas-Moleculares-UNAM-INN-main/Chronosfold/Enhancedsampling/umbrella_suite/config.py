"""Configuration objects for umbrella sampling pipelines."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional


def _default_window_centers() -> List[float]:
    return [float(x) for x in range(6, 21, 2)]


@dataclass
class UmbrellaWindowConfig:
    """Single umbrella sampling window configuration."""

    center: float
    force_constant: float = 12.0
    simulation_time_ps: float = 250.0


@dataclass
class UmbrellaPipelineConfig:
    """High level configuration for umbrella sampling orchestration."""

    structure_file: Path
    protein_selection: Optional[List[str]] = None
    ligand_selection: Optional[List[str]] = None
    window_centers: List[float] = field(default_factory=_default_window_centers)
    force_constant: float = 12.0
    simulation_time_ps: float = 250.0
    output_dir: Path = Path("umbrella_results/wnk_pilot")
    temperature_kelvin: float = 300.0
    batch_size: int = 5
    synthetic_if_missing: bool = True

    def to_windows(self) -> List[UmbrellaWindowConfig]:
        return [
            UmbrellaWindowConfig(
                center=center,
                force_constant=self.force_constant,
                simulation_time_ps=self.simulation_time_ps,
            )
            for center in self.window_centers
        ]

    def ensure_output_dir(self) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)
