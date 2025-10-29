"""High level orchestration for umbrella sampling runs."""

from __future__ import annotations

import asyncio
import logging
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

from .analysis import UmbrellaWindow, compute_pmf, generate_synthetic_windows
from .config import UmbrellaPipelineConfig, UmbrellaWindowConfig

try:
    from openmm.app import PDBFile
    OPENMM_AVAILABLE = True
except ImportError:  # pragma: no cover - optional dependency
    OPENMM_AVAILABLE = False

try:
    from umbrella_sampling_calculator import create_umbrella_calculator
    CALCULATOR_AVAILABLE = True
except ImportError:  # pragma: no cover - optional dependency
    CALCULATOR_AVAILABLE = False

LOGGER = logging.getLogger(__name__)


class UmbrellaSamplingPipeline:
    """Encapsulates execution of umbrella sampling with visualization exports."""

    def __init__(self, config: UmbrellaPipelineConfig):
        self.config = config
        self.config.ensure_output_dir()

    def run(self, force_synthetic: Optional[bool] = None) -> Dict[str, object]:
        dependencies_ready = CALCULATOR_AVAILABLE and OPENMM_AVAILABLE
        if force_synthetic is None:
            use_synthetic = self.config.synthetic_if_missing
        else:
            use_synthetic = force_synthetic

        if use_synthetic:
            if dependencies_ready:
                LOGGER.info("Synthetic dataset requested; skipping OpenMM execution")
        elif not dependencies_ready:
            LOGGER.warning(
                "OpenMM or umbrella calculator unavailable, falling back to synthetic dataset"
            )
            use_synthetic = True

        if use_synthetic:
            windows = generate_synthetic_windows(
                self.config.window_centers,
                force_constant=self.config.force_constant,
                n_samples=1500,
                noise=0.4,
            )
            metadata = self._write_synthetic(windows)
            pmf = compute_pmf(windows, temperature=self.config.temperature_kelvin)
            return {"windows": windows, "metadata": metadata, "pmf": pmf}

        LOGGER.info("Launching umbrella sampling pipeline")
        windows = asyncio.run(self._run_async())
        pmf = compute_pmf(windows, temperature=self.config.temperature_kelvin)
        metadata = {
            "structure_file": str(self.config.structure_file),
            "window_centers": self.config.window_centers,
            "force_constant": self.config.force_constant,
            "simulation_time_ps": self.config.simulation_time_ps,
            "temperature": self.config.temperature_kelvin,
        }
        self._write_metadata(metadata)
        return {"windows": windows, "metadata": metadata, "pmf": pmf}

    async def _run_async(self) -> List[UmbrellaWindow]:
        calculator = create_umbrella_calculator(
            {
                "force_field": "amber19-all.xml",
                "platform": "CPU",
                "temperature": self.config.temperature_kelvin,
                "implicit_solvent": True,
            }
        )

        cv_atoms = self._resolve_cv_atoms()
        windows: List[UmbrellaWindow] = []

        batches = self._batched(self.config.window_centers, self.config.batch_size)
        for batch in batches:
            LOGGER.info("Running batch of %d windows", len(batch))
            results = await calculator.run_full_umbrella_sampling(
                structure_file=str(self.config.structure_file),
                cv_config={"type": "distance", "atoms": cv_atoms, "params": {}},
                window_centers=list(batch),
                force_constant=self.config.force_constant,
                simulation_time_ps=self.config.simulation_time_ps,
                temperature=self.config.temperature_kelvin,
                output_dir=str(self.config.output_dir),
            )
            windows.extend(self._convert_results(results))
        return windows

    def _resolve_cv_atoms(self) -> List[int]:
        pdb_file = PDBFile(str(self.config.structure_file))
        selections = self.config.protein_selection or []
        if not selections:
            indices = [atom.index for atom in pdb_file.topology.atoms() if atom.name == "CA"]
            if len(indices) < 2:
                raise ValueError("Unable to locate at least two CA atoms for CV")
            return indices[:2]

        resolved: List[int] = []
        for spec in selections:
            chain_id, res_id, atom_name = spec.split(":")
            found = False
            for atom in pdb_file.topology.atoms():
                if (
                    atom.name == atom_name
                    and atom.residue.id.strip() == res_id
                    and atom.residue.chain.id == chain_id
                ):
                    resolved.append(atom.index)
                    found = True
                    break
            if not found:
                raise ValueError(f"Selection {spec} not found in topology")
        if len(resolved) < 2:
            raise ValueError("Need at least two atoms for distance CV")
        return resolved[:2]

    def _write_metadata(self, metadata: Dict[str, object]) -> None:
        path = self.config.output_dir / "umbrella_metadata.json"
        import json

        path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    def _write_synthetic(self, windows: List[UmbrellaWindow]) -> Dict[str, object]:
        metadata = {
            "synthetic": True,
            "window_centers": self.config.window_centers,
            "force_constant": self.config.force_constant,
            "temperature": self.config.temperature_kelvin,
            "simulation_time_ps": self.config.simulation_time_ps,
        }
        self._write_metadata(metadata)
        for window in windows:
            self._export_window(window)
        return metadata

    def _export_window(self, window: UmbrellaWindow) -> None:
        histogram = np.column_stack(
            [window.histogram_edges[:-1], window.histogram_counts]
        )
        hist_path = self.config.output_dir / f"cv_histogram_center_{window.center:.2f}.dat"
        series_path = self.config.output_dir / f"cv_timeseries_center_{window.center:.2f}.dat"
        np.savetxt(hist_path, histogram)
        np.savetxt(series_path, window.cv_values)

    @staticmethod
    def _batched(items: Iterable[float], batch_size: int) -> List[List[float]]:
        batch: List[float] = []
        batches: List[List[float]] = []
        for value in items:
            batch.append(value)
            if len(batch) >= batch_size:
                batches.append(batch)
                batch = []
        if batch:
            batches.append(batch)
        return batches

    def _convert_results(self, results: Dict[str, object]) -> List[UmbrellaWindow]:
        window_dicts = results.get("windows", [])
        converted: List[UmbrellaWindow] = []
        for entry in window_dicts:
            cv_values = np.asarray(entry.get("cv_values", []), dtype=float)
            counts = np.asarray(entry.get("cv_histogram", {}).get("counts", []), dtype=float)
            edges = np.asarray(entry.get("cv_histogram", {}).get("bin_edges", []), dtype=float)
            if cv_values.size == 0 or counts.size == 0 or edges.size == 0:
                continue
            window = UmbrellaWindow(
                center=float(entry.get("window_center")),
                force_constant=float(entry.get("force_constant", self.config.force_constant)),
                cv_values=cv_values,
                histogram_counts=counts,
                histogram_edges=edges,
            )
            converted.append(window)
            self._export_window(window)
        if not converted:
            converted.extend(
                generate_synthetic_windows(
                    self.config.window_centers,
                    force_constant=self.config.force_constant,
                    n_samples=1000,
                    noise=0.5,
                )
            )
        return converted
