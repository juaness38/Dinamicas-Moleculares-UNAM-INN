"""Unit tests for umbrella suite quick pilots."""

from pathlib import Path

import numpy as np

from Chronosfold.umbrella_suite.analysis import compute_pmf, generate_synthetic_windows
from Chronosfold.umbrella_suite.config import UmbrellaPipelineConfig
from Chronosfold.umbrella_suite.pipeline import UmbrellaSamplingPipeline


def test_generate_synthetic_windows_overlap():
    centers = [8.0, 9.0, 10.0]
    windows = generate_synthetic_windows(centers, force_constant=10.0, n_samples=500)
    assert len(windows) == len(centers)
    diffs = [abs(w.center - w.mean_cv) for w in windows]
    assert np.mean(diffs) < 1.0


def test_pipeline_returns_windows(tmp_path: Path):
    config = UmbrellaPipelineConfig(
        structure_file=Path("Chronosfold/WNK/5DRB.pdb"),
        window_centers=[8.0, 9.0],
        force_constant=12.0,
        simulation_time_ps=5.0,
        output_dir=tmp_path,
    )
    pipeline = UmbrellaSamplingPipeline(config)
    result = pipeline.run()
    windows = result["windows"]
    assert windows, "Pipeline should return at least one window"
    pmf = compute_pmf(windows, temperature=300.0)
    assert not pmf.empty
