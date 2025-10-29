# ChronosFold Umbrella Suite

Pipeline and visualization scaffolding for umbrella sampling and WHAM analysis with a focus on WNK kinases.

## Components

- `config.py`: dataclasses for pipeline configuration and window definitions.
- `analysis.py`: data loading, synthetic dataset generation, and PMF computation (MBAR/WHAM fallback).
- `pipeline.py`: orchestration layer that runs the umbrella calculator when OpenMM is available or generates synthetic pilots otherwise.
- `visualization.py`: reusable Matplotlib diagnostics analogous to the OpenMM cookbook tutorial.
- `run_wnk_pipeline.py`: CLI entry point for picosecond-scale pilots on `WNK/5DRB.pdb`.

## Usage

Generate a quick pilot and diagnostics:

```bash
python -m Chronosfold.umbrella_suite.run_wnk_pipeline --time-ps 50 --windows 6
```

Run the video wrapper:

```bash
python Chronosfold/VIDEOSUITE/umbrella/run_umbrella_visualization.py --animation
```

Outputs are placed in `umbrella_results/wnk_pilot` by default.

## Testing

Unit tests covering synthetic generation and pipeline fallbacks live in `Chronosfold/tests/test_umbrella_suite.py`.
