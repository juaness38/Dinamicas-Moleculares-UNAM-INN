# ChronosFold Video Suite (VIDEOSUITE)

This folder groups the four video pipelines referenced in the repo with simple wrappers and a minimal client:

- 2D synthetic trajectory animations (GIF/MP4)
- MD + MiniChronoGPT + optional Milvus states + timeline + Ramachandran
- PyMOL-based rendering from PDB/DCD and semantic annotations
- Alex Chen breakthrough visualization suite (ChronosFold → PCA → GIF with novelty detection)
- Umbrella sampling + WHAM diagnostics (WNK pilot)

Each wrapper delegates to existing code in the repository to avoid duplication.

## Prerequisites

- Common: Python 3.10+, pip, and ffmpeg in PATH (for MP4)
- 2D animations: matplotlib, seaborn, imageio, pillow
- MD demo: openmm, mdtraj, numpy, matplotlib, scikit-learn (optional: pymilvus for Milvus fetch)
- PyMOL: pymol-open-source and ffmpeg
- Alex Chen: mdtraj, scikit-learn, matplotlib, imageio (same as MD demo)

See `requirements.txt` for a combined (superset) list. Install only what you need.

## 1) 2D synthetic animation

- Wrapper: `2d/run_2d_animation.py`
- Output: MP4 (default `alanine_dipeptide_demo.mp4`)
- Internals: calls `chronosfold_scaffold/src/chronosfold/visualization/plotting.py` via the FastAPI router function `animate_mp4()`.

Example usage:
- Write MP4 in current folder with defaults (frames=200, fps=20):
  - `python VIDEOSUITE/2d/run_2d_animation.py`
- Custom output/fps/color:
  - `python VIDEOSUITE/2d/run_2d_animation.py --out demo.mp4 --frames 300 --fps 24 --color-by cluster`

## 2) MD + MiniChronoGPT + Milvus (optional) + animation

- Wrapper: `md_demo/run_minicronos_demo.py`
- Output: `reports/minicronos_demo.gif` (or MP4 if ffmpeg available) + `reports/minicronos_demo_meta.txt`
- Internals: delegates to `scripts/minicronos_demo_video.py` using `PYTHONPATH=chronosfold_scaffold/src`.

Examples:
- Only MD + generative continuation:
  - `python VIDEOSUITE/md_demo/run_minicronos_demo.py --steps 400 --pred-steps 60 --gif reports/minicronos_demo.gif`
- With Milvus (requires ZILLIZ_URI/ZILLIZ_TOKEN in `.env`):
  - `python VIDEOSUITE/md_demo/run_minicronos_demo.py --steps 400 --milvus --collection protein_chains --limit 20 --pred-steps 60 --gif reports/minicronos_demo.gif`

## 3) PyMOL + annotations → MP4

- Wrapper: `pymol/run_pymol_video.py`
- Output: `<out>/semantic_trajectory.mp4`
- Internals: runs `python -m ste.render.pymol_video` with provided PDB/DCD/annotations.

Example:
- `python VIDEOSUITE/pymol/run_pymol_video.py --pdb topology.pdb --traj traj.dcd --annotations semantic_store/annotations.json --out outdir --fps 30`

## 4) FastAPI client for raw MP4 endpoint

- Client: `clients/animate_mp4_client.py`
- Calls `GET /chronosfold/animate_mp4_raw` and saves MP4.

Example:
- `python VIDEOSUITE/clients/animate_mp4_client.py --host http://localhost:8000 --out demo.mp4 --frames 240 --fps 24`

## 4) Alex Chen breakthrough visualization

- Wrapper: `alex_chen/run_alex_chen_suite.py`
- Output: 4 GIF files (default in `alex_chen_suite/`)
- Internals: calls `scripts/chronos_video_pipeline.py` with `--preset alex_chen`

Generated GIFs:
- `chronos_novelty_predictions.gif`: Novel conformations with red star markers
- `chronos_density_map.gif`: Density-based conformational coloring
- `chronos_temporal_evolution.gif`: Time-based trajectory evolution
- `chronos_gpt.gif`: Original ChronosGPT breakthrough visualization

Example usage:
- Generate for SARS-CoV-2 Spike (default):
  - `python VIDEOSUITE/alex_chen/run_alex_chen_suite.py`
- Custom protein:
  - `python VIDEOSUITE/alex_chen/run_alex_chen_suite.py --pdb 1hvr.pdb --out custom_output/`

## 5) Umbrella sampling diagnostics (WNK)

- Wrapper: `umbrella/run_umbrella_visualization.py`
- Output: `umbrella_results/wnk_videos/umbrella_diagnostics.png` (+ optional MP4 if `--animation`)
- Internals: delegates to `Chronosfold.umbrella_suite` pipeline and visualization helpers.

Example usage:
- Quick pilot with synthetic fallback:
  - `python VIDEOSUITE/umbrella/run_umbrella_visualization.py`
- Forced animation output:
  - `python VIDEOSUITE/umbrella/run_umbrella_visualization.py --animation --frames 180 --fps 30`

Notes:
- If OpenMM/OpenMMTools are installed the wrapper will attempt picosecond umbrella sampling on `Chronosfold/WNK/5DRB.pdb` (distance CV between Lys245 and Glu292 Cα atoms).
- Without OpenMM it produces synthetic overlapping windows so the visualization pipeline and WHAM post-processing remain testable.

## Notes & Troubleshooting

- ffmpeg missing: install it or use GIF fallbacks in the MD demo.
- OpenMM/MDTraj not installed: run only the 2D pipeline or install the MD deps.
- Milvus optional: the MD demo continues without Milvus if env vars are absent.
- Alex Chen suite: requires the same dependencies as MD demo (mdtraj, scikit-learn)
- PyMOL headless: use `pymol -cq`; set `PYMOL_EXEC` and `FFMPEG_EXEC` env vars if needed.
