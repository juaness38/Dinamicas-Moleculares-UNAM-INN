# ChronosFold “All-in” bundle (2025-08-15)

Contents
- chronos_analytics.mp4 + chronos_analytics.json (sidecar)
- chronodiffusion_test.json + chronodiffusion_test.md (skipped: PyTorch missing)
- dual_inference_balanced.json or dual_inference_balanced_info.json (status)

To complete pending items
1) Install PyTorch (CPU is fine) in either the conda env or the repo venv:
   - Conda env (af-openmm-311):
     conda run -n af-openmm-311 python -m pip install --index-url https://download.pytorch.org/whl/cpu torch --only-binary=:all: --disable-pip-version-check
   - Or repo venv:
     .venv\Scripts\python -m pip install --index-url https://download.pytorch.org/whl/cpu torch --only-binary=:all:
2) ChronoDiffusion test:
   $env:PYTHONPATH="chronosfold_scaffold/src"; conda run -n af-openmm-311 python scripts\chronodiffusion_test.py --epochs 6 --timesteps 80 --samples 8 --dim 32 --out videos\15\08\25\chronos_allin\chronodiffusion_test.json --markdown videos\15\08\25\chronos_allin\chronodiffusion_test.md
3) Dual inference JSON:
   $env:PYTHONPATH="chronosfold_scaffold/src"; .venv\Scripts\python scripts\dual_inference_cli.py --preset balanced --seq 0,1,2,3,2,1,0 --codebook-k 0 --samples 24 --protein HP35_1YRF > videos\15\08\25\chronos_allin\dual_inference_balanced.json
4) MD demo (OpenMM OK):
   conda run -n af-openmm-311 python VIDEOSUITE\md_demo\run_minicronos_demo.py --steps 120 --pred-steps 30 --gif videos\VIDEOSUITE\md_demo\minicronos_demo.gif
   (MiniChronoGPT is optional; script now falls back gracefully if torch is missing.)
