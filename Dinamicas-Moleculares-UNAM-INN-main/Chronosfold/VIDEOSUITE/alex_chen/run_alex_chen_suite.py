#!/usr/bin/env python3
"""
Alex Chen Breakthrough Visualization Pipeline
============================================

Official VIDEOSUITE wrapper for Alex Chen breakthrough visualization.
Generates 4 GIF files using ChronosFold → PCA → GIF pipeline.

Output GIFs:
- chronos_novelty_predictions.gif: Novel conformations with red star markers
- chronos_density_map.gif: Density-based conformational coloring  
- chronos_temporal_evolution.gif: Time-based trajectory evolution
- chronos_gpt.gif: Original ChronosGPT breakthrough visualization

Usage:
    python VIDEOSUITE/alex_chen/run_alex_chen_suite.py
    python VIDEOSUITE/alex_chen/run_alex_chen_suite.py --pdb 1hvr.pdb --out custom_dir/
"""

import argparse
import subprocess
import sys
from pathlib import Path

# Get the repository root (2 levels up from VIDEOSUITE/alex_chen/)
ROOT = Path(__file__).resolve().parents[2]
PY = sys.executable


def generate_alex_chen_suite(pdb_file: str, output_dir: Path, verbose: bool = False):
    """Generate the 4 Alex Chen breakthrough GIFs"""
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if verbose:
        print(f"Alex Chen Breakthrough Suite")
        print(f"PDB: {pdb_file}")
        print(f"Output: {output_dir}")
        print("-" * 50)
    
    # Alex Chen breakthrough configurations
    configs = [
        ("chronos_novelty_predictions", "novelty"),
        ("chronos_density_map", "density"),
        ("chronos_temporal_evolution", "time"),
        ("chronos_gpt", "novelty")
    ]
    
    success_count = 0
    
    for name, color_by in configs:
        if verbose:
            print(f"Generating {name}...")
        
        output_path = output_dir / f"{name}.gif"
        
        # Use the official chronos_video_pipeline with alex_chen preset
        cmd = [
            PY, str(ROOT / "scripts" / "chronos_video_pipeline.py"),
            "--preset", "alex_chen",
            "--pdb", pdb_file,
            "--out", str(output_path),
            "--color-by", color_by
        ]
        
        try:
            result = subprocess.run(cmd, 
                                  capture_output=True, 
                                  text=True, 
                                  timeout=600,  # 10 minutes timeout
                                  cwd=str(ROOT))
            
            if result.returncode == 0 and output_path.exists():
                size_mb = output_path.stat().st_size / (1024 * 1024)
                if verbose:
                    print(f"  ✓ Success: {name} ({size_mb:.1f} MB)")
                success_count += 1
            else:
                if verbose:
                    print(f"  ✗ Failed: {name}")
                    if result.stderr:
                        print(f"    Error: {result.stderr[:200]}...")
        
        except subprocess.TimeoutExpired:
            if verbose:
                print(f"  ✗ Timeout: {name}")
        except Exception as e:
            if verbose:
                print(f"  ✗ Error: {name} - {e}")
    
    if verbose:
        print("-" * 50)
        print(f"Alex Chen Suite Complete: {success_count}/4 videos generated")
    
    return success_count


def main():
    parser = argparse.ArgumentParser(description="Alex Chen Breakthrough Visualization Suite")
    parser.add_argument("--pdb", default="6VXX.pdb", 
                       help="PDB file to process (default: 6VXX.pdb)")
    parser.add_argument("--out", type=Path, default=Path("alex_chen_suite"),
                       help="Output directory (default: alex_chen_suite)")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Verbose output")
    
    args = parser.parse_args()
    
    # Validate PDB file exists
    pdb_path = ROOT / args.pdb
    if not pdb_path.exists():
        print(f"Error: PDB file not found: {pdb_path}")
        return 1
    
    print("ChronosFold VIDEOSUITE - Alex Chen Breakthrough Pipeline")
    print("=" * 60)
    
    success_count = generate_alex_chen_suite(
        pdb_file=str(args.pdb),
        output_dir=args.out,
        verbose=args.verbose
    )
    
    if success_count == 4:
        print("\n🎉 Alex Chen Breakthrough Suite Complete!")
        print(f"📁 Output: {args.out.resolve()}")
        return 0
    else:
        print(f"\n⚠️ Partial success: {success_count}/4 videos generated")
        return 1


if __name__ == "__main__":
    sys.exit(main())
