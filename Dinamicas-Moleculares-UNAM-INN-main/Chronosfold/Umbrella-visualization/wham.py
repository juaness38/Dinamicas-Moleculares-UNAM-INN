"""CLI tool to compute PMF using umbrella_suite analysis helpers."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from Chronosfold.umbrella_suite.analysis import compute_pmf, load_umbrella_dataset
from Chronosfold.umbrella_suite.visualization import plot_umbrella_diagnostics


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute PMF from umbrella sampling outputs",
        epilog="Examples:\n"
               "  python wham.py results/ --method mbar\n"
               "  python wham.py results/ --method wham --max-iter 5000\n"
               "  python wham.py results/ --method histogram --bins 500\n",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "results",
        type=Path,
        default=Path("umbrella_results/wnk_pilot"),
        help="Directory containing umbrella outputs",
    )
    parser.add_argument(
        "--method",
        type=str,
        default="mbar",
        choices=["mbar", "wham", "histogram"],
        help="PMF calculation method: 'mbar' (default, most accurate), "
             "'wham' (classic iterative), 'histogram' (simple, fast)",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=300.0,
        help="Temperature in Kelvin (default: 300.0)",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=300,
        help="Number of histogram bins (default: 300)",
    )
    parser.add_argument(
        "--max-iter",
        type=int,
        default=10000,
        help="Max iterations for WHAM (default: 10000)",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1e-6,
        help="Convergence tolerance for WHAM (default: 1e-6)",
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=None,
        help="Diagnostics plot path (default: results/pmf_diagnostics.png)",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=None,
        help="CSV output path (default: results/pmf.csv)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    
    # Set default paths if not provided
    if args.csv is None:
        args.csv = args.results / f"pmf_{args.method}.csv"
    if args.plot is None:
        args.plot = args.results / f"pmf_diagnostics_{args.method}.png"
    
    # Load data
    print(f"Loading umbrella windows from {args.results}...")
    windows, metadata = load_umbrella_dataset(args.results)
    if not windows:
        raise SystemExit(f"No umbrella windows found in {args.results}")
    
    print(f"Found {len(windows)} windows")
    print(f"Method: {args.method.upper()}")
    print(f"Temperature: {args.temperature} K")
    
    # Compute PMF
    kwargs = {
        "bins": args.bins,
    }
    if args.method == "wham":
        kwargs["max_iter"] = args.max_iter
        kwargs["tolerance"] = args.tolerance
        print(f"WHAM parameters: max_iter={args.max_iter}, tolerance={args.tolerance}")
    
    print(f"\nComputing PMF using {args.method}...")
    pmf_df = compute_pmf(windows, temperature=args.temperature, method=args.method, **kwargs)
    
    # Save results
    if not pmf_df.empty:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        pmf_df.to_csv(args.csv, index=False)
        print(f"✓ PMF saved to {args.csv}")
        
        # Print summary statistics
        min_pmf = pmf_df["pmf"].min()
        max_pmf = pmf_df["pmf"].max()
        barrier = max_pmf - min_pmf
        print(f"\nPMF Statistics:")
        print(f"  Range: {pmf_df['cv'].min():.2f} - {pmf_df['cv'].max():.2f} Å")
        print(f"  Barrier: {barrier:.2f} kcal/mol")
        print(f"  Mean uncertainty: {pmf_df['uncertainty'].mean():.3f} kcal/mol")
    
    # Generate diagnostics plot
    print(f"\nGenerating diagnostics plot...")
    plot_umbrella_diagnostics(
        windows,
        pmf_df,
        title_suffix=f"({args.method.upper()})",
        save_path=args.plot
    )
    print(f"✓ Diagnostics saved to {args.plot}")
    print("\nDone!")


if __name__ == "__main__":
    main()
