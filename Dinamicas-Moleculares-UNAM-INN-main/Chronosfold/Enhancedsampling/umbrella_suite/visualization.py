"""Visualization utilities for umbrella sampling diagnostics."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from .analysis import UmbrellaWindow

sns.set_context("talk")
sns.set_style("whitegrid")


def _interpolated_heatmap(windows: List[UmbrellaWindow]) -> np.ndarray:
    max_bins = max(window.histogram_counts.size for window in windows)
    grid = []
    reference = windows[0].histogram_edges
    centers_ref = 0.5 * (reference[:-1] + reference[1:])
    for window in windows:
        centers = 0.5 * (window.histogram_edges[:-1] + window.histogram_edges[1:])
        grid.append(np.interp(centers_ref, centers, window.histogram_counts, left=0.0, right=0.0))
    return np.asarray(grid)


def plot_umbrella_diagnostics(
    windows: List[UmbrellaWindow],
    pmf_df,
    *,
    title_suffix: str = "",
    max_windows: Optional[int] = None,
    save_path: Optional[Path] = None,
) -> None:
    if not windows:
        raise ValueError("No umbrella windows provided")

    subset = windows
    if max_windows and len(windows) > max_windows:
        stride = max(1, len(windows) // max_windows)
        subset = windows[::stride]

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # 1) Heatmap of histogram overlap
    ax = axes[0, 0]
    heatmap = _interpolated_heatmap(subset)
    sns.heatmap(heatmap, cmap="viridis", ax=ax, cbar_kws={"label": "P(CV)"})
    ax.set_ylabel("Window index")
    ax.set_xlabel("CV")
    ax.set_title(f"Histogram overlap {title_suffix}")

    # 2) Mean vs target
    ax = axes[0, 1]
    centers = [w.center for w in subset]
    means = [w.mean_cv for w in subset]
    stds = [w.std_cv for w in subset]
    ax.errorbar(centers, means, yerr=stds, fmt="o", ecolor="gray", capsize=4)
    ax.plot(centers, centers, linestyle="--", color="black", alpha=0.4)
    ax.set_xlabel("Target CV")
    ax.set_ylabel("Observed CV")
    ax.set_title("Window tracking")

    # 3) PMF curve
    ax = axes[1, 0]
    ax.plot(pmf_df["cv"], pmf_df["pmf"], color="#1f77b4", linewidth=2)
    if "uncertainty" in pmf_df:
        ax.fill_between(
            pmf_df["cv"],
            pmf_df["pmf"] - pmf_df["uncertainty"],
            pmf_df["pmf"] + pmf_df["uncertainty"],
            alpha=0.3,
        )
    ax.set_xlabel("CV")
    ax.set_ylabel("Delta G (kcal/mol)")
    ax.set_title("Potential of mean force")

    # 4) CDF subset
    ax = axes[1, 1]
    selection = subset if len(subset) <= 8 else subset[:: max(1, len(subset) // 8)]
    for window in selection:
        sorted_cv = np.sort(window.cv_values)
        cumulative = np.linspace(0, 1, sorted_cv.size)
        ax.plot(sorted_cv, cumulative, label=f"xi={window.center:.1f}")
    ax.set_xlabel("CV")
    ax.set_ylabel("CDF")
    ax.set_title("Cumulative distributions")
    ax.legend(loc="lower right", fontsize=8)

    plt.tight_layout()
    if save_path:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300)
    else:
        plt.show()
    plt.close(fig)
