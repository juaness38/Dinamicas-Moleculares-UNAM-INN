"""Umbrella sampling data loading, PMF computation, and synthetic scaffolds."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd

try:
    import pymbar
    MBAR_AVAILABLE = True
except ImportError:  # pragma: no cover - optional dependency
    MBAR_AVAILABLE = False

K_BOLTZMANN_KCAL = 0.0019872041


@dataclass
class UmbrellaWindow:
    center: float
    force_constant: float
    cv_values: np.ndarray
    histogram_counts: np.ndarray
    histogram_edges: np.ndarray

    @property
    def mean_cv(self) -> float:
        return float(np.mean(self.cv_values))

    @property
    def std_cv(self) -> float:
        return float(np.std(self.cv_values))


def _histogram_from_series(values: np.ndarray, bins: int = 80) -> Tuple[np.ndarray, np.ndarray]:
    counts, edges = np.histogram(values, bins=bins, density=True)
    return counts, edges


def generate_synthetic_windows(
    centers: Iterable[float],
    force_constant: float = 12.0,
    n_samples: int = 4000,
    noise: float = 0.6,
) -> List[UmbrellaWindow]:
    windows: List[UmbrellaWindow] = []
    for center in centers:
        sigma = max(0.15, math.sqrt(1.0 / max(force_constant, 1e-6)))
        samples = np.random.normal(loc=center, scale=sigma + noise, size=n_samples)
        counts, edges = _histogram_from_series(samples)
        windows.append(
            UmbrellaWindow(
                center=center,
                force_constant=force_constant,
                cv_values=samples,
                histogram_counts=counts,
                histogram_edges=edges,
            )
        )
    return windows


def load_umbrella_dataset(results_dir: Path) -> Tuple[List[UmbrellaWindow], Dict]:
    """Load umbrella sampling results previously exported by the pipeline."""
    windows: List[UmbrellaWindow] = []
    metadata: Dict = {}

    if not results_dir.exists():
        return windows, metadata

    metadata_path = results_dir / "umbrella_metadata.json"
    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))

    for hist_file in sorted(results_dir.glob("cv_histogram_center_*.dat")):
        try:
            center = float(hist_file.stem.split("_")[-1])
        except ValueError:
            continue
        hist_data = np.loadtxt(hist_file)
        if hist_data.ndim != 2 or hist_data.shape[1] < 2:
            continue

        counts = hist_data[:, 1]
        step = hist_data[1, 0] - hist_data[0, 0] if hist_data.shape[0] > 1 else 0.1
        edges = np.concatenate([hist_data[:, 0], [hist_data[-1, 0] + step]])

        series_file = results_dir / f"cv_timeseries_center_{center:.2f}.dat"
        if series_file.exists():
            cv_values = np.loadtxt(series_file)
            if np.isscalar(cv_values):
                cv_values = np.array([cv_values], dtype=float)
        else:
            cv_values = np.repeat(center, counts.size)

        force_constant = metadata.get("force_constant", 12.0)
        windows.append(
            UmbrellaWindow(
                center=center,
                force_constant=force_constant,
                cv_values=np.asarray(cv_values, dtype=float),
                histogram_counts=np.asarray(counts, dtype=float),
                histogram_edges=np.asarray(edges, dtype=float),
            )
        )

    return windows, metadata


def _combined_bins(windows: List[UmbrellaWindow], bins: int = 300) -> np.ndarray:
    min_cv = min(window.cv_values.min() for window in windows)
    max_cv = max(window.cv_values.max() for window in windows)
    # Return bin edges array of length (bins + 1) so that number of bins == bins
    return np.linspace(min_cv, max_cv, bins + 1)


def wham_iterative(
    windows: List[UmbrellaWindow],
    temperature: float = 300.0,
    bins: int = 300,
    max_iter: int = 10000,
    tolerance: float = 1e-6,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    WHAM iterativo clásico (Kumar et al. 1992)
    
    Algoritmo:
    1. P_unbiased(ξ) ∝ Σ_i n_i(ξ) / Σ_j N_j exp[-β(V_j(ξ) - F_j)]
    2. F_k = -β^-1 ln[ Σ_ξ P_unbiased(ξ) exp(-βV_k(ξ)) ]
    3. Iterar hasta convergencia de F_k
    
    Returns:
        cv_values: Array de coordenadas colectivas (bin centers)
        pmf: Potencial de fuerza media sin sesgo
        uncertainty: Estimación de incertidumbre (std entre iteraciones finales)
    """
    if not windows:
        return np.array([]), np.array([]), np.array([])
    
    # Setup
    bin_edges = _combined_bins(windows, bins)
    cv_values = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    n_bins = cv_values.size
    n_windows = len(windows)
    beta = 1.0 / (K_BOLTZMANN_KCAL * temperature)
    
    # Construir histogramas n_i(ξ)
    histograms = np.zeros((n_windows, n_bins))
    N_k = np.zeros(n_windows)
    for k, window in enumerate(windows):
        hist, _ = np.histogram(window.cv_values, bins=bin_edges, density=False)
        histograms[k, :] = hist
        N_k[k] = window.cv_values.size
    
    # Matriz de bias potentials V_k(ξ)
    V_kn = np.zeros((n_windows, n_bins))
    for k, window in enumerate(windows):
        kappa = window.force_constant
        center = window.center
        V_kn[k, :] = 0.5 * kappa * (cv_values - center) ** 2
    
    # Inicializar free energies F_k = 0
    F_k = np.zeros(n_windows)
    F_k_history = []
    # Bin width for numerical integrals
    delta = cv_values[1] - cv_values[0] if n_bins > 1 else 1.0
    
    # Iteración WHAM
    for iteration in range(max_iter):
        # Calcular P_unbiased(ξ) (ecuación WHAM principal)
        # numerator = Σ_i n_i(ξ) (raw counts)
        numerator = np.sum(histograms, axis=0)

        denominator = np.zeros(n_bins)
        for k in range(n_windows):
            denominator += N_k[k] * np.exp(-beta * (V_kn[k, :] - F_k[k]))

        # Evitar división por cero
        denominator[denominator <= 0] = 1e-30
        P_unbiased = numerator / denominator

        # Normalizar P_unbiased: ∑ P_unbiased(ξ) = 1
        norm = np.sum(P_unbiased)
        if norm <= 0:
            norm = 1e-30
        P_unbiased /= norm

        # Actualizar F_k (free energies)
        F_k_new = np.zeros(n_windows)
        for k in range(n_windows):
            # Discrete integral: ∫ P(ξ) exp(-βV_k(ξ)) dξ ≈ Σ P(ξ_i) exp(-βV_k(ξ_i))
            # P_unbiased is normalized as ΣP=1 (discrete probability), so no Δξ needed
            integrand = P_unbiased * np.exp(-beta * V_kn[k, :])
            F_k_new[k] = -(1.0 / beta) * np.log(np.sum(integrand) + 1e-30)

        # Shift para estabilidad numérica (F_0 = 0)
        F_k_new -= F_k_new[0]

        # Convergencia
        max_change = np.max(np.abs(F_k_new - F_k))
        F_k = F_k_new
        F_k_history.append(F_k.copy())

        if max_change < tolerance:
            break
    
    # PMF = -kT ln P_unbiased
    P_unbiased[P_unbiased <= 0] = 1e-30
    pmf = -(1.0 / beta) * np.log(P_unbiased)
    
    # Normalize to minimum
    pmf -= pmf.min()
    
    # Estimación de incertidumbre (std de últimas 10 iteraciones)
    if len(F_k_history) > 10:
        recent_pmfs = []
        for old_F_k in F_k_history[-10:]:
            denominator = np.zeros(n_bins)
            for k in range(n_windows):
                denominator += N_k[k] * np.exp(-beta * (V_kn[k, :] - old_F_k[k]))
            P_temp = numerator / (denominator + 1e-30)
            # Normalize
            norm_temp = np.sum(P_temp)
            if norm_temp <= 0:
                norm_temp = 1e-30
            P_temp /= norm_temp
            recent_pmfs.append(-(1.0 / beta) * np.log(P_temp + 1e-30))
        uncertainty = np.std(recent_pmfs, axis=0)
    else:
        uncertainty = np.full(n_bins, 0.5)  # Default uncertainty
    
    return cv_values, pmf, uncertainty


def compute_pmf(
    windows: List[UmbrellaWindow],
    temperature: float = 300.0,
    method: str = "mbar",
    bins: int = 300,
    **kwargs
) -> pd.DataFrame:
    """
    Compute PMF from umbrella sampling windows.
    
    Args:
        windows: Lista de ventanas de umbrella sampling
        temperature: Temperatura en Kelvin
        method: Método de análisis
            - 'mbar': MBAR (Multistate Bennett Acceptance Ratio) [DEFAULT]
            - 'wham': WHAM iterativo clásico (Kumar 1992)
            - 'histogram': Simple histogram combination (no reweighting)
        bins: Número de bins para histogramas
        **kwargs: Argumentos adicionales (max_iter, tolerance para WHAM)
    
    Returns:
        DataFrame con columnas: cv, pmf, uncertainty
    
    Examples:
        >>> # Usando MBAR (más preciso, requiere pymbar)
        >>> pmf_df = compute_pmf(windows, method='mbar')
        
        >>> # Usando WHAM clásico (educativo, transparente)
        >>> pmf_df = compute_pmf(windows, method='wham', max_iter=5000)
        
        >>> # Usando histogram simple (rápido, menos preciso)
        >>> pmf_df = compute_pmf(windows, method='histogram')
    """
    if not windows:
        return pd.DataFrame(columns=["cv", "pmf", "uncertainty"])

    bin_edges = _combined_bins(windows, bins)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    beta = 1.0 / (K_BOLTZMANN_KCAL * temperature)

    # MÉTODO 1: MBAR (default, más preciso)
    if method == "mbar":
        if not MBAR_AVAILABLE:
            raise ImportError(
                "pymbar not available. Install with: pip install pymbar\n"
                "Or use method='wham' or method='histogram'"
            )
        if len(windows) < 2:
            raise ValueError("MBAR requires at least 2 windows")
        
        N_k = np.array([window.cv_values.size for window in windows])
        K = len(windows)
        u_kn = np.zeros((K, bin_centers.size))
        for k, window in enumerate(windows):
            kappa = window.force_constant
            for idx, center in enumerate(bin_centers):
                u_kn[k, idx] = beta * 0.5 * kappa * (center - window.center) ** 2
        mbar = pymbar.MBAR(u_kn, N_k)
        f_i, df_i = mbar.computePMF(u_kn, bin_centers, nbins=bin_centers.size)
        pmf = f_i - f_i.min()
        uncertainty = df_i

    # MÉTODO 2: WHAM ITERATIVO (Kumar 1992)
    elif method == "wham":
        max_iter = kwargs.get("max_iter", 10000)
        tolerance = kwargs.get("tolerance", 1e-6)
        cv_values, pmf, uncertainty = wham_iterative(
            windows, temperature, bins, max_iter, tolerance
        )
        bin_centers = cv_values

    # MÉTODO 3: HISTOGRAM SIMPLE (fallback)
    elif method == "histogram":
        combined = np.zeros(bin_centers.size)
        for window in windows:
            hist, _ = np.histogram(window.cv_values, bins=bin_edges, density=False)
            combined += hist
        combined[combined <= 0] = 1e-10
        pmf = -np.log(combined)
        pmf -= pmf.min()
        uncertainty = np.full_like(pmf, fill_value=np.std(pmf) * 0.1 if pmf.size else 0.0)
    
    else:
        raise ValueError(f"Unknown method '{method}'. Use 'mbar', 'wham', or 'histogram'")

    return pd.DataFrame({
        "cv": bin_centers,
        "pmf": pmf,
        "uncertainty": uncertainty,
    })
