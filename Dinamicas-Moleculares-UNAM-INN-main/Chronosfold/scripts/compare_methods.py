# Test simple histogram combination vs WHAM vs analytic
import importlib.util
import sys
from pathlib import Path
import numpy as np

analysis_path = Path(__file__).resolve().parents[1] / "umbrella_suite" / "analysis.py"
module_name = "Chronosfold.umbrella_suite.analysis"
spec = importlib.util.spec_from_file_location(module_name, str(analysis_path))
analysis = importlib.util.module_from_spec(spec)
sys.modules[module_name] = analysis
spec.loader.exec_module(analysis)

UmbrellaWindow = analysis.UmbrellaWindow
K_BOLTZMANN_KCAL = analysis.K_BOLTZMANN_KCAL

# Generate windows
k_true = 10.0
x0 = 12.0
k_umbrella = 15.0
temperature = 300.0

kB = K_BOLTZMANN_KCAL
beta = 1.0 / (kB * temperature)

window_centers = np.linspace(8.0, 16.0, 9)
windows = []

for xi in window_centers:
    k_eff = k_true + k_umbrella
    mean_eff = (k_true * x0 + k_umbrella * xi) / k_eff
    sigma_eff = np.sqrt(1.0 / (beta * k_eff))
    samples = np.random.normal(loc=mean_eff, scale=sigma_eff, size=10000)
    hist, edges = np.histogram(samples, bins=80, density=False)
    windows.append(UmbrellaWindow(
        center=xi,
        force_constant=k_umbrella,
        cv_values=samples,
        histogram_counts=hist,
        histogram_edges=edges,
    ))

# Define common bins
bins = 200
min_cv = min(w.cv_values.min() for w in windows)
max_cv = max(w.cv_values.max() for w in windows)
bin_edges = np.linspace(min_cv, max_cv, bins + 1)
cv_values = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# Method 1: Simple histogram combination
combined_hist = np.zeros(bins)
for window in windows:
    hist, _ = np.histogram(window.cv_values, bins=bin_edges, density=False)
    combined_hist += hist

combined_hist[combined_hist <= 0] = 1e-10
pmf_histogram = -(1.0 / beta) * np.log(combined_hist)
pmf_histogram -= pmf_histogram.min()

# Method 2: WHAM (from compute_pmf)
compute_pmf = analysis.compute_pmf
pmf_df_wham = compute_pmf(windows, temperature=temperature, method='wham', bins=bins, max_iter=5000, tolerance=1e-8)
pmf_wham = pmf_df_wham['pmf'].values

# Analytic PMF
pmf_true = 0.5 * k_true * (cv_values - x0) ** 2
pmf_true -= pmf_true.min()

# Compare
mask = (cv_values >= 9.0) & (cv_values <= 15.0)
rmsd_hist = np.sqrt(np.mean((pmf_histogram[mask] - pmf_true[mask])**2))
rmsd_wham = np.sqrt(np.mean((pmf_wham[mask] - pmf_true[mask])**2))

print(f"RMSD histogram vs analytic: {rmsd_hist:.4f} kcal/mol")
print(f"RMSD WHAM vs analytic: {rmsd_wham:.4f} kcal/mol")

# Show first 10 comparisons
print('\nFirst 10 (cv, histogram, WHAM, true):')
for i in range(10):
    print(f"{cv_values[i]:.3f}, {pmf_histogram[i]:.4f}, {pmf_wham[i]:.4f}, {pmf_true[i]:.4f}")

# Check around minimum
center_mask = (cv_values >= 11.5) & (cv_values <= 12.5)
print(f'\nAround x0=12 (should all be near 0):')
print(f"Histogram: min={pmf_histogram[center_mask].min():.4f}, max={pmf_histogram[center_mask].max():.4f}")
print(f"WHAM:      min={pmf_wham[center_mask].min():.4f}, max={pmf_wham[center_mask].max():.4f}")
print(f"True:      min={pmf_true[center_mask].min():.4f}, max={pmf_true[center_mask].max():.4f}")
