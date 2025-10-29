# Debug WHAM accuracy
import importlib.util
import sys
from pathlib import Path
import numpy as np

analysis_path = Path(__file__).resolve().parents[1] / "umbrella_suite" / "analysis.py"
module_name = "Chronosfold.umbrella_suite.analysis"
spec = importlib.util.spec_from_file_location(module_name, str(analysis_path))
analysis = importlib.util.module_from_spec(spec)
import types
sys.modules[module_name] = analysis
spec.loader.exec_module(analysis)

UmbrellaWindow = analysis.UmbrellaWindow
compute_pmf = analysis.compute_pmf
wham_iterative = analysis.wham_iterative

# generate synthetic windows (same as tests)
def generate_harmonic_oscillator_windows(k_true=10.0, x0=12.0, window_centers=None, k_umbrella=15.0, n_samples=5000, temperature=300.0):
    if window_centers is None:
        window_centers = np.linspace(8.0, 16.0, 9)
    kB = 0.0019872041
    beta = 1.0 / (kB * temperature)
    windows = []
    for xi in window_centers:
        k_eff = k_true + k_umbrella
        mean_eff = (k_true * x0 + k_umbrella * xi) / k_eff
        sigma_eff = np.sqrt(1.0 / (beta * k_eff))
        samples = np.random.normal(loc=mean_eff, scale=sigma_eff, size=n_samples)
        hist, edges = np.histogram(samples, bins=80, density=False)
        windows.append(UmbrellaWindow(center=xi, force_constant=k_umbrella, cv_values=samples, histogram_counts=hist, histogram_edges=edges))
    return windows

windows = generate_harmonic_oscillator_windows(n_samples=10000)
pmf_df = compute_pmf(windows, temperature=300.0, method='wham', bins=200)
cv = pmf_df['cv'].values
pmf_calc = pmf_df['pmf'].values
pmf_true = 0.5 * 10.0 * (cv - 12.0)**2
pmf_true -= pmf_true.min()

mask = (cv >= 8.5) & (cv <= 15.5)
rmsd = np.sqrt(np.mean((pmf_calc[mask] - pmf_true[mask])**2))
print(f"RMSD WHAM vs analytic: {rmsd:.4f} kcal/mol")
print('\nFirst 10 comparisons (cv, pmf_calc, pmf_true):')
for i in range(10):
    print(f"{cv[i]:.3f}, {pmf_calc[i]:.4f}, {pmf_true[i]:.4f}")

# Also compute PMF from wham_iterative directly
cv2, pmf2, unc2 = wham_iterative(windows, temperature=300.0, bins=200, max_iter=10000, tolerance=1e-8)
mask2 = (cv2 >= 8.5) & (cv2 <= 15.5)
pmf_true2 = 0.5*10.0*(cv2-12.0)**2
pmf_true2 -= pmf_true2.min()
rmsd2 = np.sqrt(np.mean((pmf2[mask2] - pmf_true2[mask2])**2))
print(f"RMSD wham_iterative vs analytic: {rmsd2:.4f}")

# Check center region around x0=12
center_mask = (cv2 >= 11.5) & (cv2 <= 12.5)
print(f'\nPMF values around x0=12 (should be minimum ~0):')
for i, (c, p, pt) in enumerate(zip(cv2[center_mask], pmf2[center_mask], pmf_true2[center_mask])):
    print(f"cv={c:.3f}, pmf_calc={p:.4f}, pmf_true={pt:.4f}, diff={p-pt:.4f}")
print('done')
