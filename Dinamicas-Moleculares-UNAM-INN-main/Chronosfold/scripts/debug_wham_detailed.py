# Detailed debug of WHAM algorithm step-by-step
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

# Generate 3 simple windows for easier debugging
def generate_simple_windows():
    """Generate 3 harmonic oscillator windows."""
    k_true = 10.0
    x0 = 12.0
    k_umbrella = 15.0
    temperature = 300.0
    
    kB = K_BOLTZMANN_KCAL
    beta = 1.0 / (kB * temperature)
    
    window_centers = [10.0, 12.0, 14.0]
    windows = []
    
    for xi in window_centers:
        k_eff = k_true + k_umbrella
        mean_eff = (k_true * x0 + k_umbrella * xi) / k_eff
        sigma_eff = np.sqrt(1.0 / (beta * k_eff))
        samples = np.random.normal(loc=mean_eff, scale=sigma_eff, size=5000)
        hist, edges = np.histogram(samples, bins=50, density=False)
        windows.append(UmbrellaWindow(
            center=xi,
            force_constant=k_umbrella,
            cv_values=samples,
            histogram_counts=hist,
            histogram_edges=edges,
        ))
    
    return windows

windows = generate_simple_windows()
print(f"Generated {len(windows)} windows")
for i, w in enumerate(windows):
    print(f"Window {i}: center={w.center}, N={w.cv_values.size}, mean_cv={w.mean_cv:.3f}, std_cv={w.std_cv:.3f}")

# Manual WHAM step-by-step
temperature = 300.0
bins = 100
beta = 1.0 / (K_BOLTZMANN_KCAL * temperature)

# Bin edges
min_cv = min(w.cv_values.min() for w in windows)
max_cv = max(w.cv_values.max() for w in windows)
bin_edges = np.linspace(min_cv, max_cv, bins + 1)
cv_values = 0.5 * (bin_edges[:-1] + bin_edges[1:])
n_bins = cv_values.size
n_windows = len(windows)

# Histograms
histograms = np.zeros((n_windows, n_bins))
N_k = np.zeros(n_windows)
for k, window in enumerate(windows):
    hist, _ = np.histogram(window.cv_values, bins=bin_edges, density=False)
    histograms[k, :] = hist
    N_k[k] = window.cv_values.size
    print(f"Window {k}: N_k={N_k[k]}, total_hist={np.sum(hist)}")

# Bias potentials
V_kn = np.zeros((n_windows, n_bins))
for k, window in enumerate(windows):
    kappa = window.force_constant
    center = window.center
    V_kn[k, :] = 0.5 * kappa * (cv_values - center) ** 2

print(f"\nBias potentials at cv={cv_values[n_bins//2]:.3f} (middle bin):")
for k in range(n_windows):
    print(f"Window {k}: V_k = {V_kn[k, n_bins//2]:.3f} kcal/mol")

# WHAM iteration
F_k = np.zeros(n_windows)
tolerance = 1e-6
max_iter = 1000

for iteration in range(max_iter):
    numerator = np.sum(histograms, axis=0)
    denominator = np.zeros(n_bins)
    for k in range(n_windows):
        denominator += N_k[k] * np.exp(-beta * (V_kn[k, :] - F_k[k]))
    
    denominator[denominator <= 0] = 1e-30
    P_unbiased = numerator / denominator
    norm = np.sum(P_unbiased)
    if norm <= 0:
        norm = 1e-30
    P_unbiased /= norm
    
    # Update F_k
    F_k_new = np.zeros(n_windows)
    for k in range(n_windows):
        integrand = P_unbiased * np.exp(-beta * V_kn[k, :])
        F_k_new[k] = -(1.0 / beta) * np.log(np.sum(integrand) + 1e-30)
    
    # Shift
    F_k_new -= F_k_new[0]
    
    # Check convergence
    max_change = np.max(np.abs(F_k_new - F_k))
    if iteration < 10 or iteration % 100 == 0:
        print(f"Iteration {iteration}: max_change={max_change:.2e}, F_k={F_k_new}")
    
    F_k = F_k_new
    
    if max_change < tolerance:
        print(f"Converged at iteration {iteration}")
        break

# Final PMF
P_unbiased[P_unbiased <= 0] = 1e-30
pmf = -(1.0 / beta) * np.log(P_unbiased)
pmf -= np.nanmin(pmf)

# Analytic PMF
k_true = 10.0
x0 = 12.0
pmf_true = 0.5 * k_true * (cv_values - x0) ** 2
pmf_true -= pmf_true.min()

# Compare around x0=12
mask = (cv_values >= 10.5) & (cv_values <= 13.5)
print(f"\nComparison around x0=12:")
print(f"{'CV':>8} {'PMF_calc':>10} {'PMF_true':>10} {'Diff':>10}")
for c, p, pt in zip(cv_values[mask], pmf[mask], pmf_true[mask]):
    print(f"{c:8.3f} {p:10.4f} {pt:10.4f} {(p-pt):10.4f}")

rmsd = np.sqrt(np.mean((pmf[mask] - pmf_true[mask])**2))
print(f"\nRMSD: {rmsd:.4f} kcal/mol")
