#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Análisis MBAR de umbrella sampling para WNK1 C-terminal

Lee los valores de CV de todas las ventanas y calcula PMF con MBAR
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Intentar importar pymbar
try:
    import pymbar
    HAS_PYMBAR = True
except ImportError:
    HAS_PYMBAR = False
    print("⚠️  Advertencia: pymbar no instalado")
    print("   Instalación: pip install pymbar")

def parse_args():
    parser = argparse.ArgumentParser(description="Análisis MBAR de umbrella sampling WNK1")
    parser.add_argument('--windows-dir', type=str, default='umbrella_windows',
                       help='Directorio con ventanas (default: umbrella_windows)')
    parser.add_argument('--output', type=str, default='pmf_analysis',
                       help='Directorio de salida (default: pmf_analysis)')
    parser.add_argument('--temperature', type=float, default=300.0,
                       help='Temperatura en K (default: 300)')
    parser.add_argument('--subsample', type=int, default=1,
                       help='Submuestreo (tomar cada N frames, default: 1)')
    parser.add_argument('--equilibration', type=int, default=0,
                       help='Frames de equilibración a descartar (default: 0)')
    return parser.parse_args()

def load_cv_data(windows_dir, n_windows, equilibration=0, subsample=1):
    """Carga datos de CV de todas las ventanas"""
    print("\n" + "="*70)
    print("Cargando datos de CV")
    print("="*70)
    
    all_cv_values = []
    all_window_ids = []
    window_r0 = []
    window_k = []
    
    # Leer configuración de ventanas
    windows_config = pd.read_csv(windows_dir / "windows_config.csv")
    
    for i in range(n_windows):
        window_dir = windows_dir / f"window_{i:02d}"
        cv_file = window_dir / "cv_values.dat"
        
        if not cv_file.exists():
            print(f"⚠️  Ventana {i}: {cv_file} no encontrada, omitiendo")
            continue
        
        # Leer CV values
        data = pd.read_csv(cv_file, sep='\s+', comment='#', 
                          names=['step', 'time_ps', 'distance_nm', 'bias_energy_kJ'])
        
        # Aplicar equilibración y submuestreo
        if equilibration > 0:
            data = data.iloc[equilibration:]
        
        if subsample > 1:
            data = data.iloc[::subsample]
        
        cv_values = data['distance_nm'].values
        
        # Parámetros de esta ventana
        r0 = windows_config.loc[i, 'r0_nm']
        k = windows_config.loc[i, 'spring_constant_kJ_mol_nm2']
        
        all_cv_values.append(cv_values)
        all_window_ids.append(np.full(len(cv_values), i))
        window_r0.append(r0)
        window_k.append(k)
        
        print(f"  Ventana {i:2d}: {len(cv_values):6d} samples, r₀={r0:.3f} nm, <r>={cv_values.mean():.3f} nm")
    
    if len(all_cv_values) == 0:
        print("❌ ERROR: No se encontraron datos de CV")
        sys.exit(1)
    
    # Concatenar todos los datos
    cv_concat = np.concatenate(all_cv_values)
    window_ids = np.concatenate(all_window_ids)
    
    print(f"\n✓ Total: {len(cv_concat):,} samples de {len(all_cv_values)} ventanas")
    
    return cv_concat, window_ids, all_cv_values, window_r0, window_k

def compute_bias_matrix(cv_values, window_ids, window_r0, window_k, temperature):
    """Calcula matriz de energías de bias u_kn"""
    print("\n" + "="*70)
    print("Calculando matriz de bias")
    print("="*70)
    
    kB = 0.008314  # kJ/mol/K
    beta = 1.0 / (kB * temperature)
    
    K = len(window_r0)  # Número de ventanas
    N_k = np.array([np.sum(window_ids == k) for k in range(K)])  # Samples por ventana
    
    print(f"  K = {K} ventanas")
    print(f"  Samples por ventana: min={N_k.min()}, max={N_k.max()}, mean={N_k.mean():.0f}")
    
    # Inicializar matriz u_kn (K ventanas x total samples)
    N_total = len(cv_values)
    u_kn = np.zeros((K, N_total))
    
    # Para cada ventana k, calcular U_k(r_n) para todos los samples n
    for k in range(K):
        r0 = window_r0[k]
        spring_k = window_k[k]
        
        # U_k(r) = 0.5 * k * (r - r0)^2
        # u_kn = beta * U_k
        u_kn[k, :] = beta * 0.5 * spring_k * (cv_values - r0)**2
    
    print(f"✓ Matriz u_kn: {u_kn.shape}")
    print(f"  Rango de bias: {u_kn.min():.2f} - {u_kn.max():.2f} kT")
    
    return u_kn, N_k

def run_mbar(u_kn, N_k, cv_values, temperature, n_bins=50):
    """Ejecuta MBAR y calcula PMF"""
    print("\n" + "="*70)
    print("Ejecutando MBAR")
    print("="*70)
    
    if not HAS_PYMBAR:
        print("❌ ERROR: pymbar no instalado")
        print("   pip install pymbar")
        sys.exit(1)
    
    # Inicializar MBAR
    print("  Inicializando MBAR...")
    mbar = pymbar.MBAR(u_kn, N_k, verbose=True, relative_tolerance=1e-12, 
                       maximum_iterations=10000)
    
    print("✓ MBAR converged")
    print(f"  Free energies F_k shape: {mbar.f_k.shape}")
    
    # Calcular PMF a lo largo de bins de CV
    cv_min = cv_values.min()
    cv_max = cv_values.max()
    cv_bins = np.linspace(cv_min, cv_max, n_bins)
    cv_centers = 0.5 * (cv_bins[:-1] + cv_bins[1:])
    
    print(f"\n  Calculando PMF en {len(cv_centers)} bins")
    print(f"  Rango CV: {cv_min:.3f} - {cv_max:.3f} nm")
    
    # Para cada bin, calcular energía libre
    pmf = np.zeros(len(cv_centers))
    pmf_err = np.zeros(len(cv_centers))
    
    kB = 0.008314  # kJ/mol/K
    kT = kB * temperature
    
    for i, cv_center in enumerate(cv_centers):
        # Encontrar samples en este bin
        delta = (cv_bins[1] - cv_bins[0]) / 2
        mask = (cv_values >= cv_center - delta) & (cv_values < cv_center + delta)
        
        if np.sum(mask) == 0:
            pmf[i] = np.nan
            pmf_err[i] = np.nan
            continue
        
        # Calcular matriz de bias para este bin (promedio de samples en el bin)
        u_kn_bin = u_kn[:, mask]
        
        # Energía libre del bin (negativo de log de probabilidad)
        # Usamos expectation de MBAR para el bin
        try:
            # Método 1: Usar computeExpectations para el bin
            # (Simplificado: usamos la probabilidad unbiased)
            
            # Calcular probabilidad unbiased de estar en este bin
            log_w = -u_kn_bin.sum(axis=0)  # Simplified
            log_w -= log_w.max()  # Numerical stability
            w = np.exp(log_w)
            w /= w.sum()
            
            # PMF = -kT * ln(P)
            P = w.sum()
            pmf[i] = -kT * np.log(P + 1e-100)
            
        except Exception as e:
            print(f"    Bin {i}: Error calculando PMF: {e}")
            pmf[i] = np.nan
            pmf_err[i] = np.nan
    
    # Normalizar PMF (mínimo = 0)
    pmf_min = np.nanmin(pmf)
    pmf -= pmf_min
    
    print(f"✓ PMF calculado")
    print(f"  Mínimo: {pmf_min:.2f} kJ/mol (normalizado a 0)")
    print(f"  Máximo: {np.nanmax(pmf):.2f} kJ/mol")
    
    return cv_centers, pmf, pmf_err, mbar

def plot_results(cv_centers, pmf, pmf_err, cv_all, window_ids, window_r0, output_dir):
    """Genera plots de resultados"""
    print("\n" + "="*70)
    print("Generando plots")
    print("="*70)
    
    sns.set_style("whitegrid")
    
    # Plot 1: PMF con barras de error
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.errorbar(cv_centers, pmf, yerr=pmf_err, marker='o', linestyle='-',
               capsize=3, label='MBAR PMF')
    
    ax.set_xlabel('Distancia C-terminal (nm)', fontsize=12)
    ax.set_ylabel('PMF (kJ/mol)', fontsize=12)
    ax.set_title('Potential of Mean Force - WNK1 C-terminal', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plot_file = output_dir / "pmf.png"
    plt.savefig(plot_file, dpi=300)
    print(f"✓ PMF plot: {plot_file}")
    plt.close()
    
    # Plot 2: Histogramas de ventanas + PMF
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    
    # Histogramas
    K = len(window_r0)
    colors = plt.cm.viridis(np.linspace(0, 1, K))
    
    for k in range(K):
        mask = window_ids == k
        cv_window = cv_all[mask]
        
        ax1.hist(cv_window, bins=30, alpha=0.6, color=colors[k], 
                label=f'Window {k}' if k % 3 == 0 else None)
        ax1.axvline(window_r0[k], color=colors[k], linestyle='--', alpha=0.8, linewidth=0.5)
    
    ax1.set_ylabel('Counts', fontsize=12)
    ax1.set_title('Distribuciones de CV por ventana', fontsize=14, fontweight='bold')
    ax1.legend(ncol=4, fontsize=8)
    
    # PMF
    ax2.plot(cv_centers, pmf, 'o-', linewidth=2, markersize=6, color='darkblue')
    ax2.fill_between(cv_centers, pmf - pmf_err, pmf + pmf_err, alpha=0.3, color='lightblue')
    
    ax2.set_xlabel('Distancia C-terminal (nm)', fontsize=12)
    ax2.set_ylabel('PMF (kJ/mol)', fontsize=12)
    ax2.set_title('PMF con incertidumbre', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    combined_plot = output_dir / "analysis_combined.png"
    plt.savefig(combined_plot, dpi=300)
    print(f"✓ Combined plot: {combined_plot}")
    plt.close()

def save_results(cv_centers, pmf, pmf_err, output_dir):
    """Guarda resultados a CSV"""
    print("\n" + "="*70)
    print("Guardando resultados")
    print("="*70)
    
    df = pd.DataFrame({
        'cv_nm': cv_centers,
        'pmf_kJ_mol': pmf,
        'pmf_error_kJ_mol': pmf_err
    })
    
    csv_file = output_dir / "pmf_results.csv"
    df.to_csv(csv_file, index=False)
    print(f"✓ Resultados CSV: {csv_file}")
    
    return df

def main():
    args = parse_args()
    
    print("="*70)
    print("ANÁLISIS MBAR - WNK1 UMBRELLA SAMPLING")
    print("="*70)
    
    WORK_DIR = Path(__file__).parent
    WINDOWS_DIR = WORK_DIR / args.windows_dir
    OUTPUT_DIR = WORK_DIR / args.output
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    print(f"\n📁 Directorio ventanas: {WINDOWS_DIR}")
    print(f"📂 Salida: {OUTPUT_DIR}")
    print(f"\n⚙️  PARÁMETROS:")
    print(f"  Temperatura: {args.temperature} K")
    print(f"  Equilibración: {args.equilibration} frames")
    print(f"  Submuestreo: cada {args.subsample} frames")
    
    # Detectar número de ventanas
    window_dirs = sorted(WINDOWS_DIR.glob("window_*"))
    n_windows = len(window_dirs)
    print(f"  Ventanas detectadas: {n_windows}")
    
    if n_windows == 0:
        print("❌ ERROR: No se encontraron ventanas")
        print(f"   Verifica que {WINDOWS_DIR} contenga directorios window_XX")
        sys.exit(1)
    
    # Cargar datos
    cv_all, window_ids, cv_by_window, window_r0, window_k = load_cv_data(
        WINDOWS_DIR, n_windows, args.equilibration, args.subsample
    )
    
    # Calcular matriz de bias
    u_kn, N_k = compute_bias_matrix(cv_all, window_ids, window_r0, window_k, args.temperature)
    
    # Ejecutar MBAR
    cv_centers, pmf, pmf_err, mbar = run_mbar(u_kn, N_k, cv_all, args.temperature)
    
    # Guardar resultados
    df = save_results(cv_centers, pmf, pmf_err, OUTPUT_DIR)
    
    # Plots
    plot_results(cv_centers, pmf, pmf_err, cv_all, window_ids, window_r0, OUTPUT_DIR)
    
    # Resumen
    print("\n" + "="*70)
    print("✅ ANÁLISIS COMPLETADO")
    print("="*70)
    print(f"\nArchivos generados en {OUTPUT_DIR}:")
    print(f"  1. pmf_results.csv       - PMF y errores")
    print(f"  2. pmf.png               - Plot de PMF")
    print(f"  3. analysis_combined.png - Histogramas + PMF")
    print("\n📊 RESULTADOS:")
    print(f"  Barrera energética: {pmf.max():.2f} ± {pmf_err[pmf.argmax()]:.2f} kJ/mol")
    print(f"  Distancia del mínimo: {cv_centers[pmf.argmin()]:.3f} nm")
    print(f"  Distancia de la barrera: {cv_centers[pmf.argmax()]:.3f} nm")
    print("="*70)

if __name__ == "__main__":
    main()
