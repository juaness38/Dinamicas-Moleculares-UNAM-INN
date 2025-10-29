#!/usr/bin/env python3
"""
Extrae PMF desde archivo HILLS de metadinámica.
Compatible con well-tempered metadynamics.
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

def load_hills(hills_file):
    """
    Carga archivo HILLS de metadinámica.
    
    Formato esperado (PLUMED/OpenMM):
    # time(ps)  CV(nm)  height(kJ/mol)  sigma(nm)
    
    Retorna:
    --------
    times : ndarray
        Tiempos de depósito (ps)
    cv_centers : ndarray
        Valores de CV donde se depositaron gaussianos
    heights : ndarray
        Alturas de los gaussianos (kJ/mol)
    sigmas : ndarray
        Anchos de los gaussianos (nm)
    """
    
    print(f"Cargando {hills_file}...")
    
    try:
        data = np.loadtxt(hills_file, comments='#')
    except Exception as e:
        print(f"❌ Error al leer {hills_file}: {e}")
        sys.exit(1)
    
    if data.shape[1] < 4:
        print(f"❌ Formato incorrecto. Se esperan 4 columnas: time, CV, height, sigma")
        sys.exit(1)
    
    times = data[:, 0]
    cv_centers = data[:, 1]
    heights = data[:, 2]
    sigmas = data[:, 3]
    
    print(f"✅ {len(times)} gaussianos cargados")
    print(f"   Tiempo: {times[0]:.1f} - {times[-1]:.1f} ps ({times[-1]/1e6:.2f} ns)")
    print(f"   CV rango: {cv_centers.min():.3f} - {cv_centers.max():.3f} nm")
    print(f"   Altura inicial: {heights[0]:.3f} kJ/mol")
    print(f"   Altura final: {heights[-1]:.3f} kJ/mol (well-tempered decay)")
    
    return times, cv_centers, heights, sigmas


def calculate_pmf(cv_centers, heights, sigmas, cv_min=None, cv_max=None, bins=200):
    """
    Calcula PMF desde gaussianos de metadinámica.
    
    PMF(s) = -V_bias(s) = -Σ h_i * exp(-(s - s_i)²/(2σ_i²))
    
    Parámetros:
    -----------
    cv_centers : ndarray
        Posiciones de los gaussianos
    heights : ndarray
        Alturas de los gaussianos (kJ/mol)
    sigmas : ndarray
        Anchos de los gaussianos (nm)
    cv_min, cv_max : float
        Rango de CV para el grid (auto-detecta si None)
    bins : int
        Resolución del PMF
    
    Retorna:
    --------
    cv_grid : ndarray
        Valores de CV (nm)
    pmf : ndarray
        Free energy (kJ/mol)
    """
    
    # Auto-detectar rango si no se especifica
    if cv_min is None:
        cv_min = cv_centers.min() - 0.1
    if cv_max is None:
        cv_max = cv_centers.max() + 0.1
    
    print(f"\nCalculando PMF...")
    print(f"   Rango CV: [{cv_min:.2f}, {cv_max:.2f}] nm")
    print(f"   Resolución: {bins} bins (Δs = {(cv_max-cv_min)/bins:.4f} nm)")
    
    # Crear grid
    cv_grid = np.linspace(cv_min, cv_max, bins)
    bias = np.zeros(bins)
    
    # Sumar todos los gaussianos
    print(f"   Sumando {len(cv_centers)} gaussianos...")
    
    for i, cv in enumerate(cv_grid):
        # Mostrar progreso cada 10%
        if i % (bins // 10) == 0:
            print(f"      Progreso: {100*i//bins}%")
        
        # Contribución de cada gaussiano en este punto del grid
        for j in range(len(cv_centers)):
            bias[i] += heights[j] * np.exp(
                -(cv - cv_centers[j])**2 / (2 * sigmas[j]**2)
            )
    
    # PMF = -Bias (teorema de metadinámica)
    pmf = -bias
    
    # Normalizar: mínimo = 0 kJ/mol
    pmf_min = pmf.min()
    pmf -= pmf_min
    
    # Estadísticas
    barrier = pmf.max()
    barrier_position = cv_grid[pmf.argmax()]
    
    print(f"\n📊 Resultados:")
    print(f"   PMF mínimo: {pmf_min:.2f} kJ/mol (normalizado a 0)")
    print(f"   Barrera de activación: {barrier:.2f} kJ/mol")
    print(f"   Posición del TS: {barrier_position:.3f} nm")
    
    return cv_grid, pmf


def analyze_convergence(times, cv_centers, heights, sigmas, cv_min, cv_max, bins=200):
    """
    Analiza convergencia calculando PMF en diferentes ventanas temporales.
    
    Retorna:
    --------
    time_windows : list
        Tiempos de las ventanas (ns)
    pmfs : list
        PMFs correspondientes a cada ventana
    """
    
    print("\n🔄 Analizando convergencia...")
    
    # Ventanas: 25%, 50%, 75%, 100% del tiempo total
    fractions = [0.25, 0.50, 0.75, 1.00]
    time_windows = []
    pmfs = []
    cv_grid = np.linspace(cv_min, cv_max, bins)
    
    for frac in fractions:
        n = int(len(times) * frac)
        time_ns = times[n-1] / 1000  # ps → ns
        time_windows.append(time_ns)
        
        print(f"   Calculando PMF hasta {time_ns:.1f} ns ({int(frac*100)}%)...")
        
        # Calcular PMF con subset de gaussianos
        bias = np.zeros(bins)
        for i, cv in enumerate(cv_grid):
            for j in range(n):
                bias[i] += heights[j] * np.exp(
                    -(cv - cv_centers[j])**2 / (2 * sigmas[j]**2)
                )
        
        pmf = -bias
        pmf -= pmf.min()
        pmfs.append(pmf)
    
    return time_windows, pmfs, cv_grid


def plot_pmf(cv_grid, pmf, output_file='pmf_metad.png'):
    """
    Genera gráfico publication-quality del PMF.
    """
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot principal
    ax.plot(cv_grid, pmf, 'o-', linewidth=2.5, markersize=4, 
            color='#2E86AB', label='Metadinámica')
    
    # Marcar barrera
    barrier_idx = pmf.argmax()
    ax.plot(cv_grid[barrier_idx], pmf[barrier_idx], 
            '*', markersize=20, color='#A23B72', 
            label=f'TS: {pmf[barrier_idx]:.2f} kJ/mol')
    
    # Anotar barrera
    ax.annotate(
        f'ΔG = {pmf[barrier_idx]:.2f} kJ/mol',
        xy=(cv_grid[barrier_idx], pmf[barrier_idx]),
        xytext=(cv_grid[barrier_idx] + 0.3, pmf[barrier_idx] + 3),
        arrowprops=dict(arrowstyle='->', lw=2, color='#A23B72'),
        fontsize=12, fontweight='bold'
    )
    
    # Etiquetas
    ax.set_xlabel('Distancia Cα-Cα (nm)', fontsize=14, fontweight='bold')
    ax.set_ylabel('PMF (kJ/mol)', fontsize=14, fontweight='bold')
    ax.set_title('Potential of Mean Force - Metadinámica Well-Tempered', 
                 fontsize=15, fontweight='bold')
    
    ax.legend(fontsize=11, loc='upper right')
    ax.grid(alpha=0.3, linestyle='--')
    
    # Formato
    ax.tick_params(labelsize=11)
    ax.set_xlim(cv_grid.min(), cv_grid.max())
    ax.set_ylim(-2, pmf.max() + 5)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✅ PMF guardado: {output_file}")
    plt.close()


def plot_convergence(time_windows, pmfs, cv_grid, output_file='convergence_metad.png'):
    """
    Gráfico de convergencia temporal.
    """
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Panel izquierdo: PMFs superpuestos
    colors = ['#E8B4B8', '#C77A84', '#A94064', '#2E86AB']
    
    for i, (time_ns, pmf) in enumerate(zip(time_windows, pmfs)):
        ax1.plot(cv_grid, pmf, '-', linewidth=2, 
                color=colors[i], label=f'{time_ns:.0f} ns')
    
    ax1.set_xlabel('CV (nm)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('PMF (kJ/mol)', fontsize=12, fontweight='bold')
    ax1.set_title('Convergencia Temporal', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(alpha=0.3)
    
    # Panel derecho: Diferencia respecto a PMF final
    pmf_final = pmfs[-1]
    
    for i in range(len(pmfs) - 1):
        diff = np.abs(pmfs[i] - pmf_final)
        ax2.plot(cv_grid, diff, '-', linewidth=2, 
                color=colors[i], label=f'{time_windows[i]:.0f} ns vs final')
    
    ax2.axhline(2.0, color='red', linestyle='--', linewidth=2, 
                label='Umbral ±2 kJ/mol')
    
    ax2.set_xlabel('CV (nm)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('|ΔPMF| (kJ/mol)', fontsize=12, fontweight='bold')
    ax2.set_title('Diferencia Respecto a PMF Final', fontsize=13, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✅ Convergencia guardada: {output_file}")
    plt.close()


def save_pmf_data(cv_grid, pmf, output_file='pmf_metad.dat'):
    """
    Guarda PMF en formato texto para análisis posterior.
    """
    
    header = "# PMF extraído de metadinámica\n"
    header += "# CV(nm)  PMF(kJ/mol)\n"
    
    data = np.column_stack((cv_grid, pmf))
    
    np.savetxt(output_file, data, fmt='%.6f', header=header)
    print(f"✅ Datos guardados: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Extrae PMF desde archivo HILLS de metadinámica',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:
  python extract_pmf_from_hills.py metad_wnk1_HILLS.txt
  python extract_pmf_from_hills.py HILLS --cv-min 2.0 --cv-max 4.0 --bins 200
  python extract_pmf_from_hills.py HILLS --no-convergence (más rápido)
        """
    )
    
    parser.add_argument('hills_file', help='Archivo HILLS de metadinámica')
    parser.add_argument('--cv-min', type=float, default=None,
                       help='Mínimo de CV (auto-detecta si no se especifica)')
    parser.add_argument('--cv-max', type=float, default=None,
                       help='Máximo de CV (auto-detecta si no se especifica)')
    parser.add_argument('--bins', type=int, default=200,
                       help='Resolución del PMF (default: 200)')
    parser.add_argument('--output-prefix', default='pmf_metad',
                       help='Prefijo para archivos de salida (default: pmf_metad)')
    parser.add_argument('--no-convergence', action='store_true',
                       help='Omitir análisis de convergencia (más rápido)')
    
    args = parser.parse_args()
    
    print("="*60)
    print("EXTRACCIÓN DE PMF - METADINÁMICA")
    print("="*60)
    
    # 1. Cargar HILLS
    times, cv_centers, heights, sigmas = load_hills(args.hills_file)
    
    # 2. Calcular PMF
    cv_grid, pmf = calculate_pmf(
        cv_centers, heights, sigmas,
        cv_min=args.cv_min,
        cv_max=args.cv_max,
        bins=args.bins
    )
    
    # 3. Guardar datos
    save_pmf_data(cv_grid, pmf, f'{args.output_prefix}.dat')
    
    # 4. Plot principal
    plot_pmf(cv_grid, pmf, f'{args.output_prefix}.png')
    
    # 5. Análisis de convergencia (opcional)
    if not args.no_convergence:
        time_windows, pmfs, cv_grid_conv = analyze_convergence(
            times, cv_centers, heights, sigmas,
            cv_grid.min(), cv_grid.max(), args.bins
        )
        plot_convergence(time_windows, pmfs, cv_grid_conv, 
                        f'{args.output_prefix}_convergence.png')
    
    print("\n" + "="*60)
    print("✅ ANÁLISIS COMPLETADO")
    print("="*60)
    print("\nArchivos generados:")
    print(f"  1. {args.output_prefix}.dat              ← Datos PMF (texto)")
    print(f"  2. {args.output_prefix}.png              ← Gráfico PMF")
    if not args.no_convergence:
        print(f"  3. {args.output_prefix}_convergence.png ← Análisis convergencia")
    
    print("\n📌 Próximo paso:")
    print("   python compare_umbrella_metad.py pmf_final.dat pmf_metad.dat")
    print()


if __name__ == '__main__':
    main()
