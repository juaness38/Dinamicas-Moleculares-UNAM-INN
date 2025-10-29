#!/usr/bin/env python3
"""
Compara PMFs de umbrella sampling y metadinámica.
Validación cruzada para WNK1 C-terminal.
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys

def load_pmf(pmf_file):
    """
    Carga archivo PMF (formato: CV  PMF  [Error]).
    
    Retorna:
    --------
    cv : ndarray
        Valores de CV (nm)
    pmf : ndarray
        Free energy (kJ/mol)
    error : ndarray or None
        Errores (kJ/mol) si disponibles
    """
    
    print(f"Cargando {pmf_file}...")
    
    try:
        data = np.loadtxt(pmf_file, comments='#')
    except Exception as e:
        print(f"❌ Error al leer {pmf_file}: {e}")
        sys.exit(1)
    
    cv = data[:, 0]
    pmf = data[:, 1]
    error = data[:, 2] if data.shape[1] >= 3 else None
    
    print(f"✅ {len(cv)} puntos cargados")
    print(f"   CV rango: {cv.min():.3f} - {cv.max():.3f} nm")
    print(f"   PMF rango: {pmf.min():.3f} - {pmf.max():.3f} kJ/mol")
    print(f"   Barrera: {pmf.max():.2f} kJ/mol")
    
    return cv, pmf, error


def interpolate_pmf(cv_ref, pmf_ref, cv_target):
    """
    Interpola PMF al grid de referencia (para comparación).
    
    Parámetros:
    -----------
    cv_ref : ndarray
        CV de referencia (donde interpolar)
    pmf_ref : ndarray
        PMF a interpolar
    cv_target : ndarray
        CV objetivo (grid común)
    
    Retorna:
    --------
    pmf_interp : ndarray
        PMF interpolado al grid objetivo
    """
    
    return np.interp(cv_target, cv_ref, pmf_ref)


def calculate_rmsd(pmf1, pmf2):
    """
    Calcula RMSD entre dos perfiles de PMF.
    
    RMSD = sqrt(mean((PMF1 - PMF2)^2))
    """
    
    if len(pmf1) != len(pmf2):
        print("⚠️  Los PMFs tienen diferente longitud, imposible calcular RMSD")
        return None
    
    rmsd = np.sqrt(np.mean((pmf1 - pmf2)**2))
    return rmsd


def calculate_statistics(pmf1, pmf2):
    """
    Calcula estadísticas de comparación entre PMFs.
    
    Retorna:
    --------
    stats : dict
        Diccionario con métricas de comparación
    """
    
    stats = {}
    
    # RMSD global
    stats['rmsd'] = calculate_rmsd(pmf1, pmf2)
    
    # Diferencia en barreras
    barrier1 = pmf1.max()
    barrier2 = pmf2.max()
    stats['barrier1'] = barrier1
    stats['barrier2'] = barrier2
    stats['barrier_diff'] = abs(barrier1 - barrier2)
    
    # Posición de las barreras
    stats['barrier_pos1'] = pmf1.argmax()
    stats['barrier_pos2'] = pmf2.argmax()
    
    # Diferencia máxima
    stats['max_diff'] = np.abs(pmf1 - pmf2).max()
    
    # Diferencia promedio
    stats['mean_diff'] = np.abs(pmf1 - pmf2).mean()
    
    return stats


def plot_comparison(cv, pmf_umbrella, pmf_metad, error_umbrella, stats, 
                   output_file='comparison_umbrella_metad.png'):
    """
    Genera gráfico comparativo publication-quality.
    """
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), 
                                    gridspec_kw={'height_ratios': [2, 1]})
    
    # ========== Panel Superior: PMFs Superpuestos ==========
    
    # Umbrella sampling con barras de error
    ax1.plot(cv, pmf_umbrella, 'o-', linewidth=2.5, markersize=6,
            color='#2E86AB', label='Umbrella Sampling', alpha=0.9)
    
    if error_umbrella is not None:
        ax1.fill_between(cv, 
                         pmf_umbrella - error_umbrella, 
                         pmf_umbrella + error_umbrella,
                         color='#2E86AB', alpha=0.2)
    
    # Metadinámica
    ax1.plot(cv, pmf_metad, 's-', linewidth=2.5, markersize=6,
            color='#A23B72', label='Metadinámica', alpha=0.9)
    
    # Marcar barreras
    idx_umb = pmf_umbrella.argmax()
    idx_met = pmf_metad.argmax()
    
    ax1.plot(cv[idx_umb], pmf_umbrella[idx_umb], '*', 
            markersize=20, color='#2E86AB', zorder=5)
    ax1.plot(cv[idx_met], pmf_metad[idx_met], '*', 
            markersize=20, color='#A23B72', zorder=5)
    
    # Anotaciones de barreras
    ax1.annotate(
        f'Umbrella: {stats["barrier1"]:.2f} kJ/mol',
        xy=(cv[idx_umb], pmf_umbrella[idx_umb]),
        xytext=(cv[idx_umb] - 0.4, pmf_umbrella[idx_umb] + 3),
        arrowprops=dict(arrowstyle='->', lw=2, color='#2E86AB'),
        fontsize=11, fontweight='bold', color='#2E86AB'
    )
    
    ax1.annotate(
        f'Metadinámica: {stats["barrier2"]:.2f} kJ/mol',
        xy=(cv[idx_met], pmf_metad[idx_met]),
        xytext=(cv[idx_met] + 0.2, pmf_metad[idx_met] + 3),
        arrowprops=dict(arrowstyle='->', lw=2, color='#A23B72'),
        fontsize=11, fontweight='bold', color='#A23B72'
    )
    
    # Etiquetas
    ax1.set_ylabel('PMF (kJ/mol)', fontsize=14, fontweight='bold')
    ax1.set_title('Validación Cruzada: Umbrella Sampling vs. Metadinámica', 
                 fontsize=15, fontweight='bold')
    ax1.legend(fontsize=12, loc='upper right')
    ax1.grid(alpha=0.3, linestyle='--')
    ax1.tick_params(labelsize=11)
    
    # ========== Panel Inferior: Diferencia ==========
    
    diff = pmf_metad - pmf_umbrella
    
    ax2.plot(cv, diff, 'o-', linewidth=2, markersize=5, 
            color='#F18F01', label='Metad - Umbrella')
    ax2.axhline(0, color='black', linestyle='-', linewidth=1.5, alpha=0.5)
    ax2.axhline(2, color='red', linestyle='--', linewidth=2, alpha=0.7, 
               label='Umbral ±2 kJ/mol')
    ax2.axhline(-2, color='red', linestyle='--', linewidth=2, alpha=0.7)
    
    # Sombrear región de aceptación
    ax2.fill_between(cv, -2, 2, color='green', alpha=0.1)
    
    # Etiquetas
    ax2.set_xlabel('Distancia Cα-Cα (nm)', fontsize=14, fontweight='bold')
    ax2.set_ylabel('ΔPMF (kJ/mol)', fontsize=14, fontweight='bold')
    ax2.set_title('Diferencia entre Métodos', fontsize=13, fontweight='bold')
    ax2.legend(fontsize=11, loc='upper right')
    ax2.grid(alpha=0.3, linestyle='--')
    ax2.tick_params(labelsize=11)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✅ Comparación guardada: {output_file}")
    plt.close()


def print_validation_report(stats, rmsd_threshold=2.0):
    """
    Imprime reporte de validación con criterios de aceptación.
    """
    
    print("\n" + "="*60)
    print("📊 REPORTE DE VALIDACIÓN CRUZADA")
    print("="*60)
    
    print("\n1️⃣  BARRERAS DE ACTIVACIÓN:")
    print(f"   Umbrella Sampling:  {stats['barrier1']:.2f} kJ/mol")
    print(f"   Metadinámica:       {stats['barrier2']:.2f} kJ/mol")
    print(f"   Diferencia:         {stats['barrier_diff']:.2f} kJ/mol")
    
    # Criterio 1: Diferencia en barrera < 2 kJ/mol
    criterion1 = stats['barrier_diff'] < 2.0
    status1 = "✅ PASS" if criterion1 else "❌ FAIL"
    print(f"   Criterio (<2 kJ/mol): {status1}")
    
    print("\n2️⃣  RMSD GLOBAL:")
    print(f"   RMSD:               {stats['rmsd']:.3f} kJ/mol")
    print(f"   Umbral aceptable:   <{rmsd_threshold:.1f} kJ/mol")
    
    # Criterio 2: RMSD < 2 kJ/mol
    criterion2 = stats['rmsd'] < rmsd_threshold
    status2 = "✅ PASS" if criterion2 else "❌ FAIL"
    print(f"   Criterio:           {status2}")
    
    print("\n3️⃣  DIFERENCIAS LOCALES:")
    print(f"   Diferencia máxima:  {stats['max_diff']:.2f} kJ/mol")
    print(f"   Diferencia promedio: {stats['mean_diff']:.2f} kJ/mol")
    
    # Criterio 3: Diferencia máxima < 3 kJ/mol
    criterion3 = stats['max_diff'] < 3.0
    status3 = "✅ PASS" if criterion3 else "⚠️  WARNING"
    print(f"   Criterio (<3 kJ/mol): {status3}")
    
    # Veredicto final
    print("\n" + "="*60)
    if criterion1 and criterion2:
        print("✅ VALIDACIÓN EXITOSA")
        print("="*60)
        print("\n🎯 Ambos métodos convergen al mismo PMF.")
        print("   Los resultados son estadísticamente consistentes.")
        print("   Puedes confiar en la barrera de activación calculada.")
    elif criterion1 or criterion2:
        print("⚠️  VALIDACIÓN PARCIAL")
        print("="*60)
        print("\n⚠️  Los métodos muestran acuerdo razonable, pero hay diferencias.")
        print("   Considera:")
        print("   • Extender tiempo de simulación de metadinámica")
        print("   • Verificar convergencia de ambos métodos")
        print("   • Analizar histogramas de solapamiento (umbrella)")
    else:
        print("❌ VALIDACIÓN FALLIDA")
        print("="*60)
        print("\n❌ Los métodos NO convergen al mismo PMF.")
        print("   Problemas posibles:")
        print("   • Metadinámica no ha convergido (necesita más tiempo)")
        print("   • Umbrella: ventanas mal espaciadas o poco solapamiento")
        print("   • Error en definición de CV o parámetros de bias")
        print("   • Sistema tiene múltiples pathways (necesita análisis 2D)")
    
    print()


def save_comparison_data(cv, pmf_umbrella, pmf_metad, output_file='comparison.dat'):
    """
    Guarda datos de comparación en archivo de texto.
    """
    
    header = "# Comparación Umbrella vs Metadinámica\n"
    header += "# CV(nm)  PMF_Umbrella(kJ/mol)  PMF_Metad(kJ/mol)  Diff(kJ/mol)\n"
    
    diff = pmf_metad - pmf_umbrella
    data = np.column_stack((cv, pmf_umbrella, pmf_metad, diff))
    
    np.savetxt(output_file, data, fmt='%.6f', header=header)
    print(f"✅ Datos de comparación guardados: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Compara PMFs de umbrella sampling y metadinámica',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:
  python compare_umbrella_metad.py pmf_final.dat pmf_metad.dat
  python compare_umbrella_metad.py umbrella.dat metad.dat --output comparison
  python compare_umbrella_metad.py pmf1.dat pmf2.dat --rmsd-threshold 3.0
        """
    )
    
    parser.add_argument('umbrella_file', 
                       help='Archivo PMF de umbrella sampling')
    parser.add_argument('metad_file', 
                       help='Archivo PMF de metadinámica')
    parser.add_argument('--output', default='comparison',
                       help='Prefijo para archivos de salida (default: comparison)')
    parser.add_argument('--rmsd-threshold', type=float, default=2.0,
                       help='Umbral de RMSD aceptable (default: 2.0 kJ/mol)')
    
    args = parser.parse_args()
    
    print("="*60)
    print("VALIDACIÓN CRUZADA: UMBRELLA vs METADINÁMICA")
    print("="*60)
    
    # 1. Cargar PMFs
    print("\n📂 Cargando datos...")
    cv_umb, pmf_umb, error_umb = load_pmf(args.umbrella_file)
    cv_met, pmf_met, error_met = load_pmf(args.metad_file)
    
    # 2. Interpolar al mismo grid (usar el de umbrella como referencia)
    print("\n🔄 Interpolando metadinámica al grid de umbrella...")
    
    if not np.array_equal(cv_umb, cv_met):
        print("   ⚠️  Grids diferentes, interpolando...")
        pmf_met_interp = interpolate_pmf(cv_met, pmf_met, cv_umb)
    else:
        print("   ✅ Grids idénticos, no requiere interpolación")
        pmf_met_interp = pmf_met
    
    cv_common = cv_umb
    
    # 3. Calcular estadísticas
    print("\n📊 Calculando estadísticas de comparación...")
    stats = calculate_statistics(pmf_umb, pmf_met_interp)
    
    # 4. Generar gráfico
    print("\n📈 Generando gráfico comparativo...")
    plot_comparison(cv_common, pmf_umb, pmf_met_interp, error_umb, stats,
                   f'{args.output}.png')
    
    # 5. Guardar datos
    save_comparison_data(cv_common, pmf_umb, pmf_met_interp, 
                        f'{args.output}.dat')
    
    # 6. Reporte de validación
    print_validation_report(stats, args.rmsd_threshold)
    
    print("\n📁 Archivos generados:")
    print(f"   1. {args.output}.png  ← Gráfico comparativo")
    print(f"   2. {args.output}.dat  ← Datos numéricos")
    print()


if __name__ == '__main__':
    main()
