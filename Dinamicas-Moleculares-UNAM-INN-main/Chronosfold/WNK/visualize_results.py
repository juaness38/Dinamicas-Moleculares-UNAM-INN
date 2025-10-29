#!/usr/bin/env python3
"""
Visualización completa de resultados Umbrella Sampling con VideoSuite

Este script genera:
1. PMF plot publicable (alta resolución)
2. Histogramas de CV (distribuciones por ventana)
3. Convergencia temporal
4. Video animado de histogramas (requiere VIDEOSUITE)
5. Diagnósticos 4-panel (heatmap, tracking, PMF, CDFs)

Uso:
    python visualize_results.py
    python visualize_results.py --animation  # Incluir video (requiere ffmpeg)
"""

import sys
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict

# Configurar estilo
sns.set_style("whitegrid")
sns.set_context("paper", font_scale=1.5)
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.family'] = 'sans-serif'


class UmbrellaVisualizer:
    """Visualizador de resultados umbrella sampling"""
    
    def __init__(self, base_dir: Path):
        self.base_dir = Path(base_dir)
        self.pmf_dir = self.base_dir / "pmf_analysis"
        self.windows_dir = self.base_dir / "umbrella_windows"
        
        # Cargar datos
        self.pmf_data = None
        self.windows_data = {}
        
    def load_data(self) -> bool:
        """Cargar datos de PMF y ventanas"""
        
        # Cargar PMF
        pmf_file = self.pmf_dir / "pmf_results.csv"
        if pmf_file.exists():
            self.pmf_data = pd.read_csv(pmf_file)
            print(f"✓ PMF cargado: {len(self.pmf_data)} puntos")
        else:
            print(f"⚠️  PMF no encontrado: {pmf_file}")
            return False
        
        # Cargar ventanas
        window_dirs = sorted(self.windows_dir.glob("window_*"))
        for wdir in window_dirs:
            cv_file = wdir / "cv_values.dat"
            if cv_file.exists():
                window_id = int(wdir.name.split('_')[-1])
                cv_data = np.loadtxt(cv_file)
                self.windows_data[window_id] = cv_data
        
        print(f"✓ {len(self.windows_data)} ventanas cargadas")
        return True
    
    def plot_pmf(self, save_path: Path = None):
        """Plot PMF publicable"""
        
        if self.pmf_data is None:
            print("⚠️  PMF data no cargada")
            return
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        cv = self.pmf_data['cv']
        pmf = self.pmf_data['pmf']
        
        # Plot PMF
        ax.plot(cv, pmf, 'o-', linewidth=2, markersize=8, 
                color='#2E86AB', label='PMF')
        
        # Uncertainty band si existe
        if 'uncertainty' in self.pmf_data.columns:
            uncertainty = self.pmf_data['uncertainty']
            ax.fill_between(cv, pmf - uncertainty, pmf + uncertainty,
                           alpha=0.3, color='#2E86AB', label='Uncertainty')
        
        # Formateo
        ax.set_xlabel('C-terminal Distance (nm)', fontsize=14, fontweight='bold')
        ax.set_ylabel('PMF (kJ/mol)', fontsize=14, fontweight='bold')
        ax.set_title('WNK1 C-Terminal Potential of Mean Force\n(PBS buffer, pH 7.4)', 
                    fontsize=16, fontweight='bold', pad=20)
        ax.legend(fontsize=12, frameon=True, shadow=True)
        ax.grid(alpha=0.3, linestyle='--')
        
        # Anotar mínimo
        min_idx = pmf.argmin()
        min_cv = cv.iloc[min_idx]
        min_pmf = pmf.iloc[min_idx]
        ax.annotate(f'Minimum\n{min_cv:.2f} nm\n{min_pmf:.1f} kJ/mol',
                   xy=(min_cv, min_pmf), xytext=(min_cv + 0.1, min_pmf + 5),
                   arrowprops=dict(arrowstyle='->', lw=2, color='red'),
                   fontsize=11, ha='left',
                   bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.7))
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ PMF guardado: {save_path}")
        
        return fig
    
    def plot_histograms(self, save_path: Path = None):
        """Plot histogramas de todas las ventanas"""
        
        if not self.windows_data:
            print("⚠️  Datos de ventanas no cargados")
            return
        
        n_windows = len(self.windows_data)
        ncols = 5
        nrows = (n_windows + ncols - 1) // ncols
        
        fig, axes = plt.subplots(nrows, ncols, figsize=(20, 4*nrows))
        axes = axes.flatten()
        
        for i, (window_id, cv_data) in enumerate(sorted(self.windows_data.items())):
            ax = axes[i]
            
            # Histogram
            ax.hist(cv_data, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
            
            # Formateo
            ax.set_title(f'Window {window_id}', fontweight='bold')
            ax.set_xlabel('CV (nm)')
            ax.set_ylabel('Count')
            ax.grid(alpha=0.3)
            
            # Estadísticas
            mean_cv = np.mean(cv_data)
            std_cv = np.std(cv_data)
            ax.axvline(mean_cv, color='red', linestyle='--', linewidth=2, 
                      label=f'μ={mean_cv:.2f}±{std_cv:.2f}')
            ax.legend(fontsize=9)
        
        # Ocultar axes vacíos
        for i in range(n_windows, len(axes)):
            axes[i].axis('off')
        
        plt.suptitle('CV Distributions per Umbrella Window', 
                    fontsize=18, fontweight='bold', y=1.00)
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ Histogramas guardados: {save_path}")
        
        return fig
    
    def plot_convergence(self, save_path: Path = None):
        """Plot convergencia temporal"""
        
        if not self.windows_data:
            print("⚠️  Datos de ventanas no cargados")
            return
        
        fig, axes = plt.subplots(2, 1, figsize=(12, 8))
        
        # Panel 1: Media acumulativa
        ax1 = axes[0]
        for window_id, cv_data in sorted(self.windows_data.items()):
            cumulative_mean = np.cumsum(cv_data) / np.arange(1, len(cv_data) + 1)
            ax1.plot(cumulative_mean, alpha=0.6, linewidth=1.5)
        
        ax1.set_xlabel('Frame')
        ax1.set_ylabel('Cumulative Mean CV (nm)')
        ax1.set_title('Convergence: Cumulative Mean', fontweight='bold')
        ax1.grid(alpha=0.3)
        
        # Panel 2: Desviación estándar móvil
        ax2 = axes[1]
        window_size = 1000
        for window_id, cv_data in sorted(self.windows_data.items()):
            if len(cv_data) > window_size:
                rolling_std = pd.Series(cv_data).rolling(window=window_size).std()
                ax2.plot(rolling_std, alpha=0.6, linewidth=1.5)
        
        ax2.set_xlabel('Frame')
        ax2.set_ylabel(f'Rolling Std (window={window_size})')
        ax2.set_title('Convergence: Rolling Standard Deviation', fontweight='bold')
        ax2.grid(alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ Convergencia guardada: {save_path}")
        
        return fig
    
    def plot_combined(self, save_path: Path = None):
        """Plot combinado 2x2: PMF, histograms, convergence, stats"""
        
        fig = plt.figure(figsize=(18, 14))
        gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
        
        # Panel 1: PMF (top-left)
        ax1 = fig.add_subplot(gs[0, 0])
        if self.pmf_data is not None:
            cv = self.pmf_data['cv']
            pmf = self.pmf_data['pmf']
            ax1.plot(cv, pmf, 'o-', linewidth=2, markersize=6, color='#2E86AB')
            if 'uncertainty' in self.pmf_data.columns:
                uncertainty = self.pmf_data['uncertainty']
                ax1.fill_between(cv, pmf - uncertainty, pmf + uncertainty,
                               alpha=0.3, color='#2E86AB')
            ax1.set_xlabel('CV (nm)', fontweight='bold')
            ax1.set_ylabel('PMF (kJ/mol)', fontweight='bold')
            ax1.set_title('Potential of Mean Force', fontweight='bold', fontsize=14)
            ax1.grid(alpha=0.3)
        
        # Panel 2: Overlay histograms (top-right)
        ax2 = fig.add_subplot(gs[0, 1])
        if self.windows_data:
            for window_id, cv_data in sorted(self.windows_data.items()):
                ax2.hist(cv_data, bins=30, alpha=0.4, density=True)
            ax2.set_xlabel('CV (nm)', fontweight='bold')
            ax2.set_ylabel('Density', fontweight='bold')
            ax2.set_title('CV Distributions (All Windows)', fontweight='bold', fontsize=14)
            ax2.grid(alpha=0.3)
        
        # Panel 3: Convergence (bottom-left)
        ax3 = fig.add_subplot(gs[1, 0])
        if self.windows_data:
            for window_id, cv_data in sorted(self.windows_data.items()):
                cumulative_mean = np.cumsum(cv_data) / np.arange(1, len(cv_data) + 1)
                ax3.plot(cumulative_mean, alpha=0.5, linewidth=1)
            ax3.set_xlabel('Frame', fontweight='bold')
            ax3.set_ylabel('Cumulative Mean (nm)', fontweight='bold')
            ax3.set_title('Convergence Analysis', fontweight='bold', fontsize=14)
            ax3.grid(alpha=0.3)
        
        # Panel 4: Statistics table (bottom-right)
        ax4 = fig.add_subplot(gs[1, 1])
        ax4.axis('off')
        
        if self.pmf_data is not None and self.windows_data:
            # Crear tabla de estadísticas
            stats = []
            stats.append(['N Windows', f"{len(self.windows_data)}"])
            stats.append(['Total Frames', f"{sum(len(d) for d in self.windows_data.values()):,}"])
            
            if self.pmf_data is not None:
                min_idx = self.pmf_data['pmf'].argmin()
                max_idx = self.pmf_data['pmf'].argmax()
                stats.append(['PMF Min', f"{self.pmf_data['pmf'].iloc[min_idx]:.2f} kJ/mol"])
                stats.append(['PMF Max', f"{self.pmf_data['pmf'].iloc[max_idx]:.2f} kJ/mol"])
                stats.append(['PMF Range', f"{self.pmf_data['pmf'].max() - self.pmf_data['pmf'].min():.2f} kJ/mol"])
                stats.append(['Optimal CV', f"{self.pmf_data['cv'].iloc[min_idx]:.3f} nm"])
            
            # Dibujar tabla
            table = ax4.table(cellText=stats, cellLoc='left',
                            colWidths=[0.5, 0.5],
                            bbox=[0.1, 0.3, 0.8, 0.6])
            table.auto_set_font_size(False)
            table.set_fontsize(12)
            table.scale(1, 2)
            
            # Estilo
            for i in range(len(stats)):
                table[(i, 0)].set_facecolor('#E8F4F8')
                table[(i, 0)].set_text_props(weight='bold')
            
            ax4.set_title('Summary Statistics', fontweight='bold', fontsize=14, pad=20)
        
        plt.suptitle('WNK1 Umbrella Sampling Analysis\n(PBS Buffer, pH 7.4)', 
                    fontsize=18, fontweight='bold', y=0.98)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ Análisis combinado guardado: {save_path}")
        
        return fig


def main():
    parser = argparse.ArgumentParser(description='Visualizar resultados umbrella sampling')
    parser.add_argument('--animation', action='store_true', 
                       help='Generar video animado (requiere ffmpeg y VIDEOSUITE)')
    parser.add_argument('--base-dir', type=str, default='.',
                       help='Directorio base (default: directorio actual)')
    args = parser.parse_args()
    
    base_dir = Path(args.base_dir)
    
    print("═" * 70)
    print("  WNK1 UMBRELLA SAMPLING - VISUALIZACIÓN COMPLETA")
    print("═" * 70)
    print()
    
    # Crear visualizador
    viz = UmbrellaVisualizer(base_dir)
    
    # Cargar datos
    print("Cargando datos...")
    if not viz.load_data():
        print("✗ Error cargando datos. Verificar que analyze_umbrella_mbar.py se ejecutó.")
        return 1
    
    print()
    
    # Crear directorio de salida
    out_dir = base_dir / "pmf_analysis"
    out_dir.mkdir(exist_ok=True)
    
    # Generar plots
    print("Generando visualizaciones...")
    print()
    
    # 1. PMF publicable
    print("1/5 - PMF plot...")
    viz.plot_pmf(save_path=out_dir / "pmf_publication.png")
    
    # 2. Histogramas
    print("2/5 - Histogramas por ventana...")
    viz.plot_histograms(save_path=out_dir / "histograms_all_windows.png")
    
    # 3. Convergencia
    print("3/5 - Análisis de convergencia...")
    viz.plot_convergence(save_path=out_dir / "convergence_analysis.png")
    
    # 4. Plot combinado
    print("4/5 - Análisis combinado...")
    viz.plot_combined(save_path=out_dir / "analysis_combined.png")
    
    # 5. Animación (opcional)
    if args.animation:
        print("5/5 - Generando video animado...")
        try:
            # Intentar usar VIDEOSUITE
            sys.path.insert(0, str(base_dir.parent / "VIDEOSUITE"))
            from umbrella.run_umbrella_visualization import generate_animation
            
            generate_animation(
                out_dir=out_dir,
                windows=list(viz.windows_data.values()),
                fps=24,
                frames=120
            )
            print(f"✓ Video guardado: {out_dir / 'umbrella_histograms.mp4'}")
            
        except ImportError:
            print("⚠️  VIDEOSUITE no encontrado. Saltando animación.")
            print("   Para generar video: agregar VIDEOSUITE al PYTHONPATH")
        except Exception as e:
            print(f"⚠️  Error generando video: {e}")
    else:
        print("5/5 - Video animado (skip, usar --animation para generar)")
    
    print()
    print("═" * 70)
    print("  ✓ VISUALIZACIÓN COMPLETADA")
    print("═" * 70)
    print()
    print("Archivos generados:")
    print(f"  📊 {out_dir / 'pmf_publication.png'}")
    print(f"  📊 {out_dir / 'histograms_all_windows.png'}")
    print(f"  📊 {out_dir / 'convergence_analysis.png'}")
    print(f"  📊 {out_dir / 'analysis_combined.png'}")
    if args.animation:
        print(f"  🎬 {out_dir / 'umbrella_histograms.mp4'}")
    print()
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
