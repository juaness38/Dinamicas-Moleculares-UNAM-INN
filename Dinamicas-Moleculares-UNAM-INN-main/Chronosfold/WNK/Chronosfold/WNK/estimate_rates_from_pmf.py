#!/usr/bin/env python3
"""
Estima tasas cinéticas desde PMF usando Teoría del Estado de Transición (Eyring)

Este script calcula tasas de transición aproximadas desde barreras de energía
libre obtenidas de umbrella sampling. Usa la ecuación de Eyring:

    k = (k_B * T / h) * exp(-ΔG‡ / RT)

IMPORTANTE - Limitaciones:
1. Asume no hay recrossing (transmission coefficient κ = 1)
2. No considera fricción conformacional (puede sobre-estimar k)
3. Es una estimación de ORDEN DE MAGNITUD, no cinética exacta
4. Para cinética precisa: WE, Milestoning, TPS, o MD ultra-larga

Uso apropiado:
- Comparar barreras entre mutantes (ΔΔG → fold-change en k)
- Estimar tiempo característico de transición (orden de magnitud)
- Identificar si proceso es rápido/lento a escala celular
- Diseño racional de mutaciones

Autor: Pipeline de Umbrella Sampling
Fecha: 2024
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import constants
from pathlib import Path
import sys

# Constantes físicas
k_B = constants.Boltzmann  # 1.380649e-23 J/K
h = constants.Planck       # 6.62607015e-34 J·s
R = constants.R            # 8.314462618 J/(mol·K)
N_A = constants.Avogadro   # 6.02214076e23 /mol

class EyringRateCalculator:
    """
    Calcula tasas cinéticas desde PMF usando ecuación de Eyring
    """
    
    def __init__(self, temperature=310.0):
        """
        Args:
            temperature: Temperatura en Kelvin (default: 310 K = 37°C)
        """
        self.temperature = temperature
        self.kT = k_B * temperature  # Energía térmica (J)
        self.RT = R * temperature    # Energía térmica (J/mol)
        
        print("="*70)
        print("EYRING RATE CALCULATOR")
        print("="*70)
        print(f"Temperatura: {temperature} K ({temperature-273.15:.1f}°C)")
        print(f"kT = {self.kT * 1e21:.3f} zJ (zeptojoules)")
        print(f"RT = {self.RT / 1000:.3f} kJ/mol")
        print("="*70)
    
    def prefactor(self):
        """
        Calcula factor prefactor de Eyring: k_B*T / h
        
        Returns:
            Prefactor en s^-1
        """
        pf = (k_B * self.temperature) / h
        return pf
    
    def boltzmann_factor(self, barrier_kj_mol):
        """
        Calcula factor de Boltzmann: exp(-ΔG / RT)
        
        Args:
            barrier_kj_mol: Barrera en kJ/mol
            
        Returns:
            Factor exponencial (adimensional)
        """
        barrier_j_mol = barrier_kj_mol * 1000  # kJ/mol → J/mol
        factor = np.exp(-barrier_j_mol / self.RT)
        return factor
    
    def calculate_rate(self, barrier_kj_mol, 
                      transmission_coeff=1.0,
                      friction_factor=1.0):
        """
        Calcula tasa de transición usando Eyring
        
        Args:
            barrier_kj_mol: Barrera energética (kJ/mol)
            transmission_coeff: Coeficiente de transmisión κ (default: 1.0)
            friction_factor: Factor de corrección por fricción (default: 1.0)
                            Para proteínas en agua: ~100-1000
            
        Returns:
            Tasa k en s^-1
        """
        pf = self.prefactor()
        bf = self.boltzmann_factor(barrier_kj_mol)
        
        k = transmission_coeff * (pf * bf) / friction_factor
        
        return k
    
    def half_life(self, rate_s):
        """
        Calcula tiempo de vida media desde tasa
        
        Args:
            rate_s: Tasa en s^-1
            
        Returns:
            t_1/2 en segundos
        """
        return np.log(2) / rate_s
    
    def mean_time(self, rate_s):
        """
        Calcula tiempo medio de transición
        
        Args:
            rate_s: Tasa en s^-1
            
        Returns:
            <t> en segundos
        """
        return 1.0 / rate_s
    
    def format_time(self, time_seconds):
        """
        Formatea tiempo en unidades apropiadas
        
        Args:
            time_seconds: Tiempo en segundos
            
        Returns:
            String formateado
        """
        if time_seconds < 1e-9:
            return f"{time_seconds * 1e12:.2f} ps"
        elif time_seconds < 1e-6:
            return f"{time_seconds * 1e9:.2f} ns"
        elif time_seconds < 1e-3:
            return f"{time_seconds * 1e6:.2f} μs"
        elif time_seconds < 1:
            return f"{time_seconds * 1e3:.2f} ms"
        elif time_seconds < 60:
            return f"{time_seconds:.2f} s"
        elif time_seconds < 3600:
            return f"{time_seconds / 60:.2f} min"
        elif time_seconds < 86400:
            return f"{time_seconds / 3600:.2f} h"
        else:
            return f"{time_seconds / 86400:.2f} days"
    
    def analyze_pmf(self, pmf_file):
        """
        Analiza PMF completo y calcula tasas para todas las barreras
        
        Args:
            pmf_file: Archivo con PMF (formato: CV PMF)
            
        Returns:
            Dict con análisis completo
        """
        print(f"\n{'='*70}")
        print(f"ANÁLISIS DE PMF: {pmf_file}")
        print(f"{'='*70}")
        
        # Cargar PMF
        data = np.loadtxt(pmf_file)
        
        if data.ndim == 1:
            cv = np.linspace(2.0, 4.0, len(data))
            pmf = data
        else:
            cv = data[:, 0]
            pmf = data[:, 1]
        
        pmf -= pmf.min()  # Referencia al mínimo global
        
        print(f"Datos cargados: {len(cv)} puntos")
        print(f"  Rango CV: {cv.min():.2f} - {cv.max():.2f} nm")
        print(f"  PMF: {pmf.min():.2f} - {pmf.max():.2f} kJ/mol")
        
        # Identificar características
        min_idx = np.argmin(pmf)
        max_idx = np.argmax(pmf)
        
        cv_min = cv[min_idx]
        cv_max = cv[max_idx]
        pmf_min = pmf[min_idx]
        pmf_max = pmf[max_idx]
        
        barrier = pmf_max - pmf_min
        
        print(f"\nCaracterísticas del PMF:")
        print(f"  Mínimo global: CV = {cv_min:.2f} nm, PMF = {pmf_min:.2f} kJ/mol")
        print(f"  Máximo (TS): CV = {cv_max:.2f} nm, PMF = {pmf_max:.2f} kJ/mol")
        print(f"  Barrera: ΔG‡ = {barrier:.2f} kJ/mol")
        
        # Buscar mínimo secundario (si existe)
        # Región después del máximo
        if max_idx < len(pmf) - 10:
            post_max = pmf[max_idx+5:]
            if len(post_max) > 0:
                sec_min_idx = max_idx + 5 + np.argmin(post_max)
                cv_sec = cv[sec_min_idx]
                pmf_sec = pmf[sec_min_idx]
                
                if pmf_sec < pmf_max - 5:  # Al menos 5 kJ/mol más bajo
                    print(f"  Mínimo secundario: CV = {cv_sec:.2f} nm, "
                          f"PMF = {pmf_sec:.2f} kJ/mol")
                    barrier_reverse = pmf_max - pmf_sec
                    print(f"  Barrera reversa: ΔG‡_rev = {barrier_reverse:.2f} kJ/mol")
        
        # Calcular tasas
        print(f"\n{'='*70}")
        print("CÁLCULO DE TASAS CINÉTICAS")
        print(f"{'='*70}")
        
        # Escenario 1: Sin fricción (límite superior)
        print("\n1. Sin fricción (límite superior teórico):")
        k_no_friction = self.calculate_rate(barrier, friction_factor=1.0)
        t_mean_nf = self.mean_time(k_no_friction)
        
        print(f"   k = {k_no_friction:.2e} s⁻¹")
        print(f"   <t> = {self.format_time(t_mean_nf)}")
        
        # Escenario 2: Fricción baja (proteína pequeña, poco solvente)
        print("\n2. Fricción baja (proteína pequeña, FF = 100):")
        k_low_friction = self.calculate_rate(barrier, friction_factor=100)
        t_mean_lf = self.mean_time(k_low_friction)
        
        print(f"   k = {k_low_friction:.2e} s⁻¹")
        print(f"   <t> = {self.format_time(t_mean_lf)}")
        
        # Escenario 3: Fricción media (proteína típica)
        print("\n3. Fricción media (proteína típica, FF = 500):")
        k_med_friction = self.calculate_rate(barrier, friction_factor=500)
        t_mean_mf = self.mean_time(k_med_friction)
        
        print(f"   k = {k_med_friction:.2e} s⁻¹")
        print(f"   <t> = {self.format_time(t_mean_mf)}")
        
        # Escenario 4: Fricción alta (proteína grande, mucho solvente)
        print("\n4. Fricción alta (proteína grande, FF = 1000):")
        k_high_friction = self.calculate_rate(barrier, friction_factor=1000)
        t_mean_hf = self.mean_time(k_high_friction)
        
        print(f"   k = {k_high_friction:.2e} s⁻¹")
        print(f"   <t> = {self.format_time(t_mean_hf)}")
        
        # Interpretación biológica
        print(f"\n{'='*70}")
        print("INTERPRETACIÓN BIOLÓGICA")
        print(f"{'='*70}")
        
        # Comparar con escalas de tiempo celulares
        timescales = {
            'Vibración atómica': 1e-15,
            'Movimiento backbone local': 1e-12,
            'Fluctuación loop': 1e-9,
            'Cambio conformacional dominio': 1e-6,
            'Plegamiento proteína pequeña': 1e-3,
            'Señalización celular': 1.0,
            'Ciclo celular': 3600,
        }
        
        print(f"\nEscala de tiempo (fricción media: {self.format_time(t_mean_mf)}):")
        
        for process, time_sec in timescales.items():
            if t_mean_mf < time_sec:
                print(f"  ✓ Más rápido que: {process} ({self.format_time(time_sec)})")
            else:
                print(f"  ✗ Más lento que: {process} ({self.format_time(time_sec)})")
        
        # Sensibilidad a mutaciones
        print(f"\n{'='*70}")
        print("SENSIBILIDAD A MUTACIONES")
        print(f"{'='*70}")
        print("\nCambios hipotéticos en barrera:")
        
        for delta_g in [-10, -5, -2, 2, 5, 10]:
            new_barrier = barrier + delta_g
            new_k = self.calculate_rate(new_barrier, friction_factor=500)
            fold_change = new_k / k_med_friction
            
            if delta_g < 0:
                effect = f"estabiliza TS ({abs(delta_g)} kJ/mol más bajo)"
            else:
                effect = f"desestabiliza TS ({delta_g} kJ/mol más alto)"
            
            print(f"  ΔΔG = {delta_g:+3d} kJ/mol ({effect})")
            print(f"    → k_mutante / k_WT = {fold_change:.1f}×")
            print(f"    → <t>_mutante = {self.format_time(self.mean_time(new_k))}")
        
        # Visualización
        self.plot_analysis(cv, pmf, barrier, 
                          k_med_friction, t_mean_mf,
                          output='rate_analysis.png')
        
        return {
            'cv': cv,
            'pmf': pmf,
            'barrier': barrier,
            'rates': {
                'no_friction': k_no_friction,
                'low_friction': k_low_friction,
                'med_friction': k_med_friction,
                'high_friction': k_high_friction,
            },
            'mean_times': {
                'no_friction': t_mean_nf,
                'low_friction': t_mean_lf,
                'med_friction': t_mean_mf,
                'high_friction': t_mean_hf,
            }
        }
    
    def plot_analysis(self, cv, pmf, barrier, rate, mean_time, output='rate_analysis.png'):
        """
        Genera figura de análisis
        """
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Panel 1: PMF con barrera
        ax1 = axes[0, 0]
        ax1.plot(cv, pmf, 'k-', linewidth=2.5, label='PMF')
        ax1.fill_between(cv, 0, pmf, alpha=0.2, color='gray')
        
        min_idx = np.argmin(pmf)
        max_idx = np.argmax(pmf)
        
        ax1.plot(cv[min_idx], pmf[min_idx], 'go', markersize=15, 
                label='Mínimo', zorder=5)
        ax1.plot(cv[max_idx], pmf[max_idx], 'r^', markersize=15,
                label='Estado Transición', zorder=5)
        
        # Flecha de barrera
        ax1.annotate('', xy=(cv[max_idx], pmf[max_idx]),
                    xytext=(cv[max_idx], pmf[min_idx]),
                    arrowprops=dict(arrowstyle='<->', color='red', lw=3))
        ax1.text(cv[max_idx] + 0.1, (pmf[max_idx] + pmf[min_idx])/2,
                f'ΔG‡ = {barrier:.1f} kJ/mol',
                fontsize=12, fontweight='bold', color='red')
        
        ax1.set_xlabel('Collective Variable (nm)', fontsize=12, fontweight='bold')
        ax1.set_ylabel('PMF (kJ/mol)', fontsize=12, fontweight='bold')
        ax1.set_title('Potential of Mean Force', fontsize=13, fontweight='bold')
        ax1.legend(fontsize=10)
        ax1.grid(True, alpha=0.3)
        
        # Panel 2: Dependencia de k con barrera
        ax2 = axes[0, 1]
        
        barriers = np.linspace(5, 40, 100)
        rates = [self.calculate_rate(b, friction_factor=500) for b in barriers]
        
        ax2.semilogy(barriers, rates, 'b-', linewidth=2.5)
        ax2.semilogy(barrier, rate, 'ro', markersize=15, 
                    label=f'Sistema actual\n({barrier:.1f} kJ/mol)', zorder=5)
        
        ax2.axvline(barrier, color='red', linestyle='--', alpha=0.5)
        ax2.axhline(rate, color='red', linestyle='--', alpha=0.5)
        
        ax2.set_xlabel('Barrera ΔG‡ (kJ/mol)', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Tasa k (s⁻¹)', fontsize=12, fontweight='bold')
        ax2.set_title('Dependencia de k con Barrera\n(Eyring equation)', 
                     fontsize=13, fontweight='bold')
        ax2.legend(fontsize=10)
        ax2.grid(True, alpha=0.3, which='both')
        
        # Panel 3: Tiempos característicos
        ax3 = axes[1, 0]
        
        friction_factors = [1, 10, 100, 500, 1000, 5000]
        times = [self.mean_time(self.calculate_rate(barrier, friction_factor=f)) 
                for f in friction_factors]
        
        colors = plt.cm.viridis(np.linspace(0, 1, len(friction_factors)))
        
        bars = ax3.barh(range(len(friction_factors)), 
                       [t * 1e9 for t in times],  # Convertir a ns
                       color=colors)
        
        ax3.set_yticks(range(len(friction_factors)))
        ax3.set_yticklabels([f'FF = {f}' for f in friction_factors])
        ax3.set_xlabel('Tiempo medio <t> (nanosegundos)', fontsize=12, fontweight='bold')
        ax3.set_ylabel('Factor de Fricción', fontsize=12, fontweight='bold')
        ax3.set_title('Efecto de Fricción en Tiempo de Transición', 
                     fontsize=13, fontweight='bold')
        ax3.set_xscale('log')
        ax3.grid(True, alpha=0.3, axis='x')
        
        # Añadir etiquetas
        for i, (bar, t) in enumerate(zip(bars, times)):
            ax3.text(bar.get_width(), bar.get_y() + bar.get_height()/2,
                    f' {self.format_time(t)}',
                    va='center', fontsize=9)
        
        # Panel 4: Sensibilidad a mutaciones
        ax4 = axes[1, 1]
        
        delta_gs = np.linspace(-15, 15, 31)
        fold_changes = [self.calculate_rate(barrier + dg, friction_factor=500) / rate 
                       for dg in delta_gs]
        
        ax4.semilogy(delta_gs, fold_changes, 'b-', linewidth=2.5)
        ax4.axhline(1.0, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
        ax4.axvline(0.0, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
        
        # Regiones de interés
        ax4.fill_between([-15, -5], [1e-5, 1e-5], [1e5, 1e5], 
                        alpha=0.1, color='green', label='Ganancia función')
        ax4.fill_between([5, 15], [1e-5, 1e-5], [1e5, 1e5],
                        alpha=0.1, color='red', label='Pérdida función')
        
        ax4.set_xlabel('Cambio en barrera ΔΔG (kJ/mol)', fontsize=12, fontweight='bold')
        ax4.set_ylabel('Cambio en tasa (k_mutante / k_WT)', fontsize=12, fontweight='bold')
        ax4.set_title('Sensibilidad a Mutaciones', fontsize=13, fontweight='bold')
        ax4.legend(fontsize=10)
        ax4.grid(True, alpha=0.3, which='both')
        
        plt.tight_layout()
        plt.savefig(output, dpi=300, bbox_inches='tight')
        print(f"\n✓ Figura guardada: {output}")
        plt.close()


def main():
    """Función principal"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Calcula tasas cinéticas desde PMF (Eyring equation)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:

  # Analizar PMF completo
  python estimate_rates_from_pmf.py mbar_results.txt
  
  # Especificar temperatura
  python estimate_rates_from_pmf.py mbar_results.txt --temp 298
  
  # Solo calcular tasa para barrera conocida
  python estimate_rates_from_pmf.py --barrier 25.0
  
  # Con factor de fricción específico
  python estimate_rates_from_pmf.py --barrier 25.0 --friction 500

Interpretación de resultados:
  - k > 10⁶ s⁻¹: Transición muy rápida (ns)
  - k ~ 10³-10⁶ s⁻¹: Transición rápida (μs)
  - k ~ 10⁰-10³ s⁻¹: Transición media (ms)
  - k < 1 s⁻¹: Transición lenta (s o más)

IMPORTANTE:
  Estas son estimaciones de ORDEN DE MAGNITUD.
  Para cinética precisa, usar métodos especializados.
        """
    )
    
    parser.add_argument('pmf_file', nargs='?', 
                       help='Archivo con PMF (formato: CV PMF)')
    parser.add_argument('--barrier', type=float,
                       help='Barrera energética en kJ/mol (si no hay archivo PMF)')
    parser.add_argument('--temp', type=float, default=310.0,
                       help='Temperatura en Kelvin (default: 310 K)')
    parser.add_argument('--friction', type=float, default=500.0,
                       help='Factor de fricción (default: 500)')
    
    args = parser.parse_args()
    
    # Crear calculadora
    calc = EyringRateCalculator(temperature=args.temp)
    
    if args.pmf_file:
        # Análisis completo de PMF
        if not Path(args.pmf_file).exists():
            print(f"✗ Error: Archivo no encontrado: {args.pmf_file}")
            return 1
        
        try:
            results = calc.analyze_pmf(args.pmf_file)
            print(f"\n{'='*70}")
            print("✓ Análisis completado")
            print(f"{'='*70}")
            return 0
        except Exception as e:
            print(f"✗ Error al analizar PMF: {e}")
            import traceback
            traceback.print_exc()
            return 1
    
    elif args.barrier:
        # Cálculo simple de tasa
        print(f"\n{'='*70}")
        print(f"CÁLCULO DE TASA PARA BARRERA ΔG‡ = {args.barrier} kJ/mol")
        print(f"{'='*70}")
        
        k = calc.calculate_rate(args.barrier, friction_factor=args.friction)
        t_mean = calc.mean_time(k)
        t_half = calc.half_life(k)
        
        print(f"\nResultados:")
        print(f"  Tasa: k = {k:.3e} s⁻¹")
        print(f"  Tiempo medio: <t> = {calc.format_time(t_mean)}")
        print(f"  Vida media: t₁/₂ = {calc.format_time(t_half)}")
        
        # Sensibilidad
        print(f"\nSensibilidad:")
        for delta_g in [-5, -2, 2, 5]:
            new_k = calc.calculate_rate(args.barrier + delta_g, 
                                       friction_factor=args.friction)
            fold = new_k / k
            print(f"  ΔΔG = {delta_g:+d} kJ/mol → k_new/k = {fold:.1f}×")
        
        return 0
    
    else:
        parser.print_help()
        return 1


if __name__ == '__main__':
    sys.exit(main())
