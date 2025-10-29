# Umbrella Sampling vs Metadinámica: Validación Cruzada

## ✅ TU PUNTO ES CORRECTO

**Tienes razón en TODO**:

1. ✅ **Mismo motor MD (OpenMM)** - mismos parámetros físicos
2. ✅ **Mismo objetivo** - calcular PMF (energía libre ΔG)
3. ✅ **Son compatibles y comparables** - deben dar el MISMO resultado
4. ✅ **Ambas están "sesgadas"** - solo difiere el tipo de bias
5. ✅ **NO se analizan por separado** - se COMPARAN para validación

---

## 🔬 Aclaración: Son DOS ALGORITMOS para el MISMO CÁLCULO

### MISMO Sistema Físico

```python
# ════════════════════════════════════════════════════════
# CONFIGURACIÓN BASE (IDÉNTICA EN AMBOS MÉTODOS)
# ════════════════════════════════════════════════════════
from openmm import *
from openmm.app import *

# Sistema
pdb = PDBFile('equilibrated.pdb')           # Mismo PDB
forcefield = ForceField('amber14-all.xml')  # Mismo forcefield

# Parámetros MD (IDÉNTICOS)
temperature = 310*kelvin                    # 37°C
pressure = 1.0*bar                          # 1 atm
friction = 1.0/picosecond                   # Langevin
timestep = 2.0*femtoseconds                 # 2 fs
ionic_strength = 0.163*molar                # PBS buffer

# Integrador (IDÉNTICO)
integrator = LangevinMiddleIntegrator(temperature, friction, timestep)

# Barostato (IDÉNTICO)
system.addForce(MonteCarloBarostat(pressure, temperature))

# ════════════════════════════════════════════════════════
# ÚNICA DIFERENCIA: TIPO DE BIAS AGREGADO
# ════════════════════════════════════════════════════════
```

### Umbrella Sampling: Bias ESTÁTICO Armónico

```python
# BIAS TIPO 1: Armónico (mantiene CV cerca de target)
cv = CustomCentroidBondForce(2, 'distance(g1,g2)')
cv.addGroup(kinase_atoms)
cv.addGroup(cterm_atoms)
cv.addBond([0, 1])

# Agregar potencial de resorte
umbrella_force = CustomCVForce('0.5*k*(cv-cv0)^2')
umbrella_force.addCollectiveVariable('cv', cv)
umbrella_force.addGlobalParameter('k', 1000.0*kilojoules_per_mole/nanometer**2)
umbrella_force.addGlobalParameter('cv0', target_distance)

system.addForce(umbrella_force)

# Ejecutar 20 ventanas × 100 ns = 2 μs total
simulation.step(50000000)  # 100 ns por ventana
```

**Hamiltoniano**:
$$H_{umbrella} = H_{MD} + \frac{k}{2}(CV - CV_0)^2$$

**Bias**: ESTÁTICO - siempre empuja hacia `cv0`

### Metadinámica: Bias DINÁMICO Gaussiano

```python
# BIAS TIPO 2: Gaussiano acumulativo (empuja a explorar)
from openmm import *

# Definir variable colectiva (MISMA que umbrella)
cv = CustomCentroidBondForce(2, 'distance(g1,g2)')
cv.addGroup(kinase_atoms)
cv.addGroup(cterm_atoms)

# Configurar metadinámica
metad = BiasVariable(cv, minValue=2.0*nanometers, maxValue=4.0*nanometers,
                     biasWidth=0.05*nanometers)  # σ gaussiano
meta_force = Metadynamics(system, [metad], 
                          temperature=310*kelvin,
                          biasFactor=5,           # Well-tempered
                          height=1.0*kilojoules_per_mole,
                          frequency=500)          # Cada 500 steps

# Ejecutar 1 trayectoria × 500 ns continua
simulation.step(250000000)  # 500 ns
```

**Hamiltoniano**:
$$H_{metad} = H_{MD} + \sum_{i=1}^{N} h \cdot e^{-\frac{(CV - CV_i)^2}{2\sigma^2}}$$

**Bias**: DINÁMICO - crece con el tiempo, empuja donde ya estuvo

---

## 🎯 Objetivo IDÉNTICO: Calcular PMF

### Potential of Mean Force (PMF)

**Definición termodinámica** (misma para ambos):
$$PMF(CV) = -k_B T \ln P(CV) + const$$

Donde:
- $P(CV)$ = Probabilidad de observar el valor CV
- $k_B T$ = Energía térmica (2.58 kJ/mol a 310 K)

### Umbrella → PMF via MBAR

```python
# Post-procesamiento: Combinar 20 ventanas
from pymbar import MBAR

# Cargar datos de 20 ventanas
u_kn = []  # Energías en cada ventana
N_k = []   # Número de muestras
for window in range(20):
    cv_values = np.loadtxt(f'window_{window:02d}/cv_values.dat')
    u_kn.append(compute_bias_energies(cv_values, k=1000, cv0=targets[window]))
    N_k.append(len(cv_values))

# MBAR resuelve ecuaciones auto-consistentes
mbar = MBAR(u_kn, N_k)
pmf_umbrella, uncertainties = mbar.computePMF()
```

**Output**: `pmf_umbrella.csv` con columnas `[cv, pmf, uncertainty]`

### Metadinámica → PMF via Bias Acumulado

```python
# Post-procesamiento: Sumar gaussianas
gaussians = []  # Lista de (cv_center, height, width)
with open('HILLS', 'r') as f:  # OpenMM guarda aquí
    for line in f:
        t, cv, height, sigma = line.split()
        gaussians.append((float(cv), float(height), float(sigma)))

# PMF = -1 × bias acumulado (en el límite de convergencia)
cv_grid = np.linspace(2.0, 4.0, 100)
pmf_metad = -1 * sum_gaussians(cv_grid, gaussians)
pmf_metad -= pmf_metad.min()  # Normalizar mínimo a 0
```

**Output**: `pmf_metad.csv` con columnas `[cv, pmf]`

---

## ✅ Validación Cruzada: DEBEN Coincidir

### Workflow Completo

```bash
# ═══════════════════════════════════════════════════════════
# PASO 1: Ejecutar Umbrella Sampling (Método 1)
# ═══════════════════════════════════════════════════════════
cd Chronosfold/WNK
bash run_complete_pipeline.sh hpc

# Output: pmf_analysis/pmf_umbrella.csv
#   cv (nm)    pmf (kJ/mol)    uncertainty (kJ/mol)
#   2.00       15.3            0.8
#   2.10       12.1            0.6
#   ...
#   4.00       22.7            1.2

# ═══════════════════════════════════════════════════════════
# PASO 2: Ejecutar Metadinámica (Método 2)
# ═══════════════════════════════════════════════════════════
python run_metadynamics.py --cv-min 2.0 --cv-max 4.0 --time 500

# Output: metadynamics/pmf_metad.csv
#   cv (nm)    pmf (kJ/mol)
#   2.00       14.8
#   2.10       11.9
#   ...
#   4.00       23.1

# ═══════════════════════════════════════════════════════════
# PASO 3: Comparar PMFs (Validación Cruzada)
# ═══════════════════════════════════════════════════════════
python compare_pmfs.py \
    --method1 pmf_analysis/pmf_umbrella.csv \
    --method2 metadynamics/pmf_metad.csv \
    --output comparison.png

# Output esperado:
# ✓ RMSD between PMFs: 1.4 kJ/mol
# ✓ Barrier height difference: 0.6 kJ/mol
# ✓ Minimum location agreement: 0.02 nm
# ✓ VALIDATION PASSED: Methods converged to same PMF
```

### Criterios de Convergencia

| Métrica | Buen Acuerdo | Necesita Más Muestreo |
|---------|--------------|----------------------|
| **RMSD global** | < 2 kJ/mol | > 5 kJ/mol |
| **Diferencia en barrera** | < 3 kJ/mol | > 8 kJ/mol |
| **Ubicación mínimo** | < 0.05 nm | > 0.2 nm |
| **Forma general** | Correlación > 0.95 | < 0.85 |

**Si NO coinciden** → Uno (o ambos) no convergieron:
- Umbrella: Aumentar tiempo por ventana (200 ns)
- Metadinámica: Aumentar tiempo total (1 μs) o reducir `height` (convergencia más lenta pero precisa)

---

## 🔄 ¿Por Qué Hacer Ambos?

### Ventajas de Cada Método

| Aspecto | Umbrella Sampling | Metadinámica |
|---------|-------------------|--------------|
| **Convergencia** | Predecible (cada ventana independiente) | Depende de parámetros (height, σ) |
| **Paralelización** | ✅ Excelente (20 jobs independientes) | ⚠️ Limitada (1 trayectoria continua) |
| **Tiempo de muestreo** | Eficiente (enfocado en cada región) | Puede ser lento (exploración global) |
| **Detección de estados** | Requiere conocer rango CV a priori | ✅ Descubre automáticamente |
| **Barreras altas** | ✅ Excelente (sampling forzado) | ⚠️ Puede quedarse atrapado |
| **Setup inicial** | Más complejo (definir ventanas) | Más simple (solo run) |

### Casos de Uso Recomendados

**Umbrella es mejor si**:
- ✅ Conoces el rango de CV relevante (2.0-4.0 nm)
- ✅ Tienes HPC con muchos cores (paralelizar ventanas)
- ✅ Barreras energéticas altas (>20 kJ/mol)
- ✅ Necesitas precisión alta (incertidumbres <2 kJ/mol)

**Metadinámica es mejor si**:
- ✅ NO conoces el rango de CV a priori (exploración)
- ✅ Múltiples mínimos desconocidos
- ✅ Quieres ver trayectorias continuas (transiciones)
- ✅ Testing rápido (una sola simulación)

### **MEJOR: Hacer Ambos** (Validación Cruzada)

```
┌─────────────────────────────────────────┐
│  GOLD STANDARD: Ambos Métodos Coinciden │
└─────────────────────────────────────────┘

Umbrella PMF ≈ Metadinámica PMF
    ↓
Confianza alta en resultado
    ↓
Publicar con ambos métodos como validación
```

**En papers de alto impacto**: Siempre muestran 2+ métodos convergiendo al mismo resultado.

---

## 🧪 Implementación Práctica

### Agregar Metadinámica al Pipeline

```python
# ═══════════════════════════════════════════════════════════
# run_metadynamics.py - Nueva adición al pipeline
# ═══════════════════════════════════════════════════════════
#!/usr/bin/env python3
"""
Metadinámica Well-Tempered para WNK1 C-terminal
Mismo sistema que umbrella sampling, diferente algoritmo de bias
"""

from openmm import *
from openmm.app import *
import numpy as np

# Cargar sistema (MISMO que umbrella)
pdb = PDBFile('prepared_system/equilibrated.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
system = forcefield.createSystem(pdb.topology, 
                                 nonbondedMethod=PME,
                                 nonbondedCutoff=1.0*nanometers,
                                 constraints=HBonds)

# Parámetros MD (IDÉNTICOS a umbrella)
temperature = 310*kelvin
pressure = 1.0*bar
friction = 1.0/picosecond
timestep = 2.0*femtoseconds

integrator = LangevinMiddleIntegrator(temperature, friction, timestep)
system.addForce(MonteCarloBarostat(pressure, temperature))

# Definir CV (MISMA que umbrella)
cv = CustomCentroidBondForce(2, 'distance(g1,g2)')
cv.addGroup(kinase_atoms)  # Cargar desde atom_groups.txt
cv.addGroup(cterm_atoms)
system.addForce(cv)

# METADINÁMICA: Configurar bias gaussiano
from openmm import *
metad = BiasVariable(cv, 
                     minValue=2.0*nanometers, 
                     maxValue=4.0*nanometers,
                     biasWidth=0.05*nanometers,  # σ = 0.05 nm
                     periodic=False)

meta_force = Metadynamics(system, [metad],
                          temperature=temperature,
                          biasFactor=5,            # Well-tempered (ΔT = 1240 K)
                          height=1.0*kilojoules_per_mole,
                          frequency=500,           # Agregar gaussiana cada 1 ps
                          saveFrequency=5000,      # Guardar HILLS cada 10 ps
                          biasDir='metadynamics')

# Ejecutar 500 ns
simulation = Simulation(pdb.topology, system, integrator, Platform.getPlatformByName('CPU'))
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()

print("Corriendo metadinámica (500 ns)...")
simulation.step(250000000)  # 250M steps × 2 fs = 500 ns

print("Guardando estado final...")
state = simulation.context.getState(getPositions=True)
PDBFile.writeFile(pdb.topology, state.getPositions(), 
                  open('metadynamics/final.pdb', 'w'))

# Post-procesamiento: Calcular PMF
meta_force.getFreeEnergy(metad, saveToFile='metadynamics/pmf_metad.csv')
print("✓ PMF guardado en metadynamics/pmf_metad.csv")
```

### Script de Comparación

```python
# ═══════════════════════════════════════════════════════════
# compare_pmfs.py - Validación cruzada
# ═══════════════════════════════════════════════════════════
#!/usr/bin/env python3
"""Comparar PMFs de umbrella y metadinámica"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.stats import spearmanr

def compare_pmfs(umbrella_file, metad_file, output='comparison.png'):
    """Comparar y validar dos PMFs"""
    
    # Cargar datos
    umbrella = pd.read_csv(umbrella_file)
    metad = pd.read_csv(metad_file)
    
    # Interpolar a misma grilla
    cv_grid = np.linspace(2.0, 4.0, 200)
    pmf_u_interp = interp1d(umbrella['cv'], umbrella['pmf'], 
                            kind='cubic', fill_value='extrapolate')
    pmf_m_interp = interp1d(metad['cv'], metad['pmf'], 
                            kind='cubic', fill_value='extrapolate')
    
    pmf_u = pmf_u_interp(cv_grid)
    pmf_m = pmf_m_interp(cv_grid)
    
    # Alinear (mínimo a 0)
    pmf_u -= pmf_u.min()
    pmf_m -= pmf_m.min()
    
    # Métricas de comparación
    rmsd = np.sqrt(np.mean((pmf_u - pmf_m)**2))
    max_diff = np.abs(pmf_u - pmf_m).max()
    corr, _ = spearmanr(pmf_u, pmf_m)
    
    # Diferencias en características
    barrier_u = pmf_u.max()
    barrier_m = pmf_m.max()
    min_loc_u = cv_grid[pmf_u.argmin()]
    min_loc_m = cv_grid[pmf_m.argmin()]
    
    # Plot comparativo
    fig, axes = plt.subplots(2, 1, figsize=(10, 10))
    
    # Panel 1: Overlay
    ax1 = axes[0]
    ax1.plot(cv_grid, pmf_u, 'b-', linewidth=2, label='Umbrella Sampling')
    ax1.plot(cv_grid, pmf_m, 'r--', linewidth=2, label='Metadinámica')
    if 'uncertainty' in umbrella.columns:
        unc_interp = interp1d(umbrella['cv'], umbrella['uncertainty'], 
                             kind='cubic', fill_value='extrapolate')
        unc = unc_interp(cv_grid)
        ax1.fill_between(cv_grid, pmf_u-unc, pmf_u+unc, alpha=0.3, color='blue')
    ax1.set_xlabel('Distancia C-terminal (nm)', fontweight='bold', fontsize=12)
    ax1.set_ylabel('PMF (kJ/mol)', fontweight='bold', fontsize=12)
    ax1.set_title('Validación Cruzada: Umbrella vs Metadinámica', 
                  fontweight='bold', fontsize=14)
    ax1.legend(fontsize=11)
    ax1.grid(alpha=0.3)
    
    # Anotar métricas
    textstr = f'RMSD: {rmsd:.2f} kJ/mol\n'
    textstr += f'Max diff: {max_diff:.2f} kJ/mol\n'
    textstr += f'Correlation: {corr:.3f}'
    ax1.text(0.02, 0.98, textstr, transform=ax1.transAxes,
            verticalalignment='top', bbox=dict(boxstyle='round', 
            facecolor='wheat', alpha=0.8), fontsize=10)
    
    # Panel 2: Diferencia
    ax2 = axes[1]
    diff = pmf_u - pmf_m
    ax2.plot(cv_grid, diff, 'k-', linewidth=2)
    ax2.axhline(0, color='gray', linestyle='--', alpha=0.5)
    ax2.fill_between(cv_grid, -2, 2, alpha=0.2, color='green', 
                    label='Buen acuerdo (<2 kJ/mol)')
    ax2.set_xlabel('Distancia C-terminal (nm)', fontweight='bold', fontsize=12)
    ax2.set_ylabel('Δ PMF (Umbrella - Metad) [kJ/mol]', fontweight='bold', fontsize=12)
    ax2.set_title('Diferencia entre Métodos', fontweight='bold', fontsize=14)
    ax2.legend(fontsize=10)
    ax2.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"✓ Comparación guardada: {output}")
    
    # Reporte
    print("\n" + "="*60)
    print("  VALIDACIÓN CRUZADA: Umbrella vs Metadinámica")
    print("="*60)
    print(f"\nRMSD global:               {rmsd:.2f} kJ/mol")
    print(f"Máxima diferencia:         {max_diff:.2f} kJ/mol")
    print(f"Correlación (Spearman):    {corr:.3f}")
    print(f"\nBarrera Umbrella:          {barrier_u:.2f} kJ/mol")
    print(f"Barrera Metadinámica:      {barrier_m:.2f} kJ/mol")
    print(f"Diferencia en barrera:     {abs(barrier_u - barrier_m):.2f} kJ/mol")
    print(f"\nMínimo Umbrella:           {min_loc_u:.3f} nm")
    print(f"Mínimo Metadinámica:       {min_loc_m:.3f} nm")
    print(f"Diferencia en ubicación:   {abs(min_loc_u - min_loc_m):.3f} nm")
    
    # Veredicto
    print("\n" + "-"*60)
    if rmsd < 2.0 and abs(barrier_u - barrier_m) < 3.0:
        print("✓ VALIDACIÓN EXITOSA: Métodos convergen al mismo PMF")
        print("  Confianza alta en resultado termodinámico")
    elif rmsd < 5.0:
        print("⚠️  VALIDACIÓN PARCIAL: Acuerdo razonable pero no óptimo")
        print("  Considerar aumentar tiempo de muestreo en ambos métodos")
    else:
        print("✗ VALIDACIÓN FALLIDA: Métodos no convergen")
        print("  Revisar convergencia de cada método individualmente")
    print("="*60 + "\n")
    
    return rmsd, corr

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--method1', required=True, help='PMF umbrella sampling')
    parser.add_argument('--method2', required=True, help='PMF metadinámica')
    parser.add_argument('--output', default='comparison.png', help='Output plot')
    args = parser.parse_args()
    
    compare_pmfs(args.method1, args.method2, args.output)
```

---

## 📊 Ejemplo de Validación Exitosa

### Output Esperado

```
═══════════════════════════════════════════════════════════
  VALIDACIÓN CRUZADA: Umbrella vs Metadinámica
═══════════════════════════════════════════════════════════

RMSD global:               1.4 kJ/mol
Máxima diferencia:         3.1 kJ/mol
Correlación (Spearman):    0.982

Barrera Umbrella:          18.3 kJ/mol
Barrera Metadinámica:      17.7 kJ/mol
Diferencia en barrera:     0.6 kJ/mol

Mínimo Umbrella:           2.34 nm
Mínimo Metadinámica:       2.31 nm
Diferencia en ubicación:   0.03 nm

───────────────────────────────────────────────────────────
✓ VALIDACIÓN EXITOSA: Métodos convergen al mismo PMF
  Confianza alta en resultado termodinámico
═══════════════════════════════════════════════════════════
```

### Plot Generado

El script genera un plot 2-panel:
1. **Overlay**: Ambos PMFs superpuestos (deben coincidir visualmente)
2. **Diferencia**: Δ PMF vs CV (debe oscilar cerca de 0 ± 2 kJ/mol)

---

## ✅ Conclusión

### **TIENES RAZÓN**:

1. ✅ **Mismo motor MD** (OpenMM con parámetros idénticos)
2. ✅ **Mismo objetivo** (calcular PMF del mismo sistema)
3. ✅ **Son compatibles** (deben dar mismo resultado)
4. ✅ **NO se analizan separado** - se COMPARAN para validación cruzada
5. ✅ **Ambos tienen bias** - solo difiere el algoritmo (armónico vs gaussiano)

### Workflow Recomendado

```
PASO 1: Umbrella Sampling (20 ventanas × 100 ns)
        ↓
PASO 2: Metadinámica (1 trayectoria × 500 ns)
        ↓
PASO 3: Comparar PMFs (deben coincidir RMSD < 2 kJ/mol)
        ↓
PASO 4: Si coinciden → Publicar ambos como validación
        Si NO coinciden → Aumentar muestreo y revisar
```

### En Papers

**Sección de Resultados**:
> "Para validar la robustez de nuestros cálculos de energía libre, empleamos dos métodos independientes: umbrella sampling con MBAR y metadinámica well-tempered. **Ambos métodos convergieron al mismo perfil de PMF** (RMSD = 1.4 kJ/mol), con una barrera de disociación del C-terminal de 18.0 ± 0.5 kJ/mol."

**Figura Suplementaria**:
- Panel A: PMF umbrella (20 ventanas, 2 μs total)
- Panel B: PMF metadinámica (500 ns continua)
- Panel C: Overlay de ambos (demostrando convergencia)
- Panel D: Diferencia (Δ PMF ≈ 0)

---

**Próximo paso**: ¿Quieres que implemente el script de metadinámica (`run_metadynamics.py`) para agregar al pipeline?
