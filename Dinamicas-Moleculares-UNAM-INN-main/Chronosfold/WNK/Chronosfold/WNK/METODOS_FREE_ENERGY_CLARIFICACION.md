# Clarificación: Métodos de Cálculo de Free Energy (ΔG)

**Autor**: Asistente IA - Preparación para reunión con la doctora  
**Fecha**: 2025  
**Propósito**: Corregir confusión conceptual entre enhanced sampling y métodos alquímicos

---

## 🚨 CORRECCIÓN CRÍTICA

### Tu argumento original contenía un error conceptual:

> "si quisiéramos calcular delta g **mejor que umbrella** usaríamos **métodos alquímicos en metadinámica**"

**PROBLEMA DETECTADO**:  
Esto mezcla dos categorías de métodos **completamente diferentes** que responden preguntas científicas distintas.

---

## 1. Cómo Metadinámica Calcula ΔG

### Algoritmo de Metadinámica

Metadinámica usa un **bias gaussiano acumulativo** que se añade durante la simulación:

$$V_{bias}(s, t) = \sum_{t'=0}^{t} h \cdot \exp\left(-\frac{(s(t) - s(t'))^2}{2\sigma^2}\right)$$

Donde:
- $s$ = Collective Variable (CV)
- $h$ = Altura del gaussiano (ej. 1.2 kJ/mol)
- $\sigma$ = Ancho del gaussiano (ej. 0.1 nm)
- $t'$ = Pasos de tiempo donde se depositó gaussiano

### Teorema Fundamental de Metadinámica

En el límite de convergencia ($t \rightarrow \infty$):

$$V_{bias}(s) \rightarrow -F(s) + \text{constante}$$

Donde $F(s)$ es el **Potential of Mean Force (PMF)**.

**Interpretación física**:
1. El sistema explora con bias acumulativo
2. Cuando visita región ya explorada → bias la empuja a salir
3. Eventualmente, el bias "rellena" todos los valles del PMF
4. En convergencia: **bias = -PMF** (perfil invertido)

### Extracción Práctica de ΔG

```python
import numpy as np

def calculate_pmf_from_metad(hills_file, cv_min, cv_max, bins=100):
    """
    Calcula PMF desde archivo de HILLS de metadinámica.
    
    Parámetros:
    -----------
    hills_file : str
        Archivo con gaussianos depositados (time, CV, height, sigma)
    cv_min, cv_max : float
        Rango de la collective variable
    bins : int
        Resolución del PMF
    
    Retorna:
    --------
    cv_values : ndarray
        Coordenadas de CV
    pmf : ndarray
        Free energy en kJ/mol
    """
    
    # 1. Leer HILLS (archivo de metadinámica)
    hills = np.loadtxt(hills_file)
    times = hills[:, 0]
    cv_centers = hills[:, 1]
    heights = hills[:, 2]
    sigmas = hills[:, 3]
    
    # 2. Crear grid de CV
    cv_values = np.linspace(cv_min, cv_max, bins)
    bias = np.zeros(bins)
    
    # 3. Sumar todos los gaussianos
    for i, cv in enumerate(cv_values):
        for j in range(len(hills)):
            bias[i] += heights[j] * np.exp(
                -(cv - cv_centers[j])**2 / (2 * sigmas[j]**2)
            )
    
    # 4. PMF = -V_bias (teorema de metadinámica)
    pmf = -bias
    
    # 5. Normalizar: mínimo = 0
    pmf -= pmf.min()
    
    # 6. Calcular barrera de activación
    delta_g = pmf.max() - pmf.min()
    
    print(f"Barrera de activación (ΔG): {delta_g:.2f} kJ/mol")
    
    return cv_values, pmf

# Uso:
cv, pmf = calculate_pmf_from_metad(
    'HILLS', 
    cv_min=2.0, 
    cv_max=4.0, 
    bins=200
)
```

### Well-Tempered Metadynamics (Recomendado)

En la práctica moderna se usa **well-tempered metadynamics**, que reduce la altura del gaussiano con el tiempo:

$$h(t) = h_0 \cdot \exp\left(-\frac{V_{bias}(s(t), t)}{k_B \Delta T}\right)$$

Donde $\Delta T$ es el "bias factor" (típicamente 10-15 para proteínas).

**Ventaja**: Converge más suavemente, evita sobre-muestreo.

---

## 2. Tres Clases de Métodos de Free Energy

### A. Enhanced Sampling (Umbrella, Metadinámica)

**Pregunta científica**: *"¿Cuál es el costo energético de una transición conformacional?"*

| Método | Tipo de Bias | Cálculo de ΔG |
|--------|--------------|---------------|
| **Umbrella Sampling** | Bias armónico estático: $k(s - s_0)^2$ | WHAM/MBAR sobre ventanas |
| **Metadinámica** | Bias gaussiano acumulativo | $-V_{bias}$ en convergencia |
| **ABF** | Bias adaptativo en fuerza | Integración de fuerza media |
| **Steered MD** | Bias de pulling | Jarzynski equality |

**Características comunes**:
- Aplican bias en el espacio de **collective variables (CV)**
- Calculan **PMF(CV)** = perfil de energía libre vs. coordenada de reacción
- Aceleran cruce de barreras energéticas
- **Mismo engine MD** (OpenMM, GROMACS, AMBER) con parámetros idénticos
- **Solo difiere el algoritmo de bias**

**Ejemplo WNK1**:
- CV = Distancia Cα(C-term) - Cα(N-term)
- Umbrella: 20 ventanas × 100 ns = 2 μs
- Metadinámica: 500 ns continuo con gaussianos cada 1 ps
- **Ambos calculan el mismo ΔG** (barrera conformacional)

### B. Alchemical Methods (FEP, TI, BAR)

**Pregunta científica**: *"¿Cuánto cambia ΔG si transformo químicamente el sistema?"*

| Método | Parámetro de Transformación | Cálculo de ΔG |
|--------|------------------------------|---------------|
| **FEP** (Free Energy Perturbation) | $\lambda$: 0 (estado A) → 1 (estado B) | Exponential averaging: $\langle e^{-\Delta U/k_B T} \rangle$ |
| **TI** (Thermodynamic Integration) | $\lambda$: 0 → 1 | Integración: $\int_0^1 \langle \partial U/\partial \lambda \rangle d\lambda$ |
| **BAR** (Bennett Acceptance Ratio) | $\lambda$: ventanas intermedias | Optimización de solapamiento |

**Transformaciones típicas**:
- **Mutación**: Ala → Gly (cambiar parámetros de residuo)
- **Unión de ligando**: Sin ligando → Con ligando (aparecer/desaparecer átomos)
- **Solvatación**: Vacío → Agua (calcular ΔG_solv)

**Características clave**:
- Usan parámetro **λ** (NO es una CV física)
- **Modifican el Hamiltoniano** del sistema durante la simulación
- Calculan **ΔΔG** de transformación química
- Requieren múltiples ventanas de λ (típicamente 10-20)

**Ejemplo de mutación S1261A en WNK1**:
```python
# λ = 0: Wildtype (Ser1261)
# λ = 0.5: Estado intermedio (partículas fantasma)
# λ = 1: Mutante (Ala1261)

# ΔΔG = ΔG_mutante - ΔG_wildtype
# Calculado vía TI:
ddG = integrate(dU/dλ, λ=0 to 1)
```

### C. Standard MD (drMD, MD Estándar)

**Pregunta científica**: *"¿Qué conformaciones explora el sistema espontáneamente?"*

**Características**:
- **Sin bias artificial**
- Solo explora según Hamiltoniano natural
- **No calcula ΔG cuantitativamente**
- Útil para: muestreo, clustering, análisis RMSD, visualización

**Limitaciones**:
- No puede cruzar barreras altas (>15 kJ/mol) en tiempos razonables
- No proporciona valores de ΔG directamente
- Requiere microsegundos para eventos raros

---

## 3. Comparación Directa: ¿Qué ΔG Calcula Cada Método?

| Pregunta Científica | Método Apropiado | ΔG Calculado | Unidades | Tiempo Típico |
|---------------------|------------------|--------------|----------|---------------|
| ¿Barrera conformacional abierto→cerrado? | **Umbrella / Metadinámica** | PMF(distancia) | kJ/mol | 1-5 μs total |
| ¿Costo de mutación Ser→Ala? | **FEP / TI / BAR** | ΔΔG(mutación) | kJ/mol | 50-200 ns |
| ¿Afinidad de unión ligando? | **FEP / TI / BAR** | ΔG(binding) | kJ/mol | 100-500 ns |
| ¿Exploración sin cuantificación? | **drMD / MD estándar** | Ninguno (cualitativo) | - | 100-1000 ns |
| ¿Validación cruzada de PMF? | **Umbrella + Metadinámica** | PMF(CV) × 2 | kJ/mol | 2× tiempo |

---

## 4. Corrección del Argumento para la Doctora

### ❌ Versión Original (INCORRECTA)

> "drMD no tiene umbrella, no podemos correrlo para correr delta g porque no calcula delta g como tal como umbrella, **y si quisiéramos calcular delta g mejor que umbrella usaríamos métodos alquímicos en metadinámica**"

**Errores conceptuales**:
1. ❌ "métodos alquímicos EN metadinámica" → No existe tal cosa
2. ❌ "mejor que umbrella" → Metadinámica no es "mejor", es alternativo (mismo resultado)
3. ❌ Implica que alquímicos calculan PMF conformacional → NO, calculan ΔΔG de transformaciones químicas

---

### ✅ Versión Corregida (CORRECTA)

**Parte 1 - Limitaciones de drMD** (mantener):
> "drMD es excelente para exploración conformacional y tiene validaciones robustas (pdbTriage, FirstAid, clustering automático), pero **no tiene métodos de enhanced sampling integrados** (umbrella, metadinámica, ABF). Por lo tanto, no puede cuantificar barreras de energía libre directamente."

**Parte 2 - Justificación de Umbrella Sampling** (nuevo):
> "Para cuantificar el **ΔG de la transición conformacional** del C-terminal de WNK1 (conformación abierta → cerrada), necesitamos un método de **enhanced sampling** que calcule el Potential of Mean Force (PMF). Seleccioné **umbrella sampling** por las siguientes razones:
> 
> 1. **CV bien definida**: La distancia Cα-Cα terminal es una coordenada de reacción natural
> 2. **Barreras altas**: Estimamos ~25 kJ/mol (requiere enhanced sampling)
> 3. **Paralelización eficiente**: 20 ventanas × 4 CPUs = uso óptimo de nuestro HPC de 48 cores
> 4. **Método establecido**: Extensa validación en literatura para transiciones conformacionales
> 5. **Post-análisis robusto**: MBAR proporciona errores estadísticos bien definidos"

**Parte 3 - Validación Opcional** (nuevo):
> "Como validación cruzada, podríamos implementar **metadinámica** (otro método de enhanced sampling) para verificar el mismo PMF. Ambos métodos usan el mismo engine de MD (OpenMM) con parámetros físicos idénticos, solo difieren en el algoritmo de bias (armónico vs. gaussiano acumulativo). Esperaríamos RMSD < 2 kJ/mol entre ambos perfiles."

**Parte 4 - Cuándo Usar Alquímicos** (nuevo):
> "Los **métodos alquímicos** (FEP, TI, BAR) son apropiados para preguntas científicas **diferentes**:
> 
> - ¿Cómo afecta la mutación S1261A a la estabilidad? → **TI/FEP** (ΔΔG de mutación)
> - ¿Cuál es la afinidad de unión de ATP? → **FEP** (ΔG_binding)
> - ¿Cómo cambia ΔG_solv al mutar? → **BAR** (ΔΔG_solv)
> 
> Estos métodos usan transformaciones químicas vía parámetro λ, NO calculan PMF de transiciones conformacionales. Para nuestro objetivo actual (barrera conformacional), **no son aplicables**."

---

## 5. Argumento Completo Listo para Defender

### Flujo Lógico

```
┌─────────────────────────────────────────────────────────┐
│ PREGUNTA CIENTÍFICA:                                    │
│ ¿Cuál es el costo energético (ΔG) de la transición     │
│ conformacional abierta→cerrada del C-terminal de WNK1?  │
└─────────────────────────────────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────┐
│ MÉTODO NECESARIO:                                       │
│ Enhanced Sampling (calcula PMF vs. CV)                  │
│                                                          │
│ Opciones:                                               │
│ • Umbrella Sampling ← SELECCIONADO                      │
│ • Metadinámica (validación opcional)                    │
│ • ABF (menos común en proteínas)                        │
│                                                          │
│ NO APLICABLE:                                           │
│ • Métodos alquímicos (FEP/TI/BAR) → Para mutaciones     │
│ • drMD estándar → Solo exploración cualitativa          │
└─────────────────────────────────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────┐
│ JUSTIFICACIÓN DE UMBRELLA:                              │
│                                                          │
│ 1. CV conocida (distancia Cα-Cα)                        │
│ 2. Barrera alta (~25 kJ/mol)                            │
│ 3. Paralelización HPC (48 cores, 20 ventanas)          │
│ 4. Convergencia verificable (MBAR bootstrap)           │
│ 5. Literatura establecida (>5000 papers)               │
└─────────────────────────────────────────────────────────┘
                          │
                          ▼
┌─────────────────────────────────────────────────────────┐
│ VALIDACIÓN:                                             │
│ • Convergencia: Error MBAR < 1 kJ/mol                   │
│ • Solapamiento: Histogramas adyacentes con overlap > 5% │
│ • Cross-check: Metadinámica (RMSD < 2 kJ/mol esperado) │
│ • Visualización: VideoSuite + 4-panel diagnostics      │
└─────────────────────────────────────────────────────────┘
```

---

## 6. Cuándo SÍ Usar Métodos Alquímicos

### Ejemplos de Preguntas Apropiadas

1. **Efecto de mutaciones puntuales**:
   ```
   Pregunta: ¿Cómo la mutación S1261A afecta la estabilidad de WNK1?
   Método: TI/FEP
   λ-pathway: Wildtype (λ=0) → Mutante (λ=1)
   Resultado: ΔΔG_fold = ΔG_mutante - ΔG_wildtype
   ```

2. **Afinidad de unión de ligando**:
   ```
   Pregunta: ¿Cuál es el ΔG de unión de ATP a WNK1?
   Método: FEP con doble transformación
   λ-pathway 1: ATP en agua → desaparece (ΔG_solv)
   λ-pathway 2: ATP en proteína → desaparece (ΔG_complex)
   Resultado: ΔG_binding = ΔG_complex - ΔG_solv
   ```

3. **Energía de solvatación**:
   ```
   Pregunta: ¿Cómo cambia ΔG_solv al sustituir Ser por Ala?
   Método: BAR
   λ-pathway: Ser solvatada → Ala solvatada
   Resultado: ΔΔG_solv
   ```

### Estructura de λ-Pathway

```python
# Ejemplo de FEP para mutación S1261A
lambda_windows = np.linspace(0, 1, 20)

for λ in lambda_windows:
    # Interpolación de parámetros
    charges = (1-λ) * charges_Ser + λ * charges_Ala
    vdw_params = (1-λ) * vdw_Ser + λ * vdw_Ala
    
    # Soft-core para evitar singularidades
    U_vdw = 4ε[(σ/(r^6 + α(1-λ)^2))^2 - (σ/(r^6 + α(1-λ)^2))]
    
    # Calcular dU/dλ
    dU_dlambda = U(λ+δ) - U(λ-δ) / (2δ)
    
    # Almacenar para TI
    store(λ, dU_dlambda)

# Integrar
ΔG = trapz(dU_dlambda, lambda_windows)
```

---

## 7. Resumen Ejecutivo para la Doctora

### Una Frase:
**"Umbrella sampling calcula el PMF de la transición conformacional del C-terminal de WNK1; métodos alquímicos (FEP/TI) son para transformaciones químicas (mutaciones, unión de ligandos), no para barreras conformacionales."**

### Tres Puntos Clave:

1. **Enhanced Sampling vs. Alchemical Methods**:
   - Enhanced sampling (umbrella, metadinámica): PMF(CV) para transiciones conformacionales
   - Alchemical (FEP, TI, BAR): ΔΔG para transformaciones químicas
   - **Son clases ortogonales**, no "mejor" una que otra

2. **Por Qué Umbrella para WNK1**:
   - CV bien definida (distancia Cα-Cα)
   - Paralelización óptima en HPC (48 cores)
   - Método gold standard para barreras conformacionales

3. **drMD vs. Enhanced Sampling**:
   - drMD: Exploración sin bias (no calcula ΔG)
   - Tiene validaciones excelentes (pdbTriage, FirstAid)
   - **NO es inferior**, simplemente responde pregunta diferente (cualitativa vs. cuantitativa)

---

## 8. Implementación de Metadinámica (Validación Opcional)

Si la doctora solicita validación cruzada, aquí está el código listo:

```python
# run_metadynamics_wnk.py
import openmm as mm
from openmm import app, unit
from openmm.app import PDBFile, Modeller, ForceField
import numpy as np

def setup_metadynamics(pdb_file, cv_atoms, output_prefix):
    """
    Configura well-tempered metadynamics para WNK1.
    
    Parámetros:
    -----------
    pdb_file : str
        PDB equilibrado (ej. 'wnk_pbs_equilibrated.pdb')
    cv_atoms : tuple
        (atom1_index, atom2_index) para distancia Cα-Cα
    output_prefix : str
        Prefijo para archivos de salida
    """
    
    # Cargar sistema
    pdb = PDBFile(pdb_file)
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*unit.nanometer,
        constraints=app.HBonds
    )
    
    # Definir Collective Variable (distancia)
    cv = mm.CustomBondForce('r')
    cv.addBond(cv_atoms[0], cv_atoms[1], [])
    cv_index = system.addForce(cv)
    
    # Parámetros de metadinámica well-tempered
    height = 1.2  # kJ/mol (altura inicial del gaussiano)
    sigma = 0.05  # nm (ancho del gaussiano)
    biasFactor = 10  # ΔT para well-tempered
    frequency = 1000  # Depositar cada 1 ps (500 steps)
    
    # Crear bias de metadinámica
    metad = mm.BiasVariable(cv_index, 2.0, 4.0, sigma, True)
    meta_force = mm.Metadynamics(
        system, [metad], 310*unit.kelvin, biasFactor,
        height*unit.kilojoules_per_mole, frequency
    )
    
    # Guardar HILLS cada 10 gaussianos
    meta_force.setReportInterval(10)
    meta_force.setReportFile(f'{output_prefix}_HILLS.txt')
    
    # Integrator
    integrator = mm.LangevinMiddleIntegrator(
        310*unit.kelvin,
        1.0/unit.picosecond,
        0.002*unit.picoseconds
    )
    
    # Simulation
    platform = mm.Platform.getPlatformByName('CPU')
    properties = {'Threads': '4'}
    simulation = app.Simulation(
        pdb.topology, system, integrator, platform, properties
    )
    simulation.context.setPositions(pdb.positions)
    
    # Reporters
    simulation.reporters.append(
        app.StateDataReporter(
            f'{output_prefix}_log.txt', 5000,
            step=True, time=True, potentialEnergy=True,
            temperature=True, speed=True
        )
    )
    simulation.reporters.append(
        app.DCDReporter(f'{output_prefix}_trajectory.dcd', 5000)
    )
    
    return simulation, meta_force

# Ejecución
if __name__ == '__main__':
    # Índices de Cα terminal (ajustar según tu sistema)
    CA_NTERM = 10  # Primer residuo visible
    CA_CTERM = 1250  # Último residuo antes del desorden
    
    sim, metad_force = setup_metadynamics(
        'wnk_pbs_equilibrated.pdb',
        (CA_NTERM, CA_CTERM),
        'metad_wnk'
    )
    
    # Correr 500 ns (250,000,000 steps)
    print("Starting well-tempered metadynamics (500 ns)...")
    sim.step(250000000)
    
    # Guardar estado final
    positions = sim.context.getState(getPositions=True).getPositions()
    with open('metad_wnk_final.pdb', 'w') as f:
        app.PDBFile.writeFile(sim.topology, positions, f)
    
    print("Metadynamics complete. Check metad_wnk_HILLS.txt for bias profile.")
```

### Comparación de Resultados

```python
# compare_umbrella_metad.py
import numpy as np
import matplotlib.pyplot as plt

def load_umbrella_pmf(mbar_file):
    """Cargar PMF de umbrella sampling (output de MBAR)."""
    data = np.loadtxt(mbar_file)
    return data[:, 0], data[:, 1], data[:, 2]  # cv, pmf, error

def load_metad_pmf(hills_file, cv_min=2.0, cv_max=4.0, bins=200):
    """Calcular PMF de metadinámica desde HILLS."""
    hills = np.loadtxt(hills_file)
    cv_centers = hills[:, 1]
    heights = hills[:, 2]
    sigmas = hills[:, 3]
    
    cv_values = np.linspace(cv_min, cv_max, bins)
    bias = np.zeros(bins)
    
    for i, cv in enumerate(cv_values):
        for j in range(len(hills)):
            bias[i] += heights[j] * np.exp(
                -(cv - cv_centers[j])**2 / (2 * sigmas[j]**2)
            )
    
    pmf = -bias
    pmf -= pmf.min()
    return cv_values, pmf

# Cargar ambos PMFs
cv_umb, pmf_umb, err_umb = load_umbrella_pmf('pmf_final.dat')
cv_met, pmf_met = load_metad_pmf('metad_wnk_HILLS.txt')

# Interpolar metadinámica al grid de umbrella
pmf_met_interp = np.interp(cv_umb, cv_met, pmf_met)

# Calcular RMSD
rmsd = np.sqrt(np.mean((pmf_umb - pmf_met_interp)**2))

# Plot
plt.figure(figsize=(10, 6))
plt.plot(cv_umb, pmf_umb, 'o-', label='Umbrella Sampling', linewidth=2)
plt.fill_between(cv_umb, pmf_umb-err_umb, pmf_umb+err_umb, alpha=0.3)
plt.plot(cv_met, pmf_met, 's-', label='Metadinámica', linewidth=2)
plt.xlabel('Distancia Cα-Cα (nm)', fontsize=12)
plt.ylabel('PMF (kJ/mol)', fontsize=12)
plt.title(f'Validación Cruzada: RMSD = {rmsd:.2f} kJ/mol', fontsize=14)
plt.legend(fontsize=11)
plt.grid(alpha=0.3)
plt.savefig('umbrella_metad_comparison.png', dpi=300, bbox_inches='tight')
plt.show()

print(f"\n{'='*50}")
print(f"VALIDACIÓN CRUZADA - UMBRELLA vs. METADINÁMICA")
print(f"{'='*50}")
print(f"RMSD entre perfiles:        {rmsd:.3f} kJ/mol")
print(f"Diferencia en barrera:       {abs(pmf_umb.max() - pmf_met.max()):.3f} kJ/mol")
print(f"Criterio de aceptación:      RMSD < 2.0 kJ/mol")
print(f"Resultado:                   {'✅ PASS' if rmsd < 2.0 else '❌ FAIL'}")
print(f"{'='*50}\n")
```

---

## 9. Referencias Clave

### Enhanced Sampling:
1. **Umbrella Sampling**: Torrie & Valleau (1977). *J. Comput. Phys.* 23, 187-199.
2. **MBAR**: Shirts & Chodera (2008). *J. Chem. Phys.* 129, 124105.
3. **Metadinámica**: Laio & Parrinello (2002). *PNAS* 99, 12562-12566.
4. **Well-Tempered Metad**: Barducci et al. (2008). *Phys. Rev. Lett.* 100, 020603.

### Alchemical Methods:
1. **FEP**: Zwanzig (1954). *J. Chem. Phys.* 22, 1420-1426.
2. **TI**: Kirkwood (1935). *J. Chem. Phys.* 3, 300-313.
3. **BAR**: Bennett (1976). *J. Comput. Phys.* 22, 245-268.
4. **Review**: Mobley & Klimovich (2012). *J. Comput. Aided Mol. Des.* 26, 93-109.

### Best Practices:
- Christ et al. (2010). *J. Chem. Inf. Model.* 50, 1787-1805. (Alchemical workflows)
- Barducci et al. (2011). *WIREs Comput. Mol. Sci.* 1, 826-843. (Metadynamics review)
- Hub et al. (2010). *J. Chem. Theory Comput.* 6, 3713-3720. (WHAM implementation)

---

## 10. Checklist Pre-Reunión con la Doctora

### Conceptos para Defender:

- [ ] **Entiendo las tres clases de métodos**:
  - Enhanced Sampling (umbrella, metadinámica) → PMF conformacional
  - Alchemical (FEP, TI, BAR) → ΔΔG de transformaciones químicas
  - Standard MD (drMD) → Exploración cualitativa

- [ ] **Puedo explicar por qué umbrella para WNK1**:
  - CV conocida (distancia Cα-Cα)
  - Barrera alta requiere enhanced sampling
  - Paralelización óptima en HPC (48 cores)

- [ ] **Sé cuándo usar métodos alquímicos**:
  - Mutaciones puntuales (ej. S1261A)
  - Afinidad de unión ligando
  - Energías de solvatación
  - **NO para barreras conformacionales**

- [ ] **Puedo explicar cómo metadinámica calcula ΔG**:
  - Bias gaussiano acumulativo
  - Teorema: $V_{bias} \rightarrow -PMF$
  - Extracción práctica desde archivo HILLS

- [ ] **Entiendo que metadinámica ≠ "mejor que umbrella"**:
  - Son métodos alternativos (misma pregunta científica)
  - Ambos calculan PMF
  - Umbrella: mejor para paralelización
  - Metadinámica: mejor para CVs desconocidas

### Frases Clave para Usar:

✅ **CORRECTO**:
- "Necesitamos enhanced sampling para cuantificar la barrera conformacional"
- "Umbrella sampling calcula el PMF mediante ventanas con bias armónico"
- "Metadinámica es una alternativa que usa bias gaussiano acumulativo"
- "Métodos alquímicos son para transformaciones químicas, no conformacionales"
- "drMD es excelente para exploración, pero no calcula ΔG"

❌ **EVITAR**:
- "Métodos alquímicos en metadinámica" → Mezcla conceptos incorrectamente
- "Metadinámica es mejor que umbrella" → No es mejor, es diferente
- "drMD no sirve" → SÍ sirve, pero para propósito diferente

---

## Conclusión

**Mensaje Principal**:  
Umbrella sampling y metadinámica son **ambos métodos de enhanced sampling** que calculan el mismo tipo de ΔG (PMF conformacional). Los **métodos alquímicos** (FEP/TI/BAR) son una **categoría completamente diferente** para transformaciones químicas. Para WNK1, necesitas enhanced sampling (umbrella), NO métodos alquímicos.

**Argumento Listo**:  
"Seleccioné umbrella sampling porque cuantifica la barrera conformacional del C-terminal de WNK1 mediante enhanced sampling. drMD es excelente para exploración pero no calcula ΔG. Metadinámica sería una alternativa válida (mismo resultado), pero umbrella se paraleliza mejor en nuestro HPC. Métodos alquímicos son para preguntas diferentes (mutaciones, unión de ligandos)."

---

**Éxito en tu reunión! 🚀**
