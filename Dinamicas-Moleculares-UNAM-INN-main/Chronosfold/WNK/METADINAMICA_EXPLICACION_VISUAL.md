# ¿Qué es Metadinámica? Explicación Visual

**Para**: Estudiante de 4to semestre de biología  
**Nivel**: Conceptual con matemáticas accesibles  
**Propósito**: Entender cómo metadinámica calcula PMF (no es MD estándar)

---

## 🚨 ACLARACIÓN CRÍTICA

### Metadinámica ≠ MD Estándar

**MD Estándar**:
```
Sistema → Explora según Hamiltoniano natural → Observas conformaciones
NO HAY BIAS → Puede quedarse atrapado en mínimos locales
```

**Metadinámica**:
```
Sistema → Explora con BIAS que se va sumando → "Rellena" valles energéticos
HAY BIAS ACUMULATIVO → Fuerza exploración de TODO el espacio
```

**Analogía**: 
- **MD estándar** = Dejar una pelota rodar en un valle → se queda abajo (mínimo local)
- **Metadinámica** = Ir llenando el valle con arena cada vez que la pelota visita → eventualmente la pelota DEBE salir y explorar otros valles

---

## 1. ¿Cómo Funciona Metadinámica? (Versión Visual)

### Paso a Paso

Imagina el PMF como un paisaje de colinas y valles:

```
Energía (PMF)
    ^
    |     __/\__                    __/\__
    |    /      \                  /      \
    |   /        \________________/        \
    |  /          Valle A  |  Barrera  | Valle B
    |_____________________________________________> CV (distancia)
       2.0 nm              3.0 nm            4.0 nm
```

**Problema de MD estándar**:
- Sistema empieza en Valle A (2.0 nm)
- Barrera es alta (~25 kJ/mol)
- Probabilidad de cruzar: $P \propto e^{-\Delta G/k_B T} = e^{-25/(0.008314 \times 310)} \approx 10^{-4}$
- **Tardaría microsegundos** en ver la transición espontáneamente

### Algoritmo de Metadinámica

**1. Sistema empieza explorando Valle A**:
```
t = 0 ps: Sistema en CV = 2.0 nm (Valle A)
PMF real: -20 kJ/mol (mínimo local)
Bias acumulado: 0 kJ/mol
```

**2. Cada 1 ps, depositas un "gaussiano" (pequeña colina)** en la posición actual:
```
t = 1 ps: Sistema sigue en CV ≈ 2.0 nm
Depositas: Gaussiano con altura h = 1.2 kJ/mol, ancho σ = 0.05 nm
Bias acumulado: +1.2 kJ/mol en CV = 2.0 nm
```

**3. Gaussianos se van sumando** (bias acumulativo):
```
t = 100 ps: Sistema visitó CV = 2.0-2.1 nm muchas veces
Depositas: 100 gaussianos en esa región
Bias acumulado: ≈ +20 kJ/mol en Valle A
```

**Efecto**: Valle A se "rellena" → **PMF efectivo** = PMF real - Bias:
```
PMF_efectivo = PMF_real + V_bias
             = -20 kJ/mol + 20 kJ/mol
             = 0 kJ/mol (¡aplanado!)
```

**4. Sistema DEBE salir del valle** porque ya no hay mínimo:
```
t = 500 ps: Valle A está "lleno"
Sistema explora hacia CV = 3.0 nm (barrera)
Depositas gaussianos en la barrera → la "aplanas"
```

**5. Eventualmente, TODO el espacio se rellena**:
```
t → ∞ (convergencia):
V_bias(CV) ≈ -PMF_real(CV) + constante
```

**TEOREMA FUNDAMENTAL**: En convergencia, el bias acumulado es **el negativo del PMF real**.

---

## 2. Ecuaciones (Nivel Accesible)

### Bias Gaussiano

Cada gaussiano depositado tiene forma:

$$G(s, s_0) = h \cdot \exp\left(-\frac{(s - s_0)^2}{2\sigma^2}\right)$$

Donde:
- $s$ = Collective variable (ej. distancia Cα-Cα)
- $s_0$ = Posición donde se depositó el gaussiano
- $h$ = Altura del gaussiano (ej. 1.2 kJ/mol)
- $\sigma$ = Ancho del gaussiano (ej. 0.05 nm)

**Visualización**:
```
Altura
    ^
 h  |     ___
    |   /     \
    | /         \___
    |_________________> CV
        s₀  (centro)
    |<- σ ->| (ancho)
```

### Bias Acumulativo

En cualquier tiempo $t$, el bias total es la **suma de TODOS los gaussianos depositados**:

$$V_{bias}(s, t) = \sum_{i=1}^{N(t)} h_i \cdot \exp\left(-\frac{(s - s_i)^2}{2\sigma_i^2}\right)$$

Donde:
- $N(t)$ = Número de gaussianos depositados hasta tiempo $t$
- $s_i$ = Posición del i-ésimo gaussiano

**Ejemplo numérico** (3 gaussianos):
```python
# Gaussiano 1: depositado en CV = 2.0 nm
G1 = 1.2 * exp(-(s - 2.0)^2 / (2*0.05^2))

# Gaussiano 2: depositado en CV = 2.05 nm
G2 = 1.2 * exp(-(s - 2.05)^2 / (2*0.05^2))

# Gaussiano 3: depositado en CV = 2.1 nm
G3 = 1.2 * exp(-(s - 2.1)^2 / (2*0.05^2))

# Bias total
V_bias = G1 + G2 + G3
```

### Teorema de Metadinámica

En el límite de convergencia ($t \rightarrow \infty$):

$$V_{bias}(s) \rightarrow -F(s) + C$$

Donde:
- $F(s)$ = Free energy (PMF real)
- $C$ = Constante arbitraria

**Por lo tanto**:

$$\boxed{F(s) = -V_{bias}(s) + C}$$

**Conclusión**: El PMF real es el **negativo del bias acumulado** (más una constante que podemos normalizar a 0).

---

## 3. Analogía del "Rellenado de Valle"

Imagina que el PMF es un valle en un jardín:

### MD Estándar
```
🏀 Pelota (sistema) en el fondo del valle
├─ Gravedad natural (PMF real)
├─ Pelota rueda, pero se queda atrapada en el fondo
└─ Para salir, necesita fluctuación térmica enorme (raro)
```

### Metadinámica
```
🏀 Pelota en el fondo del valle
├─ Cada segundo que pasa ahí, tiras una palada de arena 🏖️
├─ Arena se acumula → Valle se "rellena"
├─ Después de muchas paladas, el fondo está alto
└─ Pelota DEBE rodar hacia otro lado (exploración forzada)

Al final:
- Cantidad de arena en cada punto = -PMF(punto)
- Mides cuánta arena necesitaste → ¡Conoces el PMF original!
```

**Clave**: No esperas a que el evento ocurra espontáneamente (MD estándar), lo **fuerzas** rellenando donde ya visitaste.

---

## 4. Diferencia con Umbrella Sampling

### Umbrella Sampling
```
Estrategia: "Divide y conquista"
├─ Creas 20 ventanas con bias ESTÁTICO: k(s - s₀)²
├─ Cada ventana explora LOCAL (ej. 2.0-2.1 nm)
├─ Al final, combinas todas con MBAR/WHAM
└─ Bias NO cambia durante la simulación

Ventaja: Paralelizable (20 ventanas simultáneas)
Desventaja: Necesitas conocer el rango de CV de antemano
```

### Metadinámica
```
Estrategia: "Exploración adaptativa"
├─ 1 trayectoria con bias DINÁMICO que se acumula
├─ Sistema decide solo qué explorar (adaptativo)
├─ No necesitas predefinir ventanas
└─ Bias cambia cada 1 ps

Ventaja: Descubre transiciones inesperadas
Desventaja: Secuencial (no paraleliza fácilmente)
```

**Analogía**:
- **Umbrella**: 20 personas exploran 20 segmentos de montaña simultáneamente → rápido si tienes gente
- **Metadinámica**: 1 explorador inteligente que va marcando donde ya estuvo → lento pero adaptativo

---

## 5. ¿Cómo se Extrae el PMF?

### Proceso Práctico

Al final de la simulación de metadinámica (ej. 500 ns), tienes:

**Archivo HILLS** (registro de gaussianos):
```
# time(ps)   CV(nm)    height(kJ/mol)   sigma(nm)
0.0          2.05      1.200            0.050
1.0          2.06      1.195            0.050
2.0          2.07      1.190            0.050
...
500000.0     3.95      0.950            0.050
```

**Cálculo del PMF**:

```python
import numpy as np

def extract_pmf_from_hills(hills_file, cv_min=2.0, cv_max=4.0, bins=200):
    """
    Extrae PMF desde archivo HILLS de metadinámica.
    """
    # 1. Leer todos los gaussianos
    data = np.loadtxt(hills_file)
    times = data[:, 0]
    cv_centers = data[:, 1]  # Donde se depositaron
    heights = data[:, 2]      # Altura de cada gaussiano
    sigmas = data[:, 3]       # Ancho de cada gaussiano
    
    # 2. Crear grid de CV
    cv_grid = np.linspace(cv_min, cv_max, bins)
    bias = np.zeros(bins)
    
    # 3. SUMAR todos los gaussianos en cada punto del grid
    for i, cv in enumerate(cv_grid):
        for j in range(len(data)):
            # Contribución del gaussiano j en el punto cv
            bias[i] += heights[j] * np.exp(
                -(cv - cv_centers[j])**2 / (2 * sigmas[j]**2)
            )
    
    # 4. PMF = -Bias (teorema de metadinámica)
    pmf = -bias
    
    # 5. Normalizar: mínimo = 0 kJ/mol
    pmf -= pmf.min()
    
    return cv_grid, pmf

# Uso
cv, pmf = extract_pmf_from_hills('HILLS')

# Barrera de activación
delta_g = pmf.max() - pmf.min()
print(f"ΔG = {delta_g:.2f} kJ/mol")
```

**Salida típica**:
```
ΔG = 24.7 ± 1.2 kJ/mol
```

**Interpretación**:
- Costo energético de pasar de Valle A (2.0 nm) a Valle B (4.0 nm)
- **Mismo resultado que umbrella sampling** (si ambos convergen)

---

## 6. Well-Tempered Metadynamics (Versión Moderna)

### Problema de Metadinámica Clásica

Si depositas gaussianos de altura constante ($h$):
- Eventualmente **sobre-rellenas** → bias excede -PMF
- Sistema oscila sin converger
- Difícil saber cuándo parar

### Solución: Well-Tempered

La altura del gaussiano **disminuye** en regiones ya visitadas:

$$h(s, t) = h_0 \cdot \exp\left(-\frac{V_{bias}(s, t)}{k_B \Delta T}\right)$$

Donde:
- $h_0$ = Altura inicial (ej. 1.2 kJ/mol)
- $\Delta T$ = "Bias factor" (típicamente 10-15 para proteínas)
- $k_B$ = Constante de Boltzmann

**Efecto**:
```
Primera visita a CV = 2.0 nm:
├─ V_bias = 0 → h = 1.2 kJ/mol (altura completa)

Décima visita a CV = 2.0 nm:
├─ V_bias = 12 kJ/mol → h = 1.2 * exp(-12/(0.008314*310*10)) ≈ 0.6 kJ/mol

Centésima visita:
├─ V_bias = 20 kJ/mol → h ≈ 0.1 kJ/mol (muy pequeño)
```

**Ventajas**:
1. **Converge suavemente** → bias se estabiliza en -PMF
2. **Auto-regulación** → no necesitas adivinar cuándo parar
3. **Más preciso** → evita sobre-muestreo

**Es el estándar actual** (papers desde 2008).

---

## 7. Comparación Visual: Umbrella vs Metadinámica

### Umbrella Sampling

```
CV Space: [2.0 ─────────────────── 4.0 nm]

Ventana 1: [2.0-2.1]  ← Bias estático: k(s-2.0)²
Ventana 2: [2.1-2.2]  ← Bias estático: k(s-2.1)²
Ventana 3: [2.2-2.3]  ← Bias estático: k(s-2.2)²
...
Ventana 20: [3.9-4.0] ← Bias estático: k(s-4.0)²

├─ Cada ventana corre 100 ns (paralelo)
├─ 20 simulaciones simultáneas
└─ Al final: MBAR combina histogramas → PMF

Tiempo de reloj: ~10 días (48 cores)
Tiempo de simulación: 20 × 100 ns = 2 μs
```

### Metadinámica

```
CV Space: [2.0 ─────────────────── 4.0 nm]

Trayectoria única:
├─ t = 0-100 ns:   Explora 2.0-2.5 nm (rellena Valle A)
├─ t = 100-300 ns: Cruza barrera (rellena 2.5-3.5 nm)
├─ t = 300-500 ns: Explora 3.5-4.0 nm (rellena Valle B)
└─ Bias dinámico: ΣG(s, s_i) se actualiza cada 1 ps

├─ 1 simulación secuencial
└─ Al final: -V_bias → PMF

Tiempo de reloj: ~30 días (4-8 cores) o ~1 día (GPU)
Tiempo de simulación: 500 ns - 1 μs
```

---

## 8. ¿Por Qué Ambos Calculan el Mismo PMF?

### Fundamento Termodinámico

El PMF (Potential of Mean Force) es una **propiedad termodinámica**:

$$F(s) = -k_B T \ln P(s) + C$$

Donde:
- $P(s)$ = Probabilidad de observar CV = $s$ en equilibrio
- Es una propiedad del sistema, **NO del método**

**Analogía**: 
- Preguntar "¿Cuál es la temperatura de esta habitación?"
- Da igual si usas termómetro digital o de mercurio
- **La temperatura es la misma**, solo el método de medición cambia

### ¿Por Qué Funcionan?

**Umbrella**:
1. Bias armónico fuerza muestreo LOCAL en cada ventana
2. MBAR "descuenta" el bias → recupera P(s) real
3. Calcula F(s) = -kT ln P(s)

**Metadinámica**:
1. Bias gaussiano acumulativo "rellena" el PMF
2. En convergencia: V_bias = -F(s)
3. F(s) = -V_bias (inversión directa)

**Ambos recuperan la misma F(s)** porque es propiedad del sistema.

---

## 9. Ejemplo Concreto: WNK1 C-Terminal

### Setup Idéntico

**Sistema**:
- WNK1 kinase domain (residuos 220-1280)
- ~40,000 átomos (proteína + agua + iones)
- Forcefield: AMBER14 + TIP3P
- Condiciones: 310 K, 1 bar, PBS buffer

**Collective Variable**:
- Distancia Cα(N-term) - Cα(C-term)
- Rango: 2.0 - 4.0 nm
- Mide "apertura" del C-terminal

### Método 1: Umbrella Sampling

```bash
# 20 ventanas × 100 ns = 2 μs
python generate_umbrella_windows.py  # Crea configs
sbatch submit_umbrella_hpc_48cores.sh  # Corre en paralelo (10 días)
python analyze_umbrella_mbar.py       # Extrae PMF
```

**Resultado esperado**:
```
CV (nm)    PMF (kJ/mol)   Error
2.0        0.0            0.2
2.5        8.3            0.5
3.0        24.7           1.1    ← Barrera (TS)
3.5        15.2           0.7
4.0        2.1            0.3
```

### Método 2: Metadinámica (GPU)

```bash
# 1 trayectoria × 500 ns
python run_metadynamics_gpu.py  # Corre en GPU (1 día)
python extract_pmf_from_hills.py  # Extrae PMF desde HILLS
```

**Resultado esperado**:
```
CV (nm)    PMF (kJ/mol)
2.0        0.0
2.5        8.1
3.0        24.9    ← Barrera (casi idéntica)
3.5        15.4
4.0        2.0
```

### Validación Cruzada

```python
# RMSD entre ambos perfiles
rmsd = sqrt(mean((PMF_umbrella - PMF_metad)^2))
# Esperado: RMSD < 2 kJ/mol

# Diferencia en barrera
delta_barrier = |24.7 - 24.9| = 0.2 kJ/mol
# Excelente acuerdo!
```

**Conclusión**: Ambos métodos recuperan **el mismo PMF** (dentro de error estadístico).

---

## 10. ¿Entonces Cuál Usar?

### Criterios de Decisión

| Criterio | Umbrella | Metadinámica | Ganador |
|----------|----------|--------------|---------|
| **CV conocida** | ✅ Necesario | ⚠️ Ayuda pero no esencial | Umbrella |
| **Paralelización** | ✅ Excelente (20 cores) | ❌ Secuencial | Umbrella |
| **Recursos: HPC 48 cores** | ✅ Óptimo | ⚠️ Sub-usado | Umbrella |
| **Recursos: 1 GPU potente** | ⚠️ No aprovecha | ✅ Excelente | Metadinámica |
| **CV desconocida** | ❌ Problema | ✅ Adaptativo | Metadinámica |
| **Múltiples CVs** | ⚠️ Complejo (2D MBAR) | ✅ Natural | Metadinámica |
| **Tiempo de reloj (CPU)** | ✅ 10 días | ❌ 30-60 días | Umbrella |
| **Convergencia verificable** | ✅ Bootstrap MBAR | ⚠️ Recrossing analysis | Umbrella |
| **Validación cruzada** | ✅ Ambos (complementarios) | ✅ Ambos | Empate |

### Decisión para WNK1

**Tu caso**:
- ✅ CV conocida (distancia Cα-Cα)
- ✅ HPC 48 cores disponible
- ✅ GPU también disponible
- ⏰ Deadline de semestre (tiempo limitado)

**Estrategia óptima**:
1. **Primero**: Umbrella en CPU (48 cores) → 10 días → PMF principal
2. **Después**: Metadinámica en GPU (paralelo) → 1 día → Validación
3. **Comparar**: RMSD < 2 kJ/mol → Confianza en resultado
4. **Defensa**: "Usé dos métodos independientes que convergen al mismo PMF"

---

## 11. Implementación: Código Listo para WNK1

### `run_metadynamics_gpu.py`

```python
#!/usr/bin/env python3
"""
Metadinámica Well-Tempered para WNK1 C-terminal.
Optimizado para GPU (CUDA).
"""

import openmm as mm
from openmm import app, unit
from openmm.app import PDBFile, Modeller, ForceField
import sys

def setup_metadynamics_wnk1(
    pdb_file='wnk_pbs_equilibrated.pdb',
    ca_nterm_index=10,    # Ajustar según tu sistema
    ca_cterm_index=1250,  # Ajustar según tu sistema
    output_prefix='metad_wnk1',
    production_ns=500
):
    """
    Configura well-tempered metadynamics para WNK1.
    
    Parámetros:
    -----------
    pdb_file : str
        Sistema equilibrado (post NPT)
    ca_nterm_index : int
        Índice del átomo Cα del N-terminal
    ca_cterm_index : int
        Índice del átomo Cα del C-terminal
    output_prefix : str
        Prefijo para archivos de salida
    production_ns : int
        Duración de producción en nanosegundos
    """
    
    print("="*60)
    print("CONFIGURANDO WELL-TEMPERED METADYNAMICS")
    print("="*60)
    
    # 1. Cargar sistema
    print(f"\n1. Cargando {pdb_file}...")
    pdb = PDBFile(pdb_file)
    
    # 2. Forcefield (PBS buffer integrado)
    print("2. Creando sistema con AMBER14...")
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
    
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005
    )
    
    # 3. Definir Collective Variable (distancia Cα-Cα)
    print(f"3. Definiendo CV: distancia Cα({ca_nterm_index}) - Cα({ca_cterm_index})...")
    
    cv_force = mm.CustomBondForce('r')
    cv_force.addBond(ca_nterm_index, ca_cterm_index, [])
    cv_index = system.addForce(cv_force)
    
    # 4. Parámetros de Well-Tempered Metadynamics
    print("4. Configurando parámetros de metadinámica...")
    
    # Parámetros recomendados para proteínas
    height = 1.2 * unit.kilojoules_per_mole  # Altura inicial del gaussiano
    sigma = 0.05 * unit.nanometers           # Ancho del gaussiano (~1 Å en CV)
    biasFactor = 10                          # ΔT para well-tempered (10-15 típico)
    frequency = 500                          # Depositar cada 1 ps (500 steps × 2 fs)
    
    print(f"   - Altura inicial: {height}")
    print(f"   - Ancho gaussiano (σ): {sigma}")
    print(f"   - Bias factor: {biasFactor}")
    print(f"   - Frecuencia: cada {frequency} steps (1 ps)")
    
    # Rango de CV (2.0 - 4.0 nm basado en análisis previo)
    cv_min = 2.0 * unit.nanometers
    cv_max = 4.0 * unit.nanometers
    
    # 5. Crear BiasVariable y Metadynamics Force
    print("5. Creando Metadynamics Force...")
    
    metad_var = mm.BiasVariable(
        cv_index,
        cv_min,
        cv_max,
        sigma,
        True  # Periodic = False (distancia no es periódica)
    )
    
    temperature = 310 * unit.kelvin
    
    metad_force = mm.Metadynamics(
        system,
        [metad_var],
        temperature,
        biasFactor,
        height,
        frequency
    )
    
    # Guardar HILLS cada 10 gaussianos (~10 ps)
    metad_force.setReportInterval(10)
    hills_file = f'{output_prefix}_HILLS.txt'
    metad_force.setReportFile(hills_file)
    print(f"   - Archivo HILLS: {hills_file}")
    
    # 6. Integrator (Langevin para NPT)
    print("6. Configurando Langevin integrator (NPT)...")
    
    integrator = mm.LangevinMiddleIntegrator(
        temperature,
        1.0 / unit.picosecond,     # Fricción
        0.002 * unit.picoseconds   # Timestep = 2 fs
    )
    
    # Barostato para presión constante (1 bar)
    system.addForce(mm.MonteCarloBarostat(
        1.0 * unit.bar,
        temperature,
        25  # Actualizar cada 25 steps
    ))
    
    # 7. Platform (GPU CUDA)
    print("7. Seleccionando platform GPU...")
    
    try:
        platform = mm.Platform.getPlatformByName('CUDA')
        properties = {
            'DeviceIndex': '0',        # GPU 0 (ajustar si multi-GPU)
            'Precision': 'mixed'       # mixed = balance velocidad/precisión
        }
        print("   ✅ Platform: CUDA (GPU)")
    except Exception as e:
        print(f"   ⚠️  CUDA no disponible: {e}")
        print("   ℹ️  Usando CPU (será más lento)")
        platform = mm.Platform.getPlatformByName('CPU')
        properties = {'Threads': '8'}
    
    # 8. Crear Simulation
    print("8. Creando objeto Simulation...")
    
    simulation = app.Simulation(
        pdb.topology,
        system,
        integrator,
        platform,
        properties
    )
    
    simulation.context.setPositions(pdb.positions)
    
    # Si hay velocidades guardadas, cargarlas
    try:
        with open(pdb_file.replace('.pdb', '_velocities.xml')) as f:
            simulation.context.setVelocities(mm.XmlSerializer.deserialize(f.read()))
        print("   ✅ Velocidades cargadas desde XML")
    except FileNotFoundError:
        simulation.context.setVelocitiesToTemperature(temperature)
        print("   ℹ️  Velocidades inicializadas a 310 K")
    
    # 9. Reporters
    print("9. Configurando reporters...")
    
    # Log de energías y velocidad (cada 10 ps)
    simulation.reporters.append(
        app.StateDataReporter(
            f'{output_prefix}_log.txt',
            5000,  # Cada 10 ps
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            volume=True,
            speed=True,
            remainingTime=True,
            totalSteps=production_ns * 500000  # ns → steps
        )
    )
    
    # Trayectoria DCD (cada 100 ps para ahorrar espacio)
    simulation.reporters.append(
        app.DCDReporter(
            f'{output_prefix}_trajectory.dcd',
            50000  # Cada 100 ps
        )
    )
    
    # Checkpoint cada 10 ns (para reiniciar si falla)
    simulation.reporters.append(
        app.CheckpointReporter(
            f'{output_prefix}_checkpoint.chk',
            5000000  # Cada 10 ns
        )
    )
    
    print(f"   - Log: {output_prefix}_log.txt")
    print(f"   - Trayectoria: {output_prefix}_trajectory.dcd")
    print(f"   - Checkpoint: {output_prefix}_checkpoint.chk")
    
    return simulation, metad_force, hills_file


def run_metadynamics(simulation, production_ns, output_prefix):
    """
    Ejecuta la producción de metadinámica.
    """
    
    print("\n" + "="*60)
    print(f"INICIANDO PRODUCCIÓN ({production_ns} ns)")
    print("="*60)
    
    # Minimización rápida (por si acaso)
    print("\nMinimización inicial...")
    simulation.minimizeEnergy(maxIterations=1000)
    
    # Calcular steps totales (500,000 steps/ns con timestep 2 fs)
    total_steps = production_ns * 500000
    
    print(f"\nEjecutando {production_ns} ns ({total_steps:,} steps)...")
    print("Presiona Ctrl+C para detener de forma segura\n")
    
    try:
        simulation.step(total_steps)
    except KeyboardInterrupt:
        print("\n⚠️  Simulación interrumpida por usuario")
    
    # Guardar estado final
    print("\nGuardando estado final...")
    
    positions = simulation.context.getState(getPositions=True).getPositions()
    velocities = simulation.context.getState(getVelocities=True).getVelocities()
    
    with open(f'{output_prefix}_final.pdb', 'w') as f:
        app.PDBFile.writeFile(
            simulation.topology,
            positions,
            f
        )
    
    with open(f'{output_prefix}_velocities.xml', 'w') as f:
        f.write(mm.XmlSerializer.serialize(velocities))
    
    print(f"✅ Estado final guardado: {output_prefix}_final.pdb")


if __name__ == '__main__':
    
    # Configuración
    CA_NTERM = 10    # ⚠️ AJUSTAR según tu sistema WNK1
    CA_CTERM = 1250  # ⚠️ AJUSTAR según tu sistema WNK1
    PRODUCTION_NS = 500  # 500 ns (ajustar según GPU)
    
    print("\n🧬 WELL-TEMPERED METADYNAMICS - WNK1 C-TERMINAL")
    print("📌 Collective Variable: Distancia Cα-Cα terminal")
    print(f"⏱️  Duración: {PRODUCTION_NS} ns\n")
    
    # Setup
    sim, metad, hills = setup_metadynamics_wnk1(
        pdb_file='wnk_pbs_equilibrated.pdb',
        ca_nterm_index=CA_NTERM,
        ca_cterm_index=CA_CTERM,
        output_prefix='metad_wnk1',
        production_ns=PRODUCTION_NS
    )
    
    # Run
    run_metadynamics(sim, PRODUCTION_NS, 'metad_wnk1')
    
    print("\n" + "="*60)
    print("✅ METADINÁMICA COMPLETADA")
    print("="*60)
    print(f"\nArchivos generados:")
    print(f"  1. metad_wnk1_HILLS.txt      ← Gaussianos depositados")
    print(f"  2. metad_wnk1_log.txt        ← Energías y velocidad")
    print(f"  3. metad_wnk1_trajectory.dcd ← Trayectoria completa")
    print(f"  4. metad_wnk1_final.pdb      ← Estado final")
    print(f"\nPróximo paso:")
    print(f"  python extract_pmf_from_hills.py metad_wnk1_HILLS.txt")
    print()
```

---

## 12. Resumen Ejecutivo

### ¿Qué es Metadinámica?

**NO es MD estándar** - es una técnica de **enhanced sampling** que:

1. **Deposita gaussianos** (pequeñas colinas) cada 1 ps donde el sistema está
2. **Rellena valles energéticos** → fuerza exploración de TODO el espacio de CV
3. **En convergencia**: Bias acumulado = -PMF real
4. **Extrae PMF**: Suma todos los gaussianos y invierte el signo

### ¿Por Qué Calcula el Mismo PMF que Umbrella?

- Ambos son métodos que **recuperan propiedades termodinámicas**
- PMF = -kT ln P(CV) es propiedad del sistema, no del método
- **Umbrella**: Bias estático + MBAR → PMF
- **Metadinámica**: Bias dinámico acumulativo → -Bias = PMF

### Tu Argumento Refinado

**Versión original** (buena intuición):
> "no nos conviene metadinámica porque nos quedaríamos esperando a que pase el evento"

**Versión técnica** (para la doctora):
> "Elegí umbrella porque **paraleliza mejor** en nuestro HPC (20 ventanas × 48 cores = ~10 días). Metadinámica es secuencial (~30-60 días en CPU), pero con GPU sería ~1 día. **Ambos convergen al mismo PMF** (misma calidad termodinámica), la diferencia es eficiencia computacional. Podemos hacer ambos para validación cruzada."

---

## 13. Conclusión

### Tres Puntos Clave

1. **Metadinámica ≠ MD estándar**:
   - MD estándar: explora naturalmente (lento para barreras altas)
   - Metadinámica: aplica bias acumulativo que "rellena" valles (fuerza exploración completa)

2. **Cálculo de PMF**:
   - Deposita gaussianos donde visita
   - En convergencia: V_bias → -PMF
   - PMF = -V_bias (inversión directa)

3. **Umbrella vs Metadinámica**:
   - **NO es "mejor" o "peor"** - es diferencia de estrategia
   - Umbrella: paralelo (rápido con HPC)
   - Metadinámica: secuencial (rápido con GPU, adaptativo)
   - **Ambos calculan el mismo ΔG** (propiedad termodinámica)

### Analogía Final

**PMF = Mapa topográfico de montañas**

- **MD estándar** = Caminar aleatoriamente y esperar cruzar la montaña (muy lento)
- **Umbrella** = 20 equipos exploran 20 segmentos simultáneamente, luego combinas mapas
- **Metadinámica** = 1 explorador que va rellenando donde ya estuvo hasta que TODO está plano, mides cuánto rellenaste = altura original

**Todos miden la misma montaña, solo difieren en la estrategia de exploración.**

---

¿Quedó claro ahora? 🎯 Metadinámica **SÍ es una técnica sofisticada** (no MD estándar) que calcula PMF mediante bias acumulativo que "invierte" el perfil energético. ¡Es tan válida como umbrella, solo que diferente estrategia computacional!
