# Comparación: Umbrella Sampling vs MD Estándar (drMD)

## 📊 Resumen Ejecutivo

**TL;DR**: 
- **Umbrella Sampling**: 20 ventanas × 100 ns = **2 μs TOTAL** (pero en trayectorias separadas con bias)
- **drMD**: 1 simulación × 50 ns = **50 ns** (trayectoria continua libre)
- **NO son equivalentes** - diferentes propósitos y análisis

---

## 🔬 Comparación de Longitudes de Simulación

### Configuración Actual

| Método | Simulaciones | Tiempo/Sim | Tiempo Total | Tipo de Trayectoria |
|--------|--------------|------------|--------------|---------------------|
| **Umbrella Sampling** | 20 ventanas | 100 ns | **2000 ns (2 μs)** | Múltiples, con bias armónico |
| **drMD Production** | 1 (o 4 paralelas) | 50 ns | **50 ns (o 200 ns)** | Continua(s), sin bias |

### Detalles Técnicos

**Umbrella Sampling** (`submit_umbrella_hpc_48cores.sh`):
```bash
python run_umbrella_window.py \
    --window $WINDOW_ID \
    --steps 50000000 \      # 50M steps
    --dt 0.002 \            # 2 fs timestep
    --platform CPU
# 50,000,000 × 2 fs = 100 ns por ventana
# 20 ventanas × 100 ns = 2000 ns TOTAL
```

**drMD** (`drMD_wnk_config_hpc.yaml`):
```yaml
- name: "production"
  type: "md"
  steps: 25000000         # 25M steps
  timestep: 0.002         # 2 fs
  # 25,000,000 × 2 fs = 50 ns por simulación
  
parallelization:
  parallelCPU: 4          # Opcionalmente 4 en paralelo
  # 4 × 50 ns = 200 ns total (si se ejecutan 4)
```

---

## 🎯 Diferencias Fundamentales

### 1. **Estructura de Outputs**

#### Umbrella Sampling
```
umbrella_windows/
├── window_00/
│   ├── trajectory.dcd    # 100 ns con bias en ~2.0 nm
│   ├── cv_values.dat     # Serie temporal de distancia C-terminal
│   └── final.pdb
├── window_01/
│   ├── trajectory.dcd    # 100 ns con bias en ~2.1 nm
│   └── ...
└── window_19/
    ├── trajectory.dcd    # 100 ns con bias en ~4.0 nm
    └── ...
```

✅ **20 archivos .dcd SEPARADOS**  
✅ **Cada ventana explora región restringida** (bias armónico mantiene distancia cerca del target)  
❌ **NO se pueden concatenar** (diferentes condiciones de bias)  
✅ **Total: 2 μs de tiempo de simulación**

#### drMD / MD Estándar
```
drMD_simulations/
└── production/
    ├── production.dcd    # 50 ns sin bias (exploración libre)
    ├── state.xml
    └── checkpoint.chk
```

✅ **1 archivo .dcd CONTINUO**  
✅ **Exploración libre** (sin restricciones artificiales)  
✅ **Puede extenderse** (cargar checkpoint y continuar)  
❌ **Menos tiempo total** (50 ns vs 2000 ns)

---

### 2. **Naturaleza de las Trayectorias**

| Aspecto | Umbrella Sampling | MD Estándar (drMD) |
|---------|-------------------|-------------------|
| **Continuidad** | 20 trayectorias independientes | 1 trayectoria continua |
| **Bias** | ✅ Sí (armónico en distancia C-term) | ❌ No (evolución natural) |
| **Exploración** | Restringida a región de ventana | Libre en todo el espacio |
| **Transiciones** | ❌ No ve transiciones entre estados | ✅ Puede observar transiciones |
| **Tiempo efectivo** | Alto (2 μs total) | Moderado (50 ns) |
| **Muestreo** | Exhaustivo en 1D (distancia) | Muestreo libre multidimensional |

---

### 3. **Métodos de Análisis**

#### Umbrella Sampling → MBAR/WHAM

**Propósito**: Calcular **Potential of Mean Force (PMF)** - energía libre vs distancia

**Pipeline**:
```python
# 1. Cada ventana genera histograma de distancia
for window in range(20):
    cv_data = np.loadtxt(f"window_{window:02d}/cv_values.dat")
    histogram, bins = np.histogram(cv_data, bins=50)

# 2. MBAR combina histogramas considerando bias
from pymbar import MBAR
mbar = MBAR(u_kn, N_k)
pmf, uncertainties = mbar.computePMF()

# 3. Output: PMF(distancia)
plot(distance, pmf)  # kJ/mol vs nm
```

**Outputs**:
- ✅ **PMF**: Perfil de energía libre 1D
- ✅ **Barreras energéticas**: Altura de transiciones
- ✅ **Estados estables**: Mínimos en PMF
- ❌ **NO da trayectorias realistas** (son biased)

#### MD Estándar → Análisis Estructural

**Propósito**: Explorar **dinámica natural** sin perturbaciones

**Pipeline**:
```python
import mdtraj as md

# 1. Cargar trayectoria continua
traj = md.load('production.dcd', top='equilibrated.pdb')

# 2. Análisis estándar
rmsd = md.rmsd(traj, traj[0])          # Estabilidad
rmsf = md.rmsf(traj)                   # Fluctuaciones por residuo
contacts = md.compute_contacts(traj)   # Mapa de contactos
clusters = md.cluster(traj)            # Clustering conformacional

# 3. Propiedades termodinámicas
rg = md.compute_rg(traj)               # Radio de giro
sasa = md.shrake_rupley(traj)          # Superficie accesible
```

**Outputs**:
- ✅ **Trayectorias realistas**: Sin bias artificial
- ✅ **Transiciones observadas**: Si ocurren en 50 ns
- ✅ **Clustering**: Conformaciones representativas
- ❌ **Puede no converger PMF** (tiempo insuficiente)

---

## ⚖️ ¿Se Pueden Analizar Igual?

### ❌ **NO** - Son complementarias, no equivalentes

| Pregunta | Umbrella Sampling | MD Estándar |
|----------|-------------------|-------------|
| **¿Calcula PMF?** | ✅ Sí (propósito principal) | ⚠️ Difícil (necesita mucho tiempo) |
| **¿Muestra trayectorias realistas?** | ❌ No (biased) | ✅ Sí (natural) |
| **¿Ve transiciones entre estados?** | ❌ No (restringido) | ✅ Sí (si ocurren) |
| **¿Requiere post-procesamiento especial?** | ✅ Sí (MBAR/WHAM) | ❌ No (análisis estándar) |
| **¿Cuánto tiempo necesita?** | Moderado (100 ns/ventana) | Mucho (μs para PMF) |

---

## 🔄 ¿Se Pueden Usar Réplicas?

### ✅ **SÍ** - Ambos métodos soportan réplicas

#### Umbrella Sampling con Réplicas

**Implementación**:
```bash
# Modificar SLURM array para 3 réplicas × 20 ventanas = 60 jobs
#SBATCH --array=0-59

# En el script
WINDOW_ID=$((SLURM_ARRAY_TASK_ID % 20))
REPLICA_ID=$((SLURM_ARRAY_TASK_ID / 20))

# Salvar en directorios separados
OUT_DIR="umbrella_windows/window_${WINDOW_ID}_replica_${REPLICA_ID}"
```

**Ventajas**:
- Mejor estadística en histogramas
- Detectar histéresis
- Reducir incertidumbre en PMF

**Costo**:
- 3 réplicas × 20 ventanas × 100 ns = **6 μs total**
- 3× más tiempo de cómputo

#### MD Estándar con Réplicas

**drMD automático**:
```yaml
# En drMD, simplemente crear múltiples PDBs
inputDir/
├── wnk1_replica_1.pdb
├── wnk1_replica_2.pdb
└── wnk1_replica_3.pdb

# drMD ejecuta automáticamente 3 simulaciones
```

**Ventajas**:
- Explorar múltiples trayectorias independientes
- Mejor muestreo de espacio conformacional
- Estadística robusta en clustering

**Costo**:
- 3 réplicas × 50 ns = **150 ns total**

---

## 🎓 ¿Cuándo Usar Cada Método?

### Usar **Umbrella Sampling** cuando:

✅ **Quieres calcular energía libre** (PMF) de una reacción/transición  
✅ **Conoces la coordenada de reacción** (ej. distancia C-terminal)  
✅ **Barreras energéticas altas** (>15-20 kJ/mol) - MD libre no cruza en tiempo razonable  
✅ **Necesitas cuantificar termodinámica** (ΔG, barreras)  
✅ **Comparar múltiples condiciones** (PBS vs estándar, mutantes)  

**Ejemplo**: 
> "¿Cuál es la energía libre de disociación del C-terminal de WNK1?"

### Usar **MD Estándar** cuando:

✅ **Exploración conformacional** sin hipótesis previa  
✅ **Ver dinámica natural** sin perturbaciones  
✅ **Caracterizar fluctuaciones** (RMSD, RMSF, contactos)  
✅ **Clustering de estados** conformacionales  
✅ **Estudiar efectos alostéricos** a largo plazo  

**Ejemplo**: 
> "¿Cómo se mueve WNK1 naturalmente en solución? ¿Qué regiones son más flexibles?"

---

## 🔬 Validaciones Cruzadas

### Estrategia Recomendada: Usar Ambos Métodos

```
1. MD Estándar (exploración inicial)
   ├── 3 réplicas × 50 ns = 150 ns
   ├── Identificar conformaciones relevantes
   └── Definir coordenada de reacción (CV)

2. Umbrella Sampling (cuantificación)
   ├── Usar CV identificado en MD
   ├── 20 ventanas × 100 ns = 2 μs
   └── Calcular PMF con MBAR

3. Validación
   ├── Comparar mínimos PMF con clusters MD
   ├── Verificar barreras energéticas
   └── PMF debe explicar transiciones vistas en MD
```

**Benchmark**: Si PMF predice barrera de 25 kJ/mol, MD debería:
- Mostrar transiciones raras (bajo Boltzmann: e^(-25/RT) ≈ 0.00005)
- Pasar ~20,000 ns para ver 1 transición

---

## 🛠️ drMD: ¿Solo un Wrapper?

### ❌ **NO** - drMD hace validaciones automáticas que nosotros NO hacemos

| Validación | drMD | Nuestro Pipeline Manual |
|------------|------|-------------------------|
| **pdbTriage** | ✅ Valida PDB antes de simular | ❌ Asumimos PDB válido |
| **FirstAid** | ✅ Auto-recupera crashes (10 reintentos) | ❌ Si crashea, se pierde job |
| **Protonation states** | ✅ Automático según pH | ✅ Manual con ProPKa |
| **Clustering** | ✅ Automático post-MD | ❌ Manual con MDTraj |
| **Methods generation** | ✅ Auto-genera sección para paper | ❌ Manual |
| **Error handling** | ✅ Logging detallado + recovery | ⚠️ Básico (stderr SLURM) |
| **QC reports** | ✅ Reportes automáticos | ❌ No generamos |

### Validaciones Específicas de drMD

#### 1. **pdbTriage** (Pre-simulación)
```python
# drMD verifica automáticamente:
- Missing atoms
- Non-standard residues
- Chain breaks
- Clash detection
- Protonation consistency
```

#### 2. **FirstAid** (Durante simulación)
```python
# Si la simulación crashea, drMD automáticamente:
- Reduce timestep (2 fs → 1 fs)
- Aumenta frecuencia de constraints
- Reduce temperatura temporalmente
- Reintenta hasta 10 veces
```

#### 3. **Automatic Clustering**
```yaml
analysis:
  clustering: true
  cluster_method: "kmeans"
  n_clusters: 10
  # Genera automáticamente:
  # - cluster_centers.pdb
  # - cluster_populations.csv
  # - cluster_rmsd_matrix.npy
```

### ¿Deberíamos Usar drMD para Umbrella?

**Actualmente**: Usamos OpenMM directo (más control)

**Ventajas de drMD**:
- ✅ FirstAid podría rescatar ventanas que crashean
- ✅ pdbTriage detectaría problemas en estructuras iniciales
- ✅ Clustering automático de cada ventana

**Desventajas de drMD**:
- ❌ No tiene soporte nativo para umbrella sampling
- ❌ Tendríamos que implementar bias manualmente
- ❌ Menos control sobre reporters customizados

**Recomendación**: 
- Mantener OpenMM para umbrella (más flexible)
- Usar drMD para MD estándar (aprovecha validaciones)

---

## 📈 Aumentar Tiempo de Simulación

### Si Necesitas Más Tiempo...

#### Umbrella Sampling

**Opción 1: Aumentar tiempo por ventana**
```bash
# En submit_umbrella_hpc_48cores.sh
python run_umbrella_window.py \
    --steps 100000000 \     # 100M → 200 ns por ventana
    --window $WINDOW_ID
# Total: 20 × 200 ns = 4 μs
```

**Opción 2: Agregar réplicas**
```bash
#SBATCH --array=0-59    # 3 réplicas × 20 ventanas
# Total: 60 × 100 ns = 6 μs
```

**Opción 3: Más ventanas (mejor resolución)**
```bash
#SBATCH --array=0-39    # 40 ventanas (spacing 0.05 nm)
# Total: 40 × 100 ns = 4 μs
```

#### MD Estándar (drMD)

**Opción 1: Extender producción**
```yaml
- name: "production"
  steps: 100000000       # 100M → 200 ns
```

**Opción 2: Múltiples stages**
```yaml
- name: "production_1"
  steps: 25000000        # 50 ns
- name: "production_2"
  steps: 25000000        # 50 ns (continúa desde production_1)
- name: "production_3"
  steps: 25000000        # 50 ns
# Total: 150 ns continua
```

**Opción 3: Réplicas independientes**
```bash
# Crear múltiples PDBs en inputDir
inputDir/
├── wnk1_1.pdb
├── wnk1_2.pdb
├── wnk1_3.pdb
└── wnk1_4.pdb
# drMD ejecuta 4 × 50 ns = 200 ns total
```

---

## 🧪 Ejemplo Práctico: Análisis Completo

### Workflow Completo Recomendado

```bash
# ═══════════════════════════════════════════════════════════
# FASE 1: EXPLORACIÓN INICIAL (MD Estándar)
# ═══════════════════════════════════════════════════════════

# 1.1. Ejecutar 3 réplicas con drMD (50 ns cada una)
cd Chronosfold/WNK
python run_drMD_pipeline.py --config drMD_wnk_config_hpc.yaml

# 1.2. Analizar trayectorias
python << EOF
import mdtraj as md
trajs = [md.load(f'drMD_simulations/replica_{i}/production.dcd') 
         for i in range(1, 4)]

# RMSD, RMSF, clustering
combined = md.join(trajs)
rmsd = md.rmsd(combined, combined[0])
clusters = md.cluster.KMeans(n_clusters=5).fit(combined)

# Identificar CV relevante (distancia C-terminal)
cterm_distance = md.compute_distances(combined, [[kinase_ca, cterm_ca]])
print(f"Rango observado: {cterm_distance.min():.2f} - {cterm_distance.max():.2f} nm")
EOF

# Output: "Rango observado: 2.1 - 3.8 nm"
# → Definir ventanas umbrella en este rango

# ═══════════════════════════════════════════════════════════
# FASE 2: CUANTIFICACIÓN (Umbrella Sampling)
# ═══════════════════════════════════════════════════════════

# 2.1. Generar 20 ventanas entre 2.0-4.0 nm
python generate_umbrella_windows.py

# 2.2. Ejecutar umbrella (2 μs total)
bash run_complete_pipeline.sh hpc

# 2.3. Calcular PMF
python analyze_umbrella_mbar.py
# Output: pmf_analysis/pmf_results.csv

# 2.4. Visualizar
python visualize_results.py --animation

# ═══════════════════════════════════════════════════════════
# FASE 3: VALIDACIÓN CRUZADA
# ═══════════════════════════════════════════════════════════

# 3.1. Comparar mínimos PMF con clusters MD
python << EOF
import pandas as pd
pmf = pd.read_csv('pmf_analysis/pmf_results.csv')
min_pmf_cv = pmf.loc[pmf['pmf'].idxmin(), 'cv']

# Verificar que clusters de MD están cerca de mínimo PMF
cluster_centers_cv = [compute_cv(center) for center in cluster_centers]
print(f"PMF mínimo en: {min_pmf_cv:.2f} nm")
print(f"Clusters MD en: {cluster_centers_cv}")
# Deberían estar correlacionados!
EOF

# 3.2. Verificar tiempo de transición predicho
barrier = pmf['pmf'].max() - pmf['pmf'].min()  # kJ/mol
tau = compute_transition_time(barrier, T=310)  # Kramers theory
print(f"Barrera: {barrier:.1f} kJ/mol")
print(f"Tiempo transición esperado: {tau:.1f} ns")

# Si tau ~ 1000 ns y MD solo es 50 ns → OK no vimos transiciones
# Si tau ~ 10 ns y MD es 50 ns → Debería haber visto transiciones
```

---

## 📚 Referencias y Lecturas Adicionales

### Umbrella Sampling
- **Teoría**: Torrie & Valleau (1977) - Original umbrella sampling paper
- **MBAR**: Shirts & Chodera (2008) - Modern reweighting method
- **Best Practices**: Hub et al. (2010) - Choosing CVs and windows

### drMD
- **GitHub**: https://github.com/wells-wood-research/drMD
- **Paper**: [Pendiente - check repository]

### Comparaciones
- **Umbrella vs Metadynamics**: Barducci et al. (2011)
- **Free Energy Methods Review**: Christ & Fox (2014)

---

## ✅ Checklist de Decisión

**¿Qué método usar?**

```
[ ] ¿Quieres calcular energía libre (PMF)?
    → SÍ: Umbrella Sampling
    → NO: ↓

[ ] ¿Conoces la coordenada de reacción (CV)?
    → NO: Primero explorar con MD estándar
    → SÍ: ↓

[ ] ¿Esperas barreras energéticas altas (>20 kJ/mol)?
    → SÍ: Umbrella Sampling (MD no las cruza)
    → NO: MD estándar puede ser suficiente

[ ] ¿Quieres ver dinámica natural sin bias?
    → SÍ: MD estándar
    → NO: ↓

[ ] ¿Tienes tiempo de cómputo limitado?
    → SÍ: Umbrella (más eficiente para PMF)
    → NO: Hacer ambos (exploración + cuantificación)
```

---

## 🎯 Conclusiones Clave

1. **Longitudes de Simulación**:
   - Umbrella: **2 μs total** (20 × 100 ns), pero trayectorias cortas con bias
   - drMD: **50 ns** continua libre (o 200 ns con 4 réplicas)

2. **NO Son Equivalentes**:
   - Umbrella → PMF (energía libre 1D)
   - MD estándar → Dinámica natural (exploración multidimensional)

3. **Se Pueden Analizar Diferente**:
   - Umbrella necesita MBAR/WHAM
   - MD usa análisis estándar (RMSD, clustering, etc.)

4. **Réplicas Soportadas**:
   - Ambos métodos permiten réplicas
   - Umbrella: Mejor estadística en histogramas
   - MD: Mejor muestreo conformacional

5. **drMD NO es Solo un Wrapper**:
   - pdbTriage, FirstAid, clustering automático
   - Validaciones que nuestro pipeline manual NO tiene

6. **Uso Complementario Recomendado**:
   ```
   MD estándar → Explorar → Identificar CV
            ↓
   Umbrella Sampling → Cuantificar → Calcular PMF
            ↓
   Validación Cruzada → Verificar consistencia
   ```

---

**Última actualización**: Octubre 2024  
**Configuraciones validadas**: 
- Umbrella: 20 ventanas × 100 ns
- drMD: 50 ns production (HPC optimizado)
