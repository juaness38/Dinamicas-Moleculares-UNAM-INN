# Guía de Optimización para HPC (48 cores)

## 🖥️ Hardware Target

**Sistema**: `yoltla1` HPC Cluster

```
Architecture:          x86_64
CPU(s):                48
  On-line CPU(s):      0-47
Thread(s) per core:    2
Core(s) per socket:    12
Socket(s):             2
NUMA node(s):          2
  NUMA node0 CPU(s):   0-11,24-35
  NUMA node1 CPU(s):   12-23,36-47

Model:                 Intel(R) Xeon(R) CPU E5-2695 v2 @ 2.40GHz
CPU MHz:               1200.000 (min) - 3200.000 (max turbo)
L1d cache:             32 KiB (per core)
L1i cache:             32 KiB (per core)
L2 cache:              256 KiB (per core)
L3 cache:              30 MiB (per socket)
```

### Características Clave

- **24 cores físicos** (12 per socket × 2 sockets)
- **48 threads lógicos** (hyperthreading)
- **NUMA architecture** (2 nodos de memoria)
- **AVX, SSE4.2** (optimizaciones SIMD)
- **30 MB L3 cache** compartida por socket

---

## 🎯 Estrategias de Optimización

### 1. OpenMM Umbrella Sampling (CPU)

**Script optimizado**: `submit_umbrella_hpc_48cores.sh`

#### Configuración SLURM

```bash
#SBATCH --array=0-19              # 20 ventanas
#SBATCH --cpus-per-task=4         # 4 CPUs por ventana
#SBATCH --mem=8G                  # 8 GB por ventana
```

#### Paralelización

- **12 ventanas simultáneas** (máximo que cabe en 48 cores con 4 CPUs cada una)
- **4 threads por ventana** (óptimo para OpenMM CPU)
- **Total: 12 × 4 = 48 cores utilizados**

#### Variables de Entorno

```bash
export OPENMM_CPU_THREADS=4
export OMP_NUM_THREADS=4
export OMP_PLACES=cores           # Bind threads a cores físicos
export OMP_PROC_BIND=close        # Mantener threads en mismo socket
```

#### Estimación de Performance

- **Velocidad**: ~50-100 ns/day por ventana (depende del sistema)
- **20 ventanas × 50 ns**: ~10-20 horas wall time (con 12 jobs paralelos)
- **Memoria**: ~6-8 GB por ventana (total ~96 GB)

#### Uso

```bash
# Crear logs
mkdir -p logs

# Submit job array
sbatch submit_umbrella_hpc_48cores.sh

# Monitorear
squeue -u $USER
tail -f logs/umbrella_*.out
```

---

### 2. drMD Parallel Simulations

**Config optimizado**: `drMD_wnk_config_hpc.yaml`

#### Configuración Paralela

```yaml
parallelization:
  parallelCPU: 4              # 4 simulaciones en paralelo
  subprocessCpus: 12          # 12 cores por simulación
  numa_aware: true
  numa_strategy: "socket"     # Una simulación por socket
```

#### Distribución de Recursos

- **4 simulaciones paralelas**
- **12 cores por simulación**
- **Total: 4 × 12 = 48 cores**
- **NUMA binding**: Cada simulación en un socket diferente (minimiza latencia)

#### Estimación de Performance

- **Velocidad single**: ~100 ns/day (12 cores, CPU optimizado)
- **4 paralelas**: ~400 ns/day total throughput
- **50 ns production**: ~12 hours wall time por simulación

#### Uso

```bash
# Verificar drMD instalado
python -c "import drmd; print(drmd.__version__)"

# Ejecutar pipeline
python run_drMD_pipeline.py --config drMD_wnk_config_hpc.yaml

# Monitorear (drMD crea logs automáticamente)
tail -f drMD_simulations/*/drmd.log
```

---

## 🧪 Comparación de Estrategias

| Aspecto | OpenMM Umbrella | drMD Parallel |
|---------|-----------------|---------------|
| **Propósito** | PMF (free energy) | MD production múltiple |
| **Cores usados** | 48 (12 jobs × 4 threads) | 48 (4 sims × 12 cores) |
| **Paralelización** | Job array (SLURM) | Python multiprocessing |
| **Throughput** | 20 ventanas en ~10-20h | 4×50ns en ~12h |
| **Memoria** | ~96 GB total | ~64 GB total |
| **Output** | PMF + trayectorias | 4 trayectorias independientes |
| **Cuándo usar** | Calcular barreras energéticas | Explorar espacio conformacional |

---

## 📊 Pipeline Completo Optimizado

### Script Master: `run_complete_pipeline.sh`

```bash
# Uso HPC
bash run_complete_pipeline.sh hpc

# Uso local (solo testing)
bash run_complete_pipeline.sh local
```

### Flujo de Trabajo

```
1. Preparación (1-2 min)
   ├── prepare_system.py (PBS buffer, ProPKa)
   └── analyze_propka.py

2. Generación de Ventanas (1 min)
   └── generate_umbrella_windows.py (20 windows)

3. Simulaciones MD (10-20 horas HPC) ⭐ OPTIMIZADO
   ├── submit_umbrella_hpc_48cores.sh
   │   ├── 12 jobs simultáneos (48 cores)
   │   ├── 4 threads/job (OpenMM optimizado)
   │   └── ~50 ns por ventana
   └── Espera automática hasta completar

4. Análisis MBAR (5-10 min)
   └── analyze_umbrella_mbar.py (PMF + uncertainty)

5. Visualización (2-3 min) ⭐ NUEVO
   ├── visualize_results.py
   │   ├── PMF publicable (alta resolución)
   │   ├── Histogramas por ventana
   │   ├── Convergencia temporal
   │   ├── Análisis combinado 4-panel
   │   └── Video animado (opcional, --animation)
   └── Outputs:
       ├── pmf_publication.png
       ├── histograms_all_windows.png
       ├── convergence_analysis.png
       ├── analysis_combined.png
       └── umbrella_histograms.mp4 (si --animation)
```

---

## 🔧 Optimizaciones Avanzadas

### NUMA Binding (Opcional)

Para control fino de NUMA:

```bash
# En submit_umbrella_hpc_48cores.sh (descomentar)
export NUMA_NODE=$((SLURM_ARRAY_TASK_ID % 2))
numactl --cpunodebind=$NUMA_NODE --membind=$NUMA_NODE python run_umbrella_window.py ...
```

**Efecto**:
- Ventanas pares (0,2,4...) → NUMA node 0
- Ventanas impares (1,3,5...) → NUMA node 1
- Reduce latencia de memoria ~10-15%

### Hyperthreading

**Actual**: 4 threads/job (usa hyperthreading)
**Alternativa**: 2 threads/job (solo cores físicos)

```bash
#SBATCH --cpus-per-task=2
export OPENMM_CPU_THREADS=2
export OMP_NUM_THREADS=2
```

**Trade-off**:
- ✅ Menos contención de cache
- ❌ Solo 24 jobs paralelos (vs 12)
- ❌ Más overhead de comunicación

**Recomendación**: Mantener 4 threads (usa hyperthreading). OpenMM está optimizado para esto.

### Memory Bandwidth

Si hay contención de memoria (verificar con `htop`):

```bash
# Reducir jobs paralelos
#SBATCH --cpus-per-task=6  # Solo 8 jobs en vez de 12
export OPENMM_CPU_THREADS=6
```

---

## 📈 Monitoreo y Debugging

### Verificar Uso de CPUs

```bash
# En el nodo HPC
htop

# O desde login node
squeue -u $USER -o "%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6D %R %C"
```

### Verificar NUMA Locality

```bash
# Dentro del job
numactl --hardware
numactl --show

# Ver binding de procesos
ps -eo pid,psr,comm | grep python
```

### Logs de Performance

Cada job guarda:
```
logs/umbrella_<JOB_ID>_<WINDOW_ID>.out
logs/umbrella_<JOB_ID>_<WINDOW_ID>.err
```

Buscar:
- ✅ "Speed: XXX ns/day" (velocidad OpenMM)
- ⚠️ "Memory allocation failed" (RAM insuficiente)
- ⚠️ "CUDA_ERROR" (si accidentalmente se usa GPU)

---

## 🎬 Visualización Completa

### Generar Todas las Visualizaciones

```bash
# Después de MBAR
python visualize_results.py --animation

# Output:
# - pmf_publication.png (PMF alta resolución)
# - histograms_all_windows.png (20 histogramas)
# - convergence_analysis.png (convergencia temporal)
# - analysis_combined.png (4-panel diagnóstico)
# - umbrella_histograms.mp4 (video animado) ⭐ requiere ffmpeg
```

### Requisitos para Video

```bash
# Instalar ffmpeg (si no está)
conda install ffmpeg

# O en HPC
module load ffmpeg
```

### Integración con VIDEOSUITE

El script `visualize_results.py` intentará usar:
```python
from VIDEOSUITE.umbrella.run_umbrella_visualization import generate_animation
```

Si no encuentra VIDEOSUITE:
```bash
# Agregar al PYTHONPATH
export PYTHONPATH="/path/to/Chronosfold:$PYTHONPATH"
```

---

## ✅ Checklist de Deployment

Antes de ejecutar en HPC:

- [ ] Verificar módulos disponibles: `module avail`
- [ ] Cargar Python 3.x: `module load python/3.11`
- [ ] Activar environment: `source venv/bin/activate`
- [ ] Verificar OpenMM: `python -c "import openmm; print(openmm.__version__)"`
- [ ] Crear directorio logs: `mkdir -p logs`
- [ ] Verificar PBS system listo: `ls prepared_system/equilibrated.pdb`
- [ ] Ajustar partition en SLURM: `#SBATCH --partition=<tu_partition>`
- [ ] Verificar permisos: `chmod +x *.sh`
- [ ] Test local (3 ventanas): `bash run_complete_pipeline.sh local`
- [ ] Deploy HPC: `bash run_complete_pipeline.sh hpc`

---

## 🐛 Troubleshooting

### Problema: Jobs muy lentos

**Diagnóstico**:
```bash
# Ver velocidad en logs
grep "Speed:" logs/umbrella_*.out
```

**Soluciones**:
- Reducir `--cpus-per-task` si hay contención
- Verificar que no se está usando GPU (debe ser CPU)
- Aumentar `--cpus-per-task` a 6 (reduce jobs paralelos)

### Problema: Out of Memory

**Síntomas**: Jobs cancelados con "OOM" en stderr

**Soluciones**:
```bash
#SBATCH --mem=12G  # Aumentar de 8G a 12G
```

### Problema: NUMA performance pobre

**Verificar**:
```bash
# Dentro del job
numastat
```

**Si hay muchos "numa_miss"**, activar NUMA binding (ver sección avanzada).

---

## 📚 Referencias

- [OpenMM CPU Platform](http://docs.openmm.org/latest/userguide/library/02_running_sims.html#platforms)
- [SLURM Job Arrays](https://slurm.schedmd.com/job_array.html)
- [Intel Xeon E5-2695 v2 Specs](https://ark.intel.com/content/www/us/en/ark/products/75283/intel-xeon-processor-e52695-v2-30m-cache-2-40-ghz.html)
- [NUMA Best Practices](https://documentation.suse.com/sles/15-SP1/html/SLES-all/cha-tuning-numactl.html)

---

## 🚀 Quick Start

```bash
# 1. Preparar todo
cd Chronosfold/WNK

# 2. Ejecutar pipeline completo (HPC)
bash run_complete_pipeline.sh hpc

# 3. Esperar completar (10-20h)
squeue -u $USER

# 4. Ver resultados
display pmf_analysis/pmf_publication.png
display pmf_analysis/analysis_combined.png

# 5. (Opcional) Generar video
python visualize_results.py --animation
```

---

**Última actualización**: 2024
**Hardware validado**: Intel Xeon E5-2695 v2 (48 cores)
**Software**: OpenMM 8.0+, Python 3.11+
