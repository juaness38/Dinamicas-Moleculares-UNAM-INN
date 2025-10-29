# WNK1 C-Terminal Umbrella Sampling Pipeline

Pipeline completo de umbrella sampling para explorar cambios conformacionales del C-terminal de WNK1 kinasa usando MBAR en **condiciones fisiológicas PBS**.

## 📋 Contenido

- **Sistema**: WNK1 kinasa (rat, PDB: 5DRB, residuos 194-483)
- **Variable colectiva (CV)**: Distancia entre centros de masa del dominio kinasa y C-terminal
- **Condiciones**: PBS buffer (137 mM NaCl, 2.7 mM KCl, fosfatos 10/1.8 mM), pH 7.4
- **Método de análisis**: MBAR (Multistate Bennett Acceptance Ratio)
- **Plataforma**: HPC con SLURM job arrays + drMD pipeline paralelo

## 🗂️ Estructura de archivos

```
WNK/
├── 5DRB.pdb                       # Estructura cristalográfica de WNK1
│
├── PIPELINE PRINCIPAL (Chronosfold)
├── prepare_system.py              # 1. Preparación con ProPKa + PBS buffer
├── generate_umbrella_windows.py   # 2. Generación de 20 ventanas
├── run_umbrella_window.py         # 3. MD producción con bias armónico
├── analyze_umbrella_mbar.py       # 4. Análisis MBAR y PMF
├── submit_umbrella_hpc.sh         # SLURM job array
├── deploy_to_hpc.sh               # SSH deployment automation
│
├── PIPELINE PARALELO (drMD)
├── drMD_wnk_config.yaml           # Configuración drMD con PBS
├── run_drMD_pipeline.py           # Script para ejecutar drMD
│
├── DOCUMENTACIÓN
├── README.md                      # Este archivo
├── PROTONACION_GUIDE.md           # Guía ProPKa y estados HIS
├── PBS_BUFFER_IMPLEMENTATION.md   # Implementación PBS y limitaciones
├── PIPELINE_DIAGRAM.md            # Diagrama visual completo
└── PROPKA_INTEGRATION_SUMMARY.md  # Resumen integración ProPKa

tests/
├── test_wnk_umbrella_setup.py     # Tests pre-HPC (sistema, windows, bias)
└── test_pbs_conditions.py         # Tests PBS (fuerza iónica, pH, iones)
```

## 🚀 Uso Rápido

### Opción 1: Pipeline drMD (Automatizado, PBS real) ⭐ NUEVO

drMD ofrece automatización completa con manejo de errores y clustering.

1. **Verificar instalación drMD**:
   ```bash
   # drMD ya está en SciToolAgent/external/drMD
   cd ../../SciToolAgent/external/drMD
   pip install -e .
   cd -
   ```

2. **Ejecutar pipeline drMD**:
   ```bash
   python run_drMD_pipeline.py
   ```

   **Ventajas**:
   - ✅ PBS buffer con pH 7.4 automático
   - ✅ FirstAid (manejo automático de errores)
   - ✅ Clustering de conformaciones
   - ✅ Generación de sección de métodos
   - ✅ Logging detallado

   **Output**: `drMD_output/5DRB/production/`

3. **Revisar resultados**:
   ```bash
   # Ver logs
   cat drMD_output/00_drMD_logs/aftercare.log
   
   # Visualizar
   pymol drMD_output/5DRB/production/production_*.pdb
   
   # Ver clusters
   ls drMD_output/5DRB/05_cluster_analysis/clusters/
   ```

### Opción 2: Despliegue automático a HPC (Chronosfold)

Pipeline principal con umbrella sampling completo.

1. **Configurar credenciales HPC**:
   ```bash
   # Editar deploy_to_hpc.sh
   HPC_USER="tu_usuario"
   HPC_HOST="nombre_cluster.edu"
   HPC_WORK_DIR="/home/tu_usuario/wnk_umbrella"
   HPC_PYTHON_ENV="/path/to/venv"
   ```

2. **Ejecutar despliegue**:
   ```bash
   bash deploy_to_hpc.sh
   # Seleccionar opción 1 (Deploy completo)
   ```

3. **Monitorear jobs**:
   ```bash
   bash deploy_to_hpc.sh
   # Seleccionar opción 7 (Monitorear)
   ```

4. **Descargar resultados**:
   ```bash
   bash deploy_to_hpc.sh
   # Seleccionar opción 6 (Descargar resultados)
   ```

### Opción 2: Ejecución manual paso a paso

#### Paso 1: Preparar sistema

```bash
python prepare_system.py
```

**Output**:
- `prepared_system/propka_results.pka` - Análisis de pKa de residuos (ProPKa)
- `prepared_system/system_solvated.pdb` - Sistema completo solvatado
- `prepared_system/equilibrated.pdb` - Sistema equilibrado
- `prepared_system/equilibrated_state.xml` - Estado listo para producción

**Pasos internos**:
1. Cargar estructura 5DRB.pdb
2. Limpiar estructura (solo proteína)
3. **Determinar estados de protonación con ProPKa (pH 7.4)** ← NUEVO
4. Agregar hidrógenos con estados correctos
5. **Solvatar en PBS buffer** (137 mM NaCl equiv., pH 7.4) ← NUEVO
6. Minimización de energía
7. Equilibración NVT (100 ps)
8. Equilibración NPT (100 ps)

**⚠️ CONDICIONES PBS**:
- 137 mM NaCl
- 2.7 mM KCl (aproximado como Na+, ver `PBS_BUFFER_IMPLEMENTATION.md`)
- 10 mM Na₂HPO₄ (aproximado como Cl⁻, forcefield limitation)
- 1.8 mM KH₂PO₄ (aproximado como Cl⁻, forcefield limitation)
- Fuerza iónica total: ~163 mM
- pH: 7.4 (fisiológico)

**⚠️ IMPORTANTE**: 
- Revisar `propka_results.pka` para verificar estados de protonación de HIS (histidina)
- Ver `PROTONACION_GUIDE.md` para detalles de protonación
- Ver `PBS_BUFFER_IMPLEMENTATION.md` para limitaciones de forcefield

**Tiempo estimado**: 10-30 minutos

#### Paso 2: Generar ventanas

```bash
python generate_umbrella_windows.py
```

**Output**:
- `umbrella_windows/windows_config.csv` - Configuración de 20 ventanas
- `umbrella_windows/atom_groups.txt` - Índices de átomos para CV
- `umbrella_windows/window_XX/` - 20 directorios de ventanas

**Ventanas generadas**: 
- Rango: 1.5 - 4.0 nm
- N = 20 ventanas
- k = 1000 kJ/mol/nm²

#### Paso 3: Ejecutar simulaciones

**Localmente (prueba corta)**:
```bash
# Ventana 0, 100 ps
python run_umbrella_window.py --window 0 --steps 50000 --platform CPU
```

**En HPC (producción)**:
```bash
# Crear directorio de logs
mkdir -p logs

# Submit job array (20 ventanas en paralelo)
sbatch submit_umbrella_hpc.sh
```

**Monitorear**:
```bash
squeue -u $USER
```

**Output por ventana**:
- `window_XX/trajectory.dcd` - Trayectoria MD
- `window_XX/cv_values.dat` - Valores de CV en cada frame
- `window_XX/production.log` - Log termodinámico
- `window_XX/final.pdb` - Estructura final

**Tiempo estimado por ventana**: 2-8 horas (depende de GPU/CPU)

#### Paso 4: Análisis MBAR

```bash
python analyze_umbrella_mbar.py
```

**Output**:
- `pmf_analysis/pmf_results.csv` - PMF y errores
- `pmf_analysis/pmf.png` - Plot de PMF
- `pmf_analysis/analysis_combined.png` - Histogramas + PMF

**Requisito**: pymbar instalado (`pip install pymbar`)

## 🧪 Tests Pre-HPC

Antes de enviar a HPC, ejecutar tests de validación:

```bash
cd ..
python -m pytest tests/test_wnk_umbrella_setup.py -v
```

**Tests incluidos**:
- ✓ PDB válido y cargable
- ✓ Sistema preparado tiene agua e iones
- ✓ Ventanas configuradas correctamente
- ✓ Bias armónico funciona
- ✓ Simulación corta (100 pasos) funciona

## 📊 Parámetros de Simulación

### Sistema
- **Forcefield**: AMBER14
- **Agua**: TIP3P
- **Fuerza iónica**: 0.15 M (fisiológica)
- **Caja**: Padding 1.0 nm

### MD
- **Temperatura**: 300 K
- **Presión**: 1 bar (NPT en equilibración)
- **Integrador**: Langevin Middle
- **Timestep**: 2 fs
- **Duración**: 10 ns por ventana (5M pasos)

### Umbrella Sampling
- **CV**: Distancia entre COM(dominio kinasa) y COM(C-terminal)
- **Rango**: 1.5 - 4.0 nm
- **Ventanas**: 20
- **k (spring)**: 1000 kJ/mol/nm²
- **Frecuencia guardado**: cada 5000 pasos (10 ps)

### MBAR
- **Bins**: 50 (ajustable con `--bins`)
- **Equilibración**: 0 frames por defecto (ajustable con `--equilibration`)
- **Submuestreo**: 1 (cada frame, ajustable con `--subsample`)

## 🔧 Configuración HPC

### Recursos SLURM (submit_umbrella_hpc.sh)

```bash
#SBATCH --partition=gpu        # Partición GPU
#SBATCH --gres=gpu:1           # 1 GPU por ventana
#SBATCH --cpus-per-task=4      # 4 CPUs
#SBATCH --mem=16G              # 16 GB RAM
#SBATCH --time=48:00:00        # 48 horas máximo
#SBATCH --array=0-19           # 20 ventanas (job array)
```

### Modificar para tu cluster

Editar `submit_umbrella_hpc.sh`:

```bash
# Cargar módulos específicos de tu cluster
module load cuda/11.8
module load python/3.11
module load openmm/8.0

# O activar virtualenv
source /path/to/your/venv/bin/activate
```

## 📦 Dependencias

### Instalación local

```bash
pip install numpy pandas matplotlib seaborn openmm mdtraj pytest
pip install pymbar  # Para análisis MBAR
pip install propka  # Para análisis de protonación
```

### Instalación en HPC

```bash
# Crear virtualenv
python3 -m venv ~/venv_umbrella
source ~/venv_umbrella/bin/activate

# Instalar dependencias
pip install --upgrade pip
pip install numpy pandas matplotlib seaborn
pip install openmm  # o usar module load si está disponible
pip install mdtraj pymbar pytest
pip install propka  # Para análisis de protonación
```

## 🐛 Troubleshooting

### Error: "Platform CUDA not available"

Si no hay GPU disponible:
```bash
python run_umbrella_window.py --window 0 --platform CPU
```

En SLURM, cambiar:
```bash
#SBATCH --partition=cpu
# Eliminar: #SBATCH --gres=gpu:1
```

### Error: "pymbar not installed"

```bash
pip install pymbar
```

### Error: "OpenMM cannot find CUDA compiler"

```bash
export OPENMM_CUDA_COMPILER=/usr/local/cuda/bin/nvcc
```

### Simulación muy lenta

1. Verificar que está usando GPU:
   ```python
   # En el log debería decir:
   # "✓ Plataforma: CUDA"
   ```

2. Reducir sistema (menos agua):
   ```python
   # En prepare_system.py, cambiar:
   BOX_PADDING = 0.8 * unit.nanometer  # En vez de 1.0
   ```

3. Timestep más grande (solo si no hay constraints):
   ```bash
   python run_umbrella_window.py --dt 0.004  # 4 fs en vez de 2 fs
   ```

### Jobs HPC fallan inmediatamente

1. Verificar logs:
   ```bash
   cat logs/window_0.err
   ```

2. Tests pre-HPC:
   ```bash
   python -m pytest tests/test_wnk_umbrella_setup.py -v
   ```

3. Verificar módulos cargados:
   ```bash
   module list
   which python3
   ```

## 📈 Análisis de Resultados

### Verificar convergencia

```bash
# En cada ventana, verificar que <r> ≈ r₀
# Ejemplo para ventana 5:
awk 'NR>1 {sum+=$3; n++} END {print "Mean r:", sum/n, "nm"}' \
    umbrella_windows/window_05/cv_values.dat

# Comparar con r₀ en:
head -n 7 umbrella_windows/windows_config.csv | tail -n 1
```

### Plot de overlap entre ventanas

```python
import pandas as pd
import matplotlib.pyplot as plt

for i in range(20):
    data = pd.read_csv(f'umbrella_windows/window_{i:02d}/cv_values.dat',
                       sep='\s+', comment='#',
                       names=['step', 'time', 'r', 'bias'])
    plt.hist(data['r'], bins=30, alpha=0.5, label=f'W{i}')

plt.xlabel('Distance (nm)')
plt.ylabel('Counts')
plt.title('Window overlap')
plt.legend(ncol=4, fontsize=6)
plt.savefig('window_overlap.png', dpi=150)
```

### Interpretar PMF

```python
import pandas as pd
import matplotlib.pyplot as plt

pmf = pd.read_csv('pmf_analysis/pmf_results.csv')

# Mínimo global
min_idx = pmf['pmf_kJ_mol'].idxmin()
print(f"Mínimo global: {pmf.loc[min_idx, 'cv_nm']:.3f} nm")

# Barrera energética
max_pmf = pmf['pmf_kJ_mol'].max()
print(f"Barrera: {max_pmf:.2f} kJ/mol ({max_pmf / 2.479:.2f} kcal/mol)")

# Plot
plt.errorbar(pmf['cv_nm'], pmf['pmf_kJ_mol'], 
             yerr=pmf['pmf_error_kJ_mol'], marker='o')
plt.xlabel('Distance (nm)')
plt.ylabel('PMF (kJ/mol)')
plt.title('WNK1 C-terminal PMF')
plt.grid(True, alpha=0.3)
plt.savefig('pmf_interpreted.png', dpi=300)
```

## 📚 Referencias

### Variable Colectiva
- **Dominio kinasa**: Residuos 194-450 (estructura catalítica conservada)
- **C-terminal**: Residuos 451-483 (región regulatoria)
- **Hipótesis**: C-terminal puede adoptar conformaciones compactas vs extendidas que modulan actividad

### Método MBAR
- Shirts & Chodera (2008). *J. Chem. Phys.* 129, 124105
- Precisión típica: <0.5 kcal/mol
- Mejor que WHAM para sistemas con overlap variable

### Sistema WNK1
- PDB: 5DRB (Min et al., 2004)
- Resolución: 1.65 Å
- Especie: Rattus norvegicus (rat)
- Mutaciones oncogénicas frecuentes en C-terminal

## 📧 Soporte

Para preguntas o issues:
1. Verificar tests: `pytest tests/test_wnk_umbrella_setup.py -v`
2. Revisar logs en `logs/window_XX.err`
3. Consultar documentación de OpenMM: http://openmm.org

## 📝 To-Do / Extensiones

- [ ] SMD pulling para generar conformaciones iniciales extendidas
- [ ] Análisis de trayectoria con mdtraj (RMSD, Rg, RMSF)
- [ ] Visualización VMD/PyMOL de trayectorias
- [ ] Análisis de puentes de hidrógeno C-terminal ↔ kinase domain
- [ ] Comparación con mutantes patogénicos
- [ ] Free energy surface 2D (distancia + ángulo)

## ✅ Checklist de Ejecución

**Pre-HPC**:
- [ ] PDB 5DRB.pdb presente
- [ ] Dependencies instaladas
- [ ] Tests pasan: `pytest tests/test_wnk_umbrella_setup.py`
- [ ] Credenciales HPC configuradas en `deploy_to_hpc.sh`

**Ejecución**:
- [ ] Sistema preparado: `prepare_system.py` completado
- [ ] Ventanas generadas: 20 directorios en `umbrella_windows/`
- [ ] Jobs submitted: `sbatch submit_umbrella_hpc.sh`
- [ ] Monitorear: `squeue -u $USER`

**Post-HPC**:
- [ ] 20/20 ventanas completaron (verificar `trajectory.dcd`)
- [ ] CV data presente: `cv_values.dat` en cada ventana
- [ ] Análisis MBAR ejecutado: `analyze_umbrella_mbar.py`
- [ ] PMF interpretado: mínimo, barrera, rango

---

**Creado**: 2025
**Pipeline**: Umbrella Sampling + MBAR
**Sistema**: WNK1 C-terminal conformational change
**Para**: Presentación INN / Publicación
