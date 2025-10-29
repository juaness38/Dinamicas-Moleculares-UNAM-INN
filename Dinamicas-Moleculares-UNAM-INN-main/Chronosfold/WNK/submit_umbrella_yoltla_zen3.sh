#!/bin/bash
#SBATCH --job-name=wnk_umbrella_zen3
#SBATCH --array=0-19               # 20 ventanas (0-19)
#SBATCH --cpus-per-task=4          # 4 CPUs por ventana (16 ventanas paralelas en 64 cores)
#SBATCH --mem=8G                   # 8 GB por ventana
#SBATCH --time=48:00:00            # 48 horas máximo
#SBATCH --output=logs/umbrella_%A_%a.out
#SBATCH --error=logs/umbrella_%A_%a.err
#SBATCH --partition=qz2d-128p      # Partición AMD Zen3 (2 días, 128 cores)
#SBATCH --nodes=1                  # 1 nodo (64 cores suficiente)

# ============================================================================
# WNK1 Umbrella Sampling - SLURM Job Array para Yoltla HPC
# OPTIMIZADO PARA: 2x AMD EPYC 7003 (Zen3) - 64 cores totales
# ============================================================================
# Hardware Target:
#   - 2 sockets × 32 cores/socket = 64 cores físicos
#   - RAM: 512 GB
#   - L3 cache: 256 MB total (32 MB × 8 CCDs)
#   - NUMA architecture (2 nodes, 8 CCDs)
#
# Optimization Strategy:
#   - 4 CPUs/task × 16 tasks simultáneas = 64 cores utilizados
#   - Cada ventana usa 4 threads (optimal para OpenMM CPU)
#   - NUMA binding para minimizar latencia de memoria
#   - AVX2/AVX512 vectorization enabled
# ============================================================================

# ============================================================================
# CONFIGURACIÓN DEL ENTORNO - YOLTLA SPECIFIC
# ============================================================================

# NO cargar módulos system python (GLIBC 2.17 issues)
# Usar Anaconda3 con entorno mdgraphemb

# Activar entorno conda mdgraphemb (ya tiene OpenMM 8.0)
source ~/anaconda3/bin/activate
conda activate mdgraphemb

# Variables de entorno OpenMM - OPTIMIZADO PARA CPU AMD Zen3
export OPENMM_DEFAULT_PLATFORM=CPU
export OPENMM_CPU_THREADS=4        # 4 threads por ventana

# OpenMP settings para NUMA y AVX2
export OMP_NUM_THREADS=4
export OMP_PLACES=cores            # Bind threads a cores físicos
export OMP_PROC_BIND=close         # Keep threads close (mismo CCD)
export OMP_WAIT_POLICY=ACTIVE      # Active spinning para baja latencia

# AMD Zen3 specific optimizations
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

# Directorio de trabajo
WORK_DIR="/LUSTRE/home/lancad/2025/l.100066/Dinamicas-Moleculares-UNAM-INN-main/Chronosfold/WNK"
cd $WORK_DIR

# Crear directorio de logs si no existe
mkdir -p logs
mkdir -p umbrella_windows

# ============================================================================
# INFORMACIÓN DEL JOB
# ============================================================================

WINDOW_ID=$SLURM_ARRAY_TASK_ID

echo "==============================================="
echo "SLURM Job Array: WNK1 Umbrella Sampling"
echo "==============================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID (ventana $WINDOW_ID)"
echo "Node: $SLURMD_NODENAME"
echo "CPUs asignados: $SLURM_CPUS_PER_TASK"
echo "Memoria asignada: 8 GB"
echo "HPC: 2x AMD EPYC 7003 Zen3 (64 cores total)"
echo "Platform: CPU (4 threads, AVX2 optimized)"
echo "Conda env: mdgraphemb (OpenMM 8.0)"
echo "==============================================="
echo ""

# ============================================================================
# VERIFICACIÓN DE DEPENDENCIAS
# ============================================================================

echo "Verificando dependencias..."

# Check Python
python --version || { echo "ERROR: Python no disponible"; exit 1; }

# Check OpenMM
python -c "import openmm; print(f'OpenMM: {openmm.__version__}')" || \
    { echo "ERROR: OpenMM no disponible"; exit 1; }

# Check MDAnalysis
python -c "import MDAnalysis; print(f'MDAnalysis: {MDAnalysis.__version__}')" || \
    { echo "ERROR: MDAnalysis no disponible"; exit 1; }

# Check pymbar
python -c "import pymbar; print(f'pymbar: {pymbar.__version__}')" 2>/dev/null || \
    { echo "ADVERTENCIA: pymbar no instalado (necesario para análisis MBAR)"; }

echo "✓ Dependencias verificadas"
echo ""

# ============================================================================
# VERIFICACIÓN DE ARCHIVOS DE ENTRADA
# ============================================================================

echo "Verificando archivos de entrada..."

# Check si existe el sistema preparado
if [ ! -f "system_prepared/solvated_ionized.pdb" ]; then
    echo "ERROR: Sistema no preparado. Ejecuta primero:"
    echo "  python prepare_system.py"
    exit 1
fi

# Check si existen las ventanas generadas
if [ ! -d "umbrella_windows" ] || [ ! -f "umbrella_windows/window_${WINDOW_ID}.xml" ]; then
    echo "ERROR: Ventanas no generadas. Ejecuta primero:"
    echo "  python generate_umbrella_windows.py"
    exit 1
fi

echo "✓ Archivos de entrada verificados"
echo ""

# ============================================================================
# EJECUCIÓN DE UMBRELLA SAMPLING
# ============================================================================

echo "Iniciando umbrella sampling para ventana $WINDOW_ID..."
echo "Distancia inicial: $(grep 'distance' umbrella_windows/window_${WINDOW_ID}.xml | head -1)"
echo ""

# Ejecutar simulación de producción con bias armónico
python run_umbrella_window.py \
    --window-id $WINDOW_ID \
    --input-pdb system_prepared/solvated_ionized.pdb \
    --window-state umbrella_windows/window_${WINDOW_ID}.xml \
    --output-dir umbrella_windows \
    --steps 50000000 \
    --temperature 310.0 \
    --pressure 1.0 \
    --report-interval 5000 \
    --dcd-interval 5000 \
    --checkpoint-interval 25000 \
    --platform CPU \
    --cpu-threads 4

EXIT_CODE=$?

# ============================================================================
# VERIFICACIÓN DE RESULTADOS
# ============================================================================

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "✓ Ventana $WINDOW_ID completada exitosamente"
    
    # Verificar que se generó el archivo de trayectoria
    if [ -f "umbrella_windows/window_${WINDOW_ID}_prod.dcd" ]; then
        DCD_SIZE=$(stat -f%z "umbrella_windows/window_${WINDOW_ID}_prod.dcd" 2>/dev/null || \
                   stat -c%s "umbrella_windows/window_${WINDOW_ID}_prod.dcd" 2>/dev/null)
        echo "  Trayectoria: $(($DCD_SIZE / 1024 / 1024)) MB"
    fi
    
    # Verificar checkpoint final
    if [ -f "umbrella_windows/window_${WINDOW_ID}_checkpoint.chk" ]; then
        echo "  Checkpoint: ✓ guardado"
    fi
    
else
    echo ""
    echo "✗ ERROR: Ventana $WINDOW_ID falló con código $EXIT_CODE"
    exit $EXIT_CODE
fi

echo ""
echo "==============================================="
echo "Job completado: $(date)"
echo "==============================================="
