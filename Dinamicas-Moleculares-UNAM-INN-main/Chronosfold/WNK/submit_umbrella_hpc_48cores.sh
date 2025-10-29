#!/bin/bash
#SBATCH --job-name=wnk_umbrella
#SBATCH --array=0-19               # 20 ventanas (0-19)
#SBATCH --cpus-per-task=4          # 4 CPUs por ventana (optimizado para 48 cores)
#SBATCH --mem=8G                   # 8 GB por ventana
#SBATCH --time=48:00:00            # 48 horas máximo
#SBATCH --output=logs/umbrella_%A_%a.out
#SBATCH --error=logs/umbrella_%A_%a.err
#SBATCH --partition=default        # Ajustar según tu cluster

# ============================================================================
# WNK1 Umbrella Sampling - SLURM Job Array 
# OPTIMIZADO PARA: 2x Intel Xeon E5-2695 v2 (48 cores total)
# ============================================================================
# Hardware Target:
#   - 2 sockets × 12 cores/socket × 2 threads = 48 logical CPUs
#   - L3 cache: 30 MB por socket
#   - NUMA architecture (2 nodes)
#
# Optimization Strategy:
#   - 4 CPUs/task × 12 tasks simultáneas = 48 cores utilizados
#   - Cada ventana usa 4 threads (optimal para OpenMM CPU)
#   - NUMA binding para minimizar latencia de memoria
# ============================================================================

# ============================================================================
# CONFIGURACIÓN DEL ENTORNO
# ============================================================================

# Cargar módulos del HPC (AJUSTAR SEGÚN TU CLUSTER)
# module load python/3.11
# module load openmm/8.0

# Activar entorno virtual (AJUSTAR PATH)
# source /path/to/your/venv/bin/activate

# Variables de entorno OpenMM - OPTIMIZADO PARA CPU
export OPENMM_DEFAULT_PLATFORM=CPU
export OPENMM_CPU_THREADS=4        # 4 threads por ventana

# OpenMP settings para NUMA
export OMP_NUM_THREADS=4
export OMP_PLACES=cores            # Bind threads a cores físicos
export OMP_PROC_BIND=close         # Keep threads close (mismo socket)

# NUMA binding (opcional, descomentar si es necesario)
# export NUMA_NODE=$((SLURM_ARRAY_TASK_ID % 2))  # Alternar entre NUMA nodes
# numactl --cpunodebind=$NUMA_NODE --membind=$NUMA_NODE

# Directorio de trabajo
WORK_DIR=$(dirname $(readlink -f $0))
cd $WORK_DIR

# Crear directorio de logs si no existe
mkdir -p logs

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
echo "HPC: 2x Intel Xeon E5-2695 v2 (48 cores total)"
echo "Platform: CPU (4 threads)"
echo "==============================================="
echo ""

# ============================================================================
# VERIFICACIÓN DE DEPENDENCIAS
# ============================================================================

echo "Verificando dependencias..."

# Verificar Python
if ! command -v python &> /dev/null; then
    echo "ERROR: Python no encontrado"
    exit 1
fi

# Verificar OpenMM
python -c "import openmm" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "ERROR: OpenMM no instalado"
    exit 1
fi

echo "✓ Dependencias OK"
echo ""

# ============================================================================
# EJECUTAR SIMULACIÓN
# ============================================================================

echo "Iniciando simulación de ventana $WINDOW_ID..."
echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Ejecutar simulación
python run_umbrella_window.py \
    --window $WINDOW_ID \
    --steps 50000000 \
    --platform CPU \
    --threads 4

EXIT_CODE=$?

echo ""
echo "==============================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ Ventana $WINDOW_ID completada exitosamente"
    
    # Verificar archivos de salida
    OUT_DIR="umbrella_windows/window_$(printf "%02d" $WINDOW_ID)"
    if [ -f "$OUT_DIR/cv_values.dat" ] && [ -f "$OUT_DIR/final.pdb" ]; then
        echo "✓ Archivos de salida verificados"
        
        # Mostrar estadísticas básicas
        N_FRAMES=$(wc -l < "$OUT_DIR/cv_values.dat")
        echo "  - Frames guardados: $N_FRAMES"
        echo "  - Trayectoria: $OUT_DIR/trajectory.dcd"
        echo "  - Estado final: $OUT_DIR/final.pdb"
    else
        echo "⚠️  ADVERTENCIA: Archivos de salida incompletos"
        EXIT_CODE=2
    fi
else
    echo "✗ ERROR: Ventana $WINDOW_ID falló con código $EXIT_CODE"
fi
echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
echo "==============================================="

exit $EXIT_CODE
