#!/bin/bash
#SBATCH --job-name=wnk_umbrella
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --time=48:00:00
#SBATCH --mem=16G
#SBATCH --array=0-19
#SBATCH --output=logs/window_%a.out
#SBATCH --error=logs/window_%a.err

# ============================================================================
# SLURM Job Array para Umbrella Sampling WNK1
#
# Este script ejecuta 20 ventanas en paralelo usando job array
# Cada trabajo corre una ventana independiente
#
# Uso:
#   sbatch submit_umbrella_hpc.sh
#
# ============================================================================

echo "================================================================================"
echo "WNK1 UMBRELLA SAMPLING - VENTANA ${SLURM_ARRAY_TASK_ID}"
echo "================================================================================"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURM_NODELIST}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "GPU: $CUDA_VISIBLE_DEVICES"
echo "================================================================================"

# ============================================================================
# CONFIGURACIÓN DEL ENTORNO
# ============================================================================

# Cargar módulos del HPC (AJUSTAR SEGÚN TU CLUSTER)
# module load cuda/11.8
# module load python/3.11
# module load openmm/8.0

# Activar entorno virtual (AJUSTAR PATH)
# source /path/to/your/venv/bin/activate

# Variables de entorno OpenMM
export OPENMM_DEFAULT_PLATFORM=CUDA
export OPENMM_CUDA_COMPILER=/usr/local/cuda/bin/nvcc

# Directorio de trabajo (donde están los scripts)
WORK_DIR="${SLURM_SUBMIT_DIR}"
cd "${WORK_DIR}"

echo ""
echo "Directorio de trabajo: ${WORK_DIR}"
echo "Python: $(which python3)"
echo "OpenMM platform: ${OPENMM_DEFAULT_PLATFORM}"
echo ""

# ============================================================================
# VERIFICACIÓN PREVIA
# ============================================================================

WINDOW_ID=${SLURM_ARRAY_TASK_ID}

echo "Verificando archivos necesarios..."

# Verificar que existen los archivos de entrada
if [ ! -f "WNK/run_umbrella_window.py" ]; then
    echo "ERROR: No se encuentra WNK/run_umbrella_window.py"
    exit 1
fi

if [ ! -d "WNK/umbrella_windows/window_$(printf "%02d" ${WINDOW_ID})" ]; then
    echo "ERROR: No se encuentra directorio de ventana ${WINDOW_ID}"
    exit 1
fi

echo "✓ Archivos verificados"
echo ""

# ============================================================================
# EJECUCIÓN DE LA SIMULACIÓN
# ============================================================================

echo "================================================================================"
echo "INICIANDO MD DE PRODUCCIÓN - VENTANA ${WINDOW_ID}"
echo "================================================================================"
echo "Inicio: $(date)"
echo ""

# Parámetros de simulación
STEPS=5000000        # 10 ns con dt=2fs
TEMP=300             # K
PLATFORM=CUDA        # o CPU si no hay GPU
SAVE_INTERVAL=5000   # Guardar cada 10 ps
LOG_INTERVAL=10000   # Log cada 20 ps

# Ejecutar simulación
python3 WNK/run_umbrella_window.py \
    --window ${WINDOW_ID} \
    --steps ${STEPS} \
    --temp ${TEMP} \
    --platform ${PLATFORM} \
    --save-interval ${SAVE_INTERVAL} \
    --log-interval ${LOG_INTERVAL}

EXIT_CODE=$?

echo ""
echo "Fin: $(date)"
echo "Exit code: ${EXIT_CODE}"

# ============================================================================
# VERIFICACIÓN POST-EJECUCIÓN
# ============================================================================

if [ ${EXIT_CODE} -eq 0 ]; then
    echo ""
    echo "================================================================================"
    echo "✅ VENTANA ${WINDOW_ID} COMPLETADA EXITOSAMENTE"
    echo "================================================================================"
    
    # Verificar archivos de salida
    WINDOW_DIR="WNK/umbrella_windows/window_$(printf "%02d" ${WINDOW_ID})"
    
    echo "Archivos generados:"
    ls -lh ${WINDOW_DIR}/trajectory.dcd 2>/dev/null && echo "  ✓ trajectory.dcd" || echo "  ✗ trajectory.dcd FALTANTE"
    ls -lh ${WINDOW_DIR}/cv_values.dat 2>/dev/null && echo "  ✓ cv_values.dat" || echo "  ✗ cv_values.dat FALTANTE"
    ls -lh ${WINDOW_DIR}/production.log 2>/dev/null && echo "  ✓ production.log" || echo "  ✗ production.log FALTANTE"
    
    # Mostrar últimas líneas del log
    if [ -f "${WINDOW_DIR}/production.log" ]; then
        echo ""
        echo "Últimas 5 líneas del log:"
        tail -n 5 ${WINDOW_DIR}/production.log
    fi
    
    # Mostrar estadísticas de CV
    if [ -f "${WINDOW_DIR}/cv_values.dat" ]; then
        echo ""
        echo "Estadísticas de CV:"
        # Usar awk para calcular promedio y desviación
        awk 'NR>1 {sum+=$3; sumsq+=$3*$3; n++} END {
            if (n>0) {
                mean=sum/n; 
                std=sqrt(sumsq/n - mean*mean); 
                printf "  Mean distance: %.4f nm\n  Std dev: %.4f nm\n  Samples: %d\n", mean, std, n
            }
        }' ${WINDOW_DIR}/cv_values.dat
    fi
    
else
    echo ""
    echo "================================================================================"
    echo "❌ VENTANA ${WINDOW_ID} FALLÓ"
    echo "================================================================================"
    echo "Revisa el archivo de error: logs/window_${WINDOW_ID}.err"
    exit ${EXIT_CODE}
fi

echo "================================================================================"
