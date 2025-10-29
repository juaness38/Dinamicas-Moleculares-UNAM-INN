#!/bin/bash
#SBATCH --job-name=wnk_test_100ns
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=23:00:00
#SBATCH --output=logs/test_100ns_%j.out
#SBATCH --error=logs/test_100ns_%j.err
#SBATCH --partition=q1d-20p
#SBATCH --nodes=1

echo "========================================================================"
echo "WNK1 - TEST PRODUCCIÓN 100 ns"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"
echo ""

# Crear directorio de logs
mkdir -p logs production_test_100ns

# Activar entorno
source ~/anaconda3/bin/activate
conda activate mdgraphemb

# Configurar OpenMM para CPU
export OPENMM_DEFAULT_PLATFORM=CPU
export OPENMM_CPU_THREADS=4
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# Ejecutar simulación con unbuffered output
echo "Ejecutando test_production_100ns.py..."
python -u test_production_100ns.py

echo ""
echo "End: $(date)"
echo "========================================================================"
