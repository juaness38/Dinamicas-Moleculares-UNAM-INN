#!/bin/bash
#SBATCH --job-name=wnk_prep
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=23:00:00
#SBATCH --output=logs/prepare_%j.out
#SBATCH --error=logs/prepare_%j.err
#SBATCH --partition=q1d-20p
#SBATCH --nodes=1

echo "========================================================================"
echo "WNK1 SYSTEM PREPARATION"
echo "========================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start: $(date)"
echo ""

# Crear directorio de logs
mkdir -p logs

# Activar entorno
source ~/anaconda3/bin/activate
conda activate mdgraphemb

# Configurar OpenMM para CPU
export OPENMM_DEFAULT_PLATFORM=CPU
export OPENMM_CPU_THREADS=4

# Ejecutar preparación
echo "Ejecutando prepare_system.py..."
python -u prepare_system.py  # -u = unbuffered output

echo ""
echo "End: $(date)"
echo "========================================================================"
