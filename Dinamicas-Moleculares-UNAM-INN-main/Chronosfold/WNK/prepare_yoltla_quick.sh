#!/bin/bash
# ============================================================================
# Script de preparación rápida para WNK Umbrella Sampling en Yoltla
# ============================================================================

set -e  # Exit on error

WORK_DIR="/LUSTRE/home/lancad/2025/l.100066/Dinamicas-Moleculares-UNAM-INN-main/Chronosfold/WNK"

echo "========================================================================"
echo "  WNK Umbrella Sampling - Preparación Rápida (Yoltla HPC)"
echo "========================================================================"
echo ""

# Activar entorno conda
echo "[1/6] Activando entorno mdgraphemb..."
source ~/anaconda3/bin/activate
conda activate mdgraphemb
echo "✓ Entorno activado"
echo "  Python: $(python --version)"
echo "  OpenMM: $(python -c 'import openmm; print(openmm.__version__)')"
echo ""

# Ir al directorio de trabajo
cd $WORK_DIR

# Verificar estructura PDB
echo "[2/6] Verificando estructura 5DRB.pdb..."
if [ ! -f "5DRB.pdb" ]; then
    echo "✗ ERROR: 5DRB.pdb no encontrado"
    exit 1
fi
echo "✓ Estructura encontrada ($(stat -c%s 5DRB.pdb) bytes)"
echo ""

# Instalar dependencias faltantes si es necesario
echo "[3/6] Verificando dependencias..."

# Check pymbar (necesario para MBAR analysis)
python -c "import pymbar" 2>/dev/null || {
    echo "  Instalando pymbar..."
    pip install pymbar --quiet
}

# Check propka (necesario para protonación)
python -c "import propka" 2>/dev/null || {
    echo "  Instalando propka..."
    pip install propka --quiet
}

# Check pdbfixer (necesario para preparación)
python -c "import pdbfixer" 2>/dev/null || {
    echo "  Instalando pdbfixer..."
    pip install pdbfixer --quiet
}

echo "✓ Todas las dependencias instaladas"
echo ""

# Preparar sistema (si no existe)
echo "[4/6] Preparando sistema molecular..."
if [ -d "system_prepared" ] && [ -f "system_prepared/solvated_ionized.pdb" ]; then
    echo "  Sistema ya preparado, saltando..."
else
    echo "  Ejecutando prepare_system.py..."
    python prepare_system.py
    
    if [ $? -eq 0 ]; then
        echo "✓ Sistema preparado exitosamente"
    else
        echo "✗ ERROR: Falló la preparación del sistema"
        exit 1
    fi
fi
echo ""

# Generar ventanas de umbrella (si no existen)
echo "[5/6] Generando ventanas de umbrella sampling..."
if [ -d "umbrella_windows" ] && [ -f "umbrella_windows/window_0.xml" ]; then
    WINDOW_COUNT=$(ls umbrella_windows/window_*.xml 2>/dev/null | wc -l)
    echo "  Ya existen $WINDOW_COUNT ventanas, saltando..."
else
    echo "  Ejecutando generate_umbrella_windows.py..."
    python generate_umbrella_windows.py
    
    if [ $? -eq 0 ]; then
        WINDOW_COUNT=$(ls umbrella_windows/window_*.xml 2>/dev/null | wc -l)
        echo "✓ $WINDOW_COUNT ventanas generadas exitosamente"
    else
        echo "✗ ERROR: Falló la generación de ventanas"
        exit 1
    fi
fi
echo ""

# Test rápido de una ventana (opcional)
echo "[6/6] Test rápido de umbrella window..."
read -p "¿Ejecutar test de 1000 steps en ventana 10? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "  Ejecutando test corto..."
    python run_umbrella_window.py \
        --window-id 10 \
        --input-pdb system_prepared/solvated_ionized.pdb \
        --window-state umbrella_windows/window_10.xml \
        --output-dir test_output \
        --steps 1000 \
        --temperature 310.0 \
        --platform CPU \
        --cpu-threads 4 \
        --report-interval 100
    
    if [ $? -eq 0 ]; then
        echo "✓ Test completado exitosamente"
        ls -lh test_output/
    else
        echo "✗ Test falló"
    fi
else
    echo "  Test saltado"
fi
echo ""

# Resumen final
echo "========================================================================"
echo "  PREPARACIÓN COMPLETADA"
echo "========================================================================"
echo ""
echo "Sistema listo para umbrella sampling en Yoltla HPC"
echo ""
echo "Archivos generados:"
echo "  - system_prepared/solvated_ionized.pdb  (sistema solvatado + PBS)"
echo "  - umbrella_windows/window_*.xml         (20 ventanas)"
echo ""
echo "Para ejecutar umbrella sampling completo:"
echo "  1. Transferir submit_umbrella_yoltla_zen3.sh a Yoltla"
echo "  2. Ejecutar: sbatch submit_umbrella_yoltla_zen3.sh"
echo ""
echo "Esto lanzará 20 jobs paralelos (16 simultáneos en 64 cores)"
echo "Tiempo estimado: ~24-48 horas"
echo ""
echo "========================================================================"
