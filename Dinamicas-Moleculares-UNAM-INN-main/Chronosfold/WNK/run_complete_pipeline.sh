#!/bin/bash
# ============================================================================
# WNK1 Complete Pipeline with Visualization
# ============================================================================
# Ejecuta el pipeline completo:
#   1. Preparación del sistema (ProPKa + PBS)
#   2. Generación de ventanas
#   3. Simulaciones MD (20 ventanas)
#   4. Análisis MBAR
#   5. Visualización completa (VideoSuite)
#
# Uso:
#   bash run_complete_pipeline.sh [local|hpc]
# ============================================================================

set -e  # Exit on error

PLATFORM=${1:-local}  # Default: local
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║         WNK1 Complete Pipeline with Visualization                 ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""
echo "Platform: $PLATFORM"
echo "Directory: $SCRIPT_DIR"
echo ""

# ============================================================================
# FASE 1: Preparación del Sistema
# ============================================================================
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║  FASE 1: Preparación del Sistema (PBS + ProPKa)                   ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"

if [ ! -f "prepared_system/equilibrated.pdb" ]; then
    echo "→ Ejecutando prepare_system.py..."
    python prepare_system.py
    echo "✓ Sistema preparado"
else
    echo "✓ Sistema ya preparado (skip)"
fi

# Analizar ProPKa
if [ -f "prepared_system/propka_results.pka" ]; then
    echo "→ Analizando estados de protonación..."
    python analyze_propka.py
fi

echo ""

# ============================================================================
# FASE 2: Generación de Ventanas
# ============================================================================
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║  FASE 2: Generación de Ventanas (20 windows)                      ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"

if [ ! -f "umbrella_windows/windows_config.csv" ]; then
    echo "→ Ejecutando generate_umbrella_windows.py..."
    python generate_umbrella_windows.py
    echo "✓ 20 ventanas generadas"
else
    echo "✓ Ventanas ya generadas (skip)"
fi

echo ""

# ============================================================================
# FASE 3: Simulaciones MD
# ============================================================================
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║  FASE 3: Simulaciones MD                                          ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"

if [ "$PLATFORM" = "hpc" ]; then
    echo "→ Submitting SLURM job array..."
    
    # Crear directorio de logs
    mkdir -p logs
    
    # Submit job
    JOB_ID=$(sbatch submit_umbrella_hpc.sh | awk '{print $NF}')
    echo "✓ Job submitted: $JOB_ID"
    echo ""
    echo "Monitorear con:"
    echo "  squeue -u \$USER"
    echo "  tail -f logs/umbrella_${JOB_ID}_*.out"
    echo ""
    echo "Esperando a que terminen los jobs..."
    
    # Esperar a que terminen todos
    while squeue -j $JOB_ID &> /dev/null; do
        RUNNING=$(squeue -j $JOB_ID -h -t RUNNING | wc -l)
        PENDING=$(squeue -j $JOB_ID -h -t PENDING | wc -l)
        echo "  Running: $RUNNING, Pending: $PENDING"
        sleep 60
    done
    
    echo "✓ Todos los jobs completados"
    
elif [ "$PLATFORM" = "local" ]; then
    echo "→ Ejecutando simulaciones localmente (solo primeras 3 ventanas para test)..."
    echo "  ADVERTENCIA: Esto es solo para testing. Usar HPC para producción."
    echo ""
    
    for WINDOW in 0 1 2; do
        if [ ! -f "umbrella_windows/window_$(printf "%02d" $WINDOW)/cv_values.dat" ]; then
            echo "  Ventana $WINDOW..."
            python run_umbrella_window.py --window $WINDOW --steps 50000 --platform CPU
        else
            echo "  Ventana $WINDOW ya existe (skip)"
        fi
    done
    
    echo ""
    echo "⚠️  Solo 3 ventanas ejecutadas (test). Para pipeline completo usar HPC:"
    echo "    bash run_complete_pipeline.sh hpc"
    echo ""
fi

echo ""

# ============================================================================
# FASE 4: Análisis MBAR
# ============================================================================
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║  FASE 4: Análisis MBAR y Cálculo de PMF                           ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"

# Verificar que hay suficientes ventanas
N_WINDOWS=$(ls -d umbrella_windows/window_* 2>/dev/null | wc -l)
N_COMPLETE=$(find umbrella_windows/window_*/cv_values.dat 2>/dev/null | wc -l)

echo "Ventanas totales: $N_WINDOWS"
echo "Ventanas completas: $N_COMPLETE"

if [ $N_COMPLETE -lt 10 ]; then
    echo "⚠️  Insuficientes ventanas completadas ($N_COMPLETE < 10)"
    echo "   Saltando análisis MBAR. Ejecutar más simulaciones primero."
else
    echo "→ Ejecutando análisis MBAR..."
    python analyze_umbrella_mbar.py
    echo "✓ PMF calculado"
fi

echo ""

# ============================================================================
# FASE 5: Visualización con VideoSuite
# ============================================================================
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║  FASE 5: Visualización (VideoSuite + Diagnostics)                 ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"

if [ $N_COMPLETE -ge 10 ]; then
    echo "→ Generando visualizaciones..."
    
    # Ejecutar script de visualización
    if [ -f "visualize_results.py" ]; then
        python visualize_results.py
    else
        echo "⚠️  visualize_results.py no encontrado"
        echo "   Creando visualización básica..."
        
        # Crear visualización básica con matplotlib
        python << 'EOF'
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Cargar PMF
pmf_file = Path("pmf_analysis/pmf_results.csv")
if pmf_file.exists():
    pmf = pd.read_csv(pmf_file)
    
    # Plot PMF
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(pmf['cv'], pmf['pmf'], 'o-', linewidth=2, markersize=8)
    if 'uncertainty' in pmf.columns:
        ax.fill_between(pmf['cv'], 
                        pmf['pmf'] - pmf['uncertainty'],
                        pmf['pmf'] + pmf['uncertainty'],
                        alpha=0.3)
    ax.set_xlabel('Distance C-terminal (nm)', fontsize=14)
    ax.set_ylabel('PMF (kJ/mol)', fontsize=14)
    ax.set_title('WNK1 C-Terminal PMF (PBS buffer, pH 7.4)', fontsize=16)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig('pmf_analysis/pmf_publication.png', dpi=300, bbox_inches='tight')
    print("✓ PMF plot saved: pmf_analysis/pmf_publication.png")
else:
    print("⚠️  pmf_results.csv not found")
EOF
    fi
    
    echo "✓ Visualización completada"
    echo ""
    echo "Archivos generados:"
    echo "  - pmf_analysis/pmf.png"
    echo "  - pmf_analysis/analysis_combined.png"
    echo "  - pmf_analysis/pmf_publication.png"
    
else
    echo "⚠️  Insuficientes ventanas para visualización ($N_COMPLETE < 10)"
fi

echo ""

# ============================================================================
# RESUMEN FINAL
# ============================================================================
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║  📊 RESUMEN DEL PIPELINE                                           ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""
echo "Fase 1 - Preparación:        ✓"
echo "Fase 2 - Ventanas:           ✓"
echo "Fase 3 - Simulaciones:       $([ $N_COMPLETE -ge 20 ] && echo '✓ (20/20)' || echo "⏳ ($N_COMPLETE/20)")"
echo "Fase 4 - MBAR:               $([ -f pmf_analysis/pmf_results.csv ] && echo '✓' || echo '⏳')"
echo "Fase 5 - Visualización:      $([ -f pmf_analysis/pmf.png ] && echo '✓' || echo '⏳')"
echo ""

if [ $N_COMPLETE -ge 20 ] && [ -f pmf_analysis/pmf_results.csv ]; then
    echo "╔═══════════════════════════════════════════════════════════════════╗"
    echo "║                 🎉 PIPELINE COMPLETADO 🎉                          ║"
    echo "╚═══════════════════════════════════════════════════════════════════╝"
    echo ""
    echo "Resultados finales:"
    echo "  📁 prepared_system/         - Sistema preparado (PBS, ProPKa)"
    echo "  📁 umbrella_windows/        - 20 ventanas + trayectorias"
    echo "  📁 pmf_analysis/            - PMF, plots, análisis"
    echo ""
    echo "Próximos pasos:"
    echo "  1. Revisar PMF: display pmf_analysis/pmf_publication.png"
    echo "  2. Analizar trayectorias: pymol umbrella_windows/window_*/final.pdb"
    echo "  3. Comparar con PBS vs estándar (si se ejecutó ambos)"
    echo ""
else
    echo "╔═══════════════════════════════════════════════════════════════════╗"
    echo "║                 ⏳ PIPELINE EN PROGRESO                            ║"
    echo "╚═══════════════════════════════════════════════════════════════════╝"
    echo ""
    echo "Ventanas completadas: $N_COMPLETE/20"
    echo ""
    if [ "$PLATFORM" = "local" ]; then
        echo "Para ejecutar pipeline completo:"
        echo "  bash run_complete_pipeline.sh hpc"
    else
        echo "Esperando simulaciones HPC..."
        echo "  Monitorear: squeue -u \$USER"
    fi
fi

echo ""
echo "═══════════════════════════════════════════════════════════════════"
