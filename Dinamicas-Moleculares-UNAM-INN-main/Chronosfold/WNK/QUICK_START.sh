#!/bin/bash
# ============================================================================
# Quick Start Guide for WNK1 PBS + drMD Pipeline
# ============================================================================
# Este script NO se ejecuta automáticamente. Es una guía paso a paso.
# Copia y pega los comandos según necesites.
# ============================================================================

echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║     WNK1 PBS Buffer + drMD Pipeline - Quick Start Guide            ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""

# ============================================================================
# PASO 0: Verificar dependencias
# ============================================================================
echo "PASO 0: Verificar dependencias"
echo "────────────────────────────────────────────────────────────────────"

# Verificar Python
python --version  # Debe ser 3.8+

# Verificar OpenMM
python -c "import openmm; print(f'OpenMM {openmm.__version__}')"

# Verificar ProPKa
python -c "import propka; print('ProPKa OK')" || echo "⚠️  Instalar: pip install propka"

# Verificar MDAnalysis
python -c "import MDAnalysis; print('MDAnalysis OK')" || echo "⚠️  Instalar: pip install MDAnalysis"

# Verificar pymbar
python -c "import pymbar; print('pymbar OK')" || echo "⚠️  Instalar: pip install pymbar"

echo ""

# ============================================================================
# PASO 1: Verificar si WNK1 tiene K+ cristalográfico
# ============================================================================
echo "PASO 1: Verificar K+ en estructura"
echo "────────────────────────────────────────────────────────────────────"

grep "^HETATM.*  K  " 5DRB.pdb

if [ $? -eq 0 ]; then
    echo ""
    echo "⚠️  ENCONTRADO K+ en estructura!"
    echo "   Considerar usar CHARMM36 forcefield en lugar de amber14"
    echo "   Ver: PBS_BUFFER_IMPLEMENTATION.md sección OPCIÓN 2"
    echo ""
else
    echo ""
    echo "✅ No hay K+ cristalográfico"
    echo "   PBS aproximado (Na+/Cl-) es adecuado"
    echo ""
fi

# ============================================================================
# OPCIÓN A: Pipeline drMD (Automático, Recomendado)
# ============================================================================
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║  OPCIÓN A: drMD Pipeline (Automático) ⭐ RECOMENDADO               ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""

echo "A.1: Instalar drMD (una vez)"
echo "────────────────────────────────────────────────────────────────────"
echo "cd ../../SciToolAgent/external/drMD"
echo "pip install -e ."
echo "cd -"
echo ""

echo "A.2: Ejecutar pipeline drMD"
echo "────────────────────────────────────────────────────────────────────"
echo "python run_drMD_pipeline.py"
echo ""
echo "Esto ejecutará automáticamente:"
echo "  - Minimización"
echo "  - Equilibración NVT"
echo "  - Equilibración NPT"
echo "  - Producción 10 ns"
echo "  - Clustering de conformaciones"
echo "  - Generación de sección de métodos"
echo ""

echo "A.3: Revisar resultados drMD"
echo "────────────────────────────────────────────────────────────────────"
echo "# Ver logs"
echo "cat drMD_output/00_drMD_logs/aftercare.log"
echo ""
echo "# Visualizar trayectoria"
echo "pymol drMD_output/5DRB/production/production_*.pdb"
echo ""
echo "# Ver clusters"
echo "ls drMD_output/5DRB/05_cluster_analysis/clusters/"
echo ""

# ============================================================================
# OPCIÓN B: Pipeline Manual (Control fino)
# ============================================================================
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║  OPCIÓN B: Pipeline Manual (Control fino)                          ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""

echo "B.1: Preparar sistema con PBS"
echo "────────────────────────────────────────────────────────────────────"
echo "python prepare_system.py"
echo ""
echo "Output:"
echo "  prepared_system/propka_results.pka      ← Estados de protonación"
echo "  prepared_system/system_solvated.pdb     ← Sistema completo"
echo "  prepared_system/equilibrated.pdb        ← Listo para producción"
echo ""

echo "B.2: Analizar estados de protonación (ProPKa)"
echo "────────────────────────────────────────────────────────────────────"
echo "python analyze_propka.py"
echo ""
echo "Revisar HIS (histidinas) con pKa cerca de 7.4:"
echo "  - HID: δ-protonada"
echo "  - HIE: ε-protonada"
echo "  - HIP: ambas protonadas (+1)"
echo ""

echo "B.3: Validar condiciones PBS"
echo "────────────────────────────────────────────────────────────────────"
echo "pytest tests/test_pbs_conditions.py -v"
echo ""
echo "Tests críticos:"
echo "  ✓ Fuerza iónica 150-170 mM"
echo "  ✓ Sistema neutralizado"
echo "  ✓ Box size >1 nm padding"
echo "  ✓ pH 7.4"
echo ""

echo "B.4: Generar ventanas de umbrella sampling"
echo "────────────────────────────────────────────────────────────────────"
echo "python generate_umbrella_windows.py"
echo ""
echo "Output:"
echo "  umbrella_windows/windows_config.csv"
echo "  umbrella_windows/window_XX/initial.pdb  (20 ventanas)"
echo ""

echo "B.5: Test local (una ventana, CPU)"
echo "────────────────────────────────────────────────────────────────────"
echo "python run_umbrella_window.py --window 0 --steps 10000 --platform CPU"
echo ""
echo "Esto corre 10,000 pasos (20 ps) para verificar que funciona"
echo ""

# ============================================================================
# OPCIÓN C: Despliegue HPC (Producción)
# ============================================================================
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║  OPCIÓN C: HPC Deployment (Producción)                             ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""

echo "C.1: Configurar credenciales HPC"
echo "────────────────────────────────────────────────────────────────────"
echo "Editar deploy_to_hpc.sh:"
echo "  HPC_USER='tu_usuario'"
echo "  HPC_HOST='cluster.universidad.edu'"
echo "  HPC_WORK_DIR='/home/tu_usuario/wnk_umbrella'"
echo ""

echo "C.2: Ejecutar despliegue automático"
echo "────────────────────────────────────────────────────────────────────"
echo "bash deploy_to_hpc.sh"
echo ""
echo "Seleccionar:"
echo "  [1] Deploy completo (transferir + preparar + tests + submit)"
echo ""

echo "C.3: Monitorear jobs en HPC"
echo "────────────────────────────────────────────────────────────────────"
echo "bash deploy_to_hpc.sh"
echo ""
echo "Seleccionar:"
echo "  [7] Monitorear jobs"
echo ""

echo "C.4: Descargar resultados"
echo "────────────────────────────────────────────────────────────────────"
echo "bash deploy_to_hpc.sh"
echo ""
echo "Seleccionar:"
echo "  [6] Descargar resultados"
echo ""

echo "C.5: Analizar PMF con MBAR"
echo "────────────────────────────────────────────────────────────────────"
echo "python analyze_umbrella_mbar.py"
echo ""
echo "Output:"
echo "  pmf_analysis/pmf_results.csv"
echo "  pmf_analysis/pmf.png"
echo "  pmf_analysis/analysis_combined.png"
echo ""

# ============================================================================
# VALIDACIÓN FINAL
# ============================================================================
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║  VALIDACIÓN FINAL                                                  ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""

echo "Tests completos"
echo "────────────────────────────────────────────────────────────────────"
echo "# Pre-HPC tests (sistema, windows, bias)"
echo "pytest tests/test_wnk_umbrella_setup.py -v"
echo ""
echo "# PBS tests (fuerza iónica, pH)"
echo "pytest tests/test_pbs_conditions.py -v"
echo ""

echo "Verificar PBS vs Estándar"
echo "────────────────────────────────────────────────────────────────────"
echo "Comparar PMF con:"
echo "  1. PBS (163 mM) - IMPLEMENTADO"
echo "  2. Estándar (150 mM) - Para control"
echo ""
echo "Si Δ(PMF) < 2 kJ/mol → PBS aproximado es válido ✅"
echo "Si Δ(PMF) > 2 kJ/mol → Considerar CHARMM36 para K+ real"
echo ""

# ============================================================================
# DOCUMENTACIÓN
# ============================================================================
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║  DOCUMENTACIÓN DISPONIBLE                                          ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""
echo "README.md                          - Guía general del pipeline"
echo "PBS_BUFFER_IMPLEMENTATION.md       - Detalles técnicos PBS"
echo "PBS_DRMD_INTEGRATION_SUMMARY.md    - Resumen ejecutivo"
echo "PROTONACION_GUIDE.md               - Estados de protonación HIS"
echo "PIPELINE_DIAGRAM.md                - Diagrama visual ASCII completo"
echo "COMPLETION_REPORT.txt              - Reporte de implementación"
echo ""

# ============================================================================
# TROUBLESHOOTING
# ============================================================================
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║  TROUBLESHOOTING                                                   ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""

echo "Error: 'K+ not in forcefield'"
echo "→ Normal. Ver PBS_BUFFER_IMPLEMENTATION.md"
echo "  PBS aproxima K+ como Na+ (fuerza iónica equivalente)"
echo ""

echo "Error: 'ProPKa not found'"
echo "→ pip install propka"
echo ""

echo "Error: 'drMD not found'"
echo "→ cd ../../SciToolAgent/external/drMD && pip install -e ."
echo ""

echo "Error: 'Ionic strength too low/high'"
echo "→ pytest tests/test_pbs_conditions.py::test_ionic_strength_pbs"
echo "  Debe estar entre 150-170 mM"
echo ""

echo "Pregunta: '¿Necesito K+ real para WNK1?'"
echo "→ Verificar: grep 'K ' 5DRB.pdb | grep HETATM"
echo "  Si hay K+ en sitio activo → Considerar CHARMM36"
echo "  Si no → PBS aproximado OK"
echo ""

# ============================================================================
# SIGUIENTE PASO
# ============================================================================
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║  ⏭️  SIGUIENTE PASO RECOMENDADO                                    ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
echo ""
echo "1. Verificar K+ en 5DRB.pdb (arriba)"
echo ""
echo "2. Elegir workflow:"
echo "   A) drMD automático → python run_drMD_pipeline.py"
echo "   B) Manual          → python prepare_system.py"
echo "   C) HPC             → bash deploy_to_hpc.sh"
echo ""
echo "3. Validar PBS:"
echo "   pytest tests/test_pbs_conditions.py -v"
echo ""
echo "4. Continuar con umbrella sampling:"
echo "   python generate_umbrella_windows.py"
echo ""
echo "5. Analizar resultados:"
echo "   python analyze_umbrella_mbar.py"
echo ""
echo "╔═══════════════════════════════════════════════════════════════════╗"
echo "║                    🎉 SISTEMA LISTO 🎉                             ║"
echo "╚═══════════════════════════════════════════════════════════════════╝"
