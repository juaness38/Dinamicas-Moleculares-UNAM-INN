# Pipeline Completo: WNK1 Umbrella Sampling con ProPKa

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                    WNK1 C-TERMINAL UMBRELLA SAMPLING                         │
│                         Complete Workflow                                    │
└─────────────────────────────────────────────────────────────────────────────┘

INPUT
  │
  ├─► 5DRB.pdb (WNK1 crystal structure, rat, 1.65Å)
  │
  └─► environment.yml (dependencies: openmm, propka, pymbar, mdtraj, ...)

═══════════════════════════════════════════════════════════════════════════════
FASE 1: PREPARACIÓN DEL SISTEMA
═══════════════════════════════════════════════════════════════════════════════

┌──────────────────────┐
│ prepare_system.py    │
└──────────────────────┘
          │
          ├─► PASO 1: Cargar PDB
          │      └─► Topology: 2905 átomos (proteína + agua cristalográfica)
          │
          ├─► PASO 2: Limpiar estructura
          │      └─► Eliminar agua, heteroátomos, solo proteína
          │
          ├─► PASO 3: ProPKa Analysis ⭐ NUEVO
          │      │
          │      ├─► Ejecutar: propka.run.single(pdb, pH=7.0)
          │      │
          │      ├─► Output: propka_results.pka
          │      │      └─► pKa predicho para ASP, GLU, HIS, LYS, ARG, CYS, TYR
          │      │
          │      └─► Reporte: Estados de Histidinas (HID/HIE/HIP)
          │
          ├─► PASO 4: Agregar hidrógenos
          │      └─► modeller.addHydrogens(forcefield, pH=7.0)
          │           └─► Usa estados informados por ProPKa
          │
          ├─► PASO 5: Solvatar
          │      ├─► Agua TIP3P
          │      ├─► Padding 1.0 nm
          │      ├─► Iones Na+/Cl- (0.15 M)
          │      └─► ~50,000-80,000 átomos totales
          │
          ├─► PASO 6: Guardar sistema
          │      ├─► system_solvated.pdb
          │      ├─► system.xml (OpenMM serialized)
          │      └─► initial_state.xml
          │
          ├─► PASO 7: Minimización (1000 steps)
          │      └─► minimized.pdb
          │
          ├─► PASO 8: NVT equilibration (100 ps, 300 K)
          │      └─► nvt_equilibrated.pdb
          │
          └─► PASO 9: NPT equilibration (100 ps, 300 K, 1 bar)
                 └─► equilibrated.pdb ← READY FOR PRODUCTION
                 └─► equilibrated_state.xml

OUTPUT: prepared_system/
  ├─ propka_results.pka ⭐ NUEVO
  ├─ system_solvated.pdb
  ├─ system.xml
  ├─ minimized.pdb
  ├─ nvt_equilibrated.pdb
  ├─ equilibrated.pdb ← INPUT FOR NEXT STEP
  └─ equilibrated_state.xml

┌──────────────────────┐
│ analyze_propka.py    │ ⭐ NUEVO
└──────────────────────┘
          │
          ├─► Parse propka_results.pka
          ├─► Identificar HIS críticas (Catalytic loop, C-terminal)
          ├─► Reportar pKa perturbados (|ΔpKa| > 1.0)
          └─► Recomendaciones para validación

═══════════════════════════════════════════════════════════════════════════════
FASE 2: GENERACIÓN DE VENTANAS
═══════════════════════════════════════════════════════════════════════════════

┌────────────────────────────┐
│ generate_umbrella_windows.py│
└────────────────────────────┘
          │
          ├─► INPUT: prepared_system/equilibrated.pdb
          │
          ├─► Definir grupos para CV:
          │      ├─► Grupo 1: Dominio kinasa (residuos 194-450, CA atoms)
          │      └─► Grupo 2: C-terminal (residuos 451-483, CA atoms)
          │
          ├─► Calcular distancia inicial COM-COM
          │      └─► d_initial ≈ 2.5-3.0 nm (compacto)
          │
          ├─► Crear N=20 ventanas:
          │      ├─► r₀ = linspace(1.5, 4.0, 20) nm
          │      └─► k = 1000 kJ/mol/nm²
          │
          └─► Generar archivos

OUTPUT: umbrella_windows/
  ├─ atom_groups.txt (índices de átomos)
  ├─ windows_config.csv (parámetros de 20 ventanas)
  └─ window_XX/ (20 directorios)
       └─ initial.pdb
       └─ params.txt

═══════════════════════════════════════════════════════════════════════════════
FASE 3: SIMULACIONES MD (HPC)
═══════════════════════════════════════════════════════════════════════════════

OPCIÓN A: Local (prueba)
┌──────────────────────────┐
│ run_umbrella_window.py   │
└──────────────────────────┘
     │
     └─► python run_umbrella_window.py --window 0 --steps 50000 --platform CPU

OPCIÓN B: HPC (producción) ⭐ RECOMENDADO
┌──────────────────────────┐
│ deploy_to_hpc.sh         │
└──────────────────────────┘
          │
          ├─► MENÚ INTERACTIVO:
          │   [1] Deploy completo (transferir + preparar + tests + submit)
          │   [2] Solo transferir
          │   [3] Solo preparación
          │   [4] Solo tests
          │   [5] Solo submit jobs
          │   [6] Descargar resultados
          │   [7] Monitorear jobs
          │
          ├─► Configurar:
          │      ├─ HPC_USER="tu_usuario"
          │      ├─ HPC_HOST="cluster.edu"
          │      └─ HPC_WORK_DIR="/home/user/wnk"
          │
          ├─► Transferir archivos (scp/rsync)
          │
          ├─► Ejecutar en HPC:
          │      ├─► prepare_system.py
          │      ├─► generate_umbrella_windows.py
          │      └─► pytest test_wnk_umbrella_setup.py
          │
          └─► Submit SLURM job array

┌──────────────────────────┐
│ submit_umbrella_hpc.sh   │ (SLURM batch script)
└──────────────────────────┘
          │
          ├─► #SBATCH --array=0-19 (20 jobs en paralelo)
          ├─► #SBATCH --gres=gpu:1
          ├─► #SBATCH --time=48:00:00
          │
          └─► Para cada ventana:
                 │
                 └─► run_umbrella_window.py
                       │
                       ├─► Cargar sistema + bias armónico
                       │      U_bias = 0.5 * k * (r - r₀)²
                       │
                       ├─► MD producción (5M steps = 10 ns)
                       │      ├─ Langevin integrator (300 K)
                       │      ├─ dt = 2 fs
                       │      └─ CUDA/GPU platform
                       │
                       └─► Output:
                            ├─ trajectory.dcd
                            ├─ cv_values.dat (r, E_bias vs time)
                            ├─ production.log
                            └─ final.pdb

OUTPUT: umbrella_windows/window_XX/ (para cada ventana)
  ├─ trajectory.dcd (trayectoria completa)
  ├─ cv_values.dat (CV y bias energy vs tiempo)
  ├─ production.log (termodinámica)
  └─ final.pdb (estructura final)

═══════════════════════════════════════════════════════════════════════════════
FASE 4: ANÁLISIS MBAR
═══════════════════════════════════════════════════════════════════════════════

┌──────────────────────────┐
│ analyze_umbrella_mbar.py │
└──────────────────────────┘
          │
          ├─► Cargar cv_values.dat de 20 ventanas
          │      └─► Concatenar: ~1M-2M samples total
          │
          ├─► Calcular matriz de bias u_kn:
          │      u_kn[k,n] = β * 0.5 * k * (r_n - r₀_k)²
          │      donde k=ventana, n=sample
          │
          ├─► Ejecutar MBAR:
          │      ├─► pymbar.MBAR(u_kn, N_k)
          │      ├─► Converge free energies F_k
          │      └─► Iterative solver (~200 iterations)
          │
          ├─► Calcular PMF en bins:
          │      ├─► Bins: 50 puntos en rango 1.5-4.0 nm
          │      ├─► PMF(r) = -kT * ln P(r)
          │      └─► Normalizar: min(PMF) = 0
          │
          └─► Plots:
                 ├─► pmf.png (PMF vs distancia)
                 └─► analysis_combined.png (histogramas + PMF)

OUTPUT: pmf_analysis/
  ├─ pmf_results.csv (r, PMF, error)
  ├─ pmf.png
  └─ analysis_combined.png

═══════════════════════════════════════════════════════════════════════════════
VALIDACIÓN Y TESTS
═══════════════════════════════════════════════════════════════════════════════

┌──────────────────────────────┐
│ test_wnk_umbrella_setup.py   │
└──────────────────────────────┘
          │
          ├─► TestSystemPreparation
          │      ├─ test_pdb_file_exists
          │      ├─ test_pdb_is_valid
          │      └─ test_forcefield_loads
          │
          ├─► TestPreparedSystem
          │      ├─ test_prepared_system_has_water
          │      ├─ test_prepared_system_has_ions
          │      ├─ test_system_has_periodic_box
          │      └─ test_system_can_create_context
          │
          ├─► TestUmbrellaWindows
          │      ├─ test_windows_directory_exists
          │      ├─ test_atom_groups_are_valid
          │      └─ test_windows_config_is_valid
          │
          ├─► TestBiasForce
          │      ├─ test_can_add_centroid_force
          │      └─ test_bias_force_gives_reasonable_energy
          │
          └─► TestShortSimulation
                 └─ test_can_run_100_steps

Ejecutar: pytest test_wnk_umbrella_setup.py -v

═══════════════════════════════════════════════════════════════════════════════
DOCUMENTACIÓN
═══════════════════════════════════════════════════════════════════════════════

📄 README.md
   └─► Guía completa de uso, parámetros, troubleshooting

📄 PROTONACION_GUIDE.md ⭐ NUEVO
   └─► Estados de protonación, HIS críticas, workflow de validación

📄 PROPKA_INTEGRATION_SUMMARY.md ⭐ NUEVO
   └─► Resumen de la integración de ProPKa, por qué es importante

📄 PIPELINE_DIAGRAM.md (este archivo)
   └─► Diagrama visual del workflow completo

═══════════════════════════════════════════════════════════════════════════════
RESULTADOS ESPERADOS
═══════════════════════════════════════════════════════════════════════════════

1. PMF del C-terminal de WNK1:
   
   ┌────────────────────────────────────────┐
   │  PMF                                    │
   │   ↑                                     │
   │ 15│         ╱‾‾╲  ← Barrera            │
   │   │        ╱    ╲                      │
   │ 10│       ╱      ╲                     │
   │   │      │        ╲___                 │
   │  5│     ╱              ‾‾╲             │
   │   │    ╱                   ╲           │
   │  0├───╯─────────────────────╲─────────►│
   │   1.5   2.0   2.5   3.0   3.5   4.0   │
   │          Distancia C-term (nm)         │
   └────────────────────────────────────────┘
   
   Interpretación:
   - Mínimo global: Estado compacto (r ≈ 2.0 nm)
   - Barrera: Transición conformacional (r ≈ 2.8 nm)
   - Estado extendido: Menos favorable (r > 3.5 nm)

2. Métricas:
   - Barrera energética: 10-20 kJ/mol (esperado para cambios conformacionales)
   - Precisión MBAR: <0.5 kcal/mol (<2 kJ/mol)
   - Overlap entre ventanas: >50% para convergencia

3. Insights biológicos:
   - ¿C-terminal preferentemente compacto o extendido?
   - ¿Barrera accesible a kT (2.5 kJ/mol a 300K)?
   - ¿Múltiples mínimos? → Estados funcionales distintos

═══════════════════════════════════════════════════════════════════════════════
CRONOLOGÍA DE EJECUCIÓN
═══════════════════════════════════════════════════════════════════════════════

Local (testing):
  prepare_system.py         : 10-30 min
  generate_umbrella_windows : <1 min
  run_umbrella_window (1)   : 10-30 min (100k steps, CPU)
  analyze_propka            : <1 min
  TOTAL                     : ~1 hora

HPC (production):
  deploy_to_hpc.sh (setup)  : 5-10 min
  prepare_system.py (HPC)   : 10-30 min
  submit_umbrella_hpc.sh    : Envío inmediato
  
  Simulaciones (20 ventanas en paralelo):
    - Con GPU: 2-4 horas cada ventana
    - Con CPU: 8-12 horas cada ventana
  
  analyze_umbrella_mbar.py  : 5-10 min
  
  TOTAL (GPU)               : ~4-6 horas
  TOTAL (CPU)               : ~12-16 horas

═══════════════════════════════════════════════════════════════════════════════
ARCHIVOS FINALES
═══════════════════════════════════════════════════════════════════════════════

WNK/
├── Scripts de preparación:
│   ├── prepare_system.py ⭐ Integra ProPKa
│   ├── generate_umbrella_windows.py
│   └── analyze_propka.py ⭐ NUEVO
│
├── Scripts de simulación:
│   ├── run_umbrella_window.py
│   ├── submit_umbrella_hpc.sh
│   └── deploy_to_hpc.sh
│
├── Scripts de análisis:
│   └── analyze_umbrella_mbar.py
│
├── Datos:
│   ├── 5DRB.pdb
│   ├── prepared_system/
│   │   ├── propka_results.pka ⭐ NUEVO
│   │   ├── equilibrated.pdb
│   │   └── equilibrated_state.xml
│   ├── umbrella_windows/
│   │   ├── windows_config.csv
│   │   └── window_XX/
│   │       ├── trajectory.dcd
│   │       └── cv_values.dat
│   └── pmf_analysis/
│       ├── pmf_results.csv
│       └── pmf.png
│
└── Documentación:
    ├── README.md
    ├── PROTONACION_GUIDE.md ⭐ NUEVO
    ├── PROPKA_INTEGRATION_SUMMARY.md ⭐ NUEVO
    └── PIPELINE_DIAGRAM.md (este archivo)

═══════════════════════════════════════════════════════════════════════════════
NOTAS IMPORTANTES
═══════════════════════════════════════════════════════════════════════════════

⚠️  SIEMPRE revisar propka_results.pka antes de simulaciones largas

⚠️  Histidinas (HIS) son críticas - verificar estados de protonación

⚠️  Si HIS en C-terminal tiene pKa cerca de 7.0 → Considerar múltiples estados

✓  ProPKa mejora precisión de protonación en residuos buried/interface

✓  Pipeline completo es reproducible y automatizable

✓  HPC deployment automatizado con deploy_to_hpc.sh

═══════════════════════════════════════════════════════════════════════════════
```
