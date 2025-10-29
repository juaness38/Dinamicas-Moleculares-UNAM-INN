#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulación de producción de prueba - 100 ns
WNK1 kinase en PBS buffer
"""

import sys
import time
from pathlib import Path
from openmm import app, unit
from openmm import *
from datetime import datetime

print("="*70)
print("WNK1 - SIMULACIÓN DE PRODUCCIÓN (100 ns)")
print("="*70)

# Configuración
WORK_DIR = Path(__file__).parent
INPUT_PDB = WORK_DIR / "prepared_system" / "minimized.pdb"
OUTPUT_DIR = WORK_DIR / "production_test_100ns"
OUTPUT_DIR.mkdir(exist_ok=True)

# Parámetros de simulación
TEMPERATURE = 310 * unit.kelvin  # 37°C fisiológico
PRESSURE = 1.0 * unit.bar
TIMESTEP = 2.0 * unit.femtoseconds  # 2 fs con constraints
TOTAL_STEPS = 50_000_000  # 100 ns (50M pasos × 2 fs)
REPORT_INTERVAL = 5000  # Reportar cada 10 ps (5000 × 2 fs)
CHECKPOINT_INTERVAL = 500000  # Checkpoint cada 1 ns

# Calcular ns totales
total_ns = (TOTAL_STEPS * TIMESTEP).value_in_unit(unit.nanosecond)

print(f"\n📁 Directorio: {WORK_DIR}")
print(f"📄 Input: {INPUT_PDB}")
print(f"📂 Output: {OUTPUT_DIR}")
print(f"\n⏱️  Configuración:")
print(f"  • Temperatura: {TEMPERATURE}")
print(f"  • Presión: {PRESSURE}")
print(f"  • Timestep: {TIMESTEP}")
print(f"  • Pasos totales: {TOTAL_STEPS:,}")
print(f"  • Simulación: {total_ns:.1f} ns")
print(f"  • Reportes cada: {REPORT_INTERVAL * 2} fs = {REPORT_INTERVAL * 2 / 1000:.1f} ps")

# ============================================================================
# PASO 1: Cargar sistema minimizado
# ============================================================================
print("\n" + "="*70)
print("PASO 1: Cargando sistema minimizado")
print("="*70)

if not INPUT_PDB.exists():
    print(f"❌ ERROR: No se encuentra {INPUT_PDB}")
    sys.exit(1)

pdb = app.PDBFile(str(INPUT_PDB))
print(f"✓ Sistema cargado: {len(pdb.positions)} átomos")

# ============================================================================
# PASO 2: Crear sistema con forcefield
# ============================================================================
print("\n" + "="*70)
print("PASO 2: Creando sistema OpenMM")
print("="*70)

forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*unit.nanometer,
    constraints=app.HBonds,
    rigidWater=True,
    ewaldErrorTolerance=0.0005
)

# Agregar barostato para NPT
system.addForce(MonteCarloBarostat(PRESSURE, TEMPERATURE, 25))

print(f"✓ Sistema creado")
print(f"  • Método no-enlazado: PME")
print(f"  • Cutoff: 1.0 nm")
print(f"  • Constraints: H-bonds")
print(f"  • Barostato: Monte Carlo (cada 25 pasos)")

# ============================================================================
# PASO 3: Configurar integrador y plataforma
# ============================================================================
print("\n" + "="*70)
print("PASO 3: Configurando integrador")
print("="*70)

integrator = LangevinMiddleIntegrator(
    TEMPERATURE,
    1.0/unit.picosecond,  # Fricción
    TIMESTEP
)

# Usar CPU (Yoltla no tiene GPUs modernas)
platform = Platform.getPlatformByName('CPU')
properties = {'Threads': '4'}  # 4 threads como en el job

print(f"✓ Integrador: Langevin Middle")
print(f"  • Temperatura: {TEMPERATURE}")
print(f"  • Fricción: 1.0 ps⁻¹")
print(f"  • Timestep: {TIMESTEP}")
print(f"✓ Plataforma: CPU (4 threads)")

# ============================================================================
# PASO 4: Crear simulación
# ============================================================================
print("\n" + "="*70)
print("PASO 4: Inicializando simulación")
print("="*70)

simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)

# Asignar velocidades según temperatura
simulation.context.setVelocitiesToTemperature(TEMPERATURE)

# Energía inicial
state = simulation.context.getState(getEnergy=True)
print(f"✓ Simulación inicializada")
print(f"  • Energía potencial: {state.getPotentialEnergy()}")
print(f"  • Energía cinética: {state.getKineticEnergy()}")

# ============================================================================
# PASO 5: Configurar reportes
# ============================================================================
print("\n" + "="*70)
print("PASO 5: Configurando salidas")
print("="*70)

# DCD trajectory
dcd_file = str(OUTPUT_DIR / "production.dcd")
simulation.reporters.append(app.DCDReporter(dcd_file, REPORT_INTERVAL))
print(f"✓ Trayectoria DCD: {dcd_file}")

# Log de energía y estado
log_file = str(OUTPUT_DIR / "production.log")
simulation.reporters.append(
    app.StateDataReporter(
        log_file,
        REPORT_INTERVAL,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True,
        density=True,
        progress=True,
        remainingTime=True,
        speed=True,
        totalSteps=TOTAL_STEPS,
        separator='\t'
    )
)
print(f"✓ Log de energía: {log_file}")

# Log a stdout (consola) con progreso
simulation.reporters.append(
    app.StateDataReporter(
        sys.stdout,
        REPORT_INTERVAL * 10,  # Cada 100 ps a consola
        step=True,
        time=True,
        speed=True,
        progress=True,
        remainingTime=True,
        totalSteps=TOTAL_STEPS,
        separator='\t'
    )
)
print(f"✓ Progreso a consola cada 100 ps")

# Checkpoint
checkpoint_file = str(OUTPUT_DIR / "checkpoint.chk")
simulation.reporters.append(app.CheckpointReporter(checkpoint_file, CHECKPOINT_INTERVAL))
print(f"✓ Checkpoints cada 1 ns: {checkpoint_file}")

# ============================================================================
# PASO 6: Ejecutar simulación
# ============================================================================
print("\n" + "="*70)
print("PASO 6: INICIANDO SIMULACIÓN DE PRODUCCIÓN")
print("="*70)
print(f"\n🚀 Simulando {total_ns:.1f} ns ({TOTAL_STEPS:,} pasos)")
print(f"⏰ Inicio: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"\n{'Step':<12} {'Time (ps)':<12} {'Speed (ns/day)':<15} {'Progress':<10} {'ETA'}")
print("-"*70)

start_time = time.time()

try:
    simulation.step(TOTAL_STEPS)
except Exception as e:
    print(f"\n❌ ERROR durante la simulación: {e}")
    # Guardar estado actual
    emergency_state = str(OUTPUT_DIR / "emergency_state.xml")
    simulation.saveState(emergency_state)
    print(f"💾 Estado de emergencia guardado: {emergency_state}")
    sys.exit(1)

end_time = time.time()
elapsed_time = end_time - start_time

# ============================================================================
# PASO 7: Guardar estado final
# ============================================================================
print("\n" + "="*70)
print("PASO 7: Guardando estado final")
print("="*70)

final_state = str(OUTPUT_DIR / "final_state.xml")
simulation.saveState(final_state)
print(f"✓ Estado final guardado: {final_state}")

# Guardar PDB final
final_positions = simulation.context.getState(getPositions=True).getPositions()
final_pdb = str(OUTPUT_DIR / "final.pdb")
with open(final_pdb, 'w') as f:
    app.PDBFile.writeFile(simulation.topology, final_positions, f)
print(f"✓ Estructura final: {final_pdb}")

# ============================================================================
# RESUMEN
# ============================================================================
print("\n" + "="*70)
print("✅ SIMULACIÓN COMPLETADA")
print("="*70)

hours = elapsed_time / 3600
ns_per_day = (total_ns / elapsed_time) * 86400

print(f"\n⏱️  Estadísticas:")
print(f"  • Tiempo total: {hours:.2f} horas ({elapsed_time/60:.1f} minutos)")
print(f"  • ns simulados: {total_ns:.1f} ns")
print(f"  • Rendimiento: {ns_per_day:.2f} ns/day")
print(f"  • Velocidad promedio: {TOTAL_STEPS/elapsed_time:.1f} pasos/segundo")

print(f"\n📊 Archivos generados:")
print(f"  • {dcd_file}")
print(f"  • {log_file}")
print(f"  • {final_pdb}")
print(f"  • {final_state}")
print(f"  • {checkpoint_file}")

print(f"\n💡 Para analizar:")
print(f"  python analyze_trajectory.py {dcd_file}")
print("="*70)
