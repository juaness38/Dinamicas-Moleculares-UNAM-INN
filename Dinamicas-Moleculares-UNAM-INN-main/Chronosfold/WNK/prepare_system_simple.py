#!/usr/bin/env python3
"""
Preparación SIMPLIFICADA del sistema WNK1 para umbrella sampling
Sin dependencias externas (solo OpenMM)
"""

import os
from pathlib import Path
from openmm import app, unit
from openmm import *
import numpy as np

print("="*70)
print("WNK1 - PREPARACIÓN SIMPLE (Solo OpenMM)")
print("="*70)

# Configuración
WORK_DIR = Path(__file__).parent
PDB_FILE = WORK_DIR / "5DRB.pdb"
OUTPUT_DIR = WORK_DIR / "system_prepared"
OUTPUT_DIR.mkdir(exist_ok=True)

# Parámetros
FORCEFIELD_FILES = ['amber14-all.xml', 'amber14/tip3pfb.xml']
BOX_PADDING = 1.0 * unit.nanometer
IONIC_STRENGTH = 0.15 * unit.molar
TEMPERATURE = 310 * unit.kelvin
PRESSURE = 1.0 * unit.bar

print(f"\nArchivo PDB: {PDB_FILE}")
print(f"Salida: {OUTPUT_DIR}")

# ============================================================================
# PASO 1: Cargar y limpiar estructura
# ============================================================================
print("\n" + "="*70)
print("PASO 1: Cargando estructura")
print("="*70)

pdb = app.PDBFile(str(PDB_FILE))
print(f"✓ Cargado: {len(pdb.positions)} átomos")

# Crear modeller
modeller = app.Modeller(pdb.topology, pdb.positions)

# Eliminar agua cristalográfica
modeller.deleteWater()

# Eliminar ligandos y heteroátomos (solo proteína)
standard_residues = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
                     'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
                     'THR', 'TRP', 'TYR', 'VAL', 'HID', 'HIE', 'HIP'}

non_standard = [res for res in modeller.topology.residues() 
                if res.name not in standard_residues]

if non_standard:
    print(f"Eliminando {len(non_standard)} residuos no-estándar")
    modeller.delete(non_standard)

print(f"✓ Proteína limpia: {len(list(modeller.topology.atoms()))} átomos")

# ============================================================================
# PASO 2: Agregar hidrógenos
# ============================================================================
print("\n" + "="*70)
print("PASO 2: Agregando hidrógenos (pH 7.4)")
print("="*70)

forcefield = app.ForceField(*FORCEFIELD_FILES)

# Agregar hidrógenos sin pH-dependent protonation (simplificado)
# Usa estados de protonación por defecto
try:
    modeller.addHydrogens(forcefield)
    print(f"✓ Hidrógenos agregados: {len(list(modeller.topology.atoms()))} átomos totales")
except Exception as e:
    print(f"⚠️  Problema agregando hidrógenos: {e}")
    print("Intentando con variante neutra...")
    # Si falla, guardar y recargar (OpenMM a veces tiene problemas con caps)
    temp_pdb = OUTPUT_DIR / "temp_protein.pdb"
    with open(temp_pdb, 'w') as f:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
    
    # Recargar
    pdb_temp = app.PDBFile(str(temp_pdb))
    modeller = app.Modeller(pdb_temp.topology, pdb_temp.positions)
    modeller.addHydrogens(forcefield, pH=7.0, variants=['neutral'])
    print(f"✓ Hidrógenos agregados (intento 2): {len(list(modeller.topology.atoms()))} átomos")

# ============================================================================
# PASO 3: Solvatar en PBS buffer
# ============================================================================
print("\n" + "="*70)
print("PASO 3: Solvatando en caja TIP3P + PBS buffer")
print("="*70)

# Crear caja de agua
modeller.addSolvent(
    forcefield,
    model='tip3p',
    padding=BOX_PADDING,
    ionicStrength=IONIC_STRENGTH
)

n_atoms = len(list(modeller.topology.atoms()))
n_residues = len(list(modeller.topology.residues()))
print(f"✓ Sistema solvatado:")
print(f"  Total átomos: {n_atoms}")
print(f"  Total residuos: {n_residues}")

# Guardar sistema solvatado
solvated_pdb = OUTPUT_DIR / "solvated_ionized.pdb"
with open(solvated_pdb, 'w') as f:
    app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"✓ Guardado: {solvated_pdb}")

# ============================================================================
# PASO 4: Crear sistema y minimizar
# ============================================================================
print("\n" + "="*70)
print("PASO 4: Minimización de energía")
print("="*70)

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*unit.nanometer,
    constraints=app.HBonds
)

integrator = LangevinMiddleIntegrator(
    TEMPERATURE,
    1/unit.picosecond,
    0.004*unit.picoseconds
)

platform = Platform.getPlatformByName('CPU')
properties = {'Threads': '4'}

simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
simulation.context.setPositions(modeller.positions)

# Minimizar
print("Minimizando energía...")
initial_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
print(f"  Energía inicial: {initial_energy}")

simulation.minimizeEnergy(maxIterations=1000)

final_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
print(f"  Energía final: {final_energy}")
print(f"  ΔE: {final_energy - initial_energy}")

# Guardar estructura minimizada
minimized_positions = simulation.context.getState(getPositions=True).getPositions()
minimized_pdb = OUTPUT_DIR / "minimized.pdb"
with open(minimized_pdb, 'w') as f:
    app.PDBFile.writeFile(modeller.topology, minimized_positions, f)
print(f"✓ Guardado: {minimized_pdb}")

# Guardar estado completo para usar en umbrella sampling
state_xml = OUTPUT_DIR / "minimized_state.xml"
simulation.saveState(str(state_xml))
print(f"✓ Estado guardado: {state_xml}")

# ============================================================================
# RESUMEN FINAL
# ============================================================================
print("\n" + "="*70)
print("✓ PREPARACIÓN COMPLETADA")
print("="*70)
print(f"\nArchivos generados en: {OUTPUT_DIR}/")
print(f"  - solvated_ionized.pdb    (sistema completo)")
print(f"  - minimized.pdb           (minimizado)")
print(f"  - minimized_state.xml     (estado OpenMM)")
print(f"\nSiguiente paso:")
print(f"  python generate_umbrella_windows.py")
print("="*70)
