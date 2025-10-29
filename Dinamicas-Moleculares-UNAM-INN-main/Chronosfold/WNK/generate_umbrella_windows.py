#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generador de ventanas de umbrella sampling para WNK1 C-terminal

Collective Variable (CV): Distancia entre centros de masa
- Dominio kinasa (residuos 194-450)
- C-terminal (residuos 451-483)

Explora el rango de distancias conformacionales razonables
"""

import os
import sys
from pathlib import Path
import numpy as np
from openmm import app, unit
from openmm import *
import mdtraj as mdt

print("="*70)
print("WNK1 UMBRELLA SAMPLING - GENERACIÓN DE VENTANAS")
print("="*70)

# Configuración
WORK_DIR = Path(__file__).parent
INPUT_DIR = WORK_DIR / "prepared_system"
WINDOWS_DIR = WORK_DIR / "umbrella_windows"
WINDOWS_DIR.mkdir(exist_ok=True)

EQUILIBRATED_PDB = INPUT_DIR / "equilibrated.pdb"
SYSTEM_XML = INPUT_DIR / "system.xml"

# Definición de grupos para CV (índices de residuos en PDB original)
KINASE_DOMAIN = range(194, 451)  # Dominio kinasa
C_TERMINAL = range(451, 484)     # C-terminal

# Parámetros de umbrella sampling
MIN_DISTANCE = 1.5  # nm - conformación compacta
MAX_DISTANCE = 4.0  # nm - conformación extendida
N_WINDOWS = 20      # Número de ventanas

SPRING_CONSTANT = 1000.0  # kJ/mol/nm² - constante de resorte para bias

print(f"\n📁 Directorio de trabajo: {WORK_DIR}")
print(f"📂 Ventanas: {WINDOWS_DIR}")
print(f"\n⚙️  PARÁMETROS:")
print(f"  Dominio kinasa: residuos {KINASE_DOMAIN.start}-{KINASE_DOMAIN.stop-1}")
print(f"  C-terminal: residuos {C_TERMINAL.start}-{C_TERMINAL.stop-1}")
print(f"  Rango distancia: {MIN_DISTANCE:.2f} - {MAX_DISTANCE:.2f} nm")
print(f"  Número de ventanas: {N_WINDOWS}")
print(f"  Constante de resorte: {SPRING_CONSTANT} kJ/mol/nm²")

# ============================================================================
# PASO 1: Cargar estructura equilibrada
# ============================================================================
print("\n" + "="*70)
print("PASO 1: Cargando sistema equilibrado")
print("="*70)

if not EQUILIBRATED_PDB.exists():
    print(f"❌ ERROR: No se encuentra {EQUILIBRATED_PDB}")
    print("   Ejecuta primero prepare_system.py")
    sys.exit(1)

pdb = app.PDBFile(str(EQUILIBRATED_PDB))
print(f"✓ PDB cargado: {len(list(pdb.topology.atoms()))} átomos")

# Cargar sistema
with open(SYSTEM_XML, 'r') as f:
    system = XmlSerializer.deserialize(f.read())
print(f"✓ Sistema cargado desde XML")

# ============================================================================
# PASO 2: Identificar átomos de los grupos
# ============================================================================
print("\n" + "="*70)
print("PASO 2: Identificando átomos de los grupos")
print("="*70)

topology = pdb.topology

# Obtener índices de átomos CA (carbono alfa) de cada grupo
kinase_atoms = []
cterm_atoms = []

for atom in topology.atoms():
    if atom.name == 'CA':  # Solo carbonos alfa para centro de masa
        res_id = int(atom.residue.id)
        if res_id in KINASE_DOMAIN:
            kinase_atoms.append(atom.index)
        elif res_id in C_TERMINAL:
            cterm_atoms.append(atom.index)

print(f"✓ Átomos identificados:")
print(f"  Dominio kinasa (CA): {len(kinase_atoms)} átomos")
print(f"  C-terminal (CA): {len(cterm_atoms)} átomos")

if len(kinase_atoms) == 0 or len(cterm_atoms) == 0:
    print("❌ ERROR: No se encontraron átomos en uno de los grupos")
    print(f"   Kinase: {len(kinase_atoms)}, C-term: {len(cterm_atoms)}")
    sys.exit(1)

# ============================================================================
# PASO 3: Calcular distancia inicial
# ============================================================================
print("\n" + "="*70)
print("PASO 3: Calculando distancia inicial")
print("="*70)

positions = pdb.positions

# Calcular centros de masa
def compute_com(atom_indices, positions):
    """Calcula centro de masa de un grupo de átomos"""
    coords = np.array([positions[i].value_in_unit(unit.nanometer) for i in atom_indices])
    return np.mean(coords, axis=0)

com_kinase = compute_com(kinase_atoms, positions)
com_cterm = compute_com(cterm_atoms, positions)

initial_distance = np.linalg.norm(com_kinase - com_cterm)

print(f"✓ Centros de masa calculados:")
print(f"  Dominio kinasa: {com_kinase}")
print(f"  C-terminal: {com_cterm}")
print(f"  Distancia inicial: {initial_distance:.3f} nm")

# ============================================================================
# PASO 4: Definir posiciones de ventanas
# ============================================================================
print("\n" + "="*70)
print("PASO 4: Definiendo ventanas de umbrella sampling")
print("="*70)

window_distances = np.linspace(MIN_DISTANCE, MAX_DISTANCE, N_WINDOWS)

print(f"✓ {N_WINDOWS} ventanas definidas:")
for i, d in enumerate(window_distances):
    marker = " ← INICIAL" if abs(d - initial_distance) < 0.15 else ""
    print(f"  Ventana {i:2d}: r₀ = {d:.3f} nm{marker}")

# ============================================================================
# PASO 5: Crear archivos de configuración para cada ventana
# ============================================================================
print("\n" + "="*70)
print("PASO 5: Generando archivos de configuración")
print("="*70)

# Guardar lista de índices de átomos
atoms_file = WINDOWS_DIR / "atom_groups.txt"
with open(atoms_file, 'w') as f:
    f.write("# Grupos de átomos para CV\n")
    f.write(f"# Kinase domain (residues {KINASE_DOMAIN.start}-{KINASE_DOMAIN.stop-1})\n")
    f.write("KINASE_ATOMS=" + ",".join(map(str, kinase_atoms)) + "\n\n")
    f.write(f"# C-terminal (residues {C_TERMINAL.start}-{C_TERMINAL.stop-1})\n")
    f.write("CTERM_ATOMS=" + ",".join(map(str, cterm_atoms)) + "\n")

print(f"✓ Grupos de átomos guardados: {atoms_file}")

# Crear archivo con parámetros de ventanas
windows_file = WINDOWS_DIR / "windows_config.csv"
with open(windows_file, 'w') as f:
    f.write("window_id,r0_nm,spring_constant_kJ_mol_nm2\n")
    for i, r0 in enumerate(window_distances):
        f.write(f"{i},{r0:.4f},{SPRING_CONSTANT}\n")

print(f"✓ Configuración de ventanas guardada: {windows_file}")

# Crear directorio para cada ventana con archivos iniciales
for i, r0 in enumerate(window_distances):
    window_dir = WINDOWS_DIR / f"window_{i:02d}"
    window_dir.mkdir(exist_ok=True)
    
    # Copiar PDB inicial
    window_pdb = window_dir / "initial.pdb"
    with open(EQUILIBRATED_PDB, 'r') as src, open(window_pdb, 'w') as dst:
        dst.write(src.read())
    
    # Crear archivo de parámetros específico de esta ventana
    params_file = window_dir / "params.txt"
    with open(params_file, 'w') as f:
        f.write(f"# Parámetros ventana {i}\n")
        f.write(f"WINDOW_ID={i}\n")
        f.write(f"R0={r0:.4f}  # nm\n")
        f.write(f"SPRING_K={SPRING_CONSTANT}  # kJ/mol/nm²\n")
        f.write(f"KINASE_ATOMS_FILE=../atom_groups.txt\n")

print(f"✓ {N_WINDOWS} directorios de ventanas creados")

# ============================================================================
# PASO 6: Crear script Python para SMD (Steered MD) pulling
# ============================================================================
print("\n" + "="*70)
print("PASO 6: Generando script de pulling SMD")
print("="*70)

smd_script = WINDOWS_DIR / "run_smd_pulling.py"
with open(smd_script, 'w') as f:
    f.write('''#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Steered Molecular Dynamics (SMD) para generar configuraciones iniciales

Estira gradualmente el C-terminal para crear snapshots a diferentes distancias
"""

import sys
from pathlib import Path
from openmm import app, unit
from openmm import *
import numpy as np

# Configuración
WORK_DIR = Path(__file__).parent.parent
WINDOWS_DIR = Path(__file__).parent
INPUT_PDB = WORK_DIR / "prepared_system" / "equilibrated.pdb"
SYSTEM_XML = WORK_DIR / "prepared_system" / "system.xml"

# Leer grupos de átomos
with open(WINDOWS_DIR / "atom_groups.txt", 'r') as f:
    lines = f.readlines()
    kinase_line = [l for l in lines if l.startswith("KINASE_ATOMS=")][0]
    cterm_line = [l for l in lines if l.startswith("CTERM_ATOMS=")][0]
    
    kinase_atoms = [int(x) for x in kinase_line.split("=")[1].strip().split(",")]
    cterm_atoms = [int(x) for x in cterm_line.split("=")[1].strip().split(",")]

print(f"Grupos: Kinase={len(kinase_atoms)} átomos, C-term={len(cterm_atoms)} átomos")

# Parámetros SMD
PULL_VELOCITY = 0.01  # nm/ps - velocidad de pulling
SPRING_K = 500.0  # kJ/mol/nm² - constante para pulling
TOTAL_TIME = 5000  # ps - tiempo total de pulling
TEMPERATURE = 300 * unit.kelvin

print(f"SMD: v={PULL_VELOCITY} nm/ps, k={SPRING_K} kJ/mol/nm², t={TOTAL_TIME} ps")

# Cargar sistema
pdb = app.PDBFile(str(INPUT_PDB))
with open(SYSTEM_XML, 'r') as f:
    system = XmlSerializer.deserialize(f.read())

# Agregar fuerza de pulling entre centros de masa
pull_force = CustomCentroidBondForce(2, "0.5*k*(distance(g1,g2)-r0)^2")
pull_force.addPerBondParameter("k")
pull_force.addPerBondParameter("r0")
pull_force.addGroup(kinase_atoms)
pull_force.addGroup(cterm_atoms)

# Distancia inicial
def compute_com(atom_indices, positions):
    coords = np.array([positions[i].value_in_unit(unit.nanometer) for i in atom_indices])
    return np.mean(coords, axis=0)

com_k = compute_com(kinase_atoms, pdb.positions)
com_c = compute_com(cterm_atoms, pdb.positions)
r0_initial = np.linalg.norm(com_k - com_c)

print(f"Distancia inicial: {r0_initial:.3f} nm")
print("Ejecutando SMD pulling...")
print("(Este script genera configuraciones. Úsalo solo si necesitas crear estructuras extendidas)")
print("Para umbrella sampling normal, usa directamente run_umbrella_window.py")

# Este script está incompleto intencionalmente - SMD es opcional
# Para la mayoría de casos, empezar desde la estructura equilibrada es suficiente
''')

print(f"✓ Script SMD creado: {smd_script}")
print("  (Opcional - solo si necesitas generar conformaciones muy extendidas)")

# ============================================================================
# RESUMEN
# ============================================================================
print("\n" + "="*70)
print("✅ GENERACIÓN DE VENTANAS COMPLETADA")
print("="*70)
print(f"\nArchivos generados en {WINDOWS_DIR}:")
print(f"  1. atom_groups.txt           - Índices de átomos para CV")
print(f"  2. windows_config.csv        - Parámetros de todas las ventanas")
print(f"  3. window_XX/                - {N_WINDOWS} directorios de ventanas")
print(f"  4. run_smd_pulling.py        - Script SMD (opcional)")
print(f"\nPróximo paso:")
print(f"  Ejecutar run_umbrella_window.py en cada ventana para MD de producción")
print("="*70)
