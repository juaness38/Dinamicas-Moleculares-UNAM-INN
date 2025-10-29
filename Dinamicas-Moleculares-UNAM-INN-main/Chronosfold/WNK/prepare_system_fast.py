#!/usr/bin/env python3
"""
Preparación RÁPIDA del sistema WNK1
Solo lo esencial: estructura + solvente + minimización ligera
La equilibración se hará en cada ventana de umbrella
"""

import sys
from pathlib import Path
from openmm import app, unit
from openmm import *

print("="*70)
print("WNK1 - PREPARACIÓN RÁPIDA (sin equilibración)")
print("="*70)

WORK_DIR = Path(__file__).parent
PDB_FILE = WORK_DIR / "5DRB.pdb"
OUTPUT_DIR = WORK_DIR / "system_prepared"
OUTPUT_DIR.mkdir(exist_ok=True)

print(f"\n[1/5] Cargando {PDB_FILE}...")
pdb = app.PDBFile(str(PDB_FILE))
modeller = app.Modeller(pdb.topology, pdb.positions)

# Limpiar
print("[2/5] Limpiando estructura...")
modeller.deleteWater()
standard = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
            'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
non_standard = [r for r in modeller.topology.residues() if r.name not in standard]
if non_standard:
    modeller.delete(non_standard)
print(f"  ✓ {len(list(modeller.topology.residues()))} residuos")

# Usar PDBFixer si está disponible
print("[3/5] Agregando hidrógenos y átomos faltantes...")
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

try:
    from pdbfixer import PDBFixer
    temp_pdb = OUTPUT_DIR / "temp.pdb"
    with open(temp_pdb, 'w') as f:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
    
    fixer = PDBFixer(str(temp_pdb))
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)
    
    modeller.topology = fixer.topology
    modeller.positions = fixer.positions
    temp_pdb.unlink()
    print("  ✓ PDBFixer: terminales corregidos")
except:
    modeller.addHydrogens(forcefield, pH=7.0)
    print("  ✓ Hidrógenos agregados (sin PDBFixer)")

print(f"  Total: {len(list(modeller.topology.atoms()))} átomos")

# Solvatar
print("[4/5] Solvatando (PBS 150 mM)...")
modeller.addSolvent(forcefield, model='tip3p', padding=1.0*unit.nanometer,
                    ionicStrength=0.15*unit.molar)

n_atoms = len(list(modeller.topology.atoms()))
n_waters = sum(1 for r in modeller.topology.residues() if r.name == 'HOH')
print(f"  ✓ {n_atoms} átomos totales ({n_waters} aguas)")

# Guardar sistema solvatado
output_pdb = OUTPUT_DIR / "solvated_ionized.pdb"
with open(output_pdb, 'w') as f:
    app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"  ✓ Guardado: {output_pdb}")

# Minimización RÁPIDA (solo 100 pasos)
print("[5/5] Minimización rápida (100 pasos)...")
system = forcefield.createSystem(modeller.topology, 
                                  nonbondedMethod=app.PME,
                                  nonbondedCutoff=1.0*unit.nanometer,
                                  constraints=app.HBonds)

integrator = LangevinMiddleIntegrator(310*unit.kelvin, 1/unit.picosecond, 
                                       0.002*unit.picoseconds)
platform = Platform.getPlatformByName('CPU')
simulation = app.Simulation(modeller.topology, system, integrator, platform,
                            {'Threads': '4'})
simulation.context.setPositions(modeller.positions)

E_init = simulation.context.getState(getEnergy=True).getPotentialEnergy()
simulation.minimizeEnergy(maxIterations=100)
E_final = simulation.context.getState(getEnergy=True).getPotentialEnergy()

print(f"  E inicial: {E_init}")
print(f"  E final: {E_final}")

# Guardar minimizado
positions = simulation.context.getState(getPositions=True).getPositions()
minimized_pdb = OUTPUT_DIR / "minimized.pdb"
with open(minimized_pdb, 'w') as f:
    app.PDBFile.writeFile(modeller.topology, positions, f)
print(f"  ✓ {minimized_pdb}")

# Guardar estado
state_xml = OUTPUT_DIR / "minimized_state.xml"
simulation.saveState(str(state_xml))
print(f"  ✓ {state_xml}")

print("\n" + "="*70)
print("✅ PREPARACIÓN COMPLETADA (modo rápido)")
print("="*70)
print(f"\nArchivos en {OUTPUT_DIR}/:")
print("  - solvated_ionized.pdb    (sistema completo)")
print("  - minimized.pdb           (minimizado)")
print("  - minimized_state.xml     (estado OpenMM)")
print("\n⚠️  NOTA: Equilibración NVT/NPT se hará en cada ventana")
print("="*70)
