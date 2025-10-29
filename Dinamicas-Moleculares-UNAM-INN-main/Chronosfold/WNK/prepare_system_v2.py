#!/usr/bin/env python3
"""
Preparación WNK1 con manejo explícito de terminales
"""

from pathlib import Path
from openmm import app, unit
from openmm import *

print("="*70)
print("WNK1 - PREPARACIÓN CON TERMINAL CAPS")
print("="*70)

WORK_DIR = Path(__file__).parent
PDB_FILE = WORK_DIR / "5DRB.pdb"
OUTPUT_DIR = WORK_DIR / "system_prepared"
OUTPUT_DIR.mkdir(exist_ok=True)

FORCEFIELD_FILES = ['amber14-all.xml', 'amber14/tip3pfb.xml']

print(f"\nCargando: {PDB_FILE}")

# Cargar PDB
pdb = app.PDBFile(str(PDB_FILE))
modeller = app.Modeller(pdb.topology, pdb.positions)

# Eliminar agua
modeller.deleteWater()

# Eliminar heteroátomos
standard_residues = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
                     'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
                     'THR', 'TRP', 'TYR', 'VAL'}

non_standard = [res for res in modeller.topology.residues() 
                if res.name not in standard_residues]
if non_standard:
    print(f"Eliminando {len(non_standard)} ligandos/heteroátomos")
    modeller.delete(non_standard)

print(f"✓ Proteína: {len(list(modeller.topology.residues()))} residuos")

# Guardar PDB limpio (sin hidrógenos aún)
clean_pdb = OUTPUT_DIR / "protein_clean.pdb"
with open(clean_pdb, 'w') as f:
    app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"✓ PDB limpio guardado: {clean_pdb}")

# Crear forcefield
forcefield = app.ForceField(*FORCEFIELD_FILES)

# Agregar hidrógenos - OpenMM manejará los caps automáticamente
print("\nAgregando hidrógenos...")
try:
    # Intentar con addMissingHydrogens (más robusto)
    modeller.addHydrogens(forcefield, pH=7.4)
    print(f"✓ Hidrógenos agregados: {len(list(modeller.topology.atoms()))} átomos totales")
except ValueError as e:
    if "terminal" in str(e).lower() or "bond" in str(e).lower():
        print(f"⚠️  Problema con terminales: {e}")
        print("Intentando método alternativo...")
        
        # Recargar sin modeller
        pdb_clean = app.PDBFile(str(clean_pdb))
        
        # Crear topology con caps explícitos
        from openmm.app import Modeller, Topology
        
        new_topology = Topology()
        new_positions = []
        
        # Copiar átomos pero marcar primera/última cadena
        for chain in pdb_clean.topology.chains():
            new_chain = new_topology.addChain(chain.id)
            residues = list(chain.residues())
            
            for i, residue in enumerate(residues):
                # Marcar si es N-terminal o C-terminal
                is_first = (i == 0)
                is_last = (i == len(residues) - 1)
                
                new_residue = new_topology.addResidue(residue.name, new_chain)
                
                for atom in residue.atoms():
                    new_topology.addAtom(atom.name, atom.element, new_residue)
                    new_positions.append(pdb_clean.positions[atom.index])
        
        new_modeller = Modeller(new_topology, new_positions)
        
        # Ahora agregar hidrógenos sin especificar pH (usa defaults)
        new_modeller.addHydrogens(forcefield)
        modeller = new_modeller
        print(f"✓ Hidrógenos agregados (método 2): {len(list(modeller.topology.atoms()))} átomos")

# Solvatar
print("\nSolvatando en TIP3P + iones (150 mM)...")
modeller.addSolvent(
    forcefield,
    model='tip3p',
    padding=1.0*unit.nanometer,
    ionicStrength=0.15*unit.molar
)

n_atoms = len(list(modeller.topology.atoms()))
print(f"✓ Sistema solvatado: {n_atoms} átomos totales")

# Guardar
solvated_pdb = OUTPUT_DIR / "solvated_ionized.pdb"
with open(solvated_pdb, 'w') as f:
    app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"✓ Guardado: {solvated_pdb}")

# Minimizar
print("\nMinimizando energía...")
system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*unit.nanometer,
    constraints=app.HBonds
)

integrator = LangevinMiddleIntegrator(310*unit.kelvin, 1/unit.picosecond, 0.004*unit.picoseconds)
platform = Platform.getPlatformByName('CPU')
simulation = app.Simulation(modeller.topology, system, integrator, platform, {'Threads': '4'})
simulation.context.setPositions(modeller.positions)

initial_E = simulation.context.getState(getEnergy=True).getPotentialEnergy()
print(f"  E inicial: {initial_E}")

simulation.minimizeEnergy(maxIterations=1000)

final_E = simulation.context.getState(getEnergy=True).getPotentialEnergy()
print(f"  E final: {final_E}")

# Guardar minimizado
minimized_positions = simulation.context.getState(getPositions=True).getPositions()
minimized_pdb = OUTPUT_DIR / "minimized.pdb"
with open(minimized_pdb, 'w') as f:
    app.PDBFile.writeFile(modeller.topology, minimized_positions, f)
print(f"✓ Guardado: {minimized_pdb}")

# Guardar estado
state_xml = OUTPUT_DIR / "minimized_state.xml"
simulation.saveState(str(state_xml))
print(f"✓ Estado: {state_xml}")

print("\n" + "="*70)
print("✓ PREPARACIÓN COMPLETADA")
print("="*70)
print(f"\nArchivos en: {OUTPUT_DIR}/")
print("  - protein_clean.pdb")
print("  - solvated_ionized.pdb")
print("  - minimized.pdb")
print("  - minimized_state.xml")
print("\nSiguiente: python generate_umbrella_windows.py")
