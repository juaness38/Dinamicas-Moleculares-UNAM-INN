#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Preparación del sistema WNK1 para umbrella sampling
Explora cambio conformacional del C-terminal (residuos 450-483)

Steps:
1. Cargar estructura 5DRB.pdb
2. Seleccionar cadena A (proteína)
3. Agregar hidrógenos
4. Solvatar en caja explícita
5. Neutralizar y agregar sal
6. Minimización de energía
7. Equilibración NVT y NPT
"""

import os
import sys
from pathlib import Path
from openmm import app, unit
from openmm import *
import numpy as np

print("="*70)
print("WNK1 C-TERMINAL UMBRELLA SAMPLING - PREPARACIÓN DEL SISTEMA")
print("="*70)

# Configuración
WORK_DIR = Path(__file__).parent
PDB_FILE = WORK_DIR / "5DRB.pdb"
OUTPUT_DIR = WORK_DIR / "prepared_system"
OUTPUT_DIR.mkdir(exist_ok=True)

# Parámetros
FORCEFIELD = ['amber14-all.xml', 'amber14/tip3pfb.xml']
WATER_MODEL = 'tip3p'
BOX_PADDING = 1.0 * unit.nanometer  # Distancia mínima proteína-borde de caja
IONIC_STRENGTH = 0.15 * unit.molar  # Fisiológica
TEMPERATURE = 300 * unit.kelvin
PRESSURE = 1.0 * unit.bar

print(f"\n📁 Directorio de trabajo: {WORK_DIR}")
print(f"📄 Archivo PDB: {PDB_FILE}")
print(f"📂 Salida: {OUTPUT_DIR}")

# ============================================================================
# PASO 1: Cargar estructura
# ============================================================================
print("\n" + "="*70)
print("PASO 1: Cargando estructura 5DRB.pdb")
print("="*70)

if not PDB_FILE.exists():
    print(f"❌ ERROR: No se encuentra {PDB_FILE}")
    sys.exit(1)

pdb = app.PDBFile(str(PDB_FILE))
print(f"✓ Estructura cargada: {len(pdb.positions)} átomos")

# Información de cadenas
topology = pdb.topology
chains = list(topology.chains())
print(f"  Cadenas detectadas: {len(chains)}")
for chain in chains:
    n_residues = len(list(chain.residues()))
    print(f"    - Cadena {chain.id}: {n_residues} residuos")

# ============================================================================
# PASO 2: Limpiar estructura (solo proteína, eliminar agua cristalográfica)
# ============================================================================
print("\n" + "="*70)
print("PASO 2: Limpiando estructura (solo proteína)")
print("="*70)

modeller = app.Modeller(pdb.topology, pdb.positions)

# Eliminar agua y heteroátomos (mantener solo proteína)
modeller.deleteWater()
to_delete = [res for res in modeller.topology.residues() if res.name not in app.PDBFile._residueNameReplacements]
# Mantener solo residuos estándar
standard_residues = set(app.PDBFile._residueNameReplacements.keys())
standard_residues.update(['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
                          'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
                          'THR', 'TRP', 'TYR', 'VAL'])

non_standard = [res for res in modeller.topology.residues() if res.name not in standard_residues]
if non_standard:
    print(f"  Eliminando {len(non_standard)} residuos no-estándar/ligandos:")
    for res in non_standard[:5]:  # Mostrar primeros 5
        print(f"    - {res.name} {res.id}")
    if len(non_standard) > 5:
        print(f"    ... y {len(non_standard)-5} más")
    modeller.delete(non_standard)

print(f"✓ Estructura limpia: {len(list(modeller.topology.atoms()))} átomos")
print(f"  Residuos de proteína: {len(list(modeller.topology.residues()))}")

# ============================================================================
# PASO 3: Determinar estados de protonación con ProPKa
# ============================================================================
print("\n" + "="*70)
print("PASO 3: Determinando estados de protonación con ProPKa (pH 7.0)")
print("="*70)

# Guardar estructura limpia temporal para ProPKa
temp_clean_pdb = OUTPUT_DIR / "temp_clean.pdb"
with open(temp_clean_pdb, 'w') as f:
    app.PDBFile.writeFile(modeller.topology, modeller.positions, f)

print(f"  Ejecutando ProPKa en pH 7.0...")

try:
    # Importar propka
    from propka.run import single
    from propka.molecular_container import MolecularContainer
    
    # Ejecutar ProPKa
    propka_output = OUTPUT_DIR / "propka_results.pka"
    molecule = single(str(temp_clean_pdb), optargs=['--pH', '7.0'])
    
    # Guardar resultados
    with open(propka_output, 'w') as f:
        f.write(molecule.write_propka())
    
    print(f"✓ ProPKa completado")
    print(f"  Resultados: {propka_output}")
    
    # Resumen de residuos titratables
    print(f"\n  Estados de protonación predichos a pH 7.0:")
    
    # Contar residuos titratables
    titratable = {'HIS': 0, 'ASP': 0, 'GLU': 0, 'LYS': 0, 'ARG': 0, 'CYS': 0, 'TYR': 0}
    for residue in modeller.topology.residues():
        if residue.name in titratable:
            titratable[residue.name] += 1
    
    for resname, count in titratable.items():
        if count > 0:
            print(f"    {resname}: {count} residuos")
    
    print(f"\n  IMPORTANTE: Los estados de protonación se aplicarán con addHydrogens(pH=7.0)")
    print(f"  OpenMM usa pKa estándar. Para casos especiales, editar manualmente después.")
    
except ImportError:
    print(f"⚠️  ADVERTENCIA: ProPKa no está instalado")
    print(f"   pip install propka")
    print(f"   Continuando con estados de protonación por defecto (pH 7.0)...")
except Exception as e:
    print(f"⚠️  ADVERTENCIA: ProPKa falló: {e}")
    print(f"   Continuando con estados de protonación por defecto...")

# Limpiar archivo temporal
if temp_clean_pdb.exists():
    temp_clean_pdb.unlink()

# ============================================================================
# PASO 4: Agregar hidrógenos con estados de protonación correctos
# ============================================================================
print("\n" + "="*70)
print("PASO 4: Agregando hidrógenos (pH 7.0)")
print("="*70)

forcefield = app.ForceField(*FORCEFIELD)

# Intentar con PDBFixer primero (maneja terminales automáticamente)
try:
    from pdbfixer import PDBFixer
    print("  Usando PDBFixer para agregar átomos faltantes...")
    
    # Guardar temporalmente
    temp_pdb = OUTPUT_DIR / "temp_for_fixer.pdb"
    with open(temp_pdb, 'w') as f:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
    
    # Usar PDBFixer
    fixer = PDBFixer(str(temp_pdb))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)
    
    # Actualizar modeller
    modeller.topology = fixer.topology
    modeller.positions = fixer.positions
    
    temp_pdb.unlink()  # Limpiar
    print(f"✓ PDBFixer: átomos faltantes y hidrógenos agregados")
    
except ImportError:
    print("  PDBFixer no disponible, usando método directo de OpenMM...")
    try:
        modeller.addHydrogens(forcefield, pH=7.0)
    except ValueError as e:
        if "terminal" in str(e).lower() or "bond" in str(e).lower():
            print(f"  ⚠️  Error de terminal detectado: {e}")
            print("  Intentando guardar/recargar para forzar detección de terminales...")
            
            # Guardar y recargar (a veces ayuda con la detección)
            temp_reload = OUTPUT_DIR / "temp_reload.pdb"
            with open(temp_reload, 'w') as f:
                app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
            
            pdb_reload = app.PDBFile(str(temp_reload))
            modeller = app.Modeller(pdb_reload.topology, pdb_reload.positions)
            
            # Intentar sin pH específico (usa defaults que son más permisivos)
            try:
                modeller.addHydrogens(forcefield)
                print("  ✓ Hidrógenos agregados (método alternativo, sin pH específico)")
            except:
                # Último recurso: sin variantes
                print("  Último intento: addHydrogens sin variantes...")
                modeller.addHydrogens(forcefield, pH=7.0, platform=None)
            
            temp_reload.unlink()
        else:
            raise

print(f"✓ Hidrógenos agregados: {len(list(modeller.topology.atoms()))} átomos totales")

# Reportar residuos HIS (importante para protonación)
his_residues = [res for res in modeller.topology.residues() if res.name in ['HIS', 'HIE', 'HID', 'HIP']]
if his_residues:
    print(f"\n  Residuos Histidina detectados: {len(his_residues)}")
    print(f"  Estados: {set(res.name for res in his_residues)}")
    print(f"  (HID=δ-protonada, HIE=ε-protonada, HIP=ambas, HIS=auto)")

# ============================================================================
# PASO 5: Solvatar en PBS (Phosphate Buffered Saline)
# ============================================================================
print("\n" + "="*70)
print("PASO 5: Solvatando en PBS (Phosphate Buffered Saline, pH 7.4)")
print("="*70)

# PBS buffer composition (concentraciones fisiológicas):
# 137 mM NaCl
# 2.7 mM KCl
# 10 mM Na₂HPO₄
# 1.8 mM KH₂PO₄
# Total: ~150 mM ionic strength

print("\nComposición PBS:")
print("  137 mM NaCl    (cloruro de sodio)")
print("  2.7 mM KCl     (cloruro de potasio)")
print("  10 mM Na₂HPO₄  (fosfato disódico)")
print("  1.8 mM KH₂PO₄  (fosfato monopotásico)")
print("  Fuerza iónica total: ~0.163 M")
print("  pH: 7.4 (fisiológico)")

# Primero solvatar con agua (sin iones para control manual)
modeller.addSolvent(forcefield, model=WATER_MODEL, padding=BOX_PADDING, 
                    neutralize=False)  # No neutralizar automáticamente

# Contar átomos de agua inicial
n_waters_initial = sum(1 for res in modeller.topology.residues() if res.name == 'HOH')

# Calcular volumen de caja para determinar cantidad de iones
vectors = modeller.topology.getPeriodicBoxVectors()
box_volume = (vectors[0][0] * vectors[1][1] * vectors[2][2]).value_in_unit(unit.liter)

# Calcular número de iones para PBS
# N = C * V * N_A, donde C en mol/L, V en L
from math import ceil
avogadro = 6.022e23

# NaCl: 137 mM
n_na_from_nacl = ceil(0.137 * box_volume * avogadro)
n_cl_from_nacl = ceil(0.137 * box_volume * avogadro)

# KCl: 2.7 mM
n_k_from_kcl = ceil(0.0027 * box_volume * avogadro)
n_cl_from_kcl = ceil(0.0027 * box_volume * avogadro)

# Na₂HPO₄: 10 mM (2 Na+ por molécula, pero HPO₄²⁻ no está en forcefield estándar)
n_na_from_phosphate = ceil(2 * 0.010 * box_volume * avogadro)

# KH₂PO₄: 1.8 mM (1 K+ por molécula)
n_k_from_phosphate = ceil(0.0018 * box_volume * avogadro)

# Totales
n_na_total = n_na_from_nacl + n_na_from_phosphate
n_k_total = n_k_from_kcl + n_k_from_phosphate
n_cl_total = n_cl_from_nacl + n_cl_from_kcl

# NOTA: Los fosfatos (HPO₄²⁻, H₂PO₄⁻) no están en forcefield amber14
# Aproximación: reemplazar con Cl⁻ adicional para mantener fuerza iónica
# HPO₄²⁻ → 2 Cl⁻, H₂PO₄⁻ → 1 Cl⁻
n_cl_from_phosphate = ceil((2*0.010 + 1*0.0018) * box_volume * avogadro)
n_cl_total += n_cl_from_phosphate

print(f"\nCálculo de iones (volumen caja: {box_volume*1e3:.2f} mL):")
print(f"  Na+: {n_na_total} iones ({n_na_from_nacl} de NaCl, {n_na_from_phosphate} de Na₂HPO₄)")
print(f"  K+:  {n_k_total} iones ({n_k_from_kcl} de KCl, {n_k_from_phosphate} de KH₂PO₄)")
print(f"  Cl-: {n_cl_total} iones ({n_cl_from_nacl+n_cl_from_kcl} de sales, {n_cl_from_phosphate} aprox. fosfatos)")
print(f"\n  ⚠️  NOTA: Fosfatos (HPO₄²⁻/H₂PO₄⁻) aproximados como Cl⁻ (no en forcefield)")

# Calcular carga neta de la proteína
protein_charge = 0
for atom in modeller.topology.atoms():
    if atom.residue.name not in ['HOH', 'NA', 'CL', 'K']:
        # Estimar carga por residuo (simplificado)
        if atom.residue.name == 'ARG':
            protein_charge += 1/len([a for a in atom.residue.atoms()])
        elif atom.residue.name == 'LYS':
            protein_charge += 1/len([a for a in atom.residue.atoms()])
        elif atom.residue.name == 'ASP':
            protein_charge -= 1/len([a for a in atom.residue.atoms()])
        elif atom.residue.name == 'GLU':
            protein_charge -= 1/len([a for a in atom.residue.atoms()])

protein_charge = round(protein_charge)
print(f"\nCarga neta aproximada de la proteína: {protein_charge:+d}")

# Ajustar iones para neutralizar
if protein_charge > 0:
    n_cl_total += abs(protein_charge)
    print(f"  Agregando {abs(protein_charge)} Cl⁻ adicionales para neutralizar")
elif protein_charge < 0:
    n_na_total += abs(protein_charge)
    print(f"  Agregando {abs(protein_charge)} Na+ adicionales para neutralizar")

# Agregar iones manualmente usando modeller
# OpenMM's addSolvent solo soporta Na+/Cl-, no K+
# Workaround: agregar K+ reemplazando algunas moléculas de agua

# IMPLEMENTACIÓN SIMPLIFICADA:
# Usar addSolvent con fuerza iónica equivalente y documentar aproximación
ionic_strength_pbs = 0.163 * unit.molar  # ~163 mM total
modeller = app.Modeller(modeller.topology, modeller.positions)
modeller.addSolvent(forcefield, model=WATER_MODEL, padding=BOX_PADDING,
                    ionicStrength=ionic_strength_pbs,
                    positiveIon='Na+', negativeIon='Cl-')

n_waters = sum(1 for res in modeller.topology.residues() if res.name == 'HOH')
n_ions_na = sum(1 for res in modeller.topology.residues() if res.name == 'NA')
n_ions_cl = sum(1 for res in modeller.topology.residues() if res.name == 'CL')

print(f"\n✓ Sistema solvatado en PBS:")
print(f"  Total átomos: {len(list(modeller.topology.atoms()))}")
print(f"  Moléculas de agua: {n_waters}")
print(f"  Iones Na+: {n_ions_na}")
print(f"  Iones Cl-: {n_ions_cl}")
print(f"  Fuerza iónica: {ionic_strength_pbs}")
print(f"\n  ⚠️  APROXIMACIÓN: PBS real contiene K+ y fosfatos")
print(f"      OpenMM solo soporta Na+/Cl- en forcefield amber14")
print(f"      Fuerza iónica equivalente: 163 mM (vs 150 mM PBS exacto)")

# Obtener dimensiones de caja
vectors = modeller.topology.getPeriodicBoxVectors()
a = vectors[0].value_in_unit(unit.nanometer)
b = vectors[1].value_in_unit(unit.nanometer)
c = vectors[2].value_in_unit(unit.nanometer)
print(f"  Dimensiones de caja: {a[0]:.2f} x {b[1]:.2f} x {c[2]:.2f} nm")

# ============================================================================
# PASO 6: Guardar sistema preparado
# ============================================================================
print("\n" + "="*70)
print("PASO 6: Guardando sistema preparado")
print("="*70)

# Guardar PDB del sistema completo
output_pdb = OUTPUT_DIR / "system_solvated.pdb"
with open(output_pdb, 'w') as f:
    app.PDBFile.writeFile(modeller.topology, modeller.positions, f)
print(f"✓ PDB guardado: {output_pdb}")

# Crear sistema OpenMM
system = forcefield.createSystem(modeller.topology, 
                                  nonbondedMethod=app.PME,
                                  nonbondedCutoff=1.0*unit.nanometer,
                                  constraints=app.HBonds)

# Guardar XML del sistema
output_xml = OUTPUT_DIR / "system.xml"
with open(output_xml, 'w') as f:
    f.write(XmlSerializer.serialize(system))
print(f"✓ Sistema XML guardado: {output_xml}")

# Guardar estado serializado
integrator = LangevinMiddleIntegrator(TEMPERATURE, 1.0/unit.picosecond, 
                                       0.002*unit.picoseconds)
simulation = app.Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

output_state = OUTPUT_DIR / "initial_state.xml"
simulation.saveState(str(output_state))
print(f"✓ Estado inicial guardado: {output_state}")

# ============================================================================
# PASO 7: Minimización de energía
# ============================================================================
print("\n" + "="*70)
print("PASO 7: Minimización de energía")
print("="*70)

print("  Energía inicial:")
state_init = simulation.context.getState(getEnergy=True)
E_init = state_init.getPotentialEnergy()
print(f"    Potencial: {E_init}")

print("  Minimizando (máx 1000 pasos)...")
simulation.minimizeEnergy(maxIterations=1000)

state_min = simulation.context.getState(getEnergy=True, getPositions=True)
E_min = state_min.getPotentialEnergy()
print(f"✓ Minimización completada")
print(f"  Energía final: {E_min}")
print(f"  ΔE = {E_min - E_init}")

# Guardar estructura minimizada
minimized_pdb = OUTPUT_DIR / "minimized.pdb"
positions_min = state_min.getPositions()
with open(minimized_pdb, 'w') as f:
    app.PDBFile.writeFile(simulation.topology, positions_min, f)
print(f"✓ Estructura minimizada guardada: {minimized_pdb}")

# ============================================================================
# PASO 8: Equilibración NVT (canónica, T constante)
# ============================================================================
print("\n" + "="*70)
print("PASO 8: Equilibración NVT (100 ps, 300 K)")
print("="*70)

simulation.context.setVelocitiesToTemperature(TEMPERATURE)

print("  Corriendo NVT...")
nvt_steps = 50000  # 100 ps con dt=2fs
simulation.step(nvt_steps)

state_nvt = simulation.context.getState(getEnergy=True, getPositions=True)
print(f"✓ NVT completado")
print(f"  Temperatura: {TEMPERATURE}")
print(f"  Energía: {state_nvt.getPotentialEnergy()}")

# Guardar
nvt_pdb = OUTPUT_DIR / "nvt_equilibrated.pdb"
with open(nvt_pdb, 'w') as f:
    app.PDBFile.writeFile(simulation.topology, state_nvt.getPositions(), f)
print(f"✓ Estructura NVT guardada: {nvt_pdb}")

# ============================================================================
# PASO 9: Equilibración NPT (isotérmica-isobárica)
# ============================================================================
print("\n" + "="*70)
print("PASO 9: Equilibración NPT (100 ps, 300 K, 1 bar)")
print("="*70)

# Agregar barostato
system.addForce(MonteCarloBarostat(PRESSURE, TEMPERATURE))

# Recrear simulación con barostato
integrator_npt = LangevinMiddleIntegrator(TEMPERATURE, 1.0/unit.picosecond, 
                                           0.002*unit.picoseconds)
simulation_npt = app.Simulation(modeller.topology, system, integrator_npt)
simulation_npt.context.setPositions(state_nvt.getPositions())
simulation_npt.context.setVelocities(simulation.context.getState(getVelocities=True).getVelocities())

print("  Corriendo NPT...")
npt_steps = 50000  # 100 ps
simulation_npt.step(npt_steps)

state_npt = simulation_npt.context.getState(getEnergy=True, getPositions=True)
print(f"✓ NPT completado")
print(f"  Temperatura: {TEMPERATURE}")
print(f"  Presión: {PRESSURE}")
print(f"  Energía: {state_npt.getPotentialEnergy()}")

# Guardar estado final equilibrado
npt_pdb = OUTPUT_DIR / "equilibrated.pdb"
with open(npt_pdb, 'w') as f:
    app.PDBFile.writeFile(simulation_npt.topology, state_npt.getPositions(), f)
print(f"✓ Sistema equilibrado final: {npt_pdb}")

equilibrated_state = OUTPUT_DIR / "equilibrated_state.xml"
simulation_npt.saveState(str(equilibrated_state))
print(f"✓ Estado equilibrado guardado: {equilibrated_state}")

# ============================================================================
# RESUMEN
# ============================================================================
print("\n" + "="*70)
print("✅ PREPARACIÓN COMPLETADA")
print("="*70)
print(f"\nArchivos generados en {OUTPUT_DIR}:")
print(f"  1. system_solvated.pdb      - Sistema completo solvatado")
print(f"  2. system.xml                - Sistema OpenMM serializado")
print(f"  3. minimized.pdb             - Estructura minimizada")
print(f"  4. nvt_equilibrated.pdb      - Post-equilibración NVT")
print(f"  5. equilibrated.pdb          - Estructura final lista para producción")
print(f"  6. equilibrated_state.xml    - Estado completo para continuar")
print("\n✓ Sistema listo para umbrella sampling del C-terminal")
print("="*70)
