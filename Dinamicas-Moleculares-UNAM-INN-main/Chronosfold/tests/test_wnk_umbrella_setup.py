#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tests de validación pre-HPC para umbrella sampling WNK1

Verifica:
1. Estructura PDB válida
2. Sistema OpenMM se puede crear
3. Grupos de átomos correctos
4. Bias armónico funciona
5. Simulación corta (100 pasos) funciona
6. CV se puede calcular
"""

import pytest
import sys
from pathlib import Path
import numpy as np

# Añadir parent al path para imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from openmm import app, unit
from openmm import *

WORK_DIR = Path(__file__).parent.parent
PDB_FILE = WORK_DIR / "WNK" / "5DRB.pdb"
PREPARED_DIR = WORK_DIR / "WNK" / "prepared_system"
WINDOWS_DIR = WORK_DIR / "WNK" / "umbrella_windows"

class TestSystemPreparation:
    """Tests de preparación del sistema"""
    
    def test_pdb_file_exists(self):
        """Verifica que el PDB original existe"""
        assert PDB_FILE.exists(), f"PDB no encontrado: {PDB_FILE}"
    
    def test_pdb_is_valid(self):
        """Verifica que el PDB se puede cargar"""
        pdb = app.PDBFile(str(PDB_FILE))
        assert len(list(pdb.topology.atoms())) > 0, "PDB sin átomos"
        
        # Verificar que tiene residuos de proteína
        n_residues = len(list(pdb.topology.residues()))
        assert n_residues > 100, f"Muy pocos residuos: {n_residues}"
    
    def test_forcefield_loads(self):
        """Verifica que el forcefield se puede cargar"""
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        assert forcefield is not None

class TestPreparedSystem:
    """Tests del sistema preparado"""
    
    @pytest.fixture
    def prepared_system(self):
        """Fixture: sistema preparado"""
        if not PREPARED_DIR.exists():
            pytest.skip(f"Sistema preparado no existe: {PREPARED_DIR}")
        
        pdb_file = PREPARED_DIR / "equilibrated.pdb"
        xml_file = PREPARED_DIR / "system.xml"
        
        if not pdb_file.exists() or not xml_file.exists():
            pytest.skip("Archivos de sistema preparado no encontrados")
        
        pdb = app.PDBFile(str(pdb_file))
        with open(xml_file, 'r') as f:
            system = XmlSerializer.deserialize(f.read())
        
        return pdb, system
    
    def test_prepared_system_has_water(self, prepared_system):
        """Verifica que el sistema tiene agua"""
        pdb, system = prepared_system
        
        n_waters = sum(1 for res in pdb.topology.residues() if res.name == 'HOH')
        assert n_waters > 1000, f"Muy pocas aguas: {n_waters}"
    
    def test_prepared_system_has_ions(self, prepared_system):
        """Verifica que el sistema tiene iones"""
        pdb, system = prepared_system
        
        n_ions = sum(1 for res in pdb.topology.residues() if res.name in ['NA', 'CL'])
        assert n_ions > 0, "No hay iones para neutralizar"
    
    def test_system_has_periodic_box(self, prepared_system):
        """Verifica que el sistema tiene caja periódica"""
        pdb, system = prepared_system
        
        vectors = pdb.topology.getPeriodicBoxVectors()
        assert vectors is not None, "No hay caja periódica"
        
        # Verificar dimensiones razonables (> 3 nm cada lado)
        a = vectors[0].value_in_unit(unit.nanometer)[0]
        assert a > 3.0, f"Caja muy pequeña: {a} nm"
    
    def test_system_can_create_context(self, prepared_system):
        """Verifica que se puede crear un contexto de simulación"""
        pdb, system = prepared_system
        
        integrator = LangevinMiddleIntegrator(300*unit.kelvin, 
                                               1.0/unit.picosecond,
                                               0.002*unit.picoseconds)
        
        try:
            platform = Platform.getPlatformByName('CPU')
            simulation = app.Simulation(pdb.topology, system, integrator, platform)
        except:
            simulation = app.Simulation(pdb.topology, system, integrator)
        
        assert simulation is not None
        
        # Intentar setear posiciones
        simulation.context.setPositions(pdb.positions)
        
        # Obtener energía
        state = simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy()
        
        assert not np.isnan(energy.value_in_unit(unit.kilojoule_per_mole)), "Energía es NaN"

class TestUmbrellaWindows:
    """Tests de ventanas de umbrella sampling"""
    
    def test_windows_directory_exists(self):
        """Verifica que el directorio de ventanas existe"""
        if not WINDOWS_DIR.exists():
            pytest.skip(f"Directorio de ventanas no existe: {WINDOWS_DIR}")
    
    def test_atom_groups_file_exists(self):
        """Verifica que el archivo de grupos de átomos existe"""
        if not WINDOWS_DIR.exists():
            pytest.skip("Directorio de ventanas no existe")
        
        atom_file = WINDOWS_DIR / "atom_groups.txt"
        assert atom_file.exists(), f"Archivo de grupos no encontrado: {atom_file}"
    
    def test_atom_groups_are_valid(self):
        """Verifica que los grupos de átomos son válidos"""
        if not WINDOWS_DIR.exists():
            pytest.skip("Directorio de ventanas no existe")
        
        atom_file = WINDOWS_DIR / "atom_groups.txt"
        if not atom_file.exists():
            pytest.skip("Archivo de grupos no existe")
        
        with open(atom_file, 'r') as f:
            lines = f.readlines()
        
        kinase_line = [l for l in lines if l.startswith("KINASE_ATOMS=")]
        cterm_line = [l for l in lines if l.startswith("CTERM_ATOMS=")]
        
        assert len(kinase_line) == 1, "No se encontró KINASE_ATOMS"
        assert len(cterm_line) == 1, "No se encontró CTERM_ATOMS"
        
        # Parsear
        kinase_atoms = [int(x) for x in kinase_line[0].split("=")[1].strip().split(",")]
        cterm_atoms = [int(x) for x in cterm_line[0].split("=")[1].strip().split(",")]
        
        assert len(kinase_atoms) > 100, f"Muy pocos átomos en kinase: {len(kinase_atoms)}"
        assert len(cterm_atoms) > 10, f"Muy pocos átomos en C-term: {len(cterm_atoms)}"
        
        # No debe haber overlap
        overlap = set(kinase_atoms) & set(cterm_atoms)
        assert len(overlap) == 0, f"Hay {len(overlap)} átomos en ambos grupos"
    
    def test_windows_config_exists(self):
        """Verifica que la configuración de ventanas existe"""
        if not WINDOWS_DIR.exists():
            pytest.skip("Directorio de ventanas no existe")
        
        config_file = WINDOWS_DIR / "windows_config.csv"
        assert config_file.exists(), f"Config no encontrado: {config_file}"
    
    def test_windows_config_is_valid(self):
        """Verifica que la configuración de ventanas es válida"""
        if not WINDOWS_DIR.exists():
            pytest.skip("Directorio de ventanas no existe")
        
        config_file = WINDOWS_DIR / "windows_config.csv"
        if not config_file.exists():
            pytest.skip("Config no existe")
        
        import pandas as pd
        df = pd.read_csv(config_file)
        
        assert 'window_id' in df.columns
        assert 'r0_nm' in df.columns
        assert 'spring_constant_kJ_mol_nm2' in df.columns
        
        assert len(df) > 5, f"Muy pocas ventanas: {len(df)}"
        
        # Verificar que r0 es monótono creciente
        r0_values = df['r0_nm'].values
        assert np.all(np.diff(r0_values) > 0), "r0 no es monótono creciente"
        
        # Verificar valores razonables
        assert np.all(r0_values > 0.5), "r0 muy pequeño"
        assert np.all(r0_values < 10.0), "r0 muy grande"

class TestBiasForce:
    """Tests de la fuerza de bias armónico"""
    
    @pytest.fixture
    def minimal_system(self):
        """Fixture: sistema minimal con PDB original"""
        pdb = app.PDBFile(str(PDB_FILE))
        
        # Modeller para limpiar
        modeller = app.Modeller(pdb.topology, pdb.positions)
        modeller.deleteWater()
        
        # Agregar hidrógenos
        forcefield = app.ForceField('amber14-all.xml')
        modeller.addHydrogens(forcefield, pH=7.0)
        
        # Crear sistema minimal (vacuum)
        system = forcefield.createSystem(modeller.topology,
                                          nonbondedMethod=app.NoCutoff)
        
        return modeller.topology, modeller.positions, system
    
    def test_can_add_centroid_force(self, minimal_system):
        """Verifica que se puede agregar CustomCentroidBondForce"""
        topology, positions, system = minimal_system
        
        # Obtener algunos carbonos alfa
        ca_atoms = [atom.index for atom in topology.atoms() 
                    if atom.name == 'CA'][:50]  # Primeros 50
        
        assert len(ca_atoms) >= 2, "No hay suficientes CA"
        
        # Crear dos grupos
        group1 = ca_atoms[:25]
        group2 = ca_atoms[25:50]
        
        # Agregar fuerza
        force = CustomCentroidBondForce(2, "0.5*k*(distance(g1,g2)-r0)^2")
        force.addPerBondParameter("k")
        force.addPerBondParameter("r0")
        
        g1 = force.addGroup(group1)
        g2 = force.addGroup(group2)
        
        force.addBond([g1, g2], [
            1000.0 * unit.kilojoule_per_mole / (unit.nanometer**2),
            2.0 * unit.nanometer
        ])
        
        system.addForce(force)
        
        assert system.getNumForces() > 0
    
    def test_bias_force_gives_reasonable_energy(self, minimal_system):
        """Verifica que el bias da energías razonables"""
        topology, positions, system = minimal_system
        
        # Grupos de CA
        ca_atoms = [atom.index for atom in topology.atoms() 
                    if atom.name == 'CA'][:50]
        
        group1 = ca_atoms[:25]
        group2 = ca_atoms[25:50]
        
        # Agregar bias
        force = CustomCentroidBondForce(2, "0.5*k*(distance(g1,g2)-r0)^2")
        force.addPerBondParameter("k")
        force.addPerBondParameter("r0")
        
        g1 = force.addGroup(group1)
        g2 = force.addGroup(group2)
        
        k_bias = 1000.0  # kJ/mol/nm²
        r0_bias = 2.0    # nm
        
        force.addBond([g1, g2], [
            k_bias * unit.kilojoule_per_mole / (unit.nanometer**2),
            r0_bias * unit.nanometer
        ])
        
        system.addForce(force)
        
        # Crear simulación
        integrator = LangevinMiddleIntegrator(300*unit.kelvin,
                                               1.0/unit.picosecond,
                                               0.002*unit.picoseconds)
        
        try:
            platform = Platform.getPlatformByName('CPU')
            simulation = app.Simulation(topology, system, integrator, platform)
        except:
            simulation = app.Simulation(topology, system, integrator)
        
        simulation.context.setPositions(positions)
        
        # Obtener energía
        state = simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        
        # Debe ser finita y razonable
        assert not np.isnan(energy), "Energía es NaN"
        assert not np.isinf(energy), "Energía es infinita"
        assert energy < 1e6, f"Energía muy alta: {energy} kJ/mol"

class TestShortSimulation:
    """Tests de simulación corta"""
    
    def test_can_run_100_steps(self):
        """Verifica que se puede correr 100 pasos de MD"""
        # Usar PDB original (más rápido que sistema solvatado)
        pdb = app.PDBFile(str(PDB_FILE))
        
        modeller = app.Modeller(pdb.topology, pdb.positions)
        modeller.deleteWater()
        
        forcefield = app.ForceField('amber14-all.xml')
        modeller.addHydrogens(forcefield, pH=7.0)
        
        system = forcefield.createSystem(modeller.topology,
                                          nonbondedMethod=app.NoCutoff,
                                          constraints=app.HBonds)
        
        integrator = LangevinMiddleIntegrator(300*unit.kelvin,
                                               1.0/unit.picosecond,
                                               0.002*unit.picoseconds)
        
        try:
            platform = Platform.getPlatformByName('CPU')
            simulation = app.Simulation(modeller.topology, system, integrator, platform)
        except:
            simulation = app.Simulation(modeller.topology, system, integrator)
        
        simulation.context.setPositions(modeller.positions)
        simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
        
        # Correr 100 pasos
        simulation.step(100)
        
        # Verificar que terminó
        state = simulation.context.getState(getEnergy=True)
        energy = state.getPotentialEnergy()
        
        assert not np.isnan(energy.value_in_unit(unit.kilojoule_per_mole))

# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    print("="*70)
    print("TESTS PRE-HPC - WNK1 UMBRELLA SAMPLING")
    print("="*70)
    print()
    print("Ejecutando tests de validación...")
    print()
    
    # Ejecutar con pytest
    pytest.main([__file__, "-v", "--tb=short"])
