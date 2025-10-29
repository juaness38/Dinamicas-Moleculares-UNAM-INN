#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MD de producción con bias armónico para una ventana de umbrella sampling

Usage:
    python run_umbrella_window.py --window 0
    
Este script:
1. Carga la estructura inicial
2. Agrega bias armónico en la distancia entre dominios
3. Corre MD de producción (configurable)
4. Guarda trayectoria y valores de CV
"""

import argparse
import sys
from pathlib import Path
import numpy as np
from openmm import app, unit
from openmm import *

def parse_args():
    parser = argparse.ArgumentParser(description="Umbrella sampling MD para WNK1")
    parser.add_argument('--window', type=int, required=True,
                       help='Índice de ventana (0-19)')
    parser.add_argument('--steps', type=int, default=5000000,
                       help='Pasos de simulación (default: 5M = 10 ns)')
    parser.add_argument('--dt', type=float, default=0.002,
                       help='Timestep en ps (default: 2 fs)')
    parser.add_argument('--temp', type=float, default=300,
                       help='Temperatura en K (default: 300)')
    parser.add_argument('--platform', type=str, default='CPU',
                       choices=['CPU', 'CUDA', 'OpenCL'],
                       help='Plataforma de cálculo (default: CPU)')
    parser.add_argument('--save-interval', type=int, default=5000,
                       help='Intervalo para guardar frames (default: 5000 steps = 10 ps)')
    parser.add_argument('--log-interval', type=int, default=10000,
                       help='Intervalo para reportar progreso (default: 10000 steps)')
    return parser.parse_args()

def load_atom_groups(atom_groups_file):
    """Carga índices de átomos desde atom_groups.txt"""
    with open(atom_groups_file, 'r') as f:
        lines = f.readlines()
    
    kinase_line = [l for l in lines if l.startswith("KINASE_ATOMS=")][0]
    cterm_line = [l for l in lines if l.startswith("CTERM_ATOMS=")][0]
    
    kinase_atoms = [int(x) for x in kinase_line.split("=")[1].strip().split(",")]
    cterm_atoms = [int(x) for x in cterm_line.split("=")[1].strip().split(",")]
    
    return kinase_atoms, cterm_atoms

def compute_com_distance(positions, kinase_atoms, cterm_atoms):
    """Calcula distancia entre centros de masa"""
    kinase_coords = np.array([positions[i].value_in_unit(unit.nanometer) for i in kinase_atoms])
    cterm_coords = np.array([positions[i].value_in_unit(unit.nanometer) for i in cterm_atoms])
    
    com_kinase = np.mean(kinase_coords, axis=0)
    com_cterm = np.mean(cterm_coords, axis=0)
    
    return np.linalg.norm(com_kinase - com_cterm)

def main():
    args = parse_args()
    
    print("="*70)
    print(f"UMBRELLA SAMPLING - VENTANA {args.window}")
    print("="*70)
    
    # Paths
    WORK_DIR = Path(__file__).parent
    WINDOWS_DIR = WORK_DIR / "umbrella_windows"
    WINDOW_DIR = WINDOWS_DIR / f"window_{args.window:02d}"
    
    if not WINDOW_DIR.exists():
        print(f"❌ ERROR: Ventana {args.window} no existe")
        print(f"   Esperado en: {WINDOW_DIR}")
        print("   Ejecuta primero generate_umbrella_windows.py")
        sys.exit(1)
    
    INPUT_PDB = WINDOW_DIR / "initial.pdb"
    SYSTEM_XML = WORK_DIR / "prepared_system" / "system.xml"
    ATOM_GROUPS_FILE = WINDOWS_DIR / "atom_groups.txt"
    WINDOWS_CONFIG = WINDOWS_DIR / "windows_config.csv"
    
    print(f"\n📁 Directorio de ventana: {WINDOW_DIR}")
    
    # Leer parámetros de esta ventana
    import pandas as pd
    windows_df = pd.read_csv(WINDOWS_CONFIG)
    window_params = windows_df.iloc[args.window]
    
    r0 = window_params['r0_nm']
    spring_k = window_params['spring_constant_kJ_mol_nm2']
    
    print(f"\n⚙️  PARÁMETROS:")
    print(f"  Ventana: {args.window}")
    print(f"  r₀ = {r0:.4f} nm")
    print(f"  k = {spring_k:.1f} kJ/mol/nm²")
    print(f"  Temperatura: {args.temp} K")
    print(f"  Pasos: {args.steps:,} ({args.steps * args.dt / 1000:.1f} ns)")
    print(f"  Timestep: {args.dt} ps")
    print(f"  Plataforma: {args.platform}")
    
    # ============================================================================
    # PASO 1: Cargar sistema
    # ============================================================================
    print("\n" + "="*70)
    print("PASO 1: Cargando sistema")
    print("="*70)
    
    pdb = app.PDBFile(str(INPUT_PDB))
    print(f"✓ PDB cargado: {len(list(pdb.topology.atoms()))} átomos")
    
    with open(SYSTEM_XML, 'r') as f:
        system = XmlSerializer.deserialize(f.read())
    print(f"✓ Sistema cargado")
    
    # Cargar grupos de átomos
    kinase_atoms, cterm_atoms = load_atom_groups(ATOM_GROUPS_FILE)
    print(f"✓ Grupos de átomos: Kinase={len(kinase_atoms)}, C-term={len(cterm_atoms)}")
    
    # Calcular distancia inicial
    d_initial = compute_com_distance(pdb.positions, kinase_atoms, cterm_atoms)
    print(f"  Distancia inicial: {d_initial:.3f} nm (target: {r0:.3f} nm)")
    
    # ============================================================================
    # PASO 2: Agregar bias armónico
    # ============================================================================
    print("\n" + "="*70)
    print("PASO 2: Agregando bias armónico")
    print("="*70)
    
    # Fuerza: U = 0.5 * k * (r - r0)²
    # donde r = distancia entre centros de masa
    bias_force = CustomCentroidBondForce(2, "0.5*k*(distance(g1,g2)-r0)^2")
    bias_force.addPerBondParameter("k")
    bias_force.addPerBondParameter("r0")
    
    # Agregar grupos (centros de masa)
    g1 = bias_force.addGroup(kinase_atoms)
    g2 = bias_force.addGroup(cterm_atoms)
    
    # Agregar bond con parámetros de esta ventana
    bias_force.addBond([g1, g2], [
        spring_k * unit.kilojoule_per_mole / (unit.nanometer**2),
        r0 * unit.nanometer
    ])
    
    system.addForce(bias_force)
    print(f"✓ Bias agregado: U = 0.5 * {spring_k} * (r - {r0})²")
    
    # ============================================================================
    # PASO 3: Crear simulación
    # ============================================================================
    print("\n" + "="*70)
    print("PASO 3: Configurando simulación")
    print("="*70)
    
    # Integrador Langevin
    temperature = args.temp * unit.kelvin
    friction = 1.0 / unit.picosecond
    timestep = args.dt * unit.picoseconds
    
    integrator = LangevinMiddleIntegrator(temperature, friction, timestep)
    
    # Plataforma
    try:
        platform = Platform.getPlatformByName(args.platform)
        print(f"✓ Plataforma: {args.platform}")
    except Exception as e:
        print(f"⚠️  Advertencia: Plataforma {args.platform} no disponible")
        print(f"   Usando plataforma por defecto")
        platform = None
    
    # Crear simulación
    if platform:
        simulation = app.Simulation(pdb.topology, system, integrator, platform)
    else:
        simulation = app.Simulation(pdb.topology, system, integrator)
    
    simulation.context.setPositions(pdb.positions)
    simulation.context.setVelocitiesToTemperature(temperature)
    
    print(f"✓ Simulación creada")
    
    # Energía inicial
    state = simulation.context.getState(getEnergy=True)
    E_pot = state.getPotentialEnergy()
    E_kin = state.getKineticEnergy()
    print(f"  Energía potencial inicial: {E_pot}")
    print(f"  Energía cinética inicial: {E_kin}")
    
    # ============================================================================
    # PASO 4: Configurar reporters
    # ============================================================================
    print("\n" + "="*70)
    print("PASO 4: Configurando salida")
    print("="*70)
    
    # DCD trajectory
    dcd_file = WINDOW_DIR / "trajectory.dcd"
    simulation.reporters.append(
        app.DCDReporter(str(dcd_file), args.save_interval)
    )
    print(f"✓ Trayectoria DCD: {dcd_file}")
    print(f"  Frecuencia: cada {args.save_interval} pasos ({args.save_interval * args.dt} ps)")
    
    # Estado termodinámico
    log_file = WINDOW_DIR / "production.log"
    simulation.reporters.append(
        app.StateDataReporter(
            str(log_file), args.log_interval,
            step=True, time=True, potentialEnergy=True, kineticEnergy=True,
            totalEnergy=True, temperature=True, speed=True
        )
    )
    print(f"✓ Log termodinámico: {log_file}")
    
    # Reportero personalizado para CV
    cv_file = WINDOW_DIR / "cv_values.dat"
    cv_output = open(cv_file, 'w')
    cv_output.write("# step time_ps distance_nm bias_energy_kJ\n")
    
    class CVReporter:
        def __init__(self, file, interval, kinase_atoms, cterm_atoms, r0, k):
            self.file = file
            self.interval = interval
            self.kinase_atoms = kinase_atoms
            self.cterm_atoms = cterm_atoms
            self.r0 = r0
            self.k = k
            self._step = 0
        
        def __del__(self):
            self.file.close()
        
        def report(self, simulation, state):
            self._step += 1
            if self._step % self.interval != 0:
                return
            
            positions = state.getPositions()
            d = compute_com_distance(positions, self.kinase_atoms, self.cterm_atoms)
            
            # Energía de bias
            bias_E = 0.5 * self.k * (d - self.r0)**2
            
            step = state.getStepCount()
            time = state.getTime().value_in_unit(unit.picosecond)
            
            self.file.write(f"{step} {time:.3f} {d:.6f} {bias_E:.6f}\n")
            self.file.flush()
        
        def describeNextReport(self, simulation):
            steps_left = self.interval - (self._step % self.interval)
            return (steps_left, False, False, False, False, None)
    
    cv_reporter = CVReporter(cv_output, args.save_interval, kinase_atoms, cterm_atoms, r0, spring_k)
    simulation.reporters.append(cv_reporter)
    print(f"✓ CV data: {cv_file}")
    
    # ============================================================================
    # PASO 5: Producción MD
    # ============================================================================
    print("\n" + "="*70)
    print("PASO 5: CORRIENDO MD DE PRODUCCIÓN")
    print("="*70)
    print(f"  {args.steps:,} pasos = {args.steps * args.dt / 1000:.2f} ns")
    print(f"  Esto puede tomar un tiempo...")
    print()
    
    try:
        simulation.step(args.steps)
    except KeyboardInterrupt:
        print("\n⚠️  Simulación interrumpida por usuario")
    except Exception as e:
        print(f"\n❌ ERROR durante simulación: {e}")
        raise
    
    # ============================================================================
    # PASO 6: Guardar estado final
    # ============================================================================
    print("\n" + "="*70)
    print("PASO 6: Guardando estado final")
    print("="*70)
    
    final_state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
    
    # PDB final
    final_pdb = WINDOW_DIR / "final.pdb"
    with open(final_pdb, 'w') as f:
        app.PDBFile.writeFile(simulation.topology, final_state.getPositions(), f)
    print(f"✓ Estructura final: {final_pdb}")
    
    # Estado completo
    final_state_xml = WINDOW_DIR / "final_state.xml"
    simulation.saveState(str(final_state_xml))
    print(f"✓ Estado completo: {final_state_xml}")
    
    # Estadísticas finales
    E_pot_final = final_state.getPotentialEnergy()
    E_kin_final = final_state.getKineticEnergy()
    d_final = compute_com_distance(final_state.getPositions(), kinase_atoms, cterm_atoms)
    
    print(f"\n📊 ESTADÍSTICAS FINALES:")
    print(f"  Energía potencial: {E_pot_final}")
    print(f"  Energía cinética: {E_kin_final}")
    print(f"  Distancia CV: {d_final:.3f} nm (target: {r0:.3f} nm)")
    print(f"  Desviación: {abs(d_final - r0):.3f} nm")
    
    # ============================================================================
    # RESUMEN
    # ============================================================================
    print("\n" + "="*70)
    print(f"✅ VENTANA {args.window} COMPLETADA")
    print("="*70)
    print(f"\nArchivos generados en {WINDOW_DIR}:")
    print(f"  1. trajectory.dcd        - Trayectoria MD")
    print(f"  2. production.log        - Log termodinámico")
    print(f"  3. cv_values.dat         - Valores de CV en cada frame")
    print(f"  4. final.pdb             - Estructura final")
    print(f"  5. final_state.xml       - Estado completo final")
    print("="*70)

if __name__ == "__main__":
    main()
