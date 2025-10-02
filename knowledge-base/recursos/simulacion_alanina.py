#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simulación de Alanina Dipéptido con OpenMM
==========================================

Este script realiza una dinámica molecular básica de alanina dipéptido,
un sistema clásico para aprender simulaciones MD.

Autor: Equipo Dinámicas Moleculares UNAM-INN
Fecha: Octubre 2025
"""

import openmm as mm
from openmm import app, unit
import numpy as np
from sys import stdout

def descargar_pdb(pdb_id="1ALA", output_file="alanine_dipeptide.pdb"):
    """
    Descarga una estructura PDB desde la base de datos RCSB.
    
    Parámetros:
    -----------
    pdb_id : str
        Código PDB de 4 caracteres
    output_file : str
        Nombre del archivo de salida
    """
    import urllib.request
    
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"📥 Descargando {pdb_id} desde RCSB PDB...")
    
    try:
        urllib.request.urlretrieve(url, output_file)
        print(f"✅ Estructura guardada en: {output_file}")
        return output_file
    except Exception as e:
        print(f"❌ Error al descargar: {e}")
        return None

def configurar_sistema(pdb_file):
    """
    Configura el sistema molecular para la simulación.
    
    Parámetros:
    -----------
    pdb_file : str
        Ruta al archivo PDB
    
    Retorna:
    --------
    system : openmm.System
        Sistema configurado
    integrator : openmm.Integrator
        Integrador de Langevin
    simulation : openmm.app.Simulation
        Objeto de simulación
    """
    print("\n🔧 Configurando sistema molecular...")
    
    # Cargar estructura
    pdb = app.PDBFile(pdb_file)
    
    # Seleccionar campo de fuerza (AMBER14 con agua implícita)
    forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    
    # Crear el sistema
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.NoCutoff,  # Sin corte para sistema pequeño
        constraints=app.HBonds,         # Restringir enlaces con H
        rigidWater=True                 # Agua rígida
    )
    
    # Configurar integrador de Langevin (NVT)
    temperature = 300 * unit.kelvin
    friction = 1.0 / unit.picosecond
    timestep = 2.0 * unit.femtosecond
    
    integrator = mm.LangevinIntegrator(temperature, friction, timestep)
    
    # Crear simulación
    simulation = app.Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    
    print("✅ Sistema configurado correctamente")
    return system, integrator, simulation

def minimizar_energia(simulation):
    """
    Minimiza la energía del sistema.
    
    Parámetros:
    -----------
    simulation : openmm.app.Simulation
        Objeto de simulación
    """
    print("\n⚡ Minimizando energía...")
    
    # Energía inicial
    state = simulation.context.getState(getEnergy=True)
    energia_inicial = state.getPotentialEnergy()
    print(f"   Energía inicial: {energia_inicial}")
    
    # Minimización
    simulation.minimizeEnergy(maxIterations=100)
    
    # Energía final
    state = simulation.context.getState(getEnergy=True)
    energia_final = state.getPotentialEnergy()
    print(f"   Energía final: {energia_final}")
    print(f"   Cambio: {energia_final - energia_inicial}")

def equilibrar_sistema(simulation, pasos=1000):
    """
    Equilibra el sistema a la temperatura deseada.
    
    Parámetros:
    -----------
    simulation : openmm.app.Simulation
        Objeto de simulación
    pasos : int
        Número de pasos de equilibración
    """
    print(f"\n🌡️  Equilibrando sistema ({pasos} pasos)...")
    simulation.step(pasos)
    print("✅ Equilibración completada")

def ejecutar_produccion(simulation, pasos=10000, output_dcd="trajectory.dcd", 
                        output_pdb="output.pdb", intervalo_guardado=100):
    """
    Ejecuta la simulación de producción.
    
    Parámetros:
    -----------
    simulation : openmm.app.Simulation
        Objeto de simulación
    pasos : int
        Número de pasos de producción
    output_dcd : str
        Archivo de trayectoria
    output_pdb : str
        Archivo PDB de salida
    intervalo_guardado : int
        Guardar cada N pasos
    """
    print(f"\n🚀 Ejecutando simulación de producción ({pasos} pasos)...")
    
    # Configurar reporteros
    simulation.reporters.append(
        app.DCDReporter(output_dcd, intervalo_guardado)
    )
    
    simulation.reporters.append(
        app.StateDataReporter(
            stdout, 
            intervalo_guardado, 
            step=True,
            potentialEnergy=True,
            temperature=True,
            speed=True,
            progress=True,
            remainingTime=True,
            totalSteps=pasos
        )
    )
    
    # Ejecutar simulación
    simulation.step(pasos)
    
    # Guardar estructura final
    positions = simulation.context.getState(getPositions=True).getPositions()
    app.PDBFile.writeFile(
        simulation.topology,
        positions,
        open(output_pdb, 'w')
    )
    
    print(f"\n✅ Simulación completada")
    print(f"   Trayectoria guardada en: {output_dcd}")
    print(f"   Estructura final en: {output_pdb}")

def main():
    """Función principal."""
    print("="*60)
    print("🧬 SIMULACIÓN DE ALANINA DIPÉPTIDO CON OPENMM")
    print("="*60)
    
    # Paso 1: Descargar estructura
    pdb_file = descargar_pdb("1ALA", "alanine_dipeptide.pdb")
    
    if pdb_file is None:
        print("❌ No se pudo descargar la estructura. Abortando.")
        return
    
    # Paso 2: Configurar sistema
    system, integrator, simulation = configurar_sistema(pdb_file)
    
    # Paso 3: Minimizar energía
    minimizar_energia(simulation)
    
    # Paso 4: Equilibración
    equilibrar_sistema(simulation, pasos=1000)
    
    # Paso 5: Producción
    ejecutar_produccion(
        simulation,
        pasos=10000,
        output_dcd="alanine_trajectory.dcd",
        output_pdb="alanine_final.pdb",
        intervalo_guardado=100
    )
    
    print("\n" + "="*60)
    print("🎉 ¡Simulación completada exitosamente!")
    print("="*60)
    print("\n📊 Próximos pasos:")
    print("   1. Analiza la trayectoria con el notebook de análisis")
    print("   2. Visualiza con: alanine_analysis.ipynb")
    print("   3. Calcula diagramas de Ramachandran")

if __name__ == "__main__":
    main()
