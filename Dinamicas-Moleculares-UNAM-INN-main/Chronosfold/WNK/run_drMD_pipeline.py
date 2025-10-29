#!/usr/bin/env python3
"""
Pipeline paralelo usando drMD para WNK1 Umbrella Sampling
===========================================================

Este script ejecuta drMD para preparar y equilibrar WNK1 en condiciones PBS,
complementando el pipeline principal de Chronosfold.

Ventajas de drMD:
  - Manejo automático de errores (FirstAid)
  - Logging detallado
  - Clustering automático
  - Paralelización nativa
  - Generación automática de sección de métodos

Uso:
  python run_drMD_pipeline.py

Requiere:
  - drMD instalado (pip install drMD)
  - AmberTools 23
  - OpenMM
  - drMD_wnk_config.yaml en mismo directorio
"""

import sys
import os
from pathlib import Path
import yaml

# Verificar que drMD esté instalado
try:
    # Importar módulo drMD desde SciToolAgent/external
    drmd_path = Path(__file__).parent.parent.parent / "SciToolAgent" / "external" / "drMD" / "src"
    sys.path.insert(0, str(drmd_path))
    
    import drMD
    print("✓ drMD encontrado en SciToolAgent/external/drMD")
except ImportError:
    print("⚠️  ERROR: drMD no encontrado")
    print("\nInstalación:")
    print("  1. Instalar desde pip:")
    print("     pip install drMD")
    print("\n  2. O usar versión en SciToolAgent/external:")
    print("     cd ../../SciToolAgent/external/drMD")
    print("     pip install -e .")
    sys.exit(1)

# Verificar que AmberTools esté disponible
try:
    import parmed
    print("✓ AmberTools disponible (parmed)")
except ImportError:
    print("⚠️  ADVERTENCIA: AmberTools no detectado")
    print("   conda install -c conda-forge ambertools=23")

# Verificar OpenMM
try:
    import openmm
    print(f"✓ OpenMM {openmm.__version__} disponible")
except ImportError:
    print("⚠️  ERROR: OpenMM no disponible")
    print("   conda install -c conda-forge openmm")
    sys.exit(1)


def check_config_file(config_path: Path) -> dict:
    """Verificar que archivo de configuración exista y sea válido"""
    if not config_path.exists():
        print(f"⚠️  ERROR: {config_path} no encontrado")
        print("\nEjecutar primero:")
        print("  cp drMD_wnk_config.yaml.template drMD_wnk_config.yaml")
        sys.exit(1)
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    print(f"✓ Configuración válida: {config_path}")
    return config


def check_pdb_file(pdb_path: Path):
    """Verificar que estructura PDB exista"""
    if not pdb_path.exists():
        print(f"⚠️  ERROR: {pdb_path} no encontrado")
        print("\nAsegúrate que 5DRB.pdb esté en el directorio actual")
        sys.exit(1)
    
    print(f"✓ Estructura encontrada: {pdb_path}")


def modify_config_for_pbs(config: dict) -> dict:
    """
    Modificar configuración para PBS buffer
    
    NOTA: drMD/OpenMM no soportan K+ y fosfatos directamente
    Esta función documenta la limitación y usa aproximación
    """
    print("\n" + "="*70)
    print("CONFIGURACIÓN PBS")
    print("="*70)
    
    # Verificar que pH sea 7.4
    if config.get('miscInfo', {}).get('pH') != 7.4:
        print("⚠️  Ajustando pH a 7.4 (PBS)")
        if 'miscInfo' not in config:
            config['miscInfo'] = {}
        config['miscInfo']['pH'] = 7.4
    
    print("\nPBS buffer (aproximado):")
    print("  pH: 7.4")
    print("  Fuerza iónica: ~0.163 M (Na+/Cl-)")
    print("\nNOTA: K+ y fosfatos no soportados en forcefield amber14")
    print("      Usando fuerza iónica equivalente")
    
    return config


def run_drmd_pipeline(config_path: Path):
    """Ejecutar pipeline completo de drMD"""
    print("\n" + "="*70)
    print("EJECUTANDO drMD PIPELINE")
    print("="*70)
    
    # Verificar archivos
    config = check_config_file(config_path)
    pdb_file = Path("5DRB.pdb")
    check_pdb_file(pdb_file)
    
    # Ajustar configuración para PBS
    config = modify_config_for_pbs(config)
    
    # Guardar configuración modificada
    temp_config = Path("drMD_wnk_config_pbs.yaml")
    with open(temp_config, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)
    
    print(f"\n✓ Configuración PBS guardada: {temp_config}")
    
    # Ejecutar drMD
    print("\nIniciando drMD...")
    print("Esto puede tardar 1-2 horas dependiendo del hardware\n")
    
    try:
        drMD.main(str(temp_config))
        print("\n" + "="*70)
        print("✓ PIPELINE drMD COMPLETADO")
        print("="*70)
        
        # Reportar outputs
        output_dir = Path(config['pathInfo']['outputDir'])
        print(f"\nResultados en: {output_dir}")
        print("\nArchivos importantes:")
        print(f"  - {output_dir}/00_drMD_logs/aftercare.log")
        print(f"  - {output_dir}/5DRB/production/")
        print(f"  - {output_dir}/5DRB/production_backbone.dcd")
        
        # Si hay clustering
        cluster_dir = output_dir / "5DRB" / "05_cluster_analysis"
        if cluster_dir.exists():
            print(f"  - {cluster_dir}/clusters/")
        
        print("\nPróximos pasos:")
        print("  1. Revisar logs: cat drMD_output/00_drMD_logs/aftercare.log")
        print("  2. Visualizar: pymol drMD_output/5DRB/production/production_*.pdb")
        print("  3. Analizar trayectoria para definir umbrella windows")
        
    except Exception as e:
        print(f"\n⚠️  ERROR en drMD: {e}")
        print("\nRevisar logs en: drMD_output/00_drMD_logs/")
        sys.exit(1)


def main():
    """Punto de entrada principal"""
    print("="*70)
    print("WNK1 PIPELINE PARALELO CON drMD")
    print("="*70)
    print("\nEste pipeline complementa Chronosfold umbrella sampling")
    print("Usa drMD para preparación y equilibración automatizada\n")
    
    # Archivo de configuración
    config_file = Path("drMD_wnk_config.yaml")
    
    # Ejecutar
    run_drmd_pipeline(config_file)


if __name__ == "__main__":
    main()
