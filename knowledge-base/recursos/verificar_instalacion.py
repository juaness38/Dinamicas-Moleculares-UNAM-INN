#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script de Verificación de Instalación
======================================

Este script verifica que todas las dependencias necesarias
para realizar dinámicas moleculares estén correctamente instaladas.

Uso:
    python verificar_instalacion.py
"""

import sys

def verificar_paquete(nombre_paquete, import_name=None):
    """Verifica si un paquete está instalado y muestra su versión."""
    if import_name is None:
        import_name = nombre_paquete
    
    try:
        modulo = __import__(import_name)
        version = getattr(modulo, '__version__', 'desconocida')
        print(f"✅ {nombre_paquete:20s} v{version}")
        return True
    except ImportError:
        print(f"❌ {nombre_paquete:20s} NO INSTALADO")
        return False

def verificar_openmm_gpu():
    """Verifica las plataformas disponibles de OpenMM."""
    try:
        from openmm import Platform
        print("\n🔍 Plataformas OpenMM disponibles:")
        for i in range(Platform.getNumPlatforms()):
            platform = Platform.getPlatform(i)
            nombre = platform.getName()
            if nombre == "CUDA":
                print(f"   🚀 {nombre} (GPU NVIDIA - ¡Excelente!)")
            elif nombre == "OpenCL":
                print(f"   ⚡ {nombre} (GPU - Bueno)")
            elif nombre == "CPU":
                print(f"   💻 {nombre} (Optimizado)")
            else:
                print(f"   📌 {nombre}")
        return True
    except Exception as e:
        print(f"⚠️  Error al verificar plataformas: {e}")
        return False

def main():
    """Función principal de verificación."""
    print("="*60)
    print("🧬 VERIFICACIÓN DE INSTALACIÓN - Dinámicas Moleculares")
    print("="*60)
    print()
    
    # Lista de paquetes a verificar
    paquetes = [
        ("Python", "sys"),
        ("NumPy", "numpy"),
        ("SciPy", "scipy"),
        ("Matplotlib", "matplotlib"),
        ("Pandas", "pandas"),
        ("OpenMM", "openmm"),
        ("MDTraj", "mdtraj"),
        ("MDAnalysis", "MDAnalysis"),
        ("Jupyter", "jupyter"),
        ("NGLView", "nglview"),
        ("BioPython", "Bio"),
        ("Seaborn", "seaborn"),
    ]
    
    print("📦 Verificando paquetes instalados:")
    print("-" * 60)
    
    resultados = []
    for nombre, import_name in paquetes:
        resultado = verificar_paquete(nombre, import_name)
        resultados.append((nombre, resultado))
    
    # Verificar versión de Python
    print(f"\n🐍 Python: {sys.version}")
    
    # Verificar plataformas de OpenMM
    verificar_openmm_gpu()
    
    # Resumen
    print("\n" + "="*60)
    exitosos = sum(1 for _, r in resultados if r)
    total = len(resultados)
    porcentaje = (exitosos / total) * 100
    
    print(f"📊 RESUMEN: {exitosos}/{total} paquetes instalados ({porcentaje:.1f}%)")
    
    if exitosos == total:
        print("🎉 ¡Instalación completa! Estás listo para comenzar.")
    elif exitosos >= total * 0.8:
        print("⚠️  Instalación mayormente completa. Considera instalar los paquetes faltantes.")
    else:
        print("❌ Instalación incompleta. Revisa la guía de instalación.")
    
    print("="*60)
    
    # Mensaje de ayuda
    if exitosos < total:
        print("\n💡 Para instalar paquetes faltantes:")
        print("   Con conda: conda install -c conda-forge <paquete>")
        print("   Con pip:   pip install <paquete>")
        print("\n📖 Consulta: knowledge-base/GUIA_INSTALACION.md")

if __name__ == "__main__":
    main()
