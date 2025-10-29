#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Análisis rápido de resultados de ProPKa para WNK1

Lee propka_results.pka y reporta residuos críticos
"""

import sys
import re
from pathlib import Path

# Configuración
WORK_DIR = Path(__file__).parent
PROPKA_FILE = WORK_DIR / "prepared_system" / "propka_results.pka"

# Regiones funcionales de WNK1 (basado en PDB 5DRB, residuos 194-483)
REGIONS = {
    "N-lobe": (194, 239),
    "Catalytic_loop": (240, 250),
    "Hinge": (251, 270),
    "C-lobe": (271, 419),
    "Activation_segment": (420, 450),
    "C-terminal": (451, 483)  # Región de umbrella sampling
}

print("="*70)
print("ANÁLISIS DE ESTADOS DE PROTONACIÓN - WNK1")
print("="*70)

if not PROPKA_FILE.exists():
    print(f"\n❌ ERROR: No se encuentra {PROPKA_FILE}")
    print("   Ejecuta primero: python prepare_system.py")
    sys.exit(1)

print(f"\nArchivo: {PROPKA_FILE}")
print()

# Leer archivo ProPKa
with open(PROPKA_FILE, 'r') as f:
    content = f.read()
    lines = content.split('\n')

# Buscar sección de summary
summary_start = None
for i, line in enumerate(lines):
    if 'SUMMARY OF THIS PREDICTION' in line:
        summary_start = i
        break

if summary_start is None:
    print("⚠️  No se encontró SUMMARY en el archivo ProPKa")
    sys.exit(1)

# Parsear residuos
residues = []
for line in lines[summary_start:]:
    # Formato: ASP  45 A     3.65       3.80     20 %
    match = re.search(r'(\w{3})\s+(\d+)\s+([A-Z])\s+([\d.]+)\s+([\d.]+)\s+(\d+)\s*%', line)
    if match:
        resname = match.group(1)
        resid = int(match.group(2))
        chain = match.group(3)
        pka_pred = float(match.group(4))
        pka_model = float(match.group(5))
        buried = int(match.group(6))
        
        residues.append({
            'name': resname,
            'id': resid,
            'chain': chain,
            'pka_pred': pka_pred,
            'pka_model': pka_model,
            'buried': buried
        })

print(f"✓ {len(residues)} residuos titratables detectados")
print()

# ============================================================================
# ANÁLISIS POR TIPO DE RESIDUO
# ============================================================================

print("="*70)
print("RESUMEN POR TIPO DE RESIDUO")
print("="*70)
print()

res_types = {}
for res in residues:
    if res['name'] not in res_types:
        res_types[res['name']] = []
    res_types[res['name']].append(res)

for resname in sorted(res_types.keys()):
    count = len(res_types[resname])
    avg_pka = sum(r['pka_pred'] for r in res_types[resname]) / count
    print(f"{resname:3s}: {count:3d} residuos, pKa promedio = {avg_pka:.2f}")

# ============================================================================
# HISTIDINAS CRÍTICAS
# ============================================================================

if 'HIS' in res_types:
    print()
    print("="*70)
    print("HISTIDINAS (CRÍTICO PARA PROTONACIÓN)")
    print("="*70)
    print()
    print(f"{'Res':>6s} {'pKa':>6s} {'Buried':>7s} {'Estado pH 7.0':>20s} {'Región':>25s}")
    print("-"*70)
    
    for his in res_types['HIS']:
        resid = his['id']
        pka = his['pka_pred']
        buried = his['buried']
        
        # Predicción de estado
        if pka > 7.5:
            state = "HIP (protonated +1)"
        elif pka > 6.5:
            state = "HID/HIE (ambiguous)"
        else:
            state = "Neutral (HIE/HID)"
        
        # Región funcional
        region = "Unknown"
        for region_name, (start, end) in REGIONS.items():
            if start <= resid <= end:
                region = region_name
                break
        
        # Resaltar críticos
        marker = ""
        if region in ["Catalytic_loop", "Activation_segment", "C-terminal"]:
            marker = " ⚠️"
        
        print(f"HIS {resid:3d} {pka:6.2f} {buried:6d}% {state:>20s} {region:>25s}{marker}")
    
    print()
    print("⚠️  = Región crítica para función/umbrella sampling")

# ============================================================================
# RESIDUOS CON pKa PERTURBADO
# ============================================================================

print()
print("="*70)
print("RESIDUOS CON pKa PERTURBADO (|ΔpKa| > 1.0)")
print("="*70)
print()

perturbed = [r for r in residues if abs(r['pka_pred'] - r['pka_model']) > 1.0]

if perturbed:
    print(f"{'Res':>8s} {'pKa pred':>9s} {'pKa model':>10s} {'ΔpKa':>7s} {'Buried':>7s}")
    print("-"*70)
    
    for res in sorted(perturbed, key=lambda x: abs(x['pka_pred'] - x['pka_model']), reverse=True):
        dpka = res['pka_pred'] - res['pka_model']
        print(f"{res['name']} {res['id']:3d} {res['pka_pred']:9.2f} {res['pka_model']:10.2f} "
              f"{dpka:+7.2f} {res['buried']:6d}%")
    
    print()
    print("ΔpKa > 0: pKa aumentado (menos ácido)")
    print("ΔpKa < 0: pKa disminuido (más ácido)")
else:
    print("  No hay residuos con perturbaciones significativas")

# ============================================================================
# ANÁLISIS POR REGIÓN FUNCIONAL
# ============================================================================

print()
print("="*70)
print("ANÁLISIS POR REGIÓN FUNCIONAL")
print("="*70)
print()

for region_name, (start, end) in REGIONS.items():
    region_residues = [r for r in residues if start <= r['id'] <= end]
    
    if region_residues:
        print(f"\n{region_name} (residuos {start}-{end}):")
        print(f"  {len(region_residues)} residuos titratables:")
        
        # Contar por tipo
        types_count = {}
        for r in region_residues:
            types_count[r['name']] = types_count.get(r['name'], 0) + 1
        
        for resname, count in sorted(types_count.items()):
            print(f"    {resname}: {count}")
        
        # Mostrar HIS si hay
        his_in_region = [r for r in region_residues if r['name'] == 'HIS']
        if his_in_region:
            print(f"  HIS en esta región:")
            for his in his_in_region:
                state = "HIP" if his['pka_pred'] > 7.5 else "HID/HIE" if his['pka_pred'] > 6.5 else "Neutral"
                print(f"    HIS {his['id']}: pKa={his['pka_pred']:.2f} → {state}")

# ============================================================================
# RECOMENDACIONES
# ============================================================================

print()
print("="*70)
print("RECOMENDACIONES")
print("="*70)
print()

# HIS cerca de pH 7.0
ambiguous_his = [r for r in residues if r['name'] == 'HIS' and 6.0 < r['pka_pred'] < 8.0]

if ambiguous_his:
    print("⚠️  Histidinas ambiguas (pKa cerca de pH 7.0):")
    for his in ambiguous_his:
        print(f"    HIS {his['id']}: pKa = {his['pka_pred']:.2f}")
        
        # Región
        region = "Unknown"
        for region_name, (start, end) in REGIONS.items():
            if start <= his['id'] <= end:
                region = region_name
                break
        
        if region == "C-terminal":
            print(f"      → CRÍTICO: En región de umbrella sampling")
            print(f"      → Considerar múltiples estados de protonación (ensemble)")
        elif region == "Catalytic_loop":
            print(f"      → CRÍTICO: En loop catalítico")
            print(f"      → Verificar coordenadas de Mg²⁺ si hay")
    print()

# Residuos muy buried
very_buried = [r for r in residues if r['buried'] > 50 and abs(r['pka_pred'] - r['pka_model']) > 0.5]

if very_buried:
    print("⚠️  Residuos muy enterrados con pKa perturbado:")
    for res in very_buried[:5]:  # Primeros 5
        print(f"    {res['name']} {res['id']}: Buried {res['buried']}%, "
              f"ΔpKa = {res['pka_pred'] - res['pka_model']:+.2f}")
    if len(very_buried) > 5:
        print(f"    ... y {len(very_buried)-5} más")
    print()

# C-terminal específico
cterm_residues = [r for r in residues if 451 <= r['id'] <= 483]
if cterm_residues:
    print("✓ Región C-terminal (objetivo de umbrella sampling):")
    print(f"    {len(cterm_residues)} residuos titratables")
    
    cterm_his = [r for r in cterm_residues if r['name'] == 'HIS']
    if cterm_his:
        print(f"    {len(cterm_his)} Histidinas - Estados de protonación pueden afectar PMF")
    
    cterm_charged = [r for r in cterm_residues if r['name'] in ['ASP', 'GLU', 'LYS', 'ARG']]
    print(f"    {len(cterm_charged)} residuos cargados (ASP/GLU/LYS/ARG)")
    print()

print("="*70)
print("✓ Análisis completado")
print("="*70)
print()
print("Próximos pasos:")
print("  1. Revisar HIS críticas manualmente (PyMOL/VMD)")
print("  2. Si hay HIS ambiguas en C-terminal, considerar múltiples protonaciones")
print("  3. Continuar con: python generate_umbrella_windows.py")
print()
