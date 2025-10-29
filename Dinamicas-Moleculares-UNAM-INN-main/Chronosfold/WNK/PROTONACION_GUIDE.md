# Guía de Protonación para Simulaciones MD

## ¿Por qué es importante ProPKa?

La **protonación correcta** de residuos es crítica para:

1. **Carga neta correcta del sistema** → Afecta neutralización con iones
2. **Puentes de hidrógeno** → Estabilidad estructural
3. **Interacciones electrostáticas** → Energía y conformación
4. **Actividad catalítica** → Especialmente en kinasas como WNK1

## Residuos Titratables (pH-dependientes)

| Residuo | pKa estándar | Estado pH 7.0 | Importancia |
|---------|--------------|---------------|-------------|
| **ASP** | 3.9 | Desprotonado (-) | Carga negativa, coordinación de iones |
| **GLU** | 4.3 | Desprotonado (-) | Carga negativa, sitios catalíticos |
| **HIS** | 6.0 | Variable (±/0) | **CRÍTICO**: puede estar protonado o no |
| **CYS** | 8.3 | Protonado (SH) | Puentes disulfuro, sitios redox |
| **LYS** | 10.5 | Protonado (+) | Carga positiva, unión ATP |
| **ARG** | 12.5 | Protonado (+) | Carga positiva siempre |
| **TYR** | 10.1 | Protonado (OH) | Fosforilación en kinasas |
| **N-term** | 9.6 | Protonado (+) | Extremo amino |
| **C-term** | 2.2 | Desprotonado (-) | Extremo carboxilo |

## ⚠️ Caso Crítico: HISTIDINA (HIS)

Histidina es **especialmente problemática** porque su pKa (~6.0) está **cerca de pH fisiológico (7.0-7.4)**:

### Estados posibles de HIS:

```
      δ           ε
      |           |
   HN-C     C----NH
      ||   ||
      C----C
      |
      R

1. HID (δ-protonated): N_δ tiene H, N_ε no
2. HIE (ε-protonated): N_ε tiene H, N_δ no  
3. HIP (doubly protonated): Ambos tienen H (carga +1)
4. HIS (default): OpenMM decide automáticamente
```

### ¿Por qué importa en WNK1?

- **Sitio activo de kinasas**: Frecuentemente contienen HIS en el loop catalítico
- **Coordinación de Mg²⁺**: Iones metálicos en sitio de unión a ATP
- **Regulación alostérica**: Cambios de protonación pueden activar/desactivar

## Pipeline de Protonación en prepare_system.py

### PASO 3: ProPKa (NUEVO)

```python
from propka.run import single

# Ejecuta ProPKa a pH 7.0
molecule = single('structure.pdb', optargs=['--pH', '7.0'])

# Output: propka_results.pka
# Contiene pKa predicho para cada residuo titratable
```

**ProPKa predice pKa considerando**:
- Entorno local (buried vs exposed)
- Puentes de hidrógeno cercanos
- Interacciones electrostáticas
- Desolvation penalties

### PASO 4: addHydrogens con pH

```python
forcefield = app.ForceField('amber14-all.xml')
modeller.addHydrogens(forcefield, pH=7.0)
```

**OpenMM heuristics** (cuando ProPKa no está disponible):
- **HIS**: Compara energías de HID/HIE/HIP, elige la de menor energía
- **ASP/GLU**: Desprotonados si pH > pKa
- **LYS/ARG**: Protonados siempre (pKa >> 7)

## Interpretando resultados de ProPKa

### Archivo propka_results.pka

```
SUMMARY OF THIS PREDICTION
       Group      pKa  model-pKa   buried
   ASP  45 A     3.65       3.80     20 %
   GLU  67 A     4.45       4.50      0 %
   HIS 123 A     7.20       6.50     45 %    ← ATENCIÓN!
   LYS 189 A    10.48      10.53      5 %
```

**Interpretación**:
- **HIS 123**: pKa = 7.20 → A pH 7.0 está ~50% protonado
  - **Buried 45%** → Menos accesible, pKa shift hacia arriba
  - **Decisión**: Probablemente **HIE** (ε-protonated)

### Casos especiales en WNK1

1. **Sitio de unión a ATP**: 
   - Si hay HIS cerca del fosfato β, probablemente **HID** (interacción con Mg²⁺)

2. **Loop de activación**:
   - HIS en loop flexible → Considerar múltiples protonaciones (ensemble)

3. **Interface C-terminal ↔ dominio kinasa**:
   - HIS en la interface → Puede modular conformaciones (umbrella sampling)

## Validación de Protonación

### 1. Inspección visual (PyMOL/VMD)

```python
# Cargar prepared_system/equilibrated.pdb
# Seleccionar histidinas
select his, resn HIS+HID+HIE+HIP
show sticks, his
show spheres, his and name ND1+NE2

# Verificar puentes de hidrógeno
distance hbonds, his, all, 3.5
```

### 2. Verificar carga neta

```python
# En el output de prepare_system.py
# Buscar: "Iones (Na+/Cl-)"
# 
# Carga neta proteína = N_HIP - N_ASP - N_GLU + N_LYS + N_ARG + 1 (N-term) - 1 (C-term)
# 
# Debe ser neutralizada por: N_CL - N_NA
```

### 3. Energía después de minimización

```python
# E_potencial muy alta (>0 kJ/mol) → Posible mala protonación
# E_potencial razonable (<-1e5 kJ/mol) → Probablemente OK
```

## Casos donde ProPKa puede fallar

1. **Sitios metálicos**: ProPKa no considera Mg²⁺/Zn²⁺ explícitamente
   - **Solución**: Revisar manualmente HIS cerca de iones

2. **Ligandos covalentes**: Si hay inhibidores unidos
   - **Solución**: Usar PDB2PQR con ligandos parametrizados

3. **Mutaciones patogénicas**: Cambios de pKa no predichos
   - **Solución**: Comparar WT vs mutante con ProPKa

## Workflow Recomendado

```bash
# 1. Ejecutar prepare_system.py (incluye ProPKa automáticamente)
python prepare_system.py

# 2. Revisar propka_results.pka
cat prepared_system/propka_results.pka | grep HIS

# 3. Si hay HIS críticas, re-ejecutar con estados específicos:
# (Actualmente no implementado - requiere edición manual de PDB)

# 4. Validar con MD corta
python run_umbrella_window.py --window 0 --steps 10000 --platform CPU

# 5. Verificar estabilidad (RMSD, energía)
```

## Recursos

- **ProPKa documentation**: https://github.com/jensengroup/propka
- **PDB2PQR**: https://www.poissonboltzmann.org/
- **H++ server** (alternativa web): http://biophysics.cs.vt.edu/H++

## Para WNK1 específicamente

### Residuos clave a verificar:

1. **Catalytic loop** (residuos ~240-250):
   - HIS en el bolsillo de ATP
   - Protonación afecta afinidad por ATP/Mg²⁺

2. **Activation segment** (residuos ~420-440):
   - Fosforilación de TYR/SER/THR
   - Protonación de HIS puede estabilizar estado activado

3. **C-terminal** (residuos 451-483):
   - **Objetivo de umbrella sampling**
   - HIS en interface C-term ↔ kinase domain → Modulan conformaciones
   - Protonación puede cambiar preferencias compacto vs extendido

### Script de análisis post-ProPKa

```python
# analyze_propka_wnk.py
import re

with open('prepared_system/propka_results.pka', 'r') as f:
    lines = f.readlines()

print("Residuos HIS en WNK1:")
for line in lines:
    if 'HIS' in line:
        match = re.search(r'HIS\s+(\d+)\s+A\s+([\d.]+)', line)
        if match:
            resid = int(match.group(1))
            pka = float(match.group(2))
            
            # Predicción de estado
            if pka > 7.5:
                state = "HIP (protonated, +1)"
            elif pka > 6.5:
                state = "HID/HIE (neutral, protonation state ambiguous)"
            else:
                state = "Neutral (likely HIE or HID)"
            
            # Contexto funcional
            if 240 <= resid <= 250:
                context = "** CATALYTIC LOOP **"
            elif 420 <= resid <= 440:
                context = "** ACTIVATION SEGMENT **"
            elif 451 <= resid <= 483:
                context = "** C-TERMINAL (umbrella CV region) **"
            else:
                context = ""
            
            print(f"  HIS {resid}: pKa={pka:.2f} → {state} {context}")
```

---

**IMPORTANTE**: Después de ejecutar `prepare_system.py`, **SIEMPRE** revisar `propka_results.pka` antes de lanzar simulaciones largas en HPC.
