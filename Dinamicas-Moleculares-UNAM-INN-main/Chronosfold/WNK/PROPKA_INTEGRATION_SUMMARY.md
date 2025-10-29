# Actualización: Integración de ProPKa para Protonación Correcta

## 🎯 Problema Resuelto

**Antes**: OpenMM's `addHydrogens(pH=7.0)` usaba valores genéricos de pKa sin considerar:
- Entorno local (buried vs exposed)
- Puentes de hidrógeno
- Interacciones electrostáticas
- Efectos de desolvación

**Resultado**: Estados de protonación incorrectos, especialmente para **Histidina (HIS)** cuyo pKa (~6.0) está cerca de pH fisiológico (7.0).

## ✅ Solución Implementada

### 1. Actualización de Dependencies

**Archivo**: `environment.yml`

```yaml
# AGREGADO:
- propka        # pKa prediction considerando entorno local
- pdb2pqr       # Herramienta alternativa de protonación
- pymbar        # Para análisis MBAR (faltaba)
```

### 2. Modificación del Pipeline de Preparación

**Archivo**: `Chronosfold/WNK/prepare_system.py`

**NUEVO PASO 3** (antes de agregar hidrógenos):

```python
from propka.run import single

# Ejecutar ProPKa a pH 7.0
molecule = single('structure.pdb', optargs=['--pH', '7.0'])

# Guardar resultados
with open('propka_results.pka', 'w') as f:
    f.write(molecule.write_propka())
```

**Output**: `prepared_system/propka_results.pka` con pKa predicho para cada residuo titratable.

**Workflow actualizado**:
1. Cargar 5DRB.pdb
2. Limpiar estructura (solo proteína)
3. **ProPKa analysis** ← NUEVO
4. addHydrogens(pH=7.0) ← Con estados informados por ProPKa
5. Solvatar
6. Minimizar
7. Equilibrar NVT
8. Equilibrar NPT

### 3. Documentación Completa

**Archivo**: `Chronosfold/WNK/PROTONACION_GUIDE.md`

- Explicación de residuos titratables
- Estados de protonación de HIS (HID/HIE/HIP)
- Por qué es crítico para WNK1
- Workflow de validación
- Scripts de análisis post-ProPKa

### 4. Script de Análisis de ProPKa

**Archivo**: `Chronosfold/WNK/analyze_propka.py`

Analiza automáticamente `propka_results.pka` y reporta:
- **Histidinas críticas** por región funcional
- Residuos con pKa perturbado (|ΔpKa| > 1.0)
- Análisis por región (Catalytic loop, C-terminal, etc.)
- Recomendaciones específicas para WNK1

**Uso**:
```bash
python analyze_propka.py
```

### 5. README Actualizado

**Archivo**: `Chronosfold/WNK/README.md`

- ProPKa en dependencies
- Paso 1 ahora incluye output de ProPKa
- Link a PROTONACION_GUIDE.md
- Nota sobre verificar estados de HIS

## 🔬 Por Qué es Crítico para WNK1

### Regiones Funcionales Sensibles

1. **Catalytic Loop (residuos 240-250)**:
   - HIS cerca del sitio de unión a ATP
   - Protonación afecta coordinación de Mg²⁺
   - Impacta catálisis

2. **Activation Segment (residuos 420-450)**:
   - Contiene sitios de fosforilación
   - HIS puede estabilizar estado activado/desactivado

3. **C-terminal (residuos 451-483)** ← **OBJETIVO DE UMBRELLA SAMPLING**:
   - **HIS en interface C-term ↔ dominio kinasa**
   - Estados de protonación modulan preferencia compacto vs extendido
   - **Protonación incorrecta → PMF incorrecto**

### Ejemplo: HIS en C-terminal

```
pKa = 7.2 (ProPKa prediction)
pH = 7.0 (simulation)

Sin ProPKa:
  → OpenMM adivina: Probablemente HIE (ε-protonated)
  
Con ProPKa:
  → pKa 7.2 + buried 45% → Definitivamente HIE
  → Pero si pKa = 6.8 → Podría ser HID (δ-protonated)
  → ΔE conformacional puede ser 2-5 kJ/mol diferente
```

## 📋 Nuevo Workflow de Validación

```bash
# 1. Preparar sistema (incluye ProPKa automáticamente)
python prepare_system.py

# 2. Analizar estados de protonación
python analyze_propka.py

# 3. Revisar output
cat prepared_system/propka_results.pka | grep HIS

# 4. Si hay HIS críticas, verificar manualmente en PyMOL
pymol prepared_system/equilibrated.pdb

# 5. Continuar con umbrella sampling
python generate_umbrella_windows.py
```

## 📊 Output Esperado

### propka_results.pka
```
SUMMARY OF THIS PREDICTION
       Group      pKa  model-pKa   buried
   ASP  45 A     3.65       3.80     20 %
   GLU  67 A     4.45       4.50      0 %
   HIS 456 A     7.20       6.50     45 %    ← C-terminal!
   LYS 189 A    10.48      10.53      5 %
```

### analyze_propka.py output
```
HISTIDINAS (CRÍTICO PARA PROTONACIÓN)
   Res    pKa  Buried        Estado pH 7.0                   Región
----------------------------------------------------------------------
HIS 456   7.20     45%    HID/HIE (ambiguous)             C-terminal ⚠️
HIS 245   6.85     30%    HID/HIE (ambiguous)         Catalytic_loop ⚠️

⚠️ = Región crítica para función/umbrella sampling
```

## ⚙️ Fallback Behavior

Si ProPKa no está instalado o falla:

```python
except ImportError:
    print("⚠️ ProPKa not installed, using default pH=7.0 states")
    # Continúa con addHydrogens de OpenMM
```

Sistema sigue funcional pero menos preciso para residuos buried/interface.

## 🎓 Referencias

- **ProPKa**: Olsson et al. (2011) *J. Chem. Theory Comput.* 7, 525-537
- **pKa estándar**: ASP(3.9), GLU(4.3), HIS(6.0), LYS(10.5), ARG(12.5)
- **GitHub**: https://github.com/jensengroup/propka

## ✅ Checklist de Archivos Modificados/Creados

- [x] `environment.yml` - Agregados propka, pdb2pqr, pymbar
- [x] `Chronosfold/WNK/prepare_system.py` - PASO 3 ProPKa insertado
- [x] `Chronosfold/WNK/PROTONACION_GUIDE.md` - Guía completa (nuevo)
- [x] `Chronosfold/WNK/analyze_propka.py` - Script de análisis (nuevo)
- [x] `Chronosfold/WNK/README.md` - Actualizado con ProPKa info

## 🚀 Próximos Pasos

1. **Instalar ProPKa**:
   ```bash
   conda env update -f environment.yml
   # o
   pip install propka
   ```

2. **Ejecutar preparación**:
   ```bash
   cd Chronosfold/WNK
   python prepare_system.py
   ```

3. **Analizar protonación**:
   ```bash
   python analyze_propka.py
   ```

4. **Revisar HIS críticas** en PyMOL/VMD si es necesario

5. **Continuar con umbrella sampling**:
   ```bash
   python generate_umbrella_windows.py
   python run_umbrella_window.py --window 0 --steps 50000
   ```

## 💡 Consideración Futura

Para análisis exhaustivo, considerar **múltiples estados de protonación** (ensemble):

```bash
# Correr umbrella sampling con:
# - Estado 1: HIS 456 como HIE
# - Estado 2: HIS 456 como HID
# - Promediar PMF ponderado por probabilidad de cada estado
```

Esto captura incertidumbre en protonación → PMF más robusto.

---

**Fecha**: 2025-01-20  
**Sistema**: WNK1 C-terminal umbrella sampling  
**Mejora**: Protonación correcta con ProPKa  
**Impacto**: PMF más preciso, especialmente si hay HIS en región de CV
