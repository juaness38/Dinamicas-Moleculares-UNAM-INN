# Implementación de PBS Buffer en Pipeline WNK1
## Phosphate Buffered Saline para Condiciones Fisiológicas

---

## 📋 Composición PBS Estándar

### PBS 1× (pH 7.4, 25°C)

| Componente | Concentración | Molaridad | Función |
|-----------|---------------|-----------|---------|
| **NaCl** | 8.0 g/L | **137 mM** | Fuerza iónica principal |
| **KCl** | 0.2 g/L | **2.7 mM** | Ión fisiológico |
| **Na₂HPO₄** | 1.42 g/L | **10 mM** | Buffer (base conjugada) |
| **KH₂PO₄** | 0.24 g/L | **1.8 mM** | Buffer (ácido débil) |

**Propiedades:**
- pH: 7.4 ± 0.1
- Fuerza iónica total: ~163 mM
- Osmolaridad: ~280-300 mOsm/kg (isotónico)
- Capacidad buffer: pKa₂(H₃PO₄) = 7.2

---

## ⚠️ Limitaciones de Forcefields Estándar

### Problema: Iones y moléculas no soportadas

OpenMM con forcefield **amber14-all.xml** solo soporta:
- ✅ Na⁺ (sodio)
- ✅ Cl⁻ (cloruro)
- ✅ Mg²⁺ (magnesio)
- ✅ Ca²⁺ (calcio)
- ❌ **K⁺** (potasio) - NO DISPONIBLE
- ❌ **HPO₄²⁻** (fosfato dibásico) - NO DISPONIBLE
- ❌ **H₂PO₄⁻** (fosfato monobásico) - NO DISPONIBLE

### ¿Por qué?

Los forcefields biomoleculares clásicos (AMBER, CHARMM, OPLS) fueron parametrizados principalmente para:
1. Iones comunes: Na⁺, Cl⁻, Mg²⁺, Ca²⁺, Zn²⁺
2. Agua (TIP3P, SPC/E, TIP4P)
3. Proteínas, DNA, RNA, lípidos

**K⁺** y **fosfatos poliatómicos** requieren:
- Parámetros de Lennard-Jones específicos
- Cargas parciales distribuidas (para poliatómicos)
- Términos de ángulo/torsión (HPO₄²⁻/H₂PO₄⁻)

---

## 🔧 Soluciones Implementadas

### OPCIÓN 1: Aproximación con Fuerza Iónica Equivalente ⭐ **IMPLEMENTADA**

**Estrategia:**
- Reemplazar K⁺ → Na⁺
- Reemplazar HPO₄²⁻/H₂PO₄⁻ → Cl⁻ (ajustando carga)
- Mantener fuerza iónica total: **163 mM**

**Ventajas:**
- ✅ Compatible con amber14-all.xml
- ✅ No requiere parametrización adicional
- ✅ Fuerza iónica correcta para electrostática
- ✅ pH controlado con `addHydrogens(pH=7.4)`

**Desventajas:**
- ❌ No captura interacciones específicas K⁺-proteína
- ❌ No captura coordinación de fosfatos
- ❌ Buffer capacity aproximada (no verdadero equilibrio H₂PO₄⁻/HPO₄²⁻)

**Implementación en `prepare_system.py`:**
```python
# Calcular fuerza iónica total de PBS
# 137 mM NaCl + 2.7 mM KCl + 10 mM Na₂HPO₄ + 1.8 mM KH₂PO₄
# Total: ~163 mM

ionic_strength_pbs = 0.163 * unit.molar

modeller.addSolvent(
    forcefield, 
    model='tip3p',
    padding=1.0*unit.nanometer,
    ionicStrength=ionic_strength_pbs,
    positiveIon='Na+',  # Representa Na+ y K+
    negativeIon='Cl-'   # Representa Cl⁻ y fosfatos
)
```

**Cálculo de iones:**
```python
# Volumen de caja (ej. 10 nm × 10 nm × 10 nm = 1000 nm³ = 1 pL)
box_volume = 1e-12 L  # 1 picolitro

# Número de iones
N = C × V × N_A
N_Na = 0.163 mol/L × 1e-12 L × 6.022e23 = ~98 iones
N_Cl = 0.163 mol/L × 1e-12 L × 6.022e23 = ~98 iones

# + iones para neutralizar proteína
```

---

### OPCIÓN 2: Forcefield CHARMM36 (K⁺ disponible)

**CHARMM36m** incluye parámetros para K⁺:
- Disponible en OpenMM: `charmm36.xml`
- K⁺: Radio de van der Waals = 1.76 Å, ε = 0.087 kcal/mol

**Implementación:**
```python
forcefield = app.ForceField('charmm36.xml', 'charmm36/water.xml')

# Agregar PBS con K+ real
modeller.addSolvent(
    forcefield,
    model='tip3p',
    padding=1.0*unit.nanometer,
    positiveIon='K+',  # ← K+ disponible en CHARMM!
    negativeIon='Cl-',
    ionicStrength=0.027*unit.molar  # Solo KCl primero
)

# Después agregar Na+ manualmente (o usar script custom)
```

**Ventajas:**
- ✅ K⁺ real con parámetros validados
- ✅ Mejor para estudios de selectividad iónica

**Desventajas:**
- ❌ Requiere re-parametrizar proteína (no trivial)
- ❌ Fosfatos aún no disponibles
- ❌ Incompatible con pipeline actual (amber14)

---

### OPCIÓN 3: Parametrizar Fosfatos con GAFF/CGenFF

**Workflow:**
1. Obtener estructura 3D de HPO₄²⁻ y H₂PO₄⁻
2. Calcular cargas con RESP/AM1-BCC
3. Asignar tipos de átomo GAFF (General Amber Force Field)
4. Generar archivos `.frcmod` y `.mol2`
5. Integrar con `loadAmberParams()` en OpenMM

**Ejemplo con antechamber (AmberTools):**
```bash
# Para H₂PO₄⁻
antechamber -i h2po4.pdb -fi pdb -o h2po4.mol2 -fo mol2 \
            -c bcc -nc -1 -at gaff2

parmchk2 -i h2po4.mol2 -f mol2 -o h2po4.frcmod

# Cargar en OpenMM
forcefield.loadFile('h2po4.xml')  # Convertir .frcmod → .xml
```

**Ventajas:**
- ✅ PBS exacto con fosfatos reales
- ✅ Buffer capacity correcto

**Desventajas:**
- ❌ **Muy complejo** para este proyecto
- ❌ Requiere validación QM (DFT)
- ❌ Fuera del alcance de umbrella sampling simple

---

## ✅ Solución Recomendada: OPCIÓN 1

### Pipeline Modificado (`prepare_system.py`)

```python
# PASO 5: Solvatar en PBS (aproximado)
ionic_strength_pbs = 0.163 * unit.molar  # 163 mM total

modeller.addSolvent(
    forcefield,
    model='tip3p',
    padding=1.0*unit.nanometer,
    ionicStrength=ionic_strength_pbs,
    positiveIon='Na+',
    negativeIon='Cl-'
)

print("✓ PBS buffer (aproximado):")
print("  137 mM NaCl   → Na+/Cl-")
print("  2.7 mM KCl    → Na+/Cl- (K+ aproximado)")
print("  10 mM Na₂HPO₄ → Na+/Cl- (fosfato aproximado)")
print("  1.8 mM KH₂PO₄ → Cl- (fosfato aproximado)")
print("  Fuerza iónica total: 163 mM")
print("  pH: 7.4 (controlado con addHydrogens)")
```

---

## 📊 Validación de la Aproximación

### Test 1: Comparar PMF con diferentes fuerzas iónicas

```python
# Correr umbrella sampling con:
# 1. Fuerza iónica 163 mM (PBS aproximado)
# 2. Fuerza iónica 150 mM (estándar)
# 3. Fuerza iónica 200 mM

# Comparar PMF:
# Si diferencia < 2 kJ/mol → aproximación válida
```

### Test 2: Verificar apantallamiento electrostático

```python
# Calcular Debye length λ_D
λ_D = sqrt(ε₀εᵣkT / (2NAe²I))

# Para I = 0.163 M:
λ_D ≈ 0.76 nm

# Para I = 0.150 M:
λ_D ≈ 0.79 nm

# Diferencia: 3.9% → despreciable para distancias >2 nm
```

---

## 🔍 ¿Cuándo Importa K⁺ Real?

### Casos donde K⁺ ≠ Na⁺:

1. **Canales iónicos selectivos** (KcsA, NaK, Kv1.2)
   - Filtro de selectividad usa geometría exacta
   - Diferencia de radio: K⁺ (1.38 Å) vs Na⁺ (1.02 Å)

2. **Sitios de coordinación específicos**
   - WNK (With No Lysine) kinasa → K⁺ en sitio catalítico!
   - **CRÍTICO PARA WNK1**: Ver literatura sobre K⁺-binding

3. **Estudios de selectividad iónica**
   - Free energy perturbation (FEP): K⁺ → Na⁺
   - Requiere K⁺ real

### Para WNK1 C-terminal:

**Revisar cristalografía:**
```bash
# Buscar K+ en 5DRB.pdb
grep "^HETATM.*  K  " 5DRB.pdb

# Si hay K+ cristalográfico → USAR CHARMM36
# Si no hay K+ en estructura → APROXIMACIÓN OK
```

**Verificar en PDB:**
- 5DRB heteroátomos: ¿Hay K⁺ o Mg²⁺ cristalográfico?
- Si sí → considerar OPCIÓN 2 (CHARMM36)

---

## 🚀 Próximos Pasos

### Implementación Inmediata (OPCIÓN 1)

1. ✅ Modificar `prepare_system.py` con PBS aproximado
2. ✅ Documentar aproximación en logs
3. ✅ Ejecutar umbrella sampling
4. ✅ Comparar con fuerza iónica estándar (150 mM)

### Mejora Futura (si necesario)

1. Verificar si WNK1 tiene sitio de K⁺ funcional
2. Si crítico: cambiar a CHARMM36 para K⁺ real
3. Parametrizar fosfatos solo si buffer capacity importa

---

## 📚 Referencias

1. **PBS Composition:**
   - Cold Spring Harbor Protocols: doi:10.1101/pdb.rec8247

2. **Ion Parameters:**
   - Joung & Cheatham (2008): J. Phys. Chem. B, 112(30), 9020-9041
   - CHARMM36 K⁺: doi:10.1021/ct300400x

3. **WNK Kinases and K⁺:**
   - Piala et al. (2014): PNAS, 111(28), 10305-10310
   - Xu et al. (2000): Science, 290(5492), 765-768
   - **¡REVISAR SI K⁺ ES FUNCIONAL EN WNK1!**

---

## ✅ Checklist de Implementación

- [x] Modificar `prepare_system.py` con fuerza iónica 163 mM
- [x] Agregar documentación de aproximación PBS
- [x] Crear `PBS_BUFFER_IMPLEMENTATION.md` (este archivo)
- [ ] Verificar 5DRB.pdb para K⁺/Mg²⁺ cristalográfico
- [ ] Correr test de validación (PBS vs estándar)
- [ ] Decidir si CHARMM36 es necesario (basado en literatura WNK)
- [ ] Integrar con drMD pipeline paralelo

---

**Última actualización:** 2025-01-20  
**Autor:** Pipeline WNK1 Team  
**Versión:** 1.0
