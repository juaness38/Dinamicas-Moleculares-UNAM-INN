# Umbrella Sampling: Trayectorias Sesgadas, Termodinámica Correcta

**Pregunta clave**: "¿Si las trayectorias están sesgadas, cómo puedo analizarlas? ¿Cómo obtengo PMF correcto?"

**Respuesta corta**: Las trayectorias SÍ están sesgadas (eso es intencional), pero **MBAR/WHAM descuentan el bias matemáticamente** para recuperar termodinámica correcta. **NO recuperas cinética real**.

---

## 🚨 La Paradoja de Umbrella Sampling

### Tu Intuición es Correcta

> "Si aplico bias artificial → trayectoria no representa dinámica real → ¿cómo confío en el resultado?"

**Esta preocupación es 100% válida**. La clave está en entender:

1. **"Sesgado" ≠ "Incorrecto"** (son conceptos diferentes)
2. **Termodinámica ≠ Cinética** (umbrella calcula solo la primera)
3. **MBAR es magia matemática** (recupera propiedades sin sesgo)

---

## 1. ¿Qué Significa "Sesgado"?

### MD Estándar (Sin Sesgo)

```
Sistema explora según:
P(s) ∝ exp(-F(s) / kT)

Donde F(s) = PMF real (lo que queremos conocer)

Problema: Si F(barrera) = 25 kJ/mol
         → P(barrera) ≈ 10⁻⁴
         → Sistema evita la barrera (tardaría μs en cruzar)
```

**Analogía**: Dejas una pelota en un valle profundo. Esperarías años a que salte espontáneamente por fluctuación térmica.

### Umbrella Sampling (Con Sesgo)

```
Sistema explora según:
P_i(s) ∝ exp(-(F(s) + U_i(s)) / kT)

Donde U_i(s) = k(s - s_i)² = Bias armónico artificial

Efecto: Altera artificialmente P(s) para muestrear región
```

**Analogía**: Pones la pelota en una cuerda elástica centrada en la barrera. La pelota DEBE explorar la barrera (forzado por el resorte).

### Visualización

```
PMF Real (sin sesgo):
    
    25 kJ/mol  ___
              /   \
             /     \___
    0 kJ/mol/  A      B

P(s) sin sesgo: 99% en A, 0.01% en barrera, 1% en B


PMF Efectivo en Ventana i (con bias):

          ___/▲\___  ← Bias armónico centrado en barrera
    25   /   |   \
        /    |    \___
    0  /  A  i  B

P_i(s) con sesgo: 80% en i (barrera), 10% en A, 10% en B

¡Sistema FORZADO a explorar barrera!
```

---

## 2. ¿Cómo MBAR "Descuenta" el Bias?

### El Truco Matemático

MBAR (Multistate Bennett Acceptance Ratio) usa **teoría de re-weighting**:

Si conoces:
1. El bias que aplicaste: $U_i(s)$
2. Las probabilidades observadas (sesgadas): $P_i^{obs}(s)$

Puedes calcular:
3. Las probabilidades reales (sin sesgo): $P^{real}(s)$

### Ecuación de Re-weighting

$$P^{unbiased}(s) = \frac{\sum_{i=1}^{K} N_i \, P_i^{obs}(s) \cdot \exp(+U_i(s)/k_B T)}{\sum_{i=1}^{K} N_i \cdot \exp(+U_i(s)/k_B T + f_i)}$$

Donde:
- $K$ = Número de ventanas (20 en nuestro caso)
- $N_i$ = Número de muestras en ventana $i$
- $U_i(s)$ = Bias aplicado en ventana $i$
- $f_i$ = "Free energy" de cada ventana (calculado iterativamente por MBAR)

**En palabras**:
1. Observas histograma sesgado en cada ventana
2. MBAR calcula "cuánto bias agregaste" ($U_i$)
3. **"Resta" ese bias exponencialmente** ($\exp(+U_i/kT)$ compensa $\exp(-U_i/kT)$ del muestreo)
4. Recuperas $P(s)$ sin sesgo
5. Extraes PMF: $F(s) = -k_B T \ln P(s)$

### Analogía Fotográfica

```
Proceso:
1. Tomas 20 fotos de montaña con FILTROS de color diferentes
2. Sabes exactamente qué filtro usaste en cada foto
3. Photoshop "invierte" los filtros digitalmente
4. Recuperas COLOR REAL de la montaña

Umbrella = Fotos con filtros (trayectorias sesgadas)
MBAR = Photoshop (invierte filtros matemáticamente)
PMF = Montaña con color real (termodinámica correcta)
```

**Clave**: No necesitas trayectoria sin filtro. Basta conocer el filtro para corregirlo.

---

## 3. ¿Qué Puedes Analizar de Trayectorias Sesgadas?

### ✅ ANÁLISIS VÁLIDOS (Propiedades Termodinámicas)

#### A. PMF (Free Energy Profile)

**LO MÁS IMPORTANTE**: Después de MBAR, obtienes PMF **sin sesgo**.

```python
# analyze_umbrella_mbar.py
pmf, pmf_error = mbar.compute_pmf()

# Este PMF es CORRECTO (bias descontado)
barrier = pmf.max()  # ΔG‡ = barrera de activación
```

**Validez**: ✅ 100% correcto termodinámicamente

#### B. Estructuras Representativas

Puedes extraer conformaciones de cada ventana:

```python
# Ventana 10 (CV = 3.0 nm, la barrera)
traj = md.load('window_10_trajectory.dcd', top='system.pdb')

# Clusterizar para obtener representantes
from sklearn.cluster import KMeans
clusters = KMeans(n_clusters=5).fit(traj.xyz)
representative_TS = traj[clusters.labels_ == 0]

# Guardar estructura del estado de transición
representative_TS[0].save_pdb('transition_state_structure.pdb')
```

**Validez**: ✅ Estructuras SON reales
- El bias solo afecta **cuánto tiempo pasas en CV=3.0 nm**
- NO afecta **qué estructura tiene el sistema cuando CV=3.0 nm**
- Las conformaciones observadas son físicamente correctas

#### C. Contactos y Propiedades Locales

```python
# Analizar contactos en la barrera (ventana 10)
contacts = md.compute_contacts(traj, contacts=[[5, 120], [8, 115]])

# Estos contactos SON válidos
# El bias te LLEVÓ a CV=3.0, pero NO te dice cómo lograrlo
# El sistema decide naturalmente cómo configurarse
```

**Análisis válidos**:
- ✅ Puentes de hidrógeno en cada ventana
- ✅ RMSD respecto a cristal
- ✅ Radio de giro (Rg)
- ✅ Ángulos diedros (Ramachandran)
- ✅ SASA (superficie accesible al solvente)
- ✅ Distancias específicas (que NO sean la CV)

**Por qué son válidos**: El bias solo restringe **1 coordenada** (la CV). Las otras ~120,000 coordenadas (todos los átomos) evolucionan naturalmente.

#### D. Poblaciones Relativas

```python
# Probabilidad de estar en conformación abierta vs cerrada
P_open = integrate(P(s), s=2.0 to 2.5 nm)
P_closed = integrate(P(s), s=3.5 to 4.0 nm)

ratio = P_open / P_closed
```

**Validez**: ✅ Correcto después de MBAR
- $P(s)$ es sin sesgo (MBAR lo garantiza)
- Puedes calcular poblaciones, entropías, etc.

---

### ❌ ANÁLISIS INVÁLIDOS (Propiedades Cinéticas)

#### A. Tiempos de Residencia

```python
# ❌ INCORRECTO:
residence_time = measure_time_in_region(traj, cv_min=2.9, cv_max=3.1)

# Este tiempo está ARTIFICIALMENTE EXTENDIDO por el bias
# El sistema "debería" salir rápido, pero el resorte lo retiene
```

**Problema**: El bias $k(s - s_i)^2$ crea una **trampa artificial**. El tiempo que pasas en cada ventana NO refleja cinética real.

#### B. Frecuencia de Transiciones

```python
# ❌ INCORRECTO:
num_crossings = count_barrier_crossings(traj, barrier_cv=3.0)
k_on = num_crossings / simulation_time

# Esta frecuencia es ARTIFICIAL
# El bias facilita/dificulta cruces dependiendo de la ventana
```

**Problema**: En ventanas adyacentes a la barrera, el bias **empuja** hacia ella. Ves muchas transiciones que no ocurrirían naturalmente.

#### C. Constantes de Velocidad (k_on, k_off)

```python
# ❌ INCORRECTO:
k_off = 1 / mean_residence_time_in_A

# Esta tasa es COMPLETAMENTE SESGADA
```

**Problema**: Umbrella NO conserva la cinética. Solo conserva termodinámica (PMF).

#### D. Coeficientes de Difusión

```python
# ❌ INCORRECTO:
D = calculate_diffusion_coefficient(traj)

# El bias artificial afecta cómo el sistema difunde en el espacio de CV
```

**Problema**: El resorte armónico ralentiza/acelera artificialmente el movimiento.

---

## 4. Termodinámica vs. Cinética: La Distinción Clave

### Tabla Comparativa

| Propiedad | Tipo | Umbrella + MBAR | MD Estándar (μs) | Necesario para |
|-----------|------|-----------------|------------------|----------------|
| **PMF (ΔG)** | Termodinámica | ✅ Correcto | ✅ Correcto | Estabilidad relativa |
| **Poblaciones** | Termodinámica | ✅ Correcto | ✅ Correcto | Equilibrio |
| **Estructuras** | Termodinámica | ✅ Correcto | ✅ Correcto | Modelado |
| **Constante k** | Cinética | ❌ Sesgado | ✅ Correcto | Velocidad de reacción |
| **Tiempo t₁/₂** | Cinética | ❌ Sesgado | ✅ Correcto | Farmacología |
| **Mecanismo** | Cinética | ⚠️ Parcial | ✅ Completo | Pathway detallado |

### Pregunta: ¿Qué Necesitas para Tu Proyecto?

**Para WNK1 C-terminal**:

**Pregunta biológica**: *"¿Cuál es el costo energético de la apertura del C-terminal?"*

- **Respuesta**: ΔG‡ = 24.7 kJ/mol (barrera de activación)
- **Método apropiado**: ✅ Umbrella sampling (suficiente)
- **Paper reporta**: PMF, estructuras del TS, contactos clave
- **NO necesitas**: k_off, t₁/₂ (eso es otro paper completo)

**Si la pregunta fuera**: *"¿Cuánto tiempo tarda la apertura del C-terminal?"*

- **Respuesta**: t = 1/k_on (necesitas cinética)
- **Método apropiado**: ❌ Umbrella NO sirve
- **Alternativas**: Weighted Ensemble, Milestoning, MD ultra-largo

---

## 5. Alternativas Si Necesitas Dinámica No Sesgada

### Opción 1: MD Estándar Ultra-Largo (Fuerza Bruta)

```python
# Correr microsegundos sin bias
simulation.step(5_000_000_000)  # 10 μs en GPU

# Esperar a observar transiciones espontáneas
transitions = detect_transitions(trajectory)
k_off = len(transitions) / total_time
```

**Ventajas**:
- ✅ Dinámica 100% real (sin bias)
- ✅ Cinética correcta (k, t₁/₂)
- ✅ Mecanismo completo

**Desventajas**:
- ❌ Extremadamente costoso (semanas en GPU)
- ❌ Puede no ver transiciones (barrera muy alta)
- ❌ Requiere múltiples réplicas para estadística

**Para WNK1**: Barrera ~25 kJ/mol → P(cruzar) ~ 10⁻⁴ → necesitarías >100 μs

---

### Opción 2: Weighted Ensemble (WESTPA)

```python
# Divide espacio en bins
# Corre múltiples trayectorias cortas (sin bias)
# Replica/elimina trayectorias para mantener poblaciones uniformes

import westpa
sim = westpa.Simulation()
sim.run(n_iterations=1000)

# Obtiene k directamente
k_off = sim.calculate_rate()
```

**Ventajas**:
- ✅ Dinámica sin bias (trayectorias naturales)
- ✅ Calcula k_on, k_off directamente
- ✅ Más eficiente que MD puro (factor 10-100×)

**Desventajas**:
- ⚠️ Complejo de configurar (curva de aprendizaje)
- ⚠️ Requiere buenos criterios de binning
- ⚠️ Análisis más sofisticado

**Para WNK1**: Factible, pero proyecto completo (3-6 meses)

---

### Opción 3: Milestoning

```
Divide pathway en "hitos" (milestones)
Corre trayectorias cortas entre hitos adyacentes
Reconstruye cinética global desde trayectorias locales
```

**Ventajas**:
- ✅ Obtiene k sin bias
- ✅ Puede usar información de umbrella (posiciones de hitos)
- ✅ Menos costoso que MD puro

**Desventajas**:
- ⚠️ Requiere definir hitos (¿cuántos? ¿dónde?)
- ⚠️ Asume difusión entre hitos (puede ser incorrecto)
- ⚠️ Implementación no trivial

**Para WNK1**: Posible como extensión (después de umbrella)

---

### Opción 4: Transition Path Sampling (TPS)

```
Muestrea directamente rutas reactivas (A → B)
Usa Monte Carlo en el espacio de trayectorias
Obtiene mecanismo Y cinética
```

**Ventajas**:
- ✅ Mecanismo detallado (pathway más probable)
- ✅ Cinética correcta
- ✅ No necesita CV predefinida

**Desventajas**:
- ❌ MUY costoso computacionalmente
- ❌ Requiere expertise avanzado (papers de métodos)
- ❌ Análisis complejo

**Para WNK1**: Proyecto de doctorado completo

---

### Tabla de Decisión

| Objetivo | Método Recomendado | Tiempo de Proyecto | Dificultad |
|----------|-------------------|-------------------|------------|
| Solo ΔG (barrera) | **Umbrella / Metadinámica** | 2-4 semanas | ⭐⭐ |
| ΔG + validación | Umbrella + Metadinámica | 4-6 semanas | ⭐⭐⭐ |
| k_on, k_off (cinética) | Weighted Ensemble | 3-6 meses | ⭐⭐⭐⭐ |
| Mecanismo detallado | Milestoning / TPS | 6-12 meses | ⭐⭐⭐⭐⭐ |
| Todo lo anterior | MD ultra-largo (fuerza bruta) | Años | ⭐⭐⭐ (simple pero lento) |

---

## 6. Ejemplo Práctico: WNK1 C-Terminal

### Lo Que Umbrella TE DA

**Después de correr el pipeline completo**:

```bash
# 1. PMF (barrera de activación)
python analyze_umbrella_mbar.py
# → ΔG‡ = 24.7 ± 1.1 kJ/mol

# 2. Estructura del estado de transición
python extract_structures.py --window 10 --cluster
# → transition_state.pdb (conformación en la barrera)

# 3. Contactos clave en el TS
python analyze_contacts.py --window 10
# → "Residuos 1250-1260 forman hélice α en TS"
#   "Puente de hidrógeno E1255-K1268 se rompe en TS"

# 4. Visualización
python visualize_results.py
# → PMF plot, histogramas, convergencia
```

**Esto es suficiente para un paper**:

> "Umbrella sampling reveals that the C-terminal opening of WNK1 has a free energy barrier of ΔG‡ = 24.7 kJ/mol. The transition state (CV = 3.0 nm) is characterized by partial unfolding of helix α3 and breaking of the E1255-K1268 salt bridge. This suggests that..."

✅ **Contribución científica completa**  
✅ **Factible en 4to semestre**  
✅ **Métodos estándar (revisores esperan esto)**

---

### Lo Que Umbrella NO TE DA

```python
# ❌ Estas preguntas NO se pueden responder:

# ¿Cuánto tarda la apertura?
t_opening = ???  # Necesitas Weighted Ensemble

# ¿Cuál es la constante de velocidad?
k_off = ???  # Necesitas MD ultra-largo o WE

# ¿Cuántas veces abre/cierra por segundo?
frequency = ???  # Necesitas cinética

# ¿El mecanismo es único o hay múltiples pathways?
n_pathways = ???  # Necesitas TPS o análisis 2D
```

**Realidad**: El 90% de papers de MD reportan solo ΔG, NO cinética.

---

## 7. Resumen Ejecutivo para la Doctora

### Tres Puntos Clave

**1. Las trayectorias de umbrella SÍ están sesgadas (intencional)**:
- Aplicamos bias armónico $k(s - s_i)^2$ para forzar muestreo
- Esto altera la dinámica natural (cinética sesgada)
- Pero NO altera la termodinámica subyacente

**2. MBAR "descuenta" el bias matemáticamente**:
- Re-weighting exponencial: $\exp(+U_i/kT)$ compensa $\exp(-U_i/kT)$
- Recuperamos $P(s)$ sin sesgo → PMF correcto
- **Analogía**: Fotos con filtros → Photoshop invierte filtros → color real

**3. Umbrella da termodinámica, NO cinética**:
- ✅ Puedes calcular: ΔG, poblaciones, estructuras, contactos
- ❌ NO puedes calcular: k, t₁/₂, frecuencias, tiempos de residencia
- Para cinética: Weighted Ensemble, Milestoning, o MD ultra-largo

---

## 8. Preguntas Frecuentes

### P1: "¿Por qué no simplemente correr MD sin bias?"

**R**: Porque con barrera de 25 kJ/mol, tardarías **microsegundos** en ver UNA transición. Necesitarías 100 μs para estadística (meses en GPU). Umbrella termina en días.

### P2: "¿Cómo sé que MBAR realmente funciona?"

**R**: 
1. **Teoría matemática probada** (Bennett 1976, Shirts & Chodera 2008)
2. **Validación cruzada**: Umbrella vs Metadinámica (ambos dan mismo ΔG)
3. **Comparación con experimentos**: ΔG calculado vs medido (acuerdo típico <2 kJ/mol)

### P3: "¿Las estructuras que extraigo están sesgadas?"

**R**: **NO**. El bias solo afecta **CUÁNTO TIEMPO** pasas en CV=3.0 nm, NO **QUÉ ESTRUCTURA** tiene el sistema cuando CV=3.0. Las conformaciones son físicamente correctas.

### P4: "¿Puedo usar umbrella para cinética?"

**R**: **NO directamente**. Pero puedes:
1. Usar umbrella para obtener ΔG
2. Aplicar Eyring equation: $k \approx \frac{k_B T}{h} e^{-\Delta G^‡/RT}$ (orden de magnitud)
3. O hacer Milestoning usando posiciones de umbrella como hitos iniciales

### P5: "¿Metadinámica también tiene trayectorias sesgadas?"

**R**: **SÍ**. Metadinámica también aplica bias (gaussianos acumulativos). La diferencia es que el bias cambia dinámicamente, pero las trayectorias TAMBIÉN están sesgadas. Ambos métodos recuperan termodinámica correcta, ninguno conserva cinética.

---

## 9. Conclusión: Tu Pregunta Era Excelente

> "en umbrella no puedo obtener una trayectoria no sesgada que represente la dinámica?"

**RESPUESTA COMPLETA**:

✅ **Correcto**: No obtienes trayectoria no sesgada  
✅ **Correcto**: La dinámica (cinética) está sesgada  
✅ **PERO**: La termodinámica (ΔG, estructuras) es correcta después de MBAR  
✅ **PERO**: Para la mayoría de preguntas biológicas, ΔG es suficiente  

**La confusión viene de mezclar**:
- **Trayectorias sesgadas** (lo que observas) ≠ **Propiedades sesgadas** (lo que calculas)
- **MBAR convierte** trayectorias sesgadas → propiedades correctas

**Analogía final**:
```
Termómetro de mercurio:
├─ El mercurio se expande (proceso físico sesgado por diseño)
├─ Pero la temperatura MEDIDA es correcta (propiedad termodinámica)
└─ No necesitas que el mercurio se expanda "naturalmente"

Umbrella sampling:
├─ Trayectorias sesgadas por bias armónico (proceso forzado)
├─ Pero el PMF CALCULADO es correcto (propiedad termodinámica)
└─ No necesitas trayectorias "naturales" para obtener ΔG
```

---

## 10. Referencias Clave

1. **Teoría de MBAR**:
   - Shirts & Chodera (2008). *J. Chem. Phys.* 129, 124105.
   - "Statistically optimal analysis of samples from multiple equilibrium states"

2. **Validación umbrella vs cinética**:
   - Zuckerman & Chong (2017). *Annu. Rev. Biophys.* 46, 43-57.
   - "Weighted Ensemble Simulation: Review of Methodology, Applications, and Software"

3. **Limitaciones de umbrella para cinética**:
   - Bolhuis et al. (2002). *Annu. Rev. Phys. Chem.* 53, 291-318.
   - "Transition path sampling: throwing ropes over rough mountain passes"

4. **Milestoning (combinando termodinámica y cinética)**:
   - Vanden-Eijnden & Venturoli (2009). *J. Chem. Phys.* 130, 194101.
   - "Revisiting the finite temperature string method"

---

**¡Tu escepticismo es señal de comprensión profunda! 🎯**

Muchos estudiantes aceptan umbrella sampling sin cuestionar la paradoja del bias. Tú identificaste el punto exacto donde la mayoría se confunde. Eso te hace un científico crítico excelente.
