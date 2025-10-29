# 🧬 Umbrella Sampling y la Quinasa WNK: De los Fundamentos Físicos a la Implementación Computacional

**Una Guía Exhaustiva para Investigadores sin Conocimientos Previos de Python**

---

## 📋 Tabla de Contenidos

### PARTE I: FUNDAMENTOS TEÓRICOS
1. [¿Qué es la Energía Libre y Por Qué Importa?](#1-qué-es-la-energía-libre-y-por-qué-importa)
2. [El Problema: Eventos Raros en Biología Molecular](#2-el-problema-eventos-raros-en-biología-molecular)
3. [La Solución: Umbrella Sampling](#3-la-solución-umbrella-sampling)
4. [WHAM: Reconstruyendo el Panorama Energético Real](#4-wham-reconstruyendo-el-panorama-energético-real)
5. [Coordenadas Colectivas: Eligiendo Qué Medir](#5-coordenadas-colectivas-eligiendo-qué-medir)

### PARTE II: IMPLEMENTACIÓN TÉCNICA
6. [Arquitectura del Sistema: Un Recorrido Visual](#6-arquitectura-del-sistema-un-recorrido-visual)
7. [El Motor: umbrella_sampling_calculator.py](#7-el-motor-umbrella_sampling_calculatorpy)
8. [El Orquestador: umbrella_suite/](#8-el-orquestador-umbrella_suite)
9. [Visualización y Análisis](#9-visualización-y-análisis)
10. [Scripts de Lanzamiento Multiplataforma](#10-scripts-de-lanzamiento-multiplataforma)

### PARTE III: APLICACIÓN A WNK
11. [La Quinasa WNK: Un Sensor Molecular Extraordinario](#11-la-quinasa-wnk-un-sensor-molecular-extraordinario)
12. [El Motivo RFxV y el Dominio CCT: Una Historia de Reconocimiento Molecular](#12-el-motivo-rfxv-y-el-dominio-cct)
13. [Nuestro Objetivo: Mapear el Panorama Termodinámico de WNK](#13-nuestro-objetivo-mapear-el-panorama-termodinámico-de-wnk)
14. [Del Demo Sintético a las Simulaciones Reales](#14-del-demo-sintético-a-las-simulaciones-reales)

### APÉNDICES
- [A. Glosario de Términos](#apéndice-a-glosario-de-términos)
- [B. Referencias Bibliográficas](#apéndice-b-referencias-bibliográficas)
- [C. Recursos Computacionales Requeridos](#apéndice-c-recursos-computacionales-requeridos)

---

## PARTE I: FUNDAMENTOS TEÓRICOS

### 1. ¿Qué es la Energía Libre y Por Qué Importa?

#### 1.1 La Termodinámica en 5 Palabras (para Doctores)

Imagina que eres un director de orquesta y tu partitura tiene dos elementos:

1. **Energía (Entalpía, ΔH)**: ¿Qué tan "cómodos" están los músicos con sus asientos? (¿las moléculas están en configuraciones estables?)
2. **Desorden (Entropía, ΔS)**: ¿Cuántas formas diferentes pueden sentarse y seguir tocando bien? (¿cuántas configuraciones equivalentes existen?)

La **energía libre de Gibbs (ΔG)** es la combinación de ambas:

```
ΔG = ΔH - T·ΔS
```

Donde:
- **ΔG < 0**: El proceso ocurre espontáneamente (como una pelota rodando cuesta abajo)
- **ΔG > 0**: Necesitas empujar (energía externa) para que ocurra
- **T**: Temperatura (cuánto "agita" el sistema térmicamente)

**¿Por qué importa en biología?**

Las proteínas no son estructuras rígidas. Son como edificios que constantemente tiemblan, se flexionan y cambian de forma. Algunas de estas formas son:
- **Activas**: La proteína hace su trabajo (ej. fosforila otra proteína)
- **Inactivas**: La proteína está "apagada"

El **paisaje de energía libre** es un mapa que muestra:
- Qué formas son más probables (valles = bajo ΔG = estables)
- Qué tan difícil es pasar de una forma a otra (montañas = alto ΔG = barreras)

```
                    PAISAJE DE ENERGÍA LIBRE
                    
Energía Libre (kcal/mol)
    ↑
 30 │                      ╱╲               
    │                     ╱  ╲              ← Barrera (Estado de Transición)
 20 │                    ╱    ╲             
    │                   ╱      ╲            
 10 │      ╱╲          ╱        ╲___        
    │     ╱  ╲        ╱             ╲       
  0 │____╱    ╲______╱               ╲____  
    │    A      B           C            D
    └────────────────────────────────────→ Coordenada de Reacción
    
    A: Estado Inactivo (dímero)
    B: Barrera conformacional pequeña  
    C: Estado de Transición (alta energía)
    D: Estado Activo (monómero)
```

**Interpretación**:
- WNK pasa el 90% del tiempo en **A** (valle profundo)
- Cruzar a **D** requiere ~25 kcal/mol (raro, ocurre pocas veces por segundo)
- Si bajamos la barrera (ej. con osmolitos) → más activación

---

#### 1.2 El Potencial de Fuerza Media (PMF): El Mapa Termodinámico

El **PMF** es el equivalente a un mapa topográfico de montañas y valles, pero para energía:

**Definición matemática** (no te asustes):
```
PMF(ξ) = -kT ln⟨P(ξ)⟩
```

**Traducción para humanos**:
- **ξ** (xi): Es la coordenada que medimos (ej. distancia entre dos dominios)
- **P(ξ)**: Probabilidad de observar el sistema en esa coordenada
- **ln⟨P(ξ)⟩**: Logaritmo natural de la probabilidad
- **kT**: Constante de Boltzmann × Temperatura (energía térmica disponible)

**En palabras simples**:
> *"El PMF te dice cuánta energía libre cuesta mover el sistema a cada posición de la coordenada"*

**Ejemplo concreto con WNK**:
```
Distancia entre dominios (Å)  |  PMF (kcal/mol)  |  Interpretación
─────────────────────────────────────────────────────────────────
        8.0                    |       15.2       |  Muy cerca → comprimido → inestable
       10.0                    |        0.0       |  ← Distancia ideal (mínimo de energía)
       12.0                    |       12.8       |  Alejándose → rompiendo contactos
       14.0                    |       25.3       |  Muy lejos → completamente separado
```

Si trazas estos puntos en una gráfica, obtienes el perfil del PMF que viste arriba.

---

#### 1.3 Conceptos Clave: Entalpía vs. Entropía

##### Caso 1: Proceso Impulsado por Entalpía (ΔH domina)
```
Ejemplo: Hielo derritiéndose
    
H₂O(sólido) → H₂O(líquido)

ΔH = +6 kcal/mol  (necesitas romper enlaces → endotérmico)
ΔS = +22 cal/(mol·K)  (el líquido es más desordenado → favorable)

A T = 0°C (273 K):
ΔG = ΔH - T·ΔS = 6 - 273×(0.022) = 6 - 6 = 0 kcal/mol
             ↑                 ↑
          Desfavorable      Favorable
          
→ Equilibrio: hielo y agua coexisten
```

##### Caso 2: Proceso Impulsado por Entropía (ΔS domina)
```
Ejemplo: WNK1 dímero disociándose

                  ┌─────────────┐
    Dímero        │   H₂O H₂O   │  ← Agua ordenada en la interfaz
   (Inactivo)     │   H₂O H₂O   │
                  └─────────────┘
                        ↓
                   Osmolitos
                   (PEG, sacarosa)
                        ↓
    Monómero         [H₂O]  [H₂O]   ← Agua libre (desordenada)
    (Activo)         [H₂O]  [H₂O]

ΔH = +5 kcal/mol  (romper contactos proteína-proteína → desfavorable)
ΔS = +50 cal/(mol·K)  (liberar ~200 moléculas de agua → MUY favorable)

A T = 300 K:
ΔG = 5 - 300×(0.050) = 5 - 15 = -10 kcal/mol
                                  ↑
                              ¡Espontáneo!
```

**Mensaje clave**: WNK se activa porque liberar agua es termodinámicamente irresistible, no porque los monómeros sean más estables intrínsecamente.

---

### 2. El Problema: Eventos Raros en Biología Molecular

#### 2.1 La Barrera del Tiempo en Simulaciones

**Dinámicas Moleculares Clásicas (MD sin bias)**:
```
Paso de tiempo: 2 femtosegundos (2×10⁻¹⁵ s)
Simulación típica: 100 nanosegundos (10⁻⁷ s)
Pasos totales: 50,000,000 iteraciones
```

**El problema**:
| Evento Biológico | Escala de Tiempo | Accesible con MD Clásica |
|---|---|---|
| Vibración de enlace | femtosegundos | ✅ Sí |
| Fluctuación de cadena lateral | picosegundos | ✅ Sí |
| Movimiento de bucles | nanosegundos | ✅ Sí (apenas) |
| Apertura/cierre de dominios | microsegundos | ❌ No |
| Plegamiento de proteínas | milisegundos | ❌ No |
| Conformación activa ↔ inactiva WNK | microsegundos-milisegundos | ❌ No |

**Analogía del Túnel**:
```
Imagina que quieres medir cuántos autos cruzan un túnel de montaña cada día.

Método 1 (MD sin bias): 
    Parate en la entrada y cuenta los autos que pasan.
    Problema: Si solo 1 auto cruza cada hora, necesitarás esperar 24 horas 
              para obtener una muestra decente.
              
    En MD: Si el evento ocurre cada 10 µs, necesitarías simular ¡meses de 
           tiempo real de computadora!
```

#### 2.2 El Paisaje Energético y las Barreras

La razón por la cual estos eventos son raros es que están separados por **barreras de energía libre**:

```
         ╱‾‾‾‾‾╲
        ╱       ╲
       ╱  25 kcal ╲     ← Barrera
      ╱   /mol    ╲
─────╱             ╲────
Estado A        Estado B
(Inactivo)       (Activo)
```

**Probabilidad de cruzar la barrera** (ecuación de Boltzmann):
```
P ∝ exp(-ΔG‡/kT)

Donde:
- ΔG‡: Altura de la barrera
- kT ≈ 0.6 kcal/mol a 300 K

Ejemplo:
ΔG‡ = 15 kcal/mol → P ∝ exp(-15/0.6) = exp(-25) ≈ 10⁻¹¹
                                             ↑
                               ¡1 en 100 mil millones de intentos!
```

**Frecuencia de intentos** en una proteína:
- Fluctuaciones térmicas: ~10¹² intentos/segundo
- Tasa real de cruce: 10¹² × 10⁻¹¹ = **10 eventos/segundo**

Para observar 100 eventos (estadística decente):
```
Tiempo requerido = 100 eventos / 10 s⁻¹ = 10 segundos
```

¡Pero esto en tiempo de simulación MD = 10⁷ × tiempo real! = **Imposible**.

---

#### 2.3 Referencias Fundamentales

Los conceptos aquí presentados están basados en trabajos seminales:

**[1] Kumar S, Rosenberg JM, Bouzida D, Swendsen RH, Kollman PA. (1992)**  
*"THE weighted histogram analysis method for free-energy calculations on biomolecules. I. The method"*  
Journal of Computational Chemistry, 13(8):1011-1021.  
doi: 10.1002/jcc.540130812  
**Citaciones: 5,876** (paper fundacional de WHAM)

Este artículo estableció el método WHAM que usaremos para combinar datos de múltiples ventanas de umbrella sampling.

**[2] Oshima H, Re S, Sugita Y. (2019)**  
*"Replica-exchange umbrella sampling combined with Gaussian accelerated molecular dynamics for free-energy calculation of biomolecules"*  
Journal of Chemical Theory and Computation, 15(10):5199-5208.  
doi: 10.1021/acs.jctc.9b00501  
**Citaciones: 63**

Presenta mejoras modernas a umbrella sampling que aceleran la convergencia, relevantes para nuestro caso de WNK.

---

### 3. La Solución: Umbrella Sampling

#### 3.1 La Idea Central: Forzar el Muestreo con Potenciales de Sesgo

**Analogía del Túnel (continuación)**:
```
Método 2 (Umbrella Sampling):
    No esperes pasivamente. En su lugar:
    
    1. Coloca "estaciones de peaje" a lo largo del túnel
    2. En cada estación, ofreces un descuento para que los autos
       se detengan ahí (potencial de sesgo)
    3. Cuenta cuántos autos se detienen en cada estación
    4. Luego, usa matemáticas para "quitar" el efecto del descuento
       y calcular cuántos habrían pasado naturalmente
```

**En términos de simulación**:
```
Sin Umbrella:
    Sistema explora libremente → se queda atascado en valles profundos
    
Con Umbrella:
    Añades un potencial artificial ("umbrella") que:
    - Empuja el sistema hacia regiones de alta energía
    - Permite muestrear barreras que serían inaccesibles
    - Luego se "resta" matemáticamente para obtener el PMF real
```

---

#### 3.2 El Potencial de Umbrella: Matemáticas Simple

El potencial de umbrella es típicamente un **resorte armónico**:

```
V_umbrella(ξ) = ½ k (ξ - ξ₀)²
```

**Componentes**:
- **ξ**: Coordenada colectiva actual (ej. distancia = 10.5 Å)
- **ξ₀**: Centro de la ventana (ej. distancia objetivo = 10.0 Å)
- **k**: Constante de fuerza del resorte (ej. 12 kcal/mol/Å²)

**Interpretación física**:
```
Imagina un resorte invisible conectado entre dos átomos:

    Átomo A ═══════⊙~~~~~~~~⊙═══════ Átomo B
                   resorte
                   
- Si los átomos se alejan más de ξ₀ → el resorte tira (penalidad energética)
- Si se acercan más de ξ₀ → el resorte empuja (penalidad energética)
- El sistema prefiere estar cerca de ξ₀
```

**Ejemplo numérico**:
```
Ventana centrada en ξ₀ = 10.0 Å con k = 12 kcal/mol/Å²

Si el sistema está en ξ = 11.0 Å:
    V_umbrella = ½ × 12 × (11.0 - 10.0)²
               = 6 × 1.0²
               = 6.0 kcal/mol
               
Si el sistema está en ξ = 12.0 Å:
    V_umbrella = ½ × 12 × (12.0 - 10.0)²
               = 6 × 4.0
               = 24.0 kcal/mol  ← Penalidad fuerte!
```

El sistema tenderá a fluctuar alrededor de ξ₀ = 10.0 Å, pero explorará ~9.5-10.5 Å debido a las fluctuaciones térmicas.

---

#### 3.3 Sistema de Ventanas Superpuestas

**Clave del éxito**: Necesitas **muchas ventanas** que se **solapen**.

```
EJEMPLO: Medir PMF desde 8.0 Å hasta 14.0 Å

Configuración de 6 ventanas:
    
Energía del Umbrella
     ↑
  30 │    ╱‾╲         ╱‾╲         ╱‾╲         ╱‾╲         ╱‾╲         ╱‾╲
     │   ╱   ╲       ╱   ╲       ╱   ╲       ╱   ╲       ╱   ╲       ╱   ╲
  20 │  ╱     ╲     ╱     ╲     ╱     ╲     ╱     ╲     ╱     ╲     ╱     ╲
     │ ╱       ╲   ╱       ╲   ╱       ╲   ╱       ╲   ╱       ╲   ╱       ╲
  10 │╱         ╲ ╱         ╲ ╱         ╲ ╱         ╲ ╱         ╲ ╱         ╲
     │           X           X           X           X           X
   0 └──────────────────────────────────────────────────────────────────────→
     8.0       9.2       10.4      11.6      12.8      14.0     Distancia (Å)
     W1        W2        W3        W4        W5        W6
```

**Histogramas generados por cada ventana**:
```
Ventana 1 (centrada en 8.0 Å):
Frecuencia
    │     ╱‾╲
    │    ╱   ╲
    │   ╱     ╲___
    └──────────────→ 7.5  8.0  8.5  9.0  9.5  Distancia
    
Ventana 2 (centrada en 9.2 Å):
Frecuencia
    │              ╱‾╲
    │         ____╱   ╲
    │        ╱         ╲___
    └──────────────→ 8.5  9.0  9.5 10.0 10.5

...y así sucesivamente
```

**El solapamiento es crítico**:
```
Ventana 1 samplea: [7.5 - 9.0 Å]    ─────────
Ventana 2 samplea: [8.5 - 10.0 Å]         ─────────
                                   ↑
                              Overlap necesario!
```

Sin overlap, WHAM no puede "coser" las ventanas para formar un PMF continuo.

---

#### 3.4 Flujo de Trabajo Completo

```
PASO 1: Definir la Coordenada Colectiva (CV)
    ┌─────────────────────────────────────┐
    │ ¿Qué quieres medir?                  │
    │                                      │
    │ Opciones comunes:                   │
    │ • Distancia entre dos átomos/grupos │
    │ • Ángulo diedro                     │
    │ • RMSD desde una estructura         │
    │ • Número de contactos nativos       │
    └─────────────────────────────────────┘
                 ↓
PASO 2: Elegir Centros de Ventanas
    ┌─────────────────────────────────────┐
    │ Rango: [ξ_min, ξ_max]               │
    │ Número de ventanas: 30-60           │
    │ Espaciado: ~0.5-1.0 Å para distancias│
    └─────────────────────────────────────┘
                 ↓
PASO 3: Generar Configuraciones Iniciales
    ┌─────────────────────────────────────┐
    │ Método A: Steered MD (SMD)          │
    │   - Tira del sistema de A → B       │
    │   - Extrae snapshots                │
    │                                      │
    │ Método B: Interpolación             │
    │   - Mezcla geometrías de A y B      │
    └─────────────────────────────────────┘
                 ↓
PASO 4: Correr Simulaciones de Umbrella
    ┌─────────────────────────────────────┐
    │ Para cada ventana i:                │
    │   1. Coloca sistema en ξ₀ᵢ          │
    │   2. Añade V_umbrella = ½k(ξ-ξ₀ᵢ)²  │
    │   3. Equilibra 1-5 ns               │
    │   4. Produce 20-50 ns               │
    │   5. Guarda trayectoria             │
    └─────────────────────────────────────┘
                 ↓
PASO 5: Análisis con WHAM
    ┌─────────────────────────────────────┐
    │ Input: Histogramas de todas ventanas│
    │ Output: PMF(ξ)                      │
    │                                      │
    │ PMF te da:                          │
    │ • Energía libre relativa            │
    │ • Altura de barreras                │
    │ • Posición de estados estables      │
    └─────────────────────────────────────┘
```

---

### 4. WHAM: Reconstruyendo el Panorama Energético Real

#### 4.1 El Problema de Combinar Histogramas Sesgados

Después de correr todas las ventanas, tienes:
```
Ventana 1: Histograma H₁(ξ) con potencial V₁(ξ) = ½k(ξ - ξ₀¹)²
Ventana 2: Histograma H₂(ξ) con potencial V₂(ξ) = ½k(ξ - ξ₀²)²
...
Ventana N: Histograma Hₙ(ξ) con potencial Vₙ(ξ) = ½k(ξ - ξ₀ⁿ)²
```

**El desafío**: Cada histograma está "deformado" por su umbrella. ¿Cómo combinarlos para obtener el PMF **sin sesgo**?

**Ejemplo visual**:
```
Datos Crudos de Umbrella (SESGADOS):
    
Ventana en ξ₀ = 8 Å:
    Frecuencia
         │  ╱‾╲          ← Pico artificial debido al resorte
         │ ╱   ╲
         └───────→ 7  8  9  10 Å
         
Ventana en ξ₀ = 10 Å:
    Frecuencia
         │        ╱‾╲    ← Pico en otro lugar
         │       ╱   ╲
         └───────→ 9  10  11  12 Å

Pregunta: ¿Cuál es la distribución REAL sin los resortes?
```

---

#### 4.2 La Ecuación WHAM (sin pánico)

El método WHAM resuelve iterativamente dos ecuaciones acopladas:

**Ecuación 1: Probabilidad sin sesgo**
```
P_unbiased(ξ) ∝ N(ξ) / Σᵢ nᵢ exp[-βVᵢ(ξ) + βFᵢ]
```

**Ecuación 2: Energías libres de las ventanas**
```
exp[-βFᵢ] = Σ_ξ P_unbiased(ξ) exp[-βVᵢ(ξ)]
```

**Traducción a español**:
- **P_unbiased(ξ)**: Probabilidad verdadera de encontrar el sistema en ξ (sin resortes)
- **N(ξ)**: Número total de veces que se observó ξ en TODAS las ventanas
- **nᵢ**: Número de observaciones en la ventana i
- **Vᵢ(ξ)**: Potencial del umbrella en la ventana i
- **Fᵢ**: Energía libre del umbrella i (constante de normalización)
- **β = 1/(kT)**: Inverso de la energía térmica

**Algoritmo iterativo**:
```
1. Inicia con un guess para Fᵢ (ej. todos = 0)

2. REPITE hasta convergencia:
    a) Calcula P_unbiased(ξ) usando Ecuación 1 y los Fᵢ actuales
    b) Calcula nuevos Fᵢ usando Ecuación 2 y P_unbiased actual
    c) Checa si los valores cambiaron < tolerancia (ej. 0.0001)
    
3. Una vez convergido:
    PMF(ξ) = -kT ln[P_unbiased(ξ)] + constante
```

**Intuición**:
> WHAM es como resolver un rompecabezas donde cada pieza (ventana) está ligeramente deformada. El algoritmo busca la deformación que, al deshacerla, hace que todas las piezas encajen perfectamente.

---

#### 4.3 MBAR: La Evolución de WHAM

**MBAR** (Multistate Bennett Acceptance Ratio) es una generalización de WHAM desarrollada en 2008:

**Diferencias clave**:
| Aspecto | WHAM | MBAR |
|---|---|---|
| **Datos requeridos** | Histogramas de ξ | Configuraciones completas |
| **Precisión** | Buena | Óptima (más eficiente estadísticamente) |
| **Errores** | Estimación bootstrap | Analíticos (más confiables) |
| **Aplicabilidad** | 1D (una CV) principalmente | Multidimensional (varias CVs) |

**Cuándo usar cada uno**:
- **WHAM**: Tu análisis inicial, debugging, visualización rápida
- **MBAR**: Paper final, cálculos precisos de ΔΔG, comparación con experimento

Nuestro código soporta ambos métodos vía la librería `pymbar`.

---

#### 4.4 Convergencia y Validación

**¿Cómo sabes que tus resultados son confiables?**

##### Criterio 1: Solapamiento de Histogramas
```
Bueno:
    Ventana i:    ────────
    Ventana i+1:       ────────
                      ↑
                  Overlap > 20%
                  
Malo:
    Ventana i:    ────────
    Ventana i+1:              ────────
                              ↑
                      Gap sin sampleo!
```

##### Criterio 2: Convergencia Temporal
```
Divide tu simulación en bloques:
    
    Bloque 1: 0-10 ns   →  PMF₁(ξ)
    Bloque 2: 10-20 ns  →  PMF₂(ξ)
    Bloque 3: 20-30 ns  →  PMF₃(ξ)
    
    Calcula diferencia:
    RMSD(PMF₁, PMF₂) < 1 kcal/mol  →  ✅ Converged
    RMSD(PMF₂, PMF₃) < 0.5 kcal/mol  →  ✅ Bien converged
```

##### Criterio 3: Bootstrapping para Errores
```
WHAM/MBAR pueden estimar incertidumbre:
    
    PMF(ξ=10 Å) = 5.2 ± 0.8 kcal/mol
                       ↑
                  Error estándar
    
Si los errores son > 2 kcal/mol en regiones importantes → simula más tiempo
```

**Referencias para WHAM/MBAR**:

**[3] Shirts MR, Chodera JD. (2008)**  
*"Statistically optimal analysis of samples from multiple equilibrium states"*  
Journal of Chemical Physics, 129:124105.  
doi: 10.1063/1.2978177  
**Citaciones: 1,847** (paper fundacional de MBAR)

---

### 5. Coordenadas Colectivas: Eligiendo Qué Medir

#### 5.1 El Arte de Seleccionar CVs

**Pregunta central**: ¿Qué propiedad del sistema captura el proceso biológico de interés?

**Mala CV → Simulación inútil**  
**Buena CV → Descubrimiento científico**

```
ANALOGÍA: Estudiar tráfico de una ciudad

Mala CV: Temperatura del asfalto
    - Varía mucho, pero no te dice sobre flujo de autos
    
Mala CV: Color de los autos
    - Interesante, pero irrelevante para congestión
    
Buena CV: Número de autos por kilómetro (densidad)
    - Captura directamente el fenómeno de congestión
    
Buena CV: Velocidad promedio
    - Correlaciona con tráfico fluido vs. atascado
```

---

#### 5.2 Tipos Comunes de CVs en Proteínas

##### Tipo 1: Distancias
```
Definición: Distancia euclidiana entre dos grupos de átomos

    Ejemplo: Distancia entre Cα del residuo 245 y Cα del residuo 292
    
    d = √[(x₂-x₁)² + (y₂-y₁)² + (z₂-z₁)²]
    
Pros:
    ✅ Fácil de interpretar
    ✅ Rápida de calcular
    ✅ Directamente observable en estructuras cristalográficas
    
Contras:
    ❌ No captura rotaciones
    ❌ Puede no reflejar cambios funcionales complejos
    
Casos de uso:
    • Apertura/cierre de dominios
    • Distancia ligando-sitio activo
    • Formación de dímeros
```

**Nuestra aplicación en WNK**:
```
CV = distancia(CA_245, CA_292)
         ↑           ↑
    Residuo β3   Residuo αC
    
Esta distancia cambia dramáticamente entre:
    - Estado inactivo (dímero): ~8 Å
    - Estado activo (monómero): ~14 Å
```

##### Tipo 2: Ángulos Diedros (Dihedral Angles)
```
Definición: Ángulo de torsión entre 4 átomos

        A ─── B
               │
               C ─── D
               
    Ángulo φ = rotación alrededor del eje B-C
    
Ejemplos en proteínas:
    • Phi (φ): Ángulo backbone N-Cα
    • Psi (ψ): Ángulo backbone Cα-C
    • Chi (χ): Ángulo de cadena lateral
    
Casos de uso:
    • Transiciones de estructura secundaria (hélice ↔ hoja)
    • Rotámeros de cadenas laterales (Arg, Lys)
```

##### Tipo 3: RMSD (Root Mean Square Deviation)
```
Definición: Desviación promedio de la estructura actual vs. una referencia

    RMSD = √[Σᵢ (rᵢ - rᵢ_ref)² / N]
    
Donde:
    • rᵢ: Posición del átomo i en la estructura actual
    • rᵢ_ref: Posición en la estructura de referencia
    • N: Número de átomos
    
Interpretación:
    RMSD = 1 Å  →  Muy similar a referencia
    RMSD = 3 Å  →  Cambio conformacional moderado
    RMSD = 8 Å  →  Conformación completamente diferente
    
Casos de uso:
    • Plegamiento de proteínas (referencia = estado nativo)
    • Transiciones entre estados cristalográficos conocidos
```

##### Tipo 4: Número de Contactos
```
Definición: Cuenta cuántas parejas de átomos están dentro de un umbral

    N_contacts = Σ_ij Θ(r_cutoff - d_ij)
    
    Donde Θ es función escalón:
        Θ(x) = 1 si x > 0
        Θ(x) = 0 si x < 0
    
Ejemplo:
    r_cutoff = 4.5 Å
    
    Conformación A: 45 contactos  ←  Compacta
    Conformación B: 12 contactos  ←  Extendida
    
Casos de uso:
    • Colapso hidrofóbico en plegamiento
    • Estabilidad de interfaces proteína-proteína
```

##### Tipo 5: CVs Avanzadas (Path Collective Variables)
```
Para procesos complejos donde no hay una CV simple:

Path CV: Progreso a lo largo de una trayectoria de transición predefinida

    Ejemplo: Plegamiento de un péptido
        Estado A (desplega) ──[ruta]──> Estado B (plegado)
        
        s(t) ∈ [0, 1]
        s = 0: Completamente desplegado
        s = 1: Completamente plegado
        
    OpenMM no tiene Path CVs nativas (necesitas PLUMED)
```

---

#### 5.3 Validación de tu CV: Proyección de la Transición

**Test ácido**: Si corres una MD libre que atraviesa la transición, ¿tu CV lo detecta?

```
Experimento mental:

1. Simula WNK dímero disociándose naturalmente (tardaría microsegundos)

2. Proyecta la trayectoria sobre tu CV:

    ┌────────────────────────────────────────┐
    │ CV (distancia β3-αC)                   │
    │   ↑                                    │
    │14 │                        ┌───────────│  ← Monómero
    │   │                    ___╱            │
    │12 │                ___╱                │
    │   │            ___╱                    │
    │10 │        ___╱                        │
    │   │    ___╱                            │
    │ 8 │───╱                                │  ← Dímero
    │   └────────────────────────────────────→
    │   0 µs    1 µs     2 µs     3 µs  Tiempo
    └────────────────────────────────────────┘
    
✅ Buena CV: Cambio monotónico y suave
❌ Mala CV: Fluctúa aleatoriamente sin tendencia clara
```

**Criterios de una buena CV**:
1. **Separa estados**: Valores diferentes para A y B
2. **Continua**: Transición suave, no jumps abruptos
3. **Interpretable**: Conexión clara con mecanismo molecular
4. **Computacionalmente barata**: No ralentiza la simulación

---

#### 5.4 CVs para WNK: Nuestras Decisiones

**Objetivo**: Capturar el cambio conformacional asociado con la transición dímero → monómero y la activación de WNK.

**CV principal**:
```python
CV = distance(CA_residue_245, CA_residue_292)
```

**Justificación**:
1. **β3-strand (residuo 245)** y **αC-helix (residuo 292)** son elementos de estructura secundaria críticos en el lóbulo N-terminal.

2. En estructuras cristalográficas:
   ```
   Dímero inactivo (PDB: 6CN9): d ≈ 8.2 Å
   Monómero activo (inferido):  d ≈ 13.5 Å
   ```

3. Esta distancia correlaciona con:
   - **Compactación del sitio activo**
   - **Accesibilidad del sitio de fosforilación** (loop de activación)
   - **Posición relativa del C-helix** (marca de activación en kinasas)

**Diagrama estructural**:
```
Vista lateral de WNK1 dominio kinasa:

        Lóbulo N                    Lóbulo C
         ╱╲                           ╱╲
        ╱  ╲                         ╱  ╲
       ╱αC  ╲                       ╱    ╲
      ╱  ·   ╲                     ╱      ╲
     ╱   ·    ╲___________________╱        ╲
    ╱    ·   β3                            ╲
   ╱─────●────→                             ╲
  ╱ 245  ↓  292                              ╲
 ╱       d = CV                               ╲
╱                                              ╲

Cuando d aumenta:
    - αC se aleja de β3
    - Loop de activación se libera
    - ATP se posiciona correctamente
    → ACTIVACIÓN
```

**Referencias sobre CVs**:

**[4] Thiede EH, Van Koten B, Weare J, Dinner AR. (2016)**  
*"Eigenvector method for umbrella sampling enables error analysis"*  
arXiv:1607.03722  
(Método para validar CVs usando eigenvectores de matrices de transición)

**[5] Awasthi S, Nair NN. (2015)**  
*"Exploring high-dimensional free energy landscapes: Temperature accelerated sliced sampling"*  
arXiv:1508.05181  
(Técnicas para manejar múltiples CVs simultáneamente)

---

**FIN DE LA PARTE I: FUNDAMENTOS TEÓRICOS**

---

## PARTE II: IMPLEMENTACIÓN TÉCNICA

### 6. Arquitectura del Sistema: Un Recorrido Visual

#### 6.1 Vista de 10,000 Pies: Los Componentes

Nuestro sistema está organizado como una fábrica de análisis termodinámica:

```
┌──────────────────────────────────────────────────────────────────┐
│                     UMBRELLA SAMPLING SYSTEM                      │
└──────────────────────────────────────────────────────────────────┘
           │
           ├─── ① INPUT LAYER (Configuración)
           │    └─ umbrella_suite/config.py
           │       • UmbrellaPipelineConfig: Parámetros globales
           │       • UmbrellaWindowConfig: Config por ventana
           │
           ├─── ② ENGINE LAYER (Motor de Simulación)
           │    └─ Chronosfold/umbrella_sampling_calculator.py
           │       • run_full_umbrella_sampling()
           │       • create_umbrella_potential()
           │       • Interfaz con OpenMM
           │
           ├─── ③ ORCHESTRATION LAYER (Coordinación)
           │    └─ umbrella_suite/pipeline.py
           │       • UmbrellaSamplingPipeline class
           │       • Maneja batch de ventanas
           │       • Exporta resultados (.dat files)
           │
           ├─── ④ ANALYSIS LAYER (Post-Procesamiento)
           │    └─ umbrella_suite/analysis.py
           │       • compute_pmf() → WHAM/MBAR
           │       • generate_synthetic_windows() → Fallback
           │
           ├─── ⑤ VISUALIZATION LAYER (Gráficas)
           │    └─ umbrella_suite/visualization.py
           │       • plot_umbrella_diagnostics()
           │       • 3-panel figure system
           │
           └─── ⑥ CLI LAYER (Interfaz de Usuario)
                ├─ run_wnk_pipeline.py (Python CLI)
                ├─ run_umbrella.ps1 (Windows launcher)
                └─ run_umbrella.sh (Linux launcher)
```

**Flujo de datos típico**:
```
Usuario ejecuta script
        ↓
run_umbrella.ps1 detecta Conda
        ↓
Activa environment bsm-lancad-env
        ↓
Llama run_wnk_pipeline.py --synthetic
        ↓
Pipeline crea UmbrellaPipelineConfig
        ↓
├─ Modo REAL: Llama umbrella_sampling_calculator
│               ↓
│            OpenMM simula cada ventana
│               ↓
│            Genera trayectorias .dcd
│
└─ Modo SYNTHETIC: generate_synthetic_windows()
                ↓
             Crea datos gaussianos
                ↓
         UmbrellaWindow objects
                ↓
          Exporta .dat files
                ↓
          compute_pmf() → PMF array
                ↓
          plot_umbrella_diagnostics()
                ↓
          Guarda umbrella_diagnostics.png
                ↓
          ✅ Resultados en umbrella_results/
```

---

#### 6.2 Anatomía de un Objeto de Datos: Sin Miedo a las Clases

**Pregunta**: ¿Qué es una "clase" en Python? ¿Por qué usar `dataclass`?

**Analogía: Formulario vs. Formulario Impreso**
```
Clase = Template de formulario en blanco
    
    ┌─────────────────────────┐
    │  FORMULARIO DE PEDIDO   │
    │                         │
    │ Nombre: ____________    │
    │ Cantidad: __________    │
    │ Precio: ____________    │
    └─────────────────────────┘
    
Instancia = Formulario llenado
    
    ┌─────────────────────────┐
    │  FORMULARIO DE PEDIDO   │
    │                         │
    │ Nombre: "Laptop"        │
    │ Cantidad: 5             │
    │ Precio: 1200.00         │
    └─────────────────────────┘
```

**Nuestro código usa `dataclass`** (desde Python 3.7), que es una forma simplificada de crear clases:

```python
# SIN dataclass (tedioso):
class UmbrellaWindowConfig:
    def __init__(self, center, force_constant, simulation_time_ps):
        self.center = center
        self.force_constant = force_constant
        self.simulation_time_ps = simulation_time_ps
    
    def __repr__(self):
        return f"UmbrellaWindowConfig(center={self.center},...)"
    
    # ...más código boilerplate

# CON dataclass (automático):
@dataclass
class UmbrellaWindowConfig:
    center: float
    force_constant: float
    simulation_time_ps: int
    
    # Python genera __init__, __repr__, __eq__ automáticamente!
```

**Beneficio**: Declaras qué campos necesitas, y Python hace el resto.

---

#### 6.3 Ejemplo Concreto: Crear una Configuración

**Tarea**: Configurar un pipeline para simular 6 ventanas de 10-14 Å, 5 ps cada una.

**Paso a paso (sin código, solo conceptos)**:
```
1. CREAR CONFIG GLOBAL:
    Pipeline Config = Formulario Maestro
        ├─ Lista de centros de ventanas: [10.0, 10.8, 11.6, 12.4, 13.2, 14.0]
        ├─ Constante de fuerza: 12.0 kcal/mol/Å²
        ├─ Tiempo por ventana: 5 ps
        ├─ Temperatura: 300 K
        └─ Archivo de estructura: "wnk_structure.pdb"

2. PIPELINE LEE CONFIG:
    Para cada centro en la lista:
        Crea un "Window Config" individual:
            ├─ Centro: (valor de la lista)
            ├─ Fuerza: (hereda del config global)
            └─ Tiempo: (hereda del config global)

3. PIPELINE EJECUTA:
    Para cada Window Config:
        ├─ Llama al calculator con esos parámetros
        ├─ Calculator devuelve datos de CV vs. tiempo
        └─ Pipeline empaqueta en UmbrellaWindow object

4. PIPELINE EXPORTA:
    Para cada UmbrellaWindow:
        Escribe archivo "window_0.dat":
            8.2  0
            8.3  1
            8.5  2
            ...
            (valores de CV y tiempos)

5. ANÁLISIS LEE ARCHIVOS:
    Carga todos los .dat
        ↓
    Construye histogramas
        ↓
    Llama WHAM
        ↓
    Devuelve PMF(ξ)
```

**Lo importante**: No necesitas saber Python para entender la LÓGICA. Es como seguir una receta de cocina:
- Ingredientes = Config
- Pasos = Pipeline methods
- Plato final = PMF plot

---

### 7. El Motor: `umbrella_sampling_calculator.py`

#### 7.1 Propósito y Responsabilidades

Este archivo es el **único que habla directamente con OpenMM**. Piénsalo como el "conductor de autobús":
- Recibe pasajeros (configuraciones de ventanas)
- Conduce el autobús (corre la simulación MD)
- Deja a los pasajeros en sus destinos (devuelve resultados)

**Responsabilidades**:
1. ✅ Crear el potencial de umbrella (resorte armónico)
2. ✅ Configurar el sistema OpenMM (forcefield, solvente, integradores)
3. ✅ Ejecutar la simulación en lotes (batches) para paralelizar
4. ✅ Guardar trayectorias (opcional, .dcd files)
5. ✅ Devolver datos de la CV vs. tiempo

**NO hace**:
- ❌ Decidir cuántas ventanas correr (eso es del Pipeline)
- ❌ Calcular el PMF (eso es del Analysis)
- ❌ Generar plots (eso es del Visualization)

---

#### 7.2 Función Clave: `create_umbrella_potential()`

**Firma de la función** (no te asustes):
```python
def create_umbrella_potential(
    system: openmm.System,
    cv_atoms: Tuple[int, int],
    k: float,
    r0: float
) -> int:
```

**Traducción**:
- **`system`**: El "universo" de OpenMM donde vive tu proteína (contiene todos los átomos, fuerzas, etc.)
- **`cv_atoms`**: Tupla de 2 números = índices de los átomos entre los que mides distancia
  - Ejemplo: `(2450, 2920)` = átomo 2450 y 2920
- **`k`**: Constante del resorte en kcal/mol/Å²
- **`r0`**: Distancia objetivo (centro de la ventana) en Å
- **`-> int`**: Devuelve un número de identificación de la fuerza (OpenMM lo usa internamente)

**Lo que hace por dentro** (conceptual):
```
1. Crea un objeto "CustomBondForce":
    - Es una "fuerza personalizada" de OpenMM
    - Tú defines la fórmula matemática: "0.5*k*(r-r0)^2"

2. Especifica los parámetros:
    - 'k' = constante de fuerza
    - 'r0' = distancia objetivo

3. Añade un "bond virtual" entre los dos átomos:
    - No es un enlace químico real
    - Solo le dice a OpenMM: "calcula la distancia entre estos dos átomos"

4. Agrega esta fuerza al sistema:
    - system.addForce(umbrella_force)
    
5. Devuelve el ID:
    - Por si después quieres modificar o eliminar esta fuerza
```

**Diagrama de flujo visual**:
```
Input: Átomos A y B, k=12, r0=10.0
            ↓
┌───────────────────────────────────┐
│ Crear CustomBondForce             │
│   Fórmula: "0.5*k*(r-r0)^2"       │
└───────────────────────────────────┘
            ↓
┌───────────────────────────────────┐
│ Definir parámetros:               │
│   k = 12 kcal/mol/Å²              │
│   r0 = 10.0 Å                     │
└───────────────────────────────────┘
            ↓
┌───────────────────────────────────┐
│ Añadir bond entre átomos A y B    │
│   Índices: (2450, 2920)           │
└───────────────────────────────────┘
            ↓
┌───────────────────────────────────┐
│ Agregar fuerza al sistema OpenMM  │
└───────────────────────────────────┘
            ↓
        Listo! ✅
```

---

#### 7.3 Función Clave: `run_full_umbrella_sampling()`

Esta es la función "maestra" que orquesta toda una campaña de umbrella sampling.

**Inputs principales**:
```python
run_full_umbrella_sampling(
    pdb_path: str,               # Ruta al archivo .pdb de tu proteína
    window_centers: List[float], # [8.0, 9.2, 10.4, ...] en Å
    k: float,                    # 12.0 kcal/mol/Å²
    simulation_time_ps: int,     # 5000 ps = 5 ns
    cv_atoms: Tuple[int, int],   # (2450, 2920)
    temperature: float = 300.0,  # Kelvin
    batch_size: int = 4          # Cuántas ventanas simular en paralelo
)
```

**Output**:
```python
{
    'windows': [
        {
            'center': 8.0,
            'cv_values': [8.1, 8.2, 8.0, 7.9, ...],  # Array de distancias
            'time_ps': [0, 2, 4, 6, ...],             # Timestamps
            'k': 12.0
        },
        {
            'center': 9.2,
            ...
        },
        ...
    ],
    'pmf_estimate': [0.0, 1.5, 3.2, ...]  # PMF preliminar (opcional)
}
```

**Proceso interno** (simplificado):
```
PASO 1: PREPARACIÓN
    ├─ Cargar estructura PDB
    ├─ Añadir campo de fuerzas (AMBER, CHARMM, etc.)
    ├─ Añadir solvente (caja de agua)
    └─ Crear integrador (Langevin, para temperatura constante)

PASO 2: CREAR BATCHES
    Si tienes 24 ventanas y batch_size=4:
        Batch 1: ventanas 0-3
        Batch 2: ventanas 4-7
        ...
        Batch 6: ventanas 20-23

PASO 3: PARA CADA BATCH
    ├─ Para cada ventana en el batch:
    │   ├─ Clonar el sistema base
    │   ├─ Añadir umbrella potential con create_umbrella_potential()
    │   ├─ Establecer posiciones iniciales
    │   └─ Crear simulador
    │
    ├─ EJECUTAR EN PARALELO (asyncio):
    │   Ventana 0 corre en Core 0
    │   Ventana 1 corre en Core 1
    │   Ventana 2 corre en Core 2
    │   Ventana 3 corre en Core 3
    │
    ├─ Para cada paso de simulación:
    │   ├─ OpenMM calcula fuerzas (incluyendo umbrella)
    │   ├─ OpenMM integra ecuaciones de movimiento
    │   ├─ OpenMM actualiza posiciones
    │   └─ Guarda distancia CV cada 1 ps
    │
    └─ Cuando termina:
        └─ Devuelve cv_values, time_ps para cada ventana

PASO 4: RETORNAR RESULTADOS
    └─ Empaqueta todos los resultados en un diccionario
```

**Concepto clave: Paralelización**
```
SIN paralelización (serial):
    Ventana 1 → 5 ns → ████████████░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
    Ventana 2 → 5 ns → ░░░░░░░░░░░░████████████░░░░░░░░░░░░░░░░░░
    Ventana 3 → 5 ns → ░░░░░░░░░░░░░░░░░░░░░░░░████████████░░░░░░
    Ventana 4 → 5 ns → ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░████████
    
    Tiempo total: 20 ns de tiempo de reloj

CON paralelización (batch_size=4):
    Ventana 1 → 5 ns → ████████████
    Ventana 2 → 5 ns → ████████████
    Ventana 3 → 5 ns → ████████████
    Ventana 4 → 5 ns → ████████████
    
    Tiempo total: 5 ns de tiempo de reloj (¡4x más rápido!)
```

---

### 8. El Orquestador: `umbrella_suite/`

#### 8.1 Visión General de los Módulos

El paquete `umbrella_suite` es una colección de "especialistas":

```
umbrella_suite/
├── __init__.py         →  Define qué exportar al exterior
├── config.py           →  📋 Formularios de configuración
├── analysis.py         →  🧮 Matemáticas (WHAM, PMF, synthetic data)
├── visualization.py    →  📊 Gráficas y diagnósticos
├── pipeline.py         →  🎭 Director de orquesta
└── run_wnk_pipeline.py →  🖥️ Interfaz de línea de comandos
```

---

#### 8.2 `config.py`: Los Formularios

**Objetivo**: Centralizar TODOS los parámetros en un solo lugar para evitar errores.

##### `UmbrellaWindowConfig`
```python
@dataclass
class UmbrellaWindowConfig:
    center: float              # Centro de la ventana (Å)
    force_constant: float      # k del resorte (kcal/mol/Å²)
    simulation_time_ps: int    # Duración (picosegundos)
```

**Uso conceptual**:
```
Crear configuración para ventana centrada en 10.0 Å:
    Config = UmbrellaWindowConfig
        center = 10.0
        force_constant = 12.0
        simulation_time_ps = 5000
```

##### `UmbrellaPipelineConfig`
```python
@dataclass
class UmbrellaPipelineConfig:
    window_centers: List[float]       # [8.0, 9.2, 10.4, ...]
    force_constant: float             # 12.0
    simulation_time_ps: int           # 5000
    protein_selection: List[str]      # ["A:245:CA", "A:292:CA"]
    pdb_path: Optional[str]           # "wnk_structure.pdb"
    temperature: float                # 300.0 K
    output_dir: str                   # "umbrella_results/"
    num_windows: int                  # Calculado automáticamente
```

**Validación automática**:
```
Si creas un config con window_centers vacío:
    → Python levanta error: "¡Necesitas al menos 2 ventanas!"
    
Si pones force_constant = -5:
    → Error: "¡La constante de fuerza debe ser positiva!"
```

**Métodos helper** (funciones auxiliares):
```python
config.get_window_config(i) 
    → Devuelve UmbrellaWindowConfig para la ventana i
    
config.estimate_total_time_hours()
    → Calcula: (num_windows × simulation_time_ps) / velocidad_simulación
    → Ejemplo: "Tiempo estimado: 6.5 horas"
```

---

#### 8.3 `analysis.py`: El Matemático

**Tres funciones principales**:

##### 1. `compute_pmf()`
```python
def compute_pmf(
    windows: List[UmbrellaWindow],
    method: str = 'WHAM',  # o 'MBAR'
    temperature: float = 300.0
) -> np.ndarray
```

**Input**: Lista de objetos `UmbrellaWindow`, cada uno con:
```python
UmbrellaWindow:
    center: 10.0 Å
    cv_values: [9.8, 9.9, 10.1, 10.0, 10.2, ...]  # 5000 puntos
    time_ps: [0, 1, 2, 3, ..., 4999]
    force_constant: 12.0
```

**Output**: Array del PMF:
```python
[
    (8.0, 15.2),    # (ξ, PMF en kcal/mol)
    (8.5, 12.1),
    (9.0, 8.3),
    (9.5, 4.1),
    (10.0, 0.0),    # ← Mínimo (referencia)
    (10.5, 3.8),
    ...
]
```

**Proceso interno**:
```
1. CONSTRUIR HISTOGRAMAS:
    Para cada ventana:
        ├─ Divide el rango de CV en bins (ej. 0.1 Å)
        ├─ Cuenta cuántas veces el sistema visitó cada bin
        └─ Guarda: histogram[ventana_i][bin_j] = count
    
2. PREPARAR INPUTS PARA PYMBAR:
    ├─ u_kn: Matriz de energías del umbrella
    │        (filas = ventanas, columnas = snapshots)
    ├─ N_k: Número de snapshots por ventana
    └─ Parámetros de bias: [k₁, r₀₁, k₂, r₀₂, ...]
    
3. LLAMAR A PYMBAR.WHAM:
    pymbar_instance = WHAM(u_kn, N_k, ...)
    pmf, uncertainties = pymbar_instance.compute_free_energy()
    
4. POST-PROCESAMIENTO:
    ├─ Restar el mínimo (para que PMF(min) = 0)
    ├─ Suavizar si hay ruido excesivo (opcional)
    └─ Retornar array de PMF
```

**Ejemplo de output visual**:
```
PMF Computation Results:
    
    Minimum at ξ = 10.2 Å (PMF = 0.0 kcal/mol)
    Barrier at ξ = 13.5 Å (PMF = 25.3 kcal/mol)
    
    ΔG(8.0 → 10.2) = -15.2 kcal/mol (favorable)
    ΔG(10.2 → 14.0) = +25.3 kcal/mol (barrera)
```

---

##### 2. `generate_synthetic_windows()`
```python
def generate_synthetic_windows(
    window_centers: List[float],
    force_constant: float,
    n_samples: int = 5000,
    pmf_shape: str = 'double_well'  # o 'barrier', 'harmonic'
) -> List[UmbrellaWindow]
```

**¿Por qué existe?**
- Testing: Validar que WHAM funciona sin correr MD
- Demos: Mostrar flujo completo sin esperar horas
- Debugging: Aislar problemas de análisis vs. simulación

**PMF shapes disponibles**:
```
1. 'double_well' (dos valles):
    PMF
     ↑
  20 │        ╱‾╲
     │       ╱   ╲
  10 │  ╱‾╲ ╱     ╲
     │ ╱   V       ╲___
   0 │╱             Valle₂
     └──────────────────→ ξ
     Valle₁  Barrera
     
2. 'barrier' (un valle + barrera):
    PMF
     ↑
  30 │           ╱╲
     │          ╱  ╲
  15 │         ╱    ╲___
     │    ____╱
   0 │───╱
     └──────────────────→ ξ
     Valle    Barrera
     
3. 'harmonic' (cuenco parabólico):
    PMF
     ↑
  20 │            ╱╲
     │          ╱    ╲
  10 │        ╱        ╲
     │      ╱            ╲
   0 │────╱──────────────╲──
     └──────────────────────→ ξ
```

**Algoritmo de generación**:
```
Para cada ventana i con centro ξ₀ᵢ:
    
    1. Calcular PMF verdadero en ξ₀ᵢ:
        PMF_true(ξ₀ᵢ) = función_shape(ξ₀ᵢ)
        Ejemplo: double_well → evaluar fórmula matemática
    
    2. Crear distribución Gaussiana biased:
        P(ξ | ventana i) ∝ exp[-β(PMF_true(ξ) + V_umbrella(ξ))]
        
        Donde V_umbrella = 0.5 * k * (ξ - ξ₀ᵢ)²
    
    3. Samplear n_samples puntos desde esta distribución:
        cv_values = random.gaussian(mean=ξ₀ᵢ, std=σ)
        
        σ calculada desde la curvatura del potencial total
    
    4. Crear timestamps uniformes:
        time_ps = [0, 1, 2, ..., n_samples-1]
    
    5. Empaquetar en UmbrellaWindow:
        return UmbrellaWindow(
            center=ξ₀ᵢ,
            cv_values=cv_values,
            time_ps=time_ps,
            force_constant=k
        )
```

**Validación del método**:
```
Test: ¿WHAM recupera el PMF que usamos para generar los datos?

    PMF original (conocido): f(ξ)
           ↓
    Generar datos sintéticos con f(ξ)
           ↓
    Correr WHAM sobre esos datos
           ↓
    PMF recuperado: f'(ξ)
           ↓
    Calcular error: RMSD[f(ξ) - f'(ξ)]
    
    Si RMSD < 0.5 kcal/mol → ✅ WHAM funciona correctamente
```

---

##### 3. `UmbrellaWindow` Dataclass
```python
@dataclass
class UmbrellaWindow:
    center: float               # Centro de la ventana
    cv_values: np.ndarray       # Serie temporal de CV
    time_ps: np.ndarray         # Timestamps
    force_constant: float       # k del umbrella
    trajectory_path: Optional[str] = None  # Ruta al .dcd (si existe)
```

**Métodos útiles**:
```python
window.get_histogram(bins=50)
    → Devuelve histogram, bin_edges
    → Útil para visualizar distribución
    
window.mean_cv()
    → Promedio de cv_values
    → Debería estar cerca de 'center' si hay buen sampleo
    
window.std_cv()
    → Desviación estándar
    → Mide qué tan restringido está el sampleo
    → std pequeño → ventana muy rígida
    → std grande → mucho solapamiento con vecinos
```

---

#### 8.4 `visualization.py`: El Artista

**Función principal**: `plot_umbrella_diagnostics()`

**Objetivo**: Generar un reporte visual completo en una sola imagen.

**Layout de 3 paneles**:
```
┌──────────────────────────────────────────────────────────────┐
│  PANEL A: Histogramas Superpuestos                          │
│                                                              │
│  Frecuencia                                                  │
│      ↑                                                       │
│      │  ╱╲   ╱╲   ╱╲   ╱╲   ╱╲   ╱╲                         │
│      │ ╱  ╲ ╱  ╲ ╱  ╲ ╱  ╲ ╱  ╲ ╱  ╲                        │
│      └────────────────────────────────→ CV (Å)              │
│         8   9  10  11  12  13  14                           │
│                                                              │
│  Colores: Una curva por ventana                             │
│  Objetivo: Ver solapamiento                                 │
├──────────────────────────────────────────────────────────────┤
│  PANEL B: Potencial de Fuerza Media (PMF)                   │
│                                                              │
│  ΔG (kcal/mol)                                               │
│      ↑                                                       │
│   30 │              ╱╲                                       │
│   20 │          ___╱  ╲___                                   │
│   10 │     ____╱          ╲____                              │
│    0 │────╱                    ╲──                           │
│      └────────────────────────────────→ CV (Å)              │
│         8   9  10  11  12  13  14                           │
│                                                              │
│  + Barras de error (si están disponibles)                   │
│  + Anotaciones de barreras y mínimos                        │
├──────────────────────────────────────────────────────────────┤
│  PANEL C: Series de Tiempo por Ventana                      │
│                                                              │
│  CV (Å)                                                      │
│   14 │ ········ Window 6 ················                    │
│   13 │ ········ Window 5 ················                    │
│   12 │ ········ Window 4 ················                    │
│   11 │ ········ Window 3 ················                    │
│   10 │ ········ Window 2 ················                    │
│    9 │ ········ Window 1 ················                    │
│    8 │ ········ Window 0 ················                    │
│      └────────────────────────────────→ Time (ps)           │
│        0        1000      2000     3000                      │
│                                                              │
│  Muestra fluctuaciones y estabilidad                        │
└──────────────────────────────────────────────────────────────┘
```

**Interpretación diagnóstica**:

| Observación | Interpretación | Acción |
|---|---|---|
| Gaps en Panel A | Ventanas no se solapan | ✅ Añadir más ventanas intermedias |
| PMF ruidoso (Panel B) | Muestreo insuficiente | ✅ Simular más tiempo |
| Series tiempo con drift (Panel C) | No equilibrado | ✅ Añadir fase de equilibración |
| Barras de error > 2 kcal/mol | Poca estadística | ✅ Más tiempo o más réplicas |

---

#### 8.5 `pipeline.py`: El Director de Orquesta

**Clase `UmbrellaSamplingPipeline`**: Coordina todos los componentes.

**Método principal: `run()`**
```python
def run(self, force_synthetic: bool = False) -> Dict:
```

**Parámetros**:
- `force_synthetic=False`: Correr simulaciones reales con OpenMM
- `force_synthetic=True`: Usar datos sintéticos (para testing)

**Flujo completo**:
```
┌─────────────────────────────────────┐
│ 1. VALIDAR CONFIGURACIÓN            │
│    ✓ ¿Existe pdb_path?              │
│    ✓ ¿protein_selection válido?    │
│    ✓ ¿Parámetros físicos razonables?│
└─────────────────────────────────────┘
          ↓
┌─────────────────────────────────────┐
│ 2. PARSEAR PROTEIN_SELECTION        │
│    Input: ["A:245:CA", "A:292:CA"]  │
│    Output: (atom_idx1, atom_idx2)   │
│              ↓                       │
│    Usa topology del PDB             │
└─────────────────────────────────────┘
          ↓
┌─────────────────────────────────────┐
│ 3. DECIDIR MODO                     │
│    if force_synthetic:              │
│        → generate_synthetic_windows()│
│    else:                            │
│        → _run_async()               │
└─────────────────────────────────────┘
          ↓
┌─────────────────────────────────────┐
│ 4a. MODO REAL (_run_async)          │
│     ├─ Importar calculator          │
│     ├─ Llamar run_full_umbrella_    │
│     │  sampling()                   │
│     └─ Convertir resultados →       │
│        UmbrellaWindow objects       │
└─────────────────────────────────────┘
          OR
┌─────────────────────────────────────┐
│ 4b. MODO SYNTHETIC                  │
│     ├─ Llamar generate_synthetic_   │
│     │  windows()                    │
│     └─ Ya tienes UmbrellaWindow     │
│        objects directamente         │
└─────────────────────────────────────┘
          ↓
┌─────────────────────────────────────┐
│ 5. EXPORTAR DATOS                   │
│    Para cada ventana:               │
│    ├─ Crear output_dir/window_i.dat │
│    └─ Escribir cv_values, time_ps   │
└─────────────────────────────────────┘
          ↓
┌─────────────────────────────────────┐
│ 6. CALCULAR PMF                     │
│    pmf = compute_pmf(windows)       │
└─────────────────────────────────────┘
          ↓
┌─────────────────────────────────────┐
│ 7. VISUALIZAR                       │
│    plot_umbrella_diagnostics(...)   │
│    Guardar umbrella_diagnostics.png │
└─────────────────────────────────────┘
          ↓
┌─────────────────────────────────────┐
│ 8. RETORNAR RESULTADOS              │
│    {                                │
│      'windows': [...],              │
│      'pmf': [...],                  │
│      'output_dir': "...",           │
│      'diagnostics_plot': "..."      │
│    }                                │
└─────────────────────────────────────┘
```

**Manejo de errores**:
```python
# Ejemplo de protección contra errores comunes:

try:
    atom_indices = self._resolve_cv_atoms()
except ValueError as e:
    print(f"❌ Error: {e}")
    print("💡 Tip: Verifica que protein_selection tenga formato:")
    print("   [\"cadena:residuo:átomo\", \"cadena:residuo:átomo\"]")
    sys.exit(1)

try:
    windows = self._run_async()
except ImportError:
    print("⚠️  OpenMM no disponible. Usando modo sintético...")
    windows = generate_synthetic_windows(...)
```

---

#### 8.6 `run_wnk_pipeline.py`: La Interfaz de Usuario

**Propósito**: CLI (Command-Line Interface) para ejecutar pipelines sin escribir código Python.

**Uso básico**:
```bash
# Sintético (demo rápido):
python -m Chronosfold.umbrella_suite.run_wnk_pipeline --synthetic

# Real (simulación completa):
python -m Chronosfold.umbrella_suite.run_wnk_pipeline \
    --pdb-path mi_proteina.pdb \
    --windows 30 \
    --time-ps 50000 \
    --output-dir resultados_wnk/
```

**Argumentos disponibles**:
```
--synthetic / --no-synthetic
    Control de modo de ejecución
    
--windows INTEGER
    Número de ventanas (default: 6)
    
--time-ps INTEGER
    Tiempo de simulación por ventana en picosegundos (default: 5)
    
--pdb-path TEXT
    Ruta al archivo PDB de entrada
    
--output-dir TEXT
    Directorio para guardar resultados (default: "umbrella_results/wnk_pilot")
    
--cv-selection TEXT TEXT
    Selección de átomos para CV (default: "A:245:CA" "A:292:CA")
```

**Ejemplo de output**:
```
🧬 WNK Umbrella Sampling Pipeline
==================================

📋 Configuration:
   - Mode: SYNTHETIC (demo with Gaussian data)
   - Windows: 6 (8.0 to 14.0 Å)
   - Force constant: 12.0 kcal/mol/Å²
   - Simulation time: 5 ps per window
   - Output: umbrella_results/wnk_pilot/

🚀 Running pipeline...

✅ Generated 6 synthetic windows
✅ Computed PMF (range: 0.0 - 27.2 kcal/mol)
✅ Created diagnostics plot

📁 Results saved to:
   - Data: umbrella_results/wnk_pilot/window_*.dat
   - Plot: umbrella_results/wnk_pilot/umbrella_diagnostics.png
   - PMF:  umbrella_results/wnk_pilot/pmf.dat

🎉 Pipeline completed successfully!
```

---

### 9. Visualización y Análisis

#### 9.1 Notebook Interactivo: `umbrella_wham_visualization.ipynb`

**Propósito**: Tutorial paso a paso que combina teoría, código, y análisis.

**Estructura**:
```
Sección 1: Introducción Teórica
    ├─ Conceptos de energía libre
    ├─ Motivación para umbrella sampling
    └─ Ecuaciones WHAM simplificadas

Sección 2: Generar Datos Sintéticos
    ├─ Función para crear PMF de prueba
    ├─ Simulación de histogramas con bias
    └─ Visualización de ventanas individuales

Sección 3: Análisis con WHAM
    ├─ Implementación de WHAM (versión simplificada)
    ├─ Convergencia del algoritmo iterativo
    └─ Comparación PMF recuperado vs. verdadero

Sección 4: Caso de Estudio: WNK
    ├─ Cargar datos reales (si están disponibles)
    ├─ Análisis de convergencia temporal
    ├─ Cálculo de barreras y constantes de tasa
    └─ Interpretación biológica
```

**Ejemplo de celda ejecutable**:
```markdown
### Celda 1: Importar librerías

```python
import numpy as np
import matplotlib.pyplot as plt
from Chronosfold.umbrella_suite import (
    generate_synthetic_windows,
    compute_pmf,
    plot_umbrella_diagnostics
)
```

### Celda 2: Crear datos de prueba

```python
# Configuración
window_centers = np.linspace(8.0, 14.0, 6)
k = 12.0  # kcal/mol/Å²
n_samples = 5000

# Generar
windows = generate_synthetic_windows(
    window_centers=window_centers.tolist(),
    force_constant=k,
    n_samples=n_samples,
    pmf_shape='double_well'
)

print(f"✅ Generadas {len(windows)} ventanas")
```

### Celda 3: Visualizar histogramas

```python
fig, ax = plt.subplots(figsize=(10, 6))

for window in windows:
    hist, bins = window.get_histogram(bins=50)
    ax.plot(bins[:-1], hist, alpha=0.7, label=f'{window.center:.1f} Å')

ax.set_xlabel('CV (Å)')
ax.set_ylabel('Frecuencia')
ax.legend()
plt.show()
```

### Celda 4: Calcular PMF

```python
pmf = compute_pmf(windows, method='WHAM', temperature=300.0)

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(pmf[:, 0], pmf[:, 1], 'o-', linewidth=2)
ax.set_xlabel('CV (Å)')
ax.set_ylabel('PMF (kcal/mol)')
ax.set_title('Potencial de Fuerza Media')
ax.grid(True, alpha=0.3)
plt.show()
```
```

---

#### 9.2 Herramientas CLI Adicionales

##### `wham.py` (Umbrella-visualization/wham.py)

**Uso**:
```bash
python Umbrella-visualization/wham.py \
    --input-dir umbrella_results/wnk_pilot/ \
    --output pmf_final.png \
    --method MBAR
```

**Funcionalidad**:
- Carga archivos `window_*.dat`
- Ejecuta WHAM o MBAR
- Genera plot del PMF con barras de error
- Exporta PMF a archivo de texto

**Output típico**:
```
📊 WHAM Analysis Tool
=====================

📁 Loading windows from: umbrella_results/wnk_pilot/
   Found 6 windows

🧮 Running MBAR analysis...
   Iteration 100/500 (tolerance: 1.2e-05)
   ✅ Converged!

📈 PMF Statistics:
   - Minimum: 0.0 kcal/mol at ξ = 10.2 Å
   - Maximum: 27.3 kcal/mol at ξ = 13.8 Å
   - Barrier height: 27.3 kcal/mol
   - Mean uncertainty: ± 0.8 kcal/mol

💾 Saved:
   - Plot: pmf_final.png
   - Data: pmf_final.dat
```

---

### 10. Scripts de Lanzamiento Multiplataforma

#### 10.1 Windows: `run_umbrella.ps1` y `bootstrap_windows.ps1`

**Problema que resuelven**: Conda no siempre está en el PATH de Windows.

**Solución**: Auto-detección inteligente.

**Lugares donde busca Conda**:
```
1. Variable de entorno $env:CONDA_EXE
   Ejemplo: C:\Users\Usuario\miniconda3\Scripts\conda.exe

2. Comando 'conda' en PATH
   Verifica con: Get-Command conda

3. Rutas comunes de instalación:
   - C:\ProgramData\miniconda3\Scripts\conda.exe
   - C:\ProgramData\Anaconda3\Scripts\conda.exe
   - $env:USERPROFILE\miniconda3\Scripts\conda.exe
   - $env:USERPROFILE\Anaconda3\Scripts\conda.exe
   - $env:LOCALAPPDATA\Continuum\miniconda3\Scripts\conda.exe

4. Fallback: Preguntar al usuario
```

**Flujo de ejecución**:
```
┌─────────────────────────────────────┐
│ bootstrap_windows.ps1               │
│ (Script de entrada)                 │
└─────────────────────────────────────┘
          ↓
┌─────────────────────────────────────┐
│ ¿Usuario pasó -CondaPath?           │
│ Sí → Usar esa ruta                  │
│ No → Llamar Get-CondaCommand        │
└─────────────────────────────────────┘
          ↓
┌─────────────────────────────────────┐
│ Get-CondaCommand                    │
│ Intenta 4 métodos de detección     │
│ Devuelve ruta o $null              │
└─────────────────────────────────────┘
          ↓
┌─────────────────────────────────────┐
│ ¿Se encontró Conda?                 │
│ Sí → Continuar                      │
│ No → Mostrar error con instrucciones│
└─────────────────────────────────────┘
          ↓
┌─────────────────────────────────────┐
│ Ensure-CondaEnv "bsm-lancad-env"    │
│ ¿Existe el environment?             │
│ Sí → Activar                        │
│ No → Crear desde environment.yml    │
└─────────────────────────────────────┘
          ↓
┌─────────────────────────────────────┐
│ Ejecutar:                           │
│ python -m Chronosfold.umbrella_     │
│ suite.run_wnk_pipeline --synthetic  │
└─────────────────────────────────────┘
```

**Ejemplo de uso**:
```powershell
# Detección automática:
.\scripts\bootstrap_windows.ps1

# Especificar ruta manualmente:
.\scripts\bootstrap_windows.ps1 -CondaPath "C:\miniconda3\Scripts\conda.exe"
```

---

#### 10.2 Linux/Mac: `run_umbrella.sh` y `bootstrap_linux.sh`

**Diferencias con Windows**:
- Usa `bash` en lugar de PowerShell
- Sintaxis POSIX-compliant
- Busca en `$HOME/miniconda3/bin/conda`

**Flujo idéntico**, solo cambia la sintaxis:
```bash
#!/bin/bash

# Detectar Conda
if [ -n "$CONDA_BIN" ]; then
    CONDA_CMD="$CONDA_BIN"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
elif [ -f "$HOME/miniconda3/bin/conda" ]; then
    CONDA_CMD="$HOME/miniconda3/bin/conda"
else
    echo "❌ Conda not found!"
    exit 1
fi

# Activar environment
eval "$($CONDA_CMD shell.bash hook)"
conda activate bsm-lancad-env || {
    echo "Creating environment..."
    conda env create -f environment.yml
    conda activate bsm-lancad-env
}

# Ejecutar pipeline
python -m Chronosfold.umbrella_suite.run_wnk_pipeline --synthetic
```

**Ejemplo de uso**:
```bash
# Detección automática:
./scripts/bootstrap_linux.sh

# Especificar Conda manualmente:
export CONDA_BIN=/opt/conda/bin/conda
./scripts/bootstrap_linux.sh
```

---

**FIN DE LA PARTE II: IMPLEMENTACIÓN TÉCNICA**

---

## PARTE III: APLICACIÓN A WNK

### 11. La Quinasa WNK: Un Sensor Molecular Extraordinario

#### 11.1 Contexto Biológico: ¿Por Qué Importa WNK?

**WNK** (With-No-Lysine [K]) es una familia de serina/treonina quinasas con propiedades únicas:

**1. Homeostasis de Sal y Agua**
```
Rol fisiológico:
    WNK1/4 → Fosforila OSR1/SPAK → Regula NCC/NKCC
                                    (cotransportadores de iones)
    
Efecto neto:
    - Controla reabsorción de Na⁺/Cl⁻ en riñón
    - Mantiene presión arterial
    - Responde a volumen extracelular
```

**Evidencia clínica**:
- Mutaciones en WNK1/4 → Síndrome de Gordon (hipertensión + hipercalemia)
- Mutaciones en KLHL3/CUL3 (reguladores de WNK) → Mismo fenotipo
- Diuréticos tiazidas (fármacos) actúan downstream de WNK

**Referencias**:

**[6] Jonniya NA, Sk MF, Kar P. (2019)**  
*"Investigating phosphorylation-induced conformational changes in WNK1 kinase by molecular dynamics simulations"*  
ACS Omega, 4(17):17404-17416.  
doi: 10.1021/acsomega.9b02368  
**Citaciones: 46**

Este estudio usando MD de 200 ns mostró que la fosforilación en el loop de activación induce compactación del sitio activo de WNK1, relevante para nuestro modelo de activación.

**[7] Zhang J, Siew K, Macartney T, et al. (2015)**  
*"Critical role of the SPAK protein scaffold in regulating

 blood pressure in response to kidney injury"*  
Journal of the American Society of Nephrology, 26(10):2367-2380.  
doi: 10.1681/ASN.2014070672  
**Citaciones: 35**

Demostró que el dominio CCT de SPAK es esencial para la señalización WNK-SPAK-NCC en regulación de presión arterial.

---

**2. Sensor de Estrés Osmótico**
```
Mecanismo propuesto (Boyd-Shiwarski et al., Zhang et al.):

    Condiciones isotónicas:
        WNK1 forma dímeros auto-inhibidos
        ↓
        Inactivo (no fosforila OSR1/SPAK)
    
    Estrés hiperosmótico (↑ osmolitos):
        Dímero se disocia → Monómeros activos
        ↓
        Fosforila OSR1/SPAK → Activa NCC
        ↓
        Reabsorción de Na⁺ → Retención de agua
```

**Agentes osmóticos naturales**:
- **PEG** (polietilenglicol): Excluído del volumen, favorece monómero
- **Sacarosa**: Efecto similar
- **NaCl**: A concentraciones fisiológicas altas (>150 mM)

**Hipótesis termodinámica**:
> La transición dímero → monómero está impulsada por la liberación de moléculas de agua ordenadas en la interfaz dimérica. Los osmolitos aumentan la "recompensa" entrópica de liberar esas aguas.

---

#### 11.2 Arquitectura de Dominios de WNK

**Estructura modular**:
```
WNK1 (2382 aminoácidos):

H₂N─┬─────────┬──────────┬────────────────────┬───────────────┬─COOH
    │         │          │                    │               │
    KD        AUTOINH    COILED-COIL         CCT (C-term     PxxP
    (1-491)   (492-555)  (556-828)           tail, 1802-     motifs
                                              2097)
    
    KD: Kinase Domain (dominio catalítico)
        - Lóbulo N: β-sheets + αC-helix
        - Lóbulo C: α-helices
        - Sitio activo entre lóbulos
        - Loop de activación (T-loop): residuos 420-450
    
    AUTOINH: Región autoinhibitoria
        - Interactúa con KD para suprimir actividad basal
    
    COILED-COIL: Mediador de dimerización
        - Heptad repeats clásicos (abcdefg)
        - Residuos hidrofóbicos en posiciones 'a' y 'd'
    
    CCT (Conserved C-Terminal): Dominio de interacción
        - Reconoce motivo RFxV en sustratos
        - ESTRUCTURA: 4-5 α-helices antiparalelas
        - ~300 aminoácidos
```

**Diagrama esquemático**:
```
Vista lateral (estado dimérico):

    Monómero A                    Monómero B
    
    ╔═══════╗                    ╔═══════╗
    ║  KD   ║                    ║  KD   ║
    ║ (inactivo)                 ║ (inactivo)
    ╚═══╤═══╝                    ╚═══╤═══╝
        │                            │
        ├─── Autoinh ────────────────┤
        │    (cross-inhibición)      │
        │                            │
    ════╧════════════════════════════╧════
         Coiled-coil (interfaz)
```

---

#### 11.3 El Motivo RFxV y el Dominio CCT

##### El Código de Reconocimiento Molecular

**Motivo consenso**: **R-F-x-V/I**
- **R** (Arg): Esencial, carga positiva (+)
- **F** (Phe): Hidrofóbico aromático, ancla en bolsillo
- **x**: Cualquier aminoácido (típicamente S, T, A)
- **V/I** (Val/Ile): Hidrofóbico alifático

**Sustratos de WNK que contienen RFxV**:
```
Proteína          Secuencia RFxV      Posición
────────────────────────────────────────────────
OSR1              R-F-T-V             Residuos 430-433
SPAK              R-F-A-V             Residuos 428-431
WNK1 (sí mismo)   R-F-Q-V             Residuos 491-494 (autoinhibición)
NRBP1 (pseudokinasa) R-F-X-V          Múltiples sitios
```

**Estructura del dominio CCT** (basada en PDB 2LRU, SPAK):
```
Vista superior (mirando hacia el sitio de unión):

           Hélice α3
              │
       ╱──────┴──────╲
      ╱   Bolsillo    ╲
     │    de unión     │
     │   para RFxV     │  ← Surco hidrofóbico
     │                 │
  Hélice α2 ──────────── Hélice α4
     │                 │
      ╲               ╱
       ╲─────────────╱
         Hélice α5

Residuos clave del CCT de SPAK:
    - Arg465: Forma puente salino con Arg del motivo RFxV
    - Phe468, Leu472: Definen bolsillo hidrofóbico para Phe
    - Val493, Ile497: Acomodan Val/Ile del motivo
```

**Mecanismo de reconocimiento**:
```
PASO 1: Acercamiento inicial
    Arg del RFxV (carga +) es atraída electrostáticamente
    hacia residuos negativos en la superficie del CCT
    
PASO 2: Anclaje de Phe
    Phe del motivo se inserta en bolsillo hidrofóbico
    → Mayor superficie de contacto
    → Incrementa afinidad
    
PASO 3: Cierre conformacional
    Val/Ile final estabiliza el complejo
    → Interfaz cerrada y específica
```

**Afinidades de unión** (estimadas por SPR, Surface Plasmon Resonance):
```
WNK1-CCT + OSR1-RFxV:  Kd ~ 1-5 μM
WNK1-CCT + SPAK-RFxV:  Kd ~ 1-5 μM

Para comparación:
    - Interacción típica proteína-proteína transitoria: Kd ~ 0.1-10 μM
    - Complejo enzima-sustrato: Kd ~ 10-100 μM
    
→ WNK-CCT tiene afinidad ALTA por RFxV
```

**Referencias sobre RFxV y CCT**:

**[8] Taylor SS, Meharena HS, Kornev AP. (2024)**  
*"WNK kinases and the OSR1/SPAK kinases: Unraveling the structural basis for regulation"*  
Structure (London), 32(11):1975-1978.  
doi: 10.1016/j.str.2024.09.007  

**[9] Taylor SS, Meharena HS, Kornev AP. (2025)**  
*"WNK kinases and the OSR1/SPAK kinases: Unraveling the structural basis for regulation, Part II"*  
Structure (London), 33(1):1-4.  
doi: 10.1016/j.str.2024.11.006  

Estos trabajos recientes de Susan Taylor (experta mundial en quinasas) describen la estructuración del dominio CCT y su papel en la transmisión alostérica de señales osmóticas a través de WNK.

---

### 12. Nuestro Objetivo Científico: Mapear el Paisaje Termodinámico

#### 12.1 Hipótesis Central

**Pregunta**: ¿Cuál es el perfil de energía libre para la apertura/cierre de los lóbulos N y C del dominio quinasa de WNK1?

**Hipótesis**:
```
H1: El estado inactivo (dímero) corresponde a una configuración
    COMPACTA del dominio quinasa (distancia β3-αC pequeña).
    
H2: El estado activo (monómero) requiere APERTURA de los lóbulos
    (distancia β3-αC mayor) para permitir acceso de ATP y sustrato.
    
H3: Existe una BARRERA energética (~15-25 kcal/mol) entre ambos
    estados, consistente con la observación de que la activación
    requiere estrés osmótico.
    
H4: Los osmolitos REDUCEN la barrera al estabilizar el estado
    abierto (monómero) mediante efectos de volumen excluido.
```

**Predicción testeable**:
```
Si H1-H4 son correctas:
    
    PMF(ξ) donde ξ = distancia(β3-CA, αC-CA) mostrará:
    
    PMF
     ↑
  25 │              ╱╲          ← Barrera conformacional
     │         ____╱  ╲____
  15 │        ╱            ╲
     │   ____╱              ╲___
   0 │──╱                       ╲─
     └───────────────────────────→ ξ (Å)
     8.0  (dímero)   12.0   14.0 (monómero)
     Estado A              Estado B
```

---

#### 12.2 Diseño Experimental de la Simulación

**Coordenada Colectiva (CV)**:
```python
CV = distance(CA_residue_245, CA_residue_292)
```

**Justificación estructural**:

**Residuo 245** (en β3-strand):
- Parte del lóbulo N-terminal
- Elemento estructural estable (β-sheet)
- Adyacente al sitio de unión de ATP
- Conservado en familia de quinasas

**Residuo 292** (en αC-helix):
- Hélice regulatoria crítica
- Su posición determina estado activo/inactivo
- En quinasas activas: αC está "IN" (cerca de sitio activo)
- En quinasas inactivas: αC está "OUT" (alejada)

**Cambio esperado**:
```
Estado Inactivo (dímero):
    αC comprimida contra β3
    → distancia ~8.0-9.0 Å
    → ATP no puede posicionarse correctamente
    
Estado Activo (monómero):
    αC rota hacia afuera
    → distancia ~13.0-14.0 Å
    → Loop de activación accesible
    → ATP alineado con sustrato
```

**Diagrama de elementos estructurales**:
```
Vista top del sitio activo de WNK1:

                  Lóbulo C
                    ║
        ATP         ║
         ╱╲         ║
        ╱  ╲        ║
    ───●────●───────╫─────  ← Loop de activación (T420)
     β3│    │αC     ║
       │245 │292    ║
       │ ↕  │       ║      ↕ = Distancia CV
       │ d  │       ║
    ───┴────┴───────╫─────  ← Sitio activo
       │            ║
     Lóbulo N       ║

Cuando d aumenta:
    - αC se aleja de sitio activo
    - Loop de activación se libera
    - Sustrato puede acceder
```

---

#### 12.3 Parámetros de la Simulación

**Configuración del sistema**:
```yaml
# Estructura inicial
pdb_file: "wnk1_dimer_6cn9.pdb"
    Fuente: Protein Data Bank, estructura cristalográfica
    Resolución: 2.8 Å
    Estado: Complejo WNK1-WNK1 (dímero inactivo)

# Campo de fuerzas
forcefield: "amber14-all.xml"
    Justificación: AMBER14 optimizado para proteínas
    Alternativa: CHARMM36m (similar precisión)

# Solvente
water_model: "tip3p"
box_padding: 12 Å
    → Asegura que proteína no interactúa con sus imágenes periódicas

ions:
    Na+: Para neutralizar carga
    Cl-: Concentración fisiológica (150 mM)

# Temperatura y presión
temperature: 300 K  (27°C, fisiológico)
pressure: 1 bar
barostat: Monte Carlo (cada 25 pasos)
thermostat: Langevin (fricción = 1/ps)

# Umbrella sampling
num_windows: 30
window_range: [8.0, 14.0] Å
spacing: 0.2 Å (overlap garantizado)
force_constant: 12.0 kcal/mol/Å²

# Tiempos de simulación
equilibration: 5 ns por ventana
production: 50 ns por ventana
    Total: 55 ns × 30 ventanas = 1.65 μs tiempo agregado
```

**Estimación de recursos computacionales**:
```
Hardware: NVIDIA RTX 3090 (24 GB VRAM)
Velocidad estimada: ~150 ns/día para este sistema (~50,000 átomos)

Tiempo de reloj:
    55 ns por ventana / 150 ns/día = 0.37 días
    0.37 días × 30 ventanas = 11 días
    
    Con paralelización (4 GPUs):
        11 días / 4 = ~3 días de tiempo real

Almacenamiento:
    Trayectorias: ~500 MB por ventana × 30 = 15 GB
    Checkpoints y logs: ~5 GB
    Total: ~20 GB
```

---

#### 12.4 Análisis Biológico del PMF

Una vez obtenido el PMF, podemos extraer información cuantitativa:

##### 1. Constantes de Equilibrio
```
PMF(ξ) → ΔG°(Estado A → Estado B)

Ejemplo:
    PMF(8.0 Å) = 0.0 kcal/mol  (dímero, referencia)
    PMF(13.5 Å) = 25.3 kcal/mol (monómero)
    
    ΔG° = 25.3 kcal/mol
    
Constante de equilibrio:
    K_eq = exp(-ΔG°/RT)
    
    Con R = 1.987 cal/(mol·K), T = 300 K:
    K_eq = exp(-25300 / (1.987 × 300))
         = exp(-42.4)
         = 5.8 × 10⁻¹⁹
         
Interpretación:
    En equilibrio sin osmolitos:
    [Monómero] / [Dímero] = 5.8 × 10⁻¹⁹
    
    → Prácticamente TODO WNK está en forma dimérica
    → Consistente con baja actividad basal
```

##### 2. Constantes de Tasa (Transition State Theory)
```
Altura de la barrera: ΔG‡ = 25.3 kcal/mol (desde el mínimo)

Frecuencia de cruce de barrera:
    k = (kT/h) × exp(-ΔG‡/RT)
    
    Donde:
    - k: Constante de Boltzmann (1.381×10⁻²³ J/K)
    - T: Temperatura (300 K)
    - h: Constante de Planck (6.626×10⁻³⁴ J·s)
    - R: Constante de gases (8.314 J/(mol·K))
    
    Prefactor kT/h ~ 6.25 × 10¹² s⁻¹ a 300 K
    
    k = 6.25 × 10¹² × exp(-25300 cal/mol / (1.987 cal/(mol·K) × 300 K))
      = 6.25 × 10¹² × exp(-42.4)
      = 6.25 × 10¹² × 5.8 × 10⁻¹⁹
      = 3.6 × 10⁻⁶ s⁻¹
      
Tiempo medio de la transición:
    τ = 1/k = 2.8 × 10⁵ segundos
        ≈ 77 horas
        ≈ 3.2 días
        
Interpretación:
    Sin osmolitos, WNK necesitaría ~3 días para activarse espontáneamente
    → Biológicamente irrelevante
    → Requiere estímulo osmótico para ser fisiológicamente útil
```

##### 3. Efecto de Osmolitos (Simulación Futura)
```
Si repetimos umbrella sampling en presencia de PEG 20%:

Predicción:
    PMF_PEG(13.5 Å) = 18.0 kcal/mol (en lugar de 25.3)
    
    Reducción de barrera: ΔΔG‡ = 25.3 - 18.0 = 7.3 kcal/mol
    
Nueva constante de tasa:
    k_PEG = 6.25 × 10¹² × exp(-18000 / (1.987 × 300))
          = 6.25 × 10¹² × exp(-30.2)
          = 6.25 × 10¹² × 9.4 × 10⁻¹⁴
          = 5.9 × 10⁻¹ s⁻¹
          ≈ 0.6 eventos/segundo
          
Tiempo medio:
    τ_PEG = 1 / 0.6 s⁻¹ = 1.7 segundos
    
¡Aceleración de 10⁶ veces! (de 3 días a 2 segundos)
    → Fisiológicamente relevante
    → Explica respuesta rápida a estrés osmótico
```

---

### 13. Del Demo Sintético a las Simulaciones Reales

#### 13.1 Roadmap en 3 Fases (recordatorio)

##### **FASE 1: Proof-of-Concept Sintético** ✅ **COMPLETA**
```
Objetivo: Validar infraestructura de software sin costo computacional

Logros:
    ✅ Pipeline funcional end-to-end
    ✅ Generación de datos sintéticos (generate_synthetic_windows)
    ✅ Análisis WHAM/MBAR operativo
    ✅ Visualización de diagnósticos
    ✅ Scripts multiplataforma (Windows/Linux)
    ✅ Documentación completa

Tiempo invertido: 2 días (desarrollo + testing)
Costo computacional: Cero
```

##### **FASE 2: Piloto con OpenMM** 🔄 **EN PROGRESO**
```
Objetivo: Simulaciones cortas (5-10 ns) para validar setup físico

Tareas:
    1. Preparar estructura PDB de WNK1:
       - Descargar 6CN9 desde PDB
       - Aislar monómero A (chain A)
       - Añadir hidrógenos (pdb4amber)
       - Verificar integridad (missing loops, clashes)
       
    2. Correr 6 ventanas de prueba (8-14 Å):
       - 1 ns equilibración
       - 5 ns producción
       - Sin réplicas (single run)
       
    3. Validaciones:
       - ¿Convergencia de la CV en cada ventana?
       - ¿Overlap entre histogramas?
       - ¿PMF libre de artefactos numéricos?
       - ¿RMSD de la proteína estable (<3 Å)?
       
    4. Debugging:
       - Ajustar force_constant si ventanas no se solapan
       - Aumentar equilibración si CV no converge
       - Revisar campo de fuerzas si proteína se desnaturaliza

Tiempo estimado: 1 semana (incluyendo preparación)
Costo computacional: ~50 horas GPU
```

##### **FASE 3: Producción Completa** 📅 **FUTURO**
```
Objetivo: PMF de alta precisión para publicación

Configuración:
    - 30 ventanas (spacing 0.2 Å)
    - 5 ns equilibración por ventana
    - 50 ns producción por ventana
    - 3 réplicas independientes (para errores estadísticos)
    
Análisis avanzado:
    - Convergencia por bloques (block averaging)
    - Errores con bootstrapping (100 iteraciones)
    - Comparación WHAM vs. MBAR
    - Descomposición energética (qué fuerzas dominan la barrera)
    
Extensiones opcionales:
    - Simular con osmolitos (PEG, sacarosa)
    - Simular mutantes (ej. WNK4, mutaciones de enfermedad)
    - Calcular efectos de fosforilación (T-loop fosforilado)

Tiempo estimado: 1 mes
Costo computacional: ~1,500 horas GPU (paralelo en 4 GPUs → 2 semanas)
```

---

#### 13.2 Checklist Técnico para Fase 2

```
□ PREPARACIÓN DE ESTRUCTURA
  □ Descargar PDB 6CN9
  □ Extraer chain A (monómero)
  □ Añadir hidrógenos con pdb4amber o pdbfixer
  □ Verificar protonación de His, Asp, Glu (pH 7.0)
  □ Verificar quiralidad de aminoácidos
  □ Eliminar moléculas de agua cristalográficas (opcional)
  □ Guardar como wnk1_monomer_prepared.pdb

□ CONFIGURACIÓN DE SIMULACIÓN
  □ Crear config.py con:
      window_centers = np.linspace(8.0, 14.0, 6)
      force_constant = 12.0
      simulation_time_ps = 5000
      equilibration_time_ps = 1000
      protein_selection = ["A:245:CA", "A:292:CA"]
  □ Verificar que átomos existen en el PDB (índices válidos)
  
□ TEST RUN (1 ventana)
  □ Correr window centrada en 10.0 Å (mínimo esperado)
  □ Monitorear:
      - Energía total (debe estabilizarse tras equilibración)
      - Temperatura (debe fluctuar ~300 ± 5 K)
      - Presión (fluctuaciones grandes normales, promedio ~1 bar)
      - RMSD proteína (< 3 Å vs. estructura inicial)
      - CV (debe oscilar ~10.0 ± 0.5 Å)
  □ Visualizar trayectoria en VMD o PyMOL
      - Buscar desnaturalización, clashes, comportamiento no físico
      
□ FULL RUN (6 ventanas)
  □ Ejecutar pipeline completo
  □ Verificar que todos los jobs terminaron sin errores
  □ Cargar resultados y generar plot de diagnóstico
  □ Evaluar overlap de histogramas
  
□ ANÁLISIS PMF
  □ Calcular PMF con WHAM
  □ Calcular PMF con MBAR (comparación)
  □ Verificar que:
      - PMF es suave (sin saltos abruptos)
      - Mínimo está cerca de ξ ~ 10 Å (estado cristalográfico)
      - Barrera existe y está en rango razonable (15-30 kcal/mol)
      - Errores < 2 kcal/mol en regiones importantes
      
□ DECISIÓN: ¿PROCEDER A FASE 3?
  □ Si todas las validaciones pasan → Extender a 30 ventanas, 50 ns
  □ Si hay problemas → Iterar debugging en Fase 2
```

---

### 14. Interpretación Biológica y Perspectivas Futuras

#### 14.1 Conectando Termodinámica con Función

**Pregunta clave**: ¿Cómo el PMF obtenido explica el comportamiento fisiológico de WNK?

**Escenario 1: Condiciones Isotónicas (sin estrés)**
```
PMF muestra:
    - Valle profundo en ξ ~ 8-9 Å (dímero compacto)
    - Barrera de ~25 kcal/mol
    - Estado de monómero (ξ ~ 13-14 Å) muy desfavorable
    
Consecuencia biológica:
    → WNK permanece dimérico e inactivo
    → OSR1/SPAK no se fosforilan
    → NCC/NKCC permanecen inactivos
    → Homeostasis de Na⁺/Cl⁻ en estado basal
```

**Escenario 2: Estrés Hiperosmótico (↑ osmolitos)**
```
PMF shift prediction (basado en volumen excluido):
    - Valle del dímero: Sin cambio significativo
    - Barrera: Reducida a ~18 kcal/mol (ΔΔG ≈ -7 kcal/mol)
    - Estado de monómero: Estabilizado (más favorable)
    
Mecanismo molecular:
    Osmolitos (PEG, sacarosa) → Excluded volume effect
                              ↓
                   Favorece estado con menor superficie de solvatación
                              ↓
                   Monómero libera ~200 moléculas de agua
                              ↓
                   ΔS aumenta dramáticamente
                              ↓
                   ΔG(dímero → monómero) disminuye
                   
Consecuencia biológica:
    → WNK se disocia rápidamente (segundos a minutos)
    → Monómeros fosforilan OSR1/SPAK
    → Cascada de señalización activada
    → Reabsorción de Na⁺ → Retención de agua
    → Ajuste homeostático completado
```

---

#### 14.2 Hipótesis Adicionales para Explorar

**1. Rol del Dominio CCT en la Estabilización del Monómero**
```
Pregunta: ¿El dominio CCT (residuos 1802-2097) estabiliza el monómero
          al interactuar intramolecularmente con motivos RFxV?
          
Experimento computacional:
    - Simular WNK1 completo (con CCT) vs. truncado (sin CCT)
    - Comparar PMFs
    - Hipótesis: ΔG_barrera será menor en WNK completo
    
Predicción:
    Si el CCT estabiliza el monómero:
        PMF_completo tendrá barrera más baja que PMF_truncado
        → Validaría modelo de "pinza molecular" (CCT agarra RFxV
           del propio KD para estabilizar conformación activa)
```

**2. Efecto de Mutaciones Causantes de Enfermedad**
```
Mutaciones de interés (Síndrome de Gordon):
    - WNK1-D368A: Mutación en sitio activo (cambio de carga)
    - WNK4-Q562E: Mutación en coiled-coil (afecta dimerización)
    
Experimento computacional:
    - Modelar mutantes por mutagénesis in silico
    - Correr umbrella sampling
    - Comparar PMF_WT vs. PMF_mutante
    
Predicciones:
    a) Si mutación desestabiliza dímero:
        → Barrera más baja en mutante
        → WNK hiperactivo
        → Hipertensión (observado clínicamente ✅)
        
    b) Si mutación estabiliza dímero:
        → Barrera más alta
        → WNK hipoactivo
        → Hipotensión (menos común)
```

**3. Influencia de Fosforilación en el Loop de Activación**
```
Sitio de fosforilación: Thr420 (en el loop de activación)

Estado:
    - No fosforilado: WNK basal
    - Fosforilado: WNK máximamente activo
    
Hipótesis:
    La fosforilación en T420 introduce carga negativa (-2 tras fosforilación)
    → Repulsión electrostática con residuos negativos cercanos
    → Fuerza al loop a adoptar conformación extendida
    → Estabiliza estado de monómero activo
    
Experimento:
    - Modificar topología de T420 para simular fosforilación
      (Thr → pThr con carga -2)
    - Correr umbrella sampling
    - Comparar PMF_pT420 vs. PMF_T420
    
Predicción:
    ΔG_barrera(pT420) < ΔG_barrera(T420) por ~5-10 kcal/mol
```

---

#### 14.3 Comparación con Datos Experimentales

**Técnicas experimentales para validar PMF**:

##### 1. Single-Molecule FRET (smFRET)
```
Concepto:
    - Marca β3-strand y αC-helix con fluoróforos (donor y acceptor)
    - Mide eficiencia FRET (depende de distancia r)
    - Histograma de FRET → Distribución de distancias
    
Predicción desde PMF:
    P(r) ∝ exp[-PMF(r) / kT]
    
    Si PMF tiene dos valles (dímero y monómero):
        → Histograma de FRET mostrará dos poblaciones
        
Comparación cuantitativa:
    Experimental: Fracción de estado abierto vs. [osmolito]
    Simulación: Calcular K_eq(osmolito) desde PMF(osmolito)
```

##### 2. Analytical Ultracentrifugation (AUC)
```
Medida:
    Coeficiente de sedimentación (S) → Masa molecular aparente
    
Interpretación:
    S_dímero ≈ 6.5 S (2 × 55 kDa = 110 kDa)
    S_monómero ≈ 4.2 S (1 × 55 kDa)
    
Experimento:
    Medir S en función de concentración de osmolito
    → Extraer K_eq(dímero ⇌ monómero)
    
Comparación con simulación:
    K_eq_exp vs. K_eq_calc desde ΔG° en PMF
```

##### 3. Differential Scanning Calorimetry (DSC)
```
Medida:
    Estabilidad térmica (T_m = temperatura de melting)
    
Interpretación:
    Estado más compacto → T_m mayor
    Estado más abierto → T_m menor
    
Predicción desde PMF:
    Si PMF muestra que el estado abierto es menos estable:
        → T_m_monómero < T_m_dímero
        
Experimental:
    Medir T_m con y sin osmolitos
    → Shift en T_m correlaciona con ΔΔG del PMF
```

---

**FIN DE LA PARTE III: APLICACIÓN A WNK**

---

## APÉNDICES

### Apéndice A: Glosario de Términos

**Términos Termodinámicos**:

- **Energía Libre de Gibbs (ΔG)**: Cantidad de energía disponible para realizar trabajo útil en un proceso a temperatura y presión constantes. Combina entalpía (ΔH) y entropía (ΔS): ΔG = ΔH - TΔS.

- **Entalpía (ΔH)**: Cambio en el contenido de calor de un sistema. Refleja la formación o ruptura de enlaces.

- **Entropía (ΔS)**: Medida del desorden o número de microestados accesibles. Mayor entropía = más desorden.

- **Potencial de Fuerza Media (PMF)**: Energía libre como función de una coordenada colectiva. Representa el "paisaje energético" del sistema proyectado sobre esa coordenada.

- **Coordenada Colectiva (CV)**: Variable macroscópica que describe la configuración del sistema (ej. distancia, ángulo). En inglés: Collective Variable.

- **Barrera de Energía Libre (ΔG‡)**: Diferencia de energía entre el estado inicial y el estado de transición. Determina la velocidad de la reacción.

- **Constante de Equilibrio (K_eq)**: Relación entre concentraciones de productos y reactivos en equilibrio. K_eq = exp(-ΔG°/RT).

**Términos de Simulación**:

- **Dinámica Molecular (MD)**: Método computacional que simula el movimiento de átomos según las leyes de Newton.

- **Campo de Fuerzas (Forcefield)**: Conjunto de ecuaciones y parámetros que describen las interacciones entre átomos (enlaces, ángulos, cargas, van der Waals).

- **Integrador**: Algoritmo numérico para resolver ecuaciones de movimiento. Común: Verlet, Langevin.

- **Paso de Tiempo (Timestep)**: Incremento temporal entre cálculos sucesivos. Típicamente 1-2 femtosegundos.

- **Equilibración**: Fase inicial de simulación para relajar el sistema y alcanzar condiciones deseadas (T, P).

- **Producción**: Fase de simulación después de equilibración, donde se recolectan datos para análisis.

- **Barostat**: Algoritmo para mantener presión constante (ej. Monte Carlo, Berendsen).

- **Termostato**: Algoritmo para mantener temperatura constante (ej. Langevin, Nosé-Hoover).

**Términos de Umbrella Sampling**:

- **Potencial de Umbrella (Bias Potential)**: Potencial artificial añadido para forzar el sistema a explorar regiones específicas. Típicamente armónico: V = ½k(ξ-ξ₀)².

- **Ventana (Window)**: Simulación individual con un centro de umbrella fijo (ξ₀).

- **Constante de Fuerza (Force Constant, k)**: Rigidez del resorte del umbrella. Unidades: kcal/mol/Å².

- **Solapamiento (Overlap)**: Región donde dos ventanas consecutivas muestrean valores de CV comunes. Crítico para WHAM.

- **WHAM (Weighted Histogram Analysis Method)**: Algoritmo para combinar histogramas de múltiples ventanas y reconstruir el PMF sin sesgo.

- **MBAR (Multistate Bennett Acceptance Ratio)**: Método más general que WHAM, óptimo estadísticamente.

**Términos Biológicos**:

- **Quinasa (Kinase)**: Enzima que transfiere grupos fosfato desde ATP a residuos Ser, Thr, o Tyr en proteínas sustrato.

- **Fosforilación**: Adición de un grupo fosfato (PO₄³⁻) a una proteína. Mecanismo clave de regulación.

- **Dimerización**: Formación de un complejo de dos moléculas (dímero). Puede ser homo- (dos copias iguales) o hetero-.

- **Autoinhibición**: Mecanismo donde una región de la proteína bloquea su propio sitio activo.

- **Alosterismo**: Regulación de la actividad enzimática mediante la unión de una molécula en un sitio distante del sitio activo.

- **Coiled-Coil**: Estructura secundaria donde dos o más α-hélices se enrollan entre sí. Común en mediadores de oligomerización.

- **Motivo RFxV**: Secuencia lineal corta (Arg-Phe-X-Val) que actúa como sitio de reconocimiento molecular.

- **Dominio CCT (Conserved C-Terminal)**: Módulo de ~300 aminoácidos que reconoce motivos RFxV. Presente en WNK, OSR1, SPAK.

- **Estrés Osmótico**: Cambio en la concentración de solutos que afecta el equilibrio de agua a través de membranas.

- **Osmolito**: Molécula pequeña que contribuye a la presión osmótica (ej. PEG, sacarosa, urea).

**Términos Computacionales**:

- **Dataclass**: Clase de Python (desde 3.7) que simplifica la creación de objetos que almacenan datos.

- **Pipeline**: Secuencia automatizada de pasos computacionales.

- **Async/Await**: Paradigma de programación asíncrona en Python para ejecutar tareas en paralelo.

- **Conda**: Gestor de paquetes y entornos para Python y otras herramientas científicas.

- **CLI (Command-Line Interface)**: Interfaz basada en texto para interactuar con programas.

- **Notebook (Jupyter)**: Documento interactivo que mezcla código ejecutable, texto narrativo, y visualizaciones.

---

### Apéndice B: Referencias Bibliográficas

**Artículos Fundacionales en Umbrella Sampling y WHAM**:

**[1] Kumar S, Rosenberg JM, Bouzida D, Swendsen RH, Kollman PA. (1992)**  
"THE weighted histogram analysis method for free-energy calculations on biomolecules. I. The method"  
*Journal of Computational Chemistry*, 13(8):1011-1021.  
DOI: [10.1002/jcc.540130812](https://doi.org/10.1002/jcc.540130812)  
**Citaciones: 5,876**  
**Resumen**: Primer desarrollo completo de WHAM. Establece las ecuaciones iterativas y demuestra aplicación a rotación de enlaces en butano.

**[2] Oshima H, Re S, Sugita Y. (2019)**  
"Replica-exchange umbrella sampling combined with Gaussian accelerated molecular dynamics for free-energy calculation of biomolecules"  
*Journal of Chemical Theory and Computation*, 15(10):5199-5208.  
DOI: [10.1021/acs.jctc.9b00501](https://doi.org/10.1021/acs.jctc.9b00501)  
**Citaciones: 63**  
**Resumen**: Método GaREUS que combina umbrella sampling con aceleración gaussiana y replica-exchange. Mejora convergencia 10-100x en sistemas proteicos.

**[3] Shirts MR, Chodera JD. (2008)**  
"Statistically optimal analysis of samples from multiple equilibrium states"  
*Journal of Chemical Physics*, 129:124105.  
DOI: [10.1063/1.2978177](https://doi.org/10.1063/1.2978177)  
**Citaciones: 1,847**  
**Resumen**: Desarrollo de MBAR, generalización de WHAM que maximiza likelihood de los datos. Proporciona errores analíticos.

**Artículos sobre Collective Variables**:

**[4] Thiede EH, Van Koten B, Weare J, Dinner AR. (2016)**  
"Eigenvector method for umbrella sampling enables error analysis"  
*arXiv:1607.03722* [physics.comp-ph]  
**Resumen**: Método para validar CVs usando análisis de eigenvectores. Permite identificar CVs óptimas.

**[5] Awasthi S, Nair NN. (2015)**  
"Exploring high-dimensional free energy landscapes: Temperature accelerated sliced sampling"  
*arXiv:1508.05181* [physics.chem-ph]  
**Resumen**: Técnica para manejar múltiples CVs simultáneamente sin explosión combinatoria.

**Artículos sobre WNK Kinases**:

**[6] Jonniya NA, Sk MF, Kar P. (2019)**  
"Investigating phosphorylation-induced conformational changes in WNK1 kinase by molecular dynamics simulations"  
*ACS Omega*, 4(17):17404-17416.  
DOI: [10.1021/acsomega.9b02368](https://doi.org/10.1021/acsomega.9b02368)  
**Citaciones: 46**  
**Resumen**: MD de 200 ns mostrando que fosforilación en T420 induce compactación del sitio activo de WNK1. Relevante para mecanismo de activación.

**[7] Zhang J, Siew K, Macartney T, et al. (2015)**  
"Critical role of the SPAK protein scaffold in regulating blood pressure in response to kidney injury"  
*Journal of the American Society of Nephrology*, 26(10):2367-2380.  
DOI: [10.1681/ASN.2014070672](https://doi.org/10.1681/ASN.2014070672)  
**Citaciones: 35**  
**Resumen**: Demuestra que dominio CCT de SPAK es esencial para señalización WNK-SPAK-NCC y regulación de presión arterial.

**Artículos sobre Motivo RFxV y Dominio CCT**:

**[8] Taylor SS, Meharena HS, Kornev AP. (2024)**  
"WNK kinases and the OSR1/SPAK kinases: Unraveling the structural basis for regulation"  
*Structure (London)*, 32(11):1975-1978.  
DOI: [10.1016/j.str.2024.09.007](https://doi.org/10.1016/j.str.2024.09.007)  
**Resumen**: Review reciente de Susan Taylor sobre bases estructurales de interacciones WNK-OSR1/SPAK vía motivo RFxV.

**[9] Taylor SS, Meharena HS, Kornev AP. (2025)**  
"WNK kinases and the OSR1/SPAK kinases: Unraveling the structural basis for regulation, Part II"  
*Structure (London)*, 33(1):1-4.  
DOI: [10.1016/j.str.2024.11.006](https://doi.org/10.1016/j.str.2024.11.006)  
**Resumen**: Continuación describiendo transmisión alostérica de señales osmóticas a través de WNK.

**Libros y Recursos de Referencia**:

- **Frenkel D, Smit B. (2001)**  
  *Understanding Molecular Simulation: From Algorithms to Applications*  
  Academic Press. ISBN: 978-0122673511  
  → Capítulo 7: Free Energy Calculations (incluye umbrella sampling detallado)

- **Leach AR. (2001)**  
  *Molecular Modelling: Principles and Applications*  
  Prentice Hall. ISBN: 978-0582382107  
  → Texto introductorio excelente para MD y forcefields

- **OpenMM Documentation**  
  [http://docs.openmm.org/](http://docs.openmm.org/)  
  → User Guide, Developer Guide, API Reference

- **pymbar Documentation**  
  [https://pymbar.readthedocs.io/](https://pymbar.readthedocs.io/)  
  → Tutoriales de WHAM/MBAR con ejemplos prácticos

**Bases de Datos Estructurales**:

- **Protein Data Bank (PDB)**  
  [https://www.rcsb.org/](https://www.rcsb.org/)  
  - 6CN9: WNK1 dominio kinasa (dímero auto-inhibido)
  - 2LRU: SPAK dominio CCT en complejo con péptido RFxV

- **AlphaFold Protein Structure Database**  
  [https://alphafold.ebi.ac.uk/](https://alphafold.ebi.ac.uk/)  
  - WNK1 humano (UniProt: Q9H4A3) predicción de estructura completa

---

### Apéndice C: Recursos Computacionales Requeridos

#### C.1 Hardware Recomendado

**Para Fase 2 (Piloto, 6 ventanas × 5 ns)**:

| Componente | Mínimo | Recomendado | Óptimo |
|---|---|---|---|
| **GPU** | NVIDIA GTX 1080 (8 GB) | NVIDIA RTX 3070 (8 GB) | NVIDIA RTX 3090 (24 GB) |
| **CPU** | Intel i5 (6 cores) | Intel i7 (8 cores) | AMD Threadripper (16+ cores) |
| **RAM** | 16 GB | 32 GB | 64 GB |
| **Almacenamiento** | 50 GB SSD | 200 GB NVMe SSD | 1 TB NVMe SSD |
| **Tiempo estimado** | 5-7 días | 2-3 días | 1 día |

**Para Fase 3 (Producción, 30 ventanas × 50 ns × 3 réplicas)**:

| Componente | Configuración Cluster |
|---|---|
| **Nodos** | 4-8 nodos GPU |
| **GPU por nodo** | 2-4 × NVIDIA A100 (40 GB) o RTX 3090 |
| **CPU por nodo** | 2 × AMD EPYC (32 cores cada uno) |
| **RAM por nodo** | 256 GB |
| **Interconexión** | InfiniBand HDR (200 Gbps) |
| **Almacenamiento** | 10 TB shared storage (Lustre o BeeGFS) |
| **Tiempo estimado** | 2 semanas (paralelo completo) |

---

#### C.2 Software Dependencies

**Sistema Operativo**:
- Linux (Ubuntu 20.04+, CentOS 8+) [Recomendado]
- macOS 11+ (Intel o Apple Silicon con Rosetta)
- Windows 10/11 con WSL2 (Windows Subsystem for Linux)

**Python Environment**:
```yaml
# environment.yml
name: bsm-lancad-env
channels:
  - conda-forge
  - omnia
dependencies:
  - python=3.10
  - openmm=8.0
  - cudatoolkit=11.8  # Para GPU NVIDIA
  - numpy=1.24
  - scipy=1.10
  - pandas=2.0
  - matplotlib=3.7
  - seaborn=0.12
  - jupyter=1.0
  - pymbar=4.0
  - mdtraj=1.9
  - nglview=3.0  # Visualización en notebook
  - pytest=7.4  # Testing
  - black=23.0  # Code formatting
  - click=8.1  # CLI utilities
```

**Instalación**:
```bash
# Crear environment
conda env create -f environment.yml

# Activar
conda activate bsm-lancad-env

# Verificar OpenMM con GPU
python -m openmm.testInstallation
# Debe mostrar: "OpenMM Version: 8.0... Platform: CUDA"
```

**Alternativa: pip (si Conda no disponible)**:
```bash
python -m venv venv_umbrella
source venv_umbrella/bin/activate  # Linux/Mac
# O en Windows: venv_umbrella\Scripts\activate

pip install -r requirements.txt
```

---

#### C.3 Estimaciones de Costo y Tiempo

**Fase 2 (Piloto)**:

| Métrica | Valor |
|---|---|
| **Costo GPU cloud** (ej. AWS p3.2xlarge, $3.06/hr) | ~$150 USD |
| **Costo GPU local** (depreciación RTX 3090) | $10-20 USD equivalente |
| **Tiempo humano** (preparación + análisis) | 20-40 horas |
| **Tiempo de reloj** | 2-5 días |

**Fase 3 (Producción)**:

| Métrica | Valor |
|---|---|
| **Costo GPU cloud** (AWS p3.8xlarge cluster, $12.24/hr) | ~$4,000 USD |
| **Costo GPU local** (cluster 4 × RTX 3090) | $200-400 USD equivalente |
| **Tiempo humano** | 80-120 horas |
| **Tiempo de reloj** | 2-4 semanas |

**Costo por nanosegundo de simulación**:
```
Hardware local (RTX 3090):
    Velocidad: ~150 ns/día para sistema ~50,000 átomos
    Costo depreciación: ~$0.50/día
    → $0.003 USD por nanosegundo
    
Cloud (p3.2xlarge, Tesla V100):
    Velocidad: ~100 ns/día
    Costo: $3.06/hr × 24 hr = $73.44/día
    → $0.73 USD por nanosegundo
    
    ⚠️ Cloud es ~240x más caro para proyectos largos!
    → Invertir en hardware local si proyecto es continuo
```

---

#### C.4 Almacenamiento y Backup

**Tamaños de archivo típicos**:

| Archivo | Tamaño (1 ventana, 50 ns) | Nota |
|---|---|---|
| **Trayectoria (.dcd)** | 500-800 MB | Coordenadas cada 10 ps |
| **Topología (.pdb)** | 2-5 MB | Una vez por sistema |
| **Checkpoint (.chk)** | 50-100 MB | Para reiniciar |
| **Log files** | 1-5 MB | Energías, T, P |
| **CV timeseries (.dat)** | 500 KB | Datos de análisis |

**Total para Fase 3** (30 ventanas × 3 réplicas):
```
Trayectorias: 90 × 700 MB = 63 GB
Checkpoints: 90 × 75 MB = 6.8 GB
Análisis: 90 × 500 KB = 45 MB
Topologías y misc: 1 GB

Total: ~71 GB
```

**Estrategia de backup**:
```
Tier 1 (datos crudos de simulación):
    → Trayectorias completas .dcd
    → Backup en almacenamiento frío (ej. AWS Glacier, $0.004/GB/mes)
    → Costo: 63 GB × $0.004 = $0.25/mes
    
Tier 2 (datos procesados):
    → CV timeseries, PMF, histogramas
    → GitHub LFS o almacenamiento institucional
    → Costo: Usualmente gratuito (<5 GB)
    
Tier 3 (código y scripts):
    → Repositorio Git (GitHub, GitLab)
    → Versionado completo
    → Costo: Gratuito
```

---

#### C.5 Checklist de Infraestructura

**Antes de empezar Fase 2**:

```
□ Hardware
  □ GPU con compute capability ≥ 6.0 (CUDA 11+)
  □ Driver NVIDIA actualizado (≥ 515.x)
  □ Al menos 20 GB espacio libre en disco rápido (SSD)
  
□ Software
  □ Conda o Mamba instalado
  □ Environment bsm-lancad-env creado y testeado
  □ OpenMM testInstallation exitoso (Platform: CUDA)
  □ pymbar importa sin errores
  
□ Datos
  □ Estructura PDB preparada (wnk1_monomer_prepared.pdb)
  □ Verificada con pdbfixer o pdb4amber
  □ Átomos de CV identificados (índices correctos)
  
□ Scripts
  □ Repositorio clonado
  □ Scripts de lanzamiento testeados (run_umbrella.ps1 o .sh)
  □ Pipeline sintético ejecutado exitosamente
  
□ Documentación
  □ README.md leído
  □ Este documento (UMBRELLA_SAMPLING_EXPLICADO.md) revisado
  □ Parámetros de config.py entendidos
```

---

#### C.6 Solución de Problemas Comunes

**Error: "Platform CUDA not found"**
```
Causa: OpenMM no detecta GPU

Soluciones:
1. Verificar driver NVIDIA: nvidia-smi
   → Debe mostrar GPU y CUDA version
   
2. Reinstalar OpenMM con CUDA:
   conda install -c conda-forge openmm cudatoolkit=11.8
   
3. Forzar plataforma CPU (más lento):
   platform = Platform.getPlatformByName('CPU')
```

**Error: "Out of Memory (OOM)" en GPU**
```
Causa: Sistema demasiado grande para VRAM

Soluciones:
1. Reducir tamaño de caja de agua (box_padding = 10 Å en lugar de 15)
2. Usar precisión mixta:
   platform.setPropertyDefaultValue('Precision', 'mixed')
3. Reducir batch_size (simular menos ventanas en paralelo)
4. Usar GPU con más memoria (ej. RTX 3090 24GB)
```

**Error: "Simulation unstable / NaN energies"**
```
Causa: Paso de tiempo muy grande o clashes en estructura

Soluciones:
1. Reducir timestep:
   integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 1*femtoseconds)
   
2. Minimizar energía más agresivamente:
   simulation.minimizeEnergy(maxIterations=10000)
   
3. Aumentar fase de equilibración (más tiempo, incremento gradual de T)
   
4. Verificar estructura con VMD:
   vmd wnk1_monomer_prepared.pdb
   → Buscar clashes (átomos solapados)
```

**Ventanas no se solapan en histogramas**
```
Causa: Constante de fuerza muy alta o spacing muy grande

Soluciones:
1. Reducir force_constant:
   k = 8.0 kcal/mol/Å² (en lugar de 12.0)
   
2. Reducir spacing:
   window_centers = np.linspace(8.0, 14.0, 30)  # Más ventanas
   
3. Simular más tiempo por ventana (mejor sampleo de colas)
```

---

### Agradecimientos

Este documento fue desarrollado por el equipo del **Laboratorio de Biofísica Computacional, UNAM-INN**, como parte del proyecto de caracterización termodinámica de quinasas WNK.

**Contribuciones**:
- Diseño de pipeline: [Nombres del equipo]
- Implementación de software: [Contribuidores]
- Validación biológica: [Asesores]
- Documentación técnica: Asistencia de GitHub Copilot

**Herramientas MCP utilizadas en este documento**:
- **Semantic Scholar MCP**: Búsqueda de papers, obtención de citaciones
- **arXiv MCP**: Búsqueda de preprints en física computacional
- **Byterover MCP**: Gestión de conocimiento del proyecto

**Agradecimiento especial** a:
- OpenMM team (Stanford, Memorial Sloan Kettering)
- pymbar developers (Chodera Lab)
- Comunidad de MD en GitHub y Stack Overflow

---

### Licencia y Cita

**Licencia**: Este documento y el código asociado están bajo **MIT License**.

**Cómo citar este trabajo**:
```
@techreport{umbrella_wnk_2025,
  title={Umbrella Sampling y la Quinasa WNK: De los Fundamentos F\'isicos a la Implementaci\'on Computacional},
  author={[Equipo UNAM-INN]},
  institution={Universidad Nacional Aut\'onoma de M\'exico, Instituto Nacional de Neurolog\'ia},
  year={2025},
  url={https://github.com/[repo]/knowledge-base/UMBRELLA_SAMPLING_EXPLICADO.md}
}
```

---

### Información de Contacto

**Para preguntas técnicas sobre el código**:
- GitHub Issues: [URL del repo]
- Email técnico: [correo@institucional.mx]

**Para colaboraciones científicas**:
- PI: [Nombre del Investigador Principal]
- Email: [correo@institucional.mx]
- Lab website: [URL]

---

### Historial de Versiones

| Versión | Fecha | Cambios |
|---|---|---|
| **v1.0** | 2025-01-XX | Documento inicial completo. Incluye Partes I-III, Apéndices A-C. |
| v0.5 | 2025-01-XX | Draft interno con Parte I y II. |
| v0.1 | 2025-01-XX | Outline y estructura. |

---

### Notas Finales

**Este documento es un trabajo vivo**. Si encuentras:
- Errores técnicos o conceptuales
- Secciones poco claras
- Referencias faltantes
- Sugerencias para mejoras

Por favor abre un **GitHub Issue** o contacta al equipo directamente.

**Recursos adicionales** en el repositorio:
```
knowledge-base/
├── UMBRELLA_SAMPLING_EXPLICADO.md  ← Estás aquí
├── GUIA_INSTALACION.md            ← Setup del environment
├── WNKTHERMODYNAMICS.MD            ← Biología profunda de WNK
├── papers/                         ← PDFs de referencias (a descargar)
└── recursos/
    ├── tutorial_alanina_EXPLICADO.ipynb  ← Tutorial introductorio MD
    └── simulacion_alanina.py             ← Ejemplo simple
```

**¡Buena suerte con tus simulaciones!** 🧬🖥️🚀

---

**FIN DEL DOCUMENTO**

---

