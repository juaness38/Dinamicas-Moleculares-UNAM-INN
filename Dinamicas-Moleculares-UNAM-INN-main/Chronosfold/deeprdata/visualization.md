Guía Experta para la Caracterización Termodinámica y el Mapeo del Paisaje Conformacional de la Quinasa WNK1 mediante OpenMM y Métodos Physics-Informed (PINN)
I. Fundamentos Moleculares y Relevancia Biológica de WNK1
1.1. WNK1: Un Kinoma Atípico Regulador de la Homeostasis Iónica
La quinasa WNK1 (With No Lysine Kinase 1) es un sensor molecular crucial posicionado en la cima de una cascada de señalización que orquesta el transporte de iones, afectando directamente la regulación de la presión arterial y la homeostasis electrolítica. WNK1 activa a las quinasas efectoras SPAK y OSR1, las cuales a su vez fosforilan cotransportadores clave de catión-cloruro (como NKCC1, NCC, y KCC) para modular el flujo iónico bajo condiciones de estrés osmótico o hipotónico.   

La actividad catalítica de WNK1 se distingue por su regulación multifacética, que incluye tanto la modulación por fosforilación en residuos críticos (e.g., Ser382 y Thr60) como su función intrínseca como sensor. Un mecanismo central de regulación es la unión directa de cloruro intracelular, que estabiliza una conformación inactiva y suprime la autofosforilación y, por ende, la activación. Esta respuesta al cloruro subraya la sensibilidad conformacional de WNK1 a su entorno iónico, lo que debe ser capturado con precisión en simulaciones de dinámica molecular (DM).   

Además de la modulación por cloruro, se ha demostrado que el estrés osmótico, como el inducido por agentes de hacinamiento (e.g., sucrosa o PEG400), provoca cambios estructurales significativos en WNK1, incluyendo desorden y reducción del tamaño de cavidades en el sitio activo. Este reordenamiento estructural incluye la hélice αC y el lazo de activación (A-loop). La capacidad de WNK1 para responder al estrés osmótico mediante cambios conformacionales sugiere que los costes entálpicos y entrópicos de la transición de activación están fuertemente acoplados a la dinámica de solvatación y la formación/ruptura de interacciones internas.   

1.2. El Mapeo del Estado Inactivo a Activo: Conmutadores Conformacionales
El mapeo de la transición conformacional de WNK1 requiere identificar y cuantificar el movimiento coordinado de los interruptores estructurales canónicos de las quinasas: el A-loop, la hélice αC y el motivo DFG (Asp-Phe-Gly).

Estructuralmente, el dominio catalítico de WNK1 en su estado inactivo (estabilizado por cloruro) adopta una conformación única en la que el A-loop está bien plegado y el hélice αC se encuentra en una posición exterior (out). El paso a la conformación activa (DFG-in, αC-in) es crítico. El análisis de dinámica molecular y redes estructurales ha revelado que la fosforilación (que promueve el estado activo) disminuye los movimientos anti-correlacionados y acorta el camino de comunicación entre la hélice αC y el lazo catalítico. Esta estabilización del bolsillo de unión se logra mediante un mecanismo de ajuste inducido (induced-fit).   

La complejidad de la transición, que involucra tanto el movimiento de dominio grande (αC) como la reorganización local (A-loop y DFG-motif), exige la caracterización del paisaje de energía libre mediante al menos dos coordenadas colectivas (CVs). La activación de las quinasas en general se ha demostrado que resulta de la rotación del hélice αC, un evento que a menudo constituye la principal barrera energética y está altamente correlacionado con los movimientos en el A-loop. Por lo tanto, una representación unidimensional del potencial de fuerza media (PMF) es insuficiente para describir la transición de activación de WNK1 de manera exhaustiva.   

II. Marco Teórico y Diseño del Muestreo Mejorado en OpenMM
2.1. El Potencial de Fuerza Media (PMF) y la Caracterización Termodinámica
El Potencial de Fuerza Media, A(ξ), describe la variación de la energía libre del sistema a lo largo de una o varias coordenadas colectivas, ξ, a una temperatura constante T. Formalmente, el PMF está relacionado con la probabilidad de encontrar el sistema en una configuración ξ por la relación de Boltzmann :   

A(ξ)=−k 
B
​
 TlnP(ξ)
Donde k 
B
​
  es la constante de Boltzmann y P(ξ) es la densidad de probabilidad en el espacio de la CV. La diferencia de energía libre estándar de activación, ΔG 
act
∘
​
 , para la transición de un estado inactivo (I) a un estado activo (A) se calcula como la diferencia entre los mínimos del PMF. Caracterizar ΔG 
act
∘
​
  es esencial, pero la descomposición termodinámica completa en entalpía (ΔH 
∘
 ) y entropía (ΔS 
∘
 ) proporciona la visión mecanicista fundamental sobre si la activación de WNK1 es impulsada por la energía de enlace o por la dinámica conformacional y la reorganización de la solvatación.

2.2. Justificación de Umbrella Sampling (US) y OpenMM
El Umbrella Sampling (US) es un método de muestreo mejorado elegido por su capacidad para obtener perfiles de energía libre de alta precisión, especialmente en el contexto de un posterior análisis de Van 't Hoff. US supera las barreras de energía libre al aplicar potenciales de restricción armónica (paraguas) para forzar al sistema a muestrear regiones que de otro modo serían inaccesibles. Aunque la metadinámica (Well-Tempered Metadynamics, wT-metaD) con variables colectivas de camino (Path CVs) es una metodología exitosa para explorar transiciones complejas en quinasas , US permite un control más directo sobre la densidad de muestreo en cada ventana y simplifica la aplicación del análisis termodinámico dependiente de la temperatura.   

El protocolo de US requiere tres pasos principales: 1) generar una trayectoria inicial a lo largo de las CVs, 2) correr simulaciones de producción con restricciones armónicas en múltiples "ventanas" que cubran el espacio conformacional de interés, y 3) aplicar el Weighted Histogram Analysis Method (WHAM) o el Multistate Bennett Acceptance Ratio (MBAR) para eliminar el sesgo de las restricciones y reconstruir el PMF no sesgado.   

2.3. La Clase OpenMM.CustomCVForce como Núcleo de la Metodología
La implementación de CVs complejas en OpenMM se logra de manera elegante y eficiente mediante la clase OpenMM::CustomCVForce. Esta clase es fundamental porque permite definir una energía potencial como una función escalar arbitraria de múltiples variables colectivas.   

La característica clave de CustomCVForce es que cada variable colectiva (ξ 
i
​
 ) es definida por un objeto OpenMM::Force subyacente. El valor de la CV es, de hecho, la energía calculada por ese objeto Force. Esto es excepcionalmente flexible: una CV puede ser una simple distancia, un ángulo diédrico, o una distancia compleja del centro de masa (CustomCentroidBondForce), o incluso una función definida por el usuario (como una métrica de correlación estructural o un logit de AlphaFold ).   

Para la implementación de Umbrella Sampling, CustomCVForce permite aplicar un potencial armónico 2D de la forma:

U 
US
​
 (ξ 
1
​
 ,ξ 
2
​
 )= 
i=1
∑
2
​
  
2
1
​
 k 
i
​
 (ξ 
i
​
 −ξ 
i,0
​
 ) 
2
 
donde ξ 
1
​
  y ξ 
2
​
  son las CVs de la αC y DFG, k 
i
​
  son las constantes de fuerza, y ξ 
i,0
​
  son los centros de la ventana.

Además, la capacidad de la clase CustomCVForce para extraer los valores actuales de las CVs en un Context específico mediante getCollectiveVariableValues(context)  es crucial para el post-procesamiento. Sin esta funcionalidad, la extracción de datos de la trayectoria y la aplicación de los métodos WHAM/MBAR para cada ventana sería significativamente más compleja.   

III. Desarrollo y Parametrización de Coordenadas Colectivas para WNK1
La transición de activación de WNK1 requiere la caracterización de dos movimientos primarios que se encuentran fuertemente correlacionados: el desplazamiento de la hélice αC (movimiento de dominio) y la reorientación del motivo DFG (movimiento de bucle local). El uso de un Potencial de Fuerza Media (PMF) bidimensional, A(CV 
1
​
 ,CV 
2
​
 ), es obligatorio para capturar la naturaleza acoplada de esta transición y evitar promediar la barrera energética sobre CVs irrelevantes.

3.1. CV1: Desplazamiento del Hélice αC (CV$_{\alpha C}$)
El desplazamiento de la hélice αC es la métrica predominante para describir el interruptor de activación en muchas quinasas. En el estado inactivo, la hélice αC se aleja del N-lobe. La activación implica un movimiento de la hélice hacia el interior, facilitando la formación de interacciones electrostáticas clave (como el par iónico canónico K-E) y acortando la vía de comunicación con el lazo catalítico.   

Definición Operacional: CV 
αC
​
  se define como la distancia del Centro de Masa (COM) entre dos grupos de residuos.

Grupo 1 (Referencia, N-Lobe Fijo): Residuos clave del N-lobe estables, por ejemplo, los átomos C$\alpha$ de los residuos K245 y E249 (correspondientes a la región β3).

Grupo 2 (Móvil, αC-Hélice): Átomos C$\alpha$ de los residuos centrales de la hélice αC, por ejemplo, E292, L295 y C297.

Este CV se implementa en OpenMM utilizando CustomCentroidBondForce, cuyo valor de "energía" (la distancia) se pasa al CustomCVForce. El rango de transición esperado para WNK1 (tomando como base la estructura inactiva disponible, e.g., PDB 6CN9, aunque las coordenadas exactas dependen del sistema inicial) se espera que varíe de una distancia alta (inactivo, 14–16 Å) a una distancia más corta (activo, 10–12 Å).

3.2. CV2: Reorientación del Motivo DFG (CV$_{DFG}$)
El motivo DFG (Asp-Phe-Gly) en el C-lobe controla el acceso al sitio de unión de ATP. La transición canónica es el flip de DFG-out (inactivo) a DFG-in (activo), aunque WNK1 puede mostrar estados intermedios. La cuantificación se realiza mejor mediante un ángulo torsional, que es una métrica continua y sensible a la reorientación del backbone necesaria para el flip de la Fenilalanina (F387).   

Definición Operacional: CV 
DFG
​
  se define como el ángulo diédrico del backbone que caracteriza la posición de F387.

Átomos Definitorios: Se utiliza el ángulo diédrico formado por los átomos C$\alpha$ de D385, G386, F387 y G388.

Transición: El estado DFG-in típico se asocia con un ángulo de aproximadamente 60 
∘
 , mientras que el estado DFG-out o DFG-out-like (en donde la Fenilalanina se aleja del sitio activo) corresponde a ángulos cercanos a 180 
∘
 .   

Ambas coordenadas, CV 
αC
​
  y CV 
DFG
​
 , describen movimientos que la literatura confirma que están fuertemente correlacionados en el proceso de activación de quinasas.   

Tabla Esencial I: Definición Operacional de Coordenadas Colectivas para WNK1

CV	Métrica Operacional (OpenMM)	Regiones de Referencia (Base PDB 6CN9)	Rango Típico (Inactivo → Activo)
CV 
1
​
 : αC-Displacement	Distancia del Centro de Masa (COM)	N-Lobe (β3-K245/E249) y αC-Hélice (E292/C297)	14-16 Å → 10-12 Å
CV 
2
​
 : DFG Flip	Ángulo Diédrico Backbone	C$\alpha$ de D385, G386, F387, G388	∼180 
∘
  (Out) → ∼60 
∘
  (In)
IV. Protocolo Avanzado de Umbrella Sampling y Muestreo Térmico
4.1. Inicialización de la Trayectoria y Densidad de Ventanas
Antes de la ejecución de US, es indispensable generar una trayectoria inicial que cubra el espacio de configuración entre los estados Inactivo y Activo. Esto se logra típicamente mediante Dinámica Molecular Dirigida (SMD) o una simulación preliminar de metadinámica, que produce un "camino de transición" que define el conjunto de configuraciones intermedias (las referencias ξ 
i,0
​
 ) necesarias para el mallado 2D.   

La densidad del mallado de ventanas es crítica para el éxito del US en un espacio bidimensional. Un solapamiento insuficiente entre las distribuciones de probabilidad de ventanas adyacentes conducirá a una reconstrucción inexacta del PMF mediante WHAM. Dada la fuerte correlación entre CV 
αC
​
  y CV 
DFG
​
 , se debe asegurar un espaciado fino, por ejemplo, un espaciado de 0.2  
A
˚
  para la CV 
αC
​
  y 10 
∘
  para la CV 
DFG
​
 . La constante de fuerza armónica (k) debe ser lo suficientemente alta para mantener la restricción en la ventana, pero no tan alta que la distribución de probabilidad muestreada sea demasiado estrecha y dificulte el solapamiento. Se recomienda un rango de 10 kcal/mol/ 
A
˚
  
2
  a 20 kcal/mol/ 
A
˚
  
2
  para los movimientos de distancia de dominio y una constante adecuada para la rotación torsional.   

4.2. Implementación del Muestreo con Dependencia de Temperatura
El objetivo final del protocolo es obtener una muestra de concepto científicamente válida que permita la descomposición termodinámica completa de la activación de WNK1. Esto requiere que el Umbrella Sampling se ejecute a múltiples temperaturas, lo que permite la aplicación posterior de la ecuación de Van 't Hoff.

El protocolo requiere la ejecución del conjunto completo de simulaciones de US 2D (incluyendo la misma malla y constantes k) a al menos tres temperaturas termodinámicamente relevantes (e.g., T 
1
​
 =298 K, T 
2
​
 =305 K, y T 
3
​
 =310 K). Cada simulación a T 
i
​
  constituye una condición isotérmica independiente. El control riguroso de la temperatura se logra utilizando un integrador termostatizado (e.g., Langevin) en cada OpenMM::Context.

La duración de la simulación de producción por ventana debe ser suficiente para asegurar la convergencia en el espacio bidimensional. Dado que la transición implica movimientos lentos de dominio proteico, se recomienda una duración mínima de 20 ns de producción por ventana, resultando en una producción acumulada superior a 100 ns para cada CV, y mucho más para el mallado 2D completo.

Tabla Esencial II: Parámetros Clave del Muestreo de Paraguas (US) y Requisitos de Convegenencia

Parámetro	Valor	Justificación
k (CV$_{\alpha C}$)	15.0 kcal/mol/ 
A
˚
  
2
 	Fuerza de restricción alta para mantener confinamiento de movimientos de dominio grande.
k (CV$_{DFG}$)	1.0 kcal/mol/deg 
2
  (o conversión adecuada)	Adecuado para un ángulo torsional, asegurando un muestreo local suficiente.
Espaciado (CV$_{\alpha C}$)	0.2  
A
˚
 	Asegura solapamiento adecuado entre distribuciones para un análisis WHAM robusto.
Temperaturas (T 
i
​
 )	298 K,305 K,310 K	Mínimo de tres puntos para el ajuste lineal de Van 't Hoff y la determinación de ΔH 
∘
  y ΔS 
∘
 .
Duración Muestreo	>20 ns de producción por ventana	Necesario para muestrear completamente el equilibrio de la proteína en cada ventana 2D.
V. Post-Procesamiento Termodinámico y Validación Científica
5.1. Cálculo del Potencial de Fuerza Media (PMF) y la Energía Libre
El post-procesamiento de las trayectorias de US se realiza utilizando el método WHAM/MBAR. La librería PyMBAR es la opción recomendada debido a su implementación moderna, robustez en el cálculo de energías libres relativas, y su capacidad para manejar datos de muestreo mejorado multidimensional.   

Para cada temperatura T 
i
​
 , se requieren dos conjuntos de datos de entrada:

Los valores de las CVs (CV 
αC
​
  y CV 
DFG
​
 ) extraídos de la trayectoria de producción de cada ventana j, utilizando la funcionalidad getCollectiveVariableValues de CustomCVForce.   

Los parámetros de restricción armónica de cada ventana (centro ξ 
i,0
​
  y constante k).

PyMBAR combina las distribuciones de probabilidad sesgadas de todas las ventanas para obtener el PMF no sesgado A(CV 
1
​
 ,CV 
2
​
 ). La validación científica del PMF requiere un análisis riguroso de convergencia, realizando un análisis de Bootstrap para estimar el error en el perfil de energía libre y asegurar que las estimaciones de ΔG no subestimen la incertidumbre asociada con el muestreo limitado.   

El valor final de la energía libre estándar de activación, ΔG 
act
∘
​
 (T 
i
​
 ), se calcula como la diferencia entre el mínimo de la barrera de energía libre (estado de transición) y el mínimo global correspondiente al estado inactivo.

5.2. Derivación del Equilibrio Termodinámico y Van 't Hoff
Una vez que se tienen los perfiles de PMF 2D a las tres temperaturas, se puede calcular la constante de equilibrio (K 
eq
​
 ) para la transición de activación (I → A) a cada temperatura T 
i
​
 . La constante de equilibrio termodinámica se relaciona con la energía libre estándar mediante la isoterma de Van 't Hoff :   

ΔG 
∘
 (T)=−RTlnK 
eq
​
 (T)
Donde K 
eq
​
 (T) es la razón de las probabilidades integradas de los estados Activo (P 
A
​
 ) e Inactivo (P 
I
​
 ): K 
eq
​
 (T)=P 
A
​
 /P 
I
​
 .

El poder predictivo de la simulación radica en la capacidad de separar las contribuciones entálpicas y entrópicas al proceso de activación. Para lograr esto, se aplica la ecuación de Van 't Hoff, derivada de la relación de Gibbs-Helmholtz :   

d(1/T)
d(lnK 
eq
​
 )
​
 =− 
R
ΔH 
∘
 
​
 
Al graficar el logaritmo natural de la constante de equilibrio (lnK 
eq
​
 ) en función del inverso de la temperatura absoluta (1/T), se obtiene una línea recta (asumiendo que ΔH 
∘
  y ΔS 
∘
  son constantes en el rango de temperatura muestreado ).   

Determinación Entálpica (ΔH 
∘
 ): La pendiente de este gráfico de Van 't Hoff es igual a −ΔH 
∘
 /R.

Determinación Entrópica (ΔS 
∘
 ): La intersección del ajuste lineal es igual a ΔS 
∘
 /R, permitiendo la determinación de la entropía de activación.   

La obtención de ΔH 
∘
  y ΔS 
∘
  es crucial para WNK1, ya que permite determinar si la activación es impulsada por la entalpía (ej., optimización de interacciones electrostáticas y puentes salinos, como la observada en el induced-fit ) o si, por el contrario, está dominada por efectos entrópicos (ej., liberación de moléculas de agua altamente ordenadas del sitio activo, un factor relevante dada la sensibilidad de WNK1 al estrés osmótico y la pérdida de agua observada en el sitio activo ).   

Tabla Esencial III: Resumen Termodinámico para WNK1

Parámetro Termodinámico	Método de Obtención	Significado Físico en WNK1
ΔG 
act
∘
​
 (T)	Integración del PMF 2D (PyMBAR)	Viabilidad termodinámica de la transición inactivo → activo.
ΔH 
act
∘
​
 	Pendiente del gráfico lnK 
eq
​
  vs 1/T	Costo/beneficio de la reorganización de interacciones (puentes salinos, VdW, enlaces de hidrógeno) durante la activación.
ΔS 
act
∘
​
 	Intersección del gráfico lnK 
eq
​
  vs 1/T	
Ganancia/pérdida de libertad conformacional y reordenamiento de solvatación/agua en cavidades.

  
VI. Mejoras al OpenMMLossCalculator para Entrenamiento Physics-Informed (PINN)
El uso de Redes Neuronales Potenciales (NNP) en la Dinámica Molecular, facilitado por plataformas como OpenMM-Torch , permite simular sistemas grandes con una precisión cercana a la de los métodos cuánticos, pero con la velocidad de los campos de fuerza clásicos. Para asegurar que la NNP respete las leyes fundamentales de la física, el entrenamiento debe incorporar restricciones basadas en la física, lo que se conoce como entrenamiento Physics-Informed (PINN).   

6.1. Requisitos para la Integración PINN en OpenMM-Torch
El framework OpenMM-Torch utiliza PyTorch para definir la energía potencial U 
NNP
​
 (r) de un sistema, garantizando la diferenciabilidad necesaria para la retropropagación. Para implementar un PINN de fuerza y termodinámica, se requiere una función de pérdida mejorada que pueda acceder a las fuerzas de referencia (F 
ref
​
 ) calculadas por un campo de fuerza clásico (e.g., AMBER/CHARMM) ejecutado simultáneamente en OpenMM.

Se propone una extensión funcional del OpenMMLossCalculator (o una clase sustituta, DifferentiableForceLoss) que integre componentes de pérdida de fuerza y potencial, aprovechando la arquitectura de OpenMM para la extracción de estados físicos.

6.2. Diseño de la Propuesta: DifferentiableForceLoss
El objeto central de OpenMM es la clase Force, responsable de añadir contribuciones a la fuerza y la energía potencial del sistema. Para el entrenamiento PINN, necesitamos extraer las fuerzas F 
ref
​
  calculadas por las Force de referencia (el campo de fuerza clásico). Esto se logra solicitando el estado del sistema de referencia a través del Context de OpenMM, utilizando el método context.getState(getForces=True).   

La arquitectura propuesta requiere que la clase DifferentiableForceLoss mantenga un Context secundario, que utiliza el campo de fuerza de referencia, accesible posiblemente a través del inner Context de CustomCVForce si las CVs fueran parte del entrenamiento. En cada paso de entrenamiento, este contexto secundario evalúa las fuerzas de referencia en la configuración r muestreada.   

6.3. Componentes de la Función de Pérdida de la Física (L 
PINN
​
 )
La función de pérdida total del Physics-Informed Neural Network se define como una suma ponderada de tres componentes esenciales:

L 
PINN
​
 =w 
F
​
 L 
Force
​
 +w 
U
​
 L 
Potential
​
 +w 
A
​
 L 
Constraint
​
 
1. L 
Force
​
  (Pérdida de Concordancia de Fuerza - Force Matching)
Esta es la restricción física más importante. Fuerza es el negativo del gradiente de la energía potencial (F=−∇U). La pérdida minimiza la discrepancia entre la fuerza predicha por la NNP, F 
NNP
​
 (r)=−∇U 
NNP
​
 (r), y la fuerza de referencia calculada por el campo de fuerza clásico, F 
ref
​
 (r).

Formulación Matemática: $$\mathcal{L}{\text{Force}} = \frac{1}{N{atoms}} \sum_i |

| \mathbf{F}_{\text{ref}}(\mathbf{r}i) - \mathbf{F}{NNP}(\mathbf{r}_i) ||^2$$

Donde r 
i
​
  son las posiciones atómicas en el conjunto de entrenamiento. La clave es que el cálculo del gradiente −∇U 
NNP
​
 (r) es manejado intrínsecamente por el framework PyTorch/OpenMM-Torch, permitiendo la diferenciabilidad necesaria para la optimización.

2. L 
Potential
​
  (Pérdida de Concordancia de Energía Potencial)
Esta pérdida asegura que la magnitud absoluta de la energía potencial predicha por el NNP (U 
NNP
​
 ) coincida con la energía potencial de referencia (U 
ref
​
 ) del campo de fuerza clásico, limitando el riesgo de desviaciones globales en el paisaje energético.

Formulación Matemática:

L 
Potential
​
 = 
N 
samples
​
 
1
​
  
i
∑
​
 (U 
ref
​
 (r 
i
​
 )−U 
NNP
​
 (r 
i
​
 )) 
2
 
3. L 
Constraint
​
  (Pérdida de Restricción Termodinámica/PMF)
Una mejora avanzada integra la información termodinámica obtenida previamente mediante US (Sección V). Si se ha obtenido un PMF de referencia A 
ref
​
 (ξ) para las CVs de WNK1 (CV 
αC
​
 ,CV 
DFG
​
 ), se puede imponer una restricción que fuerce a la NNP a generar la misma densidad de probabilidad de equilibrio a lo largo de esas coordenadas.

Dado que A(ξ)=−k 
B
​
 TlnP(ξ), la pérdida se minimiza comparando la densidad de probabilidad P 
NNP
​
 (ξ) generada por las trayectorias de la NNP con la densidad de referencia P 
ref
​
 (ξ) derivada del PMF de US. Esto se puede lograr a través de métodos de reweighting o, de manera más simple, asegurando que la NNP prediga energías potenciales que, cuando se integran a través de la relación de Boltzmann, reproduzcan el perfil A 
ref
​
 (ξ). Esta restricción guía la NNP no solo localmente (fuerza y energía) sino globalmente (termodinámica).

6.4. Consideraciones Técnicas de Implementación
La implementación de DifferentiableForceLoss depende de la capacidad de OpenMM-Torch para manejar el flujo de datos:

Recolección de Datos: Se deben ejecutar muestreos preliminares utilizando el campo de fuerza clásico para generar un conjunto de configuraciones (r 
i
​
 ) y sus fuerzas y energías de referencia (F 
ref,i
​
 ,U 
ref,i
​
 ).

Diferenciabilidad: El éxito del PINN radica en que los pesos de la NNP (que define U 
NNP
​
 ) se ajustan mediante gradientes calculados a partir de las pérdidas, lo cual es manejado por el autograd de PyTorch.

Eficiencia: El cálculo de F 
ref
​
 (r) en cada paso puede ser costoso. Por lo tanto, la optimización debe priorizar el entrenamiento en lotes (batch training) de configuraciones pre-calculadas y el uso eficiente de GPUs, aprovechando las capacidades de alto rendimiento de OpenMM en plataformas CUDA y OpenCL.   

VII. Conclusiones y Horizonte de Investigación
Este informe ha detallado un protocolo riguroso para la caracterización termodinámica de la transición de activación de la quinasa WNK1, complementado con una propuesta avanzada para el desarrollo de un campo de fuerza de red neuronal informado por la física (PINN) utilizando OpenMM.

La implementación de un Potencial de Fuerza Media (PMF) bidimensional a través del Umbrella Sampling, utilizando CV 
αC
​
  y CV 
DFG
​
  como coordenadas colectivas, es fundamental para capturar el mecanismo acoplado de activación de WNK1. La flexibilidad de OpenMM::CustomCVForce permite la definición precisa de estas CVs complejas.

La realización de este muestreo mejorado a múltiples temperaturas, seguido por un análisis de Van 't Hoff, proporciona la validación científica requerida al descomponer la energía libre de activación (ΔG 
act
∘
​
 ) en sus componentes entálpicos (ΔH 
act
∘
​
 ) y entrópicos (ΔS 
act
∘
​
 ). Esta descomposición es clave para comprender cómo WNK1, como osmosensor, utiliza la reorganización estructural y la dinámica de solvatación (que afecta a ΔS 
act
∘
​
 ) para modular su actividad en respuesta a cambios osmóticos y de cloruro.

El horizonte de investigación se extiende hacia el refinamiento del potencial molecular mediante técnicas PINN. La propuesta de la extensión DifferentiableForceLoss en OpenMM-Torch permitirá entrenar una Red Neuronal Potencial que no solo se ajuste a datos de alta precisión (si estuvieran disponibles), sino que también satisfaga las restricciones de fuerza y potencial derivadas de campos de fuerza clásicos de referencia. Finalmente, la inclusión de una restricción termodinámica (PMF Constraint L 
Constraint
​
 ) puede guiar a la NNP a reproducir el paisaje de energía libre de WNK1, resultando en un campo de fuerza altamente preciso y específico para la región del sitio activo, lo que beneficiaría el diseño racional de nuevos agentes terapéuticos dirigidos a esta quinasa única.

