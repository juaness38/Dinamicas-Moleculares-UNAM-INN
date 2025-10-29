# 🎯 Guía de Umbrella Sampling con OpenMM

## 📊 ¿Qué es Umbrella Sampling?

**Umbrella Sampling** es una técnica de *enhanced sampling* para calcular **perfiles de energía libre** (PMF - Potential of Mean Force) a lo largo de una **coordenada de reacción** (collective variable).

### Diferencias con Metadinámica

| Aspecto | **Metadinámica** (drMD ya tiene) | **Umbrella Sampling** |
|---------|----------------------------------|----------------------|
| **Estrategia** | Agrega gaussianas adaptativas que "rellenan" pozos de energía | Múltiples simulaciones con restraints armónicos fijos |
| **Ventanas** | Una simulación continua | N simulaciones independientes |
| **Post-processing** | PMF directo al finalizar | Requiere WHAM o MBAR para unir ventanas |
| **Convergencia** | Automática (explora todo el espacio) | Manual (debes cubrir el rango completo) |
| **Costo computacional** | 1 simulación larga | N simulaciones medianas |
| **Mejor para** | Exploración de espacio conformacional | Cálculo preciso de barreras energéticas |

---

## 🔧 Implementación

### 1. Arquitectura del código

```
openmm_calculator.py          # Base: preparación, force fields, MD básica
         ↓
umbrella_sampling_calculator.py   # Extiende con umbrella sampling
         ↓
[tu_script_de_analisis.py]   # Usa la calculadora + WHAM
```

### 2. Componentes clave

#### A. Collective Variable (CV)
La **coordenada de reacción** que quieres estudiar:

```python
cv_config = {
    'type': 'distance',      # Tipo de CV
    'atoms': [5, 42],        # Átomos Cα entre residuos de interés
    'params': {}
}
```

**Tipos soportados:**
- `'distance'`: Distancia entre 2 átomos (e.g., apertura de binding pocket)
- `'angle'`: Ángulo entre 3 átomos
- `'dihedral'`: Ángulo dihedro entre 4 átomos (e.g., rotación de cadena lateral)
- `'rmsd'`: RMSD respecto a estructura de referencia

#### B. Ventanas de Umbrella

```python
# Ejemplo: distancia de 5 a 20 Å con 15 ventanas
window_centers = np.linspace(5.0, 20.0, 15)

# Constante de fuerza del restraint (kcal/mol/Å²)
force_constant = 10.0  # Típico: 5-20 kcal/mol/Å²
```

**Reglas de oro:**
- **Overlap**: Las ventanas deben solaparse (histogramas compartidos)
- **Número**: Más ventanas = mejor convergencia pero más costo
- **Spacing**: Típicamente 0.5-2 Å para distancias

#### C. Potencial de Umbrella

Para cada ventana `i` con centro `ξ₀ᵢ`:

```
V_umbrella(ξ) = 0.5 * k * (ξ - ξ₀ᵢ)²
```

Donde:
- `ξ` = valor de la collective variable
- `ξ₀ᵢ` = centro de la ventana
- `k` = constante de fuerza (force_constant)

---

## 💻 Ejemplo de uso

### Caso: Apertura de binding pocket

```python
import asyncio
import numpy as np
from pathlib import Path
from umbrella_sampling_calculator import create_umbrella_calculator

async def analizar_apertura_pocket():
    """
    Calcula PMF de apertura de un binding pocket
    medido como distancia entre dos residuos clave
    """
    
    # 1. Configuración
    calculator = create_umbrella_calculator({
        'platform': 'CPU',  # o 'CUDA' si tienes GPU
        'force_field': 'amber99sb.xml',
        'temperature': 300.0
    })
    
    # 2. Define collective variable
    # Ejemplo: distancia entre Cα de residuos 15 y 87
    cv_config = {
        'type': 'distance',
        'atoms': [234, 1456],  # Índices reales de átomos Cα
        'params': {}
    }
    
    # 3. Define ventanas de umbrella
    # Rango: pocket cerrado (8 Å) a abierto (18 Å)
    window_centers = np.linspace(8.0, 18.0, 20)  # 20 ventanas
    
    # 4. Configuración de simulación
    force_constant = 15.0  # kcal/mol/Å²
    simulation_time_ps = 5000.0  # 5 ns por ventana
    
    # 5. Ejecutar umbrella sampling completo
    print("🎯 Iniciando Umbrella Sampling...")
    results = await calculator.run_full_umbrella_sampling(
        structure_file="protein_equilibrated.pdb",
        cv_config=cv_config,
        window_centers=window_centers,
        force_constant=force_constant,
        simulation_time_ps=simulation_time_ps,
        temperature=300.0,
        output_dir="umbrella_results"
    )
    
    print(f"\n✅ {results['metadata']['num_windows']} ventanas completadas")
    print(f"📊 Rango CV: {results['metadata']['cv_range']}")
    print(f"⏱️  Tiempo total: {results['metadata']['total_simulation_time_ps']/1000:.1f} ns")
    
    # 6. Estimación simple de PMF (para verificación rápida)
    cv_values, pmf = calculator.estimate_pmf_simple(results['wham_data'])
    
    print("\n📈 PMF estimado (usar WHAM para precisión):")
    for cv, energy in zip(cv_values[::5], pmf[::5]):  # Cada 5 puntos
        print(f"   CV={cv:.2f} Å → ΔG={energy:.2f} kcal/mol")
    
    # 7. Guardar para análisis con WHAM externo
    np.savetxt(
        "umbrella_results/pmf_estimate.dat",
        np.column_stack([cv_values, pmf]),
        header="CV(Angstrom) PMF(kcal/mol)"
    )
    
    return results

# Ejecutar
if __name__ == "__main__":
    asyncio.run(analizar_apertura_pocket())
```

---

## 📊 Análisis con WHAM (recomendado para producción)

### Opción 1: pyWHAM (Python)

```bash
pip install pywham
```

```python
import wham

# Cargar datos de ventanas
window_data = []
for result in umbrella_results['windows']:
    window_data.append({
        'center': result['window_center'],
        'kappa': result['force_constant'],
        'data': result['cv_values']
    })

# Ejecutar WHAM
pmf = wham.wham(window_data, temperature=300)
```

### Opción 2: GROMACS WHAM

Si tienes GROMACS instalado:

```bash
# Convertir datos a formato GROMACS
# Luego:
gmx wham -it umbrella_metadata.dat -if pullf_files.dat -o pmf.xvg -hist histogram.xvg
```

### Opción 3: pymbar (más moderno)

```bash
pip install pymbar
```

```python
import pymbar

# Construir matriz de energías
u_kn = construct_energy_matrix(umbrella_results)

# MBAR
mbar = pymbar.MBAR(u_kn, N_k)
pmf, pmf_uncertainty = mbar.computePMF()
```

---

## 🎨 Visualización de resultados

```python
import matplotlib.pyplot as plt

def plot_umbrella_analysis(results):
    """Visualiza resultados de umbrella sampling"""
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. Histogramas de todas las ventanas
    ax = axes[0, 0]
    for result in results['windows']:
        hist = result['cv_histogram']
        bin_centers = (np.array(hist['bin_edges'][:-1]) + 
                      np.array(hist['bin_edges'][1:])) / 2
        ax.plot(bin_centers, hist['counts'], alpha=0.6, 
                label=f"ξ₀={result['window_center']:.1f}")
    
    ax.set_xlabel("Collective Variable (Å)")
    ax.set_ylabel("Probability")
    ax.set_title("Histogramas de ventanas (verificar overlap)")
    ax.legend(fontsize=6, ncol=2)
    
    # 2. Convergencia de cada ventana
    ax = axes[0, 1]
    for result in results['windows']:
        cv_values = result['cv_values']
        running_mean = np.cumsum(cv_values) / np.arange(1, len(cv_values)+1)
        ax.plot(running_mean, alpha=0.6)
    
    ax.set_xlabel("Simulation time (frames)")
    ax.set_ylabel("Running average CV")
    ax.set_title("Convergencia de ventanas")
    
    # 3. PMF estimado
    ax = axes[1, 0]
    # Asumiendo que calculaste pmf con estimate_pmf_simple
    # o con WHAM
    ax.plot(cv_values, pmf, 'b-', linewidth=2)
    ax.fill_between(cv_values, pmf-0.5, pmf+0.5, alpha=0.3)  # Incertidumbre estimada
    ax.set_xlabel("Collective Variable (Å)")
    ax.set_ylabel("PMF (kcal/mol)")
    ax.set_title("Potential of Mean Force")
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    
    # 4. Mapa de calor de sampling
    ax = axes[1, 1]
    centers = [r['window_center'] for r in results['windows']]
    all_cv = [r['cv_values'] for r in results['windows']]
    
    # Crear heatmap de densidad
    for i, (center, cv_data) in enumerate(zip(centers, all_cv)):
        ax.scatter([center]*len(cv_data), cv_data, alpha=0.05, s=1, c=f'C{i%10}')
    
    ax.set_xlabel("Window center (Å)")
    ax.set_ylabel("Sampled CV (Å)")
    ax.set_title("Cobertura de muestreo")
    
    plt.tight_layout()
    plt.savefig("umbrella_analysis.png", dpi=300)
    print("📊 Gráficos guardados en: umbrella_analysis.png")
```

---

## ✅ Checklist de calidad

Antes de confiar en tus resultados, verifica:

- [ ] **Overlap suficiente**: Los histogramas de ventanas adyacentes se solapan
- [ ] **Convergencia**: Las medias de CV se estabilizan en cada ventana
- [ ] **Cobertura completa**: No hay gaps en el rango de CV
- [ ] **Constante k apropiada**: No demasiado débil (poco restraint) ni fuerte (poco sampling)
- [ ] **Tiempo suficiente**: Cada ventana debe tener >1000 muestras decorreladas
- [ ] **Equilibración**: Descarta primeros 10-20% de cada trayectoria

---

## 🚀 Casos de uso científicos

### 1. Binding/Unbinding de ligando
```python
cv_config = {
    'type': 'distance',
    'atoms': [ligand_COM, pocket_COM],  # Centros de masa
}
window_centers = np.linspace(2.0, 15.0, 25)  # Contacto → disociado
```

### 2. Transición conformacional
```python
cv_config = {
    'type': 'rmsd',
    'atoms': backbone_atoms,
    'params': {'reference_positions': open_conformation}
}
window_centers = np.linspace(0.0, 5.0, 20)  # Cerrado → abierto
```

### 3. Rotación de cadena lateral
```python
cv_config = {
    'type': 'dihedral',
    'atoms': [N, CA, CB, CG],  # χ1 angle
}
window_centers = np.linspace(-180, 180, 36)  # Full rotation
```

---

## 🔬 Comparación: Cuándo usar cada método

### Usa **Metadinámica** si:
✅ Exploración de múltiples estados conformacionales  
✅ No conoces bien el landscape energético  
✅ Quieres identificar estados metastables automáticamente  
✅ Tienes pocas CVs bien definidas (<3)

### Usa **Umbrella Sampling** si:
✅ Cálculo preciso de barrera energética conocida  
✅ Quieres control fino sobre el muestreo  
✅ Análisis de binding affinities  
✅ Necesitas incertidumbres bien estimadas (con MBAR)  
✅ Puedes paralelizar las ventanas (ideal para clusters)

### Usa **ambos** si:
✅ Primero metadinámica para explorar  
✅ Luego umbrella sampling para refinar barreras importantes

---

## 📚 Referencias

1. **Umbrella Sampling**: Torrie & Valleau (1977) *J. Comput. Phys.*
2. **WHAM**: Kumar et al. (1992) *J. Comput. Chem.*
3. **MBAR**: Shirts & Chodera (2008) *J. Chem. Phys.*
4. **OpenMM Biasing**: http://docs.openmm.org/latest/userguide/theory/04_standard_forces.html

---

## 🐛 Troubleshooting

### Problema: Ventanas no convergen
**Solución**: Aumentar `simulation_time_ps` (×2-5)

### Problema: No hay overlap entre ventanas
**Solución**: Aumentar densidad de ventanas o usar `force_constant` más débil

### Problema: PMF con ruido excesivo
**Solución**: 
1. Más tiempo de simulación
2. Más ventanas
3. Usar MBAR en lugar de WHAM simple

### Problema: WHAM no converge
**Solución**:
1. Verificar que todas las ventanas tienen datos
2. Revisar que `force_constant` sea consistente
3. Aumentar tolerancia de convergencia

---

## 💡 Tips de performance

1. **Paralelización**: Cada ventana es independiente → lanza en paralelo
   ```bash
   # En cluster con SLURM
   for i in {0..19}; do
       sbatch run_window_${i}.sh
   done
   ```

2. **GPU**: Usa `platform='CUDA'` para speedup de 10-100×

3. **Checkpointing**: Guarda estados intermedios para reinicios

4. **Adaptive sampling**: Añade ventanas extra donde PMF tiene pendientes altas

---

**¿Listo para correr tu primer Umbrella Sampling? 🎯**

Revisa `umbrella_sampling_calculator.py` y adapta el ejemplo a tu sistema!
