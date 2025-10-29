# 🔧 Guía de Instalación - OpenMM y Herramientas de MD

Esta guía te ayudará a configurar tu entorno de trabajo para realizar dinámicas moleculares con OpenMM.

## 📋 Requisitos Previos

- Python 3.8 o superior
- Conda (Anaconda o Miniconda) - **RECOMENDADO**
- Git (opcional, para clonar repositorios)

## ⚡ Inicio rápido: demo de umbrella sampling

Usa los nuevos scripts para verificar tu instalación con el demo sintético de WNK.

- **Windows**: Ejecuta en PowerShell `scripts\bootstrap_windows.ps1`. Si Conda no está en `PATH`, indica la ruta explícitamente: `scripts\bootstrap_windows.ps1 -CondaPath "$env:USERPROFILE\miniconda3\Scripts\conda.exe"`.
- **Linux/macOS**: En terminal ejecuta `bash scripts/bootstrap_linux.sh`. Si Conda no está en `PATH`, define `CONDA_BIN`, por ejemplo `CONDA_BIN=$HOME/miniconda3/bin/conda bash scripts/bootstrap_linux.sh`.

Una vez creado el ambiente, puedes repetir el demo con:

```powershell
scripts\run_umbrella.ps1
```

o en Linux/macOS:

```bash
bash scripts/run_umbrella.sh
```

Los resultados (histogramas, series de tiempo y metadatos) quedan en `umbrella_results/wnk_pilot/`.

---

## 🚀 Opción 1: Instalación con Conda (RECOMENDADO)

Conda es la forma más sencilla y confiable de instalar OpenMM.

### Paso 1: Crear el ambiente virtual

```bash
# Crear ambiente con Python 3.11 según environment.yml
conda create -n bsm-lancad-env python=3.11

# Activar el ambiente
conda activate bsm-lancad-env
```

### Paso 2: Instalar OpenMM y herramientas principales

```bash
# Instalar OpenMM desde conda-forge
conda install -c conda-forge openmm

# Instalar herramientas de análisis
conda install -c conda-forge mdtraj mdanalysis

# Instalar herramientas científicas básicas
conda install -c conda-forge numpy scipy matplotlib pandas

# Instalar Jupyter para notebooks
conda install -c conda-forge jupyter jupyterlab

# Instalar herramientas de visualización
conda install -c conda-forge nglview
conda install -c conda-forge py3dmol
```

### Paso 3: Verificar la instalación

```python
# Ejecutar en Python o en un notebook
import openmm
import mdtraj
import MDAnalysis

print(f"OpenMM version: {openmm.__version__}")
print(f"MDTraj version: {mdtraj.__version__}")
print(f"MDAnalysis version: {MDAnalysis.__version__}")
print("✅ Instalación exitosa!")
```

### Paso 4: Verificar soporte GPU (opcional)

```python
from openmm import Platform

print("Plataformas disponibles:")
for i in range(Platform.getNumPlatforms()):
    platform = Platform.getPlatform(i)
    print(f"  - {platform.getName()}")
```

Plataformas esperadas:
- **Reference**: CPU básico (lento)
- **CPU**: CPU optimizado
- **CUDA**: GPU NVIDIA (más rápido)
- **OpenCL**: GPU AMD/Intel

---

## 🐍 Opción 2: Instalación con pip + venv

Si no puedes usar Conda, puedes usar pip (menos recomendado para OpenMM).

### Paso 1: Crear ambiente virtual

**En Windows (PowerShell):**
```powershell
# Crear el ambiente
python -m venv openmm-env

# Activar el ambiente
.\openmm-env\Scripts\Activate.ps1
```

**En Linux/Mac:**
```bash
# Crear el ambiente
python3 -m venv openmm-env

# Activar el ambiente
source openmm-env/bin/activate
```

### Paso 2: Instalar paquetes

```bash
# Actualizar pip
pip install --upgrade pip

# Instalar paquetes principales
pip install numpy scipy matplotlib pandas
pip install jupyter jupyterlab
pip install mdtraj MDAnalysis

# OpenMM con pip (puede requerir compilación)
pip install openmm
```

⚠️ **Nota**: La instalación de OpenMM con pip puede ser compleja en Windows. Se recomienda usar Conda.

---

## 📦 Archivo requirements.txt

Hemos creado un archivo `requirements.txt` para facilitar la instalación:

```bash
pip install -r requirements.txt
```

---

## 📦 Archivo environment.yml

Para usuarios de Conda, pueden usar el archivo `environment.yml`:

```bash
conda env create -f environment.yml
conda activate bsm-lancad-env
```

---

## 🧪 Script de Verificación

Ejecuta el script `verificar_instalacion.py` para confirmar que todo está bien:

```bash
python knowledge-base/recursos/verificar_instalacion.py
```

---

## 🐛 Problemas Comunes

### Problema 1: "conda: command not found"
**Solución**: Instala Anaconda o Miniconda desde:
- Anaconda: https://www.anaconda.com/download
- Miniconda: https://docs.conda.io/en/latest/miniconda.html

### Problema 2: OpenMM no encuentra GPU
**Solución**: 
- Verifica que tengas drivers NVIDIA actualizados
- Instala CUDA Toolkit compatible
- Reinstala OpenMM: `conda install -c conda-forge openmm`

### Problema 3: Error al importar MDAnalysis
**Solución**:
```bash
conda install -c conda-forge mdanalysis
```

### Problema 4: Jupyter no puede encontrar el kernel
**Solución**:
```bash
python -m ipykernel install --user --name bsm-lancad-env --display-name "OpenMM"
```

### Problema 5: nglview no se visualiza en Jupyter
**Solución**:
```bash
jupyter nbextension enable --py --sys-prefix nglview
jupyter labextension install @jupyter-widgets/jupyterlab-manager
```

---

## 📚 Paquetes Adicionales Útiles

### Para análisis avanzado:
```bash
conda install -c conda-forge pytraj
conda install -c conda-forge parmed
```

### Para visualización:
```bash
conda install -c conda-forge pymol-open-source
pip install py3Dmol
```

### Para manejo de estructuras:
```bash
conda install -c conda-forge biopython
conda install -c conda-forge rdkit
```

---

## 🔄 Actualizar el Ambiente

```bash
# Activar el ambiente
conda activate bsm-lancad-env

# Actualizar todos los paquetes
conda update --all

# O actualizar paquetes específicos
conda update openmm mdtraj mdanalysis
```

---

## 🗑️ Desinstalar el Ambiente

Si necesitas empezar de nuevo:

```bash
# Desactivar el ambiente
conda deactivate

# Eliminar el ambiente
conda env remove -n bsm-lancad-env
```

---

## 📝 Resumen de Comandos Rápidos

```bash
# CREAR Y ACTIVAR
conda create -n bsm-lancad-env python=3.11
conda activate bsm-lancad-env

# INSTALAR TODO
conda install -c conda-forge openmm mdtraj mdanalysis numpy scipy matplotlib pandas jupyter jupyterlab nglview

# VERIFICAR
python -c "import openmm; print(openmm.__version__)"

# USAR
jupyter lab
```

---

## 🎯 Checklist de Instalación

- [ ] Conda instalado
- [ ] Ambiente virtual creado
- [ ] OpenMM instalado
- [ ] MDTraj instalado
- [ ] MDAnalysis instalado
- [ ] Jupyter instalado
- [ ] Verificación exitosa
- [ ] GPU detectada (opcional)
- [ ] Primer notebook ejecutado
- [ ] Demo `run_wnk_pipeline --synthetic` completado

---

## 📞 Ayuda

Si tienes problemas, documenta el error en tu carpeta personal y discútelo en la próxima reunión del equipo.

**Recursos de ayuda:**
- OpenMM Forum: https://github.com/openmm/openmm/discussions
- MDAnalysis Docs: https://docs.mdanalysis.org/
- MDTraj Docs: http://mdtraj.org/

---

**Última actualización**: Octubre 20, 2025
