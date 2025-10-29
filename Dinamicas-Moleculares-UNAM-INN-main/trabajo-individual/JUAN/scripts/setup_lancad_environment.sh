#!/bin/bash
# Dinámicas Moleculares UNAM-INN - Preparación de entorno HPC (LANCAD)
# Autor: Juan Bustamante
# Fecha: 8 de octubre de 2025
# Uso: ejecutar una vez en la consola de LANCAD para preparar el entorno de trabajo

set -euo pipefail

PROJECT_NAME="Dinamicas-Moleculares-UNAM-INN"
REMOTE_ROOT="/LUSTRE/home/lancad/2025/l.100066/${PROJECT_NAME}"
WORKSPACE_DIR="${REMOTE_ROOT}/trabajo-individual/JUAN"
VENV_DIR="${WORKSPACE_DIR}/.venv"
LOG_DIR="${WORKSPACE_DIR}/logs"

underline() {
    printf '\n%s\n' "${1}"
    printf '%*s\n' "${#1}" '' | tr ' ' '-'
}

underline "🚀 Iniciando setup HPC para ${PROJECT_NAME}"

underline "📦 Configurando entorno de módulos"
module purge || true
module load python/3.11 || { echo "❌ No se pudo cargar python/3.11"; exit 1; }

underline "🧰 Creando entorno virtual"
if [[ ! -d "${VENV_DIR}" ]]; then
    python -m venv "${VENV_DIR}"
fi
source "${VENV_DIR}/bin/activate"
python -m pip install --upgrade pip

underline "📚 Instalando dependencias principales"
REQ_FILE="${WORKSPACE_DIR}/requirements.txt"
if [[ -f "${REQ_FILE}" ]]; then
    python -m pip install -r "${REQ_FILE}"
else
    echo "⚠️ No se encontró requirements.txt en ${WORKSPACE_DIR}. Instalando subconjunto base." 
    python -m pip install numpy scipy pandas matplotlib jupyter jupyterlab ipykernel MDAnalysis openmm
fi

underline "🔬 Paquetes adicionales recomendados"
python -m pip install --upgrade "nglview>=3.0.0" "py3Dmol>=1.8.0" "tqdm>=4.62.0"

TORCH_INDEX_URL="${TORCH_INDEX_URL:-https://download.pytorch.org/whl/cu121}"
PYG_INDEX_URL="${PYG_INDEX_URL:-https://pytorch-geometric.com/whl/torch-2.3.0+cu121.html}"

underline "🧠 Instalando PyTorch y PyTorch Geometric"
python -m pip install --upgrade torch torchvision torchaudio --index-url "${TORCH_INDEX_URL}"
python -m pip install --upgrade torch-scatter torch-sparse torch-cluster torch-spline-conv pyg-lib --extra-index-url "${PYG_INDEX_URL}"
python -m pip install --upgrade torch-geometric torch-geometric-temporal

REQ_BSM="${WORKSPACE_DIR}/requirements_bsm.txt"
if [[ -f "${REQ_BSM}" ]]; then
    underline "🧬 Instalando stack BSM adicional (requirements_bsm.txt)"
    python -m pip install -r "${REQ_BSM}"
fi

underline "⚡ Ajustando variables de rendimiento"
export OMP_NUM_THREADS="32"
export OPENMM_CPU_THREADS="32"
export PYTHONUNBUFFERED="1"

underline "📁 Creando estructura de carpetas"
mkdir -p "${LOG_DIR}" "${WORKSPACE_DIR}/data" "${WORKSPACE_DIR}/notebooks" "${WORKSPACE_DIR}/results"

underline "🧪 Verificación rápida de OpenMM"
python - <<'PYTEST'
import openmm
print("✅ OpenMM importado correctamente. Versión:", openmm.__version__)
PYTEST

underline "✅ Setup completado"
echo "Entorno virtual: ${VENV_DIR}"
echo "Para activar: source ${VENV_DIR}/bin/activate"
echo "Registra los cambios en tu bitácora personal (README.md)."
