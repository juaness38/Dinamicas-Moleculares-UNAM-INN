#!/usr/bin/env bash
# Launch the WNK umbrella sampling demo in synthetic mode on Linux/macOS.
# Ensures the Conda environment exists before running the pipeline.

set -euo pipefail

ENV_NAME="${1:-bsm-lancad-env}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
ENV_FILE="${REPO_ROOT}/environment.yml"

CONDA_BIN="${CONDA_BIN:-}"

if [[ -z "${CONDA_BIN}" ]]; then
  if command -v conda >/dev/null 2>&1; then
    CONDA_BIN="$(command -v conda)"
  else
    for candidate in "${HOME}/miniconda3/bin/conda" "${HOME}/anaconda3/bin/conda"; do
      if [[ -x "${candidate}" ]]; then
        CONDA_BIN="${candidate}"
        break
      fi
    done
  fi
fi

if [[ -z "${CONDA_BIN}" ]]; then
  echo "[ERROR] Conda is not available in PATH. Set CONDA_BIN to your conda executable (e.g., export CONDA_BIN=\"${HOME}/miniconda3/bin/conda\")." >&2
  exit 1
fi

if ! "${CONDA_BIN}" env list | grep -q "^${ENV_NAME}[[:space:]]"; then
  echo "Creating Conda environment '${ENV_NAME}' from ${ENV_FILE}" >&2
  "${CONDA_BIN}" env create -f "${ENV_FILE}"
fi

echo "Running umbrella sampling demo (synthetic mode)" >&2
"${CONDA_BIN}" run -n "${ENV_NAME}" python -m Chronosfold.umbrella_suite.run_wnk_pipeline --synthetic
