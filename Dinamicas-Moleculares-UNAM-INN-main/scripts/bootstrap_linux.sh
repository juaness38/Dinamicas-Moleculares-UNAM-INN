#!/usr/bin/env bash
# Prepares the Conda environment and launches the umbrella demo on Linux/macOS.

set -euo pipefail

ENV_NAME="bsm-lancad-env"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
RUNNER="${REPO_ROOT}/scripts/run_umbrella.sh"

if [[ ! -f "${RUNNER}" ]]; then
  echo "[ERROR] Unable to locate ${RUNNER}" >&2
  exit 1
fi

if [[ -z "${CONDA_BIN:-}" ]]; then
  if [[ -x "${HOME}/miniconda3/bin/conda" ]]; then
    export CONDA_BIN="${HOME}/miniconda3/bin/conda"
  elif [[ -x "${HOME}/anaconda3/bin/conda" ]]; then
    export CONDA_BIN="${HOME}/anaconda3/bin/conda"
  fi
fi

if [[ -z "${CONDA_BIN:-}" ]] && ! command -v conda >/dev/null 2>&1; then
  cat >&2 <<'EOF'
[WARN] Conda is not available in PATH.
Set CONDA_BIN to your conda executable or install Miniconda from:
  https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
(use the macOS installer on macOS) and restart your shell before rerunning this script.
You can also run:
  CONDA_BIN=$HOME/miniconda3/bin/conda bash scripts/bootstrap_linux.sh
EOF
  exit 1
fi

bash "${RUNNER}" "${ENV_NAME}"
