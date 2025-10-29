#!/bin/bash
# Install packages in drMD_wnk_umbrella environment
# Run this on Yoltla: bash install_packages.sh

echo "========================================="
echo "Installing packages for WNK umbrella sampling"
echo "========================================="

# Activate environment
source activate drMD_wnk_umbrella

echo "Environment activated: $CONDA_DEFAULT_ENV"
echo ""

# Install via conda-forge
echo "Installing OpenMM, pymbar, mdtraj..."
conda install -c conda-forge -y \
    openmm \
    pymbar \
    mdtraj \
    pdbfixer \
    numpy \
    scipy \
    pandas \
    matplotlib \
    seaborn \
    pyyaml \
    tqdm

echo ""
echo "========================================="
echo "Installation complete!"
echo "========================================="
echo ""
echo "Test the environment with:"
echo "  bash test_yoltla_env.sh"
