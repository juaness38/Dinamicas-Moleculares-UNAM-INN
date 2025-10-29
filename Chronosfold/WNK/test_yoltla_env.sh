#!/bin/bash
# Test OpenMM environment in Yoltla
# Usage: bash test_yoltla_env.sh

echo "========================================"
echo "Testing OpenMM Environment in Yoltla"
echo "========================================"

# Activate conda environment
source activate openmm-env

echo ""
echo "Python version:"
python --version

echo ""
echo "Testing OpenMM..."
python << 'EOF'
import openmm
print(f"✓ OpenMM version: {openmm.__version__}")

from openmm import Platform
n_platforms = Platform.getNumPlatforms()
platforms = [Platform.getPlatform(i).getName() for i in range(n_platforms)]
print(f"✓ Platforms: {platforms}")
if 'CUDA' in platforms:
    print("  ✓ CUDA available!")
EOF

echo ""
echo "Testing pymbar..."
python << 'EOF'
try:
    import pymbar
    print(f"✓ pymbar version: {pymbar.__version__}")
    from pymbar import timeseries
    print("✓ timeseries module available")
except ImportError as e:
    print(f"❌ pymbar NOT installed: {e}")
    print("   Install with: conda install -c conda-forge 'pymbar>=4.0'")
EOF

echo ""
echo "Testing mdtraj..."
python << 'EOF'
try:
    import mdtraj
    print(f"✓ mdtraj version: {mdtraj.__version__}")
except ImportError as e:
    print(f"❌ mdtraj NOT installed: {e}")
EOF

echo ""
echo "========================================"
echo "Test complete!"
echo "========================================"
