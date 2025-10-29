# 🚀 WNK Umbrella Sampling - Deployment Guide for Yoltla HPC

**Target System**: Yoltla Cluster (UNAM)  
**User**: l.100066  
**Mission**: Deploy validated umbrella sampling pipeline for WNK kinase simulations

---

## ⚡ Quick Start (5 Minutes)

```bash
# 1. Connect to Yoltla
ssh yoltla

# 2. Navigate to WNK directory
cd ~/mica/Chronosfold/WNK

# 3. Create conda environment
conda env create -f environment_wnk_umbrella.yml

# 4. Activate environment
conda activate drMD_wnk_umbrella

# 5. Verify installation
python test_environment.py

# 6. Test on small system (alanine dipeptide)
cd tests
sbatch submit_alanine_test.sh

# 7. Monitor job
watch squeue -u l.100066
```

---

## 📋 Pre-Flight Checklist

### ✅ Before Running Production WNK Simulations

- [ ] **1. Conda environment created** (`drMD_wnk_umbrella`)
- [ ] **2. pymbar >= 4.0 installed** (verify: `python -c "import pymbar; print(pymbar.__version__)"`)
- [ ] **3. CUDA available** (verify: `python -c "from openmm import Platform; print([Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])"`)
- [ ] **4. CORRECTED MBAR script** (`analyze_umbrella_mbar_CORRECTED.py` ready)
- [ ] **5. Test system validated** (alanine dipeptide PMF matches literature)
- [ ] **6. Storage space checked** (`df -h ~/mica` - need ~100 GB free)
- [ ] **7. Umbrella windows generated** (`python generate_umbrella_windows.py --n-windows 20`)
- [ ] **8. SLURM script configured** (`submit_umbrella_hpc_48cores.sh` ready)

---

## 🔧 Detailed Setup Instructions

### Step 1: Upload Files to Yoltla

**From local machine** (Windows PowerShell):

```powershell
# Navigate to local WNK directory
cd C:\Users\busta\Downloads\Dinamicas-Moleculares-UNAM-INN-main\Chronosfold\WNK

# Upload conda environment file
scp environment_wnk_umbrella.yml yoltla:~/mica/Chronosfold/WNK/

# Upload corrected MBAR script
scp analyze_umbrella_mbar_CORRECTED.py yoltla:~/mica/Chronosfold/WNK/

# Upload validation report
scp ../../docs/validation/WNK_UMBRELLA_VALIDATION_REPORT.md yoltla:~/mica/Chronosfold/WNK/docs/
```

**Alternative**: Files already uploaded via FileZilla (as mentioned by user).

---

### Step 2: Create Conda Environment

**On Yoltla**:

```bash
# Load Anaconda module (if not in .bashrc)
module load anaconda/2021.04

# Create environment
cd ~/mica/Chronosfold/WNK
conda env create -f environment_wnk_umbrella.yml

# This will install:
# - OpenMM 7.7.0 (upgraded from 7.4.1)
# - pymbar 4.0+ (with timeseries module)
# - openmm-plumed (for future metadynamics)
# - All analysis tools (mdtraj, matplotlib, etc.)

# Expected time: 5-10 minutes
```

---

### Step 3: Verify Installation

**Create verification script**: `test_environment.py`

```python
#!/usr/bin/env python3
"""Verify drMD_wnk_umbrella environment"""

import sys

print("="*70)
print("ENVIRONMENT VERIFICATION")
print("="*70)

# Test 1: OpenMM
try:
    import openmm
    print(f"✓ OpenMM: {openmm.__version__}")
    
    from openmm import Platform
    platforms = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
    print(f"  Platforms: {platforms}")
    
    if 'CUDA' in platforms:
        print("  ✓ CUDA available")
    else:
        print("  ⚠️  CUDA not available - will use CPU (slower)")
except ImportError as e:
    print(f"❌ OpenMM import failed: {e}")
    sys.exit(1)

# Test 2: pymbar
try:
    import pymbar
    version = pymbar.__version__
    major = int(version.split('.')[0])
    
    print(f"✓ pymbar: {version}")
    
    if major >= 4:
        from pymbar import timeseries
        print("  ✓ timeseries module available")
    else:
        print(f"  ❌ pymbar {version} is too old - need >= 4.0")
        print("     Upgrade: pip install --upgrade pymbar")
        sys.exit(1)
except ImportError as e:
    print(f"❌ pymbar import failed: {e}")
    sys.exit(1)

# Test 3: PLUMED (optional)
try:
    from openmmplumed import PlumedForce
    print("✓ OpenMM-PLUMED: installed")
except ImportError:
    print("⚠️  OpenMM-PLUMED: not installed (optional for metadynamics)")

# Test 4: Analysis tools
try:
    import mdtraj
    import pandas as pd
    import matplotlib
    import seaborn
    print("✓ Analysis tools: mdtraj, pandas, matplotlib, seaborn")
except ImportError as e:
    print(f"⚠️  Some analysis tools missing: {e}")

print("\n" + "="*70)
print("ENVIRONMENT READY FOR PRODUCTION")
print("="*70)
```

**Run verification**:

```bash
conda activate drMD_wnk_umbrella
python test_environment.py
```

**Expected output**:
```
======================================================================
ENVIRONMENT VERIFICATION
======================================================================
✓ OpenMM: 7.7.0
  Platforms: ['Reference', 'CPU', 'CUDA', 'OpenCL']
  ✓ CUDA available
✓ pymbar: 4.0.3
  ✓ timeseries module available
✓ OpenMM-PLUMED: installed
✓ Analysis tools: mdtraj, pandas, matplotlib, seaborn

======================================================================
ENVIRONMENT READY FOR PRODUCTION
======================================================================
```

---

### Step 4: Test on Alanine Dipeptide (Validation)

**Why**: Alanine dipeptide is a **benchmark system** with known PMF from literature.

**Setup test**:

```bash
cd ~/mica/Chronosfold/WNK/tests

# Create test directory
mkdir -p alanine_test
cd alanine_test

# Download alanine dipeptide structure (public PDB)
wget https://files.rcsb.org/download/1ALA.pdb

# Or use local structure if available
```

**Create test SLURM script**: `submit_alanine_test.sh`

```bash
#!/bin/bash
#SBATCH --job-name=ala_umbrella_test
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=logs/alanine_test_%j.log

# Load modules
module load anaconda/2021.04

# Activate environment
source activate drMD_wnk_umbrella

# Run quick umbrella sampling test (5 windows, 500 ps each)
cd ~/mica/Chronosfold/WNK/tests/alanine_test

# Generate windows
python ../../generate_umbrella_windows.py \
    --n-windows 5 \
    --r-min 0.3 \
    --r-max 0.7 \
    --spring-constant 500 \
    --output alanine_windows

# Run each window (serial for test)
for i in {0..4}; do
    python ../../run_umbrella_window.py \
        --window $i \
        --steps 250000 \
        --platform CPU \
        --save-interval 1000
done

# Analyze with CORRECTED MBAR
python ../../analyze_umbrella_mbar_CORRECTED.py \
    --windows-dir alanine_windows \
    --output alanine_pmf \
    --temperature 300

echo "✓ Test complete"
```

**Submit test**:

```bash
mkdir -p logs
sbatch submit_alanine_test.sh

# Monitor
watch squeue -u l.100066

# Check logs
tail -f logs/alanine_test_*.log
```

**Expected test time**: ~2 hours (CPU mode)

**Validation**: Compare PMF with literature (Kentsis et al. 2004, Garcia & Sanbonmatsu 2001)

---

### Step 5: Generate WNK Umbrella Windows

**Once test passes**:

```bash
cd ~/mica/Chronosfold/WNK

# Generate 20 windows spanning interdomain distance
python generate_umbrella_windows.py \
    --n-windows 20 \
    --r-min 2.0 \
    --r-max 6.0 \
    --spring-constant 1000 \
    --output umbrella_windows

# This creates:
# umbrella_windows/
#   ├── window_00/ (r0=2.0 nm)
#   ├── window_01/ (r0=2.2 nm)
#   ├── ...
#   ├── window_19/ (r0=6.0 nm)
#   ├── windows_config.csv
#   └── atom_groups.txt
```

---

### Step 6: Submit Production SLURM Job Array

**Create SLURM submission script**: `submit_umbrella_production.sh`

```bash
#!/bin/bash
#SBATCH --job-name=wnk_umbrella
#SBATCH --partition=gpu  # Use GPU partition if available, else cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12  # 12 cores per window
#SBATCH --gres=gpu:1  # Request 1 GPU (if available)
#SBATCH --time=24:00:00  # 24 hours should be enough for 10 ns
#SBATCH --array=0-19  # 20 windows (job array)
#SBATCH --output=logs/wnk_umbrella_%A_%a.log

# Load modules
module load anaconda/2021.04
module load cuda/11.2  # Match CUDA version

# Activate environment
source activate drMD_wnk_umbrella

# Set CUDA compiler (if using CUDA)
export OPENMM_CUDA_COMPILER=/usr/local/cuda-11.2/bin/nvcc

# Navigate to work directory
cd ~/mica/Chronosfold/WNK

# Run umbrella window for this array task
# SLURM_ARRAY_TASK_ID = window index (0-19)

python run_umbrella_window.py \
    --window $SLURM_ARRAY_TASK_ID \
    --steps 5000000 \
    --dt 0.002 \
    --temp 310 \
    --platform CUDA \
    --save-interval 5000 \
    --log-interval 10000

echo "✓ Window $SLURM_ARRAY_TASK_ID complete"
```

**Submit production jobs**:

```bash
# Create logs directory
mkdir -p logs

# Submit job array (all 20 windows in parallel)
sbatch submit_umbrella_production.sh

# Monitor overall progress
watch squeue -u l.100066

# Check individual window logs
tail -f logs/wnk_umbrella_*_5.log  # Window 5 log
```

**Expected wall time**:
- With GPU: **6-12 hours** per window (10 ns @ ~100 ns/day)
- With CPU (48 cores): **24-48 hours** per window

**Parallelization**:
- 20 windows run **simultaneously** (job array)
- **Total wall time = time for slowest window** (~6-12 hours)

---

### Step 7: Monitor Progress

**Check job status**:

```bash
# All jobs
squeue -u l.100066

# Specific job array
squeue --job JOBID

# Job details
scontrol show job JOBID
```

**Check convergence** (while running):

```bash
# Every 2 hours, check CV sampling
cd ~/mica/Chronosfold/WNK

python check_sampling_progress.py --windows-dir umbrella_windows

# This will show:
# - Samples collected per window
# - CV range sampled
# - Estimated completion time
```

---

### Step 8: Analyze Results (CORRECTED MBAR)

**When all windows complete**:

```bash
cd ~/mica/Chronosfold/WNK

# Run CORRECTED MBAR analysis
python analyze_umbrella_mbar_CORRECTED.py \
    --windows-dir umbrella_windows \
    --output pmf_final \
    --temperature 310 \
    --n-bins 50 \
    --bootstrap 100 \
    --time-blocks 10

# This produces:
# pmf_final/
#   ├── pmf.png (PMF plot with error bars)
#   ├── pmf.pdf (publication quality)
#   ├── overlap_matrix.png (convergence diagnostic)
#   ├── pmf_data.csv (numerical results)
#   └── mbar_diagnostics.txt (overlap, autocorrelation, etc.)
```

**Key outputs**:
- **PMF curve**: Free energy vs interdomain distance
- **Error bars**: Statistical uncertainties (from MBAR + autocorrelation)
- **Overlap matrix**: Convergence diagnostic (should be > 0.03 for adjacent windows)

---

### Step 9: Validation & Quality Control

**Checks to perform**:

```bash
# 1. Check overlap matrix
cat pmf_final/mbar_diagnostics.txt | grep "overlap"
# All adjacent windows should have O_ij > 0.03

# 2. Check autocorrelation
cat pmf_final/mbar_diagnostics.txt | grep "statistical inefficiency"
# g should be < 50 (otherwise need longer simulation)

# 3. Check equilibration
cat pmf_final/mbar_diagnostics.txt | grep "equilibration"
# t0 should be < 50% of total simulation time

# 4. Visual inspection
# Download plots and check:
scp yoltla:~/mica/Chronosfold/WNK/pmf_final/pmf.png .
scp yoltla:~/mica/Chronosfold/WNK/pmf_final/overlap_matrix.png .
```

**Quality criteria** (from literature):
- ✅ Adjacent overlap > 0.03 (good sampling)
- ✅ PMF smooth (no jagged artifacts)
- ✅ Error bars < 2 kJ/mol (good precision)
- ✅ Equilibration < 50% of trajectory (efficient sampling)
- ✅ Statistical inefficiency g < 50 (decorrelated samples)

---

### Step 10: Cross-Validation with Metadynamics (Optional)

**To further validate results**:

```bash
# Run metadynamics on same CV
sbatch submit_metadynamics_validation.sh

# Compare PMFs
python compare_umbrella_metad.py \
    --umbrella pmf_final/pmf_data.csv \
    --metad metad_results/fes.dat \
    --output validation_report
```

**Expected**: PMFs should agree within error bars (1-2 kJ/mol).

---

## 📊 Expected Results

### System: WNK1 C-terminal (~50,000 atoms)

**Simulation Parameters**:
- 20 umbrella windows
- 10 ns per window
- Spring constant: 1000 kJ/mol/nm²
- Temperature: 310 K (37°C)
- Platform: CUDA (GPU)

**Expected PMF Features**:
1. **Compact state** (r ~ 2-3 nm): Low free energy (0-5 kJ/mol)
2. **Transition barrier** (r ~ 3.5-4 nm): Peak (10-20 kJ/mol)
3. **Extended state** (r ~ 5-6 nm): Moderate energy (5-15 kJ/mol)

**Performance**:
- **Throughput**: ~100 ns/day per window (CUDA)
- **Wall time**: 6-12 hours (parallel, 20 windows)
- **Data size**: ~100 GB total (trajectories)
- **Analysis time**: 10-30 minutes (MBAR)

---

## 🐛 Troubleshooting

### Issue 1: pymbar version too old

**Symptom**: `ImportError: cannot import name 'timeseries' from 'pymbar'`

**Solution**:
```bash
conda activate drMD_wnk_umbrella
pip install --upgrade pymbar
python -c "import pymbar; print(pymbar.__version__)"  # Should be >= 4.0
```

---

### Issue 2: CUDA not available

**Symptom**: `Platform.getPlatformByName('CUDA')` fails

**Solution**:
```bash
# Load CUDA module
module load cuda/11.2

# Set CUDA compiler
export OPENMM_CUDA_COMPILER=/usr/local/cuda-11.2/bin/nvcc

# Verify
python -c "from openmm import Platform; print([Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())])"
```

---

### Issue 3: Poor overlap between windows

**Symptom**: Overlap matrix shows O_ij < 0.03

**Solution**:
```bash
# Option 1: Add more windows
python generate_umbrella_windows.py --n-windows 30  # Instead of 20

# Option 2: Increase sampling time
# Edit submit script: --steps 10000000  # 20 ns instead of 10 ns

# Option 3: Adjust spring constants
# For steep PMF regions: increase k (more restraint)
# For flat PMF regions: decrease k (more flexibility)
```

---

### Issue 4: Simulation crashes

**Symptom**: Job fails with "Particle coordinate is nan"

**Solution**:
```bash
# 1. Check initial structure
python -c "from openmm.app import PDBFile; pdb = PDBFile('initial.pdb'); print('OK')"

# 2. Reduce timestep
# Edit run_umbrella_window.py: --dt 0.001  # 1 fs instead of 2 fs

# 3. Add energy minimization
# Ensure simulation.minimizeEnergy() is called before MD
```

---

## 📚 References

### OpenMM Documentation
- CustomCVForce: https://openmm.github.io/openmm-cookbook/latest/cookbook/custom_forces.html
- Best Practices: https://openmm.github.io/openmm-cookbook/latest/cookbook/best_practices.html

### pymbar Documentation
- MBAR Tutorial: https://pymbar.readthedocs.io/
- Timeseries: https://pymbar.readthedocs.io/en/master/timeseries.html

### Scientific Papers
- Li et al. (2022): MBAR error analysis - *J. Chem. Phys.*
- Awasthi et al. (2015): US-MetaD hybrid - *J. Comput. Chem.*
- Chodera (2016): Equilibration detection - *arXiv:1512.07858*

---

## 🎯 Success Criteria

Before considering the umbrella sampling **production-ready**:

- ✅ Conda environment installs without errors
- ✅ Test system (alanine dipeptide) PMF matches literature
- ✅ All 20 WNK windows complete without crashes
- ✅ Overlap matrix shows O_ij > 0.03 for adjacent windows
- ✅ PMF has smooth profile with error bars < 2 kJ/mol
- ✅ Autocorrelation properly handled (g < 50)
- ✅ Cross-validation with metadynamics (optional) agrees within errors

---

**Deployment Status**: ⏳ **READY FOR TESTING**

**Next Actions**:
1. ✅ Upload files to Yoltla (DONE via FileZilla)
2. ⏳ Create conda environment
3. ⏳ Run verification script
4. ⏳ Test on alanine dipeptide
5. ⏳ Deploy WNK production

**Estimated Time to Production**: **1-2 days** (including validation)

---

**Contact**:
- **HPC Support**: Dr. Aris Thorne (autonomous simulation)
- **Analysis Support**: Dr. Yuan Chen (free energy methods)
- **Validation**: Dr. Sofia Petrov (experimental benchmarks)

**Documentation**: See `WNK_UMBRELLA_VALIDATION_REPORT.md` for full technical details.
