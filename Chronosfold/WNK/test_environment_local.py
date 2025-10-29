#!/usr/bin/env python3
"""
🔬 WNK Umbrella Sampling Environment Verification

Tests all critical components before production deployment.
Run this after creating the conda environment.

Usage:
    conda activate drMD_wnk_umbrella
    python test_environment.py
"""

import sys
import os

def print_header(title):
    print("\n" + "="*70)
    print(title)
    print("="*70)

def test_openmm():
    """Test OpenMM installation and platforms"""
    print_header("TESTING OPENMM")
    
    try:
        import openmm
        print(f"✓ OpenMM version: {openmm.__version__}")
        
        # Check if version is recent enough
        version_parts = openmm.__version__().split('.')
        major = int(version_parts[0])
        minor = int(version_parts[1])
        
        if major < 7 or (major == 7 and minor < 7):
            print(f"  ⚠️  OpenMM {openmm.__version__()} is outdated")
            print("     Recommended: >= 7.7.0 for latest features")
        
        # Test platforms
        from openmm import Platform
        n_platforms = Platform.getNumPlatforms()
        platforms = [Platform.getPlatform(i).getName() for i in range(n_platforms)]
        
        print(f"✓ Available platforms: {platforms}")
        
        # Check for GPU platforms
        if 'CUDA' in platforms:
            print("  ✓ CUDA available (GPU acceleration enabled)")
            cuda_platform = Platform.getPlatformByName('CUDA')
            print(f"    Properties: {cuda_platform.getPropertyNames()}")
        elif 'OpenCL' in platforms:
            print("  ⚠️  OpenCL available (GPU acceleration with lower performance)")
        else:
            print("  ⚠️  No GPU platforms available - will use CPU (slower)")
        
        return True
        
    except Exception as e:
        print(f"❌ OpenMM test FAILED: {e}")
        return False

def test_pymbar():
    """Test pymbar installation and version"""
    print_header("TESTING PYMBAR")
    
    try:
        import pymbar
        version = pymbar.__version__
        print(f"✓ pymbar version: {version}")
        
        # Check version
        major = int(version.split('.')[0])
        
        if major < 4:
            print(f"  ❌ pymbar {version} is TOO OLD")
            print("     CRITICAL: Need >= 4.0 for timeseries module")
            print("     Upgrade: pip install --upgrade pymbar")
            return False
        
        # Test timeseries module
        try:
            from pymbar import timeseries
            print("  ✓ timeseries module available")
            
            # Test functions
            test_data = [1.0, 1.1, 1.0, 0.9, 1.0] * 100  # Fake correlated data
            t0, g, Neff = timeseries.detectEquilibration(test_data)
            print(f"    Test: detectEquilibration() works (t0={t0}, g={g:.1f}, Neff={Neff:.0f})")
            
        except ImportError:
            print("  ❌ timeseries module NOT FOUND")
            print("     Check pymbar installation")
            return False
        
        return True
        
    except Exception as e:
        print(f"❌ pymbar test FAILED: {e}")
        return False

def test_plumed():
    """Test OpenMM-PLUMED (optional)"""
    print_header("TESTING OPENMM-PLUMED (OPTIONAL)")
    
    try:
        from openmmplumed import PlumedForce
        print("✓ OpenMM-PLUMED installed")
        print("  Can use metadynamics and advanced collective variables")
        return True
    except ImportError:
        print("⚠️  OpenMM-PLUMED not installed")
        print("   This is OPTIONAL for umbrella sampling")
        print("   Install for metadynamics: conda install -c conda-forge openmm-plumed")
        return False

def test_analysis_tools():
    """Test analysis dependencies"""
    print_header("TESTING ANALYSIS TOOLS")
    
    tools = {
        'numpy': 'Numerical arrays',
        'scipy': 'Scientific computing',
        'pandas': 'Data analysis',
        'matplotlib': 'Plotting',
        'seaborn': 'Statistical visualization',
        'mdtraj': 'Trajectory analysis',
        'pdbfixer': 'PDB structure fixing'
    }
    
    all_ok = True
    
    for tool, description in tools.items():
        try:
            exec(f"import {tool}")
            version = eval(f"{tool}.__version__")
            print(f"✓ {tool:12s} {version:10s} ({description})")
        except ImportError:
            print(f"❌ {tool:12s} NOT FOUND ({description})")
            all_ok = False
        except AttributeError:
            print(f"✓ {tool:12s} installed    ({description})")
    
    return all_ok

def test_custom_forces():
    """Test OpenMM custom forces (critical for umbrella sampling)"""
    print_header("TESTING OPENMM CUSTOM FORCES")
    
    try:
        from openmm import CustomCentroidBondForce, CustomCVForce
        print("✓ CustomCentroidBondForce available")
        print("✓ CustomCVForce available")
        
        # Create a simple test force
        force = CustomCentroidBondForce(2, "0.5*k*(distance(g1,g2)-r0)^2")
        force.addPerBondParameter("k")
        force.addPerBondParameter("r0")
        
        print("  ✓ Umbrella sampling force creation works")
        return True
        
    except Exception as e:
        print(f"❌ Custom forces test FAILED: {e}")
        return False

def test_file_access():
    """Test that we can access expected directories"""
    print_header("TESTING FILE SYSTEM ACCESS")
    
    # Expected paths
    paths_to_check = [
        ('Current directory', '.'),
        ('Home directory', os.path.expanduser('~')),
        ('MICA workspace', os.path.expanduser('~/mica')),
    ]
    
    all_ok = True
    
    for name, path in paths_to_check:
        if os.path.exists(path):
            print(f"✓ {name:20s} accessible: {path}")
        else:
            print(f"⚠️  {name:20s} NOT FOUND: {path}")
            all_ok = False
    
    # Check write permissions
    try:
        test_file = 'test_write_permissions.tmp'
        with open(test_file, 'w') as f:
            f.write('test')
        os.remove(test_file)
        print("✓ Write permissions OK")
    except Exception as e:
        print(f"❌ Cannot write to current directory: {e}")
        all_ok = False
    
    return all_ok

def print_summary(results):
    """Print summary of all tests"""
    print_header("SUMMARY")
    
    all_tests = [
        ('OpenMM', results['openmm']),
        ('pymbar', results['pymbar']),
        ('PLUMED', results['plumed']),
        ('Analysis tools', results['analysis']),
        ('Custom forces', results['custom_forces']),
        ('File system', results['filesystem']),
    ]
    
    passed = sum(1 for _, result in all_tests if result)
    total = len(all_tests)
    
    print(f"\nTests passed: {passed}/{total}")
    
    for name, result in all_tests:
        status = "✓ PASS" if result else "❌ FAIL"
        print(f"  {status:8s} {name}")
    
    print("\n" + "="*70)
    
    # Overall status
    critical_tests = ['openmm', 'pymbar', 'custom_forces']
    critical_passed = all(results[test] for test in critical_tests)
    
    if critical_passed:
        print("✓ ENVIRONMENT READY FOR PRODUCTION")
        print("\nNext steps:")
        print("  1. Generate umbrella windows: python generate_umbrella_windows.py")
        print("  2. Submit SLURM job: sbatch submit_umbrella_production.sh")
        print("  3. Analyze results: python analyze_umbrella_mbar_CORRECTED.py")
        return 0
    else:
        print("❌ ENVIRONMENT NOT READY")
        print("\nFix the failed tests before running production simulations")
        print("See YOLTLA_DEPLOYMENT_GUIDE.md for troubleshooting")
        return 1

def main():
    print("="*70)
    print("WNK UMBRELLA SAMPLING ENVIRONMENT VERIFICATION")
    print("="*70)
    print(f"Python version: {sys.version}")
    print(f"Python executable: {sys.executable}")
    
    # Run all tests
    results = {
        'openmm': test_openmm(),
        'pymbar': test_pymbar(),
        'plumed': test_plumed(),
        'analysis': test_analysis_tools(),
        'custom_forces': test_custom_forces(),
        'filesystem': test_file_access(),
    }
    
    # Print summary
    exit_code = print_summary(results)
    sys.exit(exit_code)

if __name__ == '__main__':
    main()
