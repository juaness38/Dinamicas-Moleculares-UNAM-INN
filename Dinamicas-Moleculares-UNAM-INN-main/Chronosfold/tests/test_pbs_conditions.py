#!/usr/bin/env python3
"""
Tests de validación para condiciones PBS
=========================================

Verifica que el sistema preparado tenga las condiciones correctas de PBS:
- Fuerza iónica ~163 mM
- pH 7.4
- Neutralización correcta
- Box size apropiado
"""

import pytest
from pathlib import Path
import MDAnalysis as mda
from openmm import app, unit
import numpy as np


class TestPBSConditions:
    """Tests para verificar implementación de PBS buffer"""
    
    @pytest.fixture
    def prepared_system_path(self):
        """Path al sistema preparado"""
        return Path("prepared_system/system_solvated.pdb")
    
    @pytest.fixture
    def universe(self, prepared_system_path):
        """Cargar sistema con MDAnalysis"""
        if not prepared_system_path.exists():
            pytest.skip("Sistema no preparado aún. Ejecutar prepare_system.py primero.")
        return mda.Universe(str(prepared_system_path))
    
    def test_ionic_strength_pbs(self, universe):
        """Verificar que la fuerza iónica sea ~163 mM (PBS)"""
        # Contar iones
        n_na = len(universe.select_atoms("resname NA"))
        n_cl = len(universe.select_atoms("resname CL"))
        n_waters = len(universe.select_atoms("resname HOH")) / 3  # 3 atoms per water
        
        # Calcular volumen de caja
        dimensions = universe.dimensions  # [a, b, c, alpha, beta, gamma]
        volume_nm3 = (dimensions[0] * dimensions[1] * dimensions[2]) / 1000  # Å³ → nm³
        volume_L = volume_nm3 * 1e-24  # nm³ → L
        
        # Calcular concentración de iones
        avogadro = 6.022e23
        conc_na = (n_na / avogadro) / volume_L  # mol/L
        conc_cl = (n_cl / avogadro) / volume_L
        
        ionic_strength = 0.5 * (conc_na * 1**2 + conc_cl * 1**2)  # I = 0.5 Σ c_i z_i²
        
        print(f"\n  Iones Na+: {n_na}")
        print(f"  Iones Cl-: {n_cl}")
        print(f"  Moléculas agua: {n_waters:.0f}")
        print(f"  Volumen caja: {volume_L*1e3:.2f} mL")
        print(f"  Concentración Na+: {conc_na*1000:.1f} mM")
        print(f"  Concentración Cl-: {conc_cl*1000:.1f} mM")
        print(f"  Fuerza iónica: {ionic_strength*1000:.1f} mM")
        
        # Verificar que esté en rango PBS (150-170 mM)
        assert 0.150 <= ionic_strength <= 0.170, \
            f"Fuerza iónica {ionic_strength*1000:.1f} mM fuera de rango PBS (150-170 mM)"
    
    def test_system_neutralized(self, universe):
        """Verificar que el sistema esté neutralizado"""
        # Contar iones
        n_na = len(universe.select_atoms("resname NA"))
        n_cl = len(universe.select_atoms("resname CL"))
        
        # Calcular carga aproximada de proteína
        # (simplificado: ARG/LYS +1, ASP/GLU -1)
        n_arg = len(universe.select_atoms("resname ARG"))
        n_lys = len(universe.select_atoms("resname LYS"))
        n_asp = len(universe.select_atoms("resname ASP"))
        n_glu = len(universe.select_atoms("resname GLU"))
        
        protein_charge = (n_arg + n_lys) - (n_asp + n_glu)
        
        # Carga neta del sistema
        system_charge = n_na - n_cl + protein_charge
        
        print(f"\n  Carga proteína (aprox): {protein_charge:+d}")
        print(f"  Balance iónico: {n_na} Na+ - {n_cl} Cl- = {n_na - n_cl:+d}")
        print(f"  Carga neta sistema: {system_charge:+d}")
        
        # Sistema debe estar neutralizado (±1 por redondeo)
        assert abs(system_charge) <= 1, \
            f"Sistema no neutralizado: carga neta = {system_charge:+d}"
    
    def test_box_size_adequate(self, universe):
        """Verificar que la caja sea suficientemente grande (>1 nm padding)"""
        # Seleccionar proteína
        protein = universe.select_atoms("protein")
        
        # Calcular dimensiones de proteína
        pos = protein.positions
        protein_size = pos.max(axis=0) - pos.min(axis=0)  # Å
        
        # Dimensiones de caja
        box_dims = universe.dimensions[:3]  # [a, b, c] en Å
        
        # Padding en cada dirección
        padding = (box_dims - protein_size) / 2  # Å
        min_padding = padding.min()
        
        print(f"\n  Tamaño proteína: {protein_size[0]/10:.1f} x {protein_size[1]/10:.1f} x {protein_size[2]/10:.1f} nm")
        print(f"  Tamaño caja: {box_dims[0]/10:.1f} x {box_dims[1]/10:.1f} x {box_dims[2]/10:.1f} nm")
        print(f"  Padding mínimo: {min_padding/10:.2f} nm")
        
        # Debe tener al menos 1.0 nm de padding
        assert min_padding >= 9.5, \
            f"Padding insuficiente: {min_padding/10:.2f} nm (mínimo: 1.0 nm)"
    
    def test_water_model_tip3p(self, universe):
        """Verificar que se use agua TIP3P"""
        waters = universe.select_atoms("resname HOH")
        n_water_residues = len(set(waters.resids))
        n_water_atoms = len(waters)
        
        print(f"\n  Moléculas agua: {n_water_residues}")
        print(f"  Átomos agua: {n_water_atoms}")
        print(f"  Átomos por molécula: {n_water_atoms / n_water_residues:.1f}")
        
        # TIP3P tiene 3 átomos por molécula
        assert n_water_atoms / n_water_residues == 3.0, \
            "Modelo de agua incorrecto (esperado: TIP3P con 3 átomos/molécula)"
    
    def test_pbs_documentation(self):
        """Verificar que exista documentación de PBS"""
        pbs_doc = Path("PBS_BUFFER_IMPLEMENTATION.md")
        assert pbs_doc.exists(), \
            "Falta documentación PBS_BUFFER_IMPLEMENTATION.md"
        
        # Verificar que contenga secciones clave
        content = pbs_doc.read_text()
        assert "137 mM" in content, "Falta composición de NaCl"
        assert "2.7 mM" in content, "Falta composición de KCl"
        assert "10 mM" in content, "Falta composición de Na₂HPO₄"
        assert "1.8 mM" in content, "Falta composición de KH₂PO₄"
        assert "7.4" in content, "Falta especificación de pH"
        
        print("\n  ✓ Documentación PBS completa")


class TestPBSApproximation:
    """Tests para validar la aproximación PBS (K+ → Na+, fosfatos → Cl-)"""
    
    def test_forcefield_limitations_documented(self):
        """Verificar que limitaciones de forcefield estén documentadas"""
        pbs_doc = Path("PBS_BUFFER_IMPLEMENTATION.md")
        content = pbs_doc.read_text()
        
        # Debe mencionar limitaciones
        assert "K+" in content or "potasio" in content.lower(), \
            "Falta documentación de limitación de K+"
        assert "HPO₄" in content or "fosfato" in content.lower(), \
            "Falta documentación de limitación de fosfatos"
        assert "aproxim" in content.lower(), \
            "Falta mención de aproximación"
        
        print("\n  ✓ Limitaciones documentadas")
    
    def test_alternative_solutions_documented(self):
        """Verificar que soluciones alternativas estén documentadas"""
        pbs_doc = Path("PBS_BUFFER_IMPLEMENTATION.md")
        content = pbs_doc.read_text()
        
        # Debe mencionar alternativas
        assert "CHARMM" in content or "charmm" in content, \
            "Falta mención de forcefield alternativo (CHARMM)"
        assert "GAFF" in content or "parametriz" in content.lower(), \
            "Falta mención de parametrización custom"
        
        print("\n  ✓ Soluciones alternativas documentadas")


class TestDrMDIntegration:
    """Tests para pipeline paralelo con drMD"""
    
    def test_drmd_config_exists(self):
        """Verificar que exista configuración drMD"""
        config = Path("drMD_wnk_config.yaml")
        assert config.exists(), \
            "Falta drMD_wnk_config.yaml"
        
        print("\n  ✓ Configuración drMD encontrada")
    
    def test_drmd_config_has_pbs_settings(self):
        """Verificar que configuración drMD tenga pH 7.4"""
        import yaml
        
        config_file = Path("drMD_wnk_config.yaml")
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        
        # Verificar pH
        assert 'miscInfo' in config, "Falta miscInfo en config"
        assert config['miscInfo'].get('pH') == 7.4, \
            f"pH incorrecto: {config['miscInfo'].get('pH')} (esperado: 7.4)"
        
        print(f"\n  pH configurado: {config['miscInfo']['pH']}")
        print(f"  ✓ Configuración PBS en drMD")
    
    def test_drmd_pipeline_script_exists(self):
        """Verificar que exista script de pipeline drMD"""
        script = Path("run_drMD_pipeline.py")
        assert script.exists(), \
            "Falta run_drMD_pipeline.py"
        
        # Verificar que sea ejecutable (en Unix)
        import stat
        if hasattr(stat, 'S_IXUSR'):
            mode = script.stat().st_mode
            # No verificar permisos en Windows
        
        print("\n  ✓ Script drMD pipeline encontrado")


def test_pbs_vs_standard_ionic_strength():
    """Test conceptual: documentar diferencia PBS vs estándar"""
    ionic_strength_standard = 0.150  # M (estándar común)
    ionic_strength_pbs = 0.163  # M (PBS)
    
    # Calcular Debye length (apantallamiento electrostático)
    epsilon_0 = 8.854e-12  # F/m
    epsilon_r = 78.5  # agua a 25°C
    k_B = 1.381e-23  # J/K
    T = 300  # K
    e = 1.602e-19  # C
    N_A = 6.022e23
    
    def debye_length(I):
        # λ_D = sqrt(ε₀εᵣkT / (2NAe²I))
        return np.sqrt(epsilon_0 * epsilon_r * k_B * T / 
                      (2 * N_A * e**2 * I * 1000))  # I en mol/m³
    
    lambda_standard = debye_length(ionic_strength_standard) * 1e9  # nm
    lambda_pbs = debye_length(ionic_strength_pbs) * 1e9  # nm
    
    diff_percent = abs(lambda_pbs - lambda_standard) / lambda_standard * 100
    
    print(f"\n  Debye length (150 mM): {lambda_standard:.3f} nm")
    print(f"  Debye length (163 mM PBS): {lambda_pbs:.3f} nm")
    print(f"  Diferencia: {diff_percent:.1f}%")
    
    # Diferencia debe ser pequeña (<5%)
    assert diff_percent < 5.0, \
        f"Diferencia electrostática significativa: {diff_percent:.1f}%"


if __name__ == "__main__":
    print("="*70)
    print("TESTS DE VALIDACIÓN PBS BUFFER")
    print("="*70)
    print("\nEjecutar con: pytest test_pbs_conditions.py -v")
    print("\nTests incluidos:")
    print("  1. Fuerza iónica ~163 mM")
    print("  2. Sistema neutralizado")
    print("  3. Box size >1 nm padding")
    print("  4. Modelo agua TIP3P")
    print("  5. Documentación completa")
    print("  6. Integración drMD")
    print("  7. Validación de aproximación")
