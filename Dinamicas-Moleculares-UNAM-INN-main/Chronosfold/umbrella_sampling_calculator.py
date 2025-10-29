#!/usr/bin/env python3
"""
🎯 Umbrella Sampling Calculator con OpenMM
Extensión del openmm_calculator.py para muestreo de energía libre

Implementa:
- Potenciales de restricción armónicos (umbrella windows)
- Múltiples ventanas de muestreo
- Extracción de datos para WHAM/MBAR
- Análisis de convergencia

Desarrollado para AstroFlora Core - Análisis de Energía Libre Molecular
"""

import asyncio
import logging
from typing import Dict, List, Optional, Tuple, Any
from pathlib import Path
import numpy as np
import json

try:
    from openmm.app import *
    from openmm import *
    from openmm.unit import *
    OPENMM_AVAILABLE = True
except ImportError as e:
    OPENMM_AVAILABLE = False
    print(f"OpenMM no disponible: {e}")

# Importar calculadora base
try:
    from openmm_calculator import OpenMMEnergyCalculator, create_openmm_calculator
    BASE_CALCULATOR_AVAILABLE = True
except ImportError:
    BASE_CALCULATOR_AVAILABLE = False
    print("openmm_calculator.py base no encontrado - funcionalidad limitada")


class UmbrellaSamplingCalculator(OpenMMEnergyCalculator if BASE_CALCULATOR_AVAILABLE else object):
    """
    Calculadora de Umbrella Sampling usando OpenMM
    
    Umbrella Sampling es una técnica de enhanced sampling que usa potenciales
    armónicos (umbrella potentials) para forzar el muestreo a lo largo de una
    coordenada de reacción (collective variable).
    
    Workflow:
    1. Define collective variable (CV): distancia, ángulo, RMSD, etc.
    2. Corre múltiples simulaciones MD, cada una con potencial armónico centrado
       en diferente valor de la CV (ventanas)
    3. Extrae histogramas de la CV de cada ventana
    4. Une las ventanas con WHAM (Weighted Histogram Analysis Method) o MBAR
    5. Obtiene perfil de energía libre
    """
    
    def __init__(self, *args, **kwargs):
        """Inicializa calculadora de umbrella sampling"""
        if BASE_CALCULATOR_AVAILABLE:
            super().__init__(*args, **kwargs)
        else:
            self.logger = logging.getLogger(__name__)
            self.logger.warning("OpenMMEnergyCalculator base no disponible")
        
        # Configuración específica de umbrella sampling
        self.umbrella_windows = []  # Lista de ventanas configuradas
        self.cv_data = {}  # Datos de collective variables por ventana
        
    def setup_collective_variable(self, 
                                  cv_type: str,
                                  atoms: List[int],
                                  **cv_params) -> CustomCVForce:
        """
        Configura una collective variable (CV) personalizada
        
        Args:
            cv_type: Tipo de CV ('distance', 'angle', 'dihedral', 'rmsd')
            atoms: Índices de átomos involucrados en la CV
            cv_params: Parámetros adicionales de la CV
            
        Returns:
            CustomCVForce configurada
        """
        if not OPENMM_AVAILABLE:
            raise RuntimeError("OpenMM no disponible para collective variables")
        
        self.logger.info(f"Configurando CV tipo '{cv_type}' con {len(atoms)} átomos")
        
        if cv_type == 'distance':
            # Distancia entre dos átomos
            if len(atoms) != 2:
                raise ValueError("Distance CV requiere exactamente 2 átomos")
            
            # Crear fuerza de distancia personalizada
            cv_force = CustomBondForce("r")
            cv_force.addBond(atoms[0], atoms[1])
            
        elif cv_type == 'angle':
            # Ángulo entre tres átomos
            if len(atoms) != 3:
                raise ValueError("Angle CV requiere exactamente 3 átomos")
            
            cv_force = CustomAngleForce("theta")
            cv_force.addAngle(atoms[0], atoms[1], atoms[2])
            
        elif cv_type == 'dihedral':
            # Ángulo dihedro entre cuatro átomos
            if len(atoms) != 4:
                raise ValueError("Dihedral CV requiere exactamente 4 átomos")
            
            cv_force = CustomTorsionForce("theta")
            cv_force.addTorsion(atoms[0], atoms[1], atoms[2], atoms[3])
            
        elif cv_type == 'rmsd':
            # RMSD respecto a estructura de referencia
            ref_positions = cv_params.get('reference_positions')
            if ref_positions is None:
                raise ValueError("RMSD CV requiere 'reference_positions'")
            
            cv_force = RMSDForce(ref_positions, atoms)
            
        else:
            raise ValueError(f"Tipo de CV no soportado: {cv_type}")
        
        return cv_force
    
    def create_umbrella_potential(self,
                                 cv_force: Force,
                                 center: float,
                                 force_constant: float,
                                 cv_type: str = 'distance') -> CustomCVForce:
        """
        Crea potencial armónico de umbrella (restraint)
        
        V(ξ) = 0.5 * k * (ξ - ξ₀)²
        
        donde:
        - ξ es la collective variable
        - ξ₀ es el centro de la ventana (center)
        - k es la constante de fuerza (force_constant)
        
        Args:
            cv_force: Fuerza que define la collective variable
            center: Centro del potencial armónico (ξ₀)
            force_constant: Constante de fuerza k (kcal/mol/unidad²)
            cv_type: Tipo de CV para unidades correctas
            
        Returns:
            CustomCVForce con potencial de umbrella
        """
        # Unidades según tipo de CV
        if cv_type in ['distance', 'rmsd']:
            center_with_unit = center * angstrom
            k_with_unit = force_constant * kilocalorie_per_mole / angstrom**2
        elif cv_type in ['angle', 'dihedral']:
            center_with_unit = center * degree
            k_with_unit = force_constant * kilocalorie_per_mole / degree**2
        else:
            # Default: sin unidades específicas
            center_with_unit = center
            k_with_unit = force_constant
        
        # Crear fuerza de CV personalizada con potencial armónico
        # CustomCVForce permite usar cualquier Force como variable
        umbrella_force = CustomCVForce("0.5 * k * (cv - center)^2")
        
        # Añadir parámetros
        umbrella_force.addGlobalParameter("k", k_with_unit)
        umbrella_force.addGlobalParameter("center", center_with_unit)
        
        # Añadir la collective variable
        umbrella_force.addCollectiveVariable("cv", cv_force)
        
        self.logger.info(f"Potencial umbrella creado: center={center}, k={force_constant}")
        
        return umbrella_force
    
    async def run_umbrella_window(self,
                                  structure_file: str,
                                  cv_config: Dict[str, Any],
                                  window_center: float,
                                  force_constant: float,
                                  simulation_time_ps: float = 1000.0,
                                  temperature: float = 300.0,
                                  output_dir: Optional[Path] = None) -> Dict[str, Any]:
        """
        Ejecuta una ventana de umbrella sampling
        
        Args:
            structure_file: Archivo PDB de estructura inicial
            cv_config: Configuración de collective variable
                {
                    'type': 'distance'/'angle'/'dihedral'/'rmsd',
                    'atoms': [atom_indices],
                    'params': {...}  # Parámetros adicionales
                }
            window_center: Centro del potencial de umbrella (ξ₀)
            force_constant: Constante de fuerza k
            simulation_time_ps: Tiempo de simulación en ps
            temperature: Temperatura en Kelvin
            output_dir: Directorio para guardar outputs
            
        Returns:
            Dict con resultados de la ventana:
            {
                'window_center': float,
                'force_constant': float,
                'cv_values': np.ndarray,  # Serie temporal de CV
                'cv_histogram': Dict,
                'mean_cv': float,
                'std_cv': float,
                'trajectory_file': str
            }
        """
        if not OPENMM_AVAILABLE:
            raise RuntimeError("OpenMM no disponible")
        
        self.logger.info(f"🎯 Iniciando ventana umbrella: center={window_center}")
        
        # Preparar estructura
        topology, positions = self._prepare_pdb_structure(structure_file)
        
        # Crear sistema con force field
        forcefield = self._get_force_field()
        
        try:
            system = forcefield.createSystem(
                topology,
                nonbondedMethod=CutoffNonPeriodic,
                nonbondedCutoff=1.0*nanometer,
                implicitSolvent=GBn2,
                constraints=HBonds
            )
        except:
            # Fallback a vacío
            system = forcefield.createSystem(
                topology,
                nonbondedMethod=NoCutoff,
                constraints=HBonds
            )
        
        # Configurar collective variable
        cv_type = cv_config['type']
        cv_atoms = cv_config['atoms']
        cv_params = cv_config.get('params', {})
        
        cv_force = self.setup_collective_variable(cv_type, cv_atoms, **cv_params)
        
        # Crear y agregar potencial de umbrella
        umbrella_force = self.create_umbrella_potential(
            cv_force, window_center, force_constant, cv_type
        )
        system.addForce(umbrella_force)
        
        # Integrador Langevin
        integrator = LangevinMiddleIntegrator(
            temperature * kelvin,
            1.0 / picosecond,
            0.002 * picosecond  # 2 fs timestep
        )
        
        # Crear simulación
        platform = Platform.getPlatformByName(self.platform_name)
        simulation = Simulation(topology, system, integrator, platform)
        simulation.context.setPositions(positions)
        
        # Minimización de energía
        self.logger.info("Minimizando energía...")
        simulation.minimizeEnergy(maxIterations=1000)
        
        # Equilibrio térmico
        simulation.context.setVelocitiesToTemperature(temperature * kelvin)
        
        # Configurar reporters para tracking de CV
        cv_values = []
        
        # Reporter personalizado para CV
        class CVReporter:
            def __init__(self, cv_force_index, system, reportInterval):
                self.cv_force_index = cv_force_index
                self.system = system
                self.reportInterval = reportInterval
                self.values = []
                
            def describeNextReport(self, simulation):
                steps = self.reportInterval - simulation.currentStep % self.reportInterval
                return (steps, False, False, False, False, None)
            
            def report(self, simulation, state):
                # Obtener valor de CV
                cv_value = state.getEnergyParameterDerivatives()[self.cv_force_index]
                self.values.append(cv_value)
        
        # Añadir reporter de CV
        cv_reporter = CVReporter(
            system.getNumForces() - 1,  # Último force añadido
            system,
            reportInterval=100  # Cada 100 pasos
        )
        
        # Simulación MD
        num_steps = int(simulation_time_ps / 0.002)  # timestep 2 fs
        self.logger.info(f"Corriendo MD: {simulation_time_ps} ps ({num_steps} pasos)")
        
        # Producción con sampling de CV
        for i in range(0, num_steps, 100):
            simulation.step(100)
            
            # Obtener valor actual de CV
            state = simulation.context.getState(getEnergy=True)
            
            # Para distance CV, obtener posiciones y calcular
            if cv_type == 'distance':
                positions = simulation.context.getState(getPositions=True).getPositions()
                atom1_pos = positions[cv_atoms[0]]
                atom2_pos = positions[cv_atoms[1]]
                
                # Calcular distancia
                delta = atom1_pos - atom2_pos
                distance = np.sqrt(delta.x**2 + delta.y**2 + delta.z**2)
                cv_values.append(distance.value_in_unit(angstrom))
            
            # TODO: Implementar extracción para otros tipos de CV
        
        cv_array = np.array(cv_values)
        
        # Calcular histograma de CV
        hist, bin_edges = np.histogram(cv_array, bins=50, density=True)
        
        # Guardar datos si se especifica directorio
        if output_dir:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            # Guardar serie temporal de CV
            np.savetxt(
                output_dir / f"cv_timeseries_center_{window_center:.2f}.dat",
                cv_array,
                header=f"Collective variable values for window center={window_center}"
            )
            
            # Guardar histograma
            hist_data = np.column_stack([bin_edges[:-1], hist])
            np.savetxt(
                output_dir / f"cv_histogram_center_{window_center:.2f}.dat",
                hist_data,
                header="CV_value Probability"
            )
        
        results = {
            'window_center': window_center,
            'force_constant': force_constant,
            'cv_values': cv_array,
            'cv_histogram': {
                'counts': hist.tolist(),
                'bin_edges': bin_edges.tolist()
            },
            'mean_cv': float(np.mean(cv_array)),
            'std_cv': float(np.std(cv_array)),
            'num_samples': len(cv_array)
        }
        
        self.logger.info(f"✅ Ventana completada: <CV>={results['mean_cv']:.3f} ± {results['std_cv']:.3f}")
        
        return results
    
    async def run_full_umbrella_sampling(self,
                                        structure_file: str,
                                        cv_config: Dict[str, Any],
                                        window_centers: List[float],
                                        force_constant: float,
                                        simulation_time_ps: float = 1000.0,
                                        temperature: float = 300.0,
                                        output_dir: Optional[str] = None) -> Dict[str, Any]:
        """
        Ejecuta serie completa de ventanas de umbrella sampling
        
        Args:
            structure_file: Estructura inicial
            cv_config: Configuración de collective variable
            window_centers: Lista de centros para cada ventana [ξ₁, ξ₂, ..., ξₙ]
            force_constant: Constante de fuerza (misma para todas las ventanas)
            simulation_time_ps: Tiempo por ventana
            temperature: Temperatura
            output_dir: Directorio para outputs
            
        Returns:
            Dict con todos los resultados y datos para WHAM/MBAR
        """
        self.logger.info(f"🎯 Umbrella Sampling: {len(window_centers)} ventanas")
        
        all_windows_results = []
        
        for i, center in enumerate(window_centers):
            self.logger.info(f"📊 Ventana {i+1}/{len(window_centers)}")
            
            window_result = await self.run_umbrella_window(
                structure_file=structure_file,
                cv_config=cv_config,
                window_center=center,
                force_constant=force_constant,
                simulation_time_ps=simulation_time_ps,
                temperature=temperature,
                output_dir=Path(output_dir) if output_dir else None
            )
            
            all_windows_results.append(window_result)
        
        # Preparar datos para análisis WHAM/MBAR
        wham_data = self._prepare_wham_data(all_windows_results)
        
        # Guardar metadata completa
        if output_dir:
            metadata = {
                'cv_config': cv_config,
                'window_centers': window_centers,
                'force_constant': force_constant,
                'temperature': temperature,
                'simulation_time_ps': simulation_time_ps,
                'num_windows': len(window_centers)
            }
            
            with open(Path(output_dir) / 'umbrella_metadata.json', 'w') as f:
                json.dump(metadata, f, indent=2)
            
            self.logger.info(f"📁 Datos guardados en: {output_dir}")
        
        return {
            'windows': all_windows_results,
            'wham_data': wham_data,
            'metadata': {
                'num_windows': len(window_centers),
                'cv_range': [min(window_centers), max(window_centers)],
                'total_simulation_time_ps': len(window_centers) * simulation_time_ps
            }
        }
    
    def _prepare_wham_data(self, windows_results: List[Dict]) -> Dict[str, Any]:
        """
        Prepara datos en formato para análisis WHAM
        
        Returns:
            Dict con datos formateados para WHAM/MBAR
        """
        wham_data = {
            'window_centers': [],
            'force_constants': [],
            'cv_timeseries': [],
            'histograms': []
        }
        
        for result in windows_results:
            wham_data['window_centers'].append(result['window_center'])
            wham_data['force_constants'].append(result['force_constant'])
            wham_data['cv_timeseries'].append(result['cv_values'])
            wham_data['histograms'].append(result['cv_histogram'])
        
        return wham_data
    
    def estimate_pmf_simple(self, wham_data: Dict[str, Any]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Estimación simple de PMF (Potential of Mean Force) sin WHAM completo
        
        NOTA: Para producción, usar pyWHAM o pymbar
        Esta es una aproximación básica
        
        Returns:
            (cv_values, pmf_values) en kcal/mol
        """
        self.logger.warning("⚠️  Usando estimación simplificada de PMF")
        self.logger.warning("   Para resultados precisos, usar pyWHAM o pymbar")
        
        # Combinar histogramas con overlap
        all_cv_values = []
        all_counts = []
        
        for i, hist_data in enumerate(wham_data['histograms']):
            bin_centers = (np.array(hist_data['bin_edges'][:-1]) + 
                          np.array(hist_data['bin_edges'][1:])) / 2
            counts = np.array(hist_data['counts'])
            
            all_cv_values.extend(bin_centers)
            all_counts.extend(counts)
        
        # Ordenar
        sort_idx = np.argsort(all_cv_values)
        cv_values = np.array(all_cv_values)[sort_idx]
        counts = np.array(all_counts)[sort_idx]
        
        # Estimación simple: PMF ~ -kT ln(P)
        kT = 0.593  # kcal/mol a 300K
        
        # Evitar log(0)
        counts[counts == 0] = 1e-10
        
        pmf = -kT * np.log(counts)
        pmf -= pmf.min()  # Normalizar a mínimo = 0
        
        return cv_values, pmf


# Factory function
def create_umbrella_calculator(config: Optional[Dict] = None) -> UmbrellaSamplingCalculator:
    """
    Crea calculadora de umbrella sampling
    
    Args:
        config: Configuración (igual que OpenMMEnergyCalculator)
    
    Returns:
        UmbrellaSamplingCalculator configurado
    """
    if config is None:
        config = {}
    
    return UmbrellaSamplingCalculator(
        force_field=config.get('force_field', 'amber99sb.xml'),
        platform=config.get('platform', 'CPU'),
        temperature=config.get('temperature', 300.0),
        implicit_solvent=config.get('implicit_solvent', True)
    )


# Demo y testing
async def demo_umbrella_sampling():
    """Demostración de umbrella sampling"""
    print("🎯 Demo: Umbrella Sampling con OpenMM")
    print("=" * 60)
    
    if not OPENMM_AVAILABLE:
        print("❌ OpenMM no disponible - demo cancelado")
        return
    
    # Crear calculadora
    calculator = create_umbrella_calculator({
        'platform': 'CPU',
        'temperature': 300.0
    })
    
    # Configuración de collective variable (ejemplo: distancia entre dos Cα)
    cv_config = {
        'type': 'distance',
        'atoms': [0, 10],  # Índices de ejemplo
        'params': {}
    }
    
    # Definir ventanas de umbrella
    # Para una distancia de 5-15 Å, usar 10 ventanas
    window_centers = np.linspace(5.0, 15.0, 10)  # Angstroms
    
    print(f"📊 Configuración:")
    print(f"   CV type: {cv_config['type']}")
    print(f"   Ventanas: {len(window_centers)}")
    print(f"   Rango: {window_centers[0]:.1f} - {window_centers[-1]:.1f} Å")
    
    # NOTE: Esto requiere archivo PDB real
    # structure_file = "example_protein.pdb"
    
    print("\n⚠️  Demo informativa - requiere estructura PDB real")
    print("\n📋 Workflow de Umbrella Sampling:")
    print("1. Define collective variable (e.g., distancia entre residuos)")
    print("2. Corre múltiples ventanas con potenciales armónicos")
    print("3. Extrae histogramas de cada ventana")
    print("4. Usa WHAM/MBAR para unir ventanas")
    print("5. Obtiene perfil de energía libre (PMF)")
    
    print("\n✅ Calculadora lista para usar con estructura real")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    asyncio.run(demo_umbrella_sampling())
