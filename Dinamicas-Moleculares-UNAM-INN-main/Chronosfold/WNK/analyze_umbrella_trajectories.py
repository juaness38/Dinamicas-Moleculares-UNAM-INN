#!/usr/bin/env python3
"""
Análisis de trayectorias de umbrella sampling.
Qué SÍ puedes analizar (estructuras, contactos) y qué NO (cinética).
"""

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import argparse

# ============================================================================
# ✅ ANÁLISIS VÁLIDOS (Propiedades Termodinámica/Estructurales)
# ============================================================================

def extract_representative_structures(window_id, n_clusters=5, output_prefix='repr'):
    """
    Extrae estructuras representativas de una ventana mediante clustering.
    
    ✅ VÁLIDO: Las conformaciones SON físicamente correctas.
    El bias solo afecta CUÁNTO tiempo pasas en esta ventana,
    NO QUÉ estructuras adopta el sistema cuando CV = window_center.
    
    Parámetros:
    -----------
    window_id : int
        ID de la ventana (0-19)
    n_clusters : int
        Número de clusters (conformaciones representativas)
    """
    
    print(f"\n{'='*60}")
    print(f"EXTRAYENDO ESTRUCTURAS REPRESENTATIVAS - VENTANA {window_id}")
    print(f"{'='*60}")
    
    # Cargar trayectoria
    traj_file = f'window_{window_id}_prod.dcd'
    top_file = 'system.pdb'
    
    print(f"Cargando {traj_file}...")
    traj = md.load(traj_file, top=top_file)
    print(f"✅ {len(traj)} frames cargados ({traj.time[-1]/1000:.1f} ns)")
    
    # Seleccionar solo proteína (más rápido para clustering)
    protein = traj.topology.select('protein')
    traj_protein = traj.atom_slice(protein)
    
    # Alinear (quitar rotaciones/traslaciones)
    traj_protein.superpose(traj_protein[0])
    
    # Clustering basado en RMSD
    print(f"\nRealizando clustering (k-means, k={n_clusters})...")
    
    # Matriz de distancias RMSD (reducida para velocidad)
    indices = np.random.choice(len(traj_protein), 
                               min(1000, len(traj_protein)), 
                               replace=False)
    
    rmsd_matrix = np.zeros((len(indices), len(indices)))
    for i, idx_i in enumerate(indices):
        rmsd_matrix[i] = md.rmsd(traj_protein[indices], traj_protein[idx_i])
    
    # K-means en espacio de RMSD
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    labels = kmeans.fit_predict(rmsd_matrix)
    
    # Extraer representante de cada cluster (más cercano al centroide)
    representatives = []
    
    for cluster_id in range(n_clusters):
        cluster_frames = indices[labels == cluster_id]
        
        if len(cluster_frames) == 0:
            continue
        
        # Calcular centroide
        centroid = traj_protein[cluster_frames].xyz.mean(axis=0)
        
        # Frame más cercano al centroide
        distances = np.sqrt(((traj_protein[cluster_frames].xyz - centroid)**2).sum(axis=(1,2)))
        repr_idx = cluster_frames[distances.argmin()]
        
        representatives.append(repr_idx)
        
        # Guardar estructura
        output_file = f'{output_prefix}_window{window_id}_cluster{cluster_id}.pdb'
        traj[repr_idx].save_pdb(output_file)
        
        print(f"   Cluster {cluster_id}: {len(cluster_frames)} frames "
              f"({100*len(cluster_frames)/len(indices):.1f}%) → {output_file}")
    
    print(f"\n✅ {len(representatives)} estructuras representativas guardadas")
    return representatives


def analyze_contacts(window_id, contact_pairs, cutoff=0.45):
    """
    Analiza contactos (distancias < cutoff) en una ventana.
    
    ✅ VÁLIDO: Los contactos SON físicamente correctos.
    El bias no afecta QUÉ residuos interactúan cuando CV = window_center.
    
    Parámetros:
    -----------
    window_id : int
        ID de la ventana
    contact_pairs : list of tuples
        Pares de residuos a analizar, ej. [(5, 120), (8, 115)]
    cutoff : float
        Distancia de corte para considerar contacto (nm)
    """
    
    print(f"\n{'='*60}")
    print(f"ANALIZANDO CONTACTOS - VENTANA {window_id}")
    print(f"{'='*60}")
    
    # Cargar trayectoria
    traj = md.load(f'window_{window_id}_prod.dcd', top='system.pdb')
    
    # Calcular distancias
    print(f"\nCalculando distancias para {len(contact_pairs)} pares...")
    
    distances, residue_pairs = md.compute_contacts(
        traj, 
        contacts=contact_pairs,
        scheme='closest-heavy'
    )
    
    # Analizar cada par
    for i, (res1, res2) in enumerate(contact_pairs):
        dist = distances[:, i]
        
        # Fracción de tiempo en contacto
        fraction_contact = (dist < cutoff).sum() / len(dist)
        
        # Estadísticas
        mean_dist = dist.mean()
        std_dist = dist.std()
        
        print(f"\n   Par {res1}-{res2}:")
        print(f"      Distancia media: {mean_dist:.3f} ± {std_dist:.3f} nm")
        print(f"      Fracción en contacto (<{cutoff} nm): {fraction_contact:.2%}")
        print(f"      Mín/Máx: {dist.min():.3f} / {dist.max():.3f} nm")
        
        # Clasificar
        if fraction_contact > 0.8:
            print(f"      ✅ Contacto ESTABLE")
        elif fraction_contact > 0.3:
            print(f"      ⚠️  Contacto INTERMITENTE")
        else:
            print(f"      ❌ SIN contacto")
    
    return distances


def analyze_secondary_structure(window_id):
    """
    Analiza estructura secundaria (hélices, sheets) en una ventana.
    
    ✅ VÁLIDO: La estructura secundaria ES físicamente correcta.
    """
    
    print(f"\n{'='*60}")
    print(f"ANÁLISIS DE ESTRUCTURA SECUNDARIA - VENTANA {window_id}")
    print(f"{'='*60}")
    
    # Cargar trayectoria
    traj = md.load(f'window_{window_id}_prod.dcd', top='system.pdb')
    
    # Calcular DSSP (estructura secundaria)
    print("\nCalculando DSSP...")
    dssp = md.compute_dssp(traj, simplified=True)
    
    # H = alpha helix, E = beta sheet, C = coil
    ss_types = {'H': 'Alpha-helix', 'E': 'Beta-sheet', 'C': 'Coil'}
    
    for ss_code, ss_name in ss_types.items():
        fraction = (dssp == ss_code).sum() / dssp.size
        print(f"   {ss_name}: {fraction:.1%}")
    
    # Análisis por residuo (promedio temporal)
    ss_per_residue = []
    for res_idx in range(dssp.shape[1]):
        ss_profile = dssp[:, res_idx]
        most_common = max(set(ss_profile), key=list(ss_profile).count)
        ss_per_residue.append(most_common)
    
    # Identificar regiones helicoidales
    helix_regions = []
    in_helix = False
    start = None
    
    for i, ss in enumerate(ss_per_residue):
        if ss == 'H' and not in_helix:
            start = i
            in_helix = True
        elif ss != 'H' and in_helix:
            helix_regions.append((start, i-1))
            in_helix = False
    
    if in_helix:
        helix_regions.append((start, len(ss_per_residue)-1))
    
    print(f"\n   Regiones helicoidales identificadas:")
    for start, end in helix_regions:
        print(f"      Residuos {start}-{end} ({end-start+1} residuos)")
    
    return dssp


def calculate_rg_and_sasa(window_id):
    """
    Calcula radio de giro (Rg) y superficie accesible (SASA).
    
    ✅ VÁLIDO: Propiedades geométricas correctas.
    """
    
    print(f"\n{'='*60}")
    print(f"PROPIEDADES GEOMÉTRICAS - VENTANA {window_id}")
    print(f"{'='*60}")
    
    # Cargar trayectoria
    traj = md.load(f'window_{window_id}_prod.dcd', top='system.pdb')
    
    # Solo proteína
    protein = traj.topology.select('protein')
    traj_protein = traj.atom_slice(protein)
    
    # Radio de giro
    rg = md.compute_rg(traj_protein)
    print(f"\n   Radio de giro (Rg):")
    print(f"      Media: {rg.mean():.3f} ± {rg.std():.3f} nm")
    print(f"      Rango: [{rg.min():.3f}, {rg.max():.3f}] nm")
    
    # SASA (superficie accesible al solvente)
    sasa = md.shrake_rupley(traj_protein)
    total_sasa = sasa.sum(axis=1)
    
    print(f"\n   Superficie accesible (SASA):")
    print(f"      Media: {total_sasa.mean():.2f} ± {total_sasa.std():.2f} nm²")
    print(f"      Rango: [{total_sasa.min():.2f}, {total_sasa.max():.2f}] nm²")
    
    return rg, sasa


# ============================================================================
# ❌ ANÁLISIS INVÁLIDOS (Propiedades Cinéticas) - SOLO PARA DEMOSTRACIÓN
# ============================================================================

def measure_residence_time_INVALID(window_id, cv_min, cv_max):
    """
    ❌ INCORRECTO: Mide tiempo de residencia en región de CV.
    
    PROBLEMA: El bias armónico crea una trampa artificial.
    El tiempo medido NO refleja la cinética real del sistema.
    
    Este código se incluye SOLO para demostrar por qué es inválido.
    """
    
    print(f"\n{'='*60}")
    print(f"⚠️  ANÁLISIS INVÁLIDO (SOLO DEMOSTRACIÓN)")
    print(f"{'='*60}")
    
    # Cargar trayectoria
    traj = md.load(f'window_{window_id}_prod.dcd', top='system.pdb')
    
    # Calcular CV (asumiendo distancia simple)
    # En realidad, necesitarías calcular la misma CV que en umbrella
    cv_atoms = [10, 1250]  # Ejemplo
    cv_values = md.compute_distances(traj, [cv_atoms])[:, 0]
    
    # Frames dentro de la región
    in_region = (cv_values >= cv_min) & (cv_values <= cv_max)
    fraction_in_region = in_region.sum() / len(in_region)
    
    # "Tiempo de residencia" (SESGADO)
    residence_time = fraction_in_region * traj.time[-1]
    
    print(f"\n❌ Tiempo de residencia en CV=[{cv_min:.2f}, {cv_max:.2f}]:")
    print(f"   {residence_time:.2f} ps")
    print(f"\n⚠️  ADVERTENCIA: Este valor está ARTIFICIALMENTE EXTENDIDO")
    print(f"   por el bias armónico. NO representa cinética real.")
    print(f"   El sistema 'debería' salir más rápido, pero el resorte")
    print(f"   lo retiene artificialmente.")
    
    return residence_time


def count_transitions_INVALID(window_id, barrier_cv):
    """
    ❌ INCORRECTO: Cuenta transiciones a través de barrera.
    
    PROBLEMA: El bias facilita/dificulta cruces artificialmente.
    La frecuencia medida NO refleja la cinética real.
    """
    
    print(f"\n{'='*60}")
    print(f"⚠️  ANÁLISIS INVÁLIDO (SOLO DEMOSTRACIÓN)")
    print(f"{'='*60}")
    
    # Cargar trayectoria
    traj = md.load(f'window_{window_id}_prod.dcd', top='system.pdb')
    
    # Calcular CV
    cv_atoms = [10, 1250]
    cv_values = md.compute_distances(traj, [cv_atoms])[:, 0]
    
    # Detectar cruces
    crossings = 0
    for i in range(1, len(cv_values)):
        if (cv_values[i-1] < barrier_cv) and (cv_values[i] >= barrier_cv):
            crossings += 1
        elif (cv_values[i-1] >= barrier_cv) and (cv_values[i] < barrier_cv):
            crossings += 1
    
    # "Frecuencia" (SESGADA)
    frequency = crossings / (traj.time[-1] / 1000)  # cruces/ns
    
    print(f"\n❌ Frecuencia de cruces en CV={barrier_cv:.2f} nm:")
    print(f"   {crossings} cruces en {traj.time[-1]/1000:.1f} ns")
    print(f"   Frecuencia: {frequency:.2f} cruces/ns")
    print(f"\n⚠️  ADVERTENCIA: Esta frecuencia es ARTIFICIAL.")
    print(f"   En ventanas cercanas a la barrera, el bias EMPUJA")
    print(f"   hacia ella, creando cruces que no ocurrirían naturalmente.")
    
    return crossings, frequency


# ============================================================================
# COMPARACIÓN: Análisis Válido vs Inválido
# ============================================================================

def demonstrate_valid_vs_invalid(window_id=10):
    """
    Demuestra qué análisis son válidos y cuáles no.
    """
    
    print("\n" + "="*70)
    print("DEMOSTRACIÓN: ANÁLISIS VÁLIDOS vs INVÁLIDOS")
    print("="*70)
    
    print("\n✅ ANÁLISIS VÁLIDOS (Termodinámica/Estructura):")
    print("   " + "-"*60)
    
    # 1. Estructuras representativas
    print("\n1. Extrayendo estructuras representativas...")
    representatives = extract_representative_structures(window_id, n_clusters=3)
    print("   ✅ Resultado VÁLIDO: Conformaciones físicamente correctas")
    
    # 2. Contactos
    print("\n2. Analizando contactos clave...")
    contact_pairs = [(1255, 1268), (1250, 1260)]  # Ejemplo
    contacts = analyze_contacts(window_id, contact_pairs)
    print("   ✅ Resultado VÁLIDO: Interacciones físicamente correctas")
    
    # 3. Estructura secundaria
    print("\n3. Analizando estructura secundaria...")
    dssp = analyze_secondary_structure(window_id)
    print("   ✅ Resultado VÁLIDO: Elementos estructurales correctos")
    
    # 4. Propiedades geométricas
    print("\n4. Calculando Rg y SASA...")
    rg, sasa = calculate_rg_and_sasa(window_id)
    print("   ✅ Resultado VÁLIDO: Propiedades geométricas correctas")
    
    print("\n" + "-"*70)
    print("\n❌ ANÁLISIS INVÁLIDOS (Cinética):")
    print("   " + "-"*60)
    
    # 5. Tiempo de residencia (INVÁLIDO)
    print("\n5. Midiendo tiempo de residencia (INVÁLIDO)...")
    residence = measure_residence_time_INVALID(window_id, cv_min=2.9, cv_max=3.1)
    print("   ❌ Resultado SESGADO: NO representa cinética real")
    
    # 6. Frecuencia de transiciones (INVÁLIDO)
    print("\n6. Contando transiciones (INVÁLIDO)...")
    crossings, freq = count_transitions_INVALID(window_id, barrier_cv=3.0)
    print("   ❌ Resultado SESGADO: NO representa cinética real")
    
    print("\n" + "="*70)
    print("RESUMEN:")
    print("="*70)
    print("""
    ✅ Umbrella sampling (con MBAR) proporciona:
       • PMF correcto (termodinámica)
       • Estructuras representativas (conformaciones)
       • Contactos y propiedades locales (geometría)
    
    ❌ Umbrella sampling NO proporciona:
       • Tiempos de residencia (cinética sesgada)
       • Frecuencias de transición (sesgadas por bias)
       • Constantes de velocidad k_on, k_off
       • Tiempos de vida t₁/₂
    
    Para cinética, necesitas:
       • Weighted Ensemble (WESTPA)
       • Milestoning
       • MD ultra-largo sin bias (fuerza bruta)
    """)


# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Análisis de trayectorias de umbrella sampling',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:

  # Análisis válidos (estructuras, contactos)
  python analyze_umbrella_trajectories.py --window 10 --mode structures
  python analyze_umbrella_trajectories.py --window 10 --mode contacts
  python analyze_umbrella_trajectories.py --window 10 --mode secondary
  python analyze_umbrella_trajectories.py --window 10 --mode geometry
  
  # Demostración completa (válidos vs inválidos)
  python analyze_umbrella_trajectories.py --window 10 --mode demo
  
Nota: Los análisis cinéticos (tiempos, frecuencias) están marcados
      como INVÁLIDOS porque las trayectorias están sesgadas por el
      bias armónico. Solo se incluyen para demostración educativa.
        """
    )
    
    parser.add_argument('--window', type=int, default=10,
                       help='ID de la ventana a analizar (0-19)')
    parser.add_argument('--mode', choices=['structures', 'contacts', 'secondary', 
                                           'geometry', 'demo'],
                       default='demo',
                       help='Tipo de análisis')
    
    args = parser.parse_args()
    
    if args.mode == 'structures':
        extract_representative_structures(args.window)
    elif args.mode == 'contacts':
        contact_pairs = [(1255, 1268), (1250, 1260)]  # Ajustar según sistema
        analyze_contacts(args.window, contact_pairs)
    elif args.mode == 'secondary':
        analyze_secondary_structure(args.window)
    elif args.mode == 'geometry':
        calculate_rg_and_sasa(args.window)
    elif args.mode == 'demo':
        demonstrate_valid_vs_invalid(args.window)


if __name__ == '__main__':
    main()
