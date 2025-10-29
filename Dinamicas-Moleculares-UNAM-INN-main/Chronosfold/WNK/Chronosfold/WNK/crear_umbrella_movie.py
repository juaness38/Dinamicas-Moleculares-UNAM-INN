#!/usr/bin/env python3
"""
Crea visualización animada de transición conformacional
desde resultados de umbrella sampling

Este script NO genera una trayectoria de MD clásica (eso requeriría
microsegundos de simulación), sino una ANIMACIÓN INTERPOLADA que:

1. Extrae estructuras representativas de cada ventana de umbrella
2. Interpola suavemente entre ellas
3. Genera video MP4 con overlay de PMF

Ventajas:
- Visualización clara del mecanismo estructural
- Estructuras son correctas (vienen de umbrella + MBAR)
- Relacionado con PMF (energía en cada punto)
- Duración razonable (~40 segundos)

Limitaciones:
- NO muestra cinética real (velocidad artificial)
- NO muestra fluctuaciones térmicas (suavizado)
- NO es trayectoria continua de MD

Autor: Pipeline de Umbrella Sampling
Fecha: 2024
"""

import mdtraj as md
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Backend sin display
import matplotlib.pyplot as plt
from matplotlib import animation
from sklearn.cluster import KMeans
import os
import sys
from pathlib import Path

# Configuración
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi'] = 100
plt.rcParams['font.size'] = 10

class UmbrellaMovieCreator:
    """
    Crea película animada desde resultados de umbrella sampling
    """
    
    def __init__(self, 
                 n_windows=20,
                 base_dir='.',
                 output_file='WNK1_transicion_umbrella.mp4'):
        """
        Args:
            n_windows: Número de ventanas de umbrella
            base_dir: Directorio con archivos umbrella_window_*.dcd
            output_file: Nombre del archivo de video de salida
        """
        self.n_windows = n_windows
        self.base_dir = Path(base_dir)
        self.output_file = output_file
        self.fps = 25
        self.frames_per_window = 50  # Interpolación entre ventanas
        
        # Índices de átomos Cα (se determinarán de topología)
        self.ca_nterm_index = None
        self.ca_cterm_index = None
        
        print("="*70)
        print("UMBRELLA MOVIE CREATOR")
        print("="*70)
        print(f"Ventanas: {n_windows}")
        print(f"Interpolación: {self.frames_per_window} frames/ventana")
        print(f"FPS: {self.fps}")
        print(f"Output: {output_file}")
        print("="*70)
    
    def find_ca_indices(self, topology):
        """
        Encuentra índices de Cα N-terminal y C-terminal
        
        Args:
            topology: MDTraj topology object
            
        Returns:
            (ca_nterm_index, ca_cterm_index)
        """
        # Buscar residuos terminales
        residues = list(topology.residues)
        
        # N-terminal: primer residuo con Cα
        for res in residues:
            ca_atoms = [atom.index for atom in res.atoms if atom.name == 'CA']
            if ca_atoms:
                ca_nterm = ca_atoms[0]
                break
        
        # C-terminal: último residuo con Cα
        for res in reversed(residues):
            ca_atoms = [atom.index for atom in res.atoms if atom.name == 'CA']
            if ca_atoms:
                ca_cterm = ca_atoms[0]
                break
        
        print(f"Índices Cα encontrados:")
        print(f"  N-terminal: {ca_nterm} (residuo {topology.atom(ca_nterm).residue})")
        print(f"  C-terminal: {ca_cterm} (residuo {topology.atom(ca_cterm).residue})")
        
        return ca_nterm, ca_cterm
    
    def extract_representative_structure(self, window_id, method='medoid'):
        """
        Extrae estructura más representativa de una ventana
        
        Args:
            window_id: ID de ventana (1-20)
            method: 'medoid' (más cercana al centro) o 'centroid' (K-means)
            
        Returns:
            MDTraj Trajectory con 1 frame (estructura representativa)
        """
        # Buscar archivos de trayectoria
        possible_files = [
            self.base_dir / f"umbrella_window_{window_id}.dcd",
            self.base_dir / f"window_{window_id}.dcd",
            self.base_dir / f"window{window_id}.dcd",
        ]
        
        traj_file = None
        for f in possible_files:
            if f.exists():
                traj_file = f
                break
        
        if traj_file is None:
            raise FileNotFoundError(
                f"No se encontró trayectoria para ventana {window_id}. "
                f"Buscado: {[str(f) for f in possible_files]}"
            )
        
        # Buscar topología
        possible_tops = [
            self.base_dir / "wnk1_system.pdb",
            self.base_dir / "system.pdb",
            self.base_dir / "topology.pdb",
            self.base_dir / f"umbrella_window_{window_id}_start.pdb",
        ]
        
        top_file = None
        for t in possible_tops:
            if t.exists():
                top_file = t
                break
        
        if top_file is None:
            raise FileNotFoundError(
                f"No se encontró topología. Buscado: {[str(t) for t in possible_tops]}"
            )
        
        # Cargar trayectoria
        print(f"  Cargando: {traj_file.name} (top: {top_file.name})", end='')
        traj = md.load(str(traj_file), top=str(top_file))
        print(f" - {traj.n_frames} frames")
        
        # Determinar índices Cα si es la primera vez
        if self.ca_nterm_index is None:
            self.ca_nterm_index, self.ca_cterm_index = self.find_ca_indices(traj.topology)
        
        # Alinear todas las estructuras (CA atoms)
        ca_indices = traj.topology.select('name CA')
        traj.superpose(traj[0], atom_indices=ca_indices)
        
        # Método 1: Medoid (estructura más cercana al promedio)
        if method == 'medoid':
            # Calcular RMSD de cada frame al promedio
            avg_xyz = traj.xyz.mean(axis=0)
            rmsd_to_avg = np.sqrt(((traj.xyz - avg_xyz)**2).sum(axis=(1,2)))
            representative_idx = np.argmin(rmsd_to_avg)
        
        # Método 2: K-means clustering
        elif method == 'centroid':
            # Matriz de distancias RMSD
            rmsd_matrix = np.empty((traj.n_frames, traj.n_frames))
            for i in range(traj.n_frames):
                rmsd_matrix[i] = md.rmsd(traj, traj, frame=i)
            
            # K-means con k=1
            kmeans = KMeans(n_clusters=1, random_state=42)
            kmeans.fit(rmsd_matrix)
            
            # Frame más cercano al centroide
            centroid = kmeans.cluster_centers_[0]
            distances = np.linalg.norm(rmsd_matrix - centroid, axis=1)
            representative_idx = np.argmin(distances)
        
        else:
            raise ValueError(f"Método desconocido: {method}")
        
        print(f"    → Estructura representativa: frame {representative_idx}/{traj.n_frames}")
        
        return traj[representative_idx]
    
    def interpolate_structures(self, struct1, struct2, n_frames):
        """
        Interpola linealmente entre dos estructuras
        
        Args:
            struct1: MDTraj Trajectory (1 frame)
            struct2: MDTraj Trajectory (1 frame)
            n_frames: Número de frames interpolados
            
        Returns:
            List de MDTraj Trajectories (frames interpolados)
        """
        # Alinear struct2 a struct1
        ca_indices = struct1.topology.select('name CA')
        struct2.superpose(struct1, atom_indices=ca_indices)
        
        # Coordenadas cartesianas
        xyz1 = struct1.xyz[0]  # (n_atoms, 3)
        xyz2 = struct2.xyz[0]
        
        # Interpolación lineal
        frames = []
        for alpha in np.linspace(0, 1, n_frames):
            xyz_interp = (1 - alpha) * xyz1 + alpha * xyz2
            
            # Crear nuevo frame
            frame = struct1.slice(0)
            frame.xyz[0] = xyz_interp
            frames.append(frame)
        
        return frames
    
    def calculate_cv(self, structure):
        """
        Calcula variable colectiva (distancia Cα-Cα)
        
        Args:
            structure: MDTraj Trajectory (1 frame)
            
        Returns:
            Distancia en nm
        """
        xyz = structure.xyz[0]
        ca_nterm_pos = xyz[self.ca_nterm_index]
        ca_cterm_pos = xyz[self.ca_cterm_index]
        distance = np.linalg.norm(ca_cterm_pos - ca_nterm_pos)
        return distance
    
    def load_pmf(self):
        """
        Carga PMF desde archivo MBAR
        
        Returns:
            (cv_array, pmf_array) en (nm, kJ/mol)
        """
        possible_files = [
            self.base_dir / "mbar_results.txt",
            self.base_dir / "pmf.txt",
            self.base_dir / "pmf_profile.txt",
        ]
        
        pmf_file = None
        for f in possible_files:
            if f.exists():
                pmf_file = f
                break
        
        if pmf_file is None:
            print("  ⚠ No se encontró archivo PMF. Usando PMF simulado.")
            # PMF simulado (parábola con barrera)
            cv = np.linspace(2.0, 4.0, 100)
            pmf = 25 * (cv - 2.0) * (4.0 - cv) / ((3.0-2.0)*(4.0-3.0))
            pmf -= pmf.min()
            return cv, pmf
        
        print(f"  Cargando PMF: {pmf_file.name}")
        data = np.loadtxt(pmf_file)
        
        if data.ndim == 1:
            # Formato: solo valores PMF, asumir CV equiespaciado
            cv = np.linspace(2.0, 4.0, len(data))
            pmf = data
        else:
            # Formato: columnas [CV, PMF]
            cv = data[:, 0]
            pmf = data[:, 1]
        
        print(f"    → {len(cv)} puntos, rango CV: {cv.min():.2f}-{cv.max():.2f} nm")
        print(f"    → PMF: {pmf.min():.1f} a {pmf.max():.1f} kJ/mol")
        
        return cv, pmf
    
    def render_frame(self, structure, cv_value, pmf_value, 
                    cv_pmf, pmf_pmf, frame_number, total_frames):
        """
        Renderiza un frame con proteína y overlay de PMF
        
        Args:
            structure: MDTraj Trajectory (1 frame)
            cv_value: Valor actual de CV (nm)
            pmf_value: Valor actual de PMF (kJ/mol)
            cv_pmf: Array de CV del PMF completo
            pmf_pmf: Array de PMF completo
            frame_number: Número de frame actual
            total_frames: Total de frames
            
        Returns:
            Figura matplotlib
        """
        fig = plt.figure(figsize=(14, 6))
        
        # Panel 1: Estructura 3D (simplificada - backbone)
        ax1 = fig.add_subplot(121, projection='3d')
        
        # Extraer coordenadas de Cα (backbone)
        ca_indices = structure.topology.select('name CA')
        xyz = structure.xyz[0, ca_indices, :]
        
        # Plot como ribbon
        ax1.plot(xyz[:, 0] * 10, xyz[:, 1] * 10, xyz[:, 2] * 10,
                'o-', linewidth=2.5, markersize=3, 
                color='steelblue', alpha=0.8)
        
        # Highlight N-terminal y C-terminal
        ca_nterm_idx = np.where(ca_indices == self.ca_nterm_index)[0][0]
        ca_cterm_idx = np.where(ca_indices == self.ca_cterm_index)[0][0]
        
        ax1.plot(xyz[ca_nterm_idx, 0]*10, xyz[ca_nterm_idx, 1]*10, xyz[ca_nterm_idx, 2]*10,
                'go', markersize=10, label='N-terminal', zorder=10)
        ax1.plot(xyz[ca_cterm_idx, 0]*10, xyz[ca_cterm_idx, 1]*10, xyz[ca_cterm_idx, 2]*10,
                'ro', markersize=10, label='C-terminal', zorder=10)
        
        # Línea entre N y C-terminal
        ax1.plot([xyz[ca_nterm_idx, 0]*10, xyz[ca_cterm_idx, 0]*10],
                [xyz[ca_nterm_idx, 1]*10, xyz[ca_cterm_idx, 1]*10],
                [xyz[ca_nterm_idx, 2]*10, xyz[ca_cterm_idx, 2]*10],
                'r--', linewidth=2, alpha=0.5, label=f'CV = {cv_value:.2f} nm')
        
        ax1.set_title('WNK1 Kinase C-terminal\nConformational Transition', 
                     fontsize=14, fontweight='bold')
        ax1.set_xlabel('X (Å)', fontsize=11)
        ax1.set_ylabel('Y (Å)', fontsize=11)
        ax1.set_zlabel('Z (Å)', fontsize=11)
        ax1.legend(fontsize=9, loc='upper left')
        ax1.grid(True, alpha=0.3)
        
        # Vista consistente
        ax1.view_init(elev=20, azim=45)
        
        # Panel 2: PMF con indicador de posición
        ax2 = fig.add_subplot(122)
        
        # Plot PMF completo
        ax2.plot(cv_pmf, pmf_pmf, 'k-', linewidth=2.5, label='PMF', zorder=1)
        ax2.fill_between(cv_pmf, 0, pmf_pmf, alpha=0.1, color='gray')
        
        # Indicador de posición actual
        ax2.axvline(cv_value, color='red', linestyle='--', 
                   linewidth=2, alpha=0.7, zorder=2)
        ax2.plot(cv_value, pmf_value, 'ro', markersize=20, 
                zorder=3, label='Current Position')
        
        # Anotación con valores
        ax2.annotate(f'CV = {cv_value:.2f} nm\nPMF = {pmf_value:.1f} kJ/mol',
                    xy=(cv_value, pmf_value),
                    xytext=(10, 10), textcoords='offset points',
                    fontsize=10, fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.7),
                    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0',
                                   color='red', lw=2))
        
        # Identificar mínimos y barrera
        min_idx = np.argmin(pmf_pmf)
        max_idx = np.argmax(pmf_pmf)
        
        ax2.plot(cv_pmf[min_idx], pmf_pmf[min_idx], 'go', markersize=12,
                label='Global Minimum', zorder=3)
        ax2.plot(cv_pmf[max_idx], pmf_pmf[max_idx], 'r^', markersize=12,
                label='Transition State', zorder=3)
        
        ax2.set_xlabel('Collective Variable: Cα-Cα Distance (nm)', 
                      fontsize=12, fontweight='bold')
        ax2.set_ylabel('PMF (kJ/mol)', fontsize=12, fontweight='bold')
        ax2.set_title('Potential of Mean Force\n(from Umbrella Sampling + MBAR)', 
                     fontsize=13, fontweight='bold')
        ax2.legend(fontsize=9, loc='upper right')
        ax2.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
        
        # Límites
        ax2.set_xlim(cv_pmf.min() - 0.1, cv_pmf.max() + 0.1)
        ax2.set_ylim(pmf_pmf.min() - 3, pmf_pmf.max() + 8)
        
        # Contador de frames
        fig.text(0.98, 0.02, f'Frame {frame_number}/{total_frames}',
                ha='right', va='bottom', fontsize=9, 
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        
        return fig
    
    def create_movie(self):
        """
        Pipeline completo: estructuras → interpolación → video
        """
        print("\n" + "="*70)
        print("PASO 1/5: EXTRAYENDO ESTRUCTURAS REPRESENTATIVAS")
        print("="*70)
        
        estructuras = []
        cv_values = []
        
        for i in range(1, self.n_windows + 1):
            print(f"Ventana {i}/{self.n_windows}:")
            struct = self.extract_representative_structure(i, method='medoid')
            estructuras.append(struct)
            
            # Calcular CV
            cv = self.calculate_cv(struct)
            cv_values.append(cv)
            print(f"    → CV = {cv:.3f} nm")
        
        print(f"\n✓ {len(estructuras)} estructuras extraídas")
        print(f"  Rango CV: {min(cv_values):.2f} - {max(cv_values):.2f} nm")
        
        # PASO 2: Interpolar
        print("\n" + "="*70)
        print("PASO 2/5: INTERPOLANDO ENTRE ESTRUCTURAS")
        print("="*70)
        
        all_frames = []
        all_cv_values = []
        
        for i in range(len(estructuras) - 1):
            print(f"Interpolando {i+1}→{i+2}...", end='')
            
            frames_interp = self.interpolate_structures(
                estructuras[i], 
                estructuras[i+1],
                self.frames_per_window
            )
            
            cv_interp = np.linspace(
                cv_values[i],
                cv_values[i+1],
                self.frames_per_window
            )
            
            all_frames.extend(frames_interp)
            all_cv_values.extend(cv_interp)
            
            print(f" {len(frames_interp)} frames")
        
        total_frames = len(all_frames)
        duration = total_frames / self.fps
        
        print(f"\n✓ {total_frames} frames totales")
        print(f"  Duración estimada: {duration:.1f} segundos @ {self.fps} fps")
        
        # PASO 3: Cargar PMF
        print("\n" + "="*70)
        print("PASO 3/5: CARGANDO PMF")
        print("="*70)
        
        cv_pmf, pmf_pmf = self.load_pmf()
        
        # PASO 4: Renderizar frames
        print("\n" + "="*70)
        print("PASO 4/5: RENDERIZANDO FRAMES")
        print("="*70)
        print("Este paso puede tomar 5-15 minutos dependiendo del sistema...")
        
        # Crear directorio temporal
        temp_dir = Path("temp_movie_frames")
        temp_dir.mkdir(exist_ok=True)
        print(f"Directorio temporal: {temp_dir}")
        
        for idx, (frame, cv) in enumerate(zip(all_frames, all_cv_values)):
            if idx % 50 == 0 or idx == total_frames - 1:
                progress = (idx + 1) / total_frames * 100
                print(f"  Renderizando frame {idx+1}/{total_frames} ({progress:.1f}%)...")
            
            # Interpolar PMF en CV actual
            pmf_val = np.interp(cv, cv_pmf, pmf_pmf)
            
            # Renderizar
            fig = self.render_frame(
                frame, cv, pmf_val,
                cv_pmf, pmf_pmf,
                idx + 1, total_frames
            )
            
            # Guardar
            output_path = temp_dir / f"frame_{idx:04d}.png"
            fig.savefig(output_path, dpi=100, bbox_inches='tight')
            plt.close(fig)
        
        print(f"✓ {total_frames} frames renderizados")
        
        # PASO 5: Compilar video
        print("\n" + "="*70)
        print("PASO 5/5: COMPILANDO VIDEO CON FFMPEG")
        print("="*70)
        
        # Comando FFmpeg
        ffmpeg_cmd = (
            f"ffmpeg -y -framerate {self.fps} "
            f"-i {temp_dir}/frame_%04d.png "
            f"-c:v libx264 -pix_fmt yuv420p "
            f"-preset medium -crf 23 "
            f"{self.output_file}"
        )
        
        print(f"Ejecutando: {ffmpeg_cmd[:80]}...")
        ret = os.system(ffmpeg_cmd)
        
        if ret == 0:
            print("✓ Video compilado exitosamente")
            
            # Limpiar archivos temporales
            print("\nLimpiando archivos temporales...")
            for f in temp_dir.glob("*.png"):
                f.unlink()
            temp_dir.rmdir()
            print("✓ Archivos temporales eliminados")
            
            # Verificar tamaño
            size_mb = os.path.getsize(self.output_file) / (1024 * 1024)
            print(f"\n{'='*70}")
            print("✓✓✓ PELÍCULA CREADA EXITOSAMENTE ✓✓✓")
            print(f"{'='*70}")
            print(f"Archivo: {self.output_file}")
            print(f"Tamaño: {size_mb:.1f} MB")
            print(f"Duración: {duration:.1f} segundos")
            print(f"Resolución: ~1400x600 @ {self.fps} fps")
            print(f"{'='*70}")
            
        else:
            print("✗ Error al compilar video con FFmpeg")
            print("  Los frames individuales están en:", temp_dir)
            print("  Puedes compilarlos manualmente o revisar FFmpeg")
            sys.exit(1)
        
        return self.output_file


def main():
    """Función principal"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Crea película animada desde umbrella sampling',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos de uso:

  # Uso básico (20 ventanas, directorio actual)
  python crear_umbrella_movie.py
  
  # Especificar número de ventanas
  python crear_umbrella_movie.py --windows 30
  
  # Especificar directorio de datos
  python crear_umbrella_movie.py --dir ../resultados_umbrella/
  
  # Nombre de archivo personalizado
  python crear_umbrella_movie.py --output mi_proteina.mp4
  
  # Más frames de interpolación (suavidad)
  python crear_umbrella_movie.py --frames-per-window 100

Requisitos de archivos:
  - umbrella_window_*.dcd (trayectorias de cada ventana)
  - wnk1_system.pdb o system.pdb (topología)
  - mbar_results.txt (PMF) [opcional, se simulará si falta]

Dependencias:
  - mdtraj
  - numpy
  - matplotlib
  - scikit-learn
  - ffmpeg (sistema)
        """
    )
    
    parser.add_argument('--windows', type=int, default=20,
                       help='Número de ventanas de umbrella (default: 20)')
    parser.add_argument('--dir', type=str, default='.',
                       help='Directorio con archivos de umbrella (default: .)')
    parser.add_argument('--output', type=str, default='WNK1_transicion_umbrella.mp4',
                       help='Archivo de salida (default: WNK1_transicion_umbrella.mp4)')
    parser.add_argument('--frames-per-window', type=int, default=50,
                       help='Frames interpolados entre ventanas (default: 50)')
    parser.add_argument('--fps', type=int, default=25,
                       help='Frames por segundo del video (default: 25)')
    
    args = parser.parse_args()
    
    # Crear instancia
    creator = UmbrellaMovieCreator(
        n_windows=args.windows,
        base_dir=args.dir,
        output_file=args.output
    )
    
    creator.fps = args.fps
    creator.frames_per_window = args.frames_per_window
    
    # Ejecutar
    try:
        video_file = creator.create_movie()
        print(f"\n✓ Proceso completado. Video: {video_file}")
        return 0
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
