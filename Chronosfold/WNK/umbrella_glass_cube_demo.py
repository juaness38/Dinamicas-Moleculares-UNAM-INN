#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
UMBRELLA GLASS CUBE - DEMO VERSION (Standalone)
No requiere mdtraj - usa datos sintéticos completamente

Este demo muestra cómo funcionaría el Glass Cube con datos reales
"""

import numpy as np
import plotly.graph_objects as go
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class WindowPhysicsDemo(object):
    """Métricas físicas simuladas de una ventana"""
    def __init__(self, window_id, n_windows=20):
        self.window_id = window_id
        
        # CV center
        self.cv_center = 2.0 + (window_id - 1) * (4.0 - 2.0) / (n_windows - 1)
        
        # PMF simulado (barrera de 25 kJ/mol en CV=3.0 nm)
        cv_range = np.linspace(2.0, 4.0, 100)
        pmf_curve = 25 * (cv_range - 2.0) * (4.0 - cv_range) / ((3.0-2.0)*(4.0-3.0))
        pmf_curve -= pmf_curve.min()
        self.pmf_value = np.interp(self.cv_center, cv_range, pmf_curve)
        
        # ENERGÍAS (kJ/mol) - más alto cerca del TS
        base_energy = -25000 + self.pmf_value * 1000
        self.energy_mean = base_energy
        self.energy_std = 500
        self.energy_min = base_energy - 1000
        self.energy_max = base_energy + 1000
        
        # FUERZAS (pN) - más alto cerca del TS
        self.force_mean = 30 + self.pmf_value * 2.5
        self.force_std = 10
        self.force_max = self.force_mean + 30
        
        # RMSD (nm) - aumenta con la apertura
        self.rmsd_mean = 0.10 + (self.cv_center - 2.0) * 0.08
        self.rmsd_max = self.rmsd_mean * 2
        self.rmsd_values = np.random.normal(self.rmsd_mean, 0.03, 1000)
        
        # RADIUS OF GYRATION (nm) - aumenta al abrirse
        self.rg_mean = 1.75 + (self.cv_center - 2.0) * 0.15
        self.rg_std = 0.04
        self.rg_values = np.random.normal(self.rg_mean, self.rg_std, 1000)
        
        # PCA (exploración conformacional)
        n_frames = 1000
        # Windows cerca del TS exploran más
        exploration_factor = 1.0 + (1.0 - abs(self.cv_center - 3.0))
        self.pca_x = np.random.normal(0, exploration_factor, n_frames)
        self.pca_y = np.random.normal(0, exploration_factor * 0.7, n_frames)
        self.pca_explained_variance = [0.45, 0.28]
        
        # RAMACHANDRAN (phi/psi angles)
        # Simular diferentes regiones para diferentes conformaciones
        if self.cv_center < 2.5:  # Cerrado - región beta
            phi_center, psi_center = -120, 120
        elif self.cv_center > 3.5:  # Abierto - región alpha
            phi_center, psi_center = -60, -45
        else:  # Transición - mezcla
            phi_center, psi_center = -90, 60
        
        n_residues = 50
        self.phi_angles = np.random.normal(phi_center, 30, (n_frames, n_residues))
        self.psi_angles = np.random.normal(psi_center, 30, (n_frames, n_residues))
        
        # PROTEÍNA 3D (estructura sintética)
        # Simulamos backbone de proteína que se "abre"
        n_residues_3d = 100
        t = np.linspace(0, 4*np.pi, n_residues_3d)
        
        # Apertura aumenta con CV
        opening_factor = (self.cv_center - 2.0) / 2.0  # 0 a 1
        
        # Forma de hélice que se abre
        self.protein_x = np.cos(t) * (1 + opening_factor * 0.5)
        self.protein_y = np.sin(t) * (1 + opening_factor * 0.5)
        self.protein_z = t / (4*np.pi) * 3 - 1.5
        
        logger.info(f"  Window {window_id}: CV={self.cv_center:.2f} nm, PMF={self.pmf_value:.1f} kJ/mol")


class UmbrellaGlassCubeDemo:
    """Glass Cube con datos completamente sintéticos"""
    
    def __init__(self, n_windows=20):
        self.n_windows = n_windows
        self.window_physics = {}
        self.global_ranges = {}
        
        # ChronosGPT Colors
        self.colors = {
            "background": "black",
            "text": "#00FFFF",
            "grid": "#333333",
            "protein": "#4ECDC4",
            "highlight": "#FFD93D",
            'thermometer_energy': '#FFD700',
            'thermometer_force': '#9370DB',
            'thermometer_rmsd': '#FFB6C1',
            'thermometer_rg': '#00FFFF'
        }
        
        logger.info("="*70)
        logger.info("🧊 UMBRELLA GLASS CUBE - DEMO MODE")
        logger.info("="*70)
        logger.info("Generando 20 ventanas con datos sintéticos...")
        logger.info("="*70)
        
    def generate_synthetic_data(self):
        """Genera datos sintéticos para todas las ventanas"""
        logger.info("\nGenerando datos físicos para cada ventana...\n")
        
        for window_id in range(1, self.n_windows + 1):
            physics = WindowPhysicsDemo(window_id, self.n_windows)
            self.window_physics[window_id] = physics
        
        # Calcular rangos globales
        energies = [w.energy_mean for w in self.window_physics.values()]
        forces = [w.force_mean for w in self.window_physics.values()]
        rmsds = [w.rmsd_mean for w in self.window_physics.values()]
        rgs = [w.rg_mean for w in self.window_physics.values()]
        
        self.global_ranges = {
            'energy': (min(energies), max(energies)),
            'force': (min(forces), max(forces)),
            'rmsd': (min(rmsds), max(rmsds)),
            'rg': (min(rgs), max(rgs))
        }
        
        logger.info(f"\n{'='*70}")
        logger.info("Rangos globales para normalización de thermometers:")
        logger.info(f"  🟡 Energy: {self.global_ranges['energy'][0]:.0f} - {self.global_ranges['energy'][1]:.0f} kJ/mol")
        logger.info(f"  🟣 Force: {self.global_ranges['force'][0]:.1f} - {self.global_ranges['force'][1]:.1f} pN")
        logger.info(f"  🌸 RMSD: {self.global_ranges['rmsd'][0]:.3f} - {self.global_ranges['rmsd'][1]:.3f} nm")
        logger.info(f"  🔵 Rg: {self.global_ranges['rg'][0]:.3f} - {self.global_ranges['rg'][1]:.3f} nm")
        logger.info(f"{'='*70}\n")
    
    def create_glass_cube(self, output_file="demo_glass_cube.html"):
        """Crea la visualización Glass Cube"""
        logger.info("Creando visualización Glass Cube...")
        
        fig = go.Figure()
        
        # Crear frames (1 por ventana)
        frames = []
        for window_id in range(1, self.n_windows + 1):
            logger.info(f"  Frame {window_id}/{self.n_windows}...")
            frame_data = self._create_window_frame(window_id)
            frames.append(frame_data)
        
        fig.frames = frames
        
        # Frame inicial
        initial_traces = self._get_traces_for_window(1)
        for trace in initial_traces:
            fig.add_trace(trace)
        
        # Styling
        self._apply_styling(fig)
        
        # Guardar
        fig.write_html(output_file)
        
        logger.info(f"\n{'='*70}")
        logger.info("✅ VISUALIZACIÓN CREADA EXITOSAMENTE")
        logger.info(f"{'='*70}")
        logger.info(f"📁 Archivo: {output_file}")
        logger.info(f"🪟 Windows: {len(self.window_physics)}")
        logger.info(f"\n🎮 CONTROLES:")
        logger.info(f"  ▶️  Play/Pause: Botones en pantalla")
        logger.info(f"  🎚️  Slider: Navega entre ventanas manualmente")
        logger.info(f"  🖱️  Mouse: Rota, zoom, pan del cubo 3D")
        logger.info(f"{'='*70}\n")
        
        return output_file
    
    def _create_window_frame(self, window_id):
        """Crea frame para una ventana"""
        traces = self._get_traces_for_window(window_id)
        physics = self.window_physics[window_id]
        
        return go.Frame(
            data=traces,
            name=f'Window_{window_id}',
            layout=go.Layout(
                title=dict(
                    text=f"🪟 Window {window_id}/20 | CV = {physics.cv_center:.2f} nm | PMF = {physics.pmf_value:.1f} kJ/mol",
                    font=dict(size=18, color="#00FFFF")
                )
            )
        )
    
    def _get_traces_for_window(self, window_id):
        """Genera trazas para una ventana"""
        physics = self.window_physics[window_id]
        traces = []
        
        # 1. Ramachandran floor
        traces.extend(self._create_ramachandran_floor(physics))
        
        # 2. PCA layer
        traces.extend(self._create_pca_layer(physics))
        
        # 3. Protein structure
        traces.extend(self._create_protein_structure(physics))
        
        # 4. Thermometers
        traces.extend(self._create_thermometers(physics))
        
        # 5. PMF indicator
        traces.append(self._create_pmf_indicator(physics))
        
        return traces
    
    def _create_ramachandran_floor(self, physics):
        """Ramachandran plot en el piso"""
        floor_z = -2.5
        
        # Grid
        phi_bins = np.linspace(-180, 180, 30)
        psi_bins = np.linspace(-180, 180, 30)
        
        # Histograma 2D
        phi_flat = physics.phi_angles.flatten()
        psi_flat = physics.psi_angles.flatten()
        
        H, _, _ = np.histogram2d(phi_flat, psi_flat, bins=[phi_bins, psi_bins])
        
        # Grid para plot
        phi_centers = (phi_bins[:-1] + phi_bins[1:]) / 2
        psi_centers = (psi_bins[:-1] + psi_bins[1:]) / 2
        phi_mesh, psi_mesh = np.meshgrid(phi_centers, psi_centers)
        
        # Escalar a coordenadas del cubo
        phi_scaled = phi_mesh / 180 * 3
        psi_scaled = psi_mesh / 180 * 3
        
        trace = go.Surface(
            x=phi_scaled,
            y=psi_scaled,
            z=np.full_like(phi_scaled, floor_z),
            surfacecolor=H.T,
            colorscale='Viridis',
            opacity=0.5,
            showscale=False,
            name='📊 Ramachandran',
            hovertemplate='φ: %{customdata[0]:.0f}°<br>ψ: %{customdata[1]:.0f}°<extra></extra>',
            customdata=np.dstack([phi_mesh, psi_mesh])
        )
        
        return [trace]
    
    def _create_pca_layer(self, physics):
        """PCA como firefly points"""
        pca_z = 1.0
        
        # Normalizar
        pca_x_scaled = physics.pca_x / (np.abs(physics.pca_x).max() + 0.01) * 2.5
        pca_y_scaled = physics.pca_y / (np.abs(physics.pca_y).max() + 0.01) * 2.5
        
        # Detectar conformaciones "novel"
        distances = np.sqrt(pca_x_scaled**2 + pca_y_scaled**2)
        threshold = np.percentile(distances, 75)
        is_novel = distances > threshold
        
        colors = ['#FFD93D' if n else '#4ECDC4' for n in is_novel]
        sizes = [10 if n else 5 for n in is_novel]
        
        trace = go.Scatter3d(
            x=pca_x_scaled,
            y=pca_y_scaled,
            z=np.full(len(pca_x_scaled), pca_z),
            mode='markers',
            marker=dict(color=colors, size=sizes, opacity=0.6),
            name=f'🔴 PCA ({physics.pca_explained_variance[0]:.0%} var)',
            hovertemplate='PC1: %{x:.2f}<br>PC2: %{y:.2f}<extra></extra>'
        )
        
        return [trace]
    
    def _create_protein_structure(self, physics):
        """Estructura 3D sintética"""
        protein_z_offset = 3.25
        
        # Normalizar coordenadas
        coords = np.column_stack([physics.protein_x, physics.protein_y, physics.protein_z])
        center = coords.mean(axis=0)
        coords_centered = coords - center
        max_range = np.abs(coords_centered).max()
        if max_range > 0:
            coords_scaled = coords_centered / max_range * 2.0
        else:
            coords_scaled = coords_centered
        
        # Color gradient (N-term azul → C-term rojo)
        n_res = len(coords_scaled)
        colors_gradient = np.linspace(0, 1, n_res)
        
        trace = go.Scatter3d(
            x=coords_scaled[:, 0],
            y=coords_scaled[:, 1],
            z=coords_scaled[:, 2] + protein_z_offset,
            mode='lines+markers',
            line=dict(color=self.colors['protein'], width=6),
            marker=dict(size=6, color=colors_gradient, colorscale='Turbo', showscale=False),
            name='🧬 Protein',
            hovertemplate='Residue %{pointNumber}<extra></extra>'
        )
        
        return [trace]
    
    def _create_thermometers(self, physics):
        """Edge thermometers con datos reales"""
        cube_range = 3
        traces = []
        
        # Normalizar
        def normalize(value, range_tuple):
            vmin, vmax = range_tuple
            if vmax - vmin == 0:
                return 0.5
            norm = (value - vmin) / (vmax - vmin)
            return max(0, min(1, norm))
        
        energy_norm = normalize(physics.energy_mean, self.global_ranges['energy'])
        force_norm = normalize(physics.force_mean, self.global_ranges['force'])
        rmsd_norm = normalize(physics.rmsd_mean, self.global_ranges['rmsd'])
        rg_norm = normalize(physics.rg_mean, self.global_ranges['rg'])
        
        # Alturas
        def get_height(norm):
            return -cube_range + (2 * cube_range * norm)
        
        # 🟡 ENERGY (front-left)
        energy_h = get_height(energy_norm)
        traces.append(go.Scatter3d(
            x=[-cube_range, -cube_range],
            y=[-cube_range, -cube_range],
            z=[-cube_range, energy_h],
            mode='lines+markers',
            line=dict(color=self.colors['thermometer_energy'], width=10),
            marker=dict(color=self.colors['thermometer_energy'], size=8),
            name=f'🟡 Energy: {physics.energy_mean:.0f} kJ/mol',
            hovertemplate=f'E={physics.energy_mean:.0f} kJ/mol<br>Level: {energy_norm:.0%}<extra></extra>'
        ))
        
        # 🟣 FORCE (front-right)
        force_h = get_height(force_norm)
        traces.append(go.Scatter3d(
            x=[cube_range, cube_range],
            y=[-cube_range, -cube_range],
            z=[-cube_range, force_h],
            mode='lines+markers',
            line=dict(color=self.colors['thermometer_force'], width=10),
            marker=dict(color=self.colors['thermometer_force'], size=8),
            name=f'🟣 Force: {physics.force_mean:.1f} pN',
            hovertemplate=f'F={physics.force_mean:.1f} pN<br>Level: {force_norm:.0%}<extra></extra>'
        ))
        
        # 🌸 RMSD (back-left)
        rmsd_h = get_height(rmsd_norm)
        traces.append(go.Scatter3d(
            x=[-cube_range, -cube_range],
            y=[cube_range, cube_range],
            z=[-cube_range, rmsd_h],
            mode='lines+markers',
            line=dict(color=self.colors['thermometer_rmsd'], width=10),
            marker=dict(color=self.colors['thermometer_rmsd'], size=8),
            name=f'🌸 RMSD: {physics.rmsd_mean:.3f} nm',
            hovertemplate=f'RMSD={physics.rmsd_mean:.3f} nm<br>Level: {rmsd_norm:.0%}<extra></extra>'
        ))
        
        # 🔵 Rg (back-right)
        rg_h = get_height(rg_norm)
        traces.append(go.Scatter3d(
            x=[cube_range, cube_range],
            y=[cube_range, cube_range],
            z=[-cube_range, rg_h],
            mode='lines+markers',
            line=dict(color=self.colors['thermometer_rg'], width=10),
            marker=dict(color=self.colors['thermometer_rg'], size=8),
            name=f'🔵 Rg: {physics.rg_mean:.3f} nm',
            hovertemplate=f'Rg={physics.rg_mean:.3f} nm<br>Level: {rg_norm:.0%}<extra></extra>'
        ))
        
        return traces
    
    def _create_pmf_indicator(self, physics):
        """Indicador de posición en PMF"""
        # Posición en arista bottom-front
        cv_scaled = (physics.cv_center - 2.0) / (4.0 - 2.0) * 6 - 3
        pmf_scaled = physics.pmf_value / 25.0 * 6 - 3
        
        return go.Scatter3d(
            x=[cv_scaled],
            y=[-3],
            z=[pmf_scaled],
            mode='markers',
            marker=dict(
                size=18,
                color='red',
                symbol='diamond',
                line=dict(color='white', width=3)
            ),
            name=f'📍 PMF Position',
            hovertemplate=f'CV: {physics.cv_center:.2f} nm<br>PMF: {physics.pmf_value:.1f} kJ/mol<extra></extra>'
        )
    
    def _apply_styling(self, fig):
        """Aplica styling del Glass Cube"""
        fig.update_scenes(
            bgcolor="black",
            xaxis=dict(
                title="X", gridcolor="#333333", color="#00FFFF",
                backgroundcolor="black", showbackground=True, range=[-3.5, 3.5]
            ),
            yaxis=dict(
                title="Y", gridcolor="#333333", color="#00FFFF",
                backgroundcolor="black", showbackground=True, range=[-3.5, 3.5]
            ),
            zaxis=dict(
                title="Z: Floor → PCA → Protein", gridcolor="#333333", color="#00FFFF",
                backgroundcolor="black", showbackground=True, range=[-3, 5]
            ),
            camera=dict(eye=dict(x=1.5, y=1.5, z=1.2)),
            aspectratio=dict(x=1, y=1, z=1)
        )
        
        fig.update_layout(
            title=dict(
                text="🧊 UMBRELLA GLASS CUBE - Physical Analysis per Window",
                font=dict(size=22, color="#00FFFF", family="Arial Black"),
                x=0.5
            ),
            showlegend=True,
            legend=dict(
                x=1.02, y=0.98,
                bgcolor="rgba(0,0,0,0.9)",
                bordercolor="#00FFFF",
                borderwidth=2,
                font=dict(color="#00FFFF", size=11)
            ),
            width=1600,
            height=1000,
            margin=dict(l=30, r=200, t=80, b=30),
            paper_bgcolor="black",
            plot_bgcolor="black",
            font=dict(color="#00FFFF", size=12),
            updatemenus=[{
                "type": "buttons",
                "showactive": False,
                "buttons": [
                    {
                        "label": "▶️ Play",
                        "method": "animate",
                        "args": [None, {
                            "frame": {"duration": 1000, "redraw": True},
                            "fromcurrent": True,
                            "mode": "immediate"
                        }]
                    },
                    {
                        "label": "⏸️ Pause",
                        "method": "animate",
                        "args": [[None], {
                            "frame": {"duration": 0, "redraw": False},
                            "mode": "immediate"
                        }]
                    }
                ],
                "x": 0.1, "y": 0.02,
                "xanchor": "left", "yanchor": "bottom",
                "bgcolor": "#00FFFF",
                "font": dict(size=14, color="black")
            }],
            sliders=[{
                "active": 0,
                "steps": [
                    {
                        "method": "animate",
                        "args": [[f"Window_{i}"], {
                            "frame": {"duration": 0, "redraw": True},
                            "mode": "immediate"
                        }],
                        "label": f"W{i}"
                    }
                    for i in range(1, self.n_windows + 1)
                ],
                "x": 0.1, "y": 0.0,
                "len": 0.8,
                "xanchor": "left", "yanchor": "top",
                "bgcolor": "#333333",
                "bordercolor": "#00FFFF",
                "borderwidth": 2,
                "font": dict(color="#00FFFF", size=10)
            }]
        )


def main():
    """Demo principal"""
    print("\n" + "="*70)
    print("🧊 UMBRELLA GLASS CUBE - DEMO VERSION")
    print("="*70)
    print("Visualización 3D interactiva de Umbrella Sampling")
    print("Usando datos sintéticos para demostración")
    print("="*70 + "\n")
    
    try:
        # Crear instancia
        ugc = UmbrellaGlassCubeDemo(n_windows=20)
        
        # Generar datos sintéticos
        ugc.generate_synthetic_data()
        
        # Crear visualización
        output = ugc.create_glass_cube("demo_glass_cube.html")
        
        print(f"\n🎉 ¡LISTO! Abre el archivo en tu navegador:")
        print(f"   👉 {output}\n")
        
        return 0
        
    except Exception as e:
        logger.error(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    import sys
    sys.exit(main())
