# 🛠️ Recursos para Aprendizaje

Esta carpeta contiene scripts, herramientas y recursos útiles para el aprendizaje de dinámicas moleculares.

## 📄 Archivos Disponibles

### `verificar_instalacion.py`
Script para verificar que todas las dependencias estén correctamente instaladas.

**Uso:**
```bash
python knowledge-base/recursos/verificar_instalacion.py
```

### `simulacion_alanina.py`
Script completo para ejecutar una simulación de dinámica molecular de alanina dipéptido.

**Uso:**
```bash
python knowledge-base/recursos/simulacion_alanina.py
```

**Características:**
- Descarga automática de estructura PDB
- Minimización de energía
- Equilibración del sistema
- Simulación de producción (20 ps)
- Generación de trayectoria DCD

### `tutorial_alanina_dipeptido.ipynb`
Jupyter notebook interactivo con tutorial COMPLETO de simulación y análisis.

**Contenido:**
- ✅ Simulación paso a paso con explicaciones detalladas
- ✅ Análisis de energías
- ✅ Cálculo de RMSD
- ✅ Diagrama de Ramachandran
- ✅ Análisis PCA (componentes principales)
- ✅ Visualización 3D interactiva
- ✅ Generación de gráficos de calidad

**Uso:**
```bash
jupyter lab knowledge-base/recursos/tutorial_alanina_dipeptido.ipynb
```

## 🔗 Links Útiles por Categoría

### 📚 Documentación Oficial
- [OpenMM User Guide](http://docs.openmm.org/latest/userguide/)
- [OpenMM API Reference](http://docs.openmm.org/latest/api-python/)
- [MDTraj Documentation](http://mdtraj.org/)
- [MDAnalysis Documentation](https://docs.mdanalysis.org/)

### 🎓 Tutoriales y Cursos
- [OpenMM Tutorial MSBS](https://github.com/molmod/openmm-tutorial-msbs) - Tutorial base
- [OpenMM Cookbook](http://docs.openmm.org/latest/userguide/cookbook.html)
- [MDAnalysis Tutorial](https://userguide.mdanalysis.org/stable/examples/README.html)
- [Molecular Dynamics Basics](http://www.ks.uiuc.edu/Training/Tutorials/)

### 🎬 Videos Educativos
- [Introduction to Molecular Dynamics](https://www.youtube.com/results?search_query=molecular+dynamics+simulation+tutorial)
- [OpenMM Webinars](https://openmm.org/documentation)

### 📖 Libros y Papers
- Frenkel & Smit - Understanding Molecular Simulation
- Allen & Tildesley - Computer Simulation of Liquids
- Leach - Molecular Modelling: Principles and Applications

### 🧰 Herramientas de Visualización
- [VMD](https://www.ks.uiuc.edu/Research/vmd/) - Visual Molecular Dynamics
- [PyMOL](https://pymol.org/) - Visualización de proteínas
- [NGLView](https://github.com/nglviewer/nglview) - Visualización en Jupyter
- [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) - Visualización molecular moderna

### 🗄️ Bases de Datos
- [Protein Data Bank (PDB)](https://www.rcsb.org/) - Estructuras de proteínas
- [CHARMM-GUI](https://www.charmm-gui.org/) - Preparación de sistemas
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/) - Moléculas pequeñas

### 💬 Comunidades y Foros
- [OpenMM Discussions](https://github.com/openmm/openmm/discussions)
- [CCL - Computational Chemistry List](http://www.ccl.net/)
- [r/comp_chem](https://www.reddit.com/r/comp_chem/) - Subreddit

### 🔧 Force Fields (Campos de Fuerza)
- [AMBER Force Fields](http://ambermd.org/)
- [CHARMM Force Fields](https://www.charmm.org/)
- [GROMOS](http://www.gromos.net/)

## 📝 Plantillas y Scripts

### Script Básico de Simulación
Ver: `tutoriales-openmm/` para ejemplos completos

### Análisis de Trayectorias
Próximamente: templates de análisis RMSD, RMSF, energías, etc.

## 🎯 Recursos Recomendados por Fase

### Fase 1: Fundamentos
- [ ] Leer: Introduction to Molecular Dynamics (primer capítulo de Frenkel & Smit)
- [ ] Ver: Video introductorio de MD
- [ ] Instalar: Todas las herramientas necesarias

### Fase 2: Primeros Pasos
- [ ] Tutorial: OpenMM First Steps
- [ ] Herramienta: NGLView para visualización
- [ ] Práctica: Sistema simple (agua)

### Fase 3: Intermedio
- [ ] Tutorial: Protein in Water
- [ ] Herramienta: MDTraj para análisis
- [ ] Práctica: Simulación de proteína

### Fase 4: Avanzado
- [ ] Paper: Enhanced Sampling Methods
- [ ] Herramienta: Free Energy Calculations
- [ ] Práctica: Técnicas especializadas

## 🆕 Agregar Nuevos Recursos

Cuando encuentres un recurso útil:
1. Agrégalo a este README en la categoría apropiada
2. Si es un archivo, súbelo a esta carpeta
3. Documenta por qué es útil
4. Compártelo con el equipo

---

**Última actualización**: Octubre 2, 2025
