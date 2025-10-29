# 👤 Espacio de Trabajo - JUAN

## 📝 Log de Aprendizaje

### Semana 1
- [x] Instalación de OpenMM en entorno local
- [x] Configuración del entorno de trabajo (venv + requirements)
- [ ] Primer tutorial completado

### Semana 2
- [ ] Revisión de simulación clásica (tutoría 1)
- [ ] Primer benchmark en LANCAD

## 🎯 Objetivos Personales (Q4 2025)

1. Ejecutar simulaciones de dinámica molecular en la consola LANCAD para el caso de prueba de alanina.
2. Documentar cada sesión en este README y respaldar resultados en `results/`.
3. Consolidar un flujo reproducible de transferencia HPC ↔️ local.

## 💡 Notas y Reflexiones

### 2025-10-08
- Se creó el script `scripts/setup_lancad_environment.sh` para automatizar la preparación del entorno en la consola LANCAD.
- Recordatorio: mantener nombres consistentes con la carpeta `JUAN` y limpiar entornos temporales después de cada sesión HPC.
- Nuevo `requirements_bsm.txt` con dependencias de embeddings/BSM y carpeta `package/` para empaquetar el entorno y transferirlo vía SFTP.
- Asegurar copiar `mdgraphemb/` y scripts de descarga (`download_mdcath.py`, `mdgraphemb_download.py`) al paquete antes de subir a LANCAD.

## 🐛 Problemas Encontrados y Soluciones

| Problema | Solución | Fecha |
|----------|----------|-------|
| Pendiente | | |

## 🔐 Acceso HPC LANCAD

```bash
# SFTP (transferencia interactiva)
sftp l.100066@148.206.50.61
sftp> cd LUSTRE/home/lancad/2025/l.100066/Dinamicas-Moleculares-UNAM-INN/trabajo-individual/JUAN

# Empaquetar y subir carpeta JUAN (ejemplo)
tar czf JUAN.tar.gz JUAN
put JUAN.tar.gz

# SSH para ejecutar comandos
ssh l.100066@148.206.50.61
cd /LUSTRE/home/lancad/2025/l.100066/Dinamicas-Moleculares-UNAM-INN/trabajo-individual/JUAN
```

## ⚙️ Pasos para preparar consola LANCAD (solo una vez)

1. Conectarse por SSH utilizando las instrucciones anteriores.
2. Dentro del directorio `JUAN/`, ejecutar:

	```bash
	chmod +x scripts/setup_lancad_environment.sh
	./scripts/setup_lancad_environment.sh
	```

3. Activar el entorno cuando se necesite:

	```bash
	source .venv/bin/activate
	```

4. Registrar en este README cualquier paquete adicional instalado.

## 📊 Proyectos Individuales

### Proyecto 1: Alanina en solución
**Descripción**: Reproducir el benchmark de alanina del repositorio (`tutoriales-openmm/alanina`).
**Estado**: Preparación de entorno HPC.
**Archivos**: Pendiente

### Proyecto 2: Stack BSM en LANCAD
**Descripción**: Replicar librerías BSM (embeddings + PyTorch Geometric) en el cluster para ejecutar pruebas.
**Estado**: Instalación de dependencias y verificación inicial.
**Archivos**: `environment.yml`, `requirements_bsm.txt`, `scripts/setup_lancad_environment.sh`, `package/README.md`

## 🔗 Recursos Útiles Personales

- `scripts/setup_lancad_environment.sh`
- `tutoriales-openmm/`
- Documentación de módulos LANCAD
- Carpeta `package/` con instrucciones de transferencia
