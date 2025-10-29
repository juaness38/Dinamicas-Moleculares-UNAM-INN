# Paquete de Transferencia HPC (JUAN)

Este directorio agrupa los archivos mínimos que deben copiarse vía SFTP/FileZilla hacia la consola LANCAD (`/LUSTRE/home/lancad/2025/l.100066/Dinamicas-Moleculares-UNAM-INN/trabajo-individual/JUAN`).

## Contenido recomendado

| Archivo | Descripción |
| --- | --- |
| `environment.yml` | Entorno conda completo (`bsm-lancad-env`) incluyendo PyTorch, PyTorch Geometric, OpenMM y dependencias BSM. |
| `requirements.txt` | Requerimientos base del curso (MD). Copiar desde la raíz del repositorio original si aún no existe en la carpeta `JUAN/`. |
| `requirements_bsm.txt` | Paquetes pip adicionales para BSM (embeddings, API, conectores). |
| `scripts/setup_lancad_environment.sh` | Script automatizado para preparar entorno virtual y realizar la instalación completa. |
| `README.md` | Bitácora personal (opcional) para seguir actualizando desde el cluster. |
| `mdgraphemb/` | Scripts/notebooks para extracción de ESE (copiar desde raíz si se requieren). |
| `scripts/mdcath_download.py` | Script de descarga mdCATH necesario para poblar datasets ESE (ver sección de instrucciones). |

> **Nota:** `environment.yml` y `requirements.txt` viven en la raíz del proyecto. Antes de comprimir/transferir, cópialos dentro de este directorio para enviar un paquete autocontenido.

## Flujo sugerido (local → LANCAD)

```powershell
# Dentro del repositorio local
cd ...\trabajo-individual\JUAN
robocopy .\package .\package_tmp environment.yml requirements.txt scripts\setup_lancad_environment.sh requirements_bsm.txt
Compress-Archive -Path .\package_tmp\* -DestinationPath JUAN_hpc_package.zip
# Subir `JUAN_hpc_package.zip` con FileZilla
```

En el cluster:

```bash
cd /LUSTRE/home/lancad/2025/l.100066/Dinamicas-Moleculares-UNAM-INN/trabajo-individual/JUAN
unzip ~/JUAN_hpc_package.zip
chmod +x scripts/setup_lancad_environment.sh
./scripts/setup_lancad_environment.sh
```

## Variables importantes

- `TORCH_INDEX_URL` y `PYG_INDEX_URL` pueden ajustarse antes de ejecutar el script si LANCAD cambia de versión CUDA.
- Resultados, notebooks y logs se guardan automáticamente en `JUAN/results`, `JUAN/notebooks` y `JUAN/logs`.

Mantén este paquete actualizado cada vez que agregues scripts o dependencias nuevas necesarias para ejecutar BSM en el cluster.
