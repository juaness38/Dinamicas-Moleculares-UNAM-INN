#!/bin/bash
# ============================================================================
# Script de despliegue SSH para HPC
#
# Este script automatiza:
# 1. Conexión SSH al cluster
# 2. Transferencia de archivos
# 3. Preparación del entorno
# 4. Ejecución de tests pre-HPC
# 5. Submisión del job array
#
# Uso:
#   bash deploy_to_hpc.sh
#
# AJUSTAR CONFIGURACIÓN ABAJO
# ============================================================================

# ============================================================================
# CONFIGURACIÓN - AJUSTAR ESTOS VALORES
# ============================================================================

# Credenciales HPC
HPC_USER="tu_usuario"              # Tu usuario en el cluster
HPC_HOST="nombre_cluster.edu"      # Hostname del cluster
HPC_PORT=22                         # Puerto SSH (usualmente 22)

# Directorios
LOCAL_DIR="$(pwd)"                             # Directorio local (donde está este script)
HPC_WORK_DIR="/home/${HPC_USER}/wnk_umbrella"  # Directorio de trabajo en HPC

# Entorno Python en HPC
HPC_PYTHON_ENV="/path/to/your/venv"            # Path al virtualenv en HPC
# o usar: HPC_PYTHON_ENV="module load python/3.11"  # si usas modules

# Configuración SLURM (opcional, sobrescribe valores en submit_umbrella_hpc.sh)
SLURM_PARTITION="gpu"              # Partición a usar
SLURM_TIME="48:00:00"              # Tiempo máximo
SLURM_MEM="16G"                    # Memoria por job

# ============================================================================
# COLORES PARA OUTPUT
# ============================================================================

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# ============================================================================
# FUNCIONES
# ============================================================================

print_header() {
    echo -e "${BLUE}================================================================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}================================================================================${NC}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ ERROR: $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠ WARNING: $1${NC}"
}

check_ssh_connection() {
    print_header "VERIFICANDO CONEXIÓN SSH"
    
    if ssh -p ${HPC_PORT} -o ConnectTimeout=10 ${HPC_USER}@${HPC_HOST} "echo 'Connection OK'" &>/dev/null; then
        print_success "Conexión SSH exitosa"
        return 0
    else
        print_error "No se pudo conectar a ${HPC_USER}@${HPC_HOST}:${HPC_PORT}"
        echo "Verifica:"
        echo "  1. Credenciales (HPC_USER, HPC_HOST)"
        echo "  2. VPN activa (si es necesaria)"
        echo "  3. Llaves SSH configuradas"
        return 1
    fi
}

create_remote_directory() {
    print_header "CREANDO DIRECTORIO REMOTO"
    
    ssh -p ${HPC_PORT} ${HPC_USER}@${HPC_HOST} "mkdir -p ${HPC_WORK_DIR}/logs"
    
    if [ $? -eq 0 ]; then
        print_success "Directorio creado: ${HPC_WORK_DIR}"
    else
        print_error "No se pudo crear directorio remoto"
        return 1
    fi
}

transfer_files() {
    print_header "TRANSFIRIENDO ARCHIVOS"
    
    echo "Empaquetando archivos locales..."
    
    # Crear tarball con los archivos necesarios
    tar -czf /tmp/wnk_umbrella.tar.gz \
        Chronosfold/WNK/*.py \
        Chronosfold/WNK/*.sh \
        Chronosfold/WNK/5DRB.pdb \
        Chronosfold/tests/test_wnk_umbrella_setup.py \
        2>/dev/null
    
    if [ ! -f /tmp/wnk_umbrella.tar.gz ]; then
        print_error "No se pudo crear tarball"
        return 1
    fi
    
    print_success "Tarball creado: $(du -h /tmp/wnk_umbrella.tar.gz | cut -f1)"
    
    echo "Transfiriendo a HPC..."
    
    scp -P ${HPC_PORT} /tmp/wnk_umbrella.tar.gz ${HPC_USER}@${HPC_HOST}:${HPC_WORK_DIR}/
    
    if [ $? -eq 0 ]; then
        print_success "Transferencia completada"
        
        # Desempaquetar en HPC
        echo "Desempaquetando en HPC..."
        ssh -p ${HPC_PORT} ${HPC_USER}@${HPC_HOST} "cd ${HPC_WORK_DIR} && tar -xzf wnk_umbrella.tar.gz"
        
        print_success "Archivos desempaquetados"
        
        # Limpiar tarball local
        rm /tmp/wnk_umbrella.tar.gz
    else
        print_error "Fallo en transferencia SCP"
        return 1
    fi
}

run_preparation() {
    print_header "EJECUTANDO PREPARACIÓN DEL SISTEMA (REMOTO)"
    
    print_warning "Este paso puede tomar 10-30 minutos"
    echo "Ejecutando prepare_system.py en HPC..."
    
    ssh -p ${HPC_PORT} ${HPC_USER}@${HPC_HOST} << EOF
cd ${HPC_WORK_DIR}

# Activar entorno
source ${HPC_PYTHON_ENV}/bin/activate 2>/dev/null || module load python/3.11

# Ejecutar preparación
python3 Chronosfold/WNK/prepare_system.py 2>&1 | tee logs/prepare_system.log

exit \${PIPESTATUS[0]}
EOF
    
    if [ $? -eq 0 ]; then
        print_success "Sistema preparado exitosamente"
    else
        print_error "Fallo en preparación del sistema"
        echo "Revisa logs/prepare_system.log en el HPC"
        return 1
    fi
}

run_window_generation() {
    print_header "GENERANDO VENTANAS DE UMBRELLA SAMPLING (REMOTO)"
    
    echo "Ejecutando generate_umbrella_windows.py en HPC..."
    
    ssh -p ${HPC_PORT} ${HPC_USER}@${HPC_HOST} << EOF
cd ${HPC_WORK_DIR}

# Activar entorno
source ${HPC_PYTHON_ENV}/bin/activate 2>/dev/null || module load python/3.11

# Generar ventanas
python3 Chronosfold/WNK/generate_umbrella_windows.py 2>&1 | tee logs/generate_windows.log

exit \${PIPESTATUS[0]}
EOF
    
    if [ $? -eq 0 ]; then
        print_success "Ventanas generadas exitosamente"
    else
        print_error "Fallo en generación de ventanas"
        return 1
    fi
}

run_pre_hpc_tests() {
    print_header "EJECUTANDO TESTS PRE-HPC (REMOTO)"
    
    echo "Ejecutando test_wnk_umbrella_setup.py en HPC..."
    
    ssh -p ${HPC_PORT} ${HPC_USER}@${HPC_HOST} << EOF
cd ${HPC_WORK_DIR}

# Activar entorno
source ${HPC_PYTHON_ENV}/bin/activate 2>/dev/null || module load python/3.11

# Ejecutar tests
python3 -m pytest Chronosfold/tests/test_wnk_umbrella_setup.py -v 2>&1 | tee logs/pre_hpc_tests.log

exit \${PIPESTATUS[0]}
EOF
    
    if [ $? -eq 0 ]; then
        print_success "Tests pre-HPC pasaron ✓"
    else
        print_warning "Algunos tests fallaron. Revisa logs/pre_hpc_tests.log"
        echo -n "¿Continuar de todos modos? (y/N): "
        read -r response
        if [[ ! "$response" =~ ^[Yy]$ ]]; then
            print_error "Despliegue cancelado"
            return 1
        fi
    fi
}

submit_job() {
    print_header "SUBMISIÓN DE JOB ARRAY"
    
    echo "Submitting SLURM job array..."
    
    # Modificar submit script con configuración actualizada
    ssh -p ${HPC_PORT} ${HPC_USER}@${HPC_HOST} << EOF
cd ${HPC_WORK_DIR}

# Actualizar variables en submit script
sed -i "s|#SBATCH --partition=.*|#SBATCH --partition=${SLURM_PARTITION}|" Chronosfold/WNK/submit_umbrella_hpc.sh
sed -i "s|#SBATCH --time=.*|#SBATCH --time=${SLURM_TIME}|" Chronosfold/WNK/submit_umbrella_hpc.sh
sed -i "s|#SBATCH --mem=.*|#SBATCH --mem=${SLURM_MEM}|" Chronosfold/WNK/submit_umbrella_hpc.sh

# Activar entorno en submit script
sed -i "s|# source /path/to/your/venv/bin/activate|source ${HPC_PYTHON_ENV}/bin/activate|" Chronosfold/WNK/submit_umbrella_hpc.sh

# Submit
sbatch Chronosfold/WNK/submit_umbrella_hpc.sh

EOF
    
    if [ $? -eq 0 ]; then
        print_success "Job array submitted exitosamente"
        echo ""
        echo "Para monitorear los jobs:"
        echo "  ssh ${HPC_USER}@${HPC_HOST}"
        echo "  squeue -u ${HPC_USER}"
        echo ""
        echo "Para cancelar todos los jobs:"
        echo "  scancel -u ${HPC_USER}"
        echo ""
        echo "Logs en: ${HPC_WORK_DIR}/logs/"
    else
        print_error "Fallo en submisión de job"
        return 1
    fi
}

retrieve_results() {
    print_header "RECUPERANDO RESULTADOS"
    
    echo "Verificando si las simulaciones terminaron..."
    
    # Contar cuántas ventanas completaron
    COMPLETED=$(ssh -p ${HPC_PORT} ${HPC_USER}@${HPC_HOST} \
        "ls ${HPC_WORK_DIR}/Chronosfold/WNK/umbrella_windows/window_*/trajectory.dcd 2>/dev/null | wc -l")
    
    echo "Ventanas completadas: ${COMPLETED}/20"
    
    if [ ${COMPLETED} -eq 20 ]; then
        print_success "Todas las ventanas completaron"
        
        echo "Descargando resultados..."
        
        # Crear directorio local para resultados
        mkdir -p results_from_hpc
        
        # Descargar con rsync
        rsync -avz -e "ssh -p ${HPC_PORT}" \
            ${HPC_USER}@${HPC_HOST}:${HPC_WORK_DIR}/Chronosfold/WNK/umbrella_windows/ \
            results_from_hpc/
        
        if [ $? -eq 0 ]; then
            print_success "Resultados descargados a: results_from_hpc/"
        else
            print_error "Fallo en descarga de resultados"
            return 1
        fi
    else
        print_warning "Solo ${COMPLETED}/20 ventanas completaron"
        echo "Espera a que terminen todas antes de descargar"
        echo ""
        echo "Para verificar status:"
        echo "  ssh ${HPC_USER}@${HPC_HOST} 'squeue -u ${HPC_USER}'"
    fi
}

# ============================================================================
# MAIN
# ============================================================================

main() {
    print_header "DESPLIEGUE SSH A HPC - WNK1 UMBRELLA SAMPLING"
    
    echo "Configuración:"
    echo "  HPC: ${HPC_USER}@${HPC_HOST}:${HPC_PORT}"
    echo "  Directorio remoto: ${HPC_WORK_DIR}"
    echo "  Directorio local: ${LOCAL_DIR}"
    echo ""
    
    # Verificación de configuración
    if [ "${HPC_USER}" = "tu_usuario" ] || [ "${HPC_HOST}" = "nombre_cluster.edu" ]; then
        print_error "Debes configurar HPC_USER y HPC_HOST en este script"
        echo "Edita deploy_to_hpc.sh y ajusta las variables en la sección CONFIGURACIÓN"
        exit 1
    fi
    
    # Menú interactivo
    echo "Selecciona acción:"
    echo "  1) Deploy completo (transferir + preparar + tests + submit)"
    echo "  2) Solo transferir archivos"
    echo "  3) Solo ejecutar preparación (asume archivos ya transferidos)"
    echo "  4) Solo ejecutar tests"
    echo "  5) Solo submit job array"
    echo "  6) Solo descargar resultados"
    echo "  7) Monitorear jobs activos"
    echo ""
    echo -n "Opción [1-7]: "
    read -r option
    
    case $option in
        1)
            check_ssh_connection || exit 1
            create_remote_directory || exit 1
            transfer_files || exit 1
            run_preparation || exit 1
            run_window_generation || exit 1
            run_pre_hpc_tests || exit 1
            submit_job || exit 1
            ;;
        2)
            check_ssh_connection || exit 1
            create_remote_directory || exit 1
            transfer_files || exit 1
            ;;
        3)
            check_ssh_connection || exit 1
            run_preparation || exit 1
            run_window_generation || exit 1
            ;;
        4)
            check_ssh_connection || exit 1
            run_pre_hpc_tests || exit 1
            ;;
        5)
            check_ssh_connection || exit 1
            submit_job || exit 1
            ;;
        6)
            check_ssh_connection || exit 1
            retrieve_results || exit 1
            ;;
        7)
            check_ssh_connection || exit 1
            echo "Jobs activos:"
            ssh -p ${HPC_PORT} ${HPC_USER}@${HPC_HOST} "squeue -u ${HPC_USER}"
            ;;
        *)
            print_error "Opción inválida"
            exit 1
            ;;
    esac
    
    print_header "COMPLETADO"
}

# Ejecutar
main
