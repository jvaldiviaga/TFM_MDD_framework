#!/usr/bin/env Rscript
# ==========================================================
# Script: preprocessing_master.R
# Objetivo: Detectar automáticamente el tipo de datos ómicos y ejecutar el script de preprocesamiento correspondiente
# Uso: (base) jessica@DESKTOP-UFC7TIA:~/TFM_MDD$ Rscript scripts/general/preprocessing_master.R data/transcriptomics/GSEXXXXX
# ==========================================================

# ==========================================================
# FORTALEZAS Y JUSTIFICACIÓN DEL FLUJO

# Este script actúa como punto de entrada principal al preprocesamiento del framework,
# automatizando la detección del tipo de datos ómicos y delegando la ejecución al script especializado
# correspondiente. Su diseño permite adaptar el pipeline a cualquier dataset sin intervención manual.

# - Detecta de forma automática el tipo de ómica (microarrays, RNA-seq, metilación, GWAS, neuroimagen),
#   a partir de patrones de nombre de archivo y del contenido de `series_matrix.txt.gz` si está presente.
# - Llama al script específico según el tipo detectado, sin necesidad de modificar rutas o nombres.
# - Verifica si el dataset ya fue procesado (normalización previa), y evita reprocesarlo innecesariamente.
# - Centraliza la lógica de preprocesamiento para todo el framework, lo que facilita la gestión en entornos
#   como Snakemake, clusters o ejecución por lotes.
# - Integra control de errores en todas las fases: ruta no encontrada, tipo no reconocido, doble ejecución, etc.
# - Registra logs estructurados de inicio, fin y duración del proceso, lo que permite trazabilidad
#   y auditoría en producción.
# - Escalable y compatible con futuras ómicas o pipelines adicionales, simplemente añadiendo nuevos scripts fuente.

# ------------------------------------------------------------
# LIMITACIONES Y MEJORAS PENDIENTES
# ------------------------------------------------------------

# 1. La detección del tipo de dato se basa únicamente en extensiones de archivo y contenido textual parcial.
#    -Esto puede ser impreciso si los archivos están mal nombrados o si hay múltiples tecnologías en la misma carpeta.

# 2. Si el script de preprocesamiento específico falla, el sistema no captura ni comunica el error de forma detallada.
#    -Debería implementar gestión de excepciones para registrar errores en el log y permitir continuar con otros datasets.

# 3. No valida si el script especializado correspondiente realmente existe antes de hacer `source(...)`.
#    -Si un tipo de dato está detectado pero no implementado aún (como `.fastq`), falla silenciosamente.

# 4. Las rutas a los scripts fuente están codificadas de forma fija.
#    -Esto limita la portabilidad y flexibilidad del framework. Mejora futura: cargar rutas desde una configuración externa.

# 5. No registra metainformación adicional sobre los archivos detectados (porcentaje de tipos, tamaño, n de muestras...).
#    -Esta información sería útil para estadísticas de uso, validación cruzada o prefiltrado automatizado.

# 6. El script no guarda información sobre el éxito o fallo del preprocesamiento en archivos individuales por dataset.
#    -Todo se resume en el log global, pero no hay evidencia directa por carpeta (e.g., `.success`, `.error`).

# 7. El log no incluye aún un hash de los archivos originales ni `sessionInfo()` del entorno base.
#    -Esto limita la reproducibilidad y trazabilidad completa en análisis con múltiples versiones o entornos.

# Estas limitaciones son propias de su función como controlador general del flujo,
#    y no afectan directamente a la calidad técnica del preprocesamiento, que reside en los scripts especializados.


# ==========================================================

# 1. Leer el argumento del usuario (ruta al dataset)
args <- commandArgs(trailingOnly = TRUE)

# Si no se especifica carpeta, se usa una ruta por defecto (útil para pruebas)
#input_dir <- ifelse(length(args) > 0, args[1],
#                    "~/TFM_MDD/data/transcriptomics/GSE98793")
#cat(paste0("📁 Carpeta detectada: ", input_dir, "\n"))

# Comprobar que el usuario ha proporcionado la ruta al dataset
# ⚠️ A partir de la versión final del framework se requiere pasar la ruta explícitamente.
# Esto evita errores silenciosos y hace el código compatible con entornos automatizados (Snakemake, Bash, etc.)
if (length(args) == 0) {
  stop("❌ Debes indicar la ruta completa al dataset como argumento (ej: data/transcriptomics/GSEXXXXX)")
}

# Guardamos la ruta de entrada como variable global
input_dir <- args[1]
cat(paste0("📁 Carpeta detectada: ", input_dir, "\n"))

# En versiones anteriores se permitía usar una ruta por defecto para pruebas,
# pero esto puede ocultar errores o limitar la generalización del framework.
# Por eso ahora se exige que el usuario proporcione explícitamente la ruta del dataset a procesar.


# 2. Comprobar que la carpeta existe
if (!dir.exists(input_dir)) {
  stop("❌ La carpeta especificada no existe: ", input_dir)
}

# 3. Listar archivos para detección del tipo de dato
file_list <- list.files(input_dir, full.names = TRUE, recursive = TRUE)
cat("📂 Archivos detectados:\n")
print(head(basename(file_list)))

# 4. Detección automática del tipo de dato ómico
tipo_dato <- NULL  # Inicializamos variable

# Detectamos si hay archivos .CEL → microarrays
if (any(grepl("\\.CEL$", file_list))) {
  tipo_dato <- "microarray"
  
# Si hay archivos .counts, .csv o .txt → posible RNA-seq  
} else if (any(grepl("\\.counts$|\\.txt$|\\.csv$", file_list))) {
  
  # Buscamos el archivo series_matrix para afinar si es microarray o RNA-seq
  series_matrix_path <- list.files(input_dir, pattern = "series_matrix.txt.gz$", full.names = TRUE)
  
  if (length(series_matrix_path) > 0) {
    cat("🔍 Detectado series_matrix. Analizando tecnología...\n")
    lines <- tolower(readLines(gzfile(series_matrix_path[1]), n = 200))  # Leemos primeras 200 líneas
    
    if (any(grepl("rna-seq|rnaseq|sequencing", lines))) {
      tipo_dato <- "rnaseq"
    } else if (any(grepl("expression profiling by array|array assay|agilent|microarray", lines))) {
      tipo_dato <- "microarray"
    } else {
      tipo_dato <- "rnaseq"  # Si no hay pistas claras, asumimos RNA-seq
    }
    cat(paste0("🧬 Tipo de datos detectado automáticamente desde series_matrix: ", tipo_dato, "\n"))
  } else {
    tipo_dato <- "rnaseq"
    cat("⚠️ No se encontró series_matrix.txt.gz. Asumiendo RNA-seq.\n")
  }

  # Detectamos archivos .idat → epigenómica (metilación)
} else if (any(grepl("\\.idat$", file_list))) {
  tipo_dato <- "metilacion"


  # Detectamos archivos .vcf, .ped, .raw, .map → genómica (GWAS)  
} else if (any(grepl("\\.vcf$|\\.ped$|\\.raw$|\\.map$", file_list))) {
  tipo_dato <- "genomica"

  # Detectamos archivos de neuroimagen → MEG o fMRI
} else if (any(grepl("\\.fif$|\\.nii.gz$", file_list))) {
  tipo_dato <- "neuroimagen"

  # Si no se detecta ningún patrón, error  
} else {
  stop("❌ No se pudo identificar el tipo de datos automáticamente.")
}


# 5. Verificar si ya fue procesado (normalización completada)
# Si existe ya el archivo 'normalized_expression.rds' en la carpeta de resultados,
# se asume que el dataset ya fue preprocesado previamente
results_dir <- gsub("data", "results", input_dir)
results_file <- file.path(results_dir, "normalized_expression.rds")
if (file.exists(results_file)) {
  cat("ℹ️ Este dataset ya ha sido procesado. Saliendo...\n")
  quit(save = "no")
}

# 6. Preparar entorno y logging
# Hacemos la variable visible globalmente (necesario si se usa en scripts fuente)
input_dir <<- input_dir  # Hacer visible globalmente para los scripts

# Extraemos el ID del dataset desde el nombre de la carpeta
dataset_id <- basename(input_dir)

# Ruta al log centralizado del preprocesamiento
log_path <- "~/TFM_MDD/results/log_preprocessing.txt"

# Guardamos hora de inicio para registrar duración total
start_time <- Sys.time()

# Línea inicial de log con información del tipo de dato
log_entry_start <- paste0("🔹 ", format(start_time, "%Y-%m-%d %H:%M:%S"),
                          " | INICIO | Dataset: ", dataset_id,
                          " | Tipo: ", toupper(tipo_dato))
cat(log_entry_start, "\n")

# 7. Ejecutar el script de preprocesamiento correspondiente
cat(paste0("🚀 Iniciando preprocesamiento para datos tipo: ", toupper(tipo_dato), "\n"))

if (tipo_dato == "microarray") {
  source("~/TFM_MDD/scripts/transcriptomics/preprocessing_microarray.R")
} else if (tipo_dato == "rnaseq") {
  source("~/TFM_MDD/scripts/transcriptomics/preprocessing_rnaseq.R")
} else if (tipo_dato == "metilacion") {
  source("~/TFM_MDD/scripts/epigenomics/preprocessing_metilacion.R")
} else if (tipo_dato == "genomica") {
  source("~/TFM_MDD/scripts/genomics/preprocessing_gwas.R")
} else if (tipo_dato == "neuroimagen") {
  source("~/TFM_MDD/scripts/neuroimaging/preprocessing_neuroimage.R")
}

# 8. Registrar finalización
end_time <- Sys.time()
duration <- round(difftime(end_time, start_time, units = "mins"), 2)

# Entrada final en el log, con duración total
log_entry_end <- paste0("✅ ", format(end_time, "%Y-%m-%d %H:%M:%S"),
                        " | FIN    | Dataset: ", dataset_id,
                        " | Duración: ", duration, " min")

# Línea separadora visual
log_separator <- paste(rep("─", 80), collapse = "")

# Escribimos log completo en archivo (append)
cat(c(log_entry_start, log_entry_end, log_separator),
    sep = "\n", file = log_path, append = TRUE)

# Mensaje final por consola
cat("✅ Preprocesamiento adaptativo completado con éxito.\n")
