#!/usr/bin/env Rscript
# ==========================================================
# Script: download_metadata_geo.R
# Objetivo: Descargar metadatos clínicos y series_matrix desde GEO
# Uso: (base) jessica@DESKTOP-UFC7TIA:~/TFM_MDD$ Rscript scripts/general/download_metadata_geo.R GSE[ID]
# Rscript scripts/general/download_metadata_geo.R GSEXXXXX
# ==========================================================


# 1. Cargar librerías necesarias
suppressPackageStartupMessages({
  # Instalamos GEOquery si no está disponible (requiere BiocManager)
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    install.packages("BiocManager")
    BiocManager::install("GEOquery")
  }
  library(GEOquery)  # GEOquery permite acceder a metadatos y datos desde GEO
})


# 2. Leer argumento del dataset
args <- commandArgs(trailingOnly = TRUE)

# Validamos que el usuario haya proporcionado un ID GEO
if (length(args) == 0) {
  stop("❌ Debes proporcionar un ID GEO (ejemplo: GSE98793)")
}

dataset_id <- args[1]  # Guardamos el ID del dataset como variable

cat(paste0("📥 Descargando metadatos para ", dataset_id, "...\n"))

# 3. Crear carpeta de salida
# Los metadatos se guardarán en una subcarpeta de transcriptómica (mejorable si se generaliza a otras ómicas)
output_dir <- file.path("~/TFM_MDD/data/transcriptomics", dataset_id)

# Creamos la carpeta destino si no existe (permite ejecución sin errores)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 4. Descargar metadatos clínicos si es posible
# Intentamos cargar el objeto GSE como ExpressionSet desde GEO
tryCatch({
  gse_object <- getGEO(dataset_id, GSEMatrix = TRUE)  # Descargamos en formato estructurado si está disponible
  
  # Verificamos si el resultado es un ExpressionSet (estructura esperada)
  if (class(gse_object)[1] == "ExpressionSet" || (is.list(gse_object) && class(gse_object[[1]])[1] == "ExpressionSet")) {
    
    # Si hay varios objetos (por ejemplo, múltiples plataformas), tomamos el primero
    if (is.list(gse_object)) {
      gse_object <- gse_object[[1]]
    }
    
    # Extraemos los metadatos (pData)
    pheno_data <- pData(gse_object)
    
    # Guardamos los metadatos en un archivo CSV
    output_file <- file.path(output_dir, paste0(dataset_id, "_metadata.csv"))
    write.csv(pheno_data, output_file, row.names = TRUE)
    
    cat(paste0("✅ Metadatos clínicos guardados en: ", output_file, "\n"))
  } else {
    cat("⚠️ No se pudo extraer metadatos clínicos estructurados (no es un ExpressionSet).\n")
  }
}, error = function(e) {
  cat("⚠️ No se pudo descargar automáticamente los metadatos clínicos.\n")
})

# 5. Descargar el archivo series_matrix completo
# Construimos la URL directa al archivo series_matrix del dataset
matrix_url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/",
                     substr(dataset_id, 1, nchar(dataset_id)-3),
                     "nnn/", dataset_id, "/matrix/", dataset_id, "_series_matrix.txt.gz")

# Aunque ya se haya descargado un ExpressionSet con getGEO(),
# descargamos también el archivo original series_matrix.txt.gz
# porque:
# - Contiene todos los metadatos sin procesar tal como fueron publicados en GEO.
# - Algunas anotaciones pueden estar ausentes o truncadas en pData().
# - Sirve como respaldo completo para parsing manual o verificación posterior.
# - Es útil si getGEO() falla o si el dataset no tiene formato GSEMatrix.


# Definimos la ruta local donde se guardará el archivo
matrix_file <- file.path(output_dir, paste0(dataset_id, "_series_matrix.txt.gz"))

cat(paste0("🌐 Descargando series_matrix desde: ", matrix_url, "\n"))

# Descargamos el archivo .gz con modo binario (wb)
tryCatch({
  download.file(matrix_url, destfile = matrix_file, mode = "wb")
  cat(paste0("✅ Archivo series_matrix guardado en: ", matrix_file, "\n"))
}, error = function(e) {
  cat("⚠️ No se pudo descargar automáticamente el series_matrix.txt.gz.\n")
  cat("🔗 Descárgalo manualmente en:\n")
  cat(paste0("   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", dataset_id, "\n"))
})

# 6. Mensaje de cierre
cat("🎯 Descarga de metadata completada.\n")

