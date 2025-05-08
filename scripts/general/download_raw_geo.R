#!/usr/bin/env Rscript
# ==========================================================
# Script: download_raw_geo.R
# Objetivo: Descargar archivos crudos (.tar) desde GEO automáticamente
# Uso: Se llama desde prepare_geo_dataset.sh
# ==========================================================


# 1. Cargar los paquetes necesarios
# Cargamos los paquetes necesarios sin mostrar mensajes en pantalla
suppressPackageStartupMessages({
  # Si el paquete GEOquery no está instalado, lo instalamos automáticamente
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    install.packages("BiocManager")  # Instalación condicional de BiocManager
    BiocManager::install("GEOquery") # Instalación si no está presente de GEOquery desde Bioconductor
  }
  library(GEOquery)  # Cargamos GEOquery, que permite descargar datos desde GEO (NCBI)
})

# 2. Capturar el ID del dataset desde el argumento de consola
args <- commandArgs(trailingOnly = TRUE)

# Si no se especifica un ID de dataset (ej. GSE98793), detenemos el script con error
dataset_id <- ifelse(length(args) > 0, args[1], stop("❌ Falta ID del dataset"))

# 3. Definir la carpeta de destino para la descarga
# Estructura: ~/TFM_MDD/data/transcriptomics/<GSE_ID>/
base_dir <- file.path("~/TFM_MDD/data/transcriptomics", dataset_id)

# ⚠️ AVISO IMPORTANTE:
# Actualmente este script guarda todos los datasets descargados en la carpeta 'data/transcriptomics'.
# Esto NO es adecuado para otros tipos de ómica (ej. metilación, GWAS, neuroimagen).
# 🛠️ MEJORA FUTURA:
# Permitir pasar como segundo argumento el tipo de ómica (transcriptomics, epigenomics, genomics...).
# Alternativamente, implementar detección automática posterior y reubicación según 'dataset_type.txt'.


# 4. Ejecutar la descarga del archivo .tar suplementario asociado al dataset
# Este archivo suele contener los archivos .CEL, .txt o .fastq crudos
getGEOSuppFiles(dataset_id, baseDir = base_dir)

# 5. Confirmar éxito por consola
cat("✅ Descarga cruda completada en:", base_dir, "\n")
