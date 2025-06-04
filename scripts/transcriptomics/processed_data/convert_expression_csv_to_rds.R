#!/usr/bin/env Rscript

# ===============================================================
# Script: convert_expression_csv_to_rds.R
# Objetivo: Convertir una matriz de expresión (genes x muestras)
#           desde formato .csv a .rds en la carpeta results/
# ===============================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Uso: Rscript convert_expression_csv_to_rds.R <entrada_csv> <dataset_id>")
}

input_csv <- args[1]
dataset_id <- args[2]

# === Ruta de salida: carpeta de resultados ===
output_dir <- file.path("~/TFM_MDD/results", dataset_id)
output_rds <- file.path(output_dir, "normalized_expression.rds")

# Crear carpeta si no existe
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Leer CSV y guardar como RDS
expr <- read.csv(input_csv, row.names = 1, check.names = FALSE)
saveRDS(expr, file = output_rds)

# Mostrar resumen
cat("✅ Matriz de expresión guardada como .rds en:\n", output_rds, "\n")
cat("🔍 Nº genes:", nrow(expr), " | Nº muestras:", ncol(expr), "\n")
cat("📈 Rango de expresión:", range(expr, na.rm = TRUE)[1], "-", range(expr, na.rm = TRUE)[2], "\n")
