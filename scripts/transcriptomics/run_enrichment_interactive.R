#!/usr/bin/env Rscript
# ==========================================================
# Script: run_enrichment_interactive.R
# Objetivo: Ejecutar enriquecimiento funcional interactivo
# Uso: Rscript scripts/run_enrichment_interactive.R
# ==========================================================

# ==========================================================
# OBJETIVO GENERAL:
# Este script funciona como un asistente interactivo para seleccionar de forma sencilla 
# los resultados de expresión diferencial (DEG) que se quieren usar en el análisis de 
# enriquecimiento funcional. Evita errores manuales al construir rutas y lanza automáticamente 
# el script `enrichment_analysis.R` con los parámetros adecuados.

# FORTALEZAS Y JUSTIFICACIÓN DEL FLUJO

# - Facilita la reutilización del pipeline para múltiples configuraciones sin duplicar scripts.
# - Permite al usuario elegir de forma asistida entre distintas comparaciones (con/sin covariables).
# - Evita errores en la ruta de entrada, ejecutando directamente el script de enriquecimiento funcional.
# - Integra la lógica dentro del flujo del framework sin exigir conocimientos de terminal.
# - Compatible con Snakemake, testing por lotes o ejecución desde entornos interactivos.

# ------------------------------------------------------------
# MEJORAS APLICADAS EN ESTA VERSIÓN
# ------------------------------------------------------------

# 1. Verifica que `DEG_results_mapped.csv` exista antes de ejecutar el script de enriquecimiento.
# 2. Registra un log por dataset cada vez que se lanza el análisis.
# 3. Controla errores si el usuario pulsa Enter sin seleccionar subcarpeta válida.

# ==========================================================

# --- Cargar librerías necesarias ---
suppressPackageStartupMessages({
  library(fs)  # Para gestión robusta de rutas y ficheros
})

cat("🧬 Asistente de enriquecimiento funcional\n")

# === Paso 1: Pedir ID del dataset ===
dataset_id <- readline("📂 Introduce el ID del dataset (ej. GSE98793): ")

# === Paso 2: Detectar subcarpetas de resultados ===
diff_dir <- file.path("~/TFM_MDD/results", dataset_id, "differential_expression")
if (!dir_exists(diff_dir)) stop("❌ No se encontró la carpeta de resultados:", diff_dir)

subcarpetas <- dir_ls(diff_dir, type = "directory")
if (length(subcarpetas) == 0) stop("❌ No se encontraron subcarpetas de expresión diferencial")

cat("\n📁 Subcarpetas de resultados detectadas:\n")
for (i in seq_along(subcarpetas)) {
  cat(sprintf("  %2d) %s\n", i, path_file(subcarpetas[i])))
}

# === Paso 3: Seleccionar una configuración ===
opcion <- as.integer(readline("\n❓ Elige la subcarpeta para análisis funcional (número): "))
if (is.na(opcion) || opcion < 1 || opcion > length(subcarpetas)) stop("❌ Opción inválida")

subcarpeta_elegida <- subcarpetas[opcion]
subcarpeta_nombre <- path_file(subcarpeta_elegida)

# === Paso 4: Verificar que existe DEG_results_mapped.csv ===
mapped_path <- file.path(diff_dir, subcarpeta_nombre, "DEG_results_mapped.csv")
if (!file.exists(mapped_path)) stop("❌ No se encontró el archivo DEG_results_mapped.csv en la carpeta seleccionada.")

# === Paso 5: Ejecutar el script de enriquecimiento ===
cmd <- sprintf("Rscript ~/TFM_MDD/scripts/transcriptomics/enrichment_analysis.R %s %s", dataset_id, subcarpeta_nombre)
cat("\n🚀 Ejecutando enriquecimiento con:\n", cmd, "\n\n")
system(cmd)

# === Paso 6: Registrar log del análisis ===
log_dir <- file.path("~/TFM_MDD/results", dataset_id, "enrichment")
dir_create(log_dir, recurse = TRUE)
log_path <- file.path(log_dir, "log_enrichment.txt")

timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_line <- sprintf("[%s] Ejecutado enriquecimiento sobre: %s / %s", timestamp, dataset_id, subcarpeta_nombre)
writeLines(log_line, log_path, sep = "\n", useBytes = TRUE, append = TRUE)

cat("\n📝 Log actualizado en:", log_path, "\n")
