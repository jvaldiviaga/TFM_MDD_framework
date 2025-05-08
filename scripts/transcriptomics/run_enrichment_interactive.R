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
#  LIMITACIONES ACTUALES Y MEJORAS PENDIENTES
# ------------------------------------------------------------

# 1. No valida si la subcarpeta seleccionada contiene los archivos necesarios (`DEG_results_mapped.csv`).
#    -El script puede fallar si el usuario elige una carpeta incompleta o mal procesada.

# 2. No permite seleccionar de forma explícita si se quiere hacer enriquecimiento sobre UP, DOWN o ALL.
#    -Esto depende del comportamiento del script de análisis.

# 3. No ofrece vista previa del contenido ni de los DEG detectados antes de lanzar el análisis.

# 4. No registra log propio de entrada/salida o éxito/fracaso de la ejecución.
#    -Debería guardar un `log_enrichment.txt` por dataset con fecha y carpeta elegida.

# 5. No está preparado aún para ser llamado en modo silencioso (por script externo o Snakemake),
#    ya que depende de `readline()` interactivamente.

# 6. No verifica que `enrichment_analysis.R` esté presente en la ruta especificada.

# 7. No permite lanzar enriquecimiento para varios análisis en cadena (modo batch).

# ==========================================================

# Carga de librerías
suppressPackageStartupMessages({
  library(fs)  # Para manejar rutas de forma robusta

})

cat("🧬 Asistente de enriquecimiento funcional\n")

# === Paso 1: Pedir ID del dataset ===
dataset_id <- readline("📂 Introduce el ID del dataset (ej. GSE98793): ")

# === Paso 2: Detectar subcarpetas de resultados ===
diff_dir <- file.path("~/TFM_MDD/results", dataset_id, "differential_expression")
if (!dir_exists(diff_dir)) stop("❌ No se encontró la carpeta de resultados:", diff_dir)

subcarpetas <- dir_ls(diff_dir, type = "directory")
if (length(subcarpetas) == 0) stop("❌ No se encontraron subcarpetas de expresión diferencial")

# Mostrar lista de configuraciones disponibles
cat("\n📁 Subcarpetas de resultados detectadas:\n")
for (i in seq_along(subcarpetas)) {
  cat(sprintf("  %2d) %s\n", i, path_file(subcarpetas[i])))
}

# === Paso 3: Seleccionar una configuración ===
opcion <- as.integer(readline("\n❓ Elige la subcarpeta para análisis funcional (número): "))
subcarpeta_elegida <- subcarpetas[opcion]
if (is.na(subcarpeta_elegida)) stop("❌ Opción inválida")

# === Paso 4: Ejecutar el script de enriquecimiento ===
cmd <- sprintf("Rscript ~/TFM_MDD/scripts/transcriptomics/enrichment_analysis.R %s %s", dataset_id, path_file(subcarpeta_elegida))
cat("\n🚀 Ejecutando enriquecimiento con:\n", cmd, "\n\n")
system(cmd)
