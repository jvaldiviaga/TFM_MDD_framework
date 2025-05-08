#!/usr/bin/env Rscript
# ==========================================================
# Script: check_metadata_columns.R
# Objetivo: Mostrar columnas, ejemplos de valores e interpretación automática
# Uso: Rscript scripts/check_metadata_columns.R <dataset_id>
# ==========================================================

args <- commandArgs(trailingOnly = TRUE)
dataset_id <- ifelse(length(args) > 0, args[1], stop("❌ Indica el ID del dataset (ej. GSE98793)"))

meta_path <- file.path("~/TFM_MDD/data/transcriptomics", dataset_id, paste0(dataset_id, "_metadata_clean.csv"))
if (!file.exists(meta_path)) stop("❌ No se encontró el archivo:", meta_path)

metadata <- read.csv(meta_path)
colnames(metadata) <- make.names(colnames(metadata))

# === Reglas de interpretación semántica ===
interpretar_columna <- function(nombre) {
  nombre_lower <- tolower(nombre)
  if (grepl("group|diagnos|status", nombre_lower)) return("🧠 Diagnóstico clínico del sujeto (ej. MDD vs Control)")
  if (grepl("anxiety", nombre_lower)) return("😰 Presencia o ausencia de ansiedad comórbida")
  if (grepl("sex|gender", nombre_lower)) return("🚻 Sexo biológico del sujeto")
  if (grepl("age", nombre_lower)) return("🎂 Edad del sujeto (numérica)")
  if (grepl("batch|platform|array", nombre_lower)) return("🧪 Lote técnico o plataforma experimental")
  if (grepl("treatment|drug|medicated", nombre_lower)) return("💊 Tratamiento recibido o estado medicado")
  if (grepl("severity|hamd|bdi|score|scale", nombre_lower)) return("📊 Escala clínica o nivel de gravedad")
  if (grepl("smoking|alcohol|habit", nombre_lower)) return("🚬 Hábito del sujeto (tabaquismo, alcohol, etc.)")
  return("⚠️ Descripción no disponible. Revisa los valores.")
}

# === Mostrar columnas con interpretación ===
cat("\n📋 Columnas en", dataset_id, "con interpretación automática:\n\n")
for (col in colnames(metadata)) {
  cat("•", col, "→", interpretar_columna(col), "\n")
}

# === Mostrar ejemplos de valores por columna ===
cat("\n🔍 Vista previa de valores por columna:\n\n")
for (col in colnames(metadata)) {
  valores <- unique(na.omit(metadata[[col]]))
  cat("•", col, ": ", paste(head(valores, 5), collapse = ", "))
  if (length(valores) > 5) cat(", ...")
  cat("\n")
}

# === Sugerencias de columnas candidatas para comparación ===
cat("\n🤖 Columnas candidatas para comparación clínica (2–3 niveles distintos):\n\n")
for (col in colnames(metadata)) {
  valores <- unique(na.omit(metadata[[col]]))
  if (length(valores) >= 2 && length(valores) <= 3) {
    cat("✔️", col, "→", paste(valores, collapse = " / "), "\n")
  }
}

cat("\n✅ Revisión finalizada. Usa esta información para definir tu variable de comparación y posibles covariables.\n")

