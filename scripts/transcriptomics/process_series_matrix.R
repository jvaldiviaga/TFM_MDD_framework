#!/usr/bin/env Rscript
# ============================================================
# Script: process_series_matrix.R
# Objetivo: extraer expresión y metadatos desde series_matrix.txt.gz
# Adaptado al framework de Jessica (TFM_MDD) con diagnóstico por pantalla
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("❌ Debes indicar el ID del dataset (ej: GSE39653)")

dataset_id <- args[1]
input_path <- paste0("data/transcriptomics/", dataset_id, "/", dataset_id, "_series_matrix.txt.gz")
output_dir <- paste0("data/transcriptomics/", dataset_id)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("📂 Procesando dataset:", dataset_id, "\n")

# ------------------ CARGA RAW ------------------
lines <- readLines(input_path)
start_line <- grep("!series_matrix_table_begin", lines)
end_line <- grep("!series_matrix_table_end", lines)

metadata_lines <- lines[1:(start_line - 1)]
expression_lines <- lines[(start_line + 1):(end_line - 1)]

# ------------------ EXPRESIÓN ------------------
cat("🔍 Extrayendo matriz de expresión...\n")
expr_matrix <- read.delim(text = expression_lines, header = TRUE, sep = "\t", check.names = FALSE)
expr_outfile <- file.path(output_dir, paste0(dataset_id, "_expression.csv"))
write.csv(expr_matrix, expr_outfile, row.names = FALSE)
cat("✅ Matriz de expresión guardada en:", expr_outfile, "\n")

# ------------------ METADATOS ------------------
cat("🔍 Reformateando metadatos clínicos...\n")
meta_lines <- grep("^!Sample_", metadata_lines, value = TRUE)

# Agrupar campos por tipo de característica
meta_list <- list()
for (line in meta_lines) {
  parts <- strsplit(line, "\t")[[1]]
  label_raw <- parts[1]
  values <- parts[-1]
  
  label <- gsub("^!Sample_characteristics_ch1\\s*", "", label_raw)
  label <- gsub("^!Sample_", "", label)  # Por si acaso
  label <- gsub(":", "", label)
  label <- trimws(label)
  
  # Determinar nombre de la variable
  if (grepl(":", values[1])) {
    # Separar por clave: valor
    values_clean <- sapply(values, function(v) sub(".*?:", "", v))
    key <- sub(":", "", strsplit(values[1], ":")[[1]][1])
    meta_list[[key]] <- values_clean
  } else {
    meta_list[[label]] <- values
  }
}

meta_df <- as.data.frame(meta_list, stringsAsFactors = FALSE)

# ------------------ SampleID ------------------
sample_ids <- colnames(expr_matrix)[-1]
if (nrow(meta_df) == length(sample_ids)) {
  meta_df$SampleID <- sample_ids
  cols <- c("SampleID", setdiff(names(meta_df), "SampleID"))
  meta_df <- meta_df[, cols]
  cat("✅ SampleID correctamente asignado.\n")
} else {
  cat("⚠️ Nº de muestras en metadata:", nrow(meta_df), "\n")
  cat("⚠️ Nº de muestras en expresión:", length(sample_ids), "\n")
  cat("🛑 No se pudo asignar correctamente los SampleID. Se guardará sin asignar.\n")
}

# ------------------ DIAGNÓSTICO ------------------
cat("🔎 Diagnóstico de metadatos:\n")
cat("• Nº muestras:", nrow(meta_df), "\n")
cat("• Nº variables clínicas:", ncol(meta_df) - 1, "\n")

for (col in colnames(meta_df)[-1]) {
  values <- meta_df[[col]]
  if (length(unique(values)) == 1) {
    cat("⚠️ Advertencia: la columna", col, "no tiene variabilidad.\n")
  } else if (any(is.na(values)) || any(values == "")) {
    cat("⚠️ Advertencia: valores vacíos o NA detectados en", col, "\n")
  } else if (grepl("group|diagnos|status|condition", tolower(col))) {
    cat("✔️ Posible columna de grupo detectada:", col, "\n")
  }
}

# ------------------ GUARDADO ------------------
meta_outfile <- file.path(output_dir, paste0(dataset_id, "_metadata.csv"))
write.csv(meta_df, meta_outfile, row.names = FALSE)
cat("✅ Metadatos clínicos guardados en:", meta_outfile, "\n")

