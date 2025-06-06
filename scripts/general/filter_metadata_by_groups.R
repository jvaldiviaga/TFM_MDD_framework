#!/usr/bin/env Rscript

# =====================================================
# Filtra metadata + (opcional) matriz de expresión .rds
# =====================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Uso: Rscript filter_metadata_and_expression.R <input.csv> <output.csv> --group_column COL --group_A A --group_B B [--expression_rds=path]")
}

# =====================
# PARÁMETROS OBLIGATORIOS
# =====================
input_file <- args[1]
output_file <- args[2]
group_column <- sub("--group_column=", "", args[grep("--group_column", args)])
group_A <- sub("--group_A=", "", args[grep("--group_A", args)])
group_B <- sub("--group_B=", "", args[grep("--group_B", args)])

# =====================
# PARÁMETRO OPCIONAL: expresión
# =====================
expression_rds <- NULL
if (any(grepl("--expression_rds=", args))) {
  expression_rds <- sub("--expression_rds=", "", args[grep("--expression_rds=", args)])
}

# =====================
# FILTRADO DEL METADATA
# =====================
df <- read.csv(input_file)

df[[group_column]] <- trimws(gsub('"', "", df[[group_column]]))
group_A_clean <- trimws(gsub('"', "", group_A))
group_B_clean <- trimws(gsub('"', "", group_B))

valores_unicos <- unique(df[[group_column]])
if (!(group_A_clean %in% valores_unicos) || !(group_B_clean %in% valores_unicos)) {
  cat("🔎 Valores únicos disponibles en", group_column, ":\n")
  print(valores_unicos)
  stop("❌ Uno o ambos grupos no existen tras limpieza.")
}

df_filtrado <- df[df[[group_column]] %in% c(group_A_clean, group_B_clean), ]
write.csv(df_filtrado, output_file, row.names = FALSE)

cat(sprintf("✅ Metadatos filtrados guardados en: %s (%d muestras)\n", output_file, nrow(df_filtrado)))

# =====================
# RESUMEN DE GRUPOS
# =====================
summary_file <- sub(".csv", "_group_summary.txt", output_file)
group_counts <- table(df_filtrado[[group_column]])

lines <- c(
  sprintf("Resumen de grupos en: %s", output_file),
  sprintf("Total: %d muestras", nrow(df_filtrado)),
  sprintf("Columna de agrupación: %s", group_column),
  paste(sprintf("• %s: %d muestras", names(group_counts), as.integer(group_counts)), collapse = "\n")
)
writeLines(lines, summary_file)

cat("📊 Resumen de grupos:\n")
for (g in names(group_counts)) {
  cat(sprintf("• %s: %d muestras\n", g, group_counts[[g]]))
}

# =====================
# FILTRADO DE MATRIZ DE EXPRESIÓN
# =====================
if (!is.null(expression_rds)) {
  if (!file.exists(expression_rds)) {
    stop(paste("❌ No se encuentra el archivo de expresión:", expression_rds))
  }
  
  cat("🧬 Cargando matriz de expresión:", expression_rds, "\n")
  expr <- readRDS(expression_rds)
  
  # Verificación: debe ser un data.frame con columnas = muestras
  if (!is.data.frame(expr)) stop("❌ El archivo .rds no es un data.frame válido.")
  
  # Determinar muestras en común
  common_samples <- intersect(colnames(expr), df_filtrado$SampleID)
  
  if (length(common_samples) == 0) {
    stop("❌ Ninguna muestra del metadata coincide con la matriz de expresión.")
  }
  
  # Filtrado
  expr_filtered <- expr[, common_samples]
  
  # Sobrescribir
  saveRDS(expr_filtered, expression_rds)
  cat(sprintf("✅ Matriz de expresión filtrada y sobrescrita (%d muestras)\n", length(common_samples)))
}

