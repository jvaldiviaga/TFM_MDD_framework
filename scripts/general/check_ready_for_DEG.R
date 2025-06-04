#!/usr/bin/env Rscript

cat("🧪 Verificación de datos antes del análisis de expresión diferencial...\n")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Uso: Rscript check_ready_for_DEG.R <expresion.csv> <metadatos.csv>")
}

exp_file <- args[1]
meta_file <- args[2]
log_file <- sub(".csv", "_qc_ready_for_DEG.txt", exp_file)

# =======================
# 1. Cargar datos
# =======================
expr <- tryCatch(read.csv(exp_file, row.names = 1), error = function(e) NULL)
meta <- tryCatch(read.csv(meta_file), error = function(e) NULL)

if (is.null(expr)) stop("❌ Error cargando matriz de expresión.")
if (is.null(meta)) stop("❌ Error cargando metadatos.")

# =======================
# 2. Verificar numérico
# =======================
if (!all(sapply(expr, is.numeric))) {
  writeLines("❌ La matriz de expresión contiene valores no numéricos.", log_file)
  stop("❌ La matriz de expresión contiene valores no numéricos.")
}

# =======================
# 3. Rango de expresión
# =======================
min_val <- min(expr, na.rm = TRUE)
max_val <- max(expr, na.rm = TRUE)
is_log2 <- max_val < 30 && min_val > 0
cat(sprintf("📈 Rango expresión: [%.2f – %.2f] → %s\n",
            min_val, max_val, ifelse(is_log2, "log2 OK", "⚠️ NO log2")))

# =======================
# 4. Coincidencia de muestras
# =======================
samples_expr <- colnames(expr)
samples_meta <- meta[,1]
overlap <- intersect(samples_expr, samples_meta)
if (length(overlap) < 0.8 * length(samples_expr)) {
  writeLines("❌ Poca coincidencia entre IDs de expresión y metadatos.", log_file)
  stop("❌ Inconsistencia de muestras entre expresión y metadatos.")
} else {
  cat(sprintf("🔗 Coincidencia de muestras: %d/%d\n", length(overlap), length(samples_expr)))
}

# =======================
# 5. Detectar variable de agrupación real
# =======================
excluded_cols <- c("status", "geo_accession", "SampleID", "title", "source_name_ch1")
meta_clinical <- meta[ , !(names(meta) %in% excluded_cols), drop = FALSE]

group_col <- NULL
for (col in names(meta_clinical)) {
  unique_vals <- unique(na.omit(meta[[col]]))
  if (length(unique_vals) >= 2 && length(unique_vals) <= 5) {
    group_col <- col
    break
  }
}

if (is.null(group_col)) {
  writeLines("❌ No se encontró una variable de agrupación válida (con 2–5 niveles).", log_file)
  stop("❌ No hay variable clínica con múltiples niveles.")
}

group_values <- meta[[group_col]]
n_groups <- length(unique(group_values))
group_names <- paste(unique(group_values), collapse = ", ")
cat(sprintf("👥 Variable de agrupación: '%s' → %d grupos: %s\n", group_col, n_groups, group_names))

# =======================
# 6. Resultado final
# =======================
writeLines(c(
  "✅ Expresión numérica: OK",
  sprintf("✅ Rango expresión: [%.2f – %.2f] → %s", min_val, max_val, ifelse(is_log2, "log2 OK", "⚠️ NO log2")),
  sprintf("✅ Coincidencia de muestras: %d/%d", length(overlap), length(samples_expr)),
  sprintf("✅ Variable de agrupación detectada: %s (%d grupos: %s)", group_col, n_groups, group_names),
  if (!is_log2) "⚠️ Requiere transformación log2(expr + 1)" else "✔️ Sin transformación adicional necesaria",
  "📄 Datos listos para análisis de expresión diferencial."
), log_file)

cat("✅ Comprobaciones completadas. Puedes continuar con el análisis.\n")

