#!/usr/bin/env Rscript
# ==========================================================
# Script: preprocessing_rnaseq.R
# Objetivo: Preprocesar RNA-seq, aplicar QC adaptativo disruptivo, normalizar y preparar datos para DEG
# ==========================================================

# -------------------------------
# 1. Cargar librerías necesarias
# -------------------------------
suppressPackageStartupMessages({
  library(edgeR)              # Normalización TMM
  library(limma)              # voom
  library(ggplot2)            # Gráficos
  library(Biobase)            # ExpressionSet para arrayQualityMetrics
  library(arrayQualityMetrics) # QC avanzado
  library(pheatmap)           # Correlación
})

# -------------------------------
# 2. Definir rutas
# -------------------------------
input_dir <- get("input_dir", envir = .GlobalEnv)
dataset_id <- basename(input_dir)
results_dir <- file.path("~/TFM_MDD/results", dataset_id)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

cat("📥 Procesando RNA-seq dataset:", dataset_id, "\n")

# -------------------------------
# 3. Leer los archivos de counts
# -------------------------------

# Buscar solo archivos válidos
count_files <- list.files(file.path(input_dir, "final"),
                          pattern = "\\.txt$",
                          full.names = TRUE)

# ⚡ Filtrar archivos basura
count_files <- count_files[!grepl("dataset_type|GPL|series_matrix", count_files, ignore.case = TRUE)]

# ⚡ Adicional: sólo mantener archivos de muestras (GSM...)
count_files <- count_files[grepl("GSM", count_files)]

if (length(count_files) == 0) {
  stop("❌ No se encontraron archivos válidos de counts en la carpeta final.")
}

# Combinar counts
combined_counts <- NULL

for (file in count_files) {
  # Siempre leer como tabulado y SIN header
  data <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)
  
  if (ncol(data) >= 2) {
    gene_col <- as.character(data[, 1])
    count_col <- suppressWarnings(as.numeric(as.character(data[, 2])))
  } else {
    stop("❌ Archivo ", basename(file), " no tiene columnas suficientes para procesar.")
  }
  
  # ⚡ Filtrar genes inválidos
  valid_idx <- !(is.na(gene_col) | is.na(count_col) | gene_col %in% c("", "*", "DATA", "TYPE"))
  gene_col <- gene_col[valid_idx]
  count_col <- count_col[valid_idx]
  
  sample_name <- tools::file_path_sans_ext(basename(file))
  
  if (is.null(combined_counts)) {
    combined_counts <- data.frame(GeneID = gene_col)
  }
  
  combined_counts[[sample_name]] <- count_col
}

# ⚡ Eliminar genes duplicados
combined_counts <- combined_counts[!duplicated(combined_counts$GeneID), ]

# ⚡ Asignar rownames
rownames(combined_counts) <- combined_counts$GeneID
combined_counts$GeneID <- NULL

cat("✅ Matriz de counts combinada: ", nrow(combined_counts), " genes x ", ncol(combined_counts), " muestras\n")


# -------------------------------
# 4. Normalización preliminar (TMM + CPM sin log)
# -------------------------------
dge <- DGEList(counts = combined_counts)
dge <- calcNormFactors(dge, method = "TMM")
norm_counts <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)

# -------------------------------
# 5. QC adaptativo disruptivo
# -------------------------------
qc_dir <- file.path(results_dir, "qc_raw")
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

# 5.1 PCA de counts normalizados
pca_data <- prcomp(t(norm_counts), scale. = TRUE)
pca_scores <- pca_data$x[, 1:2]

png(file.path(qc_dir, "pca_raw.png"), width = 1000)
plot(pca_scores, pch = 19, col = "blue",
     main = "PCA de counts normalizados (RAW)",
     xlab = "PC1", ylab = "PC2")
dev.off()

# Detectar outliers en PCA
pca_dist <- sqrt(rowSums(scale(pca_scores)^2))
pca_outliers <- names(pca_dist[pca_dist > mean(pca_dist) + 2 * sd(pca_dist)])

# 5.2 Correlación muestra-muestra
cor_matrix <- cor(norm_counts, method = "pearson")
pheatmap(cor_matrix, 
         main = "Correlación entre muestras",
         filename = file.path(qc_dir, "heatmap_correlation.png"),
         show_rownames = FALSE, show_colnames = FALSE)

cor_mean <- colMeans(cor_matrix)
cor_outliers <- names(cor_mean[cor_mean < 0.8])

# 5.3 Carga total de expresión
total_counts <- colSums(combined_counts)
median_counts <- median(total_counts)
low_counts <- names(total_counts[total_counts < 0.7 * median_counts])

# -------------------------------
# 6. Decisión combinada de outliers
# -------------------------------
all_outliers <- unique(c(pca_outliers, cor_outliers, low_counts))
decision_file <- file.path(results_dir, "decision_qc_pre.txt")

if (length(all_outliers) == 0) {
  cat("✅ No se detectaron outliers graves. Normalizando todo el dataset.\n")
  final_counts <- combined_counts
  decision_text <- "Método: CPM + TMM sobre todo el dataset. Sin outliers graves."
  
} else if (length(all_outliers) < ncol(combined_counts) / 2) {
  cat("⚠️ Se detectaron algunos outliers: ", length(all_outliers), ". Filtrándolos.\n")
  keep_samples <- setdiff(colnames(combined_counts), all_outliers)
  final_counts <- combined_counts[, keep_samples]
  decision_text <- paste0("Método: CPM + TMM tras filtrar ", length(all_outliers), " muestras outliers.")
  
} else {
  cat("❌ Alta heterogeneidad: usando voom.\n")
  final_counts <- combined_counts
  decision_text <- "Método: voom normalization (alta heterogeneidad técnica)"
}

# -------------------------------
# 7. Normalización final
# -------------------------------
dge_final <- DGEList(counts = final_counts)
dge_final <- calcNormFactors(dge_final, method = "TMM")

if (grepl("voom", decision_text)) {
  v <- voom(dge_final, plot = FALSE)
  norm_final <- v$E
} else {
  norm_final <- cpm(dge_final, normalized.lib.sizes = TRUE, log = TRUE)
}

# Guardar decisión
writeLines(c(decision_text,
             paste0("Muestras procesadas: ", ncol(final_counts)),
             paste0("Fecha: ", Sys.time())),
           decision_file)

# -------------------------------
# 8. Guardar matriz normalizada
# -------------------------------
saveRDS(norm_final, file = file.path(results_dir, "normalized_expression.rds"))
cat("💾 Matriz normalizada guardada en:", file.path(results_dir, "normalized_expression.rds"), "\n")

# -------------------------------
# 9. Ejecutar arrayQualityMetrics si es posible
# -------------------------------
if (ncol(norm_final) <= 100) {
  eset <- ExpressionSet(assayData = as.matrix(norm_final))
  aqm_dir <- file.path(results_dir, "qc_arrayqualitymetrics")
  dir.create(aqm_dir, recursive = TRUE, showWarnings = FALSE)
  
  arrayQualityMetrics(expressionset = eset,
                      outdir = aqm_dir,
                      force = TRUE,
                      do.logtransform = FALSE)
  
  cat("✅ arrayQualityMetrics ejecutado. Resultados en:", aqm_dir, "\n")
} else {
  cat("⚠️ Demasiadas muestras para arrayQualityMetrics. Saltando...\n")
}

# -------------------------------
# 10. Mensaje de cierre
# -------------------------------
cat("🎯 Preprocesamiento completo para", dataset_id, "\n")
